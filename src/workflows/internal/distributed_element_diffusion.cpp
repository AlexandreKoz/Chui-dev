#include "workflows/internal/distributed_element_diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "cosmosim/parallel/composition_exchange.hpp"
#include "cosmosim/parallel/distributed_mesh.hpp"
#include "workflows/internal/amr_migration_payload.hpp"

namespace cosmosim::workflows::internal {
namespace {

struct ClaimRecord {
  DistributedElementDiffusionCoupling::FaceKey key{};
  std::int32_t other_rank = 0;
  std::int32_t reserved = 0;
  double rho_kappa_other = 0.0;
  double area_code = 0.0;
  double distance_code = 0.0;
};

struct ConductanceRecord {
  DistributedElementDiffusionCoupling::FaceKey key{};
  double conductance = 0.0;
};

static_assert(std::is_trivially_copyable_v<ClaimRecord>);
static_assert(std::is_trivially_copyable_v<ConductanceRecord>);

// Identical expression to physics::conductance (harmonic mean of rho*kappa
// times area / centre distance).
[[nodiscard]] double faceConductance(double a, double b, double area, double distance) noexcept {
  if (!(area > 0.0) || !(distance > 0.0)) return 0.0;
  const double h = (a > 0.0 && b > 0.0) ? 2.0 * a * b / (a + b) : 0.0;
  return h * area / distance;
}

template <class Record>
[[nodiscard]] std::vector<Record> exchangeRecords(const parallel::MpiContext& mpi,
                                                  const std::vector<Record>& outbound,
                                                  const std::vector<int>& destinations) {
  const std::vector<std::byte> bytes = parallel::executeBlockingDestinationRecordExchange(
      mpi,
      std::span<const std::byte>(reinterpret_cast<const std::byte*>(outbound.data()),
                                 outbound.size() * sizeof(Record)),
      sizeof(Record), destinations, 0U);
  if (bytes.size() % sizeof(Record) != 0U) {
    throw std::runtime_error("distributed element diffusion: received a truncated record stream");
  }
  std::vector<Record> inbound(bytes.size() / sizeof(Record));
  if (!inbound.empty()) std::memcpy(inbound.data(), bytes.data(), bytes.size());
  return inbound;
}

[[nodiscard]] std::uint64_t saturatingMul(std::uint64_t a, std::uint64_t b) noexcept {
  return (a != 0U && b > std::numeric_limits<std::uint64_t>::max() / a)
      ? std::numeric_limits<std::uint64_t>::max() : a * b;
}

[[nodiscard]] std::uint64_t saturatingAdd(std::uint64_t a, std::uint64_t b) noexcept {
  return a > std::numeric_limits<std::uint64_t>::max() - b ? std::numeric_limits<std::uint64_t>::max() : a + b;
}

}  // namespace

RemoteDiffusionCellGather gatherRemoteDiffusionCells(
    const core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    core::MemoryGovernor* governor,
    double scale_factor,
    bool periodic,
    const std::array<double, 3>& periodic_box_extent_comoving,
    std::uint64_t sequence) {
  RemoteDiffusionCellGather gather;
  const int rank = mpi_context.worldRank();
  const auto local_patch_records = buildMigrationAmrPatchPayloadRecords(state, rank);
  std::vector<parallel::AmrPatchCellPayloadRecord> received;
  std::vector<core::MemoryReservation> admitted;
  const auto exchange = parallel::executeBlockingDirectedAmrPatchPayloadExchange(
      mpi_context, local_patch_records,
      [&state, rank](std::span<const parallel::AmrPatchBoundaryCellRequest> requests, std::uint64_t first_record,
                     std::size_t max_records, std::vector<parallel::AmrPatchCellPayloadRecord>& output) {
        fillMigrationAmrPatchBoundaryCellPayloadChunk(state, rank, requests, first_record, max_records, output);
      },
      [&received](int, std::span<const parallel::AmrPatchCellPayloadRecord> records) {
        received.insert(received.end(), records.begin(), records.end());
      },
      [governor, &admitted](int, std::span<const parallel::AmrPatchPayloadRecord>,
                            std::span<const parallel::AmrPatchBoundaryCellRequest>, std::uint64_t record_count) {
        if (governor == nullptr || record_count == 0U) return;
        // Received payload + decoded RemoteDiffusionCell coexist.
        const std::uint64_t bytes = saturatingMul(
            record_count, sizeof(parallel::AmrPatchCellPayloadRecord) + sizeof(RemoteDiffusionCell));
        core::MemoryReservation reservation = governor->reserve(
            core::MemoryClass::kPhaseResident, bytes, "sources.element_diffusion.remote_cells");
        reservation.commit();
        admitted.push_back(std::move(reservation));
      },
      0U, sequence,
      parallel::AmrInterfaceTopologyOptions{
          .periodic_box_extent_comoving = periodic ? periodic_box_extent_comoving : std::array<double, 3>{0.0, 0.0, 0.0},
          .ghost_layers = 1U});

  std::exception_ptr failure;
  try {
    std::unordered_map<std::uint64_t, const parallel::AmrPatchPayloadRecord*> patch_by_id;
    patch_by_id.reserve(exchange.patch_payloads_received.size());
    std::vector<std::uint64_t> parent_ids;
    for (const auto& record : exchange.patch_payloads_received) {
      patch_by_id.emplace(record.patch_id, &record);
      if (record.parent_patch_id != 0U) parent_ids.push_back(record.parent_patch_id);
    }
    for (std::uint32_t p = 0; p < state.patches.size(); ++p) {
      if (state.patches.parent_patch_id[p] != 0U) parent_ids.push_back(state.patches.parent_patch_id[p]);
    }
    std::sort(parent_ids.begin(), parent_ids.end());
    gather.cells.reserve(received.size());
    for (const auto& cell : received) {
      const auto it = patch_by_id.find(cell.patch_id);
      if (it == patch_by_id.end()) {
        throw std::runtime_error("distributed diffusion: remote boundary cell has no patch metadata (fail closed)");
      }
      // Production AMR stores leaves only; a patch that is somebody's parent
      // is covered and never a diffusion neighbour.
      if (std::binary_search(parent_ids.begin(), parent_ids.end(), cell.patch_id)) continue;
      const parallel::AmrPatchPayloadRecord& patch = *it->second;
      RemoteDiffusionCell out;
      out.gas_cell_id = cell.gas_cell_id;
      out.owner_rank = cell.owner_rank;
      out.patch_dims = {patch.cell_dim_x, patch.cell_dim_y, patch.cell_dim_z};
      if (out.patch_dims[0] == 0U || out.patch_dims[1] == 0U || out.patch_dims[2] == 0U) {
        throw std::runtime_error("distributed diffusion: remote patch has a zero dimension");
      }
      out.patch_origin_comoving = {patch.origin_x_comoving, patch.origin_y_comoving, patch.origin_z_comoving};
      out.cell_width_comoving = {patch.extent_x_comoving / static_cast<double>(out.patch_dims[0]),
                                 patch.extent_y_comoving / static_cast<double>(out.patch_dims[1]),
                                 patch.extent_z_comoving / static_cast<double>(out.patch_dims[2])};
      const std::uint32_t off = cell.local_cell_offset;
      out.index = {off % out.patch_dims[0], (off / out.patch_dims[0]) % out.patch_dims[1],
                   off / (out.patch_dims[0] * out.patch_dims[1])};
      if (out.index[2] >= out.patch_dims[2]) {
        throw std::runtime_error("distributed diffusion: remote cell offset outside its patch");
      }
      out.velocity_peculiar = {cell.velocity_x_peculiar, cell.velocity_y_peculiar, cell.velocity_z_peculiar};
      const double volume_factor = scale_factor * scale_factor * scale_factor;
      out.level = static_cast<std::int32_t>(patch.level);
      out.density_physical_code = cell.density_code > 0.0
          ? cell.density_code / volume_factor
          : 0.0;
      gather.cells.push_back(out);
    }
    // Deterministic order (gas_cell_id), duplicates (multiple images or
    // requests) collapsed: one virtual row per remote cell.
    std::sort(gather.cells.begin(), gather.cells.end(),
              [](const RemoteDiffusionCell& l, const RemoteDiffusionCell& r) { return l.gas_cell_id < r.gas_cell_id; });
    gather.cells.erase(std::unique(gather.cells.begin(), gather.cells.end(),
                                   [](const RemoteDiffusionCell& l, const RemoteDiffusionCell& r) {
                                     return l.gas_cell_id == r.gas_cell_id;
                                   }),
                       gather.cells.end());
  } catch (...) {
    failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(failure, "distributed diffusion remote cell decode");

  // Global physical envelope of all owned patches (periodic wrap planes).
  std::array<double, 6> envelope{};
  for (int c = 0; c < 3; ++c) {
    envelope[c] = std::numeric_limits<double>::infinity();
    envelope[3 + c] = std::numeric_limits<double>::infinity();  // stores -max
  }
  for (std::uint32_t p = 0; p < state.patches.size(); ++p) {
    if (state.patches.owning_rank[p] != static_cast<std::uint32_t>(rank)) continue;
    const std::array<double, 3> origin{state.patches.origin_x_comoving[p], state.patches.origin_y_comoving[p],
                                       state.patches.origin_z_comoving[p]};
    const std::array<double, 3> extent{state.patches.extent_x_comoving[p], state.patches.extent_y_comoving[p],
                                       state.patches.extent_z_comoving[p]};
    const std::array<std::uint32_t, 3> dims{state.patches.cell_dim_x[p], state.patches.cell_dim_y[p],
                                            state.patches.cell_dim_z[p]};
    for (int c = 0; c < 3; ++c) {
      if (dims[c] == 0U) continue;
      // Same arithmetic as the topology builder's envelope.
      const double o = origin[c] * scale_factor;
      const double d = extent[c] / static_cast<double>(dims[c]) * scale_factor;
      envelope[c] = std::min(envelope[c], o);
      envelope[3 + c] = std::min(envelope[3 + c], -(o + d * static_cast<double>(dims[c])));
    }
  }
  for (double& value : envelope) value = mpi_context.allreduceMinDouble(value);
  bool finite = true;
  for (const double value : envelope) finite = finite && std::isfinite(value);
  if (finite) {
    gather.global_envelope_physical = std::array<double, 6>{
        envelope[0], envelope[1], envelope[2], -envelope[3], -envelope[4], -envelope[5]};
  }
  return gather;
}

DistributedElementDiffusionCoupling::DistributedElementDiffusionCoupling(
    const parallel::MpiContext& mpi_context,
    core::MemoryGovernor* governor,
    const core::SimulationState& state,
    std::span<const RemoteDiffusionFace> remote_faces,
    std::span<const double> rho_kappa_code,
    std::size_t lane_upper_bound)
    : m_mpi(mpi_context),
      m_governor(governor),
      m_state(state),
      m_faces(remote_faces),
      m_rho_kappa(rho_kappa_code),
      m_lane_upper_bound(lane_upper_bound) {}

void DistributedElementDiffusionCoupling::reduceFailure(std::exception_ptr local_failure, const char* stage) {
  m_mpi.rethrowCollectivePreparationFailure(local_failure, stage);
}

double DistributedElementDiffusionCoupling::globalMin(double value) {
  return m_mpi.allreduceMinDouble(value);
}

void DistributedElementDiffusionCoupling::globalSum(std::span<double> values) {
  m_mpi.allreduceSumDoublesInPlace(values);
}

void DistributedElementDiffusionCoupling::setup() {
  const int rank = m_mpi.worldRank();
  std::vector<ClaimRecord> claims;
  std::vector<int> claim_destinations;
  std::exception_ptr failure;
  try {
    // P0-49-01: lane vectors travel as chunked schema-2 records, so every
    // constructible layout (ion-rich C/O/Si/Fe + dust + SGS) is admitted;
    // only a layout above the static bound is refused (impossible by
    // CompositionLayout::make).
    const std::size_t chunks_per_row = parallel::compositionLaneChunkCount(m_lane_upper_bound);
    // Admission: tables + one subcycle's state/delta send+receive buffers.
    const std::uint64_t n = m_faces.size();
    std::uint64_t bytes = saturatingMul(n, sizeof(AuthorityFace) + sizeof(DependentFace) +
                                               2U * sizeof(ClaimRecord) + 2U * sizeof(ConductanceRecord) +
                                               sizeof(std::pair<std::uint32_t, int>) +
                                               sizeof(std::pair<std::uint64_t, std::uint32_t>));
    // Per subcycle at most one state and one delta lane vector per face on
    // each side (chunked), plus the reassembled flat lanes.
    bytes = saturatingAdd(bytes, saturatingMul(2U, parallel::compositionLaneExchangeAdmissionBytes(
                                                       n, m_lane_upper_bound)));
    (void)chunks_per_row;
    if (m_governor != nullptr && bytes != 0U) {
      m_reservation = m_governor->reserve(core::MemoryClass::kPhaseResident, bytes,
                                          "sources.element_diffusion.remote");
      m_reservation.commit();
    }
    m_admitted_bytes = bytes;
    for (const RemoteDiffusionFace& face : m_faces) {
      if (face.local_row >= m_rho_kappa.size()) {
        throw std::runtime_error("distributed diffusion remote face references an invalid local row");
      }
      if (face.local_gas_cell_id == face.remote_gas_cell_id) {
        throw std::runtime_error("distributed diffusion remote face joins a cell to itself");
      }
      const std::uint8_t periodic = face.periodic ? 1U : 0U;
      if (face.local_gas_cell_id < face.remote_gas_cell_id) {
        AuthorityFace a;
        a.key = FaceKey{face.local_gas_cell_id, face.remote_gas_cell_id, face.axis, face.local_side, periodic, 0U};
        a.local_row = face.local_row;
        a.other_rank = face.remote_rank;
        a.area_code = face.area_code;
        a.distance_code = face.center_distance_code;
        m_authority.push_back(a);
      } else {
        DependentFace d;
        d.key = FaceKey{face.remote_gas_cell_id, face.local_gas_cell_id, face.axis,
                        static_cast<std::int8_t>(-face.local_side), periodic, 0U};
        d.local_row = face.local_row;
        d.authority_rank = face.remote_rank;
        m_dependent.push_back(d);
        ClaimRecord claim;
        claim.key = d.key;
        claim.other_rank = rank;
        claim.rho_kappa_other = m_rho_kappa[face.local_row];
        claim.area_code = face.area_code;
        claim.distance_code = face.center_distance_code;
        claims.push_back(claim);
        claim_destinations.push_back(face.remote_rank);
      }
    }
    std::sort(m_authority.begin(), m_authority.end(),
              [](const AuthorityFace& l, const AuthorityFace& r) { return l.key < r.key; });
    std::sort(m_dependent.begin(), m_dependent.end(),
              [](const DependentFace& l, const DependentFace& r) { return l.key < r.key; });
    for (std::size_t i = 1; i < m_authority.size(); ++i) {
      if (m_authority[i].key == m_authority[i - 1U].key) {
        throw std::runtime_error("distributed diffusion: duplicate authority interface face");
      }
    }
  } catch (...) {
    failure = std::current_exception();
  }
  reduceFailure(failure, "distributed diffusion topology admission");
  std::vector<ClaimRecord> inbound = exchangeRecords(m_mpi, claims, claim_destinations);

  std::vector<ConductanceRecord> replies;
  std::vector<int> reply_destinations;
  try {
    std::sort(inbound.begin(), inbound.end(), [](const ClaimRecord& l, const ClaimRecord& r) { return l.key < r.key; });
    if (inbound.size() != m_authority.size()) {
      throw std::runtime_error(
          "distributed diffusion: asymmetric interface topology (authority faces " +
          std::to_string(m_authority.size()) + ", claims " + std::to_string(inbound.size()) +
          "); refusing a non-conservative distributed flux");
    }
    for (std::size_t i = 0; i < inbound.size(); ++i) {
      AuthorityFace& a = m_authority[i];
      const ClaimRecord& c = inbound[i];
      if (!(c.key == a.key) || c.other_rank != a.other_rank) {
        throw std::runtime_error("distributed diffusion: interface claim does not match the authority face");
      }
      const double area_scale = std::max(std::abs(a.area_code), 1.0e-300);
      const double dist_scale = std::max(std::abs(a.distance_code), 1.0e-300);
      if (std::abs(c.area_code - a.area_code) > 1.0e-10 * area_scale ||
          std::abs(c.distance_code - a.distance_code) > 1.0e-10 * dist_scale) {
        throw std::runtime_error("distributed diffusion: interface geometry differs between owners");
      }
      a.conductance = faceConductance(m_rho_kappa[a.local_row], c.rho_kappa_other, a.area_code, a.distance_code);
      replies.push_back(ConductanceRecord{a.key, a.conductance});
      reply_destinations.push_back(a.other_rank);
    }
  } catch (...) {
    failure = std::current_exception();
  }
  reduceFailure(failure, "distributed diffusion interface validation");
  std::vector<ConductanceRecord> conductances = exchangeRecords(m_mpi, replies, reply_destinations);
  try {
    std::sort(conductances.begin(), conductances.end(),
              [](const ConductanceRecord& l, const ConductanceRecord& r) { return l.key < r.key; });
    if (conductances.size() != m_dependent.size()) {
      throw std::runtime_error("distributed diffusion: missing conductance replies for dependent faces");
    }
    for (std::size_t i = 0; i < conductances.size(); ++i) {
      if (!(conductances[i].key == m_dependent[i].key)) {
        throw std::runtime_error("distributed diffusion: conductance reply does not match a dependent face");
      }
      m_dependent[i].conductance = conductances[i].conductance;
    }
    for (const DependentFace& d : m_dependent) {
      if (d.conductance > 0.0) m_exports.emplace_back(d.local_row, d.authority_rank);
    }
    std::sort(m_exports.begin(), m_exports.end());
    m_exports.erase(std::unique(m_exports.begin(), m_exports.end()), m_exports.end());
    for (const DependentFace& d : m_dependent) {
      m_dependent_rows.emplace_back(m_state.gas_cells.gas_cell_id[d.local_row], d.local_row);
    }
    std::sort(m_dependent_rows.begin(), m_dependent_rows.end());
    m_dependent_rows.erase(std::unique(m_dependent_rows.begin(), m_dependent_rows.end()), m_dependent_rows.end());
  } catch (...) {
    failure = std::current_exception();
  }
  reduceFailure(failure, "distributed diffusion conductance agreement");
  m_setup_done = true;
}

void DistributedElementDiffusionCoupling::addRemoteConductance(std::span<double> conductance_sum) {
  if (!m_setup_done) throw std::logic_error("distributed diffusion coupling used before setup()");
  for (const AuthorityFace& a : m_authority) {
    if (a.conductance > 0.0 && a.local_row < conductance_sum.size()) conductance_sum[a.local_row] += a.conductance;
  }
  for (const DependentFace& d : m_dependent) {
    if (d.conductance > 0.0 && d.local_row < conductance_sum.size()) conductance_sum[d.local_row] += d.conductance;
  }
}

void DistributedElementDiffusionCoupling::addRemoteFaceTransfers(
    std::span<const double> values, std::span<const double> gas_mass, std::size_t lane_count, double dt_sub,
    std::span<double> delta) {
  if (!m_setup_done) throw std::logic_error("distributed diffusion coupling used before setup()");
  const int rank = m_mpi.worldRank();
  const std::size_t rows = gas_mass.size();
  std::exception_ptr failure;
  // Collective sequence: every rank calls this once per subcycle.
  const std::uint64_t sequence = ++m_exchange_sequence;
  const auto identity = [&](parallel::CompositionLaneStream stream) {
    return parallel::CompositionLaneStreamIdentity{
        .stream = stream,
        .layout_signature = m_state.gas_composition.enabled() ? m_state.gas_composition.descriptor().layout_signature
                                                              : 0U,
        .lane_count = static_cast<std::uint32_t>(lane_count),
        .exchange_sequence = sequence};
  };
  std::vector<parallel::CompositionLaneChunkRecord> states;
  std::vector<int> state_destinations;
  try {
    if (lane_count > m_lane_upper_bound || values.size() != lane_count * rows || delta.size() != lane_count * rows) {
      throw std::invalid_argument("distributed diffusion transfer: inconsistent lane/row extents");
    }
    const auto state_identity = identity(parallel::CompositionLaneStream::kDiffusionState);
    std::array<double, parallel::k_max_composition_layout_lanes> lane_values{};
    states.reserve(m_exports.size() * parallel::compositionLaneChunkCount(lane_count));
    for (const auto& [row, dest] : m_exports) {
      for (std::size_t k = 0; k < lane_count; ++k) lane_values[k] = values[k * rows + row];
      parallel::appendCompositionLaneChunks(state_identity, m_state.gas_cells.gas_cell_id[row], rank, dest,
                                            gas_mass[row], std::span<const double>(lane_values.data(), lane_count),
                                            states, state_destinations);
    }
  } catch (...) {
    failure = std::current_exception();
  }
  reduceFailure(failure, "distributed diffusion state export");
  const auto state_chunks = parallel::exchangeCompositionLaneChunks(
      m_mpi, states, state_destinations, sequence, m_governor, "sources.element_diffusion.state_exchange");
  states.clear();

  std::vector<parallel::CompositionLaneChunkRecord> delta_chunks_out;
  std::vector<int> delta_destinations;
  try {
    const parallel::AssembledCompositionLanes remote_states = parallel::assembleCompositionLaneChunks(
        state_chunks.records, identity(parallel::CompositionLaneStream::kDiffusionState), rank);
    const auto find_state = [&remote_states](std::uint64_t id) -> std::size_t {
      const auto it = std::lower_bound(remote_states.gas_cell_id.begin(), remote_states.gas_cell_id.end(), id);
      if (it == remote_states.gas_cell_id.end() || *it != id) {
        throw std::runtime_error("distributed diffusion: authority face has no fresh remote state (fail closed)");
      }
      return static_cast<std::size_t>(it - remote_states.gas_cell_id.begin());
    };
    // One outbound transfer vector per (other cell, other rank), summed in
    // authority-face key order (deterministic for a given decomposition).
    std::vector<std::pair<std::uint64_t, std::size_t>> slot;  // other id -> outbound index
    std::vector<std::uint64_t> out_id;
    std::vector<int> out_rank;
    std::vector<double> out_delta;  // flat, lane_count per outbound entry
    for (const AuthorityFace& a : m_authority) {
      if (!(a.conductance > 0.0)) continue;
      const std::size_t other = find_state(a.key.other_id);
      const double other_gas_mass = remote_states.scalar[other];
      const auto other_values = remote_states.lanesOf(other);
      if (!(other_gas_mass > 0.0) || !(gas_mass[a.local_row] > 0.0)) {
        throw std::runtime_error("distributed diffusion: interface cell without composition gas mass");
      }
      auto it = std::lower_bound(slot.begin(), slot.end(), std::make_pair(a.key.other_id, std::size_t{0}));
      std::size_t index = 0U;
      if (it != slot.end() && it->first == a.key.other_id) {
        index = it->second;
      } else {
        index = out_id.size();
        out_id.push_back(a.key.other_id);
        out_rank.push_back(a.other_rank);
        out_delta.resize(out_delta.size() + lane_count, 0.0);
        slot.insert(it, std::make_pair(a.key.other_id, index));
      }
      const std::size_t row = a.local_row;
      for (std::size_t k = 0; k < lane_count; ++k) {
        // Same expression as the local face loop (a = authority, b = other).
        const double x_a = values[k * rows + row] / gas_mass[row];
        const double x_b = other_values[k] / other_gas_mass;
        const double transfer = dt_sub * a.conductance * (x_a - x_b);
        delta[k * rows + row] -= transfer;
        out_delta[index * lane_count + k] += transfer;
      }
    }
    const auto delta_identity = identity(parallel::CompositionLaneStream::kDiffusionDelta);
    for (std::size_t i = 0; i < out_id.size(); ++i) {
      parallel::appendCompositionLaneChunks(
          delta_identity, out_id[i], rank, out_rank[i], 0.0,
          std::span<const double>(out_delta.data() + i * lane_count, lane_count), delta_chunks_out,
          delta_destinations);
    }
  } catch (...) {
    failure = std::current_exception();
  }
  reduceFailure(failure, "distributed diffusion authority transfers");
  const auto imported_chunks = parallel::exchangeCompositionLaneChunks(
      m_mpi, delta_chunks_out, delta_destinations, sequence, m_governor, "sources.element_diffusion.delta_exchange");
  try {
    const parallel::AssembledCompositionLanes imported = parallel::assembleCompositionLaneChunks(
        imported_chunks.records, identity(parallel::CompositionLaneStream::kDiffusionDelta), rank);
    for (std::size_t i = 0; i < imported.size(); ++i) {
      const auto it = std::lower_bound(m_dependent_rows.begin(), m_dependent_rows.end(),
                                       std::make_pair(imported.gas_cell_id[i], std::uint32_t{0}));
      if (it == m_dependent_rows.end() || it->first != imported.gas_cell_id[i]) {
        throw std::runtime_error("distributed diffusion: imported transfer for an unknown interface row");
      }
      const std::size_t row = it->second;
      const auto d = imported.lanesOf(i);
      for (std::size_t k = 0; k < lane_count; ++k) delta[k * rows + row] += d[k];
    }
  } catch (...) {
    failure = std::current_exception();
  }
  reduceFailure(failure, "distributed diffusion transfer import");
}

}  // namespace cosmosim::workflows::internal
