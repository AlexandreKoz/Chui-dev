// CHUI-MHD-REF-v1 workflow composition (pass 02). UNVALIDATED source.
// See workflows/internal/mhd_runtime_bridge.hpp for the ownership contract.

#include "workflows/internal/mhd_runtime_bridge.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/stage_transaction.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"
#include "cosmosim/mhd/mhd_initial_conditions.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"

namespace cosmosim::workflows::internal {
namespace {

constexpr std::uint64_t k_fnv_offset = 1469598103934665603ULL;
constexpr std::uint64_t k_fnv_prime = 1099511628211ULL;
constexpr std::size_t k_descriptor_words = 12U;
constexpr std::uint64_t k_mib = 1024ULL * 1024ULL;

void mix(std::uint64_t& hash, std::uint64_t word) noexcept {
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (word >> (8 * byte)) & 0xFFULL;
    hash *= k_fnv_prime;
  }
}

[[noreturn]] void fail(std::string_view phase, const std::string& detail) {
  throw std::runtime_error("MHD_RUNTIME (" + std::string(phase) + "): " + detail);
}

[[nodiscard]] std::uint64_t checkedBytesMultiply(std::uint64_t lhs, std::uint64_t rhs,
                                                 std::string_view context) {
  if (rhs != 0U && lhs > std::numeric_limits<std::uint64_t>::max() / rhs) {
    fail("admission", std::string(context) + " overflows 64-bit byte arithmetic");
  }
  return lhs * rhs;
}

[[nodiscard]] std::int64_t roundToLattice(double value, std::string_view what) {
  const double rounded = std::nearbyint(value);
  if (!std::isfinite(value) || std::abs(value - rounded) > 1.0e-6 ||
      std::abs(rounded) > 9.0e15) {
    std::ostringstream stream;
    stream << what << " does not lie on the integer AMR lattice (value " << value << ")";
    fail("lattice", stream.str());
  }
  return static_cast<std::int64_t>(rounded);
}

[[nodiscard]] core::StageRejectionKind rejectionKindFor(mhd::MhdErrorCode code) noexcept {
  switch (code) {
    case mhd::MhdErrorCode::kNonPositiveState:
    case mhd::MhdErrorCode::kStiffnessBudgetExceeded:
      return core::StageRejectionKind::kNumericalRetry;
    case mhd::MhdErrorCode::kMemoryBudgetExceeded:
      return core::StageRejectionKind::kResourceRefused;
    default:
      return core::StageRejectionKind::kInvalidState;
  }
}

struct Descriptor {
  std::uint64_t patch_id = 0U;
  std::int32_t level = 0;
  std::array<std::int64_t, 3> cells{};
  std::array<double, 3> origin{};
  std::array<double, 3> extent{};
  std::int32_t owner_rank = 0;
  std::size_t local_patch_row = static_cast<std::size_t>(-1);
};

}  // namespace

// ---------------------------------------------------------------------------
// Transport
// ---------------------------------------------------------------------------

int MpiMhdHaloTransport::rank() const noexcept { return std::max(m_mpi.worldRank(), 0); }

int MpiMhdHaloTransport::size() const noexcept {
  return m_mpi.isEnabled() ? std::max(m_mpi.worldSize(), 1) : 1;
}

std::vector<std::vector<std::uint8_t>> MpiMhdHaloTransport::exchange(
    const std::vector<std::vector<std::uint8_t>>& send) {
  if (size() == 1) {
    return std::vector<std::vector<std::uint8_t>>(1U);
  }
  return parallel::exchangeBoundedAlltoallBytes(m_mpi, send);
}

std::uint64_t MpiMhdHaloTransport::allreduceMaxUint64(std::uint64_t value) {
  return m_mpi.allreduceMaxUint64(value);
}

double MpiMhdHaloTransport::allreduceMinDouble(double value) {
  return m_mpi.allreduceMinDouble(value);
}

double MpiMhdHaloTransport::allreduceMaxDouble(double value) {
  return -m_mpi.allreduceMinDouble(-value);
}

void MpiMhdHaloTransport::allreduceSumDoubles(std::span<double> values) {
  m_mpi.allreduceSumDoublesInPlace(values);
}

// ---------------------------------------------------------------------------
// Bridge
// ---------------------------------------------------------------------------

MhdRuntimeBridge::MhdRuntimeBridge(const core::SimulationConfig& config,
                                   const parallel::MpiContext& mpi_context,
                                   core::MemoryGovernor* memory_governor,
                                   const core::UnitSystem& units)
    : m_config(config),
      m_mpi(mpi_context),
      m_governor(memory_governor),
      m_units(units),
      m_transport(mpi_context) {
  if (!config.mhd.enabled) {
    throw std::logic_error("MhdRuntimeBridge constructed for an MHD-disabled run");
  }
  // Pass 04: direct construction must obey the same precise admission as the
  // capability report and the HydroAmrRuntime factory (composition, observables,
  // linear-response neutrinos, regrid, rebalance, ...).
  if (const std::string refusal = mhdProductionDispatchRefusal(config, core::ModePolicy{});
      !refusal.empty()) {
    throw std::invalid_argument("MHD_UNSUPPORTED_PHYSICS: " + refusal);
  }
  m_policy = mhd::mhdNumericsPolicyFromConfig(config.mhd);
  if (config.mhd.nonideal.enabled) {
    m_provider = mhd::makeConfiguredCoefficientProvider(config.mhd);
    mhd::MhdPhysicalPropertyModel properties{};
    properties.units = &m_units;
    // Part-1 reference mean molecular weight, the runtime-wide convention of
    // temperatureKelvinPerInternalEnergyCode (composition-dependent mu is a
    // separate gated closure; see pass-02 ledger R-MU).
    properties.mean_molecular_weight = hydro::k_part1_reference_mean_molecular_weight;
    m_rate = std::make_unique<mhd::MhdNonIdealRateOperator>(
        m_provider.get(), properties, m_policy.stiffness, config.mhd.nonideal.ohmic_enabled,
        config.mhd.nonideal.hall_enabled, config.mhd.nonideal.ambipolar_enabled);
  }
}

MhdRuntimeBridge::~MhdRuntimeBridge() = default;

std::uint64_t MhdRuntimeBridge::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = m_stepper ? m_stepper->ownedCapacityBytes() : 0U;
  for (const LocalPatch& patch : m_patches) {
    bytes += patch.gas.ownedCapacityBytes() + patch.rows.capacity() * sizeof(std::uint32_t);
  }
  return bytes;
}

mhd::MhdLevelAdvanceOptions MhdRuntimeBridge::advanceOptions() const {
  mhd::MhdLevelAdvanceOptions options{};
  options.nonideal = m_rate.get();
  // Synchronized leaf AMR reconciles coarse/fine fluxes and EMFs at every
  // stage: no pending register, so neither edge nor boundary-flux integrals
  // are required by the production path (both remain available for audits).
  options.produce_edge_integrals = false;
  // Pass 04 (M03-12): accepted open/wall boundary gas and non-ideal Poynting
  // work enter the conservation budget. The integrals live in the stepper's
  // admitted level-resident boundary planes (O(patch surface)); they are
  // requested only when the global domain has a non-periodic side, so a fully
  // periodic run keeps its pass-02 memory plan.
  bool physical_boundary = false;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    physical_boundary = physical_boundary ||
        m_level_domain.lower[axis] != mhd::MhdBoundaryKind::kPeriodic ||
        m_level_domain.upper[axis] != mhd::MhdBoundaryKind::kPeriodic;
  }
  options.produce_boundary_flux_integrals = physical_boundary;
  return options;
}

std::uint64_t MhdRuntimeBridge::topologyFingerprint(const core::SimulationState& state) const {
  std::uint64_t hash = k_fnv_offset;
  const core::PatchSoa& patches = state.patches;
  mix(hash, patches.size());
  for (std::size_t p = 0; p < patches.size(); ++p) {
    if (patches.cell_count[p] == 0U) {
      continue;
    }
    mix(hash, patches.patch_id[p]);
    mix(hash, static_cast<std::uint64_t>(static_cast<std::int64_t>(patches.level[p])));
    mix(hash, patches.first_cell[p]);
    mix(hash, patches.cell_count[p]);
    mix(hash, std::bit_cast<std::uint64_t>(patches.origin_x_comoving[p]));
    mix(hash, std::bit_cast<std::uint64_t>(patches.origin_y_comoving[p]));
    mix(hash, std::bit_cast<std::uint64_t>(patches.origin_z_comoving[p]));
    mix(hash, std::bit_cast<std::uint64_t>(patches.extent_x_comoving[p]));
    mix(hash, std::bit_cast<std::uint64_t>(patches.extent_y_comoving[p]));
    mix(hash, std::bit_cast<std::uint64_t>(patches.extent_z_comoving[p]));
    mix(hash, patches.cell_dim_x[p]);
    mix(hash, patches.cell_dim_y[p]);
    mix(hash, patches.cell_dim_z[p]);
  }
  mix(hash, state.cells.size());
  return hash;
}

void MhdRuntimeBridge::buildDirectory(const core::SimulationState& state,
                                      const MhdDomainSpec& domain) {
  const int my_rank = m_transport.rank();
  std::vector<std::uint8_t> local_bytes;
  std::vector<Descriptor> local;
  std::exception_ptr failure;
  try {
    const core::PatchSoa& patches = state.patches;
    if (!patches.isConsistent()) {
      fail("directory", "PatchSoa lanes are inconsistent");
    }
    for (std::size_t p = 0; p < patches.size(); ++p) {
      if (patches.cell_count[p] == 0U) {
        continue;
      }
      if (patches.cell_dim_x[p] == 0U || patches.cell_dim_y[p] == 0U || patches.cell_dim_z[p] == 0U ||
          !(patches.extent_x_comoving[p] > 0.0) || !(patches.extent_y_comoving[p] > 0.0) ||
          !(patches.extent_z_comoving[p] > 0.0)) {
        fail("directory",
             "MHD requires explicit restart-authoritative patch geometry on every leaf patch");
      }
      Descriptor descriptor{};
      descriptor.patch_id = patches.patch_id[p];
      descriptor.level = patches.level[p];
      descriptor.cells = {patches.cell_dim_x[p], patches.cell_dim_y[p], patches.cell_dim_z[p]};
      descriptor.origin = {patches.origin_x_comoving[p], patches.origin_y_comoving[p],
                           patches.origin_z_comoving[p]};
      descriptor.extent = {patches.extent_x_comoving[p], patches.extent_y_comoving[p],
                           patches.extent_z_comoving[p]};
      descriptor.owner_rank = my_rank;
      descriptor.local_patch_row = p;
      local.push_back(descriptor);
      const std::array<std::uint64_t, k_descriptor_words> words{
          descriptor.patch_id,
          static_cast<std::uint64_t>(static_cast<std::int64_t>(descriptor.level)),
          static_cast<std::uint64_t>(descriptor.cells[0]),
          static_cast<std::uint64_t>(descriptor.cells[1]),
          static_cast<std::uint64_t>(descriptor.cells[2]),
          std::bit_cast<std::uint64_t>(descriptor.origin[0]),
          std::bit_cast<std::uint64_t>(descriptor.origin[1]),
          std::bit_cast<std::uint64_t>(descriptor.origin[2]),
          std::bit_cast<std::uint64_t>(descriptor.extent[0]),
          std::bit_cast<std::uint64_t>(descriptor.extent[1]),
          std::bit_cast<std::uint64_t>(descriptor.extent[2]),
          static_cast<std::uint64_t>(static_cast<std::int64_t>(descriptor.owner_rank))};
      const std::size_t offset = local_bytes.size();
      local_bytes.resize(offset + sizeof(words));
      std::memcpy(local_bytes.data() + offset, words.data(), sizeof(words));
    }
  } catch (...) {
    failure = std::current_exception();
  }
  m_mpi.rethrowCollectivePreparationFailure(failure, "MHD leaf directory preparation");
  // Pass 04 (M03-04): admit the gathered directory, the decoded descriptor
  // table and the replacement row maps BEFORE they are allocated. Counts are
  // agreed collectively (zero-patch ranks contribute zero), so every rank
  // admits from the same global descriptor count. The previous row maps stay
  // charged until the replacement is built (old + new coexistence).
  std::uint64_t local_cells = 0U;
  for (const Descriptor& d : local) {
    local_cells += static_cast<std::uint64_t>(d.cells[0] * d.cells[1] * d.cells[2]);
  }
  const std::uint64_t global_descriptors =
      m_transport.size() > 1 ? m_mpi.allreduceSumUint64(local.size()) : local.size();
  core::MemoryReservation directory_reservation{};
  core::MemoryReservation row_map_reservation{};
  failure = nullptr;
  try {
    if (m_governor != nullptr) {
      const std::uint64_t per_descriptor = k_descriptor_words * sizeof(std::uint64_t) +
          sizeof(Descriptor) + sizeof(mhd::MhdRemotePatchDescriptor) +
          6U * sizeof(std::uint64_t) /* canonical epoch tuple */;
      const std::uint64_t directory_bytes = core::checkedMemoryBytesAdd(
          checkedBytesMultiply(global_descriptors, per_descriptor,
                                           "MHD directory bytes"),
          local_bytes.capacity(), "MHD directory bytes");
      directory_reservation = m_governor->reserve(core::MemoryClass::kCommunication,
                                                  directory_bytes, "mhd_runtime.directory");
      const std::uint64_t row_bytes = core::checkedMemoryBytesAdd(
          checkedBytesMultiply(local_cells, sizeof(std::uint32_t),
                                           "MHD row-map bytes"),
          checkedBytesMultiply(local.size(), sizeof(LocalPatch), "MHD row-map bytes"),
          "MHD row-map bytes");
      row_map_reservation = m_governor->reserve(core::MemoryClass::kPersistentCache, row_bytes,
                                                "mhd_runtime.row_maps");
    }
  } catch (...) {
    failure = std::current_exception();
  }
  m_mpi.rethrowCollectivePreparationFailure(failure, "MHD directory and row-map admission");
  if (directory_reservation.pending()) {
    directory_reservation.commit();
  }
  if (row_map_reservation.pending()) {
    row_map_reservation.commit();
  }
  const std::vector<std::uint8_t> gathered =
      m_transport.size() > 1 ? m_mpi.allgatherBytesBounded(local_bytes) : local_bytes;
  std::uint64_t directory_digest = 0U;
  failure = nullptr;
  try {
    if (gathered.size() % (k_descriptor_words * sizeof(std::uint64_t)) != 0U) {
      fail("directory", "gathered leaf directory has a truncated record");
    }
    std::vector<Descriptor> all;
    const std::size_t count = gathered.size() / (k_descriptor_words * sizeof(std::uint64_t));
    all.reserve(count);
    for (std::size_t r = 0; r < count; ++r) {
      std::array<std::uint64_t, k_descriptor_words> words{};
      std::memcpy(words.data(), gathered.data() + r * sizeof(words), sizeof(words));
      Descriptor descriptor{};
      descriptor.patch_id = words[0];
      descriptor.level = static_cast<std::int32_t>(static_cast<std::int64_t>(words[1]));
      descriptor.cells = {static_cast<std::int64_t>(words[2]), static_cast<std::int64_t>(words[3]),
                          static_cast<std::int64_t>(words[4])};
      descriptor.origin = {std::bit_cast<double>(words[5]), std::bit_cast<double>(words[6]),
                           std::bit_cast<double>(words[7])};
      descriptor.extent = {std::bit_cast<double>(words[8]), std::bit_cast<double>(words[9]),
                           std::bit_cast<double>(words[10])};
      descriptor.owner_rank = static_cast<std::int32_t>(static_cast<std::int64_t>(words[11]));
      all.push_back(descriptor);
    }
    if (all.empty()) {
      fail("directory", "no leaf patch exists on any rank");
    }
    // Lattice: base level = coarsest leaf level; spacings must follow ratio 2.
    std::int32_t base_level = all.front().level;
    for (const Descriptor& d : all) {
      base_level = std::min(base_level, d.level);
    }
    std::array<double, 3> base_spacing{};
    bool have_base = false;
    for (const Descriptor& d : all) {
      if (d.level != base_level) {
        continue;
      }
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        const double spacing = d.extent[axis] / static_cast<double>(d.cells[axis]);
        if (!have_base) {
          base_spacing[axis] = spacing;
        } else if (std::abs(spacing - base_spacing[axis]) > 1.0e-9 * base_spacing[axis]) {
          fail("directory", "base-level leaves disagree on the cell spacing");
        }
      }
      have_base = true;
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (!(domain.extent_comoving[axis] > 0.0) || !std::isfinite(domain.origin_comoving[axis])) {
        fail("directory", "invalid global AMR hierarchy domain");
      }
      m_level_domain.cells[axis] =
          roundToLattice(domain.extent_comoving[axis] / base_spacing[axis], "domain extent");
      m_level_domain.origin_comoving[axis] = domain.origin_comoving[axis];
    }
    m_level_domain.base_level = base_level;
    const mhd::MhdPatchBoundary boundary = mhd::mhdBoundaryFromCondition(domain.boundary);
    m_level_domain.lower = boundary.lower;
    m_level_domain.upper = boundary.upper;
    m_patches.clear();
    m_remote.clear();
    // Pass 04 (M03-01): canonical GLOBAL directory digest. Integer tuples
    // (level, patch_id, lattice origin, cells, owner) sorted by
    // (level, patch_id); the allgather order and any floating bits are
    // irrelevant, so every rank computes the same value.
    std::vector<std::array<std::int64_t, 9>> canonical;
    canonical.reserve(all.size());
    std::int32_t finest_level = base_level;
    for (const Descriptor& d : all) {
      finest_level = std::max(finest_level, d.level);
    }
    m_finest_level = finest_level;
    m_base_spacing = base_spacing;
    for (const Descriptor& d : all) {
      const double factor = std::ldexp(1.0, -(d.level - base_level));
      std::array<double, 3> spacing{};
      std::array<std::int64_t, 3> origin_cell{};
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        spacing[axis] = base_spacing[axis] * factor;
        const double measured = d.extent[axis] / static_cast<double>(d.cells[axis]);
        if (std::abs(measured - spacing[axis]) > 1.0e-9 * spacing[axis]) {
          fail("directory", "leaf spacing is not base spacing / 2^(level - base)");
        }
        origin_cell[axis] = roundToLattice(
            (d.origin[axis] - domain.origin_comoving[axis]) / spacing[axis], "patch origin");
      }
      std::uint64_t generation = k_fnv_offset;
      mix(generation, static_cast<std::uint64_t>(static_cast<std::int64_t>(d.level)));
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        mix(generation, static_cast<std::uint64_t>(origin_cell[axis]));
        mix(generation, static_cast<std::uint64_t>(d.cells[axis]));
      }
      generation |= 1U;
      canonical.push_back({static_cast<std::int64_t>(d.level),
                           static_cast<std::int64_t>(d.patch_id), origin_cell[0], origin_cell[1],
                           origin_cell[2], d.cells[0], d.cells[1], d.cells[2],
                           static_cast<std::int64_t>(d.owner_rank)});
      if (d.owner_rank != my_rank) {
        mhd::MhdRemotePatchDescriptor remote{};
        remote.patch_id = d.patch_id;
        remote.geometry_generation = generation;
        remote.level = d.level;
        remote.origin_cell = origin_cell;
        remote.cells = d.cells;
        remote.owner_rank = d.owner_rank;
        m_remote.push_back(remote);
        continue;
      }
      LocalPatch patch{};
      patch.patch_id = d.patch_id;
      patch.geometry_generation = generation;
      patch.level = d.level;
      patch.origin_cell = origin_cell;
      patch.geometry.nx = static_cast<std::size_t>(d.cells[0]);
      patch.geometry.ny = static_cast<std::size_t>(d.cells[1]);
      patch.geometry.nz = static_cast<std::size_t>(d.cells[2]);
      patch.geometry.dx_comoving = spacing[0];
      patch.geometry.dy_comoving = spacing[1];
      patch.geometry.dz_comoving = spacing[2];
      patch.geometry.origin_x_comoving =
          domain.origin_comoving[0] + static_cast<double>(origin_cell[0]) * spacing[0];
      patch.geometry.origin_y_comoving =
          domain.origin_comoving[1] + static_cast<double>(origin_cell[1]) * spacing[1];
      patch.geometry.origin_z_comoving =
          domain.origin_comoving[2] + static_cast<double>(origin_cell[2]) * spacing[2];
      // Patch-local cell -> dense row from cell centres (no ordering assumption).
      const std::size_t cells = patch.geometry.cellCount();
      patch.rows.assign(cells, std::numeric_limits<std::uint32_t>::max());
      // Locate the PatchSoa row of this local descriptor.
      std::size_t patch_row = static_cast<std::size_t>(-1);
      for (const Descriptor& mine : local) {
        if (mine.patch_id == d.patch_id) {
          patch_row = mine.local_patch_row;
        }
      }
      if (patch_row == static_cast<std::size_t>(-1)) {
        fail("directory", "local leaf missing from the local descriptor list");
      }
      const std::uint32_t first = state.patches.first_cell[patch_row];
      const std::uint32_t count = state.patches.cell_count[patch_row];
      if (count != cells) {
        fail("directory", "patch cell count disagrees with its cell dimensions");
      }
      for (std::uint32_t offset = 0; offset < count; ++offset) {
        const std::uint32_t row = first + offset;
        const std::array<double, 3> center{state.cells.center_x_comoving[row],
                                           state.cells.center_y_comoving[row],
                                           state.cells.center_z_comoving[row]};
        const std::array<double, 3> origin{patch.geometry.origin_x_comoving,
                                           patch.geometry.origin_y_comoving,
                                           patch.geometry.origin_z_comoving};
        std::array<std::int64_t, 3> ijk{};
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          ijk[axis] =
              roundToLattice((center[axis] - origin[axis]) / spacing[axis] - 0.5, "cell centre");
          if (ijk[axis] < 0 || ijk[axis] >= d.cells[axis]) {
            fail("directory", "cell centre lies outside its owning patch");
          }
        }
        const std::size_t index = static_cast<std::size_t>(
            ijk[0] + d.cells[0] * (ijk[1] + d.cells[1] * ijk[2]));
        if (patch.rows[index] != std::numeric_limits<std::uint32_t>::max()) {
          fail("directory", "two gas rows map to one patch cell");
        }
        patch.rows[index] = row;
      }
      m_patches.push_back(std::move(patch));
    }
    std::sort(canonical.begin(), canonical.end(),
              [](const auto& lhs, const auto& rhs) {
                return std::tie(lhs[0], lhs[1]) < std::tie(rhs[0], rhs[1]);
              });
    for (std::size_t r = 1; r < canonical.size(); ++r) {
      if (canonical[r][1] == canonical[r - 1][1]) {
        fail("directory", "duplicate patch_id across ranks in the global leaf directory");
      }
    }
    std::uint64_t digest = k_fnv_offset;
    mix(digest, canonical.size());
    mix(digest, static_cast<std::uint64_t>(static_cast<std::int64_t>(base_level)));
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      mix(digest, static_cast<std::uint64_t>(m_level_domain.cells[axis]));
    }
    for (const auto& tuple : canonical) {
      for (const std::int64_t word : tuple) {
        mix(digest, static_cast<std::uint64_t>(word));
      }
    }
    directory_digest = digest;
  } catch (...) {
    failure = std::current_exception();
  }
  m_mpi.rethrowCollectivePreparationFailure(failure, "MHD leaf directory and lattice");
  // The gathered bytes and decode tables are released on return; only the
  // retained row maps stay charged (old maps were already cleared above).
  m_row_map_reservation = std::move(row_map_reservation);
  directory_reservation.release();
  // Collectively committed generation: buildDirectory is entered by every
  // rank the same number of times (rebuild decisions are reduced), so the
  // counter is rank-invariant. The wire epoch distinguishes successive
  // hierarchies even when their descriptors coincide.
  ++m_hierarchy_generation;
  agreeHierarchyValue(directory_digest, "global leaf directory digest");
  m_directory_digest = directory_digest;
  std::uint64_t epoch = directory_digest;
  mix(epoch, m_hierarchy_generation);
  m_hierarchy_epoch = epoch | 1U;
  agreeHierarchyValue(m_hierarchy_epoch, "hierarchy wire epoch");
}

void MhdRuntimeBridge::agreeHierarchyValue(std::uint64_t value, std::string_view what) {
  if (m_transport.size() <= 1) {
    return;
  }
  const std::uint64_t high = m_mpi.allreduceMaxUint64(value);
  const std::uint64_t high_complement = m_mpi.allreduceMaxUint64(~value);
  const std::uint64_t local_disagree = (high != value || high_complement != ~value) ? 1U : 0U;
  // Every rank learns the same verdict before anyone throws (no partial
  // mutation, no rank left waiting in a later collective).
  if (m_mpi.allreduceMaxUint64(local_disagree) != 0U) {
    fail("epoch", "ranks disagree on the " + std::string(what) +
                      "; the hierarchy is not globally consistent (stale or corrupt directory)");
  }
}

void MhdRuntimeBridge::seedFreshField(core::SimulationState& state, const MhdDomainSpec& domain,
                                      double scale_factor) {
  const core::MhdConfig::Io& io = m_config.mhd.io;
  core::MagneticFaceSoa storage{};
  core::MemoryReservation layout_reservation{};
  core::MemoryReservation seed_scratch_reservation{};
  std::uint64_t process_baseline_before = 0U;
  std::uint64_t old_owned = 0U;
  std::exception_ptr failure;
  try {
    if (io.initial_field_seed == core::MhdInitialFieldSeedConfigMode::kNone) {
      fail("initial conditions",
           "fresh start without a magnetic state: set mhd_io.initial_field_seed (uniform, "
           "field_loop, circular_alfven, stochastic) or restart from a cosmosim_restart_v26 "
           "checkpoint; cell-centred fields are never projected into face flux");
    }
    // Pass 04 (M03-04): checked upper bound of the authoritative layout (three
    // flux lanes + records + alignment slack) and of the per-patch seed edge
    // scratch, admitted BEFORE any allocation. The layout is canonical
    // persistent state; the edge scratch is a scratch arena held only while
    // seeding.
    std::uint64_t face_values = 0U;
    std::uint64_t max_edge_values = 0U;
    for (const LocalPatch& patch : m_patches) {
      const mhd::MhdPatchGeometry& g = patch.geometry;
      face_values = core::checkedMemoryBytesAdd(
          face_values, g.faceCountX() + g.faceCountY() + g.faceCountZ(), "MHD seed layout");
      const std::uint64_t edges = static_cast<std::uint64_t>(g.nx) * (g.ny + 1U) * (g.nz + 1U) +
          static_cast<std::uint64_t>(g.nx + 1U) * g.ny * (g.nz + 1U) +
          static_cast<std::uint64_t>(g.nx + 1U) * (g.ny + 1U) * g.nz;
      max_edge_values = std::max(max_edge_values, edges);
    }
    constexpr std::uint64_t k_lane_alignment_slack = 3U * 256U;
    const std::uint64_t layout_bytes = core::checkedMemoryBytesAdd(
        core::checkedMemoryBytesAdd(
            checkedBytesMultiply(face_values, sizeof(double), "MHD seed layout"),
            checkedBytesMultiply(m_patches.size(), sizeof(core::MagneticFacePatchRecord),
                                 "MHD seed layout"),
            "MHD seed layout"),
        k_lane_alignment_slack + 128U, "MHD seed layout");
    const std::uint64_t scratch_bytes = core::checkedMemoryBytesAdd(
        checkedBytesMultiply(max_edge_values, sizeof(double), "MHD seed scratch"),
        k_lane_alignment_slack, "MHD seed scratch");
    if (m_governor != nullptr) {
      process_baseline_before = m_governor->snapshot().baseline_owned_bytes;
      old_owned = core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(state));
      if (process_baseline_before < old_owned) {
        // The process baseline has not been measured yet (very first
        // allocation of the run): the replacement is then charged relative to
        // the state it replaces only.
        process_baseline_before = old_owned;
      }
      layout_reservation = m_governor->reserve(core::MemoryClass::kCanonicalPersistent,
                                               layout_bytes, "mhd_runtime.seed_layout");
      seed_scratch_reservation = m_governor->reserve(core::MemoryClass::kScratchArena,
                                                     scratch_bytes, "mhd_runtime.seed_scratch");
      layout_reservation.commit();
      seed_scratch_reservation.commit();
    }
    std::vector<mhd::MhdAdapterPatchBinding> bindings;
    bindings.reserve(m_patches.size());
    for (const LocalPatch& patch : m_patches) {
      bindings.push_back({patch.patch_id, patch.geometry_generation, patch.level, patch.geometry});
    }
    storage = mhd::buildMagneticFaceLayout(
        bindings, m_provider ? m_provider->tableSha256().value_or(std::string()) : std::string());
    const std::array<double, 3> uniform{io.initial_field_uniform_b_x_comoving,
                                        io.initial_field_uniform_b_y_comoving,
                                        io.initial_field_uniform_b_z_comoving};
    mhd::MhdEdgePotential potential{};
    switch (io.initial_field_seed) {
      case core::MhdInitialFieldSeedConfigMode::kNone:
      case core::MhdInitialFieldSeedConfigMode::kUniform:
        break;
      case core::MhdInitialFieldSeedConfigMode::kFieldLoop:
        potential = mhd::fieldLoopPotential(
            io.initial_field_amplitude_comoving, io.initial_field_length_comoving,
            {io.initial_field_center_x_comoving, io.initial_field_center_y_comoving,
             io.initial_field_center_z_comoving});
        break;
      case core::MhdInitialFieldSeedConfigMode::kCircularAlfven:
        potential = mhd::circularAlfvenPotential(io.initial_field_amplitude_comoving,
                                                 io.initial_field_length_comoving);
        break;
      case core::MhdInitialFieldSeedConfigMode::kStochastic:
        potential = mhd::stochasticSeedPotential(
            domain.extent_comoving, io.initial_field_k_max, io.initial_field_amplitude_comoving,
            io.initial_field_spectral_index, io.initial_field_random_seed);
        break;
    }
    // Pass 04 (M03-10): one common geometric edge authority on the composite
    // hierarchy. Every edge line integral is evaluated as the pairwise
    // (binary-tree) sum of the integrals over its unit sub-edges of the
    // FINEST lattice present anywhere in the global hierarchy, with sub-edge
    // start points computed from integer finest-lattice coordinates. A coarse
    // edge therefore equals the floating-point sum of its two fine halves with
    // the identical association, so the coarse interface face flux C(A)
    // equals the sum of its four fine faces to roundoff and the strict entry
    // consistency check holds for non-polynomial potentials as well.
    if (potential) {
      potential = mhd::compositeEdgePotential(std::move(potential), domain.origin_comoving,
                                              m_base_spacing,
                                              m_finest_level - m_level_domain.base_level);
    }
    double divergence = 0.0;
    for (std::size_t index = 0; index < m_patches.size(); ++index) {
      const LocalPatch& patch = m_patches[index];
      const mhd::MhdAdapterPatchBinding& binding = bindings[index];
      mhd::MhdFaceFluxMutableView view = mhd::makeMutableMagneticFaceView(storage, binding, index);
      divergence = std::max(divergence, mhd::seedFaceFlux(patch.geometry, uniform, potential,
                                                          view.phi_x, view.phi_y, view.phi_z));
    }
    if (!(divergence <= m_policy.divergence_tolerance)) {
      std::ostringstream stream;
      stream << "seeded magnetic field violates the strict-CT divergence bound (" << divergence
             << ")";
      fail("initial conditions", stream.str());
    }
    storage.divergence_max_normalized = divergence;
    (void)scale_factor;  // Seeds are B_c (comoving); B_phys = B_c / a_initial^2.
    if (m_governor != nullptr) {
      const std::uint64_t owned = static_cast<std::uint64_t>(
          (storage.face_flux_x_comoving.capacity() + storage.face_flux_y_comoving.capacity() +
           storage.face_flux_z_comoving.capacity()) * sizeof(double) +
          storage.patches.capacity() * sizeof(core::MagneticFacePatchRecord));
      if (owned > layout_reservation.bytes()) {
        fail("admission", "seeded magnetic layout capacity exceeds its pre-admitted bound");
      }
    }
  } catch (...) {
    failure = std::current_exception();
  }
  seed_scratch_reservation.release();
  // Collective success precedes mutation of the authoritative state: a seed
  // failure on any rank leaves every rank's (empty) magnetic state intact.
  m_mpi.rethrowCollectivePreparationFailure(failure, "MHD fresh-start magnetic seed");
  const std::uint64_t new_owned = static_cast<std::uint64_t>(
      (storage.face_flux_x_comoving.capacity() + storage.face_flux_y_comoving.capacity() +
       storage.face_flux_z_comoving.capacity()) * sizeof(double) +
      storage.patches.capacity() * sizeof(core::MagneticFacePatchRecord));
  const std::uint64_t replaced = static_cast<std::uint64_t>(
      (state.magnetic_faces.face_flux_x_comoving.capacity() +
       state.magnetic_faces.face_flux_y_comoving.capacity() +
       state.magnetic_faces.face_flux_z_comoving.capacity()) * sizeof(double) +
      state.magnetic_faces.patches.capacity() * sizeof(core::MagneticFacePatchRecord));
  state.magnetic_faces = std::move(storage);  // nothrow move of the vectors
  if (layout_reservation.valid() && layout_reservation.committed()) {
    // Hand the admitted bytes to the process baseline without a transient
    // double count (the old, empty state is retired in the same step).
    const std::uint64_t base = process_baseline_before >= replaced
        ? process_baseline_before - replaced
        : 0U;
    layout_reservation.reconcileBaselineOwnedAndRelease(base + new_owned);
  }
  (void)old_owned;
}

void MhdRuntimeBridge::bindFaces(core::SimulationState& state) {
  mhd::validateMagneticFaceSoa(state.magnetic_faces);
  // Provider identity is part of restart admission (REF 8.1/13.2).
  const std::string expected_sha =
      m_provider ? m_provider->tableSha256().value_or(std::string()) : std::string();
  if (state.magnetic_faces.coefficient_table_sha256 != expected_sha) {
    fail("bind", "magnetic state carries coefficient table SHA-256 '" +
                     state.magnetic_faces.coefficient_table_sha256 +
                     "' but the configured provider is '" + expected_sha + "'");
  }
  if (state.magnetic_faces.patches.size() != m_patches.size()) {
    fail("bind", "magnetic face records do not correspond one-to-one with the local leaves "
                 "(regrid or migration without magnetic transfer)");
  }
  for (LocalPatch& patch : m_patches) {
    const std::size_t record = mhd::findMagneticFacePatch(state.magnetic_faces, patch.patch_id,
                                                          patch.geometry_generation);
    if (record == state.magnetic_faces.patches.size()) {
      fail("bind", "no magnetic face record for leaf patch " + std::to_string(patch.patch_id) +
                       " at its current geometry generation (stale magnetic topology)");
    }
    const core::MagneticFacePatchRecord& r = state.magnetic_faces.patches[record];
    // Pass 04 (M03-13): full three-axis geometry (dimensions, widths and
    // origins) - the geometry generation alone is a lattice digest and does
    // not cover the comoving widths/origins recorded in the restart record.
    const std::array<double, 3> record_width{r.cell_width_x_comoving, r.cell_width_y_comoving,
                                             r.cell_width_z_comoving};
    const std::array<double, 3> record_origin{r.origin_x_comoving, r.origin_y_comoving,
                                              r.origin_z_comoving};
    const std::array<double, 3> width{patch.geometry.dx_comoving, patch.geometry.dy_comoving,
                                      patch.geometry.dz_comoving};
    const std::array<double, 3> origin{patch.geometry.origin_x_comoving,
                                       patch.geometry.origin_y_comoving,
                                       patch.geometry.origin_z_comoving};
    bool geometry_matches = r.level == patch.level && r.nx == patch.geometry.nx &&
        r.ny == patch.geometry.ny && r.nz == patch.geometry.nz;
    for (std::size_t axis = 0; axis < 3U && geometry_matches; ++axis) {
      geometry_matches = std::isfinite(record_width[axis]) && std::isfinite(record_origin[axis]) &&
          std::abs(record_width[axis] - width[axis]) <= 1.0e-12 * width[axis] &&
          std::abs(record_origin[axis] - origin[axis]) <=
              1.0e-9 * std::max(width[axis], std::abs(origin[axis]));
    }
    if (!geometry_matches) {
      fail("bind", "magnetic face record geometry (level, dims, x/y/z widths or origins) "
                   "disagrees with leaf patch " + std::to_string(patch.patch_id));
    }
    patch.face_record = record;
  }
}

void MhdRuntimeBridge::ensureMagneticState(core::SimulationState& state,
                                           const core::IntegratorState& integrator,
                                           const MhdDomainSpec& domain, bool allow_seed) {
  const bool has_leaves = [&]() {
    for (std::size_t p = 0; p < state.patches.size(); ++p) {
      if (state.patches.cell_count[p] != 0U) {
        return true;
      }
    }
    return false;
  }();
  const std::uint64_t local_empty = has_leaves && state.magnetic_faces.empty() ? 1U : 0U;
  const std::uint64_t local_present = !state.magnetic_faces.empty() ? 1U : 0U;
  const std::uint64_t ranks_empty =
      m_transport.size() > 1 ? m_mpi.allreduceSumUint64(local_empty) : local_empty;
  const std::uint64_t ranks_present =
      m_transport.size() > 1 ? m_mpi.allreduceSumUint64(local_present) : local_present;
  if (ranks_empty != 0U && ranks_present != 0U) {
    fail("admission", "some ranks hold a magnetic state and others do not");
  }
  if (ranks_empty != 0U) {
    if (!allow_seed) {
      fail("admission",
           "the magnetic state was not initialized at the legal pre-step boundary "
           "(HydroAmrRuntime::initializeMagneticState); a timestep candidate or kHydroUpdate "
           "never seeds the authoritative field");
    }
    if (integrator.step_index != 0U) {
      fail("admission",
           "magnetized run past step 0 without an authoritative magnetic state (restart from a "
           "non-magnetic checkpoint?); fields are never reseeded after the fresh start");
    }
    buildDirectory(state, domain);
    seedFreshField(state, domain, integrator.current_scale_factor);
    m_stepper.reset();
  }
  bind(state, domain);
}

void MhdRuntimeBridge::bind(core::SimulationState& state, const MhdDomainSpec& domain) {
  const std::uint64_t fingerprint = topologyFingerprint(state);
  const bool same = m_stepper != nullptr && fingerprint == m_fingerprint &&
      m_bound_face_data == state.magnetic_faces.face_flux_x_comoving.data() &&
      m_bound_identity_generation == state.gas_cell_identity.generation();
  const std::uint64_t rebuild_local = same ? 0U : 1U;
  const std::uint64_t rebuild =
      m_transport.size() > 1 ? m_mpi.allreduceMaxUint64(rebuild_local) : rebuild_local;
  if (rebuild == 0U) {
    return;
  }
  m_stepper.reset();
  m_staging_reservation.release();
  buildDirectory(state, domain);
  std::exception_ptr failure;
  try {
    bindFaces(state);
    // Gas staging (kPhaseResident): six comoving conserved lanes per local
    // cell, admitted BEFORE allocation.
    std::uint64_t staging = 0U;
    for (const LocalPatch& patch : m_patches) {
      staging += 6U * mhd::mhdAlignedLaneBytes(patch.geometry.cellCount());
    }
    if (m_governor != nullptr && staging != 0U) {
      m_staging_reservation =
          m_governor->reserve(core::MemoryClass::kPhaseResident, staging, "mhd_runtime.gas_staging");
    }
    for (LocalPatch& patch : m_patches) {
      patch.gas.resize(patch.geometry.cellCount());
    }
    if (m_staging_reservation.pending()) {
      m_staging_reservation.commit();
    }
  } catch (...) {
    failure = std::current_exception();
  }
  m_mpi.rethrowCollectivePreparationFailure(failure, "MHD face binding and gas staging admission");
  failure = nullptr;
  try {
    std::vector<mhd::MhdLevelPatchBinding> bindings;
    bindings.reserve(m_patches.size());
    for (std::size_t index = 0; index < m_patches.size(); ++index) {
      LocalPatch& patch = m_patches[index];
      mhd::MhdAdapterPatchBinding adapter{patch.patch_id, patch.geometry_generation, patch.level,
                                          patch.geometry};
      mhd::MhdLevelPatchBinding binding{};
      binding.patch_id = patch.patch_id;
      binding.geometry_generation = patch.geometry_generation;
      binding.level = patch.level;
      binding.origin_cell = patch.origin_cell;
      binding.geometry = patch.geometry;
      binding.gas = &patch.gas;
      binding.faces = mhd::makeMutableMagneticFaceView(state.magnetic_faces, adapter,
                                                       patch.face_record);
      // Stable address: m_patches is never resized while the stepper lives.
      binding.faces.geometry = &patch.geometry;
      binding.owner_rank = m_transport.rank();
      bindings.push_back(binding);
    }
    mhd::MhdLevelStepperConfig stepper_config{};
    stepper_config.batch_patches =
        static_cast<std::uint32_t>(m_config.mhd.performance.active_batch_target_patches);
    stepper_config.retry_mode = mhd::MhdRetryMode::kExternal;
    stepper_config.transport = &m_transport;
    stepper_config.require_memory_admission = m_governor != nullptr;
    // Pass 04 (M03-01): globally agreed epoch, never the local fingerprint.
    stepper_config.hierarchy_epoch = m_hierarchy_epoch;
    stepper_config.topology_governor = m_governor;
    stepper_config.comm_budget_bytes = m_config.mhd.performance.mhd_comm_budget_mib_per_rank * k_mib;
    m_stepper = std::make_unique<mhd::MhdLevelStepper>(m_policy, m_level_domain, std::move(bindings),
                                                       m_remote, stepper_config);
  } catch (...) {
    failure = std::current_exception();
  }
  m_mpi.rethrowCollectivePreparationFailure(failure, "MHD hierarchy topology");
  // The stepper's own digest of its global directory must agree on every rank
  // (it includes the owner ranks the bridge assigned).
  agreeHierarchyValue(m_stepper->directoryDigest(), "stepper directory digest");
  failure = nullptr;
  try {
    if (m_governor != nullptr) {
      mhd::MhdResourceLimits limits{};
      limits.scratch_bytes_per_active_cell =
          m_config.mhd.performance.mhd_scratch_max_bytes_per_active_cell;
      limits.comm_budget_bytes_per_rank =
          m_config.mhd.performance.mhd_comm_budget_mib_per_rank * k_mib;
      limits.comm_max_bytes_per_rank = m_config.mhd.performance.mhd_comm_max_mib_per_rank * k_mib;
      limits.io_chunk_bytes_per_rank = m_config.mhd.performance.mhd_io_chunk_mib_per_rank * k_mib;
      limits.active_batch_target_patches =
          static_cast<std::uint32_t>(m_config.mhd.performance.active_batch_target_patches);
      m_stepper->admitMemory(*m_governor, advanceOptions(), limits);
    }
  } catch (...) {
    failure = std::current_exception();
  }
  m_mpi.rethrowCollectivePreparationFailure(failure, "MHD memory admission");
  m_fingerprint = fingerprint;
  m_bound_face_data = state.magnetic_faces.face_flux_x_comoving.data();
  m_bound_identity_generation = state.gas_cell_identity.generation();
}

void MhdRuntimeBridge::loadGas(const core::SimulationState& state, double scale_factor) {
  std::exception_ptr failure;
  try {
    for (LocalPatch& patch : m_patches) {
      mhd::MhdAdapterPatchBinding adapter{patch.patch_id, patch.geometry_generation, patch.level,
                                          patch.geometry};
      const mhd::MhdFaceFluxView view =
          mhd::makeMagneticFaceView(state.magnetic_faces, adapter, patch.face_record);
      mhd::loadMhdGasFromState(state, patch.geometry, patch.rows, view, scale_factor, patch.gas);
    }
  } catch (...) {
    failure = std::current_exception();
  }
  m_mpi.rethrowCollectivePreparationFailure(failure, "MHD gas staging load");
}

double MhdRuntimeBridge::stableTimeStep(core::SimulationState& state,
                                        const core::IntegratorState& integrator,
                                        const MhdDomainSpec& domain, bool cosmological) {
  ensureMagneticState(state, integrator, domain);
  const double a = integrator.current_scale_factor > 0.0 ? integrator.current_scale_factor : 1.0;
  loadGas(state, a);
  mhd::MhdStageContext interval{};
  interval.time_code = integrator.current_time_code;
  interval.scale_factor = a;
  interval.scale_factor_end = a;
  interval.comoving_coordinates = cosmological;
  try {
    return m_stepper->estimateStableTimeStep(interval, advanceOptions());
  } catch (const mhd::MhdError& error) {
    throw core::StageTransactionRejected(rejectionKindFor(error.code()),
                                         std::string("MHD timestep candidate: ") + error.what());
  }
}

void MhdRuntimeBridge::advance(core::StepContext& context, const MhdDomainSpec& domain,
                               double temperature_k_per_internal_energy_code) {
  core::SimulationState& state = context.state;
  ensureMagneticState(state, context.integrator_state, domain);
  const core::CosmologicalStepFactors& step = context.timeline_step;
  mhd::MhdStageContext interval{};
  interval.dt_code = context.integrator_state.dt_time_code;
  interval.time_code = context.integrator_state.current_time_code;
  interval.comoving_coordinates = step.cosmological;
  if (step.cosmological) {
    interval.scale_factor = step.scale_factor_begin;
    interval.scale_factor_end = step.scale_factor_end;
    interval.hubble_rate_code = step.hubble_begin_code;
    interval.hubble_rate_end_code = step.hubble_end_code;
    // REF 4.3/8.5 expansion ownership (pass-02 R-EXP): the MHD ideal stage is
    // the single owner of -H m_c and -H(2K_c + 3p_c + M_c); the gas KDK kick
    // of magnetized runs applies no drag or adiabatic scaling.
    interval.apply_expansion_source = true;
  } else {
    interval.scale_factor = 1.0;
    interval.scale_factor_end = 1.0;
    interval.hubble_rate_code = 0.0;
    interval.hubble_rate_end_code = 0.0;
    interval.apply_expansion_source = false;
  }
  loadGas(state, interval.scale_factor);
  ++m_telemetry.attempted_intervals;
  try {
    m_last = m_stepper->advance(interval, advanceOptions());
  } catch (const mhd::MhdError& error) {
    ++m_telemetry.rejected_intervals;
    throw core::StageTransactionRejected(rejectionKindFor(error.code()),
                                         std::string("MHD interval: ") + error.what());
  }
  m_telemetry.ideal_stage_evaluations_attempted += m_last.ideal_stage_evaluations_attempted;
  m_telemetry.mpi_bytes_sent += m_last.mpi_bytes_sent;
  m_telemetry.mpi_bytes_received += m_last.mpi_bytes_received;
  m_telemetry.exchange_rounds += m_last.exchange_rounds;
  m_telemetry.exchange_phases += m_last.exchange_phases;
  if (!m_last.accepted) {
    ++m_telemetry.rejected_intervals;
    std::ostringstream stream;
    stream << "MHD interval rejected (retryable, nothing committed): " << m_last.last_rejection
           << "; suggested dt " << m_last.dt_suggested_code;
    throw core::StageTransactionRejected(core::StageRejectionKind::kNumericalRetry, stream.str());
  }
  // Pass 04 (M03-12): accepted physical-boundary budget of this interval.
  // Lane 0 is mass, lane 4 the ideal TOTAL (gas + magnetic) energy flux of
  // the HLLD/HLLE face solver. Gas: the record integral I = sum_stage w F / a updates a cell by
  // -/+ I / dx, so the comoving amount leaving through a face is +/- I * A
  // (upper side positive outward, lower side negative). Non-ideal Poynting
  // energy uses the same convention with the energy-record integral.
  // Rank-local sums are reduced collectively (zero-patch ranks contribute 0).
  {
    std::array<double, 3> boundary{0.0, 0.0, 0.0};  // mass, gas energy, Poynting energy out
    const mhd::MhdLevelAdvanceOptions options = advanceOptions();
    if (options.produce_boundary_flux_integrals) {
      for (std::size_t index = 0; index < m_patches.size(); ++index) {
        const mhd::MhdPatchBoundary& sides = m_stepper->patchBoundary(index);
        const mhd::MhdBoundaryFaceFluxes& gas = m_stepper->boundaryGasFluxIntegral(index);
        const mhd::MhdBoundaryFaceEnergyFluxes& energy =
            m_stepper->boundaryEnergyFluxIntegral(index);
        const mhd::MhdPatchGeometry& g = m_patches[index].geometry;
        const std::array<std::size_t, 3> n{g.nx, g.ny, g.nz};
        const std::array<double, 3> area{g.faceAreaX(), g.faceAreaY(), g.faceAreaZ()};
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          const std::size_t t1 = (axis + 1U) % 3U;
          const std::size_t t2 = (axis + 2U) % 3U;
          for (std::size_t side = 0; side < 2U; ++side) {
            const mhd::MhdBoundaryKind kind = side == 0U ? sides.lower[axis] : sides.upper[axis];
            if (kind != mhd::MhdBoundaryKind::kOpen &&
                kind != mhd::MhdBoundaryKind::kReflectingConductingWall) {
              continue;
            }
            const double outward = side == 0U ? -1.0 : 1.0;
            for (std::size_t b = 0; b < n[t2]; ++b) {
              for (std::size_t a = 0; a < n[t1]; ++a) {
                const double* record = gas.at(axis, side, a, b);
                boundary[0] += outward * record[0] * area[axis];
                boundary[1] += outward * record[4] * area[axis];
                if (m_rate != nullptr) {
                  boundary[2] += outward * energy.at(axis, side, a, b) * area[axis];
                }
              }
            }
          }
        }
      }
    }
    if (m_transport.size() > 1) {
      m_transport.allreduceSumDoubles(std::span<double>(boundary.data(), boundary.size()));
    }
    m_telemetry.boundary_mass_out_comoving += boundary[0];
    m_telemetry.boundary_gas_energy_out_comoving += boundary[1];
    m_telemetry.boundary_poynting_energy_out_comoving += boundary[2];
    m_telemetry.last_boundary_mass_out_comoving = boundary[0];
    m_telemetry.last_boundary_energy_out_comoving = boundary[1] + boundary[2];
  }
  std::exception_ptr failure;
  try {
    hydro::GasThermodynamicRefreshPolicy policy{};
    policy.adiabatic_index = m_policy.gamma_adiabatic;
    policy.density_floor_code = 0.0;
    policy.pressure_floor_code = 0.0;
    policy.scale_factor = interval.scale_factor_end;
    policy.redshift = 1.0 / interval.scale_factor_end - 1.0;
    policy.closure = nullptr;
    policy.temperature_k_per_internal_energy_code = temperature_k_per_internal_energy_code;
    for (LocalPatch& patch : m_patches) {
      mhd::MhdAdapterPatchBinding adapter{patch.patch_id, patch.geometry_generation, patch.level,
                                          patch.geometry};
      const mhd::MhdFaceFluxView view =
          mhd::makeMagneticFaceView(state.magnetic_faces, adapter, patch.face_record);
      mhd::storeMhdGasToState(state, patch.geometry, patch.rows, patch.gas, view,
                              interval.scale_factor_end, m_policy.gamma_adiabatic);
      for (const std::uint32_t row : patch.rows) {
        (void)hydro::refreshGasThermodynamicCell(state, row, policy);
      }
    }
    state.magnetic_faces.divergence_max_normalized = m_last.divergence_max_normalized;
  } catch (...) {
    failure = std::current_exception();
  }
  m_mpi.rethrowCollectivePreparationFailure(failure, "MHD accepted-state store");
  ++m_telemetry.accepted_intervals;
  m_telemetry.ideal_stage_evaluations_accepted += m_last.ideal_stage_evaluations_accepted;
  m_telemetry.hlld_fallback_faces += m_last.ideal.hlle_fallback_faces;
  m_telemetry.first_order_faces += m_last.ideal.first_order_faces;
  m_telemetry.hall_subcycles += m_last.transport.hall_subcycle_count;
  m_telemetry.rkl_stages += m_last.transport.rkl_stage_count;
  m_telemetry.rkl_cycles += m_last.transport.rkl_cycle_count;
  m_telemetry.ghost_prolongation_slope_fallbacks += m_last.ghost_prolongation_slope_fallbacks;
  m_telemetry.scratch_capacity_bytes = m_last.scratch_capacity_bytes;
  m_telemetry.staging_capacity_bytes = ownedCapacityBytes() - m_last.scratch_capacity_bytes;
  m_telemetry.divergence_max_normalized = m_last.divergence_max_normalized;
  m_telemetry.last_dt_accepted_code = m_last.dt_accepted_code;
  m_telemetry.last_magnetic_energy = m_last.budget_end.energy_magnetic;
  m_telemetry.last_total_energy = m_last.budget_end.energy_total;
}

}  // namespace cosmosim::workflows::internal
