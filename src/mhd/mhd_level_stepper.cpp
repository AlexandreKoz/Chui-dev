// Synchronized, transactional, topology-aware MHD interval owner.
//
// CHUI-MHD-REF v1.0 sections 5.2, 6.1, 7.3-7.5, 8.5, 15 and CHUI-MHD-PERF v1.0
// sections 5-9, 13. UNVALIDATED source (pass 02; never compiled or run in
// the source-only campaign). See
// docs/architecture/chui_mhd_implementation_pass02_status.md.

#include "cosmosim/mhd/mhd_level_stepper.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <sstream>
#include <tuple>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/openmp_runtime.hpp"

#if COSMOSIM_HAVE_OPENMP
#include <omp.h>
#endif

namespace cosmosim::mhd {
namespace {

constexpr std::uint64_t k_payload_magic = 0x434855494D484432ULL;  // "CHUIMHD2"
constexpr std::size_t k_header_words = 4U;
constexpr std::size_t k_max_batch_slots = 64U;

MhdFailureContext makeContext(std::string_view owner, std::string detail) {
  MhdFailureContext context{};
  context.owner = std::string(owner);
  context.detail = std::move(detail);
  return context;
}

// Failures that a shorter interval can legitimately cure. Everything else is
// a contract violation that must not masquerade as a numerical retry.
[[nodiscard]] bool retryable(MhdErrorCode code) noexcept {
  return code == MhdErrorCode::kNonPositiveState ||
      code == MhdErrorCode::kStiffnessBudgetExceeded;
}

// Collective failure word: non-retryable codes dominate retryable ones so all
// ranks take the same rollback/retry decision.
[[nodiscard]] std::uint64_t encodeFailure(MhdErrorCode code) noexcept {
  const std::uint64_t priority = retryable(code) ? 1U : 2U;
  return (priority << 8U) | (static_cast<std::uint64_t>(code) + 1U);
}

[[nodiscard]] MhdErrorCode decodeFailure(std::uint64_t word) noexcept {
  const std::uint64_t code = (word & 0xFFU);
  if (code == 0U || code > static_cast<std::uint64_t>(MhdErrorCode::kMpiEpochMismatch) + 1U) {
    return MhdErrorCode::kInvalidField;
  }
  return static_cast<MhdErrorCode>(code - 1U);
}

[[nodiscard]] std::size_t faceIndex(const MhdPatchGeometry& geometry, std::size_t axis,
                                    std::int64_t i, std::int64_t j, std::int64_t k) noexcept {
  const auto ui = static_cast<std::size_t>(i);
  const auto uj = static_cast<std::size_t>(j);
  const auto uk = static_cast<std::size_t>(k);
  return axis == 0U ? geometry.faceIndexX(ui, uj, uk)
                    : (axis == 1U ? geometry.faceIndexY(ui, uj, uk) : geometry.faceIndexZ(ui, uj, uk));
}

[[nodiscard]] std::size_t edgeIndex(const MhdPatchGeometry& geometry, std::size_t axis,
                                    std::int64_t i, std::int64_t j, std::int64_t k) noexcept {
  const auto ui = static_cast<std::size_t>(i);
  const auto uj = static_cast<std::size_t>(j);
  const auto uk = static_cast<std::size_t>(k);
  return axis == 0U ? geometry.edgeIndexX(ui, uj, uk)
                    : (axis == 1U ? geometry.edgeIndexY(ui, uj, uk) : geometry.edgeIndexZ(ui, uj, uk));
}

[[nodiscard]] std::span<double> emfLane(MhdEdgeEmfStorage& emf, std::size_t axis) noexcept {
  return axis == 0U ? emf.emfX() : (axis == 1U ? emf.emfY() : emf.emfZ());
}

[[nodiscard]] std::span<const double> emfLane(const MhdEdgeEmfStorage& emf,
                                              std::size_t axis) noexcept {
  return axis == 0U ? emf.emfX() : (axis == 1U ? emf.emfY() : emf.emfZ());
}

[[nodiscard]] double faceArea(const MhdPatchGeometry& geometry, std::size_t axis) noexcept {
  return axis == 0U ? geometry.faceAreaX()
                    : (axis == 1U ? geometry.faceAreaY() : geometry.faceAreaZ());
}

[[nodiscard]] std::size_t extentAlong(const MhdPatchGeometry& geometry, std::size_t axis) noexcept {
  return axis == 0U ? geometry.nx : (axis == 1U ? geometry.ny : geometry.nz);
}

void copyLane(std::span<double> to, std::span<const double> from) noexcept {
  std::copy(from.begin(), from.end(), to.begin());
}

void copyEdgeField(MhdEdgeEmfStorage& to, const MhdEdgeEmfStorage& from) noexcept {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    copyLane(emfLane(to, axis), emfLane(from, axis));
  }
}

void copyRecords(MhdBoundaryFaceFluxes& to, const MhdBoundaryFaceFluxes& from) noexcept {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    for (std::size_t side = 0; side < 2U; ++side) {
      copyLane(to.values[axis][side].span(), from.values[axis][side].span());
    }
  }
}

void copyRecords(MhdBoundaryFaceEnergyFluxes& to, const MhdBoundaryFaceEnergyFluxes& from) noexcept {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    for (std::size_t side = 0; side < 2U; ++side) {
      copyLane(to.values[axis][side].span(), from.values[axis][side].span());
    }
  }
}

void zeroRecords(MhdBoundaryFaceFluxes& records) noexcept {
  for (auto& axis : records.values) {
    for (MhdDoubleArray& lane : axis) {
      lane.assignFill(0.0);
    }
  }
}

void zeroRecords(MhdBoundaryFaceEnergyFluxes& records) noexcept {
  for (auto& axis : records.values) {
    for (MhdDoubleArray& lane : axis) {
      lane.assignFill(0.0);
    }
  }
}

// target = alpha * base + beta * target + weight * add, lane-wise.
template <class Records>
void combineRecords(Records& target, const Records& base, double alpha, double beta,
                    const Records* add, double weight) noexcept {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    for (std::size_t side = 0; side < 2U; ++side) {
      std::span<double> t = target.values[axis][side].span();
      const std::span<const double> b = base.values[axis][side].span();
      for (std::size_t index = 0; index < t.size(); ++index) {
        double value = alpha * b[index] + beta * t[index];
        if (add != nullptr) {
          value += weight * add->values[axis][side][index];
        }
        t[index] = value;
      }
    }
  }
}

// Maximum oriented flux imbalance per cell (no allocation).
[[nodiscard]] double maxFluxImbalance(const MhdPatchGeometry& geometry, std::span<const double> phi_x,
                                      std::span<const double> phi_y,
                                      std::span<const double> phi_z) noexcept {
  double maximum = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const double imbalance =
            phi_x[geometry.faceIndexX(i + 1U, j, k)] - phi_x[geometry.faceIndexX(i, j, k)] +
            phi_y[geometry.faceIndexY(i, j + 1U, k)] - phi_y[geometry.faceIndexY(i, j, k)] +
            phi_z[geometry.faceIndexZ(i, j, k + 1U)] - phi_z[geometry.faceIndexZ(i, j, k)];
        if (!std::isfinite(imbalance)) {
          return std::numeric_limits<double>::infinity();
        }
        maximum = std::max(maximum, std::abs(imbalance));
      }
    }
  }
  return maximum;
}

[[nodiscard]] double minmod(double a, double b) noexcept {
  if (a * b <= 0.0) {
    return 0.0;
  }
  return a > 0.0 ? std::min(a, b) : std::max(a, b);
}

[[nodiscard]] std::array<double, 9> packPrimitive(const MhdPrimitiveState& w) noexcept {
  return {w.rho_comoving, w.vel_x_peculiar, w.vel_y_peculiar, w.vel_z_peculiar,
          w.pressure_comoving, w.bx_comoving, w.by_comoving, w.bz_comoving,
          w.metallicity_mass_fraction};
}

[[nodiscard]] MhdPrimitiveState unpackPrimitive(const double* v) noexcept {
  MhdPrimitiveState w{};
  w.rho_comoving = v[0];
  w.vel_x_peculiar = v[1];
  w.vel_y_peculiar = v[2];
  w.vel_z_peculiar = v[3];
  w.pressure_comoving = v[4];
  w.bx_comoving = v[5];
  w.by_comoving = v[6];
  w.bz_comoving = v[7];
  w.metallicity_mass_fraction = v[8];
  return w;
}

// RAII interval transaction: restores the rollback image unless committed.
template <class Restore>
class RollbackGuard {
 public:
  explicit RollbackGuard(Restore restore) noexcept : m_restore(restore) {}
  RollbackGuard(const RollbackGuard&) = delete;
  RollbackGuard& operator=(const RollbackGuard&) = delete;
  ~RollbackGuard() {
    if (m_active) {
      m_restore();
    }
  }
  void commit() noexcept { m_active = false; }

 private:
  Restore m_restore;
  bool m_active = true;
};

}  // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

MhdLevelStepper::MhdLevelStepper(MhdNumericsPolicy policy, MhdLevelDomain domain,
                                 std::vector<MhdLevelPatchBinding> patches,
                                 std::vector<MhdRemotePatchDescriptor> remote_patches,
                                 MhdLevelStepperConfig config)
    : m_policy(policy), m_operator(policy), m_domain(domain), m_config(config) {
  m_transport = m_config.transport != nullptr ? m_config.transport : &m_serial_transport;
  if (m_config.batch_patches == 0U || m_config.batch_patches > k_max_batch_slots) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_level_stepper", "batch_patches must lie in [1, 64]"));
  }
  if (m_config.comm_budget_bytes == 0U ||
      m_config.comm_budget_bytes > k_mhd_max_comm_budget_bytes_per_rank) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_level_stepper",
                               "communication budget must lie in (0, 128 MiB] (PERF 8.1)"));
  }
  m_local.reserve(patches.size());
  for (MhdLevelPatchBinding& binding : patches) {
    PatchRuntime runtime{};
    runtime.binding = binding;
    m_local.push_back(std::move(runtime));
  }
  buildTopology(std::move(remote_patches));
  const std::size_t slots =
      std::min<std::size_t>(m_config.batch_patches, std::max<std::size_t>(m_local.size(), 1U));
  m_slots.resize(m_local.empty() ? 0U : slots);
  const auto world = static_cast<std::size_t>(m_transport->size());
  m_send_values.resize(world);
  m_send_bytes.resize(world);
  m_recv_values.resize(world);
}

std::span<double> MhdLevelStepper::faceLane(PatchRuntime& patch, std::size_t axis) noexcept {
  return axis == 0U ? patch.binding.faces.phi_x
                    : (axis == 1U ? patch.binding.faces.phi_y : patch.binding.faces.phi_z);
}

std::span<const double> MhdLevelStepper::faceLane(const PatchRuntime& patch,
                                                  std::size_t axis) const noexcept {
  return axis == 0U ? patch.binding.faces.phi_x
                    : (axis == 1U ? patch.binding.faces.phi_y : patch.binding.faces.phi_z);
}

const MhdEdgeEmfStorage& MhdLevelStepper::edgeIntegral(std::size_t index) const {
  return m_local.at(index).integral;
}

const MhdBoundaryFaceFluxes& MhdLevelStepper::boundaryGasFluxIntegral(std::size_t index) const {
  return m_local.at(index).gas_integral;
}

const MhdBoundaryFaceEnergyFluxes& MhdLevelStepper::boundaryEnergyFluxIntegral(
    std::size_t index) const {
  return m_local.at(index).energy_integral;
}

// ---------------------------------------------------------------------------
// Collective error agreement and exchange framework
// ---------------------------------------------------------------------------

void MhdLevelStepper::agree(const std::exception_ptr& local_failure, std::string_view phase) {
  std::uint64_t local_word = 0U;
  std::string local_detail;
  if (local_failure) {
    try {
      std::rethrow_exception(local_failure);
    } catch (const MhdError& error) {
      local_word = encodeFailure(error.code());
      local_detail = error.what();
    } catch (const std::bad_alloc&) {
      local_word = encodeFailure(MhdErrorCode::kMemoryBudgetExceeded);
      local_detail = "ungoverned allocation failure (std::bad_alloc)";
    } catch (const std::exception& error) {
      local_word = encodeFailure(MhdErrorCode::kInvalidField);
      local_detail = std::string("foreign failure: ") + error.what();
    } catch (...) {
      local_word = encodeFailure(MhdErrorCode::kInvalidField);
      local_detail = "foreign failure of unknown type";
    }
  }
  const std::uint64_t global_word =
      m_transport->size() > 1 ? m_transport->allreduceMaxUint64(local_word) : local_word;
  if (global_word == 0U) {
    return;
  }
  const MhdErrorCode code = decodeFailure(global_word);
  if (local_failure && local_word == global_word) {
    try {
      std::rethrow_exception(local_failure);
    } catch (const MhdError&) {
      throw;
    } catch (...) {
      // Foreign or allocation failure: rethrown as the agreed MHD code so all
      // ranks take the same transactional decision.
    }
  }
  MhdFailureContext context = makeContext(
      "mhd_level_stepper",
      std::string("collective failure during ") + std::string(phase) +
          (local_detail.empty() ? std::string(" (raised on another rank)")
                                : std::string(": ") + local_detail));
  context.rank = m_transport->rank();
  throw MhdError(code, std::move(context));
}

template <class Fn>
void MhdLevelStepper::collective(std::string_view phase, Fn&& fn) {
  std::exception_ptr failure;
  try {
    fn();
  } catch (...) {
    failure = std::current_exception();
  }
  agree(failure, phase);
}

template <class Packer, class Applier>
void MhdLevelStepper::runItems(std::uint64_t phase_id, std::span<const Item> items, Packer&& pack,
                               Applier&& apply, bool remote_only) {
  const int world = m_transport->size();
  const int me = m_transport->rank();
  const std::uint64_t phase_word = (phase_id << 32U) ^ (m_phase_counter++ & 0xFFFFFFFFULL);
  if (world > 1) {
    std::exception_ptr failure;
    try {
      for (std::vector<double>& values : m_send_values) {
        values.clear();
      }
      for (const Item& item : items) {
        const TopoPatch& src = m_topo[item.src];
        const TopoPatch& dst = m_topo[item.dst];
        if (src.local == k_npos || dst.owner_rank == me) {
          continue;
        }
        std::vector<double>& out = m_send_values[static_cast<std::size_t>(dst.owner_rank)];
        const std::size_t before = out.size();
        pack(item, out);
        if (out.size() - before != itemValueCount(item)) {
          throw MhdError(MhdErrorCode::kMpiEpochMismatch,
                         makeContext("mhd_level_exchange", "packed value count mismatch"));
        }
      }
      for (std::size_t rank = 0; rank < m_send_values.size(); ++rank) {
        std::vector<std::uint8_t>& bytes = m_send_bytes[rank];
        bytes.clear();
        const std::vector<double>& values = m_send_values[rank];
        if (values.empty() || static_cast<int>(rank) == me) {
          continue;
        }
        const std::array<std::uint64_t, k_header_words> header{
            k_payload_magic, phase_word, m_config.hierarchy_epoch,
            static_cast<std::uint64_t>(values.size())};
        bytes.resize(sizeof(header) + values.size() * sizeof(double));
        std::memcpy(bytes.data(), header.data(), sizeof(header));
        std::memcpy(bytes.data() + sizeof(header), values.data(), values.size() * sizeof(double));
      }
    } catch (...) {
      failure = std::current_exception();
    }
    agree(failure, "exchange packing");
    const std::vector<std::vector<std::uint8_t>> received = m_transport->exchange(m_send_bytes);
    failure = nullptr;
    try {
      if (received.size() != static_cast<std::size_t>(world)) {
        throw MhdError(MhdErrorCode::kMpiEpochMismatch,
                       makeContext("mhd_level_exchange", "transport returned a wrong rank count"));
      }
      for (std::size_t rank = 0; rank < received.size(); ++rank) {
        m_bytes_sent += m_send_bytes[rank].size();
        m_bytes_received += received[rank].size();
        std::vector<double>& values = m_recv_values[rank];
        values.clear();
        const std::vector<std::uint8_t>& bytes = received[rank];
        if (bytes.empty()) {
          continue;
        }
        std::array<std::uint64_t, k_header_words> header{};
        if (bytes.size() < sizeof(header)) {
          throw MhdError(MhdErrorCode::kMpiEpochMismatch,
                         makeContext("mhd_level_exchange", "truncated payload header"));
        }
        std::memcpy(header.data(), bytes.data(), sizeof(header));
        if (header[0] != k_payload_magic || header[1] != phase_word ||
            header[2] != m_config.hierarchy_epoch ||
            bytes.size() != sizeof(header) + header[3] * sizeof(double)) {
          MhdFailureContext context = makeContext(
              "mhd_level_exchange",
              "payload magic, phase, hierarchy epoch or count disagree with the receiver plan");
          context.rank = static_cast<std::int64_t>(rank);
          throw MhdError(MhdErrorCode::kMpiEpochMismatch, std::move(context));
        }
        values.resize(static_cast<std::size_t>(header[3]));
        std::memcpy(values.data(), bytes.data() + sizeof(header), values.size() * sizeof(double));
      }
    } catch (...) {
      failure = std::current_exception();
    }
    agree(failure, "exchange receive");
  }
  std::exception_ptr failure;
  try {
    std::vector<std::size_t> cursor(static_cast<std::size_t>(world), 0U);
    for (const Item& item : items) {
      const TopoPatch& src = m_topo[item.src];
      const TopoPatch& dst = m_topo[item.dst];
      if (dst.local == k_npos) {
        continue;
      }
      const std::size_t count = itemValueCount(item);
      if (src.local != k_npos) {
        if (remote_only) {
          continue;
        }
        m_exchange_scratch.clear();
        pack(item, m_exchange_scratch);
        apply(item, m_exchange_scratch.data());
        continue;
      }
      const auto rank = static_cast<std::size_t>(src.owner_rank);
      const std::vector<double>& values = m_recv_values[rank];
      if (cursor[rank] + count > values.size()) {
        throw MhdError(MhdErrorCode::kMpiEpochMismatch,
                       makeContext("mhd_level_exchange",
                                   "received fewer values than the receiver plan requires"));
      }
      apply(item, values.data() + cursor[rank]);
      cursor[rank] += count;
    }
    for (std::size_t rank = 0; rank < cursor.size(); ++rank) {
      if (static_cast<int>(rank) != me && cursor[rank] != m_recv_values[rank].size()) {
        throw MhdError(MhdErrorCode::kMpiEpochMismatch,
                       makeContext("mhd_level_exchange",
                                   "received more values than the receiver plan consumes"));
      }
    }
  } catch (...) {
    failure = std::current_exception();
  }
  agree(failure, "exchange apply");
  ++m_exchange_phases;
}

// Values of one item computed from its (local) source patch canonical or
// stage state, appended in (k, j, i) order of the destination box.
void MhdLevelStepper::computeItemValues(const Item& item, double scale_factor,
                                        std::vector<double>& out) const {
  const PatchRuntime& source = m_local[m_topo[item.src].local];
  const MhdPatchGeometry& g = source.binding.geometry;
  const MhdGasStateSoa& gas = *source.binding.gas;
  const std::span<const double> phi_x = source.binding.faces.phi_x;
  const std::span<const double> phi_y = source.binding.faces.phi_y;
  const std::span<const double> phi_z = source.binding.faces.phi_z;
  const std::array<std::span<const double>, 3> phi{phi_x, phi_y, phi_z};
  const double inverse_sqrt_a = 1.0 / std::sqrt(scale_factor);
  const std::size_t axis = item.axis;
  const std::size_t t1 = (axis + 1U) % 3U;
  const std::size_t t2 = (axis + 2U) % 3U;
  const auto fail = [&](std::string detail) {
    MhdFailureContext context = makeContext("mhd_level_halo", std::move(detail));
    context.patch_id = source.binding.patch_id;
    context.geometry_generation = source.binding.geometry_generation;
    context.scale_factor = scale_factor;
    throw MhdError(MhdErrorCode::kNonPositiveState, std::move(context));
  };
  const auto direct_cell = [&](std::int64_t i, std::int64_t j, std::int64_t k) {
    const auto ui = static_cast<std::size_t>(i);
    const auto uj = static_cast<std::size_t>(j);
    const auto uk = static_cast<std::size_t>(k);
    MhdPrimitiveState w{};
    if (!m_operator.recoverCellPrimitive(
            gas.loadGas(g.cellIndex(ui, uj, uk)),
            MhdStepOperator::cellFieldFromFaces(g, phi_x, phi_y, phi_z, ui, uj, uk), scale_factor,
            w)) {
      fail("inadmissible source cell for a ghost/sample value (no floor applied)");
    }
    const std::array<double, 9> v = packPrimitive(w);
    out.insert(out.end(), v.begin(), v.end());
  };
  const auto restricted_cell = [&](std::int64_t i0, std::int64_t j0, std::int64_t k0) {
    MhdConservedState sum{};
    for (std::int64_t dk = 0; dk < 2; ++dk) {
      for (std::int64_t dj = 0; dj < 2; ++dj) {
        for (std::int64_t di = 0; di < 2; ++di) {
          sum += gas.loadGas(g.cellIndex(static_cast<std::size_t>(i0 + di),
                                         static_cast<std::size_t>(j0 + dj),
                                         static_cast<std::size_t>(k0 + dk)));
        }
      }
    }
    MhdConservedState mean{};
    mean.rho_comoving = sum.rho_comoving / 8.0;
    mean.momentum_density_x_comoving = sum.momentum_density_x_comoving / 8.0;
    mean.momentum_density_y_comoving = sum.momentum_density_y_comoving / 8.0;
    mean.momentum_density_z_comoving = sum.momentum_density_z_comoving / 8.0;
    mean.total_energy_density_comoving = sum.total_energy_density_comoving / 8.0;
    mean.metal_mass_density_comoving = sum.metal_mass_density_comoving / 8.0;
    // Coarse face fluxes = sums of the four fine faces; B_c = mean of the two.
    const std::array<std::int64_t, 3> base{i0, j0, k0};
    std::array<double, 3> b{};
    for (std::size_t n = 0; n < 3U; ++n) {
      const std::size_t a1 = (n + 1U) % 3U;
      const std::size_t a2 = (n + 2U) % 3U;
      double lower = 0.0;
      double upper = 0.0;
      for (std::int64_t d2 = 0; d2 < 2; ++d2) {
        for (std::int64_t d1 = 0; d1 < 2; ++d1) {
          std::array<std::int64_t, 3> f = base;
          f[a1] += d1;
          f[a2] += d2;
          lower += phi[n][faceIndex(g, n, f[0], f[1], f[2])];
          f[n] += 2;
          upper += phi[n][faceIndex(g, n, f[0], f[1], f[2])];
        }
      }
      b[n] = 0.5 * (lower + upper) / (4.0 * faceArea(g, n));
    }
    MhdPrimitiveState w{};
    if (!m_operator.recoverCellPrimitive(mean, b, scale_factor, w)) {
      fail("inadmissible restricted fine state for a ghost/sample value");
    }
    const std::array<double, 9> v = packPrimitive(w);
    out.insert(out.end(), v.begin(), v.end());
  };
  const Box& box = item.dst_box;
  const auto scale = static_cast<std::int64_t>(item.scale);
  for (std::int64_t k = box.lo[2]; k < box.hi[2]; ++k) {
    for (std::int64_t j = box.lo[1]; j < box.hi[1]; ++j) {
      for (std::int64_t i = box.lo[0]; i < box.hi[0]; ++i) {
        const std::array<std::int64_t, 3> s{scale * i + item.src_offset[0],
                                            scale * j + item.src_offset[1],
                                            scale * k + item.src_offset[2]};
        switch (item.kind) {
          case ItemKind::kGhostCellDirect:
          case ItemKind::kSampleCellDirect:
            direct_cell(s[0], s[1], s[2]);
            break;
          case ItemKind::kGhostCellRestricted:
          case ItemKind::kSampleCellRestricted:
            restricted_cell(s[0], s[1], s[2]);
            break;
          case ItemKind::kGhostFaceDirect:
          case ItemKind::kSampleFaceDirect:
            out.push_back(phi[axis][faceIndex(g, axis, s[0], s[1], s[2])] /
                          faceArea(g, axis) * inverse_sqrt_a);
            break;
          case ItemKind::kGhostFaceRestricted: {
            double sum = 0.0;
            for (std::int64_t d2 = 0; d2 < 2; ++d2) {
              for (std::int64_t d1 = 0; d1 < 2; ++d1) {
                std::array<std::int64_t, 3> f = s;
                f[t1] += d1;
                f[t2] += d2;
                sum += phi[axis][faceIndex(g, axis, f[0], f[1], f[2])];
              }
            }
            out.push_back(sum / (4.0 * faceArea(g, axis)) * inverse_sqrt_a);
            break;
          }
          case ItemKind::kFaceOwner:
            out.push_back(phi[axis][faceIndex(g, axis, s[0], s[1], s[2])]);
            break;
          case ItemKind::kFaceFine: {
            double sum = 0.0;
            for (std::int64_t d2 = 0; d2 < 2; ++d2) {
              for (std::int64_t d1 = 0; d1 < 2; ++d1) {
                std::array<std::int64_t, 3> f = s;
                f[t1] += d1;
                f[t2] += d2;
                sum += phi[axis][faceIndex(g, axis, f[0], f[1], f[2])];
              }
            }
            out.push_back(sum);
            break;
          }
          case ItemKind::kEdgeOwner:
            out.push_back(emfLane(source.emf, axis)[edgeIndex(g, axis, s[0], s[1], s[2])]);
            break;
          case ItemKind::kEdgeFine: {
            std::array<std::int64_t, 3> second = s;
            second[axis] += 1;
            const std::span<const double> lane = emfLane(source.emf, axis);
            out.push_back(lane[edgeIndex(g, axis, s[0], s[1], s[2])] +
                          lane[edgeIndex(g, axis, second[0], second[1], second[2])]);
            break;
          }
          case ItemKind::kGasRecordOwner:
          case ItemKind::kGasRecordFine:
          case ItemKind::kEnergyRecordOwner:
          case ItemKind::kEnergyRecordFine: {
            const auto node = static_cast<std::size_t>(s[axis]);
            const std::size_t n_axis = extentAlong(g, axis);
            if (node != 0U && node != n_axis) {
              throw MhdError(MhdErrorCode::kStaleGeometry,
                             makeContext("mhd_level_exchange",
                                         "boundary record item does not lie on a patch boundary"));
            }
            const std::size_t side = node == 0U ? 0U : 1U;
            const bool gas_record = item.kind == ItemKind::kGasRecordOwner ||
                item.kind == ItemKind::kGasRecordFine;
            const bool fine = item.kind == ItemKind::kGasRecordFine ||
                item.kind == ItemKind::kEnergyRecordFine;
            const std::int64_t span = fine ? 2 : 1;
            std::array<double, MhdBoundaryFaceFluxes::k_values_per_face> sum{};
            double energy_sum = 0.0;
            for (std::int64_t d2 = 0; d2 < span; ++d2) {
              for (std::int64_t d1 = 0; d1 < span; ++d1) {
                const auto c1 = static_cast<std::size_t>(s[t1] + d1);
                const auto c2 = static_cast<std::size_t>(s[t2] + d2);
                if (gas_record) {
                  const double* record = source.records.at(axis, side, c1, c2);
                  for (std::size_t lane = 0; lane < sum.size(); ++lane) {
                    sum[lane] += record[lane];
                  }
                } else {
                  energy_sum += source.energy_records.at(axis, side, c1, c2);
                }
              }
            }
            const double inverse = fine ? 0.25 : 1.0;
            if (gas_record) {
              for (const double value : sum) {
                out.push_back(value * inverse);
              }
            } else {
              out.push_back(energy_sum * inverse);
            }
            break;
          }
        }
      }
    }
  }
}

void MhdLevelStepper::exchangeHalo(double scale_factor) {
  if (m_transport->size() <= 1) {
    return;  // every halo source is local and is read directly at batch time
  }
  runItems(
      1U, std::span<const Item>(m_halo_items),
      [&](const Item& item, std::vector<double>& out) {
        computeItemValues(item, scale_factor, out);
      },
      [&](const Item& item, const double* values) {
        PatchRuntime& patch = m_local[m_topo[item.dst].local];
        const std::size_t count = itemValueCount(item);
        if (item.import_offset + count > patch.halo_import.size()) {
          throw MhdError(MhdErrorCode::kMpiEpochMismatch,
                         makeContext("mhd_level_halo", "import offset outside the admitted buffer"));
        }
        std::copy(values, values + count, patch.halo_import.data() + item.import_offset);
      },
      true);
}

// ---------------------------------------------------------------------------
// Batch processing
// ---------------------------------------------------------------------------

template <class Fn>
void MhdLevelStepper::forEachLocalPatchBatched(Fn&& fn) {
  const std::size_t count = m_local.size();
  const std::size_t slots = m_slots.size();
  if (count == 0U) {
    return;
  }
  for (std::size_t begin = 0; begin < count; begin += slots) {
    const std::size_t chunk = std::min(slots, count - begin);
    std::array<std::exception_ptr, k_max_batch_slots> failures{};
    const auto chunk_signed = static_cast<std::ptrdiff_t>(chunk);
#if COSMOSIM_HAVE_OPENMP
    const int threads = std::max(
        1, std::min(core::configuredOpenMpThreadCount(), static_cast<int>(chunk)));
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (std::ptrdiff_t offset = 0; offset < chunk_signed; ++offset) {
      const auto slot = static_cast<std::size_t>(offset);
      try {
        fn(begin + slot, m_slots[slot]);
      } catch (...) {
        failures[slot] = std::current_exception();
      }
    }
    for (std::size_t slot = 0; slot < chunk; ++slot) {
      if (failures[slot]) {
        std::rethrow_exception(failures[slot]);
      }
    }
  }
}

void MhdLevelStepper::fillExternalGhosts(std::size_t patch_index, BatchSlot& slot,
                                         double scale_factor) {
  PatchRuntime& patch = m_local[patch_index];
  MhdPrimitivePatch& primitives = slot.workspace.primitives;
  MhdStaggeredFieldPatch& staggered = slot.workspace.staggered;
  const HaloPlan& plan = patch.halo;
  for (const Item& item : plan.items) {
    const double* values = nullptr;
    if (m_topo[item.src].local != k_npos) {
      slot.values.clear();
      computeItemValues(item, scale_factor, slot.values);
      values = slot.values.data();
    } else {
      values = patch.halo_import.data() + item.import_offset;
    }
    switch (item.kind) {
      case ItemKind::kGhostCellDirect:
      case ItemKind::kGhostCellRestricted: {
        std::size_t cursor = 0U;
        for (std::int64_t k = item.dst_box.lo[2]; k < item.dst_box.hi[2]; ++k) {
          for (std::int64_t j = item.dst_box.lo[1]; j < item.dst_box.hi[1]; ++j) {
            for (std::int64_t i = item.dst_box.lo[0]; i < item.dst_box.hi[0]; ++i) {
              primitives.store({static_cast<std::ptrdiff_t>(i), static_cast<std::ptrdiff_t>(j),
                                static_cast<std::ptrdiff_t>(k)},
                               unpackPrimitive(values + cursor));
              cursor += 9U;
            }
          }
        }
        break;
      }
      case ItemKind::kGhostFaceDirect:
      case ItemKind::kGhostFaceRestricted: {
        std::size_t cursor = 0U;
        const MhdAxis normal = static_cast<MhdAxis>(item.axis);
        for (std::int64_t k = item.dst_box.lo[2]; k < item.dst_box.hi[2]; ++k) {
          for (std::int64_t j = item.dst_box.lo[1]; j < item.dst_box.hi[1]; ++j) {
            for (std::int64_t i = item.dst_box.lo[0]; i < item.dst_box.hi[0]; ++i) {
              staggered.at(normal, {static_cast<std::ptrdiff_t>(i), static_cast<std::ptrdiff_t>(j),
                                    static_cast<std::ptrdiff_t>(k)}) = values[cursor++];
            }
          }
        }
        break;
      }
      case ItemKind::kSampleCellDirect:
      case ItemKind::kSampleCellRestricted:
        std::copy(values, values + 9, slot.cell_samples.data() + 9U * item.slot);
        break;
      case ItemKind::kSampleFaceDirect:
        slot.face_samples[item.slot] = values[0];
        break;
      default:
        throw MhdError(MhdErrorCode::kStaleGeometry,
                       makeContext("mhd_level_halo", "non-halo item in a halo plan"));
    }
  }
  // Coarse -> fine ghost cells: minmod-limited linear primitives (each axis
  // contributes at most a quarter of a one-sided difference, so density and
  // pressure stay positive by construction); a missing neighbour sample
  // (physical boundary or a two-level jump) gives a counted zero slope.
  for (const ProlongCell& recipe : plan.prolong_cells) {
    const double* center = slot.cell_samples.data() + 9U * static_cast<std::size_t>(recipe.center);
    std::array<double, 9> value{};
    for (std::size_t lane = 0; lane < 9U; ++lane) {
      value[lane] = center[lane];
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const std::int32_t minus = recipe.neighbour[2U * axis];
      const std::int32_t plus = recipe.neighbour[2U * axis + 1U];
      if (minus < 0 || plus < 0) {
        ++slot.slope_fallbacks;
        continue;
      }
      const double* lower = slot.cell_samples.data() + 9U * static_cast<std::size_t>(minus);
      const double* upper = slot.cell_samples.data() + 9U * static_cast<std::size_t>(plus);
      const double offset = 0.25 * static_cast<double>(recipe.half[axis]);
      for (std::size_t lane = 0; lane < 9U; ++lane) {
        value[lane] += offset * minmod(center[lane] - lower[lane], upper[lane] - center[lane]);
      }
    }
    MhdPrimitiveState state = unpackPrimitive(value.data());
    if (!(state.rho_comoving > 0.0) || !(state.pressure_comoving > 0.0) ||
        !std::isfinite(state.rho_comoving) || !std::isfinite(state.pressure_comoving)) {
      state = unpackPrimitive(center);
      ++slot.slope_fallbacks;
    }
    primitives.store(recipe.dst, state);
  }
  for (const ProlongFace& recipe : plan.prolong_faces) {
    const double first = slot.face_samples[static_cast<std::size_t>(recipe.first)];
    const double value = recipe.second < 0
        ? first
        : 0.5 * (first + slot.face_samples[static_cast<std::size_t>(recipe.second)]);
    staggered.at(static_cast<MhdAxis>(recipe.axis), recipe.dst) = value;
  }
}

void MhdLevelStepper::prepareSlot(std::size_t patch_index, BatchSlot& slot, double scale_factor) {
  PatchRuntime& patch = m_local[patch_index];
  const MhdFaceFluxView view{&patch.binding.geometry,
                             patch.binding.faces.phi_x,
                             patch.binding.faces.phi_y,
                             patch.binding.faces.phi_z,
                             patch.binding.patch_id,
                             patch.binding.geometry_generation,
                             patch.binding.level};
  const std::size_t bad = m_operator.recoverPrimitives(patch.binding.geometry, *patch.binding.gas,
                                                       view, scale_factor, slot.workspace);
  if (bad != 0U) {
    MhdFailureContext failure = makeContext(
        "mhd_level_stepper",
        std::to_string(bad) + " inadmissible cell(s) at stage preparation (no floor applied)");
    failure.patch_id = patch.binding.patch_id;
    failure.geometry_generation = patch.binding.geometry_generation;
    failure.scale_factor = scale_factor;
    throw MhdError(MhdErrorCode::kNonPositiveState, std::move(failure));
  }
  fillExternalGhosts(patch_index, slot, scale_factor);
  MhdStepOperator::fillPhysicalGhosts(patch.boundary, slot.workspace);
}

// ---------------------------------------------------------------------------
// Owner reconciliation and coarse <- fine overrides
// ---------------------------------------------------------------------------

void MhdLevelStepper::applyItemsWrite(std::span<const Item> items, std::uint64_t phase_id) {
  const auto apply = [&](const Item& item, const double* values) {
    PatchRuntime& patch = m_local[m_topo[item.dst].local];
    const MhdPatchGeometry& g = patch.binding.geometry;
    const std::size_t axis = item.axis;
    const std::size_t t1 = (axis + 1U) % 3U;
    const std::size_t t2 = (axis + 2U) % 3U;
    std::size_t cursor = 0U;
    for (std::int64_t k = item.dst_box.lo[2]; k < item.dst_box.hi[2]; ++k) {
      for (std::int64_t j = item.dst_box.lo[1]; j < item.dst_box.hi[1]; ++j) {
        for (std::int64_t i = item.dst_box.lo[0]; i < item.dst_box.hi[0]; ++i) {
          switch (item.kind) {
            case ItemKind::kEdgeOwner:
            case ItemKind::kEdgeFine:
              emfLane(patch.emf, axis)[edgeIndex(g, axis, i, j, k)] = values[cursor++];
              break;
            case ItemKind::kGasRecordOwner:
            case ItemKind::kGasRecordFine:
            case ItemKind::kEnergyRecordOwner:
            case ItemKind::kEnergyRecordFine: {
              const std::array<std::int64_t, 3> d{i, j, k};
              const auto node = static_cast<std::size_t>(d[axis]);
              if (node != 0U && node != extentAlong(g, axis)) {
                throw MhdError(MhdErrorCode::kStaleGeometry,
                               makeContext("mhd_level_exchange",
                                           "record destination is not a patch boundary face"));
              }
              const std::size_t side = node == 0U ? 0U : 1U;
              const auto c1 = static_cast<std::size_t>(d[t1]);
              const auto c2 = static_cast<std::size_t>(d[t2]);
              if (item.kind == ItemKind::kGasRecordOwner || item.kind == ItemKind::kGasRecordFine) {
                double* record = patch.records.at(axis, side, c1, c2);
                std::copy(values + cursor, values + cursor + MhdBoundaryFaceFluxes::k_values_per_face,
                          record);
                cursor += MhdBoundaryFaceFluxes::k_values_per_face;
              } else {
                patch.energy_records.at(axis, side, c1, c2) = values[cursor++];
              }
              break;
            }
            case ItemKind::kFaceOwner:
            case ItemKind::kFaceFine:
              faceLane(patch, axis)[faceIndex(g, axis, i, j, k)] = values[cursor++];
              break;
            default:
              throw MhdError(MhdErrorCode::kStaleGeometry,
                             makeContext("mhd_level_exchange", "unexpected item kind"));
          }
        }
      }
    }
  };
  runItems(
      phase_id, items,
      [&](const Item& item, std::vector<double>& out) { computeItemValues(item, 1.0, out); },
      apply, false);
}

void MhdLevelStepper::reconcileEdges() {
  collective("edges.self", [&]() {
    for (const Item& item : m_edge_self_items) {
      PatchRuntime& patch = m_local[m_topo[item.dst].local];
      const MhdPatchGeometry& g = patch.binding.geometry;
      std::span<double> lane = emfLane(patch.emf, item.axis);
      for (std::int64_t k = item.dst_box.lo[2]; k < item.dst_box.hi[2]; ++k) {
        for (std::int64_t j = item.dst_box.lo[1]; j < item.dst_box.hi[1]; ++j) {
          for (std::int64_t i = item.dst_box.lo[0]; i < item.dst_box.hi[0]; ++i) {
            lane[edgeIndex(g, item.axis, i, j, k)] =
                lane[edgeIndex(g, item.axis, i + item.src_offset[0], j + item.src_offset[1],
                               k + item.src_offset[2])];
          }
        }
      }
    }
  });
  applyItemsWrite(m_edge_cross_items, 10U);
  // Coarse <- fine, finest coarse level first (a coarse edge may itself be the
  // fine source of the next coarser override).
  for (std::size_t slot = m_edge_fine_items.size(); slot-- > 0U;) {
    applyItemsWrite(m_edge_fine_items[slot], 11U);
  }
}

void MhdLevelStepper::reconcileGasRecords() {
  collective("records.self", [&]() {
    for (const Item& item : m_record_self_items) {
      PatchRuntime& patch = m_local[m_topo[item.dst].local];
      const std::size_t axis = item.axis;
      const std::size_t t1 = (axis + 1U) % 3U;
      const std::size_t t2 = (axis + 2U) % 3U;
      for (std::int64_t k = item.dst_box.lo[2]; k < item.dst_box.hi[2]; ++k) {
        for (std::int64_t j = item.dst_box.lo[1]; j < item.dst_box.hi[1]; ++j) {
          for (std::int64_t i = item.dst_box.lo[0]; i < item.dst_box.hi[0]; ++i) {
            const std::array<std::int64_t, 3> d{i, j, k};
            const std::array<std::int64_t, 3> s{i + item.src_offset[0], j + item.src_offset[1],
                                                k + item.src_offset[2]};
            const std::size_t dst_side = d[axis] == 0 ? 0U : 1U;
            const std::size_t src_side = s[axis] == 0 ? 0U : 1U;
            const double* from = patch.records.at(axis, src_side, static_cast<std::size_t>(s[t1]),
                                                  static_cast<std::size_t>(s[t2]));
            double* to = patch.records.at(axis, dst_side, static_cast<std::size_t>(d[t1]),
                                          static_cast<std::size_t>(d[t2]));
            std::copy(from, from + MhdBoundaryFaceFluxes::k_values_per_face, to);
          }
        }
      }
    }
  });
  applyItemsWrite(m_record_cross_items, 20U);
  applyItemsWrite(m_record_fine_items, 21U);
}

void MhdLevelStepper::reconcileEnergyRecords() {
  collective("energy_records.self", [&]() {
    for (const Item& item : m_energy_self_items) {
      PatchRuntime& patch = m_local[m_topo[item.dst].local];
      const std::size_t axis = item.axis;
      const std::size_t t1 = (axis + 1U) % 3U;
      const std::size_t t2 = (axis + 2U) % 3U;
      for (std::int64_t k = item.dst_box.lo[2]; k < item.dst_box.hi[2]; ++k) {
        for (std::int64_t j = item.dst_box.lo[1]; j < item.dst_box.hi[1]; ++j) {
          for (std::int64_t i = item.dst_box.lo[0]; i < item.dst_box.hi[0]; ++i) {
            const std::array<std::int64_t, 3> d{i, j, k};
            const std::array<std::int64_t, 3> s{i + item.src_offset[0], j + item.src_offset[1],
                                                k + item.src_offset[2]};
            patch.energy_records.at(axis, d[axis] == 0 ? 0U : 1U, static_cast<std::size_t>(d[t1]),
                                    static_cast<std::size_t>(d[t2])) =
                patch.energy_records.at(axis, s[axis] == 0 ? 0U : 1U,
                                        static_cast<std::size_t>(s[t1]),
                                        static_cast<std::size_t>(s[t2]));
          }
        }
      }
    }
  });
  applyItemsWrite(m_energy_cross_items, 30U);
  applyItemsWrite(m_energy_fine_items, 31U);
}

void MhdLevelStepper::reconcileFaces(bool compare_only, double* max_mismatch) {
  m_compare_max = 0.0;
  const auto self_pass = [&]() {
    for (const Item& item : m_face_self_items) {
      PatchRuntime& patch = m_local[m_topo[item.dst].local];
      const MhdPatchGeometry& g = patch.binding.geometry;
      std::span<double> lane = faceLane(patch, item.axis);
      for (std::int64_t k = item.dst_box.lo[2]; k < item.dst_box.hi[2]; ++k) {
        for (std::int64_t j = item.dst_box.lo[1]; j < item.dst_box.hi[1]; ++j) {
          for (std::int64_t i = item.dst_box.lo[0]; i < item.dst_box.hi[0]; ++i) {
            const double value = lane[faceIndex(g, item.axis, i + item.src_offset[0],
                                                j + item.src_offset[1], k + item.src_offset[2])];
            double& target = lane[faceIndex(g, item.axis, i, j, k)];
            if (compare_only) {
              m_compare_max = std::max(m_compare_max, std::abs(target - value));
            } else {
              target = value;
            }
          }
        }
      }
    }
  };
  collective("faces.self", self_pass);
  if (!compare_only) {
    applyItemsWrite(m_face_cross_items, 40U);
    applyItemsWrite(m_face_fine_items, 41U);
    return;
  }
  const auto compare = [&](const Item& item, const double* values) {
    PatchRuntime& patch = m_local[m_topo[item.dst].local];
    const MhdPatchGeometry& g = patch.binding.geometry;
    const std::span<const double> lane = faceLane(patch, item.axis);
    std::size_t cursor = 0U;
    for (std::int64_t k = item.dst_box.lo[2]; k < item.dst_box.hi[2]; ++k) {
      for (std::int64_t j = item.dst_box.lo[1]; j < item.dst_box.hi[1]; ++j) {
        for (std::int64_t i = item.dst_box.lo[0]; i < item.dst_box.hi[0]; ++i) {
          m_compare_max = std::max(
              m_compare_max, std::abs(lane[faceIndex(g, item.axis, i, j, k)] - values[cursor++]));
        }
      }
    }
  };
  const auto pack = [&](const Item& item, std::vector<double>& out) {
    computeItemValues(item, 1.0, out);
  };
  runItems(42U, std::span<const Item>(m_face_cross_items), pack, compare, false);
  runItems(43U, std::span<const Item>(m_face_fine_items), pack, compare, false);
  if (max_mismatch != nullptr) {
    *max_mismatch = m_transport->size() > 1 ? m_transport->allreduceMaxDouble(m_compare_max)
                                            : m_compare_max;
  }
}

double MhdLevelStepper::validateCanonicalConsistency() {
  double scale = 0.0;
  collective("validate.scale", [&]() {
    for (const PatchRuntime& patch : m_local) {
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        for (const double value : faceLane(patch, axis)) {
          if (!std::isfinite(value)) {
            MhdFailureContext context =
                makeContext("mhd_level_stepper", "non-finite canonical face flux");
            context.patch_id = patch.binding.patch_id;
            throw MhdError(MhdErrorCode::kInvalidField, std::move(context));
          }
          scale = std::max(scale, std::abs(value));
        }
      }
    }
  });
  if (m_transport->size() > 1) {
    scale = m_transport->allreduceMaxDouble(scale);
  }
  double mismatch = 0.0;
  reconcileFaces(true, &mismatch);
  const double normalizer = std::max(scale, std::numeric_limits<double>::min());
  if (mismatch > 1.0e-12 * normalizer) {
    std::ostringstream stream;
    stream << "shared face mirrors or coarse interface faces disagree with their owner beyond "
              "roundoff (max |difference| "
           << mismatch << ", flux scale " << scale << "); the canonical state is not consistent";
    throw MhdError(MhdErrorCode::kInvalidField, makeContext("mhd_level_stepper", stream.str()));
  }
  return mismatch / normalizer;
}

// ---------------------------------------------------------------------------
// Allocation, snapshots and rollback
// ---------------------------------------------------------------------------

void MhdLevelStepper::requireModeAdmitted(const MhdLevelAdvanceOptions& options) const {
  const bool admitted = m_admission.admitted();
  if (m_config.require_memory_admission && !admitted) {
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded,
                   makeContext("mhd_level_stepper",
                               "admitMemory() must succeed before any MHD workspace is allocated"));
  }
  if (!admitted) {
    return;
  }
  const bool nonideal = options.nonideal != nullptr;
  const bool super = nonideal && options.nonideal->diffusionEnabled();
  if ((nonideal && !m_admitted_nonideal) || (super && !m_admitted_super) ||
      (options.produce_edge_integrals && !m_admitted_mode.produce_edge_integrals) ||
      (options.produce_boundary_flux_integrals &&
       !m_admitted_mode.produce_boundary_flux_integrals)) {
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded,
                   makeContext("mhd_level_stepper",
                               "interval mode needs workspaces that were not admitted"));
  }
}

void MhdLevelStepper::ensureAllocated(const MhdLevelAdvanceOptions& options) {
  const bool nonideal = options.nonideal != nullptr;
  const bool super = nonideal && options.nonideal->diffusionEnabled();
  const bool integral = options.produce_edge_integrals;
  const bool flux_integrals = options.produce_boundary_flux_integrals;
  const bool first = m_local.empty() ? false : m_local.front().rollback_gas.size() == 0U;
  const bool grow = first || (nonideal && !m_allocated_nonideal) || (super && !m_allocated_super) ||
      (integral && !m_allocated_integral) || (flux_integrals && !m_allocated_flux_integrals);
  if (!grow) {
    return;
  }
  for (PatchRuntime& patch : m_local) {
    const MhdPatchGeometry& g = patch.binding.geometry;
    const std::size_t cells = g.cellCount();
    const std::array<std::size_t, 3> faces{g.faceCountX(), g.faceCountY(), g.faceCountZ()};
    if (first) {
      patch.rollback_gas.resize(cells);
      patch.stage0_gas.resize(cells);
      patch.increment.resize(cells);
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        patch.rollback_faces[axis].resize(faces[axis], 0.0);
        patch.stage0_faces[axis].resize(faces[axis], 0.0);
      }
      patch.emf.resize(g);
      patch.records.resize(g);
      patch.halo_import.resize(static_cast<std::size_t>(patch.halo.import_doubles), 0.0);
    }
    if (nonideal && !m_allocated_nonideal) {
      patch.energy_rate.resize(cells, 0.0);
      patch.energy_records.resize(g);
    }
    if (super && !m_allocated_super) {
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        patch.yjm2_faces[axis].resize(faces[axis], 0.0);
      }
      patch.yjm2_energy.resize(cells, 0.0);
      patch.energy_rate_y0.resize(cells, 0.0);
      patch.emf_y0.resize(g);
    }
    if (integral && !m_allocated_integral) {
      patch.integral.resize(g);
      patch.stage0_integral.resize(g);
      if (super || m_allocated_super) {
        patch.yjm2_integral.resize(g);
      }
    }
    if (integral && super && patch.yjm2_integral.emfX().size() != g.edgeCountX()) {
      patch.yjm2_integral.resize(g);
    }
    if (flux_integrals && !m_allocated_flux_integrals) {
      patch.gas_integral.resize(g);
      patch.stage0_gas_integral.resize(g);
      patch.rollback_gas_integral.resize(g);
    }
    if (flux_integrals && nonideal && patch.energy_integral.values[0][0].size() == 0U &&
        g.ny * g.nz != 0U) {
      patch.energy_integral.resize(g);
      patch.stage0_energy_integral.resize(g);
      patch.rollback_energy_integral.resize(g);
    }
    if (flux_integrals && super && patch.energy_records_y0.values[0][0].size() == 0U &&
        g.ny * g.nz != 0U) {
      patch.yjm2_energy_integral.resize(g);
      patch.energy_records_y0.resize(g);
    }
  }
  // Batch slots grow to the per-lane maximum over the local patch shapes and
  // keep that capacity (MhdDoubleArray retains capacity across shapes).
  std::uint32_t max_cell_samples = 0U;
  std::uint32_t max_face_samples = 0U;
  std::size_t max_item_values = 0U;
  for (const PatchRuntime& patch : m_local) {
    max_cell_samples = std::max(max_cell_samples, patch.halo.cell_sample_count);
    max_face_samples = std::max(max_face_samples, patch.halo.face_sample_count);
    for (const Item& item : patch.halo.items) {
      max_item_values = std::max(max_item_values, itemValueCount(item));
    }
  }
  for (BatchSlot& slot : m_slots) {
    for (const PatchRuntime& patch : m_local) {
      slot.workspace.resize(patch.binding.geometry, m_policy.reconstruction);
      if (nonideal) {
        slot.nonideal.resize(patch.binding.geometry, false);
      }
    }
    if (slot.cell_samples.size() < 9U * static_cast<std::size_t>(max_cell_samples)) {
      slot.cell_samples.resize(9U * static_cast<std::size_t>(max_cell_samples), 0.0);
    }
    if (slot.face_samples.size() < static_cast<std::size_t>(max_face_samples)) {
      slot.face_samples.resize(static_cast<std::size_t>(max_face_samples), 0.0);
    }
    slot.values.reserve(max_item_values);
  }
  std::size_t max_exchange_item = 0U;
  for (const std::vector<Item>* list :
       {&m_face_cross_items, &m_edge_cross_items, &m_record_cross_items, &m_face_fine_items,
        &m_record_fine_items, &m_energy_cross_items, &m_energy_fine_items}) {
    for (const Item& item : *list) {
      max_exchange_item = std::max(max_exchange_item, itemValueCount(item));
    }
  }
  m_exchange_scratch.reserve(std::max(max_exchange_item, max_item_values));
  m_allocated_nonideal = m_allocated_nonideal || nonideal;
  m_allocated_super = m_allocated_super || super;
  m_allocated_integral = m_allocated_integral || integral;
  m_allocated_flux_integrals = m_allocated_flux_integrals || flux_integrals;
  if (m_admission.admitted()) {
    if (!m_admission_committed) {
      // Verifies measured retained capacity <= admitted bytes, then commits.
      m_admission.commit(ownedCapacityBytes());
      if (m_comm_reservation.pending()) {
        m_comm_reservation.commit();
      }
      m_admission_committed = true;
    } else if (ownedCapacityBytes() > m_admitted_plan.phasePlusScratchBytes()) {
      MhdFailureContext failure = makeContext(
          "mhd_level_stepper",
          "retained MHD capacity grew beyond the admitted plan (unexpected allocation growth)");
      failure.requested_bytes = ownedCapacityBytes();
      failure.admitted_bytes = m_admitted_plan.phasePlusScratchBytes();
      throw MhdError(MhdErrorCode::kMemoryBudgetExceeded, std::move(failure));
    }
  }
}

void MhdLevelStepper::snapshotRollback() {
  for (PatchRuntime& patch : m_local) {
    for (std::size_t lane = 0; lane < MhdGasStateSoa::k_lane_count; ++lane) {
      copyLane(patch.rollback_gas.lane(lane), patch.binding.gas->lane(lane));
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      copyLane(patch.rollback_faces[axis].span(), faceLane(patch, axis));
    }
  }
}

// No-throw restoration from preallocated images (sizes fixed at allocation).
void MhdLevelStepper::restoreRollback() noexcept {
  for (PatchRuntime& patch : m_local) {
    for (std::size_t lane = 0; lane < MhdGasStateSoa::k_lane_count; ++lane) {
      copyLane(patch.binding.gas->lane(lane), patch.rollback_gas.lane(lane));
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      copyLane(faceLane(patch, axis), patch.rollback_faces[axis].span());
    }
    if (m_track_integral) {
      patch.integral.zero();
    }
    if (m_track_flux_integrals) {
      zeroRecords(patch.gas_integral);
      if (m_nonideal != nullptr) {
        zeroRecords(patch.energy_integral);
      }
    }
  }
}

void MhdLevelStepper::snapshotStage0(bool gas_all_lanes) {
  for (PatchRuntime& patch : m_local) {
    if (gas_all_lanes) {
      for (std::size_t lane = 0; lane < MhdGasStateSoa::k_lane_count; ++lane) {
        copyLane(patch.stage0_gas.lane(lane), patch.binding.gas->lane(lane));
      }
    } else {
      copyLane(patch.stage0_gas.totalEnergy(), patch.binding.gas->totalEnergy());
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      copyLane(patch.stage0_faces[axis].span(), faceLane(patch, axis));
    }
    if (m_track_integral) {
      copyEdgeField(patch.stage0_integral, patch.integral);
    }
    if (m_track_flux_integrals) {
      copyRecords(patch.stage0_gas_integral, patch.gas_integral);
      if (m_nonideal != nullptr) {
        copyRecords(patch.stage0_energy_integral, patch.energy_integral);
      }
    }
  }
}

// magnetic := alpha * stage0 + beta * magnetic   (faces and tracked integral)
void MhdLevelStepper::combineMagnetic(PatchRuntime& patch, double alpha, double beta) {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    std::span<double> lane = faceLane(patch, axis);
    const MhdDoubleArray& base = patch.stage0_faces[axis];
    for (std::size_t index = 0; index < lane.size(); ++index) {
      lane[index] = alpha * base[index] + beta * lane[index];
    }
  }
  if (m_track_integral) {
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      std::span<double> lane = emfLane(patch.integral, axis);
      const std::span<const double> base = emfLane(patch.stage0_integral, axis);
      for (std::size_t index = 0; index < lane.size(); ++index) {
        lane[index] = alpha * base[index] + beta * lane[index];
      }
    }
  }
}

// Phi -= weight * C(emf); the tracked time integral gains weight * emf so that
// Phi_end = Phi_begin - C(I) holds exactly for every linear stage combination.
void MhdLevelStepper::applyCt(PatchRuntime& patch, const MhdEdgeEmfStorage& emf, double weight) {
  applyInductionUpdate(patch.binding.geometry, patch.binding.faces.phi_x,
                       patch.binding.faces.phi_y, patch.binding.faces.phi_z, emf, weight);
  if (m_track_integral) {
    accumulateEdgeEmfIntegral(patch.binding.geometry, patch.integral, emf, weight);
  }
}

// ---------------------------------------------------------------------------
// Ideal SSPRK2 (operator I)
// ---------------------------------------------------------------------------

void MhdLevelStepper::idealStage(double alpha, double beta, double weight, double scale_factor,
                                 double hubble, bool expansion, MhdStageCounters& counters,
                                 MhdLevelAdvanceResult& result) {
  exchangeHalo(scale_factor);
  collective("ideal.evaluate", [&]() {
    forEachLocalPatchBatched([&](std::size_t index, BatchSlot& slot) {
      PatchRuntime& patch = m_local[index];
      prepareSlot(index, slot, scale_factor);
      for (std::size_t lane = 0; lane < MhdGasStateSoa::k_lane_count; ++lane) {
        std::span<double> values = patch.increment.lane(lane);
        std::fill(values.begin(), values.end(), 0.0);
      }
      m_operator.evaluateIdealStage(patch.binding.geometry, patch.boundary, scale_factor, hubble,
                                    expansion, weight, slot.workspace, patch.increment,
                                    slot.counters);
      copyEdgeField(patch.emf, slot.workspace.emf);
      copyRecords(patch.records, slot.workspace.boundary_fluxes);
    });
    for (BatchSlot& slot : m_slots) {
      counters += slot.counters;
      slot.counters = MhdStageCounters{};
      result.ghost_prolongation_slope_fallbacks += slot.slope_fallbacks;
      slot.slope_fallbacks = 0U;
    }
  });
  ++result.ideal_stage_evaluations_attempted;
  reconcileGasRecords();
  reconcileEdges();
  collective("ideal.apply", [&]() {
    for (PatchRuntime& patch : m_local) {
      MhdStepOperator::applyBoundaryFaceFluxes(patch.binding.geometry, patch.records, scale_factor,
                                               weight, patch.increment);
      if (m_track_flux_integrals) {
        combineRecords(patch.gas_integral, patch.stage0_gas_integral, alpha, beta, &patch.records,
                       weight / scale_factor);
      }
      for (std::size_t lane = 0; lane < MhdGasStateSoa::k_lane_count; ++lane) {
        std::span<double> values = patch.binding.gas->lane(lane);
        const std::span<const double> base = patch.stage0_gas.lane(lane);
        const std::span<const double> delta = patch.increment.lane(lane);
        for (std::size_t cell = 0; cell < values.size(); ++cell) {
          values[cell] = alpha * base[cell] + beta * values[cell] + delta[cell];
        }
      }
      combineMagnetic(patch, alpha, beta);
      applyCt(patch, patch.emf, weight);
    }
  });
  reconcileFaces(false, nullptr);
}

void MhdLevelStepper::idealSsprk2(const MhdStageContext& interval, MhdStageCounters& counters,
                                  MhdLevelAdvanceResult& result) {
  const double dt = interval.dt_code;
  const bool expansion = interval.comoving_coordinates && interval.apply_expansion_source;
  snapshotStage0(true);
  // Stage 1: Y1 = Y0 + dt L(Y0, t_n).
  idealStage(0.0, 1.0, dt, interval.scale_factor, interval.hubble_rate_code, expansion, counters,
             result);
  // Stage 2: Y = Y0/2 + (Y1 + dt L(Y1, t_n + dt))/2.
  idealStage(0.5, 0.5, 0.5 * dt, interval.scale_factor_end, interval.hubble_rate_end_code,
             expansion, counters, result);
}

// ---------------------------------------------------------------------------
// Non-ideal operators D and H
// ---------------------------------------------------------------------------

double MhdLevelStepper::evaluateNonIdealRate(MhdNonIdealOperatorKind kind, double scale_factor,
                                             MhdTransportStageDiagnostics& diagnostics,
                                             MhdLevelAdvanceResult& result) {
  exchangeHalo(scale_factor);
  double local_stable = std::numeric_limits<double>::infinity();
  collective("nonideal.emf", [&]() {
    for (BatchSlot& slot : m_slots) {
      slot.stable_dt = std::numeric_limits<double>::infinity();
    }
    forEachLocalPatchBatched([&](std::size_t index, BatchSlot& slot) {
      PatchRuntime& patch = m_local[index];
      prepareSlot(index, slot, scale_factor);
      slot.nonideal.resize(patch.binding.geometry, false);
      m_nonideal->evaluateCoefficients(patch.binding.geometry, slot.workspace.primitives,
                                       scale_factor, slot.nonideal, slot.transport);
      slot.stable_dt = std::min(
          slot.stable_dt,
          kind == MhdNonIdealOperatorKind::kDiffusion
              ? m_nonideal->diffusionStableTimeStep(patch.binding.geometry, slot.nonideal)
              : m_nonideal->hallStableTimeStep(patch.binding.geometry, slot.nonideal));
      m_nonideal->computeEmf(kind, patch.binding.geometry, patch.boundary, slot.workspace.staggered,
                             scale_factor, slot.nonideal, slot.workspace.emf, slot.transport);
      copyEdgeField(patch.emf, slot.workspace.emf);
    });
    for (BatchSlot& slot : m_slots) {
      local_stable = std::min(local_stable, slot.stable_dt);
      diagnostics += slot.transport;
      slot.transport = MhdTransportStageDiagnostics{};
      result.ghost_prolongation_slope_fallbacks += slot.slope_fallbacks;
      slot.slope_fallbacks = 0U;
    }
  });
  const double stable =
      m_transport->size() > 1 ? m_transport->allreduceMinDouble(local_stable) : local_stable;
  reconcileEdges();
  collective("nonideal.energy", [&]() {
    forEachLocalPatchBatched([&](std::size_t index, BatchSlot& slot) {
      PatchRuntime& patch = m_local[index];
      prepareSlot(index, slot, scale_factor);
      slot.nonideal.resize(patch.binding.geometry, false);
      std::fill(patch.energy_rate.begin(), patch.energy_rate.end(), 0.0);
      m_nonideal->applyEnergyFlux(patch.binding.geometry, slot.workspace.primitives, patch.emf,
                                  scale_factor, 1.0, slot.nonideal, patch.energy_rate.span());
      copyRecords(patch.energy_records, slot.nonideal.boundary_energy);
    });
    for (BatchSlot& slot : m_slots) {
      slot.slope_fallbacks = 0U;  // already counted by the EMF pass
    }
  });
  reconcileEnergyRecords();
  collective("nonideal.energy_apply", [&]() {
    for (PatchRuntime& patch : m_local) {
      applyBoundaryEnergyFluxes(patch.binding.geometry, patch.energy_records, 1.0,
                                patch.energy_rate.span());
    }
  });
  return stable;
}

void MhdLevelStepper::explicitNonIdealSubstep(MhdNonIdealOperatorKind kind,
                                              const MhdStageContext& sub,
                                              MhdTransportStageDiagnostics& diagnostics,
                                              MhdLevelAdvanceResult& result) {
  const double h = sub.dt_code;
  const auto require_stable = [&](double stable) {
    if (h > stable * (1.0 + 1.0e-12)) {
      MhdFailureContext failure = makeContext(
          kind == MhdNonIdealOperatorKind::kDiffusion ? "mhd_diffusion" : "mhd_hall",
          "stage coefficients grew beyond the forecast explicit stability bound");
      failure.time_code = sub.time_code;
      failure.scale_factor = sub.scale_factor;
      throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
    }
  };
  snapshotStage0(false);
  require_stable(evaluateNonIdealRate(kind, sub.scale_factor, diagnostics, result));
  collective("nonideal.substep1", [&]() {
    for (PatchRuntime& patch : m_local) {
      std::span<double> energy = patch.binding.gas->totalEnergy();
      for (std::size_t cell = 0; cell < energy.size(); ++cell) {
        energy[cell] += h * patch.energy_rate[cell];
      }
      if (m_track_flux_integrals) {
        combineRecords(patch.energy_integral, patch.stage0_energy_integral, 0.0, 1.0,
                       &patch.energy_records, h);
      }
      applyCt(patch, patch.emf, h);
    }
  });
  reconcileFaces(false, nullptr);
  require_stable(evaluateNonIdealRate(kind, sub.scale_factor_end, diagnostics, result));
  collective("nonideal.substep2", [&]() {
    for (PatchRuntime& patch : m_local) {
      std::span<double> energy = patch.binding.gas->totalEnergy();
      const std::span<const double> energy0 = patch.stage0_gas.totalEnergy();
      for (std::size_t cell = 0; cell < energy.size(); ++cell) {
        energy[cell] = 0.5 * energy0[cell] + 0.5 * (energy[cell] + h * patch.energy_rate[cell]);
      }
      if (m_track_flux_integrals) {
        combineRecords(patch.energy_integral, patch.stage0_energy_integral, 0.5, 0.5,
                       &patch.energy_records, 0.5 * h);
      }
      combineMagnetic(patch, 0.5, 0.5);
      applyCt(patch, patch.emf, 0.5 * h);
    }
  });
  reconcileFaces(false, nullptr);
}

void MhdLevelStepper::rkl2Cycle(const MhdStageContext& sub, std::int32_t stages,
                                MhdTransportStageDiagnostics& diagnostics,
                                MhdLevelAdvanceResult& result) {
  const MhdRkl2Coefficients coefficients = makeRkl2Coefficients(stages);
  const double tau = sub.dt_code;
  // All stages of one cycle use the cycle-midpoint scale factor for the
  // comoving coefficient conversion (second order in the slowly varying a);
  // physical properties are re-evaluated from every stage state.
  const double a_mid = 0.5 * (sub.scale_factor + sub.scale_factor_end);
  const double bound = rkl2StabilityFactor(stages);
  const auto require_stable = [&](double stable) {
    if (tau > stable * bound * (1.0 + 1.0e-12)) {
      MhdFailureContext failure =
          makeContext("mhd_rkl2", "stage coefficients grew beyond the RKL2 stability interval");
      failure.predicted_stages = static_cast<std::uint64_t>(stages);
      failure.time_code = sub.time_code;
      throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
    }
  };
  snapshotStage0(false);  // Y_0 (energy + faces + integrals)
  require_stable(
      evaluateNonIdealRate(MhdNonIdealOperatorKind::kDiffusion, a_mid, diagnostics, result));
  collective("rkl2.stage1", [&]() {
    for (PatchRuntime& patch : m_local) {
      // Keep L(Y_0) and Y_{j-2} := Y_0.
      copyEdgeField(patch.emf_y0, patch.emf);
      copyLane(patch.energy_rate_y0.span(), patch.energy_rate.span());
      copyLane(patch.yjm2_energy.span(), patch.stage0_gas.totalEnergy());
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        copyLane(patch.yjm2_faces[axis].span(), patch.stage0_faces[axis].span());
      }
      if (m_track_integral) {
        copyEdgeField(patch.yjm2_integral, patch.stage0_integral);
      }
      if (m_track_flux_integrals) {
        copyRecords(patch.energy_records_y0, patch.energy_records);
        copyRecords(patch.yjm2_energy_integral, patch.stage0_energy_integral);
        combineRecords(patch.energy_integral, patch.stage0_energy_integral, 0.0, 1.0,
                       &patch.energy_records_y0, coefficients.mu_tilde[1] * tau);
      }
      // Y_1 = Y_0 + mu~_1 tau L(Y_0).
      std::span<double> energy = patch.binding.gas->totalEnergy();
      for (std::size_t cell = 0; cell < energy.size(); ++cell) {
        energy[cell] += coefficients.mu_tilde[1] * tau * patch.energy_rate_y0[cell];
      }
      applyCt(patch, patch.emf_y0, coefficients.mu_tilde[1] * tau);
    }
  });
  reconcileFaces(false, nullptr);
  ++diagnostics.rkl_stage_count;
  for (std::int32_t j = 2; j <= stages; ++j) {
    const auto index = static_cast<std::size_t>(j);
    const double mu = coefficients.mu[index];
    const double nu = coefficients.nu[index];
    const double y0_weight = 1.0 - mu - nu;
    const double mu_tilde_tau = coefficients.mu_tilde[index] * tau;
    const double gamma_tilde_tau = coefficients.gamma_tilde[index] * tau;
    require_stable(
        evaluateNonIdealRate(MhdNonIdealOperatorKind::kDiffusion, a_mid, diagnostics, result));
    collective("rkl2.stage", [&]() {
      for (PatchRuntime& patch : m_local) {
        // Y_j = mu Y_{j-1} + nu Y_{j-2} + (1-mu-nu) Y_0 + mu~ tau L(Y_{j-1})
        //       + gamma~ tau L(Y_0); rotate Y_{j-2}.
        std::span<double> energy = patch.binding.gas->totalEnergy();
        const std::span<const double> energy0 = patch.stage0_gas.totalEnergy();
        for (std::size_t cell = 0; cell < energy.size(); ++cell) {
          const double previous = energy[cell];
          energy[cell] = mu * previous + nu * patch.yjm2_energy[cell] + y0_weight * energy0[cell] +
              mu_tilde_tau * patch.energy_rate[cell] + gamma_tilde_tau * patch.energy_rate_y0[cell];
          patch.yjm2_energy[cell] = previous;
        }
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          std::span<double> lane = faceLane(patch, axis);
          MhdDoubleArray& yjm2 = patch.yjm2_faces[axis];
          const MhdDoubleArray& y0 = patch.stage0_faces[axis];
          for (std::size_t face = 0; face < lane.size(); ++face) {
            const double previous = lane[face];
            lane[face] = mu * previous + nu * yjm2[face] + y0_weight * y0[face];
            yjm2[face] = previous;
          }
        }
        if (m_track_integral) {
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            std::span<double> lane = emfLane(patch.integral, axis);
            std::span<double> yjm2 = emfLane(patch.yjm2_integral, axis);
            const std::span<const double> y0 = emfLane(patch.stage0_integral, axis);
            for (std::size_t edge = 0; edge < lane.size(); ++edge) {
              const double previous = lane[edge];
              lane[edge] = mu * previous + nu * yjm2[edge] + y0_weight * y0[edge];
              yjm2[edge] = previous;
            }
          }
        }
        if (m_track_flux_integrals) {
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            for (std::size_t side = 0; side < 2U; ++side) {
              std::span<double> lane = patch.energy_integral.values[axis][side].span();
              std::span<double> yjm2 = patch.yjm2_energy_integral.values[axis][side].span();
              const std::span<const double> y0 =
                  patch.stage0_energy_integral.values[axis][side].span();
              const std::span<const double> rate = patch.energy_records.values[axis][side].span();
              const std::span<const double> rate0 =
                  patch.energy_records_y0.values[axis][side].span();
              for (std::size_t face = 0; face < lane.size(); ++face) {
                const double previous = lane[face];
                lane[face] = mu * previous + nu * yjm2[face] + y0_weight * y0[face] +
                    mu_tilde_tau * rate[face] + gamma_tilde_tau * rate0[face];
                yjm2[face] = previous;
              }
            }
          }
        }
        // Magnetic parts of the L terms (L_Phi = -C(E)).
        applyCt(patch, patch.emf, mu_tilde_tau);
        applyCt(patch, patch.emf_y0, gamma_tilde_tau);
      }
    });
    reconcileFaces(false, nullptr);
    ++diagnostics.rkl_stage_count;
  }
  ++diagnostics.rkl_cycle_count;
}

void MhdLevelStepper::diffusionHalf(const MhdStageContext& sub,
                                    MhdTransportStageDiagnostics& diagnostics,
                                    MhdLevelAdvanceResult& result) {
  if (!m_nonideal->diffusionEnabled() || !(sub.dt_code > 0.0)) {
    return;
  }
  const double stable = m_forecast_dt_diffusion;
  if (!std::isfinite(stable)) {
    return;  // all diffusion coefficients exactly zero: exact no-op
  }
  const double ratio = sub.dt_code / stable;
  if (ratio <= static_cast<double>(m_nonideal->limits().rkl_switch_ratio)) {
    // Explicit SSPRK2 with the least number of equal substeps h <= dt_D.
    const auto substeps =
        static_cast<std::int32_t>(std::max(1.0, std::ceil(ratio * (1.0 - 1.0e-12))));
    for (std::int32_t step = 0; step < substeps; ++step) {
      const double begin = static_cast<double>(step) / substeps;
      const double end = static_cast<double>(step + 1) / substeps;
      explicitNonIdealSubstep(MhdNonIdealOperatorKind::kDiffusion,
                              subIntervalContext(sub, begin, end), diagnostics, result);
      ++diagnostics.explicit_diffusion_substeps;
    }
    return;
  }
  const MhdRkl2Plan plan = planRkl2(ratio, m_nonideal->limits().rkl_max_stages);
  if (m_interval_rkl_stages + plan.totalStages() > m_nonideal->limits().rkl_stage_cost_cap) {
    MhdFailureContext failure =
        makeContext("mhd_rkl2", "RKL2 stage forecast exceeds the per-interval cost cap");
    failure.predicted_stages =
        static_cast<std::uint64_t>(m_interval_rkl_stages + plan.totalStages());
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
  }
  m_interval_rkl_stages += plan.totalStages();
  for (std::int32_t cycle = 0; cycle < plan.cycles; ++cycle) {
    const double begin = static_cast<double>(cycle) / plan.cycles;
    const double end = static_cast<double>(cycle + 1) / plan.cycles;
    rkl2Cycle(subIntervalContext(sub, begin, end), plan.stages_per_cycle, diagnostics, result);
  }
}

void MhdLevelStepper::hallHalf(const MhdStageContext& sub,
                               MhdTransportStageDiagnostics& diagnostics,
                               MhdLevelAdvanceResult& result) {
  if (!m_nonideal->hallEnabled() || !(sub.dt_code > 0.0)) {
    return;
  }
  const double stable = m_forecast_dt_hall;
  if (!std::isfinite(stable)) {
    return;  // eta_H exactly zero everywhere: exact no-op
  }
  const double ratio = sub.dt_code / stable;
  if (!std::isfinite(ratio) ||
      ratio > static_cast<double>(std::numeric_limits<std::int32_t>::max())) {
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded,
                   makeContext("mhd_hall", "Hall subcycle forecast is not representable"));
  }
  const auto subcycles =
      static_cast<std::int64_t>(std::max(1.0, std::ceil(ratio * (1.0 - 1.0e-12))));
  if (m_interval_hall_subcycles + subcycles > m_nonideal->limits().hall_subcycle_cost_cap) {
    MhdFailureContext failure =
        makeContext("mhd_hall", "Hall subcycle forecast exceeds the per-interval cost cap");
    failure.predicted_stages = static_cast<std::uint64_t>(m_interval_hall_subcycles + subcycles);
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
  }
  m_interval_hall_subcycles += subcycles;
  for (std::int64_t step = 0; step < subcycles; ++step) {
    const double begin = static_cast<double>(step) / static_cast<double>(subcycles);
    const double end = static_cast<double>(step + 1) / static_cast<double>(subcycles);
    explicitNonIdealSubstep(MhdNonIdealOperatorKind::kHall, subIntervalContext(sub, begin, end),
                            diagnostics, result);
    ++diagnostics.hall_subcycle_count;
  }
}

// Pre-mutation stiffness forecast at the interval-start state (rank-global
// minima): the explicit bounds and, when requested, the per-interval Hall/RKL
// cost caps for both halves. Only scratch is written.
void MhdLevelStepper::forecastNonIdeal(const MhdStageContext& interval,
                                       MhdTransportStageDiagnostics& diagnostics,
                                       MhdLevelAdvanceResult& result, bool enforce_caps) {
  exchangeHalo(interval.scale_factor);
  double dt_diffusion = std::numeric_limits<double>::infinity();
  double dt_hall = std::numeric_limits<double>::infinity();
  collective("nonideal.forecast", [&]() {
    std::array<double, k_max_batch_slots> slot_d{};
    std::array<double, k_max_batch_slots> slot_h{};
    slot_d.fill(std::numeric_limits<double>::infinity());
    slot_h.fill(std::numeric_limits<double>::infinity());
    forEachLocalPatchBatched([&](std::size_t index, BatchSlot& slot) {
      PatchRuntime& patch = m_local[index];
      prepareSlot(index, slot, interval.scale_factor);
      slot.nonideal.resize(patch.binding.geometry, false);
      m_nonideal->evaluateCoefficients(patch.binding.geometry, slot.workspace.primitives,
                                       interval.scale_factor, slot.nonideal, slot.transport);
      const auto s = static_cast<std::size_t>(&slot - m_slots.data());
      slot_d[s] = std::min(slot_d[s], m_nonideal->diffusionStableTimeStep(patch.binding.geometry,
                                                                           slot.nonideal));
      slot_h[s] = std::min(slot_h[s],
                           m_nonideal->hallStableTimeStep(patch.binding.geometry, slot.nonideal));
    });
    for (std::size_t s = 0; s < m_slots.size(); ++s) {
      dt_diffusion = std::min(dt_diffusion, slot_d[s]);
      dt_hall = std::min(dt_hall, slot_h[s]);
      diagnostics += m_slots[s].transport;
      m_slots[s].transport = MhdTransportStageDiagnostics{};
      result.ghost_prolongation_slope_fallbacks += m_slots[s].slope_fallbacks;
      m_slots[s].slope_fallbacks = 0U;
    }
  });
  if (m_transport->size() > 1) {
    dt_diffusion = m_transport->allreduceMinDouble(dt_diffusion);
    dt_hall = m_transport->allreduceMinDouble(dt_hall);
  }
  m_forecast_dt_diffusion =
      m_nonideal->diffusionEnabled() ? dt_diffusion : std::numeric_limits<double>::infinity();
  m_forecast_dt_hall =
      m_nonideal->hallEnabled() ? dt_hall : std::numeric_limits<double>::infinity();
  diagnostics.min_diffusion_dt_code =
      std::isfinite(m_forecast_dt_diffusion) ? m_forecast_dt_diffusion : 0.0;
  diagnostics.min_hall_dt_code = std::isfinite(m_forecast_dt_hall) ? m_forecast_dt_hall : 0.0;
  if (!enforce_caps) {
    return;
  }
  const double half = 0.5 * interval.dt_code;
  if (std::isfinite(m_forecast_dt_hall)) {
    const double subcycles = 2.0 * std::ceil(half / m_forecast_dt_hall);
    if (subcycles > static_cast<double>(m_nonideal->limits().hall_subcycle_cost_cap)) {
      MhdFailureContext failure =
          makeContext("mhd_hall", "forecast Hall subcycles exceed the per-interval cost cap");
      failure.predicted_stages = static_cast<std::uint64_t>(std::min(subcycles, 1.0e18));
      throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
    }
  }
  if (std::isfinite(m_forecast_dt_diffusion)) {
    const double ratio = half / m_forecast_dt_diffusion;
    if (ratio > static_cast<double>(m_nonideal->limits().rkl_switch_ratio)) {
      const MhdRkl2Plan plan = planRkl2(ratio, m_nonideal->limits().rkl_max_stages);
      if (2 * plan.totalStages() > m_nonideal->limits().rkl_stage_cost_cap) {
        MhdFailureContext failure =
            makeContext("mhd_rkl2", "forecast RKL2 stages exceed the per-interval cost cap");
        failure.predicted_stages = static_cast<std::uint64_t>(2 * plan.totalStages());
        throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
      }
    }
  }
}

void MhdLevelStepper::validateAccepted(double scale_factor, MhdLevelAdvanceResult& result) {
  double divergence = 0.0;
  collective("validate.accepted", [&]() {
    for (PatchRuntime& patch : m_local) {
      const MhdPatchGeometry& g = patch.binding.geometry;
      for (std::size_t k = 0; k < g.nz; ++k) {
        for (std::size_t j = 0; j < g.ny; ++j) {
          for (std::size_t i = 0; i < g.nx; ++i) {
            MhdPrimitiveState w{};
            if (!m_operator.recoverCellPrimitive(
                    patch.binding.gas->loadGas(g.cellIndex(i, j, k)),
                    MhdStepOperator::cellFieldFromFaces(g, patch.binding.faces.phi_x,
                                                        patch.binding.faces.phi_y,
                                                        patch.binding.faces.phi_z, i, j, k),
                    scale_factor, w)) {
              MhdFailureContext failure =
                  makeContext("mhd_level_stepper", "inadmissible cell(s) in the candidate end state");
              failure.patch_id = patch.binding.patch_id;
              failure.scale_factor = scale_factor;
              throw MhdError(MhdErrorCode::kNonPositiveState, std::move(failure));
            }
          }
        }
      }
      const double normalizer = magneticDivergenceNormalizer(
          g, patch.binding.faces.phi_x, patch.binding.faces.phi_y, patch.binding.faces.phi_z);
      const double imbalance = maxFluxImbalance(g, patch.binding.faces.phi_x,
                                                patch.binding.faces.phi_y,
                                                patch.binding.faces.phi_z);
      const double norm = normalizer > 0.0 && std::isfinite(normalizer)
          ? imbalance / normalizer
          : (imbalance > 0.0 ? std::numeric_limits<double>::infinity() : 0.0);
      divergence = std::max(divergence, norm);
      MhdFailureContext context = makeContext("mhd_level_stepper", "post-interval strict CT");
      context.patch_id = patch.binding.patch_id;
      context.geometry_generation = patch.binding.geometry_generation;
      context.scale_factor = scale_factor;
      requireDivergenceContract(norm, m_policy.divergence_tolerance, context);
    }
  });
  result.divergence_max_normalized =
      m_transport->size() > 1 ? m_transport->allreduceMaxDouble(divergence) : divergence;
}

// ---------------------------------------------------------------------------
// Memory plan and governed admission
// ---------------------------------------------------------------------------

std::uint64_t MhdLevelStepper::topologyCapacityBytes() const noexcept {
  std::uint64_t bytes = m_topo.capacity() * sizeof(TopoPatch) +
      m_local.capacity() * sizeof(PatchRuntime) + m_slots.capacity() * sizeof(BatchSlot);
  for (const auto& bins : m_bins) {
    bytes += bins.capacity() * sizeof(bins.front());
  }
  for (const std::vector<Item>* list :
       {&m_halo_items, &m_face_self_items, &m_face_cross_items, &m_edge_self_items,
        &m_edge_cross_items, &m_record_self_items, &m_record_cross_items, &m_face_fine_items,
        &m_record_fine_items, &m_energy_self_items, &m_energy_cross_items, &m_energy_fine_items}) {
    bytes += list->capacity() * sizeof(Item);
  }
  for (const std::vector<Item>& list : m_edge_fine_items) {
    bytes += list.capacity() * sizeof(Item);
  }
  for (const PatchRuntime& patch : m_local) {
    bytes += patch.halo.items.capacity() * sizeof(Item) +
        patch.halo.prolong_cells.capacity() * sizeof(ProlongCell) +
        patch.halo.prolong_faces.capacity() * sizeof(ProlongFace);
  }
  return bytes;
}

std::uint64_t MhdLevelStepper::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (const PatchRuntime& patch : m_local) {
    bytes += patch.halo_import.ownedCapacityBytes() + patch.rollback_gas.ownedCapacityBytes() +
        patch.stage0_gas.ownedCapacityBytes() + patch.increment.ownedCapacityBytes() +
        patch.emf.ownedCapacityBytes() + patch.records.ownedCapacityBytes() +
        patch.energy_rate.ownedCapacityBytes() + patch.energy_records.ownedCapacityBytes() +
        patch.yjm2_energy.ownedCapacityBytes() + patch.energy_rate_y0.ownedCapacityBytes() +
        patch.emf_y0.ownedCapacityBytes() + patch.integral.ownedCapacityBytes() +
        patch.stage0_integral.ownedCapacityBytes() + patch.yjm2_integral.ownedCapacityBytes() +
        patch.gas_integral.ownedCapacityBytes() + patch.stage0_gas_integral.ownedCapacityBytes() +
        patch.rollback_gas_integral.ownedCapacityBytes() +
        patch.energy_integral.ownedCapacityBytes() +
        patch.stage0_energy_integral.ownedCapacityBytes() +
        patch.yjm2_energy_integral.ownedCapacityBytes() +
        patch.energy_records_y0.ownedCapacityBytes() +
        patch.rollback_energy_integral.ownedCapacityBytes();
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      bytes += patch.rollback_faces[axis].ownedCapacityBytes() +
          patch.stage0_faces[axis].ownedCapacityBytes() + patch.yjm2_faces[axis].ownedCapacityBytes();
    }
  }
  for (const BatchSlot& slot : m_slots) {
    bytes += slot.workspace.ownedCapacityBytes() + slot.nonideal.ownedCapacityBytes() +
        slot.cell_samples.ownedCapacityBytes() + slot.face_samples.ownedCapacityBytes() +
        slot.values.capacity() * sizeof(double);
  }
  bytes += m_exchange_scratch.capacity() * sizeof(double);
  return bytes;
}

MhdLevelMemoryPlan MhdLevelStepper::planMemory(const MhdLevelAdvanceOptions& mode,
                                               const MhdResourceLimits& limits) const {
  (void)limits;
  MhdLevelMemoryPlan plan{};
  const bool nonideal = mode.nonideal != nullptr;
  const bool super = nonideal && mode.nonideal->diffusionEnabled();
  const auto lane = [](std::size_t count) { return mhdAlignedLaneBytes(count); };
  const auto record_bytes = [&](const MhdPatchGeometry& g, std::size_t per_face) {
    std::uint64_t bytes = 0U;
    const std::array<std::size_t, 3> plane{g.ny * g.nz, g.nx * g.nz, g.nx * g.ny};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      bytes += 2U * lane(per_face * plane[axis]);
    }
    return bytes;
  };
  // Per-lane maxima of the batch workspace over local shapes (retained
  // capacity of a reusable slot: every MhdDoubleArray keeps its largest size).
  constexpr std::size_t k_slot_lanes = 9U + 3U + 6U + 3U + 6U + 3U + 1U + 6U;
  std::array<std::size_t, k_slot_lanes> slot_counts{};
  std::uint32_t max_cell_samples = 0U;
  std::uint32_t max_face_samples = 0U;
  std::size_t max_item_values = 0U;
  const std::size_t g = m_policy.ghostCells();
  const std::size_t gt = k_mhd_staggered_tangential_ghost;
  for (const PatchRuntime& patch : m_local) {
    const MhdPatchGeometry& geo = patch.binding.geometry;
    const std::size_t cells = geo.cellCount();
    const std::uint64_t gas6 = 6U * lane(cells);
    const std::uint64_t faces3 =
        lane(geo.faceCountX()) + lane(geo.faceCountY()) + lane(geo.faceCountZ());
    const std::uint64_t edges3 =
        lane(geo.edgeCountX()) + lane(geo.edgeCountY()) + lane(geo.edgeCountZ());
    const std::uint64_t gas_records = record_bytes(geo, MhdBoundaryFaceFluxes::k_values_per_face);
    const std::uint64_t energy_records = record_bytes(geo, 1U);
    std::uint64_t resident = 2U * (gas6 + faces3) + gas6 + edges3 + gas_records +
        lane(static_cast<std::size_t>(patch.halo.import_doubles));
    if (nonideal) {
      resident += lane(cells) + energy_records;
    }
    if (super) {
      resident += faces3 + 2U * lane(cells) + edges3;
    }
    if (mode.produce_edge_integrals) {
      resident += 2U * edges3 + (super ? edges3 : 0U);
    }
    if (mode.produce_boundary_flux_integrals) {
      resident += 3U * gas_records;
      if (nonideal) {
        resident += 3U * energy_records;
      }
      if (super) {
        resident += 2U * energy_records;
      }
    }
    plan.level_resident_bytes += resident;
    plan.active_cells += cells;
    const std::array<std::size_t, 3> n{geo.nx, geo.ny, geo.nz};
    std::array<std::size_t, k_slot_lanes> counts{};
    std::size_t cursor = 0U;
    for (std::size_t lane_index = 0; lane_index < 9U; ++lane_index) {
      counts[cursor++] = (n[0] + 2U * g) * (n[1] + 2U * g) * (n[2] + 2U * g);
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      counts[cursor++] =
          (n[axis] + 1U) * (n[(axis + 1U) % 3U] + 2U * gt) * (n[(axis + 2U) % 3U] + 2U * gt);
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const std::size_t count =
          (n[axis] + 1U) * (n[(axis + 1U) % 3U] + 2U) * (n[(axis + 2U) % 3U] + 2U);
      counts[cursor++] = count;
      counts[cursor++] = count;
    }
    counts[cursor++] = geo.edgeCountX();
    counts[cursor++] = geo.edgeCountY();
    counts[cursor++] = geo.edgeCountZ();
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const std::size_t plane = n[(axis + 1U) % 3U] * n[(axis + 2U) % 3U];
      counts[cursor++] = MhdBoundaryFaceFluxes::k_values_per_face * plane;
      counts[cursor++] = MhdBoundaryFaceFluxes::k_values_per_face * plane;
    }
    for (std::size_t lane_index = 0; lane_index < 3U; ++lane_index) {
      counts[cursor++] = nonideal ? (n[0] + 2U) * (n[1] + 2U) * (n[2] + 2U) : 0U;
    }
    counts[cursor++] = nonideal ? cells : 0U;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const std::size_t plane = nonideal ? n[(axis + 1U) % 3U] * n[(axis + 2U) % 3U] : 0U;
      counts[cursor++] = plane;
      counts[cursor++] = plane;
    }
    for (std::size_t index = 0; index < k_slot_lanes; ++index) {
      slot_counts[index] = std::max(slot_counts[index], counts[index]);
    }
    max_cell_samples = std::max(max_cell_samples, patch.halo.cell_sample_count);
    max_face_samples = std::max(max_face_samples, patch.halo.face_sample_count);
    for (const Item& item : patch.halo.items) {
      max_item_values = std::max(max_item_values, itemValueCount(item));
    }
  }
  std::uint64_t per_slot = 0U;
  for (const std::size_t count : slot_counts) {
    per_slot += lane(count);
  }
  per_slot += lane(9U * static_cast<std::size_t>(max_cell_samples)) +
      lane(static_cast<std::size_t>(max_face_samples)) + max_item_values * sizeof(double);
  plan.batch_patches = static_cast<std::uint32_t>(m_slots.size());
  plan.batch_scratch_bytes = per_slot * m_slots.size() + max_item_values * sizeof(double);
  plan.topology_bytes = topologyCapacityBytes();
  // Communication: largest single exchange phase (send + receive staging, in
  // both the double and the byte representation).
  const int me = m_transport->rank();
  std::uint64_t max_phase = 0U;
  const auto phase_bytes = [&](const std::vector<Item>& items, std::size_t per_value_scale) {
    std::uint64_t bytes = 0U;
    for (const Item& item : items) {
      const bool src_local = m_topo[item.src].local != k_npos;
      const bool dst_local = m_topo[item.dst].local != k_npos;
      if (src_local == dst_local) {
        continue;
      }
      if (src_local && m_topo[item.dst].owner_rank == me) {
        continue;
      }
      bytes += 2U * per_value_scale * itemValueCount(item) * sizeof(double);
    }
    return bytes;
  };
  for (const std::vector<Item>* list :
       {&m_halo_items, &m_face_cross_items, &m_edge_cross_items, &m_record_cross_items,
        &m_face_fine_items, &m_record_fine_items, &m_energy_cross_items, &m_energy_fine_items}) {
    max_phase = std::max(max_phase, phase_bytes(*list, 1U));
  }
  for (const std::vector<Item>& list : m_edge_fine_items) {
    max_phase = std::max(max_phase, phase_bytes(list, 1U));
  }
  const auto world = static_cast<std::uint64_t>(m_transport->size());
  plan.communication_bytes =
      world > 1U ? max_phase + world * 2U * k_header_words * sizeof(std::uint64_t) : 0U;
  for (const PatchRuntime& patch : m_local) {
    const MhdPatchGeometry& g = patch.binding.geometry;
    plan.canonical_face_flux_bytes += 8U * g.faceCountTotal();
  }
  return plan;
}

void MhdLevelStepper::admitMemory(core::MemoryGovernor& governor, const MhdLevelAdvanceOptions& mode,
                                  const MhdResourceLimits& limits) {
  const MhdLevelMemoryPlan plan = planMemory(mode, limits);
  if (plan.active_cells != 0U &&
      plan.bytesPerActiveCell() > static_cast<double>(limits.scratch_bytes_per_active_cell)) {
    std::ostringstream stream;
    stream << "MHD phase+scratch per active cell (" << plan.bytesPerActiveCell()
           << " B: level-resident " << plan.level_resident_bytes << " B + batch scratch "
           << plan.batch_scratch_bytes << " B over " << plan.active_cells
           << " cells) exceeds the configured cap (MP02); reduce edge/boundary integral "
              "production, shrink batch_patches or use larger patches";
    MhdFailureContext failure = makeContext("mhd_level_stepper", stream.str());
    failure.requested_bytes = plan.phasePlusScratchBytes();
    failure.safe_limit_bytes = limits.scratch_bytes_per_active_cell * plan.active_cells;
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded, std::move(failure));
  }
  const std::uint64_t comm_limit = std::min(m_config.comm_budget_bytes, limits.comm_max_bytes_per_rank);
  if (plan.communication_bytes > comm_limit) {
    MhdFailureContext failure = makeContext(
        "mhd_level_stepper",
        "one MHD exchange phase exceeds the per-rank communication budget (MP03); repartition "
        "or select the explicit high-memory communication profile");
    failure.requested_bytes = plan.communication_bytes;
    failure.safe_limit_bytes = comm_limit;
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded, std::move(failure));
  }
  // Release a previous admission only after the new plan was validated.
  m_admission.release();
  m_comm_reservation.release();
  m_topology_reservation.release();
  MhdPatchMemoryEstimate aggregate{};
  aggregate.rollback_image_bytes = plan.level_resident_bytes;
  aggregate.primitive_cache_bytes = plan.batch_scratch_bytes;
  m_admission = MhdMemoryAdmission(governor, aggregate, 1U, false, "mhd_level_stepper");
  try {
    if (plan.topology_bytes != 0U) {
      m_topology_reservation = governor.reserve(core::MemoryClass::kPhaseResident,
                                                plan.topology_bytes, "mhd_level_stepper.topology");
      // Topology metadata already exists (built at construction).
      m_topology_reservation.commit();
    }
    if (plan.communication_bytes != 0U) {
      m_comm_reservation = governor.reserve(core::MemoryClass::kCommunication,
                                            plan.communication_bytes, "mhd_level_stepper.comm");
    }
  } catch (const core::MemoryAdmissionError& error) {
    m_admission.release();
    m_topology_reservation.release();
    m_comm_reservation.release();
    MhdFailureContext failure = makeContext(
        "mhd_level_stepper", std::string("MHD topology/communication admission refused: ") +
                                 error.what());
    failure.requested_bytes = plan.topology_bytes + plan.communication_bytes;
    failure.safe_limit_bytes = governor.policy().hard_limit_bytes;
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded, std::move(failure));
  }
  m_admitted_mode = mode;
  m_admitted_mode.timeline = {};
  m_admitted_mode.nonideal = nullptr;  // never dereferenced after admission
  m_admitted_nonideal = mode.nonideal != nullptr;
  m_admitted_super = m_admitted_nonideal && mode.nonideal->diffusionEnabled();
  m_admitted_plan = plan;
  m_admission_committed = false;
}

// ---------------------------------------------------------------------------
// Interval driver
// ---------------------------------------------------------------------------

MhdLevelAdvanceResult MhdLevelStepper::advance(const MhdStageContext& interval,
                                               const MhdLevelAdvanceOptions& options) {
  if (!(interval.dt_code > 0.0) || !std::isfinite(interval.dt_code) ||
      !(interval.scale_factor > 0.0) || !(interval.scale_factor_end > 0.0) ||
      !std::isfinite(interval.hubble_rate_code) || !std::isfinite(interval.hubble_rate_end_code)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_level_stepper", "invalid interval (dt, a or H)"));
  }
  if (!interval.comoving_coordinates &&
      (interval.scale_factor != 1.0 || interval.scale_factor_end != 1.0 ||
       interval.hubble_rate_code != 0.0 || interval.hubble_rate_end_code != 0.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_level_stepper",
                               "isolated (non-comoving) runs require a = 1 and H = 0 exactly"));
  }
  // Collective: an admission/allocation refusal on any rank stops every rank
  // before the first exchange (no rank may enter a stage alone).
  collective("allocation", [&]() {
    requireModeAdmitted(options);
    ensureAllocated(options);
  });
  m_nonideal = options.nonideal;
  m_track_integral = options.produce_edge_integrals;
  m_track_flux_integrals = options.produce_boundary_flux_integrals;
  MhdLevelAdvanceResult result{};
  result.dt_requested_code = interval.dt_code;
  const std::uint64_t bytes_sent_before = m_bytes_sent;
  const std::uint64_t bytes_received_before = m_bytes_received;
  const std::uint64_t phases_before = m_exchange_phases;
  // (1) Read-only validation BEFORE anything is written; (2) the rollback
  // image is the exact caller state; (3) normalization happens inside the
  // transaction.
  (void)validateCanonicalConsistency();
  result.budget_begin = levelBudget(interval.scale_factor);
  snapshotRollback();
  const bool constant_timeline = interval.scale_factor == interval.scale_factor_end &&
      interval.hubble_rate_code == 0.0 && interval.hubble_rate_end_code == 0.0;
  const bool internal_retry = m_config.retry_mode == MhdRetryMode::kInternalTimeline &&
      (constant_timeline || static_cast<bool>(options.timeline));
  MhdStageContext attempt = interval;
  for (std::int32_t retry = 0; retry <= m_policy.max_step_retries; ++retry) {
    RollbackGuard guard([this]() noexcept { restoreRollback(); });
    MhdStageCounters attempt_ideal{};
    MhdTransportStageDiagnostics attempt_transport{};
    const std::uint64_t stages_before = result.ideal_stage_evaluations_attempted;
    std::exception_ptr failure;
    MhdErrorCode code = MhdErrorCode::kInvalidField;
    try {
      reconcileFaces(false, nullptr);
      for (PatchRuntime& patch : m_local) {
        if (m_track_integral) {
          patch.integral.zero();
        }
        if (m_track_flux_integrals) {
          zeroRecords(patch.gas_integral);
          if (m_nonideal != nullptr) {
            zeroRecords(patch.energy_integral);
          }
        }
      }
      m_interval_rkl_stages = 0;
      m_interval_hall_subcycles = 0;
      if (m_nonideal != nullptr) {
        forecastNonIdeal(attempt, attempt_transport, result, true);
        const MhdStageContext first_half = subIntervalContext(attempt, 0.0, 0.5);
        const MhdStageContext second_half = subIntervalContext(attempt, 0.5, 1.0);
        diffusionHalf(first_half, attempt_transport, result);
        hallHalf(first_half, attempt_transport, result);
        idealSsprk2(attempt, attempt_ideal, result);
        hallHalf(second_half, attempt_transport, result);
        diffusionHalf(second_half, attempt_transport, result);
      } else {
        idealSsprk2(attempt, attempt_ideal, result);
      }
      validateAccepted(attempt.scale_factor_end, result);
      result.budget_end = levelBudget(attempt.scale_factor_end);
      guard.commit();
      result.accepted = true;
      result.dt_accepted_code = attempt.dt_code;
      result.retry_count = retry;
      result.ideal = attempt_ideal;
      result.transport = attempt_transport;
      result.ideal_attempted += attempt_ideal;
      result.transport_attempted += attempt_transport;
      result.ideal_stage_evaluations_accepted =
          result.ideal_stage_evaluations_attempted - stages_before;
      const double faces =
          static_cast<double>(std::max<std::uint64_t>(1U, result.ideal.faces_solved));
      result.fallback_fraction =
          static_cast<double>(result.ideal.hlle_fallback_faces + result.ideal.first_order_faces) /
          faces;
      result.fallback_fraction_gate_exceeded =
          result.fallback_fraction > m_policy.fallback_fraction_limit;
      result.scratch_capacity_bytes = ownedCapacityBytes();
      result.mpi_bytes_sent = m_bytes_sent - bytes_sent_before;
      result.mpi_bytes_received = m_bytes_received - bytes_received_before;
      result.exchange_phases = m_exchange_phases - phases_before;
      return result;
    } catch (const MhdError& error) {
      failure = std::current_exception();
      code = error.code();
      result.last_rejection = error.what();
    } catch (const std::bad_alloc&) {
      failure = std::current_exception();
      code = MhdErrorCode::kMemoryBudgetExceeded;
      result.last_rejection = "ungoverned allocation failure inside an MHD interval";
    } catch (const std::exception& error) {
      failure = std::current_exception();
      code = MhdErrorCode::kInvalidField;
      result.last_rejection = error.what();
    }
    // Atomic rejection: every local patch returns to the interval-start image
    // (the guard runs on every exit path below). Attempted counters only.
    result.ideal_attempted += attempt_ideal;
    result.transport_attempted += attempt_transport;
    ++result.rejected_attempts;
    if (!retryable(code) || retry == m_policy.max_step_retries) {
      guard.commit();
      restoreRollback();
      std::rethrow_exception(failure);
    }
    const double shorter = 0.5 * attempt.dt_code;
    if (!internal_retry) {
      // The scheduler owns the timeline: return the rejection; nothing of
      // this interval is committed.
      guard.commit();
      restoreRollback();
      result.accepted = false;
      result.retry_count = retry + 1;
      result.dt_suggested_code = shorter;
      result.scratch_capacity_bytes = ownedCapacityBytes();
      result.mpi_bytes_sent = m_bytes_sent - bytes_sent_before;
      result.mpi_bytes_received = m_bytes_received - bytes_received_before;
      result.exchange_phases = m_exchange_phases - phases_before;
      return result;
    }
    attempt.dt_code = shorter;
    if (options.timeline) {
      const std::pair<double, double> end = options.timeline(attempt.time_code + shorter);
      attempt.scale_factor_end = end.first;
      attempt.hubble_rate_end_code = end.second;
    }
    result.retry_count = retry + 1;
  }
  return result;
}

double MhdLevelStepper::estimateStableTimeStep(const MhdStageContext& interval,
                                               const MhdLevelAdvanceOptions& options) {
  double dt = std::numeric_limits<double>::infinity();
  collective("estimate.ideal", [&]() {
    for (const PatchRuntime& patch : m_local) {
      const MhdFaceFluxView view{&patch.binding.geometry,
                                 patch.binding.faces.phi_x,
                                 patch.binding.faces.phi_y,
                                 patch.binding.faces.phi_z,
                                 patch.binding.patch_id,
                                 patch.binding.geometry_generation,
                                 patch.binding.level};
      dt = std::min(dt, m_operator.estimateIdealTimeStep(patch.binding.geometry,
                                                         *patch.binding.gas, view,
                                                         interval.scale_factor, m_policy.cfl_ideal));
    }
  });
  if (options.nonideal != nullptr) {
    collective("estimate.allocation", [&]() {
      requireModeAdmitted(options);
      ensureAllocated(options);
    });
    m_nonideal = options.nonideal;
    MhdTransportStageDiagnostics scratch_diagnostics{};
    MhdLevelAdvanceResult scratch{};
    forecastNonIdeal(interval, scratch_diagnostics, scratch, false);
    const MhdStiffnessLimits& limits = m_nonideal->limits();
    if (m_nonideal->hallEnabled() && std::isfinite(m_forecast_dt_hall)) {
      // 2 * ceil(dt/2 / dt_H) <= cap  <=>  dt <= cap * dt_H (cap even part).
      dt = std::min(dt, static_cast<double>(limits.hall_subcycle_cost_cap / 2) * 2.0 *
                            m_forecast_dt_hall);
    }
    if (m_nonideal->diffusionEnabled() && std::isfinite(m_forecast_dt_diffusion)) {
      // Largest per-half ratio whose RKL2 plan fits half of the stage cap.
      const std::int32_t stages_per_half = std::max(2, limits.rkl_stage_cost_cap / 2);
      const std::int32_t cycles = std::max(1, stages_per_half / limits.rkl_max_stages);
      const std::int32_t per_cycle = std::min(limits.rkl_max_stages, stages_per_half);
      const double ratio = static_cast<double>(cycles) * rkl2StabilityFactor(per_cycle);
      dt = std::min(dt, 2.0 * ratio * m_forecast_dt_diffusion);
    }
  }
  return m_transport->size() > 1 ? m_transport->allreduceMinDouble(dt) : dt;
}

MhdConservationBudget MhdLevelStepper::levelBudget(double scale_factor) {
  MhdConservationBudget total{};
  collective("budget", [&]() {
    for (const PatchRuntime& patch : m_local) {
      const MhdFaceFluxView view{&patch.binding.geometry,
                                 patch.binding.faces.phi_x,
                                 patch.binding.faces.phi_y,
                                 patch.binding.faces.phi_z,
                                 patch.binding.patch_id,
                                 patch.binding.geometry_generation,
                                 patch.binding.level};
      const MhdConservationBudget budget =
          m_operator.patchBudget(patch.binding.geometry, *patch.binding.gas, view, scale_factor);
      total.mass += budget.mass;
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        total.momentum[axis] += budget.momentum[axis];
      }
      total.energy_thermal += budget.energy_thermal;
      total.energy_kinetic += budget.energy_kinetic;
      total.energy_magnetic += budget.energy_magnetic;
      total.energy_total += budget.energy_total;
      total.divergence_max_normalized =
          std::max(total.divergence_max_normalized, budget.divergence_max_normalized);
    }
  });
  if (m_transport->size() > 1) {
    std::array<double, 8> sums{total.mass,           total.momentum[0],     total.momentum[1],
                               total.momentum[2],    total.energy_thermal, total.energy_kinetic,
                               total.energy_magnetic, total.energy_total};
    m_transport->allreduceSumDoubles(sums);
    total.mass = sums[0];
    total.momentum = {sums[1], sums[2], sums[3]};
    total.energy_thermal = sums[4];
    total.energy_kinetic = sums[5];
    total.energy_magnetic = sums[6];
    total.energy_total = sums[7];
    total.divergence_max_normalized =
        m_transport->allreduceMaxDouble(total.divergence_max_normalized);
  }
  return total;
}

}  // namespace cosmosim::mhd
