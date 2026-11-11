#include "workflows/internal/amr_geometry_assurance.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/amr/amr_geometry_coverage.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"

namespace cosmosim::workflows::internal {
namespace {

constexpr std::uint64_t k_fnv_offset = 1469598103934665603ULL;
constexpr std::uint64_t k_fnv_prime = 1099511628211ULL;

void hashWord(std::uint64_t& hash, std::uint64_t word) noexcept {
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (word >> (8 * byte)) & 0xFFULL;
    hash *= k_fnv_prime;
  }
}

[[nodiscard]] amr::AmrPatchGeometryRecord recordForPatch(
    const core::SimulationState& state, std::size_t patch) noexcept {
  const core::PatchSoa& p = state.patches;
  amr::AmrPatchGeometryRecord record;
  record.patch_id = p.patch_id[patch];
  record.level = p.level[patch];
  record.owning_rank = p.owning_rank[patch];
  record.cell_count = p.cell_count[patch];
  record.origin_comoving = {p.origin_x_comoving[patch], p.origin_y_comoving[patch], p.origin_z_comoving[patch]};
  record.extent_comoving = {p.extent_x_comoving[patch], p.extent_y_comoving[patch], p.extent_z_comoving[patch]};
  record.cell_dims = {p.cell_dim_x[patch], p.cell_dim_y[patch], p.cell_dim_z[patch]};
  return record;
}

[[nodiscard]] std::string phaseLabel(std::string_view phase, std::string_view step) {
  return std::string(phase) + " (" + std::string(step) + ")";
}

[[nodiscard]] std::string formatDouble(double value) { return std::to_string(value); }

}  // namespace

std::uint64_t amrPopulatedGeometryFingerprint(
    const core::SimulationState& state, const AmrGeometryAssurancePolicy& policy) noexcept {
  std::uint64_t hash = k_fnv_offset;
  hashWord(hash, policy.require_explicit_geometry ? 1U : 0U);
  hashWord(hash, policy.periodic ? 1U : 0U);
  for (const double extent : policy.periodic_box_extent_comoving) hashWord(hash, std::bit_cast<std::uint64_t>(extent));
  const core::PatchSoa& p = state.patches;
  const std::size_t count = p.isConsistent() ? p.size() : 0U;
  hashWord(hash, p.isConsistent() ? 1U : 0U);
  hashWord(hash, static_cast<std::uint64_t>(p.size()));
  for (std::size_t patch = 0; patch < count; ++patch) {
    if (p.cell_count[patch] == 0U) continue;
    hashWord(hash, p.patch_id[patch]);
    hashWord(hash, static_cast<std::uint64_t>(static_cast<std::uint32_t>(p.level[patch])));
    hashWord(hash, (static_cast<std::uint64_t>(p.owning_rank[patch]) << 32U) | p.cell_count[patch]);
    hashWord(hash, std::bit_cast<std::uint64_t>(p.origin_x_comoving[patch]));
    hashWord(hash, std::bit_cast<std::uint64_t>(p.origin_y_comoving[patch]));
    hashWord(hash, std::bit_cast<std::uint64_t>(p.origin_z_comoving[patch]));
    hashWord(hash, std::bit_cast<std::uint64_t>(p.extent_x_comoving[patch]));
    hashWord(hash, std::bit_cast<std::uint64_t>(p.extent_y_comoving[patch]));
    hashWord(hash, std::bit_cast<std::uint64_t>(p.extent_z_comoving[patch]));
    hashWord(hash, (static_cast<std::uint64_t>(p.cell_dim_x[patch]) << 32U) |
                       (static_cast<std::uint64_t>(p.cell_dim_y[patch]) << 16U) | p.cell_dim_z[patch]);
  }
  return hash;
}

AmrGeometryAssuranceResult AmrGeometryAssurance::assure(
    const core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    core::MemoryGovernor* memory_governor,
    const AmrGeometryAssurancePolicy& policy,
    std::string_view phase) {
  // C1: cache gate. A rank whose local geometry is unchanged and whose cache
  // holds an accepted proof votes 0; any change anywhere re-runs the proof on
  // every rank. The cache is only ever set after a collectively accepted
  // proof, so all ranks agree on its validity.
  const std::uint64_t fingerprint = amrPopulatedGeometryFingerprint(state, policy);
  const std::uint64_t local_changed = (!m_valid || fingerprint != m_fingerprint) ? 1U : 0U;
  if (mpi_context.allreduceMaxUint64(local_changed) == 0U) {
    AmrGeometryAssuranceResult cached = m_result;
    cached.revalidated = false;
    return cached;
  }
  m_valid = false;

  const double inf = std::numeric_limits<double>::infinity();
  std::array<double, 3> local_lower{inf, inf, inf};
  std::array<double, 3> local_negative_upper{inf, inf, inf};
  std::array<double, 3> local_finest{inf, inf, inf};
  // Reservations are declared before the buffers they cover so the buffers
  // are destroyed (memory returned) before the lease is released.
  core::MemoryReservation scratch_reservation;
  core::MemoryReservation gather_reservation;
  std::vector<amr::AmrPatchGeometryRecord> records;
  std::vector<amr::AmrSignedLatticeCorner> corners;

  // V1: local geometry admission.
  std::exception_ptr failure;
  try {
    if (!state.patches.isConsistent()) {
      throw std::runtime_error("AMR patch SoA lanes are inconsistent (lane sizes differ)");
    }
    std::size_t populated = 0U;
    for (std::size_t patch = 0; patch < state.patches.size(); ++patch) {
      if (state.patches.cell_count[patch] != 0U) ++populated;
    }
    if (memory_governor != nullptr && populated != 0U) {
      const std::uint64_t bytes = static_cast<std::uint64_t>(core::checkedSizeMultiply(
          populated, sizeof(amr::AmrPatchGeometryRecord) + 8U * sizeof(amr::AmrSignedLatticeCorner),
          "AMR geometry proof local scratch bytes"));
      scratch_reservation = memory_governor->reserve(
          core::MemoryClass::kScratchArena, bytes, "hydro.amr.geometry_proof.local_scratch");
      scratch_reservation.commit();
    }
    records.reserve(populated);
    for (std::size_t patch = 0; patch < state.patches.size(); ++patch) {
      if (state.patches.cell_count[patch] == 0U) continue;  // inactive parent / metadata-only
      const amr::AmrPatchGeometryRecord record = recordForPatch(state, patch);
      if (amr::isLegacyUndeclaredAmrPatchGeometry(record)) {
        if (policy.require_explicit_geometry) {
          throw std::runtime_error(
              "AMR patch " + std::to_string(record.patch_id) + " (rank " + std::to_string(record.owning_rank) +
              ", level " + std::to_string(record.level) +
              ") is populated but declares no geometry (legacy all-zero descriptor); production AMR hydro "
              "requires explicit geometry on every populated leaf");
        }
        continue;  // non-AMR legacy row: not part of any geometric tiling
      }
      const std::string issue = amr::diagnoseAmrPatchGeometry(record);
      if (!issue.empty()) throw std::runtime_error(issue);
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        local_lower[axis] = std::min(local_lower[axis], record.origin_comoving[axis]);
        local_negative_upper[axis] =
            std::min(local_negative_upper[axis], -(record.origin_comoving[axis] + record.extent_comoving[axis]));
        local_finest[axis] = std::min(
            local_finest[axis], record.extent_comoving[axis] / static_cast<double>(record.cell_dims[axis]));
      }
      records.push_back(record);
    }
  } catch (...) {
    failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(failure, phaseLabel(phase, "AMR patch geometry admission"));

  // C2: populated geometric leaves anywhere?
  const std::uint64_t global_populated =
      mpi_context.allreduceSumUint64(static_cast<std::uint64_t>(records.size()));
  if (global_populated == 0U) {
    m_result = AmrGeometryAssuranceResult{};
    m_result.revalidated = true;
    m_fingerprint = fingerprint;
    m_valid = true;
    return m_result;
  }

  // C3: global bounding box and finest populated cell width.
  std::array<double, 3> lower{};
  std::array<double, 3> upper{};
  std::array<double, 3> finest{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    lower[axis] = mpi_context.allreduceMinDouble(local_lower[axis]);
    upper[axis] = -mpi_context.allreduceMinDouble(local_negative_upper[axis]);
    finest[axis] = mpi_context.allreduceMinDouble(local_finest[axis]);
  }

  // V2: lattice, periodic-box agreement, local corners + local cancellation.
  amr::AmrCoverageLattice lattice;
  double local_volume = 0.0;
  failure = nullptr;
  try {
    lattice = amr::makeAmrCoverageLattice(lower, upper, finest);
    if (policy.periodic) {
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        const double length = policy.periodic_box_extent_comoving[axis];
        if (!(length > 0.0) || !std::isfinite(length)) continue;  // non-periodic axis of the frame
        const double unit = lattice.unit_comoving[axis];
        const double tol = amr::amrCoverageSnapTolerance(static_cast<double>(lattice.cells[axis]));
        if (std::abs(lower[axis] / unit) > tol || std::abs((upper[axis] - length) / unit) > tol) {
          throw std::runtime_error(
              "periodic AMR hydro requires the populated leaf tiling to cover the periodic box [0, " +
              formatDouble(length) + ") on axis " + std::to_string(axis) + "; the leaves span [" +
              formatDouble(lower[axis]) + ", " + formatDouble(upper[axis]) +
              "). A periodic sub-box hydro domain has no defined wrap and is not supported");
        }
      }
    }
    corners.reserve(records.size() * 8U);
    for (const amr::AmrPatchGeometryRecord& record : records) {
      amr::appendAmrPatchSignedCorners(lattice, record, corners);
      local_volume += static_cast<double>(amr::amrPatchLatticeVolume(lattice, record));
    }
    records.clear();
    records.shrink_to_fit();
    amr::aggregateAmrSignedCorners(corners);
  } catch (...) {
    failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(failure, phaseLabel(phase, "AMR lattice representability"));

  // C4: residual size and (diagnostic) summed lattice volume.
  const std::uint64_t global_residual =
      mpi_context.allreduceSumUint64(static_cast<std::uint64_t>(corners.size()));
  const double global_volume = mpi_context.allreduceSumDouble(local_volume);

  // V3: gather admission. Root holds the gathered bytes and their decoded
  // copy; other ranks hold only their own (already reserved) residual.
  const int root = 0;
  const bool is_root = mpi_context.worldRank() == root;
  failure = nullptr;
  try {
    if (memory_governor != nullptr && is_root && global_residual != 0U) {
      const std::uint64_t bytes = static_cast<std::uint64_t>(core::checkedSizeMultiply(
          core::checkedIntegralNarrow<std::size_t>(global_residual, "AMR geometry proof gathered corner count"),
          2U * sizeof(amr::AmrSignedLatticeCorner), "AMR geometry proof gather bytes"));
      gather_reservation = memory_governor->reserve(
          core::MemoryClass::kCommunication, bytes, "hydro.amr.geometry_proof.root_gather");
      gather_reservation.commit();
    }
  } catch (...) {
    failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(failure, phaseLabel(phase, "AMR coverage gather admission"));

  // C5: bounded gather of the residual measure.
  const std::span<const std::uint8_t> local_bytes(
      reinterpret_cast<const std::uint8_t*>(corners.data()), corners.size() * sizeof(amr::AmrSignedLatticeCorner));
  std::vector<std::uint8_t> gathered = mpi_context.gatherBytesToRoot(local_bytes, root);

  // Root: exact verdict. Any root-side failure becomes the verdict text so the
  // broadcast below always happens (no rank is left blocked).
  std::string verdict_text;
  if (is_root) {
    try {
      if (gathered.size() % sizeof(amr::AmrSignedLatticeCorner) != 0U) {
        throw std::runtime_error("AMR coverage gather returned a truncated corner record");
      }
      corners.resize(gathered.size() / sizeof(amr::AmrSignedLatticeCorner));
      if (!gathered.empty()) std::memcpy(corners.data(), gathered.data(), gathered.size());
      gathered.clear();
      gathered.shrink_to_fit();
      amr::aggregateAmrSignedCorners(corners);
      const amr::AmrTilingVerdict verdict =
          amr::verifyAmrExactTiling(lattice, corners, static_cast<long double>(global_volume));
      if (!verdict.exact_tiling) verdict_text = verdict.diagnostic;
    } catch (const std::exception& error) {
      verdict_text = std::string("AMR coverage proof failed on the root rank: ") + error.what();
    } catch (...) {
      verdict_text = "AMR coverage proof failed on the root rank with a non-standard exception";
    }
  }
  gathered.clear();
  corners.clear();
  corners.shrink_to_fit();

  // C6: identical verdict on every rank.
  const std::vector<std::uint8_t> verdict_bytes = mpi_context.broadcastBytesFromRoot(
      std::span<const std::uint8_t>(reinterpret_cast<const std::uint8_t*>(verdict_text.data()), verdict_text.size()),
      root);
  if (!verdict_bytes.empty()) {
    throw std::runtime_error(
        std::string(phase) + ": " + std::string(verdict_bytes.begin(), verdict_bytes.end()));
  }

  m_result = AmrGeometryAssuranceResult{};
  m_result.declared = true;
  m_result.populated_patches = global_populated;
  m_result.gathered_residual_points = global_residual;
  m_result.revalidated = true;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    m_result.origin_comoving[axis] = lower[axis];
    m_result.extent_comoving[axis] = upper[axis] - lower[axis];
  }
  m_fingerprint = fingerprint;
  m_valid = true;
  return m_result;
}

}  // namespace cosmosim::workflows::internal
