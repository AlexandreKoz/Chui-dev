#include "cosmosim/mhd/mhd_memory.hpp"

#include <algorithm>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <utility>

namespace cosmosim::mhd {
namespace {

[[nodiscard]] std::uint64_t lane(std::uint64_t count) {
  if (count > std::numeric_limits<std::size_t>::max()) {
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded,
                   MhdFailureContext{.owner = "mhd_memory", .detail = "lane count overflow"});
  }
  return mhdAlignedLaneBytes(static_cast<std::size_t>(count));
}

[[nodiscard]] std::uint64_t checkedProduct(std::initializer_list<std::uint64_t> factors) {
  std::uint64_t product = 1U;
  for (const std::uint64_t factor : factors) {
    if (factor != 0U && product > std::numeric_limits<std::uint64_t>::max() / factor) {
      throw MhdError(MhdErrorCode::kMemoryBudgetExceeded,
                     MhdFailureContext{.owner = "mhd_memory", .detail = "count product overflow"});
    }
    product *= factor;
  }
  return product;
}

[[nodiscard]] std::uint64_t checkedSum(std::uint64_t lhs, std::uint64_t rhs) {
  if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded,
                   MhdFailureContext{.owner = "mhd_memory", .detail = "byte sum overflow"});
  }
  return lhs + rhs;
}

}  // namespace

std::uint64_t mhdFaceCount(const MhdPatchGeometry& geometry) {
  const std::uint64_t nx = geometry.nx;
  const std::uint64_t ny = geometry.ny;
  const std::uint64_t nz = geometry.nz;
  return checkedSum(checkedSum(checkedProduct({nx + 1U, ny, nz}), checkedProduct({nx, ny + 1U, nz})),
                    checkedProduct({nx, ny, nz + 1U}));
}

std::uint64_t mhdEdgeCount(const MhdPatchGeometry& geometry) {
  const std::uint64_t nx = geometry.nx;
  const std::uint64_t ny = geometry.ny;
  const std::uint64_t nz = geometry.nz;
  return checkedSum(
      checkedSum(checkedProduct({nx, ny + 1U, nz + 1U}), checkedProduct({nx + 1U, ny, nz + 1U})),
      checkedProduct({nx + 1U, ny + 1U, nz}));
}

MhdPatchMemoryEstimate estimatePatchMemory(const MhdMemoryRequest& request) {
  const MhdPatchGeometry& geometry = request.geometry;
  if (!geometry.valid()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_memory", .detail = "invalid patch geometry"});
  }
  const std::array<std::uint64_t, 3> n{geometry.nx, geometry.ny, geometry.nz};
  const std::uint64_t g = mhdGhostCells(request.reconstruction);
  const std::uint64_t gt = k_mhd_staggered_tangential_ghost;
  const std::uint64_t cells = checkedProduct({n[0], n[1], n[2]});
  const std::array<std::uint64_t, 3> faces{checkedProduct({n[0] + 1U, n[1], n[2]}),
                                           checkedProduct({n[0], n[1] + 1U, n[2]}),
                                           checkedProduct({n[0], n[1], n[2] + 1U})};
  const std::array<std::uint64_t, 3> edges{checkedProduct({n[0], n[1] + 1U, n[2] + 1U}),
                                           checkedProduct({n[0] + 1U, n[1], n[2] + 1U}),
                                           checkedProduct({n[0] + 1U, n[1] + 1U, n[2]})};
  const auto edge_storage = [&]() {
    return lane(edges[0]) + lane(edges[1]) + lane(edges[2]);
  };
  const auto face_lanes = [&]() {
    return lane(faces[0]) + lane(faces[1]) + lane(faces[2]);
  };
  MhdPatchMemoryEstimate estimate{};
  estimate.canonical_face_flux_bytes = face_lanes();
  estimate.primitive_cache_bytes =
      9U * lane(checkedProduct({n[0] + 2U * g, n[1] + 2U * g, n[2] + 2U * g}));
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::uint64_t t1 = n[(axis + 1U) % 3U];
    const std::uint64_t t2 = n[(axis + 2U) % 3U];
    estimate.staggered_cache_bytes += lane(checkedProduct({n[axis] + 1U, t1 + 2U * gt, t2 + 2U * gt}));
    estimate.face_signal_bytes += 2U * lane(checkedProduct({n[axis] + 1U, t1 + 2U, t2 + 2U}));
    estimate.boundary_flux_bytes += 2U * lane(checkedProduct({MhdBoundaryFaceFluxes::k_values_per_face, t1, t2}));
  }
  estimate.stage_emf_bytes = edge_storage();
  estimate.rollback_image_bytes = MhdGasStateSoa::k_lane_count * lane(cells) + face_lanes();
  estimate.stage0_image_bytes = MhdGasStateSoa::k_lane_count * lane(cells) + face_lanes();
  if (request.nonideal) {
    std::uint64_t bytes = 3U * lane(checkedProduct({n[0] + 2U, n[1] + 2U, n[2] + 2U}));
    bytes += lane(cells);  // energy_rate_stage
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      bytes += 2U * lane(checkedProduct({n[(axis + 1U) % 3U], n[(axis + 2U) % 3U]}));
    }
    if (request.super_stepping) {
      bytes += face_lanes() + 2U * lane(cells) + edge_storage();
    }
    estimate.nonideal_bytes = bytes;
  }
  if (request.edge_integrals) {
    estimate.edge_integral_bytes = 3U * edge_storage();
  }
  return estimate;
}

double scratchBytesPerActiveCell(const MhdMemoryRequest& request, std::size_t patch_count) {
  const MhdPatchMemoryEstimate estimate = estimatePatchMemory(request);
  const double cells = static_cast<double>(request.geometry.cellCount()) *
      static_cast<double>(std::max<std::size_t>(patch_count, 1U));
  return static_cast<double>(estimate.phasePlusScratchBytes()) *
      static_cast<double>(std::max<std::size_t>(patch_count, 1U)) / cells;
}

std::size_t selectMhdBatchPatches(const MhdMemoryRequest& request,
                                  std::uint64_t available_bytes,
                                  const MhdResourceLimits& limits) {
  const MhdPatchMemoryEstimate estimate = estimatePatchMemory(request);
  const std::uint64_t per_patch = estimate.phasePlusScratchBytes();
  const std::uint64_t cap_per_patch =
      checkedProduct({limits.scratch_bytes_per_active_cell, request.geometry.cellCount()});
  if (per_patch > cap_per_patch || per_patch > available_bytes || per_patch == 0U) {
    return 0U;
  }
  const std::uint64_t by_memory = available_bytes / per_patch;
  return static_cast<std::size_t>(
      std::min<std::uint64_t>(by_memory, std::max<std::uint32_t>(1U, limits.active_batch_target_patches)));
}

void requireHostEnvelope(std::span<const std::uint64_t> rank_hard_limits_on_host,
                         const MhdResourceLimits& limits,
                         bool performance_qualified) {
  std::uint64_t total = 0U;
  for (const std::uint64_t limit : rank_hard_limits_on_host) {
    if (limit == 0U && performance_qualified) {
      throw MhdError(MhdErrorCode::kMemoryBudgetExceeded,
                     MhdFailureContext{.owner = "mhd_memory",
                                       .detail = "performance-qualified MHD run has an unlimited "
                                                 "(hard_limit_bytes=0) rank governor (MP11)"});
    }
    total = checkedSum(total, limit);
  }
  if (total > limits.host_safe_bytes) {
    MhdFailureContext context{};
    context.owner = "mhd_memory";
    context.detail = "sum of rank hard limits on one host exceeds the whole-host safe ceiling";
    context.requested_bytes = total;
    context.safe_limit_bytes = limits.host_safe_bytes;
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded, std::move(context));
  }
}

MhdMemoryAdmission::MhdMemoryAdmission(core::MemoryGovernor& governor,
                                       const MhdPatchMemoryEstimate& per_patch,
                                       std::size_t patch_count,
                                       bool charge_canonical_faces,
                                       std::string owner)
    : m_owner(std::move(owner)) {
  const std::uint64_t count = patch_count;
  m_phase_bytes = checkedProduct({per_patch.phaseResidentBytes(), count});
  m_scratch_bytes = checkedProduct({per_patch.scratchArenaBytes(), count});
  m_canonical_bytes =
      charge_canonical_faces ? checkedProduct({per_patch.canonical_face_flux_bytes, count}) : 0U;
  try {
    if (m_canonical_bytes != 0U) {
      m_canonical = governor.reserve(core::MemoryClass::kCanonicalPersistent, m_canonical_bytes,
                                     m_owner + ".face_flux");
    }
    m_phase = governor.reserve(core::MemoryClass::kPhaseResident, m_phase_bytes,
                               m_owner + ".rollback_stage0");
    m_scratch = governor.reserve(core::MemoryClass::kScratchArena, m_scratch_bytes,
                                 m_owner + ".stage_scratch");
  } catch (const core::MemoryAdmissionError& error) {
    m_canonical.release();
    m_phase.release();
    m_scratch.release();
    MhdFailureContext context{};
    context.owner = m_owner;
    context.detail = std::string("MHD workspace admission refused: ") + error.what();
    context.requested_bytes = m_phase_bytes + m_scratch_bytes + m_canonical_bytes;
    context.safe_limit_bytes = governor.policy().hard_limit_bytes;
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded, std::move(context));
  }
  m_admitted = true;
}

void MhdMemoryAdmission::commit(std::uint64_t actual_phase_plus_scratch_bytes) {
  if (!m_admitted) {
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded,
                   MhdFailureContext{.owner = m_owner, .detail = "commit without admission"});
  }
  if (actual_phase_plus_scratch_bytes > m_phase_bytes + m_scratch_bytes) {
    MhdFailureContext context{};
    context.owner = m_owner;
    context.detail = "actual retained MHD scratch capacity exceeds the admitted estimate "
                     "(unexpected allocation growth)";
    context.requested_bytes = actual_phase_plus_scratch_bytes;
    context.admitted_bytes = m_phase_bytes + m_scratch_bytes;
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded, std::move(context));
  }
  for (core::MemoryReservation* reservation : {&m_canonical, &m_phase, &m_scratch}) {
    if (reservation->pending()) {
      reservation->commit();
    }
  }
}

void MhdMemoryAdmission::release() noexcept {
  m_canonical.release();
  m_phase.release();
  m_scratch.release();
  m_admitted = false;
}

}  // namespace cosmosim::mhd
