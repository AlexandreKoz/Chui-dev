#include "cosmosim/mhd/mhd_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>

namespace cosmosim::mhd {
namespace {

constexpr std::uint64_t k_fnv_offset_basis = 1469598103934665603ULL;
constexpr std::uint64_t k_fnv_prime = 1099511628211ULL;

void hashMix(std::uint64_t& hash, std::uint64_t value) noexcept {
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & 0xFFULL;
    hash *= k_fnv_prime;
  }
}

std::uint64_t signedToUnsigned(std::int64_t value) noexcept {
  return static_cast<std::uint64_t>(value);
}

}  // namespace

// ---------------------------------------------------------------------------
// MhdEdgeEmfStorage
// ---------------------------------------------------------------------------

MhdEdgeEmfStorage::MhdEdgeEmfStorage(const MhdPatchGeometry& geometry) {
  resize(geometry);
}

void MhdEdgeEmfStorage::resize(const MhdPatchGeometry& geometry) {
  m_geometry = geometry;
  m_emf_x.resize(geometry.edgeCountX(), 0.0);
  m_emf_y.resize(geometry.edgeCountY(), 0.0);
  m_emf_z.resize(geometry.edgeCountZ(), 0.0);
}

void MhdEdgeEmfStorage::zero() noexcept {
  m_emf_x.assignFill(0.0);
  m_emf_y.assignFill(0.0);
  m_emf_z.assignFill(0.0);
}

// ---------------------------------------------------------------------------
// MhdCellFieldCache (ephemeral derived lane)
// ---------------------------------------------------------------------------

void MhdCellFieldCache::refresh(const MhdPatchGeometry& patch_geometry,
                                std::span<const double> phi_x,
                                std::span<const double> phi_y,
                                std::span<const double> phi_z) {
  geometry = patch_geometry;
  const std::size_t cells = patch_geometry.cellCount();
  bx_comoving.resize(cells, 0.0);
  by_comoving.resize(cells, 0.0);
  bz_comoving.resize(cells, 0.0);
  if (phi_x.size() != patch_geometry.faceCountX() ||
      phi_y.size() != patch_geometry.faceCountY() ||
      phi_z.size() != patch_geometry.faceCountZ()) {
    return;
  }
  const double area_x = patch_geometry.faceAreaX();
  const double area_y = patch_geometry.faceAreaY();
  const double area_z = patch_geometry.faceAreaZ();
  const double inv_area_x = area_x > 0.0 ? 1.0 / area_x : 0.0;
  const double inv_area_y = area_y > 0.0 ? 1.0 / area_y : 0.0;
  const double inv_area_z = area_z > 0.0 ? 1.0 / area_z : 0.0;
  for (std::size_t k = 0; k < patch_geometry.nz; ++k) {
    for (std::size_t j = 0; j < patch_geometry.ny; ++j) {
      for (std::size_t i = 0; i < patch_geometry.nx; ++i) {
        const std::size_t cell = patch_geometry.cellIndex(i, j, k);
        bx_comoving[cell] = 0.5 * (phi_x[patch_geometry.faceIndexX(i, j, k)] +
                                   phi_x[patch_geometry.faceIndexX(i + 1U, j, k)]) *
            inv_area_x;
        by_comoving[cell] = 0.5 * (phi_y[patch_geometry.faceIndexY(i, j, k)] +
                                   phi_y[patch_geometry.faceIndexY(i, j + 1U, k)]) *
            inv_area_y;
        bz_comoving[cell] = 0.5 * (phi_z[patch_geometry.faceIndexZ(i, j, k)] +
                                   phi_z[patch_geometry.faceIndexZ(i, j, k + 1U)]) *
            inv_area_z;
      }
    }
  }
}

std::array<double, 3> MhdCellFieldCache::field(std::size_t i, std::size_t j,
                                               std::size_t k) const {
  const std::size_t cell = geometry.cellIndex(i, j, k);
  if (cell >= bx_comoving.size()) {
    return {0.0, 0.0, 0.0};
  }
  return {bx_comoving[cell], by_comoving[cell], bz_comoving[cell]};
}

// ---------------------------------------------------------------------------
// Discrete divergence of face fluxes
// ---------------------------------------------------------------------------

std::vector<double> faceFluxDivergence(const MhdPatchGeometry& geometry,
                                       std::span<const double> phi_x,
                                       std::span<const double> phi_y,
                                       std::span<const double> phi_z) {
  std::vector<double> divergence(geometry.cellCount(), 0.0);
  if (phi_x.size() != geometry.faceCountX() ||
      phi_y.size() != geometry.faceCountY() ||
      phi_z.size() != geometry.faceCountZ()) {
    return divergence;
  }
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const double dx_term = phi_x[geometry.faceIndexX(i + 1U, j, k)] -
            phi_x[geometry.faceIndexX(i, j, k)];
        const double dy_term = phi_y[geometry.faceIndexY(i, j + 1U, k)] -
            phi_y[geometry.faceIndexY(i, j, k)];
        const double dz_term = phi_z[geometry.faceIndexZ(i, j, k + 1U)] -
            phi_z[geometry.faceIndexZ(i, j, k)];
        divergence[geometry.cellIndex(i, j, k)] = dx_term + dy_term + dz_term;
      }
    }
  }
  return divergence;
}

double maxNormalizedFaceDivergence(const MhdPatchGeometry& geometry,
                                   std::span<const double> phi_x,
                                   std::span<const double> phi_y,
                                   std::span<const double> phi_z,
                                   double normalizer) {
  const std::vector<double> divergence =
      faceFluxDivergence(geometry, phi_x, phi_y, phi_z);
  double max_abs = 0.0;
  for (const double value : divergence) {
    max_abs = std::max(max_abs, std::abs(value));
  }
  if (!(normalizer > 0.0) || !std::isfinite(normalizer)) {
    return max_abs > 0.0 ? std::numeric_limits<double>::infinity() : 0.0;
  }
  return max_abs / normalizer;
}

double faceFluxMagnitudeSum(std::span<const double> phi_x,
                            std::span<const double> phi_y,
                            std::span<const double> phi_z) {
  double sum = 0.0;
  for (const double value : phi_x) {
    sum += std::abs(value);
  }
  for (const double value : phi_y) {
    sum += std::abs(value);
  }
  for (const double value : phi_z) {
    sum += std::abs(value);
  }
  return sum;
}

// Domain interior volume V = sum_cell dx*dy*dz. A face flux Phi has units of
// B (in the divergence normalizer) * area, so D(Phi)_i = sum over the six
// oriented faces of the stored flux. The discrete identity D . C = 0 holds on
// the index lattice; no periodic wrap is applied here because wrap identity is
// a boundary-condition property owned by the caller that exchanged the faces.
void applyConstrainedTransportUpdate(const MhdPatchGeometry& geometry,
                                     std::span<double> phi_x,
                                     std::span<double> phi_y,
                                     std::span<double> phi_z,
                                     std::span<const double> emf_x,
                                     std::span<const double> emf_y,
                                     std::span<const double> emf_z,
                                     double dt_code) {
  if (phi_x.size() != geometry.faceCountX() ||
      phi_y.size() != geometry.faceCountY() ||
      phi_z.size() != geometry.faceCountZ() ||
      emf_x.size() != geometry.edgeCountX() ||
      emf_y.size() != geometry.edgeCountY() ||
      emf_z.size() != geometry.edgeCountZ()) {
    return;
  }
  // dPhi_x/dt = -[ (Ez_j+1 - Ez_j) - (Ey_k+1 - Ey_k) ]
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i <= geometry.nx; ++i) {
        const double curl_term =
            (emf_z[geometry.edgeIndexZ(i, j + 1U, k)] -
             emf_z[geometry.edgeIndexZ(i, j, k)]) -
            (emf_y[geometry.edgeIndexY(i, j, k + 1U)] -
             emf_y[geometry.edgeIndexY(i, j, k)]);
        phi_x[geometry.faceIndexX(i, j, k)] -= dt_code * curl_term;
      }
    }
  }
  // dPhi_y/dt = -[ (Ex_k+1 - Ex_k) - (Ez_i+1 - Ez_i) ]
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j <= geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const double curl_term =
            (emf_x[geometry.edgeIndexX(i, j, k + 1U)] -
             emf_x[geometry.edgeIndexX(i, j, k)]) -
            (emf_z[geometry.edgeIndexZ(i + 1U, j, k)] -
             emf_z[geometry.edgeIndexZ(i, j, k)]);
        phi_y[geometry.faceIndexY(i, j, k)] -= dt_code * curl_term;
      }
    }
  }
  // dPhi_z/dt = -[ (Ey_i+1 - Ey_i) - (Ex_j+1 - Ex_j) ]
  for (std::size_t k = 0; k <= geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const double curl_term =
            (emf_y[geometry.edgeIndexY(i + 1U, j, k)] -
             emf_y[geometry.edgeIndexY(i, j, k)]) -
            (emf_x[geometry.edgeIndexX(i, j + 1U, k)] -
             emf_x[geometry.edgeIndexX(i, j, k)]);
        phi_z[geometry.faceIndexZ(i, j, k)] -= dt_code * curl_term;
      }
    }
  }
}

double discreteComplexResidualX(const MhdPatchGeometry& geometry) {
  MhdFaceFluxStorage flux(geometry);
  MhdEdgeEmfStorage emf(geometry);
  auto ex = emf.emfX();
  for (std::size_t index = 0; index < ex.size(); ++index) {
    ex[index] = 1.0 + 0.5 * static_cast<double>(index % 7U);
  }
  applyConstrainedTransportUpdate(geometry, flux.fluxX(), flux.fluxY(), flux.fluxZ(),
                                  emf.emfX(), emf.emfY(), emf.emfZ(), 1.0);
  const std::vector<double> divergence =
      faceFluxDivergence(geometry, flux.fluxX(), flux.fluxY(), flux.fluxZ());
  double residual = 0.0;
  for (const double value : divergence) {
    residual = std::max(residual, std::abs(value));
  }
  return residual;
}

double discreteComplexResidualY(const MhdPatchGeometry& geometry) {
  MhdFaceFluxStorage flux(geometry);
  MhdEdgeEmfStorage emf(geometry);
  auto ey = emf.emfY();
  for (std::size_t index = 0; index < ey.size(); ++index) {
    ey[index] = -0.25 + 0.125 * static_cast<double>(index % 11U);
  }
  applyConstrainedTransportUpdate(geometry, flux.fluxX(), flux.fluxY(), flux.fluxZ(),
                                  emf.emfX(), emf.emfY(), emf.emfZ(), 1.0);
  const std::vector<double> divergence =
      faceFluxDivergence(geometry, flux.fluxX(), flux.fluxY(), flux.fluxZ());
  double residual = 0.0;
  for (const double value : divergence) {
    residual = std::max(residual, std::abs(value));
  }
  return residual;
}

double discreteComplexResidualZ(const MhdPatchGeometry& geometry) {
  MhdFaceFluxStorage flux(geometry);
  MhdEdgeEmfStorage emf(geometry);
  auto ez = emf.emfZ();
  for (std::size_t index = 0; index < ez.size(); ++index) {
    ez[index] = 2.0 - 0.375 * static_cast<double>(index % 13U);
  }
  applyConstrainedTransportUpdate(geometry, flux.fluxX(), flux.fluxY(), flux.fluxZ(),
                                  emf.emfX(), emf.emfY(), emf.emfZ(), 1.0);
  const std::vector<double> divergence =
      faceFluxDivergence(geometry, flux.fluxX(), flux.fluxY(), flux.fluxZ());
  double residual = 0.0;
  for (const double value : divergence) {
    residual = std::max(residual, std::abs(value));
  }
  return residual;
}

// ---------------------------------------------------------------------------
// Deterministic geometric identities
// ---------------------------------------------------------------------------

bool stableFaceKeyBefore(const MhdFaceKey& lhs, const MhdFaceKey& rhs) noexcept {
  if (lhs.refinement_level != rhs.refinement_level) {
    return lhs.refinement_level < rhs.refinement_level;
  }
  if (lhs.patch_id != rhs.patch_id) {
    return lhs.patch_id < rhs.patch_id;
  }
  if (lhs.axis != rhs.axis) {
    return static_cast<int>(lhs.axis) < static_cast<int>(rhs.axis);
  }
  if (lhs.i != rhs.i) {
    return lhs.i < rhs.i;
  }
  if (lhs.j != rhs.j) {
    return lhs.j < rhs.j;
  }
  return lhs.k < rhs.k;
}

bool stableEdgeKeyBefore(const MhdEdgeKey& lhs, const MhdEdgeKey& rhs) noexcept {
  if (lhs.refinement_level != rhs.refinement_level) {
    return lhs.refinement_level < rhs.refinement_level;
  }
  if (lhs.patch_id != rhs.patch_id) {
    return lhs.patch_id < rhs.patch_id;
  }
  if (lhs.axis != rhs.axis) {
    return static_cast<int>(lhs.axis) < static_cast<int>(rhs.axis);
  }
  if (lhs.i != rhs.i) {
    return lhs.i < rhs.i;
  }
  if (lhs.j != rhs.j) {
    return lhs.j < rhs.j;
  }
  return lhs.k < rhs.k;
}

std::uint64_t faceKeyHash(const MhdFaceKey& key) noexcept {
  std::uint64_t hash = k_fnv_offset_basis;
  hashMix(hash, key.patch_id);
  hashMix(hash, key.geometry_generation);
  hashMix(hash, signedToUnsigned(key.refinement_level));
  hashMix(hash, static_cast<std::uint64_t>(key.axis));
  hashMix(hash, signedToUnsigned(key.i));
  hashMix(hash, signedToUnsigned(key.j));
  hashMix(hash, signedToUnsigned(key.k));
  return hash;
}

std::uint64_t edgeKeyHash(const MhdEdgeKey& key) noexcept {
  std::uint64_t hash = k_fnv_offset_basis;
  hashMix(hash, key.patch_id);
  hashMix(hash, key.geometry_generation);
  hashMix(hash, signedToUnsigned(key.refinement_level));
  hashMix(hash, static_cast<std::uint64_t>(key.axis));
  // Edge keys are distinguished from face keys by an explicit domain tag.
  hashMix(hash, 0xE0E0E0E0ULL);
  hashMix(hash, signedToUnsigned(key.i));
  hashMix(hash, signedToUnsigned(key.j));
  hashMix(hash, signedToUnsigned(key.k));
  return hash;
}

std::int64_t globalFaceCoordinate(const MhdPatchDescriptor& patch,
                                  MhdAxis axis,
                                  std::size_t local_index_along_axis) noexcept {
  const std::int64_t local = static_cast<std::int64_t>(local_index_along_axis);
  switch (axis) {
    case MhdAxis::kX:
      return patch.origin_i + local;
    case MhdAxis::kY:
      return patch.origin_j + local;
    case MhdAxis::kZ:
      return patch.origin_k + local;
  }
  return local;
}

std::int64_t wrapPeriodic(std::int64_t index, std::int64_t extent) noexcept {
  if (extent <= 0) {
    return index;
  }
  std::int64_t wrapped = index % extent;
  if (wrapped < 0) {
    wrapped += extent;
  }
  return wrapped;
}

void validateBoundarySupport(const MhdPatchBoundary& boundary,
                             const MhdFailureContext& context) {
  const auto reject = [&](std::string detail) {
    MhdFailureContext failure = context;
    failure.detail = std::move(detail);
    throw MhdError(MhdErrorCode::kUnsupportedPhysics, std::move(failure));
  };

  for (int axis_index = 0; axis_index < 3; ++axis_index) {
    const std::size_t axis = static_cast<std::size_t>(axis_index);
    const char* name = axis_index == 0 ? "x" : (axis_index == 1 ? "y" : "z");
    const bool lower_periodic = boundary.lower[axis] == MhdBoundaryKind::kPeriodic;
    const bool upper_periodic = boundary.upper[axis] == MhdBoundaryKind::kPeriodic;
    if (lower_periodic != upper_periodic) {
      reject(std::string("one-sided periodic boundary on axis ") + name +
             " is not a supported MHD topology");
    }
    if (boundary.lower[axis] == MhdBoundaryKind::kImportedMpi ||
        boundary.upper[axis] == MhdBoundaryKind::kImportedMpi) {
      reject(std::string("remote-rank magnetic halo on axis ") + name +
             " requested, but no MPI face/edge halo protocol is implemented");
    }
  }
}

}  // namespace cosmosim::mhd
