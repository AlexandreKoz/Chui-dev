#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/mhd/mhd_state.hpp"

namespace cosmosim::mhd {

// Virtual (non-Cartesian tomography) edges are not part of REF-v1.
enum class MhdAxis : std::uint8_t { kX = 0, kY = 1, kZ = 2 };
enum class MhdFaceSide : std::uint8_t { kLower = 0, kUpper = 1 };

// Physical boundary classes supported by the reference. Shearing-periodic
// domains are explicitly unsupported (section 7.5) and have no enumerator.
enum class MhdBoundaryKind : std::uint8_t {
  kPeriodic = 0,
  kOpen = 1,
  kReflectingConductingWall = 2,
  // Remote-rank neighbour. No MPI magnetic halo protocol exists in this source
  // pass; any patch side of this kind is rejected by validateBoundarySupport.
  kImportedMpi = 3,
  // Same-rank, same-level neighbour patch. Ghost cells and padded staggered
  // faces on this side are filled by MhdLevelStepper from the neighbour's
  // stage state; shared faces/edges are reconciled from the deterministic
  // owner (lowest (level, patch_id)).
  kNeighborPatch = 4,
};

struct MhdPatchBoundary {
  std::array<MhdBoundaryKind, 3> lower{
      MhdBoundaryKind::kPeriodic, MhdBoundaryKind::kPeriodic, MhdBoundaryKind::kPeriodic};
  std::array<MhdBoundaryKind, 3> upper{
      MhdBoundaryKind::kPeriodic, MhdBoundaryKind::kPeriodic, MhdBoundaryKind::kPeriodic};
};

// Edge-integrated EMF storage. Values are oriented line integrals along the
// positive axis:
//   emf_x: (nx,   ny+1, nz+1)
//   emf_y: (nx+1, ny,   nz+1)
//   emf_z: (nx+1, ny+1, nz  )
// They are the single time-integrated authority synchronized by AMR/MPI.
class MhdEdgeEmfStorage final {
 public:
  MhdEdgeEmfStorage() = default;
  explicit MhdEdgeEmfStorage(const MhdPatchGeometry& geometry);

  void resize(const MhdPatchGeometry& geometry);
  void zero() noexcept;
  [[nodiscard]] const MhdPatchGeometry& geometry() const noexcept { return m_geometry; }
  [[nodiscard]] std::span<double> emfX() noexcept { return m_emf_x.span(); }
  [[nodiscard]] std::span<const double> emfX() const noexcept { return m_emf_x.span(); }
  [[nodiscard]] std::span<double> emfY() noexcept { return m_emf_y.span(); }
  [[nodiscard]] std::span<const double> emfY() const noexcept { return m_emf_y.span(); }
  [[nodiscard]] std::span<double> emfZ() noexcept { return m_emf_z.span(); }
  [[nodiscard]] std::span<const double> emfZ() const noexcept { return m_emf_z.span(); }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept {
    return m_emf_x.ownedCapacityBytes() + m_emf_y.ownedCapacityBytes() +
        m_emf_z.ownedCapacityBytes();
  }

 private:
  MhdPatchGeometry m_geometry{};
  MhdDoubleArray m_emf_x;
  MhdDoubleArray m_emf_y;
  MhdDoubleArray m_emf_z;
};

// Cell-centered derived field cache (ephemeral, never an authority).
struct MhdCellFieldCache {
  MhdPatchGeometry geometry{};
  MhdDoubleArray bx_comoving;
  MhdDoubleArray by_comoving;
  MhdDoubleArray bz_comoving;

  void refresh(const MhdPatchGeometry& patch_geometry,
               std::span<const double> phi_x,
               std::span<const double> phi_y,
               std::span<const double> phi_z);
  [[nodiscard]] std::array<double, 3> field(std::size_t i, std::size_t j, std::size_t k) const;
};

// Discrete face-to-cell divergence operator D. Returns the signed oriented
// flux imbalance per cell (no division by volume).
[[nodiscard]] std::vector<double> faceFluxDivergence(
    const MhdPatchGeometry& geometry,
    std::span<const double> phi_x,
    std::span<const double> phi_y,
    std::span<const double> phi_z);

// Maximum |D(Phi)| over cells divided by `normalizer` (see mhd_validation).
[[nodiscard]] double maxNormalizedFaceDivergence(
    const MhdPatchGeometry& geometry,
    std::span<const double> phi_x,
    std::span<const double> phi_y,
    std::span<const double> phi_z,
    double normalizer);

// Sum of |Phi| over all faces, used as one divergence normalizer candidate.
[[nodiscard]] double faceFluxMagnitudeSum(std::span<const double> phi_x,
                                          std::span<const double> phi_y,
                                          std::span<const double> phi_z);

// Apply the discrete Faraday update Phi <- Phi - dt * C(E) on one patch.
// `emf_*` hold time-integrated edge values (already multiplied by the interval
// they represent); dt_code must be 1.0 when the integrals are pre-integrated.
void applyConstrainedTransportUpdate(const MhdPatchGeometry& geometry,
                                     std::span<double> phi_x,
                                     std::span<double> phi_y,
                                     std::span<double> phi_z,
                                     std::span<const double> emf_x,
                                     std::span<const double> emf_y,
                                     std::span<const double> emf_z,
                                     double dt_code);

// Discrete on-patch tests used by unit tests and preflight: verifies
// D*C = 0 for every interior cell on the periodic patch without requiring the
// full solver. Returns the maximum residual of the composed operator.
[[nodiscard]] double discreteComplexResidualX(const MhdPatchGeometry& geometry);
[[nodiscard]] double discreteComplexResidualY(const MhdPatchGeometry& geometry);
[[nodiscard]] double discreteComplexResidualZ(const MhdPatchGeometry& geometry);

// Face and edge identities. Identity is geometric, never a dense gas row.
struct MhdFaceKey {
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::int32_t refinement_level = 0;
  MhdAxis axis = MhdAxis::kX;
  std::int64_t i = 0;
  std::int64_t j = 0;
  std::int64_t k = 0;
};

struct MhdEdgeKey {
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::int32_t refinement_level = 0;
  MhdAxis axis = MhdAxis::kX;
  std::int64_t i = 0;
  std::int64_t j = 0;
  std::int64_t k = 0;
};

struct MhdPatchDescriptor {
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::int32_t refinement_level = 0;
  MhdPatchGeometry geometry{};
  // Lower corner in the level's integer lattice.
  std::int64_t origin_i = 0;
  std::int64_t origin_j = 0;
  std::int64_t origin_k = 0;
};

// Deterministic stable ordering: lowest (refinement_level, patch_id, axis,
// i, j, k) tuple owns a shared face/edge. Returns true when lhs precedes rhs.
[[nodiscard]] bool stableFaceKeyBefore(const MhdFaceKey& lhs, const MhdFaceKey& rhs) noexcept;
[[nodiscard]] bool stableEdgeKeyBefore(const MhdEdgeKey& lhs, const MhdEdgeKey& rhs) noexcept;

// Compact deterministic hash of a geometric key (64-bit FNV-1a over the
// canonical field order). Stable across ranks and builds.
[[nodiscard]] std::uint64_t faceKeyHash(const MhdFaceKey& key) noexcept;
[[nodiscard]] std::uint64_t edgeKeyHash(const MhdEdgeKey& key) noexcept;

// Global lattice coordinate of a local face index.
[[nodiscard]] std::int64_t globalFaceCoordinate(const MhdPatchDescriptor& patch,
                                                MhdAxis axis,
                                                std::size_t local_index_along_axis) noexcept;

// Periodic wrap helper for one axis index within [0, extent).
[[nodiscard]] std::int64_t wrapPeriodic(std::int64_t index, std::int64_t extent) noexcept;

// Refuse unsupported boundary combinations before allocation.
void validateBoundarySupport(const MhdPatchBoundary& boundary, const MhdFailureContext& context);

}  // namespace cosmosim::mhd
