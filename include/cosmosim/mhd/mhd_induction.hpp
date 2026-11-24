#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_riemann.hpp"

namespace cosmosim::mhd {

// Per-face upwind signals produced by the Riemann sweep, stored axis-major for
// contiguous traversal. `mass_flux` is owner->neighbour along the face normal.
// `electric_t1`/`electric_t2` are the two tangential ideal electric-field
// components in the local right-handed basis (n, t1, t2) described in
// mhd_riemann.hpp; they are NOT advected scalars.
struct MhdFaceSignalPatch {
  MhdPatchGeometry geometry{};
  std::vector<double> mass_flux_x;
  std::vector<double> electric_t1_x;
  std::vector<double> electric_t2_x;
  std::vector<double> face_speed_x;

  std::vector<double> mass_flux_y;
  std::vector<double> electric_t1_y;
  std::vector<double> electric_t2_y;
  std::vector<double> face_speed_y;

  std::vector<double> mass_flux_z;
  std::vector<double> electric_t1_z;
  std::vector<double> electric_t2_z;
  std::vector<double> face_speed_z;

  void resize(const MhdPatchGeometry& patch_geometry);
};

// Upwind UCT edge EMF construction (CHUI-MHD-REF v1.0 section 6.1). The edge
// value is the inverse-fast-wave-speed weighted combination of the four
// adjacent face electric-field estimates. This is an HLL-qualified upwind
// average, not an unqualified arithmetic four-face average: a face whose
// signal fan is faster receives a proportionally smaller weight, and a
// degenerate (all-zero-speed) configuration uses the arithmetic mean only
// after proving the four estimates agree within tolerance.
class MhdEmfBuilder {
 public:
  static constexpr double k_degenerate_weight_floor = 1.0e-30;
  static constexpr double k_disagreement_relative_tolerance = 1.0e-10;

  // Build all three edge arrays from the face signals. `emf` is overwritten.
  static void buildUctHll(const MhdPatchGeometry& geometry,
                          const MhdFaceSignalPatch& signals,
                          MhdEdgeEmfStorage& emf,
                          std::uint64_t* degenerate_edge_count = nullptr);

  // Compose the ideal electric field at a face from a Riemann flux in the
  // local basis. Exposed so the solver and tests share one sign convention.
  [[nodiscard]] static std::array<double, 3> faceElectricFieldFromFlux(
      const MhdRiemannFlux& flux) noexcept;
};

// Inductive (curl) update of the magnetic face fluxes for one stage interval.
// `emf` holds time-integrated edge values for this interval.
void applyInductionUpdate(const MhdPatchGeometry& geometry,
                          std::span<double> phi_x,
                          std::span<double> phi_y,
                          std::span<double> phi_z,
                          const MhdEdgeEmfStorage& emf);

// Accumulate an interval into a running edge-integral authority. This is the
// operation AMR/MPI synchronizes; no per-substep edge history is retained.
void accumulateEdgeEmfIntegral(const MhdPatchGeometry& geometry,
                               MhdEdgeEmfStorage& accumulator,
                               const MhdEdgeEmfStorage& interval_emf);

// Clamp-free magnetic-magnitude normalizer. Returns max(eps_scale, max|Phi|),
// where eps_scale is the specification's 1e-30 * B_unit_code * A_ref surrogate
// supplied by the caller. It never fabricates a finite field from B=0.
[[nodiscard]] double magneticDivergenceNormalizer(
    std::span<const double> phi_x,
    std::span<const double> phi_y,
    std::span<const double> phi_z,
    double eps_scale);

}  // namespace cosmosim::mhd
