#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <utility>

#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_induction.hpp"
#include "cosmosim/mhd/mhd_riemann.hpp"
#include "cosmosim/mhd/mhd_state.hpp"
#include "cosmosim/mhd/mhd_transport.hpp"
#include "cosmosim/mhd/mhd_validation.hpp"

namespace cosmosim::mhd {

enum class MhdReconstructionMode : std::uint8_t { kPlmMc = 0, kPpm = 1 };
enum class MhdTimeIntegrator : std::uint8_t { kSsprk2 = 0 };
enum class MhdDivergencePolicy : std::uint8_t { kStrictCt = 0 };
enum class MhdEmfScheme : std::uint8_t { kUctHll = 0 };
enum class MhdSpatialMethod : std::uint8_t { kFiniteVolumeUnsplit = 0 };
enum class MhdMagneticRepresentation : std::uint8_t { kFaceFluxCt = 0 };

[[nodiscard]] std::string_view mhdReconstructionName(MhdReconstructionMode mode) noexcept;

// Ghost depth required by each reconstruction (CHUI-MHD-REF v1.0 section 6.2):
// PLM-MC needs cells i-1..i+1 for the face states of cell i and the face sweep
// reconstructs cells -1..n, hence 2 ghost layers; bounded PPM needs i-2..i+2,
// hence 3. The UCT-HLL builder needs 2 (MC corner slopes at cell -1).
[[nodiscard]] constexpr std::size_t mhdGhostCells(MhdReconstructionMode mode) noexcept {
  return mode == MhdReconstructionMode::kPpm ? 3U : 2U;
}
inline constexpr std::size_t k_mhd_staggered_tangential_ghost = 2U;

struct MhdNumericsPolicy {
  MhdReconstructionMode reconstruction = MhdReconstructionMode::kPpm;
  MhdRiemann::Solver riemann = MhdRiemann::Solver::kHlld;
  MhdEmfScheme emf_scheme = MhdEmfScheme::kUctHll;
  MhdTimeIntegrator time_integrator = MhdTimeIntegrator::kSsprk2;
  MhdDivergencePolicy divergence_policy = MhdDivergencePolicy::kStrictCt;
  // Calorically perfect gas. REF-v1 certification uses 5/3; the runtime
  // adapter rejects any EOS whose MHD-compatible closure is not established.
  double gamma_adiabatic = k_mhd_gamma_default;
  double cfl_ideal = 0.30;
  std::int32_t max_step_retries = 8;
  double fallback_fraction_limit = 0.001;
  double floor_energy_fraction_limit = 1.0e-8;
  // Post-step strict-CT divergence acceptance (G03 normalization).
  double divergence_tolerance = 1.0e-11;
  bool deterministic = true;
  MhdStiffnessLimits stiffness{};

  [[nodiscard]] std::size_t ghostCells() const noexcept { return mhdGhostCells(reconstruction); }
};

// Gas patch state (comoving conserved). In MHD mode the energy lane is the
// single authoritative total energy E_c = rho_c e + K_c + B_c^2/(2a); the
// cell-centred B_c entering M_c is always re-derived from the authoritative
// face fluxes, never stored here.
class MhdGasStateSoa final {
 public:
  static constexpr std::size_t k_lane_count = 6U;

  MhdGasStateSoa() = default;
  explicit MhdGasStateSoa(std::size_t cell_count);
  void resize(std::size_t cell_count);
  [[nodiscard]] std::size_t size() const noexcept { return m_rho.size(); }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;

  [[nodiscard]] std::span<double> rho() noexcept { return m_rho.span(); }
  [[nodiscard]] std::span<const double> rho() const noexcept { return m_rho.span(); }
  [[nodiscard]] std::span<double> momentumX() noexcept { return m_mx.span(); }
  [[nodiscard]] std::span<const double> momentumX() const noexcept { return m_mx.span(); }
  [[nodiscard]] std::span<double> momentumY() noexcept { return m_my.span(); }
  [[nodiscard]] std::span<const double> momentumY() const noexcept { return m_my.span(); }
  [[nodiscard]] std::span<double> momentumZ() noexcept { return m_mz.span(); }
  [[nodiscard]] std::span<const double> momentumZ() const noexcept { return m_mz.span(); }
  [[nodiscard]] std::span<double> totalEnergy() noexcept { return m_energy.span(); }
  [[nodiscard]] std::span<const double> totalEnergy() const noexcept { return m_energy.span(); }
  [[nodiscard]] std::span<double> metalMass() noexcept { return m_metal.span(); }
  [[nodiscard]] std::span<const double> metalMass() const noexcept { return m_metal.span(); }
  [[nodiscard]] std::span<double> lane(std::size_t index) noexcept;
  [[nodiscard]] std::span<const double> lane(std::size_t index) const noexcept;

  // Gas lanes only (rho, m, E, metal density); the B lanes of the returned
  // MhdConservedState are zero and must be filled from face fluxes.
  [[nodiscard]] MhdConservedState loadGas(std::size_t cell) const noexcept;
  void storeGas(std::size_t cell, const MhdConservedState& state) noexcept;
  void copyFrom(const MhdGasStateSoa& other);
  // this = alpha * lhs + beta * rhs (all lanes, cell-wise).
  void assignLinearCombination(double alpha, const MhdGasStateSoa& lhs,
                               double beta, const MhdGasStateSoa& rhs);

 private:
  MhdDoubleArray m_rho;
  MhdDoubleArray m_mx;
  MhdDoubleArray m_my;
  MhdDoubleArray m_mz;
  MhdDoubleArray m_energy;
  MhdDoubleArray m_metal;
};

// Numerical gas fluxes on the six patch-boundary face planes. A boundary face
// flux is not applied by the patch sweep: the level stepper reconciles each
// shared logical face from its deterministic owner (identical flux on both
// sides, exact conservation) and then applies it. Storage is O(patch surface).
struct MhdBoundaryFaceFluxes {
  static constexpr std::size_t k_values_per_face = MhdGasStateSoa::k_lane_count;
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  // [axis][side]: side 0 = lower plane (face index 0), 1 = upper (index n).
  std::array<std::array<MhdDoubleArray, 2>, 3> values{};

  void resize(const MhdPatchGeometry& geometry);
  [[nodiscard]] std::size_t faceCount(std::size_t axis) const noexcept;
  // Tangential coordinates (t1, t2) are the cyclic successors of `axis`.
  [[nodiscard]] double* at(std::size_t axis, std::size_t side, std::size_t t1,
                           std::size_t t2) noexcept;
  [[nodiscard]] const double* at(std::size_t axis, std::size_t side, std::size_t t1,
                                 std::size_t t2) const noexcept;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

struct MhdStageCounters {
  std::uint64_t faces_solved = 0U;
  std::uint64_t hlle_fallback_faces = 0U;
  std::uint64_t first_order_faces = 0U;
  std::uint64_t ppm_characteristic_fallback_cells = 0U;
  std::uint64_t edges_built = 0U;

  MhdStageCounters& operator+=(const MhdStageCounters& rhs) noexcept;
};

// Patch-local phase workspace. One instance per active patch; reused across
// SSPRK2 stages, RKL stages and Hall subcycles (stage count never multiplies
// residency). All lanes are MhdDoubleArray with retained capacity accounted
// by ownedCapacityBytes().
struct MhdStageWorkspace {
  MhdPrimitivePatch primitives{};
  MhdStaggeredFieldPatch staggered{};
  MhdFaceSignalPatch signals{};
  MhdEdgeEmfStorage emf{};
  MhdBoundaryFaceFluxes boundary_fluxes{};

  void resize(const MhdPatchGeometry& geometry, MhdReconstructionMode reconstruction);
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

// Bounded face states of one cell along one axis (characteristic PPM or
// primitive PLM-MC). `lower`/`upper` are the cell's values at its lower and
// upper faces in the global frame; magnetic lanes are B_eff.
struct MhdCellFaceStates {
  MhdPrimitiveState lower{};
  MhdPrimitiveState upper{};
};

// Ideal-MHD reference patch kernels (CHUI-MHD-REF v1.0 section 6). There is no
// per-face virtual dispatch or allocation; algorithm selection happens before
// the face loops. These kernels never commit a step on their own: the
// interval-level owner (MhdLevelStepper) holds the rollback image, the retry
// policy and the shared-face/edge reconciliation.
class MhdStepOperator {
 public:
  explicit MhdStepOperator(MhdNumericsPolicy policy = {});

  [[nodiscard]] const MhdNumericsPolicy& policy() const noexcept { return m_policy; }

  // Recover interior primitives (B lanes = B_eff) and the interior staggered
  // normal field B_eff from the gas state and authoritative face fluxes at the
  // stage scale factor. Returns the number of inadmissible cells (rho<=0,
  // thermal<=0, Z outside [0,1], or non-finite); no floor is applied.
  [[nodiscard]] std::size_t recoverPrimitives(const MhdPatchGeometry& geometry,
                                              const MhdGasStateSoa& gas,
                                              const MhdFaceFluxView& faces,
                                              double scale_factor,
                                              MhdStageWorkspace& workspace) const;

  // Single-cell primitive recovery with exactly the arithmetic of
  // recoverPrimitives (used by the level stepper to sample neighbour, finer
  // and coarser leaves for ghost data). Returns false for an inadmissible
  // state; never floors.
  [[nodiscard]] bool recoverCellPrimitive(const MhdConservedState& gas,
                                          const std::array<double, 3>& b_comoving,
                                          double scale_factor,
                                          MhdPrimitiveState& out) const noexcept;
  // Cell-centred comoving field B_c from authoritative face fluxes
  // (arithmetic mean of the two bordering face fields per axis).
  [[nodiscard]] static std::array<double, 3> cellFieldFromFaces(
      const MhdPatchGeometry& geometry, std::span<const double> phi_x,
      std::span<const double> phi_y, std::span<const double> phi_z, std::size_t i,
      std::size_t j, std::size_t k) noexcept;

  // Fill ghost cells and staggered ghost rows for every patch side whose kind
  // is physical (periodic self-wrap, open, conducting wall). kNeighborPatch
  // sides must already have been filled by the level stepper.
  static void fillPhysicalGhosts(const MhdPatchBoundary& boundary, MhdStageWorkspace& workspace);

  // One SSPRK2 operator evaluation L(Y) at stage scale factor `a_stage` and
  // Hubble rate `hubble_stage`, applied with `flux_weight` (= dt for stage 1,
  // dt/2 for stage 2) into `gas_target`:
  //   gas_target[cell] += flux_weight * ( -div(F)/a + S_expansion )
  // for all INTERIOR faces; patch-boundary faces are recorded in
  // workspace.boundary_fluxes for owner reconciliation. The UCT-HLL edge EMF
  // rate is left in workspace.emf (not yet applied). Workspace primitives and
  // staggered fields must hold the stage state with ghosts filled.
  // Throws MhdError(kNonPositiveState) when no admissible flux exists at a
  // face (HLLD -> HLLE -> first-order cell-average Rusanov all failed).
  void evaluateIdealStage(const MhdPatchGeometry& geometry,
                          const MhdPatchBoundary& boundary,
                          double a_stage,
                          double hubble_stage,
                          bool apply_expansion_source,
                          double flux_weight,
                          MhdStageWorkspace& workspace,
                          MhdGasStateSoa& gas_target,
                          MhdStageCounters& counters) const;

  // Apply the (reconciled) boundary-face gas fluxes:
  //   lower plane f=0  -> cell 0       gains +w F/(a dx)
  //   upper plane f=n  -> cell n-1     loses -w F/(a dx)
  static void applyBoundaryFaceFluxes(const MhdPatchGeometry& geometry,
                                      const MhdBoundaryFaceFluxes& fluxes,
                                      double a_stage,
                                      double flux_weight,
                                      MhdGasStateSoa& gas_target);

  // Fast-magnetosonic CFL candidate: dt <= cfl / max_cell sum_axis
  // (|u_axis| + c_f,axis) / (a dx_axis) with physical signal speeds from
  // (rho_c, p_c, B_eff). Throws MHD_NONPOSITIVE_STATE for an inadmissible cell
  // (never silently skipped).
  [[nodiscard]] double estimateIdealTimeStep(const MhdPatchGeometry& geometry,
                                             const MhdGasStateSoa& gas,
                                             const MhdFaceFluxView& faces,
                                             double scale_factor,
                                             double cfl) const;

  // Conservative budget of a patch (comoving volume integrals; magnetic energy
  // uses B_c^2/(2a) at the supplied scale factor).
  [[nodiscard]] MhdConservationBudget patchBudget(const MhdPatchGeometry& geometry,
                                                  const MhdGasStateSoa& gas,
                                                  const MhdFaceFluxView& faces,
                                                  double scale_factor) const;

  // Characteristic PPM (bounded, Colella-Woodward limiting in the MHD
  // primitive eigenbasis of the cell) or primitive PLM-MC face states of cell
  // `cell` along `axis`. Returns false only for an inadmissible cell; a
  // characteristic projection that is singular or yields a non-positive face
  // density/pressure falls back to PLM-MC for that cell (counted).
  [[nodiscard]] bool reconstructCell(const MhdPrimitivePatch& primitives,
                                     const MhdIndex3& cell,
                                     MhdAxis axis,
                                     MhdCellFaceStates& states,
                                     std::uint64_t* characteristic_fallbacks) const;

  // Scalar one-dimensional limiters, exposed for unit tests. Both return
  // {value at lower face, value at upper face} of the centre cell. PLM-MC uses
  // only (c_minus, c, c_plus); scalar PPM uses all five values.
  [[nodiscard]] static std::pair<double, double> reconstructPlmMc(
      double c_minus, double c, double c_plus) noexcept;
  [[nodiscard]] static std::pair<double, double> reconstructPpmScalar(
      double c_minus2, double c_minus, double c, double c_plus, double c_plus2) noexcept;

  // Apply physical ghost rules to the padded primitive cache (exposed for
  // tests): periodic wrap; open zero-gradient; conducting wall with odd
  // wall-normal velocity, even wall-normal field and odd wall-tangential field.
  static void fillGhostPrimitives(const MhdPatchBoundary& boundary, MhdPrimitivePatch& primitives);

 private:
  MhdNumericsPolicy m_policy;
};

}  // namespace cosmosim::mhd
