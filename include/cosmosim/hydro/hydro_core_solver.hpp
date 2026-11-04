#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace cosmosim::hydro {

constexpr std::size_t k_invalid_cell_index = static_cast<std::size_t>(-1);
constexpr std::size_t k_invalid_ghost_cell_slot = static_cast<std::size_t>(-1);
constexpr std::size_t k_hydro_automatic_active_batch_max_cells = 16384U;
constexpr std::size_t k_hydro_active_batch_alignment_cells = 8U;
constexpr std::uint64_t k_hydro_runtime_batch_scratch_budget_bytes_per_cell = 4096ULL;

// Cell-centered primitive variables in comoving coordinates.
// rho_comoving: comoving mass density.
// vel_*_peculiar: peculiar velocity components.
// pressure_comoving: comoving-frame thermal pressure.
struct HydroPrimitiveState {
  double rho_comoving = 0.0;
  double vel_x_peculiar = 0.0;
  double vel_y_peculiar = 0.0;
  double vel_z_peculiar = 0.0;
  double pressure_comoving = 0.0;
  // Specific internal energy carried by the conserved total-energy lane. This
  // is explicit because a barotropic/effective closure need not satisfy
  // P=(gamma-1)rho*u.
  double specific_internal_energy_code = 0.0;
  // Characteristic signal speed squared used by reconstruction/Riemann/CFL.
  // For ideal gas this is gamma*P/rho; effective closures provide dP/drho.
  double signal_speed_squared_code = 0.0;
  bool uses_effective_ism = false;
  // Passive gas metal mass fraction, reconstructed with the hydrodynamic state.
  double metallicity_mass_fraction = 0.0;
};

// Cell-centered conserved variables in comoving coordinates.
// mass_density_comoving: rho.
// momentum_density_*_comoving: rho * u.
// total_energy_density_comoving: rho * (e_internal + 0.5 |u|^2).
struct HydroConservedState {
  double mass_density_comoving = 0.0;
  double momentum_density_x_comoving = 0.0;
  double momentum_density_y_comoving = 0.0;
  double momentum_density_z_comoving = 0.0;
  double total_energy_density_comoving = 0.0;
  // Conserved passive scalar rho * Z.
  double metal_mass_density_comoving = 0.0;

  HydroConservedState& operator+=(const HydroConservedState& rhs);
  HydroConservedState& operator-=(const HydroConservedState& rhs);
};

[[nodiscard]] HydroConservedState operator+(HydroConservedState lhs, const HydroConservedState& rhs);
[[nodiscard]] HydroConservedState operator-(HydroConservedState lhs, const HydroConservedState& rhs);
[[nodiscard]] HydroConservedState operator*(double scalar, HydroConservedState state);

class HydroConservedStateSoa {
 public:
  explicit HydroConservedStateSoa(std::size_t cell_count = 0);

  void resize(std::size_t cell_count);
  [[nodiscard]] std::size_t size() const;
  [[nodiscard]] std::uint64_t currentSizeBytes() const;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const;

  [[nodiscard]] HydroConservedState loadCell(std::size_t cell_index) const;
  void storeCell(std::size_t cell_index, const HydroConservedState& cell_state);

  [[nodiscard]] std::span<double> massDensityComoving();
  [[nodiscard]] std::span<const double> massDensityComoving() const;

  [[nodiscard]] std::span<double> momentumDensityXComoving();
  [[nodiscard]] std::span<const double> momentumDensityXComoving() const;

  [[nodiscard]] std::span<double> momentumDensityYComoving();
  [[nodiscard]] std::span<const double> momentumDensityYComoving() const;

  [[nodiscard]] std::span<double> momentumDensityZComoving();
  [[nodiscard]] std::span<const double> momentumDensityZComoving() const;

  [[nodiscard]] std::span<double> totalEnergyDensityComoving();
  [[nodiscard]] std::span<const double> totalEnergyDensityComoving() const;

  [[nodiscard]] std::span<double> metalMassDensityComoving();
  [[nodiscard]] std::span<const double> metalMassDensityComoving() const;

  // Optional dual-energy lane: auxiliary internal-energy density rho*e
  // (comoving). It is NOT a conserved flux lane: HydroCoreSolver evolves it
  // with the non-conservative internal-energy equation and writes the
  // selected thermal energy back into it. Callers initialize it from the
  // canonical specific internal energy (GasCellSidecar::internal_energy_code)
  // before an advance; NaN marks "unknown" (the solver then falls back to
  // E - K for that row and counts the fallback). Enabling/resizing fills new
  // rows with NaN; disabling releases the lane.
  void enableAuxiliaryInternalEnergy(bool enabled);
  [[nodiscard]] bool auxiliaryInternalEnergyEnabled() const noexcept;
  [[nodiscard]] std::span<double> auxiliaryInternalEnergyDensityComoving();
  [[nodiscard]] std::span<const double> auxiliaryInternalEnergyDensityComoving() const;

 private:
  std::vector<double> m_mass_density_comoving;
  std::vector<double> m_momentum_density_x_comoving;
  std::vector<double> m_momentum_density_y_comoving;
  std::vector<double> m_momentum_density_z_comoving;
  std::vector<double> m_total_energy_density_comoving;
  std::vector<double> m_metal_mass_density_comoving;
  std::vector<double> m_auxiliary_internal_energy_density_comoving;
  bool m_auxiliary_internal_energy_enabled = false;
};

class HydroPrimitiveCacheSoa {
 public:
  explicit HydroPrimitiveCacheSoa(std::size_t cell_count = 0);

  void resize(std::size_t cell_count);
  [[nodiscard]] std::size_t size() const;

  [[nodiscard]] HydroPrimitiveState loadCell(std::size_t cell_index) const;
  void storeCell(std::size_t cell_index, const HydroPrimitiveState& primitive_state);

 private:
  std::vector<double> m_rho_comoving;
  std::vector<double> m_vel_x_peculiar;
  std::vector<double> m_vel_y_peculiar;
  std::vector<double> m_vel_z_peculiar;
  std::vector<double> m_pressure_comoving;
  std::vector<double> m_specific_internal_energy_code;
  std::vector<double> m_signal_speed_squared_code;
  std::vector<std::uint8_t> m_uses_effective_ism;
  std::vector<double> m_metallicity_mass_fraction;
};

struct HydroPatchColdData {
  std::uint64_t patch_id = 0;
  int refinement_level = 0;
};

enum class HydroFaceAxis {
  kX,
  kY,
  kZ,
};

enum class HydroFaceSide {
  kLower,
  kUpper,
};

enum class HydroFluxRegisterFaceRole {
  kNone,
  kCoarse,
  kFine,
};

struct HydroFluxRegisterFace {
  HydroFluxRegisterFaceRole role = HydroFluxRegisterFaceRole::kNone;
  std::uint64_t register_key = 0;
  std::uint64_t coarse_patch_id = 0;
  std::uint64_t coarse_gas_cell_id = 0;
  std::size_t coarse_cell_index = k_invalid_cell_index;
  int level = 0;
  HydroFaceAxis axis = HydroFaceAxis::kX;
  HydroFaceSide orientation = HydroFaceSide::kLower;
  double coarse_orientation_sign = 1.0;
};

enum class HydroBoundaryKind {
  kInterior,
  kPeriodic,
  kOpen,
  kReflective,
  kImportedMpi,
};

enum class HydroGhostMutationRights {
  kReadOnlyImported,
  kWritablePhysicalBoundaryScratch,
};

struct HydroGhostCell {
  std::size_t owner_real_cell = k_invalid_cell_index;
  std::size_t source_real_cell = k_invalid_cell_index;
  std::size_t ghost_cell = k_invalid_cell_index;
  std::size_t ghost_slot = k_invalid_ghost_cell_slot;
  HydroBoundaryKind boundary_kind = HydroBoundaryKind::kInterior;
  HydroFaceAxis axis = HydroFaceAxis::kX;
  HydroFaceSide side = HydroFaceSide::kLower;
  HydroGhostMutationRights mutation_rights = HydroGhostMutationRights::kWritablePhysicalBoundaryScratch;
  std::uint64_t origin_gas_cell_id = 0;
  int origin_rank = -1;
  std::uint64_t hydro_sync_epoch = 0;
  std::uint64_t decomposition_epoch = 0;
};

// One Cartesian face. Reconstruction stencil cells are ordered ALONG THE
// FACE NORMAL n (owner -> neighbor is +n):
//   owner_minus2_cell -> owner_minus_cell -> owner_cell | face |
//   neighbor_cell -> neighbor_plus_cell -> neighbor_plus2_cell
// k_invalid_cell_index marks an unavailable stencil cell. Interior faces of
// makeCartesianPatchGeometry have n = +axis; physical-boundary ghost faces on
// the lower side have n = -axis and no stencil beyond owner/neighbor.
// owner_minus2_cell / neighbor_plus2_cell are consumed only by PPM; they are
// appended after `axis` so existing designated initializers stay valid.
struct HydroFace {
  std::size_t owner_cell = k_invalid_cell_index;
  std::size_t neighbor_cell = k_invalid_cell_index;
  std::size_t owner_minus_cell = k_invalid_cell_index;
  std::size_t neighbor_plus_cell = k_invalid_cell_index;
  std::size_t ghost_cell_slot = k_invalid_ghost_cell_slot;
  double area_comoving = 0.0;
  double normal_x = 0.0;
  double normal_y = 0.0;
  double normal_z = 0.0;
  HydroFaceAxis axis = HydroFaceAxis::kX;
  std::size_t owner_minus2_cell = k_invalid_cell_index;
  std::size_t neighbor_plus2_cell = k_invalid_cell_index;
};

struct HydroPatchGeometry {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  double origin_x_comoving = 0.0;
  double origin_y_comoving = 0.0;
  double origin_z_comoving = 0.0;
  double cell_width_x_comoving = 0.0;
  double cell_width_y_comoving = 0.0;
  double cell_width_z_comoving = 0.0;
  double cell_volume_comoving = 0.0;
  std::vector<HydroFace> faces;
  std::vector<HydroGhostCell> ghost_cells;
  std::vector<HydroFluxRegisterFace> flux_register_faces;

  [[nodiscard]] std::size_t cellCount() const;
  [[nodiscard]] std::size_t totalCellStorageCount() const;
  [[nodiscard]] std::size_t linearCellIndex(std::size_t i, std::size_t j, std::size_t k) const;
  [[nodiscard]] std::array<std::size_t, 3> cellIjk(std::size_t row) const;
  [[nodiscard]] std::size_t neighborCell(
      std::size_t row,
      int di,
      int dj,
      int dk) const;
};

struct HydroFluxRegisterRecord {
  HydroFluxRegisterFaceRole role = HydroFluxRegisterFaceRole::kNone;
  std::uint64_t register_key = 0;
  std::uint64_t coarse_patch_id = 0;
  std::uint64_t coarse_gas_cell_id = 0;
  std::size_t coarse_cell_index = k_invalid_cell_index;
  int level = 0;
  HydroFaceAxis axis = HydroFaceAxis::kX;
  HydroFaceSide orientation = HydroFaceSide::kLower;
  double face_area_comoving = 0.0;
  double dt_code = 0.0;
  HydroConservedState flux_code;
};

class HydroFluxRegisterSink {
 public:
  virtual ~HydroFluxRegisterSink() = default;

  virtual void recordFaceFlux(const HydroFluxRegisterRecord& record) = 0;
};

struct HydroActiveSetView {
  std::span<const std::size_t> active_cells;
  std::span<const std::size_t> active_faces;
};

struct HydroUpdateContext {
  double dt_code = 0.0;
  // Explicit source/transport interval. The scalar scale_factor/hubble_rate
  // are the evaluation epoch (normally the interval midpoint).
  double time_begin_code = 0.0;
  double time_end_code = 0.0;
  double scale_factor_begin = 1.0;
  double scale_factor = 1.0;
  double scale_factor_end = 1.0;
  // H(a) in inverse code time, owned by the cosmological timeline.
  double hubble_rate_begin_code = 0.0;
  double hubble_rate_code = 0.0;
  double hubble_rate_end_code = 0.0;
  bool comoving_coordinates = false;
};

struct HydroThermodynamicClosureResult {
  double pressure_comoving = 0.0;
  double signal_speed_squared_code = 0.0;
  double target_specific_internal_energy_code = 0.0;
  bool uses_effective_ism = false;
  bool valid = true;
};

class HydroThermodynamicClosure {
 public:
  virtual ~HydroThermodynamicClosure() = default;
  [[nodiscard]] virtual HydroThermodynamicClosureResult evaluate(
      std::size_t cell_index,
      const HydroConservedState& conserved,
      const HydroPrimitiveState& ideal_primitive,
      double scale_factor,
      double redshift) const = 0;
};

struct HydroCellSourceProperties {
  double mass_density_physical_cgs = 0.0;
  double hydrogen_number_density_cgs = 0.0;
  double metallicity_mass_fraction = 0.0;
  double temperature_k = 0.0;
};

class HydroCellSourcePropertyProvider {
 public:
  virtual ~HydroCellSourcePropertyProvider() = default;
  [[nodiscard]] virtual HydroCellSourceProperties sourcePropertiesForCell(
      std::size_t cell_index) const = 0;
};

// CHUI Part 2: per-face mass transfer of the FINAL (conservative) flux of
// one advance, used to transport composition lanes with the identical mass
// flux (consistent multi-fluid advection: species flux = F_rho * X_donor).
// mass_transfer_code = (dt/a) * A_f * F_rho, i.e. code MASS moved from
// owner_cell to neighbor_cell across the face (negative = reverse). Cells are
// solver geometry rows; neighbor_cell is k_invalid_cell_index for outflow /
// physical-boundary faces without a real or ghost target. Called once per
// active face per advance, single-threaded, in active-face order.
class HydroFaceMassFluxSink {
 public:
  virtual ~HydroFaceMassFluxSink() = default;
  // Called by the RUNTIME (fixed grid / AMR orchestrator) around each patch
  // advance so the sink can resolve solver rows into stable identities.
  // dense_row_by_real_cell maps real geometry rows [0, cellCount) to
  // SimulationState rows; ghost_origin_gas_cell_id_by_slot gives the stable
  // id of the cell that supplied each ghost slot (0 = unknown/physical).
  virtual void beginAdvance(const HydroPatchGeometry& geometry,
                            std::span<const std::uint32_t> dense_row_by_real_cell,
                            std::span<const std::uint64_t> ghost_origin_gas_cell_id_by_slot) {
    (void)geometry;
    (void)dense_row_by_real_cell;
    (void)ghost_origin_gas_cell_id_by_slot;
  }
  virtual void endAdvance() noexcept {}
  virtual void recordFaceMassTransfer(std::size_t face_index, std::size_t owner_cell,
                                      std::size_t neighbor_cell, double mass_transfer_code) = 0;
};

struct HydroSourceContext {
  HydroUpdateContext update;
  // Scale-free comoving gravitational kernel A from TreePM. The comoving
  // gravity source applies the A/a^2 peculiar-velocity response exactly once.
  std::span<const double> gravity_accel_x_peculiar;
  std::span<const double> gravity_accel_y_peculiar;
  std::span<const double> gravity_accel_z_peculiar;
  // Physical mass density in g cm^-3 at the source-evaluation epoch.
  // Cooling must not infer this quantity from a comoving code-density lane.
  std::span<const double> mass_density_physical_cgs;
  std::span<const double> hydrogen_number_density_cgs;
  std::span<const double> metallicity_mass_fraction;
  std::span<const double> temperature_k;
  // Optional on-demand provider. When present, source terms prefer it over the
  // legacy dense property spans so fixed-grid production hydro does not need
  // population-scale derived staging arrays.
  const HydroCellSourcePropertyProvider* source_property_provider = nullptr;
  const HydroThermodynamicClosure* thermodynamic_closure = nullptr;
  double redshift = 0.0;
  // Authoritative comoving volume of the cells being advanced (uniform per
  // Cartesian patch). HydroCoreSolver overwrites this with
  // HydroPatchGeometry::cell_volume_comoving before invoking source terms, so
  // source ledgers can integrate energy densities into energies (audit F7).
  // Direct callers of HydroSourceTerm::sourceForCell must set it explicitly.
  double cell_volume_comoving = 0.0;
  // Optional composition-transport hook (CHUI Part 2). Null = no recording.
  HydroFaceMassFluxSink* face_mass_flux_sink = nullptr;
  // CHUI Part 2 ISM regime gate. Empty = no gate (legacy). Otherwise indexed
  // by SOLVER row for HydroCoreSolver (the AMR orchestrator receives a
  // dense-row span and builds patch-local spans): 0 = the effective
  // multiphase closure must NOT be applied (resolved regime owns the cell),
  // 1 = the effective closure may apply. Keeps hydro closure, relaxation and
  // thermochemistry mutually exclusive per cell.
  std::span<const std::uint8_t> effective_closure_allowed;
};

[[nodiscard]] inline bool effectiveClosureAllowed(const HydroSourceContext& context, std::size_t cell_index) noexcept {
  return context.effective_closure_allowed.empty() ||
      (cell_index < context.effective_closure_allowed.size() && context.effective_closure_allowed[cell_index] != 0U);
}

// Volume-integrated conserved diagnostics over a selected set of real cells.
// Density-like conserved lanes are multiplied by HydroPatchGeometry::cell_volume_comoving.
// internal_energy is derived as total energy minus kinetic energy, not stored as an
// independent conserved authority.
struct HydroConservationTotals {
  double mass = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;
  double total_energy = 0.0;
  double internal_energy = 0.0;
  double metal_mass = 0.0;

  HydroConservationTotals& operator+=(const HydroConservationTotals& rhs);
  HydroConservationTotals& operator-=(const HydroConservationTotals& rhs);
};

[[nodiscard]] HydroConservationTotals operator+(
    HydroConservationTotals lhs,
    const HydroConservationTotals& rhs);
[[nodiscard]] HydroConservationTotals operator-(
    HydroConservationTotals lhs,
    const HydroConservationTotals& rhs);

struct HydroConservationReport {
  HydroConservationTotals before;
  HydroConservationTotals after;
  HydroConservationTotals flux_delta;
  HydroConservationTotals source_delta;
  HydroConservationTotals floor_delta;
  // Dual energy: volume-integrated thermal correction
  // sum (rho e_selected - (E - K)) V over cells where the auxiliary internal
  // energy (or the thermal floor) replaced the conservative thermal energy.
  // E itself is not modified by the solver; this correction enters the
  // canonical state when the runtime stores e from the auxiliary lane. It is
  // the complete non-conservative energy ledger of the formalism (only the
  // total_energy/internal_energy entries are non-zero).
  HydroConservationTotals dual_energy_sync_delta;
  HydroConservationTotals residual;
  std::uint64_t cell_count = 0;
  std::uint64_t internal_energy_floor_count = 0;
  std::uint64_t dual_energy_sync_count = 0;
};

class HydroSourceTerm {
 public:
  virtual ~HydroSourceTerm() = default;

  [[nodiscard]] virtual HydroConservedState sourceForCell(
      std::size_t cell_index,
      const HydroConservedState& conserved,
      const HydroPrimitiveState& primitive,
      const HydroSourceContext& context) const = 0;
};

class ComovingGravityExpansionSource final : public HydroSourceTerm {
 public:
  [[nodiscard]] HydroConservedState sourceForCell(
      std::size_t cell_index,
      const HydroConservedState& conserved,
      const HydroPrimitiveState& primitive,
      const HydroSourceContext& context) const override;
};

struct HydroProfileEvent {
  std::uint64_t bytes_moved = 0;
  std::uint64_t face_count = 0;
  std::uint64_t limiter_clip_count = 0;
  std::uint64_t positivity_fallback_count = 0;
  std::uint64_t riemann_fallback_count = 0;
  // Total Riemann evaluations (MUSCL-Hancock: 1 per face; CTU: 2 per face).
  std::uint64_t riemann_solve_count = 0;
  std::uint64_t stencil_degradation_count = 0;
  // CTU transverse-correction outcome per face side.
  std::uint64_t ctu_transverse_corrected_states = 0;
  // Side state whose cell had no transverse information (non-periodic
  // ghost / imported MPI ghost): normal-only predictor used.
  std::uint64_t ctu_transverse_degraded_states = 0;
  // Imported AMR/MPI ghost side corrected with the adjacent owner's
  // transverse increment. This one-sided CTU closure preserves temporal
  // centering at a distributed/coarse-fine interface without pretending the
  // remote ghost owns a local transverse-flux stencil. Qualification still
  // requires the boundary convergence tests.
  std::uint64_t ctu_transverse_boundary_proxy_states = 0;
  // Corrected state was inadmissible; normal-only predictor used.
  std::uint64_t ctu_transverse_inadmissible_states = 0;
  // Cells whose final thermal energy came from the auxiliary update.
  std::uint64_t dual_energy_auxiliary_cells = 0;
  // Cells whose auxiliary lane was NaN/non-positive at step start (e.g.
  // imported ghosts or an uninitialized lane): E - K was used as (rho e)^n.
  std::uint64_t dual_energy_unknown_auxiliary_cells = 0;
  // Cells where neither E - K nor the auxiliary value was positive.
  std::uint64_t dual_energy_floor_cells = 0;
  double reconstruct_ms = 0.0;
  double riemann_ms = 0.0;
  double accumulate_ms = 0.0;
  double source_ms = 0.0;
  double total_ms = 0.0;
  HydroConservationReport conservation;
};

struct HydroActiveBatchPolicy {
  // Maximum number of active real cells represented by one reconstruction /
  // Riemann workspace batch. Zero selects the deterministic repository
  // automatic maximum. This controls transient derived state only; all face
  // flux deltas are accumulated against the common step-start state before
  // any authoritative cell is committed.
  std::size_t max_active_cells = 0;
  // Optional sparse per-touched-cell lanes that must be admitted before the
  // advance. Set them from HydroCoreSolver::requiredWorkspaceFeatures() so
  // workspaceBytes/prepareForActiveSet reserve them up front; the batched
  // solver refuses to grow unplanned lanes.
  bool ctu_transverse_workspace = false;
  bool dual_energy_workspace = false;
};

// Time-integration family (orthogonal to reconstruction and Riemann solver).
enum class HydroIntegrationScheme : std::uint8_t {
  // Directionally unsplit MUSCL-Hancock: one Riemann solve per face using
  // the normal-only Hancock predictor of the selected reconstruction.
  kMusclHancock = 0,
  // Corner-transport-upwind (Colella 1990), simplified 3D variant: pass 1
  // solves every face with the normal-predicted states; per-cell transverse
  // conservative increments are accumulated from those fluxes; pass 2 adds
  // half of the transverse increments (all axes except the face axis) to each
  // face state and re-solves. Two Riemann solves per face (6 per cell in 3D,
  // 4 in 2D). Stability: Courant <= 1 in 2D; <= 1/2 in 3D for this variant
  // (no corner-coupled 12-solve terms of Saltzman 1994).
  kCornerTransportUpwind = 1,
};

enum class HydroDualEnergyMode : std::uint8_t {
  kOff = 0,
  // Bryan et al. (1995)/Enzo-style switch: an auxiliary internal-energy lane
  // (HydroConservedStateSoa, initialized from the canonical specific internal
  // energy) supplies the thermal state where E - K < eta E; E is never
  // overwritten by the solver (see HydroCoreSolver).
  kThermalFraction = 1,
};

struct HydroDualEnergyPolicy {
  HydroDualEnergyMode mode = HydroDualEnergyMode::kOff;
  // eta: use the auxiliary internal energy when (E - K) < eta * E after the
  // conservative update, or when E - K is non-positive.
  double thermal_fraction_switch = 1.0e-3;
};

struct HydroNumericsPolicy {
  HydroIntegrationScheme integration_scheme = HydroIntegrationScheme::kMusclHancock;
  HydroDualEnergyPolicy dual_energy{};
};

[[nodiscard]] std::string_view hydroIntegrationSchemeToString(HydroIntegrationScheme scheme);
[[nodiscard]] HydroIntegrationScheme hydroIntegrationSchemeFromString(std::string_view name);
[[nodiscard]] std::string_view hydroDualEnergyModeToString(HydroDualEnergyMode mode);
[[nodiscard]] HydroDualEnergyMode hydroDualEnergyModeFromString(std::string_view name);

struct HydroScratchHighWater {
  std::size_t active_cell_batch_capacity = 0;
  std::size_t face_batch_capacity = 0;
  std::size_t primitive_stencil_capacity = 0;
  std::size_t touched_cell_capacity = 0;
  std::uint64_t capacity_bytes = 0;
};

struct HydroScratchBuffers {
  // Transient per-step workspaces. None of these arrays are restart-authoritative.
  std::vector<HydroConservedState> cell_delta;
  std::vector<HydroPrimitiveState> left_states;
  std::vector<HydroPrimitiveState> right_states;
  std::vector<HydroConservedState> fluxes;
  std::vector<std::size_t> touched_cells;
  std::vector<std::size_t> source_active_cells;
  std::vector<std::size_t> primitive_stencil_cells;
  std::vector<HydroPrimitiveState> primitive_stencil_states;
  std::vector<std::size_t> full_active_cells;
  std::vector<std::size_t> full_active_faces;
  // CTU only: per touched cell, the conservative increment produced by the
  // pass-1 fluxes of the faces normal to x, y and z (same ordering as
  // touched_cells). Transient, batch-scoped, never restart truth.
  std::vector<std::array<HydroConservedState, 3>> ctu_transverse_delta;
  // Dual energy only: per touched cell, the advective internal-energy flux
  // divergence and the face-velocity divergence (both already multiplied by
  // dt/(a V) and the face area), consumed by the auxiliary update.
  std::vector<double> aux_internal_energy_advection;
  std::vector<double> aux_velocity_divergence;
  HydroScratchHighWater high_water{};

  // Admit and physically prepare every vector used by the sparse solver before
  // the first numerical update. The returned capacity is the physical charge;
  // workspaceBytes is the checked upper bound for an initially empty scratch.
  [[nodiscard]] static std::uint64_t workspaceBytes(
      std::size_t active_cell_count, std::size_t active_face_count,
      std::size_t real_cell_count, std::size_t total_cell_count,
      const HydroActiveBatchPolicy& batch_policy);
  [[nodiscard]] std::uint64_t prepareForActiveSet(
      std::size_t active_cell_count, std::size_t active_face_count,
      std::size_t real_cell_count, std::size_t total_cell_count,
      const HydroActiveBatchPolicy& batch_policy);
  void resizeFaceBatch(std::size_t active_cell_batch_count, std::size_t active_face_batch_count);
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const;
};

// Conservative finite-volume update (comoving form). With rho the comoving
// mass density, u the peculiar velocity, x comoving and a the scale factor at
// the flux epoch (HydroUpdateContext::scale_factor, normally the interval
// midpoint),
//   U_i^{n+1} = U_i^n - dt/(a V_i) sum_faces A_f F_f(W_L, W_R) + dt S_i
// with the Hancock/CTU predictors using the same dt/(a dx) factor
// (hydroTransportCourantFactors). Non-expanding runs use a = 1.
class HydroCoreSolver {
 public:
  explicit HydroCoreSolver(double adiabatic_index, HydroNumericsPolicy numerics = {});

  [[nodiscard]] double adiabaticIndex() const;
  [[nodiscard]] const HydroNumericsPolicy& numericsPolicy() const noexcept;
  // Batch policy with the optional workspace lanes this solver needs.
  [[nodiscard]] HydroActiveBatchPolicy requiredWorkspaceFeatures(
      HydroActiveBatchPolicy batch_policy) const noexcept;

  [[nodiscard]] static HydroConservedState conservedFromPrimitive(
      const HydroPrimitiveState& primitive,
      double adiabatic_index);

  [[nodiscard]] static HydroPrimitiveState primitiveFromConserved(
      const HydroConservedState& conserved,
      double adiabatic_index);

  [[nodiscard]] static HydroConservationTotals conservationTotals(
      const HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      std::span<const std::size_t> cells);

  void advancePatch(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroProfileEvent* profile = nullptr,
      HydroFluxRegisterSink* flux_register_sink = nullptr) const;

  void advancePatchWithScratch(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroScratchBuffers& scratch,
      HydroPrimitiveCacheSoa* primitive_cache,
      HydroProfileEvent* profile = nullptr,
      HydroFluxRegisterSink* flux_register_sink = nullptr) const;

  void advancePatchActiveSet(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active_set,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroProfileEvent* profile = nullptr,
      HydroFluxRegisterSink* flux_register_sink = nullptr) const;

  void advancePatchActiveSetWithScratch(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active_set,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroScratchBuffers& scratch,
      HydroPrimitiveCacheSoa* primitive_cache,
      HydroProfileEvent* profile = nullptr,
      HydroFluxRegisterSink* flux_register_sink = nullptr) const;

  void advancePatchActiveSetBatchedWithScratch(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active_set,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroScratchBuffers& scratch,
      HydroPrimitiveCacheSoa* primitive_cache,
      const HydroActiveBatchPolicy& batch_policy,
      HydroProfileEvent* profile = nullptr,
      HydroFluxRegisterSink* flux_register_sink = nullptr) const;

 private:
  double m_adiabatic_index = 5.0 / 3.0;
  HydroNumericsPolicy m_numerics{};
};

}  // namespace cosmosim::hydro
