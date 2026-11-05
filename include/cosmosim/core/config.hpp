#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/particle_species.hpp"

namespace cosmosim::core {

enum class SimulationMode {
  kCosmoCube,
  kZoomIn,
  kIsolatedGalaxy,
  kIsolatedCluster,
};

enum class GravitySolver {
  kTreePm,
};

enum class HydroSolver {
  kGodunovFv,
};

// Hydro numerical-method selectors. These are core-owned configuration
// enums (core must not include hydro headers); workflows map them onto
// cosmosim::hydro types. Defaults preserve the historical PLM/MC/MUSCL-Hancock/HLLC path.
enum class HydroReconstructionMode {
  kPiecewiseConstant,
  kPiecewiseLinear,
  kPiecewiseParabolic,
};

enum class HydroSlopeLimiterMode {
  kMinmod,
  kMonotonizedCentral,
  kVanLeer,
};

enum class HydroIntegrationMode {
  kMusclHancock,
  kCornerTransportUpwind,
};

enum class HydroRiemannSolverMode {
  kHllc,
  kHlle,
};

enum class HydroDualEnergyConfigMode {
  kOff,
  kThermalFraction,
};

// CHUI Part 3: stellar wind / pre-SN mechanical channel.
enum class SspWindModelConfig {
  kNone,
  kTableV1,
};

// CHUI Part 3 Milestone I: gas AMR time stepping. kSynchronized (default)
// advances every level with the coordinator interval; kRecursiveSubcycling
// advances level L with dt / r^(L - L_min) (Berger-Oliger, nested temporal
// histories, per-level reflux).
enum class AmrTimeSteppingConfigMode {
  kSynchronized,
  kRecursiveSubcycling,
};

// CHUI Part 3: production AMR regrid. kOff keeps the historical static
// hierarchy; kAutomatic refines/derefines at synchronized hierarchy
// boundaries every amr_regrid_interval_steps accepted steps.
enum class AmrRegridConfigMode {
  kOff,
  kAutomatic,
};

// A03 force-resolution contract for gas AMR. The reference mode preserves
// historical behavior and records the unqualified regime. The production
// policy limits optional refinement below gas softening and rejects a
// mandatory request that gravity cannot resolve.
enum class AmrForceResolutionPolicy {
  kAllowCountedReference,
  kLimitOptionalRejectMandatory,
};

// V02 research controller. Fixed reference preserves independent chemistry
// and deterministic max-operation AMR policies. Joint v1 coordinates their
// retained state and RAM headroom without suppressing mandatory physics.
enum class V02ControlMode {
  kFixedReference,
  kJointErrorResourceV1,
};

// Node-level memory admission policy. kOff preserves historical semantics:
// only parallel.process_memory_budget_bytes (0 = unlimited) bounds a rank.
enum class NodeMemoryAdmissionMode {
  kOff,
  // Determine a node limit (explicit node budget, cgroup/job limit, physical
  // RAM; minimum of the known values), subtract the node reserve once and
  // give each co-located rank an equal share. Deterministic admission.
  kSharedEqual,
  // Same node limit and single reserve, but one node-shared accounting pool
  // (MPI-3 shared window, lock-free counter): ranks may use unequal amounts,
  // the node-wide sum is bounded. Admission under contention is not
  // run-to-run reproducible.
  kSharedPool,
};

enum class TreePmAssignmentScheme {
  kCic,
  kTsc,
};

enum class TreePmForceGradient {
  // Reference: three spectral derivatives (one inverse FFT per axis).
  kSpectral,
  // Opt-in, unqualified: one inverse FFT plus a 4th-order centered difference.
  kFiniteDifference4,
};

// Long-range PM field reuse between refreshes (numerics.treepm_pm_reuse_policy).
// kExactSourceState (reference): a field is reused only for the identical
//   source state it was built from; treepm_update_cadence_steps must be 1.
// kStaleWithinCadence (experimental, unqualified): with cadence N > 1, the
//   field built at a refresh step is interpolated at current target positions
//   for the next N-1 global KDK force evaluations (first-order in N dt in the
//   long-range force); decomposition epoch, mesh and split must be unchanged.
// kMultirateOuterKick (experimental, unqualified; hierarchical KDK only): an
//   r-RESPA outer PM level. treepm_update_cadence_steps is the TARGET number
//   of coarse blocks per PM cycle; the cycle length T is fixed at opening from
//   that target, the PM displacement criterion and the next output/end
//   boundary, and inner blocks are re-planned by the reference per-block
//   criteria (an interior PM displacement violation rolls the cycle back and
//   retries it shorter). PM is evaluated only at cycle endpoints and applied
//   as two half-cycle kicks (core::HierarchicalPmCyclePlan). Formally second
//   order; adaptive stepping breaks strict reversibility (docs/time_integration.md).
enum class TreePmPmReusePolicy {
  kExactSourceState,
  kStaleWithinCadence,
  kMultirateOuterKick,
};

// Response to a MemoryGovernor refusal (on any rank) of the multirate PM cycle
// rollback snapshot (numerics.treepm_multirate_rollback_memory_policy).
// kSingleBlockFallback (default): every rank opens that cycle as one coarse
//   block (cannot be rejected; all timestep and PM displacement limits still
//   apply). The schedule then depends on the memory budget; each fallback is
//   recorded (report/counters/block records).
// kFail: the run fails collectively instead (reproducibility mode: the
//   numerical schedule never depends on the memory budget).
enum class MultirateRollbackMemoryPolicy {
  kSingleBlockFallback,
  kFail,
};

enum class TreePmOpeningCriterion {
  kGeometric,
  kComDistance,
  kRelativeForceError,
};

// Displacement limit on the hierarchical coarse (PM) interval
// (see docs/time_integration.md).
enum class HierarchicalCoarseDisplacementCriterion {
  kMaxSpeedMesh,
  kRmsSplitScale,
};

// Adaptive TreePM error estimate (see docs/tree_pm_coupling.md).
enum class TreePmAdaptiveErrorEstimator {
  kResidualDerivative,
  kRelativeMonopole,
};

// Adaptive TreePM treatment of nodes straddling the short-range cutoff.
enum class TreePmCutoffNodePolicy {
  kContain,
  kComDistance,
};

enum class PmDecompositionMode {
  kSlab,
  kPencil,
};

enum class CoordinateFrame {
  kComoving,
  kPhysical,
};

enum class ModeHydroBoundary {
  kAuto,
  kPeriodic,
  kOpen,
  kReflective,
};

enum class ModeGravityBoundary {
  kAuto,
  kPeriodic,
  kIsolatedMonopole,
};

enum class ZoomLongRangeStrategy {
  kDisabled,
  kGlobalCoarsePlusFocusedHighResCorrection,
};

enum class InitialConditionConvention {
  kGenerated,
  kChuiCanonicalV1,
  kGadgetArepoBridgeV1,
  kManifestV1,
};

enum class InitialConditionCoordinateFrame {
  kUnspecified,
  kComoving,
  kPhysical,
};

enum class InitialConditionVelocityConvention {
  kUnspecified,
  kPhysicalPeculiar,
  // Legacy/custom convention: stored = sqrt(a) * peculiar, so import multiplies by 1/sqrt(a).
  kSqrtAScaledPeculiar,
  // Standard GADGET/AREPO cosmological snapshot convention: peculiar = stored * sqrt(a).
  kGadgetArepoStoredPeculiar,
  kComovingCoordinateRate,
};

enum class InitialConditionSpeciesPolicy {
  kReject,
  kDarkMatter,
  kStar,
  kBlackHole,
  kTracer,
};

enum class InitialConditionMissingFieldPolicy {
  kReject,
  kReconstruct,
  kUseConfigValue,
  kDialectDefinedDefault,
  kPreserveUnavailable,
};

enum class StarFormationModelKind {
  kLegacySchmidtThreshold,
  kAdaptiveBoundJeans,
  kEffectiveMultiphaseTngLike,
};

enum class EffectiveIsmEosRelaxation {
  kInstantaneous,
  kFiniteTimescale,
};

enum class EffectiveIsmBirthMassConvention {
  kInitialSspMass,
  kLongLivedMass,
};

enum class EffectiveIsmFeedbackCoupling {
  kExternalFeedbackCalibrated,
  kEffectiveEosOnly,
};

enum class StarFormationCollapseTimescale {
  kFreeFall,
  kMinimumFreeFallOrCompression,
};

enum class StarParticleMassPolicy {
  kFixed,
  kGasResolutionFraction,
};

enum class FeedbackMode {
  kThermal,
  kKinetic,
  kMomentum,
  kThermalKineticMomentum,
};

enum class FeedbackVariant {
  kNone,
  kDelayedCooling,
  kStochastic,
};

enum class UvBackgroundModel {
  kNone,
  // Legacy analytic proxies of an overall heating scale (pre-Part-2); refused
  // by every Part 2 thermochemistry backend.
  kHm12,
  kFg20,
  // Production: TREECOOL-format table at physics.uv_background_table_path.
  kTreecoolTable,
};

enum class SelfShieldingModel {
  kNone,
  kRahmati13Like,
};

enum class CoolingModel {
  kPrimordial,
  kPrimordialMetalLine,
};


// CHUI Part 2 baryonic-physics selectors (see docs/configuration.md).
enum class ThermochemistryBackendMode : std::uint8_t {
  kLegacyCooling = 0,
  kAnalyticReference = 1,
  kEquilibriumPrimordial = 2,
  kFixedHybrid = 3,
  kFixedDetailed = 4,
};

enum class MetalCoolingModeConfig : std::uint8_t {
  kNone = 0,
  kFineStructureCold = 1,
  kElementTable = 2,
  kElementTableAndFineStructure = 3,
};

enum class MetalCoolingTableRangeConfig : std::uint8_t {
  kReject = 0,
  kClampDensityAndRedshift = 1,
};

enum class DustModelConfig : std::uint8_t {
  kNone = 0,
  kTwoReservoirV1 = 1,
};

enum class IsmRegimeModelConfig : std::uint8_t {
  kLegacy = 0,
  kEffectiveOrResolved = 1,
};

enum class SgsTurbulenceModelConfig : std::uint8_t {
  kNone = 0,
  kBudgetedV1 = 1,
};

enum class ChemistryAdaptiveSelectionConfig : std::uint8_t {
  kOff = 0,
  kLocalV1 = 1,
};

enum class StellarPopulationModelConfig : std::uint8_t {
  kLegacyTable = 0,
  kChabrierReferenceV1 = 1,
  kChabrierYieldTableV1 = 2,
};

enum class SupernovaFeedbackModelConfig : std::uint8_t {
  kLegacy = 0,
  kMechanicalResolvedUnresolved = 1,
};

enum class MetalSpeciesMode {
  kTotalOnly,
  kCoreElements,
};

enum class MetalDiffusionModel {
  kNone,
  kSmagorinsky,
};

enum class MetalDiffusionTimeIntegrator {
  kExplicitSubcycling,
  kRkl2,
};

enum class IntegratorTimeVariable {
  kScaleFactor,
  kLogScaleFactor,
  kCodeTime,
  kPhysicalTime,
};

struct CosmologyConfig {
  double omega_matter = 0.315;
  double omega_lambda = 0.685;
  double omega_baryon = 0.049;
  double hubble_param = 0.674;
  double sigma8 = 0.811;
  double scalar_index_ns = 0.965;
  double box_size_x_mpc_comoving = 50.0;
  double box_size_y_mpc_comoving = 50.0;
  double box_size_z_mpc_comoving = 50.0;
  // Legacy scalar compatibility lane; canonical runtime geometry must use axis-aware fields above.
  double box_size_mpc_comoving = 50.0;
};

struct NumericsConfig {
  double a_begin = 1.0;
  double a_end = 1.0;
  double z_begin = 0.0;
  double z_end = 0.0;
  double t_code_begin = 0.0;
  double t_code_end = 1.0;
  double t_phys_begin = 0.0;
  double t_phys_end = 0.0;
  IntegratorTimeVariable integrator_time_variable = IntegratorTimeVariable::kScaleFactor;
  double cosmology_max_delta_ln_a = 1.0e-2;
  double cosmology_max_hubble_time_fraction = 1.0e-2;
  double source_max_fractional_change = 0.1;
  double hydro_density_floor_code = 1.0e-10;
  double hydro_pressure_floor_code = 1.0e-10;
  // Maximum real-cell count represented by one hydro reconstruction/Riemann
  // scratch batch. 0 selects the deterministic runtime automatic maximum,
  // further bounded by current MemoryGovernor headroom.
  std::uint64_t hydro_active_batch_max_cells = 0;
  HydroReconstructionMode hydro_reconstruction = HydroReconstructionMode::kPiecewiseLinear;
  HydroSlopeLimiterMode hydro_slope_limiter = HydroSlopeLimiterMode::kMonotonizedCentral;
  HydroIntegrationMode hydro_integration_scheme = HydroIntegrationMode::kMusclHancock;
  HydroRiemannSolverMode hydro_riemann_solver = HydroRiemannSolverMode::kHllc;
  HydroDualEnergyConfigMode hydro_dual_energy = HydroDualEnergyConfigMode::kOff;
  // Dual-energy switch eta in (0, 1): use the auxiliary internal energy when
  // the conservative thermal energy falls below eta times the total energy.
  double hydro_dual_energy_eta = 1.0e-3;
  int max_global_steps = 1024;
  // Zero retains the safe global KDK default. Nonzero selects the provisional
  // synchronized collisionless DMO hierarchy; scientific qualification pending.
  int hierarchical_max_rung = 0;
  HierarchicalCoarseDisplacementCriterion hierarchical_coarse_displacement_criterion =
      HierarchicalCoarseDisplacementCriterion::kMaxSpeedMesh;
  double hierarchical_coarse_rms_displacement_factor = 0.25;
  // Softened-encounter rung criterion dt <= eta sqrt(a^3 (r_nn^2+eps^2)^{3/2} /
  // (G m_max)) from each row's nearest direct pair; bounds the timestep where
  // the acceleration criterion sees |a| -> 0 during a pair crossing.
  bool hierarchical_encounter_timestep_enabled = true;
  // Dimensionless accuracy coefficient of the gravitational timestep
  // criterion dt = eta sqrt(eps a^3 / |A|) (scale-free comoving acceleration A)
  // and of the softened-encounter rung bound dt <= eta sqrt(a^3/T).  Note this
  // is not AREPO's ErrTolIntAccuracy, whose criterion is sqrt(2 eta eps/|a|).
  // 0.2 is the historical hardcoded compatibility value.
  double gravity_timestep_eta = 0.2;
  int amr_max_level = 10;
  // CHUI Part 3 automatic AMR regrid (see docs/configuration.md). A criterion
  // value of 0 disables that criterion; kAutomatic requires at least one.
  AmrRegridConfigMode amr_regrid_mode = AmrRegridConfigMode::kOff;
  AmrForceResolutionPolicy amr_force_resolution_policy =
      AmrForceResolutionPolicy::kLimitOptionalRejectMandatory;
  V02ControlMode v02_control_mode = V02ControlMode::kFixedReference;
  // Fraction of governor headroom exposed to optional regrid candidates in
  // joint mode after reserving observed adaptive-chemistry growth.
  double v02_optional_memory_headroom_fraction = 0.5;
  // CHUI Part 3 collective global-step rollback (rung zero): 0 disables.
  int global_step_max_retries = 0;
  AmrTimeSteppingConfigMode amr_time_stepping = AmrTimeSteppingConfigMode::kSynchronized;
  // CHUI Part 3 (A04 RESEARCH): predictive event-aware refinement of cells
  // hosting stars younger than this age (pending CC SN / wind feedback);
  // 0 disables. Requires amr_regrid_mode = automatic.
  double amr_predictive_feedback_age_years = 0.0;
  // R2 (A04, RESEARCH/HEURISTIC in expectation mode): refine host cells of
  // stars whose canonical event chronology yields a CC/Ia event within this
  // lookahead (also at the drifted position). 0 disables.
  double amr_predictive_event_lookahead_years = 0.0;
  int amr_regrid_interval_steps = 1;
  int amr_regrid_max_operations = 64;
  double amr_refine_cell_mass_code = 0.0;
  double amr_refine_jump_threshold = 0.0;
  double amr_refine_jeans_cells = 4.0;
  double amr_derefine_hysteresis = 0.5;
  // R1 (A03) physics-aware indicators; 0 disables each. MANDATORY: cooling
  // length in cells (legacy in-solver cooling only). OPTIONAL (budgeted):
  // metallicity-interface jump, species-interface jump (composition network
  // backends only), vorticity |curl v| dx / c_s.
  double amr_refine_cooling_length_cells = 0.0;
  double amr_refine_metallicity_jump = 0.0;
  double amr_refine_species_jump = 0.0;
  double amr_refine_vorticity = 0.0;
  double gravity_softening_kpc_comoving = 1.0;
  double gravity_softening_gas_kpc_comoving = -1.0;
  double gravity_softening_dark_matter_kpc_comoving = -1.0;
  double gravity_softening_star_kpc_comoving = -1.0;
  double gravity_softening_black_hole_kpc_comoving = -1.0;
  double gravity_softening_tracer_kpc_comoving = -1.0;
  GravitySolver gravity_solver = GravitySolver::kTreePm;
  HydroSolver hydro_solver = HydroSolver::kGodunovFv;
  int treepm_pm_grid_nx = 16;
  int treepm_pm_grid_ny = 16;
  int treepm_pm_grid_nz = 16;
  // Legacy scalar compatibility lane; canonical runtime geometry must use axis-aware fields above.
  int treepm_pm_grid = 16;
  double treepm_asmth_cells = 1.25;
  double treepm_rcut_cells = 6.25;
  TreePmOpeningCriterion treepm_tree_opening_criterion = TreePmOpeningCriterion::kComDistance;
  double treepm_tree_opening_theta = 0.7;
  double treepm_tree_relative_force_tolerance = 0.005;
  double treepm_tree_relative_force_acceleration_floor = 1.0e-30;
  bool treepm_gaussian_pair_lookup_enabled = false;
  bool treepm_full_residual_lookup_enabled = false;
  bool treepm_leaf_simd_enabled = false;
  bool treepm_radial_coefficients_enabled = false;
  int treepm_residual_block_size = 64;
  int treepm_tree_max_leaf_size = 16;

  bool treepm_full_mac_diagnostics = true;
  bool treepm_adaptive_acceptance_enabled = false;
  double treepm_adaptive_maximum_opening_angle = 0.25;
  TreePmAdaptiveErrorEstimator treepm_adaptive_error_estimator =
      TreePmAdaptiveErrorEstimator::kResidualDerivative;
  TreePmCutoffNodePolicy treepm_cutoff_node_policy = TreePmCutoffNodePolicy::kContain;
  bool treepm_group_walk_enabled = false;
  // 0 keeps legacy group-walk grouping (active targets of one source leaf);
  // N in [2, treepm_residual_block_size] groups up to N tree-order targets,
  // decoupled from treepm_tree_max_leaf_size.
  int treepm_group_walk_target_size = 0;
  bool treepm_identical_source_tree_reuse_enabled = false;
  bool treepm_topology_refit_enabled = false;
  TreePmAssignmentScheme treepm_assignment_scheme = TreePmAssignmentScheme::kTsc;
  TreePmForceGradient treepm_pm_force_gradient = TreePmForceGradient::kSpectral;
  bool treepm_enable_window_deconvolution = true;
  int treepm_update_cadence_steps = 1;
  TreePmPmReusePolicy treepm_pm_reuse_policy = TreePmPmReusePolicy::kExactSourceState;
  MultirateRollbackMemoryPolicy treepm_multirate_rollback_memory_policy =
      MultirateRollbackMemoryPolicy::kSingleBlockFallback;
  // Explicitly test/diagnostic-only. Production workflows fail closed when
  // TreePM is requested without a production FFT backend.
  bool treepm_allow_diagnostic_naive_dft = false;
  PmDecompositionMode treepm_pm_decomposition_mode = PmDecompositionMode::kSlab;
  std::uint64_t treepm_tree_exchange_batch_bytes = 4ULL * 1024ULL * 1024ULL;
  std::uint64_t treepm_pm_exchange_batch_bytes = 16ULL * 1024ULL * 1024ULL;
};

struct PhysicsConfig {
  bool enable_cooling = true;
  bool enable_star_formation = true;
  bool enable_feedback = true;
  bool enable_stellar_evolution = true;
  std::string reionization_model = "hm12";
  UvBackgroundModel uv_background_model = UvBackgroundModel::kHm12;
  SelfShieldingModel self_shielding_model = SelfShieldingModel::kNone;
  CoolingModel cooling_model = CoolingModel::kPrimordial;
  std::string metal_line_table_path;
  double temperature_floor_k = 100.0;
  StarFormationModelKind star_formation_model = StarFormationModelKind::kLegacySchmidtThreshold;
  double sf_density_threshold_code = 10.0;
  double sf_temperature_threshold_k = 1.0e4;
  double sf_min_converging_flow_rate_code = 0.0;
  double sf_epsilon_ff = 0.01;
  double sf_bound_alpha_vir_max = 1.0;
  bool sf_require_converging_flow = true;
  StarFormationCollapseTimescale sf_collapse_timescale =
      StarFormationCollapseTimescale::kFreeFall;
  double sf_jeans_mass_floor_code = 100.0;
  StarParticleMassPolicy sf_star_particle_mass_policy =
      StarParticleMassPolicy::kFixed;
  double sf_target_star_particle_mass_code = 0.1;
  double sf_target_star_particle_mass_fraction = 0.25;
  double sf_min_star_particle_mass_code = 0.1;
  double sf_max_star_particle_mass_code = 1.0e30;
  std::uint32_t sf_max_spawn_particles_per_cell_step = 8;
  double sf_max_fractional_mass_conversion = 0.25;
  double sf_min_remaining_gas_fraction = 0.01;
  double sf_min_remaining_gas_mass_code = 1.0e-12;
  double sf_temperature_safety_ceiling_k = 0.0;
  bool sf_stochastic_spawning = true;
  std::uint64_t sf_random_seed = 123456789ull;

  // SH03/TNG-inspired unresolved multiphase-ISM closure. Values are a
  // versioned CHUI reference parameter set, not an exact IllustrisTNG match.
  std::string sf_effective_parameter_set = "chui_sh03_tng_like_v1";
  double sf_effective_n_h_threshold_cgs = 0.13;
  double sf_effective_min_baryon_overdensity = 0.0;
  double sf_effective_t_star_at_threshold_code = 1.0;
  double sf_effective_evaporation_factor_at_threshold = 1000.0;
  double sf_effective_evaporation_density_exponent = -0.8;
  double sf_effective_cold_phase_temperature_k = 1000.0;
  double sf_effective_supernova_specific_energy_code = 1.0e5;
  double sf_effective_massive_star_fraction = 0.1;
  double sf_effective_q_eos = 0.3;
  double sf_effective_isothermal_temperature_k = 1.0e4;
  double sf_effective_hot_excess_tolerance = 1.0;
  EffectiveIsmEosRelaxation sf_effective_eos_relaxation =
      EffectiveIsmEosRelaxation::kInstantaneous;
  double sf_effective_eos_relaxation_timescale_code = 0.0;
  std::uint32_t sf_effective_eos_table_bins = 256;
  double sf_effective_eos_max_density_ratio = 1.0e6;
  EffectiveIsmBirthMassConvention sf_effective_birth_mass_convention =
      EffectiveIsmBirthMassConvention::kInitialSspMass;
  EffectiveIsmFeedbackCoupling sf_effective_feedback_coupling =
      EffectiveIsmFeedbackCoupling::kExternalFeedbackCalibrated;

  FeedbackMode fb_mode = FeedbackMode::kThermalKineticMomentum;
  FeedbackVariant fb_variant = FeedbackVariant::kNone;
  bool fb_use_returned_mass_budget = true;
  double fb_epsilon_thermal = 0.6;
  double fb_epsilon_kinetic = 0.3;
  double fb_epsilon_momentum = 0.1;
  double fb_sn_energy_erg_per_mass_code = 1.0e49;
  double fb_momentum_code_per_mass_code = 3.0e3;
  std::uint32_t fb_neighbor_count = 8;
  double fb_delayed_cooling_time_code = 0.0;
  double fb_stochastic_event_probability = 0.25;
  std::uint64_t fb_random_seed = 42424242ull;
  std::string stellar_evolution_table_path;
  double stellar_evolution_hubble_time_years = 1.44e10;
  bool stellar_evolution_require_production_table = false;

  // Total metal mass remains authoritative. Element mode is reserved until a
  // complete sidecar/yield/cooling schema is selected; validation fails closed.
  MetalSpeciesMode metal_species_mode = MetalSpeciesMode::kTotalOnly;
  bool enable_metal_diffusion = false;
  MetalDiffusionModel metal_diffusion_model = MetalDiffusionModel::kNone;
  MetalDiffusionTimeIntegrator metal_diffusion_time_integrator =
      MetalDiffusionTimeIntegrator::kExplicitSubcycling;
  double metal_diffusion_coefficient = 0.05;
  double metal_diffusion_cfl = 0.4;
  std::uint32_t metal_diffusion_max_subcycles = 128;
  std::uint32_t metal_diffusion_max_rkl_stages = 64;
  double metal_diffusion_coefficient_floor_code = 0.0;
  double metal_diffusion_coefficient_ceiling_code = 1.0e30;

  bool enable_black_hole_agn = false;
  // Seeding requires an authoritative halo/candidate provider. The reference
  // workflow currently has none, so this is explicit and fail-closed.
  bool bh_enable_seeding = false;
  double bh_seed_halo_mass_threshold_code = 1.0e3;
  double bh_seed_mass_code = 1.0;
  std::uint32_t bh_seed_max_per_cell = 1;
  double bh_alpha_bondi = 1.0;
  bool bh_use_eddington_cap = true;
  double bh_epsilon_r = 0.1;
  double bh_epsilon_f = 0.05;
  double bh_feedback_coupling_efficiency = 1.0;
  double bh_duty_cycle_active_edd_ratio_threshold = 0.01;
  double bh_proton_mass_si = 1.67262192369e-27;
  double bh_thomson_cross_section_si = 6.6524587321e-29;
  double bh_newton_g_si = 6.67430e-11;
  double bh_speed_of_light_si = 2.99792458e8;
  bool enable_tracers = false;
  bool tracer_track_mass = true;
  double tracer_min_host_mass_code = 0.0;

  // ---- CHUI Part 2 baryonic physics (defaults keep Part 1 behaviour) ----
  ThermochemistryBackendMode thermochemistry_backend = ThermochemistryBackendMode::kLegacyCooling;
  MetalCoolingModeConfig metal_cooling_mode = MetalCoolingModeConfig::kNone;
  MetalCoolingTableRangeConfig metal_cooling_table_range = MetalCoolingTableRangeConfig::kReject;
  DustModelConfig dust_model = DustModelConfig::kNone;
  IsmRegimeModelConfig ism_regime_model = IsmRegimeModelConfig::kLegacy;
  SgsTurbulenceModelConfig sgs_turbulence_model = SgsTurbulenceModelConfig::kNone;
  ChemistryAdaptiveSelectionConfig chemistry_adaptive_selection = ChemistryAdaptiveSelectionConfig::kOff;
  StellarPopulationModelConfig stellar_population_model = StellarPopulationModelConfig::kLegacyTable;
  SupernovaFeedbackModelConfig fb_supernova_model = SupernovaFeedbackModelConfig::kLegacy;
  double thermochemistry_relative_tolerance = 1.0e-4;
  double thermochemistry_abundance_absolute_tolerance = 1.0e-12;
  double thermochemistry_analytic_relaxation_time_years = 1.0e6;
  double thermochemistry_analytic_equilibrium_temperature_k = 1.0e4;
  double isrf_g0_habing = 0.0;
  double chemistry_shielding_length_max_pc = 100.0;
  double metal_cooling_switch_temperature_k = 1.0e4;
  double chemistry_primordial_helium_mass_fraction = 0.2449;
  double dust_accretion_timescale_ref_years = 2.0e7;
  double dust_sputtering_timescale_ref_years = 1.7e8;
  double dust_sn_swept_mass_msun = 6800.0;
  double dust_sn_destruction_efficiency = 0.3;
  double dust_agb_condensation_efficiency = 0.2;
  double dust_cc_condensation_efficiency = 0.1;
  double ism_jeans_resolution_cells = 4.0;
  double ism_regime_hysteresis_fraction = 0.2;
  double ism_regime_min_dwell_years = 1.0e6;
  double sgs_dissipation_coefficient = 1.58;
  double sgs_eddy_viscosity_coefficient = 0.05;
  double sgs_feedback_energy_fraction = 0.0;
  double chemistry_adaptive_hysteresis_fraction = 0.3;
  double chemistry_adaptive_min_dwell_years = 1.0e6;
  double chemistry_adaptive_cold_temperature_k = 3.0e3;
  // R3: promote to detailed chemistry when t_cool < ratio x interval (0 off).
  double chemistry_adaptive_cooling_time_ratio = 0.0;
  // C03: normalized previous accepted solver error and conservation-residual
  // fraction thresholds (0 disables each history criterion).
  double chemistry_adaptive_error_estimate_threshold = 0.0;
  double chemistry_adaptive_residual_fraction_threshold = 0.0;
  double chemistry_adaptive_dense_nh_cm3 = 1.0;
  double chemistry_shadow_sample_fraction = 0.0;
  double ssp_imf_mass_min_msun = 0.1;
  double ssp_imf_mass_max_msun = 100.0;
  double ssp_cc_mass_min_msun = 8.0;
  double ssp_cc_mass_max_msun = 40.0;
  double ssp_snia_per_msun = 1.3e-3;
  double ssp_snia_delay_min_years = 4.0e7;
  double ssp_snia_delay_max_years = 1.37e10;
  double ssp_energy_per_sn_erg = 1.0e51;
  double ssp_reference_net_metal_yield = 0.0;
  double fb_presn_wind_momentum_lsun_per_c = 0.0;
  std::uint32_t thermochemistry_max_substeps = 4096U;
  std::uint32_t thermochemistry_max_newton_iterations = 12U;
  std::uint32_t chemistry_shadow_max_cells = 0U;
  std::string uv_background_table_path = "";
  std::string uv_background_table_sha256 = "";
  std::string metal_cooling_table_path = "";
  std::string metal_cooling_table_sha256 = "";
  // C02 fixed detailed metal-ion reference: "none" or a comma-separated
  // subset of C,O,Si,Fe (normalized to catalog order). Requires
  // thermochemistry_backend=fixed_detailed and an external rate dataset
  // (schema chui_metal_ion_rates_v1) bound by SHA-256; no rates are bundled.
  std::string chemistry_metal_ion_elements = "none";
  std::string metal_ion_rate_table_path = "";
  std::string metal_ion_rate_table_sha256 = "";
  std::string ssp_yield_table_path = "";
  std::string ssp_yield_table_sha256 = "";
  bool chemistry_metallicity_scaled_dust_proxy = false;
  // Explicit opt-in to UNQUALIFIED (unverified) rate coefficients; see
  // docs/chemistry_rate_provenance.md. Default refuses those backends.
  bool thermochemistry_allow_unqualified_rates = false;
  // Part 2.5: explicit acceptance of a metal cooling table whose embedded UV
  // background does not match the configured CHUI UVB (see
  // MetalCoolingTable::checkUvBackgroundCompatibility). Default refuses.
  bool metal_cooling_accept_uv_background_mismatch = false;
  bool ssp_stochastic_sn_events = true;
  bool ssp_clustered_associations = false;
  // CHUI Part 3 (docs/configuration.md): correlated associations (RESEARCH)
  // and the table-driven wind/pre-SN mechanical channel.
  double ssp_association_radius_code = 0.0;
  double ssp_association_time_window_years = 3.0e6;
  SspWindModelConfig ssp_wind_model = SspWindModelConfig::kNone;
  std::string ssp_wind_table_path;
  std::string ssp_wind_table_sha256;
};

enum class ScienceSnapshotLayout : std::uint8_t {
  kAuto = 0,
  kSingle = 1,
  kAggregated = 2,
  kSharded = 3,
};

// Write-time verification of a non-collective science snapshot member.
// kFullReadback materializes the member and compares every field against the
// authoritative state (strict validation). kStructural verifies schema,
// file kind, header cosmology/time, provenance hash and per-type counts from
// header metadata only; the root-rank SHA-256 completion marker and the
// independent snapshot-set validation still run in both modes. Collective
// single-file snapshots always use the bounded rank-partition readback.
enum class SnapshotVerificationPolicy : std::uint8_t {
  kFullReadback = 0,
  kStructural = 1,
};

struct OutputConfig {
  std::string run_name = "cosmosim_run";
  std::string output_directory = "outputs";
  std::string output_stem = "snapshot";
  std::string restart_stem = "restart";
  int snapshot_interval_steps = 64;
  // Zero disables code-time cadence.  When positive, the workflow clips a
  // step at each ordered event and persists the next event in restart state.
  double snapshot_interval_time_code = 0.0;
  ScienceSnapshotLayout snapshot_layout = ScienceSnapshotLayout::kAuto;
  // Zero selects backend policy. This is intentionally independent of MPI
  // world size so future aggregated science I/O is not rank-count coupled.
  std::uint32_t snapshot_num_files = 0U;
  bool write_restarts = true;
  // 0 retains the legacy coupling: with write_restarts=true every science
  // snapshot event also writes a checkpoint. N > 0 selects an independent
  // checkpoint cadence every N completed global steps (coarse blocks under
  // hierarchical KDK); snapshot events then no longer imply checkpoints.
  // The next checkpoint index is derived from the restart-persisted step
  // index, so resumed and uninterrupted runs share one checkpoint timeline.
  int restart_interval_steps = 0;
  // Zero disables. Positive: checkpoint at the first restart-safe boundary at
  // or after each multiple of this code-time interval from t_code_begin
  // (no step clipping; checkpoints are not science products). Combines with
  // restart_interval_steps (either triggers) and implies the independent
  // checkpoint policy.
  double restart_interval_time_code = 0.0;
  SnapshotVerificationPolicy snapshot_verification = SnapshotVerificationPolicy::kFullReadback;
};


struct AnalysisConfig {
  enum class DiagnosticsExecutionPolicy : std::uint8_t {
    kRunHealthOnly = 0,
    kRunHealthAndLightScience = 1,
    kAllIncludingProvisional = 2,
  };

  bool enable_diagnostics = true;
  bool enable_halo_workflow = false;
  bool halo_on_the_fly = false;
  DiagnosticsExecutionPolicy diagnostics_execution_policy =
      DiagnosticsExecutionPolicy::kRunHealthAndLightScience;
  int run_health_interval_steps = 1;
  int science_light_interval_steps = 8;
  int science_heavy_interval_steps = 64;
  int retention_bundle_count = 8;
  int power_spectrum_mesh_n = 16;
  int power_spectrum_bin_count = 12;
  int sf_history_bin_count = 16;
  int quicklook_grid_n = 32;
  std::string diagnostics_stem = "diagnostics";
  std::string halo_catalog_stem = "halo_catalog";
  std::string merger_tree_stem = "merger_tree_plan";
  double halo_fof_linking_length_factor = 0.2;
  int halo_fof_min_group_size = 16;
  bool halo_include_gas = true;
  bool halo_include_stars = true;
  bool halo_include_black_holes = true;
};

struct ParallelConfig {
  int mpi_ranks_expected = 1;
  // OpenMP team request: 0 uses the runtime default, 1 requests serial execution,
  // and N > 1 requests N threads when this binary was built with OpenMP.
  // Builds without OpenMP reject requests above one during validation.
  int omp_threads = 1;
  int gpu_devices = 0;
  bool deterministic_reduction = true;
  double decomposition_particle_count_weight = 1.0;
  double decomposition_gas_cell_weight = 1.5;
  double decomposition_tree_interaction_weight = 1.0;
  double decomposition_pm_mesh_weight = 0.25;
  double decomposition_amr_patch_weight = 1.0;
  double decomposition_active_fraction_weight = 2.0;
  double decomposition_memory_pressure_weight = 1.0 / (1024.0 * 1024.0);
  double decomposition_gpu_occupancy_weight = 0.0;
  double decomposition_generic_work_weight = 0.5;
  bool decomposition_runtime_rebalance_enabled = true;
  bool decomposition_spatial_work_enabled = false;
  bool decomposition_debug_exact_ownership_audit = false;
  double decomposition_rebalance_imbalance_trigger = 1.25;
  double decomposition_rebalance_memory_trigger = 1.50;
  double decomposition_rebalance_max_migrated_load_fraction = 0.25;
  double decomposition_measured_tree_pair_weight = 1.0;
  double decomposition_measured_pm_cell_weight = 1.0;
  double decomposition_measured_amr_cell_weight = 1.0;
  double decomposition_measured_hydro_face_weight = 1.0;
  double decomposition_measured_wall_ms_weight = 1.0;
  std::uint64_t isolated_pm_root_workspace_limit_bytes = 256ULL * 1024ULL * 1024ULL;
  std::uint64_t zoom_high_res_allgather_limit_bytes = 256ULL * 1024ULL * 1024ULL;
  // Optional per-rank gravity peak budget. Zero disables the budget gate but
  // never disables the estimator/report.
  std::uint64_t gravity_memory_budget_bytes = 0U;
  // Additional reserve for backend/library allocations not observable through
  // owner containers plus a fractional margin applied to the complete pre-run
  // gravity estimate. Defaults preserve the historical budget behavior.
  std::uint64_t gravity_backend_unknown_reserve_bytes = 0U;
  double gravity_memory_safety_margin_fraction = 0.0;
  // Authoritative per-rank total-process budget for production DMO preflight.
  // Zero preserves historical parameter-file behavior and disables only this
  // process-wide gate; the estimate/report is still produced. The legacy
  // gravity_memory_budget_bytes remains a gravity-subsystem compatibility gate.
  std::uint64_t process_memory_budget_bytes = 0U;
  // Explicit reserves for library/runtime memory that CHUI does not own and
  // cannot calculate exactly from container extents.
  std::uint64_t process_mpi_unknown_reserve_bytes = 0U;
  std::uint64_t process_fftw_unknown_reserve_bytes = 0U;
  std::uint64_t process_hdf5_unknown_reserve_bytes = 0U;
  std::uint64_t process_allocator_unknown_reserve_bytes = 0U;
  // Owned output/restart staging that is permitted to overlap resident
  // canonical/PM state. IC-import staging is deliberately excluded once the
  // integrator owns the canonical state because those lifetimes do not overlap.
  std::uint64_t process_output_restart_overlap_bytes = 0U;
  double process_memory_safety_margin_fraction = 0.0;
  // Node-aware admission (see NodeMemoryAdmissionMode). The node budget is an
  // explicit node-wide byte limit (0 = discover); the node reserve is
  // subtracted exactly once at node level before the per-rank share.
  NodeMemoryAdmissionMode node_memory_admission = NodeMemoryAdmissionMode::kOff;
  std::uint64_t node_memory_budget_bytes = 0U;
  std::uint64_t node_memory_reserve_bytes = 0U;
};

struct UnitsConfig {
  std::string length_unit = "mpc";
  std::string mass_unit = "msun";
  std::string velocity_unit = "km_s";
  CoordinateFrame coordinate_frame = CoordinateFrame::kComoving;
};

struct ModeConfig {
  SimulationMode mode = SimulationMode::kZoomIn;
  std::string ic_file = "generated";
  InitialConditionConvention ic_convention =
      InitialConditionConvention::kGenerated;
  std::string ic_manifest_file;
  // A generic GADGET/AREPO bridge is scientifically valid only when this
  // complete source convention is supplied. Zero/NaN/unspecified defaults are
  // deliberate fail-closed sentinels, not physical defaults.
  double ic_bridge_source_length_unit_to_si = 0.0;
  double ic_bridge_source_mass_unit_to_si = 0.0;
  double ic_bridge_source_velocity_unit_to_si = 0.0;
  InitialConditionCoordinateFrame ic_bridge_coordinate_frame =
      InitialConditionCoordinateFrame::kUnspecified;
  InitialConditionVelocityConvention ic_bridge_velocity_convention =
      InitialConditionVelocityConvention::kUnspecified;
  double ic_bridge_length_hubble_exponent =
      std::numeric_limits<double>::quiet_NaN();
  double ic_bridge_length_scale_factor_exponent =
      std::numeric_limits<double>::quiet_NaN();
  double ic_bridge_mass_hubble_exponent =
      std::numeric_limits<double>::quiet_NaN();
  double ic_bridge_mass_scale_factor_exponent =
      std::numeric_limits<double>::quiet_NaN();
  double ic_bridge_velocity_hubble_exponent =
      std::numeric_limits<double>::quiet_NaN();
  double ic_bridge_velocity_scale_factor_exponent =
      std::numeric_limits<double>::quiet_NaN();
  std::uint64_t ic_chunk_particle_count = 65536ULL;
  std::uint64_t ic_staging_particle_count = 65536ULL;
  InitialConditionSpeciesPolicy ic_part_type2_policy =
      InitialConditionSpeciesPolicy::kReject;
  InitialConditionSpeciesPolicy ic_part_type3_policy =
      InitialConditionSpeciesPolicy::kReject;
  InitialConditionMissingFieldPolicy ic_gas_internal_energy_policy =
      InitialConditionMissingFieldPolicy::kReject;
  InitialConditionMissingFieldPolicy ic_gas_density_policy =
      InitialConditionMissingFieldPolicy::kReject;
  InitialConditionMissingFieldPolicy ic_star_formation_time_policy =
      InitialConditionMissingFieldPolicy::kReject;
  InitialConditionMissingFieldPolicy ic_star_initial_mass_policy =
      InitialConditionMissingFieldPolicy::kReject;
  InitialConditionMissingFieldPolicy ic_star_metallicity_policy =
      InitialConditionMissingFieldPolicy::kReject;
  InitialConditionMissingFieldPolicy ic_bh_mdot_policy =
      InitialConditionMissingFieldPolicy::kReject;
  double ic_gas_internal_energy_value_code = 0.0;
  double ic_gas_density_value_code = 0.0;
  double ic_star_formation_time_value = 0.0;
  double ic_star_initial_mass_value_code = 0.0;
  double ic_star_metallicity_value = 0.0;
  double ic_bh_mdot_value_code = 0.0;
  bool zoom_high_res_region = false;
  std::string zoom_region_file;
  ZoomLongRangeStrategy zoom_long_range_strategy =
      ZoomLongRangeStrategy::kDisabled;
  double zoom_region_center_x_mpc_comoving = 0.5;
  double zoom_region_center_y_mpc_comoving = 0.5;
  double zoom_region_center_z_mpc_comoving = 0.5;
  double zoom_region_radius_mpc_comoving = 0.0;
  int zoom_focused_pm_grid_nx = 0;
  int zoom_focused_pm_grid_ny = 0;
  int zoom_focused_pm_grid_nz = 0;
  double zoom_contamination_radius_mpc_comoving = 0.0;
  ModeHydroBoundary hydro_boundary = ModeHydroBoundary::kAuto;
  ModeGravityBoundary gravity_boundary = ModeGravityBoundary::kAuto;
};

struct CompatibilityConfig {
  bool allow_unknown_keys = false;
};

struct SimulationConfig final {
 private:
  struct ConstructionToken {};
  explicit SimulationConfig(ConstructionToken) {}

 public:
  SimulationConfig() = delete;
  SimulationConfig(const SimulationConfig&) = default;
  SimulationConfig(SimulationConfig&&) noexcept = default;
  SimulationConfig& operator=(const SimulationConfig&) = default;
  SimulationConfig& operator=(SimulationConfig&&) noexcept = default;
  ~SimulationConfig() = default;

  int schema_version = 1;
  CosmologyConfig cosmology;
  NumericsConfig numerics;
  PhysicsConfig physics;
  OutputConfig output;
  ParallelConfig parallel;
  AnalysisConfig analysis;
  UnitsConfig units;
  ModeConfig mode;
  CompatibilityConfig compatibility;

 private:
  friend struct FrozenConfig;
  friend SimulationConfig makeUnvalidatedSimulationConfigForTests();
};

using NormalizedConfig = SimulationConfig;

[[nodiscard]] SimulationConfig makeUnvalidatedSimulationConfigForTests();

struct ProvenanceMetadata {
  std::string source_name;
  std::uint64_t config_hash = 0;
  std::string config_hash_hex;
  std::vector<std::string> deprecation_warnings;
};

struct UserConfigEntry {
  std::string canonical_key;
  std::string value;
  std::string source_key;
  int source_line = 0;
  bool from_alias = false;
};

struct UserConfig {
  std::string source_name;
  std::vector<UserConfigEntry> entries;
  std::vector<std::string> alias_resolution_notes;
};

struct FrozenConfig {
  FrozenConfig();

  UserConfig user_config;
  NormalizedConfig config;
  std::string raw_text;
  std::string normalized_text;
  ProvenanceMetadata provenance;
};

struct DerivedRuntimeConfig {
  double t_code_begin = 0.0;
  double t_code_end = 1.0;
  std::array<double, 3> box_size_mpc_comoving{};
  std::array<int, 3> treepm_pm_grid_shape{};
  std::array<double, k_particle_species_count> gravity_softening_kpc_comoving_by_species{};
  std::uint64_t normalized_config_hash = 0;
  std::string normalized_config_hash_hex;
};

struct ParseOptions {
  bool allow_unknown_keys = false;
};

class ConfigError : public std::runtime_error {
 public:
  explicit ConfigError(const std::string& message);
};

[[nodiscard]] FrozenConfig loadFrozenConfigFromFile(
    const std::filesystem::path& path,
    const ParseOptions& options = {});

[[nodiscard]] FrozenConfig loadFrozenConfigFromString(
    const std::string& config_text,
    const std::string& source_name,
    const ParseOptions& options = {});

[[nodiscard]] DerivedRuntimeConfig deriveRuntimeConfig(const FrozenConfig& frozen_config);
[[nodiscard]] std::string serializeDerivedRuntimeConfig(const DerivedRuntimeConfig& derived_config);

void writeNormalizedConfigSnapshot(
    const FrozenConfig& frozen_config,
    const std::filesystem::path& run_directory);

[[nodiscard]] std::string modeToString(SimulationMode mode);
[[nodiscard]] std::string gravitySolverToString(GravitySolver solver);
[[nodiscard]] std::string hydroSolverToString(HydroSolver solver);
[[nodiscard]] std::string hydroReconstructionModeToString(HydroReconstructionMode mode);
[[nodiscard]] std::string hydroSlopeLimiterModeToString(HydroSlopeLimiterMode mode);
[[nodiscard]] std::string hydroIntegrationModeToString(HydroIntegrationMode mode);
[[nodiscard]] std::string hydroRiemannSolverModeToString(HydroRiemannSolverMode mode);
[[nodiscard]] std::string hydroDualEnergyConfigModeToString(HydroDualEnergyConfigMode mode);
[[nodiscard]] std::string amrRegridConfigModeToString(AmrRegridConfigMode mode);
[[nodiscard]] std::string amrForceResolutionPolicyToString(AmrForceResolutionPolicy mode);
[[nodiscard]] std::string v02ControlModeToString(V02ControlMode mode);
[[nodiscard]] std::string amrTimeSteppingConfigModeToString(AmrTimeSteppingConfigMode mode);
[[nodiscard]] std::string sspWindModelConfigToString(SspWindModelConfig mode);
[[nodiscard]] std::string nodeMemoryAdmissionModeToString(NodeMemoryAdmissionMode mode);
[[nodiscard]] std::string coordinateFrameToString(CoordinateFrame frame);
[[nodiscard]] std::string initialConditionConventionToString(
    InitialConditionConvention convention);
[[nodiscard]] std::string initialConditionCoordinateFrameToString(
    InitialConditionCoordinateFrame frame);
[[nodiscard]] std::string initialConditionVelocityConventionToString(
    InitialConditionVelocityConvention convention);
[[nodiscard]] std::string initialConditionSpeciesPolicyToString(
    InitialConditionSpeciesPolicy policy);
[[nodiscard]] std::string initialConditionMissingFieldPolicyToString(
    InitialConditionMissingFieldPolicy policy);
[[nodiscard]] std::string modeHydroBoundaryToString(ModeHydroBoundary boundary);
[[nodiscard]] std::string modeGravityBoundaryToString(ModeGravityBoundary boundary);
[[nodiscard]] std::string starFormationModelKindToString(StarFormationModelKind model);
[[nodiscard]] std::string starFormationCollapseTimescaleToString(StarFormationCollapseTimescale mode);
[[nodiscard]] std::string starParticleMassPolicyToString(StarParticleMassPolicy policy);
[[nodiscard]] std::string effectiveIsmEosRelaxationToString(EffectiveIsmEosRelaxation mode);
[[nodiscard]] std::string effectiveIsmBirthMassConventionToString(EffectiveIsmBirthMassConvention mode);
[[nodiscard]] std::string effectiveIsmFeedbackCouplingToString(EffectiveIsmFeedbackCoupling mode);
[[nodiscard]] std::string feedbackModeToString(FeedbackMode mode);
[[nodiscard]] std::string feedbackVariantToString(FeedbackVariant variant);
[[nodiscard]] std::string uvBackgroundModelToString(UvBackgroundModel model);
[[nodiscard]] std::string selfShieldingModelToString(SelfShieldingModel model);
[[nodiscard]] std::string coolingModelToString(CoolingModel model);
[[nodiscard]] std::string metalSpeciesModeToString(MetalSpeciesMode mode);
[[nodiscard]] std::string thermochemistryBackendModeToString(ThermochemistryBackendMode value);
[[nodiscard]] std::string metalCoolingModeConfigToString(MetalCoolingModeConfig value);
[[nodiscard]] std::string metalCoolingTableRangeConfigToString(MetalCoolingTableRangeConfig value);
[[nodiscard]] std::string dustModelConfigToString(DustModelConfig value);
[[nodiscard]] std::string ismRegimeModelConfigToString(IsmRegimeModelConfig value);
[[nodiscard]] std::string sgsTurbulenceModelConfigToString(SgsTurbulenceModelConfig value);
[[nodiscard]] std::string chemistryAdaptiveSelectionConfigToString(ChemistryAdaptiveSelectionConfig value);
[[nodiscard]] std::string stellarPopulationModelConfigToString(StellarPopulationModelConfig value);
[[nodiscard]] std::string supernovaFeedbackModelConfigToString(SupernovaFeedbackModelConfig value);

[[nodiscard]] std::string metalDiffusionModelToString(MetalDiffusionModel model);
[[nodiscard]] std::string metalDiffusionTimeIntegratorToString(
    MetalDiffusionTimeIntegrator integrator);
[[nodiscard]] std::string integratorTimeVariableToString(IntegratorTimeVariable variable);

}  // namespace cosmosim::core
