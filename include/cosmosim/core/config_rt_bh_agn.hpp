#pragma once

// CHUI-RT-BH-AGN-REF-v1, shared contract B04: typed, versioned options for
// radiative transfer, the advanced black-hole engine and AGN feedback.
//
// Keys live in the authoritative `.param.txt` [physics] section with the
// spellings proposed by the master design section 7 (physics.rt_*,
// physics.bh_*, physics.agn_*); see docs/configuration.md. Disabled defaults
// reproduce the predecessor behaviour exactly: the normalized dump emits these
// blocks only when the corresponding model is enabled (config-hash and restart
// admission of existing runs unchanged) and non-default values on a disabled
// model are rejected as inert keys.
//
// Compatibility mapping (amendment A1): physics.enable_black_hole_agn keeps its
// meaning. physics.bh_engine_model = legacy_host_cell (default) runs the
// unchanged BlackHoleAgnModel. The reference engine reuses bh_alpha_bondi
// (Bondi alpha), bh_epsilon_r (radiative efficiency), bh_use_eddington_cap
// (reference inflow cap) and the bh_*_si constants; bh_epsilon_f and
// bh_feedback_coupling_efficiency are legacy-only and rejected when present
// together with the reference engine.

#include <cstdint>
#include <string>

namespace cosmosim::core {

// ---------------------------------------------------------------- radiation
enum class RtTransportConfig : std::uint8_t {
  kM1 = 0,
  kM1DirectRays = 1,  // P1: implemented, rejected at admission (unqualified)
  kVet = 2,           // P1: implemented (grey), rejected at admission
};
enum class RtIntegratorConfig : std::uint8_t {
  kImexPdArs = 0,          // He, Wibking & Krumholz (2024) eqs. 28-29
  kExplicitReference = 1,  // explicit comparator; NOT the thick-limit path
};
enum class RtSpectralProfileConfig : std::uint8_t {
  kGrey = 0,
  kIonizing3 = 1,
  kBroadUvIr = 2,
  kAgnHard = 3,  // P1: rejected at admission (no qualified hard microphysics)
};
enum class RtRiemannSolverConfig : std::uint8_t { kHll = 0, kGlf = 1 };
enum class RtReconstructionConfig : std::uint8_t { kPlmMc = 0, kFirstOrder = 1 };
enum class RtOpacityModelConfig : std::uint8_t {
  kPhotoionization = 0,    // sigma_g n_species from chemistry species state
  kConstantReference = 1,  // REFERENCE ONLY (idealized tests)
};
enum class RtBoundaryConfig : std::uint8_t { kPeriodic = 0, kOutflow = 1, kReflecting = 2 };
enum class RtSourceShapeConfig : std::uint8_t {
  kSedTable = 0,            // production: versioned SED tables bound by SHA-256
  kBlackbodyReference = 1,  // REFERENCE ONLY analytic shape (tests/idealized)
  kPowerLawReference = 2,   // REFERENCE ONLY analytic shape (tests/idealized)
};
enum class RtUvbPolicyConfig : std::uint8_t {
  kRejectOverlap = 0,
  kBackgroundPlusLocalSources = 1,  // UVB in chemistry, RT carries local sources only
};

struct RadiationTransportConfig {
  bool enabled = false;
  RtTransportConfig transport = RtTransportConfig::kM1;
  RtIntegratorConfig integrator = RtIntegratorConfig::kImexPdArs;
  RtSpectralProfileConfig spectral_profile = RtSpectralProfileConfig::kIonizing3;
  double spectrum_upper_edge_ev = 0.0;  // explicit upper edge (required, > last profile edge)
  double c_fraction = 1.0;              // c_hat / c in (0, 1]
  double cfl = 0.3;
  std::uint32_t max_subcycles = 4096U;
  RtRiemannSolverConfig riemann_solver = RtRiemannSolverConfig::kHll;
  RtReconstructionConfig reconstruction = RtReconstructionConfig::kPlmMc;
  double newton_relative_tolerance = 1.0e-10;
  std::uint32_t newton_max_iterations = 32U;
  std::uint32_t max_retries = 8U;
  RtOpacityModelConfig opacity_model = RtOpacityModelConfig::kPhotoionization;
  double constant_absorption_cm2_g = 0.0;
  double constant_scattering_cm2_g = 0.0;
  bool allow_reference_models = false;
  RtBoundaryConfig boundary = RtBoundaryConfig::kPeriodic;
  bool stellar_sources = false;
  bool agn_sources = false;
  RtSourceShapeConfig source_shape = RtSourceShapeConfig::kSedTable;
  std::string stellar_sed_table_path;   // schema chui_ssp_sed_table_v1
  std::string stellar_sed_table_sha256;
  std::string agn_sed_table_path;       // schema chui_agn_sed_table_v1
  std::string agn_sed_table_sha256;
  double source_blackbody_temperature_k = 1.0e5;
  double source_power_law_index = 1.5;  // L_nu ~ nu^-index
  RtUvbPolicyConfig uvb_policy = RtUvbPolicyConfig::kRejectOverlap;
  bool operator==(const RadiationTransportConfig&) const = default;
};

// ------------------------------------------------------------- black holes
enum class BhEngineModelConfig : std::uint8_t { kLegacyHostCell = 0, kReferenceV1 = 1 };
enum class BhSupplyModelConfig : std::uint8_t {
  kBondiTorque = 0,             // hot Bondi + cold gravitational torque (partition of unity)
  kBondi = 1,                   // unboosted / alpha-boosted Bondi-Hoyle-Lyttleton
  kBondiAngularSuppressed = 2,  // Rosas-Guevara et al. (2015) comparator
  kTorque = 3,                  // Hopkins & Quataert (2011) torque only
  kResolvedFlux = 4,            // P1: implemented, rejected at admission
};
enum class BhDiskModelConfig : std::uint8_t { kInstantaneous = 0, kScalarReservoir = 1 };
enum class BhMassLedgerConfig : std::uint8_t { kLegacyNewtonian = 0, kRelativisticDefect = 1 };
enum class BhEddingtonCapConfig : std::uint8_t { kNone = 0, kReferenceRate = 1 };
enum class BhDynamicsModelConfig : std::uint8_t { kNone = 0, kUnresolvedChandrasekhar = 1 };
enum class BhMergerModelConfig : std::uint8_t { kNone = 0, kBoundPair = 1 };
enum class BhSeedModelConfig : std::uint8_t { kNone = 0, kHaloFofSerial = 1 };
enum class BhSpinModelConfig : std::uint8_t { kNone = 0, kKerrThinDisk = 1 };  // P1: rejected
enum class BhAccretionStateModelConfig : std::uint8_t { kNone = 0, kThickThinSlim = 1 };  // P1: rejected

struct BlackHoleEngineConfig {
  BhEngineModelConfig model = BhEngineModelConfig::kLegacyHostCell;
  BhSupplyModelConfig supply_model = BhSupplyModelConfig::kBondiTorque;
  double angular_suppression_c_visc = 0.0;       // Rosas-Guevara C_visc (required > 0 when selected)
  double phase_temperature_split_k = 1.0e5;      // hot/cold partition (SIMBA reference choice)
  double torque_epsilon = 0.0;                   // epsilon_T supply normalization (no default calibration)
  double kernel_radius_min_kpc_phys = 0.0;       // physical, > 0 required
  double kernel_radius_max_kpc_phys = 0.0;       // physical, >= min
  double kernel_target_gas_mass_code = 0.0;      // > 0 required
  std::uint32_t kernel_min_cells = 8U;           // adequacy flag threshold
  std::uint32_t torque_min_samples = 32U;        // kinematic decomposition adequacy
  BhDiskModelConfig disk_model = BhDiskModelConfig::kScalarReservoir;
  double disk_viscous_time_years = 0.0;          // required > 0 for scalar_reservoir
  double disk_wind_loading = 0.0;                // eta_w >= 0
  double disk_max_relative_bh_mass_change = 0.01;  // controlled substep tolerance
  BhMassLedgerConfig mass_ledger = BhMassLedgerConfig::kLegacyNewtonian;
  BhEddingtonCapConfig eddington_cap = BhEddingtonCapConfig::kReferenceRate;
  double eddington_reference_efficiency = 0.1;   // mdot_Edd,ref = L_Edd/(eps_ref c^2)
  BhDynamicsModelConfig dynamics_model = BhDynamicsModelConfig::kNone;
  std::uint32_t df_neighbor_count = 64U;
  double df_timestep_eta = 0.1;
  BhMergerModelConfig merger_model = BhMergerModelConfig::kNone;
  double merger_radius_softening_factor = 2.0;
  BhSeedModelConfig seed_model = BhSeedModelConfig::kNone;
  std::uint32_t seed_cadence_steps = 1U;
  BhSpinModelConfig spin_model = BhSpinModelConfig::kNone;
  BhAccretionStateModelConfig accretion_state_model = BhAccretionStateModelConfig::kNone;
  bool operator==(const BlackHoleEngineConfig&) const = default;
};

// -------------------------------------------------------------------- AGN
enum class AgnFeedbackModelConfig : std::uint8_t { kLegacyThermalHost = 0, kChannelsV1 = 1 };
enum class AgnThermalModeConfig : std::uint8_t { kContinuous = 0, kReservoir = 1 };
enum class AgnRadiationCouplingConfig : std::uint8_t {
  kNone = 0,             // radiative budget recorded as escaping (no coupling)
  kRtTransport = 1,      // radiative budget transported by RT (requires rt agn_sources)
  kLocalEmpirical = 2,   // named no-RT local coupling of a declared fraction
};

struct AgnFeedbackConfig {
  AgnFeedbackModelConfig model = AgnFeedbackModelConfig::kLegacyThermalHost;
  AgnThermalModeConfig thermal_mode = AgnThermalModeConfig::kContinuous;
  double thermal_efficiency = 0.0;     // fraction of dm0 c^2
  double thermal_delta_t_k = 0.0;      // reservoir target heating (required > 0 for reservoir)
  double kick_efficiency = 0.0;        // momentum-only kicks: fraction of dm0 c^2
  double wind_efficiency = 0.0;        // mass-loaded wide-angle wind
  double wind_velocity_km_s = 0.0;
  double jet_efficiency = 0.0;         // bipolar volumetric jet
  double jet_velocity_km_s = 0.0;
  double jet_half_opening_deg = 15.0;
  double max_velocity_c_fraction = 0.1;  // declared nonrelativistic domain
  AgnRadiationCouplingConfig radiation_coupling = AgnRadiationCouplingConfig::kNone;
  double local_empirical_fraction = 0.0;
  bool refinement_enabled = false;
  double refinement_injection_cells = 8.0;  // N_j: injection diameter / dx_phys
  bool magnetic_injection = false;  // P1: rejected at admission
  bool nozzle = false;              // P1: rejected at admission
  bool operator==(const AgnFeedbackConfig&) const = default;
};

[[nodiscard]] std::string rtTransportConfigToString(RtTransportConfig value);
[[nodiscard]] std::string rtIntegratorConfigToString(RtIntegratorConfig value);
[[nodiscard]] std::string rtSpectralProfileConfigToString(RtSpectralProfileConfig value);
[[nodiscard]] std::string rtRiemannSolverConfigToString(RtRiemannSolverConfig value);
[[nodiscard]] std::string rtReconstructionConfigToString(RtReconstructionConfig value);
[[nodiscard]] std::string rtOpacityModelConfigToString(RtOpacityModelConfig value);
[[nodiscard]] std::string rtBoundaryConfigToString(RtBoundaryConfig value);
[[nodiscard]] std::string rtSourceShapeConfigToString(RtSourceShapeConfig value);
[[nodiscard]] std::string rtUvbPolicyConfigToString(RtUvbPolicyConfig value);
[[nodiscard]] std::string bhEngineModelConfigToString(BhEngineModelConfig value);
[[nodiscard]] std::string bhSupplyModelConfigToString(BhSupplyModelConfig value);
[[nodiscard]] std::string bhDiskModelConfigToString(BhDiskModelConfig value);
[[nodiscard]] std::string bhMassLedgerConfigToString(BhMassLedgerConfig value);
[[nodiscard]] std::string bhEddingtonCapConfigToString(BhEddingtonCapConfig value);
[[nodiscard]] std::string bhDynamicsModelConfigToString(BhDynamicsModelConfig value);
[[nodiscard]] std::string bhMergerModelConfigToString(BhMergerModelConfig value);
[[nodiscard]] std::string bhSeedModelConfigToString(BhSeedModelConfig value);
[[nodiscard]] std::string bhSpinModelConfigToString(BhSpinModelConfig value);
[[nodiscard]] std::string bhAccretionStateModelConfigToString(BhAccretionStateModelConfig value);
[[nodiscard]] std::string agnFeedbackModelConfigToString(AgnFeedbackModelConfig value);
[[nodiscard]] std::string agnThermalModeConfigToString(AgnThermalModeConfig value);
[[nodiscard]] std::string agnRadiationCouplingConfigToString(AgnRadiationCouplingConfig value);

}  // namespace cosmosim::core
