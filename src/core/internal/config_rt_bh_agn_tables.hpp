#pragma once

// Internal spelling tables for the RT/BH/AGN typed options
// (include/cosmosim/core/config_rt_bh_agn.hpp). One table per enum is the
// single authority for both parsing and normalized serialization, so a
// spelling can never drift between the two directions.

#include <algorithm>
#include <array>
#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>

#include "cosmosim/core/config.hpp"

namespace cosmosim::core::internal {

template <typename E>
struct RtBhAgnChoice {
  std::string_view name;
  E value;
};

template <typename E, std::size_t N>
[[nodiscard]] E parseRtBhAgnChoice(const std::string& raw, const std::string& key,
                                   const std::array<RtBhAgnChoice<E>, N>& table) {
  std::string lower = raw;
  std::transform(lower.begin(), lower.end(), lower.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  for (const auto& entry : table) {
    if (entry.name == lower) return entry.value;
  }
  std::string allowed;
  for (const auto& entry : table) {
    if (!allowed.empty()) allowed += ", ";
    allowed += std::string(entry.name);
  }
  throw ConfigError("key '" + key + "': invalid value '" + raw + "' (allowed: " + allowed + ")");
}

template <typename E, std::size_t N>
[[nodiscard]] std::string rtBhAgnChoiceToString(E value, const std::array<RtBhAgnChoice<E>, N>& table,
                                                const char* type_name) {
  for (const auto& entry : table) {
    if (entry.value == value) return std::string(entry.name);
  }
  throw ConfigError(std::string("unhandled ") + type_name + " enum value during serialization");
}

inline constexpr std::array<RtBhAgnChoice<RtTransportConfig>, 3> k_rt_transport_config_names{{
    {"m1", RtTransportConfig::kM1},
    {"m1_direct_rays", RtTransportConfig::kM1DirectRays},
    {"vet", RtTransportConfig::kVet},
}};
inline constexpr std::array<RtBhAgnChoice<RtIntegratorConfig>, 2> k_rt_integrator_config_names{{
    {"imex_pd_ars", RtIntegratorConfig::kImexPdArs},
    {"explicit_reference", RtIntegratorConfig::kExplicitReference},
}};
inline constexpr std::array<RtBhAgnChoice<RtSpectralProfileConfig>, 4> k_rt_spectral_profile_config_names{{
    {"grey", RtSpectralProfileConfig::kGrey},
    {"ionizing3", RtSpectralProfileConfig::kIonizing3},
    {"broad_uv_ir", RtSpectralProfileConfig::kBroadUvIr},
    {"agn_hard", RtSpectralProfileConfig::kAgnHard},
}};
inline constexpr std::array<RtBhAgnChoice<RtRiemannSolverConfig>, 2> k_rt_riemann_solver_config_names{{
    {"hll", RtRiemannSolverConfig::kHll},
    {"glf", RtRiemannSolverConfig::kGlf},
}};
inline constexpr std::array<RtBhAgnChoice<RtReconstructionConfig>, 2> k_rt_reconstruction_config_names{{
    {"plm_mc", RtReconstructionConfig::kPlmMc},
    {"first_order", RtReconstructionConfig::kFirstOrder},
}};
inline constexpr std::array<RtBhAgnChoice<RtOpacityModelConfig>, 2> k_rt_opacity_model_config_names{{
    {"photoionization", RtOpacityModelConfig::kPhotoionization},
    {"constant_reference", RtOpacityModelConfig::kConstantReference},
}};
inline constexpr std::array<RtBhAgnChoice<RtBoundaryConfig>, 3> k_rt_boundary_config_names{{
    {"periodic", RtBoundaryConfig::kPeriodic},
    {"outflow", RtBoundaryConfig::kOutflow},
    {"reflecting", RtBoundaryConfig::kReflecting},
}};
inline constexpr std::array<RtBhAgnChoice<RtSourceShapeConfig>, 3> k_rt_source_shape_config_names{{
    {"sed_table", RtSourceShapeConfig::kSedTable},
    {"blackbody_reference", RtSourceShapeConfig::kBlackbodyReference},
    {"power_law_reference", RtSourceShapeConfig::kPowerLawReference},
}};
inline constexpr std::array<RtBhAgnChoice<RtUvbPolicyConfig>, 2> k_rt_uvb_policy_config_names{{
    {"reject_overlap", RtUvbPolicyConfig::kRejectOverlap},
    {"background_plus_local_sources", RtUvbPolicyConfig::kBackgroundPlusLocalSources},
}};
inline constexpr std::array<RtBhAgnChoice<BhEngineModelConfig>, 2> k_bh_engine_model_config_names{{
    {"legacy_host_cell", BhEngineModelConfig::kLegacyHostCell},
    {"reference_v1", BhEngineModelConfig::kReferenceV1},
}};
inline constexpr std::array<RtBhAgnChoice<BhSupplyModelConfig>, 5> k_bh_supply_model_config_names{{
    {"bondi_torque", BhSupplyModelConfig::kBondiTorque},
    {"bondi", BhSupplyModelConfig::kBondi},
    {"bondi_angular_suppressed", BhSupplyModelConfig::kBondiAngularSuppressed},
    {"torque", BhSupplyModelConfig::kTorque},
    {"resolved_flux", BhSupplyModelConfig::kResolvedFlux},
}};
inline constexpr std::array<RtBhAgnChoice<BhDiskModelConfig>, 2> k_bh_disk_model_config_names{{
    {"instantaneous", BhDiskModelConfig::kInstantaneous},
    {"scalar_reservoir", BhDiskModelConfig::kScalarReservoir},
}};
inline constexpr std::array<RtBhAgnChoice<BhMassLedgerConfig>, 2> k_bh_mass_ledger_config_names{{
    {"legacy_newtonian", BhMassLedgerConfig::kLegacyNewtonian},
    {"relativistic_defect", BhMassLedgerConfig::kRelativisticDefect},
}};
inline constexpr std::array<RtBhAgnChoice<BhEddingtonCapConfig>, 2> k_bh_eddington_cap_config_names{{
    {"none", BhEddingtonCapConfig::kNone},
    {"reference_rate", BhEddingtonCapConfig::kReferenceRate},
}};
inline constexpr std::array<RtBhAgnChoice<BhDynamicsModelConfig>, 2> k_bh_dynamics_model_config_names{{
    {"none", BhDynamicsModelConfig::kNone},
    {"unresolved_chandrasekhar", BhDynamicsModelConfig::kUnresolvedChandrasekhar},
}};
inline constexpr std::array<RtBhAgnChoice<BhMergerModelConfig>, 2> k_bh_merger_model_config_names{{
    {"none", BhMergerModelConfig::kNone},
    {"bound_pair", BhMergerModelConfig::kBoundPair},
}};
inline constexpr std::array<RtBhAgnChoice<BhSeedModelConfig>, 2> k_bh_seed_model_config_names{{
    {"none", BhSeedModelConfig::kNone},
    {"halo_fof_serial", BhSeedModelConfig::kHaloFofSerial},
}};
inline constexpr std::array<RtBhAgnChoice<BhSpinModelConfig>, 2> k_bh_spin_model_config_names{{
    {"none", BhSpinModelConfig::kNone},
    {"kerr_thin_disk", BhSpinModelConfig::kKerrThinDisk},
}};
inline constexpr std::array<RtBhAgnChoice<BhAccretionStateModelConfig>, 2>
    k_bh_accretion_state_model_config_names{{
        {"none", BhAccretionStateModelConfig::kNone},
        {"thick_thin_slim", BhAccretionStateModelConfig::kThickThinSlim},
    }};
inline constexpr std::array<RtBhAgnChoice<AgnFeedbackModelConfig>, 2> k_agn_feedback_model_config_names{{
    {"legacy_thermal_host", AgnFeedbackModelConfig::kLegacyThermalHost},
    {"channels_v1", AgnFeedbackModelConfig::kChannelsV1},
}};
inline constexpr std::array<RtBhAgnChoice<AgnThermalModeConfig>, 2> k_agn_thermal_mode_config_names{{
    {"continuous", AgnThermalModeConfig::kContinuous},
    {"reservoir", AgnThermalModeConfig::kReservoir},
}};
inline constexpr std::array<RtBhAgnChoice<AgnRadiationCouplingConfig>, 3>
    k_agn_radiation_coupling_config_names{{
        {"none", AgnRadiationCouplingConfig::kNone},
        {"rt_transport", AgnRadiationCouplingConfig::kRtTransport},
        {"local_empirical", AgnRadiationCouplingConfig::kLocalEmpirical},
    }};

}  // namespace cosmosim::core::internal
