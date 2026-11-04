#include "cosmosim/physics/baryonic_runtime_config.hpp"

#include <stdexcept>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

ThermochemistryBackend thermochemistryBackendFromConfig(const core::PhysicsConfig& physics) {
  switch (physics.thermochemistry_backend) {
    case core::ThermochemistryBackendMode::kLegacyCooling: return ThermochemistryBackend::kLegacyCooling;
    case core::ThermochemistryBackendMode::kAnalyticReference: return ThermochemistryBackend::kAnalyticReference;
    case core::ThermochemistryBackendMode::kEquilibriumPrimordial: return ThermochemistryBackend::kEquilibriumPrimordial;
    case core::ThermochemistryBackendMode::kFixedHybrid: return ThermochemistryBackend::kFixedHybrid;
    case core::ThermochemistryBackendMode::kFixedDetailed: return ThermochemistryBackend::kFixedDetailed;
  }
  throw std::invalid_argument("unhandled thermochemistry backend");
}

bool part2ThermochemistrySelected(const core::PhysicsConfig& physics) {
  return physics.thermochemistry_backend != core::ThermochemistryBackendMode::kLegacyCooling;
}

CompositionLayoutRequest compositionLayoutRequestFromConfig(const core::PhysicsConfig& physics) {
  CompositionLayoutRequest request;
  switch (physics.thermochemistry_backend) {
    case core::ThermochemistryBackendMode::kFixedHybrid:
      request.network = CompositionNetworkKind::kHybridLight;
      break;
    case core::ThermochemistryBackendMode::kFixedDetailed:
      request.network = CompositionNetworkKind::kDetailedCarbonOxygen;
      break;
    default:
      request.network = CompositionNetworkKind::kElementsOnly;
      break;
  }
  request.dust = physics.dust_model == core::DustModelConfig::kTwoReservoirV1
      ? CompositionDustKind::kTwoReservoirV1
      : CompositionDustKind::kNone;
  request.sgs_energy = physics.sgs_turbulence_model == core::SgsTurbulenceModelConfig::kBudgetedV1;
  return request;
}

PrimordialCompositionPolicy primordialPolicyFromConfig(const core::PhysicsConfig& physics) {
  return PrimordialCompositionPolicy{.helium_mass_fraction = physics.chemistry_primordial_helium_mass_fraction};
}

bool ensureGasCompositionInitialized(core::SimulationState& state, const core::PhysicsConfig& physics) {
  if (!part2ThermochemistrySelected(physics)) {
    if (state.gas_composition.enabled()) {
      throw std::runtime_error(
          "state carries a CHUI gas composition sidecar but the configuration selects legacy_cooling; "
          "resuming element-resolved state with element tracking disabled would discard it");
    }
    return false;
  }
  const CompositionLayout layout = CompositionLayout::make(compositionLayoutRequestFromConfig(physics));
  const bool was_enabled = state.gas_composition.enabled();
  if (!was_enabled && restoreGasCompositionFromModuleSidecar(state, layout)) {
    return false;  // restored from a checkpoint, not initialized
  }
  enableGasComposition(state, layout, primordialPolicyFromConfig(physics));
  return !was_enabled;
}

BaryonicDataProducts loadBaryonicDataProducts(const core::PhysicsConfig& physics) {
  BaryonicDataProducts data;
  if (physics.uv_background_model == core::UvBackgroundModel::kTreecoolTable) {
    data.uv_background = std::make_shared<const UvBackgroundTable>(
        UvBackgroundTable::loadTreecool(physics.uv_background_table_path, physics.uv_background_table_sha256));
  }
  if (physics.metal_cooling_mode == core::MetalCoolingModeConfig::kElementTable ||
      physics.metal_cooling_mode == core::MetalCoolingModeConfig::kElementTableAndFineStructure) {
    data.metal_cooling_table = std::make_shared<const MetalCoolingTable>(
        MetalCoolingTable::loadFromFile(physics.metal_cooling_table_path, physics.metal_cooling_table_sha256));
    // Throws on a refused UVB mismatch (Part 2.5 policy).
    data.uv_background_compatibility = data.metal_cooling_table->checkUvBackgroundCompatibility(
        data.uv_background != nullptr ? data.uv_background->provenance().sha256_hex : std::string(),
        physics.metal_cooling_accept_uv_background_mismatch);
  }
  return data;
}

ThermochemistryOperatorConfig makeThermochemistryOperatorConfig(
    const core::SimulationConfig& config, const core::UnitSystem& units, const BaryonicDataProducts& data) {
  const core::PhysicsConfig& p = config.physics;
  ThermochemistryOperatorConfig op;
  ThermochemistryConfig& s = op.solver;
  s.backend = thermochemistryBackendFromConfig(p);
  s.relative_tolerance = p.thermochemistry_relative_tolerance;
  s.abundance_absolute_tolerance = p.thermochemistry_abundance_absolute_tolerance;
  s.max_substeps = p.thermochemistry_max_substeps;
  s.max_newton_iterations = p.thermochemistry_max_newton_iterations;
  s.temperature_floor_k = p.temperature_floor_k;
  s.analytic_relaxation_time_s = p.thermochemistry_analytic_relaxation_time_years * k_seconds_per_year;
  s.analytic_equilibrium_temperature_k = p.thermochemistry_analytic_equilibrium_temperature_k;
  s.uv_background = data.uv_background.get();
  s.self_shielding = p.self_shielding_model != core::SelfShieldingModel::kNone;
  s.isrf_g0_habing = p.isrf_g0_habing;
  s.shielding_length_max_cm = p.chemistry_shielding_length_max_pc * 3.0856775814913673e18;
  switch (p.metal_cooling_mode) {
    case core::MetalCoolingModeConfig::kNone: s.metal_mode = MetalCoolingMode::kNone; break;
    case core::MetalCoolingModeConfig::kFineStructureCold: s.metal_mode = MetalCoolingMode::kFineStructureCold; break;
    case core::MetalCoolingModeConfig::kElementTable: s.metal_mode = MetalCoolingMode::kElementTable; break;
    case core::MetalCoolingModeConfig::kElementTableAndFineStructure:
      s.metal_mode = MetalCoolingMode::kElementTableAndFineStructure;
      break;
  }
  s.metal_table = data.metal_cooling_table.get();
  s.metal_table_range = p.metal_cooling_table_range == core::MetalCoolingTableRangeConfig::kReject
      ? MetalTableRangePolicy::kReject
      : MetalTableRangePolicy::kClampDensityAndRedshift;
  s.fine_structure_switch_temperature_k = p.metal_cooling_switch_temperature_k;
  s.allow_unqualified_rates = p.thermochemistry_allow_unqualified_rates;
  s.cosmological = config.mode.mode == core::SimulationMode::kCosmoCube ||
      config.mode.mode == core::SimulationMode::kZoomIn;
  op.units.velocity_code_to_cm_s = units.velocity_si_per_code * 100.0;
  op.units.length_code_to_cm = units.length_si_per_code * 100.0;
  op.units.mass_code_to_g = units.mass_si_per_code * 1000.0;
  op.units.density_code_to_cgs = op.units.mass_code_to_g / (op.units.length_code_to_cm * op.units.length_code_to_cm *
                                                            op.units.length_code_to_cm);
  op.units.comoving = config.units.coordinate_frame == core::CoordinateFrame::kComoving;
  op.metallicity_scaled_dust_proxy = p.chemistry_metallicity_scaled_dust_proxy;
  op.adaptive.enabled = p.chemistry_adaptive_selection == core::ChemistryAdaptiveSelectionConfig::kLocalV1;
  op.adaptive.cold_temperature_k = p.chemistry_adaptive_cold_temperature_k;
  op.adaptive.dense_n_h_cm3 = p.chemistry_adaptive_dense_nh_cm3;
  op.adaptive.hysteresis_fraction = p.chemistry_adaptive_hysteresis_fraction;
  op.adaptive.minimum_dwell_years = p.chemistry_adaptive_min_dwell_years;
  op.adaptive.cooling_time_ratio = p.chemistry_adaptive_cooling_time_ratio;
  op.adaptive.error_estimate_threshold =
      p.chemistry_adaptive_error_estimate_threshold;
  op.adaptive.residual_fraction_threshold =
      p.chemistry_adaptive_residual_fraction_threshold;
  op.shadow.sample_fraction = p.chemistry_shadow_sample_fraction;
  op.shadow.max_cells = p.chemistry_shadow_max_cells;
  return op;
}

IsmRegimeConfig makeIsmRegimeConfig(const core::PhysicsConfig& p) {
  return IsmRegimeConfig{.effective_threshold_n_h_cm3 = p.sf_effective_n_h_threshold_cgs,
                         .jeans_resolution_cells = p.ism_jeans_resolution_cells,
                         .hysteresis_fraction = p.ism_regime_hysteresis_fraction,
                         .minimum_dwell_years = p.ism_regime_min_dwell_years};
}

DustConfig makeDustConfig(const core::PhysicsConfig& p) {
  return DustConfig{.accretion_timescale_ref_s = p.dust_accretion_timescale_ref_years * k_seconds_per_year,
                    .sputtering_timescale_ref_s = p.dust_sputtering_timescale_ref_years * k_seconds_per_year,
                    .sn_swept_mass_per_1e51_erg_g = p.dust_sn_swept_mass_msun * k_solar_mass_g,
                    .sn_destruction_efficiency = p.dust_sn_destruction_efficiency,
                    .agb_condensation_efficiency = p.dust_agb_condensation_efficiency,
                    .cc_condensation_efficiency = p.dust_cc_condensation_efficiency};
}

SgsConfig makeSgsConfig(const core::PhysicsConfig& p) {
  return SgsConfig{
      .dissipation_coefficient = p.sgs_dissipation_coefficient,
      .eddy_viscosity_coefficient = p.sgs_eddy_viscosity_coefficient};
}

}  // namespace cosmosim::physics
