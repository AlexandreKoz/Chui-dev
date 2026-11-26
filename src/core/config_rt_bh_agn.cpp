#include "cosmosim/core/config_rt_bh_agn.hpp"

#include <cmath>
#include <string>

#include "cosmosim/core/config.hpp"
#include "core/internal/config_rt_bh_agn_validation.hpp"
#include "core/internal/config_rt_bh_agn_tables.hpp"

namespace cosmosim::core {

std::string rtTransportConfigToString(RtTransportConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_rt_transport_config_names, "RtTransportConfig");
}
std::string rtIntegratorConfigToString(RtIntegratorConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_rt_integrator_config_names, "RtIntegratorConfig");
}
std::string rtSpectralProfileConfigToString(RtSpectralProfileConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_rt_spectral_profile_config_names, "RtSpectralProfileConfig");
}
std::string rtRiemannSolverConfigToString(RtRiemannSolverConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_rt_riemann_solver_config_names, "RtRiemannSolverConfig");
}
std::string rtReconstructionConfigToString(RtReconstructionConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_rt_reconstruction_config_names, "RtReconstructionConfig");
}
std::string rtOpacityModelConfigToString(RtOpacityModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_rt_opacity_model_config_names, "RtOpacityModelConfig");
}
std::string rtBoundaryConfigToString(RtBoundaryConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_rt_boundary_config_names, "RtBoundaryConfig");
}
std::string rtSourceShapeConfigToString(RtSourceShapeConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_rt_source_shape_config_names, "RtSourceShapeConfig");
}
std::string rtUvbPolicyConfigToString(RtUvbPolicyConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_rt_uvb_policy_config_names, "RtUvbPolicyConfig");
}
std::string bhEngineModelConfigToString(BhEngineModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_engine_model_config_names, "BhEngineModelConfig");
}
std::string bhSupplyModelConfigToString(BhSupplyModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_supply_model_config_names, "BhSupplyModelConfig");
}
std::string bhDiskModelConfigToString(BhDiskModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_disk_model_config_names, "BhDiskModelConfig");
}
std::string bhMassLedgerConfigToString(BhMassLedgerConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_mass_ledger_config_names, "BhMassLedgerConfig");
}
std::string bhEddingtonCapConfigToString(BhEddingtonCapConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_eddington_cap_config_names, "BhEddingtonCapConfig");
}
std::string bhDynamicsModelConfigToString(BhDynamicsModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_dynamics_model_config_names, "BhDynamicsModelConfig");
}
std::string bhMergerModelConfigToString(BhMergerModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_merger_model_config_names, "BhMergerModelConfig");
}
std::string bhSeedModelConfigToString(BhSeedModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_seed_model_config_names, "BhSeedModelConfig");
}
std::string bhSpinModelConfigToString(BhSpinModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_spin_model_config_names, "BhSpinModelConfig");
}
std::string bhAccretionStateModelConfigToString(BhAccretionStateModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_bh_accretion_state_model_config_names,
                                         "BhAccretionStateModelConfig");
}
std::string agnFeedbackModelConfigToString(AgnFeedbackModelConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_agn_feedback_model_config_names, "AgnFeedbackModelConfig");
}
std::string agnThermalModeConfigToString(AgnThermalModeConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_agn_thermal_mode_config_names, "AgnThermalModeConfig");
}
std::string agnRadiationCouplingConfigToString(AgnRadiationCouplingConfig v) {
  return internal::rtBhAgnChoiceToString(v, internal::k_agn_radiation_coupling_config_names,
                                         "AgnRadiationCouplingConfig");
}

namespace internal {
namespace {

void require(bool condition, const std::string& message) {
  if (!condition) throw ConfigError(message);
}
[[nodiscard]] bool finitePositive(double v) { return std::isfinite(v) && v > 0.0; }
[[nodiscard]] bool finiteNonNegative(double v) { return std::isfinite(v) && v >= 0.0; }
[[nodiscard]] bool isCosmological(const SimulationConfig& c) {
  return c.mode.mode == SimulationMode::kCosmoCube || c.mode.mode == SimulationMode::kZoomIn;
}
[[nodiscard]] bool networkThermochemistry(const PhysicsConfig& p) {
  return p.thermochemistry_backend == ThermochemistryBackendMode::kFixedHybrid ||
      p.thermochemistry_backend == ThermochemistryBackendMode::kFixedDetailed;
}

void validateRadiation(const SimulationConfig& config) {
  const RadiationTransportConfig& rt = config.physics.rt;
  if (!rt.enabled) {
    require(rt == RadiationTransportConfig{},
            "physics.rt_enabled=false but one or more physics.rt_* keys differ from their defaults; those keys "
            "would have no effect. Remove them or set physics.rt_enabled=true");
    return;
  }
  require(rt.transport == RtTransportConfig::kM1,
          "physics.rt_transport=" + rtTransportConfigToString(rt.transport) +
              " is a P1 backend: its source implementation exists but it is NOT qualified; production admission "
              "is refused in pass 01 (use m1)");
  require(rt.spectral_profile != RtSpectralProfileConfig::kAgnHard,
          "physics.rt_spectral_profile=agn_hard requires qualified Compton/secondary-electron/X-ray microphysics "
          "that does not exist in the thermochemistry authority (FAIL_CLOSED)");
  switch (rt.spectral_profile) {
    case RtSpectralProfileConfig::kGrey:
      require(rt.spectrum_upper_edge_ev == 0.0,
              "physics.rt_spectral_profile=grey is a bolometric thermal group; physics.rt_spectrum_upper_edge_ev "
              "must stay 0");
      require(rt.opacity_model == RtOpacityModelConfig::kConstantReference,
              "physics.rt_spectral_profile=grey requires an explicit opacity model; only constant_reference "
              "(REFERENCE ONLY) exists, no production grey opacity table is provided");
      break;
    case RtSpectralProfileConfig::kIonizing3:
      require(std::isfinite(rt.spectrum_upper_edge_ev) && rt.spectrum_upper_edge_ev > 54.4,
              "physics.rt_spectral_profile=ionizing3 requires physics.rt_spectrum_upper_edge_ev > 54.4 eV "
              "(explicit upper edge; no implicit spectral truncation)");
      break;
    case RtSpectralProfileConfig::kBroadUvIr:
      require(rt.spectrum_upper_edge_ev == 0.0,
              "physics.rt_spectral_profile=broad_uv_ir has fixed edges 0.01..100 eV; "
              "physics.rt_spectrum_upper_edge_ev must stay 0");
      require(false,
              "physics.rt_spectral_profile=broad_uv_ir: the non-ionizing groups need dust/IR opacities from the dust "
              "microphysics authority; no qualified provider exists, so those groups would be untreated "
              "(FAIL_CLOSED)");
      break;
    case RtSpectralProfileConfig::kAgnHard:
      break;
  }
  require(std::isfinite(rt.c_fraction) && rt.c_fraction > 0.0 && rt.c_fraction <= 1.0,
          "physics.rt_c_fraction must lie in (0, 1]");
  require(!(isCosmological(config) && rt.c_fraction < 1.0),
          "physics.rt_c_fraction < 1 in a cosmological mode: the cosmological reduced-speed-of-light equation "
          "and its invariants are not derived/qualified (RT companion section 4); use full c");
  require(std::isfinite(rt.cfl) && rt.cfl > 0.0 && rt.cfl <= 0.5,
          "physics.rt_cfl must lie in (0, 0.5] (sum-of-speeds multidimensional CFL with MC reconstruction)");
  require(rt.max_subcycles >= 1U, "physics.rt_max_subcycles must be >= 1");
  require(std::isfinite(rt.newton_relative_tolerance) && rt.newton_relative_tolerance > 0.0 &&
              rt.newton_relative_tolerance <= 1.0e-4,
          "physics.rt_newton_relative_tolerance must lie in (0, 1e-4]");
  require(rt.newton_max_iterations >= 1U && rt.newton_max_iterations <= 1000U,
          "physics.rt_newton_max_iterations must lie in [1, 1000]");
  require(rt.max_retries <= 8U, "physics.rt_max_retries must lie in [0, 8]");
  if (rt.opacity_model == RtOpacityModelConfig::kConstantReference) {
    require(!config.physics.enable_cooling,
            "physics.rt_opacity_model=constant_reference couples gas and radiation through LTE emission/absorption "
            "(the radiative cooling of the reference RHD problem); physics.enable_cooling must be false to avoid "
            "duplicate cooling");
    require(rt.allow_reference_models,
            "physics.rt_opacity_model=constant_reference is REFERENCE ONLY; set physics.rt_allow_reference_models="
            "true to admit it for idealized tests");
    require(finiteNonNegative(rt.constant_absorption_cm2_g) && finiteNonNegative(rt.constant_scattering_cm2_g) &&
                (rt.constant_absorption_cm2_g + rt.constant_scattering_cm2_g) > 0.0,
            "physics.rt_constant_absorption_cm2_g / rt_constant_scattering_cm2_g must be >= 0 with a positive sum");
  } else {
    require(rt.constant_absorption_cm2_g == 0.0 && rt.constant_scattering_cm2_g == 0.0,
            "physics.rt_constant_* opacities are only meaningful with physics.rt_opacity_model=constant_reference");
    require(false,
            "physics.rt_opacity_model=photoionization: the RT-side frozen-species absorption, transported photo-rate "
            "provider and chemistry photoionization counters exist (amendment A5), but the per-row reconciliation "
            "of RT-absorbed photon energy against the chemistry-integrated photoionizations is not yet wired into "
            "BaryonicSourceStage; without it one absorption event would not create identical matter/radiation "
            "debit and credit (FAIL_CLOSED, pass 01)");
    require(networkThermochemistry(config.physics),
            "physics.rt_opacity_model=photoionization requires a species-evolving thermochemistry backend "
            "(fixed_hybrid or fixed_detailed): opacities are sigma_g n_species of the chemistry authority");
    require(config.physics.enable_cooling,
            "physics.rt_opacity_model=photoionization requires physics.enable_cooling=true (single thermochemical "
            "authority couples to the radiation field)");
  }
  if (rt.source_shape != RtSourceShapeConfig::kSedTable) {
    require(rt.allow_reference_models,
            "physics.rt_source_shape=" + rtSourceShapeConfigToString(rt.source_shape) +
                " is an analytic REFERENCE ONLY spectrum (not an astrophysical SED); set "
                "physics.rt_allow_reference_models=true");
    require(finitePositive(rt.source_blackbody_temperature_k) && std::isfinite(rt.source_power_law_index),
            "physics.rt_source_blackbody_temperature_k must be > 0 and physics.rt_source_power_law_index finite");
  }
  if (rt.stellar_sources) {
    require(false,
            "physics.rt_stellar_sources: no age/metallicity-dependent stellar ionizing-luminosity provider exists in "
            "the stellar-population authority (the SED-table schema carries a spectral shape only); FAIL_CLOSED");
    require(config.physics.enable_star_formation, "physics.rt_stellar_sources requires star formation");
    if (rt.source_shape == RtSourceShapeConfig::kSedTable) {
      require(!rt.stellar_sed_table_path.empty() && rt.stellar_sed_table_sha256.size() == 64U,
              "physics.rt_stellar_sources with rt_source_shape=sed_table requires physics.rt_stellar_sed_table_path "
              "and a 64-hex physics.rt_stellar_sed_table_sha256 (schema chui_ssp_sed_table_v1); no stellar SED is "
              "bundled");
    }
  }
  if (rt.agn_sources) {
    require(config.physics.agn.radiation_coupling == AgnRadiationCouplingConfig::kRtTransport,
            "physics.rt_agn_sources=true requires physics.agn_radiation_coupling=rt_transport (the AGN radiative "
            "budget must be transported exactly once)");
    if (rt.source_shape == RtSourceShapeConfig::kSedTable) {
      require(!rt.agn_sed_table_path.empty() && rt.agn_sed_table_sha256.size() == 64U,
              "physics.rt_agn_sources with rt_source_shape=sed_table requires physics.rt_agn_sed_table_path and a "
              "64-hex physics.rt_agn_sed_table_sha256 (schema chui_agn_sed_table_v1); no AGN SED is bundled");
    }
  }
  if (rt.stellar_sources || rt.agn_sources || rt.opacity_model == RtOpacityModelConfig::kPhotoionization) {
    if (config.physics.uv_background_model != UvBackgroundModel::kNone) {
      require(rt.uvb_policy == RtUvbPolicyConfig::kBackgroundPlusLocalSources,
              "a UV background is configured together with transported ionizing radiation; select "
              "physics.rt_uvb_policy=background_plus_local_sources (UVB stays in chemistry, RT transports local "
              "sources only, disjoint by origin) or disable the UVB");
    }
  }
  require(config.parallel.mpi_ranks_expected <= 1,
          "physics.rt_enabled=true with more than one MPI rank: the RT ghost/halo exchange across ranks is not "
          "implemented in pass 01 (single-rank level-synchronous transport only; FAIL_CLOSED)");
  require(config.numerics.amr_regrid_mode == AmrRegridConfigMode::kOff,
          "physics.rt_enabled=true requires numerics.amr_regrid_mode=off: radiation moments are not yet carried "
          "through regrid/migration transactions (FAIL_CLOSED)");
  require(config.numerics.amr_time_stepping == AmrTimeSteppingConfigMode::kSynchronized,
          "physics.rt_enabled=true requires numerics.amr_time_stepping=synchronized (level-synchronous RT; no "
          "falsely supported asynchronous mode)");
  require(!config.mhd.enabled,
          "physics.rt_enabled=true with mhd.enabled=true: radiation-MHD coupling (radiation pressure on "
          "magnetized gas, Lorentz-Heaviside energy reconciliation) is outside the MHD reference scope and not "
          "qualified in pass 01 (FAIL_CLOSED)");
}

void validateBlackHoleEngine(const SimulationConfig& config) {
  const PhysicsConfig& p = config.physics;
  const BlackHoleEngineConfig& bh = p.bh_engine;
  if (bh.model == BhEngineModelConfig::kLegacyHostCell) {
    require(bh == BlackHoleEngineConfig{},
            "physics.bh_engine_model=legacy_host_cell but one or more advanced physics.bh_* engine keys differ "
            "from their defaults; those keys would have no effect. Remove them or select reference_v1");
    return;
  }
  require(p.enable_black_hole_agn, "physics.bh_engine_model=reference_v1 requires physics.enable_black_hole_agn=true");
  require(bh.supply_model != BhSupplyModelConfig::kResolvedFlux,
          "physics.bh_supply_model=resolved_flux is a P1 closure: implemented, NOT qualified; admission refused");
  require(finitePositive(bh.kernel_radius_min_kpc_phys) && finitePositive(bh.kernel_radius_max_kpc_phys) &&
              bh.kernel_radius_max_kpc_phys >= bh.kernel_radius_min_kpc_phys,
          "physics.bh_kernel_radius_{min,max}_kpc_phys must be > 0 with max >= min (explicit physical support)");
  require(finitePositive(bh.kernel_target_gas_mass_code), "physics.bh_kernel_target_gas_mass_code must be > 0");
  require(bh.kernel_min_cells >= 1U, "physics.bh_kernel_min_cells must be >= 1");
  require(finitePositive(bh.phase_temperature_split_k), "physics.bh_phase_temperature_split_k must be > 0");
  const bool uses_torque = bh.supply_model == BhSupplyModelConfig::kBondiTorque ||
      bh.supply_model == BhSupplyModelConfig::kTorque;
  if (uses_torque) {
    require(finitePositive(bh.torque_epsilon),
            "the torque supply requires an explicit physics.bh_torque_epsilon > 0: SIMBA's normalization includes "
            "unresolved retention and must be re-calibrated as a disk-boundary SUPPLY (BH companion 4.2); no "
            "default calibration is provided");
    require(bh.torque_min_samples >= 1U, "physics.bh_torque_min_samples must be >= 1");
  } else {
    require(bh.torque_epsilon == 0.0, "physics.bh_torque_epsilon is only meaningful with a torque supply model");
  }
  if (bh.supply_model == BhSupplyModelConfig::kBondiAngularSuppressed) {
    require(finitePositive(bh.angular_suppression_c_visc),
            "physics.bh_supply_model=bondi_angular_suppressed requires physics.bh_angular_suppression_c_visc > 0 "
            "(Rosas-Guevara et al. 2015 viscosity parameter)");
  } else {
    require(bh.angular_suppression_c_visc == 0.0,
            "physics.bh_angular_suppression_c_visc is only meaningful with bondi_angular_suppressed");
  }
  if (bh.disk_model == BhDiskModelConfig::kScalarReservoir) {
    require(finitePositive(bh.disk_viscous_time_years),
            "physics.bh_disk_model=scalar_reservoir requires physics.bh_disk_viscous_time_years > 0");
  } else {
    require(bh.disk_viscous_time_years == 0.0,
            "physics.bh_disk_viscous_time_years is only meaningful with the scalar reservoir");
    require(bh.disk_wind_loading == 0.0,
            "physics.bh_disk_wind_loading requires the scalar reservoir (the instantaneous comparator has no disk)");
  }
  require(finiteNonNegative(bh.disk_wind_loading), "physics.bh_disk_wind_loading must be >= 0");
  require(finitePositive(bh.disk_max_relative_bh_mass_change) && bh.disk_max_relative_bh_mass_change <= 0.5,
          "physics.bh_disk_max_relative_bh_mass_change must lie in (0, 0.5]");
  require(finitePositive(bh.eddington_reference_efficiency) && bh.eddington_reference_efficiency < 1.0,
          "physics.bh_eddington_reference_efficiency must lie in (0, 1)");
  require(p.bh_use_eddington_cap == (bh.eddington_cap == BhEddingtonCapConfig::kReferenceRate),
          "compatibility mapping: physics.bh_use_eddington_cap must agree with physics.bh_eddington_cap "
          "(true <-> reference_rate)");
  if (bh.dynamics_model == BhDynamicsModelConfig::kUnresolvedChandrasekhar) {
    require(bh.df_neighbor_count >= 8U, "physics.bh_df_neighbor_count must be >= 8");
    require(finitePositive(bh.df_timestep_eta) && bh.df_timestep_eta <= 1.0,
            "physics.bh_df_timestep_eta must lie in (0, 1]");
  }
  if (bh.merger_model == BhMergerModelConfig::kBoundPair) {
    require(finitePositive(bh.merger_radius_softening_factor),
            "physics.bh_merger_radius_softening_factor must be > 0");
  }
  if (bh.seed_model == BhSeedModelConfig::kHaloFofSerial) {
    require(config.parallel.mpi_ranks_expected <= 1,
            "physics.bh_seed_model=halo_fof_serial is a single-rank provider; distributed seeding has no bounded "
            "authoritative halo provider (root-global FoF is inadmissible) and stays FAIL_CLOSED");
    require(bh.seed_cadence_steps >= 1U, "physics.bh_seed_cadence_steps must be >= 1");
  }
  require(p.bh_enable_seeding == (bh.seed_model != BhSeedModelConfig::kNone),
          "compatibility mapping: with the reference engine physics.bh_enable_seeding must agree with "
          "physics.bh_seed_model (true <-> a seed provider is selected)");
  require(bh.spin_model == BhSpinModelConfig::kNone,
          "physics.bh_spin_model=kerr_thin_disk is a P1 model: Kerr ISCO/vector-J source exists but published warp "
          "alignment coefficients are not pinned and nothing is qualified; admission refused");
  require(bh.accretion_state_model == BhAccretionStateModelConfig::kNone,
          "physics.bh_accretion_state_model=thick_thin_slim requires a versioned published efficiency/photon-"
          "trapping model that is not pinned (no invented efficiency fits); admission refused");
}

void validateAgn(const SimulationConfig& config) {
  const PhysicsConfig& p = config.physics;
  const AgnFeedbackConfig& agn = p.agn;
  if (agn.model == AgnFeedbackModelConfig::kLegacyThermalHost) {
    require(agn == AgnFeedbackConfig{},
            "physics.agn_feedback_model=legacy_thermal_host but one or more physics.agn_* keys differ from their "
            "defaults; those keys would have no effect. Remove them or select channels_v1");
    require(p.bh_engine.model == BhEngineModelConfig::kLegacyHostCell,
            "physics.bh_engine_model=reference_v1 requires physics.agn_feedback_model=channels_v1 (the legacy "
            "thermal host deposit belongs to the legacy callback only)");
    return;
  }
  require(p.bh_engine.model == BhEngineModelConfig::kReferenceV1,
          "physics.agn_feedback_model=channels_v1 requires physics.bh_engine_model=reference_v1 (budgets come from "
          "the BH/disk engine)");
  const double sum = agn.thermal_efficiency + agn.kick_efficiency + agn.wind_efficiency + agn.jet_efficiency;
  require(finiteNonNegative(agn.thermal_efficiency) && finiteNonNegative(agn.kick_efficiency) &&
              finiteNonNegative(agn.wind_efficiency) && finiteNonNegative(agn.jet_efficiency),
          "physics.agn_*_efficiency values must be >= 0");
  require(sum + p.bh_epsilon_r < 1.0,
          "sum of AGN mechanical efficiencies plus physics.bh_epsilon_r must be < 1 (fraction of dm0 c^2)");
  if (agn.thermal_mode == AgnThermalModeConfig::kReservoir) {
    require(finitePositive(agn.thermal_delta_t_k),
            "physics.agn_thermal_mode=reservoir requires physics.agn_thermal_delta_t_k > 0");
  } else {
    require(agn.thermal_delta_t_k == 0.0, "physics.agn_thermal_delta_t_k is only meaningful in reservoir mode");
  }
  require(finitePositive(agn.max_velocity_c_fraction) && agn.max_velocity_c_fraction <= 0.3,
          "physics.agn_max_velocity_c_fraction must lie in (0, 0.3] (declared nonrelativistic domain)");
  const double c_km_s = 299792.458;
  if (agn.wind_efficiency > 0.0) {
    require(p.bh_engine.disk_wind_loading > 0.0,
            "physics.agn_wind_efficiency > 0 requires a mass origin: physics.bh_disk_wind_loading > 0 (winds are "
            "loaded with disk-wind mass; entrainment of ambient gas is the kick/jet channel, never new mass)");
    require(finitePositive(agn.wind_velocity_km_s), "a wind channel requires physics.agn_wind_velocity_km_s > 0");
    require(agn.wind_velocity_km_s <= agn.max_velocity_c_fraction * c_km_s,
            "physics.agn_wind_velocity_km_s exceeds the declared nonrelativistic domain (rejected, never capped)");
  }
  if (agn.jet_efficiency > 0.0) {
    require(finitePositive(agn.jet_velocity_km_s), "a jet channel requires physics.agn_jet_velocity_km_s > 0");
    require(agn.jet_velocity_km_s <= agn.max_velocity_c_fraction * c_km_s,
            "physics.agn_jet_velocity_km_s exceeds the declared nonrelativistic domain (rejected, never capped)");
    require(agn.jet_half_opening_deg > 0.0 && agn.jet_half_opening_deg < 90.0,
            "physics.agn_jet_half_opening_deg must lie in (0, 90)");
  }
  switch (agn.radiation_coupling) {
    case AgnRadiationCouplingConfig::kNone:
      require(agn.local_empirical_fraction == 0.0, "physics.agn_local_empirical_fraction requires local_empirical");
      break;
    case AgnRadiationCouplingConfig::kRtTransport:
      require(p.rt.enabled && p.rt.agn_sources,
              "physics.agn_radiation_coupling=rt_transport requires physics.rt_enabled and physics.rt_agn_sources");
      require(agn.local_empirical_fraction == 0.0,
              "the no-RT local empirical coupling is mutually exclusive with transporting the same photons");
      break;
    case AgnRadiationCouplingConfig::kLocalEmpirical:
      require(!(p.rt.enabled && p.rt.agn_sources),
              "the no-RT local empirical coupling is mutually exclusive with RT transport of AGN photons");
      require(std::isfinite(agn.local_empirical_fraction) && agn.local_empirical_fraction > 0.0 &&
                  agn.local_empirical_fraction <= 1.0,
              "physics.agn_local_empirical_fraction must lie in (0, 1]");
      break;
  }
  if (agn.refinement_enabled) {
    require(finitePositive(agn.refinement_injection_cells),
            "physics.agn_refinement_injection_cells must be > 0 (chosen by convergence, not a physical constant)");
  }
  require(!agn.magnetic_injection,
          "physics.agn_magnetic_injection is a P1 CT-compatible channel: source exists but is not qualified; "
          "admission refused");
  require(!agn.nozzle, "physics.agn_nozzle is a P1 resolved-nozzle channel: not qualified; admission refused");
  // MHD runs: AGN gas transactions act on gas lanes only and never touch the
  // CT face flux (existing MHD source-owner policy, mhd pass 02 B2); the BH
  // sink keeps the magnetic field (retained-field sink policy, recorded in
  // the BH engine ledger as retained magnetic energy, not removed).
}

}  // namespace

void validateRtBhAgnConfig(const SimulationConfig& config) {
  validateRadiation(config);
  validateBlackHoleEngine(config);
  validateAgn(config);
}

}  // namespace internal
}  // namespace cosmosim::core
