#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/physics/cooling_heating.hpp"

namespace cosmosim::physics {

inline constexpr std::uint32_t kEffectiveMultiphaseEosSchemaVersion = 1;

struct EffectiveMultiphaseEosConfig {
  std::string parameter_set = "chui_sh03_tng_like_v1";
  double hydrogen_mass_fraction = 0.76;
  double n_h_threshold_cgs = 0.13;
  double min_baryon_overdensity = 0.0;
  double t_star_at_threshold_code = 1.0;
  double evaporation_factor_at_threshold = 1000.0;
  double evaporation_density_exponent = -0.8;
  double cold_phase_temperature_k = 1000.0;
  double supernova_specific_energy_code = 1.0e5;
  double massive_star_fraction = 0.1;
  double q_eos = 0.3;
  double isothermal_temperature_k = 1.0e4;
  double hot_excess_tolerance = 1.0;
  core::EffectiveIsmEosRelaxation relaxation = core::EffectiveIsmEosRelaxation::kInstantaneous;
  double relaxation_timescale_code = 0.0;
  std::uint32_t table_bins = 256;
  double max_density_ratio = 1.0e6;
  core::EffectiveIsmBirthMassConvention birth_mass_convention =
      core::EffectiveIsmBirthMassConvention::kInitialSspMass;
  core::EffectiveIsmFeedbackCoupling feedback_coupling =
      core::EffectiveIsmFeedbackCoupling::kExternalFeedbackCalibrated;
  double adiabatic_index = 5.0 / 3.0;
  double mean_molecular_weight_cold = 1.22;
  double mean_molecular_weight_ionized = 0.59;
};

struct EffectiveMultiphaseEosEntry {
  double density_phys_code = 0.0;
  double density_ratio = 0.0;
  double cold_mass_fraction = 0.0;
  double specific_internal_energy_hot_code = 0.0;
  double specific_internal_energy_eff_code = 0.0;
  double pressure_phys_code = 0.0;
  double signal_speed_squared_code = 0.0;
  double star_formation_timescale_code = 0.0;
  double cooling_timescale_code = 0.0;
  bool valid = false;
};

struct EffectiveMultiphaseEosLookup {
  bool above_threshold = false;
  bool valid = false;
  EffectiveMultiphaseEosEntry entry;
};

class EffectiveMultiphaseEosTable {
 public:
  EffectiveMultiphaseEosTable(
      EffectiveMultiphaseEosConfig config,
      core::UnitSystem units,
      CoolingRateProvider cooling_provider);

  [[nodiscard]] const EffectiveMultiphaseEosConfig& config() const noexcept;
  [[nodiscard]] const core::UnitSystem& units() const noexcept;
  [[nodiscard]] double thresholdDensityPhysCode() const noexcept;
  [[nodiscard]] double thresholdDensityCgs() const noexcept;
  [[nodiscard]] double specificInternalEnergyFromTemperatureCode(
      double temperature_k,
      double mean_molecular_weight) const;
  [[nodiscard]] EffectiveMultiphaseEosEntry evaluateDirect(double density_phys_code) const;
  [[nodiscard]] EffectiveMultiphaseEosLookup lookup(double density_phys_code) const;
  [[nodiscard]] std::span<const EffectiveMultiphaseEosEntry> entries() const noexcept;
  [[nodiscard]] std::uint64_t tableHash() const noexcept;
  [[nodiscard]] std::string tableHashHex() const;
  [[nodiscard]] std::string coolingReferenceDescription() const;

 private:
  void build();
  void normalizeThresholdContinuity();
  void computeSignalSpeedDerivative();
  void computeHash();

  EffectiveMultiphaseEosConfig m_config;
  core::UnitSystem m_units;
  CoolingRateProvider m_cooling_provider;
  double m_threshold_density_phys_code = 0.0;
  double m_threshold_density_cgs = 0.0;
  std::vector<EffectiveMultiphaseEosEntry> m_entries;
  std::uint64_t m_table_hash = 0;
};

class EffectiveIsmThermodynamicClosure final : public hydro::HydroThermodynamicClosure {
 public:
  explicit EffectiveIsmThermodynamicClosure(EffectiveMultiphaseEosTable table);
  explicit EffectiveIsmThermodynamicClosure(
      std::shared_ptr<const EffectiveMultiphaseEosTable> table);

  [[nodiscard]] hydro::HydroThermodynamicClosureResult evaluate(
      std::size_t cell_index,
      const hydro::HydroConservedState& conserved,
      const hydro::HydroPrimitiveState& ideal_primitive,
      double scale_factor,
      double redshift) const override;

  [[nodiscard]] const EffectiveMultiphaseEosTable& table() const noexcept;

 private:
  std::shared_ptr<const EffectiveMultiphaseEosTable> m_table;
};

// Integrated code-energy ledger (energy, not energy density): every entry is
// sum(rho_comoving * delta_u * V_comoving) over adjusted cells. Persisted in the
// "effective_multiphase_ism" sidecar with energy_ledger_version=2; sidecars
// without that key carry the pre-repair energy-density sums, which are
// dimensionally invalid and are discarded (not reinterpreted) on restart.
inline constexpr int kEffectiveIsmEnergyLedgerVersion = 2;

struct EffectiveIsmEnergyLedger {
  double energy_added_code = 0.0;
  double energy_removed_code = 0.0;
  double net_energy_adjustment_code = 0.0;
  std::uint64_t adjusted_cell_count = 0;
};

class EffectiveIsmEnergyRelaxationSource final : public hydro::HydroSourceTerm {
 public:
  explicit EffectiveIsmEnergyRelaxationSource(const EffectiveIsmThermodynamicClosure& closure);

  [[nodiscard]] hydro::HydroConservedState sourceForCell(
      std::size_t cell_index,
      const hydro::HydroConservedState& conserved,
      const hydro::HydroPrimitiveState& primitive,
      const hydro::HydroSourceContext& context) const override;

  void resetLedger() const noexcept;
  [[nodiscard]] EffectiveIsmEnergyLedger ledger() const noexcept;

 private:
  const EffectiveIsmThermodynamicClosure& m_closure;
  mutable EffectiveIsmEnergyLedger m_ledger;
};

[[nodiscard]] EffectiveMultiphaseEosConfig makeEffectiveMultiphaseEosConfig(
    const core::PhysicsConfig& physics_config);
[[nodiscard]] CoolingRateProvider makeEffectiveIsmReferenceCoolingProvider(
    const core::PhysicsConfig& physics_config);

}  // namespace cosmosim::physics
