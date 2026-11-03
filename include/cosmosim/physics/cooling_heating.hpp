#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::physics {

struct CoolingTableProvenance {
  std::string source_path;
  std::string source_format;
  std::string table_label;
  std::uint64_t row_count = 0;
};

struct CoolingRateQuery {
  double temperature_k = 1.0e4;
  double hydrogen_number_density_cgs = 1.0e-4;
  double metallicity_mass_fraction = 0.0;
  double redshift = 0.0;
};

struct CoolingRateResult {
  double cooling_rate_erg_cm3_s = 0.0;
  double heating_rate_erg_cm3_s = 0.0;

  [[nodiscard]] double netCoolingRateErgCm3S() const;
};

class MetalLineCoolingTable {
 public:
  [[nodiscard]] static MetalLineCoolingTable loadFromTextFile(
      const std::string& path,
      const std::string& table_label = "metal_line");

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] double lookupCoolingRateErgCm3S(double temperature_k) const;
  [[nodiscard]] const CoolingTableProvenance& provenance() const noexcept;

 private:
  CoolingTableProvenance m_provenance;
  std::vector<double> m_log10_temperature_k;
  std::vector<double> m_log10_cooling_rate_erg_cm3_s;
};

struct CoolingModelConfig {
  core::UvBackgroundModel uv_background_model = core::UvBackgroundModel::kHm12;
  core::SelfShieldingModel self_shielding_model = core::SelfShieldingModel::kNone;
  double temperature_floor_k = 100.0;
  double primordial_helium_mass_fraction = 0.24;
  bool enable_metal_line_cooling = false;
  double max_fractional_energy_change_per_substep = 0.1;
  int max_subcycles = 16;
  double temperature_per_internal_energy_k = 1.0e6;
};

class CoolingRateProvider {
 public:
  explicit CoolingRateProvider(CoolingModelConfig config);

  void setMetalLineTable(MetalLineCoolingTable table);
  [[nodiscard]] bool hasMetalLineTable() const noexcept;
  [[nodiscard]] CoolingRateResult lookupRates(const CoolingRateQuery& query) const;
  [[nodiscard]] const CoolingModelConfig& config() const noexcept;
  [[nodiscard]] std::optional<CoolingTableProvenance> metalLineTableProvenance() const;

 private:
  CoolingModelConfig m_config;
  std::optional<MetalLineCoolingTable> m_metal_line_table;
};

struct CoolingIntegrationDiagnostics {
  // Step attempts (accepted + rejected); bounded by max_subcycles.
  int subcycles_used = 0;
  double suggested_next_dt_code = 0.0;
  // True iff the attempt budget ran out before the interval closed (F2).
  bool hit_subcycle_limit = false;
  int accepted_substeps = 0;
  int rejected_attempts = 0;
  int floor_activations = 0;
  double max_error_estimate = 0.0;
};

// CHUI Part 2 (audit F2): `completed` is true only when the whole requested
// interval was integrated. An incomplete result returns the UNCHANGED input
// energy and must be rejected by the caller.
struct CoolingIntegrationResult {
  double specific_internal_energy_code = 0.0;
  bool completed = false;
  double requested_dt_code = 0.0;
  double completed_dt_code = 0.0;
  CoolingIntegrationDiagnostics diagnostics;
};

class CoolingSourceIntegrator {
 public:
  explicit CoolingSourceIntegrator(
      double energy_floor_code,
      double specific_energy_rate_cgs_to_code = 1.0);

  [[nodiscard]] CoolingIntegrationResult integrateSpecificInternalEnergy(
      double specific_internal_energy_code,
      double mass_density_physical_cgs,
      double dt_code,
      const CoolingRateQuery& query,
      const CoolingRateProvider& rate_provider) const;

 private:
  double m_energy_floor_code = 0.0;
  double m_specific_energy_rate_cgs_to_code = 1.0;
};

class CoolingHeatingSource final : public hydro::HydroSourceTerm {
 public:
  CoolingHeatingSource(CoolingRateProvider provider, CoolingSourceIntegrator integrator);

  [[nodiscard]] hydro::HydroConservedState sourceForCell(
      std::size_t cell_index,
      const hydro::HydroConservedState& conserved,
      const hydro::HydroPrimitiveState& primitive,
      const hydro::HydroSourceContext& context) const override;

 private:
  CoolingRateProvider m_rate_provider;
  CoolingSourceIntegrator m_integrator;
};

}  // namespace cosmosim::physics
