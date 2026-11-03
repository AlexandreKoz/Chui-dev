#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/chemical_network.hpp"
#include "cosmosim/physics/metal_cooling_table.hpp"
#include "cosmosim/physics/radiation_background.hpp"
#include "cosmosim/physics/source_interval.hpp"

namespace cosmosim::physics {

// ---------------------------------------------------------------------------
// CHUI Part 2 one-zone thermochemistry (isochoric source step).
//
// Equations (proper cgs, density rho fixed over the source interval):
//   dn_i/dt = sum_r nu_ir R_r(n, T)                       (network backends)
//   du/dt   = [Gamma(n, T) - Lambda_rad(n, T) - Lambda_chem(n, T)] / rho
//   T       = (gamma - 1) rho u / (k_B n_particles),
//   n_particles = sum_i n_i + n_e + n_metal_particles
// with n_e from charge neutrality. The evolved energy is the THERMAL specific
// internal energy u (canonical CHUI authority, GasCellSidecar::
// internal_energy_code). Chemical binding energy (ionization/dissociation) is
// NOT part of u: collisional ionization/dissociation move energy from u into
// the chemical reservoir (Lambda_chem) and the reservoir's change is ledgered
// exactly from the species change (chemical_binding_delta_erg). Radiative
// losses leave the gas. gamma = 5/3 for all particles (H2 rotational modes
// not represented — stated approximation consistent with the hydro EOS).
//
// Integrator: backward (implicit) Euler with a full Newton iteration
// (finite-difference Jacobian, dense LU with partial pivoting on <= 11
// unknowns), local error from step doubling (one step h vs two steps h/2,
// Richardson-corrected solution when positive). Per-component scaled
// tolerances: abundance atol per H nucleus, energy relative. The interval is
// closed EXACTLY (long double clock, last step = remainder). Exhausting
// max_substeps, Newton failure at the minimum step, a non-finite rate or a
// negative species beyond tolerance returns kRetryShorterInterval with the
// completed sub-interval as the legal shorter-step request — never a partial
// success (audit F2). Linear element/charge invariants are preserved by
// backward Euler; a final projection removes round-off and is ledgered.
// ---------------------------------------------------------------------------

enum class MetalCoolingMode : std::uint8_t {
  kNone = 0,
  // C+ 158 um + O I 63 um only (cold gas, T < switch temperature).
  kFineStructureCold = 1,
  // Element table only (BLOCKED BY DATA until supplied).
  kElementTable = 2,
  // Table above the switch temperature, fine structure below, smooth blend
  // over +-0.15 dex (weights sum to one: no double counting).
  kElementTableAndFineStructure = 3,
};

struct ThermochemistryConfig {
  ThermochemistryBackend backend = ThermochemistryBackend::kFixedHybrid;
  double relative_tolerance = 1.0e-4;
  double abundance_absolute_tolerance = 1.0e-12;  // per H nucleus
  std::uint32_t max_substeps = 4096U;
  std::uint32_t max_newton_iterations = 12U;
  double minimum_step_fraction = 1.0e-14;  // of the requested interval
  double temperature_floor_k = 10.0;
  double adiabatic_index = 5.0 / 3.0;
  // Analytic reference (REFERENCE ONLY): du/dt = -(u - u_eq)/t_relax + h0.
  double analytic_relaxation_time_s = 3.15576e13;  // 1 Myr
  double analytic_equilibrium_temperature_k = 1.0e4;
  double analytic_constant_heating_erg_g_s = 0.0;
  // Radiation.
  const UvBackgroundTable* uv_background = nullptr;  // null => no UVB
  bool self_shielding = true;
  double isrf_g0_habing = 0.0;
  double shielding_length_max_cm = 3.0857e20;  // 100 pc cap on the Jeans length
  // Metals.
  MetalCoolingMode metal_mode = MetalCoolingMode::kNone;
  const MetalCoolingTable* metal_table = nullptr;
  MetalTableRangePolicy metal_table_range = MetalTableRangePolicy::kReject;
  double fine_structure_switch_temperature_k = 1.0e4;
  bool cosmological = false;  // enables Compton exchange with the CMB
};

struct ThermochemistryCellInput {
  double mass_density_cgs = 1.0e-24;               // proper, gas + dust
  double specific_internal_energy_erg_g = 1.0e12;  // thermal
  // Gas-phase element mass fractions (relative to the total cell mass, dust
  // excluded from these numerators).
  std::array<double, k_chemical_element_count> gas_element_mass_fraction{};
  // Species number densities (network backends); ignored otherwise.
  NetworkState species_cm3{};
  double dust_to_gas_relative_mw = 0.0;
  double redshift = 0.0;
};

struct ThermochemistryCellOutput {
  SourceAdvanceStatus status = SourceAdvanceStatus::kCompleted;
  double specific_internal_energy_erg_g = 0.0;
  NetworkState species_cm3{};
  double temperature_k = 0.0;
  double mean_molecular_weight = 0.0;
  double electron_per_hydrogen = 0.0;
  double completed_seconds = 0.0;
  double suggested_retry_seconds = 0.0;
  std::string diagnostic;
};

class ThermochemistrySolver {
 public:
  explicit ThermochemistrySolver(ThermochemistryConfig config);

  [[nodiscard]] const ThermochemistryConfig& config() const noexcept { return m_config; }
  [[nodiscard]] std::size_t speciesCount() const noexcept { return m_species_count; }

  // Advance one cell over exactly dt_seconds. Ledger energies are per unit
  // volume [erg cm^-3]; the caller multiplies by the proper cell volume.
  void advanceCell(const ThermochemistryCellInput& input, double dt_seconds,
                   ThermochemistryCellOutput& output, SourceLedger& ledger_per_volume) const;

  // Ionization-equilibrium H/He state at fixed T (KWH96 iteration).
  [[nodiscard]] NetworkState equilibriumPrimordialState(double temperature_k, double n_h, double n_he,
                                                        const RateEnvironment& env) const;
  [[nodiscard]] static double meanMolecularWeight(double mass_density_cgs, double n_particles) noexcept;
  // Total particle number density of a state (species + electrons + metal
  // particles not represented by stored species).
  [[nodiscard]] double particleDensity(const ThermochemistryCellInput& input, const NetworkState& n,
                                       double n_e) const;
  [[nodiscard]] RateEnvironment environment(const ThermochemistryCellInput& input, const NetworkState& n,
                                            double temperature_k) const;
  // Volumetric net heating-cooling for diagnostics/tests.
  [[nodiscard]] CoolingHeatingTerms coolingHeating(const ThermochemistryCellInput& input,
                                                   const NetworkState& n, double temperature_k) const;

 private:
  ThermochemistryConfig m_config;
  const NetworkManifest* m_manifest = nullptr;
  std::size_t m_species_count = 0;
};

}  // namespace cosmosim::physics
