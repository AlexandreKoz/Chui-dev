#pragma once

// CHUI-RT-BH-AGN-REF-v1, RT companion section 6 "Local nonlinear solve"
// (RT04): local stiff radiation-matter exchange S(U) used by both implicit
// IMEX stages, solving  U - alpha S(U) = U_rhs  cell by cell.
//
// Four-force convention (one authority): matter receives +c G^0 and +G,
// radiation -c_hat G^0 and -c c_hat G (variables E, G = F c_hat/c), so the
// closed exchange invariants are E_gas + (c/c_hat) E and p_gas + F/(c c_hat)
// (full-c balance only for c_hat = c).
//
// Models:
//  * grey thermal REFERENCE (constant opacities, LTE emission a T^4):
//    mixed-frame O(v/c) terms following Krumholz et al. (2007, ApJ 667, 626)
//    as used by He, Wibking & Krumholz (2024):
//      G^0 = rho kP (E - a T^4) + rho (kF - 2 kP) v.F / c^2
//      G   = rho kF F/c - rho kP a T^4 v/c - rho kF (v E + v.P)/c
//    (transcribed; audit 02 must verify every coefficient against the paper).
//    F is eliminated linearly (Eddington tensor D = P/E frozen at the stage
//    right-hand side), leaving a 2x2 damped Newton on (T, E) with backtracking
//    that keeps T > 0 and E >= 0; residuals are scaled per conserved quantity;
//    bounded iterations; on failure the caller rolls back and retries with a
//    smaller dt. No energy floor is ever applied silently.
//  * photoionization (frozen absorber densities over the RT interval):
//    chi_g = sum_s n_s sigma_E,g,s;  E' = E_rhs/(1 + alpha c_hat chi_g),
//    G' = G_rhs/(1 + alpha c_hat chi_g); radiation pressure on matter
//    dp = alpha chi_g F'/c (exact for pure absorption). The absorbed energy
//    alpha c chi E' is LEDGERED per absorber for the chemistry authority (it is
//    NOT deposited as heat here; chemistry owns photoheating). Production
//    admission is FAIL_CLOSED until the per-row reconciliation is wired.

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/physics/thermochemistry_operator.hpp"
#include "cosmosim/radiation/radiation_spectrum.hpp"

namespace cosmosim::radiation {

struct RadiationMatterCell {
  double density = 0.0;                   // proper code density
  double specific_internal_energy = 0.0;  // code
  double temperature_per_specific_energy = 0.0;  // K per code specific energy (cell EOS ratio, frozen)
  std::array<double, 3> velocity{};
  std::array<double, k_radiation_absorber_count> absorber_number_density_cgs{};
};

struct RadiationCouplingConstants {
  double c_code = 0.0;
  double c_hat_code = 0.0;
  double kappa_absorption_code = 0.0;   // area / mass (code)
  double kappa_scattering_code = 0.0;
  double radiation_constant_code = 0.0; // a: energy density code per K^4
  double length_cm_per_code = 1.0;      // proper code length in cm (photoionization opacity)
  double energy_density_cgs_per_code = 1.0;
  double density_cgs_per_code = 1.0;    // proper code density -> g cm^-3
  double relative_tolerance = 1.0e-10;
  std::uint32_t max_iterations = 32U;
};

struct RadiationLocalResult {
  bool converged = true;
  std::uint32_t iterations = 0U;
  double residual = 0.0;
  double gas_energy_change = 0.0;           // proper code energy density (thermal + kinetic)
  std::array<double, 3> gas_momentum_change{};  // proper code momentum density
  std::array<double, k_radiation_absorber_count> absorbed_energy_density{};  // photoionization ledger
};

// One cell, all groups: E[g], G[g] (proper code energy density and reduced
// flux) are updated in place from their right-hand-side values.
[[nodiscard]] RadiationLocalResult solveGreyThermalImplicit(RadiationMatterCell& cell, double& e_rad,
                                                            std::array<double, 3>& g_rad, double alpha,
                                                            const RadiationCouplingConstants& k);

[[nodiscard]] RadiationLocalResult solvePhotoionizationImplicit(const RadiationMatterCell& cell,
                                                                const RadiationSpectralDefinition& spectrum,
                                                                std::span<double> e_rad,
                                                                std::span<std::array<double, 3>> g_rad, double alpha,
                                                                const RadiationCouplingConstants& k);

// Transported photo-rates of one cell from its group moments (proper cgs
// energy densities): Gamma_s = sum_g c sigma_N,g,s E_g/<h nu>_g,
// heat_s = sum_g c sigma_N,g,s (E_g/<h nu>_g)(eps_g,s - I_s).
[[nodiscard]] physics::TransportedPhotoRates transportedPhotoRates(const RadiationSpectralDefinition& spectrum,
                                                                   std::span<const double> group_energy_density_cgs);

// Per-row cache of interval-averaged transported rates (implements the
// chemistry provider interface). Non-persistent derived cache: rebuilt from
// the canonical moments; refreshed at every RT commit.
class RadiationPhotoRateCache final : public physics::TransportedPhotoRateProvider {
 public:
  void resize(std::size_t rows) { m_rates.assign(rows, physics::TransportedPhotoRates{}); }
  void set(std::size_t row, const physics::TransportedPhotoRates& r) { m_rates.at(row) = r; }
  [[nodiscard]] physics::TransportedPhotoRates ratesForRow(const core::SimulationState& state,
                                                           std::uint32_t row) const override;
  [[nodiscard]] std::size_t size() const noexcept { return m_rates.size(); }

 private:
  std::vector<physics::TransportedPhotoRates> m_rates;
};

}  // namespace cosmosim::radiation
