// CHUI-EXT-PHYS-REF-v1 gates V00/VN/VC (S01, N01): shared cosmological and
// radiation background authority. Expected values are derived independently
// of the implementation (closed-form Fermi-Dirac moments, Einstein-de Sitter
// distances, number-density formula for the nonrelativistic limit).
// NOT RUN in pass 01 (source-only campaign).
#include "../support/test_require.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/cosmology.hpp"

namespace c = cosmosim::core;
namespace k = cosmosim::core::constants;

namespace {

bool closeRel(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

// EXT-N01-T01: Fermi-Dirac moments, analytic limits and table accuracy.
void testFermiDiracMoments() {
  const double pi4 = std::pow(k::k_pi, 4);
  CHUI_REQUIRE(closeRel(c::fermiDiracEnergyIntegralQuadrature(0.0), 7.0 * pi4 / 120.0, 1.0e-13));
  CHUI_REQUIRE(closeRel(c::fermiDiracPressureIntegralQuadrature(0.0), 7.0 * pi4 / 360.0, 1.0e-13));
  CHUI_REQUIRE(closeRel(c::fermiDiracNumberIntegral(), 1.5 * k::k_zeta_3, 1.0e-15));
  // Table vs direct quadrature across the transition (relative 1e-9).
  for (double y = 1.1e-5; y < 9.0e5; y *= 1.37) {
    CHUI_REQUIRE(closeRel(c::fermiDiracEnergyIntegral(y), c::fermiDiracEnergyIntegralQuadrature(y), 1.0e-9));
    CHUI_REQUIRE(closeRel(c::fermiDiracPressureIntegral(y), c::fermiDiracPressureIntegralQuadrature(y), 1.0e-9));
    // Physical bounds: 0 < P <= rho/3.
    CHUI_REQUIRE(c::fermiDiracPressureIntegral(y) <= c::fermiDiracEnergyIntegral(y) / 3.0 * (1.0 + 1.0e-12));
  }
  // Asymptotic branches agree with quadrature at the table edges.
  CHUI_REQUIRE(closeRel(c::fermiDiracEnergyIntegral(0.5e-5), c::fermiDiracEnergyIntegralQuadrature(0.5e-5), 1.0e-12));
  CHUI_REQUIRE(closeRel(c::fermiDiracEnergyIntegral(2.0e6), c::fermiDiracEnergyIntegralQuadrature(2.0e6), 1.0e-10));
  CHUI_REQUIRE(closeRel(c::fermiDiracPressureIntegral(2.0e6), c::fermiDiracPressureIntegralQuadrature(2.0e6), 1.0e-8));
  // Nonrelativistic limit rho -> m n: I_rho(y)/y -> (3/2) zeta(3).
  CHUI_REQUIRE(closeRel(c::fermiDiracEnergyIntegral(1.0e8) / 1.0e8, 1.5 * k::k_zeta_3, 1.0e-12));
  // Mean momentum <x> = (7 pi^4/120)/((3/2) zeta 3) = 3.1513...
  CHUI_REQUIRE(closeRel(c::fermiDiracMeanMomentumOverTemperature(), (7.0 * pi4 / 120.0) / (1.5 * k::k_zeta_3), 1.0e-15));
  bool threw = false;
  try {
    (void)c::fermiDiracEnergyIntegral(-1.0);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

c::CosmologyBackgroundConfig planckLike(double n_ur, std::vector<double> masses, double ratio) {
  c::CosmologyBackgroundConfig cfg;
  cfg.hubble_param = 0.674;
  cfg.omega_matter = 0.315;
  cfg.radiation_species.enabled = true;
  cfg.radiation_species.t_cmb0_k = 2.7255;
  cfg.radiation_species.n_ur = n_ur;
  cfg.radiation_species.ncdm_masses_ev = std::move(masses);
  cfg.radiation_species.ncdm_temperature_ratio = cfg.radiation_species.ncdm_masses_ev.empty() ? 0.0 : ratio;
  // Independent photon density: a_rad T^4 / c^2 / rho_crit0.
  const double h0 = cfg.hubble_param * k::k_hubble_100_km_s_mpc_si;
  const double rho_crit0 = 3.0 * h0 * h0 / (8.0 * k::k_pi * k::k_newton_g_si);
  const double omega_gamma = k::k_radiation_constant_si * std::pow(2.7255, 4) /
                             (k::k_speed_of_light_si * k::k_speed_of_light_si * rho_crit0);
  const double omega_ur = n_ur * 7.0 / 8.0 * std::pow(4.0 / 11.0, 4.0 / 3.0) * omega_gamma;
  cfg.omega_lambda = 1.0 - cfg.omega_matter - omega_gamma - omega_ur;
  return cfg;
}

// EXT-S01-T01: photon density and CMB number density against literature values
// for T = 2.7255 K (Omega_gamma h^2 = 2.4728e-5; n_gamma = 410.7 cm^-3).
void testPhotonBackground() {
  const c::LambdaCdmBackground bg(planckLike(0.0, {}, 0.0));
  const double h2 = 0.674 * 0.674;
  CHUI_REQUIRE(closeRel(bg.photonDensityParameter0() * h2, 2.4728e-5, 2.0e-4));
  CHUI_REQUIRE(closeRel(bg.cmbPhotonNumberDensitySi(1.0) * 1.0e-6, 410.7, 1.0e-3));
  CHUI_REQUIRE(closeRel(bg.cmbTemperatureK(0.25), 4.0 * 2.7255, 1.0e-15));
  CHUI_REQUIRE(closeRel(bg.cmbEnergyDensitySi(0.5), 16.0 * bg.cmbEnergyDensitySi(1.0), 1.0e-14));
  CHUI_REQUIRE(closeRel(bg.eFactor(1.0), 1.0, 1.0e-10));
  // Radiation domination: E^2 a^4 -> Omega_gamma as a -> 0.
  const double a = 1.0e-7;
  CHUI_REQUIRE(closeRel(bg.eFactor(a) * bg.eFactor(a) * std::pow(a, 4), bg.photonDensityParameter0(), 1.0e-3));
  // Disabled background refuses a CMB temperature instead of inventing one.
  bool threw = false;
  try {
    (void)c::LambdaCdmBackground().cmbTemperatureK(1.0);
  } catch (const std::logic_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

// EXT-N01-T02: massive eigenstates — massless and nonrelativistic limits,
// no double counting in Omega_matter, CLASS 93.14 eV reference relation.
void testMassiveNeutrinos() {
  const double ratio = 0.71611;  // CLASS default T_ncdm/T_cmb
  const c::LambdaCdmBackground bg(planckLike(2.0328, {0.06}, ratio));
  CHUI_REQUIRE(bg.massiveNeutrinoCount() == 1U);
  // Massless limit at early times: (7/8) r^4 Omega_gamma a^-4.
  const double a_early = 1.0e-9;
  const double massless = 7.0 / 8.0 * std::pow(ratio, 4) * bg.photonDensityParameter0() / std::pow(a_early, 4);
  CHUI_REQUIRE(closeRel(bg.massiveNeutrinoDensityParameter(0, a_early), massless, 1.0e-6));
  CHUI_REQUIRE(closeRel(bg.massiveNeutrinoPressureParameter(0, a_early), massless / 3.0, 1.0e-6));
  // Nonrelativistic limit today from the number density, independent formula:
  // rho = m * g/(2 pi^2) (k T_nu/(hbar c))^3 (3/2) zeta(3), g = 2.
  const double kt = k::k_boltzmann_si * ratio * 2.7255;
  const double inv_len = kt / (k::k_reduced_planck_si * k::k_speed_of_light_si);
  const double n_nu = 2.0 / (2.0 * k::k_pi * k::k_pi) * inv_len * inv_len * inv_len * 1.5 * k::k_zeta_3;
  const double rho_nr = 0.06 * k::k_electron_volt_si * n_nu / (k::k_speed_of_light_si * k::k_speed_of_light_si);
  CHUI_REQUIRE(closeRel(bg.massiveNeutrinoDensityParameter(0, 1.0), rho_nr / bg.criticalDensity0Si(), 1.0e-5));
  // CLASS convention reference: Omega_ncdm h^2 = m / 93.14 eV (0.1% level).
  CHUI_REQUIRE(closeRel(bg.massiveNeutrinoDensityParameter(0, 1.0) * 0.674 * 0.674, 0.06 / 93.14, 1.0e-3));
  // Omega_cb + Omega_ncdm,0 == Omega_matter exactly (no double counting).
  CHUI_REQUIRE(closeRel(bg.coldMatterDensityParameter0() + bg.totalMassiveNeutrinoDensityParameter(1.0), 0.315, 1.0e-14));
  CHUI_REQUIRE(closeRel(bg.eFactor(1.0), 1.0, 1.0e-10));
  // Validity diagnostic: <p>/(m c) = 3.1514 k T_nu / (m c^2).
  const double beta = bg.massiveNeutrinoMeanMomentumOverMass(0, 1.0);
  CHUI_REQUIRE(closeRel(beta, c::fermiDiracMeanMomentumOverTemperature() * kt / (0.06 * k::k_electron_volt_si), 1.0e-12));
  // Massless-only configuration: a_eq = (Omega_gamma + Omega_ur) / Omega_cb exactly.
  const c::LambdaCdmBackground massless_bg(planckLike(3.046, {}, 0.0));
  const double expected_aeq = (massless_bg.photonDensityParameter0() + massless_bg.masslessNeutrinoDensityParameter0()) /
                              massless_bg.coldMatterDensityParameter0();
  CHUI_REQUIRE(closeRel(massless_bg.matterRadiationEqualityScaleFactor(), expected_aeq, 1.0e-10));
}

// EXT-S01-T02 / V00: disabled species reproduce the legacy formula bit-for-bit
// and refuse inert values; wrong closure names the required Omega_lambda.
void testLegacyEquivalenceAndRefusals() {
  c::CosmologyBackgroundConfig legacy;
  legacy.hubble_param = 0.7;
  legacy.omega_matter = 0.3;
  legacy.omega_lambda = 0.7;
  const c::LambdaCdmBackground bg(legacy);
  for (double a : {0.01, 0.1, 0.5, 1.0}) {
    const double expected = std::sqrt(0.3 / (a * a * a) + 0.7);
    CHUI_REQUIRE(bg.eFactor(a) == expected);
  }
  c::CosmologyBackgroundConfig inert = legacy;
  inert.radiation_species.t_cmb0_k = 2.7255;  // enabled == false
  bool threw = false;
  try {
    const c::LambdaCdmBackground rejected(inert);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  c::CosmologyBackgroundConfig open_closure = planckLike(3.046, {}, 0.0);
  open_closure.omega_lambda = 0.685;  // ignores radiation
  std::string message;
  try {
    const c::LambdaCdmBackground rejected(open_closure);
  } catch (const std::invalid_argument& error) {
    message = error.what();
  }
  CHUI_REQUIRE(message.find("strict closure requires omega_lambda") != std::string::npos);
  c::CosmologyBackgroundConfig double_count = planckLike(3.046, {}, 0.0);
  double_count.omega_radiation = 9.0e-5;
  threw = false;
  try {
    const c::LambdaCdmBackground rejected(double_count);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

// EXT-S01-T03: distances. Einstein-de Sitter chi(a) = 2 c/H0 (1 - sqrt a).
void testDistances() {
  c::CosmologyBackgroundConfig eds;
  eds.hubble_param = 0.7;
  eds.omega_matter = 1.0;
  eds.omega_lambda = 0.0;
  const c::LambdaCdmBackground bg(eds);
  const double dh = k::k_speed_of_light_si / bg.hubble0Si();
  for (double a : {0.9, 0.5, 0.2, 0.05}) {
    const double expected = 2.0 * dh * (1.0 - std::sqrt(a));
    CHUI_REQUIRE(closeRel(bg.comovingDistanceSi(a), expected, 1.0e-9));
    CHUI_REQUIRE(closeRel(bg.angularDiameterDistanceSi(a), a * expected, 1.0e-9));
    CHUI_REQUIRE(closeRel(bg.luminosityDistanceSi(a), expected / a, 1.0e-9));
    CHUI_REQUIRE(closeRel(bg.scaleFactorAtComovingDistanceSi(expected, 0.01), a, 1.0e-9));
  }
  CHUI_REQUIRE(bg.comovingDistanceSi(1.0) == 0.0);
  // EdS conformal time: eta(a) = 2 sqrt(a) / H0.
  CHUI_REQUIRE(closeRel(bg.conformalTimeIntervalSi(0.25, 1.0), 2.0 * (1.0 - 0.5) / bg.hubble0Si(), 1.0e-9));
}

// EXT-S01-T04 / VI: .param.txt admission, normalized round trip and refusals.
std::string cosmologicalDeck(const std::string& extra_cosmology, const std::string& extra_sections) {
  const c::CosmologyBackgroundConfig reference = planckLike(2.0328, {0.06}, 0.71611);
  char lambda[64];
  std::snprintf(lambda, sizeof(lambda), "%.17g", reference.omega_lambda);
  return "[mode]\nmode = cosmo_cube\n[numerics]\na_begin = 0.02\na_end = 0.5\n[cosmology]\nomega_matter = 0.315\n"
         "omega_baryon = 0.049\nhubble_param = 0.674\nomega_lambda = " +
         std::string(lambda) +
         "\nradiation_background = photons_neutrinos\nt_cmb0_k = 2.7255\nn_ur = 2.0328\nncdm_masses_ev = 0.06\n"
         "ncdm_temperature_ratio = 0.71611\n" +
         extra_cosmology + extra_sections;
}

bool rejects(const std::string& text, const std::string& needle) {
  try {
    (void)c::loadFrozenConfigFromString(text, "ext_background_case");
  } catch (const c::ConfigError& error) {
    return std::string(error.what()).find(needle) != std::string::npos;
  }
  return false;
}

void testConfigAdmission() {
  // Disabled default: no new normalized text.
  const auto legacy = c::loadFrozenConfigFromString("mode.mode = isolated_galaxy\n", "legacy");
  CHUI_REQUIRE(legacy.normalized_text.find("radiation_background") == std::string::npos);
  CHUI_REQUIRE(legacy.normalized_text.find("[neutrinos]") == std::string::npos);
  // Massive eigenstates without an explicit perturbation treatment are refused.
  CHUI_REQUIRE(rejects(cosmologicalDeck("", ""), "declare the perturbation treatment explicitly"));
  const auto frozen = c::loadFrozenConfigFromString(cosmologicalDeck("", "[neutrinos]\nmethod = background_only\n"), "bg");
  CHUI_REQUIRE(frozen.config.cosmology.radiation_species.enabled);
  CHUI_REQUIRE(frozen.normalized_text.find("radiation_background = photons_neutrinos") != std::string::npos);
  const auto again = c::loadFrozenConfigFromString(frozen.normalized_text, "round_trip");
  CHUI_REQUIRE(again.normalized_text == frozen.normalized_text);
  // The background constructed by the runtime factory equals the reference.
  const c::LambdaCdmBackground runtime(c::cosmologyBackgroundConfigFromSimulation(frozen.config));
  CHUI_REQUIRE(closeRel(runtime.eFactor(0.1), c::LambdaCdmBackground(planckLike(2.0328, {0.06}, 0.71611)).eFactor(0.1), 1.0e-12));
  // Inert keys and unsupported regimes.
  CHUI_REQUIRE(rejects("mode.mode = isolated_galaxy\ncosmology.t_cmb0_k = 2.7\n", "would have no effect"));
  CHUI_REQUIRE(rejects("mode.mode = isolated_galaxy\ncosmology.radiation_background = photons_neutrinos\n"
                       "cosmology.t_cmb0_k = 2.7\n", "cosmological comoving mode"));
  std::string early = cosmologicalDeck("", "[neutrinos]\nmethod = background_only\n");
  early.replace(early.find("a_begin = 0.02"), 14, "a_begin = 1e-05");
  CHUI_REQUIRE(rejects(early, "matter-radiation equality"));
  CHUI_REQUIRE(rejects(cosmologicalDeck("", "[neutrinos]\nmethod = delta_f\n"), "FAIL_CLOSED"));
  CHUI_REQUIRE(rejects(cosmologicalDeck("", "[neutrinos]\nmethod = transfer_ratio\n"), "requires neutrinos.table_paths"));
}

}  // namespace

int main() {
  testFermiDiracMoments();
  testPhotonBackground();
  testMassiveNeutrinos();
  testLegacyEquivalenceAndRefusals();
  testDistances();
  testConfigAdmission();
  return 0;
}
