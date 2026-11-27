// CHUI-EXT-PHYS-REF-v1 gates VL/VRA/VN (S02): Voigt/Faddeeva accuracy against
// closed forms and an independent quadrature, atomic-line consistency,
// Boltzmann transfer-table parsing/refusal and interpolation, frequency axis.
// NOT RUN in pass 01 (source-only campaign).
#include "../support/test_require.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/spectral/atomic_data.hpp"
#include "cosmosim/spectral/boltzmann_tables.hpp"
#include "cosmosim/spectral/frequency_axis.hpp"

namespace s = cosmosim::spectral;
namespace c = cosmosim::core;
namespace k = cosmosim::core::constants;

namespace {

bool closeRel(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

// EXT-L01-T01: Voigt H against w(i a) = exp(a^2) erfc(a) and the direct
// substitution quadrature over a damping/frequency grid incl. far wings.
void testVoigt() {
  for (double a : {0.25, 0.5, 1.0, 2.0, 4.0}) {
    CHUI_REQUIRE(closeRel(s::voigtH(a, 0.0), std::exp(a * a) * std::erfc(a), 1.0e-12));
  }
  CHUI_REQUIRE(s::voigtH(0.0, 1.3) == std::exp(-1.69));
  for (double a : {1.0e-4, 1.0e-2, 0.1, 1.0, 5.0}) {
    for (double x = -40.0; x <= 40.0; x += 0.73) {
      const double h = s::voigtH(a, x);
      const double ref = s::voigtHQuadrature(a, x);
      CHUI_REQUIRE(std::abs(h - ref) <= 1.0e-9 * ref + 1.0e-14);
    }
  }
  // Lorentzian far wing H -> a / (sqrt(pi) x^2) for |x| >> 1.
  CHUI_REQUIRE(closeRel(s::voigtH(1.0e-3, 300.0), 1.0e-3 / (std::sqrt(k::k_pi) * 300.0 * 300.0), 1.0e-4));
  // Normalization int H dx = sqrt(pi) (tails beyond |x| = 2000 are below 1e-4 relative for a = 0.1).
  double integral = 0.0;
  const double dx = 0.005;
  for (double x = -2000.0; x <= 2000.0; x += dx) integral += s::voigtH(0.1, x) * dx;
  CHUI_REQUIRE(closeRel(integral, std::sqrt(k::k_pi), 1.0e-4));
}

// EXT-L01-T02 / EXT-RA01-T01: line data internal consistency.
void testLines() {
  const s::AtomicLine& lya = s::hydrogenLymanAlpha();
  const double nu = lya.rest_frequency_hz;
  const double a_from_f = 2.0 * k::k_pi * k::k_elementary_charge_si * k::k_elementary_charge_si * nu * nu /
                          (k::k_vacuum_permittivity_si * k::k_electron_mass_si * std::pow(k::k_speed_of_light_si, 3)) *
                          (lya.lower_weight / lya.upper_weight) * lya.oscillator_strength;
  CHUI_REQUIRE(closeRel(a_from_f, lya.einstein_a_s, 5.0e-3));
  // Classical integrated cross-section pi e^2/(m_e c) = 0.02654 cm^2 Hz (CGS) -> 2.654e-6 m^2 Hz.
  CHUI_REQUIRE(closeRel(s::lineIntegratedCrossSectionSi(lya) / lya.oscillator_strength, 2.654e-6, 1.0e-3));
  const s::AtomicLine& hi = s::hydrogen21cm();
  const double t_star = k::k_planck_si * hi.rest_frequency_hz / k::k_boltzmann_si;
  CHUI_REQUIRE(closeRel(t_star, 0.06817, 1.0e-3));
  CHUI_REQUIRE(closeRel(hi.rest_wavelength_m, 0.211061, 1.0e-5));
  // Doppler width for T = 1e4 K: b = sqrt(2kT/m_H) = 12.85 km/s.
  const double width = s::dopplerWidthHz(lya, 1.0e4);
  CHUI_REQUIRE(closeRel(width * k::k_speed_of_light_si / nu, 12.85e3, 2.0e-3));
}

c::LambdaCdmBackground oneNeutrino() {
  c::CosmologyBackgroundConfig cfg;
  cfg.hubble_param = 0.674;
  cfg.omega_matter = 0.315;
  cfg.radiation_species.enabled = true;
  cfg.radiation_species.t_cmb0_k = 2.7255;
  cfg.radiation_species.n_ur = 2.0328;
  cfg.radiation_species.ncdm_masses_ev = {0.06};
  cfg.radiation_species.ncdm_temperature_ratio = 0.71611;
  const double h0 = cfg.hubble_param * k::k_hubble_100_km_s_mpc_si;
  const double rho_crit0 = 3.0 * h0 * h0 / (8.0 * k::k_pi * k::k_newton_g_si);
  const double omega_gamma = k::k_radiation_constant_si * std::pow(2.7255, 4) /
                             (k::k_speed_of_light_si * k::k_speed_of_light_si * rho_crit0);
  cfg.omega_lambda = 1.0 - 0.315 - omega_gamma * (1.0 + 2.0328 * 7.0 / 8.0 * std::pow(4.0 / 11.0, 4.0 / 3.0));
  return c::LambdaCdmBackground(cfg);
}

std::string num(double value) {
  std::ostringstream text;
  text.precision(17);
  text << value;
  return text.str();
}

std::string chuiTable(double z, double ratio_slope) {
  // delta_cb = 1, delta_ncdm = 0.5 + ratio_slope ln(k) -> ratio linear in ln k.
  std::string text =
      "# schema: chui_boltzmann_transfer_v1\n# source_code: synthetic_test\n# source_version: 0\n# gauge: newtonian\n"
      "# redshift: " + num(z) + "\n# k_unit: h/Mpc\n# ncdm_masses_ev: 0.06\n"
      "# columns: k delta_cb delta_ncdm theta_cb theta_ncdm\n";
  for (double kk : {0.01, 0.03, 0.1, 0.3, 1.0, 3.0}) {
    text += num(kk) + " 1.0 " + num(0.5 + ratio_slope * std::log(kk)) + " -0.2 -0.1\n";
  }
  return text;
}

// EXT-N02-T01: chui_v1/CLASS/CAMB parsing, refusals and interpolation.
void testBoltzmannTables() {
  const c::LambdaCdmBackground bg = oneNeutrino();
  s::BoltzmannLoadRequest req;
  req.path = "synthetic_z49";
  req.format = c::BoltzmannTableFormatConfig::kChuiV1;
  req.declared_gauge = c::BoltzmannGaugeConfig::kNewtonian;
  req.declared_redshift = 49.0;
  req.omega_baryon = 0.049;
  auto early = s::parseBoltzmannTransferEpoch(chuiTable(49.0, 0.0), req, bg);
  req.path = "synthetic_z9";
  req.declared_redshift = 9.0;
  auto late = s::parseBoltzmannTransferEpoch(chuiTable(9.0, 0.01), req, bg);
  CHUI_REQUIRE(early.hasVelocities() && late.k_h_mpc.size() == 6U);
  const s::BoltzmannTransferSet set({early, late});
  // At a of the first epoch: ratio 0.5 everywhere.
  CHUI_REQUIRE(closeRel(set.ncdmToColdRatio(0.2, 1.0 / 50.0), 0.5, 1.0e-12));
  // At the late epoch the ratio is linear in ln k, so ln-k interpolation is exact.
  CHUI_REQUIRE(closeRel(set.ncdmToColdRatio(0.2, 0.1), 0.5 + 0.01 * std::log(0.2), 1.0e-12));
  // Midway in ln a: average of the two epochs.
  const double a_mid = std::exp(0.5 * (std::log(1.0 / 50.0) + std::log(0.1)));
  CHUI_REQUIRE(closeRel(set.ncdmToColdRatio(1.0, a_mid), 0.5, 1.0e-12));
  CHUI_REQUIRE(closeRel(set.thetaNcdmOverDeltaCold(0.5, 0.05), -0.1, 1.0e-12));
  bool threw = false;
  try {
    (void)set.ncdmToColdRatio(10.0, 0.05);  // beyond k_max
  } catch (const std::out_of_range&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  threw = false;
  try {
    (void)set.ncdmToColdRatio(0.1, 0.5);  // beyond the last epoch
  } catch (const std::out_of_range&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  // Gauge mismatch and redshift mismatch are refused.
  req.declared_gauge = c::BoltzmannGaugeConfig::kSynchronous;
  threw = false;
  try {
    (void)s::parseBoltzmannTransferEpoch(chuiTable(9.0, 0.0), req, bg);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  req.declared_gauge = c::BoltzmannGaugeConfig::kNewtonian;
  req.declared_redshift = 10.0;
  threw = false;
  try {
    (void)s::parseBoltzmannTransferEpoch(chuiTable(9.0, 0.0), req, bg);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  // CLASS layout with one ncdm species; delta_cb is the Omega-weighted mean.
  const std::string class_text =
      "# Transfer functions T_i(k) for adiabatic (AD) mode (normalized to initial curvature=1) at redshift z=9\n"
      "#    1:k (h/Mpc)    2:d_g    3:d_b    4:d_cdm    5:d_ur    6:d_ncdm[0]    7:d_tot\n"
      "0.01 1 2 4 1 3 4\n0.1 1 2 4 1 3 4\n1 1 2 4 1 3 4\n2 1 2 4 1 3 4\n";
  s::BoltzmannLoadRequest class_req = req;
  class_req.format = c::BoltzmannTableFormatConfig::kClassTk;
  class_req.declared_gauge = c::BoltzmannGaugeConfig::kSynchronous;
  class_req.declared_redshift = 9.0;
  const auto class_epoch = s::parseBoltzmannTransferEpoch(class_text, class_req, bg);
  const double omega_cdm = bg.coldMatterDensityParameter0() - 0.049;
  const double expected_cb = (omega_cdm * 4.0 + 0.049 * 2.0) / (omega_cdm + 0.049);
  CHUI_REQUIRE(closeRel(class_epoch.delta_cb[1], expected_cb, 1.0e-14));
  CHUI_REQUIRE(closeRel(class_epoch.delta_ncdm[1], 3.0, 1.0e-14));
  CHUI_REQUIRE(!class_epoch.hasVelocities());
  // CAMB: 13 columns, delta_cb = no_nu (column 8), delta_ncdm = mass_nu (column 6).
  const std::string camb_text =
      "0.01 1 2 3 4 5 6 7 8 9 10 11 12\n0.1 1 2 3 4 5 6 7 8 9 10 11 12\n"
      "1 1 2 3 4 5 6 7 8 9 10 11 12\n2 1 2 3 4 5 6 7 8 9 10 11 12\n";
  s::BoltzmannLoadRequest camb_req = class_req;
  camb_req.format = c::BoltzmannTableFormatConfig::kCambTransfer;
  const auto camb_epoch = s::parseBoltzmannTransferEpoch(camb_text, camb_req, bg);
  CHUI_REQUIRE(camb_epoch.delta_cb[0] == 7.0 && camb_epoch.delta_ncdm[0] == 5.0);
}

// EXT-RA05-T01: frequency axis Gauss-Legendre channel averages.
void testFrequencyAxis() {
  const s::FrequencyAxis axis = s::FrequencyAxis::uniform(1.0e8, 2.0e8, 10, 4);
  CHUI_REQUIRE(axis.channelCount() == 10U && closeRel(axis.widthHz(3), 1.0e7, 1.0e-14));
  double weight_sum = 0.0;
  for (double w : axis.nodeWeights()) weight_sum += w;
  CHUI_REQUIRE(closeRel(weight_sum, 1.0, 1.0e-14));
  // 4-node rule averages a cubic exactly: mean of nu^3 over [lo, hi] = (hi^4 - lo^4) / (4 (hi - lo)).
  const auto nodes = axis.nodeFrequenciesHz(2);
  double avg = 0.0;
  for (std::size_t i = 0; i < nodes.size(); ++i) avg += axis.nodeWeights()[i] * std::pow(nodes[i] / 1.0e8, 3);
  const double lo = 1.2;
  const double hi = 1.3;
  CHUI_REQUIRE(closeRel(avg, (std::pow(hi, 4) - std::pow(lo, 4)) / (4.0 * (hi - lo)), 1.0e-13));
}

}  // namespace

int main() {
  testVoigt();
  testLines();
  testBoltzmannTables();
  testFrequencyAxis();
  return 0;
}
