// CHUI-EXT-PHYS-REF-v1 gate VN (N03-N06): free-streaming kernel, superconformal
// time, transfer-ratio multiplier, linear-response cold-limit against an
// independent RK4 integration of the joint growth ODE, persistence/replay/
// capacity refusal, delta-f weights and signed-mass refusal, nonrelativistic
// drift, FD sampling moments, species spectra normalization and noise.
// NOT RUN in pass 01 (source-only campaign).
#include "../support/test_require.hpp"

#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/neutrinos/neutrino_linear_response.hpp"
#include "cosmosim/neutrinos/neutrino_particles.hpp"
#include "cosmosim/neutrinos/species_power_spectrum.hpp"
#include "cosmosim/spectral/boltzmann_tables.hpp"

namespace c = cosmosim::core;
namespace k = cosmosim::core::constants;
namespace nu = cosmosim::neutrinos;
namespace s = cosmosim::spectral;

namespace {

bool closeRel(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

std::string num(double v) {
  std::ostringstream text;
  text.precision(17);
  text << v;
  return text.str();
}

c::CosmologyBackgroundConfig neutrinoBackground(double mass_ev) {
  c::CosmologyBackgroundConfig cfg;
  cfg.hubble_param = 0.674;
  cfg.omega_matter = 0.3;
  cfg.radiation_species.enabled = true;
  cfg.radiation_species.t_cmb0_k = 2.7255;
  cfg.radiation_species.n_ur = 0.0;
  cfg.radiation_species.ncdm_masses_ev = {mass_ev};
  cfg.radiation_species.ncdm_temperature_ratio = c::instantaneousDecouplingNeutrinoTemperatureRatio();
  const double h0 = cfg.hubble_param * k::k_hubble_100_km_s_mpc_si;
  const double rho_crit0 = 3.0 * h0 * h0 / (8.0 * k::k_pi * k::k_newton_g_si);
  const double omega_gamma = k::k_radiation_constant_si * std::pow(2.7255, 4) /
                             (k::k_speed_of_light_si * k::k_speed_of_light_si * rho_crit0);
  cfg.omega_lambda = 1.0 - cfg.omega_matter - omega_gamma;
  return cfg;
}

std::string epochText(double z, double delta_ratio, double theta_ratio_class) {
  std::string text =
      "# schema: chui_boltzmann_transfer_v1\n# source_code: synthetic_test\n# source_version: 0\n# gauge: newtonian\n"
      "# redshift: " + num(z) + "\n# k_unit: h/Mpc\n# ncdm_masses_ev: " + num(1.0) +
      "\n# columns: k delta_cb delta_ncdm theta_cb theta_ncdm\n";
  for (double kk : {1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0, 10.0}) {
    text += num(kk) + " 1 " + num(delta_ratio) + " " + num(theta_ratio_class) + " " + num(theta_ratio_class) + "\n";
  }
  return text;
}

s::BoltzmannTransferSet syntheticSet(const c::LambdaCdmBackground& bg, double z_early, double z_late,
                                     double delta_ratio, double theta_ratio_class) {
  s::BoltzmannLoadRequest req;
  req.format = c::BoltzmannTableFormatConfig::kChuiV1;
  req.declared_gauge = c::BoltzmannGaugeConfig::kNewtonian;
  req.omega_baryon = 0.05;
  req.path = "early";
  req.declared_redshift = z_early;
  auto e0 = s::parseBoltzmannTransferEpoch(epochText(z_early, delta_ratio, theta_ratio_class), req, bg);
  req.path = "late";
  req.declared_redshift = z_late;
  auto e1 = s::parseBoltzmannTransferEpoch(epochText(z_late, delta_ratio, theta_ratio_class), req, bg);
  return s::BoltzmannTransferSet({e0, e1});
}

// EXT-N03-T01: FD characteristic function.
void testCharacteristicFunction() {
  CHUI_REQUIRE(closeRel(nu::fermiDiracCharacteristicFunction(0.0), 1.0, 1.0e-12));
  const double x2_mean = 22.5 * k::k_zeta_5 / (1.5 * k::k_zeta_3);
  const double y = 0.01;
  CHUI_REQUIRE(std::abs(nu::fermiDiracCharacteristicFunction(y) - (1.0 - y * y * x2_mean / 6.0)) < 1.0e-8);
  for (double yy = 0.05; yy < 50.0; yy *= 1.21) {
    CHUI_REQUIRE(std::abs(nu::fermiDiracCharacteristicFunction(yy) - nu::fermiDiracCharacteristicFunctionQuadrature(yy)) <
                 2.0e-7);
  }
  // Asymptotic branch matches quadrature beyond the table.
  CHUI_REQUIRE(closeRel(nu::fermiDiracCharacteristicFunction(60.0), nu::fermiDiracCharacteristicFunctionQuadrature(60.0),
                        1.0e-5));
}

// EXT-N03-T02: superconformal time in Einstein-de Sitter, s = -2 a^-1/2 / H0.
void testSuperconformalTime() {
  c::CosmologyBackgroundConfig eds;
  eds.hubble_param = 0.7;
  eds.omega_matter = 1.0;
  eds.omega_lambda = 0.0;
  const c::LambdaCdmBackground bg(eds);
  const double expected = 2.0 * (1.0 / std::sqrt(0.02) - 1.0 / std::sqrt(0.5)) / bg.hubble0Si();
  CHUI_REQUIRE(closeRel(nu::superconformalTimeIntervalSi(bg, 0.02, 0.5), expected, 1.0e-12));
}

// EXT-N03-T03: transfer ratio multiplier M = 1 + (rho_nu/rho_cb) R.
void testTransferRatio() {
  const c::LambdaCdmBackground bg(neutrinoBackground(1.0));
  nu::NeutrinoSourceSettings settings;
  settings.method = c::NeutrinoMethodConfig::kTransferRatio;
  settings.k_bins = 16;
  settings.k_min_h_mpc = 1.0e-3;
  settings.k_max_h_mpc = 1.0;
  settings.hubble_param = 0.674;
  settings.omega_baryon = 0.05;
  nu::NeutrinoSpectralSource source(bg, syntheticSet(bg, 49.0, 0.0, 0.25, 0.0), settings, nullptr, nullptr);
  source.beginSolve(0.5);
  CHUI_REQUIRE(!source.needsObservation());
  source.finalizeObservation({}, {});
  const double f = source.densityFraction(0.5);
  CHUI_REQUIRE(closeRel(source.multiplier(0.674 * 0.1), 1.0 + f * 0.25, 1.0e-12));
  bool threw = false;
  try {
    source.requireWavenumberDomain(0.674 * 1.0e-4, 0.674 * 0.5);  // below k_min
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

// EXT-N03-T04: cold limit. Neutrinos and cold matter with identical initial
// density and velocity obey the same ODE d^2 delta/ds^2 = 4 pi G a^4 rho_m delta
// when free streaming is negligible (k u0 tau << 1, heavy 1 eV eigenstate,
// k <= 1e-3 h/Mpc), so R must stay 1. The cold history comes from an
// independent RK4 integration in s of (a, delta, delta').
void testLinearResponseColdLimit() {
  const c::LambdaCdmBackground bg(neutrinoBackground(1.0));
  const double rho_crit0 = bg.criticalDensity0Si();
  const auto rho_m = [&](double a) {
    return (bg.coldMatterDensityParameter0() / (a * a * a) + bg.totalMassiveNeutrinoDensityParameter(a)) * rho_crit0;
  };
  const double a_i = 1.0 / 50.0;
  // Initial data: delta = 1, d delta/ds = a^2 H delta (any data works: both species share it).
  double a = a_i;
  double delta = 1.0;
  double ddelta = a * a * bg.hubbleSi(a);
  const double theta_s_over_delta = -ddelta;                      // continuity: d delta/ds = -theta_s
  const double theta_class = theta_s_over_delta / (a_i * k::k_speed_of_light_si / k::k_megaparsec_si);
  nu::NeutrinoSourceSettings settings;
  settings.method = c::NeutrinoMethodConfig::kLinearResponse;
  settings.k_bins = 8;
  settings.history_capacity = 2048;
  settings.k_min_h_mpc = 1.0e-4;
  settings.k_max_h_mpc = 1.0e-3;
  settings.hubble_param = 0.674;
  settings.omega_baryon = 0.05;
  nu::NeutrinoSpectralSource source(bg, syntheticSet(bg, 49.0, 9.0, 1.0, theta_class), settings, nullptr, nullptr);
  const double total_s = nu::superconformalTimeIntervalSi(bg, a_i, 0.1);
  constexpr int steps = 800;
  const double ds = total_s / steps;
  const auto rhs = [&](double aa, double dd, double dp, double out[3]) {
    out[0] = aa * aa * aa * bg.hubbleSi(aa);                         // da/ds
    out[1] = dp;                                                     // d delta/ds
    out[2] = 4.0 * k::k_pi * k::k_newton_g_si * std::pow(aa, 4) * rho_m(aa) * dd;
  };
  const std::vector<double> modes(8, 1.0);
  for (int step = 0; step <= steps; ++step) {
    if (step % 4 == 0) {  // every 4th RK4 step is a "PM solve"
      source.beginSolve(a);
      const std::vector<double> power(8, delta * delta);
      source.finalizeObservation(power, modes);
    }
    double k1[3], k2[3], k3[3], k4[3];
    rhs(a, delta, ddelta, k1);
    rhs(a + 0.5 * ds * k1[0], delta + 0.5 * ds * k1[1], ddelta + 0.5 * ds * k1[2], k2);
    rhs(a + 0.5 * ds * k2[0], delta + 0.5 * ds * k2[1], ddelta + 0.5 * ds * k2[2], k3);
    rhs(a + ds * k3[0], delta + ds * k3[1], ddelta + ds * k3[2], k4);
    a += ds / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
    delta += ds / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
    ddelta += ds / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);
  }
  for (const double ratio : source.responseRatio()) CHUI_REQUIRE(std::abs(ratio - 1.0) < 2.0e-3);
  // EXT-N03-T05: persistence round trip and deterministic replay.
  const c::ModuleSidecarBlock block = source.serialize();
  nu::NeutrinoSpectralSource restored(bg, syntheticSet(bg, 49.0, 9.0, 1.0, theta_class), settings, nullptr, nullptr);
  restored.restore(block);
  CHUI_REQUIRE(restored.historySize() == source.historySize());
  CHUI_REQUIRE(restored.fingerprint() == source.fingerprint());
  const std::size_t before = source.historySize();
  source.beginSolve(0.05);  // earlier epoch: later entries are discarded
  CHUI_REQUIRE(source.historySize() < before);
}

// EXT-N03-T06: bounded history refuses instead of dropping entries.
void testHistoryCapacity() {
  const c::LambdaCdmBackground bg(neutrinoBackground(1.0));
  nu::NeutrinoSourceSettings settings;
  settings.method = c::NeutrinoMethodConfig::kLinearResponse;
  settings.k_bins = 8;
  settings.history_capacity = 16;
  settings.k_min_h_mpc = 1.0e-4;
  settings.k_max_h_mpc = 1.0e-3;
  settings.hubble_param = 0.674;
  settings.omega_baryon = 0.05;
  nu::NeutrinoSpectralSource source(bg, syntheticSet(bg, 49.0, 9.0, 1.0, 0.0), settings, nullptr, nullptr);
  const std::vector<double> ones(8, 1.0);
  bool threw = false;
  try {
    for (int i = 0; i < 17; ++i) {
      source.beginSolve(0.02 + 0.001 * i);
      source.finalizeObservation(ones, ones);
    }
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw && source.historySize() == 16U);
}

// Pass 04 (M04-01): transactional rollback of the live linear-response
// history. A rejected interval containing PM refreshes must leave the live
// history, current ratio/multiplier and fingerprint exactly equal to the
// opening accepted state; the retried (shortened) trajectory must equal a
// fresh direct run with the same accepted refresh epochs. Covers first use
// (empty history), repeated retries, and the replay truncation of an accepted
// node at the opening scale factor. NOT_RUN (source-only pass 04).
void testLinearResponseTransactionRollback() {
  const c::LambdaCdmBackground bg(neutrinoBackground(1.0));
  nu::NeutrinoSourceSettings settings;
  settings.method = c::NeutrinoMethodConfig::kLinearResponse;
  settings.k_bins = 8;
  settings.history_capacity = 64;
  settings.k_min_h_mpc = 1.0e-4;
  settings.k_max_h_mpc = 1.0e-3;
  settings.hubble_param = 0.674;
  settings.omega_baryon = 0.05;
  const auto make = [&]() {
    return nu::NeutrinoSpectralSource(bg, syntheticSet(bg, 49.0, 9.0, 1.0, 0.0), settings, nullptr, nullptr);
  };
  const std::vector<double> modes(8, 1.0);
  const auto refresh = [&](nu::NeutrinoSpectralSource& source, double a, double amplitude) {
    source.beginSolve(a);
    const std::vector<double> power(8, amplitude * amplitude);
    source.finalizeObservation(power, modes);
  };
  const auto same_state = [](const nu::NeutrinoSpectralSource& lhs, const nu::NeutrinoSpectralSource& rhs) {
    CHUI_REQUIRE(lhs.historySize() == rhs.historySize());
    CHUI_REQUIRE(lhs.fingerprint() == rhs.fingerprint());
    const c::ModuleSidecarBlock left = lhs.serialize();
    const c::ModuleSidecarBlock right = rhs.serialize();
    CHUI_REQUIRE(left.payload == right.payload);  // bitwise: history, ratios, initial rows
  };

  // (1) First use inside a rejected attempt: rollback restores the empty
  //     history (including the initial-epoch status).
  {
    nu::NeutrinoSpectralSource source = make();
    const nu::NeutrinoSpectralSource pristine = make();
    source.beginTransaction();
    refresh(source, 0.020, 1.0);
    refresh(source, 0.021, 1.05);
    CHUI_REQUIRE(source.historySize() == 2U);
    source.rollbackTransaction();
    CHUI_REQUIRE(!source.transactionOpen());
    same_state(source, pristine);
  }
  // (2) Rejected attempt with refreshes AFTER accepted history, retried twice
  //     with shorter intervals; compare with a direct run of the accepted
  //     epochs only.
  {
    nu::NeutrinoSpectralSource source = make();
    nu::NeutrinoSpectralSource direct = make();
    refresh(source, 0.020, 1.0);
    refresh(direct, 0.020, 1.0);
    for (int retry = 0; retry < 2; ++retry) {
      source.beginTransaction();
      refresh(source, 0.026 - 0.002 * retry, 1.3);  // rejected epochs beyond the retry
      source.rollbackTransaction();
    }
    source.beginTransaction();
    refresh(source, 0.022, 1.1);
    source.commitTransaction();
    refresh(direct, 0.022, 1.1);
    same_state(source, direct);
    // Subsequent accepted solves stay identical (no rejected node survived).
    refresh(source, 0.025, 1.2);
    refresh(direct, 0.025, 1.2);
    same_state(source, direct);
  }
  // (3) A replay at the opening scale factor pops an ACCEPTED node; a rollback
  //     must restore it exactly.
  {
    nu::NeutrinoSpectralSource source = make();
    refresh(source, 0.020, 1.0);
    refresh(source, 0.022, 1.1);
    const c::ModuleSidecarBlock opening = source.serialize();
    nu::NeutrinoSpectralSource reference = make();
    reference.restore(opening);
    source.beginTransaction();
    refresh(source, 0.022, 1.4);  // truncates the accepted node at a = 0.022
    refresh(source, 0.023, 1.5);
    source.rollbackTransaction();
    same_state(source, reference);
  }
  // (4) Commit keeps the attempt; a later rollback of a NEW transaction does
  //     not touch committed nodes.
  {
    nu::NeutrinoSpectralSource source = make();
    source.beginTransaction();
    refresh(source, 0.020, 1.0);
    source.commitTransaction();
    const c::ModuleSidecarBlock committed = source.serialize();
    source.beginTransaction();
    refresh(source, 0.030, 2.0);
    source.rollbackTransaction();
    CHUI_REQUIRE(source.serialize().payload == committed.payload);
  }
}

// EXT-N04-T01: delta-f weights, signed-mass refusal, NR drift, sampling moments.
void testParticleKernels() {
  const nu::NeutrinoSpeciesScales scales{1.0, c::instantaneousDecouplingNeutrinoTemperatureRatio() * 2.7255};
  nu::NeutrinoMarker marker = nu::sampleThermalMarker(scales, 1234U, 7U, {0.0, 0.0, 0.0}, 1.0);
  CHUI_REQUIRE(std::abs(nu::deltaFWeight(marker, scales)) < 1.0e-12);  // sampled from f_ref
  // A kick that lowers |q| raises f_ref(q) above f_i: weight becomes negative (not clipped).
  const double q = std::hypot(marker.momentum_si[0], marker.momentum_si[1], marker.momentum_si[2]);
  for (std::size_t i = 0; i < 3U; ++i) marker.momentum_si[i] *= 0.5;
  CHUI_REQUIRE(nu::deltaFWeight(marker, scales) < 0.0);
  const double signed_mass = nu::deltaFGravitationalMass(marker, scales, 0.5);
  bool threw = false;
  try {
    const std::array<double, 2> masses{1.0, signed_mass};
    nu::requireNonNegativeTreeMasses(masses);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  // Nonrelativistic drift: dx = (q/m) int dt/a^2 = (q/m) Delta s for a heavy marker at late times.
  c::CosmologyBackgroundConfig eds;
  eds.hubble_param = 0.7;
  eds.omega_matter = 1.0;
  eds.omega_lambda = 0.0;
  const c::LambdaCdmBackground bg(eds);
  const nu::NeutrinoSpeciesScales heavy{1.0e3, 1.0};
  nu::NeutrinoMarker slow;
  slow.momentum_si = {1.0e-3 * heavy.momentumScaleSi(), 0.0, 0.0};
  nu::driftMarker(slow, heavy, bg, 0.5, 1.0);
  const double expected = slow.momentum_si[0] / heavy.massKg() * nu::superconformalTimeIntervalSi(bg, 0.5, 1.0);
  CHUI_REQUIRE(closeRel(slow.position_m[0], expected, 1.0e-9));
  static_cast<void>(q);
  // Sample mean of x = |q| c / (k T) -> 7 pi^4 / 120 / (1.5 zeta 3) = 3.1514 (6 sigma tolerance).
  double sum = 0.0;
  constexpr int n = 100000;
  for (int i = 0; i < n; ++i) sum += nu::sampleFermiDiracMomentumMagnitude(nu::neutrinoCounterUniform(99U, i, 0U));
  CHUI_REQUIRE(std::abs(sum / n - c::fermiDiracMeanMomentumOverTemperature()) < 0.035);
}

// EXT-N05-T01: a single cosine mode on the mesh and noise normalizations.
void testSpeciesSpectra() {
  const std::size_t n = 16;
  nu::SpeciesMesh mesh;
  mesh.n = n;
  mesh.box_size = 100.0;
  mesh.delta.assign(n * n * n, 0.0);
  const double amplitude = 0.1;
  for (std::size_t ix = 0; ix < n; ++ix) {
    for (std::size_t iy = 0; iy < n; ++iy) {
      for (std::size_t iz = 0; iz < n; ++iz) {
        mesh.delta[(ix * n + iy) * n + iz] = amplitude * std::cos(2.0 * k::k_pi * static_cast<double>(ix) / n);
      }
    }
  }
  const auto bins = nu::crossPowerSpectrum(mesh, mesh, 12, false);
  // 12 bins of width 0.625 k_f: bin 0 holds only the 6 fundamental modes; two carry V A^2 / 4 each.
  const double volume = 1.0e6;
  CHUI_REQUIRE(bins[0].modes == 6U);
  CHUI_REQUIRE(closeRel(bins[0].power, 2.0 * volume * amplitude * amplitude / 4.0 / 6.0, 1.0e-12));
  const std::vector<double> equal(1000, 2.0);
  CHUI_REQUIRE(closeRel(nu::shotNoisePower(equal, volume), volume / 1000.0, 1.0e-14));
  const auto spectra = nu::crossSpectraFromRatio(bins, std::vector<double>(12, 0.5), 0.1);
  CHUI_REQUIRE(closeRel(spectra.total[0].power, (0.9 + 0.1 * 0.5) * (0.9 + 0.1 * 0.5) * bins[0].power, 1.0e-12));
}

}  // namespace

int main() {
  testCharacteristicFunction();
  testSuperconformalTime();
  testTransferRatio();
  testLinearResponseColdLimit();
  testHistoryCapacity();
  testLinearResponseTransactionRollback();
  testParticleKernels();
  testSpeciesSpectra();
  return 0;
}
