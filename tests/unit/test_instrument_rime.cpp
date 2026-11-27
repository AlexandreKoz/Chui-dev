// CHUI-EXT-PHYS-REF-v1 gate VO (O01-O05): point source at the phase centre,
// phase sign and amplitude of an offset source, polarization coherency
// normalization, leakage, uvw for a polar source, beam normalization,
// bandwidth smearing against the analytic sinc, noise variance and
// determinism, flags and baseline count. Expected values are analytic.
// NOT RUN in pass 01.
#include "../support/test_require.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "cosmosim/core/constants.hpp"
#include "cosmosim/instruments/interferometry.hpp"

namespace in = cosmosim::instruments;
namespace k = cosmosim::core::constants;

namespace {

bool closeAbs(double a, double b, double tol) { return std::abs(a - b) <= tol; }

std::vector<in::Station> twoStations(double east_m) {
  in::Station a;
  a.name = "A";
  in::Station b;
  b.name = "B";
  b.enu_m = {east_m, 0.0, 0.0};
  return {a, b};
}

in::ObservationSettings narrowBand(double nu, std::uint32_t subsamples) {
  in::ObservationSettings s;
  s.latitude_rad = 0.0;
  s.phase_centre_dec_rad = 0.0;
  s.hour_angle_start_rad = 0.0;
  s.integration_s = 1.0e-6;  // negligible Earth rotation
  s.time_steps = 1;
  s.frequencies = cosmosim::spectral::FrequencyAxis({nu * (1.0 - 1.0e-9), nu * (1.0 + 1.0e-9)}, subsamples);
  return s;
}

std::vector<in::VisibilityBatch> run(const std::vector<in::Station>& st, const std::vector<in::SkyComponent>& sky,
                                     const in::ObservationSettings& s, const in::FlagFunction& flags = {}) {
  std::vector<in::VisibilityBatch> out;
  (void)in::simulateVisibilitiesDirect(st, sky, s, [&out](const in::VisibilityBatch& b) { out.push_back(b); }, flags);
  return out;
}

// EXT-O05-T01: unpolarized point source at the phase centre: XX = YY = I, XY = YX = 0.
void testPhaseCentre() {
  const auto batches = run(twoStations(100.0), {{0.0, 0.0, {{2.0, 0.0, 0.0, 0.0}}}}, narrowBand(1.0e9, 1));
  const auto v = batches.at(0).data.at(0);
  CHUI_REQUIRE(closeAbs(v[0].real(), 2.0, 1e-12) && closeAbs(v[3].real(), 2.0, 1e-12));
  CHUI_REQUIRE(std::abs(v[1]) < 1e-12 && std::abs(v[2]) < 1e-12);
  // Polarized: Q = 0.5 gives XX = I + Q, YY = I - Q; U, V populate XY = U + iV.
  const auto pol = run(twoStations(100.0), {{0.0, 0.0, {{2.0, 0.5, 0.3, 0.1}}}}, narrowBand(1.0e9, 1)).at(0).data.at(0);
  CHUI_REQUIRE(closeAbs(pol[0].real(), 2.5, 1e-12) && closeAbs(pol[3].real(), 1.5, 1e-12));
  CHUI_REQUIRE(closeAbs(pol[1].real(), 0.3, 1e-12) && closeAbs(pol[1].imag(), 0.1, 1e-12));
}

// EXT-O05-T02: east-west baseline at HA = 0 and dec = 0: u = L / lambda, and an
// offset source at l picks up exp(-2 pi i u l) (w(n - 1) term included exactly).
void testPhaseSign() {
  const double nu = 1.0e9;
  const double lambda = k::k_speed_of_light_si / nu;
  const double east = 50.0;
  const auto batches = run(twoStations(east), {{0.01, 0.0, {{1.0, 0.0, 0.0, 0.0}}}}, narrowBand(nu, 1));
  const auto& b = batches.at(0);
  CHUI_REQUIRE(closeAbs(b.uvw_m[0][0], east, 1e-9) && closeAbs(b.uvw_m[0][1], 0.0, 1e-9));
  const double n = std::sqrt(1.0 - 1.0e-4);
  const double phase = -2.0 * k::k_pi * (east / lambda * 0.01 + b.uvw_m[0][2] / lambda * (n - 1.0));
  CHUI_REQUIRE(closeAbs(std::arg(b.data[0][0]), std::remainder(phase, 2.0 * k::k_pi), 1e-9));
  CHUI_REQUIRE(closeAbs(std::abs(b.data[0][0]), 1.0, 1e-12));
}

// EXT-O05-T03: bandwidth smearing of a flat-spectrum offset source equals the sinc average.
void testBandwidthSmearing() {
  const double nu = 1.0e9;
  const double width = 5.0e6;
  in::ObservationSettings s = narrowBand(nu, 32);
  s.frequencies = cosmosim::spectral::FrequencyAxis({nu - 0.5 * width, nu + 0.5 * width}, 32);
  const double east = 2000.0;
  const double l = 0.02;
  const auto v = run(twoStations(east), {{l, 0.0, {{1.0, 0.0, 0.0, 0.0}}}}, s).at(0).data.at(0)[0];
  // |<exp(-2 pi i nu tau)>| over the band = |sinc(pi width tau)|, tau = east l / c (w term ~ second order).
  const double tau = east * l / k::k_speed_of_light_si;
  const double x = k::k_pi * width * tau;
  CHUI_REQUIRE(closeAbs(std::abs(v), std::abs(std::sin(x) / x), 2e-3));
}

// EXT-O05-T04: leakage and gains enter as J_p B J_q^H.
void testJones() {
  auto st = twoStations(10.0);
  st[0].leakage_x = {0.05, 0.0};
  st[1].gain_x = {2.0, 0.0};
  const auto v = run(st, {{0.0, 0.0, {{1.0, 0.0, 0.0, 0.0}}}}, narrowBand(1.0e9, 1)).at(0).data.at(0);
  // J_p = [[1, 0.05], [0, 1]], J_q = diag(2, 1), B = I: V = J_p J_q^H = [[2, 0.05], [0, 1]].
  CHUI_REQUIRE(closeAbs(v[0].real(), 2.0, 1e-12) && closeAbs(v[1].real(), 0.05, 1e-12) && closeAbs(v[3].real(), 1.0, 1e-12));
  // Beams: unity at the centre, Gaussian half power at FWHM/2.
  in::BeamSettings g{in::BeamModel::kGaussian, 0.1, 1.0e9, 0.0};
  CHUI_REQUIRE(closeAbs(in::beamVoltage(g, 0.0, 0.0, 1.0e9), 1.0, 1e-15));
  const double half = std::sin(0.05);
  CHUI_REQUIRE(closeAbs(std::pow(in::beamVoltage(g, half, 0.0, 1.0e9), 2), 0.5, 1e-12));
  in::BeamSettings airy{in::BeamModel::kAiry, 0.0, 0.0, 25.0};
  CHUI_REQUIRE(closeAbs(in::beamVoltage(airy, 0.0, 0.0, 1.0e9), 1.0, 1e-12));
  // First Airy null at x = 3.8317 (first zero of J1).
  const double s0 = 3.8317059702075123 * k::k_speed_of_light_si / (k::k_pi * 25.0 * 1.0e9);
  CHUI_REQUIRE(std::abs(in::beamVoltage(airy, s0, 0.0, 1.0e9)) < 1e-6);
}

// EXT-O05-T05: noise variance, determinism, flags, baseline count.
void testNoiseAndFlags() {
  std::vector<in::Station> st;
  for (int i = 0; i < 12; ++i) {
    in::Station s;
    s.name = "S" + std::to_string(i);
    s.enu_m = {13.0 * i, 7.0 * i * i, 0.0};
    st.push_back(s);
  }
  CHUI_REQUIRE(in::baselineCount(12) == 66U);
  in::ObservationSettings s = narrowBand(1.0e9, 1);
  s.frequencies = cosmosim::spectral::FrequencyAxis::uniform(1.0e9, 1.001e9, 64, 1);
  s.sefd_jy = 400.0;
  s.integration_s = 10.0;
  s.time_steps = 4;
  const auto a = run(st, {}, s);
  const auto b = run(st, {}, s);
  double sum2 = 0.0;
  std::size_t count = 0;
  for (std::size_t t = 0; t < a.size(); ++t) {
    for (std::size_t i = 0; i < a[t].data.size(); ++i) {
      CHUI_REQUIRE(a[t].data[i][0] == b[t].data[i][0]);  // deterministic seeds
      sum2 += std::norm(a[t].data[i][0]) / 2.0;
      ++count;
    }
  }
  const double sigma = in::thermalNoiseSigmaJy(400.0, 1.0e6 / 64.0 * 1.0, 10.0);
  const double expected = sigma * sigma;
  CHUI_REQUIRE(std::abs(sum2 / count / expected - 1.0) < 5.0 * std::sqrt(1.0 / count));
  const auto flagged = run(st, {}, s, [](std::size_t baseline, std::size_t, std::size_t channel) {
    return baseline == 3U && channel < 8U;
  });
  CHUI_REQUIRE(flagged[0].flags[3 * 64 + 0] == 1U && flagged[0].flags[3 * 64 + 8] == 0U);
  CHUI_REQUIRE(std::abs(flagged[0].data[3 * 64 + 0][0]) == 0.0);
}

// EXT-O05-T06: a source at the celestial pole has w = baseline z and rotation-invariant |uv|.
void testPolarUvw() {
  const std::array<double, 3> l{30.0, 40.0, 5.0};
  const auto a = in::baselineUvwMetres(l, 0.0, 0.5 * k::k_pi);
  const auto b = in::baselineUvwMetres(l, 1.0, 0.5 * k::k_pi);
  CHUI_REQUIRE(closeAbs(a[2], 5.0, 1e-12) && closeAbs(b[2], 5.0, 1e-12));
  CHUI_REQUIRE(closeAbs(std::hypot(a[0], a[1]), 50.0, 1e-12) && closeAbs(std::hypot(b[0], b[1]), 50.0, 1e-12));
}

}  // namespace

int main() {
  testPhaseCentre();
  testPhaseSign();
  testBandwidthSmearing();
  testJones();
  testNoiseAndFlags();
  testPolarUvw();
  return 0;
}
