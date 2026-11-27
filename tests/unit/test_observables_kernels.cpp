// CHUI-EXT-PHYS-REF-v1 gates VL/VC/VRA (L01-L04, C01-C04, RA01-RA06):
// independent analytic expectations for line optical depths, equivalent
// widths, damping wings and velocity-averaging conservation; 21 cm ordered
// transfer and the radio-astronomy column relation; spin-temperature limits;
// uniform-slab SZ, kSZ sign and the tSZ null; ISW static/linear potentials;
// single-plane multi-plane = Born; flat-sky realization variance; free-free
// Kirchhoff, RM constant, Faraday screen, singular-K Stokes steps; synchrotron
// kernel normalization and kernel-vs-closed-form agreement; resonant Ly-a
// sampler moments and photon budget; FITS header blocking. NOT RUN in pass 01.
#include "../support/test_require.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/lightcone/ray_geometry.hpp"
#include "cosmosim/observables/cmb_secondaries.hpp"
#include "cosmosim/observables/fits_writer.hpp"
#include "cosmosim/observables/line_transfer.hpp"
#include "cosmosim/observables/radio_transfer.hpp"
#include "cosmosim/observables/resonant_lyman_alpha.hpp"
#include "cosmosim/radiation/radiation_angular.hpp"
#include "cosmosim/spectral/atomic_data.hpp"

namespace o = cosmosim::observables;
namespace k = cosmosim::core::constants;
namespace s = cosmosim::spectral;

namespace {

bool closeRel(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

o::LineSegment slab(double n, double path, double t) {
  o::LineSegment seg;
  seg.absorber_density_m3 = n;
  seg.path_proper_m = path;
  seg.temperature_k = t;
  return seg;
}

// EXT-L01-T03: Gaussian core depth, thin-line equivalent width, Lorentz wing, velocity-integral conservation.
void testLymanAlphaForest() {
  const s::AtomicLine& lya = s::hydrogenLymanAlpha();
  const double sigma = s::lineIntegratedCrossSectionSi(lya);
  const double t = 1.0e4;
  const double b = std::sqrt(2.0 * k::k_boltzmann_si * t / lya.atom_mass_kg);
  o::VelocityAxis axis{-2.0e5, 2.0e5, 4096, false};
  // Thin column N = 1e15 m^-2 (1e11 cm^-2) in a 1 m segment.
  const std::vector<o::LineSegment> thin{slab(1.0e15, 1.0, t)};
  const auto tau_d = o::lymanAlphaOpticalDepth(thin, axis, o::LineProfileModel::kDopplerOnly);
  const double tau0 = 1.0e15 * sigma * lya.rest_wavelength_m / (std::sqrt(k::k_pi) * b);
  CHUI_REQUIRE(closeRel(*std::max_element(tau_d.begin(), tau_d.end()), tau0, 2.0e-3));
  const auto stats = o::forestStatistics(tau_d, axis, 10, 0.0, thin);
  const double ew_thin = sigma * lya.rest_wavelength_m * lya.rest_wavelength_m * 1.0e15 / k::k_speed_of_light_si;
  CHUI_REQUIRE(closeRel(stats.equivalent_width_m, ew_thin, 1.0e-2));
  CHUI_REQUIRE(closeRel(stats.neutral_column_m2, 1.0e15, 1.0e-14));
  // Damping wing of a DLA (N = 1e25 m^-2) at 1000 km/s: tau = N sigma lambda^2 A / (4 pi^2 dv^2).
  o::VelocityAxis wing{9.99e5, 1.001e6, 1, false};
  const std::vector<o::LineSegment> dla{slab(1.0e25, 1.0, t)};
  const double dv = 1.0e6;
  const double expected = 1.0e25 * sigma * lya.rest_wavelength_m * lya.rest_wavelength_m * lya.einstein_a_s /
                          (4.0 * k::k_pi * k::k_pi * dv * dv);
  CHUI_REQUIRE(closeRel(o::lymanAlphaOpticalDepth(dla, wing, o::LineProfileModel::kVoigt)[0], expected, 1.0e-2));
  CHUI_REQUIRE(o::lymanAlphaOpticalDepth(dla, wing, o::LineProfileModel::kDopplerOnly)[0] < 1.0e-30);
  // A segment spread over 50 km/s of Hubble flow keeps int tau du = N sigma lambda0.
  o::LineSegment spread = slab(1.0e15, 1.0, t);
  spread.velocity_begin_m_s = -2.5e4;
  spread.velocity_end_m_s = 2.5e4;
  const auto tau_s = o::lymanAlphaOpticalDepth({spread}, axis, o::LineProfileModel::kDopplerOnly);
  double integral = 0.0;
  for (const double v : tau_s) integral += v * axis.width();
  CHUI_REQUIRE(closeRel(integral, 1.0e15 * sigma * lya.rest_wavelength_m, 1.0e-3));
}

// EXT-RA01-T02: ordered 21 cm transfer, column relation, spin limits.
void testHi21cm() {
  o::VelocityAxis axis{-1.0e5, 1.0e5, 2048, false};
  o::LineSegment cold = slab(1.0e6, 3.0e19, 100.0);  // n = 1 cm^-3, l = 1 kpc
  cold.spin_temperature_k = 100.0;
  const auto tb = o::hiBrightnessTemperature({cold}, axis, o::HiTransferModel::kFullOpticalDepth, 0.0);
  const auto thin = o::hiBrightnessTemperature({cold}, axis, o::HiTransferModel::kOpticallyThin, 0.0);
  // Recover tau(v) = -ln(1 - T_b/T_s); N_HI [cm^-2] = 1.823e18 T_s int tau dv [K km/s] (radio-astronomy
  // relation; 1% tolerance covers the A_10 literature spread).
  double int_tau = 0.0;
  for (std::size_t c = 0; c < tb.size(); ++c) int_tau += -std::log1p(-tb[c] / 100.0) * axis.width() / 1.0e3;
  CHUI_REQUIRE(closeRel(1.0e6 * 3.0e19 * 1.0e-4, 1.823e18 * 100.0 * int_tau, 1.0e-2));
  CHUI_REQUIRE(thin[1024] >= tb[1024]);  // the thin shortcut overestimates a saturated line
  // Ordered two-layer solution: far layer (T_s = 50) behind near layer (T_s = 200).
  o::LineSegment near = cold;
  near.spin_temperature_k = 200.0;
  o::LineSegment far = cold;
  far.spin_temperature_k = 50.0;
  const auto two = o::hiBrightnessTemperature({near, far}, axis, o::HiTransferModel::kFullOpticalDepth, 2.7);
  const auto t_near = o::hiBrightnessTemperature({near}, axis, o::HiTransferModel::kFullOpticalDepth, 0.0);
  const double tau_near = -std::log1p(-t_near[1024] / 200.0);
  const auto t_far_tau = o::hiBrightnessTemperature({far}, axis, o::HiTransferModel::kFullOpticalDepth, 0.0);
  const double tau_far = -std::log1p(-t_far_tau[1024] / 50.0);
  const double behind = 2.7 * std::exp(-tau_far) + 50.0 * (1.0 - std::exp(-tau_far));
  CHUI_REQUIRE(closeRel(two[1024], behind * std::exp(-tau_near) + 200.0 * (1.0 - std::exp(-tau_near)), 1.0e-12));
  // Spin-temperature limits.
  CHUI_REQUIRE(closeRel(o::spinTemperature({100.0, 30.0, 1.0e12, 0.0, 0.0}), 100.0, 1.0e-9));
  CHUI_REQUIRE(closeRel(o::spinTemperature({100.0, 30.0, 0.0, 0.0, 0.0}), 30.0, 1.0e-14));
}

// EXT-C03-T01: SZ slab, kSZ sign, tSZ null near 217 GHz, ISW.
void testSz() {
  o::SzSegment seg{3.0857e22, 1.0e3, 1.0e8, 3.0e5};
  const auto sz = o::integrateSz({seg});
  const double expected_y = k::k_thomson_cross_section_si * 1.0e3 * k::k_boltzmann_si * 1.0e8 * 3.0857e22 /
                            (k::k_electron_mass_si * k::k_speed_of_light_si * k::k_speed_of_light_si);
  CHUI_REQUIRE(closeRel(sz.compton_y, expected_y, 1.0e-14));
  CHUI_REQUIRE(sz.ksz_delta_t_over_t < 0.0);  // receding gas -> decrement
  double lo = 150.0e9;
  double hi = 300.0e9;
  for (int i = 0; i < 100; ++i) {
    const double mid = 0.5 * (lo + hi);
    (o::tszSpectralFactor(mid, 2.7255) < 0.0 ? lo : hi) = mid;
  }
  CHUI_REQUIRE(std::abs(0.5 * (lo + hi) - 217.5e9) < 1.5e9);
  cosmosim::core::CosmologyBackgroundConfig eds;
  eds.hubble_param = 0.7;
  eds.omega_matter = 1.0;
  eds.omega_lambda = 0.0;
  const cosmosim::core::LambdaCdmBackground bg(eds);
  CHUI_REQUIRE(o::integrateIsw({{1.0e24, 0.5, 1.0e4, 1.0e4}}, 0.4, 0.6, bg) == 0.0);
  const double isw = o::integrateIsw({{1.0e24, 0.5, 1.0e4, 2.0e4}}, 0.4, 0.6, bg);
  const double expected_isw = 2.0 / (k::k_speed_of_light_si * k::k_speed_of_light_si) * (1.0e4 / 0.2) * 0.5 *
                              bg.hubbleSi(0.5) * 0.5 * 1.0e24 / k::k_speed_of_light_si;
  CHUI_REQUIRE(closeRel(isw, expected_isw, 1.0e-12));
}

// EXT-C02-T01: one lens plane: multi-plane ray shooting equals Born (weak sinusoid).
void testLensing() {
  o::LensPlane plane;
  plane.chi_m = 1.0e25;
  plane.n = 256;  // 64 cells per wavelength keeps bilinear interpolation error ~1e-3
  plane.side_m = 1.0e23;
  plane.sigma_m.resize(256 * 256);
  for (std::size_t iy = 0; iy < 256; ++iy) {
    for (std::size_t ix = 0; ix < 256; ++ix) {
      plane.sigma_m[iy * 256 + ix] = 1.0e19 * std::cos(2.0 * k::k_pi * 4.0 * (ix + 0.5) / 256.0);
    }
  }
  const double h0 = 0.7 * k::k_hubble_100_km_s_mpc_si;
  const auto lens = o::rayShootLensPlanes({plane}, 2.0e25, h0, 0.3, 128, 1.0e23 / 1.0e25);
  double diff = 0.0;
  double norm = 0.0;
  for (std::size_t i = 0; i < lens.kappa_born.size(); ++i) {
    const std::size_t rx = i % 128;
    const std::size_t ry = i / 128;
    if (rx == 0 || ry == 0 || rx == 127 || ry == 127) continue;  // one-sided Jacobian at the edges
    diff += std::pow(lens.kappa_multiplane[i] - lens.kappa_born[i], 2);
    norm += lens.kappa_born[i] * lens.kappa_born[i];
  }
  CHUI_REQUIRE(std::sqrt(diff / norm) < 0.02);
}

// EXT-C01-T01: flat-sky realization of a white spectrum: pixel variance = C (N^2 - 1)/A.
void testFlatSkyRealization() {
  o::CmbSpectraTable table;
  table.ell = {1.0, 1.0e6};
  table.cl_tt = {2.0, 2.0};
  table.cl_ee = {0.0, 0.0};
  table.cl_te = {0.0, 0.0};
  const std::size_t n = 128;
  const double side = 0.05;
  const auto map = o::gaussianRealizationFlatSky(table, n, side, 42U);
  double mean = 0.0;
  for (const double v : map.values) mean += v / (n * n);
  double var = 0.0;
  for (const double v : map.values) var += (v - mean) * (v - mean) / (n * n);
  CHUI_REQUIRE(closeRel(var, 2.0 * (n * n - 1.0) / (side * side), 0.05));
}

// EXT-RA03/RA04-T01: Kirchhoff, RM constant, Faraday screen, singular-K Stokes.
void testRadio() {
  const double nu = 1.4e9;
  CHUI_REQUIRE(closeRel(o::freeFreeEmissivity(nu, 8.0e3, 1.0e6, 1.0e6) / o::freeFreeAbsorption(nu, 8.0e3, 1.0e6, 1.0e6),
                        o::planckBnu(nu, 8.0e3), 1.0e-12));
  // 0.812 rad m^-2 per (cm^-3 microgauss pc).
  CHUI_REQUIRE(closeRel(o::rotationMeasureCoefficient() * 1.0e6 * 1.0e-10 * 3.0857e16, 0.812, 2.0e-3));
  const double rho = 3.0;
  const auto rot = o::propagationMatrix(0, 0, 0, 0, 0, 0, rho);
  const auto out = o::stokesStep({1.0, 1.0, 0.0, 0.0}, {0, 0, 0, 0}, rot, 0.7);
  CHUI_REQUIRE(std::abs(out[1] - std::cos(rho * 0.7)) < 1.0e-12 && std::abs(out[2] - std::sin(rho * 0.7)) < 1.0e-12);
  const auto emit = o::stokesStep({0, 0, 0, 0}, {2.0, 0, 0, 0}, o::propagationMatrix(0, 0, 0, 0, 0, 0, 0), 5.0);
  CHUI_REQUIRE(closeRel(emit[0], 10.0, 1.0e-13));  // singular K: I = j L
  const auto absorb = o::stokesStep({0, 0, 0, 0}, {2.0, 0, 0, 0}, o::propagationMatrix(0.5, 0, 0, 0, 0, 0, 0), 3.0);
  CHUI_REQUIRE(closeRel(absorb[0], 4.0 * (1.0 - std::exp(-1.5)), 1.0e-12));
  // DM of a uniform screen.
  const auto dm = o::integrateDispersion({{3.0857e22, 1.0e4, 0.0, 1.0}});
  CHUI_REQUIRE(closeRel(dm.dm_m2, 1.0e4 * 3.0857e22 / 2.0, 1.0e-14));
}

// EXT-RA03-T02: synchrotron kernel normalization and kernel vs closed forms.
void testSynchrotron() {
  // int_0^inf F(x) dx = 8 pi / (9 sqrt 3); trapezoid in ln x.
  double integral = 0.0;
  for (double lx = std::log(1.0e-7); lx < std::log(60.0); lx += 1.0e-3) {
    const double x0 = std::exp(lx);
    const double x1 = std::exp(lx + 1.0e-3);
    integral += 0.5 * (x0 * o::synchrotronKernelF(x0) + x1 * o::synchrotronKernelF(x1)) * 1.0e-3;
  }
  CHUI_REQUIRE(closeRel(integral, 8.0 * k::k_pi / (9.0 * std::sqrt(3.0)), 2.0e-3));
  const double p = 2.5;
  const double c_m3 = 1.0e2;
  o::ElectronSpectrumSample electrons;
  for (double lg = 1.0; lg <= 6.0; lg += 0.01) {
    const double g = std::pow(10.0, lg);
    electrons.gamma.push_back(g);
    electrons.number_per_m3.push_back(c_m3 * std::pow(g, -p));
  }
  const double b = 1.0e-9;
  const auto kernel = o::synchrotronFromSpectrum(electrons, 1.4e9, b, 1.0);
  CHUI_REQUIRE(closeRel(kernel.j_i, o::synchrotronPowerLawEmissivity(c_m3, p, 1.4e9, b, 1.0), 2.0e-2));
  CHUI_REQUIRE(closeRel(kernel.alpha_i, o::synchrotronPowerLawAbsorption(c_m3, p, 1.4e9, b, 1.0), 3.0e-2));
  CHUI_REQUIRE(closeRel(kernel.linear_fraction, o::synchrotronPowerLawPolarization(p), 1.0e-2));
  // Spectral slopes: j ~ nu^-(p-1)/2, alpha ~ nu^-(p+4)/2.
  const double j2 = o::synchrotronPowerLawEmissivity(c_m3, p, 2.8e9, b, 1.0);
  CHUI_REQUIRE(closeRel(j2 / o::synchrotronPowerLawEmissivity(c_m3, p, 1.4e9, b, 1.0), std::pow(2.0, -(p - 1.0) / 2.0), 1.0e-12));
}

// EXT-L02-T01: exact parallel-velocity sampler (moments vs quadrature) and photon budget.
void testResonantLymanAlpha() {
  for (const double x : {0.0, 3.0}) {
    const double a = 0.01;
    // Target moments by direct quadrature of exp(-u^2)/((x-u)^2 + a^2).
    double z = 0.0;
    double m1 = 0.0;
    double m2 = 0.0;
    for (double u = -8.0; u <= 8.0; u += 1.0e-5) {
      const double f = std::exp(-u * u) / ((x - u) * (x - u) + a * a);
      z += f;
      m1 += u * f;
      m2 += u * u * f;
    }
    m1 /= z;
    m2 /= z;
    std::uint64_t counter = 0;
    const auto uniform = [&]() {
      const std::uint64_t bits = cosmosim::radiation::radiationCounterRandom(11U, 3U, counter++);
      return (static_cast<double>(bits >> 11U) + 0.5) * (1.0 / 9007199254740992.0);
    };
    double s1 = 0.0;
    double s2 = 0.0;
    const int n = 200000;
    for (int i = 0; i < n; ++i) {
      const double u = o::sampleParallelAtomVelocity(x, a, uniform);
      s1 += u;
      s2 += u * u;
    }
    const double mean = s1 / n;
    const double var = s2 / n - mean * mean;
    const double se = std::sqrt(std::max(m2 - m1 * m1, 1.0e-12) / n);
    CHUI_REQUIRE(std::abs(mean - m1) < 5.0 * se);
    CHUI_REQUIRE(closeRel(var, m2 - m1 * m1, 0.03));
  }
  // Photon budget in a small optically thick static cube with dust.
  const auto geometry = cosmosim::lightcone::LeafGeometry::uniformGrid({4, 4, 4}, {0, 0, 0}, {1, 1, 1},
                                                                     cosmosim::lightcone::RayBoundary::kOpen);
  const std::size_t cells = geometry.cellCount();
  std::vector<double> n_hi(cells, 1.0), temp(cells, 1.0e4), zero(cells, 0.0), dust(cells, 1.0e-16);
  o::ResonantMedium medium{n_hi, temp, zero, zero, zero, dust, 1.0e16, 0.0};
  o::ResonantSettings settings;
  settings.packets = 2000;
  const auto result = o::runResonantLymanAlpha(geometry, medium, {0.5, 0.5, 0.5}, settings);
  CHUI_REQUIRE(std::abs(result.escape_fraction + result.dust_absorbed_fraction + result.truncated_fraction - 1.0) < 1.0e-12);
  CHUI_REQUIRE(result.mean_scatterings > 1.0);
}

// EXT-O01-T01: FITS header cards and 2880-byte blocking.
void testFitsHeader() {
  o::FitsImage image;
  image.axes.push_back({4, "VRAD", "m/s", 1.0, 0.0, 1.0e3});
  image.data.assign(4, 1.0);
  image.history.push_back(std::string(200, 'a'));
  const std::string header = o::fitsHeader(image);
  CHUI_REQUIRE(header.size() % 2880U == 0U);
  CHUI_REQUIRE(header.rfind("SIMPLE  =                    T", 0) == 0U);
  CHUI_REQUIRE(header.find("END     ") != std::string::npos);
}

}  // namespace

int main() {
  testLymanAlphaForest();
  testHi21cm();
  testSz();
  testLensing();
  testFlatSkyRealization();
  testRadio();
  testSynchrotron();
  testResonantLymanAlpha();
  testFitsHeader();
  return 0;
}
