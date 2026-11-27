// CHUI-EXT-PHYS-REF-v1 N04 neutrino marker kernels (reference only; source
// only, NOT compiled or executed in pass 01).

#include "cosmosim/neutrinos/neutrino_particles.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/constants.hpp"
#include "neutrinos/internal/gauss_legendre.hpp"

namespace cosmosim::neutrinos {
namespace {

namespace k = core::constants;

constexpr std::size_t k_cdf_nodes = 2048U;
constexpr double k_cdf_x_max = 40.0;

[[nodiscard]] std::uint64_t splitMix64(std::uint64_t x) noexcept {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31U);
}

[[nodiscard]] double fd(double x) { return x > 40.0 ? std::exp(-x) / (1.0 + std::exp(-x)) : 1.0 / (std::exp(x) + 1.0); }

// Cumulative distribution of x^2 f(x) on a uniform grid (Simpson panels of
// 32 sub-steps per node interval; normalized by the last value, which equals
// (3/2) zeta(3) to < 1e-12).
[[nodiscard]] const std::vector<double>& fdCdf() {
  static const std::vector<double> cdf = [] {
    std::vector<double> values(k_cdf_nodes + 1U, 0.0);
    const double dx = k_cdf_x_max / static_cast<double>(k_cdf_nodes);
    constexpr int sub = 32;
    for (std::size_t i = 0; i < k_cdf_nodes; ++i) {
      const double x0 = dx * static_cast<double>(i);
      const double h = dx / sub;
      double s = 0.0;
      for (int j = 0; j < sub; ++j) {
        const double a = x0 + h * j;
        const double m = a + 0.5 * h;
        const double b = a + h;
        s += h / 6.0 * (a * a * fd(a) + 4.0 * m * m * fd(m) + b * b * fd(b));
      }
      values[i + 1U] = values[i] + s;
    }
    const double total = values.back();
    for (double& v : values) v /= total;
    return values;
  }();
  return cdf;
}

[[nodiscard]] double integrateDrift(const NeutrinoSpeciesScales& scales, double q, const core::LambdaCdmBackground& bg,
                                    double a0, double a1) {
  // int da / (a^2 H eps) = int d ln a / (a H eps), 16-point Gauss-Legendre
  // (computed nodes) on panels of <= 0.01 in ln a.
  const internal::GaussLegendreRule& rule = internal::gaussLegendre16();
  const double l0 = std::log(a0);
  const double l1 = std::log(a1);
  const auto panels = static_cast<std::size_t>(std::max(1.0, std::ceil(std::abs(l1 - l0) / 0.01)));
  const double width = (l1 - l0) / static_cast<double>(panels);
  double sum = 0.0;
  for (std::size_t p = 0; p < panels; ++p) {
    const double mid = l0 + width * (static_cast<double>(p) + 0.5);
    for (std::size_t i = 0; i < rule.nodes.size(); ++i) {
      const double a = std::exp(mid + 0.5 * width * rule.nodes[i]);
      sum += 0.5 * width * rule.weights[i] / (a * bg.hubbleSi(a) * comovingEnergyMass(scales, q, a));
    }
  }
  return sum;
}

}  // namespace

double neutrinoCounterUniform(std::uint64_t seed, std::uint64_t particle_id, std::uint64_t counter) noexcept {
  const std::uint64_t bits = splitMix64(splitMix64(seed ^ splitMix64(particle_id)) ^ counter);
  // 53 random mantissa bits, mapped to the open interval (0, 1).
  return (static_cast<double>(bits >> 11U) + 0.5) * (1.0 / 9007199254740992.0);
}

double sampleFermiDiracMomentumMagnitude(double u) {
  if (!(u > 0.0) || !(u < 1.0)) throw std::invalid_argument("FD sampling needs u in (0, 1)");
  const std::vector<double>& cdf = fdCdf();
  const auto it = std::upper_bound(cdf.begin(), cdf.end(), u);
  const auto hi = static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(std::distance(cdf.begin(), it), 1,
                                                                      static_cast<std::ptrdiff_t>(cdf.size() - 1U)));
  const double dx = k_cdf_x_max / static_cast<double>(k_cdf_nodes);
  const double t = (u - cdf[hi - 1U]) / (cdf[hi] - cdf[hi - 1U]);
  return dx * (static_cast<double>(hi - 1U) + t);
}

double NeutrinoSpeciesScales::momentumScaleSi() const {
  if (!(temperature_k > 0.0)) throw std::invalid_argument("neutrino temperature must be > 0");
  return k::k_boltzmann_si * temperature_k / k::k_speed_of_light_si;
}

double NeutrinoSpeciesScales::massKg() const {
  if (!(mass_ev > 0.0)) throw std::invalid_argument("neutrino mass must be > 0");
  return mass_ev * k::k_electron_volt_si / (k::k_speed_of_light_si * k::k_speed_of_light_si);
}

NeutrinoMarker sampleThermalMarker(const NeutrinoSpeciesScales& scales, std::uint64_t seed, std::uint64_t particle_id,
                                   const std::array<double, 3>& position_m, double rest_mass_particle_kg) {
  NeutrinoMarker marker;
  marker.particle_id = particle_id;
  marker.position_m = position_m;
  marker.rest_mass_particle_kg = rest_mass_particle_kg;
  const double x = sampleFermiDiracMomentumMagnitude(neutrinoCounterUniform(seed, particle_id, 0U));
  const double mu = 2.0 * neutrinoCounterUniform(seed, particle_id, 1U) - 1.0;
  const double phi = 2.0 * k::k_pi * neutrinoCounterUniform(seed, particle_id, 2U);
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  const double q = x * scales.momentumScaleSi();
  marker.momentum_si = {q * sin_theta * std::cos(phi), q * sin_theta * std::sin(phi), q * mu};
  marker.initial_reference_density = fd(x);
  return marker;
}

double comovingEnergyMass(const NeutrinoSpeciesScales& scales, double q, double a) {
  const double c = k::k_speed_of_light_si;
  const double m = scales.massKg();
  return std::sqrt(q * q * c * c + m * m * c * c * c * c * a * a) / (c * c);
}

void driftMarker(NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales, const core::LambdaCdmBackground& background,
                 double a0, double a1) {
  if (!(a0 > 0.0) || !(a1 > 0.0)) throw std::invalid_argument("drift needs a > 0");
  const double q = std::hypot(marker.momentum_si[0], marker.momentum_si[1], marker.momentum_si[2]);
  const double factor = integrateDrift(scales, q, background, a0, a1);
  for (std::size_t i = 0; i < 3U; ++i) marker.position_m[i] += marker.momentum_si[i] * factor;
}

void kickMarker(NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales, const std::array<double, 3>& grad_psi,
                const core::LambdaCdmBackground& background, double a0, double a1) {
  if (!(a0 > 0.0) || !(a1 > 0.0)) throw std::invalid_argument("kick needs a > 0");
  const double c = k::k_speed_of_light_si;
  const double q = std::hypot(marker.momentum_si[0], marker.momentum_si[1], marker.momentum_si[2]);
  // dq = -(eps/a)(1 + v^2/c^2) grad Psi dt; dt = da / (a H). With q frozen the
  // a-dependence is eps(q, a)(1 + v^2)/a^2/H; integrate with the midpoint of
  // ln a refined by 16 sub-panels (explicit first-order kick in q).
  const double l0 = std::log(a0);
  const double l1 = std::log(a1);
  constexpr int panels = 16;
  const double width = (l1 - l0) / panels;
  double factor = 0.0;
  for (int p = 0; p < panels; ++p) {
    const double a = std::exp(l0 + width * (p + 0.5));
    const double eps = comovingEnergyMass(scales, q, a);
    const double v_over_c = q / (eps * c);
    factor += width * eps * (1.0 + v_over_c * v_over_c) / (a * background.hubbleSi(a));
  }
  for (std::size_t i = 0; i < 3U; ++i) marker.momentum_si[i] -= factor * grad_psi[i];
}

double deltaFWeight(const NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales) {
  if (!(marker.initial_reference_density > 0.0)) throw std::invalid_argument("delta-f marker has no initial density");
  const double q = std::hypot(marker.momentum_si[0], marker.momentum_si[1], marker.momentum_si[2]);
  return 1.0 - fd(q / scales.momentumScaleSi()) / marker.initial_reference_density;
}

double markerEnergyMass(const NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales, double a) {
  const double q = std::hypot(marker.momentum_si[0], marker.momentum_si[1], marker.momentum_si[2]);
  const double gamma = comovingEnergyMass(scales, q, a) / (scales.massKg() * a);
  return marker.rest_mass_particle_kg * gamma;
}

double deltaFGravitationalMass(const NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales, double a) {
  return deltaFWeight(marker, scales) * markerEnergyMass(marker, scales, a);
}

void requireNonNegativeTreeMasses(std::span<const double> masses) {
  for (const double m : masses) {
    if (!(m >= 0.0)) {
      throw std::invalid_argument(
          "signed (delta-f) estimator mass on the positive-mass tree path: multipoles and opening criteria assume "
          "positive masses; use the PM-only scale-limited neutrino force model");
    }
  }
}

}  // namespace cosmosim::neutrinos
