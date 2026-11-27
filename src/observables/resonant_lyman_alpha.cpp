// CHUI-EXT-PHYS-REF-v1 L02/L03 resonant Ly-a reference Monte Carlo (source
// only, NOT compiled or executed in pass 01).

#include "cosmosim/observables/resonant_lyman_alpha.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"
#include "cosmosim/radiation/radiation_angular.hpp"
#include "cosmosim/spectral/atomic_data.hpp"

namespace cosmosim::observables {
namespace {

namespace k = core::constants;
constexpr double k_pi = 3.141592653589793238462643383279502884;

struct Rng {
  std::uint64_t seed;
  std::uint64_t stream;
  std::uint64_t counter = 0;
  double operator()() {
    const std::uint64_t bits = radiation::radiationCounterRandom(seed, stream, counter++);
    return (static_cast<double>(bits >> 11U) + 0.5) * (1.0 / 9007199254740992.0);
  }
};

[[nodiscard]] std::array<double, 3> isotropicDirection(Rng& rng) {
  const double mu = 2.0 * rng() - 1.0;
  const double phi = 2.0 * k_pi * rng();
  const double s = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  return {s * std::cos(phi), s * std::sin(phi), mu};
}

void basis(const std::array<double, 3>& d, std::array<double, 3>& e1, std::array<double, 3>& e2) {
  std::array<double, 3> t = std::abs(d[2]) < 0.9 ? std::array<double, 3>{0.0, 0.0, 1.0} : std::array<double, 3>{1.0, 0.0, 0.0};
  e1 = {d[1] * t[2] - d[2] * t[1], d[2] * t[0] - d[0] * t[2], d[0] * t[1] - d[1] * t[0]};
  const double n = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  for (double& v : e1) v /= n;
  e2 = {d[1] * e1[2] - d[2] * e1[1], d[2] * e1[0] - d[0] * e1[2], d[0] * e1[1] - d[1] * e1[0]};
}

[[nodiscard]] double dot(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

}  // namespace

ResonantResult runResonantLymanAlpha(const lightcone::LeafGeometry& geometry, const ResonantMedium& m,
                                     const std::array<double, 3>& source, const ResonantSettings& s) {
  if (geometry.boundary() != lightcone::RayBoundary::kOpen) {
    throw std::invalid_argument("resonant Ly-a reference needs an open domain (escape boundary)");
  }
  const std::size_t cells = geometry.cellCount();
  if (m.n_hi_m3.size() < cells || m.temperature_k.size() < cells || m.vx.size() < cells || m.vy.size() < cells ||
      m.vz.size() < cells || (!m.dust_absorption_per_m.empty() && m.dust_absorption_per_m.size() < cells)) {
    throw std::invalid_argument("resonant medium spans smaller than the geometry");
  }
  if (s.packets == 0U || s.spectrum_bins == 0U || !(s.spectrum_w_max > s.spectrum_w_min)) {
    throw std::invalid_argument("invalid resonant Ly-a settings");
  }
  const spectral::AtomicLine& line = spectral::hydrogenLymanAlpha();
  const double sigma_lambda = spectral::lineIntegratedCrossSectionSi(line) * line.rest_wavelength_m;
  const double photon_momentum = k::k_planck_si * line.rest_frequency_hz / k::k_speed_of_light_si;
  const double rate_per_packet = s.photon_rate_s / static_cast<double>(s.packets);
  ResonantResult result;
  result.escaped_spectrum.assign(s.spectrum_bins, 0.0);
  result.scattering_rate_per_cell.assign(cells, 0.0);
  result.momentum_rate_per_cell.assign(cells, {0.0, 0.0, 0.0});
  const auto extent = geometry.domainExtent();
  const double chunk = std::sqrt(extent[0] * extent[0] + extent[1] * extent[1] + extent[2] * extent[2]);
  std::uint64_t total_scatterings = 0;
  std::uint64_t escaped = 0;
  std::uint64_t dusted = 0;
  std::uint64_t truncated = 0;
  const auto bulk = [&](std::size_t cell, const std::array<double, 3>& r) {
    std::array<double, 3> v{m.vx[cell], m.vy[cell], m.vz[cell]};
    if (m.hubble_rate_s != 0.0) {
      for (std::size_t a = 0; a < 3U; ++a) v[a] += m.hubble_rate_s * (r[a] - source[a]) * m.length_unit_m;
    }
    return v;
  };
  for (std::uint64_t id = 0; id < s.packets; ++id) {
    Rng rng{s.seed, id};
    std::array<double, 3> pos = source;
    std::array<double, 3> dir = isotropicDirection(rng);
    double w = s.initial_w_m_s;
    std::uint64_t scatterings = 0;
    bool alive = true;
    bool escaped_now = false;
    while (alive) {
      double tau_target = -std::log(rng());
      bool interacted = false;
      while (!interacted) {
        const auto segments = geometry.traceStraightRay(pos, dir, chunk, 1U << 20U);
        double covered = 0.0;
        for (const auto& seg : segments) {
          const std::size_t c = seg.cell_row;
          const double t = m.temperature_k[c];
          const double b = std::sqrt(2.0 * k::k_boltzmann_si * t / line.atom_mass_kg);
          const double a = line.einstein_a_s * line.rest_wavelength_m / (4.0 * k::k_pi * b);
          const double chi_d = m.dust_absorption_per_m.empty() ? 0.0 : m.dust_absorption_per_m[c];
          // Sub-steps: Hubble flow changes x along the segment; static media use one piece.
          const double length = seg.length();
          const std::size_t pieces = m.hubble_rate_s == 0.0
              ? 1U
              : std::max<std::size_t>(1U, static_cast<std::size_t>(std::ceil(
                    std::abs(m.hubble_rate_s) * length * m.length_unit_m / (0.1 * b))));
          const double piece = length / static_cast<double>(pieces);
          for (std::size_t q = 0; q < pieces && !interacted; ++q) {
            const double s0 = seg.s_begin + piece * static_cast<double>(q);
            std::array<double, 3> mid{pos[0] + (s0 + 0.5 * piece) * dir[0], pos[1] + (s0 + 0.5 * piece) * dir[1],
                                      pos[2] + (s0 + 0.5 * piece) * dir[2]};
            const auto v = bulk(c, mid);
            const double x = (w - dot(v, dir)) / b;
            const double chi = m.n_hi_m3[c] * sigma_lambda * spectral::voigtH(a, x) / (std::sqrt(k::k_pi) * b) + chi_d;
            const double dtau = chi * piece * m.length_unit_m;
            if (dtau < tau_target || !(chi > 0.0)) {
              tau_target -= dtau;
              covered = s0 + piece - 0.0;
              continue;
            }
            const double l = tau_target / chi / m.length_unit_m;
            for (std::size_t ax = 0; ax < 3U; ++ax) pos[ax] += (s0 + l) * dir[ax];
            interacted = true;
            if (rng() * chi < chi_d) {
              ++dusted;
              alive = false;
              break;
            }
            // Resonant scattering in the atom frame.
            const double u_par = sampleParallelAtomVelocity(x, a, rng);
            const double r = std::sqrt(-std::log(rng()));
            const double ph = 2.0 * k_pi * rng();
            std::array<double, 3> e1{};
            std::array<double, 3> e2{};
            basis(dir, e1, e2);
            std::array<double, 3> u{};
            for (std::size_t ax = 0; ax < 3U; ++ax) {
              u[ax] = u_par * dir[ax] + r * std::cos(ph) * e1[ax] + r * std::sin(ph) * e2[ax];
            }
            std::array<double, 3> out{};
            if (s.phase == LymanAlphaPhaseFunction::kIsotropic) {
              out = isotropicDirection(rng);
            } else {
              double mu = 0.0;
              do {
                mu = 2.0 * rng() - 1.0;
              } while (rng() * 2.0 > 1.0 + mu * mu);
              const double phi = 2.0 * k_pi * rng();
              const double st = std::sqrt(std::max(0.0, 1.0 - mu * mu));
              for (std::size_t ax = 0; ax < 3U; ++ax) {
                out[ax] = mu * dir[ax] + st * (std::cos(phi) * e1[ax] + std::sin(phi) * e2[ax]);
              }
            }
            double x_out = x - dot(u, dir) + dot(u, out);
            if (s.recoil) {
              const double g = k::k_planck_si * line.rest_frequency_hz /
                               (k::k_speed_of_light_si * std::sqrt(2.0 * k::k_boltzmann_si * t * line.atom_mass_kg));
              x_out -= g * (1.0 - dot(dir, out));
            }
            const auto v_here = bulk(c, pos);
            w = x_out * b + dot(v_here, out);
            for (std::size_t ax = 0; ax < 3U; ++ax) {
              result.momentum_rate_per_cell[c][ax] += photon_momentum * (dir[ax] - out[ax]) * rate_per_packet;
            }
            result.scattering_rate_per_cell[c] += rate_per_packet;
            dir = out;
            ++scatterings;
            if (scatterings >= s.max_scatterings) {
              ++truncated;
              alive = false;
            }
          }
          if (interacted) break;
          covered = seg.s_end;
        }
        if (interacted) break;
        if (covered < chunk * (1.0 - 1.0e-12)) {  // left the open domain
          escaped_now = true;
          alive = false;
          break;
        }
        for (std::size_t ax = 0; ax < 3U; ++ax) pos[ax] += covered * dir[ax];
      }
    }
    total_scatterings += scatterings;
    if (escaped_now) {
      ++escaped;
      const double f = (w - s.spectrum_w_min) / (s.spectrum_w_max - s.spectrum_w_min);
      if (f >= 0.0 && f < 1.0) result.escaped_spectrum[static_cast<std::size_t>(f * static_cast<double>(s.spectrum_bins))] += 1.0;
    }
  }
  const double n = static_cast<double>(s.packets);
  result.escape_fraction = static_cast<double>(escaped) / n;
  result.dust_absorbed_fraction = static_cast<double>(dusted) / n;
  result.truncated_fraction = static_cast<double>(truncated) / n;
  result.mean_scatterings = static_cast<double>(total_scatterings) / n;
  return result;
}

}  // namespace cosmosim::observables
