// CHUI-EXT-PHYS-REF-v1 C01-C03 CMB primary inputs and secondaries (source
// only, NOT compiled or executed in pass 01).

#include "cosmosim/observables/cmb_secondaries.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <map>
#include <sstream>
#include <stdexcept>

#include "core/internal/sha256.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/radiation/radiation_angular.hpp"
#include "cosmosim/spectral/reference_fft.hpp"

namespace cosmosim::observables {
namespace {

namespace k = core::constants;

[[nodiscard]] std::string trim(const std::string& s) {
  std::size_t lo = 0;
  std::size_t hi = s.size();
  while (lo < hi && std::isspace(static_cast<unsigned char>(s[lo]))) ++lo;
  while (hi > lo && std::isspace(static_cast<unsigned char>(s[hi - 1]))) --hi;
  return s.substr(lo, hi - lo);
}

[[nodiscard]] double uniform(std::uint64_t seed, std::uint64_t stream, std::uint64_t counter) {
  const std::uint64_t bits = radiation::radiationCounterRandom(seed, stream, counter);
  return (static_cast<double>(bits >> 11U) + 0.5) * (1.0 / 9007199254740992.0);
}

[[nodiscard]] double gaussian(std::uint64_t seed, std::uint64_t stream, std::uint64_t counter) {
  const double r = std::sqrt(-2.0 * std::log(uniform(seed, stream, 2U * counter)));
  return r * std::cos(2.0 * k::k_pi * uniform(seed, stream, 2U * counter + 1U));
}

// Bilinear periodic interpolation on an n x n plane of side L centred on 0.
[[nodiscard]] double samplePeriodic(const std::vector<double>& field, std::size_t n, double side, double x, double y) {
  const double fx = (x / side + 0.5) * static_cast<double>(n) - 0.5;
  const double fy = (y / side + 0.5) * static_cast<double>(n) - 0.5;
  const double ix = std::floor(fx);
  const double iy = std::floor(fy);
  const double tx = fx - ix;
  const double ty = fy - iy;
  const auto wrap = [n](double i) {
    const long m = static_cast<long>(n);
    return static_cast<std::size_t>(((static_cast<long>(i) % m) + m) % m);
  };
  const std::size_t x0 = wrap(ix);
  const std::size_t x1 = wrap(ix + 1.0);
  const std::size_t y0 = wrap(iy);
  const std::size_t y1 = wrap(iy + 1.0);
  return (1.0 - tx) * (1.0 - ty) * field[y0 * n + x0] + tx * (1.0 - ty) * field[y0 * n + x1] +
         (1.0 - tx) * ty * field[y1 * n + x0] + tx * ty * field[y1 * n + x1];
}

}  // namespace

double CmbSpectraTable::clTT(double l) const {
  if (!(l >= ell.front()) || !(l <= ell.back())) throw std::out_of_range("multipole outside the C_l table");
  const auto it = std::upper_bound(ell.begin(), ell.end(), l);
  std::size_t hi = static_cast<std::size_t>(std::distance(ell.begin(), it));
  hi = std::clamp<std::size_t>(hi, 1U, ell.size() - 1U);
  const std::size_t lo = hi - 1U;
  const double t = (l - ell[lo]) / (ell[hi] - ell[lo]);
  return (1.0 - t) * cl_tt[lo] + t * cl_tt[hi];
}

CmbSpectraTable parseCmbSpectraTable(const std::string& bytes, const std::string& origin) {
  CmbSpectraTable table;
  std::map<std::string, std::string> meta;
  std::istringstream stream(bytes);
  std::string line;
  std::vector<std::vector<double>> rows;
  while (std::getline(stream, line)) {
    const std::string t = trim(line);
    if (t.empty()) continue;
    if (t.front() == '#') {
      const auto colon = t.find(':');
      if (colon != std::string::npos) meta[trim(t.substr(1, colon - 1))] = trim(t.substr(colon + 1));
      continue;
    }
    std::istringstream row(t);
    std::vector<double> values;
    double v = 0.0;
    while (row >> v) values.push_back(v);
    if (values.size() < 4U) throw std::runtime_error("C_l table '" + origin + "': rows need ell TT EE TE");
    rows.push_back(values);
  }
  const auto need = [&](const char* key) {
    const auto it = meta.find(key);
    if (it == meta.end()) throw std::runtime_error("C_l table '" + origin + "': metadata '" + key + "' missing");
    return it->second;
  };
  if (need("schema") != "chui_cmb_cl_v1") throw std::runtime_error("C_l table '" + origin + "': wrong schema");
  if (need("units") != "muK2") throw std::runtime_error("C_l table '" + origin + "': units must be muK2");
  const std::string norm = need("normalization");
  if (norm != "dl" && norm != "cl") throw std::runtime_error("C_l table '" + origin + "': normalization dl|cl");
  if (need("columns").rfind("ell TT EE TE", 0) != 0) {
    throw std::runtime_error("C_l table '" + origin + "': columns must start with 'ell TT EE TE'");
  }
  table.source_code = need("source_code");
  table.source_version = need("source_version");
  table.lensed = need("lensed") == "yes";
  for (const auto& r : rows) {
    const double l = r[0];
    if (!(l >= 1.0) || (!table.ell.empty() && !(l > table.ell.back()))) {
      throw std::runtime_error("C_l table '" + origin + "': ell must start at >= 1 and increase");
    }
    const double factor = norm == "dl" ? 2.0 * k::k_pi / (l * (l + 1.0)) : 1.0;
    table.ell.push_back(l);
    table.cl_tt.push_back(r[1] * factor);
    table.cl_ee.push_back(r[2] * factor);
    table.cl_te.push_back(r[3] * factor);
  }
  if (table.ell.size() < 2U) throw std::runtime_error("C_l table '" + origin + "': needs >= 2 rows");
  table.sha256 = core::internal::sha256Hex(bytes);
  return table;
}

FlatSkyMap gaussianRealizationFlatSky(const CmbSpectraTable& table, std::size_t n, double side, std::uint64_t seed) {
  if (!spectral::isPowerOfTwo(n) || !(side > 0.0)) throw std::invalid_argument("flat-sky map needs n = 2^m, side > 0");
  const double area = side * side;
  std::vector<std::complex<double>> a(n * n, {0.0, 0.0});
  for (std::size_t iy = 0; iy < n; ++iy) {
    for (std::size_t ix = 0; ix < n; ++ix) {
      const long mx = spectral::fftFrequencyIndex(ix, n);
      const long my = spectral::fftFrequencyIndex(iy, n);
      if (mx == 0 && my == 0) continue;  // monopole not part of the anisotropy field
      // Canonical half: generate only for (my > 0) or (my == 0 and mx > 0), plus self-conjugate modes.
      const std::size_t cx = (n - ix) % n;
      const std::size_t cy = (n - iy) % n;
      const bool self = cx == ix && cy == iy;
      const bool canonical = self || iy < cy || (iy == cy && ix < cx);
      if (!canonical) continue;
      const double l = 2.0 * k::k_pi / side * std::sqrt(static_cast<double>(mx * mx + my * my));
      const double cl = table.clTT(l);
      const std::uint64_t stream = static_cast<std::uint64_t>(iy * n + ix);
      if (self) {
        a[iy * n + ix] = {std::sqrt(area * cl) * gaussian(seed, stream, 0U), 0.0};
      } else {
        const std::complex<double> z(gaussian(seed, stream, 0U), gaussian(seed, stream, 1U));
        a[iy * n + ix] = std::sqrt(0.5 * area * cl) * z;
        a[cy * n + cx] = std::conj(a[iy * n + ix]);
      }
    }
  }
  spectral::fft2d(a, n, true);  // includes 1/n^2
  FlatSkyMap map;
  map.n = n;
  map.side_rad = side;
  map.values.resize(n * n);
  const double scale = static_cast<double>(n * n) / area;
  for (std::size_t i = 0; i < n * n; ++i) map.values[i] = scale * a[i].real();
  return map;
}

SzIntegrals integrateSz(const std::vector<SzSegment>& segments) {
  SzIntegrals out;
  const double me_c2 = k::k_electron_mass_si * k::k_speed_of_light_si * k::k_speed_of_light_si;
  for (const SzSegment& s : segments) {
    if (!(s.path_proper_m >= 0.0) || !(s.electron_density_m3 >= 0.0) || !(s.electron_temperature_k >= 0.0)) {
      throw std::invalid_argument("SZ segment needs non-negative path, density and temperature");
    }
    const double dtau = k::k_thomson_cross_section_si * s.electron_density_m3 * s.path_proper_m;
    out.optical_depth += dtau;
    out.compton_y += dtau * k::k_boltzmann_si * s.electron_temperature_k / me_c2;
    out.ksz_delta_t_over_t -= dtau * s.radial_velocity_m_s / k::k_speed_of_light_si;
  }
  return out;
}

double tszSpectralFactor(double frequency_hz, double t_cmb_k) {
  if (!(frequency_hz > 0.0) || !(t_cmb_k > 0.0)) throw std::invalid_argument("tSZ factor needs nu, T > 0");
  const double x = k::k_planck_si * frequency_hz / (k::k_boltzmann_si * t_cmb_k);
  return x / std::tanh(0.5 * x) - 4.0;
}

double integrateIsw(const std::vector<IswSegment>& segments, double a0, double a1,
                    const core::LambdaCdmBackground& background) {
  if (!(a1 > a0)) throw std::invalid_argument("ISW needs two potential snapshots a0 < a1");
  const double c = k::k_speed_of_light_si;
  double sum = 0.0;
  for (const IswSegment& s : segments) {
    if (s.scale_factor < a0 || s.scale_factor > a1) {
      throw std::out_of_range("ISW segment epoch outside the potential snapshots (no temporal extrapolation)");
    }
    const double dphi_dt = (s.potential_a1 - s.potential_a0) / (a1 - a0) * s.scale_factor *
                           background.hubbleSi(s.scale_factor);
    sum += 2.0 / (c * c) * dphi_dt * s.scale_factor * s.comoving_length_m / c;
  }
  return sum;
}

LensingResult rayShootLensPlanes(const std::vector<LensPlane>& planes_in, double chi_s, double h0, double omega_m,
                                 std::size_t n_rays, double field) {
  if (planes_in.empty() || n_rays < 3U || !(field > 0.0) || !(chi_s > 0.0)) {
    throw std::invalid_argument("lens ray shooting needs planes, >= 3 rays per side, field and source distance");
  }
  std::vector<LensPlane> planes;
  for (const LensPlane& p : planes_in) {
    if (p.chi_m < chi_s) planes.push_back(p);
  }
  std::sort(planes.begin(), planes.end(), [](const LensPlane& a, const LensPlane& b) { return a.chi_m < b.chi_m; });
  const double c = k::k_speed_of_light_si;
  const double factor = 3.0 * h0 * h0 * omega_m / (c * c);
  // Deflection fields per plane from the 2-D periodic Poisson equation.
  std::vector<std::vector<double>> alpha_x(planes.size());
  std::vector<std::vector<double>> alpha_y(planes.size());
  for (std::size_t j = 0; j < planes.size(); ++j) {
    const LensPlane& p = planes[j];
    if (!spectral::isPowerOfTwo(p.n) || p.sigma_m.size() != p.n * p.n || !(p.side_m > 0.0)) {
      throw std::invalid_argument("lens plane needs n = 2^m and a matching sigma field");
    }
    std::vector<std::complex<double>> s(p.sigma_m.begin(), p.sigma_m.end());
    spectral::fft2d(s, p.n, false);
    std::vector<std::complex<double>> ax(p.n * p.n);
    std::vector<std::complex<double>> ay(p.n * p.n);
    for (std::size_t iy = 0; iy < p.n; ++iy) {
      for (std::size_t ix = 0; ix < p.n; ++ix) {
        const double kx = 2.0 * k::k_pi / p.side_m * static_cast<double>(spectral::fftFrequencyIndex(ix, p.n));
        const double ky = 2.0 * k::k_pi / p.side_m * static_cast<double>(spectral::fftFrequencyIndex(iy, p.n));
        const double k2 = kx * kx + ky * ky;
        const std::size_t i = iy * p.n + ix;
        if (k2 == 0.0) {
          ax[i] = ay[i] = {0.0, 0.0};
          continue;
        }
        const std::complex<double> psi = -factor * s[i] / k2;
        ax[i] = std::complex<double>(0.0, kx) * psi;
        ay[i] = std::complex<double>(0.0, ky) * psi;
      }
    }
    spectral::fft2d(ax, p.n, true);
    spectral::fft2d(ay, p.n, true);
    alpha_x[j].resize(p.n * p.n);
    alpha_y[j].resize(p.n * p.n);
    for (std::size_t i = 0; i < p.n * p.n; ++i) {
      alpha_x[j][i] = ax[i].real();
      alpha_y[j][i] = ay[i].real();
    }
  }
  LensingResult out;
  out.n = n_rays;
  out.field_rad = field;
  const std::size_t nr = n_rays * n_rays;
  out.kappa_born.assign(nr, 0.0);
  out.kappa_multiplane.assign(nr, 0.0);
  out.deflection_x_rad.assign(nr, 0.0);
  out.deflection_y_rad.assign(nr, 0.0);
  const double dtheta = field / static_cast<double>(n_rays);
  std::vector<double> beta_x(nr);
  std::vector<double> beta_y(nr);
  for (std::size_t ry = 0; ry < n_rays; ++ry) {
    for (std::size_t rx = 0; rx < n_rays; ++rx) {
      const double t0x = -0.5 * field + (static_cast<double>(rx) + 0.5) * dtheta;
      const double t0y = -0.5 * field + (static_cast<double>(ry) + 0.5) * dtheta;
      double tx = t0x;
      double ty = t0y;
      double x = 0.0;
      double y = 0.0;
      double chi_prev = 0.0;
      double born = 0.0;
      for (std::size_t j = 0; j < planes.size(); ++j) {
        const LensPlane& p = planes[j];
        x += (p.chi_m - chi_prev) * tx;
        y += (p.chi_m - chi_prev) * ty;
        chi_prev = p.chi_m;
        born += 0.5 * factor * p.chi_m * (chi_s - p.chi_m) / chi_s *
                samplePeriodic(p.sigma_m, p.n, p.side_m, p.chi_m * t0x, p.chi_m * t0y);
        tx -= samplePeriodic(alpha_x[j], p.n, p.side_m, x, y);
        ty -= samplePeriodic(alpha_y[j], p.n, p.side_m, x, y);
      }
      x += (chi_s - chi_prev) * tx;
      y += (chi_s - chi_prev) * ty;
      const std::size_t i = ry * n_rays + rx;
      beta_x[i] = x / chi_s;
      beta_y[i] = y / chi_s;
      out.kappa_born[i] = born;
      out.deflection_x_rad[i] = t0x - beta_x[i];
      out.deflection_y_rad[i] = t0y - beta_y[i];
    }
  }
  for (std::size_t ry = 0; ry < n_rays; ++ry) {
    for (std::size_t rx = 0; rx < n_rays; ++rx) {
      const std::size_t xm = rx == 0 ? 0 : rx - 1;
      const std::size_t xp = rx + 1 == n_rays ? rx : rx + 1;
      const std::size_t ym = ry == 0 ? 0 : ry - 1;
      const std::size_t yp = ry + 1 == n_rays ? ry : ry + 1;
      const double dbx = (beta_x[ry * n_rays + xp] - beta_x[ry * n_rays + xm]) / (static_cast<double>(xp - xm) * dtheta);
      const double dby = (beta_y[yp * n_rays + rx] - beta_y[ym * n_rays + rx]) / (static_cast<double>(yp - ym) * dtheta);
      out.kappa_multiplane[ry * n_rays + rx] = 1.0 - 0.5 * (dbx + dby);
    }
  }
  return out;
}

}  // namespace cosmosim::observables
