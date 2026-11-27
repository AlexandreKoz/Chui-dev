// CHUI-EXT-PHYS-REF-v1 N05 species statistics (source only, NOT compiled in pass 01).

#include "cosmosim/neutrinos/species_power_spectrum.hpp"

#include <cmath>
#include <stdexcept>

#include "cosmosim/spectral/reference_fft.hpp"

namespace cosmosim::neutrinos {
namespace {

constexpr double k_two_pi = 6.283185307179586476925286766559005768;

[[nodiscard]] double sincHalf(double x) {
  // sin(x/2) / (x/2)
  const double h = 0.5 * x;
  return std::abs(h) < 1.0e-8 ? 1.0 : std::sin(h) / h;
}

}  // namespace

SpeciesMesh depositCloudInCell(std::span<const std::array<double, 3>> positions, std::span<const double> masses,
                               std::size_t n, double box_size, double mean_density_per_volume) {
  if (positions.size() != masses.size()) throw std::invalid_argument("CIC deposit: positions/masses size mismatch");
  if (!spectral::isPowerOfTwo(n) || !(box_size > 0.0) || !(mean_density_per_volume > 0.0)) {
    throw std::invalid_argument("CIC deposit: n must be a power of two, box and mean density > 0");
  }
  SpeciesMesh mesh;
  mesh.n = n;
  mesh.box_size = box_size;
  mesh.delta.assign(n * n * n, 0.0);
  const double dx = box_size / static_cast<double>(n);
  const double cell_volume = dx * dx * dx;
  for (std::size_t p = 0; p < positions.size(); ++p) {
    std::array<std::size_t, 3> i0{};
    std::array<double, 3> frac{};
    for (std::size_t d = 0; d < 3U; ++d) {
      double u = positions[p][d] / dx - 0.5;  // cell-centered CIC
      u -= std::floor(u / static_cast<double>(n)) * static_cast<double>(n);
      const double cell = std::floor(u);
      frac[d] = u - cell;
      i0[d] = static_cast<std::size_t>(cell) % n;
    }
    for (int cx = 0; cx < 2; ++cx) {
      for (int cy = 0; cy < 2; ++cy) {
        for (int cz = 0; cz < 2; ++cz) {
          const double w = (cx ? frac[0] : 1.0 - frac[0]) * (cy ? frac[1] : 1.0 - frac[1]) *
                           (cz ? frac[2] : 1.0 - frac[2]);
          const std::size_t ix = (i0[0] + static_cast<std::size_t>(cx)) % n;
          const std::size_t iy = (i0[1] + static_cast<std::size_t>(cy)) % n;
          const std::size_t iz = (i0[2] + static_cast<std::size_t>(cz)) % n;
          mesh.delta[(ix * n + iy) * n + iz] += w * masses[p];
        }
      }
    }
  }
  for (double& value : mesh.delta) value = value / (cell_volume * mean_density_per_volume) - 1.0;
  return mesh;
}

std::vector<SpectrumBin> crossPowerSpectrum(const SpeciesMesh& a, const SpeciesMesh& b, std::size_t bins,
                                            bool deconvolve_cic) {
  if (a.n != b.n || a.box_size != b.box_size || a.delta.size() != b.delta.size() || bins == 0U) {
    throw std::invalid_argument("cross spectrum needs matching meshes and >= 1 bin");
  }
  const std::size_t n = a.n;
  std::vector<std::complex<double>> fa(a.delta.begin(), a.delta.end());
  std::vector<std::complex<double>> fb(b.delta.begin(), b.delta.end());
  spectral::fft3d(fa, n, false);
  spectral::fft3d(fb, n, false);
  const double k_fund = k_two_pi / a.box_size;
  const double k_nyq = k_fund * static_cast<double>(n / 2U);
  const double volume = a.box_size * a.box_size * a.box_size;
  const double norm = volume / (static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n) *
                                static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n));
  const double dx = a.box_size / static_cast<double>(n);
  std::vector<SpectrumBin> out(bins);
  std::vector<double> k_sum(bins, 0.0);
  for (std::size_t ix = 0; ix < n; ++ix) {
    for (std::size_t iy = 0; iy < n; ++iy) {
      for (std::size_t iz = 0; iz < n; ++iz) {
        const double kx = k_fund * static_cast<double>(spectral::fftFrequencyIndex(ix, n));
        const double ky = k_fund * static_cast<double>(spectral::fftFrequencyIndex(iy, n));
        const double kz = k_fund * static_cast<double>(spectral::fftFrequencyIndex(iz, n));
        const double kk = std::sqrt(kx * kx + ky * ky + kz * kz);
        if (kk == 0.0 || kk > k_nyq) continue;
        const auto bin = std::min(static_cast<std::size_t>((kk - 0.5 * k_fund) / (k_nyq - 0.5 * k_fund) *
                                                           static_cast<double>(bins)),
                                  bins - 1U);
        const std::size_t index = (ix * n + iy) * n + iz;
        double value = norm * (fa[index] * std::conj(fb[index])).real();
        if (deconvolve_cic) {
          const double w = sincHalf(kx * dx) * sincHalf(ky * dx) * sincHalf(kz * dx);
          value /= std::pow(w, 4);  // CIC window W = sinc^2 per axis, power ~ W^2
        }
        out[bin].power += value;
        out[bin].modes += 1U;
        k_sum[bin] += kk;
      }
    }
  }
  for (std::size_t i = 0; i < bins; ++i) {
    if (out[i].modes > 0U) {
      out[i].power /= static_cast<double>(out[i].modes);
      out[i].k_center = k_sum[i] / static_cast<double>(out[i].modes);
    }
  }
  return out;
}

SpeciesSpectra combineSpeciesSpectra(std::vector<SpectrumBin> cold_cold, std::vector<SpectrumBin> cold_nu,
                                     std::vector<SpectrumBin> nu_nu, double f_nu) {
  if (cold_cold.size() != cold_nu.size() || cold_cold.size() != nu_nu.size()) {
    throw std::invalid_argument("species spectra must share bins");
  }
  if (!(f_nu >= 0.0) || !(f_nu < 1.0)) throw std::invalid_argument("f_nu must lie in [0, 1)");
  SpeciesSpectra out;
  out.total = cold_cold;
  const double f_cb = 1.0 - f_nu;
  for (std::size_t i = 0; i < cold_cold.size(); ++i) {
    out.total[i].power = f_cb * f_cb * cold_cold[i].power + 2.0 * f_cb * f_nu * cold_nu[i].power +
                         f_nu * f_nu * nu_nu[i].power;
  }
  out.cold_cold = std::move(cold_cold);
  out.cold_nu = std::move(cold_nu);
  out.nu_nu = std::move(nu_nu);
  return out;
}

SpeciesSpectra crossSpectraFromRatio(const std::vector<SpectrumBin>& cold_cold, std::span<const double> ratio_per_bin,
                                     double f_nu) {
  if (ratio_per_bin.size() != cold_cold.size()) throw std::invalid_argument("ratio bins must match spectrum bins");
  std::vector<SpectrumBin> cn = cold_cold;
  std::vector<SpectrumBin> nn = cold_cold;
  for (std::size_t i = 0; i < cold_cold.size(); ++i) {
    cn[i].power = ratio_per_bin[i] * cold_cold[i].power;
    nn[i].power = ratio_per_bin[i] * ratio_per_bin[i] * cold_cold[i].power;
  }
  return combineSpeciesSpectra(cold_cold, std::move(cn), std::move(nn), f_nu);
}

double shotNoisePower(std::span<const double> masses, double box_volume) {
  double sum = 0.0;
  double sum2 = 0.0;
  for (const double m : masses) {
    sum += m;
    sum2 += m * m;
  }
  if (!(sum > 0.0)) throw std::invalid_argument("shot noise needs positive total mass");
  return box_volume * sum2 / (sum * sum);
}

double deltaFNoisePower(std::span<const double> signed_masses, double background_mass, double box_volume) {
  if (!(background_mass > 0.0)) throw std::invalid_argument("delta-f noise needs the background mass");
  double sum2 = 0.0;
  for (const double m : signed_masses) sum2 += m * m;
  return box_volume * sum2 / (background_mass * background_mass);
}

}  // namespace cosmosim::neutrinos
