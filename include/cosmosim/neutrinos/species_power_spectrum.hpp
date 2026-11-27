#pragma once

// CHUI-EXT-PHYS-REF-v1 N05: species-resolved density statistics on a periodic
// reference mesh (serial post-processing; uses spectral::fft3d).
//
//   delta_X(x) = rho_X(x) / rho_bar_X - 1 on an n^3 CIC mesh, with rho_bar_X
//   the BACKGROUND mean of the species (for delta-f markers the analytic
//   background density, never the signed marker sum).
//   P_XY(k) = V <Re(delta_X,k delta_Y,k^*)> / N^2 averaged in |k| shells,
//   each mode divided by the CIC window W(k)^2 = prod_i sinc^4(k_i dx / 2)
//   when deconvolution is requested.
//   Total matter: P_m = f_cb^2 P_cb,cb + 2 f_cb f_nu P_cb,nu + f_nu^2 P_nu,nu,
//   f_X = rho_bar_X / (rho_bar_cb + rho_bar_nu).
//   Noise: ordinary particles P_shot = V sum m_i^2 / (sum m_i)^2; delta-f
//   markers P_noise = V sum (w_i m_i)^2 / M_bg^2 with M_bg = rho_bar_nu V.
// Linear-response / transfer-ratio neutrino fields are not sampled: their
// spectra follow exactly from P_nu,nu = R^2 P_cb,cb and P_cb,nu = R P_cb,cb
// (crossSpectraFromRatio), so the estimator family is selected explicitly
// and never substituted silently.

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace cosmosim::neutrinos {

struct SpeciesMesh {
  std::size_t n = 0;
  double box_size = 0.0;                       // same length unit as positions
  std::vector<double> delta;                   // n^3 row-major density contrast
};

// CIC deposit of (signed) masses; mean_density_per_volume is the background
// mean used to form the contrast. Positions are wrapped periodically.
[[nodiscard]] SpeciesMesh depositCloudInCell(std::span<const std::array<double, 3>> positions,
                                             std::span<const double> masses, std::size_t n, double box_size,
                                             double mean_density_per_volume);

struct SpectrumBin {
  double k_center = 0.0;   // 2 pi / length unit
  double power = 0.0;      // length^3
  std::uint64_t modes = 0;
};

[[nodiscard]] std::vector<SpectrumBin> crossPowerSpectrum(const SpeciesMesh& a, const SpeciesMesh& b,
                                                          std::size_t bins, bool deconvolve_cic);

struct SpeciesSpectra {
  std::vector<SpectrumBin> cold_cold;
  std::vector<SpectrumBin> cold_nu;
  std::vector<SpectrumBin> nu_nu;
  std::vector<SpectrumBin> total;
};

[[nodiscard]] SpeciesSpectra combineSpeciesSpectra(std::vector<SpectrumBin> cold_cold, std::vector<SpectrumBin> cold_nu,
                                                   std::vector<SpectrumBin> nu_nu, double f_nu);
// Linear-response/transfer-ratio family: P_cn = R P_cc, P_nn = R^2 P_cc with
// R evaluated at each bin center through `ratio_at_k`.
[[nodiscard]] SpeciesSpectra crossSpectraFromRatio(const std::vector<SpectrumBin>& cold_cold,
                                                   std::span<const double> ratio_per_bin, double f_nu);

[[nodiscard]] double shotNoisePower(std::span<const double> masses, double box_volume);
[[nodiscard]] double deltaFNoisePower(std::span<const double> signed_masses, double background_mass, double box_volume);

}  // namespace cosmosim::neutrinos
