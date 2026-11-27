#pragma once

// CHUI-EXT-PHYS-REF-v1 L02/L03: frequency-dependent resonant Ly-a packet
// transport, UNACCELERATED reference (no core skipping, no diffusion hybrid).
//
// Frequency: each packet carries its LAB-frame velocity offset w = c (nu -
// nu0)/nu0 [m/s]; in a cell with bulk velocity V and Doppler parameter
// b = sqrt(2kT/m_H) the local frequency is x = (w - V . k) / b. Opacity
// chi = n_HI sigma_int lambda0 H(a, x)/(sqrt(pi) b) + chi_dust (absorption only).
// Scattering (derived): the parallel atom velocity u_par (units of b) is
// sampled EXACTLY from f(u) ~ exp(-u^2) / ((x - u)^2 + a^2) by a two-component
// rejection envelope split at u_s = x - d, d = max(1, |x|/2) (x >= 0; mirrored
// for x < 0): a Gaussian component on u < u_s accepted with
// d^2 / ((x-u)^2 + a^2), and a truncated Cauchy component on u >= u_s
// accepted with exp(-u^2) / max(exp(-u_s^2) [u_s > 0], 1). Perpendicular
// components are Gaussian (variance 1/2). New direction: isotropic or dipole
// (P ~ 1 + mu^2), named options (no unverified core/wing mixtures). Atom-frame
// coherent scattering: x_out = x_in - u . k_in + u . k_out, optional recoil
// x_out -= g (1 - k_in . k_out), g = h nu0 / (c sqrt(2 k T m_H)) (derived).
// Cosmological/Hubble redshifting inside the medium: optional uniform Hubble
// rate H_local adds the bulk velocity field V = H r about the source.
// Estimators (L03), each separately reported: escaped spectrum (lab w),
// escape fraction, scatterings per cell (scattering-rate estimator
// = scatterings x photon rate per packet), momentum deposition per cell
// (h nu0 / c (k_in - k_out) per scattering), dust-absorbed energy, and a
// truncation ledger (packets reaching max_scatterings; never dropped). A
// valid escape spectrum does NOT validate the force or Wouthuysen-Field
// estimators; dynamical coupling and spin coupling stay FAIL_CLOSED.

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/lightcone/ray_geometry.hpp"

namespace cosmosim::observables {

enum class LymanAlphaPhaseFunction : unsigned char { kIsotropic = 0, kDipole = 1 };

struct ResonantMedium {
  std::span<const double> n_hi_m3;
  std::span<const double> temperature_k;
  std::span<const double> vx, vy, vz;   // bulk velocity m/s
  std::span<const double> dust_absorption_per_m;  // may be empty (no dust)
  double length_unit_m = 1.0;          // geometry length -> proper metres
  double hubble_rate_s = 0.0;          // optional uniform Hubble flow about the source
};

struct ResonantSettings {
  std::uint64_t seed = 7U;
  std::uint64_t packets = 10000U;
  std::uint64_t max_scatterings = 100000000ULL;
  LymanAlphaPhaseFunction phase = LymanAlphaPhaseFunction::kIsotropic;
  bool recoil = false;
  double initial_w_m_s = 0.0;     // emission frequency offset (line centre by default)
  double photon_rate_s = 1.0;     // source photon emission rate (for rate estimators)
  double spectrum_w_min = -1.0e6;
  double spectrum_w_max = 1.0e6;
  std::size_t spectrum_bins = 200;
};

struct ResonantResult {
  std::vector<double> escaped_spectrum;   // packets per bin (lab w)
  double escape_fraction = 0.0;
  double dust_absorbed_fraction = 0.0;
  double truncated_fraction = 0.0;
  std::vector<double> scattering_rate_per_cell;      // scatterings / s
  std::vector<std::array<double, 3>> momentum_rate_per_cell;  // N (kg m s^-2)
  double mean_scatterings = 0.0;
};

// Exact sampler of the parallel atom velocity (units of b), uniform draws from `next`.
template <typename Uniform>
[[nodiscard]] double sampleParallelAtomVelocity(double x, double a, Uniform&& next);

[[nodiscard]] ResonantResult runResonantLymanAlpha(const lightcone::LeafGeometry& geometry, const ResonantMedium& medium,
                                                   const std::array<double, 3>& source, const ResonantSettings& settings);

}  // namespace cosmosim::observables

#include "cosmosim/observables/resonant_lyman_alpha_sampler.hpp"
