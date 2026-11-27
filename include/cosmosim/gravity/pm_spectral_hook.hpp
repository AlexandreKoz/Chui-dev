#pragma once

// CHUI-EXT-PHYS-REF-v1 N03: optional isotropic spectral source hook of the
// periodic PM Poisson operator. It lets a species that is NOT represented by
// deposited particles (massive neutrinos in the transfer-ratio and
// linear-response models) contribute to the long-range potential without a
// duplicate PM engine:
//
//   Psi_k = -4 pi G W(k)^-2 S_split(k) M(|k|) delta_rho_cb,k / k^2,   k != 0,
//   M(|k|) = 1 + (rho_bar_nu / rho_bar_cb) R(|k|),  R = delta_nu / delta_cb.
//
// Contract:
//   * The deposited density is the cold (CDM + baryon) density only; the hook
//     multiplies the already mean-subtracted mode, so k = 0 stays exactly 0
//     and the neutrino mean density enters only H(a) (S01 authority).
//   * The multiplier commutes with the assignment-window deconvolution and
//     the TreePM Gaussian long-range split; short-range tree forces carry no
//     neutrino contribution (the model's declared scale limit: the error is
//     bounded by M(k_split) - 1, reported by the hook).
//   * beginSolve/finalizeObservation are COLLECTIVE: every rank of the PM
//     communicator calls them exactly once per periodic Poisson solve, in the
//     same order, including ranks that own no spectral modes.
//   * Observation is deterministic: the PM solver visits every locally owned
//     logical mode once, SERIALLY in storage order, adds
//     multiplicity * |delta_rho_k|^2 and multiplicity to bin `binIndex(|k|)`
//     (multiplicity 1 or 2 for the r2c half spectrum), then passes the
//     rank-local bin sums to finalizeObservation, which reduces them across
//     ranks. k is in 1/(PM length unit) (comoving Mpc for the TreePM runtime).
//   * multiplier() is const and thread-safe (called inside OpenMP loops).
//   * Isolated/open and CUDA PM solves refuse a non-null hook.

#include <cstddef>
#include <cstdint>
#include <span>

namespace cosmosim::gravity {

class PmSpectralSourceHook {
 public:
  static constexpr std::size_t k_no_bin = static_cast<std::size_t>(-1);
  virtual ~PmSpectralSourceHook() = default;
  virtual void beginSolve(double scale_factor) = 0;
  [[nodiscard]] virtual bool needsObservation() const noexcept = 0;
  [[nodiscard]] virtual std::size_t binCount() const noexcept = 0;
  // Bin of |k| or k_no_bin (modes outside the observed range are not binned).
  [[nodiscard]] virtual std::size_t binIndex(double k_magnitude) const noexcept = 0;
  // Collective. Spans are the rank-local sums (size binCount(), empty when
  // needsObservation() is false).
  virtual void finalizeObservation(std::span<const double> local_bin_power,
                                   std::span<const double> local_bin_modes) = 0;
  // Throws unless [k_min, k_max] (smallest nonzero and largest mesh |k|) lies
  // inside the tabulated multiplier domain. Called before any mode loop, so
  // the noexcept multiplier below never meets an untabulated wavenumber.
  virtual void requireWavenumberDomain(double k_min, double k_max) const = 0;
  // Isotropic multiplier for |k| > 0 (domain pre-verified).
  [[nodiscard]] virtual double multiplier(double k_magnitude) const noexcept = 0;
  // Rank-uniform identity of the hook state (collective fingerprint lane).
  [[nodiscard]] virtual std::uint64_t fingerprint() const noexcept = 0;
};

}  // namespace cosmosim::gravity
