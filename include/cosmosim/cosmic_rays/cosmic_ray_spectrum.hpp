#pragma once

// CHUI-EXT-PHYS-REF-v1 CR04-CR07: momentum-resolved cosmic-ray spectra.
//
// Momentum variable p~ = p / (m c) of the species (proton, electron,
// positron); user-selected fixed logarithmic or explicit bin edges (no
// adaptive groups in pass 01). Distribution f(p~) with
//   n = 4 pi int p~^2 f dp~,   e = 4 pi int p~^2 f T dp~,   T = (sqrt(1+p~^2) - 1) m c^2,
//   P = (4 pi / 3) int p~^2 f (p v) dp~,   p v = m c^2 p~^2 / sqrt(1 + p~^2).
// Representations:
//   kFixedBins          — reference: number per bin, f constant in p~ within a bin.
//   kPiecewisePowerLaw  — advanced parity: number AND energy per bin are the
//     conserved variables; f = f0 (p~/p~_lo)^-q is reconstructed per bin by
//     bisection on the monotone mean kinetic energy <T>(q) = e/n with q in the
//     declared [q_min, q_max]. A ratio outside the admissible range keeps the
//     bound slope and is COUNTED (moments stay authoritative; the slope is a
//     derived reconstruction only). Recalled lead: Girichidis et al. 2020
//     (MNRAS 491, 993) for the n/e power-law scheme; Miniati 2001 (COMPASS).
// Momentum-space transport for dp~/dt = A(p~) (adiabatic -div u / 3 p~,
// losses): edge number flux Phi_N = 4 pi p~_e^2 f_up(p~_e) A(p~_e) from the
// upwind bin, energy flux Phi_E = T(p~_e) Phi_N, plus the in-bin work
// W_i = 4 pi int_bin p~^2 f dT/dp~ A dp~ (24-node Gauss-Legendre in ln p~):
//   dN_i/dt = Phi_N(lo) - Phi_N(hi),  dE_i/dt = Phi_E(lo) - Phi_E(hi) + W_i.
// Fluxes through the outer edges p~_0, p~_N are ledgered (spectral boundary
// losses), inflow from outside the grid is zero. Explicit, positivity refused.
// Pressure from the reconstructed spectrum (no grey P = e/3 assumption).

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"

namespace cosmosim::cosmic_rays {

enum class CrSpecies : std::uint8_t { kProton = 0, kElectron = 1, kPositron = 2 };
enum class CrSpectralRepresentation : std::uint8_t { kFixedBins = 0, kPiecewisePowerLaw = 1 };

[[nodiscard]] double crSpeciesRestEnergyJ(CrSpecies species) noexcept;

struct CrMomentumGrid {
  std::vector<double> edges;  // p~ = p/(m c), strictly increasing, > 0
  [[nodiscard]] std::size_t bins() const noexcept { return edges.empty() ? 0U : edges.size() - 1U; }
  [[nodiscard]] static CrMomentumGrid logarithmic(double p_min, double p_max, std::size_t bins);
  void validate() const;
};

struct CrPowerLawBin {
  double q = 0.0;
  double f0 = 0.0;   // f at the lower edge
  bool clamped = false;
};

struct CrSpectralSettings {
  CrSpecies species = CrSpecies::kProton;
  CrSpectralRepresentation representation = CrSpectralRepresentation::kPiecewisePowerLaw;
  double q_min = -10.0;
  double q_max = 15.0;
};

// Moments of f = f0 (p/p_lo)^-q on [p_lo, p_hi] (number; kinetic energy and
// pressure in units of m c^2 per unit f0).
[[nodiscard]] double powerLawNumber(double p_lo, double p_hi, double q, double f0);
[[nodiscard]] double powerLawEnergy(double p_lo, double p_hi, double q, double f0, double rest_energy);
[[nodiscard]] double powerLawPressure(double p_lo, double p_hi, double q, double f0, double rest_energy);
[[nodiscard]] CrPowerLawBin reconstructPowerLaw(double p_lo, double p_hi, double number, double energy,
                                                double rest_energy, double q_min, double q_max);

// Per-cell spectral state: lanes [bin * cells + cell] (bin-major SoA).
class CrSpectralState {
 public:
  CrSpectralState(CrMomentumGrid grid, CrSpectralSettings settings, std::size_t cells,
                  core::MemoryGovernor* governor);
  [[nodiscard]] static std::uint64_t bytesFor(std::size_t bins, std::size_t cells, CrSpectralRepresentation rep);
  [[nodiscard]] const CrMomentumGrid& grid() const noexcept { return m_grid; }
  [[nodiscard]] const CrSpectralSettings& settings() const noexcept { return m_settings; }
  [[nodiscard]] std::size_t cells() const noexcept { return m_cells; }
  [[nodiscard]] std::span<double> number(std::size_t bin);
  [[nodiscard]] std::span<double> energy(std::size_t bin);  // empty for kFixedBins
  [[nodiscard]] std::span<const double> number(std::size_t bin) const;
  [[nodiscard]] std::span<const double> energy(std::size_t bin) const;
  [[nodiscard]] CrPowerLawBin reconstruct(std::size_t bin, std::size_t cell) const;
  [[nodiscard]] double binEnergy(std::size_t bin, std::size_t cell) const;  // derived for fixed bins
  [[nodiscard]] double totalEnergy(std::size_t cell) const;
  [[nodiscard]] double totalNumber(std::size_t cell) const;
  [[nodiscard]] double pressure(std::size_t cell) const;
  [[nodiscard]] double distribution(std::size_t cell, double p) const;  // f(p~) from the reconstruction

 private:
  CrMomentumGrid m_grid;
  CrSpectralSettings m_settings;
  std::size_t m_cells = 0;
  core::MemoryReservation m_reservation;
  std::vector<double> m_number;
  std::vector<double> m_energy;
};

// dp~/dt as a function of p~ for one cell (adiabatic + losses).
struct CrMomentumRate {
  double adiabatic_coefficient = 0.0;  // A_ad = -div u / 3 (1/T): dp~/dt = A_ad p~
  double radiative_coefficient = 0.0;  // C (1/T): dp~/dt = -C p~ gamma (synchrotron + inverse Compton)
  [[nodiscard]] double operator()(double p) const noexcept {
    return adiabatic_coefficient * p - radiative_coefficient * p * std::sqrt(1.0 + p * p);
  }
};

struct CrSpectralLedger {
  double number_lost_low = 0.0;
  double number_lost_high = 0.0;
  double energy_lost_low = 0.0;
  double energy_lost_high = 0.0;
  double work = 0.0;              // sum of W_i dt (energy change from dp/dt)
  std::uint64_t clamped_slopes = 0U;
};

[[nodiscard]] double momentumTimestep(const CrMomentumGrid& grid, const CrMomentumRate& rate, double cfl);
// Advances one cell's spectrum by dt; returns false (cell untouched) on positivity failure.
[[nodiscard]] bool advanceMomentumSpace(CrSpectralState& state, std::size_t cell, const CrMomentumRate& rate,
                                        double dt, CrSpectralLedger& ledger);

}  // namespace cosmosim::cosmic_rays
