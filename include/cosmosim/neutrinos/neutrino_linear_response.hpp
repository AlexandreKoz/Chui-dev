#pragma once

// CHUI-EXT-PHYS-REF-v1 N03: massive-neutrino long-range source models for
// the periodic PM operator (gravity::PmSpectralSourceHook). Two DISTINCT,
// separately admitted models:
//
// (a) kTransferRatio — instantaneous ratio R(k, a) = delta_ncdm / delta_cb
//     interpolated from versioned Boltzmann tables (spectral::
//     BoltzmannTransferSet). The neutrino perturbation is slaved to the
//     current cold mode with a PRESCRIBED linear-theory correlation; it does
//     not respond to the simulated potential. Never labelled linear response.
//
// (b) kLinearResponse — Newtonian collisionless linear response to the
//     evolving cold + neutrino potential (derivation in
//     docs/architecture/chui_extended_physics_numerical_contracts_pass01.md,
//     section N03). In superconformal time ds = dt / a^2, comoving x, with
//     u = q / m the conserved canonical velocity of the nonrelativistic
//     Fermi-Dirac distribution and F(y) = <sinc(y x)>_FD its characteristic
//     function (x = q c / (k T_nu,0), y = k u_0 tau, u_0 = k T_nu,0 / (m c)):
//
//       delta_nu(k, s) = [delta_nu,i - theta_nu,i (s - s_i)] F(k u_0 (s - s_i))
//         + int_{s_i}^{s} ds' (s - s') F(k u_0 (s - s')) 4 pi G a(s')^4
//               [rho_cb(s') delta_cb(k, s') + rho_nu(s') delta_nu(k, s')],
//
//     theta_s = i k . dx/ds (continuity d delta/ds = -theta_s; theta_s = a
//     theta_conformal). The cold history enters through the PHASE-COHERENT
//     amplitude approximation delta_cb(k, s') = delta_cb(k, s) A(k, s')/A(k, s)
//     with A the measured rms amplitude of the deposited cold density in |k|
//     bins (Ali-Haimoud & Bird 2013 style; recalled lead). The (s - s')
//     weight vanishes at s' = s, so the update is explicit. Domain: all
//     eigenstates degenerate (one FD velocity scale), <|p|>/(m c) below the
//     admitted bound at a_begin, sub-horizon PM wavenumbers, periodic box.
//
// State: per |k| bin a bounded history (s, a, A, R). It is authoritative run
// state: persisted as the whole-block module sidecar `chui_neutrino_response_v1`
// (restart round trip, replicated identically on every rank after the
// collective reduction) and replayed deterministically on retry: a solve at a
// scale factor <= the newest history entry truncates the later entries
// before recomputing (rollback-safe without a separate journal). A full
// history REFUSES (never drops entries).

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/gravity/pm_spectral_hook.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/spectral/boltzmann_tables.hpp"

namespace cosmosim::neutrinos {

inline constexpr const char* k_neutrino_response_module = "chui_neutrino_response_v1";
inline constexpr std::uint32_t k_neutrino_response_schema = 1U;

// Characteristic function of the nonrelativistic FD velocity distribution,
// F(y) = int x^2 f(x) sinc(x y) dx / int x^2 f(x) dx, f = 1/(e^x + 1).
// Table on y in [0, 50] (step 0.01, 4-point Lagrange) and the asymptotic
// series F ~ (1/(2y^4) + 1/(2y^6) + 3/(2y^8)) / ((3/2) zeta(3)) beyond.
[[nodiscard]] double fermiDiracCharacteristicFunction(double y);
// Direct oscillatory quadrature (reference for tests; panel count grows with y).
[[nodiscard]] double fermiDiracCharacteristicFunctionQuadrature(double y);

// Superconformal-time interval s(a1) - s(a0) = int da / (a^3 H) [SI s].
[[nodiscard]] double superconformalTimeIntervalSi(const core::LambdaCdmBackground& background, double a0, double a1);

struct NeutrinoSourceSettings {
  core::NeutrinoMethodConfig method = core::NeutrinoMethodConfig::kTransferRatio;
  std::size_t k_bins = 64U;
  std::size_t history_capacity = 4096U;
  double k_min_h_mpc = 0.0;
  double k_max_h_mpc = 0.0;
  double hubble_param = 0.0;           // converts PM k [1/Mpc] -> [h/Mpc]
  double omega_baryon = 0.0;
};

class NeutrinoSpectralSource final : public gravity::PmSpectralSourceHook {
 public:
  NeutrinoSpectralSource(const core::LambdaCdmBackground& background, spectral::BoltzmannTransferSet tables,
                         NeutrinoSourceSettings settings, const parallel::MpiContext* mpi,
                         core::MemoryGovernor* governor);

  // Builds the engine from a frozen configuration (loads and verifies the
  // tables). Returns null for methods that need no PM hook.
  [[nodiscard]] static std::unique_ptr<NeutrinoSpectralSource> fromConfig(const core::SimulationConfig& config,
                                                                          const parallel::MpiContext* mpi,
                                                                          core::MemoryGovernor* governor);

  // gravity::PmSpectralSourceHook
  void beginSolve(double scale_factor) override;
  [[nodiscard]] bool needsObservation() const noexcept override;
  [[nodiscard]] std::size_t binCount() const noexcept override { return m_bin_k_h_mpc.size(); }
  [[nodiscard]] std::size_t binIndex(double k_pm) const noexcept override;
  void finalizeObservation(std::span<const double> local_bin_power, std::span<const double> local_bin_modes) override;
  void requireWavenumberDomain(double k_min_pm, double k_max_pm) const override;
  [[nodiscard]] double multiplier(double k_pm) const noexcept override;
  [[nodiscard]] std::uint64_t fingerprint() const noexcept override;

  // Current R(k) = delta_nu / delta_cb per bin (after the last solve).
  [[nodiscard]] std::span<const double> responseRatio() const noexcept { return m_current_ratio; }
  [[nodiscard]] std::span<const double> binCentersHMpc() const noexcept { return m_bin_k_h_mpc; }
  [[nodiscard]] std::size_t historySize() const noexcept { return m_history_s.size(); }
  // rho_nu(a) / rho_cb(a) from the background authority.
  [[nodiscard]] double densityFraction(double scale_factor) const;
  // Upper bound of the neglected short-range neutrino force: M(k) - 1 at k.
  [[nodiscard]] double shortRangeOmissionBound(double k_pm) const noexcept { return multiplier(k_pm) - 1.0; }

  // Persistence (whole-block module sidecar).
  [[nodiscard]] core::ModuleSidecarBlock serialize() const;
  void restore(const core::ModuleSidecarBlock& block);
  void storeTo(core::ModuleSidecarRegistry& registry) const { registry.upsert(serialize()); }
  // Restores when the sidecar exists; returns whether it did.
  bool restoreFrom(const core::ModuleSidecarRegistry& registry);

  // Pure evaluation of the linear-response ratio for one bin given an
  // explicit history (exposed for analytic tests): s nodes ascending,
  // amplitude ratios A(s_j)/A(s_now), ratios R(s_j), scale factors a_j.
  [[nodiscard]] double evaluateLinearResponse(std::size_t bin, double s_now, double a_now,
                                              std::span<const double> s_nodes, std::span<const double> a_nodes,
                                              std::span<const double> amplitude_over_now,
                                              std::span<const double> ratio_nodes) const;

 private:
  void computeTransferRatio(double scale_factor);
  void computeLinearResponse(double scale_factor, std::span<const double> amplitude);
  [[nodiscard]] double superconformalTimeSi(double scale_factor) const;

  core::LambdaCdmBackground m_background;
  spectral::BoltzmannTransferSet m_tables;
  NeutrinoSourceSettings m_settings;
  const parallel::MpiContext* m_mpi = nullptr;
  core::MemoryReservation m_reservation;
  double m_velocity_scale_m_s = 0.0;     // u_0 = k T_nu,0 / (m c)
  double m_pm_k_to_si = 0.0;             // 1/Mpc -> 1/m
  double m_initial_scale_factor = 0.0;
  std::vector<double> m_initial_delta_ratio;   // delta_nu,i / delta_cb,i per bin
  std::vector<double> m_initial_theta_ratio_s; // theta_nu,i (s-time, 1/s) / delta_cb,i per bin
  std::vector<double> m_bin_k_h_mpc;
  std::vector<double> m_bin_ln_edges;          // ln(k [h/Mpc]) edges, size bins + 1
  // History (bounded): flattened [entry * bins + bin].
  std::vector<double> m_history_s;
  std::vector<double> m_history_a;
  std::vector<double> m_history_amplitude;
  std::vector<double> m_history_ratio;
  std::vector<double> m_current_ratio;
  std::vector<double> m_current_multiplier;
  double m_solve_scale_factor = 0.0;
  bool m_solve_open = false;
};

}  // namespace cosmosim::neutrinos
