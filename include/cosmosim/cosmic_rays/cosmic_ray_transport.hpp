#pragma once

// CHUI-EXT-PHYS-REF-v1 CR01-CR03, CR07: cosmic-ray transport kernels on one
// uniform Cartesian patch. Lazy state: nothing is allocated unless a CR
// model is constructed (no CR memory in DMO or CR-disabled runs).
//
// Conventions. Energy densities e [any consistent unit], velocities u and
// streaming speeds v_st [L/T], diffusion coefficients kappa [L^2/T]. B is any
// cell-centred field sample of the CT authority (only its direction and
// magnitude enter; b = B/|B|). Grey closure: P = (gamma_cr - 1) e with the
// user-declared adiabatic index (4/3 for the ultrarelativistic grey
// approximation; spectral populations compute P from the spectrum instead).
//
// CR01 grey energy equation (one-moment, Pfrommer et al. 2017 form; recalled lead):
//   d e/dt + div(e u) = -P div u - div F_diff - div F_st + v_st . grad P + Q - Lambda,
// gas exchange: momentum -grad P, energy -u . grad P (work) + H_st,
// H_st = -v_st . grad P >= 0 (streaming heating), so CR + gas energy is
// conserved up to boundary fluxes, sources and losses (explicit ledger).
// Discretization: first-order donor-cell advection of e with face velocity
// (u_L + u_R)/2 (conservative); P div u and u . grad P with central
// differences on the same cells (the pair sums to -div(P u)'s discrete
// analogue only to truncation error — the ledger reports the residual).
//
// CR02 anisotropic diffusion: F = -kappa_par b (b . grad e) - kappa_perp (grad e - b (b . grad e)).
// Face fluxes use the normal gradient and transverse gradients limited with
// the monotonized-central limiter of the four adjacent central differences
// (Sharma & Hammett 2007 "limited" scheme; recalled lead) so no new extrema
// appear for any field obliquity. Weak field: when |B| < weak_field_threshold
// the declared weak-field model applies (isotropic kappa_iso) — never an
// undefined direction. Explicit stability: dt <= cfl dx_min^2 / (2 d kappa_max);
// RKL2 super-time-stepping (mhd::makeRkl2Coefficients, Meyer, Balsara & Aslam 2014)
// covers larger steps.
// Streaming (regularized, Sharma, Colella & Martin 2010 form; recalled lead):
//   v_st = -v_s b tanh(l_reg (b . grad e) / e_floor_safe) [smooth sign],
//   F_st = (e + P) v_st; explicit with dt <= cfl dx / v_s * min(1, dx / l_reg).
//
// CR03 two-moment (DERIVED contract, not a copy of photon M1):
//   d e/dt + div F = -(u + v_st) . f_gas,
//   d F/dt + (c_r^2 / 3) grad e = -(c_r^2 / 3) K^-1 (F - F_eq),
//   F_eq = (e + P)(u + v_st),  K = kappa_par b b + kappa_perp (I - b b),
//   f_gas = (1/3) K^-1 (F - F_eq)  (force density on the gas),
// i.e. dF/dt = -(F - F_target) / tau with F_target = F_eq - K grad e and
// tau = 3 K / c_r^2. Its relaxed limit is F = F_eq - K grad e and the energy
// equation reduces to the one-moment CR01 form. Discretization: STAGGERED
// face-normal fluxes (state f[c][a] = normal flux on the + face of cell c
// along axis a); per face the target uses the CR02 Sharma-Hammett limited
// diffusive flux and a donor-cell enthalpy flux, tau uses the normal
// diffusivity kappa_n = b_a^2 kappa_par + (1 - b_a^2) kappa_perp, and the
// relaxation F <- F_target + (F - F_target) exp(-dt / tau) is exact for the
// frozen target (kappa_n -> 0: instantaneous; kappa_n -> inf: the wave term
// (c_r^2/3) grad e survives through F_target / tau). The stiff limit is the
// CR02 explicit diffusion stencil with no numerical diffusion from the wave
// speed (asymptotic preserving in space); the CFL dt <= cfl / sum c_r/dx is
// sufficient in both limits because c_r dx / kappa >> 1 (stiff) implies the
// explicit diffusion bound exceeds the wave bound. c_r is an explicit
// reduced signal speed (validity: c_r >> |u| + v_s, c_r >> kappa / L). The
// gas impulse per face is the time integral of f_gas over the substep with
// the exact exponential average of F; the transferred energy
// (u + v_st)_face . impulse is split equally between the two cells (CR loses,
// gas gains). Boundaries: periodic or reflecting (zero normal flux);
// outflow is refused. Positivity: a step producing e < 0 is refused
// (caller retries; state untouched).

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace cosmosim::cosmic_rays {

enum class CrBoundary : std::uint8_t { kPeriodic = 0, kOutflow = 1, kReflecting = 2 };

struct CrGrid {
  std::array<std::size_t, 3> n{1, 1, 1};
  std::array<double, 3> dx{1.0, 1.0, 1.0};
  CrBoundary boundary = CrBoundary::kPeriodic;
  [[nodiscard]] std::size_t cells() const noexcept { return n[0] * n[1] * n[2]; }
  [[nodiscard]] std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return i + n[0] * (j + n[1] * k);
  }
  [[nodiscard]] double cellVolume() const noexcept { return dx[0] * dx[1] * dx[2]; }
};

struct CrMediumView {
  std::span<const double> ux, uy, uz;     // gas velocity
  std::span<const double> bx, by, bz;     // cell-centred magnetic field sample
  std::span<const double> streaming_speed;  // v_s >= 0 (e.g. Alfven speed); empty = no streaming
};

struct CrDiffusionModel {
  double kappa_parallel = 0.0;
  double kappa_perpendicular = 0.0;
  double weak_field_threshold = 0.0;   // |B| below -> isotropic kappa_isotropic_weak
  double kappa_isotropic_weak = 0.0;
};

struct CrGasExchange {
  std::vector<double> momentum_x, momentum_y, momentum_z;  // force density x dt (impulse per volume)
  std::vector<double> energy;                              // energy per volume deposited in gas
  void resize(std::size_t n);
};

struct CrStepLedger {
  double cr_energy_change = 0.0;       // volume integrated
  double gas_energy_change = 0.0;      // from work + streaming heating
  double boundary_flux_out = 0.0;      // energy leaving through open faces
  double streaming_heating = 0.0;
  double adiabatic_work = 0.0;         // -int P div u dV dt (to CR)
  [[nodiscard]] double residual() const noexcept {
    return cr_energy_change + gas_energy_change + boundary_flux_out;
  }
};

[[nodiscard]] double mcLimiter(double a, double b, double c, double d) noexcept;

// CR01: advection + adiabatic work + gas exchange (explicit, first order).
void advectGrey(const CrGrid& grid, std::vector<double>& e, const CrMediumView& medium, double gamma_cr, double dt,
                CrGasExchange& gas, CrStepLedger& ledger);
[[nodiscard]] double advectionTimestep(const CrGrid& grid, const CrMediumView& medium, double cfl);

// CR02: anisotropic diffusion operator L(e) = -div F (per unit volume rate).
void anisotropicDiffusionRate(const CrGrid& grid, std::span<const double> e, const CrMediumView& medium,
                              const CrDiffusionModel& model, std::vector<double>& rate,
                              double* boundary_flux_out_rate = nullptr);
[[nodiscard]] double diffusionExplicitTimestep(const CrGrid& grid, const CrDiffusionModel& model, double cfl);
void diffuseExplicit(const CrGrid& grid, std::vector<double>& e, const CrMediumView& medium,
                     const CrDiffusionModel& model, double dt, CrStepLedger& ledger);
// RKL2 super-step of length dt with the minimal stage count for stability.
std::int32_t diffuseRkl2(const CrGrid& grid, std::vector<double>& e, const CrMediumView& medium,
                         const CrDiffusionModel& model, double dt, double cfl, CrStepLedger& ledger);

// Regularized streaming with gas heating (exactly balanced in the ledger).
void streamGrey(const CrGrid& grid, std::vector<double>& e, const CrMediumView& medium, double gamma_cr,
                double regularization_length, double dt, CrGasExchange& gas, CrStepLedger& ledger);
[[nodiscard]] double streamingTimestep(const CrGrid& grid, const CrMediumView& medium, double regularization_length,
                                       double cfl);

// CR03 two-moment state and step.
struct CrTwoMomentState {
  std::vector<double> e;
  std::vector<std::array<double, 3>> f;  // face-normal fluxes on the + faces of each cell
};
struct CrTwoMomentSettings {
  double reduced_speed = 0.0;   // c_r (explicit, validity-bounded)
  double gamma_cr = 4.0 / 3.0;
  double cfl = 0.4;
  CrDiffusionModel diffusion{};
};
struct CrTwoMomentReport {
  std::size_t substeps = 0;
  double max_speed_ratio = 0.0;  // max(|u| + v_s) / c_r (validity diagnostic)
  double gas_energy_gain = 0.0;  // volume integrated transfer to the gas
  double cr_energy_change = 0.0;
};
// Gas impulse/energy of the last successful call (per cell, per volume).
// Passing a null exchange skips the bookkeeping.
// Returns false (state untouched) if a substep would produce e < 0.
[[nodiscard]] bool advanceTwoMoment(const CrGrid& grid, CrTwoMomentState& state, const CrMediumView& medium,
                                    const CrTwoMomentSettings& settings, double dt, CrTwoMomentReport& report,
                                    CrGasExchange* gas = nullptr);

}  // namespace cosmosim::cosmic_rays
