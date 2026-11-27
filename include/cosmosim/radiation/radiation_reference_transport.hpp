#pragma once

// CHUI-EXT-PHYS-REF-v1 R03: grey reference transport backends on one uniform
// Cartesian patch (static medium, lab frame = fluid frame; frame terms
// O(v/c) are a declared omission). These are DISTINCT closures/backends, not
// variants of the M1 Riemann path, and they are NOT dispatched by
// RadiationRuntime in pass 01 (`physics.rt_transport = vet` stays refused at
// admission until AMR records, multigroup coupling and the runtime adapter
// exist). Units are any consistent set (length L, time T, c in L/T,
// opacities chi in 1/L, emissivity eta in energy / (L^3 T sr)).
//
// 1. Discrete ordinates formal solver (independent angular reference):
//    first-order upwind "step" discretization of n . grad I + chi_t I = eta_tot,
//      I_c = (eta_tot + sum_a |n_a| I_up,a / dx_a) / (chi_t + sum_a |n_a| / dx_a),
//    swept in upwind order per direction (positivity preserving). Isotropic
//    scattering by source iteration eta_tot = eta + chi_s J with the mean
//    intensity J = (1/4 pi) int I dOmega, stopped on a relative residual. Periodic boundaries iterate the inflow to
//    convergence; open boundaries impose vacuum (I_in = 0) inflow.
//    Quadrature: N equal-weight Fibonacci-sphere directions (w = 4 pi / N).
// 2. VET: normalized Eddington tensor f_ij = sum w I n_i n_j / sum w I per
//    cell from (1); realizable by construction (convex combination of
//    projectors n n^T: PSD, trace 1). Moment transport
//      dE/dt + div F = c chi_a (aT^4 - E),  dF/dt + c^2 div (f E) = -c chi_t F
//    with first-order local Lax-Friedrichs fluxes at signal speed c (stable
//    because the eigenvalues of f along any axis lie in [0, 1]) and the
//    stiff sources treated implicitly (backward Euler). The closure is
//    refreshed every `closure_cadence` transport steps from a steady formal
//    solution (time-independent VET closure, declared). A post-update
//    realizability projection |F| <= c E is counted, never silent.
// 3. FLD comparator: dE/dt = div(D grad E) + c chi_a (aT^4 - E) + S,
//    D = c lambda(R) / chi_t, R = |grad E| / (chi_t E), Levermore-Pomraning
//    limiter lambda = (2 + R) / (6 + 3R + R^2) (lambda -> 1/3 thick,
//    -> 1/R free streaming). Backward Euler with face D (harmonic mean),
//    lagged in R; SPD system solved by Jacobi-preconditioned CG.
//    Boundaries: periodic, reflecting (zero flux) or vacuum (E_ghost = 0).

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace cosmosim::radiation {

enum class ReferenceBoundary : std::uint8_t { kPeriodic = 0, kOpenVacuum = 1, kReflecting = 2 };

struct UniformPatch {
  std::array<std::size_t, 3> n{1, 1, 1};
  std::array<double, 3> dx{1.0, 1.0, 1.0};
  ReferenceBoundary boundary = ReferenceBoundary::kPeriodic;
  [[nodiscard]] std::size_t cells() const noexcept { return n[0] * n[1] * n[2]; }
  [[nodiscard]] std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return i + n[0] * (j + n[1] * k);
  }
};

[[nodiscard]] std::vector<std::array<double, 3>> fibonacciDirections(std::size_t count);

struct OrdinatesSettings {
  std::size_t directions = 96;
  std::size_t max_sweeps = 200;        // periodic inflow / source iterations
  double relative_tolerance = 1.0e-10;
};

struct OrdinatesSolution {
  std::vector<std::array<double, 3>> directions;
  std::vector<double> intensity;      // [direction * cells + cell]
  std::vector<double> mean_intensity;  // J = (1/4 pi) sum w I
  std::size_t sweeps = 0;
  double residual = 0.0;
  bool converged = false;
};

// eta: emissivity per cell (isotropic, true emission only); chi_a, chi_s per cell.
[[nodiscard]] OrdinatesSolution solveDiscreteOrdinates(const UniformPatch& patch, const std::vector<double>& eta,
                                                       const std::vector<double>& chi_a,
                                                       const std::vector<double>& chi_s,
                                                       const OrdinatesSettings& settings);

// f_ij per cell (row-major 9 values); cells with zero intensity get the isotropic 1/3 I.
[[nodiscard]] std::vector<std::array<double, 9>> eddingtonTensors(const UniformPatch& patch,
                                                                  const OrdinatesSolution& solution);

struct VetState {
  std::vector<double> energy;
  std::vector<std::array<double, 3>> flux;
  std::vector<std::array<double, 9>> eddington;
};

struct VetSettings {
  double c = 1.0;
  double cfl = 0.4;
  double radiation_constant = 1.0;
  std::size_t closure_cadence = 1;
  OrdinatesSettings ordinates{};
};

struct VetStepReport {
  std::size_t substeps = 0;
  std::size_t closure_updates = 0;
  std::size_t realizability_projections = 0;
  double absorbed_energy = 0.0;      // volume-integrated
  double emitted_energy = 0.0;
  double boundary_outflow_energy = 0.0;
};

class VetTransport {
 public:
  VetTransport(UniformPatch patch, VetSettings settings);
  // Advances (E, F) by dt with medium (chi_a, chi_s, temperature T); the
  // closure is rebuilt from the formal solution of the current emission.
  [[nodiscard]] VetStepReport advance(VetState& state, const std::vector<double>& chi_a, const std::vector<double>& chi_s,
                                      const std::vector<double>& temperature, double dt);
  [[nodiscard]] double maxStableTimestep() const noexcept;
  [[nodiscard]] const UniformPatch& patch() const noexcept { return m_patch; }

 private:
  void refreshClosure(VetState& state, const std::vector<double>& chi_a, const std::vector<double>& chi_s,
                      const std::vector<double>& temperature);
  UniformPatch m_patch;
  VetSettings m_settings;
  std::size_t m_steps_since_closure = 0;
};

struct FldSettings {
  double c = 1.0;
  double radiation_constant = 1.0;
  std::size_t max_cg_iterations = 10000;
  double cg_relative_tolerance = 1.0e-12;
};

struct FldStepReport {
  std::size_t cg_iterations = 0;
  double cg_residual = 0.0;
  double absorbed_energy = 0.0;
  double source_energy = 0.0;
};

[[nodiscard]] double levermorePomraningLimiter(double r) noexcept;

// Advances E by one backward-Euler step dt (volume-integrated energies in the report).
[[nodiscard]] FldStepReport advanceFldImplicit(const UniformPatch& patch, std::vector<double>& energy,
                                               const std::vector<double>& chi_a, const std::vector<double>& chi_s,
                                               const std::vector<double>& temperature,
                                               const std::vector<double>& source_rate, double dt,
                                               const FldSettings& settings);

}  // namespace cosmosim::radiation
