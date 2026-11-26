#pragma once

// CHUI-RT-BH-AGN-REF-v1, BH companion section 5 (BH05): unresolved scalar
// accretion-disk reservoir.
//
//   dM_d/dt = S - (1 + eta_w) mdot_in,   mdot_in = min(M_d / t_visc, mdot_cap)
//   wind     = eta_w mdot_in
//
// Uncapped (M_d <= M* = mdot_cap t_visc): exact constant-coefficient solution
//   tau = t_visc/(1+eta_w),  M' = M e^{-dt/tau} + S tau (1 - e^{-dt/tau})
// using expm1 (no cancellation at dt << tau). Capped (M_d > M*): the inner
// inflow saturates at the cap and the excess REMAINS in the disk,
//   M' = M + (S - (1+eta_w) mdot_cap) dt   (linear, exact),
// i.e. the Eddington limit changes the ODE instead of clipping after the drain.
// The integrator is piecewise analytic across the cap threshold (exact crossing
// time), so the integrated ledger D = M + S dt - M' is exact to round-off and
// M' >= 0 for S >= 0. The cap is held fixed over one call; time-varying
// coefficients (BH mass growth changes mdot_cap) are handled by the caller's
// controlled substepping (BlackHoleEngine, relative mass-change tolerance).
// The legacy instantaneous comparator has no disk (requested supply is capped
// before withdrawal, as in BlackHoleAgnModel).

#include <cstdint>
#include <limits>

namespace cosmosim::physics {

struct BhDiskStepInputs {
  double disk_mass_code = 0.0;
  double supply_rate_code = 0.0;   // accepted supply S >= 0 (constant over dt)
  double viscous_time_code = 0.0;  // > 0
  double wind_loading = 0.0;       // eta_w >= 0
  double inner_rate_cap_code = std::numeric_limits<double>::infinity();
  double dt_code = 0.0;
};

struct BhDiskStepResult {
  double disk_mass_end_code = 0.0;
  double drained_mass_code = 0.0;  // D = M + S dt - M'
  double inner_mass_code = 0.0;    // dm0 reaching the inner engine (after winds)
  double wind_mass_code = 0.0;
  double capped_time_code = 0.0;   // time spent at the inner-rate cap
  std::uint32_t pieces = 0U;
};

// 1 - exp(-x) without cancellation for small x >= 0.
[[nodiscard]] double oneMinusExpNeg(double x) noexcept;

[[nodiscard]] BhDiskStepResult integrateScalarDisk(const BhDiskStepInputs& inputs);

}  // namespace cosmosim::physics
