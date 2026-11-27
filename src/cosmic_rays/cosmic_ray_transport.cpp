// CHUI-EXT-PHYS-REF-v1 CR01-CR03 transport kernels (source only, NOT compiled
// or executed in pass 01).

#include "cosmosim/cosmic_rays/cosmic_ray_transport.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cosmosim/mhd/mhd_transport.hpp"

namespace cosmosim::cosmic_rays {
namespace {

using Ijk = std::array<long, 3>;

struct Access {
  const CrGrid& grid;
  [[nodiscard]] bool inside(const Ijk& p) const {
    for (std::size_t a = 0; a < 3U; ++a) {
      if (p[a] < 0 || p[a] >= static_cast<long>(grid.n[a])) return false;
    }
    return true;
  }
  [[nodiscard]] Ijk wrap(Ijk p) const {
    for (std::size_t a = 0; a < 3U; ++a) {
      const long n = static_cast<long>(grid.n[a]);
      p[a] = ((p[a] % n) + n) % n;
    }
    return p;
  }
  [[nodiscard]] Ijk clampIndex(Ijk p) const {
    for (std::size_t a = 0; a < 3U; ++a) p[a] = std::clamp<long>(p[a], 0, static_cast<long>(grid.n[a]) - 1);
    return p;
  }
  [[nodiscard]] std::size_t idx(const Ijk& p) const {
    return grid.index(static_cast<std::size_t>(p[0]), static_cast<std::size_t>(p[1]), static_cast<std::size_t>(p[2]));
  }
  // Scalar ghost rule: periodic wrap, outflow -> 0 (escape), reflecting -> mirror.
  [[nodiscard]] double scalar(std::span<const double> f, const Ijk& p) const {
    if (inside(p)) return f[idx(p)];
    switch (grid.boundary) {
      case CrBoundary::kPeriodic:
        return f[idx(wrap(p))];
      case CrBoundary::kOutflow:
        return 0.0;
      case CrBoundary::kReflecting:
        return f[idx(clampIndex(p))];
    }
    return 0.0;
  }
  // Medium ghost rule: periodic wrap, otherwise zero-gradient copy.
  [[nodiscard]] double medium(std::span<const double> f, const Ijk& p) const {
    if (inside(p)) return f[idx(p)];
    return grid.boundary == CrBoundary::kPeriodic ? f[idx(wrap(p))] : f[idx(clampIndex(p))];
  }
};

[[nodiscard]] std::array<double, 3> fieldAt(const Access& acc, const CrMediumView& m, const Ijk& p) {
  return {acc.medium(m.bx, p), acc.medium(m.by, p), acc.medium(m.bz, p)};
}

[[nodiscard]] std::array<double, 3> velocityAt(const Access& acc, const CrMediumView& m, const Ijk& p) {
  return {acc.medium(m.ux, p), acc.medium(m.uy, p), acc.medium(m.uz, p)};
}

void requireMedium(const CrGrid& grid, const CrMediumView& m, bool need_streaming) {
  const std::size_t n = grid.cells();
  if (m.ux.size() != n || m.uy.size() != n || m.uz.size() != n || m.bx.size() != n || m.by.size() != n ||
      m.bz.size() != n) {
    throw std::invalid_argument("CR medium view does not match the grid");
  }
  if (need_streaming && m.streaming_speed.size() != n) throw std::invalid_argument("CR streaming speed missing");
}

[[nodiscard]] bool activeAxis(const CrGrid& grid, std::size_t a) {
  return !(grid.n[a] == 1U && grid.boundary == CrBoundary::kPeriodic);
}

// Central gradient of a scalar at a cell.
[[nodiscard]] std::array<double, 3> centralGradient(const Access& acc, std::span<const double> f, const Ijk& p) {
  std::array<double, 3> g{};
  for (std::size_t a = 0; a < 3U; ++a) {
    if (!activeAxis(acc.grid, a)) continue;
    Ijk plus = p;
    Ijk minus = p;
    plus[a] += 1;
    minus[a] -= 1;
    g[a] = (acc.scalar(f, plus) - acc.scalar(f, minus)) / (2.0 * acc.grid.dx[a]);
  }
  return g;
}

// Face gradient on the + face of p along axis a: normal from the two cells,
// transverse from the MC-limited four one-sided differences.
[[nodiscard]] std::array<double, 3> faceGradient(const Access& acc, std::span<const double> f, const Ijk& p,
                                                 std::size_t a) {
  Ijk q = p;
  q[a] += 1;
  std::array<double, 3> g{};
  g[a] = (acc.scalar(f, q) - acc.scalar(f, p)) / acc.grid.dx[a];
  for (std::size_t t = 0; t < 3U; ++t) {
    if (t == a || !activeAxis(acc.grid, t)) continue;
    Ijk pp = p, pm = p, qp = q, qm = q;
    pp[t] += 1;
    pm[t] -= 1;
    qp[t] += 1;
    qm[t] -= 1;
    const double h = acc.grid.dx[t];
    g[t] = mcLimiter((acc.scalar(f, pp) - acc.scalar(f, p)) / h, (acc.scalar(f, p) - acc.scalar(f, pm)) / h,
                     (acc.scalar(f, qp) - acc.scalar(f, q)) / h, (acc.scalar(f, q) - acc.scalar(f, qm)) / h);
  }
  return g;
}

[[nodiscard]] std::array<double, 3> unitOrZero(const std::array<double, 3>& v, double threshold, bool& weak) {
  const double m = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  weak = !(m > threshold) || m == 0.0;
  if (weak) return {0.0, 0.0, 0.0};
  return {v[0] / m, v[1] / m, v[2] / m};
}

template <typename FaceFlux>
void forEachFace(const CrGrid& grid, FaceFlux&& visit) {
  // Visits every + face once; for non-periodic grids also the low boundary faces
  // (reported with `low_boundary = true`, flux oriented along +axis).
  for (std::size_t a = 0; a < 3U; ++a) {
    if (!activeAxis(grid, a)) continue;
    for (std::size_t k = 0; k < grid.n[2]; ++k) {
      for (std::size_t j = 0; j < grid.n[1]; ++j) {
        for (std::size_t i = 0; i < grid.n[0]; ++i) {
          const Ijk p{static_cast<long>(i), static_cast<long>(j), static_cast<long>(k)};
          visit(p, a, false);
          if (grid.boundary != CrBoundary::kPeriodic && p[a] == 0) {
            Ijk ghost = p;
            ghost[a] = -1;
            visit(ghost, a, true);
          }
        }
      }
    }
  }
}

// Scatter a +face flux (oriented +axis) between the cell p and its + neighbour.
void scatterFlux(const Access& acc, const Ijk& p, std::size_t a, double flux_times_dt_over_dx,
                 std::vector<double>& delta, double& outflow) {
  Ijk q = p;
  q[a] += 1;
  const bool p_in = acc.inside(p);
  const bool q_in = acc.inside(q);
  const double area_dt_factor = acc.grid.cellVolume();  // flux*dt/dx * V = flux*dt*area
  if (acc.grid.boundary == CrBoundary::kPeriodic) {
    delta[acc.idx(p)] -= flux_times_dt_over_dx;
    delta[acc.idx(acc.wrap(q))] += flux_times_dt_over_dx;
    return;
  }
  if (p_in) {
    delta[acc.idx(p)] -= flux_times_dt_over_dx;
  } else {
    outflow -= flux_times_dt_over_dx * area_dt_factor;  // flux into the domain through the low face
  }
  if (q_in) {
    delta[acc.idx(q)] += flux_times_dt_over_dx;
  } else {
    outflow += flux_times_dt_over_dx * area_dt_factor;
  }
}

}  // namespace

void CrGasExchange::resize(std::size_t n) {
  momentum_x.assign(n, 0.0);
  momentum_y.assign(n, 0.0);
  momentum_z.assign(n, 0.0);
  energy.assign(n, 0.0);
}

double mcLimiter(double a, double b, double c, double d) noexcept {
  const bool all_pos = a > 0.0 && b > 0.0 && c > 0.0 && d > 0.0;
  const bool all_neg = a < 0.0 && b < 0.0 && c < 0.0 && d < 0.0;
  if (!all_pos && !all_neg) return 0.0;
  const double sign = all_pos ? 1.0 : -1.0;
  const double m = std::min({2.0 * std::abs(a), 2.0 * std::abs(b), 2.0 * std::abs(c), 2.0 * std::abs(d),
                             std::abs(a + b + c + d) / 4.0});
  return sign * m;
}

double advectionTimestep(const CrGrid& grid, const CrMediumView& m, double cfl) {
  double rate = 0.0;
  for (std::size_t c = 0; c < grid.cells(); ++c) {
    rate = std::max(rate, std::abs(m.ux[c]) / grid.dx[0] + std::abs(m.uy[c]) / grid.dx[1] + std::abs(m.uz[c]) / grid.dx[2]);
  }
  return rate > 0.0 ? cfl / rate : 1.0e300;
}

void advectGrey(const CrGrid& grid, std::vector<double>& e, const CrMediumView& m, double gamma_cr, double dt,
                CrGasExchange& gas, CrStepLedger& ledger) {
  requireMedium(grid, m, false);
  if (e.size() != grid.cells()) throw std::invalid_argument("CR energy size mismatch");
  if (!(gamma_cr > 1.0)) throw std::invalid_argument("gamma_cr must exceed 1");
  const Access acc{grid};
  const std::size_t n = grid.cells();
  if (gas.energy.size() != n) gas.resize(n);
  std::vector<double> delta(n, 0.0);
  double outflow = 0.0;
  forEachFace(grid, [&](const Ijk& p, std::size_t a, bool) {
    Ijk q = p;
    q[a] += 1;
    const std::array<std::span<const double>, 3> u{m.ux, m.uy, m.uz};
    double u_face = 0.5 * (acc.medium(u[a], p) + acc.medium(u[a], q));
    if (grid.boundary == CrBoundary::kReflecting && (!acc.inside(p) || !acc.inside(q))) u_face = 0.0;
    const double upwind = u_face > 0.0 ? acc.scalar(e, p) : acc.scalar(e, q);
    scatterFlux(acc, p, a, u_face * upwind * dt / grid.dx[a], delta, outflow);
  });
  const double volume = grid.cellVolume();
  for (std::size_t k = 0; k < grid.n[2]; ++k) {
    for (std::size_t j = 0; j < grid.n[1]; ++j) {
      for (std::size_t i = 0; i < grid.n[0]; ++i) {
        const Ijk p{static_cast<long>(i), static_cast<long>(j), static_cast<long>(k)};
        const std::size_t c = grid.index(i, j, k);
        double div_u = 0.0;
        for (std::size_t a = 0; a < 3U; ++a) {
          if (!activeAxis(grid, a)) continue;
          Ijk plus = p;
          Ijk minus = p;
          plus[a] += 1;
          minus[a] -= 1;
          const std::array<std::span<const double>, 3> u{m.ux, m.uy, m.uz};
          div_u += (acc.medium(u[a], plus) - acc.medium(u[a], minus)) / (2.0 * grid.dx[a]);
        }
        const double pressure = (gamma_cr - 1.0) * e[c];
        const double work = -pressure * div_u * dt;  // to CR
        delta[c] += work;
        ledger.adiabatic_work += work * volume;
        // Gas exchange: force -grad P, work -u . grad P (P from the pre-step state).
        std::array<double, 3> grad_p{};
        for (std::size_t a = 0; a < 3U; ++a) {
          if (!activeAxis(grid, a)) continue;
          Ijk plus = p;
          Ijk minus = p;
          plus[a] += 1;
          minus[a] -= 1;
          grad_p[a] = (gamma_cr - 1.0) * (acc.scalar(e, plus) - acc.scalar(e, minus)) / (2.0 * grid.dx[a]);
        }
        const auto u = velocityAt(acc, m, p);
        gas.momentum_x[c] -= grad_p[0] * dt;
        gas.momentum_y[c] -= grad_p[1] * dt;
        gas.momentum_z[c] -= grad_p[2] * dt;
        const double gas_work = -(u[0] * grad_p[0] + u[1] * grad_p[1] + u[2] * grad_p[2]) * dt;
        gas.energy[c] += gas_work;
        ledger.gas_energy_change += gas_work * volume;
      }
    }
  }
  for (std::size_t c = 0; c < n; ++c) {
    e[c] += delta[c];
    ledger.cr_energy_change += delta[c] * volume;
    if (e[c] < 0.0) throw std::runtime_error("CR advection produced negative energy (reduce dt)");
  }
  ledger.boundary_flux_out += outflow;
}

void anisotropicDiffusionRate(const CrGrid& grid, std::span<const double> e, const CrMediumView& m,
                              const CrDiffusionModel& model, std::vector<double>& rate, double* boundary_out_rate) {
  requireMedium(grid, m, false);
  const Access acc{grid};
  rate.assign(grid.cells(), 0.0);
  double outflow = 0.0;
  forEachFace(grid, [&](const Ijk& p, std::size_t a, bool) {
    Ijk q = p;
    q[a] += 1;
    const auto bp = fieldAt(acc, m, p);
    const auto bq = fieldAt(acc, m, q);
    const std::array<double, 3> bf{0.5 * (bp[0] + bq[0]), 0.5 * (bp[1] + bq[1]), 0.5 * (bp[2] + bq[2])};
    bool weak = false;
    const auto b = unitOrZero(bf, model.weak_field_threshold, weak);
    const auto g = faceGradient(acc, e, p, a);
    double flux = 0.0;
    if (weak) {
      flux = -model.kappa_isotropic_weak * g[a];
    } else {
      const double bg = b[0] * g[0] + b[1] * g[1] + b[2] * g[2];
      flux = -model.kappa_parallel * b[a] * bg - model.kappa_perpendicular * (g[a] - b[a] * bg);
    }
    if (grid.boundary == CrBoundary::kReflecting && (!acc.inside(p) || !acc.inside(q))) flux = 0.0;
    scatterFlux(acc, p, a, flux / grid.dx[a], rate, outflow);
  });
  if (boundary_out_rate != nullptr) *boundary_out_rate = outflow;
}

double diffusionExplicitTimestep(const CrGrid& grid, const CrDiffusionModel& model, double cfl) {
  const double kappa = std::max({model.kappa_parallel, model.kappa_perpendicular, model.kappa_isotropic_weak});
  if (!(kappa > 0.0)) return 1.0e300;
  double inv = 0.0;
  for (std::size_t a = 0; a < 3U; ++a) {
    if (activeAxis(grid, a)) inv += 1.0 / (grid.dx[a] * grid.dx[a]);
  }
  return cfl / (2.0 * kappa * inv);
}

void diffuseExplicit(const CrGrid& grid, std::vector<double>& e, const CrMediumView& m, const CrDiffusionModel& model,
                     double dt, CrStepLedger& ledger) {
  std::vector<double> rate;
  double out_rate = 0.0;
  anisotropicDiffusionRate(grid, e, m, model, rate, &out_rate);
  const double volume = grid.cellVolume();
  for (std::size_t c = 0; c < e.size(); ++c) {
    e[c] += dt * rate[c];
    ledger.cr_energy_change += dt * rate[c] * volume;
  }
  ledger.boundary_flux_out += dt * out_rate;
}

std::int32_t diffuseRkl2(const CrGrid& grid, std::vector<double>& e, const CrMediumView& m,
                         const CrDiffusionModel& model, double dt, double cfl, CrStepLedger& ledger) {
  const double dt_explicit = diffusionExplicitTimestep(grid, model, cfl);
  const double ratio = dt / dt_explicit;
  if (ratio <= 1.0) {
    diffuseExplicit(grid, e, m, model, dt, ledger);
    return 1;
  }
  const std::int32_t stages = mhd::selectRkl2StageCount(ratio, 4096);
  if (stages <= 0 || stages > 4096) {
    throw std::runtime_error("CR RKL2: required stage count exceeds 4096 (split the interval)");
  }
  const mhd::MhdRkl2Coefficients co = mhd::makeRkl2Coefficients(stages);
  const std::size_t n = e.size();
  std::vector<double> y0 = e;
  std::vector<double> l0;
  double out0 = 0.0;
  anisotropicDiffusionRate(grid, y0, m, model, l0, &out0);
  std::vector<double> ym2 = y0;
  std::vector<double> ym1(n);
  for (std::size_t c = 0; c < n; ++c) ym1[c] = y0[c] + co.mu_tilde[1] * dt * l0[c];
  std::vector<double> lj;
  for (std::int32_t j = 2; j <= stages; ++j) {
    const auto jj = static_cast<std::size_t>(j);
    double out = 0.0;
    anisotropicDiffusionRate(grid, ym1, m, model, lj, &out);
    std::vector<double> yj(n);
    for (std::size_t c = 0; c < n; ++c) {
      yj[c] = co.mu[jj] * ym1[c] + co.nu[jj] * ym2[c] + (1.0 - co.mu[jj] - co.nu[jj]) * y0[c] +
              co.mu_tilde[jj] * dt * lj[c] + co.gamma_tilde[jj] * dt * l0[c];
    }
    ym2 = std::move(ym1);
    ym1 = std::move(yj);
  }
  const double volume = grid.cellVolume();
  double before = 0.0;
  double after = 0.0;
  for (std::size_t c = 0; c < n; ++c) {
    before += e[c];
    after += ym1[c];
  }
  e = std::move(ym1);
  ledger.cr_energy_change += (after - before) * volume;
  // RKL2 conserves the discrete total exactly for periodic/reflecting grids;
  // with outflow faces the loss is the total change.
  ledger.boundary_flux_out += grid.boundary == CrBoundary::kOutflow ? -(after - before) * volume : 0.0;
  return stages;
}

double streamingTimestep(const CrGrid& grid, const CrMediumView& m, double l_reg, double cfl) {
  double vmax = 0.0;
  for (const double v : m.streaming_speed) vmax = std::max(vmax, v);
  if (!(vmax > 0.0)) return 1.0e300;
  double dxmin = 1.0e300;
  for (std::size_t a = 0; a < 3U; ++a) {
    if (activeAxis(grid, a)) dxmin = std::min(dxmin, grid.dx[a]);
  }
  return cfl * dxmin / vmax * std::min(1.0, dxmin / std::max(l_reg, 1.0e-300));
}

void streamGrey(const CrGrid& grid, std::vector<double>& e, const CrMediumView& m, double gamma_cr, double l_reg,
                double dt, CrGasExchange& gas, CrStepLedger& ledger) {
  requireMedium(grid, m, true);
  if (!(l_reg > 0.0)) throw std::invalid_argument("streaming regularization length must be > 0");
  const Access acc{grid};
  const std::size_t n = grid.cells();
  if (gas.energy.size() != n) gas.resize(n);
  const double tiny = 1.0e-300;
  std::vector<double> delta(n, 0.0);
  double outflow = 0.0;
  forEachFace(grid, [&](const Ijk& p, std::size_t a, bool) {
    Ijk q = p;
    q[a] += 1;
    const auto bp = fieldAt(acc, m, p);
    const auto bq = fieldAt(acc, m, q);
    bool weak = false;
    const auto b = unitOrZero({0.5 * (bp[0] + bq[0]), 0.5 * (bp[1] + bq[1]), 0.5 * (bp[2] + bq[2])}, 0.0, weak);
    if (weak) return;
    const double e_face = 0.5 * (acc.scalar(e, p) + acc.scalar(e, q));
    if (!(e_face > tiny)) return;
    const auto g = faceGradient(acc, e, p, a);
    const double bg = b[0] * g[0] + b[1] * g[1] + b[2] * g[2];
    const double vs = 0.5 * (acc.medium(m.streaming_speed, p) + acc.medium(m.streaming_speed, q));
    const double v_st_a = -vs * b[a] * std::tanh(l_reg * bg / e_face);
    double flux = gamma_cr * e_face * v_st_a;  // (e + P) v_st with P = (gamma - 1) e
    if (grid.boundary == CrBoundary::kReflecting && (!acc.inside(p) || !acc.inside(q))) flux = 0.0;
    scatterFlux(acc, p, a, flux * dt / grid.dx[a], delta, outflow);
  });
  const double volume = grid.cellVolume();
  for (std::size_t k = 0; k < grid.n[2]; ++k) {
    for (std::size_t j = 0; j < grid.n[1]; ++j) {
      for (std::size_t i = 0; i < grid.n[0]; ++i) {
        const Ijk p{static_cast<long>(i), static_cast<long>(j), static_cast<long>(k)};
        const std::size_t c = grid.index(i, j, k);
        bool weak = false;
        const auto b = unitOrZero(fieldAt(acc, m, p), 0.0, weak);
        if (weak || !(e[c] > tiny)) continue;
        const auto g = centralGradient(acc, e, p);
        const double bg = b[0] * g[0] + b[1] * g[1] + b[2] * g[2];
        // H = -v_st . grad P = v_s tanh(l b.grad e / e) (gamma - 1) b.grad e >= 0.
        const double heating = m.streaming_speed[c] * std::tanh(l_reg * bg / e[c]) * (gamma_cr - 1.0) * bg * dt;
        delta[c] -= heating;
        gas.energy[c] += heating;
        ledger.streaming_heating += heating * volume;
        ledger.gas_energy_change += heating * volume;
      }
    }
  }
  for (std::size_t c = 0; c < n; ++c) {
    e[c] += delta[c];
    ledger.cr_energy_change += delta[c] * volume;
    if (e[c] < 0.0) throw std::runtime_error("CR streaming produced negative energy (reduce dt)");
  }
  ledger.boundary_flux_out += outflow;
}

bool advanceTwoMoment(const CrGrid& grid, CrTwoMomentState& state, const CrMediumView& m,
                      const CrTwoMomentSettings& s, double dt, CrTwoMomentReport& report, CrGasExchange* gas) {
  const bool streaming = !m.streaming_speed.empty();
  requireMedium(grid, m, streaming);
  const std::size_t n = grid.cells();
  if (state.e.size() != n || state.f.size() != n) throw std::invalid_argument("two-moment state size mismatch");
  if (!(s.reduced_speed > 0.0)) throw std::invalid_argument("two-moment CR needs an explicit reduced speed c_r > 0");
  if (grid.boundary == CrBoundary::kOutflow) {
    throw std::invalid_argument("two-moment CR: outflow boundaries are not supported (periodic or reflecting)");
  }
  const Access acc{grid};
  double rate = 0.0;
  for (std::size_t a = 0; a < 3U; ++a) {
    if (activeAxis(grid, a)) rate += s.reduced_speed / grid.dx[a];
  }
  const std::size_t substeps = std::max<std::size_t>(1U, static_cast<std::size_t>(std::ceil(dt * rate / s.cfl)));
  const double h = dt / static_cast<double>(substeps);
  const double cr2 = s.reduced_speed * s.reduced_speed;
  const double volume = grid.cellVolume();
  CrTwoMomentState work = state;
  CrGasExchange exchange;
  exchange.resize(n);
  double gas_gain = 0.0;
  double cr_change = 0.0;
  for (std::size_t step = 0; step < substeps; ++step) {
    std::vector<double> transfer(n, 0.0);  // energy moved CR -> gas per cell (per volume)
    std::array<std::vector<double>, 3> impulse{std::vector<double>(n, 0.0), std::vector<double>(n, 0.0),
                                               std::vector<double>(n, 0.0)};
    // 1. Face fluxes: exact relaxation toward the frozen target.
    for (std::size_t a = 0; a < 3U; ++a) {
      if (!activeAxis(grid, a)) continue;
      const std::array<std::span<const double>, 3> u_lanes{m.ux, m.uy, m.uz};
      for (std::size_t k = 0; k < grid.n[2]; ++k) {
        for (std::size_t j = 0; j < grid.n[1]; ++j) {
          for (std::size_t i = 0; i < grid.n[0]; ++i) {
            const Ijk p{static_cast<long>(i), static_cast<long>(j), static_cast<long>(k)};
            Ijk q = p;
            q[a] += 1;
            const std::size_t c = grid.index(i, j, k);
            if (!acc.inside(q) && grid.boundary == CrBoundary::kReflecting) {
              work.f[c][a] = 0.0;  // wall: zero normal flux
              continue;
            }
            const auto bp = fieldAt(acc, m, p);
            const auto bq = fieldAt(acc, m, q);
            bool weak = false;
            const auto b = unitOrZero({0.5 * (bp[0] + bq[0]), 0.5 * (bp[1] + bq[1]), 0.5 * (bp[2] + bq[2])},
                                      s.diffusion.weak_field_threshold, weak);
            const auto g = faceGradient(acc, work.e, p, a);
            double diffusive = 0.0;
            double kappa_n = 0.0;
            if (weak) {
              kappa_n = s.diffusion.kappa_isotropic_weak;
              diffusive = -kappa_n * g[a];
            } else {
              const double bg = b[0] * g[0] + b[1] * g[1] + b[2] * g[2];
              diffusive = -s.diffusion.kappa_parallel * b[a] * bg - s.diffusion.kappa_perpendicular * (g[a] - b[a] * bg);
              kappa_n = b[a] * b[a] * s.diffusion.kappa_parallel + (1.0 - b[a] * b[a]) * s.diffusion.kappa_perpendicular;
            }
            const double u_face = 0.5 * (acc.medium(u_lanes[a], p) + acc.medium(u_lanes[a], q));
            double v_st = 0.0;
            if (streaming && !weak) {
              const double bg = b[0] * g[0] + b[1] * g[1] + b[2] * g[2];
              const double vs = 0.5 * (acc.medium(m.streaming_speed, p) + acc.medium(m.streaming_speed, q));
              v_st = -vs * b[a] * (bg > 0.0 ? 1.0 : (bg < 0.0 ? -1.0 : 0.0));
            }
            const double v_total = u_face + v_st;
            const double e_up = v_total > 0.0 ? acc.scalar(work.e, p) : acc.scalar(work.e, q);
            const double f_eq = s.gamma_cr * e_up * v_total;
            const double f_target = f_eq + diffusive;
            const double f_pre = work.f[c][a];
            double f_post = f_target;
            double face_impulse = 0.0;
            if (kappa_n > 0.0) {
              const double x = h * cr2 / (3.0 * kappa_n);
              const double decay = std::exp(-x);
              f_post = f_target + (f_pre - f_target) * decay;
              const double avg_factor = x > 1.0e-12 ? -std::expm1(-x) / x : 1.0 - 0.5 * x;
              const double f_avg = f_target + (f_pre - f_target) * avg_factor;
              face_impulse = h * (f_avg - f_eq) / (3.0 * kappa_n);
            } else {
              // Perfect coupling: the normal force is the CR pressure gradient.
              face_impulse = -h * (s.gamma_cr - 1.0) * g[a];
            }
            work.f[c][a] = f_post;
            const double energy_moved = v_total * face_impulse;
            const std::size_t cq = acc.idx(acc.inside(q) ? q : acc.wrap(q));
            transfer[c] += 0.5 * energy_moved;
            transfer[cq] += 0.5 * energy_moved;
            impulse[a][c] += 0.5 * face_impulse;
            impulse[a][cq] += 0.5 * face_impulse;
            report.max_speed_ratio = std::max(report.max_speed_ratio, std::abs(v_total) / s.reduced_speed);
          }
        }
      }
    }
    // 2. Energy update with the new face fluxes and the gas transfer.
    std::vector<double> e_new = work.e;
    for (std::size_t k = 0; k < grid.n[2]; ++k) {
      for (std::size_t j = 0; j < grid.n[1]; ++j) {
        for (std::size_t i = 0; i < grid.n[0]; ++i) {
          const Ijk p{static_cast<long>(i), static_cast<long>(j), static_cast<long>(k)};
          const std::size_t c = grid.index(i, j, k);
          double div = 0.0;
          for (std::size_t a = 0; a < 3U; ++a) {
            if (!activeAxis(grid, a)) continue;
            Ijk lo = p;
            lo[a] -= 1;
            double f_low = 0.0;
            if (acc.inside(lo)) {
              f_low = work.f[acc.idx(lo)][a];
            } else if (grid.boundary == CrBoundary::kPeriodic) {
              f_low = work.f[acc.idx(acc.wrap(lo))][a];
            }  // reflecting wall: 0
            div += (work.f[c][a] - f_low) / grid.dx[a];
          }
          e_new[c] -= h * div + transfer[c];
          if (e_new[c] < 0.0) return false;  // positivity refusal: state untouched
        }
      }
    }
    for (std::size_t c = 0; c < n; ++c) {
      cr_change += (e_new[c] - work.e[c]) * volume;
      gas_gain += transfer[c] * volume;
      exchange.energy[c] += transfer[c];
      exchange.momentum_x[c] += impulse[0][c];
      exchange.momentum_y[c] += impulse[1][c];
      exchange.momentum_z[c] += impulse[2][c];
    }
    work.e = std::move(e_new);
    ++report.substeps;
  }
  state = std::move(work);
  report.gas_energy_gain += gas_gain;
  report.cr_energy_change += cr_change;
  if (gas != nullptr) *gas = std::move(exchange);
  return true;
}

}  // namespace cosmosim::cosmic_rays
