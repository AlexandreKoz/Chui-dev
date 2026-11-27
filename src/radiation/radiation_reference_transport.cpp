// CHUI-EXT-PHYS-REF-v1 R03 reference transport backends (source only, NOT
// compiled or executed in pass 01).

#include "cosmosim/radiation/radiation_reference_transport.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace cosmosim::radiation {
namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

[[nodiscard]] std::size_t wrapIndex(long i, std::size_t n) {
  const long m = static_cast<long>(n);
  return static_cast<std::size_t>(((i % m) + m) % m);
}

void requireField(const UniformPatch& patch, const std::vector<double>& field, const char* name) {
  if (field.size() != patch.cells()) throw std::invalid_argument(std::string(name) + " has the wrong cell count");
  for (const double v : field) {
    if (!std::isfinite(v) || v < 0.0) throw std::invalid_argument(std::string(name) + " must be finite and >= 0");
  }
}

}  // namespace

std::vector<std::array<double, 3>> fibonacciDirections(std::size_t count) {
  if (count < 4U) throw std::invalid_argument("at least 4 directions are required");
  std::vector<std::array<double, 3>> dirs(count);
  const double golden = k_pi * (3.0 - std::sqrt(5.0));
  for (std::size_t i = 0; i < count; ++i) {
    const double z = 1.0 - (2.0 * static_cast<double>(i) + 1.0) / static_cast<double>(count);
    const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    const double phi = golden * static_cast<double>(i);
    dirs[i] = {r * std::cos(phi), r * std::sin(phi), z};
  }
  return dirs;
}

OrdinatesSolution solveDiscreteOrdinates(const UniformPatch& patch, const std::vector<double>& eta,
                                         const std::vector<double>& chi_a, const std::vector<double>& chi_s,
                                         const OrdinatesSettings& settings) {
  requireField(patch, eta, "eta");
  requireField(patch, chi_a, "chi_a");
  requireField(patch, chi_s, "chi_s");
  if (patch.boundary == ReferenceBoundary::kReflecting) {
    throw std::invalid_argument("discrete ordinates: reflecting boundaries need mirrored ordinates (not implemented)");
  }
  OrdinatesSolution out;
  out.directions = fibonacciDirections(settings.directions);
  const std::size_t nd = out.directions.size();
  const std::size_t nc = patch.cells();
  const double weight = 4.0 * k_pi / static_cast<double>(nd);
  out.intensity.assign(nd * nc, 0.0);
  out.mean_intensity.assign(nc, 0.0);
  std::vector<double> source(nc);
  const auto& n = patch.n;
  for (std::size_t sweep = 0; sweep < settings.max_sweeps; ++sweep) {
    for (std::size_t c = 0; c < nc; ++c) source[c] = eta[c] + chi_s[c] * out.mean_intensity[c];
    double max_change = 0.0;
    double max_value = 0.0;
    for (std::size_t d = 0; d < nd; ++d) {
      const auto& dir = out.directions[d];
      double* intensity = out.intensity.data() + d * nc;
      const std::array<long, 3> step{dir[0] >= 0.0 ? 1L : -1L, dir[1] >= 0.0 ? 1L : -1L, dir[2] >= 0.0 ? 1L : -1L};
      const std::array<double, 3> w{std::abs(dir[0]) / patch.dx[0], std::abs(dir[1]) / patch.dx[1],
                                    std::abs(dir[2]) / patch.dx[2]};
      for (std::size_t kk = 0; kk < n[2]; ++kk) {
        const long k = step[2] > 0 ? static_cast<long>(kk) : static_cast<long>(n[2] - 1U - kk);
        for (std::size_t jj = 0; jj < n[1]; ++jj) {
          const long j = step[1] > 0 ? static_cast<long>(jj) : static_cast<long>(n[1] - 1U - jj);
          for (std::size_t ii = 0; ii < n[0]; ++ii) {
            const long i = step[0] > 0 ? static_cast<long>(ii) : static_cast<long>(n[0] - 1U - ii);
            const std::array<long, 3> ijk{i, j, k};
            double inflow = 0.0;
            for (std::size_t a = 0; a < 3U; ++a) {
              if (w[a] == 0.0) continue;
              std::array<long, 3> up = ijk;
              up[a] -= step[a];
              double value = 0.0;
              if (up[a] < 0 || up[a] >= static_cast<long>(n[a])) {
                if (patch.boundary == ReferenceBoundary::kPeriodic) {
                  up[a] = static_cast<long>(wrapIndex(up[a], n[a]));
                  value = intensity[patch.index(static_cast<std::size_t>(up[0]), static_cast<std::size_t>(up[1]),
                                                static_cast<std::size_t>(up[2]))];
                }  // open vacuum: no inflow
              } else {
                value = intensity[patch.index(static_cast<std::size_t>(up[0]), static_cast<std::size_t>(up[1]),
                                              static_cast<std::size_t>(up[2]))];
              }
              inflow += w[a] * value;
            }
            const std::size_t c = patch.index(static_cast<std::size_t>(i), static_cast<std::size_t>(j),
                                              static_cast<std::size_t>(k));
            const double denominator = chi_a[c] + chi_s[c] + w[0] + w[1] + w[2];
            const double updated = (source[c] + inflow) / denominator;
            max_change = std::max(max_change, std::abs(updated - intensity[c]));
            max_value = std::max(max_value, std::abs(updated));
            intensity[c] = updated;
          }
        }
      }
    }
    for (std::size_t c = 0; c < nc; ++c) {
      double sum = 0.0;
      for (std::size_t d = 0; d < nd; ++d) sum += out.intensity[d * nc + c];
      out.mean_intensity[c] = weight * sum / (4.0 * k_pi);
    }
    out.sweeps = sweep + 1U;
    out.residual = max_value > 0.0 ? max_change / max_value : 0.0;
    if (out.residual <= settings.relative_tolerance) {
      out.converged = true;
      break;
    }
  }
  return out;
}

std::vector<std::array<double, 9>> eddingtonTensors(const UniformPatch& patch, const OrdinatesSolution& solution) {
  const std::size_t nc = patch.cells();
  const std::size_t nd = solution.directions.size();
  std::vector<std::array<double, 9>> f(nc);
  for (std::size_t c = 0; c < nc; ++c) {
    std::array<double, 9> sum{};
    double norm = 0.0;
    for (std::size_t d = 0; d < nd; ++d) {
      const double i = solution.intensity[d * nc + c];
      const auto& n = solution.directions[d];
      norm += i;
      for (std::size_t a = 0; a < 3U; ++a) {
        for (std::size_t b = 0; b < 3U; ++b) sum[a * 3U + b] += i * n[a] * n[b];
      }
    }
    if (norm > 0.0) {
      for (double& v : sum) v /= norm;
      f[c] = sum;
    } else {
      f[c] = {1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0};
    }
  }
  return f;
}

VetTransport::VetTransport(UniformPatch patch, VetSettings settings) : m_patch(patch), m_settings(settings) {
  if (!(m_settings.c > 0.0) || !(m_settings.cfl > 0.0) || !(m_settings.cfl <= 1.0) || m_settings.closure_cadence == 0U) {
    throw std::invalid_argument("invalid VET settings");
  }
  m_steps_since_closure = m_settings.closure_cadence;  // force a closure on the first step
}

double VetTransport::maxStableTimestep() const noexcept {
  double rate = 0.0;
  for (std::size_t a = 0; a < 3U; ++a) {
    if (m_patch.n[a] > 1U) rate += m_settings.c / m_patch.dx[a];
  }
  return rate > 0.0 ? m_settings.cfl / rate : 0.0;
}

void VetTransport::refreshClosure(VetState& state, const std::vector<double>& chi_a, const std::vector<double>& chi_s,
                                  const std::vector<double>& temperature) {
  const std::size_t nc = m_patch.cells();
  std::vector<double> eta(nc);
  for (std::size_t c = 0; c < nc; ++c) {
    // LTE emissivity per steradian: chi_a B, B = c a T^4 / (4 pi).
    eta[c] = chi_a[c] * m_settings.c * m_settings.radiation_constant * std::pow(temperature[c], 4) / (4.0 * k_pi);
  }
  const OrdinatesSolution solution = solveDiscreteOrdinates(m_patch, eta, chi_a, chi_s, m_settings.ordinates);
  if (!solution.converged) throw std::runtime_error("VET closure: formal solution did not converge");
  state.eddington = eddingtonTensors(m_patch, solution);
}

VetStepReport VetTransport::advance(VetState& state, const std::vector<double>& chi_a, const std::vector<double>& chi_s,
                                    const std::vector<double>& temperature, double dt) {
  const std::size_t nc = m_patch.cells();
  requireField(m_patch, chi_a, "chi_a");
  requireField(m_patch, chi_s, "chi_s");
  requireField(m_patch, temperature, "temperature");
  if (state.energy.size() != nc || state.flux.size() != nc) throw std::invalid_argument("VET state size mismatch");
  if (!(dt > 0.0)) throw std::invalid_argument("VET step needs dt > 0");
  VetStepReport report;
  const double dt_max = maxStableTimestep();
  const auto substeps = static_cast<std::size_t>(std::ceil(dt / dt_max));
  const double h = dt / static_cast<double>(std::max<std::size_t>(1U, substeps));
  const double c = m_settings.c;
  const double volume = m_patch.dx[0] * m_patch.dx[1] * m_patch.dx[2];
  std::vector<std::array<double, 4>> u(nc);
  std::vector<std::array<double, 4>> du(nc);
  for (std::size_t s = 0; s < std::max<std::size_t>(1U, substeps); ++s) {
    if (m_steps_since_closure >= m_settings.closure_cadence || state.eddington.size() != nc) {
      refreshClosure(state, chi_a, chi_s, temperature);
      ++report.closure_updates;
      m_steps_since_closure = 0;
    }
    for (std::size_t i = 0; i < nc; ++i) {
      u[i] = {state.energy[i], state.flux[i][0], state.flux[i][1], state.flux[i][2]};
      du[i] = {0.0, 0.0, 0.0, 0.0};
    }
    const auto physicalFlux = [&](const std::array<double, 4>& q, const std::array<double, 9>& f, std::size_t a) {
      return std::array<double, 4>{q[1 + a], c * c * f[a * 3U + 0U] * q[0], c * c * f[a * 3U + 1U] * q[0],
                                   c * c * f[a * 3U + 2U] * q[0]};
    };
    for (std::size_t a = 0; a < 3U; ++a) {
      if (m_patch.n[a] == 1U && m_patch.boundary == ReferenceBoundary::kPeriodic) continue;  // degenerate axis
      for (std::size_t k = 0; k < m_patch.n[2]; ++k) {
        for (std::size_t j = 0; j < m_patch.n[1]; ++j) {
          for (std::size_t i = 0; i < m_patch.n[0]; ++i) {
            // Face on the + side of cell (i, j, k) along axis a.
            std::array<long, 3> right{static_cast<long>(i), static_cast<long>(j), static_cast<long>(k)};
            right[a] += 1;
            const std::size_t lc = m_patch.index(i, j, k);
            const bool at_edge = right[a] >= static_cast<long>(m_patch.n[a]);
            std::array<double, 4> ur{};
            std::array<double, 9> fr = state.eddington[lc];
            std::size_t rc = lc;
            bool has_right_cell = true;
            if (at_edge) {
              if (m_patch.boundary == ReferenceBoundary::kPeriodic) {
                right[a] = 0;
              } else {
                has_right_cell = false;
                if (m_patch.boundary == ReferenceBoundary::kReflecting) {
                  ur = u[lc];
                  ur[1 + a] = -ur[1 + a];
                }  // vacuum ghost: ur = 0
              }
            }
            if (has_right_cell) {
              rc = m_patch.index(static_cast<std::size_t>(right[0]), static_cast<std::size_t>(right[1]),
                                 static_cast<std::size_t>(right[2]));
              ur = u[rc];
              fr = state.eddington[rc];
            }
            const auto gl = physicalFlux(u[lc], state.eddington[lc], a);
            const auto gr = physicalFlux(ur, fr, a);
            std::array<double, 4> face{};
            for (std::size_t q = 0; q < 4U; ++q) face[q] = 0.5 * (gl[q] + gr[q]) - 0.5 * c * (ur[q] - u[lc][q]);
            for (std::size_t q = 0; q < 4U; ++q) du[lc][q] -= h / m_patch.dx[a] * face[q];
            if (has_right_cell) {
              for (std::size_t q = 0; q < 4U; ++q) du[rc][q] += h / m_patch.dx[a] * face[q];
            } else {
              report.boundary_outflow_energy += face[0] * h * volume / m_patch.dx[a];
            }
            // Low-side boundary face of the first cell (non-periodic).
            if (m_patch.boundary != ReferenceBoundary::kPeriodic) {
              std::array<long, 3> ijk{static_cast<long>(i), static_cast<long>(j), static_cast<long>(k)};
              if (ijk[a] == 0) {
                std::array<double, 4> ul{};
                if (m_patch.boundary == ReferenceBoundary::kReflecting) {
                  ul = u[lc];
                  ul[1 + a] = -ul[1 + a];
                }
                const auto gll = physicalFlux(ul, state.eddington[lc], a);
                const auto grr = physicalFlux(u[lc], state.eddington[lc], a);
                std::array<double, 4> low{};
                for (std::size_t q = 0; q < 4U; ++q) low[q] = 0.5 * (gll[q] + grr[q]) - 0.5 * c * (u[lc][q] - ul[q]);
                for (std::size_t q = 0; q < 4U; ++q) du[lc][q] += h / m_patch.dx[a] * low[q];
                report.boundary_outflow_energy -= low[0] * h * volume / m_patch.dx[a];
              }
            }
          }
        }
      }
    }
    for (std::size_t i = 0; i < nc; ++i) {
      double e = u[i][0] + du[i][0];
      std::array<double, 3> f{u[i][1] + du[i][1], u[i][2] + du[i][2], u[i][3] + du[i][3]};
      // Stiff local sources, backward Euler.
      const double emission = m_settings.radiation_constant * std::pow(temperature[i], 4);
      const double ka = h * c * chi_a[i];
      const double e_new = (e + ka * emission) / (1.0 + ka);
      report.absorbed_energy += ka * e_new * volume;
      report.emitted_energy += ka * emission * volume;
      e = std::max(0.0, e_new);
      const double damp = 1.0 / (1.0 + h * c * (chi_a[i] + chi_s[i]));
      for (double& v : f) v *= damp;
      const double fmag = std::sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
      if (fmag > c * e) {
        const double scale = fmag > 0.0 ? c * e / fmag : 0.0;
        for (double& v : f) v *= scale;
        ++report.realizability_projections;
      }
      state.energy[i] = e;
      state.flux[i] = f;
    }
    ++m_steps_since_closure;
    ++report.substeps;
  }
  return report;
}

double levermorePomraningLimiter(double r) noexcept { return (2.0 + r) / (6.0 + 3.0 * r + r * r); }

FldStepReport advanceFldImplicit(const UniformPatch& patch, std::vector<double>& energy, const std::vector<double>& chi_a,
                                 const std::vector<double>& chi_s, const std::vector<double>& temperature,
                                 const std::vector<double>& source_rate, double dt, const FldSettings& settings) {
  const std::size_t nc = patch.cells();
  requireField(patch, energy, "energy");
  requireField(patch, chi_a, "chi_a");
  requireField(patch, chi_s, "chi_s");
  requireField(patch, temperature, "temperature");
  if (source_rate.size() != nc) throw std::invalid_argument("source_rate has the wrong cell count");
  if (!(dt > 0.0)) throw std::invalid_argument("FLD step needs dt > 0");
  const double c = settings.c;
  // Lagged diffusion coefficient per cell from the central-gradient ratio R.
  std::vector<double> diffusion(nc);
  const auto neighbor = [&](std::size_t i, std::size_t j, std::size_t k, std::size_t a, long offset,
                            bool& exists) -> std::size_t {
    std::array<long, 3> ijk{static_cast<long>(i), static_cast<long>(j), static_cast<long>(k)};
    ijk[a] += offset;
    exists = true;
    if (ijk[a] < 0 || ijk[a] >= static_cast<long>(patch.n[a])) {
      if (patch.boundary != ReferenceBoundary::kPeriodic) {
        exists = false;
        return 0U;
      }
      ijk[a] = static_cast<long>(wrapIndex(ijk[a], patch.n[a]));
    }
    return patch.index(static_cast<std::size_t>(ijk[0]), static_cast<std::size_t>(ijk[1]),
                       static_cast<std::size_t>(ijk[2]));
  };
  for (std::size_t k = 0; k < patch.n[2]; ++k) {
    for (std::size_t j = 0; j < patch.n[1]; ++j) {
      for (std::size_t i = 0; i < patch.n[0]; ++i) {
        const std::size_t cc = patch.index(i, j, k);
        const double chi_t = chi_a[cc] + chi_s[cc];
        if (!(chi_t > 0.0)) throw std::invalid_argument("FLD needs chi_a + chi_s > 0 in every cell");
        double grad2 = 0.0;
        for (std::size_t a = 0; a < 3U; ++a) {
          bool ep = false;
          bool em = false;
          const std::size_t p = neighbor(i, j, k, a, 1, ep);
          const std::size_t m = neighbor(i, j, k, a, -1, em);
          const double ep_value = ep ? energy[p] : (patch.boundary == ReferenceBoundary::kReflecting ? energy[cc] : 0.0);
          const double em_value = em ? energy[m] : (patch.boundary == ReferenceBoundary::kReflecting ? energy[cc] : 0.0);
          const double g = (ep_value - em_value) / (2.0 * patch.dx[a]);
          grad2 += g * g;
        }
        const double r = energy[cc] > 0.0 ? std::sqrt(grad2) / (chi_t * energy[cc]) : 0.0;
        diffusion[cc] = c * levermorePomraningLimiter(r) / chi_t;
      }
    }
  }
  // Matrix-free SPD operator.
  std::vector<double> diag(nc, 1.0);
  for (std::size_t cc = 0; cc < nc; ++cc) diag[cc] += dt * c * chi_a[cc];
  const auto faceD = [&](std::size_t a_cell, std::size_t b_cell) {
    const double da = diffusion[a_cell];
    const double db = diffusion[b_cell];
    return (da + db) > 0.0 ? 2.0 * da * db / (da + db) : 0.0;
  };
  const auto apply = [&](const std::vector<double>& x, std::vector<double>& y) {
    for (std::size_t k = 0; k < patch.n[2]; ++k) {
      for (std::size_t j = 0; j < patch.n[1]; ++j) {
        for (std::size_t i = 0; i < patch.n[0]; ++i) {
          const std::size_t cc = patch.index(i, j, k);
          double value = (1.0 + dt * c * chi_a[cc]) * x[cc];
          for (std::size_t a = 0; a < 3U; ++a) {
            if (patch.n[a] == 1U && patch.boundary == ReferenceBoundary::kPeriodic) continue;
            const double inv_dx2 = 1.0 / (patch.dx[a] * patch.dx[a]);
            for (const long off : {-1L, 1L}) {
              bool exists = false;
              const std::size_t nb = neighbor(i, j, k, a, off, exists);
              if (exists) {
                value += dt * faceD(cc, nb) * inv_dx2 * (x[cc] - x[nb]);
              } else if (patch.boundary == ReferenceBoundary::kOpenVacuum) {
                value += dt * diffusion[cc] * inv_dx2 * x[cc];  // ghost E = 0
              }
            }
          }
          y[cc] = value;
        }
      }
    }
  };
  for (std::size_t k = 0; k < patch.n[2]; ++k) {
    for (std::size_t j = 0; j < patch.n[1]; ++j) {
      for (std::size_t i = 0; i < patch.n[0]; ++i) {
        const std::size_t cc = patch.index(i, j, k);
        for (std::size_t a = 0; a < 3U; ++a) {
          if (patch.n[a] == 1U && patch.boundary == ReferenceBoundary::kPeriodic) continue;
          const double inv_dx2 = 1.0 / (patch.dx[a] * patch.dx[a]);
          for (const long off : {-1L, 1L}) {
            bool exists = false;
            const std::size_t nb = neighbor(i, j, k, a, off, exists);
            if (exists) {
              diag[cc] += dt * faceD(cc, nb) * inv_dx2;
            } else if (patch.boundary == ReferenceBoundary::kOpenVacuum) {
              diag[cc] += dt * diffusion[cc] * inv_dx2;
            }
          }
        }
      }
    }
  }
  FldStepReport report;
  const double volume = patch.dx[0] * patch.dx[1] * patch.dx[2];
  std::vector<double> b(nc);
  for (std::size_t cc = 0; cc < nc; ++cc) {
    const double emission = settings.radiation_constant * std::pow(temperature[cc], 4);
    b[cc] = energy[cc] + dt * (c * chi_a[cc] * emission + source_rate[cc]);
    report.source_energy += dt * (c * chi_a[cc] * emission + source_rate[cc]) * volume;
  }
  // Jacobi-preconditioned conjugate gradient from the previous energy.
  std::vector<double> x = energy;
  std::vector<double> r(nc);
  std::vector<double> z(nc);
  std::vector<double> p(nc);
  std::vector<double> ap(nc);
  apply(x, ap);
  double b_norm = 0.0;
  for (std::size_t cc = 0; cc < nc; ++cc) {
    r[cc] = b[cc] - ap[cc];
    z[cc] = r[cc] / diag[cc];
    p[cc] = z[cc];
    b_norm += b[cc] * b[cc];
  }
  b_norm = std::sqrt(b_norm);
  double rz = 0.0;
  for (std::size_t cc = 0; cc < nc; ++cc) rz += r[cc] * z[cc];
  double r_norm = 0.0;
  for (std::size_t it = 0; it < settings.max_cg_iterations; ++it) {
    r_norm = 0.0;
    for (const double v : r) r_norm += v * v;
    r_norm = std::sqrt(r_norm);
    report.cg_iterations = it;
    if (r_norm <= settings.cg_relative_tolerance * std::max(b_norm, 1.0e-300)) break;
    apply(p, ap);
    double pap = 0.0;
    for (std::size_t cc = 0; cc < nc; ++cc) pap += p[cc] * ap[cc];
    if (!(pap > 0.0)) throw std::runtime_error("FLD CG: operator is not positive definite");
    const double alpha = rz / pap;
    double rz_new = 0.0;
    for (std::size_t cc = 0; cc < nc; ++cc) {
      x[cc] += alpha * p[cc];
      r[cc] -= alpha * ap[cc];
      z[cc] = r[cc] / diag[cc];
      rz_new += r[cc] * z[cc];
    }
    const double beta = rz_new / rz;
    rz = rz_new;
    for (std::size_t cc = 0; cc < nc; ++cc) p[cc] = z[cc] + beta * p[cc];
  }
  report.cg_residual = b_norm > 0.0 ? r_norm / b_norm : 0.0;
  if (report.cg_residual > settings.cg_relative_tolerance) {
    throw std::runtime_error("FLD CG did not converge within max_cg_iterations (step refused)");
  }
  for (std::size_t cc = 0; cc < nc; ++cc) {
    energy[cc] = x[cc];
    report.absorbed_energy += dt * c * chi_a[cc] * x[cc] * volume;
  }
  return report;
}

}  // namespace cosmosim::radiation
