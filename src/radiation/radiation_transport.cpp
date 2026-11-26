#include "cosmosim/radiation/radiation_transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <map>
#include <utility>

#include "cosmosim/radiation/radiation_closure.hpp"

namespace cosmosim::radiation {
namespace {

using Vec3 = std::array<double, 3>;

[[nodiscard]] std::uint32_t dim(const core::RadiationPatchRecord& r, int a) noexcept {
  return a == 0 ? r.nx : (a == 1 ? r.ny : r.nz);
}
[[nodiscard]] double width(const core::RadiationPatchRecord& r, int a) noexcept {
  return a == 0 ? r.cell_width_x_comoving : (a == 1 ? r.cell_width_y_comoving : r.cell_width_z_comoving);
}
[[nodiscard]] double origin(const core::RadiationPatchRecord& r, int a) noexcept {
  return a == 0 ? r.origin_x_comoving : (a == 1 ? r.origin_y_comoving : r.origin_z_comoving);
}
[[nodiscard]] std::uint32_t cellIndex(const core::RadiationPatchRecord& r, std::uint32_t i, std::uint32_t j,
                                      std::uint32_t k) noexcept {
  return i + r.nx * (j + r.ny * k);
}
[[nodiscard]] std::uint32_t ncell(const core::RadiationPatchRecord& r) noexcept { return r.nx * r.ny * r.nz; }

// Transverse axes of a slab normal to `a`.
void transverse(int a, int& u, int& v) noexcept {
  u = (a + 1) % 3;
  v = (a + 2) % 3;
  if (u > v) std::swap(u, v);
}

[[nodiscard]] double mcSlope(double left, double centre, double right) noexcept {
  const double a = centre - left;
  const double b = right - centre;
  if (a * b <= 0.0) return 0.0;
  const double s = a > 0.0 ? 1.0 : -1.0;
  return s * std::min({2.0 * std::abs(a), 2.0 * std::abs(b), 0.5 * std::abs(a + b)});
}

[[nodiscard]] std::uint64_t faceKey(std::uint32_t record, std::uint32_t cell, int axis, int side) noexcept {
  return (static_cast<std::uint64_t>(record) << 34U) | (static_cast<std::uint64_t>(cell) << 3U) |
      (static_cast<std::uint64_t>(axis) << 1U) | static_cast<std::uint64_t>(side);
}

}  // namespace

RadiationTopology::RadiationTopology(const core::RadiationMomentSoa& moments, const RadiationTransportSettings& s)
    : m_settings(s) {
  if (!moments.isConsistent()) throw std::invalid_argument("RadiationTopology: inconsistent radiation moments");
  m_records.reserve(moments.patches.size());
  for (const core::RadiationPatchRecord& rec : moments.patches) {
    RecordInfo info;
    info.record = rec;
    m_records.push_back(info);
  }
  const bool periodic = s.boundary == core::RtBoundaryConfig::kPeriodic;
  std::uint64_t cursor = 0U;
  for (RecordInfo& info : m_records) {
    for (int a = 0; a < 3; ++a) {
      int u = 0;
      int v = 0;
      transverse(a, u, v);
      info.slab_offset[a] = cursor;
      cursor += 4ULL * dim(info.record, u) * dim(info.record, v);  // 2 sides x 2 layers
    }
  }
  m_ghosts.resize(cursor);
  for (std::uint32_t ri = 0; ri < m_records.size(); ++ri) {
    const core::RadiationPatchRecord& r = m_records[ri].record;
    for (int a = 0; a < 3; ++a) {
      int ua = 0;
      int va = 0;
      transverse(a, ua, va);
      for (int side = 0; side < 2; ++side) {
        for (int layer = 0; layer < 2; ++layer) {
          for (std::uint32_t iu = 0; iu < dim(r, ua); ++iu) {
            for (std::uint32_t iv = 0; iv < dim(r, va); ++iv) {
              const std::int64_t ia = side == 0 ? -1 - layer : static_cast<std::int64_t>(dim(r, a)) + layer;
              Vec3 p{};
              p[a] = origin(r, a) + (static_cast<double>(ia) + 0.5) * width(r, a);
              p[ua] = origin(r, ua) + (iu + 0.5) * width(r, ua);
              p[va] = origin(r, va) + (iv + 0.5) * width(r, va);
              GhostSource g;
              bool outside = false;
              for (int d = 0; d < 3; ++d) {
                const double lo = s.domain_origin_comoving[d];
                const double hi = lo + s.domain_extent_comoving[d];
                if (s.domain_extent_comoving[d] > 0.0 && (p[d] < lo || p[d] > hi)) {
                  if (periodic) {
                    p[d] -= s.domain_extent_comoving[d] * std::floor((p[d] - lo) / s.domain_extent_comoving[d]);
                  } else {
                    outside = true;
                  }
                }
              }
              auto interior_cell = [&](std::int64_t along) {
                std::uint32_t idx[3] = {0U, 0U, 0U};
                idx[a] = static_cast<std::uint32_t>(along);
                idx[ua] = iu;
                idx[va] = iv;
                return cellIndex(r, idx[0], idx[1], idx[2]);
              };
              if (outside) {
                if (s.boundary == core::RtBoundaryConfig::kReflecting) {
                  g.kind = GhostSource::Kind::kReflecting;
                  g.reflect_axis = static_cast<std::int8_t>(a);
                  g.record = ri;
                  g.cells[0] = interior_cell(side == 0 ? layer : static_cast<std::int64_t>(dim(r, a)) - 1 - layer);
                } else {
                  g.kind = GhostSource::Kind::kOutflow;
                  g.record = ri;
                  g.cells[0] = interior_cell(side == 0 ? 0 : static_cast<std::int64_t>(dim(r, a)) - 1);
                }
              } else {
                const std::int64_t q = findContainingRecord(p, r.level);
                if (q < 0) throw std::runtime_error("RadiationTopology: leaf tiling has a gap at a ghost cell");
                const core::RadiationPatchRecord& rq = m_records[static_cast<std::size_t>(q)].record;
                auto cell_of = [&](const core::RadiationPatchRecord& rr, const Vec3& pt) {
                  std::uint32_t idx[3];
                  for (int d = 0; d < 3; ++d) {
                    const double t = (pt[d] - origin(rr, d)) / width(rr, d);
                    idx[d] = static_cast<std::uint32_t>(std::clamp<std::int64_t>(
                        static_cast<std::int64_t>(std::floor(t)), 0, static_cast<std::int64_t>(dim(rr, d)) - 1));
                  }
                  return cellIndex(rr, idx[0], idx[1], idx[2]);
                };
                if (rq.level == r.level) {
                  g.kind = GhostSource::Kind::kSameLevel;
                  g.record = static_cast<std::uint32_t>(q);
                  g.cells[0] = cell_of(rq, p);
                } else if (rq.level < r.level) {
                  if (r.level - rq.level != 1) throw std::runtime_error("RadiationTopology: tiling is not 2:1 balanced");
                  g.kind = GhostSource::Kind::kCoarser;
                  g.record = static_cast<std::uint32_t>(q);
                  g.cells[0] = cell_of(rq, p);
                } else {
                  if (rq.level - r.level != 1) throw std::runtime_error("RadiationTopology: tiling is not 2:1 balanced");
                  // Eight ratio-2 children of the ghost cell; all must lie in
                  // ONE finer record for this pass-01 representation.
                  g.kind = GhostSource::Kind::kFiner;
                  g.record = static_cast<std::uint32_t>(q);
                  g.cell_count = 8U;
                  std::uint32_t n = 0U;
                  for (int cx = -1; cx <= 1; cx += 2) {
                    for (int cy = -1; cy <= 1; cy += 2) {
                      for (int cz = -1; cz <= 1; cz += 2) {
                        const Vec3 pc{p[0] + 0.25 * cx * width(r, 0), p[1] + 0.25 * cy * width(r, 1),
                                      p[2] + 0.25 * cz * width(r, 2)};
                        if (findContainingRecord(pc, rq.level) != q) {
                          throw std::runtime_error(
                              "RadiationTopology: ratio-2 children of a coarse ghost span several fine records");
                        }
                        g.cells[n++] = cell_of(rq, pc);
                      }
                    }
                  }
                  if (layer == 0) {
                    std::uint32_t idx[3];
                    idx[a] = side == 0 ? 0U : dim(r, a) - 1U;
                    idx[ua] = iu;
                    idx[va] = iv;
                    m_coarse_fine_faces.push_back(faceKey(ri, cellIndex(r, idx[0], idx[1], idx[2]), a, side));
                  }
                }
              }
              const std::uint64_t slab = m_records[ri].slab_offset[a] +
                  (static_cast<std::uint64_t>(side * 2 + layer) * dim(r, ua) + iu) * dim(r, va) + iv;
              m_ghosts[slab] = g;
            }
          }
        }
      }
    }
  }
  std::sort(m_coarse_fine_faces.begin(), m_coarse_fine_faces.end());
}

std::int64_t RadiationTopology::findContainingRecord(const std::array<double, 3>& p, int preferred_level) const {
  std::int64_t best = -1;
  for (std::size_t i = 0; i < m_records.size(); ++i) {
    const core::RadiationPatchRecord& r = m_records[i].record;
    bool inside = true;
    for (int d = 0; d < 3 && inside; ++d) {
      const double lo = origin(r, d);
      const double hi = lo + width(r, d) * dim(r, d);
      inside = p[d] > lo && p[d] < hi;
    }
    if (inside) {
      if (r.level == preferred_level) return static_cast<std::int64_t>(i);
      best = static_cast<std::int64_t>(i);
    }
  }
  return best;
}

const RadiationTopology::GhostSource& RadiationTopology::ghost(std::uint32_t record, int axis, int side, int layer,
                                                               std::uint32_t u, std::uint32_t v) const {
  const core::RadiationPatchRecord& r = m_records[record].record;
  int ua = 0;
  int va = 0;
  transverse(axis, ua, va);
  const std::uint64_t slab = m_records[record].slab_offset[axis] +
      (static_cast<std::uint64_t>(side * 2 + layer) * dim(r, ua) + u) * dim(r, va) + v;
  return m_ghosts[slab];
}

bool RadiationTopology::faceIsCoarseFine(std::uint32_t record, std::uint32_t cell, int axis, int side) const {
  return std::binary_search(m_coarse_fine_faces.begin(), m_coarse_fine_faces.end(), faceKey(record, cell, axis, side));
}

RadiationLevelSolver::RadiationLevelSolver(const core::RadiationMomentSoa& moments, const RadiationTopology& topology,
                                           RadiationTransportSettings settings)
    : m_moments(moments), m_topology(topology), m_settings(settings), m_groups(moments.group_count) {
  if (!(m_settings.c_hat_code > 0.0)) throw std::invalid_argument("RadiationLevelSolver: c_hat must be > 0");
}

double RadiationLevelSolver::stableTimeStep(double cfl, double a) const {
  double best = std::numeric_limits<double>::infinity();
  for (const core::RadiationPatchRecord& r : m_moments.patches) {
    const double s = m_settings.c_hat_code *
        (1.0 / (a * r.cell_width_x_comoving) + 1.0 / (a * r.cell_width_y_comoving) + 1.0 / (a * r.cell_width_z_comoving));
    best = std::min(best, cfl / s);
  }
  return best;
}

bool RadiationLevelSolver::admissible(const RadiationWorkState& st, std::uint64_t& bad) const {
  bad = 0U;
  for (std::size_t i = 0; i < st.size(); ++i) {
    const M1State m{st.e[i], {st.gx[i], st.gy[i], st.gz[i]}};
    if (!m1Admissible(m, m_settings.c_hat_code, 1.0e-12)) ++bad;
  }
  return bad == 0U;
}

bool RadiationLevelSolver::transportRate(const RadiationWorkState& st, double inv_a, bool first_order,
                                         RadiationWorkState& rate, RadiationStageDiagnostics& diag) const {
  const double c = m_settings.c_hat_code;
  const auto& recs = m_moments.patches;
  rate.resize(st.size());
  using Kind = RadiationTopology::GhostSource::Kind;
  // Value of (record, group, cell) or of a ghost slab entry.
  auto load = [&](std::uint32_t ri, std::uint32_t g, std::uint32_t cell) {
    const core::RadiationPatchRecord& r = recs[ri];
    const std::size_t k = r.value_offset + static_cast<std::size_t>(g) * ncell(r) + cell;
    return M1State{st.e[k], {st.gx[k], st.gy[k], st.gz[k]}};
  };
  auto ghostValue = [&](const RadiationTopology::GhostSource& gs, std::uint32_t g) {
    if (gs.kind == Kind::kFiner) {
      M1State sum{};
      for (std::uint32_t n = 0; n < gs.cell_count; ++n) {
        const M1State m = load(gs.record, g, gs.cells[n]);
        sum.energy += m.energy / gs.cell_count;
        for (int d = 0; d < 3; ++d) sum.flux[d] += m.flux[d] / gs.cell_count;
      }
      return sum;
    }
    M1State m = load(gs.record, g, gs.cells[0]);
    if (gs.kind == Kind::kReflecting) m.flux[gs.reflect_axis] = -m.flux[gs.reflect_axis];
    return m;
  };
  // Cell along axis a at signed index ia (others fixed), in or out of record.
  auto valueAlong = [&](std::uint32_t ri, std::uint32_t g, int a, std::int64_t ia, std::uint32_t iu, std::uint32_t iv,
                        int ua, int va) {
    const core::RadiationPatchRecord& r = recs[ri];
    const std::int64_t n = dim(r, a);
    if (ia >= 0 && ia < n) {
      std::uint32_t idx[3];
      idx[a] = static_cast<std::uint32_t>(ia);
      idx[ua] = iu;
      idx[va] = iv;
      return load(ri, g, cellIndex(r, idx[0], idx[1], idx[2]));
    }
    const int side = ia < 0 ? 0 : 1;
    const int layer = ia < 0 ? static_cast<int>(-1 - ia) : static_cast<int>(ia - n);
    return ghostValue(m_topology.ghost(ri, a, side, std::min(layer, 1), iu, iv), g);
  };
  auto reconstruct = [&](const M1State& l, const M1State& m, const M1State& r, M1State& minus, M1State& plus) {
    if (first_order) {
      minus = m;
      plus = m;
      return;
    }
    M1State half{};
    half.energy = 0.5 * mcSlope(l.energy, m.energy, r.energy);
    for (int d = 0; d < 3; ++d) half.flux[d] = 0.5 * mcSlope(l.flux[d], m.flux[d], r.flux[d]);
    M1State neg{-half.energy, {-half.flux[0], -half.flux[1], -half.flux[2]}};
    const double theta = std::min(m1AdmissibleFraction(m, half, c), m1AdmissibleFraction(m, neg, c));
    plus.energy = m.energy + theta * half.energy;
    minus.energy = m.energy - theta * half.energy;
    for (int d = 0; d < 3; ++d) {
      plus.flux[d] = m.flux[d] + theta * half.flux[d];
      minus.flux[d] = m.flux[d] - theta * half.flux[d];
    }
  };
  bool finite = true;
  auto riemann = [&](const M1State& left, const M1State& right, int a) {
    const auto fl = m1PhysicalFlux(left, a, c);
    const auto fr = m1PhysicalFlux(right, a, c);
    const std::array<double, 4> ul{left.energy, left.flux[0], left.flux[1], left.flux[2]};
    const std::array<double, 4> ur{right.energy, right.flux[0], right.flux[1], right.flux[2]};
    std::array<double, 4> f{};
    const bool use_hll = !first_order && m_settings.riemann == core::RtRiemannSolverConfig::kHll;
    if (use_hll) {
      const M1WaveSpeeds sl = m1WaveSpeeds(left, a, c);
      const M1WaveSpeeds sr = m1WaveSpeeds(right, a, c);
      const double s_l = std::min(sl.min, sr.min);
      const double s_r = std::max(sl.max, sr.max);
      if (s_l >= 0.0) return fl;
      if (s_r <= 0.0) return fr;
      bool ok = s_r > s_l;
      for (int q = 0; q < 4 && ok; ++q) {
        f[q] = (s_r * fl[q] - s_l * fr[q] + s_l * s_r * (ur[q] - ul[q])) / (s_r - s_l);
        ok = std::isfinite(f[q]);
      }
      if (ok) return f;
      ++diag.hll_fallback_faces;
    }
    for (int q = 0; q < 4; ++q) {
      f[q] = 0.5 * (fl[q] + fr[q]) - 0.5 * c * (ur[q] - ul[q]);
      if (!std::isfinite(f[q])) finite = false;
    }
    return f;
  };
  // Face flux between signed cells (ia-1, ia) along a.
  auto faceFlux = [&](std::uint32_t ri, std::uint32_t g, int a, std::int64_t ia, std::uint32_t iu, std::uint32_t iv,
                      int ua, int va) {
    const M1State lm2 = valueAlong(ri, g, a, ia - 2, iu, iv, ua, va);
    const M1State lm1 = valueAlong(ri, g, a, ia - 1, iu, iv, ua, va);
    const M1State l0 = valueAlong(ri, g, a, ia, iu, iv, ua, va);
    const M1State lp1 = valueAlong(ri, g, a, ia + 1, iu, iv, ua, va);
    M1State left_minus, left_plus, right_minus, right_plus;
    reconstruct(lm2, lm1, l0, left_minus, left_plus);
    reconstruct(lm1, l0, lp1, right_minus, right_plus);
    return riemann(left_plus, right_minus, a);
  };

  // Pre-pass: fine fluxes on coarse/fine faces, accumulated into the coarse
  // register (area ratio 1/4 for ratio-2 refinement).
  std::map<std::pair<std::uint64_t, std::uint32_t>, std::array<double, 4>> reg;
  for (std::uint32_t ri = 0; ri < recs.size(); ++ri) {
    const core::RadiationPatchRecord& r = recs[ri];
    for (int a = 0; a < 3; ++a) {
      int ua = 0;
      int va = 0;
      transverse(a, ua, va);
      for (int side = 0; side < 2; ++side) {
        for (std::uint32_t iu = 0; iu < dim(r, ua); ++iu) {
          for (std::uint32_t iv = 0; iv < dim(r, va); ++iv) {
            const auto& gs = m_topology.ghost(ri, a, side, 0, iu, iv);
            if (gs.kind != Kind::kCoarser) continue;
            for (std::uint32_t g = 0; g < m_groups; ++g) {
              const std::int64_t face = side == 0 ? 0 : static_cast<std::int64_t>(dim(r, a));
              const auto f = faceFlux(ri, g, a, face, iu, iv, ua, va);
              const auto key = std::make_pair(faceKey(gs.record, gs.cells[0], a, 1 - side), g);
              auto& acc = reg[key];
              for (int q = 0; q < 4; ++q) acc[q] += 0.25 * f[q];
            }
          }
        }
      }
    }
  }
  // Main pass: divergence per record with shared/replaced faces.
  for (std::uint32_t ri = 0; ri < recs.size(); ++ri) {
    const core::RadiationPatchRecord& r = recs[ri];
    for (std::uint32_t g = 0; g < m_groups; ++g) {
      const std::size_t base = r.value_offset + static_cast<std::size_t>(g) * ncell(r);
      for (std::uint32_t c0 = 0; c0 < ncell(r); ++c0) {
        rate.e[base + c0] = 0.0;
        rate.gx[base + c0] = 0.0;
        rate.gy[base + c0] = 0.0;
        rate.gz[base + c0] = 0.0;
      }
      for (int a = 0; a < 3; ++a) {
        int ua = 0;
        int va = 0;
        transverse(a, ua, va);
        const double factor = inv_a / width(r, a);
        for (std::uint32_t iu = 0; iu < dim(r, ua); ++iu) {
          for (std::uint32_t iv = 0; iv < dim(r, va); ++iv) {
            for (std::int64_t face = 0; face <= static_cast<std::int64_t>(dim(r, a)); ++face) {
              std::array<double, 4> f{};
              const int side = face == 0 ? 0 : 1;
              bool replaced = false;
              if (face == 0 || face == static_cast<std::int64_t>(dim(r, a))) {
                std::uint32_t idx[3];
                idx[a] = face == 0 ? 0U : dim(r, a) - 1U;
                idx[ua] = iu;
                idx[va] = iv;
                const std::uint32_t cell = cellIndex(r, idx[0], idx[1], idx[2]);
                if (m_topology.faceIsCoarseFine(ri, cell, a, side)) {
                  const auto it = reg.find(std::make_pair(faceKey(ri, cell, a, side), g));
                  if (it == reg.end()) throw std::logic_error("RadiationLevelSolver: missing fine flux register");
                  f = it->second;
                  replaced = true;
                }
              }
              if (!replaced) f = faceFlux(ri, g, a, face, iu, iv, ua, va);
              if (face == 0 || face == static_cast<std::int64_t>(dim(r, a))) {
                const auto& gs = m_topology.ghost(ri, a, side, 0, iu, iv);
                if (gs.kind == Kind::kOutflow || gs.kind == Kind::kReflecting) {
                  const double area = width(r, ua) * width(r, va);
                  diag.boundary_energy_outflow_code += (side == 1 ? f[0] : -f[0]) * area;
                }
              }
              auto apply = [&](std::int64_t cell_along, double sign) {
                if (cell_along < 0 || cell_along >= static_cast<std::int64_t>(dim(r, a))) return;
                std::uint32_t idx[3];
                idx[a] = static_cast<std::uint32_t>(cell_along);
                idx[ua] = iu;
                idx[va] = iv;
                const std::size_t k = base + cellIndex(r, idx[0], idx[1], idx[2]);
                rate.e[k] += sign * factor * f[0];
                rate.gx[k] += sign * factor * f[1];
                rate.gy[k] += sign * factor * f[2];
                rate.gz[k] += sign * factor * f[3];
              };
              apply(face - 1, -1.0);  // flux leaves the left cell
              apply(face, +1.0);      // and enters the right cell
            }
          }
        }
      }
    }
  }
  return finite;
}

}  // namespace cosmosim::radiation
