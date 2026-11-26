#include "cosmosim/radiation/radiation_angular.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>

namespace cosmosim::radiation {
namespace {

using Vec3 = std::array<double, 3>;
constexpr double k_pi = 3.141592653589793238462643383279502884;

[[nodiscard]] double dot(const Vec3& a, const Vec3& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
[[nodiscard]] Vec3 normalize(const Vec3& a) noexcept {
  const double n = std::sqrt(dot(a, a));
  return n > 0.0 ? Vec3{a[0] / n, a[1] / n, a[2] / n} : Vec3{0.0, 0.0, 1.0};
}
[[nodiscard]] Vec3 cross(const Vec3& a, const Vec3& b) noexcept {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

[[nodiscard]] double recordWidth(const core::RadiationPatchRecord& r, int a) noexcept {
  return a == 0 ? r.cell_width_x_comoving : (a == 1 ? r.cell_width_y_comoving : r.cell_width_z_comoving);
}
[[nodiscard]] double recordOrigin(const core::RadiationPatchRecord& r, int a) noexcept {
  return a == 0 ? r.origin_x_comoving : (a == 1 ? r.origin_y_comoving : r.origin_z_comoving);
}
[[nodiscard]] std::uint32_t recordDim(const core::RadiationPatchRecord& r, int a) noexcept {
  return a == 0 ? r.nx : (a == 1 ? r.ny : r.nz);
}

[[nodiscard]] std::vector<std::size_t> cellBases(const core::RadiationMomentSoa& m) {
  std::vector<std::size_t> base(m.patches.size() + 1U, 0U);
  for (std::size_t r = 0; r < m.patches.size(); ++r) {
    base[r + 1U] = base[r] + static_cast<std::size_t>(m.patches[r].nx) * m.patches[r].ny * m.patches[r].nz;
  }
  return base;
}

// Jacobi eigen-decomposition of a symmetric 3x3 matrix.
void jacobiEigen(std::array<double, 9> a, std::array<double, 3>& eval, std::array<double, 9>& evec) {
  evec = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  for (int sweep = 0; sweep < 50; ++sweep) {
    const double off = a[1] * a[1] + a[2] * a[2] + a[5] * a[5];
    if (off < 1.0e-30) break;
    for (int p = 0; p < 2; ++p) {
      for (int q = p + 1; q < 3; ++q) {
        const double apq = a[3 * p + q];
        if (std::abs(apq) < 1.0e-300) continue;
        const double theta = 0.5 * (a[3 * q + q] - a[3 * p + p]) / apq;
        const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0);
        const double s = t * c;
        for (int k = 0; k < 3; ++k) {
          const double akp = a[3 * k + p];
          const double akq = a[3 * k + q];
          a[3 * k + p] = c * akp - s * akq;
          a[3 * k + q] = s * akp + c * akq;
        }
        for (int k = 0; k < 3; ++k) {
          const double apk = a[3 * p + k];
          const double aqk = a[3 * q + k];
          a[3 * p + k] = c * apk - s * aqk;
          a[3 * q + k] = s * apk + c * aqk;
        }
        for (int k = 0; k < 3; ++k) {
          const double vkp = evec[3 * k + p];
          const double vkq = evec[3 * k + q];
          evec[3 * k + p] = c * vkp - s * vkq;
          evec[3 * k + q] = s * vkp + c * vkq;
        }
      }
    }
  }
  eval = {a[0], a[4], a[8]};
}

}  // namespace

RadiationRecordLocator::RadiationRecordLocator(const core::RadiationMomentSoa& m, const physics::PeriodicBox& box,
                                               std::array<double, 3> origin)
    : m_moments(m), m_box(box), m_origin(origin) {}

bool RadiationRecordLocator::insideDomain(const Vec3& p) const {
  for (int a = 0; a < 3; ++a) {
    if (m_box.extent[a] > 0.0 && (p[a] < m_origin[a] || p[a] > m_origin[a] + m_box.extent[a])) return false;
  }
  return true;
}

bool RadiationRecordLocator::locate(const Vec3& p, RecordCellRef& out) const {
  for (std::uint32_t r = 0; r < m_moments.patches.size(); ++r) {
    const core::RadiationPatchRecord& rr = m_moments.patches[r];
    std::uint32_t idx[3];
    bool inside = true;
    for (int a = 0; a < 3 && inside; ++a) {
      const double t = (p[a] - recordOrigin(rr, a)) / recordWidth(rr, a);
      if (!(t >= 0.0) || t >= recordDim(rr, a)) {
        inside = false;
      } else {
        idx[a] = static_cast<std::uint32_t>(t);
      }
    }
    if (inside) {
      out.record = r;
      out.cell = idx[0] + rr.nx * (idx[1] + rr.ny * idx[2]);
      return true;
    }
  }
  return false;
}

double RadiationRecordLocator::exitDistance(const RecordCellRef& ref, const Vec3& p, const Vec3& d) const {
  const core::RadiationPatchRecord& rr = m_moments.patches[ref.record];
  const std::uint32_t ix = ref.cell % rr.nx;
  const std::uint32_t iy = (ref.cell / rr.nx) % rr.ny;
  const std::uint32_t iz = ref.cell / (rr.nx * rr.ny);
  const std::uint32_t idx[3] = {ix, iy, iz};
  double t = std::numeric_limits<double>::infinity();
  for (int a = 0; a < 3; ++a) {
    if (d[a] == 0.0) continue;
    const double lo = recordOrigin(rr, a) + idx[a] * recordWidth(rr, a);
    const double hi = lo + recordWidth(rr, a);
    const double ta = d[a] > 0.0 ? (hi - p[a]) / d[a] : (lo - p[a]) / d[a];
    t = std::min(t, std::max(ta, 0.0));
  }
  return t;
}

double RadiationRecordLocator::cellArea(const RecordCellRef& ref) const {
  const core::RadiationPatchRecord& rr = m_moments.patches[ref.record];
  return std::cbrt(rr.cell_width_x_comoving * rr.cell_width_y_comoving * rr.cell_width_z_comoving) *
      std::cbrt(rr.cell_width_x_comoving * rr.cell_width_y_comoving * rr.cell_width_z_comoving);
}

RayTracingResult traceDirectRays(const RadiationRecordLocator& loc, const core::RadiationMomentSoa& m, const Vec3& src,
                                 double e0, const RadiationOpacityFunction& opacity,
                                 const physics::ProperFrameContext& frame, const RayTracingSettings& set) {
  RayTracingResult out;
  const auto base = cellBases(m);
  out.absorbed_energy.assign(base.back(), 0.0);
  struct Ray {
    std::uint64_t id;
    Vec3 pos;
    Vec3 dir;
    double energy;
    double solid_angle;
    double travelled;  // stored length from the source
  };
  std::deque<Ray> queue;
  const std::uint32_t n0 = std::max<std::uint32_t>(set.initial_rays, 1U);
  const double golden = k_pi * (3.0 - std::sqrt(5.0));
  for (std::uint32_t k = 0; k < n0; ++k) {
    const double z = 1.0 - (2.0 * k + 1.0) / n0;
    const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    const double phi = golden * k;
    queue.push_back(Ray{k + 1U, src, {r * std::cos(phi), r * std::sin(phi), z}, e0 / n0, 4.0 * k_pi / n0, 0.0});
  }
  const double e_min = set.min_energy_fraction * e0;
  while (!queue.empty()) {
    Ray ray = queue.front();
    queue.pop_front();
    ++out.rays_traced;
    for (std::uint64_t guard = 0; guard < 100000000ULL; ++guard) {
      if (ray.energy <= e_min) {
        out.terminated_energy += ray.energy;
        break;
      }
      if (!loc.insideDomain(ray.pos)) {
        out.escaped_energy += ray.energy;
        break;
      }
      RecordCellRef ref;
      if (!loc.locate(ray.pos, ref)) {
        out.escaped_energy += ray.energy;
        break;
      }
      // Adaptive splitting (Abel & Wandelt): Omega r^2 > A_cell / q.
      if (ray.solid_angle * ray.travelled * ray.travelled > loc.cellArea(ref) / set.split_quality) {
        if (queue.size() + 4U > set.max_rays) {
          ++out.queue_overflows;
          out.terminated_energy += ray.energy;
          break;
        }
        const Vec3 helper = std::abs(ray.dir[2]) < 0.9 ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
        const Vec3 e1 = normalize(cross(ray.dir, helper));
        const Vec3 e2 = cross(ray.dir, e1);
        const double off = 0.5 * std::sqrt(ray.solid_angle / 4.0);
        for (int sx = -1; sx <= 1; sx += 2) {
          for (int sy = -1; sy <= 1; sy += 2) {
            const Vec3 d = normalize({ray.dir[0] + off * (sx * e1[0] + sy * e2[0]),
                                      ray.dir[1] + off * (sx * e1[1] + sy * e2[1]),
                                      ray.dir[2] + off * (sx * e1[2] + sy * e2[2])});
            const std::uint64_t child = ray.id * 4ULL + static_cast<std::uint64_t>((sx + 1) + (sy + 1) / 2);
            queue.push_back(Ray{child, ray.pos, d, 0.25 * ray.energy, 0.25 * ray.solid_angle, ray.travelled});
          }
        }
        ++out.splits;
        break;
      }
      const double s = loc.exitDistance(ref, ray.pos, ray.dir);
      const double step = s + 1.0e-9 * loc.cellArea(ref) / std::max(std::sqrt(loc.cellArea(ref)), 1.0e-300);
      const double tau = opacity(ref.record, ref.cell) * frame.properLength(s);
      const double absorbed = ray.energy * -std::expm1(-tau);
      out.absorbed_energy[base[ref.record] + ref.cell] += absorbed;
      ray.energy -= absorbed;
      for (int a = 0; a < 3; ++a) ray.pos[a] += step * ray.dir[a];
      ray.travelled += step;
    }
  }
  return out;
}

VetTensor vetFromAngularSamples(const std::vector<Vec3>& dirs, const std::vector<double>& intensity,
                                const std::vector<double>& weights) {
  VetTensor out;
  if (dirs.size() != intensity.size() || dirs.size() != weights.size()) {
    throw std::invalid_argument("vetFromAngularSamples: sample arrays differ in length");
  }
  std::array<double, 9> p{};
  double e = 0.0;
  for (std::size_t k = 0; k < dirs.size(); ++k) {
    if (intensity[k] < 0.0 || !(weights[k] > 0.0)) throw std::invalid_argument("vetFromAngularSamples: invalid sample");
    const Vec3 n = normalize(dirs[k]);
    const double wi = weights[k] * intensity[k];
    e += wi;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) p[3 * i + j] += wi * n[i] * n[j];
    }
  }
  if (!(e > 0.0)) return out;
  for (double& v : p) v /= e;
  for (int i = 0; i < 3; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      const double s = 0.5 * (p[3 * i + j] + p[3 * j + i]);
      p[3 * i + j] = s;
      p[3 * j + i] = s;
    }
  }
  std::array<double, 3> ev{};
  std::array<double, 9> vec{};
  jacobiEigen(p, ev, vec);
  double trace = 0.0;
  for (double& l : ev) {
    l = std::max(l, 0.0);
    trace += l;
  }
  if (!(trace > 0.0)) return out;
  for (double& l : ev) l /= trace;
  double diff = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double v = 0.0;
      for (int k = 0; k < 3; ++k) v += vec[3 * i + k] * ev[k] * vec[3 * j + k];
      diff += (v - p[3 * i + j]) * (v - p[3 * i + j]);
      out.f[3 * i + j] = v;
    }
  }
  out.projection_size = std::sqrt(diff);
  out.defined = true;
  return out;
}

std::uint64_t radiationCounterRandom(std::uint64_t seed, std::uint64_t stream, std::uint64_t counter) noexcept {
  std::uint64_t z = seed ^ (stream * 0x9e3779b97f4a7c15ULL) ^ (counter * 0xd1b54a32d192ed03ULL);
  z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31U);
}

MonteCarloResult runMonteCarloReference(const RadiationRecordLocator& loc, const core::RadiationMomentSoa& m,
                                        const Vec3& src, double e0, const RadiationOpacityFunction& opacity,
                                        const physics::ProperFrameContext& frame, const MonteCarloSettings& set) {
  MonteCarloResult out;
  const auto base = cellBases(m);
  out.absorbed_energy.assign(base.back(), 0.0);
  std::vector<double> sum_sq(base.back(), 0.0);
  if (set.packets == 0U) return out;
  const double eps = e0 / static_cast<double>(set.packets);
  auto uniform = [&](std::uint64_t packet, std::uint64_t& counter) {
    const std::uint64_t r = radiationCounterRandom(set.seed, packet, counter++);
    return (static_cast<double>(r >> 11U) + 0.5) * (1.0 / 9007199254740992.0);
  };
  auto isotropic = [&](std::uint64_t packet, std::uint64_t& counter) {
    const double z = 2.0 * uniform(packet, counter) - 1.0;
    const double phi = 2.0 * k_pi * uniform(packet, counter);
    const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    return Vec3{r * std::cos(phi), r * std::sin(phi), z};
  };
  std::uint64_t escaped = 0U;
  for (std::uint64_t batch_begin = 0; batch_begin < set.packets; batch_begin += std::max<std::uint64_t>(set.batch, 1U)) {
    const std::uint64_t batch_end = std::min(set.packets, batch_begin + std::max<std::uint64_t>(set.batch, 1U));
    for (std::uint64_t pk = batch_begin; pk < batch_end; ++pk) {
      std::uint64_t counter = 0U;
      Vec3 pos = src;
      Vec3 dir = isotropic(pk, counter);
      double tau_left = -std::log(uniform(pk, counter));
      bool done = false;
      for (std::uint32_t ev = 0; ev < set.max_events_per_packet && !done; ++ev) {
        if (!loc.insideDomain(pos)) {
          ++escaped;
          out.escaped_energy += eps;
          done = true;
          break;
        }
        RecordCellRef ref;
        if (!loc.locate(pos, ref)) {
          ++escaped;
          out.escaped_energy += eps;
          done = true;
          break;
        }
        const double s = loc.exitDistance(ref, pos, dir);
        const double chi = opacity(ref.record, ref.cell);
        const double tau_cell = chi * frame.properLength(s);
        if (chi > 0.0 && tau_cell >= tau_left) {
          const double s_int = tau_left / (chi * frame.lengthScale());
          for (int a = 0; a < 3; ++a) pos[a] += s_int * dir[a];
          if (uniform(pk, counter) < set.scattering_albedo) {
            dir = isotropic(pk, counter);
            tau_left = -std::log(uniform(pk, counter));
          } else {
            const std::size_t k = base[ref.record] + ref.cell;
            out.absorbed_energy[k] += eps;
            sum_sq[k] += eps * eps;
            done = true;
          }
        } else {
          tau_left -= tau_cell;
          const double step = s * (1.0 + 1.0e-9) + 1.0e-12 * std::sqrt(loc.cellArea(ref));
          for (int a = 0; a < 3; ++a) pos[a] += step * dir[a];
        }
      }
      if (!done) ++out.truncated_packets;
      ++out.packets;
    }
  }
  const double n = static_cast<double>(out.packets);
  out.absorbed_energy_variance.assign(base.back(), 0.0);
  for (std::size_t k = 0; k < base.back(); ++k) {
    const double mean = out.absorbed_energy[k] / n;
    const double mean_sq = sum_sq[k] / n;
    out.absorbed_energy_variance[k] = n * std::max(mean_sq - mean * mean, 0.0);  // variance of the total estimator
  }
  const double f = static_cast<double>(escaped) / n;
  out.escaped_fraction_standard_error = std::sqrt(std::max(f * (1.0 - f), 0.0) / n);
  return out;
}

}  // namespace cosmosim::radiation
