// CHUI-EXT-PHYS-REF-v1 S03 leaf ray service (source only, NOT compiled in pass 01).

#include "cosmosim/lightcone/ray_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cosmosim::lightcone {
namespace {

[[nodiscard]] bool contains(const LeafPatch& leaf, const std::array<double, 3>& p) noexcept {
  for (std::size_t a = 0; a < 3U; ++a) {
    if (!(p[a] >= leaf.origin[a]) || !(p[a] < leaf.origin[a] + leaf.extent[a])) return false;
  }
  return true;
}

}  // namespace

LeafGeometry::LeafGeometry(std::vector<LeafPatch> leaves, std::array<double, 3> domain_origin,
                           std::array<double, 3> domain_extent, RayBoundary boundary, std::size_t buckets_per_axis)
    : m_leaves(std::move(leaves)), m_origin(domain_origin), m_extent(domain_extent), m_boundary(boundary) {
  if (m_leaves.empty()) throw std::invalid_argument("LeafGeometry needs at least one leaf patch");
  for (std::size_t a = 0; a < 3U; ++a) {
    if (!std::isfinite(m_origin[a]) || !(m_extent[a] > 0.0)) throw std::invalid_argument("invalid ray domain");
  }
  m_finest_width = std::numeric_limits<double>::infinity();
  for (const LeafPatch& leaf : m_leaves) {
    std::size_t cells = 1;
    for (std::size_t a = 0; a < 3U; ++a) {
      if (leaf.dims[a] == 0U || !(leaf.extent[a] > 0.0)) throw std::invalid_argument("leaf patch without geometry");
      if (leaf.origin[a] < m_origin[a] - 1.0e-12 * m_extent[a] ||
          leaf.origin[a] + leaf.extent[a] > m_origin[a] + m_extent[a] * (1.0 + 1.0e-12)) {
        throw std::invalid_argument("leaf patch " + std::to_string(leaf.patch_id) + " lies outside the ray domain");
      }
      m_finest_width = std::min(m_finest_width, leaf.extent[a] / static_cast<double>(leaf.dims[a]));
      cells *= leaf.dims[a];
    }
    m_cell_count += cells;
  }
  if (buckets_per_axis == 0U) {
    const double target = std::cbrt(static_cast<double>(m_leaves.size())) * 2.0;
    buckets_per_axis = static_cast<std::size_t>(std::clamp(target, 1.0, 256.0));
  }
  m_buckets = buckets_per_axis;
  const std::size_t total = m_buckets * m_buckets * m_buckets;
  std::vector<std::uint32_t> counts(total + 1U, 0U);
  const auto range = [&](const LeafPatch& leaf, std::size_t a, std::size_t& lo, std::size_t& hi) {
    const double width = m_extent[a] / static_cast<double>(m_buckets);
    const double l = (leaf.origin[a] - m_origin[a]) / width;
    const double h = (leaf.origin[a] + leaf.extent[a] - m_origin[a]) / width;
    lo = static_cast<std::size_t>(std::clamp(std::floor(l), 0.0, static_cast<double>(m_buckets - 1U)));
    hi = static_cast<std::size_t>(std::clamp(std::ceil(h) - 1.0, 0.0, static_cast<double>(m_buckets - 1U)));
  };
  for (int pass = 0; pass < 2; ++pass) {
    if (pass == 1) {
      for (std::size_t i = 1; i <= total; ++i) counts[i] += counts[i - 1U];
      m_bucket_begin = counts;
      m_bucket_leaves.assign(counts.back(), 0U);
      std::fill(counts.begin(), counts.end(), 0U);
    }
    for (std::size_t l = 0; l < m_leaves.size(); ++l) {
      std::array<std::size_t, 3> lo{};
      std::array<std::size_t, 3> hi{};
      for (std::size_t a = 0; a < 3U; ++a) range(m_leaves[l], a, lo[a], hi[a]);
      for (std::size_t i = lo[0]; i <= hi[0]; ++i) {
        for (std::size_t j = lo[1]; j <= hi[1]; ++j) {
          for (std::size_t k = lo[2]; k <= hi[2]; ++k) {
            const std::size_t b = (i * m_buckets + j) * m_buckets + k;
            if (pass == 0) {
              ++counts[b + 1U];
            } else {
              m_bucket_leaves[m_bucket_begin[b] + counts[b]++] = static_cast<std::uint32_t>(l);
            }
          }
        }
      }
    }
  }
}

LeafGeometry LeafGeometry::fromPatches(const core::PatchSoa& patches, std::array<double, 3> domain_origin,
                                       std::array<double, 3> domain_extent, RayBoundary boundary) {
  std::vector<LeafPatch> leaves;
  for (std::size_t p = 0; p < patches.size(); ++p) {
    if (patches.cell_count[p] == 0U) continue;  // metadata-only parent: never a leaf
    LeafPatch leaf;
    leaf.patch_id = patches.patch_id[p];
    leaf.level = patches.level[p];
    leaf.origin = {patches.origin_x_comoving[p], patches.origin_y_comoving[p], patches.origin_z_comoving[p]};
    leaf.extent = {patches.extent_x_comoving[p], patches.extent_y_comoving[p], patches.extent_z_comoving[p]};
    leaf.dims = {patches.cell_dim_x[p], patches.cell_dim_y[p], patches.cell_dim_z[p]};
    leaf.first_cell = patches.first_cell[p];
    if (leaf.dims[0] == 0U || leaf.extent[0] == 0.0) {
      throw std::invalid_argument("patch " + std::to_string(leaf.patch_id) +
                                  " carries the legacy undeclared geometry; rays need restart-authoritative patch "
                                  "geometry");
    }
    if (static_cast<std::uint64_t>(leaf.dims[0]) * leaf.dims[1] * leaf.dims[2] != patches.cell_count[p]) {
      throw std::invalid_argument("patch cell_count disagrees with its cell dimensions");
    }
    leaves.push_back(leaf);
  }
  return LeafGeometry(std::move(leaves), domain_origin, domain_extent, boundary);
}

LeafGeometry LeafGeometry::uniformGrid(std::array<std::uint32_t, 3> dims, std::array<double, 3> domain_origin,
                                       std::array<double, 3> domain_extent, RayBoundary boundary) {
  LeafPatch leaf;
  leaf.patch_id = 0;
  leaf.origin = domain_origin;
  leaf.extent = domain_extent;
  leaf.dims = dims;
  leaf.first_cell = 0;
  return LeafGeometry({leaf}, domain_origin, domain_extent, boundary, 1U);
}

std::size_t LeafGeometry::bucketOf(const std::array<double, 3>& p) const noexcept {
  std::array<std::size_t, 3> b{};
  for (std::size_t a = 0; a < 3U; ++a) {
    const double u = (p[a] - m_origin[a]) / m_extent[a] * static_cast<double>(m_buckets);
    b[a] = static_cast<std::size_t>(std::clamp(std::floor(u), 0.0, static_cast<double>(m_buckets - 1U)));
  }
  return (b[0] * m_buckets + b[1]) * m_buckets + b[2];
}

bool LeafGeometry::locate(const std::array<double, 3>& point, LeafCellRef& out) const {
  const std::size_t b = bucketOf(point);
  for (std::uint32_t i = m_bucket_begin[b]; i < m_bucket_begin[b + 1U]; ++i) {
    const std::uint32_t l = m_bucket_leaves[i];
    const LeafPatch& leaf = m_leaves[l];
    if (!contains(leaf, point)) continue;
    out.leaf = l;
    for (std::size_t a = 0; a < 3U; ++a) {
      const double width = leaf.extent[a] / static_cast<double>(leaf.dims[a]);
      const double u = std::floor((point[a] - leaf.origin[a]) / width);
      out.ijk[a] = static_cast<std::uint32_t>(std::clamp(u, 0.0, static_cast<double>(leaf.dims[a] - 1U)));
    }
    out.cell_row = leaf.first_cell + out.ijk[0] + leaf.dims[0] * (out.ijk[1] + leaf.dims[1] * out.ijk[2]);
    return true;
  }
  return false;
}

std::array<double, 3> LeafGeometry::wrap(const std::array<double, 3>& p) const noexcept {
  std::array<double, 3> q = p;
  for (std::size_t a = 0; a < 3U; ++a) {
    const double rel = (q[a] - m_origin[a]) / m_extent[a];
    q[a] = m_origin[a] + (rel - std::floor(rel)) * m_extent[a];
    if (q[a] >= m_origin[a] + m_extent[a]) q[a] = m_origin[a];
  }
  return q;
}

std::array<double, 3> LeafGeometry::cellCenter(const LeafCellRef& ref) const {
  const LeafPatch& leaf = m_leaves.at(ref.leaf);
  std::array<double, 3> c{};
  for (std::size_t a = 0; a < 3U; ++a) {
    c[a] = leaf.origin[a] + (static_cast<double>(ref.ijk[a]) + 0.5) * leaf.extent[a] / static_cast<double>(leaf.dims[a]);
  }
  return c;
}

double LeafGeometry::cellWidth(std::uint32_t leaf, std::size_t axis) const {
  const LeafPatch& patch = m_leaves.at(leaf);
  return patch.extent.at(axis) / static_cast<double>(patch.dims.at(axis));
}

std::uint64_t LeafGeometry::indexBytes() const noexcept {
  return (m_bucket_begin.size() + m_bucket_leaves.size()) * sizeof(std::uint32_t) + m_leaves.size() * sizeof(LeafPatch);
}

std::vector<RaySegment> LeafGeometry::traceStraightRay(const std::array<double, 3>& start,
                                                       const std::array<double, 3>& d, double max_length,
                                                       std::size_t max_segments) const {
  const double norm = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  if (!(std::abs(norm - 1.0) < 1.0e-12)) throw std::invalid_argument("ray direction must be a unit vector");
  if (!(max_length > 0.0) || max_segments == 0U) throw std::invalid_argument("ray needs positive length and capacity");
  std::vector<RaySegment> segments;
  const double nudge = 1.0e-6 * m_finest_width;
  std::array<double, 3> pos = start;
  if (m_boundary == RayBoundary::kPeriodic) pos = wrap(pos);
  double s = 0.0;
  while (s < max_length) {
    std::array<double, 3> probe{pos[0] + nudge * d[0], pos[1] + nudge * d[1], pos[2] + nudge * d[2]};
    if (m_boundary == RayBoundary::kPeriodic) {
      probe = wrap(probe);
    } else {
      bool outside = false;
      for (std::size_t a = 0; a < 3U; ++a) {
        outside = outside || probe[a] < m_origin[a] || probe[a] >= m_origin[a] + m_extent[a];
      }
      if (outside) break;  // open boundary: the ray leaves the domain
    }
    LeafCellRef ref;
    if (!locate(probe, ref)) {
      throw std::runtime_error("ray probe fell into a hole of the leaf tiling (geometry contract violated)");
    }
    const LeafPatch& leaf = m_leaves[ref.leaf];
    double t_exit = std::numeric_limits<double>::infinity();
    for (std::size_t a = 0; a < 3U; ++a) {
      if (d[a] == 0.0) continue;
      const double width = leaf.extent[a] / static_cast<double>(leaf.dims[a]);
      // Cell faces in the (possibly wrapped) frame of the probe; shift so the
      // distance is measured from `pos` consistently with the probe image.
      const double image_shift = probe[a] - (pos[a] + nudge * d[a]);
      const double lower = leaf.origin[a] + static_cast<double>(ref.ijk[a]) * width - image_shift;
      const double upper = lower + width;
      const double t = d[a] > 0.0 ? (upper - pos[a]) / d[a] : (lower - pos[a]) / d[a];
      t_exit = std::min(t_exit, t);
    }
    if (!(t_exit > 0.0)) t_exit = nudge;  // degenerate tie: advance deterministically
    const double t = std::min(t_exit, max_length - s);
    if (segments.size() >= max_segments) {
      throw std::runtime_error("ray exceeded max_segments = " + std::to_string(max_segments) +
                               " (bounded ray buffer; raise the capacity or shorten the ray)");
    }
    RaySegment segment;
    segment.patch_id = leaf.patch_id;
    segment.cell_row = ref.cell_row;
    segment.leaf = ref.leaf;
    segment.s_begin = s;
    segment.s_end = s + t;
    std::array<double, 3> mid{pos[0] + 0.5 * t * d[0], pos[1] + 0.5 * t * d[1], pos[2] + 0.5 * t * d[2]};
    segment.midpoint = m_boundary == RayBoundary::kPeriodic ? wrap(mid) : mid;
    segments.push_back(segment);
    for (std::size_t a = 0; a < 3U; ++a) pos[a] += t * d[a];
    s += t;
  }
  return segments;
}

std::vector<RaySegment> mergeRankSegments(std::vector<std::vector<RaySegment>> per_rank) {
  std::vector<RaySegment> merged;
  for (auto& list : per_rank) merged.insert(merged.end(), list.begin(), list.end());
  std::sort(merged.begin(), merged.end(), [](const RaySegment& a, const RaySegment& b) {
    if (a.s_begin != b.s_begin) return a.s_begin < b.s_begin;
    if (a.patch_id != b.patch_id) return a.patch_id < b.patch_id;
    return a.cell_row < b.cell_row;
  });
  const double scale = merged.empty() ? 0.0 : merged.back().s_end;
  for (std::size_t i = 1; i < merged.size(); ++i) {
    if (merged[i].s_begin < merged[i - 1U].s_end - 1.0e-12 * scale) {
      throw std::runtime_error("ray segments from different ranks overlap (two owners of one path interval)");
    }
  }
  return merged;
}

std::array<double, 3> unitVectorFromAngles(double theta, double phi) {
  return {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};
}

}  // namespace cosmosim::lightcone
