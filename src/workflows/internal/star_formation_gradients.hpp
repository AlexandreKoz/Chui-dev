#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::workflows::internal {

// ---------------------------------------------------------------------------
// CHUI Part 2 audit repair (supplementary risk "star-formation boundary
// gradients"). The former per-patch derivative assumed dense storage order
// equalled patch-local order and, at patch edges, silently used one-sided
// in-patch differences even when a same-level neighbour patch, a periodic
// image or a coarser cell existed.
//
// This index locates the true face neighbour of any local gas cell by
// GEOMETRY: every local leaf cell is keyed by (level, integer lattice
// coordinates of its centre at that level's cell width). The neighbour
// along +-axis is searched at the same level, then at the next coarser level
// (coarse/fine interface). Periodic modes wrap the query point into the box.
// Differences use the actual centre separation, so unequal spacings across
// coarse/fine interfaces are handled. Outcomes per side:
//   found                -> used;
//   outside a physical (open/reflective) boundary -> one-sided difference;
//   missing otherwise (remote rank / unsupported topology) -> the gradient is
//   reported UNAVAILABLE (NaN); it is never treated as zero.
// Memory: one sorted (key,row) vector of 24 B per local cell per stage.
// ---------------------------------------------------------------------------
class StarFormationNeighborIndex {
 public:
  struct Entry {
    std::uint64_t key = 0;
    std::uint32_t row = 0;
    std::uint32_t level = 0;
  };

  StarFormationNeighborIndex(const core::SimulationState& state, bool periodic, std::array<double, 3> box_extent,
                             std::span<const std::uint8_t> leaf_mask)
      : m_state(state), m_periodic(periodic), m_box(box_extent) {
    m_entries.reserve(state.cells.size());
    for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
      if (!leaf_mask.empty() && (row >= leaf_mask.size() || leaf_mask[row] == 0U)) continue;
      const auto width = cellWidth(row);
      if (!(width[0] > 0.0)) continue;
      const std::uint32_t level = levelOf(row);
      const std::array<double, 3> center{state.cells.center_x_comoving[row], state.cells.center_y_comoving[row],
                                         state.cells.center_z_comoving[row]};
      m_entries.push_back(Entry{latticeKey(level, center, width), row, level});
      m_width_by_level_x.resize(std::max<std::size_t>(m_width_by_level_x.size(), level + 1U), 0.0);
      m_width_by_level_x[level] = width[0];
      m_width_by_level_y.resize(m_width_by_level_x.size(), 0.0);
      m_width_by_level_z.resize(m_width_by_level_x.size(), 0.0);
      m_width_by_level_y[level] = width[1];
      m_width_by_level_z[level] = width[2];
    }
    std::sort(m_entries.begin(), m_entries.end(),
              [](const Entry& a, const Entry& b) { return a.key < b.key || (a.key == b.key && a.row < b.row); });
  }

  struct Gradient {
    std::array<std::array<double, 3>, 3> velocity{};
    bool valid = false;
    bool used_coarser_neighbor = false;
    bool used_one_sided = false;
  };

  // Velocity gradient (physical spacing = comoving * length_to_physical).
  [[nodiscard]] Gradient velocityGradient(std::uint32_t row, double length_to_physical) const {
    Gradient g;
    const auto width = cellWidth(row);
    if (!(width[0] > 0.0)) return g;
    const std::array<std::span<const double>, 3> v{m_state.gas_cells.velocity_x_peculiar,
                                                   m_state.gas_cells.velocity_y_peculiar,
                                                   m_state.gas_cells.velocity_z_peculiar};
    const std::array<double, 3> center{m_state.cells.center_x_comoving[row], m_state.cells.center_y_comoving[row],
                                       m_state.cells.center_z_comoving[row]};
    const std::uint32_t level = levelOf(row);
    for (int axis = 0; axis < 3; ++axis) {
      double x_minus = 0.0, x_plus = 0.0;
      std::uint32_t r_minus = core::kInvalidGasCellRow, r_plus = core::kInvalidGasCellRow;
      bool minus_physical = false, plus_physical = false;
      for (int side = -1; side <= 1; side += 2) {
        std::array<double, 3> query = center;
        query[axis] += side * width[axis];
        bool outside_physical = false;
        if (!wrapOrFlag(query, axis, outside_physical)) {
          (side < 0 ? minus_physical : plus_physical) = outside_physical;
          continue;
        }
        bool coarser = false;
        const std::uint32_t nb = find(query, level, coarser);
        if (nb == core::kInvalidGasCellRow) continue;
        g.used_coarser_neighbor = g.used_coarser_neighbor || coarser;
        // Separation along the axis with the periodic minimum image.
        double dxc = coordinate(nb, axis) - center[axis];
        if (m_periodic && m_box[axis] > 0.0) dxc -= m_box[axis] * std::round(dxc / m_box[axis]);
        if (side < 0) { r_minus = nb; x_minus = dxc; } else { r_plus = nb; x_plus = dxc; }
      }
      for (int c = 0; c < 3; ++c) {
        double derivative = std::numeric_limits<double>::quiet_NaN();
        const double f0 = v[c][row];
        if (r_minus != core::kInvalidGasCellRow && r_plus != core::kInvalidGasCellRow && x_plus > x_minus) {
          derivative = (v[c][r_plus] - v[c][r_minus]) / ((x_plus - x_minus) * length_to_physical);
        } else if (r_plus != core::kInvalidGasCellRow && minus_physical && x_plus > 0.0) {
          derivative = (v[c][r_plus] - f0) / (x_plus * length_to_physical);
          g.used_one_sided = true;
        } else if (r_minus != core::kInvalidGasCellRow && plus_physical && x_minus < 0.0) {
          derivative = (f0 - v[c][r_minus]) / (x_minus * length_to_physical);
          g.used_one_sided = true;
        } else if (axisCollapsed(axis)) {
          derivative = 0.0;  // genuinely lower-dimensional domain along this axis
        }
        g.velocity[c][axis] = derivative;
      }
    }
    g.valid = true;
    for (const auto& row_g : g.velocity) {
      for (const double d : row_g) g.valid = g.valid && std::isfinite(d);
    }
    return g;
  }

  [[nodiscard]] std::size_t entryBytes() const noexcept { return m_entries.capacity() * sizeof(Entry); }

 private:
  [[nodiscard]] std::array<double, 3> cellWidth(std::uint32_t row) const {
    if (row >= m_state.cells.patch_index.size()) return {0.0, 0.0, 0.0};
    const std::uint32_t p = m_state.cells.patch_index[row];
    if (p >= m_state.patches.size()) return {0.0, 0.0, 0.0};
    const auto& patches = m_state.patches;
    if (patches.cell_dim_x[p] == 0U || patches.cell_dim_y[p] == 0U || patches.cell_dim_z[p] == 0U) return {0.0, 0.0, 0.0};
    return {patches.extent_x_comoving[p] / patches.cell_dim_x[p], patches.extent_y_comoving[p] / patches.cell_dim_y[p],
            patches.extent_z_comoving[p] / patches.cell_dim_z[p]};
  }
  [[nodiscard]] std::uint32_t levelOf(std::uint32_t row) const {
    const std::uint32_t p = m_state.cells.patch_index[row];
    return p < m_state.patches.size() ? static_cast<std::uint32_t>(std::max(m_state.patches.level[p], 0)) : 0U;
  }
  [[nodiscard]] double coordinate(std::uint32_t row, int axis) const {
    return axis == 0 ? m_state.cells.center_x_comoving[row]
                     : (axis == 1 ? m_state.cells.center_y_comoving[row] : m_state.cells.center_z_comoving[row]);
  }
  [[nodiscard]] bool axisCollapsed(int axis) const {
    // A periodic box axis covered by a single cell at level 0 has no gradient.
    return m_periodic && axis < 3 && !m_width_by_level_x.empty() &&
        std::abs((axis == 0 ? m_width_by_level_x[0] : axis == 1 ? m_width_by_level_y[0] : m_width_by_level_z[0]) -
                 m_box[axis]) <= 1.0e-12 * m_box[axis];
  }
  // Periodic: wrap into [0, L). Non-periodic: false when the query is
  // outside the box (physical boundary) — flagged so one-sided is legal.
  [[nodiscard]] bool wrapOrFlag(std::array<double, 3>& q, int axis, bool& outside_physical) const {
    if (m_box[axis] <= 0.0) return true;
    if (m_periodic) {
      q[axis] -= m_box[axis] * std::floor(q[axis] / m_box[axis]);
      return true;
    }
    if (q[axis] < 0.0 || q[axis] >= m_box[axis]) {
      outside_physical = true;
      return false;
    }
    return true;
  }
  [[nodiscard]] static std::uint64_t pack(std::uint32_t level, std::int64_t i, std::int64_t j, std::int64_t k) {
    const auto u = [](std::int64_t x) { return static_cast<std::uint64_t>(x + (1LL << 19)) & 0xFFFFFULL; };
    return (static_cast<std::uint64_t>(level & 0xFU) << 60U) | (u(i) << 40U) | (u(j) << 20U) | u(k);
  }
  [[nodiscard]] std::uint64_t latticeKey(std::uint32_t level, const std::array<double, 3>& p,
                                         const std::array<double, 3>& w) const {
    return pack(level, static_cast<std::int64_t>(std::floor(p[0] / w[0])),
                static_cast<std::int64_t>(std::floor(p[1] / w[1])), static_cast<std::int64_t>(std::floor(p[2] / w[2])));
  }
  [[nodiscard]] std::uint32_t lookup(std::uint64_t key) const {
    const auto it = std::lower_bound(m_entries.begin(), m_entries.end(), key,
                                     [](const Entry& e, std::uint64_t k) { return e.key < k; });
    return (it != m_entries.end() && it->key == key) ? it->row : core::kInvalidGasCellRow;
  }
  [[nodiscard]] std::uint32_t find(const std::array<double, 3>& q, std::uint32_t level, bool& coarser) const {
    coarser = false;
    for (std::int64_t lv = level; lv >= 0 && lv + 2 > static_cast<std::int64_t>(level); --lv) {
      const std::size_t l = static_cast<std::size_t>(lv);
      if (l >= m_width_by_level_x.size() || !(m_width_by_level_x[l] > 0.0)) continue;
      const std::array<double, 3> w{m_width_by_level_x[l], m_width_by_level_y[l], m_width_by_level_z[l]};
      const std::uint32_t row = lookup(latticeKey(static_cast<std::uint32_t>(l), q, w));
      if (row != core::kInvalidGasCellRow) {
        coarser = l < level;
        return row;
      }
    }
    return core::kInvalidGasCellRow;
  }

  const core::SimulationState& m_state;
  bool m_periodic = false;
  std::array<double, 3> m_box{};
  std::vector<Entry> m_entries;
  std::vector<double> m_width_by_level_x, m_width_by_level_y, m_width_by_level_z;
};

}  // namespace cosmosim::workflows::internal
