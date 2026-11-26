#include "cosmosim/physics/black_hole/bh_neighborhood.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <type_traits>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

static_assert(std::is_trivially_copyable_v<BhSupportQuery>, "BhSupportQuery is a wire record");
static_assert(std::is_trivially_copyable_v<BhSupportReply>, "BhSupportReply is a wire record");

bool bhCellHalfWidthStored(const core::SimulationState& state, std::uint32_t row,
                           std::array<double, 3>& hw) noexcept {
  if (row >= state.cells.size()) return false;
  if (row < state.cells.patch_index.size()) {
    const std::uint32_t p = state.cells.patch_index[row];
    if (p < state.patches.size()) {
      const std::uint32_t nx = state.patches.cell_dim_x[p];
      const std::uint32_t ny = state.patches.cell_dim_y[p];
      const std::uint32_t nz = state.patches.cell_dim_z[p];
      const double ex = state.patches.extent_x_comoving[p];
      const double ey = state.patches.extent_y_comoving[p];
      const double ez = state.patches.extent_z_comoving[p];
      if (nx != 0U && ny != 0U && nz != 0U && ex > 0.0 && ey > 0.0 && ez > 0.0) {
        hw = {0.5 * ex / nx, 0.5 * ey / ny, 0.5 * ez / nz};
        return true;
      }
    }
  }
  const double rho = row < state.gas_cells.size() ? state.gas_cells.density_code[row] : 0.0;
  const double m = state.cells.mass_code[row];
  if (rho > 0.0 && m > 0.0) {
    const double side = std::cbrt(m / rho);
    hw = {0.5 * side, 0.5 * side, 0.5 * side};
    return true;
  }
  return false;
}

std::vector<BhSupportQuery> LocalBhSupportExchange::routeQueries(std::span<const BhSupportQuery> q) {
  std::vector<BhSupportQuery> out(q.begin(), q.end());
  std::sort(out.begin(), out.end(), [](const BhSupportQuery& a, const BhSupportQuery& b) {
    return std::tie(a.origin_rank, a.query_index) < std::tie(b.origin_rank, b.query_index);
  });
  return out;
}

std::vector<BhSupportReply> LocalBhSupportExchange::returnReplies(std::span<const BhSupportReply> r) {
  std::vector<BhSupportReply> out(r.begin(), r.end());
  std::sort(out.begin(), out.end(), [](const BhSupportReply& a, const BhSupportReply& b) {
    return std::tie(a.query_index, a.responder_rank) < std::tie(b.query_index, b.responder_rank);
  });
  return out;
}

BhProximityFilter::BhProximityFilter(std::span<const BhSupportQuery> queries, const PeriodicBox& box)
    : m_queries(queries), m_box(box) {
  double r_max = 0.0;
  for (const BhSupportQuery& q : queries) {
    if (!(q.support_radius_stored >= 0.0) || !std::isfinite(q.support_radius_stored)) {
      throw std::invalid_argument("BhProximityFilter: invalid query support radius");
    }
    r_max = std::max(r_max, q.support_radius_stored);
  }
  m_bucket = r_max > 0.0 ? r_max : 1.0;
  if (m_box.periodic) {
    for (std::size_t a = 0; a < 3; ++a) {
      if (!(m_box.extent[a] > 0.0)) throw std::invalid_argument("BhProximityFilter: periodic box needs extents");
      // Bucket count per axis; a single bucket when the support spans the box.
      m_bucket_count[a] = std::max<std::int64_t>(1, static_cast<std::int64_t>(std::floor(m_box.extent[a] / m_bucket)));
    }
  }
  for (std::uint32_t i = 0; i < queries.size(); ++i) {
    const auto b = bucketOf(m_box.wrap(queries[i].position_stored));
    m_buckets[bucketKey(b[0], b[1], b[2])].push_back(i);
  }
}

std::array<std::int64_t, 3> BhProximityFilter::bucketOf(const std::array<double, 3>& p) const noexcept {
  std::array<std::int64_t, 3> b{};
  for (std::size_t a = 0; a < 3; ++a) {
    if (m_box.periodic) {
      const double cell = m_box.extent[a] / static_cast<double>(m_bucket_count[a]);
      b[a] = std::clamp<std::int64_t>(static_cast<std::int64_t>(std::floor(p[a] / cell)), 0, m_bucket_count[a] - 1);
    } else {
      b[a] = static_cast<std::int64_t>(std::floor(p[a] / m_bucket));
    }
  }
  return b;
}

std::int64_t BhProximityFilter::bucketKey(std::int64_t ix, std::int64_t iy, std::int64_t iz) const noexcept {
  if (m_box.periodic) {
    ix = ((ix % m_bucket_count[0]) + m_bucket_count[0]) % m_bucket_count[0];
    iy = ((iy % m_bucket_count[1]) + m_bucket_count[1]) % m_bucket_count[1];
    iz = ((iz % m_bucket_count[2]) + m_bucket_count[2]) % m_bucket_count[2];
  }
  // 21-bit signed fields: sufficient for any domain / bucket ratio below 2^20.
  constexpr std::int64_t k_mask = (1LL << 21) - 1;
  return ((ix & k_mask) << 42) | ((iy & k_mask) << 21) | (iz & k_mask);
}

void BhProximityFilter::matches(const std::array<double, 3>& position, std::vector<std::uint32_t>& out) const {
  out.clear();
  if (m_queries.empty()) return;
  const auto b = bucketOf(m_box.wrap(position));
  // With fewer than 3 periodic buckets on an axis the 27-neighbourhood wraps
  // onto itself; deduplicate below.
  for (std::int64_t dx = -1; dx <= 1; ++dx) {
    for (std::int64_t dy = -1; dy <= 1; ++dy) {
      for (std::int64_t dz = -1; dz <= 1; ++dz) {
        const auto it = m_buckets.find(bucketKey(b[0] + dx, b[1] + dy, b[2] + dz));
        if (it == m_buckets.end()) continue;
        for (const std::uint32_t qi : it->second) {
          const BhSupportQuery& q = m_queries[qi];
          const auto d = m_box.delta(q.position_stored, position);
          if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= q.support_radius_stored * q.support_radius_stored) {
            out.push_back(qi);
          }
        }
      }
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
}

}  // namespace cosmosim::physics
