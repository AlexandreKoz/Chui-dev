#include "cosmosim/radiation/radiation_amr.hpp"

#include <algorithm>
#include <cmath>

namespace cosmosim::radiation {
namespace {

[[nodiscard]] double mc(double l, double m, double r) noexcept {
  const double a = m - l;
  const double b = r - m;
  if (a * b <= 0.0) return 0.0;
  const double s = a > 0.0 ? 1.0 : -1.0;
  return s * std::min({2.0 * std::abs(a), 2.0 * std::abs(b), 0.5 * std::abs(a + b)});
}

}  // namespace

M1State restrictRadiationChildren(const std::array<M1State, 8>& ch) noexcept {
  M1State p{};
  for (const M1State& c : ch) {
    p.energy += c.energy / 8.0;
    for (int d = 0; d < 3; ++d) p.flux[d] += c.flux[d] / 8.0;
  }
  return p;
}

std::array<M1State, 8> prolongRadiationParent(const M1State& parent,
                                              const std::array<std::array<M1State, 2>, 3>& nb, double c) noexcept {
  // Per-axis slopes (cell-width units) of (E, F).
  std::array<M1State, 3> slope{};
  for (int a = 0; a < 3; ++a) {
    slope[a].energy = mc(nb[a][0].energy, parent.energy, nb[a][1].energy);
    for (int d = 0; d < 3; ++d) slope[a].flux[d] = mc(nb[a][0].flux[d], parent.flux[d], nb[a][1].flux[d]);
  }
  std::array<M1State, 8> out{};
  std::array<M1State, 8> delta{};
  for (int cz = 0; cz < 2; ++cz) {
    for (int cy = 0; cy < 2; ++cy) {
      for (int cx = 0; cx < 2; ++cx) {
        const int idx = cx + 2 * (cy + 2 * cz);
        const double sx = cx == 0 ? -0.25 : 0.25;
        const double sy = cy == 0 ? -0.25 : 0.25;
        const double sz = cz == 0 ? -0.25 : 0.25;
        delta[idx].energy = sx * slope[0].energy + sy * slope[1].energy + sz * slope[2].energy;
        for (int d = 0; d < 3; ++d) {
          delta[idx].flux[d] = sx * slope[0].flux[d] + sy * slope[1].flux[d] + sz * slope[2].flux[d];
        }
      }
    }
  }
  double theta = 1.0;
  for (const M1State& dl : delta) theta = std::min(theta, m1AdmissibleFraction(parent, dl, c));
  for (int i = 0; i < 8; ++i) {
    out[i].energy = parent.energy + theta * delta[i].energy;
    for (int d = 0; d < 3; ++d) out[i].flux[d] = parent.flux[d] + theta * delta[i].flux[d];
  }
  return out;  // symmetric deltas sum to zero: children average to the parent exactly
}

}  // namespace cosmosim::radiation
