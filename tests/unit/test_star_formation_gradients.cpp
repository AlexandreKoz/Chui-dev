// CHUI Part 2 corrective gate — star-formation velocity gradients across
// patch boundaries, periodic boundaries, coarse/fine interfaces and missing
// remote neighbours (NOT RUN — USER-RESERVED).
//
// Reference: a linear field v_x = 3 x has dv_x/dx = 3 exactly for any
// centred or one-sided difference with the true centre separations; the
// periodic case uses the minimum-image separation. A neighbour that is
// missing (not owned locally) and not a physical boundary must make the
// gradient UNAVAILABLE (NaN), never zero.
#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "workflows/internal/star_formation_gradients.hpp"

namespace {

using cosmosim::core::SimulationState;
using cosmosim::workflows::internal::StarFormationNeighborIndex;

struct PatchSpec {
  int level;
  double x0;
  double width;   // cell width along x
  std::uint32_t nx;
};

// Builds patches stacked along x; y and z are one cell of width 1.
SimulationState makeState(const std::vector<PatchSpec>& patches, double (*vx)(double)) {
  std::uint32_t cells = 0;
  for (const auto& p : patches) cells += p.nx;
  SimulationState state;
  state.resizeCells(cells);
  state.resizePatches(patches.size());
  std::uint32_t row = 0;
  for (std::size_t i = 0; i < patches.size(); ++i) {
    const auto& p = patches[i];
    state.patches.patch_id[i] = 100U + i;
    state.patches.level[i] = p.level;
    state.patches.first_cell[i] = row;
    state.patches.cell_count[i] = p.nx;
    state.patches.origin_x_comoving[i] = p.x0;
    state.patches.extent_x_comoving[i] = p.width * p.nx;
    state.patches.extent_y_comoving[i] = 1.0;
    state.patches.extent_z_comoving[i] = 1.0;
    state.patches.cell_dim_x[i] = p.nx;
    state.patches.cell_dim_y[i] = 1;
    state.patches.cell_dim_z[i] = 1;
    for (std::uint32_t c = 0; c < p.nx; ++c, ++row) {
      const double x = p.x0 + (c + 0.5) * p.width;
      state.cells.patch_index[row] = static_cast<std::uint32_t>(i);
      state.cells.center_x_comoving[row] = x;
      state.cells.center_y_comoving[row] = 0.5;
      state.cells.center_z_comoving[row] = 0.5;
      state.gas_cells.velocity_x_peculiar[row] = vx(x);
    }
  }
  return state;
}

double linear(double x) { return 3.0 * x; }
double wave(double x) { return std::sin(2.0 * 3.14159265358979323846 * x / 4.0); }

void testAcrossSameLevelPatchBoundary() {
  const auto state = makeState({{0, 0.0, 1.0, 2}, {0, 2.0, 1.0, 2}}, linear);
  const StarFormationNeighborIndex index(state, false, {4.0, 1.0, 1.0}, {});
  const auto g = index.velocityGradient(1, 1.0);  // x = 1.5, right neighbour in the other patch
  assert(g.valid);
  assert(std::abs(g.velocity[0][0] - 3.0) < 1.0e-12);
  assert(!g.used_one_sided);
  assert(g.velocity[0][1] == 0.0 && g.velocity[0][2] == 0.0);  // one-cell physical dims
  // Physical (non-periodic) box edge: one-sided, still exact for a linear field.
  const auto edge = index.velocityGradient(0, 1.0);
  assert(edge.valid && edge.used_one_sided);
  assert(std::abs(edge.velocity[0][0] - 3.0) < 1.0e-12);
  // Physical spacing factor (comoving -> physical) divides the derivative.
  const auto scaled = index.velocityGradient(1, 0.5);
  assert(std::abs(scaled.velocity[0][0] - 6.0) < 1.0e-12);
}

void testPeriodicWrap() {
  const auto state = makeState({{0, 0.0, 1.0, 4}}, wave);
  const StarFormationNeighborIndex index(state, true, {4.0, 1.0, 1.0}, {});
  const auto g = index.velocityGradient(0, 1.0);  // x = 0.5: neighbours 3.5 (wrapped) and 1.5
  assert(g.valid && !g.used_one_sided);
  const double expected = (wave(1.5) - wave(3.5)) / 2.0;
  assert(std::abs(g.velocity[0][0] - expected) < 1.0e-12);
}

void testCoarseFineInterface() {
  // Level-0 cells at 0.5 and 1.5 (width 1); level-1 cells at 2.25, 2.75 (width 0.5).
  const auto state = makeState({{0, 0.0, 1.0, 2}, {1, 2.0, 0.5, 2}}, linear);
  const StarFormationNeighborIndex index(state, false, {3.0, 1.0, 1.0}, {});
  const auto g = index.velocityGradient(2, 1.0);  // fine cell at 2.25
  assert(g.valid && g.used_coarser_neighbor);
  // (v(2.75) - v(1.5)) / (2.75 - 1.5) = 3 with the true unequal spacing.
  assert(std::abs(g.velocity[0][0] - 3.0) < 1.0e-12);
}

void testMissingRemoteNeighbourIsUnavailable() {
  const auto state = makeState({{0, 0.0, 1.0, 2}, {0, 2.0, 1.0, 2}}, linear);
  // Patch B is "remote": excluded from the owned-leaf mask.
  const std::vector<std::uint8_t> owned{1U, 1U, 0U, 0U};
  const StarFormationNeighborIndex index(state, false, {4.0, 1.0, 1.0}, owned);
  const auto g = index.velocityGradient(1, 1.0);
  assert(!g.valid);
  assert(std::isnan(g.velocity[0][0]));  // never silently zero
}

}  // namespace

int main() {
  testAcrossSameLevelPatchBoundary();
  testPeriodicWrap();
  testCoarseFineInterface();
  testMissingRemoteNeighbourIsUnavailable();
  return 0;
}
