// I01 single ISM-regime authority fixture.
// NOT COMPILED OR RUN in this handoff (the owner reserves build/test stages).
//
// v10 §4.3: the effective-EOS closure (pressure, signal speed, relaxation)
// and resolved thermochemistry must never act on the same gas. With an
// enabled composition sidecar the per-row ISM regime tag is the single
// authority. This fixture checks the production helper used by the fixed
// grid and the AMR orchestrator (hydro::buildIsmRegimeGate), its agreement
// with the derived P/c_s refresh that feeds the CFL, and that the hydro
// solver honours the gate (gate 0 == no closure, gate 1 == legacy closure).
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/physics/baryon_composition.hpp"

namespace {

using namespace cosmosim;
namespace h = cosmosim::hydro;
namespace p = cosmosim::physics;

constexpr double k_gamma = 5.0 / 3.0;
constexpr std::size_t k_n = 8U;
constexpr std::size_t k_cells = k_n * k_n * k_n;

// Always-effective stub: EOS pressure 4x the ideal value.
class StubEffectiveClosure final : public h::HydroThermodynamicClosure {
 public:
  [[nodiscard]] h::HydroThermodynamicClosureResult evaluate(
      std::size_t, const h::HydroConservedState&, const h::HydroPrimitiveState& ideal, double, double) const override {
    return h::HydroThermodynamicClosureResult{
        .pressure_comoving = 4.0 * ideal.pressure_comoving,
        .signal_speed_squared_code = 4.0 * ideal.signal_speed_squared_code,
        .target_specific_internal_energy_code = ideal.specific_internal_energy_code,
        .uses_effective_ism = true,
        .valid = true};
  }
};

[[nodiscard]] h::HydroPatchGeometry periodicGeometry() {
  auto geometry = h::makeCartesianPatchGeometry({.nx = k_n, .ny = k_n, .nz = k_n,
                                                 .cell_width_x_comoving = 0.125,
                                                 .cell_width_y_comoving = 0.125,
                                                 .cell_width_z_comoving = 0.125});
  h::appendCartesianBoundaryGhostFaces(geometry, h::HydroBoundaryKind::kPeriodic);
  (void)h::completeReconstructionStencils(geometry);
  (void)h::completePeriodicReconstructionStencils(geometry);
  return geometry;
}

[[nodiscard]] core::SimulationState gasState(const h::HydroPatchGeometry& geometry, bool composition) {
  core::SimulationState state;
  state.resizeCells(k_cells);
  for (std::size_t r = 0; r < k_cells; ++r) {
    const auto ijk = geometry.cellIjk(r);
    const double rho = 1.0 + 0.1 * std::sin(0.785398 * static_cast<double>(ijk[0]));
    state.gas_cells.gas_cell_id[r] = 100U + r;
    state.gas_cells.density_code[r] = rho;
    state.cells.mass_code[r] = rho * geometry.cell_volume_comoving;
    state.gas_cells.internal_energy_code[r] = 10.0 + static_cast<double>(ijk[1]);
    state.gas_cells.velocity_x_peculiar[r] = 0.2;
  }
  if (composition) {
    p::enableGasComposition(state, p::CompositionLayout::make({.network = p::CompositionNetworkKind::kHybridLight}),
                            {});
    // Effective regime on even x columns only.
    for (std::size_t r = 0; r < k_cells; ++r) {
      state.gas_composition.ism_regime_tag[r] = (geometry.cellIjk(r)[0] % 2U == 0U)
          ? core::k_gas_ism_regime_effective
          : static_cast<std::uint8_t>(0U);
    }
  }
  return state;
}

void checkGateHelpers() {
  const auto geometry = periodicGeometry();
  std::vector<std::uint32_t> dense(k_cells);
  std::iota(dense.begin(), dense.end(), 0U);
  std::vector<std::uint8_t> gate{7U};

  // No sidecar: no tag, legacy closure decision, no gate.
  const auto legacy = gasState(geometry, false);
  assert(!h::ismRegimeGateActive(legacy));
  assert(h::ismRegimeAllowsEffectiveClosure(legacy, 3U) == 1U);
  assert(!h::buildIsmRegimeGate(legacy, geometry, dense, gate));
  assert(gate.empty());

  // Sidecar: the tag decides for real rows and periodic ghosts.
  const auto state = gasState(geometry, true);
  assert(h::ismRegimeGateActive(state));
  assert(h::ismRegimeAllowsEffectiveClosure(state, k_cells + 5U) == 0U);  // out of range: never effective
  assert(h::buildIsmRegimeGate(state, geometry, dense, gate));
  assert(gate.size() == geometry.totalCellStorageCount());
  for (std::size_t r = 0; r < k_cells; ++r) {
    const bool effective = state.gas_composition.ism_regime_tag[r] == core::k_gas_ism_regime_effective;
    assert(gate[r] == (effective ? 1U : 0U));
  }
  std::size_t checked_ghosts = 0U;
  for (const auto& ghost : geometry.ghost_cells) {
    if (ghost.source_real_cell < k_cells && ghost.ghost_cell < gate.size() &&
        ghost.boundary_kind == h::HydroBoundaryKind::kPeriodic) {
      assert(gate[ghost.ghost_cell] == gate[ghost.source_real_cell]);
      ++checked_ghosts;
    }
  }
  assert(checked_ghosts != 0U);
}

void checkRefreshAgreesWithGate() {
  const auto geometry = periodicGeometry();
  auto state = gasState(geometry, true);
  const StubEffectiveClosure closure;
  const h::GasThermodynamicRefreshPolicy policy{.adiabatic_index = k_gamma, .density_floor_code = 1.0e-14,
                                                .pressure_floor_code = 1.0e-14, .scale_factor = 1.0,
                                                .redshift = 0.0, .closure = &closure};
  for (const std::size_t row : {std::size_t{0U}, std::size_t{1U}}) {  // x = 0 effective, x = 1 resolved
    (void)h::refreshGasThermodynamicCell(state, row, policy);
    const double ideal = (k_gamma - 1.0) * state.gas_cells.density_code[row] * state.gas_cells.internal_energy_code[row];
    const bool effective = h::ismRegimeAllowsEffectiveClosure(state, row) != 0U;
    assert(std::abs(state.gas_cells.pressure_code[row] - (effective ? 4.0 : 1.0) * ideal) <= 1.0e-12 * ideal);
  }
}

[[nodiscard]] h::HydroConservedStateSoa advance(const h::HydroPatchGeometry& geometry,
                                                const core::SimulationState& state,
                                                const h::HydroThermodynamicClosure* closure,
                                                std::span<const std::uint8_t> gate) {
  h::HydroConservedStateSoa soa(geometry.totalCellStorageCount());
  for (std::size_t r = 0; r < k_cells; ++r) {
    const double rho = state.gas_cells.density_code[r];
    const double vx = state.gas_cells.velocity_x_peculiar[r];
    soa.storeCell(r, h::HydroConservedState{
        .mass_density_comoving = rho,
        .momentum_density_x_comoving = rho * vx,
        .total_energy_density_comoving = rho * (state.gas_cells.internal_energy_code[r] + 0.5 * vx * vx)});
  }
  h::fillHydroBoundaryGhostCells(soa, geometry, k_gamma);
  h::HydroNumericsPolicy numerics;
  numerics.integration_scheme = h::HydroIntegrationScheme::kCornerTransportUpwind;
  const h::HydroCoreSolver solver(k_gamma, numerics);
  const h::HydroUpdateContext update{.dt_code = 2.0e-3, .time_begin_code = 0.0, .time_end_code = 2.0e-3};
  const auto courant = h::hydroTransportCourantFactors(update.dt_code, 1.0, 0.125, 0.125, 0.125);
  const h::PpmReconstruction ppm(h::HydroReconstructionPolicy{
      .dt_over_dx_code = courant[0], .dt_over_cell_width_code = courant, .adiabatic_index = k_gamma});
  const h::HllcRiemannSolver riemann;
  std::vector<std::size_t> cells(k_cells);
  std::iota(cells.begin(), cells.end(), 0U);
  std::vector<std::size_t> faces(geometry.faces.size());
  std::iota(faces.begin(), faces.end(), 0U);
  h::HydroScratchBuffers scratch;
  h::HydroSourceContext context{.update = update};
  context.thermodynamic_closure = closure;
  context.effective_closure_allowed = gate;
  solver.advancePatchActiveSetBatchedWithScratch(soa, geometry, {.active_cells = cells, .active_faces = faces},
                                                 update, ppm, riemann, {}, context, scratch, nullptr,
                                                 solver.requiredWorkspaceFeatures({}));
  return soa;
}

[[nodiscard]] bool sameRealCells(const h::HydroConservedStateSoa& a, const h::HydroConservedStateSoa& b) {
  for (std::size_t r = 0; r < k_cells; ++r) {
    const auto x = a.loadCell(r);
    const auto y = b.loadCell(r);
    if (x.mass_density_comoving != y.mass_density_comoving ||
        x.momentum_density_x_comoving != y.momentum_density_x_comoving ||
        x.total_energy_density_comoving != y.total_energy_density_comoving) {
      return false;
    }
  }
  return true;
}

void checkSolverHonoursGate() {
  const auto geometry = periodicGeometry();
  const auto state = gasState(geometry, true);
  const StubEffectiveClosure closure;
  const std::vector<std::uint8_t> none(geometry.totalCellStorageCount(), 0U);
  const std::vector<std::uint8_t> all(geometry.totalCellStorageCount(), 1U);
  // Gate 0 everywhere == no closure; gate 1 everywhere == legacy (no gate).
  assert(sameRealCells(advance(geometry, state, &closure, none), advance(geometry, state, nullptr, {})));
  assert(sameRealCells(advance(geometry, state, &closure, all), advance(geometry, state, &closure, {})));
  // The production gate from the tags differs from both extremes.
  std::vector<std::uint32_t> dense(k_cells);
  std::iota(dense.begin(), dense.end(), 0U);
  std::vector<std::uint8_t> gate;
  assert(h::buildIsmRegimeGate(state, geometry, dense, gate));
  const auto tagged = advance(geometry, state, &closure, gate);
  assert(!sameRealCells(tagged, advance(geometry, state, &closure, {})));
  assert(!sameRealCells(tagged, advance(geometry, state, nullptr, {})));
}

}  // namespace

int main() {
  checkGateHelpers();
  checkRefreshAgreesWithGate();
  checkSolverHonoursGate();
  return 0;
}
