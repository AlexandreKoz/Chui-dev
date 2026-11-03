// CHUI Part 2 — conservative composition advection (NOT RUN — USER-RESERVED).
//
// Invariants:
//  * uniform-velocity / uniform-pressure contact on a periodic 4-cell patch
//    (HLLC exact => donor cell): every element lane moves with dm * X_donor,
//    crossing the periodic box face exactly like an interior face;
//  * every lane is globally conserved to round-off;
//  * after commit the composition mass equals the hydro cell mass
//    (consistency: sum_s X_s = 1) — no reconciliation mass needed;
//  * application is independent of record order (step-start donors);
//  * an imported MPI ghost refuses (Part 3 boundary) instead of dropping lanes;
//  * coarse/fine: a kCoarse face does not update the coarse cell; the fine
//    face mirrors exactly its own transfer onto the coarse cell (reflux).
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include "cosmosim/amr/amr_hydro_geometry.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/composition_transport.hpp"
#include "../support/test_require.hpp"

namespace {

namespace hy = cosmosim::hydro;
namespace ph = cosmosim::physics;
constexpr double k_gamma = 5.0 / 3.0;
constexpr std::uint32_t k_nx = 4U;

cosmosim::core::SimulationState makeState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(k_nx);
  state.resizePatches(1);
  state.patches.patch_id[0] = 7;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = k_nx;
  state.patches.morton_key[0] = 7;
  state.patches.extent_x_comoving[0] = 1.0;
  state.patches.extent_y_comoving[0] = 0.25;
  state.patches.extent_z_comoving[0] = 0.25;
  state.patches.cell_dim_x[0] = k_nx;
  state.patches.cell_dim_y[0] = 1;
  state.patches.cell_dim_z[0] = 1;
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  for (std::uint32_t row = 0; row < k_nx; ++row) {
    const double rho = 1.0 + 0.25 * static_cast<double>(row);
    state.cells.center_x_comoving[row] = (static_cast<double>(row) + 0.5) * 0.25;
    state.cells.center_y_comoving[row] = 0.125;
    state.cells.center_z_comoving[row] = 0.125;
    state.cells.mass_code[row] = rho * 0.25 * 0.25 * 0.25;
    state.gas_cells.gas_cell_id[row] = 500U + row;
    state.gas_cells.density_code[row] = rho;
    state.gas_cells.internal_energy_code[row] = 1.0 / ((k_gamma - 1.0) * rho);
    state.gas_cells.pressure_code[row] = 1.0;
    state.gas_cells.velocity_x_peculiar[row] = 0.5;
    // Strongly varying metallicity so donor selection is visible.
    state.gas_cells.metal_mass_code[row] = (0.001 + 0.01 * row) * state.cells.mass_code[row];
    records.push_back({.gas_cell_id = 500U + row, .parent_particle_id = std::nullopt,
                       .owning_patch_id = 7, .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

cosmosim::amr::PatchDescriptor descriptor() {
  cosmosim::amr::PatchDescriptor patch;
  patch.patch_id = 7;
  patch.origin_comov = {0.0, 0.0, 0.0};
  patch.extent_comov = {1.0, 0.25, 0.25};
  patch.cell_dims = {static_cast<std::uint16_t>(k_nx), 1, 1};
  return patch;
}

void testPeriodicDonorCellCompositionAdvection() {
  auto state = makeState();
  const auto layout = ph::CompositionLayout::make(
      {.network = ph::CompositionNetworkKind::kHybridLight});
  ph::enableGasComposition(state, layout, {});

  cosmosim::amr::AmrHydroGeometryOptions options;
  options.physical_boundary_kind = hy::HydroBoundaryKind::kPeriodic;
  const auto patch_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(state, descriptor(), options);
  const hy::HydroPatchGeometry& geometry = patch_geometry.geometry;
  hy::HydroConservedStateSoa conserved = cosmosim::amr::loadAmrHydroConservedState(state, patch_geometry, k_gamma);
  hy::fillHydroBoundaryGhostCells(conserved, geometry, k_gamma);
  std::vector<std::size_t> active_cells;
  for (std::size_t c = 0; c < geometry.cellCount(); ++c) active_cells.push_back(c);
  std::vector<std::size_t> active_faces;
  for (std::size_t f = 0; f < geometry.faces.size(); ++f) {
    const hy::HydroFace& face = geometry.faces[f];
    if (face.ghost_cell_slot != hy::k_invalid_ghost_cell_slot &&
        geometry.ghost_cells[face.ghost_cell_slot].side == hy::HydroFaceSide::kLower) {
      continue;  // box face evaluated once (upper copy)
    }
    active_faces.push_back(f);
  }
  std::vector<std::uint32_t> dense(k_nx);
  for (std::uint32_t r = 0; r < k_nx; ++r) dense[r] = r;

  const std::uint32_t o_lane = layout.elementLane(ph::ChemicalElement::kO);
  std::vector<double> o_before(k_nx), mass_before(k_nx);
  double o_total_before = 0.0;
  for (std::uint32_t r = 0; r < k_nx; ++r) {
    o_before[r] = state.gas_composition.value(o_lane, r);
    mass_before[r] = state.cells.mass_code[r];
    o_total_before += o_before[r];
  }

  ph::CompositionFluxAccumulator accumulator(state);
  hy::HydroSourceContext source_context{.update = {.dt_code = 0.1, .scale_factor = 1.0}};
  source_context.face_mass_flux_sink = &accumulator;
  accumulator.beginAdvance(geometry, dense, {});
  const hy::HydroCoreSolver solver(k_gamma);
  const hy::PiecewiseConstantReconstruction pcm;
  hy::HllcRiemannSolver riemann;
  hy::HydroScratchBuffers scratch;
  solver.advancePatchActiveSetWithScratch(
      conserved, geometry, hy::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
      source_context.update, pcm, riemann, {}, source_context, scratch, nullptr, nullptr);
  accumulator.endAdvance();
  CHUI_REQUIRE(accumulator.recordCount() == k_nx);  // 3 interior + 1 periodic box face

  // Store hydro mass back (what the runtime commit does).
  for (std::uint32_t r = 0; r < k_nx; ++r) {
    state.cells.mass_code[r] = conserved.loadCell(r).mass_density_comoving * geometry.cell_volume_comoving;
  }
  accumulator.prepare(state);
  const auto report = accumulator.commit(state);
  CHUI_REQUIRE(report.max_relative_mismatch < 1.0e-12);

  const double nu = 0.5 * 0.1 / 0.25;
  double o_total_after = 0.0;
  for (std::uint32_t r = 0; r < k_nx; ++r) {
    const std::uint32_t up = (r + k_nx - 1U) % k_nx;
    // dm_in = nu * M_up, dm_out = nu * M_r; donor fractions of step start.
    const double expected = o_before[r] - nu * mass_before[r] * (o_before[r] / mass_before[r]) +
        nu * mass_before[up] * (o_before[up] / mass_before[up]);
    CHUI_REQUIRE(std::abs(state.gas_composition.value(o_lane, r) - expected) < 1.0e-14);
    o_total_after += state.gas_composition.value(o_lane, r);
    CHUI_REQUIRE(std::abs(ph::compositionGasMass(state.gas_composition, r) - state.cells.mass_code[r]) <
                 1.0e-13 * state.cells.mass_code[r]);
  }
  CHUI_REQUIRE(std::abs(o_total_after - o_total_before) < 1.0e-16);
}

// Hand-built 2-cell geometry with one ghost for the boundary cases.
hy::HydroPatchGeometry makeTwoCellGeometry(hy::HydroBoundaryKind ghost_kind,
                                           hy::HydroFluxRegisterFaceRole role) {
  hy::HydroPatchGeometry g;
  g.nx = 2; g.ny = 1; g.nz = 1;
  g.cell_volume_comoving = 1.0;
  g.faces.push_back(hy::HydroFace{.owner_cell = 0, .neighbor_cell = 1, .area_comoving = 1.0});
  g.faces.push_back(hy::HydroFace{.owner_cell = 1, .neighbor_cell = 2, .ghost_cell_slot = 0, .area_comoving = 1.0});
  g.ghost_cells.push_back(hy::HydroGhostCell{.owner_real_cell = 1, .source_real_cell = 1, .ghost_cell = 2,
                                            .ghost_slot = 0, .boundary_kind = ghost_kind});
  g.flux_register_faces.resize(2);
  g.flux_register_faces[1].role = role;
  g.flux_register_faces[1].coarse_gas_cell_id = 502U;
  return g;
}

void testImportedGhostRefuses() {
  auto state = makeState();
  ph::enableGasComposition(state, ph::CompositionLayout::make({}), {});
  const auto g = makeTwoCellGeometry(hy::HydroBoundaryKind::kImportedMpi, hy::HydroFluxRegisterFaceRole::kNone);
  std::vector<std::uint32_t> dense{0U, 1U};
  ph::CompositionFluxAccumulator accumulator(state);
  accumulator.beginAdvance(g, dense, {});
  bool threw = false;
  try { accumulator.recordFaceMassTransfer(1, 1, 2, 1.0e-3); } catch (const std::runtime_error&) { threw = true; }
  CHUI_REQUIRE(threw);
}

void testFineFaceMirrorsOntoCoarseCell() {
  auto state = makeState();
  const auto layout = ph::CompositionLayout::make({});
  ph::enableGasComposition(state, layout, {});
  // Fine patch: cells 0,1; ghost slot 0 sourced from coarse cell 502 (row 2).
  const auto g = makeTwoCellGeometry(hy::HydroBoundaryKind::kInterior, hy::HydroFluxRegisterFaceRole::kFine);
  std::vector<std::uint32_t> dense{0U, 1U};
  std::vector<std::uint64_t> origins{502U};
  const std::uint32_t o_lane = layout.elementLane(ph::ChemicalElement::kO);
  double total_before = 0.0;
  for (std::uint32_t r = 0; r < k_nx; ++r) total_before += state.gas_composition.value(o_lane, r);
  const double x1 = state.gas_composition.value(o_lane, 1) / ph::compositionGasMass(state.gas_composition, 1);
  const double o2 = state.gas_composition.value(o_lane, 2);
  ph::CompositionFluxAccumulator accumulator(state);
  accumulator.beginAdvance(g, dense, origins);
  const double dm = 1.0e-4;  // fine cell 1 -> coarse cell 2
  accumulator.recordFaceMassTransfer(1, 1, hy::k_invalid_cell_index, dm);
  accumulator.endAdvance();
  // Coarse side on the coarse patch: role kCoarse must not update.
  auto coarse = makeTwoCellGeometry(hy::HydroBoundaryKind::kInterior, hy::HydroFluxRegisterFaceRole::kCoarse);
  std::vector<std::uint32_t> coarse_dense{2U, 3U};
  const std::vector<std::uint64_t> coarse_origins{501U};
  accumulator.beginAdvance(coarse, coarse_dense, coarse_origins);
  accumulator.recordFaceMassTransfer(1, 1, hy::k_invalid_cell_index, -2.0 * dm);
  accumulator.endAdvance();
  CHUI_REQUIRE(accumulator.recordCount() == 1U);
  state.cells.mass_code[1] -= dm;
  state.cells.mass_code[2] += dm;
  accumulator.prepare(state);
  (void)accumulator.commit(state);
  CHUI_REQUIRE(std::abs(state.gas_composition.value(o_lane, 2) - (o2 + dm * x1)) < 1.0e-16);
  double total_after = 0.0;
  for (std::uint32_t r = 0; r < k_nx; ++r) total_after += state.gas_composition.value(o_lane, r);
  CHUI_REQUIRE(std::abs(total_after - total_before) < 1.0e-17);
}

void testRecordOrderIndependence() {
  auto a = makeState();
  auto b = makeState();
  const auto layout = ph::CompositionLayout::make({});
  ph::enableGasComposition(a, layout, {});
  ph::enableGasComposition(b, layout, {});
  hy::HydroPatchGeometry g;
  g.nx = 4; g.ny = 1; g.nz = 1; g.cell_volume_comoving = 1.0;
  for (std::size_t c = 0; c + 1 < 4; ++c) {
    g.faces.push_back(hy::HydroFace{.owner_cell = c, .neighbor_cell = c + 1, .area_comoving = 1.0});
  }
  std::vector<std::uint32_t> dense{0U, 1U, 2U, 3U};
  ph::CompositionFluxAccumulator fa(a), fb(b);
  fa.beginAdvance(g, dense, {});
  fb.beginAdvance(g, dense, {});
  const double transfers[3] = {1.0e-3, -2.0e-3, 5.0e-4};
  for (std::size_t f = 0; f < 3; ++f) fa.recordFaceMassTransfer(f, f, f + 1, transfers[f]);
  for (std::size_t f = 3; f-- > 0;) fb.recordFaceMassTransfer(f, f, f + 1, transfers[f]);
  fa.endAdvance(); fb.endAdvance();
  for (auto* s : {&a, &b}) {
    s->cells.mass_code[0] -= transfers[0];
    s->cells.mass_code[1] += transfers[0] - transfers[1];
    s->cells.mass_code[2] += transfers[1] - transfers[2];
    s->cells.mass_code[3] += transfers[2];
  }
  fa.prepare(a); (void)fa.commit(a);
  fb.prepare(b); (void)fb.commit(b);
  for (std::size_t l = 0; l < a.gas_composition.laneCount(); ++l) {
    for (std::uint32_t r = 0; r < 4; ++r) {
      CHUI_REQUIRE(a.gas_composition.value(l, r) == b.gas_composition.value(l, r));
    }
  }
}

}  // namespace

int main() {
  testPeriodicDonorCellCompositionAdvection();
  testImportedGhostRefuses();
  testFineFaceMirrorsOntoCoarseCell();
  testRecordOrderIndependence();
  return 0;
}
