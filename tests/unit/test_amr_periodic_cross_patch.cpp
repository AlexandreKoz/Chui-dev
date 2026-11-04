// CHUI Part 3 Checkpoint A — periodic cross-patch / cross-rank AMR topology.
// STATUS: written as source only. NOT RUN — USER-RESERVED (the authoring agent
// executed no configure, build or test command).
//
// Physical invariants exercised (all exact or roundoff-level, no tuned
// tolerances):
//  1. Periodic constant state: a uniform (rho, u, P) state tiled by 2, 3 or 4
//     root patches across the periodic x box remains uniform after one step
//     with PCM, PLM and PPM (any spurious wall, missing ghost or wrong wrap
//     perturbs the cells next to the box face). Tolerance 1e-12 relative.
//  2. Translated periodic contact (uniform u, P; varying rho and Z) with PCM:
//     HLLC is exact for a contact, so one step is the donor-cell update
//     rho_i' = rho_i - nu (rho_i - rho_{i-1}) with i-1 wrapping across the box
//     face INTO ANOTHER PATCH; total mass and metal mass are conserved to
//     roundoff (same-level patch faces are evaluated from both sides with the
//     identical Riemann problem; see limitations in the Part 3 handoff).
//  3. Decomposition invariance: the 2-patch and 4-patch tilings of the same
//     8-cell state give identical cell values (tolerance 1e-13 relative).
//  4. Deep ghost layers: with ghost_layers = 3 every same-level boundary face
//     has neighbor_plus and neighbor_plus2 stencil cells (full PPM order).
//  5. Fail closed: a hierarchy leaving a periodic box face without any
//     wrapped neighbour patch is refused (std::runtime_error), never closed.
//  6. Distributed planner: two patches on different ranks that meet only
//     across the periodic wrap produce a boundary request only when the
//     periodic topology is supplied; the requested depth scales with the
//     ghost-layer count (clamped to the source patch dimension).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/parallel/distributed_mesh.hpp"
#include "../support/test_require.hpp"

namespace {

namespace amr = cosmosim::amr;
namespace core = cosmosim::core;
namespace hy = cosmosim::hydro;

constexpr double k_gamma = 5.0 / 3.0;
constexpr std::uint32_t k_cells = 8U;  // total cells along periodic x

struct CellInit {
  double rho = 1.0;
  double u = 0.0;
  double pressure = 1.0;
  double metal_fraction = 0.0;
};

// patch_count root patches of k_cells / patch_count cells each, tiling x in
// [0, 1); y, z are single-cell spanning axes (in-patch wrap).
core::SimulationState makeTiledState(std::uint32_t patch_count, const std::vector<CellInit>& cells) {
  const std::uint32_t per_patch = k_cells / patch_count;
  const double dx = 1.0 / static_cast<double>(k_cells);
  core::SimulationState state;
  state.resizeCells(k_cells);
  state.resizePatches(patch_count);
  std::vector<core::GasCellIdentityRecord> records;
  for (std::uint32_t p = 0; p < patch_count; ++p) {
    const amr::PatchDescriptor descriptor{
        .patch_id = 100U + p,
        .parent_patch_id = 0U,
        .level = 0U,
        .morton_key = 100U + p,
        .origin_comov = {static_cast<double>(p * per_patch) * dx, 0.0, 0.0},
        .extent_comov = {static_cast<double>(per_patch) * dx, 1.0, 1.0},
        .cell_dims = {static_cast<std::uint16_t>(per_patch), 1U, 1U}};
    amr::writePatchDescriptorToStateRow(state, p, descriptor);
    state.patches.first_cell[p] = p * per_patch;
    state.patches.cell_count[p] = per_patch;
    state.patches.owning_rank[p] = 0;
  }
  for (std::uint32_t row = 0; row < k_cells; ++row) {
    const std::uint32_t patch = row / per_patch;
    const CellInit& init = cells.at(row);
    state.cells.center_x_comoving[row] = (static_cast<double>(row) + 0.5) * dx;
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.patch_index[row] = patch;
    state.cells.time_bin[row] = 0;
    state.cells.mass_code[row] = init.rho * dx;
    state.gas_cells.gas_cell_id[row] = 7000U + row;  // identity independent of tiling
    state.gas_cells.parent_particle_id[row] = 0;
    state.gas_cells.density_code[row] = init.rho;
    state.gas_cells.pressure_code[row] = init.pressure;
    state.gas_cells.internal_energy_code[row] = init.pressure / ((k_gamma - 1.0) * init.rho);
    state.gas_cells.velocity_x_peculiar[row] = init.u;
    state.gas_cells.velocity_y_peculiar[row] = 0.0;
    state.gas_cells.velocity_z_peculiar[row] = 0.0;
    state.gas_cells.metal_mass_code[row] = init.metal_fraction * state.cells.mass_code[row];
    state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
    state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * init.pressure / init.rho);
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = 7000U + row, .parent_particle_id = std::nullopt,
        .owning_patch_id = 100U + patch, .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

amr::ProductionAmrHydroOptions periodicOptions(hy::HydroReconstructionScheme scheme) {
  amr::ProductionAmrHydroOptions options;
  options.physical_boundary_kind = hy::HydroBoundaryKind::kPeriodic;
  options.adiabatic_index = k_gamma;
  options.density_floor = 1.0e-12;
  options.pressure_floor = 1.0e-12;
  options.reconstruction_scheme = scheme;
  options.periodic_box_extent_comoving = {1.0, 1.0, 1.0};
  return options;
}

void advance(core::SimulationState& state, const amr::ProductionAmrHydroOptions& options, double dt) {
  const hy::HydroCoreSolver solver(k_gamma);
  hy::HllcRiemannSolver riemann;
  const hy::HydroUpdateContext update{.dt_code = dt, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  const hy::HydroSourceContext source_context{.update = update};
  (void)amr::advanceProductionAmrHydro(state, {}, update, source_context, solver, riemann, {}, options);
}

double relativeDifference(double a, double b) {
  return std::abs(a - b) / std::max({1.0, std::abs(a), std::abs(b)});
}

void testPeriodicConstantStateUnderTilings() {
  const std::vector<CellInit> uniform(k_cells, CellInit{.rho = 1.3, .u = 0.4, .pressure = 0.9, .metal_fraction = 0.02});
  for (const std::uint32_t patches : {2U, 4U, 8U}) {
    for (const auto scheme : {hy::HydroReconstructionScheme::kPiecewiseConstant,
                              hy::HydroReconstructionScheme::kPiecewiseLinear,
                              hy::HydroReconstructionScheme::kPiecewiseParabolic}) {
      if (scheme == hy::HydroReconstructionScheme::kPiecewiseParabolic && patches == 8U) {
        continue;  // single-cell patches: PPM deep layers cross two patches (covered by 4)
      }
      core::SimulationState state = makeTiledState(patches, uniform);
      advance(state, periodicOptions(scheme), 0.01);
      for (std::uint32_t row = 0; row < k_cells; ++row) {
        CHUI_REQUIRE(relativeDifference(state.gas_cells.density_code[row], 1.3) < 1.0e-12);
        CHUI_REQUIRE(relativeDifference(state.gas_cells.velocity_x_peculiar[row], 0.4) < 1.0e-12);
        CHUI_REQUIRE(relativeDifference(state.gas_cells.metal_mass_code[row] / state.cells.mass_code[row], 0.02) < 1.0e-12);
      }
    }
  }
}

std::vector<CellInit> contactCells() {
  std::vector<CellInit> cells;
  for (std::uint32_t row = 0; row < k_cells; ++row) {
    cells.push_back(CellInit{.rho = 1.0 + 0.1 * static_cast<double>(row), .u = 0.5, .pressure = 1.0,
                             .metal_fraction = 0.01 + 0.002 * static_cast<double>(row)});
  }
  return cells;
}

void testTranslatedContactCrossesPatchWrap() {
  const std::vector<CellInit> cells = contactCells();
  const double dt = 0.02;
  const double dx = 1.0 / static_cast<double>(k_cells);
  const double nu = 0.5 * dt / dx;
  std::vector<std::vector<double>> results;
  for (const std::uint32_t patches : {2U, 4U}) {
    core::SimulationState state = makeTiledState(patches, cells);
    double mass_before = 0.0;
    double metal_before = 0.0;
    for (std::uint32_t row = 0; row < k_cells; ++row) {
      mass_before += state.gas_cells.density_code[row] * dx;
      metal_before += state.gas_cells.metal_mass_code[row];
    }
    advance(state, periodicOptions(hy::HydroReconstructionScheme::kPiecewiseConstant), dt);
    double mass_after = 0.0;
    double metal_after = 0.0;
    std::vector<double> rho_after;
    for (std::uint32_t row = 0; row < k_cells; ++row) {
      mass_after += state.gas_cells.density_code[row] * dx;
      metal_after += state.gas_cells.metal_mass_code[row];
      rho_after.push_back(state.gas_cells.density_code[row]);
      const std::uint32_t upwind = (row + k_cells - 1U) % k_cells;  // wraps into the last patch
      const double expected = cells[row].rho - nu * (cells[row].rho - cells[upwind].rho);
      CHUI_REQUIRE(relativeDifference(state.gas_cells.density_code[row], expected) < 1.0e-12);
    }
    CHUI_REQUIRE(relativeDifference(mass_after, mass_before) < 1.0e-13);
    CHUI_REQUIRE(relativeDifference(metal_after, metal_before) < 1.0e-12);
    results.push_back(rho_after);
  }
  for (std::uint32_t row = 0; row < k_cells; ++row) {
    CHUI_REQUIRE(relativeDifference(results[0][row], results[1][row]) < 1.0e-13);
  }
}

void testDeepGhostLayersGiveFullPpmStencils() {
  const core::SimulationState state = makeTiledState(2U, contactCells());
  const auto descriptors = amr::buildProductionAmrPatchDescriptors(state);
  amr::AmrHydroGeometryOptions options;
  options.physical_boundary_kind = hy::HydroBoundaryKind::kPeriodic;
  options.ghost_layers = 3U;
  // x faces are patch interfaces (one direct, one across the wrap); y/z span.
  options.boundary_classes = {amr::AmrHydroBoundaryClass::kSameLevel, amr::AmrHydroBoundaryClass::kSameLevel,
                              amr::AmrHydroBoundaryClass::kPhysical, amr::AmrHydroBoundaryClass::kPhysical,
                              amr::AmrHydroBoundaryClass::kPhysical, amr::AmrHydroBoundaryClass::kPhysical};
  const auto geometry = amr::buildAmrHydroPatchGeometry(state, descriptors.front(), options);
  std::size_t deep = 0U;
  for (const auto& ghost : geometry.ghosts) {
    if (ghost.ghost_layer > 1U) {
      ++deep;
      CHUI_REQUIRE(ghost.boundary_class == amr::AmrHydroBoundaryClass::kSameLevel);
    }
  }
  CHUI_REQUIRE(deep == 4U);  // 2 x-faces * (layers 2 and 3); y/z periodic in-patch need none
  for (const hy::HydroFace& face : geometry.geometry.faces) {
    if (face.axis != hy::HydroFaceAxis::kX) continue;
    CHUI_REQUIRE(face.owner_minus_cell != hy::k_invalid_cell_index);
    CHUI_REQUIRE(face.neighbor_plus_cell != hy::k_invalid_cell_index);
    CHUI_REQUIRE(face.owner_minus2_cell != hy::k_invalid_cell_index);
    CHUI_REQUIRE(face.neighbor_plus2_cell != hy::k_invalid_cell_index);
  }
  CHUI_REQUIRE(geometry.ownedCapacityBytes() <=
               amr::amrHydroGeometryCapacity(descriptors.front(), 3U).retained_bytes);
}

void testUntiledPeriodicBoxFailsClosed() {
  std::vector<CellInit> cells(k_cells);
  core::SimulationState state = makeTiledState(2U, cells);
  // Drop the upper patch from the hierarchy: the box face of patch 100 at
  // x = 0 now has no wrapped neighbour.
  state.resizeCells(k_cells / 2U);
  state.resizePatches(1U);
  std::vector<core::GasCellIdentityRecord> records;
  for (std::uint32_t row = 0; row < k_cells / 2U; ++row) {
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = 7000U + row, .parent_particle_id = std::nullopt,
        .owning_patch_id = 100U, .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
  bool refused = false;
  try {
    advance(state, periodicOptions(hy::HydroReconstructionScheme::kPiecewiseLinear), 0.01);
  } catch (const std::runtime_error&) {
    refused = true;
  }
  CHUI_REQUIRE(refused);
}

void testDistributedPlannerSeesWrappedPeer() {
  namespace parallel = cosmosim::parallel;
  const parallel::AmrPatchPayloadRecord lower{
      .patch_id = 11U, .owner_rank = 0, .cell_count = 4U,
      .origin_x_comoving = 0.0, .origin_y_comoving = 0.0, .origin_z_comoving = 0.0,
      .extent_x_comoving = 0.25, .extent_y_comoving = 1.0, .extent_z_comoving = 1.0,
      .cell_dim_x = 4U, .cell_dim_y = 1U, .cell_dim_z = 1U};
  const parallel::AmrPatchPayloadRecord upper{
      .patch_id = 22U, .owner_rank = 1, .cell_count = 4U,
      .origin_x_comoving = 0.75, .origin_y_comoving = 0.0, .origin_z_comoving = 0.0,
      .extent_x_comoving = 0.25, .extent_y_comoving = 1.0, .extent_z_comoving = 1.0,
      .cell_dim_x = 4U, .cell_dim_y = 1U, .cell_dim_z = 1U};
  const auto none = parallel::planDirectedAmrPatchBoundaryCellRequests(
      std::span<const parallel::AmrPatchPayloadRecord>(&lower, 1U),
      std::span<const parallel::AmrPatchPayloadRecord>(&upper, 1U));
  CHUI_REQUIRE(none.empty());  // historical non-periodic plan unchanged
  const parallel::AmrInterfaceTopologyOptions topology{
      .periodic_box_extent_comoving = {1.0, 1.0, 1.0}, .ghost_layers = 3U};
  const auto wrapped = parallel::planDirectedAmrPatchBoundaryCellRequests(
      std::span<const parallel::AmrPatchPayloadRecord>(&lower, 1U),
      std::span<const parallel::AmrPatchPayloadRecord>(&upper, 1U), topology);
  CHUI_REQUIRE(wrapped.size() == 1U);
  CHUI_REQUIRE((wrapped[0].boundary_face_mask &
                static_cast<std::uint8_t>(parallel::AmrPatchBoundaryFace::kXLower)) != 0U);
  CHUI_REQUIRE(wrapped[0].boundary_face_depths[0] == 3U);  // three same-level ghost layers
  const auto reverse = parallel::planDirectedAmrPatchBoundaryCellRequests(
      std::span<const parallel::AmrPatchPayloadRecord>(&upper, 1U),
      std::span<const parallel::AmrPatchPayloadRecord>(&lower, 1U), topology);
  CHUI_REQUIRE(reverse.size() == 1U);
  CHUI_REQUIRE((reverse[0].boundary_face_mask &
                static_cast<std::uint8_t>(parallel::AmrPatchBoundaryFace::kXUpper)) != 0U);
}

}  // namespace

int main() {
  testPeriodicConstantStateUnderTilings();
  testTranslatedContactCrossesPatchWrap();
  testDeepGhostLayersGiveFullPpmStencils();
  testUntiledPeriodicBoxFailsClosed();
  testDistributedPlannerSeesWrappedPeer();
  return 0;
}
