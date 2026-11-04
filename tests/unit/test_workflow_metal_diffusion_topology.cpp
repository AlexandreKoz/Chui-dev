#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <optional>
#include <array>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"
#include "workflows/internal/metal_diffusion_topology.hpp"

namespace {

void setPatch(
    cosmosim::core::SimulationState& state,
    std::uint32_t patch,
    std::uint32_t first,
    std::uint32_t count,
    std::uint32_t nx,
    std::uint32_t ny,
    std::uint32_t nz,
    double ox,
    double oy,
    double oz,
    double ex,
    double ey,
    double ez) {
  state.patches.patch_id[patch] = 100U + patch;
  state.patches.first_cell[patch] = first;
  state.patches.cell_count[patch] = count;
  state.patches.owning_rank[patch] = 0U;
  state.patches.cell_dim_x[patch] = nx;
  state.patches.cell_dim_y[patch] = ny;
  state.patches.cell_dim_z[patch] = nz;
  state.patches.origin_x_comoving[patch] = ox;
  state.patches.origin_y_comoving[patch] = oy;
  state.patches.origin_z_comoving[patch] = oz;
  state.patches.extent_x_comoving[patch] = ex;
  state.patches.extent_y_comoving[patch] = ey;
  state.patches.extent_z_comoving[patch] = ez;
}

std::vector<cosmosim::physics::MetalDiffusionCell> makeDiffusionCells(
    std::size_t count) {
  std::vector<cosmosim::physics::MetalDiffusionCell> cells(count);
  for (std::size_t row = 0; row < count; ++row) {
    cells[row].gas_cell_id = 1000U + row;
    cells[row].gas_mass_code = 1.0;
    cells[row].metal_mass_code = row == 0U ? 0.5 : 0.0;
    cells[row].density_code = 1.0;
    cells[row].volume_code = 1.0;
    cells[row].filter_length_code = 1.0;
    cells[row].is_owned_leaf = true;
  }
  return cells;
}

void testSameLevelPatchBoundaryTransfersConservatively() {
  cosmosim::core::SimulationState state;
  state.resizeCells(2U);
  state.resizePatches(2U);
  setPatch(state, 0U, 0U, 1U, 1U, 1U, 1U, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
  setPatch(state, 1U, 1U, 1U, 1U, 1U, 1U, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0);
  state.cells.center_x_comoving[0] = 0.5;
  state.cells.center_x_comoving[1] = 1.5;
  state.cells.center_y_comoving[0] = state.cells.center_y_comoving[1] = 0.5;
  state.cells.center_z_comoving[0] = state.cells.center_z_comoving[1] = 0.5;
  state.gas_cells.velocity_x_peculiar[0] = 0.0;
  state.gas_cells.velocity_x_peculiar[1] = 1.0;

  std::vector<std::uint8_t> owned_leaf{1U, 1U};
  auto cells = makeDiffusionCells(2U);
  const auto topology = cosmosim::workflows::internal::buildMetalDiffusionTopology(
      state,
      owned_leaf,
      0U,
      1.0,
      cosmosim::core::BoundaryCondition::kOpen);
  assert(topology.interior_patch_face_count == 0U);
  assert(topology.cross_patch_face_count == 1U);
  assert(topology.periodic_wrap_face_count == 0U);
  assert(topology.faces.size() == 1U);
  assert(std::abs(topology.faces[0].area_code - 1.0) < 1.0e-14);
  assert(std::abs(topology.faces[0].center_distance_code - 1.0) < 1.0e-14);

  cosmosim::physics::MetalDiffusionConfig config;
  config.enabled = true;
  config.model = cosmosim::core::MetalDiffusionModel::kSmagorinsky;
  config.smagorinsky_coefficient = 0.05;
  config.parabolic_cfl = 0.4;
  cosmosim::physics::MetalDiffusionModel model(config);
  std::vector<double> gas_mass{1.0, 1.0};
  std::vector<double> metal_mass{0.5, 0.0};
  std::vector<double> rho_kappa(2U, 0.0);
  for (std::size_t i = 0U; i < rho_kappa.size(); ++i) {
    rho_kappa[i] = cosmosim::physics::smagorinskyRhoKappaCode(
        config, 1.0, 1.0, topology.strain_magnitude_code[i]);
  }
  cosmosim::physics::MetalDiffusionWorkspace workspace;
  const double before = metal_mass[0] + metal_mass[1];
  const auto report = model.advanceFromView(
      cosmosim::physics::MetalDiffusionFieldView{
          .gas_mass_code = gas_mass,
          .metal_mass_code = metal_mass,
          .rho_kappa_code = rho_kappa,
      },
      topology.faces, 1.0e-2, workspace);
  assert(report.faces_evaluated > 0U);
  assert(metal_mass[0] < 0.5);
  assert(metal_mass[1] > 0.0);
  assert(std::abs(metal_mass[0] + metal_mass[1] - before) < 1.0e-14);
}

void testCoarseFineOverlapSplitsFaceArea() {
  cosmosim::core::SimulationState state;
  state.resizeCells(5U);
  state.resizePatches(2U);
  setPatch(state, 0U, 0U, 1U, 1U, 1U, 1U, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0);
  state.patches.level[0] = 0U;
  setPatch(state, 1U, 1U, 4U, 1U, 2U, 2U, 1.0, 0.0, 0.0, 1.0, 2.0, 2.0);
  state.patches.level[1] = 1U;
  state.cells.center_x_comoving[0] = 0.5;
  state.cells.center_y_comoving[0] = 1.0;
  state.cells.center_z_comoving[0] = 1.0;
  std::uint32_t row = 1U;
  for (std::uint32_t k = 0U; k < 2U; ++k) {
    for (std::uint32_t j = 0U; j < 2U; ++j) {
      state.cells.center_x_comoving[row] = 1.5;
      state.cells.center_y_comoving[row] = 0.5 + static_cast<double>(j);
      state.cells.center_z_comoving[row] = 0.5 + static_cast<double>(k);
      ++row;
    }
  }
  std::vector<std::uint8_t> owned_leaf(5U, 1U);
  auto cells = makeDiffusionCells(5U);
  const auto topology = cosmosim::workflows::internal::buildMetalDiffusionTopology(
      state,
      owned_leaf,
      0U,
      1.0,
      cosmosim::core::BoundaryCondition::kOpen);
  assert(topology.cross_patch_face_count == 4U);
  double coarse_interface_area = 0.0;
  std::uint64_t coarse_interface_faces = 0U;
  for (const auto& face : topology.faces) {
    if (face.left_cell == 0U || face.right_cell == 0U) {
      coarse_interface_area += face.area_code;
      ++coarse_interface_faces;
    }
  }
  assert(coarse_interface_faces == 4U);
  assert(std::abs(coarse_interface_area - 4.0) < 1.0e-14);
}

void testPeriodicWrapIsARealDiffusionInterface() {
  cosmosim::core::SimulationState state;
  state.resizeCells(2U);
  state.resizePatches(1U);
  setPatch(state, 0U, 0U, 2U, 2U, 1U, 1U, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0);
  state.cells.center_x_comoving[0] = 0.5;
  state.cells.center_x_comoving[1] = 1.5;
  state.cells.center_y_comoving[0] = state.cells.center_y_comoving[1] = 0.5;
  state.cells.center_z_comoving[0] = state.cells.center_z_comoving[1] = 0.5;
  state.gas_cells.velocity_x_peculiar[0] = -1.0;
  state.gas_cells.velocity_x_peculiar[1] = 1.0;
  std::vector<std::uint8_t> owned_leaf{1U, 1U};
  auto cells = makeDiffusionCells(2U);
  const auto topology = cosmosim::workflows::internal::buildMetalDiffusionTopology(
      state,
      owned_leaf,
      0U,
      1.0,
      cosmosim::core::BoundaryCondition::kPeriodic);
  assert(topology.interior_patch_face_count == 1U);
  assert(topology.periodic_wrap_face_count == 1U);
  assert(topology.faces.size() == 2U);
  assert(std::count_if(
             topology.faces.begin(), topology.faces.end(),
             [](const auto& face) {
               return face.boundary_kind ==
                   cosmosim::physics::MetalDiffusionBoundaryKind::kPeriodic;
             }) == 1);
}

// G2.1/G0.5: a remote neighbour cell yields exactly one remote face with the
// owner's geometry, contributes its velocity to the strain, never creates a
// local face, and the pre-construction bound covers the measured scratch.
void testRemoteNeighbourCreatesOneRemoteFace() {
  cosmosim::core::SimulationState state;
  state.resizeCells(1U);
  state.resizePatches(1U);
  setPatch(state, 0U, 0U, 1U, 1U, 1U, 1U, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
  state.cells.center_x_comoving[0] = 0.5;
  state.cells.center_y_comoving[0] = 0.5;
  state.cells.center_z_comoving[0] = 0.5;
  state.gas_cells.gas_cell_id[0] = 7U;
  state.gas_cells.velocity_x_peculiar[0] = 0.0;
  std::vector<std::uint8_t> owned_leaf{1U};
  cosmosim::workflows::internal::RemoteDiffusionCell remote;
  remote.gas_cell_id = 3U;
  remote.owner_rank = 1;
  remote.patch_origin_comoving = {1.0, 0.0, 0.0};
  remote.cell_width_comoving = {1.0, 1.0, 1.0};
  remote.index = {0U, 0U, 0U};
  remote.patch_dims = {1U, 1U, 1U};
  remote.velocity_peculiar = {2.0, 0.0, 0.0};
  const std::vector<cosmosim::workflows::internal::RemoteDiffusionCell> remotes{remote};
  const auto topology = cosmosim::workflows::internal::buildMetalDiffusionTopology(
      state, owned_leaf, 0U, 1.0, cosmosim::core::BoundaryCondition::kOpen, remotes);
  assert(topology.faces.empty());
  assert(topology.remote_faces.size() == 1U);
  const auto& face = topology.remote_faces[0];
  assert(face.local_row == 0U && face.remote_gas_cell_id == 3U && face.local_gas_cell_id == 7U);
  assert(face.axis == 0U && face.local_side == 1 && !face.periodic && face.remote_rank == 1);
  assert(std::abs(face.area_code - 1.0) < 1.0e-14 && std::abs(face.center_distance_code - 1.0) < 1.0e-14);
  // One-sided derivative du/dx = (2 - 0) / 1 contributes to the strain.
  assert(topology.strain_magnitude_code[0] > 0.0);
  const std::uint64_t bound = cosmosim::workflows::internal::metalDiffusionTopologyConstructionBoundBytes(
      state, owned_leaf, 0U, remotes.size());
  assert(bound >= topology.construction_scratch_high_water_bytes);
  // Periodic adjacency across a split axis is a distinct second face when the
  // global envelope says the box is [0, 2): same pair, other side.
  const std::optional<std::array<double, 6>> envelope = std::array<double, 6>{0.0, 0.0, 0.0, 2.0, 1.0, 1.0};
  const auto periodic = cosmosim::workflows::internal::buildMetalDiffusionTopology(
      state, owned_leaf, 0U, 1.0, cosmosim::core::BoundaryCondition::kPeriodic, remotes, envelope);
  assert(periodic.remote_faces.size() == 2U);
  bool saw_direct = false, saw_wrap = false;
  for (const auto& f : periodic.remote_faces) {
    if (f.periodic) { saw_wrap = true; assert(f.local_side == -1); }
    else { saw_direct = true; assert(f.local_side == 1); }
  }
  assert(saw_direct && saw_wrap);
}

}  // namespace


int main() {
  testSameLevelPatchBoundaryTransfersConservatively();
  testCoarseFineOverlapSplitsFaceArea();
  testPeriodicWrapIsARealDiffusionInterface();
  testRemoteNeighbourCreatesOneRemoteFace();
  return 0;
}
