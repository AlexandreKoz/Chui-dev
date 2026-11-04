// G2.1 regression (NOT RUN — USER-RESERVED): distributed element-resolved
// diffusion on 2 and 4 ranks.
//
//  * contiguous, interleaved and empty-rank decompositions of the SAME
//    8-cell line (open and x-periodic) reproduce the single-domain operator
//    on every owned row to 1e-12 relative (roundoff: the per-row summation
//    order of face contributions may differ between decompositions);
//  * every mixed lane total is conserved globally to 1e-13 relative;
//  * the unique flux authority evaluates every interface exactly once
//    (authority faces summed over ranks == interface count);
//  * a collective refusal (interval needs more than max_subcycles) throws on
//    EVERY rank without deadlock and leaves every sidecar bitwise unchanged.
//
// Acceptance thresholds: 1e-12 is ~1e4 ulp of O(1) mass fractions after a
// few hundred subcycles of an explicit convex update; it is not widened to
// accommodate a failure.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/element_diffusion.hpp"
#include "workflows/internal/distributed_element_diffusion.hpp"
#include "workflows/internal/metal_diffusion_topology.hpp"

namespace {

namespace core = cosmosim::core;
namespace ph = cosmosim::physics;
namespace wi = cosmosim::workflows::internal;

constexpr std::uint32_t k_global_cells = 8U;

int g_failures = 0;

void require(bool condition, const std::string& what, int rank) {
  if (!condition) {
    ++g_failures;
    std::cerr << "[rank " << rank << "] FAILED: " << what << "\n";
  }
}

enum class Decomposition { kContiguous, kInterleaved, kEmptyLastRank };

int ownerOf(std::uint32_t global_cell, int world_size, Decomposition decomposition) {
  switch (decomposition) {
    case Decomposition::kContiguous:
      return static_cast<int>(global_cell * static_cast<std::uint32_t>(world_size) / k_global_cells);
    case Decomposition::kInterleaved:
      return static_cast<int>(global_cell % static_cast<std::uint32_t>(world_size));
    case Decomposition::kEmptyLastRank: {
      const int used = std::max(world_size - 1, 1);
      return static_cast<int>(global_cell * static_cast<std::uint32_t>(used) / k_global_cells);
    }
  }
  return 0;
}

// Builds the rows of `global_cells` owned by `rank` (rank < 0: everything,
// owned by rank 0, the single-domain reference).
core::SimulationState makeState(int rank, int world_size, Decomposition decomposition) {
  std::vector<std::uint32_t> mine;
  for (std::uint32_t g = 0; g < k_global_cells; ++g) {
    if (rank < 0 || ownerOf(g, world_size, decomposition) == rank) mine.push_back(g);
  }
  const std::uint32_t owner_rank = rank < 0 ? 0U : static_cast<std::uint32_t>(rank);
  core::SimulationState state;
  state.resizeCells(mine.size());
  state.resizePatches(mine.size());
  std::vector<core::GasCellIdentityRecord> identities;
  for (std::uint32_t row = 0; row < mine.size(); ++row) {
    const std::uint32_t g = mine[row];
    state.patches.patch_id[row] = 100U + g;
    state.patches.parent_patch_id[row] = 0U;
    state.patches.level[row] = 0U;
    state.patches.morton_key[row] = 1U + g;
    state.patches.first_cell[row] = row;
    state.patches.cell_count[row] = 1U;
    state.patches.owning_rank[row] = owner_rank;
    state.patches.cell_dim_x[row] = 1U;
    state.patches.cell_dim_y[row] = 1U;
    state.patches.cell_dim_z[row] = 1U;
    state.patches.origin_x_comoving[row] = static_cast<double>(g);
    state.patches.origin_y_comoving[row] = 0.0;
    state.patches.origin_z_comoving[row] = 0.0;
    state.patches.extent_x_comoving[row] = 1.0;
    state.patches.extent_y_comoving[row] = 1.0;
    state.patches.extent_z_comoving[row] = 1.0;
    state.cells.center_x_comoving[row] = static_cast<double>(g) + 0.5;
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.mass_code[row] = 1.0;
    state.gas_cells.gas_cell_id[row] = 500U + g;
    state.gas_cells.density_code[row] = 1.0;
    state.gas_cells.internal_energy_code[row] = 1.0;
    state.gas_cells.metal_mass_code[row] = 0.001 + 0.004 * static_cast<double>(g);
    identities.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = 500U + g, .parent_particle_id = std::nullopt, .owning_patch_id = 100U + g,
        .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(identities));
  ph::enableGasComposition(state, ph::CompositionLayout::make({
      .network = ph::CompositionNetworkKind::kDetailedCarbonOxygen,
      .dust = ph::CompositionDustKind::kTwoReservoirV1, .sgs_energy = true}), {});
  return state;
}

struct RunResult {
  ph::ElementDiffusionReport report{};
  std::uint64_t authority_faces = 0U;
};

RunResult runDistributed(core::SimulationState& state, const cosmosim::parallel::MpiContext& mpi, bool periodic,
                         std::uint32_t max_subcycles, double dt) {
  const std::array<double, 3> box{static_cast<double>(k_global_cells), 1.0, 1.0};
  std::vector<std::uint8_t> mask(state.cells.size(), 1U);
  const auto gather = wi::gatherRemoteDiffusionCells(state, mpi, nullptr, 1.0, periodic, box, 7U);
  const auto topology = wi::buildMetalDiffusionTopology(
      state, mask, static_cast<std::uint32_t>(mpi.worldRank()), 1.0,
      periodic ? core::BoundaryCondition::kPeriodic : core::BoundaryCondition::kOpen, gather.cells,
      gather.global_envelope_physical);
  std::vector<double> rho_kappa(state.cells.size(), 0.05);
  wi::DistributedElementDiffusionCoupling coupling(mpi, nullptr, state, topology.remote_faces, rho_kappa,
                                                   state.gas_composition.laneCount());
  coupling.setup();
  std::vector<double> workspace;
  RunResult result;
  result.report = ph::advanceElementResolvedDiffusion(state.gas_composition, topology.faces, rho_kappa, dt, 0.4,
                                                     max_subcycles, workspace, &coupling);
  result.authority_faces = coupling.authorityFaceCount();
  return result;
}

void checkDecomposition(const cosmosim::parallel::MpiContext& mpi, Decomposition decomposition, bool periodic) {
  const int rank = mpi.worldRank();
  const int size = mpi.worldSize();
  const std::string label = std::string(periodic ? "periodic" : "open") + "/decomposition " +
      std::to_string(static_cast<int>(decomposition));
  // Single-domain reference on every rank.
  core::SimulationState reference = makeState(-1, size, decomposition);
  {
    std::vector<std::uint8_t> mask(reference.cells.size(), 1U);
    const auto topology = wi::buildMetalDiffusionTopology(
        reference, mask, 0U, 1.0, periodic ? core::BoundaryCondition::kPeriodic : core::BoundaryCondition::kOpen);
    std::vector<double> rho_kappa(reference.cells.size(), 0.05);
    std::vector<double> workspace;
    (void)ph::advanceElementResolvedDiffusion(reference.gas_composition, topology.faces, rho_kappa, 0.5, 0.4, 1000U,
                                              workspace);
  }
  core::SimulationState state = makeState(rank, size, decomposition);
  std::vector<double> before_totals(state.gas_composition.laneCount(), 0.0);
  for (std::size_t lane = 0; lane < before_totals.size(); ++lane) {
    for (std::size_t row = 0; row < state.cells.size(); ++row) before_totals[lane] += state.gas_composition.value(lane, row);
  }
  mpi.allreduceSumDoublesInPlace(before_totals);
  const RunResult run = runDistributed(state, mpi, periodic, 1000U, 0.5);
  require(run.report.subcycles > 0U, label + ": operator subcycled", rank);
  // Interfaces of the line: 7 (open) or 8 (periodic). Every interface that
  // crosses ranks has exactly one authority.
  std::uint64_t remote_interfaces = 0U;
  for (std::uint32_t g = 0; g < k_global_cells; ++g) {
    const std::uint32_t right = (g + 1U) % k_global_cells;
    if (!periodic && g + 1U == k_global_cells) continue;
    if (ownerOf(g, size, decomposition) != ownerOf(right, size, decomposition)) ++remote_interfaces;
  }
  const std::uint64_t authority_total = mpi.allreduceSumUint64(run.authority_faces);
  require(authority_total == remote_interfaces, label + ": one authority per remote interface", rank);
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const std::uint64_t id = state.gas_cells.gas_cell_id[row];
    const std::uint32_t ref_row = static_cast<std::uint32_t>(id - 500U);
    for (std::size_t lane = 0; lane < state.gas_composition.laneCount(); ++lane) {
      const double a = state.gas_composition.value(lane, row);
      const double b = reference.gas_composition.value(lane, ref_row);
      require(std::abs(a - b) <= 1.0e-12 * std::max({std::abs(a), std::abs(b), 1.0e-300}),
              label + ": lane " + std::to_string(lane) + " row id " + std::to_string(id) + " matches reference",
              rank);
    }
  }
  std::vector<double> after_totals(state.gas_composition.laneCount(), 0.0);
  for (std::size_t lane = 0; lane < after_totals.size(); ++lane) {
    for (std::size_t row = 0; row < state.cells.size(); ++row) after_totals[lane] += state.gas_composition.value(lane, row);
  }
  mpi.allreduceSumDoublesInPlace(after_totals);
  for (std::size_t lane = 0; lane < after_totals.size(); ++lane) {
    require(std::abs(after_totals[lane] - before_totals[lane]) <=
                1.0e-13 * std::max(std::abs(before_totals[lane]), 1.0e-300),
            label + ": global lane total conserved", rank);
  }
}

void checkCollectiveRefusal(const cosmosim::parallel::MpiContext& mpi) {
  const int rank = mpi.worldRank();
  core::SimulationState state = makeState(rank, mpi.worldSize(), Decomposition::kContiguous);
  const core::SimulationState before = state;
  bool refused = false;
  try {
    (void)runDistributed(state, mpi, false, 1U, 1.0e6);
  } catch (const std::runtime_error&) {
    refused = true;
  }
  const std::uint64_t refusals = mpi.allreduceSumUint64(refused ? 1U : 0U);
  require(refusals == static_cast<std::uint64_t>(mpi.worldSize()), "refusal is collective", rank);
  for (std::size_t lane = 0; lane < state.gas_composition.laneCount(); ++lane) {
    for (std::size_t row = 0; row < state.cells.size(); ++row) {
      require(state.gas_composition.value(lane, row) == before.gas_composition.value(lane, row),
              "refused interval leaves the sidecar bitwise unchanged", rank);
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(&argc, &argv);
  {
    const cosmosim::parallel::MpiContext mpi{};
    for (const bool periodic : {false, true}) {
      checkDecomposition(mpi, Decomposition::kContiguous, periodic);
      checkDecomposition(mpi, Decomposition::kInterleaved, periodic);
      if (mpi.worldSize() > 1) checkDecomposition(mpi, Decomposition::kEmptyLastRank, periodic);
    }
    checkCollectiveRefusal(mpi);
    const std::uint64_t failures = mpi.allreduceSumUint64(static_cast<std::uint64_t>(g_failures));
    g_failures = static_cast<int>(failures);
  }
  MPI_Finalize();
  return g_failures == 0 ? 0 : 1;
#else
  (void)argc;
  (void)argv;
  return 77;
#endif
}
