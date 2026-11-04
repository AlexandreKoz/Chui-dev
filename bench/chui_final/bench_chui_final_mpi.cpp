// CHUI final-design MPI benchmarks (V5). SOURCE ONLY — NOT RUN — USER-RESERVED.
// Usage: mpiexec -n <ranks> bench_chui_final_mpi [cells_per_rank=4096] [repeats=5]
// A periodic 1D line of single-cell patches, contiguous over ranks. Metrics
// (rank-0 median of the per-repeat MAX over ranks):
//   distributed_diffusion_gather_seconds: directed remote-cell exchange;
//   distributed_diffusion_setup_seconds: claim/conductance handshake;
//   distributed_diffusion_advance_seconds: full element-resolved operator
//     (subcycles x 2 exchanges); distributed_diffusion_subcycles;
//   family_colocation_seconds: octet co-location census + all-gather.
// No result is claimed; numbers come only from the owner's runs.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <vector>

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/element_diffusion.hpp"
#include "workflows/internal/amr_family_colocation.hpp"
#include "workflows/internal/distributed_element_diffusion.hpp"

namespace {
namespace core = cosmosim::core;
namespace ph = cosmosim::physics;
namespace wi = cosmosim::workflows::internal;
using Clock = std::chrono::steady_clock;

core::SimulationState makeLine(int rank, std::uint32_t per_rank) {
  core::SimulationState state;
  state.resizeCells(per_rank);
  state.resizePatches(per_rank);
  std::vector<core::GasCellIdentityRecord> ids;
  for (std::uint32_t r = 0; r < per_rank; ++r) {
    const std::uint64_t g = static_cast<std::uint64_t>(rank) * per_rank + r;
    state.patches.patch_id[r] = 1U + g;
    state.patches.level[r] = 0;
    state.patches.morton_key[r] = 1U + g;
    state.patches.first_cell[r] = r;
    state.patches.cell_count[r] = 1U;
    state.patches.owning_rank[r] = static_cast<std::uint32_t>(rank);
    state.patches.cell_dim_x[r] = state.patches.cell_dim_y[r] = state.patches.cell_dim_z[r] = 1U;
    state.patches.origin_x_comoving[r] = static_cast<double>(g);
    state.patches.extent_x_comoving[r] = state.patches.extent_y_comoving[r] = state.patches.extent_z_comoving[r] = 1.0;
    state.cells.center_x_comoving[r] = static_cast<double>(g) + 0.5;
    state.cells.center_y_comoving[r] = state.cells.center_z_comoving[r] = 0.5;
    state.cells.mass_code[r] = 1.0;
    state.gas_cells.gas_cell_id[r] = 1U + g;
    state.gas_cells.density_code[r] = 1.0;
    state.gas_cells.metal_mass_code[r] = (g % 2U == 0U) ? 0.02 : 0.001;
    ids.push_back({.gas_cell_id = 1U + g, .parent_particle_id = std::nullopt, .owning_patch_id = 1U + g,
                   .local_cell_row = r});
  }
  state.gas_cell_identity.assign(std::move(ids));
  ph::enableGasComposition(state, ph::CompositionLayout::make({.network = ph::CompositionNetworkKind::kHybridLight,
      .dust = ph::CompositionDustKind::kTwoReservoirV1, .sgs_energy = false}), {});
  return state;
}
}  // namespace

int main(int argc, char** argv) {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(&argc, &argv);
  {
    const cosmosim::parallel::MpiContext mpi{};
    const std::uint32_t per_rank = argc > 1 ? static_cast<std::uint32_t>(std::strtoul(argv[1], nullptr, 10)) : 4096U;
    const int repeats = argc > 2 ? std::atoi(argv[2]) : 5;
    const double box = static_cast<double>(per_rank) * mpi.worldSize();
    std::vector<double> gather_t, setup_t, advance_t, family_t;
    std::uint32_t subcycles = 0U;
    for (int rep = 0; rep <= repeats; ++rep) {
      core::SimulationState state = makeLine(mpi.worldRank(), per_rank);
      std::vector<std::uint8_t> mask(state.cells.size(), 1U);
      auto t0 = Clock::now();
      const auto gather = wi::gatherRemoteDiffusionCells(state, mpi, nullptr, 1.0, true, {box, 1.0, 1.0}, 11U);
      const auto topology = wi::buildMetalDiffusionTopology(state, mask, static_cast<std::uint32_t>(mpi.worldRank()),
          1.0, core::BoundaryCondition::kPeriodic, gather.cells, gather.global_envelope_physical);
      const double tg = std::chrono::duration<double>(Clock::now() - t0).count();
      std::vector<double> rho_kappa(state.cells.size(), 0.05);
      t0 = Clock::now();
      wi::DistributedElementDiffusionCoupling coupling(mpi, nullptr, state, topology.remote_faces, rho_kappa,
                                                       state.gas_composition.laneCount());
      coupling.setup();
      const double ts = std::chrono::duration<double>(Clock::now() - t0).count();
      std::vector<double> workspace;
      t0 = Clock::now();
      const auto report = ph::advanceElementResolvedDiffusion(state.gas_composition, topology.faces, rho_kappa,
                                                             1.0, 0.4, 100000U, workspace, &coupling);
      const double ta = std::chrono::duration<double>(Clock::now() - t0).count();
      subcycles = report.subcycles;
      cosmosim::parallel::RuntimeRebalancePlan plan;
      t0 = Clock::now();
      (void)wi::enforceAmrFamilyColocation(state, mpi, nullptr, mpi.worldRank(), plan);
      const double tf = std::chrono::duration<double>(Clock::now() - t0).count();
      if (rep == 0) continue;  // warm-up
      gather_t.push_back(-mpi.allreduceMinDouble(-tg));
      setup_t.push_back(-mpi.allreduceMinDouble(-ts));
      advance_t.push_back(-mpi.allreduceMinDouble(-ta));
      family_t.push_back(-mpi.allreduceMinDouble(-tf));
    }
    const auto median = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2U]; };
    if (mpi.isRoot()) {
      std::printf("ranks %d\ncells_per_rank %u\n", mpi.worldSize(), per_rank);
      std::printf("distributed_diffusion_gather_seconds %.6e\n", median(gather_t));
      std::printf("distributed_diffusion_setup_seconds %.6e\n", median(setup_t));
      std::printf("distributed_diffusion_advance_seconds %.6e\n", median(advance_t));
      std::printf("distributed_diffusion_subcycles %u\n", subcycles);
      std::printf("family_colocation_seconds %.6e\n", median(family_t));
    }
  }
  MPI_Finalize();
  return 0;
#else
  (void)argc;
  (void)argv;
  std::printf("bench_chui_final_mpi requires COSMOSIM_ENABLE_MPI\n");
  return 77;
#endif
}
