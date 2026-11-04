// G2.2 regression (NOT RUN — USER-RESERVED): octet-family co-location of
// the runtime migration plan on 2 and 4 ranks.
//  * a level-1 family split across ranks is consolidated onto the owner of
//    its minimum-id child; every rank reports the same global change count;
//  * a family already co-located is untouched (no update emitted);
//  * a load-balance update that moves the minimum-id child carries the whole
//    family with it (family owner = POST-plan owner of the min-id child);
//  * level-0 patches are never touched by the rule.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "workflows/internal/amr_family_colocation.hpp"

namespace {

namespace core = cosmosim::core;
namespace par = cosmosim::parallel;

int g_failures = 0;
void require(bool ok, const std::string& what, int rank) {
  if (!ok) {
    ++g_failures;
    std::cerr << "[rank " << rank << "] FAILED: " << what << "\n";
  }
}

// Family A: parent 1000, children 10..17 dealt round-robin over ranks.
// Family B: parent 2000, children 20..27 all on the last rank.
// One level-0 patch per rank (id 900 + rank).
core::SimulationState makeState(int rank, int size) {
  std::vector<std::uint64_t> ids;
  std::vector<std::uint64_t> parents;
  std::vector<int> levels;
  for (std::uint64_t c = 0; c < 8U; ++c) {
    if (static_cast<int>(c % static_cast<std::uint64_t>(size)) == rank) {
      ids.push_back(10U + c); parents.push_back(1000U); levels.push_back(1);
    }
    if (rank == size - 1) { ids.push_back(20U + c); parents.push_back(2000U); levels.push_back(1); }
  }
  ids.push_back(900U + static_cast<std::uint64_t>(rank)); parents.push_back(0U); levels.push_back(0);
  core::SimulationState state;
  state.resizePatches(ids.size());
  for (std::size_t p = 0; p < ids.size(); ++p) {
    state.patches.patch_id[p] = ids[p];
    state.patches.parent_patch_id[p] = parents[p];
    state.patches.level[p] = levels[p];
    state.patches.owning_rank[p] = static_cast<std::uint32_t>(rank);
    state.patches.cell_count[p] = 0U;
  }
  return state;
}

int newOwner(const par::RuntimeRebalancePlan& plan, std::uint64_t id, int fallback) {
  for (const auto& u : plan.amr_patch_ownership_updates) if (u.patch_id == id) return u.new_owner_rank;
  return fallback;
}

}  // namespace

int main(int argc, char** argv) {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(&argc, &argv);
  {
    const par::MpiContext mpi{};
    const int rank = mpi.worldRank();
    const int size = mpi.worldSize();
    if (size >= 2) {
      // 1. Pure consolidation onto the owner of child 10 (rank 0).
      {
        const core::SimulationState state = makeState(rank, size);
        par::RuntimeRebalancePlan plan;
        const std::uint64_t changes =
            cosmosim::workflows::internal::enforceAmrFamilyColocation(state, mpi, nullptr, rank, plan);
        // Family A: 8 children, those not on rank 0 move (8 - ceil(8/size)).
        const std::uint64_t expected = 8U - (8U + static_cast<std::uint64_t>(size) - 1U) / static_cast<std::uint64_t>(size);
        require(changes == expected, "global change count of the split family", rank);
        require(plan.should_rebalance, "consolidation forces a migration", rank);
        for (std::size_t p = 0; p < state.patches.size(); ++p) {
          const std::uint64_t id = state.patches.patch_id[p];
          const int owner = newOwner(plan, id, rank);
          if (state.patches.parent_patch_id[p] == 1000U) require(owner == 0, "family A on rank 0", rank);
          if (state.patches.parent_patch_id[p] == 2000U) require(owner == size - 1, "family B untouched", rank);
          if (state.patches.level[p] == 0) require(owner == rank, "level-0 patch untouched", rank);
        }
      }
      // 2. The load balancer moves child 10 (on rank 0) to the last rank: the
      //    whole family follows it.
      {
        const core::SimulationState state = makeState(rank, size);
        par::RuntimeRebalancePlan plan;
        if (rank == 0) {
          plan.amr_patch_ownership_updates.push_back(
              par::AmrPatchOwnershipUpdate{.patch_id = 10U, .old_owner_rank = 0, .new_owner_rank = size - 1});
        }
        (void)cosmosim::workflows::internal::enforceAmrFamilyColocation(state, mpi, nullptr, rank, plan);
        for (std::size_t p = 0; p < state.patches.size(); ++p) {
          if (state.patches.parent_patch_id[p] != 1000U) continue;
          require(newOwner(plan, state.patches.patch_id[p], rank) == size - 1,
                  "family A follows its minimum-id child", rank);
        }
        if (rank == size - 1) {
          require(plan.amr_patch_ownership_updates.empty(), "children already on the family owner stay put", rank);
        }
      }
      // 3. Already co-located everywhere: no change, no forced migration.
      {
        core::SimulationState state = makeState(rank, size);
        for (std::size_t p = 0; p < state.patches.size(); ++p) {
          if (state.patches.parent_patch_id[p] == 1000U) state.patches.parent_patch_id[p] = 3000U + static_cast<std::uint64_t>(rank);
        }
        par::RuntimeRebalancePlan plan;
        const std::uint64_t changes =
            cosmosim::workflows::internal::enforceAmrFamilyColocation(state, mpi, nullptr, rank, plan);
        require(changes == 0U && !plan.should_rebalance, "co-located families need no migration", rank);
      }
    }
    g_failures = static_cast<int>(mpi.allreduceSumUint64(static_cast<std::uint64_t>(g_failures)));
  }
  MPI_Finalize();
  return g_failures == 0 ? 0 : 1;
#else
  (void)argc;
  (void)argv;
  return 77;
#endif
}
