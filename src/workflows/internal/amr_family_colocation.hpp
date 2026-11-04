#pragma once

// G2.2: octet-family co-location rule applied to every runtime migration
// plan (see migration_balance_runtime.cpp). Collective over the world.

#include <cstdint>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows::internal {

// Rewrites `plan` so that, after migration, every set of level > 0 leaf
// siblings (same parent_patch_id) lives on the post-plan owner of its
// minimum-id child. Returns the global number of changed assignments and
// sets plan.should_rebalance when it is non-zero.
[[nodiscard]] std::uint64_t enforceAmrFamilyColocation(
    const core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    core::MemoryGovernor* governor,
    int world_rank,
    parallel::RuntimeRebalancePlan& plan);

}  // namespace cosmosim::workflows::internal
