#include "cosmosim/workflows/source_runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <utility>
#include <unordered_map>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <memory_resource>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/adaptive_resource_control.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/governed_scratch_arena.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/retained_capacity_transaction.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/baryonic_runtime_config.hpp"
#include "cosmosim/physics/thermochemistry_operator.hpp"
#include "cosmosim/physics/baryonic_source_stage.hpp"
#include "cosmosim/physics/stellar_population.hpp"
#include "cosmosim/physics/supernova_feedback.hpp"
#include "cosmosim/core/stage_transaction.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"
#include "cosmosim/physics/element_diffusion.hpp"
#include "cosmosim/physics/stellar_wind.hpp"
#include "cosmosim/parallel/composition_exchange.hpp"
#include "cosmosim/parallel/distributed_mesh.hpp"
#include "workflows/internal/amr_migration_payload.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/mpi_stage_collective_decision.hpp"
#include "workflows/internal/particle_id_registry.hpp"
#include "workflows/internal/distributed_element_diffusion.hpp"
#include "workflows/internal/metal_diffusion_topology.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/retained_population_growth.hpp"
#include "workflows/internal/star_formation_geometry.hpp"
#include "workflows/internal/star_formation_gradients.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::workflows {
namespace {

constexpr std::size_t k_feedback_event_batch_max = 4096U;
constexpr std::uint64_t k_feedback_batch_bytes_per_star =
    static_cast<std::uint64_t>(sizeof(physics::StellarEvolutionStarBudget)) +
    static_cast<std::uint64_t>(sizeof(physics::StellarFeedbackEvent)) +
    static_cast<std::uint64_t>(sizeof(std::uint32_t));

template <class T>
void replaceSourceVectorCapacity(std::vector<T>& values, std::size_t capacity,
                                 std::string_view owner) {
  if (capacity <= values.capacity()) return;
  // Allocate and validate the replacement before relinquishing the old owner.
  // An over-grant or copy failure destroys only the replacement, leaving the
  // canonical input and its existing retained-capacity lease intact.
  std::vector<T> replacement;
  replacement.reserve(capacity);
  if (replacement.capacity() != capacity) {
    throw std::length_error(std::string(owner) + " allocator exceeded admitted capacity");
  }
  replacement.insert(replacement.end(), values.begin(), values.end());
  values.swap(replacement);
}

[[nodiscard]] double newtonGCodeFromUnits(const core::UnitSystem& units) {
  return core::newtonGravitationalConstantCode(units);
}

[[nodiscard]] physics::StarFormationConfig makeRuntimeStarFormationConfig(
    const core::SimulationConfig& config,
    const core::UnitSystem& units) {
  physics::StarFormationConfig result = physics::makeStarFormationConfig(config.physics);
  result.newton_g_code = newtonGCodeFromUnits(units);
  result.density_is_comoving = config.units.coordinate_frame == core::CoordinateFrame::kComoving;
  result.geometry_is_comoving = config.units.coordinate_frame == core::CoordinateFrame::kComoving;
  return result;
}

[[nodiscard]] physics::StellarFeedbackConfig makeRuntimeStellarFeedbackConfig(
    const core::SimulationConfig& config,
    const core::UnitSystem& units) {
  physics::StellarFeedbackConfig result =
      physics::makeStellarFeedbackConfig(config.physics);
  // Mass and metal return is part of stellar evolution even when energetic
  // feedback is disabled. Keep the deposition transaction enabled and zero only
  // the energy/momentum coupling in that case.
  result.enabled = config.physics.enable_stellar_evolution;
  if (!config.physics.enable_feedback) {
    result.epsilon_thermal = 0.0;
    result.epsilon_kinetic = 0.0;
    result.epsilon_momentum = 0.0;
  }
  constexpr double k_joule_per_erg = 1.0e-7;
  result.total_energy_code_per_erg = k_joule_per_erg /
      (units.mass_si_per_code * units.velocity_si_per_code *
       units.velocity_si_per_code);
  return result;
}

[[nodiscard]] std::shared_ptr<const physics::EffectiveMultiphaseEosTable>
makeRuntimeEffectiveEosTable(
    const core::SimulationConfig& config,
    const core::UnitSystem& units) {
  if (!config.physics.enable_star_formation ||
      config.physics.star_formation_model !=
          core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
    return nullptr;
  }
  return std::make_shared<const physics::EffectiveMultiphaseEosTable>(
      physics::makeEffectiveMultiphaseEosConfig(config.physics),
      units,
      physics::makeEffectiveIsmReferenceCoolingProvider(config.physics));
}

// Audit F4: records every gas cell mutated by feedback deposition together
// with its specific internal energy before the first mutation, so the single
// thermodynamic refresh contract can update P / c_s / T before the next CFL
// decision. Capacity is admitted through the MemoryGovernor BEFORE each event
// batch (old + new storage coexist during growth); recordGasCellMutation never
// allocates inside a deposition, so a refusal cannot leave a half-deposited
// batch.
class FeedbackThermodynamicTouchLog final
    : public physics::StellarFeedbackThermodynamicTouchSink {
 public:
  explicit FeedbackThermodynamicTouchLog(core::MemoryGovernor* governor) noexcept
      : m_governor(governor) {}

  void reserveAdditional(std::size_t additional_touches) {
    const std::size_t required = core::checkedSizeAdd(
        m_touches.size(), additional_touches, "feedback thermodynamic touch log size");
    if (required <= m_touches.capacity()) {
      return;
    }
    const std::size_t new_capacity = std::max(required, core::checkedSizeMultiply(
        m_touches.capacity(), 2U, "feedback thermodynamic touch log growth"));
    const std::uint64_t bytes = static_cast<std::uint64_t>(core::checkedSizeMultiply(
        new_capacity, sizeof(hydro::GasThermodynamicTouch), "feedback thermodynamic touch log bytes"));
    core::MemoryReservation replacement_reservation;
    if (m_governor != nullptr) {
      replacement_reservation = m_governor->reserve(
          core::MemoryClass::kScratchArena, bytes, "sources.feedback.thermodynamic_touch_log");
    }
    std::vector<hydro::GasThermodynamicTouch> replacement;
    replacement.reserve(new_capacity);
    replacement.insert(replacement.end(), m_touches.begin(), m_touches.end());
    m_touches.swap(replacement);
    if (m_governor != nullptr) {
      replacement_reservation.commit();
      m_reservation.release();
      m_reservation = std::move(replacement_reservation);
    }
  }

  void recordGasCellMutation(std::uint32_t cell_index, double internal_energy_before_code) override {
    if (m_touches.size() == m_touches.capacity()) {
      throw std::logic_error(
          "feedback thermodynamic touch log exceeded its pre-admitted capacity; "
          "reserveAdditional() must bound every deposition batch");
    }
    m_touches.push_back(hydro::GasThermodynamicTouch{
        .cell_row = cell_index, .internal_energy_before_code = internal_energy_before_code});
  }

  // Sort by row and keep the FIRST record of each row (its pre-stage energy).
  [[nodiscard]] std::span<const hydro::GasThermodynamicTouch> collapse() {
    std::stable_sort(m_touches.begin(), m_touches.end(),
        [](const hydro::GasThermodynamicTouch& a, const hydro::GasThermodynamicTouch& b) {
          return a.cell_row < b.cell_row;
        });
    const auto last = std::unique(m_touches.begin(), m_touches.end(),
        [](const hydro::GasThermodynamicTouch& a, const hydro::GasThermodynamicTouch& b) {
          return a.cell_row == b.cell_row;
        });
    m_touches.erase(last, m_touches.end());
    return m_touches;
  }

  [[nodiscard]] bool empty() const noexcept { return m_touches.empty(); }

 private:
  core::MemoryGovernor* m_governor = nullptr;
  core::MemoryReservation m_reservation;
  std::vector<hydro::GasThermodynamicTouch> m_touches;
};

template <class TFunction>
void prepareIdCollectively(const parallel::MpiContext& mpi_context,
                           std::string_view phase, TFunction&& prepare) {
  std::exception_ptr failure;
  try { prepare(); } catch (...) { failure = std::current_exception(); }
  mpi_context.rethrowCollectivePreparationFailure(failure, phase);
}

[[nodiscard]] std::uint64_t sourceParentIndexBytes(std::size_t patch_count) {
  return core::checkedMemoryBytesAdd(
      static_cast<std::uint64_t>(core::checkedSizeMultiply(
          patch_count, sizeof(std::uint64_t), "source parent index")),
      256U, "source parent index alignment");
}

[[nodiscard]] std::pmr::vector<std::uint64_t> buildSourceParentIndex(
    const core::SimulationState& state, std::pmr::memory_resource* resource) {
  std::pmr::vector<std::uint64_t> parents(resource);
  parents.reserve(state.patches.size());
  for (const std::uint64_t parent_id : state.patches.parent_patch_id) {
    if (parent_id != 0U) parents.push_back(parent_id);
  }
  std::sort(parents.begin(), parents.end());
  parents.erase(std::unique(parents.begin(), parents.end()), parents.end());
  return parents;
}

using IdVector = std::pmr::vector<std::uint64_t>;
using IdBuckets = std::pmr::vector<IdVector>;

struct ShardedUint64Exchange {
  explicit ShardedUint64Exchange(std::pmr::memory_resource* resource)
      : values(resource), recv_counts(resource), recv_displacements(resource) {}
  IdVector values;
  std::pmr::vector<std::size_t> recv_counts;
  std::pmr::vector<std::size_t> recv_displacements;
};

// The transport has O(ranks + configured round bytes) physical scratch.
// Logical received records are bounded by the caller's admitted source batch.
// No population-sized MPI round-plan object is retained.
[[nodiscard]] ShardedUint64Exchange exchangeShardedUint64Records(
    const parallel::MpiContext& mpi_context,
    const IdBuckets& values_by_rank,
    std::size_t record_width,
    std::pmr::memory_resource* resource) {
  const int world_size = mpi_context.worldSize();
  std::exception_ptr local_failure;
  try {
    if (record_width == 0U || world_size <= 0 ||
        values_by_rank.size() != static_cast<std::size_t>(world_size)) {
      throw std::invalid_argument("sharded uint64 exchange framing/rank extent mismatch");
    }
    for (const auto& values : values_by_rank) {
      if (values.size() % record_width != 0U) {
        throw std::invalid_argument("sharded uint64 exchange received a partial record");
      }
    }
  } catch (...) { local_failure = std::current_exception(); }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "sharded uint64 local framing validation");

  ShardedUint64Exchange result(resource);
  if (!mpi_context.isEnabled()) {
    result.values.assign(values_by_rank.front().begin(), values_by_rank.front().end());
    result.recv_counts.push_back(result.values.size());
    result.recv_displacements.push_back(0U);
    return result;
  }
#if COSMOSIM_ENABLE_MPI
  std::pmr::vector<std::uint64_t> send_counts64(resource), recv_counts64(resource);
  local_failure = nullptr;
  try {
    send_counts64.resize(static_cast<std::size_t>(world_size), 0U);
    recv_counts64.resize(static_cast<std::size_t>(world_size), 0U);
    for (std::size_t rank = 0; rank < values_by_rank.size(); ++rank) {
      send_counts64[rank] = core::checkedIntegralNarrow<std::uint64_t>(
          values_by_rank[rank].size(), "sharded uint64 logical send count");
    }
  } catch (...) { local_failure = std::current_exception(); }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "sharded uint64 count-buffer preparation");
  if (MPI_Alltoall(send_counts64.data(), 1, MPI_UINT64_T,
                   recv_counts64.data(), 1, MPI_UINT64_T,
                   MPI_COMM_WORLD) != MPI_SUCCESS) {
    throw std::runtime_error("sharded uint64 MPI_Alltoall count exchange failed");
  }

  std::pmr::vector<int> send_counts(resource), recv_counts(resource);
  std::pmr::vector<int> send_displacements(resource), recv_displacements(resource);
  std::pmr::vector<std::uint64_t> send_round_buffer(resource), recv_round_buffer(resource);
  std::uint64_t per_peer_cap = 0U;
  std::uint64_t global_round_count = 0U;
  local_failure = nullptr;
  try {
    const std::size_t ranks = static_cast<std::size_t>(world_size);
    result.recv_counts.resize(ranks, 0U);
    result.recv_displacements.resize(ranks, 0U);
    std::size_t total = 0U;
    for (std::size_t peer = 0; peer < ranks; ++peer) {
      result.recv_counts[peer] = core::checkedIntegralNarrow<std::size_t>(
          recv_counts64[peer], "sharded uint64 receive count");
      if (result.recv_counts[peer] % record_width != 0U) {
        throw std::runtime_error("sharded uint64 receive record framing mismatch");
      }
      result.recv_displacements[peer] = total;
      total = core::checkedSizeAdd(total, result.recv_counts[peer],
                                   "sharded uint64 receive total");
    }
    result.values.resize(total);
    // One element per peer is the minimum legal round. Both the aggregate
    // count and every displacement remain representable by MPI int.
    const std::size_t round_elements = std::max<std::size_t>(
        ranks, parallel::mpiTransportRoundLimitBytes() / sizeof(std::uint64_t));
    std::uint64_t local_max_peer = 1U;
    for (std::size_t peer = 0; peer < ranks; ++peer) {
      local_max_peer = std::max({local_max_peer, send_counts64[peer], recv_counts64[peer]});
    }
    per_peer_cap = std::min<std::uint64_t>(
        static_cast<std::uint64_t>(round_elements / ranks),
        static_cast<std::uint64_t>(std::numeric_limits<int>::max()) / ranks);
    // The common cap is finalized after the collective preparation gate.
    // Do not allocate round buffers until all ranks agree on this maximum.
    global_round_count = local_max_peer;
    const std::size_t capacity = core::checkedSizeMultiply(
        ranks, static_cast<std::size_t>(per_peer_cap), "sharded uint64 round capacity");
    send_counts.resize(ranks, 0); recv_counts.resize(ranks, 0);
    send_displacements.resize(ranks, 0); recv_displacements.resize(ranks, 0);
    (void)capacity;
  } catch (...) { local_failure = std::current_exception(); }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "sharded uint64 bounded payload preparation");
  const std::uint64_t global_max_peer = mpi_context.allreduceMaxUint64(global_round_count);
  per_peer_cap = std::max<std::uint64_t>(1U, std::min(per_peer_cap, global_max_peer));
  const std::size_t ranks = static_cast<std::size_t>(world_size);
  local_failure = nullptr;
  try {
    const std::size_t capacity = core::checkedSizeMultiply(ranks,
        static_cast<std::size_t>(per_peer_cap), "sharded uint64 round buffer capacity");
    send_round_buffer.resize(capacity); recv_round_buffer.resize(capacity);
    std::uint64_t local_round_count = 0U;
    for (std::size_t peer = 0; peer < ranks; ++peer) {
      const auto rounds = [per_peer_cap](std::uint64_t count) {
        return count / per_peer_cap + (count % per_peer_cap != 0U ? 1U : 0U);
      };
      local_round_count = std::max({local_round_count,
          rounds(send_counts64[peer]), rounds(recv_counts64[peer])});
    }
    global_round_count = local_round_count;
  } catch (...) { local_failure = std::current_exception(); }
  mpi_context.rethrowCollectivePreparationFailure(local_failure,
      "sharded uint64 round-buffer preparation");
  global_round_count = mpi_context.allreduceMaxUint64(global_round_count);
  for (std::uint64_t round = 0U; round < global_round_count; ++round) {
    const std::uint64_t offset = round * per_peer_cap;
    std::size_t send_total = 0U, recv_total = 0U;
    for (std::size_t peer = 0; peer < ranks; ++peer) {
      const std::uint64_t sc = offset < send_counts64[peer]
          ? std::min(per_peer_cap, send_counts64[peer] - offset) : 0U;
      const std::uint64_t rc = offset < recv_counts64[peer]
          ? std::min(per_peer_cap, recv_counts64[peer] - offset) : 0U;
      send_counts[peer] = static_cast<int>(sc);
      recv_counts[peer] = static_cast<int>(rc);
      send_displacements[peer] = static_cast<int>(send_total);
      recv_displacements[peer] = static_cast<int>(recv_total);
      if (sc != 0U) {
        std::memcpy(send_round_buffer.data() + send_total,
                    values_by_rank[peer].data() + static_cast<std::size_t>(offset),
                    static_cast<std::size_t>(sc) * sizeof(std::uint64_t));
      }
      send_total += static_cast<std::size_t>(sc);
      recv_total += static_cast<std::size_t>(rc);
    }
    if (MPI_Alltoallv(send_round_buffer.data(), send_counts.data(),
                      send_displacements.data(), MPI_UINT64_T,
                      recv_round_buffer.data(), recv_counts.data(),
                      recv_displacements.data(), MPI_UINT64_T,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("sharded uint64 bounded MPI_Alltoallv failed");
    }
    for (std::size_t peer = 0; peer < ranks; ++peer) {
      const std::size_t count = static_cast<std::size_t>(recv_counts[peer]);
      if (count != 0U) {
        std::memcpy(result.values.data() + result.recv_displacements[peer] +
                        static_cast<std::size_t>(offset),
                    recv_round_buffer.data() + static_cast<std::size_t>(recv_displacements[peer]),
                    count * sizeof(std::uint64_t));
      }
    }
  }
  return result;
#else
  throw std::runtime_error("distributed sharded particle-ID registry requires an MPI build");
#endif
}

// Merge two sorted, disjoint sequences into already-admitted vector capacity.
// Reverse traversal avoids an additional population-sized merge buffer.
template <class TContainer>
void appendSortedIds(TContainer& values, std::span<const std::uint64_t> incoming) {
  const std::size_t old_size = values.size();
  if (incoming.size() > values.capacity() - old_size) {
    throw std::length_error("sharded ID merge exceeds admitted capacity");
  }
  values.resize(old_size + incoming.size());
  std::size_t i = old_size, j = incoming.size(), k = values.size();
  while (j != 0U) {
    if (i != 0U && values[i - 1U] > incoming[j - 1U]) {
      values[--k] = values[--i];
    } else {
      values[--k] = incoming[--j];
    }
  }
}

[[nodiscard]] std::uint64_t checkedIdBytes(std::uint64_t count,
                                            std::uint64_t width,
                                            std::string_view label) {
  if (width != 0U && count > std::numeric_limits<std::uint64_t>::max() / width) {
    throw std::overflow_error(std::string(label) + ": byte multiplication overflow");
  }
  return count * width;
}

// The arena is a hard physical cap, not an optimistic peak estimate. The
// coefficients cover the explicitly reserved PMR arrays and their alignment.
// If an implementation requests more, null-upstream allocation fails safely.
[[nodiscard]] std::uint64_t idRoundWorkspaceBytes(std::uint64_t global_records,
                                                    std::uint64_t local_records,
                                                    std::uint64_t ranks,
                                                    std::uint64_t record_width) {
  std::uint64_t bytes = checkedIdBytes(global_records,
      record_width * sizeof(std::uint64_t) * 2U + 128U, "ID round records");
  bytes = core::checkedMemoryBytesAdd(bytes,
      checkedIdBytes(local_records, record_width * sizeof(std::uint64_t) * 2U + 64U,
                     "ID round local records"), "ID round workspace");
  bytes = core::checkedMemoryBytesAdd(bytes,
      checkedIdBytes(ranks, 512U, "ID round rank metadata"), "ID round workspace");
  const std::uint64_t transport = std::max<std::uint64_t>(
      checkedIdBytes(ranks, sizeof(std::uint64_t), "ID round minimum transport"),
      std::min<std::uint64_t>(parallel::mpiTransportRoundLimitBytes(),
          checkedIdBytes(checkedIdBytes(global_records, record_width,
              "ID round wire records"),
              checkedIdBytes(ranks, sizeof(std::uint64_t), "ID round peer width"),
              "ID round maximum transport")));
  bytes = core::checkedMemoryBytesAdd(bytes, checkedIdBytes(transport, 2U,
      "ID send/receive workspace"), "ID round workspace");
  return core::checkedMemoryBytesAdd(bytes, 4096U, "ID arena alignment allowance");
}

[[nodiscard]] std::size_t shardForUint64(std::uint64_t value, int world_size) {
  if (world_size <= 0) {
    throw std::invalid_argument("particle-ID sharding requires positive world size");
  }
  // Multiplicative mixing keeps sequential imported IDs from concentrating on
  // a subset of shards while remaining deterministic for any rank count.
  value ^= value >> 33U;
  value *= 0xff51afd7ed558ccdULL;
  value ^= value >> 33U;
  value *= 0xc4ceb9fe1a85ec53ULL;
  value ^= value >> 33U;
  return static_cast<std::size_t>(value % static_cast<std::uint64_t>(world_size));
}

class DistributedParticleIdRegistry final : public physics::ParticleIdPrecommit {
 public:
  explicit DistributedParticleIdRegistry(const parallel::MpiContext& mpi_context)
      : m_mpi_context(mpi_context) {}

  using GrowthHandler = void (*)(void*, core::SimulationState&, std::size_t, core::ParticleSpecies);
  void setGrowthHandler(void* owner, GrowthHandler handler) noexcept {
    m_growth_owner = owner;
    m_growth_handler = handler;
  }
  void setMemoryGovernor(core::MemoryGovernor* governor) noexcept {
    m_memory_governor = governor;
  }
  void finishPrecommit() noexcept override { m_result_reservation.release(); }

  // G1.6 rollback entry point (the only public invalidation path). The
  // whole-step coordinator restores SimulationState on every rank, so every
  // rank calls this together; the occupied shard is then rebuilt by the
  // collective initializeOccupiedShard() at the next precommit, from the
  // restored authoritative particle IDs. Any precommit result lease of the
  // rejected trial has no surviving consumer and is released here so the
  // next precommit's "no live consumer" guard holds on every rank.
  void invalidateForRollback() noexcept {
    m_result_reservation.release();
    invalidateOccupancyCache();
  }

  void preparePopulationGrowth(core::SimulationState& state, std::size_t count,
                               core::ParticleSpecies species) override {
    if (m_growth_handler != nullptr) m_growth_handler(m_growth_owner, state, count, species);
  }

  [[nodiscard]] std::uint64_t occupiedCapacityBytes() const {
    return checkedIdBytes(m_occupied_shard.capacity(), sizeof(std::uint64_t),
                          "ID occupied shard capacity");
  }
  [[nodiscard]] bool occupiedCapacityGoverned() const noexcept {
    return m_occupied_reservation.committed();
  }

  [[nodiscard]] std::vector<std::uint64_t> precommit(
      const core::SimulationState& state,
      std::span<const std::uint64_t> birth_keys) override {
    prepareIdCollectively(m_mpi_context, "particle-ID precommit entry", [&]() {
      if (m_result_reservation.committed()) {
        throw std::logic_error("particle-ID precommit result still has a live consumer");
      }
    });
    if (!m_initialized) initializeOccupiedShard(state.particle_sidecar.particle_id);
    validateBirthKeysDistributed(birth_keys);

    struct LocalCandidate {
      std::uint64_t particle_id = 0U;
      std::uint64_t birth_key = 0U;
      std::uint32_t collision_ordinal = 0U;
      bool resolved = false;
    };
    struct ReceivedCandidate {
      std::uint64_t particle_id = 0U;
      std::uint64_t birth_key = 0U;
      std::uint64_t origin_index = 0U;
      std::uint32_t ordinal = 0U;
      int origin_rank = 0;
    };
    const int world_size = m_mpi_context.worldSize();
    std::uint64_t local_count = 0U;
    prepareIdCollectively(m_mpi_context, "particle-ID count preparation", [&]() {
      local_count = core::checkedIntegralNarrow<std::uint64_t>(
          birth_keys.size(), "ID local precommit count");
    });
    const std::uint64_t global_count = m_mpi_context.allreduceSumUint64(local_count);
    if (global_count == 0U) return {};
    std::unique_ptr<core::GovernedScratchArena> outer;
    std::pmr::vector<LocalCandidate> candidates;
    IdVector batch_reserved_on_shard;
    prepareIdCollectively(m_mpi_context, "particle-ID outer workspace preparation", [&]() {
    const std::uint64_t outer_bytes = core::checkedMemoryBytesAdd(
        checkedIdBytes(local_count, sizeof(LocalCandidate) + sizeof(std::uint64_t) * 2U,
                       "ID outer local workspace"),
        core::checkedMemoryBytesAdd(
            checkedIdBytes(global_count, sizeof(std::uint64_t), "ID reserved candidates"),
            4096U, "ID outer alignment"), "ID outer workspace");
    outer = std::make_unique<core::GovernedScratchArena>(m_memory_governor, core::MemoryClass::kScratchArena,
                                      outer_bytes, "sources.ids.precommit");
    candidates = std::pmr::vector<LocalCandidate>(outer->resource());
    candidates.reserve(birth_keys.size());
    for (const std::uint64_t birth_key : birth_keys) {
      candidates.push_back(LocalCandidate{
          .particle_id = physics::starFormationParticleIdFromBirthKey(birth_key, 0U),
          .birth_key = birth_key});
    }
    batch_reserved_on_shard = IdVector(outer->resource());
    batch_reserved_on_shard.reserve(core::checkedIntegralNarrow<std::size_t>(
        global_count, "ID batch reservation capacity"));

    });
    for (std::uint32_t pass = 0U; pass < 1024U; ++pass) {
      // The arena is released after each collision pass; no allocator or
      // round-plan high-water can accumulate over the retry sequence.
      std::unique_ptr<core::GovernedScratchArena> scratch;
      IdBuckets requests;
      prepareIdCollectively(m_mpi_context, "particle-ID collision request preparation", [&]() {
      scratch = std::make_unique<core::GovernedScratchArena>(m_memory_governor, core::MemoryClass::kScratchArena,
          idRoundWorkspaceBytes(global_count, local_count,
                                static_cast<std::uint64_t>(world_size), 4U),
          "sources.ids.collision_round");
      requests = IdBuckets(scratch->resource());
      requests.resize(static_cast<std::size_t>(world_size));
      std::pmr::vector<std::size_t> request_counts(scratch->resource());
      request_counts.resize(requests.size(), 0U);
      for (const auto& candidate : candidates) {
        if (!candidate.resolved) ++request_counts[shardForUint64(candidate.particle_id, world_size)];
      }
      for (std::size_t rank = 0; rank < requests.size(); ++rank) {
        requests[rank].reserve(core::checkedSizeMultiply(request_counts[rank], 4U,
            "ID collision request records"));
      }
      for (std::size_t local_index = 0; local_index < candidates.size(); ++local_index) {
        const LocalCandidate& candidate = candidates[local_index];
        if (candidate.resolved) continue;
        auto& target = requests[shardForUint64(candidate.particle_id, world_size)];
        target.push_back(candidate.particle_id);
        target.push_back(candidate.birth_key);
        target.push_back(static_cast<std::uint64_t>(local_index));
        target.push_back(candidate.collision_ordinal);
      }
      });
      const ShardedUint64Exchange received = exchangeShardedUint64Records(
          m_mpi_context, requests, 4U, scratch->resource());
      IdBuckets decisions(scratch->resource());
      prepareIdCollectively(m_mpi_context, "particle-ID collision decision preparation", [&]() {
      std::pmr::vector<ReceivedCandidate> shard_candidates(scratch->resource());
      shard_candidates.reserve(received.values.size() / 4U);
      for (int origin_rank = 0; origin_rank < world_size; ++origin_rank) {
        const std::size_t begin = received.recv_displacements[static_cast<std::size_t>(origin_rank)];
        const std::size_t count = received.recv_counts[static_cast<std::size_t>(origin_rank)];
        for (std::size_t offset = 0; offset < count; offset += 4U) {
          const std::size_t index = begin + offset;
          const std::uint64_t ordinal64 = received.values[index + 3U];
          if (ordinal64 > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("particle-ID collision ordinal wire value overflow");
          }
          shard_candidates.push_back(ReceivedCandidate{
              .particle_id = received.values[index], .birth_key = received.values[index + 1U],
              .origin_index = received.values[index + 2U],
              .ordinal = static_cast<std::uint32_t>(ordinal64), .origin_rank = origin_rank});
        }
      }
      std::sort(shard_candidates.begin(), shard_candidates.end(),
                [](const ReceivedCandidate& lhs, const ReceivedCandidate& rhs) {
                  if (lhs.particle_id != rhs.particle_id) return lhs.particle_id < rhs.particle_id;
                  if (lhs.birth_key != rhs.birth_key) return lhs.birth_key < rhs.birth_key;
                  if (lhs.origin_rank != rhs.origin_rank) return lhs.origin_rank < rhs.origin_rank;
                  return lhs.origin_index < rhs.origin_index;
                });
      decisions.resize(static_cast<std::size_t>(world_size));
      std::pmr::vector<std::size_t> decision_counts(scratch->resource());
      decision_counts.resize(decisions.size(), 0U);
      for (const auto& candidate : shard_candidates) ++decision_counts[static_cast<std::size_t>(candidate.origin_rank)];
      for (std::size_t rank = 0; rank < decisions.size(); ++rank) {
        decisions[rank].reserve(core::checkedSizeMultiply(decision_counts[rank], 2U,
            "ID collision decision records"));
      }
      IdVector new_reserved(scratch->resource());
      new_reserved.reserve(shard_candidates.size());
      std::size_t begin = 0U;
      while (begin < shard_candidates.size()) {
        std::size_t end_group = begin + 1U;
        while (end_group < shard_candidates.size() &&
               shard_candidates[end_group].particle_id == shard_candidates[begin].particle_id) ++end_group;
        const std::uint64_t particle_id = shard_candidates[begin].particle_id;
        const bool collides_reserved = containsOccupied(particle_id) ||
            std::binary_search(batch_reserved_on_shard.begin(), batch_reserved_on_shard.end(), particle_id);
        const std::size_t first_rehash = collides_reserved ? begin : begin + 1U;
        if (!collides_reserved) {
          const ReceivedCandidate& winner = shard_candidates[begin];
          auto& result = decisions[static_cast<std::size_t>(winner.origin_rank)];
          result.push_back(winner.origin_index); result.push_back(0U);
          new_reserved.push_back(particle_id);
        }
        for (std::size_t index = first_rehash; index < end_group; ++index) {
          const ReceivedCandidate& collision = shard_candidates[index];
          auto& result = decisions[static_cast<std::size_t>(collision.origin_rank)];
          result.push_back(collision.origin_index); result.push_back(1U);
        }
        begin = end_group;
      }
      // Both sequences are sorted, so the merge is linear and allocation-free.
      appendSortedIds(batch_reserved_on_shard, new_reserved);
      });
      const ShardedUint64Exchange returned = exchangeShardedUint64Records(
          m_mpi_context, decisions, 2U, scratch->resource());
      prepareIdCollectively(m_mpi_context, "particle-ID collision response validation", [&]() {
      for (std::size_t index = 0; index < returned.values.size(); index += 2U) {
        const std::uint64_t local_index64 = returned.values[index];
        const std::uint64_t action = returned.values[index + 1U];
        if (local_index64 >= candidates.size() || action > 1U) {
          throw std::runtime_error("particle-ID shard returned an invalid collision decision");
        }
        LocalCandidate& candidate = candidates[static_cast<std::size_t>(local_index64)];
        if (action == 0U) { candidate.resolved = true; continue; }
        if (candidate.collision_ordinal == std::numeric_limits<std::uint32_t>::max()) {
          throw std::runtime_error("ParticleIdRegistry: deterministic collision ordinal overflow");
        }
        ++candidate.collision_ordinal;
        candidate.particle_id = physics::starFormationParticleIdFromBirthKey(
            candidate.birth_key, candidate.collision_ordinal);
      }
      });
      std::uint64_t local_unresolved = 0U;
      for (const auto& candidate : candidates) local_unresolved += candidate.resolved ? 0U : 1U;
      if (m_mpi_context.allreduceSumUint64(local_unresolved) == 0U) {
        // Admit the returned IDs before committing the occupied shard. The
        // caller releases this result lease after its birth transaction ends.
        core::MemoryReservation result_lease;
        std::vector<std::uint64_t> result;
        prepareIdCollectively(m_mpi_context,
            "particle-ID occupied-shard commit preparation", [&]() {
          if (m_memory_governor != nullptr) {
            result_lease = m_memory_governor->reserve(core::MemoryClass::kPhaseResident,
                checkedIdBytes(local_count, sizeof(std::uint64_t), "ID precommit result"),
                "sources.ids.precommit_result");
          }
          result.reserve(birth_keys.size());
          if (result.capacity() > birth_keys.size()) {
            throw std::length_error("ID precommit result exceeds admitted capacity");
          }
          for (const auto& candidate : candidates) result.push_back(candidate.particle_id);
          if (result_lease.pending()) result_lease.commit();
          reserveOccupiedCapacity(core::checkedSizeAdd(m_occupied_shard.size(),
              batch_reserved_on_shard.size(), "ID occupied append capacity"));
        });
        appendSortedIds(m_occupied_shard, batch_reserved_on_shard);
        m_result_reservation = std::move(result_lease);
        return result;
      }
    }
    throw std::runtime_error("ParticleIdRegistry: deterministic sharded collision resolution exhausted");
  }

 private:
  [[nodiscard]] bool containsOccupied(std::uint64_t id) const noexcept {
    return std::binary_search(m_occupied_shard.begin(), m_occupied_shard.end(), id);
  }

  void reserveOccupiedCapacity(std::size_t target) {
    if (target <= m_occupied_shard.capacity()) return;
    const std::size_t old_capacity = m_occupied_shard.capacity();
    const std::size_t growth = old_capacity / 2U;
    const std::size_t capacity = std::max(target, core::checkedSizeAdd(
        old_capacity, std::max<std::size_t>(growth, 1U), "ID occupied capacity growth"));
    const std::uint64_t bytes = checkedIdBytes(capacity, sizeof(std::uint64_t),
                                               "ID occupied replacement");
    auto replacement_lease = m_memory_governor == nullptr ? core::MemoryReservation{} :
        m_memory_governor->reserve(core::MemoryClass::kPersistentCache, bytes,
                                   "sources.ids.occupied_shard");
    std::vector<std::uint64_t> replacement;
    replacement.reserve(capacity);
    if (replacement.capacity() > capacity) {
      throw std::length_error("ID occupied allocator exceeded admitted replacement capacity");
    }
    replacement.insert(replacement.end(), m_occupied_shard.begin(), m_occupied_shard.end());
    if (replacement_lease.pending()) replacement_lease.commit();
    m_occupied_shard.swap(replacement);
    m_occupied_reservation = std::move(replacement_lease);
  }

  // Initial population validation is streamed in fixed local batches. The
  // persistent shard is O(owned hash shard), never replicated on every rank.
  // G1.6: a rolled-back step may have occupied IDs that no longer exist in
  // the restored state; the occupancy cache is rebuilt (collectively) from
  // the authoritative particle IDs at its next use.
  void invalidateOccupancyCache() noexcept { m_initialized = false; }

  void initializeOccupiedShard(std::span<const std::uint64_t> local_ids) {
    constexpr std::size_t k_initial_chunk = 4096U;
    const int world_size = m_mpi_context.worldSize();
    const std::uint64_t local_rounds = local_ids.size() / k_initial_chunk +
        (local_ids.size() % k_initial_chunk != 0U ? 1U : 0U);
    const std::uint64_t global_rounds = m_mpi_context.allreduceMaxUint64(local_rounds);
    m_occupied_shard.clear();
    try {
      for (std::uint64_t round = 0U; round < global_rounds; ++round) {
        const std::size_t begin = round < local_rounds
            ? static_cast<std::size_t>(round) * k_initial_chunk : local_ids.size();
        const std::size_t count = begin < local_ids.size()
            ? std::min(k_initial_chunk, local_ids.size() - begin) : 0U;
        const std::uint64_t global_count = m_mpi_context.allreduceSumUint64(count);
        std::unique_ptr<core::GovernedScratchArena> scratch;
        IdBuckets ids_by_rank;
        prepareIdCollectively(m_mpi_context, "particle-ID initial shard preparation", [&]() {
        scratch = std::make_unique<core::GovernedScratchArena>(m_memory_governor, core::MemoryClass::kScratchArena,
            idRoundWorkspaceBytes(global_count, count,
                                  static_cast<std::uint64_t>(world_size), 1U),
            "sources.ids.initial_shard");
        ids_by_rank = IdBuckets(scratch->resource());
        ids_by_rank.resize(static_cast<std::size_t>(world_size));
        std::pmr::vector<std::size_t> counts(scratch->resource());
        counts.resize(ids_by_rank.size(), 0U);
        for (std::size_t i = 0U; i < count; ++i) ++counts[shardForUint64(local_ids[begin + i], world_size)];
        for (std::size_t rank = 0; rank < counts.size(); ++rank) ids_by_rank[rank].reserve(counts[rank]);
        for (std::size_t i = 0U; i < count; ++i) ids_by_rank[shardForUint64(local_ids[begin + i], world_size)].push_back(local_ids[begin + i]);
        });
        const ShardedUint64Exchange received = exchangeShardedUint64Records(
            m_mpi_context, ids_by_rank, 1U, scratch->resource());
        IdVector shard_ids(scratch->resource());
        prepareIdCollectively(m_mpi_context, "particle-ID initial shard decode", [&]() {
          shard_ids.assign(received.values.begin(), received.values.end());
          std::sort(shard_ids.begin(), shard_ids.end());
        });
        bool invalid_local = (!shard_ids.empty() && shard_ids.front() == 0U) ||
            std::adjacent_find(shard_ids.begin(), shard_ids.end()) != shard_ids.end();
        if (!invalid_local) {
          for (const auto id : shard_ids) {
            if (containsOccupied(id)) { invalid_local = true; break; }
          }
        }
        if (m_mpi_context.allreduceSumUint64(invalid_local ? 1ULL : 0ULL) != 0U) {
          throw std::runtime_error("ParticleIdRegistry: zero or duplicate existing ID during sharded initialization");
        }
        std::exception_ptr local_failure;
        try { reserveOccupiedCapacity(core::checkedSizeAdd(m_occupied_shard.size(),
            shard_ids.size(), "ID initial shard append")); }
        catch (...) { local_failure = std::current_exception(); }
        m_mpi_context.rethrowCollectivePreparationFailure(local_failure,
            "particle-ID initial shard capacity preparation");
        appendSortedIds(m_occupied_shard, shard_ids);
      }
      m_initialized = true;
    } catch (...) {
      m_occupied_shard.clear();
      m_initialized = false;
      throw;
    }
  }

  void validateBirthKeysDistributed(std::span<const std::uint64_t> birth_keys) const {
    const int world_size = m_mpi_context.worldSize();
    std::uint64_t local_count = 0U;
    prepareIdCollectively(m_mpi_context, "particle-ID birth-key count preparation", [&]() {
      local_count = core::checkedIntegralNarrow<std::uint64_t>(
          birth_keys.size(), "ID birth-key count");
    });
    const std::uint64_t global_count = m_mpi_context.allreduceSumUint64(local_count);
    if (global_count == 0U) return;
    std::unique_ptr<core::GovernedScratchArena> scratch;
    IdBuckets keys_by_rank;
    prepareIdCollectively(m_mpi_context, "particle-ID birth-key workspace preparation", [&]() {
      scratch = std::make_unique<core::GovernedScratchArena>(m_memory_governor, core::MemoryClass::kScratchArena,
        idRoundWorkspaceBytes(global_count, local_count,
                              static_cast<std::uint64_t>(world_size), 1U),
        "sources.ids.birth_key_validation");
    keys_by_rank = IdBuckets(scratch->resource());
    keys_by_rank.resize(static_cast<std::size_t>(world_size));
    std::pmr::vector<std::size_t> counts(scratch->resource());
    counts.resize(keys_by_rank.size(), 0U);
    for (const auto key : birth_keys) ++counts[shardForUint64(key, world_size)];
    for (std::size_t rank = 0; rank < counts.size(); ++rank) keys_by_rank[rank].reserve(counts[rank]);
    for (const auto key : birth_keys) keys_by_rank[shardForUint64(key, world_size)].push_back(key);
    });
    const ShardedUint64Exchange received = exchangeShardedUint64Records(
        m_mpi_context, keys_by_rank, 1U, scratch->resource());
    IdVector shard_keys(scratch->resource());
    prepareIdCollectively(m_mpi_context, "particle-ID birth-key decode", [&]() {
      shard_keys.assign(received.values.begin(), received.values.end());
      std::sort(shard_keys.begin(), shard_keys.end());
    });
    const bool invalid_local = (!shard_keys.empty() && shard_keys.front() == 0U) ||
        std::adjacent_find(shard_keys.begin(), shard_keys.end()) != shard_keys.end();
    if (m_mpi_context.allreduceSumUint64(invalid_local ? 1ULL : 0ULL) != 0U) {
      throw std::runtime_error("ParticleIdRegistry: zero or duplicate immutable birth key across owner ranks");
    }
  }

  const parallel::MpiContext& m_mpi_context;
  core::MemoryGovernor* m_memory_governor = nullptr;
  void* m_growth_owner = nullptr;
  GrowthHandler m_growth_handler = nullptr;
  bool m_initialized = false;
  std::vector<std::uint64_t> m_occupied_shard;
  core::MemoryReservation m_occupied_reservation;
  core::MemoryReservation m_result_reservation;
};

[[nodiscard]] physics::BlackHoleAgnConfig makeRuntimeBlackHoleAgnConfig(
    const core::PhysicsConfig& physics_config,
    const core::UnitSystem& units) {
  physics::BlackHoleAgnConfig config = physics::makeBlackHoleAgnConfig(physics_config);
  config.proton_mass_code = config.proton_mass_si / units.mass_si_per_code;
  config.thomson_cross_section_code = config.thomson_cross_section_si /
      (units.length_si_per_code * units.length_si_per_code);
  config.newton_g_code = config.newton_g_si * units.mass_si_per_code *
      units.timeSiPerCode() * units.timeSiPerCode() /
      (units.length_si_per_code * units.length_si_per_code * units.length_si_per_code);
  config.speed_of_light_code = config.speed_of_light_si / units.velocity_si_per_code;
  return config;
}

using PatchCellGeometry = internal::StarFormationPatchCellGeometry;
using internal::starFormationDerivativeAtCell;
using internal::starFormationPatchCellGeometry;
using internal::starFormationPatchIsLeaf;

class FeedbackCellVolumeProvider final
    : public physics::StellarFeedbackCellVolumeProvider {
 public:
  explicit FeedbackCellVolumeProvider(const core::SimulationState& state)
      : m_state(state) {}

  [[nodiscard]] double cellVolumeCode(std::uint32_t cell_index) const override {
    if (cell_index >= m_state.cells.size() || cell_index >= m_state.gas_cells.size()) {
      return 0.0;
    }
    const PatchCellGeometry geometry =
        starFormationPatchCellGeometry(m_state, cell_index);
    if (geometry.valid) {
      return geometry.dx_stored * geometry.dy_stored * geometry.dz_stored;
    }
    const double density = m_state.gas_cells.density_code[cell_index];
    return density > 0.0 ? m_state.cells.mass_code[cell_index] / density : 0.0;
  }

 private:
  const core::SimulationState& m_state;
};

// T01 source-interval contract. One plan is derived per source dispatch from
// the integrator-owned CosmologicalStepFactors (coarse/block interval) and the
// dispatch's SourceSplitPhase; every source operator reads its clock from
// here, never from an ad-hoc local. Code time is linear in cosmic time
// (time_si_per_code), so code-time halves are physical-time halves; the
// half-interval scale factors are obtained by inverting the same background
// integral. Proper-coordinate (non-cosmological) runs use a == 1 throughout.
struct SourceStagePlan {
  core::SourceSplitPhase phase = core::SourceSplitPhase::kUnsplitCompatibility;
  // Coarse (or hierarchical-block) interval owned by the time integrator.
  double coarse_time_begin_code = 0.0;
  double coarse_time_end_code = 0.0;
  double coarse_dt_code = 0.0;
  // Smooth sub-interval advanced by this dispatch.
  double smooth_fraction_begin = 0.0;
  double smooth_fraction_end = 1.0;
  double smooth_time_begin_code = 0.0;
  double smooth_time_end_code = 0.0;
  double smooth_dt_code = 0.0;
  double smooth_scale_factor_begin = 1.0;
  double smooth_scale_factor_end = 1.0;
  // Physical-time midpoint of the smooth sub-interval (compatibility: end).
  double smooth_evaluation_scale_factor = 1.0;
  // Discrete event epoch: absolute end of the coarse interval, full duration.
  double event_time_code = 0.0;
  double event_scale_factor = 1.0;
  double event_dt_code = 0.0;
  // Ownership of this dispatch.
  bool global_smooth_thermochemistry = true;  // false: AMR levels own it
  bool smooth_transport = true;               // element diffusion + SGS shear
  bool discrete_events = false;
  bool opens_interval = false;
  bool closes_interval = false;
};

// T02 (patch 52) physical-time ledger of the discrete stochastic events of
// one source interval. Physical event time (model-defined threshold
// crossing, see physics::realizedEventChronology) is distinct from the
// publication epoch (interval end, the only legal synchronized deposition
// point). The lag = age_end - a_k >= 0 of every timed event is recorded;
// when physics.ssp_max_event_publication_lag_years > 0 it is enforced
// (fail closed) -- the time authority clips steps so it cannot be exceeded.
// Catch-up events of a legacy record (threshold already crossed at age_begin)
// are counted and exempt: their physical time precedes the interval.
struct DiscreteEventChronologyLedger {
  double interval_begin_code = 0.0;
  double interval_end_code = 0.0;  // publication epoch
  double lag_bound_years = 0.0;    // 0: bounded by the interval only
  // P1-TIME-01: physics.ssp_event_publication_mode = physical_event_barrier.
  // The time authority ended this interval at the next global physical event
  // time (+ k_ssp_event_barrier_tolerance_years); lag_bound_years then holds
  // that barrier resolution and is enforced for every forecast-covered event.
  bool physical_event_barrier = false;
  std::uint64_t timed_events = 0U;
  std::uint64_t catch_up_events = 0U;
  double earliest_event_fraction = 1.0;  // (a_k - age_b) / (age_e - age_b)
  double latest_event_fraction = 0.0;
  double max_publication_lag_years = 0.0;
  double sum_publication_lag_years = 0.0;
  // Patch 53: timed events of records that the step-start forecast did not
  // cover (cursor adopted in this interval); recorded, exempt from the bound.
  std::uint64_t unforecast_timed_events = 0U;

  // Policy label of this interval (patch 53, documentation in
  // docs/time_integration_event_chronology.md): the default 0 means "lag
  // bounded only by the accepted coarse interval", NOT "zero lag required".
  [[nodiscard]] const char* lagPolicy() const noexcept {
    if (physical_event_barrier) return "physical_event_barrier";
    return lag_bound_years > 0.0 ? "bounded_by_time_authority_clip" : "bounded_by_accepted_interval_only";
  }
  [[nodiscard]] double latestEventTimeCode() const noexcept {
    return timed_events == 0U ? interval_end_code
                              : interval_begin_code +
            latest_event_fraction * (interval_end_code - interval_begin_code);
  }

  [[nodiscard]] double earliestEventTimeCode() const noexcept {
    return timed_events == 0U ? interval_end_code
                              : interval_begin_code +
            earliest_event_fraction * (interval_end_code - interval_begin_code);
  }

  void accumulate(const physics::RealizedEventChronology& chronology, std::uint64_t particle_id, double age_begin_years,
                  double age_end_years, bool enforce_lag_bound) {
    catch_up_events += chronology.catch_up_count;
    const std::uint32_t timed = chronology.count - chronology.catch_up_count;
    if (timed == 0U) return;
    const double span = age_end_years - age_begin_years;
    if (!(span > 0.0) || chronology.earliest_age_years < age_begin_years ||
        chronology.latest_age_years > age_end_years) {
      throw std::logic_error("discrete stellar event physical age lies outside its own source interval");
    }
    timed_events += timed;
    const double lag = age_end_years - chronology.earliest_age_years;
    max_publication_lag_years = std::max(max_publication_lag_years, lag);
    sum_publication_lag_years += static_cast<double>(timed) * age_end_years - chronology.sum_age_years;
    earliest_event_fraction =
        std::min(earliest_event_fraction, (chronology.earliest_age_years - age_begin_years) / span);
    latest_event_fraction =
        std::max(latest_event_fraction, (chronology.latest_age_years - age_begin_years) / span);
    if (!enforce_lag_bound) unforecast_timed_events += timed;
    // Relative 1e-6 covers the code-time <-> stellar-years conversion
    // round-off between the forecast and the accepted interval.
    if (enforce_lag_bound && lag_bound_years > 0.0 && lag > lag_bound_years * (1.0 + 1.0e-6)) {
      throw std::runtime_error(
          "discrete SN event of star " + std::to_string(particle_id) + " would be published " +
          std::to_string(lag) + " yr after its physical time, above " +
          (physical_event_barrier
               ? std::string("the physical_event_barrier resolution = ")
               : std::string("physics.ssp_max_event_publication_lag_years = ")) +
          std::to_string(lag_bound_years) + " yr (time-authority event clip not honoured)");
    }
  }
};

// T02 (patch 52) canonical identity of one coarse/fine interface face as
// seen by one rank: the two stable gas_cell_ids (lo < hi), the face axis,
// the side of the face relative to the lo cell and the periodic-image flag.
// The hierarchy, the leaf identities and the rank ownership are frozen
// inside one coarse interval, so the exact sorted set a rank selects in the
// level pre half must equal the set it selects in the post-children closure.
struct CanonicalInterfaceFace {
  std::uint64_t lo_id = 0U;
  std::uint64_t hi_id = 0U;
  std::uint8_t axis = 0U;
  std::int8_t side = 0;
  std::uint8_t periodic = 0U;
  std::uint8_t remote = 0U;  // 1: MPI-mirrored face (the peer owner holds the other cell)

  [[nodiscard]] friend bool operator<(const CanonicalInterfaceFace& l, const CanonicalInterfaceFace& r) noexcept {
    if (l.lo_id != r.lo_id) return l.lo_id < r.lo_id;
    if (l.hi_id != r.hi_id) return l.hi_id < r.hi_id;
    if (l.axis != r.axis) return l.axis < r.axis;
    if (l.side != r.side) return l.side < r.side;
    if (l.periodic != r.periodic) return l.periodic < r.periodic;
    return l.remote < r.remote;
  }
  [[nodiscard]] friend bool operator==(const CanonicalInterfaceFace& l, const CanonicalInterfaceFace& r) noexcept {
    return l.lo_id == r.lo_id && l.hi_id == r.hi_id && l.axis == r.axis && l.side == r.side &&
        l.periodic == r.periodic && l.remote == r.remote;
  }
};

// Accumulates the observations of one Strang source interval (both global
// halves and every recursive level half) for the single controller
// publication at interval close. Plain values: rollback restores it.
struct SourceIntervalAccumulator {
  bool open = false;
  std::uint64_t sequence = 0U;            // intervals opened so far
  std::uint64_t step_index = std::numeric_limits<std::uint64_t>::max();
  bool recursive_levels = false;          // AMR levels own smooth thermochemistry
  std::uint64_t global_smooth_halves = 0U;
  std::uint64_t level_pre_halves = 0U;
  std::uint64_t level_post_halves = 0U;
  // T02: level halves that integrated level-clock diffusion/SGS faces.
  std::uint64_t level_transport_halves = 0U;
  // T02 (patch 51): per level, 0 = no pre half yet; 1 = coarse/fine split
  // closed (or none opened); n + 1 = n global coarse/fine faces opened by the
  // level pre half and awaiting the post-children synchronization.
  std::vector<std::uint64_t> pending_interface_sync;
  // T02 (patch 52): exact per-rank ledger of the coarse/fine faces each level
  // pre half opened (sorted canonical identities), consumed and compared by
  // the post-children closure. Bounded by the rank's coarse/fine face count
  // and admitted (m_interface_face_ledger_leases); empty outside an open
  // split, so the rollback copy taken at the step start is empty.
  std::vector<std::vector<CanonicalInterfaceFace>> pending_interface_faces;
  std::uint64_t level_interface_syncs = 0U;
  // T02 exactly-once discrete-event dispatch per interval.
  std::uint64_t discrete_event_dispatches = 0U;
  physics::ChemistryShadowReport shadow{};
  physics::ChemistryAdaptiveReport adaptive{};
};

class SourceRuntimeImpl final : public SourceRuntime {
 public:
  SourceRuntimeImpl(
      const core::SimulationConfig& config,
      const core::ModePolicy& mode_policy,
      const core::UnitSystem& units,
      std::uint32_t world_rank,
      const parallel::MpiContext& mpi_context,
      std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table,
      const RuntimeServices* runtime_services)
      : m_units(units),
        m_effective_eos_table(
            effective_eos_table != nullptr
                ? std::move(effective_eos_table)
                : makeRuntimeEffectiveEosTable(config, units)),
        m_star_formation(makeRuntimeStarFormationConfig(config, units), m_effective_eos_table),
        m_stellar_evolution(
            physics::makeStellarEvolutionConfig(config.physics),
            physics::loadStellarEvolutionTable(config.physics)),
        m_stellar_feedback(makeRuntimeStellarFeedbackConfig(config, units)),
        m_metal_diffusion(physics::makeMetalDiffusionConfig(config.physics)),
        m_black_hole(makeRuntimeBlackHoleAgnConfig(config.physics, units)),
        m_physics_config(config.physics),
        m_full_config(config),
        m_particle_id_registry(mpi_context),
        m_mpi_context(mpi_context),
        m_runtime_services(runtime_services),
        m_memory_governor(
            runtime_services != nullptr ? runtime_services->memory_governor : nullptr),
        m_world_rank(world_rank),
        m_coordinate_frame(config.units.coordinate_frame),
        m_hydro_boundary(mode_policy.hydro_boundary),
        m_bh_enabled(config.physics.enable_black_hole_agn),
        m_bh_seeding_requested(config.physics.bh_enable_seeding),
        m_is_cosmological(
            config.mode.mode == core::SimulationMode::kCosmoCube ||
            config.mode.mode == core::SimulationMode::kZoomIn) {
    m_particle_id_registry.setMemoryGovernor(m_memory_governor);
    m_particle_id_registry.setGrowthHandler(this,
        &SourceRuntimeImpl::preparePopulationGrowthCallback);
    if (m_bh_seeding_requested) {
      if (!m_bh_enabled) {
        throw std::invalid_argument(
            "physics.bh_enable_seeding=true requires physics.enable_black_hole_agn=true");
      }
      throw std::runtime_error(
          "BH seeding requested, but ReferenceWorkflow has no authoritative halo/candidate provider; "
          "disable bh_enable_seeding or supply a production candidate provider before claiming seeding support");
    }

    const double hubble_si = config.cosmology.hubble_param *
        core::constants::k_hubble_100_km_s_mpc_si;
    const double hubble_code = hubble_si * units.timeSiPerCode();
    const double g_code = newtonGCodeFromUnits(units);
    m_mean_baryon_density0_code = 3.0 * hubble_code * hubble_code /
        (8.0 * core::constants::k_pi * g_code) * config.cosmology.omega_baryon;
    m_hydro_density_floor_code = config.numerics.hydro_density_floor_code;
    m_hydro_pressure_floor_code = config.numerics.hydro_pressure_floor_code;
    // CHUI Part 2: resolve the thermochemistry model, its data products and
    // the composition layout ONCE (never in cell loops). Missing production
    // data stops construction with SourceModelUnavailableError.
    m_box_extent_comoving = {config.cosmology.box_size_x_mpc_comoving, config.cosmology.box_size_y_mpc_comoving,
                             config.cosmology.box_size_z_mpc_comoving};
    if (config.physics.enable_stellar_evolution &&
        config.physics.stellar_population_model != core::StellarPopulationModelConfig::kLegacyTable) {
      physics::SspConfig ssp;
      ssp.imf_mass_min_msun = config.physics.ssp_imf_mass_min_msun;
      ssp.imf_mass_max_msun = config.physics.ssp_imf_mass_max_msun;
      ssp.cc_mass_min_msun = config.physics.ssp_cc_mass_min_msun;
      ssp.cc_mass_max_msun = config.physics.ssp_cc_mass_max_msun;
      ssp.snia_per_msun = config.physics.ssp_snia_per_msun;
      ssp.snia_delay_min_years = config.physics.ssp_snia_delay_min_years;
      ssp.snia_delay_max_years = config.physics.ssp_snia_delay_max_years;
      ssp.energy_per_sn_erg = config.physics.ssp_energy_per_sn_erg;
      ssp.reference_net_metal_yield = config.physics.ssp_reference_net_metal_yield;
      ssp.stochastic_events = config.physics.ssp_stochastic_sn_events;
      if (config.physics.stellar_population_model == core::StellarPopulationModelConfig::kChabrierYieldTableV1) {
        ssp.yield_mode = physics::SspYieldMode::kTable;
        m_ssp_yield_table = std::make_shared<const physics::SspYieldTable>(physics::SspYieldTable::loadFromFile(
            config.physics.ssp_yield_table_path, config.physics.ssp_yield_table_sha256));
      }
      m_ssp_model = std::make_unique<physics::SspModel>(ssp, m_ssp_yield_table.get());
      if (config.physics.ssp_wind_model == core::SspWindModelConfig::kTableV1) {
        m_wind_table = std::make_shared<const physics::StellarWindTable>(physics::StellarWindTable::loadFromFile(
            config.physics.ssp_wind_table_path, config.physics.ssp_wind_table_sha256));
      }
    }
    if (physics::part2ThermochemistrySelected(config.physics)) {
      m_baryonic_data = physics::loadBaryonicDataProducts(config.physics);
      if (config.physics.ism_regime_model == core::IsmRegimeModelConfig::kEffectiveOrResolved) {
        m_source_stage_config.regime = physics::makeIsmRegimeConfig(config.physics);
      }
      if (config.physics.dust_model != core::DustModelConfig::kNone) {
        m_source_stage_config.dust = physics::makeDustConfig(config.physics);
      }
      if (config.physics.sgs_turbulence_model != core::SgsTurbulenceModelConfig::kNone) {
        m_source_stage_config.sgs = physics::makeSgsConfig(config.physics);
      }
      m_thermochemistry = std::make_unique<physics::ThermochemistryOperator>(
          physics::makeThermochemistryOperatorConfig(config, units, m_baryonic_data),
          physics::CompositionLayout::make(physics::compositionLayoutRequestFromConfig(config.physics)));
    }
  }

  [[nodiscard]] bool fillAmrCoolingTimes(
      core::StepContext& context, std::span<double> cooling_time_code) override {
    if (m_thermochemistry == nullptr) return false;
    if (cooling_time_code.size() != context.state.cells.size()) {
      throw std::invalid_argument(
          "AMR thermochemistry cooling-time output must cover every dense gas row");
    }
    (void)physics::ensureGasCompositionInitialized(context.state, m_physics_config);
    const double evaluation_scale_factor = m_is_cosmological
        ? std::max(context.timeline_step.scale_factor_end, 1.0e-12)
        : 1.0;
    const physics::SourceInterval interval =
        makeSourceInterval(context, evaluation_scale_factor);
    const double seconds_per_code = m_units.timeSiPerCode();
    if (!(seconds_per_code > 0.0) || !std::isfinite(seconds_per_code)) {
      throw std::runtime_error("invalid source time unit while evaluating AMR cooling length");
    }
    for (std::uint32_t row = 0; row < context.state.cells.size(); ++row) {
      const double seconds = m_thermochemistry->coolingTimeSeconds(
          context.state, row, interval);
      cooling_time_code[row] = std::isfinite(seconds)
          ? seconds / seconds_per_code
          : std::numeric_limits<double>::infinity();
    }
    return true;
  }

  // Patch 53 (C5/V01): one runtime record of every loaded external data
  // product with its checksum, declared verification state and declared
  // model limitations (IMF/metallicity coverage of the wind table, unverified
  // C02 rates, UVB compatibility decision of the metal table). The normalized
  // configuration (paths + SHA-256, opt-ins) remains the restart/provenance
  // authority; this event makes the scientific assumptions explicit in the
  // run log. Diagnostic only: never mutates state.
  void reportDataAssumptionsOnce(const core::StepContext& context) {
    if (m_data_assumptions_reported || context.profiler_session == nullptr) return;
    m_data_assumptions_reported = true;
    std::unordered_map<std::string, std::string> payload;
    if (m_ssp_yield_table != nullptr) {
      payload.emplace("ssp_yield_table_sha256", m_ssp_yield_table->sha256());
      payload.emplace("ssp_yield_table_citation", m_ssp_yield_table->citation());
      payload.emplace("ssp_yield_table_license", m_ssp_yield_table->license());
    }
    if (m_wind_table != nullptr) {
      payload.emplace("ssp_wind_table_source", m_wind_table->source());
      payload.emplace("ssp_wind_table_metallicities", std::to_string(m_wind_table->metallicityCount()));
      payload.emplace("ssp_wind_table_limitations", m_wind_table->limitations());
    }
    if (m_baryonic_data.uv_background != nullptr) {
      payload.emplace("uv_background_table_sha256", m_baryonic_data.uv_background->provenance().sha256_hex);
      payload.emplace("uv_background_metal_photorates", "none (TREECOOL supplies H/He rates only)");
    }
    if (m_baryonic_data.metal_cooling_table != nullptr) {
      payload.emplace("metal_cooling_uv_background_compatible",
                           m_baryonic_data.uv_background_compatibility.compatible ? "true" : "false");
      payload.emplace("metal_cooling_uv_background_mismatch_accepted",
                           m_baryonic_data.uv_background_compatibility.accepted_mismatch ? "true" : "false");
    }
    if (m_baryonic_data.metal_ion_rate_table != nullptr) {
      const physics::MetalIonRateTable& ions = *m_baryonic_data.metal_ion_rate_table;
      payload.emplace("metal_ion_rate_table_provenance", ions.provenance());
      payload.emplace("metal_ion_rate_table_license", ions.license());
      payload.emplace("metal_ion_rate_table_verification",
                           ions.declaredVerified() ? "verified_against_source (declared)" : "unverified");
      payload.emplace("metal_ion_energy_convention", ions.energyConvention());
    }
    if (payload.empty()) return;
    context.profiler_session->recordEvent(core::RuntimeEvent{
        .event_kind = "sources.data_assumptions",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "physics.data",
        .step_index = context.integrator_state.step_index,
        .simulation_time_code = context.integrator_state.current_time_code,
        .scale_factor = context.integrator_state.current_scale_factor,
        .message = "external physics data products, checksums, verification state and declared limitations",
        .payload = std::move(payload)});
  }

  // T02 (patch 52): event-aware time-authority input (see SourceRuntime).
  // Deterministic legacy population record (stars created before Part 2 or
  // restored without a record): primordial H/He + solar-pattern metals at the
  // stored birth metallicity, unknown birth time, immutable birth key stream.
  // Shared by the SSP stage and the time-authority forecast so both see the
  // identical thresholds.
  [[nodiscard]] physics::StellarPopulationRecord legacyStellarPopulationRecord(
      const core::SimulationState& state, std::uint32_t s) const {
    const double z = std::clamp(state.star_particles.metallicity_mass_fraction[s], 0.0, 0.5);
    const auto masses = physics::initialElementMasses(1.0, z, physics::primordialPolicyFromConfig(m_physics_config));
    physics::StellarPopulationRecord record =
        physics::makeStellarPopulationRecord(masses, 0.0, state.star_particles.birth_key[s]);
    record.chronology_initialized = 0U;
    return record;
  }

  // P1-TIME-01: the enforced publication-lag bound of the configured
  // contract. coarse_endpoint: the optional configured lag (0 = bounded only
  // by the accepted interval). physical_event_barrier: the fixed barrier
  // resolution, so the global step ends at the next physical event time.
  [[nodiscard]] bool physicalEventBarrier() const noexcept {
    return m_physics_config.ssp_event_publication_mode ==
        core::SspEventPublicationModeConfig::kPhysicalEventBarrier;
  }
  [[nodiscard]] double enforcedEventLagBoundYears() const noexcept {
    return physicalEventBarrier() ? core::k_ssp_event_barrier_tolerance_years
                                  : m_physics_config.ssp_max_event_publication_lag_years;
  }

  [[nodiscard]] double discreteEventStepEndLimitTimeCode(const core::SimulationState& state,
                                                         double current_time_code) const override {
    const double lag_years = enforcedEventLagBoundYears();
    double limit = std::numeric_limits<double>::infinity();
    if (!(lag_years > 0.0) || m_ssp_model == nullptr || !m_ssp_model->config().stochastic_events ||
        state.star_particles.size() == 0U) {
      return limit;
    }
    constexpr double k_seconds_per_year = 31557600.0;
    const double years_per_code = m_units.timeSiPerCode() / k_seconds_per_year;
    const double msun_code = m_units.massSiToCode(physics::k_solar_mass_g * 1.0e-3);
    if (!(std::isfinite(years_per_code) && years_per_code > 0.0 && std::isfinite(msun_code) && msun_code > 0.0)) {
      throw std::runtime_error("discrete event step limit requires finite positive time and mass units");
    }
    const physics::StellarPopulationStore store(state);
    for (std::uint32_t s = 0; s < state.star_particles.size(); ++s) {
      const std::uint32_t p = state.star_particles.particle_index[s];
      const std::uint64_t pid = state.particle_sidecar.particle_id[p];
      // P1-TIME-01 (post patch 53): legacy stars are forecast too. A star
      // without a record gets, in its next interval, exactly the record the
      // stage constructs deterministically from immutable birth data
      // (legacyStellarPopulationRecord); a record without an absolute cursor
      // adopts the next interval begin (= current_time_code) as its cursor.
      // Thresholds already crossed at age_now are catch-up events (counted,
      // exempt: their physical time precedes the interval); they are
      // advanced past here so the forecast is the first TIMED event.
      physics::StellarPopulationRecord record;
      bool legacy = false;
      if (store.has(pid)) {
        record = store.get(pid);
      }
      if (record.initialized == 0U) {
        record = legacyStellarPopulationRecord(state, s);
        legacy = true;
      } else if (record.chronology_initialized == 0U) {
        legacy = true;
      }
      const double age_last = std::max(state.star_particles.stellar_age_years_last[s], 0.0);
      const double age_now = legacy
          ? age_last
          : age_last + std::max(current_time_code - record.chronology_time_code, 0.0) * years_per_code;
      const double birth_msun = state.star_particles.birth_mass_code[s] / msun_code;
      if (legacy) {
        double z = 0.0;
        for (std::size_t e = 2; e < physics::k_chemical_element_count; ++e) z += record.birth_element_mass_fraction[e];
        (void)physics::realizeEvents(record, 0U, m_ssp_model->ccCumulativeNumber(age_now, z) * birth_msun);
        (void)physics::realizeEvents(record, 1U, m_ssp_model->iaCumulativeNumber(age_now) * birth_msun);
      }
      const double next_age = physics::nextStellarEventAgeYears(
          *m_ssp_model, record, age_now, birth_msun);
      if (!std::isfinite(next_age)) continue;
      const double event_time_code = current_time_code + (next_age - age_now) / years_per_code;
      limit = std::min(limit, event_time_code + lag_years / years_per_code);
    }
    return limit;
  }

  // Drops every face-identity ledger and its admission (storage first).
  void resetInterfaceFaceLedger() noexcept {
    for (auto& ledger : m_source_interval.pending_interface_faces) {
      std::vector<CanonicalInterfaceFace>().swap(ledger);
    }
    for (auto& lease : m_interface_face_ledger_leases) lease.release();
    std::vector<CanonicalInterfaceFace>().swap(m_last_coarse_fine_faces);
    m_last_coarse_fine_faces_lease.release();
    m_face_ledger_failure = nullptr;
  }

  void captureRollbackCheckpoint() override {
    ControllerCheckpoint checkpoint;
    if (m_ssp_model != nullptr) {
      checkpoint.ssp_table_clamps = m_ssp_model->table_clamps;
      checkpoint.ssp_negative_ejecta_clamps = m_ssp_model->negative_ejecta_clamps;
      checkpoint.ssp_untracked_mass_residual_msun = m_ssp_model->untracked_mass_residual_msun;
      checkpoint.ssp_nonmonotone_return_clamps = m_ssp_model->nonmonotone_return_clamps;
    }
    checkpoint.wind_metallicity_clamps = m_wind_metallicity_clamps;
    checkpoint.sf_gradient_unavailable_cells = m_sf_gradient_unavailable_cells;
    checkpoint.sf_gradient_coarse_fine_cells = m_sf_gradient_coarse_fine_cells;
    checkpoint.last_thermochemistry_ledger = m_last_thermochemistry_ledger;
    checkpoint.last_feedback_ledger = m_last_feedback_ledger;
    checkpoint.last_event_chronology = m_last_event_chronology;
    checkpoint.last_ism_regime_report = m_last_ism_regime_report;
    checkpoint.last_chemistry_adaptive_report = m_last_chemistry_adaptive_report;
    checkpoint.last_dust_ledger = m_last_dust_ledger;
    checkpoint.last_sgs_ledger = m_last_sgs_ledger;
    checkpoint.last_shadow_report = m_last_shadow_report;
    checkpoint.source_interval = m_source_interval;
    checkpoint.last_skipped_effective_rows = m_last_skipped_effective_rows;
    checkpoint.last_element_diffusion = m_last_element_diffusion;
    m_rollback_checkpoint = std::move(checkpoint);
  }

  void restoreRollbackCheckpoint() override {
    if (!m_rollback_checkpoint.has_value()) {
      throw std::logic_error("source runtime rollback restore without a captured checkpoint");
    }
    const ControllerCheckpoint& c = *m_rollback_checkpoint;
    if (m_ssp_model != nullptr) {
      m_ssp_model->table_clamps = c.ssp_table_clamps;
      m_ssp_model->negative_ejecta_clamps = c.ssp_negative_ejecta_clamps;
      m_ssp_model->untracked_mass_residual_msun = c.ssp_untracked_mass_residual_msun;
      m_ssp_model->nonmonotone_return_clamps = c.ssp_nonmonotone_return_clamps;
    }
    m_wind_metallicity_clamps = c.wind_metallicity_clamps;
    m_sf_gradient_unavailable_cells = c.sf_gradient_unavailable_cells;
    m_sf_gradient_coarse_fine_cells = c.sf_gradient_coarse_fine_cells;
    m_last_thermochemistry_ledger = c.last_thermochemistry_ledger;
    m_last_feedback_ledger = c.last_feedback_ledger;
    m_last_event_chronology = c.last_event_chronology;
    m_last_ism_regime_report = c.last_ism_regime_report;
    m_last_chemistry_adaptive_report = c.last_chemistry_adaptive_report;
    m_last_dust_ledger = c.last_dust_ledger;
    m_last_sgs_ledger = c.last_sgs_ledger;
    m_last_shadow_report = c.last_shadow_report;
    // A rejected trial may have died between the pre and post halves with an
    // open interval; the checkpoint was taken at the step rollback point, so
    // restoring it also discards every partial shadow/adaptive observation.
    m_source_interval = c.source_interval;
    resetInterfaceFaceLedger();  // the restored interval holds no open split
    m_last_skipped_effective_rows = c.last_skipped_effective_rows;
    m_last_element_diffusion = c.last_element_diffusion;
    m_particle_id_registry.invalidateForRollback();
    // The chemistry promotion bias is re-derived from the restored sidecar
    // before the next fidelity decision (G0.4); nothing to restore here.
  }

  void discardRollbackCheckpoint() noexcept override { m_rollback_checkpoint.reset(); }

  [[nodiscard]] core::MemoryReport memoryReport() const override {
    core::MemoryReportBuilder builder;
    const auto add_container = [&builder](
        core::MemorySubsystem subsystem,
        core::MemoryClass memory_class,
        std::string_view label,
        const auto& container,
        bool governed) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = subsystem,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = memory_class,
          .label = std::string(label),
          .current_size_bytes = core::currentSizeBytesForContainer(container),
          .owned_capacity_bytes = core::ownedCapacityBytesForContainer(container),
          .high_water_bytes = core::ownedCapacityBytesForContainer(container),
          .estimate_only = false,
          .governed_commitment = governed,
      });
    };
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.star_formation.contiguous_cell_batch",
        m_contiguous_cell_batch,
        m_contiguous_cell_batch_reservation.committed());
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.star_formation.input_batch",
        m_star_formation_inputs,
        m_star_formation_input_reservation.committed());
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.stellar_evolution.contiguous_star_batch",
        m_contiguous_star_batch,
        m_contiguous_star_batch_reservation.committed());
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.stellar_feedback.event_batch",
        m_feedback_events,
        m_feedback_event_reservation.committed());
    // Gate-4 step-local governed budgets are admitted per source execution and
    // released at stage end; report the batch peak as governed transient.
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kScratch,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kPhaseResident,
        .label = "sources.stellar_feedback.spatial_index",
        .current_size_bytes = m_feedback_spatial_index.ownedCapacityBytes(),
        .owned_capacity_bytes = m_feedback_spatial_index.ownedCapacityBytes(),
        .high_water_bytes = m_feedback_spatial_index.highWaterBytes(),
        .estimate_only = false,
        .governed_commitment = m_feedback_index_reservation.committed(),
    });
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kParticles,
        .lifetime = core::MemoryLifetime::kPersistent,
        .memory_class = core::MemoryClass::kPersistentCache,
        .label = "sources.ids.occupied_shard",
        .owned_capacity_bytes = m_particle_id_registry.occupiedCapacityBytes(),
        .high_water_bytes = m_particle_id_registry.occupiedCapacityBytes(),
        .estimate_only = false,
        .governed_commitment = m_particle_id_registry.occupiedCapacityGoverned(),
        .uncertainty_note = "one sorted immutable-ID cache per hash shard; no rank-global replication"});
    const bool smooth_transport_governed =
        m_diffusion_phase_reservation.committed();
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.metal_diffusion.owned_leaf_mask",
        m_owned_leaf_mask,
        smooth_transport_governed);
    // These two containers are move-replaced by each topology rebuild. Report
    // historical retained-capacity high-water explicitly rather than letting a
    // smaller later topology erase the evidence of a larger residency peak.
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kScratch,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kPhaseResident,
        .label = "sources.metal_diffusion.rho_kappa",
        .current_size_bytes = core::currentSizeBytesForContainer(
            m_diffusion_rho_kappa_code),
        .owned_capacity_bytes = core::ownedCapacityBytesForContainer(
            m_diffusion_rho_kappa_code),
        .high_water_bytes = std::max(
            m_diffusion_rho_kappa_high_water_bytes,
            core::ownedCapacityBytesForContainer(m_diffusion_rho_kappa_code)),
        .estimate_only = false,
        .governed_commitment = smooth_transport_governed,
    });
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kScratch,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kPhaseResident,
        .label = "sources.metal_diffusion.faces",
        .current_size_bytes = core::currentSizeBytesForContainer(m_diffusion_faces),
        .owned_capacity_bytes = core::ownedCapacityBytesForContainer(m_diffusion_faces),
        .high_water_bytes = std::max(
            m_diffusion_faces_high_water_bytes,
            core::ownedCapacityBytesForContainer(m_diffusion_faces)),
        .estimate_only = false,
        .governed_commitment = m_diffusion_faces_reservation.committed(),
    });
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kScratch,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kPhaseResident,
        .label = "sources.metal_diffusion.workspace",
        .current_size_bytes = m_diffusion_workspace.ownedCapacityBytes(),
        .owned_capacity_bytes = m_diffusion_workspace.ownedCapacityBytes(),
        .high_water_bytes = m_diffusion_workspace.highWaterBytes(),
        .estimate_only = false,
        .governed_commitment = smooth_transport_governed,
    });
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kScratch,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kPhaseResident,
        .label = "sources.sgs_shear.staging",
        .current_size_bytes = m_sgs_shear_workspace.ownedCapacityBytes(),
        .owned_capacity_bytes = m_sgs_shear_workspace.ownedCapacityBytes(),
        .high_water_bytes = m_sgs_shear_workspace.ownedCapacityBytes(),
        .estimate_only = false,
        .governed_commitment = smooth_transport_governed,
    });
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kScratch,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kScratchArena,
        .label = "sources.metal_diffusion.topology_construction_scratch",
        .owned_capacity_bytes = 0U,
        .high_water_bytes = m_diffusion_topology_scratch_high_water_bytes,
        .estimate_only = true,
        .governed_commitment = false,
        .uncertainty_note =
            "measured container-capacity high-water for released topology construction scratch",
    });
    return std::move(builder).finish();
  }

  [[nodiscard]] RuntimeTaskMemoryEstimate estimateMemory(
      const core::SimulationState& state) const override {
    if (!m_star_formation.config().enabled &&
        !m_stellar_evolution.config().enabled &&
        !m_metal_diffusion.config().enabled &&
        !m_source_stage_config.sgs.has_value() && !m_bh_enabled) {
      return {0U, true, {}};
    }
    const std::uint64_t cell_count = static_cast<std::uint64_t>(state.cells.size());
    const std::uint64_t star_count = static_cast<std::uint64_t>(state.star_particles.size());
    const auto multiply = [](std::uint64_t count, std::uint64_t width,
                             std::string_view label) {
      if (width != 0U && count > std::numeric_limits<std::uint64_t>::max() / width) {
        throw std::overflow_error(std::string(label) + ": source memory estimate overflow");
      }
      return count * width;
    };
    const auto growth = [&](std::uint64_t count, std::uint64_t width,
                            std::uint64_t retained, std::string_view label) {
      return count > retained ? multiply(count, width, label) : 0U;
    };
    std::uint64_t peak = 0U;
    if (m_star_formation.config().enabled) {
      constexpr std::uint64_t k_batch = 4096U;
      const std::uint64_t cells = std::min(cell_count, k_batch);
      const std::uint64_t input_bytes = growth(
          cells, sizeof(physics::StarFormationCellInput),
          m_star_formation_inputs.capacity(), "source input batch");
      const std::uint64_t row_bytes = growth(
          cells, sizeof(std::uint32_t), m_contiguous_cell_batch.capacity(),
          "source contiguous cell batch");
      peak = core::checkedMemoryBytesAdd(input_bytes, row_bytes,
                                         "source formation staging");
    }
    if (m_stellar_evolution.config().enabled) {
      const std::uint64_t batch = std::min(star_count,
          static_cast<std::uint64_t>(k_feedback_event_batch_max));
      // Gate-4: step-local governed budget buffer (admitted per execution,
      // reused across batches within the execution, released at stage end).
      // Counters-only feedback path avoids a second per-star report vector.
      std::uint64_t evolution = multiply(batch,
          sizeof(physics::StellarEvolutionStarBudget), "source evolution budgets");
      evolution = core::checkedMemoryBytesAdd(evolution,
          growth(batch, sizeof(physics::StellarFeedbackEvent),
                 m_feedback_events.capacity(), "source feedback events"),
          "source feedback staging");
      evolution = core::checkedMemoryBytesAdd(evolution,
          growth(batch, sizeof(std::uint32_t), m_contiguous_star_batch.capacity(),
                 "source contiguous star batch"), "source feedback staging");
      evolution = core::checkedMemoryBytesAdd(evolution,
          multiply(cell_count, sizeof(std::uint32_t), "source spatial index rebuild"),
          "source feedback staging");
      // Sparse active rows have their own physical arena admission.
      peak = std::max(peak, evolution);
    }
    if ((m_metal_diffusion.config().enabled ||
         m_source_stage_config.sgs.has_value()) && cell_count != 0U) {
      std::uint64_t transport = core::checkedMemoryBytesAdd(
          m_metal_diffusion.config().enabled
              ? m_metal_diffusion.requiredWorkspaceBytes(state.cells.size())
              : 0U,
          multiply(cell_count, sizeof(std::uint8_t), "source transport mask"),
          "source smooth transport phase");
      transport = core::checkedMemoryBytesAdd(
          transport,
          multiply(cell_count, sizeof(double), "source transport strain/rho-kappa"),
          "source smooth transport phase");
      if (m_source_stage_config.sgs.has_value()) {
        transport = core::checkedMemoryBytesAdd(
            transport,
            physics::SgsShearTransferWorkspace::requiredBytes(state.cells.size()),
            "source SGS shear staging");
      }
      peak = std::max(peak, transport);
    }
    return {peak, false,
            "source staging modeled; population growth, ID coordination and topology scratch require owner-local admission"};
  }

  void setAmrLevelSourceCouplingEnabled(bool enabled) noexcept override {
    m_amr_level_source_coupling_enabled = enabled;
  }

  void executeAmrLevelSmoothSources(
      core::StepContext& parent_context,
      std::span<const std::uint32_t> level_rows,
      const hydro::HydroUpdateContext& level_update,
      int level,
      core::SourceSplitPhase phase) override {
    if (!m_amr_level_source_coupling_enabled ||
        m_full_config.numerics.amr_time_stepping !=
            core::AmrTimeSteppingConfigMode::kRecursiveSubcycling) {
      throw std::logic_error(
          "recursive AMR smooth source execution was not enabled by composition");
    }
    if (phase == core::SourceSplitPhase::kUnsplitCompatibility) {
      throw std::invalid_argument(
          "recursive AMR smooth source execution requires a Strang half phase");
    }
    if (!(std::isfinite(level_update.dt_code) && level_update.dt_code > 0.0 &&
          std::isfinite(level_update.time_begin_code) &&
          std::isfinite(level_update.time_end_code) &&
          level_update.time_end_code > level_update.time_begin_code &&
          std::isfinite(level_update.scale_factor_begin) &&
          std::isfinite(level_update.scale_factor) &&
          std::isfinite(level_update.scale_factor_end) &&
          level_update.scale_factor_begin > 0.0 &&
          level_update.scale_factor > 0.0 &&
          level_update.scale_factor_end > 0.0)) {
      throw std::invalid_argument(
          "recursive AMR smooth source execution received an invalid level interval");
    }

    core::IntegratorState level_integrator = parent_context.integrator_state;
    level_integrator.current_time_code =
        level_update.time_begin_code + 0.5 * level_update.dt_code;
    level_integrator.current_scale_factor = level_update.scale_factor;
    level_integrator.current_redshift = m_is_cosmological
        ? std::max(0.0, 1.0 / level_update.scale_factor - 1.0)
        : 0.0;
    level_integrator.current_hubble_rate_code = level_update.hubble_rate_code;
    level_integrator.dt_time_code = level_update.dt_code;

    core::StepContext level_context{parent_context.state, level_integrator};
    level_context.active_set = core::ActiveSetDescriptor{
        .particle_indices = {},
        .cell_indices = level_rows,
        .particles_are_subset = true,
        .cells_are_subset = true,
        .particles_from_scheduler = false,
        .cells_from_scheduler = false,
        .has_generation_metadata = true,
        .source_particle_index_generation =
            parent_context.state.particleIndexGeneration(),
        .source_cell_index_generation = parent_context.state.cellIndexGeneration(),
        .source_scheduler_tick = parent_context.active_set.source_scheduler_tick,
        .has_global_synchronization_metadata = true,
        .globally_complete_active_set = false,
    };
    level_context.workspace = parent_context.workspace;
    level_context.cosmology_background = parent_context.cosmology_background;
    level_context.mode_policy = parent_context.mode_policy;
    level_context.profiler_session = parent_context.profiler_session;
    level_context.timeline_step = parent_context.timeline_step;
    level_context.timeline_step.time_begin_code = level_update.time_begin_code;
    level_context.timeline_step.time_end_code = level_update.time_end_code;
    level_context.timeline_step.dt_time_code = level_update.dt_code;
    level_context.timeline_step.scale_factor_begin = level_update.scale_factor_begin;
    level_context.timeline_step.scale_factor_midpoint = level_update.scale_factor;
    level_context.timeline_step.scale_factor_end = level_update.scale_factor_end;
    level_context.timeline_step.redshift_begin = m_is_cosmological
        ? std::max(0.0, 1.0 / level_update.scale_factor_begin - 1.0)
        : 0.0;
    level_context.timeline_step.redshift_end = m_is_cosmological
        ? std::max(0.0, 1.0 / level_update.scale_factor_end - 1.0)
        : 0.0;
    level_context.timeline_step.hubble_begin_code =
        level_update.hubble_rate_begin_code;
    level_context.timeline_step.hubble_end_code = level_update.hubble_rate_end_code;
    level_context.timeline_step.dt_time_si =
        m_is_cosmological && parent_context.cosmology_background != nullptr
        ? parent_context.cosmology_background->cosmicTimeIntervalSi(
              level_update.scale_factor_begin, level_update.scale_factor_end)
        : level_update.dt_code * m_units.timeSiPerCode();
    level_context.boundary = parent_context.boundary;
    level_context.stage = core::IntegrationStage::kSourceTerms;
    level_context.source_split_phase = phase;

    if (!m_source_interval.open || !m_source_interval.recursive_levels) {
      throw std::logic_error(
          "recursive AMR level source half invoked outside an open coarse source "
          "interval owned by AMR levels (composition/scheduler mismatch)");
    }
    if (phase == core::SourceSplitPhase::kStrangPreHalf) {
      if (m_source_interval.level_pre_halves ==
          std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("recursive AMR level pre-half counter overflow");
      }
      ++m_source_interval.level_pre_halves;
    } else {
      if (m_source_interval.level_post_halves >= m_source_interval.level_pre_halves) {
        throw std::logic_error(
            "recursive AMR level post-half without a matching pre-half");
      }
      ++m_source_interval.level_post_halves;
    }
    // The level half interval is [t_begin, t_mid] (pre) or [t_mid, t_end]
    // (post) of THIS level clock, evaluated at its physical midpoint. The
    // controller observation is accumulated; publication happens once at the
    // coarse post-half close (C03/C04: never per level, never per half).
    (void)physics::ensureGasCompositionInitialized(level_context.state, m_physics_config);
    const physics::SourceInterval interval = makeSourceInterval(
        level_context, sourceEvaluationScaleFactor(level_context));
    // Palindromic level split, identical to the global halves: pre = smooth
    // transport then thermochemistry, post = thermochemistry then transport.
    // Transport faces follow the causal owner rule documented at
    // restrictDiffusionFacesToLevelSelection: the pre half (t_begin, every
    // finer level synchronized) integrates same-level AND coarse/fine faces
    // owned by this level; the post half (finer levels still at t_begin)
    // integrates same-level faces only. Every rank enters.
    const double half_dt_code = 0.5 * level_update.dt_code;
    const double transport_scale_factor = phase == core::SourceSplitPhase::kStrangPreHalf
        ? 0.5 * (level_update.scale_factor_begin + level_update.scale_factor)
        : 0.5 * (level_update.scale_factor + level_update.scale_factor_end);
    const auto level_index = static_cast<std::size_t>(std::max(level, 0));
    if (phase == core::SourceSplitPhase::kStrangPreHalf) {
      if (m_source_interval.pending_interface_sync.size() <= level_index) {
        m_source_interval.pending_interface_sync.resize(level_index + 1U, 0U);
      }
      if (m_source_interval.pending_interface_sync[level_index] > 1U) {
        throw std::logic_error(
            "recursive AMR level pre half opened while the previous coarse/fine interface split of the level "
            "was never closed (missing post-children synchronization)");
      }
      executeDiffusionAndSgsShear(
          level_context, transport_scale_factor, half_dt_code,
          LevelFaceSelection{.level = static_cast<std::int32_t>(level), .same_level = true, .coarse_fine = true,
                             .phase_tag = levelPhaseTag(parent_context.integrator_state.step_index, level,
                                                        level_update.time_begin_code, 0U)});
      // T02 (patch 52): park the exact per-rank identity ledger of the opened
      // coarse/fine faces (and its lease) until the post-children closure.
      // Any local ledger failure is made collective before the count reduce.
      {
        std::exception_ptr ledger_failure = std::exchange(m_face_ledger_failure, nullptr);
        try {
          if (ledger_failure != nullptr) std::rethrow_exception(ledger_failure);
          if (m_source_interval.pending_interface_faces.size() <= level_index) {
            m_source_interval.pending_interface_faces.resize(level_index + 1U);
          }
          if (m_interface_face_ledger_leases.size() <= level_index) {
            m_interface_face_ledger_leases.resize(level_index + 1U);
          }
          if (!m_source_interval.pending_interface_faces[level_index].empty()) {
            throw std::logic_error("recursive AMR interface face ledger of the level was never consumed");
          }
        } catch (...) {
          ledger_failure = std::current_exception();
        }
        m_mpi_context.rethrowCollectivePreparationFailure(ledger_failure, "recursive AMR interface face ledger open");
        // noexcept hand-over (vector move + reservation move).
        m_source_interval.pending_interface_faces[level_index] = std::move(m_last_coarse_fine_faces);
        m_interface_face_ledger_leases[level_index] = std::move(m_last_coarse_fine_faces_lease);
        m_last_coarse_fine_faces = {};
      }
      // Global count of coarse/fine faces opened by this pre half; stored +1
      // so that "pre half seen, no coarse/fine faces" is distinguishable.
      const std::uint64_t opened = m_mpi_context.allreduceSumUint64(m_last_level_face_selection.coarse_fine);
      m_source_interval.pending_interface_sync[level_index] = opened + 1U;
    }
    executeThermochemistryInterval(
        level_context, interval, level_rows, false, level);
    if (phase == core::SourceSplitPhase::kStrangPostHalf) {
      executeDiffusionAndSgsShear(
          level_context, transport_scale_factor, half_dt_code,
          LevelFaceSelection{.level = static_cast<std::int32_t>(level), .same_level = true, .coarse_fine = false,
                             .phase_tag = levelPhaseTag(parent_context.integrator_state.step_index, level,
                                                        level_update.time_begin_code, 1U)});
    }
    if (m_source_interval.level_transport_halves == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("recursive AMR level transport counter overflow");
    }
    ++m_source_interval.level_transport_halves;
  }

  void executeAmrLevelInterfaceSynchronization(
      core::StepContext& parent_context,
      const hydro::HydroUpdateContext& level_update,
      int level) override {
    if (!m_amr_level_source_coupling_enabled || !m_source_interval.open || !m_source_interval.recursive_levels) {
      throw std::logic_error(
          "recursive AMR interface synchronization invoked outside an open level-owned source interval");
    }
    if (!(std::isfinite(level_update.dt_code) && level_update.dt_code > 0.0 &&
          level_update.scale_factor_end > 0.0 && level_update.scale_factor > 0.0)) {
      throw std::invalid_argument("recursive AMR interface synchronization received an invalid level interval");
    }
    const auto level_index = static_cast<std::size_t>(std::max(level, 0));
    if (level_index >= m_source_interval.pending_interface_sync.size() ||
        m_source_interval.pending_interface_sync[level_index] == 0U) {
      throw std::logic_error("recursive AMR interface synchronization without a matching level pre half");
    }
    // Level clock of the closing half [t_mid, t_end] (all levels >= L at t_end).
    core::IntegratorState level_integrator = parent_context.integrator_state;
    level_integrator.current_time_code = level_update.time_end_code;
    level_integrator.current_scale_factor = level_update.scale_factor_end;
    level_integrator.current_redshift = m_is_cosmological
        ? std::max(0.0, 1.0 / level_update.scale_factor_end - 1.0)
        : 0.0;
    level_integrator.dt_time_code = level_update.dt_code;
    core::StepContext level_context{parent_context.state, level_integrator};
    level_context.active_set = parent_context.active_set;
    level_context.workspace = parent_context.workspace;
    level_context.cosmology_background = parent_context.cosmology_background;
    level_context.mode_policy = parent_context.mode_policy;
    level_context.profiler_session = parent_context.profiler_session;
    level_context.timeline_step = parent_context.timeline_step;
    level_context.boundary = parent_context.boundary;
    level_context.stage = core::IntegrationStage::kSourceTerms;
    level_context.source_split_phase = core::SourceSplitPhase::kStrangPostHalf;
    executeDiffusionAndSgsShear(
        level_context, 0.5 * (level_update.scale_factor + level_update.scale_factor_end), 0.5 * level_update.dt_code,
        LevelFaceSelection{.level = static_cast<std::int32_t>(level), .same_level = false, .coarse_fine = true,
                           .phase_tag = levelPhaseTag(parent_context.integrator_state.step_index, level,
                                                      level_update.time_begin_code, 2U)});
    // T02 (patch 52): exact coverage. Every rank compares the sorted
    // canonical identities it integrated now with those its pre half opened
    // (identical set <=> no face missed, none duplicated, none foreign; a
    // missing face balanced by a duplicate elsewhere is detected, which the
    // global count alone cannot do). Mirrored remote faces are additionally
    // matched key-for-key across the two owners by the per-phase claim
    // handshake (DistributedElementDiffusionCoupling::setup, phase-tagged).
    {
      std::exception_ptr ledger_failure = std::exchange(m_face_ledger_failure, nullptr);
      try {
        if (ledger_failure != nullptr) std::rethrow_exception(ledger_failure);
        if (level_index >= m_source_interval.pending_interface_faces.size()) {
          throw std::logic_error("recursive AMR interface closure without an opened face ledger");
        }
        const auto& opened_faces = m_source_interval.pending_interface_faces[level_index];
        if (opened_faces.size() != m_last_coarse_fine_faces.size() ||
            !std::equal(opened_faces.begin(), opened_faces.end(), m_last_coarse_fine_faces.begin())) {
          const auto mismatch = std::mismatch(opened_faces.begin(), opened_faces.end(),
                                              m_last_coarse_fine_faces.begin(), m_last_coarse_fine_faces.end());
          const CanonicalInterfaceFace& face = mismatch.first != opened_faces.end()
              ? *mismatch.first
              : *mismatch.second;
          throw std::logic_error(
              "recursive AMR coarse/fine interface closure on level " + std::to_string(level) +
              " does not integrate exactly the faces its pre half opened (opened " +
              std::to_string(opened_faces.size()) + ", closing " + std::to_string(m_last_coarse_fine_faces.size()) +
              "; first differing face cells " + std::to_string(face.lo_id) + "/" + std::to_string(face.hi_id) +
              " axis " + std::to_string(face.axis) + ")");
        }
      } catch (...) {
        ledger_failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(ledger_failure, "recursive AMR interface face ledger close");
      // Consume: storage first, then its admission.
      std::vector<CanonicalInterfaceFace>().swap(m_source_interval.pending_interface_faces[level_index]);
      m_interface_face_ledger_leases[level_index].release();
      std::vector<CanonicalInterfaceFace>().swap(m_last_coarse_fine_faces);
      m_last_coarse_fine_faces_lease.release();
    }
    const std::uint64_t closed = m_mpi_context.allreduceSumUint64(m_last_level_face_selection.coarse_fine);
    if (closed + 1U != m_source_interval.pending_interface_sync[level_index]) {
      // The hierarchy is frozen inside one coarse interval (regrid only at
      // synchronized coarse boundaries), so the closing half must integrate
      // exactly the coarse/fine faces its pre half opened.
      throw std::logic_error(
          "recursive AMR coarse/fine interface synchronization integrated " + std::to_string(closed) +
          " faces but its pre half opened " +
          std::to_string(m_source_interval.pending_interface_sync[level_index] - 1U) +
          " (face missed or double integrated)");
    }
    m_source_interval.pending_interface_sync[level_index] = 1U;  // closed; pre half seen
    if (m_source_interval.level_interface_syncs == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("recursive AMR interface synchronization counter overflow");
    }
    ++m_source_interval.level_interface_syncs;
  }

  void execute(SourceMutationStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::sourceContext(
        view,
        {{RuntimeResourceKey::kSourceMutationState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleIdentity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleSpeciesIndex, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroPrimitiveState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kAmrPatchState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kEffectiveIsmThermodynamics, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead}});
    if (context.stage != core::IntegrationStage::kSourceTerms) {
      throw std::logic_error("source runtime received a non-source stage");
    }

    // T01: one explicit plan per dispatch (see SourceStagePlan). Smooth
    // operators advance only this dispatch's half interval; discrete events
    // are owned by the post half (or the compatibility dispatch) alone.
    const SourceStagePlan plan =
        makeSourceStagePlan(context, recursiveAmrLevelSourcesActive(context));
    if (plan.opens_interval) openSourceInterval(context, plan);
    requireOpenSourceInterval(context, plan);

    // Only the current source stage may borrow the scheduler. The hook itself
    // was installed once at owner construction and performs no allocation here.
    m_growth_scheduler = context.particle_scheduler;
    struct GrowthSchedulerReset {
      core::HierarchicalTimeBinScheduler*& slot;
      ~GrowthSchedulerReset() { slot = nullptr; }
    } growth_scheduler_reset{m_growth_scheduler};

    // Symmetric smooth split: the pre half applies transport then
    // thermochemistry, the post half thermochemistry then transport, so the
    // composition S_pre H S_post is palindromic in its smooth sub-operators.
    // Cooling is never also applied inside the hydro solver.
    if (plan.phase == core::SourceSplitPhase::kStrangPreHalf) {
      if (plan.smooth_transport) {
        executeDiffusionAndSgsShear(
            context, plan.smooth_evaluation_scale_factor, plan.smooth_dt_code);
      }
      if (plan.global_smooth_thermochemistry) executeThermochemistry(context, plan);
    } else {
      if (plan.global_smooth_thermochemistry) executeThermochemistry(context, plan);
      if (plan.smooth_transport) {
        executeDiffusionAndSgsShear(
            context, plan.smooth_evaluation_scale_factor, plan.smooth_dt_code);
      }
    }

    if (plan.discrete_events) {
      executeDiscreteSourceEvents(context, plan);
    }
    // The controller publishes exactly once per Strang interval, after every
    // half (and every AMR level half) contributed its observations.
    if (plan.closes_interval) closeSourceInterval(context);
    reportDataAssumptionsOnce(context);
  }

 private:
  // Discrete/stochastic source chronology (S01/S02/F01/F02/BH): executed once
  // per source interval at the absolute end-of-interval epoch with the full
  // interval duration. Ordering is fixed: star birth -> SSP/stellar returns,
  // winds/pre-SN and SN feedback -> BH seeding/accretion/feedback.
  void executeDiscreteSourceEvents(core::StepContext& context,
                                   const SourceStagePlan& plan) {
    if (!plan.discrete_events) {
      throw std::logic_error("discrete source events requested by a smooth-only plan");
    }
    // T02 discrete-event chronology (patch 51, specification-permitted
    // stochastic-event separation, route v10 T01/T02 + section 4.6): births,
    // SSP returns, winds/pre-SN, SN and BH events are NOT integrated by the
    // smooth multirate level clocks. Their physical event time is the
    // absolute end of the accepted coarse (or hierarchical-block) interval,
    // which is also their only legal publication epoch: every AMR level and
    // particle bin of the interval is synchronized there, no flux register or
    // temporal history is open, and new particles may be created without
    // violating KDK/regrid ownership. Stellar returns are integrated exactly
    // over absolute stellar age [age_last, age_now] (cumulative cursors,
    // gap/replay refusal), so the event budget is invariant under level
    // cadence, refinement, rank count and retry; unresolved deposits are
    // carried in the persistent enrichment-carry lanes. Exactly-once: one
    // dispatch per source interval (counter below, restored by rollback).
    std::exception_ptr epoch_failure;
    try {
      if (m_source_interval.discrete_event_dispatches != 0U) {
        throw std::logic_error("discrete source events dispatched twice inside one source interval");
      }
      if (plan.event_time_code != plan.coarse_time_end_code) {
        throw std::logic_error("discrete source event epoch is not the coarse interval end");
      }
      if (m_source_interval.recursive_levels) {
        if (!context.state.pending_flux_registers.empty() ||
            !context.state.amr_temporal_boundary_history.empty()) {
          throw std::logic_error(
              "discrete source events reached with open AMR flux registers or temporal histories (hierarchy not "
              "synchronized: illegal event publication epoch)");
        }
        for (const std::uint64_t pending : m_source_interval.pending_interface_sync) {
          if (pending > 1U) {
            throw std::logic_error("discrete source events reached with an open coarse/fine interface split");
          }
        }
      }
    } catch (...) {
      epoch_failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(epoch_failure, "discrete source event synchronization epoch");
    ++m_source_interval.discrete_event_dispatches;
    const double source_evaluation_scale_factor = plan.event_scale_factor;
    if (!std::isfinite(source_evaluation_scale_factor) ||
        source_evaluation_scale_factor <= 0.0) {
      throw std::runtime_error(
          "source runtime requires a finite positive discrete-event scale factor");
    }
    const std::size_t cell_count = context.state.cells.size();
    const bool local_has_mutable_gravity_sources =
        cell_count > 0U || context.state.star_particles.size() > 0U ||
        context.state.black_holes.size() > 0U;
    const std::uint64_t gravity_source_mutation_rank_count =
        m_mpi_context.allreduceSumUint64(local_has_mutable_gravity_sources ? 1ULL : 0ULL);
    bool star_formation_mutated_gas = false;
    std::vector<hydro::GasThermodynamicTouch> star_formation_touches;
    std::vector<std::uint32_t> star_formation_touched_rows;
    if (m_star_formation.config().enabled) {
      constexpr std::size_t k_star_formation_cell_batch_max = 4096U;
      const bool all_cells_active =
          !context.active_set.cells_are_subset &&
          context.active_set.cell_indices.empty();
      const std::span<const std::uint32_t> active_cells =
          context.active_set.cell_indices;
      const std::size_t local_active_count = all_cells_active
          ? cell_count : active_cells.size();
      const std::uint64_t local_batch_count = local_active_count == 0U
          ? 0U
          : 1U + static_cast<std::uint64_t>(
              (local_active_count - 1U) / k_star_formation_cell_batch_max);
      const std::uint64_t global_batch_count =
          m_mpi_context.allreduceMaxUint64(local_batch_count);
      const std::size_t required_input_capacity = std::min<std::size_t>(
          local_active_count, k_star_formation_cell_batch_max);

      core::MemoryReservation replacement_input_reservation;
      core::MemoryReservation replacement_cell_rows_reservation;
      std::exception_ptr reservation_failure;
      try {
        if (m_memory_governor != nullptr &&
            required_input_capacity > m_star_formation_inputs.capacity()) {
          const std::size_t input_bytes = core::checkedSizeMultiply(
              required_input_capacity,
              sizeof(physics::StarFormationCellInput),
              "star-formation bounded input batch");
          replacement_input_reservation = m_memory_governor->reserve(
              core::MemoryClass::kPhaseResident,
              static_cast<std::uint64_t>(input_bytes),
              "sources.star_formation.input_batch");
        }
        if (m_memory_governor != nullptr &&
            required_input_capacity > m_contiguous_cell_batch.capacity()) {
          replacement_cell_rows_reservation = m_memory_governor->reserve(
              core::MemoryClass::kPhaseResident,
              static_cast<std::uint64_t>(core::checkedSizeMultiply(
                  required_input_capacity, sizeof(std::uint32_t),
                  "source contiguous cell batch")),
              "sources.star_formation.contiguous_cell_batch");
        }
      } catch (...) {
        reservation_failure = std::current_exception();
      }
      if (m_runtime_services != nullptr) {
        FailureCoordinator(*m_runtime_services).rethrowCollectiveFailure(
            reservation_failure,
            "star-formation memory preflight");
      } else if (reservation_failure != nullptr) {
        std::rethrow_exception(reservation_failure);
      }
      if (replacement_input_reservation.pending()) {
        replacement_input_reservation.commit();
      }
      if (required_input_capacity > m_star_formation_inputs.capacity()) {
        replaceSourceVectorCapacity(m_star_formation_inputs, required_input_capacity,
                                    "source input");
        if (m_memory_governor != nullptr) {
          m_star_formation_input_reservation.release();
          m_star_formation_input_reservation =
              std::move(replacement_input_reservation);
        }
      }
      if (replacement_cell_rows_reservation.pending()) {
        replacement_cell_rows_reservation.commit();
      }
      if (required_input_capacity > m_contiguous_cell_batch.capacity()) {
        // The replacement coexists with the old retained allocation. Preserve
        // the old lease if allocation fails and reconcile only after success.
        replaceSourceVectorCapacity(m_contiguous_cell_batch, required_input_capacity,
                                    "source cell-row");
        if (m_memory_governor != nullptr) {
          m_contiguous_cell_batch_reservation.release();
          m_contiguous_cell_batch_reservation =
              std::move(replacement_cell_rows_reservation);
        }
      }

      std::unique_ptr<core::GovernedScratchArena> parent_arena;
      std::unique_ptr<std::pmr::vector<std::uint64_t>> parent_index;

      std::exception_ptr parent_failure;
      try {
        parent_arena = std::make_unique<core::GovernedScratchArena>(
            m_memory_governor, core::MemoryClass::kScratchArena,
            sourceParentIndexBytes(context.state.patches.size()),
            "sources.star_formation.parent_index");
        parent_index = std::make_unique<std::pmr::vector<std::uint64_t>>(
            buildSourceParentIndex(context.state, parent_arena->resource()));
      } catch (...) { parent_failure = std::current_exception(); }
      m_mpi_context.rethrowCollectivePreparationFailure(
          parent_failure, "source parent metadata preparation");
      std::exception_ptr neighbor_failure;
      try {
        m_sf_gradient_unavailable_cells = 0U;
        m_sf_gradient_coarse_fine_cells = 0U;
        buildOwnedLeafCellMetadata(context);
        m_sf_neighbor_index = std::make_unique<internal::StarFormationNeighborIndex>(
            context.state, m_hydro_boundary == core::BoundaryCondition::kPeriodic, m_box_extent_comoving,
            m_owned_leaf_mask);
      } catch (...) { neighbor_failure = std::current_exception(); }
      m_mpi_context.rethrowCollectivePreparationFailure(
          neighbor_failure, "star-formation neighbour index preparation");
      if (m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1) {
        // CHUI Part 3 Checkpoint C: remote neighbour halo (collective).
        importRemoteStarFormationNeighbors(context);
      }

      const std::size_t particle_count_before_birth =
          context.state.particles.size();
      std::uint64_t local_spawned_particles = 0U;
      // Star formation removes gas mass/density from its host cells.
      const auto mark_star_formation_gas_mutation = [&]() {
        if (local_spawned_particles > 0U) {
          star_formation_mutated_gas = true;
        }
      };
      for (std::uint64_t batch_round = 0U;
           batch_round < global_batch_count; ++batch_round) {
        const std::size_t batch_begin = batch_round < local_batch_count
            ? core::checkedSizeMultiply(
                core::checkedIntegralNarrow<std::size_t>(
                    batch_round,
                    "star-formation batch round"),
                k_star_formation_cell_batch_max,
                "star-formation batch begin")
            : local_active_count;
        const std::size_t batch_size = batch_begin < local_active_count
            ? std::min<std::size_t>(
                k_star_formation_cell_batch_max,
                local_active_count - batch_begin)
            : 0U;
        std::span<const std::uint32_t> cell_batch;
        if (batch_size == 0U) {
          m_contiguous_cell_batch.clear();
          cell_batch = {};
        } else if (all_cells_active) {
          m_contiguous_cell_batch.resize(batch_size);
          for (std::size_t i = 0U; i < batch_size; ++i) {
            m_contiguous_cell_batch[i] =
                core::checkedIntegralNarrow<std::uint32_t>(
                    batch_begin + i,
                    "star-formation contiguous active gas row");
          }
          cell_batch = m_contiguous_cell_batch;
        } else {
          cell_batch = active_cells.subspan(batch_begin, batch_size);
        }
        buildStarFormationInputs(
            context,
            cell_batch,
            source_evaluation_scale_factor,
            *parent_index);
        // Bound source-report and birth-plan metadata separately from the ID
        // registry's own precommit/result leases. The two report arrays are
        // exact-reserved by applyFromInputs before any birth mutation. The
        // arena owns only plans and immutable birth keys; no heap fallback.
        std::unique_ptr<core::GovernedScratchArena> metadata_arena;
        core::MemoryReservation report_reservation;
        std::exception_ptr metadata_failure;
        try {
          const std::size_t max_births = core::checkedSizeMultiply(
              batch_size,
              static_cast<std::size_t>(m_star_formation.config().max_spawn_particles_per_cell_step),
              "source metadata maximum birth count");
          const std::uint64_t plan_bytes = static_cast<std::uint64_t>(
              core::checkedSizeMultiply(batch_size,
                  physics::starFormationBirthPlanBytes(), "source metadata birth plans"));
          const std::uint64_t key_bytes = static_cast<std::uint64_t>(
              core::checkedSizeMultiply(max_births, sizeof(std::uint64_t),
                  "source metadata birth keys"));
          const std::uint64_t report_bytes = static_cast<std::uint64_t>(
              core::checkedSizeMultiply(max_births,
                  sizeof(std::uint32_t) + sizeof(std::uint64_t),
                  "source metadata report rows"));
          metadata_arena = std::make_unique<core::GovernedScratchArena>(
              m_memory_governor, core::MemoryClass::kScratchArena,
              core::checkedMemoryBytesAdd(
                  core::checkedMemoryBytesAdd(plan_bytes, key_bytes,
                      "source metadata arena"), 256U, "source metadata alignment"),
              "sources.star_formation.birth_metadata");
          if (m_memory_governor != nullptr && report_bytes != 0U) {
            report_reservation = m_memory_governor->reserve(
                core::MemoryClass::kDiagnostic, report_bytes,
                "sources.star_formation.report_metadata");
          }
        } catch (...) { metadata_failure = std::current_exception(); }
        m_mpi_context.rethrowCollectivePreparationFailure(
            metadata_failure, "source birth/report metadata admission");
        if (report_reservation.pending()) report_reservation.commit();
        m_star_formation.setBirthTimeCode(plan.event_time_code);
        const physics::StarFormationStepReport batch_report =
            m_star_formation.applyFromInputs(
                context.state,
                m_star_formation_inputs,
                plan.event_dt_code,
                source_evaluation_scale_factor,
                context.integrator_state.step_index,
                &m_particle_id_registry,
                metadata_arena->resource());
        if (batch_report.counters.spawned_particles >
            std::numeric_limits<std::uint64_t>::max() - local_spawned_particles) {
          throw std::overflow_error(
              "source runtime star-formation batch spawn count overflows uint64");
        }
        local_spawned_particles += batch_report.counters.spawned_particles;
        // Host rows whose gas mass/density the batch reduced (F4 refresh).
        // Bounded by the number of births, which is itself admitted above.
        star_formation_touched_rows.insert(
            star_formation_touched_rows.end(),
            batch_report.spawned_from_cells.begin(),
            batch_report.spawned_from_cells.end());
      }

      if (local_spawned_particles > 0U) {
        const std::size_t particle_count_after_birth =
            context.state.particles.size();
        if (particle_count_after_birth < particle_count_before_birth ||
            particle_count_after_birth - particle_count_before_birth !=
                static_cast<std::size_t>(local_spawned_particles)) {
          throw std::runtime_error(
              "source runtime star-formation report disagrees with appended particle rows");
        }
        for (std::size_t particle_index = particle_count_before_birth;
             particle_index < particle_count_after_birth; ++particle_index) {
          context.state.particle_sidecar.last_drift_time_code[particle_index] =
              context.timeline_step.time_end_code;
          context.state.particle_sidecar.last_drift_scale_factor[particle_index] =
              source_evaluation_scale_factor;
        }
        internal::synchronizeParentParticleCompatibilityMirrors(
            context.state, m_world_rank, "SourceRuntime star-formation batch");
        if (context.newly_created_particle_ids != nullptr) {
          context.newly_created_particle_ids->insert(
              context.newly_created_particle_ids->end(),
              context.state.particle_sidecar.particle_id.begin() +
                  static_cast<std::ptrdiff_t>(particle_count_before_birth),
              context.state.particle_sidecar.particle_id.end());
        }
      }
      mark_star_formation_gas_mutation();
      if (local_spawned_particles > 0U &&
          context.particle_scheduler != nullptr) {
        const std::uint64_t current_tick =
            context.particle_scheduler->currentTick();
        if (current_tick == std::numeric_limits<std::uint64_t>::max()) {
          throw std::overflow_error(
              "source runtime newborn activation tick overflows uint64");
        }
        context.particle_scheduler->appendElements(
            core::checkedIntegralNarrow<std::uint32_t>(
                local_spawned_particles,
                "source runtime newborn scheduler append count"),
            0U,
            current_tick + 1U);
      }
    }

    m_sf_neighbor_index.reset();  // stage-local; release before the next phase
    // One owner per stellar source: the Part 2 SSP/event/mechanical-feedback
    // path replaces (never adds to) the legacy table evolution + deposition.
    if (m_ssp_model != nullptr) {
      executePart2StellarPopulation(context, source_evaluation_scale_factor);
    } else {
      executeStellarEvolutionAndEnrichment(context, source_evaluation_scale_factor);
    }
    physics::BlackHoleAgnStepReport black_hole_report{};
    if (m_bh_enabled) {
      black_hole_report = m_black_hole.apply(
          context.state,
          {},
          plan.event_dt_code,
          source_evaluation_scale_factor,
          m_coordinate_frame == core::CoordinateFrame::kComoving,
          context.integrator_state.step_index,
          &m_particle_id_registry);
    }
    if (black_hole_report.counters.gas_mass_removed_code > 0.0) {
      // Generic gas particles are compatibility mirrors only. Keep them coherent
      // for I/O/legacy consumers after the authoritative gas->BH transaction.
      internal::synchronizeParentParticleCompatibilityMirrors(
          context.state, m_world_rank, "SourceRuntime BH accretion/seeding batch");
    }
    if (!black_hole_report.seeded_particle_ids.empty()) {
      if (context.newly_created_particle_ids != nullptr) {
        context.newly_created_particle_ids->insert(
            context.newly_created_particle_ids->end(),
            black_hole_report.seeded_particle_ids.begin(),
            black_hole_report.seeded_particle_ids.end());
      }
      if (context.particle_scheduler != nullptr) {
        const std::uint64_t current_tick = context.particle_scheduler->currentTick();
        if (current_tick == std::numeric_limits<std::uint64_t>::max()) {
          throw std::overflow_error("source runtime BH activation tick overflows uint64");
        }
        context.particle_scheduler->appendElements(
            static_cast<std::uint32_t>(black_hole_report.seeded_particle_ids.size()),
            0U,
            current_tick + 1U);
      }
    }
    if (star_formation_mutated_gas || black_hole_report.counters.gas_mass_removed_code > 0.0) {
      // Audit F4 (density-changing sources): star formation and BH accretion
      // change gas density, which changes the effective-EOS pressure and
      // signal speed (and the ideal-gas pressure cache). Star formation
      // reports its host rows, so only those are refreshed; the BH model does
      // not expose its accretion cells, so a BH gas removal refreshes the
      // whole rank (O(N_gas), only on steps with accretion).
      context.state.notifyGasThermodynamicMutation();
      std::unique_ptr<physics::EffectiveIsmThermodynamicClosure> closure;
      if (m_effective_eos_table != nullptr) {
        closure = std::make_unique<physics::EffectiveIsmThermodynamicClosure>(m_effective_eos_table);
      }
      const hydro::GasThermodynamicRefreshPolicy policy{
          .adiabatic_index = 5.0 / 3.0,
          .density_floor_code = m_hydro_density_floor_code,
          .pressure_floor_code = m_hydro_pressure_floor_code,
          .scale_factor = source_evaluation_scale_factor,
          .redshift = std::max(0.0, 1.0 / source_evaluation_scale_factor - 1.0),
          .closure = closure.get(),
          .temperature_k_per_internal_energy_code =
              hydro::temperatureKelvinPerInternalEnergyCode(5.0 / 3.0, m_units.velocity_si_per_code)};
      if (black_hole_report.counters.gas_mass_removed_code > 0.0) {
        (void)hydro::refreshAllGasThermodynamics(context.state, policy);
      } else {
        std::sort(star_formation_touched_rows.begin(), star_formation_touched_rows.end());
        star_formation_touched_rows.erase(
            std::unique(star_formation_touched_rows.begin(), star_formation_touched_rows.end()),
            star_formation_touched_rows.end());
        star_formation_touches.reserve(star_formation_touched_rows.size());
        for (const std::uint32_t row : star_formation_touched_rows) {
          // Star formation does not change specific internal energy.
          star_formation_touches.push_back(hydro::GasThermodynamicTouch{
              .cell_row = row,
              .internal_energy_before_code = context.state.gas_cells.internal_energy_code[row]});
        }
        (void)hydro::refreshGasThermodynamics(context.state, star_formation_touches, policy);
      }
    }
    if (gravity_source_mutation_rank_count > 0U) {
      // Source models can change gas/star/BH gravitating mass and membership.
      // Advance once per globally ordered source stage on every rank so PM/tree
      // validity follows authoritative state rather than PM refresh sequencing.
      context.state.bumpGravitySourceGeneration();
    }
  }

  [[nodiscard]] bool recursiveAmrLevelSourcesActive(
      const core::StepContext& context) const {
    if (!m_amr_level_source_coupling_enabled ||
        m_full_config.numerics.amr_time_stepping !=
            core::AmrTimeSteppingConfigMode::kRecursiveSubcycling ||
        context.source_split_phase ==
            core::SourceSplitPhase::kUnsplitCompatibility) {
      return false;
    }
    const std::uint64_t local_has_production_amr =
        amr::hasProductionAmrHydroCoverage(context.state) ? 1ULL : 0ULL;
    return m_mpi_context.allreduceSumUint64(local_has_production_amr) > 0U;
  }

  [[nodiscard]] SourceStagePlan makeSourceStagePlan(
      const core::StepContext& context, bool recursive_levels_active) const {
    SourceStagePlan plan;
    plan.phase = context.source_split_phase;
    plan.coarse_time_begin_code = context.timeline_step.time_begin_code;
    plan.coarse_time_end_code = context.timeline_step.time_end_code;
    plan.coarse_dt_code = plan.coarse_time_end_code - plan.coarse_time_begin_code;
    if (!(std::isfinite(plan.coarse_dt_code) && plan.coarse_dt_code > 0.0)) {
      // Legacy direct callers may omit explicit endpoints; the integrator's
      // accepted dt is then the only interval authority.
      plan.coarse_time_begin_code = context.integrator_state.current_time_code;
      plan.coarse_dt_code = context.integrator_state.dt_time_code;
      plan.coarse_time_end_code = plan.coarse_time_begin_code + plan.coarse_dt_code;
    }
    if (!(std::isfinite(plan.coarse_dt_code) && plan.coarse_dt_code > 0.0)) {
      throw std::invalid_argument("source stage requires a finite positive interval");
    }
    switch (plan.phase) {
      case core::SourceSplitPhase::kStrangPreHalf:
        plan.smooth_fraction_begin = 0.0;
        plan.smooth_fraction_end = 0.5;
        plan.opens_interval = true;
        break;
      case core::SourceSplitPhase::kStrangPostHalf:
        plan.smooth_fraction_begin = 0.5;
        plan.smooth_fraction_end = 1.0;
        plan.discrete_events = true;
        plan.closes_interval = true;
        break;
      case core::SourceSplitPhase::kUnsplitCompatibility:
        plan.smooth_fraction_begin = 0.0;
        plan.smooth_fraction_end = 1.0;
        plan.discrete_events = true;
        plan.opens_interval = true;
        plan.closes_interval = true;
        break;
    }
    plan.smooth_time_begin_code =
        plan.coarse_time_begin_code + plan.smooth_fraction_begin * plan.coarse_dt_code;
    plan.smooth_time_end_code =
        plan.coarse_time_begin_code + plan.smooth_fraction_end * plan.coarse_dt_code;
    plan.smooth_dt_code = plan.smooth_time_end_code - plan.smooth_time_begin_code;
    plan.smooth_scale_factor_begin =
        scaleFactorAtPhysicalFraction(context, plan.smooth_fraction_begin);
    plan.smooth_scale_factor_end =
        scaleFactorAtPhysicalFraction(context, plan.smooth_fraction_end);
    plan.smooth_evaluation_scale_factor = sourceEvaluationScaleFactor(context);
    plan.event_time_code = plan.coarse_time_end_code;
    plan.event_scale_factor = m_is_cosmological
        ? context.timeline_step.scale_factor_end
        : 1.0;
    plan.event_dt_code = context.integrator_state.dt_time_code;
    if (!(std::isfinite(plan.event_dt_code) && plan.event_dt_code > 0.0)) {
      plan.event_dt_code = plan.coarse_dt_code;
    }
    plan.global_smooth_thermochemistry = !recursive_levels_active;
    // T02: with recursive AMR levels the element diffusion and SGS shear are
    // integrated on the level clocks (executeAmrLevelSmoothSources), never
    // additionally on the coarse interval.
    plan.smooth_transport = !recursive_levels_active;
    if (!(std::isfinite(plan.smooth_dt_code) && plan.smooth_dt_code > 0.0 &&
          std::isfinite(plan.smooth_evaluation_scale_factor) &&
          plan.smooth_evaluation_scale_factor > 0.0)) {
      throw std::runtime_error(
          "source stage plan has an invalid smooth half interval or evaluation epoch");
    }
    return plan;
  }

  void openSourceInterval(const core::StepContext& context,
                          const SourceStagePlan& plan) {
    if (m_source_interval.open) {
      // A pre half after an unclosed interval means a post half was skipped:
      // the previous interval's events and controller update never ran.
      throw std::logic_error(
          "source interval opened while the previous Strang interval is still open");
    }
    if (m_source_interval.sequence == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("source interval sequence overflow");
    }
    const std::uint64_t sequence = m_source_interval.sequence + 1U;
    m_source_interval = SourceIntervalAccumulator{};
    resetInterfaceFaceLedger();
    m_source_interval.open = true;
    m_source_interval.sequence = sequence;
    m_source_interval.step_index = context.integrator_state.step_index;
    m_source_interval.recursive_levels = !plan.global_smooth_thermochemistry;
  }

  void requireOpenSourceInterval(const core::StepContext& context,
                                 const SourceStagePlan& plan) const {
    if (!m_source_interval.open ||
        m_source_interval.step_index != context.integrator_state.step_index) {
      throw std::logic_error(
          "source half dispatched outside its Strang interval (missing pre half)");
    }
    if (m_source_interval.recursive_levels != !plan.global_smooth_thermochemistry) {
      throw std::logic_error(
          "AMR-level smooth-source ownership changed inside one Strang interval");
    }
  }

  void accumulateSourceIntervalObservation(
      const physics::ChemistryShadowReport& shadow,
      const physics::ChemistryAdaptiveReport& adaptive) noexcept {
    const auto add = [](std::uint64_t& value, std::uint64_t increment) {
      value = increment > std::numeric_limits<std::uint64_t>::max() - value
          ? std::numeric_limits<std::uint64_t>::max()
          : value + increment;
    };
    physics::ChemistryShadowReport& total = m_source_interval.shadow;
    const std::uint64_t previous_samples = total.samples;
    const long double weighted_cost =
        static_cast<long double>(total.mean_substep_cost_ratio) *
            static_cast<long double>(previous_samples) +
        static_cast<long double>(shadow.mean_substep_cost_ratio) *
            static_cast<long double>(shadow.samples);
    add(total.samples, shadow.samples);
    total.max_relative_species_difference = std::max(
        total.max_relative_species_difference, shadow.max_relative_species_difference);
    total.max_relative_electron_difference = std::max(
        total.max_relative_electron_difference, shadow.max_relative_electron_difference);
    total.max_relative_energy_difference = std::max(
        total.max_relative_energy_difference, shadow.max_relative_energy_difference);
    add(total.ion_samples, shadow.ion_samples);
    total.max_ion_stage_fraction_difference = std::max(
        total.max_ion_stage_fraction_difference, shadow.max_ion_stage_fraction_difference);
    total.max_relative_ion_electron_difference = std::max(
        total.max_relative_ion_electron_difference, shadow.max_relative_ion_electron_difference);
    total.max_relative_ion_temperature_difference = std::max(
        total.max_relative_ion_temperature_difference, shadow.max_relative_ion_temperature_difference);
    total.max_relative_ion_potential_difference = std::max(
        total.max_relative_ion_potential_difference, shadow.max_relative_ion_potential_difference);
    total.mean_substep_cost_ratio = total.samples == 0U
        ? 0.0
        : static_cast<double>(weighted_cost / static_cast<long double>(total.samples));
    physics::ChemistryAdaptiveReport& a = m_source_interval.adaptive;
    add(a.promoted, adaptive.promoted);
    add(a.demoted, adaptive.demoted);
    add(a.held_by_dwell, adaptive.held_by_dwell);
    add(a.promoted_by_cooling_time, adaptive.promoted_by_cooling_time);
    add(a.promoted_by_error_history, adaptive.promoted_by_error_history);
    // Cell-population counts are state snapshots, not increments: keep the
    // latest committed population rather than summing halves/levels.
    a.detailed_cells = adaptive.detailed_cells;
    a.hybrid_cells = adaptive.hybrid_cells;
  }

  // Exactly one collective controller publication per Strang interval. All
  // ranks reach this point together (the post-half dispatch is rank-uniform).
  void closeSourceInterval(core::StepContext& context) {
    if (!m_source_interval.open) {
      throw std::logic_error("source interval closed without being opened");
    }
    if (m_source_interval.recursive_levels &&
        m_source_interval.level_pre_halves != m_source_interval.level_post_halves) {
      throw std::logic_error(
          "recursive AMR level sources left an unmatched pre half at interval close");
    }
    if (m_source_interval.recursive_levels &&
        m_source_interval.level_transport_halves !=
            m_source_interval.level_pre_halves + m_source_interval.level_post_halves) {
      // T02: every level half must have integrated its level-clock transport
      // faces exactly once (no face clock skipped or doubled).
      throw std::logic_error(
          "recursive AMR level transport halves do not match the level source halves");
    }
    for (const std::uint64_t pending : m_source_interval.pending_interface_sync) {
      if (pending > 1U) {
        throw std::logic_error(
            "recursive AMR coarse/fine interface split left open at interval close (faces integrated for dt/2 only)");
      }
    }
    for (const auto& ledger : m_source_interval.pending_interface_faces) {
      if (!ledger.empty()) {
        throw std::logic_error("recursive AMR interface face ledger left open at interval close");
      }
    }
    resetInterfaceFaceLedger();
    m_last_shadow_report = m_source_interval.shadow;
    updateGlobalChemistryControl(
        context, m_source_interval.shadow, m_source_interval.adaptive);
    m_source_interval.open = false;
    m_source_interval.shadow = {};
    m_source_interval.adaptive = {};
  }

  [[nodiscard]] double scaleFactorAtPhysicalFraction(
      const core::StepContext& context, double fraction) const {
    if (!m_is_cosmological) return 1.0;
    if (!(fraction >= 0.0 && fraction <= 1.0)) {
      throw std::invalid_argument("source physical-time fraction must lie in [0,1]");
    }
    const double begin = context.timeline_step.scale_factor_begin;
    const double end = context.timeline_step.scale_factor_end;
    if (!(std::isfinite(begin) && std::isfinite(end) && begin > 0.0 && end >= begin)) {
      throw std::runtime_error("source interval has invalid cosmological scale-factor bounds");
    }
    if (fraction == 0.0 || begin == end) return begin;
    if (fraction == 1.0) return end;
    if (context.cosmology_background == nullptr) {
      return begin + fraction * (end - begin);
    }
    const double total_seconds =
        context.cosmology_background->cosmicTimeIntervalSi(begin, end);
    if (!(std::isfinite(total_seconds) && total_seconds >= 0.0)) {
      throw std::runtime_error("source interval has invalid cosmological elapsed time");
    }
    const double target_seconds = fraction * total_seconds;
    double lower = begin;
    double upper = end;
    // Fixed-iteration inversion of the same background integral used by the
    // source interval keeps the physical-time split deterministic on all ranks.
    for (std::uint32_t iteration = 0U; iteration < 64U; ++iteration) {
      const double midpoint = lower + 0.5 * (upper - lower);
      const double elapsed =
          context.cosmology_background->cosmicTimeIntervalSi(begin, midpoint);
      if (elapsed < target_seconds) {
        lower = midpoint;
      } else {
        upper = midpoint;
      }
    }
    return lower + 0.5 * (upper - lower);
  }

  [[nodiscard]] double sourceEvaluationScaleFactor(
      const core::StepContext& context) const {
    switch (context.source_split_phase) {
      case core::SourceSplitPhase::kStrangPreHalf:
        return scaleFactorAtPhysicalFraction(context, 0.25);
      case core::SourceSplitPhase::kStrangPostHalf:
        return scaleFactorAtPhysicalFraction(context, 0.75);
      case core::SourceSplitPhase::kUnsplitCompatibility:
        return m_is_cosmological ? context.timeline_step.scale_factor_end : 1.0;
    }
    throw std::logic_error("unknown source split phase");
  }

  [[nodiscard]] physics::SourceInterval makeSourceInterval(
      const core::StepContext& context, double source_evaluation_scale_factor) const {
    physics::SourceInterval interval;
    const bool pre_half =
        context.source_split_phase == core::SourceSplitPhase::kStrangPreHalf;
    const bool post_half =
        context.source_split_phase == core::SourceSplitPhase::kStrangPostHalf;
    const double fraction_begin = post_half ? 0.5 : 0.0;
    const double fraction_end = pre_half ? 0.5 : 1.0;
    const double code_begin = context.timeline_step.time_begin_code;
    const double code_end = context.timeline_step.time_end_code;
    interval.time_begin_code = code_begin + fraction_begin * (code_end - code_begin);
    interval.time_end_code = code_begin + fraction_end * (code_end - code_begin);
    interval.scale_factor_begin = scaleFactorAtPhysicalFraction(context, fraction_begin);
    interval.scale_factor_end = scaleFactorAtPhysicalFraction(context, fraction_end);
    if (m_is_cosmological && context.cosmology_background != nullptr) {
      interval.physical_seconds = context.cosmology_background->cosmicTimeIntervalSi(
          interval.scale_factor_begin, interval.scale_factor_end);
    } else {
      interval.physical_seconds =
          (interval.time_end_code - interval.time_begin_code) * m_units.timeSiPerCode();
      if (!(interval.physical_seconds > 0.0)) {
        const double full_seconds =
            elapsedStellarEvolutionYears(context) * physics::k_seconds_per_year;
        interval.physical_seconds = (pre_half || post_half)
            ? 0.5 * full_seconds
            : full_seconds;
      }
    }
    interval.evaluation_scale_factor = source_evaluation_scale_factor;
    interval.evaluation_redshift = m_is_cosmological
        ? std::max(0.0, 1.0 / source_evaluation_scale_factor - 1.0)
        : 0.0;
    interval.cosmological = m_is_cosmological;
    interval.accepted_step_index = context.integrator_state.step_index;
    return interval;
  }

  [[nodiscard]] physics::SourceInterval makeSourceInterval(
      const core::StepContext& context,
      const hydro::HydroUpdateContext& update) const {
    physics::SourceInterval interval;
    interval.time_begin_code = update.time_begin_code;
    interval.time_end_code = update.time_end_code;
    if (!(interval.time_end_code > interval.time_begin_code)) {
      interval.time_begin_code = context.timeline_step.time_begin_code;
      interval.time_end_code = interval.time_begin_code + update.dt_code;
    }
    interval.scale_factor_begin = m_is_cosmological
        ? std::max(update.scale_factor_begin, 1.0e-12) : 1.0;
    interval.scale_factor_end = m_is_cosmological
        ? std::max(update.scale_factor_end, interval.scale_factor_begin) : 1.0;
    double elapsed_si = update.dt_code * m_units.timeSiPerCode();
    if (m_is_cosmological && context.cosmology_background != nullptr &&
        interval.scale_factor_end >= interval.scale_factor_begin) {
      elapsed_si = context.cosmology_background->cosmicTimeIntervalSi(
          interval.scale_factor_begin, interval.scale_factor_end);
    }
    interval.physical_seconds = std::max(elapsed_si, 0.0);
    interval.evaluation_scale_factor = m_is_cosmological
        ? std::max(update.scale_factor_end, 1.0e-12) : 1.0;
    interval.evaluation_redshift = m_is_cosmological
        ? std::max(0.0, 1.0 / interval.evaluation_scale_factor - 1.0) : 0.0;
    interval.cosmological = m_is_cosmological;
    interval.accepted_step_index = context.integrator_state.step_index;
    return interval;
  }

  // CHUI Part 2 thermochemistry stage. Staged inside a StageTransaction with
  // one collective vote: a cell that cannot integrate the full interval
  // rejects the stage on EVERY rank (kNumericalRetry) and nothing is
  // committed (audit F2 contract). Touched rows are refreshed before any CFL
  // or source consumer reads P, T or c_s (audit F4 contract).
  void executeThermochemistry(core::StepContext& context, const SourceStagePlan& plan) {
    (void)physics::ensureGasCompositionInitialized(context.state, m_physics_config);
    if (m_thermochemistry == nullptr) return;
    const physics::SourceInterval interval = makeSourceInterval(context, plan);
    const bool all_cells_active =
        !context.active_set.cells_are_subset && context.active_set.cell_indices.empty();
    const std::span<const std::uint32_t> rows = all_cells_active
        ? std::span<const std::uint32_t>{}
        : std::span<const std::uint32_t>(context.active_set.cell_indices);
    if (m_source_interval.global_smooth_halves ==
        std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("global smooth source half counter overflow");
    }
    ++m_source_interval.global_smooth_halves;
    executeThermochemistryInterval(context, interval, rows, all_cells_active, -1);
  }

  // Builds the physical sub-interval of one global smooth dispatch from the
  // T01 plan (half interval for Strang phases, full for compatibility).
  [[nodiscard]] physics::SourceInterval makeSourceInterval(
      const core::StepContext& context, const SourceStagePlan& plan) const {
    physics::SourceInterval interval;
    interval.time_begin_code = plan.smooth_time_begin_code;
    interval.time_end_code = plan.smooth_time_end_code;
    interval.scale_factor_begin = plan.smooth_scale_factor_begin;
    interval.scale_factor_end = plan.smooth_scale_factor_end;
    if (m_is_cosmological && context.cosmology_background != nullptr) {
      interval.physical_seconds = context.cosmology_background->cosmicTimeIntervalSi(
          interval.scale_factor_begin, interval.scale_factor_end);
    } else {
      interval.physical_seconds = plan.smooth_dt_code * m_units.timeSiPerCode();
    }
    interval.evaluation_scale_factor = plan.smooth_evaluation_scale_factor;
    interval.evaluation_redshift = m_is_cosmological
        ? std::max(0.0, 1.0 / plan.smooth_evaluation_scale_factor - 1.0)
        : 0.0;
    interval.cosmological = m_is_cosmological;
    interval.accepted_step_index = context.integrator_state.step_index;
    return interval;
  }

  void executeThermochemistryInterval(
      core::StepContext& context, const physics::SourceInterval& interval,
      std::span<const std::uint32_t> rows, bool all_cells_active,
      int recursive_level) {
    if (m_thermochemistry == nullptr) return;
    const std::size_t row_count =
        all_cells_active ? context.state.cells.size() : rows.size();
    // Local source-stage atomicity (corrective gate): regime decisions,
    // fidelity tags, dwell counters, chemistry, dust and SGS are evaluated
    // from the step-start state, staged per row, and written together only
    // after the collective accept. A rejection leaves state, metadata and
    // the published diagnostics exactly as they were.
    // G0.4: the committed promotion bias has exactly one authority, the
    // canonical sidecar "chui_adaptive_chemistry_control_v1" carried by
    // SimulationState (restart payload and rollback copy). The operator's
    // cached bias is re-derived from it before EVERY fidelity decision, so a
    // rejected/rolled-back step, a restart, or a retry cannot reuse a value
    // that was never committed to the canonical state.
    loadCommittedChemistryControl(context.state);
    physics::BaryonicSourceStage stage(*m_thermochemistry, m_source_stage_config);
    const internal::MpiStageCollectiveDecision decision(m_mpi_context);
    core::StageTransaction transaction("sources.baryonic_stage", m_memory_governor, &decision);
    (void)transaction.reserve(core::MemoryClass::kPhaseResident,
                              physics::BaryonicSourceStage::stagingBytes(
                                  row_count, context.state.gas_composition.laneCount()),
                              "sources.baryonic_stage.staged_rows");
    const bool local_work = row_count > 0U && interval.valid();
    if (local_work) {
      (void)transaction.prepare([&] {
        stage.prepare(context.state, interval, rows);
        auto& composition = context.state.gas_composition;
        core::RetainedCapacityTransaction sparse_capacity(
            [&composition] { return composition.ownedCapacityBytes(); });
        composition.planSparseRowCapacity(
            sparse_capacity,
            core::checkedSizeAdd(composition.sparseActiveRowCount(),
                                 stage.sparseRowsToActivate(),
                                 "adaptive chemistry sparse row capacity"));
        sparse_capacity.execute(
            m_memory_governor, core::MemoryClass::kCanonicalPersistent,
            "sources.baryonic_stage.adaptive_chemistry_rows");
      });
    } else {
      (void)transaction.prepare([] {});
    }
    if (!transaction.validate(true, core::StageRejectionKind::kNone, {})) {
      stage.discard();
      throw core::StageTransactionRejected(
          transaction.outcome().rejection,
          "baryonic source stage rejected (" +
              std::string(core::stageRejectionKindLabel(transaction.outcome().rejection)) +
              "): " + transaction.outcome().reason);
    }
    std::vector<std::uint32_t> touched;
    (void)transaction.commit([&] {
      if (local_work) touched = stage.commit(context.state);
    });
    if (local_work) {
      const physics::BaryonicSourceStageReport& report = stage.lastCommittedReport();
      m_last_thermochemistry_ledger = report.chemistry;
      m_last_ism_regime_report = report.regime;
      m_last_chemistry_adaptive_report = report.adaptive;
      m_last_dust_ledger = report.dust;
      m_last_sgs_ledger = report.sgs;
      m_last_shadow_report = report.shadow;
      m_last_skipped_effective_rows = report.skipped_effective_rows;
    }
    // Ranks without local work contribute zero current samples (never a stale
    // report from a previous stage). Observations are only ACCUMULATED here;
    // the collective controller update happens once at interval close
    // (closeSourceInterval), so neither a half step nor a level clock nor a
    // rejected trial can publish a controller state on its own.
    if (local_work) {
      accumulateSourceIntervalObservation(m_last_shadow_report,
                                          m_last_chemistry_adaptive_report);
    }
    if (!touched.empty()) {
      context.state.notifyGasThermodynamicMutation();
      std::unique_ptr<physics::EffectiveIsmThermodynamicClosure> closure;
      if (m_effective_eos_table != nullptr) {
        closure = std::make_unique<physics::EffectiveIsmThermodynamicClosure>(m_effective_eos_table);
      }
      const hydro::GasThermodynamicRefreshPolicy policy{
          .adiabatic_index = 5.0 / 3.0,
          .density_floor_code = m_hydro_density_floor_code,
          .pressure_floor_code = m_hydro_pressure_floor_code,
          .scale_factor = interval.evaluation_scale_factor,
          .redshift = interval.evaluation_redshift,
          .closure = closure.get(),
          .temperature_k_per_internal_energy_code =
              hydro::temperatureKelvinPerInternalEnergyCode(5.0 / 3.0, m_units.velocity_si_per_code)};
      std::vector<hydro::GasThermodynamicTouch> touches;
      touches.reserve(touched.size());
      for (const std::uint32_t row : touched) {
        touches.push_back(hydro::GasThermodynamicTouch{
            .cell_row = row, .internal_energy_before_code = context.state.gas_cells.internal_energy_code[row]});
      }
      (void)hydro::refreshGasThermodynamics(context.state, touches, policy);
    }
    if (context.profiler_session != nullptr) {
      const physics::SourceLedger& l = m_last_thermochemistry_ledger;
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "sources.thermochemistry",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "physics.thermochemistry",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "thermochemistry source ledger for the accepted interval",
          .payload = {{"cells", std::to_string(l.cells)},
                      {"skipped_effective_regime_cells", std::to_string(m_last_skipped_effective_rows)},
                      {"regime_changes", std::to_string(m_last_ism_regime_report.changed_rows.size())},
                      {"fidelity_promoted", std::to_string(m_last_chemistry_adaptive_report.promoted)},
                      {"fidelity_demoted", std::to_string(m_last_chemistry_adaptive_report.demoted)},
                      {"fidelity_promoted_by_cooling_time",
                       std::to_string(m_last_chemistry_adaptive_report.promoted_by_cooling_time)},
                      {"fidelity_promoted_by_error_history",
                       std::to_string(m_last_chemistry_adaptive_report.promoted_by_error_history)},
                      {"fidelity_held_by_dwell", std::to_string(m_last_chemistry_adaptive_report.held_by_dwell)},
                      {"dust_gas_to_dust_code", std::to_string(m_last_dust_ledger.gas_to_dust_code)},
                      {"dust_dust_to_gas_code", std::to_string(m_last_dust_ledger.dust_to_gas_code)},
                      {"sgs_dissipated_code", std::to_string(m_last_sgs_ledger.dissipated_to_heat_code)},
                      {"sgs_produced_from_resolved_code",
                       std::to_string(m_last_sgs_ledger.produced_from_resolved_code)},
                      {"sgs_resolved_stress_work_code",
                       std::to_string(m_last_sgs_ledger.resolved_stress_work_code)},
                      {"shadow_samples", std::to_string(m_last_shadow_report.samples)},
                      {"physical_seconds", std::to_string(interval.physical_seconds)},
                      {"recursive_level", std::to_string(recursive_level)},
                      {"radiative_loss_erg", std::to_string(l.radiative_loss_erg)},
                      {"heating_gain_erg", std::to_string(l.heating_gain_erg)},
                      {"chemical_binding_delta_erg", std::to_string(l.chemical_binding_delta_erg)},
                      {"thermal_energy_delta_erg", std::to_string(l.thermal_energy_delta_erg)},
                      {"floor_injection_erg", std::to_string(l.floor_injection_erg)},
                      {"accepted_substeps", std::to_string(l.accepted_substeps)},
                      {"rejected_substeps", std::to_string(l.rejected_substeps)},
                      {"newton_iterations", std::to_string(l.newton_iterations)},
                      {"metal_table_clamps", std::to_string(l.metal_table_clamps)},
                      {"max_error_estimate", std::to_string(l.max_error_estimate)},
                      {"element_projection_residual", std::to_string(l.element_mass_residual_code)},
                      {"charge_residual", std::to_string(l.charge_residual)},
                      {"metal_ion_potential_delta_erg", std::to_string(l.metal_ion_potential_delta_erg)},
                      {"metal_ion_co_transfer_nuclei", std::to_string(l.metal_ion_co_transfer_nuclei)},
                      {"metal_ion_normalization_residual", std::to_string(l.metal_ion_normalization_residual)},
                      {"metal_ion_inventory_residual", std::to_string(l.metal_ion_inventory_residual)},
                      {"shadow_ion_samples", std::to_string(m_last_shadow_report.ion_samples)},
                      {"shadow_max_ion_stage_fraction_difference",
                       std::to_string(m_last_shadow_report.max_ion_stage_fraction_difference)}}});
    }
  }

  // CHUI Part 3 Checkpoint F (C03/C04): global error-aware fidelity control.
  // Every rank reduces the sampled shadow disagreement (max) and sample count
  // (sum), applies the same deterministic control law, and persists the bias
  // in module sidecar "chui_adaptive_chemistry_control_v1" so a restart
  // continues with the identical decision history.
  // Parses the committed controller state from the canonical sidecar. Absent
  // sidecar => the documented initial state (bias 0). A malformed payload is
  // an invalid-state error, never a silent reset.
  [[nodiscard]] static double committedPromotionBias(const core::SimulationState& state) {
    const auto* stored = state.sidecars.find("chui_adaptive_chemistry_control_v1");
    if (stored == nullptr) return 0.0;
    const std::string text(reinterpret_cast<const char*>(stored->payload.data()), stored->payload.size());
    const auto pos = text.find("promotion_bias=");
    if (pos == std::string::npos) {
      throw std::runtime_error("chui_adaptive_chemistry_control_v1 sidecar has no promotion_bias entry");
    }
    std::size_t consumed = 0U;
    const double bias = std::stod(text.substr(pos + 15), &consumed);
    if (consumed == 0U || !std::isfinite(bias) || bias < 0.0) {
      throw std::runtime_error("chui_adaptive_chemistry_control_v1 sidecar has an invalid promotion_bias");
    }
    return bias;
  }

  void loadCommittedChemistryControl(const core::SimulationState& state) {
    if (m_thermochemistry == nullptr || !m_thermochemistry->config().adaptive.enabled) return;
    double bias = 0.0;
    std::exception_ptr load_failure;
    try {
      bias = m_config.numerics.v02_control_mode ==
                  core::V02ControlMode::kJointErrorResourceV1
          ? core::loadV02ControlState(state).chemistry_promotion_bias
          : committedPromotionBias(state);
    } catch (...) {
      load_failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        load_failure, "adaptive chemistry controller load");
    m_thermochemistry->setGlobalPromotionBias(bias);
  }

  void updateGlobalChemistryControl(
      core::StepContext& context,
      const physics::ChemistryShadowReport& current_shadow,
      const physics::ChemistryAdaptiveReport& interval_adaptive) {
    if (m_thermochemistry == nullptr || !m_thermochemistry->config().adaptive.enabled) return;
    // Next value is a pure function of the COMMITTED bias and this stage's
    // collective shadow reduction; it becomes authoritative only by being
    // written to the canonical sidecar (restored on rollback/restart).
    const bool joint =
        m_config.numerics.v02_control_mode == core::V02ControlMode::kJointErrorResourceV1;
    core::V02ControlState joint_state;
    double committed_bias = 0.0;
    std::exception_ptr load_failure;
    try {
      if (joint) {
        joint_state = core::loadV02ControlState(context.state);
        committed_bias = joint_state.chemistry_promotion_bias;
      } else {
        committed_bias = committedPromotionBias(context.state);
      }
    } catch (...) {
      load_failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        load_failure, "adaptive chemistry controller update load");
    // P1-49-07: the controlled disagreement covers the light network AND the
    // C02 ladders (stage fractions, ion electrons). It remains a LOCAL
    // sampled shadow disagreement, never a certified global error bar.
    const double local_ion_disagreement = std::max(
        current_shadow.max_ion_stage_fraction_difference, current_shadow.max_relative_ion_electron_difference);
    const double global_species_max =
        -m_mpi_context.allreduceMinDouble(-current_shadow.max_relative_species_difference);
    const double global_ion_max = -m_mpi_context.allreduceMinDouble(-local_ion_disagreement);
    const double global_max = std::max(global_species_max, global_ion_max);
    const std::uint64_t global_samples = m_mpi_context.allreduceSumUint64(current_shadow.samples);
    const std::uint64_t global_ion_samples = m_mpi_context.allreduceSumUint64(current_shadow.ion_samples);
    const double bias = physics::nextChemistryPromotionBias(committed_bias, global_max, global_samples);
    m_thermochemistry->setGlobalPromotionBias(bias);
    if (joint) {
      std::uint64_t local_detailed = 0U;
      for (const std::uint8_t tag : context.state.gas_composition.chemistry_model_tag) {
        local_detailed += tag == physics::k_chemistry_tag_detailed ? 1U : 0U;
      }
      joint_state.chemistry_promotion_bias = bias;
      joint_state.last_shadow_species_error = std::max(global_species_max, 0.0);
      joint_state.last_shadow_samples = global_samples;
      joint_state.last_shadow_ion_error = std::max(global_ion_max, 0.0);
      joint_state.last_shadow_ion_samples = global_ion_samples;
      joint_state.chemistry_detailed_rows =
          m_mpi_context.allreduceSumUint64(local_detailed);
      joint_state.chemistry_total_rows =
          m_mpi_context.allreduceSumUint64(
              static_cast<std::uint64_t>(context.state.gas_composition.rowCount()));
      joint_state.chemistry_promoted_last_interval =
          m_mpi_context.allreduceSumUint64(interval_adaptive.promoted);
      std::exception_ptr publication_failure;
      try {
        if (joint_state.epoch == std::numeric_limits<std::uint64_t>::max()) {
          throw std::overflow_error("V02 controller epoch overflow");
        }
        joint_state.epoch += 1U;
        core::storeV02ControlState(context.state, joint_state);
      } catch (...) {
        publication_failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(
          publication_failure, "V02 chemistry controller publication");
      return;
    }
    std::exception_ptr publication_failure;
    try {
      std::ostringstream text;
      text.precision(17);
      text << "module=chui_adaptive_chemistry_control\npromotion_bias=" << bias
           << "\nlast_global_shadow_max_relative_species_difference=" << global_species_max
           << "\nlast_global_shadow_samples=" << global_samples
           << "\nlast_global_shadow_ion_disagreement=" << global_ion_max
           << "\nlast_global_shadow_ion_samples=" << global_ion_samples << "\n";
      const std::string payload = text.str();
      core::ModuleSidecarBlock block;
      block.module_name = "chui_adaptive_chemistry_control_v1";
      block.schema_version = 1U;
      block.payload.resize(payload.size());
      for (std::size_t i = 0; i < payload.size(); ++i) {
        block.payload[i] = static_cast<std::byte>(payload[i]);
      }
      context.state.sidecars.upsert(std::move(block));
    } catch (...) {
      publication_failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        publication_failure, "fixed chemistry controller publication");
  }

  // CHUI Part 3 Checkpoint C: remote star-formation neighbours. The directed
  // AMR boundary exchange (same planner as hydro; one face layer; periodic
  // images) delivers every remote leaf cell that is a face neighbour of a
  // local patch, with its step-start centre and velocity; the patch metadata
  // gives its level and cell width. Ranks without patches still participate.
  void importRemoteStarFormationNeighbors(core::StepContext& context) {
    // Patch 53 (P0-52-01): local descriptor build voted before the exchange.
    std::vector<parallel::AmrPatchPayloadRecord> local_patch_records;
    {
      std::exception_ptr build_failure;
      try {
        local_patch_records =
            internal::buildMigrationAmrPatchPayloadRecords(context.state, m_mpi_context.worldRank());
      } catch (...) {
        build_failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(build_failure, "star-formation patch descriptor build");
    }
    // Received boundary cells (2x: geometric growth of the staging vector)
    // and their decoded RemoteCell coexist; admitted per peer inside the
    // stream's guarded admission phase (reservations declared first).
    std::vector<core::MemoryReservation> admitted;
    std::vector<parallel::AmrPatchCellPayloadRecord> received;
    const auto exchange = parallel::executeBlockingDirectedAmrPatchPayloadExchange(
        m_mpi_context, local_patch_records,
        [&context, this](std::span<const parallel::AmrPatchBoundaryCellRequest> requests, std::uint64_t first_record,
                         std::size_t max_records, std::vector<parallel::AmrPatchCellPayloadRecord>& output) {
          internal::fillMigrationAmrPatchBoundaryCellPayloadChunk(
              context.state, m_mpi_context.worldRank(), requests, first_record, max_records, output);
        },
        [&received](int, std::span<const parallel::AmrPatchCellPayloadRecord> records) {
          received.insert(received.end(), records.begin(), records.end());
        },
        [this, &admitted](int, std::span<const parallel::AmrPatchPayloadRecord>,
                          std::span<const parallel::AmrPatchBoundaryCellRequest>, std::uint64_t record_count) {
          if (m_memory_governor == nullptr || record_count == 0U) return;
          core::MemoryReservation reservation = m_memory_governor->reserve(
              core::MemoryClass::kPhaseResident,
              core::checkedSizeMultiply(
                  core::checkedIntegralNarrow<std::size_t>(record_count, "star-formation remote cell count"),
                  2U * sizeof(parallel::AmrPatchCellPayloadRecord) +
                      sizeof(internal::StarFormationNeighborIndex::RemoteCell),
                  "star-formation remote cell bytes"),
              "sources.star_formation.remote_neighbours");
          reservation.commit();
          admitted.push_back(std::move(reservation));
        },
        0U, context.integrator_state.step_index,
        parallel::AmrInterfaceTopologyOptions{
            .periodic_box_extent_comoving = m_hydro_boundary == core::BoundaryCondition::kPeriodic
                ? m_box_extent_comoving : std::array<double, 3>{0.0, 0.0, 0.0},
            .ghost_layers = 1U},
        m_memory_governor);
    // Patch 52 (P1-51-01): the decode and index import are a guarded phase
    // with a collective vote, because the star-formation batch rounds that
    // follow enter collectives; a local failure must not strand peers there.
    std::exception_ptr failure;
    try {
      std::unordered_map<std::uint64_t, parallel::AmrPatchPayloadRecord> patch_by_id;
      for (const auto& record : exchange.patch_payloads_received) patch_by_id.emplace(record.patch_id, record);
      std::vector<internal::StarFormationNeighborIndex::RemoteCell> remote;
      remote.reserve(received.size());
      for (const auto& cell : received) {
        const auto it = patch_by_id.find(cell.patch_id);
        if (it == patch_by_id.end()) {
          throw std::runtime_error("star-formation remote neighbour has no patch metadata (fail closed)");
        }
        const auto& patch = it->second;
        remote.push_back(internal::StarFormationNeighborIndex::RemoteCell{
            .center = {cell.center_x_comoving, cell.center_y_comoving, cell.center_z_comoving},
            .width = {patch.extent_x_comoving / patch.cell_dim_x, patch.extent_y_comoving / patch.cell_dim_y,
                      patch.extent_z_comoving / patch.cell_dim_z},
            .level = patch.level,
            .velocity = {cell.velocity_x_peculiar, cell.velocity_y_peculiar, cell.velocity_z_peculiar}});
      }
      if (m_sf_neighbor_index != nullptr) m_sf_neighbor_index->addRemoteCells(std::move(remote));
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(failure, "star-formation remote neighbour import");
  }

  // CHUI Part 3 Checkpoint C: remote feedback halo. For every staged source
  // the local k-NN radius bounds the global stencil; each other rank returns
  // its own k nearest owned leaf cells within that radius (step-start
  // properties), so the source owner sees the exact global k-NN stencil and
  // normalizes its weights over local AND remote cells (never a rank-local
  // renormalization). Query volume: O(active sources x ranks) small records.
  //
  // Patch 52 (P1-51-01/02): every rank-local preparation (query build,
  // response build, candidate import) runs in a guarded phase followed by a
  // collective vote, and both streams are received directly into governed
  // typed buffers (no ungoverned byte buffer + second typed copy). Patch 53:
  // the candidate storage AND its lease are handed to `deposition`, which
  // owns both (and admits its own index growth) for its whole lifetime.
  void importRemoteFeedbackCandidates(
      const core::SimulationState& state,
      physics::MechanicalFeedbackDeposition& deposition,
      std::span<const physics::FeedbackSourceBudget> budgets,
      std::uint64_t sequence) {
    struct Query {
      std::array<double, 3> position{};
      double radius2 = 0.0;
      int requester_rank = 0;
      int pad = 0;
    };
    static_assert(std::is_trivially_copyable_v<Query>);
    const int rank = m_mpi_context.worldRank();
    const int ranks = m_mpi_context.worldSize();
    const auto active = [](const physics::FeedbackSourceBudget& budget) {
      // Same activity predicate as MechanicalFeedbackDeposition::stage: a
      // wind-only or dispersion-only source still deposits, so it must see the
      // global stencil (previously skipped => rank-local stencil).
      return !(budget.ejecta_mass_code == 0.0 && budget.supernova_energy_erg == 0.0 &&
               budget.presn_momentum_code == 0.0 && budget.wind_energy_erg == 0.0 &&
               budget.wind_momentum_code == 0.0 && budget.association_dispersion_energy_code == 0.0);
    };
    core::MemoryReservation query_lease;  // outbound queries + destinations
    std::vector<Query> queries;
    std::vector<int> destinations;
    std::exception_ptr failure;
    try {
      std::size_t active_sources = 0U;
      for (const auto& budget : budgets) active_sources += active(budget) ? 1U : 0U;
      const std::size_t query_count = core::checkedSizeMultiply(
          active_sources, static_cast<std::size_t>(std::max(ranks - 1, 0)), "remote feedback halo query count");
      if (m_memory_governor != nullptr && query_count != 0U) {
        query_lease = m_memory_governor->reserve(
            core::MemoryClass::kCommunication,
            core::checkedSizeMultiply(query_count, sizeof(Query) + sizeof(int), "remote feedback halo query bytes"),
            "sources.feedback.remote_halo_queries");
        query_lease.commit();
      }
      queries.reserve(query_count);
      destinations.reserve(query_count);
      for (const auto& budget : budgets) {
        if (!active(budget)) continue;
        const double r2 = deposition.localStencilRadius2(state, budget.position_comoving);
        for (int peer = 0; peer < ranks; ++peer) {
          if (peer == rank) continue;
          queries.push_back(Query{.position = budget.position_comoving, .radius2 = r2, .requester_rank = rank});
          destinations.push_back(peer);
        }
      }
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(failure, "remote feedback halo queries");
    auto inbound = parallel::exchangeGovernedRecords<Query>(
        m_mpi_context, queries, destinations, m_memory_governor, "sources.feedback.remote_halo_query_exchange");
    core::MemoryReservation response_lease;  // outbound candidates + destinations
    std::vector<physics::FeedbackRemoteCandidate> responses;
    std::vector<int> response_destinations;
    try {
      // At most k (neighbor_count) local candidates answer one query.
      const std::size_t response_bound = core::checkedSizeMultiply(
          inbound.records.size(), static_cast<std::size_t>(std::max<std::uint32_t>(m_physics_config.fb_neighbor_count, 1U)),
          "remote feedback halo response bound");
      // P2-MEMORY-01: one phase-reused k-NN heap and candidate-row buffer
      // (capacity k each) serve every inbound query; they are admitted with
      // the responses and never reallocate inside the query loop.
      const std::size_t knn_capacity =
          static_cast<std::size_t>(std::max<std::uint32_t>(m_physics_config.fb_neighbor_count, 1U));
      if (m_memory_governor != nullptr && response_bound != 0U) {
        response_lease = m_memory_governor->reserve(
            core::MemoryClass::kCommunication,
            core::checkedSizeAdd(
                core::checkedSizeMultiply(response_bound, sizeof(physics::FeedbackRemoteCandidate) + sizeof(int),
                                          "remote feedback halo response bytes"),
                core::checkedSizeMultiply(knn_capacity,
                                          sizeof(std::pair<double, std::uint32_t>) + sizeof(std::uint32_t),
                                          "remote feedback halo k-NN scratch bytes"),
                "remote feedback halo response and scratch bytes"),
            "sources.feedback.remote_halo_responses");
        response_lease.commit();
      }
      std::vector<std::pair<double, std::uint32_t>> knn_heap_scratch;
      std::vector<std::uint32_t> knn_rows;
      if (!inbound.records.empty()) {
        knn_heap_scratch.reserve(knn_capacity);
        knn_rows.reserve(knn_capacity);
      }
      for (const Query& query : inbound.records) {
        if (query.requester_rank < 0 || query.requester_rank >= ranks || query.requester_rank == rank) {
          throw std::runtime_error("remote feedback halo: query from an invalid requester rank");
        }
        deposition.localCandidatesWithin(state, query.position, query.radius2, knn_heap_scratch, knn_rows);
        for (const std::uint32_t row : knn_rows) {
          physics::FeedbackRemoteCandidate candidate;
          candidate.gas_cell_id = state.gas_cells.gas_cell_id[row];
          candidate.owner_rank = rank;
          candidate.center_comoving = {state.cells.center_x_comoving[row], state.cells.center_y_comoving[row],
                                       state.cells.center_z_comoving[row]};
          candidate.mass_code = state.cells.mass_code[row];
          candidate.momentum_code = {candidate.mass_code * state.gas_cells.velocity_x_peculiar[row],
                                     candidate.mass_code * state.gas_cells.velocity_y_peculiar[row],
                                     candidate.mass_code * state.gas_cells.velocity_z_peculiar[row]};
          candidate.density_code = state.gas_cells.density_code[row];
          candidate.hydrogen_mass_fraction = state.gas_composition.enabled()
              ? state.gas_composition.value(state.gas_composition.descriptor().element_lane_begin, row) /
                    std::max(physics::compositionGasMass(state.gas_composition, row), 1.0e-300)
              : 0.76;
          responses.push_back(candidate);
          response_destinations.push_back(query.requester_rank);
        }
      }
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(failure, "remote feedback halo responses");
    inbound.releaseStorage();
    auto remote = parallel::exchangeGovernedRecords<physics::FeedbackRemoteCandidate>(
        m_mpi_context, responses, response_destinations, m_memory_governor,
        "sources.feedback.remote_halo_candidate_exchange");
    try {
      deposition.setRemoteCandidates(state, std::move(remote.records), std::move(remote.lease));
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(failure, "remote feedback halo import");
  }

  // G0.3: decomposition-independent association discovery and merge.
  //  1. every rank contributes compact descriptors of its ACTIVE feedback
  //     sources (non-zero budget this interval); one bounded all-gather
  //     gives every rank the identical global descriptor set (count reduced
  //     and admitted first);
  //  2. every rank runs the same deterministic friends-of-friends discovery
  //     (physics::discoverAssociationLeaders) => identical leader ids;
  //  3. member budgets are routed to the rank owning the leader (minimum id)
  //     and merged there in ascending source-id order, so the merged budget
  //     is bitwise independent of the MPI partition;
  //  4. inactive (zero-budget) sources are dropped exactly as before.
  // Every collective is entered by every rank (failures are reduced first).
  std::vector<physics::FeedbackSourceBudget> mergeAssociationsDistributed(
      const std::vector<physics::FeedbackSourceBudget>& budgets,
      const std::vector<double>& birth_years,
      const physics::FeedbackUnits& units,
      bool distributed,
      std::uint64_t sequence) {
    if (budgets.size() != birth_years.size()) {
      throw std::logic_error("association merge: budget/birth-time lists differ in length");
    }
    const auto active = [](const physics::FeedbackSourceBudget& b) {
      return b.ejecta_mass_code != 0.0 || b.supernova_energy_erg != 0.0 || b.presn_momentum_code != 0.0 ||
          b.wind_energy_erg != 0.0 || b.wind_momentum_code != 0.0 || b.association_dispersion_energy_code != 0.0;
    };
    const int rank = distributed ? m_mpi_context.worldRank() : 0;
    std::vector<physics::AssociationSourceDescriptor> local;
    std::vector<std::size_t> local_budget_index;
    std::exception_ptr failure;
    try {
      for (std::size_t i = 0; i < budgets.size(); ++i) {
        if (!active(budgets[i])) continue;
        physics::AssociationSourceDescriptor d;
        d.source_id = budgets[i].source_id;
        d.position_comoving = budgets[i].position_comoving;
        d.birth_years = birth_years[i];
        d.owner_rank = rank;
        local.push_back(d);
        local_budget_index.push_back(i);
      }
    } catch (...) {
      failure = std::current_exception();
    }
    if (distributed) m_mpi_context.rethrowCollectivePreparationFailure(failure, "association descriptors");
    else if (failure != nullptr) std::rethrow_exception(failure);

    std::vector<physics::AssociationSourceDescriptor> global;
    core::MemoryReservation global_reservation;
    if (distributed) {
      const std::uint64_t global_count = m_mpi_context.allreduceSumUint64(local.size());
      try {
        if (m_memory_governor != nullptr) {
          // Gathered bytes, decoded descriptors and discovery scratch
          // (entries, parent, leader ids, sort keys) coexist.
          const std::uint64_t per_source = 2U * sizeof(physics::AssociationSourceDescriptor) +
              sizeof(std::array<std::int64_t, 3>) + 4U * sizeof(std::uint64_t) + 2U * sizeof(std::size_t);
          global_reservation = m_memory_governor->reserve(
              core::MemoryClass::kPhaseResident, global_count * per_source, "sources.associations.global_discovery");
          global_reservation.commit();
        }
      } catch (...) {
        failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(failure, "association discovery admission");
      const auto bytes = m_mpi_context.allgatherBytesBounded(std::span<const std::uint8_t>(
          reinterpret_cast<const std::uint8_t*>(local.data()),
          local.size() * sizeof(physics::AssociationSourceDescriptor)));
      // Patch 52: the decode is a guarded phase (the next collective is the
      // discovery vote below, so a local failure here is reported there).
      try {
        if (bytes.size() != global_count * sizeof(physics::AssociationSourceDescriptor)) {
          throw std::runtime_error("association discovery: gathered descriptor size mismatch");
        }
        global.resize(static_cast<std::size_t>(global_count));
        if (!bytes.empty()) std::memcpy(global.data(), bytes.data(), bytes.size());
      } catch (...) {
        failure = std::current_exception();
      }
    } else {
      global = local;
    }
    std::vector<std::uint64_t> leaders;
    try {
      if (failure != nullptr) std::rethrow_exception(failure);
      leaders = physics::discoverAssociationLeaders(
          global, m_physics_config.ssp_association_radius_code,
          m_physics_config.ssp_association_time_window_years, units.periodic, units.box_extent_comoving);
    } catch (...) {
      failure = std::current_exception();
    }
    if (distributed) m_mpi_context.rethrowCollectivePreparationFailure(failure, "association discovery");
    else if (failure != nullptr) std::rethrow_exception(failure);

    struct MemberRecord {
      std::uint64_t leader_id = 0U;
      physics::FeedbackSourceBudget budget{};
    };
    static_assert(std::is_trivially_copyable_v<MemberRecord>);
    // Patch 52 (P1-51-01/02): routing preparation is a guarded phase with a
    // collective vote before the member exchange; members are received into
    // a governed typed buffer and appended to `kept` in a guarded phase
    // (the final merge has no further collective in this function).
    // Patch 53 (P1-52-02): every coexisting routing/merge buffer is admitted
    // before it grows, with its lease declared before the storage it covers:
    //   routing_lease: outbound + destinations + kept (reserved at local
    //     size, so routing pushes never reallocate);
    //   routing_scratch_lease (inside the routing scope): owner_of,
    //     id_to_leader, local_leader;
    //   merge_lease: kept regrowth to local + received (old and new capacity
    //     coexist during the reallocation), members and merged (each <= the
    //     member count). The received member records carry their own lease
    //     (governed exchange) and are not charged again.
    core::MemoryReservation routing_lease;
    core::MemoryReservation merge_lease;
    std::vector<MemberRecord> outbound;
    std::vector<int> destinations;
    std::vector<MemberRecord> kept;
    try {
      const std::size_t n_local = local.size();
      if (m_memory_governor != nullptr && n_local != 0U) {
        routing_lease = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            core::checkedSizeMultiply(n_local, 2U * sizeof(MemberRecord) + sizeof(int),
                                      "association routing buffer bytes"),
            "sources.associations.routing");
        routing_lease.commit();
      }
      outbound.reserve(n_local);
      destinations.reserve(n_local);
      kept.reserve(n_local);
      core::MemoryReservation routing_scratch_lease;
      if (m_memory_governor != nullptr && (global.size() + n_local) != 0U) {
        routing_scratch_lease = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            core::checkedSizeAdd(
                core::checkedSizeMultiply(global.size(),
                                          sizeof(std::pair<std::uint64_t, int>) +
                                              sizeof(std::pair<std::uint64_t, std::uint64_t>),
                                          "association owner/leader index bytes"),
                core::checkedSizeMultiply(n_local, sizeof(std::uint64_t), "association local leader bytes"),
                "association routing scratch bytes"),
            "sources.associations.routing_scratch");
        routing_scratch_lease.commit();
      }
      // Leader id -> owner rank (the descriptor whose id equals the leader).
      std::vector<std::pair<std::uint64_t, int>> owner_of;
      owner_of.reserve(global.size());
      for (const auto& d : global) owner_of.emplace_back(d.source_id, d.owner_rank);
      std::sort(owner_of.begin(), owner_of.end());
      const auto owner_rank_of = [&owner_of](std::uint64_t id) {
        const auto it = std::lower_bound(owner_of.begin(), owner_of.end(), std::make_pair(id, -1));
        if (it == owner_of.end() || it->first != id) {
          throw std::logic_error("association discovery: leader id has no descriptor");
        }
        return it->second;
      };
      std::vector<std::uint64_t> local_leader(local.size());
      {
        std::vector<std::pair<std::uint64_t, std::uint64_t>> id_to_leader(global.size());
        for (std::size_t g = 0; g < global.size(); ++g) id_to_leader[g] = {global[g].source_id, leaders[g]};
        std::sort(id_to_leader.begin(), id_to_leader.end());
        for (std::size_t l = 0; l < local.size(); ++l) {
          const auto it = std::lower_bound(id_to_leader.begin(), id_to_leader.end(),
                                           std::make_pair(local[l].source_id, std::uint64_t{0}));
          if (it == id_to_leader.end() || it->first != local[l].source_id) {
            throw std::logic_error("association discovery: local source missing from the global descriptor set");
          }
          local_leader[l] = it->second;
        }
      }
      for (std::size_t l = 0; l < local.size(); ++l) {
        MemberRecord rec{local_leader[l], budgets[local_budget_index[l]]};
        const int dest = owner_rank_of(local_leader[l]);
        if (dest == rank) {
          kept.push_back(rec);
        } else {
          outbound.push_back(rec);
          destinations.push_back(dest);
        }
      }
      if (!distributed && !outbound.empty()) {
        throw std::logic_error("association merge: single-rank merge routed a member off-rank");
      }
    } catch (...) {
      failure = std::current_exception();
    }
    if (distributed) m_mpi_context.rethrowCollectivePreparationFailure(failure, "association member routing");
    else if (failure != nullptr) std::rethrow_exception(failure);
    if (distributed) {
      auto received = parallel::exchangeGovernedRecords<MemberRecord>(
          m_mpi_context, outbound, destinations, m_memory_governor, "sources.associations.member_exchange");
      try {
        const std::size_t members_total = core::checkedSizeAdd(kept.size(), received.records.size(),
                                                               "association member count");
        if (m_memory_governor != nullptr && members_total != 0U) {
          merge_lease = m_memory_governor->reserve(
              core::MemoryClass::kPhaseResident,
              core::checkedSizeAdd(
                  core::checkedSizeMultiply(members_total, sizeof(MemberRecord), "association kept regrowth bytes"),
                  core::checkedSizeMultiply(members_total, 2U * sizeof(physics::FeedbackSourceBudget),
                                            "association members/merged bytes"),
                  "association merge bytes"),
              "sources.associations.merge");
          merge_lease.commit();
        }
        kept.reserve(members_total);
        kept.insert(kept.end(), received.records.begin(), received.records.end());
      } catch (...) {
        failure = std::current_exception();
      }
    } else if (m_memory_governor != nullptr && !kept.empty()) {
      try {
        merge_lease = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            core::checkedSizeMultiply(kept.size(), 2U * sizeof(physics::FeedbackSourceBudget),
                                      "association members/merged bytes"),
            "sources.associations.merge");
        merge_lease.commit();
      } catch (...) {
        failure = std::current_exception();
      }
    }
    std::vector<physics::FeedbackSourceBudget> merged;
    try {
      if (failure != nullptr) std::rethrow_exception(failure);
      merged.reserve(kept.size());
      std::sort(kept.begin(), kept.end(), [](const MemberRecord& l, const MemberRecord& r) {
        return l.leader_id != r.leader_id ? l.leader_id < r.leader_id : l.budget.source_id < r.budget.source_id;
      });
      std::vector<physics::FeedbackSourceBudget> members;
      members.reserve(kept.size());
      std::size_t begin = 0U;
      while (begin < kept.size()) {
        std::size_t end = begin + 1U;
        while (end < kept.size() && kept[end].leader_id == kept[begin].leader_id) ++end;
        members.clear();
        for (std::size_t k = begin; k < end; ++k) members.push_back(kept[k].budget);
        merged.push_back(physics::mergeAssociationMembers(members, units.periodic, units.box_extent_comoving));
        begin = end;
      }
    } catch (...) {
      failure = std::current_exception();
    }
    // The caller's next step is a collective (remote halo query vote), so a
    // local merge failure is reported collectively here (patch 52).
    if (distributed) m_mpi_context.rethrowCollectivePreparationFailure(failure, "association member merge");
    else if (failure != nullptr) std::rethrow_exception(failure);
    return merged;
  }

  // CHUI Part 2 stellar populations: age-resolved SSP mass/element return
  // (AGB, CC continuous; Ia per event), discrete SN events from the
  // restart-safe cumulative-Poisson chronology, and conservative mechanical
  // feedback. Staged in a StageTransaction: event ordinals, star masses/ages,
  // population records and gas deltas commit together or not at all.
  void executePart2StellarPopulation(core::StepContext& context, double source_evaluation_scale_factor) {
    core::SimulationState& state = context.state;
    const std::uint64_t ranks_with_stars =
        m_mpi_context.allreduceSumUint64(state.star_particles.size() > 0U ? 1ULL : 0ULL);
    // CHUI Part 3 Checkpoint C: cross-rank deposition (Part 2 refused here).
    // Every rank participates when any rank has stars, so the collective
    // halo and delta exchanges never deadlock.
    const bool distributed_feedback =
        ranks_with_stars > 0U && m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1;
    if (!distributed_feedback && state.star_particles.size() == 0U) return;
    if (ranks_with_stars == 0U) return;
    const double elapsed_years = elapsedStellarEvolutionYears(context);
    const bool explicit_event_interval =
        std::isfinite(context.timeline_step.time_begin_code) &&
        std::isfinite(context.timeline_step.time_end_code) &&
        context.timeline_step.time_end_code > context.timeline_step.time_begin_code;
    const double event_interval_begin_code = explicit_event_interval
        ? context.timeline_step.time_begin_code
        : context.integrator_state.current_time_code;
    const double event_interval_end_code = explicit_event_interval
        ? context.timeline_step.time_end_code
        : context.integrator_state.current_time_code +
              context.integrator_state.dt_time_code;
    const physics::ThermochemistryUnitBridge bridge =
        physics::makeThermochemistryOperatorConfig(m_full_config, m_units, m_baryonic_data).units;
    physics::FeedbackUnits units;
    units.mass_code_to_g = bridge.mass_code_to_g;
    units.velocity_code_to_cm_s = bridge.velocity_code_to_cm_s;
    units.length_code_to_cm = bridge.length_code_to_cm;
    units.density_code_to_cgs = bridge.density_code_to_cgs;
    units.scale_factor = source_evaluation_scale_factor;
    units.comoving = bridge.comoving;
    units.periodic = m_hydro_boundary == core::BoundaryCondition::kPeriodic;
    units.box_extent_comoving = m_box_extent_comoving;
    physics::MechanicalFeedbackConfig fb_config;
    fb_config.neighbor_count = std::max<std::uint32_t>(m_physics_config.fb_neighbor_count, 1U);
    fb_config.sgs_energy_fraction = m_physics_config.sgs_feedback_energy_fraction;
    const physics::DustConfig dust_config = physics::makeDustConfig(m_physics_config);
    fb_config.dust_sn_swept_mass_per_1e51_erg_g = dust_config.sn_swept_mass_per_1e51_erg_g;
    fb_config.dust_sn_destruction_efficiency = dust_config.sn_destruction_efficiency;
    // Decomposition-invariant deposition (step-start target state + per-cell
    // kinetic closure); identical on 1 and N ranks.
    fb_config.step_start_reference = true;

    // Patch 53: `deposition` owns its candidate indices, the remote candidate
    // storage and their leases (setMemoryGovernor / setRemoteCandidates).
    physics::MechanicalFeedbackDeposition deposition(fb_config, units);  // rank-uniform config check only
    deposition.setMemoryGovernor(m_memory_governor);
    // Patch 52 (P1-51-01): every rank-local allocation of the stage prologue
    // (owned-leaf metadata, candidate index, population store, sorted star
    // rows) precedes one collective vote, so a local failure cannot strand
    // peers in the first collective of the distributed feedback path.
    std::optional<physics::StellarPopulationStore> store_storage;
    std::vector<std::uint32_t> star_rows;
    {
      std::exception_ptr prologue_failure;
      try {
        buildOwnedLeafCellMetadata(context);
        std::vector<std::uint32_t> owned_leaf;
        owned_leaf.reserve(state.cells.size());
        for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
          if (row < m_owned_leaf_mask.size() && m_owned_leaf_mask[row] != 0U) owned_leaf.push_back(row);
        }
        deposition.setCandidates(state, std::move(owned_leaf));
        store_storage.emplace(state);
        if (!context.active_set.particles_are_subset || context.active_set.particle_indices.empty()) {
          star_rows.reserve(state.star_particles.size());
          for (std::uint32_t s = 0; s < state.star_particles.size(); ++s) star_rows.push_back(s);
        } else {
          std::pmr::vector<std::uint32_t> rows(std::pmr::get_default_resource());
          buildActiveStarRows(context, rows, std::pmr::get_default_resource());
          star_rows.assign(rows.begin(), rows.end());
        }
        // Deterministic source order: stable particle id (not row / rank order).
        std::sort(star_rows.begin(), star_rows.end(), [&](std::uint32_t a, std::uint32_t b) {
          return state.particle_sidecar.particle_id[state.star_particles.particle_index[a]] <
              state.particle_sidecar.particle_id[state.star_particles.particle_index[b]];
        });
      } catch (...) {
        prologue_failure = std::current_exception();
      }
      if (distributed_feedback) {
        m_mpi_context.rethrowCollectivePreparationFailure(prologue_failure, "stellar population stage prologue");
      } else if (prologue_failure != nullptr) {
        std::rethrow_exception(prologue_failure);
      }
    }
    physics::StellarPopulationStore& store = *store_storage;

    struct StarUpdate {
      std::uint32_t star_index = 0;
      std::uint32_t particle_index = 0;
      std::uint64_t particle_id = 0;
      double returned_mass_code = 0.0;
      double age_end_years = 0.0;
      double energy_erg = 0.0;
      std::array<double, 3> channel_mass{};
      std::array<double, 3> channel_energy{};
      physics::StellarPopulationRecord record{};
    };
    // Patch 53 (M01): per-star stage staging (updates, budgets, association
    // birth times) admitted at the active star count before it grows; the
    // lease is declared before the storage it covers.
    core::MemoryReservation ssp_stage_lease;
    std::vector<StarUpdate> updates;
    const double msun_code = physics::k_solar_mass_g / bridge.mass_code_to_g;
    const bool energy_coupled = m_physics_config.enable_feedback;
    std::uint64_t newborn_skipped = 0U, legacy_initialized = 0U;
    const internal::MpiStageCollectiveDecision decision(m_mpi_context);
    core::StageTransaction transaction("sources.stellar_population", m_memory_governor, &decision);
    // SSP diagnostics are restored if the stage is rejected.
    const std::uint64_t yield_clamps_before = m_ssp_model->table_clamps;
    const std::uint64_t negative_clamps_before = m_ssp_model->negative_ejecta_clamps;
    const double untracked_residual_before = m_ssp_model->untracked_mass_residual_msun;
    const std::uint64_t nonmonotone_before = m_ssp_model->nonmonotone_return_clamps;
    std::vector<physics::FeedbackSourceBudget> staged_budgets;
    std::vector<double> association_birth_years;
    DiscreteEventChronologyLedger event_chronology;
    event_chronology.interval_begin_code = event_interval_begin_code;
    event_chronology.interval_end_code = event_interval_end_code;
    event_chronology.lag_bound_years = enforcedEventLagBoundYears();
    event_chronology.physical_event_barrier = physicalEventBarrier();
    (void)transaction.prepare([&] {
      std::exception_ptr loop_failure;
      try {
      if (m_memory_governor != nullptr && !star_rows.empty()) {
        ssp_stage_lease = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            core::checkedSizeMultiply(star_rows.size(),
                                      sizeof(StarUpdate) + sizeof(physics::FeedbackSourceBudget) + sizeof(double),
                                      "SSP stage per-star staging bytes"),
            "sources.stellar_population.stage");
        ssp_stage_lease.commit();
      }
      updates.reserve(star_rows.size());
      staged_budgets.reserve(star_rows.size());
      association_birth_years.reserve(star_rows.size());
      for (const std::uint32_t s : star_rows) {
        const std::uint32_t p = state.star_particles.particle_index[s];
        const std::uint64_t pid = state.particle_sidecar.particle_id[p];
        // Audit F5: born at the end of this interval => zero age in it.
        if (physics::stellarElapsedYearsInInterval(state.star_particles.birth_tick[s],
                                                   context.integrator_state.step_index, elapsed_years) <= 0.0) {
          ++newborn_skipped;
          continue;
        }
        physics::StellarPopulationRecord record = store.get(pid);
        if (record.initialized == 0U) {
          // Explicit legacy initialization for stars created before Part 2:
          // birth mixture = primordial H/He + solar-pattern metals at the
          // stored birth metallicity; birth time unknown (0); stream key =
          // the immutable birth key. Counted, never silent.
          record = legacyStellarPopulationRecord(state, s);
          // A legacy star has no absolute event cursor. Its cumulative age
          // lanes/ordinals remain authoritative for past evolution; the
          // current interval establishes the first absolute boundary.
          record.chronology_initialized = 0U;
          ++legacy_initialized;
        }
        // P1-TIME-01: the time authority forecasts every star, legacy records
        // included (discreteEventStepEndLimitTimeCode builds the same legacy
        // record and cursor adoption), so every TIMED event is held to the
        // enforced bound; catch-up events stay counted and exempt.
        const bool event_forecast_covered = true;
        physics::advanceStellarPopulationChronology(
            record,
            event_interval_begin_code,
            event_interval_end_code,
            // Patch 52: discrete events run once per accepted coarse interval
            // (never on a level clock); patch 51 referenced an identifier
            // that is not in scope here.
            context.integrator_state.step_index);
        const double birth_mass = state.star_particles.birth_mass_code[s];
        const double birth_msun = birth_mass / msun_code;
        const double age_b = std::max(state.star_particles.stellar_age_years_last[s], 0.0);
        const double age_e = age_b + elapsed_years;
        const auto& x = record.birth_element_mass_fraction;
        double z = 0.0;
        for (std::size_t e = 2; e < physics::k_chemical_element_count; ++e) z += x[e];
        const auto agb = m_ssp_model->massiveReturn(0, age_b, age_e, x);
        const auto cc = m_ssp_model->massiveReturn(1, age_b, age_e, x);
        double n_cc = 0.0, n_ia = 0.0;
        const double lam_cc_end = m_ssp_model->ccCumulativeNumber(age_e, z) * birth_msun;
        const double lam_ia_end = m_ssp_model->iaCumulativeNumber(age_e) * birth_msun;
        if (m_ssp_model->config().stochastic_events) {
          // T02 (patch 52): each realized event keeps its model-defined
          // physical age (threshold crossing of Lambda_c); it is published at
          // the synchronized interval end with a recorded, bounded lag.
          const physics::StellarPopulationRecord record_before = record;
          const auto cc_events = physics::realizeEvents(record, 0U, lam_cc_end);
          const auto ia_events = physics::realizeEvents(record, 1U, lam_ia_end);
          n_cc = static_cast<double>(cc_events.count);
          n_ia = static_cast<double>(ia_events.count);
          for (std::uint32_t channel = 0; channel < 2U; ++channel) {
            event_chronology.accumulate(
                physics::realizedEventChronology(*m_ssp_model, record_before, channel,
                                                 channel == 0U ? cc_events.count : ia_events.count, z, birth_msun,
                                                 age_b, age_e),
                pid, age_b, age_e, event_forecast_covered);
          }
        } else {
          n_cc = lam_cc_end - m_ssp_model->ccCumulativeNumber(age_b, z) * birth_msun;
          n_ia = lam_ia_end - m_ssp_model->iaCumulativeNumber(age_b) * birth_msun;
        }
        const auto ia = m_ssp_model->iaReturnPerEvent(x);
        std::array<double, physics::k_chemical_element_count> agb_el{}, cc_el{}, ia_el{};
        for (std::size_t e = 0; e < physics::k_chemical_element_count; ++e) {
          agb_el[e] = agb.element_fraction[e] * birth_mass;
          cc_el[e] = cc.element_fraction[e] * birth_mass;
          ia_el[e] = ia.element_fraction[e] * n_ia * msun_code;
        }
        double returned = agb.mass_fraction * birth_mass + cc.mass_fraction * birth_mass +
            ia.mass_fraction * n_ia * msun_code;
        // Cumulative returned mass never exceeds the star's current mass.
        const double current_mass = state.particles.mass_code[p];
        double scale = 1.0;
        if (returned > current_mass && returned > 0.0) scale = current_mass / returned;
        returned *= scale;
        const auto agb_part = physics::partitionEjectaIntoDust(agb_el, dust_config.agb_condensation_efficiency);
        const auto cc_part = physics::partitionEjectaIntoDust(cc_el, dust_config.cc_condensation_efficiency);
        physics::FeedbackSourceBudget budget;
        budget.source_id = pid;
        budget.position_comoving = {state.particles.position_x_comoving[p], state.particles.position_y_comoving[p],
                                    state.particles.position_z_comoving[p]};
        budget.velocity_peculiar = {state.particles.velocity_x_peculiar[p], state.particles.velocity_y_peculiar[p],
                                    state.particles.velocity_z_peculiar[p]};
        budget.ejecta_mass_code = returned;
        const bool dust_lanes = state.gas_composition.enabled() &&
            state.gas_composition.descriptor().dust_lane_count > 0U;
        for (std::size_t e = 0; e < physics::k_chemical_element_count; ++e) {
          // Without dust lanes the condensable mass stays in the gas.
          budget.ejecta_elements_code[e] = scale *
              (dust_lanes ? agb_part.gas_elements[e] + cc_part.gas_elements[e] : agb_el[e] + cc_el[e]) +
              scale * ia_el[e];
        }
        if (dust_lanes) {
          for (std::size_t r = 0; r < physics::k_dust_reservoir_count; ++r) {
            budget.ejecta_dust_code[r] = scale * (agb_part.dust[r] + cc_part.dust[r]);
          }
        }
        // Corrective gate: the deposited mass is DEFINED by the element and
        // dust budgets so the gas mass and the composition lanes change by
        // exactly the same amount (the star loses the same mass).
        {
          double element_sum = 0.0;
          for (const double m : budget.ejecta_elements_code) element_sum += m;
          for (const double m : budget.ejecta_dust_code) element_sum += m;
          if (element_sum > current_mass && element_sum > 0.0) {
            // Never return more than the star holds (exact, all lanes scaled).
            const double cap = current_mass / element_sum;
            for (double& m : budget.ejecta_elements_code) m *= cap;
            for (double& m : budget.ejecta_dust_code) m *= cap;
            element_sum = current_mass;
          }
          returned = element_sum;
          budget.ejecta_mass_code = element_sum;
        }
        const double e_sn = energy_coupled ? (n_cc + n_ia) * m_ssp_model->config().energy_per_sn_erg : 0.0;
        budget.supernova_count = energy_coupled ? n_cc + n_ia : 0.0;
        budget.supernova_energy_erg = e_sn;
        if (energy_coupled && m_physics_config.fb_presn_wind_momentum_lsun_per_c > 0.0 &&
            age_b < m_ssp_model->lifetimeYears(m_ssp_model->config().cc_mass_min_msun, z)) {
          // EXPERIMENTAL user closure: p_dot = q M_* L_sun / c for ages
          // before the last core-collapse SN (no SSP luminosity table).
          const double p_dot_cgs = m_physics_config.fb_presn_wind_momentum_lsun_per_c * birth_msun * 3.828e33 / 2.99792458e10;
          budget.presn_momentum_code = p_dot_cgs * elapsed_years * physics::k_seconds_per_year /
              (bridge.mass_code_to_g * bridge.velocity_code_to_cm_s);
        }
        association_birth_years.push_back(
            record.birth_time_code * (m_units.length_si_per_code / m_units.velocity_si_per_code) / physics::k_seconds_per_year);
        if (energy_coupled && m_wind_table != nullptr) {
          // CHUI Part 3 wind/pre-SN channel: energy and momentum only (the
          // wind mass is part of the CC continuous return above).
          const physics::StellarWindIntervalBudget wind = m_wind_table->interval(age_b, age_e, z);
          budget.wind_energy_erg = wind.energy_erg_per_msun * birth_msun;
          budget.wind_momentum_code = wind.momentum_msun_kms_per_msun * birth_msun *
              (physics::k_solar_mass_g / bridge.mass_code_to_g) * (1.0e5 / bridge.velocity_code_to_cm_s);
          if (wind.metallicity_clamped) ++m_wind_metallicity_clamps;
        }
        staged_budgets.push_back(budget);
        StarUpdate update;
        update.star_index = s;
        update.particle_index = p;
        update.particle_id = pid;
        update.returned_mass_code = returned;
        update.age_end_years = age_e;
        update.energy_erg = e_sn;
        update.channel_mass = {scale * agb.mass_fraction * birth_mass, scale * cc.mass_fraction * birth_mass,
                               scale * ia.mass_fraction * n_ia * msun_code};
        update.channel_energy = {0.0, energy_coupled ? n_cc * m_ssp_model->config().energy_per_sn_erg : 0.0,
                                 energy_coupled ? n_ia * m_ssp_model->config().energy_per_sn_erg : 0.0};
        record.returned_mass_code += returned;
        update.record = record;
        updates.push_back(update);
      }
      } catch (...) {
        loop_failure = std::current_exception();
      }
      if (distributed_feedback) {
        m_mpi_context.rethrowCollectivePreparationFailure(loop_failure, "stellar population budgets");
      } else if (loop_failure != nullptr) {
        std::rethrow_exception(loop_failure);
      }
      if (m_physics_config.ssp_clustered_associations) {
        // RESEARCH: correlated associations / superbubbles. G0.3 ordering:
        // associations are discovered GLOBALLY and merged BEFORE any remote
        // candidate query, because merging moves the deposition point. The
        // merged budgets conserve mass, elements, dust, energy (including the
        // internal velocity-dispersion kinetic energy) and vector momentum.
        staged_budgets = mergeAssociationsDistributed(
            staged_budgets, association_birth_years, units, distributed_feedback,
            context.integrator_state.step_index);
      }
      if (distributed_feedback) {
        importRemoteFeedbackCandidates(state, deposition, staged_budgets, context.integrator_state.step_index);
      }
      std::exception_ptr stage_failure;
      try {
        deposition.reserveStagedRows(core::checkedSizeMultiply(
            staged_budgets.size(), static_cast<std::size_t>(fb_config.neighbor_count), "feedback staged row bound"));
        for (const physics::FeedbackSourceBudget& budget : staged_budgets) deposition.stage(state, budget);
      } catch (...) {
        stage_failure = std::current_exception();
      }
      if (distributed_feedback) {
        m_mpi_context.rethrowCollectivePreparationFailure(stage_failure, "mechanical feedback staging");
        // Patch 52 (P1-51-01/02): outbound delta preparation and the inbound
        // merge are guarded phases with collective votes; the inbound deltas
        // are received directly into a governed typed buffer whose lease
        // lives until they are merged.
        std::vector<physics::FeedbackCellDelta> outbound;
        std::vector<int> destinations;
        try {
          outbound = deposition.remoteDeltas(m_mpi_context.worldRank());
          destinations.reserve(outbound.size());
          for (const auto& d : outbound) destinations.push_back(d.owner_rank);
        } catch (...) {
          stage_failure = std::current_exception();
        }
        m_mpi_context.rethrowCollectivePreparationFailure(stage_failure, "mechanical feedback remote delta export");
        auto inbound = parallel::exchangeGovernedRecords<physics::FeedbackCellDelta>(
            m_mpi_context, outbound, destinations, m_memory_governor, "sources.feedback.remote_delta_exchange");
        try {
          // Deterministic merge order (cell id, then source rank).
          std::sort(inbound.records.begin(), inbound.records.end(), [](const auto& l, const auto& r) {
            return l.gas_cell_id != r.gas_cell_id ? l.gas_cell_id < r.gas_cell_id : l.source_rank < r.source_rank;
          });
          deposition.reserveStagedRows(inbound.records.size());
          for (const auto& d : inbound.records) deposition.addImportedDelta(state, d);
        } catch (...) {
          stage_failure = std::current_exception();
        }
        m_mpi_context.rethrowCollectivePreparationFailure(stage_failure, "mechanical feedback remote delta import");
      } else if (stage_failure != nullptr) {
        std::rethrow_exception(stage_failure);
      }
    });
    if (!transaction.validate(true, core::StageRejectionKind::kNone, {})) {
      m_ssp_model->table_clamps = yield_clamps_before;
      m_ssp_model->negative_ejecta_clamps = negative_clamps_before;
      m_ssp_model->untracked_mass_residual_msun = untracked_residual_before;
      m_ssp_model->nonmonotone_return_clamps = nonmonotone_before;
      throw core::StageTransactionRejected(transaction.outcome().rejection,
          "stellar population stage rejected: " + transaction.outcome().reason);
    }
    std::vector<std::uint32_t> touched;
    (void)transaction.commit([&] {
      touched = deposition.commit(state);
      for (const StarUpdate& u : updates) {
        state.particles.mass_code[u.particle_index] =
            std::max(state.particles.mass_code[u.particle_index] - u.returned_mass_code, 0.0);
        state.star_particles.stellar_age_years_last[u.star_index] = u.age_end_years;
        state.star_particles.stellar_returned_mass_cumulative_code[u.star_index] += u.returned_mass_code;
        state.star_particles.stellar_feedback_energy_cumulative_erg[u.star_index] += u.energy_erg;
        for (std::size_t c = 0; c < 3U; ++c) {
          state.star_particles.stellar_returned_mass_channel_cumulative_code[c][u.star_index] += u.channel_mass[c];
          state.star_particles.stellar_feedback_energy_channel_cumulative_erg[c][u.star_index] += u.channel_energy[c];
        }
        store.put(u.particle_id, u.record);  // event ordinals consumed only here
      }
      store.store(state);
      if (!updates.empty()) state.notifyParticleMassMutation();
    });
    m_last_feedback_ledger = deposition.ledger();
    m_last_event_chronology = event_chronology;  // published only with the commit
    internal::synchronizeParentParticleCompatibilityMirrors(state, m_world_rank, "SourceRuntime Part 2 stellar population");
    if (!touched.empty()) {
      state.notifyGasThermodynamicMutation();
      std::unique_ptr<physics::EffectiveIsmThermodynamicClosure> closure;
      if (m_effective_eos_table != nullptr) {
        closure = std::make_unique<physics::EffectiveIsmThermodynamicClosure>(m_effective_eos_table);
      }
      const hydro::GasThermodynamicRefreshPolicy policy{
          .adiabatic_index = 5.0 / 3.0,
          .density_floor_code = m_hydro_density_floor_code,
          .pressure_floor_code = m_hydro_pressure_floor_code,
          .scale_factor = source_evaluation_scale_factor,
          .redshift = m_is_cosmological ? std::max(0.0, 1.0 / source_evaluation_scale_factor - 1.0) : 0.0,
          .closure = closure.get(),
          .temperature_k_per_internal_energy_code =
              hydro::temperatureKelvinPerInternalEnergyCode(5.0 / 3.0, m_units.velocity_si_per_code)};
      std::vector<hydro::GasThermodynamicTouch> touches;
      for (const std::uint32_t row : touched) {
        touches.push_back({.cell_row = row, .internal_energy_before_code = state.gas_cells.internal_energy_code[row]});
      }
      (void)hydro::refreshGasThermodynamics(state, touches, policy);
    }
    if (context.profiler_session != nullptr) {
      const physics::FeedbackLedger& l = m_last_feedback_ledger;
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "sources.stellar_population",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "physics.stellar_feedback.mechanical",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "Part 2 SSP return and mechanical feedback ledger",
          .payload = {{"sources", std::to_string(l.sources)},
                      {"newborn_skipped_f5", std::to_string(newborn_skipped)},
                      {"legacy_initialized_records", std::to_string(legacy_initialized)},
                      {"mass_injected_code", std::to_string(l.mass_injected_code)},
                      {"energy_budget_code", std::to_string(l.energy_budget_code)},
                      {"kinetic_energy_injected_code", std::to_string(l.kinetic_energy_injected_code)},
                      {"thermal_energy_injected_code", std::to_string(l.thermal_energy_injected_code)},
                      {"sgs_energy_injected_code", std::to_string(l.sgs_energy_injected_code)},
                      {"radial_momentum_requested_code", std::to_string(l.radial_momentum_requested_code)},
                      {"radial_momentum_delivered_code", std::to_string(l.radial_momentum_delivered_code)},
                      {"net_vector_momentum_residual_code", std::to_string(l.net_vector_momentum_residual_code)},
                      {"resolved_cells", std::to_string(l.resolved_cells)},
                      {"unresolved_cells", std::to_string(l.unresolved_cells)},
                      {"dust_injected_code", std::to_string(l.dust_injected_code)},
                      {"dust_destroyed_code", std::to_string(l.dust_destroyed_code)},
                      {"yield_table_clamps", std::to_string(m_ssp_model->table_clamps)},
                      {"ssp_nonmonotone_return_clamps", std::to_string(m_ssp_model->nonmonotone_return_clamps)},
                      {"event_publication_time_code", std::to_string(m_last_event_chronology.interval_end_code)},
                      {"timed_discrete_events", std::to_string(m_last_event_chronology.timed_events)},
                      {"catch_up_discrete_events", std::to_string(m_last_event_chronology.catch_up_events)},
                      {"earliest_event_time_code",
                       std::to_string(m_last_event_chronology.earliestEventTimeCode())},
                      {"max_event_publication_lag_years",
                       std::to_string(m_last_event_chronology.max_publication_lag_years)},
                      {"sum_event_publication_lag_years",
                       std::to_string(m_last_event_chronology.sum_publication_lag_years)},
                      {"event_publication_lag_bound_years",
                       std::to_string(m_last_event_chronology.lag_bound_years)},
                      {"latest_event_time_code", std::to_string(m_last_event_chronology.latestEventTimeCode())},
                      {"unforecast_timed_discrete_events",
                       std::to_string(m_last_event_chronology.unforecast_timed_events)},
                      {"event_lag_policy", m_last_event_chronology.lagPolicy()},
                      {"event_deposition_semantics",
                       m_last_event_chronology.physical_event_barrier
                           ? "deposited at the synchronized interval end, which the time authority placed at the "
                             "next global physical event time (barrier resolution event_publication_lag_bound_years); "
                             "unforecast legacy-record events are counted and exempt"
                           : "deposited at the synchronized interval end (publication epoch), not at the physical "
                             "event time; physical times are model-defined threshold crossings"}}});
    }
  }

  void preparePopulationGrowth(core::SimulationState& state, std::size_t count,
                               core::ParticleSpecies species) {
    auto* scheduler = m_growth_scheduler;
      std::exception_ptr local_failure;
      try {
        if (count != 0U) {
          const std::size_t particle_target = core::checkedLocalCountAdd(
              state.particles.size(), count, core::kMaxLocalParticleCount,
              "particle", "source population growth");
          const std::size_t species_target = core::checkedSizeAdd(
              species == core::ParticleSpecies::kStar ? state.star_particles.size()
                  : state.black_holes.size(), count, "source species growth");
          const auto baseline = [&]() {
            std::uint64_t bytes = core::memoryReportBaselineOwnedBytes(
                core::collectSimulationMemoryReport(state));
            if (scheduler != nullptr) bytes = core::checkedMemoryBytesAdd(
                bytes, scheduler->ownedCapacityBytes(), "source scheduler baseline");
            return bytes;
          };
          core::RetainedCapacityTransaction plan(baseline);
          internal::planParticlePopulationGrowth(plan, state, particle_target,
                                                 species_target, species);
          if (scheduler != nullptr) {
            scheduler->planAppendCapacity(plan,
                core::checkedIntegralNarrow<std::uint32_t>(count,
                    "source scheduler growth count"), 0U);
          }
          plan.execute(m_memory_governor, core::MemoryClass::kCanonicalPersistent,
                       "sources.population_growth");
        }
      } catch (...) {
        local_failure = std::current_exception();
      }
      if (m_runtime_services != nullptr) {
        FailureCoordinator(*m_runtime_services).rethrowCollectiveFailure(
            local_failure, "source population growth admission");
      } else if (local_failure != nullptr) {
        std::rethrow_exception(local_failure);
      }
  }

  static void preparePopulationGrowthCallback(
      void* owner, core::SimulationState& state, std::size_t count,
      core::ParticleSpecies species) {
    static_cast<SourceRuntimeImpl*>(owner)->preparePopulationGrowth(state, count, species);
  }

  void buildOwnedLeafCellMetadata(const core::StepContext& context) {
    const core::SimulationState& state = context.state;
    const std::size_t cell_count = state.cells.size();
    m_owned_leaf_mask.assign(cell_count, 0U);
    core::GovernedScratchArena parent_arena(
        m_memory_governor, core::MemoryClass::kScratchArena,
        sourceParentIndexBytes(state.patches.size()),
        "sources.metal_diffusion.parent_index");
    const auto patch_ids_with_children = buildSourceParentIndex(
        state, parent_arena.resource());
    for (std::uint32_t cell_index = 0U; cell_index < cell_count; ++cell_index) {
      const PatchCellGeometry geometry = starFormationPatchCellGeometry(state, cell_index);
      bool owned_leaf = true;
      if (geometry.valid) {
        owned_leaf = state.patches.owning_rank[geometry.patch_index] == m_world_rank &&
            starFormationPatchIsLeaf(
                state, geometry.patch_index, patch_ids_with_children);
      }
      m_owned_leaf_mask[cell_index] = owned_leaf ? 1U : 0U;
    }
  }

  void buildActiveStarRows(
      const core::StepContext& context,
      std::pmr::vector<std::uint32_t>& active_star_rows,
      std::pmr::memory_resource* resource) {
    const core::SimulationState& state = context.state;
    active_star_rows.clear();
    if (!context.active_set.particles_are_subset ||
        context.active_set.particle_indices.empty()) return;
    std::pmr::vector<std::uint32_t> active_particles(resource);
    active_particles.assign(context.active_set.particle_indices.begin(),
                            context.active_set.particle_indices.end());
    std::sort(active_particles.begin(), active_particles.end());
    active_particles.erase(std::unique(active_particles.begin(), active_particles.end()),
                           active_particles.end());
    active_star_rows.reserve(state.star_particles.size());
    // Preserve canonical star-row ordering independently of scheduler order.
    for (std::uint32_t star_index = 0;
         star_index < state.star_particles.size(); ++star_index) {
      if (std::binary_search(active_particles.begin(), active_particles.end(),
              state.star_particles.particle_index[star_index])) {
        active_star_rows.push_back(star_index);
      }
    }
  }

  [[nodiscard]] double elapsedStellarEvolutionYears(
      const core::StepContext& context) const {
    constexpr double k_seconds_per_year = 31557600.0;
    double elapsed_si = context.timeline_step.dt_time_si;
    if (m_is_cosmological && context.cosmology_background != nullptr &&
        context.timeline_step.scale_factor_begin > 0.0 &&
        context.timeline_step.scale_factor_end >=
            context.timeline_step.scale_factor_begin) {
      elapsed_si = context.cosmology_background->cosmicTimeIntervalSi(
          context.timeline_step.scale_factor_begin,
          context.timeline_step.scale_factor_end);
    } else if (!(elapsed_si > 0.0)) {
      elapsed_si = context.integrator_state.dt_time_code * m_units.timeSiPerCode();
    }
    return std::max(elapsed_si / k_seconds_per_year, 0.0);
  }

  void executeStellarEvolutionAndEnrichment(
      core::StepContext& context, double source_evaluation_scale_factor) {
    if (!m_stellar_evolution.config().enabled) {
      return;
    }
    FeedbackThermodynamicTouchLog thermodynamic_touches(m_memory_governor);
    const bool all_stars_active =
        !context.active_set.particles_are_subset ||
        context.active_set.particle_indices.empty();
    std::unique_ptr<core::GovernedScratchArena> active_rows_arena;
    std::unique_ptr<std::pmr::vector<std::uint32_t>> active_star_rows;
    std::exception_ptr active_rows_failure;
    try {
      if (!all_stars_active) {
        const std::uint64_t bytes = core::checkedMemoryBytesAdd(
            static_cast<std::uint64_t>(core::checkedSizeMultiply(
                context.active_set.particle_indices.size(), sizeof(std::uint32_t),
                "source active particle lookup")),
            static_cast<std::uint64_t>(core::checkedSizeMultiply(
                context.state.star_particles.size(), sizeof(std::uint32_t),
                "source active star rows")), "source active rows");
        active_rows_arena = std::make_unique<core::GovernedScratchArena>(
            m_memory_governor, core::MemoryClass::kScratchArena,
            core::checkedMemoryBytesAdd(bytes, 256U, "source active alignment"),
            "sources.stellar_evolution.active_rows");
        active_star_rows = std::make_unique<std::pmr::vector<std::uint32_t>>(
            active_rows_arena->resource());
        buildActiveStarRows(context, *active_star_rows, active_rows_arena->resource());
      }
    } catch (...) { active_rows_failure = std::current_exception(); }
    m_mpi_context.rethrowCollectivePreparationFailure(
        active_rows_failure, "source active-star metadata preparation");
    const std::size_t active_star_count = all_stars_active
        ? context.state.star_particles.size()
        : active_star_rows->size();
    const std::size_t feedback_index_bytes_size = active_star_count == 0U
        ? 0U
        : core::checkedSizeMultiply(
              context.state.cells.size(),
              sizeof(std::uint32_t),
              "stellar-feedback spatial-index staging");
    std::size_t feedback_batch_max = std::min<std::size_t>(
        active_star_count, k_feedback_event_batch_max);
    if (m_memory_governor != nullptr && feedback_batch_max != 0U) {
      // Only new physical allocations consume headroom. Reused event/index
      // capacities remain charged to the governor; a vector replacement must
      // coexist with its old storage, so growth costs its full new capacity.
      const std::array<core::DeterministicBatchWorkspace, 3> workspaces{{
          {sizeof(physics::StellarEvolutionStarBudget), 0U},
          {sizeof(physics::StellarFeedbackEvent),
           static_cast<std::uint64_t>(m_feedback_events.capacity())},
          {all_stars_active ? sizeof(std::uint32_t) : 0U,
           static_cast<std::uint64_t>(m_contiguous_star_batch.capacity())},
      }};
      const core::DeterministicBatchSizingResult sizing =
          core::selectDeterministicBatchSizeForWorkspaces(
              m_memory_governor->snapshot(),
              core::DeterministicBatchSizingPolicy{
                  .requested_max_items = static_cast<std::uint64_t>(feedback_batch_max),
                  .bytes_per_item = 1U,
                  .fixed_reserve_bytes = static_cast<std::uint64_t>(feedback_index_bytes_size),
                  .minimum_items = 1U,
                  .alignment_items = 1U,
                  .headroom_use_basis_points = 10000U,
              }, workspaces);
      if (sizing.selected_items == 0U) {
        throw std::runtime_error(
            "stellar-feedback headroom cannot admit one evolution/event item after spatial-index staging");
      }
      feedback_batch_max = core::checkedIntegralNarrow<std::size_t>(
          sizing.selected_items, "stellar-feedback selected batch size");
    }
    const std::size_t required_event_capacity = feedback_batch_max;
    const std::size_t event_bytes_size =
        required_event_capacity > m_feedback_events.capacity()
        ? core::checkedSizeMultiply(
              required_event_capacity,
              sizeof(physics::StellarFeedbackEvent),
              "stellar-feedback event batch")
        : 0U;
    const std::size_t contiguous_star_bytes_size =
        required_event_capacity > m_contiguous_star_batch.capacity()
        ? core::checkedSizeMultiply(
              required_event_capacity, sizeof(std::uint32_t),
              "stellar-feedback contiguous star batch")
        : 0U;
    core::MemoryReservation feedback_index_rebuild_reservation;
    core::MemoryReservation replacement_event_reservation;
    core::MemoryReservation replacement_contiguous_star_reservation;
    std::exception_ptr reservation_failure;
    try {
      if (m_memory_governor != nullptr && feedback_index_bytes_size != 0U) {
        feedback_index_rebuild_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            static_cast<std::uint64_t>(feedback_index_bytes_size),
            "sources.stellar_feedback.spatial_index_rebuild");
      }
      if (m_memory_governor != nullptr && event_bytes_size != 0U) {
        replacement_event_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            static_cast<std::uint64_t>(event_bytes_size),
            "sources.stellar_feedback.event_batch");
      }
      if (m_memory_governor != nullptr && contiguous_star_bytes_size != 0U) {
        replacement_contiguous_star_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            static_cast<std::uint64_t>(contiguous_star_bytes_size),
            "sources.stellar_feedback.contiguous_star_batch");
      }
    } catch (...) {
      reservation_failure = std::current_exception();
    }
    if (m_runtime_services != nullptr) {
      FailureCoordinator(*m_runtime_services).rethrowCollectiveFailure(
          reservation_failure,
          "stellar-feedback memory preflight");
    } else if (reservation_failure != nullptr) {
      std::rethrow_exception(reservation_failure);
    }
    if (active_star_count == 0U) {
      return;
    }
    if (feedback_index_rebuild_reservation.pending()) {
      feedback_index_rebuild_reservation.commit();
    }
    if (replacement_event_reservation.pending()) {
      replacement_event_reservation.commit();
    }
    if (replacement_contiguous_star_reservation.pending()) {
      replacement_contiguous_star_reservation.commit();
    }
    const double elapsed_years = elapsedStellarEvolutionYears(context);

    const std::uint64_t leaf_bytes = sourceParentIndexBytes(context.state.patches.size());
    core::GovernedScratchArena leaf_arena(
        m_memory_governor, core::MemoryClass::kScratchArena,
        leaf_bytes, "sources.stellar_feedback.leaf_metadata");
    const auto patch_ids_with_children = buildSourceParentIndex(
        context.state, leaf_arena.resource());
    std::vector<std::uint32_t> owned_leaf_cells;
    owned_leaf_cells.reserve(context.state.cells.size());
    for (std::uint32_t cell_index = 0;
         cell_index < context.state.cells.size(); ++cell_index) {
      const PatchCellGeometry geometry =
          starFormationPatchCellGeometry(context.state, cell_index);
      bool owned_leaf = true;
      if (geometry.valid) {
        owned_leaf =
            context.state.patches.owning_rank[geometry.patch_index] == m_world_rank &&
            starFormationPatchIsLeaf(
                context.state, geometry.patch_index, patch_ids_with_children);
      }
      if (owned_leaf) {
        owned_leaf_cells.push_back(cell_index);
      }
    }

    const physics::StellarFeedbackGeometryView geometry_view{
        .particle_position_x_comoving =
            context.state.particles.position_x_comoving,
        .particle_position_y_comoving =
            context.state.particles.position_y_comoving,
        .particle_position_z_comoving =
            context.state.particles.position_z_comoving,
        .cell_center_x_comoving = context.state.cells.center_x_comoving,
        .cell_center_y_comoving = context.state.cells.center_y_comoving,
        .cell_center_z_comoving = context.state.cells.center_z_comoving,
        .gas_cell_id = context.state.gas_cells.gas_cell_id,
    };
    m_feedback_spatial_index.rebuildFromCellIndices(
        geometry_view, std::move(owned_leaf_cells));
    if (m_memory_governor != nullptr) {
      // The moved candidate allocation has become the retained index and the
      // old index buffer has been released by vector move-assignment. Transfer
      // the committed rebuild reservation only after that ownership change.
      m_feedback_index_reservation.release();
      m_feedback_index_reservation = std::move(feedback_index_rebuild_reservation);
    }
    FeedbackCellVolumeProvider volume_provider(context.state);
    physics::StellarFeedbackDepositionView deposition_view{
        .cell_mass_code = context.state.cells.mass_code,
        .gas_density_code = context.state.gas_cells.density_code,
        .gas_internal_energy_code =
            context.state.gas_cells.internal_energy_code,
        .gas_metal_mass_code = context.state.gas_cells.metal_mass_code,
        .cell_volume_provider = &volume_provider,
        .thermodynamic_touch_sink = &thermodynamic_touches,
    };
    m_feedback_events.clear();
    if (required_event_capacity > m_feedback_events.capacity()) {
      replaceSourceVectorCapacity(m_feedback_events, required_event_capacity,
                                  "sources.stellar_feedback.event_batch");
      if (m_memory_governor != nullptr) {
        m_feedback_event_reservation.release();
        m_feedback_event_reservation = std::move(replacement_event_reservation);
      }
    }
    if (required_event_capacity > m_contiguous_star_batch.capacity()) {
      replaceSourceVectorCapacity(m_contiguous_star_batch, required_event_capacity,
                                  "sources.stellar_feedback.contiguous_star_batch");
      if (m_memory_governor != nullptr) {
        m_contiguous_star_batch_reservation.release();
        m_contiguous_star_batch_reservation =
            std::move(replacement_contiguous_star_reservation);
      }
    }
    physics::StellarEvolutionBatchWorkspace evolution_workspace;
    for (std::size_t batch_begin = 0U;
         batch_begin < active_star_count;
         batch_begin += feedback_batch_max) {
      physics::StellarFeedbackStepReport feedback_step_report;
      const std::size_t batch_size = std::min<std::size_t>(
          feedback_batch_max, active_star_count - batch_begin);
      std::span<const std::uint32_t> star_batch;
      if (all_stars_active) {
        m_contiguous_star_batch.resize(batch_size);
        for (std::size_t i = 0U; i < batch_size; ++i) {
          m_contiguous_star_batch[i] = core::checkedIntegralNarrow<std::uint32_t>(
              batch_begin + i,
              "stellar-feedback contiguous active star row");
        }
        star_batch = m_contiguous_star_batch;
      } else {
        star_batch = std::span<const std::uint32_t>(*active_star_rows)
            .subspan(batch_begin, batch_size);
      }
      // Gate-4 production: caller-owned reusable governed budget buffer.
      // evaluateElapsedYearsGoverned admits capacity BEFORE growth, fills
      // within admitted capacity, and reconciles actual capacity() for reuse.
      // The old uncontrolled evaluateElapsedYears path is no longer used by
      // production; it remains only for diagnostics/tests.
      const physics::StellarEvolutionStepCounters evolution_counters =
          m_stellar_evolution.evaluateElapsedYearsGoverned(
              context.state, star_batch, elapsed_years, m_memory_governor,
              evolution_workspace,
              // Audit F5: stars born by THIS source stage (birth tick = this
              // step's tick, t_birth = t_end) have zero age in the interval.
              context.integrator_state.step_index);
      m_feedback_events.clear();
      for (const physics::StellarEvolutionStarBudget& budget :
           evolution_workspace.budgets) {
        m_feedback_events.push_back(physics::StellarFeedbackEvent{
            .star_index = budget.star_index,
            .returned_mass_code = budget.interval.returned_mass_code,
            .returned_metals_code = budget.interval.returned_metals_code,
            .feedback_energy_erg = budget.interval.feedback_energy_erg,
        });
      }
      if (!m_feedback_events.empty()) {
        // Admit the touch-log capacity for this batch before any deposition:
        // each event deposits into at most neighbor_count target cells.
        thermodynamic_touches.reserveAdditional(core::checkedSizeMultiply(
            m_feedback_events.size(),
            static_cast<std::size_t>(m_stellar_feedback.config().neighbor_count),
            "feedback thermodynamic touch batch bound"));
        m_stellar_feedback.applyEventsCountersOnly(
            context.state,
            m_stellar_feedback_state,
            geometry_view,
            &m_feedback_spatial_index,
            deposition_view,
            m_feedback_events,
            context.integrator_state.dt_time_code,
            feedback_step_report);
      }

      // Each bounded event batch is deposited or durably carried before its
      // matching stellar-evolution ledger advances. No processed event history
      // survives beyond this batch. Commit moves (not copies) the governed
      // workspace budgets into the report so no second population-scale copy
      // coexists with the admitted buffer; ownership moves back afterwards to
      // preserve retained capacity for the next batch.
      physics::StellarEvolutionStepReport commit_report;
      commit_report.counters = evolution_counters;
      commit_report.budgets = std::move(evolution_workspace.budgets);
      m_stellar_evolution.commitBudgets(context.state, commit_report);
      evolution_workspace.budgets = std::move(commit_report.budgets);
    }
    internal::synchronizeParentParticleCompatibilityMirrors(
        context.state, m_world_rank,
        "SourceRuntime stellar-evolution enrichment batch");
    if (!thermodynamic_touches.empty()) {
      // Audit F4: feedback changed density / internal energy / metals. Mark
      // the derived caches stale and refresh them through the single
      // closure contract at the source epoch, so the next CFL decision (and
      // any output) sees the post-feedback signal speed.
      context.state.notifyGasThermodynamicMutation();
      std::unique_ptr<physics::EffectiveIsmThermodynamicClosure> closure;
      if (m_effective_eos_table != nullptr) {
        closure = std::make_unique<physics::EffectiveIsmThermodynamicClosure>(m_effective_eos_table);
      }
      const hydro::GasThermodynamicRefreshReport refresh_report = hydro::refreshGasThermodynamics(
          context.state,
          thermodynamic_touches.collapse(),
          hydro::GasThermodynamicRefreshPolicy{
              .adiabatic_index = 5.0 / 3.0,
              .density_floor_code = m_hydro_density_floor_code,
              .pressure_floor_code = m_hydro_pressure_floor_code,
              .scale_factor = source_evaluation_scale_factor,
              .redshift = std::max(0.0, 1.0 / source_evaluation_scale_factor - 1.0),
              .closure = closure.get(),
              .temperature_k_per_internal_energy_code =
                  hydro::temperatureKelvinPerInternalEnergyCode(5.0 / 3.0, m_units.velocity_si_per_code)});
      m_last_thermodynamic_refresh = refresh_report;
      if (context.profiler_session != nullptr) {
        context.profiler_session->recordEvent(core::RuntimeEvent{
            .event_kind = "sources.thermodynamic_refresh",
            .severity = core::RuntimeEventSeverity::kInfo,
            .subsystem = "sources.stellar_feedback",
            .step_index = context.integrator_state.step_index,
            .simulation_time_code = context.integrator_state.current_time_code,
            .scale_factor = source_evaluation_scale_factor,
            .message = "refreshed derived gas thermodynamics after feedback deposition (audit F4)",
            .payload = {
                {"refreshed_cells", std::to_string(refresh_report.refreshed_cells)},
                {"effective_closure_cells", std::to_string(refresh_report.effective_closure_cells)},
                {"pressure_floor_cells", std::to_string(refresh_report.pressure_floor_cells)},
                {"temperature_unscaled_cells", std::to_string(refresh_report.temperature_unscaled_cells)},
                {"max_sound_speed_increase_ratio",
                 std::to_string(refresh_report.max_sound_speed_increase_ratio)},
                {"max_sound_speed_code", std::to_string(refresh_report.max_sound_speed_code)}}});
      }
    }
  }

  // T02 causal multilevel face ownership (patch 51; replaces the patch-50
  // finest-leaf clock, which read a coarse neighbour that had already
  // completed its interval while the fine side was at a substep start).
  //
  // Owner of a two-sided face f between leaves a and b:
  //   owner level L(f) = min(level(a), level(b))  (the COARSER side).
  // Integration of f over one level-L interval [t_b, t_e] (dt = t_e - t_b):
  //   * same-level face (level(a) == level(b) == L): dt/2 in the level-L
  //     pre half at t_b and dt/2 in the level-L post half at t_e (both
  //     sides advance in the same level sweep, so they are contemporaneous);
  //   * coarse/fine face (max level > L): dt/2 in the level-L pre half at
  //     t_b (every finer level is still synchronized with L at t_b; the
  //     transfer precedes the level-L history-start capture, so finer
  //     children interpolate from an already-updated coarse boundary) and
  //     dt/2 in executeAmrLevelInterfaceSynchronization at t_e (after every
  //     child substep, the complete hydro reflux and the retirement of the
  //     level-L history, when all levels >= L are synchronized at t_e).
  // Hence each face is integrated on both of its sides for exactly the same
  // physical interval [t_b, t_e], with states read at the same physical
  // instant on both sides (equal/opposite transfers applied immediately: no
  // flux register or coarse/fine mismatch exists), in a palindromic
  // S(dt/2) [subtree advance] S(dt/2) split around the level-L subtree.
  // Never applied in the post half or in any finer-level half, so no face is
  // integrated twice. Remote faces use the same symmetric min/max on both
  // owners, so the distributed claim/authority handshake stays one-to-one.
  // One-sided (open/reflective) faces belong to their single cell's level.
  struct LevelFaceSelection {
    std::int32_t level = 0;
    bool same_level = true;   // faces with min == max == level
    bool coarse_fine = false; // faces with min == level < max
    std::uint64_t phase_tag = 0U;  // exchange identity of this phase
  };
  struct LevelFaceSelectionCount {
    std::uint64_t same_level = 0U;
    std::uint64_t coarse_fine = 0U;
  };

  [[nodiscard]] static std::int32_t leafLevelOfRow(const core::SimulationState& state, std::uint32_t row) {
    if (row >= state.cells.patch_index.size()) return 0;
    const std::uint32_t patch = state.cells.patch_index[row];
    return patch < state.patches.level.size() ? state.patches.level[patch] : 0;
  }

  [[nodiscard]] static std::uint64_t levelPhaseTag(std::uint64_t step_index, int level, double time_begin_code,
                                                   std::uint32_t phase) noexcept {
    std::uint64_t bits = 0U;
    static_assert(sizeof(bits) == sizeof(time_begin_code));
    std::memcpy(&bits, &time_begin_code, sizeof(bits));
    std::uint64_t h = 1469598103934665603ULL;
    const auto mix = [&h](std::uint64_t v) { h ^= v; h *= 1099511628211ULL; };
    mix(step_index);
    mix(static_cast<std::uint64_t>(static_cast<std::uint32_t>(level)));
    mix(bits);
    mix(phase);
    return h == 0U ? 1U : h;
  }

  LevelFaceSelectionCount restrictDiffusionFacesToLevelSelection(
      const core::SimulationState& state,
      std::vector<physics::MetalDiffusionFace>& faces,
      std::vector<internal::RemoteDiffusionFace>& remote_faces,
      const LevelFaceSelection& selection) {
    LevelFaceSelectionCount count;
    const auto keep = [&](std::int32_t lo, std::int32_t hi) {
      if (lo != selection.level) return false;
      if (hi == lo) {
        if (selection.same_level) ++count.same_level;
        return selection.same_level;
      }
      if (selection.coarse_fine) ++count.coarse_fine;
      return selection.coarse_fine;
    };
    faces.erase(std::remove_if(faces.begin(), faces.end(),
                               [&](const physics::MetalDiffusionFace& face) {
                                 const std::int32_t left = leafLevelOfRow(state, face.left_cell);
                                 const bool two_sided =
                                     face.boundary_kind == physics::MetalDiffusionBoundaryKind::kInternal ||
                                     face.boundary_kind == physics::MetalDiffusionBoundaryKind::kPeriodic;
                                 const std::int32_t right = two_sided ? leafLevelOfRow(state, face.right_cell) : left;
                                 return !keep(std::min(left, right), std::max(left, right));
                               }),
                faces.end());
    remote_faces.erase(
        std::remove_if(remote_faces.begin(), remote_faces.end(),
                       [&](const internal::RemoteDiffusionFace& f) {
                         const std::int32_t local = leafLevelOfRow(state, f.local_row);
                         return !keep(std::min(local, f.remote_level), std::max(local, f.remote_level));
                       }),
        remote_faces.end());
    m_last_level_face_selection = count;
    // T02 (patch 52): exact identity ledger of the kept coarse/fine faces.
    // Built here (the remote face list is phase-local), admitted before it
    // is allocated; any local failure is parked in m_face_ledger_failure and
    // reported by the collective vote of the calling level phase (never
    // thrown between the diffusion/SGS collectives of this dispatch).
    std::vector<CanonicalInterfaceFace>().swap(m_last_coarse_fine_faces);
    m_last_coarse_fine_faces_lease.release();
    if (selection.coarse_fine && m_face_ledger_failure == nullptr) {
      try {
        const std::size_t bound = core::checkedSizeAdd(faces.size(), remote_faces.size(), "interface face ledger");
        if (m_memory_governor != nullptr && bound != 0U) {
          m_last_coarse_fine_faces_lease = m_memory_governor->reserve(
              core::MemoryClass::kPhaseResident,
              core::checkedSizeMultiply(bound, sizeof(CanonicalInterfaceFace), "interface face ledger bytes"),
              "sources.amr_interface_face_ledger");
          m_last_coarse_fine_faces_lease.commit();
        }
        m_last_coarse_fine_faces.reserve(bound);
        const auto gid = [&state](std::uint32_t row) {
          if (row >= state.gas_cells.gas_cell_id.size()) {
            throw std::runtime_error("interface face ledger: face row outside the gas-cell identity lane");
          }
          return state.gas_cells.gas_cell_id[row];
        };
        for (const physics::MetalDiffusionFace& face : faces) {
          const bool two_sided = face.boundary_kind == physics::MetalDiffusionBoundaryKind::kInternal ||
              face.boundary_kind == physics::MetalDiffusionBoundaryKind::kPeriodic;
          if (!two_sided) continue;  // one-sided faces belong to their cell's level (same-level)
          const std::int32_t left = leafLevelOfRow(state, face.left_cell);
          const std::int32_t right = leafLevelOfRow(state, face.right_cell);
          if (left == right) continue;  // same-level face
          const std::uint64_t l = gid(face.left_cell);
          const std::uint64_t r = gid(face.right_cell);
          m_last_coarse_fine_faces.push_back(CanonicalInterfaceFace{
              .lo_id = std::min(l, r), .hi_id = std::max(l, r), .axis = face.axis,
              .side = static_cast<std::int8_t>(l < r ? 1 : -1),
              .periodic = static_cast<std::uint8_t>(
                  face.boundary_kind == physics::MetalDiffusionBoundaryKind::kPeriodic ? 1U : 0U),
              .remote = 0U});
        }
        for (const internal::RemoteDiffusionFace& f : remote_faces) {
          if (leafLevelOfRow(state, f.local_row) == f.remote_level) continue;
          const bool local_is_lo = f.local_gas_cell_id < f.remote_gas_cell_id;
          m_last_coarse_fine_faces.push_back(CanonicalInterfaceFace{
              .lo_id = std::min(f.local_gas_cell_id, f.remote_gas_cell_id),
              .hi_id = std::max(f.local_gas_cell_id, f.remote_gas_cell_id), .axis = f.axis,
              .side = static_cast<std::int8_t>(local_is_lo ? f.local_side : -f.local_side),
              .periodic = static_cast<std::uint8_t>(f.periodic ? 1U : 0U), .remote = 1U});
        }
        std::sort(m_last_coarse_fine_faces.begin(), m_last_coarse_fine_faces.end());
        if (std::adjacent_find(m_last_coarse_fine_faces.begin(), m_last_coarse_fine_faces.end()) !=
            m_last_coarse_fine_faces.end()) {
          throw std::logic_error(
              "recursive AMR coarse/fine phase selected the same canonical interface face twice (double "
              "integration)");
        }
        if (m_last_coarse_fine_faces.size() != count.coarse_fine) {
          throw std::logic_error("interface face ledger size differs from the coarse/fine selection count");
        }
      } catch (...) {
        m_face_ledger_failure = std::current_exception();
      }
    }
    return count;
  }

  void executeDiffusionAndSgsShear(
      core::StepContext& context,
      double source_evaluation_scale_factor,
      double source_dt_code,
      std::optional<LevelFaceSelection> level_selection = std::nullopt) {
    m_last_level_face_selection = {};
    std::vector<CanonicalInterfaceFace>().swap(m_last_coarse_fine_faces);
    m_last_coarse_fine_faces_lease.release();
    const bool diffusion_enabled = m_metal_diffusion.config().enabled;
    // Configuration is rank-uniform. Do not derive collective participation
    // from local row count/sidecar population (empty ranks must still enter).
    const bool sgs_shear_enabled = m_source_stage_config.sgs.has_value();
    if (!diffusion_enabled && !sgs_shear_enabled) {
      return;
    }
    if (!(std::isfinite(source_dt_code) && source_dt_code > 0.0)) {
      throw std::invalid_argument(
          "smooth diffusion/SGS transport requires a finite positive source subinterval");
    }
    if (sgs_shear_enabled) {
      (void)physics::ensureGasCompositionInitialized(
          context.state, m_physics_config);
    }
    // G2.1: on more than one rank every rank enters every collective of the
    // distributed element-resolved operator, including ranks without gas.
    const bool distributed = m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1;
    if (!distributed && context.state.cells.size() == 0U) {
      return;
    }
    if (distributed && diffusion_enabled) {
      const std::uint64_t local_composition = context.state.gas_composition.enabled() ? 1U : 0U;
      const std::uint64_t any_composition = m_mpi_context.allreduceMaxUint64(local_composition);
      const std::uint64_t all_composition = m_mpi_context.allreduceMinUint64(local_composition);
      if (any_composition != all_composition) {
        throw std::runtime_error(
            "distributed metal diffusion: composition sidecar enabled on some ranks only (invalid state)");
      }
      if (any_composition == 0U) {
        throw std::runtime_error(
            "metal diffusion on multiple MPI ranks requires the element-resolved composition sidecar; the "
            "legacy total-Z scalar operator is single-rank only (refusing a decomposition-dependent result)");
      }
    }

    const std::size_t cell_count = context.state.cells.size();
    const std::uint64_t workspace_bytes = diffusion_enabled
        ? m_metal_diffusion.requiredWorkspaceBytes(cell_count)
        : 0U;
    const std::size_t mask_bytes = core::checkedSizeMultiply(
        cell_count, sizeof(std::uint8_t), "metal-diffusion owned-leaf mask");
    const std::size_t rho_kappa_bytes = core::checkedSizeMultiply(
        cell_count, sizeof(double), "metal-diffusion rho-kappa field");
    std::uint64_t phase_bytes = core::checkedMemoryBytesAdd(
        workspace_bytes,
        static_cast<std::uint64_t>(mask_bytes),
        "metal-diffusion phase workspace plus ownership mask");
    phase_bytes = core::checkedMemoryBytesAdd(
        phase_bytes,
        static_cast<std::uint64_t>(rho_kappa_bytes),
        "metal-diffusion phase plus rho-kappa field");
    if (sgs_shear_enabled) {
      phase_bytes = core::checkedMemoryBytesAdd(
          phase_bytes,
          physics::SgsShearTransferWorkspace::requiredBytes(cell_count),
          "smooth transport phase plus SGS shear staging");
    }

    core::MemoryReservation replacement_phase_reservation;
    std::exception_ptr reservation_failure;
    try {
      if (m_memory_governor != nullptr &&
          (m_owned_leaf_mask.capacity() < cell_count ||
           m_diffusion_rho_kappa_code.capacity() < cell_count ||
           m_diffusion_workspace.ownedCapacityBytes() < workspace_bytes ||
           (sgs_shear_enabled &&
            m_sgs_shear_workspace.ownedCapacityBytes() <
                physics::SgsShearTransferWorkspace::requiredBytes(cell_count)))) {
        replacement_phase_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            phase_bytes,
            "sources.metal_diffusion.cell_phase");
      }
    } catch (...) {
      reservation_failure = std::current_exception();
    }
    if (m_runtime_services != nullptr) {
      FailureCoordinator(*m_runtime_services).rethrowCollectiveFailure(
          reservation_failure,
          "metal-diffusion memory preflight");
    } else if (reservation_failure != nullptr) {
      std::rethrow_exception(reservation_failure);
    }
    if (replacement_phase_reservation.pending()) {
      replacement_phase_reservation.commit();
    }

    buildOwnedLeafCellMetadata(context);
    const double scale_factor = m_coordinate_frame == core::CoordinateFrame::kComoving
        ? source_evaluation_scale_factor
        : 1.0;
    // G2.1: remote owned leaf neighbours (one ghost layer) and the global
    // periodic envelope; collective on every rank.
    internal::RemoteDiffusionCellGather remote_gather;
    if (distributed) {
      remote_gather = internal::gatherRemoteDiffusionCells(
          context.state, m_mpi_context, m_memory_governor, scale_factor,
          m_hydro_boundary == core::BoundaryCondition::kPeriodic, m_box_extent_comoving,
          context.integrator_state.step_index);
    }
    // G0.5: admit the topology builder's complete construction working set
    // (new faces/remote faces/strain coexisting with the retained graph, the
    // boundary list, plane map, unique sets and directed refs) BEFORE it is
    // allocated. Released when this stage returns.
    core::MemoryReservation topology_reservation;
    internal::MetalDiffusionTopologyResult topology;
    std::exception_ptr topology_failure;
    try {
      if (m_memory_governor != nullptr) {
        topology_reservation = m_memory_governor->reserve(
            core::MemoryClass::kScratchArena,
            internal::metalDiffusionTopologyConstructionBoundBytes(
                context.state, m_owned_leaf_mask, m_world_rank, remote_gather.cells.size()),
            "sources.metal_diffusion.topology_construction");
        topology_reservation.commit();
      }
      topology = internal::buildMetalDiffusionTopology(
          context.state,
          m_owned_leaf_mask,
          m_world_rank,
          scale_factor,
          m_hydro_boundary,
          remote_gather.cells,
          distributed ? remote_gather.global_envelope_physical : std::nullopt);
    } catch (...) {
      topology_failure = std::current_exception();
    }
    if (distributed) {
      m_mpi_context.rethrowCollectivePreparationFailure(topology_failure, "metal-diffusion topology construction");
    } else if (topology_failure != nullptr) {
      std::rethrow_exception(topology_failure);
    }
    m_diffusion_topology_scratch_high_water_bytes = std::max(
        m_diffusion_topology_scratch_high_water_bytes,
        topology.construction_scratch_high_water_bytes);

    const std::uint64_t current_faces_capacity =
        core::ownedCapacityBytesForContainer(m_diffusion_faces);
    const std::uint64_t replacement_faces_capacity =
        core::ownedCapacityBytesForContainer(topology.faces);
    const std::uint64_t current_rho_kappa_capacity =
        core::ownedCapacityBytesForContainer(m_diffusion_rho_kappa_code);
    const std::uint64_t replacement_rho_kappa_capacity =
        core::ownedCapacityBytesForContainer(topology.strain_magnitude_code);
    m_diffusion_faces_high_water_bytes = std::max(
        m_diffusion_faces_high_water_bytes,
        std::max(current_faces_capacity, replacement_faces_capacity));
    m_diffusion_rho_kappa_high_water_bytes = std::max(
        m_diffusion_rho_kappa_high_water_bytes,
        std::max(current_rho_kappa_capacity, replacement_rho_kappa_capacity));

    // The retained diffusion face graph is scientifically meaningful O(Nface)
    // phase topology, but it still consumes real process memory. Admit the new
    // retained allocation before replacing the previous graph so the governor
    // accounts for old/new coexistence during the handoff. Construction-local
    // scratch remains separately visible through the measured scratch high-water.
    core::MemoryReservation replacement_faces_reservation;
    std::exception_ptr faces_reservation_failure;
    try {
      if (m_memory_governor != nullptr && replacement_faces_capacity != 0U) {
        replacement_faces_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            replacement_faces_capacity,
            "sources.metal_diffusion.faces");
      }
    } catch (...) {
      faces_reservation_failure = std::current_exception();
    }
    if (m_runtime_services != nullptr) {
      FailureCoordinator(*m_runtime_services).rethrowCollectiveFailure(
          faces_reservation_failure,
          "metal-diffusion topology residency preflight");
    } else if (faces_reservation_failure != nullptr) {
      std::rethrow_exception(faces_reservation_failure);
    }
    if (replacement_faces_reservation.pending()) {
      replacement_faces_reservation.commit();
    }

    m_diffusion_faces = std::move(topology.faces);
    m_diffusion_rho_kappa_code = std::move(topology.strain_magnitude_code);
    if (level_selection.has_value()) {
      // T02: this dispatch is one AMR level phase; only the faces owned by
      // this level and phase are integrated (restrictDiffusionFacesToLevelSelection).
      (void)restrictDiffusionFacesToLevelSelection(context.state, m_diffusion_faces, topology.remote_faces,
                                                   *level_selection);
    }
    if (m_memory_governor != nullptr) {
      m_diffusion_faces_reservation.release();
      if (replacement_faces_reservation.committed()) {
        m_diffusion_faces_reservation = std::move(replacement_faces_reservation);
      }
    }
    if (m_diffusion_rho_kappa_code.size() != cell_count) {
      throw std::runtime_error(
          "metal diffusion topology returned the wrong strain-field extent");
    }

    // I02: the same authoritative face graph now drives the resolved SGS
    // stress. Local faces receive equal/opposite momentum and stress-work
    // updates. For a remote face, both owners reconstruct the same canonical
    // lower/upper traction from exchanged step-start state and publish only
    // their local half. The existing distributed topology handshake is used
    // with zero scalar lanes solely to prove one-to-one interface geometry
    // before any rank commits a momentum/energy mutation.
    if (sgs_shear_enabled) {
      physics::SgsLedger shear_ledger;
      if (distributed) {
        internal::DistributedElementDiffusionCoupling topology_agreement(
            m_mpi_context, m_memory_governor, context.state,
            topology.remote_faces, m_diffusion_rho_kappa_code, 0U,
            level_selection.has_value() ? level_selection->phase_tag : 0U);
        topology_agreement.setup();
      }
      const internal::MpiStageCollectiveDecision decision(m_mpi_context);
      core::StageTransaction transaction(
          "sources.sgs_resolved_shear", m_memory_governor, &decision);
      (void)transaction.prepare([&] {
        shear_ledger = physics::prepareSgsShearTransfer(
            context.state,
            m_owned_leaf_mask,
            m_diffusion_rho_kappa_code,
            m_diffusion_faces,
            topology.remote_faces,
            scale_factor,
            source_dt_code,
            *m_source_stage_config.sgs,
            m_sgs_shear_workspace);
      });
      if (!transaction.validate(true, core::StageRejectionKind::kNone, {})) {
        m_sgs_shear_workspace.clear();
        throw core::StageTransactionRejected(
            transaction.outcome().rejection,
            "SGS resolved-shear stage rejected (" +
                std::string(core::stageRejectionKindLabel(
                    transaction.outcome().rejection)) +
                "): " + transaction.outcome().reason);
      }

      const bool touched_thermodynamics = std::any_of(
          m_sgs_shear_workspace.touched.begin(),
          m_sgs_shear_workspace.touched.end(),
          [](std::uint8_t value) { return value != 0U; });
      (void)transaction.commit([&] {
        physics::commitSgsShearTransfer(
            context.state, m_sgs_shear_workspace);
      });
      if (touched_thermodynamics) {
        context.state.notifyGasThermodynamicMutation();
        std::unique_ptr<physics::EffectiveIsmThermodynamicClosure> closure;
        if (m_effective_eos_table != nullptr) {
          closure = std::make_unique<physics::EffectiveIsmThermodynamicClosure>(
              m_effective_eos_table);
        }
        const hydro::GasThermodynamicRefreshPolicy refresh_policy{
            .adiabatic_index = 5.0 / 3.0,
            .density_floor_code = m_hydro_density_floor_code,
            .pressure_floor_code = m_hydro_pressure_floor_code,
            .scale_factor = source_evaluation_scale_factor,
            .redshift = m_is_cosmological
                ? std::max(0.0, 1.0 / source_evaluation_scale_factor - 1.0)
                : 0.0,
            .closure = closure.get(),
            .temperature_k_per_internal_energy_code =
                hydro::temperatureKelvinPerInternalEnergyCode(
                    5.0 / 3.0, m_units.velocity_si_per_code)};
        for (std::uint32_t row = 0U; row < cell_count; ++row) {
          if (m_sgs_shear_workspace.touched[row] != 0U) {
            (void)hydro::refreshGasThermodynamicCell(
                context.state, row, refresh_policy);
          }
        }
        context.state.markGasThermodynamicsRefreshed();
      }
      m_last_sgs_ledger.produced_from_resolved_code +=
          shear_ledger.produced_from_resolved_code;
      m_last_sgs_ledger.resolved_stress_work_code +=
          shear_ledger.resolved_stress_work_code;
      m_last_sgs_ledger.momentum_impulse_l1_code +=
          shear_ledger.momentum_impulse_l1_code;
      m_last_sgs_ledger.shear_faces += shear_ledger.shear_faces;
      m_last_sgs_ledger.shear_cells += shear_ledger.shear_cells;
      if (context.profiler_session != nullptr) {
        context.profiler_session->recordEvent(core::RuntimeEvent{
            .event_kind = "sources.sgs_shear",
            .severity = core::RuntimeEventSeverity::kInfo,
            .subsystem = "physics.sgs",
            .step_index = context.integrator_state.step_index,
            .simulation_time_code = context.integrator_state.current_time_code,
            .scale_factor = source_evaluation_scale_factor,
            .message = "conservative resolved-stress to SGS transfer",
            .payload = {
                {"faces", std::to_string(shear_ledger.shear_faces)},
                {"cells", std::to_string(shear_ledger.shear_cells)},
                {"produced_from_resolved_code",
                 std::to_string(shear_ledger.produced_from_resolved_code)},
                {"resolved_stress_work_code",
                 std::to_string(shear_ledger.resolved_stress_work_code)},
                {"momentum_impulse_l1_code",
                 std::to_string(shear_ledger.momentum_impulse_l1_code)},
                {"dt_code", std::to_string(source_dt_code)}}});
      }
    }

    if (!diffusion_enabled) {
      if (m_memory_governor != nullptr &&
          replacement_phase_reservation.committed()) {
        m_diffusion_phase_reservation.release();
        m_diffusion_phase_reservation =
            std::move(replacement_phase_reservation);
      }
      return;
    }

    for (std::uint32_t cell_index = 0U; cell_index < cell_count; ++cell_index) {
      if (m_owned_leaf_mask[cell_index] == 0U) {
        m_diffusion_rho_kappa_code[cell_index] = 0.0;
        continue;
      }
      const PatchCellGeometry geometry =
          starFormationPatchCellGeometry(context.state, cell_index);
      double volume_stored = context.state.gas_cells.density_code[cell_index] > 0.0
          ? context.state.cells.mass_code[cell_index] /
                context.state.gas_cells.density_code[cell_index]
          : 0.0;
      if (geometry.valid) {
        volume_stored = geometry.dx_stored * geometry.dy_stored *
            geometry.dz_stored;
      }
      const double volume_phys = volume_stored * scale_factor * scale_factor *
          scale_factor;
      const double gas_mass = context.state.cells.mass_code[cell_index];
      const double density = volume_phys > 0.0 ? gas_mass / volume_phys : 0.0;
      const double filter_length = volume_phys > 0.0
          ? std::cbrt(volume_phys) : 0.0;
      m_diffusion_rho_kappa_code[cell_index] =
          physics::smagorinskyRhoKappaCode(
              m_metal_diffusion.config(),
              density,
              filter_length,
              m_diffusion_rho_kappa_code[cell_index]);
    }

    if (m_diffusion_faces.empty() && !distributed) {
      if (m_memory_governor != nullptr &&
          replacement_phase_reservation.committed()) {
        m_diffusion_phase_reservation.release();
        m_diffusion_phase_reservation = std::move(replacement_phase_reservation);
      }
      return;
    }
    if (context.state.gas_composition.enabled()) {
      // CHUI Part 3 Checkpoint C: element-resolved conservative mixing of
      // every element/species/dust lane with the SAME diffusivity and face
      // graph (replaces the Part 2 refusal). The total-Z mirror is then
      // re-derived from the evolved lanes; it is never diffused on its own.
      const auto& sidecar = context.state.gas_composition;
      core::MemoryReservation element_reservation;
      std::exception_ptr element_failure;
      physics::ElementDiffusionReport element_report;
      // G0.5: elementDiffusionWorkspaceBytes is the single estimate of the
      // operator's scratch INCLUDING the per-subcycle delta buffer; the delta
      // is no longer added here (it was the only extra term, so the previous
      // sum double-counted it once the estimator became complete).
      try {
        if (m_memory_governor != nullptr) {
          element_reservation = m_memory_governor->reserve(
              core::MemoryClass::kScratchArena,
              physics::elementDiffusionWorkspaceBytes(sidecar.laneCount(), sidecar.rowCount()),
              "sources.element_diffusion.workspace");
          element_reservation.commit();
        }
      } catch (...) {
        element_failure = std::current_exception();
      }
      // A local refusal must be agreed before any rank enters the operator's
      // collectives (otherwise peers would block in the face exchange).
      if (distributed) {
        m_mpi_context.rethrowCollectivePreparationFailure(element_failure, "element diffusion admission");
      }
      if (element_failure == nullptr) {
        try {
          std::optional<internal::DistributedElementDiffusionCoupling> coupling;
          if (distributed) {
            coupling.emplace(m_mpi_context, m_memory_governor, context.state, topology.remote_faces,
                             m_diffusion_rho_kappa_code, sidecar.laneCount(),
                             level_selection.has_value() ? level_selection->phase_tag : 0U);
            coupling->setup();
          }
          element_report = physics::advanceElementResolvedDiffusion(
              context.state.gas_composition, m_diffusion_faces, m_diffusion_rho_kappa_code,
              source_dt_code, m_metal_diffusion.config().parabolic_cfl,
              m_metal_diffusion.config().max_subcycles, m_element_diffusion_workspace,
              coupling.has_value() ? &*coupling : nullptr, m_memory_governor);
        } catch (...) {
          element_failure = std::current_exception();
        }
      }
      m_element_diffusion_workspace.clear();
      m_element_diffusion_workspace.shrink_to_fit();
      if (m_runtime_services != nullptr) {
        FailureCoordinator(*m_runtime_services).rethrowCollectiveFailure(
            element_failure, "element-resolved metal diffusion");
      } else if (element_failure != nullptr) {
        std::rethrow_exception(element_failure);
      }
      (void)physics::synchronizeMetalMirrorFromComposition(context.state, std::span<const std::uint32_t>{});
      m_last_element_diffusion = element_report;
      if (m_memory_governor != nullptr && replacement_phase_reservation.committed()) {
        m_diffusion_phase_reservation.release();
        m_diffusion_phase_reservation = std::move(replacement_phase_reservation);
      }
      return;
    }
    const physics::MetalDiffusionStepReport report =
        m_metal_diffusion.advanceFromView(
            physics::MetalDiffusionFieldView{
                .gas_mass_code = context.state.cells.mass_code,
                .metal_mass_code = context.state.gas_cells.metal_mass_code,
                .rho_kappa_code = m_diffusion_rho_kappa_code,
            },
            m_diffusion_faces,
            source_dt_code,
            m_diffusion_workspace);
    if (m_memory_governor != nullptr &&
        replacement_phase_reservation.committed()) {
      m_diffusion_phase_reservation.release();
      m_diffusion_phase_reservation = std::move(replacement_phase_reservation);
    }

    std::ostringstream metadata;
    metadata << "module=metal_diffusion\n";
    metadata << "model=" << core::metalDiffusionModelToString(
        m_metal_diffusion.config().model) << "\n";
    metadata << "time_integrator=" <<
        core::metalDiffusionTimeIntegratorToString(
            m_metal_diffusion.config().time_integrator) << "\n";
    metadata << "stable_dt_code=" << report.stable_dt_code << "\n";
    metadata << "subcycles=" << report.subcycles << "\n";
    metadata << "rkl_stages=" << report.rkl_stages << "\n";
    metadata << "conservation_residual_code=" <<
        report.conservation_residual_code << "\n";
    const std::string text = metadata.str();
    core::ModuleSidecarBlock sidecar;
    sidecar.module_name = "metal_diffusion";
    sidecar.schema_version = 1U;
    sidecar.payload.resize(text.size());
    for (std::size_t i = 0U; i < text.size(); ++i) {
      sidecar.payload[i] = static_cast<std::byte>(text[i]);
    }
    context.state.sidecars.upsert(std::move(sidecar));
  }

  void buildStarFormationInputs(
      const core::StepContext& context,
      std::span<const std::uint32_t> active_cells,
      double source_evaluation_scale_factor,
      std::span<const std::uint64_t> patch_ids_with_children) {
    const core::SimulationState& state = context.state;
    m_star_formation_inputs.clear();

    const double scale_factor = source_evaluation_scale_factor;
    const double length_to_physical =
        m_coordinate_frame == core::CoordinateFrame::kComoving ? scale_factor : 1.0;

    for (const std::uint32_t cell_index : active_cells) {
      physics::StarFormationCellInput input;
      input.cell_index = cell_index;
      input.owning_rank = m_world_rank;
      input.is_active = true;
      if (cell_index >= state.cells.size() || cell_index >= state.gas_cells.size()) {
        input.is_active = false;
        m_star_formation_inputs.push_back(input);
        continue;
      }

      input.gas_cell_id = state.gas_cells.gas_cell_id[cell_index];
      if (input.gas_cell_id == 0U) {
        input.gas_cell_id = state.gasCellIdForLocalRow(cell_index).value_or(0U);
      }
      input.gas_mass_code = state.cells.mass_code[cell_index];
      input.gas_density_code = state.gas_cells.density_code[cell_index];
      input.gas_temperature_k = state.gas_cells.temperature_code[cell_index];
      input.gas_sound_speed_code = state.gas_cells.sound_speed_code[cell_index];
      input.gas_specific_internal_energy_code =
          state.gas_cells.internal_energy_code[cell_index];
      input.is_cosmological = m_is_cosmological;
      if (m_is_cosmological && m_mean_baryon_density0_code > 0.0 && scale_factor > 0.0) {
        const double density_phys = m_coordinate_frame == core::CoordinateFrame::kComoving
            ? input.gas_density_code / (scale_factor * scale_factor * scale_factor)
            : input.gas_density_code;
        const double mean_density_phys = m_mean_baryon_density0_code /
            (scale_factor * scale_factor * scale_factor);
        input.baryon_overdensity = density_phys / std::max(mean_density_phys, 1.0e-30);
      }
      input.velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[cell_index];
      input.velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[cell_index];
      input.velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[cell_index];
      input.gas_metal_mass_code = state.gas_cells.metal_mass_code[cell_index];
      input.metallicity_mass_fraction = input.gas_mass_code > 0.0
          ? input.gas_metal_mass_code / input.gas_mass_code
          : 0.0;
      input.center_x_comoving = state.cells.center_x_comoving[cell_index];
      input.center_y_comoving = state.cells.center_y_comoving[cell_index];
      input.center_z_comoving = state.cells.center_z_comoving[cell_index];

      const PatchCellGeometry geometry = starFormationPatchCellGeometry(state, cell_index);
      if (geometry.valid) {
        input.cell_volume_code = geometry.dx_stored * geometry.dy_stored * geometry.dz_stored;
        input.is_owned = state.patches.owning_rank[geometry.patch_index] == m_world_rank;
        input.is_leaf = starFormationPatchIsLeaf(state, geometry.patch_index, patch_ids_with_children);
        input.is_ghost = !input.is_owned;

        // CHUI Part 2: geometric face neighbours across patches, periodic
        // images and coarse/fine interfaces; an unavailable neighbour makes
        // the gradient NaN (never zero), which the adaptive collapse model
        // treats as ineligible.
        std::array<std::array<double, 3>, 3> gradient{};
        bool gradient_valid = false;
        if (m_sf_neighbor_index != nullptr) {
          const auto g = m_sf_neighbor_index->velocityGradient(cell_index, length_to_physical);
          gradient = g.velocity;
          gradient_valid = g.valid;
          if (!g.valid) ++m_sf_gradient_unavailable_cells;
          if (g.used_coarser_neighbor) ++m_sf_gradient_coarse_fine_cells;
        }
        if (gradient_valid) {
          input.velocity_divergence_code =
              gradient[0][0] + gradient[1][1] + gradient[2][2];
          input.velocity_gradient_frobenius_sq_code = 0.0;
          for (const auto& row : gradient) {
            for (const double value : row) {
              input.velocity_gradient_frobenius_sq_code += value * value;
            }
          }
        } else {
          input.velocity_divergence_code = std::numeric_limits<double>::quiet_NaN();
          input.velocity_gradient_frobenius_sq_code =
              std::numeric_limits<double>::quiet_NaN();
        }
      } else {
        input.cell_volume_code = input.gas_density_code > 0.0
            ? input.gas_mass_code / input.gas_density_code
            : 0.0;
        input.is_owned = true;
        input.is_leaf = true;
        input.is_ghost = false;
        // Adaptive production requires a real patch gradient. Legacy mode keeps the
        // historical converging-flow value of zero for non-AMR compatibility states.
        if (m_star_formation.config().model ==
            core::StarFormationModelKind::kAdaptiveBoundJeans) {
          input.velocity_divergence_code = std::numeric_limits<double>::quiet_NaN();
          input.velocity_gradient_frobenius_sq_code =
              std::numeric_limits<double>::quiet_NaN();
        } else {
          input.velocity_divergence_code = 0.0;
          input.velocity_gradient_frobenius_sq_code = 0.0;
        }
      }
      m_star_formation_inputs.push_back(input);
    }
  }

  core::UnitSystem m_units;
  std::shared_ptr<const physics::EffectiveMultiphaseEosTable> m_effective_eos_table;
  physics::StarFormationModel m_star_formation;
  physics::StellarEvolutionBookkeeper m_stellar_evolution;
  physics::StellarFeedbackModel m_stellar_feedback;
  physics::MetalDiffusionModel m_metal_diffusion;
  physics::StellarFeedbackModuleState m_stellar_feedback_state;
  physics::BlackHoleAgnModel m_black_hole;
  core::PhysicsConfig m_physics_config{};
  physics::BaryonicDataProducts m_baryonic_data{};
  std::unique_ptr<physics::ThermochemistryOperator> m_thermochemistry;
  physics::SourceLedger m_last_thermochemistry_ledger{};
  std::shared_ptr<const physics::SspYieldTable> m_ssp_yield_table;
  std::unique_ptr<physics::SspModel> m_ssp_model;
  physics::FeedbackLedger m_last_feedback_ledger{};
  DiscreteEventChronologyLedger m_last_event_chronology{};
  bool m_data_assumptions_reported = false;
  std::unique_ptr<internal::StarFormationNeighborIndex> m_sf_neighbor_index;
  std::array<double, 3> m_box_extent_comoving{};
  std::uint64_t m_sf_gradient_unavailable_cells = 0U;
  std::uint64_t m_sf_gradient_coarse_fine_cells = 0U;
  // SimulationConfig has no default constructor; initialized from the
  // constructor's member-initializer list (declaration order preserved).
  core::SimulationConfig m_full_config;
  physics::IsmRegimeReport m_last_ism_regime_report{};
  physics::ChemistryAdaptiveReport m_last_chemistry_adaptive_report{};
  physics::DustUpdateLedger m_last_dust_ledger{};
  physics::SgsLedger m_last_sgs_ledger{};
  physics::ChemistryShadowReport m_last_shadow_report{};
  bool m_amr_level_source_coupling_enabled = false;
  // T01/T02 single authority for the open Strang source interval. Plain
  // value state: captured/restored by the whole-step rollback checkpoint.
  SourceIntervalAccumulator m_source_interval{};
  std::uint64_t m_last_skipped_effective_rows = 0U;
  physics::BaryonicSourceStageConfig m_source_stage_config{};
  struct ControllerCheckpoint {
    std::uint64_t ssp_table_clamps = 0U;
    std::uint64_t ssp_negative_ejecta_clamps = 0U;
    double ssp_untracked_mass_residual_msun = 0.0;
    std::uint64_t ssp_nonmonotone_return_clamps = 0U;
    std::uint64_t wind_metallicity_clamps = 0U;
    std::uint64_t sf_gradient_unavailable_cells = 0U;
    std::uint64_t sf_gradient_coarse_fine_cells = 0U;
    physics::SourceLedger last_thermochemistry_ledger{};
    physics::FeedbackLedger last_feedback_ledger{};
    DiscreteEventChronologyLedger last_event_chronology{};
    physics::IsmRegimeReport last_ism_regime_report{};
    physics::ChemistryAdaptiveReport last_chemistry_adaptive_report{};
    physics::DustUpdateLedger last_dust_ledger{};
    physics::SgsLedger last_sgs_ledger{};
    physics::ChemistryShadowReport last_shadow_report{};
    SourceIntervalAccumulator source_interval{};
    std::uint64_t last_skipped_effective_rows = 0U;
    physics::ElementDiffusionReport last_element_diffusion{};
  };
  std::optional<ControllerCheckpoint> m_rollback_checkpoint;
  DistributedParticleIdRegistry m_particle_id_registry;
  const parallel::MpiContext& m_mpi_context;
  const RuntimeServices* m_runtime_services = nullptr;
  core::MemoryGovernor* m_memory_governor = nullptr;
  core::HierarchicalTimeBinScheduler* m_growth_scheduler = nullptr;
  core::MemoryReservation m_feedback_index_reservation;
  core::MemoryReservation m_feedback_event_reservation;
  core::MemoryReservation m_contiguous_star_batch_reservation;
  core::MemoryReservation m_star_formation_input_reservation;
  core::MemoryReservation m_contiguous_cell_batch_reservation;
  core::MemoryReservation m_diffusion_phase_reservation;
  core::MemoryReservation m_diffusion_faces_reservation;
  // Face counts selected by the last level-phase transport dispatch (local;
  // used for the coarse/fine pre/sync pairing check).
  LevelFaceSelectionCount m_last_level_face_selection{};
  // T02 (patch 52) face-identity ledger plumbing (see CanonicalInterfaceFace).
  // Each lease is declared before the storage it covers.
  core::MemoryReservation m_last_coarse_fine_faces_lease;
  std::vector<CanonicalInterfaceFace> m_last_coarse_fine_faces;
  std::exception_ptr m_face_ledger_failure;
  std::vector<core::MemoryReservation> m_interface_face_ledger_leases;  // per level, covers pending_interface_faces
  std::uint32_t m_world_rank = 0;
  core::CoordinateFrame m_coordinate_frame = core::CoordinateFrame::kComoving;
  core::BoundaryCondition m_hydro_boundary = core::BoundaryCondition::kOpen;
  bool m_bh_enabled = false;
  bool m_bh_seeding_requested = false;
  bool m_is_cosmological = false;
  double m_mean_baryon_density0_code = 0.0;
  double m_hydro_density_floor_code = 0.0;
  double m_hydro_pressure_floor_code = 0.0;
  hydro::GasThermodynamicRefreshReport m_last_thermodynamic_refresh{};
  std::vector<std::uint32_t> m_contiguous_cell_batch;
  std::vector<physics::StarFormationCellInput> m_star_formation_inputs;
  std::vector<std::uint32_t> m_contiguous_star_batch;
  std::vector<physics::StellarFeedbackEvent> m_feedback_events;
  // Gate-4: step-local governed budget workspace is created per source
  // execution (not retained across stages) so the retained budget capacity
  // does not compete with the hydro sweep workspace under the tight process
  // envelope. Within one execution the buffer is reused across batches with
  // admission-before-growth and actual-capacity reconciliation.
  physics::StellarFeedbackSpatialIndex m_feedback_spatial_index;
  std::vector<std::uint8_t> m_owned_leaf_mask;
  std::vector<double> m_diffusion_rho_kappa_code;
  std::vector<physics::MetalDiffusionFace> m_diffusion_faces;
  std::shared_ptr<const physics::StellarWindTable> m_wind_table;
  std::uint64_t m_wind_metallicity_clamps = 0U;
  std::vector<double> m_element_diffusion_workspace;
  physics::ElementDiffusionReport m_last_element_diffusion{};
  physics::MetalDiffusionWorkspace m_diffusion_workspace;
  physics::SgsShearTransferWorkspace m_sgs_shear_workspace;
  std::uint64_t m_diffusion_rho_kappa_high_water_bytes = 0U;
  std::uint64_t m_diffusion_faces_high_water_bytes = 0U;
  std::uint64_t m_diffusion_topology_scratch_high_water_bytes = 0U;
};

 }  // namespace

namespace internal {
std::unique_ptr<physics::ParticleIdPrecommit> makeDistributedParticleIdRegistry(
    const parallel::MpiContext& mpi_context, core::MemoryGovernor* governor) {
  auto registry = std::make_unique<DistributedParticleIdRegistry>(mpi_context);
  registry->setMemoryGovernor(governor);
  return registry;
}
}  // namespace internal

std::unique_ptr<SourceRuntime> makeSourceRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const core::UnitSystem& units,
    std::uint32_t world_rank,
    const parallel::MpiContext& mpi_context,
    std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table,
    const RuntimeServices* runtime_services) {
  return std::make_unique<SourceRuntimeImpl>(
      config,
      mode_policy,
      units,
      world_rank,
      mpi_context,
      std::move(effective_eos_table),
      runtime_services);
}

}  // namespace cosmosim::workflows
