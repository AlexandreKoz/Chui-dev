#include "cosmosim/workflows/source_runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <memory_resource>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/governed_scratch_arena.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/particle_id_registry.hpp"
#include "workflows/internal/metal_diffusion_topology.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/retained_population_growth.hpp"
#include "workflows/internal/star_formation_geometry.hpp"

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
  }

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
    const bool diffusion_governed = m_diffusion_phase_reservation.committed();
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.metal_diffusion.owned_leaf_mask",
        m_owned_leaf_mask,
        diffusion_governed);
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
        .governed_commitment = diffusion_governed,
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
        .governed_commitment = diffusion_governed,
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
        !m_metal_diffusion.config().enabled && !m_bh_enabled) {
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
    if (m_metal_diffusion.config().enabled && cell_count != 0U) {
      std::uint64_t diffusion = core::checkedMemoryBytesAdd(
          m_metal_diffusion.requiredWorkspaceBytes(state.cells.size()),
          multiply(cell_count, sizeof(std::uint8_t), "source diffusion mask"),
          "source diffusion phase");
      diffusion = core::checkedMemoryBytesAdd(diffusion,
          multiply(cell_count, sizeof(double), "source diffusion rho-kappa"),
          "source diffusion phase");
      peak = std::max(peak, diffusion);
    }
    return {peak, false,
            "source staging modeled; population growth, ID coordination and topology scratch require owner-local admission"};
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

    // SourceTerms executes after drift/hydro on the step-end state while the
    // committed IntegratorState remains at the step beginning until commitStep().
    // Make the physical evaluation epoch explicit so every source model sees
    // the state and cosmological conversion from the same instant.
    const double source_evaluation_scale_factor = m_is_cosmological
        ? context.timeline_step.scale_factor_end
        : 1.0;
    if (!std::isfinite(source_evaluation_scale_factor) ||
        source_evaluation_scale_factor <= 0.0) {
      throw std::runtime_error(
          "source runtime requires a finite positive source-stage scale factor");
    }

    // Only the current source stage may borrow the scheduler. The hook itself
    // was installed once at owner construction and performs no allocation here.
    m_growth_scheduler = context.particle_scheduler;
    struct GrowthSchedulerReset {
      core::HierarchicalTimeBinScheduler*& slot;
      ~GrowthSchedulerReset() { slot = nullptr; }
    } growth_scheduler_reset{m_growth_scheduler};

    const std::size_t cell_count = context.state.cells.size();
    const bool local_has_mutable_gravity_sources =
        cell_count > 0U || context.state.star_particles.size() > 0U ||
        context.state.black_holes.size() > 0U;
    const std::uint64_t gravity_source_mutation_rank_count =
        m_mpi_context.allreduceSumUint64(local_has_mutable_gravity_sources ? 1ULL : 0ULL);
    bool star_formation_mutated_gas = false;
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
        const physics::StarFormationStepReport batch_report =
            m_star_formation.applyFromInputs(
                context.state,
                m_star_formation_inputs,
                context.integrator_state.dt_time_code,
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

    executeStellarEvolutionAndEnrichment(context, source_evaluation_scale_factor);
    executeMetalDiffusion(context, source_evaluation_scale_factor);

    physics::BlackHoleAgnStepReport black_hole_report{};
    if (m_bh_enabled) {
      black_hole_report = m_black_hole.apply(
          context.state,
          {},
          context.integrator_state.dt_time_code,
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
      // signal speed (and the ideal-gas pressure cache). These models do not
      // expose touched-cell lists, so refresh every local gas cell once
      // (O(N_gas), only on steps where they removed gas). Temperature is
      // derived with the runtime-wide factor.
      context.state.notifyGasThermodynamicMutation();
      std::unique_ptr<physics::EffectiveIsmThermodynamicClosure> closure;
      if (m_effective_eos_table != nullptr) {
        closure = std::make_unique<physics::EffectiveIsmThermodynamicClosure>(m_effective_eos_table);
      }
      (void)hydro::refreshAllGasThermodynamics(
          context.state,
          hydro::GasThermodynamicRefreshPolicy{
              .adiabatic_index = 5.0 / 3.0,
              .density_floor_code = m_hydro_density_floor_code,
              .pressure_floor_code = m_hydro_pressure_floor_code,
              .scale_factor = source_evaluation_scale_factor,
              .redshift = std::max(0.0, 1.0 / source_evaluation_scale_factor - 1.0),
              .closure = closure.get(),
              .temperature_k_per_internal_energy_code =
                  hydro::temperatureKelvinPerInternalEnergyCode(5.0 / 3.0, m_units.velocity_si_per_code)});
    }
    if (gravity_source_mutation_rank_count > 0U) {
      // Source models can change gas/star/BH gravitating mass and membership.
      // Advance once per globally ordered source stage on every rank so PM/tree
      // validity follows authoritative state rather than PM refresh sequencing.
      context.state.bumpGravitySourceGeneration();
    }
  }

 private:
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
              evolution_workspace);
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

  void executeMetalDiffusion(
      core::StepContext& context,
      double source_evaluation_scale_factor) {
    if (!m_metal_diffusion.config().enabled || context.state.cells.size() == 0U) {
      return;
    }
    if (m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1) {
      throw std::runtime_error(
          "metal diffusion is not permitted on multiple MPI ranks until the directed gas-cell interface exchange "
          "and equal-and-opposite remote flux commit are available; refusing a decomposition-dependent result");
    }

    const std::size_t cell_count = context.state.cells.size();
    const std::uint64_t workspace_bytes =
        m_metal_diffusion.requiredWorkspaceBytes(cell_count);
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

    core::MemoryReservation replacement_phase_reservation;
    std::exception_ptr reservation_failure;
    try {
      if (m_memory_governor != nullptr &&
          (m_owned_leaf_mask.capacity() < cell_count ||
           m_diffusion_rho_kappa_code.capacity() < cell_count ||
           m_diffusion_workspace.ownedCapacityBytes() < workspace_bytes)) {
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
    internal::MetalDiffusionTopologyResult topology =
        internal::buildMetalDiffusionTopology(
            context.state,
            m_owned_leaf_mask,
            m_world_rank,
            scale_factor,
            m_hydro_boundary);
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

    if (m_diffusion_faces.empty()) {
      if (m_memory_governor != nullptr &&
          replacement_phase_reservation.committed()) {
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
            context.integrator_state.dt_time_code,
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

        const double dx_phys = geometry.dx_stored * length_to_physical;
        const double dy_phys = geometry.dy_stored * length_to_physical;
        const double dz_phys = geometry.dz_stored * length_to_physical;
        const std::array<std::span<const double>, 3> velocity_fields{
            state.gas_cells.velocity_x_peculiar,
            state.gas_cells.velocity_y_peculiar,
            state.gas_cells.velocity_z_peculiar};
        std::array<std::array<double, 3>, 3> gradient{};
        bool gradient_valid = true;
        for (std::size_t component = 0; component < velocity_fields.size(); ++component) {
          gradient[component][0] = starFormationDerivativeAtCell(
              velocity_fields[component], geometry, 0, dx_phys, cell_index);
          gradient[component][1] = starFormationDerivativeAtCell(
              velocity_fields[component], geometry, 1, dy_phys, cell_index);
          gradient[component][2] = starFormationDerivativeAtCell(
              velocity_fields[component], geometry, 2, dz_phys, cell_index);
          for (std::size_t axis = 0; axis < gradient[component].size(); ++axis) {
            gradient_valid = gradient_valid && std::isfinite(gradient[component][axis]);
          }
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
  physics::MetalDiffusionWorkspace m_diffusion_workspace;
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
