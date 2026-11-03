#include "cosmosim/io/restart_checkpoint.hpp"

#include <algorithm>
#include <bit>
#include <array>
#include <cmath>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <optional>
#include <memory_resource>
#include <sstream>
#include <streambuf>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <utility>
#include <unordered_set>
#include <span>
#include <type_traits>

#include "core/internal/sha256.hpp"
#include "cosmosim/core/governed_scratch_arena.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/io/io_contract.hpp"
#include "io/internal/transactional_file.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <fcntl.h>
#include <hdf5.h>
#include <unistd.h>
#endif

namespace cosmosim::io {
namespace {

constexpr std::uint64_t k_offset_basis = 14695981039346656037ull;
constexpr std::uint64_t k_prime = 1099511628211ull;
constexpr std::uint32_t k_restart_schema_v14 = 14;
constexpr std::uint32_t k_restart_schema_v15 = 15;
constexpr std::uint32_t k_restart_schema_v16 = 16;
constexpr std::uint32_t k_restart_schema_v17 = 17;
constexpr std::uint32_t k_restart_schema_v18 = 18;
constexpr std::uint32_t k_restart_schema_v19 = 19;
constexpr std::uint32_t k_restart_schema_v20 = 20;
constexpr std::uint32_t k_restart_schema_v21 = 21;
constexpr std::uint32_t k_restart_schema_v22 = 22;
constexpr std::uint32_t k_restart_schema_v23 = 23;
constexpr std::uint32_t k_restart_schema_v24 = 24;
constexpr std::string_view k_restart_schema_name_v14 = "cosmosim_restart_v14";
constexpr std::string_view k_restart_schema_name_v15 = "cosmosim_restart_v15";
constexpr std::string_view k_restart_schema_name_v16 = "cosmosim_restart_v16";
constexpr std::string_view k_restart_schema_name_v17 = "cosmosim_restart_v17";
constexpr std::string_view k_restart_schema_name_v18 = "cosmosim_restart_v18";
constexpr std::string_view k_restart_schema_name_v19 = "cosmosim_restart_v19";
constexpr std::string_view k_restart_schema_name_v20 = "cosmosim_restart_v20";
constexpr std::string_view k_restart_schema_name_v21 = "cosmosim_restart_v21";
constexpr std::string_view k_restart_schema_name_v22 = "cosmosim_restart_v22";
constexpr std::string_view k_restart_schema_name_v23 = "cosmosim_restart_v23";
constexpr std::string_view k_hierarchical_endpoint_group = "/gravity_force_cache/hierarchical_endpoint";
constexpr std::string_view k_gas_identity_row_policy = "explicit_dense_local_cell_row";
constexpr std::string_view k_gas_cell_scheduler_identity_key = "gas_cell_id";

[[nodiscard]] std::string hexU64(std::uint64_t value) {
  constexpr char k_hex[] = "0123456789abcdef";
  std::string text(16, '0');
  for (int i = 15; i >= 0; --i) {
    text[static_cast<std::size_t>(i)] = k_hex[value & 0xF];
    value >>= 4;
  }
  return text;
}

[[nodiscard]] std::uint64_t fnv1aAppend(std::uint64_t seed, std::span<const std::byte> bytes) {
  std::uint64_t hash = seed;
  for (const std::byte value : bytes) {
    hash ^= static_cast<std::uint64_t>(std::to_integer<unsigned char>(value));
    hash *= k_prime;
  }
  return hash;
}


struct RestartIntegrityDigests {
  std::uint64_t legacy_fnv1a = 0;
  std::string sha256_hex;
};

[[nodiscard]] std::string sha256DigestHex(const std::array<std::uint8_t, 32>& digest) {
  static constexpr char HEX[] = "0123456789abcdef";
  std::string out(64U, '0');
  for (std::size_t i = 0; i < digest.size(); ++i) {
    out[i * 2U] = HEX[digest[i] >> 4U];
    out[i * 2U + 1U] = HEX[digest[i] & 0x0fU];
  }
  return out;
}

template <typename T>
void appendCanonicalSha256(core::internal::Sha256& hash, T value) {
  using RawT = std::remove_cv_t<T>;
  static_assert(std::is_arithmetic_v<RawT> || std::is_enum_v<RawT>);
  if constexpr (std::is_floating_point_v<RawT>) {
    if constexpr (sizeof(RawT) == sizeof(std::uint64_t)) {
      appendCanonicalSha256(hash, std::bit_cast<std::uint64_t>(value));
    } else if constexpr (sizeof(RawT) == sizeof(std::uint32_t)) {
      appendCanonicalSha256(hash, std::bit_cast<std::uint32_t>(value));
    }
  } else if constexpr (std::is_enum_v<RawT>) {
    using ValueT = std::underlying_type_t<RawT>;
    using UnsignedT = std::make_unsigned_t<ValueT>;
    UnsignedT bits = static_cast<UnsignedT>(value);
    std::array<std::uint8_t, sizeof(UnsignedT)> bytes{};
    for (std::size_t i = 0; i < bytes.size(); ++i) {
      bytes[i] = static_cast<std::uint8_t>(bits & static_cast<UnsignedT>(0xffU));
      bits >>= 8U;
    }
    hash.update(bytes.data(), bytes.size());
  } else {
    using UnsignedT = std::make_unsigned_t<RawT>;
    UnsignedT bits = static_cast<UnsignedT>(value);
    std::array<std::uint8_t, sizeof(UnsignedT)> bytes{};
    for (std::size_t i = 0; i < bytes.size(); ++i) {
      bytes[i] = static_cast<std::uint8_t>(bits & static_cast<UnsignedT>(0xffU));
      bits >>= 8U;
    }
    hash.update(bytes.data(), bytes.size());
  }
}

template <typename VectorLike>
[[nodiscard]] std::span<const std::byte> asBytesSpan(const VectorLike& values) {
  using ValueType = typename VectorLike::value_type;
  return {reinterpret_cast<const std::byte*>(values.data()), values.size() * sizeof(ValueType)};
}

void validateSchedulerPersistentStateForRestart(
    const core::TimeBinPersistentState& scheduler_state,
    std::size_t expected_element_count,
    std::string_view context,
    std::string_view scheduler_path) {
  if (scheduler_state.bin_index.size() != scheduler_state.next_activation_tick.size() ||
      scheduler_state.bin_index.size() != scheduler_state.active_flag.size() ||
      scheduler_state.bin_index.size() != scheduler_state.pending_bin_index.size()) {
    throw std::invalid_argument(
        std::string(context) + ": " + std::string(scheduler_path) +
        " persistent arrays must have matching sizes");
  }
  if (scheduler_state.bin_index.size() != expected_element_count) {
    throw std::invalid_argument(
        std::string(context) + ": " + std::string(scheduler_path) +
        " element count does not match the authoritative state extent");
  }
  for (std::size_t i = 0; i < scheduler_state.bin_index.size(); ++i) {
    if (scheduler_state.bin_index[i] > scheduler_state.max_bin) {
      throw std::invalid_argument(
          std::string(context) + ": " + std::string(scheduler_path) + "/bin_index exceeds max_bin");
    }
    if (scheduler_state.active_flag[i] > 1U) {
      throw std::invalid_argument(
          std::string(context) + ": " + std::string(scheduler_path) + "/active_flag must contain only 0 or 1");
    }
    const bool pending_is_unset =
        scheduler_state.pending_bin_index[i] == core::HierarchicalTimeBinScheduler::k_unset_pending_bin;
    if (!pending_is_unset && scheduler_state.pending_bin_index[i] > scheduler_state.max_bin) {
      throw std::invalid_argument(
          std::string(context) + ": " + std::string(scheduler_path) + "/pending_bin_index exceeds max_bin");
    }
  }
}

void validateGravityForceCacheForRestart(
    const GravityForceCachePersistentView& cache,
    const core::SimulationState& state,
    std::string_view context) {
  const auto validate_triplet = [&](std::span<const double> x,
                                    std::span<const double> y,
                                    std::span<const double> z,
                                    std::size_t expected,
                                    std::string_view label) {
    if (x.size() != y.size() || x.size() != z.size()) {
      throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                  std::string(label) + " component extents differ");
    }
    if (cache.valid && x.size() != expected) {
      throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                  std::string(label) + " extent does not match authoritative state");
    }
    if (!cache.valid && !x.empty()) {
      throw std::invalid_argument(std::string(context) +
                                  ": invalid gravity force cache must not retain stale component lanes");
    }
    for (const double value : x) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache contains non-finite values");
      }
    }
    for (const double value : y) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache contains non-finite values");
      }
    }
    for (const double value : z) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache contains non-finite values");
      }
    }
  };
  const auto validate_identity = [&](std::span<const std::uint64_t> cache_ids,
                                     std::span<const std::uint64_t> state_ids,
                                     std::string_view label) {
    if (!cache.valid && !cache_ids.empty()) {
      throw std::invalid_argument(std::string(context) + ": invalid gravity force cache must not retain stale " +
                                  std::string(label) + " identity lanes");
    }
    if (!cache.valid) {
      return;
    }
    if (cache_ids.size() != state_ids.size()) {
      throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                  std::string(label) + " identity extent does not match authoritative state");
    }
    std::unordered_set<std::uint64_t> seen;
    seen.reserve(cache_ids.size());
    for (std::size_t index = 0; index < cache_ids.size(); ++index) {
      if (cache_ids[index] == 0U || !seen.insert(cache_ids[index]).second) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                    std::string(label) + " identity lane contains a zero or duplicate ID");
      }
      if (cache_ids[index] != state_ids[index]) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                    std::string(label) + " identity lane does not match authoritative dense row");
      }
    }
  };
  validate_triplet(cache.particle_accel_x_comoving, cache.particle_accel_y_comoving,
                   cache.particle_accel_z_comoving, state.particles.size(), "particles");
  validate_triplet(cache.cell_accel_x_comoving, cache.cell_accel_y_comoving,
                   cache.cell_accel_z_comoving, state.cells.size(), "gas cells");
  validate_identity(cache.particle_id, state.particle_sidecar.particle_id, "particle");
  validate_identity(cache.gas_cell_id, state.gas_cells.gas_cell_id, "gas-cell");
  // v24 hierarchical endpoint: all-or-nothing, rows follow particle_id.
  const std::array<std::span<const double>, 7> endpoint_lanes{
      cache.hierarchical_tree_accel_x_comoving, cache.hierarchical_tree_accel_y_comoving,
      cache.hierarchical_tree_accel_z_comoving, cache.hierarchical_pm_accel_x_comoving,
      cache.hierarchical_pm_accel_y_comoving, cache.hierarchical_pm_accel_z_comoving,
      cache.hierarchical_nearest_pair_r2_comoving};
  if (!cache.hierarchical_endpoint_valid) {
    for (const auto lane : endpoint_lanes) {
      if (!lane.empty()) {
        throw std::invalid_argument(std::string(context) +
                                    ": invalid hierarchical endpoint must not retain stale lanes");
      }
    }
    return;
  }
  if (!cache.valid) {
    throw std::invalid_argument(std::string(context) +
                                ": hierarchical endpoint requires a valid gravity force cache");
  }
  for (std::size_t lane = 0; lane < endpoint_lanes.size(); ++lane) {
    if (endpoint_lanes[lane].size() != state.particles.size()) {
      throw std::invalid_argument(std::string(context) +
                                  ": hierarchical endpoint lane extent does not match authoritative state");
    }
    for (const double value : endpoint_lanes[lane]) {
      // The nearest-pair lane is +infinity for a target without a direct pair.
      const bool legal = lane + 1U == endpoint_lanes.size()
          ? !std::isnan(value) && value >= 0.0
          : std::isfinite(value);
      if (!legal) {
        throw std::invalid_argument(std::string(context) +
                                    ": hierarchical endpoint lane contains an illegal value");
      }
    }
  }
}

void validateRestartTimeBinMirrorsAgainstScheduler(
    const core::SimulationState& state,
    const core::TimeBinPersistentState& scheduler_state,
    std::string_view context) {
  validateSchedulerPersistentStateForRestart(
      scheduler_state, state.particles.size(), context, "/scheduler");
  if (state.particles.size() != scheduler_state.bin_index.size()) {
    throw std::invalid_argument(
        std::string(context) + ": logical particle time_bin length must match /scheduler/bin_index length");
  }
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    if (state.particleTimeBin(i) != scheduler_state.bin_index[i]) {
      throw std::invalid_argument(
          std::string(context) + ": /state/particles/time_bin is stale relative to /scheduler/bin_index");
    }
  }
}

void validateGasCellTimeBinMirrorsAgainstScheduler(
    const core::SimulationState& state,
    const core::TimeBinPersistentState& scheduler_state,
    std::span<const std::uint64_t> gas_cell_ids,
    std::string_view context) {
  validateSchedulerPersistentStateForRestart(
      scheduler_state, state.cells.size(), context, "/gas_cell_scheduler");
  state.requireGasCellIdentityMapCoversDenseRows(context);
  if (state.cells.time_bin.size() != scheduler_state.bin_index.size()) {
    throw std::invalid_argument(
        std::string(context) + ": /state/cells/time_bin length must match /gas_cell_scheduler/bin_index length");
  }
  if (gas_cell_ids.size() != state.cells.size()) {
    throw std::invalid_argument(
        std::string(context) + ": /gas_cell_scheduler/gas_cell_id length must match CellSoa extent");
  }
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const core::GasCellIdentityRecord* identity = state.gas_cell_identity.findByLocalRow(cell_index);
    if (identity == nullptr) {
      throw std::invalid_argument(std::string(context) + ": gas-cell identity map is missing a dense row");
    }
    if (gas_cell_ids[cell_index] != identity->gas_cell_id) {
      throw std::invalid_argument(
          std::string(context) + ": /gas_cell_scheduler/gas_cell_id is stale relative to GasCellIdentityMap");
    }
    if (state.cells.time_bin[cell_index] != scheduler_state.bin_index[cell_index]) {
      throw std::invalid_argument(
          std::string(context) + ": /state/cells/time_bin is stale relative to /gas_cell_scheduler/bin_index");
    }
  }
}

void rebuildRestartTimeBinMirrorsFromScheduler(
    core::SimulationState& state,
    const core::TimeBinPersistentState& scheduler_state) {
  if (state.particles.time_bin.size() != scheduler_state.bin_index.size()) {
    throw std::invalid_argument("restart reader: particle time_bin mirror length must match scheduler bin_index length");
  }
  for (std::size_t i = 0; i < state.particles.time_bin.size(); ++i) {
    state.particles.time_bin[i] = scheduler_state.bin_index[i];
  }
}

void rebuildGasCellTimeBinMirrorsFromScheduler(
    core::SimulationState& state,
    const core::TimeBinPersistentState& scheduler_state,
    std::span<const std::uint64_t> gas_cell_ids) {
  validateGasCellTimeBinMirrorsAgainstScheduler(state, scheduler_state, gas_cell_ids, "restart reader");
  for (std::size_t i = 0; i < state.cells.time_bin.size(); ++i) {
    state.cells.time_bin[i] = scheduler_state.bin_index[i];
  }
}

[[nodiscard]] std::vector<std::uint64_t> gasCellIdsByDenseLocalRow(
    const core::SimulationState& state,
    std::string_view context) {
  state.requireGasCellIdentityMapCoversDenseRows(context);
  std::vector<std::uint64_t> gas_cell_ids(state.cells.size(), 0U);
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const core::GasCellIdentityRecord* identity = state.gas_cell_identity.findByLocalRow(cell_index);
    if (identity == nullptr) {
      throw std::runtime_error(std::string(context) + ": gas-cell identity map is missing a dense row");
    }
    gas_cell_ids[cell_index] = identity->gas_cell_id;
  }
  return gas_cell_ids;
}

[[nodiscard]] core::TimeBinPersistentState legacyGasCellSchedulerStateFromMirrors(
    const core::SimulationState& state,
    const core::TimeBinPersistentState& particle_scheduler_state) {
  // Compatibility only for older direct API callers.  The reference workflow
  // always supplies a dedicated gas-cell scheduler; v19 then persists that
  // authoritative state explicitly.  This branch never maps cells through
  // parent particles and therefore remains valid for parentless cells.
  core::TimeBinPersistentState gas_scheduler_state;
  gas_scheduler_state.current_tick = particle_scheduler_state.current_tick;
  gas_scheduler_state.max_bin = particle_scheduler_state.max_bin;
  gas_scheduler_state.bin_index.assign(state.cells.time_bin.begin(), state.cells.time_bin.end());
  gas_scheduler_state.next_activation_tick.resize(state.cells.size(), particle_scheduler_state.current_tick);
  gas_scheduler_state.active_flag.assign(state.cells.size(), 0U);
  gas_scheduler_state.pending_bin_index.assign(
      state.cells.size(), core::HierarchicalTimeBinScheduler::k_unset_pending_bin);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint8_t bin = gas_scheduler_state.bin_index[cell_index];
    if (bin > gas_scheduler_state.max_bin) {
      throw std::invalid_argument("restart legacy gas-cell scheduler mirror has a bin above max_bin");
    }
    const std::uint64_t interval = 1ULL << bin;
    gas_scheduler_state.next_activation_tick[cell_index] =
        ((particle_scheduler_state.current_tick / interval) + 1U) * interval;
  }
  return gas_scheduler_state;
}

void validateHydroGeometryStateForRestart(
    const core::SimulationState& state,
    std::string_view context) {
  if (state.cells.size() != state.gas_cells.size()) {
    throw std::invalid_argument(std::string(context) + ": /state/cells and /state/gas_cells sizes must match");
  }
  if (!state.cells.isConsistent()) {
    throw std::invalid_argument(std::string(context) + ": /state/cells lanes must have matching lengths");
  }
  if (!state.gas_cells.isConsistent()) {
    throw std::invalid_argument(std::string(context) + ": /state/gas_cells lanes must have matching lengths");
  }
  if (!state.patches.isConsistent()) {
    throw std::invalid_argument(std::string(context) + ": /state/patches lanes must have matching lengths");
  }
  if (state.cells.size() == 0) {
    if (state.patches.size() != 0) {
      throw std::invalid_argument(std::string(context) + ": /state/patches must be empty when no gas cells exist");
    }
    return;
  }
  // PatchSoa is mandatory for AMR-backed gas geometry, but legacy/non-AMR
  // Cartesian states may intentionally have no patch descriptors.  In that
  // case `patch_index == 0` and `owning_patch_id == 0` are sentinel metadata;
  // they are not a reason to reject a restart whose cell identity and sidecar
  // lanes are otherwise complete.  Do not manufacture a fake patch merely to
  // satisfy restart I/O.
  if (state.patches.size() == 0U) {
    for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
      if (state.cells.patch_index[cell_index] != 0U) {
        throw std::invalid_argument(
            std::string(context) + ": /state/cells/patch_index must be zero when PatchSoa is absent");
      }
    }
  } else {
    std::vector<std::uint8_t> cell_covered_by_patch(state.cells.size(), 0U);
    for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      const std::uint64_t first_cell = state.patches.first_cell[patch_index];
      const std::uint64_t cell_count = state.patches.cell_count[patch_index];
      if (cell_count == 0U) {
        throw std::invalid_argument(std::string(context) + ": /state/patches/cell_count entries must be non-zero");
      }
      if (first_cell > state.cells.size() || cell_count > state.cells.size() - first_cell) {
        throw std::invalid_argument(std::string(context) + ": /state/patches cell range exceeds /state/cells");
      }
      for (std::uint64_t offset = 0; offset < cell_count; ++offset) {
        const std::size_t cell_index = static_cast<std::size_t>(first_cell + offset);
        if (cell_covered_by_patch[cell_index] != 0U) {
          throw std::invalid_argument(std::string(context) + ": /state/patches cell ranges overlap");
        }
        cell_covered_by_patch[cell_index] = 1U;
        if (state.cells.patch_index[cell_index] != patch_index) {
          throw std::invalid_argument(
              std::string(context) + ": /state/cells/patch_index does not match /state/patches cell ranges");
        }
      }
    }
    for (std::size_t cell_index = 0; cell_index < cell_covered_by_patch.size(); ++cell_index) {
      if (cell_covered_by_patch[cell_index] == 0U) {
        throw std::invalid_argument(std::string(context) + ": /state/patches do not cover every gas cell");
      }
    }
  }
  if (!state.gas_cell_identity.isConsistent() ||
      !state.gas_cell_identity.coversDenseLocalRows(state.cells.size())) {
    throw std::invalid_argument(
        std::string(context) + ": /state/gas_cell_identity must cover dense local cell rows");
  }
  if (!state.gasCellIdentityMapMatchesSidecarLanes()) {
    throw std::invalid_argument(
        std::string(context) + ": /state/gas_cell_identity records must match gas-cell mirror lanes and patch ownership");
  }
  for (const core::GasCellIdentityRecord& record : state.gas_cell_identity.records()) {
    if (record.parent_particle_id.has_value() && *record.parent_particle_id == 0U) {
      throw std::invalid_argument(
          std::string(context) + ": /state/gas_cell_identity/parent_particle_id must be nonzero when has_parent_particle is true");
    }
  }
}

void validateOutputCadenceStateForRestart(
    const OutputCadencePersistentState& output_state,
    const core::IntegratorState& integrator_state,
    std::string_view context) {
  if (!output_state.output_enabled) {
    return;
  }
  if (output_state.snapshot_interval_steps == 0 &&
      output_state.snapshot_interval_time_code == 0.0) {
    throw std::invalid_argument(
        std::string(context) +
        ": /output_cadence requires a positive step or code-time interval when output is enabled");
  }
  if (!std::isfinite(output_state.snapshot_interval_time_code) ||
      output_state.snapshot_interval_time_code < 0.0) {
    throw std::invalid_argument(
        std::string(context) +
        ": /output_cadence/@snapshot_interval_time_code must be finite and >= 0");
  }
  if (output_state.last_completed_step_index != integrator_state.step_index) {
    throw std::invalid_argument(
        std::string(context) + ": /output_cadence/@last_completed_step_index must match /integrator/@step_index");
  }
  if (output_state.snapshot_interval_steps > 0 &&
      output_state.next_snapshot_step_index <= output_state.last_completed_step_index) {
    throw std::invalid_argument(
        std::string(context) + ": /output_cadence/@next_snapshot_step_index must be strictly after the restart step");
  }
  if (output_state.snapshot_interval_time_code > 0.0 &&
      (!std::isfinite(output_state.next_snapshot_time_code) ||
       output_state.next_snapshot_time_code <= integrator_state.current_time_code)) {
    throw std::invalid_argument(
        std::string(context) +
        ": /output_cadence/@next_snapshot_time_code must be finite and strictly after the restart time");
  }
}

void validateStochasticStateForRestart(
    const StochasticPersistentState& stochastic_state,
    const core::IntegratorState& integrator_state,
    std::string_view context) {
  std::vector<std::string> seen_modules;
  seen_modules.reserve(stochastic_state.modules.size());
  for (const StochasticModulePersistentState& module_state : stochastic_state.modules) {
    if (module_state.module_name.empty()) {
      throw std::invalid_argument(std::string(context) + ": /stochastic_state module_name must be non-empty");
    }
    if (module_state.module_name.find('/') != std::string::npos) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name +
          " module_name must not contain '/'");
    }
    if (std::find(seen_modules.begin(), seen_modules.end(), module_state.module_name) != seen_modules.end()) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state duplicate module_name '" + module_state.module_name + "'");
    }
    seen_modules.push_back(module_state.module_name);
    if (module_state.schema_version == 0) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name + "/@schema_version must be > 0");
    }
    if (module_state.rng_policy.empty()) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name + "/@rng_policy must be non-empty");
    }
    if (!module_state.deterministic_from_serialized_inputs) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name +
          " declares non-deterministic RNG state, but stateful RNG engine serialization is not implemented");
    }
    if (module_state.last_committed_step_index != integrator_state.step_index) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name +
          "/@last_committed_step_index must match /integrator/@step_index");
    }
  }
}

[[nodiscard]] std::vector<StochasticModulePersistentState> sortedStochasticModules(
    const StochasticPersistentState& stochastic_state) {
  std::vector<StochasticModulePersistentState> modules = stochastic_state.modules;
  std::sort(
      modules.begin(),
      modules.end(),
      [](const StochasticModulePersistentState& lhs, const StochasticModulePersistentState& rhs) {
        return lhs.module_name < rhs.module_name;
      });
  return modules;
}

[[nodiscard]] std::vector<core::GasCellIdentityRecord> gasIdentityRecordsSortedByLocalRow(
    const core::GasCellIdentityMap& identity_map) {
  std::vector<core::GasCellIdentityRecord> records(
      identity_map.records().begin(), identity_map.records().end());
  std::sort(
      records.begin(),
      records.end(),
      [](const core::GasCellIdentityRecord& lhs, const core::GasCellIdentityRecord& rhs) {
        return lhs.local_cell_row < rhs.local_cell_row;
      });
  return records;
}

#if COSMOSIM_ENABLE_HDF5

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t handle = -1) : m_handle(handle) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept : m_handle(other.m_handle) { other.m_handle = -1; }
  Hdf5Handle& operator=(Hdf5Handle&& other) noexcept {
    if (this != &other) {
      close();
      m_handle = other.m_handle;
      other.m_handle = -1;
    }
    return *this;
  }
  ~Hdf5Handle() { close(); }

  [[nodiscard]] hid_t get() const { return m_handle; }
  [[nodiscard]] bool valid() const { return m_handle >= 0; }

 private:
  void close() {
    if (m_handle < 0) {
      return;
    }
    const H5I_type_t type = H5Iget_type(m_handle);
    if (type == H5I_FILE) {
      H5Fclose(m_handle);
    } else if (type == H5I_GROUP) {
      H5Gclose(m_handle);
    } else if (type == H5I_DATASET) {
      H5Dclose(m_handle);
    } else if (type == H5I_DATASPACE) {
      H5Sclose(m_handle);
    } else if (type == H5I_ATTR) {
      H5Aclose(m_handle);
    } else if (type == H5I_DATATYPE) {
      H5Tclose(m_handle);
    }
    m_handle = -1;
  }

  hid_t m_handle = -1;
};

void writeScalarStringAttribute(hid_t location, std::string_view key, const std::string& value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle string_type(H5Tcopy(H5T_C_S1));
  if (!scalar_space.valid() || !string_type.valid()) {
    throw std::runtime_error("failed to create string attribute type");
  }
  if (H5Tset_size(string_type.get(), std::max<std::size_t>(std::size_t{1}, value.size() + 1)) < 0 || H5Tset_strpad(string_type.get(), H5T_STR_NULLTERM) < 0) {
    throw std::runtime_error("failed to configure string attribute type");
  }
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), string_type.get(), scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), string_type.get(), value.c_str()) < 0) {
    throw std::runtime_error("failed to write string attribute: " + std::string(key));
  }
}

[[nodiscard]] std::string readScalarStringAttribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    throw std::runtime_error("missing string attribute: " + std::string(key));
  }
  Hdf5Handle attr_type(H5Aget_type(attr.get()));
  const std::size_t length = H5Tget_size(attr_type.get());
  std::string value(length, '\0');
  if (H5Aread(attr.get(), attr_type.get(), value.data()) < 0) {
    throw std::runtime_error("failed reading string attribute: " + std::string(key));
  }
  while (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return value;
}

void writeScalarU32Attribute(hid_t location, std::string_view key, std::uint32_t value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), H5T_STD_U32LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error("failed writing u32 attribute: " + std::string(key));
  }
}

[[nodiscard]] std::uint32_t readScalarU32Attribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  std::uint32_t value = 0;
  if (!attr.valid() || H5Aread(attr.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error("failed reading u32 attribute: " + std::string(key));
  }
  return value;
}

void writeScalarU64Attribute(hid_t location, std::string_view key, std::uint64_t value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), H5T_STD_U64LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_UINT64, &value) < 0) {
    throw std::runtime_error("failed writing u64 attribute: " + std::string(key));
  }
}

[[nodiscard]] std::uint64_t readScalarU64Attribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  std::uint64_t value = 0;
  if (!attr.valid() || H5Aread(attr.get(), H5T_NATIVE_UINT64, &value) < 0) {
    throw std::runtime_error("failed reading u64 attribute: " + std::string(key));
  }
  return value;
}

void writeScalarF64Attribute(hid_t location, std::string_view key, double value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), H5T_IEEE_F64LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error("failed writing f64 attribute: " + std::string(key));
  }
}

[[nodiscard]] double readScalarF64Attribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  double value = 0.0;
  if (!attr.valid() || H5Aread(attr.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error("failed reading f64 attribute: " + std::string(key));
  }
  return value;
}

[[nodiscard]] bool hdf5LinkExists(hid_t root, std::string_view path) {
  return H5Lexists(root, std::string(path).c_str(), H5P_DEFAULT) > 0;
}

[[nodiscard]] bool hdf5AttributeExists(hid_t location, std::string_view key) {
  return H5Aexists(location, std::string(key).c_str()) > 0;
}

void requireHdf5Link(hid_t root, std::string_view path) {
  if (!hdf5LinkExists(root, path)) {
    throw std::runtime_error("restart schema validation missing required path: " + std::string(path));
  }
}

Hdf5Handle openRequiredGroup(hid_t root, std::string_view path) {
  requireHdf5Link(root, path);
  Hdf5Handle group(H5Gopen2(root, std::string(path).c_str(), H5P_DEFAULT));
  if (!group.valid()) {
    throw std::runtime_error("restart schema validation expected group at path: " + std::string(path));
  }
  return group;
}

void requireHdf5Attribute(hid_t location, std::string_view location_path, std::string_view key) {
  if (!hdf5AttributeExists(location, key)) {
    throw std::runtime_error(
        "restart schema validation missing required attribute: " + std::string(location_path) + "/@" +
        std::string(key));
  }
}

void requireHdf5Dataset1d(hid_t root, std::string_view path) {
  requireHdf5Link(root, path);
  Hdf5Handle dataset(H5Dopen2(root, std::string(path).c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("restart schema validation expected dataset at path: " + std::string(path));
  }
  Hdf5Handle space(H5Dget_space(dataset.get()));
  if (!space.valid() || H5Sget_simple_extent_ndims(space.get()) != 1) {
    throw std::runtime_error("restart schema validation expected 1D dataset at path: " + std::string(path));
  }
}

void validateFileKindAttribute(
    hid_t file,
    std::string_view expected_kind,
    std::string_view reader_label) {
  const auto& shared_names = sharedIoContractNames();
  if (!hdf5AttributeExists(file, shared_names.file_kind_attribute)) {
    std::string hint;
    if (hdf5LinkExists(file, "/Header")) {
      hint = "; file appears to be a science snapshot, not a restart checkpoint";
    }
    throw std::runtime_error(
        std::string(reader_label) + " missing required file-kind attribute /@" +
        std::string(shared_names.file_kind_attribute) + hint);
  }
  const std::string file_kind = readScalarStringAttribute(file, shared_names.file_kind_attribute);
  if (file_kind != expected_kind) {
    throw std::runtime_error(
        std::string(reader_label) + " rejected HDF5 file-kind '" + file_kind + "'; expected '" +
        std::string(expected_kind) + "'");
  }
}

void validateRestartCheckpointSchema(hid_t file, std::uint32_t schema_version) {
  const auto& shared_names = sharedIoContractNames();
  validateFileKindAttribute(file, shared_names.restart_checkpoint_file_kind, "restart reader");

  requireHdf5Attribute(file, "/", "restart_schema_name");
  requireHdf5Attribute(file, "/", "restart_schema_version");
  requireHdf5Attribute(file, "/", "normalized_config_hash_hex");
  requireHdf5Attribute(file, "/", "payload_integrity_hash_hex");
  requireHdf5Attribute(file, "/", "payload_integrity_hash");
  if (schema_version >= k_restart_schema_v23) {
    requireHdf5Attribute(file, "/", "payload_integrity_algorithm");
    requireHdf5Attribute(file, "/", "payload_integrity_sha256_hex");
  }
  requireHdf5Dataset1d(file, "/normalized_config_text");
  requireHdf5Dataset1d(file, "/provenance_record");

  requireHdf5Dataset1d(file, "/state/state_metadata");
  requireHdf5Dataset1d(file, "/state/particles/position_x_comoving");
  requireHdf5Dataset1d(file, "/state/particles/position_y_comoving");
  requireHdf5Dataset1d(file, "/state/particles/position_z_comoving");
  requireHdf5Dataset1d(file, "/state/particles/velocity_x_peculiar");
  requireHdf5Dataset1d(file, "/state/particles/velocity_y_peculiar");
  requireHdf5Dataset1d(file, "/state/particles/velocity_z_peculiar");
  requireHdf5Dataset1d(file, "/state/particles/mass_code");
  requireHdf5Dataset1d(file, "/state/particles/time_bin");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/particle_id");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/sfc_key");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/species_tag");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/particle_flags");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/owning_rank");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/last_drift_time_code");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/last_drift_scale_factor");
  requireHdf5Dataset1d(file, "/state/cells/center_x_comoving");
  requireHdf5Dataset1d(file, "/state/cells/center_y_comoving");
  requireHdf5Dataset1d(file, "/state/cells/center_z_comoving");
  requireHdf5Dataset1d(file, "/state/cells/mass_code");
  requireHdf5Dataset1d(file, "/state/cells/time_bin");
  requireHdf5Dataset1d(file, "/state/cells/patch_index");
  requireHdf5Dataset1d(file, "/state/gas_cells/gas_cell_id");
  requireHdf5Dataset1d(file, "/state/gas_cells/parent_particle_id");
  requireHdf5Dataset1d(file, "/state/gas_cells/velocity_x_peculiar");
  requireHdf5Dataset1d(file, "/state/gas_cells/velocity_y_peculiar");
  requireHdf5Dataset1d(file, "/state/gas_cells/velocity_z_peculiar");
  requireHdf5Dataset1d(file, "/state/gas_cells/density_code");
  requireHdf5Dataset1d(file, "/state/gas_cells/pressure_code");
  requireHdf5Dataset1d(file, "/state/gas_cells/internal_energy_code");
  if (schema_version >= k_restart_schema_v22) {
    requireHdf5Dataset1d(file, "/state/gas_cells/metal_mass_code");
  }
  requireHdf5Dataset1d(file, "/state/gas_cells/temperature_code");
  requireHdf5Dataset1d(file, "/state/gas_cells/sound_speed_code");
  if (schema_version >= k_restart_schema_v18) {
    Hdf5Handle identity_group = openRequiredGroup(file, "/state/gas_cell_identity");
    requireHdf5Attribute(identity_group.get(), "/state/gas_cell_identity", "local_row_reconstruction_policy");
    requireHdf5Attribute(identity_group.get(), "/state/gas_cell_identity", "identity_generation_at_write");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/gas_cell_id");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/has_parent_particle");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/parent_particle_id");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/owning_patch_id");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/local_cell_row");
  }
  requireHdf5Dataset1d(file, "/state/patches/patch_id");
  requireHdf5Dataset1d(file, "/state/patches/level");
  requireHdf5Dataset1d(file, "/state/patches/first_cell");
  requireHdf5Dataset1d(file, "/state/patches/cell_count");
  if (schema_version >= k_restart_schema_v18) {
    requireHdf5Dataset1d(file, "/state/patches/parent_patch_id");
    requireHdf5Dataset1d(file, "/state/patches/morton_key");
    requireHdf5Dataset1d(file, "/state/patches/origin_x_comoving");
    requireHdf5Dataset1d(file, "/state/patches/origin_y_comoving");
    requireHdf5Dataset1d(file, "/state/patches/origin_z_comoving");
    requireHdf5Dataset1d(file, "/state/patches/extent_x_comoving");
    requireHdf5Dataset1d(file, "/state/patches/extent_y_comoving");
    requireHdf5Dataset1d(file, "/state/patches/extent_z_comoving");
    requireHdf5Dataset1d(file, "/state/patches/cell_dim_x");
    requireHdf5Dataset1d(file, "/state/patches/cell_dim_y");
    requireHdf5Dataset1d(file, "/state/patches/cell_dim_z");
  }
  if (schema_version >= k_restart_schema_v18) {
    Hdf5Handle pending_group = openRequiredGroup(file, "/state/amr_pending_flux_registers");
    requireHdf5Attribute(pending_group.get(), "/state/amr_pending_flux_registers", "schema_version");
    for (std::string_view dataset : {"register_key", "coarse_patch_id", "coarse_gas_cell_id",
                                     "coarse_cell_index", "level", "axis", "orientation",
                                     "expected_area_comov", "coarse_area_accumulated_comov",
                                     "fine_area_accumulated_comov", "interval_start_code",
                                     "interval_end_code", "coarse_dt_code", "expected_fine_substeps",
                                     "completed_fine_substeps", "fine_substep_coverage_mask",
                                     "coarse_face_count", "fine_face_count",
                                     "gas_cell_identity_generation", "patch_geometry_generation",
                                     "coarse_mass_flux_integral_code",
                                     "coarse_momentum_x_flux_integral_code",
                                     "coarse_momentum_y_flux_integral_code",
                                     "coarse_momentum_z_flux_integral_code",
                                     "coarse_total_energy_flux_integral_code",
                                     "fine_mass_flux_integral_code", "fine_momentum_x_flux_integral_code",
                                     "fine_momentum_y_flux_integral_code",
                                     "fine_momentum_z_flux_integral_code",
                                     "fine_total_energy_flux_integral_code"}) {
      requireHdf5Dataset1d(file, std::string("/state/amr_pending_flux_registers/") + std::string(dataset));
    }
    if (schema_version >= k_restart_schema_v22) {
      requireHdf5Dataset1d(file, "/state/amr_pending_flux_registers/coarse_metal_mass_flux_integral_code");
      requireHdf5Dataset1d(file, "/state/amr_pending_flux_registers/fine_metal_mass_flux_integral_code");
    }
  }
  if (schema_version >= k_restart_schema_v18) {
    Hdf5Handle temporal_group = openRequiredGroup(file, "/state/amr_temporal_boundary_history");
    requireHdf5Attribute(temporal_group.get(), "/state/amr_temporal_boundary_history", "schema_version");
    for (std::string_view dataset : {"patch_id", "patch_level", "patch_geometry_fingerprint",
                                     "gas_cell_identity_generation", "interval_start_code",
                                     "interval_end_code", "end_state_valid", "cell_offset", "cell_count",
                                     "gas_cell_id", "patch_local_cell", "start_mass_density_comoving",
                                     "start_momentum_density_x_comoving", "start_momentum_density_y_comoving",
                                     "start_momentum_density_z_comoving", "start_total_energy_density_comoving",
                                     "end_mass_density_comoving", "end_momentum_density_x_comoving",
                                     "end_momentum_density_y_comoving", "end_momentum_density_z_comoving",
                                     "end_total_energy_density_comoving"}) {
      requireHdf5Dataset1d(file, std::string("/state/amr_temporal_boundary_history/") + std::string(dataset));
    }
    if (schema_version >= k_restart_schema_v22) {
      requireHdf5Dataset1d(file, "/state/amr_temporal_boundary_history/start_metal_mass_density_comoving");
      requireHdf5Dataset1d(file, "/state/amr_temporal_boundary_history/end_metal_mass_density_comoving");
    }
  }
  requireHdf5Dataset1d(file, "/state/species_count_by_species");
  requireHdf5Link(file, "/state/module_sidecars");
  requireHdf5Dataset1d(file, "/state/module_sidecar_names");

  Hdf5Handle integrator_group = openRequiredGroup(file, "/integrator");
  for (std::string_view attr : {"current_time_code", "current_scale_factor", "current_redshift",
                                "current_hubble_rate_code", "time_si_per_code", "dt_time_code",
                                "last_drift_factor_code", "last_first_kick_factor_code",
                                "last_second_kick_factor_code", "last_first_hubble_drag_factor",
                                "last_second_hubble_drag_factor", "step_index", "scheme",
                                "current_boundary_kind", "last_completed_boundary_kind", "inside_kdk_step",
                                "last_completed_restart_safe", "time_bins_hierarchical", "time_bins_active_bin",
                                "time_bins_max_bin", "pm_long_range_field_valid", "pm_cadence_steps",
                                "pm_gravity_kick_opportunity", "pm_last_refresh_opportunity",
                                "pm_field_version", "pm_last_refresh_step_index",
                                "pm_last_refresh_scale_factor", "pm_refresh_commit_pending",
                                "pm_pending_refresh_opportunity", "pm_pending_refresh_field_version"}) {
    requireHdf5Attribute(integrator_group.get(), "/integrator", attr);
  }
  if (schema_version >= k_restart_schema_v20) {
    requireHdf5Attribute(integrator_group.get(), "/integrator", "pm_refresh_enabled");
  }

  if (schema_version >= k_restart_schema_v20) {
    Hdf5Handle force_cache_group = openRequiredGroup(file, "/gravity_force_cache");
    requireHdf5Attribute(force_cache_group.get(), "/gravity_force_cache", "valid");
    for (std::string_view dataset : {"particle_id", "gas_cell_id",
                                     "particle_accel_x_comoving", "particle_accel_y_comoving",
                                     "particle_accel_z_comoving", "cell_accel_x_comoving",
                                     "cell_accel_y_comoving", "cell_accel_z_comoving"}) {
      requireHdf5Dataset1d(file, std::string("/gravity_force_cache/") + std::string(dataset));
    }
  }
  if (schema_version >= k_restart_schema_v24) {
    Hdf5Handle endpoint_group = openRequiredGroup(file, std::string(k_hierarchical_endpoint_group));
    requireHdf5Attribute(endpoint_group.get(), std::string(k_hierarchical_endpoint_group), "valid");
    requireHdf5Attribute(endpoint_group.get(), std::string(k_hierarchical_endpoint_group), "pm_field_version");
    for (std::string_view dataset : {"tree_accel_x_comoving", "tree_accel_y_comoving", "tree_accel_z_comoving",
                                     "pm_accel_x_comoving", "pm_accel_y_comoving", "pm_accel_z_comoving",
                                     "nearest_pair_r2_comoving"}) {
      requireHdf5Dataset1d(file, std::string(k_hierarchical_endpoint_group) + "/" + std::string(dataset));
    }
  }

  Hdf5Handle scheduler_group = openRequiredGroup(file, "/scheduler");
  requireHdf5Attribute(scheduler_group.get(), "/scheduler", "current_tick");
  requireHdf5Attribute(scheduler_group.get(), "/scheduler", "max_bin");
  requireHdf5Dataset1d(file, "/scheduler/bin_index");
  requireHdf5Dataset1d(file, "/scheduler/next_activation_tick");
  requireHdf5Dataset1d(file, "/scheduler/active_flag");
  requireHdf5Dataset1d(file, "/scheduler/pending_bin_index");
  if (schema_version >= k_restart_schema_v19) {
    Hdf5Handle gas_scheduler_group = openRequiredGroup(file, "/gas_cell_scheduler");
    requireHdf5Attribute(gas_scheduler_group.get(), "/gas_cell_scheduler", "identity_key");
    requireHdf5Attribute(gas_scheduler_group.get(), "/gas_cell_scheduler", "current_tick");
    requireHdf5Attribute(gas_scheduler_group.get(), "/gas_cell_scheduler", "max_bin");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/gas_cell_id");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/bin_index");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/next_activation_tick");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/active_flag");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/pending_bin_index");
  }

  Hdf5Handle output_group = openRequiredGroup(file, "/output_cadence");
  for (std::string_view attr : {"output_enabled", "write_restarts", "snapshot_due", "checkpoint_due",
                                "last_completed_step_index", "snapshot_interval_steps",
                                "next_snapshot_step_index", "snapshot_stem", "restart_stem"}) {
    requireHdf5Attribute(output_group.get(), "/output_cadence", attr);
  }
  if (schema_version >= k_restart_schema_v21) {
    requireHdf5Attribute(output_group.get(), "/output_cadence", "snapshot_interval_time_code");
    requireHdf5Attribute(output_group.get(), "/output_cadence", "next_snapshot_time_code");
  }

  Hdf5Handle stochastic_group = openRequiredGroup(file, "/stochastic_state");
  requireHdf5Attribute(stochastic_group.get(), "/stochastic_state", "module_count");
  requireHdf5Dataset1d(file, "/stochastic_state/module_names");

  Hdf5Handle diagnostics_group = openRequiredGroup(file, "/restart_diagnostics");
  for (std::string_view attr : {"restart_schema_name", "restart_schema_version", "current_boundary_kind",
                                "last_completed_boundary_kind", "restart_safe", "step_index",
                                "scheduler_current_tick", "scheduler_max_bin", "scheduler_element_count",
                                "scheduler_active_count", "scheduler_pending_transition_count",
                                "pm_cadence_steps", "pm_gravity_kick_opportunity", "pm_field_version",
                                "pm_last_refresh_opportunity", "pm_last_refresh_step_index",
                                "pm_refresh_commit_pending", "pm_long_range_field_valid", "output_enabled",
                                "output_snapshot_due", "output_checkpoint_due",
                                "output_last_completed_step_index", "output_next_snapshot_step_index",
                                "stochastic_module_count"}) {
    requireHdf5Attribute(diagnostics_group.get(), "/restart_diagnostics", attr);
  }
  if (schema_version >= k_restart_schema_v21) {
    requireHdf5Attribute(
        diagnostics_group.get(), "/restart_diagnostics", "output_next_snapshot_time_code");
  }

  if (schema_version >= k_restart_schema_v19) {
    for (std::string_view attr : {"gas_cell_scheduler_current_tick", "gas_cell_scheduler_max_bin",
                                  "gas_cell_scheduler_element_count", "gas_cell_scheduler_active_count",
                                  "gas_cell_scheduler_pending_transition_count"}) {
      requireHdf5Attribute(diagnostics_group.get(), "/restart_diagnostics", attr);
    }
  }

  requireHdf5Dataset1d(file, "/distributed_gravity/state");
}

template <typename T>
void writeDataset1d(
    hid_t group,
    std::string_view name,
    hid_t file_type,
    hid_t memory_type,
    std::span<const T> values) {
  hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group, std::string(name).c_str(), file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed creating dataset: " + std::string(name));
  }
  if (!values.empty() && H5Dwrite(dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    throw std::runtime_error("failed writing dataset: " + std::string(name));
  }
}

template <typename T, typename Allocator>
void writeDataset1d(
    hid_t group,
    std::string_view name,
    hid_t file_type,
    hid_t memory_type,
    const std::vector<T, Allocator>& values) {
  writeDataset1d<T>(group, name, file_type, memory_type, std::span<const T>(values.data(), values.size()));
}


template <typename T>
void writeUniformDataset1d(
    hid_t group, std::string_view name, hid_t file_type, hid_t memory_type,
    std::size_t count, T value) {
  hsize_t dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group, std::string(name).c_str(), file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!dataset.valid()) throw std::runtime_error("failed creating dataset: " + std::string(name));
  constexpr std::size_t k_chunk_elements = 64U * 1024U;
  std::vector<T> chunk(std::min(count, k_chunk_elements), value);
  for (std::size_t offset = 0; offset < count; offset += chunk.size()) {
    const std::size_t n = std::min(chunk.size(), count - offset);
    hsize_t start[1] = {static_cast<hsize_t>(offset)};
    hsize_t extent[1] = {static_cast<hsize_t>(n)};
    Hdf5Handle file_space(H5Dget_space(dataset.get()));
    if (!file_space.valid() || H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, start, nullptr, extent, nullptr) < 0) {
      throw std::runtime_error("failed selecting uniform dataset hyperslab: " + std::string(name));
    }
    Hdf5Handle memory_space(H5Screate_simple(1, extent, nullptr));
    if (!memory_space.valid() || H5Dwrite(dataset.get(), memory_type, memory_space.get(), file_space.get(), H5P_DEFAULT, chunk.data()) < 0) {
      throw std::runtime_error("failed writing uniform dataset: " + std::string(name));
    }
  }
}

// Read directly into the destination allocator. In particular, an aligned
// canonical lane must not first materialize a second full std::vector.
template <typename T, typename Allocator>
void readDataset1dInto(hid_t group, std::string_view name, hid_t memory_type,
                       std::vector<T, Allocator>& values) {
  const std::string dataset_name(name);
  const htri_t exists = H5Lexists(group, dataset_name.c_str(), H5P_DEFAULT);
  if (exists < 0) {
    throw std::runtime_error("failed checking dataset presence: " + dataset_name);
  }
  if (exists == 0) {
    throw std::runtime_error("missing dataset: " + dataset_name);
  }
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("missing dataset: " + dataset_name);
  }
  Hdf5Handle space(H5Dget_space(dataset.get()));
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_ndims(space.get()) != 1 ||
      H5Sget_simple_extent_dims(space.get(), dims, nullptr) != 1) {
    throw std::runtime_error("unexpected rank for dataset: " + dataset_name);
  }
  const std::size_t count = core::checkedIntegralNarrow<std::size_t>(
      dims[0], "restart dataset element count");
  (void)core::checkedSizeMultiply(count, sizeof(T), "restart dataset byte count");
  if (count > values.max_size()) {
    throw std::length_error("restart dataset exceeds destination capacity: " + dataset_name);
  }
  values.resize(count);
  if (!values.empty() && H5Dread(dataset.get(), memory_type, H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT, values.data()) < 0) {
    throw std::runtime_error("failed reading dataset: " + dataset_name);
  }
}

template <typename T>
[[nodiscard]] std::vector<T> readDataset1d(hid_t group, std::string_view name, hid_t memory_type) {
  std::vector<T> values;
  readDataset1dInto(group, name, memory_type, values);
  return values;
}

template <typename T>
[[nodiscard]] core::AlignedVector<T> readDataset1dAligned(
    hid_t group, std::string_view name, hid_t memory_type) {
  core::AlignedVector<T> values;
  readDataset1dInto(group, name, memory_type, values);
  return values;
}

void writeStringDataset(hid_t group, std::string_view name, const std::string& value) {
  // HDF5's byte dataset accepts the canonical string directly; no second
  // string-sized byte vector is necessary. The on-disk schema is unchanged.
  writeDataset1d<std::uint8_t>(group, name, H5T_STD_U8LE, H5T_NATIVE_UINT8,
      std::span<const std::uint8_t>(
          reinterpret_cast<const std::uint8_t*>(value.data()), value.size()));
}

[[nodiscard]] std::string readStringDataset(hid_t group, std::string_view name) {
  const std::string dataset_name(name);
  const htri_t exists = H5Lexists(group, dataset_name.c_str(), H5P_DEFAULT);
  if (exists < 0) throw std::runtime_error("failed checking dataset presence: " + dataset_name);
  if (exists == 0) throw std::runtime_error("missing dataset: " + dataset_name);
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) throw std::runtime_error("missing dataset: " + dataset_name);
  Hdf5Handle space(H5Dget_space(dataset.get()));
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_ndims(space.get()) != 1 ||
      H5Sget_simple_extent_dims(space.get(), dims, nullptr) != 1) {
    throw std::runtime_error("unexpected rank for dataset: " + dataset_name);
  }
  const std::size_t count = core::checkedIntegralNarrow<std::size_t>(
      dims[0], "restart string byte count");
  std::string value;
  if (count > value.max_size()) throw std::length_error("restart string exceeds maximum size");
  value.resize(count);
  if (!value.empty() && H5Dread(dataset.get(), H5T_NATIVE_UINT8, H5S_ALL,
                                H5S_ALL, H5P_DEFAULT, value.data()) < 0) {
    throw std::runtime_error("failed reading dataset: " + dataset_name);
  }
  return value;
}

[[nodiscard]] hid_t openOrCreateGroup(hid_t parent, const std::string& name) {
  if (H5Lexists(parent, name.c_str(), H5P_DEFAULT) > 0) {
    return H5Gopen2(parent, name.c_str(), H5P_DEFAULT);
  }
  return H5Gcreate2(parent, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void writeStarSidecarGroup(hid_t state_group, const core::StarParticleSidecar& stars) {
  Hdf5Handle star_group(openOrCreateGroup(state_group, "star_particles"));
  writeDataset1d(star_group.get(), "particle_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, stars.particle_index);
  writeDataset1d(
      star_group.get(), "formation_scale_factor", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.formation_scale_factor);
  writeDataset1d(star_group.get(), "birth_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.birth_mass_code);
  writeDataset1d(
      star_group.get(),
      "metallicity_mass_fraction",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.metallicity_mass_fraction);
  writeDataset1d(star_group.get(), "birth_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, stars.birth_key);
  writeDataset1d(
      star_group.get(), "parent_gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, stars.parent_gas_cell_id);
  writeDataset1d(star_group.get(), "birth_tick", H5T_STD_U64LE, H5T_NATIVE_UINT64, stars.birth_tick);
  writeDataset1d(star_group.get(), "birth_ordinal", H5T_STD_U32LE, H5T_NATIVE_UINT32, stars.birth_ordinal);
  writeDataset1d(
      star_group.get(),
      "stellar_age_years_last",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.stellar_age_years_last);
  writeDataset1d(
      star_group.get(),
      "stellar_returned_mass_cumulative_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.stellar_returned_mass_cumulative_code);
  writeDataset1d(
      star_group.get(),
      "stellar_returned_metals_cumulative_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.stellar_returned_metals_cumulative_code);
  writeDataset1d(star_group.get(), "stellar_newly_synthesized_metals_cumulative_code",
      H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.stellar_newly_synthesized_metals_cumulative_code);
  writeDataset1d(star_group.get(), "enrichment_carry_mass_code",
      H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.enrichment_carry_mass_code);
  writeDataset1d(star_group.get(), "enrichment_carry_metals_code",
      H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.enrichment_carry_metals_code);
  writeDataset1d(star_group.get(), "enrichment_carry_feedback_energy_erg",
      H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.enrichment_carry_feedback_energy_erg);
  writeDataset1d(star_group.get(), "enrichment_carry_momentum_code",
      H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.enrichment_carry_momentum_code);
  writeDataset1d(star_group.get(), "stellar_deposited_mass_cumulative_code",
      H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.stellar_deposited_mass_cumulative_code);
  writeDataset1d(star_group.get(), "stellar_deposited_metals_cumulative_code",
      H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.stellar_deposited_metals_cumulative_code);
  writeDataset1d(star_group.get(), "stellar_deposited_feedback_energy_cumulative_erg",
      H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.stellar_deposited_feedback_energy_cumulative_erg);
  writeDataset1d(
      star_group.get(),
      "stellar_feedback_energy_cumulative_erg",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.stellar_feedback_energy_cumulative_erg);

  for (std::size_t channel = 0; channel < stars.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    const std::string suffix = std::to_string(channel);
    writeDataset1d(
        star_group.get(),
        "stellar_returned_mass_channel_cumulative_code_" + suffix,
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        stars.stellar_returned_mass_channel_cumulative_code[channel]);
    writeDataset1d(
        star_group.get(),
        "stellar_returned_metals_channel_cumulative_code_" + suffix,
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        stars.stellar_returned_metals_channel_cumulative_code[channel]);
    writeDataset1d(
        star_group.get(),
        "stellar_feedback_energy_channel_cumulative_erg_" + suffix,
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        stars.stellar_feedback_energy_channel_cumulative_erg[channel]);
  }
}

void readStarSidecarGroup(hid_t state_group, core::StarParticleSidecar& stars) {
  Hdf5Handle star_group(H5Gopen2(state_group, "star_particles", H5P_DEFAULT));
  stars.particle_index = readDataset1dAligned<std::uint32_t>(star_group.get(), "particle_index", H5T_NATIVE_UINT32);
  stars.formation_scale_factor =
      readDataset1dAligned<double>(star_group.get(), "formation_scale_factor", H5T_NATIVE_DOUBLE);
  stars.birth_mass_code = readDataset1dAligned<double>(star_group.get(), "birth_mass_code", H5T_NATIVE_DOUBLE);
  stars.metallicity_mass_fraction =
      readDataset1dAligned<double>(star_group.get(), "metallicity_mass_fraction", H5T_NATIVE_DOUBLE);
  if (H5Lexists(star_group.get(), "birth_key", H5P_DEFAULT) > 0) {
    stars.birth_key = readDataset1dAligned<std::uint64_t>(star_group.get(), "birth_key", H5T_NATIVE_UINT64);
    stars.parent_gas_cell_id =
        readDataset1dAligned<std::uint64_t>(star_group.get(), "parent_gas_cell_id", H5T_NATIVE_UINT64);
    stars.birth_tick = readDataset1dAligned<std::uint64_t>(star_group.get(), "birth_tick", H5T_NATIVE_UINT64);
    stars.birth_ordinal =
        readDataset1dAligned<std::uint32_t>(star_group.get(), "birth_ordinal", H5T_NATIVE_UINT32);
  } else {
    stars.birth_key.assign(stars.particle_index.size(), 0U);
    stars.parent_gas_cell_id.assign(stars.particle_index.size(), 0U);
    stars.birth_tick.assign(stars.particle_index.size(), 0U);
    stars.birth_ordinal.assign(stars.particle_index.size(), 0U);
  }
  stars.stellar_age_years_last =
      readDataset1dAligned<double>(star_group.get(), "stellar_age_years_last", H5T_NATIVE_DOUBLE);
  stars.stellar_returned_mass_cumulative_code = readDataset1dAligned<double>(
      star_group.get(), "stellar_returned_mass_cumulative_code", H5T_NATIVE_DOUBLE);
  stars.stellar_returned_metals_cumulative_code = readDataset1dAligned<double>(
      star_group.get(), "stellar_returned_metals_cumulative_code", H5T_NATIVE_DOUBLE);
  const std::size_t star_count = stars.particle_index.size();
  auto read_optional_double_lane = [&](const char* name) {
    if (H5Lexists(star_group.get(), name, H5P_DEFAULT) > 0) {
      return readDataset1dAligned<double>(star_group.get(), name, H5T_NATIVE_DOUBLE);
    }
    return core::AlignedVector<double>(star_count, 0.0);
  };
  stars.stellar_newly_synthesized_metals_cumulative_code =
      read_optional_double_lane("stellar_newly_synthesized_metals_cumulative_code");
  stars.enrichment_carry_mass_code = read_optional_double_lane("enrichment_carry_mass_code");
  stars.enrichment_carry_metals_code = read_optional_double_lane("enrichment_carry_metals_code");
  stars.enrichment_carry_feedback_energy_erg =
      read_optional_double_lane("enrichment_carry_feedback_energy_erg");
  stars.enrichment_carry_momentum_code =
      read_optional_double_lane("enrichment_carry_momentum_code");
  stars.stellar_deposited_mass_cumulative_code =
      read_optional_double_lane("stellar_deposited_mass_cumulative_code");
  stars.stellar_deposited_metals_cumulative_code =
      read_optional_double_lane("stellar_deposited_metals_cumulative_code");
  stars.stellar_deposited_feedback_energy_cumulative_erg =
      read_optional_double_lane("stellar_deposited_feedback_energy_cumulative_erg");
  stars.stellar_feedback_energy_cumulative_erg = readDataset1dAligned<double>(
      star_group.get(), "stellar_feedback_energy_cumulative_erg", H5T_NATIVE_DOUBLE);

  for (std::size_t channel = 0; channel < stars.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    const std::string suffix = std::to_string(channel);
    stars.stellar_returned_mass_channel_cumulative_code[channel] = readDataset1dAligned<double>(
        star_group.get(),
        "stellar_returned_mass_channel_cumulative_code_" + suffix,
        H5T_NATIVE_DOUBLE);
    stars.stellar_returned_metals_channel_cumulative_code[channel] = readDataset1dAligned<double>(
        star_group.get(),
        "stellar_returned_metals_channel_cumulative_code_" + suffix,
        H5T_NATIVE_DOUBLE);
    stars.stellar_feedback_energy_channel_cumulative_erg[channel] = readDataset1dAligned<double>(
        star_group.get(),
        "stellar_feedback_energy_channel_cumulative_erg_" + suffix,
        H5T_NATIVE_DOUBLE);
  }
}

void writeGasCellIdentityGroup(hid_t state_group, const core::SimulationState& state) {
  Hdf5Handle identity_group(openOrCreateGroup(state_group, "gas_cell_identity"));
  writeScalarStringAttribute(identity_group.get(), "local_row_reconstruction_policy", std::string(k_gas_identity_row_policy));
  writeScalarU64Attribute(identity_group.get(), "identity_generation_at_write", state.gasCellIdentityGeneration());

  const auto records = gasIdentityRecordsSortedByLocalRow(state.gas_cell_identity);
  std::vector<std::uint64_t> gas_cell_id;
  std::vector<std::uint8_t> has_parent_particle;
  std::vector<std::uint64_t> parent_particle_id;
  std::vector<std::uint64_t> owning_patch_id;
  std::vector<std::uint32_t> local_cell_row;
  gas_cell_id.reserve(records.size());
  has_parent_particle.reserve(records.size());
  parent_particle_id.reserve(records.size());
  owning_patch_id.reserve(records.size());
  local_cell_row.reserve(records.size());
  for (const core::GasCellIdentityRecord& record : records) {
    gas_cell_id.push_back(record.gas_cell_id);
    has_parent_particle.push_back(record.parent_particle_id.has_value() ? 1U : 0U);
    parent_particle_id.push_back(record.parent_particle_id.value_or(0U));
    owning_patch_id.push_back(record.owning_patch_id);
    local_cell_row.push_back(record.local_cell_row);
  }

  writeDataset1d(identity_group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_id);
  writeDataset1d(identity_group.get(), "has_parent_particle", H5T_STD_U8LE, H5T_NATIVE_UINT8, has_parent_particle);
  writeDataset1d(identity_group.get(), "parent_particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, parent_particle_id);
  writeDataset1d(identity_group.get(), "owning_patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, owning_patch_id);
  writeDataset1d(identity_group.get(), "local_cell_row", H5T_STD_U32LE, H5T_NATIVE_UINT32, local_cell_row);
}

void materializeLegacyGasCellIdentityMapFromMirrors(core::SimulationState& state) {
  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint64_t gas_cell_id = state.gas_cells.gas_cell_id[cell_index];
    const std::uint64_t parent_particle_id = state.gas_cells.parent_particle_id[cell_index];
    if (gas_cell_id == 0U || parent_particle_id == 0U || gas_cell_id != parent_particle_id) {
      throw std::runtime_error(
          "legacy v14 restart gas-cell identity requires nonzero gas_cell_id == parent_particle_id");
    }
    std::uint64_t owning_patch_id = 0U;
    if (!state.patches.patch_id.empty()) {
      const std::uint32_t patch_index = state.cells.patch_index[cell_index];
      if (patch_index >= state.patches.size()) {
        throw std::runtime_error("legacy v14 restart gas-cell identity has invalid patch_index");
      }
      owning_patch_id = state.patches.patch_id[patch_index];
    }
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = parent_particle_id,
        .owning_patch_id = owning_patch_id,
        .local_cell_row = static_cast<std::uint32_t>(cell_index),
    });
  }
  state.replaceGasCellIdentityRecords(std::move(records));
}

void readGasCellIdentityGroup(hid_t state_group, core::SimulationState& state, std::uint32_t schema_version) {
  if (schema_version == k_restart_schema_v14 && H5Lexists(state_group, "gas_cell_identity", H5P_DEFAULT) <= 0) {
    materializeLegacyGasCellIdentityMapFromMirrors(state);
    return;
  }

  Hdf5Handle identity_group(H5Gopen2(state_group, "gas_cell_identity", H5P_DEFAULT));
  if (!identity_group.valid()) {
    throw std::runtime_error("restart is missing /state/gas_cell_identity");
  }
  const std::string row_policy =
      readScalarStringAttribute(identity_group.get(), "local_row_reconstruction_policy");
  if (row_policy != k_gas_identity_row_policy) {
    throw std::runtime_error("unsupported /state/gas_cell_identity local_row_reconstruction_policy");
  }
  const std::uint64_t identity_generation_at_write =
      readScalarU64Attribute(identity_group.get(), "identity_generation_at_write");

  const auto gas_cell_id =
      readDataset1d<std::uint64_t>(identity_group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
  const auto has_parent_particle =
      readDataset1d<std::uint8_t>(identity_group.get(), "has_parent_particle", H5T_NATIVE_UINT8);
  const auto parent_particle_id =
      readDataset1d<std::uint64_t>(identity_group.get(), "parent_particle_id", H5T_NATIVE_UINT64);
  const auto owning_patch_id =
      readDataset1d<std::uint64_t>(identity_group.get(), "owning_patch_id", H5T_NATIVE_UINT64);
  const auto local_cell_row =
      readDataset1d<std::uint32_t>(identity_group.get(), "local_cell_row", H5T_NATIVE_UINT32);
  if (gas_cell_id.size() != has_parent_particle.size() ||
      gas_cell_id.size() != parent_particle_id.size() ||
      gas_cell_id.size() != owning_patch_id.size() ||
      gas_cell_id.size() != local_cell_row.size()) {
    throw std::runtime_error("/state/gas_cell_identity datasets must have matching lengths");
  }

  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(gas_cell_id.size());
  for (std::size_t i = 0; i < gas_cell_id.size(); ++i) {
    if (has_parent_particle[i] > 1U) {
      throw std::runtime_error("/state/gas_cell_identity/has_parent_particle must contain only 0 or 1");
    }
    if (has_parent_particle[i] == 0U && parent_particle_id[i] != 0U) {
      throw std::runtime_error(
          "/state/gas_cell_identity/parent_particle_id must be 0 when has_parent_particle is false");
    }
    if (has_parent_particle[i] != 0U && parent_particle_id[i] == 0U) {
      throw std::runtime_error(
          "/state/gas_cell_identity/parent_particle_id must be nonzero when has_parent_particle is true");
    }
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = gas_cell_id[i],
        .parent_particle_id =
            has_parent_particle[i] != 0U ? std::optional<std::uint64_t>(parent_particle_id[i]) : std::nullopt,
        .owning_patch_id = owning_patch_id[i],
        .local_cell_row = local_cell_row[i],
    });
  }
  // Validate persisted compatibility mirrors before normalizing them from the
  // canonical map.  Otherwise a corrupt /state/cells/patch_index lane could be
  // silently overwritten during restore and evade restart integrity checks.
  if (records.size() != state.cells.size()) {
    throw std::runtime_error("/state/gas_cell_identity record count must match /state/cells extent");
  }
  for (const core::GasCellIdentityRecord& record : records) {
    if (record.local_cell_row >= state.cells.size()) {
      throw std::runtime_error("/state/gas_cell_identity/local_cell_row exceeds /state/cells extent");
    }
    const std::uint32_t row = record.local_cell_row;
    if (state.gas_cells.gas_cell_id[row] != record.gas_cell_id ||
        state.gas_cells.parent_particle_id[row] != record.parent_particle_id.value_or(0U)) {
      throw std::runtime_error("/state/gas_cell_identity records disagree with gas-cell compatibility mirrors");
    }
    if (record.owning_patch_id != 0U) {
      if (state.cells.patch_index[row] >= state.patches.size() ||
          state.patches.patch_id[state.cells.patch_index[row]] != record.owning_patch_id) {
        throw std::runtime_error("/state/cells/patch_index disagrees with /state/gas_cell_identity owning_patch_id");
      }
    }
  }
  state.restoreGasCellIdentityRecords(std::move(records), identity_generation_at_write);
}

void writePendingFluxRegisterGroup(hid_t state_group, const core::PendingFluxRegisterStore& store) {
  Hdf5Handle group(openOrCreateGroup(state_group, "amr_pending_flux_registers"));
  writeScalarU32Attribute(group.get(), "schema_version", 2U);
  const auto records = store.records();
  std::vector<std::uint64_t> register_key;
  std::vector<std::uint64_t> coarse_patch_id;
  std::vector<std::uint64_t> coarse_gas_cell_id;
  std::vector<std::uint64_t> coarse_cell_index;
  std::vector<std::uint8_t> level;
  std::vector<std::uint8_t> axis;
  std::vector<std::uint8_t> orientation;
  std::vector<double> expected_area_comov;
  std::vector<double> coarse_area_accumulated_comov;
  std::vector<double> fine_area_accumulated_comov;
  std::vector<double> interval_start_code;
  std::vector<double> interval_end_code;
  std::vector<double> coarse_dt_code;
  std::vector<std::uint32_t> expected_fine_substeps;
  std::vector<std::uint32_t> completed_fine_substeps;
  std::vector<std::uint64_t> fine_substep_coverage_mask;
  std::vector<std::uint32_t> coarse_face_count;
  std::vector<std::uint32_t> fine_face_count;
  std::vector<std::uint64_t> gas_cell_identity_generation;
  std::vector<std::uint64_t> patch_geometry_generation;
  std::vector<double> coarse_mass_flux_integral_code;
  std::vector<double> coarse_momentum_x_flux_integral_code;
  std::vector<double> coarse_momentum_y_flux_integral_code;
  std::vector<double> coarse_momentum_z_flux_integral_code;
  std::vector<double> coarse_total_energy_flux_integral_code;
  std::vector<double> coarse_metal_mass_flux_integral_code;
  std::vector<double> fine_mass_flux_integral_code;
  std::vector<double> fine_momentum_x_flux_integral_code;
  std::vector<double> fine_momentum_y_flux_integral_code;
  std::vector<double> fine_momentum_z_flux_integral_code;
  std::vector<double> fine_total_energy_flux_integral_code;
  std::vector<double> fine_metal_mass_flux_integral_code;
  register_key.reserve(records.size());
  for (const core::PendingFluxRegisterRecord& record : records) {
    register_key.push_back(record.register_key);
    coarse_patch_id.push_back(record.coarse_patch_id);
    coarse_gas_cell_id.push_back(record.coarse_gas_cell_id);
    coarse_cell_index.push_back(static_cast<std::uint64_t>(record.coarse_cell_index));
    level.push_back(record.level);
    axis.push_back(record.axis);
    orientation.push_back(record.orientation);
    expected_area_comov.push_back(record.expected_area_comov);
    coarse_area_accumulated_comov.push_back(record.coarse_area_accumulated_comov);
    fine_area_accumulated_comov.push_back(record.fine_area_accumulated_comov);
    interval_start_code.push_back(record.interval_start_code);
    interval_end_code.push_back(record.interval_end_code);
    coarse_dt_code.push_back(record.coarse_dt_code);
    expected_fine_substeps.push_back(record.expected_fine_substeps);
    completed_fine_substeps.push_back(record.completed_fine_substeps);
    fine_substep_coverage_mask.push_back(record.fine_substep_coverage_mask);
    coarse_face_count.push_back(record.coarse_face_count);
    fine_face_count.push_back(record.fine_face_count);
    gas_cell_identity_generation.push_back(record.gas_cell_identity_generation);
    patch_geometry_generation.push_back(record.patch_geometry_generation);
    coarse_mass_flux_integral_code.push_back(record.coarse_mass_flux_integral_code);
    coarse_momentum_x_flux_integral_code.push_back(record.coarse_momentum_x_flux_integral_code);
    coarse_momentum_y_flux_integral_code.push_back(record.coarse_momentum_y_flux_integral_code);
    coarse_momentum_z_flux_integral_code.push_back(record.coarse_momentum_z_flux_integral_code);
    coarse_total_energy_flux_integral_code.push_back(record.coarse_total_energy_flux_integral_code);
    coarse_metal_mass_flux_integral_code.push_back(record.coarse_metal_mass_flux_integral_code);
    fine_mass_flux_integral_code.push_back(record.fine_mass_flux_integral_code);
    fine_momentum_x_flux_integral_code.push_back(record.fine_momentum_x_flux_integral_code);
    fine_momentum_y_flux_integral_code.push_back(record.fine_momentum_y_flux_integral_code);
    fine_momentum_z_flux_integral_code.push_back(record.fine_momentum_z_flux_integral_code);
    fine_total_energy_flux_integral_code.push_back(record.fine_total_energy_flux_integral_code);
    fine_metal_mass_flux_integral_code.push_back(record.fine_metal_mass_flux_integral_code);
  }
  writeDataset1d(group.get(), "register_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, register_key);
  writeDataset1d(group.get(), "coarse_patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, coarse_patch_id);
  writeDataset1d(group.get(), "coarse_gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, coarse_gas_cell_id);
  writeDataset1d(group.get(), "coarse_cell_index", H5T_STD_U64LE, H5T_NATIVE_UINT64, coarse_cell_index);
  writeDataset1d(group.get(), "level", H5T_STD_U8LE, H5T_NATIVE_UINT8, level);
  writeDataset1d(group.get(), "axis", H5T_STD_U8LE, H5T_NATIVE_UINT8, axis);
  writeDataset1d(group.get(), "orientation", H5T_STD_U8LE, H5T_NATIVE_UINT8, orientation);
  writeDataset1d(group.get(), "expected_area_comov", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, expected_area_comov);
  writeDataset1d(group.get(), "coarse_area_accumulated_comov", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_area_accumulated_comov);
  writeDataset1d(group.get(), "fine_area_accumulated_comov", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_area_accumulated_comov);
  writeDataset1d(group.get(), "interval_start_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, interval_start_code);
  writeDataset1d(group.get(), "interval_end_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, interval_end_code);
  writeDataset1d(group.get(), "coarse_dt_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_dt_code);
  writeDataset1d(group.get(), "expected_fine_substeps", H5T_STD_U32LE, H5T_NATIVE_UINT32, expected_fine_substeps);
  writeDataset1d(group.get(), "completed_fine_substeps", H5T_STD_U32LE, H5T_NATIVE_UINT32, completed_fine_substeps);
  writeDataset1d(group.get(), "fine_substep_coverage_mask", H5T_STD_U64LE, H5T_NATIVE_UINT64, fine_substep_coverage_mask);
  writeDataset1d(group.get(), "coarse_face_count", H5T_STD_U32LE, H5T_NATIVE_UINT32, coarse_face_count);
  writeDataset1d(group.get(), "fine_face_count", H5T_STD_U32LE, H5T_NATIVE_UINT32, fine_face_count);
  writeDataset1d(group.get(), "gas_cell_identity_generation", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_identity_generation);
  writeDataset1d(group.get(), "patch_geometry_generation", H5T_STD_U64LE, H5T_NATIVE_UINT64, patch_geometry_generation);
  writeDataset1d(group.get(), "coarse_mass_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_mass_flux_integral_code);
  writeDataset1d(group.get(), "coarse_momentum_x_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_momentum_x_flux_integral_code);
  writeDataset1d(group.get(), "coarse_momentum_y_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_momentum_y_flux_integral_code);
  writeDataset1d(group.get(), "coarse_momentum_z_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_momentum_z_flux_integral_code);
  writeDataset1d(group.get(), "coarse_total_energy_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_total_energy_flux_integral_code);
  writeDataset1d(group.get(), "coarse_metal_mass_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_metal_mass_flux_integral_code);
  writeDataset1d(group.get(), "fine_mass_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_mass_flux_integral_code);
  writeDataset1d(group.get(), "fine_momentum_x_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_momentum_x_flux_integral_code);
  writeDataset1d(group.get(), "fine_momentum_y_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_momentum_y_flux_integral_code);
  writeDataset1d(group.get(), "fine_momentum_z_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_momentum_z_flux_integral_code);
  writeDataset1d(group.get(), "fine_total_energy_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_total_energy_flux_integral_code);
  writeDataset1d(group.get(), "fine_metal_mass_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_metal_mass_flux_integral_code);
}

void readPendingFluxRegisterGroup(hid_t state_group, core::SimulationState& state, std::uint32_t schema_version) {
  if (schema_version < k_restart_schema_v19 && H5Lexists(state_group, "amr_pending_flux_registers", H5P_DEFAULT) <= 0) {
    state.pending_flux_registers.clear();
    return;
  }
  Hdf5Handle group(H5Gopen2(state_group, "amr_pending_flux_registers", H5P_DEFAULT));
  if (!group.valid()) {
    state.pending_flux_registers.clear();
    return;
  }
  const auto register_key = readDataset1d<std::uint64_t>(group.get(), "register_key", H5T_NATIVE_UINT64);
  const auto coarse_patch_id = readDataset1d<std::uint64_t>(group.get(), "coarse_patch_id", H5T_NATIVE_UINT64);
  const auto coarse_gas_cell_id = readDataset1d<std::uint64_t>(group.get(), "coarse_gas_cell_id", H5T_NATIVE_UINT64);
  const auto coarse_cell_index = readDataset1d<std::uint64_t>(group.get(), "coarse_cell_index", H5T_NATIVE_UINT64);
  const auto level = readDataset1d<std::uint8_t>(group.get(), "level", H5T_NATIVE_UINT8);
  const auto axis = readDataset1d<std::uint8_t>(group.get(), "axis", H5T_NATIVE_UINT8);
  const auto orientation = readDataset1d<std::uint8_t>(group.get(), "orientation", H5T_NATIVE_UINT8);
  const auto expected_area_comov = readDataset1d<double>(group.get(), "expected_area_comov", H5T_NATIVE_DOUBLE);
  const auto coarse_area_accumulated_comov = readDataset1d<double>(group.get(), "coarse_area_accumulated_comov", H5T_NATIVE_DOUBLE);
  const auto fine_area_accumulated_comov = readDataset1d<double>(group.get(), "fine_area_accumulated_comov", H5T_NATIVE_DOUBLE);
  const auto interval_start_code = readDataset1d<double>(group.get(), "interval_start_code", H5T_NATIVE_DOUBLE);
  const auto interval_end_code = readDataset1d<double>(group.get(), "interval_end_code", H5T_NATIVE_DOUBLE);
  const auto coarse_dt_code = readDataset1d<double>(group.get(), "coarse_dt_code", H5T_NATIVE_DOUBLE);
  const auto expected_fine_substeps = readDataset1d<std::uint32_t>(group.get(), "expected_fine_substeps", H5T_NATIVE_UINT32);
  const auto completed_fine_substeps = readDataset1d<std::uint32_t>(group.get(), "completed_fine_substeps", H5T_NATIVE_UINT32);
  const auto fine_substep_coverage_mask = readDataset1d<std::uint64_t>(group.get(), "fine_substep_coverage_mask", H5T_NATIVE_UINT64);
  const auto coarse_face_count = readDataset1d<std::uint32_t>(group.get(), "coarse_face_count", H5T_NATIVE_UINT32);
  const auto fine_face_count = readDataset1d<std::uint32_t>(group.get(), "fine_face_count", H5T_NATIVE_UINT32);
  const auto gas_cell_identity_generation = readDataset1d<std::uint64_t>(group.get(), "gas_cell_identity_generation", H5T_NATIVE_UINT64);
  const auto patch_geometry_generation = readDataset1d<std::uint64_t>(group.get(), "patch_geometry_generation", H5T_NATIVE_UINT64);
  const auto coarse_mass_flux_integral_code = readDataset1d<double>(group.get(), "coarse_mass_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto coarse_momentum_x_flux_integral_code = readDataset1d<double>(group.get(), "coarse_momentum_x_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto coarse_momentum_y_flux_integral_code = readDataset1d<double>(group.get(), "coarse_momentum_y_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto coarse_momentum_z_flux_integral_code = readDataset1d<double>(group.get(), "coarse_momentum_z_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto coarse_total_energy_flux_integral_code = readDataset1d<double>(group.get(), "coarse_total_energy_flux_integral_code", H5T_NATIVE_DOUBLE);
  const std::size_t n = register_key.size();
  const auto coarse_metal_mass_flux_integral_code =
      H5Lexists(group.get(), "coarse_metal_mass_flux_integral_code", H5P_DEFAULT) > 0
          ? readDataset1d<double>(group.get(), "coarse_metal_mass_flux_integral_code", H5T_NATIVE_DOUBLE)
          : std::vector<double>(n, 0.0);
  const auto fine_mass_flux_integral_code = readDataset1d<double>(group.get(), "fine_mass_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_momentum_x_flux_integral_code = readDataset1d<double>(group.get(), "fine_momentum_x_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_momentum_y_flux_integral_code = readDataset1d<double>(group.get(), "fine_momentum_y_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_momentum_z_flux_integral_code = readDataset1d<double>(group.get(), "fine_momentum_z_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_total_energy_flux_integral_code = readDataset1d<double>(group.get(), "fine_total_energy_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_metal_mass_flux_integral_code =
      H5Lexists(group.get(), "fine_metal_mass_flux_integral_code", H5P_DEFAULT) > 0
          ? readDataset1d<double>(group.get(), "fine_metal_mass_flux_integral_code", H5T_NATIVE_DOUBLE)
          : std::vector<double>(n, 0.0);
  const bool sizes_match = coarse_patch_id.size() == n && coarse_gas_cell_id.size() == n && coarse_cell_index.size() == n &&
      level.size() == n && axis.size() == n && orientation.size() == n && expected_area_comov.size() == n &&
      coarse_area_accumulated_comov.size() == n && fine_area_accumulated_comov.size() == n && interval_start_code.size() == n &&
      interval_end_code.size() == n && coarse_dt_code.size() == n && expected_fine_substeps.size() == n &&
      completed_fine_substeps.size() == n && fine_substep_coverage_mask.size() == n && coarse_face_count.size() == n &&
      fine_face_count.size() == n && gas_cell_identity_generation.size() == n && patch_geometry_generation.size() == n &&
      coarse_mass_flux_integral_code.size() == n && coarse_momentum_x_flux_integral_code.size() == n &&
      coarse_momentum_y_flux_integral_code.size() == n && coarse_momentum_z_flux_integral_code.size() == n &&
      coarse_total_energy_flux_integral_code.size() == n && coarse_metal_mass_flux_integral_code.size() == n &&
      fine_mass_flux_integral_code.size() == n &&
      fine_momentum_x_flux_integral_code.size() == n && fine_momentum_y_flux_integral_code.size() == n &&
      fine_momentum_z_flux_integral_code.size() == n && fine_total_energy_flux_integral_code.size() == n &&
      fine_metal_mass_flux_integral_code.size() == n;
  if (!sizes_match) {
    throw std::runtime_error("/state/amr_pending_flux_registers datasets must have matching lengths");
  }
  std::vector<core::PendingFluxRegisterRecord> records;
  records.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    records.push_back(core::PendingFluxRegisterRecord{
        .register_key = register_key[i],
        .coarse_patch_id = coarse_patch_id[i],
        .coarse_gas_cell_id = coarse_gas_cell_id[i],
        .coarse_cell_index = static_cast<std::size_t>(coarse_cell_index[i]),
        .level = level[i],
        .axis = axis[i],
        .orientation = orientation[i],
        .expected_area_comov = expected_area_comov[i],
        .coarse_area_accumulated_comov = coarse_area_accumulated_comov[i],
        .fine_area_accumulated_comov = fine_area_accumulated_comov[i],
        .interval_start_code = interval_start_code[i],
        .interval_end_code = interval_end_code[i],
        .coarse_dt_code = coarse_dt_code[i],
        .expected_fine_substeps = expected_fine_substeps[i],
        .completed_fine_substeps = completed_fine_substeps[i],
        .fine_substep_coverage_mask = fine_substep_coverage_mask[i],
        .coarse_face_count = coarse_face_count[i],
        .fine_face_count = fine_face_count[i],
        .gas_cell_identity_generation = gas_cell_identity_generation[i],
        .patch_geometry_generation = patch_geometry_generation[i],
        .coarse_mass_flux_integral_code = coarse_mass_flux_integral_code[i],
        .coarse_momentum_x_flux_integral_code = coarse_momentum_x_flux_integral_code[i],
        .coarse_momentum_y_flux_integral_code = coarse_momentum_y_flux_integral_code[i],
        .coarse_momentum_z_flux_integral_code = coarse_momentum_z_flux_integral_code[i],
        .coarse_total_energy_flux_integral_code = coarse_total_energy_flux_integral_code[i],
        .coarse_metal_mass_flux_integral_code = coarse_metal_mass_flux_integral_code[i],
        .fine_mass_flux_integral_code = fine_mass_flux_integral_code[i],
        .fine_momentum_x_flux_integral_code = fine_momentum_x_flux_integral_code[i],
        .fine_momentum_y_flux_integral_code = fine_momentum_y_flux_integral_code[i],
        .fine_momentum_z_flux_integral_code = fine_momentum_z_flux_integral_code[i],
        .fine_total_energy_flux_integral_code = fine_total_energy_flux_integral_code[i],
        .fine_metal_mass_flux_integral_code = fine_metal_mass_flux_integral_code[i]});
  }
  state.pending_flux_registers.assign(std::move(records));
}


void writeAmrTemporalBoundaryHistoryGroup(
    hid_t state_group,
    const core::AmrTemporalBoundaryHistoryStore& store) {
  Hdf5Handle group(openOrCreateGroup(state_group, "amr_temporal_boundary_history"));
  writeScalarU32Attribute(group.get(), "schema_version", 2U);
  std::vector<std::uint64_t> patch_id;
  std::vector<std::uint8_t> patch_level;
  std::vector<std::uint64_t> patch_geometry_fingerprint;
  std::vector<std::uint64_t> gas_cell_identity_generation;
  std::vector<double> interval_start_code;
  std::vector<double> interval_end_code;
  std::vector<std::uint8_t> end_state_valid;
  std::vector<std::uint64_t> cell_offset;
  std::vector<std::uint64_t> cell_count;
  std::vector<std::uint64_t> gas_cell_id;
  std::vector<std::uint64_t> patch_local_cell;
  std::vector<double> start_mass_density_comoving;
  std::vector<double> start_momentum_density_x_comoving;
  std::vector<double> start_momentum_density_y_comoving;
  std::vector<double> start_momentum_density_z_comoving;
  std::vector<double> start_total_energy_density_comoving;
  std::vector<double> start_metal_mass_density_comoving;
  std::vector<double> end_mass_density_comoving;
  std::vector<double> end_momentum_density_x_comoving;
  std::vector<double> end_momentum_density_y_comoving;
  std::vector<double> end_momentum_density_z_comoving;
  std::vector<double> end_total_energy_density_comoving;
  std::vector<double> end_metal_mass_density_comoving;
  const auto records = store.records();
  patch_id.reserve(records.size());
  std::uint64_t offset = 0U;
  for (const core::AmrTemporalBoundaryHistoryRecord& record : records) {
    patch_id.push_back(record.patch_id);
    patch_level.push_back(record.patch_level);
    patch_geometry_fingerprint.push_back(record.patch_geometry_fingerprint);
    gas_cell_identity_generation.push_back(record.gas_cell_identity_generation);
    interval_start_code.push_back(record.interval_start_code);
    interval_end_code.push_back(record.interval_end_code);
    end_state_valid.push_back(record.end_state_valid ? 1U : 0U);
    cell_offset.push_back(offset);
    cell_count.push_back(static_cast<std::uint64_t>(record.cells.size()));
    for (const core::AmrTemporalBoundaryHistoryCellRecord& cell : record.cells) {
      gas_cell_id.push_back(cell.gas_cell_id);
      patch_local_cell.push_back(static_cast<std::uint64_t>(cell.patch_local_cell));
      start_mass_density_comoving.push_back(cell.start_mass_density_comoving);
      start_momentum_density_x_comoving.push_back(cell.start_momentum_density_x_comoving);
      start_momentum_density_y_comoving.push_back(cell.start_momentum_density_y_comoving);
      start_momentum_density_z_comoving.push_back(cell.start_momentum_density_z_comoving);
      start_total_energy_density_comoving.push_back(cell.start_total_energy_density_comoving);
      start_metal_mass_density_comoving.push_back(cell.start_metal_mass_density_comoving);
      end_mass_density_comoving.push_back(cell.end_mass_density_comoving);
      end_momentum_density_x_comoving.push_back(cell.end_momentum_density_x_comoving);
      end_momentum_density_y_comoving.push_back(cell.end_momentum_density_y_comoving);
      end_momentum_density_z_comoving.push_back(cell.end_momentum_density_z_comoving);
      end_total_energy_density_comoving.push_back(cell.end_total_energy_density_comoving);
      end_metal_mass_density_comoving.push_back(cell.end_metal_mass_density_comoving);
      ++offset;
    }
  }
  writeDataset1d(group.get(), "patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, patch_id);
  writeDataset1d(group.get(), "patch_level", H5T_STD_U8LE, H5T_NATIVE_UINT8, patch_level);
  writeDataset1d(group.get(), "patch_geometry_fingerprint", H5T_STD_U64LE, H5T_NATIVE_UINT64, patch_geometry_fingerprint);
  writeDataset1d(group.get(), "gas_cell_identity_generation", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_identity_generation);
  writeDataset1d(group.get(), "interval_start_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, interval_start_code);
  writeDataset1d(group.get(), "interval_end_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, interval_end_code);
  writeDataset1d(group.get(), "end_state_valid", H5T_STD_U8LE, H5T_NATIVE_UINT8, end_state_valid);
  writeDataset1d(group.get(), "cell_offset", H5T_STD_U64LE, H5T_NATIVE_UINT64, cell_offset);
  writeDataset1d(group.get(), "cell_count", H5T_STD_U64LE, H5T_NATIVE_UINT64, cell_count);
  writeDataset1d(group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_id);
  writeDataset1d(group.get(), "patch_local_cell", H5T_STD_U64LE, H5T_NATIVE_UINT64, patch_local_cell);
  writeDataset1d(group.get(), "start_mass_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_mass_density_comoving);
  writeDataset1d(group.get(), "start_momentum_density_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_momentum_density_x_comoving);
  writeDataset1d(group.get(), "start_momentum_density_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_momentum_density_y_comoving);
  writeDataset1d(group.get(), "start_momentum_density_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_momentum_density_z_comoving);
  writeDataset1d(group.get(), "start_total_energy_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_total_energy_density_comoving);
  writeDataset1d(group.get(), "start_metal_mass_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_metal_mass_density_comoving);
  writeDataset1d(group.get(), "end_mass_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_mass_density_comoving);
  writeDataset1d(group.get(), "end_momentum_density_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_momentum_density_x_comoving);
  writeDataset1d(group.get(), "end_momentum_density_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_momentum_density_y_comoving);
  writeDataset1d(group.get(), "end_momentum_density_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_momentum_density_z_comoving);
  writeDataset1d(group.get(), "end_total_energy_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_total_energy_density_comoving);
  writeDataset1d(group.get(), "end_metal_mass_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_metal_mass_density_comoving);
}

void readAmrTemporalBoundaryHistoryGroup(
    hid_t state_group,
    core::SimulationState& state,
    std::uint32_t schema_version) {
  if (schema_version < k_restart_schema_v19) {
    // Older checkpoints cannot prove a mid-subcycle continuation has the
    // time-aligned coarse history required by temporal coarse-to-fine fills.
    // A legacy file is therefore safe only when no deferred AMR synchronization
    // state remains. Pending registers are the conservative observable marker
    // for an open local subcycling interval in pre-v18 checkpoints.
    if (!state.pending_flux_registers.empty()) {
      throw std::runtime_error(
          "restart reader: legacy restart with pending AMR flux registers cannot resume "
          "without v18 temporal boundary history");
    }
    state.amr_temporal_boundary_history.clear();
    return;
  }
  Hdf5Handle group(H5Gopen2(state_group, "amr_temporal_boundary_history", H5P_DEFAULT));
  if (!group.valid()) {
    throw std::runtime_error("restart reader: v18 restart is missing /state/amr_temporal_boundary_history");
  }
  const std::uint32_t temporal_schema_version = readScalarU32Attribute(group.get(), "schema_version");
  if (temporal_schema_version != 1U && temporal_schema_version != 2U) {
    throw std::runtime_error("restart reader: unsupported /state/amr_temporal_boundary_history schema_version");
  }
  const auto patch_id = readDataset1d<std::uint64_t>(group.get(), "patch_id", H5T_NATIVE_UINT64);
  const auto patch_level = readDataset1d<std::uint8_t>(group.get(), "patch_level", H5T_NATIVE_UINT8);
  const auto patch_geometry_fingerprint = readDataset1d<std::uint64_t>(group.get(), "patch_geometry_fingerprint", H5T_NATIVE_UINT64);
  const auto gas_cell_identity_generation = readDataset1d<std::uint64_t>(group.get(), "gas_cell_identity_generation", H5T_NATIVE_UINT64);
  const auto interval_start_code = readDataset1d<double>(group.get(), "interval_start_code", H5T_NATIVE_DOUBLE);
  const auto interval_end_code = readDataset1d<double>(group.get(), "interval_end_code", H5T_NATIVE_DOUBLE);
  const auto end_state_valid = readDataset1d<std::uint8_t>(group.get(), "end_state_valid", H5T_NATIVE_UINT8);
  const auto cell_offset = readDataset1d<std::uint64_t>(group.get(), "cell_offset", H5T_NATIVE_UINT64);
  const auto cell_count = readDataset1d<std::uint64_t>(group.get(), "cell_count", H5T_NATIVE_UINT64);
  const auto gas_cell_id = readDataset1d<std::uint64_t>(group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
  const auto patch_local_cell = readDataset1d<std::uint64_t>(group.get(), "patch_local_cell", H5T_NATIVE_UINT64);
  const auto start_mass_density_comoving = readDataset1d<double>(group.get(), "start_mass_density_comoving", H5T_NATIVE_DOUBLE);
  const auto start_momentum_density_x_comoving = readDataset1d<double>(group.get(), "start_momentum_density_x_comoving", H5T_NATIVE_DOUBLE);
  const auto start_momentum_density_y_comoving = readDataset1d<double>(group.get(), "start_momentum_density_y_comoving", H5T_NATIVE_DOUBLE);
  const auto start_momentum_density_z_comoving = readDataset1d<double>(group.get(), "start_momentum_density_z_comoving", H5T_NATIVE_DOUBLE);
  const auto start_total_energy_density_comoving = readDataset1d<double>(group.get(), "start_total_energy_density_comoving", H5T_NATIVE_DOUBLE);
  const std::size_t nc = gas_cell_id.size();
  const auto start_metal_mass_density_comoving =
      H5Lexists(group.get(), "start_metal_mass_density_comoving", H5P_DEFAULT) > 0
          ? readDataset1d<double>(group.get(), "start_metal_mass_density_comoving", H5T_NATIVE_DOUBLE)
          : std::vector<double>(nc, 0.0);
  const auto end_mass_density_comoving = readDataset1d<double>(group.get(), "end_mass_density_comoving", H5T_NATIVE_DOUBLE);
  const auto end_momentum_density_x_comoving = readDataset1d<double>(group.get(), "end_momentum_density_x_comoving", H5T_NATIVE_DOUBLE);
  const auto end_momentum_density_y_comoving = readDataset1d<double>(group.get(), "end_momentum_density_y_comoving", H5T_NATIVE_DOUBLE);
  const auto end_momentum_density_z_comoving = readDataset1d<double>(group.get(), "end_momentum_density_z_comoving", H5T_NATIVE_DOUBLE);
  const auto end_total_energy_density_comoving = readDataset1d<double>(group.get(), "end_total_energy_density_comoving", H5T_NATIVE_DOUBLE);
  const auto end_metal_mass_density_comoving =
      H5Lexists(group.get(), "end_metal_mass_density_comoving", H5P_DEFAULT) > 0
          ? readDataset1d<double>(group.get(), "end_metal_mass_density_comoving", H5T_NATIVE_DOUBLE)
          : std::vector<double>(nc, 0.0);
  const std::size_t n = patch_id.size();
  const bool record_sizes = patch_level.size() == n && patch_geometry_fingerprint.size() == n &&
      gas_cell_identity_generation.size() == n && interval_start_code.size() == n && interval_end_code.size() == n &&
      end_state_valid.size() == n && cell_offset.size() == n && cell_count.size() == n;
  const bool cell_sizes = patch_local_cell.size() == nc && start_mass_density_comoving.size() == nc &&
      start_momentum_density_x_comoving.size() == nc && start_momentum_density_y_comoving.size() == nc &&
      start_momentum_density_z_comoving.size() == nc && start_total_energy_density_comoving.size() == nc &&
      start_metal_mass_density_comoving.size() == nc && end_mass_density_comoving.size() == nc && end_momentum_density_x_comoving.size() == nc &&
      end_momentum_density_y_comoving.size() == nc && end_momentum_density_z_comoving.size() == nc &&
      end_total_energy_density_comoving.size() == nc && end_metal_mass_density_comoving.size() == nc;
  if (!record_sizes || !cell_sizes) {
    throw std::runtime_error("/state/amr_temporal_boundary_history datasets must have matching lengths");
  }
  std::vector<core::AmrTemporalBoundaryHistoryRecord> records;
  records.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (end_state_valid[i] > 1U || cell_offset[i] > nc || cell_count[i] > nc - cell_offset[i] ||
        !std::isfinite(interval_start_code[i]) || !std::isfinite(interval_end_code[i]) ||
        interval_end_code[i] <= interval_start_code[i]) {
      throw std::runtime_error("/state/amr_temporal_boundary_history has invalid interval or record bounds");
    }
    core::AmrTemporalBoundaryHistoryRecord record;
    record.patch_id = patch_id[i];
    record.patch_level = patch_level[i];
    record.patch_geometry_fingerprint = patch_geometry_fingerprint[i];
    record.gas_cell_identity_generation = gas_cell_identity_generation[i];
    record.interval_start_code = interval_start_code[i];
    record.interval_end_code = interval_end_code[i];
    record.end_state_valid = end_state_valid[i] != 0U;
    const std::size_t begin = static_cast<std::size_t>(cell_offset[i]);
    const std::size_t count = static_cast<std::size_t>(cell_count[i]);
    record.cells.reserve(count);
    for (std::size_t j = 0; j < count; ++j) {
      const std::size_t k = begin + j;
      record.cells.push_back(core::AmrTemporalBoundaryHistoryCellRecord{
          .gas_cell_id = gas_cell_id[k],
          .patch_local_cell = static_cast<std::size_t>(patch_local_cell[k]),
          .start_mass_density_comoving = start_mass_density_comoving[k],
          .start_momentum_density_x_comoving = start_momentum_density_x_comoving[k],
          .start_momentum_density_y_comoving = start_momentum_density_y_comoving[k],
          .start_momentum_density_z_comoving = start_momentum_density_z_comoving[k],
          .start_total_energy_density_comoving = start_total_energy_density_comoving[k],
          .start_metal_mass_density_comoving = start_metal_mass_density_comoving[k],
          .end_mass_density_comoving = end_mass_density_comoving[k],
          .end_momentum_density_x_comoving = end_momentum_density_x_comoving[k],
          .end_momentum_density_y_comoving = end_momentum_density_y_comoving[k],
          .end_momentum_density_z_comoving = end_momentum_density_z_comoving[k],
          .end_total_energy_density_comoving = end_total_energy_density_comoving[k],
          .end_metal_mass_density_comoving = end_metal_mass_density_comoving[k]});
    }
    records.push_back(std::move(record));
  }
  state.amr_temporal_boundary_history.assign(std::move(records));
}

void validateAmrTemporalBoundaryHistoryForRestart(const core::SimulationState& state) {
  for (const core::AmrTemporalBoundaryHistoryRecord& record : state.amr_temporal_boundary_history.records()) {
    if (record.patch_id == 0U || record.cells.empty() || !record.end_state_valid ||
        record.gas_cell_identity_generation != state.gasCellIdentityGeneration() ||
        !std::isfinite(record.interval_start_code) || !std::isfinite(record.interval_end_code) ||
        record.interval_end_code <= record.interval_start_code) {
      throw std::runtime_error("restart reader: active AMR temporal history is incomplete or stale");
    }
    for (std::size_t i = 0; i < record.cells.size(); ++i) {
      const core::AmrTemporalBoundaryHistoryCellRecord& cell = record.cells[i];
      const auto finite_state = [](double value) { return std::isfinite(value); };
      if (cell.gas_cell_id == 0U || !finite_state(cell.start_mass_density_comoving) ||
          !finite_state(cell.start_momentum_density_x_comoving) ||
          !finite_state(cell.start_momentum_density_y_comoving) ||
          !finite_state(cell.start_momentum_density_z_comoving) ||
          !finite_state(cell.start_total_energy_density_comoving) ||
          !finite_state(cell.start_metal_mass_density_comoving) ||
          !finite_state(cell.end_mass_density_comoving) ||
          !finite_state(cell.end_momentum_density_x_comoving) ||
          !finite_state(cell.end_momentum_density_y_comoving) ||
          !finite_state(cell.end_momentum_density_z_comoving) ||
          !finite_state(cell.end_total_energy_density_comoving) ||
          !finite_state(cell.end_metal_mass_density_comoving) ||
          cell.start_mass_density_comoving <= 0.0 || cell.end_mass_density_comoving <= 0.0 ||
          cell.start_metal_mass_density_comoving < 0.0 ||
          cell.end_metal_mass_density_comoving < 0.0 ||
          cell.start_metal_mass_density_comoving > cell.start_mass_density_comoving ||
          cell.end_metal_mass_density_comoving > cell.end_mass_density_comoving) {
        throw std::runtime_error("restart reader: temporal history contains inadmissible conserved state");
      }
      for (std::size_t j = 0; j < i; ++j) {
        if (record.cells[j].gas_cell_id == cell.gas_cell_id ||
            record.cells[j].patch_local_cell == cell.patch_local_cell) {
          throw std::runtime_error("restart reader: temporal history has duplicate stable cell identity or patch-local index");
        }
      }
      const auto row = state.rowForGasCellId(cell.gas_cell_id);
      if (!row.has_value() || state.cells.patch_index[*row] >= state.patches.size() ||
          state.patches.patch_id[state.cells.patch_index[*row]] != record.patch_id) {
        throw std::runtime_error("restart reader: AMR temporal history references a missing or moved coarse gas cell");
      }
    }
  }
}

void writeStateGroup(hid_t root, const core::SimulationState& state) {
  Hdf5Handle state_group(openOrCreateGroup(root, "/state"));
  Hdf5Handle particles_group(openOrCreateGroup(state_group.get(), "particles"));
  writeDataset1d(particles_group.get(), "position_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.position_x_comoving);
  writeDataset1d(particles_group.get(), "position_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.position_y_comoving);
  writeDataset1d(particles_group.get(), "position_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.position_z_comoving);
  writeDataset1d(particles_group.get(), "velocity_x_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.velocity_x_peculiar);
  writeDataset1d(particles_group.get(), "velocity_y_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.velocity_y_peculiar);
  writeDataset1d(particles_group.get(), "velocity_z_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.velocity_z_peculiar);
  writeDataset1d(particles_group.get(), "mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.mass_code);
  if (state.hasHomogeneousDmoMetadata() && state.particles.time_bin.empty()) {
    writeUniformDataset1d(particles_group.get(), "time_bin", H5T_STD_U8LE, H5T_NATIVE_UINT8, state.particles.size(), static_cast<std::uint8_t>(0U));
  } else {
    writeDataset1d(particles_group.get(), "time_bin", H5T_STD_U8LE, H5T_NATIVE_UINT8, state.particles.time_bin);
  }

  Hdf5Handle particle_sidecar_group(openOrCreateGroup(state_group.get(), "particle_sidecar"));
  writeDataset1d(particle_sidecar_group.get(), "particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.particle_sidecar.particle_id);
  if (state.hasHomogeneousDmoMetadata()) {
    const std::size_t n = state.particles.size();
    const std::size_t probe = 0U;
    if (state.particle_sidecar.sfc_key.empty()) {
      writeUniformDataset1d(particle_sidecar_group.get(), "sfc_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, n, n == 0U ? 0ULL : state.particleSfcKey(probe));
    } else {
      writeDataset1d(particle_sidecar_group.get(), "sfc_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.particle_sidecar.sfc_key);
    }
    writeUniformDataset1d(particle_sidecar_group.get(), "species_tag", H5T_STD_U32LE, H5T_NATIVE_UINT32, n, n == 0U ? static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter) : state.particleSpeciesTag(probe));
    writeUniformDataset1d(particle_sidecar_group.get(), "particle_flags", H5T_STD_U32LE, H5T_NATIVE_UINT32, n, n == 0U ? 0U : state.particleFlags(probe));
    writeUniformDataset1d(particle_sidecar_group.get(), "owning_rank", H5T_STD_U32LE, H5T_NATIVE_UINT32, n, n == 0U ? 0U : state.particleOwningRank(probe));
    writeUniformDataset1d(particle_sidecar_group.get(), "last_drift_time_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, n, n == 0U ? 0.0 : state.particleLastDriftTimeCode(probe));
    writeUniformDataset1d(particle_sidecar_group.get(), "last_drift_scale_factor", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, n, n == 0U ? 1.0 : state.particleLastDriftScaleFactor(probe));
  } else {
    writeDataset1d(particle_sidecar_group.get(), "sfc_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.particle_sidecar.sfc_key);
    writeDataset1d(particle_sidecar_group.get(), "species_tag", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.particle_sidecar.species_tag);
    writeDataset1d(particle_sidecar_group.get(), "particle_flags", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.particle_sidecar.particle_flags);
    writeDataset1d(particle_sidecar_group.get(), "owning_rank", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.particle_sidecar.owning_rank);
    writeDataset1d(particle_sidecar_group.get(), "last_drift_time_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particle_sidecar.last_drift_time_code);
    writeDataset1d(particle_sidecar_group.get(), "last_drift_scale_factor", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particle_sidecar.last_drift_scale_factor);
  }
  writeDataset1d(
      particle_sidecar_group.get(),
      "gravity_softening_comoving",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.particle_sidecar.gravity_softening_comoving);
  writeDataset1d(
      particle_sidecar_group.get(),
      "has_gravity_softening_override",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      state.particle_sidecar.has_gravity_softening_override);

  Hdf5Handle cells_group(openOrCreateGroup(state_group.get(), "cells"));
  writeDataset1d(cells_group.get(), "center_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.center_x_comoving);
  writeDataset1d(cells_group.get(), "center_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.center_y_comoving);
  writeDataset1d(cells_group.get(), "center_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.center_z_comoving);
  writeDataset1d(cells_group.get(), "mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.mass_code);
  writeDataset1d(cells_group.get(), "time_bin", H5T_STD_U8LE, H5T_NATIVE_UINT8, state.cells.time_bin);
  writeDataset1d(cells_group.get(), "patch_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.cells.patch_index);

  Hdf5Handle gas_group(openOrCreateGroup(state_group.get(), "gas_cells"));
  writeDataset1d(gas_group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.gas_cells.gas_cell_id);
  writeDataset1d(gas_group.get(), "parent_particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.gas_cells.parent_particle_id);
  writeDataset1d(gas_group.get(), "velocity_x_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.velocity_x_peculiar);
  writeDataset1d(gas_group.get(), "velocity_y_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.velocity_y_peculiar);
  writeDataset1d(gas_group.get(), "velocity_z_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.velocity_z_peculiar);
  writeDataset1d(gas_group.get(), "density_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.density_code);
  writeDataset1d(gas_group.get(), "pressure_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.pressure_code);
  writeDataset1d(gas_group.get(), "internal_energy_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.internal_energy_code);
  writeDataset1d(gas_group.get(), "metal_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.metal_mass_code);
  writeDataset1d(gas_group.get(), "temperature_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.temperature_code);
  writeDataset1d(gas_group.get(), "sound_speed_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.sound_speed_code);

  writeGasCellIdentityGroup(state_group.get(), state);

  Hdf5Handle patches_group(openOrCreateGroup(state_group.get(), "patches"));
  writeDataset1d(patches_group.get(), "patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.patches.patch_id);
  writeDataset1d(patches_group.get(), "level", H5T_STD_I32LE, H5T_NATIVE_INT32, state.patches.level);
  writeDataset1d(patches_group.get(), "first_cell", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.patches.first_cell);
  writeDataset1d(patches_group.get(), "cell_count", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.patches.cell_count);
  writeDataset1d(patches_group.get(), "parent_patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.patches.parent_patch_id);
  writeDataset1d(patches_group.get(), "morton_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.patches.morton_key);
  writeDataset1d(patches_group.get(), "origin_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.origin_x_comoving);
  writeDataset1d(patches_group.get(), "origin_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.origin_y_comoving);
  writeDataset1d(patches_group.get(), "origin_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.origin_z_comoving);
  writeDataset1d(patches_group.get(), "extent_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.extent_x_comoving);
  writeDataset1d(patches_group.get(), "extent_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.extent_y_comoving);
  writeDataset1d(patches_group.get(), "extent_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.extent_z_comoving);
  writeDataset1d(patches_group.get(), "cell_dim_x", H5T_STD_U16LE, H5T_NATIVE_UINT16, state.patches.cell_dim_x);
  writeDataset1d(patches_group.get(), "cell_dim_y", H5T_STD_U16LE, H5T_NATIVE_UINT16, state.patches.cell_dim_y);
  writeDataset1d(patches_group.get(), "cell_dim_z", H5T_STD_U16LE, H5T_NATIVE_UINT16, state.patches.cell_dim_z);
  writeDataset1d(patches_group.get(), "owning_rank", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.patches.owning_rank);

  writePendingFluxRegisterGroup(state_group.get(), state.pending_flux_registers);
  writeAmrTemporalBoundaryHistoryGroup(state_group.get(), state.amr_temporal_boundary_history);

  writeDataset1d(state_group.get(), "species_count_by_species", H5T_STD_U64LE, H5T_NATIVE_UINT64, std::vector<std::uint64_t>(state.species.count_by_species.begin(), state.species.count_by_species.end()));

  writeStarSidecarGroup(state_group.get(), state.star_particles);

  Hdf5Handle bh_group(openOrCreateGroup(state_group.get(), "black_holes"));
  writeDataset1d(bh_group.get(), "particle_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.black_holes.particle_index);
  writeDataset1d(bh_group.get(), "host_cell_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.black_holes.host_cell_index);
  writeDataset1d(bh_group.get(), "subgrid_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.subgrid_mass_code);
  writeDataset1d(bh_group.get(), "accretion_rate_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.accretion_rate_code);
  writeDataset1d(bh_group.get(), "feedback_energy_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.feedback_energy_code);
  writeDataset1d(bh_group.get(), "eddington_ratio", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.eddington_ratio);
  writeDataset1d(
      bh_group.get(),
      "cumulative_accreted_mass_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.cumulative_accreted_mass_code);
  writeDataset1d(
      bh_group.get(),
      "cumulative_feedback_energy_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.cumulative_feedback_energy_code);
  writeDataset1d(
      bh_group.get(),
      "duty_cycle_active_time_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.duty_cycle_active_time_code);
  writeDataset1d(
      bh_group.get(),
      "duty_cycle_total_time_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.duty_cycle_total_time_code);

  Hdf5Handle tracer_group(openOrCreateGroup(state_group.get(), "tracers"));
  writeDataset1d(tracer_group.get(), "particle_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.tracers.particle_index);
  writeDataset1d(tracer_group.get(), "parent_particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.tracers.parent_particle_id);
  writeDataset1d(tracer_group.get(), "injection_step", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.tracers.injection_step);
  writeDataset1d(tracer_group.get(), "host_cell_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.tracers.host_cell_index);
  writeDataset1d(tracer_group.get(), "mass_fraction_of_host", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.tracers.mass_fraction_of_host);
  writeDataset1d(tracer_group.get(), "last_host_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.tracers.last_host_mass_code);
  writeDataset1d(
      tracer_group.get(),
      "cumulative_exchanged_mass_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.tracers.cumulative_exchanged_mass_code);

  writeStringDataset(state_group.get(), "state_metadata", state.metadata.serialize());

  Hdf5Handle sidecars_group(openOrCreateGroup(state_group.get(), "module_sidecars"));
  const auto sidecar_blocks = state.sidecars.blocksSortedByName();
  std::string module_names_text;
  for (const core::ModuleSidecarBlock* block : sidecar_blocks) {
    Hdf5Handle module_group(openOrCreateGroup(sidecars_group.get(), block->module_name));
    writeScalarU32Attribute(module_group.get(), "schema_version", block->schema_version);
    writeScalarU32Attribute(module_group.get(), "particle_indexed", block->particle_indexed ? 1U : 0U);
    writeScalarU32Attribute(module_group.get(), "row_stride_bytes", block->row_stride_bytes);
    writeScalarU32Attribute(module_group.get(), "required_species_mask", block->required_species_mask);
    writeScalarU32Attribute(module_group.get(), "requirement_kind", static_cast<std::uint32_t>(block->requirement.kind));
    writeScalarU32Attribute(module_group.get(), "requirement_species_mask", block->requirement.species_mask);
    writeScalarU32Attribute(module_group.get(), "requirement_particle_flags_mask", block->requirement.particle_flags_mask);
    writeScalarF64Attribute(module_group.get(), "requirement_threshold_code", block->requirement.threshold_code);
    // std::byte and uint8_t are both one-byte object representations. Write
    // the canonical payload directly; a second full sidecar is unnecessary.
    writeDataset1d(module_group.get(), "payload", H5T_STD_U8LE, H5T_NATIVE_UINT8,
                   std::span<const std::uint8_t>(
                       reinterpret_cast<const std::uint8_t*>(block->payload.data()),
                       block->payload.size()));
    if (block->isParticleIndexed()) {
      writeDataset1d(
          module_group.get(),
          "particle_id_by_row",
          H5T_STD_U64LE,
          H5T_NATIVE_UINT64,
          block->particle_id_by_row);
    }
    module_names_text += block->module_name;
    module_names_text.push_back('\n');
  }
  writeStringDataset(state_group.get(), "module_sidecar_names", module_names_text);
}

void readStateGroup(hid_t root, core::SimulationState& state, std::uint32_t schema_version) {
  Hdf5Handle state_group(H5Gopen2(root, "/state", H5P_DEFAULT));
  Hdf5Handle particles_group(H5Gopen2(state_group.get(), "particles", H5P_DEFAULT));
  state.particles.position_x_comoving = readDataset1dAligned<double>(particles_group.get(), "position_x_comoving", H5T_NATIVE_DOUBLE);
  state.particles.position_y_comoving = readDataset1dAligned<double>(particles_group.get(), "position_y_comoving", H5T_NATIVE_DOUBLE);
  state.particles.position_z_comoving = readDataset1dAligned<double>(particles_group.get(), "position_z_comoving", H5T_NATIVE_DOUBLE);
  state.particles.velocity_x_peculiar = readDataset1dAligned<double>(particles_group.get(), "velocity_x_peculiar", H5T_NATIVE_DOUBLE);
  state.particles.velocity_y_peculiar = readDataset1dAligned<double>(particles_group.get(), "velocity_y_peculiar", H5T_NATIVE_DOUBLE);
  state.particles.velocity_z_peculiar = readDataset1dAligned<double>(particles_group.get(), "velocity_z_peculiar", H5T_NATIVE_DOUBLE);
  state.particles.mass_code = readDataset1dAligned<double>(particles_group.get(), "mass_code", H5T_NATIVE_DOUBLE);
  state.particles.time_bin = readDataset1dAligned<std::uint8_t>(particles_group.get(), "time_bin", H5T_NATIVE_UINT8);

  Hdf5Handle particle_sidecar_group(H5Gopen2(state_group.get(), "particle_sidecar", H5P_DEFAULT));
  state.particle_sidecar.particle_id =
      readDataset1dAligned<std::uint64_t>(particle_sidecar_group.get(), "particle_id", H5T_NATIVE_UINT64);
  state.particle_sidecar.sfc_key =
      readDataset1dAligned<std::uint64_t>(particle_sidecar_group.get(), "sfc_key", H5T_NATIVE_UINT64);
  state.particle_sidecar.species_tag =
      readDataset1dAligned<std::uint32_t>(particle_sidecar_group.get(), "species_tag", H5T_NATIVE_UINT32);
  state.particle_sidecar.particle_flags =
      readDataset1dAligned<std::uint32_t>(particle_sidecar_group.get(), "particle_flags", H5T_NATIVE_UINT32);
  state.particle_sidecar.owning_rank =
      readDataset1dAligned<std::uint32_t>(particle_sidecar_group.get(), "owning_rank", H5T_NATIVE_UINT32);
  state.particle_sidecar.last_drift_time_code =
      readDataset1dAligned<double>(particle_sidecar_group.get(), "last_drift_time_code", H5T_NATIVE_DOUBLE);
  state.particle_sidecar.last_drift_scale_factor =
      readDataset1dAligned<double>(particle_sidecar_group.get(), "last_drift_scale_factor", H5T_NATIVE_DOUBLE);
  state.particle_sidecar.gravity_softening_comoving =
      readDataset1dAligned<double>(particle_sidecar_group.get(), "gravity_softening_comoving", H5T_NATIVE_DOUBLE);
  state.particle_sidecar.has_gravity_softening_override =
      readDataset1dAligned<std::uint8_t>(
          particle_sidecar_group.get(), "has_gravity_softening_override", H5T_NATIVE_UINT8);

  Hdf5Handle cells_group(H5Gopen2(state_group.get(), "cells", H5P_DEFAULT));
  state.cells.center_x_comoving = readDataset1dAligned<double>(cells_group.get(), "center_x_comoving", H5T_NATIVE_DOUBLE);
  state.cells.center_y_comoving = readDataset1dAligned<double>(cells_group.get(), "center_y_comoving", H5T_NATIVE_DOUBLE);
  state.cells.center_z_comoving = readDataset1dAligned<double>(cells_group.get(), "center_z_comoving", H5T_NATIVE_DOUBLE);
  state.cells.mass_code = readDataset1dAligned<double>(cells_group.get(), "mass_code", H5T_NATIVE_DOUBLE);
  state.cells.time_bin = readDataset1dAligned<std::uint8_t>(cells_group.get(), "time_bin", H5T_NATIVE_UINT8);
  state.cells.patch_index = readDataset1dAligned<std::uint32_t>(cells_group.get(), "patch_index", H5T_NATIVE_UINT32);

  Hdf5Handle gas_group(H5Gopen2(state_group.get(), "gas_cells", H5P_DEFAULT));
  if (H5Lexists(gas_group.get(), "gas_cell_id", H5P_DEFAULT) > 0 &&
      H5Lexists(gas_group.get(), "parent_particle_id", H5P_DEFAULT) > 0) {
    state.gas_cells.gas_cell_id =
        readDataset1dAligned<std::uint64_t>(gas_group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
    state.gas_cells.parent_particle_id =
        readDataset1dAligned<std::uint64_t>(gas_group.get(), "parent_particle_id", H5T_NATIVE_UINT64);
  } else {
    throw std::runtime_error("restart is missing gas-cell identity datasets gas_cell_id/parent_particle_id");
  }
  state.gas_cells.density_code = readDataset1dAligned<double>(gas_group.get(), "density_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.velocity_x_peculiar =
      readDataset1dAligned<double>(gas_group.get(), "velocity_x_peculiar", H5T_NATIVE_DOUBLE);
  state.gas_cells.velocity_y_peculiar =
      readDataset1dAligned<double>(gas_group.get(), "velocity_y_peculiar", H5T_NATIVE_DOUBLE);
  state.gas_cells.velocity_z_peculiar =
      readDataset1dAligned<double>(gas_group.get(), "velocity_z_peculiar", H5T_NATIVE_DOUBLE);
  state.gas_cells.pressure_code = readDataset1dAligned<double>(gas_group.get(), "pressure_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.internal_energy_code = readDataset1dAligned<double>(gas_group.get(), "internal_energy_code", H5T_NATIVE_DOUBLE);
  if (H5Lexists(gas_group.get(), "metal_mass_code", H5P_DEFAULT) > 0) {
    state.gas_cells.metal_mass_code =
        readDataset1dAligned<double>(gas_group.get(), "metal_mass_code", H5T_NATIVE_DOUBLE);
  } else {
    state.gas_cells.metal_mass_code.assign(state.cells.size(), 0.0);
  }
  state.gas_cells.temperature_code = readDataset1dAligned<double>(gas_group.get(), "temperature_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.sound_speed_code = readDataset1dAligned<double>(gas_group.get(), "sound_speed_code", H5T_NATIVE_DOUBLE);

  Hdf5Handle patches_group(H5Gopen2(state_group.get(), "patches", H5P_DEFAULT));
  state.patches.patch_id = readDataset1dAligned<std::uint64_t>(patches_group.get(), "patch_id", H5T_NATIVE_UINT64);
  state.patches.level = readDataset1dAligned<std::int32_t>(patches_group.get(), "level", H5T_NATIVE_INT32);
  state.patches.first_cell = readDataset1dAligned<std::uint32_t>(patches_group.get(), "first_cell", H5T_NATIVE_UINT32);
  state.patches.cell_count = readDataset1dAligned<std::uint32_t>(patches_group.get(), "cell_count", H5T_NATIVE_UINT32);
  if (H5Lexists(patches_group.get(), "parent_patch_id", H5P_DEFAULT) > 0) {
    state.patches.parent_patch_id = readDataset1dAligned<std::uint64_t>(patches_group.get(), "parent_patch_id", H5T_NATIVE_UINT64);
    state.patches.morton_key = readDataset1dAligned<std::uint64_t>(patches_group.get(), "morton_key", H5T_NATIVE_UINT64);
    state.patches.origin_x_comoving = readDataset1dAligned<double>(patches_group.get(), "origin_x_comoving", H5T_NATIVE_DOUBLE);
    state.patches.origin_y_comoving = readDataset1dAligned<double>(patches_group.get(), "origin_y_comoving", H5T_NATIVE_DOUBLE);
    state.patches.origin_z_comoving = readDataset1dAligned<double>(patches_group.get(), "origin_z_comoving", H5T_NATIVE_DOUBLE);
    state.patches.extent_x_comoving = readDataset1dAligned<double>(patches_group.get(), "extent_x_comoving", H5T_NATIVE_DOUBLE);
    state.patches.extent_y_comoving = readDataset1dAligned<double>(patches_group.get(), "extent_y_comoving", H5T_NATIVE_DOUBLE);
    state.patches.extent_z_comoving = readDataset1dAligned<double>(patches_group.get(), "extent_z_comoving", H5T_NATIVE_DOUBLE);
    state.patches.cell_dim_x = readDataset1dAligned<std::uint16_t>(patches_group.get(), "cell_dim_x", H5T_NATIVE_UINT16);
    state.patches.cell_dim_y = readDataset1dAligned<std::uint16_t>(patches_group.get(), "cell_dim_y", H5T_NATIVE_UINT16);
    state.patches.cell_dim_z = readDataset1dAligned<std::uint16_t>(patches_group.get(), "cell_dim_z", H5T_NATIVE_UINT16);
  } else {
    state.patches.parent_patch_id.assign(state.patches.patch_id.size(), 0U);
    state.patches.morton_key.assign(state.patches.patch_id.size(), 0U);
    state.patches.origin_x_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.origin_y_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.origin_z_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.extent_x_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.extent_y_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.extent_z_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.cell_dim_x.assign(state.patches.patch_id.size(), 0U);
    state.patches.cell_dim_y.assign(state.patches.patch_id.size(), 0U);
    state.patches.cell_dim_z.assign(state.patches.patch_id.size(), 0U);
  }
  if (H5Lexists(patches_group.get(), "owning_rank", H5P_DEFAULT) > 0) {
    state.patches.owning_rank = readDataset1dAligned<std::uint32_t>(patches_group.get(), "owning_rank", H5T_NATIVE_UINT32);
  } else {
    state.patches.owning_rank.assign(state.patches.patch_id.size(), 0U);
  }

  readGasCellIdentityGroup(state_group.get(), state, schema_version);
  readPendingFluxRegisterGroup(state_group.get(), state, schema_version);
  readAmrTemporalBoundaryHistoryGroup(state_group.get(), state, schema_version);

  const auto species_count = readDataset1d<std::uint64_t>(state_group.get(), "species_count_by_species", H5T_NATIVE_UINT64);
  if (species_count.size() != state.species.count_by_species.size()) {
    throw std::runtime_error("invalid species_count_by_species size in restart");
  }
  std::copy(species_count.begin(), species_count.end(), state.species.count_by_species.begin());

  readStarSidecarGroup(state_group.get(), state.star_particles);

  Hdf5Handle bh_group(H5Gopen2(state_group.get(), "black_holes", H5P_DEFAULT));
  state.black_holes.particle_index = readDataset1dAligned<std::uint32_t>(bh_group.get(), "particle_index", H5T_NATIVE_UINT32);
  state.black_holes.host_cell_index = readDataset1dAligned<std::uint32_t>(bh_group.get(), "host_cell_index", H5T_NATIVE_UINT32);
  state.black_holes.subgrid_mass_code = readDataset1dAligned<double>(bh_group.get(), "subgrid_mass_code", H5T_NATIVE_DOUBLE);
  state.black_holes.accretion_rate_code =
      readDataset1dAligned<double>(bh_group.get(), "accretion_rate_code", H5T_NATIVE_DOUBLE);
  state.black_holes.feedback_energy_code =
      readDataset1dAligned<double>(bh_group.get(), "feedback_energy_code", H5T_NATIVE_DOUBLE);
  state.black_holes.eddington_ratio = readDataset1dAligned<double>(bh_group.get(), "eddington_ratio", H5T_NATIVE_DOUBLE);
  state.black_holes.cumulative_accreted_mass_code =
      readDataset1dAligned<double>(bh_group.get(), "cumulative_accreted_mass_code", H5T_NATIVE_DOUBLE);
  state.black_holes.cumulative_feedback_energy_code =
      readDataset1dAligned<double>(bh_group.get(), "cumulative_feedback_energy_code", H5T_NATIVE_DOUBLE);
  state.black_holes.duty_cycle_active_time_code =
      readDataset1dAligned<double>(bh_group.get(), "duty_cycle_active_time_code", H5T_NATIVE_DOUBLE);
  state.black_holes.duty_cycle_total_time_code =
      readDataset1dAligned<double>(bh_group.get(), "duty_cycle_total_time_code", H5T_NATIVE_DOUBLE);

  Hdf5Handle tracer_group(H5Gopen2(state_group.get(), "tracers", H5P_DEFAULT));
  state.tracers.particle_index = readDataset1dAligned<std::uint32_t>(tracer_group.get(), "particle_index", H5T_NATIVE_UINT32);
  state.tracers.parent_particle_id =
      readDataset1dAligned<std::uint64_t>(tracer_group.get(), "parent_particle_id", H5T_NATIVE_UINT64);
  state.tracers.injection_step =
      readDataset1dAligned<std::uint64_t>(tracer_group.get(), "injection_step", H5T_NATIVE_UINT64);
  state.tracers.host_cell_index =
      readDataset1dAligned<std::uint32_t>(tracer_group.get(), "host_cell_index", H5T_NATIVE_UINT32);
  state.tracers.mass_fraction_of_host =
      readDataset1dAligned<double>(tracer_group.get(), "mass_fraction_of_host", H5T_NATIVE_DOUBLE);
  state.tracers.last_host_mass_code =
      readDataset1dAligned<double>(tracer_group.get(), "last_host_mass_code", H5T_NATIVE_DOUBLE);
  state.tracers.cumulative_exchanged_mass_code =
      readDataset1dAligned<double>(tracer_group.get(), "cumulative_exchanged_mass_code", H5T_NATIVE_DOUBLE);

  state.metadata = core::StateMetadata::deserialize(readStringDataset(state_group.get(), "state_metadata"));
  state.rebuildSpeciesIndex();

  Hdf5Handle sidecars_group(H5Gopen2(state_group.get(), "module_sidecars", H5P_DEFAULT));
  const std::string module_names_text = readStringDataset(state_group.get(), "module_sidecar_names");
  std::istringstream module_names_stream(module_names_text);
  std::string module_name;
  while (std::getline(module_names_stream, module_name)) {
    if (module_name.empty()) {
      continue;
    }
    Hdf5Handle module_group(H5Gopen2(sidecars_group.get(), module_name.c_str(), H5P_DEFAULT));
    core::ModuleSidecarBlock block;
    block.module_name = module_name;
    block.schema_version = readScalarU32Attribute(module_group.get(), "schema_version");
    if (hdf5AttributeExists(module_group.get(), "particle_indexed")) {
      block.particle_indexed = readScalarU32Attribute(module_group.get(), "particle_indexed") != 0U;
    }
    if (hdf5AttributeExists(module_group.get(), "row_stride_bytes")) {
      block.row_stride_bytes = readScalarU32Attribute(module_group.get(), "row_stride_bytes");
    }
    if (hdf5AttributeExists(module_group.get(), "required_species_mask")) {
      block.required_species_mask = readScalarU32Attribute(module_group.get(), "required_species_mask");
    }
    if (hdf5AttributeExists(module_group.get(), "requirement_kind")) {
      block.requirement.kind = static_cast<core::ModuleSidecarRequirementKind>(
          readScalarU32Attribute(module_group.get(), "requirement_kind"));
    }
    if (hdf5AttributeExists(module_group.get(), "requirement_species_mask")) {
      block.requirement.species_mask = readScalarU32Attribute(module_group.get(), "requirement_species_mask");
    }
    if (hdf5AttributeExists(module_group.get(), "requirement_particle_flags_mask")) {
      block.requirement.particle_flags_mask = readScalarU32Attribute(module_group.get(), "requirement_particle_flags_mask");
    }
    if (hdf5AttributeExists(module_group.get(), "requirement_threshold_code")) {
      block.requirement.threshold_code = readScalarF64Attribute(module_group.get(), "requirement_threshold_code");
    }
    block.payload = readDataset1d<std::byte>(
        module_group.get(), "payload", H5T_NATIVE_UINT8);
    if (block.particle_indexed || block.row_stride_bytes != 0U ||
        hdf5LinkExists(module_group.get(), "particle_id_by_row")) {
      block.particle_indexed = true;
      block.particle_id_by_row =
          readDataset1d<std::uint64_t>(module_group.get(), "particle_id_by_row", H5T_NATIVE_UINT64);
    }
    state.sidecars.upsert(std::move(block));
  }

  validateHydroGeometryStateForRestart(state, "restart reader");
  if (!state.validateOwnershipInvariants()) {
    throw std::runtime_error("restart state failed ownership invariant validation");
  }
}

// Preserve the existing contiguous uint8 dataset schema while keeping the
// serialization staging independent of the ownership-table population.
// HDF5's internal memory is governed by the existing external-runtime reserve.
class DistributedRestartHdf5Streambuf final : public std::streambuf {
 public:
  DistributedRestartHdf5Streambuf(hid_t dataset, std::size_t expected_bytes)
      : m_dataset(dataset), m_expected_bytes(expected_bytes) {}

  void finish() {
    flush();
    if (m_written != m_expected_bytes) {
      throw std::logic_error("distributed restart serialized length changed during write");
    }
  }

 protected:
  std::streamsize xsputn(const char* bytes, std::streamsize count) override {
    if (count < 0) throw std::length_error("negative distributed restart write length");
    const std::size_t size = static_cast<std::size_t>(count);
    if (size > m_expected_bytes - m_written - m_used) {
      throw std::length_error("distributed restart exceeds admitted dataset extent");
    }
    std::size_t consumed = 0U;
    while (consumed < size) {
      const std::size_t take = std::min(size - consumed, m_buffer.size() - m_used);
      std::memcpy(m_buffer.data() + m_used, bytes + consumed, take);
      m_used += take;
      consumed += take;
      if (m_used == m_buffer.size()) flush();
    }
    return count;
  }

  int_type overflow(int_type ch) override {
    if (!traits_type::eq_int_type(ch, traits_type::eof())) {
      const char value = traits_type::to_char_type(ch);
      (void)xsputn(&value, 1);
    }
    return traits_type::not_eof(ch);
  }

 private:
  void flush() {
    if (m_used == 0U) return;
    const hsize_t offset[1] = {static_cast<hsize_t>(m_written)};
    const hsize_t count[1] = {static_cast<hsize_t>(m_used)};
    Hdf5Handle file_space(H5Dget_space(m_dataset));
    Hdf5Handle memory_space(H5Screate_simple(1, count, nullptr));
    if (!file_space.valid() || !memory_space.valid() ||
        H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, offset, nullptr, count, nullptr) < 0 ||
        H5Dwrite(m_dataset, H5T_NATIVE_UINT8, memory_space.get(), file_space.get(),
                 H5P_DEFAULT, m_buffer.data()) < 0) {
      throw std::runtime_error("failed streaming distributed restart metadata");
    }
    m_written += m_used;
    m_used = 0U;
  }

  hid_t m_dataset;
  std::size_t m_expected_bytes;
  std::size_t m_written = 0U;
  std::size_t m_used = 0U;
  std::array<std::uint8_t, 65536U> m_buffer{};
};

void writeDistributedGravityGroup(hid_t root, const parallel::DistributedRestartState& distributed_state) {
  Hdf5Handle group(openOrCreateGroup(root, "/distributed_gravity"));
  const std::size_t byte_count = distributed_state.serializedSizeBytes();
  const hsize_t dims[1] = {core::checkedIntegralNarrow<hsize_t>(
      byte_count, "distributed restart dataset extent")};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group.get(), "state", H5T_STD_U8LE, space.get(),
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!dataset.valid()) throw std::runtime_error("failed creating distributed restart dataset");
  DistributedRestartHdf5Streambuf buffer(dataset.get(), byte_count);
  std::ostream stream(&buffer);
  stream.exceptions(std::ios::badbit | std::ios::failbit);
  distributed_state.serializeTo(stream);
  buffer.finish();
}

void readDistributedGravityGroup(hid_t root, parallel::DistributedRestartState& distributed_state) {
  Hdf5Handle group(H5Gopen2(root, "/distributed_gravity", H5P_DEFAULT));
  distributed_state = parallel::DistributedRestartState::deserialize(readStringDataset(group.get(), "state"));
}


void writeOutputCadenceGroup(hid_t root, const OutputCadencePersistentState& output_state) {
  Hdf5Handle group(openOrCreateGroup(root, "/output_cadence"));
  writeScalarU32Attribute(group.get(), "output_enabled", output_state.output_enabled ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "write_restarts", output_state.write_restarts ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "snapshot_due", output_state.snapshot_due ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "checkpoint_due", output_state.checkpoint_due ? 1U : 0U);
  writeScalarU64Attribute(group.get(), "last_completed_step_index", output_state.last_completed_step_index);
  writeScalarU64Attribute(group.get(), "snapshot_interval_steps", output_state.snapshot_interval_steps);
  writeScalarU64Attribute(group.get(), "next_snapshot_step_index", output_state.next_snapshot_step_index);
  writeScalarF64Attribute(
      group.get(), "snapshot_interval_time_code", output_state.snapshot_interval_time_code);
  writeScalarF64Attribute(group.get(), "next_snapshot_time_code", output_state.next_snapshot_time_code);
  writeScalarStringAttribute(group.get(), "snapshot_stem", output_state.snapshot_stem);
  writeScalarStringAttribute(group.get(), "restart_stem", output_state.restart_stem);
}

[[nodiscard]] OutputCadencePersistentState readOutputCadenceGroup(
    hid_t root,
    std::uint32_t schema_version) {
  Hdf5Handle group(H5Gopen2(root, "/output_cadence", H5P_DEFAULT));
  OutputCadencePersistentState output_state;
  output_state.output_enabled = readScalarU32Attribute(group.get(), "output_enabled") != 0U;
  output_state.write_restarts = readScalarU32Attribute(group.get(), "write_restarts") != 0U;
  output_state.snapshot_due = readScalarU32Attribute(group.get(), "snapshot_due") != 0U;
  output_state.checkpoint_due = readScalarU32Attribute(group.get(), "checkpoint_due") != 0U;
  output_state.last_completed_step_index = readScalarU64Attribute(group.get(), "last_completed_step_index");
  output_state.snapshot_interval_steps = readScalarU64Attribute(group.get(), "snapshot_interval_steps");
  output_state.next_snapshot_step_index = readScalarU64Attribute(group.get(), "next_snapshot_step_index");
  if (schema_version >= k_restart_schema_v21) {
    output_state.snapshot_interval_time_code =
        readScalarF64Attribute(group.get(), "snapshot_interval_time_code");
    output_state.next_snapshot_time_code =
        readScalarF64Attribute(group.get(), "next_snapshot_time_code");
  }
  output_state.snapshot_stem = readScalarStringAttribute(group.get(), "snapshot_stem");
  output_state.restart_stem = readScalarStringAttribute(group.get(), "restart_stem");
  return output_state;
}

void writeStochasticStateGroup(hid_t root, const StochasticPersistentState& stochastic_state) {
  Hdf5Handle group(openOrCreateGroup(root, "/stochastic_state"));
  const auto modules = sortedStochasticModules(stochastic_state);
  writeScalarU32Attribute(group.get(), "module_count", static_cast<std::uint32_t>(modules.size()));
  std::string module_names;
  for (const StochasticModulePersistentState& module_state : modules) {
    if (!module_names.empty()) {
      module_names.push_back('\n');
    }
    module_names += module_state.module_name;
    Hdf5Handle module_group(openOrCreateGroup(root, "/stochastic_state/" + module_state.module_name));
    writeScalarU32Attribute(module_group.get(), "schema_version", module_state.schema_version);
    writeScalarStringAttribute(module_group.get(), "rng_policy", module_state.rng_policy);
    writeScalarU64Attribute(module_group.get(), "random_seed", module_state.random_seed);
    writeScalarU32Attribute(module_group.get(), "rank_local_seed_offset", module_state.rank_local_seed_offset);
    writeScalarU64Attribute(module_group.get(), "last_committed_step_index", module_state.last_committed_step_index);
    writeScalarU32Attribute(
        module_group.get(),
        "deterministic_from_serialized_inputs",
        module_state.deterministic_from_serialized_inputs ? 1U : 0U);
  }
  writeStringDataset(group.get(), "module_names", module_names);
}

[[nodiscard]] StochasticPersistentState readStochasticStateGroup(hid_t root) {
  Hdf5Handle group(H5Gopen2(root, "/stochastic_state", H5P_DEFAULT));
  if (!group.valid()) {
    throw std::runtime_error("restart schema validation missing required group: /stochastic_state");
  }
  const std::string module_names_text = readStringDataset(group.get(), "module_names");
  StochasticPersistentState stochastic_state;
  std::istringstream module_names_stream(module_names_text);
  std::string module_name;
  while (std::getline(module_names_stream, module_name)) {
    if (module_name.empty()) {
      continue;
    }
    Hdf5Handle module_group(H5Gopen2(group.get(), module_name.c_str(), H5P_DEFAULT));
    if (!module_group.valid()) {
      throw std::runtime_error("restart schema validation missing required group: /stochastic_state/" + module_name);
    }
    StochasticModulePersistentState module_state;
    module_state.module_name = module_name;
    module_state.schema_version = readScalarU32Attribute(module_group.get(), "schema_version");
    module_state.rng_policy = readScalarStringAttribute(module_group.get(), "rng_policy");
    module_state.random_seed = readScalarU64Attribute(module_group.get(), "random_seed");
    module_state.rank_local_seed_offset = readScalarU32Attribute(module_group.get(), "rank_local_seed_offset");
    module_state.last_committed_step_index = readScalarU64Attribute(module_group.get(), "last_committed_step_index");
    module_state.deterministic_from_serialized_inputs =
        readScalarU32Attribute(module_group.get(), "deterministic_from_serialized_inputs") != 0U;
    stochastic_state.modules.push_back(std::move(module_state));
  }
  const std::uint32_t module_count = readScalarU32Attribute(group.get(), "module_count");
  if (module_count != stochastic_state.modules.size()) {
    throw std::runtime_error("restart schema validation invalid /stochastic_state/@module_count");
  }
  return stochastic_state;
}

[[nodiscard]] RestartDiagnosticsSummary makeRestartDiagnosticsSummary(
    const core::IntegratorState& integrator_state,
    const core::TimeBinPersistentState& scheduler_state,
    const core::TimeBinPersistentState& gas_cell_scheduler_state,
    const OutputCadencePersistentState& output_state,
    const StochasticPersistentState& stochastic_state) {
  RestartDiagnosticsSummary diagnostics;
  diagnostics.restart_schema_name = restartSchema().name;
  diagnostics.restart_schema_version = restartSchema().version;
  diagnostics.current_boundary_kind = std::string(core::stepBoundaryKindName(integrator_state.current_boundary_kind));
  diagnostics.last_completed_boundary_kind =
      std::string(core::stepBoundaryKindName(integrator_state.last_completed_boundary_kind));
  diagnostics.restart_safe = !integrator_state.inside_kdk_step && integrator_state.last_completed_restart_safe &&
      core::isRestartSafeBoundary(integrator_state.last_completed_boundary_kind);
  diagnostics.step_index = integrator_state.step_index;
  diagnostics.scheduler_current_tick = scheduler_state.current_tick;
  diagnostics.scheduler_max_bin = scheduler_state.max_bin;
  diagnostics.scheduler_element_count = scheduler_state.bin_index.size();
  diagnostics.scheduler_active_count = static_cast<std::uint64_t>(
      std::count(scheduler_state.active_flag.begin(), scheduler_state.active_flag.end(), static_cast<std::uint8_t>(1U)));
  diagnostics.scheduler_pending_transition_count = static_cast<std::uint64_t>(std::count_if(
      scheduler_state.pending_bin_index.begin(),
      scheduler_state.pending_bin_index.end(),
      [](std::uint8_t pending_bin) {
        return pending_bin != core::HierarchicalTimeBinScheduler::k_unset_pending_bin;
      }));
  diagnostics.gas_cell_scheduler_current_tick = gas_cell_scheduler_state.current_tick;
  diagnostics.gas_cell_scheduler_max_bin = gas_cell_scheduler_state.max_bin;
  diagnostics.gas_cell_scheduler_element_count = gas_cell_scheduler_state.bin_index.size();
  diagnostics.gas_cell_scheduler_active_count = static_cast<std::uint64_t>(
      std::count(gas_cell_scheduler_state.active_flag.begin(), gas_cell_scheduler_state.active_flag.end(), static_cast<std::uint8_t>(1U)));
  diagnostics.gas_cell_scheduler_pending_transition_count = static_cast<std::uint64_t>(std::count_if(
      gas_cell_scheduler_state.pending_bin_index.begin(),
      gas_cell_scheduler_state.pending_bin_index.end(),
      [](std::uint8_t pending_bin) {
        return pending_bin != core::HierarchicalTimeBinScheduler::k_unset_pending_bin;
      }));

  const core::PmSynchronizationPersistentState pm_sync_state =
      integrator_state.pm_sync_state.exportPersistentState();
  diagnostics.pm_cadence_steps = pm_sync_state.cadence_steps;
  diagnostics.pm_gravity_kick_opportunity = pm_sync_state.gravity_kick_opportunity;
  diagnostics.pm_field_version = pm_sync_state.field_version;
  diagnostics.pm_last_refresh_opportunity = pm_sync_state.last_refresh_opportunity;
  diagnostics.pm_last_refresh_step_index = pm_sync_state.last_refresh_step_index;
  diagnostics.pm_refresh_commit_pending = pm_sync_state.refresh_commit_pending;
  diagnostics.pm_long_range_field_valid = integrator_state.pm_long_range_field_valid;

  diagnostics.output_enabled = output_state.output_enabled;
  diagnostics.output_snapshot_due = output_state.snapshot_due;
  diagnostics.output_checkpoint_due = output_state.checkpoint_due;
  diagnostics.output_last_completed_step_index = output_state.last_completed_step_index;
  diagnostics.output_next_snapshot_step_index = output_state.next_snapshot_step_index;
  diagnostics.output_next_snapshot_time_code = output_state.next_snapshot_time_code;
  diagnostics.stochastic_module_count = stochastic_state.modules.size();
  return diagnostics;
}

void writeRestartDiagnosticsGroup(hid_t root, const RestartDiagnosticsSummary& diagnostics) {
  Hdf5Handle group(openOrCreateGroup(root, "/restart_diagnostics"));
  writeScalarStringAttribute(group.get(), "restart_schema_name", diagnostics.restart_schema_name);
  writeScalarU32Attribute(group.get(), "restart_schema_version", diagnostics.restart_schema_version);
  writeScalarStringAttribute(group.get(), "current_boundary_kind", diagnostics.current_boundary_kind);
  writeScalarStringAttribute(group.get(), "last_completed_boundary_kind", diagnostics.last_completed_boundary_kind);
  writeScalarU32Attribute(group.get(), "restart_safe", diagnostics.restart_safe ? 1U : 0U);
  writeScalarU64Attribute(group.get(), "step_index", diagnostics.step_index);
  writeScalarU64Attribute(group.get(), "scheduler_current_tick", diagnostics.scheduler_current_tick);
  writeScalarU32Attribute(group.get(), "scheduler_max_bin", diagnostics.scheduler_max_bin);
  writeScalarU64Attribute(group.get(), "scheduler_element_count", diagnostics.scheduler_element_count);
  writeScalarU64Attribute(group.get(), "scheduler_active_count", diagnostics.scheduler_active_count);
  writeScalarU64Attribute(
      group.get(), "scheduler_pending_transition_count", diagnostics.scheduler_pending_transition_count);
  writeScalarU64Attribute(
      group.get(), "gas_cell_scheduler_current_tick", diagnostics.gas_cell_scheduler_current_tick);
  writeScalarU32Attribute(
      group.get(), "gas_cell_scheduler_max_bin", diagnostics.gas_cell_scheduler_max_bin);
  writeScalarU64Attribute(
      group.get(), "gas_cell_scheduler_element_count", diagnostics.gas_cell_scheduler_element_count);
  writeScalarU64Attribute(
      group.get(), "gas_cell_scheduler_active_count", diagnostics.gas_cell_scheduler_active_count);
  writeScalarU64Attribute(
      group.get(),
      "gas_cell_scheduler_pending_transition_count",
      diagnostics.gas_cell_scheduler_pending_transition_count);
  writeScalarU64Attribute(group.get(), "pm_cadence_steps", diagnostics.pm_cadence_steps);
  writeScalarU64Attribute(group.get(), "pm_gravity_kick_opportunity", diagnostics.pm_gravity_kick_opportunity);
  writeScalarU64Attribute(group.get(), "pm_field_version", diagnostics.pm_field_version);
  writeScalarU64Attribute(group.get(), "pm_last_refresh_opportunity", diagnostics.pm_last_refresh_opportunity);
  writeScalarU64Attribute(group.get(), "pm_last_refresh_step_index", diagnostics.pm_last_refresh_step_index);
  writeScalarU32Attribute(group.get(), "pm_refresh_commit_pending", diagnostics.pm_refresh_commit_pending ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "pm_long_range_field_valid", diagnostics.pm_long_range_field_valid ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "output_enabled", diagnostics.output_enabled ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "output_snapshot_due", diagnostics.output_snapshot_due ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "output_checkpoint_due", diagnostics.output_checkpoint_due ? 1U : 0U);
  writeScalarU64Attribute(group.get(), "output_last_completed_step_index", diagnostics.output_last_completed_step_index);
  writeScalarU64Attribute(group.get(), "output_next_snapshot_step_index", diagnostics.output_next_snapshot_step_index);
  writeScalarF64Attribute(
      group.get(), "output_next_snapshot_time_code", diagnostics.output_next_snapshot_time_code);
  writeScalarU64Attribute(group.get(), "stochastic_module_count", diagnostics.stochastic_module_count);
}

[[nodiscard]] RestartDiagnosticsSummary readRestartDiagnosticsGroup(
    hid_t root,
    std::uint32_t schema_version) {
  Hdf5Handle group(H5Gopen2(root, "/restart_diagnostics", H5P_DEFAULT));
  if (!group.valid()) {
    throw std::runtime_error("restart schema validation missing required group: /restart_diagnostics");
  }
  RestartDiagnosticsSummary diagnostics;
  diagnostics.restart_schema_name = readScalarStringAttribute(group.get(), "restart_schema_name");
  diagnostics.restart_schema_version = readScalarU32Attribute(group.get(), "restart_schema_version");
  diagnostics.current_boundary_kind = readScalarStringAttribute(group.get(), "current_boundary_kind");
  diagnostics.last_completed_boundary_kind = readScalarStringAttribute(group.get(), "last_completed_boundary_kind");
  diagnostics.restart_safe = readScalarU32Attribute(group.get(), "restart_safe") != 0U;
  diagnostics.step_index = readScalarU64Attribute(group.get(), "step_index");
  diagnostics.scheduler_current_tick = readScalarU64Attribute(group.get(), "scheduler_current_tick");
  diagnostics.scheduler_max_bin = readScalarU32Attribute(group.get(), "scheduler_max_bin");
  diagnostics.scheduler_element_count = readScalarU64Attribute(group.get(), "scheduler_element_count");
  diagnostics.scheduler_active_count = readScalarU64Attribute(group.get(), "scheduler_active_count");
  diagnostics.scheduler_pending_transition_count =
      readScalarU64Attribute(group.get(), "scheduler_pending_transition_count");
  if (schema_version >= k_restart_schema_v19) {
    diagnostics.gas_cell_scheduler_current_tick =
        readScalarU64Attribute(group.get(), "gas_cell_scheduler_current_tick");
    diagnostics.gas_cell_scheduler_max_bin =
        readScalarU32Attribute(group.get(), "gas_cell_scheduler_max_bin");
    diagnostics.gas_cell_scheduler_element_count =
        readScalarU64Attribute(group.get(), "gas_cell_scheduler_element_count");
    diagnostics.gas_cell_scheduler_active_count =
        readScalarU64Attribute(group.get(), "gas_cell_scheduler_active_count");
    diagnostics.gas_cell_scheduler_pending_transition_count =
        readScalarU64Attribute(group.get(), "gas_cell_scheduler_pending_transition_count");
  }
  diagnostics.pm_cadence_steps = readScalarU64Attribute(group.get(), "pm_cadence_steps");
  diagnostics.pm_gravity_kick_opportunity = readScalarU64Attribute(group.get(), "pm_gravity_kick_opportunity");
  diagnostics.pm_field_version = readScalarU64Attribute(group.get(), "pm_field_version");
  diagnostics.pm_last_refresh_opportunity = readScalarU64Attribute(group.get(), "pm_last_refresh_opportunity");
  diagnostics.pm_last_refresh_step_index = readScalarU64Attribute(group.get(), "pm_last_refresh_step_index");
  diagnostics.pm_refresh_commit_pending = readScalarU32Attribute(group.get(), "pm_refresh_commit_pending") != 0U;
  diagnostics.pm_long_range_field_valid = readScalarU32Attribute(group.get(), "pm_long_range_field_valid") != 0U;
  diagnostics.output_enabled = readScalarU32Attribute(group.get(), "output_enabled") != 0U;
  diagnostics.output_snapshot_due = readScalarU32Attribute(group.get(), "output_snapshot_due") != 0U;
  diagnostics.output_checkpoint_due = readScalarU32Attribute(group.get(), "output_checkpoint_due") != 0U;
  diagnostics.output_last_completed_step_index = readScalarU64Attribute(group.get(), "output_last_completed_step_index");
  diagnostics.output_next_snapshot_step_index = readScalarU64Attribute(group.get(), "output_next_snapshot_step_index");
  if (schema_version >= k_restart_schema_v21) {
    diagnostics.output_next_snapshot_time_code =
        readScalarF64Attribute(group.get(), "output_next_snapshot_time_code");
  }
  diagnostics.stochastic_module_count = readScalarU64Attribute(group.get(), "stochastic_module_count");
  return diagnostics;
}
#endif

}  // namespace

const RestartSchema& restartSchema() {
  static const RestartSchema schema{};
  return schema;
}

GravityForceCachePersistentView makeGravityForceCachePersistentView(
    const GravityForceCachePersistentState& cache) noexcept {
  return GravityForceCachePersistentView{
      .valid = cache.valid,
      .particle_id = cache.particle_id,
      .gas_cell_id = cache.gas_cell_id,
      .particle_accel_x_comoving = cache.particle_accel_x_comoving,
      .particle_accel_y_comoving = cache.particle_accel_y_comoving,
      .particle_accel_z_comoving = cache.particle_accel_z_comoving,
      .cell_accel_x_comoving = cache.cell_accel_x_comoving,
      .cell_accel_y_comoving = cache.cell_accel_y_comoving,
      .cell_accel_z_comoving = cache.cell_accel_z_comoving,
      .hierarchical_endpoint_valid = cache.hierarchical_endpoint_valid,
      .hierarchical_pm_field_version = cache.hierarchical_pm_field_version,
      .hierarchical_tree_accel_x_comoving = cache.hierarchical_tree_accel_x_comoving,
      .hierarchical_tree_accel_y_comoving = cache.hierarchical_tree_accel_y_comoving,
      .hierarchical_tree_accel_z_comoving = cache.hierarchical_tree_accel_z_comoving,
      .hierarchical_pm_accel_x_comoving = cache.hierarchical_pm_accel_x_comoving,
      .hierarchical_pm_accel_y_comoving = cache.hierarchical_pm_accel_y_comoving,
      .hierarchical_pm_accel_z_comoving = cache.hierarchical_pm_accel_z_comoving,
      .hierarchical_nearest_pair_r2_comoving = cache.hierarchical_nearest_pair_r2_comoving};
}

bool isRestartSchemaCompatible(std::uint32_t file_schema_version) {
  return file_schema_version == restartSchema().version ||
      file_schema_version == k_restart_schema_v23 ||
      file_schema_version == k_restart_schema_v22 ||
      file_schema_version == k_restart_schema_v21 ||
      file_schema_version == k_restart_schema_v20 ||
      file_schema_version == k_restart_schema_v19 ||
      file_schema_version == k_restart_schema_v18 ||
      file_schema_version == k_restart_schema_v17 || file_schema_version == k_restart_schema_v16 ||
      file_schema_version == k_restart_schema_v15 || file_schema_version == k_restart_schema_v14;
}

const std::vector<std::string_view>& exactRestartCompletenessChecklist() {
  static const std::vector<std::string_view> checklist = {
      "simulation_state_lanes_and_metadata",
      "particle_identity_softening_and_drift_epoch_lanes",
      "gas_cell_identity_lanes",
      "gas_metal_mass_and_star_birth_identity_lanes",
      "hydro_geometry_patch_state",
      "amr_pending_flux_register_state",
      "amr_temporal_boundary_history_state",
      "species_specific_sidecars",
      "module_sidecars_with_schema_versions",
      "integrator_state",
      "integrator_owned_pm_sync_state",
      "gravity_force_cache_at_kdk_boundary",
      "scheduler_persistent_state",
      "gas_cell_scheduler_persistent_state_keyed_by_gas_cell_id",
      "output_cadence_persistent_state",
      "stochastic_module_persistent_state",
      "restart_diagnostics_summary",
      "distributed_gravity_state",
      "normalized_config_text_and_hash",
      "provenance_record",
      "payload_integrity_hash_and_hex"};
  return checklist;
}

RestartIntegrityDigests restartPayloadIntegrityDigestsImpl(
    const RestartWritePayload& payload,
    bool include_gas_identity_records,
    bool include_pending_flux_registers,
    bool include_temporal_boundary_history,
    bool include_gas_cell_scheduler,
    bool include_gravity_force_cache,
    bool include_output_time_cadence,
    bool include_star_formation_v22,
    bool include_hierarchical_endpoint_v24,
    core::MemoryGovernor* memory_governor = nullptr) {
  if (payload.persistent_state.simulation_state == nullptr || payload.integrator_state == nullptr || payload.scheduler == nullptr) {
    throw std::invalid_argument("restart payload must provide state, integrator_state, and scheduler");
  }
  (void)memory_governor;  // Integrity traversal is allocation-free in the production payload path.
  validateContinuationMetadata(
      payload.normalized_config_text,
      payload.normalized_config_hash_hex,
      payload.provenance,
      "restart payload");
  validateHydroGeometryStateForRestart(*payload.persistent_state.simulation_state, "restart payload");
  if (!payload.persistent_state.simulation_state->validateOwnershipInvariants()) {
    throw std::invalid_argument("restart payload state failed ownership invariant validation");
  }
  const GravityForceCachePersistentView empty_force_cache{};
  const GravityForceCachePersistentView owning_force_cache_view =
      payload.gravity_force_cache != nullptr
          ? makeGravityForceCachePersistentView(*payload.gravity_force_cache)
          : empty_force_cache;
  const GravityForceCachePersistentView& gravity_force_cache =
      payload.gravity_force_cache_view != nullptr
          ? *payload.gravity_force_cache_view
          : owning_force_cache_view;
  if (include_gravity_force_cache) {
    validateGravityForceCacheForRestart(
        gravity_force_cache, *payload.persistent_state.simulation_state, "restart payload");
  }
  const bool uniform_particle_scheduler =
      payload.scheduler->representation() == core::SchedulerRepresentation::kUniformRungZero;
  const core::SimulationState& validation_state = *payload.persistent_state.simulation_state;
  core::assertCanWriteCheckpointAtBoundary(*payload.integrator_state, payload.scheduler->currentTick());
  if (payload.scheduler->elementCount() != validation_state.particles.size()) {
    throw std::invalid_argument("restart payload particle scheduler extent does not match particle state");
  }
  if (uniform_particle_scheduler && payload.scheduler->maxBin() != 0U) {
    throw std::invalid_argument("restart payload uniform particle scheduler rung contract is invalid");
  }
  for (std::size_t i = 0; i < validation_state.particles.size(); ++i) {
    if (validation_state.particleTimeBin(i) !=
        payload.scheduler->binIndex(static_cast<std::uint32_t>(i))) {
      throw std::invalid_argument(
          "restart payload particle time_bin mirror is stale relative to scheduler authority");
    }
  }
  if (payload.gas_cell_scheduler != nullptr) {
    if (payload.gas_cell_scheduler->currentTick() != payload.scheduler->currentTick()) {
      throw std::invalid_argument("restart payload particle and gas-cell schedulers must share current_tick");
    }
    if (payload.gas_cell_scheduler->elementCount() != validation_state.cells.size()) {
      throw std::invalid_argument("restart payload gas-cell scheduler extent does not match cell state");
    }
    validation_state.requireGasCellIdentityMapCoversDenseRows("restart payload");
    for (std::uint32_t cell_index = 0; cell_index < validation_state.cells.size(); ++cell_index) {
      if (validation_state.gas_cell_identity.findByLocalRow(cell_index) == nullptr) {
        throw std::invalid_argument("restart payload gas-cell identity map is missing a dense row");
      }
      if (validation_state.cells.time_bin[cell_index] !=
          payload.gas_cell_scheduler->binIndex(cell_index)) {
        throw std::invalid_argument(
            "restart payload cell time_bin mirror is stale relative to gas-cell scheduler authority");
      }
    }
  } else if (validation_state.cells.size() != 0U) {
    throw std::invalid_argument(
        "restart payload with gas cells must provide the authoritative gas-cell scheduler");
  }
  if (payload.distributed_gravity_state.world_size <= 0) {
    throw std::invalid_argument("restart payload distributed_gravity_state.world_size must be positive");
  }
  if (payload.distributed_gravity_state.owningRankItemCount() !=
      payload.persistent_state.simulation_state->particles.size()) {
    throw std::invalid_argument(
        "restart payload distributed_gravity_state logical ownership count must match particle count");
  }
  if (payload.distributed_gravity_state.pm_slab_begin_x_by_rank.size() !=
          static_cast<std::size_t>(payload.distributed_gravity_state.world_size) ||
      payload.distributed_gravity_state.pm_slab_end_x_by_rank.size() !=
          static_cast<std::size_t>(payload.distributed_gravity_state.world_size)) {
    throw std::invalid_argument(
        "restart payload distributed_gravity_state PM slab ownership must match world_size");
  }
  if (payload.distributed_gravity_state.long_range_restart_policy != "deterministic_rebuild") {
    throw std::invalid_argument(
        "restart payload distributed_gravity_state.long_range_restart_policy must be deterministic_rebuild");
  }
  validateOutputCadenceStateForRestart(
      payload.output_cadence_state, *payload.integrator_state, "restart payload");
  validateStochasticStateForRestart(
      payload.stochastic_state, *payload.integrator_state, "restart payload");

  std::uint64_t hash = k_offset_basis;
  core::internal::Sha256 strong_hash;
  strong_hash.update("cosmosim-restart-payload-v23");

  const auto append_u64 = [&hash, &strong_hash](std::uint64_t value) {
    const std::array<std::byte, sizeof(std::uint64_t)> bytes =
        std::bit_cast<std::array<std::byte, sizeof(std::uint64_t)>>(value);
    hash = fnv1aAppend(hash, bytes);
    appendCanonicalSha256(strong_hash, value);
  };

  const auto append_string = [&hash, &strong_hash, &append_u64](const std::string& value) {
    append_u64(static_cast<std::uint64_t>(value.size()));
    hash = fnv1aAppend(hash, {reinterpret_cast<const std::byte*>(value.data()), value.size()});
    strong_hash.update(reinterpret_cast<const std::uint8_t*>(value.data()), value.size());
  };

  append_string(payload.normalized_config_text);
  append_string(payload.normalized_config_hash_hex);
  append_string(core::serializeProvenanceRecord(payload.provenance));
  append_string(payload.persistent_state.simulation_state->metadata.serialize());

  const auto append_any_vec = [&hash, &strong_hash, &append_u64](const auto& values) {
    const auto bytes = asBytesSpan(values);
    append_u64(static_cast<std::uint64_t>(values.size()));
    append_u64(static_cast<std::uint64_t>(bytes.size()));
    hash = fnv1aAppend(hash, bytes);
    for (const auto value : values) {
      appendCanonicalSha256(strong_hash, value);
    }
  };


  const auto append_logical_lane = [&hash, &strong_hash, &append_u64](std::size_t count, auto getter) {
    using Value = std::decay_t<decltype(getter(std::size_t{}))>;
    append_u64(static_cast<std::uint64_t>(count));
    append_u64(static_cast<std::uint64_t>(count) * sizeof(Value));
    for (std::size_t i = 0; i < count; ++i) {
      const Value value = getter(i);
      hash = fnv1aAppend(hash, {reinterpret_cast<const std::byte*>(&value), sizeof(Value)});
      appendCanonicalSha256(strong_hash, value);
    }
  };

  const core::SimulationState& state = *payload.persistent_state.simulation_state;
  append_any_vec(state.particles.position_x_comoving);
  append_any_vec(state.particles.position_y_comoving);
  append_any_vec(state.particles.position_z_comoving);
  append_any_vec(state.particles.velocity_x_peculiar);
  append_any_vec(state.particles.velocity_y_peculiar);
  append_any_vec(state.particles.velocity_z_peculiar);
  append_any_vec(state.particles.mass_code);
  if (state.hasHomogeneousDmoMetadata()) {
    append_logical_lane(state.particles.size(), [&](std::size_t i) { return state.particleTimeBin(i); });
  } else {
    append_any_vec(state.particles.time_bin);
  }
  append_any_vec(state.particle_sidecar.particle_id);
  if (state.hasHomogeneousDmoMetadata()) {
    append_logical_lane(state.particles.size(), [&](std::size_t i) { return state.particleSfcKey(i); });
    append_logical_lane(state.particles.size(), [&](std::size_t i) { return state.particleSpeciesTag(i); });
    append_logical_lane(state.particles.size(), [&](std::size_t i) { return state.particleFlags(i); });
    append_logical_lane(state.particles.size(), [&](std::size_t i) { return state.particleOwningRank(i); });
    append_logical_lane(state.particles.size(), [&](std::size_t i) { return state.particleLastDriftTimeCode(i); });
    append_logical_lane(state.particles.size(), [&](std::size_t i) { return state.particleLastDriftScaleFactor(i); });
  } else {
    append_any_vec(state.particle_sidecar.sfc_key);
    append_any_vec(state.particle_sidecar.species_tag);
    append_any_vec(state.particle_sidecar.particle_flags);
    append_any_vec(state.particle_sidecar.owning_rank);
    append_any_vec(state.particle_sidecar.last_drift_time_code);
    append_any_vec(state.particle_sidecar.last_drift_scale_factor);
  }
  append_any_vec(state.particle_sidecar.gravity_softening_comoving);
  append_any_vec(state.particle_sidecar.has_gravity_softening_override);

  append_any_vec(state.cells.center_x_comoving);
  append_any_vec(state.cells.center_y_comoving);
  append_any_vec(state.cells.center_z_comoving);
  append_any_vec(state.cells.mass_code);
  append_any_vec(state.cells.time_bin);
  append_any_vec(state.cells.patch_index);

  append_any_vec(state.gas_cells.gas_cell_id);
  append_any_vec(state.gas_cells.parent_particle_id);
  append_any_vec(state.gas_cells.velocity_x_peculiar);
  append_any_vec(state.gas_cells.velocity_y_peculiar);
  append_any_vec(state.gas_cells.velocity_z_peculiar);
  append_any_vec(state.gas_cells.density_code);
  append_any_vec(state.gas_cells.pressure_code);
  append_any_vec(state.gas_cells.internal_energy_code);
  if (include_star_formation_v22) {
    append_any_vec(state.gas_cells.metal_mass_code);
  }
  append_any_vec(state.gas_cells.temperature_code);
  append_any_vec(state.gas_cells.sound_speed_code);
  if (include_gas_identity_records) {
    append_string(std::string(k_gas_identity_row_policy));
    state.requireGasCellIdentityMapCoversDenseRows("restart integrity");
    const std::size_t identity_count = state.cells.size();
    const auto identity_at = [&](std::size_t i) -> const core::GasCellIdentityRecord& {
      const auto* record = state.gas_cell_identity.findByLocalRow(static_cast<std::uint32_t>(i));
      if (record == nullptr) throw std::logic_error("restart integrity gas identity row missing");
      return *record;
    };
    append_logical_lane(identity_count, [&](std::size_t i) { return identity_at(i).gas_cell_id; });
    append_logical_lane(identity_count, [&](std::size_t i) {
      return static_cast<std::uint8_t>(identity_at(i).parent_particle_id.has_value() ? 1U : 0U);
    });
    append_logical_lane(identity_count, [&](std::size_t i) {
      return identity_at(i).parent_particle_id.value_or(0U);
    });
    append_logical_lane(identity_count, [&](std::size_t i) { return identity_at(i).owning_patch_id; });
    append_logical_lane(identity_count, [&](std::size_t i) { return identity_at(i).local_cell_row; });
  }

  append_any_vec(state.patches.patch_id);
  append_any_vec(state.patches.level);
  append_any_vec(state.patches.first_cell);
  append_any_vec(state.patches.cell_count);
  append_any_vec(state.patches.parent_patch_id);
  append_any_vec(state.patches.morton_key);
  append_any_vec(state.patches.origin_x_comoving);
  append_any_vec(state.patches.origin_y_comoving);
  append_any_vec(state.patches.origin_z_comoving);
  append_any_vec(state.patches.extent_x_comoving);
  append_any_vec(state.patches.extent_y_comoving);
  append_any_vec(state.patches.extent_z_comoving);
  append_any_vec(state.patches.cell_dim_x);
  append_any_vec(state.patches.cell_dim_y);
  append_any_vec(state.patches.cell_dim_z);
  append_any_vec(state.patches.owning_rank);
  if (include_pending_flux_registers) {
    const auto pending_flux_records = state.pending_flux_registers.records();
    append_u64(static_cast<std::uint64_t>(pending_flux_records.size()));
    std::uint64_t previous_key = 0U;
    bool have_previous = false;
    for (std::size_t ordered = 0; ordered < pending_flux_records.size(); ++ordered) {
      const core::PendingFluxRegisterRecord* selected = nullptr;
      for (const core::PendingFluxRegisterRecord& candidate : pending_flux_records) {
        if (have_previous && candidate.register_key <= previous_key) continue;
        if (selected == nullptr || candidate.register_key < selected->register_key) selected = &candidate;
      }
      if (selected == nullptr) {
        throw std::invalid_argument("restart integrity pending-flux register keys must be unique");
      }
      const core::PendingFluxRegisterRecord& record = *selected;
      previous_key = record.register_key;
      have_previous = true;
      append_u64(record.register_key);
      append_u64(record.coarse_patch_id);
      append_u64(record.coarse_gas_cell_id);
      append_u64(static_cast<std::uint64_t>(record.coarse_cell_index));
      append_u64(record.level);
      append_u64(record.axis);
      append_u64(record.orientation);
      append_u64(std::bit_cast<std::uint64_t>(record.expected_area_comov));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_area_accumulated_comov));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_area_accumulated_comov));
      append_u64(std::bit_cast<std::uint64_t>(record.interval_start_code));
      append_u64(std::bit_cast<std::uint64_t>(record.interval_end_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_dt_code));
      append_u64(record.expected_fine_substeps);
      append_u64(record.completed_fine_substeps);
      append_u64(record.fine_substep_coverage_mask);
      append_u64(record.coarse_face_count);
      append_u64(record.fine_face_count);
      append_u64(record.gas_cell_identity_generation);
      append_u64(record.patch_geometry_generation);
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_mass_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_momentum_x_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_momentum_y_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_momentum_z_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_total_energy_flux_integral_code));
      if (include_star_formation_v22) append_u64(std::bit_cast<std::uint64_t>(record.coarse_metal_mass_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_mass_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_momentum_x_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_momentum_y_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_momentum_z_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_total_energy_flux_integral_code));
      if (include_star_formation_v22) append_u64(std::bit_cast<std::uint64_t>(record.fine_metal_mass_flux_integral_code));
    }
  }

  if (include_temporal_boundary_history) {
    const auto temporal_source = state.amr_temporal_boundary_history.records();
    append_u64(static_cast<std::uint64_t>(temporal_source.size()));
    std::uint64_t previous_patch_id = 0U;
    bool have_previous_patch = false;
    for (std::size_t ordered = 0; ordered < temporal_source.size(); ++ordered) {
      const core::AmrTemporalBoundaryHistoryRecord* selected = nullptr;
      for (const auto& candidate : temporal_source) {
        if (have_previous_patch && candidate.patch_id <= previous_patch_id) continue;
        if (selected == nullptr || candidate.patch_id < selected->patch_id) selected = &candidate;
      }
      if (selected == nullptr) {
        throw std::invalid_argument("restart integrity temporal-history patch IDs must be unique");
      }
      const auto& record = *selected;
      previous_patch_id = record.patch_id;
      have_previous_patch = true;
      append_u64(record.patch_id);
      append_u64(record.patch_level);
      append_u64(record.patch_geometry_fingerprint);
      append_u64(record.gas_cell_identity_generation);
      append_u64(std::bit_cast<std::uint64_t>(record.interval_start_code));
      append_u64(std::bit_cast<std::uint64_t>(record.interval_end_code));
      append_u64(record.end_state_valid ? 1U : 0U);
      append_u64(static_cast<std::uint64_t>(record.cells.size()));
      std::uint64_t previous_local = 0U;
      std::uint64_t previous_id = 0U;
      bool have_previous_cell = false;
      for (std::size_t cell_ordered = 0; cell_ordered < record.cells.size(); ++cell_ordered) {
        const core::AmrTemporalBoundaryHistoryCellRecord* selected_cell = nullptr;
        for (const auto& candidate : record.cells) {
          const std::uint64_t local = static_cast<std::uint64_t>(candidate.patch_local_cell);
          const bool after_previous = !have_previous_cell || local > previous_local ||
              (local == previous_local && candidate.gas_cell_id > previous_id);
          if (!after_previous) continue;
          if (selected_cell == nullptr || local < static_cast<std::uint64_t>(selected_cell->patch_local_cell) ||
              (local == static_cast<std::uint64_t>(selected_cell->patch_local_cell) &&
               candidate.gas_cell_id < selected_cell->gas_cell_id)) {
            selected_cell = &candidate;
          }
        }
        if (selected_cell == nullptr) {
          throw std::invalid_argument("restart integrity temporal-history cell keys must be unique");
        }
        const auto& cell = *selected_cell;
        previous_local = static_cast<std::uint64_t>(cell.patch_local_cell);
        previous_id = cell.gas_cell_id;
        have_previous_cell = true;
        append_u64(cell.gas_cell_id);
        append_u64(static_cast<std::uint64_t>(cell.patch_local_cell));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_mass_density_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_momentum_density_x_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_momentum_density_y_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_momentum_density_z_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_total_energy_density_comoving));
        if (include_star_formation_v22) append_u64(std::bit_cast<std::uint64_t>(cell.start_metal_mass_density_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_mass_density_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_momentum_density_x_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_momentum_density_y_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_momentum_density_z_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_total_energy_density_comoving));
        if (include_star_formation_v22) append_u64(std::bit_cast<std::uint64_t>(cell.end_metal_mass_density_comoving));
      }
    }
  }

  append_any_vec(state.star_particles.particle_index);
  append_any_vec(state.star_particles.formation_scale_factor);
  append_any_vec(state.star_particles.birth_mass_code);
  append_any_vec(state.star_particles.metallicity_mass_fraction);
  if (include_star_formation_v22) {
    append_any_vec(state.star_particles.birth_key);
    append_any_vec(state.star_particles.parent_gas_cell_id);
    append_any_vec(state.star_particles.birth_tick);
    append_any_vec(state.star_particles.birth_ordinal);
  }
  append_any_vec(state.star_particles.stellar_age_years_last);
  append_any_vec(state.star_particles.stellar_returned_mass_cumulative_code);
  append_any_vec(state.star_particles.stellar_returned_metals_cumulative_code);
  append_any_vec(state.star_particles.stellar_newly_synthesized_metals_cumulative_code);
  append_any_vec(state.star_particles.stellar_feedback_energy_cumulative_erg);
  append_any_vec(state.star_particles.enrichment_carry_mass_code);
  append_any_vec(state.star_particles.enrichment_carry_metals_code);
  append_any_vec(state.star_particles.enrichment_carry_feedback_energy_erg);
  append_any_vec(state.star_particles.enrichment_carry_momentum_code);
  append_any_vec(state.star_particles.stellar_deposited_mass_cumulative_code);
  append_any_vec(state.star_particles.stellar_deposited_metals_cumulative_code);
  append_any_vec(state.star_particles.stellar_deposited_feedback_energy_cumulative_erg);
  for (std::size_t channel = 0; channel < state.star_particles.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    append_any_vec(state.star_particles.stellar_returned_mass_channel_cumulative_code[channel]);
    append_any_vec(state.star_particles.stellar_returned_metals_channel_cumulative_code[channel]);
    append_any_vec(state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel]);
  }
  append_any_vec(state.black_holes.particle_index);
  append_any_vec(state.black_holes.host_cell_index);
  append_any_vec(state.black_holes.subgrid_mass_code);
  append_any_vec(state.black_holes.accretion_rate_code);
  append_any_vec(state.black_holes.feedback_energy_code);
  append_any_vec(state.black_holes.eddington_ratio);
  append_any_vec(state.black_holes.cumulative_accreted_mass_code);
  append_any_vec(state.black_holes.cumulative_feedback_energy_code);
  append_any_vec(state.black_holes.duty_cycle_active_time_code);
  append_any_vec(state.black_holes.duty_cycle_total_time_code);
  append_any_vec(state.tracers.particle_index);
  append_any_vec(state.tracers.parent_particle_id);
  append_any_vec(state.tracers.injection_step);
  append_any_vec(state.tracers.host_cell_index);
  append_any_vec(state.tracers.mass_fraction_of_host);
  append_any_vec(state.tracers.last_host_mass_code);
  append_any_vec(state.tracers.cumulative_exchanged_mass_code);
  append_any_vec(state.species.count_by_species);

  const auto ordered_sidecars = state.sidecars.blocksSortedByName();
  for (const core::ModuleSidecarBlock* block : ordered_sidecars) {
    append_string(block->module_name);
    append_u64(static_cast<std::uint64_t>(block->schema_version));
    append_u64(block->particle_indexed ? 1ULL : 0ULL);
    append_u64(static_cast<std::uint64_t>(block->row_stride_bytes));
    append_u64(static_cast<std::uint64_t>(block->required_species_mask));
    append_u64(static_cast<std::uint64_t>(block->requirement.kind));
    append_u64(static_cast<std::uint64_t>(block->requirement.species_mask));
    append_u64(static_cast<std::uint64_t>(block->requirement.particle_flags_mask));
    append_u64(std::bit_cast<std::uint64_t>(block->requirement.threshold_code));
    append_u64(static_cast<std::uint64_t>(block->particle_id_by_row.size()));
    for (const std::uint64_t particle_id : block->particle_id_by_row) {
      append_u64(particle_id);
    }
    append_u64(static_cast<std::uint64_t>(block->payload.size()));
    hash = fnv1aAppend(hash, std::span<const std::byte>(block->payload.data(), block->payload.size()));
    strong_hash.update(reinterpret_cast<const std::uint8_t*>(block->payload.data()), block->payload.size());
  }

  append_u64(payload.integrator_state->step_index);
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_time_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_scale_factor));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_redshift));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_hubble_rate_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->time_si_per_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->dt_time_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_drift_factor_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_first_kick_factor_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_second_kick_factor_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_first_hubble_drag_factor));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_second_hubble_drag_factor));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->scheme));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->current_boundary_kind));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->last_completed_boundary_kind));
  append_u64(payload.integrator_state->inside_kdk_step ? 1ull : 0ull);
  append_u64(payload.integrator_state->last_completed_restart_safe ? 1ull : 0ull);
  append_u64(payload.integrator_state->time_bins.hierarchical_enabled ? 1ull : 0ull);
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->time_bins.active_bin));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->time_bins.max_bin));
  append_u64(payload.integrator_state->pm_long_range_field_valid ? 1ull : 0ull);
  if (include_gravity_force_cache) {
    append_u64(payload.integrator_state->pm_refresh_enabled ? 1ull : 0ull);
  }
  const core::PmSynchronizationPersistentState pm_sync_state = payload.integrator_state->pm_sync_state.exportPersistentState();
  append_u64(pm_sync_state.cadence_steps);
  append_u64(pm_sync_state.gravity_kick_opportunity);
  append_u64(pm_sync_state.last_refresh_opportunity);
  append_u64(pm_sync_state.field_version);
  append_u64(pm_sync_state.last_refresh_step_index);
  append_u64(std::bit_cast<std::uint64_t>(pm_sync_state.last_refresh_scale_factor));
  append_u64(pm_sync_state.refresh_commit_pending ? 1ull : 0ull);
  append_u64(pm_sync_state.pending_refresh_opportunity);
  append_u64(pm_sync_state.pending_refresh_field_version);
  if (include_gravity_force_cache) {
    append_u64(gravity_force_cache.valid ? 1ULL : 0ULL);
    append_any_vec(gravity_force_cache.particle_id);
    append_any_vec(gravity_force_cache.gas_cell_id);
    append_any_vec(gravity_force_cache.particle_accel_x_comoving);
    append_any_vec(gravity_force_cache.particle_accel_y_comoving);
    append_any_vec(gravity_force_cache.particle_accel_z_comoving);
    append_any_vec(gravity_force_cache.cell_accel_x_comoving);
    append_any_vec(gravity_force_cache.cell_accel_y_comoving);
    append_any_vec(gravity_force_cache.cell_accel_z_comoving);
  }
  if (include_hierarchical_endpoint_v24) {
    append_u64(gravity_force_cache.hierarchical_endpoint_valid ? 1ULL : 0ULL);
    append_u64(gravity_force_cache.hierarchical_pm_field_version);
    append_any_vec(gravity_force_cache.hierarchical_tree_accel_x_comoving);
    append_any_vec(gravity_force_cache.hierarchical_tree_accel_y_comoving);
    append_any_vec(gravity_force_cache.hierarchical_tree_accel_z_comoving);
    append_any_vec(gravity_force_cache.hierarchical_pm_accel_x_comoving);
    append_any_vec(gravity_force_cache.hierarchical_pm_accel_y_comoving);
    append_any_vec(gravity_force_cache.hierarchical_pm_accel_z_comoving);
    append_any_vec(gravity_force_cache.hierarchical_nearest_pair_r2_comoving);
  }

  append_u64(payload.scheduler->currentTick());
  append_u64(static_cast<std::uint64_t>(payload.scheduler->maxBin()));
  const std::size_t scheduler_count = payload.scheduler->elementCount();
  append_logical_lane(scheduler_count, [&](std::size_t i) {
    return payload.scheduler->binIndex(static_cast<std::uint32_t>(i));
  });
  append_logical_lane(scheduler_count, [&](std::size_t i) {
    return payload.scheduler->nextActivationTick(static_cast<std::uint32_t>(i));
  });
  append_logical_lane(scheduler_count, [&](std::size_t i) {
    return static_cast<std::uint8_t>(
        payload.scheduler->isElementActive(static_cast<std::uint32_t>(i)) ? 1U : 0U);
  });
  append_logical_lane(scheduler_count, [&](std::size_t i) {
    return payload.scheduler->pendingBinIndex(static_cast<std::uint32_t>(i));
  });
  if (include_gas_cell_scheduler) {
    append_string(std::string(k_gas_cell_scheduler_identity_key));
    const std::size_t gas_scheduler_count = validation_state.cells.size();
    const std::uint64_t gas_current_tick =
        payload.gas_cell_scheduler != nullptr ? payload.gas_cell_scheduler->currentTick()
                                              : payload.scheduler->currentTick();
    const std::uint8_t gas_max_bin =
        payload.gas_cell_scheduler != nullptr ? payload.gas_cell_scheduler->maxBin()
                                              : payload.scheduler->maxBin();
    append_u64(gas_current_tick);
    append_u64(static_cast<std::uint64_t>(gas_max_bin));
    append_logical_lane(gas_scheduler_count, [&](std::size_t i) {
      const auto* identity = validation_state.gas_cell_identity.findByLocalRow(
          static_cast<std::uint32_t>(i));
      if (identity == nullptr) throw std::logic_error("gas scheduler integrity identity row missing");
      return identity->gas_cell_id;
    });
    append_logical_lane(gas_scheduler_count, [&](std::size_t i) {
      return payload.gas_cell_scheduler != nullptr
          ? payload.gas_cell_scheduler->binIndex(static_cast<std::uint32_t>(i))
          : validation_state.cells.time_bin[i];
    });
    append_logical_lane(gas_scheduler_count, [&](std::size_t i) {
      return payload.gas_cell_scheduler != nullptr
          ? payload.gas_cell_scheduler->nextActivationTick(static_cast<std::uint32_t>(i))
          : gas_current_tick;
    });
    append_logical_lane(gas_scheduler_count, [&](std::size_t i) {
      return static_cast<std::uint8_t>(
          payload.gas_cell_scheduler != nullptr &&
                  payload.gas_cell_scheduler->isElementActive(static_cast<std::uint32_t>(i))
              ? 1U
              : 0U);
    });
    append_logical_lane(gas_scheduler_count, [&](std::size_t i) {
      return payload.gas_cell_scheduler != nullptr
          ? payload.gas_cell_scheduler->pendingBinIndex(static_cast<std::uint32_t>(i))
          : validation_state.cells.time_bin[i];
    });
  }
  append_u64(payload.output_cadence_state.output_enabled ? 1ull : 0ull);
  append_u64(payload.output_cadence_state.write_restarts ? 1ull : 0ull);
  append_u64(payload.output_cadence_state.snapshot_due ? 1ull : 0ull);
  append_u64(payload.output_cadence_state.checkpoint_due ? 1ull : 0ull);
  append_u64(payload.output_cadence_state.last_completed_step_index);
  append_u64(payload.output_cadence_state.snapshot_interval_steps);
  append_u64(payload.output_cadence_state.next_snapshot_step_index);
  if (include_output_time_cadence) {
    append_u64(std::bit_cast<std::uint64_t>(
        payload.output_cadence_state.snapshot_interval_time_code));
    append_u64(std::bit_cast<std::uint64_t>(
        payload.output_cadence_state.next_snapshot_time_code));
  }
  append_string(payload.output_cadence_state.snapshot_stem);
  append_string(payload.output_cadence_state.restart_stem);
  for (const StochasticModulePersistentState& module_state : sortedStochasticModules(payload.stochastic_state)) {
    append_string(module_state.module_name);
    append_u64(module_state.schema_version);
    append_string(module_state.rng_policy);
    append_u64(module_state.random_seed);
    append_u64(module_state.rank_local_seed_offset);
    append_u64(module_state.last_committed_step_index);
    append_u64(module_state.deterministic_from_serialized_inputs ? 1ull : 0ull);
  }
  // Hash the exact same canonical text as before, including its length
  // prefix, without retaining another population-scale serialization string.
  const std::size_t distributed_bytes =
      payload.distributed_gravity_state.serializedSizeBytes();
  append_u64(static_cast<std::uint64_t>(distributed_bytes));
  class DistributedRestartHashStreambuf final : public std::streambuf {
   public:
    DistributedRestartHashStreambuf(std::uint64_t& legacy_hash,
                                   core::internal::Sha256& strong_hash,
                                   std::size_t expected_bytes)
        : m_legacy_hash(legacy_hash), m_strong_hash(strong_hash),
          m_expected_bytes(expected_bytes) {}
    void finish() const {
      if (m_written != m_expected_bytes) {
        throw std::logic_error("distributed restart changed length during integrity hashing");
      }
    }
   protected:
    std::streamsize xsputn(const char* bytes, std::streamsize count) override {
      if (count < 0) throw std::length_error("negative distributed restart hash length");
      const std::size_t size = static_cast<std::size_t>(count);
      if (size > m_expected_bytes - m_written) {
        throw std::length_error("distributed restart exceeds integrity hash extent");
      }
      const auto* data = reinterpret_cast<const std::uint8_t*>(bytes);
      m_legacy_hash = fnv1aAppend(m_legacy_hash,
          {reinterpret_cast<const std::byte*>(bytes), size});
      m_strong_hash.update(data, size);
      m_written += size;
      return count;
    }
    int_type overflow(int_type ch) override {
      if (!traits_type::eq_int_type(ch, traits_type::eof())) {
        const char value = traits_type::to_char_type(ch);
        (void)xsputn(&value, 1);
      }
      return traits_type::not_eof(ch);
    }
   private:
    std::uint64_t& m_legacy_hash;
    core::internal::Sha256& m_strong_hash;
    std::size_t m_expected_bytes;
    std::size_t m_written = 0U;
  };
  DistributedRestartHashStreambuf distributed_buffer(hash, strong_hash, distributed_bytes);
  std::ostream distributed_stream(&distributed_buffer);
  distributed_stream.exceptions(std::ios::badbit | std::ios::failbit);
  payload.distributed_gravity_state.serializeTo(distributed_stream);
  distributed_buffer.finish();

  return RestartIntegrityDigests{hash, sha256DigestHex(strong_hash.finish())};
}

std::uint64_t restartPayloadIntegrityHash(
    const RestartWritePayload& payload, core::MemoryGovernor* governor) {
  return restartPayloadIntegrityDigestsImpl(payload, true, true, true, true, true, true, true, true, governor).legacy_fnv1a;
}

std::string restartPayloadIntegrityHashHex(
    const RestartWritePayload& payload, core::MemoryGovernor* governor) {
  return hexU64(restartPayloadIntegrityHash(payload, governor));
}


namespace {
// CHUI Part 2: true when every Part 2 physics key holds its default value
// (the only configuration a pre-Part-2 checkpoint can be migrated into).
[[nodiscard]] bool part2PhysicsKeysAtDefaults(const core::PhysicsConfig& physics) {
  const core::PhysicsConfig defaults{};
  return physics.thermochemistry_backend == defaults.thermochemistry_backend &&
      physics.metal_cooling_mode == defaults.metal_cooling_mode &&
      physics.metal_cooling_table_range == defaults.metal_cooling_table_range &&
      physics.dust_model == defaults.dust_model &&
      physics.ism_regime_model == defaults.ism_regime_model &&
      physics.sgs_turbulence_model == defaults.sgs_turbulence_model &&
      physics.chemistry_adaptive_selection == defaults.chemistry_adaptive_selection &&
      physics.stellar_population_model == defaults.stellar_population_model &&
      physics.fb_supernova_model == defaults.fb_supernova_model &&
      physics.thermochemistry_relative_tolerance == defaults.thermochemistry_relative_tolerance &&
      physics.thermochemistry_abundance_absolute_tolerance == defaults.thermochemistry_abundance_absolute_tolerance &&
      physics.thermochemistry_analytic_relaxation_time_years == defaults.thermochemistry_analytic_relaxation_time_years &&
      physics.thermochemistry_analytic_equilibrium_temperature_k == defaults.thermochemistry_analytic_equilibrium_temperature_k &&
      physics.isrf_g0_habing == defaults.isrf_g0_habing &&
      physics.chemistry_shielding_length_max_pc == defaults.chemistry_shielding_length_max_pc &&
      physics.metal_cooling_switch_temperature_k == defaults.metal_cooling_switch_temperature_k &&
      physics.chemistry_primordial_helium_mass_fraction == defaults.chemistry_primordial_helium_mass_fraction &&
      physics.dust_accretion_timescale_ref_years == defaults.dust_accretion_timescale_ref_years &&
      physics.dust_sputtering_timescale_ref_years == defaults.dust_sputtering_timescale_ref_years &&
      physics.dust_sn_swept_mass_msun == defaults.dust_sn_swept_mass_msun &&
      physics.dust_sn_destruction_efficiency == defaults.dust_sn_destruction_efficiency &&
      physics.dust_agb_condensation_efficiency == defaults.dust_agb_condensation_efficiency &&
      physics.dust_cc_condensation_efficiency == defaults.dust_cc_condensation_efficiency &&
      physics.ism_jeans_resolution_cells == defaults.ism_jeans_resolution_cells &&
      physics.ism_regime_hysteresis_fraction == defaults.ism_regime_hysteresis_fraction &&
      physics.ism_regime_min_dwell_years == defaults.ism_regime_min_dwell_years &&
      physics.sgs_dissipation_coefficient == defaults.sgs_dissipation_coefficient &&
      physics.sgs_eddy_viscosity_coefficient == defaults.sgs_eddy_viscosity_coefficient &&
      physics.sgs_feedback_energy_fraction == defaults.sgs_feedback_energy_fraction &&
      physics.chemistry_adaptive_hysteresis_fraction == defaults.chemistry_adaptive_hysteresis_fraction &&
      physics.chemistry_adaptive_min_dwell_years == defaults.chemistry_adaptive_min_dwell_years &&
      physics.chemistry_adaptive_cold_temperature_k == defaults.chemistry_adaptive_cold_temperature_k &&
      physics.chemistry_adaptive_dense_nh_cm3 == defaults.chemistry_adaptive_dense_nh_cm3 &&
      physics.chemistry_shadow_sample_fraction == defaults.chemistry_shadow_sample_fraction &&
      physics.ssp_imf_mass_min_msun == defaults.ssp_imf_mass_min_msun &&
      physics.ssp_imf_mass_max_msun == defaults.ssp_imf_mass_max_msun &&
      physics.ssp_cc_mass_min_msun == defaults.ssp_cc_mass_min_msun &&
      physics.ssp_cc_mass_max_msun == defaults.ssp_cc_mass_max_msun &&
      physics.ssp_snia_per_msun == defaults.ssp_snia_per_msun &&
      physics.ssp_snia_delay_min_years == defaults.ssp_snia_delay_min_years &&
      physics.ssp_snia_delay_max_years == defaults.ssp_snia_delay_max_years &&
      physics.ssp_energy_per_sn_erg == defaults.ssp_energy_per_sn_erg &&
      physics.ssp_reference_net_metal_yield == defaults.ssp_reference_net_metal_yield &&
      physics.fb_presn_wind_momentum_lsun_per_c == defaults.fb_presn_wind_momentum_lsun_per_c &&
      physics.thermochemistry_max_substeps == defaults.thermochemistry_max_substeps &&
      physics.thermochemistry_max_newton_iterations == defaults.thermochemistry_max_newton_iterations &&
      physics.chemistry_shadow_max_cells == defaults.chemistry_shadow_max_cells &&
      physics.uv_background_table_path == defaults.uv_background_table_path &&
      physics.uv_background_table_sha256 == defaults.uv_background_table_sha256 &&
      physics.metal_cooling_table_path == defaults.metal_cooling_table_path &&
      physics.metal_cooling_table_sha256 == defaults.metal_cooling_table_sha256 &&
      physics.ssp_yield_table_path == defaults.ssp_yield_table_path &&
      physics.ssp_yield_table_sha256 == defaults.ssp_yield_table_sha256 &&
      physics.chemistry_metallicity_scaled_dust_proxy == defaults.chemistry_metallicity_scaled_dust_proxy &&
      physics.thermochemistry_allow_unqualified_rates == defaults.thermochemistry_allow_unqualified_rates &&
      physics.ssp_stochastic_sn_events == defaults.ssp_stochastic_sn_events &&
      physics.ssp_clustered_associations == defaults.ssp_clustered_associations;
}
}  // namespace

std::string restartConfigurationCompatibilityRecord(
    const RestartReadResult& restart, const core::FrozenConfig& current) {
  const auto fail = [](std::string_view detail) -> void {
    throw std::runtime_error("restart normalized configuration compatibility refused: " + std::string(detail));
  };
  validateContinuationMetadata(restart.normalized_config_text,
      restart.normalized_config_hash_hex, restart.provenance, "restart configuration admission");
  const std::string begin_marker = "restart_config_compatibility_begin=treepm_normalized_defaults_v2\n";
  const std::string end_marker = "restart_config_compatibility_end=treepm_normalized_defaults_v2\n";
  if (restart.normalized_config_hash_hex == current.provenance.config_hash_hex) {
    if (restart.normalized_config_text != current.normalized_text) fail("hash agrees but normalized text differs");
    // Preserve a previous migration record across subsequent exact resumes.
    for (const std::string revision : {"treepm_normalized_defaults_v2", "treepm_normalized_defaults_v1"}) {
      const std::string historical_begin = "restart_config_compatibility_begin="+revision+"\n";
      const std::string historical_end = "restart_config_compatibility_end="+revision+"\n";
      const auto begin = restart.provenance.derived_runtime_state.find(historical_begin);
      if (begin == std::string::npos) continue;
      const auto end = restart.provenance.derived_runtime_state.find(historical_end,begin);
      if (end == std::string::npos) fail("unterminated historical compatibility record");
      return restart.provenance.derived_runtime_state.substr(begin,end+historical_end.size()-begin);
    }
    return {};
  }
  // This exact historical default projection is also valid for the current
  // fully verified v24 checkpoint. v24 adds output-cadence integrity but does
  // not relax normalized-config provenance or global-KDK admission. Refuse
  // older or unknown schemas rather than treating them as equivalent.
  const auto restart_version = restart.diagnostics.restart_schema_version;
  if (!restart.payload_integrity_verified ||
      (restart_version != k_restart_schema_v23 && restart_version != k_restart_schema_v24) ||
      restart.payload_integrity_algorithm != "sha256-canonical-le-v1" || restart.payload_integrity_sha256_hex.size() != 64U) {
    fail("migration requires a complete, integrity-verified v23/v24 checkpoint");
  }
  const auto& original = restart.provenance;
  if (original.schema_version != "provenance_v7" || original.config_schema_name != "cosmosim_config" ||
      original.config_schema_version != "1" || original.integrity_digest_algorithm != "sha256" ||
      original.normalized_config != restart.normalized_config_text ||
      original.config_hash_hex != restart.normalized_config_hash_hex ||
      original.normalized_config_hash_hex != restart.normalized_config_hash_hex ||
      original.normalized_config_sha256_hex != core::strongConfigHashSha256Hex(restart.normalized_config_text)) {
    fail("complete original normalized configuration and v7 provenance are required");
  }
  if (current.config.numerics.treepm_radial_coefficients_enabled ||
      current.config.numerics.treepm_full_residual_lookup_enabled ||
      current.config.numerics.treepm_leaf_simd_enabled ||
      current.config.numerics.treepm_residual_block_size != 64 ||
      current.config.numerics.treepm_tree_max_leaf_size != 16 ||
      current.config.numerics.treepm_gaussian_pair_lookup_enabled ||
      !current.config.numerics.treepm_full_mac_diagnostics || current.config.numerics.hierarchical_max_rung != 0 ||
      // CHUI Part 1 keys must hold their historical-equivalent defaults.
      current.config.numerics.hydro_reconstruction != core::HydroReconstructionMode::kPiecewiseLinear ||
      current.config.numerics.hydro_slope_limiter != core::HydroSlopeLimiterMode::kMonotonizedCentral ||
      current.config.numerics.hydro_integration_scheme != core::HydroIntegrationMode::kMusclHancock ||
      current.config.numerics.hydro_riemann_solver != core::HydroRiemannSolverMode::kHllc ||
      current.config.numerics.hydro_dual_energy != core::HydroDualEnergyConfigMode::kOff ||
      current.config.numerics.hydro_dual_energy_eta != 1.0e-3 ||
      current.config.parallel.node_memory_admission != core::NodeMemoryAdmissionMode::kOff ||
      current.config.parallel.node_memory_budget_bytes != 0U ||
      current.config.parallel.node_memory_reserve_bytes != 0U ||
      // CHUI Part 2 keys must hold their (physics-disabled) defaults.
      !part2PhysicsKeysAtDefaults(current.config.physics) ||
      current.config.schema_version != 1 ||
      current.provenance.config_hash_hex != core::stableConfigHashHex(current.normalized_text)) {
    fail("new options must retain legacy defaults and global KDK must remain selected");
  }
  if (restart.integrator_state.inside_kdk_step || !restart.integrator_state.last_completed_restart_safe ||
      !core::isRestartSafeBoundary(restart.integrator_state.last_completed_boundary_kind) ||
      restart.integrator_state.time_bins.max_bin != 0U ||
      restart.scheduler_state.max_bin != 0U || restart.gas_cell_scheduler_state.max_bin != 0U ||
      restart.scheduler_state.current_tick != restart.gas_cell_scheduler_state.current_tick ||
      restart.state.metadata.normalized_config_hash != core::stableConfigHash(restart.normalized_config_text) ||
      restart.state.metadata.normalized_config_hash_hex != restart.normalized_config_hash_hex) {
    fail("closed global KDK state and original state configuration identity required");
  }
  // This is a versioned literal projection, not permissive re-normalization:
  // it cannot accept another missing key, a changed scientific value, or an
  // unsupported historical default. Roundtrip also verifies today's defaults.
  if (current.normalized_text.size() > 1024U*1024U || restart.normalized_config_text.size() > 1024U*1024U) {
    fail("normalized configuration exceeds migration metadata budget");
  }
  std::string projected = current.normalized_text;
  const auto removeDefaults = [&](std::initializer_list<std::string_view> lines) {
    for (const auto line : lines) {
      const auto position = projected.find(line);
      if (position == std::string::npos || projected.find(line, position+line.size()) != std::string::npos)
        fail("current normalized default representation is not the supported revision");
      projected.erase(position, line.size());
    }
  };
  // CHUI Part 1 (20261103) added nine keys whose defaults reproduce the
  // historical hydro method and per-process memory semantics. Their current
  // values were verified above to be those defaults, so remove each complete
  // "key = value" line exactly once (value text comes from the same
  // normalizer, avoiding a second formatting rule). Every predecessor
  // checkpoint lacks this group.
  const auto removeKeyLines = [&](std::initializer_list<std::string_view> keys) {
    for (const auto key : keys) {
      const std::string prefix = "\n" + std::string(key) + " = ";
      const auto position = projected.find(prefix);
      if (position == std::string::npos || projected.find(prefix, position + prefix.size()) != std::string::npos)
        fail("current normalized Part 1 default representation is not the supported revision");
      const auto line_end = projected.find('\n', position + 1U);
      if (line_end == std::string::npos) fail("unterminated Part 1 normalized default line");
      projected.erase(position + 1U, line_end - position);
    }
  };
  // CHUI Part 2 (20261103) added the baryonic-physics key group below; all
  // were verified above to hold their defaults (no Part 2 physics active).
  if (projected.find("\nthermochemistry_backend = ") != std::string::npos &&
      restart.normalized_config_text.find("\nthermochemistry_backend = ") == std::string::npos) {
    removeKeyLines({
        "thermochemistry_backend",
        "metal_cooling_mode",
        "metal_cooling_table_range",
        "dust_model",
        "ism_regime_model",
        "sgs_turbulence_model",
        "chemistry_adaptive_selection",
        "stellar_population_model",
        "fb_supernova_model",
        "thermochemistry_relative_tolerance",
        "thermochemistry_abundance_absolute_tolerance",
        "thermochemistry_analytic_relaxation_time_years",
        "thermochemistry_analytic_equilibrium_temperature_k",
        "isrf_g0_habing",
        "chemistry_shielding_length_max_pc",
        "metal_cooling_switch_temperature_k",
        "chemistry_primordial_helium_mass_fraction",
        "dust_accretion_timescale_ref_years",
        "dust_sputtering_timescale_ref_years",
        "dust_sn_swept_mass_msun",
        "dust_sn_destruction_efficiency",
        "dust_agb_condensation_efficiency",
        "dust_cc_condensation_efficiency",
        "ism_jeans_resolution_cells",
        "ism_regime_hysteresis_fraction",
        "ism_regime_min_dwell_years",
        "sgs_dissipation_coefficient",
        "sgs_eddy_viscosity_coefficient",
        "sgs_feedback_energy_fraction",
        "chemistry_adaptive_hysteresis_fraction",
        "chemistry_adaptive_min_dwell_years",
        "chemistry_adaptive_cold_temperature_k",
        "chemistry_adaptive_dense_nh_cm3",
        "chemistry_shadow_sample_fraction",
        "ssp_imf_mass_min_msun",
        "ssp_imf_mass_max_msun",
        "ssp_cc_mass_min_msun",
        "ssp_cc_mass_max_msun",
        "ssp_snia_per_msun",
        "ssp_snia_delay_min_years",
        "ssp_snia_delay_max_years",
        "ssp_energy_per_sn_erg",
        "ssp_reference_net_metal_yield",
        "fb_presn_wind_momentum_lsun_per_c",
        "thermochemistry_max_substeps",
        "thermochemistry_max_newton_iterations",
        "chemistry_shadow_max_cells",
        "uv_background_table_path",
        "uv_background_table_sha256",
        "metal_cooling_table_path",
        "metal_cooling_table_sha256",
        "ssp_yield_table_path",
        "ssp_yield_table_sha256",
        "chemistry_metallicity_scaled_dust_proxy",
        "thermochemistry_allow_unqualified_rates",
        "ssp_stochastic_sn_events",
        "ssp_clustered_associations"});
  }
  if (projected.find("\nhydro_reconstruction = ") != std::string::npos &&
      restart.normalized_config_text.find("\nhydro_reconstruction = ") == std::string::npos) {
    removeKeyLines({"hydro_reconstruction", "hydro_slope_limiter", "hydro_integration_scheme",
        "hydro_riemann_solver", "hydro_dual_energy", "hydro_dual_energy_eta",
        "node_memory_admission", "node_memory_budget_bytes", "node_memory_reserve_bytes"});
  }
  const std::string after_part1_projection = projected;
  // Support exact predecessor groups, not arbitrary missing keys.
  if (projected != restart.normalized_config_text) {
    removeDefaults({"treepm_radial_coefficients_enabled = false\n", "treepm_full_residual_lookup_enabled = false\n", "treepm_leaf_simd_enabled = false\n",
        "treepm_residual_block_size = 64\n", "treepm_tree_max_leaf_size = 16\n"});
  }
  if (projected != restart.normalized_config_text) {
    removeDefaults({"treepm_gaussian_pair_lookup_enabled = false\n", "treepm_full_mac_diagnostics = true\n"});
  }
  if (projected != restart.normalized_config_text) {
    projected = after_part1_projection;
    removeDefaults({"treepm_gaussian_pair_lookup_enabled = false\n", "treepm_full_mac_diagnostics = true\n"});
  }
  if (projected != restart.normalized_config_text) fail("scientific/settings equivalence not established by exact legacy projection");
  if (core::loadFrozenConfigFromString(projected, "restart_treepm_normalized_defaults_v2").normalized_text != current.normalized_text) {
    fail("legacy default normalization does not reproduce current configuration");
  }
  const auto hex = [](std::string_view text) {
    constexpr char digits[] = "0123456789abcdef";
    std::string result;
    result.reserve(text.size()*2U);
    for (const char value : text) {
      const auto byte = static_cast<unsigned char>(value);
      result.push_back(digits[byte >> 4U]); result.push_back(digits[byte & 15U]);
    }
    return result;
  };
  const auto original_provenance = core::serializeProvenanceRecord(original);
  if (original_provenance.size() > 4U*1024U*1024U) fail("original provenance exceeds migration metadata budget");
  return begin_marker +
      "original_normalized_hash_hex=" + restart.normalized_config_hash_hex + "\n" +
      "original_normalized_sha256_hex=" + original.normalized_config_sha256_hex + "\n" +
      "original_payload_sha256_hex=" + restart.payload_integrity_sha256_hex + "\n" +
      "derived_normalized_hash_hex=" + current.provenance.config_hash_hex + "\n" +
      "derived_normalized_sha256_hex=" + core::strongConfigHashSha256Hex(current.normalized_text) + "\n" +
      "original_normalized_text_utf8_hex=" + hex(restart.normalized_config_text) + "\n" +
      "original_provenance_utf8_hex=" + hex(original_provenance) + "\n" + end_marker;
}

void writeRestartCheckpointHdf5(
    const std::filesystem::path& output_path,
    const RestartWritePayload& payload,
    const RestartWritePolicy& policy) {
#if !COSMOSIM_ENABLE_HDF5
  (void)output_path;
  (void)payload;
  (void)policy;
  throw std::runtime_error("restart checkpoint requires COSMOSIM_ENABLE_HDF5=ON");
#else
  if (payload.persistent_state.simulation_state == nullptr || payload.integrator_state == nullptr || payload.scheduler == nullptr) {
    throw std::invalid_argument("restart write payload must include state, integrator_state, and scheduler");
  }
  validateContinuationMetadata(
      payload.normalized_config_text,
      payload.normalized_config_hash_hex,
      payload.provenance,
      "restart writer");
  if (!payload.persistent_state.simulation_state->validateOwnershipInvariants()) {
    throw std::invalid_argument("cannot checkpoint invalid simulation state");
  }
  const bool uniform_particle_scheduler =
      payload.scheduler->representation() == core::SchedulerRepresentation::kUniformRungZero;
  core::TimeBinPersistentState scheduler_state_for_validation;
  scheduler_state_for_validation.current_tick = payload.scheduler->currentTick();
  scheduler_state_for_validation.max_bin = payload.scheduler->maxBin();
  if (!uniform_particle_scheduler) {
    scheduler_state_for_validation = payload.scheduler->exportPersistentState();
  }
  const core::TimeBinPersistentState gas_cell_scheduler_state_for_validation =
      payload.gas_cell_scheduler != nullptr
          ? payload.gas_cell_scheduler->exportPersistentState()
          : legacyGasCellSchedulerStateFromMirrors(
                *payload.persistent_state.simulation_state, scheduler_state_for_validation);
  const std::vector<std::uint64_t> gas_cell_scheduler_ids =
      gasCellIdsByDenseLocalRow(*payload.persistent_state.simulation_state, "restart writer");
  core::assertCanWriteCheckpointAtBoundary(*payload.integrator_state, scheduler_state_for_validation.current_tick);
  if (gas_cell_scheduler_state_for_validation.current_tick != scheduler_state_for_validation.current_tick) {
    throw std::invalid_argument("restart writer particle and gas-cell schedulers must share current_tick");
  }
  if (uniform_particle_scheduler) {
    if (payload.scheduler->elementCount() != payload.persistent_state.simulation_state->particles.size() ||
        payload.scheduler->maxBin() != 0U) {
      throw std::invalid_argument("restart writer uniform particle scheduler extent/rung contract is invalid");
    }
    for (std::size_t i = 0; i < payload.persistent_state.simulation_state->particles.size(); ++i) {
      if (payload.persistent_state.simulation_state->particleTimeBin(i) != payload.scheduler->binIndex(static_cast<std::uint32_t>(i))) {
        throw std::invalid_argument("restart writer particle time_bin mirror is stale relative to uniform scheduler authority");
      }
    }
  } else {
    validateRestartTimeBinMirrorsAgainstScheduler(
        *payload.persistent_state.simulation_state,
        scheduler_state_for_validation,
        "restart writer");
  }
  validateGasCellTimeBinMirrorsAgainstScheduler(
      *payload.persistent_state.simulation_state,
      gas_cell_scheduler_state_for_validation,
      gas_cell_scheduler_ids,
      "restart writer");
  validateOutputCadenceStateForRestart(
      payload.output_cadence_state, *payload.integrator_state, "restart writer");
  validateStochasticStateForRestart(
      payload.stochastic_state, *payload.integrator_state, "restart writer");
  const GravityForceCachePersistentView empty_force_cache{};
  const GravityForceCachePersistentView owning_force_cache_view =
      payload.gravity_force_cache != nullptr
          ? makeGravityForceCachePersistentView(*payload.gravity_force_cache)
          : empty_force_cache;
  const GravityForceCachePersistentView& gravity_force_cache =
      payload.gravity_force_cache_view != nullptr
          ? *payload.gravity_force_cache_view
          : owning_force_cache_view;
  validateGravityForceCacheForRestart(
      gravity_force_cache, *payload.persistent_state.simulation_state, "restart writer");

  internal::validateTemporarySuffix(policy.temporary_suffix);
  internal::TransactionalFileTarget transaction(
      output_path,
      policy.temporary_suffix,
      policy.enable_fsync_finalize
          ? internal::FileDurability::kDurablePublication
          : internal::FileDurability::kAtomicVisibility);

  Hdf5Handle file(H5Fcreate(
      transaction.temporaryPath().string().c_str(),
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed to create temporary restart file: " + transaction.temporaryPath().string());
  }

  const RestartIntegrityDigests integrity =
      restartPayloadIntegrityDigestsImpl(payload, true, true, true, true, true, true, true, true,
                                        policy.memory_governor);

  const auto& shared_names = sharedIoContractNames();
  writeScalarStringAttribute(
      file.get(), shared_names.file_kind_attribute, std::string(shared_names.restart_checkpoint_file_kind));
  writeScalarStringAttribute(file.get(), "restart_schema_name", restartSchema().name);
  writeScalarU32Attribute(file.get(), "restart_schema_version", restartSchema().version);
  writeScalarStringAttribute(file.get(), "normalized_config_hash_hex", payload.normalized_config_hash_hex);
  writeScalarStringAttribute(file.get(), "payload_integrity_hash_hex", hexU64(integrity.legacy_fnv1a));
  writeScalarU64Attribute(file.get(), "payload_integrity_hash", integrity.legacy_fnv1a);
  writeScalarStringAttribute(file.get(), "payload_integrity_algorithm", "sha256-canonical-le-v1");
  writeScalarStringAttribute(file.get(), "payload_integrity_sha256_hex", integrity.sha256_hex);
  writeStringDataset(
      file.get(),
      std::string(shared_names.normalized_config_text_dataset),
      payload.normalized_config_text);
  writeStringDataset(
      file.get(),
      std::string(shared_names.provenance_record_dataset),
      core::serializeProvenanceRecord(payload.provenance));

  writeStateGroup(file.get(), *payload.persistent_state.simulation_state);

  Hdf5Handle integrator_group(openOrCreateGroup(file.get(), "/integrator"));
  writeScalarF64Attribute(integrator_group.get(), "current_time_code", payload.integrator_state->current_time_code);
  writeScalarF64Attribute(integrator_group.get(), "current_scale_factor", payload.integrator_state->current_scale_factor);
  writeScalarF64Attribute(integrator_group.get(), "current_redshift", payload.integrator_state->current_redshift);
  writeScalarF64Attribute(integrator_group.get(), "current_hubble_rate_code", payload.integrator_state->current_hubble_rate_code);
  writeScalarF64Attribute(integrator_group.get(), "time_si_per_code", payload.integrator_state->time_si_per_code);
  writeScalarF64Attribute(integrator_group.get(), "dt_time_code", payload.integrator_state->dt_time_code);
  writeScalarF64Attribute(integrator_group.get(), "last_drift_factor_code", payload.integrator_state->last_drift_factor_code);
  writeScalarF64Attribute(integrator_group.get(), "last_first_kick_factor_code", payload.integrator_state->last_first_kick_factor_code);
  writeScalarF64Attribute(integrator_group.get(), "last_second_kick_factor_code", payload.integrator_state->last_second_kick_factor_code);
  writeScalarF64Attribute(integrator_group.get(), "last_first_hubble_drag_factor", payload.integrator_state->last_first_hubble_drag_factor);
  writeScalarF64Attribute(integrator_group.get(), "last_second_hubble_drag_factor", payload.integrator_state->last_second_hubble_drag_factor);
  writeScalarU64Attribute(integrator_group.get(), "step_index", payload.integrator_state->step_index);
  writeScalarU32Attribute(integrator_group.get(), "scheme", static_cast<std::uint32_t>(payload.integrator_state->scheme));
  writeScalarU32Attribute(integrator_group.get(), "current_boundary_kind", static_cast<std::uint32_t>(payload.integrator_state->current_boundary_kind));
  writeScalarU32Attribute(integrator_group.get(), "last_completed_boundary_kind", static_cast<std::uint32_t>(payload.integrator_state->last_completed_boundary_kind));
  writeScalarU32Attribute(integrator_group.get(), "inside_kdk_step", payload.integrator_state->inside_kdk_step ? 1U : 0U);
  writeScalarU32Attribute(integrator_group.get(), "last_completed_restart_safe", payload.integrator_state->last_completed_restart_safe ? 1U : 0U);
  writeScalarU32Attribute(integrator_group.get(), "time_bins_hierarchical", payload.integrator_state->time_bins.hierarchical_enabled ? 1U : 0U);
  writeScalarU32Attribute(integrator_group.get(), "time_bins_active_bin", payload.integrator_state->time_bins.active_bin);
  writeScalarU32Attribute(integrator_group.get(), "time_bins_max_bin", payload.integrator_state->time_bins.max_bin);
  writeScalarU32Attribute(integrator_group.get(), "pm_long_range_field_valid", payload.integrator_state->pm_long_range_field_valid ? 1U : 0U);
  writeScalarU32Attribute(integrator_group.get(), "pm_refresh_enabled", payload.integrator_state->pm_refresh_enabled ? 1U : 0U);
  const core::PmSynchronizationPersistentState pm_sync_state = payload.integrator_state->pm_sync_state.exportPersistentState();
  writeScalarU64Attribute(integrator_group.get(), "pm_cadence_steps", pm_sync_state.cadence_steps);
  writeScalarU64Attribute(integrator_group.get(), "pm_gravity_kick_opportunity", pm_sync_state.gravity_kick_opportunity);
  writeScalarU64Attribute(integrator_group.get(), "pm_last_refresh_opportunity", pm_sync_state.last_refresh_opportunity);
  writeScalarU64Attribute(integrator_group.get(), "pm_field_version", pm_sync_state.field_version);
  writeScalarU64Attribute(integrator_group.get(), "pm_last_refresh_step_index", pm_sync_state.last_refresh_step_index);
  writeScalarF64Attribute(integrator_group.get(), "pm_last_refresh_scale_factor", pm_sync_state.last_refresh_scale_factor);
  writeScalarU32Attribute(integrator_group.get(), "pm_refresh_commit_pending", pm_sync_state.refresh_commit_pending ? 1U : 0U);
  writeScalarU64Attribute(integrator_group.get(), "pm_pending_refresh_opportunity", pm_sync_state.pending_refresh_opportunity);
  writeScalarU64Attribute(integrator_group.get(), "pm_pending_refresh_field_version", pm_sync_state.pending_refresh_field_version);

  Hdf5Handle force_cache_group(openOrCreateGroup(file.get(), "/gravity_force_cache"));
  writeScalarU32Attribute(force_cache_group.get(), "valid", gravity_force_cache.valid ? 1U : 0U);
  writeDataset1d(force_cache_group.get(), "particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gravity_force_cache.particle_id);
  writeDataset1d(force_cache_group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gravity_force_cache.gas_cell_id);
  writeDataset1d(force_cache_group.get(), "particle_accel_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.particle_accel_x_comoving);
  writeDataset1d(force_cache_group.get(), "particle_accel_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.particle_accel_y_comoving);
  writeDataset1d(force_cache_group.get(), "particle_accel_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.particle_accel_z_comoving);
  writeDataset1d(force_cache_group.get(), "cell_accel_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.cell_accel_x_comoving);
  writeDataset1d(force_cache_group.get(), "cell_accel_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.cell_accel_y_comoving);
  writeDataset1d(force_cache_group.get(), "cell_accel_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.cell_accel_z_comoving);
  Hdf5Handle endpoint_group(openOrCreateGroup(file.get(), std::string(k_hierarchical_endpoint_group)));
  writeScalarU32Attribute(endpoint_group.get(), "valid", gravity_force_cache.hierarchical_endpoint_valid ? 1U : 0U);
  writeScalarU64Attribute(endpoint_group.get(), "pm_field_version", gravity_force_cache.hierarchical_pm_field_version);
  writeDataset1d(endpoint_group.get(), "tree_accel_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.hierarchical_tree_accel_x_comoving);
  writeDataset1d(endpoint_group.get(), "tree_accel_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.hierarchical_tree_accel_y_comoving);
  writeDataset1d(endpoint_group.get(), "tree_accel_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.hierarchical_tree_accel_z_comoving);
  writeDataset1d(endpoint_group.get(), "pm_accel_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.hierarchical_pm_accel_x_comoving);
  writeDataset1d(endpoint_group.get(), "pm_accel_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.hierarchical_pm_accel_y_comoving);
  writeDataset1d(endpoint_group.get(), "pm_accel_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.hierarchical_pm_accel_z_comoving);
  writeDataset1d(endpoint_group.get(), "nearest_pair_r2_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.hierarchical_nearest_pair_r2_comoving);

  const core::TimeBinPersistentState scheduler_state = scheduler_state_for_validation;
  Hdf5Handle scheduler_group(openOrCreateGroup(file.get(), "/scheduler"));
  writeScalarU64Attribute(scheduler_group.get(), "current_tick", payload.scheduler->currentTick());
  writeScalarU32Attribute(scheduler_group.get(), "max_bin", payload.scheduler->maxBin());
  if (uniform_particle_scheduler) {
    const std::size_t scheduler_count = payload.scheduler->elementCount();
    writeUniformDataset1d<std::uint8_t>(scheduler_group.get(), "bin_index", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_count, 0U);
    writeUniformDataset1d<std::uint64_t>(scheduler_group.get(), "next_activation_tick", H5T_STD_U64LE, H5T_NATIVE_UINT64, scheduler_count, payload.scheduler->currentTick());
    writeUniformDataset1d<std::uint8_t>(scheduler_group.get(), "active_flag", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_count, 0U);
    const std::uint8_t uniform_pending = scheduler_count == 0U
        ? core::HierarchicalTimeBinScheduler::k_unset_pending_bin
        : payload.scheduler->pendingBinIndex(0U);
    writeUniformDataset1d<std::uint8_t>(scheduler_group.get(), "pending_bin_index", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_count, uniform_pending);
  } else {
    writeDataset1d(scheduler_group.get(), "bin_index", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_state.bin_index);
    writeDataset1d(scheduler_group.get(), "next_activation_tick", H5T_STD_U64LE, H5T_NATIVE_UINT64, scheduler_state.next_activation_tick);
    writeDataset1d(scheduler_group.get(), "active_flag", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_state.active_flag);
    writeDataset1d(scheduler_group.get(), "pending_bin_index", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_state.pending_bin_index);
  }

  Hdf5Handle gas_scheduler_group(openOrCreateGroup(file.get(), "/gas_cell_scheduler"));
  writeScalarStringAttribute(
      gas_scheduler_group.get(), "identity_key", std::string(k_gas_cell_scheduler_identity_key));
  writeScalarU64Attribute(
      gas_scheduler_group.get(), "current_tick", gas_cell_scheduler_state_for_validation.current_tick);
  writeScalarU32Attribute(
      gas_scheduler_group.get(), "max_bin", gas_cell_scheduler_state_for_validation.max_bin);
  writeDataset1d(
      gas_scheduler_group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_scheduler_ids);
  writeDataset1d(
      gas_scheduler_group.get(),
      "bin_index",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      gas_cell_scheduler_state_for_validation.bin_index);
  writeDataset1d(
      gas_scheduler_group.get(),
      "next_activation_tick",
      H5T_STD_U64LE,
      H5T_NATIVE_UINT64,
      gas_cell_scheduler_state_for_validation.next_activation_tick);
  writeDataset1d(
      gas_scheduler_group.get(),
      "active_flag",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      gas_cell_scheduler_state_for_validation.active_flag);
  writeDataset1d(
      gas_scheduler_group.get(),
      "pending_bin_index",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      gas_cell_scheduler_state_for_validation.pending_bin_index);
  writeOutputCadenceGroup(file.get(), payload.output_cadence_state);
  writeStochasticStateGroup(file.get(), payload.stochastic_state);
  RestartDiagnosticsSummary restart_diagnostics = makeRestartDiagnosticsSummary(
      *payload.integrator_state,
      scheduler_state,
      gas_cell_scheduler_state_for_validation,
      payload.output_cadence_state,
      payload.stochastic_state);
  if (uniform_particle_scheduler) {
    restart_diagnostics.scheduler_current_tick = payload.scheduler->currentTick();
    restart_diagnostics.scheduler_max_bin = payload.scheduler->maxBin();
    restart_diagnostics.scheduler_element_count = payload.scheduler->elementCount();
    restart_diagnostics.scheduler_active_count = 0U;
    restart_diagnostics.scheduler_pending_transition_count = 0U;
  }
  writeRestartDiagnosticsGroup(file.get(), restart_diagnostics);
  writeDistributedGravityGroup(file.get(), payload.distributed_gravity_state);

  if (H5Fflush(file.get(), H5F_SCOPE_GLOBAL) < 0) {
    throw std::runtime_error("failed to flush temporary restart file");
  }
  file = Hdf5Handle();
  transaction.publish();
#endif
}

std::uint64_t restartReadCandidateStagingBytes(
    const RestartReadCandidateDimensions& dimensions) {
  // Conservative per-row widths for the dominant candidate families. These are
  // physical capacity upper bounds admitted BEFORE any population vector is
  // constructed. Checked arithmetic rejects overflow before allocation.
  const auto mul = [](std::uint64_t count, std::uint64_t width, std::string_view label) {
    if (width != 0U && count > std::numeric_limits<std::uint64_t>::max() / width) {
      throw std::overflow_error(std::string(label) + ": restart candidate extent overflows uint64");
    }
    return count * width;
  };
  // Particles: 7x f64 SoA + 1x u8 time_bin + sidecar (2x u64 + 2x u32 + 3x f64)
  // + gravity force cache (u64 id + 3x f64) + v24 hierarchical endpoint (7x f64).
  constexpr std::uint64_t k_particle_bytes =
      7U * sizeof(double) + sizeof(std::uint8_t) + 2U * sizeof(std::uint64_t) +
      2U * sizeof(std::uint32_t) + 3U * sizeof(double) +
      sizeof(std::uint64_t) + 3U * sizeof(double) + 7U * sizeof(double);
  // Cells: 4x f64 + u8 + u32 + gas sidecar (2x u64 + 9x f64).
  constexpr std::uint64_t k_cell_bytes =
      4U * sizeof(double) + sizeof(std::uint8_t) + sizeof(std::uint32_t) +
      2U * sizeof(std::uint64_t) + 9U * sizeof(double);
  // Patches: coarse patch SoA estimate.
  constexpr std::uint64_t k_patch_bytes =
      3U * sizeof(std::uint64_t) + sizeof(std::int32_t) + 3U * sizeof(std::uint32_t) +
      6U * sizeof(double) + 3U * sizeof(std::uint16_t) + sizeof(std::uint32_t);
  // Pending flux register record + temporal cell record estimates.
  constexpr std::uint64_t k_pending_bytes = 32U * sizeof(double) + 8U * sizeof(std::uint64_t);
  constexpr std::uint64_t k_temporal_bytes = 16U * sizeof(double) + 4U * sizeof(std::uint64_t);
  std::uint64_t total = mul(dimensions.particle_count, k_particle_bytes, "restart particle candidate");
  total = core::checkedMemoryBytesAdd(
      total, mul(dimensions.cell_count, k_cell_bytes, "restart cell candidate"),
      "restart candidate cells");
  total = core::checkedMemoryBytesAdd(
      total, mul(dimensions.patch_count, k_patch_bytes, "restart patch candidate"),
      "restart candidate patches");
  total = core::checkedMemoryBytesAdd(
      total, mul(dimensions.pending_flux_count, k_pending_bytes, "restart pending candidate"),
      "restart candidate pending flux");
  total = core::checkedMemoryBytesAdd(
      total,
      mul(dimensions.temporal_cell_count, k_temporal_bytes, "restart temporal candidate"),
      "restart candidate temporal history");
  return total;
}

#if COSMOSIM_ENABLE_HDF5
[[nodiscard]] std::uint64_t hdf5DatasetLength1dChecked(hid_t loc, const char* path, std::string_view label) {
  Hdf5Handle dataset(H5Dopen2(loc, path, H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("restart preflight missing dataset: " + std::string(path));
  }
  Hdf5Handle space(H5Dget_space(dataset.get()));
  if (!space.valid() || H5Sget_simple_extent_ndims(space.get()) != 1) {
    throw std::runtime_error("restart preflight expected 1D dataset: " + std::string(path));
  }
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_dims(space.get(), dims, nullptr) != 1) {
    throw std::runtime_error("restart preflight cannot read extent: " + std::string(path));
  }
  // hsize_t is 64-bit; reject values that cannot fit size_t allocation.
  const std::uint64_t length = static_cast<std::uint64_t>(dims[0]);
  if (length > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw std::length_error(std::string(label) + ": restart dimension exceeds addressable size_t");
  }
  // Impossible-dimension guard before any allocation: cap at 2^40 rows (~1T).
  constexpr std::uint64_t k_max_plausible_rows = std::uint64_t{1} << 40U;
  if (length > k_max_plausible_rows) {
    throw std::length_error(std::string(label) + ": restart dimension is implausible before allocation");
  }
  return length;
}

[[nodiscard]] RestartReadCandidateDimensions preflightRestartCandidateDimensions(
    hid_t file, const RestartReadPolicy& policy) {
  RestartReadCandidateDimensions dims;
  dims.particle_count =
      hdf5DatasetLength1dChecked(file, "/state/particles/position_x_comoving", "particles");
  dims.cell_count =
      hdf5DatasetLength1dChecked(file, "/state/cells/center_x_comoving", "cells");
  dims.patch_count =
      hdf5DatasetLength1dChecked(file, "/state/patches/patch_id", "patches");
  // Optional groups: missing links mean zero count (older schemas).
  if (H5Lexists(file, "/state/amr_pending_flux_registers/register_key", H5P_DEFAULT) > 0) {
    dims.pending_flux_count = hdf5DatasetLength1dChecked(
        file, "/state/amr_pending_flux_registers/register_key", "pending flux");
  }
  if (H5Lexists(file, "/state/amr_temporal_boundary_history/cell_count", H5P_DEFAULT) > 0) {
    // Temporal history cell count is per-record cells; use record count as bound.
    dims.temporal_cell_count = hdf5DatasetLength1dChecked(
        file, "/state/amr_temporal_boundary_history/patch_id", "temporal history");
  }
  if (policy.enforce_dimension_consistency) {
    // Cross-dataset consistency BEFORE allocation: all lanes of a family must
    // agree, otherwise reject instead of truncating.
    const std::uint64_t py =
        hdf5DatasetLength1dChecked(file, "/state/particles/position_y_comoving", "particles.y");
    const std::uint64_t pid =
        hdf5DatasetLength1dChecked(file, "/state/particle_sidecar/particle_id", "particle_id");
    if (py != dims.particle_count || pid != dims.particle_count) {
      throw std::runtime_error("restart dimension inconsistency: particle lanes disagree");
    }
    const std::uint64_t gas =
        hdf5DatasetLength1dChecked(file, "/state/gas_cells/gas_cell_id", "gas_cells");
    if (gas != dims.cell_count) {
      throw std::runtime_error("restart dimension inconsistency: cell/gas-cell lanes disagree");
    }
  }
  return dims;
}
#endif


namespace {

#if COSMOSIM_ENABLE_HDF5
template <typename T, typename Getter>
void verifyDataset1dChunked(
    hid_t file,
    std::string_view path,
    std::size_t expected_count,
    Getter getter,
    hid_t memory_type,
    std::uint64_t workspace_limit_bytes,
    RestartVerificationResult& result) {
  Hdf5Handle dataset(H5Dopen2(file, std::string(path).c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("restart verification missing dataset: " + std::string(path));
  }
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  if (!file_space.valid() || H5Sget_simple_extent_ndims(file_space.get()) != 1) {
    throw std::runtime_error("restart verification expected 1D dataset: " + std::string(path));
  }
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_dims(file_space.get(), dims, nullptr) != 1 ||
      dims[0] != static_cast<hsize_t>(expected_count)) {
    throw std::runtime_error("restart verification extent mismatch: " + std::string(path));
  }
  const std::size_t capacity = std::max<std::size_t>(
      1U, static_cast<std::size_t>(workspace_limit_bytes / std::max<std::size_t>(sizeof(T), 1U)));
  std::vector<T> chunk(std::min(expected_count, capacity));
  result.workspace_high_water_bytes = std::max<std::uint64_t>(
      result.workspace_high_water_bytes,
      static_cast<std::uint64_t>(chunk.capacity()) * sizeof(T));
  for (std::size_t offset = 0; offset < expected_count;) {
    const std::size_t count = std::min(chunk.size(), expected_count - offset);
    hsize_t start[1] = {static_cast<hsize_t>(offset)};
    hsize_t extent[1] = {static_cast<hsize_t>(count)};
    Hdf5Handle selected_space(H5Dget_space(dataset.get()));
    Hdf5Handle memory_space(H5Screate_simple(1, extent, nullptr));
    if (!selected_space.valid() || !memory_space.valid() ||
        H5Sselect_hyperslab(selected_space.get(), H5S_SELECT_SET, start, nullptr, extent, nullptr) < 0 ||
        H5Dread(dataset.get(), memory_type, memory_space.get(), selected_space.get(), H5P_DEFAULT, chunk.data()) < 0) {
      throw std::runtime_error("restart verification failed reading dataset: " + std::string(path));
    }
    for (std::size_t local = 0; local < count; ++local) {
      const T observed = chunk[local];
      const T expected = static_cast<T>(getter(offset + local));
      bool equal = false;
      if constexpr (std::is_floating_point_v<T>) {
        equal = std::bit_cast<std::uint64_t>(static_cast<double>(observed)) ==
            std::bit_cast<std::uint64_t>(static_cast<double>(expected));
      } else {
        equal = observed == expected;
      }
      if (!equal) {
        result.failed_field = std::string(path);
        result.detail = "first mismatch at row " + std::to_string(offset + local);
        throw std::runtime_error(
            "restart verification mismatch: " + std::string(path) +
            " row=" + std::to_string(offset + local));
      }
    }
    result.verified_bytes += static_cast<std::uint64_t>(count) * sizeof(T);
    offset += count;
  }
  ++result.verified_dataset_count;
}

template <typename T>
void verifyDatasetSpanChunked(
    hid_t file, std::string_view path, std::span<const T> expected, hid_t memory_type,
    std::uint64_t workspace_limit_bytes, RestartVerificationResult& result) {
  verifyDataset1dChunked<T>(
      file, path, expected.size(),
      [&](std::size_t i) { return expected[i]; }, memory_type,
      workspace_limit_bytes, result);
}

void verifyStringDatasetChunked(
    hid_t file, std::string_view path, std::string_view expected,
    std::uint64_t workspace_limit_bytes, RestartVerificationResult& result) {
  verifyDataset1dChunked<std::uint8_t>(
      file, path, expected.size(),
      [&](std::size_t i) { return static_cast<std::uint8_t>(expected[i]); }, H5T_NATIVE_UINT8,
      workspace_limit_bytes, result);
}

void verifyDistributedRestartStateChunked(
    hid_t file,
    const parallel::DistributedRestartState& expected,
    std::uint64_t workspace_limit_bytes,
    RestartVerificationResult& result) {
  Hdf5Handle dataset(H5Dopen2(file, "/distributed_gravity/state", H5P_DEFAULT));
  if (!dataset.valid()) throw std::runtime_error("restart verification missing /distributed_gravity/state");
  Hdf5Handle space(H5Dget_space(dataset.get()));
  hsize_t dims[1] = {0};
  if (!space.valid() || H5Sget_simple_extent_ndims(space.get()) != 1 ||
      H5Sget_simple_extent_dims(space.get(), dims, nullptr) != 1 ||
      dims[0] != static_cast<hsize_t>(expected.serializedSizeBytes())) {
    throw std::runtime_error("restart verification distributed state extent mismatch");
  }
  const std::size_t capacity = std::max<std::size_t>(1U, static_cast<std::size_t>(workspace_limit_bytes));
  std::vector<char> chunk(std::min<std::size_t>(capacity, expected.serializedSizeBytes()));
  result.workspace_high_water_bytes = std::max<std::uint64_t>(result.workspace_high_water_bytes, chunk.capacity());
  class ComparingStreambuf final : public std::streambuf {
   public:
    ComparingStreambuf(hid_t dataset, std::vector<char>& chunk, RestartVerificationResult& result)
        : m_dataset(dataset), m_chunk(chunk), m_result(result) {}
    void finish(std::size_t expected_bytes) const {
      if (m_offset != expected_bytes) throw std::runtime_error("distributed restart serialization length changed during verification");
    }
   protected:
    std::streamsize xsputn(const char* bytes, std::streamsize count) override {
      if (count < 0) throw std::length_error("negative distributed restart verification write");
      std::size_t remaining = static_cast<std::size_t>(count);
      const char* expected_bytes = bytes;
      while (remaining != 0U) {
        const std::size_t n = std::min(remaining, m_chunk.size());
        hsize_t start[1] = {static_cast<hsize_t>(m_offset)};
        hsize_t extent[1] = {static_cast<hsize_t>(n)};
        Hdf5Handle file_space(H5Dget_space(m_dataset));
        Hdf5Handle memory_space(H5Screate_simple(1, extent, nullptr));
        if (!file_space.valid() || !memory_space.valid() ||
            H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, start, nullptr, extent, nullptr) < 0 ||
            H5Dread(m_dataset, H5T_NATIVE_UINT8, memory_space.get(), file_space.get(), H5P_DEFAULT, m_chunk.data()) < 0) {
          throw std::runtime_error("failed streaming /distributed_gravity/state");
        }
        if (std::memcmp(m_chunk.data(), expected_bytes, n) != 0) {
          m_result.failed_field = "/distributed_gravity/state";
          m_result.detail = "serialized byte mismatch near offset " + std::to_string(m_offset);
          throw std::runtime_error("restart verification distributed state byte mismatch");
        }
        m_offset += n;
        expected_bytes += n;
        remaining -= n;
        m_result.verified_bytes += n;
      }
      return count;
    }
    int_type overflow(int_type ch) override {
      if (!traits_type::eq_int_type(ch, traits_type::eof())) {
        const char value = traits_type::to_char_type(ch);
        (void)xsputn(&value, 1);
      }
      return traits_type::not_eof(ch);
    }
   private:
    hid_t m_dataset;
    std::vector<char>& m_chunk;
    RestartVerificationResult& m_result;
    std::size_t m_offset = 0U;
  };
  ComparingStreambuf buffer(dataset.get(), chunk, result);
  std::ostream stream(&buffer);
  stream.exceptions(std::ios::badbit | std::ios::failbit);
  expected.serializeTo(stream);
  buffer.finish(expected.serializedSizeBytes());
  ++result.verified_dataset_count;
}
#endif

}  // namespace

RestartVerificationResult verifyRestartCheckpointHdf5(
    const std::filesystem::path& input_path,
    const RestartWritePayload& expected,
    const RestartVerificationPolicy& policy) {
  RestartVerificationResult result;
#if !COSMOSIM_ENABLE_HDF5
  (void)input_path; (void)expected; (void)policy;
  throw std::runtime_error("restart checkpoint verification requires COSMOSIM_ENABLE_HDF5=ON");
#else
  if (expected.persistent_state.simulation_state == nullptr || expected.integrator_state == nullptr ||
      expected.scheduler == nullptr) {
    throw std::invalid_argument("restart verifier requires state, integrator, and scheduler authority");
  }
  if (policy.workspace_limit_bytes < sizeof(std::uint64_t)) {
    throw std::invalid_argument("restart verification workspace limit is too small");
  }
  core::MemoryReservation workspace_reservation;
  if (policy.memory_governor != nullptr) {
    workspace_reservation = policy.memory_governor->reserve(
        core::MemoryClass::kDiagnostic, policy.workspace_limit_bytes, "io.restart.streaming_verify");
    workspace_reservation.commit();
  }
  try {
    Hdf5Handle file(H5Fopen(input_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    if (!file.valid()) throw std::runtime_error("failed to open restart checkpoint for streaming verification: " + input_path.string());
    validateFileKindAttribute(file.get(), sharedIoContractNames().restart_checkpoint_file_kind, "restart verifier");
    result.schema_name = readScalarStringAttribute(file.get(), "restart_schema_name");
    result.schema_version = readScalarU32Attribute(file.get(), "restart_schema_version");
    if (result.schema_name != restartSchema().name || result.schema_version != restartSchema().version) {
      throw std::runtime_error("write-time restart verifier requires the current restart schema");
    }
    validateRestartCheckpointSchema(file.get(), result.schema_version);
    const auto require_u32_attr = [&](hid_t loc, std::string_view key, std::uint32_t expected_value, std::string_view field) {
      if (readScalarU32Attribute(loc, key) != expected_value) { result.failed_field = std::string(field); throw std::runtime_error("restart verification attribute mismatch: " + std::string(field)); }
    };
    const auto require_u64_attr = [&](hid_t loc, std::string_view key, std::uint64_t expected_value, std::string_view field) {
      if (readScalarU64Attribute(loc, key) != expected_value) { result.failed_field = std::string(field); throw std::runtime_error("restart verification attribute mismatch: " + std::string(field)); }
    };
    const auto require_f64_attr = [&](hid_t loc, std::string_view key, double expected_value, std::string_view field) {
      if (std::bit_cast<std::uint64_t>(readScalarF64Attribute(loc, key)) != std::bit_cast<std::uint64_t>(expected_value)) { result.failed_field = std::string(field); throw std::runtime_error("restart verification attribute mismatch: " + std::string(field)); }
    };
    const auto require_string_attr = [&](hid_t loc, std::string_view key, const std::string& expected_value, std::string_view field) {
      if (readScalarStringAttribute(loc, key) != expected_value) { result.failed_field = std::string(field); throw std::runtime_error("restart verification attribute mismatch: " + std::string(field)); }
    };

    require_string_attr(file.get(), "normalized_config_hash_hex", expected.normalized_config_hash_hex, "/@normalized_config_hash_hex");

    Hdf5Handle integrator_group(openRequiredGroup(file.get(), "/integrator"));
    const auto& integ = *expected.integrator_state;
    require_f64_attr(integrator_group.get(), "current_time_code", integ.current_time_code, "/integrator/@current_time_code");
    require_f64_attr(integrator_group.get(), "current_scale_factor", integ.current_scale_factor, "/integrator/@current_scale_factor");
    require_f64_attr(integrator_group.get(), "current_redshift", integ.current_redshift, "/integrator/@current_redshift");
    require_f64_attr(integrator_group.get(), "current_hubble_rate_code", integ.current_hubble_rate_code, "/integrator/@current_hubble_rate_code");
    require_f64_attr(integrator_group.get(), "time_si_per_code", integ.time_si_per_code, "/integrator/@time_si_per_code");
    require_f64_attr(integrator_group.get(), "dt_time_code", integ.dt_time_code, "/integrator/@dt_time_code");
    require_f64_attr(integrator_group.get(), "last_drift_factor_code", integ.last_drift_factor_code, "/integrator/@last_drift_factor_code");
    require_f64_attr(integrator_group.get(), "last_first_kick_factor_code", integ.last_first_kick_factor_code, "/integrator/@last_first_kick_factor_code");
    require_f64_attr(integrator_group.get(), "last_second_kick_factor_code", integ.last_second_kick_factor_code, "/integrator/@last_second_kick_factor_code");
    require_f64_attr(integrator_group.get(), "last_first_hubble_drag_factor", integ.last_first_hubble_drag_factor, "/integrator/@last_first_hubble_drag_factor");
    require_f64_attr(integrator_group.get(), "last_second_hubble_drag_factor", integ.last_second_hubble_drag_factor, "/integrator/@last_second_hubble_drag_factor");
    require_u64_attr(integrator_group.get(), "step_index", integ.step_index, "/integrator/@step_index");
    require_u32_attr(integrator_group.get(), "scheme", static_cast<std::uint32_t>(integ.scheme), "/integrator/@scheme");
    require_u32_attr(integrator_group.get(), "current_boundary_kind", static_cast<std::uint32_t>(integ.current_boundary_kind), "/integrator/@current_boundary_kind");
    require_u32_attr(integrator_group.get(), "last_completed_boundary_kind", static_cast<std::uint32_t>(integ.last_completed_boundary_kind), "/integrator/@last_completed_boundary_kind");
    require_u32_attr(integrator_group.get(), "inside_kdk_step", integ.inside_kdk_step ? 1U : 0U, "/integrator/@inside_kdk_step");
    require_u32_attr(integrator_group.get(), "last_completed_restart_safe", integ.last_completed_restart_safe ? 1U : 0U, "/integrator/@last_completed_restart_safe");
    require_u32_attr(integrator_group.get(), "time_bins_hierarchical", integ.time_bins.hierarchical_enabled ? 1U : 0U, "/integrator/@time_bins_hierarchical");
    require_u32_attr(integrator_group.get(), "time_bins_active_bin", integ.time_bins.active_bin, "/integrator/@time_bins_active_bin");
    require_u32_attr(integrator_group.get(), "time_bins_max_bin", integ.time_bins.max_bin, "/integrator/@time_bins_max_bin");
    require_u32_attr(integrator_group.get(), "pm_long_range_field_valid", integ.pm_long_range_field_valid ? 1U : 0U, "/integrator/@pm_long_range_field_valid");
    require_u32_attr(integrator_group.get(), "pm_refresh_enabled", integ.pm_refresh_enabled ? 1U : 0U, "/integrator/@pm_refresh_enabled");
    const core::PmSynchronizationPersistentState pm_sync = integ.pm_sync_state.exportPersistentState();
    require_u64_attr(integrator_group.get(), "pm_cadence_steps", pm_sync.cadence_steps, "/integrator/@pm_cadence_steps");
    require_u64_attr(integrator_group.get(), "pm_gravity_kick_opportunity", pm_sync.gravity_kick_opportunity, "/integrator/@pm_gravity_kick_opportunity");
    require_u64_attr(integrator_group.get(), "pm_last_refresh_opportunity", pm_sync.last_refresh_opportunity, "/integrator/@pm_last_refresh_opportunity");
    require_u64_attr(integrator_group.get(), "pm_field_version", pm_sync.field_version, "/integrator/@pm_field_version");
    require_u64_attr(integrator_group.get(), "pm_last_refresh_step_index", pm_sync.last_refresh_step_index, "/integrator/@pm_last_refresh_step_index");
    require_f64_attr(integrator_group.get(), "pm_last_refresh_scale_factor", pm_sync.last_refresh_scale_factor, "/integrator/@pm_last_refresh_scale_factor");
    require_u32_attr(integrator_group.get(), "pm_refresh_commit_pending", pm_sync.refresh_commit_pending ? 1U : 0U, "/integrator/@pm_refresh_commit_pending");
    require_u64_attr(integrator_group.get(), "pm_pending_refresh_opportunity", pm_sync.pending_refresh_opportunity, "/integrator/@pm_pending_refresh_opportunity");
    require_u64_attr(integrator_group.get(), "pm_pending_refresh_field_version", pm_sync.pending_refresh_field_version, "/integrator/@pm_pending_refresh_field_version");

    Hdf5Handle output_group(openRequiredGroup(file.get(), "/output_cadence"));
    const auto& out = expected.output_cadence_state;
    require_u32_attr(output_group.get(), "output_enabled", out.output_enabled ? 1U : 0U, "/output_cadence/@output_enabled");
    require_u32_attr(output_group.get(), "write_restarts", out.write_restarts ? 1U : 0U, "/output_cadence/@write_restarts");
    require_u32_attr(output_group.get(), "snapshot_due", out.snapshot_due ? 1U : 0U, "/output_cadence/@snapshot_due");
    require_u32_attr(output_group.get(), "checkpoint_due", out.checkpoint_due ? 1U : 0U, "/output_cadence/@checkpoint_due");
    require_u64_attr(output_group.get(), "last_completed_step_index", out.last_completed_step_index, "/output_cadence/@last_completed_step_index");
    require_u64_attr(output_group.get(), "snapshot_interval_steps", out.snapshot_interval_steps, "/output_cadence/@snapshot_interval_steps");
    require_u64_attr(output_group.get(), "next_snapshot_step_index", out.next_snapshot_step_index, "/output_cadence/@next_snapshot_step_index");
    require_f64_attr(output_group.get(), "snapshot_interval_time_code", out.snapshot_interval_time_code, "/output_cadence/@snapshot_interval_time_code");
    require_f64_attr(output_group.get(), "next_snapshot_time_code", out.next_snapshot_time_code, "/output_cadence/@next_snapshot_time_code");
    require_string_attr(output_group.get(), "snapshot_stem", out.snapshot_stem, "/output_cadence/@snapshot_stem");
    require_string_attr(output_group.get(), "restart_stem", out.restart_stem, "/output_cadence/@restart_stem");

    const core::SimulationState& state = *expected.persistent_state.simulation_state;
    const std::uint64_t w = policy.workspace_limit_bytes;
#define VERIFY_SPAN(PATH, TYPE, MEMTYPE, VALUES) \
    verifyDatasetSpanChunked<TYPE>(file.get(), PATH, std::span<const TYPE>((VALUES).data(), (VALUES).size()), MEMTYPE, w, result)
    VERIFY_SPAN("/state/particles/position_x_comoving", double, H5T_NATIVE_DOUBLE, state.particles.position_x_comoving);
    VERIFY_SPAN("/state/particles/position_y_comoving", double, H5T_NATIVE_DOUBLE, state.particles.position_y_comoving);
    VERIFY_SPAN("/state/particles/position_z_comoving", double, H5T_NATIVE_DOUBLE, state.particles.position_z_comoving);
    VERIFY_SPAN("/state/particles/velocity_x_peculiar", double, H5T_NATIVE_DOUBLE, state.particles.velocity_x_peculiar);
    VERIFY_SPAN("/state/particles/velocity_y_peculiar", double, H5T_NATIVE_DOUBLE, state.particles.velocity_y_peculiar);
    VERIFY_SPAN("/state/particles/velocity_z_peculiar", double, H5T_NATIVE_DOUBLE, state.particles.velocity_z_peculiar);
    VERIFY_SPAN("/state/particles/mass_code", double, H5T_NATIVE_DOUBLE, state.particles.mass_code);
    verifyDataset1dChunked<std::uint8_t>(file.get(), "/state/particles/time_bin", state.particles.size(),
        [&](std::size_t i){ return state.particleTimeBin(i); }, H5T_NATIVE_UINT8, w, result);
    VERIFY_SPAN("/state/particle_sidecar/particle_id", std::uint64_t, H5T_NATIVE_UINT64, state.particle_sidecar.particle_id);
    verifyDataset1dChunked<std::uint64_t>(file.get(), "/state/particle_sidecar/sfc_key", state.particles.size(),
        [&](std::size_t i){ return state.particleSfcKey(i); }, H5T_NATIVE_UINT64, w, result);
    verifyDataset1dChunked<std::uint32_t>(file.get(), "/state/particle_sidecar/species_tag", state.particles.size(),
        [&](std::size_t i){ return state.particleSpeciesTag(i); }, H5T_NATIVE_UINT32, w, result);
    verifyDataset1dChunked<std::uint32_t>(file.get(), "/state/particle_sidecar/particle_flags", state.particles.size(),
        [&](std::size_t i){ return state.particleFlags(i); }, H5T_NATIVE_UINT32, w, result);
    verifyDataset1dChunked<std::uint32_t>(file.get(), "/state/particle_sidecar/owning_rank", state.particles.size(),
        [&](std::size_t i){ return state.particleOwningRank(i); }, H5T_NATIVE_UINT32, w, result);
    verifyDataset1dChunked<double>(file.get(), "/state/particle_sidecar/last_drift_time_code", state.particles.size(),
        [&](std::size_t i){ return state.particleLastDriftTimeCode(i); }, H5T_NATIVE_DOUBLE, w, result);
    verifyDataset1dChunked<double>(file.get(), "/state/particle_sidecar/last_drift_scale_factor", state.particles.size(),
        [&](std::size_t i){ return state.particleLastDriftScaleFactor(i); }, H5T_NATIVE_DOUBLE, w, result);
    VERIFY_SPAN("/state/particle_sidecar/gravity_softening_comoving", double, H5T_NATIVE_DOUBLE, state.particle_sidecar.gravity_softening_comoving);
    VERIFY_SPAN("/state/particle_sidecar/has_gravity_softening_override", std::uint8_t, H5T_NATIVE_UINT8, state.particle_sidecar.has_gravity_softening_override);

    VERIFY_SPAN("/state/cells/center_x_comoving", double, H5T_NATIVE_DOUBLE, state.cells.center_x_comoving);
    VERIFY_SPAN("/state/cells/center_y_comoving", double, H5T_NATIVE_DOUBLE, state.cells.center_y_comoving);
    VERIFY_SPAN("/state/cells/center_z_comoving", double, H5T_NATIVE_DOUBLE, state.cells.center_z_comoving);
    VERIFY_SPAN("/state/cells/mass_code", double, H5T_NATIVE_DOUBLE, state.cells.mass_code);
    VERIFY_SPAN("/state/cells/time_bin", std::uint8_t, H5T_NATIVE_UINT8, state.cells.time_bin);
    VERIFY_SPAN("/state/cells/patch_index", std::uint32_t, H5T_NATIVE_UINT32, state.cells.patch_index);

    VERIFY_SPAN("/state/gas_cells/gas_cell_id", std::uint64_t, H5T_NATIVE_UINT64, state.gas_cells.gas_cell_id);
    VERIFY_SPAN("/state/gas_cells/parent_particle_id", std::uint64_t, H5T_NATIVE_UINT64, state.gas_cells.parent_particle_id);
    VERIFY_SPAN("/state/gas_cells/velocity_x_peculiar", double, H5T_NATIVE_DOUBLE, state.gas_cells.velocity_x_peculiar);
    VERIFY_SPAN("/state/gas_cells/velocity_y_peculiar", double, H5T_NATIVE_DOUBLE, state.gas_cells.velocity_y_peculiar);
    VERIFY_SPAN("/state/gas_cells/velocity_z_peculiar", double, H5T_NATIVE_DOUBLE, state.gas_cells.velocity_z_peculiar);
    VERIFY_SPAN("/state/gas_cells/density_code", double, H5T_NATIVE_DOUBLE, state.gas_cells.density_code);
    VERIFY_SPAN("/state/gas_cells/pressure_code", double, H5T_NATIVE_DOUBLE, state.gas_cells.pressure_code);
    VERIFY_SPAN("/state/gas_cells/internal_energy_code", double, H5T_NATIVE_DOUBLE, state.gas_cells.internal_energy_code);
    VERIFY_SPAN("/state/gas_cells/metal_mass_code", double, H5T_NATIVE_DOUBLE, state.gas_cells.metal_mass_code);
    VERIFY_SPAN("/state/gas_cells/temperature_code", double, H5T_NATIVE_DOUBLE, state.gas_cells.temperature_code);
    VERIFY_SPAN("/state/gas_cells/sound_speed_code", double, H5T_NATIVE_DOUBLE, state.gas_cells.sound_speed_code);

    VERIFY_SPAN("/state/patches/patch_id", std::uint64_t, H5T_NATIVE_UINT64, state.patches.patch_id);
    VERIFY_SPAN("/state/patches/level", std::int32_t, H5T_NATIVE_INT32, state.patches.level);
    VERIFY_SPAN("/state/patches/first_cell", std::uint32_t, H5T_NATIVE_UINT32, state.patches.first_cell);
    VERIFY_SPAN("/state/patches/cell_count", std::uint32_t, H5T_NATIVE_UINT32, state.patches.cell_count);
    VERIFY_SPAN("/state/patches/parent_patch_id", std::uint64_t, H5T_NATIVE_UINT64, state.patches.parent_patch_id);
    VERIFY_SPAN("/state/patches/morton_key", std::uint64_t, H5T_NATIVE_UINT64, state.patches.morton_key);
    VERIFY_SPAN("/state/patches/origin_x_comoving", double, H5T_NATIVE_DOUBLE, state.patches.origin_x_comoving);
    VERIFY_SPAN("/state/patches/origin_y_comoving", double, H5T_NATIVE_DOUBLE, state.patches.origin_y_comoving);
    VERIFY_SPAN("/state/patches/origin_z_comoving", double, H5T_NATIVE_DOUBLE, state.patches.origin_z_comoving);
    VERIFY_SPAN("/state/patches/extent_x_comoving", double, H5T_NATIVE_DOUBLE, state.patches.extent_x_comoving);
    VERIFY_SPAN("/state/patches/extent_y_comoving", double, H5T_NATIVE_DOUBLE, state.patches.extent_y_comoving);
    VERIFY_SPAN("/state/patches/extent_z_comoving", double, H5T_NATIVE_DOUBLE, state.patches.extent_z_comoving);
    VERIFY_SPAN("/state/patches/cell_dim_x", std::uint16_t, H5T_NATIVE_UINT16, state.patches.cell_dim_x);
    VERIFY_SPAN("/state/patches/cell_dim_y", std::uint16_t, H5T_NATIVE_UINT16, state.patches.cell_dim_y);
    VERIFY_SPAN("/state/patches/cell_dim_z", std::uint16_t, H5T_NATIVE_UINT16, state.patches.cell_dim_z);
    VERIFY_SPAN("/state/patches/owning_rank", std::uint32_t, H5T_NATIVE_UINT32, state.patches.owning_rank);

    verifyDataset1dChunked<std::uint64_t>(file.get(), "/state/gas_cell_identity/gas_cell_id", state.gas_cell_identity.size(),
        [&](std::size_t i){ const auto* r = state.gas_cell_identity.findByLocalRow(static_cast<std::uint32_t>(i)); if (r == nullptr) throw std::runtime_error("missing gas identity local row during restart verification"); return r->gas_cell_id; }, H5T_NATIVE_UINT64, w, result);
    verifyDataset1dChunked<std::uint8_t>(file.get(), "/state/gas_cell_identity/has_parent_particle", state.gas_cell_identity.size(),
        [&](std::size_t i){ const auto* r = state.gas_cell_identity.findByLocalRow(static_cast<std::uint32_t>(i)); return static_cast<std::uint8_t>(r != nullptr && r->parent_particle_id.has_value() ? 1U : 0U); }, H5T_NATIVE_UINT8, w, result);
    verifyDataset1dChunked<std::uint64_t>(file.get(), "/state/gas_cell_identity/parent_particle_id", state.gas_cell_identity.size(),
        [&](std::size_t i){ const auto* r = state.gas_cell_identity.findByLocalRow(static_cast<std::uint32_t>(i)); return r == nullptr ? 0ULL : r->parent_particle_id.value_or(0ULL); }, H5T_NATIVE_UINT64, w, result);
    verifyDataset1dChunked<std::uint64_t>(file.get(), "/state/gas_cell_identity/owning_patch_id", state.gas_cell_identity.size(),
        [&](std::size_t i){ const auto* r = state.gas_cell_identity.findByLocalRow(static_cast<std::uint32_t>(i)); return r == nullptr ? 0ULL : r->owning_patch_id; }, H5T_NATIVE_UINT64, w, result);
    verifyDataset1dChunked<std::uint32_t>(file.get(), "/state/gas_cell_identity/local_cell_row", state.gas_cell_identity.size(),
        [&](std::size_t i){ return static_cast<std::uint32_t>(i); }, H5T_NATIVE_UINT32, w, result);

    const auto pending_records = state.pending_flux_registers.records();
#define VERIFY_PENDING(NAME, TYPE, MEMTYPE, EXPR) \
    verifyDataset1dChunked<TYPE>(file.get(), "/state/amr_pending_flux_registers/" NAME, pending_records.size(), \
        [&](std::size_t i){ return static_cast<TYPE>(pending_records[i].EXPR); }, MEMTYPE, w, result)
    VERIFY_PENDING("register_key", std::uint64_t, H5T_NATIVE_UINT64, register_key);
    VERIFY_PENDING("coarse_patch_id", std::uint64_t, H5T_NATIVE_UINT64, coarse_patch_id);
    VERIFY_PENDING("coarse_gas_cell_id", std::uint64_t, H5T_NATIVE_UINT64, coarse_gas_cell_id);
    VERIFY_PENDING("coarse_cell_index", std::uint64_t, H5T_NATIVE_UINT64, coarse_cell_index);
    VERIFY_PENDING("level", std::uint8_t, H5T_NATIVE_UINT8, level);
    VERIFY_PENDING("axis", std::uint8_t, H5T_NATIVE_UINT8, axis);
    VERIFY_PENDING("orientation", std::uint8_t, H5T_NATIVE_UINT8, orientation);
    VERIFY_PENDING("expected_area_comov", double, H5T_NATIVE_DOUBLE, expected_area_comov);
    VERIFY_PENDING("coarse_area_accumulated_comov", double, H5T_NATIVE_DOUBLE, coarse_area_accumulated_comov);
    VERIFY_PENDING("fine_area_accumulated_comov", double, H5T_NATIVE_DOUBLE, fine_area_accumulated_comov);
    VERIFY_PENDING("interval_start_code", double, H5T_NATIVE_DOUBLE, interval_start_code);
    VERIFY_PENDING("interval_end_code", double, H5T_NATIVE_DOUBLE, interval_end_code);
    VERIFY_PENDING("coarse_dt_code", double, H5T_NATIVE_DOUBLE, coarse_dt_code);
    VERIFY_PENDING("expected_fine_substeps", std::uint32_t, H5T_NATIVE_UINT32, expected_fine_substeps);
    VERIFY_PENDING("completed_fine_substeps", std::uint32_t, H5T_NATIVE_UINT32, completed_fine_substeps);
    VERIFY_PENDING("fine_substep_coverage_mask", std::uint64_t, H5T_NATIVE_UINT64, fine_substep_coverage_mask);
    VERIFY_PENDING("coarse_face_count", std::uint32_t, H5T_NATIVE_UINT32, coarse_face_count);
    VERIFY_PENDING("fine_face_count", std::uint32_t, H5T_NATIVE_UINT32, fine_face_count);
    VERIFY_PENDING("gas_cell_identity_generation", std::uint64_t, H5T_NATIVE_UINT64, gas_cell_identity_generation);
    VERIFY_PENDING("patch_geometry_generation", std::uint64_t, H5T_NATIVE_UINT64, patch_geometry_generation);
    VERIFY_PENDING("coarse_mass_flux_integral_code", double, H5T_NATIVE_DOUBLE, coarse_mass_flux_integral_code);
    VERIFY_PENDING("coarse_momentum_x_flux_integral_code", double, H5T_NATIVE_DOUBLE, coarse_momentum_x_flux_integral_code);
    VERIFY_PENDING("coarse_momentum_y_flux_integral_code", double, H5T_NATIVE_DOUBLE, coarse_momentum_y_flux_integral_code);
    VERIFY_PENDING("coarse_momentum_z_flux_integral_code", double, H5T_NATIVE_DOUBLE, coarse_momentum_z_flux_integral_code);
    VERIFY_PENDING("coarse_total_energy_flux_integral_code", double, H5T_NATIVE_DOUBLE, coarse_total_energy_flux_integral_code);
    VERIFY_PENDING("coarse_metal_mass_flux_integral_code", double, H5T_NATIVE_DOUBLE, coarse_metal_mass_flux_integral_code);
    VERIFY_PENDING("fine_mass_flux_integral_code", double, H5T_NATIVE_DOUBLE, fine_mass_flux_integral_code);
    VERIFY_PENDING("fine_momentum_x_flux_integral_code", double, H5T_NATIVE_DOUBLE, fine_momentum_x_flux_integral_code);
    VERIFY_PENDING("fine_momentum_y_flux_integral_code", double, H5T_NATIVE_DOUBLE, fine_momentum_y_flux_integral_code);
    VERIFY_PENDING("fine_momentum_z_flux_integral_code", double, H5T_NATIVE_DOUBLE, fine_momentum_z_flux_integral_code);
    VERIFY_PENDING("fine_total_energy_flux_integral_code", double, H5T_NATIVE_DOUBLE, fine_total_energy_flux_integral_code);
    VERIFY_PENDING("fine_metal_mass_flux_integral_code", double, H5T_NATIVE_DOUBLE, fine_metal_mass_flux_integral_code);
#undef VERIFY_PENDING

    const auto temporal_records = state.amr_temporal_boundary_history.records();
#define VERIFY_TEMPORAL_RECORD(NAME, TYPE, MEMTYPE, EXPR) \
    verifyDataset1dChunked<TYPE>(file.get(), "/state/amr_temporal_boundary_history/" NAME, temporal_records.size(), \
        [&](std::size_t i){ return static_cast<TYPE>(temporal_records[i].EXPR); }, MEMTYPE, w, result)
    VERIFY_TEMPORAL_RECORD("patch_id", std::uint64_t, H5T_NATIVE_UINT64, patch_id);
    VERIFY_TEMPORAL_RECORD("patch_level", std::uint8_t, H5T_NATIVE_UINT8, patch_level);
    VERIFY_TEMPORAL_RECORD("patch_geometry_fingerprint", std::uint64_t, H5T_NATIVE_UINT64, patch_geometry_fingerprint);
    VERIFY_TEMPORAL_RECORD("gas_cell_identity_generation", std::uint64_t, H5T_NATIVE_UINT64, gas_cell_identity_generation);
    VERIFY_TEMPORAL_RECORD("interval_start_code", double, H5T_NATIVE_DOUBLE, interval_start_code);
    VERIFY_TEMPORAL_RECORD("interval_end_code", double, H5T_NATIVE_DOUBLE, interval_end_code);
    verifyDataset1dChunked<std::uint8_t>(file.get(), "/state/amr_temporal_boundary_history/end_state_valid", temporal_records.size(),
        [&](std::size_t i){ return static_cast<std::uint8_t>(temporal_records[i].end_state_valid ? 1U : 0U); }, H5T_NATIVE_UINT8, w, result);
    std::size_t temporal_offset_next = 0U;
    std::uint64_t temporal_offset = 0U;
    verifyDataset1dChunked<std::uint64_t>(file.get(), "/state/amr_temporal_boundary_history/cell_offset", temporal_records.size(),
        [&](std::size_t i){
          if (i != temporal_offset_next) throw std::logic_error("temporal verification accessor order changed");
          const std::uint64_t value = temporal_offset;
          temporal_offset += static_cast<std::uint64_t>(temporal_records[i].cells.size());
          ++temporal_offset_next;
          return value;
        }, H5T_NATIVE_UINT64, w, result);
    verifyDataset1dChunked<std::uint64_t>(file.get(), "/state/amr_temporal_boundary_history/cell_count", temporal_records.size(),
        [&](std::size_t i){ return static_cast<std::uint64_t>(temporal_records[i].cells.size()); }, H5T_NATIVE_UINT64, w, result);
    std::size_t temporal_cell_total = 0U;
    for (const auto& record : temporal_records) temporal_cell_total += record.cells.size();
    const auto verify_temporal_cells = [&](std::string_view name, auto member, hid_t memory_type) {
      using T = std::decay_t<decltype(std::declval<core::AmrTemporalBoundaryHistoryCellRecord>().*member)>;
      std::size_t record_index = 0U;
      std::size_t cell_index = 0U;
      verifyDataset1dChunked<T>(file.get(), std::string("/state/amr_temporal_boundary_history/") + std::string(name), temporal_cell_total,
          [&](std::size_t) {
            while (record_index < temporal_records.size() && cell_index == temporal_records[record_index].cells.size()) {
              ++record_index; cell_index = 0U;
            }
            if (record_index >= temporal_records.size()) throw std::logic_error("temporal verification flattened accessor overflow");
            return temporal_records[record_index].cells[cell_index++].*member;
          }, memory_type, w, result);
    };
    verify_temporal_cells("gas_cell_id", &core::AmrTemporalBoundaryHistoryCellRecord::gas_cell_id, H5T_NATIVE_UINT64);
    verify_temporal_cells("patch_local_cell", &core::AmrTemporalBoundaryHistoryCellRecord::patch_local_cell, H5T_NATIVE_UINT64);
    verify_temporal_cells("start_mass_density_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::start_mass_density_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("start_momentum_density_x_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::start_momentum_density_x_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("start_momentum_density_y_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::start_momentum_density_y_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("start_momentum_density_z_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::start_momentum_density_z_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("start_total_energy_density_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::start_total_energy_density_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("start_metal_mass_density_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::start_metal_mass_density_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("end_mass_density_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::end_mass_density_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("end_momentum_density_x_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::end_momentum_density_x_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("end_momentum_density_y_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::end_momentum_density_y_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("end_momentum_density_z_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::end_momentum_density_z_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("end_total_energy_density_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::end_total_energy_density_comoving, H5T_NATIVE_DOUBLE);
    verify_temporal_cells("end_metal_mass_density_comoving", &core::AmrTemporalBoundaryHistoryCellRecord::end_metal_mass_density_comoving, H5T_NATIVE_DOUBLE);
#undef VERIFY_TEMPORAL_RECORD

    VERIFY_SPAN("/state/star_particles/particle_index", std::uint32_t, H5T_NATIVE_UINT32, state.star_particles.particle_index);
    VERIFY_SPAN("/state/star_particles/formation_scale_factor", double, H5T_NATIVE_DOUBLE, state.star_particles.formation_scale_factor);
    VERIFY_SPAN("/state/star_particles/birth_mass_code", double, H5T_NATIVE_DOUBLE, state.star_particles.birth_mass_code);
    VERIFY_SPAN("/state/star_particles/metallicity_mass_fraction", double, H5T_NATIVE_DOUBLE, state.star_particles.metallicity_mass_fraction);
    VERIFY_SPAN("/state/star_particles/birth_key", std::uint64_t, H5T_NATIVE_UINT64, state.star_particles.birth_key);
    VERIFY_SPAN("/state/star_particles/parent_gas_cell_id", std::uint64_t, H5T_NATIVE_UINT64, state.star_particles.parent_gas_cell_id);
    VERIFY_SPAN("/state/star_particles/birth_tick", std::uint64_t, H5T_NATIVE_UINT64, state.star_particles.birth_tick);
    VERIFY_SPAN("/state/star_particles/birth_ordinal", std::uint32_t, H5T_NATIVE_UINT32, state.star_particles.birth_ordinal);
    VERIFY_SPAN("/state/star_particles/stellar_age_years_last", double, H5T_NATIVE_DOUBLE, state.star_particles.stellar_age_years_last);
    VERIFY_SPAN("/state/star_particles/stellar_returned_mass_cumulative_code", double, H5T_NATIVE_DOUBLE, state.star_particles.stellar_returned_mass_cumulative_code);
    VERIFY_SPAN("/state/star_particles/stellar_returned_metals_cumulative_code", double, H5T_NATIVE_DOUBLE, state.star_particles.stellar_returned_metals_cumulative_code);
    VERIFY_SPAN("/state/star_particles/stellar_newly_synthesized_metals_cumulative_code", double, H5T_NATIVE_DOUBLE, state.star_particles.stellar_newly_synthesized_metals_cumulative_code);
    VERIFY_SPAN("/state/star_particles/enrichment_carry_mass_code", double, H5T_NATIVE_DOUBLE, state.star_particles.enrichment_carry_mass_code);
    VERIFY_SPAN("/state/star_particles/enrichment_carry_metals_code", double, H5T_NATIVE_DOUBLE, state.star_particles.enrichment_carry_metals_code);
    VERIFY_SPAN("/state/star_particles/enrichment_carry_feedback_energy_erg", double, H5T_NATIVE_DOUBLE, state.star_particles.enrichment_carry_feedback_energy_erg);
    VERIFY_SPAN("/state/star_particles/enrichment_carry_momentum_code", double, H5T_NATIVE_DOUBLE, state.star_particles.enrichment_carry_momentum_code);
    VERIFY_SPAN("/state/star_particles/stellar_deposited_mass_cumulative_code", double, H5T_NATIVE_DOUBLE, state.star_particles.stellar_deposited_mass_cumulative_code);
    VERIFY_SPAN("/state/star_particles/stellar_deposited_metals_cumulative_code", double, H5T_NATIVE_DOUBLE, state.star_particles.stellar_deposited_metals_cumulative_code);
    VERIFY_SPAN("/state/star_particles/stellar_deposited_feedback_energy_cumulative_erg", double, H5T_NATIVE_DOUBLE, state.star_particles.stellar_deposited_feedback_energy_cumulative_erg);
    VERIFY_SPAN("/state/star_particles/stellar_feedback_energy_cumulative_erg", double, H5T_NATIVE_DOUBLE, state.star_particles.stellar_feedback_energy_cumulative_erg);
    for (std::size_t channel = 0; channel < state.star_particles.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
      const std::string suffix = std::to_string(channel);
      verifyDatasetSpanChunked<double>(file.get(), "/state/star_particles/stellar_returned_mass_channel_cumulative_code_" + suffix, state.star_particles.stellar_returned_mass_channel_cumulative_code[channel], H5T_NATIVE_DOUBLE, w, result);
      verifyDatasetSpanChunked<double>(file.get(), "/state/star_particles/stellar_returned_metals_channel_cumulative_code_" + suffix, state.star_particles.stellar_returned_metals_channel_cumulative_code[channel], H5T_NATIVE_DOUBLE, w, result);
      verifyDatasetSpanChunked<double>(file.get(), "/state/star_particles/stellar_feedback_energy_channel_cumulative_erg_" + suffix, state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel], H5T_NATIVE_DOUBLE, w, result);
    }

    VERIFY_SPAN("/state/black_holes/particle_index", std::uint32_t, H5T_NATIVE_UINT32, state.black_holes.particle_index);
    VERIFY_SPAN("/state/black_holes/host_cell_index", std::uint32_t, H5T_NATIVE_UINT32, state.black_holes.host_cell_index);
    VERIFY_SPAN("/state/black_holes/subgrid_mass_code", double, H5T_NATIVE_DOUBLE, state.black_holes.subgrid_mass_code);
    VERIFY_SPAN("/state/black_holes/accretion_rate_code", double, H5T_NATIVE_DOUBLE, state.black_holes.accretion_rate_code);
    VERIFY_SPAN("/state/black_holes/feedback_energy_code", double, H5T_NATIVE_DOUBLE, state.black_holes.feedback_energy_code);
    VERIFY_SPAN("/state/black_holes/eddington_ratio", double, H5T_NATIVE_DOUBLE, state.black_holes.eddington_ratio);
    VERIFY_SPAN("/state/black_holes/cumulative_accreted_mass_code", double, H5T_NATIVE_DOUBLE, state.black_holes.cumulative_accreted_mass_code);
    VERIFY_SPAN("/state/black_holes/cumulative_feedback_energy_code", double, H5T_NATIVE_DOUBLE, state.black_holes.cumulative_feedback_energy_code);
    VERIFY_SPAN("/state/black_holes/duty_cycle_active_time_code", double, H5T_NATIVE_DOUBLE, state.black_holes.duty_cycle_active_time_code);
    VERIFY_SPAN("/state/black_holes/duty_cycle_total_time_code", double, H5T_NATIVE_DOUBLE, state.black_holes.duty_cycle_total_time_code);
    VERIFY_SPAN("/state/tracers/particle_index", std::uint32_t, H5T_NATIVE_UINT32, state.tracers.particle_index);
    VERIFY_SPAN("/state/tracers/parent_particle_id", std::uint64_t, H5T_NATIVE_UINT64, state.tracers.parent_particle_id);
    VERIFY_SPAN("/state/tracers/injection_step", std::uint64_t, H5T_NATIVE_UINT64, state.tracers.injection_step);
    VERIFY_SPAN("/state/tracers/host_cell_index", std::uint32_t, H5T_NATIVE_UINT32, state.tracers.host_cell_index);
    VERIFY_SPAN("/state/tracers/mass_fraction_of_host", double, H5T_NATIVE_DOUBLE, state.tracers.mass_fraction_of_host);
    VERIFY_SPAN("/state/tracers/last_host_mass_code", double, H5T_NATIVE_DOUBLE, state.tracers.last_host_mass_code);
    VERIFY_SPAN("/state/tracers/cumulative_exchanged_mass_code", double, H5T_NATIVE_DOUBLE, state.tracers.cumulative_exchanged_mass_code);
    verifyDatasetSpanChunked<std::uint64_t>(file.get(), "/state/species_count_by_species",
        std::span<const std::uint64_t>(state.species.count_by_species.data(), state.species.count_by_species.size()),
        H5T_NATIVE_UINT64, w, result);
#undef VERIFY_SPAN

    verifyStringDatasetChunked(file.get(), "/state/state_metadata", state.metadata.serialize(), w, result);
    verifyStringDatasetChunked(file.get(), std::string(sharedIoContractNames().normalized_config_text_dataset), expected.normalized_config_text, w, result);
    verifyStringDatasetChunked(file.get(), std::string(sharedIoContractNames().provenance_record_dataset), core::serializeProvenanceRecord(expected.provenance), w, result);
    std::string expected_module_names;
    for (const core::ModuleSidecarBlock* block : state.sidecars.blocksSortedByName()) {
      expected_module_names += block->module_name;
      expected_module_names.push_back('\n');
      const std::string base = "/state/module_sidecars/" + block->module_name;
      Hdf5Handle module_group(openRequiredGroup(file.get(), base));
      require_u32_attr(module_group.get(), "schema_version", block->schema_version, base + "/@schema_version");
      require_u32_attr(module_group.get(), "particle_indexed", block->particle_indexed ? 1U : 0U, base + "/@particle_indexed");
      require_u32_attr(module_group.get(), "row_stride_bytes", block->row_stride_bytes, base + "/@row_stride_bytes");
      require_u32_attr(module_group.get(), "required_species_mask", block->required_species_mask, base + "/@required_species_mask");
      require_u32_attr(module_group.get(), "requirement_kind", static_cast<std::uint32_t>(block->requirement.kind), base + "/@requirement_kind");
      require_u32_attr(module_group.get(), "requirement_species_mask", block->requirement.species_mask, base + "/@requirement_species_mask");
      require_u32_attr(module_group.get(), "requirement_particle_flags_mask", block->requirement.particle_flags_mask, base + "/@requirement_particle_flags_mask");
      require_f64_attr(module_group.get(), "requirement_threshold_code", block->requirement.threshold_code, base + "/@requirement_threshold_code");
      verifyDataset1dChunked<std::uint8_t>(file.get(), base + "/payload", block->payload.size(),
          [&](std::size_t i){ return std::to_integer<std::uint8_t>(block->payload[i]); }, H5T_NATIVE_UINT8, w, result);
      if (block->isParticleIndexed()) {
        verifyDatasetSpanChunked<std::uint64_t>(file.get(), base + "/particle_id_by_row", block->particle_id_by_row, H5T_NATIVE_UINT64, w, result);
      }
    }
    verifyStringDatasetChunked(file.get(), "/state/module_sidecar_names", expected_module_names, w, result);

    Hdf5Handle stochastic_group(openRequiredGroup(file.get(), "/stochastic_state"));
    const auto stochastic_modules = sortedStochasticModules(expected.stochastic_state);
    require_u32_attr(stochastic_group.get(), "module_count", static_cast<std::uint32_t>(stochastic_modules.size()), "/stochastic_state/@module_count");
    std::string stochastic_names;
    for (const auto& module : stochastic_modules) {
      if (!stochastic_names.empty()) stochastic_names.push_back('\n');
      stochastic_names += module.module_name;
      const std::string base = "/stochastic_state/" + module.module_name;
      Hdf5Handle module_group(openRequiredGroup(file.get(), base));
      require_u32_attr(module_group.get(), "schema_version", module.schema_version, base + "/@schema_version");
      require_string_attr(module_group.get(), "rng_policy", module.rng_policy, base + "/@rng_policy");
      require_u64_attr(module_group.get(), "random_seed", module.random_seed, base + "/@random_seed");
      require_u32_attr(module_group.get(), "rank_local_seed_offset", module.rank_local_seed_offset, base + "/@rank_local_seed_offset");
      require_u64_attr(module_group.get(), "last_committed_step_index", module.last_committed_step_index, base + "/@last_committed_step_index");
      require_u32_attr(module_group.get(), "deterministic_from_serialized_inputs", module.deterministic_from_serialized_inputs ? 1U : 0U, base + "/@deterministic_from_serialized_inputs");
    }
    verifyStringDatasetChunked(file.get(), "/stochastic_state/module_names", stochastic_names, w, result);

    const GravityForceCachePersistentView empty_cache{};
    const GravityForceCachePersistentView owning_cache_view =
        expected.gravity_force_cache != nullptr
            ? makeGravityForceCachePersistentView(*expected.gravity_force_cache)
            : empty_cache;
    const GravityForceCachePersistentView& cache =
        expected.gravity_force_cache_view != nullptr ? *expected.gravity_force_cache_view : owning_cache_view;
    if ((readScalarU32Attribute(openRequiredGroup(file.get(), "/gravity_force_cache").get(), "valid") != 0U) != cache.valid) {
      throw std::runtime_error("restart verification gravity force-cache valid mismatch");
    }
    verifyDatasetSpanChunked<std::uint64_t>(file.get(), "/gravity_force_cache/particle_id", cache.particle_id, H5T_NATIVE_UINT64, w, result);
    verifyDatasetSpanChunked<std::uint64_t>(file.get(), "/gravity_force_cache/gas_cell_id", cache.gas_cell_id, H5T_NATIVE_UINT64, w, result);
    verifyDatasetSpanChunked<double>(file.get(), "/gravity_force_cache/particle_accel_x_comoving", cache.particle_accel_x_comoving, H5T_NATIVE_DOUBLE, w, result);
    verifyDatasetSpanChunked<double>(file.get(), "/gravity_force_cache/particle_accel_y_comoving", cache.particle_accel_y_comoving, H5T_NATIVE_DOUBLE, w, result);
    verifyDatasetSpanChunked<double>(file.get(), "/gravity_force_cache/particle_accel_z_comoving", cache.particle_accel_z_comoving, H5T_NATIVE_DOUBLE, w, result);
    verifyDatasetSpanChunked<double>(file.get(), "/gravity_force_cache/cell_accel_x_comoving", cache.cell_accel_x_comoving, H5T_NATIVE_DOUBLE, w, result);
    verifyDatasetSpanChunked<double>(file.get(), "/gravity_force_cache/cell_accel_y_comoving", cache.cell_accel_y_comoving, H5T_NATIVE_DOUBLE, w, result);
    verifyDatasetSpanChunked<double>(file.get(), "/gravity_force_cache/cell_accel_z_comoving", cache.cell_accel_z_comoving, H5T_NATIVE_DOUBLE, w, result);
    {
      const std::string endpoint(k_hierarchical_endpoint_group);
      Hdf5Handle endpoint_group(openRequiredGroup(file.get(), endpoint));
      require_u32_attr(endpoint_group.get(), "valid", cache.hierarchical_endpoint_valid ? 1U : 0U, endpoint + "/@valid");
      require_u64_attr(endpoint_group.get(), "pm_field_version", cache.hierarchical_pm_field_version, endpoint + "/@pm_field_version");
      verifyDatasetSpanChunked<double>(file.get(), endpoint + "/tree_accel_x_comoving", cache.hierarchical_tree_accel_x_comoving, H5T_NATIVE_DOUBLE, w, result);
      verifyDatasetSpanChunked<double>(file.get(), endpoint + "/tree_accel_y_comoving", cache.hierarchical_tree_accel_y_comoving, H5T_NATIVE_DOUBLE, w, result);
      verifyDatasetSpanChunked<double>(file.get(), endpoint + "/tree_accel_z_comoving", cache.hierarchical_tree_accel_z_comoving, H5T_NATIVE_DOUBLE, w, result);
      verifyDatasetSpanChunked<double>(file.get(), endpoint + "/pm_accel_x_comoving", cache.hierarchical_pm_accel_x_comoving, H5T_NATIVE_DOUBLE, w, result);
      verifyDatasetSpanChunked<double>(file.get(), endpoint + "/pm_accel_y_comoving", cache.hierarchical_pm_accel_y_comoving, H5T_NATIVE_DOUBLE, w, result);
      verifyDatasetSpanChunked<double>(file.get(), endpoint + "/pm_accel_z_comoving", cache.hierarchical_pm_accel_z_comoving, H5T_NATIVE_DOUBLE, w, result);
      verifyDatasetSpanChunked<double>(file.get(), endpoint + "/nearest_pair_r2_comoving", cache.hierarchical_nearest_pair_r2_comoving, H5T_NATIVE_DOUBLE, w, result);
    }

    Hdf5Handle scheduler_group(openRequiredGroup(file.get(), "/scheduler"));
    require_u64_attr(scheduler_group.get(), "current_tick", expected.scheduler->currentTick(), "/scheduler/@current_tick");
    require_u32_attr(scheduler_group.get(), "max_bin", expected.scheduler->maxBin(), "/scheduler/@max_bin");
    verifyDataset1dChunked<std::uint8_t>(file.get(), "/scheduler/bin_index", expected.scheduler->elementCount(),
        [&](std::size_t i){ return expected.scheduler->binIndex(static_cast<std::uint32_t>(i)); }, H5T_NATIVE_UINT8, w, result);
    verifyDataset1dChunked<std::uint64_t>(file.get(), "/scheduler/next_activation_tick", expected.scheduler->elementCount(),
        [&](std::size_t i){ return expected.scheduler->nextActivationTick(static_cast<std::uint32_t>(i)); }, H5T_NATIVE_UINT64, w, result);
    verifyDataset1dChunked<std::uint8_t>(file.get(), "/scheduler/active_flag", expected.scheduler->elementCount(),
        [&](std::size_t i){ return static_cast<std::uint8_t>(expected.scheduler->isElementActive(static_cast<std::uint32_t>(i)) ? 1U : 0U); }, H5T_NATIVE_UINT8, w, result);
    verifyDataset1dChunked<std::uint8_t>(file.get(), "/scheduler/pending_bin_index", expected.scheduler->elementCount(),
        [&](std::size_t i){ return expected.scheduler->pendingBinIndex(static_cast<std::uint32_t>(i)); }, H5T_NATIVE_UINT8, w, result);
    if (expected.gas_cell_scheduler != nullptr) {
      Hdf5Handle gas_scheduler_group(openRequiredGroup(file.get(), "/gas_cell_scheduler"));
      require_string_attr(gas_scheduler_group.get(), "identity_key", std::string(k_gas_cell_scheduler_identity_key), "/gas_cell_scheduler/@identity_key");
      require_u64_attr(gas_scheduler_group.get(), "current_tick", expected.gas_cell_scheduler->currentTick(), "/gas_cell_scheduler/@current_tick");
      require_u32_attr(gas_scheduler_group.get(), "max_bin", expected.gas_cell_scheduler->maxBin(), "/gas_cell_scheduler/@max_bin");
      verifyDatasetSpanChunked<std::uint64_t>(file.get(), "/gas_cell_scheduler/gas_cell_id", state.gas_cells.gas_cell_id, H5T_NATIVE_UINT64, w, result);
      verifyDataset1dChunked<std::uint8_t>(file.get(), "/gas_cell_scheduler/bin_index", expected.gas_cell_scheduler->elementCount(),
          [&](std::size_t i){ return expected.gas_cell_scheduler->binIndex(static_cast<std::uint32_t>(i)); }, H5T_NATIVE_UINT8, w, result);
      verifyDataset1dChunked<std::uint64_t>(file.get(), "/gas_cell_scheduler/next_activation_tick", expected.gas_cell_scheduler->elementCount(),
          [&](std::size_t i){ return expected.gas_cell_scheduler->nextActivationTick(static_cast<std::uint32_t>(i)); }, H5T_NATIVE_UINT64, w, result);
      verifyDataset1dChunked<std::uint8_t>(file.get(), "/gas_cell_scheduler/active_flag", expected.gas_cell_scheduler->elementCount(),
          [&](std::size_t i){ return static_cast<std::uint8_t>(expected.gas_cell_scheduler->isElementActive(static_cast<std::uint32_t>(i)) ? 1U : 0U); }, H5T_NATIVE_UINT8, w, result);
      verifyDataset1dChunked<std::uint8_t>(file.get(), "/gas_cell_scheduler/pending_bin_index", expected.gas_cell_scheduler->elementCount(),
          [&](std::size_t i){ return expected.gas_cell_scheduler->pendingBinIndex(static_cast<std::uint32_t>(i)); }, H5T_NATIVE_UINT8, w, result);
    }
    verifyDistributedRestartStateChunked(file.get(), expected.distributed_gravity_state, w, result);

    const RestartIntegrityDigests expected_integrity = restartPayloadIntegrityDigestsImpl(
        expected, true, true, true, true, true, true, true, true, policy.memory_governor);
    result.legacy_payload_hash = readScalarU64Attribute(file.get(), "payload_integrity_hash");
    result.payload_sha256_hex = readScalarStringAttribute(file.get(), "payload_integrity_sha256_hex");
    if (result.legacy_payload_hash != expected_integrity.legacy_fnv1a ||
        result.payload_sha256_hex != expected_integrity.sha256_hex ||
        readScalarStringAttribute(file.get(), "payload_integrity_algorithm") != "sha256-canonical-le-v1") {
      result.failed_field = "/@payload_integrity";
      throw std::runtime_error("restart verification integrity metadata mismatch");
    }
    result.diagnostics = readRestartDiagnosticsGroup(file.get(), result.schema_version);
    std::uint64_t scheduler_active_count = 0U;
    std::uint64_t scheduler_pending_count = 0U;
    for (std::size_t i = 0; i < expected.scheduler->elementCount(); ++i) {
      scheduler_active_count += expected.scheduler->isElementActive(static_cast<std::uint32_t>(i)) ? 1ULL : 0ULL;
      scheduler_pending_count += expected.scheduler->pendingBinIndex(static_cast<std::uint32_t>(i)) != core::HierarchicalTimeBinScheduler::k_unset_pending_bin ? 1ULL : 0ULL;
    }
    if (result.diagnostics.restart_schema_name != restartSchema().name ||
        result.diagnostics.restart_schema_version != restartSchema().version ||
        result.diagnostics.step_index != expected.integrator_state->step_index ||
        result.diagnostics.scheduler_current_tick != expected.scheduler->currentTick() ||
        result.diagnostics.scheduler_max_bin != expected.scheduler->maxBin() ||
        result.diagnostics.scheduler_element_count != expected.scheduler->elementCount() ||
        result.diagnostics.scheduler_active_count != scheduler_active_count ||
        result.diagnostics.scheduler_pending_transition_count != scheduler_pending_count ||
        result.diagnostics.output_enabled != expected.output_cadence_state.output_enabled ||
        result.diagnostics.output_snapshot_due != expected.output_cadence_state.snapshot_due ||
        result.diagnostics.output_checkpoint_due != expected.output_cadence_state.checkpoint_due ||
        result.diagnostics.output_last_completed_step_index != expected.output_cadence_state.last_completed_step_index ||
        result.diagnostics.output_next_snapshot_step_index != expected.output_cadence_state.next_snapshot_step_index ||
        std::bit_cast<std::uint64_t>(result.diagnostics.output_next_snapshot_time_code) !=
            std::bit_cast<std::uint64_t>(expected.output_cadence_state.next_snapshot_time_code) ||
        result.diagnostics.stochastic_module_count != expected.stochastic_state.modules.size()) {
      result.failed_field = "/restart_diagnostics";
      throw std::runtime_error("restart verification diagnostics summary mismatch");
    }
    result.ok = true;
    result.detail = "streaming field-exact verification complete";
    workspace_reservation.release();
    return result;
  } catch (...) {
    workspace_reservation.release();
    throw;
  }
#endif
}

RestartReadResult readRestartCheckpointHdf5(
    const std::filesystem::path& input_path,
    const RestartReadPolicy& policy) {
#if !COSMOSIM_ENABLE_HDF5
  (void)input_path;
  (void)policy;
  throw std::runtime_error("restart checkpoint requires COSMOSIM_ENABLE_HDF5=ON");
#else
  Hdf5Handle file(H5Fopen(input_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed to open restart checkpoint: " + input_path.string());
  }

  validateFileKindAttribute(
      file.get(),
      sharedIoContractNames().restart_checkpoint_file_kind,
      "restart reader");
  const std::string schema_name = readScalarStringAttribute(file.get(), "restart_schema_name");
  const std::uint32_t schema_version = readScalarU32Attribute(file.get(), "restart_schema_version");
  const bool current_schema =
      schema_name == restartSchema().name && schema_version == restartSchema().version;
  const bool legacy_v14_schema =
      schema_name == k_restart_schema_name_v14 && schema_version == k_restart_schema_v14;
  const bool legacy_v15_schema =
      schema_name == k_restart_schema_name_v15 && schema_version == k_restart_schema_v15;
  const bool legacy_v16_schema =
      schema_name == k_restart_schema_name_v16 && schema_version == k_restart_schema_v16;
  const bool legacy_v17_schema =
      schema_name == k_restart_schema_name_v17 && schema_version == k_restart_schema_v17;
  const bool legacy_v18_schema =
      schema_name == k_restart_schema_name_v18 && schema_version == k_restart_schema_v18;
  const bool legacy_v19_schema =
      schema_name == k_restart_schema_name_v19 && schema_version == k_restart_schema_v19;
  const bool legacy_v20_schema =
      schema_name == k_restart_schema_name_v20 && schema_version == k_restart_schema_v20;
  const bool legacy_v21_schema =
      schema_name == k_restart_schema_name_v21 && schema_version == k_restart_schema_v21;
  const bool legacy_v23_schema =
      schema_name == k_restart_schema_name_v23 && schema_version == k_restart_schema_v23;
  if ((!current_schema && !legacy_v14_schema && !legacy_v15_schema && !legacy_v16_schema &&
       !legacy_v17_schema && !legacy_v18_schema && !legacy_v19_schema && !legacy_v20_schema &&
       !legacy_v21_schema && !legacy_v23_schema) ||
      !isRestartSchemaCompatible(schema_version)) {
    throw std::runtime_error(
        "restart schema is not compatible: file='" + schema_name + "' v" + std::to_string(schema_version) +
        ", expected='" + restartSchema().name + "' v" + std::to_string(restartSchema().version));
  }
  validateRestartCheckpointSchema(file.get(), schema_version);

  // Gate-5 real readback transaction: preflight HDF5 dimensions BEFORE any
  // population allocation, reject overflow/implausible extents, then admit the
  // ACTUAL candidate storage via policy.memory_governor. The live simulation
  // already occupies committed memory; the candidate is incremental and
  // coexists with it. The reader is the physical owner of this reservation
  // (workflow does NOT double-reserve). Candidate vectors below read directly
  // into final storage (no temp-then-copy) and verification runs before the
  // candidate could ever become canonical (caller decides publication).
  const RestartReadCandidateDimensions candidate_dims =
      preflightRestartCandidateDimensions(file.get(), policy);
  const std::uint64_t candidate_bytes = restartReadCandidateStagingBytes(candidate_dims);
  core::MemoryReservation readback_reservation;
  if (policy.memory_governor != nullptr && candidate_bytes > 0U) {
    readback_reservation = policy.memory_governor->reserve(
        core::MemoryClass::kPhaseResident, candidate_bytes,
        "restart.readback.candidate");
    readback_reservation.commit();
  }

  RestartReadResult result;
  result.normalized_config_hash_hex = readScalarStringAttribute(file.get(), "normalized_config_hash_hex");
  result.payload_hash_hex = readScalarStringAttribute(file.get(), "payload_integrity_hash_hex");
  result.payload_hash = readScalarU64Attribute(file.get(), "payload_integrity_hash");
  if (schema_version >= k_restart_schema_v23) {
    result.payload_integrity_algorithm = readScalarStringAttribute(file.get(), "payload_integrity_algorithm");
    result.payload_integrity_sha256_hex = readScalarStringAttribute(file.get(), "payload_integrity_sha256_hex");
    if (result.payload_integrity_algorithm != "sha256-canonical-le-v1") {
      throw std::runtime_error("unsupported restart payload integrity algorithm: " + result.payload_integrity_algorithm);
    }
  } else {
    result.payload_integrity_algorithm = "fnv1a64-native-legacy";
  }

  const auto& shared_names = sharedIoContractNames();
  result.normalized_config_text =
      readStringDataset(file.get(), std::string(shared_names.normalized_config_text_dataset));
  result.provenance = core::deserializeProvenanceRecord(
      readStringDataset(file.get(), std::string(shared_names.provenance_record_dataset)));

  readStateGroup(file.get(), result.state, schema_version);
  validateHydroGeometryStateForRestart(result.state, "restart reader");
  validateAmrTemporalBoundaryHistoryForRestart(result.state);

  Hdf5Handle integrator_group(H5Gopen2(file.get(), "/integrator", H5P_DEFAULT));
  result.integrator_state.current_time_code = readScalarF64Attribute(integrator_group.get(), "current_time_code");
  result.integrator_state.current_scale_factor = readScalarF64Attribute(integrator_group.get(), "current_scale_factor");
  result.integrator_state.current_redshift = readScalarF64Attribute(integrator_group.get(), "current_redshift");
  result.integrator_state.current_hubble_rate_code = readScalarF64Attribute(integrator_group.get(), "current_hubble_rate_code");
  result.integrator_state.time_si_per_code = readScalarF64Attribute(integrator_group.get(), "time_si_per_code");
  result.integrator_state.dt_time_code = readScalarF64Attribute(integrator_group.get(), "dt_time_code");
  result.integrator_state.last_drift_factor_code = readScalarF64Attribute(integrator_group.get(), "last_drift_factor_code");
  result.integrator_state.last_first_kick_factor_code = readScalarF64Attribute(integrator_group.get(), "last_first_kick_factor_code");
  result.integrator_state.last_second_kick_factor_code = readScalarF64Attribute(integrator_group.get(), "last_second_kick_factor_code");
  result.integrator_state.last_first_hubble_drag_factor = readScalarF64Attribute(integrator_group.get(), "last_first_hubble_drag_factor");
  result.integrator_state.last_second_hubble_drag_factor = readScalarF64Attribute(integrator_group.get(), "last_second_hubble_drag_factor");
  result.integrator_state.step_index = readScalarU64Attribute(integrator_group.get(), "step_index");
  result.integrator_state.scheme =
      static_cast<core::TimeStepScheme>(readScalarU32Attribute(integrator_group.get(), "scheme"));
  result.integrator_state.current_boundary_kind =
      static_cast<core::StepBoundaryKind>(readScalarU32Attribute(integrator_group.get(), "current_boundary_kind"));
  result.integrator_state.last_completed_boundary_kind =
      static_cast<core::StepBoundaryKind>(readScalarU32Attribute(integrator_group.get(), "last_completed_boundary_kind"));
  result.integrator_state.inside_kdk_step = readScalarU32Attribute(integrator_group.get(), "inside_kdk_step") != 0U;
  result.integrator_state.last_completed_restart_safe =
      readScalarU32Attribute(integrator_group.get(), "last_completed_restart_safe") != 0U;
  result.integrator_state.time_bins.hierarchical_enabled =
      readScalarU32Attribute(integrator_group.get(), "time_bins_hierarchical") != 0;
  result.integrator_state.time_bins.active_bin =
      static_cast<std::uint8_t>(readScalarU32Attribute(integrator_group.get(), "time_bins_active_bin"));
  result.integrator_state.time_bins.max_bin =
      static_cast<std::uint8_t>(readScalarU32Attribute(integrator_group.get(), "time_bins_max_bin"));
  result.integrator_state.pm_long_range_field_valid =
      readScalarU32Attribute(integrator_group.get(), "pm_long_range_field_valid") != 0U;
  result.integrator_state.pm_refresh_enabled =
      schema_version >= k_restart_schema_v20 &&
      readScalarU32Attribute(integrator_group.get(), "pm_refresh_enabled") != 0U;
  core::PmSynchronizationPersistentState pm_sync_state;
  pm_sync_state.cadence_steps = readScalarU64Attribute(integrator_group.get(), "pm_cadence_steps");
  pm_sync_state.gravity_kick_opportunity = readScalarU64Attribute(integrator_group.get(), "pm_gravity_kick_opportunity");
  pm_sync_state.last_refresh_opportunity = readScalarU64Attribute(integrator_group.get(), "pm_last_refresh_opportunity");
  pm_sync_state.field_version = readScalarU64Attribute(integrator_group.get(), "pm_field_version");
  pm_sync_state.last_refresh_step_index = readScalarU64Attribute(integrator_group.get(), "pm_last_refresh_step_index");
  pm_sync_state.last_refresh_scale_factor = readScalarF64Attribute(integrator_group.get(), "pm_last_refresh_scale_factor");
  pm_sync_state.refresh_commit_pending = readScalarU32Attribute(integrator_group.get(), "pm_refresh_commit_pending") != 0U;
  pm_sync_state.pending_refresh_opportunity = readScalarU64Attribute(integrator_group.get(), "pm_pending_refresh_opportunity");
  pm_sync_state.pending_refresh_field_version = readScalarU64Attribute(integrator_group.get(), "pm_pending_refresh_field_version");
  result.integrator_state.pm_sync_state.importPersistentState(pm_sync_state);
  if (schema_version >= k_restart_schema_v20) {
    Hdf5Handle force_cache_group(H5Gopen2(file.get(), "/gravity_force_cache", H5P_DEFAULT));
    result.gravity_force_cache.valid = readScalarU32Attribute(force_cache_group.get(), "valid") != 0U;
    result.gravity_force_cache.particle_id = readDataset1d<std::uint64_t>(force_cache_group.get(), "particle_id", H5T_NATIVE_UINT64);
    result.gravity_force_cache.gas_cell_id = readDataset1d<std::uint64_t>(force_cache_group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
    result.gravity_force_cache.particle_accel_x_comoving = readDataset1d<double>(force_cache_group.get(), "particle_accel_x_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.particle_accel_y_comoving = readDataset1d<double>(force_cache_group.get(), "particle_accel_y_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.particle_accel_z_comoving = readDataset1d<double>(force_cache_group.get(), "particle_accel_z_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.cell_accel_x_comoving = readDataset1d<double>(force_cache_group.get(), "cell_accel_x_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.cell_accel_y_comoving = readDataset1d<double>(force_cache_group.get(), "cell_accel_y_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.cell_accel_z_comoving = readDataset1d<double>(force_cache_group.get(), "cell_accel_z_comoving", H5T_NATIVE_DOUBLE);
    if (schema_version >= k_restart_schema_v24) {
      Hdf5Handle endpoint_group(H5Gopen2(file.get(), k_hierarchical_endpoint_group.data(), H5P_DEFAULT));
      auto& cache = result.gravity_force_cache;
      cache.hierarchical_endpoint_valid = readScalarU32Attribute(endpoint_group.get(), "valid") != 0U;
      cache.hierarchical_pm_field_version = readScalarU64Attribute(endpoint_group.get(), "pm_field_version");
      cache.hierarchical_tree_accel_x_comoving = readDataset1d<double>(endpoint_group.get(), "tree_accel_x_comoving", H5T_NATIVE_DOUBLE);
      cache.hierarchical_tree_accel_y_comoving = readDataset1d<double>(endpoint_group.get(), "tree_accel_y_comoving", H5T_NATIVE_DOUBLE);
      cache.hierarchical_tree_accel_z_comoving = readDataset1d<double>(endpoint_group.get(), "tree_accel_z_comoving", H5T_NATIVE_DOUBLE);
      cache.hierarchical_pm_accel_x_comoving = readDataset1d<double>(endpoint_group.get(), "pm_accel_x_comoving", H5T_NATIVE_DOUBLE);
      cache.hierarchical_pm_accel_y_comoving = readDataset1d<double>(endpoint_group.get(), "pm_accel_y_comoving", H5T_NATIVE_DOUBLE);
      cache.hierarchical_pm_accel_z_comoving = readDataset1d<double>(endpoint_group.get(), "pm_accel_z_comoving", H5T_NATIVE_DOUBLE);
      cache.hierarchical_nearest_pair_r2_comoving = readDataset1d<double>(endpoint_group.get(), "nearest_pair_r2_comoving", H5T_NATIVE_DOUBLE);
    }
    validateGravityForceCacheForRestart(
        makeGravityForceCachePersistentView(result.gravity_force_cache), result.state, "restart reader");
  }

  Hdf5Handle scheduler_group(H5Gopen2(file.get(), "/scheduler", H5P_DEFAULT));
  result.scheduler_state.current_tick = readScalarU64Attribute(scheduler_group.get(), "current_tick");
  result.scheduler_state.max_bin = static_cast<std::uint8_t>(readScalarU32Attribute(scheduler_group.get(), "max_bin"));
  result.scheduler_state.bin_index = readDataset1d<std::uint8_t>(scheduler_group.get(), "bin_index", H5T_NATIVE_UINT8);
  result.scheduler_state.next_activation_tick =
      readDataset1d<std::uint64_t>(scheduler_group.get(), "next_activation_tick", H5T_NATIVE_UINT64);
  result.scheduler_state.active_flag = readDataset1d<std::uint8_t>(scheduler_group.get(), "active_flag", H5T_NATIVE_UINT8);
  result.scheduler_state.pending_bin_index =
      readDataset1d<std::uint8_t>(scheduler_group.get(), "pending_bin_index", H5T_NATIVE_UINT8);
  validateRestartTimeBinMirrorsAgainstScheduler(result.state, result.scheduler_state, "restart reader");
  rebuildRestartTimeBinMirrorsFromScheduler(result.state, result.scheduler_state);
  if (schema_version >= k_restart_schema_v19) {
    Hdf5Handle gas_scheduler_group(H5Gopen2(file.get(), "/gas_cell_scheduler", H5P_DEFAULT));
    const std::string identity_key = readScalarStringAttribute(gas_scheduler_group.get(), "identity_key");
    if (identity_key != k_gas_cell_scheduler_identity_key) {
      throw std::runtime_error("restart reader: /gas_cell_scheduler identity_key must be gas_cell_id");
    }
    result.gas_cell_scheduler_state.current_tick =
        readScalarU64Attribute(gas_scheduler_group.get(), "current_tick");
    result.gas_cell_scheduler_state.max_bin =
        static_cast<std::uint8_t>(readScalarU32Attribute(gas_scheduler_group.get(), "max_bin"));
    result.gas_cell_scheduler_ids =
        readDataset1d<std::uint64_t>(gas_scheduler_group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
    result.gas_cell_scheduler_state.bin_index =
        readDataset1d<std::uint8_t>(gas_scheduler_group.get(), "bin_index", H5T_NATIVE_UINT8);
    result.gas_cell_scheduler_state.next_activation_tick =
        readDataset1d<std::uint64_t>(gas_scheduler_group.get(), "next_activation_tick", H5T_NATIVE_UINT64);
    result.gas_cell_scheduler_state.active_flag =
        readDataset1d<std::uint8_t>(gas_scheduler_group.get(), "active_flag", H5T_NATIVE_UINT8);
    result.gas_cell_scheduler_state.pending_bin_index =
        readDataset1d<std::uint8_t>(gas_scheduler_group.get(), "pending_bin_index", H5T_NATIVE_UINT8);
    if (result.gas_cell_scheduler_state.current_tick != result.scheduler_state.current_tick) {
      throw std::runtime_error("restart reader: particle and gas-cell schedulers have different current_tick values");
    }
    validateGasCellTimeBinMirrorsAgainstScheduler(
        result.state,
        result.gas_cell_scheduler_state,
        result.gas_cell_scheduler_ids,
        "restart reader");
    rebuildGasCellTimeBinMirrorsFromScheduler(
        result.state, result.gas_cell_scheduler_state, result.gas_cell_scheduler_ids);
  } else {
    result.gas_cell_scheduler_state = legacyGasCellSchedulerStateFromMirrors(
        result.state, result.scheduler_state);
    result.gas_cell_scheduler_ids = gasCellIdsByDenseLocalRow(result.state, "legacy restart reader");
  }
  result.output_cadence_state = readOutputCadenceGroup(file.get(), schema_version);
  validateOutputCadenceStateForRestart(result.output_cadence_state, result.integrator_state, "restart reader");
  result.stochastic_state = readStochasticStateGroup(file.get());
  validateStochasticStateForRestart(result.stochastic_state, result.integrator_state, "restart reader");
  result.diagnostics = readRestartDiagnosticsGroup(file.get(), schema_version);
  if (result.diagnostics.restart_schema_name != schema_name ||
      result.diagnostics.restart_schema_version != schema_version ||
      result.diagnostics.step_index != result.integrator_state.step_index ||
      result.diagnostics.scheduler_current_tick != result.scheduler_state.current_tick ||
      result.diagnostics.scheduler_element_count != result.scheduler_state.bin_index.size() ||
      (schema_version >= k_restart_schema_v19 &&
       (result.diagnostics.gas_cell_scheduler_current_tick != result.gas_cell_scheduler_state.current_tick ||
        result.diagnostics.gas_cell_scheduler_element_count != result.gas_cell_scheduler_state.bin_index.size())) ||
      result.diagnostics.output_last_completed_step_index != result.output_cadence_state.last_completed_step_index ||
      (schema_version >= k_restart_schema_v21 &&
       result.diagnostics.output_next_snapshot_time_code !=
           result.output_cadence_state.next_snapshot_time_code) ||
      result.diagnostics.stochastic_module_count != result.stochastic_state.modules.size()) {
    throw std::runtime_error("restart diagnostics summary is inconsistent with authoritative restart state");
  }
  readDistributedGravityGroup(file.get(), result.distributed_gravity_state);

  RestartWritePayload verify_payload;
  verify_payload.persistent_state.simulation_state = &result.state;
  verify_payload.integrator_state = &result.integrator_state;
  core::HierarchicalTimeBinScheduler verify_scheduler(result.scheduler_state.max_bin);
  verify_scheduler.importPersistentState(result.scheduler_state);
  core::HierarchicalTimeBinScheduler verify_gas_cell_scheduler(result.gas_cell_scheduler_state.max_bin);
  verify_gas_cell_scheduler.importPersistentState(result.gas_cell_scheduler_state);
  verify_payload.scheduler = &verify_scheduler;
  // RestartWritePayload borrows this span-only view during integrity traversal;
  // keep the view object alive until restartPayloadIntegrityDigestsImpl returns.
  std::optional<GravityForceCachePersistentView> verify_force_cache_view;
  if (schema_version >= k_restart_schema_v20) {
    verify_force_cache_view.emplace(makeGravityForceCachePersistentView(result.gravity_force_cache));
    verify_payload.gravity_force_cache_view = &*verify_force_cache_view;
  }
  if (schema_version >= k_restart_schema_v19) {
    verify_payload.gas_cell_scheduler = &verify_gas_cell_scheduler;
  }
  verify_payload.normalized_config_hash_hex = result.normalized_config_hash_hex;
  verify_payload.normalized_config_text = result.normalized_config_text;
  verify_payload.provenance = result.provenance;
  verify_payload.distributed_gravity_state = result.distributed_gravity_state;
  verify_payload.output_cadence_state = result.output_cadence_state;
  verify_payload.stochastic_state = result.stochastic_state;

  const RestartIntegrityDigests computed_integrity =
      restartPayloadIntegrityDigestsImpl(
          verify_payload,
          schema_version >= k_restart_schema_v15,
          schema_version >= k_restart_schema_v17,
          schema_version >= k_restart_schema_v18,
          schema_version >= k_restart_schema_v19,
          schema_version >= k_restart_schema_v20,
          schema_version >= k_restart_schema_v21,
          schema_version >= k_restart_schema_v22,
          schema_version >= k_restart_schema_v24);
  if (computed_integrity.legacy_fnv1a != result.payload_hash ||
      hexU64(computed_integrity.legacy_fnv1a) != result.payload_hash_hex) {
    throw std::runtime_error("restart legacy payload integrity hash mismatch");
  }
  if (schema_version >= k_restart_schema_v23 &&
      computed_integrity.sha256_hex != result.payload_integrity_sha256_hex) {
    throw std::runtime_error("restart canonical SHA-256 payload integrity mismatch");
  }

  result.payload_integrity_verified = true;

  // The checkpoint stores the original flux-register geometry epoch for its
  // integrity digest. The reconstructed SimulationState has a *new* in-memory
  // dense-cell-index epoch; the original epoch is not authoritative in this
  // process. After SHA-256/FNV verification, validate every stable target and
  // rebase only the transient geometry marker. Never rebase before the digest
  // or a corrupted checkpoint could be accepted as an ordinary restart.
  for (core::PendingFluxRegisterRecord& pending :
       result.state.pending_flux_registers.mutableRecords()) {
    const auto row = result.state.rowForGasCellId(pending.coarse_gas_cell_id);
    const auto owning_patch =
        result.state.owningPatchIdForGasCellId(pending.coarse_gas_cell_id);
    if (!row.has_value() || *row != pending.coarse_cell_index ||
        !owning_patch.has_value() || *owning_patch != pending.coarse_patch_id ||
        pending.gas_cell_identity_generation != result.state.gasCellIdentityGeneration()) {
      throw std::runtime_error(
          "restart pending AMR flux register has an invalid stable gas/patch target or identity generation");
    }
    pending.patch_geometry_generation = result.state.cellIndexGeneration();
  }

  // The on-disk schema remains fully explicit. Only after complete logical
  // validation and integrity verification do we release redundant homogeneous
  // DMO metadata capacities. Nonuniform imported SFC or other metadata remains
  // materialized exactly because compactHomogeneousDmoMetadata rejects it or
  // preserves the exact SFC exception lane.
  if (result.state.particles.size() != 0U &&
      result.state.particle_sidecar.owning_rank.size() == result.state.particles.size()) {
    const std::uint32_t restored_local_rank = result.state.particle_sidecar.owning_rank.front();
    (void)result.state.compactHomogeneousDmoMetadata(restored_local_rank);
  }

  return result;
#endif
}

}  // namespace cosmosim::io
