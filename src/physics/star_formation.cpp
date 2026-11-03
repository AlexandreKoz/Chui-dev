#include "cosmosim/physics/star_formation.hpp"

#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/stellar_population.hpp"

#include <algorithm>
#include <optional>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cosmosim::physics {
namespace {

constexpr double kDensityFloor = 1.0e-30;
constexpr double kMassFloor = 1.0e-30;
constexpr double kTimeFloor = 1.0e-30;
constexpr std::uint64_t kGeneratedParticleIdTag = 0x8000000000000000ull;
constexpr std::uint64_t kParticleIdDomain = 0x31d9a7c5f47b2e11ull;
constexpr double kU01Scale = 0x1.0p-53;

[[nodiscard]] std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ull;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ull;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebull;
  return value ^ (value >> 31U);
}

[[nodiscard]] bool finitePositive(double value) {
  return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] bool finiteNonNegative(double value) {
  return std::isfinite(value) && value >= 0.0;
}

[[nodiscard]] double safeMetallicity(const StarFormationCellInput& cell) {
  if (finitePositive(cell.gas_mass_code) && finiteNonNegative(cell.gas_metal_mass_code)) {
    return std::clamp(cell.gas_metal_mass_code / cell.gas_mass_code, 0.0, 1.0);
  }
  return std::clamp(
      std::isfinite(cell.metallicity_mass_fraction) ? cell.metallicity_mass_fraction : 0.0,
      0.0,
      1.0);
}

void incrementRejectionCounter(
    StarFormationCounters& counters,
    StarFormationRejectionReason reason) {
  switch (reason) {
    case StarFormationRejectionReason::kNone:
      break;
    case StarFormationRejectionReason::kInactive:
      ++counters.rejected_inactive;
      break;
    case StarFormationRejectionReason::kNotOwned:
      ++counters.rejected_not_owned;
      break;
    case StarFormationRejectionReason::kNonLeaf:
      ++counters.rejected_non_leaf;
      break;
    case StarFormationRejectionReason::kGhost:
      ++counters.rejected_ghost;
      break;
    case StarFormationRejectionReason::kInvalidState:
      ++counters.rejected_invalid_state;
      break;
    case StarFormationRejectionReason::kNonPositiveMass:
      ++counters.rejected_non_positive_mass;
      break;
    case StarFormationRejectionReason::kNonPositiveVolume:
      ++counters.rejected_non_positive_volume;
      break;
    case StarFormationRejectionReason::kNotConverging:
      ++counters.rejected_not_converging;
      break;
    case StarFormationRejectionReason::kUnbound:
      ++counters.rejected_unbound;
      break;
    case StarFormationRejectionReason::kJeansStable:
      ++counters.rejected_jeans_stable;
      break;
    case StarFormationRejectionReason::kLegacyDensity:
      ++counters.rejected_legacy_density_threshold;
      break;
    case StarFormationRejectionReason::kLegacyTemperature:
      ++counters.rejected_legacy_temperature_threshold;
      break;
    case StarFormationRejectionReason::kAdaptiveMassFloor:
      ++counters.rejected_adaptive_mass_floor;
      break;
    case StarFormationRejectionReason::kEffectiveBelowDensityThreshold:
      ++counters.rejected_effective_below_density_threshold;
      break;
    case StarFormationRejectionReason::kEffectiveBelowOverdensityThreshold:
      ++counters.rejected_effective_below_overdensity_threshold;
      break;
    case StarFormationRejectionReason::kEffectiveHotAboveEos:
      ++counters.rejected_effective_hot_above_eos;
      break;
    case StarFormationRejectionReason::kEffectiveInvalidEquilibrium:
      ++counters.rejected_effective_invalid_equilibrium;
      break;
  }
}


void accumulateStarFormationCounters(StarFormationCounters& cumulative, const StarFormationCounters& latest) {
  cumulative.scanned_cells += latest.scanned_cells;
  cumulative.rejected_inactive += latest.rejected_inactive;
  cumulative.rejected_not_owned += latest.rejected_not_owned;
  cumulative.rejected_non_leaf += latest.rejected_non_leaf;
  cumulative.rejected_ghost += latest.rejected_ghost;
  cumulative.rejected_invalid_state += latest.rejected_invalid_state;
  cumulative.rejected_non_positive_mass += latest.rejected_non_positive_mass;
  cumulative.rejected_non_positive_volume += latest.rejected_non_positive_volume;
  cumulative.rejected_not_converging += latest.rejected_not_converging;
  cumulative.rejected_unbound += latest.rejected_unbound;
  cumulative.rejected_jeans_stable += latest.rejected_jeans_stable;
  cumulative.rejected_legacy_density_threshold += latest.rejected_legacy_density_threshold;
  cumulative.rejected_legacy_temperature_threshold += latest.rejected_legacy_temperature_threshold;
  cumulative.rejected_adaptive_mass_floor += latest.rejected_adaptive_mass_floor;
  cumulative.rejected_effective_below_density_threshold += latest.rejected_effective_below_density_threshold;
  cumulative.rejected_effective_below_overdensity_threshold += latest.rejected_effective_below_overdensity_threshold;
  cumulative.rejected_effective_hot_above_eos += latest.rejected_effective_hot_above_eos;
  cumulative.rejected_effective_invalid_equilibrium += latest.rejected_effective_invalid_equilibrium;
  cumulative.eligible_cells += latest.eligible_cells;
  cumulative.effective_cells_scanned += latest.effective_cells_scanned;
  cumulative.effective_cells_above_threshold += latest.effective_cells_above_threshold;
  cumulative.effective_cells_on_eos += latest.effective_cells_on_eos;
  cumulative.effective_cells_hot_above_eos += latest.effective_cells_hot_above_eos;
  cumulative.effective_cells_invalid_equilibrium += latest.effective_cells_invalid_equilibrium;
  cumulative.spawn_events += latest.spawn_events;
  cumulative.spawned_particles += latest.spawned_particles;
  cumulative.effective_cold_gas_mass_code += latest.effective_cold_gas_mass_code;
  cumulative.effective_instantaneous_sfr_code += latest.effective_instantaneous_sfr_code;
  cumulative.expected_spawn_mass_code += latest.expected_spawn_mass_code;
  cumulative.spawned_mass_code += latest.spawned_mass_code;
  cumulative.gas_mass_removed_code += latest.gas_mass_removed_code;
  cumulative.star_mass_created_code += latest.star_mass_created_code;
  cumulative.mass_residual_code += latest.mass_residual_code;
  cumulative.gas_momentum_removed_x_code += latest.gas_momentum_removed_x_code;
  cumulative.gas_momentum_removed_y_code += latest.gas_momentum_removed_y_code;
  cumulative.gas_momentum_removed_z_code += latest.gas_momentum_removed_z_code;
  cumulative.star_momentum_created_x_code += latest.star_momentum_created_x_code;
  cumulative.star_momentum_created_y_code += latest.star_momentum_created_y_code;
  cumulative.star_momentum_created_z_code += latest.star_momentum_created_z_code;
  cumulative.momentum_residual_norm_code += latest.momentum_residual_norm_code;
  cumulative.gas_metal_mass_removed_code += latest.gas_metal_mass_removed_code;
  cumulative.star_metal_mass_created_code += latest.star_metal_mass_created_code;
  cumulative.metal_mass_residual_code += latest.metal_mass_residual_code;
  cumulative.gas_internal_energy_removed_code += latest.gas_internal_energy_removed_code;
  cumulative.star_kinetic_energy_created_code += latest.star_kinetic_energy_created_code;
  cumulative.star_formation_internal_energy_sink_code += latest.star_formation_internal_energy_sink_code;
  if (latest.effective_min_cold_fraction > 0.0 && (cumulative.effective_min_cold_fraction == 0.0 || latest.effective_min_cold_fraction < cumulative.effective_min_cold_fraction)) cumulative.effective_min_cold_fraction = latest.effective_min_cold_fraction;
  if (latest.effective_min_pressure_code > 0.0 && (cumulative.effective_min_pressure_code == 0.0 || latest.effective_min_pressure_code < cumulative.effective_min_pressure_code)) cumulative.effective_min_pressure_code = latest.effective_min_pressure_code;
  if (latest.minimum_free_fall_time_code > 0.0 && (cumulative.minimum_free_fall_time_code == 0.0 || latest.minimum_free_fall_time_code < cumulative.minimum_free_fall_time_code)) cumulative.minimum_free_fall_time_code = latest.minimum_free_fall_time_code;
  if (latest.minimum_collapse_time_code > 0.0 && (cumulative.minimum_collapse_time_code == 0.0 || latest.minimum_collapse_time_code < cumulative.minimum_collapse_time_code)) cumulative.minimum_collapse_time_code = latest.minimum_collapse_time_code;
  cumulative.effective_max_cold_fraction = std::max(cumulative.effective_max_cold_fraction, latest.effective_max_cold_fraction);
  cumulative.effective_max_pressure_code = std::max(cumulative.effective_max_pressure_code, latest.effective_max_pressure_code);
  cumulative.maximum_spawn_probability = std::max(cumulative.maximum_spawn_probability, latest.maximum_spawn_probability);
  cumulative.maximum_fractional_mass_conversion = std::max(cumulative.maximum_fractional_mass_conversion, latest.maximum_fractional_mass_conversion);
}

void writeStarFormationCounters(
    std::ostringstream& stream,
    std::string_view prefix,
    const StarFormationCounters& counters) {
  stream << prefix << "scanned_cells=" << counters.scanned_cells << "\n";
  stream << prefix << "rejected_inactive=" << counters.rejected_inactive << "\n";
  stream << prefix << "rejected_not_owned=" << counters.rejected_not_owned << "\n";
  stream << prefix << "rejected_non_leaf=" << counters.rejected_non_leaf << "\n";
  stream << prefix << "rejected_ghost=" << counters.rejected_ghost << "\n";
  stream << prefix << "rejected_invalid_state=" << counters.rejected_invalid_state << "\n";
  stream << prefix << "rejected_non_positive_mass=" << counters.rejected_non_positive_mass << "\n";
  stream << prefix << "rejected_non_positive_volume=" << counters.rejected_non_positive_volume << "\n";
  stream << prefix << "rejected_not_converging=" << counters.rejected_not_converging << "\n";
  stream << prefix << "rejected_unbound=" << counters.rejected_unbound << "\n";
  stream << prefix << "rejected_jeans_stable=" << counters.rejected_jeans_stable << "\n";
  stream << prefix << "rejected_legacy_density_threshold=" << counters.rejected_legacy_density_threshold << "\n";
  stream << prefix << "rejected_legacy_temperature_threshold=" << counters.rejected_legacy_temperature_threshold << "\n";
  stream << prefix << "rejected_adaptive_mass_floor=" << counters.rejected_adaptive_mass_floor << "\n";
  stream << prefix << "rejected_effective_below_density_threshold=" << counters.rejected_effective_below_density_threshold << "\n";
  stream << prefix << "rejected_effective_below_overdensity_threshold=" << counters.rejected_effective_below_overdensity_threshold << "\n";
  stream << prefix << "rejected_effective_hot_above_eos=" << counters.rejected_effective_hot_above_eos << "\n";
  stream << prefix << "rejected_effective_invalid_equilibrium=" << counters.rejected_effective_invalid_equilibrium << "\n";
  stream << prefix << "eligible_cells=" << counters.eligible_cells << "\n";
  stream << prefix << "effective_cells_scanned=" << counters.effective_cells_scanned << "\n";
  stream << prefix << "effective_cells_above_threshold=" << counters.effective_cells_above_threshold << "\n";
  stream << prefix << "effective_cells_on_eos=" << counters.effective_cells_on_eos << "\n";
  stream << prefix << "effective_cells_hot_above_eos=" << counters.effective_cells_hot_above_eos << "\n";
  stream << prefix << "effective_cells_invalid_equilibrium=" << counters.effective_cells_invalid_equilibrium << "\n";
  stream << prefix << "spawn_events=" << counters.spawn_events << "\n";
  stream << prefix << "spawned_particles=" << counters.spawned_particles << "\n";
  stream << prefix << "effective_cold_gas_mass_code=" << counters.effective_cold_gas_mass_code << "\n";
  stream << prefix << "effective_instantaneous_sfr_code=" << counters.effective_instantaneous_sfr_code << "\n";
  stream << prefix << "expected_spawn_mass_code=" << counters.expected_spawn_mass_code << "\n";
  stream << prefix << "spawned_mass_code=" << counters.spawned_mass_code << "\n";
  stream << prefix << "gas_mass_removed_code=" << counters.gas_mass_removed_code << "\n";
  stream << prefix << "star_mass_created_code=" << counters.star_mass_created_code << "\n";
  stream << prefix << "mass_residual_code=" << counters.mass_residual_code << "\n";
  stream << prefix << "gas_momentum_removed_x_code=" << counters.gas_momentum_removed_x_code << "\n";
  stream << prefix << "gas_momentum_removed_y_code=" << counters.gas_momentum_removed_y_code << "\n";
  stream << prefix << "gas_momentum_removed_z_code=" << counters.gas_momentum_removed_z_code << "\n";
  stream << prefix << "star_momentum_created_x_code=" << counters.star_momentum_created_x_code << "\n";
  stream << prefix << "star_momentum_created_y_code=" << counters.star_momentum_created_y_code << "\n";
  stream << prefix << "star_momentum_created_z_code=" << counters.star_momentum_created_z_code << "\n";
  stream << prefix << "momentum_residual_norm_code=" << counters.momentum_residual_norm_code << "\n";
  stream << prefix << "gas_metal_mass_removed_code=" << counters.gas_metal_mass_removed_code << "\n";
  stream << prefix << "star_metal_mass_created_code=" << counters.star_metal_mass_created_code << "\n";
  stream << prefix << "metal_mass_residual_code=" << counters.metal_mass_residual_code << "\n";
  stream << prefix << "gas_internal_energy_removed_code=" << counters.gas_internal_energy_removed_code << "\n";
  stream << prefix << "star_kinetic_energy_created_code=" << counters.star_kinetic_energy_created_code << "\n";
  stream << prefix << "star_formation_internal_energy_sink_code=" << counters.star_formation_internal_energy_sink_code << "\n";
  stream << prefix << "effective_min_cold_fraction=" << counters.effective_min_cold_fraction << "\n";
  stream << prefix << "effective_min_pressure_code=" << counters.effective_min_pressure_code << "\n";
  stream << prefix << "minimum_free_fall_time_code=" << counters.minimum_free_fall_time_code << "\n";
  stream << prefix << "minimum_collapse_time_code=" << counters.minimum_collapse_time_code << "\n";
  stream << prefix << "effective_max_cold_fraction=" << counters.effective_max_cold_fraction << "\n";
  stream << prefix << "effective_max_pressure_code=" << counters.effective_max_pressure_code << "\n";
  stream << prefix << "maximum_spawn_probability=" << counters.maximum_spawn_probability << "\n";
  stream << prefix << "maximum_fractional_mass_conversion=" << counters.maximum_fractional_mass_conversion << "\n";
}

[[nodiscard]] StarFormationCounters cumulativeCountersFromState(const core::SimulationState& state) {
  StarFormationCounters counters;
  const core::ModuleSidecarBlock* block = state.sidecars.find("star_formation");
  if (block == nullptr || block->payload.empty()) return counters;
  std::string text;
  text.reserve(block->payload.size());
  for (const std::byte value : block->payload) text.push_back(static_cast<char>(value));
  std::unordered_map<std::string, std::string> values;
  std::istringstream stream(text);
  std::string line;
  while (std::getline(stream, line)) {
    const std::size_t split = line.find('=');
    if (split != std::string::npos) values.emplace(line.substr(0, split), line.substr(split + 1U));
  }
  const auto read_u64 = [&](std::string_view name, std::uint64_t& target) {
    const auto it = values.find("cumulative." + std::string(name));
    if (it != values.end()) target = static_cast<std::uint64_t>(std::stoull(it->second));
  };
  const auto read_double = [&](std::string_view name, double& target) {
    const auto it = values.find("cumulative." + std::string(name));
    if (it != values.end()) target = std::stod(it->second);
  };
  try {
    read_u64("scanned_cells", counters.scanned_cells);
    read_u64("rejected_inactive", counters.rejected_inactive);
    read_u64("rejected_not_owned", counters.rejected_not_owned);
    read_u64("rejected_non_leaf", counters.rejected_non_leaf);
    read_u64("rejected_ghost", counters.rejected_ghost);
    read_u64("rejected_invalid_state", counters.rejected_invalid_state);
    read_u64("rejected_non_positive_mass", counters.rejected_non_positive_mass);
    read_u64("rejected_non_positive_volume", counters.rejected_non_positive_volume);
    read_u64("rejected_not_converging", counters.rejected_not_converging);
    read_u64("rejected_unbound", counters.rejected_unbound);
    read_u64("rejected_jeans_stable", counters.rejected_jeans_stable);
    read_u64("rejected_legacy_density_threshold", counters.rejected_legacy_density_threshold);
    read_u64("rejected_legacy_temperature_threshold", counters.rejected_legacy_temperature_threshold);
    read_u64("rejected_adaptive_mass_floor", counters.rejected_adaptive_mass_floor);
    read_u64("rejected_effective_below_density_threshold", counters.rejected_effective_below_density_threshold);
    read_u64("rejected_effective_below_overdensity_threshold", counters.rejected_effective_below_overdensity_threshold);
    read_u64("rejected_effective_hot_above_eos", counters.rejected_effective_hot_above_eos);
    read_u64("rejected_effective_invalid_equilibrium", counters.rejected_effective_invalid_equilibrium);
    read_u64("eligible_cells", counters.eligible_cells);
    read_u64("effective_cells_scanned", counters.effective_cells_scanned);
    read_u64("effective_cells_above_threshold", counters.effective_cells_above_threshold);
    read_u64("effective_cells_on_eos", counters.effective_cells_on_eos);
    read_u64("effective_cells_hot_above_eos", counters.effective_cells_hot_above_eos);
    read_u64("effective_cells_invalid_equilibrium", counters.effective_cells_invalid_equilibrium);
    read_u64("spawn_events", counters.spawn_events);
    read_u64("spawned_particles", counters.spawned_particles);
    read_double("effective_cold_gas_mass_code", counters.effective_cold_gas_mass_code);
    read_double("effective_instantaneous_sfr_code", counters.effective_instantaneous_sfr_code);
    read_double("expected_spawn_mass_code", counters.expected_spawn_mass_code);
    read_double("spawned_mass_code", counters.spawned_mass_code);
    read_double("gas_mass_removed_code", counters.gas_mass_removed_code);
    read_double("star_mass_created_code", counters.star_mass_created_code);
    read_double("mass_residual_code", counters.mass_residual_code);
    read_double("gas_momentum_removed_x_code", counters.gas_momentum_removed_x_code);
    read_double("gas_momentum_removed_y_code", counters.gas_momentum_removed_y_code);
    read_double("gas_momentum_removed_z_code", counters.gas_momentum_removed_z_code);
    read_double("star_momentum_created_x_code", counters.star_momentum_created_x_code);
    read_double("star_momentum_created_y_code", counters.star_momentum_created_y_code);
    read_double("star_momentum_created_z_code", counters.star_momentum_created_z_code);
    read_double("momentum_residual_norm_code", counters.momentum_residual_norm_code);
    read_double("gas_metal_mass_removed_code", counters.gas_metal_mass_removed_code);
    read_double("star_metal_mass_created_code", counters.star_metal_mass_created_code);
    read_double("metal_mass_residual_code", counters.metal_mass_residual_code);
    read_double("gas_internal_energy_removed_code", counters.gas_internal_energy_removed_code);
    read_double("star_kinetic_energy_created_code", counters.star_kinetic_energy_created_code);
    read_double("star_formation_internal_energy_sink_code", counters.star_formation_internal_energy_sink_code);
    read_double("effective_min_cold_fraction", counters.effective_min_cold_fraction);
    read_double("effective_min_pressure_code", counters.effective_min_pressure_code);
    read_double("minimum_free_fall_time_code", counters.minimum_free_fall_time_code);
    read_double("minimum_collapse_time_code", counters.minimum_collapse_time_code);
    read_double("effective_max_cold_fraction", counters.effective_max_cold_fraction);
    read_double("effective_max_pressure_code", counters.effective_max_pressure_code);
    read_double("maximum_spawn_probability", counters.maximum_spawn_probability);
    read_double("maximum_fractional_mass_conversion", counters.maximum_fractional_mass_conversion);
  } catch (const std::exception&) {
    return StarFormationCounters{};
  }
  return counters;
}


struct StarBirthPlan {
  std::uint64_t gas_cell_id = 0;
  std::uint32_t local_cell_row = 0;
  std::uint32_t owning_rank = 0;
  std::uint64_t birth_tick = 0;
  std::uint32_t spawn_count = 0;
  double total_spawn_mass_code = 0.0;
  double parent_mass_code = 0.0;
  double parent_density_code = 0.0;
  double parent_pressure_code = 0.0;
  double parent_specific_internal_energy_code = 0.0;
  double parent_metal_mass_code = 0.0;
  double parent_metallicity_mass_fraction = 0.0;
  double parent_velocity_x_peculiar = 0.0;
  double parent_velocity_y_peculiar = 0.0;
  double parent_velocity_z_peculiar = 0.0;
  double parent_center_x_comoving = 0.0;
  double parent_center_y_comoving = 0.0;
  double parent_center_z_comoving = 0.0;
};

}  // namespace

std::uint64_t starFormationBirthKey(
    std::uint64_t gas_cell_id,
    std::uint64_t global_integration_tick,
    std::uint32_t birth_attempt_ordinal,
    std::uint32_t model_schema_version) {
  std::uint64_t mixed = splitmix64(gas_cell_id ^ 0x53f09d9d7a0f43b1ull);
  mixed ^= splitmix64(global_integration_tick ^ 0xc2b2ae3d27d4eb4full);
  mixed ^= splitmix64(static_cast<std::uint64_t>(birth_attempt_ordinal) ^ 0x165667b19e3779f9ull);
  mixed ^= splitmix64(static_cast<std::uint64_t>(model_schema_version) ^ 0x85ebca77c2b2ae63ull);
  const std::uint64_t result = splitmix64(mixed);
  return result == 0U ? 1U : result;
}

double starFormationUniform01(
    std::uint64_t global_seed,
    std::uint64_t gas_cell_id,
    std::uint64_t global_integration_tick,
    std::uint32_t birth_attempt_ordinal,
    std::uint32_t rng_key_schema_version) {
  std::uint64_t mixed = splitmix64(global_seed ^ 0xd6e8feb86659fd93ull);
  mixed ^= splitmix64(gas_cell_id ^ 0xa0761d6478bd642full);
  mixed ^= splitmix64(global_integration_tick ^ 0xe7037ed1a0b428dbull);
  mixed ^= splitmix64(static_cast<std::uint64_t>(birth_attempt_ordinal) ^ 0x8ebc6af09c88c6e3ull);
  mixed ^= splitmix64(static_cast<std::uint64_t>(rng_key_schema_version) ^ 0x589965cc75374cc3ull);
  const std::uint64_t bits = splitmix64(mixed);
  return static_cast<double>(bits >> 11U) * kU01Scale;
}

std::uint64_t starFormationParticleIdFromBirthKey(
    std::uint64_t birth_key,
    std::uint32_t collision_ordinal) {
  std::uint64_t particle_id = kGeneratedParticleIdTag |
      (splitmix64(birth_key ^ kParticleIdDomain ^
          splitmix64(static_cast<std::uint64_t>(collision_ordinal))) & ~kGeneratedParticleIdTag);
  if (particle_id == 0U) particle_id = kGeneratedParticleIdTag | 1U;
  return particle_id;
}

std::vector<std::uint64_t> precommitStarParticleIdsExact(
    std::span<const std::uint64_t> existing_particle_ids,
    std::span<const std::uint64_t> birth_keys) {
  std::vector<std::uint64_t> occupied(existing_particle_ids.begin(), existing_particle_ids.end());
  std::sort(occupied.begin(), occupied.end());
  if ((!occupied.empty() && occupied.front() == 0U) ||
      std::adjacent_find(occupied.begin(), occupied.end()) != occupied.end()) {
    throw std::runtime_error("ParticleIdRegistry: existing particle IDs are zero or duplicated");
  }

  std::vector<std::uint64_t> sorted_birth_keys(birth_keys.begin(), birth_keys.end());
  std::sort(sorted_birth_keys.begin(), sorted_birth_keys.end());
  if ((!sorted_birth_keys.empty() && sorted_birth_keys.front() == 0U) ||
      std::adjacent_find(sorted_birth_keys.begin(), sorted_birth_keys.end()) !=
          sorted_birth_keys.end()) {
    throw std::runtime_error(
        "ParticleIdRegistry: zero or duplicate immutable birth key in precommit batch");
  }

  struct CandidateRecord {
    std::uint64_t particle_id = 0U;
    std::uint64_t birth_key = 0U;
    std::size_t original_index = 0U;
    std::uint32_t collision_ordinal = 0U;
  };
  std::vector<CandidateRecord> candidates;
  candidates.reserve(birth_keys.size());
  for (std::size_t index = 0; index < birth_keys.size(); ++index) {
    candidates.push_back(CandidateRecord{
        .particle_id = starFormationParticleIdFromBirthKey(birth_keys[index], 0U),
        .birth_key = birth_keys[index],
        .original_index = index,
        .collision_ordinal = 0U,
    });
  }

  for (std::uint32_t pass = 0U; pass < 1024U; ++pass) {
    std::sort(candidates.begin(), candidates.end(), [](const CandidateRecord& lhs,
                                                       const CandidateRecord& rhs) {
      if (lhs.particle_id != rhs.particle_id) return lhs.particle_id < rhs.particle_id;
      return lhs.birth_key < rhs.birth_key;
    });
    bool had_collision = false;
    std::size_t begin = 0U;
    while (begin < candidates.size()) {
      std::size_t end_group = begin + 1U;
      while (end_group < candidates.size() &&
             candidates[end_group].particle_id == candidates[begin].particle_id) {
        ++end_group;
      }
      const bool collides_existing = std::binary_search(
          occupied.begin(), occupied.end(), candidates[begin].particle_id);
      const std::size_t first_to_rehash = collides_existing ? begin : begin + 1U;
      for (std::size_t index = first_to_rehash; index < end_group; ++index) {
        CandidateRecord& candidate = candidates[index];
        if (candidate.collision_ordinal == std::numeric_limits<std::uint32_t>::max()) {
          throw std::runtime_error(
              "ParticleIdRegistry: deterministic collision ordinal overflow");
        }
        ++candidate.collision_ordinal;
        candidate.particle_id = starFormationParticleIdFromBirthKey(
            candidate.birth_key, candidate.collision_ordinal);
        had_collision = true;
      }
      begin = end_group;
    }
    if (!had_collision) {
      std::vector<std::uint64_t> result(birth_keys.size(), 0U);
      for (const CandidateRecord& candidate : candidates) {
        result[candidate.original_index] = candidate.particle_id;
      }
      return result;
    }
  }
  throw std::runtime_error(
      "ParticleIdRegistry: deterministic collision resolution exhausted");
}

std::vector<std::uint64_t> LocalParticleIdRegistry::precommit(
    const core::SimulationState& state,
    std::span<const std::uint64_t> birth_keys) {
  return precommitStarParticleIdsExact(state.particle_sidecar.particle_id, birth_keys);
}

std::size_t starFormationBirthPlanBytes() noexcept {
  return sizeof(StarBirthPlan);
}

StarFormationModel::StarFormationModel(
    StarFormationConfig config,
    std::shared_ptr<const EffectiveMultiphaseEosTable> effective_eos_table)
    : m_config(std::move(config)), m_effective_eos_table(std::move(effective_eos_table)) {
  if (!finitePositive(m_config.newton_g_code)) {
    throw std::invalid_argument("StarFormationModel: newton_g_code must be finite and > 0");
  }
  if (!finiteNonNegative(m_config.epsilon_ff) || m_config.epsilon_ff > 1.0) {
    throw std::invalid_argument("StarFormationModel: epsilon_ff must be in [0, 1]");
  }
  if (!finitePositive(m_config.min_star_particle_mass_code) ||
      !finitePositive(m_config.max_star_particle_mass_code) ||
      m_config.max_star_particle_mass_code < m_config.min_star_particle_mass_code) {
    throw std::invalid_argument("StarFormationModel: invalid star-particle mass bounds");
  }
  if (m_config.max_spawn_particles_per_cell_step == 0U) {
    throw std::invalid_argument("StarFormationModel: max spawn count must be > 0");
  }
  if (!(m_config.max_fractional_mass_conversion > 0.0 &&
        m_config.max_fractional_mass_conversion < 1.0)) {
    throw std::invalid_argument("StarFormationModel: max fractional conversion must be in (0, 1)");
  }
  if (!(m_config.min_remaining_gas_fraction >= 0.0 &&
        m_config.min_remaining_gas_fraction < 1.0) ||
      m_config.min_remaining_gas_mass_code < 0.0) {
    throw std::invalid_argument("StarFormationModel: invalid remaining-gas constraints");
  }
  if (m_config.model == core::StarFormationModelKind::kLegacySchmidtThreshold) {
    if (!finitePositive(m_config.density_threshold_code) ||
        !finitePositive(m_config.temperature_threshold_k)) {
      throw std::invalid_argument("StarFormationModel: legacy thresholds must be finite and > 0");
    }
  } else if (m_config.model == core::StarFormationModelKind::kAdaptiveBoundJeans) {
    if (!finitePositive(m_config.bound_alpha_vir_max) ||
        !finiteNonNegative(m_config.jeans_mass_floor_code) ||
        !finitePositive(m_config.target_star_particle_mass_code) ||
        !(m_config.target_star_particle_mass_fraction > 0.0 &&
          m_config.target_star_particle_mass_fraction <= 1.0)) {
      throw std::invalid_argument("StarFormationModel: invalid adaptive_bound_jeans parameters");
    }
  } else {
    if (m_effective_eos_table == nullptr ||
        !finiteNonNegative(m_config.effective_min_baryon_overdensity) ||
        !finiteNonNegative(m_config.effective_hot_excess_tolerance) ||
        m_config.effective_massive_star_fraction < 0.0 ||
        m_config.effective_massive_star_fraction >= 1.0) {
      throw std::invalid_argument("StarFormationModel: effective multiphase model requires a valid EOS table");
    }
  }
}

const StarFormationConfig& StarFormationModel::config() const noexcept { return m_config; }

double StarFormationModel::physicalDensityCode(double stored_density_code, double scale_factor) const {
  if (!finitePositive(stored_density_code) || !finitePositive(scale_factor)) {
    return 0.0;
  }
  if (!m_config.density_is_comoving) {
    return stored_density_code;
  }
  return stored_density_code / (scale_factor * scale_factor * scale_factor);
}

double StarFormationModel::physicalCellScaleCode(double stored_volume_code, double scale_factor) const {
  if (!finitePositive(stored_volume_code) || !finitePositive(scale_factor)) {
    return 0.0;
  }
  const double length_stored = std::cbrt(stored_volume_code);
  return m_config.geometry_is_comoving ? length_stored * scale_factor : length_stored;
}

double StarFormationModel::virialParameter(
    double physical_density_code,
    double physical_cell_scale_code,
    double sound_speed_code,
    double velocity_gradient_frobenius_sq_code) const {
  if (!finitePositive(physical_density_code) || !finitePositive(physical_cell_scale_code) ||
      !finiteNonNegative(sound_speed_code) ||
      !finiteNonNegative(velocity_gradient_frobenius_sq_code)) {
    return std::numeric_limits<double>::infinity();
  }
  const double inverse_scale = 1.0 / physical_cell_scale_code;
  const double support_rate_sq = sound_speed_code * sound_speed_code * inverse_scale * inverse_scale;
  const double denominator = 8.0 * core::constants::k_pi * m_config.newton_g_code * physical_density_code;
  return (velocity_gradient_frobenius_sq_code + 2.0 * support_rate_sq) /
      std::max(denominator, kDensityFloor);
}

double StarFormationModel::jeansMassCode(double physical_density_code, double sound_speed_code) const {
  if (!finitePositive(physical_density_code) || !finiteNonNegative(sound_speed_code)) {
    return std::numeric_limits<double>::infinity();
  }
  const double coefficient = std::pow(core::constants::k_pi, 2.5) / 6.0;
  const double denominator = std::pow(m_config.newton_g_code, 1.5) * std::sqrt(physical_density_code);
  return coefficient * sound_speed_code * sound_speed_code * sound_speed_code /
      std::max(denominator, kDensityFloor);
}

double StarFormationModel::freeFallTimeCode(double gas_density_code) const {
  if (!finitePositive(gas_density_code)) {
    return std::numeric_limits<double>::infinity();
  }
  return std::sqrt(
      3.0 * core::constants::k_pi /
      (32.0 * m_config.newton_g_code * std::max(gas_density_code, kDensityFloor)));
}

double StarFormationModel::compressionTimeCode(double velocity_divergence_code) const {
  if (!std::isfinite(velocity_divergence_code) || velocity_divergence_code >= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  return -1.0 / velocity_divergence_code;
}

double StarFormationModel::sourceTimeStepLimitCode(double collapse_time_code) const {
  if (!finitePositive(collapse_time_code) || m_config.epsilon_ff <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  return -collapse_time_code * std::log1p(-m_config.max_fractional_mass_conversion) /
      m_config.epsilon_ff;
}

double StarFormationModel::sfrDensityRateCode(double gas_density_code) const {
  if (!finitePositive(gas_density_code)) {
    return 0.0;
  }
  const double t_ff = freeFallTimeCode(gas_density_code);
  return m_config.epsilon_ff * gas_density_code / std::max(t_ff, kTimeFloor);
}

StarFormationCellOutcome StarFormationModel::evaluateCell(
    const StarFormationCellInput& cell,
    double dt_code,
    double scale_factor,
    std::uint64_t global_integration_tick) const {
  StarFormationCellOutcome outcome;
  if (!m_config.enabled) {
    outcome.rejection_reason = StarFormationRejectionReason::kInactive;
    return outcome;
  }
  if (!cell.is_active) {
    outcome.rejection_reason = StarFormationRejectionReason::kInactive;
    return outcome;
  }
  if (!cell.is_owned) {
    outcome.rejection_reason = StarFormationRejectionReason::kNotOwned;
    return outcome;
  }
  if (!cell.is_leaf) {
    outcome.rejection_reason = StarFormationRejectionReason::kNonLeaf;
    return outcome;
  }
  if (cell.is_ghost) {
    outcome.rejection_reason = StarFormationRejectionReason::kGhost;
    return outcome;
  }
  if (!finitePositive(cell.gas_mass_code)) {
    outcome.rejection_reason = StarFormationRejectionReason::kNonPositiveMass;
    return outcome;
  }
  if (!finitePositive(cell.gas_density_code) || !finiteNonNegative(cell.gas_temperature_k) ||
      !finiteNonNegative(cell.gas_sound_speed_code) ||
      !std::isfinite(cell.velocity_divergence_code) ||
      !finiteNonNegative(cell.velocity_gradient_frobenius_sq_code) ||
      !std::isfinite(cell.velocity_x_peculiar) || !std::isfinite(cell.velocity_y_peculiar) ||
      !std::isfinite(cell.velocity_z_peculiar) || !finiteNonNegative(cell.gas_metal_mass_code)) {
    outcome.rejection_reason = StarFormationRejectionReason::kInvalidState;
    return outcome;
  }

  if (m_config.model == core::StarFormationModelKind::kLegacySchmidtThreshold) {
    if (cell.gas_density_code < m_config.density_threshold_code) {
      outcome.rejection_reason = StarFormationRejectionReason::kLegacyDensity;
      return outcome;
    }
    if (cell.gas_temperature_k > m_config.temperature_threshold_k) {
      outcome.rejection_reason = StarFormationRejectionReason::kLegacyTemperature;
      return outcome;
    }
    if (cell.velocity_divergence_code > m_config.min_converging_flow_rate_code) {
      outcome.rejection_reason = StarFormationRejectionReason::kNotConverging;
      return outcome;
    }
    outcome.physical_density_code = cell.gas_density_code;
    outcome.physical_cell_scale_code =
        finitePositive(cell.cell_volume_code) ? std::cbrt(cell.cell_volume_code) : 0.0;
  } else if (m_config.model == core::StarFormationModelKind::kAdaptiveBoundJeans) {
    if (cell.gas_cell_id == 0U) {
      outcome.rejection_reason = StarFormationRejectionReason::kInvalidState;
      return outcome;
    }
    if (!finitePositive(cell.cell_volume_code)) {
      outcome.rejection_reason = StarFormationRejectionReason::kNonPositiveVolume;
      return outcome;
    }
    outcome.physical_density_code = physicalDensityCode(cell.gas_density_code, scale_factor);
    outcome.physical_cell_scale_code = physicalCellScaleCode(cell.cell_volume_code, scale_factor);
    if (!finitePositive(outcome.physical_density_code) ||
        !finitePositive(outcome.physical_cell_scale_code)) {
      outcome.rejection_reason = StarFormationRejectionReason::kInvalidState;
      return outcome;
    }
    if (m_config.temperature_safety_ceiling_k > 0.0 &&
        cell.gas_temperature_k > m_config.temperature_safety_ceiling_k) {
      outcome.rejection_reason = StarFormationRejectionReason::kLegacyTemperature;
      return outcome;
    }
    if (m_config.require_converging_flow && cell.velocity_divergence_code >= 0.0) {
      outcome.rejection_reason = StarFormationRejectionReason::kNotConverging;
      return outcome;
    }
    outcome.virial_parameter = virialParameter(
        outcome.physical_density_code,
        outcome.physical_cell_scale_code,
        cell.gas_sound_speed_code,
        cell.velocity_gradient_frobenius_sq_code);
    if (!(outcome.virial_parameter < m_config.bound_alpha_vir_max)) {
      outcome.rejection_reason = StarFormationRejectionReason::kUnbound;
      return outcome;
    }
    outcome.jeans_mass_code = jeansMassCode(outcome.physical_density_code, cell.gas_sound_speed_code);
    const double instability_mass_code = std::max(m_config.jeans_mass_floor_code, cell.gas_mass_code);
    if (!(outcome.jeans_mass_code < instability_mass_code)) {
      outcome.rejection_reason = StarFormationRejectionReason::kJeansStable;
      return outcome;
    }
  } else {
    if (cell.gas_cell_id == 0U || m_effective_eos_table == nullptr ||
        !finitePositive(cell.gas_specific_internal_energy_code)) {
      outcome.rejection_reason = StarFormationRejectionReason::kEffectiveInvalidEquilibrium;
      return outcome;
    }
    outcome.physical_density_code = physicalDensityCode(cell.gas_density_code, scale_factor);
    outcome.physical_cell_scale_code = finitePositive(cell.cell_volume_code)
        ? physicalCellScaleCode(cell.cell_volume_code, scale_factor)
        : 0.0;
    const EffectiveMultiphaseEosLookup equilibrium =
        m_effective_eos_table->lookup(outcome.physical_density_code);
    if (!equilibrium.above_threshold) {
      outcome.rejection_reason = StarFormationRejectionReason::kEffectiveBelowDensityThreshold;
      return outcome;
    }
    if (cell.is_cosmological && m_config.effective_min_baryon_overdensity > 0.0 &&
        cell.baryon_overdensity < m_config.effective_min_baryon_overdensity) {
      outcome.rejection_reason = StarFormationRejectionReason::kEffectiveBelowOverdensityThreshold;
      return outcome;
    }
    if (!equilibrium.valid || !finitePositive(equilibrium.entry.star_formation_timescale_code) ||
        !finiteNonNegative(equilibrium.entry.cold_mass_fraction)) {
      outcome.rejection_reason = StarFormationRejectionReason::kEffectiveInvalidEquilibrium;
      return outcome;
    }
    outcome.cold_mass_fraction = equilibrium.entry.cold_mass_fraction;
    outcome.effective_specific_internal_energy_code =
        equilibrium.entry.specific_internal_energy_eff_code;
    outcome.effective_pressure_code = (m_effective_eos_table->config().adiabatic_index - 1.0) *
        cell.gas_density_code * outcome.effective_specific_internal_energy_code;
    outcome.effective_signal_speed_squared_code =
        equilibrium.entry.signal_speed_squared_code;
    if (cell.gas_specific_internal_energy_code >
        outcome.effective_specific_internal_energy_code *
            (1.0 + m_config.effective_hot_excess_tolerance)) {
      outcome.rejection_reason = StarFormationRejectionReason::kEffectiveHotAboveEos;
      return outcome;
    }
    outcome.collapse_time_code = equilibrium.entry.star_formation_timescale_code;
    outcome.free_fall_time_code = freeFallTimeCode(outcome.physical_density_code);
    outcome.compression_time_code = compressionTimeCode(cell.velocity_divergence_code);
  }

  outcome.eligible = true;
  outcome.rejection_reason = StarFormationRejectionReason::kNone;
  if (m_config.model != core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
    outcome.free_fall_time_code = freeFallTimeCode(outcome.physical_density_code);
    outcome.compression_time_code = compressionTimeCode(cell.velocity_divergence_code);
    outcome.collapse_time_code = outcome.free_fall_time_code;
    if (m_config.collapse_timescale ==
        core::StarFormationCollapseTimescale::kMinimumFreeFallOrCompression) {
      outcome.collapse_time_code = std::min(outcome.free_fall_time_code, outcome.compression_time_code);
    }
    outcome.sfr_density_rate_code = m_config.epsilon_ff * outcome.physical_density_code /
        std::max(outcome.collapse_time_code, kTimeFloor);
  } else {
    const double long_lived_factor =
        m_config.effective_birth_mass_convention == core::EffectiveIsmBirthMassConvention::kLongLivedMass
        ? (1.0 - m_config.effective_massive_star_fraction)
        : 1.0;
    outcome.sfr_density_rate_code = long_lived_factor * outcome.cold_mass_fraction *
        outcome.physical_density_code / std::max(outcome.collapse_time_code, kTimeFloor);
  }

  if (!(dt_code > 0.0) || !finitePositive(outcome.collapse_time_code)) {
    return outcome;
  }

  double expected_mass_code = 0.0;
  if (m_config.model == core::StarFormationModelKind::kLegacySchmidtThreshold) {
    if (m_config.epsilon_ff <= 0.0) return outcome;
    expected_mass_code = m_config.epsilon_ff * cell.gas_mass_code * dt_code /
        std::max(outcome.collapse_time_code, kTimeFloor);
  } else if (m_config.model == core::StarFormationModelKind::kAdaptiveBoundJeans) {
    if (m_config.epsilon_ff <= 0.0) return outcome;
    const double exponent = -m_config.epsilon_ff * dt_code /
        std::max(outcome.collapse_time_code, kTimeFloor);
    expected_mass_code = cell.gas_mass_code * (-std::expm1(exponent));
  } else {
    const double long_lived_factor =
        m_config.effective_birth_mass_convention == core::EffectiveIsmBirthMassConvention::kLongLivedMass
        ? (1.0 - m_config.effective_massive_star_fraction)
        : 1.0;
    const double rate_per_mass = long_lived_factor * outcome.cold_mass_fraction /
        std::max(outcome.collapse_time_code, kTimeFloor);
    expected_mass_code = cell.gas_mass_code * (-std::expm1(-rate_per_mass * dt_code));
  }

  const double max_by_step = cell.gas_mass_code * m_config.max_fractional_mass_conversion;
  const double max_by_fraction = cell.gas_mass_code * (1.0 - m_config.min_remaining_gas_fraction);
  const double max_by_absolute = std::max(0.0, cell.gas_mass_code - m_config.min_remaining_gas_mass_code);
  const double maximum_transfer = std::max(
      0.0,
      std::min({cell.gas_mass_code, max_by_step, max_by_fraction, max_by_absolute}));
  outcome.expected_spawn_mass_code = std::clamp(expected_mass_code, 0.0, maximum_transfer);
  if (!(outcome.expected_spawn_mass_code > 0.0) ||
      maximum_transfer < m_config.min_star_particle_mass_code) {
    return outcome;
  }
  // Deterministic births cannot produce a star lighter than the configured
  // minimum.  Return before the stochastic target/count feasibility search:
  // a subminimum expected mass is a legitimate no-birth step, not an error.
  if (!m_config.stochastic_spawning &&
      outcome.expected_spawn_mass_code < m_config.min_star_particle_mass_code) {
    return outcome;
  }

  if (m_config.model == core::StarFormationModelKind::kLegacySchmidtThreshold ||
      m_config.star_particle_mass_policy == core::StarParticleMassPolicy::kFixed) {
    outcome.target_particle_mass_code = m_config.target_star_particle_mass_code;
    if (m_config.model == core::StarFormationModelKind::kLegacySchmidtThreshold) {
      outcome.target_particle_mass_code = m_config.min_star_particle_mass_code;
    }
  } else {
    outcome.target_particle_mass_code = cell.gas_mass_code * m_config.target_star_particle_mass_fraction;
  }
  outcome.target_particle_mass_code = std::clamp(
      outcome.target_particle_mass_code,
      m_config.min_star_particle_mass_code,
      std::min(m_config.max_star_particle_mass_code, maximum_transfer));

  // Find the closest feasible effective target mass for which both possible
  // integer-plus-Bernoulli outcomes fit inside the available gas mass. This
  // avoids the biased post-draw clipping that occurs when ceil(lambda) times a
  // nominal target mass exceeds maximum_transfer.
  const double nominal_target_mass_code = outcome.target_particle_mass_code;
  double feasible_target_mass_code = 0.0;
  double best_target_distance = std::numeric_limits<double>::infinity();
  for (std::uint32_t maximum_count = 1U;
       maximum_count <= m_config.max_spawn_particles_per_cell_step;
       ++maximum_count) {
    const double count = static_cast<double>(maximum_count);
    const double lower = std::max(
        m_config.min_star_particle_mass_code,
        outcome.expected_spawn_mass_code / count);
    double upper = std::min(
        m_config.max_star_particle_mass_code,
        maximum_transfer / count);
    if (maximum_count > 1U) {
      upper = std::min(
          upper,
          std::nextafter(
              outcome.expected_spawn_mass_code /
                  static_cast<double>(maximum_count - 1U),
              0.0));
    }
    if (!(lower <= upper) || !finitePositive(lower) || !finitePositive(upper)) {
      continue;
    }
    const double candidate = std::clamp(nominal_target_mass_code, lower, upper);
    const double distance = std::abs(candidate - nominal_target_mass_code);
    if (distance < best_target_distance) {
      best_target_distance = distance;
      feasible_target_mass_code = candidate;
    }
  }
  if (!finitePositive(feasible_target_mass_code)) {
    std::ostringstream error;
    error << "StarFormationModel: configured particle-mass/count bounds cannot represent bounded expected mass"
          << " (expected_mass_code=" << outcome.expected_spawn_mass_code
          << ", maximum_transfer_code=" << maximum_transfer
          << ", minimum_particle_mass_code=" << m_config.min_star_particle_mass_code
          << ", maximum_particle_mass_code=" << m_config.max_star_particle_mass_code
          << ", max_spawn_particles_per_cell_step=" << m_config.max_spawn_particles_per_cell_step
          << ", stochastic=" << m_config.stochastic_spawning << ')';
    throw std::runtime_error(error.str());
  }
  outcome.target_particle_mass_code = feasible_target_mass_code;

  if (!m_config.stochastic_spawning) {
    if (outcome.expected_spawn_mass_code < m_config.min_star_particle_mass_code) {
      return outcome;
    }
    outcome.spawned_particle_count = 1U;
    outcome.spawned_mass_code = outcome.expected_spawn_mass_code;
    outcome.spawn_probability = 1.0;
    return outcome;
  }

  const double lambda = outcome.expected_spawn_mass_code /
      std::max(outcome.target_particle_mass_code, kMassFloor);
  std::uint64_t base_count = static_cast<std::uint64_t>(std::floor(lambda));
  const double remainder_probability = std::clamp(lambda - static_cast<double>(base_count), 0.0, 1.0);
  outcome.spawn_probability = remainder_probability;
  outcome.random_u01 = starFormationUniform01(
      m_config.random_seed,
      cell.gas_cell_id != 0U ? cell.gas_cell_id : static_cast<std::uint64_t>(cell.cell_index) + 1U,
      global_integration_tick,
      0U);
  if (outcome.random_u01 < remainder_probability) {
    ++base_count;
  }
  if (base_count > m_config.max_spawn_particles_per_cell_step) {
    throw std::runtime_error(
        "StarFormationModel: target-mass/count policy cannot represent bounded expected mass");
  }
  if (base_count == 0U) {
    return outcome;
  }
  outcome.spawned_particle_count = static_cast<std::uint32_t>(base_count);
  outcome.spawned_mass_code =
      static_cast<double>(outcome.spawned_particle_count) * outcome.target_particle_mass_code;
  if (outcome.spawned_mass_code > maximum_transfer * (1.0 + 8.0 * std::numeric_limits<double>::epsilon())) {
    throw std::runtime_error(
        "StarFormationModel: feasible stochastic outcome exceeded available gas mass");
  }
  if (outcome.spawned_mass_code < m_config.min_star_particle_mass_code) {
    outcome.spawned_particle_count = 0U;
    outcome.spawned_mass_code = 0.0;
  }
  return outcome;
}

bool StarFormationModel::isEligible(const StarFormationCellInput& cell) const {
  return evaluateCell(cell, 0.0, 1.0, 0U).eligible;
}

double StarFormationModel::expectedSpawnMassCode(
    const StarFormationCellInput& cell,
    double dt_code) const {
  return evaluateCell(cell, dt_code, 1.0, 0U).expected_spawn_mass_code;
}

StarFormationCellOutcome StarFormationModel::sampleCellOutcome(
    const StarFormationCellInput& cell,
    double dt_code,
    std::uint64_t step_index,
    std::uint32_t ignored_rank_local_seed_offset,
    double scale_factor) const {
  (void)ignored_rank_local_seed_offset;
  return evaluateCell(cell, dt_code, scale_factor, step_index);
}

StarFormationStepReport StarFormationModel::applyFromInputs(
    core::SimulationState& state,
    std::span<const StarFormationCellInput> cell_inputs,
    double dt_code,
    double scale_factor,
    std::uint64_t global_integration_tick,
    ParticleIdPrecommit* id_precommit,
    std::pmr::memory_resource* metadata_scratch) const {
  StarFormationStepReport report;
  if (!m_config.enabled || !(dt_code > 0.0)) {
    return report;
  }

  const std::uint64_t identity_generation_before = state.gasCellIdentityGeneration();
  // Birth plans are event-local. Do not reserve one 136-byte record for every
  // scanned gas cell before eligibility is known.
  std::pmr::vector<StarBirthPlan> plans(metadata_scratch == nullptr
      ? std::pmr::get_default_resource() : metadata_scratch);
  if (metadata_scratch != nullptr) plans.reserve(cell_inputs.size());

  for (const StarFormationCellInput& cell : cell_inputs) {
    ++report.counters.scanned_cells;
    if (m_config.model == core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
      ++report.counters.effective_cells_scanned;
    }
    const StarFormationCellOutcome outcome = evaluateCell(
        cell, dt_code, scale_factor, global_integration_tick);
    if (!outcome.eligible) {
      incrementRejectionCounter(report.counters, outcome.rejection_reason);
      if (outcome.rejection_reason == StarFormationRejectionReason::kEffectiveHotAboveEos) {
        ++report.counters.effective_cells_hot_above_eos;
      } else if (outcome.rejection_reason == StarFormationRejectionReason::kEffectiveInvalidEquilibrium) {
        ++report.counters.effective_cells_invalid_equilibrium;
      }
      continue;
    }
    ++report.counters.eligible_cells;
    if (m_config.model == core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
      ++report.counters.effective_cells_above_threshold;
      ++report.counters.effective_cells_on_eos;
      report.counters.effective_cold_gas_mass_code += outcome.cold_mass_fraction * cell.gas_mass_code;
      report.counters.effective_instantaneous_sfr_code +=
          outcome.sfr_density_rate_code * cell.cell_volume_code;
      report.counters.effective_min_cold_fraction = std::min(
          report.counters.effective_min_cold_fraction, outcome.cold_mass_fraction);
      report.counters.effective_max_cold_fraction = std::max(
          report.counters.effective_max_cold_fraction, outcome.cold_mass_fraction);
      if (report.counters.effective_min_pressure_code == 0.0) {
        report.counters.effective_min_pressure_code = outcome.effective_pressure_code;
      } else {
        report.counters.effective_min_pressure_code = std::min(
            report.counters.effective_min_pressure_code, outcome.effective_pressure_code);
      }
      report.counters.effective_max_pressure_code = std::max(
          report.counters.effective_max_pressure_code, outcome.effective_pressure_code);
    }
    report.counters.expected_spawn_mass_code += outcome.expected_spawn_mass_code;
    if (outcome.free_fall_time_code > 0.0 && std::isfinite(outcome.free_fall_time_code)) {
      if (report.counters.minimum_free_fall_time_code == 0.0) {
        report.counters.minimum_free_fall_time_code = outcome.free_fall_time_code;
      } else {
        report.counters.minimum_free_fall_time_code =
            std::min(report.counters.minimum_free_fall_time_code, outcome.free_fall_time_code);
      }
    }
    if (outcome.collapse_time_code > 0.0 && std::isfinite(outcome.collapse_time_code)) {
      if (report.counters.minimum_collapse_time_code == 0.0) {
        report.counters.minimum_collapse_time_code = outcome.collapse_time_code;
      } else {
        report.counters.minimum_collapse_time_code =
            std::min(report.counters.minimum_collapse_time_code, outcome.collapse_time_code);
      }
    }
    report.counters.maximum_spawn_probability = std::max(
        report.counters.maximum_spawn_probability,
        outcome.spawn_probability);
    if (outcome.spawned_particle_count == 0U || outcome.spawned_mass_code <= kMassFloor) {
      continue;
    }
    if (cell.cell_index >= state.cells.size() || cell.cell_index >= state.gas_cells.size()) {
      ++report.counters.rejected_invalid_state;
      continue;
    }

    const double authoritative_mass = state.cells.mass_code[cell.cell_index];
    if (!finitePositive(authoritative_mass) ||
        std::abs(authoritative_mass - cell.gas_mass_code) >
            1.0e-12 * std::max(1.0, std::abs(authoritative_mass))) {
      throw std::runtime_error("StarFormationModel: stale birth plan gas mass");
    }

    StarBirthPlan plan;
    plan.gas_cell_id = cell.gas_cell_id != 0U
        ? cell.gas_cell_id
        : static_cast<std::uint64_t>(cell.cell_index) + 1U;
    plan.local_cell_row = cell.cell_index;
    plan.owning_rank = cell.owning_rank;
    plan.birth_tick = global_integration_tick;
    plan.spawn_count = outcome.spawned_particle_count;
    plan.total_spawn_mass_code = std::min(outcome.spawned_mass_code, authoritative_mass);
    plan.parent_mass_code = authoritative_mass;
    plan.parent_density_code = state.gas_cells.density_code[cell.cell_index];
    plan.parent_pressure_code = state.gas_cells.pressure_code[cell.cell_index];
    plan.parent_specific_internal_energy_code = state.gas_cells.internal_energy_code[cell.cell_index];
    const double state_metal_mass = state.gas_cells.metal_mass_code[cell.cell_index];
    plan.parent_metal_mass_code = state_metal_mass > 0.0
        ? state_metal_mass
        : safeMetallicity(cell) * authoritative_mass;
    plan.parent_metallicity_mass_fraction = std::clamp(
        plan.parent_metal_mass_code / std::max(authoritative_mass, kMassFloor), 0.0, 1.0);
    plan.parent_velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[cell.cell_index];
    plan.parent_velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[cell.cell_index];
    plan.parent_velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[cell.cell_index];
    plan.parent_center_x_comoving = state.cells.center_x_comoving[cell.cell_index];
    plan.parent_center_y_comoving = state.cells.center_y_comoving[cell.cell_index];
    plan.parent_center_z_comoving = state.cells.center_z_comoving[cell.cell_index];
    plans.push_back(plan);
  }

  std::sort(plans.begin(), plans.end(), [](const StarBirthPlan& lhs, const StarBirthPlan& rhs) {
    if (lhs.gas_cell_id != rhs.gas_cell_id) {
      return lhs.gas_cell_id < rhs.gas_cell_id;
    }
    return lhs.local_cell_row < rhs.local_cell_row;
  });
  for (std::size_t i = 1; i < plans.size(); ++i) {
    if (plans[i - 1].gas_cell_id == plans[i].gas_cell_id) {
      throw std::runtime_error("StarFormationModel: duplicate gas-cell birth plan");
    }
  }
  if (state.gasCellIdentityGeneration() != identity_generation_before) {
    throw std::runtime_error("StarFormationModel: gas-cell identity changed between evaluation and apply");
  }

  std::uint64_t total_new_particles_u64 = 0U;
  for (const StarBirthPlan& plan : plans) {
    total_new_particles_u64 += plan.spawn_count;
  }
  if (total_new_particles_u64 > std::numeric_limits<std::uint32_t>::max()) {
    throw std::overflow_error("StarFormationModel: source step creates too many particles");
  }
  const std::size_t total_new_particles = static_cast<std::size_t>(total_new_particles_u64);

  // IDs are deterministic functions of immutable birth identity. Validate the
  // complete local birth batch without a full scan of the existing particle set
  // or one heap allocation per ID. The workflow ownership gate performs the exact
  // global duplicate-ID check after legal source mutation and before acceptance.
  std::pmr::vector<std::uint64_t> new_birth_keys(metadata_scratch == nullptr
      ? std::pmr::get_default_resource() : metadata_scratch);
  new_birth_keys.reserve(total_new_particles);
  // The production owner has admitted both report arrays before entering the
  // distributed ID precommit. Exact reserves prevent allocator growth after
  // canonical state begins its birth transaction.
  if (metadata_scratch != nullptr) {
    report.spawned_from_cells.reserve(total_new_particles);
    report.birth_keys.reserve(total_new_particles);
  }
  for (const StarBirthPlan& plan : plans) {
    for (std::uint32_t ordinal = 0U; ordinal < plan.spawn_count; ++ordinal) {
      new_birth_keys.push_back(starFormationBirthKey(
          plan.gas_cell_id, plan.birth_tick, ordinal, m_config.metadata_schema_version));
    }
  }
  LocalParticleIdRegistry local_registry;
  ParticleIdPrecommit& registry = id_precommit != nullptr
      ? *id_precommit
      : static_cast<ParticleIdPrecommit&>(local_registry);
  const std::vector<std::uint64_t> new_particle_ids = registry.precommit(state, new_birth_keys);
  struct PrecommitResultScope {
    ParticleIdPrecommit& registry;
    ~PrecommitResultScope() { registry.finishPrecommit(); }
  } precommit_result_scope{registry};
  if (new_particle_ids.size() != new_birth_keys.size()) {
    throw std::runtime_error("StarFormationModel: particle-ID precommit returned the wrong batch size");
  }

  registry.preparePopulationGrowth(state, total_new_particles, core::ParticleSpecies::kStar);

  // Distributed precommit is collective. Empty-birth ranks must reach it for
  // every globally ordered source batch before returning.
  if (plans.empty()) {
    StarFormationCounters cumulative = cumulativeCountersFromState(state);
    accumulateStarFormationCounters(cumulative, report.counters);
    state.sidecars.upsert(buildMetadataSidecar(
        report.counters, state.metadata.normalized_config_hash_hex, &cumulative));
    return report;
  }

  const std::size_t old_particle_count = state.particles.size();
  const std::size_t old_star_count = state.star_particles.size();
  state.resizeParticles(old_particle_count + total_new_particles);
  state.star_particles.resize(old_star_count + total_new_particles);

  // CHUI Part 2: element-resolved birth composition. Stars are a material
  // reservoir: the removed fraction of every gas element, dust solid and SGS
  // energy lane leaves the cell (exact proportional removal) and the star
  // records its birth element mixture in chui_stellar_population_v1.
  const bool composition = state.gas_composition.enabled();
  std::optional<CompositionLayout> composition_layout;
  std::optional<StellarPopulationStore> population;
  if (composition) {
    composition_layout = CompositionLayout::fromDescriptor(state.gas_composition.descriptor());
    population.emplace(state);
  }
  std::size_t append_offset = 0U;
  for (const StarBirthPlan& plan : plans) {
    const double removal_fraction = plan.total_spawn_mass_code /
        std::max(plan.parent_mass_code, kMassFloor);
    const double remaining_fraction = std::clamp(1.0 - removal_fraction, 0.0, 1.0);
    const double removed_metal_mass = plan.parent_metallicity_mass_fraction * plan.total_spawn_mass_code;

    state.cells.mass_code[plan.local_cell_row] = plan.parent_mass_code - plan.total_spawn_mass_code;
    state.gas_cells.density_code[plan.local_cell_row] = plan.parent_density_code * remaining_fraction;
    state.gas_cells.pressure_code[plan.local_cell_row] = plan.parent_pressure_code * remaining_fraction;
    state.gas_cells.metal_mass_code[plan.local_cell_row] =
        std::max(0.0, plan.parent_metal_mass_code - removed_metal_mass);
    std::array<double, k_chemical_element_count> birth_x{};
    if (composition) {
      const auto removed = removeGasFractionFromComposition(
          state.gas_composition, *composition_layout, plan.local_cell_row, std::clamp(removal_fraction, 0.0, 1.0));
      double removed_total = 0.0;
      for (const double m : removed) removed_total += m;
      for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
        birth_x[e] = removed_total > 0.0 ? removed[e] / removed_total : 0.0;
      }
      const std::uint32_t row = plan.local_cell_row;
      (void)synchronizeMetalMirrorFromComposition(state, std::span<const std::uint32_t>(&row, 1U));
    }
    // velocity, temperature, sound speed, and specific internal energy remain unchanged
    // under proportional removal. Conserved momentum/energy are therefore reduced by
    // exactly the mass-proportional ledgers below when hydro rebuilds its conserved view.

    const double per_particle_mass = plan.total_spawn_mass_code /
        static_cast<double>(plan.spawn_count);
    for (std::uint32_t ordinal = 0U; ordinal < plan.spawn_count; ++ordinal) {
      const std::size_t particle_index = old_particle_count + append_offset;
      const std::size_t star_index = old_star_count + append_offset;
      state.particles.position_x_comoving[particle_index] = plan.parent_center_x_comoving;
      state.particles.position_y_comoving[particle_index] = plan.parent_center_y_comoving;
      state.particles.position_z_comoving[particle_index] = plan.parent_center_z_comoving;
      state.particles.velocity_x_peculiar[particle_index] = plan.parent_velocity_x_peculiar;
      state.particles.velocity_y_peculiar[particle_index] = plan.parent_velocity_y_peculiar;
      state.particles.velocity_z_peculiar[particle_index] = plan.parent_velocity_z_peculiar;
      state.particles.mass_code[particle_index] = per_particle_mass;
      state.notifyParticleMassMutation();
      state.particles.time_bin[particle_index] = 0U;
      state.particle_sidecar.particle_id[particle_index] = new_particle_ids[append_offset];
      state.particle_sidecar.sfc_key[particle_index] = new_particle_ids[append_offset];
      state.particle_sidecar.species_tag[particle_index] =
          static_cast<std::uint32_t>(core::ParticleSpecies::kStar);
      state.particle_sidecar.particle_flags[particle_index] = 0U;
      state.particle_sidecar.owning_rank[particle_index] = plan.owning_rank;
      state.particle_sidecar.last_drift_time_code[particle_index] = 0.0;
      state.particle_sidecar.last_drift_scale_factor[particle_index] = scale_factor;

      state.star_particles.particle_index[star_index] = static_cast<std::uint32_t>(particle_index);
      state.star_particles.formation_scale_factor[star_index] = scale_factor;
      state.star_particles.birth_mass_code[star_index] = per_particle_mass;
      state.star_particles.metallicity_mass_fraction[star_index] =
          plan.parent_metallicity_mass_fraction;
      state.star_particles.birth_key[star_index] = new_birth_keys[append_offset];
      state.star_particles.parent_gas_cell_id[star_index] = plan.gas_cell_id;
      state.star_particles.birth_tick[star_index] = plan.birth_tick;
      state.star_particles.birth_ordinal[star_index] = ordinal;
      state.star_particles.stellar_age_years_last[star_index] = 0.0;
      state.star_particles.stellar_returned_mass_cumulative_code[star_index] = 0.0;
      state.star_particles.stellar_returned_metals_cumulative_code[star_index] = 0.0;
      state.star_particles.stellar_newly_synthesized_metals_cumulative_code[star_index] = 0.0;
      state.star_particles.stellar_feedback_energy_cumulative_erg[star_index] = 0.0;
      state.star_particles.enrichment_carry_mass_code[star_index] = 0.0;
      state.star_particles.enrichment_carry_metals_code[star_index] = 0.0;
      state.star_particles.enrichment_carry_feedback_energy_erg[star_index] = 0.0;
      state.star_particles.enrichment_carry_momentum_code[star_index] = 0.0;
      state.star_particles.stellar_deposited_mass_cumulative_code[star_index] = 0.0;
      state.star_particles.stellar_deposited_metals_cumulative_code[star_index] = 0.0;
      state.star_particles.stellar_deposited_feedback_energy_cumulative_erg[star_index] = 0.0;
      for (std::size_t channel = 0;
           channel < state.star_particles.stellar_returned_mass_channel_cumulative_code.size();
           ++channel) {
        state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][star_index] = 0.0;
        state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][star_index] = 0.0;
        state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][star_index] = 0.0;
      }
      if (population.has_value()) {
        // Birth epoch: the end of the accepted interval (t_birth = t_end;
        // F5 chronology). Event stream keyed by the immutable birth key.
        population->put(new_particle_ids[append_offset],
                        makeStellarPopulationRecord(birth_x, m_birth_time_code, new_birth_keys[append_offset]));
      }
      report.birth_keys.push_back(new_birth_keys[append_offset]);
      ++append_offset;
    }

    ++report.counters.spawn_events;
    report.counters.spawned_particles += plan.spawn_count;
    report.counters.spawned_mass_code += plan.total_spawn_mass_code;
    report.spawned_from_cells.push_back(plan.local_cell_row);

    const double momentum_x = plan.total_spawn_mass_code * plan.parent_velocity_x_peculiar;
    const double momentum_y = plan.total_spawn_mass_code * plan.parent_velocity_y_peculiar;
    const double momentum_z = plan.total_spawn_mass_code * plan.parent_velocity_z_peculiar;
    const double speed_sq =
        plan.parent_velocity_x_peculiar * plan.parent_velocity_x_peculiar +
        plan.parent_velocity_y_peculiar * plan.parent_velocity_y_peculiar +
        plan.parent_velocity_z_peculiar * plan.parent_velocity_z_peculiar;
    const double internal_energy_removed =
        plan.total_spawn_mass_code * plan.parent_specific_internal_energy_code;
    const double kinetic_energy_created = 0.5 * plan.total_spawn_mass_code * speed_sq;

    report.counters.gas_mass_removed_code += plan.total_spawn_mass_code;
    report.counters.star_mass_created_code += plan.total_spawn_mass_code;
    report.counters.gas_momentum_removed_x_code += momentum_x;
    report.counters.gas_momentum_removed_y_code += momentum_y;
    report.counters.gas_momentum_removed_z_code += momentum_z;
    report.counters.star_momentum_created_x_code += momentum_x;
    report.counters.star_momentum_created_y_code += momentum_y;
    report.counters.star_momentum_created_z_code += momentum_z;
    report.counters.gas_metal_mass_removed_code += removed_metal_mass;
    report.counters.star_metal_mass_created_code += removed_metal_mass;
    report.counters.gas_internal_energy_removed_code += internal_energy_removed;
    report.counters.star_kinetic_energy_created_code += kinetic_energy_created;
    report.counters.star_formation_internal_energy_sink_code += internal_energy_removed;
    report.counters.maximum_fractional_mass_conversion = std::max(
        report.counters.maximum_fractional_mass_conversion, removal_fraction);
  }

  report.counters.mass_residual_code =
      report.counters.gas_mass_removed_code - report.counters.star_mass_created_code;
  const double residual_px = report.counters.gas_momentum_removed_x_code -
      report.counters.star_momentum_created_x_code;
  const double residual_py = report.counters.gas_momentum_removed_y_code -
      report.counters.star_momentum_created_y_code;
  const double residual_pz = report.counters.gas_momentum_removed_z_code -
      report.counters.star_momentum_created_z_code;
  report.counters.momentum_residual_norm_code =
      std::sqrt(residual_px * residual_px + residual_py * residual_py + residual_pz * residual_pz);
  report.counters.metal_mass_residual_code =
      report.counters.gas_metal_mass_removed_code - report.counters.star_metal_mass_created_code;

  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)] +=
      report.counters.spawned_particles;
  state.rebuildSpeciesIndex();
  if (population.has_value()) population->store(state);
  StarFormationCounters cumulative = cumulativeCountersFromState(state);
  accumulateStarFormationCounters(cumulative, report.counters);
  state.sidecars.upsert(buildMetadataSidecar(
      report.counters, state.metadata.normalized_config_hash_hex, &cumulative));
  return report;
}

StarFormationStepReport StarFormationModel::applyFromView(
    core::SimulationState& state,
    StarFormationRuntimeView view,
    double dt_code,
    double scale_factor,
    std::uint64_t step_index,
    std::uint32_t owning_rank) const {
  std::vector<StarFormationCellInput> inputs;
  inputs.reserve(view.active_cell_indices.size());
  for (const std::uint32_t cell_index : view.active_cell_indices) {
    if (cell_index >= state.cells.size() || cell_index >= state.gas_cells.size()) {
      continue;
    }
    StarFormationCellInput input;
    input.cell_index = cell_index;
    input.gas_cell_id = cell_index < state.gas_cells.gas_cell_id.size()
        ? state.gas_cells.gas_cell_id[cell_index]
        : 0U;
    input.owning_rank = owning_rank;
    input.gas_mass_code = state.cells.mass_code[cell_index];
    input.gas_density_code = state.gas_cells.density_code[cell_index];
    input.cell_volume_code = finitePositive(input.gas_density_code)
        ? input.gas_mass_code / input.gas_density_code
        : 0.0;
    input.gas_temperature_k = state.gas_cells.temperature_code[cell_index];
    input.gas_sound_speed_code = state.gas_cells.sound_speed_code[cell_index];
    input.velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[cell_index];
    input.velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[cell_index];
    input.velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[cell_index];
    input.velocity_divergence_code = cell_index < view.velocity_divergence_code.size()
        ? view.velocity_divergence_code[cell_index]
        : 0.0;
    input.velocity_gradient_frobenius_sq_code =
        input.velocity_divergence_code * input.velocity_divergence_code / 3.0;
    input.metallicity_mass_fraction = cell_index < view.metallicity_mass_fraction.size()
        ? view.metallicity_mass_fraction[cell_index]
        : 0.0;
    input.gas_metal_mass_code = state.gas_cells.metal_mass_code[cell_index] > 0.0
        ? state.gas_cells.metal_mass_code[cell_index]
        : input.metallicity_mass_fraction * input.gas_mass_code;
    input.center_x_comoving = state.cells.center_x_comoving[cell_index];
    input.center_y_comoving = state.cells.center_y_comoving[cell_index];
    input.center_z_comoving = state.cells.center_z_comoving[cell_index];
    inputs.push_back(input);
  }
  return applyFromInputs(state, inputs, dt_code, scale_factor, step_index);
}

StarFormationStepReport StarFormationModel::apply(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_indices,
    std::span<const double> velocity_divergence_code,
    std::span<const double> metallicity_mass_fraction,
    double dt_code,
    double scale_factor,
    std::uint64_t step_index,
    std::uint32_t owning_rank) const {
  StarFormationRuntimeView view{
      .active_cell_indices = active_cell_indices,
      .center_x_comoving = state.cells.center_x_comoving,
      .center_y_comoving = state.cells.center_y_comoving,
      .center_z_comoving = state.cells.center_z_comoving,
      .gas_mass_code = state.cells.mass_code,
      .gas_density_code = state.gas_cells.density_code,
      .gas_temperature_k = state.gas_cells.temperature_code,
      .velocity_divergence_code = velocity_divergence_code,
      .metallicity_mass_fraction = metallicity_mass_fraction,
  };
  return applyFromView(state, view, dt_code, scale_factor, step_index, owning_rank);
}

core::ModuleSidecarBlock StarFormationModel::buildMetadataSidecar(
    const StarFormationCounters& counters,
    std::string_view configuration_hash,
    const StarFormationCounters* cumulative_counters) const {
  std::ostringstream stream;
  stream << "module_name=star_formation\n";
  stream << "model_name=" << core::starFormationModelKindToString(m_config.model) << "\n";
  stream << "model_schema_version=" << m_config.metadata_schema_version << "\n";
  stream << "rng_algorithm=splitmix64_counter_keyed\n";
  stream << "rng_key_schema_version=" << kStarFormationRngKeySchemaVersion << "\n";
  stream << "particle_id_collision_detection=exact_existing_local_and_distributed_precommit\n";
  stream << "configuration_hash=" << configuration_hash << "\n";
  stream << "counter_scope=latest_and_cumulative_source_stages\n";
  stream << "physical_density_convention="
         << (m_config.density_is_comoving ? "rho_stored_div_a_cubed" : "stored_density_is_physical") << "\n";
  stream << "velocity_gradient_convention=physical_inverse_time_patch_finite_difference\n";
  stream << "jeans_mass_convention=pi_pow_5_over_2_over_6\n";
  stream << "particle_mass_policy="
         << core::starParticleMassPolicyToString(m_config.star_particle_mass_policy) << "\n";
  stream << "epsilon_ff=" << m_config.epsilon_ff << "\n";
  stream << "bound_alpha_vir_max=" << m_config.bound_alpha_vir_max << "\n";
  stream << "jeans_mass_floor_code=" << m_config.jeans_mass_floor_code << "\n";
  stream << "random_seed=" << m_config.random_seed << "\n";
  // Flat latest-step compatibility fields remain available while the complete
  // namespaced latest/cumulative ledgers provide restart-persistent accounting.
  stream << "scanned_cells=" << counters.scanned_cells << "\n";
  stream << "eligible_cells=" << counters.eligible_cells << "\n";
  stream << "spawn_events=" << counters.spawn_events << "\n";
  stream << "spawned_particles=" << counters.spawned_particles << "\n";
  stream << "expected_spawn_mass_code=" << counters.expected_spawn_mass_code << "\n";
  stream << "spawned_mass_code=" << counters.spawned_mass_code << "\n";
  stream << "mass_residual_code=" << counters.mass_residual_code << "\n";
  stream << "momentum_residual_norm_code=" << counters.momentum_residual_norm_code << "\n";
  stream << "metal_mass_residual_code=" << counters.metal_mass_residual_code << "\n";
  stream << "gas_internal_energy_removed_code=" << counters.gas_internal_energy_removed_code << "\n";
  stream << "star_formation_internal_energy_sink_code="
         << counters.star_formation_internal_energy_sink_code << "\n";
  writeStarFormationCounters(stream, "latest.", counters);
  writeStarFormationCounters(
      stream, "cumulative.", cumulative_counters != nullptr ? *cumulative_counters : counters);
  if (m_effective_eos_table != nullptr) {
    stream << "effective_eos_schema_version=" << kEffectiveMultiphaseEosSchemaVersion << "\n";
    stream << "effective_eos_parameter_set=" << m_effective_eos_table->config().parameter_set << "\n";
    stream << "effective_eos_table_hash=" << m_effective_eos_table->tableHashHex() << "\n";
    stream << "effective_eos_table_bins=" << m_effective_eos_table->entries().size() << "\n";
    stream << "effective_eos_density_threshold_code=" << m_effective_eos_table->thresholdDensityPhysCode() << "\n";
    stream << "effective_eos_cooling_reference=" << m_effective_eos_table->coolingReferenceDescription() << "\n";
  }

  const std::string text = stream.str();
  core::ModuleSidecarBlock block;
  block.module_name = "star_formation";
  block.schema_version = m_config.metadata_schema_version;
  block.payload.resize(text.size());
  for (std::size_t i = 0; i < text.size(); ++i) {
    block.payload[i] = static_cast<std::byte>(text[i]);
  }
  return block;
}

StarFormationConfig makeStarFormationConfig(const core::PhysicsConfig& physics_config) {
  StarFormationConfig config;
  config.enabled = physics_config.enable_star_formation;
  config.model = physics_config.star_formation_model;
  config.density_threshold_code = physics_config.sf_density_threshold_code;
  config.temperature_threshold_k = physics_config.sf_temperature_threshold_k;
  config.min_converging_flow_rate_code = physics_config.sf_min_converging_flow_rate_code;
  config.epsilon_ff = physics_config.sf_epsilon_ff;
  config.bound_alpha_vir_max = physics_config.sf_bound_alpha_vir_max;
  config.require_converging_flow = physics_config.sf_require_converging_flow;
  config.collapse_timescale = physics_config.sf_collapse_timescale;
  config.jeans_mass_floor_code = physics_config.sf_jeans_mass_floor_code;
  config.star_particle_mass_policy = physics_config.sf_star_particle_mass_policy;
  config.target_star_particle_mass_code = physics_config.sf_target_star_particle_mass_code;
  config.target_star_particle_mass_fraction = physics_config.sf_target_star_particle_mass_fraction;
  config.min_star_particle_mass_code = physics_config.sf_min_star_particle_mass_code;
  config.max_star_particle_mass_code = physics_config.sf_max_star_particle_mass_code;
  config.max_spawn_particles_per_cell_step = physics_config.sf_max_spawn_particles_per_cell_step;
  config.max_fractional_mass_conversion = physics_config.sf_max_fractional_mass_conversion;
  config.min_remaining_gas_fraction = physics_config.sf_min_remaining_gas_fraction;
  config.min_remaining_gas_mass_code = physics_config.sf_min_remaining_gas_mass_code;
  config.temperature_safety_ceiling_k = physics_config.sf_temperature_safety_ceiling_k;
  config.stochastic_spawning = physics_config.sf_stochastic_spawning;
  config.random_seed = physics_config.sf_random_seed;
  config.effective_min_baryon_overdensity = physics_config.sf_effective_min_baryon_overdensity;
  config.effective_hot_excess_tolerance = physics_config.sf_effective_hot_excess_tolerance;
  config.effective_massive_star_fraction = physics_config.sf_effective_massive_star_fraction;
  config.effective_birth_mass_convention = physics_config.sf_effective_birth_mass_convention;
  return config;
}

StarFormationCallback::StarFormationCallback(StarFormationModel model, std::uint32_t owning_rank)
    : m_model(std::move(model)), m_owning_rank(owning_rank) {}

std::string_view StarFormationCallback::callbackName() const { return "star_formation_callback"; }

std::span<const core::IntegrationStage> StarFormationCallback::integrationStages() const {
  static constexpr std::array stages{core::IntegrationStage::kSourceTerms};
  return stages;
}

std::span<const core::StageContract> StarFormationCallback::stageContracts() const {
  static constexpr std::array contracts{core::StageContract{
      .stage = core::IntegrationStage::kSourceTerms,
      .required_inputs = core::StageDataDomain::kParticles | core::StageDataDomain::kGasCells,
      .mutated_state = core::StageDataDomain::kParticles | core::StageDataDomain::kGasCells |
          core::StageDataDomain::kDiagnostics,
      .produced_outputs = core::StageDataDomain::kParticles | core::StageDataDomain::kGasCells |
          core::StageDataDomain::kDiagnostics,
      .allowed_side_effects = core::StageDataDomain::kDiagnostics,
      .sync_requirements = core::StageSyncRequirement::kLocalOnly,
      .active_set_family = core::StageActiveSetFamily::kGasCells,
      .restart_safety = core::StageSafety::kUnsafe,
      .output_safety = core::StageSafety::kUnsafe,
      .owner = core::StageSubsystem::kSources,
  }};
  return contracts;
}

void StarFormationCallback::onStage(core::StepContext& context) {
  if (context.stage != core::IntegrationStage::kSourceTerms) {
    throw std::logic_error("star formation handler received an unregistered stage");
  }
  const std::size_t cell_count = context.state.cells.size();
  if (cell_count == 0U) {
    m_last_step_report = {};
    return;
  }
  ensureFieldSizes(cell_count);
  std::span<const std::uint32_t> active_cells = context.active_set.cell_indices;
  if (!context.active_set.cells_are_subset && active_cells.empty()) {
    m_full_cell_indices.resize(cell_count);
    std::iota(m_full_cell_indices.begin(), m_full_cell_indices.end(), 0U);
    active_cells = m_full_cell_indices;
  }
  m_last_step_report = m_model.apply(
      context.state,
      active_cells,
      m_velocity_divergence_code,
      m_metallicity_mass_fraction,
      context.integrator_state.dt_time_code,
      context.integrator_state.current_scale_factor,
      context.integrator_state.step_index,
      m_owning_rank);
}

void StarFormationCallback::setVelocityDivergenceCode(
    std::span<const double> velocity_divergence_code) {
  m_velocity_divergence_code.assign(
      velocity_divergence_code.begin(), velocity_divergence_code.end());
}

void StarFormationCallback::setMetallicityMassFraction(
    std::span<const double> metallicity_mass_fraction) {
  m_metallicity_mass_fraction.assign(
      metallicity_mass_fraction.begin(), metallicity_mass_fraction.end());
}

void StarFormationCallback::setRankLocalSeedOffset(std::uint32_t owning_rank) {
  m_owning_rank = owning_rank;
}

const StarFormationStepReport& StarFormationCallback::lastStepReport() const noexcept {
  return m_last_step_report;
}

void StarFormationCallback::ensureFieldSizes(std::size_t cell_count) {
  if (m_model.config().model != core::StarFormationModelKind::kLegacySchmidtThreshold &&
      (m_velocity_divergence_code.size() < cell_count ||
       m_metallicity_mass_fraction.size() < cell_count)) {
    throw std::runtime_error(
        "StarFormationCallback: adaptive/effective models require explicitly injected authoritative fields; use SourceRuntime in production");
  }
  if (m_velocity_divergence_code.size() < cell_count) {
    m_velocity_divergence_code.resize(cell_count, 0.0);
  }
  if (m_metallicity_mass_fraction.size() < cell_count) {
    m_metallicity_mass_fraction.resize(cell_count, 0.0);
  }
}

}  // namespace cosmosim::physics
