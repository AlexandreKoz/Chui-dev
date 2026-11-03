#include "cosmosim/physics/stellar_feedback.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace cosmosim::physics {
namespace {

constexpr double k_mass_floor = 1.0e-20;
constexpr double k_energy_floor = 1.0e-30;
constexpr double k_distance_floor = 1.0e-12;
constexpr double k_u01_norm =
    1.0 / static_cast<double>(std::numeric_limits<std::uint64_t>::max());

struct WeightedCellDistance {
  std::uint32_t cell_index = 0;
  std::uint64_t gas_cell_id = 0;
  double distance2 = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
};

[[nodiscard]] std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ull;
  x = (x ^ (x >> 30u)) * 0xbf58476d1ce4e5b9ull;
  x = (x ^ (x >> 27u)) * 0x94d049bb133111ebull;
  return x ^ (x >> 31u);
}

[[nodiscard]] double uniform01(
    std::uint64_t seed, std::uint32_t star_index, std::uint64_t step_seed) {
  std::uint64_t mixed = seed;
  mixed ^= splitmix64(step_seed + 0x9e3779b97f4a7c15ull);
  mixed ^= splitmix64(
      static_cast<std::uint64_t>(star_index) + 0xbf58476d1ce4e5b9ull);
  return static_cast<double>(splitmix64(mixed)) * k_u01_norm;
}

[[nodiscard]] bool deterministicDistanceOrder(
    const WeightedCellDistance& lhs,
    const WeightedCellDistance& rhs) noexcept {
  if (lhs.distance2 != rhs.distance2) {
    return lhs.distance2 < rhs.distance2;
  }
  if (lhs.gas_cell_id != rhs.gas_cell_id) {
    return lhs.gas_cell_id < rhs.gas_cell_id;
  }
  return lhs.cell_index < rhs.cell_index;
}

void insertBoundedNearest(
    std::vector<WeightedCellDistance>& nearest,
    WeightedCellDistance candidate,
    std::size_t keep_count) {
  if (keep_count == 0U) {
    return;
  }
  const auto position = std::lower_bound(
      nearest.begin(), nearest.end(), candidate, deterministicDistanceOrder);
  if (nearest.size() < keep_count) {
    nearest.insert(position, std::move(candidate));
    return;
  }
  if (position == nearest.end()) {
    return;
  }
  nearest.insert(position, std::move(candidate));
  nearest.pop_back();
}

}  // namespace

void StellarFeedbackSpatialIndex::rebuild(
    const StellarFeedbackGeometryView& geometry_view) {
  const std::size_t cell_count = geometry_view.cell_center_x_comoving.size();
  if (geometry_view.cell_center_y_comoving.size() != cell_count ||
      geometry_view.cell_center_z_comoving.size() != cell_count) {
    throw std::invalid_argument(
        "stellar-feedback spatial index received inconsistent geometry");
  }
  if (!geometry_view.gas_cell_id.empty() &&
      geometry_view.gas_cell_id.size() != cell_count) {
    throw std::invalid_argument(
        "stellar-feedback spatial index gas-cell IDs do not cover geometry");
  }
  if (!geometry_view.is_owned_leaf.empty() &&
      geometry_view.is_owned_leaf.size() != cell_count) {
    throw std::invalid_argument(
        "stellar-feedback spatial index ownership mask does not cover geometry");
  }

  std::vector<std::uint32_t> candidates;
  const std::size_t candidate_count = geometry_view.candidate_cell_indices.empty()
      ? cell_count
      : geometry_view.candidate_cell_indices.size();
  candidates.reserve(candidate_count);
  const auto append = [&](std::uint32_t cell_index) {
    if (cell_index >= cell_count ||
        (!geometry_view.is_owned_leaf.empty() &&
         geometry_view.is_owned_leaf[cell_index] == 0U)) {
      return;
    }
    const double x = geometry_view.cell_center_x_comoving[cell_index];
    const double y = geometry_view.cell_center_y_comoving[cell_index];
    const double z = geometry_view.cell_center_z_comoving[cell_index];
    if (std::isfinite(x) && std::isfinite(y) && std::isfinite(z)) {
      candidates.push_back(cell_index);
    }
  };
  if (geometry_view.candidate_cell_indices.empty()) {
    for (std::uint32_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      append(cell_index);
    }
  } else {
    for (const std::uint32_t cell_index : geometry_view.candidate_cell_indices) {
      append(cell_index);
    }
  }

  rebuildFromCellIndices(geometry_view, std::move(candidates));
}

void StellarFeedbackSpatialIndex::rebuildFromCellIndices(
    const StellarFeedbackGeometryView& geometry_view,
    std::vector<std::uint32_t> cell_indices) {
  const std::size_t cell_count = geometry_view.cell_center_x_comoving.size();
  if (geometry_view.cell_center_y_comoving.size() != cell_count ||
      geometry_view.cell_center_z_comoving.size() != cell_count ||
      (!geometry_view.gas_cell_id.empty() &&
       geometry_view.gas_cell_id.size() != cell_count)) {
    throw std::invalid_argument(
        "stellar-feedback spatial index received inconsistent geometry");
  }
  for (const std::uint32_t cell_index : cell_indices) {
    if (cell_index >= cell_count) {
      throw std::invalid_argument(
          "stellar-feedback spatial index candidate row is out of range");
    }
  }
  m_sorted_cell_indices = std::move(cell_indices);

  std::sort(
      m_sorted_cell_indices.begin(), m_sorted_cell_indices.end(),
      [&](std::uint32_t lhs, std::uint32_t rhs) {
        const double lhs_x = geometry_view.cell_center_x_comoving[lhs];
        const double rhs_x = geometry_view.cell_center_x_comoving[rhs];
        if (lhs_x != rhs_x) {
          return lhs_x < rhs_x;
        }
        const std::uint64_t lhs_id = geometry_view.gas_cell_id.empty()
            ? static_cast<std::uint64_t>(lhs) + 1U
            : geometry_view.gas_cell_id[lhs];
        const std::uint64_t rhs_id = geometry_view.gas_cell_id.empty()
            ? static_cast<std::uint64_t>(rhs) + 1U
            : geometry_view.gas_cell_id[rhs];
        return lhs_id != rhs_id ? lhs_id < rhs_id : lhs < rhs;
      });
  m_high_water_bytes = std::max(m_high_water_bytes, ownedCapacityBytes());
}

void StellarFeedbackSpatialIndex::clear() noexcept {
  m_sorted_cell_indices.clear();
}

std::span<const std::uint32_t>
StellarFeedbackSpatialIndex::sortedCellIndices() const noexcept {
  return m_sorted_cell_indices;
}

std::uint64_t StellarFeedbackSpatialIndex::ownedCapacityBytes() const noexcept {
  return static_cast<std::uint64_t>(m_sorted_cell_indices.capacity()) *
      static_cast<std::uint64_t>(sizeof(std::uint32_t));
}

std::uint64_t StellarFeedbackSpatialIndex::highWaterBytes() const noexcept {
  return m_high_water_bytes;
}

void StellarFeedbackModuleState::ensureStarCapacity(std::size_t) noexcept {
  // Intentionally no-op. StarParticleSidecar is the single persistent carry
  // authority; M2C forbids a duplicate O(N_star) runtime mirror.
}

StellarFeedbackModel::StellarFeedbackModel(StellarFeedbackConfig config)
    : m_config(std::move(config)) {
  if (m_config.epsilon_thermal < 0.0 || m_config.epsilon_kinetic < 0.0 ||
      m_config.epsilon_momentum < 0.0 ||
      !std::isfinite(m_config.total_energy_code_per_erg) ||
      !(m_config.total_energy_code_per_erg > 0.0)) {
    throw std::invalid_argument(
        "StellarFeedbackModel: efficiencies and energy conversion are invalid");
  }
  if (m_config.sn_energy_erg_per_mass_code <= 0.0 ||
      m_config.momentum_code_per_mass_code < 0.0) {
    throw std::invalid_argument(
        "StellarFeedbackModel: budget scales must be physically positive");
  }
  if (m_config.stochastic_event_probability <= 0.0 ||
      m_config.stochastic_event_probability > 1.0) {
    throw std::invalid_argument(
        "StellarFeedbackModel: stochastic_event_probability must be in (0,1]");
  }
  if (m_config.neighbor_count == 0) {
    throw std::invalid_argument("StellarFeedbackModel: neighbor_count must be > 0");
  }
}

const StellarFeedbackConfig& StellarFeedbackModel::config() const noexcept {
  return m_config;
}

StellarFeedbackBudget StellarFeedbackModel::computeBudget(
    double source_mass_code,
    double returned_mass_code,
    double returned_metals_code) const {
  return computeBudgetFromEnergy(
      source_mass_code,
      returned_mass_code,
      returned_metals_code,
      std::max(source_mass_code, 0.0) * m_config.sn_energy_erg_per_mass_code);
}

StellarFeedbackBudget StellarFeedbackModel::computeBudgetFromEnergy(
    double source_mass_code,
    double returned_mass_code,
    double returned_metals_code,
    double feedback_energy_erg) const {
  if (!std::isfinite(source_mass_code) || !std::isfinite(returned_mass_code) ||
      !std::isfinite(returned_metals_code) || !std::isfinite(feedback_energy_erg) ||
      source_mass_code < 0.0 || returned_mass_code < 0.0 ||
      returned_metals_code < 0.0 || feedback_energy_erg < 0.0 ||
      returned_metals_code > returned_mass_code +
          64.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, returned_mass_code)) {
    throw std::invalid_argument(
        "StellarFeedbackModel: invalid returned mass, metal, or energy budget");
  }
  StellarFeedbackBudget budget;
  budget.source_mass_code = source_mass_code;
  budget.returned_mass_code = returned_mass_code;
  budget.returned_metals_code = returned_metals_code;
  budget.total_energy_erg = feedback_energy_erg;
  budget.thermal_energy_erg = feedback_energy_erg * m_config.epsilon_thermal;
  budget.kinetic_energy_erg = feedback_energy_erg * m_config.epsilon_kinetic;
  budget.momentum_budget_code = source_mass_code *
      m_config.momentum_code_per_mass_code * m_config.epsilon_momentum;

  if (m_config.mode == StellarFeedbackMode::kThermal) {
    budget.kinetic_energy_erg = 0.0;
    budget.momentum_budget_code = 0.0;
  } else if (m_config.mode == StellarFeedbackMode::kKinetic) {
    budget.thermal_energy_erg = 0.0;
    budget.momentum_budget_code = 0.0;
  } else if (m_config.mode == StellarFeedbackMode::kMomentum) {
    budget.thermal_energy_erg = 0.0;
    budget.kinetic_energy_erg = 0.0;
  }
  return budget;
}

std::vector<StellarFeedbackTarget> StellarFeedbackModel::selectTargets(
    const StellarFeedbackGeometryView& geometry_view,
    std::uint32_t particle_index) const {
  const std::size_t cell_count = geometry_view.cell_center_x_comoving.size();
  if (particle_index >= geometry_view.particle_position_x_comoving.size() ||
      geometry_view.particle_position_y_comoving.size() <= particle_index ||
      geometry_view.particle_position_z_comoving.size() <= particle_index ||
      geometry_view.cell_center_y_comoving.size() != cell_count ||
      geometry_view.cell_center_z_comoving.size() != cell_count || cell_count == 0U) {
    return {};
  }
  if (!geometry_view.gas_cell_id.empty() &&
      geometry_view.gas_cell_id.size() != cell_count) {
    throw std::invalid_argument("stellar-feedback gas-cell IDs do not cover geometry");
  }
  if (!geometry_view.is_owned_leaf.empty() &&
      geometry_view.is_owned_leaf.size() != cell_count) {
    throw std::invalid_argument("stellar-feedback ownership mask does not cover geometry");
  }

  const double px = geometry_view.particle_position_x_comoving[particle_index];
  const double py = geometry_view.particle_position_y_comoving[particle_index];
  const double pz = geometry_view.particle_position_z_comoving[particle_index];
  const std::size_t candidate_count = geometry_view.candidate_cell_indices.empty()
      ? cell_count : geometry_view.candidate_cell_indices.size();
  const std::size_t keep_count =
      std::min<std::size_t>(m_config.neighbor_count, candidate_count);
  std::vector<WeightedCellDistance> distances;
  distances.reserve(keep_count);

  const auto append_candidate = [&](std::uint32_t cell_index) {
    if (cell_index >= cell_count ||
        (!geometry_view.is_owned_leaf.empty() &&
         geometry_view.is_owned_leaf[cell_index] == 0U)) {
      return;
    }
    const double dx = geometry_view.cell_center_x_comoving[cell_index] - px;
    const double dy = geometry_view.cell_center_y_comoving[cell_index] - py;
    const double dz = geometry_view.cell_center_z_comoving[cell_index] - pz;
    const double d2 = dx * dx + dy * dy + dz * dz;
    if (!std::isfinite(d2)) {
      return;
    }
    insertBoundedNearest(distances, WeightedCellDistance{
        .cell_index = cell_index,
        .gas_cell_id = geometry_view.gas_cell_id.empty()
            ? static_cast<std::uint64_t>(cell_index) + 1U
            : geometry_view.gas_cell_id[cell_index],
        .distance2 = d2,
        .dx = dx,
        .dy = dy,
        .dz = dz,
    }, keep_count);
  };
  if (geometry_view.candidate_cell_indices.empty()) {
    for (std::uint32_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      append_candidate(cell_index);
    }
  } else {
    for (const std::uint32_t cell_index : geometry_view.candidate_cell_indices) {
      append_candidate(cell_index);
    }
  }

  if (distances.empty()) {
    return {};
  }

  double weight_sum = 0.0;
  for (const WeightedCellDistance& distance : distances) {
    weight_sum += 1.0 /
        std::sqrt(std::max(distance.distance2, k_distance_floor));
  }
  if (!(weight_sum > 0.0) || !std::isfinite(weight_sum)) {
    return {};
  }

  std::vector<StellarFeedbackTarget> targets;
  targets.reserve(distances.size());
  for (const WeightedCellDistance& distance : distances) {
    const double inv_r = 1.0 /
        std::sqrt(std::max(distance.distance2, k_distance_floor));
    targets.push_back(StellarFeedbackTarget{
        .cell_index = distance.cell_index,
        .weight = inv_r / weight_sum,
        .radial_dx_comoving = distance.dx,
        .radial_dy_comoving = distance.dy,
        .radial_dz_comoving = distance.dz,
    });
  }
  return targets;
}

std::vector<StellarFeedbackTarget> StellarFeedbackModel::selectTargets(
    const StellarFeedbackGeometryView& geometry_view,
    const StellarFeedbackSpatialIndex& spatial_index,
    std::uint32_t particle_index) const {
  const std::size_t cell_count = geometry_view.cell_center_x_comoving.size();
  if (particle_index >= geometry_view.particle_position_x_comoving.size() ||
      geometry_view.particle_position_y_comoving.size() <= particle_index ||
      geometry_view.particle_position_z_comoving.size() <= particle_index ||
      geometry_view.cell_center_y_comoving.size() != cell_count ||
      geometry_view.cell_center_z_comoving.size() != cell_count ||
      spatial_index.sortedCellIndices().empty()) {
    return {};
  }
  const double px = geometry_view.particle_position_x_comoving[particle_index];
  const double py = geometry_view.particle_position_y_comoving[particle_index];
  const double pz = geometry_view.particle_position_z_comoving[particle_index];
  const auto sorted = spatial_index.sortedCellIndices();
  const auto right_begin = std::lower_bound(
      sorted.begin(), sorted.end(), px,
      [&](std::uint32_t cell_index, double x) {
        return geometry_view.cell_center_x_comoving[cell_index] < x;
      });
  std::size_t right = static_cast<std::size_t>(right_begin - sorted.begin());
  std::size_t left = right;
  const std::size_t keep_count =
      std::min<std::size_t>(m_config.neighbor_count, sorted.size());
  std::vector<WeightedCellDistance> nearest;
  nearest.reserve(keep_count);

  const auto x_distance2 = [&](std::uint32_t cell_index) {
    const double dx = geometry_view.cell_center_x_comoving[cell_index] - px;
    return dx * dx;
  };
  const auto evaluate = [&](std::uint32_t cell_index) {
    if (cell_index >= cell_count) {
      throw std::runtime_error("stellar-feedback spatial index contains an invalid cell row");
    }
    const double dx = geometry_view.cell_center_x_comoving[cell_index] - px;
    const double dy = geometry_view.cell_center_y_comoving[cell_index] - py;
    const double dz = geometry_view.cell_center_z_comoving[cell_index] - pz;
    const double d2 = dx * dx + dy * dy + dz * dz;
    if (!std::isfinite(d2)) {
      return;
    }
    insertBoundedNearest(nearest, WeightedCellDistance{
        .cell_index = cell_index,
        .gas_cell_id = geometry_view.gas_cell_id.empty()
            ? static_cast<std::uint64_t>(cell_index) + 1U
            : geometry_view.gas_cell_id[cell_index],
        .distance2 = d2,
        .dx = dx,
        .dy = dy,
        .dz = dz,
    }, keep_count);
  };

  while (left > 0U || right < sorted.size()) {
    const double left_dx2 = left > 0U
        ? x_distance2(sorted[left - 1U])
        : std::numeric_limits<double>::infinity();
    const double right_dx2 = right < sorted.size()
        ? x_distance2(sorted[right])
        : std::numeric_limits<double>::infinity();
    if (nearest.size() == keep_count) {
      const double worst_distance2 = nearest.back().distance2;
      if (left_dx2 > worst_distance2 && right_dx2 > worst_distance2) {
        break;
      }
    }
    if (left_dx2 <= right_dx2 && left > 0U) {
      --left;
      evaluate(sorted[left]);
    } else if (right < sorted.size()) {
      evaluate(sorted[right]);
      ++right;
    } else {
      break;
    }
  }

  if (nearest.empty()) {
    return {};
  }
  double weight_sum = 0.0;
  for (const WeightedCellDistance& distance : nearest) {
    weight_sum += 1.0 /
        std::sqrt(std::max(distance.distance2, k_distance_floor));
  }
  if (!(weight_sum > 0.0) || !std::isfinite(weight_sum)) {
    return {};
  }
  std::vector<StellarFeedbackTarget> targets;
  targets.reserve(nearest.size());
  for (const WeightedCellDistance& distance : nearest) {
    const double inv_r = 1.0 /
        std::sqrt(std::max(distance.distance2, k_distance_floor));
    targets.push_back(StellarFeedbackTarget{
        .cell_index = distance.cell_index,
        .weight = inv_r / weight_sum,
        .radial_dx_comoving = distance.dx,
        .radial_dy_comoving = distance.dy,
        .radial_dz_comoving = distance.dz,
    });
  }
  return targets;
}

std::vector<StellarFeedbackTarget> StellarFeedbackModel::selectTargets(
    const core::SimulationState& state,
    std::uint32_t particle_index) const {
  const StellarFeedbackGeometryView geometry_view{
      .particle_position_x_comoving = state.particles.position_x_comoving,
      .particle_position_y_comoving = state.particles.position_y_comoving,
      .particle_position_z_comoving = state.particles.position_z_comoving,
      .cell_center_x_comoving = state.cells.center_x_comoving,
      .cell_center_y_comoving = state.cells.center_y_comoving,
      .cell_center_z_comoving = state.cells.center_z_comoving,
      .gas_cell_id = state.gas_cells.gas_cell_id,
  };
  return selectTargets(geometry_view, particle_index);
}

StellarFeedbackStepReport StellarFeedbackModel::applyWithViews(
    core::SimulationState& state,
    StellarFeedbackModuleState& module_state,
    const StellarFeedbackGeometryView& geometry_view,
    StellarFeedbackDepositionView deposition_view,
    std::span<const std::uint32_t> active_star_indices,
    std::span<const double> returned_mass_delta_code,
    std::span<const double> returned_metals_delta_code,
    double dt_code,
    std::span<const double> feedback_energy_delta_erg) const {
  std::vector<StellarFeedbackEvent> events;
  events.reserve(active_star_indices.size());
  for (const std::uint32_t star_index : active_star_indices) {
    const double returned_mass = star_index < returned_mass_delta_code.size()
        ? returned_mass_delta_code[star_index]
        : 0.0;
    const double returned_metals = star_index < returned_metals_delta_code.size()
        ? returned_metals_delta_code[star_index]
        : 0.0;
    double energy = 0.0;
    if (star_index < feedback_energy_delta_erg.size()) {
      energy = feedback_energy_delta_erg[star_index];
    } else if (star_index < state.star_particles.size()) {
      const double source_mass = m_config.use_returned_mass_budget
          ? returned_mass
          : std::max(
                state.star_particles.birth_mass_code[star_index] * dt_code,
                0.0);
      energy = std::max(source_mass, 0.0) * m_config.sn_energy_erg_per_mass_code;
    }
    events.push_back(StellarFeedbackEvent{
        .star_index = star_index,
        .returned_mass_code = returned_mass,
        .returned_metals_code = returned_metals,
        .feedback_energy_erg = energy,
    });
  }
  return applyEventsWithViews(
      state,
      module_state,
      geometry_view,
      nullptr,
      deposition_view,
      events,
      dt_code);
}

StellarFeedbackStepReport StellarFeedbackModel::applyEventsWithViews(
    core::SimulationState& state,
    StellarFeedbackModuleState& module_state,
    const StellarFeedbackGeometryView& geometry_view,
    const StellarFeedbackSpatialIndex* spatial_index,
    StellarFeedbackDepositionView deposition_view,
    std::span<const StellarFeedbackEvent> events,
    double dt_code) const {
  StellarFeedbackStepReport report;
  const std::size_t cell_count = geometry_view.cell_center_x_comoving.size();
  if (geometry_view.cell_center_y_comoving.size() != cell_count ||
      geometry_view.cell_center_z_comoving.size() != cell_count ||
      deposition_view.cell_mass_code.size() < cell_count ||
      deposition_view.gas_density_code.size() < cell_count ||
      deposition_view.gas_internal_energy_code.size() < cell_count ||
      deposition_view.gas_metal_mass_code.size() < cell_count ||
      (!deposition_view.cell_volume_code.empty() &&
       deposition_view.cell_volume_code.size() < cell_count)) {
    throw std::invalid_argument(
        "stellar-feedback deposition views do not cover every gas cell");
  }
  if (!m_config.enabled || dt_code <= 0.0 || state.star_particles.size() == 0U) {
    state.sidecars.upsert(buildMetadataSidecar(report));
    return report;
  }
  if (!state.star_particles.isConsistent() ||
      state.star_particles.enrichment_carry_mass_code.size() !=
          state.star_particles.size()) {
    throw std::runtime_error(
        "stellar-feedback persistent star sidecar is inconsistent");
  }

  module_state.ensureStarCapacity(state.star_particles.size());
  const std::uint64_t step_seed = state.metadata.step_index;

  for (const StellarFeedbackEvent& event : events) {
    const std::uint32_t star_index = event.star_index;
    ++report.counters.scanned_stars;
    if (star_index >= state.star_particles.size()) {
      continue;
    }
    const std::uint32_t particle_index =
        state.star_particles.particle_index[star_index];
    if (particle_index >= state.particles.size()) {
      continue;
    }

    const double returned_mass = event.returned_mass_code;
    const double returned_metals = event.returned_metals_code;
    const double source_mass = m_config.use_returned_mass_budget
        ? returned_mass
        : std::max(state.star_particles.birth_mass_code[star_index] * dt_code, 0.0);
    const double energy = event.feedback_energy_erg;

    StellarFeedbackStarReport star_report;
    star_report.star_index = star_index;
    star_report.particle_index = particle_index;
    star_report.budget = computeBudgetFromEnergy(
        source_mass, returned_mass, returned_metals, energy);

    const double carry_mass =
        state.star_particles.enrichment_carry_mass_code[star_index];
    const double carry_metals =
        state.star_particles.enrichment_carry_metals_code[star_index];
    const double carry_energy =
        state.star_particles.enrichment_carry_feedback_energy_erg[star_index];
    const double carry_momentum =
        state.star_particles.enrichment_carry_momentum_code[star_index];

    star_report.budget.returned_mass_code += carry_mass;
    star_report.budget.returned_metals_code += carry_metals;
    const double thermal_fraction = star_report.budget.total_energy_erg > 0.0
        ? star_report.budget.thermal_energy_erg /
              star_report.budget.total_energy_erg : 1.0;
    star_report.budget.total_energy_erg += carry_energy;
    star_report.budget.thermal_energy_erg += carry_energy * thermal_fraction;
    star_report.budget.kinetic_energy_erg += carry_energy * (1.0 - thermal_fraction);
    star_report.budget.momentum_budget_code += carry_momentum;
    if (star_report.budget.returned_metals_code >
        star_report.budget.returned_mass_code +
            64.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, star_report.budget.returned_mass_code)) {
      throw std::runtime_error(
          "stellar-feedback carried metals exceed carried returned mass");
    }

    state.star_particles.enrichment_carry_mass_code[star_index] = 0.0;
    state.star_particles.enrichment_carry_metals_code[star_index] = 0.0;
    state.star_particles.enrichment_carry_feedback_energy_erg[star_index] = 0.0;
    state.star_particles.enrichment_carry_momentum_code[star_index] = 0.0;

    if (star_report.budget.returned_mass_code <= k_mass_floor &&
        star_report.budget.thermal_energy_erg <= k_energy_floor &&
        star_report.budget.kinetic_energy_erg <= k_energy_floor &&
        star_report.budget.momentum_budget_code <= k_mass_floor) {
      continue;
    }

    star_report.stochastic_event_fired =
        m_config.variant != StellarFeedbackVariant::kStochastic ||
        stochasticEventFires(star_index, step_seed);
    std::vector<StellarFeedbackTarget> targets = spatial_index != nullptr
        ? selectTargets(geometry_view, *spatial_index, particle_index)
        : selectTargets(geometry_view, particle_index);
    star_report.target_count = targets.size();
    report.counters.target_cells_visited += targets.size();

    if (!star_report.stochastic_event_fired || targets.empty()) {
      star_report.unresolved_mass_code = star_report.budget.returned_mass_code;
      star_report.unresolved_metals_code = star_report.budget.returned_metals_code;
      star_report.unresolved_thermal_energy_erg =
          star_report.budget.thermal_energy_erg;
      star_report.unresolved_kinetic_energy_erg =
          star_report.budget.kinetic_energy_erg;
      star_report.unresolved_momentum_code =
          star_report.budget.momentum_budget_code;
    } else {
      // Momentum injection remains durably carried until a velocity-conservative
      // deposition contract is available. It is never discarded or misapplied to u.
      star_report.unresolved_momentum_code =
          star_report.budget.momentum_budget_code;
      for (const StellarFeedbackTarget& target : targets) {
        const std::uint32_t cell_index = target.cell_index;
        const double weight = target.weight;
        if (cell_index >= cell_count || !(weight >= 0.0) || !std::isfinite(weight)) {
          throw std::runtime_error("stellar-feedback produced an invalid target");
        }

        const double mass_add =
            star_report.budget.returned_mass_code * weight;
        const double metals_add =
            star_report.budget.returned_metals_code * weight;
        const double thermal_add =
            star_report.budget.thermal_energy_erg * weight;
        const double kinetic_add =
            star_report.budget.kinetic_energy_erg * weight;
        const double old_mass = deposition_view.cell_mass_code[cell_index];
        const double old_density = deposition_view.gas_density_code[cell_index];
        const double old_u = deposition_view.gas_internal_energy_code[cell_index];
        if (deposition_view.thermodynamic_touch_sink != nullptr) {
          // Audit F4: report the mutation before it happens so the caller
          // can refresh derived P / c_s / T before the next CFL decision.
          deposition_view.thermodynamic_touch_sink->recordGasCellMutation(cell_index, old_u);
        }
        double volume = 0.0;
        if (!deposition_view.cell_volume_code.empty()) {
          volume = deposition_view.cell_volume_code[cell_index];
        } else if (deposition_view.cell_volume_provider != nullptr) {
          volume = deposition_view.cell_volume_provider->cellVolumeCode(cell_index);
        } else {
          volume = old_density > 0.0 ? old_mass / old_density : 0.0;
        }
        if (!(volume > 0.0) || !std::isfinite(volume) || !(old_mass >= 0.0) ||
            !std::isfinite(old_u)) {
          throw std::runtime_error(
              "stellar-feedback target lacks a valid cell volume or energy state");
        }
        const double new_mass = old_mass + mass_add;
        const double old_internal_total_code = old_mass * old_u;

        deposition_view.cell_mass_code[cell_index] = new_mass;
        deposition_view.gas_density_code[cell_index] = new_mass / volume;
        deposition_view.gas_metal_mass_code[cell_index] += metals_add;
        if (deposition_view.gas_metal_mass_code[cell_index] > new_mass +
            64.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, new_mass)) {
          throw std::runtime_error(
              "stellar-feedback deposition would create metal mass above gas mass");
        }

        if (m_config.variant == StellarFeedbackVariant::kDelayedCooling) {
          star_report.delayed_cooling_applied = true;
          star_report.unresolved_thermal_energy_erg += thermal_add;
        } else {
          star_report.deposited_thermal_energy_erg += thermal_add;
        }
        // Until a total-energy/momentum injection operator is present, the named
        // kinetic energy channel is thermalized conservatively in total energy.
        star_report.deposited_kinetic_energy_erg += kinetic_add;
        const double energy_deposited_erg =
            (m_config.variant == StellarFeedbackVariant::kDelayedCooling
                 ? kinetic_add
                 : thermal_add + kinetic_add);
        const double new_internal_total_code = old_internal_total_code +
            energy_deposited_erg * m_config.total_energy_code_per_erg;
        deposition_view.gas_internal_energy_code[cell_index] = new_mass > 0.0
            ? new_internal_total_code / new_mass : 0.0;

        star_report.deposited_mass_code += mass_add;
        star_report.deposited_metals_code += metals_add;
      }
    }

    const double unresolved_energy =
        star_report.unresolved_thermal_energy_erg +
        star_report.unresolved_kinetic_energy_erg;
    state.star_particles.enrichment_carry_mass_code[star_index] =
        star_report.unresolved_mass_code;
    state.star_particles.enrichment_carry_metals_code[star_index] =
        star_report.unresolved_metals_code;
    state.star_particles.enrichment_carry_feedback_energy_erg[star_index] =
        unresolved_energy;
    state.star_particles.enrichment_carry_momentum_code[star_index] =
        star_report.unresolved_momentum_code;
    state.star_particles.stellar_deposited_mass_cumulative_code[star_index] +=
        star_report.deposited_mass_code;
    state.star_particles.stellar_deposited_metals_cumulative_code[star_index] +=
        star_report.deposited_metals_code;
    state.star_particles.stellar_deposited_feedback_energy_cumulative_erg[star_index] +=
        star_report.deposited_thermal_energy_erg +
        star_report.deposited_kinetic_energy_erg;

    ++report.counters.feedback_stars;
    report.counters.source_mass_code += star_report.budget.source_mass_code;
    report.counters.deposited_mass_code += star_report.deposited_mass_code;
    report.counters.deposited_metals_code += star_report.deposited_metals_code;
    report.counters.deposited_thermal_energy_erg +=
        star_report.deposited_thermal_energy_erg;
    report.counters.deposited_kinetic_energy_erg +=
        star_report.deposited_kinetic_energy_erg;
    report.counters.deposited_momentum_code +=
        star_report.deposited_momentum_code;
    report.counters.unresolved_mass_code += star_report.unresolved_mass_code;
    report.counters.unresolved_metals_code += star_report.unresolved_metals_code;
    report.counters.unresolved_thermal_energy_erg +=
        star_report.unresolved_thermal_energy_erg;
    report.counters.unresolved_kinetic_energy_erg +=
        star_report.unresolved_kinetic_energy_erg;
    report.counters.unresolved_momentum_code +=
        star_report.unresolved_momentum_code;
    report.star_reports.push_back(star_report);
  }

  state.sidecars.upsert(buildMetadataSidecar(report));
  return report;
}

void StellarFeedbackModel::applyEventsCountersOnly(
    core::SimulationState& state,
    StellarFeedbackModuleState& module_state,
    const StellarFeedbackGeometryView& geometry_view,
    const StellarFeedbackSpatialIndex* spatial_index,
    StellarFeedbackDepositionView deposition_view,
    std::span<const StellarFeedbackEvent> events,
    double dt_code,
    StellarFeedbackStepReport& report_out) const {
  const std::size_t cell_count = geometry_view.cell_center_x_comoving.size();
  if (!m_config.enabled || dt_code <= 0.0 || state.star_particles.size() == 0U) {
    return;
  }
  if (!state.star_particles.isConsistent() ||
      state.star_particles.enrichment_carry_mass_code.size() !=
          state.star_particles.size()) {
    throw std::runtime_error(
        "stellar-feedback persistent star sidecar is inconsistent");
  }
  module_state.ensureStarCapacity(state.star_particles.size());
  const std::uint64_t step_seed = state.metadata.step_index;

  for (const StellarFeedbackEvent& event : events) {
    const std::uint32_t star_index = event.star_index;
    ++report_out.counters.scanned_stars;
    if (star_index >= state.star_particles.size()) {
      continue;
    }
    const std::uint32_t particle_index =
        state.star_particles.particle_index[star_index];
    if (particle_index >= state.particles.size()) {
      continue;
    }

    const double returned_mass = event.returned_mass_code;
    const double returned_metals = event.returned_metals_code;
    const double source_mass = m_config.use_returned_mass_budget
        ? returned_mass
        : std::max(state.star_particles.birth_mass_code[star_index] * dt_code, 0.0);
    const double energy = event.feedback_energy_erg;

    StellarFeedbackStarReport star_report;
    star_report.star_index = star_index;
    star_report.particle_index = particle_index;
    star_report.budget = computeBudgetFromEnergy(
        source_mass, returned_mass, returned_metals, energy);

    const double carry_mass =
        state.star_particles.enrichment_carry_mass_code[star_index];
    const double carry_metals =
        state.star_particles.enrichment_carry_metals_code[star_index];
    const double carry_energy =
        state.star_particles.enrichment_carry_feedback_energy_erg[star_index];
    const double carry_momentum =
        state.star_particles.enrichment_carry_momentum_code[star_index];

    star_report.budget.returned_mass_code += carry_mass;
    star_report.budget.returned_metals_code += carry_metals;
    const double thermal_fraction = star_report.budget.total_energy_erg > 0.0
        ? star_report.budget.thermal_energy_erg /
              star_report.budget.total_energy_erg : 1.0;
    star_report.budget.total_energy_erg += carry_energy;
    star_report.budget.thermal_energy_erg += carry_energy * thermal_fraction;
    star_report.budget.kinetic_energy_erg += carry_energy * (1.0 - thermal_fraction);
    star_report.budget.momentum_budget_code += carry_momentum;
    if (star_report.budget.returned_metals_code >
        star_report.budget.returned_mass_code +
            64.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, star_report.budget.returned_mass_code)) {
      throw std::runtime_error(
          "stellar-feedback carried metals exceed carried returned mass");
    }

    state.star_particles.enrichment_carry_mass_code[star_index] = 0.0;
    state.star_particles.enrichment_carry_metals_code[star_index] = 0.0;
    state.star_particles.enrichment_carry_feedback_energy_erg[star_index] = 0.0;
    state.star_particles.enrichment_carry_momentum_code[star_index] = 0.0;

    if (star_report.budget.returned_mass_code <= k_mass_floor &&
        star_report.budget.thermal_energy_erg <= k_energy_floor &&
        star_report.budget.kinetic_energy_erg <= k_energy_floor &&
        star_report.budget.momentum_budget_code <= k_mass_floor) {
      continue;
    }

    star_report.stochastic_event_fired =
        m_config.variant != StellarFeedbackVariant::kStochastic ||
        stochasticEventFires(star_index, step_seed);
    std::vector<StellarFeedbackTarget> targets = spatial_index != nullptr
        ? selectTargets(geometry_view, *spatial_index, particle_index)
        : selectTargets(geometry_view, particle_index);
    report_out.counters.target_cells_visited += targets.size();

    if (!star_report.stochastic_event_fired || targets.empty()) {
      star_report.unresolved_mass_code = star_report.budget.returned_mass_code;
      star_report.unresolved_metals_code = star_report.budget.returned_metals_code;
      star_report.unresolved_thermal_energy_erg =
          star_report.budget.thermal_energy_erg;
      star_report.unresolved_kinetic_energy_erg =
          star_report.budget.kinetic_energy_erg;
      star_report.unresolved_momentum_code =
          star_report.budget.momentum_budget_code;
    } else {
      star_report.unresolved_momentum_code =
          star_report.budget.momentum_budget_code;
      for (const StellarFeedbackTarget& target : targets) {
        const std::uint32_t cell_index = target.cell_index;
        const double weight = target.weight;
        if (cell_index >= cell_count || !(weight >= 0.0) || !std::isfinite(weight)) {
          throw std::runtime_error("stellar-feedback produced an invalid target");
        }

        const double mass_add =
            star_report.budget.returned_mass_code * weight;
        const double metals_add =
            star_report.budget.returned_metals_code * weight;
        const double thermal_add =
            star_report.budget.thermal_energy_erg * weight;
        const double kinetic_add =
            star_report.budget.kinetic_energy_erg * weight;
        const double old_mass = deposition_view.cell_mass_code[cell_index];
        const double old_density = deposition_view.gas_density_code[cell_index];
        const double old_u = deposition_view.gas_internal_energy_code[cell_index];
        if (deposition_view.thermodynamic_touch_sink != nullptr) {
          // Audit F4: report the mutation before it happens so the caller
          // can refresh derived P / c_s / T before the next CFL decision.
          deposition_view.thermodynamic_touch_sink->recordGasCellMutation(cell_index, old_u);
        }
        double volume = 0.0;
        if (!deposition_view.cell_volume_code.empty()) {
          volume = deposition_view.cell_volume_code[cell_index];
        } else if (deposition_view.cell_volume_provider != nullptr) {
          volume = deposition_view.cell_volume_provider->cellVolumeCode(cell_index);
        } else {
          volume = old_density > 0.0 ? old_mass / old_density : 0.0;
        }
        if (!(volume > 0.0) || !std::isfinite(volume) || !(old_mass >= 0.0) ||
            !std::isfinite(old_u)) {
          throw std::runtime_error(
              "stellar-feedback target lacks a valid cell volume or energy state");
        }
        const double new_mass = old_mass + mass_add;
        const double old_internal_total_code = old_mass * old_u;

        deposition_view.cell_mass_code[cell_index] = new_mass;
        deposition_view.gas_density_code[cell_index] = new_mass / volume;
        deposition_view.gas_metal_mass_code[cell_index] += metals_add;
        if (deposition_view.gas_metal_mass_code[cell_index] > new_mass +
            64.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, new_mass)) {
          throw std::runtime_error(
              "stellar-feedback deposition would create metal mass above gas mass");
        }

        if (m_config.variant == StellarFeedbackVariant::kDelayedCooling) {
          star_report.delayed_cooling_applied = true;
          star_report.unresolved_thermal_energy_erg += thermal_add;
        } else {
          star_report.deposited_thermal_energy_erg += thermal_add;
        }
        star_report.deposited_kinetic_energy_erg += kinetic_add;
        const double energy_deposited_erg =
            (m_config.variant == StellarFeedbackVariant::kDelayedCooling
                 ? kinetic_add
                 : thermal_add + kinetic_add);
        const double new_internal_total_code = old_internal_total_code +
            energy_deposited_erg * m_config.total_energy_code_per_erg;
        deposition_view.gas_internal_energy_code[cell_index] = new_mass > 0.0
            ? new_internal_total_code / new_mass : 0.0;

        star_report.deposited_mass_code += mass_add;
        star_report.deposited_metals_code += metals_add;
      }
    }

    const double unresolved_energy =
        star_report.unresolved_thermal_energy_erg +
        star_report.unresolved_kinetic_energy_erg;
    state.star_particles.enrichment_carry_mass_code[star_index] =
        star_report.unresolved_mass_code;
    state.star_particles.enrichment_carry_metals_code[star_index] =
        star_report.unresolved_metals_code;
    state.star_particles.enrichment_carry_feedback_energy_erg[star_index] =
        unresolved_energy;
    state.star_particles.enrichment_carry_momentum_code[star_index] =
        star_report.unresolved_momentum_code;
    state.star_particles.stellar_deposited_mass_cumulative_code[star_index] +=
        star_report.deposited_mass_code;
    state.star_particles.stellar_deposited_metals_cumulative_code[star_index] +=
        star_report.deposited_metals_code;
    state.star_particles.stellar_deposited_feedback_energy_cumulative_erg[star_index] +=
        star_report.deposited_thermal_energy_erg +
        star_report.deposited_kinetic_energy_erg;

    ++report_out.counters.feedback_stars;
    report_out.counters.source_mass_code += star_report.budget.source_mass_code;
    report_out.counters.deposited_mass_code += star_report.deposited_mass_code;
    report_out.counters.deposited_metals_code += star_report.deposited_metals_code;
    report_out.counters.deposited_thermal_energy_erg +=
        star_report.deposited_thermal_energy_erg;
    report_out.counters.deposited_kinetic_energy_erg +=
        star_report.deposited_kinetic_energy_erg;
    report_out.counters.deposited_momentum_code +=
        star_report.deposited_momentum_code;
    report_out.counters.unresolved_mass_code += star_report.unresolved_mass_code;
    report_out.counters.unresolved_metals_code += star_report.unresolved_metals_code;
    report_out.counters.unresolved_thermal_energy_erg +=
        star_report.unresolved_thermal_energy_erg;
    report_out.counters.unresolved_kinetic_energy_erg +=
        star_report.unresolved_kinetic_energy_erg;
    report_out.counters.unresolved_momentum_code +=
        star_report.unresolved_momentum_code;
    // No star_reports.push_back — this is the production path that avoids
    // materializing a per-star report vector.
  }

  state.sidecars.upsert(buildMetadataSidecar(report_out));
}

StellarFeedbackStepReport StellarFeedbackModel::apply(
    core::SimulationState& state,
    StellarFeedbackModuleState& module_state,
    std::span<const std::uint32_t> active_star_indices,
    std::span<const double> returned_mass_delta_code,
    std::span<const double> returned_metals_delta_code,
    double dt_code,
    std::span<const double> feedback_energy_delta_erg) const {
  const StellarFeedbackGeometryView geometry_view{
      .particle_position_x_comoving = state.particles.position_x_comoving,
      .particle_position_y_comoving = state.particles.position_y_comoving,
      .particle_position_z_comoving = state.particles.position_z_comoving,
      .cell_center_x_comoving = state.cells.center_x_comoving,
      .cell_center_y_comoving = state.cells.center_y_comoving,
      .cell_center_z_comoving = state.cells.center_z_comoving,
      .gas_cell_id = state.gas_cells.gas_cell_id,
  };
  StellarFeedbackDepositionView deposition_view{
      .cell_mass_code = state.cells.mass_code,
      .gas_density_code = state.gas_cells.density_code,
      .gas_internal_energy_code = state.gas_cells.internal_energy_code,
      .gas_metal_mass_code = state.gas_cells.metal_mass_code,
  };
  return applyWithViews(
      state, module_state, geometry_view, deposition_view,
      active_star_indices, returned_mass_delta_code,
      returned_metals_delta_code, dt_code, feedback_energy_delta_erg);
}

core::ModuleSidecarBlock StellarFeedbackModel::buildMetadataSidecar(
    const StellarFeedbackStepReport& report) const {
  std::ostringstream stream;
  stream << "module=stellar_feedback\n";
  stream << "schema_version=" << m_config.metadata_schema_version << "\n";
  stream << "mode=" << modeToString(m_config.mode) << "\n";
  stream << "variant=" << variantToString(m_config.variant) << "\n";
  stream << "carry_authority=star_particle_persistent_sidecar\n";
  stream << "kinetic_energy_policy=thermalized_total_energy\n";
  stream << "momentum_policy=persistent_carry_until_velocity_operator\n";
  stream << "total_energy_code_per_erg=" << m_config.total_energy_code_per_erg << "\n";
  stream << "scanned_stars=" << report.counters.scanned_stars << "\n";
  stream << "feedback_stars=" << report.counters.feedback_stars << "\n";
  stream << "target_cells_visited=" << report.counters.target_cells_visited << "\n";
  stream << "source_mass_code=" << report.counters.source_mass_code << "\n";
  stream << "deposited_mass_code=" << report.counters.deposited_mass_code << "\n";
  stream << "deposited_metals_code=" << report.counters.deposited_metals_code << "\n";
  stream << "deposited_thermal_energy_erg="
         << report.counters.deposited_thermal_energy_erg << "\n";
  stream << "deposited_kinetic_energy_erg="
         << report.counters.deposited_kinetic_energy_erg << "\n";
  stream << "unresolved_mass_code=" << report.counters.unresolved_mass_code << "\n";
  stream << "unresolved_metals_code=" << report.counters.unresolved_metals_code << "\n";
  stream << "unresolved_feedback_energy_erg="
         << report.counters.unresolved_thermal_energy_erg +
                report.counters.unresolved_kinetic_energy_erg << "\n";
  stream << "unresolved_momentum_code="
         << report.counters.unresolved_momentum_code << "\n";

  const std::string text = stream.str();
  core::ModuleSidecarBlock block;
  block.module_name = "stellar_feedback";
  block.schema_version = m_config.metadata_schema_version;
  block.payload.resize(text.size());
  for (std::size_t i = 0; i < text.size(); ++i) {
    block.payload[i] = static_cast<std::byte>(text[i]);
  }
  return block;
}

std::string StellarFeedbackModel::modeToString(StellarFeedbackMode mode) {
  switch (mode) {
    case StellarFeedbackMode::kThermal: return "thermal";
    case StellarFeedbackMode::kKinetic: return "kinetic";
    case StellarFeedbackMode::kMomentum: return "momentum";
    case StellarFeedbackMode::kThermalKineticMomentum:
      return "thermal_kinetic_momentum";
  }
  return "thermal_kinetic_momentum";
}

std::string StellarFeedbackModel::variantToString(
    StellarFeedbackVariant variant) {
  switch (variant) {
    case StellarFeedbackVariant::kNone: return "none";
    case StellarFeedbackVariant::kDelayedCooling: return "delayed_cooling";
    case StellarFeedbackVariant::kStochastic: return "stochastic";
  }
  return "none";
}

bool StellarFeedbackModel::stochasticEventFires(
    std::uint32_t star_index, std::uint64_t step_seed) const {
  return uniform01(m_config.random_seed, star_index, step_seed) <
      m_config.stochastic_event_probability;
}

StellarFeedbackConfig makeStellarFeedbackConfig(
    const core::PhysicsConfig& physics_config) {
  StellarFeedbackConfig config;
  config.enabled = physics_config.enable_feedback;
  config.use_returned_mass_budget = physics_config.fb_use_returned_mass_budget;
  config.epsilon_thermal = physics_config.fb_epsilon_thermal;
  config.epsilon_kinetic = physics_config.fb_epsilon_kinetic;
  config.epsilon_momentum = physics_config.fb_epsilon_momentum;
  config.sn_energy_erg_per_mass_code =
      physics_config.fb_sn_energy_erg_per_mass_code;
  config.momentum_code_per_mass_code =
      physics_config.fb_momentum_code_per_mass_code;
  config.neighbor_count = physics_config.fb_neighbor_count;
  config.delayed_cooling_time_code =
      physics_config.fb_delayed_cooling_time_code;
  config.stochastic_event_probability =
      physics_config.fb_stochastic_event_probability;
  config.random_seed = physics_config.fb_random_seed;
  switch (physics_config.fb_mode) {
    case core::FeedbackMode::kThermal:
      config.mode = StellarFeedbackMode::kThermal; break;
    case core::FeedbackMode::kKinetic:
      config.mode = StellarFeedbackMode::kKinetic; break;
    case core::FeedbackMode::kMomentum:
      config.mode = StellarFeedbackMode::kMomentum; break;
    case core::FeedbackMode::kThermalKineticMomentum:
      config.mode = StellarFeedbackMode::kThermalKineticMomentum; break;
  }
  switch (physics_config.fb_variant) {
    case core::FeedbackVariant::kDelayedCooling:
      config.variant = StellarFeedbackVariant::kDelayedCooling; break;
    case core::FeedbackVariant::kStochastic:
      config.variant = StellarFeedbackVariant::kStochastic; break;
    case core::FeedbackVariant::kNone:
      config.variant = StellarFeedbackVariant::kNone; break;
  }
  return config;
}

}  // namespace cosmosim::physics
