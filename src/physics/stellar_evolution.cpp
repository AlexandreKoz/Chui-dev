#include "cosmosim/physics/stellar_evolution.hpp"

#include "cosmosim/core/checked_arithmetic.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace cosmosim::physics {

std::uint64_t stellarEvolutionBatchStagingBytes(std::size_t batch_capacity) {
  return static_cast<std::uint64_t>(core::checkedSizeMultiply(
      batch_capacity, sizeof(StellarEvolutionStarBudget),
      "stellar-evolution governed budget staging bytes"));
}

namespace {

constexpr double k_mass_floor = 1.0e-20;
constexpr double k_seconds_per_year = 31557600.0;
constexpr double k_validation_tolerance = 1.0e-10;

[[nodiscard]] bool finiteNonnegative(double value) {
  return std::isfinite(value) && value >= 0.0;
}

[[nodiscard]] std::string trim(const std::string& input) {
  const auto begin = input.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) {
    return {};
  }
  const auto end = input.find_last_not_of(" \t\r\n");
  return input.substr(begin, end - begin + 1);
}

[[nodiscard]] double interpolateLinear(
    double x0,
    double x1,
    double y0,
    double y1,
    double x) {
  if (!(x1 > x0)) {
    return y0;
  }
  const double fraction = std::clamp((x - x0) / (x1 - x0), 0.0, 1.0);
  return y0 + fraction * (y1 - y0);
}

struct Bracket {
  std::size_t left = 0;
  std::size_t right = 0;
  double fraction = 0.0;
};

[[nodiscard]] Bracket linearBracket(std::span<const double> grid, double value) {
  if (grid.empty()) {
    return {};
  }
  const double clamped = std::clamp(value, grid.front(), grid.back());
  const auto upper = std::lower_bound(grid.begin(), grid.end(), clamped);
  if (upper == grid.begin()) {
    return {0U, 0U, 0.0};
  }
  if (upper == grid.end()) {
    const std::size_t last = grid.size() - 1U;
    return {last, last, 0.0};
  }
  const std::size_t right = static_cast<std::size_t>(upper - grid.begin());
  const std::size_t left = right - 1U;
  return {left, right, (clamped - grid[left]) / (grid[right] - grid[left])};
}

[[nodiscard]] Bracket ageBracket(std::span<const double> grid, double age_years) {
  if (grid.empty()) {
    return {};
  }
  const double clamped = std::clamp(age_years, grid.front(), grid.back());
  const auto upper = std::lower_bound(grid.begin(), grid.end(), clamped);
  if (upper == grid.begin()) {
    return {0U, 0U, 0.0};
  }
  if (upper == grid.end()) {
    const std::size_t last = grid.size() - 1U;
    return {last, last, 0.0};
  }
  const std::size_t right = static_cast<std::size_t>(upper - grid.begin());
  const std::size_t left = right - 1U;
  // Zero-age is treated explicitly; positive intervals use log-age interpolation.
  if (grid[left] <= 0.0 || clamped <= 0.0) {
    return {left, right, (clamped - grid[left]) / (grid[right] - grid[left])};
  }
  const double log_left = std::log(grid[left]);
  const double log_right = std::log(grid[right]);
  const double log_value = std::log(clamped);
  return {left, right, (log_value - log_left) / (log_right - log_left)};
}

[[nodiscard]] double bilinear(
    const StellarEvolutionTable& table,
    const std::vector<double>& values,
    const Bracket& age,
    const Bracket& metallicity) {
  const auto sample = [&](std::size_t z_index, std::size_t age_index) {
    return values[table.flatIndex(z_index, age_index)];
  };
  const double low_z = interpolateLinear(
      0.0, 1.0, sample(metallicity.left, age.left), sample(metallicity.left, age.right),
      age.fraction);
  const double high_z = interpolateLinear(
      0.0, 1.0, sample(metallicity.right, age.left), sample(metallicity.right, age.right),
      age.fraction);
  return interpolateLinear(0.0, 1.0, low_z, high_z, metallicity.fraction);
}

[[nodiscard]] const std::vector<double>& totalEjectedValues(
    const StellarEvolutionTable& table) {
  return table.total_ejected_metal_fraction_total.empty()
      ? table.metal_yield_fraction_total
      : table.total_ejected_metal_fraction_total;
}

[[nodiscard]] const std::vector<double>& totalEjectedChannelValues(
    const StellarEvolutionTable& table,
    std::size_t channel) {
  return table.total_ejected_metal_fraction_channel[channel].empty()
      ? table.metal_yield_fraction_channel[channel]
      : table.total_ejected_metal_fraction_channel[channel];
}

[[nodiscard]] bool monotonicNondecreasing(
    const StellarEvolutionTable& table,
    const std::vector<double>& values) {
  const std::size_t age_count = table.age_yr.size();
  const std::size_t z_count = table.birth_metallicity_mass_fraction.empty()
      ? 1U
      : table.birth_metallicity_mass_fraction.size();
  if (values.size() != age_count * z_count) {
    return false;
  }
  for (std::size_t z_index = 0; z_index < z_count; ++z_index) {
    for (std::size_t age_index = 0; age_index < age_count; ++age_index) {
      const double value = values[table.flatIndex(z_index, age_index)];
      if (!finiteNonnegative(value)) {
        return false;
      }
      if (age_index > 0U) {
        const double previous = values[table.flatIndex(z_index, age_index - 1U)];
        if (value + k_validation_tolerance < previous) {
          return false;
        }
      }
    }
  }
  return true;
}

[[nodiscard]] double nonnegativeDifference(double end, double begin) {
  const double difference = end - begin;
  if (difference < -k_validation_tolerance) {
    throw std::runtime_error("stellar-evolution cumulative table decreased over an interval");
  }
  return std::max(difference, 0.0);
}

void parseMetadataLine(StellarEvolutionTable& table, const std::string& metadata) {
  const std::size_t equal_pos = metadata.find('=');
  if (equal_pos == std::string::npos) {
    return;
  }
  const std::string key = trim(metadata.substr(0, equal_pos));
  const std::string value = trim(metadata.substr(equal_pos + 1U));
  if (key == "table_id") table.table_id = value;
  else if (key == "table_version") table.table_version = value;
  else if (key == "source_papers") table.source_papers = value;
  else if (key == "source_repository") table.source_repository = value;
  else if (key == "redistribution_license") table.redistribution_license = value;
  else if (key == "sha256") table.sha256 = value;
  else if (key == "imf") table.imf = value;
  else if (key == "stellar_mass_range") table.stellar_mass_range = value;
  else if (key == "solar_abundance_reference") table.solar_abundance_reference = value;
  else if (key == "production_calibrated") table.production_calibrated = value == "true";
}

}  // namespace

std::size_t StellarEvolutionTable::flatIndex(
    std::size_t metallicity_index,
    std::size_t age_index) const {
  return metallicity_index * age_yr.size() + age_index;
}

bool StellarEvolutionTable::isConsistent() const noexcept {
  try {
    if (age_yr.size() < 2U) return false;
    for (std::size_t index = 0; index < age_yr.size(); ++index) {
      if (!finiteNonnegative(age_yr[index])) return false;
      if (index > 0U && !(age_yr[index] > age_yr[index - 1U])) return false;
    }

    const std::size_t effective_z_count = birth_metallicity_mass_fraction.empty()
        ? 1U
        : birth_metallicity_mass_fraction.size();
    for (std::size_t index = 0; index < birth_metallicity_mass_fraction.size(); ++index) {
      const double metallicity = birth_metallicity_mass_fraction[index];
      if (!std::isfinite(metallicity) || metallicity < 0.0 || metallicity > 1.0) return false;
      if (index > 0U && !(metallicity > birth_metallicity_mass_fraction[index - 1U])) return false;
    }

    const std::size_t expected = age_yr.size() * effective_z_count;
    if (return_fraction_total.size() != expected ||
        totalEjectedValues(*this).size() != expected ||
        energy_erg_per_initial_mass_code.size() != expected) {
      return false;
    }
    if (!newly_synthesized_metal_fraction_total.empty() &&
        newly_synthesized_metal_fraction_total.size() != expected) return false;
    if (!event_count_per_initial_mass_code_total.empty() &&
        event_count_per_initial_mass_code_total.size() != expected) return false;

    if (!monotonicNondecreasing(*this, return_fraction_total) ||
        !monotonicNondecreasing(*this, totalEjectedValues(*this)) ||
        !monotonicNondecreasing(*this, energy_erg_per_initial_mass_code)) return false;
    if (!newly_synthesized_metal_fraction_total.empty() &&
        !monotonicNondecreasing(*this, newly_synthesized_metal_fraction_total)) return false;
    if (!event_count_per_initial_mass_code_total.empty() &&
        !monotonicNondecreasing(*this, event_count_per_initial_mass_code_total)) return false;

    for (std::size_t flat_index = 0; flat_index < expected; ++flat_index) {
      const double returned = return_fraction_total[flat_index];
      const double ejected_metals = totalEjectedValues(*this)[flat_index];
      if (returned > 1.0 + k_validation_tolerance ||
          ejected_metals > returned + k_validation_tolerance) return false;
      if (!newly_synthesized_metal_fraction_total.empty() &&
          newly_synthesized_metal_fraction_total[flat_index] >
              ejected_metals + k_validation_tolerance) return false;
    }

    for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
      if (return_fraction_channel[channel].size() != expected ||
          totalEjectedChannelValues(*this, channel).size() != expected ||
          energy_erg_per_initial_mass_channel[channel].size() != expected) return false;
      if (!newly_synthesized_metal_fraction_channel[channel].empty() &&
          newly_synthesized_metal_fraction_channel[channel].size() != expected) return false;
      if (!event_count_per_initial_mass_code_channel[channel].empty() &&
          event_count_per_initial_mass_code_channel[channel].size() != expected) return false;
      if (!monotonicNondecreasing(*this, return_fraction_channel[channel]) ||
          !monotonicNondecreasing(*this, totalEjectedChannelValues(*this, channel)) ||
          !monotonicNondecreasing(*this, energy_erg_per_initial_mass_channel[channel])) return false;
    }

    for (std::size_t flat_index = 0; flat_index < expected; ++flat_index) {
      double returned_sum = 0.0;
      double metals_sum = 0.0;
      double energy_sum = 0.0;
      for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
        returned_sum += return_fraction_channel[channel][flat_index];
        metals_sum += totalEjectedChannelValues(*this, channel)[flat_index];
        energy_sum += energy_erg_per_initial_mass_channel[channel][flat_index];
      }
      const auto close = [](double lhs, double rhs) {
        return std::abs(lhs - rhs) <= 1.0e-8 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
      };
      if (!close(returned_sum, return_fraction_total[flat_index]) ||
          !close(metals_sum, totalEjectedValues(*this)[flat_index]) ||
          !close(energy_sum, energy_erg_per_initial_mass_code[flat_index])) return false;
    }
    return true;
  } catch (...) {
    return false;
  }
}

void StellarEvolutionTable::requireConsistent() const {
  if (!isConsistent()) {
    throw std::runtime_error("StellarEvolutionTable: invalid or non-conservative cumulative table");
  }
}

StellarEvolutionCumulative StellarEvolutionTable::evaluate(
    double age_yr_value,
    double birth_metallicity_mass_fraction_value) const {
  requireConsistent();
  if (!std::isfinite(age_yr_value) || !std::isfinite(birth_metallicity_mass_fraction_value)) {
    throw std::invalid_argument("stellar-evolution interpolation inputs must be finite");
  }
  const Bracket age = ageBracket(age_yr, std::max(age_yr_value, 0.0));
  const Bracket metallicity = birth_metallicity_mass_fraction.empty()
      ? Bracket{0U, 0U, 0.0}
      : linearBracket(
            birth_metallicity_mass_fraction,
            std::clamp(birth_metallicity_mass_fraction_value, 0.0, 1.0));

  StellarEvolutionCumulative result;
  result.return_fraction_total = bilinear(*this, return_fraction_total, age, metallicity);
  result.total_ejected_metal_fraction_total =
      bilinear(*this, totalEjectedValues(*this), age, metallicity);
  result.newly_synthesized_metal_fraction_total = newly_synthesized_metal_fraction_total.empty()
      ? 0.0
      : bilinear(*this, newly_synthesized_metal_fraction_total, age, metallicity);
  result.event_count_per_initial_mass_code_total = event_count_per_initial_mass_code_total.empty()
      ? 0.0
      : bilinear(*this, event_count_per_initial_mass_code_total, age, metallicity);
  result.energy_erg_per_initial_mass_code =
      bilinear(*this, energy_erg_per_initial_mass_code, age, metallicity);
  result.metal_yield_fraction_total = result.total_ejected_metal_fraction_total;

  for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
    result.return_fraction_channel[channel] =
        bilinear(*this, return_fraction_channel[channel], age, metallicity);
    result.total_ejected_metal_fraction_channel[channel] =
        bilinear(*this, totalEjectedChannelValues(*this, channel), age, metallicity);
    result.newly_synthesized_metal_fraction_channel[channel] =
        newly_synthesized_metal_fraction_channel[channel].empty()
        ? 0.0
        : bilinear(*this, newly_synthesized_metal_fraction_channel[channel], age, metallicity);
    result.event_count_per_initial_mass_code_channel[channel] =
        event_count_per_initial_mass_code_channel[channel].empty()
        ? 0.0
        : bilinear(*this, event_count_per_initial_mass_code_channel[channel], age, metallicity);
    result.energy_erg_per_initial_mass_channel[channel] =
        bilinear(*this, energy_erg_per_initial_mass_channel[channel], age, metallicity);
    result.metal_yield_fraction_channel[channel] =
        result.total_ejected_metal_fraction_channel[channel];
  }
  return result;
}

StellarEvolutionCumulative StellarEvolutionTable::evaluateAtAgeYears(
    double age_yr_value) const {
  const double reference_metallicity = birth_metallicity_mass_fraction.empty()
      ? 0.0
      : birth_metallicity_mass_fraction.front();
  return evaluate(age_yr_value, reference_metallicity);
}

StellarEvolutionIntervalBudget StellarEvolutionTable::integrateInterval(
    double age_begin_yr,
    double age_end_yr,
    double initial_mass_code,
    double birth_metallicity_mass_fraction_value) const {
  if (!std::isfinite(initial_mass_code) || initial_mass_code < 0.0) {
    throw std::invalid_argument("initial stellar mass must be finite and nonnegative");
  }
  const double begin = std::max(age_begin_yr, 0.0);
  const double end = std::max(age_end_yr, begin);
  const StellarEvolutionCumulative begin_state =
      evaluate(begin, birth_metallicity_mass_fraction_value);
  const StellarEvolutionCumulative end_state =
      evaluate(end, birth_metallicity_mass_fraction_value);

  StellarEvolutionIntervalBudget interval;
  interval.returned_mass_code = nonnegativeDifference(
      end_state.return_fraction_total, begin_state.return_fraction_total) * initial_mass_code;
  interval.returned_metals_code = nonnegativeDifference(
      end_state.total_ejected_metal_fraction_total,
      begin_state.total_ejected_metal_fraction_total) * initial_mass_code;
  interval.newly_synthesized_metals_code = nonnegativeDifference(
      end_state.newly_synthesized_metal_fraction_total,
      begin_state.newly_synthesized_metal_fraction_total) * initial_mass_code;
  interval.event_count = nonnegativeDifference(
      end_state.event_count_per_initial_mass_code_total,
      begin_state.event_count_per_initial_mass_code_total) * initial_mass_code;
  interval.feedback_energy_erg = nonnegativeDifference(
      end_state.energy_erg_per_initial_mass_code,
      begin_state.energy_erg_per_initial_mass_code) * initial_mass_code;

  if (interval.returned_metals_code > interval.returned_mass_code + k_validation_tolerance) {
    throw std::runtime_error("stellar-evolution interval ejects more metals than total mass");
  }
  for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
    interval.returned_mass_channel_code[channel] = nonnegativeDifference(
        end_state.return_fraction_channel[channel],
        begin_state.return_fraction_channel[channel]) * initial_mass_code;
    interval.returned_metals_channel_code[channel] = nonnegativeDifference(
        end_state.total_ejected_metal_fraction_channel[channel],
        begin_state.total_ejected_metal_fraction_channel[channel]) * initial_mass_code;
    interval.newly_synthesized_metals_channel_code[channel] = nonnegativeDifference(
        end_state.newly_synthesized_metal_fraction_channel[channel],
        begin_state.newly_synthesized_metal_fraction_channel[channel]) * initial_mass_code;
    interval.event_count_channel[channel] = nonnegativeDifference(
        end_state.event_count_per_initial_mass_code_channel[channel],
        begin_state.event_count_per_initial_mass_code_channel[channel]) * initial_mass_code;
    interval.feedback_energy_channel_erg[channel] = nonnegativeDifference(
        end_state.energy_erg_per_initial_mass_channel[channel],
        begin_state.energy_erg_per_initial_mass_channel[channel]) * initial_mass_code;
  }
  return interval;
}

StellarEvolutionIntervalBudget StellarEvolutionTable::integrateInterval(
    double age_begin_yr,
    double age_end_yr,
    double initial_mass_code) const {
  const double reference_metallicity = birth_metallicity_mass_fraction.empty()
      ? 0.0
      : birth_metallicity_mass_fraction.front();
  return integrateInterval(age_begin_yr, age_end_yr, initial_mass_code, reference_metallicity);
}

StellarEvolutionTable StellarEvolutionTable::loadFromTextFile(
    const std::string& path,
    const std::string& source_tag) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("StellarEvolutionTable: failed to open '" + path + "'");
  }

  struct Row {
    double metallicity = 0.0;
    double age = 0.0;
    double returned = 0.0;
    double total_metals = 0.0;
    double new_metals = 0.0;
    double events = 0.0;
    double energy = 0.0;
    std::array<double, 3> returned_channel{};
    std::array<double, 3> total_metals_channel{};
    std::array<double, 3> new_metals_channel{};
    std::array<double, 3> events_channel{};
    std::array<double, 3> energy_channel{};
  };

  StellarEvolutionTable table;
  table.source_path = path;
  std::vector<Row> rows;
  std::string line;
  std::size_t line_number = 0U;
  while (std::getline(input, line)) {
    ++line_number;
    const std::string normalized = trim(line);
    if (normalized.empty()) continue;
    if (normalized.front() == '#') {
      parseMetadataLine(table, trim(normalized.substr(1U)));
      continue;
    }
    Row row;
    std::istringstream stream(normalized);
    stream >> row.metallicity >> row.age >> row.returned >> row.total_metals >>
        row.new_metals >> row.events >> row.energy;
    for (double& value : row.returned_channel) stream >> value;
    for (double& value : row.total_metals_channel) stream >> value;
    for (double& value : row.new_metals_channel) stream >> value;
    for (double& value : row.events_channel) stream >> value;
    for (double& value : row.energy_channel) stream >> value;
    if (!stream) {
      throw std::runtime_error(
          "StellarEvolutionTable: malformed v2 row at " + path + ":" +
          std::to_string(line_number) + " (" + source_tag + ")");
    }
    std::string trailing_token;
    if (stream >> trailing_token) {
      throw std::runtime_error(
          "StellarEvolutionTable: unexpected trailing token at " + path + ":" +
          std::to_string(line_number) + " (" + source_tag + ")");
    }
    rows.push_back(row);
  }
  if (rows.empty()) {
    throw std::runtime_error("StellarEvolutionTable: no data rows in '" + path + "'");
  }

  std::sort(rows.begin(), rows.end(), [](const Row& lhs, const Row& rhs) {
    return std::tie(lhs.metallicity, lhs.age) < std::tie(rhs.metallicity, rhs.age);
  });
  for (const Row& row : rows) {
    if (table.birth_metallicity_mass_fraction.empty() ||
        row.metallicity != table.birth_metallicity_mass_fraction.back()) {
      table.birth_metallicity_mass_fraction.push_back(row.metallicity);
    }
  }
  for (const Row& row : rows) {
    if (row.metallicity == table.birth_metallicity_mass_fraction.front()) {
      table.age_yr.push_back(row.age);
    }
  }
  const std::size_t expected = table.birth_metallicity_mass_fraction.size() * table.age_yr.size();
  if (rows.size() != expected) {
    throw std::runtime_error("StellarEvolutionTable: table is not a complete tensor-product grid");
  }
  for (std::size_t index = 0; index < rows.size(); ++index) {
    const Row& row = rows[index];
    const std::size_t age_index = index % table.age_yr.size();
    if (row.age != table.age_yr[age_index]) {
      throw std::runtime_error("StellarEvolutionTable: age grid differs between metallicity planes");
    }
    table.return_fraction_total.push_back(row.returned);
    table.total_ejected_metal_fraction_total.push_back(row.total_metals);
    table.newly_synthesized_metal_fraction_total.push_back(row.new_metals);
    table.event_count_per_initial_mass_code_total.push_back(row.events);
    table.energy_erg_per_initial_mass_code.push_back(row.energy);
    for (std::size_t channel = 0; channel < 3U; ++channel) {
      table.return_fraction_channel[channel].push_back(row.returned_channel[channel]);
      table.total_ejected_metal_fraction_channel[channel].push_back(
          row.total_metals_channel[channel]);
      table.newly_synthesized_metal_fraction_channel[channel].push_back(
          row.new_metals_channel[channel]);
      table.event_count_per_initial_mass_code_channel[channel].push_back(
          row.events_channel[channel]);
      table.energy_erg_per_initial_mass_channel[channel].push_back(
          row.energy_channel[channel]);
    }
  }
  table.requireConsistent();
  return table;
}

StellarEvolutionTable StellarEvolutionTable::makeBuiltinReference() {
  // Deliberately non-calibrated and zero-yield. This safe compatibility table
  // prevents invented astrophysical constants from entering production runs.
  StellarEvolutionTable table;
  table.table_id = "builtin_zero_yield";
  table.table_version = "v2_safe_compatibility";
  table.source_path = "builtin";
  table.redistribution_license = "CHUI source license";
  table.imf = "none";
  table.stellar_mass_range = "not_applicable";
  table.solar_abundance_reference = "not_applicable";
  table.production_calibrated = false;
  table.age_yr = {0.0, 1.4e10};
  table.birth_metallicity_mass_fraction = {0.0, 0.04};
  const std::size_t count = table.age_yr.size() * table.birth_metallicity_mass_fraction.size();
  table.return_fraction_total.assign(count, 0.0);
  table.total_ejected_metal_fraction_total.assign(count, 0.0);
  table.newly_synthesized_metal_fraction_total.assign(count, 0.0);
  table.event_count_per_initial_mass_code_total.assign(count, 0.0);
  table.energy_erg_per_initial_mass_code.assign(count, 0.0);
  for (std::size_t channel = 0; channel < 3U; ++channel) {
    table.return_fraction_channel[channel].assign(count, 0.0);
    table.total_ejected_metal_fraction_channel[channel].assign(count, 0.0);
    table.newly_synthesized_metal_fraction_channel[channel].assign(count, 0.0);
    table.event_count_per_initial_mass_code_channel[channel].assign(count, 0.0);
    table.energy_erg_per_initial_mass_channel[channel].assign(count, 0.0);
  }
  table.requireConsistent();
  return table;
}

StellarEvolutionBookkeeper::StellarEvolutionBookkeeper(
    StellarEvolutionConfig config,
    StellarEvolutionTable table)
    : m_config(std::move(config)), m_table(std::move(table)) {
  m_table.requireConsistent();
  if (!(m_config.hubble_time_years > 0.0)) {
    throw std::invalid_argument("stellar_evolution_hubble_time_years must be positive");
  }
  if (m_config.require_production_calibrated_table && !m_table.production_calibrated) {
    throw std::invalid_argument(
        "production stellar evolution requires a calibrated table with provenance");
  }
}

const StellarEvolutionConfig& StellarEvolutionBookkeeper::config() const noexcept {
  return m_config;
}

const StellarEvolutionTable& StellarEvolutionBookkeeper::table() const noexcept {
  return m_table;
}

double StellarEvolutionBookkeeper::evaluateStarAgeYears(
    double formation_scale_factor,
    double current_scale_factor,
    const core::LambdaCdmBackground& background) const {
  if (current_scale_factor < formation_scale_factor) {
    throw std::invalid_argument("current scale factor precedes stellar formation");
  }
  return background.cosmicTimeIntervalSi(formation_scale_factor, current_scale_factor) /
      k_seconds_per_year;
}

double StellarEvolutionBookkeeper::evaluateStarAgeYears(
    double formation_scale_factor,
    double current_scale_factor) const {
  // Preserve source compatibility without retaining the old logarithmic proxy.
  const core::LambdaCdmBackground reference_background;
  return evaluateStarAgeYears(formation_scale_factor, current_scale_factor, reference_background);
}

StellarEvolutionRuntimeView StellarEvolutionBookkeeper::makeRuntimeView(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_star_indices) const {
  return StellarEvolutionRuntimeView{
      .active_star_indices = active_star_indices,
      .particle_index = state.star_particles.particle_index,
      .birth_mass_code = state.star_particles.birth_mass_code,
      .formation_scale_factor = state.star_particles.formation_scale_factor,
      .birth_metallicity_mass_fraction = state.star_particles.metallicity_mass_fraction,
      .stellar_age_years_last = state.star_particles.stellar_age_years_last,
      .stellar_returned_mass_cumulative_code =
          state.star_particles.stellar_returned_mass_cumulative_code,
      .stellar_returned_metals_cumulative_code =
          state.star_particles.stellar_returned_metals_cumulative_code,
      .stellar_newly_synthesized_metals_cumulative_code =
          state.star_particles.stellar_newly_synthesized_metals_cumulative_code,
      .stellar_feedback_energy_cumulative_erg =
          state.star_particles.stellar_feedback_energy_cumulative_erg,
      .returned_mass_channel_cumulative_code = {
          state.star_particles.stellar_returned_mass_channel_cumulative_code[0],
          state.star_particles.stellar_returned_mass_channel_cumulative_code[1],
          state.star_particles.stellar_returned_mass_channel_cumulative_code[2]},
      .returned_metals_channel_cumulative_code = {
          state.star_particles.stellar_returned_metals_channel_cumulative_code[0],
          state.star_particles.stellar_returned_metals_channel_cumulative_code[1],
          state.star_particles.stellar_returned_metals_channel_cumulative_code[2]},
      .feedback_energy_channel_cumulative_erg = {
          state.star_particles.stellar_feedback_energy_channel_cumulative_erg[0],
          state.star_particles.stellar_feedback_energy_channel_cumulative_erg[1],
          state.star_particles.stellar_feedback_energy_channel_cumulative_erg[2]},
      .particle_mass_code = state.particles.mass_code,
      .birth_tick = state.star_particles.birth_tick,
  };
}

StellarEvolutionStepReport StellarEvolutionBookkeeper::evaluateElapsedYears(
    const core::SimulationState& state,
    std::span<const std::uint32_t> active_star_indices,
    double elapsed_years) const {
  // The evaluation view is logically read-only; spans are mutable only because
  // the same narrow ABI is shared with commitBudgetsFromView.
  core::SimulationState& mutable_state = const_cast<core::SimulationState&>(state);
  return evaluateElapsedYearsFromView(
      makeRuntimeView(mutable_state, active_star_indices), elapsed_years);
}

StellarEvolutionStepReport StellarEvolutionBookkeeper::evaluateElapsedYearsFromView(
    StellarEvolutionRuntimeView view,
    double elapsed_years) const {
  StellarEvolutionStepReport report;
  if (!m_config.enabled || elapsed_years <= 0.0 || view.particle_index.empty()) {
    return report;
  }
  if (!std::isfinite(elapsed_years)) {
    throw std::invalid_argument("stellar elapsed time must be finite");
  }

  for (const std::uint32_t star_index : view.active_star_indices) {
    ++report.counters.scanned_stars;
    if (star_index >= view.particle_index.size() || star_index >= view.birth_mass_code.size() ||
        star_index >= view.birth_metallicity_mass_fraction.size() ||
        star_index >= view.stellar_age_years_last.size()) {
      continue;
    }
    const std::uint32_t particle_index = view.particle_index[star_index];
    if (particle_index >= view.particle_mass_code.size()) continue;
    const double birth_mass = view.birth_mass_code[star_index];
    const double current_mass = view.particle_mass_code[particle_index];
    if (!(birth_mass > k_mass_floor) || !(current_mass >= 0.0)) continue;

    const double age_begin = std::max(view.stellar_age_years_last[star_index], 0.0);
    const double star_elapsed = stellarElapsedYearsInInterval(
        star_index < view.birth_tick.size() ? view.birth_tick[star_index] : 0U, view.newborn_birth_tick,
        elapsed_years);
    if (!(star_elapsed > 0.0)) continue;  // F5
    const double age_end = age_begin + star_elapsed;
    StellarEvolutionIntervalBudget interval = m_table.integrateInterval(
        age_begin,
        age_end,
        birth_mass,
        view.birth_metallicity_mass_fraction[star_index]);
    const double returned_mass = std::min(interval.returned_mass_code, current_mass);
    if (interval.returned_mass_code > 0.0 && returned_mass < interval.returned_mass_code) {
      const double scale = returned_mass / interval.returned_mass_code;
      interval.returned_mass_code = returned_mass;
      interval.returned_metals_code *= scale;
      interval.newly_synthesized_metals_code *= scale;
      interval.feedback_energy_erg *= scale;
      interval.event_count *= scale;
      for (std::size_t channel = 0; channel < 3U; ++channel) {
        interval.returned_mass_channel_code[channel] *= scale;
        interval.returned_metals_channel_code[channel] *= scale;
        interval.newly_synthesized_metals_channel_code[channel] *= scale;
        interval.feedback_energy_channel_erg[channel] *= scale;
        interval.event_count_channel[channel] *= scale;
      }
    }
    if (interval.returned_metals_code > interval.returned_mass_code + k_validation_tolerance) {
      throw std::runtime_error("stellar interval violates metal <= returned mass");
    }

    StellarEvolutionStarBudget budget{
        .star_index = star_index,
        .particle_index = particle_index,
        .star_age_begin_years = age_begin,
        .star_age_end_years = age_end,
        .mass_old_code = current_mass,
        .mass_new_code = std::max(current_mass - interval.returned_mass_code, 0.0),
        .interval = interval,
    };
    ++report.counters.evolved_stars;
    report.counters.returned_mass_code += interval.returned_mass_code;
    report.counters.returned_metals_code += interval.returned_metals_code;
    report.counters.newly_synthesized_metals_code += interval.newly_synthesized_metals_code;
    report.counters.event_count += interval.event_count;
    report.counters.feedback_energy_erg += interval.feedback_energy_erg;
    report.budgets.push_back(std::move(budget));
  }
  return report;
}

double stellarElapsedYearsInInterval(std::uint64_t birth_tick, std::uint64_t newborn_birth_tick,
                                     double interval_years) noexcept {
  if (newborn_birth_tick != std::numeric_limits<std::uint64_t>::max() && birth_tick == newborn_birth_tick) {
    return 0.0;  // born at t_end of this interval (F5)
  }
  return std::max(interval_years, 0.0);
}

StellarEvolutionStepCounters StellarEvolutionBookkeeper::evaluateElapsedYearsGoverned(
    const core::SimulationState& state,
    std::span<const std::uint32_t> active_star_indices,
    double elapsed_years,
    core::MemoryGovernor* governor,
    StellarEvolutionBatchWorkspace& workspace,
    std::uint64_t newborn_birth_tick) const {
  // Physical owner: SourceRuntime's workspace.budgets owns the bytes.
  // Capacity admitted via governor BEFORE growth; retained capacity()
  // reconciled so the buffer is reused next batch. Growth occurs only after
  // admission, tight headroom rejects before any ledger mutation, and retry
  // with sufficient headroom reproduces identical budgets/counters because
  // evaluation is a pure function of (state, active set, elapsed_years).
  workspace.budgets.clear();
  StellarEvolutionStepCounters counters;
  if (!m_config.enabled || elapsed_years <= 0.0) {
    return counters;
  }
  if (!std::isfinite(elapsed_years)) {
    throw std::invalid_argument("stellar elapsed time must be finite");
  }
  if (active_star_indices.empty()) {
    return counters;
  }
  const std::size_t required_capacity = active_star_indices.size();
  const std::uint64_t required_bytes = stellarEvolutionBatchStagingBytes(required_capacity);
  const std::uint64_t retained_bytes = stellarEvolutionBatchStagingBytes(workspace.budgets.capacity());
  if (required_capacity > workspace.budgets.capacity()) {
    // Admit growth BEFORE vector growth. Old/new coexistence: the old retained
    // buffer stays live while the replacement is admitted.
    core::MemoryReservation growth;
    if (governor != nullptr && required_bytes > retained_bytes) {
      // Reserve the full replacement (not just the delta) so old/new
      // coexistence is charged during vector reallocation.
      growth = governor->reserve(core::MemoryClass::kPhaseResident, required_bytes,
                                 "sources.stellar_evolution.governed_budgets");
      growth.commit();
    } else if (governor != nullptr && required_bytes > 0U && workspace.admitted_bytes == 0U) {
      growth = governor->reserve(core::MemoryClass::kPhaseResident, required_bytes,
                                 "sources.stellar_evolution.governed_budgets");
      growth.commit();
    }
    workspace.budgets.reserve(required_capacity);
    if (governor != nullptr) {
      growth.release();
      workspace.reservation.release();
      const std::uint64_t actual_bytes =
          stellarEvolutionBatchStagingBytes(workspace.budgets.capacity());
      if (actual_bytes > 0U) {
        workspace.reservation = governor->reserve(core::MemoryClass::kPhaseResident,
                                                  actual_bytes,
                                                  "sources.stellar_evolution.governed_budgets");
        workspace.reservation.commit();
      }
      workspace.admitted_bytes = actual_bytes;
    }
  }
  // Fill within admitted capacity: no further allocation can occur because
  // size <= capacity. Deterministic input order is preserved.
  core::SimulationState& mutable_state = const_cast<core::SimulationState&>(state);
  StellarEvolutionRuntimeView view = makeRuntimeView(mutable_state, active_star_indices);
  if (view.particle_index.empty()) {
    return counters;
  }
  for (const std::uint32_t star_index : view.active_star_indices) {
    ++counters.scanned_stars;
    if (star_index >= view.particle_index.size() || star_index >= view.birth_mass_code.size() ||
        star_index >= view.birth_metallicity_mass_fraction.size() ||
        star_index >= view.stellar_age_years_last.size()) {
      continue;
    }
    const std::uint32_t particle_index = view.particle_index[star_index];
    if (particle_index >= view.particle_mass_code.size()) continue;
    const double birth_mass = view.birth_mass_code[star_index];
    const double current_mass = view.particle_mass_code[particle_index];
    if (!(birth_mass > k_mass_floor) || !(current_mass >= 0.0)) continue;
    const double age_begin = std::max(view.stellar_age_years_last[star_index], 0.0);
    const double star_elapsed = stellarElapsedYearsInInterval(
        star_index < view.birth_tick.size() ? view.birth_tick[star_index] : 0U, newborn_birth_tick, elapsed_years);
    if (!(star_elapsed > 0.0)) {
      ++counters.scanned_newborn_stars;
      continue;  // F5: zero age in the interval that created it
    }
    const double age_end = age_begin + star_elapsed;
    StellarEvolutionIntervalBudget interval = m_table.integrateInterval(
        age_begin, age_end, birth_mass,
        view.birth_metallicity_mass_fraction[star_index]);
    const double returned_mass = std::min(interval.returned_mass_code, current_mass);
    if (interval.returned_mass_code > 0.0 && returned_mass < interval.returned_mass_code) {
      const double scale = returned_mass / interval.returned_mass_code;
      interval.returned_mass_code = returned_mass;
      interval.returned_metals_code *= scale;
      interval.newly_synthesized_metals_code *= scale;
      interval.feedback_energy_erg *= scale;
      interval.event_count *= scale;
      for (std::size_t channel = 0; channel < 3U; ++channel) {
        interval.returned_mass_channel_code[channel] *= scale;
        interval.returned_metals_channel_code[channel] *= scale;
        interval.newly_synthesized_metals_channel_code[channel] *= scale;
        interval.feedback_energy_channel_erg[channel] *= scale;
        interval.event_count_channel[channel] *= scale;
      }
    }
    if (interval.returned_metals_code > interval.returned_mass_code + k_validation_tolerance) {
      throw std::runtime_error("stellar interval violates metal <= returned mass");
    }
    StellarEvolutionStarBudget budget{
        .star_index = star_index,
        .particle_index = particle_index,
        .star_age_begin_years = age_begin,
        .star_age_end_years = age_end,
        .mass_old_code = current_mass,
        .mass_new_code = std::max(current_mass - interval.returned_mass_code, 0.0),
        .interval = interval,
    };
    ++counters.evolved_stars;
    counters.returned_mass_code += interval.returned_mass_code;
    counters.returned_metals_code += interval.returned_metals_code;
    counters.newly_synthesized_metals_code += interval.newly_synthesized_metals_code;
    counters.event_count += interval.event_count;
    counters.feedback_energy_erg += interval.feedback_energy_erg;
    workspace.budgets.push_back(std::move(budget));
  }
  return counters;
}

void StellarEvolutionBookkeeper::commitBudgets(
    core::SimulationState& state,
    const StellarEvolutionStepReport& report) const {
  commitBudgetsFromView(makeRuntimeView(state, {}), report);
  // View-only kernels cannot mutate SimulationState's generation; the state-
  // owning commit publishes this mutation before any subsequent mass cache use.
  if (!report.budgets.empty()) state.notifyParticleMassMutation();
  state.sidecars.upsert(buildMetadataSidecar(report));
}

void StellarEvolutionBookkeeper::commitBudgetsFromView(
    StellarEvolutionRuntimeView view,
    const StellarEvolutionStepReport& report) const {
  for (const StellarEvolutionStarBudget& budget : report.budgets) {
    const std::size_t star_index = budget.star_index;
    if (star_index >= view.particle_index.size() ||
        budget.particle_index >= view.particle_mass_code.size()) {
      throw std::runtime_error("stellar-evolution commit references a stale row");
    }
    view.particle_mass_code[budget.particle_index] = budget.mass_new_code;
    view.stellar_age_years_last[star_index] = budget.star_age_end_years;
    view.stellar_returned_mass_cumulative_code[star_index] +=
        budget.interval.returned_mass_code;
    view.stellar_returned_metals_cumulative_code[star_index] +=
        budget.interval.returned_metals_code;
    if (!view.stellar_newly_synthesized_metals_cumulative_code.empty()) {
      view.stellar_newly_synthesized_metals_cumulative_code[star_index] +=
          budget.interval.newly_synthesized_metals_code;
    }
    view.stellar_feedback_energy_cumulative_erg[star_index] +=
        budget.interval.feedback_energy_erg;
    for (std::size_t channel = 0; channel < 3U; ++channel) {
      view.returned_mass_channel_cumulative_code[channel][star_index] +=
          budget.interval.returned_mass_channel_code[channel];
      view.returned_metals_channel_cumulative_code[channel][star_index] +=
          budget.interval.returned_metals_channel_code[channel];
      view.feedback_energy_channel_cumulative_erg[channel][star_index] +=
          budget.interval.feedback_energy_channel_erg[channel];
    }
  }
}

StellarEvolutionStepReport StellarEvolutionBookkeeper::applyElapsedYears(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_star_indices,
    double elapsed_years) const {
  StellarEvolutionStepReport report =
      evaluateElapsedYears(state, active_star_indices, elapsed_years);
  commitBudgets(state, report);
  return report;
}

StellarEvolutionStepReport StellarEvolutionBookkeeper::apply(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_star_indices,
    double current_scale_factor,
    double dt_code) const {
  (void)current_scale_factor;
  return applyElapsedYears(
      state, active_star_indices, std::max(dt_code * m_config.hubble_time_years, 0.0));
}

StellarEvolutionStepReport StellarEvolutionBookkeeper::applyFromView(
    StellarEvolutionRuntimeView view,
    double current_scale_factor,
    double dt_code) const {
  (void)current_scale_factor;
  StellarEvolutionStepReport report = evaluateElapsedYearsFromView(
      view, std::max(dt_code * m_config.hubble_time_years, 0.0));
  commitBudgetsFromView(view, report);
  return report;
}

core::ModuleSidecarBlock StellarEvolutionBookkeeper::buildMetadataSidecar(
    const StellarEvolutionStepReport& report) const {
  std::ostringstream stream;
  stream << "module=stellar_evolution\n";
  stream << "schema_version=" << m_config.metadata_schema_version << "\n";
  stream << "table_id=" << m_table.table_id << "\n";
  stream << "table_version=" << m_table.table_version << "\n";
  stream << "table_source_path=" << m_table.source_path << "\n";
  stream << "table_sha256=" << m_table.sha256 << "\n";
  stream << "production_calibrated=" << (m_table.production_calibrated ? "true" : "false") << "\n";
  stream << "imf=" << m_table.imf << "\n";
  stream << "solar_abundance_reference=" << m_table.solar_abundance_reference << "\n";
  stream << "scanned_stars=" << report.counters.scanned_stars << "\n";
  stream << "evolved_stars=" << report.counters.evolved_stars << "\n";
  stream << "returned_mass_code=" << report.counters.returned_mass_code << "\n";
  stream << "returned_metals_code=" << report.counters.returned_metals_code << "\n";
  stream << "newly_synthesized_metals_code="
         << report.counters.newly_synthesized_metals_code << "\n";
  stream << "event_count=" << report.counters.event_count << "\n";
  stream << "feedback_energy_erg=" << report.counters.feedback_energy_erg << "\n";

  const std::string text = stream.str();
  core::ModuleSidecarBlock block;
  block.module_name = "stellar_evolution";
  block.schema_version = m_config.metadata_schema_version;
  block.payload.resize(text.size());
  for (std::size_t index = 0; index < text.size(); ++index) {
    block.payload[index] = static_cast<std::byte>(text[index]);
  }
  return block;
}

StellarEvolutionConfig makeStellarEvolutionConfig(
    const core::PhysicsConfig& physics_config) {
  StellarEvolutionConfig config;
  config.enabled = physics_config.enable_stellar_evolution;
  config.table_path = physics_config.stellar_evolution_table_path;
  config.hubble_time_years = physics_config.stellar_evolution_hubble_time_years;
  config.require_production_calibrated_table =
      physics_config.stellar_evolution_require_production_table;
  return config;
}

StellarEvolutionTable loadStellarEvolutionTable(
    const core::PhysicsConfig& physics_config) {
  if (physics_config.stellar_evolution_table_path.empty()) {
    return StellarEvolutionTable::makeBuiltinReference();
  }
  return StellarEvolutionTable::loadFromTextFile(
      physics_config.stellar_evolution_table_path,
      "physics.stellar_evolution_table_path");
}

}  // namespace cosmosim::physics
