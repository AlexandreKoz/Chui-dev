#include "cosmosim/core/gas_composition_sidecar.hpp"

#include <cmath>
#include <stdexcept>
#include <string_view>

namespace cosmosim::core {

bool GasCompositionLayoutDescriptor::isConsistent() const noexcept {
  if (!enabled()) {
    return lane_count == 0U;
  }
  const auto within = [this](std::uint32_t begin, std::uint32_t count) {
    return count == 0U || (begin <= lane_count && count <= lane_count - begin);
  };
  if (!within(element_lane_begin, element_lane_count) ||
      !within(species_lane_begin, species_lane_count) ||
      !within(dust_lane_begin, dust_lane_count)) {
    return false;
  }
  if (sgs_energy_lane != k_no_lane && sgs_energy_lane >= lane_count) {
    return false;
  }
  return element_lane_count != 0U && layout_signature != 0U;
}

bool GasCompositionLayoutDescriptor::operator==(
    const GasCompositionLayoutDescriptor& rhs) const noexcept {
  return layout_version == rhs.layout_version && layout_signature == rhs.layout_signature &&
      layout_name == rhs.layout_name && lane_count == rhs.lane_count &&
      element_lane_begin == rhs.element_lane_begin && element_lane_count == rhs.element_lane_count &&
      species_lane_begin == rhs.species_lane_begin && species_lane_count == rhs.species_lane_count &&
      dust_lane_begin == rhs.dust_lane_begin && dust_lane_count == rhs.dust_lane_count &&
      sgs_energy_lane == rhs.sgs_energy_lane;
}

void GasCompositionSidecar::configure(
    const GasCompositionLayoutDescriptor& descriptor, std::size_t row_count) {
  if (!descriptor.isConsistent()) {
    throw std::invalid_argument(
        "GasCompositionSidecar::configure: inconsistent layout descriptor '" +
        descriptor.layout_name + "'");
  }
  clear();
  if (!descriptor.enabled()) {
    return;
  }
  m_descriptor = descriptor;
  m_lanes.resize(descriptor.lane_count);
  resizeRows(row_count);
}

void GasCompositionSidecar::clear() noexcept {
  m_descriptor = GasCompositionLayoutDescriptor{};
  m_lanes.clear();
  m_lanes.shrink_to_fit();
  gas_cell_id = {};
  chemistry_model_tag = {};
  ism_regime_tag = {};
  model_dwell_years = {};
  ism_regime_dwell_years = {};
  mean_molecular_weight = {};
}

void GasCompositionSidecar::resizeRows(std::size_t row_count) {
  if (!enabled()) {
    if (row_count != 0U && !m_lanes.empty()) {
      throw std::logic_error("GasCompositionSidecar: lanes present without an enabled layout");
    }
    return;
  }
  for (auto& values : m_lanes) {
    values.resize(row_count, 0.0);
  }
  gas_cell_id.resize(row_count, 0U);
  chemistry_model_tag.resize(row_count, 0U);
  ism_regime_tag.resize(row_count, 0U);
  model_dwell_years.resize(row_count, 0.0);
  ism_regime_dwell_years.resize(row_count, 0.0);
  mean_molecular_weight.resize(row_count, std::numeric_limits<double>::quiet_NaN());
}

std::span<double> GasCompositionSidecar::lane(std::size_t lane_index) {
  if (lane_index >= m_lanes.size()) {
    throw std::out_of_range("GasCompositionSidecar::lane index out of range");
  }
  return m_lanes[lane_index];
}

std::span<const double> GasCompositionSidecar::lane(std::size_t lane_index) const {
  if (lane_index >= m_lanes.size()) {
    throw std::out_of_range("GasCompositionSidecar::lane index out of range");
  }
  return m_lanes[lane_index];
}

double GasCompositionSidecar::value(std::size_t lane_index, std::size_t row) const {
  return lane(lane_index)[row];
}

void GasCompositionSidecar::setValue(std::size_t lane_index, std::size_t row, double value) {
  lane(lane_index)[row] = value;
}

void GasCompositionSidecar::addValue(std::size_t lane_index, std::size_t row, double delta) {
  lane(lane_index)[row] += delta;
}

double GasCompositionSidecar::elementMassSum(std::size_t row) const {
  double sum = 0.0;
  for (std::uint32_t l = 0; l < m_descriptor.element_lane_count; ++l) {
    sum += m_lanes[m_descriptor.element_lane_begin + l][row];
  }
  return sum;
}

void GasCompositionSidecar::scaleRow(std::size_t row, double factor) {
  if (!(factor >= 0.0) || factor > 1.0 || !std::isfinite(factor)) {
    throw std::invalid_argument("GasCompositionSidecar::scaleRow factor must be in [0, 1]");
  }
  for (auto& values : m_lanes) {
    values[row] *= factor;
  }
}

std::vector<double> GasCompositionSidecar::packRow(std::size_t row) const {
  std::vector<double> payload;
  if (!enabled()) {
    return payload;
  }
  if (row >= rowCount()) {
    throw std::out_of_range("GasCompositionSidecar::packRow row out of range");
  }
  payload.reserve(packedRowWidth());
  for (const auto& values : m_lanes) {
    payload.push_back(values[row]);
  }
  payload.push_back(static_cast<double>(chemistry_model_tag[row]));
  payload.push_back(static_cast<double>(ism_regime_tag[row]));
  payload.push_back(model_dwell_years[row]);
  payload.push_back(ism_regime_dwell_years[row]);
  return payload;
}

void GasCompositionSidecar::unpackRow(
    std::size_t row, std::span<const double> payload, std::uint64_t gas_cell_id_value) {
  if (!enabled()) {
    if (!payload.empty()) {
      throw std::runtime_error(
          "GasCompositionSidecar: received a composition payload while composition is disabled");
    }
    return;
  }
  if (payload.size() != packedRowWidth()) {
    throw std::runtime_error(
        "GasCompositionSidecar: migration payload width " + std::to_string(payload.size()) +
        " does not match layout '" + m_descriptor.layout_name + "' width " +
        std::to_string(packedRowWidth()) +
        " (a gas cell would lose or invent composition lanes)");
  }
  if (row >= rowCount()) {
    throw std::out_of_range("GasCompositionSidecar::unpackRow row out of range");
  }
  for (std::size_t l = 0; l < m_lanes.size(); ++l) {
    m_lanes[l][row] = payload[l];
  }
  const std::size_t meta = m_lanes.size();
  chemistry_model_tag[row] = static_cast<std::uint8_t>(payload[meta]);
  ism_regime_tag[row] = static_cast<std::uint8_t>(payload[meta + 1U]);
  model_dwell_years[row] = payload[meta + 2U];
  ism_regime_dwell_years[row] = payload[meta + 3U];
  mean_molecular_weight[row] = std::numeric_limits<double>::quiet_NaN();
  gas_cell_id[row] = gas_cell_id_value;
}

GasCompositionSidecar GasCompositionSidecar::emptyLike(std::size_t row_count) const {
  GasCompositionSidecar result;
  if (enabled()) {
    result.configure(m_descriptor, row_count);
  }
  return result;
}

void GasCompositionSidecar::requireAlignedWith(
    std::span<const std::uint64_t> expected_gas_cell_ids, const char* caller) const {
  if (!enabled()) {
    return;
  }
  if (expected_gas_cell_ids.size() != rowCount()) {
    throw std::runtime_error(std::string(caller) +
        ": composition sidecar row count " + std::to_string(rowCount()) +
        " differs from gas row count " + std::to_string(expected_gas_cell_ids.size()));
  }
  for (std::size_t row = 0; row < rowCount(); ++row) {
    if (gas_cell_id[row] != expected_gas_cell_ids[row]) {
      throw std::runtime_error(std::string(caller) +
          ": composition row " + std::to_string(row) + " carries gas_cell_id " +
          std::to_string(gas_cell_id[row]) + " but the gas row is " +
          std::to_string(expected_gas_cell_ids[row]) +
          " (row permutation would apply another cell's chemistry)");
    }
  }
}

bool GasCompositionSidecar::isConsistent(std::size_t expected_rows) const noexcept {
  if (!enabled()) {
    return m_lanes.empty() && gas_cell_id.empty();
  }
  if (!m_descriptor.isConsistent() || m_lanes.size() != m_descriptor.lane_count ||
      gas_cell_id.size() != expected_rows || chemistry_model_tag.size() != expected_rows ||
      ism_regime_tag.size() != expected_rows || model_dwell_years.size() != expected_rows ||
      ism_regime_dwell_years.size() != expected_rows ||
      mean_molecular_weight.size() != expected_rows) {
    return false;
  }
  for (const auto& values : m_lanes) {
    if (values.size() != expected_rows) return false;
    for (const double v : values) {
      if (!std::isfinite(v) || v < 0.0) return false;
    }
  }
  return true;
}

std::uint64_t GasCompositionSidecar::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (const auto& values : m_lanes) bytes += values.capacity() * sizeof(double);
  bytes += gas_cell_id.capacity() * sizeof(std::uint64_t);
  bytes += chemistry_model_tag.capacity() + ism_regime_tag.capacity();
  bytes += (model_dwell_years.capacity() + ism_regime_dwell_years.capacity() + mean_molecular_weight.capacity()) *
      sizeof(double);
  return bytes;
}

std::uint64_t GasCompositionSidecar::currentSizeBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (const auto& values : m_lanes) bytes += values.size() * sizeof(double);
  bytes += gas_cell_id.size() * sizeof(std::uint64_t);
  bytes += chemistry_model_tag.size() + ism_regime_tag.size();
  bytes += (model_dwell_years.size() + ism_regime_dwell_years.size() + mean_molecular_weight.size()) * sizeof(double);
  return bytes;
}

std::uint64_t fnv1a64(std::string_view text) noexcept {
  std::uint64_t hash = 14695981039346656037ULL;
  for (const char c : text) {
    hash ^= static_cast<std::uint8_t>(c);
    hash *= 1099511628211ULL;
  }
  return hash;
}

}  // namespace cosmosim::core
