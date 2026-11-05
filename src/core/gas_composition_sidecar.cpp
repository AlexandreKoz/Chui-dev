#include "cosmosim/core/gas_composition_sidecar.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string_view>

#include "cosmosim/core/retained_capacity_transaction.hpp"

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
    const GasCompositionLayoutDescriptor& descriptor,
    std::size_t row_count,
    std::uint32_t sparse_lane_begin,
    std::span<const std::uint32_t> sparse_default_source_lanes) {
  if (!descriptor.isConsistent()) {
    throw std::invalid_argument(
        "GasCompositionSidecar::configure: inconsistent layout descriptor '" +
        descriptor.layout_name + "'");
  }
  clear();
  if (!descriptor.enabled()) {
    return;
  }
  if (sparse_lane_begin == GasCompositionLayoutDescriptor::k_no_lane) {
    if (!sparse_default_source_lanes.empty()) {
      throw std::invalid_argument(
          "GasCompositionSidecar::configure: sparse defaults supplied without sparse lanes");
    }
  } else {
    if (sparse_lane_begin >= descriptor.lane_count || sparse_default_source_lanes.empty() ||
        sparse_default_source_lanes.size() > descriptor.lane_count - sparse_lane_begin) {
      throw std::invalid_argument(
          "GasCompositionSidecar::configure: invalid contiguous sparse lane range");
    }
    for (const std::uint32_t source : sparse_default_source_lanes) {
      if (source != GasCompositionLayoutDescriptor::k_no_lane &&
          (source >= sparse_lane_begin || source >= descriptor.lane_count)) {
        throw std::invalid_argument(
            "GasCompositionSidecar::configure: sparse implicit source must be an earlier dense lane");
      }
    }
  }
  m_descriptor = descriptor;
  m_lanes.resize(descriptor.lane_count);
  m_sparse_lane_begin = sparse_lane_begin;
  m_sparse_lane_count = static_cast<std::uint32_t>(sparse_default_source_lanes.size());
  m_sparse_default_source_lanes.assign(
      sparse_default_source_lanes.begin(), sparse_default_source_lanes.end());
  resizeRows(row_count);
}

void GasCompositionSidecar::clear() noexcept {
  m_descriptor = GasCompositionLayoutDescriptor{};
  m_lanes.clear();
  m_lanes.shrink_to_fit();
  m_sparse_lane_begin = GasCompositionLayoutDescriptor::k_no_lane;
  m_sparse_lane_count = 0U;
  m_sparse_default_source_lanes.clear();
  m_sparse_default_source_lanes.shrink_to_fit();
  m_sparse_gas_cell_id = {};
  m_sparse_values = {};
  gas_cell_id = {};
  chemistry_model_tag = {};
  ism_regime_tag = {};
  model_dwell_years = {};
  ism_regime_dwell_years = {};
  chemistry_last_temperature_k = {};
  chemistry_last_cooling_time_seconds = {};
  chemistry_last_error_estimate = {};
  chemistry_last_element_fraction_residual = {};
  chemistry_last_charge_residual = {};
  chemistry_cumulative_radiative_loss_erg = {};
  chemistry_cumulative_heating_gain_erg = {};
  mean_molecular_weight = {};
}

void GasCompositionSidecar::resizeRows(std::size_t row_count) {
  if (!enabled()) {
    if (row_count != 0U && !m_lanes.empty()) {
      throw std::logic_error("GasCompositionSidecar: lanes present without an enabled layout");
    }
    return;
  }
  if (row_count < gas_cell_id.size() && sparseStorageEnabled()) {
    std::size_t entry = 0U;
    while (entry < m_sparse_gas_cell_id.size()) {
      const std::uint64_t id = m_sparse_gas_cell_id[entry];
      const bool retained = std::find(gas_cell_id.begin(), gas_cell_id.begin() + row_count, id) !=
          gas_cell_id.begin() + row_count;
      if (retained) {
        ++entry;
        continue;
      }
      m_sparse_gas_cell_id.erase(m_sparse_gas_cell_id.begin() + entry);
      const auto first = m_sparse_values.begin() + entry * m_sparse_lane_count;
      m_sparse_values.erase(first, first + m_sparse_lane_count);
    }
  }
  for (std::size_t lane_index = 0; lane_index < m_lanes.size(); ++lane_index) {
    if (isSparseLane(lane_index)) {
      m_lanes[lane_index].clear();
    } else {
      m_lanes[lane_index].resize(row_count, 0.0);
    }
  }
  gas_cell_id.resize(row_count, 0U);
  chemistry_model_tag.resize(row_count, 0U);
  ism_regime_tag.resize(row_count, 0U);
  model_dwell_years.resize(row_count, 0.0);
  ism_regime_dwell_years.resize(row_count, 0.0);
  chemistry_last_temperature_k.resize(row_count, 0.0);
  chemistry_last_cooling_time_seconds.resize(row_count, 0.0);
  chemistry_last_error_estimate.resize(row_count, 0.0);
  chemistry_last_element_fraction_residual.resize(row_count, 0.0);
  chemistry_last_charge_residual.resize(row_count, 0.0);
  chemistry_cumulative_radiative_loss_erg.resize(row_count, 0.0);
  chemistry_cumulative_heating_gain_erg.resize(row_count, 0.0);
  mean_molecular_weight.resize(row_count, std::numeric_limits<double>::quiet_NaN());
}

std::span<double> GasCompositionSidecar::lane(std::size_t lane_index) {
  if (lane_index >= m_lanes.size()) {
    throw std::out_of_range("GasCompositionSidecar::lane index out of range");
  }
  if (isSparseLane(lane_index)) {
    throw std::logic_error(
        "GasCompositionSidecar::lane cannot expose a non-contiguous sparse lane; use value/setValue");
  }
  return m_lanes[lane_index];
}

std::span<const double> GasCompositionSidecar::lane(std::size_t lane_index) const {
  if (lane_index >= m_lanes.size()) {
    throw std::out_of_range("GasCompositionSidecar::lane index out of range");
  }
  if (isSparseLane(lane_index)) {
    throw std::logic_error(
        "GasCompositionSidecar::lane cannot expose a non-contiguous sparse lane; use value/setValue");
  }
  return m_lanes[lane_index];
}

bool GasCompositionSidecar::isSparseLane(std::size_t lane_index) const noexcept {
  return sparseStorageEnabled() && lane_index >= m_sparse_lane_begin &&
      lane_index < static_cast<std::size_t>(m_sparse_lane_begin) + m_sparse_lane_count;
}

bool GasCompositionSidecar::laneStoredDensely(std::size_t lane_index) const noexcept {
  return lane_index < m_lanes.size() && !isSparseLane(lane_index);
}

std::size_t GasCompositionSidecar::requireRow(std::size_t row) const {
  if (row >= rowCount()) {
    throw std::out_of_range("GasCompositionSidecar row out of range");
  }
  return row;
}

std::size_t GasCompositionSidecar::sparseEntryForId(std::uint64_t id) const noexcept {
  const auto it = std::lower_bound(m_sparse_gas_cell_id.begin(), m_sparse_gas_cell_id.end(), id);
  return it != m_sparse_gas_cell_id.end() && *it == id
      ? static_cast<std::size_t>(it - m_sparse_gas_cell_id.begin())
      : m_sparse_gas_cell_id.size();
}

bool GasCompositionSidecar::hasExplicitSparseRow(std::size_t row) const {
  requireRow(row);
  return sparseStorageEnabled() &&
      sparseEntryForId(gas_cell_id[row]) != m_sparse_gas_cell_id.size();
}

double GasCompositionSidecar::implicitSparseValue(
    std::size_t lane_index, std::size_t row) const {
  requireRow(row);
  if (!isSparseLane(lane_index)) {
    throw std::invalid_argument("GasCompositionSidecar::implicitSparseValue requires a sparse lane");
  }
  const std::uint32_t source =
      m_sparse_default_source_lanes[lane_index - m_sparse_lane_begin];
  return source == GasCompositionLayoutDescriptor::k_no_lane ? 0.0 : m_lanes[source][row];
}

std::uint32_t GasCompositionSidecar::sparseDefaultSourceLane(
    std::size_t lane_index) const {
  if (!isSparseLane(lane_index)) {
    throw std::invalid_argument(
        "GasCompositionSidecar::sparseDefaultSourceLane requires a sparse lane");
  }
  return m_sparse_default_source_lanes[lane_index - m_sparse_lane_begin];
}

double GasCompositionSidecar::value(std::size_t lane_index, std::size_t row) const {
  requireRow(row);
  if (lane_index >= m_lanes.size()) {
    throw std::out_of_range("GasCompositionSidecar::value lane out of range");
  }
  if (!isSparseLane(lane_index)) return m_lanes[lane_index][row];
  const std::size_t entry = sparseEntryForId(gas_cell_id[row]);
  return entry == m_sparse_gas_cell_id.size()
      ? implicitSparseValue(lane_index, row)
      : m_sparse_values[entry * m_sparse_lane_count + lane_index - m_sparse_lane_begin];
}

void GasCompositionSidecar::setValue(std::size_t lane_index, std::size_t row, double value) {
  requireRow(row);
  if (lane_index >= m_lanes.size()) {
    throw std::out_of_range("GasCompositionSidecar::setValue lane out of range");
  }
  if (!isSparseLane(lane_index)) {
    m_lanes[lane_index][row] = value;
    return;
  }
  std::size_t entry = sparseEntryForId(gas_cell_id[row]);
  if (entry == m_sparse_gas_cell_id.size()) {
    if (value == implicitSparseValue(lane_index, row)) return;
    activateSparseRow(row);
    entry = sparseEntryForId(gas_cell_id[row]);
  }
  m_sparse_values[entry * m_sparse_lane_count + lane_index - m_sparse_lane_begin] = value;
}

void GasCompositionSidecar::addValue(std::size_t lane_index, std::size_t row, double delta) {
  setValue(lane_index, row, value(lane_index, row) + delta);
}

void GasCompositionSidecar::activateSparseRow(std::size_t row) {
  requireRow(row);
  if (!sparseStorageEnabled() || hasExplicitSparseRow(row)) return;
  if (m_sparse_gas_cell_id.size() == m_sparse_gas_cell_id.capacity() ||
      m_sparse_values.capacity() < m_sparse_lane_count ||
      m_sparse_values.size() > m_sparse_values.capacity() - m_sparse_lane_count) {
    throw std::logic_error(
        "GasCompositionSidecar sparse activation exceeds pre-admitted capacity");
  }
  const std::uint64_t id = gas_cell_id[row];
  if (id == 0U) {
    throw std::logic_error(
        "GasCompositionSidecar cannot activate sparse chemistry without a stable gas_cell_id");
  }
  const auto id_it = std::lower_bound(m_sparse_gas_cell_id.begin(), m_sparse_gas_cell_id.end(), id);
  const std::size_t entry = static_cast<std::size_t>(id_it - m_sparse_gas_cell_id.begin());
  m_sparse_gas_cell_id.insert(id_it, id);
  m_sparse_values.insert(
      m_sparse_values.begin() + entry * m_sparse_lane_count,
      m_sparse_lane_count, 0.0);
  for (std::size_t sparse = 0; sparse < m_sparse_lane_count; ++sparse) {
    m_sparse_values[entry * m_sparse_lane_count + sparse] =
        implicitSparseValue(m_sparse_lane_begin + sparse, row);
  }
}

void GasCompositionSidecar::planSparseRowCapacity(
    RetainedCapacityTransaction& plan, std::size_t target_rows) {
  if (!sparseStorageEnabled() || target_rows <= sparseRowCapacity()) return;
  if (target_rows > m_sparse_gas_cell_id.max_size() ||
      target_rows > m_sparse_values.max_size() / m_sparse_lane_count) {
    throw std::length_error("GasCompositionSidecar sparse capacity target is too large");
  }
  plan.add(m_sparse_gas_cell_id, target_rows);
  plan.add(m_sparse_values, target_rows * m_sparse_lane_count);
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
  for (std::size_t lane_index = 0; lane_index < m_lanes.size(); ++lane_index) {
    if (!isSparseLane(lane_index)) m_lanes[lane_index][row] *= factor;
  }
  if (sparseStorageEnabled()) {
    const std::size_t entry = sparseEntryForId(gas_cell_id[row]);
    if (entry != m_sparse_gas_cell_id.size()) {
      for (std::size_t sparse = 0; sparse < m_sparse_lane_count; ++sparse) {
        m_sparse_values[entry * m_sparse_lane_count + sparse] *= factor;
      }
    }
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
  for (std::size_t lane_index = 0; lane_index < m_lanes.size(); ++lane_index) {
    payload.push_back(value(lane_index, row));
  }
  payload.push_back(static_cast<double>(chemistry_model_tag[row]));
  payload.push_back(static_cast<double>(ism_regime_tag[row]));
  payload.push_back(model_dwell_years[row]);
  payload.push_back(ism_regime_dwell_years[row]);
  payload.push_back(chemistry_last_temperature_k[row]);
  payload.push_back(chemistry_last_cooling_time_seconds[row]);
  payload.push_back(chemistry_last_error_estimate[row]);
  payload.push_back(chemistry_last_element_fraction_residual[row]);
  payload.push_back(chemistry_last_charge_residual[row]);
  payload.push_back(chemistry_cumulative_radiative_loss_erg[row]);
  payload.push_back(chemistry_cumulative_heating_gain_erg[row]);
  payload.push_back(hasExplicitSparseRow(row) ? 1.0 : 0.0);
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
  gas_cell_id[row] = gas_cell_id_value;
  const std::size_t meta = m_lanes.size();
  if (sparseStorageEnabled() && payload[meta + 11U] != 0.0) activateSparseRow(row);
  for (std::size_t l = 0; l < m_lanes.size(); ++l) setValue(l, row, payload[l]);
  chemistry_model_tag[row] = static_cast<std::uint8_t>(payload[meta]);
  ism_regime_tag[row] = static_cast<std::uint8_t>(payload[meta + 1U]);
  model_dwell_years[row] = payload[meta + 2U];
  ism_regime_dwell_years[row] = payload[meta + 3U];
  chemistry_last_temperature_k[row] = payload[meta + 4U];
  chemistry_last_cooling_time_seconds[row] = payload[meta + 5U];
  chemistry_last_error_estimate[row] = payload[meta + 6U];
  chemistry_last_element_fraction_residual[row] = payload[meta + 7U];
  chemistry_last_charge_residual[row] = payload[meta + 8U];
  chemistry_cumulative_radiative_loss_erg[row] = payload[meta + 9U];
  chemistry_cumulative_heating_gain_erg[row] = payload[meta + 10U];
  mean_molecular_weight[row] = std::numeric_limits<double>::quiet_NaN();
}

GasCompositionSidecar GasCompositionSidecar::emptyLike(
    std::size_t row_count, bool preserve_sparse_storage) const {
  GasCompositionSidecar result;
  if (enabled()) {
    result.configure(
        m_descriptor, row_count,
        preserve_sparse_storage ? m_sparse_lane_begin : GasCompositionLayoutDescriptor::k_no_lane,
        preserve_sparse_storage
            ? std::span<const std::uint32_t>(m_sparse_default_source_lanes)
            : std::span<const std::uint32_t>{});
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
      chemistry_last_temperature_k.size() != expected_rows ||
      chemistry_last_cooling_time_seconds.size() != expected_rows ||
      chemistry_last_error_estimate.size() != expected_rows ||
      chemistry_last_element_fraction_residual.size() != expected_rows ||
      chemistry_last_charge_residual.size() != expected_rows ||
      chemistry_cumulative_radiative_loss_erg.size() != expected_rows ||
      chemistry_cumulative_heating_gain_erg.size() != expected_rows ||
      mean_molecular_weight.size() != expected_rows) {
    return false;
  }
  for (std::size_t lane_index = 0; lane_index < m_lanes.size(); ++lane_index) {
    const auto& values = m_lanes[lane_index];
    if (values.size() != (isSparseLane(lane_index) ? 0U : expected_rows)) return false;
    for (const double v : values) {
      if (!std::isfinite(v) || v < 0.0) return false;
    }
  }
  if (m_sparse_default_source_lanes.size() != m_sparse_lane_count ||
      m_sparse_values.size() != m_sparse_gas_cell_id.size() * m_sparse_lane_count ||
      !std::is_sorted(m_sparse_gas_cell_id.begin(), m_sparse_gas_cell_id.end()) ||
      std::adjacent_find(m_sparse_gas_cell_id.begin(), m_sparse_gas_cell_id.end()) !=
          m_sparse_gas_cell_id.end()) {
    return false;
  }
  for (std::size_t entry = 0; entry < m_sparse_gas_cell_id.size(); ++entry) {
    if (std::find(gas_cell_id.begin(), gas_cell_id.end(), m_sparse_gas_cell_id[entry]) ==
        gas_cell_id.end()) return false;
    for (std::size_t sparse = 0; sparse < m_sparse_lane_count; ++sparse) {
      const double v = m_sparse_values[entry * m_sparse_lane_count + sparse];
      if (!std::isfinite(v) || v < 0.0) return false;
    }
  }
  for (const auto* values : {
           &model_dwell_years,
           &ism_regime_dwell_years,
           &chemistry_last_temperature_k,
           &chemistry_last_cooling_time_seconds,
           &chemistry_last_error_estimate,
           &chemistry_last_element_fraction_residual,
           &chemistry_last_charge_residual,
           &chemistry_cumulative_radiative_loss_erg,
           &chemistry_cumulative_heating_gain_erg}) {
    for (const double value : *values) {
      if (!std::isfinite(value) || value < 0.0) return false;
    }
  }
  return true;
}

std::uint64_t GasCompositionSidecar::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (const auto& values : m_lanes) bytes += values.capacity() * sizeof(double);
  bytes += m_sparse_default_source_lanes.capacity() * sizeof(std::uint32_t);
  bytes += m_sparse_gas_cell_id.capacity() * sizeof(std::uint64_t);
  bytes += m_sparse_values.capacity() * sizeof(double);
  bytes += gas_cell_id.capacity() * sizeof(std::uint64_t);
  bytes += chemistry_model_tag.capacity() + ism_regime_tag.capacity();
  bytes += (model_dwell_years.capacity() + ism_regime_dwell_years.capacity() +
            chemistry_last_temperature_k.capacity() +
            chemistry_last_cooling_time_seconds.capacity() +
            chemistry_last_error_estimate.capacity() +
            chemistry_last_element_fraction_residual.capacity() +
            chemistry_last_charge_residual.capacity() +
            chemistry_cumulative_radiative_loss_erg.capacity() +
            chemistry_cumulative_heating_gain_erg.capacity() +
            mean_molecular_weight.capacity()) * sizeof(double);
  return bytes;
}

std::uint64_t GasCompositionSidecar::currentSizeBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (const auto& values : m_lanes) bytes += values.size() * sizeof(double);
  bytes += m_sparse_default_source_lanes.size() * sizeof(std::uint32_t);
  bytes += m_sparse_gas_cell_id.size() * sizeof(std::uint64_t);
  bytes += m_sparse_values.size() * sizeof(double);
  bytes += gas_cell_id.size() * sizeof(std::uint64_t);
  bytes += chemistry_model_tag.size() + ism_regime_tag.size();
  bytes += (model_dwell_years.size() + ism_regime_dwell_years.size() +
            chemistry_last_temperature_k.size() +
            chemistry_last_cooling_time_seconds.size() +
            chemistry_last_error_estimate.size() +
            chemistry_last_element_fraction_residual.size() +
            chemistry_last_charge_residual.size() +
            chemistry_cumulative_radiative_loss_erg.size() +
            chemistry_cumulative_heating_gain_erg.size() +
            mean_molecular_weight.size()) * sizeof(double);
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
