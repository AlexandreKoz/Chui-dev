#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/soa_storage.hpp"

namespace cosmosim::core {

class RetainedCapacityTransaction;

// CHUI Part 2: optional per-gas-cell composition sidecar.
//
// Ownership contract
// ------------------
// * The sidecar is GENERIC storage owned by SimulationState. Its physical
//   meaning (which lane is which element / species / dust reservoir) is fixed
//   by a physics-side manifest (cosmosim::physics::CompositionLayout) whose
//   signature is stored in the descriptor. core never interprets lanes.
// * Every lane stores an EXTENSIVE mass in code mass units (the same unit as
//   CellSoa::mass_code). Extensive storage makes advection, migration,
//   star-formation removal and feedback injection additive and exactly
//   conservative; intensive fractions are derived by dividing by cell mass.
//   The single exception is the optional SGS lane, which stores the
//   extensive SGS turbulent energy in code energy units (mass_code *
//   velocity_code^2).
// * Rows are aligned with gas rows (CellSoa / GasCellSidecar). The
//   gas_cell_id lane is an identity mirror used to prove row alignment
//   before every physics consumer (audit risk: cooling/composition row
//   permutation). A mismatch is a hard error, never a silent remap.
// * Disabled (layout_version == 0) means ZERO allocation: DMO and plain hydro
//   runs carry no lanes.
// * chemistry/ISM tags and dwell counters plus the accepted local
//   temperature, cooling-time, integration-error, conservation-residual and
//   radiative/heating ledgers are persistent history for model selection.
// * mean_molecular_weight is a DERIVED cache (NaN = unknown) refreshed by the
//   thermochemistry commit and consumed by the gas thermodynamic refresh to
//   convert specific internal energy into kelvin. It is not restart truth;
//   it is rebuilt from the lanes after restart.
// ISM regime tags (GasCompositionSidecar::ism_regime_tag). Owned by the
// physics ISM regime manager; consumed by hydro closures, thermodynamic
// refresh, thermochemistry and star formation so exactly one regime owns a
// cell's thermal physics.
inline constexpr std::uint8_t k_gas_ism_regime_resolved = 0U;
inline constexpr std::uint8_t k_gas_ism_regime_effective = 1U;

struct GasCompositionLayoutDescriptor {
  // 0 = disabled. Bumped by the physics manifest when lane meaning changes.
  std::uint32_t layout_version = 0;
  // Stable 64-bit signature over the ordered lane names/units.
  std::uint64_t layout_signature = 0;
  std::string layout_name;
  std::uint32_t lane_count = 0;
  std::uint32_t element_lane_begin = 0;
  std::uint32_t element_lane_count = 0;
  std::uint32_t species_lane_begin = 0;
  std::uint32_t species_lane_count = 0;
  std::uint32_t dust_lane_begin = 0;
  std::uint32_t dust_lane_count = 0;
  static constexpr std::uint32_t k_no_lane = std::numeric_limits<std::uint32_t>::max();
  std::uint32_t sgs_energy_lane = k_no_lane;

  [[nodiscard]] bool enabled() const noexcept { return layout_version != 0U && lane_count != 0U; }
  [[nodiscard]] bool isConsistent() const noexcept;
  [[nodiscard]] bool operator==(const GasCompositionLayoutDescriptor& rhs) const noexcept;
};

class GasCompositionSidecar {
 public:
  // Number of trailing metadata values appended by packRow(). The first four
  // retain the schema-1 order; schema-2 adaptive-history values and the
  // schema-3 explicit sparse-history marker follow.
  static constexpr std::size_t k_packed_metadata_count = 12U;

  // sparse_default_source_lanes selects a contiguous sparse lane suffix-like
  // range beginning at sparse_lane_begin. Each entry is either a dense lane
  // whose current row value is the implicit value of an inactive sparse row,
  // or k_no_lane for implicit zero. This keeps adaptive detailed chemistry
  // sparse without discarding its physically defined C+/CO projection.
  void configure(
      const GasCompositionLayoutDescriptor& descriptor,
      std::size_t row_count,
      std::uint32_t sparse_lane_begin = GasCompositionLayoutDescriptor::k_no_lane,
      std::span<const std::uint32_t> sparse_default_source_lanes = {});
  void clear() noexcept;
  [[nodiscard]] bool enabled() const noexcept { return m_descriptor.enabled(); }
  [[nodiscard]] const GasCompositionLayoutDescriptor& descriptor() const noexcept {
    return m_descriptor;
  }
  [[nodiscard]] std::size_t rowCount() const noexcept { return gas_cell_id.size(); }
  [[nodiscard]] std::size_t laneCount() const noexcept { return m_lanes.size(); }
  // New rows are zero-mass with unknown mu (NaN) and gas_cell_id 0.
  void resizeRows(std::size_t row_count);

  [[nodiscard]] std::span<double> lane(std::size_t lane_index);
  [[nodiscard]] std::span<const double> lane(std::size_t lane_index) const;
  [[nodiscard]] bool laneStoredDensely(std::size_t lane_index) const noexcept;
  [[nodiscard]] double value(std::size_t lane_index, std::size_t row) const;
  void setValue(std::size_t lane_index, std::size_t row, double value);
  void addValue(std::size_t lane_index, std::size_t row, double delta);

  // Element-lane total (sum of the element lanes) for one row.
  [[nodiscard]] double elementMassSum(std::size_t row) const;
  // Scale every extensive lane of a row by factor (uniform mass removal,
  // e.g. star formation or accretion). factor must be in [0, 1].
  void scaleRow(std::size_t row, double factor);

  // Migration payload: lane values followed by metadata values
  // (chemistry_model_tag, ism_regime_tag, model_dwell_years,
  // ism_regime_dwell_years). Empty when the
  // sidecar is disabled.
  [[nodiscard]] std::vector<double> packRow(std::size_t row) const;
  void unpackRow(std::size_t row, std::span<const double> payload, std::uint64_t gas_cell_id_value);
  [[nodiscard]] std::size_t packedRowWidth() const noexcept {
    return enabled() ? laneCount() + k_packed_metadata_count : 0U;
  }

  // Same layout, row_count zero-initialized rows (used by rebuilds).
  [[nodiscard]] GasCompositionSidecar emptyLike(
      std::size_t row_count, bool preserve_sparse_storage = true) const;

  // Adaptive detailed-state storage. Entries are sorted by stable gas_cell_id
  // and share one packed row-major value allocation (never a map/vector per
  // gas cell). Capacity growth must be planned before a no-fail commit.
  [[nodiscard]] bool sparseStorageEnabled() const noexcept {
    return m_sparse_lane_count != 0U;
  }
  [[nodiscard]] std::uint32_t sparseLaneBegin() const noexcept {
    return m_sparse_lane_begin;
  }
  [[nodiscard]] std::uint32_t sparseLaneCount() const noexcept {
    return m_sparse_lane_count;
  }
  [[nodiscard]] std::size_t sparseActiveRowCount() const noexcept {
    return m_sparse_gas_cell_id.size();
  }
  [[nodiscard]] std::size_t sparseRowCapacity() const noexcept {
    if (m_sparse_lane_count == 0U) return 0U;
    const std::size_t value_rows = m_sparse_values.capacity() / m_sparse_lane_count;
    return m_sparse_gas_cell_id.capacity() < value_rows
        ? m_sparse_gas_cell_id.capacity()
        : value_rows;
  }
  [[nodiscard]] bool hasExplicitSparseRow(std::size_t row) const;
  [[nodiscard]] double implicitSparseValue(std::size_t lane_index, std::size_t row) const;
  [[nodiscard]] std::uint32_t sparseDefaultSourceLane(std::size_t lane_index) const;
  void activateSparseRow(std::size_t row);
  void planSparseRowCapacity(RetainedCapacityTransaction& plan, std::size_t target_rows);

  // Row identity guard: throws std::runtime_error with a precise diagnostic
  // when gas_cell_id mirrors differ from expected (row permutation defect).
  void requireAlignedWith(std::span<const std::uint64_t> expected_gas_cell_ids,
                          const char* caller) const;
  // Finite, non-negative lanes (SGS lane included) and size coherence.
  [[nodiscard]] bool isConsistent(std::size_t expected_rows) const noexcept;

  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
  [[nodiscard]] std::uint64_t currentSizeBytes() const noexcept;

  AlignedVector<std::uint64_t> gas_cell_id;
  AlignedVector<std::uint8_t> chemistry_model_tag;
  AlignedVector<std::uint8_t> ism_regime_tag;
  AlignedVector<double> model_dwell_years;
  // Physical time spent in the current ISM regime (hysteresis / dwell rule).
  AlignedVector<double> ism_regime_dwell_years;
  // C03 local, restart-authoritative solver history. Zero means unavailable
  // for the last-* fields. Cumulative radiative/heating energies are
  // extensive and follow conservative AMR prolongation/restriction.
  AlignedVector<double> chemistry_last_temperature_k;
  AlignedVector<double> chemistry_last_cooling_time_seconds;
  AlignedVector<double> chemistry_last_error_estimate;
  AlignedVector<double> chemistry_last_element_fraction_residual;
  AlignedVector<double> chemistry_last_charge_residual;
  AlignedVector<double> chemistry_cumulative_radiative_loss_erg;
  AlignedVector<double> chemistry_cumulative_heating_gain_erg;
  AlignedVector<double> mean_molecular_weight;

 private:
  [[nodiscard]] bool isSparseLane(std::size_t lane_index) const noexcept;
  [[nodiscard]] std::size_t sparseEntryForId(std::uint64_t id) const noexcept;
  [[nodiscard]] std::size_t requireRow(std::size_t row) const;
  GasCompositionLayoutDescriptor m_descriptor;
  // One allocation per DENSE lane, never per cell. Sparse lanes retain an
  // empty vector at their manifest index so lane numbers never change.
  std::vector<AlignedVector<double>> m_lanes;
  std::uint32_t m_sparse_lane_begin = GasCompositionLayoutDescriptor::k_no_lane;
  std::uint32_t m_sparse_lane_count = 0U;
  std::vector<std::uint32_t> m_sparse_default_source_lanes;
  AlignedVector<std::uint64_t> m_sparse_gas_cell_id;
  AlignedVector<double> m_sparse_values;
};

// FNV-1a 64-bit over a byte string: used for layout signatures only (not a
// security or data-integrity hash; data files use SHA-256).
[[nodiscard]] std::uint64_t fnv1a64(std::string_view text) noexcept;

}  // namespace cosmosim::core
