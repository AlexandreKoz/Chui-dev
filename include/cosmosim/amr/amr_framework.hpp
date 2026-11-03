#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory_resource>
#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::amr {

struct ConservedState {
  double mass_code = 0.0;
  double momentum_x_code = 0.0;
  double momentum_y_code = 0.0;
  double momentum_z_code = 0.0;
  double total_energy_code = 0.0;
  double metal_mass_code = 0.0;

  ConservedState& operator+=(const ConservedState& rhs);
  ConservedState& operator-=(const ConservedState& rhs);
  ConservedState& operator*=(double factor);
};

[[nodiscard]] ConservedState operator+(ConservedState lhs, const ConservedState& rhs);
[[nodiscard]] ConservedState operator-(ConservedState lhs, const ConservedState& rhs);
[[nodiscard]] ConservedState operator*(ConservedState lhs, double factor);

struct CellMetrics {
  double density_code = 0.0;
  double pressure_code = 0.0;
  double sound_speed_code = 0.0;
  double gradient_indicator = 0.0;
  std::uint32_t particle_count = 0;
};

enum class RefinementDecision : std::uint8_t {
  kKeep,
  kRefine,
  kDerefine,
};

struct RefinementCriteria {
  bool use_mass_threshold = true;
  bool use_gradient_indicator = true;
  bool use_particle_count = true;
  bool use_jeans_condition = true;

  double mass_threshold_code = 1.0;
  double gradient_threshold = 0.5;
  std::uint32_t particle_threshold = 8;
  double jeans_resolution_cells = 4.0;
  double gravitational_constant_code = 1.0;
  double derefine_hysteresis = 0.5;
};

class RefinementEvaluator {
 public:
  [[nodiscard]] static RefinementDecision evaluateCell(
      const CellMetrics& metrics,
      double cell_width_comov,
      const RefinementCriteria& criteria);
};

struct PatchDescriptor {
  std::uint64_t patch_id = 0;
  std::uint64_t parent_patch_id = 0;
  std::uint8_t level = 0;
  std::uint64_t morton_key = 0;
  std::array<double, 3> origin_comov = {0.0, 0.0, 0.0};
  std::array<double, 3> extent_comov = {1.0, 1.0, 1.0};
  std::array<std::uint16_t, 3> cell_dims = {4, 4, 4};
};

class AmrPatch {
 public:
  explicit AmrPatch(PatchDescriptor descriptor);

  [[nodiscard]] const PatchDescriptor& descriptor() const;
  [[nodiscard]] std::size_t cellCount() const;
  [[nodiscard]] double cellVolumeComov() const;
  [[nodiscard]] double cellWidthComov() const;

  [[nodiscard]] std::span<ConservedState> conservedView();
  [[nodiscard]] std::span<const ConservedState> conservedView() const;
  [[nodiscard]] std::span<CellMetrics> metricsView();
  [[nodiscard]] std::span<const CellMetrics> metricsView() const;
  [[nodiscard]] std::span<std::uint64_t> gasCellIdView();
  [[nodiscard]] std::span<const std::uint64_t> gasCellIdView() const;

  [[nodiscard]] ConservedState totalConserved() const;
  // Capacity truth for the patch-owned physical state lanes only. This excludes
  // allocator bookkeeping and hierarchy/index container overhead.
  [[nodiscard]] std::size_t ownedStorageCapacityBytes() const;
  [[nodiscard]] bool isLeaf() const;
  void setLeaf(bool is_leaf);

 private:
  PatchDescriptor m_descriptor;
  std::vector<ConservedState> m_conserved;
  std::vector<CellMetrics> m_metrics;
  std::vector<std::uint64_t> m_gas_cell_ids;
  bool m_is_leaf = true;
};

class PatchHierarchy {
 public:
  [[nodiscard]] std::uint64_t createRootPatch(const PatchDescriptor& root);
  [[nodiscard]] std::array<std::uint64_t, 8> refinePatch(std::uint64_t patch_id);
  bool derefinePatch(std::uint64_t parent_patch_id);

  [[nodiscard]] AmrPatch* findPatch(std::uint64_t patch_id);
  [[nodiscard]] const AmrPatch* findPatch(std::uint64_t patch_id) const;

  [[nodiscard]] std::span<AmrPatch> levelView(std::size_t level);
  [[nodiscard]] std::span<const AmrPatch> levelView(std::size_t level) const;

  [[nodiscard]] std::size_t levelCount() const;
  [[nodiscard]] std::size_t patchCount() const;
  [[nodiscard]] std::size_t allocatedCellCount() const;
  [[nodiscard]] std::size_t patchStorageCapacityBytes() const;
  [[nodiscard]] std::span<const std::uint64_t> retiredGasCellIds() const;

 private:
  std::vector<std::vector<AmrPatch>> m_levels;
  std::unordered_map<std::uint64_t, std::pair<std::size_t, std::size_t>> m_patch_index;
  std::uint64_t m_next_patch_id = 1;
  std::uint64_t m_next_gas_cell_id = 1;
  std::vector<std::uint64_t> m_retired_gas_cell_ids;

  static void assignStableGasCellIds(AmrPatch& patch, std::uint64_t& next_gas_cell_id);
  void rebuildPatchIndex();
};

class ConservativeTransfer {
 public:
  [[nodiscard]] static std::vector<ConservedState> prolongateFromCoarse(
      const ConservedState& coarse,
      std::size_t fine_cell_count);

  [[nodiscard]] static ConservedState restrictToCoarse(std::span<const ConservedState> fine_cells);
};

// Local synchronized-sweep flux-register record. Incomplete records are safety
// rejects in applyFluxRegistersToSimulationState(): they are not persisted and
// not replayed later. A future Berger-Colella subcycling implementation needs a
// restart-safe pending-register owner with patch/gas identity metadata, area
// coverage, and dt coverage before deferred synchronization can be claimed.
struct FluxRegisterEntry {
  std::uint64_t register_key = 0;
  std::uint64_t coarse_patch_id = 0;
  std::uint64_t coarse_gas_cell_id = 0;
  std::size_t coarse_cell_index = 0;
  std::uint8_t level = 0;
  hydro::HydroFaceAxis axis = hydro::HydroFaceAxis::kX;
  hydro::HydroFaceSide orientation = hydro::HydroFaceSide::kLower;
  ConservedState coarse_face_flux_code;
  ConservedState fine_face_flux_code;
  double face_area_comov = 1.0;
  double coarse_area_comov = 0.0;
  double fine_area_comov = 0.0;
  double dt_code = 0.0;
  std::uint32_t coarse_face_count = 0;
  std::uint32_t fine_face_count = 0;
  [[nodiscard]] bool isComplete() const noexcept { return coarse_face_count > 0U && fine_face_count > 0U; }
};

class FluxRegisterAccumulator final : public hydro::HydroFluxRegisterSink {
 public:
  explicit FluxRegisterAccumulator(std::pmr::memory_resource* resource =
                                       std::pmr::get_default_resource());
  // Production may predeclare the complete set of legal register keys. The
  // fixed lookup and accumulators are then allocated before any flux update;
  // recordFaceFlux performs no allocation and preserves per-key accumulation
  // order. The ordinary dynamically keyed compatibility path is unchanged.
  [[nodiscard]] static std::uint64_t fixedWorkspaceBytes(std::size_t key_count);
  void prepareFixedKeys(std::span<const std::uint64_t> sorted_unique_keys);
  [[nodiscard]] std::size_t entryCount() const noexcept;
  void recordFaceFlux(const hydro::HydroFluxRegisterRecord& record) override;
  [[nodiscard]] std::vector<FluxRegisterEntry> entries() const;
  void clear();

 private:
  struct AccumulatedEntry {
    FluxRegisterEntry entry;
    ConservedState coarse_area_weighted_flux;
    ConservedState fine_area_weighted_flux;
    double coarse_area_comov = 0.0;
    double fine_area_comov = 0.0;
  };

  std::pmr::unordered_map<std::uint64_t, std::size_t> m_slot_by_key;
  std::pmr::vector<AccumulatedEntry> m_entries;
  std::pmr::vector<std::uint64_t> m_fixed_keys;
  bool m_fixed = false;
};

struct RefluxDiagnostics {
  std::size_t complete_register_count = 0;
  std::size_t skipped_incomplete_register_count = 0;
  std::size_t skipped_area_mismatch_count = 0;
  std::size_t skipped_missing_target_count = 0;
  std::size_t corrected_cells = 0;
  double corrected_mass_code = 0.0;
  double corrected_momentum_x_code = 0.0;
  double corrected_momentum_y_code = 0.0;
  double corrected_momentum_z_code = 0.0;
  double corrected_total_energy_code = 0.0;
  double corrected_metal_mass_code = 0.0;
  double corrected_energy_code = 0.0;
  double corrected_internal_energy_code = 0.0;
  // Dual energy: corrected cells whose post-reflux thermal fraction was below
  // eta kept their pre-reflux specific internal energy; the resulting
  // non-conservative thermal correction sum (rho e_kept - (E - K)) V is
  // ledgered here (code energy).
  std::size_t dual_energy_reflux_cells = 0;
  double dual_energy_reflux_thermal_correction_code = 0.0;
};

class RefluxSynchronizer {
 public:
  [[nodiscard]] static RefluxDiagnostics apply(
      PatchHierarchy& hierarchy,
      std::span<const FluxRegisterEntry> entries);
};

struct RegridDiagnostics {
  std::size_t refined_patch_count = 0;
  std::size_t derefined_patch_count = 0;
  std::size_t touched_leaf_patch_count = 0;
  std::uint64_t elapsed_microseconds = 0;
};

class RefineDerefineManager {
 public:
  explicit RefineDerefineManager(RefinementCriteria criteria);

  [[nodiscard]] RegridDiagnostics regrid(PatchHierarchy& hierarchy) const;

 private:
  RefinementCriteria m_criteria;

  [[nodiscard]] RefinementDecision evaluatePatchDecision(const AmrPatch& patch) const;
};

}  // namespace cosmosim::amr
