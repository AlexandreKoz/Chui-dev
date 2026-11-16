#pragma once

// C04: optional, non-authoritative finite-volume hydrodynamics comparison.
// Every sample is a private conservative state trajectory. References use
// the existing CTU/PPM/HLLC(+HLLE fallback) solver, not production output
// replay. A two-half-step reference is used when production itself is CTU/PPM;
// CTU/PLM production is compared with one higher-order CTU/PPM step.
// Only cells strictly inside the numerical domain of dependence of the
// captured step-start stencil are compared. This measures local *transport*
// disagreement, never global physical error or the accuracy of source terms.

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {

struct HydroShadowOptions {
  bool enabled = false;
  std::uint64_t cadence = 16U;
  std::uint64_t seed = 0U;
  // Hard local-rank budgets; no population-scale persistent shadow state.
  std::size_t max_patches_per_rank = 2U;
  std::size_t max_storage_cells_per_rank = 8192U;
  std::uint64_t max_workspace_bytes_per_rank = 64U * 1024U * 1024U;
};

struct HydroShadowMetric {
  double weighted_abs_difference = 0.0;
  double weighted_reference_magnitude = 0.0;
  double linf_difference = 0.0;
};

struct HydroShadowReport {
  std::uint64_t eligible = 0U;
  std::uint64_t selected = 0U;
  std::uint64_t completed = 0U;
  std::uint64_t compared_cells = 0U;
  std::uint64_t compared_auxiliary_cells = 0U;
  std::uint64_t skipped_cadence = 0U;
  std::uint64_t skipped_selection = 0U;
  std::uint64_t skipped_geometry = 0U;
  std::uint64_t skipped_sources_or_time = 0U;
  std::uint64_t skipped_resource = 0U;
  std::uint64_t skipped_numerical = 0U;
  std::uint64_t reference_substeps = 0U;
  std::uint64_t fallback_count = 0U;
  std::uint64_t workspace_high_water_bytes = 0U; // largest admitted reservation
  std::uint64_t observed_capacity_high_water_bytes = 0U;
  double compared_volume = 0.0;
  std::array<HydroShadowMetric, 7> components{}; // rho,mx,my,mz,E,rhoZ,rho e
  std::array<double, 6> reference_delta{};
  std::array<double, 6> production_delta{};
};

class HydroShadowReference final {
 public:
  class Sample final {
   public:
    Sample() = default;
    Sample(const Sample&) = delete;
    Sample& operator=(const Sample&) = delete;
    Sample(Sample&&) noexcept = default;
    Sample& operator=(Sample&&) noexcept = default;
    [[nodiscard]] explicit operator bool() const noexcept { return valid; }
   private:
    friend class HydroShadowReference;
    // Lease declared before state: state destructs before the lease is released.
    core::MemoryReservation lease;
    HydroConservedStateSoa initial;
    std::uint64_t admitted_bytes = 0U;
    bool valid = false;
  };

  explicit HydroShadowReference(HydroShadowOptions options, core::MemoryGovernor* governor = nullptr);
  void beginEpoch(std::uint64_t step_index);
  [[nodiscard]] Sample capture(
      std::uint64_t patch_id, int level,
      const HydroConservedStateSoa& before,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active,
      const HydroUpdateContext& update,
      const HydroNumericsPolicy& production_numerics,
      HydroReconstructionScheme production_reconstruction,
      std::size_t source_count,
      bool effective_closure_active);
  void compare(
      Sample&& sample,
      const HydroConservedStateSoa& production_after,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active,
      const HydroUpdateContext& update,
      const HydroNumericsPolicy& production_numerics,
      HydroReconstructionScheme production_reconstruction,
      double adiabatic_index,
      double density_floor,
      double pressure_floor);
  [[nodiscard]] const HydroShadowReport& report() const noexcept { return m_report; }
  [[nodiscard]] const HydroShadowOptions& options() const noexcept { return m_options; }

 private:
  HydroShadowOptions m_options;
  core::MemoryGovernor* m_governor = nullptr;
  std::uint64_t m_step_index = 0U;
  std::size_t m_storage_cells_charged = 0U;
  std::uint64_t m_workspace_bytes_charged = 0U;
  HydroShadowReport m_report{};
};

}  // namespace cosmosim::hydro
