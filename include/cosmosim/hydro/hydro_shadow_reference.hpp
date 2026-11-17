#pragma once

// C04: sampled, non-authoritative finite-volume comparison. A sample is a
// separately advanced conservative trajectory, never committed to production.
// CTU/PPM vs CTU/PPM halfsteps changes only temporal subdivision; other
// production schemes vs CTU/PPM are explicitly mixed-scheme comparisons.
// Local smooth source operators must provide a private/replay-safe adapter;
// discrete events, external composition transport and unavailable AMR temporal
// ghosts are independently counted ineligibilities, not error=0 samples.

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <memory>
#include <span>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {

struct HydroShadowOptions {
  bool enabled = false;
  std::uint64_t cadence = 16U;
  std::uint64_t seed = 0U;
  // Legacy config key hydro_shadow_max_patches_per_rank is interpreted as
  // a GLOBAL top-K bound by the v10 runtime. Per-rank storage/bytes remain local.
  std::size_t max_patches_per_rank = 2U;
  std::size_t max_storage_cells_per_rank = 8192U;
  std::uint64_t max_workspace_bytes_per_rank = 64U * 1024U * 1024U;
};

// Stable physics identity; independent of MPI rank and patch-vector order.
struct HydroShadowCandidate {
  std::uint64_t patch_id = 0;
  int level = 0;
};

// Global deterministic top-K comparator. Ranks submit at most K local
// candidates; taking top K again after allgather preserves global top K.
[[nodiscard]] bool hydroShadowCandidateLess(
    const HydroShadowCandidate& lhs, const HydroShadowCandidate& rhs,
    std::uint64_t step_index, std::uint64_t seed) noexcept;
[[nodiscard]] std::vector<HydroShadowCandidate> rankHydroShadowCandidates(
    std::span<const HydroShadowCandidate> candidates, std::size_t k,
    std::uint64_t step_index, std::uint64_t seed);

struct HydroShadowMetric {
  double weighted_abs_difference = 0.0;
  double weighted_reference_magnitude = 0.0;
  double linf_difference = 0.0;
};

struct HydroShadowReport {
  std::uint64_t eligible = 0U;
  std::uint64_t selected = 0U;
  std::uint64_t completed = 0U;
  std::uint64_t completed_closed_boundary_patches = 0U;
  std::uint64_t completed_source_active_patches = 0U;
  std::uint64_t completed_evolving_comoving_patches = 0U;
  std::uint64_t compared_cells = 0U;
  std::uint64_t compared_auxiliary_cells = 0U;
  std::uint64_t skipped_cadence = 0U;
  std::uint64_t skipped_selection = 0U;
  std::uint64_t skipped_geometry = 0U;
  std::uint64_t skipped_sources_or_time = 0U;
  std::uint64_t skipped_resource = 0U;
  std::uint64_t skipped_numerical = 0U;
  std::uint64_t skipped_unsupported_source = 0U;
  std::uint64_t skipped_external_passives = 0U;
  std::uint64_t skipped_temporal_ghost = 0U;
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
    HydroReconstructionPolicy reconstruction_policy{};
    HydroRiemannSolverKind riemann_kind = HydroRiemannSolverKind::kHllc;
    std::vector<std::unique_ptr<HydroSourceTerm>> shadow_sources;
    std::vector<HydroCellSourceProperties> captured_properties;
    std::vector<std::uint8_t> effective_closure_allowed;
    const HydroThermodynamicClosure* closure = nullptr;
    double source_redshift = 0.0;
    bool source_properties_captured = false;
    bool closed_physical_boundary = false;
    bool source_active = false;
    bool evolving_comoving = false;
    std::uint64_t admitted_bytes = 0U;
    bool valid = false;
  };

  explicit HydroShadowReference(HydroShadowOptions options, core::MemoryGovernor* governor = nullptr);
  void beginEpoch(std::uint64_t step_index);
  // The coordinator must call this once at the rank-uniform hydro boundary
  // with global candidates computed by a bounded collective. No MPI inside
  // conditional patch advances.
  void setGlobalCandidates(std::vector<HydroShadowCandidate> ranked);
  void noteGlobalSelectionResourceRefusal() noexcept { ++m_report.skipped_resource; }
  [[nodiscard]] Sample capture(
      std::uint64_t patch_id, int level,
      const HydroConservedStateSoa& before,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active,
      const HydroUpdateContext& update,
      const HydroNumericsPolicy& production_numerics,
      HydroReconstructionScheme production_reconstruction,
      const HydroReconstructionPolicy& production_policy,
      HydroRiemannSolverKind production_riemann,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context);
  void compare(
      Sample&& sample,
      const HydroConservedStateSoa& production_after,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active,
      const HydroUpdateContext& update,
      const HydroNumericsPolicy& production_numerics,
      HydroReconstructionScheme production_reconstruction);
  [[nodiscard]] const HydroShadowReport& report() const noexcept { return m_report; }
  [[nodiscard]] const HydroShadowOptions& options() const noexcept { return m_options; }

 private:
  HydroShadowOptions m_options;
  core::MemoryGovernor* m_governor = nullptr;
  std::uint64_t m_step_index = 0U;
  std::size_t m_storage_cells_charged = 0U;
  std::uint64_t m_workspace_bytes_charged = 0U;
  HydroShadowReport m_report{};
  std::vector<HydroShadowCandidate> m_global_candidates;
  std::vector<HydroShadowCandidate> m_selected_once;
};

}  // namespace cosmosim::hydro
