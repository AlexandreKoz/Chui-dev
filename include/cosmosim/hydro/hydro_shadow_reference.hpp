#pragma once

// C04: sampled, non-authoritative finite-volume comparison. A sample is a
// separately advanced conservative trajectory, never committed to production.
//
// Comparison classes (HydroShadowComparisonClass):
//  * kTimeOnlyCtuPpm: production CTU/PPM vs two private CTU/PPM halfsteps.
//    Only the temporal subdivision differs.
//  * kMixedSpatialTemporal: another production reconstruction/integrator vs
//    one private CTU/PPM step (unlike spatial AND temporal schemes).
//  * kCoupledFullPhysicsStrang: the accepted T01 interval
//    S(dt/2) H(dt) S(dt/2) (production smooth operators, composition
//    transport, captured discrete-event forcing) vs the refined chronology
//    S(dt/4) H(dt/2) S(dt/4) S(dt/4) H(dt/2) S(dt/4) on a private region.
//  * kBoundedProjection: a comparison whose comparable cell set was reduced
//    by declared dependency-cone exclusions (stale ghost authority).
//
// Ghost authority (C04-03). Every ghost slot of a sampled patch is classified
// by provenance (self-reconstructible physical/periodic, imported same-level,
// remote MPI, edge/diagonal, coarse->fine, fine->coarse, closure-mapped) and
// by time authority for the private halfsteps:
//  * kSelfRefreshed: recomputed from the PRIVATE state before each halfstep
//    (in-patch periodic wrap, open, reflective) by hydro boundary rules;
//  * kEndpointInterpolated: the caller supplies the ghost at the production
//    interval end through provideEndpointGhosts(), formed by the SAME
//    production ghost filler (owner, prolongation, restriction, temporal
//    history, MPI payload). The halfstep-2 ghost is the linear time
//    interpolation of the step-start and endpoint values, i.e. the production
//    AMR temporal-ghost rule. No flux register or real halo is touched;
//  * kStepStartOnly: no endpoint authority exists at a legal point (e.g. a
//    recursive coarse level whose finer neighbour has not advanced). The
//    step-start value is kept and every real cell inside the halfstep
//    dependency cone of that slot (Chebyshev radius
//    k_hydro_shadow_halfstep_dependency_radius around the owner cell) is
//    excluded from the comparison and counted. No value is invented.
//
// Sampling identity (C04-04). Candidates carry a physical region key derived
// from (level, global lattice origin, extent). Ranking, de-duplication and the
// single GLOBAL K bound use only that key, the step index and the seed, so the
// selected physical regions are independent of the MPI decomposition. Patch
// IDs are rank-local lookup handles only.
//
// Resource contract (C04-05). Every private copy (state, ghosts, endpoint
// ghosts, production result, geometry/active copies, composition rows,
// chronology snapshots and scratch) is charged to one MemoryGovernor lease
// with checked arithmetic BEFORE allocation. A refusal never vetoes
// production evolution. Shadow work never writes authoritative state,
// registers, ledgers, RNG cursors, event queues or controller state.

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include "cosmosim/core/gas_composition_sidecar.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::core {
class LambdaCdmBackground;
}  // namespace cosmosim::core

namespace cosmosim::hydro {

// One CTU/PPM halfstep: PPM stencil (3 cells per side of a face) plus the
// CTU transverse increment of the adjacent cell (1 cell). Real cells within
// this Chebyshev distance of the owner of a stale ghost are not comparable.
inline constexpr int k_hydro_shadow_halfstep_dependency_radius = 4;
// Donor-cell composition transport moves one cell per halfstep (Courant <= 1).
inline constexpr int k_hydro_shadow_composition_dependency_radius = 1;
// Fixed-grid physical sampling tiles (cells per axis), equal to the default
// production AMR patch width. Not a config key: it defines region identity.
inline constexpr std::uint32_t k_hydro_shadow_fixed_grid_region_extent = 8U;

struct HydroShadowOptions {
  bool enabled = false;
  std::uint64_t cadence = 16U;
  std::uint64_t seed = 0U;
  // Legacy config key hydro_shadow_max_patches_per_rank is interpreted as
  // a GLOBAL top-K bound on physical regions. Per-rank storage/bytes stay local.
  std::size_t max_patches_per_rank = 2U;
  std::size_t max_storage_cells_per_rank = 8192U;
  std::uint64_t max_workspace_bytes_per_rank = 64U * 1024U * 1024U;
};

// Candidate identity. region_key is the physical identity; patch_id/level are
// the rank-local handle of the owning patch (patch_id is NOT used for ranking
// when region_key != 0). region_key == 0 is the legacy form (key derived from
// patch_id and level) kept for standalone callers.
struct HydroShadowCandidate {
  std::uint64_t patch_id = 0;
  int level = 0;
  std::uint64_t region_key = 0;
};

// Physical region key of a Cartesian block at refinement level `level` whose
// lower corner is the global lattice index `lattice_origin` (cells of that
// level, measured from the hierarchy/domain origin) with `extent` cells per
// axis. Never zero. Independent of rank, patch ID and container order.
[[nodiscard]] std::uint64_t hydroShadowPhysicalRegionKey(
    int level, const std::array<std::int64_t, 3>& lattice_origin,
    const std::array<std::uint32_t, 3>& extent) noexcept;
// Global lattice index of a block origin; nullopt when the origin is not an
// integer number of cell widths from the domain origin (within 1e-6 cells) or
// a width is not positive.
[[nodiscard]] std::optional<std::array<std::int64_t, 3>> hydroShadowLatticeOrigin(
    const std::array<double, 3>& origin_comoving,
    const std::array<double, 3>& domain_origin_comoving,
    const std::array<double, 3>& cell_width_comoving) noexcept;
[[nodiscard]] std::uint64_t hydroShadowEffectiveKey(const HydroShadowCandidate& candidate) noexcept;

// Global deterministic top-K comparator. Ranks submit at most K local
// candidates; taking top K again after allgather preserves global top K.
[[nodiscard]] bool hydroShadowCandidateLess(
    const HydroShadowCandidate& lhs, const HydroShadowCandidate& rhs,
    std::uint64_t step_index, std::uint64_t seed) noexcept;
[[nodiscard]] std::vector<HydroShadowCandidate> rankHydroShadowCandidates(
    std::span<const HydroShadowCandidate> candidates, std::size_t k,
    std::uint64_t step_index, std::uint64_t seed);

// Mirrors the deterministic fixed-iteration FLRW inversion of the source
// owner (64 bisection steps of the exact cosmic-time integral): the scale
// factor at physical-time fraction f of [a_begin, a_end]. Linear in a when
// background == nullptr; returns a_begin when a_begin == a_end.
[[nodiscard]] double hydroShadowScaleFactorAtPhysicalFraction(
    const core::LambdaCdmBackground* background, double scale_factor_begin,
    double scale_factor_end, double fraction);

enum class HydroShadowGhostProvenance : std::uint8_t {
  kSelfPhysical = 0,          // open/reflective owner-derived ghost of this patch
  kSelfPeriodicWrap = 1,      // in-patch periodic wrap of this patch
  kImportedSameLevel = 2,     // same-level ghost from a local neighbour patch
  kImportedSameLevelRemote = 3,  // same-level ghost from a remote MPI owner
  kImportedEdgeDiagonal = 4,  // edge/corner (diagonal) stencil ghost
  kCoarseToFine = 5,          // spatial prolongation from a coarser patch
  kCoarseToFineTemporal = 6,  // prolongation from a temporal boundary history
  kFineToCoarse = 7,          // restriction from finer patches
  kPhysicalClosureMapped = 8, // deep reflective/open closure resolved by the AMR filler
  kUnknownImported = 9,       // imported ghost without further metadata
};
inline constexpr std::size_t k_hydro_shadow_ghost_provenance_count = 10U;

enum class HydroShadowGhostTimeAuthority : std::uint8_t {
  kSelfRefreshed = 0,
  kEndpointInterpolated = 1,
  kStepStartOnly = 2,
};

enum class HydroShadowComparisonClass : std::uint8_t {
  kTimeOnlyCtuPpm = 0,
  kMixedSpatialTemporal = 1,
  kCoupledFullPhysicsStrang = 2,
  kBoundedProjection = 3,
};
inline constexpr std::size_t k_hydro_shadow_comparison_class_count = 4U;
[[nodiscard]] std::string_view hydroShadowComparisonClassToString(HydroShadowComparisonClass value) noexcept;
[[nodiscard]] std::string_view hydroShadowGhostProvenanceToString(HydroShadowGhostProvenance value) noexcept;

// Default provenance from solver geometry metadata only (fixed grid and
// standalone callers): writable physical ghosts are self-reconstructible;
// a periodic ghost whose wrapped source lies in the patch is a periodic wrap;
// edge ghosts are diagonal imports; every other imported ghost is unknown.
[[nodiscard]] HydroShadowGhostProvenance hydroShadowDefaultGhostProvenance(
    const HydroPatchGeometry& geometry, std::size_t ghost_slot) noexcept;

// Real cells of `geometry` at Chebyshev distance <= radius from any owner
// cell of a ghost slot with stale_slot[slot] != 0 (non-suppressed axes only).
// Returns one flag per real cell.
[[nodiscard]] std::vector<std::uint8_t> hydroShadowDependencyConeMask(
    const HydroPatchGeometry& geometry, std::span<const std::uint8_t> stale_slot, int radius);

// A comparison region inside a sample domain: the physical key and the real
// geometry rows it covers (empty = every real cell of the domain).
struct HydroShadowComparisonRegion {
  HydroShadowCandidate identity{};
  std::vector<std::size_t> cells;
};

// Private region state used by the source chronology. Real cells only, in
// sample-domain geometry order. Lanes are comoving code densities and the
// selected specific internal energy (dual-energy authority when present).
struct HydroShadowPrivateRegion {
  std::size_t cell_count = 0U;
  double cell_volume_comoving = 0.0;
  std::vector<std::uint64_t> gas_cell_ids;
  std::vector<double> mass_density;
  std::vector<double> momentum_x;
  std::vector<double> momentum_y;
  std::vector<double> momentum_z;
  std::vector<double> total_energy;
  std::vector<double> metal_density;
  std::vector<double> specific_internal_energy;
  core::GasCompositionSidecar composition;  // disabled when composition is off

  void resize(std::size_t count);
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
  // Checked upper bound of the lanes for `count` cells and `composition_lanes`
  // packed composition doubles per row (without the sidecar object overhead).
  [[nodiscard]] static std::uint64_t estimatedBytes(std::size_t count, std::size_t composition_lanes);
};

enum class HydroShadowReplayStatus : std::uint8_t {
  kCompleted = 0,
  kModelUnavailable = 1,       // operator cannot run privately (fail closed)
  kIntervalIncomplete = 2,     // the private integrator would require a retry
  kInvalidState = 3,
};

// Private, non-authoritative replay of the production SMOOTH source operator
// (thermochemistry, ISM regime, dust, SGS dissipation) over the physical
// fraction [fraction_begin, fraction_end] of the accepted coarse interval.
// Implemented by the source owner (workflows). It must use the production
// operator and data tables, must never touch SimulationState, ledgers,
// controller state, RNG or event queues, and must evaluate the interval with
// the production cosmological epoch rules.
class HydroShadowSmoothSourceReplay {
 public:
  virtual ~HydroShadowSmoothSourceReplay() = default;
  [[nodiscard]] virtual HydroShadowReplayStatus advanceSmooth(
      HydroShadowPrivateRegion& region, double fraction_begin, double fraction_end) const = 0;
  [[nodiscard]] virtual std::string_view modelLabel() const noexcept = 0;
  // Re-derives the production derived mirrors of a private region after a
  // private hydro substep (composition-authoritative runs overwrite the
  // passive metal lane from the element/dust inventory after transport).
  // Default: no derived mirror.
  virtual void synchronizeDerivedMirrors(HydroShadowPrivateRegion& region) const { (void)region; }
};

// Chronology points of one accepted T01 interval at which the runtime records
// the PRODUCTION state of a sampled region (read-only snapshots).
enum class HydroShadowChronologyPoint : std::uint8_t {
  kPreHalfBegin = 0,           // before S1 (state U0)
  kPreHalfAfterTransport = 1,  // after S1 element diffusion / SGS shear
  kPreHalfEnd = 2,             // after S1 thermochemistry (hydro input)
  kPostHalfBegin = 3,          // after H (+ composition commit, reflux)
  kPostHalfAfterThermochemistry = 4,
  kPostHalfAfterTransport = 5, // before discrete events
  kPostHalfEnd = 6,            // after discrete events (endpoint W)
};
inline constexpr std::size_t k_hydro_shadow_chronology_point_count = 7U;

struct HydroShadowMetric {
  double weighted_abs_difference = 0.0;
  double weighted_reference_magnitude = 0.0;
  double linf_difference = 0.0;
};

// Deterministic per-region result. Aggregation over ranks sorts these by
// region key so published totals do not depend on rank order.
struct HydroShadowRegionResult {
  std::uint64_t region_key = 0U;
  int level = 0;
  HydroShadowComparisonClass comparison_class = HydroShadowComparisonClass::kTimeOnlyCtuPpm;
  std::uint64_t compared_cells = 0U;
  std::uint64_t compared_auxiliary_cells = 0U;
  std::uint64_t cone_excluded_cells = 0U;
  std::uint64_t compared_composition_cells = 0U;
  double compared_volume = 0.0;
  std::array<HydroShadowMetric, 7> components{};  // rho,mx,my,mz,E,rhoZ,rho e
  HydroShadowMetric composition{};                // all composition lanes (mass)
  std::array<double, 6> reference_delta{};
  std::array<double, 6> production_delta{};
};
// Fixed wire width of one HydroShadowRegionResult (bytes).
inline constexpr std::size_t k_hydro_shadow_region_result_wire_bytes = 8U * 64U;
void appendHydroShadowRegionResultWire(const HydroShadowRegionResult& result, std::vector<std::uint8_t>& wire);
[[nodiscard]] HydroShadowRegionResult readHydroShadowRegionResultWire(std::span<const std::uint8_t> bytes);

struct HydroShadowReport {
  std::uint64_t eligible = 0U;
  std::uint64_t selected = 0U;
  std::uint64_t globally_selected_regions = 0U;
  std::uint64_t attempted = 0U;
  std::uint64_t completed = 0U;
  std::uint64_t completed_closed_boundary_patches = 0U;
  std::uint64_t completed_source_active_patches = 0U;
  std::uint64_t completed_evolving_comoving_patches = 0U;
  std::uint64_t completed_endpoint_ghost_patches = 0U;
  std::uint64_t completed_partial_active_patches = 0U;
  std::uint64_t completed_composition_patches = 0U;
  std::uint64_t completed_chronology_regions = 0U;
  std::array<std::uint64_t, k_hydro_shadow_comparison_class_count> completed_by_class{};
  std::uint64_t compared_cells = 0U;
  std::uint64_t compared_auxiliary_cells = 0U;
  std::uint64_t compared_composition_cells = 0U;
  std::uint64_t cone_excluded_cells = 0U;
  std::uint64_t skipped_cadence = 0U;
  std::uint64_t skipped_selection = 0U;
  std::uint64_t skipped_geometry = 0U;
  std::uint64_t skipped_sources_or_time = 0U;  // legacy aggregate (kept for compatibility)
  std::uint64_t skipped_resource = 0U;
  std::uint64_t skipped_numerical = 0U;
  std::uint64_t skipped_unsupported_source = 0U;
  std::uint64_t skipped_external_passives = 0U;
  std::uint64_t skipped_temporal_ghost = 0U;
  std::uint64_t skipped_endpoint_unavailable = 0U;
  std::uint64_t skipped_endpoint_geometry_mismatch = 0U;
  std::uint64_t skipped_region_split_across_ranks = 0U;
  std::uint64_t skipped_composition_endpoint_mismatch = 0U;
  std::uint64_t discarded_rejected_attempt = 0U;
  // Chronology (full-physics Strang interval) reasons.
  std::uint64_t chronology_started = 0U;
  std::uint64_t chronology_skipped_recursive_level_sources = 0U;
  std::uint64_t chronology_skipped_unsplit_compatibility = 0U;
  std::uint64_t chronology_skipped_no_replay_operator = 0U;
  std::uint64_t chronology_skipped_replay_unavailable = 0U;
  std::uint64_t chronology_skipped_replay_incomplete = 0U;
  std::uint64_t chronology_skipped_region_changed = 0U;
  std::uint64_t chronology_skipped_partial_activity = 0U;
  std::uint64_t chronology_skipped_missing_point = 0U;
  std::uint64_t chronology_skipped_hydro_sample_missing = 0U;
  std::uint64_t chronology_skipped_resource = 0U;
  double chronology_replay_identity_linf = 0.0;  // private S1 vs production S1
  std::uint64_t composition_remote_donor_faces = 0U;
  std::uint64_t reference_substeps = 0U;
  std::uint64_t fallback_count = 0U;
  std::array<std::uint64_t, k_hydro_shadow_ghost_provenance_count> ghost_slots_by_provenance{};
  std::uint64_t ghost_slots_self_refreshed = 0U;
  std::uint64_t ghost_slots_endpoint_interpolated = 0U;
  std::uint64_t ghost_slots_step_start_only = 0U;
  std::uint64_t workspace_high_water_bytes = 0U; // largest admitted reservation
  std::uint64_t observed_capacity_high_water_bytes = 0U;
  double compared_volume = 0.0;
  std::array<HydroShadowMetric, 7> components{}; // rho,mx,my,mz,E,rhoZ,rho e
  std::array<double, 6> reference_delta{};
  std::array<double, 6> production_delta{};
};

// Complete description of one sampled patch advance (preferred API).
struct HydroShadowPatchCapture {
  HydroShadowCandidate identity{};
  // Comparison regions within this domain; empty = one region equal to
  // `identity` covering every active real cell.
  std::span<const HydroShadowComparisonRegion> regions{};
  const HydroConservedStateSoa* before = nullptr;   // real cells + step-start ghosts
  const HydroPatchGeometry* geometry = nullptr;
  HydroActiveSetView active{};
  HydroUpdateContext update{};
  HydroNumericsPolicy production_numerics{};
  HydroReconstructionScheme production_reconstruction = HydroReconstructionScheme::kPiecewiseLinear;
  HydroReconstructionPolicy production_policy{};
  HydroRiemannSolverKind production_riemann = HydroRiemannSolverKind::kHllc;
  std::span<const HydroSourceTerm* const> source_terms{};
  const HydroSourceContext* source_context = nullptr;
  // Per ghost slot; empty = hydroShadowDefaultGhostProvenance().
  std::span<const HydroShadowGhostProvenance> ghost_provenance{};
  // True when the caller will call provideEndpointGhosts() for this sample at
  // a legal rank-uniform point after the production advance.
  bool endpoint_ghosts_will_be_provided = false;
  // Private composition transport (replaces the production face sink).
  const core::GasCompositionSidecar* composition = nullptr;  // read-only production sidecar
  std::span<const std::uint32_t> composition_real_rows{};     // dense row per real geometry cell
  // Dense composition row of the LOCAL origin of each ghost slot;
  // core::kInvalidGasCellRow = remote/unknown (donor fractions unavailable).
  std::span<const std::uint32_t> composition_ghost_rows{};
  // Composition-dependent temperature closure (T = factor * mu/0.59 * u).
  double temperature_k_per_internal_energy_code = 0.0;
  // Physical cgs density per comoving code density at a = 1.
  double density_cgs_per_code = 0.0;
  bool comoving_density = false;
  // Exact private halfstep epochs.
  const core::LambdaCdmBackground* cosmology = nullptr;
  double time_si_per_code = 0.0;
  // Retain this sample after its hydro comparison for the full-physics
  // chronology (set by the runtime when a chronology record is open).
  bool retain_for_chronology = false;
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
    HydroShadowCandidate identity{};
    std::vector<HydroShadowComparisonRegion> regions;
    HydroConservedStateSoa initial;
    HydroPatchGeometry geometry;
    std::vector<std::size_t> active_cells;
    std::vector<std::size_t> active_faces;
    HydroUpdateContext update{};
    HydroNumericsPolicy numerics{};
    HydroReconstructionScheme production_reconstruction = HydroReconstructionScheme::kPiecewiseLinear;
    HydroReconstructionPolicy reconstruction_policy{};
    HydroRiemannSolverKind riemann_kind = HydroRiemannSolverKind::kHllc;
    std::vector<std::unique_ptr<HydroSourceTerm>> shadow_sources;
    std::vector<HydroCellSourceProperties> captured_properties;
    std::vector<std::uint8_t> effective_closure_allowed;
    const HydroThermodynamicClosure* closure = nullptr;
    double source_redshift = 0.0;
    std::vector<HydroShadowGhostProvenance> provenance;
    std::vector<HydroShadowGhostTimeAuthority> authority;
    std::vector<HydroConservedState> endpoint_ghosts;
    std::vector<double> endpoint_ghost_auxiliary;
    std::vector<std::uint8_t> endpoint_slot_valid;
    std::vector<HydroConservedState> production_after;
    std::vector<double> production_after_auxiliary;
    // Private composition (real rows) and ghost donor fractions.
    core::GasCompositionSidecar composition;
    std::vector<double> ghost_donor_fractions;  // slot-major x lanes, NaN = unavailable
    std::vector<std::uint64_t> gas_cell_ids;
    std::array<double, 5> epoch_scale_factor{};
    std::array<double, 5> epoch_hubble_code{};
    double temperature_k_per_internal_energy_code = 0.0;
    double density_cgs_per_code = 0.0;
    bool comoving_density = false;
    bool source_properties_captured = false;
    bool closed_physical_boundary = false;
    bool source_active = false;
    bool evolving_comoving = false;
    bool temporal_refined = false;
    bool needs_endpoint = false;
    bool endpoint_received = false;
    bool production_observed = false;
    bool partial_active = false;
    bool retain_for_chronology = false;
    bool hydro_compared = false;
    bool composition_checked = false;
    std::uint64_t admitted_bytes = 0U;
    bool valid = false;
  };

  explicit HydroShadowReference(HydroShadowOptions options, core::MemoryGovernor* governor = nullptr);
  // Opens a new sampling epoch (one accepted-step ATTEMPT). Drops every
  // pending sample/chronology of a previous attempt (counted as rejected
  // when unfinished) and clears the report.
  void beginEpoch(std::uint64_t step_index);
  [[nodiscard]] std::uint64_t epochStepIndex() const noexcept { return m_step_index; }
  [[nodiscard]] bool cadenceStep() const noexcept {
    return m_options.enabled && (m_step_index % m_options.cadence) == 0U;
  }
  // The coordinator must call this once at a rank-uniform boundary with
  // global candidates computed by a bounded collective. No MPI inside
  // conditional patch advances.
  void setGlobalCandidates(std::vector<HydroShadowCandidate> ranked);
  [[nodiscard]] std::span<const HydroShadowCandidate> globalCandidates() const noexcept {
    return m_global_candidates;
  }
  [[nodiscard]] bool isGloballySelected(std::uint64_t region_key, int level) const noexcept;
  void noteGlobalSelectionResourceRefusal() noexcept { ++m_report.skipped_resource; }
  void noteRegionSplitAcrossRanks(std::uint64_t count) noexcept {
    m_report.skipped_region_split_across_ranks += count;
  }

  // Preferred capture. Returns an invalid sample when not selected, not
  // eligible or not admitted (each counted by reason).
  [[nodiscard]] Sample capture(const HydroShadowPatchCapture& spec);
  // Legacy capture (standalone callers): one region, default provenance,
  // no endpoint ghosts, no composition.
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

  // Records the production post-advance real-cell state (pre-reflux) and
  // either runs the comparison immediately (all ghosts self-refreshed or a
  // single-step reference) or keeps the sample pending its endpoint ghosts.
  void observeProduction(Sample&& sample, const HydroConservedStateSoa& production_after);
  // Legacy one-shot API: observeProduction() followed by finalizePendingEndpoints().
  void compare(
      Sample&& sample,
      const HydroConservedStateSoa& production_after,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active,
      const HydroUpdateContext& update,
      const HydroNumericsPolicy& production_numerics,
      HydroReconstructionScheme production_reconstruction);

  struct EndpointRequest {
    HydroShadowCandidate identity{};
    double time_end_code = 0.0;
    std::size_t ghost_slot_count = 0U;
  };
  [[nodiscard]] std::vector<EndpointRequest> pendingEndpointRequests() const;
  [[nodiscard]] std::size_t pendingEndpointCount() const noexcept;
  // Endpoint ghost values formed by the production ghost filler at the end
  // of the sampled interval. slot_valid[slot] == 0 marks a slot without
  // endpoint authority (it becomes kStepStartOnly). Size mismatches are a
  // counted geometry mismatch, never an exception into production.
  void provideEndpointGhosts(
      std::uint64_t region_key, double time_end_code,
      std::span<const HydroConservedState> ghost_states,
      std::span<const double> ghost_auxiliary,
      std::span<const std::uint8_t> slot_valid) noexcept;
  // Runs every pending hydro comparison. Samples without endpoint values use
  // kStepStartOnly for imported slots (cone exclusion). Never throws.
  void finalizePendingEndpoints() noexcept;

  // ---- Full-physics chronology (kCoupledFullPhysicsStrang) ----
  // Opens a chronology record for a locally owned selected domain. Rows are
  // dense SimulationState rows in domain geometry order; the record charges
  // its snapshot storage to the governor before allocation. Returns false
  // (counted) on refusal.
  // level_clock: a recursive-AMR level window S_L H_L S_L (no discrete
  // events inside; composition transport is committed by production only at
  // the coarse end, so the private window defers it identically).
  [[nodiscard]] bool openChronology(const HydroShadowCandidate& domain,
                                    std::span<const std::uint64_t> gas_cell_ids,
                                    std::size_t composition_lanes, bool level_clock = false);
  [[nodiscard]] bool chronologyOpen(std::uint64_t domain_key) const noexcept;
  [[nodiscard]] std::vector<HydroShadowCandidate> openChronologyDomains() const;
  [[nodiscard]] std::span<const std::uint64_t> chronologyGasCellIds(std::uint64_t domain_key) const noexcept;
  // Stores a read-only production snapshot at `point`. The snapshot ids must
  // equal the record ids (otherwise the record is dropped, counted).
  void recordChronologyPoint(std::uint64_t domain_key, HydroShadowChronologyPoint point,
                             const HydroShadowPrivateRegion& snapshot) noexcept;
  void dropChronology(std::uint64_t domain_key, std::uint64_t HydroShadowReport::*reason) noexcept;
  // Hydro-only composition check (any stepping mode): retained compared
  // samples with private composition, their gas-cell ids (domain order) and
  // the comparison against the committed post-hydro production snapshot.
  [[nodiscard]] std::vector<HydroShadowCandidate> compositionCheckDomains() const;
  [[nodiscard]] std::span<const std::uint64_t> sampleGasCellIds(std::uint64_t domain_key) const noexcept;
  void compareCompositionAfterCommit(std::uint64_t domain_key, const HydroShadowPrivateRegion& snapshot,
                                     double endpoint_time_code) noexcept;
  void noteChronologyUnsupported(std::uint64_t HydroShadowReport::*reason, std::uint64_t count) noexcept;
  // Runs every complete chronology with the production smooth-operator
  // replay (nullptr = counted unavailability). Never throws.
  // level_filter: finalize only the level-clock records of that level (at
  // the end of its first post half). Without a filter (global interval
  // close) level-clock records still open are dropped as missing points.
  void finalizeChronologies(const HydroShadowSmoothSourceReplay* replay,
                            std::optional<int> level_filter = std::nullopt) noexcept;

  // Every pending sample/chronology is dropped (rejected attempt), counted.
  void discardPending() noexcept;

  [[nodiscard]] const HydroShadowReport& report() const noexcept { return m_report; }
  [[nodiscard]] std::span<const HydroShadowRegionResult> regionResults() const noexcept {
    return m_region_results;
  }
  [[nodiscard]] const HydroShadowOptions& options() const noexcept { return m_options; }

 private:
  struct ChronologyRecord {
    core::MemoryReservation lease;
    HydroShadowCandidate domain{};
    std::vector<std::uint64_t> gas_cell_ids;
    std::array<std::optional<HydroShadowPrivateRegion>, k_hydro_shadow_chronology_point_count> points;
    std::size_t composition_lanes = 0U;
    bool level_clock = false;
  };

  void runHydroComparison(Sample& sample) noexcept;
  void runChronology(ChronologyRecord& record, Sample& sample,
                     const HydroShadowSmoothSourceReplay& replay) noexcept;
  void accumulateRegionResult(const HydroShadowRegionResult& result);

  HydroShadowOptions m_options;
  core::MemoryGovernor* m_governor = nullptr;
  std::uint64_t m_step_index = 0U;
  std::size_t m_storage_cells_charged = 0U;
  std::uint64_t m_workspace_bytes_charged = 0U;
  HydroShadowReport m_report{};
  std::vector<HydroShadowCandidate> m_global_candidates;
  std::vector<HydroShadowCandidate> m_selected_once;
  std::vector<Sample> m_pending;
  std::vector<ChronologyRecord> m_chronologies;
  std::vector<HydroShadowRegionResult> m_region_results;
};

}  // namespace cosmosim::hydro
