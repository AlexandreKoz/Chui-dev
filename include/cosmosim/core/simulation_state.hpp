#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/particle_species.hpp"
#include "cosmosim/core/soa_storage.hpp"

namespace cosmosim::core {

class HierarchicalTimeBinScheduler;

enum class ParticleMetadataRepresentation : std::uint8_t {
  kMaterializedGeneric = 0,
  kHomogeneousDmo = 1,
};


struct ParticleSoa {
  // Authoritative owner for persistent gravity-hot particle truth in SimulationState.
  // Ownership lanes: position_*_comoving, velocity_*_peculiar, mass_code.
  // Deliberately no persistent acceleration lane: force accelerations are transient
  // kernel/workspace outputs and must not become a second persistent truth vector.
  AlignedVector<double> position_x_comoving;
  AlignedVector<double> position_y_comoving;
  AlignedVector<double> position_z_comoving;
  AlignedVector<double> velocity_x_peculiar;
  AlignedVector<double> velocity_y_peculiar;
  AlignedVector<double> velocity_z_peculiar;
  AlignedVector<double> mass_code;
  // Derived mirror of scheduler bin_index for diagnostics/I/O compatibility;
  // exact restart uses scheduler state.
  AlignedVector<std::uint8_t> time_bin;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct ParticleSidecar {
  // Authoritative persistent metadata owner for particle IDs/species/rank and
  // ownership-adjacent particle metadata outside gravity-hot lanes.
  AlignedVector<std::uint64_t> particle_id;
  // Space-filling-curve key used for locality-preserving reorder/grouping.
  AlignedVector<std::uint64_t> sfc_key;
  // Authoritative species-tag ownership lane.
  AlignedVector<std::uint32_t> species_tag;
  AlignedVector<std::uint32_t> particle_flags;
  // Authoritative owning-rank lane for distributed ownership truth.
  AlignedVector<std::uint32_t> owning_rank;
  // Authoritative particle-position epoch used to predict inactive TreePM sources
  // without mutating persistent particle state during local active-bin refreshes.
  AlignedVector<double> last_drift_time_code;
  AlignedVector<double> last_drift_scale_factor;
  // Optional per-particle gravity softening value sidecar (comoving code units).
  // has_gravity_softening_override is the authoritative override mask. A populated
  // value vector without a populated mask is a materialized diagnostic/default mirror,
  // not independent per-particle override truth.
  AlignedVector<double> gravity_softening_comoving;
  AlignedVector<std::uint8_t> has_gravity_softening_override;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
  [[nodiscard]] bool hasGravitySofteningOverride(std::size_t particle_index) const;
  [[nodiscard]] double gravitySofteningOverride(std::size_t particle_index) const;
  void setGravitySofteningOverride(std::size_t particle_index, double epsilon_comoving);
  void clearGravitySofteningOverride(std::size_t particle_index);
};

struct CellSoa {
  // Gravity-facing gas-cell skeleton kept separate from hydro thermodynamics.
  AlignedVector<double> center_x_comoving;
  AlignedVector<double> center_y_comoving;
  AlignedVector<double> center_z_comoving;
  AlignedVector<double> mass_code;
  // Derived mirror of scheduler bin_index for diagnostics/I/O compatibility;
  // not persistent authority.
  AlignedVector<std::uint8_t> time_bin;
  AlignedVector<std::uint32_t> patch_index;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct GasCellSidecar {
  // Persistent gas-cell identity and thermodynamic state excluded from gravity-hot paths.
  // Reconstruction gradients are transient hydro/workspace scratch, not restart truth.
  // gas_cell_id is the stable nonzero identity mirror. parent_particle_id is legacy
  // lineage metadata: a value of zero mirrors a parentless GasCellIdentityRecord.
  // Cell-local velocity lanes are persistent hydro restart truth for parentless or
  // split cells; particle velocities are compatibility mirrors only when a unique
  // local parent particle is present.
  AlignedVector<std::uint64_t> gas_cell_id;
  AlignedVector<std::uint64_t> parent_particle_id;
  AlignedVector<double> velocity_x_peculiar;
  AlignedVector<double> velocity_y_peculiar;
  AlignedVector<double> velocity_z_peculiar;
  AlignedVector<double> density_code;
  // Persistent compatibility/restart caches. Pressure, temperature, and sound
  // speed are not independent conserved authority; hydro/source/scheduler
  // synchronization refreshes them from the accepted thermodynamic state.
  // Removing these lanes requires an explicit restart/schema migration.
  AlignedVector<double> pressure_code;
  AlignedVector<double> internal_energy_code;
  // Conserved gas metal mass. Metallicity is derived as metal_mass_code / cell mass.
  AlignedVector<double> metal_mass_code;
  AlignedVector<double> temperature_code;
  AlignedVector<double> sound_speed_code;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct StarParticleSidecar {
  // Stellar-formation metadata decoupled from common particle skeleton.
  AlignedVector<std::uint32_t> particle_index;
  AlignedVector<double> formation_scale_factor;
  AlignedVector<double> birth_mass_code;
  AlignedVector<double> metallicity_mass_fraction;
  // Immutable decomposition-invariant birth identity.
  AlignedVector<std::uint64_t> birth_key;
  AlignedVector<std::uint64_t> parent_gas_cell_id;
  AlignedVector<std::uint64_t> birth_tick;
  AlignedVector<std::uint32_t> birth_ordinal;
  // Explicit stellar evolution bookkeeping lanes for auditable enrichment.
  AlignedVector<double> stellar_age_years_last;
  AlignedVector<double> stellar_returned_mass_cumulative_code;
  AlignedVector<double> stellar_returned_metals_cumulative_code;
  // Newly synthesized metals are tracked separately from total ejected metals.
  AlignedVector<double> stellar_newly_synthesized_metals_cumulative_code;
  AlignedVector<double> stellar_feedback_energy_cumulative_erg;
  // Restart- and migration-authoritative unresolved enrichment budgets.
  AlignedVector<double> enrichment_carry_mass_code;
  AlignedVector<double> enrichment_carry_metals_code;
  AlignedVector<double> enrichment_carry_feedback_energy_erg;
  AlignedVector<double> enrichment_carry_momentum_code;
  AlignedVector<double> stellar_deposited_mass_cumulative_code;
  AlignedVector<double> stellar_deposited_metals_cumulative_code;
  AlignedVector<double> stellar_deposited_feedback_energy_cumulative_erg;
  std::array<AlignedVector<double>, 3> stellar_returned_mass_channel_cumulative_code;
  std::array<AlignedVector<double>, 3> stellar_returned_metals_channel_cumulative_code;
  std::array<AlignedVector<double>, 3> stellar_feedback_energy_channel_cumulative_erg;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct BlackHoleParticleSidecar {
  // Black-hole subgrid metadata sidecar.
  AlignedVector<std::uint32_t> particle_index;
  AlignedVector<std::uint32_t> host_cell_index;
  AlignedVector<double> subgrid_mass_code;
  AlignedVector<double> accretion_rate_code;
  AlignedVector<double> feedback_energy_code;
  AlignedVector<double> eddington_ratio;
  AlignedVector<double> cumulative_accreted_mass_code;
  AlignedVector<double> cumulative_feedback_energy_code;
  AlignedVector<double> duty_cycle_active_time_code;
  AlignedVector<double> duty_cycle_total_time_code;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct TracerParticleSidecar {
  // Tracer attachment metadata sidecar.
  AlignedVector<std::uint32_t> particle_index;
  AlignedVector<std::uint64_t> parent_particle_id;
  AlignedVector<std::uint64_t> injection_step;
  AlignedVector<std::uint32_t> host_cell_index;
  AlignedVector<double> mass_fraction_of_host;
  AlignedVector<double> last_host_mass_code;
  AlignedVector<double> cumulative_exchanged_mass_code;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

inline constexpr std::uint32_t kInvalidGasCellRow = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::size_t kMaxLocalParticleCount = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::size_t kMaxLocalCellCount = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::size_t kMaxLocalPatchCount = std::numeric_limits<std::uint32_t>::max();


[[nodiscard]] inline std::size_t checkedLocalCount(
    std::size_t count,
    std::size_t max_count,
    std::string_view kind,
    std::string_view context) {
  if (count > max_count) {
    throw std::length_error(
        std::string(context) + ": local " + std::string(kind) +
        " count exceeds uint32 index capacity; repartition/redistribute the domain");
  }
  return count;
}

[[nodiscard]] inline std::size_t checkedLocalCountAdd(
    std::size_t lhs,
    std::size_t rhs,
    std::size_t max_count,
    std::string_view kind,
    std::string_view context) {
  const std::size_t total = checkedSizeAdd(lhs, rhs, context);
  return checkedLocalCount(total, max_count, kind, context);
}

[[nodiscard]] inline std::uint32_t checkedLocalRow(
    std::size_t row,
    std::size_t max_count,
    std::string_view kind,
    std::string_view context) {
  if (row >= max_count) {
    throw std::length_error(
        std::string(context) + ": local " + std::string(kind) +
        " row exceeds uint32 row space; repartition/redistribute the domain");
  }
  return checkedIntegralNarrow<std::uint32_t>(row, context);
}

[[nodiscard]] inline std::uint32_t checkedLocalParticleRow(
    std::size_t row,
    std::string_view context) {
  return checkedLocalRow(row, kMaxLocalParticleCount, "particle", context);
}

[[nodiscard]] inline std::uint32_t checkedLocalCellRow(
    std::size_t row,
    std::string_view context) {
  return checkedLocalRow(row, kMaxLocalCellCount, "gas-cell", context);
}

[[nodiscard]] inline std::uint32_t checkedLocalPatchRow(
    std::size_t row,
    std::string_view context) {
  return checkedLocalRow(row, kMaxLocalPatchCount, "AMR patch", context);
}

[[nodiscard]] inline std::uint32_t checkedLocalPatchCellCount(
    std::size_t count,
    std::string_view context) {
  (void)checkedLocalCount(count, kMaxLocalCellCount, "AMR patch cell", context);
  return checkedIntegralNarrow<std::uint32_t>(count, context);
}

struct GasCellIdentityRecord {
  // Design-seam record for future AMR/moving-mesh gas ownership.
  // gas_cell_id is the stable identity and must not be derived from a particle row index.
  // A zero gas_cell_id is reserved as an uninitialized sentinel and is rejected by GasCellIdentityMap.
  std::uint64_t gas_cell_id = 0;
  std::optional<std::uint64_t> parent_particle_id;
  std::uint64_t owning_patch_id = 0;
  std::uint32_t local_cell_row = 0;
};

class GasCellIdentityMap {
 public:
  // Authoritative production row-space for gas-cell identity. Dense local rows remain transient
  // indexing; gas_cell_id is the stable key across reorder, restart, split/merge, and migration.
  // Lookup indices are rebuilt atomically on
  // assign(), and generation() changes after every successful mutation so future production views can
  // reject stale local-row mappings.
  void assign(std::vector<GasCellIdentityRecord> records);
  // Restart-authoritative import path: restore validated records together with
  // the generation value persisted at write time, so stale view checks see the
  // same identity epoch after reload. Normal runtime mutations must keep using
  // assign()/clear(), which bump the generation.
  void assignWithGeneration(std::vector<GasCellIdentityRecord> records, std::uint64_t generation);
  void clear() noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::uint64_t generation() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
  [[nodiscard]] bool coversDenseLocalRows(std::size_t cell_count) const noexcept;
  void requireCoversDenseLocalRows(std::size_t cell_count, std::string_view caller) const;
  [[nodiscard]] std::span<const GasCellIdentityRecord> records() const noexcept;
  [[nodiscard]] const GasCellIdentityRecord* findByGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] const GasCellIdentityRecord* findByLocalRow(std::uint32_t local_cell_row) const noexcept;
  [[nodiscard]] std::optional<std::uint32_t> rowForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> gasCellIdForLocalRow(std::uint32_t local_cell_row) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> parentParticleIdForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> owningPatchIdForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::vector<std::uint32_t> rowsForParentParticleId(std::uint64_t parent_particle_id) const;
  [[nodiscard]] std::vector<std::uint32_t> rowsForPatch(std::uint64_t owning_patch_id) const;
  // Exact owned bytes for the canonical record vector plus a conservative
  // bucket/node model for the two derived unordered lookup tables. The latter
  // is intentionally conservative so MemoryGovernor baseline reconciliation
  // does not omit identity-index residency hidden behind std::unordered_map.
  [[nodiscard]] std::uint64_t recordsCurrentSizeBytes() const;
  [[nodiscard]] std::uint64_t recordsOwnedCapacityBytes() const;
  [[nodiscard]] std::uint64_t lookupOwnedBytesConservative() const;

 private:
  bool rebuildLookupTables();

  std::vector<GasCellIdentityRecord> m_records;
  std::unordered_map<std::uint64_t, std::size_t> m_index_by_gas_cell_id;
  std::unordered_map<std::uint32_t, std::size_t> m_index_by_local_row;
  std::uint64_t m_generation = 0;
};

// Build a dense new-row -> old-row remapping by stable gas_cell_id. New cells that are absent from
// old_map are marked with kInvalidGasCellRow so callers must deliberately initialize their hydro state.
[[nodiscard]] std::vector<std::uint32_t> buildGasCellNewToOldRowMap(
    const GasCellIdentityMap& old_map,
    const GasCellIdentityMap& new_map);

struct PatchSoa {
  // AMR patch descriptors and contiguous cell ranges [first_cell, first_cell + cell_count).
  AlignedVector<std::uint64_t> patch_id;
  AlignedVector<std::int32_t> level;
  AlignedVector<std::uint32_t> first_cell;
  AlignedVector<std::uint32_t> cell_count;
  // Restart-authoritative patch geometry. Zero extents/dimensions mark a legacy
  // non-AMR descriptor; production AMR hydro requires these lanes to be populated.
  AlignedVector<std::uint64_t> parent_patch_id;
  AlignedVector<std::uint64_t> morton_key;
  AlignedVector<double> origin_x_comoving;
  AlignedVector<double> origin_y_comoving;
  AlignedVector<double> origin_z_comoving;
  AlignedVector<double> extent_x_comoving;
  AlignedVector<double> extent_y_comoving;
  AlignedVector<double> extent_z_comoving;
  AlignedVector<std::uint16_t> cell_dim_x;
  AlignedVector<std::uint16_t> cell_dim_y;
  AlignedVector<std::uint16_t> cell_dim_z;
  // Authoritative AMR patch owner for distributed rebalancing. Patch payloads remain
  // contiguous in CellSoa; this lane is the ownership contract used by load-balance
  // planning and restart, not transient exchange scratch.
  AlignedVector<std::uint32_t> owning_rank;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct SpeciesContainer {
  // Explicit species counts; used as an auditable ownership ledger.
  std::array<std::uint64_t, k_particle_species_count> count_by_species{};

  [[nodiscard]] std::uint64_t totalCount() const noexcept;
  [[nodiscard]] bool isConsistentWith(const ParticleSidecar& sidecar) const noexcept;
};

struct ParticleSpeciesIndex {
  enum class Representation : std::uint8_t { kMaterialized = 0, kHomogeneousDmoIdentity = 1 };
  // Explicit species-local to global particle index mapping for generic states.
  std::array<AlignedVector<std::uint32_t>, k_particle_species_count> global_index_by_species;
  AlignedVector<std::uint32_t> local_index_by_global;
  Representation representation = Representation::kMaterialized;
  std::uint32_t identity_count = 0;

  void rebuild(const ParticleSidecar& sidecar);
  void setHomogeneousDmoIdentity(std::uint32_t particle_count);
  [[nodiscard]] std::size_t count(ParticleSpecies species) const;
  [[nodiscard]] std::span<const std::uint32_t> globalIndices(ParticleSpecies species) const;
  [[nodiscard]] std::uint32_t localIndex(std::uint32_t global_index) const;
  [[nodiscard]] std::uint32_t globalIndex(ParticleSpecies species, std::uint32_t local_index) const;
};

struct StateMetadata {
  // Schema/provenance fields that must remain stable across restart/snapshot workflows.
  std::uint32_t schema_version = 3;
  std::string run_name = "cosmosim_run";
  std::uint64_t normalized_config_hash = 0;
  std::string normalized_config_hash_hex;
  std::uint64_t step_index = 0;
  double scale_factor = 1.0;
  std::string snapshot_stem = "snapshot";
  std::string restart_stem = "restart";

  [[nodiscard]] std::string serialize() const;
  [[nodiscard]] static StateMetadata deserialize(std::string_view text);
};

enum class ModuleSidecarRequirementKind : std::uint8_t {
  kSparse = 0,
  kSpeciesMask = 1,
  kGasDensityAtLeast = 2,
  kBlackHoleAccretionAtLeast = 3,
  kParticleFlagMask = 4,
};

struct ModuleSidecarRequirement {
  // Structured coverage contract for particle-indexed module sidecars.
  // kSparse means optional rows.  kSpeciesMask reproduces the legacy mask.
  // Other predicates are evaluated against the authoritative particle/gas/BH
  // state during pack/commit so modules can require rows only for physically
  // active subsets instead of every particle of a species.
  ModuleSidecarRequirementKind kind = ModuleSidecarRequirementKind::kSparse;
  std::uint32_t species_mask = 0;
  std::uint32_t particle_flags_mask = 0;
  double threshold_code = 0.0;
};

struct ModuleSidecarBlock {
  // Opaque module payload with an independent schema version. Non-particle-indexed
  // blocks are run/module persistent state and are preserved as whole blocks.
  std::string module_name;
  std::uint32_t schema_version = 1;
  std::vector<std::byte> payload;

  // Optional particle-indexed row layout for module-owned persistent particle
  // sidecars. When enabled, each row is keyed by stable particle_id and the
  // corresponding byte slice must migrate with that particle.
  bool particle_indexed = false;
  std::uint32_t row_stride_bytes = 0;
  std::vector<std::uint64_t> particle_id_by_row;

  // Legacy unconditional species mask retained for compatibility. Structured
  // requirement is ORed with this mask when deciding whether a row is required.
  std::uint32_t required_species_mask = 0;
  ModuleSidecarRequirement requirement{};

  [[nodiscard]] bool isParticleIndexed() const noexcept;
  [[nodiscard]] bool requiresSpecies(ParticleSpecies species) const noexcept;
  [[nodiscard]] std::size_t rowCount() const noexcept;
  [[nodiscard]] std::span<const std::byte> rowPayload(std::size_t row) const;
};

struct ModuleParticleSidecarPayload {
  std::string module_name;
  std::uint32_t schema_version = 1;
  std::uint32_t row_stride_bytes = 0;
  std::uint32_t required_species_mask = 0;
  ModuleSidecarRequirement requirement{};
  std::vector<std::byte> payload;
};

class ModuleSidecarRegistry {
 public:
  // Insert or replace sidecar payload for a module.
  void upsert(ModuleSidecarBlock block);
  void clear() noexcept;
  [[nodiscard]] const ModuleSidecarBlock* find(std::string_view module_name) const;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] std::vector<const ModuleSidecarBlock*> blocksSortedByName() const;

 private:
  std::unordered_map<std::string, ModuleSidecarBlock> m_sidecars;
};

struct ParticleTransferPacket {
  // Species-specific transfer packet for MPI or host-device staging.
  ParticleSpecies species = ParticleSpecies::kDarkMatter;
  AlignedVector<std::uint64_t> particle_id;
  AlignedVector<double> position_x_comoving;
  AlignedVector<double> position_y_comoving;
  AlignedVector<double> position_z_comoving;
  AlignedVector<double> velocity_x_peculiar;
  AlignedVector<double> velocity_y_peculiar;
  AlignedVector<double> velocity_z_peculiar;
  AlignedVector<double> mass_code;
  // Derived timestep mirror copied from scheduler authority for transfer locality only.
  // Receiver-side scheduler remap/rebuild must validate or overwrite this lane
  // before exact continuation.
  AlignedVector<std::uint8_t> time_bin;
  AlignedVector<std::uint32_t> owning_rank;
  AlignedVector<double> last_drift_time_code;
  AlignedVector<double> last_drift_scale_factor;
};



struct PendingFluxRegisterRecord {
  std::uint64_t register_key = 0;
  std::uint64_t coarse_patch_id = 0;
  std::uint64_t coarse_gas_cell_id = 0;
  std::size_t coarse_cell_index = 0;
  std::uint8_t level = 0;
  std::uint8_t axis = 0;
  std::uint8_t orientation = 0;
  double expected_area_comov = 0.0;
  double coarse_area_accumulated_comov = 0.0;
  double fine_area_accumulated_comov = 0.0;
  double interval_start_code = 0.0;
  double interval_end_code = 0.0;
  double coarse_dt_code = 0.0;
  std::uint32_t expected_fine_substeps = 1;
  std::uint32_t completed_fine_substeps = 0;
  std::uint64_t fine_substep_coverage_mask = 0;
  std::uint32_t coarse_face_count = 0;
  std::uint32_t fine_face_count = 0;
  std::uint64_t gas_cell_identity_generation = 0;
  std::uint64_t patch_geometry_generation = 0;
  double coarse_mass_flux_integral_code = 0.0;
  double coarse_momentum_x_flux_integral_code = 0.0;
  double coarse_momentum_y_flux_integral_code = 0.0;
  double coarse_momentum_z_flux_integral_code = 0.0;
  double coarse_total_energy_flux_integral_code = 0.0;
  double coarse_metal_mass_flux_integral_code = 0.0;
  double fine_mass_flux_integral_code = 0.0;
  double fine_momentum_x_flux_integral_code = 0.0;
  double fine_momentum_y_flux_integral_code = 0.0;
  double fine_momentum_z_flux_integral_code = 0.0;
  double fine_total_energy_flux_integral_code = 0.0;
  double fine_metal_mass_flux_integral_code = 0.0;

  [[nodiscard]] bool hasCoarseContribution() const noexcept { return coarse_face_count > 0U; }
  [[nodiscard]] bool hasCompleteFineSubstepCoverage() const noexcept {
    return expected_fine_substeps > 0U && completed_fine_substeps >= expected_fine_substeps;
  }
  [[nodiscard]] bool isComplete() const noexcept {
    return hasCoarseContribution() && hasCompleteFineSubstepCoverage();
  }
};

struct AmrTemporalBoundaryHistoryCellRecord {
  // Stable, row-order-independent source identity for one coarse patch cell.
  std::uint64_t gas_cell_id = 0;
  std::size_t patch_local_cell = 0;
  double start_mass_density_comoving = 0.0;
  double start_momentum_density_x_comoving = 0.0;
  double start_momentum_density_y_comoving = 0.0;
  double start_momentum_density_z_comoving = 0.0;
  double start_total_energy_density_comoving = 0.0;
  double start_metal_mass_density_comoving = 0.0;
  double end_mass_density_comoving = 0.0;
  double end_momentum_density_x_comoving = 0.0;
  double end_momentum_density_y_comoving = 0.0;
  double end_momentum_density_z_comoving = 0.0;
  double end_total_energy_density_comoving = 0.0;
  double end_metal_mass_density_comoving = 0.0;
};

struct AmrTemporalBoundaryHistoryRecord {
  // Persistent coarse-to-fine temporal boundary interval.  Geometry and cell
  // identities are stored explicitly so a restart or dense-row reorder cannot
  // turn a history entry into a positional lookup.
  std::uint64_t patch_id = 0;
  std::uint8_t patch_level = 0;
  std::uint64_t patch_geometry_fingerprint = 0;
  std::uint64_t gas_cell_identity_generation = 0;
  double interval_start_code = 0.0;
  double interval_end_code = 0.0;
  bool end_state_valid = false;
  std::vector<AmrTemporalBoundaryHistoryCellRecord> cells;
};

class AmrTemporalBoundaryHistoryStore {
 public:
  [[nodiscard]] bool empty() const noexcept { return m_records.empty(); }
  [[nodiscard]] std::size_t size() const noexcept { return m_records.size(); }
  [[nodiscard]] std::uint64_t currentSizeBytes() const;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const;
  [[nodiscard]] std::span<const AmrTemporalBoundaryHistoryRecord> records() const noexcept { return m_records; }
  [[nodiscard]] std::span<AmrTemporalBoundaryHistoryRecord> mutableRecords() noexcept { return m_records; }
  void clear() noexcept { m_records.clear(); }
  void assign(std::vector<AmrTemporalBoundaryHistoryRecord> records) { m_records = std::move(records); }
  [[nodiscard]] const AmrTemporalBoundaryHistoryRecord* findByPatchId(std::uint64_t patch_id) const noexcept {
    for (const auto& record : m_records) {
      if (record.patch_id == patch_id) {
        return &record;
      }
    }
    return nullptr;
  }
  [[nodiscard]] AmrTemporalBoundaryHistoryRecord* findByPatchId(std::uint64_t patch_id) noexcept {
    for (auto& record : m_records) {
      if (record.patch_id == patch_id) {
        return &record;
      }
    }
    return nullptr;
  }
  AmrTemporalBoundaryHistoryRecord& upsertByPatchId(const AmrTemporalBoundaryHistoryRecord& record) {
    if (auto* existing = findByPatchId(record.patch_id); existing != nullptr) {
      *existing = record;
      return *existing;
    }
    m_records.push_back(record);
    return m_records.back();
  }

 private:
  std::vector<AmrTemporalBoundaryHistoryRecord> m_records;
};

class PendingFluxRegisterStore {
 public:
  [[nodiscard]] bool empty() const noexcept { return m_records.empty(); }
  [[nodiscard]] std::size_t size() const noexcept { return m_records.size(); }
  [[nodiscard]] std::uint64_t currentSizeBytes() const;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const;
  [[nodiscard]] std::span<const PendingFluxRegisterRecord> records() const noexcept { return m_records; }
  [[nodiscard]] std::span<PendingFluxRegisterRecord> mutableRecords() noexcept { return m_records; }
  void clear() noexcept { m_records.clear(); }
  void assign(std::vector<PendingFluxRegisterRecord> records) { m_records = std::move(records); }
  [[nodiscard]] const PendingFluxRegisterRecord* findByRegisterKey(std::uint64_t register_key) const noexcept {
    for (const auto& record : m_records) {
      if (record.register_key == register_key) {
        return &record;
      }
    }
    return nullptr;
  }
  [[nodiscard]] PendingFluxRegisterRecord* findByRegisterKey(std::uint64_t register_key) noexcept {
    for (auto& record : m_records) {
      if (record.register_key == register_key) {
        return &record;
      }
    }
    return nullptr;
  }
  PendingFluxRegisterRecord& upsertByRegisterKey(const PendingFluxRegisterRecord& record) {
    if (auto* existing = findByRegisterKey(record.register_key); existing != nullptr) {
      return *existing;
    }
    m_records.push_back(record);
    return m_records.back();
  }
  void eraseCompletedByKey(std::span<const std::uint64_t> register_keys) {
    m_records.erase(
        std::remove_if(
            m_records.begin(),
            m_records.end(),
            [register_keys](const PendingFluxRegisterRecord& record) {
              return std::find(register_keys.begin(), register_keys.end(), record.register_key) != register_keys.end();
            }),
        m_records.end());
  }
  // Erase records whose register_key has been zeroed (successfully applied).
  // This avoids materializing a separate keys vector during pending reflux
  // application.
  void eraseCompleted() noexcept {
    m_records.erase(
        std::remove_if(
            m_records.begin(),
            m_records.end(),
            [](const PendingFluxRegisterRecord& record) {
              return record.register_key == 0U;
            }),
        m_records.end());
  }

 private:
  std::vector<PendingFluxRegisterRecord> m_records;
};

struct GasCellMigrationFields {
  std::uint64_t gas_cell_id = 0;
  std::uint8_t has_parent_particle = 0;
  std::uint64_t parent_particle_id = 0;
  std::uint64_t owning_patch_id = 0;
  std::uint32_t destination_local_cell_row = kInvalidGasCellRow;
  std::uint64_t gas_cell_identity_generation = 0;
  std::uint64_t ghost_hydro_epoch = 0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
  double cell_mass_code = 0.0;
  std::uint8_t cell_time_bin = 0;
  std::uint32_t patch_index = 0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double density_code = 0.0;
  double pressure_code = 0.0;
  double internal_energy_code = 0.0;
  double metal_mass_code = 0.0;
  double temperature_code = 0.0;
  double sound_speed_code = 0.0;
};

struct GasCellMigrationRecord {
  // Standalone gas-cell ownership payload keyed by gas_cell_id. This is separate
  // from ParticleMigrationRecord so parentless and split/merged cells can migrate
  // without inventing a particle parent.
  std::uint32_t owning_rank = 0;
  GasCellMigrationFields fields{};
};

struct GasCellStaleGhostRecord {
  std::uint64_t gas_cell_id = 0;
  std::uint32_t local_cell_row = kInvalidGasCellRow;
  std::uint64_t gas_cell_identity_generation = 0;
  std::uint64_t ghost_hydro_epoch = 0;
};

struct GasCellMigrationCommit {
  int world_rank = 0;
  std::uint64_t expected_gas_cell_identity_generation = 0;
  std::uint64_t expected_ghost_hydro_epoch = 0;
  std::vector<std::uint32_t> outbound_local_cell_rows;
  std::vector<GasCellMigrationRecord> inbound_records;
  std::vector<GasCellStaleGhostRecord> stale_local_ghost_records;
};

struct AmrPatchMigrationFields {
  std::uint64_t patch_id = 0;
  std::uint64_t parent_patch_id = 0;
  std::uint64_t morton_key = 0;
  std::int32_t level = 0;
  std::uint32_t owning_rank = 0;
  std::uint32_t cell_count = 0;
  double origin_x_comoving = 0.0;
  double origin_y_comoving = 0.0;
  double origin_z_comoving = 0.0;
  double extent_x_comoving = 0.0;
  double extent_y_comoving = 0.0;
  double extent_z_comoving = 0.0;
  std::uint16_t cell_dim_x = 0;
  std::uint16_t cell_dim_y = 0;
  std::uint16_t cell_dim_z = 0;
};

struct GasCellSchedulerMigrationRecord {
  std::uint64_t gas_cell_id = 0;
  std::uint8_t bin_index = 0;
  std::uint64_t next_activation_tick = 0;
  std::uint8_t pending_bin_index = 0xFF;
};

struct AmrPatchMigrationRecord {
  // Atomic patch ownership payload. Patch descriptors, every authoritative
  // gas-cell row, scheduler identity, and AMR synchronization state owned by
  // the patch migrate together. Pending reflux records are keyed by the
  // coarse patch/gas-cell identity; temporal histories are keyed by patch_id.
  // Dense-row indices/generations are remapped on destination commit.
  AmrPatchMigrationFields patch{};
  std::vector<GasCellMigrationRecord> gas_cell_records;
  std::vector<GasCellSchedulerMigrationRecord> gas_cell_scheduler_records;
  std::vector<PendingFluxRegisterRecord> pending_flux_register_records;
  std::vector<AmrTemporalBoundaryHistoryRecord> temporal_boundary_history_records;
};

struct AmrPatchMigrationCommit {
  int world_rank = 0;
  std::uint64_t expected_gas_cell_identity_generation = 0;
  std::uint64_t expected_ghost_hydro_epoch = 0;
  std::vector<std::uint32_t> outbound_local_patch_indices;
  std::vector<std::uint32_t> stale_local_ghost_patch_indices;
  std::vector<AmrPatchMigrationRecord> inbound_records;
  std::vector<GasCellStaleGhostRecord> stale_local_ghost_records;
};

struct StarParticleMigrationFields {
  double formation_scale_factor = 0.0;
  double birth_mass_code = 0.0;
  double metallicity_mass_fraction = 0.0;
  std::uint64_t birth_key = 0;
  std::uint64_t parent_gas_cell_id = 0;
  std::uint64_t birth_tick = 0;
  std::uint32_t birth_ordinal = 0;
  double stellar_age_years_last = 0.0;
  double stellar_returned_mass_cumulative_code = 0.0;
  double stellar_returned_metals_cumulative_code = 0.0;
  double stellar_newly_synthesized_metals_cumulative_code = 0.0;
  double stellar_feedback_energy_cumulative_erg = 0.0;
  double enrichment_carry_mass_code = 0.0;
  double enrichment_carry_metals_code = 0.0;
  double enrichment_carry_feedback_energy_erg = 0.0;
  double enrichment_carry_momentum_code = 0.0;
  double stellar_deposited_mass_cumulative_code = 0.0;
  double stellar_deposited_metals_cumulative_code = 0.0;
  double stellar_deposited_feedback_energy_cumulative_erg = 0.0;
  std::array<double, 3> stellar_returned_mass_channel_cumulative_code{};
  std::array<double, 3> stellar_returned_metals_channel_cumulative_code{};
  std::array<double, 3> stellar_feedback_energy_channel_cumulative_erg{};
};

struct BlackHoleParticleMigrationFields {
  std::uint32_t host_cell_index = kInvalidGasCellRow;
  double subgrid_mass_code = 0.0;
  double accretion_rate_code = 0.0;
  double feedback_energy_code = 0.0;
  double eddington_ratio = 0.0;
  double cumulative_accreted_mass_code = 0.0;
  double cumulative_feedback_energy_code = 0.0;
  double duty_cycle_active_time_code = 0.0;
  double duty_cycle_total_time_code = 0.0;
};

struct TracerParticleMigrationFields {
  std::uint64_t parent_particle_id = 0;
  std::uint64_t injection_step = 0;
  std::uint32_t host_cell_index = kInvalidGasCellRow;
  double mass_fraction_of_host = 0.0;
  double last_host_mass_code = 0.0;
  double cumulative_exchanged_mass_code = 0.0;
};

struct SchedulerMigrationFields {
  // Scheduler-authoritative lanes for an element keyed by stable particle_id.
  // time_bin in ParticleSoa remains a derived mirror; these fields are the
  // migration payload needed to rebuild exact hierarchical-timestep state.
  std::uint8_t bin_index = 0;
  std::uint64_t next_activation_tick = 0;
  std::uint8_t pending_bin_index = 0xFF;
};

struct ParticleMigrationRecord {
  std::uint64_t particle_id = 0;
  std::uint64_t sfc_key = 0;
  std::uint32_t species_tag = 0;
  std::uint32_t particle_flags = 0;
  std::uint32_t owning_rank = 0;
  double position_x_comoving = 0.0;
  double position_y_comoving = 0.0;
  double position_z_comoving = 0.0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double mass_code = 0.0;
  // Derived timestep mirror copied with the migrating particle for diagnostics
  // and legacy transfer consumers. Exact continuation must use scheduler_fields
  // when has_scheduler_fields is true.
  std::uint8_t time_bin = 0;
  bool has_scheduler_fields = false;
  SchedulerMigrationFields scheduler_fields{};
  double last_drift_time_code = 0.0;
  double last_drift_scale_factor = 1.0;
  bool has_gravity_softening_value = false;
  bool has_gravity_softening_override = false;
  double gravity_softening_comoving = 0.0;
  bool has_gas_cell_fields = false;
  GasCellMigrationFields gas_cell_fields{};
  bool has_star_fields = false;
  StarParticleMigrationFields star_fields{};
  bool has_black_hole_fields = false;
  BlackHoleParticleMigrationFields black_hole_fields{};
  bool has_tracer_fields = false;
  TracerParticleMigrationFields tracer_fields{};
  std::vector<ModuleParticleSidecarPayload> module_sidecar_payloads;
};

struct ParticleMigrationCommit {
  int world_rank = 0;
  std::vector<std::uint32_t> outbound_local_indices;
  std::vector<ParticleMigrationRecord> inbound_records;
  std::vector<std::uint32_t> stale_local_ghost_indices;
  // Distributed gas/AMR migration paths move gas cells explicitly by
  // gas_cell_id/patch_id. In that mode particle migration must not infer gas
  // authority from optional parent_particle_id compatibility mirrors.
  bool preserve_gas_cell_state = false;
};

struct ParticleIdValidationResult {
  bool unique = true;
  bool nonzero = true;
};

struct ParticleIdentityValidationCertificate {
  std::uint64_t identity_generation = std::numeric_limits<std::uint64_t>::max();
  std::size_t particle_count = 0U;
  bool local_uniqueness_proven = false;
  bool nonzero_ids_proven = false;

  [[nodiscard]] bool validFor(
      std::uint64_t generation, std::size_t count, bool require_nonzero) const noexcept {
    return identity_generation == generation && particle_count == count &&
        local_uniqueness_proven && (!require_nonzero || nonzero_ids_proven);
  }
};

// Reusable scratch for exact ownership/ID validation. Particle identity uses a
// fixed-capacity radix bucket rather than population-sized scratch. Cold
// sidecar marker lanes are allocated only when the corresponding sidecar is
// present, and cell_owner exists only when AMR patch ownership must be checked.
struct OwnershipValidationWorkspace {
  static constexpr std::uint64_t k_id_validation_workspace_limit_bytes =
      32ULL * 1024ULL * 1024ULL;

  std::vector<std::uint64_t> id_bucket;
  std::vector<std::uint8_t> star_rows;
  std::vector<std::uint8_t> bh_rows;
  std::vector<std::uint8_t> tracer_rows;
  std::vector<std::uint32_t> cell_owner;

  void resize(
      std::size_t particle_count, std::size_t cell_count,
      bool need_star_rows = false, bool need_bh_rows = false,
      bool need_tracer_rows = false, bool need_cell_owner = false);
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const;
  [[nodiscard]] static std::uint64_t requiredBytes(
      std::size_t particle_count, std::size_t cell_count,
      bool need_star_rows = false, bool need_bh_rows = false,
      bool need_tracer_rows = false, bool need_cell_owner = false,
      bool include_id_validation_workspace = true);
};

// Exact, deterministic particle-ID validation with a population-independent
// CHUI-owned workspace. Callers that already hold a valid identity certificate
// should reuse that evidence instead of invoking this fallback.
[[nodiscard]] ParticleIdValidationResult validateParticleIdsExact(
    std::span<const std::uint64_t> particle_ids,
    OwnershipValidationWorkspace& scratch);

class SimulationState {
 public:
  // Single ownership root for persistent run state.
  ParticleSoa particles;
  ParticleSidecar particle_sidecar;
  CellSoa cells;
  GasCellSidecar gas_cells;
  // Authoritative in-memory gas-cell identity map. During the H2 legacy transition,
  // gas_cells.gas_cell_id and gas_cells.parent_particle_id are compatibility mirrors
  // synchronized from this map for particle-bound import/restart paths.
  GasCellIdentityMap gas_cell_identity;
  PatchSoa patches;
  // Restart-authoritative AMR hydro deferred synchronization state.
  PendingFluxRegisterStore pending_flux_registers;
  // Restart-authoritative coarse temporal boundary intervals retained while
  // local AMR subcycling has fine work outstanding.
  AmrTemporalBoundaryHistoryStore amr_temporal_boundary_history;
  SpeciesContainer species;
  ParticleSpeciesIndex particle_species_index;
  StarParticleSidecar star_particles;
  BlackHoleParticleSidecar black_holes;
  TracerParticleSidecar tracers;
  StateMetadata metadata;
  ModuleSidecarRegistry sidecars;

  void resizeParticles(std::size_t count);
  void resizeCells(std::size_t count);
  void resizePatches(std::size_t count);

  [[nodiscard]] bool validateOwnershipInvariants() const;
  [[nodiscard]] bool validateOwnershipInvariants(OwnershipValidationWorkspace& scratch) const;
  [[nodiscard]] bool validateUniqueParticleIds() const;
  [[nodiscard]] bool validateUniqueParticleIds(OwnershipValidationWorkspace& scratch) const;
  [[nodiscard]] bool validatePersistentParticleIds() const;
  [[nodiscard]] bool validatePersistentParticleIds(OwnershipValidationWorkspace& scratch) const;
  [[nodiscard]] const ParticleIdentityValidationCertificate&
  particleIdentityValidationCertificate() const noexcept {
    return m_particle_identity_validation_certificate;
  }
  [[nodiscard]] bool hasCertifiedUniqueParticleIds() const noexcept {
    return m_particle_identity_validation_certificate.validFor(
        m_particle_identity_generation, particles.size(), false);
  }
  // Callers that rewrite particle_sidecar.particle_id in place, without a
  // resize/reorder/migration commit, must call this so cached identity
  // validation certificates are not reused for the mutated lane.
  void notifyParticleIdentityMutation() noexcept { bumpParticleIdentityGeneration(); }
  [[nodiscard]] bool hasCertifiedPersistentParticleIds() const noexcept {
    return m_particle_identity_validation_certificate.validFor(
        m_particle_identity_generation, particles.size(), true);
  }
  void rebuildSpeciesIndex();
  [[nodiscard]] ParticleMetadataRepresentation particleMetadataRepresentation() const noexcept { return m_particle_metadata_representation; }
  [[nodiscard]] bool hasHomogeneousDmoMetadata() const noexcept { return m_particle_metadata_representation == ParticleMetadataRepresentation::kHomogeneousDmo; }
  [[nodiscard]] bool compactHomogeneousDmoMetadata(std::uint32_t local_rank);
  void materializeParticleMetadata();
  [[nodiscard]] std::uint32_t homogeneousDmoParticleFlags() const;
  [[nodiscard]] double homogeneousDmoLastDriftTimeCode() const;
  [[nodiscard]] double homogeneousDmoLastDriftScaleFactor() const;
  [[nodiscard]] bool homogeneousDmoSfcKeyIsUniform() const;
  [[nodiscard]] std::uint64_t homogeneousDmoUniformSfcKey() const;
  // Strong-commit seam for the compact homogeneous-DMO migration transaction.
  // Candidate lanes are fully prepared/validated before this no-allocation swap.
  void commitCompactHomogeneousDmoCandidate(
      ParticleSoa&& candidate_particles,
      ParticleSidecar&& candidate_sidecar,
      std::uint32_t local_rank,
      std::uint32_t uniform_particle_flags,
      double common_last_drift_time_code,
      double common_last_drift_scale_factor,
      bool persistent_particle_ids_certified = false);
  [[nodiscard]] std::uint32_t particleSpeciesTag(std::size_t particle_index) const;
  [[nodiscard]] std::uint32_t particleFlags(std::size_t particle_index) const;
  [[nodiscard]] std::uint32_t particleOwningRank(std::size_t particle_index) const;
  [[nodiscard]] std::uint64_t particleSfcKey(std::size_t particle_index) const;
  [[nodiscard]] double particleLastDriftTimeCode(std::size_t particle_index) const;
  [[nodiscard]] double particleLastDriftScaleFactor(std::size_t particle_index) const;
  [[nodiscard]] std::uint8_t particleTimeBin(std::size_t particle_index) const;
  void updateAllParticleDriftEpoch(double time_code, double scale_factor);
  void refreshGasCellIdentityFromParticleOrder();
  [[nodiscard]] bool gasCellIdentityMatchesParticleOrder() const;
  // Legacy/bootstrap imports only: these construct the authoritative map from
  // old particle-bound or sidecar representations. Normal production mutation
  // must flow in the opposite direction through replaceGasCellIdentityRecords().
  void refreshGasCellIdentityMapFromParticleBoundState();
  void refreshGasCellIdentityMapFromSidecarLanes();

  // Authoritative production mutation path. Records are validated and installed
  // atomically in GasCellIdentityMap, then compatibility lanes are rebuilt from
  // the map. Dense rows remain storage positions only.
  void replaceGasCellIdentityRecords(std::vector<GasCellIdentityRecord> records);
  // Restart-only variant preserving the serialized identity generation.
  void restoreGasCellIdentityRecords(
      std::vector<GasCellIdentityRecord> records,
      std::uint64_t generation);
  // Rebuild compatibility mirror lanes from GasCellIdentityMap. This never
  // imports identity from mirrors and is safe at production boundaries.
  void synchronizeGasCellIdentityCompatibilityMirrors();

  [[nodiscard]] bool gasCellIdentityMapMatchesParticleBoundState() const;
  [[nodiscard]] std::uint64_t gasCellIdentityGeneration() const noexcept;
  void requireGasCellIdentityMapCoversDenseRows(std::string_view caller) const;
  void requireGasCellIdentityMapFresh(std::uint64_t expected_generation, std::string_view caller) const;
  [[nodiscard]] std::optional<std::uint32_t> rowForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> gasCellIdForLocalRow(std::uint32_t local_cell_row) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> parentParticleIdForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> owningPatchIdForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] bool gasCellIdentityMapMatchesSidecarLanes() const;
  void requireParticleBoundGasCellContract(std::string_view caller) const;
  void legacyRequireParticleBoundGasCellContract(std::string_view caller) const;
  [[nodiscard]] std::optional<std::uint64_t> parentParticleIdForGasCellRow(std::uint32_t cell_index) const;
  [[nodiscard]] std::uint32_t gasCellRowForParticleId(std::uint64_t particle_id) const;
  [[nodiscard]] std::uint32_t gasParticleIndexForCellRow(std::uint32_t cell_index) const;

  [[nodiscard]] ParticleTransferPacket packSpeciesTransferPacket(ParticleSpecies species_tag) const;
  [[nodiscard]] std::vector<ParticleMigrationRecord> packParticleMigrationRecords(
      std::span<const std::uint32_t> local_indices,
      const HierarchicalTimeBinScheduler& scheduler) const;
  void commitParticleMigration(const ParticleMigrationCommit& commit);
  [[nodiscard]] std::vector<GasCellMigrationRecord> packGasCellMigrationRecords(
      std::span<const std::uint32_t> local_cell_rows,
      std::uint64_t ghost_hydro_epoch = 0) const;
  void commitGasCellMigration(const GasCellMigrationCommit& commit);
  [[nodiscard]] std::vector<AmrPatchMigrationRecord> packAmrPatchMigrationRecords(
      std::span<const std::uint32_t> local_patch_indices,
      std::uint64_t ghost_hydro_epoch = 0) const;
  void commitAmrPatchMigration(const AmrPatchMigrationCommit& commit);
  [[nodiscard]] std::uint64_t particleIndexGeneration() const noexcept;
  [[nodiscard]] std::uint64_t cellIndexGeneration() const noexcept;
  // Monotonic identity of authoritative gravity-relevant source state. This
  // is owned by SimulationState, not synthesized from PM refresh sequencing
  // or rank-local dense-row/index generations. Physical source mutations must
  // bump this explicitly at a communicator-consistent workflow boundary.
  [[nodiscard]] std::uint64_t gravitySourceGeneration() const noexcept;
  void bumpGravitySourceGeneration() noexcept;
  // Dense index generations are rank-local cache/view invalidation epochs.
  // Reorder, migration, or empty-rank topology can advance them differently
  // across ranks without changing the communicator-global gravity source epoch.
  void bumpParticleIndexGeneration() noexcept;
  void bumpCellIndexGeneration() noexcept;
  // Rank-local particle mass epoch for mass-derived caches (e.g. the
  // hierarchical encounter criterion's local maximum mass). Advanced by every
  // dense-row change (bumpParticleIndexGeneration: reorder, resize,
  // compaction, migration, population replacement) and by every in-place
  // particle mass write, which must call notifyParticleMassMutation(). Never
  // key a mass cache on particleIndexGeneration() alone: masses can change in
  // place (star formation, black holes, tracers, gas-parent aggregation).
  [[nodiscard]] std::uint64_t particleMassGeneration() const noexcept;
  void notifyParticleMassMutation() noexcept;
  // Gas thermodynamic cache epochs (rank-local, non-persistent; audit F4).
  // GasCellSidecar::pressure_code / sound_speed_code / temperature_code are
  // derived caches of (density, internal energy, composition, closure). A
  // stage that mutates the authoritative thermodynamic inputs outside the
  // hydro primitive store (e.g. feedback deposition) must call
  // notifyGasThermodynamicMutation(); the single refresh contract
  // hydro::refreshGasThermodynamics() calls markGasThermodynamicsRefreshed().
  // Timestep selection refuses to read stale caches.
  [[nodiscard]] std::uint64_t gasThermodynamicMutationEpoch() const noexcept;
  [[nodiscard]] std::uint64_t gasThermodynamicRefreshEpoch() const noexcept;
  void notifyGasThermodynamicMutation() noexcept;
  void markGasThermodynamicsRefreshed() noexcept;
  [[nodiscard]] bool gasThermodynamicCachesFresh() const noexcept;

 private:
  [[nodiscard]] bool validateOwnershipInvariantsImpl(
      OwnershipValidationWorkspace& scratch, bool use_bounded_id_scratch) const;
  [[nodiscard]] std::vector<ParticleMigrationRecord> packParticleMigrationRecordsCore(
      std::span<const std::uint32_t> local_indices) const;

  ParticleMetadataRepresentation m_particle_metadata_representation = ParticleMetadataRepresentation::kMaterializedGeneric;
  std::uint32_t m_uniform_species_tag = static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter);
  std::uint32_t m_uniform_particle_flags = 0U;
  std::uint32_t m_uniform_owning_rank = 0U;
  std::uint64_t m_uniform_sfc_key = 0U;
  bool m_sfc_key_is_uniform = false;
  double m_common_last_drift_time_code = 0.0;
  double m_common_last_drift_scale_factor = 1.0;
  void bumpParticleIdentityGeneration() noexcept;

  std::uint64_t m_particle_index_generation = 0;
  std::uint64_t m_particle_mass_generation = 1;
  std::uint64_t m_particle_identity_generation = 1;
  std::uint64_t m_cell_index_generation = 0;
  std::uint64_t m_gravity_source_generation = 1;
  std::uint64_t m_gas_thermodynamic_mutation_epoch = 0;
  std::uint64_t m_gas_thermodynamic_refresh_epoch = 0;
  mutable ParticleIdentityValidationCertificate m_particle_identity_validation_certificate{};
};

template <typename T>
inline constexpr bool k_is_canonical_particle_state_owner_v = false;

template <>
inline constexpr bool k_is_canonical_particle_state_owner_v<SimulationState> = true;

}  // namespace cosmosim::core

#include "cosmosim/core/simulation_state_runtime.hpp"
