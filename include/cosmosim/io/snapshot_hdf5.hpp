#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::core { class MemoryGovernor; }

namespace cosmosim::io {

// Science-snapshot semantics are explicit. CHUI-native files may contain
// project extensions; AREPO/GADGET exports use their documented cosmological
// storage scalings and never rely on a fictional merged dialect.
enum class SnapshotDialect {
  kAuto,
  kChuiNative,
  kArepoFormat3,
  kGadget4Hdf5,
};

enum class SnapshotMissingFieldPolicy {
  kReject,
  kMarkUnavailable,
  kFillDocumentedDefault,
};

struct SnapshotFieldAliases {
  std::string_view canonical_name;
  std::vector<std::string_view> read_aliases;
};

struct ScienceSnapshotSchemaMap {
  std::string_view schema_name = "chui_science_snapshot_v6";
  std::uint32_t schema_version = 6;

  std::string_view header_group = "/Header";
  std::string_view config_group = "/Config";
  std::string_view parameters_group = "/Parameters";
  std::string_view provenance_group = "/Provenance";
  std::string_view units_group = "/Units";
  std::string_view config_normalized_attribute = "normalized";

  std::array<std::string_view, 6> part_type_group = {
      "/PartType0", "/PartType1", "/PartType2", "/PartType3", "/PartType4", "/PartType5"};

  SnapshotFieldAliases coordinates{"Coordinates", {"Coordinates", "Position", "POS"}};
  SnapshotFieldAliases velocities{"Velocities", {"Velocities", "Velocity", "VEL"}};
  SnapshotFieldAliases masses{"Masses", {"Masses", "Mass"}};
  SnapshotFieldAliases particle_ids{"ParticleIDs", {"ParticleIDs", "ParticleID", "ID"}};
};

// Source compatibility aliases. New code should use ScienceSnapshot* names.
using GadgetArepoFieldAliases = SnapshotFieldAliases;
using GadgetArepoSchemaMap = ScienceSnapshotSchemaMap;

[[nodiscard]] const ScienceSnapshotSchemaMap& scienceSnapshotSchemaMap();
[[nodiscard]] const GadgetArepoSchemaMap& gadgetArepoSchemaMap();

struct SnapshotReadBudget {
  std::uint64_t max_particles = 1ULL << 34U;
  std::uint64_t max_gas_cells = 1ULL << 34U;
  std::uint64_t max_sidecar_rows = 1ULL << 34U;
  std::uint64_t max_materialized_bytes = 1ULL << 40U;
  std::uint64_t max_dataset_bytes = 1ULL << 38U;
  std::uint64_t max_attribute_bytes = 1ULL << 24U;
  std::uint64_t max_validation_id_bytes = 1ULL << 30U;
  std::uint32_t max_members = 4096U;
};

struct SnapshotSetMemberInfo {
  std::uint32_t member_index = 0;
  std::uint32_t num_files_per_snapshot = 1;
  std::array<std::uint64_t, 6> global_part_count{};
  bool has_global_part_count = false;
  std::string generation_id;
};

struct SnapshotIoPolicy {
  bool enable_compression = false;
  int compression_level = 4;
  std::size_t chunk_particle_count = 1u << 15;
  bool write_particle_type_alias_groups = false;
  bool durable_publication = false;
  SnapshotDialect dialect = SnapshotDialect::kAuto;
  bool write_optional_pressure = true;
  // Collective Parallel-HDF5 callers manage the .partial -> final publication
  // after distributed readback and failure agreement.
  bool caller_managed_publication = false;
};

struct SnapshotWritePayload {
  // Optional owner authority. Standalone callers retain the legacy ungoverned
  // API; production workflow supplies its existing process governor.
  core::MemoryGovernor* memory_governor = nullptr;
  const core::SimulationState* state = nullptr;
  const core::SimulationConfig* config = nullptr;
  std::string normalized_config_text;
  core::ProvenanceRecord provenance;
  std::string git_sha = "unknown";
  SnapshotSetMemberInfo set_member;
  // Parallel-HDF5 one-file topology. Global datasets use global_part_count;
  // each rank writes only its non-overlapping file_row_offset range.
  bool collective_single_file = false;
  std::array<std::uint64_t, 6> collective_file_row_offset{};
  bool collective_write_particle_softening = false;
  bool collective_write_particle_softening_override = false;
  // chui_science_snapshot_v7 /MHD/Diagnostics (magnetized states only). A
  // collective single file requires rank-agreed global values supplied by the
  // caller (mhd_diagnostics_global=true); sharded members write their own
  // rank-local values with DiagnosticsScope="file_member".
  bool mhd_diagnostics_global = false;
  double mhd_divergence_max_norm = 0.0;
  double mhd_magnetic_energy_comoving = 0.0;
  double mhd_scale_factor = 1.0;
};

// CHUI-MHD amendment R6: magnetized science output schema (v6 + derived
// PartType0/MagneticField [B_comov, normalized code units], PartType0/
// MagneticPressure and /MHD/Diagnostics). Never exact-restart authority.
inline constexpr std::string_view k_science_snapshot_mhd_schema_name = "chui_science_snapshot_v7";
inline constexpr std::uint32_t k_science_snapshot_mhd_schema_version = 7U;

// Derived cell-centred B_comov (arithmetic mean of the two bordering face
// fields per axis) of every gas row of `state`, evaluated from the
// authoritative MagneticFaceSoa through the leaf patch geometry; written in
// bounded chunks by the snapshot writer (no full-population mirror).
// Returns the rank-local comoving magnetic energy sum B_c^2/(2a) V.
[[nodiscard]] double mhdLocalMagneticEnergyComoving(const core::SimulationState& state,
                                                    double scale_factor);

struct SnapshotReadOptions {
  bool require_ids = true;
  bool allow_generated_particle_ids = false;
  bool require_velocities = true;
  bool allow_mass_table_fallback = true;
  bool require_complete_chui_set = true;
  bool verify_snapshot_set_member_hashes = true;
  // Generic offline validation proves global ID uniqueness by materializing
  // the logical ID set. Distributed production output may disable only this
  // redundant O(N_global) proof after exact rank-partition readback has already
  // established file IDs == authoritative globally unique runtime IDs.
  bool validate_global_id_uniqueness = true;
  SnapshotDialect dialect = SnapshotDialect::kAuto;
  SnapshotMissingFieldPolicy missing_field_policy =
      SnapshotMissingFieldPolicy::kReject;
  SnapshotReadBudget budget{};
  std::array<std::optional<core::ParticleSpecies>, 6> part_type_species_map{};
};

struct SnapshotIoReport {
  std::string schema_name = "unknown";
  std::uint32_t schema_version = 0;
  SnapshotDialect dialect = SnapshotDialect::kAuto;
  bool storage_policy_known = false;
  bool storage_policy_mixed = false;
  bool compression_enabled = false;
  int compression_level = 0;
  std::size_t chunk_particle_count = 0;
  std::uint64_t peak_staging_bytes = 0;
  std::uint64_t logical_bytes_written = 0;
  std::uint64_t chunk_write_count = 0;
  std::vector<std::string> present_aliases;
  std::vector<std::string> defaulted_fields;
  std::vector<std::string> unavailable_fields;
  std::string file_kind = "unknown";
  bool restart_compatible = false;
  bool analysis_ready = false;
  bool evolution_ready = false;
  std::vector<std::string> evolution_readiness_reasons;
  std::uint64_t materialized_state_bytes = 0;
  double header_time = 1.0;
  double header_redshift = 0.0;
  double header_box_size_x = 0.0;
  double header_box_size_y = 0.0;
  double header_box_size_z = 0.0;
  double header_omega_matter = 0.0;
  double header_omega_lambda = 0.0;
  double header_omega_baryon = 0.0;
  double header_hubble_param = 0.0;
  std::uint32_t member_index = 0;
  std::uint32_t num_files_per_snapshot = 1;
  std::array<std::uint64_t, 6> local_part_count{};
  std::array<std::uint64_t, 6> global_part_count{};
  std::string generation_id;
};

struct SnapshotReadResult {
  core::SimulationState state;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  SnapshotIoReport report;

  void requireEvolutionReady() const;
};

struct SnapshotSetInspection {
  SnapshotDialect dialect = SnapshotDialect::kAuto;
  std::string file_kind;
  std::string schema_name;
  std::uint32_t schema_version = 0;
  std::uint32_t num_files_per_snapshot = 0;
  std::array<std::uint64_t, 6> global_part_count{};
  double scale_factor = 1.0;
  double redshift = 0.0;
  double box_size_x = 0.0;
  double box_size_y = 0.0;
  double box_size_z = 0.0;
  double omega_matter = 0.0;
  double omega_lambda = 0.0;
  double omega_baryon = 0.0;
  double hubble_param = 0.0;
  std::string unit_length;
  std::string unit_mass;
  std::string unit_velocity;
  std::string coordinate_frame;
  std::string velocity_storage_convention;
  std::string config_hash_hex;
  std::string naming_rules_version;
  std::string file_naming_rules_version;
  std::string generation_id;
  std::vector<std::filesystem::path> member_paths;
  std::uint64_t total_member_file_bytes = 0;
  bool complete = false;
};

// Header-only identity of one snapshot member file. Used by the structural
// write-time verification policy; reads no particle datasets.
struct SnapshotMemberInspection {
  std::string file_kind;
  std::string schema_name;
  std::uint32_t schema_version = 0;
  SnapshotDialect dialect = SnapshotDialect::kAuto;
  std::uint32_t num_files_per_snapshot = 0;
  std::uint32_t member_index = 0;
  bool has_member_index = false;
  std::array<std::uint64_t, 6> local_part_count{};
  std::array<std::uint64_t, 6> global_part_count{};
  double scale_factor = 1.0;
  double redshift = 0.0;
  double box_size_x = 0.0;
  double box_size_y = 0.0;
  double box_size_z = 0.0;
  double omega_matter = 0.0;
  double omega_lambda = 0.0;
  double omega_baryon = 0.0;
  double hubble_param = 0.0;
  std::string config_hash_hex;
  std::string generation_id;
  std::uint64_t file_bytes = 0;
};

struct SnapshotValidationReport {
  bool valid = false;
  SnapshotSetInspection inspection;
  std::uint64_t datasets_checked = 0;
  std::uint64_t scalar_values_checked = 0;
  std::uint64_t particle_ids_checked = 0;
  std::uint64_t validation_peak_id_bytes = 0;
  std::vector<std::string> diagnostics;

  void requireValid() const;
};

void writeScienceSnapshotHdf5(
    const std::filesystem::path& output_path,
    const SnapshotWritePayload& payload,
    const SnapshotIoPolicy& policy = {});

[[nodiscard]] SnapshotReadResult readCosmoSimScienceSnapshotHdf5(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    const SnapshotReadOptions& options = {});

[[nodiscard]] SnapshotReadResult importExternalSnapshotHdf5(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    const SnapshotReadOptions& options);

[[nodiscard]] SnapshotSetInspection inspectSnapshotSet(
    const std::filesystem::path& input_path,
    const SnapshotReadOptions& options = {});

[[nodiscard]] SnapshotMemberInspection inspectSnapshotMemberHdf5(
    const std::filesystem::path& member_path,
    const SnapshotReadOptions& options = {});

[[nodiscard]] SnapshotValidationReport validateSnapshotSetHdf5(
    const std::filesystem::path& input_path,
    const SnapshotReadOptions& options = {});

enum class SnapshotCompletionIntegrityMode : std::uint8_t {
  kMemberSha256 = 0,
  kDistributedScienceReadback = 1,
};

void writeSnapshotSetCompletionMarker(
    const std::filesystem::path& snapshot_directory,
    std::string_view generation_id,
    std::uint32_t num_files_per_snapshot,
    const std::array<std::uint64_t, 6>& global_part_count,
    bool durable_publication = false,
    SnapshotCompletionIntegrityMode integrity_mode =
        SnapshotCompletionIntegrityMode::kMemberSha256);

// Compatibility entry points retained for existing callers. They now use the
// explicit CHUI science-snapshot implementation rather than a merged dialect.

// Verify only the rank-owned hyperslabs of a collective single-file science
// snapshot against the authoritative local state. The check is bounded and
// does not materialize global state on rank 0.
void verifySingleFileScienceSnapshotPartitionHdf5(
    const std::filesystem::path& input_path,
    const SnapshotWritePayload& payload,
    const SnapshotIoPolicy& policy = {});

// Publish a fully written and distributed-readback-verified science file.
void publishSingleFileScienceSnapshot(
    const std::filesystem::path& temporary_path,
    const std::filesystem::path& final_path,
    const SnapshotSetMemberInfo& member,
    bool durable_publication = false);

void writeGadgetArepoSnapshotHdf5(
    const std::filesystem::path& output_path,
    const SnapshotWritePayload& payload,
    const SnapshotIoPolicy& policy = {});

[[nodiscard]] SnapshotReadResult readGadgetArepoSnapshotHdf5(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    const SnapshotReadOptions& options = {});

}  // namespace cosmosim::io
