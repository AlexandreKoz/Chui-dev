#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/retained_capacity_transaction.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

namespace cosmosim::core {
namespace {

[[nodiscard]] bool migrationPayloadHasExplicitSparseRow(
    const GasCellMigrationFields& fields,
    const GasCompositionSidecar& composition) {
  if (!composition.sparseStorageEnabled()) return false;
  if (fields.composition_payload.size() != composition.packedRowWidth()) {
    throw std::runtime_error(
        "gas-cell migration composition payload width cannot preserve sparse chemistry history");
  }
  return fields.composition_payload.back() != 0.0;
}

void prepareAdmittedSparseCompositionCandidate(
    GasCompositionSidecar& composition,
    std::size_t active_rows,
    std::string_view owner) {
  if (!composition.sparseStorageEnabled()) return;
  RetainedCapacityTransaction capacity(
      [&composition] { return composition.ownedCapacityBytes(); });
  composition.planSparseRowCapacity(capacity, active_rows);
  // Migration/rebuild callers have already admitted the complete candidate
  // live set. This transaction enforces exact capacities without charging a
  // second overlapping reservation.
  capacity.execute(nullptr, MemoryClass::kCanonicalPersistent, owner);
}

[[nodiscard]] std::size_t findRequiredSidecarRow(
    std::span<const std::uint32_t> particle_indices,
    std::uint32_t particle_index,
    const char* sidecar_name) {
  for (std::size_t row = 0; row < particle_indices.size(); ++row) {
    if (particle_indices[row] == particle_index) {
      return row;
    }
  }
  throw std::runtime_error(std::string("missing required sidecar row for particle migration: ") + sidecar_name);
}

[[nodiscard]] std::optional<std::size_t> findModuleSidecarRowByParticleId(
    const ModuleSidecarBlock& block,
    std::uint64_t particle_id) {
  for (std::size_t row = 0; row < block.particle_id_by_row.size(); ++row) {
    if (block.particle_id_by_row[row] == particle_id) {
      return row;
    }
  }
  return std::nullopt;
}

void appendModuleRowPayload(ModuleSidecarBlock* block, std::uint64_t particle_id, std::span<const std::byte> payload) {
  if (!block->particle_indexed || block->row_stride_bytes == 0U) {
    throw std::invalid_argument("appendModuleRowPayload: destination module sidecar is not particle-indexed");
  }
  if (payload.size() != static_cast<std::size_t>(block->row_stride_bytes)) {
    throw std::invalid_argument("appendModuleRowPayload: payload row size does not match module sidecar stride");
  }
  block->particle_id_by_row.push_back(particle_id);
  block->payload.insert(block->payload.end(), payload.begin(), payload.end());
}

[[nodiscard]] ParticleIdValidationResult validateParticleIdsBoundedExactImpl(
    std::span<const std::uint64_t> particle_ids,
    OwnershipValidationWorkspace& scratch) {
  ParticleIdValidationResult result;
  for (const std::uint64_t id : particle_ids) {
    result.nonzero = result.nonzero && id != 0U;
  }
  if (particle_ids.size() < 2U) {
    return result;
  }

  const std::size_t bucket_capacity = static_cast<std::size_t>(
      OwnershipValidationWorkspace::k_id_validation_workspace_limit_bytes /
      sizeof(std::uint64_t));
  if (bucket_capacity < 2U) {
    throw std::logic_error("particle-ID validation workspace is too small");
  }
  scratch.id_bucket.clear();
  if (scratch.id_bucket.capacity() > bucket_capacity) {
    std::vector<std::uint64_t>{}.swap(scratch.id_bucket);
  }
  scratch.id_bucket.reserve(std::min(bucket_capacity, particle_ids.size()));

  struct Prefix {
    std::uint8_t bits = 0U;
    std::uint64_t value = 0U;
  };
  std::array<Prefix, 65U> pending{};
  std::size_t pending_size = 1U;
  pending[0] = Prefix{};

  const auto matches = [](std::uint64_t id, const Prefix& prefix) noexcept {
    if (prefix.bits == 0U) return true;
    if (prefix.bits == 64U) return id == prefix.value;
    const std::uint64_t mask = (1ULL << prefix.bits) - 1ULL;
    return (id & mask) == prefix.value;
  };

  while (pending_size != 0U) {
    const Prefix prefix = pending[--pending_size];
    std::size_t count = 0U;
    for (const std::uint64_t id : particle_ids) {
      if (matches(id, prefix)) {
        ++count;
        if (count > bucket_capacity) break;
      }
    }
    if (count > bucket_capacity) {
      if (prefix.bits == 64U) {
        result.unique = false;
        scratch.id_bucket.clear();
        return result;
      }
      if (pending_size + 2U > pending.size()) {
        throw std::logic_error("particle-ID validation prefix stack overflow");
      }
      const std::uint8_t child_bits = static_cast<std::uint8_t>(prefix.bits + 1U);
      pending[pending_size++] = Prefix{child_bits, prefix.value | (1ULL << prefix.bits)};
      pending[pending_size++] = Prefix{child_bits, prefix.value};
      continue;
    }

    scratch.id_bucket.clear();
    for (const std::uint64_t id : particle_ids) {
      if (matches(id, prefix)) scratch.id_bucket.push_back(id);
    }
    std::sort(scratch.id_bucket.begin(), scratch.id_bucket.end());
    if (std::adjacent_find(scratch.id_bucket.begin(), scratch.id_bucket.end()) !=
        scratch.id_bucket.end()) {
      result.unique = false;
      scratch.id_bucket.clear();
      return result;
    }
  }
  scratch.id_bucket.clear();
  return result;
}

[[nodiscard]] ParticleSpecies speciesFromTagOrThrow(std::uint32_t species_tag, const char* caller) {
  if (!isValidParticleSpeciesTag(species_tag)) {
    throw std::invalid_argument(std::string(caller) + ": invalid species tag");
  }
  return static_cast<ParticleSpecies>(species_tag);
}

[[nodiscard]] bool moduleRequiresSpeciesMask(
    const ModuleSidecarBlock& block,
    std::uint32_t species_tag,
    const char* caller) {
  return block.requiresSpecies(speciesFromTagOrThrow(species_tag, caller));
}

[[nodiscard]] bool moduleRequirementMatchesRecord(
    const ModuleSidecarRequirement& requirement,
    const ParticleMigrationRecord& record) {
  switch (requirement.kind) {
    case ModuleSidecarRequirementKind::kSparse:
      return false;
    case ModuleSidecarRequirementKind::kSpeciesMask: {
      if (!isValidParticleSpeciesTag(record.species_tag)) {
        return false;
      }
      const std::uint32_t bit = record.species_tag;
      return (requirement.species_mask & (1U << bit)) != 0U;
    }
    case ModuleSidecarRequirementKind::kGasDensityAtLeast:
      return record.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kGas) &&
          record.has_gas_cell_fields &&
          record.gas_cell_fields.density_code >= requirement.threshold_code;
    case ModuleSidecarRequirementKind::kBlackHoleAccretionAtLeast:
      return record.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole) &&
          record.has_black_hole_fields &&
          record.black_hole_fields.accretion_rate_code >= requirement.threshold_code;
    case ModuleSidecarRequirementKind::kParticleFlagMask:
      return (record.particle_flags & requirement.particle_flags_mask) == requirement.particle_flags_mask;
  }
  return false;
}

[[nodiscard]] bool moduleRequiresRecord(
    const ModuleSidecarBlock& block,
    const ParticleMigrationRecord& record,
    const char* caller) {
  return moduleRequiresSpeciesMask(block, record.species_tag, caller) ||
      moduleRequirementMatchesRecord(block.requirement, record);
}

void requireValidGasCellMigrationFields(const GasCellMigrationFields& fields, const char* caller) {
  if (fields.gas_cell_id == 0U) {
    throw std::invalid_argument(std::string(caller) + ": gas_cell_id must be nonzero");
  }
  if (fields.has_parent_particle > 1U) {
    throw std::invalid_argument(std::string(caller) + ": has_parent_particle must be 0 or 1");
  }
  if (fields.has_parent_particle == 0U && fields.parent_particle_id != 0U) {
    throw std::invalid_argument(std::string(caller) + ": parent_particle_id must be 0 for parentless cells");
  }
  if (fields.has_parent_particle != 0U && fields.parent_particle_id == 0U) {
    throw std::invalid_argument(std::string(caller) + ": parented cells require a nonzero parent_particle_id");
  }
}

void requireValidAmrPatchMigrationFields(const AmrPatchMigrationFields& fields, const char* caller) {
  if (fields.patch_id == 0U) {
    throw std::invalid_argument(std::string(caller) + ": patch_id must be nonzero");
  }
  if (fields.cell_count == 0U) {
    throw std::invalid_argument(std::string(caller) + ": patch cell_count must be nonzero");
  }
  const bool has_geometry = fields.extent_x_comoving > 0.0 || fields.extent_y_comoving > 0.0 ||
      fields.extent_z_comoving > 0.0 || fields.cell_dim_x != 0U || fields.cell_dim_y != 0U ||
      fields.cell_dim_z != 0U;
  if (has_geometry) {
    if (!(fields.extent_x_comoving > 0.0 && fields.extent_y_comoving > 0.0 &&
          fields.extent_z_comoving > 0.0) ||
        fields.cell_dim_x == 0U || fields.cell_dim_y == 0U || fields.cell_dim_z == 0U) {
      throw std::invalid_argument(std::string(caller) + ": patch geometry fields must be complete");
    }
    const std::size_t dims_product = static_cast<std::size_t>(fields.cell_dim_x) *
        static_cast<std::size_t>(fields.cell_dim_y) * static_cast<std::size_t>(fields.cell_dim_z);
    if (dims_product != fields.cell_count) {
      throw std::invalid_argument(std::string(caller) + ": patch cell dimensions must match cell_count");
    }
  }
}

[[nodiscard]] GasCellMigrationFields gasCellFieldsFromLocalRow(
    const SimulationState& state,
    std::uint32_t row,
    std::uint64_t ghost_hydro_epoch,
    const char* caller) {
  if (row >= state.cells.size()) {
    throw std::out_of_range(std::string(caller) + ": local gas-cell row is out of range");
  }
  // Callers validate the dense identity map once at the operation boundary.
  // Repeating that O(N) validation for every row makes migration packing/rebuild O(N^2).
  const auto* identity = state.gas_cell_identity.findByLocalRow(row);
  if (identity == nullptr) {
    throw std::runtime_error(std::string(caller) + ": gas-cell identity map is missing local row");
  }

  GasCellMigrationFields fields;
  fields.gas_cell_id = identity->gas_cell_id;
  fields.has_parent_particle = identity->parent_particle_id.has_value() ? 1U : 0U;
  fields.parent_particle_id = identity->parent_particle_id.value_or(0U);
  fields.owning_patch_id = identity->owning_patch_id;
  fields.destination_local_cell_row = row;
  fields.gas_cell_identity_generation = state.gasCellIdentityGeneration();
  fields.ghost_hydro_epoch = ghost_hydro_epoch;
  fields.center_x_comoving = state.cells.center_x_comoving[row];
  fields.center_y_comoving = state.cells.center_y_comoving[row];
  fields.center_z_comoving = state.cells.center_z_comoving[row];
  fields.cell_mass_code = state.cells.mass_code[row];
  fields.cell_time_bin = state.cells.time_bin[row];
  fields.patch_index = state.cells.patch_index[row];
  fields.velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[row];
  fields.velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[row];
  fields.velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[row];
  fields.density_code = state.gas_cells.density_code[row];
  fields.pressure_code = state.gas_cells.pressure_code[row];
  fields.internal_energy_code = state.gas_cells.internal_energy_code[row];
  fields.metal_mass_code = state.gas_cells.metal_mass_code[row];
  fields.temperature_code = state.gas_cells.temperature_code[row];
  fields.sound_speed_code = state.gas_cells.sound_speed_code[row];
  if (state.gas_composition.enabled()) {
    // CHUI Part 2: composition travels with the cell by stable identity.
    if (state.gas_composition.gas_cell_id[row] != fields.gas_cell_id) {
      throw std::runtime_error(std::string(caller) +
          ": composition sidecar row is not aligned with the gas-cell identity map");
    }
    fields.composition_payload = state.gas_composition.packRow(row);
  }
  requireValidGasCellMigrationFields(fields, caller);
  return fields;
}

void writeGasCellFieldsToRow(
    const GasCellMigrationFields& fields,
    std::uint32_t row,
    CellSoa& cells,
    GasCellSidecar& gas_cells,
    GasCompositionSidecar& composition) {
  requireValidGasCellMigrationFields(fields, "writeGasCellFieldsToRow");
  composition.unpackRow(row, fields.composition_payload, fields.gas_cell_id);
  cells.center_x_comoving[row] = fields.center_x_comoving;
  cells.center_y_comoving[row] = fields.center_y_comoving;
  cells.center_z_comoving[row] = fields.center_z_comoving;
  cells.mass_code[row] = fields.cell_mass_code;
  cells.time_bin[row] = fields.cell_time_bin;
  cells.patch_index[row] = fields.patch_index;
  gas_cells.gas_cell_id[row] = fields.gas_cell_id;
  gas_cells.parent_particle_id[row] = fields.has_parent_particle != 0U ? fields.parent_particle_id : 0U;
  gas_cells.velocity_x_peculiar[row] = fields.velocity_x_peculiar;
  gas_cells.velocity_y_peculiar[row] = fields.velocity_y_peculiar;
  gas_cells.velocity_z_peculiar[row] = fields.velocity_z_peculiar;
  gas_cells.density_code[row] = fields.density_code;
  gas_cells.pressure_code[row] = fields.pressure_code;
  gas_cells.internal_energy_code[row] = fields.internal_energy_code;
  gas_cells.metal_mass_code[row] = fields.metal_mass_code;
  gas_cells.temperature_code[row] = fields.temperature_code;
  gas_cells.sound_speed_code[row] = fields.sound_speed_code;
}

[[nodiscard]] AmrPatchMigrationFields patchFieldsFromLocalRow(
    const SimulationState& state,
    std::uint32_t patch_index,
    const char* caller) {
  if (patch_index >= state.patches.size()) {
    throw std::out_of_range(std::string(caller) + ": local patch index is out of range");
  }
  AmrPatchMigrationFields fields;
  fields.patch_id = state.patches.patch_id[patch_index];
  fields.parent_patch_id = state.patches.parent_patch_id[patch_index];
  fields.morton_key = state.patches.morton_key[patch_index];
  fields.level = state.patches.level[patch_index];
  fields.owning_rank = state.patches.owning_rank[patch_index];
  fields.cell_count = state.patches.cell_count[patch_index];
  fields.origin_x_comoving = state.patches.origin_x_comoving[patch_index];
  fields.origin_y_comoving = state.patches.origin_y_comoving[patch_index];
  fields.origin_z_comoving = state.patches.origin_z_comoving[patch_index];
  fields.extent_x_comoving = state.patches.extent_x_comoving[patch_index];
  fields.extent_y_comoving = state.patches.extent_y_comoving[patch_index];
  fields.extent_z_comoving = state.patches.extent_z_comoving[patch_index];
  fields.cell_dim_x = state.patches.cell_dim_x[patch_index];
  fields.cell_dim_y = state.patches.cell_dim_y[patch_index];
  fields.cell_dim_z = state.patches.cell_dim_z[patch_index];
  requireValidAmrPatchMigrationFields(fields, caller);
  return fields;
}

void writePatchFieldsToRow(
    const AmrPatchMigrationFields& fields,
    std::uint32_t first_cell,
    std::uint32_t cell_count,
    std::uint32_t patch_index,
    PatchSoa& patches) {
  requireValidAmrPatchMigrationFields(fields, "writePatchFieldsToRow");
  patches.patch_id[patch_index] = fields.patch_id;
  patches.level[patch_index] = fields.level;
  patches.first_cell[patch_index] = first_cell;
  patches.cell_count[patch_index] = cell_count;
  patches.parent_patch_id[patch_index] = fields.parent_patch_id;
  patches.morton_key[patch_index] = fields.morton_key;
  patches.origin_x_comoving[patch_index] = fields.origin_x_comoving;
  patches.origin_y_comoving[patch_index] = fields.origin_y_comoving;
  patches.origin_z_comoving[patch_index] = fields.origin_z_comoving;
  patches.extent_x_comoving[patch_index] = fields.extent_x_comoving;
  patches.extent_y_comoving[patch_index] = fields.extent_y_comoving;
  patches.extent_z_comoving[patch_index] = fields.extent_z_comoving;
  patches.cell_dim_x[patch_index] = fields.cell_dim_x;
  patches.cell_dim_y[patch_index] = fields.cell_dim_y;
  patches.cell_dim_z[patch_index] = fields.cell_dim_z;
  patches.owning_rank[patch_index] = fields.owning_rank;
}

[[nodiscard]] std::vector<GasCellIdentityRecord> gasIdentityRecordsFromSidecars(
    const CellSoa& cells,
    const GasCellSidecar& gas_cells,
    const PatchSoa& patches) {
  std::vector<GasCellIdentityRecord> records;
  records.reserve(gas_cells.size());
  for (std::uint32_t row = 0; row < gas_cells.size(); ++row) {
    const std::uint32_t patch_row = row < cells.patch_index.size() ? cells.patch_index[row] : 0U;
    std::uint64_t owning_patch_id = 0U;
    if (patch_row < patches.patch_id.size()) {
      owning_patch_id = patches.patch_id[patch_row];
    }
    const std::uint64_t mirrored_parent = gas_cells.parent_particle_id[row];
    records.push_back(GasCellIdentityRecord{
        .gas_cell_id = gas_cells.gas_cell_id[row],
        .parent_particle_id = mirrored_parent == 0U
            ? std::optional<std::uint64_t>{}
            : std::optional<std::uint64_t>{mirrored_parent},
        .owning_patch_id = owning_patch_id,
        .local_cell_row = row,
    });
  }
  return records;
}
}  // namespace

ParticleIdValidationResult validateParticleIdsExact(
    std::span<const std::uint64_t> particle_ids,
    OwnershipValidationWorkspace& scratch) {
  return validateParticleIdsBoundedExactImpl(particle_ids, scratch);
}

std::uint64_t SpeciesContainer::totalCount() const noexcept {
  std::uint64_t total = 0;
  for (const auto count : count_by_species) {
    total += count;
  }
  return total;
}

bool SpeciesContainer::isConsistentWith(const ParticleSidecar& sidecar) const noexcept {
  std::array<std::uint64_t, k_particle_species_count> measured{};
  for (const auto tag : sidecar.species_tag) {
    if (!isValidParticleSpeciesTag(tag)) {
      return false;
    }
    ++measured.at(tag);
  }
  return measured == count_by_species;
}

void ParticleSpeciesIndex::rebuild(const ParticleSidecar& sidecar) {
  representation = Representation::kMaterialized;
  identity_count = 0;
  for (auto& indices : global_index_by_species) {
    indices.clear();
  }

  local_index_by_global.resize(sidecar.size());
  for (std::size_t global_index = 0; global_index < sidecar.size(); ++global_index) {
    const auto tag = sidecar.species_tag[global_index];
    if (!isValidParticleSpeciesTag(tag)) {
      throw std::invalid_argument("ParticleSpeciesIndex.rebuild: invalid species tag");
    }
    auto& species_indices = global_index_by_species[tag];
    local_index_by_global[global_index] = checkedLocalParticleRow(
        species_indices.size(), "ParticleSpeciesIndex::rebuild species-local row");
    species_indices.push_back(checkedLocalParticleRow(
        global_index, "ParticleSpeciesIndex::rebuild global particle row"));
  }
}


void ParticleSpeciesIndex::setHomogeneousDmoIdentity(std::uint32_t particle_count) {
  for (auto& indices : global_index_by_species) {
    AlignedVector<std::uint32_t>{}.swap(indices);
  }
  AlignedVector<std::uint32_t>{}.swap(local_index_by_global);
  representation = Representation::kHomogeneousDmoIdentity;
  identity_count = particle_count;
}

std::size_t ParticleSpeciesIndex::count(ParticleSpecies species) const {
  if (representation == Representation::kHomogeneousDmoIdentity) {
    return species == ParticleSpecies::kDarkMatter ? identity_count : 0U;
  }
  return global_index_by_species[particleSpeciesIndex(species)].size();
}

std::span<const std::uint32_t> ParticleSpeciesIndex::globalIndices(ParticleSpecies species) const {
  if (representation == Representation::kHomogeneousDmoIdentity) {
    if (species == ParticleSpecies::kDarkMatter && identity_count != 0U) {
      throw std::logic_error("ParticleSpeciesIndex.globalIndices: homogeneous DMO identity is implicit; iterate [0,count) instead of materializing an identity vector");
    }
    return {};
  }
  return global_index_by_species[particleSpeciesIndex(species)];
}

std::uint32_t ParticleSpeciesIndex::localIndex(std::uint32_t global_index) const {
  if (representation == Representation::kHomogeneousDmoIdentity) {
    if (global_index >= identity_count) throw std::out_of_range("ParticleSpeciesIndex.localIndex: global index out of range");
    return global_index;
  }
  if (global_index >= local_index_by_global.size()) {
    throw std::out_of_range("ParticleSpeciesIndex.localIndex: global index out of range");
  }
  return local_index_by_global[global_index];
}

std::uint32_t ParticleSpeciesIndex::globalIndex(ParticleSpecies species, std::uint32_t local_index) const {
  if (representation == Representation::kHomogeneousDmoIdentity) {
    if (species != ParticleSpecies::kDarkMatter || local_index >= identity_count) {
      throw std::out_of_range("ParticleSpeciesIndex.globalIndex: local index out of range");
    }
    return local_index;
  }
  const auto& species_indices = global_index_by_species[particleSpeciesIndex(species)];
  if (local_index >= species_indices.size()) {
    throw std::out_of_range("ParticleSpeciesIndex.globalIndex: local index out of range");
  }
  return species_indices[local_index];
}


bool SimulationState::compactHomogeneousDmoMetadata(std::uint32_t local_rank) {
  const std::size_t count = particles.size();
  if (!particles.isConsistent() || !particle_sidecar.isConsistent() ||
      particle_sidecar.particle_id.size() != count ||
      (!particle_sidecar.sfc_key.empty() && particle_sidecar.sfc_key.size() != count) ||
      cells.size() != 0U || star_particles.size() != 0U || black_holes.size() != 0U || tracers.size() != 0U ||
      particle_sidecar.species_tag.size() != count || particle_sidecar.particle_flags.size() != count ||
      particle_sidecar.owning_rank.size() != count || particles.time_bin.size() != count ||
      particle_sidecar.last_drift_time_code.size() != count || particle_sidecar.last_drift_scale_factor.size() != count) {
    return false;
  }
  const std::uint32_t dm = static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter);
  if (count == 0U) {
    m_uniform_species_tag = dm;
    m_uniform_particle_flags = 0U;
    m_uniform_owning_rank = local_rank;
    m_uniform_sfc_key = 0U;
    m_sfc_key_is_uniform = true;
    m_common_last_drift_time_code = 0.0;
    m_common_last_drift_scale_factor = 1.0;
  } else {
    m_uniform_particle_flags = particle_sidecar.particle_flags[0];
    m_common_last_drift_time_code = particle_sidecar.last_drift_time_code[0];
    m_common_last_drift_scale_factor = particle_sidecar.last_drift_scale_factor[0];
    m_uniform_sfc_key = particle_sidecar.sfc_key.empty() ? 0U : particle_sidecar.sfc_key[0];
    m_sfc_key_is_uniform = particle_sidecar.sfc_key.empty() ||
        std::all_of(particle_sidecar.sfc_key.begin(), particle_sidecar.sfc_key.end(), [&](std::uint64_t key) { return key == m_uniform_sfc_key; });
  }
  for (std::size_t i = 0; i < count; ++i) {
    if (particle_sidecar.species_tag[i] != dm || particle_sidecar.owning_rank[i] != local_rank ||
        particle_sidecar.particle_flags[i] != m_uniform_particle_flags ||
        particle_sidecar.last_drift_time_code[i] != m_common_last_drift_time_code ||
        particle_sidecar.last_drift_scale_factor[i] != m_common_last_drift_scale_factor) {
      return false;
    }
  }
  m_particle_metadata_representation = ParticleMetadataRepresentation::kHomogeneousDmo;
  m_uniform_species_tag = dm;
  m_uniform_owning_rank = local_rank;
  // Bin mirrors are independent of uniform species/owner/drift metadata.
  // Keep the one-byte lane only when a heterogeneous scheduler requires it.
  if (std::all_of(particles.time_bin.begin(), particles.time_bin.end(), [](std::uint8_t bin) { return bin == 0U; })) {
    AlignedVector<std::uint8_t>{}.swap(particles.time_bin);
  }
  AlignedVector<std::uint32_t>{}.swap(particle_sidecar.species_tag);
  AlignedVector<std::uint32_t>{}.swap(particle_sidecar.particle_flags);
  AlignedVector<std::uint32_t>{}.swap(particle_sidecar.owning_rank);
  AlignedVector<double>{}.swap(particle_sidecar.last_drift_time_code);
  AlignedVector<double>{}.swap(particle_sidecar.last_drift_scale_factor);
  if (m_sfc_key_is_uniform) {
    AlignedVector<std::uint64_t>{}.swap(particle_sidecar.sfc_key);
  }
  particle_species_index.setHomogeneousDmoIdentity(checkedIntegralNarrow<std::uint32_t>(count, "compactHomogeneousDmoMetadata particle count"));
  species.count_by_species.fill(0U);
  species.count_by_species[particleSpeciesIndex(ParticleSpecies::kDarkMatter)] = count;
  return true;
}

std::uint32_t SimulationState::homogeneousDmoParticleFlags() const {
  if (!hasHomogeneousDmoMetadata()) throw std::logic_error("homogeneousDmoParticleFlags requires compact DMO metadata");
  return m_uniform_particle_flags;
}

double SimulationState::homogeneousDmoLastDriftTimeCode() const {
  if (!hasHomogeneousDmoMetadata()) throw std::logic_error("homogeneousDmoLastDriftTimeCode requires compact DMO metadata");
  return m_common_last_drift_time_code;
}

double SimulationState::homogeneousDmoLastDriftScaleFactor() const {
  if (!hasHomogeneousDmoMetadata()) throw std::logic_error("homogeneousDmoLastDriftScaleFactor requires compact DMO metadata");
  return m_common_last_drift_scale_factor;
}

bool SimulationState::homogeneousDmoSfcKeyIsUniform() const {
  if (!hasHomogeneousDmoMetadata()) throw std::logic_error("homogeneousDmoSfcKeyIsUniform requires compact DMO metadata");
  return m_sfc_key_is_uniform;
}

std::uint64_t SimulationState::homogeneousDmoUniformSfcKey() const {
  if (!hasHomogeneousDmoMetadata()) throw std::logic_error("homogeneousDmoUniformSfcKey requires compact DMO metadata");
  return m_uniform_sfc_key;
}

void SimulationState::commitCompactHomogeneousDmoCandidate(
    ParticleSoa&& candidate_particles,
    ParticleSidecar&& candidate_sidecar,
    std::uint32_t local_rank,
    std::uint32_t uniform_particle_flags,
    double common_last_drift_time_code,
    double common_last_drift_scale_factor,
    bool persistent_particle_ids_certified) {
  const std::size_t count = candidate_particles.size();
  if (!hasHomogeneousDmoMetadata() || cells.size() != 0U || star_particles.size() != 0U ||
      black_holes.size() != 0U || tracers.size() != 0U || !candidate_particles.isConsistent() ||
      !candidate_sidecar.isConsistent() || candidate_sidecar.particle_id.size() != count ||
      !candidate_particles.time_bin.empty() || !candidate_sidecar.species_tag.empty() ||
      !candidate_sidecar.particle_flags.empty() || !candidate_sidecar.owning_rank.empty() ||
      !candidate_sidecar.last_drift_time_code.empty() || !candidate_sidecar.last_drift_scale_factor.empty() ||
      !candidate_sidecar.gravity_softening_comoving.empty() || !candidate_sidecar.has_gravity_softening_override.empty() ||
      (!candidate_sidecar.sfc_key.empty() && candidate_sidecar.sfc_key.size() != count) ||
      !std::isfinite(common_last_drift_time_code) || !std::isfinite(common_last_drift_scale_factor) ||
      common_last_drift_scale_factor <= 0.0) {
    throw std::invalid_argument("compact DMO migration candidate violates canonical representation invariants");
  }
  std::uint64_t uniform_sfc_key = 0U;
  bool sfc_uniform = true;
  if (!candidate_sidecar.sfc_key.empty()) {
    uniform_sfc_key = candidate_sidecar.sfc_key.front();
    sfc_uniform = std::all_of(candidate_sidecar.sfc_key.begin(), candidate_sidecar.sfc_key.end(),
                              [&](std::uint64_t key) { return key == uniform_sfc_key; });
    if (sfc_uniform) AlignedVector<std::uint64_t>{}.swap(candidate_sidecar.sfc_key);
  } else if (count != 0U) {
    uniform_sfc_key = m_uniform_sfc_key;
  }

  using std::swap;
  swap(particles, candidate_particles);
  swap(particle_sidecar, candidate_sidecar);
  m_particle_metadata_representation = ParticleMetadataRepresentation::kHomogeneousDmo;
  m_uniform_species_tag = static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter);
  m_uniform_particle_flags = uniform_particle_flags;
  m_uniform_owning_rank = local_rank;
  m_uniform_sfc_key = uniform_sfc_key;
  m_sfc_key_is_uniform = sfc_uniform;
  m_common_last_drift_time_code = common_last_drift_time_code;
  m_common_last_drift_scale_factor = common_last_drift_scale_factor;
  particle_species_index.setHomogeneousDmoIdentity(
      checkedIntegralNarrow<std::uint32_t>(count, "commitCompactHomogeneousDmoCandidate particle count"));
  species.count_by_species.fill(0U);
  species.count_by_species[particleSpeciesIndex(ParticleSpecies::kDarkMatter)] = count;
  bumpParticleIdentityGeneration();
  if (persistent_particle_ids_certified) {
    m_particle_identity_validation_certificate = ParticleIdentityValidationCertificate{
        .identity_generation = m_particle_identity_generation,
        .particle_count = particle_sidecar.particle_id.size(),
        .local_uniqueness_proven = true,
        .nonzero_ids_proven = true,
    };
  }
  bumpParticleIndexGeneration();
  bumpGravitySourceGeneration();
}

void SimulationState::materializeParticleMetadata() {
  if (!hasHomogeneousDmoMetadata()) return;
  const std::size_t count = particles.size();
  if (particles.time_bin.empty()) particles.time_bin.assign(count, 0U);
  particle_sidecar.species_tag.assign(count, m_uniform_species_tag);
  particle_sidecar.particle_flags.assign(count, m_uniform_particle_flags);
  particle_sidecar.owning_rank.assign(count, m_uniform_owning_rank);
  particle_sidecar.last_drift_time_code.assign(count, m_common_last_drift_time_code);
  particle_sidecar.last_drift_scale_factor.assign(count, m_common_last_drift_scale_factor);
  if (particle_sidecar.sfc_key.empty()) particle_sidecar.sfc_key.assign(count, m_uniform_sfc_key);
  m_particle_metadata_representation = ParticleMetadataRepresentation::kMaterializedGeneric;
  particle_species_index.rebuild(particle_sidecar);
}

std::uint32_t SimulationState::particleSpeciesTag(std::size_t particle_index) const {
  if (particle_index >= particles.size()) throw std::out_of_range("particleSpeciesTag: particle index out of range");
  return hasHomogeneousDmoMetadata() ? m_uniform_species_tag : particle_sidecar.species_tag[particle_index];
}
std::uint32_t SimulationState::particleFlags(std::size_t particle_index) const {
  if (particle_index >= particles.size()) throw std::out_of_range("particleFlags: particle index out of range");
  return hasHomogeneousDmoMetadata() ? m_uniform_particle_flags : particle_sidecar.particle_flags[particle_index];
}
std::uint32_t SimulationState::particleOwningRank(std::size_t particle_index) const {
  if (particle_index >= particles.size()) throw std::out_of_range("particleOwningRank: particle index out of range");
  return hasHomogeneousDmoMetadata() ? m_uniform_owning_rank : particle_sidecar.owning_rank[particle_index];
}
std::uint64_t SimulationState::particleSfcKey(std::size_t particle_index) const {
  if (particle_index >= particles.size()) throw std::out_of_range("particleSfcKey: particle index out of range");
  if (hasHomogeneousDmoMetadata() && m_sfc_key_is_uniform) return m_uniform_sfc_key;
  return particle_sidecar.sfc_key[particle_index];
}
double SimulationState::particleLastDriftTimeCode(std::size_t particle_index) const {
  if (particle_index >= particles.size()) throw std::out_of_range("particleLastDriftTimeCode: particle index out of range");
  return hasHomogeneousDmoMetadata() ? m_common_last_drift_time_code : particle_sidecar.last_drift_time_code[particle_index];
}
double SimulationState::particleLastDriftScaleFactor(std::size_t particle_index) const {
  if (particle_index >= particles.size()) throw std::out_of_range("particleLastDriftScaleFactor: particle index out of range");
  return hasHomogeneousDmoMetadata() ? m_common_last_drift_scale_factor : particle_sidecar.last_drift_scale_factor[particle_index];
}
std::uint8_t SimulationState::particleTimeBin(std::size_t particle_index) const {
  if (particle_index >= particles.size()) throw std::out_of_range("particleTimeBin: particle index out of range");
  return particles.time_bin.empty() && hasHomogeneousDmoMetadata() ? 0U : particles.time_bin[particle_index];
}
void SimulationState::updateAllParticleDriftEpoch(double time_code, double scale_factor) {
  if (!std::isfinite(time_code) || !std::isfinite(scale_factor) || scale_factor <= 0.0) {
    throw std::invalid_argument("updateAllParticleDriftEpoch: invalid drift epoch");
  }
  if (hasHomogeneousDmoMetadata()) {
    m_common_last_drift_time_code = time_code;
    m_common_last_drift_scale_factor = scale_factor;
    return;
  }
  std::fill(particle_sidecar.last_drift_time_code.begin(), particle_sidecar.last_drift_time_code.end(), time_code);
  std::fill(particle_sidecar.last_drift_scale_factor.begin(), particle_sidecar.last_drift_scale_factor.end(), scale_factor);
}

void SimulationState::rebuildSpeciesIndex() {
  if (hasHomogeneousDmoMetadata()) {
    particle_species_index.setHomogeneousDmoIdentity(checkedIntegralNarrow<std::uint32_t>(particles.size(), "rebuildSpeciesIndex DMO identity"));
    return;
  }
  particle_species_index.rebuild(particle_sidecar);
  const auto gas_count = particle_species_index.count(ParticleSpecies::kGas);
  const bool identity_uninitialized = gas_cells.gas_cell_id.empty() ||
      gas_cells.parent_particle_id.empty() ||
      (std::all_of(gas_cells.gas_cell_id.begin(), gas_cells.gas_cell_id.end(), [](std::uint64_t id) { return id == 0; }) &&
       std::all_of(gas_cells.parent_particle_id.begin(), gas_cells.parent_particle_id.end(), [](std::uint64_t id) { return id == 0; }));
  if (gas_count > 0 && gas_count == cells.size() && gas_count == gas_cells.size() && identity_uninitialized) {
    refreshGasCellIdentityFromParticleOrder();
  }
}

bool SimulationState::validateUniqueParticleIds() const {
  if (hasCertifiedUniqueParticleIds()) return true;
  OwnershipValidationWorkspace scratch;
  return validateUniqueParticleIds(scratch);
}

bool SimulationState::validateUniqueParticleIds(
    OwnershipValidationWorkspace& scratch) const {
  if (hasCertifiedUniqueParticleIds()) return true;
  const ParticleIdValidationResult validation =
      validateParticleIdsExact(particle_sidecar.particle_id, scratch);
  if (!validation.unique) return false;
  m_particle_identity_validation_certificate = ParticleIdentityValidationCertificate{
      .identity_generation = m_particle_identity_generation,
      .particle_count = particle_sidecar.particle_id.size(),
      .local_uniqueness_proven = true,
      .nonzero_ids_proven = validation.nonzero,
  };
  return true;
}

bool SimulationState::validatePersistentParticleIds() const {
  if (hasCertifiedPersistentParticleIds()) return true;
  OwnershipValidationWorkspace scratch;
  return validatePersistentParticleIds(scratch);
}

bool SimulationState::validatePersistentParticleIds(
    OwnershipValidationWorkspace& scratch) const {
  if (hasCertifiedPersistentParticleIds()) return true;
  const ParticleIdValidationResult validation =
      validateParticleIdsExact(particle_sidecar.particle_id, scratch);
  if (!validation.unique || !validation.nonzero) return false;
  m_particle_identity_validation_certificate = ParticleIdentityValidationCertificate{
      .identity_generation = m_particle_identity_generation,
      .particle_count = particle_sidecar.particle_id.size(),
      .local_uniqueness_proven = true,
      .nonzero_ids_proven = true,
  };
  return true;
}

ParticleTransferPacket SimulationState::packSpeciesTransferPacket(ParticleSpecies species_tag) const {
  const auto indices = particle_species_index.globalIndices(species_tag);

  ParticleTransferPacket packet;
  packet.species = species_tag;
  packet.particle_id.resize(indices.size());
  packet.position_x_comoving.resize(indices.size());
  packet.position_y_comoving.resize(indices.size());
  packet.position_z_comoving.resize(indices.size());
  packet.velocity_x_peculiar.resize(indices.size());
  packet.velocity_y_peculiar.resize(indices.size());
  packet.velocity_z_peculiar.resize(indices.size());
  packet.mass_code.resize(indices.size());
  packet.time_bin.resize(indices.size());
  packet.owning_rank.resize(indices.size());
  packet.last_drift_time_code.resize(indices.size());
  packet.last_drift_scale_factor.resize(indices.size());

  for (std::size_t i = 0; i < indices.size(); ++i) {
    const auto source = indices[i];
    packet.particle_id[i] = particle_sidecar.particle_id[source];
    packet.position_x_comoving[i] = particles.position_x_comoving[source];
    packet.position_y_comoving[i] = particles.position_y_comoving[source];
    packet.position_z_comoving[i] = particles.position_z_comoving[source];
    packet.velocity_x_peculiar[i] = particles.velocity_x_peculiar[source];
    packet.velocity_y_peculiar[i] = particles.velocity_y_peculiar[source];
    packet.velocity_z_peculiar[i] = particles.velocity_z_peculiar[source];
    packet.mass_code[i] = particles.mass_code[source];
    packet.time_bin[i] = particleTimeBin(source);
    packet.owning_rank[i] = particleOwningRank(source);
    packet.last_drift_time_code[i] = particleLastDriftTimeCode(source);
    packet.last_drift_scale_factor[i] = particleLastDriftScaleFactor(source);
  }

  return packet;
}

std::vector<ParticleMigrationRecord> SimulationState::packParticleMigrationRecordsCore(
    std::span<const std::uint32_t> local_indices) const {
  if (cells.size() != 0U || gas_cells.size() != 0U) {
    this->requireGasCellIdentityMapCoversDenseRows("packParticleMigrationRecords");
  }
  std::vector<ParticleMigrationRecord> records;
  records.reserve(local_indices.size());
  for (const auto index : local_indices) {
    if (index >= particles.size()) {
      throw std::out_of_range("packParticleMigrationRecords: local index out of range");
    }
  }

  for (const auto index : local_indices) {
    ParticleMigrationRecord record;
    record.particle_id = particle_sidecar.particle_id[index];
    record.sfc_key = particleSfcKey(index);
    record.species_tag = particleSpeciesTag(index);
    record.particle_flags = particleFlags(index);
    record.owning_rank = particleOwningRank(index);
    record.last_drift_time_code = particleLastDriftTimeCode(index);
    record.last_drift_scale_factor = particleLastDriftScaleFactor(index);
    record.position_x_comoving = particles.position_x_comoving[index];
    record.position_y_comoving = particles.position_y_comoving[index];
    record.position_z_comoving = particles.position_z_comoving[index];
    record.velocity_x_peculiar = particles.velocity_x_peculiar[index];
    record.velocity_y_peculiar = particles.velocity_y_peculiar[index];
    record.velocity_z_peculiar = particles.velocity_z_peculiar[index];
    record.mass_code = particles.mass_code[index];
    record.time_bin = particleTimeBin(index);
    if (!particle_sidecar.gravity_softening_comoving.empty()) {
      record.has_gravity_softening_value = true;
      record.gravity_softening_comoving = particle_sidecar.gravity_softening_comoving[index];
    }
    if (particle_sidecar.hasGravitySofteningOverride(index)) {
      record.has_gravity_softening_value = true;
      record.has_gravity_softening_override = true;
    }

    if (record.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kGas)) {
      if (cells.size() != 0 || gas_cells.size() != 0) {
        const std::vector<std::uint32_t> rows = gas_cell_identity.rowsForParentParticleId(record.particle_id);
        if (rows.size() == 1U) {
          record.has_gas_cell_fields = true;
          record.gas_cell_fields = gasCellFieldsFromLocalRow(
              *this,
              rows.front(),
              0U,
              "packParticleMigrationRecords");
        }
      }
    } else if (record.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
      const std::size_t row =
          findRequiredSidecarRow(star_particles.particle_index, index, "star_particles");
      record.has_star_fields = true;
      record.star_fields.formation_scale_factor = star_particles.formation_scale_factor[row];
      record.star_fields.birth_mass_code = star_particles.birth_mass_code[row];
      record.star_fields.metallicity_mass_fraction = star_particles.metallicity_mass_fraction[row];
      record.star_fields.birth_key = star_particles.birth_key[row];
      record.star_fields.parent_gas_cell_id = star_particles.parent_gas_cell_id[row];
      record.star_fields.birth_tick = star_particles.birth_tick[row];
      record.star_fields.birth_ordinal = star_particles.birth_ordinal[row];
      record.star_fields.stellar_age_years_last = star_particles.stellar_age_years_last[row];
      record.star_fields.stellar_returned_mass_cumulative_code =
          star_particles.stellar_returned_mass_cumulative_code[row];
      record.star_fields.stellar_returned_metals_cumulative_code =
          star_particles.stellar_returned_metals_cumulative_code[row];
      record.star_fields.stellar_newly_synthesized_metals_cumulative_code =
          star_particles.stellar_newly_synthesized_metals_cumulative_code[row];
      record.star_fields.stellar_feedback_energy_cumulative_erg =
          star_particles.stellar_feedback_energy_cumulative_erg[row];
      record.star_fields.enrichment_carry_mass_code = star_particles.enrichment_carry_mass_code[row];
      record.star_fields.enrichment_carry_metals_code = star_particles.enrichment_carry_metals_code[row];
      record.star_fields.enrichment_carry_feedback_energy_erg =
          star_particles.enrichment_carry_feedback_energy_erg[row];
      record.star_fields.enrichment_carry_momentum_code =
          star_particles.enrichment_carry_momentum_code[row];
      record.star_fields.stellar_deposited_mass_cumulative_code =
          star_particles.stellar_deposited_mass_cumulative_code[row];
      record.star_fields.stellar_deposited_metals_cumulative_code =
          star_particles.stellar_deposited_metals_cumulative_code[row];
      record.star_fields.stellar_deposited_feedback_energy_cumulative_erg =
          star_particles.stellar_deposited_feedback_energy_cumulative_erg[row];
      for (std::size_t channel = 0; channel < record.star_fields.stellar_returned_mass_channel_cumulative_code.size();
           ++channel) {
        record.star_fields.stellar_returned_mass_channel_cumulative_code[channel] =
            star_particles.stellar_returned_mass_channel_cumulative_code[channel][row];
        record.star_fields.stellar_returned_metals_channel_cumulative_code[channel] =
            star_particles.stellar_returned_metals_channel_cumulative_code[channel][row];
        record.star_fields.stellar_feedback_energy_channel_cumulative_erg[channel] =
            star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][row];
      }
    } else if (record.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
      const std::size_t row =
          findRequiredSidecarRow(black_holes.particle_index, index, "black_holes");
      record.has_black_hole_fields = true;
      record.black_hole_fields.host_cell_index = black_holes.host_cell_index[row];
      record.black_hole_fields.subgrid_mass_code = black_holes.subgrid_mass_code[row];
      record.black_hole_fields.accretion_rate_code = black_holes.accretion_rate_code[row];
      record.black_hole_fields.feedback_energy_code = black_holes.feedback_energy_code[row];
      record.black_hole_fields.eddington_ratio = black_holes.eddington_ratio[row];
      record.black_hole_fields.cumulative_accreted_mass_code = black_holes.cumulative_accreted_mass_code[row];
      record.black_hole_fields.cumulative_feedback_energy_code = black_holes.cumulative_feedback_energy_code[row];
      record.black_hole_fields.duty_cycle_active_time_code = black_holes.duty_cycle_active_time_code[row];
      record.black_hole_fields.duty_cycle_total_time_code = black_holes.duty_cycle_total_time_code[row];
    } else if (record.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
      const std::size_t row = findRequiredSidecarRow(tracers.particle_index, index, "tracers");
      record.has_tracer_fields = true;
      record.tracer_fields.parent_particle_id = tracers.parent_particle_id[row];
      record.tracer_fields.injection_step = tracers.injection_step[row];
      record.tracer_fields.host_cell_index = tracers.host_cell_index[row];
      record.tracer_fields.mass_fraction_of_host = tracers.mass_fraction_of_host[row];
      record.tracer_fields.last_host_mass_code = tracers.last_host_mass_code[row];
      record.tracer_fields.cumulative_exchanged_mass_code = tracers.cumulative_exchanged_mass_code[row];
    }

    for (const ModuleSidecarBlock* block : sidecars.blocksSortedByName()) {
      if (!block->isParticleIndexed()) {
        continue;
      }
      const auto row = findModuleSidecarRowByParticleId(*block, record.particle_id);
      if (!row.has_value()) {
        if (moduleRequiresRecord(*block, record, "packParticleMigrationRecords")) {
          throw std::invalid_argument(
              "packParticleMigrationRecords: required module sidecar row is missing for migrating particle");
        }
        continue;
      }
      const auto row_payload = block->rowPayload(*row);
      ModuleParticleSidecarPayload module_payload;
      module_payload.module_name = block->module_name;
      module_payload.schema_version = block->schema_version;
      module_payload.row_stride_bytes = block->row_stride_bytes;
      module_payload.required_species_mask = block->required_species_mask;
      module_payload.requirement = block->requirement;
      module_payload.payload.assign(row_payload.begin(), row_payload.end());
      record.module_sidecar_payloads.push_back(std::move(module_payload));
    }

    records.push_back(record);
  }
  return records;
}

void SimulationState::commitParticleMigration(const ParticleMigrationCommit& commit) {
  if (commit.world_rank < 0) {
    throw std::invalid_argument("commitParticleMigration: world_rank must be non-negative");
  }
  const std::size_t particle_count = particles.size();
  std::vector<std::uint8_t> outbound_mask(particle_count, 0U);
  for (const auto index : commit.outbound_local_indices) {
    if (index >= particle_count) {
      throw std::out_of_range("commitParticleMigration: outbound index out of range");
    }
    if (outbound_mask[index] != 0U) {
      throw std::invalid_argument("commitParticleMigration: duplicate outbound index");
    }
    if (particle_sidecar.owning_rank[index] != static_cast<std::uint32_t>(commit.world_rank)) {
      throw std::invalid_argument("commitParticleMigration: outbound index is not owned by committing rank");
    }
    outbound_mask[index] = 1U;
  }

  std::vector<std::uint8_t> stale_ghost_mask(particle_count, 0U);
  for (const auto index : commit.stale_local_ghost_indices) {
    if (index >= particle_count) {
      throw std::out_of_range("commitParticleMigration: stale ghost index out of range");
    }
    if (stale_ghost_mask[index] != 0U) {
      throw std::invalid_argument("commitParticleMigration: duplicate stale ghost index");
    }
    if (particle_sidecar.owning_rank[index] == static_cast<std::uint32_t>(commit.world_rank)) {
      throw std::invalid_argument("commitParticleMigration: stale ghost index is owned by committing rank");
    }
    stale_ghost_mask[index] = 1U;
  }

  std::vector<std::uint8_t> remove_mask = outbound_mask;
  for (std::size_t i = 0; i < remove_mask.size(); ++i) {
    if (stale_ghost_mask[i] != 0U) {
      remove_mask[i] = 1U;
    }
  }

  const bool destination_expects_softening_value =
      !particle_sidecar.gravity_softening_comoving.empty() ||
      !particle_sidecar.has_gravity_softening_override.empty();
  const bool state_has_gas_state = cells.size() != 0 || gas_cells.size() != 0;
  if (state_has_gas_state) {
    this->requireGasCellIdentityMapCoversDenseRows("commitParticleMigration");
  }

  auto require_inbound_sidecar_contract = [destination_expects_softening_value, &commit](
                                           const ParticleMigrationRecord& inbound) {
    if (!isValidParticleSpeciesTag(inbound.species_tag)) {
      throw std::invalid_argument("commitParticleMigration: inbound record has invalid species tag");
    }

    const auto particle_species = static_cast<ParticleSpecies>(inbound.species_tag);
    const bool requires_star_fields = particle_species == ParticleSpecies::kStar;
    const bool requires_black_hole_fields = particle_species == ParticleSpecies::kBlackHole;
    const bool requires_tracer_fields = particle_species == ParticleSpecies::kTracer;

    if (inbound.has_gas_cell_fields && particle_species != ParticleSpecies::kGas) {
      throw std::invalid_argument(
          "commitParticleMigration: inbound gas-cell hydro fields do not match species tag");
    }
    if (commit.preserve_gas_cell_state && inbound.has_gas_cell_fields) {
      throw std::invalid_argument(
          "commitParticleMigration: gas-cell payloads are forbidden when preserve_gas_cell_state is set");
    }
    if (inbound.has_gas_cell_fields) {
      requireValidGasCellMigrationFields(inbound.gas_cell_fields, "commitParticleMigration");
      if (inbound.gas_cell_fields.has_parent_particle == 0U ||
          inbound.gas_cell_fields.parent_particle_id != inbound.particle_id) {
        throw std::invalid_argument(
            "commitParticleMigration: particle-embedded gas-cell payload must name the migrating parent particle");
      }
    }
    if (inbound.has_star_fields != requires_star_fields) {
      throw std::invalid_argument(
          "commitParticleMigration: inbound star sidecar fields do not match species tag");
    }
    if (inbound.has_black_hole_fields != requires_black_hole_fields) {
      throw std::invalid_argument(
          "commitParticleMigration: inbound black-hole sidecar fields do not match species tag");
    }
    if (inbound.has_tracer_fields != requires_tracer_fields) {
      throw std::invalid_argument(
          "commitParticleMigration: inbound tracer sidecar fields do not match species tag");
    }
    if (destination_expects_softening_value && !inbound.has_gravity_softening_value) {
      throw std::invalid_argument(
          "commitParticleMigration: inbound record is missing required softening value for destination sidecar");
    }
    if (inbound.has_gravity_softening_override && !inbound.has_gravity_softening_value) {
      throw std::invalid_argument(
          "commitParticleMigration: inbound softening override requires an explicit softening value");
    }
    if (inbound.has_scheduler_fields && inbound.scheduler_fields.bin_index != inbound.time_bin) {
      throw std::invalid_argument(
          "commitParticleMigration: inbound time_bin mirror is stale relative to scheduler migration fields");
    }
  };

  std::unordered_set<std::uint64_t> final_particle_ids;
  final_particle_ids.reserve(particle_count + commit.inbound_records.size());
  for (std::size_t i = 0; i < particle_count; ++i) {
    if (remove_mask[i] != 0U) {
      continue;
    }
    if (particle_sidecar.particle_id[i] == 0U) {
      throw std::invalid_argument("commitParticleMigration: kept particle has zero persistent ID");
    }
    if (!final_particle_ids.insert(particle_sidecar.particle_id[i]).second) {
      throw std::invalid_argument("commitParticleMigration: kept particles contain duplicate IDs");
    }
  }
  std::unordered_map<std::string, std::unordered_map<std::uint64_t, const ModuleParticleSidecarPayload*>>
      inbound_module_payloads_by_module;
  for (const auto& inbound : commit.inbound_records) {
    require_inbound_sidecar_contract(inbound);
    if (inbound.owning_rank != static_cast<std::uint32_t>(commit.world_rank)) {
      throw std::invalid_argument("commitParticleMigration: inbound record ownership must equal commit world rank");
    }
    if (inbound.particle_id == 0U) {
      throw std::invalid_argument("commitParticleMigration: inbound record has zero persistent particle ID");
    }
    if (!final_particle_ids.insert(inbound.particle_id).second) {
      throw std::invalid_argument("commitParticleMigration: inbound record would create duplicate particle ID");
    }

    std::unordered_set<std::string> modules_seen_for_particle;
    for (const ModuleParticleSidecarPayload& payload : inbound.module_sidecar_payloads) {
      if (payload.module_name.empty()) {
        throw std::invalid_argument("commitParticleMigration: module sidecar payload has empty module_name");
      }
      if (payload.row_stride_bytes == 0U ||
          payload.payload.size() != static_cast<std::size_t>(payload.row_stride_bytes)) {
        throw std::invalid_argument("commitParticleMigration: module sidecar payload row shape is invalid");
      }
      if (!modules_seen_for_particle.insert(payload.module_name).second) {
        throw std::invalid_argument("commitParticleMigration: duplicate module sidecar payload for one particle");
      }
      if (const ModuleSidecarBlock* local_block = sidecars.find(payload.module_name)) {
        if (!local_block->isParticleIndexed()) {
          throw std::invalid_argument("commitParticleMigration: inbound particle module sidecar targets a non-particle-indexed local block");
        }
        if (local_block->schema_version != payload.schema_version ||
            local_block->row_stride_bytes != payload.row_stride_bytes ||
            local_block->required_species_mask != payload.required_species_mask ||
            local_block->requirement.kind != payload.requirement.kind ||
            local_block->requirement.species_mask != payload.requirement.species_mask ||
            local_block->requirement.particle_flags_mask != payload.requirement.particle_flags_mask ||
            local_block->requirement.threshold_code != payload.requirement.threshold_code) {
          throw std::invalid_argument("commitParticleMigration: inbound module sidecar schema, stride, or required-coverage mismatch");
        }
      }
      auto& rows_by_particle = inbound_module_payloads_by_module[payload.module_name];
      if (!rows_by_particle.emplace(inbound.particle_id, &payload).second) {
        throw std::invalid_argument("commitParticleMigration: duplicate inbound module sidecar payload for particle ID");
      }
    }
  }

  std::vector<std::uint32_t> old_to_new(particle_count, std::numeric_limits<std::uint32_t>::max());
  std::size_t kept_count = 0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    if (remove_mask[i] == 0U) {
      old_to_new[i] = checkedLocalParticleRow(
          kept_count++, "commitParticleMigration kept particle row");
    }
  }

  const std::size_t final_count = kept_count + commit.inbound_records.size();
  ParticleSoa new_particles;
  new_particles.resize(final_count);
  ParticleSidecar new_sidecar;
  new_sidecar.resize(final_count);

  bool has_softening_sidecar = !particle_sidecar.gravity_softening_comoving.empty();
  bool has_softening_override_mask = !particle_sidecar.has_gravity_softening_override.empty();
  for (const auto& inbound : commit.inbound_records) {
    if (inbound.has_gravity_softening_value) {
      has_softening_sidecar = true;
    }
    if (inbound.has_gravity_softening_override) {
      has_softening_sidecar = true;
      has_softening_override_mask = true;
    }
  }
  if (has_softening_sidecar) {
    new_sidecar.gravity_softening_comoving.resize(final_count, 0.0);
  }
  if (has_softening_override_mask) {
    new_sidecar.has_gravity_softening_override.resize(final_count, 0U);
  }

  std::size_t write_index = 0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    if (remove_mask[i] != 0U) {
      continue;
    }
    new_particles.position_x_comoving[write_index] = particles.position_x_comoving[i];
    new_particles.position_y_comoving[write_index] = particles.position_y_comoving[i];
    new_particles.position_z_comoving[write_index] = particles.position_z_comoving[i];
    new_particles.velocity_x_peculiar[write_index] = particles.velocity_x_peculiar[i];
    new_particles.velocity_y_peculiar[write_index] = particles.velocity_y_peculiar[i];
    new_particles.velocity_z_peculiar[write_index] = particles.velocity_z_peculiar[i];
    new_particles.mass_code[write_index] = particles.mass_code[i];
    new_particles.time_bin[write_index] = particles.time_bin[i];
    new_sidecar.particle_id[write_index] = particle_sidecar.particle_id[i];
    new_sidecar.sfc_key[write_index] = particle_sidecar.sfc_key[i];
    new_sidecar.species_tag[write_index] = particle_sidecar.species_tag[i];
    new_sidecar.particle_flags[write_index] = particle_sidecar.particle_flags[i];
    new_sidecar.owning_rank[write_index] = particle_sidecar.owning_rank[i];
    new_sidecar.last_drift_time_code[write_index] = particle_sidecar.last_drift_time_code[i];
    new_sidecar.last_drift_scale_factor[write_index] = particle_sidecar.last_drift_scale_factor[i];
    if (has_softening_sidecar) {
      new_sidecar.gravity_softening_comoving[write_index] =
          particle_sidecar.gravity_softening_comoving.empty() ? 0.0 : particle_sidecar.gravity_softening_comoving[i];
    }
    if (has_softening_override_mask && particle_sidecar.hasGravitySofteningOverride(i)) {
      new_sidecar.has_gravity_softening_override[write_index] = 1U;
    }
    ++write_index;
  }

  for (const auto& inbound : commit.inbound_records) {
    new_particles.position_x_comoving[write_index] = inbound.position_x_comoving;
    new_particles.position_y_comoving[write_index] = inbound.position_y_comoving;
    new_particles.position_z_comoving[write_index] = inbound.position_z_comoving;
    new_particles.velocity_x_peculiar[write_index] = inbound.velocity_x_peculiar;
    new_particles.velocity_y_peculiar[write_index] = inbound.velocity_y_peculiar;
    new_particles.velocity_z_peculiar[write_index] = inbound.velocity_z_peculiar;
    new_particles.mass_code[write_index] = inbound.mass_code;
    new_particles.time_bin[write_index] = inbound.time_bin;
    new_sidecar.particle_id[write_index] = inbound.particle_id;
    new_sidecar.sfc_key[write_index] = inbound.sfc_key;
    new_sidecar.species_tag[write_index] = inbound.species_tag;
    new_sidecar.particle_flags[write_index] = inbound.particle_flags;
    new_sidecar.owning_rank[write_index] = inbound.owning_rank;
    new_sidecar.last_drift_time_code[write_index] = inbound.last_drift_time_code;
    new_sidecar.last_drift_scale_factor[write_index] = inbound.last_drift_scale_factor;
    if (has_softening_sidecar && inbound.has_gravity_softening_value) {
      new_sidecar.gravity_softening_comoving[write_index] = inbound.gravity_softening_comoving;
    }
    if (has_softening_override_mask && inbound.has_gravity_softening_override) {
      new_sidecar.has_gravity_softening_override[write_index] = 1U;
    }
    ++write_index;
  }

  std::unordered_map<std::uint64_t, const ParticleMigrationRecord*> inbound_record_by_particle_id;
  inbound_record_by_particle_id.reserve(commit.inbound_records.size());
  bool has_inbound_gas_cell_fields = false;
  for (const auto& inbound : commit.inbound_records) {
    inbound_record_by_particle_id.emplace(inbound.particle_id, &inbound);
    has_inbound_gas_cell_fields = has_inbound_gas_cell_fields || inbound.has_gas_cell_fields;
  }

  CellSoa rebuilt_cells;
  GasCellSidecar rebuilt_gas_cells;
  // Composition rows are rebuilt alongside gas rows (identity-carried payload).
  GasCompositionSidecar rebuilt_composition = gas_composition.emptyLike(0U);
  std::vector<std::uint32_t> old_cell_to_new(cells.size(), kInvalidGasCellRow);
  const bool rebuild_gas_state =
      !commit.preserve_gas_cell_state && (state_has_gas_state || has_inbound_gas_cell_fields);
  if (rebuild_gas_state) {
    std::unordered_set<std::uint64_t> removed_particle_ids;
    removed_particle_ids.reserve(particle_count);
    for (std::size_t i = 0; i < particle_count; ++i) {
      if (remove_mask[i] != 0U) {
        removed_particle_ids.insert(particle_sidecar.particle_id[i]);
      }
    }

    std::vector<GasCellMigrationFields> final_gas_fields;
    final_gas_fields.reserve(checkedSizeAdd(
        cells.size(), commit.inbound_records.size(), "commitParticleMigration gas-cell reserve"));
    std::unordered_set<std::uint64_t> final_gas_cell_ids;
    final_gas_cell_ids.reserve(checkedSizeAdd(
        cells.size(), commit.inbound_records.size(), "migration gas-cell id reserve"));
    for (std::uint32_t row = 0; row < cells.size(); ++row) {
      GasCellMigrationFields fields = gasCellFieldsFromLocalRow(*this, row, 0U, "commitParticleMigration");
      if (fields.has_parent_particle != 0U && removed_particle_ids.contains(fields.parent_particle_id)) {
        continue;
      }
      const std::uint32_t new_row = checkedLocalCellRow(
          final_gas_fields.size(), "commitParticleMigration kept gas-cell row");
      old_cell_to_new[row] = new_row;
      fields.destination_local_cell_row = new_row;
      if (!final_gas_cell_ids.insert(fields.gas_cell_id).second) {
        throw std::invalid_argument("commitParticleMigration: duplicate kept gas_cell_id");
      }
      final_gas_fields.push_back(fields);
    }
    for (const auto& inbound : commit.inbound_records) {
      if (!inbound.has_gas_cell_fields) {
        continue;
      }
      GasCellMigrationFields fields = inbound.gas_cell_fields;
      fields.destination_local_cell_row = checkedLocalCellRow(
          final_gas_fields.size(), "commitParticleMigration inbound gas-cell row");
      if (!final_gas_cell_ids.insert(fields.gas_cell_id).second) {
        throw std::invalid_argument("commitParticleMigration: inbound gas_cell_id duplicates local gas-cell state");
      }
      final_gas_fields.push_back(fields);
    }

    rebuilt_cells.resize(final_gas_fields.size());
    rebuilt_gas_cells.resize(final_gas_fields.size());
    rebuilt_composition = gas_composition.emptyLike(final_gas_fields.size());
    std::size_t active_sparse_rows = 0U;
    for (const auto& fields : final_gas_fields) {
      active_sparse_rows +=
          migrationPayloadHasExplicitSparseRow(fields, rebuilt_composition) ? 1U : 0U;
    }
    prepareAdmittedSparseCompositionCandidate(
        rebuilt_composition, active_sparse_rows,
        "particle_migration.adaptive_chemistry_rows");
    for (std::uint32_t row = 0; row < final_gas_fields.size(); ++row) {
      writeGasCellFieldsToRow(final_gas_fields[row], row, rebuilt_cells, rebuilt_gas_cells, rebuilt_composition);
    }
  }

  auto rebuildStarSidecar = [&](StarParticleSidecar* destination) {
    std::vector<std::uint32_t> kept_rows;
    for (std::size_t row = 0; row < star_particles.size(); ++row) {
      const auto old_particle = star_particles.particle_index[row];
      if (old_particle >= remove_mask.size()) {
        throw std::runtime_error("commitParticleMigration: stale star sidecar index");
      }
      if (remove_mask[old_particle] == 0U) {
        kept_rows.push_back(checkedLocalParticleRow(
            row, "commitParticleMigration sidecar row"));
      }
    }
    const std::size_t inbound_count = static_cast<std::size_t>(std::count_if(
        commit.inbound_records.begin(), commit.inbound_records.end(),
        [](const ParticleMigrationRecord& record) { return record.has_star_fields; }));
    destination->resize(checkedSizeAdd(kept_rows.size(), inbound_count, "commitParticleMigration star rows"));
    for (std::size_t row = 0; row < kept_rows.size(); ++row) {
      const std::size_t source = kept_rows[row];
      const auto old_particle = star_particles.particle_index[source];
      destination->particle_index[row] = old_to_new[old_particle];
      destination->formation_scale_factor[row] = star_particles.formation_scale_factor[source];
      destination->birth_mass_code[row] = star_particles.birth_mass_code[source];
      destination->metallicity_mass_fraction[row] = star_particles.metallicity_mass_fraction[source];
      destination->birth_key[row] = star_particles.birth_key[source];
      destination->parent_gas_cell_id[row] = star_particles.parent_gas_cell_id[source];
      destination->birth_tick[row] = star_particles.birth_tick[source];
      destination->birth_ordinal[row] = star_particles.birth_ordinal[source];
      destination->stellar_age_years_last[row] = star_particles.stellar_age_years_last[source];
      destination->stellar_returned_mass_cumulative_code[row] =
          star_particles.stellar_returned_mass_cumulative_code[source];
      destination->stellar_returned_metals_cumulative_code[row] =
          star_particles.stellar_returned_metals_cumulative_code[source];
      destination->stellar_newly_synthesized_metals_cumulative_code[row] =
          star_particles.stellar_newly_synthesized_metals_cumulative_code[source];
      destination->stellar_feedback_energy_cumulative_erg[row] =
          star_particles.stellar_feedback_energy_cumulative_erg[source];
      destination->enrichment_carry_mass_code[row] = star_particles.enrichment_carry_mass_code[source];
      destination->enrichment_carry_metals_code[row] = star_particles.enrichment_carry_metals_code[source];
      destination->enrichment_carry_feedback_energy_erg[row] =
          star_particles.enrichment_carry_feedback_energy_erg[source];
      destination->enrichment_carry_momentum_code[row] =
          star_particles.enrichment_carry_momentum_code[source];
      destination->stellar_deposited_mass_cumulative_code[row] =
          star_particles.stellar_deposited_mass_cumulative_code[source];
      destination->stellar_deposited_metals_cumulative_code[row] =
          star_particles.stellar_deposited_metals_cumulative_code[source];
      destination->stellar_deposited_feedback_energy_cumulative_erg[row] =
          star_particles.stellar_deposited_feedback_energy_cumulative_erg[source];
      for (std::size_t channel = 0; channel < destination->stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
        destination->stellar_returned_mass_channel_cumulative_code[channel][row] =
            star_particles.stellar_returned_mass_channel_cumulative_code[channel][source];
        destination->stellar_returned_metals_channel_cumulative_code[channel][row] =
            star_particles.stellar_returned_metals_channel_cumulative_code[channel][source];
        destination->stellar_feedback_energy_channel_cumulative_erg[channel][row] =
            star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][source];
      }
    }
    std::size_t row = kept_rows.size();
    for (std::size_t inbound_index = 0; inbound_index < commit.inbound_records.size(); ++inbound_index) {
      const auto& inbound = commit.inbound_records[inbound_index];
      if (!inbound.has_star_fields) {
        continue;
      }
      destination->particle_index[row] = checkedIntegralNarrow<std::uint32_t>(
          checkedSizeAdd(kept_count, inbound_index, "commitParticleMigration star particle index"),
          "commitParticleMigration star particle index");
      destination->formation_scale_factor[row] = inbound.star_fields.formation_scale_factor;
      destination->birth_mass_code[row] = inbound.star_fields.birth_mass_code;
      destination->metallicity_mass_fraction[row] = inbound.star_fields.metallicity_mass_fraction;
      destination->birth_key[row] = inbound.star_fields.birth_key;
      destination->parent_gas_cell_id[row] = inbound.star_fields.parent_gas_cell_id;
      destination->birth_tick[row] = inbound.star_fields.birth_tick;
      destination->birth_ordinal[row] = inbound.star_fields.birth_ordinal;
      destination->stellar_age_years_last[row] = inbound.star_fields.stellar_age_years_last;
      destination->stellar_returned_mass_cumulative_code[row] =
          inbound.star_fields.stellar_returned_mass_cumulative_code;
      destination->stellar_returned_metals_cumulative_code[row] =
          inbound.star_fields.stellar_returned_metals_cumulative_code;
      destination->stellar_newly_synthesized_metals_cumulative_code[row] =
          inbound.star_fields.stellar_newly_synthesized_metals_cumulative_code;
      destination->stellar_feedback_energy_cumulative_erg[row] =
          inbound.star_fields.stellar_feedback_energy_cumulative_erg;
      destination->enrichment_carry_mass_code[row] = inbound.star_fields.enrichment_carry_mass_code;
      destination->enrichment_carry_metals_code[row] = inbound.star_fields.enrichment_carry_metals_code;
      destination->enrichment_carry_feedback_energy_erg[row] =
          inbound.star_fields.enrichment_carry_feedback_energy_erg;
      destination->enrichment_carry_momentum_code[row] =
          inbound.star_fields.enrichment_carry_momentum_code;
      destination->stellar_deposited_mass_cumulative_code[row] =
          inbound.star_fields.stellar_deposited_mass_cumulative_code;
      destination->stellar_deposited_metals_cumulative_code[row] =
          inbound.star_fields.stellar_deposited_metals_cumulative_code;
      destination->stellar_deposited_feedback_energy_cumulative_erg[row] =
          inbound.star_fields.stellar_deposited_feedback_energy_cumulative_erg;
      for (std::size_t channel = 0; channel < inbound.star_fields.stellar_returned_mass_channel_cumulative_code.size();
           ++channel) {
        destination->stellar_returned_mass_channel_cumulative_code[channel][row] =
            inbound.star_fields.stellar_returned_mass_channel_cumulative_code[channel];
        destination->stellar_returned_metals_channel_cumulative_code[channel][row] =
            inbound.star_fields.stellar_returned_metals_channel_cumulative_code[channel];
        destination->stellar_feedback_energy_channel_cumulative_erg[channel][row] =
            inbound.star_fields.stellar_feedback_energy_channel_cumulative_erg[channel];
      }
      ++row;
    }
  };

  auto rebuildBlackHoleSidecar = [&](BlackHoleParticleSidecar* destination) {
    std::vector<std::uint32_t> kept_rows;
    for (std::size_t row = 0; row < black_holes.size(); ++row) {
      const auto old_particle = black_holes.particle_index[row];
      if (old_particle >= remove_mask.size()) {
        throw std::runtime_error("commitParticleMigration: stale black-hole sidecar index");
      }
      if (remove_mask[old_particle] == 0U) {
        kept_rows.push_back(checkedLocalParticleRow(
            row, "commitParticleMigration sidecar row"));
      }
    }
    const std::size_t inbound_count = static_cast<std::size_t>(std::count_if(
        commit.inbound_records.begin(), commit.inbound_records.end(),
        [](const ParticleMigrationRecord& record) { return record.has_black_hole_fields; }));
    destination->resize(checkedSizeAdd(kept_rows.size(), inbound_count, "commitParticleMigration black-hole rows"));
    for (std::size_t row = 0; row < kept_rows.size(); ++row) {
      const std::size_t source = kept_rows[row];
      const auto old_particle = black_holes.particle_index[source];
      destination->particle_index[row] = old_to_new[old_particle];
      destination->host_cell_index[row] = black_holes.host_cell_index[source];
      destination->subgrid_mass_code[row] = black_holes.subgrid_mass_code[source];
      destination->accretion_rate_code[row] = black_holes.accretion_rate_code[source];
      destination->feedback_energy_code[row] = black_holes.feedback_energy_code[source];
      destination->eddington_ratio[row] = black_holes.eddington_ratio[source];
      destination->cumulative_accreted_mass_code[row] = black_holes.cumulative_accreted_mass_code[source];
      destination->cumulative_feedback_energy_code[row] = black_holes.cumulative_feedback_energy_code[source];
      destination->duty_cycle_active_time_code[row] = black_holes.duty_cycle_active_time_code[source];
      destination->duty_cycle_total_time_code[row] = black_holes.duty_cycle_total_time_code[source];
    }
    std::size_t row = kept_rows.size();
    for (std::size_t inbound_index = 0; inbound_index < commit.inbound_records.size(); ++inbound_index) {
      const auto& inbound = commit.inbound_records[inbound_index];
      if (!inbound.has_black_hole_fields) {
        continue;
      }
      destination->particle_index[row] = checkedIntegralNarrow<std::uint32_t>(
          checkedSizeAdd(kept_count, inbound_index, "commitParticleMigration black-hole particle index"),
          "commitParticleMigration black-hole particle index");
      destination->host_cell_index[row] = inbound.black_hole_fields.host_cell_index;
      destination->subgrid_mass_code[row] = inbound.black_hole_fields.subgrid_mass_code;
      destination->accretion_rate_code[row] = inbound.black_hole_fields.accretion_rate_code;
      destination->feedback_energy_code[row] = inbound.black_hole_fields.feedback_energy_code;
      destination->eddington_ratio[row] = inbound.black_hole_fields.eddington_ratio;
      destination->cumulative_accreted_mass_code[row] = inbound.black_hole_fields.cumulative_accreted_mass_code;
      destination->cumulative_feedback_energy_code[row] =
          inbound.black_hole_fields.cumulative_feedback_energy_code;
      destination->duty_cycle_active_time_code[row] = inbound.black_hole_fields.duty_cycle_active_time_code;
      destination->duty_cycle_total_time_code[row] = inbound.black_hole_fields.duty_cycle_total_time_code;
      ++row;
    }
  };

  auto rebuildTracerSidecar = [&](TracerParticleSidecar* destination) {
    std::vector<std::uint32_t> kept_rows;
    for (std::size_t row = 0; row < tracers.size(); ++row) {
      const auto old_particle = tracers.particle_index[row];
      if (old_particle >= remove_mask.size()) {
        throw std::runtime_error("commitParticleMigration: stale tracer sidecar index");
      }
      if (remove_mask[old_particle] == 0U) {
        kept_rows.push_back(checkedLocalParticleRow(
            row, "commitParticleMigration sidecar row"));
      }
    }
    const std::size_t inbound_count = static_cast<std::size_t>(std::count_if(
        commit.inbound_records.begin(), commit.inbound_records.end(),
        [](const ParticleMigrationRecord& record) { return record.has_tracer_fields; }));
    destination->resize(checkedSizeAdd(kept_rows.size(), inbound_count, "commitParticleMigration tracer rows"));
    for (std::size_t row = 0; row < kept_rows.size(); ++row) {
      const std::size_t source = kept_rows[row];
      const auto old_particle = tracers.particle_index[source];
      destination->particle_index[row] = old_to_new[old_particle];
      destination->parent_particle_id[row] = tracers.parent_particle_id[source];
      destination->injection_step[row] = tracers.injection_step[source];
      destination->host_cell_index[row] = tracers.host_cell_index[source];
      destination->mass_fraction_of_host[row] = tracers.mass_fraction_of_host[source];
      destination->last_host_mass_code[row] = tracers.last_host_mass_code[source];
      destination->cumulative_exchanged_mass_code[row] = tracers.cumulative_exchanged_mass_code[source];
    }
    std::size_t row = kept_rows.size();
    for (std::size_t inbound_index = 0; inbound_index < commit.inbound_records.size(); ++inbound_index) {
      const auto& inbound = commit.inbound_records[inbound_index];
      if (!inbound.has_tracer_fields) {
        continue;
      }
      destination->particle_index[row] = checkedIntegralNarrow<std::uint32_t>(
          checkedSizeAdd(kept_count, inbound_index, "commitParticleMigration tracer particle index"),
          "commitParticleMigration tracer particle index");
      destination->parent_particle_id[row] = inbound.tracer_fields.parent_particle_id;
      destination->injection_step[row] = inbound.tracer_fields.injection_step;
      destination->host_cell_index[row] = inbound.tracer_fields.host_cell_index;
      destination->mass_fraction_of_host[row] = inbound.tracer_fields.mass_fraction_of_host;
      destination->last_host_mass_code[row] = inbound.tracer_fields.last_host_mass_code;
      destination->cumulative_exchanged_mass_code[row] = inbound.tracer_fields.cumulative_exchanged_mass_code;
      ++row;
    }
  };

  StarParticleSidecar rebuilt_star;
  BlackHoleParticleSidecar rebuilt_black_holes;
  TracerParticleSidecar rebuilt_tracers;
  rebuildStarSidecar(&rebuilt_star);
  rebuildBlackHoleSidecar(&rebuilt_black_holes);
  rebuildTracerSidecar(&rebuilt_tracers);
  if (rebuild_gas_state) {
    const auto remap_host_cell = [&](std::uint32_t old_cell_index) {
      if (old_cell_index == kInvalidGasCellRow) {
        return kInvalidGasCellRow;
      }
      if (old_cell_index >= old_cell_to_new.size()) {
        throw std::runtime_error("commitParticleMigration: sidecar host_cell_index is outside old CellSoa");
      }
      const std::uint32_t new_cell_index = old_cell_to_new[old_cell_index];
      if (new_cell_index == kInvalidGasCellRow) {
        throw std::runtime_error("commitParticleMigration: sidecar host gas cell was removed during gas migration");
      }
      return new_cell_index;
    };
    for (std::size_t row = 0; row < rebuilt_black_holes.size(); ++row) {
      rebuilt_black_holes.host_cell_index[row] = remap_host_cell(rebuilt_black_holes.host_cell_index[row]);
    }
    for (std::size_t row = 0; row < rebuilt_tracers.size(); ++row) {
      rebuilt_tracers.host_cell_index[row] = remap_host_cell(rebuilt_tracers.host_cell_index[row]);
    }
  }

  auto module_requires_final_particle = [&](const ModuleSidecarBlock& block, std::size_t final_particle) {
    ParticleMigrationRecord synthetic;
    synthetic.particle_id = new_sidecar.particle_id[final_particle];
    synthetic.species_tag = new_sidecar.species_tag[final_particle];
    synthetic.particle_flags = new_sidecar.particle_flags[final_particle];
    if (synthetic.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kGas)) {
      for (std::size_t row = 0; row < rebuilt_gas_cells.size(); ++row) {
        if (rebuilt_gas_cells.parent_particle_id[row] == synthetic.particle_id) {
          synthetic.has_gas_cell_fields = true;
          synthetic.gas_cell_fields.parent_particle_id = synthetic.particle_id;
          synthetic.gas_cell_fields.gas_cell_id = rebuilt_gas_cells.gas_cell_id[row];
          synthetic.gas_cell_fields.density_code = rebuilt_gas_cells.density_code[row];
          break;
        }
      }
    } else if (synthetic.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
      for (std::size_t row = 0; row < rebuilt_black_holes.size(); ++row) {
        if (rebuilt_black_holes.particle_index[row] == final_particle) {
          synthetic.has_black_hole_fields = true;
          synthetic.black_hole_fields.accretion_rate_code = rebuilt_black_holes.accretion_rate_code[row];
          break;
        }
      }
    }
    return moduleRequiresRecord(block, synthetic, "commitParticleMigration");
  };

  ModuleSidecarRegistry rebuilt_module_sidecars;
  std::unordered_set<std::string> rebuilt_module_names;
  auto append_particle_indexed_module_rows = [&](const ModuleSidecarBlock& source_block, ModuleSidecarBlock* destination) {
    destination->module_name = source_block.module_name;
    destination->schema_version = source_block.schema_version;
    destination->particle_indexed = true;
    destination->row_stride_bytes = source_block.row_stride_bytes;
    destination->required_species_mask = source_block.required_species_mask;
    destination->requirement = source_block.requirement;
    destination->payload.clear();
    destination->particle_id_by_row.clear();

    const auto inbound_module_it = inbound_module_payloads_by_module.find(source_block.module_name);
    for (std::size_t final_particle = 0; final_particle < final_count; ++final_particle) {
      const std::uint64_t particle_id = new_sidecar.particle_id[final_particle];
      if (inbound_module_it != inbound_module_payloads_by_module.end()) {
        const auto payload_it = inbound_module_it->second.find(particle_id);
        if (payload_it != inbound_module_it->second.end()) {
          appendModuleRowPayload(destination, particle_id, payload_it->second->payload);
          continue;
        }
      }
      const auto source_row = findModuleSidecarRowByParticleId(source_block, particle_id);
      if (source_row.has_value()) {
        appendModuleRowPayload(destination, particle_id, source_block.rowPayload(*source_row));
        continue;
      }
      if (module_requires_final_particle(source_block, final_particle)) {
        throw std::invalid_argument(
            "commitParticleMigration: required module sidecar row is missing after particle migration");
      }
    }
  };

  for (const ModuleSidecarBlock* block : sidecars.blocksSortedByName()) {
    if (!block->isParticleIndexed()) {
      rebuilt_module_sidecars.upsert(*block);
      rebuilt_module_names.insert(block->module_name);
      continue;
    }
    ModuleSidecarBlock rebuilt_block;
    append_particle_indexed_module_rows(*block, &rebuilt_block);
    rebuilt_module_sidecars.upsert(std::move(rebuilt_block));
    rebuilt_module_names.insert(block->module_name);
  }

  for (const auto& [module_name, payloads_by_particle] : inbound_module_payloads_by_module) {
    if (rebuilt_module_names.contains(module_name)) {
      continue;
    }
    if (payloads_by_particle.empty()) {
      continue;
    }
    const ModuleParticleSidecarPayload* prototype = payloads_by_particle.begin()->second;
    ModuleSidecarBlock rebuilt_block;
    rebuilt_block.module_name = module_name;
    rebuilt_block.schema_version = prototype->schema_version;
    rebuilt_block.particle_indexed = true;
    rebuilt_block.row_stride_bytes = prototype->row_stride_bytes;
    rebuilt_block.required_species_mask = prototype->required_species_mask;
    rebuilt_block.requirement = prototype->requirement;
    for (const auto& [particle_id, payload] : payloads_by_particle) {
      if (payload->schema_version != rebuilt_block.schema_version ||
          payload->row_stride_bytes != rebuilt_block.row_stride_bytes ||
          payload->required_species_mask != rebuilt_block.required_species_mask ||
          payload->requirement.kind != rebuilt_block.requirement.kind ||
          payload->requirement.species_mask != rebuilt_block.requirement.species_mask ||
          payload->requirement.particle_flags_mask != rebuilt_block.requirement.particle_flags_mask ||
          payload->requirement.threshold_code != rebuilt_block.requirement.threshold_code) {
        throw std::invalid_argument("commitParticleMigration: inbound module sidecar payloads disagree on schema, stride, or required coverage");
      }
    }
    for (std::size_t final_particle = 0; final_particle < final_count; ++final_particle) {
      const std::uint64_t particle_id = new_sidecar.particle_id[final_particle];
      const auto payload_it = payloads_by_particle.find(particle_id);
      if (payload_it != payloads_by_particle.end()) {
        appendModuleRowPayload(&rebuilt_block, particle_id, payload_it->second->payload);
        continue;
      }
      if (module_requires_final_particle(rebuilt_block, final_particle)) {
        throw std::invalid_argument(
            "commitParticleMigration: required inbound module sidecar row is missing after particle migration");
      }
    }
    rebuilt_module_sidecars.upsert(std::move(rebuilt_block));
    rebuilt_module_names.insert(module_name);
  }

  particles = std::move(new_particles);
  particle_sidecar = std::move(new_sidecar);
  if (rebuild_gas_state) {
    cells = std::move(rebuilt_cells);
    gas_cells = std::move(rebuilt_gas_cells);
    gas_composition = std::move(rebuilt_composition);
    replaceGasCellIdentityRecords(gasIdentityRecordsFromSidecars(cells, gas_cells, patches));
  }
  star_particles = std::move(rebuilt_star);
  black_holes = std::move(rebuilt_black_holes);
  tracers = std::move(rebuilt_tracers);
  sidecars = std::move(rebuilt_module_sidecars);

  species.count_by_species.fill(0);
  for (const auto tag : particle_sidecar.species_tag) {
    if (!isValidParticleSpeciesTag(tag)) {
      throw std::invalid_argument("commitParticleMigration: invalid species tag in rebuilt sidecar");
    }
    ++species.count_by_species[tag];
  }
  rebuildSpeciesIndex();
  bumpParticleIdentityGeneration();
  m_particle_identity_validation_certificate = ParticleIdentityValidationCertificate{
      .identity_generation = m_particle_identity_generation,
      .particle_count = particle_sidecar.particle_id.size(),
      .local_uniqueness_proven = true,
      .nonzero_ids_proven = true,
  };
  bumpParticleIndexGeneration();
}

std::vector<GasCellMigrationRecord> SimulationState::packGasCellMigrationRecords(
    std::span<const std::uint32_t> local_cell_rows,
    std::uint64_t ghost_hydro_epoch) const {
  this->requireGasCellIdentityMapCoversDenseRows("packGasCellMigrationRecords");
  std::vector<GasCellMigrationRecord> records;
  records.reserve(local_cell_rows.size());
  for (const std::uint32_t row : local_cell_rows) {
    const GasCellMigrationFields fields =
        gasCellFieldsFromLocalRow(*this, row, ghost_hydro_epoch, "packGasCellMigrationRecords");
    std::uint32_t owner_rank = 0U;
    if (fields.patch_index < patches.owning_rank.size()) {
      owner_rank = patches.owning_rank[fields.patch_index];
    }
    records.push_back(GasCellMigrationRecord{
        .owning_rank = owner_rank,
        .fields = fields,
    });
  }
  return records;
}

void SimulationState::commitGasCellMigration(const GasCellMigrationCommit& commit) {
  if (commit.world_rank < 0) {
    throw std::invalid_argument("commitGasCellMigration: world_rank must be non-negative");
  }
  this->requireGasCellIdentityMapCoversDenseRows("commitGasCellMigration");

  std::vector<std::uint8_t> remove_mask(cells.size(), 0U);
  const auto mark_remove = [&](std::uint32_t row, std::string_view label) {
    if (row >= cells.size()) {
      throw std::out_of_range(std::string(label) + ": local gas-cell row is out of range");
    }
    if (remove_mask[row] != 0U) {
      throw std::invalid_argument(std::string(label) + ": duplicate gas-cell row removal");
    }
    remove_mask[row] = 1U;
  };

  for (const std::uint32_t row : commit.outbound_local_cell_rows) {
    mark_remove(row, "commitGasCellMigration outbound");
  }
  for (const GasCellStaleGhostRecord& stale : commit.stale_local_ghost_records) {
    if (stale.gas_cell_identity_generation != commit.expected_gas_cell_identity_generation ||
        stale.ghost_hydro_epoch != commit.expected_ghost_hydro_epoch) {
      throw std::invalid_argument("commitGasCellMigration: stale ghost generation or hydro epoch mismatch");
    }
    if (stale.local_cell_row >= cells.size()) {
      throw std::out_of_range("commitGasCellMigration: stale ghost local row is out of range");
    }
    const auto* identity = gas_cell_identity.findByLocalRow(stale.local_cell_row);
    if (identity == nullptr || identity->gas_cell_id != stale.gas_cell_id) {
      throw std::invalid_argument("commitGasCellMigration: stale ghost gas_cell_id does not match local row");
    }
    mark_remove(stale.local_cell_row, "commitGasCellMigration stale ghost");
  }

  std::size_t removed_count = 0U;
  for (const std::uint8_t removed : remove_mask) {
    removed_count += removed != 0U ? 1U : 0U;
  }
  const std::size_t kept_count = cells.size() - removed_count;
  const std::size_t final_count = checkedLocalCountAdd(
      kept_count,
      commit.inbound_records.size(),
      kMaxLocalCellCount,
      "gas-cell",
      "commitGasCellMigration final cell count");

  std::unordered_set<std::uint64_t> final_gas_cell_ids;
  final_gas_cell_ids.reserve(final_count);
  std::vector<std::uint32_t> old_cell_to_new(cells.size(), kInvalidGasCellRow);
  CellSoa rebuilt_cells;
  GasCellSidecar rebuilt_gas_cells;
  rebuilt_cells.resize(final_count);
  rebuilt_gas_cells.resize(final_count);
  GasCompositionSidecar rebuilt_composition = gas_composition.emptyLike(final_count);
  std::size_t active_sparse_rows = 0U;
  if (gas_composition.sparseStorageEnabled()) {
    for (std::size_t row = 0; row < cells.size(); ++row) {
      if (remove_mask[row] == 0U && gas_composition.hasExplicitSparseRow(row)) {
        ++active_sparse_rows;
      }
    }
    for (const auto& inbound : commit.inbound_records) {
      active_sparse_rows += migrationPayloadHasExplicitSparseRow(
          inbound.fields, rebuilt_composition) ? 1U : 0U;
    }
  }
  prepareAdmittedSparseCompositionCandidate(
      rebuilt_composition, active_sparse_rows,
      "gas_cell_migration.adaptive_chemistry_rows");

  std::size_t write_row = 0U;
  for (std::uint32_t row = 0; row < cells.size(); ++row) {
    if (remove_mask[row] != 0U) {
      continue;
    }
    GasCellMigrationFields fields = gasCellFieldsFromLocalRow(*this, row, 0U, "commitGasCellMigration");
    const std::uint32_t destination_row = checkedLocalCellRow(
        write_row, "commitGasCellMigration kept destination row");
    fields.destination_local_cell_row = destination_row;
    old_cell_to_new[row] = destination_row;
    if (!final_gas_cell_ids.insert(fields.gas_cell_id).second) {
      throw std::invalid_argument("commitGasCellMigration: duplicate kept gas_cell_id");
    }
    writeGasCellFieldsToRow(fields, destination_row, rebuilt_cells, rebuilt_gas_cells, rebuilt_composition);
    ++write_row;
  }
  for (const GasCellMigrationRecord& inbound : commit.inbound_records) {
    if (inbound.owning_rank != static_cast<std::uint32_t>(commit.world_rank)) {
      throw std::invalid_argument("commitGasCellMigration: inbound gas-cell owner rank does not match commit rank");
    }
    requireValidGasCellMigrationFields(inbound.fields, "commitGasCellMigration");
    GasCellMigrationFields fields = inbound.fields;
    const std::uint32_t destination_row = checkedLocalCellRow(
        write_row, "commitGasCellMigration inbound destination row");
    fields.destination_local_cell_row = destination_row;
    if (!final_gas_cell_ids.insert(fields.gas_cell_id).second) {
      throw std::invalid_argument("commitGasCellMigration: inbound gas_cell_id duplicates local gas-cell state");
    }
    writeGasCellFieldsToRow(fields, destination_row, rebuilt_cells, rebuilt_gas_cells, rebuilt_composition);
    ++write_row;
  }
  if (write_row != final_count) {
    throw std::runtime_error("commitGasCellMigration: prepared cell count does not match preflight final count");
  }

  const auto remap_host_cell = [&](std::uint32_t old_cell_index) {
    if (old_cell_index == kInvalidGasCellRow) {
      return kInvalidGasCellRow;
    }
    if (old_cell_index >= old_cell_to_new.size()) {
      throw std::runtime_error("commitGasCellMigration: sidecar host_cell_index is outside old CellSoa");
    }
    const std::uint32_t new_cell_index = old_cell_to_new[old_cell_index];
    if (new_cell_index == kInvalidGasCellRow) {
      throw std::runtime_error("commitGasCellMigration: sidecar host gas cell was removed during migration");
    }
    return new_cell_index;
  };
  auto rebuilt_black_hole_hosts = black_holes.host_cell_index;
  auto rebuilt_tracer_hosts = tracers.host_cell_index;
  for (std::size_t row = 0; row < rebuilt_black_hole_hosts.size(); ++row) {
    rebuilt_black_hole_hosts[row] = remap_host_cell(rebuilt_black_hole_hosts[row]);
  }
  for (std::size_t row = 0; row < rebuilt_tracer_hosts.size(); ++row) {
    rebuilt_tracer_hosts[row] = remap_host_cell(rebuilt_tracer_hosts[row]);
  }

  auto rebuilt_identity_records = gasIdentityRecordsFromSidecars(rebuilt_cells, rebuilt_gas_cells, patches);
  if (gas_cell_identity.generation() == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("commitGasCellMigration: gas-cell identity generation overflow");
  }
  GasCellIdentityMap rebuilt_identity;
  rebuilt_identity.assignWithGeneration(
      std::move(rebuilt_identity_records), gas_cell_identity.generation() + 1U);

  cells = std::move(rebuilt_cells);
  gas_cells = std::move(rebuilt_gas_cells);
  gas_composition = std::move(rebuilt_composition);
  black_holes.host_cell_index = std::move(rebuilt_black_hole_hosts);
  tracers.host_cell_index = std::move(rebuilt_tracer_hosts);
  gas_cell_identity = std::move(rebuilt_identity);
  bumpCellIndexGeneration();
}

std::vector<AmrPatchMigrationRecord> SimulationState::packAmrPatchMigrationRecords(
    std::span<const std::uint32_t> local_patch_indices,
    std::uint64_t ghost_hydro_epoch) const {
  if (!patches.isConsistent()) {
    throw std::runtime_error("packAmrPatchMigrationRecords: PatchSoa is inconsistent");
  }
  this->requireGasCellIdentityMapCoversDenseRows("packAmrPatchMigrationRecords");

  std::vector<AmrPatchMigrationRecord> records;
  records.reserve(local_patch_indices.size());
  for (const std::uint32_t patch_index : local_patch_indices) {
    if (patch_index >= patches.size()) {
      throw std::out_of_range("packAmrPatchMigrationRecords: local patch index is out of range");
    }
    const std::uint32_t first_cell = patches.first_cell[patch_index];
    const std::uint32_t cell_count = patches.cell_count[patch_index];
    if (static_cast<std::uint64_t>(first_cell) + static_cast<std::uint64_t>(cell_count) > cells.size()) {
      throw std::runtime_error("packAmrPatchMigrationRecords: patch cell range is outside CellSoa");
    }

    AmrPatchMigrationRecord record;
    record.patch = patchFieldsFromLocalRow(*this, patch_index, "packAmrPatchMigrationRecords");
    record.gas_cell_records.reserve(cell_count);
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      const std::uint32_t cell_row = first_cell + offset;
      GasCellMigrationFields fields =
          gasCellFieldsFromLocalRow(*this, cell_row, ghost_hydro_epoch, "packAmrPatchMigrationRecords");
      fields.patch_index = 0U;
      fields.owning_patch_id = record.patch.patch_id;
      fields.destination_local_cell_row = offset;
      record.gas_cell_records.push_back(GasCellMigrationRecord{
          .owning_rank = record.patch.owning_rank,
          .fields = fields,
      });
    }
    if (record.gas_cell_records.size() != record.patch.cell_count) {
      throw std::runtime_error("packAmrPatchMigrationRecords: patch payload did not cover every gas-cell row");
    }
    for (const PendingFluxRegisterRecord& pending : pending_flux_registers.records()) {
      if (pending.coarse_patch_id == record.patch.patch_id) {
        record.pending_flux_register_records.push_back(pending);
      }
    }
    for (const AmrTemporalBoundaryHistoryRecord& history :
         amr_temporal_boundary_history.records()) {
      if (history.patch_id == record.patch.patch_id) {
        record.temporal_boundary_history_records.push_back(history);
      }
    }
    records.push_back(std::move(record));
  }
  return records;
}

void SimulationState::commitAmrPatchMigration(const AmrPatchMigrationCommit& commit) {
  if (commit.world_rank < 0) {
    throw std::invalid_argument("commitAmrPatchMigration: world_rank must be non-negative");
  }
  if (!patches.isConsistent()) {
    throw std::runtime_error("commitAmrPatchMigration: PatchSoa is inconsistent");
  }
  this->requireGasCellIdentityMapCoversDenseRows("commitAmrPatchMigration");

  std::vector<std::uint8_t> remove_patch(patches.size(), 0U);
  for (const std::uint32_t patch_index : commit.outbound_local_patch_indices) {
    if (patch_index >= patches.size()) {
      throw std::out_of_range("commitAmrPatchMigration: outbound patch index is out of range");
    }
    if (remove_patch[patch_index] != 0U) {
      throw std::invalid_argument("commitAmrPatchMigration: duplicate outbound patch index");
    }
    if (patches.owning_rank[patch_index] != static_cast<std::uint32_t>(commit.world_rank)) {
      throw std::invalid_argument("commitAmrPatchMigration: outbound patch is not owned by committing rank");
    }
    remove_patch[patch_index] = 1U;
  }
  for (const std::uint32_t patch_index : commit.stale_local_ghost_patch_indices) {
    if (patch_index >= patches.size()) {
      throw std::out_of_range("commitAmrPatchMigration: stale ghost patch index is out of range");
    }
    if (remove_patch[patch_index] != 0U) {
      throw std::invalid_argument("commitAmrPatchMigration: duplicate patch removal");
    }
    if (patches.owning_rank[patch_index] == static_cast<std::uint32_t>(commit.world_rank)) {
      throw std::invalid_argument("commitAmrPatchMigration: stale ghost patch is owned by committing rank");
    }
    remove_patch[patch_index] = 1U;
  }

  std::vector<std::uint8_t> remove_cell(cells.size(), 0U);
  for (std::uint32_t cell_row = 0; cell_row < cells.size(); ++cell_row) {
    const std::uint32_t patch_index = cells.patch_index[cell_row];
    if (patch_index >= patches.size()) {
      throw std::runtime_error("commitAmrPatchMigration: cell patch_index is outside PatchSoa");
    }
    if (remove_patch[patch_index] != 0U) {
      remove_cell[cell_row] = 1U;
    }
  }

  for (const GasCellStaleGhostRecord& stale : commit.stale_local_ghost_records) {
    if (stale.gas_cell_identity_generation != commit.expected_gas_cell_identity_generation ||
        stale.ghost_hydro_epoch != commit.expected_ghost_hydro_epoch) {
      throw std::invalid_argument("commitAmrPatchMigration: stale gas ghost generation or hydro epoch mismatch");
    }
    if (stale.local_cell_row >= cells.size()) {
      throw std::out_of_range("commitAmrPatchMigration: stale gas ghost local row is out of range");
    }
    const auto* identity = gas_cell_identity.findByLocalRow(stale.local_cell_row);
    if (identity == nullptr || identity->gas_cell_id != stale.gas_cell_id) {
      throw std::invalid_argument("commitAmrPatchMigration: stale gas ghost gas_cell_id does not match local row");
    }
    if (remove_cell[stale.local_cell_row] != 0U) {
      throw std::invalid_argument("commitAmrPatchMigration: duplicate gas-cell removal");
    }
    remove_cell[stale.local_cell_row] = 1U;
  }

  std::unordered_set<std::uint64_t> final_patch_ids;
  final_patch_ids.reserve(checkedSizeAdd(
      patches.size(), commit.inbound_records.size(), "commitAmrPatchMigration patch-id reserve"));
  for (std::uint32_t patch_index = 0; patch_index < patches.size(); ++patch_index) {
    if (remove_patch[patch_index] != 0U) {
      continue;
    }
    if (!final_patch_ids.insert(patches.patch_id[patch_index]).second) {
      throw std::invalid_argument("commitAmrPatchMigration: kept patches contain duplicate patch_id");
    }
  }
  const std::unordered_set<std::uint64_t> kept_original_patch_ids = final_patch_ids;
  for (const AmrPatchMigrationRecord& inbound : commit.inbound_records) {
    requireValidAmrPatchMigrationFields(inbound.patch, "commitAmrPatchMigration");
    if (inbound.patch.owning_rank != static_cast<std::uint32_t>(commit.world_rank)) {
      throw std::invalid_argument("commitAmrPatchMigration: inbound patch owner rank does not match commit rank");
    }
    if (inbound.gas_cell_records.size() != inbound.patch.cell_count) {
      throw std::invalid_argument("commitAmrPatchMigration: inbound patch gas-cell payload count mismatch");
    }
    std::unordered_set<std::uint64_t> inbound_gas_cell_ids;
    inbound_gas_cell_ids.reserve(inbound.gas_cell_records.size());
    for (const GasCellMigrationRecord& gas_record : inbound.gas_cell_records) {
      if (!inbound_gas_cell_ids.insert(gas_record.fields.gas_cell_id).second) {
        throw std::invalid_argument("commitAmrPatchMigration: inbound patch contains duplicate gas_cell_id");
      }
    }
    std::unordered_set<std::uint64_t> inbound_flux_keys;
    inbound_flux_keys.reserve(inbound.pending_flux_register_records.size());
    for (const PendingFluxRegisterRecord& pending : inbound.pending_flux_register_records) {
      if (pending.register_key == 0U || pending.coarse_patch_id != inbound.patch.patch_id ||
          !inbound_gas_cell_ids.contains(pending.coarse_gas_cell_id)) {
        throw std::invalid_argument(
            "commitAmrPatchMigration: inbound pending flux register does not belong to the migrated patch");
      }
      if (!inbound_flux_keys.insert(pending.register_key).second) {
        throw std::invalid_argument("commitAmrPatchMigration: inbound patch contains duplicate flux-register key");
      }
    }
    if (inbound.temporal_boundary_history_records.size() > 1U) {
      throw std::invalid_argument(
          "commitAmrPatchMigration: inbound patch contains multiple temporal-boundary histories");
    }
    for (const AmrTemporalBoundaryHistoryRecord& history : inbound.temporal_boundary_history_records) {
      if (history.patch_id != inbound.patch.patch_id || history.cells.size() != inbound.patch.cell_count) {
        throw std::invalid_argument(
            "commitAmrPatchMigration: inbound temporal-boundary history does not cover the migrated patch");
      }
      std::vector<std::uint8_t> seen_local_cell(history.cells.size(), 0U);
      for (const AmrTemporalBoundaryHistoryCellRecord& history_cell : history.cells) {
        if (history_cell.patch_local_cell >= inbound.patch.cell_count ||
            !inbound_gas_cell_ids.contains(history_cell.gas_cell_id)) {
          throw std::invalid_argument(
              "commitAmrPatchMigration: inbound temporal-boundary history references a foreign gas cell");
        }
        if (seen_local_cell[history_cell.patch_local_cell] != 0U) {
          throw std::invalid_argument(
              "commitAmrPatchMigration: inbound temporal-boundary history duplicates a patch-local cell");
        }
        seen_local_cell[history_cell.patch_local_cell] = 1U;
      }
    }
    if (!final_patch_ids.insert(inbound.patch.patch_id).second) {
      throw std::invalid_argument("commitAmrPatchMigration: inbound patch_id duplicates local patch state");
    }
  }

  std::vector<AmrPatchMigrationFields> final_patches;
  std::vector<std::vector<GasCellMigrationFields>> final_cells_by_patch;
  final_patches.reserve(checkedSizeAdd(
      patches.size(), commit.inbound_records.size(), "commitAmrPatchMigration patch reserve"));
  final_cells_by_patch.reserve(checkedSizeAdd(
      patches.size(), commit.inbound_records.size(), "commitAmrPatchMigration cell-vector reserve"));
  std::unordered_set<std::uint64_t> final_gas_cell_ids;
  final_gas_cell_ids.reserve(cells.size());

  for (std::uint32_t patch_index = 0; patch_index < patches.size(); ++patch_index) {
    if (remove_patch[patch_index] != 0U) {
      continue;
    }
    AmrPatchMigrationFields patch = patchFieldsFromLocalRow(*this, patch_index, "commitAmrPatchMigration");
    std::vector<GasCellMigrationFields> patch_cells;
    const std::uint32_t first_cell = patches.first_cell[patch_index];
    const std::uint32_t cell_count = patches.cell_count[patch_index];
    if (static_cast<std::uint64_t>(first_cell) + static_cast<std::uint64_t>(cell_count) > cells.size()) {
      throw std::runtime_error("commitAmrPatchMigration: kept patch cell range is outside CellSoa");
    }
    patch_cells.reserve(cell_count);
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      const std::uint32_t old_cell = first_cell + offset;
      if (remove_cell[old_cell] != 0U) {
        continue;
      }
      GasCellMigrationFields fields =
          gasCellFieldsFromLocalRow(*this, old_cell, 0U, "commitAmrPatchMigration");
      fields.owning_patch_id = patch.patch_id;
      patch_cells.push_back(fields);
    }
    patch.cell_count = checkedLocalPatchCellCount(
        patch_cells.size(), "commitAmrPatchMigration kept patch cell count");
    if (patch.cell_dim_x != 0U || patch.cell_dim_y != 0U || patch.cell_dim_z != 0U) {
      const std::size_t dims_product = static_cast<std::size_t>(patch.cell_dim_x) *
          static_cast<std::size_t>(patch.cell_dim_y) * static_cast<std::size_t>(patch.cell_dim_z);
      if (dims_product != patch.cell_count) {
        throw std::runtime_error("commitAmrPatchMigration: kept AMR patch lost cells without a descriptor remap");
      }
    }
    final_patches.push_back(patch);
    final_cells_by_patch.push_back(std::move(patch_cells));
  }

  for (const AmrPatchMigrationRecord& inbound : commit.inbound_records) {
    std::vector<GasCellMigrationFields> patch_cells;
    patch_cells.reserve(inbound.gas_cell_records.size());
    for (const GasCellMigrationRecord& gas_record : inbound.gas_cell_records) {
      if (gas_record.owning_rank != static_cast<std::uint32_t>(commit.world_rank)) {
        throw std::invalid_argument("commitAmrPatchMigration: inbound gas-cell owner rank does not match commit rank");
      }
      requireValidGasCellMigrationFields(gas_record.fields, "commitAmrPatchMigration");
      if (gas_record.fields.owning_patch_id != inbound.patch.patch_id) {
        throw std::invalid_argument("commitAmrPatchMigration: inbound gas-cell owning_patch_id does not match patch");
      }
      patch_cells.push_back(gas_record.fields);
    }
    final_patches.push_back(inbound.patch);
    final_cells_by_patch.push_back(std::move(patch_cells));
  }

  std::size_t total_cell_count = 0;
  for (const auto& patch_cells : final_cells_by_patch) {
    total_cell_count = checkedLocalCountAdd(
        total_cell_count,
        patch_cells.size(),
        kMaxLocalCellCount,
        "gas-cell",
        "commitAmrPatchMigration total final cells");
  }
  (void)checkedLocalCount(final_patches.size(), kMaxLocalPatchCount, "AMR patch", "commitAmrPatchMigration final patches");
  PatchSoa rebuilt_patches;
  CellSoa rebuilt_cells;
  GasCellSidecar rebuilt_gas_cells;
  rebuilt_patches.resize(final_patches.size());
  rebuilt_cells.resize(total_cell_count);
  rebuilt_gas_cells.resize(total_cell_count);
  GasCompositionSidecar rebuilt_composition = gas_composition.emptyLike(total_cell_count);
  std::size_t active_sparse_rows = 0U;
  if (rebuilt_composition.sparseStorageEnabled()) {
    for (const auto& patch_cells : final_cells_by_patch) {
      for (const auto& fields : patch_cells) {
        active_sparse_rows += migrationPayloadHasExplicitSparseRow(
            fields, rebuilt_composition) ? 1U : 0U;
      }
    }
  }
  prepareAdmittedSparseCompositionCandidate(
      rebuilt_composition, active_sparse_rows,
      "amr_patch_migration.adaptive_chemistry_rows");

  std::uint32_t write_cell = 0;
  for (std::size_t patch_row = 0; patch_row < final_patches.size(); ++patch_row) {
    const std::uint32_t patch_index = checkedLocalPatchRow(
        patch_row, "commitAmrPatchMigration final patch row");
    AmrPatchMigrationFields patch = final_patches[patch_row];
    patch.cell_count = checkedLocalPatchCellCount(
        final_cells_by_patch[patch_row].size(), "commitAmrPatchMigration final patch cell count");
    writePatchFieldsToRow(patch, write_cell, patch.cell_count, patch_index, rebuilt_patches);
    for (GasCellMigrationFields fields : final_cells_by_patch[patch_row]) {
      fields.patch_index = patch_index;
      fields.destination_local_cell_row = write_cell;
      if (!final_gas_cell_ids.insert(fields.gas_cell_id).second) {
        throw std::invalid_argument("commitAmrPatchMigration: duplicate final gas_cell_id");
      }
      writeGasCellFieldsToRow(fields, write_cell, rebuilt_cells, rebuilt_gas_cells, rebuilt_composition);
      ++write_cell;
    }
  }

  std::unordered_map<std::uint64_t, std::uint32_t> new_row_by_gas_cell_id;
  new_row_by_gas_cell_id.reserve(rebuilt_gas_cells.size());
  for (std::uint32_t row = 0; row < rebuilt_gas_cells.size(); ++row) {
    new_row_by_gas_cell_id.emplace(rebuilt_gas_cells.gas_cell_id[row], row);
  }
  const auto remap_host_cell = [&](std::uint32_t old_cell_index) {
    if (old_cell_index == kInvalidGasCellRow) {
      return kInvalidGasCellRow;
    }
    if (old_cell_index >= cells.size()) {
      throw std::runtime_error("commitAmrPatchMigration: sidecar host_cell_index is outside old CellSoa");
    }
    const std::uint64_t gas_cell_id = gas_cells.gas_cell_id[old_cell_index];
    const auto row_it = new_row_by_gas_cell_id.find(gas_cell_id);
    if (row_it == new_row_by_gas_cell_id.end()) {
      throw std::runtime_error("commitAmrPatchMigration: sidecar host gas cell was removed during patch migration");
    }
    return row_it->second;
  };
  auto rebuilt_black_hole_hosts = black_holes.host_cell_index;
  auto rebuilt_tracer_hosts = tracers.host_cell_index;
  for (std::size_t row = 0; row < rebuilt_black_hole_hosts.size(); ++row) {
    rebuilt_black_hole_hosts[row] = remap_host_cell(rebuilt_black_hole_hosts[row]);
  }
  for (std::size_t row = 0; row < rebuilt_tracer_hosts.size(); ++row) {
    rebuilt_tracer_hosts[row] = remap_host_cell(rebuilt_tracer_hosts[row]);
  }

  auto rebuilt_identity_records =
      gasIdentityRecordsFromSidecars(rebuilt_cells, rebuilt_gas_cells, rebuilt_patches);
  if (gas_cell_identity.generation() == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("commitAmrPatchMigration: gas-cell identity generation overflow");
  }
  GasCellIdentityMap rebuilt_identity;
  rebuilt_identity.assignWithGeneration(
      std::move(rebuilt_identity_records), gas_cell_identity.generation() + 1U);
  if (cellIndexGeneration() == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("commitAmrPatchMigration: cell-index generation overflow");
  }
  const std::uint64_t rebuilt_cell_index_generation = cellIndexGeneration() + 1U;

  std::vector<PendingFluxRegisterRecord> rebuilt_pending_flux_records;
  std::size_t incoming_pending_count = 0U;
  std::size_t incoming_history_count = 0U;
  for (const AmrPatchMigrationRecord& inbound : commit.inbound_records) {
    incoming_pending_count = checkedSizeAdd(
        incoming_pending_count,
        inbound.pending_flux_register_records.size(),
        "commitAmrPatchMigration pending flux count");
    incoming_history_count = checkedSizeAdd(
        incoming_history_count,
        inbound.temporal_boundary_history_records.size(),
        "commitAmrPatchMigration temporal history count");
  }
  rebuilt_pending_flux_records.reserve(checkedSizeAdd(
      pending_flux_registers.size(), incoming_pending_count,
      "commitAmrPatchMigration pending flux reserve"));
  std::unordered_set<std::uint64_t> final_flux_keys;
  final_flux_keys.reserve(checkedSizeAdd(
      pending_flux_registers.size(), incoming_pending_count,
      "commitAmrPatchMigration flux-key reserve"));
  for (const PendingFluxRegisterRecord& pending : pending_flux_registers.records()) {
    if (!kept_original_patch_ids.contains(pending.coarse_patch_id)) {
      continue;
    }
    if (!final_flux_keys.insert(pending.register_key).second) {
      throw std::invalid_argument("commitAmrPatchMigration: kept pending flux registers contain duplicate key");
    }
    rebuilt_pending_flux_records.push_back(pending);
  }
  for (const AmrPatchMigrationRecord& inbound : commit.inbound_records) {
    for (const PendingFluxRegisterRecord& pending : inbound.pending_flux_register_records) {
      if (!final_flux_keys.insert(pending.register_key).second) {
        throw std::invalid_argument(
            "commitAmrPatchMigration: inbound pending flux-register key duplicates destination state");
      }
      rebuilt_pending_flux_records.push_back(pending);
    }
  }

  std::vector<AmrTemporalBoundaryHistoryRecord> rebuilt_temporal_history_records;
  rebuilt_temporal_history_records.reserve(checkedSizeAdd(
      amr_temporal_boundary_history.size(), incoming_history_count,
      "commitAmrPatchMigration temporal history reserve"));
  std::unordered_set<std::uint64_t> final_history_patch_ids;
  final_history_patch_ids.reserve(checkedSizeAdd(
      amr_temporal_boundary_history.size(), incoming_history_count,
      "commitAmrPatchMigration temporal history key reserve"));
  for (const AmrTemporalBoundaryHistoryRecord& history : amr_temporal_boundary_history.records()) {
    if (!kept_original_patch_ids.contains(history.patch_id)) {
      continue;
    }
    if (!final_history_patch_ids.insert(history.patch_id).second) {
      throw std::invalid_argument("commitAmrPatchMigration: kept temporal histories contain duplicate patch_id");
    }
    rebuilt_temporal_history_records.push_back(history);
  }
  for (const AmrPatchMigrationRecord& inbound : commit.inbound_records) {
    for (const AmrTemporalBoundaryHistoryRecord& history : inbound.temporal_boundary_history_records) {
      if (!final_history_patch_ids.insert(history.patch_id).second) {
        throw std::invalid_argument(
            "commitAmrPatchMigration: inbound temporal history duplicates destination patch state");
      }
      rebuilt_temporal_history_records.push_back(history);
    }
  }

  std::unordered_map<std::uint64_t, std::uint32_t> final_patch_row_by_id;
  final_patch_row_by_id.reserve(rebuilt_patches.size());
  for (std::uint32_t patch_row = 0; patch_row < rebuilt_patches.size(); ++patch_row) {
    final_patch_row_by_id.emplace(rebuilt_patches.patch_id[patch_row], patch_row);
  }
  for (PendingFluxRegisterRecord& pending : rebuilt_pending_flux_records) {
    const auto patch_it = final_patch_row_by_id.find(pending.coarse_patch_id);
    const auto gas_it = new_row_by_gas_cell_id.find(pending.coarse_gas_cell_id);
    if (patch_it == final_patch_row_by_id.end() || gas_it == new_row_by_gas_cell_id.end()) {
      throw std::invalid_argument(
          "commitAmrPatchMigration: pending flux register lost its stable patch/gas target");
    }
    const std::uint32_t gas_row = gas_it->second;
    if (rebuilt_cells.patch_index[gas_row] != patch_it->second) {
      throw std::invalid_argument(
          "commitAmrPatchMigration: pending flux register gas target belongs to a different patch");
    }
    pending.coarse_cell_index = gas_row;
    pending.gas_cell_identity_generation = rebuilt_identity.generation();
    pending.patch_geometry_generation = rebuilt_cell_index_generation;
  }
  for (AmrTemporalBoundaryHistoryRecord& history : rebuilt_temporal_history_records) {
    const auto patch_it = final_patch_row_by_id.find(history.patch_id);
    if (patch_it == final_patch_row_by_id.end()) {
      throw std::invalid_argument("commitAmrPatchMigration: temporal history lost its patch target");
    }
    const std::uint32_t patch_row = patch_it->second;
    if (history.cells.size() != rebuilt_patches.cell_count[patch_row]) {
      throw std::invalid_argument(
          "commitAmrPatchMigration: temporal history no longer covers its patch after migration");
    }
    history.gas_cell_identity_generation = rebuilt_identity.generation();
    for (const AmrTemporalBoundaryHistoryCellRecord& history_cell : history.cells) {
      if (history_cell.patch_local_cell >= rebuilt_patches.cell_count[patch_row]) {
        throw std::invalid_argument(
            "commitAmrPatchMigration: temporal history patch-local cell is out of range after migration");
      }
      const std::uint32_t expected_row = rebuilt_patches.first_cell[patch_row] +
          checkedIntegralNarrow<std::uint32_t>(
              history_cell.patch_local_cell,
              "commitAmrPatchMigration temporal history local-cell row");
      if (expected_row >= rebuilt_gas_cells.size() ||
          rebuilt_gas_cells.gas_cell_id[expected_row] != history_cell.gas_cell_id) {
        throw std::invalid_argument(
            "commitAmrPatchMigration: temporal history gas-cell identity changed during migration");
      }
    }
  }

  PendingFluxRegisterStore rebuilt_pending_flux_store;
  rebuilt_pending_flux_store.assign(std::move(rebuilt_pending_flux_records));
  AmrTemporalBoundaryHistoryStore rebuilt_temporal_history_store;
  rebuilt_temporal_history_store.assign(std::move(rebuilt_temporal_history_records));

  patches = std::move(rebuilt_patches);
  cells = std::move(rebuilt_cells);
  gas_cells = std::move(rebuilt_gas_cells);
  gas_composition = std::move(rebuilt_composition);
  black_holes.host_cell_index = std::move(rebuilt_black_hole_hosts);
  tracers.host_cell_index = std::move(rebuilt_tracer_hosts);
  gas_cell_identity = std::move(rebuilt_identity);
  pending_flux_registers = std::move(rebuilt_pending_flux_store);
  amr_temporal_boundary_history = std::move(rebuilt_temporal_history_store);
  bumpCellIndexGeneration();
}

}  // namespace cosmosim::core
