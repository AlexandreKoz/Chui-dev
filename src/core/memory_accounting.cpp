#include "cosmosim/core/memory_accounting.hpp"

#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_scheduler.hpp"

namespace cosmosim::core {
namespace {

[[nodiscard]] std::uint64_t checkedU64Add(std::uint64_t lhs, std::uint64_t rhs, std::string_view context) {
  if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
    throw std::overflow_error(std::string(context) + ": uint64 addition overflow");
  }
  return lhs + rhs;
}

[[nodiscard]] std::uint64_t checkedU64Multiply(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::string_view context) {
  if (lhs != 0U && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
    throw std::overflow_error(std::string(context) + ": uint64 multiplication overflow");
  }
  return lhs * rhs;
}

[[nodiscard]] MemoryClass defaultClassFor(
    MemorySubsystem subsystem,
    MemoryLifetime lifetime) noexcept {
  if (lifetime == MemoryLifetime::kUnknown) {
    return MemoryClass::kExternalRuntime;
  }
  switch (subsystem) {
    case MemorySubsystem::kParticles:
    case MemorySubsystem::kGasHydro:
    case MemorySubsystem::kSidecars:
      return lifetime == MemoryLifetime::kPersistent
          ? MemoryClass::kCanonicalPersistent
          : MemoryClass::kPhaseResident;
    case MemorySubsystem::kTree:
    case MemorySubsystem::kPmMesh:
      return lifetime == MemoryLifetime::kPersistent
          ? MemoryClass::kPersistentCache
          : MemoryClass::kPhaseResident;
    case MemorySubsystem::kActiveSets:
    case MemorySubsystem::kOutputBuffers:
      return MemoryClass::kPhaseResident;
    case MemorySubsystem::kMpiBuffers:
      return MemoryClass::kCommunication;
    case MemorySubsystem::kScratch:
      return MemoryClass::kScratchArena;
    case MemorySubsystem::kCount:
      return MemoryClass::kExternalRuntime;
  }
  return MemoryClass::kExternalRuntime;
}

void addTotals(MemoryReport& report, const MemoryEntry& entry) {
  const std::size_t index = memorySubsystemIndex(entry.subsystem);
  if (index >= static_cast<std::size_t>(MemorySubsystem::kCount)) {
    throw std::invalid_argument("memory report entry has invalid subsystem");
  }
  const std::uint64_t bytes = entry.owned_capacity_bytes;
  if (entry.lifetime == MemoryLifetime::kPersistent) {
    report.totals.persistent_by_subsystem[index] = checkedU64Add(
        report.totals.persistent_by_subsystem[index], bytes, "persistent subsystem memory total");
    report.totals.persistent_total_bytes = checkedU64Add(
        report.totals.persistent_total_bytes, bytes, "persistent memory total");
  } else if (entry.lifetime == MemoryLifetime::kTransient) {
    report.totals.transient_by_subsystem[index] = checkedU64Add(
        report.totals.transient_by_subsystem[index], bytes, "transient subsystem memory total");
    report.totals.transient_total_bytes = checkedU64Add(
        report.totals.transient_total_bytes, bytes, "transient memory total");
  } else if (entry.lifetime == MemoryLifetime::kUnknown) {
    report.totals.unknown_by_subsystem[index] = checkedU64Add(
        report.totals.unknown_by_subsystem[index], bytes, "unknown subsystem memory total");
    report.totals.unknown_total_bytes = checkedU64Add(
        report.totals.unknown_total_bytes, bytes, "unknown memory total");
  } else {
    throw std::invalid_argument("memory report entry has invalid lifetime");
  }
}

template <typename T>
void addOwnedWithClass(
    MemoryReportBuilder& builder,
    MemorySubsystem subsystem,
    MemoryLifetime lifetime,
    MemoryClass memory_class,
    std::string_view label,
    const T& container,
    std::string_view note = {}) {
  const std::uint64_t bytes = ownedCapacityBytesForContainer(container);
  builder.addEntry(MemoryEntry{.subsystem = subsystem,
                               .lifetime = lifetime,
                               .memory_class = memory_class,
                               .label = std::string(label),
                               .owned_capacity_bytes = bytes,
                               .referenced_bytes = 0,
                               .high_water_bytes = bytes,
                               .estimate_only = false,
                               .uncertainty_note = std::string(note)});
}

template <typename T>
void addOwned(
    MemoryReportBuilder& builder,
    MemorySubsystem subsystem,
    MemoryLifetime lifetime,
    std::string_view label,
    const T& container) {
  const std::uint64_t bytes = ownedCapacityBytesForContainer(container);
  builder.addEntry(MemoryEntry{.subsystem = subsystem,
                               .lifetime = lifetime,
                               .memory_class = defaultClassFor(subsystem, lifetime),
                               .label = std::string(label),
                               .owned_capacity_bytes = bytes,
                               .referenced_bytes = 0,
                               .high_water_bytes = bytes,
                               .estimate_only = false,
                               .uncertainty_note = {}});
}

void addEstimate(
    MemoryReportBuilder& builder,
    MemorySubsystem subsystem,
    MemoryLifetime lifetime,
    std::string_view label,
    std::uint64_t bytes,
    std::string_view note) {
  builder.addEntry(MemoryEntry{.subsystem = subsystem,
                               .lifetime = lifetime,
                               .memory_class = defaultClassFor(subsystem, lifetime),
                               .label = std::string(label),
                               .owned_capacity_bytes = bytes,
                               .high_water_bytes = bytes,
                               .estimate_only = true,
                               .uncertainty_note = std::string(note)});
}

void accountParticleSoa(MemoryReportBuilder& builder, const ParticleSoa& particles) {
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.position_x_comoving", particles.position_x_comoving);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.position_y_comoving", particles.position_y_comoving);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.position_z_comoving", particles.position_z_comoving);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.velocity_x_peculiar", particles.velocity_x_peculiar);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.velocity_y_peculiar", particles.velocity_y_peculiar);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.velocity_z_peculiar", particles.velocity_z_peculiar);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.mass_code", particles.mass_code);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.time_bin_mirror", particles.time_bin);
}

void accountCellSoa(MemoryReportBuilder& builder, const CellSoa& cells) {
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "cells.center_x_comoving", cells.center_x_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "cells.center_y_comoving", cells.center_y_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "cells.center_z_comoving", cells.center_z_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "cells.mass_code", cells.mass_code);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "cells.time_bin_mirror", cells.time_bin);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "cells.patch_index", cells.patch_index);
}

void accountGasSidecar(MemoryReportBuilder& builder, const GasCellSidecar& gas) {
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "gas_cells.gas_cell_id", gas.gas_cell_id);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "gas_cells.parent_particle_id", gas.parent_particle_id);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "gas_cells.velocity_x_peculiar", gas.velocity_x_peculiar);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "gas_cells.velocity_y_peculiar", gas.velocity_y_peculiar);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "gas_cells.velocity_z_peculiar", gas.velocity_z_peculiar);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "gas_cells.density_code", gas.density_code);
  addOwnedWithClass(
      builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent,
      MemoryClass::kPersistentCache, "gas_cells.pressure_code", gas.pressure_code,
      "restart/schema compatibility mirror; hydro conserved updates refresh it");
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "gas_cells.internal_energy_code", gas.internal_energy_code);
  addOwnedWithClass(
      builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent,
      MemoryClass::kPersistentCache, "gas_cells.temperature_code", gas.temperature_code,
      "restart/source compatibility cache derived from thermodynamic state");
  addOwnedWithClass(
      builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent,
      MemoryClass::kPersistentCache, "gas_cells.sound_speed_code", gas.sound_speed_code,
      "restart/scheduler compatibility cache derived from thermodynamic closure");
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "gas_cells.metal_mass_code", gas.metal_mass_code);
}

void accountParticleSidecar(MemoryReportBuilder& builder, const ParticleSidecar& sidecar) {
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.particle_id", sidecar.particle_id);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.sfc_key", sidecar.sfc_key);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.species_tag", sidecar.species_tag);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.particle_flags", sidecar.particle_flags);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.owning_rank", sidecar.owning_rank);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.last_drift_time_code", sidecar.last_drift_time_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.last_drift_scale_factor", sidecar.last_drift_scale_factor);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.gravity_softening_comoving", sidecar.gravity_softening_comoving);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.has_gravity_softening_override", sidecar.has_gravity_softening_override);
}

void accountSpeciesIndex(MemoryReportBuilder& builder, const ParticleSpeciesIndex& species_index) {
  for (std::size_t i = 0; i < species_index.global_index_by_species.size(); ++i) {
    addOwned(builder,
             MemorySubsystem::kSidecars,
             MemoryLifetime::kPersistent,
             "particle_species_index.global_index_by_species." + std::to_string(i),
             species_index.global_index_by_species[i]);
  }
  addOwned(builder,
           MemorySubsystem::kSidecars,
           MemoryLifetime::kPersistent,
           "particle_species_index.local_index_by_global",
           species_index.local_index_by_global);
}

void accountStarSidecar(MemoryReportBuilder& builder, const StarParticleSidecar& stars) {
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.particle_index", stars.particle_index);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.birth_key", stars.birth_key);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.parent_gas_cell_id", stars.parent_gas_cell_id);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.birth_tick", stars.birth_tick);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.birth_ordinal", stars.birth_ordinal);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.formation_scale_factor", stars.formation_scale_factor);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.birth_mass_code", stars.birth_mass_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.metallicity_mass_fraction", stars.metallicity_mass_fraction);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.stellar_age_years_last", stars.stellar_age_years_last);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.stellar_returned_mass_cumulative_code", stars.stellar_returned_mass_cumulative_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.stellar_returned_metals_cumulative_code", stars.stellar_returned_metals_cumulative_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.stellar_newly_synthesized_metals_cumulative_code", stars.stellar_newly_synthesized_metals_cumulative_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.stellar_feedback_energy_cumulative_erg", stars.stellar_feedback_energy_cumulative_erg);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.enrichment_carry_mass_code", stars.enrichment_carry_mass_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.enrichment_carry_metals_code", stars.enrichment_carry_metals_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.enrichment_carry_feedback_energy_erg", stars.enrichment_carry_feedback_energy_erg);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.enrichment_carry_momentum_code", stars.enrichment_carry_momentum_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.stellar_deposited_mass_cumulative_code", stars.stellar_deposited_mass_cumulative_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.stellar_deposited_metals_cumulative_code", stars.stellar_deposited_metals_cumulative_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.stellar_deposited_feedback_energy_cumulative_erg", stars.stellar_deposited_feedback_energy_cumulative_erg);
  for (std::size_t i = 0; i < stars.stellar_returned_mass_channel_cumulative_code.size(); ++i) {
    addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.returned_mass_channel." + std::to_string(i), stars.stellar_returned_mass_channel_cumulative_code[i]);
    addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.returned_metals_channel." + std::to_string(i), stars.stellar_returned_metals_channel_cumulative_code[i]);
    addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.feedback_energy_channel." + std::to_string(i), stars.stellar_feedback_energy_channel_cumulative_erg[i]);
  }
}

void accountBlackHoleSidecar(MemoryReportBuilder& builder, const BlackHoleParticleSidecar& black_holes) {
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.particle_index", black_holes.particle_index);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.host_cell_index", black_holes.host_cell_index);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.subgrid_mass_code", black_holes.subgrid_mass_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.accretion_rate_code", black_holes.accretion_rate_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.feedback_energy_code", black_holes.feedback_energy_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.eddington_ratio", black_holes.eddington_ratio);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.cumulative_accreted_mass_code", black_holes.cumulative_accreted_mass_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.cumulative_feedback_energy_code", black_holes.cumulative_feedback_energy_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.duty_cycle_active_time_code", black_holes.duty_cycle_active_time_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.duty_cycle_total_time_code", black_holes.duty_cycle_total_time_code);
}

void accountTracerSidecar(MemoryReportBuilder& builder, const TracerParticleSidecar& tracers) {
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "tracers.particle_index", tracers.particle_index);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "tracers.parent_particle_id", tracers.parent_particle_id);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "tracers.injection_step", tracers.injection_step);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "tracers.host_cell_index", tracers.host_cell_index);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "tracers.mass_fraction_of_host", tracers.mass_fraction_of_host);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "tracers.last_host_mass_code", tracers.last_host_mass_code);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "tracers.cumulative_exchanged_mass_code", tracers.cumulative_exchanged_mass_code);
}

void accountPatchSoa(MemoryReportBuilder& builder, const PatchSoa& patches) {
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.patch_id", patches.patch_id);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.level", patches.level);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.first_cell", patches.first_cell);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.cell_count", patches.cell_count);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.parent_patch_id", patches.parent_patch_id);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.morton_key", patches.morton_key);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.origin_x_comoving", patches.origin_x_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.origin_y_comoving", patches.origin_y_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.origin_z_comoving", patches.origin_z_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.extent_x_comoving", patches.extent_x_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.extent_y_comoving", patches.extent_y_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.extent_z_comoving", patches.extent_z_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.cell_dim_x", patches.cell_dim_x);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.cell_dim_y", patches.cell_dim_y);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.cell_dim_z", patches.cell_dim_z);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "patches.owning_rank", patches.owning_rank);
}

void accountAmrIdentityAndSynchronization(
    MemoryReportBuilder& builder,
    const SimulationState& state) {
  builder.addEntry(MemoryEntry{
      .subsystem = MemorySubsystem::kGasHydro,
      .lifetime = MemoryLifetime::kPersistent,
      .memory_class = MemoryClass::kCanonicalPersistent,
      .label = "gas_cell_identity.records",
      .current_size_bytes = state.gas_cell_identity.recordsCurrentSizeBytes(),
      .owned_capacity_bytes = state.gas_cell_identity.recordsOwnedCapacityBytes(),
      .high_water_bytes = state.gas_cell_identity.recordsOwnedCapacityBytes(),
      .estimate_only = false,
      .uncertainty_note = {}});
  const std::uint64_t lookup_bytes =
      state.gas_cell_identity.lookupOwnedBytesConservative();
  builder.addEntry(MemoryEntry{
      .subsystem = MemorySubsystem::kSidecars,
      .lifetime = MemoryLifetime::kPersistent,
      .memory_class = MemoryClass::kPersistentCache,
      .label = "gas_cell_identity.lookup_tables",
      .current_size_bytes = lookup_bytes,
      .owned_capacity_bytes = lookup_bytes,
      .high_water_bytes = lookup_bytes,
      .estimate_only = false,
      .uncertainty_note =
          "Conservative unordered-map bucket/node ownership model; allocator bookkeeping is covered by process safety margin/RSS reconciliation."});
  builder.addEntry(MemoryEntry{
      .subsystem = MemorySubsystem::kGasHydro,
      .lifetime = MemoryLifetime::kPersistent,
      .memory_class = MemoryClass::kCanonicalPersistent,
      .label = "amr.pending_flux_registers",
      .current_size_bytes = state.pending_flux_registers.currentSizeBytes(),
      .owned_capacity_bytes = state.pending_flux_registers.ownedCapacityBytes(),
      .high_water_bytes = state.pending_flux_registers.ownedCapacityBytes(),
      .estimate_only = false,
      .uncertainty_note = {}});
  builder.addEntry(MemoryEntry{
      .subsystem = MemorySubsystem::kGasHydro,
      .lifetime = MemoryLifetime::kPersistent,
      .memory_class = MemoryClass::kCanonicalPersistent,
      .label = "amr.temporal_boundary_history",
      .current_size_bytes = state.amr_temporal_boundary_history.currentSizeBytes(),
      .owned_capacity_bytes = state.amr_temporal_boundary_history.ownedCapacityBytes(),
      .high_water_bytes = state.amr_temporal_boundary_history.ownedCapacityBytes(),
      .estimate_only = false,
      .uncertainty_note =
          "Includes outer record capacity and nested per-patch temporal-cell vector capacity."});
}

void accountModulePayloads(MemoryReportBuilder& builder, const ModuleSidecarRegistry& registry) {
  for (const ModuleSidecarBlock* block : registry.blocksSortedByName()) {
    if (block == nullptr) {
      continue;
    }
    addOwned(builder,
             MemorySubsystem::kSidecars,
             MemoryLifetime::kPersistent,
             "module_sidecar." + block->module_name + ".payload",
             block->payload);
    addOwned(builder,
             MemorySubsystem::kSidecars,
             MemoryLifetime::kPersistent,
             "module_sidecar." + block->module_name + ".particle_id_by_row",
             block->particle_id_by_row);
  }
}

void accountWorkspace(MemoryReportBuilder& builder, const TransientStepWorkspace& workspace) {
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_id", workspace.particle_id);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_species_tag", workspace.particle_species_tag);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_position_x_comoving", workspace.particle_position_x_comoving);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_position_y_comoving", workspace.particle_position_y_comoving);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_position_z_comoving", workspace.particle_position_z_comoving);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_velocity_x_peculiar", workspace.particle_velocity_x_peculiar);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_velocity_y_peculiar", workspace.particle_velocity_y_peculiar);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_velocity_z_peculiar", workspace.particle_velocity_z_peculiar);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_mass_code", workspace.particle_mass_code);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.gravity_particle_index", workspace.gravity_particle_index);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.hydro_cell_index", workspace.hydro_cell_index);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.hydro_gas_cell_id", workspace.hydro_gas_cell_id);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.hydro_local_cell_row", workspace.hydro_local_cell_row);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.hydro_cell_center_x_comoving", workspace.hydro_cell_center_x_comoving);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.hydro_cell_center_y_comoving", workspace.hydro_cell_center_y_comoving);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.hydro_cell_center_z_comoving", workspace.hydro_cell_center_z_comoving);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.hydro_cell_mass_code", workspace.hydro_cell_mass_code);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.hydro_cell_density_code", workspace.hydro_cell_density_code);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.hydro_cell_pressure_code", workspace.hydro_cell_pressure_code);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.cell_center_x_comoving", workspace.cell_center_x_comoving);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.cell_center_y_comoving", workspace.cell_center_y_comoving);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.cell_center_z_comoving", workspace.cell_center_z_comoving);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.cell_mass_code", workspace.cell_mass_code);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.cell_patch_index", workspace.cell_patch_index);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.cell_density_code", workspace.cell_density_code);
  addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.cell_pressure_code", workspace.cell_pressure_code);
  addOwned(builder, MemorySubsystem::kScratch, MemoryLifetime::kTransient, "workspace.hydro_recon_gradient_x", workspace.hydro_recon_gradient_x);
  addOwned(builder, MemorySubsystem::kScratch, MemoryLifetime::kTransient, "workspace.hydro_recon_gradient_y", workspace.hydro_recon_gradient_y);
  addOwned(builder, MemorySubsystem::kScratch, MemoryLifetime::kTransient, "workspace.hydro_recon_gradient_z", workspace.hydro_recon_gradient_z);
  builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kScratch,
                               .lifetime = MemoryLifetime::kTransient,
                               .memory_class = MemoryClass::kScratchArena,
                               .label = "workspace.scratch_arena",
                               .current_size_bytes = static_cast<std::uint64_t>(workspace.scratch.usedBytes()),
                               .owned_capacity_bytes = static_cast<std::uint64_t>(workspace.scratch.capacityBytes()),
                               .referenced_bytes = 0,
                               .high_water_bytes = static_cast<std::uint64_t>(workspace.scratch.capacityHighWaterBytes()),
                               .estimate_only = false,
                               .governed_commitment = workspace.scratch.governed(),
                               .uncertainty_note = "logical_high_water_bytes=" +
                                   std::to_string(workspace.scratch.logicalHighWaterBytes())});
}

}  // namespace

void MemoryReportBuilder::addEntry(MemoryEntry entry) {
  if (entry.label.empty()) {
    throw std::invalid_argument("memory report entry label must not be empty");
  }
  if (memorySubsystemIndex(entry.subsystem) >= static_cast<std::size_t>(MemorySubsystem::kCount)) {
    throw std::invalid_argument("memory report entry has invalid subsystem");
  }
  if (entry.lifetime != MemoryLifetime::kPersistent && entry.lifetime != MemoryLifetime::kTransient &&
      entry.lifetime != MemoryLifetime::kUnknown) {
    throw std::invalid_argument("memory report entry has invalid lifetime");
  }
  if (!entry.memory_class.has_value()) {
    entry.memory_class = defaultClassFor(entry.subsystem, entry.lifetime);
  }
  if (memoryClassIndex(*entry.memory_class) >=
      static_cast<std::size_t>(MemoryClass::kCount)) {
    throw std::invalid_argument("memory report entry has invalid memory class");
  }
  addTotals(m_report, entry);
  m_report.entries.push_back(std::move(entry));
}

MemoryReport MemoryReportBuilder::finish() && {
  for (std::size_t i = 0; i < static_cast<std::size_t>(MemorySubsystem::kCount); ++i) {
    if (m_report.totals.persistent_by_subsystem[i] == 0 && m_report.totals.transient_by_subsystem[i] == 0 &&
        m_report.totals.unknown_by_subsystem[i] == 0) {
      m_report.entries.push_back(MemoryEntry{.subsystem = static_cast<MemorySubsystem>(i),
                                             .lifetime = MemoryLifetime::kUnknown,
                                             .label = "category_present",
                                             .owned_capacity_bytes = 0,
                                             .referenced_bytes = 0,
                                             .high_water_bytes = 0,
                                             .estimate_only = false,
                                             .uncertainty_note = {}});
    }
  }
  return std::move(m_report);
}

std::size_t memorySubsystemIndex(MemorySubsystem subsystem) noexcept { return static_cast<std::size_t>(subsystem); }

std::string_view memorySubsystemLabel(MemorySubsystem subsystem) noexcept {
  switch (subsystem) {
    case MemorySubsystem::kParticles:
      return "particles";
    case MemorySubsystem::kGasHydro:
      return "gas_hydro";
    case MemorySubsystem::kTree:
      return "tree";
    case MemorySubsystem::kPmMesh:
      return "pm_mesh";
    case MemorySubsystem::kSidecars:
      return "sidecars";
    case MemorySubsystem::kActiveSets:
      return "active_sets";
    case MemorySubsystem::kMpiBuffers:
      return "mpi_buffers";
    case MemorySubsystem::kScratch:
      return "scratch";
    case MemorySubsystem::kOutputBuffers:
      return "output_buffers";
    case MemorySubsystem::kCount:
      return "invalid";
  }
  return "invalid";
}
std::string_view memoryLifetimeLabel(MemoryLifetime lifetime) noexcept {
  switch (lifetime) {
    case MemoryLifetime::kPersistent:
      return "persistent";
    case MemoryLifetime::kTransient:
      return "transient";
    case MemoryLifetime::kUnknown:
      return "unknown";
  }
  return "unknown";
}

MemoryReport collectSimulationMemoryReport(const SimulationState& state, const TransientStepWorkspace* workspace) {
  MemoryReportBuilder builder;
  accountParticleSoa(builder, state.particles);
  accountCellSoa(builder, state.cells);
  accountGasSidecar(builder, state.gas_cells);
  accountParticleSidecar(builder, state.particle_sidecar);
  accountSpeciesIndex(builder, state.particle_species_index);
  accountStarSidecar(builder, state.star_particles);
  accountBlackHoleSidecar(builder, state.black_holes);
  accountTracerSidecar(builder, state.tracers);
  accountPatchSoa(builder, state.patches);
  accountAmrIdentityAndSynchronization(builder, state);
  accountModulePayloads(builder, state.sidecars);
  // M01 (patch 50): the element/species/ion/dust/SGS composition sidecar is
  // canonical persistent state (82 lanes per cell for the largest C02
  // layout, plus dense or sparse detailed rows and chemistry history). It
  // was previously absent from this report, so the governor baseline, the
  // global-step rollback state-copy admission and the retained-capacity
  // baseline reconciliation all undercounted it.
  if (state.gas_composition.enabled() || state.gas_composition.ownedCapacityBytes() != 0U) {
    const std::uint64_t composition_capacity = state.gas_composition.ownedCapacityBytes();
    builder.addEntry(MemoryEntry{
        .subsystem = MemorySubsystem::kSidecars,
        .lifetime = MemoryLifetime::kPersistent,
        .memory_class = MemoryClass::kCanonicalPersistent,
        .label = "gas_composition.lanes_sparse_rows_and_history",
        .current_size_bytes = state.gas_composition.currentSizeBytes(),
        .owned_capacity_bytes = composition_capacity,
        .referenced_bytes = 0U,
        .high_water_bytes = composition_capacity,
        .estimated_next_step_bytes = composition_capacity,
        .estimate_only = false,
        .governed_commitment = false,
        .uncertainty_note = "capacity of every composition lane, sparse detailed-row store and per-row "
                            "chemistry history vector (GasCompositionSidecar::ownedCapacityBytes)"});
  }
  // Representation policy is scalar authority. Population-scale compact DMO
  // metadata is accounted only when its backing vectors actually exist above.
  constexpr std::uint64_t particle_metadata_policy_scalar_bytes =
      sizeof(ParticleMetadataRepresentation) + sizeof(std::uint32_t) * 3U +
      sizeof(std::uint64_t) + sizeof(bool) + sizeof(double) * 2U;
  builder.addEntry(MemoryEntry{
      .subsystem = MemorySubsystem::kSidecars,
      .lifetime = MemoryLifetime::kPersistent,
      .memory_class = MemoryClass::kCanonicalPersistent,
      .label = "particle_metadata.representation_policy_scalars",
      .current_size_bytes = particle_metadata_policy_scalar_bytes,
      .owned_capacity_bytes = particle_metadata_policy_scalar_bytes,
      .high_water_bytes = particle_metadata_policy_scalar_bytes,
      .estimated_next_step_bytes = particle_metadata_policy_scalar_bytes,
      .estimate_only = false,
      .uncertainty_note = state.hasHomogeneousDmoMetadata()
          ? "homogeneous DMO logical metadata is scalar/implicit; only materialized exception vectors are charged separately"
          : "generic particle metadata representation policy scalars; dense lanes are charged separately"});

  builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kTree,
                               .lifetime = MemoryLifetime::kUnknown,
                               .label = "tree.external_allocations",
                               .estimate_only = true,
                               .uncertainty_note = "Tree node storage is subsystem-internal unless a tree solver report is merged into this snapshot."});
  builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kPmMesh,
                               .lifetime = MemoryLifetime::kUnknown,
                               .label = "pm_mesh.external_allocations",
                               .estimate_only = true,
                               .uncertainty_note = "PM mesh/FFTW/GPU allocations are reported only by explicit solver hooks; third-party internals remain unknown."});
  builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kMpiBuffers,
                               .lifetime = MemoryLifetime::kUnknown,
                               .label = "mpi_buffers.external_allocations",
                               .estimate_only = true,
                               .uncertainty_note = "MPI implementation buffers are outside owned host-state accounting."});
  builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kOutputBuffers,
                               .lifetime = MemoryLifetime::kUnknown,
                               .label = "output_buffers.external_allocations",
                               .estimate_only = true,
                               .uncertainty_note = "HDF5/library staging allocations are intentionally not counted as restart truth."});

  if (workspace != nullptr) {
    accountWorkspace(builder, *workspace);
  }

  MemoryReport report = std::move(builder).finish();
  report.notes.push_back("Owned bytes use container capacity*sizeof(T); spans/views are non-owning and report zero owned bytes.");
  report.notes.push_back("Persistent/transient classification follows restart ownership: scratch and active buffers are transient diagnostics, not restart truth.");
  report.notes.push_back("Unknown/estimate entries intentionally avoid claiming allocator, MPI, GPU, FFTW, or HDF5 internal precision.");
  return report;
}

MemoryReport collectSchedulerMemoryReport(
    const HierarchicalTimeBinScheduler& particle_scheduler,
    const HierarchicalTimeBinScheduler& gas_cell_scheduler) {
  MemoryReportBuilder builder;
  const auto add_scheduler = [&builder](
      std::string label,
      const HierarchicalTimeBinScheduler& scheduler) {
    builder.addEntry(MemoryEntry{
        .subsystem = MemorySubsystem::kSidecars,
        .lifetime = MemoryLifetime::kPersistent,
        .memory_class = MemoryClass::kCanonicalPersistent,
        .label = std::move(label),
        .current_size_bytes = scheduler.logicalSizeBytes(),
        .owned_capacity_bytes = scheduler.ownedCapacityBytes(),
        .referenced_bytes = 0U,
        .high_water_bytes = scheduler.ownedCapacityHighWaterBytes(),
        .estimated_next_step_bytes = scheduler.ownedCapacityBytes(),
        .estimate_only = false,
        .governed_commitment = false,
        .uncertainty_note =
            scheduler.representation() == SchedulerRepresentation::kUniformRungZero
                ? "uniform rung-zero scheduler: one uint32 identity active-index lane plus scalar authority; absent generic lanes are not charged"
                : "generic hierarchical scheduler logical bytes, retained capacity, and retained-capacity high-water"});
  };
  add_scheduler("scheduler.particle_owned_state", particle_scheduler);
  add_scheduler("scheduler.gas_cell_owned_state", gas_cell_scheduler);
  MemoryReport report = std::move(builder).finish();
  report.notes.push_back(
      "Scheduler retained capacities are explicit persistent ownership entries and participate directly in governor baseline reconciliation.");
  return report;
}

MemoryReport mergeMemoryReports(std::span<const MemoryReport> reports) {
  MemoryReportBuilder builder;
  std::vector<std::string> notes;
  std::unordered_set<std::string> ownership_keys;
  for (const MemoryReport& report : reports) {
    for (const MemoryEntry& entry : report.entries) {
      if (entry.label == "category_present") {
        continue;
      }
      const std::string key = std::to_string(memorySubsystemIndex(entry.subsystem)) + ":" +
                              std::to_string(static_cast<unsigned int>(entry.lifetime)) + ":" + entry.label;
      if (!ownership_keys.insert(key).second) {
        throw std::invalid_argument("duplicate memory ownership entry while merging reports: " + entry.label);
      }
      builder.addEntry(entry);
    }
    notes.insert(notes.end(), report.notes.begin(), report.notes.end());
  }
  MemoryReport merged = std::move(builder).finish();
  merged.notes = std::move(notes);
  merged.notes.push_back("Merged runtime memory report; duplicate owners are avoided by subsystem-level append hooks.");
  return merged;
}

std::uint64_t memoryReportBaselineOwnedBytes(const MemoryReport& report) {
  std::uint64_t total = 0U;
  for (const MemoryEntry& entry : report.entries) {
    if (entry.estimate_only || entry.governed_commitment ||
        entry.lifetime == MemoryLifetime::kUnknown) {
      continue;
    }
    total = checkedU64Add(
        total, entry.owned_capacity_bytes,
        "memory report baseline owned byte reconciliation");
  }
  return total;
}

void attachMemoryGovernorSnapshot(
    MemoryReport& report,
    const MemoryGovernor& governor) {
  report.governor_snapshot = governor.snapshot();
}

void attachProcessMemoryObservation(
    MemoryReport& report,
    ProcessMemoryObservation observation) {
  report.process_memory_observation = observation;
  ProcessMemoryReconciliation reconciliation;
  if (report.governor_snapshot.has_value()) {
    const MemoryGovernorSnapshot& governor = *report.governor_snapshot;
    reconciliation.current_declared_residency_bytes = checkedMemoryBytesAdd(
        governor.baseline_owned_bytes,
        governor.external_runtime_reserve_bytes,
        "process current declared baseline/external residency");
    reconciliation.current_declared_residency_bytes = checkedMemoryBytesAdd(
        reconciliation.current_declared_residency_bytes,
        governor.committed_bytes,
        "process current declared committed residency");
  } else {
    reconciliation.current_declared_residency_bytes =
        memoryReportBaselineOwnedBytes(report);
  }
  reconciliation.known_accounted_bytes =
      reconciliation.current_declared_residency_bytes;
  reconciliation.observed_rss_bytes = observation.current_rss_bytes;
  reconciliation.observed_peak_rss_bytes = observation.peak_rss_bytes;
  reconciliation.observed_pss_bytes = observation.pss_bytes;
  if (observation.current_rss_bytes.has_value()) {
    reconciliation.unexplained_resident_bytes =
        *observation.current_rss_bytes > reconciliation.known_accounted_bytes
        ? *observation.current_rss_bytes - reconciliation.known_accounted_bytes
        : 0U;
    if (reconciliation.known_accounted_bytes != 0U) {
      reconciliation.observed_to_known_ratio =
          static_cast<double>(*observation.current_rss_bytes) /
          static_cast<double>(reconciliation.known_accounted_bytes);
    }
  }
  report.process_memory_reconciliation = reconciliation;
}

std::uint64_t memoryReportCommunicationHighWaterBytes(
    const MemoryReport& report) {
  std::uint64_t total = 0U;
  for (const MemoryEntry& entry : report.entries) {
    if (entry.subsystem != MemorySubsystem::kMpiBuffers ||
        entry.label == "category_present") {
      continue;
    }
    total = checkedU64Add(
        total,
        entry.high_water_bytes,
        "memory report communication high-water aggregation");
  }
  return total;
}

MemoryReport estimatePreRunMemoryBudget(const MemoryBudgetEstimateInput& input) {
  MemoryReportBuilder builder;
  const auto bytes = [](std::uint64_t count, std::uint64_t elem_bytes) {
    return checkedU64Multiply(count, elem_bytes, "pre-run memory estimate");
  };

  addEstimate(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "estimate.particles.hot_soa", bytes(input.particle_capacity, 7 * sizeof(double) + sizeof(std::uint8_t)), "capacity estimate from configured particle count; allocator overhead not included");
  addEstimate(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "estimate.particle_sidecars", bytes(input.particle_capacity, 3 * sizeof(std::uint64_t) + 3 * sizeof(std::uint32_t) + 3 * sizeof(double) + sizeof(std::uint8_t)), "particle metadata/softening sidecar capacity estimate");
  addEstimate(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "estimate.gas_hydro.cells", bytes(input.gas_cell_capacity, 13 * sizeof(double) + 2 * sizeof(std::uint64_t) + 2 * sizeof(std::uint32_t) + sizeof(std::uint8_t)), "cell and gas thermodynamic capacity estimate; hydro conserved solver storage must be reported separately by hydro module");
  addEstimate(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "estimate.star_sidecars", bytes(input.star_capacity, 2 * sizeof(std::uint32_t) + 16 * sizeof(double)), "star sidecar estimate includes channel cumulative fields");
  addEstimate(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "estimate.black_hole_sidecars", bytes(input.black_hole_capacity, 2 * sizeof(std::uint32_t) + 8 * sizeof(double)), "black-hole sidecar capacity estimate");
  addEstimate(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "estimate.tracer_sidecars", bytes(input.tracer_capacity, 2 * sizeof(std::uint64_t) + 2 * sizeof(std::uint32_t) + 3 * sizeof(double)), "tracer sidecar capacity estimate");
  addEstimate(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "estimate.active_particle_views", bytes(input.active_particle_capacity, sizeof(std::uint32_t) + 7 * sizeof(double)), "compact active particle view estimate");
  addEstimate(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "estimate.active_cell_views", bytes(input.active_cell_capacity, sizeof(std::uint32_t) + 7 * sizeof(double)), "compact active cell view estimate");
  addEstimate(builder, MemorySubsystem::kTree, MemoryLifetime::kTransient, "estimate.tree_nodes", bytes(input.tree_node_capacity, 96), "tree node byte estimate uses a conservative host-side node footprint until tree hooks report exact capacity");
  addEstimate(builder, MemorySubsystem::kPmMesh, MemoryLifetime::kTransient, "estimate.pm_mesh", bytes(input.pm_grid_cells, 4 * sizeof(double)), "PM mesh estimate covers owned real-space fields only; FFTW/GPU work buffers remain external unknowns");
  addEstimate(builder, MemorySubsystem::kMpiBuffers, MemoryLifetime::kTransient, "estimate.mpi_exchange_particles", bytes(input.mpi_exchange_particle_capacity, 128), "MPI payload estimate excludes MPI library internal buffers");
  addEstimate(builder, MemorySubsystem::kOutputBuffers, MemoryLifetime::kTransient, "estimate.output_buffers", input.output_buffer_bytes, "configured output staging estimate; HDF5 internals are unknown");

  MemoryReport report = std::move(builder).finish();
  report.notes.push_back("Pre-run entries are estimates derived from requested capacities, not measured allocations.");
  report.notes.push_back("Third-party and device allocator internals are intentionally labelled as unknown/estimate-only.");
  return report;
}

std::string formatMemoryReportHumanReadable(const MemoryReport& report) {
  std::ostringstream out;
  out << "memory_report persistent_total_bytes=" << report.totals.persistent_total_bytes
      << " transient_total_bytes=" << report.totals.transient_total_bytes
      << " unknown_total_bytes=" << report.totals.unknown_total_bytes << "\n";
  if (report.governor_snapshot.has_value()) {
    const MemoryGovernorSnapshot& governor = *report.governor_snapshot;
    out << " governor hard_limit_bytes=" << governor.hard_limit_bytes
        << " baseline_owned_bytes=" << governor.baseline_owned_bytes
        << " committed_bytes=" << governor.committed_bytes
        << " reserved_bytes=" << governor.reserved_bytes
        << " accounted_bytes=" << governor.accounted_bytes
        << " headroom_bytes=" << governor.headroom_bytes
        << " pressure=" << memoryPressureLabel(governor.pressure)
        << " rejection_count=" << governor.rejection_count << "\n";
  }
  if (report.process_memory_reconciliation.has_value()) {
    const ProcessMemoryReconciliation& process = *report.process_memory_reconciliation;
    out << " process_memory current_declared_residency_bytes="
        << process.current_declared_residency_bytes
        << " known_accounted_bytes=" << process.known_accounted_bytes;
    if (process.observed_rss_bytes.has_value()) {
      out << " observed_rss_bytes=" << *process.observed_rss_bytes;
    }
    if (process.observed_peak_rss_bytes.has_value()) {
      out << " observed_peak_rss_bytes=" << *process.observed_peak_rss_bytes;
    }
    if (process.observed_pss_bytes.has_value()) {
      out << " observed_pss_bytes=" << *process.observed_pss_bytes;
    }
    if (process.unexplained_resident_bytes.has_value()) {
      out << " unexplained_resident_bytes=" << *process.unexplained_resident_bytes;
    }
    if (process.observed_to_known_ratio.has_value()) {
      out << " observed_to_known_ratio=" << *process.observed_to_known_ratio;
    }
    out << "\n";
  }
  for (const MemoryEntry& entry : report.entries) {
    out << " - subsystem=" << memorySubsystemLabel(entry.subsystem)
        << " lifetime=" << memoryLifetimeLabel(entry.lifetime)
        << " class=" << (entry.memory_class.has_value()
            ? memoryClassLabel(*entry.memory_class) : std::string_view("unclassified"))
        << " label=" << entry.label;
    if (entry.current_size_bytes > 0) {
      out << " current_size_bytes=" << entry.current_size_bytes;
    }
    out << " owned_capacity_bytes=" << entry.owned_capacity_bytes;
    if (entry.high_water_bytes > 0) {
      out << " high_water_bytes=" << entry.high_water_bytes;
    }
    if (entry.estimated_next_step_bytes > 0) {
      out << " estimated_next_step_bytes=" << entry.estimated_next_step_bytes;
    }
    if (entry.referenced_bytes > 0) {
      out << " referenced_bytes=" << entry.referenced_bytes;
    }
    if (entry.estimate_only) {
      out << " estimate_only=true";
    }
    if (entry.governed_commitment) {
      out << " governed_commitment=true";
    }
    if (!entry.uncertainty_note.empty()) {
      out << " note=\"" << entry.uncertainty_note << "\"";
    }
    out << "\n";
  }
  return out.str();
}

}  // namespace cosmosim::core
