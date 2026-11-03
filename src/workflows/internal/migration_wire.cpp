#include "workflows/internal/migration_wire.hpp"

#include "cosmosim/core/checked_arithmetic.hpp"

#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace cosmosim::workflows::internal::migration_wire {
namespace {

class WireWriter {
 public:
  void u8(std::uint8_t value) { m_bytes.push_back(value); }

  void u16(std::uint16_t value) {
    for (unsigned shift = 0; shift < 16U; shift += 8U) {
      m_bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
    }
  }

  void u32(std::uint32_t value) {
    for (unsigned shift = 0; shift < 32U; shift += 8U) {
      m_bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
    }
  }

  void i32(std::int32_t value) { u32(std::bit_cast<std::uint32_t>(value)); }

  void u64(std::uint64_t value) {
    for (unsigned shift = 0; shift < 64U; shift += 8U) {
      m_bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffULL));
    }
  }

  void f64(double value) { u64(std::bit_cast<std::uint64_t>(value)); }

  void string(std::string_view value) {
    u64(static_cast<std::uint64_t>(value.size()));
    bytes(std::span<const std::uint8_t>(
        reinterpret_cast<const std::uint8_t*>(value.data()), value.size()));
  }

  void bytePayload(std::span<const std::byte> value) {
    u64(static_cast<std::uint64_t>(value.size()));
    bytes(std::span<const std::uint8_t>(
        reinterpret_cast<const std::uint8_t*>(value.data()), value.size()));
  }

  void bytes(std::span<const std::uint8_t> value) {
    m_bytes.insert(m_bytes.end(), value.begin(), value.end());
  }

  [[nodiscard]] std::vector<std::uint8_t> take() && { return std::move(m_bytes); }

 private:
  std::vector<std::uint8_t> m_bytes;
};

class WireReader {
 public:
  explicit WireReader(std::span<const std::uint8_t> bytes) : m_bytes(bytes) {}

  [[nodiscard]] std::uint8_t u8(std::string_view label) {
    require(sizeof(std::uint8_t), label);
    return m_bytes[m_offset++];
  }

  [[nodiscard]] std::uint16_t u16(std::string_view label) {
    require(sizeof(std::uint16_t), label);
    std::uint16_t value = 0U;
    for (unsigned shift = 0; shift < 16U; shift += 8U) {
      value |= static_cast<std::uint16_t>(m_bytes[m_offset++]) << shift;
    }
    return value;
  }

  [[nodiscard]] std::uint32_t u32(std::string_view label) {
    require(sizeof(std::uint32_t), label);
    std::uint32_t value = 0U;
    for (unsigned shift = 0; shift < 32U; shift += 8U) {
      value |= static_cast<std::uint32_t>(m_bytes[m_offset++]) << shift;
    }
    return value;
  }

  [[nodiscard]] std::int32_t i32(std::string_view label) {
    return std::bit_cast<std::int32_t>(u32(label));
  }

  [[nodiscard]] std::uint64_t u64(std::string_view label) {
    require(sizeof(std::uint64_t), label);
    std::uint64_t value = 0U;
    for (unsigned shift = 0; shift < 64U; shift += 8U) {
      value |= static_cast<std::uint64_t>(m_bytes[m_offset++]) << shift;
    }
    return value;
  }

  [[nodiscard]] double f64(std::string_view label) {
    return std::bit_cast<double>(u64(label));
  }

  [[nodiscard]] std::string string(std::string_view label) {
    const std::uint64_t size64 = u64(label);
    const std::size_t size = core::checkedIntegralNarrow<std::size_t>(size64, label);
    require(size, label);
    std::string result(reinterpret_cast<const char*>(m_bytes.data() + m_offset), size);
    m_offset += size;
    return result;
  }

  [[nodiscard]] std::vector<std::byte> bytePayload(std::string_view label) {
    const std::uint64_t size64 = u64(label);
    const std::size_t size = core::checkedIntegralNarrow<std::size_t>(size64, label);
    require(size, label);
    std::vector<std::byte> result(size);
    if (size != 0U) {
      std::memcpy(result.data(), m_bytes.data() + m_offset, size);
    }
    m_offset += size;
    return result;
  }

  [[nodiscard]] std::span<const std::uint8_t> remainingPayload(std::size_t size, std::string_view label) {
    require(size, label);
    const auto result = m_bytes.subspan(m_offset, size);
    m_offset += size;
    return result;
  }

  [[nodiscard]] bool atEnd() const noexcept { return m_offset == m_bytes.size(); }
  [[nodiscard]] std::size_t offset() const noexcept { return m_offset; }

 private:
  void require(std::size_t size, std::string_view label) const {
    if (size > m_bytes.size() - std::min(m_offset, m_bytes.size())) {
      throw std::runtime_error(
          "migration wire packet truncated while reading " + std::string(label));
    }
  }

  std::span<const std::uint8_t> m_bytes;
  std::size_t m_offset = 0U;
};

void appendRequirement(WireWriter& out, const core::ModuleSidecarRequirement& requirement) {
  out.u32(static_cast<std::uint32_t>(requirement.kind));
  out.u32(requirement.species_mask);
  out.u32(requirement.particle_flags_mask);
  out.f64(requirement.threshold_code);
}

core::ModuleSidecarRequirement readRequirement(WireReader& in) {
  core::ModuleSidecarRequirement requirement;
  requirement.kind = static_cast<core::ModuleSidecarRequirementKind>(in.u32("module_requirement_kind"));
  requirement.species_mask = in.u32("module_requirement_species_mask");
  requirement.particle_flags_mask = in.u32("module_requirement_particle_flags_mask");
  requirement.threshold_code = in.f64("module_requirement_threshold_code");
  return requirement;
}

void appendModulePayload(WireWriter& out, const core::ModuleParticleSidecarPayload& payload) {
  out.string(payload.module_name);
  out.u32(payload.schema_version);
  out.u32(payload.row_stride_bytes);
  out.u32(payload.required_species_mask);
  appendRequirement(out, payload.requirement);
  out.bytePayload(std::span<const std::byte>(payload.payload.data(), payload.payload.size()));
}

core::ModuleParticleSidecarPayload readModulePayload(WireReader& in) {
  core::ModuleParticleSidecarPayload payload;
  payload.module_name = in.string("module_name");
  payload.schema_version = in.u32("module_schema_version");
  payload.row_stride_bytes = in.u32("module_row_stride_bytes");
  payload.required_species_mask = in.u32("module_required_species_mask");
  payload.requirement = readRequirement(in);
  payload.payload = in.bytePayload("module_payload");
  return payload;
}

void appendSchedulerFields(WireWriter& out, const core::SchedulerMigrationFields& fields) {
  out.u8(fields.bin_index);
  out.u64(fields.next_activation_tick);
  out.u8(fields.pending_bin_index);
}

core::SchedulerMigrationFields readSchedulerFields(WireReader& in) {
  return core::SchedulerMigrationFields{
      .bin_index = in.u8("scheduler_bin_index"),
      .next_activation_tick = in.u64("scheduler_next_activation_tick"),
      .pending_bin_index = in.u8("scheduler_pending_bin_index"),
  };
}

void appendGasFields(WireWriter& out, const core::GasCellMigrationFields& fields) {
  out.u64(fields.gas_cell_id);
  out.u8(fields.has_parent_particle);
  out.u64(fields.parent_particle_id);
  out.u64(fields.owning_patch_id);
  out.u32(fields.destination_local_cell_row);
  out.u64(fields.gas_cell_identity_generation);
  out.u64(fields.ghost_hydro_epoch);
  out.f64(fields.center_x_comoving);
  out.f64(fields.center_y_comoving);
  out.f64(fields.center_z_comoving);
  out.f64(fields.cell_mass_code);
  out.u8(fields.cell_time_bin);
  out.u32(fields.patch_index);
  out.f64(fields.velocity_x_peculiar);
  out.f64(fields.velocity_y_peculiar);
  out.f64(fields.velocity_z_peculiar);
  out.f64(fields.density_code);
  out.f64(fields.pressure_code);
  out.f64(fields.internal_energy_code);
  out.f64(fields.metal_mass_code);
  out.f64(fields.temperature_code);
  out.f64(fields.sound_speed_code);
  // CHUI Part 2 (particle wire v3 / AMR patch wire v4): composition payload.
  out.u32(static_cast<std::uint32_t>(fields.composition_payload.size()));
  for (const double value : fields.composition_payload) out.f64(value);
}

core::GasCellMigrationFields readGasFields(WireReader& in) {
  core::GasCellMigrationFields fields;
  fields.gas_cell_id = in.u64("gas_cell_id");
  fields.has_parent_particle = in.u8("gas_has_parent_particle");
  fields.parent_particle_id = in.u64("gas_parent_particle_id");
  fields.owning_patch_id = in.u64("gas_owning_patch_id");
  fields.destination_local_cell_row = in.u32("gas_destination_local_cell_row");
  fields.gas_cell_identity_generation = in.u64("gas_cell_identity_generation");
  fields.ghost_hydro_epoch = in.u64("gas_ghost_hydro_epoch");
  fields.center_x_comoving = in.f64("gas_center_x_comoving");
  fields.center_y_comoving = in.f64("gas_center_y_comoving");
  fields.center_z_comoving = in.f64("gas_center_z_comoving");
  fields.cell_mass_code = in.f64("gas_cell_mass_code");
  fields.cell_time_bin = in.u8("gas_cell_time_bin");
  fields.patch_index = in.u32("gas_patch_index");
  fields.velocity_x_peculiar = in.f64("gas_velocity_x_peculiar");
  fields.velocity_y_peculiar = in.f64("gas_velocity_y_peculiar");
  fields.velocity_z_peculiar = in.f64("gas_velocity_z_peculiar");
  fields.density_code = in.f64("gas_density_code");
  fields.pressure_code = in.f64("gas_pressure_code");
  fields.internal_energy_code = in.f64("gas_internal_energy_code");
  fields.metal_mass_code = in.f64("gas_metal_mass_code");
  fields.temperature_code = in.f64("gas_temperature_code");
  fields.sound_speed_code = in.f64("gas_sound_speed_code");
  const std::uint32_t composition_width = in.u32("gas_composition_payload_width");
  // Bounded: the packed width is lanes + metadata of one manifest (< 4096).
  if (composition_width > 4096U) {
    throw std::runtime_error("gas composition payload width exceeds the wire bound");
  }
  fields.composition_payload.resize(composition_width);
  for (std::uint32_t i = 0; i < composition_width; ++i) {
    fields.composition_payload[i] = in.f64("gas_composition_payload_value");
  }
  return fields;
}

void appendStarFields(WireWriter& out, const core::StarParticleMigrationFields& fields) {
  out.f64(fields.formation_scale_factor);
  out.f64(fields.birth_mass_code);
  out.f64(fields.metallicity_mass_fraction);
  out.u64(fields.birth_key);
  out.u64(fields.parent_gas_cell_id);
  out.u64(fields.birth_tick);
  out.u32(fields.birth_ordinal);
  out.f64(fields.stellar_age_years_last);
  out.f64(fields.stellar_returned_mass_cumulative_code);
  out.f64(fields.stellar_returned_metals_cumulative_code);
  out.f64(fields.stellar_newly_synthesized_metals_cumulative_code);
  out.f64(fields.stellar_feedback_energy_cumulative_erg);
  out.f64(fields.enrichment_carry_mass_code);
  out.f64(fields.enrichment_carry_metals_code);
  out.f64(fields.enrichment_carry_feedback_energy_erg);
  out.f64(fields.enrichment_carry_momentum_code);
  out.f64(fields.stellar_deposited_mass_cumulative_code);
  out.f64(fields.stellar_deposited_metals_cumulative_code);
  out.f64(fields.stellar_deposited_feedback_energy_cumulative_erg);
  for (double value : fields.stellar_returned_mass_channel_cumulative_code) out.f64(value);
  for (double value : fields.stellar_returned_metals_channel_cumulative_code) out.f64(value);
  for (double value : fields.stellar_feedback_energy_channel_cumulative_erg) out.f64(value);
}

core::StarParticleMigrationFields readStarFields(WireReader& in) {
  core::StarParticleMigrationFields fields;
  fields.formation_scale_factor = in.f64("star_formation_scale_factor");
  fields.birth_mass_code = in.f64("star_birth_mass_code");
  fields.metallicity_mass_fraction = in.f64("star_metallicity_mass_fraction");
  fields.birth_key = in.u64("star_birth_key");
  fields.parent_gas_cell_id = in.u64("star_parent_gas_cell_id");
  fields.birth_tick = in.u64("star_birth_tick");
  fields.birth_ordinal = in.u32("star_birth_ordinal");
  fields.stellar_age_years_last = in.f64("star_age_years_last");
  fields.stellar_returned_mass_cumulative_code = in.f64("star_returned_mass_cumulative_code");
  fields.stellar_returned_metals_cumulative_code = in.f64("star_returned_metals_cumulative_code");
  fields.stellar_newly_synthesized_metals_cumulative_code = in.f64("star_new_metals_cumulative_code");
  fields.stellar_feedback_energy_cumulative_erg = in.f64("star_feedback_energy_cumulative_erg");
  fields.enrichment_carry_mass_code = in.f64("star_enrichment_carry_mass_code");
  fields.enrichment_carry_metals_code = in.f64("star_enrichment_carry_metals_code");
  fields.enrichment_carry_feedback_energy_erg = in.f64("star_enrichment_carry_feedback_energy_erg");
  fields.enrichment_carry_momentum_code = in.f64("star_enrichment_carry_momentum_code");
  fields.stellar_deposited_mass_cumulative_code = in.f64("star_deposited_mass_cumulative_code");
  fields.stellar_deposited_metals_cumulative_code = in.f64("star_deposited_metals_cumulative_code");
  fields.stellar_deposited_feedback_energy_cumulative_erg = in.f64("star_deposited_feedback_energy_cumulative_erg");
  for (double& value : fields.stellar_returned_mass_channel_cumulative_code) value = in.f64("star_returned_mass_channel");
  for (double& value : fields.stellar_returned_metals_channel_cumulative_code) value = in.f64("star_returned_metals_channel");
  for (double& value : fields.stellar_feedback_energy_channel_cumulative_erg) value = in.f64("star_feedback_energy_channel");
  return fields;
}

void appendBlackHoleFields(WireWriter& out, const core::BlackHoleParticleMigrationFields& fields) {
  out.u32(fields.host_cell_index);
  out.f64(fields.subgrid_mass_code);
  out.f64(fields.accretion_rate_code);
  out.f64(fields.feedback_energy_code);
  out.f64(fields.eddington_ratio);
  out.f64(fields.cumulative_accreted_mass_code);
  out.f64(fields.cumulative_feedback_energy_code);
  out.f64(fields.duty_cycle_active_time_code);
  out.f64(fields.duty_cycle_total_time_code);
}

core::BlackHoleParticleMigrationFields readBlackHoleFields(WireReader& in) {
  core::BlackHoleParticleMigrationFields fields;
  fields.host_cell_index = in.u32("bh_host_cell_index");
  fields.subgrid_mass_code = in.f64("bh_subgrid_mass_code");
  fields.accretion_rate_code = in.f64("bh_accretion_rate_code");
  fields.feedback_energy_code = in.f64("bh_feedback_energy_code");
  fields.eddington_ratio = in.f64("bh_eddington_ratio");
  fields.cumulative_accreted_mass_code = in.f64("bh_cumulative_accreted_mass_code");
  fields.cumulative_feedback_energy_code = in.f64("bh_cumulative_feedback_energy_code");
  fields.duty_cycle_active_time_code = in.f64("bh_duty_cycle_active_time_code");
  fields.duty_cycle_total_time_code = in.f64("bh_duty_cycle_total_time_code");
  return fields;
}

void appendTracerFields(WireWriter& out, const core::TracerParticleMigrationFields& fields) {
  out.u64(fields.parent_particle_id);
  out.u64(fields.injection_step);
  out.u32(fields.host_cell_index);
  out.f64(fields.mass_fraction_of_host);
  out.f64(fields.last_host_mass_code);
  out.f64(fields.cumulative_exchanged_mass_code);
}

core::TracerParticleMigrationFields readTracerFields(WireReader& in) {
  core::TracerParticleMigrationFields fields;
  fields.parent_particle_id = in.u64("tracer_parent_particle_id");
  fields.injection_step = in.u64("tracer_injection_step");
  fields.host_cell_index = in.u32("tracer_host_cell_index");
  fields.mass_fraction_of_host = in.f64("tracer_mass_fraction_of_host");
  fields.last_host_mass_code = in.f64("tracer_last_host_mass_code");
  fields.cumulative_exchanged_mass_code = in.f64("tracer_cumulative_exchanged_mass_code");
  return fields;
}

void appendAmrPatchFields(WireWriter& out, const core::AmrPatchMigrationFields& fields) {
  out.u64(fields.patch_id);
  out.u64(fields.parent_patch_id);
  out.u64(fields.morton_key);
  out.i32(fields.level);
  out.u32(fields.owning_rank);
  out.u32(fields.cell_count);
  out.f64(fields.origin_x_comoving);
  out.f64(fields.origin_y_comoving);
  out.f64(fields.origin_z_comoving);
  out.f64(fields.extent_x_comoving);
  out.f64(fields.extent_y_comoving);
  out.f64(fields.extent_z_comoving);
  out.u16(fields.cell_dim_x);
  out.u16(fields.cell_dim_y);
  out.u16(fields.cell_dim_z);
}

core::AmrPatchMigrationFields readAmrPatchFields(WireReader& in) {
  core::AmrPatchMigrationFields fields;
  fields.patch_id = in.u64("patch_id");
  fields.parent_patch_id = in.u64("patch_parent_id");
  fields.morton_key = in.u64("patch_morton_key");
  fields.level = in.i32("patch_level");
  fields.owning_rank = in.u32("patch_owning_rank");
  fields.cell_count = in.u32("patch_cell_count");
  fields.origin_x_comoving = in.f64("patch_origin_x_comoving");
  fields.origin_y_comoving = in.f64("patch_origin_y_comoving");
  fields.origin_z_comoving = in.f64("patch_origin_z_comoving");
  fields.extent_x_comoving = in.f64("patch_extent_x_comoving");
  fields.extent_y_comoving = in.f64("patch_extent_y_comoving");
  fields.extent_z_comoving = in.f64("patch_extent_z_comoving");
  fields.cell_dim_x = in.u16("patch_cell_dim_x");
  fields.cell_dim_y = in.u16("patch_cell_dim_y");
  fields.cell_dim_z = in.u16("patch_cell_dim_z");
  return fields;
}

void appendGasSchedulerRecord(WireWriter& out, const core::GasCellSchedulerMigrationRecord& record) {
  out.u64(record.gas_cell_id);
  out.u8(record.bin_index);
  out.u64(record.next_activation_tick);
  out.u8(record.pending_bin_index);
}

core::GasCellSchedulerMigrationRecord readGasSchedulerRecord(WireReader& in) {
  return core::GasCellSchedulerMigrationRecord{
      .gas_cell_id = in.u64("gas_scheduler_cell_id"),
      .bin_index = in.u8("gas_scheduler_bin_index"),
      .next_activation_tick = in.u64("gas_scheduler_next_activation_tick"),
      .pending_bin_index = in.u8("gas_scheduler_pending_bin_index"),
  };
}


void appendPendingFluxRegisterRecord(
    WireWriter& out,
    const core::PendingFluxRegisterRecord& record) {
  out.u64(record.register_key);
  out.u64(record.coarse_patch_id);
  out.u64(record.coarse_gas_cell_id);
  out.u64(core::checkedIntegralNarrow<std::uint64_t>(
      record.coarse_cell_index, "pending_flux_coarse_cell_index"));
  out.u8(record.level);
  out.u8(record.axis);
  out.u8(record.orientation);
  out.f64(record.expected_area_comov);
  out.f64(record.coarse_area_accumulated_comov);
  out.f64(record.fine_area_accumulated_comov);
  out.f64(record.interval_start_code);
  out.f64(record.interval_end_code);
  out.f64(record.coarse_dt_code);
  out.u32(record.expected_fine_substeps);
  out.u32(record.completed_fine_substeps);
  out.u64(record.fine_substep_coverage_mask);
  out.u32(record.coarse_face_count);
  out.u32(record.fine_face_count);
  out.u64(record.gas_cell_identity_generation);
  out.u64(record.patch_geometry_generation);
  out.f64(record.coarse_mass_flux_integral_code);
  out.f64(record.coarse_momentum_x_flux_integral_code);
  out.f64(record.coarse_momentum_y_flux_integral_code);
  out.f64(record.coarse_momentum_z_flux_integral_code);
  out.f64(record.coarse_total_energy_flux_integral_code);
  out.f64(record.coarse_metal_mass_flux_integral_code);
  out.f64(record.fine_mass_flux_integral_code);
  out.f64(record.fine_momentum_x_flux_integral_code);
  out.f64(record.fine_momentum_y_flux_integral_code);
  out.f64(record.fine_momentum_z_flux_integral_code);
  out.f64(record.fine_total_energy_flux_integral_code);
  out.f64(record.fine_metal_mass_flux_integral_code);
}

core::PendingFluxRegisterRecord readPendingFluxRegisterRecord(WireReader& in) {
  core::PendingFluxRegisterRecord record;
  record.register_key = in.u64("pending_flux_register_key");
  record.coarse_patch_id = in.u64("pending_flux_coarse_patch_id");
  record.coarse_gas_cell_id = in.u64("pending_flux_coarse_gas_cell_id");
  record.coarse_cell_index = core::checkedIntegralNarrow<std::size_t>(
      in.u64("pending_flux_coarse_cell_index"), "pending_flux_coarse_cell_index");
  record.level = in.u8("pending_flux_level");
  record.axis = in.u8("pending_flux_axis");
  record.orientation = in.u8("pending_flux_orientation");
  record.expected_area_comov = in.f64("pending_flux_expected_area");
  record.coarse_area_accumulated_comov = in.f64("pending_flux_coarse_area");
  record.fine_area_accumulated_comov = in.f64("pending_flux_fine_area");
  record.interval_start_code = in.f64("pending_flux_interval_start");
  record.interval_end_code = in.f64("pending_flux_interval_end");
  record.coarse_dt_code = in.f64("pending_flux_coarse_dt");
  record.expected_fine_substeps = in.u32("pending_flux_expected_fine_substeps");
  record.completed_fine_substeps = in.u32("pending_flux_completed_fine_substeps");
  record.fine_substep_coverage_mask = in.u64("pending_flux_coverage_mask");
  record.coarse_face_count = in.u32("pending_flux_coarse_face_count");
  record.fine_face_count = in.u32("pending_flux_fine_face_count");
  record.gas_cell_identity_generation = in.u64("pending_flux_identity_generation");
  record.patch_geometry_generation = in.u64("pending_flux_geometry_generation");
  record.coarse_mass_flux_integral_code = in.f64("pending_flux_coarse_mass");
  record.coarse_momentum_x_flux_integral_code = in.f64("pending_flux_coarse_momentum_x");
  record.coarse_momentum_y_flux_integral_code = in.f64("pending_flux_coarse_momentum_y");
  record.coarse_momentum_z_flux_integral_code = in.f64("pending_flux_coarse_momentum_z");
  record.coarse_total_energy_flux_integral_code = in.f64("pending_flux_coarse_energy");
  record.coarse_metal_mass_flux_integral_code = in.f64("pending_flux_coarse_metal");
  record.fine_mass_flux_integral_code = in.f64("pending_flux_fine_mass");
  record.fine_momentum_x_flux_integral_code = in.f64("pending_flux_fine_momentum_x");
  record.fine_momentum_y_flux_integral_code = in.f64("pending_flux_fine_momentum_y");
  record.fine_momentum_z_flux_integral_code = in.f64("pending_flux_fine_momentum_z");
  record.fine_total_energy_flux_integral_code = in.f64("pending_flux_fine_energy");
  record.fine_metal_mass_flux_integral_code = in.f64("pending_flux_fine_metal");
  return record;
}

void appendTemporalHistoryCell(
    WireWriter& out,
    const core::AmrTemporalBoundaryHistoryCellRecord& cell) {
  out.u64(cell.gas_cell_id);
  out.u64(core::checkedIntegralNarrow<std::uint64_t>(
      cell.patch_local_cell, "temporal_history_patch_local_cell"));
  out.f64(cell.start_mass_density_comoving);
  out.f64(cell.start_momentum_density_x_comoving);
  out.f64(cell.start_momentum_density_y_comoving);
  out.f64(cell.start_momentum_density_z_comoving);
  out.f64(cell.start_total_energy_density_comoving);
  out.f64(cell.start_metal_mass_density_comoving);
  out.f64(cell.end_mass_density_comoving);
  out.f64(cell.end_momentum_density_x_comoving);
  out.f64(cell.end_momentum_density_y_comoving);
  out.f64(cell.end_momentum_density_z_comoving);
  out.f64(cell.end_total_energy_density_comoving);
  out.f64(cell.end_metal_mass_density_comoving);
}

core::AmrTemporalBoundaryHistoryCellRecord readTemporalHistoryCell(WireReader& in) {
  core::AmrTemporalBoundaryHistoryCellRecord cell;
  cell.gas_cell_id = in.u64("temporal_history_gas_cell_id");
  cell.patch_local_cell = core::checkedIntegralNarrow<std::size_t>(
      in.u64("temporal_history_patch_local_cell"), "temporal_history_patch_local_cell");
  cell.start_mass_density_comoving = in.f64("temporal_history_start_mass");
  cell.start_momentum_density_x_comoving = in.f64("temporal_history_start_momentum_x");
  cell.start_momentum_density_y_comoving = in.f64("temporal_history_start_momentum_y");
  cell.start_momentum_density_z_comoving = in.f64("temporal_history_start_momentum_z");
  cell.start_total_energy_density_comoving = in.f64("temporal_history_start_energy");
  cell.start_metal_mass_density_comoving = in.f64("temporal_history_start_metal");
  cell.end_mass_density_comoving = in.f64("temporal_history_end_mass");
  cell.end_momentum_density_x_comoving = in.f64("temporal_history_end_momentum_x");
  cell.end_momentum_density_y_comoving = in.f64("temporal_history_end_momentum_y");
  cell.end_momentum_density_z_comoving = in.f64("temporal_history_end_momentum_z");
  cell.end_total_energy_density_comoving = in.f64("temporal_history_end_energy");
  cell.end_metal_mass_density_comoving = in.f64("temporal_history_end_metal");
  return cell;
}

void appendTemporalHistoryRecord(
    WireWriter& out,
    const core::AmrTemporalBoundaryHistoryRecord& record) {
  out.u64(record.patch_id);
  out.u8(record.patch_level);
  out.u64(record.patch_geometry_fingerprint);
  out.u64(record.gas_cell_identity_generation);
  out.f64(record.interval_start_code);
  out.f64(record.interval_end_code);
  out.u8(record.end_state_valid ? 1U : 0U);
  out.u64(static_cast<std::uint64_t>(record.cells.size()));
  for (const auto& cell : record.cells) {
    appendTemporalHistoryCell(out, cell);
  }
}

core::AmrTemporalBoundaryHistoryRecord readTemporalHistoryRecord(WireReader& in) {
  core::AmrTemporalBoundaryHistoryRecord record;
  record.patch_id = in.u64("temporal_history_patch_id");
  record.patch_level = in.u8("temporal_history_patch_level");
  record.patch_geometry_fingerprint = in.u64("temporal_history_geometry_fingerprint");
  record.gas_cell_identity_generation = in.u64("temporal_history_identity_generation");
  record.interval_start_code = in.f64("temporal_history_interval_start");
  record.interval_end_code = in.f64("temporal_history_interval_end");
  const std::uint8_t end_state_valid = in.u8("temporal_history_end_state_valid");
  if (end_state_valid > 1U) {
    throw std::runtime_error("AMR temporal history wire boolean is not canonical");
  }
  record.end_state_valid = end_state_valid != 0U;
  const std::uint64_t cell_count = in.u64("temporal_history_cell_count");
  if (cell_count > 100'000'000ULL) {
    throw std::runtime_error("AMR temporal history wire has unreasonable cell count");
  }
  record.cells.reserve(core::checkedIntegralNarrow<std::size_t>(
      cell_count, "temporal_history_cell_count"));
  for (std::uint64_t i = 0U; i < cell_count; ++i) {
    record.cells.push_back(readTemporalHistoryCell(in));
  }
  return record;
}

[[nodiscard]] std::size_t moduleWireUpperBound(const core::ModuleSidecarBlock& block) {
  constexpr std::size_t k_requirement_bytes = sizeof(std::uint32_t) * 3U + sizeof(double);
  const std::size_t metadata = core::checkedSizeAdd(
      sizeof(std::uint64_t), block.module_name.size(), "migration module name wire bytes");
  const std::size_t fixed = sizeof(std::uint32_t) * 3U + k_requirement_bytes + sizeof(std::uint64_t);
  return core::checkedSizeAdd(
      core::checkedSizeAdd(metadata, fixed, "migration module metadata wire bytes"),
      block.row_stride_bytes,
      "migration module payload wire bytes");
}

[[nodiscard]] std::size_t encodedSize(const core::ParticleMigrationRecord& record) {
  return encodeParticleMigrationRecord(record).size();
}

}  // namespace

PacketCapacityPlan planPacketCapacity(
    std::size_t transport_round_limit_bytes,
    std::size_t rank_count) {
  if (transport_round_limit_bytes == 0U || rank_count == 0U) {
    throw std::invalid_argument("migration packet capacity requires positive round/rank extents");
  }
  const std::size_t per_peer_packet_bytes = transport_round_limit_bytes / rank_count;
  if (per_peer_packet_bytes <= k_fragment_header_bytes) {
    throw std::invalid_argument(
        "migration transport round limit is too small for one fragment header per rank");
  }
  return PacketCapacityPlan{
      .per_peer_packet_bytes = per_peer_packet_bytes,
      .fragment_payload_bytes = per_peer_packet_bytes - k_fragment_header_bytes,
      .aggregate_packet_bytes = core::checkedSizeMultiply(
          per_peer_packet_bytes,
          rank_count,
          "migration aggregate packet capacity"),
  };
}

std::vector<std::uint8_t> encodeDmoParticleMigrationRecord(
    const DmoParticleMigrationRecord& record) {
  WireWriter out;
  out.u32(k_dmo_particle_record_wire_version);
  out.u64(record.particle_id);
  out.u64(record.sfc_key);
  out.f64(record.position_x_comoving);
  out.f64(record.position_y_comoving);
  out.f64(record.position_z_comoving);
  out.f64(record.velocity_x_peculiar);
  out.f64(record.velocity_y_peculiar);
  out.f64(record.velocity_z_peculiar);
  out.f64(record.mass_code);
  auto bytes = std::move(out).take();
  if (bytes.size() != dmoParticleMigrationWireBytes()) {
    throw std::logic_error("compact DMO migration wire size contract mismatch");
  }
  return bytes;
}

DmoParticleMigrationRecord decodeDmoParticleMigrationRecord(
    std::span<const std::uint8_t> bytes) {
  if (bytes.size() != dmoParticleMigrationWireBytes()) {
    throw std::runtime_error("compact DMO migration wire extent mismatch");
  }
  WireReader in(bytes);
  if (in.u32("dmo_particle_wire_version") != k_dmo_particle_record_wire_version) {
    throw std::runtime_error("unsupported compact DMO migration wire version");
  }
  DmoParticleMigrationRecord record;
  record.particle_id = in.u64("dmo_particle_id");
  record.sfc_key = in.u64("dmo_sfc_key");
  record.position_x_comoving = in.f64("dmo_position_x_comoving");
  record.position_y_comoving = in.f64("dmo_position_y_comoving");
  record.position_z_comoving = in.f64("dmo_position_z_comoving");
  record.velocity_x_peculiar = in.f64("dmo_velocity_x_peculiar");
  record.velocity_y_peculiar = in.f64("dmo_velocity_y_peculiar");
  record.velocity_z_peculiar = in.f64("dmo_velocity_z_peculiar");
  record.mass_code = in.f64("dmo_mass_code");
  if (!in.atEnd()) throw std::runtime_error("compact DMO migration wire has trailing bytes");
  return record;
}

std::vector<std::uint8_t> encodeParticleMigrationRecord(
    const core::ParticleMigrationRecord& record) {
  WireWriter out;
  out.u32(k_particle_record_wire_version);
  out.u64(record.particle_id);
  out.u64(record.sfc_key);
  out.u32(record.species_tag);
  out.u32(record.particle_flags);
  out.u32(record.owning_rank);
  out.f64(record.position_x_comoving);
  out.f64(record.position_y_comoving);
  out.f64(record.position_z_comoving);
  out.f64(record.velocity_x_peculiar);
  out.f64(record.velocity_y_peculiar);
  out.f64(record.velocity_z_peculiar);
  out.f64(record.mass_code);
  out.u8(record.time_bin);

  std::uint8_t common_flags = 0U;
  if (record.has_scheduler_fields) common_flags |= 1U << 0U;
  if (record.has_gravity_softening_value) common_flags |= 1U << 1U;
  if (record.has_gravity_softening_override) common_flags |= 1U << 2U;
  out.u8(common_flags);
  if (record.has_scheduler_fields) appendSchedulerFields(out, record.scheduler_fields);
  out.f64(record.last_drift_time_code);
  out.f64(record.last_drift_scale_factor);
  if (record.has_gravity_softening_value) out.f64(record.gravity_softening_comoving);

  std::uint8_t species_flags = 0U;
  if (record.has_gas_cell_fields) species_flags |= 1U << 0U;
  if (record.has_star_fields) species_flags |= 1U << 1U;
  if (record.has_black_hole_fields) species_flags |= 1U << 2U;
  if (record.has_tracer_fields) species_flags |= 1U << 3U;
  out.u8(species_flags);
  if (record.has_gas_cell_fields) appendGasFields(out, record.gas_cell_fields);
  if (record.has_star_fields) appendStarFields(out, record.star_fields);
  if (record.has_black_hole_fields) appendBlackHoleFields(out, record.black_hole_fields);
  if (record.has_tracer_fields) appendTracerFields(out, record.tracer_fields);

  out.u64(static_cast<std::uint64_t>(record.module_sidecar_payloads.size()));
  for (const auto& payload : record.module_sidecar_payloads) appendModulePayload(out, payload);
  return std::move(out).take();
}

core::ParticleMigrationRecord decodeParticleMigrationRecord(
    std::span<const std::uint8_t> bytes) {
  WireReader in(bytes);
  if (in.u32("particle_wire_version") != k_particle_record_wire_version) {
    throw std::runtime_error("unsupported particle migration wire version");
  }
  core::ParticleMigrationRecord record;
  record.particle_id = in.u64("particle_id");
  record.sfc_key = in.u64("sfc_key");
  record.species_tag = in.u32("species_tag");
  record.particle_flags = in.u32("particle_flags");
  record.owning_rank = in.u32("owning_rank");
  record.position_x_comoving = in.f64("position_x_comoving");
  record.position_y_comoving = in.f64("position_y_comoving");
  record.position_z_comoving = in.f64("position_z_comoving");
  record.velocity_x_peculiar = in.f64("velocity_x_peculiar");
  record.velocity_y_peculiar = in.f64("velocity_y_peculiar");
  record.velocity_z_peculiar = in.f64("velocity_z_peculiar");
  record.mass_code = in.f64("mass_code");
  record.time_bin = in.u8("time_bin");
  const std::uint8_t common_flags = in.u8("common_flags");
  if ((common_flags & 0xf8U) != 0U) throw std::runtime_error("particle migration wire common flags contain unknown bits");
  record.has_scheduler_fields = (common_flags & (1U << 0U)) != 0U;
  record.has_gravity_softening_value = (common_flags & (1U << 1U)) != 0U;
  record.has_gravity_softening_override = (common_flags & (1U << 2U)) != 0U;
  if (record.has_gravity_softening_override && !record.has_gravity_softening_value) {
    throw std::runtime_error("particle migration wire softening override lacks a softening value");
  }
  if (record.has_scheduler_fields) record.scheduler_fields = readSchedulerFields(in);
  record.last_drift_time_code = in.f64("last_drift_time_code");
  record.last_drift_scale_factor = in.f64("last_drift_scale_factor");
  if (record.has_gravity_softening_value) record.gravity_softening_comoving = in.f64("gravity_softening_comoving");

  const std::uint8_t species_flags = in.u8("species_flags");
  if ((species_flags & 0xf0U) != 0U) throw std::runtime_error("particle migration wire species flags contain unknown bits");
  record.has_gas_cell_fields = (species_flags & (1U << 0U)) != 0U;
  record.has_star_fields = (species_flags & (1U << 1U)) != 0U;
  record.has_black_hole_fields = (species_flags & (1U << 2U)) != 0U;
  record.has_tracer_fields = (species_flags & (1U << 3U)) != 0U;
  if (record.has_gas_cell_fields) record.gas_cell_fields = readGasFields(in);
  if (record.has_star_fields) record.star_fields = readStarFields(in);
  if (record.has_black_hole_fields) record.black_hole_fields = readBlackHoleFields(in);
  if (record.has_tracer_fields) record.tracer_fields = readTracerFields(in);

  const std::uint64_t module_count = in.u64("module_count");
  if (module_count > 1'000'000ULL) {
    throw std::runtime_error("particle migration wire packet has unreasonable module payload count");
  }
  record.module_sidecar_payloads.reserve(core::checkedIntegralNarrow<std::size_t>(module_count, "module_count"));
  for (std::uint64_t i = 0U; i < module_count; ++i) {
    record.module_sidecar_payloads.push_back(readModulePayload(in));
  }
  if (!in.atEnd()) throw std::runtime_error("particle migration wire record has trailing bytes");
  return record;
}

std::vector<std::uint8_t> encodeAmrPatchMigrationRecord(
    const core::AmrPatchMigrationRecord& record) {
  WireWriter out;
  out.u32(k_amr_patch_record_wire_version);
  appendAmrPatchFields(out, record.patch);
  out.u64(static_cast<std::uint64_t>(record.gas_cell_records.size()));
  for (const auto& gas_record : record.gas_cell_records) {
    out.u32(gas_record.owning_rank);
    appendGasFields(out, gas_record.fields);
  }
  out.u64(static_cast<std::uint64_t>(record.gas_cell_scheduler_records.size()));
  for (const auto& scheduler_record : record.gas_cell_scheduler_records) {
    appendGasSchedulerRecord(out, scheduler_record);
  }
  out.u64(static_cast<std::uint64_t>(record.pending_flux_register_records.size()));
  for (const auto& pending : record.pending_flux_register_records) {
    appendPendingFluxRegisterRecord(out, pending);
  }
  out.u64(static_cast<std::uint64_t>(record.temporal_boundary_history_records.size()));
  for (const auto& history : record.temporal_boundary_history_records) {
    appendTemporalHistoryRecord(out, history);
  }
  return std::move(out).take();
}

core::AmrPatchMigrationRecord decodeAmrPatchMigrationRecord(
    std::span<const std::uint8_t> bytes) {
  WireReader in(bytes);
  if (in.u32("amr_patch_wire_version") != k_amr_patch_record_wire_version) {
    throw std::runtime_error("unsupported AMR patch migration wire version");
  }
  core::AmrPatchMigrationRecord record;
  record.patch = readAmrPatchFields(in);
  const std::uint64_t gas_count = in.u64("amr_patch_gas_cell_count");
  if (gas_count > 100'000'000ULL) throw std::runtime_error("AMR patch migration wire has unreasonable gas-cell count");
  record.gas_cell_records.reserve(core::checkedIntegralNarrow<std::size_t>(gas_count, "amr_patch_gas_cell_count"));
  for (std::uint64_t i = 0U; i < gas_count; ++i) {
    core::GasCellMigrationRecord gas_record;
    gas_record.owning_rank = in.u32("gas_cell_owning_rank");
    gas_record.fields = readGasFields(in);
    record.gas_cell_records.push_back(std::move(gas_record));
  }
  const std::uint64_t scheduler_count = in.u64("amr_patch_scheduler_record_count");
  if (scheduler_count != gas_count) {
    throw std::runtime_error("AMR patch migration wire scheduler record count does not match gas-cell count");
  }
  record.gas_cell_scheduler_records.reserve(core::checkedIntegralNarrow<std::size_t>(scheduler_count, "amr_patch_scheduler_record_count"));
  for (std::uint64_t i = 0U; i < scheduler_count; ++i) {
    record.gas_cell_scheduler_records.push_back(readGasSchedulerRecord(in));
  }
  const std::uint64_t pending_count = in.u64("amr_patch_pending_flux_count");
  if (pending_count > 100'000'000ULL) {
    throw std::runtime_error("AMR patch migration wire has unreasonable pending flux-register count");
  }
  record.pending_flux_register_records.reserve(core::checkedIntegralNarrow<std::size_t>(
      pending_count, "amr_patch_pending_flux_count"));
  for (std::uint64_t i = 0U; i < pending_count; ++i) {
    record.pending_flux_register_records.push_back(readPendingFluxRegisterRecord(in));
  }
  const std::uint64_t history_count = in.u64("amr_patch_temporal_history_count");
  if (history_count > 1'000'000ULL) {
    throw std::runtime_error("AMR patch migration wire has unreasonable temporal-history count");
  }
  record.temporal_boundary_history_records.reserve(core::checkedIntegralNarrow<std::size_t>(
      history_count, "amr_patch_temporal_history_count"));
  for (std::uint64_t i = 0U; i < history_count; ++i) {
    record.temporal_boundary_history_records.push_back(readTemporalHistoryRecord(in));
  }
  if (!in.atEnd()) throw std::runtime_error("AMR patch migration wire record has trailing bytes");
  return record;
}

std::vector<std::uint8_t> encodeFragment(
    std::uint64_t record_sequence,
    std::uint64_t record_total_bytes,
    std::uint64_t fragment_offset,
    std::span<const std::uint8_t> payload) {
  if (fragment_offset > record_total_bytes ||
      payload.size() > record_total_bytes - fragment_offset) {
    throw std::invalid_argument("migration fragment lies outside encoded record extent");
  }
  WireWriter out;
  out.u32(k_fragment_wire_version);
  out.u64(record_sequence);
  out.u64(record_total_bytes);
  out.u64(fragment_offset);
  out.u32(core::checkedIntegralNarrow<std::uint32_t>(payload.size(), "migration fragment payload bytes"));
  out.bytes(payload);
  return std::move(out).take();
}

FragmentView decodeFragment(std::span<const std::uint8_t> bytes) {
  WireReader in(bytes);
  if (in.u32("fragment_wire_version") != k_fragment_wire_version) {
    throw std::runtime_error("unsupported migration fragment wire version");
  }
  FragmentView fragment;
  fragment.record_sequence = in.u64("fragment_record_sequence");
  fragment.record_total_bytes = in.u64("fragment_record_total_bytes");
  fragment.fragment_offset = in.u64("fragment_offset");
  const std::size_t payload_bytes = in.u32("fragment_payload_bytes");
  if (fragment.fragment_offset > fragment.record_total_bytes ||
      payload_bytes > fragment.record_total_bytes - fragment.fragment_offset) {
    throw std::runtime_error("migration fragment metadata exceeds encoded record extent");
  }
  fragment.payload = in.remainingPayload(payload_bytes, "fragment_payload");
  if (!in.atEnd()) throw std::runtime_error("migration fragment has trailing bytes");
  return fragment;
}

std::size_t estimateParticleMigrationWireUpperBoundBytes(
    const core::SimulationState& state,
    std::uint32_t local_particle_index) {
  if (local_particle_index >= state.particles.size()) {
    throw std::out_of_range("particle migration wire estimate local index out of range");
  }
  core::ParticleMigrationRecord representative;
  representative.species_tag = state.particleSpeciesTag(local_particle_index);
  representative.has_scheduler_fields = true;
  representative.has_gravity_softening_value = !state.particle_sidecar.gravity_softening_comoving.empty();
  representative.has_gravity_softening_override =
      state.particle_sidecar.hasGravitySofteningOverride(local_particle_index);
  const auto species = static_cast<core::ParticleSpecies>(representative.species_tag);
  representative.has_gas_cell_fields = species == core::ParticleSpecies::kGas;
  representative.has_star_fields = species == core::ParticleSpecies::kStar;
  representative.has_black_hole_fields = species == core::ParticleSpecies::kBlackHole;
  representative.has_tracer_fields = species == core::ParticleSpecies::kTracer;
  std::size_t bytes = encodedSize(representative);
  for (const core::ModuleSidecarBlock* block : state.sidecars.blocksSortedByName()) {
    if (block == nullptr || !block->isParticleIndexed()) continue;
    bytes = core::checkedSizeAdd(bytes, moduleWireUpperBound(*block), "particle migration module wire upper bound");
  }
  return bytes;
}

std::size_t estimateAmrPatchMigrationWireUpperBoundBytes(
    const core::SimulationState& state,
    std::uint32_t local_patch_index) {
  if (local_patch_index >= state.patches.size()) {
    throw std::out_of_range("AMR migration wire estimate patch index out of range");
  }
  core::AmrPatchMigrationRecord representative;
  representative.gas_cell_records.resize(state.patches.cell_count[local_patch_index]);
  representative.gas_cell_scheduler_records.resize(state.patches.cell_count[local_patch_index]);
  std::size_t bytes = encodeAmrPatchMigrationRecord(representative).size();
  const std::uint64_t patch_id = state.patches.patch_id[local_patch_index];
  constexpr std::size_t k_pending_flux_wire_bytes = 219U;
  constexpr std::size_t k_temporal_history_fixed_wire_bytes = 50U;
  constexpr std::size_t k_temporal_history_cell_wire_bytes = 112U;
  for (const core::PendingFluxRegisterRecord& pending : state.pending_flux_registers.records()) {
    if (pending.coarse_patch_id == patch_id) {
      bytes = core::checkedSizeAdd(
          bytes, k_pending_flux_wire_bytes,
          "AMR migration pending flux-register wire upper bound");
    }
  }
  for (const core::AmrTemporalBoundaryHistoryRecord& history :
       state.amr_temporal_boundary_history.records()) {
    if (history.patch_id != patch_id) continue;
    const std::size_t cell_bytes = core::checkedSizeMultiply(
        history.cells.size(), k_temporal_history_cell_wire_bytes,
        "AMR migration temporal-history cell wire upper bound");
    bytes = core::checkedSizeAdd(
        bytes,
        core::checkedSizeAdd(
            k_temporal_history_fixed_wire_bytes, cell_bytes,
            "AMR migration temporal-history wire upper bound"),
        "AMR migration wire upper bound");
  }
  return bytes;
}


std::size_t estimateParticleMigrationDynamicHeapUpperBoundBytes(
    const core::SimulationState& state,
    std::uint32_t local_particle_index) {
  if (local_particle_index >= state.particles.size()) {
    throw std::out_of_range("particle migration heap estimate local index out of range");
  }
  std::size_t bytes = 0U;
  std::size_t module_count = 0U;
  for (const core::ModuleSidecarBlock* block : state.sidecars.blocksSortedByName()) {
    if (block == nullptr || !block->isParticleIndexed()) continue;
    ++module_count;
    bytes = core::checkedSizeAdd(
        bytes,
        core::checkedSizeAdd(
            block->module_name.size(), block->row_stride_bytes,
            "particle migration module dynamic bytes"),
        "particle migration dynamic heap bytes");
  }
  bytes = core::checkedSizeAdd(
      bytes,
      core::checkedSizeMultiply(
          module_count, sizeof(core::ModuleParticleSidecarPayload),
          "particle migration module payload vector bytes"),
      "particle migration dynamic heap bytes");
  return bytes;
}

std::size_t estimateAmrPatchMigrationDynamicHeapUpperBoundBytes(
    const core::SimulationState& state,
    std::uint32_t local_patch_index) {
  if (local_patch_index >= state.patches.size()) {
    throw std::out_of_range("AMR migration heap estimate patch index out of range");
  }
  const std::size_t cell_count = state.patches.cell_count[local_patch_index];
  std::size_t bytes = core::checkedSizeAdd(
      core::checkedSizeMultiply(
          cell_count, sizeof(core::GasCellMigrationRecord),
          "AMR migration gas-cell vector bytes"),
      core::checkedSizeMultiply(
          cell_count, sizeof(core::GasCellSchedulerMigrationRecord),
          "AMR migration scheduler vector bytes"),
      "AMR migration dynamic heap bytes");
  const std::uint64_t patch_id = state.patches.patch_id[local_patch_index];
  for (const core::PendingFluxRegisterRecord& pending : state.pending_flux_registers.records()) {
    if (pending.coarse_patch_id == patch_id) {
      bytes = core::checkedSizeAdd(
          bytes, sizeof(core::PendingFluxRegisterRecord),
          "AMR migration pending flux-register vector bytes");
    }
  }
  for (const core::AmrTemporalBoundaryHistoryRecord& history :
       state.amr_temporal_boundary_history.records()) {
    if (history.patch_id != patch_id) continue;
    bytes = core::checkedSizeAdd(
        bytes, sizeof(core::AmrTemporalBoundaryHistoryRecord),
        "AMR migration temporal-history record bytes");
    bytes = core::checkedSizeAdd(
        bytes,
        core::checkedSizeMultiply(
            history.cells.size(), sizeof(core::AmrTemporalBoundaryHistoryCellRecord),
            "AMR migration temporal-history cell vector bytes"),
        "AMR migration temporal-history dynamic heap bytes");
  }
  return bytes;
}

ParticleWireReferenceWidths particleWireReferenceWidths() {
  core::ParticleMigrationRecord dm;
  dm.species_tag = static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
  dm.has_scheduler_fields = true;
  const std::size_t dm_bytes = encodeParticleMigrationRecord(dm).size();

  core::ParticleMigrationRecord gas = dm;
  gas.species_tag = static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
  gas.has_gas_cell_fields = true;
  core::ParticleMigrationRecord star = dm;
  star.species_tag = static_cast<std::uint32_t>(core::ParticleSpecies::kStar);
  star.has_star_fields = true;
  core::ParticleMigrationRecord black_hole = dm;
  black_hole.species_tag = static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);
  black_hole.has_black_hole_fields = true;
  core::ParticleMigrationRecord tracer = dm;
  tracer.species_tag = static_cast<std::uint32_t>(core::ParticleSpecies::kTracer);
  tracer.has_tracer_fields = true;

  return ParticleWireReferenceWidths{
      .dark_matter_bytes = dm_bytes,
      .gas_extra_bytes = encodeParticleMigrationRecord(gas).size() - dm_bytes,
      .star_extra_bytes = encodeParticleMigrationRecord(star).size() - dm_bytes,
      .black_hole_extra_bytes = encodeParticleMigrationRecord(black_hole).size() - dm_bytes,
      .tracer_extra_bytes = encodeParticleMigrationRecord(tracer).size() - dm_bytes,
  };
}

}  // namespace cosmosim::workflows::internal::migration_wire
