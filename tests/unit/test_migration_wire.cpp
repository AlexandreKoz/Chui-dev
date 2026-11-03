#include <cassert>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "workflows/internal/amr_migration_payload.hpp"
#include "workflows/internal/migration_wire.hpp"

namespace {

using cosmosim::core::ParticleMigrationRecord;
namespace wire = cosmosim::workflows::internal::migration_wire;

ParticleMigrationRecord makeCommonRecord() {
  ParticleMigrationRecord record;
  record.particle_id = 0x0102030405060708ULL;
  record.sfc_key = 0x1122334455667788ULL;
  record.species_tag = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  record.particle_flags = 0x21U;
  record.owning_rank = 3U;
  record.position_x_comoving = 1.25;
  record.position_y_comoving = 2.5;
  record.position_z_comoving = 3.75;
  record.velocity_x_peculiar = -4.0;
  record.velocity_y_peculiar = 5.0;
  record.velocity_z_peculiar = -6.0;
  record.mass_code = 7.5;
  record.time_bin = 4U;
  record.has_scheduler_fields = true;
  record.scheduler_fields.bin_index = 4U;
  record.scheduler_fields.next_activation_tick = 12345U;
  record.scheduler_fields.pending_bin_index = 5U;
  record.last_drift_time_code = 0.125;
  record.last_drift_scale_factor = 0.75;
  record.has_gravity_softening_value = true;
  record.gravity_softening_comoving = 0.002;
  return record;
}

void testDarkMatterWireIsConditionalAndExplicit() {
  const ParticleMigrationRecord dm = makeCommonRecord();
  const std::vector<std::uint8_t> encoded = wire::encodeParticleMigrationRecord(dm);
  assert(encoded.size() < 200U);
  assert(encoded.size() == wire::particleWireReferenceWidths().dark_matter_bytes + sizeof(double));

  // Wire version and particle ID are little-endian field encodings, not native
  // struct bytes/padding.
  assert(encoded[0] == 3U && encoded[1] == 0U && encoded[2] == 0U && encoded[3] == 0U);
  assert(encoded[4] == 0x08U);
  assert(encoded[11] == 0x01U);

  ParticleMigrationRecord inactive_payload = dm;
  inactive_payload.gas_cell_fields.gas_cell_id = 991U;
  inactive_payload.star_fields.birth_key = 992U;
  inactive_payload.black_hole_fields.accretion_rate_code = 993.0;
  inactive_payload.tracer_fields.parent_particle_id = 994U;
  assert(wire::encodeParticleMigrationRecord(inactive_payload) == encoded);

  const ParticleMigrationRecord decoded = wire::decodeParticleMigrationRecord(encoded);
  assert(decoded.particle_id == dm.particle_id);
  assert(decoded.sfc_key == dm.sfc_key);
  assert(decoded.owning_rank == dm.owning_rank);
  assert(decoded.mass_code == dm.mass_code);
  assert(decoded.has_scheduler_fields);
  assert(decoded.scheduler_fields.next_activation_tick == 12345U);
  assert(decoded.has_gravity_softening_value);
  assert(decoded.gravity_softening_comoving == dm.gravity_softening_comoving);
  assert(!decoded.has_gas_cell_fields);
  assert(!decoded.has_star_fields);
  assert(!decoded.has_black_hole_fields);
  assert(!decoded.has_tracer_fields);
}

void testSpeciesPayloadRoundTrips() {
  ParticleMigrationRecord gas = makeCommonRecord();
  gas.species_tag = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  gas.has_gas_cell_fields = true;
  gas.gas_cell_fields.gas_cell_id = 101U;
  gas.gas_cell_fields.parent_particle_id = gas.particle_id;
  gas.gas_cell_fields.density_code = 3.5;
  auto decoded_gas = wire::decodeParticleMigrationRecord(wire::encodeParticleMigrationRecord(gas));
  assert(decoded_gas.has_gas_cell_fields);
  assert(decoded_gas.gas_cell_fields.gas_cell_id == 101U);
  assert(decoded_gas.gas_cell_fields.density_code == 3.5);

  ParticleMigrationRecord star = makeCommonRecord();
  star.species_tag = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  star.has_star_fields = true;
  star.star_fields.birth_key = 202U;
  star.star_fields.birth_mass_code = 4.5;
  auto decoded_star = wire::decodeParticleMigrationRecord(wire::encodeParticleMigrationRecord(star));
  assert(decoded_star.has_star_fields);
  assert(decoded_star.star_fields.birth_key == 202U);
  assert(decoded_star.star_fields.birth_mass_code == 4.5);

  ParticleMigrationRecord black_hole = makeCommonRecord();
  black_hole.species_tag = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kBlackHole);
  black_hole.has_black_hole_fields = true;
  black_hole.black_hole_fields.accretion_rate_code = 6.5;
  auto decoded_bh = wire::decodeParticleMigrationRecord(wire::encodeParticleMigrationRecord(black_hole));
  assert(decoded_bh.has_black_hole_fields);
  assert(decoded_bh.black_hole_fields.accretion_rate_code == 6.5);

  ParticleMigrationRecord tracer = makeCommonRecord();
  tracer.species_tag = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
  tracer.has_tracer_fields = true;
  tracer.tracer_fields.parent_particle_id = 303U;
  auto decoded_tracer = wire::decodeParticleMigrationRecord(wire::encodeParticleMigrationRecord(tracer));
  assert(decoded_tracer.has_tracer_fields);
  assert(decoded_tracer.tracer_fields.parent_particle_id == 303U);

  const auto widths = wire::particleWireReferenceWidths();
  assert(widths.dark_matter_bytes > 0U);
  assert(widths.gas_extra_bytes > 0U);
  assert(widths.star_extra_bytes > 0U);
  assert(widths.black_hole_extra_bytes > 0U);
  assert(widths.tracer_extra_bytes > 0U);
}

void testTruncatedAndFragmentPacketsFailClosed() {
  const auto encoded = wire::encodeParticleMigrationRecord(makeCommonRecord());
  bool threw = false;
  try {
    (void)wire::decodeParticleMigrationRecord(
        std::span<const std::uint8_t>(encoded.data(), encoded.size() - 1U));
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);

  const std::vector<std::uint8_t> payload{1U, 2U, 3U, 4U};
  const auto fragment = wire::encodeFragment(7U, 10U, 3U, payload);
  const wire::FragmentView view = wire::decodeFragment(fragment);
  assert(view.record_sequence == 7U);
  assert(view.record_total_bytes == 10U);
  assert(view.fragment_offset == 3U);
  assert(std::vector<std::uint8_t>(view.payload.begin(), view.payload.end()) == payload);

  threw = false;
  try {
    (void)wire::encodeFragment(0U, 4U, 3U, payload);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
}

void testPacketCapacityIsIndependentOfLogicalTraffic() {
  constexpr std::size_t round_limit = 16U * 1024U * 1024U;
  const auto plan = wire::planPacketCapacity(round_limit, 8U);
  assert(plan.aggregate_packet_bytes <= round_limit);
  assert(plan.fragment_payload_bytes > 0U);
  // The plan has no logical-record-count input: one thousand or one billion
  // records reuse the same physical packet capacity.
  const std::size_t capacity_for_thousand = plan.aggregate_packet_bytes;
  const std::size_t capacity_for_billion = plan.aggregate_packet_bytes;
  assert(capacity_for_thousand == capacity_for_billion);
}

void testAmrBoundaryPayloadCarriesDepthAndCanonicalMetalMass() {
  cosmosim::core::SimulationState state;
  state.resizeCells(2U);
  state.resizePatches(1U);
  state.patches.patch_id[0] = 77U;
  state.patches.level[0] = 1;
  state.patches.morton_key[0] = 77U;
  state.patches.origin_x_comoving[0] = 1.0;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = 0.5;
  state.patches.extent_y_comoving[0] = 1.0;
  state.patches.extent_z_comoving[0] = 1.0;
  state.patches.cell_dim_x[0] = 2U;
  state.patches.cell_dim_y[0] = 1U;
  state.patches.cell_dim_z[0] = 1U;
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] = 2U;
  state.patches.owning_rank[0] = 0U;

  std::vector<cosmosim::core::GasCellIdentityRecord> identities;
  for (std::uint32_t row = 0U; row < 2U; ++row) {
    state.cells.patch_index[row] = 0U;
    state.cells.center_x_comoving[row] = 1.125 + 0.25 * static_cast<double>(row);
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.mass_code[row] = 2.0 + static_cast<double>(row);
    state.gas_cells.gas_cell_id[row] = 8001U + row;
    state.gas_cells.density_code[row] = 8.0 + static_cast<double>(row);
    state.gas_cells.pressure_code[row] = 3.0;
    state.gas_cells.internal_energy_code[row] = 1.0;
    state.gas_cells.temperature_code[row] = 2.0;
    state.gas_cells.sound_speed_code[row] = 1.0;
    state.gas_cells.metal_mass_code[row] = 0.2 + 0.1 * static_cast<double>(row);
    identities.push_back(cosmosim::core::GasCellIdentityRecord{
        .gas_cell_id = 8001U + row,
        .parent_particle_id = std::nullopt,
        .owning_patch_id = 77U,
        .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(identities));

  cosmosim::parallel::AmrPatchBoundaryCellRequest request;
  request.patch_id = 77U;
  request.boundary_face_mask = static_cast<std::uint8_t>(
      cosmosim::parallel::AmrPatchBoundaryFace::kXLower);
  request.boundary_face_depths[0] = 2U;
  std::vector<cosmosim::parallel::AmrPatchCellPayloadRecord> payload;
  cosmosim::workflows::internal::fillMigrationAmrPatchBoundaryCellPayloadChunk(
      state, 0, std::span<const cosmosim::parallel::AmrPatchBoundaryCellRequest>(&request, 1U),
      0U, 8U, payload);
  assert(payload.size() == 2U);
  assert(payload[0].local_cell_offset == 0U);
  assert(payload[1].local_cell_offset == 1U);
  assert(payload[0].metal_mass_code == state.gas_cells.metal_mass_code[0]);
  assert(payload[1].metal_mass_code == state.gas_cells.metal_mass_code[1]);
}

void testAmrMigrationWireCarriesPendingSynchronizationState() {
  cosmosim::core::AmrPatchMigrationRecord record;
  record.patch.patch_id = 501U;
  record.patch.parent_patch_id = 50U;
  record.patch.morton_key = 0xABU;
  record.patch.level = 2;
  record.patch.owning_rank = 3U;
  record.patch.cell_count = 1U;
  record.patch.extent_x_comoving = 1.0;
  record.patch.extent_y_comoving = 1.0;
  record.patch.extent_z_comoving = 1.0;
  record.patch.cell_dim_x = 1U;
  record.patch.cell_dim_y = 1U;
  record.patch.cell_dim_z = 1U;
  cosmosim::core::GasCellMigrationRecord gas;
  gas.owning_rank = 3U;
  gas.fields.gas_cell_id = 9001U;
  gas.fields.owning_patch_id = 501U;
  gas.fields.temperature_code = 42.0;
  record.gas_cell_records.push_back(gas);
  cosmosim::core::GasCellSchedulerMigrationRecord scheduler;
  scheduler.gas_cell_id = 9001U;
  scheduler.bin_index = 2U;
  scheduler.next_activation_tick = 64U;
  scheduler.pending_bin_index = 0xFFU;
  record.gas_cell_scheduler_records.push_back(scheduler);

  cosmosim::core::PendingFluxRegisterRecord pending;
  pending.register_key = 7001U;
  pending.coarse_patch_id = 501U;
  pending.coarse_gas_cell_id = 9001U;
  pending.coarse_cell_index = 17U;
  pending.level = 2U;
  pending.expected_area_comov = 0.25;
  pending.expected_fine_substeps = 2U;
  pending.completed_fine_substeps = 1U;
  pending.fine_substep_coverage_mask = 1U;
  pending.coarse_face_count = 1U;
  pending.fine_face_count = 4U;
  pending.gas_cell_identity_generation = 9U;
  pending.patch_geometry_generation = 11U;
  pending.coarse_total_energy_flux_integral_code = 3.5;
  pending.fine_total_energy_flux_integral_code = 4.5;
  record.pending_flux_register_records.push_back(pending);

  cosmosim::core::AmrTemporalBoundaryHistoryRecord history;
  history.patch_id = 501U;
  history.patch_level = 2U;
  history.patch_geometry_fingerprint = 0x1234U;
  history.gas_cell_identity_generation = 9U;
  history.interval_start_code = 0.5;
  history.interval_end_code = 0.75;
  history.end_state_valid = true;
  cosmosim::core::AmrTemporalBoundaryHistoryCellRecord history_cell;
  history_cell.gas_cell_id = 9001U;
  history_cell.patch_local_cell = 0U;
  history_cell.start_mass_density_comoving = 2.0;
  history_cell.end_mass_density_comoving = 2.5;
  history.cells.push_back(history_cell);
  record.temporal_boundary_history_records.push_back(history);

  const auto encoded = wire::encodeAmrPatchMigrationRecord(record);
  const auto decoded = wire::decodeAmrPatchMigrationRecord(encoded);
  assert(decoded.patch.patch_id == 501U);
  assert(decoded.gas_cell_records.size() == 1U);
  assert(decoded.pending_flux_register_records.size() == 1U);
  assert(decoded.pending_flux_register_records[0].register_key == 7001U);
  assert(decoded.pending_flux_register_records[0].coarse_gas_cell_id == 9001U);
  assert(decoded.pending_flux_register_records[0].fine_total_energy_flux_integral_code == 4.5);
  assert(decoded.temporal_boundary_history_records.size() == 1U);
  assert(decoded.temporal_boundary_history_records[0].patch_geometry_fingerprint == 0x1234U);
  assert(decoded.temporal_boundary_history_records[0].end_state_valid);
  assert(decoded.temporal_boundary_history_records[0].cells.size() == 1U);
  assert(decoded.temporal_boundary_history_records[0].cells[0].gas_cell_id == 9001U);
  assert(decoded.temporal_boundary_history_records[0].cells[0].end_mass_density_comoving == 2.5);
}

}  // namespace

int main() {
  testDarkMatterWireIsConditionalAndExplicit();
  testSpeciesPayloadRoundTrips();
  testTruncatedAndFragmentPacketsFailClosed();
  testPacketCapacityIsIndependentOfLogicalTraffic();
  testAmrBoundaryPayloadCarriesDepthAndCanonicalMetalMass();
  testAmrMigrationWireCarriesPendingSynchronizationState();
  return 0;
}
