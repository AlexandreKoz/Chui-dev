#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/workflows/gravity_source_ownership.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"

#include "../support/test_temp_workspace.hpp"

namespace {

[[nodiscard]] std::string buildConfig(std::string_view run_name) {
  std::ostringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[units]\n";
  stream << "length_unit = kpc\n";
  stream << "mass_unit = msun\n";
  stream << "velocity_unit = km_s\n";
  stream << "coordinate_frame = physical\n\n";
  stream << "[mode]\n";
  stream << "mode = isolated_galaxy\n";
  stream << "ic_file = generated\n";
  stream << "zoom_high_res_region = false\n";
  stream << "hydro_boundary = open\n";
  stream << "gravity_boundary = isolated_monopole\n\n";
  stream << "[numerics]\n";
  stream << "time_begin_code = 0.0\n";
  stream << "time_end_code = 0.000000001\n";
  stream << "max_global_steps = 1\n";
  stream << "hierarchical_max_rung = 0\n";
  stream << "treepm_pm_grid = 4\n";
  stream << "treepm_enable_window_deconvolution = false\n";
  stream << "treepm_update_cadence_steps = 1\n";
  stream << "treepm_allow_diagnostic_naive_dft = true\n\n";
  stream << "[physics]\n";
  stream << "enable_star_formation = true\n";
  stream << "star_formation_model = adaptive_bound_jeans\n";
  stream << "sf_epsilon_ff = 1.0\n";
  stream << "sf_bound_alpha_vir_max = 1.0\n";
  stream << "sf_require_converging_flow = true\n";
  stream << "sf_collapse_timescale = minimum_free_fall_or_compression\n";
  stream << "sf_jeans_mass_floor_code = 100.0\n";
  stream << "sf_star_particle_mass_policy = fixed\n";
  stream << "sf_target_star_particle_mass_code = 100.0\n";
  stream << "sf_min_star_particle_mass_code = 1.0\n";
  stream << "sf_max_star_particle_mass_code = 1000000.0\n";
  stream << "sf_max_spawn_particles_per_cell_step = 8\n";
  stream << "sf_max_fractional_mass_conversion = 0.25\n";
  stream << "sf_min_remaining_gas_fraction = 0.01\n";
  stream << "sf_min_remaining_gas_mass_code = 1.0\n";
  stream << "sf_temperature_safety_ceiling_k = 0.0\n";
  stream << "sf_stochastic_spawning = false\n";
  stream << "sf_random_seed = 987654321\n";
  stream << "sf_effective_t_star_at_threshold = 1.0e-9 code\n\n";
  stream << "[output]\n";
  stream << "run_name = " << run_name << "\n";
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n";
  stream << "snapshot_interval_steps = 1\n";
  stream << "write_restarts = true\n\n";
  stream << "[parallel]\n";
  stream << "mpi_ranks_expected = 1\n";
  stream << "omp_threads = 1\n";
  stream << "gpu_devices = 0\n";
  stream << "deterministic_reduction = true\n";
  return stream.str();
}

[[nodiscard]] cosmosim::core::SimulationState makeState(bool converging) {
  cosmosim::core::SimulationState state;
  state.resizeParticles(3);
  state.resizeCells(3);
  state.resizePatches(1);

  state.patches.patch_id[0] = 7001;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 3;
  state.patches.owning_rank[0] = 0;
  state.patches.origin_x_comoving[0] = 0.0;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = 0.003;
  state.patches.extent_y_comoving[0] = 0.001;
  state.patches.extent_z_comoving[0] = 0.001;
  state.patches.cell_dim_x[0] = 3;
  state.patches.cell_dim_y[0] = 1;
  state.patches.cell_dim_z[0] = 1;

  const std::array<double, 3> velocity_x = converging
      ? std::array<double, 3>{1.0, 0.0, -1.0}
      : std::array<double, 3>{-1.0, 0.0, 1.0};
  std::vector<cosmosim::core::GasCellIdentityRecord> identity_records;
  identity_records.reserve(3);
  for (std::size_t row = 0; row < 3; ++row) {
    const std::uint64_t particle_id = 5001 + row;
    state.particles.position_x_comoving[row] = 0.0005 + 0.001 * row;
    state.particles.position_y_comoving[row] = 0.0005;
    state.particles.position_z_comoving[row] = 0.0005;
    state.particles.velocity_x_peculiar[row] = velocity_x[row];
    state.particles.mass_code[row] = 1.0e6;
    state.particle_sidecar.particle_id[row] = particle_id;
    state.particle_sidecar.sfc_key[row] = particle_id;
    state.particle_sidecar.species_tag[row] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[row] = 0;

    state.cells.center_x_comoving[row] = state.particles.position_x_comoving[row];
    state.cells.center_y_comoving[row] = state.particles.position_y_comoving[row];
    state.cells.center_z_comoving[row] = state.particles.position_z_comoving[row];
    state.cells.mass_code[row] = 1.0e6;
    state.cells.patch_index[row] = 0;
    state.cells.time_bin[row] = 0;
    state.gas_cells.gas_cell_id[row] = 8001 + row;
    state.gas_cells.parent_particle_id[row] = particle_id;
    state.gas_cells.density_code[row] = 1.0e15;
    state.gas_cells.pressure_code[row] = 1.0e3;
    state.gas_cells.internal_energy_code[row] = 1.0e-3;
    state.gas_cells.temperature_code[row] = 100.0;
    state.gas_cells.sound_speed_code[row] = 1.0e-3;
    state.gas_cells.velocity_x_peculiar[row] = velocity_x[row];
    state.gas_cells.metal_mass_code[row] = 2.0e4;
    identity_records.push_back({
        .gas_cell_id = 8001 + row,
        .parent_particle_id = particle_id,
        .owning_patch_id = 7001,
        .local_cell_row = static_cast<std::uint32_t>(row),
    });
  }
  state.species.count_by_species.fill(0);
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 3;
  state.rebuildSpeciesIndex();
  state.replaceGasCellIdentityRecords(std::move(identity_records));
  return state;
}

[[nodiscard]] cosmosim::core::SimulationState makeBatchScalingState(
    std::size_t cell_count) {
  assert(cell_count > 0U);
  assert(cell_count <= std::numeric_limits<std::uint32_t>::max());
  cosmosim::core::SimulationState state;
  state.resizeParticles(cell_count);
  state.resizeCells(cell_count);
  state.resizePatches(1U);

  constexpr double k_cell_width = 0.001;
  state.patches.patch_id[0] = 7901U;
  state.patches.level[0] = 0U;
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] =
      static_cast<std::uint32_t>(cell_count);
  state.patches.owning_rank[0] = 0U;
  state.patches.origin_x_comoving[0] = 0.0;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] =
      k_cell_width * static_cast<double>(cell_count);
  state.patches.extent_y_comoving[0] = k_cell_width;
  state.patches.extent_z_comoving[0] = k_cell_width;
  state.patches.cell_dim_x[0] =
      static_cast<std::uint32_t>(cell_count);
  state.patches.cell_dim_y[0] = 1U;
  state.patches.cell_dim_z[0] = 1U;

  std::vector<cosmosim::core::GasCellIdentityRecord> identity_records;
  identity_records.reserve(cell_count);
  for (std::size_t row = 0U; row < cell_count; ++row) {
    const std::uint64_t particle_id = 900000U + row;
    const std::uint64_t gas_cell_id = 1900000U + row;
    const double x = (static_cast<double>(row) + 0.5) * k_cell_width;
    state.particles.position_x_comoving[row] = x;
    state.particles.position_y_comoving[row] = 0.5 * k_cell_width;
    state.particles.position_z_comoving[row] = 0.5 * k_cell_width;
    state.particles.mass_code[row] = 1.0e6;
    state.particle_sidecar.particle_id[row] = particle_id;
    state.particle_sidecar.sfc_key[row] = particle_id;
    state.particle_sidecar.species_tag[row] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[row] = 0U;

    state.cells.center_x_comoving[row] = x;
    state.cells.center_y_comoving[row] = 0.5 * k_cell_width;
    state.cells.center_z_comoving[row] = 0.5 * k_cell_width;
    state.cells.mass_code[row] = 1.0e6;
    state.cells.patch_index[row] = 0U;
    state.cells.time_bin[row] = 0U;
    state.gas_cells.gas_cell_id[row] = gas_cell_id;
    state.gas_cells.parent_particle_id[row] = particle_id;
    state.gas_cells.density_code[row] = 1.0e15;
    state.gas_cells.pressure_code[row] = 1.0e3;
    state.gas_cells.internal_energy_code[row] = 1.0e-3;
    state.gas_cells.temperature_code[row] = 100.0;
    state.gas_cells.sound_speed_code[row] = 1.0e-3;
    // Zero velocity gives zero divergence, so this fixture exercises the
    // all-active scan/batching path without producing a large star population.
    state.gas_cells.velocity_x_peculiar[row] = 0.0;
    state.gas_cells.velocity_y_peculiar[row] = 0.0;
    state.gas_cells.velocity_z_peculiar[row] = 0.0;
    state.gas_cells.metal_mass_code[row] = 2.0e4;
    identity_records.push_back({
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = particle_id,
        .owning_patch_id = 7901U,
        .local_cell_row = static_cast<std::uint32_t>(row),
    });
  }
  state.species.count_by_species.fill(0U);
  state.species.count_by_species[
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] =
      static_cast<std::uint64_t>(cell_count);
  state.rebuildSpeciesIndex();
  state.replaceGasCellIdentityRecords(std::move(identity_records));
  assert(state.validateOwnershipInvariants());
  return state;
}

[[nodiscard]] std::string readTextFile(const std::filesystem::path& path) {
  std::ifstream input(path);
  assert(input.good());
  return std::string(
      (std::istreambuf_iterator<char>(input)),
      std::istreambuf_iterator<char>());
}

[[nodiscard]] std::uint64_t memoryEntryOwnedCapacity(
    std::string_view profiler_json,
    std::string_view label) {
  const std::string label_token = "\"label\": \"" + std::string(label) + "\"";
  const std::size_t label_pos = profiler_json.find(label_token);
  assert(label_pos != std::string_view::npos);
  constexpr std::string_view k_capacity_key = "\"owned_capacity_bytes\":";
  const std::size_t capacity_key_pos = profiler_json.find(k_capacity_key, label_pos);
  assert(capacity_key_pos != std::string_view::npos);
  std::size_t value_begin = capacity_key_pos + k_capacity_key.size();
  while (value_begin < profiler_json.size() && profiler_json[value_begin] == ' ') {
    ++value_begin;
  }
  std::size_t value_end = value_begin;
  while (value_end < profiler_json.size() &&
         profiler_json[value_end] >= '0' && profiler_json[value_end] <= '9') {
    ++value_end;
  }
  assert(value_end > value_begin);
  return static_cast<std::uint64_t>(std::stoull(
      std::string(profiler_json.substr(value_begin, value_end - value_begin))));
}

[[nodiscard]] std::uint64_t memoryEntryHighWater(
    std::string_view profiler_json,
    std::string_view label) {
  const std::string label_token = "\"label\": \"" + std::string(label) + "\"";
  const std::size_t label_pos = profiler_json.find(label_token);
  assert(label_pos != std::string_view::npos);
  constexpr std::string_view k_high_water_key = "\"high_water_bytes\":";
  const std::size_t key_pos = profiler_json.find(k_high_water_key, label_pos);
  assert(key_pos != std::string_view::npos);
  std::size_t value_begin = key_pos + k_high_water_key.size();
  while (value_begin < profiler_json.size() && profiler_json[value_begin] == ' ') {
    ++value_begin;
  }
  std::size_t value_end = value_begin;
  while (value_end < profiler_json.size() &&
         profiler_json[value_end] >= '0' && profiler_json[value_end] <= '9') {
    ++value_end;
  }
  assert(value_end > value_begin);
  return static_cast<std::uint64_t>(std::stoull(
      std::string(profiler_json.substr(value_begin, value_end - value_begin))));
}

[[nodiscard]] bool memoryEntryGoverned(
    std::string_view profiler_json,
    std::string_view label) {
  const std::string label_token = "\"label\": \"" + std::string(label) + "\"";
  const std::size_t label_pos = profiler_json.find(label_token);
  assert(label_pos != std::string_view::npos);
  constexpr std::string_view k_governed_key = "\"governed_commitment\":";
  const std::size_t key_pos = profiler_json.find(k_governed_key, label_pos);
  assert(key_pos != std::string_view::npos);
  std::size_t value_begin = key_pos + k_governed_key.size();
  while (value_begin < profiler_json.size() && profiler_json[value_begin] == ' ') {
    ++value_begin;
  }
  return profiler_json.substr(value_begin, 4U) == "true";
}

[[nodiscard]] cosmosim::core::SimulationState makeParentlessGasState() {
  auto state = makeState(true);
  std::vector<cosmosim::core::GasCellIdentityRecord> identities(
      state.gas_cell_identity.records().begin(),
      state.gas_cell_identity.records().end());
  identities[1].parent_particle_id.reset();
  state.replaceGasCellIdentityRecords(std::move(identities));
  assert(state.validateOwnershipInvariants());
  return state;
}

[[nodiscard]] cosmosim::core::SimulationState makeSharedLineageGasState() {
  auto state = makeState(true);
  std::vector<cosmosim::core::GasCellIdentityRecord> identities(
      state.gas_cell_identity.records().begin(),
      state.gas_cell_identity.records().end());
  identities[1].parent_particle_id = identities[0].parent_particle_id;
  state.replaceGasCellIdentityRecords(std::move(identities));
  assert(state.validateOwnershipInvariants());
  return state;
}

[[nodiscard]] cosmosim::core::SimulationState withoutCompatibilityMirrors(
    cosmosim::core::SimulationState state) {
  state.resizeParticles(0);
  state.species.count_by_species.fill(0);
  state.rebuildSpeciesIndex();
  std::vector<cosmosim::core::GasCellIdentityRecord> identities(
      state.gas_cell_identity.records().begin(),
      state.gas_cell_identity.records().end());
  for (auto& identity : identities) {
    identity.parent_particle_id.reset();
  }
  state.replaceGasCellIdentityRecords(std::move(identities));
  assert(state.validateOwnershipInvariants());
  return state;
}

[[nodiscard]] cosmosim::core::SimulationState makeRefinementGravityState(bool refined) {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(refined ? 3 : 1);
  state.resizePatches(refined ? 2 : 1);

  state.particles.position_x_comoving[0] = 0.004;
  state.particles.position_y_comoving[0] = 0.0005;
  state.particles.position_z_comoving[0] = 0.0005;
  state.particles.mass_code[0] = 0.0;
  state.particle_sidecar.particle_id[0] = 9701;
  state.particle_sidecar.sfc_key[0] = 9701;
  state.particle_sidecar.species_tag[0] =
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  state.particle_sidecar.owning_rank[0] = 0;
  state.species.count_by_species.fill(0);
  state.species.count_by_species[
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = 1;
  state.rebuildSpeciesIndex();

  state.patches.patch_id[0] = 9801;
  state.patches.parent_patch_id[0] = 0;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 1;
  state.patches.owning_rank[0] = 0;
  state.patches.origin_x_comoving[0] = 0.0005;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = 0.001;
  state.patches.extent_y_comoving[0] = 0.001;
  state.patches.extent_z_comoving[0] = 0.001;
  state.patches.cell_dim_x[0] = 1;
  state.patches.cell_dim_y[0] = 1;
  state.patches.cell_dim_z[0] = 1;

  state.cells.center_x_comoving[0] = 0.001;
  state.cells.center_y_comoving[0] = 0.0005;
  state.cells.center_z_comoving[0] = 0.0005;
  state.cells.mass_code[0] = 2.0e6;
  state.cells.patch_index[0] = 0;
  state.gas_cells.gas_cell_id[0] = 9811;
  state.gas_cells.parent_particle_id[0] = 0;

  std::vector<cosmosim::core::GasCellIdentityRecord> identities;
  identities.push_back({
      .gas_cell_id = 9811,
      .parent_particle_id = std::nullopt,
      .owning_patch_id = 9801,
      .local_cell_row = 0,
  });

  if (refined) {
    state.patches.patch_id[1] = 9802;
    state.patches.parent_patch_id[1] = 9801;
    state.patches.level[1] = 1;
    state.patches.first_cell[1] = 1;
    state.patches.cell_count[1] = 2;
    state.patches.owning_rank[1] = 0;
    state.patches.origin_x_comoving[1] = 0.0005;
    state.patches.origin_y_comoving[1] = 0.0;
    state.patches.origin_z_comoving[1] = 0.0;
    state.patches.extent_x_comoving[1] = 0.001;
    state.patches.extent_y_comoving[1] = 0.001;
    state.patches.extent_z_comoving[1] = 0.001;
    state.patches.cell_dim_x[1] = 2;
    state.patches.cell_dim_y[1] = 1;
    state.patches.cell_dim_z[1] = 1;

    for (std::uint32_t row = 1; row < 3; ++row) {
      const double x = row == 1 ? 0.00075 : 0.00125;
      state.cells.center_x_comoving[row] = x;
      state.cells.center_y_comoving[row] = 0.0005;
      state.cells.center_z_comoving[row] = 0.0005;
      state.cells.mass_code[row] = 1.0e6;
      state.cells.patch_index[row] = 1;
      state.gas_cells.gas_cell_id[row] = 9811 + row;
      state.gas_cells.parent_particle_id[row] = 0;
      identities.push_back({
          .gas_cell_id = 9811 + row,
          .parent_particle_id = std::nullopt,
          .owning_patch_id = 9802,
          .local_cell_row = row,
      });
    }
  }

  state.replaceGasCellIdentityRecords(std::move(identities));
  assert(state.validateOwnershipInvariants());
  return state;
}

struct AuthoritativeGravityNumericalResult {
  std::vector<double> x, y, z, mass;
  std::vector<double> accel_x, accel_y, accel_z;
  std::size_t particle_source_count = 0;
  std::size_t gas_source_count = 0;
};

[[nodiscard]] AuthoritativeGravityNumericalResult solveAuthoritativeGravityNumerically(
    const cosmosim::core::SimulationState& state,
    double epsilon_comoving = 1.0e-5) {
  const auto rows = cosmosim::workflows::internal::selectAuthoritativeGravitySourceRows(
      state, 0U, 1U, "authoritative gas numerical acceptance");
  AuthoritativeGravityNumericalResult result;
  result.particle_source_count = rows.particle_rows.size();
  result.gas_source_count = rows.gas_cell_rows.size();
  const std::size_t source_count = result.particle_source_count + result.gas_source_count;
  result.x.reserve(source_count);
  result.y.reserve(source_count);
  result.z.reserve(source_count);
  result.mass.reserve(source_count);
  for (const std::uint32_t row : rows.particle_rows) {
    result.x.push_back(state.particles.position_x_comoving[row]);
    result.y.push_back(state.particles.position_y_comoving[row]);
    result.z.push_back(state.particles.position_z_comoving[row]);
    result.mass.push_back(state.particles.mass_code[row]);
  }
  for (const std::uint32_t row : rows.gas_cell_rows) {
    const auto source = cosmosim::workflows::internal::authoritativeGasGravitySource(
        state, row, "authoritative gas numerical acceptance");
    result.x.push_back(source.x_comoving);
    result.y.push_back(source.y_comoving);
    result.z.push_back(source.z_comoving);
    result.mass.push_back(source.mass_code);
  }

  std::vector<cosmosim::gravity::TreeLocalIndex> active(source_count);
  for (std::size_t i = 0; i < source_count; ++i) {
    active[i] = cosmosim::gravity::checkedTreeLocalIndex(
        i, "authoritative gas numerical acceptance active source");
  }
  result.accel_x.assign(source_count, 0.0);
  result.accel_y.assign(source_count, 0.0);
  result.accel_z.assign(source_count, 0.0);
  cosmosim::gravity::TreeGravityOptions options;
  options.opening_theta = 1.0e-8;
  options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric;
  options.multipole_order = cosmosim::gravity::TreeMultipoleOrder::kMonopole;
  options.max_leaf_size = 1;
  options.gravitational_constant_code = 1.0;
  options.softening.epsilon_comoving = epsilon_comoving;
  cosmosim::gravity::TreeGravitySolver solver;
  const cosmosim::gravity::GravitySourceGeneration generation{1};
  solver.build(result.x, result.y, result.z, result.mass, options, nullptr, {}, generation);
  solver.evaluateActiveSet(
      result.x, result.y, result.z, result.mass, active,
      result.accel_x, result.accel_y, result.accel_z,
      options, nullptr, {}, {}, generation);

  for (std::size_t target = 0; target < source_count; ++target) {
    std::array<double, 3> direct{};
    for (std::size_t source = 0; source < source_count; ++source) {
      if (source == target) continue;
      const double dx = result.x[source] - result.x[target];
      const double dy = result.y[source] - result.y[target];
      const double dz = result.z[source] - result.z[target];
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double inv_r3 = cosmosim::gravity::softenedInvR3(r2, epsilon_comoving);
      direct[0] += result.mass[source] * dx * inv_r3;
      direct[1] += result.mass[source] * dy * inv_r3;
      direct[2] += result.mass[source] * dz * inv_r3;
    }
    const auto close = [](double actual, double expected) {
      const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
      return std::abs(actual - expected) <= 5.0e-12 * scale;
    };
    assert(close(result.accel_x[target], direct[0]));
    assert(close(result.accel_y[target], direct[1]));
    assert(close(result.accel_z[target], direct[2]));
  }
  return result;
}

[[nodiscard]] cosmosim::core::SimulationState makeCoveredCoarseHierarchyState(bool reverse_rows) {
  cosmosim::core::SimulationState state;
  state.resizeParticles(2);
  state.resizeCells(2);
  state.resizePatches(2);

  state.patches.patch_id[0] = 9101;
  state.patches.parent_patch_id[0] = 0;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = reverse_rows ? 1 : 0;
  state.patches.cell_count[0] = 1;
  state.patches.owning_rank[0] = 0;
  state.patches.origin_x_comoving[0] = 0.0;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = 0.001;
  state.patches.extent_y_comoving[0] = 0.001;
  state.patches.extent_z_comoving[0] = 0.001;
  state.patches.cell_dim_x[0] = 1;
  state.patches.cell_dim_y[0] = 1;
  state.patches.cell_dim_z[0] = 1;

  state.patches.patch_id[1] = 9102;
  state.patches.parent_patch_id[1] = 9101;
  state.patches.level[1] = 1;
  state.patches.first_cell[1] = reverse_rows ? 0 : 1;
  state.patches.cell_count[1] = 1;
  state.patches.owning_rank[1] = 0;
  state.patches.origin_x_comoving[1] = 0.0;
  state.patches.origin_y_comoving[1] = 0.0;
  state.patches.origin_z_comoving[1] = 0.0;
  state.patches.extent_x_comoving[1] = 0.0005;
  state.patches.extent_y_comoving[1] = 0.0005;
  state.patches.extent_z_comoving[1] = 0.0005;
  state.patches.cell_dim_x[1] = 1;
  state.patches.cell_dim_y[1] = 1;
  state.patches.cell_dim_z[1] = 1;

  std::vector<cosmosim::core::GasCellIdentityRecord> identities;
  for (std::uint32_t logical = 0; logical < 2; ++logical) {
    const bool child = logical == 1U;
    const std::uint32_t row = reverse_rows ? 1U - logical : logical;
    const std::uint64_t particle_id = 9201U + logical;
    const std::uint64_t gas_cell_id = 9301U + logical;
    state.particles.position_x_comoving[row] = child ? 0.00025 : 0.0005;
    state.particles.position_y_comoving[row] = child ? 0.00025 : 0.0005;
    state.particles.position_z_comoving[row] = child ? 0.00025 : 0.0005;
    state.particles.mass_code[row] = 1.0e6;
    state.particle_sidecar.particle_id[row] = particle_id;
    state.particle_sidecar.sfc_key[row] = particle_id;
    state.particle_sidecar.species_tag[row] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[row] = 0;
    state.cells.center_x_comoving[row] = state.particles.position_x_comoving[row];
    state.cells.center_y_comoving[row] = state.particles.position_y_comoving[row];
    state.cells.center_z_comoving[row] = state.particles.position_z_comoving[row];
    state.cells.mass_code[row] = 1.0e6;
    state.cells.patch_index[row] = child ? 1U : 0U;
    state.gas_cells.gas_cell_id[row] = gas_cell_id;
    state.gas_cells.parent_particle_id[row] = particle_id;
    state.gas_cells.density_code[row] = 1.0e15;
    state.gas_cells.pressure_code[row] = 1.0e3;
    state.gas_cells.internal_energy_code[row] = 1.0e-3;
    state.gas_cells.temperature_code[row] = 100.0;
    state.gas_cells.sound_speed_code[row] = 1.0e-3;
    state.gas_cells.metal_mass_code[row] = 2.0e4;
    identities.push_back({
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = particle_id,
        .owning_patch_id = child ? 9102U : 9101U,
        .local_cell_row = row,
    });
  }
  state.species.count_by_species.fill(0);
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 2;
  state.rebuildSpeciesIndex();
  state.replaceGasCellIdentityRecords(std::move(identities));
  return state;
}

[[nodiscard]] cosmosim::core::SimulationState makeSingleLevelEffectiveState(
    std::uint32_t level,
    double extent_code,
    std::uint64_t gas_cell_id) {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(1);
  state.resizePatches(1);
  state.patches.patch_id[0] = 9400U + level;
  state.patches.level[0] = level;
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] = 1U;
  state.patches.owning_rank[0] = 0U;
  state.patches.origin_x_comoving[0] = 0.0;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = extent_code;
  state.patches.extent_y_comoving[0] = extent_code;
  state.patches.extent_z_comoving[0] = extent_code;
  state.patches.cell_dim_x[0] = 1U;
  state.patches.cell_dim_y[0] = 1U;
  state.patches.cell_dim_z[0] = 1U;
  const std::uint64_t particle_id = 9500U + level;
  state.particles.position_x_comoving[0] = 0.5 * extent_code;
  state.particles.position_y_comoving[0] = 0.5 * extent_code;
  state.particles.position_z_comoving[0] = 0.5 * extent_code;
  state.particles.mass_code[0] = 1.0e6;
  state.particle_sidecar.particle_id[0] = particle_id;
  state.particle_sidecar.sfc_key[0] = particle_id;
  state.particle_sidecar.species_tag[0] =
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  state.particle_sidecar.owning_rank[0] = 0U;
  state.cells.center_x_comoving[0] = state.particles.position_x_comoving[0];
  state.cells.center_y_comoving[0] = state.particles.position_y_comoving[0];
  state.cells.center_z_comoving[0] = state.particles.position_z_comoving[0];
  state.cells.mass_code[0] = 1.0e6;
  state.cells.patch_index[0] = 0U;
  state.gas_cells.gas_cell_id[0] = gas_cell_id;
  state.gas_cells.parent_particle_id[0] = particle_id;
  state.gas_cells.density_code[0] = 1.0e15;
  state.gas_cells.pressure_code[0] = 1.0e3;
  state.gas_cells.internal_energy_code[0] = 1.0e-3;
  state.gas_cells.temperature_code[0] = 100.0;
  state.gas_cells.sound_speed_code[0] = 1.0e-3;
  state.gas_cells.metal_mass_code[0] = 2.0e4;
  state.species.count_by_species.fill(0U);
  state.species.count_by_species[
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 1U;
  state.rebuildSpeciesIndex();
  state.replaceGasCellIdentityRecords({{
      .gas_cell_id = gas_cell_id,
      .parent_particle_id = particle_id,
      .owning_patch_id = state.patches.patch_id[0],
      .local_cell_row = 0U,
  }});
  return state;
}

[[nodiscard]] bool stageIsBracketedBy(
    const cosmosim::workflows::ReferenceWorkflowReport& report,
    std::string_view middle,
    std::string_view bracket) {
  const auto first = std::find(report.stage_sequence.begin(), report.stage_sequence.end(), bracket);
  const auto center = std::find(report.stage_sequence.begin(), report.stage_sequence.end(), middle);
  const auto last = std::find(report.stage_sequence.rbegin(), report.stage_sequence.rend(), bracket);
  return first != report.stage_sequence.end() && center != report.stage_sequence.end() &&
      last != report.stage_sequence.rend() && first < center && center < last.base() - 1;
}

[[nodiscard]] cosmosim::workflows::ReferenceWorkflowReport runEffectiveHierarchyCase(
    const std::filesystem::path& output_root,
    std::string_view run_name,
    const cosmosim::core::SimulationState& state) {
  std::string config_text = buildConfig(run_name);
  const std::string adaptive = "star_formation_model = adaptive_bound_jeans";
  config_text.replace(config_text.find(adaptive), adaptive.size(),
                      "star_formation_model = effective_multiphase_tng_like");
  config_text.replace(config_text.find("sf_max_spawn_particles_per_cell_step = 8"),
                      std::string("sf_max_spawn_particles_per_cell_step = 8").size(),
                      "sf_max_spawn_particles_per_cell_step = 1");
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      config_text, "test_star_formation_amr_source_runtime");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  return runner.run(
      output_root,
      cosmosim::workflows::ReferenceWorkflowOptions{
          .dt_time_code = 1.0e-9,
          .write_outputs = COSMOSIM_ENABLE_HDF5 != 0,
          .initial_state_override = &state,
          .max_steps_override = 1,
      });
}

[[nodiscard]] cosmosim::workflows::ReferenceWorkflowReport runCase(
    const std::filesystem::path& output_root,
    std::string_view run_name,
    const cosmosim::core::SimulationState& state) {
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      buildConfig(run_name), "test_star_formation_source_runtime");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  return runner.run(
      output_root,
      cosmosim::workflows::ReferenceWorkflowOptions{
          .dt_time_code = 1.0e-9,
          .write_outputs = COSMOSIM_ENABLE_HDF5 != 0,
          .initial_state_override = &state,
          .max_steps_override = 1,
      });
}

[[nodiscard]] std::string moduleSidecarText(
    const cosmosim::core::SimulationState& state,
    std::string_view module_name) {
  const auto* block = state.sidecars.find(std::string(module_name));
  if (block == nullptr) return {};
  return std::string(
      reinterpret_cast<const char*>(block->payload.data()), block->payload.size());
}

[[nodiscard]] cosmosim::workflows::ReferenceWorkflowReport runEffectiveCase(
    const std::filesystem::path& output_root,
    std::string_view run_name,
    const cosmosim::core::SimulationState& state) {
  std::string config_text = buildConfig(run_name);
  const std::string adaptive = "star_formation_model = adaptive_bound_jeans";
  config_text.replace(config_text.find(adaptive), adaptive.size(),
                      "star_formation_model = effective_multiphase_tng_like");
  // Effective parameters rely on the typed reference defaults; the explicit
  // selector alone changes the unresolved-ISM closure for this fixture.
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      config_text, "test_star_formation_source_runtime_effective");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  return runner.run(
      output_root,
      cosmosim::workflows::ReferenceWorkflowOptions{
          .dt_time_code = 1.0e-9,
          .write_outputs = COSMOSIM_ENABLE_HDF5 != 0,
          .initial_state_override = &state,
          .max_steps_override = 1,
      });
}

[[nodiscard]] cosmosim::workflows::ReferenceWorkflowReport runDiffusionMemoryCase(
    const std::filesystem::path& output_root,
    std::string_view run_name,
    const cosmosim::core::SimulationState& state) {
  std::string config_text = buildConfig(run_name);
  const std::string sf_enabled = "enable_star_formation = true";
  const std::size_t sf_pos = config_text.find(sf_enabled);
  assert(sf_pos != std::string::npos);
  config_text.replace(sf_pos, sf_enabled.size(), R"(enable_star_formation = false
enable_feedback = false
metal_species_mode = total_only
enable_metal_diffusion = true
metal_diffusion_model = smagorinsky
metal_diffusion_time_integrator = explicit_subcycling
metal_diffusion_coefficient = 0.05
metal_diffusion_cfl = 0.35
metal_diffusion_max_subcycles = 64
metal_diffusion_max_rkl_stages = 48
metal_diffusion_coefficient_floor_code = 0.0
metal_diffusion_coefficient_ceiling_code = 100.0)"
  );
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      config_text, "test_star_formation_source_runtime_diffusion_memory");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  return runner.run(
      output_root,
      cosmosim::workflows::ReferenceWorkflowOptions{
          .dt_time_code = 1.0e-9,
          .write_outputs = false,
          .initial_state_override = &state,
          .max_steps_override = 1,
      });
}

}  // namespace

int main() {
  const auto workspace =
      cosmosim::test_support::TestTempWorkspace::createProcessLocal(
          "cosmosim_star_formation_source_runtime");
  const auto& output_root = workspace.root();

  const auto converging_state = makeState(true);
  const auto expanding_state = makeState(false);

  // Production ownership selection is tested directly: compatibility gas
  // particles never duplicate authoritative gas-cell mass, lineage metadata
  // does not change source membership, and covered coarse AMR cells are not
  // counted alongside their leaf child.
  const auto converging_gravity_sources =
      cosmosim::workflows::internal::selectAuthoritativeGravitySourceRows(
          converging_state, 0U, 1U, "gas gravity ownership acceptance");
  assert(converging_gravity_sources.particle_rows.empty());
  assert(converging_gravity_sources.gas_cell_rows.size() == converging_state.cells.size());

  const auto parentless_state = makeParentlessGasState();
  const auto shared_lineage_state = makeSharedLineageGasState();
  const auto parentless_gravity_sources =
      cosmosim::workflows::internal::selectAuthoritativeGravitySourceRows(
          parentless_state, 0U, 1U, "parentless gas gravity acceptance");
  const auto shared_lineage_gravity_sources =
      cosmosim::workflows::internal::selectAuthoritativeGravitySourceRows(
          shared_lineage_state, 0U, 1U, "shared-lineage gas gravity acceptance");
  assert(parentless_gravity_sources.gas_cell_rows.size() == converging_gravity_sources.gas_cell_rows.size());
  assert(shared_lineage_gravity_sources.gas_cell_rows.size() == converging_gravity_sources.gas_cell_rows.size());

  const auto parentless_force = solveAuthoritativeGravityNumerically(parentless_state);
  const auto shared_lineage_force = solveAuthoritativeGravityNumerically(shared_lineage_state);
  assert(parentless_force.gas_source_count == parentless_state.cells.size());
  assert(shared_lineage_force.gas_source_count == shared_lineage_state.cells.size());
  assert(shared_lineage_force.accel_x.size() >= 2U);
  assert(shared_lineage_force.accel_x[0] != shared_lineage_force.accel_x[1]);

  // Compatibility mirrors must have zero numerical effect on gas gravity.
  const auto without_mirrors_force = solveAuthoritativeGravityNumerically(
      withoutCompatibilityMirrors(parentless_state));
  assert(parentless_force.gas_source_count == without_mirrors_force.gas_source_count);
  assert(parentless_force.mass == without_mirrors_force.mass);
  for (std::size_t i = 0; i < parentless_force.accel_x.size(); ++i) {
    assert(std::abs(parentless_force.accel_x[i] - without_mirrors_force.accel_x[i]) < 1.0e-12);
    assert(std::abs(parentless_force.accel_y[i] - without_mirrors_force.accel_y[i]) < 1.0e-12);
    assert(std::abs(parentless_force.accel_z[i] - without_mirrors_force.accel_z[i]) < 1.0e-12);
  }

  const auto coarse_force = solveAuthoritativeGravityNumerically(makeRefinementGravityState(false));
  const auto refined_force = solveAuthoritativeGravityNumerically(makeRefinementGravityState(true));
  assert(coarse_force.gas_source_count == 1U);
  assert(refined_force.gas_source_count == 2U); // covered coarse parent excluded
  const double coarse_gas_mass = std::accumulate(
      coarse_force.mass.begin() + static_cast<std::ptrdiff_t>(coarse_force.particle_source_count),
      coarse_force.mass.end(), 0.0);
  const double refined_gas_mass = std::accumulate(
      refined_force.mass.begin() + static_cast<std::ptrdiff_t>(refined_force.particle_source_count),
      refined_force.mass.end(), 0.0);
  assert(std::abs(coarse_gas_mass - refined_gas_mass) < 1.0e-12);
  // The massless collisionless probe is first in both source arrays. Genuine
  // refinement changes the discretized force while preserving total mass/COM.
  assert(std::abs(coarse_force.accel_x[0] - refined_force.accel_x[0]) > 0.0);

  const auto covered_hierarchy_state = makeCoveredCoarseHierarchyState(false);
  const auto covered_gravity_sources =
      cosmosim::workflows::internal::selectAuthoritativeGravitySourceRows(
          covered_hierarchy_state, 0U, 1U, "covered coarse gas gravity acceptance");
  assert(covered_gravity_sources.particle_rows.empty());
  assert(covered_gravity_sources.gas_cell_rows.size() == 1U);
  assert(covered_hierarchy_state.patches.level[
      covered_hierarchy_state.cells.patch_index[covered_gravity_sources.gas_cell_rows.front()]] == 1U);
  const auto converging_report = runCase(output_root, "sf_runtime_converging", converging_state);
  const auto expanding_report = runCase(output_root, "sf_runtime_expanding", expanding_state);
  const auto parentless_report = runCase(
      output_root, "sf_runtime_parentless_gas_gravity", parentless_state);
  const auto shared_lineage_report = runCase(
      output_root, "sf_runtime_shared_lineage_gas_gravity", shared_lineage_state);
  const auto effective_report = runEffectiveCase(
      output_root, "sf_runtime_effective", converging_state);
  const auto hierarchy_report = runEffectiveHierarchyCase(
      output_root, "sf_runtime_amr_hierarchy", covered_hierarchy_state);
  const auto reordered_hierarchy_report = runEffectiveHierarchyCase(
      output_root, "sf_runtime_amr_hierarchy_reordered", makeCoveredCoarseHierarchyState(true));
  const auto level0_report = runEffectiveHierarchyCase(
      output_root, "sf_runtime_effective_level0",
      makeSingleLevelEffectiveState(0U, 0.001, 9601U));
  const auto level1_report = runEffectiveHierarchyCase(
      output_root, "sf_runtime_effective_level1",
      makeSingleLevelEffectiveState(1U, 0.0005, 9602U));
  constexpr std::size_t k_star_formation_batch_acceptance_cells = 4097U;
  const auto batch_scaling_report = runCase(
      output_root,
      "sf_runtime_batch_scaling",
      makeBatchScalingState(k_star_formation_batch_acceptance_cells));
  const auto diffusion_memory_report = runDiffusionMemoryCase(
      output_root,
      "source_runtime_diffusion_memory",
      converging_state);

  assert(converging_report.completed_steps == 1);
  assert(expanding_report.completed_steps == 1);
  assert(parentless_report.completed_steps == 1);
  assert(shared_lineage_report.completed_steps == 1);
  assert(effective_report.completed_steps == 1);
  assert(hierarchy_report.completed_steps == 1);
  assert(reordered_hierarchy_report.completed_steps == 1);
  assert(level0_report.completed_steps == 1);
  assert(level1_report.completed_steps == 1);
  assert(batch_scaling_report.completed_steps == 1);
  assert(diffusion_memory_report.completed_steps == 1);
  assert(converging_report.local_particle_count > 3);
  assert(expanding_report.local_particle_count == 3);
  assert(effective_report.local_particle_count > 3);
  // The covered coarse cell is rejected; exactly one authoritative fine child
  // creates one star even when dense rows are reordered.
  assert(hierarchy_report.local_particle_count == 3);
  assert(reordered_hierarchy_report.local_particle_count == 3);
  // The physical threshold is independent of AMR level and both synchronized
  // single-level representations produce exactly one birth.
  assert(level0_report.local_particle_count == 2);
  assert(level1_report.local_particle_count == 2);
  // Smooth sources bracket hydro in the deterministic Strang reference; star
  // formation remains a once-only endpoint transaction in the post half.
  assert(stageIsBracketedBy(hierarchy_report, "hydro_update", "source_terms"));


  // M2C source-phase ownership must be visible in the same persisted memory
  // report used by the governor/provenance path, even when optional modules are
  // disabled and therefore own zero capacity.
  {
    const std::string profiler_json = readTextFile(converging_report.profiler_json_path);
    assert(profiler_json.find("sources.star_formation.input_batch") !=
           std::string::npos);
    assert(profiler_json.find("sources.stellar_feedback.event_batch") !=
           std::string::npos);
    assert(profiler_json.find("sources.metal_diffusion.workspace") !=
           std::string::npos);
  }

  // Cross the 4096-cell production batch boundary and assert that retained
  // star-formation input capacity plateaus at one bounded batch rather than
  // scaling with the all-active gas population.
  {
    const std::string profiler_json =
        readTextFile(batch_scaling_report.profiler_json_path);
    const std::uint64_t input_capacity_bytes = memoryEntryOwnedCapacity(
        profiler_json, "sources.star_formation.input_batch");
    const std::uint64_t maximum_batch_bytes =
        4096ULL * sizeof(cosmosim::physics::StarFormationCellInput);
    assert(input_capacity_bytes <= maximum_batch_bytes);
    assert(input_capacity_bytes <
           static_cast<std::uint64_t>(k_star_formation_batch_acceptance_cells) *
               sizeof(cosmosim::physics::StarFormationCellInput));
  }

  // The scientifically required diffusion face graph is retained phase
  // topology, but its physical capacity must still be governor-admitted and
  // historical high-water must never under-report current retained capacity.
  {
    const std::string profiler_json =
        readTextFile(diffusion_memory_report.profiler_json_path);
    const std::uint64_t faces_capacity = memoryEntryOwnedCapacity(
        profiler_json, "sources.metal_diffusion.faces");
    const std::uint64_t faces_high_water = memoryEntryHighWater(
        profiler_json, "sources.metal_diffusion.faces");
    assert(faces_capacity > 0U);
    assert(faces_high_water >= faces_capacity);
    assert(memoryEntryGoverned(
        profiler_json, "sources.metal_diffusion.faces"));
    const std::uint64_t rho_capacity = memoryEntryOwnedCapacity(
        profiler_json, "sources.metal_diffusion.rho_kappa");
    const std::uint64_t rho_high_water = memoryEntryHighWater(
        profiler_json, "sources.metal_diffusion.rho_kappa");
    assert(rho_high_water >= rho_capacity);
  }

#if COSMOSIM_ENABLE_HDF5
  assert(converging_report.restart_roundtrip_ok);
  assert(effective_report.restart_roundtrip_ok);
  const auto restart = cosmosim::io::readRestartCheckpointHdf5(converging_report.restart_path);
  const std::string sf_sidecar = moduleSidecarText(restart.state, "star_formation");
  assert(sf_sidecar.find("counter_scope=latest_and_cumulative_source_stages") != std::string::npos);
  assert(sf_sidecar.find("cumulative.spawned_particles=") != std::string::npos);
  assert(sf_sidecar.find("cumulative.mass_residual_code=") != std::string::npos);
  const auto effective_restart =
      cosmosim::io::readRestartCheckpointHdf5(effective_report.restart_path);
  const std::string effective_sidecar =
      moduleSidecarText(effective_restart.state, "effective_multiphase_ism");
  assert(effective_sidecar.find("table_hash=") != std::string::npos);
  assert(effective_sidecar.find("cumulative.energy_added_code=") != std::string::npos);
  assert(restart.state.star_particles.size() == converging_report.local_particle_count - 3);
  double gas_particle_mass_code = 0.0;
  double stellar_particle_mass_code = 0.0;
  for (std::size_t particle_row = 0; particle_row < restart.state.particles.size(); ++particle_row) {
    const auto species = static_cast<cosmosim::core::ParticleSpecies>(
        restart.state.particle_sidecar.species_tag[particle_row]);
    if (species == cosmosim::core::ParticleSpecies::kGas) {
      gas_particle_mass_code += restart.state.particles.mass_code[particle_row];
    } else if (species == cosmosim::core::ParticleSpecies::kStar) {
      stellar_particle_mass_code += restart.state.particles.mass_code[particle_row];
    }
  }
  double gas_cell_mass_code = 0.0;
  for (const double mass_code : restart.state.cells.mass_code) {
    gas_cell_mass_code += mass_code;
  }
  assert(std::abs(gas_particle_mass_code - gas_cell_mass_code) < 1.0e-8);
  assert(stellar_particle_mass_code > 0.0);
  assert(std::abs(
      gas_particle_mass_code + stellar_particle_mass_code -
      (gas_cell_mass_code + stellar_particle_mass_code)) < 1.0e-8);
  for (std::size_t star_row = 0; star_row < restart.state.star_particles.size(); ++star_row) {
    assert(restart.state.star_particles.parent_gas_cell_id[star_row] >= 8001);
    assert(restart.state.star_particles.parent_gas_cell_id[star_row] <= 8003);
    assert(restart.state.star_particles.birth_key[star_row] != 0);
    assert(std::abs(restart.state.star_particles.metallicity_mass_fraction[star_row] - 0.02) < 1.0e-10);
  }
#endif

  return 0;
}
