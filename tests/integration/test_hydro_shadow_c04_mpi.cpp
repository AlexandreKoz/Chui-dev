// C04-03/C04-04/C04-05 MPI fixture for the sampled hydro shadow.
// NOT COMPILED OR RUN in this handoff (the owner reserves build/test stages).
//
// Registered for 2 and 4 ranks (and runnable on 1). Four production 8^3 AMR
// root patches tile a 4 x 1 x 1 periodic box; patch p is owned by rank
// p % world_size and its patch ID depends on the world size, so neither IDs
// nor owners are physical. With CTU/PPM, K = 2 and a fixed seed:
//   * every published epoch selects EXACTLY the regions predicted from the
//     physical keys alone (rank-decomposition independent selection, one
//     global K, no rank-local ID-zero duplicates),
//   * MPI interface ghosts receive post-advance endpoint values through the
//     directed exchange (no endpoint unavailability, no cone exclusion),
//   * enabling the shadow, or forcing a uniform resource refusal, leaves the
//     final state digest bit-identical to the shadow-disabled run.
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/hydro/hydro_shadow_reference.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

#include "../support/mpi_test_workspace.hpp"

namespace {

constexpr std::uint64_t k_seed = 424242U;
constexpr std::size_t k_global_k = 2U;
constexpr int k_patches = 4;
constexpr std::size_t k_n = 8U;

[[nodiscard]] std::string configText(int world_size, std::string_view run_name, bool shadow,
                                      std::uint64_t workspace_bytes, std::uint64_t storage_cells = 65536U) {
  std::ostringstream s;
  s << "schema_version = 1\n\n[mode]\nmode = cosmo_cube\nic_file = generated\nhydro_boundary = periodic\n\n";
  s << "[cosmology]\nbox_size_x = 4.0\nbox_size_y = 1.0\nbox_size_z = 1.0\n\n";
  s << "[numerics]\ntime_begin_code = 0.01\ntime_end_code = 0.0103\nmax_global_steps = 3\n";
  s << "hierarchical_max_rung = 0\ntreepm_pm_grid_nx = 32\ntreepm_pm_grid_ny = 8\ntreepm_pm_grid_nz = 8\n";
  s << "treepm_asmth_cells = 1.25\ntreepm_rcut_cells = 3.9\ntreepm_update_cadence_steps = 1\n";
  s << "hydro_integration_scheme = ctu\nhydro_reconstruction = piecewise_parabolic\n";
  s << "hydro_shadow_enabled = " << (shadow ? "true" : "false") << "\n";
  s << "hydro_shadow_cadence_steps = 1\nhydro_shadow_seed = " << k_seed << "\n";
  s << "hydro_shadow_max_patches_per_rank = " << k_global_k << "\n";
  s << "hydro_shadow_max_storage_cells_per_rank = " << storage_cells << "\n";
  s << "hydro_shadow_max_workspace_bytes_per_rank = " << workspace_bytes << "\n\n";
  s << "[physics]\nenable_cooling = false\nenable_star_formation = false\nenable_feedback = false\n";
  s << "enable_stellar_evolution = false\nenable_black_hole_agn = false\nenable_tracers = false\n";
  s << "enable_metal_diffusion = false\n\n";
  s << "[output]\nrun_name = " << run_name << "\noutput_directory = integration_outputs\n";
  s << "output_stem = snapshot\nrestart_stem = restart\nsnapshot_interval_steps = 1\nwrite_restarts = true\n\n";
  s << "[parallel]\nmpi_ranks_expected = " << world_size << "\ndecomposition_runtime_rebalance_enabled = false\n";
  return s.str();
}

[[nodiscard]] cosmosim::amr::PatchDescriptor descriptorFor(int patch, int world_size) {
  return cosmosim::amr::PatchDescriptor{
      .patch_id = 50000U + 97U * static_cast<std::uint64_t>(world_size) + static_cast<std::uint64_t>(patch),
      .parent_patch_id = 0U,
      .level = 0U,
      .morton_key = static_cast<std::uint64_t>(patch),
      .origin_comov = {static_cast<double>(patch), 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {static_cast<std::uint16_t>(k_n), static_cast<std::uint16_t>(k_n),
                    static_cast<std::uint16_t>(k_n)}};
}

[[nodiscard]] cosmosim::core::SimulationState makeState(int world_size) {
  namespace core = cosmosim::core;
  constexpr std::size_t per_patch = k_n * k_n * k_n;
  constexpr std::size_t cells = per_patch * static_cast<std::size_t>(k_patches);
  core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(cells);
  state.resizePatches(static_cast<std::size_t>(k_patches));
  state.metadata.run_name = "hydro_shadow_c04_mpi";
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";
  state.particles.position_x_comoving[0] = 0.5;
  state.particles.position_y_comoving[0] = 0.5;
  state.particles.position_z_comoving[0] = 0.5;
  state.particles.mass_code[0] = 0.0;
  state.particle_sidecar.particle_id[0] = 1U;
  state.particle_sidecar.sfc_key[0] = 1U;
  state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
  state.particle_sidecar.owning_rank[0] = 0;
  state.species.count_by_species.fill(0U);
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)] = 1U;
  state.rebuildSpeciesIndex();
  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(cells);
  const double width = 1.0 / static_cast<double>(k_n);
  std::uint32_t row = 0U;
  for (int patch = 0; patch < k_patches; ++patch) {
    const auto descriptor = descriptorFor(patch, world_size);
    const auto p = static_cast<std::size_t>(patch);
    state.patches.patch_id[p] = descriptor.patch_id;
    state.patches.level[p] = 0;
    state.patches.first_cell[p] = row;
    state.patches.cell_count[p] = static_cast<std::uint32_t>(per_patch);
    state.patches.owning_rank[p] = static_cast<std::uint32_t>(patch % world_size);
    state.patches.origin_x_comoving[p] = descriptor.origin_comov[0];
    state.patches.origin_y_comoving[p] = 0.0;
    state.patches.origin_z_comoving[p] = 0.0;
    state.patches.extent_x_comoving[p] = 1.0;
    state.patches.extent_y_comoving[p] = 1.0;
    state.patches.extent_z_comoving[p] = 1.0;
    state.patches.cell_dim_x[p] = static_cast<std::uint16_t>(k_n);
    state.patches.cell_dim_y[p] = static_cast<std::uint16_t>(k_n);
    state.patches.cell_dim_z[p] = static_cast<std::uint16_t>(k_n);
    for (std::size_t c = 0; c < per_patch; ++c, ++row) {
      const std::size_t i = c % k_n;
      const std::size_t j = (c / k_n) % k_n;
      const std::size_t k = c / (k_n * k_n);
      const double x = static_cast<double>(patch) + (static_cast<double>(i) + 0.5) * width;
      const double rho = 1.0 + 0.05 * static_cast<double>((i + 2U * j + 3U * k + static_cast<std::size_t>(patch)) % 5U);
      const double pressure = 1.0;
      state.cells.center_x_comoving[row] = x;
      state.cells.center_y_comoving[row] = (static_cast<double>(j) + 0.5) * width;
      state.cells.center_z_comoving[row] = (static_cast<double>(k) + 0.5) * width;
      state.cells.mass_code[row] = rho * width * width * width;
      state.cells.patch_index[row] = static_cast<std::uint32_t>(patch);
      state.cells.time_bin[row] = 0U;
      state.gas_cells.gas_cell_id[row] = 800000U + row;
      state.gas_cells.density_code[row] = rho;
      state.gas_cells.pressure_code[row] = pressure;
      state.gas_cells.internal_energy_code[row] = pressure / ((5.0 / 3.0 - 1.0) * rho);
      state.gas_cells.temperature_code[row] = pressure / rho;
      state.gas_cells.sound_speed_code[row] = std::sqrt((5.0 / 3.0) * pressure / rho);
      state.gas_cells.velocity_x_peculiar[row] = 0.05;
      records.push_back(core::GasCellIdentityRecord{
          .gas_cell_id = 800000U + row, .parent_particle_id = std::nullopt,
          .owning_patch_id = descriptor.patch_id, .local_cell_row = row});
    }
  }
  state.replaceGasCellIdentityRecords(std::move(records));
  return state;
}

// Expected selection from physical keys alone (decomposition independent).
[[nodiscard]] std::string expectedKeys(std::uint64_t step) {
  std::vector<cosmosim::hydro::HydroShadowCandidate> candidates;
  for (int patch = 0; patch < k_patches; ++patch) {
    const auto descriptor = descriptorFor(patch, 1);
    candidates.push_back(cosmosim::hydro::HydroShadowCandidate{
        .patch_id = 0U, .level = 0, .region_key = cosmosim::amr::hydroShadowRegionKeyForPatch(descriptor)});
  }
  const auto ranked = cosmosim::hydro::rankHydroShadowCandidates(candidates, k_global_k, step, k_seed);
  std::string keys;
  for (const auto& candidate : ranked) {
    if (!keys.empty()) keys += ',';
    keys += std::to_string(cosmosim::hydro::hydroShadowEffectiveKey(candidate)) + ":0";
  }
  return keys;
}

[[nodiscard]] std::vector<std::string> payloadValues(const std::string& json, std::string_view key) {
  std::vector<std::string> values;
  const std::string token = "\"" + std::string(key) + "\": \"";
  std::size_t pos = 0U;
  while ((pos = json.find(token, pos)) != std::string::npos) {
    pos += token.size();
    const std::size_t end = json.find('"', pos);
    values.push_back(json.substr(pos, end - pos));
  }
  return values;
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI && COSMOSIM_ENABLE_HDF5
  MPI_Init(nullptr, nullptr);
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 1 && world_size != 2 && world_size != 4) {
    MPI_Finalize();
    return 0;
  }
  const cosmosim::core::SimulationState initial = makeState(world_size);
  auto workspace = cosmosim::test_support::createMpiSharedWorkspace("cosmosim_hydro_shadow_c04_mpi");
  const std::filesystem::path& root = workspace.root();
  const auto run = [&](std::string_view name, bool shadow, std::uint64_t workspace_bytes,
                       std::uint64_t storage_cells = 65536U) {
    const auto frozen = cosmosim::core::loadFrozenConfigFromString(
        configText(world_size, name, shadow, workspace_bytes, storage_cells), "test_hydro_shadow_c04_mpi");
    cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
    return runner.run(root / std::string(name),
                      cosmosim::workflows::ReferenceWorkflowOptions{
                          .write_outputs = true, .initial_state_override = &initial, .max_steps_override = 3});
  };
  const auto off = run("shadow_off", false, 1ULL << 30U);
  const auto on = run("shadow_on", true, 1ULL << 30U);
  const auto refused = run("shadow_refused", true, 64U);
  // Rank-ASYMMETRIC local refusal: the per-rank storage cap admits one
  // sampled patch (8^3 + ghost shell) but not two. Ranks owning both selected
  // regions refuse the second locally while other ranks admit theirs; the
  // following collective stages (endpoint vote/exchange, publication) must
  // proceed on every rank, and the physical selection must not change.
  const auto asymmetric = run("shadow_asymmetric", true, 1ULL << 30U, 4000U);
  assert(off.final_state_digest == on.final_state_digest);
  assert(off.final_state_digest == refused.final_state_digest);
  assert(off.final_state_digest == asymmetric.final_state_digest);
  if (world_rank == 0 && std::filesystem::exists(on.profiler_json_path)) {
    std::ifstream in(on.profiler_json_path);
    const std::string json((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const auto keys = payloadValues(json, "selected_region_keys");
    const auto steps = payloadValues(json, "epoch_step_index");
    assert(!keys.empty() && keys.size() == steps.size());
    for (std::size_t e = 0; e < keys.size(); ++e) {
      assert(keys[e] == expectedKeys(std::stoull(steps[e])));
    }
    for (const auto& value : payloadValues(json, "completed_unique_physical_regions")) {
      assert(value == std::to_string(k_global_k));
    }
    for (const auto& value : payloadValues(json, "skipped_endpoint_unavailable")) assert(value == "0");
    for (const auto& value : payloadValues(json, "cone_excluded_cells")) assert(value == "0");
  }
  if (world_rank == 0 && std::filesystem::exists(refused.profiler_json_path)) {
    std::ifstream in(refused.profiler_json_path);
    const std::string json((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    for (const auto& value : payloadValues(json, "completed_patches")) assert(value == "0");
  }
  if (world_rank == 0 && std::filesystem::exists(asymmetric.profiler_json_path)) {
    std::ifstream in(asymmetric.profiler_json_path);
    const std::string json((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const auto keys = payloadValues(json, "selected_region_keys");
    const auto steps = payloadValues(json, "epoch_step_index");
    const auto attempted = payloadValues(json, "attempted_regions");
    const auto completed = payloadValues(json, "completed_patches");
    const auto refused_resource = payloadValues(json, "skipped_resource");
    assert(keys.size() == steps.size() && attempted.size() == completed.size() &&
           completed.size() == refused_resource.size());
    for (std::size_t e = 0; e < keys.size(); ++e) {
      // Identity is decided before admission: refusal never swaps regions.
      assert(keys[e] == expectedKeys(std::stoull(steps[e])));
      assert(std::stoull(completed[e]) + std::stoull(refused_resource[e]) == std::stoull(attempted[e]));
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  return 0;
}
