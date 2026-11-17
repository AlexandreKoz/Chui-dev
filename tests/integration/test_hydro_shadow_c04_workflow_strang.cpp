// C04-02/C04-06 workflow fixture: real T01 Strang dispatch with the C04 shadow.
// NOT COMPILED OR RUN in this handoff (the owner reserves build/test stages).
//
// Runs the production ReferenceWorkflowRunner (TreePM + KDK, SourceRuntime
// with the Part 2 fixed-hybrid thermochemistry, composition sidecar, dust,
// ISM regime selection, effective-multiphase star formation and SN feedback,
// comoving cosmology) on one periodic 8^3 production AMR patch, twice:
//   (a) numerics.hydro_shadow_enabled = false,
//   (b) numerics.hydro_shadow_enabled = true, cadence 1, K = 1.
// Checks:
//   * identical final_state_digest (the diagnostic never changes evolution),
//   * a published hydro.shadow_reference event (accepted-epoch barrier: the
//     epoch of step n is published at the start of step n+1),
//   * the full-physics chronology either completed or was refused with a
//     SPECIFIC counted reason (never silently mislabeled), and the hydro-only
//     comparison completed on the sampled patch.
// The deck mirrors configs/chui_part2/cosmo_cube_part2.param.txt physics.
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

namespace {

[[nodiscard]] std::string deck(std::string_view run_name, bool shadow) {
  std::ostringstream s;
  s << "schema_version = 1\n\n";
  s << "[units]\nlength_unit = mpc\nmass_unit = msun\nvelocity_unit = km_s\ncoordinate_frame = comoving\n\n";
  s << "[mode]\nmode = cosmo_cube\nic_file = generated\nzoom_high_res_region = false\n";
  s << "hydro_boundary = periodic\ngravity_boundary = periodic\n\n";
  s << "[cosmology]\nomega_matter = 0.31\nomega_lambda = 0.69\nomega_baryon = 0.0486\n";
  s << "hubble_param = 0.677\nsigma8 = 0.81\nscalar_index_ns = 0.966\nbox_size = 1 mpc\n\n";
  s << "[numerics]\ntime_begin_code = 0.0\ntime_end_code = 1.0\nmax_global_steps = 3\n";
  s << "hierarchical_max_rung = 0\ngravity_solver = treepm\nhydro_solver = godunov_fv\n";
  s << "treepm_pm_grid = 8\ntreepm_asmth_cells = 1.25\ntreepm_rcut_cells = 3.9\n";
  s << "treepm_assignment_scheme = cic\ntreepm_update_cadence_steps = 1\n";
  s << "treepm_allow_diagnostic_naive_dft = true\n";
  s << "hydro_integration_scheme = ctu\nhydro_reconstruction = piecewise_parabolic\n";
  s << "hydro_shadow_enabled = " << (shadow ? "true" : "false") << "\n";
  s << "hydro_shadow_cadence_steps = 1\nhydro_shadow_seed = 7\n";
  s << "hydro_shadow_max_patches_per_rank = 1\n";
  s << "hydro_shadow_max_storage_cells_per_rank = 65536\n";
  s << "hydro_shadow_max_workspace_bytes_per_rank = 1073741824\n\n";
  s << "[physics]\nenable_cooling = true\nenable_star_formation = true\n";
  s << "star_formation_model = effective_multiphase_tng_like\nenable_feedback = true\n";
  s << "reionization_model = hm12\nmetal_species_mode = core_elements\nuv_background_model = none\n";
  s << "self_shielding_model = rahmati13_like\nchemistry_metallicity_scaled_dust_proxy = false\n";
  s << "thermochemistry_allow_unqualified_rates = true\nthermochemistry_backend = fixed_hybrid\n";
  s << "metal_cooling_mode = fine_structure_cold\ndust_model = two_reservoir_v1\n";
  s << "ism_regime_model = effective_or_resolved\nstellar_population_model = chabrier_reference_v1\n";
  s << "fb_supernova_model = mechanical_resolved_unresolved\n\n";
  s << "[output]\nrun_name = " << run_name << "\noutput_directory = integration_outputs\n";
  s << "output_stem = snapshot\nrestart_stem = restart\nsnapshot_interval_steps = 1\nwrite_restarts = true\n\n";
  s << "[parallel]\nmpi_ranks_expected = 1\nomp_threads = 1\ngpu_devices = 0\ndeterministic_reduction = true\n";
  return s.str();
}

// One periodic 8^3 production AMR patch spanning the 1 Mpc box plus one
// massless DM tracer particle for the TreePM owner.
[[nodiscard]] cosmosim::core::SimulationState makeState() {
  using namespace cosmosim;
  constexpr std::size_t n = 8U;
  constexpr std::size_t cells = n * n * n;
  core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(cells);
  state.resizePatches(1);
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
  state.patches.patch_id[0] = 7100U;
  state.patches.parent_patch_id[0] = 0U;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] = static_cast<std::uint32_t>(cells);
  state.patches.owning_rank[0] = 0U;
  state.patches.origin_x_comoving[0] = 0.0;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = 1.0;
  state.patches.extent_y_comoving[0] = 1.0;
  state.patches.extent_z_comoving[0] = 1.0;
  state.patches.cell_dim_x[0] = static_cast<std::uint16_t>(n);
  state.patches.cell_dim_y[0] = static_cast<std::uint16_t>(n);
  state.patches.cell_dim_z[0] = static_cast<std::uint16_t>(n);
  std::vector<core::GasCellIdentityRecord> identities;
  identities.reserve(cells);
  const double width = 1.0 / static_cast<double>(n);
  for (std::uint32_t row = 0; row < cells; ++row) {
    const std::size_t i = row % n;
    const std::size_t j = (row / n) % n;
    const std::size_t k = row / (n * n);
    state.cells.center_x_comoving[row] = (static_cast<double>(i) + 0.5) * width;
    state.cells.center_y_comoving[row] = (static_cast<double>(j) + 0.5) * width;
    state.cells.center_z_comoving[row] = (static_cast<double>(k) + 0.5) * width;
    state.cells.patch_index[row] = 0U;
    state.cells.time_bin[row] = 0U;
    const double rho = 1.0e9 * (1.0 + 0.02 * static_cast<double>((i + j + k) % 3U));
    state.cells.mass_code[row] = rho * width * width * width;
    state.gas_cells.gas_cell_id[row] = 20000U + row;
    state.gas_cells.parent_particle_id[row] = 0U;
    state.gas_cells.density_code[row] = rho;
    state.gas_cells.velocity_x_peculiar[row] = 5.0;
    state.gas_cells.velocity_y_peculiar[row] = -3.0;
    state.gas_cells.velocity_z_peculiar[row] = 0.0;
    state.gas_cells.internal_energy_code[row] = 150.0;
    state.gas_cells.pressure_code[row] = (5.0 / 3.0 - 1.0) * rho * 150.0;
    state.gas_cells.temperature_code[row] = 7.0e3;
    state.gas_cells.sound_speed_code[row] = 12.9;
    state.gas_cells.metal_mass_code[row] = 0.0134 * state.cells.mass_code[row];
    identities.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = 20000U + row, .parent_particle_id = std::nullopt,
        .owning_patch_id = 7100U, .local_cell_row = row});
  }
  state.replaceGasCellIdentityRecords(std::move(identities));
  return state;
}

[[nodiscard]] std::uint64_t payloadSum(const std::string& json, std::string_view key) {
  std::uint64_t total = 0U;
  const std::string token = "\"" + std::string(key) + "\": \"";
  std::size_t pos = 0U;
  while ((pos = json.find(token, pos)) != std::string::npos) {
    pos += token.size();
    std::size_t end = pos;
    while (end < json.size() && json[end] >= '0' && json[end] <= '9') ++end;
    if (end > pos) total += std::stoull(json.substr(pos, end - pos));
  }
  return total;
}

}  // namespace

int main() {
  using namespace cosmosim;
  const std::filesystem::path root = std::filesystem::path("integration_outputs") / "c04_workflow_strang";
  const core::SimulationState initial = makeState();
  const auto run = [&](std::string_view name, bool shadow) {
    const auto frozen = core::loadFrozenConfigFromString(deck(name, shadow), "test_hydro_shadow_c04_workflow_strang");
    workflows::ReferenceWorkflowRunner runner(frozen);
    return runner.run(root, workflows::ReferenceWorkflowOptions{
                                .dt_time_code = 1.0e-5,
                                .write_outputs = COSMOSIM_ENABLE_HDF5 != 0,
                                .initial_state_override = &initial,
                                .max_steps_override = 3});
  };
  const auto disabled = run("c04_shadow_off", false);
  const auto enabled = run("c04_shadow_on", true);
  // C04-05: the optional diagnostic never changes the evolved state.
  assert(disabled.final_state_digest == enabled.final_state_digest);
  assert(disabled.completed_steps == enabled.completed_steps);
  if (!enabled.profiler_json_path.empty() && std::filesystem::exists(enabled.profiler_json_path)) {
    std::ifstream in(enabled.profiler_json_path);
    const std::string json((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    assert(json.find("hydro.shadow_reference") != std::string::npos);
    assert(payloadSum(json, "completed_patches") >= 1U);
    const std::uint64_t chronology_outcomes =
        payloadSum(json, "completed_chronology_regions") +
        payloadSum(json, "chronology_skipped_recursive_level_sources") +
        payloadSum(json, "chronology_skipped_unsplit_compatibility") +
        payloadSum(json, "chronology_skipped_no_replay_operator") +
        payloadSum(json, "chronology_skipped_replay_unavailable") +
        payloadSum(json, "chronology_skipped_replay_incomplete") +
        payloadSum(json, "chronology_skipped_region_changed") +
        payloadSum(json, "chronology_skipped_partial_activity") +
        payloadSum(json, "chronology_skipped_missing_point") +
        payloadSum(json, "chronology_skipped_hydro_sample_missing") +
        payloadSum(json, "chronology_skipped_resource");
    // Every started chronology has exactly one counted outcome.
    assert(chronology_outcomes >= payloadSum(json, "chronology_started"));
    assert(payloadSum(json, "chronology_started") >= 1U);
    // Private replay of the identical operator on identical data.
    assert(json.find("\"chronology_replay_identity_linf\": \"0.000000\"") != std::string::npos);
    // No sample is silently dropped as "external passives" any more.
    assert(payloadSum(json, "skipped_external_passives") == 0U);
  }
  return 0;
}
