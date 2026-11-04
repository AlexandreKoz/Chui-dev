// CHUI final-design kernel benchmarks (V5). SOURCE ONLY — NOT RUN —
// USER-RESERVED. No timing, speedup or memory number is claimed anywhere.
//
// Usage: bench_chui_final_kernels [root_patches_per_axis=4] [repeats=5]
// Output: one "key value" line per metric; median of `repeats` after one
// warm-up. Metrics:
//   recursive_vs_synchronized_hydro: seconds per coarse step of a 3-level
//     production hierarchy (all roots refined once, every 4th level-1 patch
//     refined again; 2:1 balanced) for both drivers, and their ratio;
//   amr_regrid_plan_seconds: planProductionAmrRegrid with every R1 indicator;
//   association_fof_seconds: discoverAssociationLeaders on N clustered sources;
//   feedback_stage_seconds: MechanicalFeedbackDeposition::stage per source;
//   rollback_copy_seconds / rollback_copy_bytes: one SimulationState copy (the
//     governed full-state retry point) and its accounted baseline bytes;
//   governor_peak_accounted_bytes: MemoryGovernor accounted high-water (baseline
//     + reserved + committed) during the hydro step (accounting, not RSS).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <random>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/amr/amr_production_regrid.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/physics/supernova_feedback.hpp"

namespace {

namespace amr = cosmosim::amr;
namespace core = cosmosim::core;
namespace hy = cosmosim::hydro;
namespace ph = cosmosim::physics;
using Clock = std::chrono::steady_clock;
constexpr double k_gamma = 5.0 / 3.0;

template <class F>
double medianSeconds(int repeats, F&& f) {
  f();  // warm-up
  std::vector<double> samples;
  for (int r = 0; r < repeats; ++r) {
    const auto t0 = Clock::now();
    f();
    samples.push_back(std::chrono::duration<double>(Clock::now() - t0).count());
  }
  std::sort(samples.begin(), samples.end());
  return samples[samples.size() / 2U];
}

amr::ProductionAmrHydroOptions options(const core::SimulationState& state, core::MemoryGovernor& governor,
                                       double box) {
  governor.setBaselineOwnedBytes(core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(state)));
  amr::ProductionAmrHydroOptions o;
  o.physical_boundary_kind = hy::HydroBoundaryKind::kPeriodic;
  o.periodic_box_extent_comoving = {box, box, box};
  o.adiabatic_index = k_gamma;
  o.density_floor = 1.0e-12;
  o.pressure_floor = 1.0e-12;
  o.refinement_ratio = 2U;
  o.regrid_memory_governor = &governor;
  return o;
}

core::SimulationState makeHierarchy(std::uint32_t roots_per_axis) {
  const std::uint32_t n_roots = roots_per_axis * roots_per_axis * roots_per_axis;
  core::SimulationState state;
  state.resizeCells(8U * n_roots);
  state.resizePatches(n_roots);
  std::vector<core::GasCellIdentityRecord> records;
  for (std::uint32_t p = 0; p < n_roots; ++p) {
    const std::array<double, 3> origin{static_cast<double>(p % roots_per_axis),
                                       static_cast<double>((p / roots_per_axis) % roots_per_axis),
                                       static_cast<double>(p / (roots_per_axis * roots_per_axis))};
    amr::writePatchDescriptorToStateRow(state, p, amr::PatchDescriptor{
        .patch_id = 1U + p, .parent_patch_id = 0U, .level = 0U, .morton_key = 1U + p,
        .origin_comov = origin, .extent_comov = {1.0, 1.0, 1.0}, .cell_dims = {2U, 2U, 2U}});
    state.patches.first_cell[p] = 8U * p;
    state.patches.cell_count[p] = 8U;
    state.patches.owning_rank[p] = 0;
    for (std::uint32_t c = 0; c < 8U; ++c) {
      const std::uint32_t row = 8U * p + c;
      const double x = origin[0] + 0.25 + 0.5 * (c % 2U);
      const double rho = 1.0 + 0.3 * std::sin(2.0 * x);
      state.cells.center_x_comoving[row] = x;
      state.cells.center_y_comoving[row] = origin[1] + 0.25 + 0.5 * ((c / 2U) % 2U);
      state.cells.center_z_comoving[row] = origin[2] + 0.25 + 0.5 * (c / 4U);
      state.cells.patch_index[row] = p;
      state.cells.mass_code[row] = rho * 0.125;
      state.gas_cells.gas_cell_id[row] = 100U + row;
      state.gas_cells.density_code[row] = rho;
      state.gas_cells.pressure_code[row] = 1.0;
      state.gas_cells.internal_energy_code[row] = 1.0 / ((k_gamma - 1.0) * rho);
      state.gas_cells.velocity_x_peculiar[row] = 0.1;
      state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma / rho);
      state.gas_cells.metal_mass_code[row] = 0.01 * state.cells.mass_code[row];
      records.push_back({.gas_cell_id = 100U + row, .parent_particle_id = std::nullopt, .owning_patch_id = 1U + p,
                         .local_cell_row = row});
    }
  }
  state.gas_cell_identity.assign(std::move(records));
  const double box = static_cast<double>(roots_per_axis);
  std::uint64_t next_patch = 10000U, next_cell = 1000000U;
  const auto refineWhere = [&](auto&& select) {
    const auto descriptors = amr::buildProductionAmrPatchDescriptors(state);
    std::size_t k = 0;
    for (const auto& d : descriptors) {
      if (!select(d, k++)) continue;
      core::MemoryGovernor governor;
      (void)amr::refineProductionPatchInSimulationState(state, d, next_patch, next_cell, options(state, governor, box));
      next_patch += 8U;
      next_cell += 64U;
    }
  };
  refineWhere([](const amr::PatchDescriptor& d, std::size_t) { return d.level == 0U; });
  refineWhere([](const amr::PatchDescriptor& d, std::size_t k) { return d.level == 1U && k % 4U == 0U; });
  return state;
}

}  // namespace

int main(int argc, char** argv) {
  const std::uint32_t roots = argc > 1 ? static_cast<std::uint32_t>(std::strtoul(argv[1], nullptr, 10)) : 4U;
  const int repeats = argc > 2 ? std::atoi(argv[2]) : 5;
  const double box = static_cast<double>(roots);
  const core::SimulationState hierarchy = makeHierarchy(roots);
  std::printf("cells %zu\npatches %zu\n", hierarchy.cells.size(), hierarchy.patches.size());

  const hy::HydroCoreSolver solver(k_gamma);
  hy::HllcRiemannSolver riemann;
  const hy::HydroUpdateContext update{.dt_code = 1.0e-3, .time_begin_code = 0.0, .time_end_code = 1.0e-3,
                                      .scale_factor = 1.0, .hubble_rate_code = 0.0};
  const hy::HydroSourceContext source_context{.update = update};
  std::uint64_t peak = 0U;
  const auto hydro = [&](bool recursive) {
    return medianSeconds(repeats, [&] {
      core::SimulationState state = hierarchy;
      core::MemoryGovernor governor;
      const auto o = options(state, governor, box);
      (void)(recursive ? amr::advanceProductionAmrHydroRecursive(state, {}, update, source_context, solver, riemann, {}, o)
                       : amr::advanceProductionAmrHydro(state, {}, update, source_context, solver, riemann, {}, o));
      peak = std::max(peak, governor.snapshot().peak_accounted_bytes);
    });
  };
  const double t_sync = hydro(false);
  const double t_rec = hydro(true);
  std::printf("synchronized_hydro_seconds %.6e\nrecursive_hydro_seconds %.6e\nrecursive_over_synchronized %.4f\n",
              t_sync, t_rec, t_rec / t_sync);
  std::printf("governor_peak_accounted_bytes %llu\n", static_cast<unsigned long long>(peak));

  amr::ProductionAmrRegridPolicy policy;
  policy.criteria.use_particle_count = false;
  policy.criteria.gradient_threshold = 0.1;
  policy.physics.metallicity_jump = 0.5;
  policy.physics.vorticity = 0.1;
  policy.max_level = 4U;
  policy.periodic = true;
  policy.periodic_box_extent_comoving = {box, box, box};
  std::printf("amr_regrid_plan_seconds %.6e\n",
              medianSeconds(repeats, [&] { (void)amr::planProductionAmrRegrid(hierarchy, {}, policy); }));

  std::mt19937_64 rng(7U);
  std::uniform_real_distribution<double> u(0.0, box);
  std::vector<ph::AssociationSourceDescriptor> sources(20000U);
  for (std::size_t i = 0; i < sources.size(); ++i) {
    sources[i].source_id = 1U + i;
    sources[i].position_comoving = {u(rng), u(rng), u(rng)};
    sources[i].birth_years = 1.0e6 * (i % 5U);
  }
  std::printf("association_fof_seconds %.6e\nassociation_sources %zu\n", medianSeconds(repeats, [&] {
    (void)ph::discoverAssociationLeaders(sources, 0.05 * box, 2.0e6, true, {box, box, box});
  }), sources.size());

  {
    ph::FeedbackUnits units;
    units.periodic = true;
    units.box_extent_comoving = {box, box, box};
    std::vector<std::uint32_t> rows(hierarchy.cells.size());
    for (std::uint32_t r = 0; r < rows.size(); ++r) rows[r] = r;
    ph::FeedbackSourceBudget budget;
    budget.ejecta_mass_code = 1.0e-3;
    budget.supernova_count = 1.0;
    budget.supernova_energy_erg = 1.0e51;
    const double t_feedback = medianSeconds(repeats, [&] {
      ph::MechanicalFeedbackDeposition deposition(ph::MechanicalFeedbackConfig{}, units);
      deposition.setCandidates(hierarchy, rows);
      for (int s = 0; s < 64; ++s) {
        budget.source_id = 1U + static_cast<std::uint64_t>(s);
        budget.position_comoving = {u(rng), u(rng), u(rng)};
        deposition.stage(hierarchy, budget);
      }
    });
    std::printf("feedback_stage_seconds_per_source %.6e\n", t_feedback / 64.0);
  }

  std::printf("rollback_copy_seconds %.6e\nrollback_copy_bytes %llu\n",
              medianSeconds(repeats, [&] { const core::SimulationState copy = hierarchy; (void)copy; }),
              static_cast<unsigned long long>(
                  core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(hierarchy))));
  return 0;
}
