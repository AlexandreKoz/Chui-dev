// CHUI Part 3 bounded microbenchmarks (NOT RUN — USER-RESERVED; no numbers
// are claimed from source inspection). Measures, on deterministic synthetic
// inputs sized by argv[1] (cells per line, default 4096):
//   * element-resolved diffusion: lane-face updates per second.
// (Composition-regrid and remote-exchange benchmarks are NOT yet written;
// they are listed as open in the Part 3 handoff.)
// Output: one "key value" line per metric (machine-readable, for the owner's
// regression tables in docs/profiling.md). Methodology: median of 5 repeats
// after one warm-up; report hardware and build preset alongside.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/element_diffusion.hpp"

namespace {

namespace core = cosmosim::core;
namespace ph = cosmosim::physics;

double medianSeconds(std::vector<double> samples) {
  std::sort(samples.begin(), samples.end());
  return samples[samples.size() / 2U];
}

core::SimulationState lineState(std::uint32_t n) {
  core::SimulationState state;
  state.resizeCells(n);
  std::vector<core::GasCellIdentityRecord> records;
  for (std::uint32_t r = 0; r < n; ++r) {
    state.cells.mass_code[r] = 1.0;
    state.gas_cells.gas_cell_id[r] = 1U + r;
    state.gas_cells.density_code[r] = 1.0;
    state.gas_cells.metal_mass_code[r] = r % 2U == 0U ? 0.02 : 0.001;
    records.push_back({.gas_cell_id = 1U + r, .parent_particle_id = std::nullopt, .owning_patch_id = 1U,
                       .local_cell_row = r});
  }
  state.gas_cell_identity.assign(std::move(records));
  ph::enableGasComposition(state, ph::CompositionLayout::make({.network = ph::CompositionNetworkKind::kHybridLight,
      .dust = ph::CompositionDustKind::kTwoReservoirV1, .sgs_energy = true}), {});
  return state;
}

}  // namespace

int main(int argc, char** argv) {
  const std::uint32_t n = argc > 1 ? static_cast<std::uint32_t>(std::strtoul(argv[1], nullptr, 10)) : 4096U;
  core::SimulationState state = lineState(n);
  std::vector<ph::MetalDiffusionFace> faces;
  for (std::uint32_t i = 0; i + 1U < n; ++i) {
    faces.push_back({.left_cell = i, .right_cell = i + 1U, .area_code = 1.0, .center_distance_code = 1.0});
  }
  const std::vector<double> rho_kappa(n, 0.1);
  std::vector<double> workspace;
  std::vector<double> samples;
  std::uint64_t updates = 0U;
  for (int repeat = 0; repeat < 6; ++repeat) {
    const auto t0 = std::chrono::steady_clock::now();
    const auto report = ph::advanceElementResolvedDiffusion(state.gas_composition, faces, rho_kappa, 1.0, 0.5, 64U, workspace);
    const auto t1 = std::chrono::steady_clock::now();
    updates = report.faces_evaluated * report.mixed_lane_count;
    if (repeat > 0) samples.push_back(std::chrono::duration<double>(t1 - t0).count());
  }
  std::printf("element_diffusion_lane_face_updates_per_s %.6e\n", static_cast<double>(updates) / medianSeconds(samples));
  std::printf("element_diffusion_cells %u\n", n);
  return 0;
}
