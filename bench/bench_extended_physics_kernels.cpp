// CHUI-EXT-PHYS-REF-v1 I04 benchmark driver (source only; NOT RUN, no numbers
// are claimed). Times representative reference kernels at user-selected sizes
// and prints the governed-memory formulas next to wall times so the execution
// campaign can relate cost to the documented scale laws:
//   leaf ray traversal (segments/s), Voigt evaluations/s, FD response update,
//   CR anisotropic diffusion step (cells/s), direct RIME (component-visibilities/s).
// Usage: bench_extended_physics_kernels [grid_n] [sightlines] [stations] [components]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "cosmosim/cosmic_rays/cosmic_ray_spectrum.hpp"
#include "cosmosim/cosmic_rays/cosmic_ray_transport.hpp"
#include "cosmosim/instruments/interferometry.hpp"
#include "cosmosim/lightcone/ray_geometry.hpp"
#include "cosmosim/radiation/radiation_runtime.hpp"
#include "cosmosim/spectral/atomic_data.hpp"

namespace {

double seconds(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}

}  // namespace

int main(int argc, char** argv) {
  const std::uint32_t n = argc > 1 ? static_cast<std::uint32_t>(std::atoi(argv[1])) : 64U;
  const std::size_t sightlines = argc > 2 ? static_cast<std::size_t>(std::atoi(argv[2])) : 256U;
  const std::size_t stations = argc > 3 ? static_cast<std::size_t>(std::atoi(argv[3])) : 16U;
  const std::size_t components = argc > 4 ? static_cast<std::size_t>(std::atoi(argv[4])) : 64U;
  namespace lc = cosmosim::lightcone;
  const auto geometry = lc::LeafGeometry::uniformGrid({n, n, n}, {0, 0, 0}, {1, 1, 1}, lc::RayBoundary::kPeriodic);
  auto t0 = std::chrono::steady_clock::now();
  std::size_t segments = 0;
  for (std::size_t s = 0; s < sightlines; ++s) {
    segments += geometry.traceStraightRay({(s + 0.5) / sightlines, 0.37, 0.0}, {0, 0, 1}, 1.0, 1U << 20U).size();
  }
  auto t1 = std::chrono::steady_clock::now();
  std::printf("ray_traversal segments=%zu seconds=%.6f index_bytes=%llu\n", segments, seconds(t0, t1),
              static_cast<unsigned long long>(geometry.indexBytes()));
  t0 = std::chrono::steady_clock::now();
  double sink = 0.0;
  for (int i = 0; i < 1000000; ++i) sink += cosmosim::spectral::voigtH(4.7e-4, -30.0 + 6.0e-5 * i);
  t1 = std::chrono::steady_clock::now();
  std::printf("voigt evaluations=1000000 seconds=%.6f checksum=%.6e\n", seconds(t0, t1), sink);
  namespace cr = cosmosim::cosmic_rays;
  cr::CrGrid grid;
  grid.n = {n, n, 1};
  grid.dx = {1.0 / n, 1.0 / n, 1.0};
  std::vector<double> zero(grid.cells(), 0.0), ones(grid.cells(), 1.0), e(grid.cells(), 1.0);
  const cr::CrMediumView medium{zero, zero, zero, ones, zero, zero, {}};
  cr::CrDiffusionModel model;
  model.kappa_parallel = 1.0;
  cr::CrStepLedger ledger;
  t0 = std::chrono::steady_clock::now();
  cr::diffuseExplicit(grid, e, medium, model, cr::diffusionExplicitTimestep(grid, model, 0.4), ledger);
  t1 = std::chrono::steady_clock::now();
  std::printf("cr_diffusion cells=%zu seconds=%.6f spectral_bytes_2x16_bins=%llu\n", grid.cells(), seconds(t0, t1),
              static_cast<unsigned long long>(2U * cr::CrSpectralState::bytesFor(16U, grid.cells(),
                                                                                  cr::CrSpectralRepresentation::kPiecewisePowerLaw)));
  namespace in = cosmosim::instruments;
  std::vector<in::Station> array;
  for (std::size_t s = 0; s < stations; ++s) {
    in::Station st;
    st.name = "S" + std::to_string(s);
    st.enu_m = {17.0 * s, 5.0 * s * s, 0.0};
    array.push_back(st);
  }
  std::vector<in::SkyComponent> sky;
  for (std::size_t c = 0; c < components; ++c) sky.push_back({1.0e-3 * c, -5.0e-4 * c, {{1.0, 0.0, 0.0, 0.0}}});
  in::ObservationSettings obs;
  obs.frequencies = cosmosim::spectral::FrequencyAxis({1.0e9, 1.001e9}, 1);
  obs.integration_s = 10.0;
  std::size_t batches = 0;
  t0 = std::chrono::steady_clock::now();
  const auto report = in::simulateVisibilitiesDirect(array, sky, obs, [&batches](const in::VisibilityBatch&) { ++batches; });
  t1 = std::chrono::steady_clock::now();
  std::printf("direct_rime baselines=%llu component_evaluations=%llu seconds=%.6f\n",
              static_cast<unsigned long long>(in::baselineCount(stations)),
              static_cast<unsigned long long>(report.component_evaluations), seconds(t0, t1));
  std::printf("rt_interval_workspace_bytes(n^3 cells, 3 groups, coupled)=%llu\n",
              static_cast<unsigned long long>(cosmosim::radiation::radiationIntervalWorkspaceBytes(
                  3U * static_cast<std::size_t>(n) * n * n, static_cast<std::size_t>(n) * n * n, true)));
  return 0;
}
