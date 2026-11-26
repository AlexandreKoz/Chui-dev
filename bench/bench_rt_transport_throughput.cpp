// CHUI-RT-BH-AGN-REF-v1 V09 benchmark driver (NOT RUN; no performance claim):
// cell-group updates per second of the level-synchronous M1 transport and the
// PD-ARS runtime on a uniform periodic patch, for 1/3/6/9 groups, plus the
// canonical-byte scale law 4 G 8 N (PERF01). Synthetic fixture.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/radiation/radiation_runtime.hpp"

namespace r = cosmosim::radiation;

int main(int argc, char** argv) {
  const std::uint32_t n = argc > 1 ? static_cast<std::uint32_t>(std::strtoul(argv[1], nullptr, 10)) : 32U;
  for (const char* profile : {"grey", "ionizing3"}) {
    cosmosim::core::SimulationState s;
    s.resizePatches(1);
    s.patches.patch_id[0] = 1U;
    s.patches.cell_count[0] = n * n * n;
    s.patches.extent_x_comoving[0] = s.patches.extent_y_comoving[0] = s.patches.extent_z_comoving[0] = 1.0;
    s.patches.cell_dim_x[0] = s.patches.cell_dim_y[0] = s.patches.cell_dim_z[0] = static_cast<std::uint16_t>(n);
    s.resizeCells(n * n * n);
    r::SpectralShape flat;
    r::RadiationRuntimeSettings rs;
    rs.spectrum = r::buildRadiationSpectralDefinition(profile, 100.0, flat);
    rs.coupling.c_code = rs.coupling.c_hat_code = rs.transport.c_hat_code = 1.0;
    rs.transport.domain_extent_comoving = {1.0, 1.0, 1.0};
    r::allocateRadiationState(s, rs.spectrum, 1.0, 1.0, 0U);
    for (double& e : s.radiation.energy_density_comoving) e = 1.0;
    const r::RadiationRuntime runtime(rs);
    r::RadiationIntervalInputs in;
    in.dt_proper_code = 0.05;
    in.frame_end = cosmosim::physics::ProperFrameContext::isolated(cosmosim::core::makeUnitSystem("kpc", "msun", "km_s"));
    const auto t0 = std::chrono::steady_clock::now();
    const auto rep = runtime.advanceInterval(s, in);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const double updates = static_cast<double>(rep.substeps) * n * n * n * rs.spectrum.groupCount();
    std::printf("profile=%s n=%u groups=%zu substeps=%u seconds=%.6e cell_group_updates_per_s=%.6e canonical_bytes=%llu\n",
                profile, n, rs.spectrum.groupCount(), rep.substeps, sec, sec > 0.0 ? updates / sec : 0.0,
                static_cast<unsigned long long>(r::radiationCanonicalBytes(static_cast<std::uint64_t>(n) * n * n,
                                                                          static_cast<std::uint32_t>(rs.spectrum.groupCount()))));
  }
  return 0;
}
