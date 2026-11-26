// CHUI-RT-BH-AGN-REF-v1 V01/V02/V03/V04 radiation core tests: M1 closure
// limits and wave-speed bounds, joint admissible limiter, spectral contract
// identities, frequency-redshift conservation, free-streaming pulse speed,
// grey implicit LTE relaxation conservation, PD-ARS runtime vacuum/energy
// invariants, AMR transfer conservation, VET projection, MC reproducibility.
// NOT RUN in pass 01.
#include "../support/test_require.hpp"

#include <cmath>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/radiation/radiation_amr.hpp"
#include "cosmosim/radiation/radiation_angular.hpp"
#include "cosmosim/radiation/radiation_closure.hpp"
#include "cosmosim/radiation/radiation_cosmology.hpp"
#include "cosmosim/radiation/radiation_matter_coupling.hpp"
#include "cosmosim/radiation/radiation_runtime.hpp"
#include "cosmosim/radiation/radiation_spectrum.hpp"

namespace r = cosmosim::radiation;
namespace c = cosmosim::core;

namespace {

bool close(double a, double b, double rel, double abs_floor = 1.0e-300) {
  return std::abs(a - b) <= std::max(rel * std::max(std::abs(a), std::abs(b)), abs_floor);
}

void testClosureLimits() {
  CHUI_REQUIRE(close(r::m1EddingtonFactor(0.0), 1.0 / 3.0, 1.0e-15));
  CHUI_REQUIRE(close(r::m1EddingtonFactor(1.0), 1.0, 1.0e-15));
  const auto iso = r::m1PressureTensor({3.0, {0.0, 0.0, 0.0}}, 1.0);
  CHUI_REQUIRE(close(iso[0], 1.0, 1.0e-15) && iso[1] == 0.0 && close(iso[8], 1.0, 1.0e-15));
  const auto beam = r::m1PressureTensor({2.0, {0.0, 2.0, 0.0}}, 1.0);
  CHUI_REQUIRE(close(beam[4], 2.0, 1.0e-12) && std::abs(beam[0]) < 1.0e-12);
  // Wave speeds: isotropic +-c/sqrt3; free-streaming beam along x travels at c.
  const auto s_iso = r::m1WaveSpeeds({1.0, {0.0, 0.0, 0.0}}, 0, 1.0);
  CHUI_REQUIRE(close(s_iso.max, 1.0 / std::sqrt(3.0), 1.0e-12));
  const auto s_beam = r::m1WaveSpeeds({1.0, {0.999999, 0.0, 0.0}}, 0, 1.0);
  CHUI_REQUIRE(s_beam.max <= 1.0 && s_beam.max > 0.99);
  // Rotational invariance of the bounds: beam along y probed along y.
  const auto s_beam_y = r::m1WaveSpeeds({1.0, {0.0, 0.999999, 0.0}}, 1, 1.0);
  CHUI_REQUIRE(close(s_beam.max, s_beam_y.max, 1.0e-12));
  // Joint limiter keeps both face states admissible.
  const r::M1State base{1.0, {0.9, 0.0, 0.0}};
  const r::M1State delta{-0.5, {0.5, 0.0, 0.0}};
  const double theta = r::m1AdmissibleFraction(base, delta, 1.0);
  CHUI_REQUIRE(theta < 1.0 && theta > 0.0);
  const r::M1State lim{base.energy + theta * delta.energy, {base.flux[0] + theta * delta.flux[0], 0.0, 0.0}};
  CHUI_REQUIRE(r::m1Admissible(lim, 1.0, 1.0e-10));
}

void testSpectralIdentities() {
  r::SpectralShape bb;
  bb.kind = r::SpectralShapeKind::kBlackbody;
  bb.blackbody_temperature_k = 1.0e5;
  const auto def = r::buildRadiationSpectralDefinition("ionizing3", 100.0, bb);
  CHUI_REQUIRE(def.groupCount() == 3U && def.sha256.size() == 64U);
  CHUI_REQUIRE(close(r::photoionizationCrossSectionCm2(r::RadiationAbsorber::kHI, r::k_ionization_potential_hi_ev),
                     6.30e-18, 1.0e-12));
  CHUI_REQUIRE(r::photoionizationCrossSectionCm2(r::RadiationAbsorber::kHI, 13.0) == 0.0);
  for (const auto& g : def.groups) {
    CHUI_REQUIRE(g.mean_photon_energy_ev > g.lower_ev && g.mean_photon_energy_ev < g.upper_ev);
    for (std::size_t a = 0; a < r::k_radiation_absorber_count; ++a) {
      // sigma_N (E/<h nu>) eps == sigma_E E (exact by construction).
      if (g.sigma_number_cm2[a] > 0.0) {
        CHUI_REQUIRE(close(g.sigma_number_cm2[a] / g.mean_photon_energy_ev * g.absorbed_energy_per_ionization_ev[a],
                           g.sigma_energy_cm2[a], 1.0e-12));
      }
    }
  }
  const auto frac = r::radiationGroupLuminosityFractions(def, bb);
  CHUI_REQUIRE(close(frac[0] + frac[1] + frac[2], 1.0, 1.0e-14));
}

void testRedshiftTransferConservesInteriorEnergy() {
  r::SpectralShape flat;
  const auto def = r::buildRadiationSpectralDefinition("ionizing3", 100.0, flat);
  std::vector<double> e{1.0, 2.0, 3.0}, gx{0.1, 0.2, 0.3}, gy(3, 0.0), gz(3, 0.0);
  r::RadiationRedshiftLedger ledger;
  r::applyFrequencyRedshiftCell(def, 1.0e-3, 1.0, e, gx, gy, gz, 1.0, ledger);
  CHUI_REQUIRE(close(e[0] + e[1] + e[2] + ledger.leaked_comoving_energy, 6.0, 1.0e-14));
  CHUI_REQUIRE(e[2] < 3.0 && ledger.leaked_comoving_energy > 0.0);
}

c::SimulationState makeRadiationState(std::uint32_t n, double extent) {
  c::SimulationState s;
  s.resizePatches(1);
  s.patches.patch_id[0] = 1U;
  s.patches.first_cell[0] = 0U;
  s.patches.cell_count[0] = n * n * n;
  s.patches.extent_x_comoving[0] = s.patches.extent_y_comoving[0] = s.patches.extent_z_comoving[0] = extent;
  s.patches.cell_dim_x[0] = s.patches.cell_dim_y[0] = s.patches.cell_dim_z[0] = static_cast<std::uint16_t>(n);
  s.patches.owning_rank[0] = 0U;
  s.resizeCells(n * n * n);
  return s;
}

void testVacuumAndPulseConservation() {
  auto s = makeRadiationState(16U, 1.0);
  r::SpectralShape flat;
  r::RadiationRuntimeSettings rs;
  rs.spectrum = r::buildRadiationSpectralDefinition("grey", 0.0, flat);
  rs.config.enabled = true;
  rs.config.spectral_profile = c::RtSpectralProfileConfig::kGrey;
  rs.coupling.c_code = 1.0;
  rs.coupling.c_hat_code = 1.0;
  rs.transport.c_hat_code = 1.0;
  rs.transport.boundary = c::RtBoundaryConfig::kPeriodic;
  rs.transport.domain_extent_comoving = {1.0, 1.0, 1.0};
  r::allocateRadiationState(s, rs.spectrum, 1.0, 1.0, 0U);
  CHUI_REQUIRE(s.radiation.isConsistent());
  const r::RadiationRuntime runtime(rs);
  r::RadiationIntervalInputs in;
  in.dt_proper_code = 0.1;
  in.frame_end = cosmosim::physics::ProperFrameContext::isolated(c::makeUnitSystem("kpc", "msun", "km_s"));
  // Vacuum stays vacuum.
  const auto vac = runtime.advanceInterval(s, in);
  CHUI_REQUIRE(vac.substeps > 0U);
  for (const double v : s.radiation.energy_density_comoving) CHUI_REQUIRE(v == 0.0);
  // Gaussian pulse in a periodic box: total energy conserved by the
  // conservative fluxes (no sources/sinks), state stays admissible.
  double e0 = 0.0;
  for (std::size_t i = 0; i < s.radiation.energy_density_comoving.size(); ++i) {
    const std::uint32_t ix = static_cast<std::uint32_t>(i % 16U);
    const double x = (ix + 0.5) / 16.0 - 0.5;
    s.radiation.energy_density_comoving[i] = 1.0e-6 + std::exp(-x * x / 0.005);
    e0 += s.radiation.energy_density_comoving[i];
  }
  (void)runtime.advanceInterval(s, in);
  double e1 = 0.0;
  for (std::size_t i = 0; i < s.radiation.energy_density_comoving.size(); ++i) {
    e1 += s.radiation.energy_density_comoving[i];
    CHUI_REQUIRE(r::m1Admissible({s.radiation.energy_density_comoving[i],
                                  {s.radiation.flux_x_comoving[i], s.radiation.flux_y_comoving[i],
                                   s.radiation.flux_z_comoving[i]}},
                                 1.0, 1.0e-10));
  }
  CHUI_REQUIRE(close(e0, e1, 1.0e-12));
}

void testGreyImplicitRelaxationConservesEnergy() {
  r::RadiationMatterCell cell;
  cell.density = 1.0;
  cell.specific_internal_energy = 2.0;
  cell.temperature_per_specific_energy = 1.0;  // T = u
  r::RadiationCouplingConstants k;
  k.c_code = 10.0;
  k.c_hat_code = 10.0;
  k.kappa_absorption_code = 5.0;
  k.radiation_constant_code = 1.0e-2;
  double e = 0.01;
  std::array<double, 3> g{0.0, 0.0, 0.0};
  const double total0 = 2.0 + e;
  const auto res = r::solveGreyThermalImplicit(cell, e, g, 10.0, k);
  CHUI_REQUIRE(res.converged);
  CHUI_REQUIRE(close(cell.specific_internal_energy + e, total0, 1.0e-9));
  // Long implicit step approaches LTE: E -> a T^4.
  const double t = cell.specific_internal_energy;
  CHUI_REQUIRE(close(e, 1.0e-2 * t * t * t * t, 1.0e-2));
}

void testAmrTransferConservative() {
  const r::M1State parent{2.0, {0.5, 0.1, 0.0}};
  std::array<std::array<r::M1State, 2>, 3> nb{};
  for (auto& axis : nb) {
    axis[0] = {1.0, {0.2, 0.0, 0.0}};
    axis[1] = {3.0, {0.9, 0.1, 0.0}};
  }
  const auto kids = r::prolongRadiationParent(parent, nb, 1.0);
  const auto back = r::restrictRadiationChildren(kids);
  CHUI_REQUIRE(close(back.energy, parent.energy, 1.0e-14));
  CHUI_REQUIRE(close(back.flux[0], parent.flux[0], 1.0e-14));
  for (const auto& k : kids) CHUI_REQUIRE(r::m1Admissible(k, 1.0, 1.0e-12));
}

void testVetAndMonteCarloReproducibility() {
  // Isotropic samples -> f = I/3.
  std::vector<std::array<double, 3>> dirs{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
  const auto vet = r::vetFromAngularSamples(dirs, std::vector<double>(6, 1.0), std::vector<double>(6, 1.0));
  CHUI_REQUIRE(vet.defined && close(vet.f[0], 1.0 / 3.0, 1.0e-12) && close(vet.f[0] + vet.f[4] + vet.f[8], 1.0, 1.0e-12));
  CHUI_REQUIRE(r::radiationCounterRandom(1U, 2U, 3U) == r::radiationCounterRandom(1U, 2U, 3U));
  CHUI_REQUIRE(r::radiationCounterRandom(1U, 2U, 3U) != r::radiationCounterRandom(1U, 2U, 4U));
}

}  // namespace

int main() {
  testClosureLimits();
  testSpectralIdentities();
  testRedshiftTransferConservesInteriorEnergy();
  testVacuumAndPulseConservation();
  testGreyImplicitRelaxationConservesEnergy();
  testAmrTransferConservative();
  testVetAndMonteCarloReproducibility();
  return 0;
}
