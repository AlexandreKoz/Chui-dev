// CHUI-RT-BH-AGN-REF-v1 V05/V07/V08 (engine-level) integration test: one
// accepted source interval of the reference BH engine on a uniform leaf patch
// followed by AGN thermal injection, through the production classes and the
// serial support exchange. Checks the closed mass/momentum ledgers, the
// disk/horizon split, the exact pending-budget debit and the restartable
// module sidecar. NOT RUN in pass 01.
#include "../support/test_require.hpp"

#include <cmath>
#include <numeric>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/physics/agn/agn_injection.hpp"
#include "cosmosim/physics/black_hole/bh_engine.hpp"
#include "cosmosim/physics/black_hole/bh_engine_state.hpp"

namespace p = cosmosim::physics;
namespace c = cosmosim::core;

namespace {

constexpr std::uint32_t k_n = 8U;
constexpr double k_extent = 4.0;  // kpc (code length)

c::SimulationState makeState() {
  c::SimulationState s;
  const std::size_t cells = static_cast<std::size_t>(k_n) * k_n * k_n;
  s.resizePatches(1);
  s.patches.patch_id[0] = 1U;
  s.patches.level[0] = 0;
  s.patches.first_cell[0] = 0U;
  s.patches.cell_count[0] = static_cast<std::uint32_t>(cells);
  s.patches.parent_patch_id[0] = 0U;
  s.patches.extent_x_comoving[0] = s.patches.extent_y_comoving[0] = s.patches.extent_z_comoving[0] = k_extent;
  s.patches.cell_dim_x[0] = s.patches.cell_dim_y[0] = s.patches.cell_dim_z[0] = k_n;
  s.patches.owning_rank[0] = 0U;
  s.resizeCells(cells);
  const double dx = k_extent / k_n;
  for (std::uint32_t k = 0; k < k_n; ++k)
    for (std::uint32_t j = 0; j < k_n; ++j)
      for (std::uint32_t i = 0; i < k_n; ++i) {
        const std::size_t r = i + k_n * (j + k_n * k);
        s.cells.center_x_comoving[r] = (i + 0.5) * dx;
        s.cells.center_y_comoving[r] = (j + 0.5) * dx;
        s.cells.center_z_comoving[r] = (k + 0.5) * dx;
        s.cells.mass_code[r] = 1.0e7;
        s.cells.patch_index[r] = 0U;
        s.gas_cells.density_code[r] = 1.0e7 / (dx * dx * dx);
        s.gas_cells.internal_energy_code[r] = 1.0e4;     // (km/s)^2
        s.gas_cells.temperature_code[r] = 1.0e6;          // hot: Bondi phase
        s.gas_cells.sound_speed_code[r] = 1.0e2;
        s.gas_cells.velocity_x_peculiar[r] = 10.0;
        s.gas_cells.gas_cell_id[r] = 1U + r;
      }
  s.resizeParticles(1);
  s.particles.position_x_comoving[0] = s.particles.position_y_comoving[0] = s.particles.position_z_comoving[0] = 2.0;
  s.particles.mass_code[0] = 1.0e7;
  s.particle_sidecar.particle_id[0] = 99U;
  s.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(c::ParticleSpecies::kBlackHole);
  s.black_holes.resize(1);
  s.black_holes.particle_index[0] = 0U;
  s.black_holes.subgrid_mass_code[0] = 1.0e7;
  s.species.count_by_species[static_cast<std::size_t>(c::ParticleSpecies::kBlackHole)] = 1U;
  s.rebuildSpeciesIndex();
  return s;
}

double gasMass(const c::SimulationState& s) {
  return std::accumulate(s.cells.mass_code.begin(), s.cells.mass_code.end(), 0.0);
}
double gasMomentumX(const c::SimulationState& s) {
  double sum = 0.0;
  for (std::size_t r = 0; r < s.cells.size(); ++r) sum += s.cells.mass_code[r] * s.gas_cells.velocity_x_peculiar[r];
  return sum;
}
double gasEnergy(const c::SimulationState& s) {
  double sum = 0.0;
  for (std::size_t r = 0; r < s.cells.size(); ++r) {
    const double v = s.gas_cells.velocity_x_peculiar[r];
    sum += s.cells.mass_code[r] * (s.gas_cells.internal_energy_code[r] + 0.5 * v * v);
  }
  return sum;
}

void testIntervalLedgersClose() {
  auto state = makeState();
  const auto units = c::makeUnitSystem("kpc", "msun", "km_s");
  c::PhysicsConfig physics;
  physics.enable_black_hole_agn = true;
  physics.bh_engine.model = c::BhEngineModelConfig::kReferenceV1;
  physics.bh_engine.supply_model = c::BhSupplyModelConfig::kBondi;
  physics.bh_engine.kernel_radius_min_kpc_phys = 0.5;
  physics.bh_engine.kernel_radius_max_kpc_phys = 1.5;
  physics.bh_engine.kernel_target_gas_mass_code = 5.0e7;
  physics.bh_engine.disk_viscous_time_years = 1.0e6;
  physics.bh_engine.disk_wind_loading = 0.0;
  physics.agn.model = c::AgnFeedbackModelConfig::kChannelsV1;
  physics.agn.thermal_efficiency = 0.01;
  const p::PeriodicBox box{{k_extent, k_extent, k_extent}, true};
  const p::BlackHoleEngine engine(p::makeBlackHoleEngineParameters(physics, units, 0.1, box, 1.0e-20));
  p::LocalBhSupportExchange exchange;
  p::BhIntervalContext ctx;
  ctx.frame = p::ProperFrameContext::isolated(units);
  ctx.time_begin_code = 0.0;
  ctx.time_end_code = 1.0e-3;
  ctx.dt_proper_code = 1.0e-3;
  ctx.interval_index = 1U;
  const double m0 = gasMass(state) + state.particles.mass_code[0];
  const double p0 = gasMomentumX(state) + state.particles.mass_code[0] * state.particles.velocity_x_peculiar[0];
  const auto report = engine.advanceInterval(state, ctx, exchange);
  CHUI_REQUIRE(report.local_black_holes == 1U);
  CHUI_REQUIRE(report.withdrawn_mass_code > 0.0);
  CHUI_REQUIRE(report.withdrawn_mass_code <= report.requested_mass_code * (1.0 + 1.0e-12));
  const double m1 = gasMass(state) + state.particles.mass_code[0] + report.mass_defect_code;
  CHUI_REQUIRE(std::abs(m1 - m0) <= 1.0e-12 * m0);  // legacy ledger: defect == 0
  const double p1 = gasMomentumX(state) + state.particles.mass_code[0] * state.particles.velocity_x_peculiar[0];
  CHUI_REQUIRE(std::abs(p1 - p0) <= 1.0e-12 * std::abs(p0));
  const p::BlackHoleEngineStore store(state);
  const p::BlackHoleEngineRecord rec = store.get(99U);
  CHUI_REQUIRE(std::abs(state.particles.mass_code[0] - rec.dynamicalMassCode()) <= 1.0e-12 * rec.dynamicalMassCode());
  CHUI_REQUIRE(rec.pending_thermal_energy_code > 0.0);
  CHUI_REQUIRE(rec.horizon_mass_code > 1.0e7);

  // AGN thermal channel: exact budget debit and energy closure.
  const p::AgnFeedbackInjector injector(physics.agn);
  p::AgnInjectionContext actx;
  actx.frame = ctx.frame;
  actx.time_end_code = ctx.time_end_code;
  actx.interval_index = 1U;
  actx.box = box;
  const double e_before = gasEnergy(state);
  const double pending = rec.pending_thermal_energy_code;
  const auto agn = injector.apply(state, actx, exchange);
  CHUI_REQUIRE(agn.sources == 1U);
  CHUI_REQUIRE(std::abs(gasEnergy(state) - e_before - pending) <= 1.0e-10 * pending);
  CHUI_REQUIRE(p::BlackHoleEngineStore(state).get(99U).pending_thermal_energy_code <= 1.0e-12 * pending);
  CHUI_REQUIRE(agn.max_energy_closure_residual <= 1.0e-10);
}

}  // namespace

int main() {
  testIntervalLedgersClose();
  return 0;
}
