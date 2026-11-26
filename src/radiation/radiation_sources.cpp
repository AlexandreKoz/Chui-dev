#include "cosmosim/radiation/radiation_sources.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cosmosim/physics/black_hole/bh_engine_state.hpp"
#include "cosmosim/physics/black_hole/bh_supply.hpp"

namespace cosmosim::radiation {

AgnRadiationSourceReport buildAgnRadiationSources(core::SimulationState& state, const std::vector<double>& frac,
                                                  const physics::ProperFrameContext& frame,
                                                  const physics::PeriodicBox& box, double dt,
                                                  std::vector<RadiationSourceDeposit>& deposits) {
  AgnRadiationSourceReport report;
  if (!(dt > 0.0)) throw std::invalid_argument("buildAgnRadiationSources: dt must be > 0");
  const core::RadiationMomentSoa& rad = state.radiation;
  if (frac.size() != rad.group_count) throw std::invalid_argument("buildAgnRadiationSources: group fraction mismatch");
  physics::BlackHoleEngineStore store(state);
  const double a = frame.scale_factor;
  const double a4 = a * a * a * a;
  for (std::size_t b = 0; b < state.black_holes.size(); ++b) {
    const std::uint32_t p = state.black_holes.particle_index[b];
    const std::uint64_t id = state.particle_sidecar.particle_id[p];
    if (!store.has(id)) continue;
    physics::BlackHoleEngineRecord rec = store.get(id);
    if (!(rec.pending_radiation_energy_code > 0.0)) continue;
    ++report.sources;
    const std::array<double, 3> x{state.particles.position_x_comoving[p], state.particles.position_y_comoving[p],
                                  state.particles.position_z_comoving[p]};
    const double h = rec.last_kernel_radius_proper_code;
    struct W {
      std::uint32_t record;
      std::uint32_t cell;
      double omega;
      double volume_comoving;
    };
    std::vector<W> weights;
    double sum = 0.0;
    for (std::uint32_t r = 0; r < rad.patches.size() && h > 0.0; ++r) {
      const core::RadiationPatchRecord& rr = rad.patches[r];
      const std::array<double, 3> w{rr.cell_width_x_comoving, rr.cell_width_y_comoving, rr.cell_width_z_comoving};
      // Coarse record rejection: skip records whose box is farther than h.
      std::array<double, 3> centre{rr.origin_x_comoving + 0.5 * w[0] * rr.nx, rr.origin_y_comoving + 0.5 * w[1] * rr.ny,
                                   rr.origin_z_comoving + 0.5 * w[2] * rr.nz};
      const auto dc = box.delta(x, centre);
      const double half_diag = 0.5 * std::sqrt(std::pow(w[0] * rr.nx, 2) + std::pow(w[1] * rr.ny, 2) + std::pow(w[2] * rr.nz, 2));
      if (frame.properLength(std::sqrt(dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2]) - half_diag) > h) continue;
      const std::array<double, 3> hw{frame.properLength(0.5 * w[0]), frame.properLength(0.5 * w[1]),
                                     frame.properLength(0.5 * w[2])};
      const double vol_c = w[0] * w[1] * w[2];
      for (std::uint32_t k = 0; k < rr.nz; ++k) {
        for (std::uint32_t j = 0; j < rr.ny; ++j) {
          for (std::uint32_t i = 0; i < rr.nx; ++i) {
            const std::array<double, 3> cc{rr.origin_x_comoving + (i + 0.5) * w[0], rr.origin_y_comoving + (j + 0.5) * w[1],
                                           rr.origin_z_comoving + (k + 0.5) * w[2]};
            const auto d = box.delta(x, cc);
            const std::array<double, 3> off{frame.properLength(d[0]), frame.properLength(d[1]), frame.properLength(d[2])};
            const double wk = physics::bhCellKernelWeight(off, hw, h);
            if (!(wk > 0.0)) continue;
            const double omega = wk * frame.properVolume(vol_c);
            weights.push_back(W{r, i + rr.nx * (j + rr.ny * k), omega, vol_c});
            sum += omega;
          }
        }
      }
    }
    if (!(sum > 0.0)) {
      ++report.held_no_support;
      rec.flags |= physics::kBhFlagPendingFeedbackHeld;
      store.put(id, rec);
      continue;
    }
    const double luminosity = rec.pending_radiation_energy_code / dt;
    for (const W& w : weights) {
      // Physical energy-density rate in the cell: L f_g omega/sum / V_phys;
      // comoving rate multiplies by a^4.
      const double v_phys = frame.properVolume(w.volume_comoving);
      for (std::uint32_t g = 0; g < rad.group_count; ++g) {
        const double rate_phys = luminosity * frac[g] * (w.omega / sum) / v_phys;
        if (rate_phys > 0.0) deposits.push_back(RadiationSourceDeposit{w.record, w.cell, g, rate_phys * a4});
      }
    }
    report.emitted_energy_code += rec.pending_radiation_energy_code;
    rec.pending_radiation_energy_code = 0.0;
    store.put(id, rec);
  }
  store.store(state);
  return report;
}

}  // namespace cosmosim::radiation
