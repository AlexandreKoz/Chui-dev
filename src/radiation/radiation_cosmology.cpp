#include "cosmosim/radiation/radiation_cosmology.hpp"

#include <cmath>
#include <stdexcept>

namespace cosmosim::radiation {

void applyFrequencyRedshiftCell(const RadiationSpectralDefinition& sp, double h, double dt, std::span<double> e,
                                std::span<double> gx, std::span<double> gy, std::span<double> gz, double volume,
                                RadiationRedshiftLedger& ledger) {
  if (sp.grey || !(h > 0.0) || !(dt > 0.0)) return;
  const std::size_t n = sp.groupCount();
  if (e.size() != n || gx.size() != n || gy.size() != n || gz.size() != n) {
    throw std::invalid_argument("applyFrequencyRedshiftCell: group count mismatch");
  }
  // Highest group first so that energy entering a group from above is itself
  // redshifted only in the following substep (upwind in frequency).
  for (std::size_t gi = n; gi-- > 0U;) {
    const RadiationGroup& g = sp.groups[gi];
    const double rate = h * g.lower_ev * g.lower_edge_density_per_ev;
    const double leave = -std::expm1(-rate * dt);  // fraction crossing the lower edge
    const double de = leave * e[gi];
    const double dgx = leave * gx[gi];
    const double dgy = leave * gy[gi];
    const double dgz = leave * gz[gi];
    e[gi] -= de;
    gx[gi] -= dgx;
    gy[gi] -= dgy;
    gz[gi] -= dgz;
    if (gi == 0U) {
      ledger.leaked_comoving_energy += de * volume;
    } else {
      e[gi - 1U] += de;
      gx[gi - 1U] += dgx;
      gy[gi - 1U] += dgy;
      gz[gi - 1U] += dgz;
    }
  }
}

}  // namespace cosmosim::radiation
