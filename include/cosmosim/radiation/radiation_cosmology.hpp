#pragma once

// CHUI-RT-BH-AGN-REF-v1, RT companion section 4 (RT05): cosmological
// representation of the radiation moments.
//
// Canonical moment densities E_comov = a^4 E_phys, F_comov = a^4 F_phys at
// PHYSICAL group edges. Group-integrated homogeneous FLRW equation
//   d_t E_g + 4 H E_g + a^-1 div F_g = H [nu E_nu]_{nu_l}^{nu_u}
// becomes, in comoving densities, a pure edge transfer
//   d_t E_comov,g = a^4 H [nu E_nu]_{nu_l}^{nu_u}
// with the conservative UPWIND edge flux (photons drift to lower nu): the
// energy crossing the lower edge of group g per unit time is
//   H E_l phi_low,g E_g     (phi_low from the declared within-group shape).
// It is integrated exactly per substep (exponential, positivity preserving)
// from the highest group down; the lowest edge loses photons out of the
// modelled range (ledgered boundary leakage). The same factor acts on the
// flux of each group. Summing groups cancels the internal transfers; the
// comoving total changes only by the bottom-edge leakage. The physical
// expansion work (redshift loss beyond dilution) is ledgered as
// H E_phys V_phys dt. The hydro a^3 matter conversion is NEVER used here.
// Grey (bolometric) runs have no edges: E_comov is invariant under expansion.

#include <span>

#include "cosmosim/radiation/radiation_spectrum.hpp"

namespace cosmosim::radiation {

[[nodiscard]] inline double radiationPhysicalFromComoving(double moment_comoving, double scale_factor) noexcept {
  const double a2 = scale_factor * scale_factor;
  return moment_comoving / (a2 * a2);
}
[[nodiscard]] inline double radiationComovingFromPhysical(double moment_physical, double scale_factor) noexcept {
  const double a2 = scale_factor * scale_factor;
  return moment_physical * a2 * a2;
}

struct RadiationRedshiftLedger {
  double leaked_comoving_energy = 0.0;  // sum over cells of leaked E_comov * V_comov
};

// Apply the exact edge transfer to one cell's group moments (comoving).
void applyFrequencyRedshiftCell(const RadiationSpectralDefinition& spectrum, double hubble_code, double dt,
                                std::span<double> e_comoving, std::span<double> gx, std::span<double> gy,
                                std::span<double> gz, double cell_volume_comoving, RadiationRedshiftLedger& ledger);

}  // namespace cosmosim::radiation
