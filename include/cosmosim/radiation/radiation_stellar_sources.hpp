#pragma once

// CHUI-EXT-PHYS-REF-v1 F02: stellar ionizing-luminosity provider for RT.
//
// Data contract (external, versioned, never bundled or synthesized):
// schema `chui_ssp_group_luminosity_v1`, a text table with metadata
//   # schema: chui_ssp_group_luminosity_v1
//   # source: <population synthesis code + version, e.g. BPASS v2.2.1>
//   # imf: <declared IMF name and mass range>
//   # group_edges_ev: e0, e1, ..., eG       (must equal the RT group edges)
//   # metallicities: Z0, Z1, ...             (mass fractions, > 0, ascending)
//   # columns: age_yr metallicity group L_erg_s_per_msun Q_photons_s_per_msun
// and one row per (age, Z, group) on the full tensor grid; ages ascend and
// start at 0 (newborn populations are covered without extrapolation).
// L is the group luminosity per solar mass of INITIAL (birth) stellar mass.
//
// Model: a star particle of birth mass M_* and metallicity Z emits, over an
// accepted interval [t0, t1] with ages tau0 = t0 - t_form, tau1 = t1 - t_form,
//   E_g = M_* int_{tau0}^{tau1} L_g(tau, Z) dtau,
// where L_g is the table's piecewise-linear interpolant in age (exact
// integral over the union of nodes) and linear in log10 Z between bracketing
// metallicities. Ages beyond the table, Z outside the table, or group edges
// that differ from the RT definition are REFUSED (no extrapolation, no
// substituted spectrum). Emission is a pure function of age, so disjoint
// accepted intervals emit each photon exactly once (no pending state).
// Deposition: the energy rate E_g / dt is placed in the radiation record cell
// containing the star (bucket-indexed leaf location), isotropic in the source
// frame, converted to the comoving a^4 storage exactly like the AGN sources.
// Photon number: RT carries energy; Q_g of the table is kept as a diagnostic
// and the ledger reports the ratio Q_table / (L / <h nu>_g) so the declared
// reference within-group shape of the RT definition is auditable.
// Formation-time authority: StarParticleSidecar::formation_scale_factor, so
// stellar RT sources are admitted only in cosmological comoving modes.

#include <string>
#include <vector>

#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/physics/rt_bh_agn_frame.hpp"
#include "cosmosim/radiation/radiation_runtime.hpp"
#include "cosmosim/radiation/radiation_spectrum.hpp"

namespace cosmosim::radiation {

inline constexpr const char* k_ssp_group_luminosity_schema = "chui_ssp_group_luminosity_v1";

struct StellarLuminosityTable {
  std::string sha256;
  std::string source;
  std::string imf;
  std::vector<double> group_edges_ev;
  std::vector<double> age_yr;
  std::vector<double> metallicity;
  std::vector<double> l_erg_s_per_msun;      // [(age * nZ + z) * nG + g]
  std::vector<double> q_photons_s_per_msun;  // same layout
  [[nodiscard]] std::size_t groupCount() const noexcept { return group_edges_ev.empty() ? 0U : group_edges_ev.size() - 1U; }
  // Group luminosity per solar mass at (age, Z); throws outside the domain.
  [[nodiscard]] double luminosity(std::size_t group, double age_yr, double metallicity) const;
  [[nodiscard]] double photonRate(std::size_t group, double age_yr, double metallicity) const;
  // Exact integral of the age-interpolant over [age0, age1] (erg per Msun).
  [[nodiscard]] double integratedEnergy(std::size_t group, double age0_yr, double age1_yr, double metallicity) const;
};

[[nodiscard]] StellarLuminosityTable parseStellarLuminosityTable(const std::string& bytes, const std::string& origin);
[[nodiscard]] StellarLuminosityTable loadStellarLuminosityTable(const std::string& path, const std::string& expected_sha256);
// Throws unless the table edges equal the RT definition edges (relative 1e-9).
void requireMatchingGroupEdges(const StellarLuminosityTable& table, const RadiationSpectralDefinition& definition);

struct StellarRadiationSourceReport {
  std::uint64_t sources = 0U;
  double emitted_energy_code = 0.0;
  double table_photons = 0.0;          // sum of the table Q integrals (diagnostic)
  double transported_photon_equivalent = 0.0;  // sum E_g / <h nu>_g of the RT definition
};

struct StellarSourceInputs {
  const core::LambdaCdmBackground* background = nullptr;
  core::UnitSystem units;
  double scale_factor_begin = 1.0;
  double scale_factor_end = 1.0;
  double dt_proper_code = 0.0;
};

[[nodiscard]] StellarRadiationSourceReport buildStellarRadiationSources(
    const core::SimulationState& state, const StellarLuminosityTable& table,
    const RadiationSpectralDefinition& definition, const StellarSourceInputs& inputs,
    const physics::ProperFrameContext& frame, const physics::PeriodicBox& box,
    std::vector<RadiationSourceDeposit>& deposits);

}  // namespace cosmosim::radiation
