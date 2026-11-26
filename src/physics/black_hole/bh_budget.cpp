#include "cosmosim/physics/black_hole/bh_budget.hpp"

#include <cmath>
#include <stdexcept>

namespace cosmosim::physics {

EddingtonConstantsCode EddingtonConstantsCode::fromFrame(const ProperFrameContext& frame, double proton_mass_si,
                                                         double thomson_cross_section_si) {
  EddingtonConstantsCode c;
  c.newton_g = frame.newtonGCode();
  c.proton_mass = proton_mass_si / frame.mass_si_per_code;
  c.thomson_cross_section = thomson_cross_section_si / (frame.length_si_per_code * frame.length_si_per_code);
  c.speed_of_light = frame.speedOfLightCode();
  if (!(c.newton_g > 0.0) || !(c.proton_mass > 0.0) || !(c.thomson_cross_section > 0.0) ||
      !(c.speed_of_light > 0.0)) {
    throw std::invalid_argument("EddingtonConstantsCode: non-positive constant in code units");
  }
  return c;
}

double eddingtonLuminosityCode(double m, const EddingtonConstantsCode& c) noexcept {
  return 4.0 * core::constants::k_pi * c.newton_g * m * c.proton_mass * c.speed_of_light / c.thomson_cross_section;
}

double eddingtonReferenceRateCode(double m, double eps_ref, const EddingtonConstantsCode& c) noexcept {
  return eddingtonLuminosityCode(m, c) / (eps_ref * c.speed_of_light * c.speed_of_light);
}

BhEmissionBudget computeBhEmissionBudget(BhMassLedger ledger, double dm0, double dt, double bh_mass,
                                         const BhChannelEfficiencies& eff, double eps_ref,
                                         const EddingtonConstantsCode& c, double spin_extracted) {
  if (!(dm0 >= 0.0) || !std::isfinite(dm0) || !(dt > 0.0) || !(bh_mass > 0.0) || !(eff.total() < 1.0) ||
      eff.radiative < 0.0 || eff.thermal < 0.0 || eff.kick < 0.0 || eff.wind < 0.0 || eff.jet < 0.0 ||
      !(spin_extracted >= 0.0)) {
    throw std::invalid_argument("computeBhEmissionBudget: invalid inner mass, interval or efficiencies");
  }
  BhEmissionBudget b;
  const double c2 = c.speed_of_light * c.speed_of_light;
  const double rest = dm0 * c2;
  b.inner_mass_code = dm0;
  b.radiative_energy_code = eff.radiative * rest;
  b.thermal_energy_code = eff.thermal * rest;
  b.kick_energy_code = eff.kick * rest;
  b.wind_energy_code = eff.wind * rest;
  b.jet_energy_code = eff.jet * rest;
  b.spin_extracted_energy_code = spin_extracted;
  if (ledger == BhMassLedger::kLegacyNewtonian) {
    if (spin_extracted != 0.0) {
      throw std::invalid_argument("legacy_newtonian ledger cannot carry spin extraction (mixed ledgers)");
    }
    b.horizon_mass_increment_code = dm0;
    b.mass_defect_code = 0.0;
  } else {
    b.mass_defect_code = eff.total() * dm0 + spin_extracted / c2;
    b.horizon_mass_increment_code = (1.0 - eff.total()) * dm0 - spin_extracted / c2;
  }
  b.luminosity_code = b.radiative_energy_code / dt;
  const double l_edd = eddingtonLuminosityCode(bh_mass, c);
  b.luminosity_eddington_ratio = l_edd > 0.0 ? b.luminosity_code / l_edd : 0.0;
  const double mdot_ref = eddingtonReferenceRateCode(bh_mass, eps_ref, c);
  b.rate_eddington_reference_ratio = mdot_ref > 0.0 ? (dm0 / dt) / mdot_ref : 0.0;
  return b;
}

}  // namespace cosmosim::physics
