#pragma once

// CHUI-RT-BH-AGN-REF-v1, BH companion section 6 (BH06) and AGN companion
// section 2 (AGN01): Eddington convention, mass ledgers and the integrated
// per-channel emission budget issued by the BH/disk engine.
//
// Eddington (ionized hydrogen, Thomson): L_Edd = 4 pi G M m_p c / sigma_T,
// reference inflow rate mdot_Edd,ref = L_Edd/(eps_ref c^2) (eps_ref = 0.1
// proposed reference convention). Rate cap and luminosity are distinct:
// both mdot/mdot_Edd,ref and L/L_Edd are reported.
//
// Ledgers (mutually exclusive, versioned):
//  legacy_newtonian    : dM_bh = dm0. Feedback energy is phenomenological and
//                        is NOT a rest-energy conservation claim.
//  relativistic_defect : dm0 is the baryonic mass reaching the inner engine
//                        AFTER disk winds (wind rest mass was never part of
//                        dm0, so it is not subtracted). With eps_acc the sum
//                        of all outgoing accretion-powered fractions,
//                          dM_bh = (1 - eps_acc) dm0 - E_spin_extracted / c^2.
//                        The released energy eps_acc dm0 c^2 leaves as the
//                        channel budgets below; it is never added as rest
//                        mass into the Newtonian hydro energy lanes.

#include <array>
#include <cstdint>

#include "cosmosim/physics/rt_bh_agn_frame.hpp"

namespace cosmosim::physics {

enum class BhMassLedger : std::uint8_t { kLegacyNewtonian = 0, kRelativisticDefect = 1 };

struct BhChannelEfficiencies {
  double radiative = 0.0;  // eps_r
  double thermal = 0.0;
  double kick = 0.0;
  double wind = 0.0;
  double jet = 0.0;
  [[nodiscard]] double total() const noexcept { return radiative + thermal + kick + wind + jet; }
};

struct EddingtonConstantsCode {
  double newton_g = 0.0;
  double proton_mass = 0.0;
  double thomson_cross_section = 0.0;
  double speed_of_light = 0.0;
  [[nodiscard]] static EddingtonConstantsCode fromFrame(const ProperFrameContext& frame, double proton_mass_si,
                                                        double thomson_cross_section_si);
};

[[nodiscard]] double eddingtonLuminosityCode(double bh_mass_code, const EddingtonConstantsCode& c) noexcept;
[[nodiscard]] double eddingtonReferenceRateCode(double bh_mass_code, double reference_efficiency,
                                                const EddingtonConstantsCode& c) noexcept;

// Integrated budget of one BH for one accepted source interval.
struct BhEmissionBudget {
  std::uint64_t bh_particle_id = 0U;
  std::uint64_t interval_index = 0U;
  double interval_begin_code = 0.0;
  double interval_end_code = 0.0;
  double inner_mass_code = 0.0;        // dm0
  double radiative_energy_code = 0.0;  // L_rad dt (accretion powered)
  double thermal_energy_code = 0.0;
  double kick_energy_code = 0.0;
  double wind_energy_code = 0.0;
  double jet_energy_code = 0.0;
  double spin_extracted_energy_code = 0.0;  // P1 (separate reserve)
  double horizon_mass_increment_code = 0.0;
  double mass_defect_code = 0.0;       // rest-mass equivalent removed by the relativistic ledger
  double luminosity_code = 0.0;
  double luminosity_eddington_ratio = 0.0;
  double rate_eddington_reference_ratio = 0.0;
};

// Pure ledger application: given dm0 and efficiencies, the budget and
// horizon increment for the selected ledger.
[[nodiscard]] BhEmissionBudget computeBhEmissionBudget(BhMassLedger ledger, double inner_mass_code, double dt_code,
                                                       double bh_mass_code, const BhChannelEfficiencies& eff,
                                                       double reference_efficiency,
                                                       const EddingtonConstantsCode& constants,
                                                       double spin_extracted_energy_code = 0.0);

}  // namespace cosmosim::physics
