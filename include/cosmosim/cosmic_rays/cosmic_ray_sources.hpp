#pragma once

// CHUI-EXT-PHYS-REF-v1 CR05/CR06: injection providers and loss providers.
//
// Injection: a caller that has ALREADY DEBITED a feedback budget (SN, shock,
// AGN) passes the injected energy E_inj; the spectrum p^-alpha (in momentum,
// f ~ p^-(alpha)) between p_min and p_max is discretized on the momentum
// grid and normalized so that sum_i e_i = E_inj EXACTLY (number follows from
// the same shape). The function never reads or writes a budget itself, so a
// budget can be debited only once by its owner.
//
// Losses:
//  * electrons, synchrotron + inverse Compton (Thomson limit, derived):
//      d gamma/dt = -(4/3) sigma_T c (U_B + U_rad) beta^2 gamma^2 / (m_e c^2)
//      => dp~/dt = -C p~ gamma,  C = (4/3) sigma_T (U_B + U_rad) / (m_e c),
//    U_B = B^2 / (2 mu0), U_rad = u_cmb(a) from the S01 background authority
//    (no independent CMB temperature). Klein-Nishina corrections are a named
//    omission (valid for gamma eps_cmb << m_e c^2).
//  * protons (grey): Coulomb and hadronic loss coefficients, and the
//    thermalized hadronic fraction, come ONLY from a versioned table
//    `chui_cr_loss_rates_v1` (no coefficient is hard-coded):
//      de/dt = -(Lambda_C + Lambda_h) n_e e,  gas heating = Lambda_C n_e e +
//      f_th Lambda_h n_e e,  escaping (gamma rays / secondaries / neutrinos)
//      = (1 - f_th) Lambda_h n_e e  (explicitly ledgered).
//    Missing table -> ProviderStatus::kUnavailable -> the kernel refuses;
//    a missing provider is never a zero loss.
//  * secondary e+- production, gamma-ray emissivity, CR ionization: need
//    cross-section tables; their providers report kUnavailable in pass 01.

#include <map>
#include <string>
#include <vector>

#include "cosmosim/cosmic_rays/cosmic_ray_spectrum.hpp"
#include "cosmosim/spectral/provider_status.hpp"

namespace cosmosim::cosmic_rays {

struct CrInjectionSpectrum {
  double alpha = 0.0;   // f ~ p^-alpha (momentum-space slope), declared by the model
  double p_min = 0.0;   // p~ bounds (inside the momentum grid)
  double p_max = 0.0;
};

struct CrInjectionLedger {
  double injected_energy = 0.0;
  double injected_number = 0.0;
};

// Adds the normalized spectrum to `cell`; throws if [p_min, p_max] is not
// inside the momentum grid or the energy is negative.
void injectCosmicRays(CrSpectralState& state, std::size_t cell, double energy, const CrInjectionSpectrum& spectrum,
                      CrInjectionLedger& ledger);

// Synchrotron + inverse-Compton coefficient C [1/s] for electrons/positrons.
[[nodiscard]] double electronRadiativeLossCoefficient(double magnetic_field_tesla, double radiation_energy_density_j_m3);

struct CrLossTable {
  std::string sha256;
  std::string reference;
  std::map<std::string, double> coefficients;  // name -> value (units in the name)
};

// Schema chui_cr_loss_rates_v1: '# schema:' and '# reference:' lines, then
// rows '<name> <value>' with names
//   proton_coulomb_cm3_s, proton_hadronic_cm3_s, hadronic_thermalized_fraction.
[[nodiscard]] CrLossTable parseCrLossTable(const std::string& bytes, const std::string& origin);
[[nodiscard]] CrLossTable loadCrLossTable(const std::string& path, const std::string& expected_sha256);

class CrLossProvider {
 public:
  CrLossProvider() = default;  // no table: every proton loss is UNAVAILABLE
  explicit CrLossProvider(CrLossTable table);
  [[nodiscard]] spectral::ProviderValue protonCoulombCoefficient() const;   // cm^3 / s
  [[nodiscard]] spectral::ProviderValue protonHadronicCoefficient() const;  // cm^3 / s
  [[nodiscard]] spectral::ProviderValue hadronicThermalizedFraction() const;
  [[nodiscard]] spectral::ProviderValue secondaryElectronSource() const;   // unavailable in pass 01
  [[nodiscard]] spectral::ProviderValue gammaRayEmissivity() const;        // unavailable in pass 01

 private:
  [[nodiscard]] spectral::ProviderValue lookup(const char* name) const;
  bool m_has_table = false;
  CrLossTable m_table;
};

struct CrGreyLossLedger {
  double cr_energy_lost = 0.0;
  double gas_heating = 0.0;
  double escaping = 0.0;
};

// Exact exponential grey proton losses over dt for electron density n_e [cm^-3];
// throws if the provider is unavailable.
void applyGreyProtonLosses(std::vector<double>& e, std::span<const double> electron_density_cm3, double dt_s,
                           const CrLossProvider& provider, std::vector<double>& gas_heating, double cell_volume,
                           CrGreyLossLedger& ledger);

}  // namespace cosmosim::cosmic_rays
