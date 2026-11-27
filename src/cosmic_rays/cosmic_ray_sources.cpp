// CHUI-EXT-PHYS-REF-v1 CR05/CR06 injection and loss providers (source only,
// NOT compiled or executed in pass 01).

#include "cosmosim/cosmic_rays/cosmic_ray_sources.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "core/internal/sha256.hpp"
#include "cosmosim/core/constants.hpp"

namespace cosmosim::cosmic_rays {
namespace {

namespace k = core::constants;

[[nodiscard]] std::string trim(const std::string& s) {
  std::size_t b = 0;
  std::size_t e = s.size();
  while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
  return s.substr(b, e - b);
}

}  // namespace

void injectCosmicRays(CrSpectralState& state, std::size_t cell, double energy, const CrInjectionSpectrum& spec,
                      CrInjectionLedger& ledger) {
  if (!(energy >= 0.0) || !std::isfinite(energy)) throw std::invalid_argument("CR injection energy must be >= 0");
  if (energy == 0.0) return;
  const CrMomentumGrid& grid = state.grid();
  if (!(spec.p_min >= grid.edges.front()) || !(spec.p_max <= grid.edges.back()) || !(spec.p_max > spec.p_min)) {
    throw std::invalid_argument("CR injection momentum range must lie inside the momentum grid");
  }
  const double rest = crSpeciesRestEnergyJ(state.settings().species);
  const std::size_t nb = grid.bins();
  std::vector<double> raw_n(nb, 0.0);
  std::vector<double> raw_e(nb, 0.0);
  double total_e = 0.0;
  for (std::size_t b = 0; b < nb; ++b) {
    const double lo = std::max(grid.edges[b], spec.p_min);
    const double hi = std::min(grid.edges[b + 1U], spec.p_max);
    if (!(hi > lo)) continue;
    // f = (p / lo)^-alpha * lo^-alpha (unit normalization p^-alpha).
    const double f0 = std::pow(lo, -spec.alpha);
    raw_n[b] = powerLawNumber(lo, hi, spec.alpha, f0);
    raw_e[b] = powerLawEnergy(lo, hi, spec.alpha, f0, rest);
    total_e += raw_e[b];
  }
  const bool power_law = state.settings().representation == CrSpectralRepresentation::kPiecewisePowerLaw;
  if (!power_law) {
    // Fixed bins carry number only and represent f as constant over the whole
    // bin, so the normalization uses the energy that representation carries.
    total_e = 0.0;
    for (std::size_t b = 0; b < nb; ++b) {
      if (raw_n[b] == 0.0) continue;
      const double lo = grid.edges[b];
      const double hi = grid.edges[b + 1U];
      raw_e[b] = powerLawEnergy(lo, hi, 0.0, raw_n[b] / powerLawNumber(lo, hi, 0.0, 1.0), rest);
      total_e += raw_e[b];
    }
  }
  if (!(total_e > 0.0)) throw std::runtime_error("CR injection spectrum carries no energy on the grid");
  const double scale = energy / total_e;
  for (std::size_t b = 0; b < nb; ++b) {
    const double dn = raw_n[b] * scale;
    state.number(b)[cell] += dn;
    ledger.injected_number += dn;
    if (power_law) state.energy(b)[cell] += raw_e[b] * scale;
  }
  ledger.injected_energy += energy;  // sum of represented bin energies equals E_inj by construction
}

double electronRadiativeLossCoefficient(double b_tesla, double u_rad) {
  if (!std::isfinite(b_tesla) || !std::isfinite(u_rad) || u_rad < 0.0) {
    throw std::invalid_argument("radiative loss coefficient needs finite B and U_rad >= 0");
  }
  const double mu0 = 1.0 / (k::k_vacuum_permittivity_si * k::k_speed_of_light_si * k::k_speed_of_light_si);
  const double u_b = b_tesla * b_tesla / (2.0 * mu0);
  return 4.0 / 3.0 * k::k_thomson_cross_section_si * (u_b + u_rad) / (k::k_electron_mass_si * k::k_speed_of_light_si);
}

CrLossTable parseCrLossTable(const std::string& bytes, const std::string& origin) {
  CrLossTable table;
  std::string schema;
  std::istringstream stream(bytes);
  std::string line;
  while (std::getline(stream, line)) {
    const std::string t = trim(line);
    if (t.empty()) continue;
    if (t.front() == '#') {
      const auto colon = t.find(':');
      if (colon == std::string::npos) continue;
      const std::string key = trim(t.substr(1, colon - 1));
      if (key == "schema") schema = trim(t.substr(colon + 1));
      if (key == "reference") table.reference = trim(t.substr(colon + 1));
      continue;
    }
    std::istringstream row(t);
    std::string name;
    double value = 0.0;
    row >> name >> value;
    if (row.fail() || !std::isfinite(value) || value < 0.0) {
      throw std::runtime_error("CR loss table '" + origin + "': malformed row '" + t + "'");
    }
    if (!table.coefficients.emplace(name, value).second) {
      throw std::runtime_error("CR loss table '" + origin + "': duplicated coefficient '" + name + "'");
    }
  }
  if (schema != "chui_cr_loss_rates_v1") throw std::runtime_error("CR loss table '" + origin + "': wrong schema");
  if (table.reference.empty()) {
    throw std::runtime_error("CR loss table '" + origin + "': a '# reference:' line citing the rate source is required");
  }
  const auto fraction = table.coefficients.find("hadronic_thermalized_fraction");
  if (fraction != table.coefficients.end() && fraction->second > 1.0) {
    throw std::runtime_error("CR loss table '" + origin + "': hadronic_thermalized_fraction must lie in [0, 1]");
  }
  table.sha256 = core::internal::sha256Hex(bytes);
  return table;
}

CrLossTable loadCrLossTable(const std::string& path, const std::string& expected_sha256) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open CR loss table '" + path + "'");
  const std::string bytes((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  const std::string sha = core::internal::sha256Hex(bytes);
  if (sha != expected_sha256) throw std::runtime_error("CR loss table '" + path + "' SHA-256 mismatch");
  return parseCrLossTable(bytes, path);
}

CrLossProvider::CrLossProvider(CrLossTable table) : m_has_table(true), m_table(std::move(table)) {}

spectral::ProviderValue CrLossProvider::lookup(const char* name) const {
  if (!m_has_table) return spectral::ProviderValue::unavailable("no chui_cr_loss_rates_v1 table configured");
  const auto it = m_table.coefficients.find(name);
  if (it == m_table.coefficients.end()) {
    return spectral::ProviderValue::unavailable(std::string("coefficient '") + name + "' absent from the loss table");
  }
  return it->second == 0.0 ? spectral::ProviderValue::physicalZero() : spectral::ProviderValue::available(it->second);
}

spectral::ProviderValue CrLossProvider::protonCoulombCoefficient() const { return lookup("proton_coulomb_cm3_s"); }
spectral::ProviderValue CrLossProvider::protonHadronicCoefficient() const { return lookup("proton_hadronic_cm3_s"); }
spectral::ProviderValue CrLossProvider::hadronicThermalizedFraction() const {
  return lookup("hadronic_thermalized_fraction");
}
spectral::ProviderValue CrLossProvider::secondaryElectronSource() const {
  return spectral::ProviderValue::unavailable("secondary e+- production needs pp cross-section tables (not in pass 01)");
}
spectral::ProviderValue CrLossProvider::gammaRayEmissivity() const {
  return spectral::ProviderValue::unavailable("gamma-ray emissivity needs pion-decay tables (not in pass 01)");
}

void applyGreyProtonLosses(std::vector<double>& e, std::span<const double> n_e, double dt, const CrLossProvider& provider,
                           std::vector<double>& gas_heating, double cell_volume, CrGreyLossLedger& ledger) {
  const double coulomb = provider.protonCoulombCoefficient().require("proton Coulomb loss");
  const double hadronic = provider.protonHadronicCoefficient().require("proton hadronic loss");
  const double f_th = provider.hadronicThermalizedFraction().require("hadronic thermalized fraction");
  if (n_e.size() != e.size() || gas_heating.size() != e.size()) throw std::invalid_argument("CR loss array sizes differ");
  const double total = coulomb + hadronic;
  for (std::size_t c = 0; c < e.size(); ++c) {
    if (!(n_e[c] >= 0.0)) throw std::invalid_argument("electron density must be >= 0");
    const double lost = e[c] * -std::expm1(-total * n_e[c] * dt);
    const double to_gas = total > 0.0 ? lost * (coulomb + f_th * hadronic) / total : 0.0;
    e[c] -= lost;
    gas_heating[c] += to_gas;
    ledger.cr_energy_lost += lost * cell_volume;
    ledger.gas_heating += to_gas * cell_volume;
    ledger.escaping += (lost - to_gas) * cell_volume;
  }
}

}  // namespace cosmosim::cosmic_rays
