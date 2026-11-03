#include "cosmosim/physics/thermochemistry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cosmosim::physics {
namespace {

using S = ChemicalSpecies;
constexpr std::size_t idx(S s) { return static_cast<std::size_t>(s); }

constexpr std::size_t k_max_unknowns = k_chemical_species_count + 1U;
using Vec = std::array<double, k_max_unknowns>;
using Mat = std::array<double, k_max_unknowns * k_max_unknowns>;

constexpr double k_newton_g_cgs = 6.67430e-8;
constexpr double k_ev = 1.602176634e-12;
constexpr double k_reference_mu_ionized = 0.59;  // analytic reference only
// N_H / A_V = 1.87e21 cm^-2 mag^-1 (Bohlin, Savage & Drake 1978 with R_V = 3.1).
constexpr double k_nh_per_av = 1.87e21;

// Dense LU with partial pivoting; solves a x = b in place. False if singular.
bool luSolve(Mat a, Vec& b, std::size_t m) {
  for (std::size_t k = 0; k < m; ++k) {
    std::size_t pivot = k;
    double best = std::abs(a[k * k_max_unknowns + k]);
    for (std::size_t i = k + 1; i < m; ++i) {
      const double v = std::abs(a[i * k_max_unknowns + k]);
      if (v > best) { best = v; pivot = i; }
    }
    if (!(best > 0.0) || !std::isfinite(best)) return false;
    if (pivot != k) {
      for (std::size_t j = 0; j < m; ++j) std::swap(a[k * k_max_unknowns + j], a[pivot * k_max_unknowns + j]);
      std::swap(b[k], b[pivot]);
    }
    const double inv = 1.0 / a[k * k_max_unknowns + k];
    for (std::size_t i = k + 1; i < m; ++i) {
      const double f = a[i * k_max_unknowns + k] * inv;
      if (f == 0.0) continue;
      for (std::size_t j = k; j < m; ++j) a[i * k_max_unknowns + j] -= f * a[k * k_max_unknowns + j];
      b[i] -= f * b[k];
    }
  }
  for (std::size_t k = m; k-- > 0;) {
    double s = b[k];
    for (std::size_t j = k + 1; j < m; ++j) s -= a[k * k_max_unknowns + j] * b[j];
    b[k] = s / a[k * k_max_unknowns + k];
    if (!std::isfinite(b[k])) return false;
  }
  return true;
}

[[nodiscard]] double smoothstep(double x) {
  x = std::clamp(x, 0.0, 1.0);
  return x * x * (3.0 - 2.0 * x);
}

// Chemical potential energy of a state relative to neutral atoms [erg cm^-3].
[[nodiscard]] double bindingEnergyDensity(const NetworkState& n, std::size_t species_count) {
  double e = 0.0;
  e += 13.598 * k_ev * n[idx(S::kHII)];
  if (species_count > idx(S::kHm)) e -= 0.754 * k_ev * n[idx(S::kHm)];
  if (species_count > idx(S::kH2)) e -= 4.478 * k_ev * n[idx(S::kH2)];
  if (species_count > idx(S::kH2p)) e += (13.598 - 2.651) * k_ev * n[idx(S::kH2p)];
  e += 24.587 * k_ev * n[idx(S::kHeII)];
  e += (24.587 + 54.418) * k_ev * n[idx(S::kHeIII)];
  if (species_count > idx(S::kCII)) e += 11.260 * k_ev * n[idx(S::kCII)];
  if (species_count > idx(S::kCO)) e -= 11.09 * k_ev * n[idx(S::kCO)];
  return e;
}

}  // namespace

ThermochemistrySolver::ThermochemistrySolver(ThermochemistryConfig config) : m_config(config) {
  if (!(m_config.relative_tolerance > 0.0) || m_config.relative_tolerance >= 0.1 ||
      !(m_config.abundance_absolute_tolerance > 0.0) || m_config.max_substeps == 0U ||
      m_config.max_newton_iterations == 0U || !(m_config.temperature_floor_k > 0.0) ||
      !(m_config.adiabatic_index > 1.0)) {
    throw std::invalid_argument("ThermochemistrySolver: invalid numerical configuration");
  }
  switch (m_config.backend) {
    case ThermochemistryBackend::kFixedHybrid:
      m_manifest = &hybridLightNetworkManifest();
      break;
    case ThermochemistryBackend::kFixedDetailed:
      m_manifest = &detailedCarbonOxygenNetworkManifest();
      break;
    case ThermochemistryBackend::kAnalyticReference:
    case ThermochemistryBackend::kEquilibriumPrimordial:
      break;
    case ThermochemistryBackend::kLegacyCooling:
      throw std::invalid_argument(
          "ThermochemistrySolver does not implement legacy_cooling (owned by the in-solver CoolingHeatingSource)");
  }
  if (m_manifest != nullptr) {
    const NetworkValidationReport report = validateNetworkManifest(*m_manifest);
    if (!report.ok) {
      throw std::logic_error("network manifest '" + m_manifest->name + "' failed validation: " + report.diagnostic);
    }
    m_species_count = m_manifest->species_count;
  }
  // Part 2.5 policy: provenance verification alone never lifts the gate.
  // The equilibrium backend's rate set is primary-verified (see
  // equilibriumPrimordialManifest) but its numerical validation is
  // USER-RESERVED and not yet recorded, so it stays behind the opt-in like
  // the network backends (whose manifests are not qualified(): R13/R15/R17/
  // R20/R21 secondary-verified, R18 and H2 photodissociation heating
  // unverified). Only the analytic reference is exempt.
  const bool uses_unqualified_rates =
      m_config.backend == ThermochemistryBackend::kEquilibriumPrimordial ||
      (m_manifest != nullptr && !m_manifest->qualified());
  if (uses_unqualified_rates && !m_config.allow_unqualified_rates) {
    throw SourceModelUnavailableError(
        "thermochemistry backend '" + std::string(thermochemistryBackendLabel(m_config.backend)) +
        "' uses rate coefficients that are UNQUALIFIED (not yet verified against their publications; see "
        "docs/chemistry_rate_provenance.md). Set physics.thermochemistry_allow_unqualified_rates = true to run it "
        "as an explicitly unqualified research model");
  }
  const bool needs_table = m_config.metal_mode == MetalCoolingMode::kElementTable ||
      m_config.metal_mode == MetalCoolingMode::kElementTableAndFineStructure;
  if (needs_table && m_config.metal_table == nullptr) {
    throw SourceModelUnavailableError(
        "metal cooling mode requires an element-resolved metal cooling table (physics.metal_cooling_table_path); "
        "no production table is bundled");
  }
}

double ThermochemistrySolver::meanMolecularWeight(double mass_density_cgs, double n_particles) noexcept {
  return n_particles > 0.0 ? mass_density_cgs / (k_amu_g * n_particles) : 0.0;
}

double ThermochemistrySolver::particleDensity(const ThermochemistryCellInput& input, const NetworkState& n,
                                              double n_e) const {
  double particles = n_e;
  const std::size_t stored = std::max<std::size_t>(m_species_count, k_hybrid_species_count);
  for (std::size_t s = 0; s < stored; ++s) particles += n[s];
  for (std::size_t e = 2; e < k_chemical_element_count; ++e) {
    particles += input.gas_element_mass_fraction[e] * input.mass_density_cgs /
        (chemicalElementTable()[e].atomic_mass_amu * k_amu_g);
  }
  if (m_species_count >= k_detailed_species_count) {
    // C+ and CO are stored; the generic count above already counted C and O
    // nuclei as one particle each. CO merges one C and one O particle.
    particles -= n[idx(S::kCO)];
  }
  return particles;
}

RateEnvironment ThermochemistrySolver::environment(const ThermochemistryCellInput& input, const NetworkState& n,
                                                   double temperature_k) const {
  RateEnvironment env;
  const double rho = input.mass_density_cgs;
  const auto nuclei = [&](ChemicalElement e) {
    return input.gas_element_mass_fraction[static_cast<std::size_t>(e)] * rho /
        (chemicalElementInfo(e).atomic_mass_amu * k_amu_g);
  };
  env.temperature_k = std::max(temperature_k, 1.0);
  env.hydrogen_nuclei_cm3 = nuclei(ChemicalElement::kH);
  env.helium_nuclei_cm3 = nuclei(ChemicalElement::kHe);
  env.oxygen_gas_nuclei_cm3 = nuclei(ChemicalElement::kO);
  env.carbon_gas_nuclei_cm3 = nuclei(ChemicalElement::kC);
  env.redshift = input.redshift;
  env.dust_to_gas_relative_mw = input.dust_to_gas_relative_mw;
  env.g0_habing = m_config.isrf_g0_habing;
  // Local shielding length: Jeans length capped by shielding_length_max_cm
  // (Schaye 2001-type local approximation; NOT radiative transfer).
  const double n_e = electronDensity(n, std::max<std::size_t>(m_species_count, k_hybrid_species_count));
  const double n_part = particleDensity(input, n, n_e);
  const double mu = meanMolecularWeight(rho, n_part);
  double length = m_config.shielding_length_max_cm;
  if (mu > 0.0 && rho > 0.0) {
    const double cs2 = m_config.adiabatic_index * k_boltzmann_erg_per_k * env.temperature_k / (mu * k_amu_g);
    length = std::min(length, std::sqrt(3.14159265358979323846 * cs2 / (k_newton_g_cgs * rho)));
  }
  env.h2_column_cm2 = m_species_count > idx(S::kH2) ? n[idx(S::kH2)] * length : 0.0;
  env.visual_extinction_mag = env.hydrogen_nuclei_cm3 * length / k_nh_per_av * env.dust_to_gas_relative_mw;
  if (m_config.uv_background != nullptr) {
    const UvBackgroundRates uvb = m_config.uv_background->ratesAtRedshift(input.redshift);
    const double shield = m_config.self_shielding
        ? rahmati2013ShieldingFactor(env.hydrogen_nuclei_cm3, env.temperature_k, uvb.gamma_hi_s)
        : 1.0;
    env.gamma_hi_s = uvb.gamma_hi_s * shield;
    env.gamma_hei_s = uvb.gamma_hei_s * shield;
    env.gamma_heii_s = uvb.gamma_heii_s * shield;
    env.heat_hi_erg_s = uvb.heat_hi_erg_s * shield;
    env.heat_hei_erg_s = uvb.heat_hei_erg_s * shield;
    env.heat_heii_erg_s = uvb.heat_heii_erg_s * shield;
  }
  return env;
}

NetworkState ThermochemistrySolver::equilibriumPrimordialState(double temperature_k, double n_h, double n_he,
                                                               const RateEnvironment& env) const {
  NetworkState n{};
  const double t = std::max(temperature_k, 1.0);
  const double a_hp = rates::recombinationHII(t);
  const double a_hep = rates::recombinationHeII(t) + rates::dielectronicHeII(t);
  const double a_hepp = rates::recombinationHeIII(t);
  const double g_h0 = rates::collisionalIonizationHI(t);
  const double g_he0 = rates::collisionalIonizationHeI(t);
  const double g_hep = rates::collisionalIonizationHeII(t);
  double n_e = n_h;
  for (int iteration = 0; iteration < 200; ++iteration) {
    const double ne = std::max(n_e, 1.0e-30 * std::max(n_h, 1.0e-300));
    const double h0 = n_h * a_hp / (a_hp + g_h0 + env.gamma_hi_s / ne);
    const double hp = n_h - h0;
    const double ion_he0 = g_he0 + env.gamma_hei_s / ne;
    const double ion_hep = g_hep + env.gamma_heii_s / ne;
    double hep = 0.0, he0 = 0.0, hepp = 0.0;
    if (ion_he0 > 0.0) {
      hep = n_he / (1.0 + a_hep / ion_he0 + ion_hep / a_hepp);
      he0 = hep * a_hep / ion_he0;
      hepp = hep * ion_hep / a_hepp;
    } else {
      he0 = n_he;
    }
    const double ne_new = hp + hep + 2.0 * hepp;
    n[idx(S::kHI)] = h0;
    n[idx(S::kHII)] = hp;
    n[idx(S::kHeI)] = he0;
    n[idx(S::kHeII)] = hep;
    n[idx(S::kHeIII)] = hepp;
    if (std::abs(ne_new - n_e) <= 1.0e-12 * std::max(n_e, 1.0e-300)) break;
    n_e = 0.5 * (n_e + ne_new);
  }
  return n;
}

CoolingHeatingTerms ThermochemistrySolver::coolingHeating(const ThermochemistryCellInput& input,
                                                          const NetworkState& n, double temperature_k) const {
  const std::size_t stored = std::max<std::size_t>(m_species_count, k_hybrid_species_count);
  const RateEnvironment env = environment(input, n, temperature_k);
  const double n_e = electronDensity(n, stored);
  CoolingHeatingTerms terms = cooling::primordial(env, n, stored, n_e);
  if (!m_config.cosmological) {
    // Compton exchange with the CMB is a cosmological-frame term only.
    const double t_gamma = 2.7255 * (1.0 + std::max(env.redshift, 0.0));
    const double compton = cooling::comptonCoolingCoefficient(env.redshift) * (env.temperature_k - t_gamma) * n_e;
    if (compton >= 0.0) terms.radiative_cooling -= compton; else terms.heating += compton;
  }
  double metal = 0.0;
  const double log_t = std::log10(env.temperature_k);
  const double w_table = smoothstep((log_t - (std::log10(m_config.fine_structure_switch_temperature_k) - 0.15)) / 0.3);
  const bool use_fs = m_config.metal_mode == MetalCoolingMode::kFineStructureCold ||
      m_config.metal_mode == MetalCoolingMode::kElementTableAndFineStructure;
  const bool use_table = m_config.metal_mode == MetalCoolingMode::kElementTable ||
      m_config.metal_mode == MetalCoolingMode::kElementTableAndFineStructure;
  if (use_fs && w_table < 1.0) {
    double n_cii = env.carbon_gas_nuclei_cm3;
    double n_oi = env.oxygen_gas_nuclei_cm3;
    if (m_species_count >= k_detailed_species_count) {
      n_cii = n[idx(S::kCII)];
      n_oi = std::max(env.oxygen_gas_nuclei_cm3 - n[idx(S::kCO)], 0.0);
    } else {
      // Hybrid: C assumed singly ionized (NL97 convention); O ionization
      // locked to H by charge exchange (x_O+ ~ x_H+).
      const double neutral_h = env.hydrogen_nuclei_cm3 > 0.0
          ? (n[idx(S::kHI)] + 2.0 * n[idx(S::kH2)]) / env.hydrogen_nuclei_cm3
          : 0.0;
      n_oi *= std::clamp(neutral_h, 0.0, 1.0);
    }
    const double n_h_collider = n[idx(S::kHI)] + n[idx(S::kH2)];
    const double fs = cooling::fineStructureCII(env.temperature_k, n_cii, n_h_collider, n_e) +
        cooling::fineStructureOI(env.temperature_k, n_oi, n_h_collider);
    metal += (1.0 - w_table) * fs;
  }
  if (use_table && (m_config.metal_mode == MetalCoolingMode::kElementTable || w_table > 0.0)) {
    MetalCoolingTable::Query query;
    query.temperature_k = env.temperature_k;
    query.hydrogen_number_density_cm3 = env.hydrogen_nuclei_cm3;
    query.redshift = env.redshift;
    query.electron_per_hydrogen = env.hydrogen_nuclei_cm3 > 0.0 ? n_e / env.hydrogen_nuclei_cm3 : -1.0;
    query.gas_mass_fraction = input.gas_element_mass_fraction;
    query.hydrogen_mass_fraction = input.gas_element_mass_fraction[0];
    const double weight = m_config.metal_mode == MetalCoolingMode::kElementTable ? 1.0 : w_table;
    const MetalCoolingTable::Result table_result = m_config.metal_table->evaluate(query, m_config.metal_table_range);
    metal += weight * table_result.net_cooling_erg_cm3_s;
    terms.metal_table_clamped = table_result.clamped;
  }
  if (metal >= 0.0) {
    terms.radiative_cooling += metal;
    terms.metal_cooling = metal;
  } else {
    terms.heating -= metal;  // net metal photoheating (table can be signed)
    terms.metal_cooling = metal;
  }
  terms.heating += cooling::photoelectricHeating(env.temperature_k, env.hydrogen_nuclei_cm3, n_e,
                                                 env.g0_habing * std::exp(-1.8 * env.visual_extinction_mag),
                                                 env.dust_to_gas_relative_mw);
  return terms;
}

void ThermochemistrySolver::advanceCell(const ThermochemistryCellInput& input, double dt_seconds,
                                        ThermochemistryCellOutput& output, SourceLedger& ledger) const {
  output = ThermochemistryCellOutput{};
  output.suggested_retry_seconds = 0.0;
  if (!(dt_seconds > 0.0) || !std::isfinite(dt_seconds) || !(input.mass_density_cgs > 0.0) ||
      !(input.specific_internal_energy_erg_g > 0.0) || !std::isfinite(input.specific_internal_energy_erg_g)) {
    output.status = SourceAdvanceStatus::kRejectedInvalid;
    output.diagnostic = "thermochemistry: invalid cell input (non-positive density/energy or interval)";
    return;
  }
  const double rho = input.mass_density_cgs;
  const double gm1 = m_config.adiabatic_index - 1.0;
  const double n_h = input.gas_element_mass_fraction[0] * rho / (chemicalElementInfo(ChemicalElement::kH).atomic_mass_amu * k_amu_g);
  const double n_he = input.gas_element_mass_fraction[1] * rho / (chemicalElementInfo(ChemicalElement::kHe).atomic_mass_amu * k_amu_g);
  const double n_c = input.gas_element_mass_fraction[2] * rho / (chemicalElementInfo(ChemicalElement::kC).atomic_mass_amu * k_amu_g);
  if (!(n_h > 0.0)) {
    output.status = SourceAdvanceStatus::kRejectedInvalid;
    output.diagnostic = "thermochemistry: cell has no hydrogen";
    return;
  }
  const double u_scale = input.specific_internal_energy_erg_g;
  const bool network = m_species_count > 0U;
  const std::size_t m = network ? m_species_count + 1U : 1U;
  const std::size_t w_index = m - 1U;

  // --- state <-> unknowns ---------------------------------------------------
  const auto temperature_of = [&](const NetworkState& n, double u) {
    const double ne = electronDensity(n, std::max<std::size_t>(m_species_count, k_hybrid_species_count));
    return gm1 * u * rho / (k_boltzmann_erg_per_k * particleDensity(input, n, ne));
  };
  // Equilibrium backend: T(u) solved with the equilibrium ionization state.
  const auto equilibrium_temperature = [&](double u, NetworkState& n_out) {
    double lo = 1.0, hi = 1.0e10;
    NetworkState n{};
    for (int it = 0; it < 80; ++it) {
      const double t = std::sqrt(lo * hi);
      const RateEnvironment env = environment(input, n, t);
      n = equilibriumPrimordialState(t, n_h, n_he, env);
      const double t_of_u = temperature_of(n, u);
      if (t_of_u > t) lo = t; else hi = t;
      if (hi / lo < 1.0 + 1.0e-10) break;
    }
    const double t = std::sqrt(lo * hi);
    n_out = equilibriumPrimordialState(t, n_h, n_he, environment(input, n, t));
    return t;
  };

  const double u_floor_reference = [&] {
    // u corresponding to the temperature floor for a neutral atomic state
    // (used only to bound the analytic/equilibrium search).
    return m_config.temperature_floor_k * k_boltzmann_erg_per_k / (gm1 * 1.22 * k_amu_g);
  }();

  bool rhs_failed = false;
  const auto rhs = [&](const Vec& x, Vec& f) {
    f.fill(0.0);
    const double u = x[w_index] * u_scale;
    if (!(u > 0.0) || !std::isfinite(u)) { rhs_failed = true; return; }
    if (m_config.backend == ThermochemistryBackend::kAnalyticReference) {
      const double u_eq = k_boltzmann_erg_per_k * m_config.analytic_equilibrium_temperature_k /
          (gm1 * k_reference_mu_ionized * k_amu_g);
      f[0] = (-(u - u_eq) / m_config.analytic_relaxation_time_s + m_config.analytic_constant_heating_erg_g_s) / u_scale;
      return;
    }
    NetworkState n{};
    double t = 0.0;
    if (network) {
      for (std::size_t s = 0; s < m_species_count; ++s) n[s] = std::max(x[s], 0.0) * n_h;
      t = temperature_of(n, u);
    } else {
      t = equilibrium_temperature(u, n);
    }
    const CoolingHeatingTerms terms = coolingHeating(input, n, t);
    double du = terms.net();
    if (!network) du = terms.heating - terms.radiative_cooling - terms.chemical_cooling;
    f[w_index] = du / (rho * u_scale);
    if (network) {
      const RateEnvironment env = environment(input, n, t);
      const double ne = electronDensity(n, m_species_count);
      for (const ReactionDefinition& reaction : m_manifest->reactions) {
        const double r = reaction.rate(env, n, ne);
        if (!std::isfinite(r)) { rhs_failed = true; return; }
        if (r == 0.0) continue;
        for (std::size_t s = 0; s < m_species_count; ++s) {
          if (reaction.net_change[s] != 0) f[s] += reaction.net_change[s] * r / n_h;
        }
      }
    }
    for (std::size_t i = 0; i < m; ++i) {
      if (!std::isfinite(f[i])) { rhs_failed = true; return; }
    }
  };

  Vec atol{};
  for (std::size_t i = 0; i < m; ++i) atol[i] = m_config.abundance_absolute_tolerance;
  atol[w_index] = 1.0e-3 * m_config.relative_tolerance;  // energy is normalized by u0

  // Backward Euler step with Newton iteration and positivity damping.
  std::uint64_t newton_iterations = 0;
  double max_newton_residual = 0.0;
  const auto be_step = [&](const Vec& x0, double h, Vec& x) -> bool {
    x = x0;
    Vec f{}, f_pert{};
    for (std::uint32_t it = 0; it < m_config.max_newton_iterations; ++it) {
      ++newton_iterations;
      rhs_failed = false;
      rhs(x, f);
      if (rhs_failed) return false;
      Vec residual{};
      for (std::size_t i = 0; i < m; ++i) residual[i] = -(x[i] - x0[i] - h * f[i]);
      Mat jac{};
      for (std::size_t j = 0; j < m; ++j) {
        Vec xp = x;
        const double typical = (j == w_index) ? 1.0 : std::max(1.0e-6, atol[j] * 1.0e4);
        const double delta = 1.0e-7 * std::max(std::abs(x[j]), typical);
        xp[j] += delta;
        rhs_failed = false;
        rhs(xp, f_pert);
        if (rhs_failed) return false;
        for (std::size_t i = 0; i < m; ++i) {
          jac[i * k_max_unknowns + j] = (i == j ? 1.0 : 0.0) - h * (f_pert[i] - f[i]) / delta;
        }
      }
      Vec step = residual;
      if (!luSolve(jac, step, m)) return false;
      double lambda = 1.0;
      for (std::size_t i = 0; i < m; ++i) {
        if (x[i] + step[i] < 0.0 && step[i] < 0.0) {
          lambda = std::min(lambda, x[i] > 0.0 ? 0.9 * x[i] / -step[i] : 0.0);
        }
      }
      if (lambda <= 0.0) lambda = 1.0;  // species already at zero: allow, positivity checked after
      double norm = 0.0;
      for (std::size_t i = 0; i < m; ++i) {
        x[i] += lambda * step[i];
        norm = std::max(norm, std::abs(lambda * step[i]) /
                                  (atol[i] + m_config.relative_tolerance * std::abs(x[i])));
      }
      max_newton_residual = std::max(max_newton_residual, norm);
      if (norm < 1.0e-2) return true;
    }
    return false;
  };

  // --- initial unknowns ---------------------------------------------------------
  Vec x{};
  NetworkState n_begin{};
  if (network) {
    for (std::size_t s = 0; s < m_species_count; ++s) {
      n_begin[s] = std::max(input.species_cm3[s], 0.0);
      x[s] = n_begin[s] / n_h;
    }
  }
  x[w_index] = 1.0;
  double t_begin = 0.0;
  if (network) {
    t_begin = temperature_of(n_begin, u_scale);
  } else if (m_config.backend == ThermochemistryBackend::kEquilibriumPrimordial) {
    t_begin = equilibrium_temperature(u_scale, n_begin);
  } else {
    t_begin = gm1 * k_reference_mu_ionized * k_amu_g * u_scale / k_boltzmann_erg_per_k;
  }
  const double binding_begin = network || m_config.backend == ThermochemistryBackend::kEquilibriumPrimordial
      ? bindingEnergyDensity(n_begin, std::max<std::size_t>(m_species_count, k_hybrid_species_count))
      : 0.0;

  // Initial step from the explicit rate scale.
  Vec f0{};
  rhs_failed = false;
  rhs(x, f0);
  if (rhs_failed) {
    output.status = SourceAdvanceStatus::kRejectedInvalid;
    output.diagnostic = "thermochemistry: non-finite rates at the interval start";
    return;
  }
  double rate_scale = 0.0;
  for (std::size_t i = 0; i < m; ++i) rate_scale = std::max(rate_scale, std::abs(f0[i]) / (std::abs(x[i]) + atol[i] * 1.0e3));
  double h = rate_scale > 0.0 ? std::min(dt_seconds, 0.1 / rate_scale) : dt_seconds;
  const double h_min = m_config.minimum_step_fraction * dt_seconds;

  long double t_done = 0.0L;
  const long double t_request = static_cast<long double>(dt_seconds);
  std::uint32_t attempts = 0;
  double radiative = 0.0, heating = 0.0, chemical = 0.0, floor_injection = 0.0;
  double max_error = 0.0;
  std::uint64_t accepted = 0, rejected = 0, floor_count = 0, metal_clamp_count = 0;

  while (!sourceIntervalClosed(t_done, t_request)) {
    if (++attempts > m_config.max_substeps) {
      output.status = SourceAdvanceStatus::kRetryShorterInterval;
      output.completed_seconds = static_cast<double>(t_done);
      output.suggested_retry_seconds = std::max(static_cast<double>(t_done), h_min);
      output.diagnostic = "thermochemistry: max_substeps exhausted before the accepted interval was integrated "
                          "(completed " + std::to_string(static_cast<double>(t_done)) + " of " +
                          std::to_string(dt_seconds) + " s); stage rejected, nothing committed (audit F2)";
      ledger.rejected_substeps += rejected;
      return;
    }
    const double remaining = static_cast<double>(t_request - t_done);
    bool final_step = false;
    if (h >= remaining * (1.0 - 1.0e-12)) { h = remaining; final_step = true; }
    Vec x_full{}, x_half{}, x_half2{};
    const bool ok = be_step(x, h, x_full) && be_step(x, 0.5 * h, x_half) && be_step(x_half, 0.5 * h, x_half2);
    bool positive = ok;
    if (ok) {
      for (std::size_t i = 0; i < m; ++i) {
        if (x_half2[i] < -atol[i] || x_full[i] < -atol[i]) positive = false;
      }
      if (!(x_half2[w_index] > 0.0)) positive = false;
    }
    if (!positive) {
      ++rejected;
      h *= 0.25;
      if (h < h_min) {
        output.status = SourceAdvanceStatus::kRetryShorterInterval;
        output.completed_seconds = static_cast<double>(t_done);
        output.suggested_retry_seconds = std::max(static_cast<double>(t_done), h_min);
        output.diagnostic = "thermochemistry: Newton/positivity failure below the minimum step";
        ledger.rejected_substeps += rejected;
        return;
      }
      continue;
    }
    double err = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
      const double scale = atol[i] + m_config.relative_tolerance * std::max(std::abs(x_half2[i]), std::abs(x[i]));
      err = std::max(err, std::abs(x_full[i] - x_half2[i]) / scale);
    }
    if (err > 1.0) {
      ++rejected;
      h *= std::max(0.2, 0.9 / std::sqrt(err));
      if (h < h_min) {
        output.status = SourceAdvanceStatus::kRetryShorterInterval;
        output.completed_seconds = static_cast<double>(t_done);
        output.suggested_retry_seconds = std::max(static_cast<double>(t_done), h_min);
        output.diagnostic = "thermochemistry: error control demands a step below the minimum";
        ledger.rejected_substeps += rejected;
        return;
      }
      continue;
    }
    max_error = std::max(max_error, err);
    // Richardson extrapolation (second order) when it stays admissible.
    Vec x_new{};
    bool richardson_ok = true;
    for (std::size_t i = 0; i < m; ++i) {
      x_new[i] = 2.0 * x_half2[i] - x_full[i];
      if (x_new[i] < 0.0) richardson_ok = false;
    }
    if (!richardson_ok || !(x_new[w_index] > 0.0)) x_new = x_half2;
    for (std::size_t i = 0; i < m; ++i) x_new[i] = std::max(x_new[i], 0.0);

    // Ledger integrals (backward-Euler consistent: rates at the new state).
    {
      NetworkState n{};
      double t = 0.0;
      const double u = x_new[w_index] * u_scale;
      if (network) {
        for (std::size_t s = 0; s < m_species_count; ++s) n[s] = x_new[s] * n_h;
        t = temperature_of(n, u);
      } else if (m_config.backend == ThermochemistryBackend::kEquilibriumPrimordial) {
        t = equilibrium_temperature(u, n);
      }
      if (m_config.backend != ThermochemistryBackend::kAnalyticReference) {
        const CoolingHeatingTerms terms = coolingHeating(input, n, t);
        if (terms.metal_table_clamped) ++metal_clamp_count;  // once per accepted substep
        radiative += h * terms.radiative_cooling;
        heating += h * terms.heating;
        chemical += h * terms.chemical_cooling;
        // Temperature floor: physical lower bound of the model, applied only
        // after a successful step, ledgered as injected energy.
        if (t < m_config.temperature_floor_k && t > 0.0) {
          const double u_floor = u * m_config.temperature_floor_k / t;
          floor_injection += rho * (u_floor - u);
          x_new[w_index] = u_floor / u_scale;
          ++floor_count;
        }
      } else if (x_new[w_index] * u_scale < u_floor_reference) {
        floor_injection += rho * (u_floor_reference - x_new[w_index] * u_scale);
        x_new[w_index] = u_floor_reference / u_scale;
        ++floor_count;
      }
    }
    x = x_new;
    ++accepted;
    t_done = final_step ? t_request : t_done + static_cast<long double>(h);
    h *= std::min(4.0, 0.9 / std::sqrt(std::max(err, 1.0e-10)));
  }

  // --- projection onto exact element totals (round-off only) ----------------
  NetworkState n_end{};
  const double u_end = x[w_index] * u_scale;
  double projection = 0.0;
  if (network) {
    for (std::size_t s = 0; s < m_species_count; ++s) n_end[s] = x[s] * n_h;
    const auto project = [&](std::initializer_list<S> group, double target, auto nuclei_of) {
      double sum = 0.0;
      for (const S s : group) sum += nuclei_of(s) * n_end[idx(s)];
      if (sum > 0.0 && target > 0.0) {
        projection = std::max(projection, std::abs(sum - target) / target);
        const double f = target / sum;
        for (const S s : group) n_end[idx(s)] *= f;
      }
    };
    project({S::kHI, S::kHII, S::kHm, S::kH2, S::kH2p}, n_h,
            [](S s) { return static_cast<double>(chemicalSpeciesInfo(s).n_h); });
    project({S::kHeI, S::kHeII, S::kHeIII}, n_he, [](S) { return 1.0; });
    if (m_species_count >= k_detailed_species_count) {
      project({S::kCII, S::kCO}, n_c, [](S) { return 1.0; });
    }
    output.temperature_k = temperature_of(n_end, u_end);
  } else if (m_config.backend == ThermochemistryBackend::kEquilibriumPrimordial) {
    output.temperature_k = equilibrium_temperature(u_end, n_end);
  } else {
    output.temperature_k = gm1 * k_reference_mu_ionized * k_amu_g * u_end / k_boltzmann_erg_per_k;
  }
  const std::size_t stored = std::max<std::size_t>(m_species_count, k_hybrid_species_count);
  const double ne_end = electronDensity(n_end, stored);
  output.status = SourceAdvanceStatus::kCompleted;
  output.completed_seconds = dt_seconds;
  output.specific_internal_energy_erg_g = u_end;
  output.species_cm3 = n_end;
  output.electron_per_hydrogen = ne_end / n_h;
  output.mean_molecular_weight = m_config.backend == ThermochemistryBackend::kAnalyticReference
      ? k_reference_mu_ionized
      : meanMolecularWeight(rho, particleDensity(input, n_end, ne_end));

  ledger.radiative_loss_erg += radiative;
  ledger.heating_gain_erg += heating;
  ledger.chemical_binding_delta_erg +=
      (network || m_config.backend == ThermochemistryBackend::kEquilibriumPrimordial)
          ? bindingEnergyDensity(n_end, stored) - binding_begin
          : 0.0;
  ledger.floor_injection_erg += floor_injection;
  ledger.thermal_energy_delta_erg += rho * (u_end - u_scale);
  ledger.element_mass_residual_code = std::max(ledger.element_mass_residual_code, projection);
  double charge = 0.0;
  for (std::size_t s = 0; s < stored; ++s) charge += chemicalSpeciesInfo(static_cast<S>(s)).charge * n_end[s];
  ledger.charge_residual = std::max(ledger.charge_residual, std::abs(charge - ne_end) / n_h);
  ledger.cells += 1U;
  ledger.accepted_substeps += accepted;
  ledger.rejected_substeps += rejected;
  ledger.newton_iterations += newton_iterations;
  ledger.floor_activations += floor_count;
  ledger.metal_table_clamps += metal_clamp_count;
  ledger.max_error_estimate = std::max(ledger.max_error_estimate, max_error);
  ledger.max_newton_residual = std::max(ledger.max_newton_residual, max_newton_residual);
  (void)t_begin;
  (void)chemical;
}

}  // namespace cosmosim::physics
