#include "cosmosim/physics/stellar_population.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/radiation_background.hpp"
#include "cosmosim/physics/source_interval.hpp"

namespace cosmosim::physics {
namespace {

constexpr double k_ln10 = 2.302585092994046;
constexpr double k_chandrasekhar_msun = 1.38;
constexpr double k_neutron_star_msun = 1.4;

// On-disk/in-flight representation written by schema 1. Keep this definition
// local to the migration reader; new code must use StellarPopulationRecord.
struct StellarPopulationRecordV1 {
  std::array<double, k_chemical_element_count> birth_element_mass_fraction{};
  double birth_time_code = 0.0;
  double cc_next_threshold = 0.0;
  double ia_next_threshold = 0.0;
  std::uint32_t cc_realized = 0;
  std::uint32_t ia_realized = 0;
  std::uint64_t event_stream_key = 0;
  double returned_mass_code = 0.0;
  std::uint32_t initialized = 0;
  std::uint32_t reserved = 0;
};

[[nodiscard]] double chabrierXi(double log_m) {
  if (log_m <= 0.0) {
    const double d = log_m - std::log10(0.079);
    return 0.158 * std::exp(-d * d / (2.0 * 0.69 * 0.69));
  }
  return 0.0443 * std::pow(10.0, -1.3 * log_m);
}

template <class F>
double simpson(F f, double a, double b, int n = 64) {
  if (!(b > a)) return 0.0;
  if (n % 2 != 0) ++n;
  const double h = (b - a) / n;
  double sum = f(a) + f(b);
  for (int i = 1; i < n; ++i) sum += (i % 2 == 1 ? 4.0 : 2.0) * f(a + i * h);
  return sum * h / 3.0;
}

[[nodiscard]] std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31U);
}

[[nodiscard]] double metallicityOf(const std::array<double, k_chemical_element_count>& x) {
  double z = 0.0;
  for (std::size_t e = 2; e < k_chemical_element_count; ++e) z += x[e];
  return z;
}

[[nodiscard]] std::string trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r");
  if (b == std::string::npos) return {};
  return s.substr(b, s.find_last_not_of(" \t\r") - b + 1U);
}

}  // namespace

// ---------------------------------------------------------------------------
// Yield table
// ---------------------------------------------------------------------------
SspYieldTable SspYieldTable::loadFromFile(const std::string& path, const std::string& expected_sha256) {
  if (path.empty()) {
    throw SourceModelUnavailableError(
        "stellar_population_model=chabrier_yield_table_v1 needs physics.ssp_yield_table_path: convert a licensed "
        "published yield set (e.g. AGB: Karakas 2010; CC: Nomoto et al. 2013 / Limongi & Chieffi 2018; Ia: "
        "Iwamoto et al. 1999 W7) into chui_ssp_yield_table_v1. No yields are bundled or synthesized");
  }
  std::ifstream input(path, std::ios::binary);
  if (!input) throw SourceModelUnavailableError("cannot open SSP yield table '" + path + "'");
  const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  SspYieldTable table = loadFromString(text, path);
  table.m_sha256 = sha256FileHex(path);
  if (!expected_sha256.empty() && expected_sha256 != table.m_sha256) {
    throw SourceModelUnavailableError("SSP yield table '" + path + "' SHA-256 mismatch");
  }
  return table;
}

SspYieldTable SspYieldTable::loadFromString(const std::string& text, const std::string& label) {
  SspYieldTable table;
  std::istringstream stream(text);
  std::string line;
  bool magic = false;
  std::vector<std::string> elements;
  int channel = -1;
  std::array<std::map<std::pair<std::size_t, std::size_t>, Row>, 3> parsed;
  const auto fail = [&](const std::string& why) {
    throw SourceModelUnavailableError("SSP yield table '" + label + "': " + why);
  };
  while (std::getline(stream, line)) {
    const std::string t = trim(line);
    if (t.empty()) continue;
    if (t[0] == '#') {
      const std::string body = trim(t.substr(1));
      if (body == "chui_ssp_yield_table_v1") magic = true;
      const auto eq = body.find('=');
      if (eq != std::string::npos) {
        const std::string key = trim(body.substr(0, eq)), value = trim(body.substr(eq + 1U));
        if (key == "citation") table.m_citation = value;
        if (key == "license") table.m_license = value;
        if (key == "net_yield_sum_tolerance") {
          try {
            table.net_yield_sum_tolerance = std::stod(value);
          } catch (const std::exception&) {
            fail("net_yield_sum_tolerance is not a number");
          }
          if (!(table.net_yield_sum_tolerance >= 0.0) || table.net_yield_sum_tolerance > 0.05) {
            fail("net_yield_sum_tolerance must be in [0, 0.05]");
          }
        }
      }
      continue;
    }
    std::istringstream row(t);
    std::string key;
    row >> key;
    if (key == "elements") {
      std::string e;
      while (row >> e) elements.push_back(e);
      if (elements.size() != k_chemical_element_count) fail("elements must list the 11 CHUI elements in catalog order");
      for (std::size_t i = 0; i < elements.size(); ++i) {
        if (elements[i] != chemicalElementTable()[i].symbol) fail("element order must match the CHUI catalog");
      }
    } else if (key == "channel") {
      std::string name;
      row >> name;
      channel = name == "agb" ? 0 : name == "cc" ? 1 : name == "ia" ? 2 : -1;
      if (channel < 0) fail("unknown channel '" + name + "'");
    } else if (key == "metallicity" || key == "mass") {
      if (channel < 0) fail("axis before channel");
      std::vector<double> values;
      double v = 0.0;
      while (row >> v) values.push_back(v);
      for (std::size_t i = 1; i < values.size(); ++i) if (!(values[i] > values[i - 1])) fail(key + " axis not increasing");
      if (values.empty()) fail("empty axis");
      (key == "metallicity" ? table.m_channels[channel].metallicity : table.m_channels[channel].mass) = values;
    } else if (key == "row") {
      if (channel < 0) fail("row before channel");
      std::size_t iz = 0, im = 0;
      row >> iz;
      if (channel != 2) row >> im;
      Row r;
      row >> r.ejecta_msun;
      for (double& y : r.net_yield_msun) row >> y;
      if (row.fail()) fail("malformed row");
      parsed[channel][{iz, im}] = r;
    } else {
      fail("unknown key '" + key + "'");
    }
  }
  if (!magic) fail("missing chui_ssp_yield_table_v1 header");
  if (table.m_citation.empty() || table.m_license.empty()) fail("citation and license are mandatory");
  if (elements.empty()) fail("missing elements line");
  for (int c = 0; c < 3; ++c) {
    Channel& ch = table.m_channels[c];
    if (ch.metallicity.empty()) fail("channel " + std::to_string(c) + " missing");
    const std::size_t nm = c == 2 ? 1U : ch.mass.size();
    if (c != 2 && nm == 0U) fail("massive channel without mass axis");
    ch.rows.resize(ch.metallicity.size() * nm);
    for (std::size_t iz = 0; iz < ch.metallicity.size(); ++iz) {
      for (std::size_t im = 0; im < nm; ++im) {
        const auto it = parsed[c].find({iz, c == 2 ? 0U : im});
        if (it == parsed[c].end()) fail("missing row for channel " + std::to_string(c));
        Row r = it->second;
        if (!(r.ejecta_msun >= 0.0)) fail("negative ejecta");
        if (c != 2) {
          // Corrective gate: yields are NEVER modified. The 11 tracked net
          // yields of a mass-conserving set sum to ~0; elements outside the
          // CHUI vector make the tracked sum non-zero. The declared
          // tolerance (header key net_yield_sum_tolerance, relative to the
          // ejecta mass, default 1e-3) bounds that untracked residual; a
          // larger residual rejects the table (inconsistent input) instead
          // of being moved silently into hydrogen.
          double sum = 0.0;
          for (const double p : r.net_yield_msun) sum += p;
          const double residual = std::abs(sum) / std::max(r.ejecta_msun, 1.0e-30);
          table.max_net_yield_residual = std::max(table.max_net_yield_residual, residual);
          if (residual > table.net_yield_sum_tolerance) {
            fail("net yields of channel " + std::to_string(c) + " (Z index " + std::to_string(iz) +
                 ", mass index " + std::to_string(im) + ") sum to " + std::to_string(sum) +
                 " Msun, i.e. " + std::to_string(residual) + " of the ejecta, above net_yield_sum_tolerance=" +
                 std::to_string(table.net_yield_sum_tolerance) + "; the table is refused (yields are never rebalanced)");
          }
        } else {
          double sum = 0.0;
          for (const double p : r.net_yield_msun) {
            if (p < 0.0) fail("Ia ejecta element masses must be non-negative");
            sum += p;
          }
          if (std::abs(sum - r.ejecta_msun) > 1.0e-3 * r.ejecta_msun) fail("Ia element masses must sum to ejecta");
        }
        ch.rows[iz * nm + im] = r;
      }
    }
  }
  return table;
}

SspYieldTable::Row SspYieldTable::massive(std::size_t channel, double mass_msun, double metallicity,
                                          bool& clamped) const {
  const Channel& ch = m_channels.at(channel);
  const auto locate = [&clamped](const std::vector<double>& axis, double x, std::size_t& lo, double& w) {
    // Part 2.5: a singleton axis covers only its one coordinate; any other
    // query uses that row and is REPORTED as clamped (never silent).
    if (axis.size() == 1U) {
      lo = 0;
      w = 0.0;
      clamped = clamped || std::abs(x - axis[0]) > 1.0e-8 * std::max(1.0, std::abs(axis[0]));
      return;
    }
    if (x <= axis.front() || x >= axis.back()) clamped = clamped || x < axis.front() || x > axis.back();
    x = std::clamp(x, axis.front(), axis.back());
    const auto up = std::upper_bound(axis.begin(), axis.end(), x);
    std::size_t hi = std::clamp<std::size_t>(static_cast<std::size_t>(up - axis.begin()), 1U, axis.size() - 1U);
    lo = hi - 1U;
    w = (x - axis[lo]) / (axis[hi] - axis[lo]);
  };
  std::size_t iz = 0, im = 0;
  double wz = 0.0, wm = 0.0;
  locate(ch.metallicity, metallicity, iz, wz);
  locate(ch.mass, mass_msun, im, wm);
  const std::size_t nm = ch.mass.size();
  Row out;
  for (int dz = 0; dz < 2; ++dz) {
    const std::size_t z = std::min(iz + dz, ch.metallicity.size() - 1U);
    const double fz = dz == 0 ? 1.0 - wz : wz;
    for (int dm = 0; dm < 2; ++dm) {
      const std::size_t m = std::min(im + dm, nm - 1U);
      const double f = fz * (dm == 0 ? 1.0 - wm : wm);
      if (f == 0.0) continue;
      const Row& r = ch.rows[z * nm + m];
      out.ejecta_msun += f * r.ejecta_msun;
      for (std::size_t e = 0; e < k_chemical_element_count; ++e) out.net_yield_msun[e] += f * r.net_yield_msun[e];
    }
  }
  return out;
}

SspYieldTable::Row SspYieldTable::snia(double metallicity, bool& clamped) const {
  const Channel& ch = m_channels[2];
  if (ch.metallicity.size() == 1U) {
    // Single-metallicity Ia table: used for every Z, reported as clamped
    // whenever Z differs from the tabulated value (Part 2.5).
    clamped = clamped || std::abs(metallicity - ch.metallicity[0]) > 1.0e-8 * std::max(1.0, std::abs(ch.metallicity[0]));
    return ch.rows[0];
  }
  const double z = std::clamp(metallicity, ch.metallicity.front(), ch.metallicity.back());
  clamped = clamped || z != metallicity;
  const auto up = std::upper_bound(ch.metallicity.begin(), ch.metallicity.end(), z);
  const std::size_t hi = std::clamp<std::size_t>(static_cast<std::size_t>(up - ch.metallicity.begin()), 1U,
                                                 ch.metallicity.size() - 1U);
  const std::size_t lo = hi - 1U;
  const double w = (z - ch.metallicity[lo]) / (ch.metallicity[hi] - ch.metallicity[lo]);
  Row out;
  out.ejecta_msun = (1.0 - w) * ch.rows[lo].ejecta_msun + w * ch.rows[hi].ejecta_msun;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    out.net_yield_msun[e] = (1.0 - w) * ch.rows[lo].net_yield_msun[e] + w * ch.rows[hi].net_yield_msun[e];
  }
  return out;
}

// ---------------------------------------------------------------------------
// SSP model
// ---------------------------------------------------------------------------
SspModel::SspModel(SspConfig config, const SspYieldTable* table) : m_config(config), m_table(table) {
  if (!(m_config.imf_mass_min_msun > 0.0) || !(m_config.imf_mass_max_msun > m_config.imf_mass_min_msun) ||
      !(m_config.cc_mass_min_msun > m_config.imf_mass_min_msun) ||
      !(m_config.cc_mass_max_msun > m_config.cc_mass_min_msun) ||
      m_config.cc_mass_max_msun > m_config.imf_mass_max_msun) {
    throw std::invalid_argument("SspModel: invalid mass limits");
  }
  if (m_config.yield_mode == SspYieldMode::kTable && m_table == nullptr) {
    throw SourceModelUnavailableError("SspModel: table yield mode selected without a loaded yield table");
  }
  const double mass = simpson([](double x) { return std::pow(10.0, x) * chabrierXi(x); },
                              std::log10(m_config.imf_mass_min_msun), std::log10(m_config.imf_mass_max_msun), 512);
  m_normalization = 1.0 / mass;
}

double SspModel::phi(double m) const { return m_normalization * chabrierXi(std::log10(m)) / (m * k_ln10); }

double SspModel::numberBetween(double m_lo, double m_hi) const {
  m_lo = std::max(m_lo, m_config.imf_mass_min_msun);
  m_hi = std::min(m_hi, m_config.imf_mass_max_msun);
  if (!(m_hi > m_lo)) return 0.0;
  return m_normalization * simpson([](double x) { return chabrierXi(x); }, std::log10(m_lo), std::log10(m_hi));
}

double SspModel::massBetween(double m_lo, double m_hi) const {
  m_lo = std::max(m_lo, m_config.imf_mass_min_msun);
  m_hi = std::min(m_hi, m_config.imf_mass_max_msun);
  if (!(m_hi > m_lo)) return 0.0;
  return m_normalization *
      simpson([](double x) { return std::pow(10.0, x) * chabrierXi(x); }, std::log10(m_lo), std::log10(m_hi));
}

double SspModel::lifetimeYears(double mass_msun, double metallicity) const {
  const double lz = std::log10(std::clamp(metallicity, 7.0e-5, 3.0e-2));
  const double a0 = 10.13 + 0.07547 * lz - 0.008084 * lz * lz;
  const double a1 = -4.424 - 0.7939 * lz - 0.1187 * lz * lz;
  const double a2 = 1.262 + 0.3385 * lz + 0.05417 * lz * lz;
  const double lm = std::log10(std::max(mass_msun, 1.0e-3));
  return std::pow(10.0, a0 + a1 * lm + a2 * lm * lm);
}

double SspModel::turnoffMassMsun(double age_years, double metallicity) const {
  if (!(age_years > 0.0)) return m_config.imf_mass_max_msun;
  const double lz = std::log10(std::clamp(metallicity, 7.0e-5, 3.0e-2));
  const double a0 = 10.13 + 0.07547 * lz - 0.008084 * lz * lz;
  const double a1 = -4.424 - 0.7939 * lz - 0.1187 * lz * lz;
  const double a2 = 1.262 + 0.3385 * lz + 0.05417 * lz * lz;
  const double c = a0 - std::log10(age_years);
  const double disc = a1 * a1 - 4.0 * a2 * c;
  if (disc < 0.0) return m_config.imf_mass_max_msun;  // younger than the shortest lifetime
  const double lm = (-a1 - std::sqrt(disc)) / (2.0 * a2);
  return std::clamp(std::pow(10.0, lm), m_config.imf_mass_min_msun, m_config.imf_mass_max_msun);
}

double SspModel::remnantMassMsun(double m) const {
  if (m < m_config.cc_mass_min_msun) return std::min(m, 0.109 * m + 0.394);
  if (m <= m_config.cc_mass_max_msun) return k_neutron_star_msun;
  return m;  // direct collapse: no ejecta (stated assumption)
}

double SspModel::ccCumulativeNumber(double age_years, double metallicity) const {
  const double m_to = turnoffMassMsun(age_years, metallicity);
  return numberBetween(std::max(m_to, m_config.cc_mass_min_msun), m_config.cc_mass_max_msun);
}

double SspModel::iaCumulativeNumber(double age_years) const {
  if (!(age_years > m_config.snia_delay_min_years)) return 0.0;
  const double f = std::log(std::min(age_years, m_config.snia_delay_max_years) / m_config.snia_delay_min_years) /
      std::log(m_config.snia_delay_max_years / m_config.snia_delay_min_years);
  return m_config.snia_per_msun * std::clamp(f, 0.0, 1.0);
}

std::array<double, k_chemical_element_count> SspModel::ejectaElements(
    std::size_t channel, double m, double metallicity, const std::array<double, k_chemical_element_count>& birth_x,
    double& ejecta) const {
  std::array<double, k_chemical_element_count> out{};
  if (m_config.yield_mode == SspYieldMode::kTable) {
    bool clamped = false;
    const SspYieldTable::Row row = m_table->massive(channel, m, metallicity, clamped);
    if (clamped) ++table_clamps;
    // Gross ejecta per element = X_birth,e m_ej + p_e (net-yield convention).
    // No renormalization: the ejected mass is DEFINED as the sum of the
    // tracked element ejecta (element-conserving); its difference from the
    // tabulated ejecta mass (untracked elements, within the table tolerance)
    // is accumulated as a diagnostic. A negative gross ejecta (a net
    // destruction exceeding the birth content, possible when interpolating
    // across metallicity) is clamped to zero and COUNTED.
    double sum = 0.0;
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
      const double gross = birth_x[e] * row.ejecta_msun + row.net_yield_msun[e];
      if (gross < 0.0) ++negative_ejecta_clamps;
      out[e] = std::max(gross, 0.0);
      sum += out[e];
    }
    untracked_mass_residual_msun += row.ejecta_msun - sum;
    ejecta = sum;
    return out;
  }
  ejecta = std::max(m - remnantMassMsun(m), 0.0);
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) out[e] = birth_x[e] * ejecta;
  const double y = m_config.reference_net_metal_yield * ejecta;
  if (y > 0.0) {
    // REFERENCE ONLY: new metals in the solar pattern, taken from H and He.
    const auto& solar = solarElementMassFractions();
    const double z_sun = solarMetalMassFraction();
    const double h0 = out[0], he0 = out[1];
    const double hhe = h0 + he0;
    if (hhe > 0.0) {
      // Subtract proportional H/He shares, add the same mass as metals.
      const double take = std::min(y, hhe);
      out[0] -= take * h0 / hhe;
      out[1] -= take * he0 / hhe;
      for (std::size_t e = 2; e < k_chemical_element_count; ++e) out[e] += take * solar[e] / z_sun;
    }
  }
  return out;
}

SspModel::ChannelReturn SspModel::cumulativeMassiveReturn(
    std::size_t channel, double age, const std::array<double, k_chemical_element_count>& birth_x) const {
  // Everything ejected by the channel up to `age`: stars with masses between
  // the turnoff mass m_to(age) and the channel's upper limit.
  //
  // G3.2 (replaces the age-dependent-node Simpson rule, whose nodes moved
  // with m_to and made C(age) structurally NON-monotone): the integrand
  // f(x) = phi(m) m ln10 * ejecta(m), x = log10 m, is sampled on a FIXED
  // grid x_0..x_K spanning the channel's whole mass range (independent of
  // age), f >= 0 at every node (gross ejecta are clamped non-negative), and
  // integrated as the piecewise-LINEAR interpolant of the node values:
  //   C(age) = int_{x_to}^{x_K} f_lin(x) dx.
  // Because f_lin >= 0 and only the lower limit depends on age, C is
  // non-decreasing in age by construction, for the total and for every
  // element; interval returns C(e) - C(b) telescope exactly under step
  // subdivision. Accuracy: second order in the fixed spacing (64 intervals
  // over the channel range, ~0.03 dex for the default 1-100 Msun IMF).
  ChannelReturn result;
  const double z = metallicityOf(birth_x);
  double channel_lo = m_config.imf_mass_min_msun;
  double channel_hi = m_config.imf_mass_max_msun;
  if (channel == 0) {
    channel_hi = std::min(channel_hi, m_config.cc_mass_min_msun);
  } else {
    channel_lo = std::max(channel_lo, m_config.cc_mass_min_msun);
    channel_hi = std::min(channel_hi, m_config.cc_mass_max_msun);
  }
  if (!(channel_hi > channel_lo)) return result;
  const double m_to = turnoffMassMsun(age, z);
  if (!(m_to < channel_hi)) return result;  // nothing of this channel has evolved yet
  constexpr int k_intervals = 64;
  const double a = std::log10(channel_lo);
  const double b = std::log10(channel_hi);
  const double h = (b - a) / k_intervals;
  const double x_to = std::log10(std::max(m_to, channel_lo));
  // First node at or above x_to (clamped into [0, K]).
  int j0 = static_cast<int>(std::floor((x_to - a) / h));
  j0 = std::clamp(j0, 0, k_intervals);
  const auto node = [&](int j, double& mass_value, std::array<double, k_chemical_element_count>& elements) {
    const double x = (j == k_intervals) ? b : a + j * h;
    const double m = std::pow(10.0, x);
    const double weight = phi(m) * m * k_ln10;  // dN/dx
    double ejecta = 0.0;
    elements = ejectaElements(channel, m, z, birth_x, ejecta);
    mass_value = weight * ejecta;
    for (double& e : elements) e *= weight;
  };
  double f_hi = 0.0;
  std::array<double, k_chemical_element_count> e_hi{};
  node(k_intervals, f_hi, e_hi);
  // Sum full intervals from the top down (fixed order), then the partial one.
  for (int j = k_intervals - 1; j >= j0; --j) {
    double f_lo = 0.0;
    std::array<double, k_chemical_element_count> e_lo{};
    node(j, f_lo, e_lo);
    const double x_lo = (j == 0) ? a : a + j * h;
    const double x_up = (j + 1 == k_intervals) ? b : a + (j + 1) * h;
    if (x_to <= x_lo) {
      const double width = x_up - x_lo;
      result.mass_fraction += 0.5 * width * (f_lo + f_hi);
      for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
        result.element_fraction[e] += 0.5 * width * (e_lo[e] + e_hi[e]);
      }
    } else if (x_to < x_up) {
      // Partial interval [x_to, x_up] of the non-negative linear interpolant.
      const double t = (x_to - x_lo) / (x_up - x_lo);
      const double width = x_up - x_to;
      result.mass_fraction += 0.5 * width * ((f_lo + t * (f_hi - f_lo)) + f_hi);
      for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
        result.element_fraction[e] += 0.5 * width * ((e_lo[e] + t * (e_hi[e] - e_lo[e])) + e_hi[e]);
      }
    }
    f_hi = f_lo;
    e_hi = e_lo;
  }
  return result;
}

SspModel::ChannelReturn SspModel::massiveReturn(std::size_t channel, double age_begin, double age_end,
                                                const std::array<double, k_chemical_element_count>& birth_x) const {
  ChannelReturn result;
  if (!(age_end > age_begin)) return result;
  const ChannelReturn end = cumulativeMassiveReturn(channel, age_end, birth_x);
  const ChannelReturn begin = cumulativeMassiveReturn(channel, age_begin, birth_x);
  // G3.2: the cumulative return is non-decreasing in age BY CONSTRUCTION
  // (fixed-node non-negative linear interpolant, see cumulativeMassiveReturn).
  // The clamp below is only a round-off guard (last-ulp differences of two
  // nearly equal partial-interval evaluations); any trigger is counted and
  // the mass is re-derived from the clamped elements so element closure
  // holds. It is a diagnostic, not part of the numerical method.
  bool nonmonotone = end.mass_fraction < begin.mass_fraction;
  double element_sum = 0.0;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    const double d = end.element_fraction[e] - begin.element_fraction[e];
    nonmonotone = nonmonotone || d < 0.0;
    result.element_fraction[e] = std::max(d, 0.0);
    element_sum += result.element_fraction[e];
  }
  if (nonmonotone) {
    ++nonmonotone_return_clamps;
    result.mass_fraction = element_sum;
  } else {
    result.mass_fraction = end.mass_fraction - begin.mass_fraction;
  }
  return result;
}

SspModel::ChannelReturn SspModel::iaReturnPerEvent(const std::array<double, k_chemical_element_count>& birth_x) const {
  ChannelReturn result;
  if (m_config.yield_mode == SspYieldMode::kTable) {
    bool clamped = false;
    const SspYieldTable::Row row = m_table->snia(metallicityOf(birth_x), clamped);
    if (clamped) ++table_clamps;
    result.mass_fraction = row.ejecta_msun;
    result.element_fraction = row.net_yield_msun;  // absolute Ia ejecta per element
    return result;
  }
  // REFERENCE ONLY: Chandrasekhar-mass ejecta with the birth mixture.
  result.mass_fraction = k_chandrasekhar_msun;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) result.element_fraction[e] = birth_x[e] * k_chandrasekhar_msun;
  return result;
}

// ---------------------------------------------------------------------------
// Event chronology and persistence
// ---------------------------------------------------------------------------
double eventExponentialDraw(std::uint64_t stream_key, std::uint32_t channel, std::uint32_t ordinal) {
  const std::uint64_t h = splitmix64(stream_key ^ splitmix64((static_cast<std::uint64_t>(channel) << 32U) | ordinal));
  const double u = (static_cast<double>(h >> 11U) + 0.5) * (1.0 / 9007199254740992.0);  // (0,1)
  return -std::log1p(-u);
}

StellarPopulationRecord makeStellarPopulationRecord(const std::array<double, k_chemical_element_count>& birth_x,
                                                    double birth_time_code, std::uint64_t stream_key) {
  StellarPopulationRecord record;
  record.birth_element_mass_fraction = birth_x;
  record.birth_time_code = birth_time_code;
  record.event_stream_key = stream_key;
  record.cc_next_threshold = eventExponentialDraw(stream_key, 0U, 0U);
  record.ia_next_threshold = eventExponentialDraw(stream_key, 1U, 0U);
  record.chronology_time_code = birth_time_code;
  record.chronology_initialized = 1U;
  record.initialized = 1U;
  return record;
}

RealizedEvents realizeEvents(StellarPopulationRecord& record, std::uint32_t channel, double lambda_end) {
  RealizedEvents events;
  double& next = channel == 0U ? record.cc_next_threshold : record.ia_next_threshold;
  std::uint32_t& realized = channel == 0U ? record.cc_realized : record.ia_realized;
  events.first_ordinal = realized;
  // Bounded: the expected count over one step is tiny for any physical star
  // particle; the cap protects against corrupt state.
  while (lambda_end >= next && events.count < 1000000U) {
    ++events.count;
    ++realized;
    next += eventExponentialDraw(record.event_stream_key, channel, realized);
  }
  return events;
}

void advanceStellarPopulationChronology(
    StellarPopulationRecord& record,
    double interval_begin_code,
    double interval_end_code,
    std::uint64_t accepted_step_index) {
  if (record.initialized == 0U) {
    throw std::invalid_argument(
        "stellar population chronology requires an initialized population record");
  }
  if (!(std::isfinite(interval_begin_code) && std::isfinite(interval_end_code) &&
        interval_end_code > interval_begin_code)) {
    throw std::invalid_argument(
        "stellar population chronology requires a finite positive interval");
  }
  const auto tolerance = [](double lhs, double rhs) {
    return 64.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(lhs), std::abs(rhs)});
  };
  if (record.chronology_initialized == 0U) {
    if (interval_begin_code + tolerance(interval_begin_code, record.birth_time_code) <
        record.birth_time_code) {
      throw std::runtime_error(
          "stellar population chronology interval begins before the recorded birth time");
    }
    // Schema-1 and explicitly imported legacy records adopt the first
    // authoritative interval boundary. Their cumulative age lane remains the
    // physical age authority for already elapsed history.
    record.chronology_time_code = interval_begin_code;
    record.chronology_initialized = 1U;
  }
  if (!std::isfinite(record.chronology_time_code) ||
      std::abs(record.chronology_time_code - interval_begin_code) >
          tolerance(record.chronology_time_code, interval_begin_code)) {
    throw std::runtime_error(
        record.chronology_time_code > interval_begin_code
            ? "stellar population chronology rejected an overlapping/replayed interval"
            : "stellar population chronology rejected a gap between accepted intervals");
  }
  record.chronology_time_code = interval_end_code;
  record.chronology_step_index = accepted_step_index;
}

StellarPopulationStore::StellarPopulationStore(const core::SimulationState& state) {
  const core::ModuleSidecarBlock* block = state.sidecars.find(k_stellar_population_module);
  if (block == nullptr) return;
  const bool current_schema =
      block->schema_version == k_stellar_population_schema &&
      block->row_stride_bytes == sizeof(StellarPopulationRecord);
  const bool legacy_schema =
      block->schema_version == 1U &&
      block->row_stride_bytes == sizeof(StellarPopulationRecordV1);
  if (!block->particle_indexed || (!current_schema && !legacy_schema)) {
    throw std::runtime_error("chui_stellar_population_v1 sidecar has an incompatible schema/stride");
  }
  std::vector<std::pair<std::uint64_t, std::size_t>> order;
  for (std::size_t r = 0; r < block->rowCount(); ++r) order.emplace_back(block->particle_id_by_row[r], r);
  std::sort(order.begin(), order.end());
  for (const auto& [id, r] : order) {
    StellarPopulationRecord record{};
    const auto bytes = block->rowPayload(r);
    if (current_schema) {
      std::memcpy(&record, bytes.data(), sizeof(record));
    } else {
      StellarPopulationRecordV1 legacy{};
      std::memcpy(&legacy, bytes.data(), sizeof(legacy));
      record.birth_element_mass_fraction = legacy.birth_element_mass_fraction;
      record.birth_time_code = legacy.birth_time_code;
      record.cc_next_threshold = legacy.cc_next_threshold;
      record.ia_next_threshold = legacy.ia_next_threshold;
      record.cc_realized = legacy.cc_realized;
      record.ia_realized = legacy.ia_realized;
      record.event_stream_key = legacy.event_stream_key;
      record.returned_mass_code = legacy.returned_mass_code;
      record.initialized = legacy.initialized;
      record.reserved = legacy.reserved;
      // The first post-restart source interval adopts its authoritative
      // begin boundary; cumulative ages/ordinals already encode older work.
      record.chronology_initialized = 0U;
    }
    m_ids.push_back(id);
    m_records.push_back(record);
  }
}

bool StellarPopulationStore::has(std::uint64_t particle_id) const {
  return std::binary_search(m_ids.begin(), m_ids.end(), particle_id);
}

StellarPopulationRecord StellarPopulationStore::get(std::uint64_t particle_id) const {
  const auto it = std::lower_bound(m_ids.begin(), m_ids.end(), particle_id);
  if (it == m_ids.end() || *it != particle_id) return {};
  return m_records[static_cast<std::size_t>(it - m_ids.begin())];
}

void StellarPopulationStore::put(std::uint64_t particle_id, const StellarPopulationRecord& record) {
  const auto it = std::lower_bound(m_ids.begin(), m_ids.end(), particle_id);
  const std::size_t index = static_cast<std::size_t>(it - m_ids.begin());
  if (it != m_ids.end() && *it == particle_id) {
    m_records[index] = record;
    return;
  }
  m_ids.insert(it, particle_id);
  m_records.insert(m_records.begin() + static_cast<std::ptrdiff_t>(index), record);
}

void StellarPopulationStore::store(core::SimulationState& state) const {
  core::ModuleSidecarBlock block;
  block.module_name = k_stellar_population_module;
  block.schema_version = k_stellar_population_schema;
  block.particle_indexed = true;
  block.row_stride_bytes = static_cast<std::uint32_t>(sizeof(StellarPopulationRecord));
  block.requirement.kind = core::ModuleSidecarRequirementKind::kSparse;
  // Keep only rows of stars that exist locally (stable particle_id order).
  std::vector<std::uint64_t> live;
  live.reserve(state.star_particles.size());
  for (std::size_t s = 0; s < state.star_particles.size(); ++s) {
    const std::uint32_t p = state.star_particles.particle_index[s];
    if (p < state.particle_sidecar.particle_id.size()) live.push_back(state.particle_sidecar.particle_id[p]);
  }
  std::sort(live.begin(), live.end());
  for (std::size_t i = 0; i < m_ids.size(); ++i) {
    if (!std::binary_search(live.begin(), live.end(), m_ids[i])) continue;
    block.particle_id_by_row.push_back(m_ids[i]);
    const auto* bytes = reinterpret_cast<const std::byte*>(&m_records[i]);
    block.payload.insert(block.payload.end(), bytes, bytes + sizeof(StellarPopulationRecord));
  }
  state.sidecars.upsert(std::move(block));
}

bool forecastStellarEvent(const SspModel& model, const StellarPopulationRecord& record, double age_years,
                          double lookahead_years, double birth_mass_msun, bool stochastic, double expected_threshold) {
  if (record.initialized == 0U || !(lookahead_years > 0.0) || !(birth_mass_msun > 0.0)) return false;
  double z = 0.0;
  for (std::size_t e = 2; e < k_chemical_element_count; ++e) z += record.birth_element_mass_fraction[e];
  const double age_now = std::max(age_years, 0.0);
  const double age_end = age_now + lookahead_years;
  const double cc_end = model.ccCumulativeNumber(age_end, z) * birth_mass_msun;
  const double ia_end = model.iaCumulativeNumber(age_end) * birth_mass_msun;
  if (stochastic) {
    StellarPopulationRecord copy = record;
    return realizeEvents(copy, 0U, cc_end).count + realizeEvents(copy, 1U, ia_end).count > 0U;
  }
  const double expected = (cc_end - model.ccCumulativeNumber(age_now, z) * birth_mass_msun) +
      (ia_end - model.iaCumulativeNumber(age_now) * birth_mass_msun);
  return expected >= expected_threshold;
}

}  // namespace cosmosim::physics
