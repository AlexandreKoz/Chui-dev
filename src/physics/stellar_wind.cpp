#include "cosmosim/physics/stellar_wind.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "cosmosim/physics/radiation_background.hpp"

namespace cosmosim::physics {

namespace {

[[nodiscard]] std::string trimText(const std::string& text) {
  const auto begin = text.find_first_not_of(" \t\r");
  if (begin == std::string::npos) return {};
  return text.substr(begin, text.find_last_not_of(" \t\r") - begin + 1U);
}

// Header keys owned by the v2 metadata contract. `license` is also tolerated
// (informational) in v1 because historical converters wrote it.
constexpr std::array<std::string_view, 8> k_v2_metadata_keys{
    "license", "imf", "imf_mass_min_msun", "imf_mass_max_msun", "normalization",
    "energy_unit", "momentum_unit", "metallicity_definition"};

[[nodiscard]] bool isV2MetadataKey(std::string_view key) noexcept {
  return std::find(k_v2_metadata_keys.begin(), k_v2_metadata_keys.end(), key) != k_v2_metadata_keys.end();
}

[[nodiscard]] bool isIdentityToken(const std::string& value) noexcept {
  if (value.empty() || value.size() > 64U) return false;
  return std::all_of(value.begin(), value.end(), [](char c) {
    return (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_';
  });
}

[[nodiscard]] std::string formatNumber(double value) {
  std::ostringstream text;
  text.precision(9);
  text << value;
  return text.str();
}

}  // namespace

const char* stellarWindCompatibilityPolicyName(StellarWindCompatibilityPolicy policy) noexcept {
  switch (policy) {
    case StellarWindCompatibilityPolicy::kStrict: return "strict";
    case StellarWindCompatibilityPolicy::kApproximateResearch: return "approximate_research";
  }
  return "unknown";
}

StellarWindTable StellarWindTable::loadFromFile(const std::string& path, const std::string& expected_sha256) {
  if (path.empty()) {
    throw std::runtime_error("stellar wind table: ssp_wind_table_path is empty");
  }
  const std::string actual = sha256FileHex(path);
  if (expected_sha256.empty() || actual != expected_sha256) {
    throw std::runtime_error("stellar wind table '" + path + "': SHA-256 " + actual +
                             " does not match the configured ssp_wind_table_sha256 (refused)");
  }
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("stellar wind table: cannot open '" + path + "'");
  }
  std::ostringstream text;
  text << input.rdbuf();
  StellarWindTable table = loadFromString(text.str(), path);
  table.m_sha256 = actual;
  return table;
}

StellarWindTable StellarWindTable::loadFromString(const std::string& text, const std::string& label) {
  StellarWindTable table;
  std::istringstream lines(text);
  std::string line;
  bool header = false;
  std::size_t source_count = 0U;
  std::size_t expected_ages = 0U;
  std::map<std::string, std::vector<std::string>> declared;
  const auto fail = [&](const std::string& why) {
    throw std::runtime_error("stellar wind table '" + label + "': " + why);
  };
  while (std::getline(lines, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') {
      const std::string body = trimText(line.substr(1));
      if (body.rfind("chui_ssp_wind_table_", 0) == 0) {
        StellarWindTableSchema schema = StellarWindTableSchema::kV1;
        if (body == "chui_ssp_wind_table_v1") {
          schema = StellarWindTableSchema::kV1;
        } else if (body == "chui_ssp_wind_table_v2") {
          schema = StellarWindTableSchema::kV2;
        } else {
          fail("unsupported schema header '" + body + "' (expected chui_ssp_wind_table_v1 or _v2)");
        }
        if (header) fail("more than one schema header (contradictory declarations)");
        header = true;
        table.m_metadata.schema = schema;
        continue;
      }
      if (line.rfind("# source =", 0) == 0) {
        table.m_source = line.substr(10);  // v1 byte-compatible provenance text
        ++source_count;
        continue;
      }
      if (line.rfind("# limitation =", 0) == 0) {
        if (!table.m_limitations.empty()) table.m_limitations += "; ";
        const std::size_t first = line.find_first_not_of(" \t", 14);
        if (first != std::string::npos) table.m_limitations += line.substr(first);
        continue;
      }
      const auto eq = body.find('=');
      if (eq != std::string::npos) {
        const std::string key = trimText(body.substr(0, eq));
        if (isV2MetadataKey(key)) declared[key].push_back(trimText(body.substr(eq + 1U)));
      }
      continue;
    }
    std::istringstream fields(line);
    std::string word;
    fields >> word;
    if (word == "metallicities") {
      std::size_t n = 0U;
      fields >> n;
      table.m_z.resize(n);
      for (double& z : table.m_z) fields >> z;
      if (!fields || n == 0U) fail("malformed metallicities line");
      continue;
    }
    if (word == "ages") {
      fields >> expected_ages;
      if (!fields || expected_ages < 2U) fail("ages must be >= 2");
      continue;
    }
    std::istringstream row(line);
    double log_age = 0.0;
    row >> log_age;
    if (!row) fail("malformed data row");
    table.m_log_age.push_back(log_age);
    for (std::size_t z = 0; z < table.m_z.size(); ++z) {
      double energy = 0.0, momentum = 0.0;
      row >> energy >> momentum;
      if (!row) fail("data row has fewer than 2 x metallicity columns");
      table.m_values.push_back(energy);
      table.m_values.push_back(momentum);
    }
  }
  if (!header) fail("missing '# chui_ssp_wind_table_v1' or '# chui_ssp_wind_table_v2' header");
  if (table.m_source.empty()) fail("missing '# source =' provenance");

  StellarWindTableMetadata& meta = table.m_metadata;
  const auto join = [](const std::vector<std::string>& values) {
    std::string out;
    for (const std::string& v : values) out += (out.empty() ? "" : "; ") + v;
    return out;
  };
  if (meta.schema == StellarWindTableSchema::kV1) {
    // Historical schema: loading semantics unchanged. Its IMF is undeclared;
    // a v2 metadata key here would make the schema ambiguous.
    for (const auto& [key, values] : declared) {
      if (key != "license") {
        fail("v2 metadata key '" + key + "' in a chui_ssp_wind_table_v1 table; declare '# chui_ssp_wind_table_v2' "
             "and the complete v2 metadata set instead");
      }
      meta.license = join(values);
    }
  } else {
    if (source_count != 1U) fail("chui_ssp_wind_table_v2 requires exactly one '# source =' declaration");
    const auto single = [&](std::string_view key) -> std::string {
      const auto it = declared.find(std::string(key));
      if (it == declared.end()) fail("chui_ssp_wind_table_v2 requires '# " + std::string(key) + " ='");
      if (it->second.size() != 1U) {
        fail("'# " + std::string(key) + " =' declared " + std::to_string(it->second.size()) +
             " times (contradictory declarations)");
      }
      if (it->second.front().empty()) fail("'# " + std::string(key) + " =' is empty");
      return it->second.front();
    };
    const auto positive = [&](std::string_view key) {
      const std::string value = single(key);
      double parsed = 0.0;
      std::size_t used = 0U;
      try {
        parsed = std::stod(value, &used);
      } catch (const std::exception&) {
        fail("'# " + std::string(key) + " = " + value + "' is not a number");
      }
      if (used != value.size() || !std::isfinite(parsed) || !(parsed > 0.0)) {
        fail("'# " + std::string(key) + " = " + value + "' must be a finite positive number");
      }
      return parsed;
    };
    const auto exact = [&](std::string_view key, std::string_view expected) {
      const std::string value = single(key);
      if (value != expected) {
        fail("'# " + std::string(key) + " = " + value + "' unsupported (the only defined value is '" +
             std::string(expected) + "')");
      }
      return value;
    };
    meta.license = single("license");
    meta.imf_identity = single("imf");
    if (!isIdentityToken(meta.imf_identity)) {
      fail("'# imf = " + meta.imf_identity + "' must be an identity token [a-z0-9_]+ (<= 64 characters)");
    }
    meta.imf_mass_min_msun = positive("imf_mass_min_msun");
    meta.imf_mass_max_msun = positive("imf_mass_max_msun");
    if (!(meta.imf_mass_max_msun > meta.imf_mass_min_msun)) {
      fail("imf_mass_max_msun must exceed imf_mass_min_msun");
    }
    meta.normalization = exact("normalization", "per_msun_initial_stellar_mass");
    meta.energy_unit = exact("energy_unit", "erg_per_msun");
    meta.momentum_unit = exact("momentum_unit", "msun_km_s_per_msun");
    meta.metallicity_definition = exact("metallicity_definition", "metal_mass_fraction");
  }
  if (table.m_z.empty() || table.m_log_age.size() != expected_ages) fail("row count differs from 'ages'");
  for (std::size_t i = 1; i < table.m_z.size(); ++i) {
    if (!(table.m_z[i] > table.m_z[i - 1])) fail("metallicities must be strictly increasing");
  }
  if (table.m_z.front() < 0.0) fail("negative metallicity");
  const std::size_t nz = table.m_z.size();
  for (std::size_t a = 0; a < table.m_log_age.size(); ++a) {
    if (a > 0 && !(table.m_log_age[a] > table.m_log_age[a - 1])) fail("ages must be strictly increasing");
    for (std::size_t z = 0; z < nz; ++z) {
      for (std::size_t q = 0; q < 2U; ++q) {
        const double v = table.m_values[(a * nz + z) * 2U + q];
        if (!std::isfinite(v) || v < 0.0) fail("cumulative budgets must be finite and non-negative");
        if (a > 0 && v < table.m_values[((a - 1) * nz + z) * 2U + q]) {
          fail("cumulative budgets must be non-decreasing in age");
        }
      }
    }
  }
  return table;
}

StellarWindCompatibilityReport assessStellarWindCompatibility(
    const StellarWindTable& table, std::string_view ssp_imf_identity, double ssp_imf_mass_min_msun,
    double ssp_imf_mass_max_msun, StellarWindCompatibilityPolicy policy) {
  StellarWindCompatibilityReport report;
  const StellarWindTableMetadata& meta = table.metadata();
  std::vector<std::string> issues;
  if (meta.schema == StellarWindTableSchema::kV1) {
    issues.push_back("IMF undeclared (chui_ssp_wind_table_v1 has no machine-readable IMF; its free-text limitations "
                     "are not scientific metadata)");
  } else {
    report.imf_declared = meta.imf_identity != "unknown";
    report.imf_identity_matches = report.imf_declared && meta.imf_identity == ssp_imf_identity;
    const auto same = [](double a, double b) { return std::abs(a - b) <= 1.0e-9 * std::max(1.0, std::abs(b)); };
    report.imf_mass_limits_match =
        same(meta.imf_mass_min_msun, ssp_imf_mass_min_msun) && same(meta.imf_mass_max_msun, ssp_imf_mass_max_msun);
    if (!report.imf_declared) {
      issues.push_back("table IMF declared 'unknown'");
    } else if (!report.imf_identity_matches) {
      issues.push_back("table IMF '" + meta.imf_identity + "' differs from the CHUI SSP IMF '" +
                       std::string(ssp_imf_identity) + "'");
    }
    if (!report.imf_mass_limits_match) {
      issues.push_back("table IMF mass limits " + formatNumber(meta.imf_mass_min_msun) + ".." +
                       formatNumber(meta.imf_mass_max_msun) + " Msun differ from the configured SSP limits " +
                       formatNumber(ssp_imf_mass_min_msun) + ".." + formatNumber(ssp_imf_mass_max_msun) +
                       " Msun (per-unit-mass normalization differs)");
    }
  }
  report.single_metallicity = table.metallicityCount() == 1U;
  if (report.single_metallicity) {
    issues.push_back("single tabulated metallicity Z=" + formatNumber(table.metallicityMin()) +
                     " (every other metallicity uses that column)");
  }
  std::string joined;
  for (const std::string& issue : issues) joined += (joined.empty() ? "" : "; ") + issue;
  if (policy == StellarWindCompatibilityPolicy::kStrict) {
    report.admitted = issues.empty();
    report.refusal = joined;
    report.approximations = issues.empty() ? "none" : "";
  } else {
    report.admitted = true;
    report.approximations = issues.empty() ? "none" : joined;
  }
  return report;
}

double StellarWindTable::cumulative(std::size_t quantity, double age_years, double metallicity, bool& clamped) const {
  if (!(age_years > 0.0)) return 0.0;
  const std::size_t nz = m_z.size();
  double z = metallicity;
  if (z < m_z.front() || z > m_z.back()) {
    // D2: a single-metallicity table (e.g. one Starburst99 track set) used at
    // any other Z is an extrapolation and is REPORTED (previously silent when
    // nz == 1). Counted by the caller (m_wind_metallicity_clamps).
    clamped = true;
    z = std::clamp(z, m_z.front(), m_z.back());
  }
  std::size_t zi = 0U;
  double zw = 0.0;
  if (nz > 1U) {
    zi = static_cast<std::size_t>(std::upper_bound(m_z.begin(), m_z.end(), z) - m_z.begin());
    zi = std::clamp<std::size_t>(zi, 1U, nz - 1U) - 1U;
    zw = (z - m_z[zi]) / (m_z[zi + 1U] - m_z[zi]);
  }
  const auto at = [&](std::size_t a, std::size_t zz) { return m_values[(a * nz + zz) * 2U + quantity]; };
  const auto atZ = [&](std::size_t a) { return nz > 1U ? (1.0 - zw) * at(a, zi) + zw * at(a, zi + 1U) : at(a, 0U); };
  const double log_age = std::log10(age_years);
  if (log_age <= m_log_age.front()) {
    return atZ(0U) * age_years / std::pow(10.0, m_log_age.front());
  }
  if (log_age >= m_log_age.back()) {
    return atZ(m_log_age.size() - 1U);
  }
  const std::size_t hi = static_cast<std::size_t>(
      std::upper_bound(m_log_age.begin(), m_log_age.end(), log_age) - m_log_age.begin());
  const std::size_t lo = hi - 1U;
  const double w = (log_age - m_log_age[lo]) / (m_log_age[hi] - m_log_age[lo]);
  // Linear in log age between non-decreasing nodes: the interpolant is
  // non-decreasing, so interval budgets are never negative.
  return (1.0 - w) * atZ(lo) + w * atZ(hi);
}

StellarWindIntervalBudget StellarWindTable::interval(double age_begin_years, double age_end_years,
                                                     double metallicity) const {
  if (!(age_end_years >= age_begin_years)) {
    throw std::invalid_argument("stellar wind interval: end age precedes begin age");
  }
  StellarWindIntervalBudget budget;
  bool clamped = false;
  budget.energy_erg_per_msun = std::max(0.0, cumulative(0U, age_end_years, metallicity, clamped) -
                                                 cumulative(0U, age_begin_years, metallicity, clamped));
  budget.momentum_msun_kms_per_msun = std::max(0.0, cumulative(1U, age_end_years, metallicity, clamped) -
                                                        cumulative(1U, age_begin_years, metallicity, clamped));
  budget.metallicity_clamped = clamped;
  return budget;
}

}  // namespace cosmosim::physics
