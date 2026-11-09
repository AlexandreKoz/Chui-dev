#include "cosmosim/physics/stellar_wind.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "cosmosim/physics/radiation_background.hpp"

namespace cosmosim::physics {

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
  return loadFromString(text.str(), path);
}

StellarWindTable StellarWindTable::loadFromString(const std::string& text, const std::string& label) {
  StellarWindTable table;
  std::istringstream lines(text);
  std::string line;
  bool header = false;
  std::size_t expected_ages = 0U;
  const auto fail = [&](const std::string& why) {
    throw std::runtime_error("stellar wind table '" + label + "': " + why);
  };
  while (std::getline(lines, line)) {
    if (line.empty()) continue;
    if (line.rfind("# chui_ssp_wind_table_v1", 0) == 0) { header = true; continue; }
    if (line.rfind("# source =", 0) == 0) { table.m_source = line.substr(10); continue; }
    if (line.rfind("# limitation =", 0) == 0) {
      if (!table.m_limitations.empty()) table.m_limitations += "; ";
      const std::size_t first = line.find_first_not_of(" \t", 14);
      if (first != std::string::npos) table.m_limitations += line.substr(first);
      continue;
    }
    if (line[0] == '#') continue;
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
  if (!header) fail("missing '# chui_ssp_wind_table_v1' header");
  if (table.m_source.empty()) fail("missing '# source =' provenance");
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
