#include "cosmosim/physics/radiation_background.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <sstream>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/physics/source_interval.hpp"

namespace cosmosim::physics {

std::string sha256FileHex(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw SourceModelUnavailableError("cannot open data file '" + path + "' for checksum");
  }
  const std::string bytes((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  return core::strongConfigHashSha256Hex(bytes);
}

UvBackgroundTable UvBackgroundTable::loadTreecool(const std::string& path,
                                                  const std::string& expected_sha256_hex) {
  if (path.empty()) {
    throw SourceModelUnavailableError(
        "UV background table path is empty. Supply a TREECOOL-format file (e.g. the Faucher-Giguere 2020 "
        "TREECOOL distributed with that paper, or HM12) via physics.uv_background_table_path; CHUI does not "
        "bundle or synthesize UV background data");
  }
  std::ifstream input(path);
  if (!input) {
    throw SourceModelUnavailableError("cannot open UV background table '" + path + "'");
  }
  std::vector<double> log_z;
  std::vector<UvBackgroundRates> rows;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') continue;
    std::istringstream stream(line);
    double lz = 0.0;
    UvBackgroundRates r;
    stream >> lz >> r.gamma_hi_s >> r.gamma_hei_s >> r.gamma_heii_s >> r.heat_hi_erg_s >> r.heat_hei_erg_s >>
        r.heat_heii_erg_s;
    if (stream.fail()) {
      throw SourceModelUnavailableError("malformed TREECOOL row " + std::to_string(line_number) + " in '" + path +
                                        "' (expected 7 numeric columns)");
    }
    for (const double v : {lz, r.gamma_hi_s, r.gamma_hei_s, r.gamma_heii_s, r.heat_hi_erg_s, r.heat_hei_erg_s,
                           r.heat_heii_erg_s}) {
      if (!std::isfinite(v) || v < 0.0) {
        throw SourceModelUnavailableError("non-finite or negative TREECOOL value at row " +
                                          std::to_string(line_number) + " in '" + path + "'");
      }
    }
    log_z.push_back(lz);
    rows.push_back(r);
  }
  UvBackgroundTable table = fromRows(std::move(log_z), std::move(rows));
  table.m_provenance.source_path = path;
  table.m_provenance.sha256_hex = sha256FileHex(path);
  if (!expected_sha256_hex.empty() && expected_sha256_hex != table.m_provenance.sha256_hex) {
    throw SourceModelUnavailableError("UV background table '" + path + "' SHA-256 " +
                                      table.m_provenance.sha256_hex + " does not match the configured " +
                                      expected_sha256_hex);
  }
  return table;
}

UvBackgroundTable UvBackgroundTable::fromRows(std::vector<double> log10_one_plus_z,
                                              std::vector<UvBackgroundRates> rows) {
  if (log10_one_plus_z.size() < 2U || log10_one_plus_z.size() != rows.size()) {
    throw SourceModelUnavailableError("UV background table needs >= 2 rows with matching columns");
  }
  for (std::size_t i = 1; i < log10_one_plus_z.size(); ++i) {
    if (!(log10_one_plus_z[i] > log10_one_plus_z[i - 1U])) {
      throw SourceModelUnavailableError("UV background table log10(1+z) must be strictly increasing");
    }
  }
  UvBackgroundTable table;
  table.m_log10_one_plus_z = std::move(log10_one_plus_z);
  table.m_rows = std::move(rows);
  table.m_provenance.row_count = table.m_rows.size();
  table.m_provenance.redshift_min = std::pow(10.0, table.m_log10_one_plus_z.front()) - 1.0;
  table.m_provenance.redshift_max = std::pow(10.0, table.m_log10_one_plus_z.back()) - 1.0;
  return table;
}

UvBackgroundRates UvBackgroundTable::ratesAtRedshift(double redshift) const {
  if (empty()) return {};
  const double lz = std::log10(1.0 + std::max(redshift, 0.0));
  if (lz > m_log10_one_plus_z.back()) return {};  // before the background switches on
  if (lz <= m_log10_one_plus_z.front()) return m_rows.front();
  const auto upper = std::upper_bound(m_log10_one_plus_z.begin(), m_log10_one_plus_z.end(), lz);
  const std::size_t hi = static_cast<std::size_t>(std::distance(m_log10_one_plus_z.begin(), upper));
  const std::size_t lo = hi - 1U;
  const double w = (lz - m_log10_one_plus_z[lo]) / (m_log10_one_plus_z[hi] - m_log10_one_plus_z[lo]);
  // Rates are positive and span decades: interpolate log10 where both ends are
  // positive, linearly otherwise (zero entries are legal and must stay zero).
  const auto interp = [w](double a, double b) {
    if (a > 0.0 && b > 0.0) return std::pow(10.0, (1.0 - w) * std::log10(a) + w * std::log10(b));
    return (1.0 - w) * a + w * b;
  };
  const UvBackgroundRates& a = m_rows[lo];
  const UvBackgroundRates& b = m_rows[hi];
  return {interp(a.gamma_hi_s, b.gamma_hi_s),       interp(a.gamma_hei_s, b.gamma_hei_s),
          interp(a.gamma_heii_s, b.gamma_heii_s),   interp(a.heat_hi_erg_s, b.heat_hi_erg_s),
          interp(a.heat_hei_erg_s, b.heat_hei_erg_s), interp(a.heat_heii_erg_s, b.heat_heii_erg_s)};
}

double rahmati2013ShieldingFactor(double hydrogen_number_density_cm3, double temperature_k, double gamma_hi_uvb_s,
                                  double effective_cross_section_cm2) {
  if (!(gamma_hi_uvb_s > 0.0) || !(hydrogen_number_density_cm3 > 0.0)) return 1.0;
  const double n0 = 6.73e-3 * std::pow(effective_cross_section_cm2 / 2.49e-18, -2.0 / 3.0) *
      std::pow(std::max(temperature_k, 1.0) / 1.0e4, 0.17) * std::pow(gamma_hi_uvb_s / 1.0e-12, 2.0 / 3.0);
  const double x = hydrogen_number_density_cm3 / n0;
  return 0.98 * std::pow(1.0 + std::pow(x, 1.64), -2.28) + 0.02 * std::pow(1.0 + x, -0.84);
}

}  // namespace cosmosim::physics
