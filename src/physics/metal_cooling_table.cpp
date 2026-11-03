#include "cosmosim/physics/metal_cooling_table.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>

#include "cosmosim/physics/radiation_background.hpp"
#include "cosmosim/physics/source_interval.hpp"

namespace cosmosim::physics {
namespace {

constexpr double k_asinh_scale = 1.0e-30;

[[nodiscard]] std::string trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r");
  if (b == std::string::npos) return {};
  const auto e = s.find_last_not_of(" \t\r");
  return s.substr(b, e - b + 1U);
}

[[nodiscard]] ChemicalElement elementFromSymbol(const std::string& symbol, const std::string& label) {
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    if (chemicalElementTable()[e].symbol == symbol) return static_cast<ChemicalElement>(e);
  }
  throw SourceModelUnavailableError("metal cooling table '" + label + "' names unknown element '" + symbol + "'");
}

void requireIncreasing(const std::vector<double>& axis, const char* name, const std::string& label) {
  if (axis.empty()) throw SourceModelUnavailableError("metal cooling table '" + label + "' lacks axis " + name);
  for (std::size_t i = 1; i < axis.size(); ++i) {
    if (!(axis[i] > axis[i - 1U])) {
      throw SourceModelUnavailableError("metal cooling table '" + label + "' axis " + name +
                                        " is not strictly increasing");
    }
  }
}

// Coordinate-match tolerance for range checks (Part 2.5 interpolation
// safety). Axis values are read from ASCII written with >= 9 significant
// digits (convert_wss09_metal_cooling.py uses %.9g: relative representation
// error <= 5e-10), and query coordinates are log10 T, log10 n_H (double
// round-off ~1e-15 relative) or z = 1/a - 1 (round-off ~1e-15 (1+z)). An
// absolute tolerance of 1e-8 x max(1, |axis value|) therefore exceeds every
// representation/round-off error by >= 20x, while remaining 1e6 times
// smaller than any physically meaningful offset (a WSS09 grid step is 0.02
// dex in T, 0.1 dex in n_H, and >= 0.05 in z). A query within this tolerance
// of a boundary is treated as ON the grid (snapped, not counted as a clamp).
[[nodiscard]] double axisTolerance(double value) noexcept { return 1.0e-8 * std::max(1.0, std::abs(value)); }

// Locate interval and weight along one axis. A SINGLETON axis is a single
// tabulated coordinate, not a wildcard: a query that does not match it is out
// of range exactly like a query beyond the ends of a longer axis (rejected,
// or clamped and reported when clamping is allowed). Before Part 2.5 a
// singleton axis accepted every query, so a table converted from WSS09
// z_0.000.hdf5 alone silently supplied z = 0 rates at any redshift.
[[nodiscard]] bool locate(const std::vector<double>& axis, double x, std::size_t& lo, double& w,
                          bool clamp_allowed, bool& clamped) {
  if (!std::isfinite(x)) return false;
  const double tol_lo = axisTolerance(axis.front());
  const double tol_hi = axisTolerance(axis.back());
  if (x < axis.front() - tol_lo || x > axis.back() + tol_hi) {
    if (!clamp_allowed) return false;
    clamped = true;
  }
  x = std::clamp(x, axis.front(), axis.back());  // snaps in-tolerance queries; clamps reported above
  if (axis.size() == 1U) {
    lo = 0U;
    w = 0.0;
    return true;
  }
  const auto upper = std::upper_bound(axis.begin(), axis.end(), x);
  std::size_t hi = static_cast<std::size_t>(std::distance(axis.begin(), upper));
  hi = std::clamp<std::size_t>(hi, 1U, axis.size() - 1U);
  lo = hi - 1U;
  w = (x - axis[lo]) / (axis[hi] - axis[lo]);
  return true;
}

}  // namespace

MetalCoolingTable MetalCoolingTable::loadFromFile(const std::string& path, const std::string& expected_sha256_hex) {
  if (path.empty()) {
    throw SourceModelUnavailableError(
        "metal cooling table path is empty. Convert a licensed element-resolved table (e.g. Wiersma, Schaye & "
        "Smith 2009, MNRAS 393, 99) into chui_metal_cooling_table_v1 and set physics.metal_cooling_table_path; "
        "CHUI never synthesizes production cooling tables");
  }
  std::ifstream input(path, std::ios::binary);
  if (!input) throw SourceModelUnavailableError("cannot open metal cooling table '" + path + "'");
  const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  MetalCoolingTable table = loadFromString(text, path);
  table.m_provenance.source_path = path;
  table.m_provenance.sha256_hex = sha256FileHex(path);
  if (!expected_sha256_hex.empty() && expected_sha256_hex != table.m_provenance.sha256_hex) {
    throw SourceModelUnavailableError("metal cooling table '" + path + "' SHA-256 mismatch: file " +
                                      table.m_provenance.sha256_hex + " vs configured " + expected_sha256_hex);
  }
  return table;
}

MetalCoolingTable MetalCoolingTable::loadFromString(const std::string& text, const std::string& label) {
  MetalCoolingTable table;
  std::istringstream stream(text);
  std::string line;
  bool magic = false;
  std::map<std::string, std::vector<double>> lambda_blocks;
  std::vector<std::string> element_symbols;
  while (std::getline(stream, line)) {
    const std::string t = trim(line);
    if (t.empty()) continue;
    if (t[0] == '#') {
      const std::string body = trim(t.substr(1));
      if (body == "chui_metal_cooling_table_v1") magic = true;
      const auto eq = body.find('=');
      if (eq != std::string::npos) {
        const std::string key = trim(body.substr(0, eq));
        const std::string value = trim(body.substr(eq + 1U));
        if (key == "citation") table.m_provenance.citation = value;
        if (key == "license") table.m_provenance.license = value;
        if (key == "solar_reference") table.m_provenance.solar_reference = value;
        if (key == "uv_background") table.m_provenance.uv_background = value;
        if (key == "uv_background_treecool_sha256") table.m_provenance.uv_background_treecool_sha256 = value;
      }
      continue;
    }
    std::istringstream row(t);
    std::string key;
    row >> key;
    if (key == "elements") {
      std::string symbol;
      while (row >> symbol) element_symbols.push_back(symbol);
      continue;
    }
    std::string block_name;
    if (key == "lambda") row >> block_name;
    std::vector<double> values;
    double v = 0.0;
    while (row >> v) values.push_back(v);
    if (!row.eof()) throw SourceModelUnavailableError("metal cooling table '" + label + "': non-numeric data in '" + key + "'");
    for (const double x : values) {
      if (!std::isfinite(x)) throw SourceModelUnavailableError("metal cooling table '" + label + "': non-finite value");
    }
    if (key == "log10_T") table.m_log_t = values;
    else if (key == "log10_nH") table.m_log_n = values;
    else if (key == "redshift") table.m_redshift = values;
    else if (key == "reference_mass_fraction") table.m_reference_mass_fraction = values;
    else if (key == "reference_hydrogen_mass_fraction") {
      if (values.size() != 1U) {
        throw SourceModelUnavailableError("metal cooling table '" + label +
                                          "': reference_hydrogen_mass_fraction must hold exactly one value");
      }
      table.m_reference_hydrogen_mass_fraction = values[0];
    }
    else if (key == "ne_over_nH") table.m_ne_over_nh = values;
    else if (key == "lambda") lambda_blocks[block_name] = values;
    else throw SourceModelUnavailableError("metal cooling table '" + label + "': unknown block '" + key + "'");
  }
  if (!magic) throw SourceModelUnavailableError("metal cooling table '" + label + "' lacks the chui_metal_cooling_table_v1 header");
  if (table.m_provenance.citation.empty() || table.m_provenance.license.empty() ||
      table.m_provenance.solar_reference.empty()) {
    throw SourceModelUnavailableError("metal cooling table '" + label +
                                      "' must declare citation, license and solar_reference (provenance is mandatory)");
  }
  if (table.m_provenance.uv_background.empty()) {
    throw SourceModelUnavailableError(
        "metal cooling table '" + label +
        "' must declare uv_background (the UV background embedded in the table, or 'none' for a collisional "
        "table); required since Part 2.5 so that UVB mismatches cannot pass silently");
  }
  if (!(table.m_reference_hydrogen_mass_fraction > 0.0) || !(table.m_reference_hydrogen_mass_fraction < 1.0)) {
    throw SourceModelUnavailableError(
        "metal cooling table '" + label +
        "' must declare reference_hydrogen_mass_fraction in (0, 1): WSS09 eq. 4 scales by abundances relative to "
        "hydrogen, (X_e/X_H)/(X_e,ref/X_H,ref)");
  }
  requireIncreasing(table.m_log_t, "log10_T", label);
  requireIncreasing(table.m_log_n, "log10_nH", label);
  requireIncreasing(table.m_redshift, "redshift", label);
  const std::size_t grid = table.m_log_t.size() * table.m_log_n.size() * table.m_redshift.size();
  table.m_provenance.grid_points = grid;
  if (element_symbols.empty() || element_symbols.size() != table.m_reference_mass_fraction.size()) {
    throw SourceModelUnavailableError("metal cooling table '" + label +
                                      "': elements and reference_mass_fraction must have equal non-zero length");
  }
  if (!table.m_ne_over_nh.empty() && table.m_ne_over_nh.size() != grid) {
    throw SourceModelUnavailableError("metal cooling table '" + label + "': ne_over_nH block has wrong size");
  }
  for (double& value : table.m_ne_over_nh) value = std::asinh(value);  // stored transformed
  for (std::size_t i = 0; i < element_symbols.size(); ++i) {
    const ChemicalElement element = elementFromSymbol(element_symbols[i], label);
    if (element == ChemicalElement::kH || element == ChemicalElement::kHe) {
      throw SourceModelUnavailableError(
          "metal cooling table '" + label +
          "' declares an H/He column: primordial cooling is owned by the non-equilibrium network and would be "
          "double counted");
    }
    if (!(table.m_reference_mass_fraction[i] > 0.0)) {
      throw SourceModelUnavailableError("metal cooling table '" + label + "': reference mass fraction must be > 0");
    }
    const auto found = lambda_blocks.find(element_symbols[i]);
    if (found == lambda_blocks.end() || found->second.size() != grid) {
      throw SourceModelUnavailableError("metal cooling table '" + label + "': missing or mis-sized lambda block for " +
                                        element_symbols[i]);
    }
    std::vector<double> transformed(grid);
    for (std::size_t g = 0; g < grid; ++g) transformed[g] = std::asinh(found->second[g] / k_asinh_scale);
    table.m_elements.push_back(element);
    table.m_lambda_transformed.push_back(std::move(transformed));
  }
  return table;
}

double MetalCoolingTable::temperatureMinK() const noexcept { return std::pow(10.0, m_log_t.front()); }
double MetalCoolingTable::temperatureMaxK() const noexcept { return std::pow(10.0, m_log_t.back()); }

double MetalCoolingTable::interpolate(const std::vector<double>& block, double z, double lt, double ln,
                                      bool& clamped, MetalTableRangePolicy policy) const {
  std::size_t iz = 0, it = 0, in = 0;
  double wz = 0.0, wt = 0.0, wn = 0.0;
  const bool clamp_dn = policy == MetalTableRangePolicy::kClampDensityAndRedshift;
  if (!locate(m_log_t, lt, it, wt, false, clamped)) {
    throw SourceModelUnavailableError("metal cooling table '" + m_provenance.source_path +
                                      "': temperature outside the table (forbidden extrapolation)");
  }
  if (!locate(m_log_n, ln, in, wn, clamp_dn, clamped) || !locate(m_redshift, z, iz, wz, clamp_dn, clamped)) {
    throw SourceModelUnavailableError("metal cooling table '" + m_provenance.source_path +
                                      "': density/redshift outside the table and range policy is reject");
  }
  const std::size_t nt = m_log_t.size(), nn = m_log_n.size(), nz = m_redshift.size();
  const auto at = [&](std::size_t a, std::size_t b, std::size_t c) { return block[(a * nt + b) * nn + c]; };
  double result = 0.0;
  for (int dz = 0; dz < 2; ++dz) {
    const std::size_t z_i = std::min(iz + static_cast<std::size_t>(dz), nz - 1U);
    const double fz = dz == 0 ? 1.0 - wz : wz;
    if (fz == 0.0) continue;
    for (int dt = 0; dt < 2; ++dt) {
      const std::size_t t_i = std::min(it + static_cast<std::size_t>(dt), nt - 1U);
      const double ft = dt == 0 ? 1.0 - wt : wt;
      if (ft == 0.0) continue;
      for (int dn = 0; dn < 2; ++dn) {
        const std::size_t n_i = std::min(in + static_cast<std::size_t>(dn), nn - 1U);
        const double fn = dn == 0 ? 1.0 - wn : wn;
        if (fn == 0.0) continue;
        result += fz * ft * fn * at(z_i, t_i, n_i);
      }
    }
  }
  return result;
}

MetalCoolingTable::Result MetalCoolingTable::evaluate(const Query& query, MetalTableRangePolicy policy) const {
  Result result;
  const double lt = std::log10(std::max(query.temperature_k, 1.0));
  const double ln = std::log10(std::max(query.hydrogen_number_density_cm3, 1.0e-300));
  const double nh = query.hydrogen_number_density_cm3;
  double electron_correction = 1.0;
  if (query.electron_per_hydrogen >= 0.0 && !m_ne_over_nh.empty()) {
    // Precomputed asinh transform of the equilibrium n_e/n_H block (no
    // per-query allocation).
    bool c = false;
    const double eq = std::sinh(interpolate(m_ne_over_nh, query.redshift, lt, ln, c, policy));
    if (eq > 0.0) electron_correction = query.electron_per_hydrogen / eq;
  }
  const double x_h = query.gas_mass_fraction[static_cast<std::size_t>(ChemicalElement::kH)];
  const double hydrogen_scale = x_h > 0.0 ? m_reference_hydrogen_mass_fraction / x_h : 0.0;
  for (std::size_t i = 0; i < m_elements.size(); ++i) {
    const double transformed = interpolate(m_lambda_transformed[i], query.redshift, lt, ln, result.clamped, policy);
    const double lambda_ref = k_asinh_scale * std::sinh(transformed);
    // WSS09 eq. 4: 10^[e/H] = (X_e/X_H) / (X_e,ref/X_H,ref).
    const double abundance_ratio = hydrogen_scale *
        query.gas_mass_fraction[static_cast<std::size_t>(m_elements[i])] / m_reference_mass_fraction[i];
    result.net_cooling_erg_cm3_s += lambda_ref * nh * nh * electron_correction * std::max(abundance_ratio, 0.0);
  }
  return result;
}

UvBackgroundCompatibility MetalCoolingTable::checkUvBackgroundCompatibility(const std::string& chui_uvb_sha256,
                                                                            bool accept_mismatch) const {
  UvBackgroundCompatibility result;
  const std::string& declared = m_provenance.uv_background;
  const std::string& declared_sha = m_provenance.uv_background_treecool_sha256;
  if (declared == "none") {
    result.compatible = chui_uvb_sha256.empty();
  } else if (!declared_sha.empty() && !chui_uvb_sha256.empty()) {
    result.compatible = declared_sha == chui_uvb_sha256;
  }
  if (result.compatible) {
    result.diagnostic = "metal table uv_background '" + declared + "' matches the configured UV background";
    return result;
  }
  result.diagnostic = "metal cooling table '" + m_provenance.source_path + "' embeds uv_background '" + declared +
      "'" + (declared_sha.empty() ? std::string() : " (TREECOOL sha256 " + declared_sha + ")") +
      " but CHUI's H/He photo-rates use " +
      (chui_uvb_sha256.empty() ? std::string("no UV background")
                               : "the TREECOOL file with sha256 " + chui_uvb_sha256) +
      ": metal photoheating and ionization balance are inconsistent with the H/He network";
  if (!accept_mismatch) {
    throw SourceModelUnavailableError(result.diagnostic +
                                      ". Use a table computed with the same UVB, or set "
                                      "physics.metal_cooling_accept_uv_background_mismatch = true to accept it "
                                      "explicitly (recorded in the run provenance)");
  }
  result.accepted_mismatch = true;
  return result;
}

}  // namespace cosmosim::physics
