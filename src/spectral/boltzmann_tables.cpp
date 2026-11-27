// CHUI-EXT-PHYS-REF-v1 S02: Boltzmann transfer-table adapters (source only,
// NOT compiled or executed in pass 01).

#include "cosmosim/spectral/boltzmann_tables.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>

#include "core/internal/sha256.hpp"

namespace cosmosim::spectral {
namespace {

[[noreturn]] void fail(const std::string& path, const std::string& message) {
  throw std::runtime_error("Boltzmann table '" + path + "': " + message);
}

[[nodiscard]] std::string trim(const std::string& raw) {
  std::size_t b = 0;
  std::size_t e = raw.size();
  while (b < e && std::isspace(static_cast<unsigned char>(raw[b]))) ++b;
  while (e > b && std::isspace(static_cast<unsigned char>(raw[e - 1]))) --e;
  return raw.substr(b, e - b);
}

[[nodiscard]] std::string lower(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

struct RawTable {
  std::vector<std::string> header_lines;  // without the leading '#'
  std::vector<std::vector<double>> rows;
};

[[nodiscard]] RawTable splitTable(const std::string& bytes, const std::string& path) {
  RawTable table;
  std::istringstream stream(bytes);
  std::string line;
  std::size_t width = 0;
  while (std::getline(stream, line)) {
    const std::string text = trim(line);
    if (text.empty()) continue;
    if (text.front() == '#') {
      table.header_lines.push_back(text.substr(1));
      continue;
    }
    std::istringstream row_stream(text);
    std::vector<double> row;
    std::string token;
    while (row_stream >> token) {
      char* end = nullptr;
      const double value = std::strtod(token.c_str(), &end);
      if (end != token.c_str() + static_cast<std::ptrdiff_t>(token.size()) || !std::isfinite(value)) {
        fail(path, "non-numeric or non-finite data token '" + token + "'");
      }
      row.push_back(value);
    }
    if (width == 0U) width = row.size();
    if (row.size() != width) fail(path, "ragged data row (column count changes)");
    table.rows.push_back(std::move(row));
  }
  if (table.rows.size() < 4U) fail(path, "fewer than 4 data rows (k interpolation needs a resolved table)");
  return table;
}

[[nodiscard]] std::map<std::string, std::string> metadata(const RawTable& table) {
  std::map<std::string, std::string> out;
  for (const std::string& line : table.header_lines) {
    const auto colon = line.find(':');
    if (colon == std::string::npos) continue;
    const std::string key = lower(trim(line.substr(0, colon)));
    if (key.empty() || key.find(' ') != std::string::npos) continue;
    out[key] = trim(line.substr(colon + 1U));
  }
  return out;
}

[[nodiscard]] std::vector<std::string> splitWords(const std::string& text, char separator = ' ') {
  std::vector<std::string> out;
  std::string current;
  for (const char c : text) {
    const bool sep = separator == ' ' ? std::isspace(static_cast<unsigned char>(c)) != 0 : c == separator;
    if (sep) {
      if (!trim(current).empty()) out.push_back(trim(current));
      current.clear();
    } else {
      current.push_back(c);
    }
  }
  if (!trim(current).empty()) out.push_back(trim(current));
  return out;
}

// CLASS header: columns announced as "1:k (h/Mpc)   2:d_g   3:d_b ...".
[[nodiscard]] std::vector<std::string> classColumnNames(const RawTable& table, const std::string& path) {
  for (const std::string& line : table.header_lines) {
    if (line.find("1:") == std::string::npos) continue;
    std::vector<std::pair<std::size_t, std::size_t>> starts;  // (column index, position after ':')
    for (std::size_t pos = 0; pos < line.size(); ++pos) {
      if (!std::isdigit(static_cast<unsigned char>(line[pos]))) continue;
      if (pos > 0 && std::isdigit(static_cast<unsigned char>(line[pos - 1]))) continue;
      std::size_t end = pos;
      while (end < line.size() && std::isdigit(static_cast<unsigned char>(line[end]))) ++end;
      if (end < line.size() && line[end] == ':' && (pos == 0 || std::isspace(static_cast<unsigned char>(line[pos - 1])))) {
        starts.emplace_back(static_cast<std::size_t>(std::stoul(line.substr(pos, end - pos))), end + 1U);
      }
    }
    if (starts.empty() || starts.front().first != 1U) continue;
    std::vector<std::string> names;
    for (std::size_t i = 0; i < starts.size(); ++i) {
      if (starts[i].first != i + 1U) fail(path, "CLASS column header numbering is not consecutive");
      std::size_t stop = line.size();
      if (i + 1U < starts.size()) {
        stop = starts[i + 1U].second - 1U;
        while (stop > 0 && std::isdigit(static_cast<unsigned char>(line[stop - 1]))) --stop;
      }
      names.push_back(trim(line.substr(starts[i].second, stop - starts[i].second)));
    }
    return names;
  }
  fail(path, "no CLASS 'N:name' column header found");
}

[[nodiscard]] bool parseHeaderRedshift(const RawTable& table, double& redshift) {
  for (const std::string& line : table.header_lines) {
    const auto pos = line.find("z=");
    if (pos == std::string::npos) continue;
    const std::string tail = line.substr(pos + 2U);
    char* end = nullptr;
    const double z = std::strtod(tail.c_str(), &end);
    if (end != tail.c_str() && std::isfinite(z)) {
      redshift = z;
      return true;
    }
  }
  return false;
}

void checkRedshift(BoltzmannTransferEpoch& epoch, bool from_header, double header_z, double declared_z) {
  epoch.redshift = declared_z;
  epoch.redshift_from_header = from_header;
  if (from_header && std::abs(header_z - declared_z) > 1.0e-6 * (1.0 + declared_z)) {
    fail(epoch.path, "header redshift " + std::to_string(header_z) + " differs from declared neutrinos.table_redshifts "
                     "value " + std::to_string(declared_z));
  }
}

[[nodiscard]] std::size_t requireColumn(const std::vector<std::string>& names, const std::string& wanted,
                                        const std::string& path) {
  const auto it = std::find(names.begin(), names.end(), wanted);
  if (it == names.end()) fail(path, "required column '" + wanted + "' is missing");
  return static_cast<std::size_t>(std::distance(names.begin(), it));
}

[[nodiscard]] long findColumn(const std::vector<std::string>& names, const std::string& wanted) {
  const auto it = std::find(names.begin(), names.end(), wanted);
  return it == names.end() ? -1L : static_cast<long>(std::distance(names.begin(), it));
}

void finalizeEpoch(BoltzmannTransferEpoch& epoch) {
  for (std::size_t i = 0; i < epoch.k_h_mpc.size(); ++i) {
    if (!(epoch.k_h_mpc[i] > 0.0)) fail(epoch.path, "k must be > 0");
    if (i > 0 && !(epoch.k_h_mpc[i] > epoch.k_h_mpc[i - 1])) fail(epoch.path, "k is not strictly increasing");
    if (epoch.delta_cb[i] == 0.0) fail(epoch.path, "delta_cb vanishes at a tabulated k (ratio undefined)");
  }
}

[[nodiscard]] double coldDarkMatterWeight(const core::LambdaCdmBackground& background, double omega_baryon,
                                          const std::string& path) {
  const double omega_cdm = background.coldMatterDensityParameter0() - omega_baryon;
  if (!(omega_cdm >= 0.0) || !(omega_baryon >= 0.0)) fail(path, "inconsistent Omega_b / Omega_cb weights");
  return omega_cdm;
}

void parseChuiV1(const RawTable& table, const BoltzmannLoadRequest& request,
                 const core::LambdaCdmBackground& background, BoltzmannTransferEpoch& epoch) {
  const auto meta = metadata(table);
  const auto get = [&](const std::string& key) {
    const auto it = meta.find(key);
    if (it == meta.end()) fail(request.path, "chui_v1 metadata '" + key + "' is missing");
    return it->second;
  };
  if (get("schema") != "chui_boltzmann_transfer_v1") fail(request.path, "schema is not chui_boltzmann_transfer_v1");
  if (lower(get("k_unit")) != "h/mpc") fail(request.path, "k_unit must be h/Mpc");
  const std::string gauge = lower(get("gauge"));
  const core::BoltzmannGaugeConfig header_gauge =
      gauge == "newtonian" ? core::BoltzmannGaugeConfig::kNewtonian
      : gauge == "synchronous" ? core::BoltzmannGaugeConfig::kSynchronous
                               : (fail(request.path, "unknown gauge '" + gauge + "'"), core::BoltzmannGaugeConfig{});
  if (header_gauge != request.declared_gauge) fail(request.path, "header gauge differs from neutrinos.table_gauge");
  epoch.gauge = header_gauge;
  epoch.source_code = get("source_code");
  epoch.source_version = get("source_version");
  char* end = nullptr;
  const std::string z_text = get("redshift");
  const double header_z = std::strtod(z_text.c_str(), &end);
  if (end == z_text.c_str() || !std::isfinite(header_z)) fail(request.path, "unparseable redshift metadata");
  checkRedshift(epoch, true, header_z, request.declared_redshift);
  // Eigenstate masses must match the background authority.
  const std::vector<std::string> masses = splitWords(get("ncdm_masses_ev"), ',');
  if (masses.size() != background.massiveNeutrinoCount()) {
    fail(request.path, "ncdm_masses_ev lists " + std::to_string(masses.size()) + " eigenstates, background has " +
                           std::to_string(background.massiveNeutrinoCount()));
  }
  for (std::size_t i = 0; i < masses.size(); ++i) {
    const double m = std::stod(masses[i]);
    const double configured = background.config().radiation_species.ncdm_masses_ev.at(i);
    if (std::abs(m - configured) > 1.0e-6 * configured) {
      fail(request.path, "ncdm mass " + masses[i] + " eV differs from the configured eigenstate mass");
    }
  }
  const std::vector<std::string> names = splitWords(get("columns"));
  if (!table.rows.empty() && names.size() != table.rows.front().size()) {
    fail(request.path, "columns metadata does not match the data width");
  }
  const std::size_t ik = requireColumn(names, "k", request.path);
  const std::size_t icb = requireColumn(names, "delta_cb", request.path);
  const std::size_t inu = requireColumn(names, "delta_ncdm", request.path);
  const long itcb = findColumn(names, "theta_cb");
  const long itnu = findColumn(names, "theta_ncdm");
  if ((itcb < 0) != (itnu < 0)) fail(request.path, "theta_cb and theta_ncdm must be provided together");
  if (itcb >= 0 && header_gauge != core::BoltzmannGaugeConfig::kNewtonian) {
    fail(request.path, "velocity columns are accepted only in Newtonian gauge (synchronous theta is frame-defined)");
  }
  for (const auto& row : table.rows) {
    epoch.k_h_mpc.push_back(row[ik]);
    epoch.delta_cb.push_back(row[icb]);
    epoch.delta_ncdm.push_back(row[inu]);
    if (itcb >= 0) {
      epoch.theta_cb.push_back(row[static_cast<std::size_t>(itcb)]);
      epoch.theta_ncdm.push_back(row[static_cast<std::size_t>(itnu)]);
    }
  }
  (void)coldDarkMatterWeight(background, request.omega_baryon, request.path);
}

void parseClassTk(const RawTable& table, const BoltzmannLoadRequest& request,
                  const core::LambdaCdmBackground& background, BoltzmannTransferEpoch& epoch) {
  const std::vector<std::string> names = classColumnNames(table, request.path);
  if (names.size() != table.rows.front().size()) fail(request.path, "CLASS header width differs from data width");
  if (names.front().rfind("k", 0) != 0 || names.front().find("h/Mpc") == std::string::npos) {
    fail(request.path, "first CLASS column must be 'k (h/Mpc)'");
  }
  epoch.source_code = "CLASS";
  epoch.source_version = "undeclared";
  for (const std::string& line : table.header_lines) {
    const auto pos = line.find("CLASS v");
    if (pos != std::string::npos) epoch.source_version = trim(line.substr(pos + 6U));
  }
  epoch.gauge = request.declared_gauge;
  double header_z = 0.0;
  const bool has_z = parseHeaderRedshift(table, header_z);
  checkRedshift(epoch, has_z, header_z, request.declared_redshift);
  const std::size_t icdm = requireColumn(names, "d_cdm", request.path);
  const std::size_t ib = requireColumn(names, "d_b", request.path);
  const std::size_t n_ncdm = background.massiveNeutrinoCount();
  std::vector<std::size_t> inu;
  for (std::size_t i = 0; i < n_ncdm; ++i) {
    inu.push_back(requireColumn(names, "d_ncdm[" + std::to_string(i) + "]", request.path));
  }
  if (findColumn(names, "d_ncdm[" + std::to_string(n_ncdm) + "]") >= 0) {
    fail(request.path, "file has more ncdm species than the configured background eigenstates");
  }
  const long itcdm = findColumn(names, "t_cdm");
  if (itcdm >= 0 && request.declared_gauge == core::BoltzmannGaugeConfig::kSynchronous) {
    fail(request.path, "t_cdm present but synchronous gauge declared (theta_cdm = 0 by construction in that gauge)");
  }
  const long itb = findColumn(names, "t_b");
  std::vector<long> itnu;
  bool velocities = itcdm >= 0 && itb >= 0;
  for (std::size_t i = 0; i < n_ncdm; ++i) {
    itnu.push_back(findColumn(names, "t_ncdm[" + std::to_string(i) + "]"));
    velocities = velocities && itnu.back() >= 0;
  }
  const double omega_cdm = coldDarkMatterWeight(background, request.omega_baryon, request.path);
  const double omega_b = request.omega_baryon;
  const double a = 1.0 / (1.0 + request.declared_redshift);
  std::vector<double> rho(n_ncdm);
  std::vector<double> rho_plus_p(n_ncdm);
  double rho_sum = 0.0;
  double rho_plus_p_sum = 0.0;
  for (std::size_t i = 0; i < n_ncdm; ++i) {
    rho[i] = background.massiveNeutrinoDensityParameter(i, a);
    rho_plus_p[i] = rho[i] + background.massiveNeutrinoPressureParameter(i, a);
    rho_sum += rho[i];
    rho_plus_p_sum += rho_plus_p[i];
  }
  for (const auto& row : table.rows) {
    epoch.k_h_mpc.push_back(row[0]);
    epoch.delta_cb.push_back((omega_cdm * row[icdm] + omega_b * row[ib]) / (omega_cdm + omega_b));
    double nu = 0.0;
    for (std::size_t i = 0; i < n_ncdm; ++i) nu += rho[i] * row[inu[i]];
    epoch.delta_ncdm.push_back(nu / rho_sum);
    if (velocities) {
      epoch.theta_cb.push_back((omega_cdm * row[static_cast<std::size_t>(itcdm)] +
                                omega_b * row[static_cast<std::size_t>(itb)]) / (omega_cdm + omega_b));
      double tnu = 0.0;
      for (std::size_t i = 0; i < n_ncdm; ++i) tnu += rho_plus_p[i] * row[static_cast<std::size_t>(itnu[i])];
      epoch.theta_ncdm.push_back(tnu / rho_plus_p_sum);
    }
  }
}

void parseCambTransfer(const RawTable& table, const BoltzmannLoadRequest& request, BoltzmannTransferEpoch& epoch) {
  if (table.rows.front().size() != 13U) {
    fail(request.path, "camb_transfer expects the 13-column CAMB >= 1.0 layout, found " +
                           std::to_string(table.rows.front().size()));
  }
  if (request.declared_gauge != core::BoltzmannGaugeConfig::kSynchronous) {
    fail(request.path, "CAMB transfer output is synchronous gauge (CDM frame); declare table_gauge = synchronous");
  }
  for (const std::string& line : table.header_lines) {
    const std::vector<std::string> words = splitWords(line);
    if (words.size() == 13U && lower(words[0]) == "k/h" && (words[5] != "mass_nu" || words[7] != "no_nu")) {
      fail(request.path, "CAMB header column names do not match the expected layout (mass_nu at 6, no_nu at 8)");
    }
  }
  epoch.source_code = "CAMB";
  epoch.source_version = "undeclared";
  epoch.gauge = core::BoltzmannGaugeConfig::kSynchronous;
  checkRedshift(epoch, false, 0.0, request.declared_redshift);
  for (const auto& row : table.rows) {
    epoch.k_h_mpc.push_back(row[0]);
    epoch.delta_cb.push_back(row[7]);    // no_nu: CDM + baryon (massive nu excluded)
    epoch.delta_ncdm.push_back(row[5]);  // mass_nu
  }
}

[[nodiscard]] std::size_t intervalIndex(const std::vector<double>& grid, double x, const std::string& what) {
  if (!(x >= grid.front()) || !(x <= grid.back())) {
    throw std::out_of_range(what + " " + std::to_string(x) + " outside the tabulated domain [" +
                            std::to_string(grid.front()) + ", " + std::to_string(grid.back()) +
                            "] (no extrapolation)");
  }
  const auto it = std::upper_bound(grid.begin(), grid.end(), x);
  std::size_t hi = static_cast<std::size_t>(std::distance(grid.begin(), it));
  hi = std::clamp<std::size_t>(hi, 1U, grid.size() - 1U);
  return hi - 1U;
}

}  // namespace

BoltzmannTransferEpoch parseBoltzmannTransferEpoch(const std::string& bytes, const BoltzmannLoadRequest& request,
                                                   const core::LambdaCdmBackground& background) {
  if (background.massiveNeutrinoCount() == 0U) {
    fail(request.path, "the background has no massive eigenstate; a neutrino transfer table is meaningless");
  }
  const RawTable table = splitTable(bytes, request.path);
  BoltzmannTransferEpoch epoch;
  epoch.path = request.path;
  epoch.format = core::boltzmannTableFormatConfigToString(request.format);
  switch (request.format) {
    case core::BoltzmannTableFormatConfig::kChuiV1:
      parseChuiV1(table, request, background, epoch);
      break;
    case core::BoltzmannTableFormatConfig::kClassTk:
      parseClassTk(table, request, background, epoch);
      break;
    case core::BoltzmannTableFormatConfig::kCambTransfer:
      parseCambTransfer(table, request, epoch);
      break;
  }
  finalizeEpoch(epoch);
  return epoch;
}

BoltzmannTransferEpoch loadBoltzmannTransferEpoch(const BoltzmannLoadRequest& request,
                                                  const core::LambdaCdmBackground& background) {
  if (request.path.empty()) fail(request.path, "empty path; no Boltzmann table is bundled");
  std::ifstream input(request.path, std::ios::binary);
  if (!input) fail(request.path, "cannot open file");
  const std::string bytes((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  const std::string sha = core::internal::sha256Hex(bytes);
  if (lower(sha) != lower(request.expected_sha256)) {
    fail(request.path, "SHA-256 mismatch (expected " + request.expected_sha256 + ", got " + sha + ")");
  }
  BoltzmannTransferEpoch epoch = parseBoltzmannTransferEpoch(bytes, request, background);
  epoch.sha256 = sha;
  return epoch;
}

BoltzmannTransferSet::BoltzmannTransferSet(std::vector<BoltzmannTransferEpoch> epochs) : m_epochs(std::move(epochs)) {
  if (m_epochs.size() < 2U) throw std::invalid_argument("BoltzmannTransferSet needs >= 2 epochs (no a extrapolation)");
  m_k_min = m_epochs.front().kMin();
  m_k_max = m_epochs.front().kMax();
  for (std::size_t i = 0; i < m_epochs.size(); ++i) {
    if (m_epochs[i].gauge != m_epochs.front().gauge) throw std::invalid_argument("epochs mix gauges");
    if (i > 0 && !(m_epochs[i].redshift < m_epochs[i - 1].redshift)) {
      throw std::invalid_argument("epochs must be ordered by strictly decreasing redshift");
    }
    m_k_min = std::max(m_k_min, m_epochs[i].kMin());
    m_k_max = std::min(m_k_max, m_epochs[i].kMax());
  }
  if (!(m_k_max > m_k_min)) throw std::invalid_argument("epochs share no common k range");
}

double BoltzmannTransferSet::scaleFactorMin() const { return 1.0 / (1.0 + m_epochs.front().redshift); }
double BoltzmannTransferSet::scaleFactorMax() const { return 1.0 / (1.0 + m_epochs.back().redshift); }

double BoltzmannTransferSet::epochValue(const BoltzmannTransferEpoch& epoch, Quantity quantity, double k) const {
  if ((quantity != Quantity::kNcdmRatio) && !epoch.hasVelocities()) {
    throw std::runtime_error("Boltzmann epoch '" + epoch.path + "' carries no velocity divergence columns");
  }
  const auto value_at = [&](std::size_t i) {
    switch (quantity) {
      case Quantity::kNcdmRatio:
        return epoch.delta_ncdm[i] / epoch.delta_cb[i];
      case Quantity::kThetaNcdm:
        return epoch.theta_ncdm[i] / epoch.delta_cb[i];
      case Quantity::kThetaCold:
        return epoch.theta_cb[i] / epoch.delta_cb[i];
    }
    return 0.0;
  };
  const std::size_t lo = intervalIndex(epoch.k_h_mpc, k, "k [h/Mpc]");
  const double t = std::log(k / epoch.k_h_mpc[lo]) / std::log(epoch.k_h_mpc[lo + 1U] / epoch.k_h_mpc[lo]);
  return (1.0 - t) * value_at(lo) + t * value_at(lo + 1U);
}

double BoltzmannTransferSet::interpolate(Quantity quantity, double k, double scale_factor) const {
  if (m_epochs.empty()) throw std::logic_error("empty BoltzmannTransferSet");
  if (!(k >= m_k_min) || !(k <= m_k_max)) {
    throw std::out_of_range("k = " + std::to_string(k) + " h/Mpc outside the common table range (no extrapolation)");
  }
  std::vector<double> ln_a(m_epochs.size());
  for (std::size_t i = 0; i < m_epochs.size(); ++i) ln_a[i] = -std::log1p(m_epochs[i].redshift);
  const double x = std::log(scale_factor);
  const std::size_t lo = intervalIndex(ln_a, x, "ln a");
  const double t = (x - ln_a[lo]) / (ln_a[lo + 1U] - ln_a[lo]);
  return (1.0 - t) * epochValue(m_epochs[lo], quantity, k) + t * epochValue(m_epochs[lo + 1U], quantity, k);
}

double BoltzmannTransferSet::ncdmToColdRatio(double k, double a) const { return interpolate(Quantity::kNcdmRatio, k, a); }
double BoltzmannTransferSet::thetaNcdmOverDeltaCold(double k, double a) const {
  return interpolate(Quantity::kThetaNcdm, k, a);
}
double BoltzmannTransferSet::thetaColdOverDeltaCold(double k, double a) const {
  return interpolate(Quantity::kThetaCold, k, a);
}

std::string BoltzmannTransferSet::provenanceDigest() const {
  std::ostringstream text;
  text.precision(17);
  for (const auto& epoch : m_epochs) {
    text << epoch.format << '|' << epoch.source_code << '|' << epoch.source_version << '|'
         << core::boltzmannGaugeConfigToString(epoch.gauge) << '|' << epoch.redshift << '|' << epoch.sha256 << '\n';
  }
  return core::internal::sha256Hex(text.str());
}

}  // namespace cosmosim::spectral
