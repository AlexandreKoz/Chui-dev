#include "cosmosim/physics/metal_ion_network.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cosmosim::physics {
namespace {

constexpr std::size_t k_process_ionization = 0U;
constexpr std::size_t k_process_recombination = 1U;
constexpr std::size_t k_process_cooling = 2U;
constexpr std::size_t k_process_photoionization = 3U;
constexpr std::size_t k_process_photoheating = 4U;
constexpr std::size_t k_process_count = 5U;
constexpr std::size_t k_max_stages = 27U;  // Fe
// CODATA 2018 exact unit conversion (definition of the electron volt), not a
// fitted physical coefficient.
constexpr double k_electron_volt_erg = 1.602176634e-12;

constexpr std::array<std::string_view, k_metal_ion_element_count> k_symbols{"C", "O", "Si", "Fe"};

[[nodiscard]] std::string roman(std::size_t value) {
  static constexpr std::array<std::pair<std::size_t, std::string_view>, 9> k_table{{
      {50U, "L"}, {40U, "XL"}, {10U, "X"}, {9U, "IX"}, {5U, "V"}, {4U, "IV"}, {1U, "I"}, {0U, ""}, {0U, ""}}};
  std::string out;
  for (const auto& [number, text] : k_table) {
    if (number == 0U) break;
    while (value >= number) {
      out += text;
      value -= number;
    }
  }
  return out;
}

[[nodiscard]] MetalIonElement elementFromSymbol(std::string_view symbol) {
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    if (k_symbols[e] == symbol) return static_cast<MetalIonElement>(e);
  }
  throw std::invalid_argument("metal-ion element '" + std::string(symbol) +
                              "' is not one of C, O, Si, Fe");
}

[[nodiscard]] double parseLog10(const std::string& token, const std::string& context) {
  if (token == "-inf") return -std::numeric_limits<double>::infinity();
  std::size_t consumed = 0U;
  double value = 0.0;
  try {
    value = std::stod(token, &consumed);
  } catch (const std::exception&) {
    throw std::runtime_error(context + ": invalid number '" + token + "'");
  }
  if (consumed != token.size() || !std::isfinite(value)) {
    throw std::runtime_error(context + ": invalid number '" + token + "'");
  }
  return value;
}

}  // namespace

ChemicalElement chemicalElementOf(MetalIonElement element) noexcept {
  switch (element) {
    case MetalIonElement::kC: return ChemicalElement::kC;
    case MetalIonElement::kO: return ChemicalElement::kO;
    case MetalIonElement::kSi: return ChemicalElement::kSi;
    case MetalIonElement::kFe: return ChemicalElement::kFe;
  }
  return ChemicalElement::kC;
}

std::string_view metalIonElementSymbol(MetalIonElement element) noexcept {
  return k_symbols[static_cast<std::size_t>(element)];
}

std::string metalIonStageName(MetalIonElement element, std::size_t stage) {
  if (stage >= metalIonStageCount(element)) {
    throw std::out_of_range("metal-ion stage exceeds Z for its element");
  }
  return std::string(metalIonElementSymbol(element)) + "_" + roman(stage + 1U);
}

MetalIonElementSet MetalIonElementSet::fromMask(std::uint8_t mask) {
  if ((mask & ~static_cast<std::uint8_t>(0x0FU)) != 0U) {
    throw std::invalid_argument("metal-ion element mask has unknown bits");
  }
  MetalIonElementSet set;
  set.m_mask = mask;
  return set;
}

MetalIonElementSet MetalIonElementSet::parse(std::string_view text) {
  if (text == "none" || text.empty()) return {};
  std::uint8_t mask = 0U;
  std::size_t begin = 0U;
  while (begin <= text.size()) {
    const std::size_t comma = text.find(',', begin);
    const std::string_view token =
        text.substr(begin, comma == std::string_view::npos ? std::string_view::npos : comma - begin);
    const MetalIonElement element = elementFromSymbol(token);
    const auto bit = static_cast<std::uint8_t>(1U << static_cast<unsigned>(element));
    if ((mask & bit) != 0U) {
      throw std::invalid_argument("metal-ion element '" + std::string(token) + "' listed twice");
    }
    mask = static_cast<std::uint8_t>(mask | bit);
    if (comma == std::string_view::npos) break;
    begin = comma + 1U;
  }
  return fromMask(mask);
}

MetalIonElementSet MetalIonElementSet::fromStageCount(std::size_t stage_count) {
  for (std::uint8_t mask = 0U; mask < 16U; ++mask) {
    if (fromMask(mask).stageCount() == stage_count) return fromMask(mask);
  }
  throw std::runtime_error("no metal-ion element subset has " + std::to_string(stage_count) + " stages");
}

bool MetalIonElementSet::containsChemicalElement(ChemicalElement element) const noexcept {
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    const auto ion_element = static_cast<MetalIonElement>(e);
    if (contains(ion_element) && chemicalElementOf(ion_element) == element) return true;
  }
  return false;
}

std::size_t MetalIonElementSet::stageCount() const noexcept {
  std::size_t count = 0U;
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    if (contains(static_cast<MetalIonElement>(e))) count += metalIonStageCount(static_cast<MetalIonElement>(e));
  }
  return count;
}

std::size_t MetalIonElementSet::stageOffset(MetalIonElement element) const {
  if (!contains(element)) throw std::out_of_range("metal-ion element is not part of the selected set");
  std::size_t offset = 0U;
  for (std::size_t e = 0; e < static_cast<std::size_t>(element); ++e) {
    if (contains(static_cast<MetalIonElement>(e))) offset += metalIonStageCount(static_cast<MetalIonElement>(e));
  }
  return offset;
}

std::string MetalIonElementSet::label() const {
  if (empty()) return "none";
  std::string out;
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    if (!contains(static_cast<MetalIonElement>(e))) continue;
    if (!out.empty()) out += ",";
    out += k_symbols[e];
  }
  return out;
}

// ---------------------------------------------------------------------------
// Dataset
// ---------------------------------------------------------------------------

MetalIonRateTable MetalIonRateTable::loadFromFile(const std::string& path, MetalIonElementSet required) {
  if (path.empty()) {
    throw std::runtime_error(
        "fixed_detailed metal-ion chemistry requires physics.metal_ion_rate_table_path "
        "(schema chui_metal_ion_rates_v1); no rate coefficients are compiled into CHUI");
  }
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open metal-ion rate table '" + path + "'");
  std::ostringstream text;
  text << input.rdbuf();
  return parse(text.str(), path, required);
}

std::size_t MetalIonRateTable::seriesIndex(MetalIonElement element, std::size_t stage, std::size_t process) const {
  return (static_cast<std::size_t>(element) * k_max_stages + stage) * k_process_count + process;
}

MetalIonRateTable MetalIonRateTable::parse(std::string_view text, std::string_view source_label,
                                           MetalIonElementSet required) {
  if (required.empty()) throw std::invalid_argument("metal-ion rate table requested for an empty element set");
  MetalIonRateTable table;
  table.m_source_label = std::string(source_label);
  table.m_series.assign(k_metal_ion_element_count * k_max_stages * k_process_count, Series{});
  table.m_ionization_potential_erg.assign(k_metal_ion_element_count * k_max_stages,
                                          std::numeric_limits<double>::quiet_NaN());
  const std::string where = "metal-ion rate table '" + table.m_source_label + "'";
  bool schema_seen = false;
  bool verification_seen = false;
  std::istringstream lines{std::string(text)};
  std::string line;
  std::size_t line_number = 0U;
  while (std::getline(lines, line)) {
    ++line_number;
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos) line.resize(hash);
    std::istringstream tokens(line);
    std::string keyword;
    if (!(tokens >> keyword)) continue;
    const std::string context = where + " line " + std::to_string(line_number);
    const auto rest_of_line = [&]() {
      std::string rest;
      std::getline(tokens, rest);
      const auto first = rest.find_first_not_of(" \t");
      return first == std::string::npos ? std::string{} : rest.substr(first);
    };
    if (keyword == "schema") {
      std::string schema;
      tokens >> schema;
      if (schema != "chui_metal_ion_rates_v1") throw std::runtime_error(context + ": unsupported schema '" + schema + "'");
      schema_seen = true;
    } else if (keyword == "provenance") {
      table.m_provenance = rest_of_line();
    } else if (keyword == "license") {
      table.m_license = rest_of_line();
    } else if (keyword == "verification") {
      std::string value;
      tokens >> value;
      if (value == "verified_against_source") {
        table.m_declared_verified = true;
      } else if (value != "unverified") {
        throw std::runtime_error(context + ": verification must be verified_against_source or unverified");
      }
      verification_seen = true;
    } else if (keyword == "energy_convention") {
      std::string value;
      tokens >> value;
      if (value != "collisional_ionization_work_in_cooling") {
        throw std::runtime_error(context + ": unsupported energy_convention '" + value +
                                 "' (supported: collisional_ionization_work_in_cooling)");
      }
      if (!table.m_energy_convention.empty()) throw std::runtime_error(context + ": energy_convention declared twice");
      table.m_energy_convention = value;
    } else if (keyword == "ionization_potential") {
      std::string symbol, token;
      std::size_t stage = 0U;
      if (!(tokens >> symbol >> stage >> token)) throw std::runtime_error(context + ": malformed ionization_potential");
      const MetalIonElement element = elementFromSymbol(symbol);
      if (stage + 1U >= metalIonStageCount(element)) {
        throw std::runtime_error(context + ": ionization_potential stage must be < Z (bare nucleus cannot ionize)");
      }
      std::size_t consumed = 0U;
      double ev = 0.0;
      try {
        ev = std::stod(token, &consumed);
      } catch (const std::exception&) {
        throw std::runtime_error(context + ": invalid ionization potential '" + token + "'");
      }
      if (consumed != token.size() || !std::isfinite(ev) || !(ev > 0.0)) {
        throw std::runtime_error(context + ": ionization potential must be a positive finite eV value");
      }
      double& slot = table.m_ionization_potential_erg[static_cast<std::size_t>(element) * k_max_stages + stage];
      if (!std::isnan(slot)) throw std::runtime_error(context + ": duplicate ionization_potential");
      slot = ev * k_electron_volt_erg;
      if (tokens >> token) throw std::runtime_error(context + ": trailing tokens after ionization_potential");
    } else if (keyword == "temperature_grid_log10_k" || keyword == "redshift_grid") {
      std::size_t count = 0U;
      if (!(tokens >> count) || count < 2U) throw std::runtime_error(context + ": grid needs a count >= 2");
      std::vector<double>& grid = keyword == "redshift_grid" ? table.m_redshift : table.m_log_t;
      if (!grid.empty()) throw std::runtime_error(context + ": grid declared twice");
      grid.reserve(count);
      std::string token;
      for (std::size_t i = 0; i < count; ++i) {
        if (!(tokens >> token)) throw std::runtime_error(context + ": grid has fewer values than declared");
        const double value = parseLog10(token, context);
        if (!std::isfinite(value) || (!grid.empty() && !(value > grid.back()))) {
          throw std::runtime_error(context + ": grid must be finite and strictly increasing");
        }
        grid.push_back(value);
      }
      if (tokens >> token) throw std::runtime_error(context + ": grid has more values than declared");
    } else if (keyword == "rate" || keyword == "photo") {
      std::string symbol, process;
      std::size_t stage = 0U;
      if (!(tokens >> symbol >> stage >> process)) throw std::runtime_error(context + ": malformed rate line");
      const MetalIonElement element = elementFromSymbol(symbol);
      if (stage >= metalIonStageCount(element)) throw std::runtime_error(context + ": stage exceeds Z");
      std::size_t process_index = 0U;
      const std::vector<double>* grid = nullptr;
      if (keyword == "rate") {
        grid = &table.m_log_t;
        if (process == "collisional_ionization") process_index = k_process_ionization;
        else if (process == "recombination") process_index = k_process_recombination;
        else if (process == "cooling") process_index = k_process_cooling;
        else throw std::runtime_error(context + ": unknown temperature process '" + process + "'");
      } else {
        grid = &table.m_redshift;
        if (process == "ionization") process_index = k_process_photoionization;
        else if (process == "heating") process_index = k_process_photoheating;
        else throw std::runtime_error(context + ": unknown photo process '" + process + "'");
      }
      if (grid->empty()) throw std::runtime_error(context + ": rate line precedes its grid declaration");
      Series& series = table.m_series[table.seriesIndex(element, stage, process_index)];
      if (!series.log10_values.empty()) throw std::runtime_error(context + ": duplicate rate line");
      series.log10_values.reserve(grid->size());
      std::string token;
      for (std::size_t i = 0; i < grid->size(); ++i) {
        if (!(tokens >> token)) throw std::runtime_error(context + ": fewer rate values than grid points");
        series.log10_values.push_back(parseLog10(token, context));
      }
      if (tokens >> token) throw std::runtime_error(context + ": more rate values than grid points");
    } else {
      throw std::runtime_error(context + ": unknown keyword '" + keyword + "'");
    }
  }
  if (!schema_seen) throw std::runtime_error(where + ": missing 'schema chui_metal_ion_rates_v1'");
  if (table.m_provenance.empty() || table.m_license.empty()) {
    throw std::runtime_error(where + ": provenance and license are mandatory");
  }
  if (!verification_seen) throw std::runtime_error(where + ": missing verification declaration");
  if (table.m_log_t.empty()) throw std::runtime_error(where + ": missing temperature grid");

  std::uint8_t present = 0U;
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    const auto element = static_cast<MetalIonElement>(e);
    const std::size_t stages = metalIonStageCount(element);
    bool any = false;
    for (std::size_t s = 0; s < stages && !any; ++s) {
      for (std::size_t p = 0; p < k_process_count; ++p) {
        any = any || !table.m_series[table.seriesIndex(element, s, p)].log10_values.empty();
      }
    }
    if (!any) continue;
    present = static_cast<std::uint8_t>(present | (1U << e));
    bool photo_any = false;
    for (std::size_t s = 0; s < stages; ++s) {
      const std::string stage_label = metalIonStageName(element, s);
      for (std::size_t p = k_process_ionization; p <= k_process_cooling; ++p) {
        const Series& series = table.m_series[table.seriesIndex(element, s, p)];
        if (series.log10_values.empty()) {
          throw std::runtime_error(where + ": " + stage_label + " is missing a temperature process");
        }
      }
      const auto zero_everywhere = [](const Series& series) {
        return std::all_of(series.log10_values.begin(), series.log10_values.end(),
                           [](double v) { return v == -std::numeric_limits<double>::infinity(); });
      };
      if (s + 1U == stages && !zero_everywhere(table.m_series[table.seriesIndex(element, s, k_process_ionization)])) {
        throw std::runtime_error(where + ": bare nucleus " + stage_label + " cannot ionize");
      }
      if (s == 0U && !zero_everywhere(table.m_series[table.seriesIndex(element, s, k_process_recombination)])) {
        throw std::runtime_error(where + ": neutral " + stage_label + " cannot recombine");
      }
      photo_any = photo_any ||
          !table.m_series[table.seriesIndex(element, s, k_process_photoionization)].log10_values.empty() ||
          !table.m_series[table.seriesIndex(element, s, k_process_photoheating)].log10_values.empty();
    }
    if (photo_any) {
      for (std::size_t s = 0; s + 1U < stages; ++s) {
        for (std::size_t p = k_process_photoionization; p <= k_process_photoheating; ++p) {
          if (table.m_series[table.seriesIndex(element, s, p)].log10_values.empty()) {
            throw std::runtime_error(where + ": photo rates for " + std::string(k_symbols[e]) +
                                     " must cover every ionizable stage (missing " +
                                     metalIonStageName(element, s) + ")");
          }
        }
      }
      table.m_has_photo[e] = true;
    }
  }
  table.m_elements = MetalIonElementSet::fromMask(present);
  if ((required.mask() & ~present) != 0U) {
    throw std::runtime_error(where + ": selected metal-ion elements '" + required.label() +
                             "' are not all provided (dataset covers '" + table.m_elements.label() + "')");
  }
  return table;
}

bool MetalIonRateTable::hasIonizationEnergyContract(MetalIonElement element) const noexcept {
  if (m_energy_convention.empty() || !m_elements.contains(element)) return false;
  const std::size_t stages = metalIonStageCount(element);
  for (std::size_t k = 0; k + 1U < stages; ++k) {
    if (std::isnan(m_ionization_potential_erg[static_cast<std::size_t>(element) * k_max_stages + k])) return false;
  }
  return true;
}

double MetalIonRateTable::ionizationPotentialErg(MetalIonElement element, std::size_t stage) const {
  if (stage + 1U >= metalIonStageCount(element)) {
    throw std::out_of_range("ionization potential requested for a bare nucleus or out-of-range stage");
  }
  const double value = m_ionization_potential_erg[static_cast<std::size_t>(element) * k_max_stages + stage];
  if (std::isnan(value)) {
    throw std::out_of_range("metal-ion dataset '" + m_source_label + "' supplies no ionization potential for " +
                            metalIonStageName(element, stage));
  }
  return value;
}

double MetalIonRateTable::storedIonizationEnergyPerNucleusErg(MetalIonElement element,
                                                              std::span<const double> fractions) const {
  const std::size_t stages = metalIonStageCount(element);
  if (fractions.size() != stages) throw std::invalid_argument("metal-ion fraction span size mismatch");
  double cumulative = 0.0;
  double stored = 0.0;
  for (std::size_t k = 1; k < stages; ++k) {
    cumulative += ionizationPotentialErg(element, k - 1U);
    stored += fractions[k] * cumulative;
  }
  return stored;
}

double MetalIonRateTable::minimumTemperatureK() const noexcept { return std::pow(10.0, m_log_t.front()); }
double MetalIonRateTable::maximumTemperatureK() const noexcept { return std::pow(10.0, m_log_t.back()); }

double MetalIonRateTable::interpolate(const std::vector<double>& grid, const Series& series, double coordinate,
                                      MetalIonRangePolicy policy, bool& clamped, const char* what) {
  if (series.log10_values.empty()) return 0.0;
  if (!std::isfinite(coordinate)) throw std::range_error(std::string("metal-ion ") + what + " is not finite");
  double x = coordinate;
  if (x < grid.front() || x > grid.back()) {
    if (policy == MetalIonRangePolicy::kReject) {
      throw std::range_error(std::string("metal-ion ") + what + " outside the dataset grid");
    }
    clamped = true;
    x = std::clamp(x, grid.front(), grid.back());
  }
  const auto upper = std::upper_bound(grid.begin(), grid.end(), x);
  const std::size_t hi = std::min<std::size_t>(static_cast<std::size_t>(upper - grid.begin()), grid.size() - 1U);
  const std::size_t lo = hi == 0U ? 0U : hi - 1U;
  const double w = hi == lo ? 0.0 : (x - grid[lo]) / (grid[hi] - grid[lo]);
  const double a = series.log10_values[lo];
  const double b = series.log10_values[hi];
  if (std::isfinite(a) && std::isfinite(b)) return std::pow(10.0, a + w * (b - a));
  // One endpoint is an exact zero: interpolate the linear coefficient.
  const double va = std::isfinite(a) ? std::pow(10.0, a) : 0.0;
  const double vb = std::isfinite(b) ? std::pow(10.0, b) : 0.0;
  return va + w * (vb - va);
}

MetalIonStageRates MetalIonRateTable::evaluate(MetalIonElement element, std::size_t stage, double temperature_k,
                                               double redshift, MetalIonRangePolicy policy, bool& clamped) const {
  if (!m_elements.contains(element) || stage >= metalIonStageCount(element)) {
    throw std::out_of_range("metal-ion rate request outside the loaded dataset");
  }
  if (!(temperature_k > 0.0)) throw std::range_error("metal-ion temperature must be positive");
  const double log_t = std::log10(temperature_k);
  MetalIonStageRates r;
  r.collisional_ionization_cm3_s = interpolate(
      m_log_t, m_series[seriesIndex(element, stage, k_process_ionization)], log_t, policy, clamped, "temperature");
  r.recombination_cm3_s = interpolate(
      m_log_t, m_series[seriesIndex(element, stage, k_process_recombination)], log_t, policy, clamped, "temperature");
  r.cooling_erg_cm3_s = interpolate(
      m_log_t, m_series[seriesIndex(element, stage, k_process_cooling)], log_t, policy, clamped, "temperature");
  if (m_has_photo[static_cast<std::size_t>(element)]) {
    r.photoionization_s = interpolate(m_redshift, m_series[seriesIndex(element, stage, k_process_photoionization)],
                                      redshift, policy, clamped, "redshift");
    r.photoheating_erg_s = interpolate(m_redshift, m_series[seriesIndex(element, stage, k_process_photoheating)],
                                       redshift, policy, clamped, "redshift");
  }
  return r;
}

// ---------------------------------------------------------------------------
// Ladder solver
// ---------------------------------------------------------------------------

MetalIonLadderSolver::MetalIonLadderSolver(const MetalIonRateTable& table, MetalIonRangePolicy policy,
                                           double relative_tolerance, double absolute_tolerance,
                                           std::uint32_t max_substeps)
    : m_table(&table),
      m_policy(policy),
      m_relative_tolerance(relative_tolerance),
      m_absolute_tolerance(absolute_tolerance),
      m_max_substeps(max_substeps) {
  if (!(relative_tolerance > 0.0) || !(absolute_tolerance > 0.0) || max_substeps == 0U) {
    throw std::invalid_argument("metal-ion solver tolerances and substep cap must be positive");
  }
}

void MetalIonLadderSolver::rates(MetalIonElement element, const MetalIonEnvironment& env, std::span<double> ionize,
                                 std::span<double> recombine, MetalIonAdvanceReport& report) const {
  const std::size_t stages = metalIonStageCount(element);
  const double ne = std::max(env.electron_cm3, 0.0);
  const double attenuation = std::clamp(env.photo_attenuation, 0.0, 1.0);
  bool clamped = false;
  for (std::size_t k = 0; k < stages; ++k) {
    const MetalIonStageRates r = m_table->evaluate(element, k, env.temperature_k, env.redshift, m_policy, clamped);
    ionize[k] = k + 1U < stages ? ne * r.collisional_ionization_cm3_s + attenuation * r.photoionization_s : 0.0;
    recombine[k] = k > 0U ? ne * r.recombination_cm3_s : 0.0;
  }
  if (clamped) ++report.range_clamps;
}

void MetalIonLadderSolver::backwardEulerStep(std::span<const double> ionize, std::span<const double> recombine,
                                             double dt, std::span<double> x) {
  // (I - dt A) y = x, A tridiagonal: sub-diagonal I_{k-1}, super R_{k+1},
  // diagonal -(I_k + R_k). Thomas algorithm (diagonally dominant M-matrix).
  const std::size_t n = x.size();
  std::array<double, k_max_stages> c_prime{};
  std::array<double, k_max_stages> d_prime{};
  for (std::size_t k = 0; k < n; ++k) {
    const double diag = 1.0 + dt * (ionize[k] + recombine[k]);
    const double lower = k > 0U ? -dt * ionize[k - 1U] : 0.0;
    const double upper = k + 1U < n ? -dt * recombine[k + 1U] : 0.0;
    const double denom = diag - (k > 0U ? lower * c_prime[k - 1U] : 0.0);
    c_prime[k] = upper / denom;
    d_prime[k] = (x[k] - (k > 0U ? lower * d_prime[k - 1U] : 0.0)) / denom;
  }
  for (std::size_t k = n; k-- > 0U;) {
    x[k] = d_prime[k] - (k + 1U < n ? c_prime[k] * x[k + 1U] : 0.0);
  }
  // The exact solution is non-negative (M-matrix); remove round-off only.
  for (std::size_t k = 0; k < n; ++k) x[k] = std::max(x[k], 0.0);
}

SourceAdvanceStatus MetalIonLadderSolver::advance(MetalIonElement element, std::span<double> fractions,
                                                  const MetalIonEnvironment& env, double dt_seconds,
                                                  MetalIonAdvanceReport& report) const {
  const std::size_t stages = metalIonStageCount(element);
  if (fractions.size() != stages) throw std::invalid_argument("metal-ion fraction span size mismatch");
  if (!(dt_seconds >= 0.0) || !std::isfinite(dt_seconds)) {
    throw std::invalid_argument("metal-ion advance requires a finite non-negative interval");
  }
  if (dt_seconds == 0.0) return SourceAdvanceStatus::kCompleted;
  std::array<double, k_max_stages> ionize{};
  std::array<double, k_max_stages> recombine{};
  rates(element, env, std::span<double>(ionize.data(), stages), std::span<double>(recombine.data(), stages), report);
  long double elapsed = 0.0L;
  double h = dt_seconds;
  std::uint32_t substeps = 0U;
  std::array<double, k_max_stages> full{};
  std::array<double, k_max_stages> half{};
  while (elapsed < static_cast<long double>(dt_seconds)) {
    if (substeps++ >= m_max_substeps) return SourceAdvanceStatus::kRetryShorterInterval;
    const double remaining = static_cast<double>(static_cast<long double>(dt_seconds) - elapsed);
    h = std::min(h, remaining);
    if (!(h > 1.0e-14 * dt_seconds)) return SourceAdvanceStatus::kRetryShorterInterval;
    std::copy(fractions.begin(), fractions.end(), full.begin());
    std::copy(fractions.begin(), fractions.end(), half.begin());
    const std::span<double> full_span(full.data(), stages);
    const std::span<double> half_span(half.data(), stages);
    backwardEulerStep(std::span<const double>(ionize.data(), stages),
                      std::span<const double>(recombine.data(), stages), h, full_span);
    backwardEulerStep(std::span<const double>(ionize.data(), stages),
                      std::span<const double>(recombine.data(), stages), 0.5 * h, half_span);
    backwardEulerStep(std::span<const double>(ionize.data(), stages),
                      std::span<const double>(recombine.data(), stages), 0.5 * h, half_span);
    double error = 0.0;
    for (std::size_t k = 0; k < stages; ++k) {
      const double scale = m_relative_tolerance * std::max(std::abs(half[k]), std::abs(fractions[k])) +
          m_absolute_tolerance;
      error = std::max(error, std::abs(half[k] - full[k]) / scale);
    }
    if (!std::isfinite(error)) return SourceAdvanceStatus::kRetryShorterInterval;
    if (error > 1.0) {
      ++report.rejected_substeps;
      h *= std::max(0.2, 0.9 / std::sqrt(error));
      continue;
    }
    // Richardson extrapolation of the first-order pair; negative components
    // fall back to the (positive) two-half-step solution.
    double sum = 0.0;
    for (std::size_t k = 0; k < stages; ++k) {
      const double extrapolated = 2.0 * half[k] - full[k];
      fractions[k] = extrapolated >= 0.0 ? extrapolated : half[k];
      sum += fractions[k];
    }
    if (!(sum > 0.0)) return SourceAdvanceStatus::kRetryShorterInterval;
    report.normalization_residual = std::max(report.normalization_residual, std::abs(sum - 1.0));
    for (std::size_t k = 0; k < stages; ++k) fractions[k] /= sum;
    report.max_error_estimate = std::max(report.max_error_estimate, error);
    ++report.accepted_substeps;
    elapsed += static_cast<long double>(h);
    h *= error > 0.0 ? std::min(4.0, 0.9 / std::sqrt(error)) : 4.0;
  }
  return SourceAdvanceStatus::kCompleted;
}

void MetalIonLadderSolver::equilibrium(MetalIonElement element, std::span<double> fractions,
                                       const MetalIonEnvironment& env, MetalIonAdvanceReport& report) const {
  const std::size_t stages = metalIonStageCount(element);
  if (fractions.size() != stages) throw std::invalid_argument("metal-ion fraction span size mismatch");
  std::array<double, k_max_stages> ionize{};
  std::array<double, k_max_stages> recombine{};
  rates(element, env, std::span<double>(ionize.data(), stages), std::span<double>(recombine.data(), stages), report);
  // Detailed balance between adjacent stages: x_{k+1} R_{k+1} = x_k I_k.
  // Log-space recurrence avoids overflow over 27 Fe stages.
  std::array<double, k_max_stages> log_x{};
  log_x[0] = 0.0;
  for (std::size_t k = 0; k + 1U < stages; ++k) {
    if (!(ionize[k] > 0.0)) {
      log_x[k + 1U] = -std::numeric_limits<double>::infinity();
    } else if (!(recombine[k + 1U] > 0.0)) {
      // No recombination out of k+1: everything below is depleted.
      for (std::size_t j = 0; j <= k; ++j) log_x[j] = -std::numeric_limits<double>::infinity();
      log_x[k + 1U] = 0.0;
    } else {
      log_x[k + 1U] = log_x[k] + std::log(ionize[k]) - std::log(recombine[k + 1U]);
    }
  }
  double max_log = -std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k < stages; ++k) max_log = std::max(max_log, log_x[k]);
  double sum = 0.0;
  for (std::size_t k = 0; k < stages; ++k) {
    fractions[k] = std::isfinite(log_x[k]) ? std::exp(log_x[k] - max_log) : 0.0;
    sum += fractions[k];
  }
  for (std::size_t k = 0; k < stages; ++k) fractions[k] /= sum;
}

double MetalIonLadderSolver::cooling(MetalIonElement element, std::span<const double> stage_cm3,
                                     const MetalIonEnvironment& env, MetalIonAdvanceReport& report) const {
  const std::size_t stages = metalIonStageCount(element);
  if (stage_cm3.size() != stages) throw std::invalid_argument("metal-ion stage span size mismatch");
  bool clamped = false;
  double loss = 0.0;
  for (std::size_t k = 0; k < stages; ++k) {
    if (!(stage_cm3[k] > 0.0)) continue;
    const MetalIonStageRates r = m_table->evaluate(element, k, env.temperature_k, env.redshift, m_policy, clamped);
    loss += stage_cm3[k] * r.cooling_erg_cm3_s;
  }
  if (clamped) ++report.range_clamps;
  return std::max(env.electron_cm3, 0.0) * loss;
}

double MetalIonLadderSolver::photoheating(MetalIonElement element, std::span<const double> stage_cm3,
                                          const MetalIonEnvironment& env, MetalIonAdvanceReport& report) const {
  if (!m_table->hasPhotoRates(element)) return 0.0;
  const std::size_t stages = metalIonStageCount(element);
  if (stage_cm3.size() != stages) throw std::invalid_argument("metal-ion stage span size mismatch");
  bool clamped = false;
  double heat = 0.0;
  for (std::size_t k = 0; k + 1U < stages; ++k) {
    if (!(stage_cm3[k] > 0.0)) continue;
    const MetalIonStageRates r = m_table->evaluate(element, k, env.temperature_k, env.redshift, m_policy, clamped);
    heat += stage_cm3[k] * r.photoheating_erg_s;
  }
  if (clamped) ++report.range_clamps;
  return std::clamp(env.photo_attenuation, 0.0, 1.0) * heat;
}

double metalIonElectronsPerNucleus(std::span<const double> fractions) noexcept {
  double electrons = 0.0;
  for (std::size_t k = 1; k < fractions.size(); ++k) electrons += static_cast<double>(k) * fractions[k];
  return electrons;
}

}  // namespace cosmosim::physics
