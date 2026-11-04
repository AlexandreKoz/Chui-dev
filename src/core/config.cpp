#include "cosmosim/core/config.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/units.hpp"

namespace cosmosim::core {
namespace {

struct ParsedEntry {
  std::string value;
  int line_number = 0;
};

struct SectionRequirement {
  std::string name;
  std::vector<std::string> required_fields;
};

[[nodiscard]] double resolveFlatLambdaCdmScaleFactorAfterTime(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double elapsed_time_si) {
  if (!std::isfinite(scale_factor_begin) || scale_factor_begin <= 0.0 ||
      !std::isfinite(elapsed_time_si) || elapsed_time_si < 0.0) {
    throw std::invalid_argument("invalid scale factor or elapsed time for FLRW endpoint resolution");
  }
  if (elapsed_time_si == 0.0) {
    return scale_factor_begin;
  }

  const auto& cosmology = background.config();
  const double omega_m = cosmology.omega_matter;
  const double omega_l = cosmology.omega_lambda;
  const double hubble0 = background.hubble0Si();
  if (omega_m > 0.0 && omega_l > 0.0) {
    const double sqrt_l = std::sqrt(omega_l);
    const double q = std::sqrt(omega_l / omega_m);
    const double y_begin = std::asinh(q * std::pow(scale_factor_begin, 1.5));
    const double y_end = y_begin + 1.5 * hubble0 * sqrt_l * elapsed_time_si;
    if (!std::isfinite(y_end) || y_end > 700.0) {
      throw std::overflow_error("FLRW code-time endpoint exceeds representable scale factor");
    }
    return std::pow(std::sinh(y_end) / q, 2.0 / 3.0);
  }
  if (omega_m > 0.0) {
    const double a32 = std::pow(scale_factor_begin, 1.5) +
        1.5 * hubble0 * std::sqrt(omega_m) * elapsed_time_si;
    return std::pow(a32, 2.0 / 3.0);
  }
  if (omega_l > 0.0) {
    const double exponent = hubble0 * std::sqrt(omega_l) * elapsed_time_si;
    if (!std::isfinite(exponent) || exponent > 700.0) {
      throw std::overflow_error("FLRW code-time endpoint exceeds representable scale factor");
    }
    return scale_factor_begin * std::exp(exponent);
  }
  throw std::invalid_argument("flat LambdaCDM endpoint resolution requires matter and/or Lambda");
}

[[nodiscard]] const std::map<std::string, std::string>& deprecatedAliasRegistry();

[[nodiscard]] std::string trim(const std::string& input) {
  const auto begin = std::find_if_not(input.begin(), input.end(), [](unsigned char c) {
    return std::isspace(c) != 0;
  });
  const auto end = std::find_if_not(input.rbegin(), input.rend(), [](unsigned char c) {
    return std::isspace(c) != 0;
  }).base();

  if (begin >= end) {
    return {};
  }

  return std::string(begin, end);
}

[[nodiscard]] std::string toLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

[[nodiscard]] std::string removeInlineComment(const std::string& line) {
  bool in_single_quotes = false;
  bool in_double_quotes = false;

  const auto is_whitespace_boundary = [&line](std::size_t pos) {
    return pos == 0 || std::isspace(static_cast<unsigned char>(line[pos - 1])) != 0;
  };

  for (std::size_t pos = 0; pos < line.size(); ++pos) {
    const char current = line[pos];
    if (current == '"' && !in_single_quotes) {
      in_double_quotes = !in_double_quotes;
      continue;
    }
    if (current == '\'' && !in_double_quotes) {
      in_single_quotes = !in_single_quotes;
      continue;
    }
    if (in_single_quotes || in_double_quotes) {
      continue;
    }

    if ((current == '#' || current == ';') && is_whitespace_boundary(pos)) {
      return line.substr(0, pos);
    }
    if (current == '/' && pos + 1 < line.size() && line[pos + 1] == '/' &&
        is_whitespace_boundary(pos)) {
      return line.substr(0, pos);
    }
  }

  return line;
}

[[nodiscard]] std::map<std::string, ParsedEntry> parseEntries(const std::string& text) {
  std::istringstream stream(text);
  std::string raw_line;
  int line_number = 0;
  std::string current_section;
  std::map<std::string, ParsedEntry> entries;

  while (std::getline(stream, raw_line)) {
    ++line_number;
    const std::string line = trim(removeInlineComment(raw_line));
    if (line.empty()) {
      continue;
    }

    if (line.front() == '[' && line.back() == ']') {
      current_section = toLower(trim(line.substr(1, line.size() - 2)));
      if (current_section.empty()) {
        throw ConfigError("line " + std::to_string(line_number) + ": empty section name");
      }
      continue;
    }

    std::string key;
    std::string value;
    const std::size_t equal_pos = line.find('=');
    if (equal_pos != std::string::npos) {
      key = trim(line.substr(0, equal_pos));
      value = trim(line.substr(equal_pos + 1));
    } else {
      std::istringstream token_stream(line);
      token_stream >> key;
      std::getline(token_stream, value);
      value = trim(value);
    }

    key = toLower(trim(key));
    const bool allow_empty_value = (equal_pos != std::string::npos);
    if (key.empty() || (value.empty() && !allow_empty_value)) {
      throw ConfigError("line " + std::to_string(line_number) + ": expected key and value");
    }

    if (!current_section.empty() && key.find('.') == std::string::npos) {
      key = current_section + "." + key;
    }

    auto [it, inserted] = entries.emplace(key, ParsedEntry{value, line_number});
    if (!inserted) {
      throw ConfigError("line " + std::to_string(line_number) + ": duplicate key '" + key +
                        "' (first at line " + std::to_string(it->second.line_number) + ")");
    }
  }

  return entries;
}

[[nodiscard]] std::set<std::string> inputSections(const std::map<std::string, ParsedEntry>& entries) {
  std::set<std::string> sections;
  for (const auto& [key, _] : entries) {
    const std::size_t dot = key.find('.');
    if (dot == std::string::npos) {
      continue;
    }
    sections.insert(key.substr(0, dot));
  }
  return sections;
}

[[nodiscard]] const std::vector<SectionRequirement>& executableSchemaRequirements() {
  static const std::vector<SectionRequirement> requirements = {
      {"mode", {"mode.mode"}},
  };
  return requirements;
}

void validateRequiredSchema(const std::map<std::string, ParsedEntry>& entries) {
  const std::set<std::string> sections = inputSections(entries);
  for (const SectionRequirement& requirement : executableSchemaRequirements()) {
    const bool has_legacy_mode_root = (requirement.name == "mode" && entries.contains("mode"));
    if (!sections.contains(requirement.name) && !has_legacy_mode_root) {
      throw ConfigError("missing required section '" + requirement.name + "'");
    }
    for (const std::string& field : requirement.required_fields) {
      const bool has_legacy_mode_field = (field == "mode.mode" && entries.contains("mode"));
      if (!entries.contains(field) && !has_legacy_mode_field) {
        throw ConfigError("missing required field '" + field + "'");
      }
    }
  }
}

void applyDeprecatedAliases(
    std::map<std::string, ParsedEntry>& entries,
    std::set<std::string>& consumed,
    FrozenConfig& frozen) {
  for (const auto& [legacy_key, canonical_key] : deprecatedAliasRegistry()) {
    const auto legacy_it = entries.find(legacy_key);
    if (legacy_it == entries.end()) {
      continue;
    }
    const auto canonical_it = entries.find(canonical_key);
    if (canonical_it != entries.end()) {
      if (trim(canonical_it->second.value) != trim(legacy_it->second.value)) {
        throw ConfigError(
            "alias conflict at path '" + canonical_key + "': keys '" + legacy_key + "' and '" +
            canonical_key + "' provide different values");
      }
      throw ConfigError("deprecated key '" + legacy_key + "' cannot be combined with '" +
                        canonical_key + "'");
    }
    entries.emplace(canonical_key, legacy_it->second);
    consumed.insert(legacy_key);
    const std::string note =
        "deprecated key '" + legacy_key + "' mapped to '" + canonical_key + "'";
    frozen.provenance.deprecation_warnings.push_back(note);
    frozen.user_config.alias_resolution_notes.push_back(note);
  }
}

[[nodiscard]] bool parseBool(const std::string& value, const std::string& key) {
  const std::string lower = toLower(trim(value));
  if (lower == "true") {
    return true;
  }
  if (lower == "false") {
    return false;
  }
  throw ConfigError("key '" + key + "': invalid boolean value '" + value + "'");
}

template <typename T>
[[nodiscard]] T parseNumber(const std::string& value, const std::string& key) {
  const std::string trimmed = trim(value);
  T number{};
  const char* begin = trimmed.data();
  const char* end = trimmed.data() + trimmed.size();
  const auto [ptr, ec] = std::from_chars(begin, end, number);
  if (ec != std::errc{} || ptr != end) {
    throw ConfigError("key '" + key + "': invalid numeric value '" + value + "'");
  }
  return number;
}

[[nodiscard]] double parseFloating(const std::string& value, const std::string& key) {
  const std::string trimmed = trim(value);
  if (trimmed.empty()) {
    throw ConfigError("key '" + key + "': invalid finite floating-point value '" + value + "'");
  }
  char* consumed = nullptr;
  const double number = std::strtod(trimmed.c_str(), &consumed);
  if (consumed != trimmed.c_str() + static_cast<std::ptrdiff_t>(trimmed.size()) ||
      !std::isfinite(number)) {
    throw ConfigError("key '" + key + "': invalid finite floating-point value '" + value + "'");
  }
  return number;
}

[[nodiscard]] double parseOptionalFinite(
    const std::string& value,
    const std::string& key) {
  if (toLower(trim(value)) == "unspecified") {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return parseFloating(value, key);
}

void requireFinite(double value, const std::string& key) {
  if (!std::isfinite(value)) {
    throw ConfigError("key '" + key + "': value must be finite");
  }
}

void requireAllFinite(std::initializer_list<std::pair<double, const char*>> values) {
  for (const auto& [value, key] : values) {
    requireFinite(value, key);
  }
}

[[nodiscard]] std::pair<double, std::string> splitMagnitudeUnit(const std::string& value) {
  std::istringstream stream(value);
  double magnitude = 0.0;
  stream >> magnitude;
  if (stream.fail() || !std::isfinite(magnitude)) {
    throw ConfigError("invalid value '" + value + "': expected finite numeric magnitude");
  }

  std::string unit;
  stream >> unit;
  stream >> std::ws;
  if (!stream.eof()) {
    std::string trailing;
    std::getline(stream, trailing);
    throw ConfigError(
        "invalid value '" + value + "': unexpected trailing token(s) after magnitude/unit");
  }
  unit = toLower(unit);
  return {magnitude, unit};
}

[[nodiscard]] double lengthToMpc(double value, const std::string& unit) {
  const std::string normalized = toLower(unit);
  if (normalized.empty() || normalized == "mpc") {
    return value;
  }
  if (normalized == "kpc") {
    return value * 1.0e-3;
  }
  throw ConfigError("unsupported length unit '" + unit + "' (supported: mpc, kpc)");
}

[[nodiscard]] double parseLengthMpc(
    const std::string& value,
    const std::string& default_unit,
    const std::string& key) {
  const auto [magnitude, provided_unit] = splitMagnitudeUnit(value);
  const std::string unit = provided_unit.empty() ? default_unit : provided_unit;
  try {
    return lengthToMpc(magnitude, unit);
  } catch (const ConfigError&) {
    throw ConfigError("key '" + key + "': unsupported unit in value '" + value + "'");
  }
}

[[nodiscard]] double parseLengthKpc(
    const std::string& value,
    const std::string& default_unit,
    const std::string& key) {
  const double in_mpc = parseLengthMpc(value, default_unit, key);
  return in_mpc * 1.0e3;
}

[[nodiscard]] double massUnitSi(const std::string& unit) {
  const std::string normalized = toLower(unit);
  if (normalized == "kg") {
    return 1.0;
  }
  if (normalized == "g") {
    return constants::k_gram_si;
  }
  if (normalized == "msun") {
    return constants::k_solar_mass_si;
  }
  throw ConfigError(
      "unsupported mass unit '" + unit + "' (supported: kg, g, msun)");
}

[[nodiscard]] double parseMassCode(
    const std::string& value,
    const std::string& default_code_mass_unit,
    const std::string& key) {
  const auto [magnitude, provided_unit] = splitMagnitudeUnit(value);
  const std::string source_unit =
      provided_unit.empty() ? default_code_mass_unit : provided_unit;
  try {
    return magnitude * massUnitSi(source_unit) /
        massUnitSi(default_code_mass_unit);
  } catch (const ConfigError&) {
    throw ConfigError(
        "key '" + key + "': unsupported unit in value '" + value + "'");
  }
}

[[nodiscard]] double timeUnitSi(const std::string& unit) {
  const std::string normalized = toLower(unit);
  if (normalized == "s" || normalized == "sec" || normalized == "second") return 1.0;
  if (normalized == "yr" || normalized == "year") return (constants::k_seconds_per_megayear / 1.0e6);
  if (normalized == "myr") return 1.0e6 * (constants::k_seconds_per_megayear / 1.0e6);
  if (normalized == "gyr") return 1.0e9 * (constants::k_seconds_per_megayear / 1.0e6);
  throw ConfigError("unsupported time unit '" + unit + "' (supported: code, s, yr, Myr, Gyr)");
}

[[nodiscard]] double parseTimeCode(
    const std::string& value,
    const UnitsConfig& unit_config,
    const std::string& key) {
  const auto [magnitude, provided_unit] = splitMagnitudeUnit(value);
  if (provided_unit.empty() || provided_unit == "code") return magnitude;
  try {
    const UnitSystem units = makeUnitSystem(
        unit_config.length_unit, unit_config.mass_unit, unit_config.velocity_unit);
    return magnitude * timeUnitSi(provided_unit) / units.timeSiPerCode();
  } catch (const ConfigError&) {
    throw ConfigError("key '" + key + "': unsupported unit in value '" + value + "'");
  }
}


[[nodiscard]] double parseNumberDensityCgs(
    const std::string& value,
    const std::string& key) {
  const auto [magnitude, provided_unit] = splitMagnitudeUnit(value);
  const std::string unit = toLower(provided_unit);
  if (unit.empty() || unit == "cm^-3" || unit == "cm-3" || unit == "1/cm3" ||
      unit == "cm^{-3}") {
    return magnitude;
  }
  if (unit == "m^-3" || unit == "m-3" || unit == "1/m3" || unit == "m^{-3}") {
    return magnitude * 1.0e-6;
  }
  throw ConfigError(
      "key '" + key + "': unsupported number-density unit in value '" + value +
      "' (supported: cm^-3, m^-3)");
}

[[nodiscard]] double parseTemperatureKelvin(
    const std::string& value,
    const std::string& key) {
  const auto [magnitude, provided_unit] = splitMagnitudeUnit(value);
  const std::string unit = toLower(provided_unit);
  if (unit.empty() || unit == "k" || unit == "kelvin") {
    return magnitude;
  }
  throw ConfigError(
      "key '" + key + "': unsupported temperature unit in value '" + value +
      "' (supported: K)");
}

[[nodiscard]] double parseSpecificEnergyCode(
    const std::string& value,
    const UnitsConfig& unit_config,
    const std::string& key) {
  const auto [magnitude, provided_unit] = splitMagnitudeUnit(value);
  const std::string unit = toLower(provided_unit);
  if (unit.empty() || unit == "code") {
    return magnitude;
  }
  const UnitSystem units = makeUnitSystem(
      unit_config.length_unit, unit_config.mass_unit, unit_config.velocity_unit);
  double specific_energy_si = 0.0;
  if (unit == "j/kg" || unit == "m2/s2" || unit == "m^2/s^2") {
    specific_energy_si = magnitude;
  } else if (unit == "erg/g" || unit == "cm2/s2" || unit == "cm^2/s^2") {
    specific_energy_si = magnitude * 1.0e-4;
  } else {
    throw ConfigError(
        "key '" + key + "': unsupported specific-energy unit in value '" + value +
        "' (supported: code, J/kg, erg/g)");
  }
  const double velocity_si_per_code = units.velocity_si_per_code;
  return specific_energy_si / (velocity_si_per_code * velocity_si_per_code);
}

[[nodiscard]] std::string requireString(
    std::map<std::string, ParsedEntry>& entries,
    std::set<std::string>& consumed,
    const std::string& key,
    const std::string& default_value) {
  const auto it = entries.find(key);
  if (it == entries.end()) {
    return default_value;
  }
  consumed.insert(key);
  return trim(it->second.value);
}

[[nodiscard]] std::string sanitizeStem(const std::string& value, const std::string& key) {
  if (value.empty()) {
    throw ConfigError("key '" + key + "': value cannot be empty");
  }
  for (const char c : value) {
    const bool valid = std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-';
    if (!valid) {
      throw ConfigError(
          "key '" + key + "': only [a-zA-Z0-9_-] are allowed for stable output naming");
    }
  }
  return value;
}

[[nodiscard]] std::string modeToLowerString(SimulationMode mode) {
  switch (mode) {
    case SimulationMode::kCosmoCube:
      return "cosmo_cube";
    case SimulationMode::kZoomIn:
      return "zoom_in";
    case SimulationMode::kIsolatedGalaxy:
      return "isolated_galaxy";
    case SimulationMode::kIsolatedCluster:
      return "isolated_cluster";
  }
  throw ConfigError("unhandled SimulationMode enum value during serialization");
}

[[nodiscard]] SimulationMode parseMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "cosmo_cube") {
    return SimulationMode::kCosmoCube;
  }
  if (lower == "zoom_in") {
    return SimulationMode::kZoomIn;
  }
  if (lower == "isolated_galaxy") {
    return SimulationMode::kIsolatedGalaxy;
  }
  if (lower == "isolated_cluster") {
    return SimulationMode::kIsolatedCluster;
  }
  throw ConfigError("key 'mode.mode': invalid mode '" + value + "'");
}

[[nodiscard]] InitialConditionConvention parseInitialConditionConvention(
    const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "generated") {
    return InitialConditionConvention::kGenerated;
  }
  if (lower == "chui_canonical_v1") {
    return InitialConditionConvention::kChuiCanonicalV1;
  }
  if (lower == "gadget_arepo_bridge_v1") {
    return InitialConditionConvention::kGadgetArepoBridgeV1;
  }
  if (lower == "manifest_v1") {
    return InitialConditionConvention::kManifestV1;
  }
  throw ConfigError(
      "key 'mode.ic_convention': invalid value '" + value +
      "' (supported: generated, chui_canonical_v1, gadget_arepo_bridge_v1, manifest_v1)");
}

[[nodiscard]] InitialConditionCoordinateFrame parseInitialConditionCoordinateFrame(
    const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "unspecified") {
    return InitialConditionCoordinateFrame::kUnspecified;
  }
  if (lower == "comoving") {
    return InitialConditionCoordinateFrame::kComoving;
  }
  if (lower == "physical") {
    return InitialConditionCoordinateFrame::kPhysical;
  }
  throw ConfigError(
      "key 'mode.ic_bridge_coordinate_frame': invalid value '" + value +
      "' (supported: comoving, physical)");
}

[[nodiscard]] InitialConditionVelocityConvention
parseInitialConditionVelocityConvention(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "unspecified") {
    return InitialConditionVelocityConvention::kUnspecified;
  }
  if (lower == "physical_peculiar") {
    return InitialConditionVelocityConvention::kPhysicalPeculiar;
  }
  if (lower == "sqrt_a_scaled_peculiar") {
    return InitialConditionVelocityConvention::kSqrtAScaledPeculiar;
  }
  if (lower == "gadget_arepo_stored_peculiar") {
    return InitialConditionVelocityConvention::kGadgetArepoStoredPeculiar;
  }
  if (lower == "comoving_coordinate_rate") {
    return InitialConditionVelocityConvention::kComovingCoordinateRate;
  }
  throw ConfigError(
      "key 'mode.ic_bridge_velocity_convention': invalid value '" + value +
      "' (supported: physical_peculiar, sqrt_a_scaled_peculiar, "
      "gadget_arepo_stored_peculiar, comoving_coordinate_rate)");
}

[[nodiscard]] InitialConditionSpeciesPolicy parseInitialConditionSpeciesPolicy(
    const std::string& value,
    const std::string& key) {
  const std::string lower = toLower(trim(value));
  if (lower == "reject") {
    return InitialConditionSpeciesPolicy::kReject;
  }
  if (lower == "dark_matter") {
    return InitialConditionSpeciesPolicy::kDarkMatter;
  }
  if (lower == "star") {
    return InitialConditionSpeciesPolicy::kStar;
  }
  if (lower == "black_hole") {
    return InitialConditionSpeciesPolicy::kBlackHole;
  }
  if (lower == "tracer") {
    return InitialConditionSpeciesPolicy::kTracer;
  }
  throw ConfigError(
      "key '" + key + "': invalid value '" + value +
      "' (supported: reject, dark_matter, star, black_hole, tracer)");
}

[[nodiscard]] InitialConditionMissingFieldPolicy
parseInitialConditionMissingFieldPolicy(
    const std::string& value,
    const std::string& key) {
  const std::string lower = toLower(trim(value));
  if (lower == "reject") {
    return InitialConditionMissingFieldPolicy::kReject;
  }
  if (lower == "reconstruct") {
    return InitialConditionMissingFieldPolicy::kReconstruct;
  }
  if (lower == "use_config_value") {
    return InitialConditionMissingFieldPolicy::kUseConfigValue;
  }
  if (lower == "dialect_defined_default") {
    return InitialConditionMissingFieldPolicy::kDialectDefinedDefault;
  }
  if (lower == "preserve_unavailable") {
    return InitialConditionMissingFieldPolicy::kPreserveUnavailable;
  }
  throw ConfigError(
      "key '" + key + "': invalid value '" + value +
      "' (supported: reject, reconstruct, use_config_value, "
      "dialect_defined_default, preserve_unavailable)");
}

[[nodiscard]] GravitySolver parseGravitySolver(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "treepm") {
    return GravitySolver::kTreePm;
  }
  throw ConfigError("key 'numerics.gravity_solver': invalid value '" + value + "'");
}

[[nodiscard]] HydroSolver parseHydroSolver(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "godunov_fv") {
    return HydroSolver::kGodunovFv;
  }
  throw ConfigError("key 'numerics.hydro_solver': invalid value '" + value + "'");
}

[[nodiscard]] HydroReconstructionMode parseHydroReconstructionMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "piecewise_constant") return HydroReconstructionMode::kPiecewiseConstant;
  if (lower == "piecewise_linear") return HydroReconstructionMode::kPiecewiseLinear;
  if (lower == "piecewise_parabolic") return HydroReconstructionMode::kPiecewiseParabolic;
  throw ConfigError("key 'numerics.hydro_reconstruction': invalid value '" + value +
                    "' (expected piecewise_constant, piecewise_linear, or piecewise_parabolic)");
}

[[nodiscard]] HydroSlopeLimiterMode parseHydroSlopeLimiterMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "minmod") return HydroSlopeLimiterMode::kMinmod;
  if (lower == "mc") return HydroSlopeLimiterMode::kMonotonizedCentral;
  if (lower == "van_leer") return HydroSlopeLimiterMode::kVanLeer;
  throw ConfigError("key 'numerics.hydro_slope_limiter': invalid value '" + value +
                    "' (expected minmod, mc, or van_leer)");
}

[[nodiscard]] HydroIntegrationMode parseHydroIntegrationMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "muscl_hancock") return HydroIntegrationMode::kMusclHancock;
  if (lower == "ctu") return HydroIntegrationMode::kCornerTransportUpwind;
  throw ConfigError("key 'numerics.hydro_integration_scheme': invalid value '" + value +
                    "' (expected muscl_hancock or ctu)");
}

[[nodiscard]] HydroRiemannSolverMode parseHydroRiemannSolverMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "hllc") return HydroRiemannSolverMode::kHllc;
  if (lower == "hlle") return HydroRiemannSolverMode::kHlle;
  throw ConfigError("key 'numerics.hydro_riemann_solver': invalid value '" + value +
                    "' (expected hllc or hlle)");
}

[[nodiscard]] HydroDualEnergyConfigMode parseHydroDualEnergyConfigMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "off") return HydroDualEnergyConfigMode::kOff;
  if (lower == "thermal_fraction") return HydroDualEnergyConfigMode::kThermalFraction;
  throw ConfigError("key 'numerics.hydro_dual_energy': invalid value '" + value +
                    "' (expected off or thermal_fraction)");
}

[[nodiscard]] SspWindModelConfig parseSspWindModelConfig(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") return SspWindModelConfig::kNone;
  if (lower == "table_v1") return SspWindModelConfig::kTableV1;
  throw ConfigError("key 'physics.ssp_wind_model': invalid value '" + value + "' (expected none or table_v1)");
}

[[nodiscard]] AmrTimeSteppingConfigMode parseAmrTimeSteppingConfigMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "synchronized") return AmrTimeSteppingConfigMode::kSynchronized;
  if (lower == "recursive_subcycling") return AmrTimeSteppingConfigMode::kRecursiveSubcycling;
  throw ConfigError("key 'numerics.amr_time_stepping': invalid value '" + value +
                    "' (expected synchronized or recursive_subcycling)");
}

[[nodiscard]] AmrRegridConfigMode parseAmrRegridConfigMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "off") return AmrRegridConfigMode::kOff;
  if (lower == "automatic") return AmrRegridConfigMode::kAutomatic;
  throw ConfigError("key 'numerics.amr_regrid_mode': invalid value '" + value +
                    "' (expected off or automatic)");
}

[[nodiscard]] NodeMemoryAdmissionMode parseNodeMemoryAdmissionMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "off") return NodeMemoryAdmissionMode::kOff;
  if (lower == "shared_equal") return NodeMemoryAdmissionMode::kSharedEqual;
  if (lower == "shared_pool") return NodeMemoryAdmissionMode::kSharedPool;
  throw ConfigError("key 'parallel.node_memory_admission': invalid value '" + value +
                    "' (expected off, shared_equal or shared_pool)");
}

[[nodiscard]] TreePmAssignmentScheme parseTreePmAssignmentScheme(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "cic") {
    return TreePmAssignmentScheme::kCic;
  }
  if (lower == "tsc") {
    return TreePmAssignmentScheme::kTsc;
  }
  throw ConfigError(
      "key 'numerics.treepm_assignment_scheme': invalid value '" + value +
      "' (supported: cic, tsc)");
}

[[nodiscard]] const char* treePmPmReusePolicyName(TreePmPmReusePolicy policy) noexcept {
  switch (policy) {
    case TreePmPmReusePolicy::kStaleWithinCadence: return "stale_within_cadence";
    case TreePmPmReusePolicy::kMultirateOuterKick: return "multirate_outer_kick";
    case TreePmPmReusePolicy::kExactSourceState: break;
  }
  return "exact";
}

[[nodiscard]] TreePmPmReusePolicy parseTreePmPmReusePolicy(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "exact") return TreePmPmReusePolicy::kExactSourceState;
  if (lower == "stale_within_cadence") return TreePmPmReusePolicy::kStaleWithinCadence;
  if (lower == "multirate_outer_kick") return TreePmPmReusePolicy::kMultirateOuterKick;
  throw ConfigError(
      "key 'numerics.treepm_pm_reuse_policy': invalid value '" + value +
      "' (supported: exact, stale_within_cadence, multirate_outer_kick)");
}

[[nodiscard]] const char* multirateRollbackMemoryPolicyName(MultirateRollbackMemoryPolicy policy) noexcept {
  return policy == MultirateRollbackMemoryPolicy::kFail ? "fail" : "single_block_fallback";
}

[[nodiscard]] MultirateRollbackMemoryPolicy parseMultirateRollbackMemoryPolicy(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "single_block_fallback") return MultirateRollbackMemoryPolicy::kSingleBlockFallback;
  if (lower == "fail") return MultirateRollbackMemoryPolicy::kFail;
  throw ConfigError(
      "key 'numerics.treepm_multirate_rollback_memory_policy': invalid value '" + value +
      "' (supported: single_block_fallback, fail)");
}

[[nodiscard]] TreePmForceGradient parseTreePmForceGradient(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "spectral") return TreePmForceGradient::kSpectral;
  if (lower == "finite_difference_4") return TreePmForceGradient::kFiniteDifference4;
  throw ConfigError(
      "key 'numerics.treepm_pm_force_gradient': invalid value '" + value +
      "' (supported: spectral, finite_difference_4)");
}

[[nodiscard]] TreePmOpeningCriterion parseTreePmOpeningCriterion(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "geometric") {
    return TreePmOpeningCriterion::kGeometric;
  }
  if (lower == "com_distance") {
    return TreePmOpeningCriterion::kComDistance;
  }
  if (lower == "relative_force_error") {
    return TreePmOpeningCriterion::kRelativeForceError;
  }
  throw ConfigError(
      "key 'numerics.treepm_tree_opening_criterion': invalid value '" + value +
      "' (supported: geometric, com_distance, relative_force_error)");
}

[[nodiscard]] std::string treePmOpeningCriterionToString(TreePmOpeningCriterion criterion) {
  switch (criterion) {
    case TreePmOpeningCriterion::kGeometric:
      return "geometric";
    case TreePmOpeningCriterion::kComDistance:
      return "com_distance";
    case TreePmOpeningCriterion::kRelativeForceError:
      return "relative_force_error";
  }
  throw ConfigError("unhandled TreePmOpeningCriterion enum value during serialization");
}

[[nodiscard]] HierarchicalCoarseDisplacementCriterion parseHierarchicalCoarseDisplacementCriterion(
    const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "max_speed_mesh") return HierarchicalCoarseDisplacementCriterion::kMaxSpeedMesh;
  if (lower == "rms_split_scale") return HierarchicalCoarseDisplacementCriterion::kRmsSplitScale;
  throw ConfigError(
      "key 'numerics.hierarchical_coarse_displacement_criterion': invalid value '" + value +
      "' (supported: max_speed_mesh, rms_split_scale)");
}

[[nodiscard]] std::string hierarchicalCoarseDisplacementCriterionToString(
    HierarchicalCoarseDisplacementCriterion criterion) {
  switch (criterion) {
    case HierarchicalCoarseDisplacementCriterion::kMaxSpeedMesh:
      return "max_speed_mesh";
    case HierarchicalCoarseDisplacementCriterion::kRmsSplitScale:
      return "rms_split_scale";
  }
  throw ConfigError("unhandled HierarchicalCoarseDisplacementCriterion enum value during serialization");
}

[[nodiscard]] TreePmAdaptiveErrorEstimator parseTreePmAdaptiveErrorEstimator(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "residual_derivative") return TreePmAdaptiveErrorEstimator::kResidualDerivative;
  if (lower == "relative_monopole") return TreePmAdaptiveErrorEstimator::kRelativeMonopole;
  throw ConfigError(
      "key 'numerics.treepm_adaptive_error_estimator': invalid value '" + value +
      "' (supported: residual_derivative, relative_monopole)");
}

[[nodiscard]] std::string treePmAdaptiveErrorEstimatorToString(TreePmAdaptiveErrorEstimator estimator) {
  switch (estimator) {
    case TreePmAdaptiveErrorEstimator::kResidualDerivative:
      return "residual_derivative";
    case TreePmAdaptiveErrorEstimator::kRelativeMonopole:
      return "relative_monopole";
  }
  throw ConfigError("unhandled TreePmAdaptiveErrorEstimator enum value during serialization");
}

[[nodiscard]] TreePmCutoffNodePolicy parseTreePmCutoffNodePolicy(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "contain") return TreePmCutoffNodePolicy::kContain;
  if (lower == "com_distance") return TreePmCutoffNodePolicy::kComDistance;
  throw ConfigError(
      "key 'numerics.treepm_cutoff_node_policy': invalid value '" + value +
      "' (supported: contain, com_distance)");
}

[[nodiscard]] std::string treePmCutoffNodePolicyToString(TreePmCutoffNodePolicy policy) {
  switch (policy) {
    case TreePmCutoffNodePolicy::kContain:
      return "contain";
    case TreePmCutoffNodePolicy::kComDistance:
      return "com_distance";
  }
  throw ConfigError("unhandled TreePmCutoffNodePolicy enum value during serialization");
}

[[nodiscard]] PmDecompositionMode parsePmDecompositionMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "slab") {
    return PmDecompositionMode::kSlab;
  }
  if (lower == "pencil") {
    return PmDecompositionMode::kPencil;
  }
  throw ConfigError(
      "key 'numerics.treepm_pm_decomposition_mode': invalid value '" + value +
      "' (supported: slab, pencil)");
}

[[nodiscard]] CoordinateFrame parseCoordinateFrame(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "comoving") {
    return CoordinateFrame::kComoving;
  }
  if (lower == "physical") {
    return CoordinateFrame::kPhysical;
  }
  throw ConfigError("key 'units.coordinate_frame': invalid value '" + value + "'");
}

[[nodiscard]] ModeHydroBoundary parseModeHydroBoundary(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "auto") {
    return ModeHydroBoundary::kAuto;
  }
  if (lower == "periodic") {
    return ModeHydroBoundary::kPeriodic;
  }
  if (lower == "open") {
    return ModeHydroBoundary::kOpen;
  }
  if (lower == "reflective") {
    return ModeHydroBoundary::kReflective;
  }
  throw ConfigError("key 'mode.hydro_boundary': invalid value '" + value + "'");
}

[[nodiscard]] ModeGravityBoundary parseModeGravityBoundary(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "auto") {
    return ModeGravityBoundary::kAuto;
  }
  if (lower == "periodic") {
    return ModeGravityBoundary::kPeriodic;
  }
  if (lower == "isolated_monopole") {
    return ModeGravityBoundary::kIsolatedMonopole;
  }
  throw ConfigError("key 'mode.gravity_boundary': invalid value '" + value + "'");
}

[[nodiscard]] ZoomLongRangeStrategy parseZoomLongRangeStrategy(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "disabled") {
    return ZoomLongRangeStrategy::kDisabled;
  }
  if (lower == "global_coarse_plus_focused_highres_correction") {
    return ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection;
  }
  throw ConfigError("key 'mode.zoom_long_range_strategy': invalid value '" + value + "'");
}

[[nodiscard]] std::string zoomLongRangeStrategyToString(ZoomLongRangeStrategy strategy) {
  switch (strategy) {
    case ZoomLongRangeStrategy::kDisabled:
      return "disabled";
    case ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection:
      return "global_coarse_plus_focused_highres_correction";
  }
  throw ConfigError("unhandled ZoomLongRangeStrategy enum value during serialization");
}

[[nodiscard]] ScienceSnapshotLayout parseScienceSnapshotLayout(const std::string& value) {
  if (value == "auto") return ScienceSnapshotLayout::kAuto;
  if (value == "single") return ScienceSnapshotLayout::kSingle;
  if (value == "aggregated") return ScienceSnapshotLayout::kAggregated;
  if (value == "sharded") return ScienceSnapshotLayout::kSharded;
  throw ConfigError("output.snapshot_layout must be one of: auto, single, aggregated, sharded");
}

[[nodiscard]] std::string_view scienceSnapshotLayoutLabel(ScienceSnapshotLayout layout) {
  switch (layout) {
    case ScienceSnapshotLayout::kAuto: return "auto";
    case ScienceSnapshotLayout::kSingle: return "single";
    case ScienceSnapshotLayout::kAggregated: return "aggregated";
    case ScienceSnapshotLayout::kSharded: return "sharded";
  }
  throw ConfigError("unhandled output.snapshot_layout enum value");
}

[[nodiscard]] SnapshotVerificationPolicy parseSnapshotVerificationPolicy(const std::string& value) {
  if (value == "full_readback") return SnapshotVerificationPolicy::kFullReadback;
  if (value == "structural") return SnapshotVerificationPolicy::kStructural;
  throw ConfigError("output.snapshot_verification must be one of: full_readback, structural");
}

[[nodiscard]] std::string_view snapshotVerificationPolicyLabel(SnapshotVerificationPolicy policy) {
  switch (policy) {
    case SnapshotVerificationPolicy::kFullReadback: return "full_readback";
    case SnapshotVerificationPolicy::kStructural: return "structural";
  }
  throw ConfigError("unhandled output.snapshot_verification enum value");
}

[[nodiscard]] StarFormationModelKind parseStarFormationModelKind(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "legacy_schmidt_threshold") {
    return StarFormationModelKind::kLegacySchmidtThreshold;
  }
  if (lower == "adaptive_bound_jeans") {
    return StarFormationModelKind::kAdaptiveBoundJeans;
  }
  if (lower == "effective_multiphase_tng_like") {
    return StarFormationModelKind::kEffectiveMultiphaseTngLike;
  }
  throw ConfigError(
      "physics.star_formation_model must be one of: legacy_schmidt_threshold, adaptive_bound_jeans, effective_multiphase_tng_like");
}

[[nodiscard]] EffectiveIsmEosRelaxation parseEffectiveIsmEosRelaxation(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "instantaneous") return EffectiveIsmEosRelaxation::kInstantaneous;
  if (lower == "finite_timescale") return EffectiveIsmEosRelaxation::kFiniteTimescale;
  throw ConfigError("physics.sf_effective_eos_relaxation must be one of: instantaneous, finite_timescale");
}

[[nodiscard]] EffectiveIsmBirthMassConvention parseEffectiveIsmBirthMassConvention(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "initial_ssp_mass") return EffectiveIsmBirthMassConvention::kInitialSspMass;
  if (lower == "long_lived_mass") return EffectiveIsmBirthMassConvention::kLongLivedMass;
  throw ConfigError("physics.sf_effective_birth_mass_convention must be one of: initial_ssp_mass, long_lived_mass");
}

[[nodiscard]] EffectiveIsmFeedbackCoupling parseEffectiveIsmFeedbackCoupling(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "external_feedback_calibrated") return EffectiveIsmFeedbackCoupling::kExternalFeedbackCalibrated;
  if (lower == "effective_eos_only") return EffectiveIsmFeedbackCoupling::kEffectiveEosOnly;
  throw ConfigError("physics.sf_effective_feedback_coupling must be one of: external_feedback_calibrated, effective_eos_only");
}

[[nodiscard]] StarFormationCollapseTimescale parseStarFormationCollapseTimescale(
    const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "free_fall") {
    return StarFormationCollapseTimescale::kFreeFall;
  }
  if (lower == "minimum_free_fall_or_compression") {
    return StarFormationCollapseTimescale::kMinimumFreeFallOrCompression;
  }
  throw ConfigError(
      "physics.sf_collapse_timescale must be one of: free_fall, minimum_free_fall_or_compression");
}

[[nodiscard]] StarParticleMassPolicy parseStarParticleMassPolicy(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "fixed") {
    return StarParticleMassPolicy::kFixed;
  }
  if (lower == "gas_resolution_fraction") {
    return StarParticleMassPolicy::kGasResolutionFraction;
  }
  throw ConfigError(
      "physics.sf_star_particle_mass_policy must be one of: fixed, gas_resolution_fraction");
}

[[nodiscard]] FeedbackMode parseFeedbackMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "thermal") {
    return FeedbackMode::kThermal;
  }
  if (lower == "kinetic") {
    return FeedbackMode::kKinetic;
  }
  if (lower == "momentum") {
    return FeedbackMode::kMomentum;
  }
  if (lower == "thermal_kinetic_momentum") {
    return FeedbackMode::kThermalKineticMomentum;
  }
  throw ConfigError(
      "physics.fb_mode must be one of: thermal, kinetic, momentum, thermal_kinetic_momentum");
}

[[nodiscard]] FeedbackVariant parseFeedbackVariant(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") {
    return FeedbackVariant::kNone;
  }
  if (lower == "delayed_cooling") {
    return FeedbackVariant::kDelayedCooling;
  }
  if (lower == "stochastic") {
    return FeedbackVariant::kStochastic;
  }
  throw ConfigError("physics.fb_variant must be one of: none, delayed_cooling, stochastic");
}

[[nodiscard]] UvBackgroundModel parseUvBackgroundModel(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") {
    return UvBackgroundModel::kNone;
  }
  if (lower == "hm12") {
    return UvBackgroundModel::kHm12;
  }
  if (lower == "fg20") {
    return UvBackgroundModel::kFg20;
  }
  if (lower == "treecool_table") {
    return UvBackgroundModel::kTreecoolTable;
  }
  throw ConfigError("key 'physics.uv_background_model': invalid value '" + value + "'");
}

[[nodiscard]] SelfShieldingModel parseSelfShieldingModel(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") {
    return SelfShieldingModel::kNone;
  }
  if (lower == "rahmati13_like") {
    return SelfShieldingModel::kRahmati13Like;
  }
  throw ConfigError("key 'physics.self_shielding_model': invalid value '" + value + "'");
}

[[nodiscard]] CoolingModel parseCoolingModel(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "primordial") {
    return CoolingModel::kPrimordial;
  }
  if (lower == "primordial_metal_line") {
    return CoolingModel::kPrimordialMetalLine;
  }
  throw ConfigError("key 'physics.cooling_model': invalid value '" + value + "'");
}

[[nodiscard]] ThermochemistryBackendMode parseThermochemistryBackendMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "legacy_cooling") return ThermochemistryBackendMode::kLegacyCooling;
  if (lower == "analytic_reference") return ThermochemistryBackendMode::kAnalyticReference;
  if (lower == "equilibrium_primordial") return ThermochemistryBackendMode::kEquilibriumPrimordial;
  if (lower == "fixed_hybrid") return ThermochemistryBackendMode::kFixedHybrid;
  if (lower == "fixed_detailed") return ThermochemistryBackendMode::kFixedDetailed;
  throw ConfigError("physics.thermochemistry_backend must be one of: legacy_cooling, analytic_reference, equilibrium_primordial, fixed_hybrid, fixed_detailed");
}

[[nodiscard]] MetalCoolingModeConfig parseMetalCoolingModeConfig(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") return MetalCoolingModeConfig::kNone;
  if (lower == "fine_structure_cold") return MetalCoolingModeConfig::kFineStructureCold;
  if (lower == "element_table") return MetalCoolingModeConfig::kElementTable;
  if (lower == "element_table_and_fine_structure") return MetalCoolingModeConfig::kElementTableAndFineStructure;
  throw ConfigError("physics.metal_cooling_mode must be one of: none, fine_structure_cold, element_table, element_table_and_fine_structure");
}

[[nodiscard]] MetalCoolingTableRangeConfig parseMetalCoolingTableRangeConfig(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "reject") return MetalCoolingTableRangeConfig::kReject;
  if (lower == "clamp_density_redshift") return MetalCoolingTableRangeConfig::kClampDensityAndRedshift;
  throw ConfigError("physics.metal_cooling_table_range must be one of: reject, clamp_density_redshift");
}

[[nodiscard]] DustModelConfig parseDustModelConfig(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") return DustModelConfig::kNone;
  if (lower == "two_reservoir_v1") return DustModelConfig::kTwoReservoirV1;
  throw ConfigError("physics.dust_model must be one of: none, two_reservoir_v1");
}

[[nodiscard]] IsmRegimeModelConfig parseIsmRegimeModelConfig(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "legacy") return IsmRegimeModelConfig::kLegacy;
  if (lower == "effective_or_resolved") return IsmRegimeModelConfig::kEffectiveOrResolved;
  throw ConfigError("physics.ism_regime_model must be one of: legacy, effective_or_resolved");
}

[[nodiscard]] SgsTurbulenceModelConfig parseSgsTurbulenceModelConfig(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") return SgsTurbulenceModelConfig::kNone;
  if (lower == "budgeted_v1") return SgsTurbulenceModelConfig::kBudgetedV1;
  throw ConfigError("physics.sgs_turbulence_model must be one of: none, budgeted_v1");
}

[[nodiscard]] ChemistryAdaptiveSelectionConfig parseChemistryAdaptiveSelectionConfig(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "off") return ChemistryAdaptiveSelectionConfig::kOff;
  if (lower == "local_v1") return ChemistryAdaptiveSelectionConfig::kLocalV1;
  throw ConfigError("physics.chemistry_adaptive_selection must be one of: off, local_v1");
}

[[nodiscard]] StellarPopulationModelConfig parseStellarPopulationModelConfig(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "legacy_table") return StellarPopulationModelConfig::kLegacyTable;
  if (lower == "chabrier_reference_v1") return StellarPopulationModelConfig::kChabrierReferenceV1;
  if (lower == "chabrier_yield_table_v1") return StellarPopulationModelConfig::kChabrierYieldTableV1;
  throw ConfigError("physics.stellar_population_model must be one of: legacy_table, chabrier_reference_v1, chabrier_yield_table_v1");
}

[[nodiscard]] SupernovaFeedbackModelConfig parseSupernovaFeedbackModelConfig(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "legacy") return SupernovaFeedbackModelConfig::kLegacy;
  if (lower == "mechanical_resolved_unresolved") return SupernovaFeedbackModelConfig::kMechanicalResolvedUnresolved;
  throw ConfigError("physics.fb_supernova_model must be one of: legacy, mechanical_resolved_unresolved");
}

[[nodiscard]] MetalSpeciesMode parseMetalSpeciesMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "total_only") return MetalSpeciesMode::kTotalOnly;
  if (lower == "core_elements") return MetalSpeciesMode::kCoreElements;
  throw ConfigError("physics.metal_species_mode must be one of: total_only, core_elements");
}

[[nodiscard]] MetalDiffusionModel parseMetalDiffusionModel(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") return MetalDiffusionModel::kNone;
  if (lower == "smagorinsky") return MetalDiffusionModel::kSmagorinsky;
  throw ConfigError("physics.metal_diffusion_model must be one of: none, smagorinsky");
}

[[nodiscard]] MetalDiffusionTimeIntegrator parseMetalDiffusionTimeIntegrator(
    const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "explicit_subcycling") return MetalDiffusionTimeIntegrator::kExplicitSubcycling;
  if (lower == "rkl2") return MetalDiffusionTimeIntegrator::kRkl2;
  throw ConfigError(
      "physics.metal_diffusion_time_integrator must be one of: explicit_subcycling, rkl2");
}

[[nodiscard]] IntegratorTimeVariable parseIntegratorTimeVariable(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "scale_factor" || lower == "a") {
    return IntegratorTimeVariable::kScaleFactor;
  }
  if (lower == "ln_a" || lower == "log_scale_factor") {
    return IntegratorTimeVariable::kLogScaleFactor;
  }
  if (lower == "t_code" || lower == "code_time") {
    return IntegratorTimeVariable::kCodeTime;
  }
  if (lower == "t_phys" || lower == "physical_time") {
    return IntegratorTimeVariable::kPhysicalTime;
  }
  throw ConfigError(
      "key 'numerics.integrator_time_variable': invalid value '" + value +
      "' (supported: scale_factor, ln_a, t_code, t_phys)");
}

[[nodiscard]] double redshiftToScaleFactor(double z, const std::string& key) {
  if (!std::isfinite(z) || z <= -1.0) {
    throw ConfigError("key '" + key + "': redshift must be finite and > -1");
  }
  return 1.0 / (1.0 + z);
}

[[nodiscard]] double scaleFactorToRedshift(double a, const std::string& key) {
  if (!std::isfinite(a) || a <= 0.0) {
    throw ConfigError("key '" + key + "': scale factor must be finite and > 0");
  }
  const double z = (1.0 / a) - 1.0;
  // For very large but finite a, IEEE-754 subtraction can round the physical
  // redshift to exactly -1. Preserve the closest representable value in the
  // valid z > -1 domain; a remains the authoritative finite endpoint.
  return z > -1.0 ? z : std::nextafter(-1.0, 0.0);
}

void validateScaleRedshiftPair(double a, double z, const std::string& a_key, const std::string& z_key) {
  if (!std::isfinite(a) || a <= 0.0) {
    throw ConfigError("key '" + a_key + "': scale factor must be finite and > 0");
  }
  if (!std::isfinite(z) || z <= -1.0) {
    throw ConfigError("key '" + z_key + "': redshift must be finite and > -1");
  }
  const double expected_z = scaleFactorToRedshift(a, a_key);
  const double tolerance = 16.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(z), std::abs(expected_z)});
  if (std::abs(expected_z - z) > tolerance) {
    throw ConfigError("keys '" + a_key + "' and '" + z_key +
                      "' are inconsistent: expected " + z_key + " = 1/" + a_key + " - 1");
  }
}

void normalizeScaleRedshiftPair(
    std::map<std::string, ParsedEntry>& entries,
    std::set<std::string>& consumed,
    const std::string& a_key,
    const std::string& z_key,
    double default_a,
    double default_z,
    double& out_a,
    double& out_z) {
  const bool has_a = entries.contains(a_key);
  const bool has_z = entries.contains(z_key);

  if (has_a && has_z) {
    out_a = parseFloating(requireString(entries, consumed, a_key, std::to_string(default_a)), a_key);
    out_z = parseFloating(requireString(entries, consumed, z_key, std::to_string(default_z)), z_key);
    validateScaleRedshiftPair(out_a, out_z, a_key, z_key);
    return;
  }

  if (has_z) {
    out_z = parseFloating(requireString(entries, consumed, z_key, std::to_string(default_z)), z_key);
    out_a = redshiftToScaleFactor(out_z, z_key);
    return;
  }

  if (has_a) {
    out_a = parseFloating(requireString(entries, consumed, a_key, std::to_string(default_a)), a_key);
    out_z = scaleFactorToRedshift(out_a, a_key);
    return;
  }

  out_a = default_a;
  out_z = default_z;
  validateScaleRedshiftPair(out_a, out_z, a_key, z_key);
}

[[nodiscard]] AnalysisConfig::DiagnosticsExecutionPolicy parseDiagnosticsExecutionPolicy(
    const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "run_health_only") {
    return AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthOnly;
  }
  if (lower == "run_health_and_light_science") {
    return AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthAndLightScience;
  }
  if (lower == "all_including_provisional") {
    return AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
  }
  throw ConfigError(
      "analysis.diagnostics_execution_policy must be one of: run_health_only, "
      "run_health_and_light_science, all_including_provisional");
}

[[nodiscard]] std::string diagnosticsExecutionPolicyToString(
    AnalysisConfig::DiagnosticsExecutionPolicy policy) {
  switch (policy) {
    case AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthOnly:
      return "run_health_only";
    case AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthAndLightScience:
      return "run_health_and_light_science";
    case AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional:
      return "all_including_provisional";
  }
  throw ConfigError("unhandled DiagnosticsExecutionPolicy enum value during serialization");
}

[[nodiscard]] std::string pmDecompositionModeToString(PmDecompositionMode mode) {
  switch (mode) {
    case PmDecompositionMode::kSlab:
      return "slab";
    case PmDecompositionMode::kPencil:
      return "pencil";
  }
  throw ConfigError("unhandled PmDecompositionMode enum value during serialization");
}

struct ConfigKeySpec {
  const char* key;
  const char* default_value;
};

[[nodiscard]] const std::vector<ConfigKeySpec>& configKeyRegistry() {
  static const std::vector<ConfigKeySpec> keys = {
      {"schema_version", "1"},
      {"units.length_unit", "mpc"},
      {"units.mass_unit", "msun"},
      {"units.velocity_unit", "km_s"},
      {"units.coordinate_frame", "comoving"},
      {"mode.mode", "zoom_in"},
      {"mode.ic_file", "generated"},
      {"mode.ic_convention", "generated"},
      {"mode.ic_manifest_file", ""},
      {"mode.ic_bridge_source_length_unit_to_si", "0"},
      {"mode.ic_bridge_source_mass_unit_to_si", "0"},
      {"mode.ic_bridge_source_velocity_unit_to_si", "0"},
      {"mode.ic_bridge_coordinate_frame", "unspecified"},
      {"mode.ic_bridge_velocity_convention", "unspecified"},
      {"mode.ic_bridge_length_hubble_exponent", "unspecified"},
      {"mode.ic_bridge_length_scale_factor_exponent", "unspecified"},
      {"mode.ic_bridge_mass_hubble_exponent", "unspecified"},
      {"mode.ic_bridge_mass_scale_factor_exponent", "unspecified"},
      {"mode.ic_bridge_velocity_hubble_exponent", "unspecified"},
      {"mode.ic_bridge_velocity_scale_factor_exponent", "unspecified"},
      {"mode.ic_chunk_particle_count", "65536"},
      {"mode.ic_staging_particle_count", "65536"},
      {"mode.ic_part_type2_policy", "reject"},
      {"mode.ic_part_type3_policy", "reject"},
      {"mode.ic_gas_internal_energy_policy", "reject"},
      {"mode.ic_gas_density_policy", "reject"},
      {"mode.ic_star_formation_time_policy", "reject"},
      {"mode.ic_star_initial_mass_policy", "reject"},
      {"mode.ic_star_metallicity_policy", "reject"},
      {"mode.ic_bh_mdot_policy", "reject"},
      {"mode.ic_gas_internal_energy_value_code", "0"},
      {"mode.ic_gas_density_value_code", "0"},
      {"mode.ic_star_formation_time_value", "0"},
      {"mode.ic_star_initial_mass_value_code", "0"},
      {"mode.ic_star_metallicity_value", "0"},
      {"mode.ic_bh_mdot_value_code", "0"},
      {"mode.zoom_high_res_region", "false"},
      {"mode.zoom_region_file", ""},
      {"mode.zoom_long_range_strategy", "disabled"},
      {"mode.zoom_region_center_x", "0.5"},
      {"mode.zoom_region_center_y", "0.5"},
      {"mode.zoom_region_center_z", "0.5"},
      {"mode.zoom_region_radius", "0.0"},
      {"mode.zoom_focused_pm_grid_nx", "0"},
      {"mode.zoom_focused_pm_grid_ny", "0"},
      {"mode.zoom_focused_pm_grid_nz", "0"},
      {"mode.zoom_contamination_radius", "0.0"},
      {"mode.hydro_boundary", "auto"},
      {"mode.gravity_boundary", "auto"},
      {"cosmology.omega_matter", "0.315"},
      {"cosmology.omega_lambda", "0.685"},
      {"cosmology.omega_baryon", "0.049"},
      {"cosmology.hubble_param", "0.674"},
      {"cosmology.sigma8", "0.811"},
      {"cosmology.scalar_index_ns", "0.965"},
      {"cosmology.box_size_x", "50.0"},
      {"cosmology.box_size_y", "50.0"},
      {"cosmology.box_size_z", "50.0"},
      {"cosmology.box_size", "50.0"},
      {"numerics.a_begin", "1.0"},
      {"numerics.a_end", "1.0"},
      {"numerics.z_begin", "0.0"},
      {"numerics.z_end", "0.0"},
      {"numerics.t_code_begin", "0.0"},
      {"numerics.t_code_end", "1.0"},
      {"numerics.t_phys_begin", "0.0"},
      {"numerics.t_phys_end", "0.0"},
      {"numerics.integrator_time_variable", "scale_factor"},
      {"numerics.cosmology_max_delta_ln_a", "1.0e-2"},
      {"numerics.cosmology_max_hubble_time_fraction", "1.0e-2"},
      {"numerics.source_max_fractional_change", "0.1"},
      {"numerics.hydro_density_floor_code", "1.0e-10"},
      {"numerics.hydro_pressure_floor_code", "1.0e-10"},
      {"numerics.hydro_active_batch_max_cells", "0"},
      {"numerics.hydro_reconstruction", "piecewise_linear"},
      {"numerics.hydro_slope_limiter", "mc"},
      {"numerics.hydro_integration_scheme", "muscl_hancock"},
      {"numerics.hydro_riemann_solver", "hllc"},
      {"numerics.hydro_dual_energy", "off"},
      {"numerics.hydro_dual_energy_eta", "1.0e-3"},
      {"numerics.max_global_steps", "1024"},
      {"numerics.hierarchical_max_rung", "0"},
      {"numerics.hierarchical_coarse_displacement_criterion", "max_speed_mesh"},
      {"numerics.hierarchical_coarse_rms_displacement_factor", "0.25"},
      {"numerics.hierarchical_encounter_timestep_enabled", "true"},
      {"numerics.gravity_timestep_eta", "0.2"},
      {"numerics.amr_max_level", "10"},
      {"numerics.amr_regrid_mode", "off"},
      {"numerics.global_step_max_retries", "0"},
      {"numerics.amr_time_stepping", "synchronized"},
      {"numerics.amr_predictive_feedback_age_years", "0.0"},
      {"numerics.amr_predictive_event_lookahead_years", "0.0"},
      {"numerics.amr_regrid_interval_steps", "1"},
      {"numerics.amr_regrid_max_operations", "64"},
      {"numerics.amr_refine_cell_mass_code", "0.0"},
      {"numerics.amr_refine_jump_threshold", "0.0"},
      {"numerics.amr_refine_jeans_cells", "4.0"},
      {"numerics.amr_derefine_hysteresis", "0.5"},
      {"numerics.amr_refine_cooling_length_cells", "0.0"},
      {"numerics.amr_refine_metallicity_jump", "0.0"},
      {"numerics.amr_refine_species_jump", "0.0"},
      {"numerics.amr_refine_vorticity", "0.0"},
      {"numerics.gravity_softening", "1.0 kpc"},
      {"numerics.gravity_softening_gas", ""},
      {"numerics.gravity_softening_dark_matter", ""},
      {"numerics.gravity_softening_star", ""},
      {"numerics.gravity_softening_black_hole", ""},
      {"numerics.gravity_softening_tracer", ""},
      {"numerics.gravity_solver", "treepm"},
      {"numerics.hydro_solver", "godunov_fv"},
      {"numerics.treepm_pm_grid_nx", "16"},
      {"numerics.treepm_pm_grid_ny", "16"},
      {"numerics.treepm_pm_grid_nz", "16"},
      {"numerics.treepm_pm_grid", "16"},
      {"numerics.treepm_asmth_cells", "1.25"},
      {"numerics.treepm_rcut_cells", "6.25"},
      {"numerics.treepm_tree_opening_criterion", "com_distance"},
      {"numerics.treepm_tree_opening_theta", "0.7"},
      {"numerics.treepm_tree_relative_force_tolerance", "0.005"},
      {"numerics.treepm_tree_relative_force_acceleration_floor", "1.0e-30"},
      {"numerics.treepm_gaussian_pair_lookup_enabled", "false"},
      {"numerics.treepm_full_residual_lookup_enabled", "false"},
      {"numerics.treepm_leaf_simd_enabled", "false"},
      {"numerics.treepm_radial_coefficients_enabled", "false"},
      {"numerics.treepm_residual_block_size", "64"},
      {"numerics.treepm_tree_max_leaf_size", "16"},

      {"numerics.treepm_full_mac_diagnostics", "true"},
      {"numerics.treepm_adaptive_acceptance_enabled", "false"},
      {"numerics.treepm_adaptive_maximum_opening_angle", "0.25"},
      {"numerics.treepm_adaptive_error_estimator", "residual_derivative"},
      {"numerics.treepm_cutoff_node_policy", "contain"},
      {"numerics.treepm_group_walk_enabled", "false"},
      {"numerics.treepm_group_walk_target_size", "0"},
      {"numerics.treepm_identical_source_tree_reuse_enabled", "false"},
      {"numerics.treepm_topology_refit_enabled", "false"},
      {"numerics.treepm_assignment_scheme", "tsc"},
      {"numerics.treepm_pm_force_gradient", "spectral"},
      {"numerics.treepm_enable_window_deconvolution", "true"},
      {"numerics.treepm_update_cadence_steps", "1"},
      {"numerics.treepm_pm_reuse_policy", "exact"},
      {"numerics.treepm_multirate_rollback_memory_policy", "single_block_fallback"},
      {"numerics.treepm_allow_diagnostic_naive_dft", "false"},
      {"numerics.treepm_pm_decomposition_mode", "slab"},
      {"numerics.treepm_tree_exchange_batch_bytes", "4194304"},
      {"numerics.treepm_pm_exchange_batch_bytes", "16777216"},
      {"physics.enable_cooling", "true"},
      {"physics.enable_star_formation", "true"},
      {"physics.enable_feedback", "true"},
      {"physics.enable_stellar_evolution", "true"},
      {"physics.reionization_model", "hm12"},
      {"physics.uv_background_model", "hm12"},
      {"physics.self_shielding_model", "none"},
      {"physics.cooling_model", "primordial"},
      {"physics.metal_line_table_path", ""},
      {"physics.temperature_floor_k", "100.0"},
      {"physics.star_formation_model", "legacy_schmidt_threshold"},
      {"physics.sf_density_threshold_code", "10.0"},
      {"physics.sf_temperature_threshold_k", "1.0e4"},
      {"physics.sf_min_converging_flow_rate_code", "0.0"},
      {"physics.sf_epsilon_ff", "0.01"},
      {"physics.sf_bound_alpha_vir_max", "1.0"},
      {"physics.sf_require_converging_flow", "true"},
      {"physics.sf_collapse_timescale", "free_fall"},
      {"physics.sf_jeans_mass_floor_code", "100.0"},
      {"physics.sf_star_particle_mass_policy", "fixed"},
      {"physics.sf_target_star_particle_mass_code", "0.1"},
      {"physics.sf_target_star_particle_mass_fraction", "0.25"},
      {"physics.sf_min_star_particle_mass_code", "0.1"},
      {"physics.sf_max_star_particle_mass_code", "1.0e30"},
      {"physics.sf_max_spawn_particles_per_cell_step", "8"},
      {"physics.sf_max_fractional_mass_conversion", "0.25"},
      {"physics.sf_min_remaining_gas_fraction", "0.01"},
      {"physics.sf_min_remaining_gas_mass_code", "1.0e-12"},
      {"physics.sf_temperature_safety_ceiling_k", "0.0"},
      {"physics.sf_stochastic_spawning", "true"},
      {"physics.sf_random_seed", "123456789"},
      {"physics.sf_effective_parameter_set", "chui_sh03_tng_like_v1"},
      {"physics.sf_effective_n_h_threshold", "0.13"},
      {"physics.sf_effective_min_baryon_overdensity", "0.0"},
      {"physics.sf_effective_t_star_at_threshold", "2.2 Gyr"},
      {"physics.sf_effective_evaporation_factor_at_threshold", "1000.0"},
      {"physics.sf_effective_evaporation_density_exponent", "-0.8"},
      {"physics.sf_effective_cold_phase_temperature", "1000.0"},
      {"physics.sf_effective_supernova_specific_energy_code", "1.0e5"},
      {"physics.sf_effective_massive_star_fraction", "0.1"},
      {"physics.sf_effective_q_eos", "0.3"},
      {"physics.sf_effective_isothermal_temperature", "1.0e4"},
      {"physics.sf_effective_hot_excess_tolerance", "1.0"},
      {"physics.sf_effective_eos_relaxation", "instantaneous"},
      {"physics.sf_effective_eos_relaxation_timescale", "0.0 code"},
      {"physics.sf_effective_eos_table_bins", "256"},
      {"physics.sf_effective_eos_max_density_ratio", "1.0e6"},
      {"physics.sf_effective_birth_mass_convention", "initial_ssp_mass"},
      {"physics.sf_effective_feedback_coupling", "external_feedback_calibrated"},
      {"physics.fb_mode", "thermal_kinetic_momentum"},
      {"physics.fb_variant", "none"},
      {"physics.fb_use_returned_mass_budget", "true"},
      {"physics.fb_epsilon_thermal", "0.6"},
      {"physics.fb_epsilon_kinetic", "0.3"},
      {"physics.fb_epsilon_momentum", "0.1"},
      {"physics.fb_sn_energy_erg_per_mass_code", "1.0e49"},
      {"physics.fb_momentum_code_per_mass_code", "3.0e3"},
      {"physics.fb_neighbor_count", "8"},
      {"physics.fb_delayed_cooling_time_code", "0.0"},
      {"physics.fb_stochastic_event_probability", "0.25"},
      {"physics.fb_random_seed", "42424242"},
      {"physics.stellar_evolution_table_path", ""},
      {"physics.stellar_evolution_hubble_time_years", "1.44e10"},
      {"physics.stellar_evolution_require_production_table", "false"},
      {"physics.metal_species_mode", "total_only"},
      {"physics.enable_metal_diffusion", "false"},
      {"physics.metal_diffusion_model", "none"},
      {"physics.metal_diffusion_time_integrator", "explicit_subcycling"},
      {"physics.metal_diffusion_coefficient", "0.05"},
      {"physics.metal_diffusion_cfl", "0.4"},
      {"physics.metal_diffusion_max_subcycles", "128"},
      {"physics.metal_diffusion_max_rkl_stages", "64"},
      {"physics.metal_diffusion_coefficient_floor_code", "0.0"},
      {"physics.metal_diffusion_coefficient_ceiling_code", "1.0e30"},
      {"physics.enable_black_hole_agn", "false"},
      {"physics.bh_enable_seeding", "false"},
      {"physics.bh_seed_halo_mass_threshold_code", "1.0e3"},
      {"physics.bh_seed_mass_code", "1.0"},
      {"physics.bh_seed_max_per_cell", "1"},
      {"physics.bh_alpha_bondi", "1.0"},
      {"physics.bh_use_eddington_cap", "true"},
      {"physics.bh_epsilon_r", "0.1"},
      {"physics.bh_epsilon_f", "0.05"},
      {"physics.bh_feedback_coupling_efficiency", "1.0"},
      {"physics.bh_duty_cycle_active_edd_ratio_threshold", "0.01"},
      {"physics.bh_proton_mass_si", "1.67262192369e-27"},
      {"physics.bh_thomson_cross_section_si", "6.6524587321e-29"},
      {"physics.bh_newton_g_si", "6.67430e-11"},
      {"physics.bh_speed_of_light_si", "2.99792458e8"},
      {"physics.enable_tracers", "false"},
      {"physics.tracer_track_mass", "true"},
      {"physics.tracer_min_host_mass_code", "0.0"},
      {"physics.thermochemistry_backend", "legacy_cooling"},
      {"physics.metal_cooling_mode", "none"},
      {"physics.metal_cooling_table_range", "reject"},
      {"physics.dust_model", "none"},
      {"physics.ism_regime_model", "legacy"},
      {"physics.sgs_turbulence_model", "none"},
      {"physics.chemistry_adaptive_selection", "off"},
      {"physics.stellar_population_model", "legacy_table"},
      {"physics.fb_supernova_model", "legacy"},
      {"physics.thermochemistry_relative_tolerance", "1.0e-4"},
      {"physics.thermochemistry_abundance_absolute_tolerance", "1.0e-12"},
      {"physics.thermochemistry_analytic_relaxation_time_years", "1.0e6"},
      {"physics.thermochemistry_analytic_equilibrium_temperature_k", "1.0e4"},
      {"physics.isrf_g0_habing", "0.0"},
      {"physics.chemistry_shielding_length_max_pc", "100.0"},
      {"physics.metal_cooling_switch_temperature_k", "1.0e4"},
      {"physics.chemistry_primordial_helium_mass_fraction", "0.2449"},
      {"physics.dust_accretion_timescale_ref_years", "2.0e7"},
      {"physics.dust_sputtering_timescale_ref_years", "1.7e8"},
      {"physics.dust_sn_swept_mass_msun", "6800.0"},
      {"physics.dust_sn_destruction_efficiency", "0.3"},
      {"physics.dust_agb_condensation_efficiency", "0.2"},
      {"physics.dust_cc_condensation_efficiency", "0.1"},
      {"physics.ism_jeans_resolution_cells", "4.0"},
      {"physics.ism_regime_hysteresis_fraction", "0.2"},
      {"physics.ism_regime_min_dwell_years", "1.0e6"},
      {"physics.sgs_dissipation_coefficient", "1.58"},
      {"physics.sgs_eddy_viscosity_coefficient", "0.05"},
      {"physics.sgs_feedback_energy_fraction", "0.0"},
      {"physics.chemistry_adaptive_hysteresis_fraction", "0.3"},
      {"physics.chemistry_adaptive_min_dwell_years", "1.0e6"},
      {"physics.chemistry_adaptive_cold_temperature_k", "3.0e3"},
      {"physics.chemistry_adaptive_cooling_time_ratio", "0.0"},
      {"physics.chemistry_adaptive_dense_nh_cm3", "1.0"},
      {"physics.chemistry_shadow_sample_fraction", "0.0"},
      {"physics.ssp_imf_mass_min_msun", "0.1"},
      {"physics.ssp_imf_mass_max_msun", "100.0"},
      {"physics.ssp_cc_mass_min_msun", "8.0"},
      {"physics.ssp_cc_mass_max_msun", "40.0"},
      {"physics.ssp_snia_per_msun", "1.3e-3"},
      {"physics.ssp_snia_delay_min_years", "4.0e7"},
      {"physics.ssp_snia_delay_max_years", "1.37e10"},
      {"physics.ssp_energy_per_sn_erg", "1.0e51"},
      {"physics.ssp_reference_net_metal_yield", "0.0"},
      {"physics.fb_presn_wind_momentum_lsun_per_c", "0.0"},
      {"physics.thermochemistry_max_substeps", "4096"},
      {"physics.thermochemistry_max_newton_iterations", "12"},
      {"physics.chemistry_shadow_max_cells", "0"},
      {"physics.uv_background_table_path", ""},
      {"physics.uv_background_table_sha256", ""},
      {"physics.metal_cooling_table_path", ""},
      {"physics.metal_cooling_table_sha256", ""},
      {"physics.metal_cooling_accept_uv_background_mismatch", "false"},
      {"physics.ssp_yield_table_path", ""},
      {"physics.ssp_yield_table_sha256", ""},
      {"physics.chemistry_metallicity_scaled_dust_proxy", "false"},
      {"physics.thermochemistry_allow_unqualified_rates", "false"},
      {"physics.ssp_stochastic_sn_events", "true"},
      {"physics.ssp_clustered_associations", "false"},
      {"physics.ssp_association_radius_code", "0.0"},
      {"physics.ssp_association_time_window_years", "3.0e6"},
      {"physics.ssp_wind_model", "none"},
      {"physics.ssp_wind_table_path", ""},
      {"physics.ssp_wind_table_sha256", ""},
      {"output.run_name", "cosmosim_run"},
      {"output.output_directory", "outputs"},
      {"output.output_stem", "snapshot"},
      {"output.restart_stem", "restart"},
      {"output.snapshot_interval_steps", "64"},
      {"output.snapshot_interval_time_code", "0.0"},
      {"output.snapshot_layout", "auto"},
      {"output.snapshot_num_files", "0"},
      {"output.write_restarts", "true"},
      {"output.restart_interval_steps", "0"},
      {"output.restart_interval_time_code", "0.0"},
      {"output.snapshot_verification", "full_readback"},
      {"parallel.mpi_ranks_expected", "1"},
      {"parallel.omp_threads", "1"},
      {"parallel.gpu_devices", "0"},
      {"parallel.deterministic_reduction", "true"},
      {"parallel.decomposition_particle_count_weight", "1.0"},
      {"parallel.decomposition_gas_cell_weight", "1.5"},
      {"parallel.decomposition_tree_interaction_weight", "1.0"},
      {"parallel.decomposition_pm_mesh_weight", "0.25"},
      {"parallel.decomposition_amr_patch_weight", "1.0"},
      {"parallel.decomposition_active_fraction_weight", "2.0"},
      {"parallel.decomposition_memory_pressure_weight", "9.5367431640625e-7"},
      {"parallel.decomposition_gpu_occupancy_weight", "0.0"},
      {"parallel.decomposition_generic_work_weight", "0.5"},
      {"parallel.decomposition_runtime_rebalance_enabled", "true"},
      {"parallel.decomposition_spatial_work_enabled", "false"},
      {"parallel.decomposition_debug_exact_ownership_audit", "false"},
      {"parallel.decomposition_rebalance_imbalance_trigger", "1.25"},
      {"parallel.decomposition_rebalance_memory_trigger", "1.50"},
      {"parallel.decomposition_rebalance_max_migrated_load_fraction", "0.25"},
      {"parallel.decomposition_measured_tree_pair_weight", "1.0"},
      {"parallel.decomposition_measured_pm_cell_weight", "1.0"},
      {"parallel.decomposition_measured_amr_cell_weight", "1.0"},
      {"parallel.decomposition_measured_hydro_face_weight", "1.0"},
      {"parallel.decomposition_measured_wall_ms_weight", "1.0"},
      {"parallel.isolated_pm_root_workspace_limit_bytes", "268435456"},
      {"parallel.zoom_high_res_allgather_limit_bytes", "268435456"},
      {"parallel.gravity_memory_budget_bytes", "0"},
      {"parallel.gravity_backend_unknown_reserve_bytes", "0"},
      {"parallel.gravity_memory_safety_margin_fraction", "0.0"},
      {"parallel.process_memory_budget_bytes", "0"},
      {"parallel.process_mpi_unknown_reserve_bytes", "0"},
      {"parallel.process_fftw_unknown_reserve_bytes", "0"},
      {"parallel.process_hdf5_unknown_reserve_bytes", "0"},
      {"parallel.process_allocator_unknown_reserve_bytes", "0"},
      {"parallel.process_output_restart_overlap_bytes", "0"},
      {"parallel.process_memory_safety_margin_fraction", "0.0"},
      {"parallel.node_memory_admission", "off"},
      {"parallel.node_memory_budget_bytes", "0"},
      {"parallel.node_memory_reserve_bytes", "0"},
      {"analysis.enable_diagnostics", "true"},
      {"analysis.enable_halo_workflow", "false"},
      {"analysis.halo_on_the_fly", "false"},
      {"analysis.diagnostics_execution_policy", "run_health_and_light_science"},
      {"analysis.run_health_interval_steps", "1"},
      {"analysis.science_light_interval_steps", "8"},
      {"analysis.science_heavy_interval_steps", "64"},
      {"analysis.retention_bundle_count", "8"},
      {"analysis.power_spectrum_mesh_n", "16"},
      {"analysis.power_spectrum_bin_count", "12"},
      {"analysis.sf_history_bin_count", "16"},
      {"analysis.quicklook_grid_n", "32"},
      {"analysis.diagnostics_stem", "diagnostics"},
      {"analysis.halo_catalog_stem", "halo_catalog"},
      {"analysis.merger_tree_stem", "merger_tree_plan"},
      {"analysis.halo_fof_linking_length_factor", "0.2"},
      {"analysis.halo_fof_min_group_size", "16"},
      {"analysis.halo_include_gas", "true"},
      {"analysis.halo_include_stars", "true"},
      {"analysis.halo_include_black_holes", "true"},
      {"compatibility.allow_unknown_keys", "false"},
  };
  return keys;
}

[[nodiscard]] const std::map<std::string, std::string>& deprecatedAliasRegistry() {
  static const std::map<std::string, std::string> aliases = {
      {"omega0", "cosmology.omega_matter"},
      {"omegalambda", "cosmology.omega_lambda"},
      {"hubbleparam", "cosmology.hubble_param"},
      {"timemax", "numerics.t_code_end"},
      {"mode", "mode.mode"},
      {"run_name", "output.run_name"},
      {"numerics.initial_scale_factor", "numerics.a_begin"},
      {"numerics.initial_redshift", "numerics.z_begin"},
      {"numerics.time_begin_code", "numerics.t_code_begin"},
      {"numerics.time_end_code", "numerics.t_code_end"},
  };
  return aliases;
}

[[nodiscard]] std::string defaultFor(const std::string& key) {
  for (const auto& entry : configKeyRegistry()) {
    if (entry.key == key) {
      return entry.default_value;
    }
  }
  throw ConfigError("internal config key registry missing default for key '" + key + "'");
}

// CHUI Part 2 typed validation of the baryonic-physics keys. Every rule is
// explicit; no Part 2 model silently degrades to a cheaper one.
void validatePart2BaryonicPhysics(const SimulationConfig& config) {
  const PhysicsConfig& p = config.physics;
  const bool part2 = p.thermochemistry_backend != ThermochemistryBackendMode::kLegacyCooling;
  if (part2 != (p.metal_species_mode == MetalSpeciesMode::kCoreElements)) {
    throw ConfigError(
        "physics.metal_species_mode=core_elements is required by (and only valid with) a Part 2 "
        "thermochemistry_backend (analytic_reference, equilibrium_primordial, fixed_hybrid, fixed_detailed)");
  }
  if ((p.thermochemistry_backend == ThermochemistryBackendMode::kEquilibriumPrimordial ||
       p.thermochemistry_backend == ThermochemistryBackendMode::kFixedHybrid ||
       p.thermochemistry_backend == ThermochemistryBackendMode::kFixedDetailed) &&
      !p.thermochemistry_allow_unqualified_rates) {
    throw ConfigError(
        "thermochemistry_backend " + thermochemistryBackendModeToString(p.thermochemistry_backend) +
        " uses rate coefficients not yet verified against their publications (UNQUALIFIED; see "
        "docs/chemistry_rate_provenance.md). Set physics.thermochemistry_allow_unqualified_rates = true to run it "
        "as an explicitly unqualified research model");
  }
  // CHUI Part 3: with an element-resolved composition sidecar,
  // physics.enable_metal_diffusion selects the element-resolved mixing
  // operator (physics/element_diffusion.hpp); the total-Z mirror is derived
  // from the mixed lanes. The Part 2 refusal is lifted.
  if (part2 && !p.enable_cooling) {
    throw ConfigError("a Part 2 thermochemistry_backend requires physics.enable_cooling=true");
  }
  if (part2 && (p.uv_background_model == UvBackgroundModel::kHm12 ||
                p.uv_background_model == UvBackgroundModel::kFg20)) {
    throw ConfigError(
        "physics.uv_background_model=hm12/fg20 are legacy analytic proxies; Part 2 backends require "
        "uv_background_model=none or treecool_table with uv_background_table_path");
  }
  if (p.uv_background_model == UvBackgroundModel::kTreecoolTable) {
    if (!part2) throw ConfigError("uv_background_model=treecool_table requires a Part 2 thermochemistry_backend");
    if (p.uv_background_table_path.empty()) {
      throw ConfigError("uv_background_model=treecool_table requires physics.uv_background_table_path");
    }
  }
  if (!part2 && (p.metal_cooling_mode != MetalCoolingModeConfig::kNone ||
                 p.dust_model != DustModelConfig::kNone ||
                 p.ism_regime_model != IsmRegimeModelConfig::kLegacy ||
                 p.sgs_turbulence_model != SgsTurbulenceModelConfig::kNone ||
                 p.chemistry_adaptive_selection != ChemistryAdaptiveSelectionConfig::kOff)) {
    throw ConfigError(
        "metal_cooling_mode, dust_model, ism_regime_model, sgs_turbulence_model and "
        "chemistry_adaptive_selection require a Part 2 thermochemistry_backend");
  }
  const bool network = p.thermochemistry_backend == ThermochemistryBackendMode::kFixedHybrid ||
      p.thermochemistry_backend == ThermochemistryBackendMode::kFixedDetailed;
  if ((p.metal_cooling_mode != MetalCoolingModeConfig::kNone || p.dust_model != DustModelConfig::kNone) &&
      !network) {
    throw ConfigError("metal_cooling_mode and dust_model require thermochemistry_backend=fixed_hybrid or fixed_detailed");
  }
  if ((p.metal_cooling_mode == MetalCoolingModeConfig::kElementTable ||
       p.metal_cooling_mode == MetalCoolingModeConfig::kElementTableAndFineStructure) &&
      p.metal_cooling_table_path.empty()) {
    throw ConfigError("element-table metal cooling requires physics.metal_cooling_table_path (no table is bundled)");
  }
  if (p.metal_cooling_accept_uv_background_mismatch &&
      p.metal_cooling_mode != MetalCoolingModeConfig::kElementTable &&
      p.metal_cooling_mode != MetalCoolingModeConfig::kElementTableAndFineStructure) {
    throw ConfigError(
        "physics.metal_cooling_accept_uv_background_mismatch = true has no effect without an element-table "
        "metal_cooling_mode (no-op options are refused)");
  }
  if (p.chemistry_adaptive_selection == ChemistryAdaptiveSelectionConfig::kLocalV1 &&
      p.thermochemistry_backend != ThermochemistryBackendMode::kFixedDetailed) {
    throw ConfigError(
        "chemistry_adaptive_selection=local_v1 needs the detailed species lanes (thermochemistry_backend=fixed_detailed)");
  }
  if (!(p.thermochemistry_relative_tolerance > 0.0) || p.thermochemistry_relative_tolerance >= 0.1) {
    throw ConfigError("physics.thermochemistry_relative_tolerance must be in (0, 0.1)");
  }
  if (!(p.thermochemistry_abundance_absolute_tolerance > 0.0) || p.thermochemistry_abundance_absolute_tolerance > 1.0e-4) {
    throw ConfigError("physics.thermochemistry_abundance_absolute_tolerance must be in (0, 1e-4]");
  }
  if (p.thermochemistry_max_substeps == 0U || p.thermochemistry_max_newton_iterations == 0U) {
    throw ConfigError("thermochemistry_max_substeps and thermochemistry_max_newton_iterations must be > 0");
  }
  if (!(p.thermochemistry_analytic_relaxation_time_years > 0.0) ||
      !(p.thermochemistry_analytic_equilibrium_temperature_k > 0.0)) {
    throw ConfigError("analytic thermochemistry reference requires positive relaxation time and temperature");
  }
  if (p.isrf_g0_habing < 0.0 || !(p.chemistry_shielding_length_max_pc > 0.0) ||
      !(p.metal_cooling_switch_temperature_k > 10.0)) {
    throw ConfigError("isrf_g0_habing >= 0, chemistry_shielding_length_max_pc > 0, metal_cooling_switch_temperature_k > 10 required");
  }
  if (!(p.chemistry_primordial_helium_mass_fraction > 0.2) || p.chemistry_primordial_helium_mass_fraction > 0.3) {
    throw ConfigError("physics.chemistry_primordial_helium_mass_fraction must be in (0.2, 0.3]");
  }
  for (const double v : {p.dust_accretion_timescale_ref_years, p.dust_sputtering_timescale_ref_years,
                         p.dust_sn_swept_mass_msun}) {
    if (!(v > 0.0)) throw ConfigError("dust timescales and SN swept mass must be > 0");
  }
  for (const double v : {p.dust_sn_destruction_efficiency, p.dust_agb_condensation_efficiency,
                         p.dust_cc_condensation_efficiency, p.sgs_feedback_energy_fraction,
                         p.chemistry_shadow_sample_fraction}) {
    if (!(v >= 0.0) || v > 1.0) throw ConfigError("dust efficiencies, sgs_feedback_energy_fraction and chemistry_shadow_sample_fraction must be in [0, 1]");
  }
  if (!(p.ism_jeans_resolution_cells >= 1.0) || !(p.ism_regime_hysteresis_fraction >= 0.0) ||
      p.ism_regime_hysteresis_fraction >= 1.0 || p.ism_regime_min_dwell_years < 0.0) {
    throw ConfigError("ism_jeans_resolution_cells >= 1, ism_regime_hysteresis_fraction in [0,1), dwell >= 0 required");
  }
  if (part2 && p.ism_regime_model == IsmRegimeModelConfig::kLegacy &&
      p.star_formation_model == StarFormationModelKind::kEffectiveMultiphaseTngLike) {
    throw ConfigError(
        "a Part 2 thermochemistry_backend with the effective multiphase ISM requires "
        "ism_regime_model=effective_or_resolved, so effective-regime cells are excluded from resolved "
        "thermochemistry (no double cooling)");
  }
  if (p.ism_regime_model == IsmRegimeModelConfig::kEffectiveOrResolved &&
      p.star_formation_model != StarFormationModelKind::kEffectiveMultiphaseTngLike) {
    throw ConfigError(
        "ism_regime_model=effective_or_resolved requires star_formation_model=effective_multiphase_tng_like "
        "(the effective closure is the unresolved branch)");
  }
  if (!(p.sgs_dissipation_coefficient > 0.0) || !(p.sgs_eddy_viscosity_coefficient > 0.0)) {
    throw ConfigError("SGS closure coefficients must be > 0");
  }
  if (!(p.chemistry_adaptive_hysteresis_fraction >= 0.0) || p.chemistry_adaptive_hysteresis_fraction >= 1.0 ||
      p.chemistry_adaptive_min_dwell_years < 0.0 || !(p.chemistry_adaptive_cold_temperature_k > 10.0) ||
      !(p.chemistry_adaptive_dense_nh_cm3 > 0.0) || !(p.chemistry_adaptive_cooling_time_ratio >= 0.0)) {
    throw ConfigError("invalid chemistry_adaptive_* parameters");
  }
  if (!(p.ssp_imf_mass_min_msun > 0.0) || !(p.ssp_imf_mass_max_msun > p.ssp_imf_mass_min_msun) ||
      !(p.ssp_cc_mass_min_msun > p.ssp_imf_mass_min_msun) || !(p.ssp_cc_mass_max_msun > p.ssp_cc_mass_min_msun) ||
      p.ssp_cc_mass_max_msun > p.ssp_imf_mass_max_msun) {
    throw ConfigError("SSP mass limits must satisfy 0 < m_min < m_cc,min < m_cc,max <= m_max");
  }
  if (!(p.ssp_snia_per_msun >= 0.0) || !(p.ssp_snia_delay_min_years > 0.0) ||
      !(p.ssp_snia_delay_max_years > p.ssp_snia_delay_min_years) || !(p.ssp_energy_per_sn_erg > 0.0) ||
      p.ssp_reference_net_metal_yield < 0.0 || p.ssp_reference_net_metal_yield > 0.1 ||
      p.fb_presn_wind_momentum_lsun_per_c < 0.0) {
    throw ConfigError("invalid SSP Ia/energy/yield or pre-SN parameters");
  }
  if (p.stellar_population_model != StellarPopulationModelConfig::kLegacyTable && !part2) {
    throw ConfigError("Part 2 stellar_population_model values require a Part 2 thermochemistry_backend (element composition)");
  }
  if (p.stellar_population_model == StellarPopulationModelConfig::kChabrierYieldTableV1 &&
      p.ssp_yield_table_path.empty()) {
    throw ConfigError("stellar_population_model=chabrier_yield_table_v1 requires physics.ssp_yield_table_path (no yields are bundled)");
  }
  // Part 2.5 reproducibility repair: an external physics data file is bound
  // to the run by its SHA-256. Without a configured checksum the normalized
  // config (and therefore the config hash and restart identity) records only
  // the path, so a replaced file at the same path would silently change the
  // physics across a restart. Every configured data path therefore requires
  // its 64-hex-digit lower-case SHA-256 (verified against the file at load).
  const auto require_checksum = [](const std::string& path, const std::string& sha, const char* key) {
    if (path.empty()) return;
    const bool valid = sha.size() == 64U && std::all_of(sha.begin(), sha.end(), [](char c) {
      return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
    });
    if (!valid) {
      throw ConfigError(std::string("physics.") + key + "_path is set, so physics." + key +
                        "_sha256 must hold the file's 64-digit lower-case SHA-256 (data files are bound to the run "
                        "by checksum; tools/chui_data/validate_treecool.py and the converters print it)");
    }
  };
  require_checksum(p.uv_background_table_path, p.uv_background_table_sha256, "uv_background_table");
  require_checksum(p.metal_cooling_table_path, p.metal_cooling_table_sha256, "metal_cooling_table");
  require_checksum(p.ssp_yield_table_path, p.ssp_yield_table_sha256, "ssp_yield_table");
  if (p.stellar_population_model != StellarPopulationModelConfig::kLegacyTable && p.enable_feedback &&
      p.fb_supernova_model != SupernovaFeedbackModelConfig::kMechanicalResolvedUnresolved) {
    throw ConfigError(
        "a Part 2 stellar_population_model deposits element/dust-resolved ejecta through "
        "fb_supernova_model=mechanical_resolved_unresolved (the legacy deposition cannot carry them)");
  }
  // CHUI Part 3: correlated associations (RESEARCH, EXPERIMENTAL/UNQUALIFIED).
  if (!std::isfinite(p.ssp_association_radius_code) || p.ssp_association_radius_code < 0.0 ||
      !std::isfinite(p.ssp_association_time_window_years) || !(p.ssp_association_time_window_years > 0.0)) {
    throw ConfigError("physics.ssp_association_radius_code must be >= 0 and ssp_association_time_window_years > 0");
  }
  if (p.ssp_clustered_associations &&
      (p.ssp_association_radius_code == 0.0 || p.stellar_population_model == StellarPopulationModelConfig::kLegacyTable)) {
    throw ConfigError(
        "physics.ssp_clustered_associations=true requires ssp_association_radius_code > 0 and a Part 2 "
        "stellar_population_model (associations group SSP feedback sources)");
  }
  if (p.ssp_wind_model == SspWindModelConfig::kTableV1) {
    require_checksum(p.ssp_wind_table_path, p.ssp_wind_table_sha256, "ssp_wind_table");
    if (p.ssp_wind_table_path.empty() || p.stellar_population_model == StellarPopulationModelConfig::kLegacyTable) {
      throw ConfigError("physics.ssp_wind_model=table_v1 requires ssp_wind_table_path (+ sha256) and a Part 2 "
                        "stellar_population_model");
    }
    if (p.fb_presn_wind_momentum_lsun_per_c > 0.0) {
      throw ConfigError("physics.ssp_wind_model=table_v1 and fb_presn_wind_momentum_lsun_per_c > 0 would count "
                        "pre-SN momentum twice; select one channel");
    }
  } else if (!p.ssp_wind_table_path.empty()) {
    throw ConfigError("physics.ssp_wind_table_path is set but ssp_wind_model = none (it would have no effect)");
  }
  if (p.fb_supernova_model == SupernovaFeedbackModelConfig::kMechanicalResolvedUnresolved &&
      p.stellar_population_model == StellarPopulationModelConfig::kLegacyTable) {
    throw ConfigError("fb_supernova_model=mechanical_resolved_unresolved consumes discrete SN events from a Part 2 stellar_population_model");
  }
}

void validateConfig(const SimulationConfig& config) {
  if (config.schema_version != 1) {
    throw ConfigError("schema_version must be 1 for this build");
  }

  requireAllFinite({
      {config.cosmology.omega_matter, "cosmology.omega_matter"},
      {config.cosmology.omega_lambda, "cosmology.omega_lambda"},
      {config.cosmology.omega_baryon, "cosmology.omega_baryon"},
      {config.cosmology.hubble_param, "cosmology.hubble_param"},
      {config.cosmology.sigma8, "cosmology.sigma8"},
      {config.cosmology.scalar_index_ns, "cosmology.scalar_index_ns"},
      {config.cosmology.box_size_x_mpc_comoving, "cosmology.box_size_x"},
      {config.cosmology.box_size_y_mpc_comoving, "cosmology.box_size_y"},
      {config.cosmology.box_size_z_mpc_comoving, "cosmology.box_size_z"},
      {config.cosmology.box_size_mpc_comoving, "cosmology.box_size"},
      {config.numerics.a_begin, "numerics.a_begin"},
      {config.numerics.a_end, "numerics.a_end"},
      {config.numerics.z_begin, "numerics.z_begin"},
      {config.numerics.z_end, "numerics.z_end"},
      {config.numerics.t_code_begin, "numerics.t_code_begin"},
      {config.numerics.t_code_end, "numerics.t_code_end"},
      {config.numerics.t_phys_begin, "numerics.t_phys_begin"},
      {config.numerics.t_phys_end, "numerics.t_phys_end"},
      {config.numerics.cosmology_max_delta_ln_a, "numerics.cosmology_max_delta_ln_a"},
      {config.numerics.cosmology_max_hubble_time_fraction, "numerics.cosmology_max_hubble_time_fraction"},
      {config.numerics.source_max_fractional_change, "numerics.source_max_fractional_change"},
      {config.numerics.hydro_density_floor_code, "numerics.hydro_density_floor_code"},
      {config.numerics.hydro_pressure_floor_code, "numerics.hydro_pressure_floor_code"},
      {config.numerics.gravity_softening_kpc_comoving, "numerics.gravity_softening"},
      {config.numerics.gravity_softening_gas_kpc_comoving, "numerics.gravity_softening_gas"},
      {config.numerics.gravity_softening_dark_matter_kpc_comoving, "numerics.gravity_softening_dark_matter"},
      {config.numerics.gravity_softening_star_kpc_comoving, "numerics.gravity_softening_star"},
      {config.numerics.gravity_softening_black_hole_kpc_comoving, "numerics.gravity_softening_black_hole"},
      {config.numerics.gravity_softening_tracer_kpc_comoving, "numerics.gravity_softening_tracer"},
      {config.numerics.treepm_asmth_cells, "numerics.treepm_asmth_cells"},
      {config.numerics.treepm_rcut_cells, "numerics.treepm_rcut_cells"},
      {config.numerics.treepm_tree_opening_theta, "numerics.treepm_tree_opening_theta"},
      {config.numerics.treepm_tree_relative_force_tolerance,
       "numerics.treepm_tree_relative_force_tolerance"},
      {config.numerics.treepm_tree_relative_force_acceleration_floor,
       "numerics.treepm_tree_relative_force_acceleration_floor"},
      {config.mode.zoom_region_center_x_mpc_comoving, "mode.zoom_region_center_x"},
      {config.mode.zoom_region_center_y_mpc_comoving, "mode.zoom_region_center_y"},
      {config.mode.zoom_region_center_z_mpc_comoving, "mode.zoom_region_center_z"},
      {config.mode.zoom_region_radius_mpc_comoving, "mode.zoom_region_radius"},
      {config.mode.zoom_contamination_radius_mpc_comoving, "mode.zoom_contamination_radius"},
      {config.physics.temperature_floor_k, "physics.temperature_floor_k"},
      {config.physics.sf_density_threshold_code, "physics.sf_density_threshold_code"},
      {config.physics.sf_temperature_threshold_k, "physics.sf_temperature_threshold_k"},
      {config.physics.sf_min_converging_flow_rate_code, "physics.sf_min_converging_flow_rate_code"},
      {config.physics.sf_epsilon_ff, "physics.sf_epsilon_ff"},
      {config.physics.sf_bound_alpha_vir_max, "physics.sf_bound_alpha_vir_max"},
      {config.physics.sf_jeans_mass_floor_code, "physics.sf_jeans_mass_floor_code"},
      {config.physics.sf_target_star_particle_mass_code, "physics.sf_target_star_particle_mass_code"},
      {config.physics.sf_target_star_particle_mass_fraction, "physics.sf_target_star_particle_mass_fraction"},
      {config.physics.sf_min_star_particle_mass_code, "physics.sf_min_star_particle_mass_code"},
      {config.physics.sf_max_star_particle_mass_code, "physics.sf_max_star_particle_mass_code"},
      {config.physics.sf_max_fractional_mass_conversion, "physics.sf_max_fractional_mass_conversion"},
      {config.physics.sf_min_remaining_gas_fraction, "physics.sf_min_remaining_gas_fraction"},
      {config.physics.sf_min_remaining_gas_mass_code, "physics.sf_min_remaining_gas_mass_code"},
      {config.physics.sf_temperature_safety_ceiling_k, "physics.sf_temperature_safety_ceiling_k"},
      {config.physics.sf_effective_n_h_threshold_cgs, "physics.sf_effective_n_h_threshold"},
      {config.physics.sf_effective_min_baryon_overdensity, "physics.sf_effective_min_baryon_overdensity"},
      {config.physics.sf_effective_t_star_at_threshold_code, "physics.sf_effective_t_star_at_threshold"},
      {config.physics.sf_effective_evaporation_factor_at_threshold, "physics.sf_effective_evaporation_factor_at_threshold"},
      {config.physics.sf_effective_evaporation_density_exponent, "physics.sf_effective_evaporation_density_exponent"},
      {config.physics.sf_effective_cold_phase_temperature_k, "physics.sf_effective_cold_phase_temperature"},
      {config.physics.sf_effective_supernova_specific_energy_code, "physics.sf_effective_supernova_specific_energy_code"},
      {config.physics.sf_effective_massive_star_fraction, "physics.sf_effective_massive_star_fraction"},
      {config.physics.sf_effective_q_eos, "physics.sf_effective_q_eos"},
      {config.physics.sf_effective_isothermal_temperature_k, "physics.sf_effective_isothermal_temperature"},
      {config.physics.sf_effective_hot_excess_tolerance, "physics.sf_effective_hot_excess_tolerance"},
      {config.physics.sf_effective_eos_relaxation_timescale_code, "physics.sf_effective_eos_relaxation_timescale"},
      {config.physics.sf_effective_eos_max_density_ratio, "physics.sf_effective_eos_max_density_ratio"},
      {config.physics.fb_epsilon_thermal, "physics.fb_epsilon_thermal"},
      {config.physics.fb_epsilon_kinetic, "physics.fb_epsilon_kinetic"},
      {config.physics.fb_epsilon_momentum, "physics.fb_epsilon_momentum"},
      {config.physics.fb_sn_energy_erg_per_mass_code, "physics.fb_sn_energy_erg_per_mass_code"},
      {config.physics.fb_momentum_code_per_mass_code, "physics.fb_momentum_code_per_mass_code"},
      {config.physics.fb_delayed_cooling_time_code, "physics.fb_delayed_cooling_time_code"},
      {config.physics.fb_stochastic_event_probability, "physics.fb_stochastic_event_probability"},
      {config.physics.stellar_evolution_hubble_time_years, "physics.stellar_evolution_hubble_time_years"},
      {config.physics.bh_seed_halo_mass_threshold_code, "physics.bh_seed_halo_mass_threshold_code"},
      {config.physics.bh_seed_mass_code, "physics.bh_seed_mass_code"},
      {config.physics.bh_alpha_bondi, "physics.bh_alpha_bondi"},
      {config.physics.bh_epsilon_r, "physics.bh_epsilon_r"},
      {config.physics.bh_epsilon_f, "physics.bh_epsilon_f"},
      {config.physics.bh_feedback_coupling_efficiency, "physics.bh_feedback_coupling_efficiency"},
      {config.physics.bh_duty_cycle_active_edd_ratio_threshold, "physics.bh_duty_cycle_active_edd_ratio_threshold"},
      {config.physics.bh_proton_mass_si, "physics.bh_proton_mass_si"},
      {config.physics.bh_thomson_cross_section_si, "physics.bh_thomson_cross_section_si"},
      {config.physics.bh_newton_g_si, "physics.bh_newton_g_si"},
      {config.physics.bh_speed_of_light_si, "physics.bh_speed_of_light_si"},
      {config.physics.tracer_min_host_mass_code, "physics.tracer_min_host_mass_code"},
      {config.physics.thermochemistry_relative_tolerance, "physics.thermochemistry_relative_tolerance"},
      {config.physics.thermochemistry_abundance_absolute_tolerance, "physics.thermochemistry_abundance_absolute_tolerance"},
      {config.physics.thermochemistry_analytic_relaxation_time_years, "physics.thermochemistry_analytic_relaxation_time_years"},
      {config.physics.thermochemistry_analytic_equilibrium_temperature_k, "physics.thermochemistry_analytic_equilibrium_temperature_k"},
      {config.physics.isrf_g0_habing, "physics.isrf_g0_habing"},
      {config.physics.chemistry_shielding_length_max_pc, "physics.chemistry_shielding_length_max_pc"},
      {config.physics.metal_cooling_switch_temperature_k, "physics.metal_cooling_switch_temperature_k"},
      {config.physics.chemistry_primordial_helium_mass_fraction, "physics.chemistry_primordial_helium_mass_fraction"},
      {config.physics.dust_accretion_timescale_ref_years, "physics.dust_accretion_timescale_ref_years"},
      {config.physics.dust_sputtering_timescale_ref_years, "physics.dust_sputtering_timescale_ref_years"},
      {config.physics.dust_sn_swept_mass_msun, "physics.dust_sn_swept_mass_msun"},
      {config.physics.dust_sn_destruction_efficiency, "physics.dust_sn_destruction_efficiency"},
      {config.physics.dust_agb_condensation_efficiency, "physics.dust_agb_condensation_efficiency"},
      {config.physics.dust_cc_condensation_efficiency, "physics.dust_cc_condensation_efficiency"},
      {config.physics.ism_jeans_resolution_cells, "physics.ism_jeans_resolution_cells"},
      {config.physics.ism_regime_hysteresis_fraction, "physics.ism_regime_hysteresis_fraction"},
      {config.physics.ism_regime_min_dwell_years, "physics.ism_regime_min_dwell_years"},
      {config.physics.sgs_dissipation_coefficient, "physics.sgs_dissipation_coefficient"},
      {config.physics.sgs_eddy_viscosity_coefficient, "physics.sgs_eddy_viscosity_coefficient"},
      {config.physics.sgs_feedback_energy_fraction, "physics.sgs_feedback_energy_fraction"},
      {config.physics.chemistry_adaptive_hysteresis_fraction, "physics.chemistry_adaptive_hysteresis_fraction"},
      {config.physics.chemistry_adaptive_min_dwell_years, "physics.chemistry_adaptive_min_dwell_years"},
      {config.physics.chemistry_adaptive_cold_temperature_k, "physics.chemistry_adaptive_cold_temperature_k"},
      {config.physics.chemistry_adaptive_cooling_time_ratio, "physics.chemistry_adaptive_cooling_time_ratio"},
      {config.physics.chemistry_adaptive_dense_nh_cm3, "physics.chemistry_adaptive_dense_nh_cm3"},
      {config.physics.chemistry_shadow_sample_fraction, "physics.chemistry_shadow_sample_fraction"},
      {config.physics.ssp_imf_mass_min_msun, "physics.ssp_imf_mass_min_msun"},
      {config.physics.ssp_imf_mass_max_msun, "physics.ssp_imf_mass_max_msun"},
      {config.physics.ssp_cc_mass_min_msun, "physics.ssp_cc_mass_min_msun"},
      {config.physics.ssp_cc_mass_max_msun, "physics.ssp_cc_mass_max_msun"},
      {config.physics.ssp_snia_per_msun, "physics.ssp_snia_per_msun"},
      {config.physics.ssp_snia_delay_min_years, "physics.ssp_snia_delay_min_years"},
      {config.physics.ssp_snia_delay_max_years, "physics.ssp_snia_delay_max_years"},
      {config.physics.ssp_energy_per_sn_erg, "physics.ssp_energy_per_sn_erg"},
      {config.physics.ssp_reference_net_metal_yield, "physics.ssp_reference_net_metal_yield"},
      {config.physics.fb_presn_wind_momentum_lsun_per_c, "physics.fb_presn_wind_momentum_lsun_per_c"},
      {config.analysis.halo_fof_linking_length_factor, "analysis.halo_fof_linking_length_factor"},
  });

  if (config.numerics.t_code_end <= config.numerics.t_code_begin) {
    throw ConfigError("numerics.t_code_end must be greater than numerics.t_code_begin");
  }
  if (config.numerics.t_phys_end < config.numerics.t_phys_begin) {
    throw ConfigError("numerics.t_phys_end must be greater than or equal to numerics.t_phys_begin");
  }
  validateScaleRedshiftPair(
      config.numerics.a_begin,
      config.numerics.z_begin,
      "numerics.a_begin",
      "numerics.z_begin");
  validateScaleRedshiftPair(
      config.numerics.a_end,
      config.numerics.z_end,
      "numerics.a_end",
      "numerics.z_end");
  if (config.numerics.a_end < config.numerics.a_begin) {
    throw ConfigError(
        "numerics.a_end must be >= numerics.a_begin for forward cosmological integration; "
        "reverse integration requires an explicit future config mode");
  }
  if (config.numerics.z_end > config.numerics.z_begin) {
    throw ConfigError(
        "numerics.z_end must be <= numerics.z_begin for forward cosmological integration; "
        "reverse integration requires an explicit future config mode");
  }
  if (config.numerics.cosmology_max_delta_ln_a <= 0.0 ||
      config.numerics.cosmology_max_hubble_time_fraction <= 0.0 ||
      config.numerics.source_max_fractional_change <= 0.0 ||
      config.numerics.source_max_fractional_change > 1.0) {
    throw ConfigError("numerics cosmology/source timestep limits must be positive; source_max_fractional_change must be <= 1");
  }
  if (config.numerics.hydro_density_floor_code < 0.0 ||
      config.numerics.hydro_pressure_floor_code < 0.0) {
    throw ConfigError("numerics hydro floors must be finite and non-negative code-unit values");
  }
  if (!std::isfinite(config.numerics.hydro_dual_energy_eta) ||
      config.numerics.hydro_dual_energy_eta <= 0.0 ||
      config.numerics.hydro_dual_energy_eta >= 1.0) {
    throw ConfigError("numerics.hydro_dual_energy_eta must be finite and in (0, 1)");
  }
  if (!std::isfinite(config.numerics.amr_predictive_feedback_age_years) ||
      config.numerics.amr_predictive_feedback_age_years < 0.0 ||
      (config.numerics.amr_predictive_feedback_age_years > 0.0 &&
       config.numerics.amr_regrid_mode != AmrRegridConfigMode::kAutomatic)) {
    throw ConfigError("numerics.amr_predictive_feedback_age_years must be >= 0 and requires amr_regrid_mode = automatic");
  }
  // R2 (A04): event-aware predictive refinement from the canonical SSP event
  // chronology; needs automatic regrid and a Part 2 population model (the
  // legacy table has no event records, which would make the key a no-op).
  if (!std::isfinite(config.numerics.amr_predictive_event_lookahead_years) ||
      config.numerics.amr_predictive_event_lookahead_years < 0.0 ||
      (config.numerics.amr_predictive_event_lookahead_years > 0.0 &&
       (config.numerics.amr_regrid_mode != AmrRegridConfigMode::kAutomatic ||
        config.physics.stellar_population_model == StellarPopulationModelConfig::kLegacyTable))) {
    throw ConfigError("numerics.amr_predictive_event_lookahead_years must be >= 0; > 0 requires "
                      "amr_regrid_mode = automatic and a Part 2 stellar_population_model");
  }
  if (config.numerics.global_step_max_retries < 0 || config.numerics.global_step_max_retries > 32) {
    throw ConfigError("numerics.global_step_max_retries must be in [0, 32]");
  }
  if (config.numerics.global_step_max_retries > 0 && config.numerics.hierarchical_max_rung > 0) {
    throw ConfigError("numerics.global_step_max_retries > 0 is implemented for the rung-zero (global KDK) driver "
                      "only; hierarchical runs use the multirate PM cycle transaction");
  }
  // CHUI Part 3 automatic AMR regrid.
  if (config.numerics.amr_regrid_interval_steps < 1 || config.numerics.amr_regrid_max_operations < 1) {
    throw ConfigError("numerics.amr_regrid_interval_steps and numerics.amr_regrid_max_operations must be >= 1");
  }
  if (!std::isfinite(config.numerics.amr_refine_cell_mass_code) || config.numerics.amr_refine_cell_mass_code < 0.0 ||
      !std::isfinite(config.numerics.amr_refine_jump_threshold) || config.numerics.amr_refine_jump_threshold < 0.0 ||
      !std::isfinite(config.numerics.amr_refine_jeans_cells) || config.numerics.amr_refine_jeans_cells < 0.0) {
    throw ConfigError("numerics.amr_refine_* criteria must be finite and >= 0 (0 disables a criterion)");
  }
  for (const double value : {config.numerics.amr_refine_cooling_length_cells, config.numerics.amr_refine_metallicity_jump,
                             config.numerics.amr_refine_species_jump, config.numerics.amr_refine_vorticity}) {
    if (!std::isfinite(value) || value < 0.0) {
      throw ConfigError("numerics.amr_refine_{cooling_length_cells,metallicity_jump,species_jump,vorticity} must be "
                        "finite and >= 0 (0 disables)");
    }
  }
  if (config.numerics.amr_refine_cooling_length_cells > 0.0 &&
      (!config.physics.enable_cooling ||
       config.physics.thermochemistry_backend != ThermochemistryBackendMode::kLegacyCooling)) {
    // Explicit refusal, not a silent no-op: the cooling time is evaluated from
    // the in-solver legacy cooling-rate provider. Part 2 backends own cooling
    // in the SourceRuntime and do not yet publish per-cell cooling times.
    throw ConfigError("numerics.amr_refine_cooling_length_cells > 0 requires physics.enable_cooling = true with "
                      "thermochemistry_backend = legacy_cooling (Part 2 backends do not publish cooling times)");
  }
  if (config.numerics.amr_refine_species_jump > 0.0 &&
      config.physics.thermochemistry_backend != ThermochemistryBackendMode::kFixedHybrid &&
      config.physics.thermochemistry_backend != ThermochemistryBackendMode::kFixedDetailed) {
    throw ConfigError("numerics.amr_refine_species_jump > 0 requires a composition network backend "
                      "(thermochemistry_backend = fixed_hybrid or fixed_detailed)");
  }
  if (!std::isfinite(config.numerics.amr_derefine_hysteresis) || config.numerics.amr_derefine_hysteresis <= 0.0 ||
      config.numerics.amr_derefine_hysteresis > 1.0) {
    throw ConfigError("numerics.amr_derefine_hysteresis must be in (0, 1]");
  }
  if (config.numerics.amr_regrid_mode == AmrRegridConfigMode::kAutomatic) {
    if (config.numerics.amr_refine_cell_mass_code == 0.0 && config.numerics.amr_refine_jump_threshold == 0.0 &&
        config.numerics.amr_refine_jeans_cells == 0.0 && config.numerics.amr_refine_cooling_length_cells == 0.0 &&
        config.numerics.amr_refine_metallicity_jump == 0.0 && config.numerics.amr_refine_species_jump == 0.0 &&
        config.numerics.amr_refine_vorticity == 0.0) {
      throw ConfigError("numerics.amr_regrid_mode = automatic requires at least one non-zero amr_refine_* criterion");
    }
    if (config.numerics.amr_max_level < 1 || config.numerics.amr_max_level > 30) {
      throw ConfigError("numerics.amr_regrid_mode = automatic requires numerics.amr_max_level in [1, 30]");
    }
  }
  if (config.parallel.node_memory_admission == NodeMemoryAdmissionMode::kOff &&
      (config.parallel.node_memory_budget_bytes != 0U ||
       config.parallel.node_memory_reserve_bytes != 0U)) {
    throw ConfigError(
        "parallel.node_memory_budget_bytes/node_memory_reserve_bytes require "
        "parallel.node_memory_admission = shared_equal; refusing to silently ignore a node budget");
  }
  if (config.parallel.node_memory_admission != NodeMemoryAdmissionMode::kOff &&
      config.parallel.node_memory_budget_bytes != 0U &&
      config.parallel.node_memory_reserve_bytes >= config.parallel.node_memory_budget_bytes) {
    throw ConfigError("parallel.node_memory_reserve_bytes must be smaller than parallel.node_memory_budget_bytes");
  }
  if (config.numerics.max_global_steps <= 0) {
    throw ConfigError("numerics.max_global_steps must be > 0");
  }
  if (config.numerics.hierarchical_max_rung < 0 || config.numerics.hierarchical_max_rung > 12) {
    throw ConfigError("numerics.hierarchical_max_rung must be in [0, 12]; 0 retains reference global KDK");
  }
  if (!std::isfinite(config.numerics.hierarchical_coarse_rms_displacement_factor) ||
      config.numerics.hierarchical_coarse_rms_displacement_factor <= 0.0 ||
      config.numerics.hierarchical_coarse_rms_displacement_factor > 1.0) {
    throw ConfigError("numerics.hierarchical_coarse_rms_displacement_factor must be finite and in (0, 1]");
  }
  // eta > 1 would let a particle cross more than its softening length per
  // step under constant acceleration (dt^2 |a| / eps = eta^2), which leaves the
  // regime the softened second-order KDK criterion is designed for.
  if (!std::isfinite(config.numerics.gravity_timestep_eta) ||
      config.numerics.gravity_timestep_eta <= 0.0 ||
      config.numerics.gravity_timestep_eta > 1.0) {
    throw ConfigError("numerics.gravity_timestep_eta must be finite and in (0, 1]");
  }
  if (config.numerics.hierarchical_max_rung > 0 &&
      (config.mode.mode != SimulationMode::kCosmoCube || config.physics.enable_cooling ||
       config.physics.enable_star_formation || config.physics.enable_feedback ||
       config.physics.enable_stellar_evolution || config.physics.enable_black_hole_agn ||
       config.physics.enable_tracers || config.physics.enable_metal_diffusion)) {
    throw ConfigError("hierarchical KDK currently requires cosmo_cube collisionless DMO with source physics disabled");
  }
  if (config.parallel.decomposition_spatial_work_enabled &&
      config.mode.mode != SimulationMode::kCosmoCube) {
    throw ConfigError("parallel.decomposition_spatial_work_enabled currently requires cosmo_cube fixed periodic domains");
  }
  if (config.output.snapshot_interval_steps < 0) {
    throw ConfigError("output.snapshot_interval_steps must be >= 0");
  }
  if (!std::isfinite(config.output.snapshot_interval_time_code) ||
      config.output.snapshot_interval_time_code < 0.0) {
    throw ConfigError("output.snapshot_interval_time_code must be finite and >= 0");
  }
  if (config.output.snapshot_interval_steps == 0 &&
      config.output.snapshot_interval_time_code == 0.0) {
    throw ConfigError(
        "output requires snapshot_interval_steps > 0 or snapshot_interval_time_code > 0");
  }
  if (config.output.snapshot_layout != ScienceSnapshotLayout::kAggregated &&
      config.output.snapshot_num_files > 1U) {
    throw ConfigError(
        "output.snapshot_num_files > 1 is only valid with snapshot_layout=aggregated");
  }
  if (config.output.restart_interval_steps < 0) {
    throw ConfigError("output.restart_interval_steps must be >= 0; 0 keeps checkpoints coupled to snapshot events");
  }
  if (config.output.restart_interval_steps > 0 && !config.output.write_restarts) {
    throw ConfigError("output.restart_interval_steps > 0 requires output.write_restarts = true");
  }
  if (!std::isfinite(config.output.restart_interval_time_code) ||
      config.output.restart_interval_time_code < 0.0) {
    throw ConfigError("output.restart_interval_time_code must be finite and >= 0");
  }
  if (config.output.restart_interval_time_code > 0.0 && !config.output.write_restarts) {
    throw ConfigError("output.restart_interval_time_code > 0 requires output.write_restarts = true");
  }
  if (config.parallel.mpi_ranks_expected <= 0) {
    throw ConfigError("parallel.mpi_ranks_expected must be positive");
  }
  if (config.parallel.omp_threads < 0) {
    throw ConfigError("parallel.omp_threads must be >= 0 (0 selects the runtime default)");
  }
#if !COSMOSIM_HAVE_OPENMP
  if (config.parallel.omp_threads > 1) {
    throw ConfigError(
        "parallel.omp_threads > 1 requires a binary built with OpenMP support");
  }
#endif
  if (config.parallel.gpu_devices < 0) {
    throw ConfigError("parallel.gpu_devices must be >= 0");
  }
  if (config.parallel.decomposition_particle_count_weight < 0.0 ||
      config.parallel.decomposition_gas_cell_weight < 0.0 ||
      config.parallel.decomposition_tree_interaction_weight < 0.0 ||
      config.parallel.decomposition_pm_mesh_weight < 0.0 ||
      config.parallel.decomposition_amr_patch_weight < 0.0 ||
      config.parallel.decomposition_active_fraction_weight < 0.0 ||
      config.parallel.decomposition_memory_pressure_weight < 0.0 ||
      config.parallel.decomposition_gpu_occupancy_weight < 0.0 ||
      config.parallel.decomposition_generic_work_weight < 0.0 ||
      config.parallel.decomposition_rebalance_imbalance_trigger < 1.0 ||
      config.parallel.decomposition_rebalance_memory_trigger < 1.0 ||
      config.parallel.decomposition_rebalance_max_migrated_load_fraction < 0.0 ||
      config.parallel.decomposition_rebalance_max_migrated_load_fraction > 1.0 ||
      config.parallel.decomposition_measured_tree_pair_weight < 0.0 ||
      config.parallel.decomposition_measured_pm_cell_weight < 0.0 ||
      config.parallel.decomposition_measured_amr_cell_weight < 0.0 ||
      config.parallel.decomposition_measured_hydro_face_weight < 0.0 ||
      config.parallel.decomposition_measured_wall_ms_weight < 0.0) {
    throw ConfigError("parallel decomposition work/rebalance weights must be non-negative; rebalance triggers must be >= 1");
  }
  if (config.parallel.isolated_pm_root_workspace_limit_bytes == 0U ||
      config.parallel.zoom_high_res_allgather_limit_bytes == 0U) {
    throw ConfigError("parallel PM gather guard byte limits must be > 0");
  }
  if (!std::isfinite(config.parallel.gravity_memory_safety_margin_fraction) ||
      config.parallel.gravity_memory_safety_margin_fraction < 0.0 ||
      config.parallel.gravity_memory_safety_margin_fraction > 1.0) {
    throw ConfigError(
        "parallel.gravity_memory_safety_margin_fraction must be finite and within [0,1]");
  }
  if (!std::isfinite(config.parallel.process_memory_safety_margin_fraction) ||
      config.parallel.process_memory_safety_margin_fraction < 0.0 ||
      config.parallel.process_memory_safety_margin_fraction > 1.0) {
    throw ConfigError(
        "parallel.process_memory_safety_margin_fraction must be finite and within [0,1]");
  }
  if (config.analysis.run_health_interval_steps <= 0 ||
      config.analysis.science_light_interval_steps <= 0 ||
      config.analysis.science_heavy_interval_steps <= 0) {
    throw ConfigError("analysis cadence intervals must be > 0");
  }
  if (config.analysis.retention_bundle_count < 1) {
    throw ConfigError("analysis.retention_bundle_count must be >= 1");
  }
  if (config.analysis.power_spectrum_mesh_n < 4 ||
      config.analysis.power_spectrum_bin_count < 1 ||
      config.analysis.sf_history_bin_count < 1 ||
      config.analysis.quicklook_grid_n < 4) {
    throw ConfigError("analysis mesh/bin settings must be within conservative minimum bounds");
  }
#if !COSMOSIM_ENABLE_FFTW
  const auto power_spectrum_mesh_n =
      static_cast<unsigned int>(config.analysis.power_spectrum_mesh_n);
  if ((power_spectrum_mesh_n & (power_spectrum_mesh_n - 1U)) != 0U) {
    throw ConfigError(
        "analysis.power_spectrum_mesh_n must be a power of two when this binary uses "
        "the built-in radix-2 FFT backend; rebuild with COSMOSIM_ENABLE_FFTW=ON "
        "to use FFTW-supported non-power-of-two mesh sizes");
  }
#endif
  if (config.analysis.halo_fof_linking_length_factor <= 0.0 ||
      config.analysis.halo_fof_linking_length_factor > 1.0) {
    throw ConfigError("analysis.halo_fof_linking_length_factor must be in (0, 1]");
  }
  if (config.analysis.halo_fof_min_group_size < 2) {
    throw ConfigError("analysis.halo_fof_min_group_size must be >= 2");
  }
  if (config.cosmology.omega_matter <= 0.0 || config.cosmology.omega_lambda < 0.0 ||
      config.cosmology.omega_baryon < 0.0 || config.cosmology.hubble_param <= 0.0 ||
      config.cosmology.sigma8 <= 0.0 || config.cosmology.scalar_index_ns <= 0.0) {
    throw ConfigError(
        "cosmology requires omega_matter > 0, omega_lambda >= 0, omega_baryon >= 0, "
        "hubble_param > 0, sigma8 > 0, and scalar_index_ns > 0");
  }
  if (config.cosmology.omega_baryon > config.cosmology.omega_matter) {
    throw ConfigError("cosmology.omega_baryon must be <= cosmology.omega_matter");
  }
  constexpr double k_cosmology_closure_tolerance = 1.0e-10;
  const double flat_density_sum =
      config.cosmology.omega_matter + config.cosmology.omega_lambda;
  if (!std::isfinite(flat_density_sum) ||
      std::abs(flat_density_sum - 1.0) > k_cosmology_closure_tolerance) {
    throw ConfigError(
        "cosmology currently supports a flat radiation-free matter+Lambda background; "
        "omega_matter + omega_lambda must equal 1 within 1e-10 so H(a=1) == H0");
  }
  auto requireOptionalPositiveSoftening = [](double value, const char* key) {
    if (value == 0.0) {
      throw ConfigError(std::string(key) + " must be > 0 when specified");
    }
  };
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_gas_kpc_comoving, "numerics.gravity_softening_gas");
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_dark_matter_kpc_comoving, "numerics.gravity_softening_dark_matter");
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_star_kpc_comoving, "numerics.gravity_softening_star");
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_black_hole_kpc_comoving, "numerics.gravity_softening_black_hole");
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_tracer_kpc_comoving, "numerics.gravity_softening_tracer");
  if (config.physics.enable_feedback && !config.physics.enable_star_formation) {
    throw ConfigError(
        "physics.enable_feedback requires physics.enable_star_formation=true or an explicit future external-source policy");
  }
  if (config.physics.enable_black_hole_agn && !config.physics.enable_feedback) {
    throw ConfigError("physics.enable_black_hole_agn requires physics.enable_feedback=true");
  }
  if (config.physics.bh_enable_seeding && !config.physics.enable_black_hole_agn) {
    throw ConfigError(
        "physics.bh_enable_seeding=true requires physics.enable_black_hole_agn=true");
  }
  if (config.physics.enable_cooling &&
      config.physics.cooling_model == CoolingModel::kPrimordialMetalLine &&
      trim(config.physics.metal_line_table_path).empty()) {
    throw ConfigError(
        "physics.cooling_model=primordial_metal_line requires physics.metal_line_table_path");
  }
  if (config.physics.temperature_floor_k <= 0.0) {
    throw ConfigError("physics.temperature_floor_k must be > 0");
  }
  if (config.physics.sf_density_threshold_code <= 0.0 ||
      config.physics.sf_temperature_threshold_k <= 0.0 ||
      config.physics.sf_min_star_particle_mass_code <= 0.0) {
    throw ConfigError("legacy star formation thresholds and min star mass must be > 0");
  }
  if (config.physics.sf_epsilon_ff < 0.0 || config.physics.sf_epsilon_ff > 1.0) {
    throw ConfigError("physics.sf_epsilon_ff must be in [0, 1]");
  }
  if (config.physics.sf_bound_alpha_vir_max <= 0.0 ||
      config.physics.sf_jeans_mass_floor_code < 0.0) {
    throw ConfigError(
        "physics.sf_bound_alpha_vir_max must be > 0 and sf_jeans_mass_floor_code must be >= 0");
  }
  if (config.physics.sf_target_star_particle_mass_code <= 0.0 ||
      config.physics.sf_target_star_particle_mass_fraction <= 0.0 ||
      config.physics.sf_target_star_particle_mass_fraction > 1.0 ||
      config.physics.sf_max_star_particle_mass_code < config.physics.sf_min_star_particle_mass_code) {
    throw ConfigError(
        "star-particle target mass must be > 0, target fraction in (0, 1], and max mass >= min mass");
  }
  if (config.physics.sf_max_spawn_particles_per_cell_step == 0U) {
    throw ConfigError("physics.sf_max_spawn_particles_per_cell_step must be > 0");
  }
  if (config.physics.sf_max_fractional_mass_conversion <= 0.0 ||
      config.physics.sf_max_fractional_mass_conversion >= 1.0 ||
      config.physics.sf_min_remaining_gas_fraction < 0.0 ||
      config.physics.sf_min_remaining_gas_fraction >= 1.0 ||
      config.physics.sf_min_remaining_gas_mass_code < 0.0 ||
      config.physics.sf_temperature_safety_ceiling_k < 0.0) {
    throw ConfigError(
        "star-formation conversion/remnant fractions must be in their documented ranges and mass/temperature floors non-negative");
  }
  if (config.physics.sf_effective_n_h_threshold_cgs <= 0.0 ||
      config.physics.sf_effective_min_baryon_overdensity < 0.0 ||
      config.physics.sf_effective_t_star_at_threshold_code <= 0.0 ||
      config.physics.sf_effective_evaporation_factor_at_threshold <= 0.0 ||
      config.physics.sf_effective_cold_phase_temperature_k <= 0.0 ||
      config.physics.sf_effective_supernova_specific_energy_code <= 0.0 ||
      config.physics.sf_effective_massive_star_fraction < 0.0 ||
      config.physics.sf_effective_massive_star_fraction >= 1.0 ||
      config.physics.sf_effective_q_eos < 0.0 ||
      config.physics.sf_effective_q_eos > 1.0 ||
      config.physics.sf_effective_isothermal_temperature_k <= 0.0 ||
      config.physics.sf_effective_hot_excess_tolerance < 0.0 ||
      config.physics.sf_effective_eos_relaxation_timescale_code < 0.0 ||
      config.physics.sf_effective_eos_table_bins < 8U ||
      config.physics.sf_effective_eos_max_density_ratio <= 1.0) {
    throw ConfigError("effective multiphase star-formation parameters are outside their supported ranges");
  }
  if (config.physics.star_formation_model == StarFormationModelKind::kEffectiveMultiphaseTngLike &&
      config.physics.sf_effective_feedback_coupling == EffectiveIsmFeedbackCoupling::kEffectiveEosOnly &&
      config.physics.enable_feedback) {
    throw ConfigError(
        "physics.enable_feedback must be false when sf_effective_feedback_coupling=effective_eos_only");
  }
  if (config.physics.fb_epsilon_thermal < 0.0 || config.physics.fb_epsilon_kinetic < 0.0 ||
      config.physics.fb_epsilon_momentum < 0.0) {
    throw ConfigError("physics feedback efficiencies must be >= 0");
  }
  if (config.physics.fb_sn_energy_erg_per_mass_code <= 0.0 ||
      config.physics.fb_momentum_code_per_mass_code < 0.0) {
    throw ConfigError("physics feedback budget scales must be physically non-negative");
  }
  if (config.physics.fb_neighbor_count == 0) {
    throw ConfigError("physics.fb_neighbor_count must be > 0");
  }
  if (config.physics.fb_stochastic_event_probability <= 0.0 ||
      config.physics.fb_stochastic_event_probability > 1.0) {
    throw ConfigError("physics.fb_stochastic_event_probability must be in (0, 1]");
  }
  if (config.physics.stellar_evolution_hubble_time_years <= 0.0) {
    throw ConfigError("physics.stellar_evolution_hubble_time_years must be > 0");
  }
  if (config.physics.stellar_evolution_require_production_table &&
      config.physics.stellar_evolution_table_path.empty()) {
    throw ConfigError(
        "physics.stellar_evolution_table_path is required when "
        "stellar_evolution_require_production_table=true");
  }
  // CHUI Part 2: core_elements is the element-resolved composition sidecar,
  // available only together with a Part 2 thermochemistry backend (which owns
  // the element/species manifest). The pairing is checked in
  // validatePart2BaryonicPhysics().
  if (config.physics.enable_metal_diffusion) {
    if (config.physics.metal_diffusion_model == MetalDiffusionModel::kNone) {
      throw ConfigError("metal diffusion is enabled but metal_diffusion_model=none");
    }
    if (!(config.physics.metal_diffusion_coefficient > 0.0) ||
        !(config.physics.metal_diffusion_cfl > 0.0) ||
        config.physics.metal_diffusion_cfl > 0.5) {
      throw ConfigError(
          "metal diffusion requires coefficient > 0 and cfl in (0, 0.5]");
    }
    if (config.physics.metal_diffusion_max_subcycles == 0U ||
        config.physics.metal_diffusion_max_rkl_stages < 2U) {
      throw ConfigError("metal diffusion stage/subcycle limits are invalid");
    }
    if (config.physics.metal_diffusion_coefficient_floor_code < 0.0 ||
        config.physics.metal_diffusion_coefficient_ceiling_code <=
            config.physics.metal_diffusion_coefficient_floor_code) {
      throw ConfigError("metal diffusion coefficient floor/ceiling are invalid");
    }
    if (config.numerics.hierarchical_max_rung != 0) {
      throw ConfigError(
          "metal diffusion currently requires numerics.hierarchical_max_rung=0");
    }
  }
  if (config.physics.bh_seed_halo_mass_threshold_code <= 0.0 || config.physics.bh_seed_mass_code <= 0.0 ||
      config.physics.bh_seed_max_per_cell == 0) {
    throw ConfigError("physics BH seeding parameters must be > 0");
  }
  if (config.physics.bh_alpha_bondi <= 0.0) {
    throw ConfigError("physics.bh_alpha_bondi must be > 0");
  }
  if (config.physics.bh_epsilon_r <= 0.0 || config.physics.bh_epsilon_r > 1.0 ||
      config.physics.bh_epsilon_f < 0.0 || config.physics.bh_epsilon_f > 1.0 ||
      config.physics.bh_feedback_coupling_efficiency < 0.0 ||
      config.physics.bh_feedback_coupling_efficiency > 1.0) {
    throw ConfigError("physics BH efficiencies must be within conservative physical bounds");
  }
  if (config.physics.bh_duty_cycle_active_edd_ratio_threshold < 0.0 ||
      config.physics.bh_duty_cycle_active_edd_ratio_threshold > 1.0) {
    throw ConfigError("physics.bh_duty_cycle_active_edd_ratio_threshold must be in [0, 1]");
  }
  if (config.physics.bh_proton_mass_si <= 0.0 || config.physics.bh_thomson_cross_section_si <= 0.0 ||
      config.physics.bh_newton_g_si <= 0.0 || config.physics.bh_speed_of_light_si <= 0.0) {
    throw ConfigError("physics BH constants must be > 0");
  }
  if (config.physics.tracer_min_host_mass_code < 0.0) {
    throw ConfigError("physics.tracer_min_host_mass_code must be >= 0");
  }
  validatePart2BaryonicPhysics(config);
  if (config.cosmology.box_size_x_mpc_comoving <= 0.0 ||
      config.cosmology.box_size_y_mpc_comoving <= 0.0 ||
      config.cosmology.box_size_z_mpc_comoving <= 0.0) {
    throw ConfigError("cosmology.box_size_{x,y,z} must all be > 0");
  }
  if (config.numerics.treepm_pm_grid_nx <= 0 ||
      config.numerics.treepm_pm_grid_ny <= 0 ||
      config.numerics.treepm_pm_grid_nz <= 0) {
    throw ConfigError("numerics.treepm_pm_grid_n{xyz} must all be > 0");
  }
  if (!std::isfinite(config.numerics.treepm_asmth_cells) ||
      config.numerics.treepm_asmth_cells <= 0.0) {
    throw ConfigError("numerics.treepm_asmth_cells must be finite and > 0");
  }
  if (!std::isfinite(config.numerics.treepm_rcut_cells) ||
      config.numerics.treepm_rcut_cells <= 0.0) {
    throw ConfigError("numerics.treepm_rcut_cells must be finite and > 0");
  }
  if (!std::isfinite(config.numerics.treepm_tree_opening_theta) ||
      config.numerics.treepm_tree_opening_theta <= 0.0) {
    throw ConfigError("numerics.treepm_tree_opening_theta must be finite and > 0");
  }
  const int block = config.numerics.treepm_residual_block_size;
  if (block != 16 && block != 32 && block != 64 && block != 128 && block != 256)
    throw ConfigError("numerics.treepm_residual_block_size must be 16, 32, 64, 128 or 256");
  if (config.numerics.treepm_tree_max_leaf_size < 1 || config.numerics.treepm_tree_max_leaf_size > 1024)
    throw ConfigError("numerics.treepm_tree_max_leaf_size must be in [1,1024]");
  if (!std::isfinite(config.numerics.treepm_adaptive_maximum_opening_angle) ||
      config.numerics.treepm_adaptive_maximum_opening_angle <= 0.0 ||
      config.numerics.treepm_adaptive_maximum_opening_angle > 1.0) {
    throw ConfigError("numerics.treepm_adaptive_maximum_opening_angle must be finite and in (0, 1]");
  }
  if (config.numerics.treepm_group_walk_enabled &&
      (!config.numerics.treepm_adaptive_acceptance_enabled ||
       config.numerics.treepm_adaptive_error_estimator != TreePmAdaptiveErrorEstimator::kRelativeMonopole ||
       !config.numerics.treepm_gaussian_pair_lookup_enabled ||
       config.numerics.treepm_leaf_simd_enabled || config.numerics.treepm_full_residual_lookup_enabled)) {
    throw ConfigError(
        "numerics.treepm_group_walk_enabled requires treepm_adaptive_acceptance_enabled=true, "
        "treepm_adaptive_error_estimator=relative_monopole, treepm_gaussian_pair_lookup_enabled=true "
        "and leaf SIMD/full-residual lookup disabled");
  }
  // Groups never cross a residual block, so sizes above the block size are
  // meaningless; 0 keeps the legacy one-source-leaf grouping.
  if (config.numerics.treepm_group_walk_target_size < 0 ||
      config.numerics.treepm_group_walk_target_size > config.numerics.treepm_residual_block_size ||
      config.numerics.treepm_group_walk_target_size == 1) {
    throw ConfigError(
        "numerics.treepm_group_walk_target_size must be 0 (one source leaf) or in [2, treepm_residual_block_size]");
  }
  if (config.numerics.treepm_group_walk_target_size > 0 && !config.numerics.treepm_group_walk_enabled) {
    throw ConfigError("numerics.treepm_group_walk_target_size > 0 requires treepm_group_walk_enabled=true");
  }
  if (!std::isfinite(config.numerics.treepm_tree_relative_force_tolerance) ||
      config.numerics.treepm_tree_relative_force_tolerance <= 0.0) {
    throw ConfigError(
        "numerics.treepm_tree_relative_force_tolerance must be finite and > 0");
  }
  if (!std::isfinite(config.numerics.treepm_tree_relative_force_acceleration_floor) ||
      config.numerics.treepm_tree_relative_force_acceleration_floor <= 0.0) {
    throw ConfigError(
        "numerics.treepm_tree_relative_force_acceleration_floor must be finite and > 0");
  }
  // Exact reuse never admits a field built from another source state, so the
  // reference policy keeps cadence 1. The experimental stale-window policy
  // validates decomposition epoch, mesh, split and field version and reuses the
  // field for at most cadence-1 later global KDK force evaluations.
  if (config.numerics.treepm_pm_reuse_policy == TreePmPmReusePolicy::kExactSourceState &&
      config.numerics.treepm_update_cadence_steps != 1) {
    throw ConfigError(
        "numerics.treepm_update_cadence_steps must be exactly 1 with treepm_pm_reuse_policy=exact; "
        "select treepm_pm_reuse_policy=stale_within_cadence for experimental PM reuse");
  }
  if (config.numerics.treepm_pm_reuse_policy == TreePmPmReusePolicy::kStaleWithinCadence &&
      (config.numerics.treepm_update_cadence_steps < 2 || config.numerics.treepm_update_cadence_steps > 64 ||
       config.numerics.hierarchical_max_rung != 0)) {
    throw ConfigError(
        "numerics.treepm_pm_reuse_policy=stale_within_cadence requires treepm_update_cadence_steps in [2,64] "
        "and global KDK (hierarchical_max_rung=0; hierarchical KDK already applies PM at coarse endpoints)");
  }
  // Multirate outer PM level: a cycle of N coarse hierarchical blocks with PM
  // at the cycle endpoints only. Requires the hierarchical integrator.
  if (config.numerics.treepm_pm_reuse_policy == TreePmPmReusePolicy::kMultirateOuterKick &&
      (config.numerics.treepm_update_cadence_steps < 2 || config.numerics.treepm_update_cadence_steps > 64 ||
       config.numerics.hierarchical_max_rung <= 0)) {
    throw ConfigError(
        "numerics.treepm_pm_reuse_policy=multirate_outer_kick requires treepm_update_cadence_steps in [2,64] "
        "(coarse blocks per PM cycle) and hierarchical KDK (hierarchical_max_rung > 0)");
  }
  if (config.numerics.treepm_pm_force_gradient == TreePmForceGradient::kFiniteDifference4) {
    const bool periodic_gravity =
        config.mode.gravity_boundary == ModeGravityBoundary::kPeriodic ||
        (config.mode.gravity_boundary == ModeGravityBoundary::kAuto &&
         (config.mode.mode == SimulationMode::kCosmoCube || config.mode.mode == SimulationMode::kZoomIn));
    if (!periodic_gravity || config.parallel.gpu_devices > 0 ||
        config.numerics.treepm_pm_grid_nx < 4 || config.numerics.treepm_pm_grid_ny < 4 ||
        config.numerics.treepm_pm_grid_nz < 4) {
      throw ConfigError(
          "numerics.treepm_pm_force_gradient=finite_difference_4 requires periodic host PM with "
          ">= 4 cells per axis");
    }
  }
  if (config.numerics.treepm_tree_exchange_batch_bytes == 0) {
    throw ConfigError("numerics.treepm_tree_exchange_batch_bytes must be > 0");
  }
  if (config.numerics.treepm_pm_exchange_batch_bytes == 0) {
    throw ConfigError("numerics.treepm_pm_exchange_batch_bytes must be > 0");
  }
  if (config.mode.ic_chunk_particle_count == 0U) {
    throw ConfigError("mode.ic_chunk_particle_count must be > 0");
  }
  if (config.mode.ic_staging_particle_count == 0U) {
    throw ConfigError("mode.ic_staging_particle_count must be > 0");
  }
  if (config.mode.ic_convention == InitialConditionConvention::kGenerated) {
    if (config.mode.ic_file != "generated") {
      throw ConfigError(
          "mode.ic_convention=generated requires mode.ic_file=generated; external HDF5 input must select an explicit versioned convention");
    }
    if (!config.mode.ic_manifest_file.empty()) {
      throw ConfigError(
          "mode.ic_manifest_file must be empty when mode.ic_convention=generated");
    }
  } else {
    const bool manifest_mode =
        config.mode.ic_convention == InitialConditionConvention::kManifestV1;
    if (!manifest_mode &&
        (config.mode.ic_file.empty() || config.mode.ic_file == "generated")) {
      throw ConfigError(
          "external initial-condition conventions require a non-generated mode.ic_file path");
    }
    if (manifest_mode && config.mode.ic_manifest_file.empty()) {
      throw ConfigError(
          "mode.ic_convention=manifest_v1 requires mode.ic_manifest_file");
    }
    if (config.mode.ic_convention ==
        InitialConditionConvention::kGadgetArepoBridgeV1) {
      std::vector<std::string> missing;
      const auto require_positive = [&](double value, std::string_view key) {
        if (!std::isfinite(value) || !(value > 0.0)) {
          missing.emplace_back(key);
        }
      };
      const auto require_finite = [&](double value, std::string_view key) {
        if (!std::isfinite(value)) {
          missing.emplace_back(key);
        }
      };
      require_positive(
          config.mode.ic_bridge_source_length_unit_to_si,
          "mode.ic_bridge_source_length_unit_to_si");
      require_positive(
          config.mode.ic_bridge_source_mass_unit_to_si,
          "mode.ic_bridge_source_mass_unit_to_si");
      require_positive(
          config.mode.ic_bridge_source_velocity_unit_to_si,
          "mode.ic_bridge_source_velocity_unit_to_si");
      if (config.mode.ic_bridge_coordinate_frame ==
          InitialConditionCoordinateFrame::kUnspecified) {
        missing.emplace_back("mode.ic_bridge_coordinate_frame");
      }
      if (config.mode.ic_bridge_velocity_convention ==
          InitialConditionVelocityConvention::kUnspecified) {
        missing.emplace_back("mode.ic_bridge_velocity_convention");
      }
      require_finite(
          config.mode.ic_bridge_length_hubble_exponent,
          "mode.ic_bridge_length_hubble_exponent");
      require_finite(
          config.mode.ic_bridge_length_scale_factor_exponent,
          "mode.ic_bridge_length_scale_factor_exponent");
      require_finite(
          config.mode.ic_bridge_mass_hubble_exponent,
          "mode.ic_bridge_mass_hubble_exponent");
      require_finite(
          config.mode.ic_bridge_mass_scale_factor_exponent,
          "mode.ic_bridge_mass_scale_factor_exponent");
      require_finite(
          config.mode.ic_bridge_velocity_hubble_exponent,
          "mode.ic_bridge_velocity_hubble_exponent");
      require_finite(
          config.mode.ic_bridge_velocity_scale_factor_exponent,
          "mode.ic_bridge_velocity_scale_factor_exponent");
      if (!missing.empty()) {
        std::ostringstream message;
        message
            << "mode.ic_convention=gadget_arepo_bridge_v1 requires a complete "
               "explicit scientific convention; missing or invalid: ";
        for (std::size_t i = 0; i < missing.size(); ++i) {
          message << (i == 0 ? "" : ", ") << missing[i];
        }
        throw ConfigError(message.str());
      }
    }
  }

  const auto validate_missing_policy = [&config](
                                           InitialConditionMissingFieldPolicy policy,
                                           double configured_value,
                                           std::string_view policy_key,
                                           std::string_view value_key,
                                           bool require_positive) {
    if (policy == InitialConditionMissingFieldPolicy::kReconstruct) {
      throw ConfigError(
          std::string(policy_key) +
          "=reconstruct is not implemented by gadget_arepo_bridge_v1; "
          "select reject or an explicit supported value policy");
    }
    if (policy == InitialConditionMissingFieldPolicy::kPreserveUnavailable) {
      throw ConfigError(
          std::string(policy_key) +
          "=preserve_unavailable is not supported because the current runtime "
          "sidecars do not carry availability masks");
    }
    if (policy == InitialConditionMissingFieldPolicy::kUseConfigValue) {
      const bool valid = std::isfinite(configured_value) &&
          (require_positive ? configured_value > 0.0 : configured_value >= 0.0);
      if (!valid) {
        throw ConfigError(
            std::string(value_key) +
            (require_positive
                 ? " must be finite and > 0 when the matching policy is use_config_value"
                 : " must be finite and >= 0 when the matching policy is use_config_value"));
      }
    }
  };
  validate_missing_policy(
      config.mode.ic_gas_internal_energy_policy,
      config.mode.ic_gas_internal_energy_value_code,
      "mode.ic_gas_internal_energy_policy",
      "mode.ic_gas_internal_energy_value_code", true);
  validate_missing_policy(
      config.mode.ic_gas_density_policy,
      config.mode.ic_gas_density_value_code,
      "mode.ic_gas_density_policy", "mode.ic_gas_density_value_code", true);
  validate_missing_policy(
      config.mode.ic_star_formation_time_policy,
      config.mode.ic_star_formation_time_value,
      "mode.ic_star_formation_time_policy",
      "mode.ic_star_formation_time_value", true);
  validate_missing_policy(
      config.mode.ic_star_initial_mass_policy,
      config.mode.ic_star_initial_mass_value_code,
      "mode.ic_star_initial_mass_policy",
      "mode.ic_star_initial_mass_value_code", true);
  validate_missing_policy(
      config.mode.ic_star_metallicity_policy,
      config.mode.ic_star_metallicity_value,
      "mode.ic_star_metallicity_policy",
      "mode.ic_star_metallicity_value", false);
  validate_missing_policy(
      config.mode.ic_bh_mdot_policy, config.mode.ic_bh_mdot_value_code,
      "mode.ic_bh_mdot_policy", "mode.ic_bh_mdot_value_code", false);
  if (config.mode.zoom_region_radius_mpc_comoving < 0.0) {
    throw ConfigError("mode.zoom_region_radius must be >= 0");
  }
  if (config.mode.zoom_contamination_radius_mpc_comoving < 0.0) {
    throw ConfigError("mode.zoom_contamination_radius must be >= 0");
  }
  if (config.mode.zoom_high_res_region &&
      config.mode.zoom_region_radius_mpc_comoving <= 0.0) {
    throw ConfigError("mode.zoom_region_radius must be > 0 when mode.zoom_high_res_region is true");
  }
  if (config.mode.zoom_long_range_strategy ==
      ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection) {
    if (!config.mode.zoom_high_res_region) {
      throw ConfigError(
          "mode.zoom_long_range_strategy requires mode.zoom_high_res_region=true");
    }
    if (config.mode.zoom_focused_pm_grid_nx <= 0 ||
        config.mode.zoom_focused_pm_grid_ny <= 0 ||
        config.mode.zoom_focused_pm_grid_nz <= 0) {
      throw ConfigError(
          "mode.zoom_focused_pm_grid_n{xyz} must be > 0 for focused zoom PM correction");
    }
  }
  const ModePolicy policy = buildModePolicy(config.mode);
  validateModePolicy(config, policy);
  if (policy.gravity_boundary == GravityBoundaryModel::kPeriodicPoisson) {
    const double mesh_spacing_x = config.cosmology.box_size_x_mpc_comoving /
        static_cast<double>(config.numerics.treepm_pm_grid_nx);
    const double mesh_spacing_y = config.cosmology.box_size_y_mpc_comoving /
        static_cast<double>(config.numerics.treepm_pm_grid_ny);
    const double mesh_spacing_z = config.cosmology.box_size_z_mpc_comoving /
        static_cast<double>(config.numerics.treepm_pm_grid_nz);
    const double representative_mesh_spacing =
        std::cbrt(mesh_spacing_x * mesh_spacing_y * mesh_spacing_z);
    const double coarsest_mesh_spacing = std::max({
        mesh_spacing_x, mesh_spacing_y, mesh_spacing_z});
    const double split_scale =
        config.numerics.treepm_asmth_cells * representative_mesh_spacing;
    const double split_scale_boundary_tolerance =
        8.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(split_scale), std::abs(coarsest_mesh_spacing)});
    if (split_scale + split_scale_boundary_tolerance < coarsest_mesh_spacing) {
      throw ConfigError(
          "periodic TreePM rectangular-grid support requires the spherical split scale "
          "treepm_asmth_cells*cbrt(dx*dy*dz) to be >= the coarsest PM cell spacing; "
          "split_scale=" + std::to_string(split_scale) +
          ", coarsest_spacing=" + std::to_string(coarsest_mesh_spacing) +
          "; use a less anisotropic PM grid/box or increase treepm_asmth_cells");
    }

    // The current Plummer-softened short-range residual has been analytically
    // composition-tested through epsilon/r_s=0.20.  Plummer softening has a
    // non-compact tail, so silently accepting larger ratios under a fixed
    // finite cutoff would make the omitted residual an uncontrolled model
    // error.  Per-particle overrides are checked again at gravity-source
    // materialization because they are not part of static configuration.
    constexpr double k_max_certified_softening_to_split_ratio = 0.20;
    const double max_configured_softening_kpc = std::max({
        config.numerics.gravity_softening_kpc_comoving,
        config.numerics.gravity_softening_gas_kpc_comoving,
        config.numerics.gravity_softening_dark_matter_kpc_comoving,
        config.numerics.gravity_softening_star_kpc_comoving,
        config.numerics.gravity_softening_black_hole_kpc_comoving,
        config.numerics.gravity_softening_tracer_kpc_comoving,
    });
    const double max_configured_softening_mpc =
        max_configured_softening_kpc * 1.0e-3;
    if (max_configured_softening_mpc >
        k_max_certified_softening_to_split_ratio * split_scale) {
      throw ConfigError(
          "periodic TreePM Plummer softening exceeds the certified epsilon/r_s <= 0.20 "
          "envelope; max_configured_softening_mpc=" +
          std::to_string(max_configured_softening_mpc) +
          ", split_scale_mpc=" + std::to_string(split_scale) +
          "; reduce gravitational softening or increase treepm_asmth_cells/PM cell size");
    }
    const double cutoff_radius =
        config.numerics.treepm_rcut_cells * representative_mesh_spacing;
    const double shortest_box_axis = std::min({
        config.cosmology.box_size_x_mpc_comoving,
        config.cosmology.box_size_y_mpc_comoving,
        config.cosmology.box_size_z_mpc_comoving,
    });
    const double half_shortest_axis = 0.5 * shortest_box_axis;
    if (cutoff_radius >= std::nextafter(half_shortest_axis, 0.0)) {
      throw ConfigError(
          "periodic TreePM requires treepm_rcut_cells times the representative mesh spacing "
          "to be < half the shortest box axis; cutoff_radius=" +
          std::to_string(cutoff_radius) + ", half_shortest_axis=" +
          std::to_string(half_shortest_axis) + ", pm_grid=" +
          std::to_string(config.numerics.treepm_pm_grid_nx) + "x" +
          std::to_string(config.numerics.treepm_pm_grid_ny) + "x" +
          std::to_string(config.numerics.treepm_pm_grid_nz) +
          "; increase the PM grid or reduce treepm_rcut_cells");
    }
  }
  if (policy.gravity_boundary ==
          GravityBoundaryModel::kIsolatedMonopoleDirichlet &&
      config.numerics.treepm_enable_window_deconvolution) {
    throw ConfigError(
        "numerics.treepm_enable_window_deconvolution must be false for "
        "isolated/open gravity boundaries because no periodic assignment "
        "window is defined for that solver path");
  }
}


[[nodiscard]] std::string roundTripDoubleToString(double value) {
  char buffer[64]{};
  const auto result = std::to_chars(
      std::begin(buffer), std::end(buffer), value, std::chars_format::general);
  if (result.ec != std::errc{}) {
    throw ConfigError("failed to serialize normalized floating-point value");
  }
  return std::string(buffer, result.ptr);
}

class NormalizedConfigStream {
 public:
  NormalizedConfigStream& operator<<(double value) {
    m_stream << roundTripDoubleToString(value);
    return *this;
  }

  template <typename T>
  NormalizedConfigStream& operator<<(const T& value) {
    m_stream << value;
    return *this;
  }

  [[nodiscard]] std::string str() const { return m_stream.str(); }

 private:
  std::ostringstream m_stream;
};

[[nodiscard]] std::string buildNormalizedText(const FrozenConfig& frozen) {
  NormalizedConfigStream stream;
  const auto bridge_number = [](double value) {
    if (!std::isfinite(value)) {
      return std::string("unspecified");
    }
    return roundTripDoubleToString(value);
  };
  stream << "schema_version = " << frozen.config.schema_version << '\n';
  stream << "\n[units]\n";
  stream << "length_unit = " << frozen.config.units.length_unit << '\n';
  stream << "mass_unit = " << frozen.config.units.mass_unit << '\n';
  stream << "velocity_unit = " << frozen.config.units.velocity_unit << '\n';
  stream << "coordinate_frame = " << coordinateFrameToString(frozen.config.units.coordinate_frame) << '\n';
  stream << "\n[mode]\n";
  stream << "mode = " << modeToLowerString(frozen.config.mode.mode) << '\n';
  stream << "ic_file = " << frozen.config.mode.ic_file << '\n';
  stream << "ic_convention = "
         << initialConditionConventionToString(frozen.config.mode.ic_convention)
         << '\n';
  stream << "ic_manifest_file = " << frozen.config.mode.ic_manifest_file << '\n';
  stream << "ic_bridge_source_length_unit_to_si = "
         << frozen.config.mode.ic_bridge_source_length_unit_to_si << '\n';
  stream << "ic_bridge_source_mass_unit_to_si = "
         << frozen.config.mode.ic_bridge_source_mass_unit_to_si << '\n';
  stream << "ic_bridge_source_velocity_unit_to_si = "
         << frozen.config.mode.ic_bridge_source_velocity_unit_to_si << '\n';
  stream << "ic_bridge_coordinate_frame = "
         << initialConditionCoordinateFrameToString(
                frozen.config.mode.ic_bridge_coordinate_frame)
         << '\n';
  stream << "ic_bridge_velocity_convention = "
         << initialConditionVelocityConventionToString(
                frozen.config.mode.ic_bridge_velocity_convention)
         << '\n';
  stream << "ic_bridge_length_hubble_exponent = "
         << bridge_number(
                frozen.config.mode.ic_bridge_length_hubble_exponent)
         << '\n';
  stream << "ic_bridge_length_scale_factor_exponent = "
         << bridge_number(
                frozen.config.mode.ic_bridge_length_scale_factor_exponent)
         << '\n';
  stream << "ic_bridge_mass_hubble_exponent = "
         << bridge_number(
                frozen.config.mode.ic_bridge_mass_hubble_exponent)
         << '\n';
  stream << "ic_bridge_mass_scale_factor_exponent = "
         << bridge_number(
                frozen.config.mode.ic_bridge_mass_scale_factor_exponent)
         << '\n';
  stream << "ic_bridge_velocity_hubble_exponent = "
         << bridge_number(
                frozen.config.mode.ic_bridge_velocity_hubble_exponent)
         << '\n';
  stream << "ic_bridge_velocity_scale_factor_exponent = "
         << bridge_number(
                frozen.config.mode.ic_bridge_velocity_scale_factor_exponent)
         << '\n';
  stream << "ic_chunk_particle_count = "
         << frozen.config.mode.ic_chunk_particle_count << '\n';
  stream << "ic_staging_particle_count = "
         << frozen.config.mode.ic_staging_particle_count << '\n';
  stream << "ic_part_type2_policy = "
         << initialConditionSpeciesPolicyToString(
                frozen.config.mode.ic_part_type2_policy)
         << '\n';
  stream << "ic_part_type3_policy = "
         << initialConditionSpeciesPolicyToString(
                frozen.config.mode.ic_part_type3_policy)
         << '\n';
  stream << "ic_gas_internal_energy_policy = "
         << initialConditionMissingFieldPolicyToString(
                frozen.config.mode.ic_gas_internal_energy_policy)
         << '\n';
  stream << "ic_gas_density_policy = "
         << initialConditionMissingFieldPolicyToString(
                frozen.config.mode.ic_gas_density_policy)
         << '\n';
  stream << "ic_star_formation_time_policy = "
         << initialConditionMissingFieldPolicyToString(
                frozen.config.mode.ic_star_formation_time_policy)
         << '\n';
  stream << "ic_star_initial_mass_policy = "
         << initialConditionMissingFieldPolicyToString(
                frozen.config.mode.ic_star_initial_mass_policy)
         << '\n';
  stream << "ic_star_metallicity_policy = "
         << initialConditionMissingFieldPolicyToString(
                frozen.config.mode.ic_star_metallicity_policy)
         << '\n';
  stream << "ic_bh_mdot_policy = "
         << initialConditionMissingFieldPolicyToString(
                frozen.config.mode.ic_bh_mdot_policy)
         << '\n';
  stream << "ic_gas_internal_energy_value_code = "
         << frozen.config.mode.ic_gas_internal_energy_value_code << '\n';
  stream << "ic_gas_density_value_code = "
         << frozen.config.mode.ic_gas_density_value_code << '\n';
  stream << "ic_star_formation_time_value = "
         << frozen.config.mode.ic_star_formation_time_value << '\n';
  stream << "ic_star_initial_mass_value_code = "
         << frozen.config.mode.ic_star_initial_mass_value_code << '\n';
  stream << "ic_star_metallicity_value = "
         << frozen.config.mode.ic_star_metallicity_value << '\n';
  stream << "ic_bh_mdot_value_code = "
         << frozen.config.mode.ic_bh_mdot_value_code << '\n';
  stream << "zoom_high_res_region = " << (frozen.config.mode.zoom_high_res_region ? "true" : "false")
         << '\n';
  stream << "zoom_region_file = " << frozen.config.mode.zoom_region_file << '\n';
  stream << "zoom_long_range_strategy = " << zoomLongRangeStrategyToString(frozen.config.mode.zoom_long_range_strategy) << '\n';
  stream << "zoom_region_center_x = " << frozen.config.mode.zoom_region_center_x_mpc_comoving << '\n';
  stream << "zoom_region_center_y = " << frozen.config.mode.zoom_region_center_y_mpc_comoving << '\n';
  stream << "zoom_region_center_z = " << frozen.config.mode.zoom_region_center_z_mpc_comoving << '\n';
  stream << "zoom_region_radius = " << frozen.config.mode.zoom_region_radius_mpc_comoving << '\n';
  stream << "zoom_focused_pm_grid_nx = " << frozen.config.mode.zoom_focused_pm_grid_nx << '\n';
  stream << "zoom_focused_pm_grid_ny = " << frozen.config.mode.zoom_focused_pm_grid_ny << '\n';
  stream << "zoom_focused_pm_grid_nz = " << frozen.config.mode.zoom_focused_pm_grid_nz << '\n';
  stream << "zoom_contamination_radius = " << frozen.config.mode.zoom_contamination_radius_mpc_comoving << '\n';
  stream << "hydro_boundary = " << modeHydroBoundaryToString(frozen.config.mode.hydro_boundary) << '\n';
  stream << "gravity_boundary = " << modeGravityBoundaryToString(frozen.config.mode.gravity_boundary) << '\n';
  stream << "\n[cosmology]\n";
  stream << "omega_matter = " << frozen.config.cosmology.omega_matter << '\n';
  stream << "omega_lambda = " << frozen.config.cosmology.omega_lambda << '\n';
  stream << "omega_baryon = " << frozen.config.cosmology.omega_baryon << '\n';
  stream << "hubble_param = " << frozen.config.cosmology.hubble_param << '\n';
  stream << "sigma8 = " << frozen.config.cosmology.sigma8 << '\n';
  stream << "scalar_index_ns = " << frozen.config.cosmology.scalar_index_ns << '\n';
  stream << "box_size_x = " << frozen.config.cosmology.box_size_x_mpc_comoving << " mpc\n";
  stream << "box_size_y = " << frozen.config.cosmology.box_size_y_mpc_comoving << " mpc\n";
  stream << "box_size_z = " << frozen.config.cosmology.box_size_z_mpc_comoving << " mpc\n";
  stream << "\n[numerics]\n";
  stream << "a_begin = " << bridge_number(frozen.config.numerics.a_begin) << '\n';
  stream << "a_end = " << bridge_number(frozen.config.numerics.a_end) << '\n';
  stream << "z_begin = " << bridge_number(frozen.config.numerics.z_begin) << '\n';
  stream << "z_end = " << bridge_number(frozen.config.numerics.z_end) << '\n';
  stream << "t_code_begin = " << bridge_number(frozen.config.numerics.t_code_begin) << '\n';
  stream << "t_code_end = " << bridge_number(frozen.config.numerics.t_code_end) << '\n';
  stream << "t_phys_begin = " << frozen.config.numerics.t_phys_begin << '\n';
  stream << "t_phys_end = " << frozen.config.numerics.t_phys_end << '\n';
  stream << "integrator_time_variable = "
         << integratorTimeVariableToString(frozen.config.numerics.integrator_time_variable) << '\n';
  stream << "cosmology_max_delta_ln_a = " << frozen.config.numerics.cosmology_max_delta_ln_a << '\n';
  stream << "cosmology_max_hubble_time_fraction = " << frozen.config.numerics.cosmology_max_hubble_time_fraction << '\n';
  stream << "source_max_fractional_change = " << frozen.config.numerics.source_max_fractional_change << '\n';
  stream << "hydro_density_floor_code = " << frozen.config.numerics.hydro_density_floor_code << '\n';
  stream << "hydro_pressure_floor_code = " << frozen.config.numerics.hydro_pressure_floor_code << '\n';
  stream << "hydro_active_batch_max_cells = " << frozen.config.numerics.hydro_active_batch_max_cells << '\n';
  stream << "hydro_reconstruction = "
         << hydroReconstructionModeToString(frozen.config.numerics.hydro_reconstruction) << '\n';
  stream << "hydro_slope_limiter = "
         << hydroSlopeLimiterModeToString(frozen.config.numerics.hydro_slope_limiter) << '\n';
  stream << "hydro_integration_scheme = "
         << hydroIntegrationModeToString(frozen.config.numerics.hydro_integration_scheme) << '\n';
  stream << "hydro_riemann_solver = "
         << hydroRiemannSolverModeToString(frozen.config.numerics.hydro_riemann_solver) << '\n';
  stream << "hydro_dual_energy = "
         << hydroDualEnergyConfigModeToString(frozen.config.numerics.hydro_dual_energy) << '\n';
  stream << "hydro_dual_energy_eta = " << frozen.config.numerics.hydro_dual_energy_eta << '\n';
  stream << "max_global_steps = " << frozen.config.numerics.max_global_steps << '\n';
  stream << "hierarchical_max_rung = " << frozen.config.numerics.hierarchical_max_rung << '\n';
  stream << "hierarchical_coarse_displacement_criterion = "
         << hierarchicalCoarseDisplacementCriterionToString(
                frozen.config.numerics.hierarchical_coarse_displacement_criterion) << '\n';
  stream << "hierarchical_coarse_rms_displacement_factor = "
         << frozen.config.numerics.hierarchical_coarse_rms_displacement_factor << '\n';
  stream << "hierarchical_encounter_timestep_enabled = "
         << (frozen.config.numerics.hierarchical_encounter_timestep_enabled ? "true" : "false") << '\n';
  stream << "gravity_timestep_eta = " << frozen.config.numerics.gravity_timestep_eta << '\n';
  stream << "amr_max_level = " << frozen.config.numerics.amr_max_level << '\n';
  stream << "amr_regrid_mode = " << amrRegridConfigModeToString(frozen.config.numerics.amr_regrid_mode) << '\n';
  stream << "global_step_max_retries = " << frozen.config.numerics.global_step_max_retries << '\n';
  stream << "amr_time_stepping = " << amrTimeSteppingConfigModeToString(frozen.config.numerics.amr_time_stepping) << '\n';
  stream << "amr_predictive_feedback_age_years = " << frozen.config.numerics.amr_predictive_feedback_age_years << '\n';
  stream << "amr_predictive_event_lookahead_years = " << frozen.config.numerics.amr_predictive_event_lookahead_years << '\n';
  stream << "amr_regrid_interval_steps = " << frozen.config.numerics.amr_regrid_interval_steps << '\n';
  stream << "amr_regrid_max_operations = " << frozen.config.numerics.amr_regrid_max_operations << '\n';
  stream << "amr_refine_cell_mass_code = " << frozen.config.numerics.amr_refine_cell_mass_code << '\n';
  stream << "amr_refine_jump_threshold = " << frozen.config.numerics.amr_refine_jump_threshold << '\n';
  stream << "amr_refine_jeans_cells = " << frozen.config.numerics.amr_refine_jeans_cells << '\n';
  stream << "amr_derefine_hysteresis = " << frozen.config.numerics.amr_derefine_hysteresis << '\n';
  stream << "amr_refine_cooling_length_cells = " << frozen.config.numerics.amr_refine_cooling_length_cells << '\n';
  stream << "amr_refine_metallicity_jump = " << frozen.config.numerics.amr_refine_metallicity_jump << '\n';
  stream << "amr_refine_species_jump = " << frozen.config.numerics.amr_refine_species_jump << '\n';
  stream << "amr_refine_vorticity = " << frozen.config.numerics.amr_refine_vorticity << '\n';
  stream << "gravity_softening = " << frozen.config.numerics.gravity_softening_kpc_comoving << " kpc\n";
  auto emitOptionalSoftening = [&](std::string_view key, double value_kpc_comoving) {
    stream << key << " = ";
    if (value_kpc_comoving > 0.0) {
      stream << value_kpc_comoving << " kpc";
    }
    stream << '\n';
  };
  emitOptionalSoftening("gravity_softening_gas", frozen.config.numerics.gravity_softening_gas_kpc_comoving);
  emitOptionalSoftening("gravity_softening_dark_matter", frozen.config.numerics.gravity_softening_dark_matter_kpc_comoving);
  emitOptionalSoftening("gravity_softening_star", frozen.config.numerics.gravity_softening_star_kpc_comoving);
  emitOptionalSoftening("gravity_softening_black_hole", frozen.config.numerics.gravity_softening_black_hole_kpc_comoving);
  emitOptionalSoftening("gravity_softening_tracer", frozen.config.numerics.gravity_softening_tracer_kpc_comoving);
  stream << "gravity_solver = " << gravitySolverToString(frozen.config.numerics.gravity_solver) << '\n';
  stream << "hydro_solver = " << hydroSolverToString(frozen.config.numerics.hydro_solver) << '\n';
  stream << "treepm_pm_grid_nx = " << frozen.config.numerics.treepm_pm_grid_nx << '\n';
  stream << "treepm_pm_grid_ny = " << frozen.config.numerics.treepm_pm_grid_ny << '\n';
  stream << "treepm_pm_grid_nz = " << frozen.config.numerics.treepm_pm_grid_nz << '\n';
  stream << "treepm_asmth_cells = " << frozen.config.numerics.treepm_asmth_cells << '\n';
  stream << "treepm_rcut_cells = " << frozen.config.numerics.treepm_rcut_cells << '\n';
  stream << "treepm_tree_opening_criterion = "
         << treePmOpeningCriterionToString(frozen.config.numerics.treepm_tree_opening_criterion)
         << '\n';
  stream << "treepm_tree_opening_theta = "
         << frozen.config.numerics.treepm_tree_opening_theta << '\n';
  stream << "treepm_tree_relative_force_tolerance = "
         << frozen.config.numerics.treepm_tree_relative_force_tolerance << '\n';
  stream << "treepm_tree_relative_force_acceleration_floor = "
         << frozen.config.numerics.treepm_tree_relative_force_acceleration_floor << '\n';
  stream << "treepm_gaussian_pair_lookup_enabled = "
         << (frozen.config.numerics.treepm_gaussian_pair_lookup_enabled ? "true" : "false") << '\n';
  stream << "treepm_full_residual_lookup_enabled = " << (frozen.config.numerics.treepm_full_residual_lookup_enabled ? "true" : "false") << '\n';
  stream << "treepm_radial_coefficients_enabled = " << (frozen.config.numerics.treepm_radial_coefficients_enabled ? "true" : "false") << '\n';
  stream << "treepm_leaf_simd_enabled = " << (frozen.config.numerics.treepm_leaf_simd_enabled ? "true" : "false") << '\n';
  stream << "treepm_residual_block_size = " << frozen.config.numerics.treepm_residual_block_size << '\n';
  stream << "treepm_tree_max_leaf_size = " << frozen.config.numerics.treepm_tree_max_leaf_size << '\n';
  stream << "treepm_full_mac_diagnostics = "
         << (frozen.config.numerics.treepm_full_mac_diagnostics ? "true" : "false") << '\n';
  stream << "treepm_adaptive_acceptance_enabled = "
         << (frozen.config.numerics.treepm_adaptive_acceptance_enabled ? "true" : "false") << '\n';
  stream << "treepm_adaptive_maximum_opening_angle = "
         << frozen.config.numerics.treepm_adaptive_maximum_opening_angle << '\n';
  stream << "treepm_adaptive_error_estimator = "
         << treePmAdaptiveErrorEstimatorToString(frozen.config.numerics.treepm_adaptive_error_estimator) << '\n';
  stream << "treepm_cutoff_node_policy = "
         << treePmCutoffNodePolicyToString(frozen.config.numerics.treepm_cutoff_node_policy) << '\n';
  stream << "treepm_group_walk_enabled = "
         << (frozen.config.numerics.treepm_group_walk_enabled ? "true" : "false") << '\n';
  stream << "treepm_group_walk_target_size = " << frozen.config.numerics.treepm_group_walk_target_size << '\n';
  stream << "treepm_identical_source_tree_reuse_enabled = "
         << (frozen.config.numerics.treepm_identical_source_tree_reuse_enabled ? "true" : "false") << '\n';
  stream << "treepm_topology_refit_enabled = "
         << (frozen.config.numerics.treepm_topology_refit_enabled ? "true" : "false") << '\n';
  stream << "treepm_assignment_scheme = "
         << (frozen.config.numerics.treepm_assignment_scheme == TreePmAssignmentScheme::kCic ? "cic" : "tsc")
         << '\n';
  stream << "treepm_pm_force_gradient = "
         << (frozen.config.numerics.treepm_pm_force_gradient == TreePmForceGradient::kFiniteDifference4
                 ? "finite_difference_4" : "spectral")
         << '\n';
  stream << "treepm_enable_window_deconvolution = "
         << (frozen.config.numerics.treepm_enable_window_deconvolution ? "true" : "false") << '\n';
  stream << "treepm_update_cadence_steps = " << frozen.config.numerics.treepm_update_cadence_steps << '\n';
  stream << "treepm_pm_reuse_policy = "
         << treePmPmReusePolicyName(frozen.config.numerics.treepm_pm_reuse_policy)
         << '\n';
  stream << "treepm_multirate_rollback_memory_policy = "
         << multirateRollbackMemoryPolicyName(frozen.config.numerics.treepm_multirate_rollback_memory_policy)
         << '\n';
  stream << "treepm_allow_diagnostic_naive_dft = "
         << (frozen.config.numerics.treepm_allow_diagnostic_naive_dft ? "true" : "false") << '\n';
  stream << "treepm_pm_decomposition_mode = "
         << pmDecompositionModeToString(frozen.config.numerics.treepm_pm_decomposition_mode) << '\n';
  stream << "treepm_tree_exchange_batch_bytes = "
         << frozen.config.numerics.treepm_tree_exchange_batch_bytes << '\n';
  stream << "treepm_pm_exchange_batch_bytes = "
         << frozen.config.numerics.treepm_pm_exchange_batch_bytes << '\n';
  stream << "\n[physics]\n";
  stream << "enable_cooling = " << (frozen.config.physics.enable_cooling ? "true" : "false") << '\n';
  stream << "enable_star_formation = "
         << (frozen.config.physics.enable_star_formation ? "true" : "false") << '\n';
  stream << "enable_feedback = " << (frozen.config.physics.enable_feedback ? "true" : "false") << '\n';
  stream << "enable_stellar_evolution = "
         << (frozen.config.physics.enable_stellar_evolution ? "true" : "false") << '\n';
  stream << "reionization_model = " << frozen.config.physics.reionization_model << '\n';
  stream << "uv_background_model = " << uvBackgroundModelToString(frozen.config.physics.uv_background_model) << '\n';
  stream << "self_shielding_model = " << selfShieldingModelToString(frozen.config.physics.self_shielding_model) << '\n';
  stream << "cooling_model = " << coolingModelToString(frozen.config.physics.cooling_model) << '\n';
  stream << "metal_line_table_path = " << frozen.config.physics.metal_line_table_path << '\n';
  stream << "temperature_floor_k = " << frozen.config.physics.temperature_floor_k << '\n';
  stream << "star_formation_model = "
         << starFormationModelKindToString(frozen.config.physics.star_formation_model) << '\n';
  stream << "sf_density_threshold_code = " << frozen.config.physics.sf_density_threshold_code << '\n';
  stream << "sf_temperature_threshold_k = " << frozen.config.physics.sf_temperature_threshold_k << '\n';
  stream << "sf_min_converging_flow_rate_code = " << frozen.config.physics.sf_min_converging_flow_rate_code
         << '\n';
  stream << "sf_epsilon_ff = " << frozen.config.physics.sf_epsilon_ff << '\n';
  stream << "sf_bound_alpha_vir_max = " << frozen.config.physics.sf_bound_alpha_vir_max << '\n';
  stream << "sf_require_converging_flow = "
         << (frozen.config.physics.sf_require_converging_flow ? "true" : "false") << '\n';
  stream << "sf_collapse_timescale = "
         << starFormationCollapseTimescaleToString(frozen.config.physics.sf_collapse_timescale) << '\n';
  stream << "sf_jeans_mass_floor_code = " << frozen.config.physics.sf_jeans_mass_floor_code << '\n';
  stream << "sf_star_particle_mass_policy = "
         << starParticleMassPolicyToString(frozen.config.physics.sf_star_particle_mass_policy) << '\n';
  stream << "sf_target_star_particle_mass_code = "
         << frozen.config.physics.sf_target_star_particle_mass_code << '\n';
  stream << "sf_target_star_particle_mass_fraction = "
         << frozen.config.physics.sf_target_star_particle_mass_fraction << '\n';
  stream << "sf_min_star_particle_mass_code = "
         << frozen.config.physics.sf_min_star_particle_mass_code << '\n';
  stream << "sf_max_star_particle_mass_code = "
         << frozen.config.physics.sf_max_star_particle_mass_code << '\n';
  stream << "sf_max_spawn_particles_per_cell_step = "
         << frozen.config.physics.sf_max_spawn_particles_per_cell_step << '\n';
  stream << "sf_max_fractional_mass_conversion = "
         << frozen.config.physics.sf_max_fractional_mass_conversion << '\n';
  stream << "sf_min_remaining_gas_fraction = "
         << frozen.config.physics.sf_min_remaining_gas_fraction << '\n';
  stream << "sf_min_remaining_gas_mass_code = "
         << frozen.config.physics.sf_min_remaining_gas_mass_code << '\n';
  stream << "sf_temperature_safety_ceiling_k = "
         << frozen.config.physics.sf_temperature_safety_ceiling_k << '\n';
  stream << "sf_stochastic_spawning = "
         << (frozen.config.physics.sf_stochastic_spawning ? "true" : "false") << '\n';
  stream << "sf_random_seed = " << frozen.config.physics.sf_random_seed << '\n';
  stream << "sf_effective_parameter_set = " << frozen.config.physics.sf_effective_parameter_set << '\n';
  stream << "sf_effective_n_h_threshold = " << frozen.config.physics.sf_effective_n_h_threshold_cgs << '\n';
  stream << "sf_effective_min_baryon_overdensity = " << frozen.config.physics.sf_effective_min_baryon_overdensity << '\n';
  stream << "sf_effective_t_star_at_threshold = " << frozen.config.physics.sf_effective_t_star_at_threshold_code << " code\n";
  stream << "sf_effective_evaporation_factor_at_threshold = " << frozen.config.physics.sf_effective_evaporation_factor_at_threshold << '\n';
  stream << "sf_effective_evaporation_density_exponent = " << frozen.config.physics.sf_effective_evaporation_density_exponent << '\n';
  stream << "sf_effective_cold_phase_temperature = " << frozen.config.physics.sf_effective_cold_phase_temperature_k << '\n';
  stream << "sf_effective_supernova_specific_energy_code = " << frozen.config.physics.sf_effective_supernova_specific_energy_code << '\n';
  stream << "sf_effective_massive_star_fraction = " << frozen.config.physics.sf_effective_massive_star_fraction << '\n';
  stream << "sf_effective_q_eos = " << frozen.config.physics.sf_effective_q_eos << '\n';
  stream << "sf_effective_isothermal_temperature = " << frozen.config.physics.sf_effective_isothermal_temperature_k << '\n';
  stream << "sf_effective_hot_excess_tolerance = " << frozen.config.physics.sf_effective_hot_excess_tolerance << '\n';
  stream << "sf_effective_eos_relaxation = " << effectiveIsmEosRelaxationToString(frozen.config.physics.sf_effective_eos_relaxation) << '\n';
  stream << "sf_effective_eos_relaxation_timescale = " << frozen.config.physics.sf_effective_eos_relaxation_timescale_code << " code\n";
  stream << "sf_effective_eos_table_bins = " << frozen.config.physics.sf_effective_eos_table_bins << '\n';
  stream << "sf_effective_eos_max_density_ratio = " << frozen.config.physics.sf_effective_eos_max_density_ratio << '\n';
  stream << "sf_effective_birth_mass_convention = " << effectiveIsmBirthMassConventionToString(frozen.config.physics.sf_effective_birth_mass_convention) << '\n';
  stream << "sf_effective_feedback_coupling = " << effectiveIsmFeedbackCouplingToString(frozen.config.physics.sf_effective_feedback_coupling) << '\n';
  stream << "fb_mode = " << feedbackModeToString(frozen.config.physics.fb_mode) << '\n';
  stream << "fb_variant = " << feedbackVariantToString(frozen.config.physics.fb_variant) << '\n';
  stream << "fb_use_returned_mass_budget = "
         << (frozen.config.physics.fb_use_returned_mass_budget ? "true" : "false") << '\n';
  stream << "fb_epsilon_thermal = " << frozen.config.physics.fb_epsilon_thermal << '\n';
  stream << "fb_epsilon_kinetic = " << frozen.config.physics.fb_epsilon_kinetic << '\n';
  stream << "fb_epsilon_momentum = " << frozen.config.physics.fb_epsilon_momentum << '\n';
  stream << "fb_sn_energy_erg_per_mass_code = " << frozen.config.physics.fb_sn_energy_erg_per_mass_code << '\n';
  stream << "fb_momentum_code_per_mass_code = " << frozen.config.physics.fb_momentum_code_per_mass_code << '\n';
  stream << "fb_neighbor_count = " << frozen.config.physics.fb_neighbor_count << '\n';
  stream << "fb_delayed_cooling_time_code = " << frozen.config.physics.fb_delayed_cooling_time_code << '\n';
  stream << "fb_stochastic_event_probability = " << frozen.config.physics.fb_stochastic_event_probability << '\n';
  stream << "fb_random_seed = " << frozen.config.physics.fb_random_seed << '\n';
  stream << "stellar_evolution_table_path = " << frozen.config.physics.stellar_evolution_table_path << '\n';
  stream << "stellar_evolution_hubble_time_years = "
         << frozen.config.physics.stellar_evolution_hubble_time_years << '\n';
  stream << "stellar_evolution_require_production_table = "
         << (frozen.config.physics.stellar_evolution_require_production_table ? "true" : "false") << '\n';
  stream << "metal_species_mode = "
         << metalSpeciesModeToString(frozen.config.physics.metal_species_mode) << '\n';
  stream << "enable_metal_diffusion = "
         << (frozen.config.physics.enable_metal_diffusion ? "true" : "false") << '\n';
  stream << "metal_diffusion_model = "
         << metalDiffusionModelToString(frozen.config.physics.metal_diffusion_model) << '\n';
  stream << "metal_diffusion_time_integrator = "
         << metalDiffusionTimeIntegratorToString(
                frozen.config.physics.metal_diffusion_time_integrator) << '\n';
  stream << "metal_diffusion_coefficient = "
         << frozen.config.physics.metal_diffusion_coefficient << '\n';
  stream << "metal_diffusion_cfl = " << frozen.config.physics.metal_diffusion_cfl << '\n';
  stream << "metal_diffusion_max_subcycles = "
         << frozen.config.physics.metal_diffusion_max_subcycles << '\n';
  stream << "metal_diffusion_max_rkl_stages = "
         << frozen.config.physics.metal_diffusion_max_rkl_stages << '\n';
  stream << "metal_diffusion_coefficient_floor_code = "
         << frozen.config.physics.metal_diffusion_coefficient_floor_code << '\n';
  stream << "metal_diffusion_coefficient_ceiling_code = "
         << frozen.config.physics.metal_diffusion_coefficient_ceiling_code << '\n';
  stream << "enable_black_hole_agn = " << (frozen.config.physics.enable_black_hole_agn ? "true" : "false")
         << '\n';
  stream << "bh_enable_seeding = " << (frozen.config.physics.bh_enable_seeding ? "true" : "false") << '\n';
  stream << "bh_seed_halo_mass_threshold_code = " << frozen.config.physics.bh_seed_halo_mass_threshold_code
         << '\n';
  stream << "bh_seed_mass_code = " << frozen.config.physics.bh_seed_mass_code << '\n';
  stream << "bh_seed_max_per_cell = " << frozen.config.physics.bh_seed_max_per_cell << '\n';
  stream << "bh_alpha_bondi = " << frozen.config.physics.bh_alpha_bondi << '\n';
  stream << "bh_use_eddington_cap = " << (frozen.config.physics.bh_use_eddington_cap ? "true" : "false")
         << '\n';
  stream << "bh_epsilon_r = " << frozen.config.physics.bh_epsilon_r << '\n';
  stream << "bh_epsilon_f = " << frozen.config.physics.bh_epsilon_f << '\n';
  stream << "bh_feedback_coupling_efficiency = "
         << frozen.config.physics.bh_feedback_coupling_efficiency << '\n';
  stream << "bh_duty_cycle_active_edd_ratio_threshold = "
         << frozen.config.physics.bh_duty_cycle_active_edd_ratio_threshold << '\n';
  stream << "bh_proton_mass_si = " << frozen.config.physics.bh_proton_mass_si << '\n';
  stream << "bh_thomson_cross_section_si = " << frozen.config.physics.bh_thomson_cross_section_si << '\n';
  stream << "bh_newton_g_si = " << frozen.config.physics.bh_newton_g_si << '\n';
  stream << "bh_speed_of_light_si = " << frozen.config.physics.bh_speed_of_light_si << '\n';
  stream << "enable_tracers = " << (frozen.config.physics.enable_tracers ? "true" : "false") << '\n';
  stream << "tracer_track_mass = " << (frozen.config.physics.tracer_track_mass ? "true" : "false") << '\n';
  stream << "tracer_min_host_mass_code = " << frozen.config.physics.tracer_min_host_mass_code << '\n';
  stream << "thermochemistry_backend = " << thermochemistryBackendModeToString(frozen.config.physics.thermochemistry_backend) << '\n';
  stream << "metal_cooling_mode = " << metalCoolingModeConfigToString(frozen.config.physics.metal_cooling_mode) << '\n';
  stream << "metal_cooling_table_range = " << metalCoolingTableRangeConfigToString(frozen.config.physics.metal_cooling_table_range) << '\n';
  stream << "dust_model = " << dustModelConfigToString(frozen.config.physics.dust_model) << '\n';
  stream << "ism_regime_model = " << ismRegimeModelConfigToString(frozen.config.physics.ism_regime_model) << '\n';
  stream << "sgs_turbulence_model = " << sgsTurbulenceModelConfigToString(frozen.config.physics.sgs_turbulence_model) << '\n';
  stream << "chemistry_adaptive_selection = " << chemistryAdaptiveSelectionConfigToString(frozen.config.physics.chemistry_adaptive_selection) << '\n';
  stream << "stellar_population_model = " << stellarPopulationModelConfigToString(frozen.config.physics.stellar_population_model) << '\n';
  stream << "fb_supernova_model = " << supernovaFeedbackModelConfigToString(frozen.config.physics.fb_supernova_model) << '\n';
  stream << "thermochemistry_relative_tolerance = " << frozen.config.physics.thermochemistry_relative_tolerance << '\n';
  stream << "thermochemistry_abundance_absolute_tolerance = " << frozen.config.physics.thermochemistry_abundance_absolute_tolerance << '\n';
  stream << "thermochemistry_analytic_relaxation_time_years = " << frozen.config.physics.thermochemistry_analytic_relaxation_time_years << '\n';
  stream << "thermochemistry_analytic_equilibrium_temperature_k = " << frozen.config.physics.thermochemistry_analytic_equilibrium_temperature_k << '\n';
  stream << "isrf_g0_habing = " << frozen.config.physics.isrf_g0_habing << '\n';
  stream << "chemistry_shielding_length_max_pc = " << frozen.config.physics.chemistry_shielding_length_max_pc << '\n';
  stream << "metal_cooling_switch_temperature_k = " << frozen.config.physics.metal_cooling_switch_temperature_k << '\n';
  stream << "chemistry_primordial_helium_mass_fraction = " << frozen.config.physics.chemistry_primordial_helium_mass_fraction << '\n';
  stream << "dust_accretion_timescale_ref_years = " << frozen.config.physics.dust_accretion_timescale_ref_years << '\n';
  stream << "dust_sputtering_timescale_ref_years = " << frozen.config.physics.dust_sputtering_timescale_ref_years << '\n';
  stream << "dust_sn_swept_mass_msun = " << frozen.config.physics.dust_sn_swept_mass_msun << '\n';
  stream << "dust_sn_destruction_efficiency = " << frozen.config.physics.dust_sn_destruction_efficiency << '\n';
  stream << "dust_agb_condensation_efficiency = " << frozen.config.physics.dust_agb_condensation_efficiency << '\n';
  stream << "dust_cc_condensation_efficiency = " << frozen.config.physics.dust_cc_condensation_efficiency << '\n';
  stream << "ism_jeans_resolution_cells = " << frozen.config.physics.ism_jeans_resolution_cells << '\n';
  stream << "ism_regime_hysteresis_fraction = " << frozen.config.physics.ism_regime_hysteresis_fraction << '\n';
  stream << "ism_regime_min_dwell_years = " << frozen.config.physics.ism_regime_min_dwell_years << '\n';
  stream << "sgs_dissipation_coefficient = " << frozen.config.physics.sgs_dissipation_coefficient << '\n';
  stream << "sgs_eddy_viscosity_coefficient = " << frozen.config.physics.sgs_eddy_viscosity_coefficient << '\n';
  stream << "sgs_feedback_energy_fraction = " << frozen.config.physics.sgs_feedback_energy_fraction << '\n';
  stream << "chemistry_adaptive_hysteresis_fraction = " << frozen.config.physics.chemistry_adaptive_hysteresis_fraction << '\n';
  stream << "chemistry_adaptive_min_dwell_years = " << frozen.config.physics.chemistry_adaptive_min_dwell_years << '\n';
  stream << "chemistry_adaptive_cold_temperature_k = " << frozen.config.physics.chemistry_adaptive_cold_temperature_k << '\n';
  stream << "chemistry_adaptive_cooling_time_ratio = " << frozen.config.physics.chemistry_adaptive_cooling_time_ratio << '\n';
  stream << "chemistry_adaptive_dense_nh_cm3 = " << frozen.config.physics.chemistry_adaptive_dense_nh_cm3 << '\n';
  stream << "chemistry_shadow_sample_fraction = " << frozen.config.physics.chemistry_shadow_sample_fraction << '\n';
  stream << "ssp_imf_mass_min_msun = " << frozen.config.physics.ssp_imf_mass_min_msun << '\n';
  stream << "ssp_imf_mass_max_msun = " << frozen.config.physics.ssp_imf_mass_max_msun << '\n';
  stream << "ssp_cc_mass_min_msun = " << frozen.config.physics.ssp_cc_mass_min_msun << '\n';
  stream << "ssp_cc_mass_max_msun = " << frozen.config.physics.ssp_cc_mass_max_msun << '\n';
  stream << "ssp_snia_per_msun = " << frozen.config.physics.ssp_snia_per_msun << '\n';
  stream << "ssp_snia_delay_min_years = " << frozen.config.physics.ssp_snia_delay_min_years << '\n';
  stream << "ssp_snia_delay_max_years = " << frozen.config.physics.ssp_snia_delay_max_years << '\n';
  stream << "ssp_energy_per_sn_erg = " << frozen.config.physics.ssp_energy_per_sn_erg << '\n';
  stream << "ssp_reference_net_metal_yield = " << frozen.config.physics.ssp_reference_net_metal_yield << '\n';
  stream << "fb_presn_wind_momentum_lsun_per_c = " << frozen.config.physics.fb_presn_wind_momentum_lsun_per_c << '\n';
  stream << "thermochemistry_max_substeps = " << frozen.config.physics.thermochemistry_max_substeps << '\n';
  stream << "thermochemistry_max_newton_iterations = " << frozen.config.physics.thermochemistry_max_newton_iterations << '\n';
  stream << "chemistry_shadow_max_cells = " << frozen.config.physics.chemistry_shadow_max_cells << '\n';
  stream << "uv_background_table_path = " << frozen.config.physics.uv_background_table_path << '\n';
  stream << "uv_background_table_sha256 = " << frozen.config.physics.uv_background_table_sha256 << '\n';
  stream << "metal_cooling_table_path = " << frozen.config.physics.metal_cooling_table_path << '\n';
  stream << "metal_cooling_table_sha256 = " << frozen.config.physics.metal_cooling_table_sha256 << '\n';
  stream << "metal_cooling_accept_uv_background_mismatch = "
         << (frozen.config.physics.metal_cooling_accept_uv_background_mismatch ? "true" : "false") << '\n';
  stream << "ssp_yield_table_path = " << frozen.config.physics.ssp_yield_table_path << '\n';
  stream << "ssp_yield_table_sha256 = " << frozen.config.physics.ssp_yield_table_sha256 << '\n';
  stream << "chemistry_metallicity_scaled_dust_proxy = " << (frozen.config.physics.chemistry_metallicity_scaled_dust_proxy ? "true" : "false") << '\n';
  stream << "thermochemistry_allow_unqualified_rates = " << (frozen.config.physics.thermochemistry_allow_unqualified_rates ? "true" : "false") << '\n';
  stream << "ssp_stochastic_sn_events = " << (frozen.config.physics.ssp_stochastic_sn_events ? "true" : "false") << '\n';
  stream << "ssp_clustered_associations = " << (frozen.config.physics.ssp_clustered_associations ? "true" : "false") << '\n';
  stream << "ssp_association_radius_code = " << frozen.config.physics.ssp_association_radius_code << '\n';
  stream << "ssp_association_time_window_years = " << frozen.config.physics.ssp_association_time_window_years << '\n';
  stream << "ssp_wind_model = " << sspWindModelConfigToString(frozen.config.physics.ssp_wind_model) << '\n';
  stream << "ssp_wind_table_path = " << frozen.config.physics.ssp_wind_table_path << '\n';
  stream << "ssp_wind_table_sha256 = " << frozen.config.physics.ssp_wind_table_sha256 << '\n';
  stream << "\n[output]\n";
  stream << "run_name = " << frozen.config.output.run_name << '\n';
  stream << "output_directory = " << frozen.config.output.output_directory << '\n';
  stream << "output_stem = " << frozen.config.output.output_stem << '\n';
  stream << "restart_stem = " << frozen.config.output.restart_stem << '\n';
  stream << "snapshot_interval_steps = " << frozen.config.output.snapshot_interval_steps << '\n';
  stream << "snapshot_interval_time_code = "
         << frozen.config.output.snapshot_interval_time_code << '\n';
  stream << "snapshot_layout = "
         << scienceSnapshotLayoutLabel(frozen.config.output.snapshot_layout) << '\n';
  stream << "snapshot_num_files = "
         << frozen.config.output.snapshot_num_files << '\n';
  stream << "write_restarts = " << (frozen.config.output.write_restarts ? "true" : "false") << '\n';
  stream << "restart_interval_steps = " << frozen.config.output.restart_interval_steps << '\n';
  stream << "restart_interval_time_code = " << frozen.config.output.restart_interval_time_code << '\n';
  stream << "snapshot_verification = "
         << snapshotVerificationPolicyLabel(frozen.config.output.snapshot_verification) << '\n';
  stream << "\n[parallel]\n";
  stream << "mpi_ranks_expected = " << frozen.config.parallel.mpi_ranks_expected << '\n';
  stream << "omp_threads = " << frozen.config.parallel.omp_threads << '\n';
  stream << "gpu_devices = " << frozen.config.parallel.gpu_devices << '\n';
  stream << "deterministic_reduction = "
         << (frozen.config.parallel.deterministic_reduction ? "true" : "false") << '\n';
  stream << "decomposition_particle_count_weight = "
         << frozen.config.parallel.decomposition_particle_count_weight << '\n';
  stream << "decomposition_gas_cell_weight = "
         << frozen.config.parallel.decomposition_gas_cell_weight << '\n';
  stream << "decomposition_tree_interaction_weight = "
         << frozen.config.parallel.decomposition_tree_interaction_weight << '\n';
  stream << "decomposition_pm_mesh_weight = "
         << frozen.config.parallel.decomposition_pm_mesh_weight << '\n';
  stream << "decomposition_amr_patch_weight = "
         << frozen.config.parallel.decomposition_amr_patch_weight << '\n';
  stream << "decomposition_active_fraction_weight = "
         << frozen.config.parallel.decomposition_active_fraction_weight << '\n';
  stream << "decomposition_memory_pressure_weight = "
         << frozen.config.parallel.decomposition_memory_pressure_weight << '\n';
  stream << "decomposition_gpu_occupancy_weight = "
         << frozen.config.parallel.decomposition_gpu_occupancy_weight << '\n';
  stream << "decomposition_generic_work_weight = "
         << frozen.config.parallel.decomposition_generic_work_weight << '\n';
  stream << "decomposition_runtime_rebalance_enabled = "
         << (frozen.config.parallel.decomposition_runtime_rebalance_enabled ? "true" : "false") << '\n';
  stream << "decomposition_spatial_work_enabled = "
         << (frozen.config.parallel.decomposition_spatial_work_enabled ? "true" : "false") << '\n';
  stream << "decomposition_debug_exact_ownership_audit = "
         << (frozen.config.parallel.decomposition_debug_exact_ownership_audit ? "true" : "false") << '\n';
  stream << "decomposition_rebalance_imbalance_trigger = "
         << frozen.config.parallel.decomposition_rebalance_imbalance_trigger << '\n';
  stream << "decomposition_rebalance_memory_trigger = "
         << frozen.config.parallel.decomposition_rebalance_memory_trigger << '\n';
  stream << "decomposition_rebalance_max_migrated_load_fraction = "
         << frozen.config.parallel.decomposition_rebalance_max_migrated_load_fraction << '\n';
  stream << "decomposition_measured_tree_pair_weight = "
         << frozen.config.parallel.decomposition_measured_tree_pair_weight << '\n';
  stream << "decomposition_measured_pm_cell_weight = "
         << frozen.config.parallel.decomposition_measured_pm_cell_weight << '\n';
  stream << "decomposition_measured_amr_cell_weight = "
         << frozen.config.parallel.decomposition_measured_amr_cell_weight << '\n';
  stream << "decomposition_measured_hydro_face_weight = "
         << frozen.config.parallel.decomposition_measured_hydro_face_weight << '\n';
  stream << "decomposition_measured_wall_ms_weight = "
         << frozen.config.parallel.decomposition_measured_wall_ms_weight << '\n';
  stream << "isolated_pm_root_workspace_limit_bytes = "
         << frozen.config.parallel.isolated_pm_root_workspace_limit_bytes << '\n';
  stream << "zoom_high_res_allgather_limit_bytes = "
         << frozen.config.parallel.zoom_high_res_allgather_limit_bytes << '\n';
  stream << "gravity_memory_budget_bytes = "
         << frozen.config.parallel.gravity_memory_budget_bytes << '\n';
  stream << "gravity_backend_unknown_reserve_bytes = "
         << frozen.config.parallel.gravity_backend_unknown_reserve_bytes << '\n';
  stream << "gravity_memory_safety_margin_fraction = "
         << frozen.config.parallel.gravity_memory_safety_margin_fraction << '\n';
  stream << "process_memory_budget_bytes = "
         << frozen.config.parallel.process_memory_budget_bytes << '\n';
  stream << "process_mpi_unknown_reserve_bytes = "
         << frozen.config.parallel.process_mpi_unknown_reserve_bytes << '\n';
  stream << "process_fftw_unknown_reserve_bytes = "
         << frozen.config.parallel.process_fftw_unknown_reserve_bytes << '\n';
  stream << "process_hdf5_unknown_reserve_bytes = "
         << frozen.config.parallel.process_hdf5_unknown_reserve_bytes << '\n';
  stream << "process_allocator_unknown_reserve_bytes = "
         << frozen.config.parallel.process_allocator_unknown_reserve_bytes << '\n';
  stream << "process_output_restart_overlap_bytes = "
         << frozen.config.parallel.process_output_restart_overlap_bytes << '\n';
  stream << "process_memory_safety_margin_fraction = "
         << frozen.config.parallel.process_memory_safety_margin_fraction << '\n';
  stream << "node_memory_admission = "
         << nodeMemoryAdmissionModeToString(frozen.config.parallel.node_memory_admission) << '\n';
  stream << "node_memory_budget_bytes = "
         << frozen.config.parallel.node_memory_budget_bytes << '\n';
  stream << "node_memory_reserve_bytes = "
         << frozen.config.parallel.node_memory_reserve_bytes << '\n';
  stream << "\n[analysis]\n";
  stream << "enable_diagnostics = " << (frozen.config.analysis.enable_diagnostics ? "true" : "false")
         << '\n';
  stream << "enable_halo_workflow = " << (frozen.config.analysis.enable_halo_workflow ? "true" : "false")
         << '\n';
  stream << "halo_on_the_fly = " << (frozen.config.analysis.halo_on_the_fly ? "true" : "false") << '\n';
  stream << "diagnostics_execution_policy = "
         << diagnosticsExecutionPolicyToString(frozen.config.analysis.diagnostics_execution_policy) << '\n';
  stream << "run_health_interval_steps = " << frozen.config.analysis.run_health_interval_steps << '\n';
  stream << "science_light_interval_steps = " << frozen.config.analysis.science_light_interval_steps << '\n';
  stream << "science_heavy_interval_steps = " << frozen.config.analysis.science_heavy_interval_steps << '\n';
  stream << "retention_bundle_count = " << frozen.config.analysis.retention_bundle_count << '\n';
  stream << "power_spectrum_mesh_n = " << frozen.config.analysis.power_spectrum_mesh_n << '\n';
  stream << "power_spectrum_bin_count = " << frozen.config.analysis.power_spectrum_bin_count << '\n';
  stream << "sf_history_bin_count = " << frozen.config.analysis.sf_history_bin_count << '\n';
  stream << "quicklook_grid_n = " << frozen.config.analysis.quicklook_grid_n << '\n';
  stream << "diagnostics_stem = " << frozen.config.analysis.diagnostics_stem << '\n';
  stream << "halo_catalog_stem = " << frozen.config.analysis.halo_catalog_stem << '\n';
  stream << "merger_tree_stem = " << frozen.config.analysis.merger_tree_stem << '\n';
  stream << "halo_fof_linking_length_factor = " << frozen.config.analysis.halo_fof_linking_length_factor << '\n';
  stream << "halo_fof_min_group_size = " << frozen.config.analysis.halo_fof_min_group_size << '\n';
  stream << "halo_include_gas = " << (frozen.config.analysis.halo_include_gas ? "true" : "false") << '\n';
  stream << "halo_include_stars = " << (frozen.config.analysis.halo_include_stars ? "true" : "false") << '\n';
  stream << "halo_include_black_holes = " << (frozen.config.analysis.halo_include_black_holes ? "true" : "false")
         << '\n';
  stream << "\n[compatibility]\n";
  stream << "allow_unknown_keys = "
         << (frozen.config.compatibility.allow_unknown_keys ? "true" : "false") << '\n';
  return stream.str();
}

[[nodiscard]] FrozenConfig normalizeValidateFreeze(
    const std::map<std::string, ParsedEntry>& parsed_entries,
    const std::string& raw_text,
    const std::string& source_name,
    const ParseOptions& options) {
  std::map<std::string, ParsedEntry> entries = parsed_entries;
  std::set<std::string> consumed;
  FrozenConfig frozen;
  frozen.raw_text = raw_text;
  frozen.user_config.source_name = source_name;
  frozen.provenance.source_name = source_name;
  validateRequiredSchema(entries);
  applyDeprecatedAliases(entries, consumed, frozen);

  frozen.config.schema_version = parseNumber<int>(
      requireString(entries, consumed, "schema_version", "1"), "schema_version");

  frozen.config.units.length_unit =
      toLower(requireString(entries, consumed, "units.length_unit", frozen.config.units.length_unit));
  frozen.config.units.mass_unit =
      toLower(requireString(entries, consumed, "units.mass_unit", frozen.config.units.mass_unit));
  frozen.config.units.velocity_unit = toLower(
      requireString(entries, consumed, "units.velocity_unit", frozen.config.units.velocity_unit));
  frozen.config.units.coordinate_frame = parseCoordinateFrame(
      requireString(entries, consumed, "units.coordinate_frame", defaultFor("units.coordinate_frame")));

  frozen.config.mode.mode = parseMode(requireString(entries, consumed, "mode.mode", defaultFor("mode.mode")));
  frozen.config.mode.ic_file = requireString(entries, consumed, "mode.ic_file", defaultFor("mode.ic_file"));
  const bool has_explicit_ic_convention = entries.contains("mode.ic_convention");
  frozen.config.mode.ic_convention = parseInitialConditionConvention(
      requireString(
          entries, consumed, "mode.ic_convention",
          defaultFor("mode.ic_convention")));
  frozen.config.mode.ic_manifest_file = requireString(
      entries, consumed, "mode.ic_manifest_file",
      defaultFor("mode.ic_manifest_file"));
  frozen.config.mode.ic_bridge_source_length_unit_to_si = parseFloating(
      requireString(
          entries, consumed, "mode.ic_bridge_source_length_unit_to_si",
          defaultFor("mode.ic_bridge_source_length_unit_to_si")),
      "mode.ic_bridge_source_length_unit_to_si");
  frozen.config.mode.ic_bridge_source_mass_unit_to_si = parseFloating(
      requireString(
          entries, consumed, "mode.ic_bridge_source_mass_unit_to_si",
          defaultFor("mode.ic_bridge_source_mass_unit_to_si")),
      "mode.ic_bridge_source_mass_unit_to_si");
  frozen.config.mode.ic_bridge_source_velocity_unit_to_si = parseFloating(
      requireString(
          entries, consumed, "mode.ic_bridge_source_velocity_unit_to_si",
          defaultFor("mode.ic_bridge_source_velocity_unit_to_si")),
      "mode.ic_bridge_source_velocity_unit_to_si");
  frozen.config.mode.ic_bridge_coordinate_frame =
      parseInitialConditionCoordinateFrame(requireString(
          entries, consumed, "mode.ic_bridge_coordinate_frame",
          defaultFor("mode.ic_bridge_coordinate_frame")));
  frozen.config.mode.ic_bridge_velocity_convention =
      parseInitialConditionVelocityConvention(requireString(
          entries, consumed, "mode.ic_bridge_velocity_convention",
          defaultFor("mode.ic_bridge_velocity_convention")));
  frozen.config.mode.ic_bridge_length_hubble_exponent = parseOptionalFinite(
      requireString(
          entries, consumed, "mode.ic_bridge_length_hubble_exponent",
          defaultFor("mode.ic_bridge_length_hubble_exponent")),
      "mode.ic_bridge_length_hubble_exponent");
  frozen.config.mode.ic_bridge_length_scale_factor_exponent = parseOptionalFinite(
      requireString(
          entries, consumed, "mode.ic_bridge_length_scale_factor_exponent",
          defaultFor("mode.ic_bridge_length_scale_factor_exponent")),
      "mode.ic_bridge_length_scale_factor_exponent");
  frozen.config.mode.ic_bridge_mass_hubble_exponent = parseOptionalFinite(
      requireString(
          entries, consumed, "mode.ic_bridge_mass_hubble_exponent",
          defaultFor("mode.ic_bridge_mass_hubble_exponent")),
      "mode.ic_bridge_mass_hubble_exponent");
  frozen.config.mode.ic_bridge_mass_scale_factor_exponent = parseOptionalFinite(
      requireString(
          entries, consumed, "mode.ic_bridge_mass_scale_factor_exponent",
          defaultFor("mode.ic_bridge_mass_scale_factor_exponent")),
      "mode.ic_bridge_mass_scale_factor_exponent");
  frozen.config.mode.ic_bridge_velocity_hubble_exponent = parseOptionalFinite(
      requireString(
          entries, consumed, "mode.ic_bridge_velocity_hubble_exponent",
          defaultFor("mode.ic_bridge_velocity_hubble_exponent")),
      "mode.ic_bridge_velocity_hubble_exponent");
  frozen.config.mode.ic_bridge_velocity_scale_factor_exponent = parseOptionalFinite(
      requireString(
          entries, consumed, "mode.ic_bridge_velocity_scale_factor_exponent",
          defaultFor("mode.ic_bridge_velocity_scale_factor_exponent")),
      "mode.ic_bridge_velocity_scale_factor_exponent");
  frozen.config.mode.ic_chunk_particle_count = parseNumber<std::uint64_t>(
      requireString(
          entries, consumed, "mode.ic_chunk_particle_count",
          defaultFor("mode.ic_chunk_particle_count")),
      "mode.ic_chunk_particle_count");
  frozen.config.mode.ic_staging_particle_count = parseNumber<std::uint64_t>(
      requireString(
          entries, consumed, "mode.ic_staging_particle_count",
          defaultFor("mode.ic_staging_particle_count")),
      "mode.ic_staging_particle_count");
  frozen.config.mode.ic_part_type2_policy =
      parseInitialConditionSpeciesPolicy(
          requireString(
              entries, consumed, "mode.ic_part_type2_policy",
              defaultFor("mode.ic_part_type2_policy")),
          "mode.ic_part_type2_policy");
  frozen.config.mode.ic_part_type3_policy =
      parseInitialConditionSpeciesPolicy(
          requireString(
              entries, consumed, "mode.ic_part_type3_policy",
              defaultFor("mode.ic_part_type3_policy")),
          "mode.ic_part_type3_policy");
  frozen.config.mode.ic_gas_internal_energy_policy =
      parseInitialConditionMissingFieldPolicy(
          requireString(
              entries, consumed, "mode.ic_gas_internal_energy_policy",
              defaultFor("mode.ic_gas_internal_energy_policy")),
          "mode.ic_gas_internal_energy_policy");
  frozen.config.mode.ic_gas_density_policy =
      parseInitialConditionMissingFieldPolicy(
          requireString(
              entries, consumed, "mode.ic_gas_density_policy",
              defaultFor("mode.ic_gas_density_policy")),
          "mode.ic_gas_density_policy");
  frozen.config.mode.ic_star_formation_time_policy =
      parseInitialConditionMissingFieldPolicy(
          requireString(
              entries, consumed, "mode.ic_star_formation_time_policy",
              defaultFor("mode.ic_star_formation_time_policy")),
          "mode.ic_star_formation_time_policy");
  frozen.config.mode.ic_star_initial_mass_policy =
      parseInitialConditionMissingFieldPolicy(
          requireString(
              entries, consumed, "mode.ic_star_initial_mass_policy",
              defaultFor("mode.ic_star_initial_mass_policy")),
          "mode.ic_star_initial_mass_policy");
  frozen.config.mode.ic_star_metallicity_policy =
      parseInitialConditionMissingFieldPolicy(
          requireString(
              entries, consumed, "mode.ic_star_metallicity_policy",
              defaultFor("mode.ic_star_metallicity_policy")),
          "mode.ic_star_metallicity_policy");
  frozen.config.mode.ic_bh_mdot_policy =
      parseInitialConditionMissingFieldPolicy(
          requireString(
              entries, consumed, "mode.ic_bh_mdot_policy",
              defaultFor("mode.ic_bh_mdot_policy")),
          "mode.ic_bh_mdot_policy");
  frozen.config.mode.ic_gas_internal_energy_value_code = parseFloating(
      requireString(
          entries, consumed, "mode.ic_gas_internal_energy_value_code",
          defaultFor("mode.ic_gas_internal_energy_value_code")),
      "mode.ic_gas_internal_energy_value_code");
  frozen.config.mode.ic_gas_density_value_code = parseFloating(
      requireString(
          entries, consumed, "mode.ic_gas_density_value_code",
          defaultFor("mode.ic_gas_density_value_code")),
      "mode.ic_gas_density_value_code");
  frozen.config.mode.ic_star_formation_time_value = parseFloating(
      requireString(
          entries, consumed, "mode.ic_star_formation_time_value",
          defaultFor("mode.ic_star_formation_time_value")),
      "mode.ic_star_formation_time_value");
  frozen.config.mode.ic_star_initial_mass_value_code = parseFloating(
      requireString(
          entries, consumed, "mode.ic_star_initial_mass_value_code",
          defaultFor("mode.ic_star_initial_mass_value_code")),
      "mode.ic_star_initial_mass_value_code");
  frozen.config.mode.ic_star_metallicity_value = parseFloating(
      requireString(
          entries, consumed, "mode.ic_star_metallicity_value",
          defaultFor("mode.ic_star_metallicity_value")),
      "mode.ic_star_metallicity_value");
  frozen.config.mode.ic_bh_mdot_value_code = parseFloating(
      requireString(
          entries, consumed, "mode.ic_bh_mdot_value_code",
          defaultFor("mode.ic_bh_mdot_value_code")),
      "mode.ic_bh_mdot_value_code");
  if (!has_explicit_ic_convention && frozen.config.mode.ic_file != "generated") {
    throw ConfigError(
        "mode.ic_convention is required for external initial conditions; the runtime will not guess unit, frame, scale-factor, velocity, or species conventions from mode.ic_file");
  }
  frozen.config.mode.zoom_high_res_region = parseBool(
      requireString(entries, consumed, "mode.zoom_high_res_region", defaultFor("mode.zoom_high_res_region")),
      "mode.zoom_high_res_region");
  frozen.config.mode.zoom_region_file =
      requireString(entries, consumed, "mode.zoom_region_file", defaultFor("mode.zoom_region_file"));
  frozen.config.mode.zoom_long_range_strategy = parseZoomLongRangeStrategy(requireString(
      entries,
      consumed,
      "mode.zoom_long_range_strategy",
      defaultFor("mode.zoom_long_range_strategy")));
  frozen.config.mode.zoom_region_center_x_mpc_comoving = parseFloating(
      requireString(entries, consumed, "mode.zoom_region_center_x", defaultFor("mode.zoom_region_center_x")),
      "mode.zoom_region_center_x");
  frozen.config.mode.zoom_region_center_y_mpc_comoving = parseFloating(
      requireString(entries, consumed, "mode.zoom_region_center_y", defaultFor("mode.zoom_region_center_y")),
      "mode.zoom_region_center_y");
  frozen.config.mode.zoom_region_center_z_mpc_comoving = parseFloating(
      requireString(entries, consumed, "mode.zoom_region_center_z", defaultFor("mode.zoom_region_center_z")),
      "mode.zoom_region_center_z");
  frozen.config.mode.zoom_region_radius_mpc_comoving = parseFloating(
      requireString(entries, consumed, "mode.zoom_region_radius", defaultFor("mode.zoom_region_radius")),
      "mode.zoom_region_radius");
  frozen.config.mode.zoom_focused_pm_grid_nx = parseNumber<int>(
      requireString(entries, consumed, "mode.zoom_focused_pm_grid_nx", defaultFor("mode.zoom_focused_pm_grid_nx")),
      "mode.zoom_focused_pm_grid_nx");
  frozen.config.mode.zoom_focused_pm_grid_ny = parseNumber<int>(
      requireString(entries, consumed, "mode.zoom_focused_pm_grid_ny", defaultFor("mode.zoom_focused_pm_grid_ny")),
      "mode.zoom_focused_pm_grid_ny");
  frozen.config.mode.zoom_focused_pm_grid_nz = parseNumber<int>(
      requireString(entries, consumed, "mode.zoom_focused_pm_grid_nz", defaultFor("mode.zoom_focused_pm_grid_nz")),
      "mode.zoom_focused_pm_grid_nz");
  frozen.config.mode.zoom_contamination_radius_mpc_comoving = parseFloating(
      requireString(
          entries,
          consumed,
          "mode.zoom_contamination_radius",
          defaultFor("mode.zoom_contamination_radius")),
      "mode.zoom_contamination_radius");
  frozen.config.mode.hydro_boundary = parseModeHydroBoundary(
      requireString(entries, consumed, "mode.hydro_boundary", defaultFor("mode.hydro_boundary")));
  frozen.config.mode.gravity_boundary = parseModeGravityBoundary(
      requireString(entries, consumed, "mode.gravity_boundary", defaultFor("mode.gravity_boundary")));

  frozen.config.cosmology.omega_matter = parseFloating(
      requireString(entries, consumed, "cosmology.omega_matter", "0.315"),
      "cosmology.omega_matter");
  frozen.config.cosmology.omega_lambda = parseFloating(
      requireString(entries, consumed, "cosmology.omega_lambda", "0.685"),
      "cosmology.omega_lambda");
  frozen.config.cosmology.omega_baryon = parseFloating(
      requireString(entries, consumed, "cosmology.omega_baryon", "0.049"),
      "cosmology.omega_baryon");
  frozen.config.cosmology.hubble_param = parseFloating(
      requireString(entries, consumed, "cosmology.hubble_param", "0.674"),
      "cosmology.hubble_param");
  frozen.config.cosmology.sigma8 =
      parseFloating(requireString(entries, consumed, "cosmology.sigma8", "0.811"), "cosmology.sigma8");
  frozen.config.cosmology.scalar_index_ns = parseFloating(
      requireString(entries, consumed, "cosmology.scalar_index_ns", "0.965"),
      "cosmology.scalar_index_ns");
  const bool has_box_size_x = entries.contains("cosmology.box_size_x");
  const bool has_box_size_y = entries.contains("cosmology.box_size_y");
  const bool has_box_size_z = entries.contains("cosmology.box_size_z");
  const bool has_box_size_scalar = entries.contains("cosmology.box_size");
  if ((has_box_size_x || has_box_size_y || has_box_size_z) &&
      !(has_box_size_x && has_box_size_y && has_box_size_z)) {
    throw ConfigError(
        "cosmology.box_size_x, cosmology.box_size_y, and cosmology.box_size_z must be specified together");
  }
  if (has_box_size_x) {
    frozen.config.cosmology.box_size_x_mpc_comoving = parseLengthMpc(
        requireString(entries, consumed, "cosmology.box_size_x", defaultFor("cosmology.box_size_x")),
        frozen.config.units.length_unit,
        "cosmology.box_size_x");
    frozen.config.cosmology.box_size_y_mpc_comoving = parseLengthMpc(
        requireString(entries, consumed, "cosmology.box_size_y", defaultFor("cosmology.box_size_y")),
        frozen.config.units.length_unit,
        "cosmology.box_size_y");
    frozen.config.cosmology.box_size_z_mpc_comoving = parseLengthMpc(
        requireString(entries, consumed, "cosmology.box_size_z", defaultFor("cosmology.box_size_z")),
        frozen.config.units.length_unit,
        "cosmology.box_size_z");
    if (has_box_size_scalar) {
      consumed.insert("cosmology.box_size");
    }
  } else {
    const double box_size_scalar_mpc = parseLengthMpc(
        requireString(entries, consumed, "cosmology.box_size", defaultFor("cosmology.box_size")),
        frozen.config.units.length_unit,
        "cosmology.box_size");
    frozen.config.cosmology.box_size_x_mpc_comoving = box_size_scalar_mpc;
    frozen.config.cosmology.box_size_y_mpc_comoving = box_size_scalar_mpc;
    frozen.config.cosmology.box_size_z_mpc_comoving = box_size_scalar_mpc;
  }
  frozen.config.cosmology.box_size_mpc_comoving = frozen.config.cosmology.box_size_x_mpc_comoving;

  normalizeScaleRedshiftPair(
      entries,
      consumed,
      "numerics.a_begin",
      "numerics.z_begin",
      parseFloating(defaultFor("numerics.a_begin"), "numerics.a_begin"),
      parseFloating(defaultFor("numerics.z_begin"), "numerics.z_begin"),
      frozen.config.numerics.a_begin,
      frozen.config.numerics.z_begin);
  normalizeScaleRedshiftPair(
      entries,
      consumed,
      "numerics.a_end",
      "numerics.z_end",
      parseFloating(defaultFor("numerics.a_end"), "numerics.a_end"),
      parseFloating(defaultFor("numerics.z_end"), "numerics.z_end"),
      frozen.config.numerics.a_end,
      frozen.config.numerics.z_end);
  frozen.config.numerics.t_code_begin = parseFloating(
      requireString(entries, consumed, "numerics.t_code_begin", defaultFor("numerics.t_code_begin")),
      "numerics.t_code_begin");
  frozen.config.numerics.t_code_end = parseFloating(
      requireString(entries, consumed, "numerics.t_code_end", defaultFor("numerics.t_code_end")),
      "numerics.t_code_end");
  frozen.config.numerics.t_phys_begin = parseFloating(
      requireString(entries, consumed, "numerics.t_phys_begin", defaultFor("numerics.t_phys_begin")),
      "numerics.t_phys_begin");
  frozen.config.numerics.t_phys_end = parseFloating(
      requireString(entries, consumed, "numerics.t_phys_end", defaultFor("numerics.t_phys_end")),
      "numerics.t_phys_end");
  frozen.config.numerics.integrator_time_variable = parseIntegratorTimeVariable(requireString(
      entries,
      consumed,
      "numerics.integrator_time_variable",
      defaultFor("numerics.integrator_time_variable")));
  frozen.config.numerics.cosmology_max_delta_ln_a = parseFloating(
      requireString(entries, consumed, "numerics.cosmology_max_delta_ln_a", defaultFor("numerics.cosmology_max_delta_ln_a")),
      "numerics.cosmology_max_delta_ln_a");
  frozen.config.numerics.cosmology_max_hubble_time_fraction = parseFloating(
      requireString(entries, consumed, "numerics.cosmology_max_hubble_time_fraction", defaultFor("numerics.cosmology_max_hubble_time_fraction")),
      "numerics.cosmology_max_hubble_time_fraction");
  frozen.config.numerics.source_max_fractional_change = parseFloating(
      requireString(entries, consumed, "numerics.source_max_fractional_change", defaultFor("numerics.source_max_fractional_change")),
      "numerics.source_max_fractional_change");
  frozen.config.numerics.hydro_density_floor_code = parseFloating(
      requireString(entries, consumed, "numerics.hydro_density_floor_code", defaultFor("numerics.hydro_density_floor_code")),
      "numerics.hydro_density_floor_code");
  frozen.config.numerics.hydro_pressure_floor_code = parseFloating(
      requireString(entries, consumed, "numerics.hydro_pressure_floor_code", defaultFor("numerics.hydro_pressure_floor_code")),
      "numerics.hydro_pressure_floor_code");
  frozen.config.numerics.hydro_active_batch_max_cells = parseNumber<std::uint64_t>(
      requireString(
          entries,
          consumed,
          "numerics.hydro_active_batch_max_cells",
          defaultFor("numerics.hydro_active_batch_max_cells")),
      "numerics.hydro_active_batch_max_cells");
  frozen.config.numerics.hydro_reconstruction = parseHydroReconstructionMode(requireString(
      entries, consumed, "numerics.hydro_reconstruction", defaultFor("numerics.hydro_reconstruction")));
  frozen.config.numerics.hydro_slope_limiter = parseHydroSlopeLimiterMode(requireString(
      entries, consumed, "numerics.hydro_slope_limiter", defaultFor("numerics.hydro_slope_limiter")));
  frozen.config.numerics.hydro_integration_scheme = parseHydroIntegrationMode(requireString(
      entries, consumed, "numerics.hydro_integration_scheme",
      defaultFor("numerics.hydro_integration_scheme")));
  frozen.config.numerics.hydro_riemann_solver = parseHydroRiemannSolverMode(requireString(
      entries, consumed, "numerics.hydro_riemann_solver", defaultFor("numerics.hydro_riemann_solver")));
  frozen.config.numerics.hydro_dual_energy = parseHydroDualEnergyConfigMode(requireString(
      entries, consumed, "numerics.hydro_dual_energy", defaultFor("numerics.hydro_dual_energy")));
  frozen.config.numerics.hydro_dual_energy_eta = parseFloating(
      requireString(entries, consumed, "numerics.hydro_dual_energy_eta",
                    defaultFor("numerics.hydro_dual_energy_eta")),
      "numerics.hydro_dual_energy_eta");
  frozen.config.numerics.max_global_steps = parseNumber<int>(
      requireString(entries, consumed, "numerics.max_global_steps", "1024"),
      "numerics.max_global_steps");
  frozen.config.numerics.hierarchical_max_rung = parseNumber<int>(
      requireString(entries, consumed, "numerics.hierarchical_max_rung", "0"),
      "numerics.hierarchical_max_rung");
  frozen.config.numerics.hierarchical_coarse_displacement_criterion =
      parseHierarchicalCoarseDisplacementCriterion(requireString(
          entries, consumed, "numerics.hierarchical_coarse_displacement_criterion",
          defaultFor("numerics.hierarchical_coarse_displacement_criterion")));
  frozen.config.numerics.hierarchical_coarse_rms_displacement_factor = parseFloating(
      requireString(entries, consumed, "numerics.hierarchical_coarse_rms_displacement_factor",
                    defaultFor("numerics.hierarchical_coarse_rms_displacement_factor")),
      "numerics.hierarchical_coarse_rms_displacement_factor");
  frozen.config.numerics.hierarchical_encounter_timestep_enabled = parseBool(
      requireString(entries, consumed, "numerics.hierarchical_encounter_timestep_enabled",
                    defaultFor("numerics.hierarchical_encounter_timestep_enabled")),
      "numerics.hierarchical_encounter_timestep_enabled");
  frozen.config.numerics.gravity_timestep_eta = parseFloating(
      requireString(entries, consumed, "numerics.gravity_timestep_eta",
                    defaultFor("numerics.gravity_timestep_eta")),
      "numerics.gravity_timestep_eta");
  frozen.config.numerics.amr_max_level = parseNumber<int>(
      requireString(entries, consumed, "numerics.amr_max_level", "10"), "numerics.amr_max_level");
  frozen.config.numerics.amr_regrid_mode = parseAmrRegridConfigMode(requireString(
      entries, consumed, "numerics.amr_regrid_mode", defaultFor("numerics.amr_regrid_mode")));
  frozen.config.numerics.amr_predictive_feedback_age_years = parseFloating(
      requireString(entries, consumed, "numerics.amr_predictive_feedback_age_years",
                    defaultFor("numerics.amr_predictive_feedback_age_years")), "numerics.amr_predictive_feedback_age_years");
  frozen.config.numerics.amr_predictive_event_lookahead_years = parseFloating(
      requireString(entries, consumed, "numerics.amr_predictive_event_lookahead_years",
                    defaultFor("numerics.amr_predictive_event_lookahead_years")),
      "numerics.amr_predictive_event_lookahead_years");
  frozen.config.numerics.amr_time_stepping = parseAmrTimeSteppingConfigMode(requireString(
      entries, consumed, "numerics.amr_time_stepping", defaultFor("numerics.amr_time_stepping")));
  frozen.config.numerics.global_step_max_retries = parseNumber<int>(
      requireString(entries, consumed, "numerics.global_step_max_retries",
                    defaultFor("numerics.global_step_max_retries")), "numerics.global_step_max_retries");
  frozen.config.numerics.amr_regrid_interval_steps = parseNumber<int>(
      requireString(entries, consumed, "numerics.amr_regrid_interval_steps",
                    defaultFor("numerics.amr_regrid_interval_steps")), "numerics.amr_regrid_interval_steps");
  frozen.config.numerics.amr_regrid_max_operations = parseNumber<int>(
      requireString(entries, consumed, "numerics.amr_regrid_max_operations",
                    defaultFor("numerics.amr_regrid_max_operations")), "numerics.amr_regrid_max_operations");
  frozen.config.numerics.amr_refine_cell_mass_code = parseFloating(
      requireString(entries, consumed, "numerics.amr_refine_cell_mass_code",
                    defaultFor("numerics.amr_refine_cell_mass_code")), "numerics.amr_refine_cell_mass_code");
  frozen.config.numerics.amr_refine_jump_threshold = parseFloating(
      requireString(entries, consumed, "numerics.amr_refine_jump_threshold",
                    defaultFor("numerics.amr_refine_jump_threshold")), "numerics.amr_refine_jump_threshold");
  frozen.config.numerics.amr_refine_jeans_cells = parseFloating(
      requireString(entries, consumed, "numerics.amr_refine_jeans_cells",
                    defaultFor("numerics.amr_refine_jeans_cells")), "numerics.amr_refine_jeans_cells");
  frozen.config.numerics.amr_derefine_hysteresis = parseFloating(
      requireString(entries, consumed, "numerics.amr_derefine_hysteresis",
                    defaultFor("numerics.amr_derefine_hysteresis")), "numerics.amr_derefine_hysteresis");
  for (auto [key, target] : std::initializer_list<std::pair<const char*, double*>>{
           {"numerics.amr_refine_cooling_length_cells", &frozen.config.numerics.amr_refine_cooling_length_cells},
           {"numerics.amr_refine_metallicity_jump", &frozen.config.numerics.amr_refine_metallicity_jump},
           {"numerics.amr_refine_species_jump", &frozen.config.numerics.amr_refine_species_jump},
           {"numerics.amr_refine_vorticity", &frozen.config.numerics.amr_refine_vorticity}}) {
    *target = parseFloating(requireString(entries, consumed, key, defaultFor(key)), key);
  }
  frozen.config.numerics.gravity_softening_kpc_comoving = parseLengthKpc(
      requireString(entries, consumed, "numerics.gravity_softening", "1.0 kpc"),
      frozen.config.units.length_unit,
      "numerics.gravity_softening");
  auto parseOptionalSofteningKpc = [&](const char* key) -> double {
    const std::string value = requireString(entries, consumed, key, defaultFor(key));
    if (trim(value).empty()) {
      return -1.0;
    }
    return parseLengthKpc(value, frozen.config.units.length_unit, key);
  };
  frozen.config.numerics.gravity_softening_gas_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_gas");
  frozen.config.numerics.gravity_softening_dark_matter_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_dark_matter");
  frozen.config.numerics.gravity_softening_star_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_star");
  frozen.config.numerics.gravity_softening_black_hole_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_black_hole");
  frozen.config.numerics.gravity_softening_tracer_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_tracer");
  frozen.config.numerics.gravity_solver = parseGravitySolver(
      requireString(entries, consumed, "numerics.gravity_solver", defaultFor("numerics.gravity_solver")));
  frozen.config.numerics.hydro_solver = parseHydroSolver(
      requireString(entries, consumed, "numerics.hydro_solver", defaultFor("numerics.hydro_solver")));
  const bool has_pm_grid_nx = entries.contains("numerics.treepm_pm_grid_nx");
  const bool has_pm_grid_ny = entries.contains("numerics.treepm_pm_grid_ny");
  const bool has_pm_grid_nz = entries.contains("numerics.treepm_pm_grid_nz");
  const bool has_pm_grid_scalar = entries.contains("numerics.treepm_pm_grid");
  if ((has_pm_grid_nx || has_pm_grid_ny || has_pm_grid_nz) &&
      !(has_pm_grid_nx && has_pm_grid_ny && has_pm_grid_nz)) {
    throw ConfigError(
        "numerics.treepm_pm_grid_nx, numerics.treepm_pm_grid_ny, and numerics.treepm_pm_grid_nz must be specified together");
  }
  if (has_pm_grid_nx) {
    frozen.config.numerics.treepm_pm_grid_nx = parseNumber<int>(
        requireString(entries, consumed, "numerics.treepm_pm_grid_nx", defaultFor("numerics.treepm_pm_grid_nx")),
        "numerics.treepm_pm_grid_nx");
    frozen.config.numerics.treepm_pm_grid_ny = parseNumber<int>(
        requireString(entries, consumed, "numerics.treepm_pm_grid_ny", defaultFor("numerics.treepm_pm_grid_ny")),
        "numerics.treepm_pm_grid_ny");
    frozen.config.numerics.treepm_pm_grid_nz = parseNumber<int>(
        requireString(entries, consumed, "numerics.treepm_pm_grid_nz", defaultFor("numerics.treepm_pm_grid_nz")),
        "numerics.treepm_pm_grid_nz");
    if (has_pm_grid_scalar) {
      consumed.insert("numerics.treepm_pm_grid");
    }
  } else {
    const int pm_grid_scalar = parseNumber<int>(
        requireString(entries, consumed, "numerics.treepm_pm_grid", defaultFor("numerics.treepm_pm_grid")),
        "numerics.treepm_pm_grid");
    frozen.config.numerics.treepm_pm_grid_nx = pm_grid_scalar;
    frozen.config.numerics.treepm_pm_grid_ny = pm_grid_scalar;
    frozen.config.numerics.treepm_pm_grid_nz = pm_grid_scalar;
  }
  frozen.config.numerics.treepm_pm_grid = frozen.config.numerics.treepm_pm_grid_nx;
  frozen.config.numerics.treepm_asmth_cells = parseFloating(
      requireString(entries, consumed, "numerics.treepm_asmth_cells", defaultFor("numerics.treepm_asmth_cells")),
      "numerics.treepm_asmth_cells");
  frozen.config.numerics.treepm_rcut_cells = parseFloating(
      requireString(entries, consumed, "numerics.treepm_rcut_cells", defaultFor("numerics.treepm_rcut_cells")),
      "numerics.treepm_rcut_cells");
  frozen.config.numerics.treepm_tree_opening_criterion = parseTreePmOpeningCriterion(requireString(
      entries,
      consumed,
      "numerics.treepm_tree_opening_criterion",
      defaultFor("numerics.treepm_tree_opening_criterion")));
  frozen.config.numerics.treepm_tree_opening_theta = parseFloating(
      requireString(
          entries,
          consumed,
          "numerics.treepm_tree_opening_theta",
          defaultFor("numerics.treepm_tree_opening_theta")),
      "numerics.treepm_tree_opening_theta");
  frozen.config.numerics.treepm_tree_relative_force_tolerance = parseFloating(
      requireString(
          entries,
          consumed,
          "numerics.treepm_tree_relative_force_tolerance",
          defaultFor("numerics.treepm_tree_relative_force_tolerance")),
      "numerics.treepm_tree_relative_force_tolerance");
  frozen.config.numerics.treepm_tree_relative_force_acceleration_floor = parseFloating(
      requireString(
          entries,
          consumed,
          "numerics.treepm_tree_relative_force_acceleration_floor",
          defaultFor("numerics.treepm_tree_relative_force_acceleration_floor")),
      "numerics.treepm_tree_relative_force_acceleration_floor");
  frozen.config.numerics.treepm_gaussian_pair_lookup_enabled = parseBool(
      requireString(entries, consumed, "numerics.treepm_gaussian_pair_lookup_enabled", defaultFor("numerics.treepm_gaussian_pair_lookup_enabled")),
      "numerics.treepm_gaussian_pair_lookup_enabled");
  frozen.config.numerics.treepm_full_residual_lookup_enabled = parseBool(
      requireString(entries, consumed, "numerics.treepm_full_residual_lookup_enabled", defaultFor("numerics.treepm_full_residual_lookup_enabled")),
      "numerics.treepm_full_residual_lookup_enabled");
  frozen.config.numerics.treepm_radial_coefficients_enabled = parseBool(
      requireString(entries, consumed, "numerics.treepm_radial_coefficients_enabled", defaultFor("numerics.treepm_radial_coefficients_enabled")),
      "numerics.treepm_radial_coefficients_enabled");
  frozen.config.numerics.treepm_leaf_simd_enabled = parseBool(
      requireString(entries, consumed, "numerics.treepm_leaf_simd_enabled", defaultFor("numerics.treepm_leaf_simd_enabled")),
      "numerics.treepm_leaf_simd_enabled");
  frozen.config.numerics.treepm_residual_block_size = parseNumber<int>(
      requireString(entries, consumed, "numerics.treepm_residual_block_size", defaultFor("numerics.treepm_residual_block_size")),
      "numerics.treepm_residual_block_size");
  frozen.config.numerics.treepm_tree_max_leaf_size = parseNumber<int>(
      requireString(entries, consumed, "numerics.treepm_tree_max_leaf_size", defaultFor("numerics.treepm_tree_max_leaf_size")),
      "numerics.treepm_tree_max_leaf_size");
  frozen.config.numerics.treepm_full_mac_diagnostics = parseBool(
      requireString(entries, consumed, "numerics.treepm_full_mac_diagnostics", defaultFor("numerics.treepm_full_mac_diagnostics")),
      "numerics.treepm_full_mac_diagnostics");
  frozen.config.numerics.treepm_adaptive_acceptance_enabled = parseBool(
      requireString(entries, consumed, "numerics.treepm_adaptive_acceptance_enabled",
                    defaultFor("numerics.treepm_adaptive_acceptance_enabled")),
      "numerics.treepm_adaptive_acceptance_enabled");
  frozen.config.numerics.treepm_adaptive_maximum_opening_angle = parseFloating(
      requireString(entries, consumed, "numerics.treepm_adaptive_maximum_opening_angle",
                    defaultFor("numerics.treepm_adaptive_maximum_opening_angle")),
      "numerics.treepm_adaptive_maximum_opening_angle");
  frozen.config.numerics.treepm_adaptive_error_estimator = parseTreePmAdaptiveErrorEstimator(
      requireString(entries, consumed, "numerics.treepm_adaptive_error_estimator",
                    defaultFor("numerics.treepm_adaptive_error_estimator")));
  frozen.config.numerics.treepm_cutoff_node_policy = parseTreePmCutoffNodePolicy(
      requireString(entries, consumed, "numerics.treepm_cutoff_node_policy",
                    defaultFor("numerics.treepm_cutoff_node_policy")));
  frozen.config.numerics.treepm_group_walk_enabled = parseBool(
      requireString(entries, consumed, "numerics.treepm_group_walk_enabled",
                    defaultFor("numerics.treepm_group_walk_enabled")),
      "numerics.treepm_group_walk_enabled");
  frozen.config.numerics.treepm_group_walk_target_size = parseNumber<int>(
      requireString(entries, consumed, "numerics.treepm_group_walk_target_size",
                    defaultFor("numerics.treepm_group_walk_target_size")),
      "numerics.treepm_group_walk_target_size");
  frozen.config.numerics.treepm_identical_source_tree_reuse_enabled = parseBool(
      requireString(entries, consumed, "numerics.treepm_identical_source_tree_reuse_enabled",
                    defaultFor("numerics.treepm_identical_source_tree_reuse_enabled")),
      "numerics.treepm_identical_source_tree_reuse_enabled");
  frozen.config.numerics.treepm_topology_refit_enabled = parseBool(
      requireString(entries, consumed, "numerics.treepm_topology_refit_enabled",
                    defaultFor("numerics.treepm_topology_refit_enabled")),
      "numerics.treepm_topology_refit_enabled");
  frozen.config.numerics.treepm_assignment_scheme = parseTreePmAssignmentScheme(requireString(
      entries,
      consumed,
      "numerics.treepm_assignment_scheme",
      defaultFor("numerics.treepm_assignment_scheme")));
  frozen.config.numerics.treepm_pm_force_gradient = parseTreePmForceGradient(requireString(
      entries, consumed, "numerics.treepm_pm_force_gradient",
      defaultFor("numerics.treepm_pm_force_gradient")));
  frozen.config.numerics.treepm_enable_window_deconvolution = parseBool(
      requireString(
          entries,
          consumed,
          "numerics.treepm_enable_window_deconvolution",
          defaultFor("numerics.treepm_enable_window_deconvolution")),
      "numerics.treepm_enable_window_deconvolution");
  frozen.config.numerics.treepm_update_cadence_steps = parseNumber<int>(
      requireString(
          entries,
          consumed,
          "numerics.treepm_update_cadence_steps",
          defaultFor("numerics.treepm_update_cadence_steps")),
      "numerics.treepm_update_cadence_steps");
  frozen.config.numerics.treepm_pm_reuse_policy = parseTreePmPmReusePolicy(requireString(
      entries, consumed, "numerics.treepm_pm_reuse_policy",
      defaultFor("numerics.treepm_pm_reuse_policy")));
  frozen.config.numerics.treepm_multirate_rollback_memory_policy = parseMultirateRollbackMemoryPolicy(requireString(
      entries, consumed, "numerics.treepm_multirate_rollback_memory_policy",
      defaultFor("numerics.treepm_multirate_rollback_memory_policy")));
  frozen.config.numerics.treepm_allow_diagnostic_naive_dft = parseBool(
      requireString(
          entries,
          consumed,
          "numerics.treepm_allow_diagnostic_naive_dft",
          defaultFor("numerics.treepm_allow_diagnostic_naive_dft")),
      "numerics.treepm_allow_diagnostic_naive_dft");
  frozen.config.numerics.treepm_pm_decomposition_mode = parsePmDecompositionMode(requireString(
      entries,
      consumed,
      "numerics.treepm_pm_decomposition_mode",
      defaultFor("numerics.treepm_pm_decomposition_mode")));
  frozen.config.numerics.treepm_tree_exchange_batch_bytes = parseNumber<std::uint64_t>(
      requireString(
          entries,
          consumed,
          "numerics.treepm_tree_exchange_batch_bytes",
          defaultFor("numerics.treepm_tree_exchange_batch_bytes")),
      "numerics.treepm_tree_exchange_batch_bytes");
  frozen.config.numerics.treepm_pm_exchange_batch_bytes = parseNumber<std::uint64_t>(
      requireString(
          entries,
          consumed,
          "numerics.treepm_pm_exchange_batch_bytes",
          defaultFor("numerics.treepm_pm_exchange_batch_bytes")),
      "numerics.treepm_pm_exchange_batch_bytes");

  frozen.config.physics.enable_cooling = parseBool(
      requireString(entries, consumed, "physics.enable_cooling", "true"), "physics.enable_cooling");
  frozen.config.physics.enable_star_formation =
      parseBool(requireString(entries, consumed, "physics.enable_star_formation", "true"),
                "physics.enable_star_formation");
  frozen.config.physics.enable_feedback = parseBool(
      requireString(entries, consumed, "physics.enable_feedback", "true"), "physics.enable_feedback");
  frozen.config.physics.enable_stellar_evolution = parseBool(
      requireString(entries, consumed, "physics.enable_stellar_evolution", "true"),
      "physics.enable_stellar_evolution");
  frozen.config.physics.reionization_model =
      requireString(entries, consumed, "physics.reionization_model", "hm12");
  frozen.config.physics.uv_background_model = parseUvBackgroundModel(
      requireString(entries, consumed, "physics.uv_background_model", defaultFor("physics.uv_background_model")));
  frozen.config.physics.self_shielding_model = parseSelfShieldingModel(
      requireString(entries, consumed, "physics.self_shielding_model", defaultFor("physics.self_shielding_model")));
  frozen.config.physics.cooling_model = parseCoolingModel(
      requireString(entries, consumed, "physics.cooling_model", defaultFor("physics.cooling_model")));
  frozen.config.physics.metal_line_table_path =
      requireString(entries, consumed, "physics.metal_line_table_path", "");
  frozen.config.physics.temperature_floor_k = parseFloating(
      requireString(entries, consumed, "physics.temperature_floor_k", "100.0"),
      "physics.temperature_floor_k");
  frozen.config.physics.star_formation_model = parseStarFormationModelKind(
      requireString(entries, consumed, "physics.star_formation_model",
                    defaultFor("physics.star_formation_model")));
  frozen.config.physics.sf_density_threshold_code = parseFloating(
      requireString(entries, consumed, "physics.sf_density_threshold_code", "10.0"),
      "physics.sf_density_threshold_code");
  frozen.config.physics.sf_temperature_threshold_k = parseFloating(
      requireString(entries, consumed, "physics.sf_temperature_threshold_k", "1.0e4"),
      "physics.sf_temperature_threshold_k");
  frozen.config.physics.sf_min_converging_flow_rate_code = parseFloating(
      requireString(entries, consumed, "physics.sf_min_converging_flow_rate_code", "0.0"),
      "physics.sf_min_converging_flow_rate_code");
  frozen.config.physics.sf_epsilon_ff = parseFloating(
      requireString(entries, consumed, "physics.sf_epsilon_ff", "0.01"),
      "physics.sf_epsilon_ff");
  frozen.config.physics.sf_bound_alpha_vir_max = parseFloating(
      requireString(entries, consumed, "physics.sf_bound_alpha_vir_max", "1.0"),
      "physics.sf_bound_alpha_vir_max");
  frozen.config.physics.sf_require_converging_flow = parseBool(
      requireString(entries, consumed, "physics.sf_require_converging_flow", "true"),
      "physics.sf_require_converging_flow");
  frozen.config.physics.sf_collapse_timescale = parseStarFormationCollapseTimescale(
      requireString(entries, consumed, "physics.sf_collapse_timescale", "free_fall"));
  frozen.config.physics.sf_jeans_mass_floor_code = parseMassCode(
      requireString(entries, consumed, "physics.sf_jeans_mass_floor_code", "100.0"),
      frozen.config.units.mass_unit,
      "physics.sf_jeans_mass_floor_code");
  frozen.config.physics.sf_star_particle_mass_policy = parseStarParticleMassPolicy(
      requireString(entries, consumed, "physics.sf_star_particle_mass_policy", "fixed"));
  frozen.config.physics.sf_target_star_particle_mass_code = parseMassCode(
      requireString(entries, consumed, "physics.sf_target_star_particle_mass_code", "0.1"),
      frozen.config.units.mass_unit,
      "physics.sf_target_star_particle_mass_code");
  frozen.config.physics.sf_target_star_particle_mass_fraction = parseFloating(
      requireString(entries, consumed, "physics.sf_target_star_particle_mass_fraction", "0.25"),
      "physics.sf_target_star_particle_mass_fraction");
  frozen.config.physics.sf_min_star_particle_mass_code = parseMassCode(
      requireString(entries, consumed, "physics.sf_min_star_particle_mass_code", "0.1"),
      frozen.config.units.mass_unit,
      "physics.sf_min_star_particle_mass_code");
  frozen.config.physics.sf_max_star_particle_mass_code = parseMassCode(
      requireString(entries, consumed, "physics.sf_max_star_particle_mass_code", "1.0e30"),
      frozen.config.units.mass_unit,
      "physics.sf_max_star_particle_mass_code");
  frozen.config.physics.sf_max_spawn_particles_per_cell_step =
      parseNumber<std::uint32_t>(
          requireString(entries, consumed, "physics.sf_max_spawn_particles_per_cell_step", "8"),
          "physics.sf_max_spawn_particles_per_cell_step");
  frozen.config.physics.sf_max_fractional_mass_conversion = parseFloating(
      requireString(entries, consumed, "physics.sf_max_fractional_mass_conversion", "0.25"),
      "physics.sf_max_fractional_mass_conversion");
  frozen.config.physics.sf_min_remaining_gas_fraction = parseFloating(
      requireString(entries, consumed, "physics.sf_min_remaining_gas_fraction", "0.01"),
      "physics.sf_min_remaining_gas_fraction");
  frozen.config.physics.sf_min_remaining_gas_mass_code = parseMassCode(
      requireString(entries, consumed, "physics.sf_min_remaining_gas_mass_code", "1.0e-12"),
      frozen.config.units.mass_unit,
      "physics.sf_min_remaining_gas_mass_code");
  frozen.config.physics.sf_temperature_safety_ceiling_k = parseFloating(
      requireString(entries, consumed, "physics.sf_temperature_safety_ceiling_k", "0.0"),
      "physics.sf_temperature_safety_ceiling_k");
  frozen.config.physics.sf_stochastic_spawning = parseBool(
      requireString(entries, consumed, "physics.sf_stochastic_spawning", "true"),
      "physics.sf_stochastic_spawning");
  frozen.config.physics.sf_random_seed = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "physics.sf_random_seed", "123456789"),
      "physics.sf_random_seed");
  frozen.config.physics.sf_effective_parameter_set = requireString(
      entries, consumed, "physics.sf_effective_parameter_set", "chui_sh03_tng_like_v1");
  frozen.config.physics.sf_effective_n_h_threshold_cgs = parseNumberDensityCgs(
      requireString(entries, consumed, "physics.sf_effective_n_h_threshold", "0.13 cm^-3"),
      "physics.sf_effective_n_h_threshold");
  frozen.config.physics.sf_effective_min_baryon_overdensity = parseFloating(
      requireString(entries, consumed, "physics.sf_effective_min_baryon_overdensity", "0.0"),
      "physics.sf_effective_min_baryon_overdensity");
  frozen.config.physics.sf_effective_t_star_at_threshold_code = parseTimeCode(
      requireString(entries, consumed, "physics.sf_effective_t_star_at_threshold", "2.2 Gyr"),
      frozen.config.units, "physics.sf_effective_t_star_at_threshold");
  frozen.config.physics.sf_effective_evaporation_factor_at_threshold = parseFloating(
      requireString(entries, consumed, "physics.sf_effective_evaporation_factor_at_threshold", "1000.0"),
      "physics.sf_effective_evaporation_factor_at_threshold");
  frozen.config.physics.sf_effective_evaporation_density_exponent = parseFloating(
      requireString(entries, consumed, "physics.sf_effective_evaporation_density_exponent", "-0.8"),
      "physics.sf_effective_evaporation_density_exponent");
  frozen.config.physics.sf_effective_cold_phase_temperature_k = parseTemperatureKelvin(
      requireString(entries, consumed, "physics.sf_effective_cold_phase_temperature", "1000 K"),
      "physics.sf_effective_cold_phase_temperature");
  frozen.config.physics.sf_effective_supernova_specific_energy_code = parseSpecificEnergyCode(
      requireString(entries, consumed, "physics.sf_effective_supernova_specific_energy_code", "1.0e5 code"),
      frozen.config.units, "physics.sf_effective_supernova_specific_energy_code");
  frozen.config.physics.sf_effective_massive_star_fraction = parseFloating(
      requireString(entries, consumed, "physics.sf_effective_massive_star_fraction", "0.1"),
      "physics.sf_effective_massive_star_fraction");
  frozen.config.physics.sf_effective_q_eos = parseFloating(
      requireString(entries, consumed, "physics.sf_effective_q_eos", "0.3"),
      "physics.sf_effective_q_eos");
  frozen.config.physics.sf_effective_isothermal_temperature_k = parseTemperatureKelvin(
      requireString(entries, consumed, "physics.sf_effective_isothermal_temperature", "1.0e4 K"),
      "physics.sf_effective_isothermal_temperature");
  frozen.config.physics.sf_effective_hot_excess_tolerance = parseFloating(
      requireString(entries, consumed, "physics.sf_effective_hot_excess_tolerance", "1.0"),
      "physics.sf_effective_hot_excess_tolerance");
  frozen.config.physics.sf_effective_eos_relaxation = parseEffectiveIsmEosRelaxation(
      requireString(entries, consumed, "physics.sf_effective_eos_relaxation", "instantaneous"));
  frozen.config.physics.sf_effective_eos_relaxation_timescale_code = parseTimeCode(
      requireString(entries, consumed, "physics.sf_effective_eos_relaxation_timescale", "0.0 code"),
      frozen.config.units, "physics.sf_effective_eos_relaxation_timescale");
  frozen.config.physics.sf_effective_eos_table_bins = parseNumber<std::uint32_t>(
      requireString(entries, consumed, "physics.sf_effective_eos_table_bins", "256"),
      "physics.sf_effective_eos_table_bins");
  frozen.config.physics.sf_effective_eos_max_density_ratio = parseFloating(
      requireString(entries, consumed, "physics.sf_effective_eos_max_density_ratio", "1.0e6"),
      "physics.sf_effective_eos_max_density_ratio");
  frozen.config.physics.sf_effective_birth_mass_convention = parseEffectiveIsmBirthMassConvention(
      requireString(entries, consumed, "physics.sf_effective_birth_mass_convention", "initial_ssp_mass"));
  frozen.config.physics.sf_effective_feedback_coupling = parseEffectiveIsmFeedbackCoupling(
      requireString(entries, consumed, "physics.sf_effective_feedback_coupling", "external_feedback_calibrated"));
  frozen.config.physics.fb_mode = parseFeedbackMode(
      requireString(entries, consumed, "physics.fb_mode", defaultFor("physics.fb_mode")));
  frozen.config.physics.fb_variant = parseFeedbackVariant(
      requireString(entries, consumed, "physics.fb_variant", defaultFor("physics.fb_variant")));
  frozen.config.physics.fb_use_returned_mass_budget = parseBool(
      requireString(entries, consumed, "physics.fb_use_returned_mass_budget", "true"),
      "physics.fb_use_returned_mass_budget");
  frozen.config.physics.fb_epsilon_thermal = parseFloating(
      requireString(entries, consumed, "physics.fb_epsilon_thermal", "0.6"),
      "physics.fb_epsilon_thermal");
  frozen.config.physics.fb_epsilon_kinetic = parseFloating(
      requireString(entries, consumed, "physics.fb_epsilon_kinetic", "0.3"),
      "physics.fb_epsilon_kinetic");
  frozen.config.physics.fb_epsilon_momentum = parseFloating(
      requireString(entries, consumed, "physics.fb_epsilon_momentum", "0.1"),
      "physics.fb_epsilon_momentum");
  frozen.config.physics.fb_sn_energy_erg_per_mass_code = parseFloating(
      requireString(entries, consumed, "physics.fb_sn_energy_erg_per_mass_code", "1.0e49"),
      "physics.fb_sn_energy_erg_per_mass_code");
  frozen.config.physics.fb_momentum_code_per_mass_code = parseFloating(
      requireString(entries, consumed, "physics.fb_momentum_code_per_mass_code", "3.0e3"),
      "physics.fb_momentum_code_per_mass_code");
  frozen.config.physics.fb_neighbor_count = parseNumber<std::uint32_t>(
      requireString(entries, consumed, "physics.fb_neighbor_count", "8"),
      "physics.fb_neighbor_count");
  frozen.config.physics.fb_delayed_cooling_time_code = parseFloating(
      requireString(entries, consumed, "physics.fb_delayed_cooling_time_code", "0.0"),
      "physics.fb_delayed_cooling_time_code");
  frozen.config.physics.fb_stochastic_event_probability = parseFloating(
      requireString(entries, consumed, "physics.fb_stochastic_event_probability", "0.25"),
      "physics.fb_stochastic_event_probability");
  frozen.config.physics.fb_random_seed = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "physics.fb_random_seed", "42424242"),
      "physics.fb_random_seed");
  frozen.config.physics.stellar_evolution_table_path =
      requireString(entries, consumed, "physics.stellar_evolution_table_path", "");
  frozen.config.physics.stellar_evolution_hubble_time_years = parseFloating(
      requireString(entries, consumed, "physics.stellar_evolution_hubble_time_years", "1.44e10"),
      "physics.stellar_evolution_hubble_time_years");
  frozen.config.physics.stellar_evolution_require_production_table = parseBool(
      requireString(entries, consumed, "physics.stellar_evolution_require_production_table", "false"),
      "physics.stellar_evolution_require_production_table");
  frozen.config.physics.metal_species_mode = parseMetalSpeciesMode(
      requireString(entries, consumed, "physics.metal_species_mode", "total_only"));
  frozen.config.physics.enable_metal_diffusion = parseBool(
      requireString(entries, consumed, "physics.enable_metal_diffusion", "false"),
      "physics.enable_metal_diffusion");
  frozen.config.physics.metal_diffusion_model = parseMetalDiffusionModel(
      requireString(entries, consumed, "physics.metal_diffusion_model", "none"));
  frozen.config.physics.metal_diffusion_time_integrator = parseMetalDiffusionTimeIntegrator(
      requireString(entries, consumed, "physics.metal_diffusion_time_integrator", "explicit_subcycling"));
  frozen.config.physics.metal_diffusion_coefficient = parseFloating(
      requireString(entries, consumed, "physics.metal_diffusion_coefficient", "0.05"),
      "physics.metal_diffusion_coefficient");
  frozen.config.physics.metal_diffusion_cfl = parseFloating(
      requireString(entries, consumed, "physics.metal_diffusion_cfl", "0.4"),
      "physics.metal_diffusion_cfl");
  frozen.config.physics.metal_diffusion_max_subcycles = parseNumber<std::uint32_t>(
      requireString(entries, consumed, "physics.metal_diffusion_max_subcycles", "128"),
      "physics.metal_diffusion_max_subcycles");
  frozen.config.physics.metal_diffusion_max_rkl_stages = parseNumber<std::uint32_t>(
      requireString(entries, consumed, "physics.metal_diffusion_max_rkl_stages", "64"),
      "physics.metal_diffusion_max_rkl_stages");
  frozen.config.physics.metal_diffusion_coefficient_floor_code = parseFloating(
      requireString(entries, consumed, "physics.metal_diffusion_coefficient_floor_code", "0.0"),
      "physics.metal_diffusion_coefficient_floor_code");
  frozen.config.physics.metal_diffusion_coefficient_ceiling_code = parseFloating(
      requireString(entries, consumed, "physics.metal_diffusion_coefficient_ceiling_code", "1.0e30"),
      "physics.metal_diffusion_coefficient_ceiling_code");
  frozen.config.physics.enable_black_hole_agn = parseBool(
      requireString(entries, consumed, "physics.enable_black_hole_agn", "false"),
      "physics.enable_black_hole_agn");
  frozen.config.physics.bh_enable_seeding = parseBool(
      requireString(entries, consumed, "physics.bh_enable_seeding", "false"),
      "physics.bh_enable_seeding");
  frozen.config.physics.bh_seed_halo_mass_threshold_code = parseFloating(
      requireString(entries, consumed, "physics.bh_seed_halo_mass_threshold_code", "1.0e3"),
      "physics.bh_seed_halo_mass_threshold_code");
  frozen.config.physics.bh_seed_mass_code = parseFloating(
      requireString(entries, consumed, "physics.bh_seed_mass_code", "1.0"),
      "physics.bh_seed_mass_code");
  frozen.config.physics.bh_seed_max_per_cell = parseNumber<std::uint32_t>(
      requireString(entries, consumed, "physics.bh_seed_max_per_cell", "1"),
      "physics.bh_seed_max_per_cell");
  frozen.config.physics.bh_alpha_bondi = parseFloating(
      requireString(entries, consumed, "physics.bh_alpha_bondi", "1.0"),
      "physics.bh_alpha_bondi");
  frozen.config.physics.bh_use_eddington_cap = parseBool(
      requireString(entries, consumed, "physics.bh_use_eddington_cap", "true"),
      "physics.bh_use_eddington_cap");
  frozen.config.physics.bh_epsilon_r = parseFloating(
      requireString(entries, consumed, "physics.bh_epsilon_r", "0.1"),
      "physics.bh_epsilon_r");
  frozen.config.physics.bh_epsilon_f = parseFloating(
      requireString(entries, consumed, "physics.bh_epsilon_f", "0.05"),
      "physics.bh_epsilon_f");
  frozen.config.physics.bh_feedback_coupling_efficiency = parseFloating(
      requireString(entries, consumed, "physics.bh_feedback_coupling_efficiency", "1.0"),
      "physics.bh_feedback_coupling_efficiency");
  frozen.config.physics.bh_duty_cycle_active_edd_ratio_threshold = parseFloating(
      requireString(entries, consumed, "physics.bh_duty_cycle_active_edd_ratio_threshold", "0.01"),
      "physics.bh_duty_cycle_active_edd_ratio_threshold");
  frozen.config.physics.bh_proton_mass_si = parseFloating(
      requireString(entries, consumed, "physics.bh_proton_mass_si", "1.67262192369e-27"),
      "physics.bh_proton_mass_si");
  frozen.config.physics.bh_thomson_cross_section_si = parseFloating(
      requireString(entries, consumed, "physics.bh_thomson_cross_section_si", "6.6524587321e-29"),
      "physics.bh_thomson_cross_section_si");
  frozen.config.physics.bh_newton_g_si = parseFloating(
      requireString(entries, consumed, "physics.bh_newton_g_si", "6.67430e-11"),
      "physics.bh_newton_g_si");
  frozen.config.physics.bh_speed_of_light_si = parseFloating(
      requireString(entries, consumed, "physics.bh_speed_of_light_si", "2.99792458e8"),
      "physics.bh_speed_of_light_si");
  frozen.config.physics.enable_tracers = parseBool(
      requireString(entries, consumed, "physics.enable_tracers", "false"),
      "physics.enable_tracers");
  frozen.config.physics.tracer_track_mass = parseBool(
      requireString(entries, consumed, "physics.tracer_track_mass", "true"),
      "physics.tracer_track_mass");
  frozen.config.physics.tracer_min_host_mass_code = parseFloating(
      requireString(entries, consumed, "physics.tracer_min_host_mass_code", "0.0"),
      "physics.tracer_min_host_mass_code");
  frozen.config.physics.thermochemistry_backend = parseThermochemistryBackendMode(
      requireString(entries, consumed, "physics.thermochemistry_backend", "legacy_cooling"));
  frozen.config.physics.metal_cooling_mode = parseMetalCoolingModeConfig(
      requireString(entries, consumed, "physics.metal_cooling_mode", "none"));
  frozen.config.physics.metal_cooling_table_range = parseMetalCoolingTableRangeConfig(
      requireString(entries, consumed, "physics.metal_cooling_table_range", "reject"));
  frozen.config.physics.dust_model = parseDustModelConfig(
      requireString(entries, consumed, "physics.dust_model", "none"));
  frozen.config.physics.ism_regime_model = parseIsmRegimeModelConfig(
      requireString(entries, consumed, "physics.ism_regime_model", "legacy"));
  frozen.config.physics.sgs_turbulence_model = parseSgsTurbulenceModelConfig(
      requireString(entries, consumed, "physics.sgs_turbulence_model", "none"));
  frozen.config.physics.chemistry_adaptive_selection = parseChemistryAdaptiveSelectionConfig(
      requireString(entries, consumed, "physics.chemistry_adaptive_selection", "off"));
  frozen.config.physics.stellar_population_model = parseStellarPopulationModelConfig(
      requireString(entries, consumed, "physics.stellar_population_model", "legacy_table"));
  frozen.config.physics.fb_supernova_model = parseSupernovaFeedbackModelConfig(
      requireString(entries, consumed, "physics.fb_supernova_model", "legacy"));
  frozen.config.physics.thermochemistry_relative_tolerance = parseFloating(
      requireString(entries, consumed, "physics.thermochemistry_relative_tolerance", "1.0e-4"), "physics.thermochemistry_relative_tolerance");
  frozen.config.physics.thermochemistry_abundance_absolute_tolerance = parseFloating(
      requireString(entries, consumed, "physics.thermochemistry_abundance_absolute_tolerance", "1.0e-12"), "physics.thermochemistry_abundance_absolute_tolerance");
  frozen.config.physics.thermochemistry_analytic_relaxation_time_years = parseFloating(
      requireString(entries, consumed, "physics.thermochemistry_analytic_relaxation_time_years", "1.0e6"), "physics.thermochemistry_analytic_relaxation_time_years");
  frozen.config.physics.thermochemistry_analytic_equilibrium_temperature_k = parseFloating(
      requireString(entries, consumed, "physics.thermochemistry_analytic_equilibrium_temperature_k", "1.0e4"), "physics.thermochemistry_analytic_equilibrium_temperature_k");
  frozen.config.physics.isrf_g0_habing = parseFloating(
      requireString(entries, consumed, "physics.isrf_g0_habing", "0.0"), "physics.isrf_g0_habing");
  frozen.config.physics.chemistry_shielding_length_max_pc = parseFloating(
      requireString(entries, consumed, "physics.chemistry_shielding_length_max_pc", "100.0"), "physics.chemistry_shielding_length_max_pc");
  frozen.config.physics.metal_cooling_switch_temperature_k = parseFloating(
      requireString(entries, consumed, "physics.metal_cooling_switch_temperature_k", "1.0e4"), "physics.metal_cooling_switch_temperature_k");
  frozen.config.physics.chemistry_primordial_helium_mass_fraction = parseFloating(
      requireString(entries, consumed, "physics.chemistry_primordial_helium_mass_fraction", "0.2449"), "physics.chemistry_primordial_helium_mass_fraction");
  frozen.config.physics.dust_accretion_timescale_ref_years = parseFloating(
      requireString(entries, consumed, "physics.dust_accretion_timescale_ref_years", "2.0e7"), "physics.dust_accretion_timescale_ref_years");
  frozen.config.physics.dust_sputtering_timescale_ref_years = parseFloating(
      requireString(entries, consumed, "physics.dust_sputtering_timescale_ref_years", "1.7e8"), "physics.dust_sputtering_timescale_ref_years");
  frozen.config.physics.dust_sn_swept_mass_msun = parseFloating(
      requireString(entries, consumed, "physics.dust_sn_swept_mass_msun", "6800.0"), "physics.dust_sn_swept_mass_msun");
  frozen.config.physics.dust_sn_destruction_efficiency = parseFloating(
      requireString(entries, consumed, "physics.dust_sn_destruction_efficiency", "0.3"), "physics.dust_sn_destruction_efficiency");
  frozen.config.physics.dust_agb_condensation_efficiency = parseFloating(
      requireString(entries, consumed, "physics.dust_agb_condensation_efficiency", "0.2"), "physics.dust_agb_condensation_efficiency");
  frozen.config.physics.dust_cc_condensation_efficiency = parseFloating(
      requireString(entries, consumed, "physics.dust_cc_condensation_efficiency", "0.1"), "physics.dust_cc_condensation_efficiency");
  frozen.config.physics.ism_jeans_resolution_cells = parseFloating(
      requireString(entries, consumed, "physics.ism_jeans_resolution_cells", "4.0"), "physics.ism_jeans_resolution_cells");
  frozen.config.physics.ism_regime_hysteresis_fraction = parseFloating(
      requireString(entries, consumed, "physics.ism_regime_hysteresis_fraction", "0.2"), "physics.ism_regime_hysteresis_fraction");
  frozen.config.physics.ism_regime_min_dwell_years = parseFloating(
      requireString(entries, consumed, "physics.ism_regime_min_dwell_years", "1.0e6"), "physics.ism_regime_min_dwell_years");
  frozen.config.physics.sgs_dissipation_coefficient = parseFloating(
      requireString(entries, consumed, "physics.sgs_dissipation_coefficient", "1.58"), "physics.sgs_dissipation_coefficient");
  frozen.config.physics.sgs_eddy_viscosity_coefficient = parseFloating(
      requireString(entries, consumed, "physics.sgs_eddy_viscosity_coefficient", "0.05"), "physics.sgs_eddy_viscosity_coefficient");
  frozen.config.physics.sgs_feedback_energy_fraction = parseFloating(
      requireString(entries, consumed, "physics.sgs_feedback_energy_fraction", "0.0"), "physics.sgs_feedback_energy_fraction");
  frozen.config.physics.chemistry_adaptive_hysteresis_fraction = parseFloating(
      requireString(entries, consumed, "physics.chemistry_adaptive_hysteresis_fraction", "0.3"), "physics.chemistry_adaptive_hysteresis_fraction");
  frozen.config.physics.chemistry_adaptive_min_dwell_years = parseFloating(
      requireString(entries, consumed, "physics.chemistry_adaptive_min_dwell_years", "1.0e6"), "physics.chemistry_adaptive_min_dwell_years");
  frozen.config.physics.chemistry_adaptive_cold_temperature_k = parseFloating(
      requireString(entries, consumed, "physics.chemistry_adaptive_cold_temperature_k", "3.0e3"), "physics.chemistry_adaptive_cold_temperature_k");
  frozen.config.physics.chemistry_adaptive_cooling_time_ratio = parseFloating(
      requireString(entries, consumed, "physics.chemistry_adaptive_cooling_time_ratio", "0.0"),
      "physics.chemistry_adaptive_cooling_time_ratio");
  frozen.config.physics.chemistry_adaptive_dense_nh_cm3 = parseFloating(
      requireString(entries, consumed, "physics.chemistry_adaptive_dense_nh_cm3", "1.0"), "physics.chemistry_adaptive_dense_nh_cm3");
  frozen.config.physics.chemistry_shadow_sample_fraction = parseFloating(
      requireString(entries, consumed, "physics.chemistry_shadow_sample_fraction", "0.0"), "physics.chemistry_shadow_sample_fraction");
  frozen.config.physics.ssp_imf_mass_min_msun = parseFloating(
      requireString(entries, consumed, "physics.ssp_imf_mass_min_msun", "0.1"), "physics.ssp_imf_mass_min_msun");
  frozen.config.physics.ssp_imf_mass_max_msun = parseFloating(
      requireString(entries, consumed, "physics.ssp_imf_mass_max_msun", "100.0"), "physics.ssp_imf_mass_max_msun");
  frozen.config.physics.ssp_cc_mass_min_msun = parseFloating(
      requireString(entries, consumed, "physics.ssp_cc_mass_min_msun", "8.0"), "physics.ssp_cc_mass_min_msun");
  frozen.config.physics.ssp_cc_mass_max_msun = parseFloating(
      requireString(entries, consumed, "physics.ssp_cc_mass_max_msun", "40.0"), "physics.ssp_cc_mass_max_msun");
  frozen.config.physics.ssp_snia_per_msun = parseFloating(
      requireString(entries, consumed, "physics.ssp_snia_per_msun", "1.3e-3"), "physics.ssp_snia_per_msun");
  frozen.config.physics.ssp_snia_delay_min_years = parseFloating(
      requireString(entries, consumed, "physics.ssp_snia_delay_min_years", "4.0e7"), "physics.ssp_snia_delay_min_years");
  frozen.config.physics.ssp_snia_delay_max_years = parseFloating(
      requireString(entries, consumed, "physics.ssp_snia_delay_max_years", "1.37e10"), "physics.ssp_snia_delay_max_years");
  frozen.config.physics.ssp_energy_per_sn_erg = parseFloating(
      requireString(entries, consumed, "physics.ssp_energy_per_sn_erg", "1.0e51"), "physics.ssp_energy_per_sn_erg");
  frozen.config.physics.ssp_reference_net_metal_yield = parseFloating(
      requireString(entries, consumed, "physics.ssp_reference_net_metal_yield", "0.0"), "physics.ssp_reference_net_metal_yield");
  frozen.config.physics.fb_presn_wind_momentum_lsun_per_c = parseFloating(
      requireString(entries, consumed, "physics.fb_presn_wind_momentum_lsun_per_c", "0.0"), "physics.fb_presn_wind_momentum_lsun_per_c");
  frozen.config.physics.thermochemistry_max_substeps = parseNumber<std::uint32_t>(
      requireString(entries, consumed, "physics.thermochemistry_max_substeps", "4096"), "physics.thermochemistry_max_substeps");
  frozen.config.physics.thermochemistry_max_newton_iterations = parseNumber<std::uint32_t>(
      requireString(entries, consumed, "physics.thermochemistry_max_newton_iterations", "12"), "physics.thermochemistry_max_newton_iterations");
  frozen.config.physics.chemistry_shadow_max_cells = parseNumber<std::uint32_t>(
      requireString(entries, consumed, "physics.chemistry_shadow_max_cells", "0"), "physics.chemistry_shadow_max_cells");
  frozen.config.physics.uv_background_table_path = requireString(entries, consumed, "physics.uv_background_table_path", "");
  frozen.config.physics.uv_background_table_sha256 = requireString(entries, consumed, "physics.uv_background_table_sha256", "");
  frozen.config.physics.metal_cooling_table_path = requireString(entries, consumed, "physics.metal_cooling_table_path", "");
  frozen.config.physics.metal_cooling_table_sha256 = requireString(entries, consumed, "physics.metal_cooling_table_sha256", "");
  frozen.config.physics.metal_cooling_accept_uv_background_mismatch = parseBool(
      requireString(entries, consumed, "physics.metal_cooling_accept_uv_background_mismatch", "false"),
      "physics.metal_cooling_accept_uv_background_mismatch");
  frozen.config.physics.ssp_yield_table_path = requireString(entries, consumed, "physics.ssp_yield_table_path", "");
  frozen.config.physics.ssp_yield_table_sha256 = requireString(entries, consumed, "physics.ssp_yield_table_sha256", "");
  frozen.config.physics.chemistry_metallicity_scaled_dust_proxy = parseBool(
      requireString(entries, consumed, "physics.chemistry_metallicity_scaled_dust_proxy", "false"), "physics.chemistry_metallicity_scaled_dust_proxy");
  frozen.config.physics.thermochemistry_allow_unqualified_rates = parseBool(
      requireString(entries, consumed, "physics.thermochemistry_allow_unqualified_rates", "false"), "physics.thermochemistry_allow_unqualified_rates");
  frozen.config.physics.ssp_stochastic_sn_events = parseBool(
      requireString(entries, consumed, "physics.ssp_stochastic_sn_events", "true"), "physics.ssp_stochastic_sn_events");
  frozen.config.physics.ssp_clustered_associations = parseBool(
      requireString(entries, consumed, "physics.ssp_clustered_associations", "false"), "physics.ssp_clustered_associations");
  frozen.config.physics.ssp_association_radius_code = parseFloating(
      requireString(entries, consumed, "physics.ssp_association_radius_code", "0.0"), "physics.ssp_association_radius_code");
  frozen.config.physics.ssp_association_time_window_years = parseFloating(
      requireString(entries, consumed, "physics.ssp_association_time_window_years", "3.0e6"),
      "physics.ssp_association_time_window_years");
  frozen.config.physics.ssp_wind_model = parseSspWindModelConfig(
      requireString(entries, consumed, "physics.ssp_wind_model", "none"));
  frozen.config.physics.ssp_wind_table_path =
      requireString(entries, consumed, "physics.ssp_wind_table_path", "");
  frozen.config.physics.ssp_wind_table_sha256 =
      requireString(entries, consumed, "physics.ssp_wind_table_sha256", "");

  frozen.config.output.run_name =
      requireString(entries, consumed, "output.run_name", frozen.config.output.run_name);
  frozen.config.output.output_directory =
      requireString(entries, consumed, "output.output_directory", frozen.config.output.output_directory);
  frozen.config.output.output_stem = sanitizeStem(
      requireString(entries, consumed, "output.output_stem", frozen.config.output.output_stem),
      "output.output_stem");
  frozen.config.output.restart_stem = sanitizeStem(
      requireString(entries, consumed, "output.restart_stem", frozen.config.output.restart_stem),
      "output.restart_stem");
  frozen.config.output.snapshot_interval_steps = parseNumber<int>(
      requireString(entries, consumed, "output.snapshot_interval_steps", "64"),
      "output.snapshot_interval_steps");
  frozen.config.output.snapshot_interval_time_code = parseFloating(
      requireString(entries, consumed, "output.snapshot_interval_time_code", "0.0"),
      "output.snapshot_interval_time_code");
  frozen.config.output.snapshot_layout = parseScienceSnapshotLayout(
      requireString(entries, consumed, "output.snapshot_layout", "auto"));
  frozen.config.output.snapshot_num_files = parseNumber<std::uint32_t>(
      requireString(entries, consumed, "output.snapshot_num_files", "0"),
      "output.snapshot_num_files");
  frozen.config.output.write_restarts = parseBool(
      requireString(entries, consumed, "output.write_restarts", "true"), "output.write_restarts");
  frozen.config.output.restart_interval_steps = parseNumber<int>(
      requireString(entries, consumed, "output.restart_interval_steps",
                    defaultFor("output.restart_interval_steps")),
      "output.restart_interval_steps");
  frozen.config.output.restart_interval_time_code = parseFloating(
      requireString(entries, consumed, "output.restart_interval_time_code",
                    defaultFor("output.restart_interval_time_code")),
      "output.restart_interval_time_code");
  frozen.config.output.snapshot_verification = parseSnapshotVerificationPolicy(
      requireString(entries, consumed, "output.snapshot_verification",
                    defaultFor("output.snapshot_verification")));

  frozen.config.parallel.mpi_ranks_expected = parseNumber<int>(
      requireString(entries, consumed, "parallel.mpi_ranks_expected", "1"),
      "parallel.mpi_ranks_expected");
  frozen.config.parallel.omp_threads = parseNumber<int>(
      requireString(entries, consumed, "parallel.omp_threads", "1"), "parallel.omp_threads");
  frozen.config.parallel.gpu_devices = parseNumber<int>(
      requireString(entries, consumed, "parallel.gpu_devices", "0"), "parallel.gpu_devices");
  frozen.config.parallel.deterministic_reduction =
      parseBool(requireString(entries, consumed, "parallel.deterministic_reduction", "true"),
                "parallel.deterministic_reduction");
  frozen.config.parallel.decomposition_particle_count_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_particle_count_weight",
                    defaultFor("parallel.decomposition_particle_count_weight")),
      "parallel.decomposition_particle_count_weight");
  frozen.config.parallel.decomposition_gas_cell_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_gas_cell_weight",
                    defaultFor("parallel.decomposition_gas_cell_weight")),
      "parallel.decomposition_gas_cell_weight");
  frozen.config.parallel.decomposition_tree_interaction_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_tree_interaction_weight",
                    defaultFor("parallel.decomposition_tree_interaction_weight")),
      "parallel.decomposition_tree_interaction_weight");
  frozen.config.parallel.decomposition_pm_mesh_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_pm_mesh_weight",
                    defaultFor("parallel.decomposition_pm_mesh_weight")),
      "parallel.decomposition_pm_mesh_weight");
  frozen.config.parallel.decomposition_amr_patch_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_amr_patch_weight",
                    defaultFor("parallel.decomposition_amr_patch_weight")),
      "parallel.decomposition_amr_patch_weight");
  frozen.config.parallel.decomposition_active_fraction_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_active_fraction_weight",
                    defaultFor("parallel.decomposition_active_fraction_weight")),
      "parallel.decomposition_active_fraction_weight");
  frozen.config.parallel.decomposition_memory_pressure_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_memory_pressure_weight",
                    defaultFor("parallel.decomposition_memory_pressure_weight")),
      "parallel.decomposition_memory_pressure_weight");
  frozen.config.parallel.decomposition_gpu_occupancy_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_gpu_occupancy_weight",
                    defaultFor("parallel.decomposition_gpu_occupancy_weight")),
      "parallel.decomposition_gpu_occupancy_weight");
  frozen.config.parallel.decomposition_generic_work_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_generic_work_weight",
                    defaultFor("parallel.decomposition_generic_work_weight")),
      "parallel.decomposition_generic_work_weight");
  frozen.config.parallel.decomposition_runtime_rebalance_enabled = parseBool(
      requireString(entries, consumed, "parallel.decomposition_runtime_rebalance_enabled",
                    defaultFor("parallel.decomposition_runtime_rebalance_enabled")),
      "parallel.decomposition_runtime_rebalance_enabled");
  frozen.config.parallel.decomposition_spatial_work_enabled = parseBool(
      requireString(entries, consumed, "parallel.decomposition_spatial_work_enabled",
                    defaultFor("parallel.decomposition_spatial_work_enabled")),
      "parallel.decomposition_spatial_work_enabled");
  frozen.config.parallel.decomposition_debug_exact_ownership_audit = parseBool(
      requireString(entries, consumed, "parallel.decomposition_debug_exact_ownership_audit",
                    defaultFor("parallel.decomposition_debug_exact_ownership_audit")),
      "parallel.decomposition_debug_exact_ownership_audit");
  frozen.config.parallel.decomposition_rebalance_imbalance_trigger = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_rebalance_imbalance_trigger",
                    defaultFor("parallel.decomposition_rebalance_imbalance_trigger")),
      "parallel.decomposition_rebalance_imbalance_trigger");
  frozen.config.parallel.decomposition_rebalance_memory_trigger = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_rebalance_memory_trigger",
                    defaultFor("parallel.decomposition_rebalance_memory_trigger")),
      "parallel.decomposition_rebalance_memory_trigger");
  frozen.config.parallel.decomposition_rebalance_max_migrated_load_fraction = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_rebalance_max_migrated_load_fraction",
                    defaultFor("parallel.decomposition_rebalance_max_migrated_load_fraction")),
      "parallel.decomposition_rebalance_max_migrated_load_fraction");
  frozen.config.parallel.decomposition_measured_tree_pair_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_tree_pair_weight",
                    defaultFor("parallel.decomposition_measured_tree_pair_weight")),
      "parallel.decomposition_measured_tree_pair_weight");
  frozen.config.parallel.decomposition_measured_pm_cell_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_pm_cell_weight",
                    defaultFor("parallel.decomposition_measured_pm_cell_weight")),
      "parallel.decomposition_measured_pm_cell_weight");
  frozen.config.parallel.decomposition_measured_amr_cell_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_amr_cell_weight",
                    defaultFor("parallel.decomposition_measured_amr_cell_weight")),
      "parallel.decomposition_measured_amr_cell_weight");
  frozen.config.parallel.decomposition_measured_hydro_face_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_hydro_face_weight",
                    defaultFor("parallel.decomposition_measured_hydro_face_weight")),
      "parallel.decomposition_measured_hydro_face_weight");
  frozen.config.parallel.decomposition_measured_wall_ms_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_wall_ms_weight",
                    defaultFor("parallel.decomposition_measured_wall_ms_weight")),
      "parallel.decomposition_measured_wall_ms_weight");
  frozen.config.parallel.isolated_pm_root_workspace_limit_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.isolated_pm_root_workspace_limit_bytes",
                    defaultFor("parallel.isolated_pm_root_workspace_limit_bytes")),
      "parallel.isolated_pm_root_workspace_limit_bytes");
  frozen.config.parallel.zoom_high_res_allgather_limit_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.zoom_high_res_allgather_limit_bytes",
                    defaultFor("parallel.zoom_high_res_allgather_limit_bytes")),
      "parallel.zoom_high_res_allgather_limit_bytes");
  frozen.config.parallel.gravity_memory_budget_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.gravity_memory_budget_bytes",
                    defaultFor("parallel.gravity_memory_budget_bytes")),
      "parallel.gravity_memory_budget_bytes");
  frozen.config.parallel.gravity_backend_unknown_reserve_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.gravity_backend_unknown_reserve_bytes",
                    defaultFor("parallel.gravity_backend_unknown_reserve_bytes")),
      "parallel.gravity_backend_unknown_reserve_bytes");
  frozen.config.parallel.gravity_memory_safety_margin_fraction = parseNumber<double>(
      requireString(entries, consumed, "parallel.gravity_memory_safety_margin_fraction",
                    defaultFor("parallel.gravity_memory_safety_margin_fraction")),
      "parallel.gravity_memory_safety_margin_fraction");
  frozen.config.parallel.process_memory_budget_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.process_memory_budget_bytes",
                    defaultFor("parallel.process_memory_budget_bytes")),
      "parallel.process_memory_budget_bytes");
  frozen.config.parallel.process_mpi_unknown_reserve_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.process_mpi_unknown_reserve_bytes",
                    defaultFor("parallel.process_mpi_unknown_reserve_bytes")),
      "parallel.process_mpi_unknown_reserve_bytes");
  frozen.config.parallel.process_fftw_unknown_reserve_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.process_fftw_unknown_reserve_bytes",
                    defaultFor("parallel.process_fftw_unknown_reserve_bytes")),
      "parallel.process_fftw_unknown_reserve_bytes");
  frozen.config.parallel.process_hdf5_unknown_reserve_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.process_hdf5_unknown_reserve_bytes",
                    defaultFor("parallel.process_hdf5_unknown_reserve_bytes")),
      "parallel.process_hdf5_unknown_reserve_bytes");
  frozen.config.parallel.process_allocator_unknown_reserve_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.process_allocator_unknown_reserve_bytes",
                    defaultFor("parallel.process_allocator_unknown_reserve_bytes")),
      "parallel.process_allocator_unknown_reserve_bytes");
  frozen.config.parallel.process_output_restart_overlap_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.process_output_restart_overlap_bytes",
                    defaultFor("parallel.process_output_restart_overlap_bytes")),
      "parallel.process_output_restart_overlap_bytes");
  frozen.config.parallel.process_memory_safety_margin_fraction = parseNumber<double>(
      requireString(entries, consumed, "parallel.process_memory_safety_margin_fraction",
                    defaultFor("parallel.process_memory_safety_margin_fraction")),
      "parallel.process_memory_safety_margin_fraction");
  frozen.config.parallel.node_memory_admission = parseNodeMemoryAdmissionMode(requireString(
      entries, consumed, "parallel.node_memory_admission", defaultFor("parallel.node_memory_admission")));
  frozen.config.parallel.node_memory_budget_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.node_memory_budget_bytes",
                    defaultFor("parallel.node_memory_budget_bytes")),
      "parallel.node_memory_budget_bytes");
  frozen.config.parallel.node_memory_reserve_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.node_memory_reserve_bytes",
                    defaultFor("parallel.node_memory_reserve_bytes")),
      "parallel.node_memory_reserve_bytes");

  frozen.config.analysis.enable_diagnostics = parseBool(
      requireString(entries, consumed, "analysis.enable_diagnostics", "true"),
      "analysis.enable_diagnostics");
  frozen.config.analysis.enable_halo_workflow = parseBool(
      requireString(entries, consumed, "analysis.enable_halo_workflow", "false"),
      "analysis.enable_halo_workflow");
  frozen.config.analysis.halo_on_the_fly = parseBool(
      requireString(entries, consumed, "analysis.halo_on_the_fly", "false"),
      "analysis.halo_on_the_fly");
  frozen.config.analysis.diagnostics_execution_policy = parseDiagnosticsExecutionPolicy(requireString(
      entries,
      consumed,
      "analysis.diagnostics_execution_policy",
      diagnosticsExecutionPolicyToString(frozen.config.analysis.diagnostics_execution_policy)));
  frozen.config.analysis.run_health_interval_steps = parseNumber<int>(
      requireString(entries, consumed, "analysis.run_health_interval_steps", "1"),
      "analysis.run_health_interval_steps");
  frozen.config.analysis.science_light_interval_steps = parseNumber<int>(
      requireString(entries, consumed, "analysis.science_light_interval_steps", "8"),
      "analysis.science_light_interval_steps");
  frozen.config.analysis.science_heavy_interval_steps = parseNumber<int>(
      requireString(entries, consumed, "analysis.science_heavy_interval_steps", "64"),
      "analysis.science_heavy_interval_steps");
  frozen.config.analysis.retention_bundle_count = parseNumber<int>(
      requireString(entries, consumed, "analysis.retention_bundle_count", "8"),
      "analysis.retention_bundle_count");
  frozen.config.analysis.power_spectrum_mesh_n = parseNumber<int>(
      requireString(entries, consumed, "analysis.power_spectrum_mesh_n", "16"),
      "analysis.power_spectrum_mesh_n");
  frozen.config.analysis.power_spectrum_bin_count = parseNumber<int>(
      requireString(entries, consumed, "analysis.power_spectrum_bin_count", "12"),
      "analysis.power_spectrum_bin_count");
  frozen.config.analysis.sf_history_bin_count = parseNumber<int>(
      requireString(entries, consumed, "analysis.sf_history_bin_count", "16"),
      "analysis.sf_history_bin_count");
  frozen.config.analysis.quicklook_grid_n = parseNumber<int>(
      requireString(entries, consumed, "analysis.quicklook_grid_n", "32"),
      "analysis.quicklook_grid_n");
  frozen.config.analysis.diagnostics_stem = sanitizeStem(
      requireString(entries, consumed, "analysis.diagnostics_stem", frozen.config.analysis.diagnostics_stem),
      "analysis.diagnostics_stem");
  frozen.config.analysis.halo_catalog_stem = sanitizeStem(
      requireString(entries, consumed, "analysis.halo_catalog_stem", frozen.config.analysis.halo_catalog_stem),
      "analysis.halo_catalog_stem");
  frozen.config.analysis.merger_tree_stem = sanitizeStem(
      requireString(entries, consumed, "analysis.merger_tree_stem", frozen.config.analysis.merger_tree_stem),
      "analysis.merger_tree_stem");
  frozen.config.analysis.halo_fof_linking_length_factor = parseFloating(
      requireString(
          entries,
          consumed,
          "analysis.halo_fof_linking_length_factor",
          "0.2"),
      "analysis.halo_fof_linking_length_factor");
  frozen.config.analysis.halo_fof_min_group_size = parseNumber<int>(
      requireString(entries, consumed, "analysis.halo_fof_min_group_size", "16"),
      "analysis.halo_fof_min_group_size");
  frozen.config.analysis.halo_include_gas = parseBool(
      requireString(entries, consumed, "analysis.halo_include_gas", "true"),
      "analysis.halo_include_gas");
  frozen.config.analysis.halo_include_stars = parseBool(
      requireString(entries, consumed, "analysis.halo_include_stars", "true"),
      "analysis.halo_include_stars");
  frozen.config.analysis.halo_include_black_holes = parseBool(
      requireString(entries, consumed, "analysis.halo_include_black_holes", "true"),
      "analysis.halo_include_black_holes");

  const bool compatible_by_file = parseBool(
      requireString(entries, consumed, "compatibility.allow_unknown_keys", "false"),
      "compatibility.allow_unknown_keys");
  frozen.config.compatibility.allow_unknown_keys = options.allow_unknown_keys || compatible_by_file;

  // Resolve one authoritative integration endpoint before validation/provenance.
  // Cosmological configurations retain the declared scale-factor default.  An
  // otherwise-unqualified code-time endpoint selects code-time authority,
  // while non-cosmological modes use code time unless explicitly rejected
  // below.  Supplying competing endpoint families without selecting an
  // authority is ambiguous and therefore fails closed.
  const ModePolicy resolved_mode_policy = buildModePolicy(frozen.config.mode);
  if (!entries.contains("numerics.integrator_time_variable")) {
    if (!resolved_mode_policy.cosmological_comoving_frame) {
      frozen.config.numerics.integrator_time_variable =
          IntegratorTimeVariable::kCodeTime;
    } else {
      const bool supplied_scale_endpoint =
          entries.contains("numerics.a_end") || entries.contains("numerics.z_end");
      const bool supplied_code_endpoint = entries.contains("numerics.t_code_end");
      if (supplied_scale_endpoint && supplied_code_endpoint) {
        throw ConfigError(
            "cosmological config supplies both scale-factor/redshift and code-time endpoints without selecting numerics.integrator_time_variable");
      }
      if (supplied_scale_endpoint) {
        frozen.config.numerics.integrator_time_variable =
            IntegratorTimeVariable::kScaleFactor;
      } else {
        // With no explicit endpoint family, preserve the historical default
        // finite code-time span rather than inventing a zero-length a=1 step.
        frozen.config.numerics.integrator_time_variable =
            IntegratorTimeVariable::kCodeTime;
      }
    }
  }
  const auto variable = frozen.config.numerics.integrator_time_variable;
  if (variable == IntegratorTimeVariable::kPhysicalTime) {
    throw ConfigError(
        "numerics.integrator_time_variable=physical_time is not yet supported; use scale_factor, ln_a, or code_time");
  }
  if (!resolved_mode_policy.cosmological_comoving_frame &&
      variable != IntegratorTimeVariable::kCodeTime) {
    throw ConfigError(
        "non-cosmological modes currently support only numerics.integrator_time_variable=code_time");
  }
  if (resolved_mode_policy.cosmological_comoving_frame) {
    CosmologyBackgroundConfig background_config;
    background_config.hubble_param = frozen.config.cosmology.hubble_param;
    background_config.omega_matter = frozen.config.cosmology.omega_matter;
    background_config.omega_lambda = frozen.config.cosmology.omega_lambda;
    const LambdaCdmBackground background = [&]() {
      try {
        return LambdaCdmBackground(background_config);
      } catch (const std::exception& error) {
        throw ConfigError(
            "invalid cosmological background while resolving integration authority in " +
            source_name + ": " + error.what());
      }
    }();
    const UnitSystem units = makeUnitSystem(
        frozen.config.units.length_unit, frozen.config.units.mass_unit,
        frozen.config.units.velocity_unit);

    if (variable == IntegratorTimeVariable::kScaleFactor ||
        variable == IntegratorTimeVariable::kLogScaleFactor) {
      const double elapsed_code = background.cosmicTimeIntervalSi(
          frozen.config.numerics.a_begin, frozen.config.numerics.a_end) /
          units.timeSiPerCode();
      if (!std::isfinite(elapsed_code) || elapsed_code <= 0.0) {
        throw ConfigError(
            "resolved FLRW scale-factor interval is not finite and positive in " + source_name + ": a_begin=" +
            std::to_string(frozen.config.numerics.a_begin) + ", a_end=" +
            std::to_string(frozen.config.numerics.a_end) + ", elapsed_code=" +
            std::to_string(elapsed_code));
      }
      const double resolved_end = frozen.config.numerics.t_code_begin + elapsed_code;
      if (entries.contains("numerics.t_code_end")) {
        const double supplied = frozen.config.numerics.t_code_end;
        const double scale = std::max({1.0, std::abs(supplied), std::abs(resolved_end)});
        if (std::abs(supplied - resolved_end) > 1.0e-10 * scale) {
          throw ConfigError(
              "numerics.t_code_end contradicts the FLRW endpoint implied by authoritative a_begin/a_end");
        }
      }
      frozen.config.numerics.t_code_end = resolved_end;
    } else if (variable == IntegratorTimeVariable::kCodeTime) {
      const double elapsed_code =
          frozen.config.numerics.t_code_end - frozen.config.numerics.t_code_begin;
      if (!std::isfinite(elapsed_code) || elapsed_code <= 0.0) {
        throw ConfigError(
            "code-time authority requires a finite positive t_code_end - t_code_begin interval");
      }
      double resolved_a_end = 0.0;
      try {
        resolved_a_end = resolveFlatLambdaCdmScaleFactorAfterTime(
            background, frozen.config.numerics.a_begin,
            elapsed_code * units.timeSiPerCode());
      } catch (const std::exception& error) {
        throw ConfigError(
            "failed to resolve code-time FLRW endpoint in " + source_name +
            ": " + error.what());
      }
      if (!std::isfinite(resolved_a_end) ||
          resolved_a_end < frozen.config.numerics.a_begin) {
        throw ConfigError(
            "code-time authority produced an invalid FLRW final scale factor");
      }
      if (entries.contains("numerics.a_end") || entries.contains("numerics.z_end")) {
        const double supplied = frozen.config.numerics.a_end;
        const double scale = std::max({1.0, std::abs(supplied), std::abs(resolved_a_end)});
        if (std::abs(supplied - resolved_a_end) > 1.0e-10 * scale) {
          throw ConfigError(
              "numerics.a_end/z_end contradicts the FLRW endpoint implied by authoritative code time");
        }
      }
      frozen.config.numerics.a_end = resolved_a_end;
      frozen.config.numerics.z_end =
          scaleFactorToRedshift(resolved_a_end, "numerics.a_end");
    }
  }


  for (const auto& [key, entry] : entries) {
    UserConfigEntry user_entry;
    user_entry.canonical_key = key;
    user_entry.value = trim(entry.value);
    user_entry.source_key = key;
    user_entry.source_line = entry.line_number;
    frozen.user_config.entries.push_back(std::move(user_entry));
  }

  validateConfig(frozen.config);

  std::vector<std::string> unknown;
  std::set<std::string> known_keys;
  for (const auto& entry : configKeyRegistry()) {
    known_keys.insert(entry.key);
  }
  for (const auto& [legacy_key, _] : deprecatedAliasRegistry()) {
    known_keys.insert(legacy_key);
  }

  for (const auto& [key, _] : entries) {
    if (consumed.contains(key)) {
      continue;
    }
    if (!known_keys.contains(key)) {
      unknown.push_back(key);
      continue;
    }
    throw ConfigError("internal config parsing error: known key '" + key + "' was not consumed");
  }

  if (!unknown.empty() && !frozen.config.compatibility.allow_unknown_keys) {
    std::ostringstream stream;
    stream << "unknown parameter keys:";
    for (const std::string& key : unknown) {
      stream << ' ' << key;
    }
    stream << " (set compatibility.allow_unknown_keys=true to bypass)";
    throw ConfigError(stream.str());
  }

  frozen.normalized_text = buildNormalizedText(frozen);
  frozen.provenance.config_hash = stableConfigHash(frozen.normalized_text);
  frozen.provenance.config_hash_hex = stableConfigHashHex(frozen.normalized_text);

  return frozen;
}

}  // namespace

ConfigError::ConfigError(const std::string& message) : std::runtime_error(message) {}

SimulationConfig makeUnvalidatedSimulationConfigForTests() {
  SimulationConfig config(SimulationConfig::ConstructionToken{});
  // Unit/integration fixtures are allowed to exercise the tiny deterministic
  // naive-DFT backend on builds without FFTW. Parsed production decks remain
  // fail-closed unless they explicitly opt into the same diagnostic override.
  config.numerics.treepm_allow_diagnostic_naive_dft = true;
  return config;
}

FrozenConfig::FrozenConfig()
    : config(SimulationConfig::ConstructionToken{}) {}

FrozenConfig loadFrozenConfigFromFile(const std::filesystem::path& path, const ParseOptions& options) {
  std::ifstream input(path);
  if (!input) {
    throw ConfigError("failed to open config file: " + path.string());
  }
  std::ostringstream contents;
  contents << input.rdbuf();
  return loadFrozenConfigFromString(contents.str(), path.string(), options);
}

FrozenConfig loadFrozenConfigFromString(
    const std::string& config_text,
    const std::string& source_name,
    const ParseOptions& options) {
  const auto entries = parseEntries(config_text);
  return normalizeValidateFreeze(entries, config_text, source_name, options);
}

DerivedRuntimeConfig deriveRuntimeConfig(const FrozenConfig& frozen_config) {
  const SimulationConfig& config = frozen_config.config;

  if (!(config.numerics.t_code_end > config.numerics.t_code_begin)) {
    throw ConfigError(
        "derived runtime config requires numerics.t_code_end > numerics.t_code_begin");
  }
  if (config.cosmology.box_size_x_mpc_comoving <= 0.0 ||
      config.cosmology.box_size_y_mpc_comoving <= 0.0 ||
      config.cosmology.box_size_z_mpc_comoving <= 0.0) {
    throw ConfigError("derived runtime config requires positive comoving box-size axes");
  }
  if (config.numerics.treepm_pm_grid_nx <= 0 ||
      config.numerics.treepm_pm_grid_ny <= 0 ||
      config.numerics.treepm_pm_grid_nz <= 0) {
    throw ConfigError("derived runtime config requires positive TreePM PM grid shape");
  }
  if (!(config.numerics.gravity_softening_kpc_comoving > 0.0)) {
    throw ConfigError("derived runtime config requires positive numerics.gravity_softening");
  }

  const auto choose_softening = [fallback = config.numerics.gravity_softening_kpc_comoving](
                                    double species_value) {
    return species_value > 0.0 ? species_value : fallback;
  };

  DerivedRuntimeConfig derived;
  derived.t_code_begin = config.numerics.t_code_begin;
  derived.t_code_end = config.numerics.t_code_end;
  derived.box_size_mpc_comoving = {
      config.cosmology.box_size_x_mpc_comoving,
      config.cosmology.box_size_y_mpc_comoving,
      config.cosmology.box_size_z_mpc_comoving};
  derived.treepm_pm_grid_shape = {
      config.numerics.treepm_pm_grid_nx,
      config.numerics.treepm_pm_grid_ny,
      config.numerics.treepm_pm_grid_nz};
  derived.gravity_softening_kpc_comoving_by_species = {
      choose_softening(config.numerics.gravity_softening_dark_matter_kpc_comoving),
      choose_softening(config.numerics.gravity_softening_gas_kpc_comoving),
      choose_softening(config.numerics.gravity_softening_star_kpc_comoving),
      choose_softening(config.numerics.gravity_softening_black_hole_kpc_comoving),
      choose_softening(config.numerics.gravity_softening_tracer_kpc_comoving)};
  derived.normalized_config_hash = frozen_config.provenance.config_hash;
  derived.normalized_config_hash_hex = frozen_config.provenance.config_hash_hex;
  return derived;
}

std::string serializeDerivedRuntimeConfig(const DerivedRuntimeConfig& derived_config) {
  std::ostringstream stream;
  stream << std::setprecision(std::numeric_limits<double>::max_digits10);
  stream << "t_code_begin=" << derived_config.t_code_begin << '\n';
  stream << "t_code_end=" << derived_config.t_code_end << '\n';
  stream << "box_size_x_mpc_comoving=" << derived_config.box_size_mpc_comoving[0] << '\n';
  stream << "box_size_y_mpc_comoving=" << derived_config.box_size_mpc_comoving[1] << '\n';
  stream << "box_size_z_mpc_comoving=" << derived_config.box_size_mpc_comoving[2] << '\n';
  stream << "treepm_pm_grid_nx=" << derived_config.treepm_pm_grid_shape[0] << '\n';
  stream << "treepm_pm_grid_ny=" << derived_config.treepm_pm_grid_shape[1] << '\n';
  stream << "treepm_pm_grid_nz=" << derived_config.treepm_pm_grid_shape[2] << '\n';
  stream << "gravity_softening_dark_matter_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[particleSpeciesIndex(ParticleSpecies::kDarkMatter)] << '\n';
  stream << "gravity_softening_gas_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[particleSpeciesIndex(ParticleSpecies::kGas)] << '\n';
  stream << "gravity_softening_star_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[particleSpeciesIndex(ParticleSpecies::kStar)] << '\n';
  stream << "gravity_softening_black_hole_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[particleSpeciesIndex(ParticleSpecies::kBlackHole)] << '\n';
  stream << "gravity_softening_tracer_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[particleSpeciesIndex(ParticleSpecies::kTracer)] << '\n';
  stream << "normalized_config_hash_hex=" << derived_config.normalized_config_hash_hex << '\n';
  return stream.str();
}

void writeNormalizedConfigSnapshot(
    const FrozenConfig& frozen_config,
    const std::filesystem::path& run_directory) {
  std::filesystem::create_directories(run_directory);
  const std::filesystem::path snapshot_path = run_directory / "normalized_config.param.txt";
  std::ofstream output(snapshot_path);
  if (!output) {
    throw ConfigError("failed to write normalized config snapshot: " + snapshot_path.string());
  }
  output << "# normalized_config generated by cosmosim\n";
  output << "# source = " << frozen_config.provenance.source_name << '\n';
  output << frozen_config.normalized_text;
}

std::string modeToString(SimulationMode mode) {
  return modeToLowerString(mode);
}

std::string gravitySolverToString(GravitySolver solver) {
  switch (solver) {
    case GravitySolver::kTreePm:
      return "treepm";
  }
  throw ConfigError("unhandled GravitySolver enum value during serialization");
}

std::string hydroSolverToString(HydroSolver solver) {
  switch (solver) {
    case HydroSolver::kGodunovFv:
      return "godunov_fv";
  }
  throw ConfigError("unhandled HydroSolver enum value during serialization");
}

std::string hydroReconstructionModeToString(HydroReconstructionMode mode) {
  switch (mode) {
    case HydroReconstructionMode::kPiecewiseConstant: return "piecewise_constant";
    case HydroReconstructionMode::kPiecewiseLinear: return "piecewise_linear";
    case HydroReconstructionMode::kPiecewiseParabolic: return "piecewise_parabolic";
  }
  throw ConfigError("unhandled HydroReconstructionMode enum value during serialization");
}

std::string hydroSlopeLimiterModeToString(HydroSlopeLimiterMode mode) {
  switch (mode) {
    case HydroSlopeLimiterMode::kMinmod: return "minmod";
    case HydroSlopeLimiterMode::kMonotonizedCentral: return "mc";
    case HydroSlopeLimiterMode::kVanLeer: return "van_leer";
  }
  throw ConfigError("unhandled HydroSlopeLimiterMode enum value during serialization");
}

std::string hydroIntegrationModeToString(HydroIntegrationMode mode) {
  switch (mode) {
    case HydroIntegrationMode::kMusclHancock: return "muscl_hancock";
    case HydroIntegrationMode::kCornerTransportUpwind: return "ctu";
  }
  throw ConfigError("unhandled HydroIntegrationMode enum value during serialization");
}

std::string hydroRiemannSolverModeToString(HydroRiemannSolverMode mode) {
  switch (mode) {
    case HydroRiemannSolverMode::kHllc: return "hllc";
    case HydroRiemannSolverMode::kHlle: return "hlle";
  }
  throw ConfigError("unhandled HydroRiemannSolverMode enum value during serialization");
}

std::string hydroDualEnergyConfigModeToString(HydroDualEnergyConfigMode mode) {
  switch (mode) {
    case HydroDualEnergyConfigMode::kOff: return "off";
    case HydroDualEnergyConfigMode::kThermalFraction: return "thermal_fraction";
  }
  throw ConfigError("unhandled HydroDualEnergyConfigMode enum value during serialization");
}

std::string sspWindModelConfigToString(SspWindModelConfig mode) {
  switch (mode) {
    case SspWindModelConfig::kNone: return "none";
    case SspWindModelConfig::kTableV1: return "table_v1";
  }
  throw ConfigError("unhandled SspWindModelConfig enum value during serialization");
}

std::string amrTimeSteppingConfigModeToString(AmrTimeSteppingConfigMode mode) {
  switch (mode) {
    case AmrTimeSteppingConfigMode::kSynchronized: return "synchronized";
    case AmrTimeSteppingConfigMode::kRecursiveSubcycling: return "recursive_subcycling";
  }
  throw ConfigError("unhandled AmrTimeSteppingConfigMode enum value during serialization");
}

std::string amrRegridConfigModeToString(AmrRegridConfigMode mode) {
  switch (mode) {
    case AmrRegridConfigMode::kOff: return "off";
    case AmrRegridConfigMode::kAutomatic: return "automatic";
  }
  throw ConfigError("unhandled AmrRegridConfigMode enum value during serialization");
}

std::string nodeMemoryAdmissionModeToString(NodeMemoryAdmissionMode mode) {
  switch (mode) {
    case NodeMemoryAdmissionMode::kOff: return "off";
    case NodeMemoryAdmissionMode::kSharedEqual: return "shared_equal";
    case NodeMemoryAdmissionMode::kSharedPool: return "shared_pool";
  }
  throw ConfigError("unhandled NodeMemoryAdmissionMode enum value during serialization");
}

std::string coordinateFrameToString(CoordinateFrame frame) {
  switch (frame) {
    case CoordinateFrame::kComoving:
      return "comoving";
    case CoordinateFrame::kPhysical:
      return "physical";
  }
  throw ConfigError("unhandled CoordinateFrame enum value during serialization");
}

std::string initialConditionConventionToString(
    InitialConditionConvention convention) {
  switch (convention) {
    case InitialConditionConvention::kGenerated:
      return "generated";
    case InitialConditionConvention::kChuiCanonicalV1:
      return "chui_canonical_v1";
    case InitialConditionConvention::kGadgetArepoBridgeV1:
      return "gadget_arepo_bridge_v1";
    case InitialConditionConvention::kManifestV1:
      return "manifest_v1";
  }
  throw ConfigError(
      "unhandled InitialConditionConvention enum value during serialization");
}

std::string initialConditionCoordinateFrameToString(
    InitialConditionCoordinateFrame frame) {
  switch (frame) {
    case InitialConditionCoordinateFrame::kUnspecified:
      return "unspecified";
    case InitialConditionCoordinateFrame::kComoving:
      return "comoving";
    case InitialConditionCoordinateFrame::kPhysical:
      return "physical";
  }
  throw ConfigError(
      "unhandled InitialConditionCoordinateFrame enum value during serialization");
}

std::string initialConditionVelocityConventionToString(
    InitialConditionVelocityConvention convention) {
  switch (convention) {
    case InitialConditionVelocityConvention::kUnspecified:
      return "unspecified";
    case InitialConditionVelocityConvention::kPhysicalPeculiar:
      return "physical_peculiar";
    case InitialConditionVelocityConvention::kSqrtAScaledPeculiar:
      return "sqrt_a_scaled_peculiar";
    case InitialConditionVelocityConvention::kGadgetArepoStoredPeculiar:
      return "gadget_arepo_stored_peculiar";
    case InitialConditionVelocityConvention::kComovingCoordinateRate:
      return "comoving_coordinate_rate";
  }
  throw ConfigError(
      "unhandled InitialConditionVelocityConvention enum value during serialization");
}

std::string initialConditionSpeciesPolicyToString(
    InitialConditionSpeciesPolicy policy) {
  switch (policy) {
    case InitialConditionSpeciesPolicy::kReject:
      return "reject";
    case InitialConditionSpeciesPolicy::kDarkMatter:
      return "dark_matter";
    case InitialConditionSpeciesPolicy::kStar:
      return "star";
    case InitialConditionSpeciesPolicy::kBlackHole:
      return "black_hole";
    case InitialConditionSpeciesPolicy::kTracer:
      return "tracer";
  }
  throw ConfigError(
      "unhandled InitialConditionSpeciesPolicy enum value during serialization");
}

std::string initialConditionMissingFieldPolicyToString(
    InitialConditionMissingFieldPolicy policy) {
  switch (policy) {
    case InitialConditionMissingFieldPolicy::kReject:
      return "reject";
    case InitialConditionMissingFieldPolicy::kReconstruct:
      return "reconstruct";
    case InitialConditionMissingFieldPolicy::kUseConfigValue:
      return "use_config_value";
    case InitialConditionMissingFieldPolicy::kDialectDefinedDefault:
      return "dialect_defined_default";
    case InitialConditionMissingFieldPolicy::kPreserveUnavailable:
      return "preserve_unavailable";
  }
  throw ConfigError(
      "unhandled InitialConditionMissingFieldPolicy enum value during serialization");
}

std::string modeHydroBoundaryToString(ModeHydroBoundary boundary) {
  switch (boundary) {
    case ModeHydroBoundary::kAuto:
      return "auto";
    case ModeHydroBoundary::kPeriodic:
      return "periodic";
    case ModeHydroBoundary::kOpen:
      return "open";
    case ModeHydroBoundary::kReflective:
      return "reflective";
  }
  throw ConfigError("unhandled ModeHydroBoundary enum value during serialization");
}

std::string modeGravityBoundaryToString(ModeGravityBoundary boundary) {
  switch (boundary) {
    case ModeGravityBoundary::kAuto:
      return "auto";
    case ModeGravityBoundary::kPeriodic:
      return "periodic";
    case ModeGravityBoundary::kIsolatedMonopole:
      return "isolated_monopole";
  }
  throw ConfigError("unhandled ModeGravityBoundary enum value during serialization");
}

std::string starFormationModelKindToString(StarFormationModelKind model) {
  switch (model) {
    case StarFormationModelKind::kLegacySchmidtThreshold:
      return "legacy_schmidt_threshold";
    case StarFormationModelKind::kAdaptiveBoundJeans:
      return "adaptive_bound_jeans";
    case StarFormationModelKind::kEffectiveMultiphaseTngLike:
      return "effective_multiphase_tng_like";
  }
  throw ConfigError("unhandled StarFormationModelKind enum value during serialization");
}

std::string starFormationCollapseTimescaleToString(StarFormationCollapseTimescale mode) {
  switch (mode) {
    case StarFormationCollapseTimescale::kFreeFall:
      return "free_fall";
    case StarFormationCollapseTimescale::kMinimumFreeFallOrCompression:
      return "minimum_free_fall_or_compression";
  }
  throw ConfigError("unhandled StarFormationCollapseTimescale enum value during serialization");
}

std::string starParticleMassPolicyToString(StarParticleMassPolicy policy) {
  switch (policy) {
    case StarParticleMassPolicy::kFixed:
      return "fixed";
    case StarParticleMassPolicy::kGasResolutionFraction:
      return "gas_resolution_fraction";
  }
  throw ConfigError("unhandled StarParticleMassPolicy enum value during serialization");
}

std::string effectiveIsmEosRelaxationToString(EffectiveIsmEosRelaxation mode) {
  switch (mode) {
    case EffectiveIsmEosRelaxation::kInstantaneous: return "instantaneous";
    case EffectiveIsmEosRelaxation::kFiniteTimescale: return "finite_timescale";
  }
  throw ConfigError("unhandled EffectiveIsmEosRelaxation enum value during serialization");
}

std::string effectiveIsmBirthMassConventionToString(EffectiveIsmBirthMassConvention mode) {
  switch (mode) {
    case EffectiveIsmBirthMassConvention::kInitialSspMass: return "initial_ssp_mass";
    case EffectiveIsmBirthMassConvention::kLongLivedMass: return "long_lived_mass";
  }
  throw ConfigError("unhandled EffectiveIsmBirthMassConvention enum value during serialization");
}

std::string effectiveIsmFeedbackCouplingToString(EffectiveIsmFeedbackCoupling mode) {
  switch (mode) {
    case EffectiveIsmFeedbackCoupling::kExternalFeedbackCalibrated: return "external_feedback_calibrated";
    case EffectiveIsmFeedbackCoupling::kEffectiveEosOnly: return "effective_eos_only";
  }
  throw ConfigError("unhandled EffectiveIsmFeedbackCoupling enum value during serialization");
}

std::string feedbackModeToString(FeedbackMode mode) {
  switch (mode) {
    case FeedbackMode::kThermal:
      return "thermal";
    case FeedbackMode::kKinetic:
      return "kinetic";
    case FeedbackMode::kMomentum:
      return "momentum";
    case FeedbackMode::kThermalKineticMomentum:
      return "thermal_kinetic_momentum";
  }
  throw ConfigError("unhandled FeedbackMode enum value during serialization");
}

std::string feedbackVariantToString(FeedbackVariant variant) {
  switch (variant) {
    case FeedbackVariant::kNone:
      return "none";
    case FeedbackVariant::kDelayedCooling:
      return "delayed_cooling";
    case FeedbackVariant::kStochastic:
      return "stochastic";
  }
  throw ConfigError("unhandled FeedbackVariant enum value during serialization");
}


std::string thermochemistryBackendModeToString(ThermochemistryBackendMode value) {
  switch (value) {
    case ThermochemistryBackendMode::kLegacyCooling: return "legacy_cooling";
    case ThermochemistryBackendMode::kAnalyticReference: return "analytic_reference";
    case ThermochemistryBackendMode::kEquilibriumPrimordial: return "equilibrium_primordial";
    case ThermochemistryBackendMode::kFixedHybrid: return "fixed_hybrid";
    case ThermochemistryBackendMode::kFixedDetailed: return "fixed_detailed";
  }
  throw ConfigError("unhandled ThermochemistryBackendMode enum value during serialization");
}

std::string metalCoolingModeConfigToString(MetalCoolingModeConfig value) {
  switch (value) {
    case MetalCoolingModeConfig::kNone: return "none";
    case MetalCoolingModeConfig::kFineStructureCold: return "fine_structure_cold";
    case MetalCoolingModeConfig::kElementTable: return "element_table";
    case MetalCoolingModeConfig::kElementTableAndFineStructure: return "element_table_and_fine_structure";
  }
  throw ConfigError("unhandled MetalCoolingModeConfig enum value during serialization");
}

std::string metalCoolingTableRangeConfigToString(MetalCoolingTableRangeConfig value) {
  switch (value) {
    case MetalCoolingTableRangeConfig::kReject: return "reject";
    case MetalCoolingTableRangeConfig::kClampDensityAndRedshift: return "clamp_density_redshift";
  }
  throw ConfigError("unhandled MetalCoolingTableRangeConfig enum value during serialization");
}

std::string dustModelConfigToString(DustModelConfig value) {
  switch (value) {
    case DustModelConfig::kNone: return "none";
    case DustModelConfig::kTwoReservoirV1: return "two_reservoir_v1";
  }
  throw ConfigError("unhandled DustModelConfig enum value during serialization");
}

std::string ismRegimeModelConfigToString(IsmRegimeModelConfig value) {
  switch (value) {
    case IsmRegimeModelConfig::kLegacy: return "legacy";
    case IsmRegimeModelConfig::kEffectiveOrResolved: return "effective_or_resolved";
  }
  throw ConfigError("unhandled IsmRegimeModelConfig enum value during serialization");
}

std::string sgsTurbulenceModelConfigToString(SgsTurbulenceModelConfig value) {
  switch (value) {
    case SgsTurbulenceModelConfig::kNone: return "none";
    case SgsTurbulenceModelConfig::kBudgetedV1: return "budgeted_v1";
  }
  throw ConfigError("unhandled SgsTurbulenceModelConfig enum value during serialization");
}

std::string chemistryAdaptiveSelectionConfigToString(ChemistryAdaptiveSelectionConfig value) {
  switch (value) {
    case ChemistryAdaptiveSelectionConfig::kOff: return "off";
    case ChemistryAdaptiveSelectionConfig::kLocalV1: return "local_v1";
  }
  throw ConfigError("unhandled ChemistryAdaptiveSelectionConfig enum value during serialization");
}

std::string stellarPopulationModelConfigToString(StellarPopulationModelConfig value) {
  switch (value) {
    case StellarPopulationModelConfig::kLegacyTable: return "legacy_table";
    case StellarPopulationModelConfig::kChabrierReferenceV1: return "chabrier_reference_v1";
    case StellarPopulationModelConfig::kChabrierYieldTableV1: return "chabrier_yield_table_v1";
  }
  throw ConfigError("unhandled StellarPopulationModelConfig enum value during serialization");
}

std::string supernovaFeedbackModelConfigToString(SupernovaFeedbackModelConfig value) {
  switch (value) {
    case SupernovaFeedbackModelConfig::kLegacy: return "legacy";
    case SupernovaFeedbackModelConfig::kMechanicalResolvedUnresolved: return "mechanical_resolved_unresolved";
  }
  throw ConfigError("unhandled SupernovaFeedbackModelConfig enum value during serialization");
}

std::string uvBackgroundModelToString(UvBackgroundModel model) {
  switch (model) {
    case UvBackgroundModel::kNone:
      return "none";
    case UvBackgroundModel::kHm12:
      return "hm12";
    case UvBackgroundModel::kFg20:
      return "fg20";
    case UvBackgroundModel::kTreecoolTable:
      return "treecool_table";
  }
  throw ConfigError("unhandled UvBackgroundModel enum value during serialization");
}

std::string selfShieldingModelToString(SelfShieldingModel model) {
  switch (model) {
    case SelfShieldingModel::kNone:
      return "none";
    case SelfShieldingModel::kRahmati13Like:
      return "rahmati13_like";
  }
  throw ConfigError("unhandled SelfShieldingModel enum value during serialization");
}

std::string coolingModelToString(CoolingModel model) {
  switch (model) {
    case CoolingModel::kPrimordial:
      return "primordial";
    case CoolingModel::kPrimordialMetalLine:
      return "primordial_metal_line";
  }
  throw ConfigError("unhandled CoolingModel enum value during serialization");
}

std::string metalSpeciesModeToString(MetalSpeciesMode mode) {
  switch (mode) {
    case MetalSpeciesMode::kTotalOnly: return "total_only";
    case MetalSpeciesMode::kCoreElements: return "core_elements";
  }
  throw ConfigError("unhandled MetalSpeciesMode enum value during serialization");
}

std::string metalDiffusionModelToString(MetalDiffusionModel model) {
  switch (model) {
    case MetalDiffusionModel::kNone: return "none";
    case MetalDiffusionModel::kSmagorinsky: return "smagorinsky";
  }
  throw ConfigError("unhandled MetalDiffusionModel enum value during serialization");
}

std::string metalDiffusionTimeIntegratorToString(
    MetalDiffusionTimeIntegrator integrator) {
  switch (integrator) {
    case MetalDiffusionTimeIntegrator::kExplicitSubcycling: return "explicit_subcycling";
    case MetalDiffusionTimeIntegrator::kRkl2: return "rkl2";
  }
  throw ConfigError(
      "unhandled MetalDiffusionTimeIntegrator enum value during serialization");
}

std::string integratorTimeVariableToString(IntegratorTimeVariable variable) {
  switch (variable) {
    case IntegratorTimeVariable::kScaleFactor:
      return "scale_factor";
    case IntegratorTimeVariable::kLogScaleFactor:
      return "ln_a";
    case IntegratorTimeVariable::kCodeTime:
      return "t_code";
    case IntegratorTimeVariable::kPhysicalTime:
      return "t_phys";
  }
  throw ConfigError("unhandled IntegratorTimeVariable enum value during serialization");
}

}  // namespace cosmosim::core
