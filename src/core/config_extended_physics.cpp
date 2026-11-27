// CHUI-EXT-PHYS-REF-v1 I02: typed parse, normalized text and admission of the
// extended-physics options, and the S01 single background construction path.
// Source-only pass 01: NOT compiled or executed by the implementing agent.

#include "core/internal/config_extended_physics.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <string_view>

#include "cosmosim/core/simulation_mode.hpp"

namespace cosmosim::core {

namespace {

template <typename E>
struct ExtendedChoice {
  std::string_view name;
  E value;
};

[[nodiscard]] std::string lowerTrimmed(const std::string& raw) {
  std::size_t begin = 0;
  std::size_t end = raw.size();
  while (begin < end && std::isspace(static_cast<unsigned char>(raw[begin]))) ++begin;
  while (end > begin && std::isspace(static_cast<unsigned char>(raw[end - 1]))) --end;
  std::string out = raw.substr(begin, end - begin);
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return out;
}

[[nodiscard]] std::string trimmed(const std::string& raw) {
  std::size_t begin = 0;
  std::size_t end = raw.size();
  while (begin < end && std::isspace(static_cast<unsigned char>(raw[begin]))) ++begin;
  while (end > begin && std::isspace(static_cast<unsigned char>(raw[end - 1]))) --end;
  return raw.substr(begin, end - begin);
}

template <typename E, std::size_t N>
[[nodiscard]] E parseChoice(const std::string& raw, const std::string& key,
                            const std::array<ExtendedChoice<E>, N>& table) {
  const std::string lower = lowerTrimmed(raw);
  for (const auto& entry : table) {
    if (entry.name == lower) return entry.value;
  }
  std::string allowed;
  for (const auto& entry : table) {
    if (!allowed.empty()) allowed += ", ";
    allowed += std::string(entry.name);
  }
  throw ConfigError("key '" + key + "': invalid value '" + raw + "' (allowed: " + allowed + ")");
}

template <typename E, std::size_t N>
[[nodiscard]] std::string choiceToString(E value, const std::array<ExtendedChoice<E>, N>& table,
                                         const char* type_name) {
  for (const auto& entry : table) {
    if (entry.value == value) return std::string(entry.name);
  }
  throw ConfigError(std::string("unhandled ") + type_name + " enum value during serialization");
}

[[nodiscard]] double parseFinite(const std::string& raw, const std::string& key) {
  const std::string text = trimmed(raw);
  if (text.empty()) throw ConfigError("key '" + key + "': empty floating-point value");
  char* end = nullptr;
  const double value = std::strtod(text.c_str(), &end);
  if (end != text.c_str() + static_cast<std::ptrdiff_t>(text.size()) || !std::isfinite(value)) {
    throw ConfigError("key '" + key + "': invalid finite floating-point value '" + raw + "'");
  }
  return value;
}

[[nodiscard]] bool parseBoolean(const std::string& raw, const std::string& key) {
  const std::string lower = lowerTrimmed(raw);
  if (lower == "true") return true;
  if (lower == "false") return false;
  throw ConfigError("key '" + key + "': invalid boolean value '" + raw + "'");
}

[[nodiscard]] std::uint32_t parseUint32(const std::string& raw, const std::string& key) {
  const std::string text = trimmed(raw);
  if (text.empty() || !std::all_of(text.begin(), text.end(), [](unsigned char c) { return std::isdigit(c); })) {
    throw ConfigError("key '" + key + "': invalid unsigned integer '" + raw + "'");
  }
  const unsigned long long value = std::strtoull(text.c_str(), nullptr, 10);
  if (value > std::numeric_limits<std::uint32_t>::max()) {
    throw ConfigError("key '" + key + "': value exceeds 32-bit range");
  }
  return static_cast<std::uint32_t>(value);
}

[[nodiscard]] std::uint64_t parseUint64(const std::string& raw, const std::string& key) {
  const std::string text = trimmed(raw);
  if (text.empty() || text.size() > 20U ||
      !std::all_of(text.begin(), text.end(), [](unsigned char c) { return std::isdigit(c); })) {
    throw ConfigError("key '" + key + "': invalid unsigned integer '" + raw + "'");
  }
  errno = 0;
  const unsigned long long value = std::strtoull(text.c_str(), nullptr, 10);
  if (errno == ERANGE) throw ConfigError("key '" + key + "': value exceeds 64-bit range");
  return static_cast<std::uint64_t>(value);
}

// Comma-separated list; an empty value is the empty list.
[[nodiscard]] std::vector<std::string> parseStringList(const std::string& raw) {
  std::vector<std::string> out;
  const std::string text = trimmed(raw);
  if (text.empty()) return out;
  std::size_t start = 0;
  while (start <= text.size()) {
    const std::size_t comma = text.find(',', start);
    const std::string item = trimmed(text.substr(start, comma == std::string::npos ? std::string::npos
                                                                                   : comma - start));
    out.push_back(item);
    if (comma == std::string::npos) break;
    start = comma + 1U;
  }
  return out;
}

[[nodiscard]] std::vector<double> parseDoubleList(const std::string& raw, const std::string& key) {
  std::vector<double> out;
  for (const std::string& item : parseStringList(raw)) out.push_back(parseFinite(item, key));
  return out;
}

void require(bool condition, const std::string& message) {
  if (!condition) throw ConfigError(message);
}

[[nodiscard]] bool isHex64(const std::string& value) {
  return value.size() == 64U &&
         std::all_of(value.begin(), value.end(), [](unsigned char c) { return std::isxdigit(c) != 0; });
}

template <typename T>
void writeList(std::ostream& stream, const std::vector<T>& values) {
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0U) stream << ", ";
    stream << values[i];
  }
}

// ----------------------------------------------------------- spelling tables
inline constexpr std::array<ExtendedChoice<NeutrinoMethodConfig>, 6> k_neutrino_method_names{{
    {"none", NeutrinoMethodConfig::kNone},
    {"background_only", NeutrinoMethodConfig::kBackgroundOnly},
    {"transfer_ratio", NeutrinoMethodConfig::kTransferRatio},
    {"linear_response", NeutrinoMethodConfig::kLinearResponse},
    {"particles", NeutrinoMethodConfig::kParticles},
    {"delta_f", NeutrinoMethodConfig::kDeltaF},
}};
inline constexpr std::array<ExtendedChoice<BoltzmannTableFormatConfig>, 3> k_boltzmann_table_format_names{{
    {"chui_v1", BoltzmannTableFormatConfig::kChuiV1},
    {"class_tk", BoltzmannTableFormatConfig::kClassTk},
    {"camb_transfer", BoltzmannTableFormatConfig::kCambTransfer},
}};
inline constexpr std::array<ExtendedChoice<BoltzmannGaugeConfig>, 2> k_boltzmann_gauge_names{{
    {"synchronous", BoltzmannGaugeConfig::kSynchronous},
    {"newtonian", BoltzmannGaugeConfig::kNewtonian},
}};

constexpr const char* k_instantaneous_decoupling_token = "instantaneous_decoupling";

}  // namespace

std::string neutrinoMethodConfigToString(NeutrinoMethodConfig value) {
  return choiceToString(value, k_neutrino_method_names, "NeutrinoMethodConfig");
}
std::string boltzmannTableFormatConfigToString(BoltzmannTableFormatConfig value) {
  return choiceToString(value, k_boltzmann_table_format_names, "BoltzmannTableFormatConfig");
}
std::string boltzmannGaugeConfigToString(BoltzmannGaugeConfig value) {
  return choiceToString(value, k_boltzmann_gauge_names, "BoltzmannGaugeConfig");
}

// --------------------------------------------------------------- S01 factory
CosmologyBackgroundConfig cosmologyBackgroundConfigFromSimulation(const SimulationConfig& config) {
  CosmologyBackgroundConfig background;
  background.hubble_param = config.cosmology.hubble_param;
  background.omega_matter = config.cosmology.omega_matter;
  background.omega_lambda = config.cosmology.omega_lambda;
  // No config key feeds omega_radiation or omega_curvature (legacy flat
  // matter+Lambda); relativistic species come only from radiation_species.
  background.omega_radiation = 0.0;
  background.omega_curvature = 0.0;
  background.radiation_species = config.cosmology.radiation_species;
  return background;
}

namespace internal {

const std::vector<ExtendedKeyDefault>& extendedPhysicsKeyDefaults() {
  static const std::vector<ExtendedKeyDefault> keys = {
      // S01/N01 background species ([cosmology]).
      {"cosmology.radiation_background", "none"},
      {"cosmology.t_cmb0_k", "0"},
      {"cosmology.n_ur", "0"},
      {"cosmology.ncdm_masses_ev", ""},
      {"cosmology.ncdm_temperature_ratio", "0"},
      // N02-N05 ([neutrinos]).
      {"neutrinos.method", "none"},
      {"neutrinos.table_format", "chui_v1"},
      {"neutrinos.table_gauge", "synchronous"},
      {"neutrinos.table_paths", ""},
      {"neutrinos.table_sha256", ""},
      {"neutrinos.table_redshifts", ""},
      {"neutrinos.response_k_bins", "64"},
      {"neutrinos.response_history_capacity", "4096"},
      {"neutrinos.response_max_mean_velocity_c", "0"},
      {"neutrinos.response_k_min_h_mpc", "0"},
      {"neutrinos.response_k_max_h_mpc", "0"},
  };
  static const std::vector<ExtendedKeyDefault> all_keys = [] {
    std::vector<ExtendedKeyDefault> merged = keys;
    const std::vector<ExtendedKeyDefault>& sections = extendedPhysicsSectionKeyDefaults();
    merged.insert(merged.end(), sections.begin(), sections.end());
    return merged;
  }();
  return all_keys;
}

void parseExtendedPhysicsConfig(const ExtendedKeyLookup& lookup, SimulationConfig& config) {
  // ---------------------------------------------------------------- [cosmology]
  {
    const std::string mode = lowerTrimmed(lookup("cosmology.radiation_background"));
    if (mode == "none") {
      config.cosmology.radiation_species.enabled = false;
    } else if (mode == "photons_neutrinos") {
      config.cosmology.radiation_species.enabled = true;
    } else {
      throw ConfigError("key 'cosmology.radiation_background': invalid value '" + mode +
                        "' (allowed: none, photons_neutrinos)");
    }
    auto& species = config.cosmology.radiation_species;
    species.t_cmb0_k = parseFinite(lookup("cosmology.t_cmb0_k"), "cosmology.t_cmb0_k");
    species.n_ur = parseFinite(lookup("cosmology.n_ur"), "cosmology.n_ur");
    species.ncdm_masses_ev = parseDoubleList(lookup("cosmology.ncdm_masses_ev"), "cosmology.ncdm_masses_ev");
    const std::string ratio = lowerTrimmed(lookup("cosmology.ncdm_temperature_ratio"));
    if (ratio == k_instantaneous_decoupling_token) {
      species.ncdm_temperature_ratio = instantaneousDecouplingNeutrinoTemperatureRatio();
      config.cosmology.ncdm_temperature_ratio_instantaneous_decoupling = true;
    } else {
      species.ncdm_temperature_ratio = parseFinite(ratio, "cosmology.ncdm_temperature_ratio");
      config.cosmology.ncdm_temperature_ratio_instantaneous_decoupling = false;
    }
  }
  // ---------------------------------------------------------------- [neutrinos]
  {
    auto& nu = config.extended.neutrinos;
    nu.method = parseChoice(lookup("neutrinos.method"), "neutrinos.method", k_neutrino_method_names);
    nu.table_format =
        parseChoice(lookup("neutrinos.table_format"), "neutrinos.table_format", k_boltzmann_table_format_names);
    nu.table_gauge = parseChoice(lookup("neutrinos.table_gauge"), "neutrinos.table_gauge", k_boltzmann_gauge_names);
    nu.table_paths = parseStringList(lookup("neutrinos.table_paths"));
    nu.table_sha256 = parseStringList(lookup("neutrinos.table_sha256"));
    nu.table_redshifts = parseDoubleList(lookup("neutrinos.table_redshifts"), "neutrinos.table_redshifts");
    nu.response_k_bins = parseUint32(lookup("neutrinos.response_k_bins"), "neutrinos.response_k_bins");
    nu.response_history_capacity =
        parseUint32(lookup("neutrinos.response_history_capacity"), "neutrinos.response_history_capacity");
    nu.response_max_mean_velocity_c =
        parseFinite(lookup("neutrinos.response_max_mean_velocity_c"), "neutrinos.response_max_mean_velocity_c");
    nu.response_k_min_h_mpc = parseFinite(lookup("neutrinos.response_k_min_h_mpc"), "neutrinos.response_k_min_h_mpc");
    nu.response_k_max_h_mpc = parseFinite(lookup("neutrinos.response_k_max_h_mpc"), "neutrinos.response_k_max_h_mpc");
  }
  parseExtendedPhysicsSectionsAfterNeutrinos(lookup, config);
}

void appendCosmologyRadiationNormalizedText(std::ostream& stream, const CosmologyConfig& cosmology) {
  const RadiationSpeciesConfig& species = cosmology.radiation_species;
  if (!species.enabled) return;
  stream << "radiation_background = photons_neutrinos\n";
  stream << "t_cmb0_k = " << species.t_cmb0_k << '\n';
  stream << "n_ur = " << species.n_ur << '\n';
  stream << "ncdm_masses_ev = ";
  writeList(stream, species.ncdm_masses_ev);
  stream << '\n';
  if (cosmology.ncdm_temperature_ratio_instantaneous_decoupling) {
    stream << "ncdm_temperature_ratio = " << k_instantaneous_decoupling_token << '\n';
  } else {
    stream << "ncdm_temperature_ratio = " << species.ncdm_temperature_ratio << '\n';
  }
}

void appendExtendedPhysicsNormalizedText(std::ostream& stream, const ExtendedPhysicsConfig& extended) {
  const NeutrinoPerturbationConfig& nu = extended.neutrinos;
  if (nu.method != NeutrinoMethodConfig::kNone) {
    stream << "\n[neutrinos]\n";
    stream << "method = " << neutrinoMethodConfigToString(nu.method) << '\n';
    stream << "table_format = " << boltzmannTableFormatConfigToString(nu.table_format) << '\n';
    stream << "table_gauge = " << boltzmannGaugeConfigToString(nu.table_gauge) << '\n';
    stream << "table_paths = ";
    writeList(stream, nu.table_paths);
    stream << "\ntable_sha256 = ";
    writeList(stream, nu.table_sha256);
    stream << "\ntable_redshifts = ";
    writeList(stream, nu.table_redshifts);
    stream << '\n';
    stream << "response_k_bins = " << nu.response_k_bins << '\n';
    stream << "response_history_capacity = " << nu.response_history_capacity << '\n';
    stream << "response_max_mean_velocity_c = " << nu.response_max_mean_velocity_c << '\n';
    stream << "response_k_min_h_mpc = " << nu.response_k_min_h_mpc << '\n';
    stream << "response_k_max_h_mpc = " << nu.response_k_max_h_mpc << '\n';
  }
  appendExtendedPhysicsSectionsAfterNeutrinos(stream, extended);
}

}  // namespace internal

// ------------------------------------------------------------------ admission
namespace {

[[nodiscard]] bool isCosmologicalRun(const SimulationConfig& config) {
  return buildModePolicy(config.mode).cosmological_comoving_frame;
}

void validateRadiationBackground(const SimulationConfig& config) {
  const CosmologyConfig& cosmology = config.cosmology;
  const RadiationSpeciesConfig& species = cosmology.radiation_species;
  if (!species.enabled) {
    require(species == RadiationSpeciesConfig{} && !cosmology.ncdm_temperature_ratio_instantaneous_decoupling,
            "cosmology.radiation_background=none but cosmology.t_cmb0_k / n_ur / ncdm_masses_ev / "
            "ncdm_temperature_ratio differ from their defaults; those keys would have no effect. Remove them or "
            "set cosmology.radiation_background=photons_neutrinos");
    return;
  }
  require(isCosmologicalRun(config),
          "cosmology.radiation_background=photons_neutrinos requires a cosmological comoving mode (the "
          "relativistic species enter H(a) only)");
  // Strict closure, positivity and Omega_cb > 0 are enforced by the same
  // constructor the runtime uses (single authority); its message is kept.
  try {
    const LambdaCdmBackground background(cosmologyBackgroundConfigFromSimulation(config));
    require(config.cosmology.omega_baryon <= background.coldMatterDensityParameter0(),
            "cosmology.omega_baryon must be <= Omega_cb = omega_matter - Omega_ncdm,0 (massive neutrinos are part "
            "of omega_matter)");
    // Unsupported early regime: the Newtonian TreePM integrator evolves no
    // relativistic perturbations, so a radiation-dominated start is refused.
    const double a_eq = background.matterRadiationEqualityScaleFactor();
    require(config.numerics.a_begin > a_eq,
            "numerics.a_begin = " + std::to_string(config.numerics.a_begin) +
                " precedes matter-radiation equality a_eq = " + std::to_string(a_eq) +
                " of the configured background; radiation-dominated starts are outside the Newtonian "
                "TreePM domain (no relativistic perturbation evolution)");
  } catch (const ConfigError&) {
    throw;
  } catch (const std::exception& error) {
    throw ConfigError(std::string("cosmology radiation background rejected: ") + error.what());
  }
}

void validateNeutrinos(const SimulationConfig& config) {
  const NeutrinoPerturbationConfig& nu = config.extended.neutrinos;
  const RadiationSpeciesConfig& species = config.cosmology.radiation_species;
  const bool has_masses = species.enabled && !species.ncdm_masses_ev.empty();
  if (nu.method == NeutrinoMethodConfig::kNone) {
    NeutrinoPerturbationConfig defaults;
    require(nu == defaults,
            "neutrinos.method=none but other [neutrinos] keys differ from their defaults; those keys would have no "
            "effect. Remove them or select a neutrinos.method");
    require(!has_masses,
            "massive neutrino eigenstates are configured (cosmology.ncdm_masses_ev) but neutrinos.method=none: "
            "declare the perturbation treatment explicitly (background_only, transfer_ratio or linear_response); "
            "neglected perturbations are never implied");
    return;
  }
  require(has_masses,
          "neutrinos.method=" + neutrinoMethodConfigToString(nu.method) +
              " requires cosmology.radiation_background=photons_neutrinos with at least one "
              "cosmology.ncdm_masses_ev entry");
  if (nu.method == NeutrinoMethodConfig::kParticles || nu.method == NeutrinoMethodConfig::kDeltaF) {
    require(false,
            "neutrinos.method=" + neutrinoMethodConfigToString(nu.method) +
                ": the sampling, relativistic drift/kick and delta-f weight kernels exist "
                "(cosmosim/neutrinos/neutrino_particles.hpp) but the canonical particle model has no neutrino "
                "species (core::ParticleSpecies, k_particle_species_count), so species-resolved softening, "
                "migration, restart and the signed-weight PM deposition edge are not wired; a signed weight would "
                "otherwise pass through the positive-mass tree (FAIL_CLOSED)");
  }
  if (nu.method == NeutrinoMethodConfig::kBackgroundOnly) {
    require(nu.table_paths.empty() && nu.table_sha256.empty() && nu.table_redshifts.empty(),
            "neutrinos.method=background_only takes no Boltzmann tables (inert keys)");
    NeutrinoPerturbationConfig defaults;
    defaults.method = nu.method;
    require(nu == defaults, "neutrinos.method=background_only: response_* / table_* keys would have no effect");
    return;
  }
  // transfer_ratio / linear_response: versioned Boltzmann tables and a PM hook.
  require(!nu.table_paths.empty(),
          "neutrinos.method=" + neutrinoMethodConfigToString(nu.method) +
              " requires neutrinos.table_paths (no Boltzmann tables are bundled; no native Einstein-Boltzmann "
              "engine exists)");
  require(nu.table_paths.size() == nu.table_sha256.size() && nu.table_paths.size() == nu.table_redshifts.size(),
          "neutrinos.table_paths, table_sha256 and table_redshifts must have the same number of entries");
  for (const std::string& digest : nu.table_sha256) {
    require(isHex64(digest), "neutrinos.table_sha256 entries must be 64 hexadecimal characters");
  }
  for (std::size_t i = 0; i < nu.table_redshifts.size(); ++i) {
    require(nu.table_redshifts[i] >= 0.0, "neutrinos.table_redshifts must be >= 0");
    if (i > 0U) {
      require(nu.table_redshifts[i] < nu.table_redshifts[i - 1U],
              "neutrinos.table_redshifts must be strictly decreasing (ascending scale factor)");
    }
  }
  require(nu.table_redshifts.size() >= 2U,
          "neutrinos.method=" + neutrinoMethodConfigToString(nu.method) +
              " needs at least two table epochs: the ratio/initial state is interpolated in a and never "
              "extrapolated");
  const double z_begin = 1.0 / config.numerics.a_begin - 1.0;
  const double z_end = 1.0 / config.numerics.a_end - 1.0;
  if (nu.method == NeutrinoMethodConfig::kTransferRatio) {
    require(nu.table_redshifts.front() >= z_begin && nu.table_redshifts.back() <= z_end,
            "neutrinos.method=transfer_ratio: the table redshifts must bracket [z_end, z_begin] (no extrapolation "
            "in a)");
  } else {
    require(nu.table_redshifts.front() >= z_begin && nu.table_redshifts.back() <= z_begin,
            "neutrinos.method=linear_response: the table epochs must bracket z_begin (initial delta_nu, theta_nu "
            "and delta_cb are interpolated at the start, never extrapolated)");
    require(std::isfinite(nu.response_max_mean_velocity_c) && nu.response_max_mean_velocity_c > 0.0 &&
                nu.response_max_mean_velocity_c < 1.0,
            "neutrinos.response_max_mean_velocity_c must be set in (0, 1): the Newtonian nonrelativistic "
            "free-streaming kernel is admitted only below this <|p|>/(m c) at a_begin (no default validity limit)");
    require(nu.response_history_capacity >= 16U && nu.response_history_capacity <= 1048576U,
            "neutrinos.response_history_capacity must lie in [16, 1048576]");
    try {
      const LambdaCdmBackground background(cosmologyBackgroundConfigFromSimulation(config));
      for (std::size_t i = 0; i < background.massiveNeutrinoCount(); ++i) {
        const double beta = background.massiveNeutrinoMeanMomentumOverMass(i, config.numerics.a_begin);
        require(beta < nu.response_max_mean_velocity_c,
                "neutrinos.method=linear_response: eigenstate " + std::to_string(i) +
                    " has <|p|>/(m c) = " + std::to_string(beta) +
                    " at a_begin, above neutrinos.response_max_mean_velocity_c; the nonrelativistic kernel is "
                    "outside its validity domain (start later or use transfer_ratio)");
      }
    } catch (const ConfigError&) {
      throw;
    } catch (const std::exception& error) {
      throw ConfigError(std::string("neutrinos.method=linear_response background check failed: ") + error.what());
    }
  }
  require(nu.response_k_bins >= 8U && nu.response_k_bins <= 4096U,
          "neutrinos.response_k_bins must lie in [8, 4096]");
  require(std::isfinite(nu.response_k_min_h_mpc) && nu.response_k_min_h_mpc > 0.0 &&
              std::isfinite(nu.response_k_max_h_mpc) && nu.response_k_max_h_mpc > nu.response_k_min_h_mpc,
          "neutrinos.response_k_min_h_mpc and response_k_max_h_mpc must satisfy 0 < k_min < k_max (h/Mpc); the "
          "PM solve refuses wavenumbers outside this domain instead of extrapolating");
  require(config.parallel.gpu_devices <= 0,
          "neutrinos.method=" + neutrinoMethodConfigToString(nu.method) +
              ": the PM spectral source hook exists only on the host periodic Poisson operator; CUDA PM "
              "(parallel.gpu_devices > 0) is refused");
  if (nu.method == NeutrinoMethodConfig::kLinearResponse) {
    require(nu.table_format != BoltzmannTableFormatConfig::kCambTransfer,
            "neutrinos.method=linear_response needs massive-neutrino velocity divergences for the initial "
            "free-streaming term; CAMB transfer output has none (use chui_v1 or Newtonian-gauge class_tk)");
    require(nu.table_gauge == BoltzmannGaugeConfig::kNewtonian,
            "neutrinos.method=linear_response requires Newtonian-gauge tables (the response equation and its "
            "initial theta are Newtonian-gauge quantities)");
    for (const double m : species.ncdm_masses_ev) {
      require(std::abs(m - species.ncdm_masses_ev.front()) <= 1.0e-12 * species.ncdm_masses_ev.front(),
              "neutrinos.method=linear_response supports degenerate eigenstate masses only (one Fermi-Dirac "
              "velocity scale); use transfer_ratio for non-degenerate hierarchies");
    }
  }
  require(buildModePolicy(config.mode).gravity_boundary == GravityBoundaryModel::kPeriodicPoisson,
          "neutrinos.method=" + neutrinoMethodConfigToString(nu.method) +
              " modifies the periodic PM Poisson operator; isolated/open gravity is not supported");
}

}  // namespace

void validateExtendedPhysicsConfig(const SimulationConfig& config) {
  validateRadiationBackground(config);
  validateNeutrinos(config);
  internal::validateExtendedPhysicsSectionsAfterNeutrinos(config);
}

}  // namespace cosmosim::core
