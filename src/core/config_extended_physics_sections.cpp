// CHUI-EXT-PHYS-REF-v1 I02: [cosmic_rays], [observables] and [instruments]
// typed options. Source-only pass 01: NOT compiled or executed.

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <string>
#include <string_view>

#include "cosmosim/core/simulation_mode.hpp"

#include "core/internal/config_extended_physics.hpp"

namespace cosmosim::core {
namespace {

template <typename E>
struct SectionChoice {
  std::string_view name;
  E value;
};

[[nodiscard]] std::string lowerTrimmedValue(const std::string& raw) {
  std::size_t b = 0;
  std::size_t e = raw.size();
  while (b < e && std::isspace(static_cast<unsigned char>(raw[b]))) ++b;
  while (e > b && std::isspace(static_cast<unsigned char>(raw[e - 1]))) --e;
  std::string out = raw.substr(b, e - b);
  std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return out;
}

template <typename E, std::size_t N>
[[nodiscard]] E parseSectionChoice(const std::string& raw, const std::string& key,
                                   const std::array<SectionChoice<E>, N>& table) {
  const std::string lower = lowerTrimmedValue(raw);
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
[[nodiscard]] std::string sectionChoiceToString(E value, const std::array<SectionChoice<E>, N>& table,
                                                const char* type_name) {
  for (const auto& entry : table) {
    if (entry.value == value) return std::string(entry.name);
  }
  throw ConfigError(std::string("unhandled ") + type_name + " enum value during serialization");
}

void requireSection(bool condition, const std::string& message) {
  if (!condition) throw ConfigError(message);
}

inline constexpr std::array<SectionChoice<CosmicRayModelConfig>, 4> k_cosmic_ray_model_names{{
    {"none", CosmicRayModelConfig::kNone},
    {"grey", CosmicRayModelConfig::kGrey},
    {"two_moment", CosmicRayModelConfig::kTwoMoment},
    {"spectral", CosmicRayModelConfig::kSpectral},
}};

inline constexpr std::array<SectionChoice<LosAxisConfig>, 3> k_los_axis_names{{
    {"x", LosAxisConfig::kX}, {"y", LosAxisConfig::kY}, {"z", LosAxisConfig::kZ}}};
inline constexpr std::array<SectionChoice<ForestProfileConfig>, 2> k_forest_profile_names{{
    {"voigt", ForestProfileConfig::kVoigt}, {"doppler_only", ForestProfileConfig::kDopplerOnly}}};
inline constexpr std::array<SectionChoice<HiSpinModelConfig>, 2> k_hi_spin_model_names{{
    {"saturated_kinetic", HiSpinModelConfig::kSaturatedKinetic},
    {"collisional_table", HiSpinModelConfig::kCollisionalTable}}};
inline constexpr std::array<SectionChoice<HiTransferConfig>, 2> k_hi_transfer_names{{
    {"full_optical_depth", HiTransferConfig::kFullOpticalDepth}, {"optically_thin", HiTransferConfig::kOpticallyThin}}};

inline constexpr std::array<SectionChoice<InstrumentBackendConfig>, 3> k_instrument_backend_names{{
    {"direct_rime", InstrumentBackendConfig::kDirectRime},
    {"oskar_export", InstrumentBackendConfig::kOskarExport},
    {"accelerated", InstrumentBackendConfig::kAccelerated}}};
inline constexpr std::array<SectionChoice<InstrumentBeamConfig>, 3> k_instrument_beam_names{{
    {"none", InstrumentBeamConfig::kNone}, {"gaussian", InstrumentBeamConfig::kGaussian}, {"airy", InstrumentBeamConfig::kAiry}}};

[[nodiscard]] bool parseSectionBool(const std::string& raw, const std::string& key) {
  const std::string v = lowerTrimmedValue(raw);
  if (v == "true") return true;
  if (v == "false") return false;
  throw ConfigError("key '" + key + "': invalid boolean value '" + raw + "'");
}

[[nodiscard]] std::uint32_t parseSectionUint(const std::string& raw, const std::string& key) {
  const std::string v = lowerTrimmedValue(raw);
  if (v.empty() || v.size() > 10U || !std::all_of(v.begin(), v.end(), [](unsigned char c) { return std::isdigit(c); })) {
    throw ConfigError("key '" + key + "': invalid unsigned integer '" + raw + "'");
  }
  const unsigned long long value = std::stoull(v);
  if (value > 0xFFFFFFFFULL) throw ConfigError("key '" + key + "': value exceeds 32-bit range");
  return static_cast<std::uint32_t>(value);
}

[[nodiscard]] double parseSectionDouble(const std::string& raw, const std::string& key) {
  const std::string v = lowerTrimmedValue(raw);
  char* end = nullptr;
  const double value = std::strtod(v.c_str(), &end);
  if (v.empty() || end != v.c_str() + static_cast<std::ptrdiff_t>(v.size()) || !std::isfinite(value)) {
    throw ConfigError("key '" + key + "': invalid finite floating-point value '" + raw + "'");
  }
  return value;
}

[[nodiscard]] std::string trimmedValue(const std::string& raw) {
  std::size_t b = 0;
  std::size_t e = raw.size();
  while (b < e && std::isspace(static_cast<unsigned char>(raw[b]))) ++b;
  while (e > b && std::isspace(static_cast<unsigned char>(raw[e - 1]))) --e;
  return raw.substr(b, e - b);
}

}  // namespace

std::string losAxisConfigToString(LosAxisConfig value) { return sectionChoiceToString(value, k_los_axis_names, "LosAxisConfig"); }
std::string forestProfileConfigToString(ForestProfileConfig value) {
  return sectionChoiceToString(value, k_forest_profile_names, "ForestProfileConfig");
}
std::string hiSpinModelConfigToString(HiSpinModelConfig value) {
  return sectionChoiceToString(value, k_hi_spin_model_names, "HiSpinModelConfig");
}
std::string hiTransferConfigToString(HiTransferConfig value) {
  return sectionChoiceToString(value, k_hi_transfer_names, "HiTransferConfig");
}

std::string instrumentBackendConfigToString(InstrumentBackendConfig value) {
  return sectionChoiceToString(value, k_instrument_backend_names, "InstrumentBackendConfig");
}
std::string instrumentBeamConfigToString(InstrumentBeamConfig value) {
  return sectionChoiceToString(value, k_instrument_beam_names, "InstrumentBeamConfig");
}

std::string cosmicRayModelConfigToString(CosmicRayModelConfig value) {
  return sectionChoiceToString(value, k_cosmic_ray_model_names, "CosmicRayModelConfig");
}

namespace internal {

const std::vector<ExtendedKeyDefault>& extendedPhysicsSectionKeyDefaults() {
  static const std::vector<ExtendedKeyDefault> keys = {
      {"cosmic_rays.model", "none"},
      {"observables.enabled", "false"},
      {"observables.interval_steps", "0"},
      {"observables.los_axis", "z"},
      {"observables.sightlines_per_side", "0"},
      {"observables.output_stem", "observables"},
      {"observables.lyman_alpha_forest", "false"},
      {"observables.forest_profile", "voigt"},
      {"observables.forest_pixels", "0"},
      {"observables.forest_pdf_bins", "20"},
      {"observables.hi_cube", "false"},
      {"observables.hi_spin_model", "saturated_kinetic"},
      {"observables.hi_transfer", "full_optical_depth"},
      {"observables.hi_channels", "0"},
      {"observables.hi_velocity_half_width_km_s", "0"},
      {"observables.hi_collision_table_path", ""},
      {"observables.hi_collision_table_sha256", ""},
      {"observables.sz_maps", "false"},
      {"observables.dispersion_map", "false"},
      {"instruments.enabled", "false"},
      {"instruments.backend", "direct_rime"},
      {"instruments.station_layout_path", ""},
      {"instruments.station_layout_sha256", ""},
      {"instruments.latitude_deg", "0"},
      {"instruments.longitude_deg", "0"},
      {"instruments.phase_centre_ra_deg", "0"},
      {"instruments.phase_centre_dec_deg", "0"},
      {"instruments.hour_angle_start_deg", "0"},
      {"instruments.integration_s", "0"},
      {"instruments.time_steps", "0"},
      {"instruments.time_subsamples", "1"},
      {"instruments.channel_subsamples", "1"},
      {"instruments.sefd_jy", "0"},
      {"instruments.noise_seed", "1"},
      {"instruments.beam", "none"},
      {"instruments.beam_fwhm_deg", "0"},
      {"instruments.beam_reference_mhz", "0"},
      {"instruments.dish_diameter_m", "0"},
      {"instruments.source_distance_mpc", "0"},
  };
  return keys;
}

void parseExtendedPhysicsSectionsAfterNeutrinos(const ExtendedKeyLookup& lookup, SimulationConfig& config) {
  config.extended.cosmic_rays.model =
      parseSectionChoice(lookup("cosmic_rays.model"), "cosmic_rays.model", k_cosmic_ray_model_names);
  ObservablesConfig& o = config.extended.observables;
  o.enabled = parseSectionBool(lookup("observables.enabled"), "observables.enabled");
  o.interval_steps = parseSectionUint(lookup("observables.interval_steps"), "observables.interval_steps");
  o.los_axis = parseSectionChoice(lookup("observables.los_axis"), "observables.los_axis", k_los_axis_names);
  o.sightlines_per_side = parseSectionUint(lookup("observables.sightlines_per_side"), "observables.sightlines_per_side");
  o.output_stem = trimmedValue(lookup("observables.output_stem"));
  o.lyman_alpha_forest = parseSectionBool(lookup("observables.lyman_alpha_forest"), "observables.lyman_alpha_forest");
  o.forest_profile =
      parseSectionChoice(lookup("observables.forest_profile"), "observables.forest_profile", k_forest_profile_names);
  o.forest_pixels = parseSectionUint(lookup("observables.forest_pixels"), "observables.forest_pixels");
  o.forest_pdf_bins = parseSectionUint(lookup("observables.forest_pdf_bins"), "observables.forest_pdf_bins");
  o.hi_cube = parseSectionBool(lookup("observables.hi_cube"), "observables.hi_cube");
  o.hi_spin_model =
      parseSectionChoice(lookup("observables.hi_spin_model"), "observables.hi_spin_model", k_hi_spin_model_names);
  o.hi_transfer = parseSectionChoice(lookup("observables.hi_transfer"), "observables.hi_transfer", k_hi_transfer_names);
  o.hi_channels = parseSectionUint(lookup("observables.hi_channels"), "observables.hi_channels");
  o.hi_velocity_half_width_km_s =
      parseSectionDouble(lookup("observables.hi_velocity_half_width_km_s"), "observables.hi_velocity_half_width_km_s");
  o.hi_collision_table_path = trimmedValue(lookup("observables.hi_collision_table_path"));
  o.hi_collision_table_sha256 = trimmedValue(lookup("observables.hi_collision_table_sha256"));
  o.sz_maps = parseSectionBool(lookup("observables.sz_maps"), "observables.sz_maps");
  o.dispersion_map = parseSectionBool(lookup("observables.dispersion_map"), "observables.dispersion_map");
  InstrumentsConfig& in = config.extended.instruments;
  const auto key = [](const char* name) { return std::string("instruments.") + name; };
  in.enabled = parseSectionBool(lookup(key("enabled")), key("enabled"));
  in.backend = parseSectionChoice(lookup(key("backend")), key("backend"), k_instrument_backend_names);
  in.station_layout_path = trimmedValue(lookup(key("station_layout_path")));
  in.station_layout_sha256 = trimmedValue(lookup(key("station_layout_sha256")));
  in.latitude_deg = parseSectionDouble(lookup(key("latitude_deg")), key("latitude_deg"));
  in.longitude_deg = parseSectionDouble(lookup(key("longitude_deg")), key("longitude_deg"));
  in.phase_centre_ra_deg = parseSectionDouble(lookup(key("phase_centre_ra_deg")), key("phase_centre_ra_deg"));
  in.phase_centre_dec_deg = parseSectionDouble(lookup(key("phase_centre_dec_deg")), key("phase_centre_dec_deg"));
  in.hour_angle_start_deg = parseSectionDouble(lookup(key("hour_angle_start_deg")), key("hour_angle_start_deg"));
  in.integration_s = parseSectionDouble(lookup(key("integration_s")), key("integration_s"));
  in.time_steps = parseSectionUint(lookup(key("time_steps")), key("time_steps"));
  in.time_subsamples = parseSectionUint(lookup(key("time_subsamples")), key("time_subsamples"));
  in.channel_subsamples = parseSectionUint(lookup(key("channel_subsamples")), key("channel_subsamples"));
  in.sefd_jy = parseSectionDouble(lookup(key("sefd_jy")), key("sefd_jy"));
  in.noise_seed = parseSectionUint(lookup(key("noise_seed")), key("noise_seed"));
  in.beam = parseSectionChoice(lookup(key("beam")), key("beam"), k_instrument_beam_names);
  in.beam_fwhm_deg = parseSectionDouble(lookup(key("beam_fwhm_deg")), key("beam_fwhm_deg"));
  in.beam_reference_mhz = parseSectionDouble(lookup(key("beam_reference_mhz")), key("beam_reference_mhz"));
  in.dish_diameter_m = parseSectionDouble(lookup(key("dish_diameter_m")), key("dish_diameter_m"));
  in.source_distance_mpc = parseSectionDouble(lookup(key("source_distance_mpc")), key("source_distance_mpc"));
}

void appendExtendedPhysicsSectionsAfterNeutrinos(std::ostream& stream, const ExtendedPhysicsConfig& extended) {
  if (extended.cosmic_rays.model != CosmicRayModelConfig::kNone) {
    stream << "\n[cosmic_rays]\n";
    stream << "model = " << cosmicRayModelConfigToString(extended.cosmic_rays.model) << '\n';
  }
  const ObservablesConfig& o = extended.observables;
  if (o.enabled) {
    stream << "\n[observables]\n";
    stream << "enabled = true\n";
    stream << "interval_steps = " << o.interval_steps << '\n';
    stream << "los_axis = " << losAxisConfigToString(o.los_axis) << '\n';
    stream << "sightlines_per_side = " << o.sightlines_per_side << '\n';
    stream << "output_stem = " << o.output_stem << '\n';
    stream << "lyman_alpha_forest = " << (o.lyman_alpha_forest ? "true" : "false") << '\n';
    stream << "forest_profile = " << forestProfileConfigToString(o.forest_profile) << '\n';
    stream << "forest_pixels = " << o.forest_pixels << '\n';
    stream << "forest_pdf_bins = " << o.forest_pdf_bins << '\n';
    stream << "hi_cube = " << (o.hi_cube ? "true" : "false") << '\n';
    stream << "hi_spin_model = " << hiSpinModelConfigToString(o.hi_spin_model) << '\n';
    stream << "hi_transfer = " << hiTransferConfigToString(o.hi_transfer) << '\n';
    stream << "hi_channels = " << o.hi_channels << '\n';
    stream << "hi_velocity_half_width_km_s = " << o.hi_velocity_half_width_km_s << '\n';
    stream << "hi_collision_table_path = " << o.hi_collision_table_path << '\n';
    stream << "hi_collision_table_sha256 = " << o.hi_collision_table_sha256 << '\n';
    stream << "sz_maps = " << (o.sz_maps ? "true" : "false") << '\n';
    stream << "dispersion_map = " << (o.dispersion_map ? "true" : "false") << '\n';
  }
  const InstrumentsConfig& in = extended.instruments;
  if (in.enabled) {
    stream << "\n[instruments]\nenabled = true\n";
    stream << "backend = " << instrumentBackendConfigToString(in.backend) << '\n';
    stream << "station_layout_path = " << in.station_layout_path << '\n';
    stream << "station_layout_sha256 = " << in.station_layout_sha256 << '\n';
    stream << "latitude_deg = " << in.latitude_deg << '\n';
    stream << "longitude_deg = " << in.longitude_deg << '\n';
    stream << "phase_centre_ra_deg = " << in.phase_centre_ra_deg << '\n';
    stream << "phase_centre_dec_deg = " << in.phase_centre_dec_deg << '\n';
    stream << "hour_angle_start_deg = " << in.hour_angle_start_deg << '\n';
    stream << "integration_s = " << in.integration_s << '\n';
    stream << "time_steps = " << in.time_steps << '\n';
    stream << "time_subsamples = " << in.time_subsamples << '\n';
    stream << "channel_subsamples = " << in.channel_subsamples << '\n';
    stream << "sefd_jy = " << in.sefd_jy << '\n';
    stream << "noise_seed = " << in.noise_seed << '\n';
    stream << "beam = " << instrumentBeamConfigToString(in.beam) << '\n';
    stream << "beam_fwhm_deg = " << in.beam_fwhm_deg << '\n';
    stream << "beam_reference_mhz = " << in.beam_reference_mhz << '\n';
    stream << "dish_diameter_m = " << in.dish_diameter_m << '\n';
    stream << "source_distance_mpc = " << in.source_distance_mpc << '\n';
  }
}

void validateExtendedPhysicsSectionsAfterNeutrinos(const SimulationConfig& config) {
  const CosmicRayConfig& cr = config.extended.cosmic_rays;
  if (cr.model != CosmicRayModelConfig::kNone) {
    requireSection(config.mhd.enabled,
                   "cosmic_rays.model=" + cosmicRayModelConfigToString(cr.model) +
                       " needs the CT magnetic-field authority (mhd.enabled=true) for field-aligned transport");
    requireSection(false,
                   "cosmic_rays.model=" + cosmicRayModelConfigToString(cr.model) +
                       ": the transport, spectral, loss and injection kernels exist (cosmosim/cosmic_rays/*) but the "
                       "MHD level stepper has no CR stage; CR state is not allocated on AMR patches, not carried "
                       "through halo/EMF exchange, regrid, migration or restart, and no SN/AGN feedback budget "
                       "debits CR energy (FAIL_CLOSED; see the extended-physics ledger, CR07)");
  }
  const InstrumentsConfig& in = config.extended.instruments;
  if (!in.enabled) {
    requireSection(in == InstrumentsConfig{},
                   "instruments.enabled=false but other [instruments] keys differ from their defaults; those keys "
                   "would have no effect. Remove them or set instruments.enabled=true");
  } else {
    requireSection(in.backend != InstrumentBackendConfig::kAccelerated,
                   "instruments.backend=accelerated: no GPU direct-sum, NUFFT or antenna-matrix backend exists and no "
                   "accuracy comparison against the direct reference has been made (FAIL_CLOSED; use direct_rime)");
    requireSection(config.extended.observables.enabled && config.extended.observables.hi_cube,
                   "instruments observe the observables.hi_cube product of the same step; enable it");
    requireSection(!in.station_layout_path.empty() && in.station_layout_sha256.size() == 64U,
                   "instruments require station_layout_path and a 64-hex station_layout_sha256 (no array is bundled)");
    requireSection(std::abs(in.latitude_deg) <= 90.0 && std::abs(in.phase_centre_dec_deg) <= 90.0,
                   "instruments latitude/declination must lie in [-90, 90] degrees");
    requireSection(in.integration_s > 0.0 && in.time_steps >= 1U && in.time_subsamples >= 1U &&
                       in.channel_subsamples >= 1U && in.channel_subsamples <= 64U,
                   "instruments need integration_s > 0, time_steps >= 1 and subsamples in [1, 64]");
    requireSection(in.sefd_jy >= 0.0, "instruments.sefd_jy must be >= 0 (0 = noiseless)");
    if (in.beam == InstrumentBeamConfig::kGaussian) {
      requireSection(in.beam_fwhm_deg > 0.0 && in.beam_reference_mhz > 0.0,
                     "instruments.beam=gaussian requires beam_fwhm_deg and beam_reference_mhz > 0 (no default beam)");
    } else if (in.beam == InstrumentBeamConfig::kAiry) {
      requireSection(in.dish_diameter_m > 0.0, "instruments.beam=airy requires dish_diameter_m > 0");
    }
    const bool cosmological_run = buildModePolicy(config.mode).cosmological_comoving_frame;
    requireSection(cosmological_run ? in.source_distance_mpc == 0.0 : in.source_distance_mpc > 0.0,
                   "instruments.source_distance_mpc must be > 0 for isolated runs (angular scale of the cube) and 0 "
                   "for cosmological runs (the S01 transverse comoving distance is used)");
  }
  const ObservablesConfig& o = config.extended.observables;
  if (!o.enabled) {
    requireSection(o == ObservablesConfig{},
                   "observables.enabled=false but other [observables] keys differ from their defaults; those keys "
                   "would have no effect. Remove them or set observables.enabled=true");
    return;
  }
  requireSection(o.interval_steps >= 1U, "observables.interval_steps must be >= 1");
  requireSection(o.sightlines_per_side >= 1U && o.sightlines_per_side <= 4096U,
                 "observables.sightlines_per_side must lie in [1, 4096]");
  requireSection(!o.output_stem.empty() && o.output_stem.find('/') == std::string::npos &&
                     o.output_stem.find("..") == std::string::npos,
                 "observables.output_stem must be a plain file stem (no path separators)");
  requireSection(o.lyman_alpha_forest || o.hi_cube || o.sz_maps || o.dispersion_map,
                 "observables.enabled=true requires at least one product (lyman_alpha_forest, hi_cube, sz_maps, "
                 "dispersion_map)");
  requireSection(config.parallel.mpi_ranks_expected <= 1,
                 "observables: the inline sightline pipeline traverses rank-local leaves only; ordered multi-rank "
                 "segment exchange is not implemented (FAIL_CLOSED for mpi_ranks_expected > 1)");
  requireSection(config.physics.thermochemistry_backend != ThermochemistryBackendMode::kLegacyCooling,
                 "observables need species-resolved composition lanes (n_HI, n_e) from a Part 2 thermochemistry "
                 "backend; the legacy cooling path carries no species and no value is substituted");
  const bool cosmological = buildModePolicy(config.mode).cosmological_comoving_frame;
  if (o.lyman_alpha_forest) {
    requireSection(cosmological, "observables.lyman_alpha_forest needs a cosmological comoving run (Hubble velocity axis)");
    requireSection(o.forest_pixels >= 16U && (o.forest_pixels & (o.forest_pixels - 1U)) == 0U,
                   "observables.forest_pixels must be a power of two >= 16 (1-D flux power spectrum)");
    requireSection(o.forest_pdf_bins >= 2U, "observables.forest_pdf_bins must be >= 2");
  }
  if (o.hi_cube) {
    requireSection(o.hi_channels >= 1U && o.hi_channels <= 65536U, "observables.hi_channels must lie in [1, 65536]");
    requireSection(o.hi_velocity_half_width_km_s > 0.0, "observables.hi_velocity_half_width_km_s must be > 0");
    if (cosmological) {
      requireSection(config.cosmology.radiation_species.enabled,
                     "observables.hi_cube in a cosmological run needs the CMB background T_cmb(a) "
                     "(cosmology.radiation_background=photons_neutrinos); no independent CMB temperature exists");
    }
    if (o.hi_spin_model == HiSpinModelConfig::kCollisionalTable) {
      requireSection(cosmological, "observables.hi_spin_model=collisional_table needs T_gamma(z) (cosmological run)");
      requireSection(!o.hi_collision_table_path.empty() && o.hi_collision_table_sha256.size() == 64U,
                     "observables.hi_spin_model=collisional_table requires hi_collision_table_path and a 64-hex "
                     "hi_collision_table_sha256 (schema chui_21cm_collision_rates_v1; no kappa_10 is bundled)");
    } else {
      requireSection(o.hi_collision_table_path.empty() && o.hi_collision_table_sha256.empty(),
                     "observables.hi_collision_table_* set with hi_spin_model=saturated_kinetic (inert keys)");
    }
  }
}

}  // namespace internal
}  // namespace cosmosim::core
