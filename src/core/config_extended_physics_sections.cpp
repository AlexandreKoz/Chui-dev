// CHUI-EXT-PHYS-REF-v1 I02: [cosmic_rays], [observables] and [instruments]
// typed options. Source-only pass 01: NOT compiled or executed.

#include <algorithm>
#include <array>
#include <cctype>
#include <string_view>

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

}  // namespace

std::string cosmicRayModelConfigToString(CosmicRayModelConfig value) {
  return sectionChoiceToString(value, k_cosmic_ray_model_names, "CosmicRayModelConfig");
}

namespace internal {

const std::vector<ExtendedKeyDefault>& extendedPhysicsSectionKeyDefaults() {
  static const std::vector<ExtendedKeyDefault> keys = {
      {"cosmic_rays.model", "none"},
  };
  return keys;
}

void parseExtendedPhysicsSectionsAfterNeutrinos(const ExtendedKeyLookup& lookup, SimulationConfig& config) {
  config.extended.cosmic_rays.model =
      parseSectionChoice(lookup("cosmic_rays.model"), "cosmic_rays.model", k_cosmic_ray_model_names);
}

void appendExtendedPhysicsSectionsAfterNeutrinos(std::ostream& stream, const ExtendedPhysicsConfig& extended) {
  if (extended.cosmic_rays.model != CosmicRayModelConfig::kNone) {
    stream << "\n[cosmic_rays]\n";
    stream << "model = " << cosmicRayModelConfigToString(extended.cosmic_rays.model) << '\n';
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
}

}  // namespace internal
}  // namespace cosmosim::core
