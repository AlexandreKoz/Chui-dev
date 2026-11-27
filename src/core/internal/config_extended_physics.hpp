#pragma once

// Internal parse/serialize hooks of the CHUI-EXT-PHYS-REF-v1 typed options
// (include/cosmosim/core/config_extended_physics.hpp). config.cpp owns the
// `.param.txt` tokenizer and the consumed-key bookkeeping; this file owns the
// key table, the typed parse, the normalized text and the admission rules of
// the extended sections, so each spelling has exactly one authority.

#include <functional>
#include <ostream>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"

namespace cosmosim::core::internal {

struct ExtendedKeyDefault {
  const char* key;
  const char* default_value;
};

// Every key of [cosmology] radiation background and the extended sections.
[[nodiscard]] const std::vector<ExtendedKeyDefault>& extendedPhysicsKeyDefaults();

// Returns the trimmed raw value of `key` (or its registry default) and marks
// the key consumed when present in the input.
using ExtendedKeyLookup = std::function<std::string(const std::string& key)>;

void parseExtendedPhysicsConfig(const ExtendedKeyLookup& lookup, SimulationConfig& config);

void appendCosmologyRadiationNormalizedText(std::ostream& stream, const CosmologyConfig& cosmology);
void appendExtendedPhysicsNormalizedText(std::ostream& stream, const ExtendedPhysicsConfig& extended);

// [cosmic_rays], [observables], [instruments] (config_extended_physics_sections.cpp).
void parseExtendedPhysicsSectionsAfterNeutrinos(const ExtendedKeyLookup& lookup, SimulationConfig& config);
void appendExtendedPhysicsSectionsAfterNeutrinos(std::ostream& stream, const ExtendedPhysicsConfig& extended);
void validateExtendedPhysicsSectionsAfterNeutrinos(const SimulationConfig& config);
// Keys of those sections, appended to extendedPhysicsKeyDefaults().
[[nodiscard]] const std::vector<ExtendedKeyDefault>& extendedPhysicsSectionKeyDefaults();

}  // namespace cosmosim::core::internal
