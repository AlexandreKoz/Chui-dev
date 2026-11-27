// CHUI-EXT-PHYS-REF-v1 I02: [cosmic_rays], [observables] and [instruments]
// typed options. Source-only pass 01: NOT compiled or executed.

#include "core/internal/config_extended_physics.hpp"

namespace cosmosim::core::internal {

const std::vector<ExtendedKeyDefault>& extendedPhysicsSectionKeyDefaults() {
  static const std::vector<ExtendedKeyDefault> keys = {};
  return keys;
}

void parseExtendedPhysicsSectionsAfterNeutrinos(const ExtendedKeyLookup& /*lookup*/, SimulationConfig& /*config*/) {}

void appendExtendedPhysicsSectionsAfterNeutrinos(std::ostream& /*stream*/, const ExtendedPhysicsConfig& /*extended*/) {}

void validateExtendedPhysicsSectionsAfterNeutrinos(const SimulationConfig& /*config*/) {}

}  // namespace cosmosim::core::internal
