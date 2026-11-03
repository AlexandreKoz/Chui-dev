#pragma once

#include <memory>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/chemical_network.hpp"
#include "cosmosim/physics/metal_cooling_table.hpp"
#include "cosmosim/physics/radiation_background.hpp"
#include "cosmosim/physics/thermochemistry_operator.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

// ---------------------------------------------------------------------------
// CHUI Part 2 translation of validated core::PhysicsConfig into physics
// objects. Resolved ONCE at runtime construction (never inside cell loops).
// ---------------------------------------------------------------------------
[[nodiscard]] ThermochemistryBackend thermochemistryBackendFromConfig(const core::PhysicsConfig& physics);
[[nodiscard]] bool part2ThermochemistrySelected(const core::PhysicsConfig& physics);
[[nodiscard]] CompositionLayoutRequest compositionLayoutRequestFromConfig(const core::PhysicsConfig& physics);
[[nodiscard]] PrimordialCompositionPolicy primordialPolicyFromConfig(const core::PhysicsConfig& physics);

// Idempotent: enables the composition sidecar with the configured layout on
// a state whose gas carries none yet (fresh initial conditions), using the
// explicit primordial + solar-pattern initialization policy; verifies the
// layout when it is already present (e.g. restored from a restart). No-op
// when Part 2 composition is not selected. Returns true if it initialized.
bool ensureGasCompositionInitialized(core::SimulationState& state, const core::PhysicsConfig& physics);

// Data products owned by the runtime (immutable, shared across cells).
struct BaryonicDataProducts {
  std::shared_ptr<const UvBackgroundTable> uv_background;
  std::shared_ptr<const MetalCoolingTable> metal_cooling_table;
};
// Loads every external data file required by the configuration. Missing or
// corrupt data throw SourceModelUnavailableError with acquisition guidance;
// nothing is synthesized.
[[nodiscard]] BaryonicDataProducts loadBaryonicDataProducts(const core::PhysicsConfig& physics);

[[nodiscard]] ThermochemistryOperatorConfig makeThermochemistryOperatorConfig(
    const core::SimulationConfig& config, const core::UnitSystem& units, const BaryonicDataProducts& data);

}  // namespace cosmosim::physics
