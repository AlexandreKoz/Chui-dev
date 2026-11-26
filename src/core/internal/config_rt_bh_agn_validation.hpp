#pragma once

#include "cosmosim/core/config.hpp"

namespace cosmosim::core::internal {

// Typed admission of the RT/BH/AGN options (config_rt_bh_agn.hpp). Throws
// ConfigError with a precise reason; never degrades a model silently.
void validateRtBhAgnConfig(const SimulationConfig& config);

}  // namespace cosmosim::core::internal
