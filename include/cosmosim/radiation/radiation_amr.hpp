#pragma once

// CHUI-RT-BH-AGN-REF-v1, RT companion section 7 (RT06): conservative,
// realizability-preserving ratio-2 transfer operators for radiation moments.
//
// Restriction: volume-weighted average of the 8 children (admissible by
// convexity of the cone E >= 0, |F| <= c E). Prolongation: parent value plus a
// jointly limited linear correction (MC slopes from the parent's axis
// neighbours, one theta in [0,1] for all eight children so every child is
// admissible); children average exactly to the parent (conservation).
// STATUS: library operators used by the level-synchronous ghost logic
// semantics and available for the regrid transaction; dynamic regrid with RT
// is FAIL_CLOSED at admission until moments migrate through
// amr_production_regrid (recorded gap).

#include <array>

#include "cosmosim/radiation/radiation_closure.hpp"

namespace cosmosim::radiation {

[[nodiscard]] M1State restrictRadiationChildren(const std::array<M1State, 8>& children) noexcept;

// neighbours[axis][0] = low neighbour, [1] = high neighbour of the parent.
// Child ordering: index = cx + 2 (cy + 2 cz), c* in {0,1} (low/high half).
[[nodiscard]] std::array<M1State, 8> prolongRadiationParent(const M1State& parent,
                                                            const std::array<std::array<M1State, 2>, 3>& neighbours,
                                                            double c) noexcept;

}  // namespace cosmosim::radiation
