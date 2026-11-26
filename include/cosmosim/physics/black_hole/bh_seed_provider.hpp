#pragma once

// CHUI-RT-BH-AGN-REF-v1, BH companion section 3.1 (BH01): halo seeding.
//
// The seeder consumes authoritative halo candidate records from a provider;
// it never builds halo truth itself. Reference rule: at most one new seed per
// eligible, UNOCCUPIED halo per seeding event. Occupation is determined from
// membership at the event time (a halo containing any BH member is occupied),
// so halo merger/split and restart never cause reseeding through changed
// indices; the immutable birth key is mixed from the decomposition-independent
// halo key (minimum member particle ID) and the integration step. The gas
// donor is an owned gas cell named by the provider; the seed mass must not
// exceed the donor mass minus the floor (failed candidates are recorded, never
// converted into a synthetic mass source). Gas -> BH conversion transfers
// exactly the seed mass, its momentum, metals/elements, and the donor's
// thermal energy (into the BH engine's captured-energy ledger) through a
// SourceTransactionBatch; IDs come from the existing particle-ID precommit
// registry; new BHs are scheduled for the next legal tick by the caller.
//
// Providers: `halo_fof_serial` (workflows, single rank, FofHaloFinder on the
// local state). Distributed seeding has no bounded authoritative provider and
// remains FAIL_CLOSED at admission (root-global FoF is inadmissible).

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/physics/black_hole/bh_engine_state.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

class ParticleIdPrecommit;

struct BhHaloSeedCandidate {
  std::uint64_t halo_key = 0U;           // minimum member particle ID (decomposition independent)
  double halo_mass_code = 0.0;
  std::array<double, 3> center_stored{};
  bool occupied = false;                 // any BH member at the event time
  std::uint64_t donor_gas_cell_id = 0U;  // densest owned gas member (provider's choice, documented)
  std::uint32_t owner_rank = 0U;
};

struct BhSeedingParameters {
  double halo_mass_threshold_code = 0.0;
  double seed_mass_code = 0.0;
  double mass_floor_code = 1.0e-20;
};

struct BhSeedingReport {
  std::uint64_t candidates = 0U;
  std::uint64_t occupied = 0U;
  std::uint64_t below_threshold = 0U;
  std::uint64_t donor_insufficient = 0U;
  std::uint64_t seeded = 0U;
  double seeded_mass_code = 0.0;
  std::vector<std::uint64_t> seeded_particle_ids;
  std::vector<std::uint32_t> touched_gas_rows;
};

// Applies the reference seeding rule to provider candidates (local rank).
[[nodiscard]] BhSeedingReport applyHaloSeeding(core::SimulationState& state,
                                               std::span<const BhHaloSeedCandidate> candidates,
                                               const BhSeedingParameters& parameters, std::uint64_t step_index,
                                               std::uint64_t interval_index, double event_time_code,
                                               ParticleIdPrecommit& id_precommit);

[[nodiscard]] std::uint64_t bhHaloSeedBirthKey(std::uint64_t halo_key, std::uint64_t step_index) noexcept;

}  // namespace cosmosim::physics
