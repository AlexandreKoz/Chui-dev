#include "cosmosim/physics/black_hole/bh_seed_provider.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/source_transaction.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace cosmosim::physics {

std::uint64_t bhHaloSeedBirthKey(std::uint64_t halo_key, std::uint64_t step_index) noexcept {
  // SplitMix64 mixing with a domain tag distinct from the legacy cell-based key.
  std::uint64_t value = halo_key ^ (step_index * 0x9e3779b97f4a7c15ULL) ^ 0x42484845414c4f31ULL;  // "BHHEALO1"
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  value ^= value >> 31U;
  return value == 0U ? 0x510e527fade682d1ULL : value;
}

BhSeedingReport applyHaloSeeding(core::SimulationState& state, std::span<const BhHaloSeedCandidate> candidates,
                                 const BhSeedingParameters& par, std::uint64_t step_index,
                                 std::uint64_t interval_index, double event_time_code,
                                 ParticleIdPrecommit& registry) {
  if (!(par.seed_mass_code > 0.0) || !(par.halo_mass_threshold_code > 0.0)) {
    throw std::invalid_argument("applyHaloSeeding: seed mass and halo threshold must be > 0");
  }
  BhSeedingReport report;
  struct Accepted {
    const BhHaloSeedCandidate* candidate;
    std::uint32_t donor_row;
  };
  std::vector<Accepted> accepted;
  std::vector<std::uint64_t> seen;
  for (const BhHaloSeedCandidate& c : candidates) {
    ++report.candidates;
    if (std::find(seen.begin(), seen.end(), c.halo_key) != seen.end()) {
      throw std::runtime_error("applyHaloSeeding: duplicated halo key in one seeding event");
    }
    seen.push_back(c.halo_key);
    if (c.occupied) {
      ++report.occupied;
      continue;
    }
    if (c.halo_mass_code < par.halo_mass_threshold_code) {
      ++report.below_threshold;
      continue;
    }
    const std::optional<std::uint32_t> row = state.rowForGasCellId(c.donor_gas_cell_id);
    if (!row.has_value() || state.cells.mass_code[*row] <= par.seed_mass_code + par.mass_floor_code) {
      ++report.donor_insufficient;  // recorded failure; no synthetic mass
      continue;
    }
    accepted.push_back(Accepted{&c, *row});
  }
  // Deterministic order by halo key.
  std::sort(accepted.begin(), accepted.end(),
            [](const Accepted& a, const Accepted& b) { return a.candidate->halo_key < b.candidate->halo_key; });
  // Two halos may name the same donor row: only the first (lowest key) seeds.
  std::vector<Accepted> unique;
  for (const Accepted& a : accepted) {
    const bool clash = std::any_of(unique.begin(), unique.end(),
                                   [&](const Accepted& u) { return u.donor_row == a.donor_row; });
    if (clash) {
      ++report.donor_insufficient;
    } else {
      unique.push_back(a);
    }
  }
  std::vector<std::uint64_t> keys;
  keys.reserve(unique.size());
  for (const Accepted& a : unique) keys.push_back(bhHaloSeedBirthKey(a.candidate->halo_key, step_index));
  const std::vector<std::uint64_t> ids = registry.precommit(state, keys);
  struct Finish {
    ParticleIdPrecommit& r;
    ~Finish() { r.finishPrecommit(); }
  } finish{registry};
  if (ids.size() != unique.size()) throw std::runtime_error("applyHaloSeeding: ID precommit returned a wrong count");
  if (unique.empty()) return report;

  // Gas side: one conservative transaction per seed (no mutation before the
  // whole batch prepared and validated).
  SourceTransactionBatch batch(TransactionFloors{par.mass_floor_code, 0.0});
  struct Transfer {
    double mass;
    std::array<double, 3> velocity;
    double thermal;
    std::array<double, 3> position;
  };
  std::vector<Transfer> transfers;
  for (std::size_t i = 0; i < unique.size(); ++i) {
    const std::uint32_t row = unique[i].donor_row;
    const double m = state.cells.mass_code[row];
    const double f = par.seed_mass_code / m;
    const std::array<double, 3> v{state.gas_cells.velocity_x_peculiar[row], state.gas_cells.velocity_y_peculiar[row],
                                  state.gas_cells.velocity_z_peculiar[row]};
    SourceTransactionPlan plan;
    plan.id = SourceTransactionId{unique[i].candidate->halo_key, interval_index, SourceChannel::kBhSeedConversion, 0U};
    plan.topology_epoch = state.cellIndexGeneration();
    plan.intended_time_code = event_time_code;
    GasCellIncrement inc;
    inc.gas_cell_id = state.gas_cells.gas_cell_id[row];
    inc.local_row = row;
    inc.withdrawal_fraction = f;
    inc.delta.mass_code = -par.seed_mass_code;
    const double thermal = par.seed_mass_code * state.gas_cells.internal_energy_code[row];
    double v2 = 0.0;
    for (std::size_t a = 0; a < 3; ++a) {
      inc.delta.momentum_code[a] = -par.seed_mass_code * v[a];
      v2 += v[a] * v[a];
    }
    inc.delta.total_energy_code = -(thermal + 0.5 * par.seed_mass_code * v2);
    plan.cells.push_back(inc);
    plan.source_increment = inc.delta.negated();
    batch.add(std::move(plan));
    transfers.push_back(Transfer{par.seed_mass_code, v, thermal,
                                 {state.cells.center_x_comoving[row], state.cells.center_y_comoving[row],
                                  state.cells.center_z_comoving[row]}});
  }
  if (!batch.prepare(state) || !batch.validate() || !batch.authorize(true)) {
    throw std::runtime_error("applyHaloSeeding: seed conversion rejected: " + batch.report().failure);
  }
  registry.preparePopulationGrowth(state, unique.size(), core::ParticleSpecies::kBlackHole);
  report.touched_gas_rows = batch.commit(state);
  BlackHoleEngineStore store(state);
  for (std::size_t i = 0; i < unique.size(); ++i) {
    const std::size_t p = state.particles.size();
    if (p >= static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
      throw std::overflow_error("applyHaloSeeding: local particle count exceeds the uint32 sidecar index");
    }
    state.resizeParticles(p + 1U);
    state.particles.position_x_comoving[p] = transfers[i].position[0];
    state.particles.position_y_comoving[p] = transfers[i].position[1];
    state.particles.position_z_comoving[p] = transfers[i].position[2];
    state.particles.velocity_x_peculiar[p] = transfers[i].velocity[0];
    state.particles.velocity_y_peculiar[p] = transfers[i].velocity[1];
    state.particles.velocity_z_peculiar[p] = transfers[i].velocity[2];
    state.particles.mass_code[p] = transfers[i].mass;
    state.particle_sidecar.particle_id[p] = ids[i];
    state.particle_sidecar.sfc_key[p] = ids[i];
    state.particle_sidecar.species_tag[p] = static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);
    state.particle_sidecar.particle_flags[p] = 0U;
    state.particle_sidecar.owning_rank[p] = unique[i].candidate->owner_rank;
    const std::size_t b = state.black_holes.size();
    state.black_holes.resize(b + 1U);
    state.black_holes.particle_index[b] = static_cast<std::uint32_t>(p);
    state.black_holes.host_cell_index[b] = unique[i].donor_row;
    state.black_holes.subgrid_mass_code[b] = transfers[i].mass;
    state.black_holes.accretion_rate_code[b] = 0.0;
    state.black_holes.feedback_energy_code[b] = 0.0;
    state.black_holes.eddington_ratio[b] = 0.0;
    state.black_holes.cumulative_accreted_mass_code[b] = 0.0;
    state.black_holes.cumulative_feedback_energy_code[b] = 0.0;
    state.black_holes.duty_cycle_active_time_code[b] = 0.0;
    state.black_holes.duty_cycle_total_time_code[b] = 0.0;
    BlackHoleEngineRecord rec{};
    rec.horizon_mass_code = transfers[i].mass;
    // The seed converts gas rest mass into the BH; the donor's thermal and
    // kinetic energy is not carried by a collisionless particle beyond its
    // momentum: thermal energy goes to the named captured-energy ledger.
    rec.cumulative_captured_gas_energy_code = transfers[i].thermal;
    rec.flags = kBhFlagInitialized | kBhFlagSeededByHalo;
    rec.seed_halo_key = unique[i].candidate->halo_key;
    rec.source_time_code = event_time_code;
    rec.interval_index = interval_index;
    store.put(ids[i], rec);
    ++report.seeded;
    report.seeded_mass_code += transfers[i].mass;
    report.seeded_particle_ids.push_back(ids[i]);
  }
  state.notifyParticleMassMutation();
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)] += report.seeded;
  state.rebuildSpeciesIndex();
  store.store(state);
  // The gravity-source generation is advanced once, collectively, by the
  // workflow at the end of the discrete source epoch (communicator-consistent).
  return report;
}

}  // namespace cosmosim::physics
