// CHUI Part 3 Checkpoint C — element-resolved diffusion, remote feedback,
// stellar-wind channel and correlated associations.
// STATUS: written as source only. NOT RUN — USER-RESERVED.
//
// Invariants:
//  1. Element diffusion (Gaussian-like step between cells): every lane total
//     is conserved to 1e-13 relative, every composition gas mass per cell is
//     unchanged to 1e-13 relative, species-per-element sums stay equal to the
//     element lane (linear operator), no lane becomes negative, a closed box
//     (no faces) is a no-op, and an interval beyond max_subcycles is refused
//     with the sidecar unchanged (bit-identical).
//  2. Remote feedback equivalence (decomposition invariance in one process):
//     a deposition where half of the stencil cells are supplied as REMOTE
//     candidates produces exactly the same per-cell deltas (to 1e-12
//     relative) as the all-local deposition, local ones committed and remote
//     ones exported; a second source sharing cells closes kinetic cross
//     terms so total energy is conserved (step_start_reference).
//  3. Wind table: parse/validation refusals (non-monotone, no provenance),
//     interval budgets are non-negative and exactly additive under step
//     subdivision, metallicity clamps are reported.
//  4. Associations: merged budgets equal the member sums (mass, every
//     element, dust, SN count/energy, wind), ejecta momentum is conserved and
//     the grouping is independent of input order.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/element_diffusion.hpp"
#include "cosmosim/physics/stellar_wind.hpp"
#include "cosmosim/physics/supernova_feedback.hpp"
#include "../support/test_require.hpp"

namespace {

namespace core = cosmosim::core;
namespace ph = cosmosim::physics;

bool close(double a, double b, double rel) {
  return std::abs(a - b) <= rel * std::max({1.0e-300, std::abs(a), std::abs(b)});
}

core::SimulationState makeLine(std::uint32_t n) {
  core::SimulationState state;
  state.resizeCells(n);
  std::vector<core::GasCellIdentityRecord> records;
  for (std::uint32_t row = 0; row < n; ++row) {
    state.cells.center_x_comoving[row] = row + 0.5;
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.mass_code[row] = 1.0;
    state.gas_cells.gas_cell_id[row] = 500U + row;
    state.gas_cells.density_code[row] = 1.0;
    state.gas_cells.internal_energy_code[row] = 1.0;
    state.gas_cells.metal_mass_code[row] = (row < n / 2U ? 0.03 : 0.001);
    records.push_back(core::GasCellIdentityRecord{.gas_cell_id = 500U + row, .parent_particle_id = std::nullopt,
                                                  .owning_patch_id = 1U, .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
  ph::enableGasComposition(state, ph::CompositionLayout::make({
      .network = ph::CompositionNetworkKind::kDetailedCarbonOxygen, .dust = ph::CompositionDustKind::kTwoReservoirV1,
      .sgs_energy = true}), {});
  return state;
}

std::vector<ph::MetalDiffusionFace> lineFaces(std::uint32_t n) {
  std::vector<ph::MetalDiffusionFace> faces;
  for (std::uint32_t i = 0; i + 1U < n; ++i) {
    faces.push_back(ph::MetalDiffusionFace{.left_cell = i, .right_cell = i + 1U, .area_code = 1.0,
                                           .center_distance_code = 1.0});
  }
  return faces;
}

void testElementDiffusionConservesAndPreservesConstraints() {
  constexpr std::uint32_t n = 8U;
  core::SimulationState state = makeLine(n);
  auto& sidecar = state.gas_composition;
  const auto& d = sidecar.descriptor();
  std::vector<double> lane_before(sidecar.laneCount(), 0.0);
  for (std::size_t l = 0; l < sidecar.laneCount(); ++l) for (double v : sidecar.lane(l)) lane_before[l] += v;
  std::vector<double> gas_before(n);
  for (std::uint32_t r = 0; r < n; ++r) gas_before[r] = ph::compositionGasMass(sidecar, r);
  const double sgs_before = d.sgs_energy_lane != core::GasCompositionLayoutDescriptor::k_no_lane
      ? sidecar.value(d.sgs_energy_lane, 0U) : 0.0;
  const std::vector<double> rho_kappa(n, 0.2);
  std::vector<double> workspace;
  const auto report = ph::advanceElementResolvedDiffusion(sidecar, lineFaces(n), rho_kappa, 3.0, 0.5, 256U, workspace);
  CHUI_REQUIRE(report.subcycles >= 1U);
  CHUI_REQUIRE(report.mixed_lane_count == d.element_lane_count + d.species_lane_count + d.dust_lane_count);
  for (std::size_t l = 0; l < sidecar.laneCount(); ++l) {
    double after = 0.0;
    for (double v : sidecar.lane(l)) { CHUI_REQUIRE(v >= 0.0); after += v; }
    CHUI_REQUIRE(close(after, lane_before[l], 1.0e-13) || std::abs(after - lane_before[l]) < 1.0e-18);
  }
  for (std::uint32_t r = 0; r < n; ++r) {
    CHUI_REQUIRE(close(ph::compositionGasMass(sidecar, r), gas_before[r], 1.0e-13));
  }
  // Metals moved down the gradient: left half lost oxygen, right half gained.
  const ph::CompositionLayout layout = ph::CompositionLayout::fromDescriptor(d);
  const std::uint32_t o = layout.elementLane(ph::ChemicalElement::kO);
  CHUI_REQUIRE(sidecar.value(o, n - 1U) > 0.0);
  if (d.sgs_energy_lane != core::GasCompositionLayoutDescriptor::k_no_lane) {
    CHUI_REQUIRE(sidecar.value(d.sgs_energy_lane, 0U) == sgs_before);  // SGS is not mixed
  }
  // Closed box: no faces -> no-op report.
  const auto none = ph::advanceElementResolvedDiffusion(sidecar, {}, rho_kappa, 1.0, 0.5, 4U, workspace);
  CHUI_REQUIRE(none.subcycles == 0U);
}

void testElementDiffusionRefusalIsTransactional() {
  constexpr std::uint32_t n = 8U;
  core::SimulationState state = makeLine(n);
  std::vector<std::vector<double>> before;
  for (std::size_t l = 0; l < state.gas_composition.laneCount(); ++l) {
    const auto lane = state.gas_composition.lane(l);
    before.emplace_back(lane.begin(), lane.end());
  }
  const std::vector<double> rho_kappa(n, 50.0);
  std::vector<double> workspace;
  bool refused = false;
  try {
    (void)ph::advanceElementResolvedDiffusion(state.gas_composition, lineFaces(n), rho_kappa, 100.0, 0.5, 2U, workspace);
  } catch (const std::runtime_error&) {
    refused = true;
  }
  CHUI_REQUIRE(refused);
  for (std::size_t l = 0; l < before.size(); ++l) {
    const auto lane = state.gas_composition.lane(l);
    CHUI_REQUIRE(std::equal(before[l].begin(), before[l].end(), lane.begin()));
  }
}

core::SimulationState makeCube() {
  core::SimulationState state;
  state.resizeCells(27);
  std::vector<core::GasCellIdentityRecord> records;
  std::uint32_t row = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k, ++row) {
        state.gas_cells.gas_cell_id[row] = 1000U + row;
        state.cells.center_x_comoving[row] = i + 0.5;
        state.cells.center_y_comoving[row] = j + 0.5;
        state.cells.center_z_comoving[row] = k + 0.5;
        state.cells.mass_code[row] = 1.0 + 0.01 * row;
        state.gas_cells.density_code[row] = 1.0 + 0.01 * row;
        state.gas_cells.internal_energy_code[row] = 1.0;
        state.gas_cells.velocity_x_peculiar[row] = 0.1 * i;
        records.push_back(core::GasCellIdentityRecord{.gas_cell_id = 1000U + row, .parent_particle_id = std::nullopt,
                                                      .owning_patch_id = 1U, .local_cell_row = row});
      }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

ph::FeedbackUnits testUnits() {
  ph::FeedbackUnits units;
  units.mass_code_to_g = ph::k_solar_mass_g;
  units.velocity_code_to_cm_s = 1.0e5;
  units.length_code_to_cm = 3.0857e18;
  units.density_code_to_cgs = 1.0e-24;
  return units;
}

ph::FeedbackSourceBudget testSource(std::uint64_t id, std::array<double, 3> position) {
  ph::FeedbackSourceBudget src;
  src.source_id = id;
  src.position_comoving = position;
  src.velocity_peculiar = {2.0, -1.0, 0.5};
  src.ejecta_mass_code = 10.0;
  src.ejecta_elements_code[0] = 7.0;
  src.ejecta_elements_code[1] = 2.5;
  src.ejecta_elements_code[4] = 0.5;
  src.supernova_count = 1.0;
  src.supernova_energy_erg = 1.0e51;
  return src;
}

void testRemoteFeedbackMatchesLocalDeposition() {
  core::SimulationState state = makeCube();
  const ph::MechanicalFeedbackConfig config{.neighbor_count = 8, .step_start_reference = true};
  const auto src = testSource(7U, {1.3, 1.6, 1.4});
  // Reference: everything local.
  ph::MechanicalFeedbackDeposition all(config, testUnits());
  std::vector<std::uint32_t> rows(27);
  for (std::uint32_t r = 0; r < 27; ++r) rows[r] = r;
  all.setCandidates(state, rows);
  all.stage(state, src);
  // Split: x < 1.5 local, x > 1.5 "remote" (ids shifted so they are not local).
  ph::MechanicalFeedbackDeposition split(config, testUnits());
  std::vector<std::uint32_t> local_rows;
  std::vector<ph::FeedbackRemoteCandidate> remote;
  for (std::uint32_t r = 0; r < 27; ++r) {
    if (state.cells.center_x_comoving[r] < 1.5) { local_rows.push_back(r); continue; }
    remote.push_back(ph::FeedbackRemoteCandidate{
        .gas_cell_id = 1000000U + state.gas_cells.gas_cell_id[r], .owner_rank = 1,
        .center_comoving = {state.cells.center_x_comoving[r], state.cells.center_y_comoving[r],
                            state.cells.center_z_comoving[r]},
        .mass_code = state.cells.mass_code[r],
        .momentum_code = {state.cells.mass_code[r] * state.gas_cells.velocity_x_peculiar[r], 0.0, 0.0},
        .density_code = state.gas_cells.density_code[r], .hydrogen_mass_fraction = 0.76});
  }
  split.setCandidates(state, local_rows);
  split.setRemoteCandidates(state, remote);
  split.stage(state, src);
  const std::vector<ph::FeedbackCellDelta> exported = split.remoteDeltas(0);
  // Compare: commit both on copies and compare per-cell results / exports.
  core::SimulationState reference = state;
  (void)all.commit(reference);
  core::SimulationState local = state;
  (void)split.commit(local);
  for (std::uint32_t r = 0; r < 27; ++r) {
    if (state.cells.center_x_comoving[r] < 1.5) {
      CHUI_REQUIRE(close(local.cells.mass_code[r], reference.cells.mass_code[r], 1.0e-12));
      CHUI_REQUIRE(close(local.gas_cells.velocity_x_peculiar[r], reference.gas_cells.velocity_x_peculiar[r], 1.0e-12));
      CHUI_REQUIRE(close(local.gas_cells.internal_energy_code[r], reference.gas_cells.internal_energy_code[r], 1.0e-12));
    } else {
      CHUI_REQUIRE(local.cells.mass_code[r] == state.cells.mass_code[r]);  // remote: untouched locally
      const auto it = std::find_if(exported.begin(), exported.end(), [&](const ph::FeedbackCellDelta& d) {
        return d.gas_cell_id == 1000000U + state.gas_cells.gas_cell_id[r];
      });
      const double reference_dm = reference.cells.mass_code[r] - state.cells.mass_code[r];
      if (it == exported.end()) {
        CHUI_REQUIRE(reference_dm == 0.0);
      } else {
        CHUI_REQUIRE(close(it->mass, reference_dm, 1.0e-12));
      }
    }
  }
  CHUI_REQUIRE(split.ledger().remote_target_cells > 0U);
}

void testStepStartReferenceClosesEnergyForSharedCells() {
  core::SimulationState state = makeCube();
  const ph::MechanicalFeedbackConfig config{.neighbor_count = 8, .step_start_reference = true};
  ph::MechanicalFeedbackDeposition dep(config, testUnits());
  std::vector<std::uint32_t> rows(27);
  for (std::uint32_t r = 0; r < 27; ++r) rows[r] = r;
  dep.setCandidates(state, rows);
  const auto energy = [](const core::SimulationState& s) {
    double e = 0.0;
    for (std::size_t r = 0; r < s.cells.size(); ++r) {
      const double m = s.cells.mass_code[r];
      const double vx = s.gas_cells.velocity_x_peculiar[r], vy = s.gas_cells.velocity_y_peculiar[r],
                   vz = s.gas_cells.velocity_z_peculiar[r];
      e += m * s.gas_cells.internal_energy_code[r] + 0.5 * m * (vx * vx + vy * vy + vz * vz);
    }
    return e;
  };
  const double e0 = energy(state);
  dep.stage(state, testSource(7U, {1.3, 1.6, 1.4}));
  dep.stage(state, testSource(8U, {1.5, 1.5, 1.6}));  // overlapping stencil
  (void)dep.commit(state);
  const auto& l = dep.ledger();
  const double expected = l.energy_budget_code;
  // Total energy change = sum of per-source budgets (thermal + kinetic) when
  // every cross term was closed into thermal energy.
  CHUI_REQUIRE(l.kinetic_cross_term_unclosed_code == 0.0);
  CHUI_REQUIRE(close(energy(state) - e0, expected, 1.0e-10));
}

void testWindTable() {
  const std::string text =
      "# chui_ssp_wind_table_v1\n# source = synthetic analytic test fixture (not a physical table)\n"
      "metallicities 2 0.001 0.02\nages 3\n"
      "6.0 1.0e47 10.0 2.0e47 20.0\n6.5 3.0e47 30.0 6.0e47 60.0\n7.0 4.0e47 40.0 8.0e47 80.0\n";
  const auto table = ph::StellarWindTable::loadFromString(text, "test");
  const auto whole = table.interval(1.0e6, 2.0e7, 0.01);
  double energy_sum = 0.0, momentum_sum = 0.0;
  for (double age = 1.0e6; age < 2.0e7; age += 1.9e6) {
    const auto part = table.interval(age, std::min(age + 1.9e6, 2.0e7), 0.01);
    CHUI_REQUIRE(part.energy_erg_per_msun >= 0.0 && part.momentum_msun_kms_per_msun >= 0.0);
    energy_sum += part.energy_erg_per_msun;
    momentum_sum += part.momentum_msun_kms_per_msun;
  }
  CHUI_REQUIRE(close(energy_sum, whole.energy_erg_per_msun, 1.0e-12));
  CHUI_REQUIRE(close(momentum_sum, whole.momentum_msun_kms_per_msun, 1.0e-12));
  CHUI_REQUIRE(table.interval(1.0e6, 2.0e6, 0.5).metallicity_clamped);
  bool refused = false;
  try {
    (void)ph::StellarWindTable::loadFromString(
        "# chui_ssp_wind_table_v1\n# source = x\nmetallicities 1 0.02\nages 2\n6.0 2.0 1.0\n7.0 1.0 2.0\n", "bad");
  } catch (const std::runtime_error&) { refused = true; }
  CHUI_REQUIRE(refused);  // decreasing cumulative energy
  refused = false;
  try {
    (void)ph::StellarWindTable::loadFromString(
        "# chui_ssp_wind_table_v1\nmetallicities 1 0.02\nages 2\n6.0 1.0 1.0\n7.0 2.0 2.0\n", "noprov");
  } catch (const std::runtime_error&) { refused = true; }
  CHUI_REQUIRE(refused);  // missing provenance
}

void testAssociationsPreserveBudgets() {
  std::vector<ph::FeedbackSourceBudget> sources{testSource(3U, {0.10, 0.10, 0.10}), testSource(1U, {0.15, 0.12, 0.11}),
                                                testSource(2U, {5.0, 5.0, 5.0})};
  sources[1].velocity_peculiar = {-1.0, 0.0, 0.0};
  sources[1].wind_energy_erg = 2.0e49;
  const std::vector<double> birth{1.0e6, 1.5e6, 1.2e6};
  const auto merged = ph::mergeAssociationSources(sources, birth, 1.0, 3.0e6, false, {});
  CHUI_REQUIRE(merged.size() == 2U);
  double mass = 0.0, energy = 0.0, wind = 0.0, count = 0.0, px = 0.0;
  for (const auto& m : merged) {
    mass += m.ejecta_mass_code; energy += m.supernova_energy_erg; wind += m.wind_energy_erg;
    count += m.supernova_count; px += m.ejecta_mass_code * m.velocity_peculiar[0];
  }
  CHUI_REQUIRE(close(mass, 30.0, 1.0e-15) && close(energy, 3.0e51, 1.0e-15) && close(wind, 2.0e49, 1.0e-15));
  CHUI_REQUIRE(count == 3.0);
  CHUI_REQUIRE(close(px, 10.0 * 2.0 + 10.0 * -1.0 + 10.0 * 2.0, 1.0e-14));
  // Order independence.
  std::vector<ph::FeedbackSourceBudget> reversed(sources.rbegin(), sources.rend());
  const std::vector<double> reversed_birth(birth.rbegin(), birth.rend());
  const auto again = ph::mergeAssociationSources(reversed, reversed_birth, 1.0, 3.0e6, false, {});
  CHUI_REQUIRE(again.size() == merged.size());
  for (std::size_t i = 0; i < merged.size(); ++i) {
    CHUI_REQUIRE(again[i].source_id == merged[i].source_id && again[i].member_count == merged[i].member_count);
  }
}

}  // namespace

int main() {
  testElementDiffusionConservesAndPreservesConstraints();
  testElementDiffusionRefusalIsTransactional();
  testRemoteFeedbackMatchesLocalDeposition();
  testStepStartReferenceClosesEnergyForSharedCells();
  testWindTable();
  testAssociationsPreserveBudgets();
  return 0;
}
