#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

// ===========================================================================
// CHUI Part 2 conservative stellar feedback deposition
// (fb_supernova_model = mechanical_resolved_unresolved).
//
// Per source s with ejecta mass m_ej, element/dust masses, N_SN supernovae of
// total energy E (erg) and an optional pre-SN momentum budget:
//  1. Stencil: the k nearest OWNED LEAF gas cells (periodic minimum image in
//     periodic modes). Weights w_i proportional to the solid angle proxy
//     V_i^(2/3) / r_i^2, normalized to sum 1 over the actually owned cells.
//  2. Mass, elements, dust: dm_i = w_i m_ej (same weights for every lane).
//  3. Momentum magnitude per cell (Hopkins et al. 2018, MNRAS 477, 1578,
//     eq. 18-19 form):  p_i = w_i * min( p_0 sqrt(1 + m_i/(w_i m_ej)),
//                                         p_term(n_i, Z_i) N_SN )
//     p_0 = sqrt(2 m_ej E); terminal momentum per SN from Kim & Ostriker
//     (2015, ApJ 802, 99) eq. 17:  p_term = 2.17e5 Msun km/s E51^0.93 n^-0.13
//     (n = n_H of the cell, proper cm^-3). A cell where the energy-conserving
//     branch is smaller than the terminal momentum is RESOLVED (Sedov phase
//     captured), else UNRESOLVED (terminal momentum). Both are counted.
//  4. Vector priority (mathematically consistent, documented): directions
//     n_i point from the source to the cell centre. The net vector sum
//     R = sum p_i n_i of an asymmetric stencil is removed by subtracting
//     R * dm'_i / sum dm' from each cell (mass-weighted), so the injected
//     NET momentum in the star frame is exactly zero; the radial scalar sum
//     before/after this correction is ledgered (radial deficit). In the lab
//     frame each cell also receives dm_i v_star (ejecta carry the source
//     velocity), conserving total momentum exactly.
//  5. Energy closure: total injected energy E + 0.5 m_ej v_star^2 is split
//     into the exact kinetic-energy change of every cell (computed from the
//     final momenta) and thermal energy (the remainder, distributed with
//     w_i). If the kinetic change would exceed the budget, the momenta are
//     scaled down (bisection) until it fits; nothing is created. A configured
//     fraction of the thermal part goes to the SGS reservoir instead.
//  6. Dust: SN shocks destroy solids in the swept mass (McKee 1989 form,
//     destroyDustBySupernova), applied once per committed source.
//  7. Staging: all contributions of all sources accumulate into per-row
//     deltas (deterministic order: sources sorted by star particle id); the
//     state is written only by commit(), after budgets validated. An event
//     is consumed only by a successful commit (no loss/duplication on retry).
// Topology: single-rank (fixed grid or local AMR, periodic or physical
// boundaries). Distributed deposition across MPI ranks requires Part 3
// remote ghost deposition; the runtime REFUSES it (fail closed) rather than
// renormalizing onto local cells.
// ===========================================================================

struct FeedbackSourceBudget {
  std::uint64_t source_id = 0;  // star particle id (deterministic order)
  std::array<double, 3> position_comoving{};
  std::array<double, 3> velocity_peculiar{};
  double ejecta_mass_code = 0.0;
  std::array<double, k_chemical_element_count> ejecta_elements_code{};
  std::array<double, k_dust_reservoir_count> ejecta_dust_code{};
  double supernova_count = 0.0;
  double supernova_energy_erg = 0.0;
  double presn_momentum_code = 0.0;  // radiation/wind momentum budget
};

struct FeedbackUnits {
  double mass_code_to_g = 1.0;
  double velocity_code_to_cm_s = 1.0;
  double length_code_to_cm = 1.0;
  double density_code_to_cgs = 1.0;
  double scale_factor = 1.0;
  bool comoving = false;
  bool periodic = false;
  std::array<double, 3> box_extent_comoving{};
};

struct MechanicalFeedbackConfig {
  std::uint32_t neighbor_count = 32U;
  double sgs_energy_fraction = 0.0;
  double dust_sn_swept_mass_per_1e51_erg_g = 6800.0 * k_solar_mass_g;
  double dust_sn_destruction_efficiency = 0.3;
};

struct FeedbackLedger {
  double mass_injected_code = 0.0;
  std::array<double, k_chemical_element_count> elements_injected_code{};
  double dust_injected_code = 0.0;
  double dust_destroyed_code = 0.0;
  double energy_budget_code = 0.0;
  double kinetic_energy_injected_code = 0.0;
  double thermal_energy_injected_code = 0.0;
  double sgs_energy_injected_code = 0.0;
  double radial_momentum_requested_code = 0.0;
  double radial_momentum_delivered_code = 0.0;
  double net_vector_momentum_residual_code = 0.0;  // |sum dp - m_ej v_star| (round-off)
  double momentum_scaled_down_sources = 0.0;
  std::uint64_t resolved_cells = 0;
  std::uint64_t unresolved_cells = 0;
  std::uint64_t sources = 0;
};

// Exposed for tests: KO15 terminal momentum per SN (Msun km/s).
[[nodiscard]] double kimOstriker2015TerminalMomentumMsunKms(double energy_per_sn_erg, double n_h_cm3);

class MechanicalFeedbackDeposition {
 public:
  MechanicalFeedbackDeposition(MechanicalFeedbackConfig config, FeedbackUnits units);

  // Candidate cells: owned leaf rows. Builds an x-sorted index.
  void setCandidates(const core::SimulationState& state, std::vector<std::uint32_t> owned_leaf_rows);
  // Stage one source (no state mutation). Throws on invalid input.
  void stage(const core::SimulationState& state, const FeedbackSourceBudget& source);
  // Apply every staged delta; returns touched rows (sorted, unique).
  std::vector<std::uint32_t> commit(core::SimulationState& state);
  [[nodiscard]] const FeedbackLedger& ledger() const noexcept { return m_ledger; }
  [[nodiscard]] std::size_t stagedRows() const noexcept { return m_rows.size(); }

 private:
  struct RowDelta {
    std::uint32_t row = 0;
    double mass = 0.0;
    std::array<double, 3> momentum{};
    double thermal_energy = 0.0;
    double sgs_energy = 0.0;
    std::array<double, k_chemical_element_count> elements{};
    std::array<double, k_dust_reservoir_count> dust{};
    double swept_mass_for_dust = 0.0;
  };
  [[nodiscard]] std::vector<std::pair<double, std::uint32_t>> nearest(const core::SimulationState& state,
                                                                      const std::array<double, 3>& p) const;
  RowDelta& delta(std::uint32_t row);

  MechanicalFeedbackConfig m_config;
  FeedbackUnits m_units;
  std::vector<std::uint32_t> m_sorted_rows;  // by x
  std::vector<double> m_sorted_x;
  std::vector<RowDelta> m_rows;
  std::vector<std::uint32_t> m_row_slot;  // dense row -> slot+1 (0 = none)
  FeedbackLedger m_ledger;
};

}  // namespace cosmosim::physics
