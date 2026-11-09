#pragma once

#include <array>
#include <type_traits>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
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
  // CHUI Part 3 wind/pre-SN mechanical channel (no mass: wind mass is in the
  // SSP CC continuous return). Energy enters the energy closure; momentum is
  // deposited radially like presn_momentum_code.
  double wind_energy_erg = 0.0;
  double wind_momentum_code = 0.0;
  // Number of member stars merged into this source (associations); 1 = single.
  std::uint32_t member_count = 1U;
  std::uint32_t reserved_padding = 0U;
  // G0.3: kinetic energy of the member ejecta in the association centre-of-
  // momentum frame, 1/2 sum_i m_i |v_i - V|^2 (code mass * code velocity^2).
  // Depositing the merged ejecta at the single velocity V would otherwise
  // discard it; the deposition adds it to the available energy so the merged
  // and unmerged models inject the same total energy. Zero for single stars.
  double association_dispersion_energy_code = 0.0;
};
static_assert(std::is_trivially_copyable_v<FeedbackSourceBudget>,
              "FeedbackSourceBudget is exchanged as a byte record across ranks");

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
  // CHUI Part 3: evaluate every source against the STEP-START target state
  // (no dependence on the order of earlier staged sources). Required for
  // results independent of the MPI decomposition; the kinetic-energy cross
  // terms of several sources sharing a cell are closed per cell at commit so
  // the total injected energy stays exact (see FeedbackLedger).
  bool step_start_reference = false;
};

// CHUI Part 3: remote (other-rank) leaf gas cell imported as a deposition
// candidate (step-start properties; owner applies the resulting delta).
struct FeedbackRemoteCandidate {
  std::uint64_t gas_cell_id = 0U;
  int owner_rank = 0;
  std::array<double, 3> center_comoving{};
  double mass_code = 0.0;
  std::array<double, 3> momentum_code{};
  double density_code = 0.0;
  double hydrogen_mass_fraction = 0.76;
};

// Owner-targeted wire record of one cell's accumulated feedback delta.
struct FeedbackCellDelta {
  std::uint64_t gas_cell_id = 0U;
  int source_rank = 0;
  int owner_rank = 0;
  double mass = 0.0;
  std::array<double, 3> momentum{};
  double thermal_energy = 0.0;
  double sgs_energy = 0.0;
  double intended_kinetic = 0.0;
  std::array<double, k_chemical_element_count> elements{};
  std::array<double, k_dust_reservoir_count> dust{};
  double swept_mass_for_dust = 0.0;
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
  // Per-cell kinetic closure (step_start_reference): sum over cells of
  // (intended kinetic change - actual kinetic change) moved into thermal
  // energy, and the part that could not be (thermal floor) — ledgered.
  double kinetic_cross_term_to_thermal_code = 0.0;
  double kinetic_cross_term_unclosed_code = 0.0;
  std::uint64_t remote_target_cells = 0;
  std::uint64_t resolved_cells = 0;
  std::uint64_t unresolved_cells = 0;
  std::uint64_t sources = 0;
};

// CHUI Part 3 RESEARCH (EXPERIMENTAL/UNQUALIFIED): correlated stellar
// associations / superbubbles (G0.3 decomposition-independent contract).
//
// Discovery: friends-of-friends over the GLOBAL set of active feedback
// sources. Sources i and j are linked when their minimum-image comoving
// separation is <= radius_code AND |t_birth,i - t_birth,j| <= window_years.
// An association is a connected component; its stable identity is the
// minimum member source_id (the "leader"). The component set depends only on
// the global source set, never on its distribution over ranks or input order.
//
// Merge: members are combined in ascending source_id order (deterministic
// floating-point summation). Summed ejecta, elements, dust, SN count/energy,
// pre-SN and wind budgets; position = ejecta-mass-weighted minimum-image
// centroid about the leader (arithmetic mean if the ejecta mass is zero),
// wrapped into the periodic box; velocity V = total ejecta momentum / total
// ejecta mass (conserves vector momentum); the internal velocity-dispersion
// kinetic energy 1/2 sum m_i |v_i - V|^2 (plus any member's own) is carried in
// association_dispersion_energy_code. Mass, elements, dust, energy and vector
// momentum are therefore those of the unmerged members.
struct AssociationSourceDescriptor {
  std::uint64_t source_id = 0U;
  std::array<double, 3> position_comoving{};
  double birth_years = 0.0;
  std::int32_t owner_rank = 0;
  std::int32_t reserved = 0;
};
static_assert(std::is_trivially_copyable_v<AssociationSourceDescriptor>);

// Returns, for each descriptor of `sources` (any order, unique source ids),
// the leader source id of its association. Throws on duplicate ids or
// non-positive radius/window.
[[nodiscard]] std::vector<std::uint64_t> discoverAssociationLeaders(
    std::span<const AssociationSourceDescriptor> sources, double radius_code, double window_years,
    bool periodic, const std::array<double, 3>& box_extent_comoving);

// Merges the members of ONE association (any order; sorted internally by id).
[[nodiscard]] FeedbackSourceBudget mergeAssociationMembers(
    std::span<const FeedbackSourceBudget> members, bool periodic, const std::array<double, 3>& box_extent_comoving);

// Single-rank convenience: discovery + merge over `sources` (output ordered by
// leader id). With all sources on one rank this equals the distributed result.
[[nodiscard]] std::vector<FeedbackSourceBudget> mergeAssociationSources(
    std::span<const FeedbackSourceBudget> sources, std::span<const double> birth_years,
    double radius_code, double window_years, bool periodic, const std::array<double, 3>& box_extent_comoving);

// Exposed for tests: KO15 terminal momentum per SN (Msun km/s).
[[nodiscard]] double kimOstriker2015TerminalMomentumMsunKms(double energy_per_sn_erg, double n_h_cm3);

class MechanicalFeedbackDeposition {
 public:
  MechanicalFeedbackDeposition(MechanicalFeedbackConfig config, FeedbackUnits units);

  // Patch 53 (P1-52-02): optional MemoryGovernor for the candidate indices
  // (kPhaseResident). Set before setCandidates; null = accounting disabled
  // (no-governor policy of the memory specification).
  void setMemoryGovernor(core::MemoryGovernor* governor) noexcept { m_governor = governor; }
  // Candidate cells: owned leaf rows. Builds an x-sorted index. The index
  // (sorted rows, sorted x, row-slot map) is admitted before it is built.
  void setCandidates(const core::SimulationState& state, std::vector<std::uint32_t> owned_leaf_rows);
  // Part 3: append remote candidates (after setCandidates). Deterministic:
  // candidates are ordered by (x, stable id) and k-NN ties by stable id.
  // Patch 53: strongly exception-safe (the new index is built beside the old
  // one and swapped in only after validation, so a failure leaves the local
  // index intact); the new index is admitted with old/new coexistence before
  // it is allocated; `remote_storage_lease` (the governed exchange lease of
  // `remote`) is adopted together with the storage it covers.
  void setRemoteCandidates(const core::SimulationState& state, std::vector<FeedbackRemoteCandidate> remote,
                           core::MemoryReservation remote_storage_lease = {});
  // Squared distance of the k-th nearest LOCAL candidate (infinity if fewer
  // than k local candidates). Any remote cell closer than this may enter the
  // global stencil, so it bounds the remote halo query of a source.
  [[nodiscard]] double localStencilRadius2(const core::SimulationState& state,
                                           const std::array<double, 3>& position) const;
  // Local candidates (dense rows) among the k nearest within r2_max.
  [[nodiscard]] std::vector<std::uint32_t> localCandidatesWithin(const core::SimulationState& state,
                                                                 const std::array<double, 3>& position,
                                                                 double r2_max) const;
  // Deltas staged on remote candidates, for delivery to their owners.
  [[nodiscard]] std::vector<FeedbackCellDelta> remoteDeltas(int source_rank) const;
  // Merge a delta delivered by another rank for a local cell (before commit).
  void addImportedDelta(const core::SimulationState& state, const FeedbackCellDelta& delta);
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
    double intended_kinetic = 0.0;
  };
  struct CandidateProps {
    std::array<double, 3> center{};
    double mass = 0.0;
    std::array<double, 3> momentum{};
    double density = 0.0;
    double x_h = 0.76;
    std::uint64_t stable_id = 0U;
  };
  [[nodiscard]] CandidateProps props(const core::SimulationState& state, std::uint32_t candidate) const;
  [[nodiscard]] std::vector<std::pair<double, std::uint32_t>> nearest(const core::SimulationState& state,
                                                                      const std::array<double, 3>& p) const;
  RowDelta& delta(std::uint32_t row);

  MechanicalFeedbackConfig m_config;
  FeedbackUnits m_units;
  core::MemoryGovernor* m_governor = nullptr;
  // Leases are declared before the storage they cover (implicit destruction
  // frees the storage first).
  core::MemoryReservation m_index_lease;   // m_sorted_rows + m_sorted_x + m_row_slot
  core::MemoryReservation m_remote_lease;  // m_remote (adopted exchange lease)
  std::vector<std::uint32_t> m_sorted_rows;  // by x
  std::vector<double> m_sorted_x;
  std::vector<RowDelta> m_rows;
  std::vector<std::uint32_t> m_row_slot;  // candidate -> slot+1 (0 = none)
  std::uint32_t m_local_count = 0U;       // candidates >= this are remote
  std::vector<FeedbackRemoteCandidate> m_remote;
  FeedbackLedger m_ledger;
};

}  // namespace cosmosim::physics
