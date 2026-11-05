#pragma once

// CHUI Part 3 Checkpoint B — automatic production refine/derefine (audit F6).
//
// Ownership: SimulationState is the only hierarchy authority (it stores leaf
// patches only; a refined patch is replaced by its eight children and a
// derefined octet by its parent). This module adds a PURE, deterministic
// planner and a transactional applier on top of the existing conservative
// single-patch operations (refineProductionPatchInSimulationState /
// derefineProductionPatchInSimulationState), which already admit old+new
// candidate storage through the MemoryGovernor before any mutation, conserve
// mass/momentum/total energy exactly (piecewise-constant prolongation, volume
// sum restriction) and, since Part 3, regrid the composition sidecar
// (elements, species, dust, SGS) with identical weights.
//
// Legality rules (all deterministic, rank-local knowledge only):
//  * Regrid runs only at a synchronized hierarchy boundary: no temporal
//    boundary history and no pending flux registers (otherwise the plan is
//    empty and `skipped_not_quiescent` is set; never an error).
//  * Passes alternate between refine-only and derefine-only (selected by the
//    caller from the step index). With current neighbour levels known (local
//    + remote interface metadata), simultaneous refinements only raise and
//    simultaneous derefinements only lower levels, so the 2:1 face balance
//    (|level difference| <= 1, periodic wrap included) is preserved across
//    ranks without exchanging decisions.
//  * Refine: any cell requests refinement, level < max_level, every face
//    neighbour has level >= patch level (else `blocked_by_balance`).
//  * Derefine: eight sibling leaf patches (same parent id, same level, all on
//    this rank) all of whose cells allow derefinement (hysteresis), no forced
//    refine request in any child, and no face neighbour of the parent volume
//    at level > child level (else blocked). Octets split across ranks are
//    skipped (`skipped_remote_sibling`) — migration consolidation is the
//    balance runtime's job.
//  * Operations are ordered by ascending patch id. Mandatory refinement is
//    emitted before optional requests and bypasses `max_operations`; only
//    optional operations are truncated by that deterministic resource budget.
//
// Indicators (typed, A03 seam): cell mass, normalized density/pressure jump
// (contact and shock), comoving Truelove-Jeans resolution
// lambda_J,com = c_s sqrt(pi a / (G rho_com)) (the evaluator is fed
// rho_com / a), and explicit forced-refine row requests from physics
// operators (predictive feedback / chemistry / gravity, Checkpoint F).
//
// IDs: new patch and gas-cell ids are drawn from caller-provided windows; the
// runtime derives disjoint rank-striped windows from global maxima so ids are
// globally unique and independent of OpenMP/row order.

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::amr {

enum class AmrRegridPassKind : std::uint8_t {
  kRefineOnly = 0,
  kDerefineOnly = 1,
};

// R1 (A03): physics-aware indicators. MANDATORY criteria (accuracy): the
// existing Jeans and cell-mass criteria and the cooling length. A cell that
// violates a mandatory criterion is refined first and is never derefined
// (hysteresis), independent of any resource budget. OPTIONAL criteria
// (structure tracking, resource-prioritized): the jump indicator, metallicity
// and species interfaces and vorticity; their refinements are scheduled only
// after every mandatory one and are the first to be dropped by the budget.
// All stencils are in-patch face neighbours (one-sided at patch faces), so
// the decision of a patch depends only on that patch's own cells: it is
// identical on any MPI decomposition (patches are never split across ranks)
// and needs no ghost exchange; interfaces between patches are resolved by the
// 2:1 balance rule.
struct PhysicsRefinementIndicators {
  // MANDATORY: refine when the cooling length c_s t_cool spans fewer than
  // cooling_length_cells cells (physical lengths: c_s and t_cool are physical
  // code units, dx_phys = a dx_com). cooling_time_code: one value per dense
  // gas row, |e / (de/dt)_cool| in code time; +inf (or non-positive) = no net
  // cooling. Derefine only if c_s t_cool >= cooling_length_cells dx / hysteresis.
  double cooling_length_cells = 0.0;
  std::span<const double> cooling_time_code{};
  // OPTIONAL: max_faces |Z_a - Z_b| / max(Z_a, Z_b, metallicity_floor), Z =
  // metal mass fraction (dimensionless; 0 disables).
  double metallicity_jump = 0.0;
  double metallicity_floor = 1.0e-6;
  // OPTIONAL: max over species lanes and faces of |x_a - x_b|, x = species
  // mass fraction of the composition gas mass (needs the composition
  // sidecar; 0 disables).
  double species_jump = 0.0;
  // OPTIONAL: |curl v| dx / c_s with comoving gradients and widths (scale-free:
  // equals the physical |omega| dx_phys / c_s). 0 disables.
  double vorticity = 0.0;
};

enum class AmrForceResolutionAction : std::uint8_t {
  kAllowCountedReference = 0,
  kLimitOptionalRejectMandatory = 1,
};

// V02 keeps unlike quantities typed. Candidates are ordered by normalized
// indicator benefit per incremental cell-step, then MPI surface and churn;
// no scalar "universal physical error" is invented. When enabled, ceilings
// are exact (UINT64_MAX means unlimited). Mandatory refinements bypass every
// optional ceiling.
struct OptionalRefinementControl {
  bool enabled = false;
  double benefit_floor = 1.0;
  std::uint64_t memory_budget_bytes = 0U;
  std::uint64_t estimated_bytes_per_cell = 0U;
  std::uint64_t cell_step_budget = 0U;
  std::uint64_t mpi_face_cell_budget = 0U;
  std::uint64_t churn_cell_budget = 0U;
};

struct ProductionAmrRegridPolicy {
  RefinementCriteria criteria{};
  PhysicsRefinementIndicators physics{};
  std::uint8_t max_level = 10;
  // Optional-refinement budget (refine pass); total operation budget for the
  // derefine pass. Mandatory accuracy refinements are never charged here.
  std::size_t max_operations = 64;
  AmrRegridPassKind pass = AmrRegridPassKind::kRefineOnly;
  double scale_factor = 1.0;
  bool periodic = false;
  std::array<double, 3> periodic_box_extent_comoving{0.0, 0.0, 0.0};
  // Sorted, unique dense rows whose owning patch must be refined (physics
  // operator requests, e.g. predictive refinement before scheduled feedback).
  std::span<const std::uint32_t> forced_refine_rows{};
  // CHUI Part 3 Checkpoint E (G01): gas gravitational softening (comoving
  // code length). Refinement is NEVER blocked by it (mandatory Jeans/hydro
  // refinement is not coarsened to fit gravity); refinements whose child cell
  // width falls below it are COUNTED as the unqualified regime "hydro
  // resolution finer than the gravitational force resolution". 0 disables.
  double gravity_softening_code = 0.0;
  AmrForceResolutionAction force_resolution_action =
      AmrForceResolutionAction::kLimitOptionalRejectMandatory;
  OptionalRefinementControl optional_control{};
};

enum class AmrRegridOperationKind : std::uint8_t {
  kRefine = 0,
  kDerefine = 1,
};

struct ProductionAmrRegridOperation {
  AmrRegridOperationKind kind = AmrRegridOperationKind::kRefine;
  // Refine: the leaf patch to refine. Derefine: the reconstructed parent.
  PatchDescriptor patch{};
  std::size_t new_patch_ids = 0;
  std::size_t new_gas_cell_ids = 0;
  double optional_benefit = 0.0;
  std::uint64_t incremental_cell_steps = 0U;
  std::uint64_t mpi_face_cells = 0U;
  std::uint64_t churn_cells = 0U;
  std::uint64_t estimated_memory_bytes = 0U;
};

struct ProductionAmrRegridPlan {
  std::vector<ProductionAmrRegridOperation> operations;
  bool skipped_not_quiescent = false;
  std::size_t refine_requests = 0;
  std::size_t derefine_requests = 0;
  std::size_t blocked_by_balance = 0;
  std::size_t blocked_by_max_level = 0;
  std::size_t skipped_remote_sibling = 0;
  std::size_t truncated_by_budget = 0;
  std::size_t refined_below_gravity_softening = 0;
  std::size_t optional_blocked_by_force_resolution = 0;
  // R1: classification of refinement requests and budget outcome.
  std::size_t mandatory_refine_requests = 0;
  std::size_t optional_refine_requests = 0;
  // Hard resolution requirements are never truncated by max_operations.
  // These counters identify the two remaining reasons a mandatory request
  // cannot be satisfied in the current pass; the runtime reports them and
  // fails closed rather than advancing below a mandatory accuracy criterion.
  std::size_t mandatory_truncated_by_budget = 0;  // retained for schema/telemetry compatibility; always zero
  std::size_t optional_truncated_by_budget = 0;
  // A01/A04/V02: admitted optional operations and their summed proxy costs
  // (cell-steps, MPI interface cells, regrid churn, memory) and the smallest
  // admitted normalized benefit (0 when none was admitted).
  std::size_t optional_admitted_operations = 0;
  std::uint64_t optional_admitted_cell_steps = 0U;
  std::uint64_t optional_admitted_mpi_face_cells = 0U;
  std::uint64_t optional_admitted_churn_cells = 0U;
  std::uint64_t optional_admitted_memory_bytes = 0U;
  double optional_min_admitted_benefit = 0.0;
  std::size_t mandatory_blocked_by_max_level = 0;
  std::size_t mandatory_blocked_by_balance = 0;
  std::size_t required_patch_ids = 0;
  std::size_t required_gas_cell_ids = 0;
};

struct ProductionAmrRegridIdWindow {
  std::uint64_t first_patch_id = 0;
  std::uint64_t patch_id_count = 0;
  std::uint64_t first_gas_cell_id = 0;
  std::uint64_t gas_cell_id_count = 0;
};

struct ProductionAmrRegridSummary {
  std::size_t applied_refinements = 0;
  std::size_t applied_derefinements = 0;
  std::size_t created_gas_cells = 0;
  std::size_t retired_gas_cells = 0;
  double mass_change_code = 0.0;
  double momentum_x_change_code = 0.0;
  double momentum_y_change_code = 0.0;
  double momentum_z_change_code = 0.0;
  double total_energy_change_code = 0.0;
};

// Per-cell indicator used by the planner (exposed for tests): max over the
// in-patch face neighbours of max(|drho| / min(rho), |dP| / min(P)).
[[nodiscard]] double amrRegridJumpIndicator(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    std::span<const std::uint32_t> patch_local_rows,
    std::size_t patch_cell);

[[nodiscard]] ProductionAmrRegridPlan planProductionAmrRegrid(
    const core::SimulationState& state,
    std::span<const PatchDescriptor> remote_neighbor_patches,
    const ProductionAmrRegridPolicy& policy);

// Disjoint per-rank windows above the global maxima. Throws on overflow.
[[nodiscard]] ProductionAmrRegridIdWindow rankStripedAmrRegridIdWindow(
    std::uint64_t global_max_patch_id,
    std::uint64_t global_max_gas_cell_id,
    int rank,
    int world_size,
    std::uint64_t patch_ids_per_rank,
    std::uint64_t gas_cell_ids_per_rank);

// Applies every operation in order. Each operation is an admitted,
// all-or-nothing candidate commit; an exception leaves the state at the last
// fully committed operation (a valid, conservative hierarchy) and propagates
// so the runtime's collective failure handling stops the run.
[[nodiscard]] ProductionAmrRegridSummary applyProductionAmrRegridPlan(
    core::SimulationState& state,
    const ProductionAmrRegridPlan& plan,
    const ProductionAmrRegridIdWindow& window,
    const ProductionAmrHydroOptions& options);

}  // namespace cosmosim::amr
