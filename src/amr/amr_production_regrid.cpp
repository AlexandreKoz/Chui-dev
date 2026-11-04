#include "cosmosim/amr/amr_production_regrid.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

#include "amr/internal/amr_periodic_topology.hpp"
#include "cosmosim/amr/amr_patch_indexing.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"

namespace cosmosim::amr {
namespace {

constexpr double k_tol = 1.0e-10;

[[nodiscard]] std::size_t cellCount(const PatchDescriptor& patch) {
  return checkedPatchCellCount(patch, "production AMR regrid");
}

[[nodiscard]] double minCellWidth(const PatchDescriptor& patch) {
  double width = std::numeric_limits<double>::infinity();
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    width = std::min(width, patch.extent_comov[axis] / static_cast<double>(patch.cell_dims[axis]));
  }
  return width;
}

// Face contact on any of the six sides (periodic minimum image on periodic
// axes), with positive-area transverse overlap.
[[nodiscard]] bool patchesTouch(
    const PatchDescriptor& a, const PatchDescriptor& b, const internal::AmrPeriodicFrame& frame) {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const double a_lo = a.origin_comov[axis];
    const double a_hi = a_lo + a.extent_comov[axis];
    const double b_lo = b.origin_comov[axis];
    const double b_hi = b_lo + b.extent_comov[axis];
    const double scale = std::max({1.0, std::abs(a_hi), std::abs(b_hi)});
    const bool contact =
        std::abs(internal::minImageDifference(frame, axis, a_lo, b_hi)) <= k_tol * scale ||
        std::abs(internal::minImageDifference(frame, axis, a_hi, b_lo)) <= k_tol * scale;
    if (!contact) {
      continue;
    }
    bool overlap = true;
    for (std::size_t other = 0; other < 3U && overlap; ++other) {
      if (other == axis) continue;
      const double lo = std::max(a.origin_comov[other], b.origin_comov[other]);
      const double hi = std::min(a.origin_comov[other] + a.extent_comov[other],
                                 b.origin_comov[other] + b.extent_comov[other]);
      overlap = hi - lo > k_tol * std::max(1.0, std::abs(hi));
    }
    if (overlap) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] bool quiescent(const core::SimulationState& state) {
  return state.amr_temporal_boundary_history.empty() && state.pending_flux_registers.empty();
}

[[nodiscard]] RefinementDecision patchDecision(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    const ProductionAmrRegridPolicy& policy,
    bool& forced) {
  const std::vector<std::uint32_t> rows =
      rowsForPatchLocalCellOrder(state, patch, "production AMR regrid indicator");
  const double width = minCellWidth(patch);
  const double a = policy.scale_factor > 0.0 ? policy.scale_factor : 1.0;
  bool refine_any = false;
  bool derefine_all = true;
  forced = false;
  for (std::size_t cell = 0; cell < rows.size(); ++cell) {
    const std::uint32_t row = rows[cell];
    if (!policy.forced_refine_rows.empty() &&
        std::binary_search(policy.forced_refine_rows.begin(), policy.forced_refine_rows.end(), row)) {
      forced = true;
      refine_any = true;
      derefine_all = false;
      continue;
    }
    const CellMetrics metrics{
        // Comoving Jeans length c_s sqrt(pi a / (G rho_com)): feed rho_com / a.
        .density_code = state.gas_cells.density_code[row] / a,
        .pressure_code = state.gas_cells.pressure_code[row],
        .sound_speed_code = state.gas_cells.sound_speed_code[row],
        .gradient_indicator = amrRegridJumpIndicator(state, patch, rows, cell),
        .particle_count = 0U};
    RefinementCriteria criteria = policy.criteria;
    criteria.use_particle_count = false;
    // Mass criterion on the comoving cell mass rho_com dx^3 (not rho/a).
    const RefinementDecision mass_free = [&] {
      RefinementCriteria no_mass = criteria;
      no_mass.use_mass_threshold = false;
      return RefinementEvaluator::evaluateCell(metrics, width, no_mass);
    }();
    RefinementDecision decision = mass_free;
    if (criteria.use_mass_threshold) {
      const double cell_mass = state.gas_cells.density_code[row] * width * width * width;
      if (cell_mass > criteria.mass_threshold_code) {
        decision = RefinementDecision::kRefine;
      } else if (decision == RefinementDecision::kDerefine &&
                 cell_mass > criteria.mass_threshold_code * criteria.derefine_hysteresis) {
        decision = RefinementDecision::kKeep;
      }
    }
    if (decision == RefinementDecision::kRefine) {
      refine_any = true;
      derefine_all = false;
    } else if (decision == RefinementDecision::kKeep) {
      derefine_all = false;
    }
  }
  if (refine_any) return RefinementDecision::kRefine;
  if (derefine_all) return RefinementDecision::kDerefine;
  return RefinementDecision::kKeep;
}

}  // namespace

double amrRegridJumpIndicator(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    std::span<const std::uint32_t> patch_local_rows,
    std::size_t patch_cell) {
  const std::size_t nx = patch.cell_dims[0];
  const std::size_t ny = patch.cell_dims[1];
  const std::size_t nz = patch.cell_dims[2];
  const std::size_t i = patch_cell % nx;
  const std::size_t j = (patch_cell / nx) % ny;
  const std::size_t k = patch_cell / (nx * ny);
  const std::uint32_t row = patch_local_rows[patch_cell];
  const double rho = state.gas_cells.density_code[row];
  const double pressure = state.gas_cells.pressure_code[row];
  double indicator = 0.0;
  const auto visit = [&](std::size_t ii, std::size_t jj, std::size_t kk) {
    const std::uint32_t other = patch_local_rows[ii + nx * (jj + ny * kk)];
    const double rho_o = state.gas_cells.density_code[other];
    const double p_o = state.gas_cells.pressure_code[other];
    const double rho_min = std::max(std::min(rho, rho_o), 1.0e-300);
    const double p_min = std::max(std::min(pressure, p_o), 1.0e-300);
    indicator = std::max({indicator, std::abs(rho - rho_o) / rho_min, std::abs(pressure - p_o) / p_min});
  };
  if (i > 0U) visit(i - 1U, j, k);
  if (i + 1U < nx) visit(i + 1U, j, k);
  if (j > 0U) visit(i, j - 1U, k);
  if (j + 1U < ny) visit(i, j + 1U, k);
  if (k > 0U) visit(i, j, k - 1U);
  if (k + 1U < nz) visit(i, j, k + 1U);
  return indicator;
}

ProductionAmrRegridPlan planProductionAmrRegrid(
    const core::SimulationState& state,
    std::span<const PatchDescriptor> remote_neighbor_patches,
    const ProductionAmrRegridPolicy& policy) {
  ProductionAmrRegridPlan plan;
  if (!hasProductionAmrHydroCoverage(state)) {
    return plan;
  }
  if (!quiescent(state)) {
    plan.skipped_not_quiescent = true;
    return plan;
  }
  if (!std::is_sorted(policy.forced_refine_rows.begin(), policy.forced_refine_rows.end())) {
    throw std::invalid_argument("production AMR regrid forced_refine_rows must be sorted");
  }
  const internal::AmrPeriodicFrame frame =
      internal::makeAmrPeriodicFrame(policy.periodic, policy.periodic_box_extent_comoving);
  std::vector<PatchDescriptor> local = buildProductionAmrPatchDescriptors(state);
  std::sort(local.begin(), local.end(),
            [](const PatchDescriptor& l, const PatchDescriptor& r) { return l.patch_id < r.patch_id; });
  std::vector<PatchDescriptor> all = local;
  all.insert(all.end(), remote_neighbor_patches.begin(), remote_neighbor_patches.end());

  const auto neighbourLevelsAtLeast = [&](const PatchDescriptor& patch, int min_level,
                                          const std::vector<std::uint64_t>& exclude) {
    for (const PatchDescriptor& other : all) {
      if (other.patch_id == patch.patch_id ||
          std::find(exclude.begin(), exclude.end(), other.patch_id) != exclude.end()) {
        continue;
      }
      if (patchesTouch(patch, other, frame) && static_cast<int>(other.level) < min_level) {
        return false;
      }
    }
    return true;
  };
  const auto neighbourLevelsAtMost = [&](const PatchDescriptor& patch, int max_level,
                                         const std::vector<std::uint64_t>& exclude) {
    for (const PatchDescriptor& other : all) {
      if (other.patch_id == patch.patch_id ||
          std::find(exclude.begin(), exclude.end(), other.patch_id) != exclude.end()) {
        continue;
      }
      if (patchesTouch(patch, other, frame) && static_cast<int>(other.level) > max_level) {
        return false;
      }
    }
    return true;
  };

  if (policy.pass == AmrRegridPassKind::kRefineOnly) {
    for (const PatchDescriptor& patch : local) {
      bool forced = false;
      if (patchDecision(state, patch, policy, forced) != RefinementDecision::kRefine) {
        continue;
      }
      ++plan.refine_requests;
      if (patch.level >= policy.max_level) {
        ++plan.blocked_by_max_level;
        continue;
      }
      if (!neighbourLevelsAtLeast(patch, static_cast<int>(patch.level), {})) {
        ++plan.blocked_by_balance;
        continue;
      }
      if (plan.operations.size() >= policy.max_operations) {
        ++plan.truncated_by_budget;
        continue;
      }
      if (policy.gravity_softening_code > 0.0 && 0.5 * minCellWidth(patch) < policy.gravity_softening_code) {
        ++plan.refined_below_gravity_softening;
      }
      const std::size_t cells = core::checkedSizeMultiply(cellCount(patch), 8U, "regrid refine cells");
      plan.operations.push_back(ProductionAmrRegridOperation{
          .kind = AmrRegridOperationKind::kRefine, .patch = patch,
          .new_patch_ids = 8U, .new_gas_cell_ids = cells});
    }
  } else {
    // Group local leaves by parent id (deterministic: std::map ordered by id).
    std::map<std::uint64_t, std::vector<PatchDescriptor>> octets;
    for (const PatchDescriptor& patch : local) {
      if (patch.parent_patch_id != 0U && patch.level > 0U) {
        octets[patch.parent_patch_id].push_back(patch);
      }
    }
    std::map<std::uint64_t, std::size_t> remote_children;
    for (const PatchDescriptor& remote : remote_neighbor_patches) {
      if (remote.parent_patch_id != 0U) ++remote_children[remote.parent_patch_id];
    }
    for (const auto& [parent_id, children] : octets) {
      bool all_derefine = true;
      for (const PatchDescriptor& child : children) {
        bool forced = false;
        if (patchDecision(state, child, policy, forced) != RefinementDecision::kDerefine) {
          all_derefine = false;
          break;
        }
      }
      if (!all_derefine) {
        continue;
      }
      ++plan.derefine_requests;
      // G2.2: runtime migration co-locates octet families (enforceAmrFamilyColocation),
      // so a split family is transient and is retried after the next migration
      // boundary rather than being a permanent derefinement barrier.
      if (children.size() != 8U || remote_children.contains(parent_id)) {
        ++plan.skipped_remote_sibling;
        continue;
      }
      const std::uint8_t child_level = children.front().level;
      bool consistent = true;
      std::array<double, 3> origin{std::numeric_limits<double>::infinity(),
                                   std::numeric_limits<double>::infinity(),
                                   std::numeric_limits<double>::infinity()};
      for (const PatchDescriptor& child : children) {
        consistent = consistent && child.level == child_level &&
            child.cell_dims == children.front().cell_dims;
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          consistent = consistent && std::abs(child.extent_comov[axis] - children.front().extent_comov[axis]) <=
              k_tol * std::max(1.0, children.front().extent_comov[axis]);
          origin[axis] = std::min(origin[axis], child.origin_comov[axis]);
        }
      }
      if (!consistent) {
        ++plan.skipped_remote_sibling;
        continue;
      }
      PatchDescriptor parent = children.front();
      parent.patch_id = parent_id;
      parent.level = static_cast<std::uint8_t>(child_level - 1U);
      parent.origin_comov = origin;
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        parent.extent_comov[axis] = 2.0 * children.front().extent_comov[axis];
      }
      // Inverse of the refine morton rule child = parent * 8 + octant + 1.
      const auto lowest = std::min_element(children.begin(), children.end(),
          [](const PatchDescriptor& l, const PatchDescriptor& r) { return l.morton_key < r.morton_key; });
      parent.morton_key = lowest->morton_key == 0U ? 0U : (lowest->morton_key - 1U) / 8U;
      // Leaves store only their immediate parent. The grandparent link of the
      // reconstructed parent is recovered from a same-level sibling of the
      // parent (same morton group, local or remote) that is still a leaf. If
      // no such sibling exists (all siblings refined) the link is unknown and
      // stays 0: this parent can then not be derefined further until the
      // link becomes recoverable (documented limitation; never a wrong id).
      parent.parent_patch_id = 0U;
      if (parent.level > 0U && parent.morton_key != 0U) {
        const std::uint64_t group = (parent.morton_key - 1U) / 8U;
        for (const PatchDescriptor& other : all) {
          if (other.level == parent.level && other.morton_key != 0U &&
              (other.morton_key - 1U) / 8U == group && other.parent_patch_id != 0U) {
            parent.parent_patch_id = other.parent_patch_id;
            break;
          }
        }
      }
      std::vector<std::uint64_t> exclude;
      for (const PatchDescriptor& child : children) exclude.push_back(child.patch_id);
      if (!neighbourLevelsAtMost(parent, static_cast<int>(child_level), exclude)) {
        ++plan.blocked_by_balance;
        continue;
      }
      if (plan.operations.size() >= policy.max_operations) {
        ++plan.truncated_by_budget;
        continue;
      }
      plan.operations.push_back(ProductionAmrRegridOperation{
          .kind = AmrRegridOperationKind::kDerefine, .patch = parent,
          .new_patch_ids = 0U, .new_gas_cell_ids = cellCount(parent)});
    }
  }
  for (const ProductionAmrRegridOperation& operation : plan.operations) {
    plan.required_patch_ids = core::checkedSizeAdd(plan.required_patch_ids, operation.new_patch_ids,
                                                   "regrid patch ids");
    plan.required_gas_cell_ids = core::checkedSizeAdd(plan.required_gas_cell_ids, operation.new_gas_cell_ids,
                                                      "regrid gas ids");
  }
  return plan;
}

ProductionAmrRegridIdWindow rankStripedAmrRegridIdWindow(
    std::uint64_t global_max_patch_id,
    std::uint64_t global_max_gas_cell_id,
    int rank,
    int world_size,
    std::uint64_t patch_ids_per_rank,
    std::uint64_t gas_cell_ids_per_rank) {
  if (rank < 0 || world_size <= 0 || rank >= world_size) {
    throw std::invalid_argument("rankStripedAmrRegridIdWindow: invalid rank");
  }
  const auto stripe = [&](std::uint64_t global_max, std::uint64_t per_rank) {
    const std::uint64_t r = static_cast<std::uint64_t>(rank);
    const std::uint64_t ranks = static_cast<std::uint64_t>(world_size);
    if (per_rank == 0U) {
      return global_max + 1U;
    }
    if (per_rank > (std::numeric_limits<std::uint64_t>::max() - global_max - 1U) / ranks) {
      throw std::overflow_error("rankStripedAmrRegridIdWindow: id space exhausted");
    }
    return global_max + 1U + r * per_rank;
  };
  return ProductionAmrRegridIdWindow{
      .first_patch_id = stripe(global_max_patch_id, patch_ids_per_rank),
      .patch_id_count = patch_ids_per_rank,
      .first_gas_cell_id = stripe(global_max_gas_cell_id, gas_cell_ids_per_rank),
      .gas_cell_id_count = gas_cell_ids_per_rank};
}

ProductionAmrRegridSummary applyProductionAmrRegridPlan(
    core::SimulationState& state,
    const ProductionAmrRegridPlan& plan,
    const ProductionAmrRegridIdWindow& window,
    const ProductionAmrHydroOptions& options) {
  ProductionAmrRegridSummary summary;
  if (plan.operations.empty()) {
    return summary;
  }
  if (plan.required_patch_ids > window.patch_id_count ||
      plan.required_gas_cell_ids > window.gas_cell_id_count) {
    throw std::runtime_error("production AMR regrid plan exceeds its identity window (fail closed)");
  }
  std::uint64_t next_patch_id = window.first_patch_id;
  std::uint64_t next_gas_cell_id = window.first_gas_cell_id;
  const auto accumulate = [&summary](const ProductionAmrRegridDiagnostics& d) {
    summary.created_gas_cells += d.created_gas_cell_count;
    summary.retired_gas_cells += d.retired_gas_cell_count;
    summary.mass_change_code += d.conserved_mass_after - d.conserved_mass_before;
    summary.momentum_x_change_code += d.conserved_momentum_x_after - d.conserved_momentum_x_before;
    summary.momentum_y_change_code += d.conserved_momentum_y_after - d.conserved_momentum_y_before;
    summary.momentum_z_change_code += d.conserved_momentum_z_after - d.conserved_momentum_z_before;
    summary.total_energy_change_code += d.conserved_total_energy_after - d.conserved_total_energy_before;
  };
  for (const ProductionAmrRegridOperation& operation : plan.operations) {
    if (operation.kind == AmrRegridOperationKind::kRefine) {
      const auto diagnostics = refineProductionPatchInSimulationState(
          state, operation.patch, next_patch_id, next_gas_cell_id, options);
      next_patch_id += operation.new_patch_ids;
      next_gas_cell_id += operation.new_gas_cell_ids;
      ++summary.applied_refinements;
      accumulate(diagnostics);
    } else {
      const auto diagnostics = derefineProductionPatchInSimulationState(
          state, operation.patch, next_gas_cell_id, options);
      next_gas_cell_id += operation.new_gas_cell_ids;
      ++summary.applied_derefinements;
      accumulate(diagnostics);
    }
  }
  return summary;
}

}  // namespace cosmosim::amr
