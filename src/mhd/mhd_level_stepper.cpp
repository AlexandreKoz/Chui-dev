#include "cosmosim/mhd/mhd_level_stepper.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <tuple>

namespace cosmosim::mhd {
namespace {

MhdFailureContext makeContext(std::string_view owner, std::string detail) {
  MhdFailureContext context{};
  context.owner = std::string(owner);
  context.detail = std::move(detail);
  return context;
}

// Failures that a shorter interval can legitimately cure. Everything else is
// a contract violation that must not masquerade as a numerical retry.
[[nodiscard]] bool retryable(MhdErrorCode code) noexcept {
  return code == MhdErrorCode::kNonPositiveState ||
      code == MhdErrorCode::kStiffnessBudgetExceeded;
}

[[nodiscard]] std::int64_t floorDiv(std::int64_t value, std::int64_t divisor) noexcept {
  std::int64_t quotient = value / divisor;
  if ((value % divisor != 0) && ((value < 0) != (divisor < 0))) {
    --quotient;
  }
  return quotient;
}

[[nodiscard]] bool physicalKind(MhdBoundaryKind kind) noexcept {
  return kind == MhdBoundaryKind::kPeriodic || kind == MhdBoundaryKind::kOpen ||
      kind == MhdBoundaryKind::kReflectingConductingWall;
}

[[nodiscard]] std::size_t edgeIndex(const MhdPatchGeometry& geometry, std::size_t axis,
                                    const std::array<std::int64_t, 3>& e) noexcept {
  const auto i = static_cast<std::size_t>(e[0]);
  const auto j = static_cast<std::size_t>(e[1]);
  const auto k = static_cast<std::size_t>(e[2]);
  return axis == 0U ? geometry.edgeIndexX(i, j, k)
                    : (axis == 1U ? geometry.edgeIndexY(i, j, k) : geometry.edgeIndexZ(i, j, k));
}

[[nodiscard]] std::size_t faceIndex(const MhdPatchGeometry& geometry, std::size_t axis,
                                    const std::array<std::int64_t, 3>& c) noexcept {
  const auto i = static_cast<std::size_t>(c[0]);
  const auto j = static_cast<std::size_t>(c[1]);
  const auto k = static_cast<std::size_t>(c[2]);
  return axis == 0U ? geometry.faceIndexX(i, j, k)
                    : (axis == 1U ? geometry.faceIndexY(i, j, k) : geometry.faceIndexZ(i, j, k));
}

[[nodiscard]] std::span<double> emfLane(MhdEdgeEmfStorage& emf, std::size_t axis) noexcept {
  return axis == 0U ? emf.emfX() : (axis == 1U ? emf.emfY() : emf.emfZ());
}

// In-place copy of an edge field (no reallocation once sized; stage loops must
// not allocate).
void copyEdgeField(MhdEdgeEmfStorage& to, const MhdEdgeEmfStorage& from) {
  const MhdPatchGeometry& geometry = from.geometry();
  if (to.geometry().nx != geometry.nx || to.geometry().ny != geometry.ny ||
      to.geometry().nz != geometry.nz) {
    to.resize(geometry);
  }
  std::copy(from.emfX().begin(), from.emfX().end(), to.emfX().begin());
  std::copy(from.emfY().begin(), from.emfY().end(), to.emfY().begin());
  std::copy(from.emfZ().begin(), from.emfZ().end(), to.emfZ().begin());
}

}  // namespace

// ---------------------------------------------------------------------------
// Construction and topology
// ---------------------------------------------------------------------------

MhdLevelStepper::MhdLevelStepper(MhdNumericsPolicy policy, MhdLevelDomain domain,
                                 std::vector<MhdLevelPatchBinding> patches)
    : m_policy(policy), m_operator(policy), m_domain(domain) {
  m_patches.reserve(patches.size());
  for (MhdLevelPatchBinding& binding : patches) {
    PatchRuntime runtime{};
    runtime.binding = binding;
    m_patches.push_back(std::move(runtime));
  }
  validateTopology();
}

void MhdLevelStepper::validateTopology() {
  if (m_patches.empty()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_level_stepper", "a level needs at least one patch"));
  }
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (!physicalKind(m_domain.lower[axis]) || !physicalKind(m_domain.upper[axis])) {
      throw MhdError(MhdErrorCode::kUnsupportedPhysics,
                     makeContext("mhd_level_stepper",
                                 "level domain boundaries must be periodic, open or conducting "
                                 "wall; remote-rank and shearing boundaries are not implemented"));
    }
    if ((m_domain.lower[axis] == MhdBoundaryKind::kPeriodic) !=
        (m_domain.upper[axis] == MhdBoundaryKind::kPeriodic)) {
      throw MhdError(MhdErrorCode::kUnsupportedPhysics,
                     makeContext("mhd_level_stepper", "one-sided periodic level domain"));
    }
  }
  const MhdLevelPatchBinding& first = m_patches.front().binding;
  m_patch_cells = {static_cast<std::int64_t>(first.geometry.nx),
                   static_cast<std::int64_t>(first.geometry.ny),
                   static_cast<std::int64_t>(first.geometry.nz)};
  std::int64_t expected = 1;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (m_patch_cells[axis] <= 0 || m_domain.cells[axis] <= 0 ||
        m_domain.cells[axis] % m_patch_cells[axis] != 0) {
      throw MhdError(MhdErrorCode::kInvalidField,
                     makeContext("mhd_level_stepper",
                                 "level domain must be an exact multiple of the uniform patch "
                                 "shape"));
    }
    m_grid_count[axis] = m_domain.cells[axis] / m_patch_cells[axis];
    expected *= m_grid_count[axis];
  }
  if (static_cast<std::int64_t>(m_patches.size()) != expected) {
    throw MhdError(MhdErrorCode::kStaleGeometry,
                   makeContext("mhd_level_stepper",
                               "patch set does not tile the level domain exactly once"));
  }
  m_grid_lookup.assign(static_cast<std::size_t>(expected), std::numeric_limits<std::size_t>::max());
  const std::size_t ghost = m_policy.ghostCells();
  for (std::size_t index = 0; index < m_patches.size(); ++index) {
    PatchRuntime& patch = m_patches[index];
    const MhdLevelPatchBinding& binding = patch.binding;
    const MhdPatchGeometry& geometry = binding.geometry;
    if (!geometry.valid() || geometry.nx != first.geometry.nx || geometry.ny != first.geometry.ny ||
        geometry.nz != first.geometry.nz || geometry.dx_comoving != first.geometry.dx_comoving ||
        geometry.dy_comoving != first.geometry.dy_comoving ||
        geometry.dz_comoving != first.geometry.dz_comoving) {
      throw MhdError(MhdErrorCode::kInvalidField,
                     makeContext("mhd_level_stepper",
                                 "level patches must share one valid shape and spacing"));
    }
    if (binding.level != first.level) {
      throw MhdError(MhdErrorCode::kUnsupportedPhysics,
                     makeContext("mhd_level_stepper",
                                 "mixed refinement levels: synchronized multi-level MHD stepping "
                                 "is not implemented in MhdLevelStepper"));
    }
    if (binding.gas == nullptr || binding.gas->size() != geometry.cellCount() ||
        !binding.faces.valid() || binding.faces.patch_id != binding.patch_id ||
        binding.faces.geometry_generation != binding.geometry_generation) {
      throw MhdError(MhdErrorCode::kStaleGeometry,
                     makeContext("mhd_level_stepper",
                                 "patch gas/face views are missing, mis-sized or carry a stale "
                                 "(patch_id, geometry_generation)"));
    }
    if (std::min({geometry.nx, geometry.ny, geometry.nz}) < ghost) {
      throw MhdError(MhdErrorCode::kInvalidField,
                     makeContext("mhd_level_stepper",
                                 "patch extent is smaller than the reconstruction ghost depth"));
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (binding.origin_cell[axis] < 0 || binding.origin_cell[axis] % m_patch_cells[axis] != 0 ||
          binding.origin_cell[axis] >= m_domain.cells[axis]) {
        throw MhdError(MhdErrorCode::kStaleGeometry,
                       makeContext("mhd_level_stepper",
                                   "patch origin is not aligned to the uniform patch lattice"));
      }
      patch.grid[axis] = binding.origin_cell[axis] / m_patch_cells[axis];
      const double spacing = axis == 0U ? geometry.dx_comoving
                                        : (axis == 1U ? geometry.dy_comoving : geometry.dz_comoving);
      const double origin = axis == 0U ? geometry.origin_x_comoving
                                       : (axis == 1U ? geometry.origin_y_comoving
                                                     : geometry.origin_z_comoving);
      const double expected_origin = static_cast<double>(binding.origin_cell[axis]) * spacing;
      if (std::abs(origin - expected_origin) > 1.0e-9 * std::max(1.0, std::abs(expected_origin))) {
        throw MhdError(MhdErrorCode::kStaleGeometry,
                       makeContext("mhd_level_stepper",
                                   "patch comoving origin disagrees with its lattice origin"));
      }
    }
    const std::size_t slot = static_cast<std::size_t>(
        patch.grid[0] + m_grid_count[0] * (patch.grid[1] + m_grid_count[1] * patch.grid[2]));
    if (m_grid_lookup[slot] != std::numeric_limits<std::size_t>::max()) {
      throw MhdError(MhdErrorCode::kStaleGeometry,
                     makeContext("mhd_level_stepper", "two patches occupy one lattice slot"));
    }
    m_grid_lookup[slot] = index;
    for (std::size_t other = 0; other < index; ++other) {
      if (m_patches[other].binding.patch_id == binding.patch_id) {
        throw MhdError(MhdErrorCode::kStaleGeometry,
                       makeContext("mhd_level_stepper", "duplicate patch_id in level"));
      }
    }
    // Patch-side boundary kinds.
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const bool periodic = m_domain.lower[axis] == MhdBoundaryKind::kPeriodic;
      const bool single = m_grid_count[axis] == 1;
      const auto side_kind = [&](bool lower) {
        const bool at_domain_edge =
            lower ? patch.grid[axis] == 0 : patch.grid[axis] == m_grid_count[axis] - 1;
        if (!at_domain_edge) {
          return MhdBoundaryKind::kNeighborPatch;
        }
        if (periodic) {
          return single ? MhdBoundaryKind::kPeriodic : MhdBoundaryKind::kNeighborPatch;
        }
        return lower ? m_domain.lower[axis] : m_domain.upper[axis];
      };
      patch.boundary.lower[axis] = side_kind(true);
      patch.boundary.upper[axis] = side_kind(false);
    }
    validateBoundarySupport(patch.boundary, makeContext("mhd_level_stepper", ""));
  }
}

std::size_t MhdLevelStepper::patchAtGrid(const std::array<std::int64_t, 3>& grid) const {
  const std::size_t slot =
      static_cast<std::size_t>(grid[0] + m_grid_count[0] * (grid[1] + m_grid_count[1] * grid[2]));
  return m_grid_lookup[slot];
}

std::span<double> MhdLevelStepper::faceLane(PatchRuntime& patch, std::size_t axis) noexcept {
  return axis == 0U ? patch.binding.faces.phi_x
                    : (axis == 1U ? patch.binding.faces.phi_y : patch.binding.faces.phi_z);
}

const MhdEdgeEmfStorage& MhdLevelStepper::edgeIntegral(std::size_t index) const {
  return m_patches.at(index).integral;
}

// ---------------------------------------------------------------------------
// Snapshots
// ---------------------------------------------------------------------------

void MhdLevelStepper::snapshotRollback() {
  for (PatchRuntime& patch : m_patches) {
    patch.rollback_gas.copyFrom(*patch.binding.gas);
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const std::span<double> lane = faceLane(patch, axis);
      if (patch.rollback_faces[axis].size() != lane.size()) {
        patch.rollback_faces[axis].resize(lane.size(), 0.0);
      }
      std::copy(lane.begin(), lane.end(), patch.rollback_faces[axis].begin());
    }
    if (m_track_integral) {
      patch.integral.resize(patch.binding.geometry);
      patch.integral.zero();
    }
  }
}

void MhdLevelStepper::restoreRollback() {
  for (PatchRuntime& patch : m_patches) {
    patch.binding.gas->copyFrom(patch.rollback_gas);
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      std::span<double> lane = faceLane(patch, axis);
      std::copy(patch.rollback_faces[axis].begin(), patch.rollback_faces[axis].end(),
                lane.begin());
    }
    if (m_track_integral) {
      patch.integral.zero();
    }
  }
}

void MhdLevelStepper::snapshotStage0(bool gas_all_lanes) {
  for (PatchRuntime& patch : m_patches) {
    if (gas_all_lanes || patch.stage0_gas.size() != patch.binding.gas->size()) {
      patch.stage0_gas.copyFrom(*patch.binding.gas);
    } else {
      const std::span<const double> energy = patch.binding.gas->totalEnergy();
      std::copy(energy.begin(), energy.end(), patch.stage0_gas.totalEnergy().begin());
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const std::span<double> lane = faceLane(patch, axis);
      if (patch.stage0_faces[axis].size() != lane.size()) {
        patch.stage0_faces[axis].resize(lane.size(), 0.0);
      }
      std::copy(lane.begin(), lane.end(), patch.stage0_faces[axis].begin());
    }
    if (m_track_integral) {
      copyEdgeField(patch.stage0_integral, patch.integral);
    }
  }
}

// magnetic := alpha * stage0 + beta * magnetic   (faces and tracked integral)
void MhdLevelStepper::combineMagnetic(PatchRuntime& patch, double alpha, double beta) {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    std::span<double> lane = faceLane(patch, axis);
    const MhdDoubleArray& base = patch.stage0_faces[axis];
    for (std::size_t index = 0; index < lane.size(); ++index) {
      lane[index] = alpha * base[index] + beta * lane[index];
    }
  }
  if (m_track_integral) {
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      std::span<double> lane = emfLane(patch.integral, axis);
      std::span<double> base = emfLane(patch.stage0_integral, axis);
      for (std::size_t index = 0; index < lane.size(); ++index) {
        lane[index] = alpha * base[index] + beta * lane[index];
      }
    }
  }
}

// Phi -= weight * C(emf); the tracked time integral gains weight * emf so that
// Phi_end = Phi_begin - C(I) holds exactly for every linear stage combination.
void MhdLevelStepper::applyCt(PatchRuntime& patch, const MhdEdgeEmfStorage& emf, double weight) {
  applyInductionUpdate(patch.binding.geometry, patch.binding.faces.phi_x,
                       patch.binding.faces.phi_y, patch.binding.faces.phi_z, emf, weight);
  if (m_track_integral) {
    accumulateEdgeEmfIntegral(patch.binding.geometry, patch.integral, emf, weight);
  }
}

// ---------------------------------------------------------------------------
// Stage preparation and halo exchange
// ---------------------------------------------------------------------------

void MhdLevelStepper::prepareStage(double scale_factor) {
  for (PatchRuntime& patch : m_patches) {
    const MhdFaceFluxView view{&patch.binding.geometry,
                               patch.binding.faces.phi_x,
                               patch.binding.faces.phi_y,
                               patch.binding.faces.phi_z,
                               patch.binding.patch_id,
                               patch.binding.geometry_generation,
                               patch.binding.level};
    const std::size_t bad = m_operator.recoverPrimitives(patch.binding.geometry, *patch.binding.gas,
                                                         view, scale_factor, patch.workspace);
    if (bad != 0U) {
      MhdFailureContext failure = makeContext(
          "mhd_level_stepper",
          std::to_string(bad) + " inadmissible cell(s) at stage preparation (no floor applied)");
      failure.patch_id = patch.binding.patch_id;
      failure.geometry_generation = patch.binding.geometry_generation;
      failure.scale_factor = scale_factor;
      throw MhdError(MhdErrorCode::kNonPositiveState, std::move(failure));
    }
  }
  fillNeighbourGhosts();
  for (PatchRuntime& patch : m_patches) {
    MhdStepOperator::fillPhysicalGhosts(patch.boundary, patch.workspace);
  }
}

void MhdLevelStepper::fillNeighbourGhosts() {
  for (PatchRuntime& patch : m_patches) {
    MhdPrimitivePatch& primitives = patch.workspace.primitives;
    const auto g = static_cast<std::int64_t>(primitives.ghost);
    const std::array<std::int64_t, 3> n = m_patch_cells;
    // Cells.
    for (std::int64_t k = -g; k < n[2] + g; ++k) {
      for (std::int64_t j = -g; j < n[1] + g; ++j) {
        for (std::int64_t i = -g; i < n[0] + g; ++i) {
          const std::array<std::int64_t, 3> local{i, j, k};
          bool ghost_cell = false;
          bool neighbour_only = true;
          std::array<std::int64_t, 3> global{};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            global[axis] = patch.binding.origin_cell[axis] + local[axis];
            if (local[axis] >= 0 && local[axis] < n[axis]) {
              continue;
            }
            ghost_cell = true;
            const MhdBoundaryKind kind =
                local[axis] < 0 ? patch.boundary.lower[axis] : patch.boundary.upper[axis];
            if (kind != MhdBoundaryKind::kNeighborPatch) {
              neighbour_only = false;
              break;
            }
            const std::int64_t extent = m_domain.cells[axis];
            global[axis] = ((global[axis] % extent) + extent) % extent;
          }
          if (!ghost_cell || !neighbour_only) {
            continue;
          }
          std::array<std::int64_t, 3> grid{};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            grid[axis] = floorDiv(global[axis], m_patch_cells[axis]);
          }
          const PatchRuntime& source = m_patches[patchAtGrid(grid)];
          MhdIndex3 source_local{};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            source_local[axis] =
                static_cast<std::ptrdiff_t>(global[axis] - source.binding.origin_cell[axis]);
          }
          primitives.store({static_cast<std::ptrdiff_t>(i), static_cast<std::ptrdiff_t>(j),
                            static_cast<std::ptrdiff_t>(k)},
                           source.workspace.primitives.load(source_local));
        }
      }
    }
    // Staggered faces (tangential ghost rows only).
    MhdStaggeredFieldPatch& staggered = patch.workspace.staggered;
    const auto gt = static_cast<std::int64_t>(staggered.tangential_ghost);
    for (std::size_t normal = 0; normal < 3U; ++normal) {
      const MhdAxis normal_axis = static_cast<MhdAxis>(normal);
      std::array<std::int64_t, 3> lo{};
      std::array<std::int64_t, 3> hi{};
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        lo[axis] = axis == normal ? 0 : -gt;
        hi[axis] = axis == normal ? n[axis] + 1 : n[axis] + gt;
      }
      for (std::int64_t k = lo[2]; k < hi[2]; ++k) {
        for (std::int64_t j = lo[1]; j < hi[1]; ++j) {
          for (std::int64_t i = lo[0]; i < hi[0]; ++i) {
            const std::array<std::int64_t, 3> local{i, j, k};
            bool ghost_face = false;
            bool neighbour_only = true;
            std::array<std::int64_t, 3> global{};
            for (std::size_t axis = 0; axis < 3U; ++axis) {
              global[axis] = patch.binding.origin_cell[axis] + local[axis];
              if (axis == normal || (local[axis] >= 0 && local[axis] < n[axis])) {
                continue;
              }
              ghost_face = true;
              const MhdBoundaryKind kind =
                  local[axis] < 0 ? patch.boundary.lower[axis] : patch.boundary.upper[axis];
              if (kind != MhdBoundaryKind::kNeighborPatch) {
                neighbour_only = false;
                break;
              }
              const std::int64_t extent = m_domain.cells[axis];
              global[axis] = ((global[axis] % extent) + extent) % extent;
            }
            if (!ghost_face || !neighbour_only) {
              continue;
            }
            std::array<std::int64_t, 3> grid{};
            for (std::size_t axis = 0; axis < 3U; ++axis) {
              if (axis == normal) {
                // Face node coordinate: pick the patch containing cell G (or the
                // last patch for the upper domain face); shared faces are
                // reconciled so either holder carries the same value.
                grid[axis] = std::min(floorDiv(global[axis], m_patch_cells[axis]),
                                      m_grid_count[axis] - 1);
              } else {
                grid[axis] = floorDiv(global[axis], m_patch_cells[axis]);
              }
            }
            const PatchRuntime& source = m_patches[patchAtGrid(grid)];
            MhdIndex3 source_local{};
            for (std::size_t axis = 0; axis < 3U; ++axis) {
              source_local[axis] =
                  static_cast<std::ptrdiff_t>(global[axis] - source.binding.origin_cell[axis]);
            }
            staggered.at(normal_axis, {static_cast<std::ptrdiff_t>(i), static_cast<std::ptrdiff_t>(j),
                                       static_cast<std::ptrdiff_t>(k)}) =
                source.workspace.staggered.at(normal_axis, source_local);
          }
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Owner reconciliation of shared faces and edges.
//
// A logical node coordinate G on a patch-lattice boundary is held by the patch
// above it (local 0) and the patch below it (local n). The owner is the holder
// with the lowest (patch_id, local coordinates); its value is never written by
// reconciliation, so the result is independent of iteration order.
// ---------------------------------------------------------------------------

namespace {

struct Holder {
  std::uint64_t patch_id = 0U;
  std::size_t patch_index = 0U;
  std::array<std::int64_t, 3> local{};
};

[[nodiscard]] bool holderBefore(const Holder& lhs, const Holder& rhs) noexcept {
  return std::tie(lhs.patch_id, lhs.local) < std::tie(rhs.patch_id, rhs.local);
}

}  // namespace

void MhdLevelStepper::reconcileEdges() {
  // Owner selection works on canonical (wrapped) global node coordinates.
  const auto owner = [&](std::size_t c, const std::array<std::int64_t, 3>& global) {
    std::array<std::array<std::pair<std::int64_t, std::int64_t>, 2>, 3> options{};
    std::array<int, 3> option_count{};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const std::int64_t extent = m_domain.cells[axis];
      const bool periodic = m_domain.lower[axis] == MhdBoundaryKind::kPeriodic;
      std::int64_t g = global[axis];
      if (axis != c && periodic) {
        g = ((g % extent) + extent) % extent;
      }
      if (axis == c) {
        options[axis][0] = {floorDiv(g, m_patch_cells[axis]),
                            g - floorDiv(g, m_patch_cells[axis]) * m_patch_cells[axis]};
        option_count[axis] = 1;
        continue;
      }
      const std::int64_t gi = floorDiv(g, m_patch_cells[axis]);
      const std::int64_t remainder = g - gi * m_patch_cells[axis];
      int count = 0;
      if (remainder != 0) {
        options[axis][count++] = {gi, remainder};
      } else {
        if (gi < m_grid_count[axis]) {
          options[axis][count++] = {gi, 0};
        }
        std::int64_t below = gi - 1;
        if (below < 0 && periodic) {
          below += m_grid_count[axis];
        }
        if (below >= 0 && below < m_grid_count[axis]) {
          options[axis][count++] = {below, m_patch_cells[axis]};
        }
      }
      option_count[axis] = count;
    }
    Holder best{};
    bool have = false;
    for (int ox = 0; ox < option_count[0]; ++ox) {
      for (int oy = 0; oy < option_count[1]; ++oy) {
        for (int oz = 0; oz < option_count[2]; ++oz) {
          const std::array<std::int64_t, 3> grid{options[0][ox].first, options[1][oy].first,
                                                 options[2][oz].first};
          Holder candidate{};
          candidate.patch_index = patchAtGrid(grid);
          candidate.patch_id = m_patches[candidate.patch_index].binding.patch_id;
          candidate.local = {options[0][ox].second, options[1][oy].second, options[2][oz].second};
          if (!have || holderBefore(candidate, best)) {
            best = candidate;
            have = true;
          }
        }
      }
    }
    return best;
  };
  for (std::size_t index = 0; index < m_patches.size(); ++index) {
    PatchRuntime& patch = m_patches[index];
    const MhdPatchGeometry& geometry = patch.binding.geometry;
    for (std::size_t c = 0; c < 3U; ++c) {
      const std::size_t a = (c + 1U) % 3U;
      const std::size_t b = (c + 2U) % 3U;
      std::span<double> lane = emfLane(patch.workspace.emf, c);
      for (std::int64_t kc = 0; kc < m_patch_cells[c]; ++kc) {
        for (std::int64_t nb = 0; nb <= m_patch_cells[b]; ++nb) {
          for (std::int64_t na = 0; na <= m_patch_cells[a]; ++na) {
            const bool on_boundary = na == 0 || na == m_patch_cells[a] || nb == 0 ||
                nb == m_patch_cells[b];
            if (!on_boundary) {
              continue;
            }
            std::array<std::int64_t, 3> local{};
            local[c] = kc;
            local[a] = na;
            local[b] = nb;
            std::array<std::int64_t, 3> global{};
            for (std::size_t axis = 0; axis < 3U; ++axis) {
              global[axis] = patch.binding.origin_cell[axis] + local[axis];
            }
            const Holder holder = owner(c, global);
            if (holder.patch_index == index && holder.local == local) {
              continue;
            }
            PatchRuntime& source = m_patches[holder.patch_index];
            lane[edgeIndex(geometry, c, local)] =
                emfLane(source.workspace.emf, c)[edgeIndex(source.binding.geometry, c, holder.local)];
          }
        }
      }
    }
  }
}

namespace {

// Holders of a face whose normal-axis node coordinate is G_n (global).
struct FaceHolders {
  std::array<std::pair<std::int64_t, std::int64_t>, 2> options{};  // (grid, local)
  int count = 0;
};

}  // namespace

double MhdLevelStepper::reconcileFaces() {
  double max_difference = 0.0;
  for (std::size_t index = 0; index < m_patches.size(); ++index) {
    PatchRuntime& patch = m_patches[index];
    for (std::size_t normal = 0; normal < 3U; ++normal) {
      const std::size_t t1 = (normal + 1U) % 3U;
      const std::size_t t2 = (normal + 2U) % 3U;
      const bool periodic = m_domain.lower[normal] == MhdBoundaryKind::kPeriodic;
      std::span<double> lane = faceLane(patch, normal);
      for (int side = 0; side < 2; ++side) {
        const std::int64_t local_n = side == 0 ? 0 : m_patch_cells[normal];
        std::int64_t g = patch.binding.origin_cell[normal] + local_n;
        if (periodic) {
          g %= m_domain.cells[normal];
        }
        const std::int64_t gi = g / m_patch_cells[normal];
        Holder best{};
        bool have = false;
        const auto consider = [&](std::int64_t grid_n, std::int64_t holder_local) {
          std::array<std::int64_t, 3> grid = patch.grid;
          grid[normal] = grid_n;
          Holder candidate{};
          candidate.patch_index = patchAtGrid(grid);
          candidate.patch_id = m_patches[candidate.patch_index].binding.patch_id;
          candidate.local = {0, 0, 0};
          candidate.local[normal] = holder_local;
          if (!have || holderBefore(candidate, best)) {
            best = candidate;
            have = true;
          }
        };
        if (gi < m_grid_count[normal]) {
          consider(gi, 0);
        }
        std::int64_t below = gi - 1;
        if (below < 0 && periodic) {
          below += m_grid_count[normal];
        }
        if (below >= 0 && below < m_grid_count[normal]) {
          consider(below, m_patch_cells[normal]);
        }
        if (!have || (best.patch_index == index && best.local[normal] == local_n)) {
          continue;
        }
        PatchRuntime& source = m_patches[best.patch_index];
        std::span<double> source_lane = faceLane(source, normal);
        for (std::int64_t c2 = 0; c2 < m_patch_cells[t2]; ++c2) {
          for (std::int64_t c1 = 0; c1 < m_patch_cells[t1]; ++c1) {
            std::array<std::int64_t, 3> local{};
            local[normal] = local_n;
            local[t1] = c1;
            local[t2] = c2;
            std::array<std::int64_t, 3> source_local = local;
            source_local[normal] = best.local[normal];
            double& target = lane[faceIndex(patch.binding.geometry, normal, local)];
            const double value = source_lane[faceIndex(source.binding.geometry, normal, source_local)];
            max_difference = std::max(max_difference, std::abs(target - value));
            target = value;
          }
        }
      }
    }
  }
  return max_difference;
}

void MhdLevelStepper::reconcileBoundaryGasFluxes() {
  for (std::size_t index = 0; index < m_patches.size(); ++index) {
    PatchRuntime& patch = m_patches[index];
    for (std::size_t normal = 0; normal < 3U; ++normal) {
      const bool periodic = m_domain.lower[normal] == MhdBoundaryKind::kPeriodic;
      for (std::size_t side = 0; side < 2U; ++side) {
        const std::int64_t local_n = side == 0U ? 0 : m_patch_cells[normal];
        std::int64_t g = patch.binding.origin_cell[normal] + local_n;
        if (periodic) {
          g %= m_domain.cells[normal];
        }
        const std::int64_t gi = g / m_patch_cells[normal];
        Holder best{};
        bool have = false;
        const auto consider = [&](std::int64_t grid_n, std::int64_t holder_local) {
          std::array<std::int64_t, 3> grid = patch.grid;
          grid[normal] = grid_n;
          Holder candidate{};
          candidate.patch_index = patchAtGrid(grid);
          candidate.patch_id = m_patches[candidate.patch_index].binding.patch_id;
          candidate.local = {0, 0, 0};
          candidate.local[normal] = holder_local;
          if (!have || holderBefore(candidate, best)) {
            best = candidate;
            have = true;
          }
        };
        if (gi < m_grid_count[normal]) {
          consider(gi, 0);
        }
        std::int64_t below = gi - 1;
        if (below < 0 && periodic) {
          below += m_grid_count[normal];
        }
        if (below >= 0 && below < m_grid_count[normal]) {
          consider(below, m_patch_cells[normal]);
        }
        if (!have || (best.patch_index == index && best.local[normal] == local_n)) {
          continue;
        }
        const std::size_t source_side = best.local[normal] == 0 ? 0U : 1U;
        const PatchRuntime& source = m_patches[best.patch_index];
        const std::size_t count = patch.workspace.boundary_fluxes.faceCount(normal);
        const std::size_t n_t1 = static_cast<std::size_t>(m_patch_cells[(normal + 1U) % 3U]);
        for (std::size_t face = 0; face < count; ++face) {
          const std::size_t c1 = face % n_t1;
          const std::size_t c2 = face / n_t1;
          const double* from = source.workspace.boundary_fluxes.at(normal, source_side, c1, c2);
          double* to = patch.workspace.boundary_fluxes.at(normal, side, c1, c2);
          std::copy(from, from + MhdBoundaryFaceFluxes::k_values_per_face, to);
        }
        // Owner's own energy record for the non-ideal path is reconciled by
        // reconcileBoundaryEnergyFluxes.
      }
    }
  }
}

void MhdLevelStepper::reconcileBoundaryEnergyFluxes() {
  for (std::size_t index = 0; index < m_patches.size(); ++index) {
    PatchRuntime& patch = m_patches[index];
    for (std::size_t normal = 0; normal < 3U; ++normal) {
      const bool periodic = m_domain.lower[normal] == MhdBoundaryKind::kPeriodic;
      for (std::size_t side = 0; side < 2U; ++side) {
        const std::int64_t local_n = side == 0U ? 0 : m_patch_cells[normal];
        std::int64_t g = patch.binding.origin_cell[normal] + local_n;
        if (periodic) {
          g %= m_domain.cells[normal];
        }
        const std::int64_t gi = g / m_patch_cells[normal];
        Holder best{};
        bool have = false;
        const auto consider = [&](std::int64_t grid_n, std::int64_t holder_local) {
          std::array<std::int64_t, 3> grid = patch.grid;
          grid[normal] = grid_n;
          Holder candidate{};
          candidate.patch_index = patchAtGrid(grid);
          candidate.patch_id = m_patches[candidate.patch_index].binding.patch_id;
          candidate.local = {0, 0, 0};
          candidate.local[normal] = holder_local;
          if (!have || holderBefore(candidate, best)) {
            best = candidate;
            have = true;
          }
        };
        if (gi < m_grid_count[normal]) {
          consider(gi, 0);
        }
        std::int64_t below = gi - 1;
        if (below < 0 && periodic) {
          below += m_grid_count[normal];
        }
        if (below >= 0 && below < m_grid_count[normal]) {
          consider(below, m_patch_cells[normal]);
        }
        if (!have || (best.patch_index == index && best.local[normal] == local_n)) {
          continue;
        }
        const std::size_t source_side = best.local[normal] == 0 ? 0U : 1U;
        const PatchRuntime& source = m_patches[best.patch_index];
        const MhdDoubleArray& from = source.nonideal.boundary_energy.values[normal][source_side];
        MhdDoubleArray& to = patch.nonideal.boundary_energy.values[normal][side];
        std::copy(from.begin(), from.end(), to.begin());
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Ideal SSPRK2 (operator I)
// ---------------------------------------------------------------------------

void MhdLevelStepper::idealSsprk2(const MhdStageContext& interval, MhdLevelAdvanceResult& result) {
  const double dt = interval.dt_code;
  const bool expansion = interval.comoving_coordinates && interval.apply_expansion_source;
  snapshotStage0(true);
  // Stage 1: Y1 = Y0 + dt L(Y0, t_n).
  prepareStage(interval.scale_factor);
  for (PatchRuntime& patch : m_patches) {
    m_operator.evaluateIdealStage(patch.binding.geometry, patch.boundary, interval.scale_factor,
                                  interval.hubble_rate_code, expansion, dt, patch.workspace,
                                  *patch.binding.gas, result.ideal);
  }
  ++result.ideal_stage_evaluations_attempted;
  reconcileBoundaryGasFluxes();
  reconcileEdges();
  for (PatchRuntime& patch : m_patches) {
    MhdStepOperator::applyBoundaryFaceFluxes(patch.binding.geometry,
                                             patch.workspace.boundary_fluxes,
                                             interval.scale_factor, dt, *patch.binding.gas);
    applyCt(patch, patch.workspace.emf, dt);
  }
  reconcileFaces();
  // Stage 2: Y = Y0/2 + (Y1 + dt L(Y1, t_n + dt))/2.
  prepareStage(interval.scale_factor_end);
  for (PatchRuntime& patch : m_patches) {
    patch.binding.gas->assignLinearCombination(0.5, patch.stage0_gas, 0.5, *patch.binding.gas);
    combineMagnetic(patch, 0.5, 0.5);
    m_operator.evaluateIdealStage(patch.binding.geometry, patch.boundary,
                                  interval.scale_factor_end, interval.hubble_rate_end_code,
                                  expansion, 0.5 * dt, patch.workspace, *patch.binding.gas,
                                  result.ideal);
  }
  ++result.ideal_stage_evaluations_attempted;
  reconcileBoundaryGasFluxes();
  reconcileEdges();
  for (PatchRuntime& patch : m_patches) {
    MhdStepOperator::applyBoundaryFaceFluxes(patch.binding.geometry,
                                             patch.workspace.boundary_fluxes,
                                             interval.scale_factor_end, 0.5 * dt,
                                             *patch.binding.gas);
    applyCt(patch, patch.workspace.emf, 0.5 * dt);
  }
  reconcileFaces();
}

// ---------------------------------------------------------------------------
// Non-ideal operators D and H
// ---------------------------------------------------------------------------

double MhdLevelStepper::evaluateNonIdealRate(MhdNonIdealOperatorKind kind, double scale_factor,
                                             MhdLevelAdvanceResult& result) {
  prepareStage(scale_factor);
  double stable = std::numeric_limits<double>::infinity();
  for (PatchRuntime& patch : m_patches) {
    m_nonideal->evaluateCoefficients(patch.binding.geometry, patch.workspace.primitives,
                                     scale_factor, patch.nonideal, result.transport);
    stable = std::min(stable, kind == MhdNonIdealOperatorKind::kDiffusion
                                  ? m_nonideal->diffusionStableTimeStep(patch.binding.geometry,
                                                                        patch.nonideal)
                                  : m_nonideal->hallStableTimeStep(patch.binding.geometry,
                                                                   patch.nonideal));
    m_nonideal->computeEmf(kind, patch.binding.geometry, patch.boundary, patch.workspace.staggered,
                           scale_factor, patch.nonideal, patch.workspace.emf, result.transport);
  }
  reconcileEdges();
  for (PatchRuntime& patch : m_patches) {
    std::fill(patch.nonideal.energy_rate_stage.begin(), patch.nonideal.energy_rate_stage.end(), 0.0);
    m_nonideal->applyEnergyFlux(patch.binding.geometry, patch.workspace.primitives,
                                patch.workspace.emf, scale_factor, 1.0, patch.nonideal,
                                patch.nonideal.energy_rate_stage.span());
  }
  reconcileBoundaryEnergyFluxes();
  for (PatchRuntime& patch : m_patches) {
    applyBoundaryEnergyFluxes(patch.binding.geometry, patch.nonideal.boundary_energy, 1.0,
                              patch.nonideal.energy_rate_stage.span());
  }
  return stable;
}

void MhdLevelStepper::explicitNonIdealSubstep(MhdNonIdealOperatorKind kind,
                                              const MhdStageContext& sub,
                                              MhdLevelAdvanceResult& result) {
  const double h = sub.dt_code;
  const auto require_stable = [&](double stable) {
    if (h > stable * (1.0 + 1.0e-12)) {
      MhdFailureContext failure = makeContext(
          kind == MhdNonIdealOperatorKind::kDiffusion ? "mhd_diffusion" : "mhd_hall",
          "stage coefficients grew beyond the forecast explicit stability bound");
      failure.time_code = sub.time_code;
      failure.scale_factor = sub.scale_factor;
      throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
    }
  };
  snapshotStage0(false);
  require_stable(evaluateNonIdealRate(kind, sub.scale_factor, result));
  for (PatchRuntime& patch : m_patches) {
    std::span<double> energy = patch.binding.gas->totalEnergy();
    for (std::size_t cell = 0; cell < energy.size(); ++cell) {
      energy[cell] += h * patch.nonideal.energy_rate_stage[cell];
    }
    applyCt(patch, patch.workspace.emf, h);
  }
  reconcileFaces();
  require_stable(evaluateNonIdealRate(kind, sub.scale_factor_end, result));
  for (PatchRuntime& patch : m_patches) {
    std::span<double> energy = patch.binding.gas->totalEnergy();
    const std::span<const double> energy0 = patch.stage0_gas.totalEnergy();
    for (std::size_t cell = 0; cell < energy.size(); ++cell) {
      energy[cell] = 0.5 * energy0[cell] +
          0.5 * (energy[cell] + h * patch.nonideal.energy_rate_stage[cell]);
    }
    combineMagnetic(patch, 0.5, 0.5);
    applyCt(patch, patch.workspace.emf, 0.5 * h);
  }
  reconcileFaces();
}

void MhdLevelStepper::rkl2Cycle(const MhdStageContext& sub, std::int32_t stages,
                                MhdLevelAdvanceResult& result) {
  const MhdRkl2Coefficients coefficients = makeRkl2Coefficients(stages);
  const double tau = sub.dt_code;
  // All stages of one cycle use the cycle-midpoint scale factor for the
  // comoving coefficient conversion (second order in the slowly varying a);
  // physical properties are re-evaluated from every stage state.
  const double a_mid = 0.5 * (sub.scale_factor + sub.scale_factor_end);
  const double bound = rkl2StabilityFactor(stages);
  const auto require_stable = [&](double stable) {
    if (tau > stable * bound * (1.0 + 1.0e-12)) {
      MhdFailureContext failure =
          makeContext("mhd_rkl2", "stage coefficients grew beyond the RKL2 stability interval");
      failure.predicted_stages = static_cast<std::uint64_t>(stages);
      failure.time_code = sub.time_code;
      throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
    }
  };
  for (PatchRuntime& patch : m_patches) {
    patch.nonideal.resize(patch.binding.geometry, true);
  }
  snapshotStage0(false);  // Y_0 (energy + faces + integral)
  require_stable(evaluateNonIdealRate(MhdNonIdealOperatorKind::kDiffusion, a_mid, result));
  for (PatchRuntime& patch : m_patches) {
    // Keep L(Y_0) and Y_{j-2} := Y_0.
    copyEdgeField(patch.nonideal.emf_y0, patch.workspace.emf);
    std::copy(patch.nonideal.energy_rate_stage.begin(), patch.nonideal.energy_rate_stage.end(),
              patch.nonideal.energy_rate_y0.begin());
    const std::span<const double> energy0 = patch.stage0_gas.totalEnergy();
    std::copy(energy0.begin(), energy0.end(), patch.nonideal.yjm2_energy.begin());
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      std::copy(patch.stage0_faces[axis].begin(), patch.stage0_faces[axis].end(),
                patch.nonideal.yjm2_faces[axis].begin());
    }
    if (m_track_integral) {
      copyEdgeField(patch.yjm2_integral, patch.stage0_integral);
    }
    // Y_1 = Y_0 + mu~_1 tau L(Y_0).
    std::span<double> energy = patch.binding.gas->totalEnergy();
    for (std::size_t cell = 0; cell < energy.size(); ++cell) {
      energy[cell] += coefficients.mu_tilde[1] * tau * patch.nonideal.energy_rate_y0[cell];
    }
    applyCt(patch, patch.nonideal.emf_y0, coefficients.mu_tilde[1] * tau);
  }
  reconcileFaces();
  ++result.transport.rkl_stage_count;
  for (std::int32_t j = 2; j <= stages; ++j) {
    const auto index = static_cast<std::size_t>(j);
    const double mu = coefficients.mu[index];
    const double nu = coefficients.nu[index];
    const double y0_weight = 1.0 - mu - nu;
    const double mu_tilde_tau = coefficients.mu_tilde[index] * tau;
    const double gamma_tilde_tau = coefficients.gamma_tilde[index] * tau;
    require_stable(evaluateNonIdealRate(MhdNonIdealOperatorKind::kDiffusion, a_mid, result));
    for (PatchRuntime& patch : m_patches) {
      // Energy: Y_j = mu Y_{j-1} + nu Y_{j-2} + (1-mu-nu) Y_0
      //             + mu~ tau L(Y_{j-1}) + gamma~ tau L(Y_0); rotate Y_{j-2}.
      std::span<double> energy = patch.binding.gas->totalEnergy();
      const std::span<const double> energy0 = patch.stage0_gas.totalEnergy();
      for (std::size_t cell = 0; cell < energy.size(); ++cell) {
        const double previous = energy[cell];
        energy[cell] = mu * previous + nu * patch.nonideal.yjm2_energy[cell] +
            y0_weight * energy0[cell] + mu_tilde_tau * patch.nonideal.energy_rate_stage[cell] +
            gamma_tilde_tau * patch.nonideal.energy_rate_y0[cell];
        patch.nonideal.yjm2_energy[cell] = previous;
      }
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        std::span<double> lane = faceLane(patch, axis);
        MhdDoubleArray& yjm2 = patch.nonideal.yjm2_faces[axis];
        const MhdDoubleArray& y0 = patch.stage0_faces[axis];
        for (std::size_t face = 0; face < lane.size(); ++face) {
          const double previous = lane[face];
          lane[face] = mu * previous + nu * yjm2[face] + y0_weight * y0[face];
          yjm2[face] = previous;
        }
      }
      if (m_track_integral) {
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          std::span<double> lane = emfLane(patch.integral, axis);
          std::span<double> yjm2 = emfLane(patch.yjm2_integral, axis);
          std::span<double> y0 = emfLane(patch.stage0_integral, axis);
          for (std::size_t edge = 0; edge < lane.size(); ++edge) {
            const double previous = lane[edge];
            lane[edge] = mu * previous + nu * yjm2[edge] + y0_weight * y0[edge];
            yjm2[edge] = previous;
          }
        }
      }
      // Magnetic parts of the L terms (L_Phi = -C(E)).
      applyCt(patch, patch.workspace.emf, mu_tilde_tau);
      applyCt(patch, patch.nonideal.emf_y0, gamma_tilde_tau);
    }
    reconcileFaces();
    ++result.transport.rkl_stage_count;
  }
  ++result.transport.rkl_cycle_count;
}

void MhdLevelStepper::diffusionHalf(const MhdStageContext& sub, MhdLevelAdvanceResult& result) {
  if (!m_nonideal->diffusionEnabled() || !(sub.dt_code > 0.0)) {
    return;
  }
  const double stable = m_forecast_dt_diffusion;
  if (!std::isfinite(stable)) {
    return;  // all diffusion coefficients exactly zero: exact no-op
  }
  const double ratio = sub.dt_code / stable;
  if (ratio <= static_cast<double>(m_nonideal->limits().rkl_switch_ratio)) {
    // Explicit SSPRK2 with the least number of equal substeps h <= dt_D.
    const auto substeps = static_cast<std::int32_t>(std::max(1.0, std::ceil(ratio * (1.0 - 1.0e-12))));
    for (std::int32_t step = 0; step < substeps; ++step) {
      const double begin = static_cast<double>(step) / substeps;
      const double end = static_cast<double>(step + 1) / substeps;
      explicitNonIdealSubstep(MhdNonIdealOperatorKind::kDiffusion,
                              subIntervalContext(sub, begin, end), result);
      ++result.transport.explicit_diffusion_substeps;
    }
    return;
  }
  const MhdRkl2Plan plan = planRkl2(ratio, m_nonideal->limits().rkl_max_stages);
  if (m_interval_rkl_stages + plan.totalStages() > m_nonideal->limits().rkl_stage_cost_cap) {
    MhdFailureContext failure =
        makeContext("mhd_rkl2", "RKL2 stage forecast exceeds the per-interval cost cap");
    failure.predicted_stages = static_cast<std::uint64_t>(m_interval_rkl_stages + plan.totalStages());
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
  }
  m_interval_rkl_stages += plan.totalStages();
  for (std::int32_t cycle = 0; cycle < plan.cycles; ++cycle) {
    const double begin = static_cast<double>(cycle) / plan.cycles;
    const double end = static_cast<double>(cycle + 1) / plan.cycles;
    rkl2Cycle(subIntervalContext(sub, begin, end), plan.stages_per_cycle, result);
  }
}

void MhdLevelStepper::hallHalf(const MhdStageContext& sub, MhdLevelAdvanceResult& result) {
  if (!m_nonideal->hallEnabled() || !(sub.dt_code > 0.0)) {
    return;
  }
  const double stable = m_forecast_dt_hall;
  if (!std::isfinite(stable)) {
    return;  // eta_H exactly zero everywhere: exact no-op
  }
  const double ratio = sub.dt_code / stable;
  if (!std::isfinite(ratio) || ratio > static_cast<double>(std::numeric_limits<std::int32_t>::max())) {
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded,
                   makeContext("mhd_hall", "Hall subcycle forecast is not representable"));
  }
  const auto subcycles = static_cast<std::int64_t>(std::max(1.0, std::ceil(ratio * (1.0 - 1.0e-12))));
  if (m_interval_hall_subcycles + subcycles > m_nonideal->limits().hall_subcycle_cost_cap) {
    MhdFailureContext failure =
        makeContext("mhd_hall", "Hall subcycle forecast exceeds the per-interval cost cap");
    failure.predicted_stages = static_cast<std::uint64_t>(m_interval_hall_subcycles + subcycles);
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
  }
  m_interval_hall_subcycles += subcycles;
  for (std::int64_t step = 0; step < subcycles; ++step) {
    const double begin = static_cast<double>(step) / static_cast<double>(subcycles);
    const double end = static_cast<double>(step + 1) / static_cast<double>(subcycles);
    explicitNonIdealSubstep(MhdNonIdealOperatorKind::kHall, subIntervalContext(sub, begin, end),
                            result);
    ++result.transport.hall_subcycle_count;
  }
}

// Pre-mutation stiffness forecast at the interval-start state: the explicit
// bounds and the per-interval Hall/RKL cost forecasts for both halves.
void MhdLevelStepper::forecastNonIdeal(const MhdStageContext& interval,
                                       MhdLevelAdvanceResult& result) {
  prepareStage(interval.scale_factor);
  double dt_diffusion = std::numeric_limits<double>::infinity();
  double dt_hall = std::numeric_limits<double>::infinity();
  for (PatchRuntime& patch : m_patches) {
    m_nonideal->evaluateCoefficients(patch.binding.geometry, patch.workspace.primitives,
                                     interval.scale_factor, patch.nonideal, result.transport);
    dt_diffusion = std::min(dt_diffusion, m_nonideal->diffusionStableTimeStep(
                                              patch.binding.geometry, patch.nonideal));
    dt_hall = std::min(dt_hall,
                       m_nonideal->hallStableTimeStep(patch.binding.geometry, patch.nonideal));
  }
  m_forecast_dt_diffusion = m_nonideal->diffusionEnabled() ? dt_diffusion
                                                          : std::numeric_limits<double>::infinity();
  m_forecast_dt_hall = m_nonideal->hallEnabled() ? dt_hall : std::numeric_limits<double>::infinity();
  result.transport.min_diffusion_dt_code =
      std::isfinite(m_forecast_dt_diffusion) ? m_forecast_dt_diffusion : 0.0;
  result.transport.min_hall_dt_code = std::isfinite(m_forecast_dt_hall) ? m_forecast_dt_hall : 0.0;
  const double half = 0.5 * interval.dt_code;
  if (std::isfinite(m_forecast_dt_hall)) {
    const double subcycles = 2.0 * std::ceil(half / m_forecast_dt_hall);
    if (subcycles > static_cast<double>(m_nonideal->limits().hall_subcycle_cost_cap)) {
      MhdFailureContext failure =
          makeContext("mhd_hall", "forecast Hall subcycles exceed the per-interval cost cap");
      failure.predicted_stages = static_cast<std::uint64_t>(std::min(subcycles, 1.0e18));
      throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
    }
  }
  if (std::isfinite(m_forecast_dt_diffusion)) {
    const double ratio = half / m_forecast_dt_diffusion;
    if (ratio > static_cast<double>(m_nonideal->limits().rkl_switch_ratio)) {
      const MhdRkl2Plan plan = planRkl2(ratio, m_nonideal->limits().rkl_max_stages);
      if (2 * plan.totalStages() > m_nonideal->limits().rkl_stage_cost_cap) {
        MhdFailureContext failure =
            makeContext("mhd_rkl2", "forecast RKL2 stages exceed the per-interval cost cap");
        failure.predicted_stages = static_cast<std::uint64_t>(2 * plan.totalStages());
        throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(failure));
      }
    }
  }
}

void MhdLevelStepper::validateAccepted(double scale_factor, MhdLevelAdvanceResult& result) {
  double divergence = 0.0;
  for (PatchRuntime& patch : m_patches) {
    const MhdFaceFluxView view{&patch.binding.geometry,
                               patch.binding.faces.phi_x,
                               patch.binding.faces.phi_y,
                               patch.binding.faces.phi_z,
                               patch.binding.patch_id,
                               patch.binding.geometry_generation,
                               patch.binding.level};
    const std::size_t bad = m_operator.recoverPrimitives(patch.binding.geometry, *patch.binding.gas,
                                                         view, scale_factor, patch.workspace);
    if (bad != 0U) {
      MhdFailureContext failure =
          makeContext("mhd_level_stepper", "inadmissible cell(s) in the candidate end state");
      failure.patch_id = patch.binding.patch_id;
      throw MhdError(MhdErrorCode::kNonPositiveState, std::move(failure));
    }
    const double normalizer = magneticDivergenceNormalizer(
        patch.binding.geometry, view.phi_x, view.phi_y, view.phi_z);
    const double norm = maxNormalizedFaceDivergence(patch.binding.geometry, view.phi_x,
                                                    view.phi_y, view.phi_z, normalizer);
    divergence = std::max(divergence, norm);
    MhdFailureContext context = makeContext("mhd_level_stepper", "post-interval strict CT");
    context.patch_id = patch.binding.patch_id;
    context.geometry_generation = patch.binding.geometry_generation;
    context.scale_factor = scale_factor;
    requireDivergenceContract(norm, m_policy.divergence_tolerance, context);
  }
  result.divergence_max_normalized = divergence;
}

// ---------------------------------------------------------------------------
// Interval driver
// ---------------------------------------------------------------------------

void MhdLevelStepper::admitMemory(core::MemoryGovernor& governor,
                                  const MhdLevelAdvanceOptions& mode,
                                  const MhdResourceLimits& limits) {
  MhdMemoryRequest request{};
  request.geometry = m_patches.front().binding.geometry;
  request.reconstruction = m_policy.reconstruction;
  request.nonideal = mode.nonideal != nullptr;
  request.super_stepping = mode.nonideal != nullptr && mode.nonideal->diffusionEnabled();
  request.edge_integrals = mode.produce_edge_integrals;
  const double per_cell = scratchBytesPerActiveCell(request, m_patches.size());
  if (per_cell > static_cast<double>(limits.scratch_bytes_per_active_cell)) {
    MhdFailureContext failure = makeContext(
        "mhd_level_stepper",
        "MHD phase+scratch per active cell exceeds the configured cap (MP02); shrink the "
        "patch batch or tile the patch");
    failure.requested_bytes = static_cast<std::uint64_t>(per_cell);
    failure.safe_limit_bytes = limits.scratch_bytes_per_active_cell;
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded, std::move(failure));
  }
  m_admission.release();
  m_admission = MhdMemoryAdmission(governor, estimatePatchMemory(request), m_patches.size(),
                                   false, "mhd_level_stepper");
  m_admitted_request = request;
  m_admission_committed = false;
}

// Mirrors of shared faces are synchronized from their owner before any stage
// work (CHUI-MHD-REF v1.0 section 5.2). A difference beyond roundoff is a
// genuine inconsistency of the authoritative state and is refused.
void MhdLevelStepper::canonicalizeMirrors() {
  double scale = 0.0;
  for (PatchRuntime& patch : m_patches) {
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      for (const double value : faceLane(patch, axis)) {
        scale = std::max(scale, std::abs(value));
      }
    }
  }
  const double mismatch = reconcileFaces();
  if (mismatch > 1.0e-12 * std::max(scale, std::numeric_limits<double>::min())) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_level_stepper",
                               "shared face mirrors disagree beyond roundoff"));
  }
}

MhdLevelAdvanceResult MhdLevelStepper::advance(const MhdStageContext& interval,
                                               const MhdLevelAdvanceOptions& options) {
  if (m_admission.admitted()) {
    const bool needs_nonideal = options.nonideal != nullptr;
    const bool needs_super = needs_nonideal && options.nonideal->diffusionEnabled();
    if ((needs_nonideal && !m_admitted_request.nonideal) ||
        (needs_super && !m_admitted_request.super_stepping) ||
        (options.produce_edge_integrals && !m_admitted_request.edge_integrals)) {
      throw MhdError(MhdErrorCode::kMemoryBudgetExceeded,
                     makeContext("mhd_level_stepper",
                                 "interval mode needs workspaces that were not admitted"));
    }
  }
  if (!(interval.dt_code > 0.0) || !std::isfinite(interval.dt_code) ||
      !(interval.scale_factor > 0.0) || !(interval.scale_factor_end > 0.0) ||
      !std::isfinite(interval.hubble_rate_code) || !std::isfinite(interval.hubble_rate_end_code)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_level_stepper", "invalid interval (dt, a or H)"));
  }
  if (!interval.comoving_coordinates &&
      (interval.scale_factor != 1.0 || interval.scale_factor_end != 1.0 ||
       interval.hubble_rate_code != 0.0 || interval.hubble_rate_end_code != 0.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_level_stepper",
                               "isolated (non-comoving) runs require a = 1 and H = 0 exactly"));
  }
  m_nonideal = options.nonideal;
  m_track_integral = options.produce_edge_integrals;
  const bool constant_timeline = interval.scale_factor == interval.scale_factor_end &&
      interval.hubble_rate_code == 0.0 && interval.hubble_rate_end_code == 0.0;
  const bool can_retry_internally = constant_timeline || static_cast<bool>(options.timeline);

  canonicalizeMirrors();
  MhdLevelAdvanceResult result{};
  result.dt_requested_code = interval.dt_code;
  result.budget_begin = levelBudget(interval.scale_factor);
  snapshotRollback();
  MhdStageContext attempt = interval;
  for (std::int32_t retry = 0; retry <= m_policy.max_step_retries; ++retry) {
    m_interval_rkl_stages = 0;
    m_interval_hall_subcycles = 0;
    const std::uint64_t stages_before = result.ideal_stage_evaluations_attempted;
    try {
      if (m_nonideal != nullptr) {
        forecastNonIdeal(attempt, result);
        const MhdStageContext first_half = subIntervalContext(attempt, 0.0, 0.5);
        const MhdStageContext second_half = subIntervalContext(attempt, 0.5, 1.0);
        diffusionHalf(first_half, result);
        hallHalf(first_half, result);
        idealSsprk2(attempt, result);
        hallHalf(second_half, result);
        diffusionHalf(second_half, result);
      } else {
        idealSsprk2(attempt, result);
      }
      validateAccepted(attempt.scale_factor_end, result);
      result.accepted = true;
      result.dt_accepted_code = attempt.dt_code;
      result.retry_count = retry;
      result.ideal_stage_evaluations_accepted =
          result.ideal_stage_evaluations_attempted - stages_before;
      const double faces = static_cast<double>(std::max<std::uint64_t>(1U, result.ideal.faces_solved));
      result.fallback_fraction =
          static_cast<double>(result.ideal.hlle_fallback_faces + result.ideal.first_order_faces) /
          faces;
      result.fallback_fraction_gate_exceeded =
          result.fallback_fraction > m_policy.fallback_fraction_limit;
      result.budget_end = levelBudget(attempt.scale_factor_end);
      result.scratch_capacity_bytes = ownedCapacityBytes();
      if (m_admission.admitted() && !m_admission_committed) {
        m_admission.commit(result.scratch_capacity_bytes);
        m_admission_committed = true;
      }
      return result;
    } catch (const MhdError& error) {
      // Atomic rejection: every patch returns to the interval-start image.
      restoreRollback();
      result.last_rejection = error.what();
      if (!retryable(error.code())) {
        throw;
      }
      if (retry == m_policy.max_step_retries) {
        throw;
      }
      const double shorter = 0.5 * attempt.dt_code;
      if (!can_retry_internally) {
        // The scheduler owns the timeline: return the rejection.
        result.accepted = false;
        result.retry_count = retry + 1;
        result.dt_suggested_code = shorter;
        result.scratch_capacity_bytes = ownedCapacityBytes();
        return result;
      }
      attempt.dt_code = shorter;
      if (options.timeline) {
        const std::pair<double, double> end = options.timeline(attempt.time_code + shorter);
        attempt.scale_factor_end = end.first;
        attempt.hubble_rate_end_code = end.second;
      }
      result.retry_count = retry + 1;
    }
  }
  return result;
}

double MhdLevelStepper::estimateStableTimeStep(const MhdStageContext& interval,
                                               const MhdLevelAdvanceOptions& options) {
  double dt = std::numeric_limits<double>::infinity();
  for (PatchRuntime& patch : m_patches) {
    const MhdFaceFluxView view{&patch.binding.geometry,
                               patch.binding.faces.phi_x,
                               patch.binding.faces.phi_y,
                               patch.binding.faces.phi_z,
                               patch.binding.patch_id,
                               patch.binding.geometry_generation,
                               patch.binding.level};
    dt = std::min(dt, m_operator.estimateIdealTimeStep(patch.binding.geometry, *patch.binding.gas,
                                                       view, interval.scale_factor,
                                                       m_policy.cfl_ideal));
  }
  if (options.nonideal != nullptr) {
    m_nonideal = options.nonideal;
    MhdLevelAdvanceResult scratch{};
    canonicalizeMirrors();
    prepareStage(interval.scale_factor);
    double dt_diffusion = std::numeric_limits<double>::infinity();
    double dt_hall = std::numeric_limits<double>::infinity();
    for (PatchRuntime& patch : m_patches) {
      m_nonideal->evaluateCoefficients(patch.binding.geometry, patch.workspace.primitives,
                                       interval.scale_factor, patch.nonideal, scratch.transport);
      dt_diffusion = std::min(dt_diffusion, m_nonideal->diffusionStableTimeStep(
                                                patch.binding.geometry, patch.nonideal));
      dt_hall = std::min(dt_hall,
                         m_nonideal->hallStableTimeStep(patch.binding.geometry, patch.nonideal));
    }
    const MhdStiffnessLimits& limits = m_nonideal->limits();
    if (m_nonideal->hallEnabled() && std::isfinite(dt_hall)) {
      // 2 * ceil(dt/2 / dt_H) <= cap  <=>  dt <= cap * dt_H (cap even part).
      dt = std::min(dt, static_cast<double>(limits.hall_subcycle_cost_cap / 2) * 2.0 * dt_hall);
    }
    if (m_nonideal->diffusionEnabled() && std::isfinite(dt_diffusion)) {
      // Largest per-half ratio whose RKL2 plan fits half of the stage cap.
      const std::int32_t stages_per_half = std::max(2, limits.rkl_stage_cost_cap / 2);
      const std::int32_t cycles = std::max(1, stages_per_half / limits.rkl_max_stages);
      const std::int32_t per_cycle = std::min(limits.rkl_max_stages, stages_per_half);
      const double ratio = static_cast<double>(cycles) * rkl2StabilityFactor(per_cycle);
      dt = std::min(dt, 2.0 * ratio * dt_diffusion);
    }
  }
  return dt;
}

MhdConservationBudget MhdLevelStepper::levelBudget(double scale_factor) const {
  MhdConservationBudget total{};
  for (const PatchRuntime& patch : m_patches) {
    const MhdFaceFluxView view{&patch.binding.geometry,
                               patch.binding.faces.phi_x,
                               patch.binding.faces.phi_y,
                               patch.binding.faces.phi_z,
                               patch.binding.patch_id,
                               patch.binding.geometry_generation,
                               patch.binding.level};
    const MhdConservationBudget budget =
        m_operator.patchBudget(patch.binding.geometry, *patch.binding.gas, view, scale_factor);
    total.mass += budget.mass;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      total.momentum[axis] += budget.momentum[axis];
    }
    total.energy_thermal += budget.energy_thermal;
    total.energy_kinetic += budget.energy_kinetic;
    total.energy_magnetic += budget.energy_magnetic;
    total.energy_total += budget.energy_total;
    total.divergence_max_normalized =
        std::max(total.divergence_max_normalized, budget.divergence_max_normalized);
  }
  return total;
}

std::uint64_t MhdLevelStepper::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (const PatchRuntime& patch : m_patches) {
    bytes += patch.workspace.ownedCapacityBytes() + patch.nonideal.ownedCapacityBytes() +
        patch.rollback_gas.ownedCapacityBytes() + patch.stage0_gas.ownedCapacityBytes() +
        patch.integral.ownedCapacityBytes() + patch.stage0_integral.ownedCapacityBytes() +
        patch.yjm2_integral.ownedCapacityBytes();
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      bytes += patch.rollback_faces[axis].ownedCapacityBytes() +
          patch.stage0_faces[axis].ownedCapacityBytes();
    }
  }
  return bytes;
}

}  // namespace cosmosim::mhd
