// Topology construction and exchange planning for MhdLevelStepper.
//
// CHUI-MHD-REF v1.0 sections 5.2, 7.2-7.5 and CHUI-MHD-PERF v1.0 sections 6, 8:
// the composite leaf hierarchy (ratio-2 levels, arbitrary rectangular patch
// extents, sparse fine coverage, local and remote patches) is described by a
// global directory that every rank builds identically. All exchange plans are
// deterministic functions of that directory, so a sender and a receiver agree
// on item order and payload counts without exchanging metadata. Planning
// uses a sparse bin index (O(P log P) setup, no dense lattice lookup and no
// O(P^2) scans). UNVALIDATED source (pass 02); see
// docs/architecture/chui_mhd_implementation_pass02_status.md.

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <tuple>

#include "cosmosim/mhd/mhd_level_stepper.hpp"

namespace cosmosim::mhd {
namespace {

MhdFailureContext topologyContext(std::string detail) {
  MhdFailureContext context{};
  context.owner = "mhd_level_topology";
  context.detail = std::move(detail);
  return context;
}

[[noreturn]] void throwTopology(MhdErrorCode code, std::string detail) {
  throw MhdError(code, topologyContext(std::move(detail)));
}

[[nodiscard]] std::int64_t floorDiv(std::int64_t value, std::int64_t divisor) noexcept {
  std::int64_t quotient = value / divisor;
  if ((value % divisor != 0) && ((value < 0) != (divisor < 0))) {
    --quotient;
  }
  return quotient;
}

[[nodiscard]] std::int64_t wrapIndex(std::int64_t value, std::int64_t extent) noexcept {
  const std::int64_t r = value % extent;
  return r < 0 ? r + extent : r;
}

[[nodiscard]] bool physicalKind(MhdBoundaryKind kind) noexcept {
  return kind == MhdBoundaryKind::kPeriodic || kind == MhdBoundaryKind::kOpen ||
      kind == MhdBoundaryKind::kReflectingConductingWall;
}

[[nodiscard]] double spacingAlong(const MhdPatchGeometry& geometry, std::size_t axis) noexcept {
  return axis == 0U ? geometry.dx_comoving
                    : (axis == 1U ? geometry.dy_comoving : geometry.dz_comoving);
}

[[nodiscard]] double originAlong(const MhdPatchGeometry& geometry, std::size_t axis) noexcept {
  return axis == 0U ? geometry.origin_x_comoving
                    : (axis == 1U ? geometry.origin_y_comoving : geometry.origin_z_comoving);
}

}  // namespace

std::vector<std::vector<std::uint8_t>> MhdSerialHaloTransport::exchange(
    const std::vector<std::vector<std::uint8_t>>& send) {
  if (send.size() != 1U) {
    throwTopology(MhdErrorCode::kMpiEpochMismatch,
                  "serial transport received a multi-rank payload set");
  }
  return std::vector<std::vector<std::uint8_t>>(1U);
}

// ---------------------------------------------------------------------------
// Level lattice helpers
// ---------------------------------------------------------------------------

std::array<std::int64_t, 3> MhdLevelStepper::levelCells(std::int32_t level) const noexcept {
  const std::int32_t depth = level - m_domain.base_level;
  std::array<std::int64_t, 3> cells = m_domain.cells;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    cells[axis] <<= depth;
  }
  return cells;
}

bool MhdLevelStepper::periodicAxis(std::size_t axis) const noexcept {
  return m_domain.lower[axis] == MhdBoundaryKind::kPeriodic;
}

MhdPatchBoundary MhdLevelStepper::boundaryForPatch(const TopoPatch& patch) const {
  const std::array<std::int64_t, 3> extent = levelCells(patch.level);
  MhdPatchBoundary boundary{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (periodicAxis(axis)) {
      const bool spans = patch.box.lo[axis] == 0 && patch.box.hi[axis] == extent[axis];
      boundary.lower[axis] = spans ? MhdBoundaryKind::kPeriodic : MhdBoundaryKind::kNeighborPatch;
      boundary.upper[axis] = boundary.lower[axis];
      continue;
    }
    boundary.lower[axis] =
        patch.box.lo[axis] == 0 ? m_domain.lower[axis] : MhdBoundaryKind::kNeighborPatch;
    boundary.upper[axis] =
        patch.box.hi[axis] == extent[axis] ? m_domain.upper[axis] : MhdBoundaryKind::kNeighborPatch;
  }
  return boundary;
}

// ---------------------------------------------------------------------------
// Directory construction and validation
// ---------------------------------------------------------------------------

void MhdLevelStepper::buildTopology(std::vector<MhdRemotePatchDescriptor> remote) {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (!physicalKind(m_domain.lower[axis]) || !physicalKind(m_domain.upper[axis])) {
      throwTopology(MhdErrorCode::kUnsupportedPhysics,
                    "domain boundaries must be periodic, open or conducting wall; shearing "
                    "boundaries are not part of REF-v1");
    }
    if ((m_domain.lower[axis] == MhdBoundaryKind::kPeriodic) !=
        (m_domain.upper[axis] == MhdBoundaryKind::kPeriodic)) {
      throwTopology(MhdErrorCode::kUnsupportedPhysics, "one-sided periodic domain");
    }
    if (m_domain.cells[axis] <= 0) {
      throwTopology(MhdErrorCode::kInvalidField, "base-level domain needs positive cell counts");
    }
  }
  if (m_domain.base_level < 0) {
    throwTopology(MhdErrorCode::kInvalidField, "negative base level");
  }
  const int my_rank = m_transport->rank();
  const int world = m_transport->size();
  const std::size_t ghost = m_policy.ghostCells();

  // Directory: local bindings first (validated), then remote descriptors.
  std::vector<TopoPatch> directory;
  directory.reserve(m_local.size() + remote.size());
  // Per-level spacing consistency over local patches (ratio-2 levels).
  std::map<std::int32_t, std::array<double, 3>> level_spacing;
  for (std::size_t index = 0; index < m_local.size(); ++index) {
    const MhdLevelPatchBinding& binding = m_local[index].binding;
    const MhdPatchGeometry& geometry = binding.geometry;
    if (!geometry.valid()) {
      throwTopology(MhdErrorCode::kInvalidField, "invalid local patch geometry");
    }
    if (binding.gas == nullptr || binding.gas->size() != geometry.cellCount() ||
        !binding.faces.valid() || binding.faces.patch_id != binding.patch_id ||
        binding.faces.geometry_generation != binding.geometry_generation ||
        binding.faces.phi_x.size() != geometry.faceCountX() ||
        binding.faces.phi_y.size() != geometry.faceCountY() ||
        binding.faces.phi_z.size() != geometry.faceCountZ()) {
      throwTopology(MhdErrorCode::kStaleGeometry,
                    "patch gas/face views are missing, mis-sized or carry a stale "
                    "(patch_id, geometry_generation)");
    }
    if (binding.owner_rank != my_rank) {
      throwTopology(MhdErrorCode::kMpiEpochMismatch,
                    "local patch binding names an owner rank other than this rank");
    }
    if (binding.level < m_domain.base_level || binding.level > m_domain.base_level + 30) {
      throwTopology(MhdErrorCode::kInvalidField, "patch level outside the supported range");
    }
    if (std::min({geometry.nx, geometry.ny, geometry.nz}) < ghost) {
      throwTopology(MhdErrorCode::kInvalidField,
                    "patch extent is smaller than the reconstruction ghost depth");
    }
    const std::array<double, 3> spacing{geometry.dx_comoving, geometry.dy_comoving,
                                        geometry.dz_comoving};
    const auto found = level_spacing.find(binding.level);
    if (found == level_spacing.end()) {
      level_spacing.emplace(binding.level, spacing);
    } else {
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        if (std::abs(found->second[axis] - spacing[axis]) > 1.0e-12 * spacing[axis]) {
          throwTopology(MhdErrorCode::kStaleGeometry,
                        "patches of one level must share one cell spacing");
        }
      }
    }
    TopoPatch entry{};
    entry.patch_id = binding.patch_id;
    entry.geometry_generation = binding.geometry_generation;
    entry.level = binding.level;
    entry.owner_rank = binding.owner_rank;
    entry.box.lo = binding.origin_cell;
    entry.box.hi = {binding.origin_cell[0] + static_cast<std::int64_t>(geometry.nx),
                    binding.origin_cell[1] + static_cast<std::int64_t>(geometry.ny),
                    binding.origin_cell[2] + static_cast<std::int64_t>(geometry.nz)};
    entry.local = index;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const double expected = m_domain.origin_comoving[axis] +
          static_cast<double>(binding.origin_cell[axis]) * spacing[axis];
      if (std::abs(originAlong(geometry, axis) - expected) >
          1.0e-9 * std::max(std::abs(expected), spacing[axis])) {
        throwTopology(MhdErrorCode::kStaleGeometry,
                      "patch comoving origin disagrees with its lattice origin");
      }
    }
    directory.push_back(entry);
  }
  // Ratio-2 spacing across the locally present levels.
  for (auto it = level_spacing.begin(); it != level_spacing.end(); ++it) {
    const auto next = std::next(it);
    if (next == level_spacing.end()) {
      break;
    }
    const double factor = std::ldexp(1.0, -(next->first - it->first));
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (std::abs(next->second[axis] - factor * it->second[axis]) >
          1.0e-12 * it->second[axis]) {
        throwTopology(MhdErrorCode::kStaleGeometry,
                      "level spacings are not related by the refinement ratio 2");
      }
    }
  }
  for (const MhdRemotePatchDescriptor& descriptor : remote) {
    if (descriptor.owner_rank == my_rank || descriptor.owner_rank < 0 ||
        descriptor.owner_rank >= world) {
      throwTopology(MhdErrorCode::kMpiEpochMismatch,
                    "remote patch descriptor names this rank or an invalid rank");
    }
    if (descriptor.level < m_domain.base_level || descriptor.level > m_domain.base_level + 30 ||
        descriptor.cells[0] < static_cast<std::int64_t>(ghost) ||
        descriptor.cells[1] < static_cast<std::int64_t>(ghost) ||
        descriptor.cells[2] < static_cast<std::int64_t>(ghost)) {
      throwTopology(MhdErrorCode::kInvalidField, "invalid remote patch descriptor");
    }
    TopoPatch entry{};
    entry.patch_id = descriptor.patch_id;
    entry.geometry_generation = descriptor.geometry_generation;
    entry.level = descriptor.level;
    entry.owner_rank = descriptor.owner_rank;
    entry.box.lo = descriptor.origin_cell;
    entry.box.hi = {descriptor.origin_cell[0] + descriptor.cells[0],
                    descriptor.origin_cell[1] + descriptor.cells[1],
                    descriptor.origin_cell[2] + descriptor.cells[2]};
    entry.local = k_npos;
    directory.push_back(entry);
  }
  if (directory.empty() && world == 1) {
    throwTopology(MhdErrorCode::kInvalidField, "a hierarchy needs at least one patch");
  }
  std::sort(directory.begin(), directory.end(), [](const TopoPatch& lhs, const TopoPatch& rhs) {
    return std::tie(lhs.level, lhs.patch_id) < std::tie(rhs.level, rhs.patch_id);
  });
  {
    std::vector<std::uint64_t> ids;
    ids.reserve(directory.size());
    for (const TopoPatch& entry : directory) {
      ids.push_back(entry.patch_id);
    }
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) {
      throwTopology(MhdErrorCode::kStaleGeometry, "duplicate patch_id in the hierarchy");
    }
  }
  m_max_level = m_domain.base_level;
  for (const TopoPatch& entry : directory) {
    m_max_level = std::max(m_max_level, entry.level);
    const std::array<std::int64_t, 3> extent = levelCells(entry.level);
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (entry.box.lo[axis] < 0 || entry.box.hi[axis] > extent[axis] ||
          entry.box.lo[axis] >= entry.box.hi[axis]) {
        throwTopology(MhdErrorCode::kStaleGeometry, "patch box outside its level domain");
      }
      if (entry.level > m_domain.base_level &&
          ((entry.box.lo[axis] & 1) != 0 || (entry.box.hi[axis] & 1) != 0)) {
        throwTopology(MhdErrorCode::kStaleGeometry,
                      "refined patch boxes must be aligned to the ratio-2 parent lattice");
      }
    }
  }
  m_level_count = m_max_level - m_domain.base_level + 1;
  m_topo = std::move(directory);
  for (std::size_t index = 0; index < m_topo.size(); ++index) {
    if (m_topo[index].local != k_npos) {
      m_local[m_topo[index].local].topo = static_cast<std::uint32_t>(index);
    }
  }
  buildSpatialIndex();

  // Disjointness (same level and adjacent levels) and exact domain coverage.
  std::vector<std::pair<std::uint32_t, std::array<std::int64_t, 3>>> hits;
  double covered = 0.0;
  for (std::uint32_t index = 0; index < m_topo.size(); ++index) {
    const TopoPatch& patch = m_topo[index];
    queryLevel(patch.level, patch.box, hits);
    for (const auto& [other, shift] : hits) {
      if (other == index && shift == std::array<std::int64_t, 3>{0, 0, 0}) {
        continue;
      }
      throwTopology(MhdErrorCode::kStaleGeometry, "overlapping leaf patches on one level");
    }
    if (patch.level > m_domain.base_level) {
      Box coarse{};
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        coarse.lo[axis] = patch.box.lo[axis] / 2;
        coarse.hi[axis] = patch.box.hi[axis] / 2;
      }
      queryLevel(patch.level - 1, coarse, hits);
      if (!hits.empty()) {
        throwTopology(MhdErrorCode::kStaleGeometry,
                      "a refined leaf overlaps a coarser leaf (not a leaf tiling)");
      }
    }
    covered += std::ldexp(static_cast<double>(patch.box.volume()),
                          -3 * (patch.level - m_domain.base_level));
  }
  const double domain_volume = static_cast<double>(m_domain.cells[0]) *
      static_cast<double>(m_domain.cells[1]) * static_cast<double>(m_domain.cells[2]);
  if (!m_topo.empty() && covered != domain_volume) {
    std::ostringstream stream;
    stream << "leaf tiling does not cover the domain exactly once (covered base-level volume "
           << covered << " versus " << domain_volume << "); the directory must list every leaf "
           << "of every rank";
    throwTopology(MhdErrorCode::kStaleGeometry, stream.str());
  }

  // Patch-side kinds of local patches.
  for (PatchRuntime& runtime : m_local) {
    runtime.boundary = boundaryForPatch(m_topo[runtime.topo]);
    validateBoundarySupport(runtime.boundary, topologyContext(""));
  }

  // Halo plans for local destinations; items whose source is local and whose
  // destination is remote are collected for exchange packing.
  m_halo_items.clear();
  std::vector<std::uint8_t> plan_needed(m_topo.size(), 0U);
  for (const PatchRuntime& runtime : m_local) {
    const TopoPatch& patch = m_topo[runtime.topo];
    // Any patch within ghost reach of a local patch, at levels L-1..L+1, may
    // need data from it.
    for (std::int32_t level = std::max(patch.level - 1, m_domain.base_level);
         level <= std::min(patch.level + 1, m_max_level); ++level) {
      Box query = patch.box;
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        std::int64_t lo = patch.box.lo[axis];
        std::int64_t hi = patch.box.hi[axis];
        if (level > patch.level) {
          lo *= 2;
          hi *= 2;
        } else if (level < patch.level) {
          lo = floorDiv(lo, 2);
          hi = floorDiv(hi + 1, 2);
        }
        const std::int64_t reach = static_cast<std::int64_t>(ghost) + 4;
        query.lo[axis] = lo - reach;
        query.hi[axis] = hi + reach;
      }
      queryLevel(level, query, hits);
      for (const auto& [other, shift] : hits) {
        (void)shift;
        plan_needed[other] = 1U;
      }
    }
  }
  for (PatchRuntime& runtime : m_local) {
    planHalo(runtime.topo, runtime.halo);
  }
  for (std::uint32_t index = 0; index < m_topo.size(); ++index) {
    if (m_topo[index].local != k_npos) {
      const HaloPlan& plan = m_local[m_topo[index].local].halo;
      for (const Item& item : plan.items) {
        if (m_topo[item.src].owner_rank != m_topo[item.dst].owner_rank) {
          m_halo_items.push_back(item);
        }
      }
      continue;
    }
    if (plan_needed[index] == 0U) {
      continue;
    }
    HaloPlan remote_plan{};
    planHalo(index, remote_plan);
    for (const Item& item : remote_plan.items) {
      if (m_topo[item.src].local != k_npos) {
        m_halo_items.push_back(item);
      }
    }
  }
  planReconciliation();
  planCrossLevel();
  finalizeItemLists();
}

// ---------------------------------------------------------------------------
// Sparse bin index
// ---------------------------------------------------------------------------

void MhdLevelStepper::buildSpatialIndex() {
  const std::size_t levels = static_cast<std::size_t>(m_max_level - m_domain.base_level + 1);
  m_bins.assign(levels, {});
  m_bin_size.assign(levels, {1, 1, 1});
  for (const TopoPatch& patch : m_topo) {
    const auto slot = static_cast<std::size_t>(patch.level - m_domain.base_level);
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      m_bin_size[slot][axis] =
          std::max(m_bin_size[slot][axis], patch.box.hi[axis] - patch.box.lo[axis]);
    }
  }
  for (std::uint32_t index = 0; index < m_topo.size(); ++index) {
    const TopoPatch& patch = m_topo[index];
    const auto slot = static_cast<std::size_t>(patch.level - m_domain.base_level);
    const std::array<std::int64_t, 3>& size = m_bin_size[slot];
    std::array<std::int64_t, 3> lo{};
    std::array<std::int64_t, 3> hi{};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      lo[axis] = floorDiv(patch.box.lo[axis], size[axis]);
      hi[axis] = floorDiv(patch.box.hi[axis] - 1, size[axis]);
    }
    for (std::int64_t bz = lo[2]; bz <= hi[2]; ++bz) {
      for (std::int64_t by = lo[1]; by <= hi[1]; ++by) {
        for (std::int64_t bx = lo[0]; bx <= hi[0]; ++bx) {
          m_bins[slot].push_back({{bx, by, bz}, index});
        }
      }
    }
  }
  for (auto& bins : m_bins) {
    std::sort(bins.begin(), bins.end());
  }
}

void MhdLevelStepper::queryLevel(
    std::int32_t level, const Box& query,
    std::vector<std::pair<std::uint32_t, std::array<std::int64_t, 3>>>& out) const {
  out.clear();
  if (level < m_domain.base_level || level > m_max_level || query.empty()) {
    return;
  }
  const auto slot = static_cast<std::size_t>(level - m_domain.base_level);
  const std::array<std::int64_t, 3> extent = levelCells(level);
  const std::array<std::int64_t, 3>& size = m_bin_size[slot];
  const auto& bins = m_bins[slot];
  // Periodic images that the query can reach.
  std::array<std::array<std::int64_t, 3>, 3> shift_options{};
  std::array<int, 3> shift_count{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    int count = 0;
    shift_options[axis][count++] = 0;
    if (periodicAxis(axis)) {
      if (query.lo[axis] < 0) {
        shift_options[axis][count++] = -extent[axis];
      }
      if (query.hi[axis] > extent[axis]) {
        shift_options[axis][count++] = extent[axis];
      }
    }
    shift_count[axis] = count;
  }
  for (int sx = 0; sx < shift_count[0]; ++sx) {
    for (int sy = 0; sy < shift_count[1]; ++sy) {
      for (int sz = 0; sz < shift_count[2]; ++sz) {
        const std::array<std::int64_t, 3> shift{shift_options[0][sx], shift_options[1][sy],
                                                shift_options[2][sz]};
        // Patch image = patch.box + shift; query in domain coordinates.
        Box local = query;
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          local.lo[axis] = std::max<std::int64_t>(query.lo[axis] - shift[axis], 0);
          local.hi[axis] = std::min<std::int64_t>(query.hi[axis] - shift[axis], extent[axis]);
        }
        if (local.empty()) {
          continue;
        }
        std::array<std::int64_t, 3> blo{};
        std::array<std::int64_t, 3> bhi{};
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          blo[axis] = floorDiv(local.lo[axis], size[axis]);
          bhi[axis] = floorDiv(local.hi[axis] - 1, size[axis]);
        }
        for (std::int64_t bz = blo[2]; bz <= bhi[2]; ++bz) {
          for (std::int64_t by = blo[1]; by <= bhi[1]; ++by) {
            for (std::int64_t bx = blo[0]; bx <= bhi[0]; ++bx) {
              const std::array<std::int64_t, 3> key{bx, by, bz};
              auto range = std::equal_range(
                  bins.begin(), bins.end(),
                  std::pair<std::array<std::int64_t, 3>, std::uint32_t>{key, 0U},
                  [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
              for (auto it = range.first; it != range.second; ++it) {
                const Box& box = m_topo[it->second].box;
                bool intersects = true;
                for (std::size_t axis = 0; axis < 3U; ++axis) {
                  intersects = intersects && box.lo[axis] < local.hi[axis] &&
                      local.lo[axis] < box.hi[axis];
                }
                if (intersects) {
                  out.emplace_back(it->second, shift);
                }
              }
            }
          }
        }
      }
    }
  }
  std::sort(out.begin(), out.end(), [&](const auto& lhs, const auto& rhs) {
    return std::tie(m_topo[lhs.first].patch_id, lhs.second) <
        std::tie(m_topo[rhs.first].patch_id, rhs.second);
  });
  out.erase(std::unique(out.begin(), out.end()), out.end());
}

std::int64_t MhdLevelStepper::findCellHolder(std::int32_t level,
                                             const std::array<std::int64_t, 3>& cell,
                                             std::array<std::int64_t, 3>& wrapped) const {
  if (level < m_domain.base_level || level > m_max_level) {
    return -1;
  }
  const std::array<std::int64_t, 3> extent = levelCells(level);
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (periodicAxis(axis)) {
      wrapped[axis] = wrapIndex(cell[axis], extent[axis]);
    } else {
      if (cell[axis] < 0 || cell[axis] >= extent[axis]) {
        return -1;
      }
      wrapped[axis] = cell[axis];
    }
  }
  const auto slot = static_cast<std::size_t>(level - m_domain.base_level);
  const std::array<std::int64_t, 3>& size = m_bin_size[slot];
  const std::array<std::int64_t, 3> key{floorDiv(wrapped[0], size[0]),
                                        floorDiv(wrapped[1], size[1]),
                                        floorDiv(wrapped[2], size[2])};
  const auto& bins = m_bins[slot];
  auto range = std::equal_range(
      bins.begin(), bins.end(), std::pair<std::array<std::int64_t, 3>, std::uint32_t>{key, 0U},
      [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
  for (auto it = range.first; it != range.second; ++it) {
    const Box& box = m_topo[it->second].box;
    if (wrapped[0] >= box.lo[0] && wrapped[0] < box.hi[0] && wrapped[1] >= box.lo[1] &&
        wrapped[1] < box.hi[1] && wrapped[2] >= box.lo[2] && wrapped[2] < box.hi[2]) {
      return static_cast<std::int64_t>(it->second);
    }
  }
  return -1;
}

// ---------------------------------------------------------------------------
// Halo planning (ghost cells, ghost staggered faces, prolongation samples)
// ---------------------------------------------------------------------------

std::size_t MhdLevelStepper::itemValueCount(const Item& item) const noexcept {
  const auto volume = static_cast<std::size_t>(item.dst_box.volume());
  switch (item.kind) {
    case ItemKind::kGhostCellDirect:
    case ItemKind::kGhostCellRestricted:
      return 9U * volume;
    case ItemKind::kSampleCellDirect:
    case ItemKind::kSampleCellRestricted:
      return 9U;
    case ItemKind::kSampleFaceDirect:
      return 1U;
    case ItemKind::kGasRecordOwner:
    case ItemKind::kGasRecordFine:
      return MhdBoundaryFaceFluxes::k_values_per_face * volume;
    default:
      return volume;
  }
}

void MhdLevelStepper::planHalo(std::uint32_t topo_index, HaloPlan& plan) const {
  plan = HaloPlan{};
  const TopoPatch& patch = m_topo[topo_index];
  const MhdPatchBoundary boundary = boundaryForPatch(patch);
  const auto g = static_cast<std::int64_t>(m_policy.ghostCells());
  const auto gt = static_cast<std::int64_t>(k_mhd_staggered_tangential_ghost);
  const std::array<std::int64_t, 3> n{patch.box.hi[0] - patch.box.lo[0],
                                      patch.box.hi[1] - patch.box.lo[1],
                                      patch.box.hi[2] - patch.box.lo[2]};
  const std::int32_t level = patch.level;
  std::vector<std::pair<std::uint32_t, std::array<std::int64_t, 3>>> hits;
  std::map<std::tuple<std::uint32_t, std::uint8_t, std::uint8_t, std::int64_t, std::int64_t,
                      std::int64_t>,
           std::int32_t>
      sample_slots;
  const auto side_external = [&](std::size_t axis, int pattern) {
    if (pattern == 0) {
      return true;
    }
    const MhdBoundaryKind kind = pattern < 0 ? boundary.lower[axis] : boundary.upper[axis];
    return kind == MhdBoundaryKind::kNeighborPatch;
  };
  const auto add_item = [&](Item item) {
    item.dst = topo_index;
    plan.items.push_back(item);
  };
  const auto sample = [&](ItemKind kind, std::uint8_t axis, std::uint32_t src,
                          const std::array<std::int64_t, 3>& src_local) {
    const auto key = std::make_tuple(src, static_cast<std::uint8_t>(kind), axis, src_local[0],
                                     src_local[1], src_local[2]);
    const auto found = sample_slots.find(key);
    if (found != sample_slots.end()) {
      return found->second;
    }
    Item item{};
    item.kind = kind;
    item.axis = axis;
    item.src = src;
    item.dst_box = Box{{0, 0, 0}, {1, 1, 1}};
    item.src_offset = src_local;
    item.slot = kind == ItemKind::kSampleFaceDirect ? plan.face_sample_count++
                                                    : plan.cell_sample_count++;
    add_item(item);
    const auto slot = static_cast<std::int32_t>(item.slot);
    sample_slots.emplace(key, slot);
    return slot;
  };
  const std::array<std::int64_t, 3> extent = levelCells(level);
  const std::array<std::int64_t, 3> coarse_extent =
      level > m_domain.base_level ? levelCells(level - 1) : extent;

  // --- Ghost cells -------------------------------------------------------
  for (int pz = -1; pz <= 1; ++pz) {
    for (int py = -1; py <= 1; ++py) {
      for (int px = -1; px <= 1; ++px) {
        const std::array<int, 3> pattern{px, py, pz};
        if (px == 0 && py == 0 && pz == 0) {
          continue;
        }
        bool external = true;
        Box sub{};
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          external = external && side_external(axis, pattern[axis]);
          sub.lo[axis] = pattern[axis] < 0 ? -g : (pattern[axis] == 0 ? 0 : n[axis]);
          sub.hi[axis] = pattern[axis] < 0 ? 0 : (pattern[axis] == 0 ? n[axis] : n[axis] + g);
        }
        if (!external) {
          continue;
        }
        Box global = sub;
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          global.lo[axis] += patch.box.lo[axis];
          global.hi[axis] += patch.box.lo[axis];
        }
        const std::array<std::int64_t, 3> sub_n{sub.hi[0] - sub.lo[0], sub.hi[1] - sub.lo[1],
                                                sub.hi[2] - sub.lo[2]};
        std::vector<std::uint8_t> marked(static_cast<std::size_t>(sub.volume()), 0U);
        const auto mark = [&](const Box& global_box) {
          for (std::int64_t k = global_box.lo[2]; k < global_box.hi[2]; ++k) {
            for (std::int64_t j = global_box.lo[1]; j < global_box.hi[1]; ++j) {
              for (std::int64_t i = global_box.lo[0]; i < global_box.hi[0]; ++i) {
                const std::int64_t li = i - global.lo[0];
                const std::int64_t lj = j - global.lo[1];
                const std::int64_t lk = k - global.lo[2];
                marked[static_cast<std::size_t>(li + sub_n[0] * (lj + sub_n[1] * lk))] = 1U;
              }
            }
          }
        };
        // (1) Same-level leaves.
        queryLevel(level, global, hits);
        for (const auto& [src, shift] : hits) {
          const TopoPatch& source = m_topo[src];
          Box overlap{};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            overlap.lo[axis] = std::max(global.lo[axis], source.box.lo[axis] + shift[axis]);
            overlap.hi[axis] = std::min(global.hi[axis], source.box.hi[axis] + shift[axis]);
          }
          if (overlap.empty()) {
            continue;
          }
          Item item{};
          item.kind = ItemKind::kGhostCellDirect;
          item.src = src;
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            item.dst_box.lo[axis] = overlap.lo[axis] - patch.box.lo[axis];
            item.dst_box.hi[axis] = overlap.hi[axis] - patch.box.lo[axis];
            item.src_offset[axis] = patch.box.lo[axis] - source.box.lo[axis] - shift[axis];
          }
          add_item(item);
          mark(overlap);
        }
        // (2) Finer leaves: restriction of eight children.
        if (level < m_max_level) {
          Box fine_query{};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            fine_query.lo[axis] = 2 * global.lo[axis];
            fine_query.hi[axis] = 2 * global.hi[axis];
          }
          queryLevel(level + 1, fine_query, hits);
          for (const auto& [src, shift] : hits) {
            const TopoPatch& source = m_topo[src];
            Box overlap{};
            for (std::size_t axis = 0; axis < 3U; ++axis) {
              overlap.lo[axis] =
                  std::max(global.lo[axis], (source.box.lo[axis] + shift[axis]) / 2);
              overlap.hi[axis] =
                  std::min(global.hi[axis], (source.box.hi[axis] + shift[axis]) / 2);
            }
            if (overlap.empty()) {
              continue;
            }
            Item item{};
            item.kind = ItemKind::kGhostCellRestricted;
            item.scale = 2U;
            item.src = src;
            for (std::size_t axis = 0; axis < 3U; ++axis) {
              item.dst_box.lo[axis] = overlap.lo[axis] - patch.box.lo[axis];
              item.dst_box.hi[axis] = overlap.hi[axis] - patch.box.lo[axis];
              item.src_offset[axis] =
                  2 * patch.box.lo[axis] - source.box.lo[axis] - shift[axis];
            }
            add_item(item);
            mark(overlap);
          }
        }
        // (3) Coarser leaves: prolongation from coarse samples.
        for (std::int64_t lk = 0; lk < sub_n[2]; ++lk) {
          for (std::int64_t lj = 0; lj < sub_n[1]; ++lj) {
            for (std::int64_t li = 0; li < sub_n[0]; ++li) {
              if (marked[static_cast<std::size_t>(li + sub_n[0] * (lj + sub_n[1] * lk))] != 0U) {
                continue;
              }
              const std::array<std::int64_t, 3> cell{global.lo[0] + li, global.lo[1] + lj,
                                                     global.lo[2] + lk};
              if (level == m_domain.base_level) {
                throwTopology(MhdErrorCode::kStaleGeometry,
                              "base-level ghost cell is not covered by any leaf (tiling gap)");
              }
              std::array<std::int64_t, 3> wrapped{};
              for (std::size_t axis = 0; axis < 3U; ++axis) {
                wrapped[axis] = periodicAxis(axis) ? wrapIndex(cell[axis], extent[axis])
                                                   : cell[axis];
              }
              const std::array<std::int64_t, 3> coarse{wrapped[0] / 2, wrapped[1] / 2,
                                                       wrapped[2] / 2};
              std::array<std::int64_t, 3> coarse_wrapped{};
              const std::int64_t holder = findCellHolder(level - 1, coarse, coarse_wrapped);
              if (holder < 0) {
                throwTopology(MhdErrorCode::kStaleGeometry,
                              "ghost cell is covered neither by a same-level, a finer nor the "
                              "parent level leaf: MHD requires a 2:1 balanced leaf tiling");
              }
              ProlongCell recipe{};
              recipe.dst = {static_cast<std::ptrdiff_t>(cell[0] - patch.box.lo[0]),
                            static_cast<std::ptrdiff_t>(cell[1] - patch.box.lo[1]),
                            static_cast<std::ptrdiff_t>(cell[2] - patch.box.lo[2])};
              const auto holder_index = static_cast<std::uint32_t>(holder);
              const Box& holder_box = m_topo[holder_index].box;
              recipe.center = sample(
                  ItemKind::kSampleCellDirect, 0U, holder_index,
                  {coarse_wrapped[0] - holder_box.lo[0], coarse_wrapped[1] - holder_box.lo[1],
                   coarse_wrapped[2] - holder_box.lo[2]});
              for (std::size_t axis = 0; axis < 3U; ++axis) {
                recipe.half[axis] = (wrapped[axis] & 1) == 0 ? static_cast<std::int8_t>(-1)
                                                             : static_cast<std::int8_t>(1);
                for (int side = 0; side < 2; ++side) {
                  std::array<std::int64_t, 3> neighbour = coarse_wrapped;
                  neighbour[axis] += side == 0 ? -1 : 1;
                  if (!periodicAxis(axis) &&
                      (neighbour[axis] < 0 || neighbour[axis] >= coarse_extent[axis])) {
                    continue;  // zero slope toward a physical boundary (counted at use)
                  }
                  std::array<std::int64_t, 3> neighbour_wrapped{};
                  const std::int64_t coarse_holder =
                      findCellHolder(level - 1, neighbour, neighbour_wrapped);
                  std::int32_t slot = -1;
                  if (coarse_holder >= 0) {
                    const Box& box = m_topo[static_cast<std::size_t>(coarse_holder)].box;
                    slot = sample(ItemKind::kSampleCellDirect, 0U,
                                  static_cast<std::uint32_t>(coarse_holder),
                                  {neighbour_wrapped[0] - box.lo[0], neighbour_wrapped[1] - box.lo[1],
                                   neighbour_wrapped[2] - box.lo[2]});
                  } else {
                    const std::array<std::int64_t, 3> child{2 * neighbour_wrapped[0],
                                                            2 * neighbour_wrapped[1],
                                                            2 * neighbour_wrapped[2]};
                    std::array<std::int64_t, 3> child_wrapped{};
                    const std::int64_t fine_holder = findCellHolder(level, child, child_wrapped);
                    if (fine_holder >= 0) {
                      const Box& box = m_topo[static_cast<std::size_t>(fine_holder)].box;
                      slot = sample(ItemKind::kSampleCellRestricted, 0U,
                                    static_cast<std::uint32_t>(fine_holder),
                                    {child_wrapped[0] - box.lo[0], child_wrapped[1] - box.lo[1],
                                     child_wrapped[2] - box.lo[2]});
                    }
                  }
                  recipe.neighbour[2U * axis + static_cast<std::size_t>(side)] = slot;
                }
              }
              plan.prolong_cells.push_back(recipe);
            }
          }
        }
      }
    }
  }

  // --- Ghost staggered faces (tangential ghost rows) ---------------------
  for (std::size_t normal = 0; normal < 3U; ++normal) {
    const std::size_t t1 = (normal + 1U) % 3U;
    const std::size_t t2 = (normal + 2U) % 3U;
    for (int p2 = -1; p2 <= 1; ++p2) {
      for (int p1 = -1; p1 <= 1; ++p1) {
        if (p1 == 0 && p2 == 0) {
          continue;
        }
        if (!side_external(t1, p1) || !side_external(t2, p2)) {
          continue;
        }
        Box sub{};
        sub.lo[normal] = 0;
        sub.hi[normal] = n[normal] + 1;
        const std::array<std::pair<std::size_t, int>, 2> tangential{{{t1, p1}, {t2, p2}}};
        for (const auto& [axis, pattern] : tangential) {
          sub.lo[axis] = pattern < 0 ? -gt : (pattern == 0 ? 0 : n[axis]);
          sub.hi[axis] = pattern < 0 ? 0 : (pattern == 0 ? n[axis] : n[axis] + gt);
        }
        Box global = sub;
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          global.lo[axis] += patch.box.lo[axis];
          global.hi[axis] += patch.box.lo[axis];
        }
        const std::array<std::int64_t, 3> sub_n{sub.hi[0] - sub.lo[0], sub.hi[1] - sub.lo[1],
                                                sub.hi[2] - sub.lo[2]};
        std::vector<std::uint8_t> marked(static_cast<std::size_t>(sub.volume()), 0U);
        const auto mark = [&](const Box& global_box) {
          for (std::int64_t k = global_box.lo[2]; k < global_box.hi[2]; ++k) {
            for (std::int64_t j = global_box.lo[1]; j < global_box.hi[1]; ++j) {
              for (std::int64_t i = global_box.lo[0]; i < global_box.hi[0]; ++i) {
                const std::int64_t li = i - global.lo[0];
                const std::int64_t lj = j - global.lo[1];
                const std::int64_t lk = k - global.lo[2];
                marked[static_cast<std::size_t>(li + sub_n[0] * (lj + sub_n[1] * lk))] = 1U;
              }
            }
          }
        };
        // Cells adjacent to the faces (for leaf queries).
        Box cell_query = global;
        cell_query.lo[normal] -= 1;
        // (1) Same-level holders (closed normal range).
        queryLevel(level, cell_query, hits);
        for (const auto& [src, shift] : hits) {
          const TopoPatch& source = m_topo[src];
          Box overlap{};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            const std::int64_t lo = source.box.lo[axis] + shift[axis];
            const std::int64_t hi = source.box.hi[axis] + shift[axis] + (axis == normal ? 1 : 0);
            overlap.lo[axis] = std::max(global.lo[axis], lo);
            overlap.hi[axis] = std::min(global.hi[axis], hi);
          }
          if (overlap.empty()) {
            continue;
          }
          Item item{};
          item.kind = ItemKind::kGhostFaceDirect;
          item.axis = static_cast<std::uint8_t>(normal);
          item.src = src;
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            item.dst_box.lo[axis] = overlap.lo[axis] - patch.box.lo[axis];
            item.dst_box.hi[axis] = overlap.hi[axis] - patch.box.lo[axis];
            item.src_offset[axis] = patch.box.lo[axis] - source.box.lo[axis] - shift[axis];
          }
          add_item(item);
          mark(overlap);
        }
        // (2) Finer holders: sum of the four fine faces of a coarse face.
        if (level < m_max_level) {
          Box fine_query{};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            fine_query.lo[axis] = 2 * cell_query.lo[axis];
            fine_query.hi[axis] = 2 * cell_query.hi[axis];
          }
          queryLevel(level + 1, fine_query, hits);
          for (const auto& [src, shift] : hits) {
            const TopoPatch& source = m_topo[src];
            Box overlap{};
            for (std::size_t axis = 0; axis < 3U; ++axis) {
              const std::int64_t lo = (source.box.lo[axis] + shift[axis]) / 2;
              const std::int64_t hi =
                  (source.box.hi[axis] + shift[axis]) / 2 + (axis == normal ? 1 : 0);
              overlap.lo[axis] = std::max(global.lo[axis], lo);
              overlap.hi[axis] = std::min(global.hi[axis], hi);
            }
            if (overlap.empty()) {
              continue;
            }
            Item item{};
            item.kind = ItemKind::kGhostFaceRestricted;
            item.axis = static_cast<std::uint8_t>(normal);
            item.scale = 2U;
            item.src = src;
            for (std::size_t axis = 0; axis < 3U; ++axis) {
              item.dst_box.lo[axis] = overlap.lo[axis] - patch.box.lo[axis];
              item.dst_box.hi[axis] = overlap.hi[axis] - patch.box.lo[axis];
              item.src_offset[axis] =
                  2 * patch.box.lo[axis] - source.box.lo[axis] - shift[axis];
            }
            add_item(item);
            mark(overlap);
          }
        }
        // (3) Coarser holders: face sampling (even node: coarse face; odd node:
        // mean of the coarse cell's two faces).
        for (std::int64_t lk = 0; lk < sub_n[2]; ++lk) {
          for (std::int64_t lj = 0; lj < sub_n[1]; ++lj) {
            for (std::int64_t li = 0; li < sub_n[0]; ++li) {
              if (marked[static_cast<std::size_t>(li + sub_n[0] * (lj + sub_n[1] * lk))] != 0U) {
                continue;
              }
              const std::array<std::int64_t, 3> face{global.lo[0] + li, global.lo[1] + lj,
                                                     global.lo[2] + lk};
              if (level == m_domain.base_level) {
                throwTopology(MhdErrorCode::kStaleGeometry,
                              "base-level ghost face has no holder leaf (tiling gap)");
              }
              std::array<std::int64_t, 3> wrapped{};
              for (std::size_t axis = 0; axis < 3U; ++axis) {
                wrapped[axis] = (axis != normal && periodicAxis(axis))
                    ? wrapIndex(face[axis], extent[axis])
                    : face[axis];
              }
              std::array<std::int64_t, 3> coarse_cell{wrapped[0] / 2, wrapped[1] / 2,
                                                      wrapped[2] / 2};
              const bool even = (wrapped[normal] & 1) == 0;
              ProlongFace recipe{};
              recipe.axis = static_cast<std::uint8_t>(normal);
              recipe.dst = {static_cast<std::ptrdiff_t>(face[0] - patch.box.lo[0]),
                            static_cast<std::ptrdiff_t>(face[1] - patch.box.lo[1]),
                            static_cast<std::ptrdiff_t>(face[2] - patch.box.lo[2])};
              std::array<std::int64_t, 3> holder_wrapped{};
              std::int64_t holder = -1;
              if (even) {
                // Coarse face node N/2 is held by the coarse cell above or below.
                coarse_cell[normal] = wrapped[normal] / 2;
                holder = findCellHolder(level - 1, coarse_cell, holder_wrapped);
                std::int64_t node_local = 0;
                if (holder < 0) {
                  std::array<std::int64_t, 3> below = coarse_cell;
                  below[normal] -= 1;
                  holder = findCellHolder(level - 1, below, holder_wrapped);
                  if (holder >= 0) {
                    node_local = 1;
                  }
                }
                if (holder < 0) {
                  throwTopology(MhdErrorCode::kStaleGeometry,
                                "ghost face has no same-level, finer or parent holder (2:1 "
                                "balance violated)");
                }
                const Box& box = m_topo[static_cast<std::size_t>(holder)].box;
                std::array<std::int64_t, 3> local{holder_wrapped[0] - box.lo[0],
                                                  holder_wrapped[1] - box.lo[1],
                                                  holder_wrapped[2] - box.lo[2]};
                local[normal] += node_local;
                recipe.first = sample(ItemKind::kSampleFaceDirect, static_cast<std::uint8_t>(normal),
                                      static_cast<std::uint32_t>(holder), local);
              } else {
                holder = findCellHolder(level - 1, coarse_cell, holder_wrapped);
                if (holder < 0) {
                  throwTopology(MhdErrorCode::kStaleGeometry,
                                "ghost face has no same-level, finer or parent holder (2:1 "
                                "balance violated)");
                }
                const Box& box = m_topo[static_cast<std::size_t>(holder)].box;
                std::array<std::int64_t, 3> local{holder_wrapped[0] - box.lo[0],
                                                  holder_wrapped[1] - box.lo[1],
                                                  holder_wrapped[2] - box.lo[2]};
                recipe.first = sample(ItemKind::kSampleFaceDirect, static_cast<std::uint8_t>(normal),
                                      static_cast<std::uint32_t>(holder), local);
                local[normal] += 1;
                recipe.second = sample(ItemKind::kSampleFaceDirect,
                                       static_cast<std::uint8_t>(normal),
                                       static_cast<std::uint32_t>(holder), local);
              }
              plan.prolong_faces.push_back(recipe);
            }
          }
        }
      }
    }
  }

  // Import layout for items whose source owner differs from the destination
  // owner (deterministic: item order).
  for (Item& item : plan.items) {
    if (m_topo[item.src].owner_rank != patch.owner_rank) {
      item.import_offset = plan.import_doubles;
      plan.import_doubles += itemValueCount(item);
    }
  }
}

// ---------------------------------------------------------------------------
// Same-level owner reconciliation plans
// ---------------------------------------------------------------------------

void MhdLevelStepper::planReconciliation() {
  m_face_self_items.clear();
  m_face_cross_items.clear();
  m_edge_self_items.clear();
  m_edge_cross_items.clear();
  m_record_self_items.clear();
  m_record_cross_items.clear();
  std::vector<std::pair<std::uint32_t, std::array<std::int64_t, 3>>> hits;
  // Shared items between destination `dst` and source image (src, shift);
  // the source must precede the destination in the owner order.
  const auto pair_items = [&](std::uint32_t dst, std::uint32_t src,
                              const std::array<std::int64_t, 3>& shift) {
    const TopoPatch& p = m_topo[dst];
    const TopoPatch& q = m_topo[src];
    const bool self = dst == src;
    if (self) {
      // Only images with a lexicographically positive shift hold a lower
      // local coordinate (the owner) of the same logical item.
      if (shift == std::array<std::int64_t, 3>{0, 0, 0}) {
        return;
      }
      const bool positive = shift[0] > 0 || (shift[0] == 0 && shift[1] > 0) ||
          (shift[0] == 0 && shift[1] == 0 && shift[2] > 0);
      if (!positive) {
        return;
      }
    } else if (!(q.patch_id < p.patch_id)) {
      return;
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      // Faces of normal `axis`: closed node range along axis, cells otherwise.
      Box face{};
      bool face_ok = true;
      Box edge{};
      bool edge_ok = true;
      for (std::size_t d = 0; d < 3U; ++d) {
        const std::int64_t q_lo = q.box.lo[d] + shift[d];
        const std::int64_t q_hi = q.box.hi[d] + shift[d];
        const std::int64_t face_extra = d == axis ? 1 : 0;
        face.lo[d] = std::max(p.box.lo[d], q_lo);
        face.hi[d] = std::min(p.box.hi[d] + face_extra, q_hi + face_extra);
        face_ok = face_ok && face.lo[d] < face.hi[d];
        const std::int64_t edge_extra = d == axis ? 0 : 1;
        edge.lo[d] = std::max(p.box.lo[d], q_lo);
        edge.hi[d] = std::min(p.box.hi[d] + edge_extra, q_hi + edge_extra);
        edge_ok = edge_ok && edge.lo[d] < edge.hi[d];
      }
      const auto make = [&](ItemKind kind, const Box& global) {
        Item item{};
        item.kind = kind;
        item.axis = static_cast<std::uint8_t>(axis);
        item.src = src;
        item.dst = dst;
        for (std::size_t d = 0; d < 3U; ++d) {
          item.dst_box.lo[d] = global.lo[d] - p.box.lo[d];
          item.dst_box.hi[d] = global.hi[d] - p.box.lo[d];
          item.src_offset[d] = p.box.lo[d] - q.box.lo[d] - shift[d];
        }
        return item;
      };
      if (face_ok) {
        // A shared face lies on a boundary plane of both holders.
        (self ? m_face_self_items : m_face_cross_items).push_back(make(ItemKind::kFaceOwner, face));
        (self ? m_record_self_items : m_record_cross_items)
            .push_back(make(ItemKind::kGasRecordOwner, face));
      }
      if (edge_ok) {
        (self ? m_edge_self_items : m_edge_cross_items).push_back(make(ItemKind::kEdgeOwner, edge));
      }
    }
  };
  for (const PatchRuntime& runtime : m_local) {
    const std::uint32_t local = runtime.topo;
    const TopoPatch& patch = m_topo[local];
    Box query = patch.box;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      query.lo[axis] -= 1;
      query.hi[axis] += 1;
    }
    queryLevel(patch.level, query, hits);
    for (const auto& [other, shift] : hits) {
      // Local destination, any source (local or remote).
      pair_items(local, other, shift);
      // Local source, remote destination (for packing). Destination image
      // shift is the inverse image.
      if (other != local && m_topo[other].local == k_npos) {
        const std::array<std::int64_t, 3> inverse{-shift[0], -shift[1], -shift[2]};
        pair_items(other, local, inverse);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Coarse <- fine interface plans (leaf tiling: interfaces lie on patch
// boundaries; the fine side is the single authority, REF 5.2/7.3)
// ---------------------------------------------------------------------------

void MhdLevelStepper::planCrossLevel() {
  m_face_fine_items.clear();
  m_record_fine_items.clear();
  m_edge_fine_items.assign(static_cast<std::size_t>(std::max(m_level_count, 1)), {});
  m_interface_face_count = 0U;
  if (m_level_count <= 1) {
    return;
  }
  std::vector<std::pair<std::uint32_t, std::array<std::int64_t, 3>>> hits;
  // Coarse destination `dst` (level L-1) and fine source image (src, fine
  // shift at level L, even).
  const auto pair_items = [&](std::uint32_t dst, std::uint32_t src,
                              const std::array<std::int64_t, 3>& fine_shift) {
    const TopoPatch& c = m_topo[dst];
    const TopoPatch& f = m_topo[src];
    std::array<std::int64_t, 3> flo{};
    std::array<std::int64_t, 3> fhi{};
    for (std::size_t d = 0; d < 3U; ++d) {
      flo[d] = (f.box.lo[d] + fine_shift[d]) / 2;
      fhi[d] = (f.box.hi[d] + fine_shift[d]) / 2;
    }
    const auto slot = static_cast<std::size_t>(c.level - m_domain.base_level);
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      Box face{};
      bool face_ok = true;
      Box edge{};
      bool edge_ok = true;
      for (std::size_t d = 0; d < 3U; ++d) {
        const std::int64_t face_extra = d == axis ? 1 : 0;
        face.lo[d] = std::max(c.box.lo[d], flo[d]);
        face.hi[d] = std::min(c.box.hi[d] + face_extra, fhi[d] + face_extra);
        face_ok = face_ok && face.lo[d] < face.hi[d];
        const std::int64_t edge_extra = d == axis ? 0 : 1;
        edge.lo[d] = std::max(c.box.lo[d], flo[d]);
        edge.hi[d] = std::min(c.box.hi[d] + edge_extra, fhi[d] + edge_extra);
        edge_ok = edge_ok && edge.lo[d] < edge.hi[d];
      }
      const auto make = [&](ItemKind kind, const Box& global) {
        Item item{};
        item.kind = kind;
        item.axis = static_cast<std::uint8_t>(axis);
        item.scale = 2U;
        item.src = src;
        item.dst = dst;
        for (std::size_t d = 0; d < 3U; ++d) {
          item.dst_box.lo[d] = global.lo[d] - c.box.lo[d];
          item.dst_box.hi[d] = global.hi[d] - c.box.lo[d];
          item.src_offset[d] = 2 * c.box.lo[d] - f.box.lo[d] - fine_shift[d];
        }
        return item;
      };
      if (face_ok) {
        m_face_fine_items.push_back(make(ItemKind::kFaceFine, face));
        m_record_fine_items.push_back(make(ItemKind::kGasRecordFine, face));
        if (c.local != k_npos) {
          m_interface_face_count += static_cast<std::uint64_t>(face.volume());
        }
      }
      if (edge_ok) {
        m_edge_fine_items[slot].push_back(make(ItemKind::kEdgeFine, edge));
      }
    }
  };
  for (const PatchRuntime& runtime : m_local) {
    const std::uint32_t local = runtime.topo;
    const TopoPatch& patch = m_topo[local];
    // Local coarse destination: fine sources at level L+1.
    if (patch.level < m_max_level) {
      Box query{};
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        query.lo[axis] = 2 * patch.box.lo[axis] - 2;
        query.hi[axis] = 2 * patch.box.hi[axis] + 2;
      }
      queryLevel(patch.level + 1, query, hits);
      for (const auto& [fine, shift] : hits) {
        pair_items(local, fine, shift);
      }
    }
    // Local fine source, remote coarse destination at level L-1.
    if (patch.level > m_domain.base_level) {
      Box query{};
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        query.lo[axis] = patch.box.lo[axis] / 2 - 1;
        query.hi[axis] = patch.box.hi[axis] / 2 + 1;
      }
      queryLevel(patch.level - 1, query, hits);
      for (const auto& [coarse, coarse_shift] : hits) {
        if (m_topo[coarse].local != k_npos) {
          continue;  // generated from the coarse side above
        }
        // Coarse image at coarse_shift <=> fine image at -2 * coarse_shift.
        const std::array<std::int64_t, 3> fine_shift{-2 * coarse_shift[0], -2 * coarse_shift[1],
                                                     -2 * coarse_shift[2]};
        pair_items(coarse, local, fine_shift);
      }
    }
  }
}

void MhdLevelStepper::finalizeItemLists() {
  const auto order = [&](const Item& lhs, const Item& rhs) {
    const std::uint64_t lhs_src = ~m_topo[lhs.src].patch_id;
    const std::uint64_t rhs_src = ~m_topo[rhs.src].patch_id;
    return std::tie(lhs.dst, lhs.kind, lhs.axis, lhs_src, lhs.src_offset, lhs.dst_box.lo) <
        std::tie(rhs.dst, rhs.kind, rhs.axis, rhs_src, rhs.src_offset, rhs.dst_box.lo);
  };
  // Self items: within one destination apply the largest shift last (its
  // source is the lexicographically smallest local holder, i.e. the owner).
  const auto self_order = [&](const Item& lhs, const Item& rhs) {
    const std::array<std::int64_t, 3> lhs_shift{-lhs.src_offset[0], -lhs.src_offset[1],
                                                -lhs.src_offset[2]};
    const std::array<std::int64_t, 3> rhs_shift{-rhs.src_offset[0], -rhs.src_offset[1],
                                                -rhs.src_offset[2]};
    return std::tie(lhs.dst, lhs.kind, lhs.axis, lhs_shift, lhs.dst_box.lo) <
        std::tie(rhs.dst, rhs.kind, rhs.axis, rhs_shift, rhs.dst_box.lo);
  };
  for (std::vector<Item>* list : {&m_face_self_items, &m_edge_self_items, &m_record_self_items}) {
    std::sort(list->begin(), list->end(), self_order);
  }
  for (std::vector<Item>* list : {&m_face_cross_items, &m_edge_cross_items, &m_record_cross_items,
                                  &m_face_fine_items, &m_record_fine_items}) {
    std::sort(list->begin(), list->end(), order);
    list->erase(std::unique(list->begin(), list->end(),
                            [](const Item& lhs, const Item& rhs) {
                              return lhs.kind == rhs.kind && lhs.axis == rhs.axis &&
                                  lhs.src == rhs.src && lhs.dst == rhs.dst &&
                                  lhs.src_offset == rhs.src_offset &&
                                  lhs.dst_box.lo == rhs.dst_box.lo &&
                                  lhs.dst_box.hi == rhs.dst_box.hi;
                            }),
                list->end());
  }
  for (std::vector<Item>& list : m_edge_fine_items) {
    std::sort(list.begin(), list.end(), order);
  }
  const auto energy_variant = [](const std::vector<Item>& gas, ItemKind kind,
                                 std::vector<Item>& out) {
    out = gas;
    for (Item& item : out) {
      item.kind = kind;
    }
  };
  energy_variant(m_record_self_items, ItemKind::kEnergyRecordOwner, m_energy_self_items);
  energy_variant(m_record_cross_items, ItemKind::kEnergyRecordOwner, m_energy_cross_items);
  energy_variant(m_record_fine_items, ItemKind::kEnergyRecordFine, m_energy_fine_items);
  // Halo exchange items are in (destination, plan order) already.
}

}  // namespace cosmosim::mhd
