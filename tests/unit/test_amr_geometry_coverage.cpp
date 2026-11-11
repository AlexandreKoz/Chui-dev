// CHUI patch 58 (P1-GEOMETRY) — exact AMR leaf-tiling proof.
// Authored without execution: NOT RUN — USER-RESERVED.
//
// Acceptance fixtures (prompt cases 1-9): valid uniform rectangle; valid
// mixed-level dyadic refinement; disjoint patches split across ranks; overlap
// compensating a hole with equal total volume (must reject — the patch-56
// volume check accepted it); plain hole; invalid/non-finite/negative/zero
// extent; a rank owning no patches; periodic and non-periodic boundaries;
// ownership changes (restart/regrid/migration) that keep or break the tiling.

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/amr/amr_geometry_coverage.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "workflows/internal/amr_geometry_assurance.hpp"

namespace {

using cosmosim::amr::AmrPatchGeometryRecord;
using cosmosim::amr::verifyAmrPatchSetTiling;

AmrPatchGeometryRecord box(std::uint64_t id, std::int32_t level, std::uint32_t rank,
                           std::array<double, 3> origin, std::array<double, 3> extent,
                           std::array<std::uint32_t, 3> dims = {2U, 2U, 2U}) {
  AmrPatchGeometryRecord r;
  r.patch_id = id;
  r.level = level;
  r.owning_rank = rank;
  r.origin_comoving = origin;
  r.extent_comoving = extent;
  r.cell_dims = dims;
  r.cell_count = dims[0] * dims[1] * dims[2];
  return r;
}

// Eight children of [o, o + 2h)^3 at the next level.
void appendOctants(std::vector<AmrPatchGeometryRecord>& out, std::uint64_t first_id, std::int32_t level,
                   std::array<double, 3> o, double h, std::uint32_t rank) {
  for (unsigned octant = 0U; octant < 8U; ++octant) {
    out.push_back(box(first_id + octant, level, rank,
                      {o[0] + ((octant & 1U) ? h : 0.0), o[1] + ((octant & 2U) ? h : 0.0),
                       o[2] + ((octant & 4U) ? h : 0.0)},
                      {h, h, h}));
  }
}

void caseUniformRectangle() {
  const std::vector<AmrPatchGeometryRecord> patches{box(1, 0, 0, {0.0, 0.0, 0.0}, {1.0, 2.0, 4.0}, {4U, 8U, 16U})};
  const auto verdict = verifyAmrPatchSetTiling(patches);
  assert(verdict.exact_tiling);
  assert(verdict.residual_points == 0U);
}

void caseMixedLevelDyadic() {
  // Root [0,1)^3 split into 8 level-1 children; child 0 refined again into
  // 8 level-2 children (fine/coarse adjacency on three faces).
  std::vector<AmrPatchGeometryRecord> patches;
  for (unsigned octant = 1U; octant < 8U; ++octant) {
    patches.push_back(box(10 + octant, 1, 0,
                          {(octant & 1U) ? 0.5 : 0.0, (octant & 2U) ? 0.5 : 0.0, (octant & 4U) ? 0.5 : 0.0},
                          {0.5, 0.5, 0.5}));
  }
  appendOctants(patches, 100, 2, {0.0, 0.0, 0.0}, 0.25, 0);
  const auto verdict = verifyAmrPatchSetTiling(patches);
  assert(verdict.exact_tiling);
}

void caseSplitAcrossRanks() {
  // Two slabs owned by different ranks; the lattice measure is rank-agnostic.
  const std::vector<AmrPatchGeometryRecord> patches{
      box(1, 0, 0, {0.0, 0.0, 0.0}, {0.5, 1.0, 1.0}),
      box(2, 0, 1, {0.5, 0.0, 0.0}, {0.5, 1.0, 1.0})};
  assert(verifyAmrPatchSetTiling(patches).exact_tiling);
}

void caseOverlapCompensatingHole() {
  // Audit counterexample: x-slabs [0,0.5), [0.4,0.7), [0.8,1): volumes sum to
  // 1 = bounding-box volume, but [0.4,0.5) is doubly covered and [0.7,0.8) is
  // empty. Cell widths 0.1 so every face is on the lattice.
  const std::vector<AmrPatchGeometryRecord> patches{
      box(1, 0, 0, {0.0, 0.0, 0.0}, {0.5, 1.0, 1.0}, {5U, 1U, 1U}),
      box(2, 0, 0, {0.4, 0.0, 0.0}, {0.3, 1.0, 1.0}, {3U, 1U, 1U}),
      box(3, 0, 1, {0.8, 0.0, 0.0}, {0.2, 1.0, 1.0}, {2U, 1U, 1U})};
  const auto verdict = verifyAmrPatchSetTiling(patches);
  assert(!verdict.exact_tiling);
  assert(verdict.residual_points > 0U);
  assert(verdict.diagnostic.find("compensating") != std::string::npos);
}

void casePlainHole() {
  std::vector<AmrPatchGeometryRecord> patches;
  for (unsigned octant = 0U; octant < 7U; ++octant) {  // octant 7 missing
    patches.push_back(box(10 + octant, 1, 0,
                          {(octant & 1U) ? 0.5 : 0.0, (octant & 2U) ? 0.5 : 0.0, (octant & 4U) ? 0.5 : 0.0},
                          {0.5, 0.5, 0.5}));
  }
  const auto verdict = verifyAmrPatchSetTiling(patches);
  assert(!verdict.exact_tiling);
  assert(verdict.diagnostic.find("hole") != std::string::npos);
}

void casePureOverlap() {
  const std::vector<AmrPatchGeometryRecord> patches{
      box(1, 0, 0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}),
      box(2, 0, 1, {0.0, 0.0, 0.0}, {0.5, 0.5, 0.5})};  // duplicate ownership of a sub-box
  const auto verdict = verifyAmrPatchSetTiling(patches);
  assert(!verdict.exact_tiling);
  assert(verdict.diagnostic.find("overlap") != std::string::npos);
}

void caseInvalidGeometry() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<std::vector<AmrPatchGeometryRecord>> bad{
      {box(7, 0, 3, {0.0, 0.0, 0.0}, {1.0, nan, 1.0})},
      {box(7, 0, 3, {0.0, 0.0, 0.0}, {1.0, 1.0, -1.0})},
      {box(7, 0, 3, {0.0, 0.0, 0.0}, {0.0, 1.0, 1.0})},
      {box(7, 0, 3, {inf, 0.0, 0.0}, {1.0, 1.0, 1.0})},
      {box(7, 0, 3, {1.0e308, 0.0, 0.0}, {1.0e308, 1.0, 1.0})},  // upper overflows
      {box(7, -1, 3, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0})}};
  for (const auto& patches : bad) {
    const auto verdict = verifyAmrPatchSetTiling(patches);
    assert(!verdict.exact_tiling);
    assert(verdict.diagnostic.find("AMR patch 7") != std::string::npos);
    assert(verdict.diagnostic.find("rank 3") != std::string::npos);
  }
  // Cell count disagreeing with the dimensions.
  AmrPatchGeometryRecord mismatch = box(8, 0, 0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  mismatch.cell_count = 7U;
  assert(!cosmosim::amr::diagnoseAmrPatchGeometry(mismatch).empty());
  // Misaligned (non-lattice) face.
  const std::vector<AmrPatchGeometryRecord> misaligned{
      box(1, 0, 0, {0.0, 0.0, 0.0}, {0.5, 1.0, 1.0}),
      box(2, 0, 0, {0.5 + 1.0e-3, 0.0, 0.0}, {0.5 - 1.0e-3, 1.0, 1.0})};
  assert(!verifyAmrPatchSetTiling(misaligned).exact_tiling);
}

cosmosim::core::SimulationState stateFromRecords(const std::vector<AmrPatchGeometryRecord>& records) {
  cosmosim::core::SimulationState state;
  state.resizePatches(records.size());
  for (std::size_t i = 0; i < records.size(); ++i) {
    const AmrPatchGeometryRecord& r = records[i];
    state.patches.patch_id[i] = r.patch_id;
    state.patches.level[i] = r.level;
    state.patches.first_cell[i] = 0U;
    state.patches.cell_count[i] = r.cell_count;
    state.patches.owning_rank[i] = r.owning_rank;
    state.patches.parent_patch_id[i] = 0U;
    state.patches.morton_key[i] = r.patch_id;
    state.patches.origin_x_comoving[i] = r.origin_comoving[0];
    state.patches.origin_y_comoving[i] = r.origin_comoving[1];
    state.patches.origin_z_comoving[i] = r.origin_comoving[2];
    state.patches.extent_x_comoving[i] = r.extent_comoving[0];
    state.patches.extent_y_comoving[i] = r.extent_comoving[1];
    state.patches.extent_z_comoving[i] = r.extent_comoving[2];
    state.patches.cell_dim_x[i] = static_cast<std::uint16_t>(r.cell_dims[0]);
    state.patches.cell_dim_y[i] = static_cast<std::uint16_t>(r.cell_dims[1]);
    state.patches.cell_dim_z[i] = static_cast<std::uint16_t>(r.cell_dims[2]);
  }
  return state;
}

void caseEmptyRankAndAssuranceCache() {
  namespace wi = cosmosim::workflows::internal;
  const cosmosim::parallel::MpiContext serial;  // MPI disabled: collectives are local identities
  wi::AmrGeometryAssurance assurance;
  wi::AmrGeometryAssurancePolicy policy;

  // A rank owning no patches: domain undeclared, no throw.
  cosmosim::core::SimulationState empty;
  const auto none = assurance.assure(empty, serial, nullptr, policy, "test empty rank");
  assert(!none.declared && none.revalidated);

  // Mixed-level valid hierarchy plus an inactive (cell_count = 0) parent.
  std::vector<AmrPatchGeometryRecord> records;
  appendOctants(records, 100, 1, {0.0, 0.0, 0.0}, 0.5, 0);
  AmrPatchGeometryRecord parent = box(1, 0, 0, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  parent.cell_count = 0U;  // metadata-only parent never enters the tiling
  records.push_back(parent);
  cosmosim::core::SimulationState state = stateFromRecords(records);
  const auto first = assurance.assure(state, serial, nullptr, policy, "test valid");
  assert(first.declared && first.revalidated && first.populated_patches == 8U);
  assert(std::abs(first.extent_comoving[0] - 1.0) < 1.0e-15);
  const auto second = assurance.assure(state, serial, nullptr, policy, "test cached");
  assert(second.declared && !second.revalidated);

  // Migration/restart ownership change keeping the geometry: re-proved, valid.
  state.patches.owning_rank[0] = 1U;
  const auto migrated = assurance.assure(state, serial, nullptr, policy, "test migrated");
  assert(migrated.revalidated && migrated.declared);

  // Regrid/restart corruption producing an overlap: refused, cache invalid.
  state.patches.origin_x_comoving[1] = 0.25;
  bool threw = false;
  try {
    (void)assurance.assure(state, serial, nullptr, policy, "test corrupted");
  } catch (const std::runtime_error& error) {
    threw = std::string(error.what()).find("test corrupted") != std::string::npos;
  }
  assert(threw);

  // Legacy all-zero populated descriptor: refused by production policy only.
  cosmosim::core::SimulationState legacy = stateFromRecords({box(9, 0, 0, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0U, 0U, 0U})});
  legacy.patches.cell_count[0] = 4U;
  wi::AmrGeometryAssurance strict;
  threw = false;
  try {
    (void)strict.assure(legacy, serial, nullptr, policy, "test legacy strict");
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
  wi::AmrGeometryAssurancePolicy lenient = policy;
  lenient.require_explicit_geometry = false;
  wi::AmrGeometryAssurance relaxed;
  assert(!relaxed.assure(legacy, serial, nullptr, lenient, "test legacy lenient").declared);
}

void casePeriodicAndNonPeriodic() {
  namespace wi = cosmosim::workflows::internal;
  const cosmosim::parallel::MpiContext serial;
  // Legitimate isolated sub-box hydro domain inside a larger gravity box.
  const std::vector<AmrPatchGeometryRecord> sub_box{box(1, 0, 0, {2.0, 2.0, 2.0}, {1.0, 1.0, 1.0})};
  cosmosim::core::SimulationState state = stateFromRecords(sub_box);
  wi::AmrGeometryAssurancePolicy open_policy;
  wi::AmrGeometryAssurance open_assurance;
  const auto open = open_assurance.assure(state, serial, nullptr, open_policy, "test open sub-box");
  assert(open.declared && std::abs(open.origin_comoving[0] - 2.0) < 1.0e-15);

  // The same sub-box under a periodic boundary: refused (no defined wrap).
  wi::AmrGeometryAssurancePolicy periodic_policy;
  periodic_policy.periodic = true;
  periodic_policy.periodic_box_extent_comoving = {5.0, 5.0, 5.0};
  wi::AmrGeometryAssurance periodic_assurance;
  bool threw = false;
  try {
    (void)periodic_assurance.assure(state, serial, nullptr, periodic_policy, "test periodic sub-box");
  } catch (const std::runtime_error& error) {
    threw = std::string(error.what()).find("periodic box") != std::string::npos;
  }
  assert(threw);

  // Full periodic box tiling: accepted.
  const std::vector<AmrPatchGeometryRecord> full{box(1, 0, 0, {0.0, 0.0, 0.0}, {5.0, 5.0, 5.0})};
  cosmosim::core::SimulationState full_state = stateFromRecords(full);
  wi::AmrGeometryAssurance full_assurance;
  assert(full_assurance.assure(full_state, serial, nullptr, periodic_policy, "test periodic full").declared);
}

}  // namespace

int main() {
  caseUniformRectangle();
  caseMixedLevelDyadic();
  caseSplitAcrossRanks();
  caseOverlapCompensatingHole();
  casePlainHole();
  casePureOverlap();
  caseInvalidGeometry();
  caseEmptyRankAndAssuranceCache();
  casePeriodicAndNonPeriodic();
  std::cout << "test_amr_geometry_coverage: all cases passed\n";
  return 0;
}
