#pragma once

// CHUI patch 58 (P1-GEOMETRY): collective, memory-admitted AMR leaf-tiling
// assurance on the production path.
//
// The exact proof (signed-corner coverage measure, see
// include/cosmosim/amr/amr_geometry_coverage.hpp) is evaluated collectively:
//
//   C1   allreduce MAX  : did any rank's populated-patch geometry fingerprint
//                         change since the last accepted proof?  If not, the
//                         cached (identical on every rank) result is returned:
//                         steady-state cost = one O(P_local) hash + one scalar
//                         collective per call.
//   V1   vote           : every populated patch is well formed (finite,
//                         positive, overflow-free, level >= 0, cell count =
//                         dims product; legacy undeclared descriptors refused
//                         when the policy requires explicit geometry).
//   C2   allreduce SUM  : populated patch count (0 => domain undeclared).
//   C3   9 x allreduce MIN : bounding box and finest cell width per axis.
//   V2   vote           : lattice representability of every local patch,
//                         periodic-box agreement, local signed-corner
//                         aggregation (interior corners cancel locally).
//   C4   allreduce SUM  : residual corner count and lattice volume.
//   V3   vote           : memory admission of the gather (root: 2 x global
//                         residual bytes, others: local residual bytes).
//   C5   bounded gather of residual corners to rank 0 (MpiContext rounds).
//   C6   broadcast of the root verdict; every rank throws the same error.
//
// Every rank executes the same collective sequence whatever it finds locally,
// so a rank that detects malformed geometry never throws before a collective
// its peers are blocked in. Memory: local scratch 8 x 32 B per local populated
// patch plus one 64 B record (reserved as kScratchArena and released on
// return); root receive 2 x 32 B per globally gathered residual point
// (kCommunication). After local cancellation a rank retains only the corners
// of the boundary of its owned region, so the gathered volume scales with
// the summed rank-region surface vertex count, never with N_cells, and is
// bounded above by 8 x P_total. No rank ever holds all patch records.

#include <array>
#include <cstdint>
#include <string_view>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows::internal {

struct AmrGeometryAssurancePolicy {
  // Production AMR hydro: a populated patch without explicit geometry (the
  // legacy all-zero descriptor) can never be a tiling leaf and is refused.
  // Callers that only want to validate whatever geometric leaves exist (and
  // treat the all-zero descriptor as a non-AMR row) pass false.
  bool require_explicit_geometry = true;
  // Periodic hydro boundary: the leaf tiling must coincide with the periodic
  // box [0, L_axis) on every axis with L_axis > 0 (box convention of
  // amr/internal/amr_periodic_topology.hpp). A periodic sub-box hydro domain
  // has no defined wrap and is refused instead of being wrapped through a gap.
  bool periodic = false;
  std::array<double, 3> periodic_box_extent_comoving{0.0, 0.0, 0.0};
};

struct AmrGeometryAssuranceResult {
  // False when no rank owns a populated geometric patch (domain undeclared:
  // origin/extent are zero, the historical "no domain" contract).
  bool declared = false;
  std::array<double, 3> origin_comoving{0.0, 0.0, 0.0};
  std::array<double, 3> extent_comoving{0.0, 0.0, 0.0};
  std::uint64_t populated_patches = 0U;
  std::uint64_t gathered_residual_points = 0U;
  // True when this call ran the full proof; false when the cache was reused.
  bool revalidated = false;
};

// Owned by one runtime (one instance per runtime object). Not thread-safe;
// call outside OpenMP regions, on every rank, at the same phase boundary.
class AmrGeometryAssurance final {
 public:
  [[nodiscard]] AmrGeometryAssuranceResult assure(
      const core::SimulationState& state,
      const parallel::MpiContext& mpi_context,
      core::MemoryGovernor* memory_governor,
      const AmrGeometryAssurancePolicy& policy,
      std::string_view phase);

  // Forces the next assure() to re-run the full proof on this rank (and
  // therefore, through C1, on every rank).
  void invalidate() noexcept { m_valid = false; }

  [[nodiscard]] const AmrGeometryAssuranceResult& lastResult() const noexcept { return m_result; }

 private:
  bool m_valid = false;
  std::uint64_t m_fingerprint = 0U;
  AmrGeometryAssuranceResult m_result{};
};

// Local, deterministic fingerprint of the populated patch geometry and the
// policy (FNV-1a over the exact bit patterns). Exposed for tests.
[[nodiscard]] std::uint64_t amrPopulatedGeometryFingerprint(
    const core::SimulationState& state, const AmrGeometryAssurancePolicy& policy) noexcept;

}  // namespace cosmosim::workflows::internal
