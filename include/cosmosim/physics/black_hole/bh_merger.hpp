#pragma once

// CHUI-RT-BH-AGN-REF-v1, BH companion section 8 (BH07): bound-pair mergers.
//
// Criteria (all required; proximity alone is prohibited):
//  1. proximity  r_proper < f_merge * eps_proper (eps = larger pair softening),
//  2. bound in the softened two-body potential consistent with the gravity
//     softening scale: E = mu |dv|^2/2 - G M1 M2 / sqrt(r^2 + eps^2) < 0
//     (Plummer-equivalent; the repository's spline softening differs inside
//     eps by an O(1) factor, documented as an approximation of the criterion),
//  3. converging (dr . dv < 0) or unresolved (r < eps).
// Local tidal context is NOT evaluated in pass 01 (recorded gap). Pairs are
// selected deterministically and disjointly: most-bound first, ties by
// (min id, max id). The survivor is the more massive member (tie: lower id);
// it receives mass, momentum (centre of mass), disk, pending budgets and
// lineage. The Newtonian merger conserves the tracked mass and momentum; it
// does not claim a gravitational-wave remnant (no GW loss, no recoil in P0).

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/physics/black_hole/bh_engine_state.hpp"
#include "cosmosim/physics/rt_bh_agn_frame.hpp"

namespace cosmosim::physics {

struct BhMergerCandidate {
  std::uint64_t particle_id = 0U;
  std::uint32_t local_bh_row = 0U;
  std::uint32_t owner_rank = 0U;
  std::array<double, 3> position_stored{};
  std::array<double, 3> velocity_peculiar{};
  double dynamical_mass_code = 0.0;
  double softening_proper = 0.0;
};

struct BhMergerPair {
  std::uint64_t survivor_id = 0U;
  std::uint64_t absorbed_id = 0U;
  std::uint32_t survivor_index = 0U;  // indices into the candidate span
  std::uint32_t absorbed_index = 0U;
  double binding_energy_code = 0.0;
  double separation_proper = 0.0;
};

struct BhMergerCriteria {
  double radius_softening_factor = 2.0;
  double newton_g_code = 0.0;
};

[[nodiscard]] bool bhPairMergeEligible(const BhMergerCandidate& a, const BhMergerCandidate& b,
                                       const BhMergerCriteria& criteria, const ProperFrameContext& frame,
                                       const PeriodicBox& box, double* binding_energy_out = nullptr,
                                       double* separation_out = nullptr) noexcept;

// O(n^2) over the local candidate set (BH counts are sparse; candidates are
// pre-filtered by the caller's proximity buckets). Deterministic, disjoint.
[[nodiscard]] std::vector<BhMergerPair> selectDisjointMergerPairs(std::span<const BhMergerCandidate> candidates,
                                                                  const BhMergerCriteria& criteria,
                                                                  const ProperFrameContext& frame,
                                                                  const PeriodicBox& box);

// Combine engine records (survivor receives everything; extensive sums).
[[nodiscard]] BlackHoleEngineRecord mergeBlackHoleEngineRecords(const BlackHoleEngineRecord& survivor,
                                                                const BlackHoleEngineRecord& absorbed,
                                                                std::uint64_t absorbed_id);

// P1-BH03 (binary delay building block): Peters (1964) circular-orbit GW
// inspiral time t = (5/256) c^5 a^4 / (G^3 m1 m2 (m1 + m2)). Eccentric
// enhancement, stellar/gas hardening and NR remnant/recoil fits are NOT
// implemented (no pinned fits; FAIL_CLOSED for binary-delay admission).
[[nodiscard]] double petersCircularInspiralTime(double separation, double m1, double m2, double newton_g,
                                                double speed_of_light) noexcept;

}  // namespace cosmosim::physics
