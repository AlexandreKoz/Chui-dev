#pragma once

// C04-02: chronology hooks between the T01 source owner and the hydro shadow
// diagnostic. The source owner calls the observer at fixed points of every
// source dispatch so the sampled-region PRODUCTION state can be snapshotted
// (read-only) at the boundaries of each smooth sub-operator, and lends its
// private smooth-operator replay at the end of the post half.
//
// Ordering of one accepted Strang interval (non-recursive AMR / fixed grid):
//   pre half : kHalfBegin -> transport -> kAfterTransport -> thermochemistry -> kHalfEnd
//   hydro    : (hydro stage; the hydro owner records the post-hydro point)
//   post half: kHalfBegin -> thermochemistry -> kAfterThermochemistry ->
//              transport -> kAfterTransport -> discrete events -> kHalfEnd
// Recursive AMR: the global halves skip the smooth operators (levels own
// them) and every level clock reports the same boundaries through
// onLevelSourceShadowBoundary (no discrete events inside a level window).
//
// Contract: the observer never mutates SimulationState, source ledgers,
// controller state, RNG or event queues. Every call is made on every rank in
// the same order (the source dispatch is rank-uniform), so the observer may
// perform rank-uniform collectives (selection, publication). Local failures
// inside the observer are contained and counted; they never veto the step.

#include <cstdint>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/hydro/hydro_shadow_reference.hpp"

namespace cosmosim::workflows {

enum class HydroShadowSourceBoundary : std::uint8_t {
  kHalfBegin = 0,
  kAfterTransport = 1,
  kAfterThermochemistry = 2,
  kHalfEnd = 3,
};

class HydroShadowChronologyObserver {
 public:
  virtual ~HydroShadowChronologyObserver() = default;
  // recursive_level_sources: the AMR level clocks own the smooth halves
  // (global halves skipped); the full-interval private chronology is then
  // not reconstructible from these boundaries and is refused with a counted
  // reason. replay is non-null only at kHalfEnd of the post half (and of the
  // unsplit compatibility dispatch) when the owner has a smooth operator; it
  // is valid for the duration of the call only.
  virtual void onSourceShadowBoundary(
      core::StepContext& context, core::SourceSplitPhase phase,
      HydroShadowSourceBoundary boundary, bool recursive_level_sources,
      const hydro::HydroShadowSmoothSourceReplay* replay) = 0;
  // T02 recursive AMR level clocks: the same boundaries for the level-local
  // Strang halves S_L(dt_L/2) H_L(dt_L) S_L(dt_L/2) of `level`, called on
  // every rank for every level clock (local work only: no collective may be
  // issued from here). level_context carries the LEVEL interval (timeline
  // begin/end/scale factors), so a replay lent at the post-half kHalfEnd
  // integrates fractions of that level interval. Default: not observed.
  virtual void onLevelSourceShadowBoundary(
      core::StepContext& level_context, int level, core::SourceSplitPhase phase,
      HydroShadowSourceBoundary boundary, const hydro::HydroShadowSmoothSourceReplay* replay) {
    (void)level_context;
    (void)level;
    (void)phase;
    (void)boundary;
    (void)replay;
  }
};

}  // namespace cosmosim::workflows
