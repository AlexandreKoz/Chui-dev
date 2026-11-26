#pragma once

// CHUI-RT-BH-AGN-REF-v1, RT companion sections 5 and 7 (RT02, RT06):
// level-synchronous finite-volume M1 transport over the authoritative leaf
// radiation records.
//
// Discretization: dU/dt = -(1/a) sum_faces A F_hat / V (comoving A, V; the
// factor 1/a converts the comoving divergence to physical; dt is proper).
// Transport variables (E, G = F c_hat/c) with speed c_hat: the reduced-speed
// system is the standard M1 system in (E, G) (RT companion section 6), so the
// physical flux F is recovered exactly as G c/c_hat.
// Reconstruction: unsplit piecewise-linear MC per axis with JOINT limiting of
// (E, G) — one theta per cell and axis so both face states remain in the
// admissible cone (no independent clipping of E or F). First-order variant.
// Riemann: HLL with exact M1 eigenvalue bounds (radiation_closure), GLF
// (global safe bound c_hat) per-face fallback when HLL is non-finite.
// Faces: each face flux is computed by every adjacent record from bitwise
// identical inputs (same code path, same ghost data), hence equal and
// opposite; at coarse/fine interfaces the coarse face flux is REPLACED by the
// area sum of the fine face fluxes (fine authority, per stage: reflux).
// Ghosts (axis slabs, width 2, no corners needed): same level copy; finer ->
// volume average of the ratio-2 children (admissible by convexity); coarser ->
// piecewise-constant prolongation (admissible); domain boundary: periodic
// wrap, outflow (zero gradient) or reflecting (normal G negated). Boundary
// conditions are the radiation's own declaration (physics.rt_boundary),
// independent of gravity.
// Stage admissibility: after a stage every cell must satisfy E >= 0,
// |G| <= c_hat E; otherwise the WHOLE stage is recomputed first-order GLF
// (realizable for the sum-of-speeds CFL <= 1/2) and counted; if still
// inadmissible the stage is rejected (caller rolls back and halves dt).
// Scope (pass 01): single-rank records (multi-rank halo exchange is not
// implemented; admission refuses rt with > 1 MPI rank).

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/config_rt_bh_agn.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::radiation {

struct RadiationTransportSettings {
  core::RtRiemannSolverConfig riemann = core::RtRiemannSolverConfig::kHll;
  core::RtReconstructionConfig reconstruction = core::RtReconstructionConfig::kPlmMc;
  core::RtBoundaryConfig boundary = core::RtBoundaryConfig::kPeriodic;
  std::array<double, 3> domain_origin_comoving{0.0, 0.0, 0.0};
  std::array<double, 3> domain_extent_comoving{0.0, 0.0, 0.0};
  double c_hat_code = 0.0;          // transport speed (c_hat), code velocity
};

// Working moments (E, G) for all records, layout identical to the canonical
// lanes (group-major per record).
struct RadiationWorkState {
  std::vector<double> e;
  std::vector<double> gx;
  std::vector<double> gy;
  std::vector<double> gz;
  void resize(std::size_t n) {
    e.assign(n, 0.0);
    gx.assign(n, 0.0);
    gy.assign(n, 0.0);
    gz.assign(n, 0.0);
  }
  [[nodiscard]] std::size_t size() const noexcept { return e.size(); }
};

struct RadiationStageDiagnostics {
  std::uint64_t hll_fallback_faces = 0U;
  std::uint64_t first_order_stage_fallbacks = 0U;
  std::uint64_t inadmissible_cells = 0U;
  double boundary_energy_outflow_code = 0.0;  // integral of E-flux leaving the domain (comoving units * dt)
};

// Ghost-source map and reflux structure of the leaf record tiling.
class RadiationTopology {
 public:
  RadiationTopology(const core::RadiationMomentSoa& moments, const RadiationTransportSettings& settings);

  struct GhostSource {
    enum class Kind : std::uint8_t { kSameLevel = 0, kCoarser = 1, kFiner = 2, kOutflow = 3, kReflecting = 4 };
    Kind kind = Kind::kOutflow;
    std::uint32_t record = 0U;
    std::array<std::uint32_t, 8> cells{};  // kFiner: 8 children; otherwise cells[0]
    std::uint8_t cell_count = 1U;
    std::int8_t reflect_axis = -1;
  };
  // Ghost slab index: record, axis (0..2), side (0 = low, 1 = high), layer
  // (0 = adjacent, 1 = second), transverse (u, v) of the slab.
  [[nodiscard]] const GhostSource& ghost(std::uint32_t record, int axis, int side, int layer, std::uint32_t u,
                                         std::uint32_t v) const;
  // Coarse faces whose neighbor is finer: (record, cell, axis, side) -> true.
  [[nodiscard]] bool faceIsCoarseFine(std::uint32_t record, std::uint32_t cell, int axis, int side) const;
  [[nodiscard]] std::size_t recordCount() const noexcept { return m_records.size(); }

 private:
  friend class RadiationLevelSolver;
  struct RecordInfo {
    core::RadiationPatchRecord record{};
    std::array<std::uint64_t, 3> slab_offset{};  // into m_ghosts per axis (both sides, both layers)
  };
  [[nodiscard]] std::int64_t findContainingRecord(const std::array<double, 3>& point, int preferred_level) const;
  std::vector<RecordInfo> m_records;
  std::vector<GhostSource> m_ghosts;
  std::vector<std::uint64_t> m_coarse_fine_faces;  // sorted keys
  RadiationTransportSettings m_settings;
};

// Explicit transport operator T(U) (flux divergence) on the full record set.
class RadiationLevelSolver {
 public:
  RadiationLevelSolver(const core::RadiationMomentSoa& moments, const RadiationTopology& topology,
                       RadiationTransportSettings settings);
  // rate = T(U) (per proper time). Returns false if any flux is non-finite.
  // `inverse_scale_factor` = 1/a at the stage time. `first_order` forces
  // first-order GLF everywhere (stage fallback).
  [[nodiscard]] bool transportRate(const RadiationWorkState& state, double inverse_scale_factor, bool first_order,
                                   RadiationWorkState& rate, RadiationStageDiagnostics& diagnostics) const;
  // Stable proper timestep: cfl / max_cells(sum_d c_hat / dx_phys_d).
  [[nodiscard]] double stableTimeStep(double cfl, double scale_factor) const;
  [[nodiscard]] bool admissible(const RadiationWorkState& state, std::uint64_t& bad_cells) const;
  [[nodiscard]] std::uint32_t groups() const noexcept { return m_groups; }

 private:
  const core::RadiationMomentSoa& m_moments;
  const RadiationTopology& m_topology;
  RadiationTransportSettings m_settings;
  std::uint32_t m_groups = 0U;
};

}  // namespace cosmosim::radiation
