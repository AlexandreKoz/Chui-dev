#pragma once

// CHUI-EXT-PHYS-REF-v1 S03 (P1): lightcone shell geometry, comoving distance
// tables, replication / missing-mode metadata and temporal interpolation.
// All distances come from the S01 background authority.
//
// Shells: snapshots i (ascending a_i, i.e. descending chi_i) own the radial
// shell between the comoving-distance midpoints to their neighbours; the
// outermost/innermost shells are closed by chi(a_min_snapshot) and 0 (or the
// configured chi_min). Inside a shell, emission scale factor along a ray is
// a(chi) (exact inversion of the table), and fields may be interpolated in
// time between the two bracketing snapshots with weights linear in a
// (TemporalWeights) — or taken from the shell's own snapshot (piecewise
// constant), an explicit, recorded choice.
// Replication: a periodic box of side L tiles space; a ray to chi_max crosses
// ceil(chi_max / L) replicas; modes with k < 2 pi / L are absent and
// structures repeat along the line of sight when replicas > 1 — both are
// recorded, never hidden. Full-sky statistics from a small replicated box are
// not claimed.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "cosmosim/core/cosmology.hpp"

namespace cosmosim::lightcone {

class ComovingDistanceTable {
 public:
  // chi(a) [m] tabulated on `points` nodes uniform in ln a over [a_min, 1].
  ComovingDistanceTable(const core::LambdaCdmBackground& background, double a_min, std::size_t points = 4096U);
  [[nodiscard]] double comovingDistanceM(double scale_factor) const;  // throws outside [a_min, 1]
  [[nodiscard]] double scaleFactorAt(double comoving_distance_m) const;
  [[nodiscard]] double aMin() const noexcept { return m_a_min; }
  [[nodiscard]] double chiMaxM() const noexcept { return m_chi.front(); }

 private:
  double m_a_min = 0.0;
  std::vector<double> m_ln_a;  // ascending
  std::vector<double> m_chi;   // descending with ln a
};

struct LightconeShell {
  std::size_t snapshot_index = 0;
  double snapshot_scale_factor = 0.0;
  double chi_inner_m = 0.0;
  double chi_outer_m = 0.0;
  double a_inner = 0.0;  // larger a (closer)
  double a_outer = 0.0;
};

struct ReplicationMetadata {
  double box_size_m = 0.0;
  double chi_max_m = 0.0;
  std::uint64_t replicas_along_ray = 0;
  bool repeated_structures = false;
  double missing_mode_k_max_per_m = 0.0;  // 2 pi / L: modes below are absent
  std::string strategy = "periodic_tiling_fixed_orientation";
};

[[nodiscard]] std::vector<LightconeShell> planLightconeShells(const ComovingDistanceTable& table,
                                                              const std::vector<double>& snapshot_scale_factors,
                                                              double chi_min_m = 0.0);
[[nodiscard]] ReplicationMetadata replicationMetadata(double box_size_m, double chi_max_m);

struct TemporalWeights {
  std::size_t lower_snapshot = 0;
  std::size_t upper_snapshot = 0;
  double upper_weight = 0.0;  // value = (1 - w) lower + w upper
};
// Bracketing snapshots of a_emit (ascending list); throws outside the range.
[[nodiscard]] TemporalWeights temporalWeights(const std::vector<double>& snapshot_scale_factors, double a_emit);

// Observer-frame redshift of a source at cosmological redshift z_cos moving
// with line-of-sight peculiar velocity v_los (positive = receding):
// 1 + z_obs = (1 + z_cos)(1 + v_los / c) (first order in v/c, declared).
[[nodiscard]] double observedRedshift(double z_cos, double v_los_m_s);
// Snapshot (fixed-time) sightline velocity coordinate: Hubble velocity of a
// comoving offset dx at scale factor a, v = a H(a) dx (proper).
[[nodiscard]] double hubbleVelocityMS(const core::LambdaCdmBackground& background, double a, double comoving_offset_m);

}  // namespace cosmosim::lightcone
