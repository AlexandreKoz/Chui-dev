#pragma once

// CHUI Part 3 Checkpoint E — zoom-in contamination diagnostic (G01 / zoom_in).
//
// A zoom region is clean only by MEASUREMENT, never by IC metadata: this
// diagnostic classifies collisionless particles by mass (high resolution =
// mass <= min_mass x (1 + high_res_mass_tolerance), where min_mass is the
// global minimum supplied by the caller after an MPI min-reduction) and
// reports, inside a sphere (periodic minimum image when requested) around the
// target-halo centre: counts and masses of high- and low-resolution
// particles, the low-resolution mass fraction and the distance of the nearest
// low-resolution particle. Rank-local partial results combine by sums (counts,
// masses) and min (nearest distance); `combine` implements exactly that.
// Gas cells are not classified (they share the high-resolution region by
// construction of the refinement criteria; contamination is a collisionless
// statement). Status: IMPLEMENTED IN SOURCE; emitted at the run-health
// cadence by AnalysisRuntime (event analysis.zoom_contamination, G3.4).
// EXPERIMENTAL/UNQUALIFIED until owner-run zoom validation.

#include <array>
#include <cstdint>
#include <limits>
#include <span>

namespace cosmosim::analysis {

struct ZoomContaminationReport {
  std::uint64_t high_res_count = 0U;
  std::uint64_t low_res_count = 0U;
  double high_res_mass = 0.0;
  double low_res_mass = 0.0;
  double nearest_low_res_distance = std::numeric_limits<double>::infinity();
  [[nodiscard]] double lowResMassFraction() const noexcept {
    const double total = high_res_mass + low_res_mass;
    return total > 0.0 ? low_res_mass / total : 0.0;
  }
};

[[nodiscard]] ZoomContaminationReport zoomContaminationDiagnostics(
    std::span<const double> x, std::span<const double> y, std::span<const double> z,
    std::span<const double> mass, const std::array<double, 3>& center, double radius,
    double global_min_mass, double high_res_mass_tolerance, bool periodic,
    const std::array<double, 3>& box_extent);

[[nodiscard]] ZoomContaminationReport combine(const ZoomContaminationReport& a,
                                              const ZoomContaminationReport& b) noexcept;

}  // namespace cosmosim::analysis
