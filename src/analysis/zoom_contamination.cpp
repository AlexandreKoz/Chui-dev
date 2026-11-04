#include "cosmosim/analysis/zoom_contamination.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cosmosim::analysis {

ZoomContaminationReport zoomContaminationDiagnostics(
    std::span<const double> x, std::span<const double> y, std::span<const double> z,
    std::span<const double> mass, const std::array<double, 3>& center, double radius,
    double global_min_mass, double high_res_mass_tolerance, bool periodic,
    const std::array<double, 3>& box_extent) {
  if (x.size() != y.size() || x.size() != z.size() || x.size() != mass.size()) {
    throw std::invalid_argument("zoomContaminationDiagnostics: lane extents differ");
  }
  if (!(radius > 0.0) || !(global_min_mass > 0.0) || !(high_res_mass_tolerance >= 0.0)) {
    throw std::invalid_argument("zoomContaminationDiagnostics: radius, min mass must be > 0, tolerance >= 0");
  }
  const double high_res_limit = global_min_mass * (1.0 + high_res_mass_tolerance);
  ZoomContaminationReport report;
  for (std::size_t i = 0; i < x.size(); ++i) {
    std::array<double, 3> d{x[i] - center[0], y[i] - center[1], z[i] - center[2]};
    if (periodic) {
      for (std::size_t c = 0; c < 3U; ++c) {
        if (box_extent[c] > 0.0) d[c] -= box_extent[c] * std::round(d[c] / box_extent[c]);
      }
    }
    const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    const bool high_res = mass[i] <= high_res_limit;
    if (!high_res) {
      report.nearest_low_res_distance = std::min(report.nearest_low_res_distance, r);
    }
    if (r > radius) continue;
    if (high_res) {
      ++report.high_res_count;
      report.high_res_mass += mass[i];
    } else {
      ++report.low_res_count;
      report.low_res_mass += mass[i];
    }
  }
  return report;
}

ZoomContaminationReport combine(const ZoomContaminationReport& a, const ZoomContaminationReport& b) noexcept {
  ZoomContaminationReport r;
  r.high_res_count = a.high_res_count + b.high_res_count;
  r.low_res_count = a.low_res_count + b.low_res_count;
  r.high_res_mass = a.high_res_mass + b.high_res_mass;
  r.low_res_mass = a.low_res_mass + b.low_res_mass;
  r.nearest_low_res_distance = std::min(a.nearest_low_res_distance, b.nearest_low_res_distance);
  return r;
}

}  // namespace cosmosim::analysis
