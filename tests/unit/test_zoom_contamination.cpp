// CHUI Part 3 Checkpoint E — zoom contamination diagnostic (two-mass test).
// STATUS: written as source only. NOT RUN — USER-RESERVED.
// Invariants: exact counts/masses of high/low-resolution particles inside the
// sphere, periodic minimum image across the box seam, nearest low-res distance,
// and decomposition invariance (combining two partitions == whole set).
#include <array>
#include <cmath>
#include <vector>

#include "cosmosim/analysis/zoom_contamination.hpp"
#include "../support/test_require.hpp"

int main() {
  namespace an = cosmosim::analysis;
  // High-res (m = 1) clustered at the box corner seam; two low-res (m = 64).
  const std::vector<double> x{0.01, 0.99, 0.02, 0.98, 0.50, 0.10};
  const std::vector<double> y{0.01, 0.99, 0.99, 0.01, 0.50, 0.00};
  const std::vector<double> z{0.00, 0.00, 0.00, 0.00, 0.50, 0.00};
  const std::vector<double> m{1.0, 1.0, 1.0, 1.0, 64.0, 64.0};
  const std::array<double, 3> center{0.0, 0.0, 0.0};
  const std::array<double, 3> box{1.0, 1.0, 1.0};
  const auto all = an::zoomContaminationDiagnostics(x, y, z, m, center, 0.2, 1.0, 1.0e-6, true, box);
  CHUI_REQUIRE(all.high_res_count == 4U);  // only correct with the periodic image
  CHUI_REQUIRE(all.low_res_count == 1U && all.low_res_mass == 64.0);
  CHUI_REQUIRE(std::abs(all.nearest_low_res_distance - 0.1) < 1.0e-15);
  CHUI_REQUIRE(std::abs(all.lowResMassFraction() - 64.0 / 68.0) < 1.0e-15);
  const std::size_t split = 3U;
  const auto a = an::zoomContaminationDiagnostics(
      std::span<const double>(x).first(split), std::span<const double>(y).first(split),
      std::span<const double>(z).first(split), std::span<const double>(m).first(split), center, 0.2, 1.0, 1.0e-6,
      true, box);
  const auto b = an::zoomContaminationDiagnostics(
      std::span<const double>(x).subspan(split), std::span<const double>(y).subspan(split),
      std::span<const double>(z).subspan(split), std::span<const double>(m).subspan(split), center, 0.2, 1.0, 1.0e-6,
      true, box);
  const auto merged = an::combine(a, b);
  CHUI_REQUIRE(merged.high_res_count == all.high_res_count && merged.low_res_count == all.low_res_count);
  CHUI_REQUIRE(merged.nearest_low_res_distance == all.nearest_low_res_distance);
  return 0;
}
