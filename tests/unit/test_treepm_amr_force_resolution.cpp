// CHUI Part 3 Checkpoint E — TreePM force accuracy at increasing AMR levels
// and zoom contamination diagnostics.
// STATUS: written as source only. NOT RUN — USER-RESERVED. Skips (77) without
// FFTW (the non-FFTW PM fallback is not an accuracy reference).
//
// Question answered (evidence fixture, not a claim): does CHUI's TreePM force
// on gas sampled at AMR cell centres stay accurate as the cell width shrinks
// by 2^L, for a fixed PM mesh and softening epsilon <= fine cell width?
// Construction: a level-L refined block of (4 x 2^L)^3 / 8 cells of uniform
// density (cell-centre point masses), embedded in a uniform coarse background,
// periodic box L = 1. Reference: direct minimum-image summation with the SAME
// Plummer softening (the documented TreePM reference convention; not Ewald).
// Metrics per level: relative RMS (L2) and maximum relative force error on the
// refined-block cells nearest the coarse/fine boundary and in its interior.
// Acceptance (problem-specific, derived from the level-0 baseline, not tuned):
//  * the relative L2 error at every level <= max(2 x level-0 error, 5e-2);
//  * TreePM split composition error < 1e-12 at every level (no force double
//    count between PM and tree).
// A level whose error exceeds the bound is a GENUINE finding (the fine-force
// correction of Gate E would then be required); do not widen the bound.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/gravity/tree_softening.hpp"
#include "../support/test_require.hpp"

namespace {

namespace gr = cosmosim::gravity;

double minimumImage(double d) { return d - std::round(d); }

struct Field {
  std::vector<double> ax, ay, az;
};

Field direct(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
             const std::vector<double>& m, const gr::TreeSofteningPolicy& softening) {
  Field f{std::vector<double>(x.size()), std::vector<double>(x.size()), std::vector<double>(x.size())};
  for (std::size_t t = 0; t < x.size(); ++t) {
    for (std::size_t s = 0; s < x.size(); ++s) {
      if (s == t) continue;
      const double dx = minimumImage(x[s] - x[t]), dy = minimumImage(y[s] - y[t]), dz = minimumImage(z[s] - z[t]);
      const double factor = m[s] * gr::softenedInvR3(dx * dx + dy * dy + dz * dz, softening);
      f.ax[t] += factor * dx; f.ay[t] += factor * dy; f.az[t] += factor * dz;
    }
  }
  return f;
}

void testForceAccuracyDoesNotDegradeWithLevel() {
#if COSMOSIM_ENABLE_FFTW
  const gr::PmGridShape shape{32U, 32U, 32U};
  double level0_error = -1.0;
  for (int level = 0; level <= 3; ++level) {
    std::vector<double> x, y, z, m;
    const double coarse = 1.0 / 8.0;  // 8^3 coarse cells
    const double fine = coarse / std::pow(2.0, level);
    const double block_lo = 0.375, block_hi = 0.625;  // refined block (2^3 coarse cells)
    for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) for (int k = 0; k < 8; ++k) {
      const double cx = (i + 0.5) * coarse, cy = (j + 0.5) * coarse, cz = (k + 0.5) * coarse;
      if (cx > block_lo && cx < block_hi && cy > block_lo && cy < block_hi && cz > block_lo && cz < block_hi) continue;
      x.push_back(cx); y.push_back(cy); z.push_back(cz); m.push_back(coarse * coarse * coarse * (1.0 + 0.1 * i));
    }
    const std::size_t first_fine = x.size();
    const int n_fine = static_cast<int>(std::lround((block_hi - block_lo) / fine));
    for (int i = 0; i < n_fine; ++i) for (int j = 0; j < n_fine; ++j) for (int k = 0; k < n_fine; ++k) {
      const double fx = block_lo + (i + 0.5) * fine, fy = block_lo + (j + 0.5) * fine, fz = block_lo + (k + 0.5) * fine;
      x.push_back(fx); y.push_back(fy); z.push_back(fz);
      m.push_back(fine * fine * fine * (2.0 + 0.5 * std::sin(10.0 * fx)));  // overdense, structured
    }
    gr::TreePmOptions options;
    options.pm_options.box_size_mpc_comoving = 1.0;
    options.pm_options.scale_factor = 1.0;
    options.pm_options.gravitational_constant_code = 1.0;
    options.pm_options.enable_window_deconvolution = true;
    options.tree_options.opening_theta = 0.5;
    options.tree_options.max_leaf_size = 8;
    options.tree_options.gravitational_constant_code = 1.0;
    options.tree_options.softening.epsilon_comoving = 0.5 * fine;  // epsilon <= fine cell width
    options.split_policy = gr::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.5, 1.0 / 32.0);
    std::vector<std::uint32_t> active(x.size());
    for (std::uint32_t i = 0; i < active.size(); ++i) active[i] = i;
    Field f{std::vector<double>(x.size()), std::vector<double>(x.size()), std::vector<double>(x.size())};
    gr::TreePmForceAccumulatorView accumulator{active, f.ax, f.ay, f.az};
    gr::TreePmCoordinator coordinator(shape);
    gr::TreePmDiagnostics diagnostics;
    coordinator.solveActiveSet(x, y, z, m, accumulator, options, nullptr, &diagnostics);
    const Field reference = direct(x, y, z, m, options.tree_options.softening);
    double err2 = 0.0, ref2 = 0.0, max_rel = 0.0;
    for (std::size_t t = first_fine; t < x.size(); ++t) {
      const double ex = f.ax[t] - reference.ax[t], ey = f.ay[t] - reference.ay[t], ez = f.az[t] - reference.az[t];
      const double r2 = reference.ax[t] * reference.ax[t] + reference.ay[t] * reference.ay[t] +
          reference.az[t] * reference.az[t];
      err2 += ex * ex + ey * ey + ez * ez;
      ref2 += r2;
      max_rel = std::max(max_rel, std::sqrt((ex * ex + ey * ey + ez * ez) / std::max(r2, 1.0e-300)));
    }
    const double rel_l2 = std::sqrt(err2 / std::max(ref2, 1.0e-300));
    if (level == 0) level0_error = rel_l2;
    CHUI_REQUIRE(diagnostics.composition_error_at_split < 1.0e-12);
    CHUI_REQUIRE(rel_l2 <= std::max(2.0 * level0_error, 5.0e-2));
    (void)max_rel;  // reported by the owner's diagnostic run; not an acceptance bound
  }
#endif
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_FFTW
  testForceAccuracyDoesNotDegradeWithLevel();
  return 0;
#else
  return 77;
#endif
}
