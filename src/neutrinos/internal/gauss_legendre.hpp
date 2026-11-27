#pragma once

// Gauss-Legendre rules computed by Newton iteration on P_n (no transcribed
// node tables). Shared by the neutrino kernels.

#include <cmath>
#include <cstddef>
#include <vector>

namespace cosmosim::neutrinos::internal {

struct GaussLegendreRule {
  std::vector<double> nodes;    // on [-1, 1]
  std::vector<double> weights;  // sum to 2
};

[[nodiscard]] inline GaussLegendreRule makeGaussLegendreRule(std::size_t n) {
  GaussLegendreRule rule;
  rule.nodes.resize(n);
  rule.weights.resize(n);
  const double pi = 3.141592653589793238462643383279502884;
  for (std::size_t i = 0; i < n; ++i) {
    double x = std::cos(pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
    double dp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0;
      double p1 = x;
      for (std::size_t k = 2; k <= n; ++k) {
        const double p2 = ((2.0 * static_cast<double>(k) - 1.0) * x * p1 - (static_cast<double>(k) - 1.0) * p0) /
                          static_cast<double>(k);
        p0 = p1;
        p1 = p2;
      }
      dp = static_cast<double>(n) * (x * p1 - p0) / (x * x - 1.0);
      const double step = p1 / dp;
      x -= step;
      if (std::abs(step) < 1.0e-16) break;
    }
    rule.nodes[i] = x;
    rule.weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
  }
  return rule;
}

[[nodiscard]] inline const GaussLegendreRule& gaussLegendre16() {
  static const GaussLegendreRule rule = makeGaussLegendreRule(16U);
  return rule;
}

}  // namespace cosmosim::neutrinos::internal
