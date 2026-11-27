#pragma once

// Template implementation of sampleParallelAtomVelocity (included by
// resonant_lyman_alpha.hpp). Derivation in that header and in the numerical
// contracts (L02).

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cosmosim::observables {

template <typename Uniform>
double sampleParallelAtomVelocity(double x_in, double a, Uniform&& next) {
  constexpr double pi = 3.141592653589793238462643383279502884;
  const double sign = x_in < 0.0 ? -1.0 : 1.0;
  const double x = std::abs(x_in);
  const double d = std::max(1.0, 0.5 * x);
  const double u_s = x - d;
  const double bound2 = u_s > 0.0 ? std::exp(-u_s * u_s) : 1.0;
  // Envelope masses: W1 = int_{-inf}^{u_s} exp(-u^2)/d^2, W2 = bound2 int_{u_s}^{inf} du/((x-u)^2+a^2).
  const double w1 = 0.5 * std::sqrt(pi) * (1.0 + std::erf(u_s)) / (d * d);
  const double w2 = bound2 / a * (0.5 * pi + std::atan(d / a));
  for (int attempt = 0; attempt < 100000000; ++attempt) {
    if (next() * (w1 + w2) < w1) {
      // Gaussian (variance 1/2) truncated to u < u_s by resampling (Box-Muller).
      double u = 0.0;
      do {
        const double r = std::sqrt(-std::log(next()));
        u = r * std::cos(2.0 * pi * next());
      } while (u >= u_s);
      const double accept = d * d / ((x - u) * (x - u) + a * a);
      if (next() < accept) return sign * u;
    } else {
      // Cauchy centred at x, width a, truncated to u >= u_s: phi uniform in [-atan(d/a), pi/2).
      const double lo = -std::atan(d / a);
      const double phi = lo + (0.5 * pi - lo) * next();
      const double u = x + a * std::tan(phi);
      if (next() < std::exp(-u * u) / bound2) return sign * u;
    }
  }
  throw std::runtime_error("parallel atom velocity sampler exceeded 1e8 rejections (a must be > 0)");
}

}  // namespace cosmosim::observables
