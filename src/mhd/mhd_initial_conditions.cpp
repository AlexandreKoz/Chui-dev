#include "cosmosim/mhd/mhd_initial_conditions.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "cosmosim/mhd/mhd_induction.hpp"

namespace cosmosim::mhd {
namespace {

// SplitMix64: deterministic, platform-independent integer stream.
struct SplitMix64 {
  std::uint64_t state = 0U;
  std::uint64_t next() noexcept {
    std::uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27U)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31U);
  }
  double uniform() noexcept {  // (0, 1]
    return (static_cast<double>(next() >> 11U) + 1.0) * (1.0 / 9007199254740992.0);
  }
  double gaussian() noexcept {  // Box-Muller, one variate per call
    const double u1 = uniform();
    const double u2 = uniform();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * std::numbers::pi * u2);
  }
};

struct FourierMode {
  std::array<double, 3> k{};
  std::array<double, 3> cos_amplitude{};
  std::array<double, 3> sin_amplitude{};
};

[[nodiscard]] std::array<double, 3> projectPerpendicular(std::array<double, 3> v,
                                                         const std::array<double, 3>& k) {
  const double k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
  const double dot = v[0] * k[0] + v[1] * k[1] + v[2] * k[2];
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    v[axis] -= dot * k[axis] / k2;
  }
  return v;
}

}  // namespace

double integrateEdgePotential(
    const std::function<std::array<double, 3>(std::array<double, 3>)>& vector_potential,
    MhdAxis axis,
    std::array<double, 3> start_comoving,
    double length) {
  // Three-point Gauss-Legendre on [0, length] (exact for quintics).
  constexpr std::array<double, 3> nodes{-0.7745966692414834, 0.0, 0.7745966692414834};
  constexpr std::array<double, 3> weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  const std::size_t direction = mhdAxisIndex(axis);
  double sum = 0.0;
  for (std::size_t q = 0; q < 3U; ++q) {
    std::array<double, 3> point = start_comoving;
    point[direction] += 0.5 * length * (1.0 + nodes[q]);
    sum += weights[q] * vector_potential(point)[direction];
  }
  return 0.5 * length * sum;
}

double seedFaceFlux(const MhdPatchGeometry& geometry,
                    std::array<double, 3> uniform_b_comoving,
                    const MhdEdgePotential& potential,
                    std::span<double> phi_x,
                    std::span<double> phi_y,
                    std::span<double> phi_z) {
  if (!geometry.valid() || phi_x.size() != geometry.faceCountX() ||
      phi_y.size() != geometry.faceCountY() || phi_z.size() != geometry.faceCountZ()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_initial_conditions",
                                     .detail = "seed target does not match the patch"});
  }
  for (const double value : uniform_b_comoving) {
    if (!std::isfinite(value)) {
      throw MhdError(MhdErrorCode::kInvalidField,
                     MhdFailureContext{.owner = "mhd_initial_conditions",
                                       .detail = "non-finite uniform seed field"});
    }
  }
  std::fill(phi_x.begin(), phi_x.end(), uniform_b_comoving[0] * geometry.faceAreaX());
  std::fill(phi_y.begin(), phi_y.end(), uniform_b_comoving[1] * geometry.faceAreaY());
  std::fill(phi_z.begin(), phi_z.end(), uniform_b_comoving[2] * geometry.faceAreaZ());
  if (potential) {
    MhdEdgeEmfStorage edges(geometry);
    const std::array<double, 3> origin{geometry.origin_x_comoving, geometry.origin_y_comoving,
                                       geometry.origin_z_comoving};
    const std::array<double, 3> spacing{geometry.dx_comoving, geometry.dy_comoving,
                                        geometry.dz_comoving};
    const auto fill = [&](MhdAxis axis, std::span<double> lane, std::size_t ni, std::size_t nj,
                          std::size_t nk, auto index) {
      for (std::size_t k = 0; k < nk; ++k) {
        for (std::size_t j = 0; j < nj; ++j) {
          for (std::size_t i = 0; i < ni; ++i) {
            const std::array<double, 3> start{origin[0] + static_cast<double>(i) * spacing[0],
                                              origin[1] + static_cast<double>(j) * spacing[1],
                                              origin[2] + static_cast<double>(k) * spacing[2]};
            const double value = potential(axis, start, spacing[mhdAxisIndex(axis)]);
            if (!std::isfinite(value)) {
              throw MhdError(MhdErrorCode::kInvalidField,
                             MhdFailureContext{.owner = "mhd_initial_conditions",
                                               .detail = "non-finite edge potential"});
            }
            lane[index(i, j, k)] = value;
          }
        }
      }
    };
    fill(MhdAxis::kX, edges.emfX(), geometry.nx, geometry.ny + 1U, geometry.nz + 1U,
         [&](std::size_t i, std::size_t j, std::size_t k) { return geometry.edgeIndexX(i, j, k); });
    fill(MhdAxis::kY, edges.emfY(), geometry.nx + 1U, geometry.ny, geometry.nz + 1U,
         [&](std::size_t i, std::size_t j, std::size_t k) { return geometry.edgeIndexY(i, j, k); });
    fill(MhdAxis::kZ, edges.emfZ(), geometry.nx + 1U, geometry.ny + 1U, geometry.nz,
         [&](std::size_t i, std::size_t j, std::size_t k) { return geometry.edgeIndexZ(i, j, k); });
    // Stokes: Phi_f = sum over boundary edges of the oriented A line integral,
    // i.e. Phi += C(A) == applyConstrainedTransportUpdate with weight -1.
    applyConstrainedTransportUpdate(geometry, phi_x, phi_y, phi_z, edges.emfX(), edges.emfY(),
                                    edges.emfZ(), -1.0);
  }
  const double normalizer = magneticDivergenceNormalizer(geometry, phi_x, phi_y, phi_z);
  return maxNormalizedFaceDivergence(geometry, phi_x, phi_y, phi_z, normalizer);
}

MhdEdgePotential fieldLoopPotential(double amplitude, double radius, std::array<double, 3> center) {
  if (!(radius > 0.0) || !std::isfinite(amplitude)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_initial_conditions",
                                     .detail = "invalid field-loop parameters"});
  }
  return [amplitude, radius, center](MhdAxis axis, std::array<double, 3> start, double length) {
    if (axis != MhdAxis::kZ) {
      return 0.0;  // A = (0, 0, A_z(x, y))
    }
    const double dx = start[0] - center[0];
    const double dy = start[1] - center[1];
    const double r = std::sqrt(dx * dx + dy * dy);
    return r < radius ? amplitude * (radius - r) * length : 0.0;
  };
}

MhdEdgePotential circularAlfvenPotential(double amplitude, double wavelength) {
  if (!(wavelength > 0.0) || !std::isfinite(amplitude)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_initial_conditions",
                                     .detail = "invalid Alfven-wave parameters"});
  }
  const double k = 2.0 * std::numbers::pi / wavelength;
  return [amplitude, k](MhdAxis axis, std::array<double, 3> start, double length) {
    // A_y, A_z depend on x only, so y- and z-edge integrals are exact.
    if (axis == MhdAxis::kY) {
      return amplitude * std::sin(k * start[0]) / k * length;
    }
    if (axis == MhdAxis::kZ) {
      return amplitude * std::cos(k * start[0]) / k * length;
    }
    return 0.0;
  };
}

MhdEdgePotential stochasticSeedPotential(std::array<double, 3> box_comoving, std::int32_t k_max,
                                         double amplitude, double spectral_index,
                                         std::uint64_t seed) {
  if (k_max < 1 || k_max > 64 || !(box_comoving[0] > 0.0) || !(box_comoving[1] > 0.0) ||
      !(box_comoving[2] > 0.0) || !std::isfinite(amplitude) || !std::isfinite(spectral_index)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_initial_conditions",
                                     .detail = "invalid stochastic seed parameters"});
  }
  SplitMix64 random{seed};
  std::vector<FourierMode> modes;
  // Half k-space (Hermitian symmetry: the real cos/sin pair represents both
  // +k and -k), deterministic lexicographic order.
  for (std::int32_t kz = 0; kz <= k_max; ++kz) {
    for (std::int32_t ky = -k_max; ky <= k_max; ++ky) {
      for (std::int32_t kx = -k_max; kx <= k_max; ++kx) {
        const bool positive_half = kz > 0 || (kz == 0 && (ky > 0 || (ky == 0 && kx > 0)));
        const double magnitude = std::sqrt(static_cast<double>(kx * kx + ky * ky + kz * kz));
        if (!positive_half || magnitude > static_cast<double>(k_max)) {
          continue;
        }
        FourierMode mode{};
        mode.k = {2.0 * std::numbers::pi * kx / box_comoving[0],
                  2.0 * std::numbers::pi * ky / box_comoving[1],
                  2.0 * std::numbers::pi * kz / box_comoving[2]};
        const double scale = amplitude * std::pow(magnitude, spectral_index);
        std::array<double, 3> c{random.gaussian(), random.gaussian(), random.gaussian()};
        std::array<double, 3> s{random.gaussian(), random.gaussian(), random.gaussian()};
        c = projectPerpendicular(c, mode.k);
        s = projectPerpendicular(s, mode.k);
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          mode.cos_amplitude[axis] = scale * c[axis];
          mode.sin_amplitude[axis] = scale * s[axis];
        }
        modes.push_back(mode);
      }
    }
  }
  const auto potential = [modes](std::array<double, 3> x) {
    std::array<double, 3> a{};
    for (const FourierMode& mode : modes) {
      const double phase = mode.k[0] * x[0] + mode.k[1] * x[1] + mode.k[2] * x[2];
      const double cosine = std::cos(phase);
      const double sine = std::sin(phase);
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        a[axis] += mode.cos_amplitude[axis] * cosine + mode.sin_amplitude[axis] * sine;
      }
    }
    return a;
  };
  return [potential](MhdAxis axis, std::array<double, 3> start, double length) {
    return integrateEdgePotential(potential, axis, start, length);
  };
}

}  // namespace cosmosim::mhd
