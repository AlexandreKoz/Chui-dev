#include "cosmosim/physics/agn/agn_magnetic.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cosmosim/physics/agn/agn_injection.hpp"

namespace cosmosim::physics {
namespace {

[[nodiscard]] double profile(const std::array<double, 3>& p, const std::array<double, 3>& c, double h) noexcept {
  const double dx = p[0] - c[0];
  const double dy = p[1] - c[1];
  const double dz = p[2] - c[2];
  const double q2 = (dx * dx + dy * dy + dz * dz) / (h * h);
  if (q2 >= 1.0) return 0.0;
  const double t = 1.0 - q2;
  return t * t * t;
}

}  // namespace

MagneticFaceFluxes buildToroidalMagneticInjection(const MagneticPatchGeometry& g, const std::array<double, 3>& c,
                                                  const std::array<double, 3>& axis_in, double h) {
  if (g.nx == 0U || g.ny == 0U || g.nz == 0U || !(h > 0.0)) {
    throw std::invalid_argument("buildToroidalMagneticInjection: invalid geometry or radius");
  }
  const double an = std::sqrt(axis_in[0] * axis_in[0] + axis_in[1] * axis_in[1] + axis_in[2] * axis_in[2]);
  if (!(an > 0.0)) throw std::invalid_argument("buildToroidalMagneticInjection: undefined axis");
  const std::array<double, 3> n{axis_in[0] / an, axis_in[1] / an, axis_in[2] / an};
  const std::size_t nx = g.nx, ny = g.ny, nz = g.nz;
  const double dx = g.width[0], dy = g.width[1], dz = g.width[2];
  // Edge potentials (edge midpoints).
  std::vector<double> ax(nx * (ny + 1) * (nz + 1)), ay((nx + 1) * ny * (nz + 1)), az((nx + 1) * (ny + 1) * nz);
  for (std::size_t k = 0; k <= nz; ++k)
    for (std::size_t j = 0; j <= ny; ++j)
      for (std::size_t i = 0; i < nx; ++i)
        ax[i + nx * (j + (ny + 1) * k)] =
            n[0] * profile({g.origin[0] + (i + 0.5) * dx, g.origin[1] + j * dy, g.origin[2] + k * dz}, c, h);
  for (std::size_t k = 0; k <= nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i <= nx; ++i)
        ay[i + (nx + 1) * (j + ny * k)] =
            n[1] * profile({g.origin[0] + i * dx, g.origin[1] + (j + 0.5) * dy, g.origin[2] + k * dz}, c, h);
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j <= ny; ++j)
      for (std::size_t i = 0; i <= nx; ++i)
        az[i + (nx + 1) * (j + (ny + 1) * k)] =
            n[2] * profile({g.origin[0] + i * dx, g.origin[1] + j * dy, g.origin[2] + (k + 0.5) * dz}, c, h);
  auto AX = [&](std::size_t i, std::size_t j, std::size_t k) { return ax[i + nx * (j + (ny + 1) * k)]; };
  auto AY = [&](std::size_t i, std::size_t j, std::size_t k) { return ay[i + (nx + 1) * (j + ny * k)]; };
  auto AZ = [&](std::size_t i, std::size_t j, std::size_t k) { return az[i + (nx + 1) * (j + (ny + 1) * k)]; };
  MagneticFaceFluxes f;
  f.x.assign((nx + 1) * ny * nz, 0.0);
  f.y.assign(nx * (ny + 1) * nz, 0.0);
  f.z.assign(nx * ny * (nz + 1), 0.0);
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i <= nx; ++i)
        f.x[i + (nx + 1) * (j + ny * k)] = dz * (AZ(i, j + 1, k) - AZ(i, j, k)) - dy * (AY(i, j, k + 1) - AY(i, j, k));
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j <= ny; ++j)
      for (std::size_t i = 0; i < nx; ++i)
        f.y[i + nx * (j + (ny + 1) * k)] = dx * (AX(i, j, k + 1) - AX(i, j, k)) - dz * (AZ(i + 1, j, k) - AZ(i, j, k));
  for (std::size_t k = 0; k <= nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i)
        f.z[i + nx * (j + ny * k)] = dy * (AY(i + 1, j, k) - AY(i, j, k)) - dx * (AX(i, j + 1, k) - AX(i, j, k));
  return f;
}

double magneticDiscreteDivergenceMax(const MagneticPatchGeometry& g, const MagneticFaceFluxes& f) {
  const std::size_t nx = g.nx, ny = g.ny, nz = g.nz;
  double worst = 0.0;
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i) {
        const double d = f.x[(i + 1) + (nx + 1) * (j + ny * k)] - f.x[i + (nx + 1) * (j + ny * k)] +
            f.y[i + nx * ((j + 1) + (ny + 1) * k)] - f.y[i + nx * (j + (ny + 1) * k)] +
            f.z[i + nx * (j + ny * (k + 1))] - f.z[i + nx * (j + ny * k)];
        worst = std::max(worst, std::abs(d));
      }
  return worst;
}

MagneticEnergyQuadratic magneticInjectionEnergyTerms(const MagneticPatchGeometry& g, const MagneticFaceFluxes& e,
                                                     const MagneticFaceFluxes& b) {
  const std::size_t nx = g.nx, ny = g.ny, nz = g.nz;
  const double ax = g.width[1] * g.width[2], ay = g.width[0] * g.width[2], az = g.width[0] * g.width[1];
  const double vol = g.width[0] * g.width[1] * g.width[2];
  MagneticEnergyQuadratic q;
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i) {
        const std::size_t x0 = i + (nx + 1) * (j + ny * k), x1 = x0 + 1;
        const std::size_t y0 = i + nx * (j + (ny + 1) * k), y1 = i + nx * ((j + 1) + (ny + 1) * k);
        const std::size_t z0 = i + nx * (j + ny * k), z1 = i + nx * (j + ny * (k + 1));
        const double bex = 0.5 * (e.x[x0] + e.x[x1]) / ax, bey = 0.5 * (e.y[y0] + e.y[y1]) / ay,
                     bez = 0.5 * (e.z[z0] + e.z[z1]) / az;
        const double bix = 0.5 * (b.x[x0] + b.x[x1]) / ax, biy = 0.5 * (b.y[y0] + b.y[y1]) / ay,
                     biz = 0.5 * (b.z[z0] + b.z[z1]) / az;
        q.linear += vol * (bex * bix + bey * biy + bez * biz);
        q.quadratic += 0.5 * vol * (bix * bix + biy * biy + biz * biz);
      }
  return q;
}

bool solveMagneticInjectionAmplitude(const MagneticEnergyQuadratic& t, double energy, double& amplitude) noexcept {
  return solveEnergyConsistentAmplitude(t.quadratic, t.linear, energy, amplitude);
}

}  // namespace cosmosim::physics
