#pragma once

// Shared deterministic fixtures for the CHUI-MHD-REF-v1 test sources.
// NOT EXECUTED in the source-only campaign; see
// docs/architecture/chui_mhd_implementation_pass01_status.md.

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <numbers>
#include <vector>

#include "../../support/test_require.hpp"
#include "cosmosim/mhd/mhd_initial_conditions.hpp"
#include "cosmosim/mhd/mhd_level_stepper.hpp"
#include "cosmosim/mhd/mhd_step_operator.hpp"

namespace cosmosim::mhd::test {

struct GasPrimitive {
  double rho = 1.0;
  std::array<double, 3> velocity{};
  double pressure = 1.0;
  double metal_fraction = 0.0;
};

using GasProfile = std::function<GasPrimitive(std::array<double, 3> cell_center)>;

// A same-level tiling of a rectangular domain by uniform patches. Storage is
// reserved up front so the bindings (raw views) stay valid.
class TestLevel {
 public:
  TestLevel(std::array<std::int64_t, 3> domain_cells, std::array<std::int64_t, 3> patch_cells,
            std::array<double, 3> spacing, MhdLevelDomain boundaries)
      : m_domain(boundaries), m_patch_cells(patch_cells), m_spacing(spacing) {
    m_domain.cells = domain_cells;
    std::array<std::int64_t, 3> grid{};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      grid[axis] = domain_cells[axis] / patch_cells[axis];
    }
    const std::size_t count = static_cast<std::size_t>(grid[0] * grid[1] * grid[2]);
    m_gas.reserve(count);
    m_faces.reserve(count);
    std::uint64_t id = 100U;
    for (std::int64_t gz = 0; gz < grid[2]; ++gz) {
      for (std::int64_t gy = 0; gy < grid[1]; ++gy) {
        for (std::int64_t gx = 0; gx < grid[0]; ++gx) {
          MhdPatchGeometry geometry{};
          geometry.nx = static_cast<std::size_t>(patch_cells[0]);
          geometry.ny = static_cast<std::size_t>(patch_cells[1]);
          geometry.nz = static_cast<std::size_t>(patch_cells[2]);
          geometry.dx_comoving = spacing[0];
          geometry.dy_comoving = spacing[1];
          geometry.dz_comoving = spacing[2];
          geometry.origin_x_comoving = static_cast<double>(gx * patch_cells[0]) * spacing[0];
          geometry.origin_y_comoving = static_cast<double>(gy * patch_cells[1]) * spacing[1];
          geometry.origin_z_comoving = static_cast<double>(gz * patch_cells[2]) * spacing[2];
          m_gas.emplace_back(geometry.cellCount());
          m_faces.emplace_back(geometry);
          MhdLevelPatchBinding binding{};
          // Deliberately non-monotone ids so owner selection is exercised.
          binding.patch_id = id;
          id = id * 7U % 1009U + 3U;
          binding.geometry_generation = 1U;
          binding.level = 0;
          binding.origin_cell = {gx * patch_cells[0], gy * patch_cells[1], gz * patch_cells[2]};
          binding.geometry = geometry;
          m_bindings.push_back(binding);
        }
      }
    }
    for (std::size_t index = 0; index < m_bindings.size(); ++index) {
      m_bindings[index].gas = &m_gas[index];
      m_bindings[index].faces = makeMutableFaceFluxView(m_faces[index], m_bindings[index].patch_id,
                                                        1U, 0);
    }
  }

  // Seed faces (uniform + curl A, evaluated at global positions) and gas with
  // the MHD total energy E_c = p/(g-1) + rho u^2/2 + B_c^2/(2a).
  void initialize(const GasProfile& gas, std::array<double, 3> uniform_b,
                  const MhdEdgePotential& potential, double scale_factor, double gamma) {
    for (std::size_t index = 0; index < m_bindings.size(); ++index) {
      MhdFaceFluxStorage& faces = m_faces[index];
      const double divergence = seedFaceFlux(faces.geometry(), uniform_b, potential,
                                             faces.fluxX(), faces.fluxY(), faces.fluxZ());
      CHUI_REQUIRE(divergence < 1.0e-12);
      const MhdPatchGeometry& geometry = faces.geometry();
      for (std::size_t k = 0; k < geometry.nz; ++k) {
        for (std::size_t j = 0; j < geometry.ny; ++j) {
          for (std::size_t i = 0; i < geometry.nx; ++i) {
            const std::array<double, 3> center{
                geometry.origin_x_comoving + (static_cast<double>(i) + 0.5) * geometry.dx_comoving,
                geometry.origin_y_comoving + (static_cast<double>(j) + 0.5) * geometry.dy_comoving,
                geometry.origin_z_comoving + (static_cast<double>(k) + 0.5) * geometry.dz_comoving};
            const GasPrimitive w = gas(center);
            const std::array<double, 3> b = faces.cellFieldComoving(i, j, k);
            MhdConservedState conserved{};
            conserved.rho_comoving = w.rho;
            conserved.momentum_density_x_comoving = w.rho * w.velocity[0];
            conserved.momentum_density_y_comoving = w.rho * w.velocity[1];
            conserved.momentum_density_z_comoving = w.rho * w.velocity[2];
            conserved.total_energy_density_comoving = w.pressure / (gamma - 1.0) +
                0.5 * w.rho *
                    (w.velocity[0] * w.velocity[0] + w.velocity[1] * w.velocity[1] +
                     w.velocity[2] * w.velocity[2]) +
                0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / scale_factor;
            conserved.metal_mass_density_comoving = w.rho * w.metal_fraction;
            m_gas[index].storeGas(geometry.cellIndex(i, j, k), conserved);
          }
        }
      }
    }
  }

  [[nodiscard]] MhdLevelStepper makeStepper(const MhdNumericsPolicy& policy) const {
    return MhdLevelStepper(policy, m_domain, m_bindings);
  }

  [[nodiscard]] std::size_t patchCount() const noexcept { return m_bindings.size(); }
  [[nodiscard]] MhdGasStateSoa& gas(std::size_t index) { return m_gas[index]; }
  [[nodiscard]] MhdFaceFluxStorage& faces(std::size_t index) { return m_faces[index]; }
  [[nodiscard]] const MhdLevelPatchBinding& binding(std::size_t index) const {
    return m_bindings[index];
  }

  // Global cell lookup (domain coordinates) of a conserved lane value.
  [[nodiscard]] double cellValue(std::size_t lane, std::array<std::int64_t, 3> cell) {
    for (std::size_t index = 0; index < m_bindings.size(); ++index) {
      const auto& origin = m_bindings[index].origin_cell;
      bool inside = true;
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        inside = inside && cell[axis] >= origin[axis] &&
            cell[axis] < origin[axis] + m_patch_cells[axis];
      }
      if (inside) {
        const MhdPatchGeometry& geometry = m_bindings[index].geometry;
        return m_gas[index].lane(lane)[geometry.cellIndex(
            static_cast<std::size_t>(cell[0] - origin[0]),
            static_cast<std::size_t>(cell[1] - origin[1]),
            static_cast<std::size_t>(cell[2] - origin[2]))];
      }
    }
    CHUI_REQUIRE(false);
    return 0.0;
  }

 private:
  MhdLevelDomain m_domain;
  std::array<std::int64_t, 3> m_patch_cells;
  std::array<double, 3> m_spacing;
  std::vector<MhdGasStateSoa> m_gas;
  std::vector<MhdFaceFluxStorage> m_faces;
  std::vector<MhdLevelPatchBinding> m_bindings;
};

[[nodiscard]] inline MhdLevelDomain periodicDomain() { return MhdLevelDomain{}; }

[[nodiscard]] inline MhdLevelDomain uniformDomain(MhdBoundaryKind kind) {
  MhdLevelDomain domain{};
  domain.lower = {kind, kind, kind};
  domain.upper = {kind, kind, kind};
  return domain;
}

[[nodiscard]] inline MhdStageContext isolatedInterval(double time, double dt) {
  MhdStageContext context{};
  context.time_code = time;
  context.dt_code = dt;
  return context;
}

[[nodiscard]] inline double relativeDifference(double lhs, double rhs, double scale) {
  return std::abs(lhs - rhs) / std::max(scale, 1.0e-300);
}

}  // namespace cosmosim::mhd::test
