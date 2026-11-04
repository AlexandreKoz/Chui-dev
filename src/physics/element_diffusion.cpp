#include "cosmosim/physics/element_diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cosmosim::physics {
namespace {

[[nodiscard]] double harmonic(double a, double b) noexcept {
  return (a > 0.0 && b > 0.0) ? 2.0 * a * b / (a + b) : 0.0;
}

[[nodiscard]] double conductance(std::span<const double> rho_kappa, const MetalDiffusionFace& face) noexcept {
  if (face.boundary_kind == MetalDiffusionBoundaryKind::kReflective ||
      face.boundary_kind == MetalDiffusionBoundaryKind::kOpen ||
      face.left_cell >= rho_kappa.size() || face.right_cell >= rho_kappa.size() ||
      face.left_cell == face.right_cell ||
      !(face.area_code > 0.0) || !(face.center_distance_code > 0.0)) {
    return 0.0;
  }
  return harmonic(rho_kappa[face.left_cell], rho_kappa[face.right_cell]) * face.area_code /
      face.center_distance_code;
}

}  // namespace

std::uint64_t elementDiffusionWorkspaceBytes(std::size_t lanes, std::size_t rows) noexcept {
  // values (lanes*rows) + gas mass + conductance sum (2*rows) + delta
  // (lanes*rows) + lane index/before/scale/after (O(lanes)).
  constexpr std::uint64_t k_max = std::numeric_limits<std::uint64_t>::max();
  const auto mul = [](std::uint64_t a, std::uint64_t b) -> std::uint64_t {
    return (a != 0U && b > k_max / a) ? k_max : a * b;
  };
  const auto add = [](std::uint64_t a, std::uint64_t b) -> std::uint64_t {
    return a > k_max - b ? k_max : a + b;
  };
  const std::uint64_t lane_rows = mul(static_cast<std::uint64_t>(lanes), static_cast<std::uint64_t>(rows));
  std::uint64_t doubles = add(mul(lane_rows, 2U), mul(static_cast<std::uint64_t>(rows), 2U));
  doubles = add(doubles, mul(static_cast<std::uint64_t>(lanes), 3U));  // before/scale + after
  return add(mul(doubles, sizeof(double)), mul(static_cast<std::uint64_t>(lanes), sizeof(std::uint32_t)));
}

ElementDiffusionReport advanceElementResolvedDiffusion(
    core::GasCompositionSidecar& sidecar,
    std::span<const MetalDiffusionFace> faces,
    std::span<const double> rho_kappa_code,
    double dt_code,
    double parabolic_cfl,
    std::uint32_t max_subcycles,
    std::vector<double>& workspace,
    ElementDiffusionRemoteCoupling* coupling) {
  ElementDiffusionReport report;
  // Arguments that are identical on every rank (config / global dt) may
  // return before any collective; data-dependent exits may not.
  if (!sidecar.enabled() || !(dt_code > 0.0)) {
    return report;
  }
  if (coupling == nullptr && faces.empty()) {
    return report;
  }
  if (!std::isfinite(dt_code) || !(parabolic_cfl > 0.0) || parabolic_cfl > 1.0 || max_subcycles == 0U) {
    throw std::invalid_argument("element diffusion: dt must be finite, 0 < cfl <= 1 and max_subcycles >= 1");
  }
  const auto reduce = [coupling](std::exception_ptr failure, const char* stage) {
    if (coupling != nullptr) {
      coupling->reduceFailure(failure, stage);
    } else if (failure != nullptr) {
      std::rethrow_exception(failure);
    }
  };
  const core::GasCompositionLayoutDescriptor& d = sidecar.descriptor();
  const std::size_t rows = sidecar.rowCount();
  // Mixed lanes: elements, species, dust (mass lanes). SGS energy excluded.
  std::vector<std::uint32_t> lanes;
  std::size_t n_lanes = 0U;
  double* values = nullptr;
  double* gas_mass = nullptr;
  double* conductance_sum = nullptr;
  std::exception_ptr failure;
  try {
    if (rho_kappa_code.size() != rows) {
      throw std::invalid_argument("element diffusion: rho_kappa extent differs from composition rows");
    }
    for (std::uint32_t l = 0; l < d.element_lane_count; ++l) lanes.push_back(d.element_lane_begin + l);
    for (std::uint32_t l = 0; l < d.species_lane_count; ++l) lanes.push_back(d.species_lane_begin + l);
    for (std::uint32_t l = 0; l < d.dust_lane_count; ++l) lanes.push_back(d.dust_lane_begin + l);
    n_lanes = lanes.size();
    workspace.assign((2U * n_lanes + 2U) * rows, 0.0);
    values = workspace.data();
    gas_mass = workspace.data() + n_lanes * rows;
    conductance_sum = gas_mass + rows;
    for (std::size_t k = 0; k < n_lanes; ++k) {
      const auto lane = sidecar.lane(lanes[k]);
      std::copy(lane.begin(), lane.end(), values + k * rows);
    }
    for (std::size_t row = 0; row < rows; ++row) {
      double m = 0.0;
      for (std::uint32_t l = 0; l < d.element_lane_count; ++l) m += sidecar.value(d.element_lane_begin + l, row);
      for (std::uint32_t l = 0; l < d.dust_lane_count; ++l) m += sidecar.value(d.dust_lane_begin + l, row);
      gas_mass[row] = m;
    }
    for (const MetalDiffusionFace& face : faces) {
      const double c = conductance(rho_kappa_code, face);
      if (!(c > 0.0)) continue;
      conductance_sum[face.left_cell] += c;
      conductance_sum[face.right_cell] += c;
    }
  } catch (...) {
    failure = std::current_exception();
  }
  reduce(failure, "element diffusion preparation");
  report.mixed_lane_count = n_lanes;
  if (coupling != nullptr) {
    coupling->addRemoteConductance(std::span<double>(conductance_sum, rows));
  }
  double stable = std::numeric_limits<double>::infinity();
  try {
    for (std::size_t row = 0; row < rows; ++row) {
      if (conductance_sum[row] > 0.0) {
        if (!(gas_mass[row] > 0.0)) {
          throw std::runtime_error("element diffusion: a mixing cell has no composition gas mass");
        }
        stable = std::min(stable, parabolic_cfl * gas_mass[row] / conductance_sum[row]);
      }
    }
  } catch (...) {
    failure = std::current_exception();
  }
  reduce(failure, "element diffusion stability bound");
  if (coupling != nullptr) stable = coupling->globalMin(stable);
  if (!std::isfinite(stable)) {
    return report;  // no active face on any rank
  }
  report.stable_dt_code = stable;
  // Identical on every rank (global stable, global dt): a consistent refusal.
  const double needed = std::ceil(dt_code / stable);
  if (needed > static_cast<double>(max_subcycles)) {
    throw std::runtime_error(
        "element diffusion: interval needs " + std::to_string(needed) + " subcycles > max_subcycles " +
        std::to_string(max_subcycles) + " (fail closed; reject and retry with a shorter interval)");
  }
  const std::uint32_t subcycles = std::max<std::uint32_t>(1U, static_cast<std::uint32_t>(needed));
  const double dt_sub = dt_code / static_cast<double>(subcycles);
  // before/scale are global lane totals (for the global conservation check).
  std::vector<double> before_scale(2U * n_lanes, 0.0);
  for (std::size_t k = 0; k < n_lanes; ++k) {
    for (std::size_t row = 0; row < rows; ++row) {
      before_scale[k] += values[k * rows + row];
      before_scale[n_lanes + k] += std::abs(values[k * rows + row]);
    }
  }
  if (coupling != nullptr) coupling->globalSum(before_scale);
  const std::span<double> delta(values + n_lanes * rows + 2U * rows, n_lanes * rows);
  for (std::uint32_t sub = 0; sub < subcycles; ++sub) {
    std::fill(delta.begin(), delta.end(), 0.0);
    for (const MetalDiffusionFace& face : faces) {
      const double c = conductance(rho_kappa_code, face);
      if (!(c > 0.0)) continue;
      const std::size_t a = face.left_cell;
      const std::size_t b = face.right_cell;
      for (std::size_t k = 0; k < n_lanes; ++k) {
        const double x_a = values[k * rows + a] / gas_mass[a];
        const double x_b = values[k * rows + b] / gas_mass[b];
        const double transfer = dt_sub * c * (x_a - x_b);
        delta[k * rows + a] -= transfer;
        delta[k * rows + b] += transfer;
      }
      ++report.faces_evaluated;
    }
    if (coupling != nullptr) {
      // Fresh remote values of THIS subcycle (values/gas_mass before update).
      coupling->addRemoteFaceTransfers(std::span<const double>(values, n_lanes * rows),
                                       std::span<const double>(gas_mass, rows), n_lanes, dt_sub, delta);
    }
    for (std::size_t i = 0; i < delta.size(); ++i) {
      values[i] += delta[i];
    }
    // Gas mass is invariant in exact arithmetic (sum x = 1); recompute from
    // the evolved element + dust lanes so roundoff never accumulates.
    for (std::size_t row = 0; row < rows; ++row) {
      double m = 0.0;
      for (std::size_t k = 0; k < d.element_lane_count; ++k) m += values[k * rows + row];
      for (std::size_t k = 0; k < d.dust_lane_count; ++k) {
        m += values[(d.element_lane_count + d.species_lane_count + k) * rows + row];
      }
      if (conductance_sum[row] > 0.0) {
        report.max_relative_cell_mass_change =
            std::max(report.max_relative_cell_mass_change, std::abs(m - gas_mass[row]) / gas_mass[row]);
      }
      gas_mass[row] = m;
    }
  }
  std::vector<double> after(n_lanes, 0.0);
  try {
    for (std::size_t k = 0; k < n_lanes; ++k) {
      for (std::size_t row = 0; row < rows; ++row) {
        const double v = values[k * rows + row];
        if (!std::isfinite(v) || v < -1.0e-14 * std::max(before_scale[n_lanes + k], 1.0e-300)) {
          throw std::runtime_error("element diffusion: lane became negative (positivity bound violated)");
        }
        after[k] += v;
      }
    }
  } catch (...) {
    failure = std::current_exception();
  }
  // Every rank commits or every rank throws (no partial distributed commit).
  reduce(failure, "element diffusion positivity");
  if (coupling != nullptr) coupling->globalSum(after);
  for (std::size_t k = 0; k < n_lanes; ++k) {
    report.max_relative_lane_residual = std::max(
        report.max_relative_lane_residual,
        std::abs(after[k] - before_scale[k]) / std::max(before_scale[n_lanes + k], 1.0e-300));
  }
  // Commit (no exception past this point).
  for (std::size_t k = 0; k < n_lanes; ++k) {
    auto lane = sidecar.lane(lanes[k]);
    for (std::size_t row = 0; row < rows; ++row) {
      lane[row] = std::max(values[k * rows + row], 0.0);
    }
  }
  report.subcycles = subcycles;
  return report;
}

}  // namespace cosmosim::physics
