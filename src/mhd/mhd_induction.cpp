#include "cosmosim/mhd/mhd_induction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace cosmosim::mhd {
namespace {

double inverseFloor(const MhdFaceSignalPatch& signals) noexcept {
  double floor = MhdEmfBuilder::k_degenerate_weight_floor;
  const auto consider = [&floor](std::span<const double> speeds) {
    for (const double speed : speeds) {
      if (std::isfinite(speed)) {
        floor = std::max(floor, std::abs(speed) * 1.0e-15);
      }
    }
  };
  consider(signals.face_speed_x);
  consider(signals.face_speed_y);
  consider(signals.face_speed_z);
  return floor;
}

// HLL-qualified weighted edge value from four face estimates. A face whose
// signal fan is faster contributes proportionally less. When every weight is
// at the floor (a genuinely static fan) the arithmetic mean is used, but only
// if the four estimates agree within the declared relative tolerance.
double blendEdge(const std::array<double, 4>& estimates,
                 const std::array<double, 4>& speeds,
                 double floor,
                 std::uint64_t* degenerate_count,
                 const char* label) {
  double weight_sum = 0.0;
  double weighted = 0.0;
  for (std::size_t index = 0; index < estimates.size(); ++index) {
    const double weight = 1.0 / (std::abs(speeds[index]) + floor);
    weight_sum += weight;
    weighted += weight * estimates[index];
  }
  if (weight_sum > 0.0 && std::isfinite(weight_sum)) {
    // Degenerate only when every signal speed is exactly zero.
    bool all_at_floor = true;
    for (const double speed : speeds) {
      if (std::abs(speed) > 0.0) {
        all_at_floor = false;
        break;
      }
    }
    if (!all_at_floor) {
      return weighted / weight_sum;
    }
    if (degenerate_count != nullptr) {
      ++(*degenerate_count);
    }
  }
  double minimum_estimate = estimates[0];
  double maximum_estimate = estimates[0];
  double mean = 0.0;
  for (const double estimate : estimates) {
    minimum_estimate = std::min(minimum_estimate, estimate);
    maximum_estimate = std::max(maximum_estimate, estimate);
    mean += estimate;
  }
  mean *= 0.25;
  const double spread = maximum_estimate - minimum_estimate;
  const double scale = std::max(std::abs(mean), 1.0);
  if (spread > MhdEmfBuilder::k_disagreement_relative_tolerance * scale) {
    std::ostringstream stream;
    stream << "degenerate edge EMF blend for " << label
           << " found disagreeing face estimates (spread=" << spread
           << ", mean=" << mean << ')';
    MhdFailureContext context{};
    context.owner = "mhd_emf_builder";
    context.detail = stream.str();
    throw MhdError(MhdErrorCode::kInvalidField, std::move(context));
  }
  return mean;
}

}  // namespace

void MhdFaceSignalPatch::resize(const MhdPatchGeometry& patch_geometry) {
  geometry = patch_geometry;
  const std::size_t count_x = patch_geometry.faceCountX();
  const std::size_t count_y = patch_geometry.faceCountY();
  const std::size_t count_z = patch_geometry.faceCountZ();
  mass_flux_x.assign(count_x, 0.0);
  electric_t1_x.assign(count_x, 0.0);
  electric_t2_x.assign(count_x, 0.0);
  face_speed_x.assign(count_x, 0.0);
  mass_flux_y.assign(count_y, 0.0);
  electric_t1_y.assign(count_y, 0.0);
  electric_t2_y.assign(count_y, 0.0);
  face_speed_y.assign(count_y, 0.0);
  mass_flux_z.assign(count_z, 0.0);
  electric_t1_z.assign(count_z, 0.0);
  electric_t2_z.assign(count_z, 0.0);
  face_speed_z.assign(count_z, 0.0);
}

std::array<double, 3> MhdEmfBuilder::faceElectricFieldFromFlux(
    const MhdRiemannFlux& flux) noexcept {
  return {flux.electric_field[0], flux.electric_field[1], flux.electric_field[2]};
}

void MhdEmfBuilder::buildUctHll(const MhdPatchGeometry& geometry,
                                const MhdFaceSignalPatch& signals,
                                MhdEdgeEmfStorage& emf,
                                std::uint64_t* degenerate_edge_count) {
  emf.resize(geometry);
  const double floor = inverseFloor(signals);
  const std::size_t nx = geometry.nx;
  const std::size_t ny = geometry.ny;
  const std::size_t nz = geometry.nz;

  // x-directed edges: (nx, ny+1, nz+1)
  for (std::size_t k = 0; k <= nz; ++k) {
    for (std::size_t j = 0; j <= ny; ++j) {
      for (std::size_t i = 0; i < nx; ++i) {
        // Adjacent y-faces: (i,j,k) and (i,j,k-1) in z.
        // Adjacent z-faces: (i,j,k) and (i,j-1,k) in y.
        const std::size_t y_face_high = geometry.faceIndexY(i, j, k);
        const std::size_t y_face_low = k > 0 ? geometry.faceIndexY(i, j, k - 1U) : y_face_high;
        const std::size_t z_face_high = geometry.faceIndexZ(i, j, k);
        const std::size_t z_face_low = j > 0 ? geometry.faceIndexZ(i, j - 1U, k) : z_face_high;
        // E_x on a y-face is the second tangential component; on a z-face it is
        // the first tangential component.
        const std::array<double, 4> estimates{
            signals.electric_t2_y[y_face_high],
            signals.electric_t2_y[y_face_low],
            signals.electric_t1_z[z_face_high],
            signals.electric_t1_z[z_face_low],
        };
        const std::array<double, 4> speeds{
            signals.face_speed_y[y_face_high],
            signals.face_speed_y[y_face_low],
            signals.face_speed_z[z_face_high],
            signals.face_speed_z[z_face_low],
        };
        const double field_component =
            blendEdge(estimates, speeds, floor, degenerate_edge_count, "E_x");
        // Edge-integrated EMF = field * edge length.
        emf.emfX()[geometry.edgeIndexX(i, j, k)] = field_component * geometry.dx_comoving;
      }
    }
  }

  // y-directed edges: (nx+1, ny, nz+1)
  for (std::size_t k = 0; k <= nz; ++k) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t i = 0; i <= nx; ++i) {
        const std::size_t x_face_high = geometry.faceIndexX(i, j, k);
        const std::size_t x_face_low = k > 0 ? geometry.faceIndexX(i, j, k - 1U) : x_face_high;
        const std::size_t z_face_high = geometry.faceIndexZ(i, j, k);
        const std::size_t z_face_low = i > 0 ? geometry.faceIndexZ(i - 1U, j, k) : z_face_high;
        // E_y on an x-face is the second tangential component; on a z-face it is
        // the second tangential component as well.
        const std::array<double, 4> estimates{
            signals.electric_t2_x[x_face_high],
            signals.electric_t2_x[x_face_low],
            signals.electric_t2_z[z_face_high],
            signals.electric_t2_z[z_face_low],
        };
        const std::array<double, 4> speeds{
            signals.face_speed_x[x_face_high],
            signals.face_speed_x[x_face_low],
            signals.face_speed_z[z_face_high],
            signals.face_speed_z[z_face_low],
        };
        const double field_component =
            blendEdge(estimates, speeds, floor, degenerate_edge_count, "E_y");
        emf.emfY()[geometry.edgeIndexY(i, j, k)] = field_component * geometry.dy_comoving;
      }
    }
  }

  // z-directed edges: (nx+1, ny+1, nz)
  for (std::size_t k = 0; k < nz; ++k) {
    for (std::size_t j = 0; j <= ny; ++j) {
      for (std::size_t i = 0; i <= nx; ++i) {
        const std::size_t x_face_high = geometry.faceIndexX(i, j, k);
        const std::size_t x_face_low = j > 0 ? geometry.faceIndexX(i, j - 1U, k) : x_face_high;
        const std::size_t y_face_high = geometry.faceIndexY(i, j, k);
        const std::size_t y_face_low = i > 0 ? geometry.faceIndexY(i - 1U, j, k) : y_face_high;
        // E_z on an x-face is the first tangential component; on a y-face it is
        // the first tangential component as well.
        const std::array<double, 4> estimates{
            signals.electric_t1_x[x_face_high],
            signals.electric_t1_x[x_face_low],
            signals.electric_t1_y[y_face_high],
            signals.electric_t1_y[y_face_low],
        };
        const std::array<double, 4> speeds{
            signals.face_speed_x[x_face_high],
            signals.face_speed_x[x_face_low],
            signals.face_speed_y[y_face_high],
            signals.face_speed_y[y_face_low],
        };
        const double field_component =
            blendEdge(estimates, speeds, floor, degenerate_edge_count, "E_z");
        emf.emfZ()[geometry.edgeIndexZ(i, j, k)] = field_component * geometry.dz_comoving;
      }
    }
  }
}

void applyInductionUpdate(const MhdPatchGeometry& geometry,
                          std::span<double> phi_x,
                          std::span<double> phi_y,
                          std::span<double> phi_z,
                          const MhdEdgeEmfStorage& emf) {
  // The edge integrals already carry the interval duration.
  applyConstrainedTransportUpdate(geometry, phi_x, phi_y, phi_z,
                                  emf.emfX(), emf.emfY(), emf.emfZ(), 1.0);
}

void accumulateEdgeEmfIntegral(const MhdPatchGeometry& geometry,
                               MhdEdgeEmfStorage& accumulator,
                               const MhdEdgeEmfStorage& interval_emf) {
  if (accumulator.geometry().cellCount() != geometry.cellCount()) {
    accumulator.resize(geometry);
  }
  const auto accumulate = [](std::span<double> target, std::span<const double> source) {
    const std::size_t count = std::min(target.size(), source.size());
    for (std::size_t index = 0; index < count; ++index) {
      target[index] += source[index];
    }
  };
  accumulate(accumulator.emfX(), interval_emf.emfX());
  accumulate(accumulator.emfY(), interval_emf.emfY());
  accumulate(accumulator.emfZ(), interval_emf.emfZ());
}

double magneticDivergenceNormalizer(std::span<const double> phi_x,
                                    std::span<const double> phi_y,
                                    std::span<const double> phi_z,
                                    double eps_scale) {
  double maximum = 0.0;
  const auto consider = [&maximum](std::span<const double> values) {
    for (const double value : values) {
      if (std::isfinite(value)) {
        maximum = std::max(maximum, std::abs(value));
      }
    }
  };
  consider(phi_x);
  consider(phi_y);
  consider(phi_z);
  const double floor = std::isfinite(eps_scale) && eps_scale > 0.0 ? eps_scale : 0.0;
  return std::max(maximum, floor);
}

}  // namespace cosmosim::mhd
