#include "cosmosim/mhd/mhd_induction.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <string>

namespace cosmosim::mhd {
namespace {

[[nodiscard]] std::array<std::size_t, 3> patchExtents(std::size_t nx, std::size_t ny,
                                                      std::size_t nz) noexcept {
  return {nx, ny, nz};
}

[[nodiscard]] double spacingAlong(const MhdPatchGeometry& geometry, std::size_t axis) noexcept {
  return axis == 0U ? geometry.dx_comoving
                    : (axis == 1U ? geometry.dy_comoving : geometry.dz_comoving);
}

[[nodiscard]] std::size_t edgeIndex(const MhdPatchGeometry& geometry, MhdAxis axis,
                                    const MhdIndex3& e) noexcept {
  const auto i = static_cast<std::size_t>(e[0]);
  const auto j = static_cast<std::size_t>(e[1]);
  const auto k = static_cast<std::size_t>(e[2]);
  switch (axis) {
    case MhdAxis::kX: return geometry.edgeIndexX(i, j, k);
    case MhdAxis::kY: return geometry.edgeIndexY(i, j, k);
    case MhdAxis::kZ: return geometry.edgeIndexZ(i, j, k);
  }
  return 0U;
}

[[nodiscard]] std::span<double> emfLane(MhdEdgeEmfStorage& emf, MhdAxis axis) noexcept {
  switch (axis) {
    case MhdAxis::kX: return emf.emfX();
    case MhdAxis::kY: return emf.emfY();
    case MhdAxis::kZ: return emf.emfZ();
  }
  return emf.emfX();
}

[[noreturn]] void throwInvalid(std::string detail) {
  MhdFailureContext context{};
  context.owner = "mhd_emf_builder";
  context.detail = std::move(detail);
  throw MhdError(MhdErrorCode::kInvalidField, std::move(context));
}

}  // namespace

// ---------------------------------------------------------------------------
// MhdPrimitivePatch
// ---------------------------------------------------------------------------

void MhdPrimitivePatch::resize(std::size_t patch_nx, std::size_t patch_ny,
                               std::size_t patch_nz, std::size_t ghost_cells) {
  nx = patch_nx;
  ny = patch_ny;
  nz = patch_nz;
  ghost = ghost_cells;
  px = nx + 2U * ghost;
  py = ny + 2U * ghost;
  pz = nz + 2U * ghost;
  const std::size_t count = px * py * pz;
  for (MhdDoubleArray* lane : {&rho, &vx, &vy, &vz, &pressure, &bx, &by, &bz, &metal}) {
    if (lane->size() != count) {
      lane->resize(count, 0.0);
    }
  }
}

bool MhdPrimitivePatch::contains(const MhdIndex3& c) const noexcept {
  const auto g = static_cast<std::ptrdiff_t>(ghost);
  return c[0] >= -g && c[0] < static_cast<std::ptrdiff_t>(nx) + g && c[1] >= -g &&
      c[1] < static_cast<std::ptrdiff_t>(ny) + g && c[2] >= -g &&
      c[2] < static_cast<std::ptrdiff_t>(nz) + g;
}

const MhdDoubleArray& MhdPrimitivePatch::velocity(MhdAxis axis) const noexcept {
  return axis == MhdAxis::kX ? vx : (axis == MhdAxis::kY ? vy : vz);
}

const MhdDoubleArray& MhdPrimitivePatch::field(MhdAxis axis) const noexcept {
  return axis == MhdAxis::kX ? bx : (axis == MhdAxis::kY ? by : bz);
}

MhdPrimitiveState MhdPrimitivePatch::load(const MhdIndex3& c) const noexcept {
  const std::size_t n = index(c);
  MhdPrimitiveState state{};
  state.rho_comoving = rho[n];
  state.vel_x_peculiar = vx[n];
  state.vel_y_peculiar = vy[n];
  state.vel_z_peculiar = vz[n];
  state.pressure_comoving = pressure[n];
  state.bx_comoving = bx[n];
  state.by_comoving = by[n];
  state.bz_comoving = bz[n];
  state.metallicity_mass_fraction = metal[n];
  return state;
}

void MhdPrimitivePatch::store(const MhdIndex3& c, const MhdPrimitiveState& state) noexcept {
  const std::size_t n = index(c);
  rho[n] = state.rho_comoving;
  vx[n] = state.vel_x_peculiar;
  vy[n] = state.vel_y_peculiar;
  vz[n] = state.vel_z_peculiar;
  pressure[n] = state.pressure_comoving;
  bx[n] = state.bx_comoving;
  by[n] = state.by_comoving;
  bz[n] = state.bz_comoving;
  metal[n] = state.metallicity_mass_fraction;
}

std::uint64_t MhdPrimitivePatch::ownedCapacityBytes() const noexcept {
  return rho.ownedCapacityBytes() + vx.ownedCapacityBytes() + vy.ownedCapacityBytes() +
      vz.ownedCapacityBytes() + pressure.ownedCapacityBytes() + bx.ownedCapacityBytes() +
      by.ownedCapacityBytes() + bz.ownedCapacityBytes() + metal.ownedCapacityBytes();
}

// ---------------------------------------------------------------------------
// MhdStaggeredFieldPatch
// ---------------------------------------------------------------------------

void MhdStaggeredFieldPatch::resize(std::size_t patch_nx, std::size_t patch_ny,
                                    std::size_t patch_nz, std::size_t ghost) {
  nx = patch_nx;
  ny = patch_ny;
  nz = patch_nz;
  tangential_ghost = ghost;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::array<std::size_t, 3> e = extents(static_cast<MhdAxis>(axis));
    const std::size_t count = e[0] * e[1] * e[2];
    if (normal_field[axis].size() != count) {
      normal_field[axis].resize(count, 0.0);
    }
  }
}

std::array<std::size_t, 3> MhdStaggeredFieldPatch::extents(MhdAxis normal) const noexcept {
  std::array<std::size_t, 3> e = patchExtents(nx, ny, nz);
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    e[axis] = axis == mhdAxisIndex(normal) ? e[axis] + 1U : e[axis] + 2U * tangential_ghost;
  }
  return e;
}

bool MhdStaggeredFieldPatch::contains(MhdAxis normal, const MhdIndex3& c) const noexcept {
  const std::array<std::size_t, 3> n = patchExtents(nx, ny, nz);
  const auto g = static_cast<std::ptrdiff_t>(tangential_ghost);
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const auto extent = static_cast<std::ptrdiff_t>(n[axis]);
    if (axis == mhdAxisIndex(normal)) {
      if (c[axis] < 0 || c[axis] > extent) {
        return false;
      }
    } else if (c[axis] < -g || c[axis] >= extent + g) {
      return false;
    }
  }
  return true;
}

std::size_t MhdStaggeredFieldPatch::index(MhdAxis normal, const MhdIndex3& c) const noexcept {
  const std::array<std::size_t, 3> e = extents(normal);
  const auto g = static_cast<std::ptrdiff_t>(tangential_ghost);
  std::array<std::size_t, 3> p{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    p[axis] = static_cast<std::size_t>(axis == mhdAxisIndex(normal) ? c[axis] : c[axis] + g);
  }
  return p[0] + e[0] * (p[1] + e[1] * p[2]);
}

std::uint64_t MhdStaggeredFieldPatch::ownedCapacityBytes() const noexcept {
  return normal_field[0].ownedCapacityBytes() + normal_field[1].ownedCapacityBytes() +
      normal_field[2].ownedCapacityBytes();
}

// ---------------------------------------------------------------------------
// MhdFaceSignalPatch
// ---------------------------------------------------------------------------

void MhdFaceSignalPatch::resize(const MhdPatchGeometry& patch_geometry) {
  nx = patch_geometry.nx;
  ny = patch_geometry.ny;
  nz = patch_geometry.nz;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::array<std::size_t, 3> e = extents(static_cast<MhdAxis>(axis));
    const std::size_t count = e[0] * e[1] * e[2];
    if (speed_plus[axis].size() != count) {
      speed_plus[axis].resize(count, 0.0);
      speed_minus[axis].resize(count, 0.0);
    }
    for (MhdDoubleArray& correction : induction_correction[axis]) {
      if (correction.size() != count) {
        correction.resize(count, 0.0);
      }
    }
  }
}

std::array<std::size_t, 3> MhdFaceSignalPatch::extents(MhdAxis normal) const noexcept {
  std::array<std::size_t, 3> e = patchExtents(nx, ny, nz);
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    e[axis] = axis == mhdAxisIndex(normal) ? e[axis] + 1U : e[axis] + 2U;
  }
  return e;
}

bool MhdFaceSignalPatch::contains(MhdAxis normal, const MhdIndex3& c) const noexcept {
  const std::array<std::size_t, 3> n = patchExtents(nx, ny, nz);
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const auto extent = static_cast<std::ptrdiff_t>(n[axis]);
    if (axis == mhdAxisIndex(normal)) {
      if (c[axis] < 0 || c[axis] > extent) {
        return false;
      }
    } else if (c[axis] < -1 || c[axis] > extent) {
      return false;
    }
  }
  return true;
}

std::size_t MhdFaceSignalPatch::index(MhdAxis normal, const MhdIndex3& c) const noexcept {
  const std::array<std::size_t, 3> e = extents(normal);
  std::array<std::size_t, 3> p{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    p[axis] = static_cast<std::size_t>(axis == mhdAxisIndex(normal) ? c[axis] : c[axis] + 1);
  }
  return p[0] + e[0] * (p[1] + e[1] * p[2]);
}

std::uint64_t MhdFaceSignalPatch::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    bytes += speed_plus[axis].ownedCapacityBytes() + speed_minus[axis].ownedCapacityBytes();
    for (const MhdDoubleArray& correction : induction_correction[axis]) {
      bytes += correction.ownedCapacityBytes();
    }
  }
  return bytes;
}

// ---------------------------------------------------------------------------
// Staggered ghost fill
// ---------------------------------------------------------------------------

void fillStaggeredGhosts(const MhdPatchBoundary& boundary, MhdStaggeredFieldPatch& staggered) {
  const std::array<std::size_t, 3> n = patchExtents(staggered.nx, staggered.ny, staggered.nz);
  const auto g = static_cast<std::ptrdiff_t>(staggered.tangential_ghost);
  for (std::size_t normal = 0; normal < 3U; ++normal) {
    const MhdAxis normal_axis = static_cast<MhdAxis>(normal);
    // Periodic exactness: both boundary planes of a self-periodic patch are one
    // logical face set and must carry identical values.
    if (boundary.lower[normal] == MhdBoundaryKind::kPeriodic) {
      const auto extent_n = static_cast<std::ptrdiff_t>(n[normal]);
      const std::size_t t1 = (normal + 1U) % 3U;
      const std::size_t t2 = (normal + 2U) % 3U;
      for (std::ptrdiff_t b = 0; b < static_cast<std::ptrdiff_t>(n[t2]); ++b) {
        for (std::ptrdiff_t a = 0; a < static_cast<std::ptrdiff_t>(n[t1]); ++a) {
          MhdIndex3 lower{};
          lower[normal] = 0;
          lower[t1] = a;
          lower[t2] = b;
          MhdIndex3 upper = lower;
          upper[normal] = extent_n;
          if (staggered.at(normal_axis, lower) != staggered.at(normal_axis, upper)) {
            throwInvalid("periodic patch carries different face fluxes on its two "
                         "boundary planes of the same logical face set");
          }
        }
      }
    }
    for (int pass = 1; pass <= 2; ++pass) {
      const std::size_t t = (normal + static_cast<std::size_t>(pass)) % 3U;
      const std::size_t other = (normal + static_cast<std::size_t>(3 - pass)) % 3U;
      const auto extent_t = static_cast<std::ptrdiff_t>(n[t]);
      const auto extent_n = static_cast<std::ptrdiff_t>(n[normal]);
      const auto extent_o = static_cast<std::ptrdiff_t>(n[other]);
      for (std::ptrdiff_t co = -g; co < extent_o + g; ++co) {
        for (std::ptrdiff_t cn = 0; cn <= extent_n; ++cn) {
          for (std::ptrdiff_t ct = -g; ct < extent_t + g; ++ct) {
            if (ct >= 0 && ct < extent_t) {
              continue;
            }
            const bool lower_side = ct < 0;
            const MhdBoundaryKind kind = lower_side ? boundary.lower[t] : boundary.upper[t];
            if (kind == MhdBoundaryKind::kNeighborPatch ||
                kind == MhdBoundaryKind::kImportedMpi) {
              continue;
            }
            std::ptrdiff_t source = ct;
            double sign = 1.0;
            switch (kind) {
              case MhdBoundaryKind::kPeriodic:
                source = ((ct % extent_t) + extent_t) % extent_t;
                break;
              case MhdBoundaryKind::kOpen:
                source = lower_side ? 0 : extent_t - 1;
                break;
              case MhdBoundaryKind::kReflectingConductingWall:
                source = lower_side ? (-1 - ct) : (2 * extent_t - 1 - ct);
                sign = -1.0;  // face component tangential to the wall is odd
                break;
              default:
                break;
            }
            MhdIndex3 target{};
            target[normal] = cn;
            target[t] = ct;
            target[other] = co;
            MhdIndex3 origin = target;
            origin[t] = source;
            staggered.at(normal_axis, target) = sign * staggered.at(normal_axis, origin);
          }
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// UCT-HLL
// ---------------------------------------------------------------------------

void MhdEmfBuilder::buildUctHll(const MhdPatchGeometry& geometry,
                                const MhdPatchBoundary& boundary,
                                const MhdPrimitivePatch& primitives,
                                const MhdStaggeredFieldPatch& staggered,
                                const MhdFaceSignalPatch& signals,
                                double scale_factor,
                                MhdEdgeEmfStorage& emf,
                                Diagnostics* diagnostics) {
  if (!(scale_factor > 0.0) || !std::isfinite(scale_factor)) {
    throwInvalid("UCT-HLL requires a finite positive stage scale factor");
  }
  if (primitives.ghost < 2U || staggered.tangential_ghost < 2U ||
      primitives.nx != geometry.nx || primitives.ny != geometry.ny ||
      primitives.nz != geometry.nz || staggered.nx != geometry.nx ||
      staggered.ny != geometry.ny || staggered.nz != geometry.nz ||
      signals.nx != geometry.nx || signals.ny != geometry.ny || signals.nz != geometry.nz) {
    throwInvalid("UCT-HLL stage buffers do not match the patch geometry or ghost depth");
  }
  if (emf.geometry().nx != geometry.nx || emf.geometry().ny != geometry.ny ||
      emf.geometry().nz != geometry.nz) {
    emf.resize(geometry);
  }
  const double inverse_sqrt_a = 1.0 / std::sqrt(scale_factor);
  const std::array<std::size_t, 3> n{geometry.nx, geometry.ny, geometry.nz};
  std::uint64_t built = 0U;

  for (std::size_t c = 0; c < 3U; ++c) {
    const MhdAxis edge_axis = static_cast<MhdAxis>(c);
    const std::size_t a = (c + 1U) % 3U;
    const std::size_t b = (c + 2U) % 3U;
    const MhdAxis axis_a = static_cast<MhdAxis>(a);
    const MhdAxis axis_b = static_cast<MhdAxis>(b);
    const double edge_length = spacingAlong(geometry, c);
    const MhdDoubleArray& velocity_a = primitives.velocity(axis_a);
    const MhdDoubleArray& velocity_b = primitives.velocity(axis_b);
    std::span<double> target = emfLane(emf, edge_axis);

    const auto cellSlope = [&](const MhdDoubleArray& lane, MhdIndex3 cell, std::size_t along) {
      MhdIndex3 minus = cell;
      MhdIndex3 plus = cell;
      minus[along] -= 1;
      plus[along] += 1;
      return mcLimitedSlope(lane[primitives.index(minus)], lane[primitives.index(cell)],
                            lane[primitives.index(plus)]);
    };
    const auto faceToward = [&](MhdAxis normal, MhdIndex3 face, std::size_t along,
                                double direction) {
      MhdIndex3 minus = face;
      MhdIndex3 plus = face;
      minus[along] -= 1;
      plus[along] += 1;
      const double center = staggered.at(normal, face);
      return center + direction * 0.5 *
          mcLimitedSlope(staggered.at(normal, minus), center, staggered.at(normal, plus));
    };

    for (std::ptrdiff_t kc = 0; kc < static_cast<std::ptrdiff_t>(n[c]); ++kc) {
      for (std::ptrdiff_t nb = 0; nb <= static_cast<std::ptrdiff_t>(n[b]); ++nb) {
        for (std::ptrdiff_t na = 0; na <= static_cast<std::ptrdiff_t>(n[a]); ++na) {
          // a-normal faces adjacent to the edge, below/above in b.
          MhdIndex3 face_a_lower{};
          face_a_lower[c] = kc;
          face_a_lower[a] = na;
          face_a_lower[b] = nb - 1;
          MhdIndex3 face_a_upper = face_a_lower;
          face_a_upper[b] = nb;
          // b-normal faces adjacent to the edge, below/above in a.
          MhdIndex3 face_b_lower{};
          face_b_lower[c] = kc;
          face_b_lower[b] = nb;
          face_b_lower[a] = na - 1;
          MhdIndex3 face_b_upper = face_b_lower;
          face_b_upper[a] = na;

          const std::size_t sa_lo = signals.index(axis_a, face_a_lower);
          const std::size_t sa_hi = signals.index(axis_a, face_a_upper);
          const std::size_t sb_lo = signals.index(axis_b, face_b_lower);
          const std::size_t sb_hi = signals.index(axis_b, face_b_upper);
          const double ap_a = std::max(signals.speed_plus[a][sa_lo], signals.speed_plus[a][sa_hi]);
          const double am_a =
              std::max(signals.speed_minus[a][sa_lo], signals.speed_minus[a][sa_hi]);
          const double ap_b = std::max(signals.speed_plus[b][sb_lo], signals.speed_plus[b][sb_hi]);
          const double am_b =
              std::max(signals.speed_minus[b][sb_lo], signals.speed_minus[b][sb_hi]);
          const double denom_a = ap_a + am_a;
          const double denom_b = ap_b + am_b;
          if (!(denom_a > 0.0) || !(denom_b > 0.0) || !std::isfinite(denom_a) ||
              !std::isfinite(denom_b)) {
            std::ostringstream stream;
            stream << "UCT-HLL edge along axis " << c << " at (" << na << ',' << nb << ','
                   << kc << ") has a degenerate signal fan (a_a+=" << ap_a << ", a_a-="
                   << am_a << ", a_b+=" << ap_b << ", a_b-=" << am_b << ')';
            throwInvalid(stream.str());
          }

          // Staggered fields reconstructed to the edge.
          const double field_a_lower = faceToward(axis_a, face_a_lower, b, +1.0);
          const double field_a_upper = faceToward(axis_a, face_a_upper, b, -1.0);
          const double field_b_lower = faceToward(axis_b, face_b_lower, a, +1.0);
          const double field_b_upper = faceToward(axis_b, face_b_upper, a, -1.0);

          // Quadrant electric fields from corner-extrapolated velocities.
          std::array<double, 4> quadrant{};  // order (--), (-+), (+-), (++) in (sa, sb)
          for (int qa = 0; qa < 2; ++qa) {
            for (int qb = 0; qb < 2; ++qb) {
              MhdIndex3 cell{};
              cell[c] = kc;
              cell[a] = na - 1 + qa;
              cell[b] = nb - 1 + qb;
              const double toward_a = qa == 0 ? 0.5 : -0.5;
              const double toward_b = qb == 0 ? 0.5 : -0.5;
              const std::size_t center = primitives.index(cell);
              const double corner_va = velocity_a[center] +
                  toward_a * cellSlope(velocity_a, cell, a) +
                  toward_b * cellSlope(velocity_a, cell, b);
              const double corner_vb = velocity_b[center] +
                  toward_a * cellSlope(velocity_b, cell, a) +
                  toward_b * cellSlope(velocity_b, cell, b);
              const double field_b = qa == 0 ? field_b_lower : field_b_upper;
              const double field_a = qb == 0 ? field_a_lower : field_a_upper;
              quadrant[static_cast<std::size_t>(2 * qa + qb)] =
                  -(corner_va * field_b - corner_vb * field_a);
            }
          }
          const double upwinded =
              (ap_a * ap_b * quadrant[0] + ap_a * am_b * quadrant[1] +
               am_a * ap_b * quadrant[2] + am_a * am_b * quadrant[3]) /
              (denom_a * denom_b);
          const double dissipation_a = ap_a * am_a / denom_a * (field_b_upper - field_b_lower);
          const double dissipation_b = ap_b * am_b / denom_b * (field_a_upper - field_a_lower);
          double electric_eff = upwinded + dissipation_a - dissipation_b;
          if (signals.induction_correction_valid) {
            // a-faces: transverse component B_b is slot (b - a - 1) % 3 == 0;
            // b-faces: B_a is slot (a - b - 1) % 3 == 1 (cyclic a = c+1, b = c+2).
            const MhdDoubleArray& correction_a = signals.induction_correction[a][0];
            const MhdDoubleArray& correction_b = signals.induction_correction[b][1];
            const double face_a_correction =
                (ap_b * correction_a[sa_lo] + am_b * correction_a[sa_hi]) / denom_b;
            const double face_b_correction =
                (ap_a * correction_b[sb_lo] + am_a * correction_b[sb_hi]) / denom_a;
            electric_eff += -face_a_correction + face_b_correction;
          }
          if (!std::isfinite(electric_eff)) {
            throwInvalid("UCT-HLL produced a non-finite edge electric field");
          }
          MhdIndex3 edge{};
          edge[c] = kc;
          edge[a] = na;
          edge[b] = nb;
          target[edgeIndex(geometry, edge_axis, edge)] =
              electric_eff * inverse_sqrt_a * edge_length;
          ++built;
        }
      }
    }
  }
  const std::uint64_t zeroed = applyEdgeBoundaryRules(geometry, boundary, emf);
  if (diagnostics != nullptr) {
    diagnostics->edges_built += built;
    diagnostics->wall_edges_zeroed += zeroed;
  }
}

std::uint64_t MhdEmfBuilder::applyEdgeBoundaryRules(const MhdPatchGeometry& geometry,
                                                    const MhdPatchBoundary& boundary,
                                                    MhdEdgeEmfStorage& emf) {
  const std::array<std::size_t, 3> n{geometry.nx, geometry.ny, geometry.nz};
  std::uint64_t zeroed = 0U;
  for (std::size_t c = 0; c < 3U; ++c) {
    const MhdAxis edge_axis = static_cast<MhdAxis>(c);
    std::span<double> lane = emfLane(emf, edge_axis);
    for (int pass = 1; pass <= 2; ++pass) {
      const std::size_t t = (c + static_cast<std::size_t>(pass)) % 3U;
      const std::size_t o = (c + static_cast<std::size_t>(3 - pass)) % 3U;
      const auto extent_t = static_cast<std::ptrdiff_t>(n[t]);
      for (std::ptrdiff_t kc = 0; kc < static_cast<std::ptrdiff_t>(n[c]); ++kc) {
        for (std::ptrdiff_t no = 0; no <= static_cast<std::ptrdiff_t>(n[o]); ++no) {
          MhdIndex3 lower{};
          lower[c] = kc;
          lower[o] = no;
          lower[t] = 0;
          MhdIndex3 upper = lower;
          upper[t] = extent_t;
          const std::size_t lower_index = edgeIndex(geometry, edge_axis, lower);
          const std::size_t upper_index = edgeIndex(geometry, edge_axis, upper);
          if (boundary.lower[t] == MhdBoundaryKind::kPeriodic) {
            lane[upper_index] = lane[lower_index];
          }
          // Perfectly conducting wall: tangential E vanishes on the wall plane.
          if (boundary.lower[t] == MhdBoundaryKind::kReflectingConductingWall) {
            lane[lower_index] = 0.0;
            ++zeroed;
          }
          if (boundary.upper[t] == MhdBoundaryKind::kReflectingConductingWall) {
            lane[upper_index] = 0.0;
            ++zeroed;
          }
        }
      }
    }
  }
  return zeroed;
}

void applyInductionUpdate(const MhdPatchGeometry& geometry,
                          std::span<double> phi_x,
                          std::span<double> phi_y,
                          std::span<double> phi_z,
                          const MhdEdgeEmfStorage& emf_rate,
                          double weight) {
  applyConstrainedTransportUpdate(geometry, phi_x, phi_y, phi_z, emf_rate.emfX(),
                                  emf_rate.emfY(), emf_rate.emfZ(), weight);
}

void accumulateEdgeEmfIntegral(const MhdPatchGeometry& geometry,
                               MhdEdgeEmfStorage& accumulator,
                               const MhdEdgeEmfStorage& interval_emf,
                               double weight) {
  if (accumulator.geometry().nx != geometry.nx || accumulator.geometry().ny != geometry.ny ||
      accumulator.geometry().nz != geometry.nz) {
    throwInvalid("edge EMF accumulator geometry does not match the patch");
  }
  if (interval_emf.emfX().size() != accumulator.emfX().size() ||
      interval_emf.emfY().size() != accumulator.emfY().size() ||
      interval_emf.emfZ().size() != accumulator.emfZ().size()) {
    throwInvalid("edge EMF interval geometry does not match the accumulator");
  }
  const auto accumulate = [weight](std::span<double> target, std::span<const double> source) {
    for (std::size_t index = 0; index < target.size(); ++index) {
      target[index] += weight * source[index];
    }
  };
  accumulate(accumulator.emfX(), interval_emf.emfX());
  accumulate(accumulator.emfY(), interval_emf.emfY());
  accumulate(accumulator.emfZ(), interval_emf.emfZ());
}

double magneticDivergenceNormalizer(const MhdPatchGeometry& geometry,
                                    std::span<const double> phi_x,
                                    std::span<const double> phi_y,
                                    std::span<const double> phi_z) {
  double maximum = 0.0;
  const auto consider = [&maximum](std::span<const double> values) {
    for (const double value : values) {
      if (!std::isfinite(value)) {
        maximum = std::numeric_limits<double>::infinity();
        return;
      }
      maximum = std::max(maximum, std::abs(value));
    }
  };
  consider(phi_x);
  consider(phi_y);
  consider(phi_z);
  const double area_ref =
      std::max({geometry.faceAreaX(), geometry.faceAreaY(), geometry.faceAreaZ()});
  const double eps_scale = k_mhd_b_zero_regularization_fraction * area_ref;
  return std::max(maximum, eps_scale);
}

}  // namespace cosmosim::mhd
