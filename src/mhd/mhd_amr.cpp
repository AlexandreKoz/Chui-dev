#include "cosmosim/mhd/mhd_amr.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace cosmosim::mhd {
namespace {

MhdFailureContext makeContext(std::string_view owner, std::string detail) {
  MhdFailureContext context{};
  context.owner = std::string(owner);
  context.detail = std::move(detail);
  return context;
}

[[noreturn]] void throwRegister(std::string detail) {
  throw MhdError(MhdErrorCode::kAmrRegisterIncomplete,
                 makeContext("mhd_emf_register", std::move(detail)));
}

[[nodiscard]] bool evenExtents(const MhdPatchGeometry& geometry) noexcept {
  return geometry.nx % 2U == 0U && geometry.ny % 2U == 0U && geometry.nz % 2U == 0U &&
      geometry.nx >= 2U && geometry.ny >= 2U && geometry.nz >= 2U;
}

[[nodiscard]] std::size_t faceIndex(const MhdPatchGeometry& geometry, std::size_t axis,
                                    std::size_t i, std::size_t j, std::size_t k) noexcept {
  return axis == 0U ? geometry.faceIndexX(i, j, k)
                    : (axis == 1U ? geometry.faceIndexY(i, j, k) : geometry.faceIndexZ(i, j, k));
}

[[nodiscard]] double faceArea(const MhdPatchGeometry& geometry, std::size_t axis) noexcept {
  return axis == 0U ? geometry.faceAreaX()
                    : (axis == 1U ? geometry.faceAreaY() : geometry.faceAreaZ());
}

[[nodiscard]] std::span<double> lane(MhdFaceFluxStorage& storage, std::size_t axis) noexcept {
  return axis == 0U ? storage.fluxX() : (axis == 1U ? storage.fluxY() : storage.fluxZ());
}

[[nodiscard]] std::span<const double> lane(const MhdFaceFluxStorage& storage,
                                           std::size_t axis) noexcept {
  return axis == 0U ? storage.fluxX() : (axis == 1U ? storage.fluxY() : storage.fluxZ());
}

[[nodiscard]] bool allFinite(std::span<const double> values) noexcept {
  return std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); });
}

// Cell-centred comoving field from face fluxes.
[[nodiscard]] std::array<double, 3> cellField(const MhdPatchGeometry& geometry,
                                              std::span<const double> phi_x,
                                              std::span<const double> phi_y,
                                              std::span<const double> phi_z, std::size_t i,
                                              std::size_t j, std::size_t k) noexcept {
  return {0.5 * (phi_x[geometry.faceIndexX(i, j, k)] + phi_x[geometry.faceIndexX(i + 1U, j, k)]) /
              geometry.faceAreaX(),
          0.5 * (phi_y[geometry.faceIndexY(i, j, k)] + phi_y[geometry.faceIndexY(i, j + 1U, k)]) /
              geometry.faceAreaY(),
          0.5 * (phi_z[geometry.faceIndexZ(i, j, k)] + phi_z[geometry.faceIndexZ(i, j, k + 1U)]) /
              geometry.faceAreaZ()};
}

[[nodiscard]] double thermalEnergy(const MhdConservedState& gas, const std::array<double, 3>& b,
                                   double scale_factor) noexcept {
  if (!(gas.rho_comoving > 0.0)) {
    return -1.0;
  }
  const double momentum2 = gas.momentum_density_x_comoving * gas.momentum_density_x_comoving +
      gas.momentum_density_y_comoving * gas.momentum_density_y_comoving +
      gas.momentum_density_z_comoving * gas.momentum_density_z_comoving;
  return gas.total_energy_density_comoving - 0.5 * momentum2 / gas.rho_comoving -
      0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / scale_factor;
}

// Face-local limited transverse slope of the parent face field (normalized to
// the half-cell coordinate in [-1, 1]); zero at the patch edge (one-sided
// lower-order closure) or when `zero_slopes`.
struct ParentFaceField {
  const MhdPatchGeometry* geometry = nullptr;
  std::span<const double> phi;
  std::size_t axis = 0;
  bool zero_slopes = false;

  [[nodiscard]] double value(std::array<std::size_t, 3> c) const noexcept {
    return phi[faceIndex(*geometry, axis, c[0], c[1], c[2])] / faceArea(*geometry, axis);
  }
  [[nodiscard]] double slope(std::array<std::size_t, 3> c, std::size_t along) const noexcept {
    if (zero_slopes) {
      return 0.0;
    }
    const std::array<std::size_t, 3> n{geometry->nx, geometry->ny, geometry->nz};
    if (c[along] == 0U || c[along] + 1U >= n[along]) {
      return 0.0;
    }
    std::array<std::size_t, 3> minus = c;
    std::array<std::size_t, 3> plus = c;
    minus[along] -= 1U;
    plus[along] += 1U;
    // MC slope per parent cell, halved to the normalized coordinate.
    return 0.5 * mcLimitedSlope(value(minus), value(c), value(plus));
  }
};

}  // namespace

std::uint64_t mhdEmfRegisterBound(const MhdPatchGeometry& coarse_geometry) {
  const std::uint64_t nx = coarse_geometry.nx;
  const std::uint64_t ny = coarse_geometry.ny;
  const std::uint64_t nz = coarse_geometry.nz;
  // Edges lying on the six boundary planes, counted once (inclusion-exclusion
  // over the patch boundary shell).
  const std::uint64_t all = coarse_geometry.edgeCountTotal();
  const std::uint64_t interior = (nx >= 2U && ny >= 2U && nz >= 2U)
      ? nx * (ny - 1U) * (nz - 1U) + (nx - 1U) * ny * (nz - 1U) + (nx - 1U) * (ny - 1U) * nz
      : 0U;
  return all - interior;
}

// ---------------------------------------------------------------------------
// MhdEmfRegister
// ---------------------------------------------------------------------------

void MhdEmfRegister::bindInterval(const MhdAmrIntervalId& interval,
                                  std::uint64_t coarse_patch_id,
                                  std::uint64_t coarse_geometry_generation,
                                  std::int32_t fine_substeps) {
  if (fine_substeps < 1 || fine_substeps > 8) {
    throwRegister("fine substep count must lie in [1, 8]");
  }
  if (!(interval.time_end_code > interval.time_begin_code) ||
      interval.coarse_geometry_generation != coarse_geometry_generation) {
    throwRegister("register interval is empty or its geometry generation is stale");
  }
  m_records.clear();
  m_interval = interval;
  m_coarse_patch_id = coarse_patch_id;
  m_coarse_geometry_generation = coarse_geometry_generation;
  m_fine_substeps = fine_substeps;
  m_bound = true;
  m_declared = false;
}

void MhdEmfRegister::declareInterfaceEdges(std::span<const MhdCoarseEdgeKey> keys) {
  if (!m_bound || m_declared) {
    throwRegister("interface edges must be declared exactly once after bindInterval");
  }
  m_records.reserve(keys.size());
  for (const MhdCoarseEdgeKey& key : keys) {
    if (key.coarse_patch_id != m_coarse_patch_id) {
      throwRegister("declared edge belongs to a different coarse patch");
    }
    Record record{};
    record.key = key;
    m_records.push_back(record);
  }
  std::sort(m_records.begin(), m_records.end(),
            [](const Record& lhs, const Record& rhs) { return lhs.key < rhs.key; });
  for (std::size_t index = 1U; index < m_records.size(); ++index) {
    if (m_records[index - 1U].key == m_records[index].key) {
      throwRegister("duplicate interface edge declaration");
    }
  }
  m_declared = true;
}

MhdEmfRegister::Record& MhdEmfRegister::locate(const MhdAmrIntervalId& interval,
                                               const MhdCoarseEdgeKey& key) {
  if (!m_bound || !m_declared) {
    throwRegister("contribution before the register interval/edge set is bound");
  }
  if (!(interval == m_interval)) {
    throwRegister("contribution for a stale or different synchronization interval");
  }
  const auto found = std::lower_bound(
      m_records.begin(), m_records.end(), key,
      [](const Record& record, const MhdCoarseEdgeKey& value) { return record.key < value; });
  if (found == m_records.end() || !(found->key == key)) {
    throwRegister("contribution for an undeclared coarse edge");
  }
  return *found;
}

void MhdEmfRegister::addCoarse(const MhdAmrIntervalId& interval, const MhdCoarseEdgeKey& key,
                               double integrated_emf) {
  if (!std::isfinite(integrated_emf)) {
    throwRegister("non-finite coarse EMF integral");
  }
  Record& record = locate(interval, key);
  if (record.coarse_present) {
    throwRegister("duplicate coarse EMF contribution");
  }
  record.integrated_coarse = integrated_emf;
  record.coarse_present = true;
}

void MhdEmfRegister::addFine(const MhdAmrIntervalId& interval, const MhdCoarseEdgeKey& key,
                             std::int32_t fine_segment, std::int32_t fine_substep,
                             double integrated_emf) {
  if (!std::isfinite(integrated_emf)) {
    throwRegister("non-finite fine EMF integral");
  }
  if (fine_segment < 0 || fine_segment > 1 || fine_substep < 0 || fine_substep >= m_fine_substeps) {
    throwRegister("fine contribution segment/substep outside the declared coverage");
  }
  Record& record = locate(interval, key);
  const std::uint32_t bit = 1U << static_cast<std::uint32_t>(fine_segment * m_fine_substeps +
                                                               fine_substep);
  if ((record.fine_coverage_mask & bit) != 0U) {
    throwRegister("duplicate fine EMF contribution");
  }
  record.fine_coverage_mask |= bit;
  record.integrated_fine += integrated_emf;
}

void MhdEmfRegister::finalize() const {
  if (!m_bound || !m_declared) {
    throwRegister("finalize before bind/declare");
  }
  const std::uint32_t full = (1U << static_cast<std::uint32_t>(2 * m_fine_substeps)) - 1U;
  for (const Record& record : m_records) {
    if (!record.coarse_present || record.fine_coverage_mask != full) {
      std::ostringstream stream;
      stream << "incomplete interval coverage at coarse edge (axis "
             << static_cast<int>(record.key.axis) << ", " << record.key.i << ',' << record.key.j
             << ',' << record.key.k << "): coarse=" << record.coarse_present
             << " fine_mask=" << record.fine_coverage_mask << " expected=" << full;
      throwRegister(stream.str());
    }
  }
}

void MhdEmfRegister::clear() noexcept {
  m_records.clear();
  m_bound = false;
  m_declared = false;
}

std::uint64_t MhdEmfRegister::ownedCapacityBytes() const noexcept {
  return static_cast<std::uint64_t>(m_records.capacity() * sizeof(Record));
}

// ---------------------------------------------------------------------------
// MhdAmrSynchronizer
// ---------------------------------------------------------------------------

MhdPatchGeometry MhdAmrSynchronizer::childGeometry(const MhdPatchGeometry& parent,
                                                   std::int32_t child_index) {
  if (child_index < 0 || child_index > 7) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_amr", "child index outside [0,7]"));
  }
  if (!parent.valid() || !evenExtents(parent)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_amr", "ratio-two transfer requires a valid even-extent patch"));
  }
  const std::int32_t di = child_index % 2;
  const std::int32_t dj = (child_index / 2) % 2;
  const std::int32_t dk = child_index / 4;
  MhdPatchGeometry child{};
  child.nx = parent.nx;
  child.ny = parent.ny;
  child.nz = parent.nz;
  child.dx_comoving = 0.5 * parent.dx_comoving;
  child.dy_comoving = 0.5 * parent.dy_comoving;
  child.dz_comoving = 0.5 * parent.dz_comoving;
  child.origin_x_comoving = parent.origin_x_comoving +
      static_cast<double>(di) * 0.5 * static_cast<double>(parent.nx) * parent.dx_comoving;
  child.origin_y_comoving = parent.origin_y_comoving +
      static_cast<double>(dj) * 0.5 * static_cast<double>(parent.ny) * parent.dy_comoving;
  child.origin_z_comoving = parent.origin_z_comoving +
      static_cast<double>(dk) * 0.5 * static_cast<double>(parent.nz) * parent.dz_comoving;
  return child;
}

void MhdAmrSynchronizer::restrictFaceFluxes(const MhdPatchGeometry& parent_geometry,
                                            std::span<double> parent_phi_x,
                                            std::span<double> parent_phi_y,
                                            std::span<double> parent_phi_z,
                                            const std::array<MhdFaceFluxStorage, 8>& children) {
  if (!evenExtents(parent_geometry) || parent_phi_x.size() != parent_geometry.faceCountX() ||
      parent_phi_y.size() != parent_geometry.faceCountY() ||
      parent_phi_z.size() != parent_geometry.faceCountZ()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_amr", "restriction target does not match the parent"));
  }
  for (std::int32_t block = 0; block < 8; ++block) {
    const MhdFaceFluxStorage& child = children[static_cast<std::size_t>(block)];
    const MhdPatchGeometry expected = childGeometry(parent_geometry, block);
    if (child.geometry().nx != expected.nx || child.geometry().ny != expected.ny ||
        child.geometry().nz != expected.nz || child.fluxX().size() != expected.faceCountX() ||
        child.fluxY().size() != expected.faceCountY() ||
        child.fluxZ().size() != expected.faceCountZ() || !allFinite(child.fluxX()) ||
        !allFinite(child.fluxY()) || !allFinite(child.fluxZ())) {
      throw MhdError(MhdErrorCode::kInvalidField,
                     makeContext("mhd_amr", "child face storage has the wrong shape or "
                                            "non-finite values"));
    }
  }
  const std::array<std::size_t, 3> n{parent_geometry.nx, parent_geometry.ny, parent_geometry.nz};
  const std::array<std::size_t, 3> half{n[0] / 2U, n[1] / 2U, n[2] / 2U};
  // Compute into a temporary first so a sibling mismatch never leaves a
  // partially written parent.
  std::array<std::vector<double>, 3> result{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::size_t t1 = (axis + 1U) % 3U;
    const std::size_t t2 = (axis + 2U) % 3U;
    result[axis].assign(axis == 0U ? parent_geometry.faceCountX()
                                   : (axis == 1U ? parent_geometry.faceCountY()
                                                 : parent_geometry.faceCountZ()),
                        0.0);
    for (std::size_t b = 0; b < n[t2]; ++b) {
      for (std::size_t a = 0; a < n[t1]; ++a) {
        for (std::size_t f = 0; f <= n[axis]; ++f) {
          std::array<std::size_t, 3> parent_face{};
          parent_face[axis] = f;
          parent_face[t1] = a;
          parent_face[t2] = b;
          const auto sum_block = [&](std::size_t block_normal, std::size_t child_normal) {
            std::array<std::size_t, 3> block{};
            block[axis] = block_normal;
            block[t1] = a / half[t1];
            block[t2] = b / half[t2];
            const std::size_t child_index = block[0] + 2U * (block[1] + 2U * block[2]);
            const MhdFaceFluxStorage& child = children[child_index];
            const std::span<const double> values = lane(child, axis);
            double sum = 0.0;
            for (std::size_t q2 = 0; q2 < 2U; ++q2) {
              for (std::size_t q1 = 0; q1 < 2U; ++q1) {
                std::array<std::size_t, 3> c{};
                c[axis] = child_normal;
                c[t1] = 2U * (a - block[t1] * half[t1]) + q1;
                c[t2] = 2U * (b - block[t2] * half[t2]) + q2;
                sum += values[faceIndex(child.geometry(), axis, c[0], c[1], c[2])];
              }
            }
            return sum;
          };
          double value = 0.0;
          if (f < half[axis]) {
            value = sum_block(0U, 2U * f);
          } else if (f > half[axis]) {
            value = sum_block(1U, 2U * (f - half[axis]));
          } else {
            // Internal sibling interface: one authoritative read, verified.
            value = sum_block(1U, 0U);
            const double mirror = sum_block(0U, n[axis]);
            const double scale = std::max({std::abs(value), std::abs(mirror),
                                           std::numeric_limits<double>::min()});
            if (std::abs(value - mirror) > 1.0e-13 * scale) {
              throw MhdError(MhdErrorCode::kInvalidField,
                             makeContext("mhd_amr",
                                         "sibling children disagree on a shared interface face"));
            }
          }
          result[axis][faceIndex(parent_geometry, axis, parent_face[0], parent_face[1],
                                 parent_face[2])] = value;
        }
      }
    }
  }
  std::copy(result[0].begin(), result[0].end(), parent_phi_x.begin());
  std::copy(result[1].begin(), result[1].end(), parent_phi_y.begin());
  std::copy(result[2].begin(), result[2].end(), parent_phi_z.begin());
}

MhdAmrSynchronizer::TransferResult MhdAmrSynchronizer::prolongFaceFluxes(
    const MhdPatchGeometry& parent_geometry,
    std::span<const double> parent_phi_x,
    std::span<const double> parent_phi_y,
    std::span<const double> parent_phi_z,
    std::array<MhdFaceFluxStorage, 8>& children,
    double divergence_tolerance,
    bool zero_slopes) {
  TransferResult result{};
  result.used_zero_slope_fallback = zero_slopes;
  if (!parent_geometry.valid() || !evenExtents(parent_geometry)) {
    result.detail = "ratio-two prolongation requires a valid even-extent parent";
    return result;
  }
  if (parent_phi_x.size() != parent_geometry.faceCountX() ||
      parent_phi_y.size() != parent_geometry.faceCountY() ||
      parent_phi_z.size() != parent_geometry.faceCountZ() || !allFinite(parent_phi_x) ||
      !allFinite(parent_phi_y) || !allFinite(parent_phi_z)) {
    result.detail = "parent face lanes are mis-sized or non-finite";
    return result;
  }
  const double parent_normalizer =
      magneticDivergenceNormalizer(parent_geometry, parent_phi_x, parent_phi_y, parent_phi_z);
  if (maxNormalizedFaceDivergence(parent_geometry, parent_phi_x, parent_phi_y, parent_phi_z,
                                  parent_normalizer) > divergence_tolerance) {
    result.detail = "parent field violates the discrete divergence contract";
    return result;
  }
  std::array<MhdFaceFluxStorage, 8> staged{};
  for (std::int32_t block = 0; block < 8; ++block) {
    staged[static_cast<std::size_t>(block)].resize(childGeometry(parent_geometry, block));
  }
  const std::array<std::size_t, 3> n{parent_geometry.nx, parent_geometry.ny, parent_geometry.nz};
  const std::array<std::size_t, 3> half{n[0] / 2U, n[1] / 2U, n[2] / 2U};
  const std::array<ParentFaceField, 3> parent{
      ParentFaceField{&parent_geometry, parent_phi_x, 0U, zero_slopes},
      ParentFaceField{&parent_geometry, parent_phi_y, 1U, zero_slopes},
      ParentFaceField{&parent_geometry, parent_phi_z, 2U, zero_slopes}};

  for (std::size_t k = 0; k < n[2]; ++k) {
    for (std::size_t j = 0; j < n[1]; ++j) {
      for (std::size_t i = 0; i < n[0]; ++i) {
        const std::array<std::size_t, 3> cell{i, j, k};
        // Polynomial coefficients per component (normal axis d, transverse
        // axes t1 = d+1, t2 = d+2): mean +- normal coefficient, slopes.
        std::array<double, 3> normal_coefficient{};   // a_x, b_y, c_z
        std::array<double, 3> mean{};                 // (B+ + B-)/2
        std::array<std::array<double, 3>, 3> slope_sum{};   // [d][t] = (S_t+ + S_t-)/2
        std::array<std::array<double, 3>, 3> slope_diff{};  // [d][t] = (S_t+ - S_t-)/2
        for (std::size_t d = 0; d < 3U; ++d) {
          std::array<std::size_t, 3> lower = cell;
          std::array<std::size_t, 3> upper = cell;
          upper[d] += 1U;
          const double b_minus = parent[d].value(lower);
          const double b_plus = parent[d].value(upper);
          mean[d] = 0.5 * (b_plus + b_minus);
          normal_coefficient[d] = 0.5 * (b_plus - b_minus);
          for (std::size_t t = 0; t < 3U; ++t) {
            if (t == d) {
              continue;
            }
            const double s_minus = parent[d].slope(lower, t);
            const double s_plus = parent[d].slope(upper, t);
            slope_sum[d][t] = 0.5 * (s_plus + s_minus);
            slope_diff[d][t] = 0.5 * (s_plus - s_minus);
          }
        }
        // Quadratic coefficients from the solenoidal constraints:
        //   2 a_xx + b_xy + c_xz = 0 (and cyclic), with b_xy = slope_diff[y][x].
        std::array<double, 3> quadratic{};
        for (std::size_t d = 0; d < 3U; ++d) {
          const std::size_t t1 = (d + 1U) % 3U;
          const std::size_t t2 = (d + 2U) % 3U;
          quadratic[d] = -0.5 * (slope_diff[t1][d] + slope_diff[t2][d]);
        }
        // Write fine faces of this parent cell.
        for (std::size_t d = 0; d < 3U; ++d) {
          const std::size_t t1 = (d + 1U) % 3U;
          const std::size_t t2 = (d + 2U) % 3U;
          const double fine_area = 0.25 * faceArea(parent_geometry, d);
          std::array<std::size_t, 3> block{};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            block[axis] = cell[axis] / half[axis];
          }
          const std::size_t child_index = block[0] + 2U * (block[1] + 2U * block[2]);
          MhdFaceFluxStorage& child = staged[child_index];
          std::span<double> values = lane(child, d);
          for (int position = -1; position <= 1; ++position) {
            for (int q2 = 0; q2 < 2; ++q2) {
              for (int q1 = 0; q1 < 2; ++q1) {
                const double y1 = q1 == 0 ? -0.5 : 0.5;
                const double y2 = q2 == 0 ? -0.5 : 0.5;
                double field_value = 0.0;
                if (position == 0) {
                  // Internal fine face: polynomial at X = 0.
                  field_value = (mean[d] - quadratic[d]) + slope_sum[d][t1] * y1 +
                      slope_sum[d][t2] * y2;
                } else {
                  // Fine subface of a parent face: face-local linear
                  // reconstruction (identical from both adjacent cells and
                  // summing exactly to the parent face over its quadrants).
                  std::array<std::size_t, 3> face = cell;
                  if (position == 1) {
                    face[d] += 1U;
                  }
                  field_value = parent[d].value(face) + parent[d].slope(face, t1) * y1 +
                      parent[d].slope(face, t2) * y2;
                }
                std::array<std::size_t, 3> c{};
                c[d] = 2U * (cell[d] - block[d] * half[d]) + static_cast<std::size_t>(position + 1);
                c[t1] = 2U * (cell[t1] - block[t1] * half[t1]) + static_cast<std::size_t>(q1);
                c[t2] = 2U * (cell[t2] - block[t2] * half[t2]) + static_cast<std::size_t>(q2);
                values[faceIndex(child.geometry(), d, c[0], c[1], c[2])] = field_value * fine_area;
              }
            }
          }
        }
      }
    }
  }
  for (std::size_t block = 0; block < 8U; ++block) {
    const MhdFaceFluxStorage& child = staged[block];
    const double normalizer = magneticDivergenceNormalizer(child.geometry(), child.fluxX(),
                                                           child.fluxY(), child.fluxZ());
    result.max_cell_divergence_residual = std::max(
        result.max_cell_divergence_residual,
        maxNormalizedFaceDivergence(child.geometry(), child.fluxX(), child.fluxY(), child.fluxZ(),
                                    normalizer));
  }
  result.ok = result.max_cell_divergence_residual <= divergence_tolerance;
  if (!result.ok) {
    std::ostringstream stream;
    stream << "ratio-two prolongation exceeded the divergence tolerance (residual="
           << result.max_cell_divergence_residual << ", tolerance=" << divergence_tolerance << ')';
    result.detail = stream.str();
    return result;
  }
  for (std::size_t block = 0; block < 8U; ++block) {
    children[block] = std::move(staged[block]);
  }
  return result;
}

MhdAmrSynchronizer::TransferResult MhdAmrSynchronizer::prolongPatch(
    const MhdPatchGeometry& parent_geometry,
    const MhdGasStateSoa& parent_gas,
    const MhdFaceFluxView& parent_faces,
    double scale_factor,
    double gamma_adiabatic,
    std::array<MhdFaceFluxStorage, 8>& children_faces,
    std::array<MhdGasStateSoa, 8>& children_gas,
    double divergence_tolerance) {
  TransferResult result{};
  if (!parent_faces.valid() || parent_gas.size() != parent_geometry.cellCount() ||
      !(scale_factor > 0.0) || !(gamma_adiabatic > 1.0)) {
    result.detail = "invalid prolongation inputs";
    return result;
  }
  const std::array<std::size_t, 3> n{parent_geometry.nx, parent_geometry.ny, parent_geometry.nz};
  const std::array<std::size_t, 3> half{n[0] / 2U, n[1] / 2U, n[2] / 2U};
  for (int attempt = 0; attempt < 2; ++attempt) {
    const bool zero_slopes = attempt == 1;
    std::array<MhdFaceFluxStorage, 8> faces{};
    TransferResult magnetic =
        prolongFaceFluxes(parent_geometry, parent_faces.phi_x, parent_faces.phi_y,
                          parent_faces.phi_z, faces, divergence_tolerance, zero_slopes);
    if (!magnetic.ok) {
      result = magnetic;
      return result;  // a divergence violation is not curable by lower order
    }
    std::array<MhdGasStateSoa, 8> gas{};
    for (std::size_t block = 0; block < 8U; ++block) {
      gas[block].resize(faces[block].geometry().cellCount());
    }
    bool admissible = true;
    double energy_change = 0.0;
    for (std::size_t k = 0; k < n[2] && admissible; ++k) {
      for (std::size_t j = 0; j < n[1] && admissible; ++j) {
        for (std::size_t i = 0; i < n[0] && admissible; ++i) {
          const std::size_t parent_cell = parent_geometry.cellIndex(i, j, k);
          const MhdConservedState conserved = parent_gas.loadGas(parent_cell);
          const std::array<double, 3> b_parent = cellField(
              parent_geometry, parent_faces.phi_x, parent_faces.phi_y, parent_faces.phi_z, i, j, k);
          const double magnetic_parent =
              0.5 * (b_parent[0] * b_parent[0] + b_parent[1] * b_parent[1] +
                     b_parent[2] * b_parent[2]) /
              scale_factor;
          const double thermal_parent = thermalEnergy(conserved, b_parent, scale_factor);
          const std::array<std::size_t, 3> cell{i, j, k};
          std::array<std::size_t, 3> block{};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            block[axis] = cell[axis] / half[axis];
          }
          const std::size_t child_index = block[0] + 2U * (block[1] + 2U * block[2]);
          const MhdFaceFluxStorage& child = faces[child_index];
          const MhdPatchGeometry& child_geometry = child.geometry();
          std::array<double, 8> magnetic_fine{};
          double magnetic_mean = 0.0;
          for (std::size_t q = 0; q < 8U; ++q) {
            const std::size_t ci = 2U * (i - block[0] * half[0]) + (q & 1U);
            const std::size_t cj = 2U * (j - block[1] * half[1]) + ((q >> 1U) & 1U);
            const std::size_t ck = 2U * (k - block[2] * half[2]) + ((q >> 2U) & 1U);
            const std::array<double, 3> b = cellField(child_geometry, child.fluxX(), child.fluxY(),
                                                      child.fluxZ(), ci, cj, ck);
            magnetic_fine[q] = 0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / scale_factor;
            magnetic_mean += 0.125 * magnetic_fine[q];
          }
          const double delta = magnetic_mean - magnetic_parent;
          const double thermal_fine = thermal_parent - delta;
          if (!(thermal_parent > 0.0) || !(thermal_fine > 0.0) || !std::isfinite(thermal_fine)) {
            admissible = false;
            break;
          }
          energy_change += delta * parent_geometry.cellVolume();
          const double kinetic = conserved.total_energy_density_comoving - thermal_parent -
              magnetic_parent;
          for (std::size_t q = 0; q < 8U; ++q) {
            const std::size_t ci = 2U * (i - block[0] * half[0]) + (q & 1U);
            const std::size_t cj = 2U * (j - block[1] * half[1]) + ((q >> 1U) & 1U);
            const std::size_t ck = 2U * (k - block[2] * half[2]) + ((q >> 2U) & 1U);
            MhdConservedState fine = conserved;
            fine.total_energy_density_comoving = thermal_fine + kinetic + magnetic_fine[q];
            gas[child_index].storeGas(child_geometry.cellIndex(ci, cj, ck), fine);
          }
        }
      }
    }
    if (!admissible) {
      if (!zero_slopes) {
        continue;  // lower-order magnetic split
      }
      result.detail =
          "prolongation produced a non-positive fine thermal energy even with the zero-slope "
          "magnetic split; refinement rejected before commit";
      return result;
    }
    for (std::size_t block = 0; block < 8U; ++block) {
      children_faces[block] = std::move(faces[block]);
      children_gas[block] = std::move(gas[block]);
    }
    result = magnetic;
    result.magnetic_energy_change = energy_change;
    result.ok = true;
    return result;
  }
  return result;
}

MhdAmrSynchronizer::TransferResult MhdAmrSynchronizer::restrictPatch(
    const MhdPatchGeometry& parent_geometry,
    const std::array<MhdFaceFluxStorage, 8>& children_faces,
    const std::array<MhdGasStateSoa, 8>& children_gas,
    double scale_factor,
    double gamma_adiabatic,
    MhdGasStateSoa& parent_gas,
    std::span<double> parent_phi_x,
    std::span<double> parent_phi_y,
    std::span<double> parent_phi_z) {
  TransferResult result{};
  if (!(scale_factor > 0.0) || !(gamma_adiabatic > 1.0)) {
    result.detail = "invalid restriction inputs";
    return result;
  }
  MhdFaceFluxStorage faces(parent_geometry);
  restrictFaceFluxes(parent_geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ(), children_faces);
  const std::array<std::size_t, 3> n{parent_geometry.nx, parent_geometry.ny, parent_geometry.nz};
  const std::array<std::size_t, 3> half{n[0] / 2U, n[1] / 2U, n[2] / 2U};
  MhdGasStateSoa gas(parent_geometry.cellCount());
  for (std::size_t block = 0; block < 8U; ++block) {
    if (children_gas[block].size() != children_faces[block].geometry().cellCount()) {
      result.detail = "child gas state does not match its geometry";
      return result;
    }
  }
  for (std::size_t k = 0; k < n[2]; ++k) {
    for (std::size_t j = 0; j < n[1]; ++j) {
      for (std::size_t i = 0; i < n[0]; ++i) {
        const std::array<std::size_t, 3> cell{i, j, k};
        std::array<std::size_t, 3> block{};
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          block[axis] = cell[axis] / half[axis];
        }
        const std::size_t child_index = block[0] + 2U * (block[1] + 2U * block[2]);
        const MhdPatchGeometry& child_geometry = children_faces[child_index].geometry();
        std::array<double, MhdGasStateSoa::k_lane_count> sum{};
        for (std::size_t q = 0; q < 8U; ++q) {
          const std::size_t ci = 2U * (i - block[0] * half[0]) + (q & 1U);
          const std::size_t cj = 2U * (j - block[1] * half[1]) + ((q >> 1U) & 1U);
          const std::size_t ck = 2U * (k - block[2] * half[2]) + ((q >> 2U) & 1U);
          const std::size_t fine_cell = child_geometry.cellIndex(ci, cj, ck);
          for (std::size_t lane_index = 0; lane_index < sum.size(); ++lane_index) {
            sum[lane_index] += children_gas[child_index].lane(lane_index)[fine_cell];
          }
        }
        const std::size_t parent_cell = parent_geometry.cellIndex(i, j, k);
        for (std::size_t lane_index = 0; lane_index < sum.size(); ++lane_index) {
          gas.lane(lane_index)[parent_cell] = 0.125 * sum[lane_index];
        }
        const std::array<double, 3> b =
            cellField(parent_geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ(), i, j, k);
        if (!(thermalEnergy(gas.loadGas(parent_cell), b, scale_factor) > 0.0)) {
          result.detail = "restriction produced a non-positive parent thermal energy; "
                          "coarsening rejected before commit";
          return result;
        }
      }
    }
  }
  parent_gas.copyFrom(gas);
  std::copy(faces.fluxX().begin(), faces.fluxX().end(), parent_phi_x.begin());
  std::copy(faces.fluxY().begin(), faces.fluxY().end(), parent_phi_y.begin());
  std::copy(faces.fluxZ().begin(), faces.fluxZ().end(), parent_phi_z.begin());
  result.ok = true;
  return result;
}

void MhdAmrSynchronizer::applyRefluxCurl(const MhdPatchGeometry& coarse_geometry,
                                         std::span<double> coarse_phi_x,
                                         std::span<double> coarse_phi_y,
                                         std::span<double> coarse_phi_z,
                                         const MhdEmfRegister& reflux_register) {
  if (coarse_phi_x.size() != coarse_geometry.faceCountX() ||
      coarse_phi_y.size() != coarse_geometry.faceCountY() ||
      coarse_phi_z.size() != coarse_geometry.faceCountZ()) {
    throwRegister("coarse face lanes do not match the coarse patch");
  }
  reflux_register.finalize();
  MhdEdgeEmfStorage correction(coarse_geometry);
  for (const MhdEmfRegister::Record& record : reflux_register.records()) {
    const MhdCoarseEdgeKey& key = record.key;
    if (key.i < 0 || key.j < 0 || key.k < 0) {
      throwRegister("register edge index outside the coarse patch");
    }
    const auto i = static_cast<std::size_t>(key.i);
    const auto j = static_cast<std::size_t>(key.j);
    const auto k = static_cast<std::size_t>(key.k);
    const double delta = record.integrated_fine - record.integrated_coarse;
    switch (key.axis) {
      case MhdAxis::kX:
        if (i >= coarse_geometry.nx || j > coarse_geometry.ny || k > coarse_geometry.nz) {
          throwRegister("register edge index outside the coarse patch");
        }
        correction.emfX()[coarse_geometry.edgeIndexX(i, j, k)] = delta;
        break;
      case MhdAxis::kY:
        if (i > coarse_geometry.nx || j >= coarse_geometry.ny || k > coarse_geometry.nz) {
          throwRegister("register edge index outside the coarse patch");
        }
        correction.emfY()[coarse_geometry.edgeIndexY(i, j, k)] = delta;
        break;
      case MhdAxis::kZ:
        if (i > coarse_geometry.nx || j > coarse_geometry.ny || k >= coarse_geometry.nz) {
          throwRegister("register edge index outside the coarse patch");
        }
        correction.emfZ()[coarse_geometry.edgeIndexZ(i, j, k)] = delta;
        break;
    }
  }
  applyConstrainedTransportUpdate(coarse_geometry, coarse_phi_x, coarse_phi_y, coarse_phi_z,
                                  correction.emfX(), correction.emfY(), correction.emfZ(), 1.0);
}

double MhdAmrSynchronizer::magneticEnergySum(const MhdPatchGeometry& geometry,
                                             std::span<const double> phi_x,
                                             std::span<const double> phi_y,
                                             std::span<const double> phi_z,
                                             double scale_factor) {
  if (phi_x.size() != geometry.faceCountX() || phi_y.size() != geometry.faceCountY() ||
      phi_z.size() != geometry.faceCountZ() || !(scale_factor > 0.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_amr", "magnetic energy request with invalid inputs"));
  }
  const double volume = geometry.cellVolume();
  double total = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::array<double, 3> b = cellField(geometry, phi_x, phi_y, phi_z, i, j, k);
        total += 0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / scale_factor * volume;
      }
    }
  }
  return total;
}

std::array<double, 3> MhdAmrSynchronizer::netBoundaryFlux(const MhdPatchGeometry& geometry,
                                                          std::span<const double> phi_x,
                                                          std::span<const double> phi_y,
                                                          std::span<const double> phi_z) {
  if (phi_x.size() != geometry.faceCountX() || phi_y.size() != geometry.faceCountY() ||
      phi_z.size() != geometry.faceCountZ()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_amr", "net flux request with mis-sized lanes"));
  }
  std::array<double, 3> net{};
  double upper = 0.0;
  double lower = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      lower += phi_x[geometry.faceIndexX(0, j, k)];
      upper += phi_x[geometry.faceIndexX(geometry.nx, j, k)];
    }
  }
  net[0] = upper - lower;
  upper = lower = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t i = 0; i < geometry.nx; ++i) {
      lower += phi_y[geometry.faceIndexY(i, 0, k)];
      upper += phi_y[geometry.faceIndexY(i, geometry.ny, k)];
    }
  }
  net[1] = upper - lower;
  upper = lower = 0.0;
  for (std::size_t j = 0; j < geometry.ny; ++j) {
    for (std::size_t i = 0; i < geometry.nx; ++i) {
      lower += phi_z[geometry.faceIndexZ(i, j, 0)];
      upper += phi_z[geometry.faceIndexZ(i, j, geometry.nz)];
    }
  }
  net[2] = upper - lower;
  return net;
}

}  // namespace cosmosim::mhd
