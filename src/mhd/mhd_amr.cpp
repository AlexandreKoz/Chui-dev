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

// Coarse/child face correspondence for ratio two. A child block `di` in one
// dimension spans child faces [di*(nx/2)*2, ...] so child face `f` maps to
// coarse face `di*(nx/2) + floor(f/2)`; odd `f` is a half-cell face that sits
// strictly inside the coarse cell `di*(nx/2) + (f-1)/2`.
std::size_t coarseFaceForChild(std::size_t block, std::size_t child_face,
                               std::size_t coarse_extent) {
  return block * (coarse_extent / 2U) + child_face / 2U;
}

double parentFaceInterpolated(std::span<const double> parent_faces,
                              const MhdPatchGeometry& parent,
                              MhdAxis axis,
                              std::size_t coarse_extent,
                              std::size_t child_face,
                              std::size_t block,
                              std::size_t other_a,
                              std::size_t other_b) {
  const std::size_t coarse_face = coarseFaceForChild(block, child_face, coarse_extent);
  const bool half_face = (child_face % 2U) == 1U;
  const auto face_at = [&](std::size_t index) {
    if (axis == MhdAxis::kX) {
      return parent_faces[parent.faceIndexX(index, other_a, other_b)];
    }
    if (axis == MhdAxis::kY) {
      return parent_faces[parent.faceIndexY(other_a, index, other_b)];
    }
    return parent_faces[parent.faceIndexZ(other_a, other_b, index)];
  };
  if (!half_face) {
    return face_at(coarse_face);
  }
  // Internal half face: opposite-parent-face mean, which keeps the child block
  // exactly divergence free for a divergence-free coarse cell.
  const std::size_t upper = std::min(coarse_face + 1U, coarse_extent);
  return 0.5 * (face_at(coarse_face) + face_at(upper));
}

// Map a child block-local normal face index to a coarse face in the *other*
// two axes: each child face pair shares one coarse face in each transverse
// direction.

}  // namespace

std::uint64_t mhdEmfRegisterBound(const MhdPatchGeometry& coarse_geometry) {
  // Conservative bound: every edge of the coarse patch shell can receive a
  // fine contribution. Computed before stepping so admission never depends on
  // runtime growth.
  return static_cast<std::uint64_t>(coarse_geometry.edgeCountTotal());
}

// ---------------------------------------------------------------------------
// MhdEmfRegister
// ---------------------------------------------------------------------------

MhdEmfRegister::MhdEmfRegister(std::uint64_t max_entries_per_interval)
    : m_max_entries(max_entries_per_interval) {}

bool MhdEmfRegister::addContribution(const MhdEmfRegisterRecord& record) {
  if (m_max_entries != 0U && m_records.size() >= m_max_entries) {
    return false;
  }
  m_records.push_back(record);
  return true;
}

void MhdEmfRegister::clear() noexcept {
  m_records.clear();
  m_has_interval = false;
}

std::uint64_t MhdEmfRegister::ownedCapacityBytes() const noexcept {
  return static_cast<std::uint64_t>(m_records.capacity() * sizeof(MhdEmfRegisterRecord));
}

void MhdEmfRegister::finalizeDeterministicOrder() {
  std::sort(m_records.begin(), m_records.end(),
            [](const MhdEmfRegisterRecord& lhs, const MhdEmfRegisterRecord& rhs) {
              if (lhs.edge_hash != rhs.edge_hash) {
                return lhs.edge_hash < rhs.edge_hash;
              }
              if (lhs.coarse_patch_id != rhs.coarse_patch_id) {
                return lhs.coarse_patch_id < rhs.coarse_patch_id;
              }
              if (lhs.source_rank != rhs.source_rank) {
                return lhs.source_rank < rhs.source_rank;
              }
              return lhs.substep_ordinal < rhs.substep_ordinal;
            });
  for (std::size_t index = 1U; index < m_records.size(); ++index) {
    const MhdEmfRegisterRecord& previous = m_records[index - 1U];
    const MhdEmfRegisterRecord& current = m_records[index];
    if (previous.edge_hash == current.edge_hash &&
        previous.coarse_patch_id == current.coarse_patch_id &&
        previous.source_rank == current.source_rank &&
        previous.substep_ordinal == current.substep_ordinal) {
      throw MhdError(MhdErrorCode::kAmrRegisterIncomplete,
                     makeContext("mhd_emf_register",
                                 "duplicate (edge, patch, rank, substep) register entry"));
    }
  }
}

void MhdEmfRegister::bindInterval(const MhdAmrIntervalId& interval) {
  m_interval = interval;
  m_has_interval = true;
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
  if (parent.nx % 2U != 0U || parent.ny % 2U != 0U || parent.nz % 2U != 0U) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_amr", "ratio-two prolongation requires even patch extents"));
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
  child.origin_x_comoving =
      parent.origin_x_comoving + static_cast<double>(di) * 0.5 * parent.nx * parent.dx_comoving;
  child.origin_y_comoving =
      parent.origin_y_comoving + static_cast<double>(dj) * 0.5 * parent.ny * parent.dy_comoving;
  child.origin_z_comoving =
      parent.origin_z_comoving + static_cast<double>(dk) * 0.5 * parent.nz * parent.dz_comoving;
  return child;
}

void MhdAmrSynchronizer::restrictFaceFluxes(
    const MhdPatchGeometry& parent_geometry,
    std::span<double> parent_phi_x,
    std::span<double> parent_phi_y,
    std::span<double> parent_phi_z,
    const std::array<MhdFaceFluxStorage, 8>& children) {
  if (parent_geometry.nx % 2U != 0U || parent_geometry.ny % 2U != 0U ||
      parent_geometry.nz % 2U != 0U) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_amr", "ratio-two restriction requires even patch extents"));
  }
  const std::size_t nx = parent_geometry.nx;
  const std::size_t ny = parent_geometry.ny;
  const std::size_t nz = parent_geometry.nz;
  const std::size_t half_x = nx / 2U;
  const std::size_t half_y = ny / 2U;
  const std::size_t half_z = nz / 2U;
  std::fill(parent_phi_x.begin(), parent_phi_x.end(), 0.0);
  std::fill(parent_phi_y.begin(), parent_phi_y.end(), 0.0);
  std::fill(parent_phi_z.begin(), parent_phi_z.end(), 0.0);

  for (std::int32_t dk = 0; dk < 2; ++dk) {
    for (std::int32_t dj = 0; dj < 2; ++dj) {
      for (std::int32_t di = 0; di < 2; ++di) {
        const std::int32_t block = di + 2 * (dj + 2 * dk);
        const MhdFaceFluxStorage& child = children[static_cast<std::size_t>(block)];
        const std::size_t coarse_offset_x = static_cast<std::size_t>(di) * half_x;
        const std::size_t coarse_offset_y = static_cast<std::size_t>(dj) * half_y;
        const std::size_t coarse_offset_z = static_cast<std::size_t>(dk) * half_z;

        // x faces: one block in x, four sub-faces in (y,z).
        for (std::size_t kc = 0; kc < half_z; ++kc) {
          for (std::size_t jc = 0; jc < half_y; ++jc) {
            for (std::size_t fc = 0; fc <= half_x; ++fc) {
              const std::size_t coarse_face = coarse_offset_x + fc;
              double sum = 0.0;
              for (std::size_t sub_k = 0; sub_k < 2U; ++sub_k) {
                for (std::size_t sub_j = 0; sub_j < 2U; ++sub_j) {
                  const std::size_t child_i = 2U * fc;
                  const std::size_t child_j = 2U * jc + sub_j;
                  const std::size_t child_k = 2U * kc + sub_k;
                  sum += child.faceFluxX(child_i, child_j, child_k);
                }
              }
              parent_phi_x[parent_geometry.faceIndexX(coarse_face,
                                                      coarse_offset_y + jc,
                                                      coarse_offset_z + kc)] += sum;
            }
          }
        }

        // y faces: one block in y, four sub-faces in (x,z).
        for (std::size_t kc = 0; kc < half_z; ++kc) {
          for (std::size_t jc = 0; jc <= half_y; ++jc) {
            for (std::size_t ic = 0; ic < half_x; ++ic) {
              const std::size_t coarse_face = coarse_offset_y + jc;
              double sum = 0.0;
              for (std::size_t sub_k = 0; sub_k < 2U; ++sub_k) {
                for (std::size_t sub_i = 0; sub_i < 2U; ++sub_i) {
                  const std::size_t child_i = 2U * ic + sub_i;
                  const std::size_t child_j = 2U * jc;
                  const std::size_t child_k = 2U * kc + sub_k;
                  sum += child.faceFluxY(child_i, child_j, child_k);
                }
              }
              parent_phi_y[parent_geometry.faceIndexY(coarse_offset_x + ic,
                                                      coarse_face,
                                                      coarse_offset_z + kc)] += sum;
            }
          }
        }

        // z faces: one block in z, four sub-faces in (x,y).
        for (std::size_t kc = 0; kc <= half_z; ++kc) {
          for (std::size_t jc = 0; jc < half_y; ++jc) {
            for (std::size_t ic = 0; ic < half_x; ++ic) {
              const std::size_t coarse_face = coarse_offset_z + kc;
              double sum = 0.0;
              for (std::size_t sub_j = 0; sub_j < 2U; ++sub_j) {
                for (std::size_t sub_i = 0; sub_i < 2U; ++sub_i) {
                  const std::size_t child_i = 2U * ic + sub_i;
                  const std::size_t child_j = 2U * jc + sub_j;
                  const std::size_t child_k = 2U * kc;
                  sum += child.faceFluxZ(child_i, child_j, child_k);
                }
              }
              parent_phi_z[parent_geometry.faceIndexZ(coarse_offset_x + ic,
                                                      coarse_offset_y + jc,
                                                      coarse_face)] += sum;
            }
          }
        }
      }
    }
  }
}

MhdAmrSynchronizer::TransferResult MhdAmrSynchronizer::prolongFaceFluxes(
    const MhdPatchGeometry& parent_geometry,
    std::span<const double> parent_phi_x,
    std::span<const double> parent_phi_y,
    std::span<const double> parent_phi_z,
    std::array<MhdFaceFluxStorage, 8>& children,
    double divergence_tolerance) {
  TransferResult result{};
  if (parent_geometry.nx % 2U != 0U || parent_geometry.ny % 2U != 0U ||
      parent_geometry.nz % 2U != 0U) {
    result.detail = "ratio-two prolongation requires even patch extents";
    return result;
  }
  if (parent_phi_x.size() != parent_geometry.faceCountX() ||
      parent_phi_y.size() != parent_geometry.faceCountY() ||
      parent_phi_z.size() != parent_geometry.faceCountZ()) {
    result.detail = "parent face lanes do not match the parent geometry";
    return result;
  }
  const std::size_t nx = parent_geometry.nx;
  const std::size_t ny = parent_geometry.ny;
  const std::size_t nz = parent_geometry.nz;
  const std::size_t half_x = nx / 2U;
  const std::size_t half_y = ny / 2U;
  const std::size_t half_z = nz / 2U;

  for (std::int32_t block = 0; block < 8; ++block) {
    const MhdPatchGeometry child = childGeometry(parent_geometry, block);
    const std::int32_t di = block % 2;
    const std::int32_t dj = (block / 2) % 2;
    const std::int32_t dk = block / 4;
    MhdFaceFluxStorage& target = children[static_cast<std::size_t>(block)];
    if (target.geometry().nx != child.nx || target.geometry().ny != child.ny ||
        target.geometry().nz != child.nz) {
      target.resize(child);
    }
    // x faces
    for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t j = 0; j < ny; ++j) {
        const std::size_t coarse_j = static_cast<std::size_t>(dj) * half_y + j / 2U;
        const std::size_t coarse_k = static_cast<std::size_t>(dk) * half_z + k / 2U;
        for (std::size_t f = 0; f <= nx; ++f) {
          const double value = parentFaceInterpolated(
              parent_phi_x, parent_geometry, MhdAxis::kX, nx, f,
              static_cast<std::size_t>(di), coarse_j, coarse_k);
          target.setFaceFluxX(f, j, k, value);
        }
      }
    }
    // y faces
    for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t j = 0; j <= ny; ++j) {
        const std::size_t coarse_k = static_cast<std::size_t>(dk) * half_z + k / 2U;
        for (std::size_t i = 0; i < nx; ++i) {
          const std::size_t coarse_i = static_cast<std::size_t>(di) * half_x + i / 2U;
          const double value = parentFaceInterpolated(
              parent_phi_y, parent_geometry, MhdAxis::kY, ny, j,
              static_cast<std::size_t>(dj), coarse_i, coarse_k);
          target.setFaceFluxY(i, j, k, value);
        }
      }
    }
    // z faces
    for (std::size_t k = 0; k <= nz; ++k) {
      for (std::size_t j = 0; j < ny; ++j) {
        const std::size_t coarse_j = static_cast<std::size_t>(dj) * half_y + j / 2U;
        for (std::size_t i = 0; i < nx; ++i) {
          const std::size_t coarse_i = static_cast<std::size_t>(di) * half_x + i / 2U;
          const double value = parentFaceInterpolated(
              parent_phi_z, parent_geometry, MhdAxis::kZ, nz, k,
              static_cast<std::size_t>(dk), coarse_i, coarse_j);
          target.setFaceFluxZ(i, j, k, value);
        }
      }
    }
    // The interior half faces must resolve to the same value on both sides of a
    // block boundary; `parentFaceInterpolated` is a pure function of the coarse
    // face, so the construction is already consistent.
    const double residual = maxNormalizedFaceDivergence(
        child, target.fluxX(), target.fluxY(), target.fluxZ(),
        std::max(faceFluxMagnitudeSum(target.fluxX(), target.fluxY(), target.fluxZ()),
                 1.0));
    result.max_cell_divergence_residual =
        std::max(result.max_cell_divergence_residual, residual);
  }

  result.used_piecewise_constant = true;
  result.ok = result.max_cell_divergence_residual <= divergence_tolerance;
  if (!result.ok) {
    std::ostringstream stream;
    stream << "ratio-two prolongation exceeded the divergence tolerance (residual="
           << result.max_cell_divergence_residual << ", tolerance=" << divergence_tolerance
           << ')';
    result.detail = stream.str();
  }
  return result;
}

void MhdAmrSynchronizer::applyRefluxCurl(const MhdPatchGeometry& coarse_geometry,
                                         std::span<double> coarse_phi_x,
                                         std::span<double> coarse_phi_y,
                                         std::span<double> coarse_phi_z,
                                         const MhdEmfRegister& register,
                                         std::uint64_t coarse_patch_id) {
  if (coarse_phi_x.size() != coarse_geometry.faceCountX() ||
      coarse_phi_y.size() != coarse_geometry.faceCountY() ||
      coarse_phi_z.size() != coarse_geometry.faceCountZ()) {
    throw MhdError(MhdErrorCode::kAmrRegisterIncomplete,
                   makeContext("mhd_reflux", "coarse face lanes do not match the coarse patch"));
  }
  MhdEdgeEmfStorage flux_correction(coarse_geometry);
  bool any_record = false;
  for (const MhdEmfRegisterRecord& record : register.records()) {
    if (record.coarse_patch_id != coarse_patch_id) {
      continue;
    }
    any_record = true;
    std::span<double> target;
    std::size_t index = 0;
    switch (record.edge_axis) {
      case MhdAxis::kX:
        target = flux_correction.emfX();
        index = coarse_geometry.edgeIndexX(static_cast<std::size_t>(record.i),
                                           static_cast<std::size_t>(record.j),
                                           static_cast<std::size_t>(record.k));
        break;
      case MhdAxis::kY:
        target = flux_correction.emfY();
        index = coarse_geometry.edgeIndexY(static_cast<std::size_t>(record.i),
                                           static_cast<std::size_t>(record.j),
                                           static_cast<std::size_t>(record.k));
        break;
      case MhdAxis::kZ:
        target = flux_correction.emfZ();
        index = coarse_geometry.edgeIndexZ(static_cast<std::size_t>(record.i),
                                           static_cast<std::size_t>(record.j),
                                           static_cast<std::size_t>(record.k));
        break;
    }
    if (index >= target.size()) {
      throw MhdError(MhdErrorCode::kAmrRegisterIncomplete,
                     makeContext("mhd_reflux", "register edge index outside the coarse patch"));
    }
    target[index] += record.integrated_emf;
  }
  if (!any_record) {
    // No fine coverage for this patch is a legitimate no-op.
    return;
  }
  applyConstrainedTransportUpdate(coarse_geometry, coarse_phi_x, coarse_phi_y,
                                  coarse_phi_z, flux_correction.emfX(),
                                  flux_correction.emfY(), flux_correction.emfZ(), 1.0);
}

double MhdAmrSynchronizer::magneticEnergySum(const MhdPatchGeometry& geometry,
                                             std::span<const double> phi_x,
                                             std::span<const double> phi_y,
                                             std::span<const double> phi_z,
                                             double scale_factor) {
  if (phi_x.size() != geometry.faceCountX() ||
      phi_y.size() != geometry.faceCountY() ||
      phi_z.size() != geometry.faceCountZ()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double volume = geometry.cellVolume();
  const double factor = scale_factor > 0.0 ? 1.0 / (2.0 * scale_factor) : 0.5;
  double total = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const double bx = 0.5 * (phi_x[geometry.faceIndexX(i, j, k)] +
                                 phi_x[geometry.faceIndexX(i + 1U, j, k)]) /
            geometry.faceAreaX();
        const double by = 0.5 * (phi_y[geometry.faceIndexY(i, j, k)] +
                                 phi_y[geometry.faceIndexY(i, j + 1U, k)]) /
            geometry.faceAreaY();
        const double bz = 0.5 * (phi_z[geometry.faceIndexZ(i, j, k)] +
                                 phi_z[geometry.faceIndexZ(i, j, k + 1U)]) /
            geometry.faceAreaZ();
        total += factor * (bx * bx + by * by + bz * bz) * volume;
      }
    }
  }
  return total;
}

std::array<double, 3> MhdAmrSynchronizer::netBoundaryFlux(
    const MhdPatchGeometry& geometry,
    std::span<const double> phi_x,
    std::span<const double> phi_y,
    std::span<const double> phi_z) {
  std::array<double, 3> net{0.0, 0.0, 0.0};
  if (phi_x.size() != geometry.faceCountX() ||
      phi_y.size() != geometry.faceCountY() ||
      phi_z.size() != geometry.faceCountZ()) {
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
  }
  double upper = 0.0;
  double lower = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      lower += phi_x[geometry.faceIndexX(0, j, k)];
      upper += phi_x[geometry.faceIndexX(geometry.nx, j, k)];
    }
  }
  net[0] = upper - lower;
  upper = 0.0;
  lower = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t i = 0; i < geometry.nx; ++i) {
      lower += phi_y[geometry.faceIndexY(i, 0, k)];
      upper += phi_y[geometry.faceIndexY(i, geometry.ny, k)];
    }
  }
  net[1] = upper - lower;
  upper = 0.0;
  lower = 0.0;
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
