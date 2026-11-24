#include "cosmosim/mhd/mhd_state.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <new>
#include <utility>

namespace cosmosim::mhd {
namespace {

// Rounded-up byte size for `count` doubles in a 64-byte aligned allocation.
std::size_t alignedBytesFor(std::size_t count) noexcept {
  if (count == 0U) {
    return 0U;
  }
  constexpr std::size_t k_alignment = static_cast<std::size_t>(k_mhd_face_alignment_bytes);
  if (count > (std::numeric_limits<std::size_t>::max() / sizeof(double))) {
    throw std::bad_array_new_length();
  }
  const std::size_t bytes = count * sizeof(double);
  const std::size_t remainder = bytes % k_alignment;
  return remainder == 0U ? bytes : bytes + (k_alignment - remainder);
}

}  // namespace

void MhdDoubleArray::AlignedDeleter::operator()(double* pointer) const noexcept {
  if (pointer == nullptr) {
    return;
  }
  ::operator delete(static_cast<void*>(pointer),
                    std::align_val_t{static_cast<std::size_t>(k_mhd_face_alignment_bytes)});
}

std::size_t MhdDoubleArray::allocatorSizeBytes(std::size_t size) noexcept {
  return alignedBytesFor(size);
}

double* MhdDoubleArray::allocateAligned(std::size_t size) {
  if (size == 0U) {
    return nullptr;
  }
  const std::size_t bytes = alignedBytesFor(size);
  void* raw = ::operator new(bytes,
                             std::align_val_t{static_cast<std::size_t>(k_mhd_face_alignment_bytes)});
  return static_cast<double*>(raw);
}

MhdDoubleArray::MhdDoubleArray(std::size_t size, double fill_value) {
  resize(size, fill_value);
}

MhdDoubleArray::MhdDoubleArray(const MhdDoubleArray& other) {
  if (other.m_size != 0U) {
    m_data = std::unique_ptr<double, AlignedDeleter>(allocateAligned(other.m_size));
    m_size = other.m_size;
    m_capacity = static_cast<std::uint64_t>(allocatorSizeBytes(other.m_size));
    std::memcpy(m_data.get(), other.m_data.get(), m_size * sizeof(double));
  }
}

MhdDoubleArray& MhdDoubleArray::operator=(const MhdDoubleArray& other) {
  if (this == &other) {
    return *this;
  }
  MhdDoubleArray copy(other);
  swap(copy);
  return *this;
}

MhdDoubleArray::MhdDoubleArray(MhdDoubleArray&& other) noexcept
    : m_data(std::move(other.m_data)),
      m_size(other.m_size),
      m_capacity(other.m_capacity) {
  other.m_size = 0U;
  other.m_capacity = 0U;
}

MhdDoubleArray& MhdDoubleArray::operator=(MhdDoubleArray&& other) noexcept {
  if (this == &other) {
    return *this;
  }
  m_data = std::move(other.m_data);
  m_size = other.m_size;
  m_capacity = other.m_capacity;
  other.m_size = 0U;
  other.m_capacity = 0U;
  return *this;
}

void MhdDoubleArray::resize(std::size_t size, double fill_value) {
  if (size == m_size) {
    if (fill_value != 0.0) {
      std::fill_n(m_data.get(), m_size, fill_value);
    }
    return;
  }
  MhdDoubleArray replacement(size, fill_value);
  swap(replacement);
}

void MhdDoubleArray::assignFill(double value) noexcept {
  std::fill_n(m_data.get(), m_size, value);
}

void MhdDoubleArray::swap(MhdDoubleArray& other) noexcept {
  m_data.swap(other.m_data);
  std::swap(m_size, other.m_size);
  std::swap(m_capacity, other.m_capacity);
}

// ---------------------------------------------------------------------------
// MhdFaceFluxStorage
// ---------------------------------------------------------------------------

MhdFaceFluxStorage::MhdFaceFluxStorage(const MhdPatchGeometry& geometry) {
  resize(geometry);
}

void MhdFaceFluxStorage::resize(const MhdPatchGeometry& geometry) {
  m_geometry = geometry;
  m_phi_x.resize(geometry.faceCountX(), 0.0);
  m_phi_y.resize(geometry.faceCountY(), 0.0);
  m_phi_z.resize(geometry.faceCountZ(), 0.0);
}

std::array<double, 3> MhdFaceFluxStorage::cellFieldComoving(
    std::size_t i, std::size_t j, std::size_t k) const noexcept {
  const double area_x = m_geometry.faceAreaX();
  const double area_y = m_geometry.faceAreaY();
  const double area_z = m_geometry.faceAreaZ();
  const double bx = area_x > 0.0
      ? 0.5 * (faceFluxX(i, j, k) + faceFluxX(i + 1U, j, k)) / area_x
      : 0.0;
  const double by = area_y > 0.0
      ? 0.5 * (faceFluxY(i, j, k) + faceFluxY(i, j + 1U, k)) / area_y
      : 0.0;
  const double bz = area_z > 0.0
      ? 0.5 * (faceFluxZ(i, j, k) + faceFluxZ(i, j, k + 1U)) / area_z
      : 0.0;
  return {bx, by, bz};
}

bool MhdFaceFluxView::valid() const noexcept {
  if (geometry == nullptr || !geometry->valid()) {
    return false;
  }
  return phi_x.size() == geometry->faceCountX() &&
      phi_y.size() == geometry->faceCountY() &&
      phi_z.size() == geometry->faceCountZ();
}

bool MhdFaceFluxMutableView::valid() const noexcept {
  if (geometry == nullptr || !geometry->valid()) {
    return false;
  }
  return phi_x.size() == geometry->faceCountX() &&
      phi_y.size() == geometry->faceCountY() &&
      phi_z.size() == geometry->faceCountZ();
}

MhdFaceFluxView makeFaceFluxView(const MhdFaceFluxStorage& storage,
                                 std::uint64_t patch_id,
                                 std::uint64_t geometry_generation,
                                 int refinement_level) noexcept {
  return MhdFaceFluxView{
      .geometry = &storage.geometry(),
      .phi_x = storage.fluxX(),
      .phi_y = storage.fluxY(),
      .phi_z = storage.fluxZ(),
      .patch_id = patch_id,
      .geometry_generation = geometry_generation,
      .refinement_level = refinement_level,
  };
}

MhdFaceFluxMutableView makeMutableFaceFluxView(MhdFaceFluxStorage& storage,
                                               std::uint64_t patch_id,
                                               std::uint64_t geometry_generation,
                                               int refinement_level) noexcept {
  return MhdFaceFluxMutableView{
      .geometry = &storage.geometry(),
      .phi_x = storage.fluxX(),
      .phi_y = storage.fluxY(),
      .phi_z = storage.fluxZ(),
      .patch_id = patch_id,
      .geometry_generation = geometry_generation,
      .refinement_level = refinement_level,
  };
}

// ---------------------------------------------------------------------------
// Exact storage budget (CHUI-MHD-PERF v1.0 section 4)
// ---------------------------------------------------------------------------

MhdPatchMemoryEstimate estimatePatchMemory(const MhdPatchGeometry& geometry,
                                           std::size_t gas_lane_count,
                                           bool nonideal_enabled,
                                           std::int32_t rkl_max_stages) {
  MhdPatchMemoryEstimate estimate{};
  const auto doubles_to_bytes = [](std::size_t count) {
    return static_cast<std::uint64_t>(alignedBytesFor(count));
  };

  const std::size_t cells = geometry.cellCount();
  const std::size_t faces = geometry.faceCountTotal();
  const std::size_t edges = geometry.edgeCountTotal();

  estimate.face_flux_capacity_bytes = doubles_to_bytes(faces);
  estimate.edge_emf_capacity_bytes = doubles_to_bytes(edges);
  // SSPRK2 trial gas copy: one full conserved-lane copy.
  estimate.gas_trial_capacity_bytes = doubles_to_bytes(cells * gas_lane_count);
  estimate.edge_current_capacity_bytes =
      nonideal_enabled ? doubles_to_bytes(cells * 3U) : 0U;

  std::uint64_t rkl_extra = 0U;
  if (nonideal_enabled) {
    const std::int32_t stages = rkl_max_stages > 1 ? rkl_max_stages : 1;
    // RKL2 reuses a single stage edge buffer plus a single stage energy lane,
    // so residency does not grow with the stage count. The bound records the
    // per-cycle reusable footprint; the stage count only affects the cost.
    (void)stages;
    rkl_extra = doubles_to_bytes(edges) + doubles_to_bytes(cells);
  }
  estimate.rkl_cycle_extra_capacity_bytes = rkl_extra;

  // Phase scratch excludes the canonical persistent face flux.
  estimate.total_phase_scratch_bytes = estimate.edge_emf_capacity_bytes +
      estimate.gas_trial_capacity_bytes + estimate.edge_current_capacity_bytes +
      estimate.rkl_cycle_extra_capacity_bytes;

  estimate.total_canonical_persistent_bytes = estimate.face_flux_capacity_bytes;
  return estimate;
}

}  // namespace cosmosim::mhd
