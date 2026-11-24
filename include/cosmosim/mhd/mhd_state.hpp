#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/mhd/mhd_validation.hpp"

namespace cosmosim::mhd {

// ---------------------------------------------------------------------------
// Canonical magnetic degrees of freedom.
//
// The authoritative persistent representation is the oriented, face-integrated
// comoving magnetic flux Phi_axis = B_axis_comov * A_face_comov in double
// precision. Shapes follow CHUI-MHD-REF v1.0 section 5.1:
//   Phi_x: (nx+1, ny,   nz  )
//   Phi_y: (nx,   ny+1, nz  )
//   Phi_z: (nx,   ny,   nz+1)
// Each array is axis-fastest: i + nx_axis*(j + ny_axis*k).
// A cell-centered field is derived and ephemeral; it is never an authority.
// ---------------------------------------------------------------------------

inline constexpr std::uint64_t k_mhd_face_alignment_bytes = 64U;
inline constexpr double k_mhd_gamma_default = 5.0 / 3.0;
// Regularization of |B| used only inside coefficients that must vanish at B=0
// (specification section 4.1). Never applied to the field itself.
inline constexpr double k_mhd_b_zero_regularization_fraction = 1.0e-30;

// 64-byte aligned owning buffer of doubles with explicit retained capacity.
class MhdDoubleArray final {
 public:
  struct AlignedDeleter {
    void operator()(double* pointer) const noexcept;
  };

  MhdDoubleArray() = default;
  explicit MhdDoubleArray(std::size_t size, double fill_value = 0.0);
  ~MhdDoubleArray() = default;

  MhdDoubleArray(const MhdDoubleArray& other);
  MhdDoubleArray& operator=(const MhdDoubleArray& other);
  MhdDoubleArray(MhdDoubleArray&& other) noexcept;
  MhdDoubleArray& operator=(MhdDoubleArray&& other) noexcept;

  void resize(std::size_t size, double fill_value = 0.0);
  void assignFill(double value) noexcept;
  void fill(double value) noexcept { assignFill(value); }
  [[nodiscard]] std::size_t size() const noexcept { return m_size; }
  [[nodiscard]] bool empty() const noexcept { return m_size == 0U; }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept { return m_capacity; }
  [[nodiscard]] std::span<double> span() noexcept { return {m_data.get(), m_size}; }
  [[nodiscard]] std::span<const double> span() const noexcept { return {m_data.get(), m_size}; }
  [[nodiscard]] double* data() noexcept { return m_data.get(); }
  [[nodiscard]] const double* data() const noexcept { return m_data.get(); }
  [[nodiscard]] double* begin() noexcept { return m_data.get(); }
  [[nodiscard]] double* end() noexcept { return m_data.get() + m_size; }
  [[nodiscard]] const double* begin() const noexcept { return m_data.get(); }
  [[nodiscard]] const double* end() const noexcept { return m_data.get() + m_size; }
  [[nodiscard]] const double* cbegin() const noexcept { return begin(); }
  [[nodiscard]] const double* cend() const noexcept { return end(); }
  [[nodiscard]] double& operator[](std::size_t index) noexcept { return m_data[index]; }
  [[nodiscard]] double operator[](std::size_t index) const noexcept { return m_data[index]; }
  void swap(MhdDoubleArray& other) noexcept;

 private:
  static std::size_t allocatorSizeBytes(std::size_t size) noexcept;
  static double* allocateAligned(std::size_t size);
  std::unique_ptr<double, AlignedDeleter> m_data{};
  std::size_t m_size = 0U;
  std::uint64_t m_capacity = 0U;
};

// Cartesian patch geometry for one MHD patch. Lengths are comoving code units.
struct MhdPatchGeometry {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  double origin_x_comoving = 0.0;
  double origin_y_comoving = 0.0;
  double origin_z_comoving = 0.0;
  double dx_comoving = 0.0;
  double dy_comoving = 0.0;
  double dz_comoving = 0.0;

  [[nodiscard]] std::size_t cellCount() const noexcept { return nx * ny * nz; }
  [[nodiscard]] std::size_t faceCountX() const noexcept { return (nx + 1U) * ny * nz; }
  [[nodiscard]] std::size_t faceCountY() const noexcept { return nx * (ny + 1U) * nz; }
  [[nodiscard]] std::size_t faceCountZ() const noexcept { return nx * ny * (nz + 1U); }
  [[nodiscard]] std::size_t faceCountTotal() const noexcept {
    return faceCountX() + faceCountY() + faceCountZ();
  }
  [[nodiscard]] std::size_t edgeCountX() const noexcept { return nx * (ny + 1U) * (nz + 1U); }
  [[nodiscard]] std::size_t edgeCountY() const noexcept { return (nx + 1U) * ny * (nz + 1U); }
  [[nodiscard]] std::size_t edgeCountZ() const noexcept { return (nx + 1U) * (ny + 1U) * nz; }
  [[nodiscard]] std::size_t edgeCountTotal() const noexcept {
    return edgeCountX() + edgeCountY() + edgeCountZ();
  }
  [[nodiscard]] double faceAreaX() const noexcept { return dy_comoving * dz_comoving; }
  [[nodiscard]] double faceAreaY() const noexcept { return dx_comoving * dz_comoving; }
  [[nodiscard]] double faceAreaZ() const noexcept { return dx_comoving * dy_comoving; }
  [[nodiscard]] double cellVolume() const noexcept {
    return dx_comoving * dy_comoving * dz_comoving;
  }
  [[nodiscard]] std::size_t cellIndex(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return i + nx * (j + ny * k);
  }
  [[nodiscard]] std::size_t faceIndexX(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return i + (nx + 1U) * (j + ny * k);
  }
  [[nodiscard]] std::size_t faceIndexY(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return i + nx * (j + (ny + 1U) * k);
  }
  [[nodiscard]] std::size_t faceIndexZ(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return i + nx * (j + ny * k);
  }
  [[nodiscard]] std::size_t edgeIndexX(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return i + nx * (j + (ny + 1U) * k);
  }
  [[nodiscard]] std::size_t edgeIndexY(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return i + (nx + 1U) * (j + ny * k);
  }
  [[nodiscard]] std::size_t edgeIndexZ(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return i + (nx + 1U) * (j + (ny + 1U) * k);
  }
  [[nodiscard]] bool valid() const noexcept {
    return nx > 0U && ny > 0U && nz > 0U && dx_comoving > 0.0 && dy_comoving > 0.0 &&
        dz_comoving > 0.0;
  }
};

// Owning patch-local face flux authority.
class MhdFaceFluxStorage final {
 public:
  MhdFaceFluxStorage() = default;
  explicit MhdFaceFluxStorage(const MhdPatchGeometry& geometry);

  void resize(const MhdPatchGeometry& geometry);
  [[nodiscard]] const MhdPatchGeometry& geometry() const noexcept { return m_geometry; }
  [[nodiscard]] std::span<double> fluxX() noexcept { return m_phi_x.span(); }
  [[nodiscard]] std::span<const double> fluxX() const noexcept { return m_phi_x.span(); }
  [[nodiscard]] std::span<double> fluxY() noexcept { return m_phi_y.span(); }
  [[nodiscard]] std::span<const double> fluxY() const noexcept { return m_phi_y.span(); }
  [[nodiscard]] std::span<double> fluxZ() noexcept { return m_phi_z.span(); }
  [[nodiscard]] std::span<const double> fluxZ() const noexcept { return m_phi_z.span(); }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept {
    return m_phi_x.ownedCapacityBytes() + m_phi_y.ownedCapacityBytes() +
        m_phi_z.ownedCapacityBytes();
  }
  [[nodiscard]] double faceFluxX(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return m_phi_x[m_geometry.faceIndexX(i, j, k)];
  }
  [[nodiscard]] double faceFluxY(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return m_phi_y[m_geometry.faceIndexY(i, j, k)];
  }
  [[nodiscard]] double faceFluxZ(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return m_phi_z[m_geometry.faceIndexZ(i, j, k)];
  }
  void setFaceFluxX(std::size_t i, std::size_t j, std::size_t k, double value) noexcept {
    m_phi_x[m_geometry.faceIndexX(i, j, k)] = value;
  }
  void setFaceFluxY(std::size_t i, std::size_t j, std::size_t k, double value) noexcept {
    m_phi_y[m_geometry.faceIndexY(i, j, k)] = value;
  }
  void setFaceFluxZ(std::size_t i, std::size_t j, std::size_t k, double value) noexcept {
    m_phi_z[m_geometry.faceIndexZ(i, j, k)] = value;
  }
  // Derived cell-centered comoving field, arithmetic mean of the two bordering
  // face fields (specification section 5.1).
  [[nodiscard]] std::array<double, 3> cellFieldComoving(
      std::size_t i, std::size_t j, std::size_t k) const noexcept;

 private:
  MhdPatchGeometry m_geometry{};
  MhdDoubleArray m_phi_x;
  MhdDoubleArray m_phi_y;
  MhdDoubleArray m_phi_z;
};

// Non-owning face flux view used by solvers. The geometry epoch is part of the
// identity: a stale view must be rebound after regrid/migration.
struct MhdFaceFluxView {
  const MhdPatchGeometry* geometry = nullptr;
  std::span<const double> phi_x;
  std::span<const double> phi_y;
  std::span<const double> phi_z;
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  int refinement_level = 0;

  [[nodiscard]] bool valid() const noexcept;
};

struct MhdFaceFluxMutableView {
  const MhdPatchGeometry* geometry = nullptr;
  std::span<double> phi_x;
  std::span<double> phi_y;
  std::span<double> phi_z;
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  int refinement_level = 0;

  [[nodiscard]] bool valid() const noexcept;
};

[[nodiscard]] MhdFaceFluxView makeFaceFluxView(const MhdFaceFluxStorage& storage,
                                               std::uint64_t patch_id,
                                               std::uint64_t geometry_generation,
                                               int refinement_level) noexcept;
[[nodiscard]] MhdFaceFluxMutableView makeMutableFaceFluxView(
    MhdFaceFluxStorage& storage,
    std::uint64_t patch_id,
    std::uint64_t geometry_generation,
    int refinement_level) noexcept;

// Runtime stage context shared by every MHD operator. Scale factor and
// Hubble rate are SUPPLIED by the authoritative cosmological timeline (MHD
// never recomputes a(t)); isolated runs use a=1, H=0 at both endpoints.
//   scale_factor / hubble_rate_code         : values at time_code (t_n)
//   scale_factor_end / hubble_rate_end_code : values at time_code + dt_code
// SSPRK2 stage 1 is evaluated at t_n and stage 2 at t_n + dt, so both
// endpoints are required. `apply_expansion_source` selects whether the MHD
// operator owns the comoving expansion source (m_dot=-H m, E_dot=-H(2K+3p+M))
// inside its stages; the workflow sets it false when another owner already
// applies the gas Hubble source so the drag is never applied twice.
struct MhdStageContext {
  double dt_code = 0.0;
  double time_code = 0.0;
  double scale_factor = 1.0;
  double hubble_rate_code = 0.0;
  double scale_factor_end = 1.0;
  double hubble_rate_end_code = 0.0;
  int refinement_level = 0;
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::uint32_t stage_ordinal = 0U;
  bool comoving_coordinates = false;
  bool apply_expansion_source = true;
};

// Linear interpolation of the timeline-supplied endpoints at a fraction of
// the interval. Used only for operator-split sub-intervals (D/H halves).
[[nodiscard]] MhdStageContext subIntervalContext(const MhdStageContext& interval,
                                                 double begin_fraction,
                                                 double end_fraction) noexcept;

// Conservation audit for one accepted MHD interval. All energies are comoving
// code-unit volume integrals. Magnetic energy uses B_c^2/(2a).
struct MhdConservationBudget {
  double mass = 0.0;
  std::array<double, 3> momentum{};
  double energy_thermal = 0.0;
  double energy_kinetic = 0.0;
  double energy_magnetic = 0.0;
  double energy_total = 0.0;
  double boundary_energy_flux = 0.0;
  double source_energy_work = 0.0;
  double floor_added_energy = 0.0;
  double divergence_max_normalized = 0.0;
  double divergence_normalizer = 0.0;
};

// Diagnostic sink. Per-face/per-edge dispatchers are banned; the sink is
// called at patch/stage granularity only.
class MhdDiagnosticSink {
 public:
  virtual ~MhdDiagnosticSink() = default;
  virtual void recordStage(const MhdStageContext& context, const MhdConservationBudget& budget) = 0;
  virtual void recordCounters(std::string_view counter, double value) = 0;
};

// Retained capacity in bytes of one MhdDoubleArray lane of `count` doubles
// (64-byte rounded, the exact allocator request). Exact allocation accounting
// for MHD stage state lives in cosmosim/mhd/mhd_memory.hpp.
[[nodiscard]] std::uint64_t mhdAlignedLaneBytes(std::size_t count);

// One-sided 512 B/active-real-cell engineering cap for MHD-owned phase-resident
// plus scratch storage. Excludes persistent face flux and MPI staging.
inline constexpr std::uint64_t k_mhd_max_scratch_bytes_per_active_cell = 512U;
inline constexpr std::uint64_t k_mhd_default_comm_budget_bytes_per_rank = 64ULL * 1024ULL * 1024ULL;
inline constexpr std::uint64_t k_mhd_max_comm_budget_bytes_per_rank = 128ULL * 1024ULL * 1024ULL;
inline constexpr std::uint64_t k_mhd_io_chunk_target_bytes_per_rank = 8ULL * 1024ULL * 1024ULL;

}  // namespace cosmosim::mhd
