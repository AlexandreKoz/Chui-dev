#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/mhd/mhd_induction.hpp"
#include "cosmosim/mhd/mhd_state.hpp"

namespace cosmosim::mhd {

// One accepted coarse/fine synchronization interval. Register entries are only
// meaningful for an exactly matching interval identity.
struct MhdAmrIntervalId {
  std::int32_t coarse_level = 0;
  double time_begin_code = 0.0;
  double time_end_code = 0.0;
  std::uint64_t coarse_geometry_generation = 0U;
  std::int32_t owner_rank = 0;

  friend bool operator==(const MhdAmrIntervalId&, const MhdAmrIntervalId&) = default;
};

// Sparse time-integrated edge-EMMF contribution at a coarse/fine interface.
// `edge_hash` is a stable geometric key; `coarse_cell_index` locates the cell
// whose face receives the curl correction.
struct MhdEmfRegisterRecord {
  std::uint64_t edge_hash = 0U;
  MhdAxis edge_axis = MhdAxis::kX;
  std::int64_t i = 0;
  std::int64_t j = 0;
  std::int64_t k = 0;
  std::uint64_t coarse_patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::int32_t coarse_level = 0;
  std::int32_t source_rank = 0;
  std::int32_t substep_ordinal = 0;
  double integrated_emf = 0.0;
  bool fine_contribution = true;
};

// Bounded, key-sorted register. Memory grows with the interface edge count and
// the number of accepted intervals, never with every face or every substep:
// substep contributions are accumulated immediately into one entry.
class MhdEmfRegister {
 public:
  explicit MhdEmfRegister(std::uint64_t max_entries_per_interval = 0U);

  // Returns false when the deterministic capacity bound would be exceeded.
  bool addContribution(const MhdEmfRegisterRecord& record);
  void clear() noexcept;
  [[nodiscard]] std::size_t size() const noexcept { return m_records.size(); }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
  [[nodiscard]] std::span<const MhdEmfRegisterRecord> records() const noexcept {
    return {m_records.data(), m_records.size()};
  }
  // Sorted, duplicate-free deterministic order (edge hash, then rank, then
  // substep ordinal). A duplicate (key, rank, ordinal) is a protocol error.
  void finalizeDeterministicOrder();
  [[nodiscard]] bool coversInterval(const MhdAmrIntervalId& interval) const noexcept {
    return m_interval == interval && !m_records.empty();
  }
  [[nodiscard]] const MhdAmrIntervalId& interval() const noexcept { return m_interval; }
  void bindInterval(const MhdAmrIntervalId& interval);

 private:
  std::vector<MhdEmfRegisterRecord> m_records;
  std::uint64_t m_max_entries = 0U;
  MhdAmrIntervalId m_interval{};
  bool m_has_interval = false;
};

// Ratio-two AMR magnetic transfer and reflux-curl (sections 7.2-7.3).
class MhdAmrSynchronizer {
 public:
  static constexpr std::int32_t k_refinement_ratio = 2;

  // Child geometry of a 2x2x2 split: child index = di + 2*(dj + 2*dk).
  [[nodiscard]] static MhdPatchGeometry childGeometry(const MhdPatchGeometry& parent,
                                                      std::int32_t child_index);

  struct TransferResult {
    bool ok = false;
    bool used_piecewise_constant = false;
    double max_cell_divergence_residual = 0.0;
    std::string detail;
  };

  // Exactly conservative oriented sum: parent = sum(child boundary fluxes).
  static void restrictFaceFluxes(const MhdPatchGeometry& parent_geometry,
                                 std::span<double> parent_phi_x,
                                 std::span<double> parent_phi_y,
                                 std::span<double> parent_phi_z,
                                 const std::array<MhdFaceFluxStorage, 8>& children);

  // Divergence-preserving ratio-two prolongation. Boundary faces come from a
  // limited bilinear split whose sum equals the parent exactly; internal faces
  // use the opposite-parent-face mean (Toth-Roe compatible 1D reduction) and a
  // bounded deterministic divergence-repair sweep. If the repair cannot reach
  // the contract tolerance the caller must reject the regrid transaction.
  static TransferResult prolongFaceFluxes(const MhdPatchGeometry& parent_geometry,
                                          std::span<const double> parent_phi_x,
                                          std::span<const double> parent_phi_y,
                                          std::span<const double> parent_phi_z,
                                          std::array<MhdFaceFluxStorage, 8>& children,
                                          double divergence_tolerance = 1.0e-12);

  // Apply Phi_coarse <- Phi_coarse - C(I_fine - I_coarse) using the register
  // entries that belong to the coarse patch. Incomplete intervals are refused.
  static void applyRefluxCurl(const MhdPatchGeometry& coarse_geometry,
                              std::span<double> coarse_phi_x,
                              std::span<double> coarse_phi_y,
                              std::span<double> coarse_phi_z,
                              const MhdEmfRegister& register,
                              std::uint64_t coarse_patch_id);

  // Sum of the magnetic energy (B_c^2/(2a)) implied by the face fluxes.
  [[nodiscard]] static double magneticEnergySum(const MhdPatchGeometry& geometry,
                                                std::span<const double> phi_x,
                                                std::span<const double> phi_y,
                                                std::span<const double> phi_z,
                                                double scale_factor);

  // Total boundary flux through the patch, per axis (net periodic flux check).
  [[nodiscard]] static std::array<double, 3> netBoundaryFlux(
      const MhdPatchGeometry& geometry,
      std::span<const double> phi_x,
      std::span<const double> phi_y,
      std::span<const double> phi_z);
};

// Upper bound on register entries for one coarse patch of the given geometry:
// every edge of the coarse cell shell adjacent to a covered face. Computed
// before stepping so admission never depends on runtime growth.
[[nodiscard]] std::uint64_t mhdEmfRegisterBound(const MhdPatchGeometry& coarse_geometry);

}  // namespace cosmosim::mhd
