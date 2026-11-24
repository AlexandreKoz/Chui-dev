#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/mhd/mhd_induction.hpp"
#include "cosmosim/mhd/mhd_state.hpp"
#include "cosmosim/mhd/mhd_step_operator.hpp"

namespace cosmosim::mhd {

// One coarse/fine synchronization interval. Register entries are meaningful
// only for an exactly matching interval identity.
struct MhdAmrIntervalId {
  std::int32_t coarse_level = 0;
  double time_begin_code = 0.0;
  double time_end_code = 0.0;
  std::uint64_t coarse_geometry_generation = 0U;
  std::uint64_t hierarchy_epoch = 0U;
  std::int32_t owner_rank = 0;

  friend bool operator==(const MhdAmrIntervalId&, const MhdAmrIntervalId&) = default;
};

// Full geometric key of a coarse interface edge: never a hash identity.
struct MhdCoarseEdgeKey {
  std::uint64_t coarse_patch_id = 0U;
  MhdAxis axis = MhdAxis::kX;
  // Coarse-patch-local edge lattice coordinates (cell index along `axis`,
  // node indices along the two tangential axes).
  std::int64_t i = 0;
  std::int64_t j = 0;
  std::int64_t k = 0;

  friend auto operator<=>(const MhdCoarseEdgeKey&, const MhdCoarseEdgeKey&) = default;
};

// Sparse, pre-declared, time-integrated EMF reflux register for one coarse
// patch and one synchronization interval (CHUI-MHD-REF v1.0 section 7.3).
//
// Every coarse edge on the coarse/fine interface is DECLARED before stepping
// (bounded capacity, computed from topology). During the interval it must
// receive exactly one time-integrated coarse contribution I_coarse and, for
// each of its two collinear fine edge segments, exactly one time-integrated
// contribution per fine substep (1 for synchronized levels, 2 for ratio-2
// subcycling). Contributions are accumulated immediately into the record
// (no N_substeps x N_edges history). finalize() rejects any gap, duplicate,
// unknown key or interval mismatch with MHD_AMR_REGISTER_INCOMPLETE.
class MhdEmfRegister {
 public:
  struct Record {
    MhdCoarseEdgeKey key{};
    double integrated_coarse = 0.0;
    double integrated_fine = 0.0;
    std::uint32_t fine_coverage_mask = 0U;  // bit (segment * substeps + substep)
    bool coarse_present = false;
  };

  MhdEmfRegister() = default;

  // Binds the interval and fine substep count (1..8). Clears prior state.
  void bindInterval(const MhdAmrIntervalId& interval, std::uint64_t coarse_patch_id,
                    std::uint64_t coarse_geometry_generation, std::int32_t fine_substeps);
  // Declare the interface edge set (sorted, unique enforced). Capacity is
  // reserved once here; contributions never grow the register.
  void declareInterfaceEdges(std::span<const MhdCoarseEdgeKey> keys);
  void addCoarse(const MhdAmrIntervalId& interval, const MhdCoarseEdgeKey& key,
                 double integrated_emf);
  void addFine(const MhdAmrIntervalId& interval, const MhdCoarseEdgeKey& key,
               std::int32_t fine_segment, std::int32_t fine_substep, double integrated_emf);
  // Verify complete coverage; throws MHD_AMR_REGISTER_INCOMPLETE.
  void finalize() const;
  void clear() noexcept;

  [[nodiscard]] std::size_t size() const noexcept { return m_records.size(); }
  [[nodiscard]] std::span<const Record> records() const noexcept {
    return {m_records.data(), m_records.size()};
  }
  [[nodiscard]] const MhdAmrIntervalId& interval() const noexcept { return m_interval; }
  [[nodiscard]] bool bound() const noexcept { return m_bound; }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;

 private:
  [[nodiscard]] Record& locate(const MhdAmrIntervalId& interval, const MhdCoarseEdgeKey& key);

  std::vector<Record> m_records;
  MhdAmrIntervalId m_interval{};
  std::uint64_t m_coarse_patch_id = 0U;
  std::uint64_t m_coarse_geometry_generation = 0U;
  std::int32_t m_fine_substeps = 1;
  bool m_bound = false;
  bool m_declared = false;
};

// Ratio-two AMR magnetic transfer and reflux-curl (sections 7.2-7.3).
//
// Child layout: a parent patch of (nx, ny, nz) cells (all even) is covered by
// eight children of the SAME (nx, ny, nz) shape at half spacing; child index
// = di + 2*(dj + 2*dk) covers parent cells [d*n/2, (d+1)*n/2) per axis.
class MhdAmrSynchronizer {
 public:
  static constexpr std::int32_t k_refinement_ratio = 2;

  [[nodiscard]] static MhdPatchGeometry childGeometry(const MhdPatchGeometry& parent,
                                                      std::int32_t child_index);

  struct TransferResult {
    bool ok = false;
    bool used_zero_slope_fallback = false;
    double max_cell_divergence_residual = 0.0;
    double magnetic_energy_change = 0.0;
    std::string detail;
  };

  // Exactly conservative oriented restriction: each parent face equals the sum
  // of the four fine subfaces covering it, each counted exactly once. A parent
  // face on an internal sibling interface is read from one child and verified
  // against the other (mismatch -> MHD_INVALID_FIELD, nothing written).
  static void restrictFaceFluxes(const MhdPatchGeometry& parent_geometry,
                                 std::span<double> parent_phi_x,
                                 std::span<double> parent_phi_y,
                                 std::span<double> parent_phi_z,
                                 const std::array<MhdFaceFluxStorage, 8>& children);

  // Divergence-free ratio-two prolongation (Balsara 2001 reconstruction, which
  // for ratio two coincides with the Toth & Roe 2002 compatible formulas):
  // inside each parent cell the field is the unique solenoidal polynomial
  //   B_x = a0 + ax x + ay y + az z + axx x^2 + axy x y + axz x z (cyclic),
  // matching on every parent face its mean and limited (MC) transverse slopes.
  // Fine face fluxes are exact integrals of that polynomial, so every fine
  // cell has zero discrete divergence when the parent cell does, fine
  // subfaces sum exactly to their parent face, and the shared subfaces of
  // sibling cells agree by construction (face-local data). With
  // `zero_slopes` the transverse slopes are zero (lower-order fallback).
  // Inputs must be finite with parent divergence <= tolerance.
  static TransferResult prolongFaceFluxes(const MhdPatchGeometry& parent_geometry,
                                          std::span<const double> parent_phi_x,
                                          std::span<const double> parent_phi_y,
                                          std::span<const double> parent_phi_z,
                                          std::array<MhdFaceFluxStorage, 8>& children,
                                          double divergence_tolerance = 1.0e-12,
                                          bool zero_slopes = false);

  // Transactional gas + magnetic prolongation. Gas mass, momentum and metals
  // are copied piecewise-constant (exactly conservative). Total energy is
  // conserved exactly per parent cell: the fine magnetic energy changes by
  // Delta M = mean(M_fine) - M_parent, and the parent thermal energy is
  // re-allocated as U' = U - Delta M uniformly over its children. If U' <= 0
  // the zero-slope magnetic split is tried; if still inadmissible the
  // transfer fails and NOTHING is written to the outputs.
  static TransferResult prolongPatch(const MhdPatchGeometry& parent_geometry,
                                     const MhdGasStateSoa& parent_gas,
                                     const MhdFaceFluxView& parent_faces,
                                     double scale_factor,
                                     double gamma_adiabatic,
                                     std::array<MhdFaceFluxStorage, 8>& children_faces,
                                     std::array<MhdGasStateSoa, 8>& children_gas,
                                     double divergence_tolerance = 1.0e-12);

  // Transactional restriction: faces by unique oriented sum, gas by volume
  // averaging of every conserved lane (total energy included, so aggregate
  // energy is exact); the parent thermal energy U = E - K - M(restricted B)
  // must stay positive or the transfer fails without writing.
  static TransferResult restrictPatch(const MhdPatchGeometry& parent_geometry,
                                      const std::array<MhdFaceFluxStorage, 8>& children_faces,
                                      const std::array<MhdGasStateSoa, 8>& children_gas,
                                      double scale_factor,
                                      double gamma_adiabatic,
                                      MhdGasStateSoa& parent_gas,
                                      std::span<double> parent_phi_x,
                                      std::span<double> parent_phi_y,
                                      std::span<double> parent_phi_z);

  // Phi_coarse <- Phi_coarse - C(I_fine - I_coarse) on the declared interface
  // edges of a finalized register (curl correction, never an overwrite).
  static void applyRefluxCurl(const MhdPatchGeometry& coarse_geometry,
                              std::span<double> coarse_phi_x,
                              std::span<double> coarse_phi_y,
                              std::span<double> coarse_phi_z,
                              const MhdEmfRegister& reflux_register);

  [[nodiscard]] static double magneticEnergySum(const MhdPatchGeometry& geometry,
                                                std::span<const double> phi_x,
                                                std::span<const double> phi_y,
                                                std::span<const double> phi_z,
                                                double scale_factor);

  // Net oriented flux through the two boundary planes of each axis.
  [[nodiscard]] static std::array<double, 3> netBoundaryFlux(
      const MhdPatchGeometry& geometry,
      std::span<const double> phi_x,
      std::span<const double> phi_y,
      std::span<const double> phi_z);
};

// Upper bound on register records for one coarse patch: every edge of the
// coarse patch boundary shell. Computed before stepping so admission never
// depends on runtime growth.
[[nodiscard]] std::uint64_t mhdEmfRegisterBound(const MhdPatchGeometry& coarse_geometry);

}  // namespace cosmosim::mhd
