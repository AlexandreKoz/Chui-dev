#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::core {
class SimulationState;
class MemoryGovernor;
}

namespace cosmosim::physics {

// ---------------------------------------------------------------------------
// Conservative composition advection (CHUI Part 2).
//
// Update law (consistent multi-fluid advection, Plewa & Mueller 1999, A&A
// 342, 179, first-order donor concentration): for every final hydro face with
// mass transfer dm (owner -> neighbor), each extensive lane s moves
//     dM_s = dm * X_s(donor),   X_s = M_s / sum_e,d M   (step-start state)
// with donor = owner if dm >= 0 else neighbor. Because sum_s X_s = 1 over the
// element+dust lanes, the composition mass changes by exactly the hydro mass
// transfer, so composition and gas mass stay consistent; every lane is
// conserved face by face (antisymmetric pairs). The SGS energy lane is
// advected the same way with X = E_sgs / M (specific SGS energy).
//
// Ownership of a face update:
//  * real/real (incl. periodic wrapped targets): both sides updated here.
//  * real/ghost: only the real side; the neighbouring patch updates its own
//    cell from its own advance (same-level faces carry identical fluxes).
//  * AMR coarse register faces (role kCoarse): the coarse side is NOT updated
//    from the coarse flux; the fine patch's kFine faces apply the mirrored
//    transfer to the coarse cell (coarse_gas_cell_id). This is the
//    composition reflux: the coarse cell receives exactly the fine-face
//    transfers, so lanes are conserved across coarse/fine interfaces.
//  * imported MPI ghosts / remote AMR sources: FAIL CLOSED (Part 3 owns
//    distributed composition ghosts); a precise diagnostic is thrown.
// Deferred application: records are applied once at stage end with the
// step-start donor fractions, so patch iteration order cannot change results.
// Positivity: an outflow never exceeds the donor mass under the hydro CFL;
// the apply step verifies non-negativity and refuses (throws) otherwise.
// ---------------------------------------------------------------------------
struct CompositionTransportCommitReport {
  double max_relative_mismatch = 0.0;   // before reconciliation
  double reconciled_mass_code = 0.0;    // signed mass added to match hydro
  std::uint64_t reconciled_rows = 0U;   // rows with mismatch > 1e-10
};

class CompositionFluxAccumulator final : public hydro::HydroFaceMassFluxSink {
 public:
  static constexpr std::uint32_t k_unresolved = std::numeric_limits<std::uint32_t>::max();

  explicit CompositionFluxAccumulator(const core::SimulationState& state);

  // Begin one patch advance. dense_row_by_real_cell maps the geometry's real
  // rows [0, cellCount) to SimulationState dense rows. Ghost rows are resolved
  // through ghost_origin_gas_cell_id (indexed by ghost slot; 0 = use owner).
  void beginAdvance(const hydro::HydroPatchGeometry& geometry,
                    std::span<const std::uint32_t> dense_row_by_real_cell,
                    std::span<const std::uint64_t> ghost_origin_gas_cell_id_by_slot) override;
  void endAdvance() noexcept override;

  void recordFaceMassTransfer(std::size_t face_index, std::size_t owner_cell,
                              std::size_t neighbor_cell, double mass_transfer_code) override;

  // Two-phase application, usable inside StageTransaction:
  //  prepare(): may throw; builds the lane deltas from the UNMODIFIED
  //    step-start sidecar and validates positivity. No state mutation.
  //  commit(): adds the deltas, then reconciles each row with the hydro cell
  //    mass (density-floor mass or round-off) by scaling with the row's own
  //    composition; the reconciled mass is ledgered in the report.
  void prepare(const core::SimulationState& state);
  CompositionTransportCommitReport commit(core::SimulationState& state) noexcept;

  [[nodiscard]] std::size_t recordCount() const noexcept { return m_records.size(); }
  void reserve(std::size_t records) { m_records.reserve(records); }
  // Records + the lanes x rows delta buffer (phase-resident transport
  // buffer; the caller admits it through the MemoryGovernor).
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept {
    return m_records.capacity() * sizeof(Record) + m_delta.capacity() * sizeof(double);
  }
  [[nodiscard]] static std::uint64_t deltaBufferBytes(std::size_t lanes, std::size_t rows) noexcept {
    return static_cast<std::uint64_t>(lanes) * rows * sizeof(double);
  }

 private:
  struct Record {
    std::uint32_t owner_row = k_unresolved;     // dense row of face owner side
    std::uint32_t neighbor_row = k_unresolved;  // dense row of neighbor side
    double mass = 0.0;                          // owner -> neighbor
    std::uint8_t update_owner = 0;
    std::uint8_t update_neighbor = 0;
  };
  [[nodiscard]] std::uint32_t resolveRow(std::size_t storage_row, bool& is_real) const;

  const core::SimulationState& m_state;
  const hydro::HydroPatchGeometry* m_geometry = nullptr;
  std::span<const std::uint32_t> m_dense_row_by_real_cell;
  std::span<const std::uint64_t> m_ghost_origin_by_slot;
  std::vector<Record> m_records;
  std::vector<double> m_delta;
  bool m_prepared = false;
};

}  // namespace cosmosim::physics
