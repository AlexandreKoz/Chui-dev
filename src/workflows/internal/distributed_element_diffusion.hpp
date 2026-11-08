#pragma once

// G2.1 — distributed element-resolved diffusion (CHUI Part 3 closure).
//
// Protocol (all steps collective over the world communicator; ranks without
// gas, patches or faces participate with empty buffers):
//
//  1. gatherRemoteDiffusionCells: one directed AMR boundary exchange (same
//     planner as hydro, one ghost layer, periodic images) delivers the
//     remote owned leaf cells that touch local patches, with step-start
//     velocity and the owner's patch geometry.
//  2. buildMetalDiffusionTopology(..., remote cells, global envelope) yields
//     the local face graph and the local<->remote interface faces.
//  3. DistributedElementDiffusionCoupling::setup: the UNIQUE FLUX AUTHORITY of
//     every remote face is the side with the smaller gas_cell_id. The other
//     side sends a claim (its rho*kappa, area, distance). The authority
//     matches claims one-to-one (a missing or extra face on either side is a
//     collective invalid-state error: an asymmetric interface would silently
//     break conservation), computes the conductance once and returns it, so
//     both sides add the identical value to their stable-dt bound.
//  4. Per subcycle (inside advanceElementResolvedDiffusion): non-authority
//     sides send FRESH lane values and composition gas mass of their
//     interface rows; the authority evaluates each face exactly once with the
//     local face formula and returns the equal-and-opposite lane transfer.
//     Imported transfers are applied in (gas_cell_id, source rank) order.
//
// Guarantees: every lane total is conserved to roundoff across ranks (each
// face transfer is computed once and applied with opposite sign on the two
// owners); positivity uses the global stable dt; subcycle count, refusal and
// commit are global decisions. Decomposition invariance holds to roundoff
// (summation order of a row's face contributions may differ between
// decompositions); bitwise identity across decompositions is NOT claimed.
// Units: faces and conductances are physical (scaled by a) exactly as in the
// single-rank operator; dt is the accepted code interval.
// Wire (P0-49-01): interface lane state and deltas are exchanged as chunked
// schema-2 parallel::CompositionLaneChunkRecord streams bound to the layout
// signature, so every constructible layout (ion-rich C/O/Si/Fe, dust, SGS;
// <= parallel::k_max_composition_layout_lanes) is supported.
// Memory: exchange records and lookup tables are admitted through the
// MemoryGovernor (owner "sources.element_diffusion.remote") before
// allocation; refusal is a collective failure before any mutation.

#include <array>
#include <cstdint>
#include <exception>
#include <optional>
#include <span>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/element_diffusion.hpp"
#include "workflows/internal/metal_diffusion_topology.hpp"

namespace cosmosim::workflows::internal {

struct RemoteDiffusionCellGather {
  std::vector<RemoteDiffusionCell> cells;
  // Global physical domain envelope (allreduced over all owned patches).
  std::optional<std::array<double, 6>> global_envelope_physical;
};

[[nodiscard]] RemoteDiffusionCellGather gatherRemoteDiffusionCells(
    const core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    core::MemoryGovernor* governor,
    double scale_factor,
    bool periodic,
    const std::array<double, 3>& periodic_box_extent_comoving,
    std::uint64_t sequence);

class DistributedElementDiffusionCoupling final : public physics::ElementDiffusionRemoteCoupling {
 public:
  DistributedElementDiffusionCoupling(const parallel::MpiContext& mpi_context,
                                      core::MemoryGovernor* governor,
                                      const core::SimulationState& state,
                                      std::span<const RemoteDiffusionFace> remote_faces,
                                      std::span<const double> rho_kappa_code,
                                      std::size_t lane_upper_bound,
                                      std::uint64_t phase_tag = 0U);

  // Collective topology validation and conductance agreement (step 3).
  void setup();

  void reduceFailure(std::exception_ptr local_failure, const char* stage) override;
  void addRemoteConductance(std::span<double> conductance_sum) override;
  [[nodiscard]] double globalMin(double value) override;
  void globalSum(std::span<double> values) override;
  void addRemoteFaceTransfers(std::span<const double> values, std::span<const double> gas_mass,
                              std::size_t lane_count, double dt_sub, std::span<double> delta) override;

  [[nodiscard]] std::uint64_t authorityFaceCount() const noexcept { return m_authority.size(); }
  [[nodiscard]] std::uint64_t dependentFaceCount() const noexcept { return m_dependent.size(); }
  [[nodiscard]] std::uint64_t admittedBytes() const noexcept { return m_admitted_bytes; }

  struct FaceKey {
    std::uint64_t authority_id = 0U;
    std::uint64_t other_id = 0U;
    std::uint8_t axis = 0U;
    std::int8_t authority_side = 0;  // +1: other cell on the +axis side of the authority cell
    std::uint8_t periodic = 0U;
    std::uint8_t reserved = 0U;
    friend bool operator<(const FaceKey& l, const FaceKey& r) noexcept {
      if (l.authority_id != r.authority_id) return l.authority_id < r.authority_id;
      if (l.other_id != r.other_id) return l.other_id < r.other_id;
      if (l.axis != r.axis) return l.axis < r.axis;
      if (l.authority_side != r.authority_side) return l.authority_side < r.authority_side;
      return l.periodic < r.periodic;
    }
    friend bool operator==(const FaceKey& l, const FaceKey& r) noexcept {
      return !(l < r) && !(r < l);
    }
  };

 private:
  struct AuthorityFace {
    FaceKey key{};
    std::uint32_t local_row = 0U;
    int other_rank = 0;
    double area_code = 0.0;
    double distance_code = 0.0;
    double conductance = 0.0;
  };
  struct DependentFace {
    FaceKey key{};
    std::uint32_t local_row = 0U;
    int authority_rank = 0;
    double conductance = 0.0;
  };

  const parallel::MpiContext& m_mpi;
  core::MemoryGovernor* m_governor = nullptr;
  const core::SimulationState& m_state;
  std::span<const RemoteDiffusionFace> m_faces;
  std::span<const double> m_rho_kappa;
  std::size_t m_lane_upper_bound = 0U;
  std::vector<AuthorityFace> m_authority;  // sorted by key
  std::vector<DependentFace> m_dependent;  // sorted by key
  // Distinct (local_row, authority_rank) exports, sorted.
  std::vector<std::pair<std::uint32_t, int>> m_exports;
  // gas_cell_id -> local row for rows receiving imported transfers (sorted).
  std::vector<std::pair<std::uint64_t, std::uint32_t>> m_dependent_rows;
  core::MemoryReservation m_reservation;
  std::uint64_t m_admitted_bytes = 0U;
  // Per-coupling collective subcycle counter bound into every chunked lane
  // record (stale/foreign records are refused on reassembly).
  std::uint64_t m_exchange_sequence = 0U;
  bool m_setup_done = false;
};

}  // namespace cosmosim::workflows::internal
