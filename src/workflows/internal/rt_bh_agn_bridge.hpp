#pragma once

// CHUI-RT-BH-AGN-REF-v1 production composition of the reference BH engine,
// halo seeding, AGN channel injection and radiative transfer inside the
// SourceRuntime discrete-event epoch (INT01).
//
// Chronology at the coarse-interval end epoch (all AMR levels synchronized,
// no open flux register, every rank participating):
//   star birth -> stellar returns/feedback (existing) ->
//   [this bridge] BH supply/withdrawal/disk/ledger/DF/mergers ->
//   halo seeding (cadence) -> AGN injection of pending budgets ->
//   RT interval (AGN radiative budget emitted exactly once, PD-ARS IMEX,
//   staged radiation+matter commit).
// The legacy BlackHoleAgnModel callback is NOT executed for an interval
// handled here (one owner per source; selected by physics.bh_engine_model).
// Caller responsibilities after the epoch: compatibility-mirror sync,
// thermodynamic cache refresh of touched rows, scheduler append for newborn
// BHs, newly-created-ID ledger, one collective gravity-source generation bump.

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/agn/agn_injection.hpp"
#include "cosmosim/physics/black_hole/bh_engine.hpp"
#include "cosmosim/physics/black_hole/bh_neighborhood.hpp"
#include "cosmosim/radiation/radiation_runtime.hpp"
#include "cosmosim/radiation/radiation_spectrum.hpp"
#include "cosmosim/radiation/radiation_stellar_sources.hpp"

namespace cosmosim::core {
class HierarchicalTimeBinScheduler;
}
namespace cosmosim::physics {
class ParticleIdPrecommit;
}

namespace cosmosim::workflows::internal {

// Distributed owner-target support exchange (bounded alltoall; rank boxes of
// owned leaf cells allgathered once per epoch: O(ranks) metadata, no global
// cell/BH catalog). Replies are matched to the queries of the immediately
// preceding routeQueries() call by position.
class MpiBhSupportExchange final : public physics::BhSupportExchange {
 public:
  MpiBhSupportExchange(const parallel::MpiContext& mpi, const core::SimulationState& state,
                       std::span<const std::uint8_t> owned_leaf_mask, const physics::PeriodicBox& box);
  [[nodiscard]] std::vector<physics::BhSupportQuery> routeQueries(
      std::span<const physics::BhSupportQuery> local_queries) override;
  [[nodiscard]] std::vector<physics::BhSupportReply> returnReplies(
      std::span<const physics::BhSupportReply> replies) override;
  [[nodiscard]] bool allAgree(bool local_ok) override;
  [[nodiscard]] std::uint32_t worldRank() const noexcept override;

 private:
  [[nodiscard]] bool sphereIntersectsRank(const physics::BhSupportQuery& q, int rank) const;
  const parallel::MpiContext& m_mpi;
  physics::PeriodicBox m_box;
  std::vector<std::array<double, 7>> m_rank_boxes;  // lo[3], hi[3], nonempty flag
  std::vector<physics::BhSupportQuery> m_last_routed;
};

struct RtBhAgnEpochInputs {
  double time_begin_code = 0.0;
  double time_end_code = 0.0;
  double scale_factor_begin = 1.0;
  double scale_factor_end = 1.0;
  double physical_seconds = 0.0;
  double hubble_begin_si = 0.0;
  double hubble_end_si = 0.0;
  bool cosmological = false;
  std::uint64_t step_index = 0U;
  std::uint64_t interval_index = 0U;
  std::span<const std::uint8_t> owned_leaf_mask;
  core::HierarchicalTimeBinScheduler* particle_scheduler = nullptr;
  physics::ParticleIdPrecommit* id_precommit = nullptr;
};

struct RtBhAgnEpochReport {
  bool gas_mutated = false;
  bool particle_mass_mutated = false;
  std::vector<std::uint32_t> touched_gas_rows;
  std::vector<std::uint64_t> seeded_particle_ids;
  std::vector<std::uint64_t> removed_particle_ids;
  physics::BhEngineIntervalReport engine{};
  physics::AgnInjectionReport agn{};
  radiation::RadiationIntervalReport radiation{};
  radiation::StellarRadiationSourceReport stellar_sources{};  // F02 emission ledger
  std::uint64_t seeded = 0U;
};

class RtBhAgnSourceBridge {
 public:
  RtBhAgnSourceBridge(const core::SimulationConfig& config, const core::UnitSystem& units,
                      const parallel::MpiContext& mpi, std::array<double, 3> box_extent_comoving, bool periodic);
  [[nodiscard]] bool engineEnabled() const noexcept { return m_engine.has_value(); }
  [[nodiscard]] bool radiationEnabled() const noexcept { return m_radiation.has_value(); }
  // F01: forwards the run's MemoryGovernor to the RT interval workspaces.
  void setMemoryGovernor(core::MemoryGovernor* governor) noexcept {
    if (m_radiation.has_value()) m_radiation->setMemoryGovernor(governor);
  }
  [[nodiscard]] RtBhAgnEpochReport executeDiscreteEpoch(core::SimulationState& state, const RtBhAgnEpochInputs& in);

 private:
  [[nodiscard]] std::vector<physics::BhHaloSeedCandidate> serialHaloCandidates(const core::SimulationState& state,
                                                                               const RtBhAgnEpochInputs& in) const;
  core::SimulationConfig m_config;
  core::UnitSystem m_units;
  const parallel::MpiContext& m_mpi;
  physics::PeriodicBox m_box;
  std::optional<physics::BlackHoleEngine> m_engine;
  std::optional<physics::AgnFeedbackInjector> m_agn;
  std::optional<radiation::RadiationRuntime> m_radiation;
  std::unique_ptr<radiation::TabulatedSpectrum> m_agn_sed;
  std::vector<double> m_agn_group_fractions;
  // CHUI-EXT-PHYS-REF-v1 F02 stellar luminosity provider (null when disabled).
  std::unique_ptr<radiation::StellarLuminosityTable> m_stellar_table;
  std::unique_ptr<core::LambdaCdmBackground> m_background;
};

}  // namespace cosmosim::workflows::internal
