#pragma once

// CHUI-RT-BH-AGN-REF-v1 reference BH engine (physics.bh_engine_model =
// reference_v1): per-source-interval coordinator of supply estimation (BH03),
// conservative gas withdrawal (BH04), disk reservoir (BH05), horizon ledger and
// emission budgets (BH06), unresolved dynamical friction (BH02) and bound-pair
// mergers (BH07). Seeding (BH01) is BlackHoleHaloSeeder (bh_seed_provider.hpp).
//
// Chronology: called ONCE per accepted source interval at the discrete-event
// epoch (coarse interval end, every AMR level synchronized, no open flux
// register), replacing — never adding to — the legacy BlackHoleAgnModel
// callback. All gas mutation goes through SourceTransactionBatch with a
// unanimous collective vote before any commit. Pending feedback budgets are
// written to the chui_bh_engine_v1 store and consumed by AgnFeedbackInjector
// in the same epoch (and by RT as time-integrated sources).
//
// Units: G, c, m_p, sigma_T in code units from ProperFrameContext; rates use
// the PROPER interval dt_proper_code = physical_seconds / time_si_per_code;
// densities/lengths/volumes proper; velocities peculiar.

#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/config_rt_bh_agn.hpp"
#include "cosmosim/physics/black_hole/bh_budget.hpp"
#include "cosmosim/physics/black_hole/bh_engine_state.hpp"
#include "cosmosim/physics/black_hole/bh_neighborhood.hpp"
#include "cosmosim/physics/rt_bh_agn_frame.hpp"

namespace cosmosim::core {
class SimulationState;
class HierarchicalTimeBinScheduler;
}

namespace cosmosim::physics {

struct BlackHoleEngineParameters {
  core::BlackHoleEngineConfig model{};
  BhChannelEfficiencies efficiencies{};
  double alpha_bondi = 1.0;
  double proton_mass_si = k_proton_mass_si;
  double thomson_cross_section_si = k_thomson_cross_section_si;
  double viscous_time_code = 0.0;          // proper code time
  double kernel_radius_min_proper_code = 0.0;
  double kernel_radius_max_proper_code = 0.0;
  double mass_floor_code = 1.0e-20;
  double specific_internal_energy_floor_code = 0.0;
  double max_withdrawal_fraction = 0.5;    // per row per interval (numerical safety limit; unmet supply recorded)
  double bh_softening_stored_code = 0.0;   // BH gravitational softening (stored/comoving)
  PeriodicBox box{};
  std::uint32_t max_disk_substeps = 4096U;
};

struct BhIntervalContext {
  ProperFrameContext frame{};
  double time_begin_code = 0.0;
  double time_end_code = 0.0;
  double dt_proper_code = 0.0;
  std::uint64_t interval_index = 0U;
  BhOwnedLeafView leaves{};
  core::HierarchicalTimeBinScheduler* particle_scheduler = nullptr;  // merger removal remap (may be null in tests)
};

struct BhEngineIntervalReport {
  std::uint64_t local_black_holes = 0U;
  double requested_mass_code = 0.0;
  double withdrawn_mass_code = 0.0;
  double unmet_mass_code = 0.0;
  double inner_mass_code = 0.0;
  double disk_wind_mass_code = 0.0;
  double mass_defect_code = 0.0;
  double radiated_energy_code = 0.0;
  double mechanical_energy_code = 0.0;
  double captured_gas_energy_code = 0.0;
  double df_reservoir_energy_code = 0.0;
  std::uint64_t kernel_inadequate = 0U;
  std::uint64_t torque_inadequate = 0U;
  std::uint64_t eddington_capped = 0U;
  std::uint64_t df_applied = 0U;
  std::uint64_t mergers = 0U;
  std::uint64_t disk_substeps = 0U;
  double max_transaction_residual = 0.0;
  std::vector<BhEmissionBudget> budgets;
  std::vector<std::uint32_t> touched_gas_rows;
  std::vector<std::uint64_t> removed_particle_ids;
  bool gas_mutated = false;
  bool particle_mass_mutated = false;
};

class BlackHoleEngine {
 public:
  explicit BlackHoleEngine(BlackHoleEngineParameters parameters);

  [[nodiscard]] const BlackHoleEngineParameters& parameters() const noexcept { return m_parameters; }

  // Collective: every rank calls it with its local state (empty ranks too).
  [[nodiscard]] BhEngineIntervalReport advanceInterval(core::SimulationState& state, const BhIntervalContext& context,
                                                       BhSupportExchange& exchange) const;

 private:
  BlackHoleEngineParameters m_parameters;
};

[[nodiscard]] BlackHoleEngineParameters makeBlackHoleEngineParameters(const core::PhysicsConfig& physics,
                                                                      const core::UnitSystem& units,
                                                                      double bh_softening_stored_code,
                                                                      const PeriodicBox& box,
                                                                      double mass_floor_code);

}  // namespace cosmosim::physics
