#pragma once

// CHUI-RT-BH-AGN-REF-v1, BH companion section 2: canonical advanced BH state.
//
// Authority (amendment A2): the reference engine's persistent per-BH truth is
// a particle-indexed module sidecar `chui_bh_engine_v1` keyed by the stable BH
// particle ID. It migrates with the particle, survives reorder/compaction and
// is persisted in native restart through the existing ModuleSidecarBlock
// machinery (independent schema version). It is allocated only for BH rows of
// runs that select physics.bh_engine_model=reference_v1; there is no cost per
// DM/gas particle. BlackHoleParticleSidecar lanes keep their legacy meaning;
// for the reference engine `subgrid_mass_code` mirrors the horizon mass and
// `accretion_rate_code` the accepted inner-engine rate (derived mirrors,
// refreshed at every commit). Particle gravitating mass is
//   m_dyn = M_bh + M_disk + M_pending_ejecta     (no double counting).
// Derived quantities (a*, L/L_Edd, r_B/dx, kernel radius) are diagnostics.

#include <array>
#include <cstdint>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

inline constexpr const char* k_bh_engine_module = "chui_bh_engine_v1";
inline constexpr std::uint32_t k_bh_engine_schema = 1U;

enum BhEngineFlag : std::uint32_t {
  kBhFlagInitialized = 1U << 0U,
  kBhFlagKernelInadequate = 1U << 1U,     // fewer than kernel_min_cells / no support
  kBhFlagBondiUnresolved = 1U << 2U,      // r_B < dx (always true for subgrid Bondi)
  kBhFlagTorqueInadequate = 1U << 3U,     // decomposition samples/axis unreliable
  kBhFlagEddingtonCapped = 1U << 4U,
  kBhFlagSupportTruncated = 1U << 5U,     // rate limited by availability/floors
  kBhFlagPendingFeedbackHeld = 1U << 6U,  // feedback held (no support / dilution / velocity domain)
  kBhFlagSeededByHalo = 1U << 7U,
  kBhFlagMergedRemnant = 1U << 8U,
};

// Trivially copyable, fixed-layout row (bytes are the persistent schema).
struct BlackHoleEngineRecord {
  double horizon_mass_code = 0.0;   // physical horizon/gravitational mass M_bh
  double disk_mass_code = 0.0;      // unresolved baryonic disk reservoir M_d
  std::array<double, k_chemical_element_count> disk_element_mass_code{};
  double disk_metal_mass_code = 0.0;  // legacy total-metal mirror of the disk
  std::array<double, 3> disk_angular_momentum_code{};  // J_d (captured orbital J; P1 evolves it)
  std::array<double, 3> spin_angular_momentum_code{};  // J_bh (P1; zero in P0)
  // Pending feedback budgets (restart truth until released or explicitly lost).
  double pending_thermal_energy_code = 0.0;
  double pending_kick_energy_code = 0.0;
  double pending_wind_energy_code = 0.0;
  double pending_wind_mass_code = 0.0;
  std::array<double, k_chemical_element_count> pending_wind_element_mass_code{};
  double pending_wind_metal_mass_code = 0.0;
  double pending_jet_energy_code = 0.0;
  double pending_radiation_energy_code = 0.0;
  double pending_magnetic_energy_code = 0.0;  // P1 (always zero in pass 01 production)
  // Cumulative ledgers (actual transfers, not requests).
  double cumulative_supply_requested_mass_code = 0.0;
  double cumulative_withdrawn_mass_code = 0.0;
  double cumulative_unmet_supply_mass_code = 0.0;
  double cumulative_inner_inflow_mass_code = 0.0;   // dm0 after disk winds
  double cumulative_disk_wind_mass_code = 0.0;
  double cumulative_mass_defect_code = 0.0;         // relativistic ledger: eps_acc dm0 (+ spin)
  double cumulative_radiated_energy_code = 0.0;
  double cumulative_mechanical_energy_code = 0.0;   // generated (thermal+kick+wind+jet)
  double cumulative_deposited_energy_code = 0.0;    // actually injected into gas
  double cumulative_escaped_radiation_energy_code = 0.0;
  double cumulative_captured_gas_energy_code = 0.0; // thermal+kinetic of withdrawn gas (unresolved reservoir)
  double cumulative_capture_dissipation_code = 0.0; // inelastic capture KE loss (unresolved reservoir)
  double cumulative_df_reservoir_energy_code = 0.0; // dynamical friction work not represented in resolved KE
  double cumulative_spin_extracted_energy_code = 0.0;
  // Last-interval diagnostics (cache; refreshed each interval).
  double last_supply_rate_code = 0.0;
  double last_inner_rate_code = 0.0;
  double last_eddington_rate_ratio = 0.0;       // mdot / mdot_Edd,ref
  double last_luminosity_eddington_ratio = 0.0; // L / L_Edd
  double last_kernel_radius_proper_code = 0.0;
  double last_bondi_radius_over_dx = 0.0;
  double last_hot_gas_mass_code = 0.0;
  double last_cold_gas_mass_code = 0.0;
  double source_time_code = 0.0;                // end time of the last accepted interval
  std::uint64_t interval_index = 0U;
  std::uint64_t event_sequence = 0U;
  std::uint64_t seed_halo_key = 0U;
  std::uint64_t merger_count = 0U;
  std::uint64_t last_merged_partner_id = 0U;
  std::uint64_t rng_counter = 0U;
  std::uint32_t flags = 0U;
  std::uint32_t reserved = 0U;

  [[nodiscard]] double pendingEjectaMassCode() const noexcept { return pending_wind_mass_code; }
  [[nodiscard]] double dynamicalMassCode() const noexcept {
    return horizon_mass_code + disk_mass_code + pending_wind_mass_code;
  }
};

class BlackHoleEngineStore {
 public:
  BlackHoleEngineStore() = default;
  // Load the module sidecar (empty store when absent). Throws on schema or
  // stride mismatch (never reinterprets foreign bytes).
  explicit BlackHoleEngineStore(const core::SimulationState& state);

  [[nodiscard]] bool has(std::uint64_t particle_id) const;
  [[nodiscard]] BlackHoleEngineRecord get(std::uint64_t particle_id) const;
  void put(std::uint64_t particle_id, const BlackHoleEngineRecord& record);
  void erase(std::uint64_t particle_id);
  // Initialize records of BH rows that have none, adopting the legacy lanes
  // (horizon = subgrid mass; disk/pending zero). Returns rows initialized.
  std::size_t adoptLegacyBlackHoles(const core::SimulationState& state);
  // Replace the module sidecar. Rows of particles that no longer exist
  // locally are dropped; a live BH row without a record is refused.
  void store(core::SimulationState& state) const;
  [[nodiscard]] std::size_t size() const noexcept { return m_ids.size(); }
  [[nodiscard]] const std::vector<std::uint64_t>& ids() const noexcept { return m_ids; }

 private:
  std::vector<std::uint64_t> m_ids;  // sorted stable particle IDs
  std::vector<BlackHoleEngineRecord> m_records;
};

}  // namespace cosmosim::physics
