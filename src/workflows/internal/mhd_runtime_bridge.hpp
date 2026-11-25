#pragma once

// Workflow-side composition of the CHUI-MHD-REF-v1 synchronized reference path
// (pass 02, UNVALIDATED source; never compiled or run in the source-only
// campaign). See docs/architecture/chui_mhd_implementation_pass02_status.md.
//
// Ownership:
//  * core::SimulationState::gas_cells (density, peculiar velocity, specific
//    internal energy, metal mass) stays the gas authority; the bridge stages
//    a comoving conserved copy (E_c = rho e + K + B_c^2/(2a)) only for the
//    duration of one kHydroUpdate and writes it back after acceptance.
//  * core::SimulationState::magnetic_faces is the single magnetic authority.
//  * The leaf AMR tiling of core::SimulationState::patches (explicit patch
//    geometry) is the topology; the bridge never creates or moves patches.
//  * The global timestep scheduler and the global-step retry transaction of
//    the TimeCoordinator own dt and retries: a rejected MHD interval is
//    returned as StageTransactionRejected(kNumericalRetry) with nothing
//    committed (MhdRetryMode::kExternal).

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/mhd/mhd_level_stepper.hpp"
#include "cosmosim/mhd/mhd_runtime_adapter.hpp"
#include "cosmosim/mhd/mhd_transport.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows::internal {

struct MhdDomainSpec {
  std::array<double, 3> origin_comoving{};
  std::array<double, 3> extent_comoving{};
  core::BoundaryCondition boundary = core::BoundaryCondition::kPeriodic;
};

// MhdHaloTransport on top of the existing bounded MPI byte all-to-all and the
// MpiContext reductions (no second MPI layer).
class MpiMhdHaloTransport final : public mhd::MhdHaloTransport {
 public:
  explicit MpiMhdHaloTransport(const parallel::MpiContext& mpi_context) noexcept
      : m_mpi(mpi_context) {}
  [[nodiscard]] int rank() const noexcept override;
  [[nodiscard]] int size() const noexcept override;
  [[nodiscard]] std::vector<std::vector<std::uint8_t>> exchange(
      const std::vector<std::vector<std::uint8_t>>& send) override;
  [[nodiscard]] std::uint64_t allreduceMaxUint64(std::uint64_t value) override;
  [[nodiscard]] double allreduceMinDouble(double value) override;
  [[nodiscard]] double allreduceMaxDouble(double value) override;
  void allreduceSumDoubles(std::span<double> values) override;

 private:
  const parallel::MpiContext& m_mpi;
};

struct MhdRuntimeTelemetry {
  std::uint64_t accepted_intervals = 0U;
  std::uint64_t rejected_intervals = 0U;
  std::uint64_t ideal_stage_evaluations_attempted = 0U;
  std::uint64_t ideal_stage_evaluations_accepted = 0U;
  std::uint64_t hlld_fallback_faces = 0U;
  std::uint64_t first_order_faces = 0U;
  std::uint64_t hall_subcycles = 0U;
  std::uint64_t rkl_stages = 0U;
  std::uint64_t rkl_cycles = 0U;
  std::uint64_t mpi_bytes_sent = 0U;
  std::uint64_t mpi_bytes_received = 0U;
  std::uint64_t ghost_prolongation_slope_fallbacks = 0U;
  std::uint64_t scratch_capacity_bytes = 0U;
  std::uint64_t staging_capacity_bytes = 0U;
  double divergence_max_normalized = 0.0;
  double last_dt_accepted_code = 0.0;
  double last_magnetic_energy = 0.0;
  double last_total_energy = 0.0;
};

class MhdRuntimeBridge {
 public:
  MhdRuntimeBridge(const core::SimulationConfig& config, const parallel::MpiContext& mpi_context,
                   core::MemoryGovernor* memory_governor, const core::UnitSystem& units);
  ~MhdRuntimeBridge();
  MhdRuntimeBridge(const MhdRuntimeBridge&) = delete;
  MhdRuntimeBridge& operator=(const MhdRuntimeBridge&) = delete;

  // COLLECTIVE. Ensure the authoritative magnetic state matches the current
  // leaf topology. Fresh start (empty magnetic state at step 0): build the
  // face layout from the real patch topology and seed it with the configured
  // Phi = Phi_uniform + C(A). Restart / later steps: validate only (never
  // reseed); an empty or stale magnetic state fails closed.
  void ensureMagneticState(core::SimulationState& state, const core::IntegratorState& integrator,
                           const MhdDomainSpec& domain);

  // COLLECTIVE scheduler candidate (rank-global): fast-magnetosonic CFL and
  // non-ideal cost-bound candidates at the current state; +infinity on an
  // empty hierarchy. Never mutates canonical state.
  [[nodiscard]] double stableTimeStep(core::SimulationState& state,
                                      const core::IntegratorState& integrator,
                                      const MhdDomainSpec& domain, bool cosmological);

  // COLLECTIVE kHydroUpdate for magnetized runs: one synchronized
  // D/H/I/H/D interval over [t, t + dt] at the timeline endpoints. On
  // rejection nothing is written and StageTransactionRejected is thrown.
  void advance(core::StepContext& context, const MhdDomainSpec& domain,
               double temperature_k_per_internal_energy_code);

  [[nodiscard]] const MhdRuntimeTelemetry& telemetry() const noexcept { return m_telemetry; }
  [[nodiscard]] const mhd::MhdLevelAdvanceResult& lastResult() const noexcept { return m_last; }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;

 private:
  struct LocalPatch {
    std::uint64_t patch_id = 0U;
    std::uint64_t geometry_generation = 0U;
    std::int32_t level = 0;
    std::array<std::int64_t, 3> origin_cell{};
    mhd::MhdPatchGeometry geometry{};
    std::vector<std::uint32_t> rows;  // patch-local cell -> dense gas row
    std::size_t face_record = 0U;
    mhd::MhdGasStateSoa gas{};
  };

  void bind(core::SimulationState& state, const MhdDomainSpec& domain);
  void buildDirectory(const core::SimulationState& state, const MhdDomainSpec& domain);
  void seedFreshField(core::SimulationState& state, const MhdDomainSpec& domain,
                      double scale_factor);
  void bindFaces(core::SimulationState& state);
  void loadGas(const core::SimulationState& state, double scale_factor);
  [[nodiscard]] mhd::MhdLevelAdvanceOptions advanceOptions() const;
  [[nodiscard]] std::uint64_t topologyFingerprint(const core::SimulationState& state) const;

  const core::SimulationConfig& m_config;
  const parallel::MpiContext& m_mpi;
  core::MemoryGovernor* m_governor = nullptr;
  core::UnitSystem m_units;
  MpiMhdHaloTransport m_transport;
  mhd::MhdNumericsPolicy m_policy{};
  std::unique_ptr<mhd::MhdCoefficientProvider> m_provider;
  std::unique_ptr<mhd::MhdNonIdealRateOperator> m_rate;
  mhd::MhdLevelDomain m_level_domain{};
  std::vector<LocalPatch> m_patches;
  std::vector<mhd::MhdRemotePatchDescriptor> m_remote;
  std::unique_ptr<mhd::MhdLevelStepper> m_stepper;
  std::uint64_t m_fingerprint = 0U;
  const double* m_bound_face_data = nullptr;
  std::uint64_t m_bound_identity_generation = 0U;
  core::MemoryReservation m_staging_reservation{};
  mhd::MhdLevelAdvanceResult m_last{};
  MhdRuntimeTelemetry m_telemetry{};
};

}  // namespace cosmosim::workflows::internal
