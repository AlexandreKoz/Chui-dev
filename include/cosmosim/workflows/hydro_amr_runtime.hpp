#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"
#include "cosmosim/workflows/runtime_module_registry.hpp"

namespace cosmosim::core {
struct MemoryReport;
}  // namespace cosmosim::core

namespace cosmosim::physics {
class EffectiveMultiphaseEosTable;
}  // namespace cosmosim::physics

namespace cosmosim::workflows {

struct RuntimeServices;
class SourceRuntime;

// Owner interface for the currently supported rung-zero hydro/AMR path.
// Communication counters are exposed as values so callers cannot reach the
// owner's ghost, patch, conserved-state, or geometry workspaces.
class HydroAmrRuntime {
 public:
  virtual ~HydroAmrRuntime() = default;

  virtual void execute(HydroAmrStageView& view) = 0;
  // Source owner is constructed after hydro/AMR in the runtime composition.
  // Binding is non-owning; the composition owns both runtimes for the full
  // workflow lifetime. It enables T02 level-local thermochemistry without a
  // second source implementation inside the AMR module.
  virtual void bindSourceRuntime(SourceRuntime* source_runtime) noexcept = 0;
  [[nodiscard]] virtual RuntimeTaskMemoryEstimate estimateMemory(
      const core::SimulationState&) const {
    return {0U, false, "hydro/AMR owner does not publish a complete memory model"};
  }

  [[nodiscard]] virtual const hydro::HydroProfileEvent&
  lastHydroProfile() const noexcept = 0;
  [[nodiscard]] virtual const core::HydroCflDiagnostics&
  lastHydroCflDiagnostics() const noexcept = 0;
  [[nodiscard]] virtual std::uint64_t
  ghostExchangeBytesRecent() const noexcept = 0;
  [[nodiscard]] virtual std::size_t
  remoteImportedGhostCount() const noexcept = 0;
  [[nodiscard]] virtual std::size_t
  remoteInterfaceFaceCount() const noexcept = 0;
  [[nodiscard]] virtual std::size_t
  remoteStaleInvalidPayloadCount() const noexcept = 0;
  // Production AMR uses directed patch/cell exchange rather than the legacy
  // fixed-Cartesian imported-ghost counters above. Expose narrow telemetry so
  // workflow qualification can assert that a real inter-rank AMR interface
  // was exercised without reaching into transport workspaces.
  [[nodiscard]] virtual std::uint64_t
  directedAmrRemoteInterfaceCount() const noexcept = 0;
  [[nodiscard]] virtual std::uint64_t
  directedAmrPatchCellRecordsReceived() const noexcept = 0;
  [[nodiscard]] virtual core::MemoryReport memoryReport() const = 0;
};

[[nodiscard]] std::unique_ptr<HydroAmrRuntime> makeHydroAmrRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const RuntimeServices& services,
    std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table = nullptr);

}  // namespace cosmosim::workflows
