#include "cosmosim/workflows/hydro_amr_runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/amr/amr_production_regrid.hpp"
#include "cosmosim/amr/amr_patch_indexing.hpp"
#include "cosmosim/parallel/composition_exchange.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/adaptive_resource_control.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/retained_capacity_transaction.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/physics/cooling_heating.hpp"
#include "cosmosim/physics/stellar_population.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/baryonic_runtime_config.hpp"
#include "cosmosim/physics/composition_transport.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "cosmosim/workflows/source_runtime.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/cartesian_gas_cell_layout.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/amr_migration_payload.hpp"
#include "workflows/internal/particle_ghost_runtime.hpp"
#include "workflows/internal/mpi_stage_collective_decision.hpp"
#include "cosmosim/core/stage_transaction.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::workflows {
namespace {

constexpr double k_gamma_adiabatic = 5.0 / 3.0;


// Final-closure predictive AMR wire record.  Positions are deliberately the
// only distributed payload: every rank resolves them against its own leaf
// patches, so a star predicted to cross a decomposition boundary refines the
// destination rank as well as the source rank.  This keeps the prediction
// independent of MPI ownership without introducing a second patch authority.
struct PredictiveRefinementPoint {
  double x_comoving = 0.0;
  double y_comoving = 0.0;
  double z_comoving = 0.0;
};
static_assert(std::is_trivially_copyable_v<PredictiveRefinementPoint>);

[[nodiscard]] bool finitePredictivePoint(const PredictiveRefinementPoint& point) noexcept {
  return std::isfinite(point.x_comoving) && std::isfinite(point.y_comoving) &&
      std::isfinite(point.z_comoving);
}

[[nodiscard]] std::vector<PredictiveRefinementPoint> allgatherPredictiveRefinementPoints(
    const parallel::MpiContext& mpi_context,
    std::span<const PredictiveRefinementPoint> local_points) {
  std::exception_ptr local_failure;
  try {
    for (const PredictiveRefinementPoint& point : local_points) {
      if (!finitePredictivePoint(point)) {
        throw std::invalid_argument("predictive AMR refinement point is non-finite");
      }
    }
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "predictive AMR local point validation");

  std::span<const std::uint8_t> local_bytes;
  local_failure = nullptr;
  try {
    const std::size_t byte_count = core::checkedSizeMultiply(
        local_points.size(), sizeof(PredictiveRefinementPoint),
        "predictive AMR point byte count");
    local_bytes = std::span<const std::uint8_t>(
        reinterpret_cast<const std::uint8_t*>(local_points.data()), byte_count);
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "predictive AMR point byte framing");

  const std::vector<std::uint8_t> gathered_bytes =
      mpi_context.allgatherBytesBounded(local_bytes);
  std::vector<PredictiveRefinementPoint> gathered;
  local_failure = nullptr;
  try {
    if (gathered_bytes.size() % sizeof(PredictiveRefinementPoint) != 0U) {
      throw std::runtime_error("predictive AMR all-gather returned a partial point record");
    }
    gathered.resize(gathered_bytes.size() / sizeof(PredictiveRefinementPoint));
    if (!gathered_bytes.empty()) {
      std::memcpy(gathered.data(), gathered_bytes.data(), gathered_bytes.size());
    }
    for (const PredictiveRefinementPoint& point : gathered) {
      if (!finitePredictivePoint(point)) {
        throw std::runtime_error("predictive AMR all-gather returned a non-finite point");
      }
    }
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "predictive AMR point reassembly validation");
  return gathered;
}

// Memory-aware OPTIONAL refinement budget.  Hard accuracy requests never use
// this limiter.  The estimate is intentionally conservative and models the
// incremental retained state created by one worst-case local refinement plus
// composition lanes.  Transactional candidate-state admission remains the
// hard memory authority when an operation is actually applied.
[[nodiscard]] std::size_t memoryAwareOptionalRegridBudget(
    const core::SimulationState& state,
    core::MemoryGovernor* governor,
    std::size_t configured_max_operations) {
  if (configured_max_operations == 0U || governor == nullptr) {
    return configured_max_operations;
  }
  const core::MemoryGovernorSnapshot snapshot = governor->snapshot();
  if (snapshot.headroom_bytes == std::numeric_limits<std::uint64_t>::max()) {
    return configured_max_operations;
  }
  std::size_t max_patch_cells = 1U;
  for (std::size_t patch = 0; patch < state.patches.size(); ++patch) {
    max_patch_cells = std::max<std::size_t>(max_patch_cells, state.patches.cell_count[patch]);
  }
  // A refinement replaces N rows by 8N rows: +7N retained rows and +7 patch
  // descriptors.  The byte model deliberately over-counts the dense/gas
  // sidecars and lookup overhead because this is only an advisory budget.
  const std::size_t composition_lane_bytes = state.gas_composition.enabled()
      ? core::checkedSizeMultiply(
            state.gas_composition.laneCount(), sizeof(double),
            "optional AMR composition lane bytes")
      : 0U;
  const std::size_t per_cell_bytes = core::checkedSizeAdd(
      24U * sizeof(double) + 6U * sizeof(std::uint64_t) + 8U * sizeof(std::uint32_t),
      composition_lane_bytes,
      "optional AMR retained per-cell bytes");
  const std::size_t incremental_cells = core::checkedSizeMultiply(
      max_patch_cells, 7U, "optional AMR incremental cells");
  const std::size_t estimated_bytes = core::checkedSizeAdd(
      core::checkedSizeMultiply(
          incremental_cells, per_cell_bytes,
          "optional AMR incremental state bytes"),
      7U * (16U * sizeof(std::uint64_t) + 16U * sizeof(double)),
      "optional AMR incremental patch bytes");
  if (estimated_bytes == 0U) {
    return configured_max_operations;
  }
  const core::DeterministicBatchSizingResult selected = core::selectDeterministicBatchSize(
      snapshot,
      core::DeterministicBatchSizingPolicy{
          .requested_max_items = static_cast<std::uint64_t>(configured_max_operations),
          .bytes_per_item = core::checkedIntegralNarrow<std::uint64_t>(
              estimated_bytes, "optional AMR budget bytes per operation"),
          .fixed_reserve_bytes = 0U,
          .minimum_items = 1U,
          .alignment_items = 1U,
          .headroom_use_basis_points = 5000U,
      });
  return core::checkedIntegralNarrow<std::size_t>(
      selected.selected_items, "optional AMR memory-aware operation budget");
}

// Mapping of the core-owned configuration selectors onto hydro types. Policy
// resolution happens once at runtime construction / per step, never inside
// per-face loops.
[[nodiscard]] hydro::HydroNumericsPolicy hydroNumericsPolicyFromConfig(const core::SimulationConfig& config) {
  hydro::HydroNumericsPolicy policy;
  policy.integration_scheme =
      config.numerics.hydro_integration_scheme == core::HydroIntegrationMode::kCornerTransportUpwind
          ? hydro::HydroIntegrationScheme::kCornerTransportUpwind
          : hydro::HydroIntegrationScheme::kMusclHancock;
  policy.dual_energy.mode =
      config.numerics.hydro_dual_energy == core::HydroDualEnergyConfigMode::kThermalFraction
          ? hydro::HydroDualEnergyMode::kThermalFraction
          : hydro::HydroDualEnergyMode::kOff;
  policy.dual_energy.thermal_fraction_switch = config.numerics.hydro_dual_energy_eta;
  return policy;
}

[[nodiscard]] hydro::HydroSlopeLimiter hydroSlopeLimiterFromConfig(const core::SimulationConfig& config) {
  switch (config.numerics.hydro_slope_limiter) {
    case core::HydroSlopeLimiterMode::kMinmod:
      return hydro::HydroSlopeLimiter::kMinmod;
    case core::HydroSlopeLimiterMode::kMonotonizedCentral:
      return hydro::HydroSlopeLimiter::kMonotonizedCentral;
    case core::HydroSlopeLimiterMode::kVanLeer:
      return hydro::HydroSlopeLimiter::kVanLeer;
  }
  return hydro::HydroSlopeLimiter::kMonotonizedCentral;
}

[[nodiscard]] hydro::HydroReconstructionScheme hydroReconstructionSchemeFromConfig(
    const core::SimulationConfig& config) {
  switch (config.numerics.hydro_reconstruction) {
    case core::HydroReconstructionMode::kPiecewiseConstant:
      return hydro::HydroReconstructionScheme::kPiecewiseConstant;
    case core::HydroReconstructionMode::kPiecewiseLinear:
      return hydro::HydroReconstructionScheme::kPiecewiseLinear;
    case core::HydroReconstructionMode::kPiecewiseParabolic:
      return hydro::HydroReconstructionScheme::kPiecewiseParabolic;
  }
  return hydro::HydroReconstructionScheme::kPiecewiseLinear;
}

[[nodiscard]] hydro::HydroRiemannSolverKind hydroRiemannSolverKindFromConfig(
    const core::SimulationConfig& config) {
  return config.numerics.hydro_riemann_solver == core::HydroRiemannSolverMode::kHlle
      ? hydro::HydroRiemannSolverKind::kHlle
      : hydro::HydroRiemannSolverKind::kHllc;
}

// Temperature cache conversion T[K] = factor * e[code], fully ionized
// primordial gas (mu = 0.59), matching the cooling module's convention
// (CoolingModelConfig::temperature_per_internal_energy_k). One rule for the
// fixed-grid store, the AMR scatter/reflux and the feedback refresh.
[[nodiscard]] double temperatureKelvinPerInternalEnergyCode(const core::UnitSystem& units) {
  return hydro::temperatureKelvinPerInternalEnergyCode(k_gamma_adiabatic, units.velocity_si_per_code);
}

// Build the selected reconstruction for one advance. courant_factors must be
// dt/(a dx_axis) at the flux-update epoch (hydroTransportCourantFactors).
[[nodiscard]] std::unique_ptr<hydro::HydroReconstruction> makeConfiguredReconstruction(
    const core::SimulationConfig& config,
    const std::array<double, 3>& courant_factors) {
  const hydro::HydroReconstructionPolicy policy{
      .limiter = hydroSlopeLimiterFromConfig(config),
      .dt_over_dx_code = courant_factors[0],
      .dt_over_cell_width_code = courant_factors,
      .rho_floor = config.numerics.hydro_density_floor_code > 0.0
          ? config.numerics.hydro_density_floor_code : 1.0e-30,
      .pressure_floor = config.numerics.hydro_pressure_floor_code > 0.0
          ? config.numerics.hydro_pressure_floor_code : 1.0e-30,
      .enable_muscl_hancock_predictor = true,
      .adiabatic_index = k_gamma_adiabatic,
  };
  switch (hydroReconstructionSchemeFromConfig(config)) {
    case hydro::HydroReconstructionScheme::kPiecewiseConstant:
      return std::make_unique<hydro::PiecewiseConstantReconstruction>();
    case hydro::HydroReconstructionScheme::kPiecewiseLinear:
      return std::make_unique<hydro::MusclHancockReconstruction>(policy);
    case hydro::HydroReconstructionScheme::kPiecewiseParabolic:
      return std::make_unique<hydro::PpmReconstruction>(policy);
  }
  throw std::invalid_argument("unhandled hydro reconstruction scheme");
}

class SimulationStateHydroSourcePropertyProvider final
    : public hydro::HydroCellSourcePropertyProvider {
 public:
  SimulationStateHydroSourcePropertyProvider(
      const core::SimulationState& state,
      const core::UnitSystem& units,
      double density_floor_code,
      double density_frame_factor,
      std::span<const std::uint32_t> dense_row_by_geometry_row)
      : m_state(state),
        m_units(units),
        m_density_floor_code(density_floor_code),
        m_density_frame_factor(density_frame_factor),
        m_dense_row_by_geometry_row(dense_row_by_geometry_row) {}

  // HydroCoreSolver passes SOLVER GEOMETRY rows. CHUI Part 2 audit repair
  // (supplementary risk "cooling geometry row permutation"): the previous
  // provider indexed SimulationState dense rows with the geometry row, so any
  // layout whose geometry order differed from dense storage order cooled each
  // cell with another cell's density/temperature/metallicity. The geometry
  // row is now translated through the authoritative layout map.
  [[nodiscard]] hydro::HydroCellSourceProperties sourcePropertiesForCell(
      std::size_t geometry_row) const override {
    if (geometry_row >= m_dense_row_by_geometry_row.size()) {
      throw std::out_of_range(
          "SimulationStateHydroSourcePropertyProvider geometry row is outside the Cartesian layout map");
    }
    const std::size_t cell_index = m_dense_row_by_geometry_row[geometry_row];
    if (cell_index >= m_state.cells.size() || cell_index >= m_state.gas_cells.size()) {
      throw std::out_of_range(
          "SimulationStateHydroSourcePropertyProvider cell index is outside gas state");
    }
    constexpr double k_hydrogen_mass_fraction = 0.76;
    constexpr double k_proton_mass_g = 1.67262192369e-24;
    const double gas_mass_code = m_state.cells.mass_code[cell_index];
    const double metallicity = gas_mass_code > 0.0
        ? std::clamp(
              m_state.gas_cells.metal_mass_code[cell_index] / gas_mass_code,
              0.0,
              1.0)
        : 0.0;
    const double rho_stored_code = std::max(
        m_state.gas_cells.density_code[cell_index], m_density_floor_code);
    const double mass_density_physical_cgs =
        m_units.densityCodeToCgs(rho_stored_code * m_density_frame_factor);
    return hydro::HydroCellSourceProperties{
        .mass_density_physical_cgs = mass_density_physical_cgs,
        .hydrogen_number_density_cgs =
            k_hydrogen_mass_fraction * mass_density_physical_cgs / k_proton_mass_g,
        .metallicity_mass_fraction = metallicity,
        .temperature_k = m_state.gas_cells.temperature_code[cell_index]};
  }

 private:
  const core::SimulationState& m_state;
  const core::UnitSystem& m_units;
  double m_density_floor_code = 0.0;
  double m_density_frame_factor = 1.0;
  std::span<const std::uint32_t> m_dense_row_by_geometry_row;
};
[[nodiscard]] internal::CartesianGasCellRowLayout requireCartesianGasCellRowLayout(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    std::string_view caller) {
  internal::CartesianGasCellLayoutBuildResult result =
      internal::buildCartesianGasCellRowLayout(state, config);
  if (!result.ok()) {
    throw std::runtime_error(
        std::string(caller) + ": fixed Cartesian hydro geometry rejected: " + result.diagnostic);
  }
  return std::move(result.layout);
}

struct HydroGhostConservedSnapshot {
  std::vector<std::uint32_t> cell_indices;
  std::vector<std::uint64_t> gas_cell_ids;
  std::vector<std::uint64_t> parent_particle_ids;
  std::vector<std::uint32_t> owner_ranks;
  std::vector<hydro::HydroConservedState> conserved_state;
};

struct HydroConservativeGhostSyncReport {
  std::size_t restored_ghost_cells = 0;
  double rejected_remote_delta_l1 = 0.0;
  std::vector<parallel::HydroConservativeFluxCorrectionRecord> correction_records;
};

struct HydroRemoteGhostBoundaryReport {
  std::size_t boundary_metadata_sent = 0;
  std::size_t boundary_metadata_received = 0;
  std::size_t requested_cells = 0;
  std::size_t received_cells = 0;
  std::size_t imported_ghosts = 0;
  std::size_t interface_faces = 0;
  std::size_t stale_or_invalid_payloads = 0;
  std::uint64_t request_bytes = 0;
  std::uint64_t payload_bytes = 0;
};

struct HydroFloorApplicationDiagnostics {
  std::uint64_t density_floor_cell_count = 0U;
  std::uint64_t pressure_floor_cell_count = 0U;
  double density_floor_correction_l1_code = 0.0;
  double pressure_floor_correction_l1_code = 0.0;
};

struct HydroBoundaryCellAdvertisement {
  std::uint64_t gas_cell_id = 0;
  int owner_rank = 0;
  std::uint64_t hydro_sync_epoch = 0;
  std::uint64_t decomposition_epoch = 0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
};

[[nodiscard]] bool finiteBoundaryAdvertisement(const HydroBoundaryCellAdvertisement& record) {
  return record.gas_cell_id != 0U &&
      record.owner_rank >= 0 &&
      std::isfinite(record.center_x_comoving) &&
      std::isfinite(record.center_y_comoving) &&
      std::isfinite(record.center_z_comoving);
}

[[nodiscard]] std::vector<HydroBoundaryCellAdvertisement> exchangeHydroBoundaryAdvertisements(
    const parallel::MpiContext& mpi_context,
    std::span<const HydroBoundaryCellAdvertisement> local_records) {
  std::exception_ptr local_failure;
  try {
    for (const HydroBoundaryCellAdvertisement& record : local_records) {
      if (!finiteBoundaryAdvertisement(record)) {
        throw std::invalid_argument("hydro boundary advertisement is invalid");
      }
      if (record.owner_rank != mpi_context.worldRank()) {
        throw std::invalid_argument("hydro boundary advertisement owner rank does not match MPI context");
      }
    }
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "hydro boundary advertisement local validation");
  if (!mpi_context.isEnabled()) {
    return std::vector<HydroBoundaryCellAdvertisement>(local_records.begin(), local_records.end());
  }
#if COSMOSIM_ENABLE_MPI
  std::span<const std::uint8_t> local_bytes;
  local_failure = nullptr;
  try {
    const std::size_t byte_count = core::checkedSizeMultiply(
        local_records.size(), sizeof(HydroBoundaryCellAdvertisement),
        "hydro boundary advertisement local byte count");
    local_bytes = std::span<const std::uint8_t>(
        reinterpret_cast<const std::uint8_t*>(local_records.data()), byte_count);
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "hydro boundary advertisement byte framing");

  const std::vector<std::uint8_t> recv_bytes =
      mpi_context.allgatherBytesBounded(local_bytes);
  std::vector<HydroBoundaryCellAdvertisement> result;
  local_failure = nullptr;
  try {
    if (recv_bytes.size() % sizeof(HydroBoundaryCellAdvertisement) != 0U) {
      throw std::runtime_error(
          "hydro boundary advertisement exchange returned a partial record");
    }
    result.resize(recv_bytes.size() / sizeof(HydroBoundaryCellAdvertisement));
    if (!recv_bytes.empty()) {
      std::memcpy(result.data(), recv_bytes.data(), recv_bytes.size());
    }
    for (const HydroBoundaryCellAdvertisement& record : result) {
      if (!finiteBoundaryAdvertisement(record) ||
          record.owner_rank >= mpi_context.worldSize()) {
        throw std::runtime_error(
            "hydro boundary advertisement exchange returned invalid metadata");
      }
    }
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "hydro boundary advertisement reassembly validation");
  return result;
#else
  throw std::runtime_error("hydro boundary advertisement exchange requires MPI support when MPI context is enabled");
#endif
}

[[nodiscard]] std::uint8_t hydroAxisCode(hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return 0U;
    case hydro::HydroFaceAxis::kY:
      return 1U;
    case hydro::HydroFaceAxis::kZ:
      return 2U;
  }
  return 0U;
}

[[nodiscard]] std::uint8_t hydroSideCode(hydro::HydroFaceSide side) {
  return side == hydro::HydroFaceSide::kLower ? 0U : 1U;
}

// Wire conversion of one flux-register half (shared by the synchronized and
// recursive distributed paths).
[[nodiscard]] parallel::AmrFluxRegisterPayloadRecord fluxRegisterPayloadFromEntry(
    const amr::FluxRegisterEntry& entry, int source_rank, int owner_rank, const core::SimulationState& state) {
    parallel::AmrFluxRegisterPayloadRecord record;
    record.register_key = entry.register_key;
    record.coarse_patch_id = entry.coarse_patch_id;
    record.coarse_gas_cell_id = entry.coarse_gas_cell_id;
    record.coarse_cell_index = static_cast<std::uint64_t>(entry.coarse_cell_index);
    record.level = entry.level;
    record.axis = hydroAxisCode(entry.axis);
    record.orientation = hydroSideCode(entry.orientation);
    record.source_rank = source_rank;
    record.owner_rank = owner_rank;
    record.gas_cell_identity_generation = state.gasCellIdentityGeneration();
    record.patch_geometry_generation = state.cellIndexGeneration();
    record.coarse_mass_flux_code = entry.coarse_face_flux_code.mass_code;
    record.coarse_momentum_x_flux_code = entry.coarse_face_flux_code.momentum_x_code;
    record.coarse_momentum_y_flux_code = entry.coarse_face_flux_code.momentum_y_code;
    record.coarse_momentum_z_flux_code = entry.coarse_face_flux_code.momentum_z_code;
    record.coarse_total_energy_flux_code = entry.coarse_face_flux_code.total_energy_code;
    record.coarse_metal_mass_flux_code = entry.coarse_face_flux_code.metal_mass_code;
    record.fine_mass_flux_code = entry.fine_face_flux_code.mass_code;
    record.fine_momentum_x_flux_code = entry.fine_face_flux_code.momentum_x_code;
    record.fine_momentum_y_flux_code = entry.fine_face_flux_code.momentum_y_code;
    record.fine_momentum_z_flux_code = entry.fine_face_flux_code.momentum_z_code;
    record.fine_total_energy_flux_code = entry.fine_face_flux_code.total_energy_code;
    record.fine_metal_mass_flux_code = entry.fine_face_flux_code.metal_mass_code;
    record.face_area_comov = entry.face_area_comov;
    record.coarse_area_comov = entry.coarse_area_comov;
    record.fine_area_comov = entry.fine_area_comov;
    record.dt_code = entry.dt_code;
    record.coarse_face_count = entry.coarse_face_count;
    record.fine_face_count = entry.fine_face_count;
  parallel::validateAmrFluxRegisterPayloadRecord(record);
  return record;
}

[[nodiscard]] amr::FluxRegisterEntry fluxRegisterEntryFromPayload(
    const parallel::AmrFluxRegisterPayloadRecord& payload) {
  return amr::FluxRegisterEntry{
        .register_key = payload.register_key,
        .coarse_patch_id = payload.coarse_patch_id,
        .coarse_gas_cell_id = payload.coarse_gas_cell_id,
        .coarse_cell_index = static_cast<std::size_t>(payload.coarse_cell_index),
        .level = payload.level,
        .axis = payload.axis == 0U ? hydro::HydroFaceAxis::kX : (payload.axis == 1U ? hydro::HydroFaceAxis::kY : hydro::HydroFaceAxis::kZ),
        .orientation = payload.orientation == 0U ? hydro::HydroFaceSide::kLower : hydro::HydroFaceSide::kUpper,
        .coarse_face_flux_code = amr::ConservedState{
            .mass_code = payload.coarse_mass_flux_code,
            .momentum_x_code = payload.coarse_momentum_x_flux_code,
            .momentum_y_code = payload.coarse_momentum_y_flux_code,
            .momentum_z_code = payload.coarse_momentum_z_flux_code,
            .total_energy_code = payload.coarse_total_energy_flux_code,
            .metal_mass_code = payload.coarse_metal_mass_flux_code},
        .fine_face_flux_code = amr::ConservedState{
            .mass_code = payload.fine_mass_flux_code,
            .momentum_x_code = payload.fine_momentum_x_flux_code,
            .momentum_y_code = payload.fine_momentum_y_flux_code,
            .momentum_z_code = payload.fine_momentum_z_flux_code,
            .total_energy_code = payload.fine_total_energy_flux_code,
            .metal_mass_code = payload.fine_metal_mass_flux_code},
        .face_area_comov = payload.face_area_comov,
        .coarse_area_comov = payload.coarse_area_comov,
        .fine_area_comov = payload.fine_area_comov,
        .dt_code = payload.dt_code,
        .coarse_face_count = payload.coarse_face_count,
        .fine_face_count = payload.fine_face_count};
}

[[nodiscard]] std::uint64_t hydroInterfaceFaceKey(
    std::uint64_t local_gas_cell_id,
    std::uint64_t remote_gas_cell_id,
    hydro::HydroFaceAxis axis) {
  const std::uint64_t lo = std::min(local_gas_cell_id, remote_gas_cell_id);
  const std::uint64_t hi = std::max(local_gas_cell_id, remote_gas_cell_id);
  std::uint64_t key = 1469598103934665603ULL;
  auto mix = [&](std::uint64_t value) {
    key ^= value;
    key *= 1099511628211ULL;
  };
  mix(lo);
  mix(hi);
  mix(static_cast<std::uint64_t>(hydroAxisCode(axis) + 1U));
  return key == 0U ? 1U : key;
}

[[nodiscard]] bool nearlySameCoordinate(double lhs, double rhs, double tolerance) {
  return std::abs(lhs - rhs) <= tolerance;
}

[[nodiscard]] HydroGhostConservedSnapshot snapshotHydroGhostConservedCells(
    const core::SimulationState& state,
    const hydro::HydroConservedStateSoa& conserved,
    const hydro::HydroPatchGeometry& geometry,
    std::span<const std::uint32_t> geometry_row_by_dense_row,
    std::uint32_t world_rank) {
  HydroGhostConservedSnapshot snapshot;
  state.requireGasCellIdentityMapCoversDenseRows("snapshotHydroGhostConservedCells");
  const auto particle_row_by_id = internal::buildParticleRowById(state);
  snapshot.cell_indices.reserve(state.cells.size());
  snapshot.gas_cell_ids.reserve(state.cells.size());
  snapshot.parent_particle_ids.reserve(state.cells.size());
  snapshot.owner_ranks.reserve(state.cells.size());
  snapshot.conserved_state.reserve(state.cells.size());
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint32_t owner_rank = internal::gasCellOwnerRankForLocalRow(
        state,
        cell_index,
        particle_row_by_id,
        "snapshotHydroGhostConservedCells");
    if (owner_rank == world_rank) {
      continue;
    }
    const core::GasCellIdentityRecord& identity = internal::gasCellIdentityRecordForLocalRow(
        state,
        cell_index,
        "snapshotHydroGhostConservedCells");
    snapshot.cell_indices.push_back(cell_index);
    snapshot.gas_cell_ids.push_back(identity.gas_cell_id);
    snapshot.parent_particle_ids.push_back(identity.parent_particle_id.value_or(0U));
    snapshot.owner_ranks.push_back(owner_rank);
    if (cell_index >= geometry_row_by_dense_row.size()) {
      throw std::out_of_range("snapshotHydroGhostConservedCells: dense row is outside Cartesian geometry map");
    }
    snapshot.conserved_state.push_back(conserved.loadCell(geometry_row_by_dense_row[cell_index]));
  }
  for (const hydro::HydroGhostCell& ghost : geometry.ghost_cells) {
    if (ghost.boundary_kind != hydro::HydroBoundaryKind::kImportedMpi) {
      continue;
    }
    if (ghost.origin_rank < 0 || ghost.origin_rank == static_cast<int>(world_rank) ||
        ghost.origin_gas_cell_id == 0U || ghost.ghost_cell >= conserved.size()) {
      throw std::runtime_error("snapshotHydroGhostConservedCells: imported MPI ghost metadata is invalid");
    }
    snapshot.cell_indices.push_back(static_cast<std::uint32_t>(ghost.ghost_cell));
    snapshot.gas_cell_ids.push_back(ghost.origin_gas_cell_id);
    snapshot.parent_particle_ids.push_back(0U);
    snapshot.owner_ranks.push_back(static_cast<std::uint32_t>(ghost.origin_rank));
    snapshot.conserved_state.push_back(conserved.loadCell(ghost.ghost_cell));
  }
  return snapshot;
}

[[nodiscard]] HydroConservativeGhostSyncReport restoreHydroGhostConservedCells(
    hydro::HydroConservedStateSoa& conserved,
    const HydroGhostConservedSnapshot& snapshot,
    std::uint32_t world_rank) {
  if (snapshot.cell_indices.size() != snapshot.conserved_state.size() ||
      snapshot.cell_indices.size() != snapshot.gas_cell_ids.size() ||
      snapshot.cell_indices.size() != snapshot.parent_particle_ids.size() ||
      snapshot.cell_indices.size() != snapshot.owner_ranks.size()) {
    throw std::invalid_argument("hydro ghost snapshot is inconsistent");
  }
  HydroConservativeGhostSyncReport report;
  report.correction_records.reserve(snapshot.cell_indices.size());
  for (std::size_t i = 0; i < snapshot.cell_indices.size(); ++i) {
    const std::uint32_t cell_index = snapshot.cell_indices[i];
    if (cell_index >= conserved.size()) {
      throw std::out_of_range("hydro ghost snapshot cell index is outside conserved state");
    }
    const hydro::HydroConservedState before = snapshot.conserved_state[i];
    const hydro::HydroConservedState after = conserved.loadCell(cell_index);
    const double dm = after.mass_density_comoving - before.mass_density_comoving;
    const double dmx = after.momentum_density_x_comoving - before.momentum_density_x_comoving;
    const double dmy = after.momentum_density_y_comoving - before.momentum_density_y_comoving;
    const double dmz = after.momentum_density_z_comoving - before.momentum_density_z_comoving;
    const double de = after.total_energy_density_comoving - before.total_energy_density_comoving;
    const double dmz_metal =
        after.metal_mass_density_comoving - before.metal_mass_density_comoving;
    report.rejected_remote_delta_l1 += std::abs(dm);
    report.rejected_remote_delta_l1 += std::abs(dmx);
    report.rejected_remote_delta_l1 += std::abs(dmy);
    report.rejected_remote_delta_l1 += std::abs(dmz);
    report.rejected_remote_delta_l1 += std::abs(de);
    report.rejected_remote_delta_l1 += std::abs(dmz_metal);
    if (dm != 0.0 || dmx != 0.0 || dmy != 0.0 || dmz != 0.0 || de != 0.0 ||
        dmz_metal != 0.0) {
      parallel::HydroConservativeFluxCorrectionRecord record;
      record.gas_cell_id = snapshot.gas_cell_ids[i];
      record.parent_particle_id = snapshot.parent_particle_ids[i];
      record.source_rank = static_cast<int>(world_rank);
      record.owner_rank = static_cast<int>(snapshot.owner_ranks[i]);
      record.delta_mass_density_comoving = dm;
      record.delta_momentum_density_x_comoving = dmx;
      record.delta_momentum_density_y_comoving = dmy;
      record.delta_momentum_density_z_comoving = dmz;
      record.delta_total_energy_density_comoving = de;
      record.delta_metal_mass_density_comoving = dmz_metal;
      parallel::validateHydroConservativeFluxCorrectionRecord(record);
      report.correction_records.push_back(record);
    }
    conserved.storeCell(cell_index, before);
    ++report.restored_ghost_cells;
  }
  return report;
}



// T02 bridge between the AMR library's level-clock seam and the workflow
// source owner. Owns no state: it forwards each level half to
// SourceRuntime::executeAmrLevelSmoothSources() with the coarse StepContext
// that owns the accepted interval. Lifetime: one hydro stage (stack object);
// both references outlive it.
class RecursiveAmrLevelSourceAdapter final : public amr::RecursiveAmrLevelSourceOperator {
 public:
  RecursiveAmrLevelSourceAdapter(SourceRuntime& source, core::StepContext& step) noexcept
      : m_source(source), m_step(step) {}

  void apply(core::SimulationState& state,
             std::span<const std::uint32_t> level_rows,
             const hydro::HydroUpdateContext& level_update,
             int level,
             amr::RecursiveAmrSourcePhase phase) override {
    if (&state != &m_step.state) {
      throw std::logic_error(
          "recursive AMR level source adapter received a foreign SimulationState");
    }
    const core::SourceSplitPhase split_phase =
        phase == amr::RecursiveAmrSourcePhase::kPreHydroHalf
        ? core::SourceSplitPhase::kStrangPreHalf
        : core::SourceSplitPhase::kStrangPostHalf;
    m_source.executeAmrLevelSmoothSources(m_step, level_rows, level_update, level, split_phase);
  }

 private:
  SourceRuntime& m_source;
  core::StepContext& m_step;
};

class HydroAmrRuntimeImpl final : public HydroAmrRuntime {
 public:
  HydroAmrRuntimeImpl(
      const core::SimulationConfig& config,
      const core::ModePolicy& mode_policy,
      const RuntimeServices& services,
      std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table)
      : m_config(config),
        m_mode_policy(mode_policy),
        m_mpi_context(services.mpi_context),
        m_memory_governor(services.memory_governor),
        m_units(core::makeUnitSystem(
            config.units.length_unit, config.units.mass_unit, config.units.velocity_unit)),
        m_temperature_k_per_internal_energy_code(temperatureKelvinPerInternalEnergyCode(m_units)),
        m_solver(k_gamma_adiabatic, hydroNumericsPolicyFromConfig(config)),
        m_reconstruction(makeConfiguredReconstruction(config, {0.0, 0.0, 0.0})),
        m_riemann_solver(hydro::makeHydroRiemannSolver(hydroRiemannSolverKindFromConfig(config))) {
    // CHUI Part 2: one owner per physical source. A Part 2 thermochemistry
    // backend is advanced by the SourceRuntime thermochemistry stage; the
    // legacy in-solver cooling source is then NOT registered (no double
    // cooling). legacy_cooling keeps the (F2-repaired) in-solver path.
    if (m_config.physics.enable_cooling && !physics::part2ThermochemistrySelected(m_config.physics)) {
      constexpr double k_mean_molecular_weight_ionized = 0.59;
      constexpr double k_proton_mass_kg = 1.67262192369e-27;
      constexpr double k_boltzmann_j_per_k = 1.380649e-23;
      const double velocity_sq_si =
          m_units.velocity_si_per_code * m_units.velocity_si_per_code;
      physics::CoolingModelConfig cooling_config{
          .uv_background_model = m_config.physics.uv_background_model,
          .self_shielding_model = m_config.physics.self_shielding_model,
          .temperature_floor_k = m_config.physics.temperature_floor_k,
          .primordial_helium_mass_fraction = 0.24,
          .enable_metal_line_cooling =
              m_config.physics.cooling_model == core::CoolingModel::kPrimordialMetalLine,
          .max_fractional_energy_change_per_substep = 0.1,
          .max_subcycles = 32,
          .temperature_per_internal_energy_k =
              (k_gamma_adiabatic - 1.0) * k_mean_molecular_weight_ionized *
              k_proton_mass_kg / k_boltzmann_j_per_k * velocity_sq_si,
      };
      physics::CoolingRateProvider provider(cooling_config);
      if (cooling_config.enable_metal_line_cooling) {
        if (m_config.physics.metal_line_table_path.empty()) {
          throw std::invalid_argument(
              "metal-line cooling is enabled but physics.metal_line_table_path is empty");
        }
        provider.setMetalLineTable(physics::MetalLineCoolingTable::loadFromTextFile(
            m_config.physics.metal_line_table_path));
      }
      const double temperature_energy_floor = std::max(
          m_config.physics.temperature_floor_k /
              cooling_config.temperature_per_internal_energy_k,
          1.0e-30);
      // Cooling rates are erg g^-1 s^-1 after division by physical CGS density.
      // 1 erg/g = 1e-4 J/kg = 1e-4 m^2/s^2.
      const double specific_energy_rate_cgs_to_code =
          1.0e-4 * m_units.timeSiPerCode() / velocity_sq_si;
      // R1: the same rate provider evaluates per-cell cooling times for the
      // mandatory cooling-length refinement criterion (read-only copy).
      m_cooling_rate_provider.emplace(provider);
      m_cooling_specific_rate_cgs_to_code = specific_energy_rate_cgs_to_code;
      m_cooling_temperature_per_energy_k = cooling_config.temperature_per_internal_energy_k;
      m_cooling_helium_mass_fraction = cooling_config.primordial_helium_mass_fraction;
      m_cooling_source = std::make_unique<physics::CoolingHeatingSource>(
          std::move(provider),
          physics::CoolingSourceIntegrator(
              temperature_energy_floor, specific_energy_rate_cgs_to_code));
    }

    if (m_config.physics.enable_star_formation &&
        m_config.physics.star_formation_model ==
            core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
      if (effective_eos_table == nullptr) {
        core::UnitSystem units = core::makeUnitSystem(
            m_config.units.length_unit,
            m_config.units.mass_unit,
            m_config.units.velocity_unit);
        effective_eos_table =
            std::make_shared<const physics::EffectiveMultiphaseEosTable>(
                physics::makeEffectiveMultiphaseEosConfig(m_config.physics),
                std::move(units),
                physics::makeEffectiveIsmReferenceCoolingProvider(m_config.physics));
      }
      m_effective_ism_closure =
          std::make_unique<physics::EffectiveIsmThermodynamicClosure>(
              std::move(effective_eos_table));
      m_effective_ism_energy_source =
          std::make_unique<physics::EffectiveIsmEnergyRelaxationSource>(*m_effective_ism_closure);
    }
  }

  [[nodiscard]] const hydro::HydroProfileEvent& lastHydroProfile() const noexcept { return m_last_hydro_profile; }
  [[nodiscard]] const internal::SolverGhostRefreshReport&
  lastGhostRefreshReport() const noexcept {
    return m_last_ghost_refresh;
  }
  [[nodiscard]] const HydroRemoteGhostBoundaryReport& lastRemoteGhostBoundaryReport() const noexcept {
    return m_last_remote_ghost_boundary_report;
  }
  [[nodiscard]] const core::HydroCflDiagnostics& lastHydroCflDiagnostics() const noexcept {
    return m_last_hydro_cfl_diagnostics;
  }
  [[nodiscard]] std::uint64_t ghostExchangeBytesRecent() const noexcept override {
    return m_last_ghost_refresh.sent_bytes + m_last_ghost_refresh.received_bytes;
  }
  [[nodiscard]] std::size_t remoteImportedGhostCount() const noexcept override {
    return m_last_remote_ghost_boundary_report.imported_ghosts;
  }
  [[nodiscard]] std::size_t remoteInterfaceFaceCount() const noexcept override {
    return m_last_remote_ghost_boundary_report.interface_faces;
  }
  [[nodiscard]] std::size_t remoteStaleInvalidPayloadCount() const noexcept override {
    return m_last_remote_ghost_boundary_report.stale_or_invalid_payloads;
  }
  [[nodiscard]] std::uint64_t directedAmrRemoteInterfaceCount() const noexcept override {
    return m_last_directed_amr_diagnostics.remote_interface_count;
  }
  [[nodiscard]] std::uint64_t directedAmrPatchCellRecordsReceived() const noexcept override {
    return m_last_directed_amr_diagnostics.directed_patch_cell_records_received;
  }
  [[nodiscard]] core::MemoryReport memoryReport() const override {
    core::MemoryReportBuilder builder;
    const std::uint64_t conserved_capacity = m_conserved.ownedCapacityBytes();
    if (conserved_capacity != 0U) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kGasHydro,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kPhaseResident,
          .label = "hydro.runtime.fixed_patch_conserved_workspace",
          .current_size_bytes = m_conserved.currentSizeBytes(),
          .owned_capacity_bytes = conserved_capacity,
          .high_water_bytes = conserved_capacity});
    }
    const std::uint64_t scratch_capacity = m_scratch.ownedCapacityBytes();
    if (scratch_capacity != 0U || m_last_amr_hydro_diagnostics.scratch_high_water_bytes != 0U) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kScratch,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kScratchArena,
          .label = "hydro.runtime.active_batch_scratch",
          .current_size_bytes = scratch_capacity,
          .owned_capacity_bytes = scratch_capacity,
          .high_water_bytes = std::max(
              m_scratch.high_water.capacity_bytes,
              m_last_amr_hydro_diagnostics.scratch_high_water_bytes),
          .governed_commitment = false,
          .uncertainty_note = "AMR per-patch scratch is released after each patch; high-water retains the observed peak"});
    }
    if (m_last_amr_hydro_diagnostics.geometry_capacity_bytes != 0U) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kGasHydro,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kPhaseResident,
          .label = "hydro.amr.geometry_retained",
          .high_water_bytes = m_last_amr_hydro_diagnostics.geometry_capacity_bytes,
          .estimate_only = false,
          .uncertainty_note = "actual nested capacity; parent object-vector storage and admission are tracked by the stage owner; released at stage exit"});
    }
    if (m_last_amr_hydro_diagnostics.prepared_ghost_capacity_bytes != 0U) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kMpiBuffers,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kCommunication,
          .label = "hydro.amr.step_start_ghost_snapshots",
          .high_water_bytes = m_last_amr_hydro_diagnostics.prepared_ghost_capacity_bytes,
          .estimate_only = true,
          .uncertainty_note = "surface-scaled synchronization snapshot released at hydro-phase end"});
    }
    if (m_last_amr_hydro_diagnostics.max_patch_conserved_bytes != 0U) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kGasHydro,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kPhaseResident,
          .label = "hydro.amr.max_patch_conserved_workspace",
          .high_water_bytes = m_last_amr_hydro_diagnostics.max_patch_conserved_bytes,
          .estimate_only = true,
          .uncertainty_note = "one-patch-at-a-time conserved staging; released before the next phase"});
    }
    const std::uint64_t active_metadata_capacity = core::checkedMemoryBytesAdd(
        core::ownedCapacityBytesForContainer(m_active_cells),
        core::ownedCapacityBytesForContainer(m_active_faces),
        "hydro active metadata capacity");
    if (active_metadata_capacity != 0U) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kActiveSets,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kPhaseResident,
          .label = "hydro.runtime.active_index_metadata",
          .current_size_bytes = core::checkedMemoryBytesAdd(
              core::currentSizeBytesForContainer(m_active_cells),
              core::currentSizeBytesForContainer(m_active_faces),
              "hydro active metadata logical bytes"),
          .owned_capacity_bytes = active_metadata_capacity,
          .high_water_bytes = active_metadata_capacity});
    }
    if (m_communication_high_water_bytes != 0U) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kMpiBuffers,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kCommunication,
          .label = "hydro.communication.total_resident_high_water",
          .high_water_bytes = m_communication_high_water_bytes,
          .estimate_only = false,
          .uncertainty_note = "capacity-based simultaneous directed-AMR transport plus sparse remote-ghost ownership; traffic volume is reported separately"});
    }
    if (m_remote_ghost_capacity_high_water_bytes != 0U) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kMpiBuffers,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kCommunication,
          .label = "hydro.remote_ghost.sparse_cell_capacity",
          .high_water_bytes = m_remote_ghost_capacity_high_water_bytes,
          .estimate_only = false,
          .uncertainty_note = "compact boundary-cell store; residency scales with advertised remote interface cells"});
    }
    if (m_remote_ghost_metadata_high_water_bytes != 0U) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kMpiBuffers,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kCommunication,
          .label = "hydro.remote_ghost.patch_metadata",
          .high_water_bytes = m_remote_ghost_metadata_high_water_bytes,
          .estimate_only = false,
          .uncertainty_note = "O(remote patches) sparse ghost descriptors; no remote solver geometry or full-volume SoA"});
    }
    return std::move(builder).finish();
  }

  void bindSourceRuntime(SourceRuntime* source_runtime) noexcept override {
    if (m_source_runtime != nullptr && m_source_runtime != source_runtime) {
      m_source_runtime->setAmrLevelSourceCouplingEnabled(false);
    }
    m_source_runtime = source_runtime;
    // T02 composition contract: under recursive subcycling the AMR level
    // clocks own the smooth baryonic halves. The same predicate installs the
    // RecursiveAmrLevelSourceAdapter in runProductionAmrHydroPath(), so the
    // source owner's global smooth halves are skipped iff levels run them.
    if (m_source_runtime != nullptr) {
      m_source_runtime->setAmrLevelSourceCouplingEnabled(
          m_config.numerics.amr_time_stepping ==
          core::AmrTimeSteppingConfigMode::kRecursiveSubcycling);
    }
  }

  [[nodiscard]] RuntimeTaskMemoryEstimate estimateMemory(
      const core::SimulationState& state) const override {
    if (state.cells.size() == 0U && state.patches.size() == 0U) {
      return {0U, true, {}};
    }
    const std::uint64_t batch = std::min<std::uint64_t>(
        static_cast<std::uint64_t>(state.cells.size()),
        m_config.numerics.hydro_active_batch_max_cells == 0U
            ? hydro::k_hydro_automatic_active_batch_max_cells
            : m_config.numerics.hydro_active_batch_max_cells);
    const std::uint64_t width = hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell;
    if (width != 0U && batch > std::numeric_limits<std::uint64_t>::max() / width) {
      throw std::overflow_error("hydro task scratch estimate overflows uint64");
    }
    return {batch * width, false,
            "active-batch scratch modeled; AMR transaction, ghost and geometry peaks require owner-local admission"};
  }

  void execute(HydroAmrStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::hydroAmrContext(
        view,
        {{RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kSourceMutationState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroPrimitiveState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kAmrPatchState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kEffectiveIsmThermodynamics, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead}});
    if (context.stage != core::IntegrationStage::kHydroUpdate) {
      throw std::logic_error("hydro handler received an unregistered stage");
    }
    const auto hydro_cell_rank_count_after_refresh = [&]() {
      const std::uint64_t local_has_hydro_cells = context.state.cells.size() == 0 ? 0ULL : 1ULL;
      return m_mpi_context.allreduceSumUint64(local_has_hydro_cells);
    };
    m_current_world_rank = static_cast<std::uint32_t>(std::max(m_mpi_context.worldRank(), 0));
    // CHUI Part 2: composition must exist before the first transport.
    (void)physics::ensureGasCompositionInitialized(context.state, m_config.physics);
    m_last_ghost_refresh = internal::refreshParticleGhostsForSolver(
        context, m_mpi_context, "hydro.godunov", &m_ghost_cache_lifecycle);
    // sent_bytes/received_bytes are traffic counters, not resident-memory
    // capacity. Do not fold them into the memory high-water metric.
    const std::uint64_t hydro_cell_rank_count = hydro_cell_rank_count_after_refresh();
    const std::uint64_t local_has_production_amr =
        amr::hasProductionAmrHydroCoverage(context.state) ? 1ULL : 0ULL;
    const std::uint64_t production_amr_rank_count =
        m_mpi_context.allreduceSumUint64(local_has_production_amr);
    const bool all_ranks_have_hydro_cells =
        !m_mpi_context.isEnabled() ||
        hydro_cell_rank_count == static_cast<std::uint64_t>(m_mpi_context.worldSize());
    m_last_hydro_profile = {};
    m_last_amr_hydro_diagnostics = {};
    m_last_floor_diagnostics = auditConfiguredHydroFloors(context.state);

    // Production AMR exchanges are collective.  Ranks with zero local cells
    // must still enter the control-plane and zero-record payload exchanges
    // whenever any peer owns a production patch.
    if (production_amr_rank_count > 0U) {
      if (context.state.cells.size() != 0U) {
        context.state.requireGasCellIdentityMapCoversDenseRows("hydro callback");
      }
      runProductionAmrHydroPath(context);
      persistEffectiveIsmDiagnostics(context.state);
      // Every rank advances the same authoritative gravity-source epoch when
      // any rank owns production gas. This covers gas mass/leaf-topology
      // mutations without an O(N) source fingerprint.
      context.state.bumpGravitySourceGeneration();
      recordHydroFloorDiagnostics(context);
      return;
    }
    if (context.state.cells.size() == 0U) {
      persistEffectiveIsmDiagnostics(context.state);
      if (hydro_cell_rank_count > 0U) {
        context.state.bumpGravitySourceGeneration();
      }
      return;
    }

    // The legacy Cartesian fallback discovers cross-rank boundaries by globally
    // replicating boundary advertisements.  That path is useful for small
    // single-rank validation problems, but it is not a scalable production MPI
    // contract and must not be selected silently.  Distributed production hydro
    // is owned by the stable-gas-cell-ID AMR exchange above.
    if (m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1) {
      throw std::runtime_error(
          "multi-rank fixed-Cartesian hydro is validation-only and is disabled in production; "
          "provide complete production AMR hydro coverage so gas-cell-ID directed exchange is used");
    }

    context.state.requireGasCellIdentityMapCoversDenseRows("hydro callback");
    const auto particle_row_by_id = internal::buildParticleRowById(context.state);
    std::vector<std::uint8_t> parent_mirror_updated(context.state.particles.size(), 0U);

    if (m_mpi_context.isEnabled()) {
      m_cached_geometry_key.reset();
    }
    rebuildGeometryIfNeeded(context.state, context.integrator_state.dt_time_code);
    verifyAcceptedHydroCfl(context, context.active_set.cell_indices);

    m_conserved.resize(m_geometry.totalCellStorageCount());
    const bool dual_energy_enabled =
        m_solver.numericsPolicy().dual_energy.mode != hydro::HydroDualEnergyMode::kOff;
    // Re-arm (NaN) the auxiliary lane every step; rows are filled below.
    m_conserved.enableAuxiliaryInternalEnergy(dual_energy_enabled);
    for (std::size_t cell_index = 0; cell_index < context.state.cells.size(); ++cell_index) {
      const double rho = std::max(
          context.state.gas_cells.density_code[cell_index],
          m_config.numerics.hydro_density_floor_code);
      const double internal_energy_floor = m_config.numerics.hydro_pressure_floor_code > 0.0
          ? m_config.numerics.hydro_pressure_floor_code /
                ((k_gamma_adiabatic - 1.0) * std::max(rho, 1.0e-30))
          : 0.0;
      // Canonical thermal authority is e; the derived pressure cache must not
      // raise the thermal state (it may be stale or hold an effective-EOS
      // pressure, which the closure re-derives inside the solver). States that
      // never set e (legacy initial conditions) fall back to the cache.
      const double canonical_e = context.state.gas_cells.internal_energy_code[cell_index];
      const double internal_energy = std::max(
          canonical_e > 0.0
              ? canonical_e
              : context.state.gas_cells.pressure_code[cell_index] /
                    ((k_gamma_adiabatic - 1.0) * std::max(rho, 1.0e-30)),
          internal_energy_floor);
      const double pressure = std::max(
          (k_gamma_adiabatic - 1.0) * rho * internal_energy,
          m_config.numerics.hydro_pressure_floor_code);
      const double sound_speed_squared_floor = m_config.numerics.hydro_pressure_floor_code > 0.0
          ? k_gamma_adiabatic * m_config.numerics.hydro_pressure_floor_code /
                std::max(rho, 1.0e-30)
          : 0.0;
      const hydro::HydroPrimitiveState primitive{
          .rho_comoving = rho,
          .vel_x_peculiar = context.state.gas_cells.velocity_x_peculiar[cell_index],
          .vel_y_peculiar = context.state.gas_cells.velocity_y_peculiar[cell_index],
          .vel_z_peculiar = context.state.gas_cells.velocity_z_peculiar[cell_index],
          .pressure_comoving = pressure,
          .specific_internal_energy_code = internal_energy,
          .signal_speed_squared_code = std::max(
              context.state.gas_cells.sound_speed_code[cell_index] *
                  context.state.gas_cells.sound_speed_code[cell_index],
              sound_speed_squared_floor),
          .uses_effective_ism = false,
          .metallicity_mass_fraction = std::clamp(
              context.state.gas_cells.metal_mass_code[cell_index] /
                  std::max(context.state.cells.mass_code[cell_index], 1.0e-30),
              0.0,
              1.0),
      };
      if (cell_index >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro callback dense row is outside Cartesian geometry map");
      }
      m_conserved.storeCell(
          m_geometry_row_by_dense_row[cell_index],
          hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma_adiabatic));
      if (dual_energy_enabled) {
        // Auxiliary thermal lane from the canonical specific internal energy
        // (exact; never recovered from E - K).
        m_conserved.auxiliaryInternalEnergyDensityComoving()[m_geometry_row_by_dense_row[cell_index]] =
            rho * internal_energy;
      }
    }
    m_last_remote_ghost_boundary_report = refreshRemoteHydroBoundaryGhosts(
        context,
        particle_row_by_id);
    hydro::fillHydroBoundaryGhostCells(m_conserved, m_geometry, k_gamma_adiabatic);
    hydro::fillHydroBoundaryGhostAuxiliaryInternalEnergy(m_conserved, m_geometry);
    const std::uint32_t world_rank = static_cast<std::uint32_t>(m_mpi_context.worldRank());
    const HydroGhostConservedSnapshot ghost_conserved_snapshot = snapshotHydroGhostConservedCells(
        context.state, m_conserved, m_geometry, m_geometry_row_by_dense_row, world_rank);
    const hydro::HydroActiveSetView active_view = buildActiveFaceView(context.active_set.cell_indices);

    // Hydro transport/source terms own the full accepted interval. Evaluate
    // background terms at its midpoint; AMR fine substeps derive their own
    // interval/evaluation epochs instead of inheriting the coarse end state.
    const double hydro_scale_factor = context.cosmology_background != nullptr
        ? context.timeline_step.scale_factor_midpoint
        : context.integrator_state.current_scale_factor;
    const bool cosmological_comoving =
        m_mode_policy.cosmological_comoving_frame &&
        m_config.units.coordinate_frame == core::CoordinateFrame::kComoving;
    const double hubble_rate_code = cosmological_comoving && context.cosmology_background != nullptr
        ? context.cosmology_background->hubbleSi(hydro_scale_factor) *
              context.timeline_step.time_si_per_code
        : 0.0;

    hydro::HydroUpdateContext update{
        .dt_code = context.integrator_state.dt_time_code,
        .time_begin_code = context.timeline_step.time_begin_code,
        .time_end_code = context.timeline_step.time_end_code,
        .scale_factor_begin = context.timeline_step.scale_factor_begin,
        .scale_factor = std::max(1.0e-12, hydro_scale_factor),
        .scale_factor_end = context.timeline_step.scale_factor_end,
        .hubble_rate_begin_code = context.timeline_step.hubble_begin_code,
        .hubble_rate_code = hubble_rate_code,
        .hubble_rate_end_code = context.timeline_step.hubble_end_code,
        .comoving_coordinates =
            m_config.units.coordinate_frame == core::CoordinateFrame::kComoving,
    };

    // Audit F1 follow-up: the Hancock/CTU predictors use the same
    // dt/(a dx) transport factor as the conservative flux update
    // dt/(a V) sum A F, evaluated at the same epoch a = update.scale_factor.
    m_reconstruction = makeConfiguredReconstruction(
        m_config,
        hydro::hydroTransportCourantFactors(
            update.dt_code,
            update.scale_factor,
            m_geometry.cell_width_x_comoving,
            m_geometry.cell_width_y_comoving,
            m_geometry.cell_width_z_comoving));

    const double density_frame_factor =
        m_config.units.coordinate_frame == core::CoordinateFrame::kComoving
            ? 1.0 / std::pow(update.scale_factor, 3.0)
            : 1.0;
    const SimulationStateHydroSourcePropertyProvider source_property_provider(
        context.state,
        m_units,
        m_config.numerics.hydro_density_floor_code,
        density_frame_factor,
        m_dense_row_by_geometry_row);
    hydro::HydroSourceContext source_context{
        .update = update,
        .gravity_accel_x_peculiar = {},
        .gravity_accel_y_peculiar = {},
        .gravity_accel_z_peculiar = {},
        .mass_density_physical_cgs = {},
        .hydrogen_number_density_cgs = {},
        .metallicity_mass_fraction = {},
        .temperature_k = {},
        .source_property_provider = &source_property_provider,
        .thermodynamic_closure = m_effective_ism_closure.get(),
        .redshift = std::max(0.0, (update.scale_factor > 0.0 ? (1.0 / update.scale_factor - 1.0) : 0.0)),
    };

    // Gravity and homogeneous expansion are integrated by the gravity-owned
    // KDK half-kicks around hydro/source evolution. Hydro owns transport plus
    // local thermochemical source terms only, preventing double application.
    std::vector<const hydro::HydroSourceTerm*> sources;
    if (m_cooling_source != nullptr) {
      sources.push_back(m_cooling_source.get());
    }
    if (m_effective_ism_energy_source != nullptr) {
      m_effective_ism_energy_source->resetLedger();
      sources.push_back(m_effective_ism_energy_source.get());
    }
    const hydro::HydroActiveBatchPolicy active_batch_policy = selectHydroActiveBatchPolicy();
    auto active_batch_reservation = reserveHydroActiveBatch(active_batch_policy);
    // CHUI Part 2: composition lanes are advected with the final hydro mass
    // flux (physics::CompositionFluxAccumulator). Deltas are staged and
    // committed inside the same StageTransaction as the hydro state.
    std::unique_ptr<physics::CompositionFluxAccumulator> composition_flux;
    core::MemoryReservation composition_transport_reservation;
    if (context.state.gas_composition.enabled()) {
      composition_flux = std::make_unique<physics::CompositionFluxAccumulator>(context.state);
      if (m_memory_governor != nullptr) {
        composition_transport_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            physics::CompositionFluxAccumulator::deltaBufferBytes(
                context.state.gas_composition.laneCount(), context.state.cells.size()),
            "hydro.composition_transport.delta_buffer");
      }
      source_context.face_mass_flux_sink = composition_flux.get();
      composition_flux->beginAdvance(m_geometry, m_dense_row_by_geometry_row, {});
    }
    // Audit F2: a source that cannot integrate the full accepted interval
    // throws StageNumericalRetryRequest from inside the solver. The solver
    // only writes the staging copy m_conserved, so the failure is converted
    // into a transaction rejection below; SimulationState stays untouched.
    std::string source_retry_failure;
    try {
      m_solver.advancePatchActiveSetBatchedWithScratch(
          m_conserved,
          m_geometry,
          active_view,
          update,
          *m_reconstruction,
          *m_riemann_solver,
          sources,
          source_context,
          m_scratch,
          nullptr,
          active_batch_policy,
          &m_last_hydro_profile);
    } catch (const core::StageNumericalRetryRequest& retry) {
      source_retry_failure = retry.what();
    }
    if (composition_flux != nullptr) {
      composition_flux->endAdvance();
    }
    if (context.profiler_session != nullptr) {
      const hydro::HydroConservationReport& conservation = m_last_hydro_profile.conservation;
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.conservation",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.godunov",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "computed volume-integrated hydro conservation totals for the local update stage",
          .payload = {{"cell_count", std::to_string(conservation.cell_count)},
                      {"before_mass", std::to_string(conservation.before.mass)},
                      {"after_mass", std::to_string(conservation.after.mass)},
                      {"flux_delta_mass", std::to_string(conservation.flux_delta.mass)},
                      {"source_delta_mass", std::to_string(conservation.source_delta.mass)},
                      {"floor_delta_mass", std::to_string(conservation.floor_delta.mass)},
                      {"residual_mass", std::to_string(conservation.residual.mass)},
                      {"residual_momentum_x", std::to_string(conservation.residual.momentum_x)},
                      {"residual_momentum_y", std::to_string(conservation.residual.momentum_y)},
                      {"residual_momentum_z", std::to_string(conservation.residual.momentum_z)},
                      {"residual_total_energy", std::to_string(conservation.residual.total_energy)},
                      {"residual_internal_energy", std::to_string(conservation.residual.internal_energy)},
                      {"flux_delta_total_energy", std::to_string(conservation.flux_delta.total_energy)},
                      {"source_delta_total_energy", std::to_string(conservation.source_delta.total_energy)},
                      {"floor_delta_total_energy", std::to_string(conservation.floor_delta.total_energy)},
                      {"internal_energy_floor_count", std::to_string(conservation.internal_energy_floor_count)},
                      {"dual_energy_sync_count", std::to_string(conservation.dual_energy_sync_count)},
                      {"dual_energy_sync_delta_total_energy", std::to_string(conservation.dual_energy_sync_delta.total_energy)},
                      {"riemann_solve_count", std::to_string(m_last_hydro_profile.riemann_solve_count)},
                      {"positivity_fallback_count", std::to_string(m_last_hydro_profile.positivity_fallback_count)},
                      {"stencil_degradation_count", std::to_string(m_last_hydro_profile.stencil_degradation_count)},
                      {"ctu_transverse_corrected_states", std::to_string(m_last_hydro_profile.ctu_transverse_corrected_states)},
                      {"ctu_transverse_degraded_states", std::to_string(m_last_hydro_profile.ctu_transverse_degraded_states)},
                      {"ctu_transverse_ghost_corrected_states", std::to_string(m_last_hydro_profile.ctu_transverse_ghost_corrected_states)},
                      {"ctu_transverse_ghost_axes", std::to_string(m_last_hydro_profile.ctu_transverse_ghost_axes)},
                      {"ctu_transverse_inadmissible_states", std::to_string(m_last_hydro_profile.ctu_transverse_inadmissible_states)},
                      {"mass_tolerance", "1e-10"},
                      {"momentum_tolerance", "1e-10"},
                      {"energy_tolerance", "1e-10"}}});
    }

    const HydroConservativeGhostSyncReport ghost_sync_report = restoreHydroGhostConservedCells(
        m_conserved, ghost_conserved_snapshot, world_rank);
    std::vector<parallel::HydroConservativeFluxCorrectionRecord> global_flux_corrections;
    if (all_ranks_have_hydro_cells) {
      global_flux_corrections =
          parallel::executeBlockingHydroConservativeFluxCorrectionExchange(
              m_mpi_context, ghost_sync_report.correction_records, context.integrator_state.step_index);
    } else if (!ghost_sync_report.correction_records.empty()) {
      throw std::runtime_error(
          "hydro conservative flux correction produced MPI correction records while at least one rank has no hydro cells");
    }
    std::size_t applied_flux_corrections = 0;
    double applied_flux_delta_l1 = 0.0;
    for (const parallel::HydroConservativeFluxCorrectionRecord& correction : global_flux_corrections) {
      if (correction.owner_rank != static_cast<int>(world_rank)) {
        continue;
      }
      const auto correction_row = context.state.rowForGasCellId(correction.gas_cell_id);
      if (!correction_row.has_value()) {
        throw std::runtime_error("hydro conservative flux correction targeted an unknown local gas_cell_id");
      }
      const std::uint32_t cell_index = *correction_row;
      const std::uint32_t owner_rank = internal::gasCellOwnerRankForLocalRow(
          context.state,
          cell_index,
          particle_row_by_id,
          "hydro conservative flux correction");
      if (owner_rank != world_rank) {
        throw std::runtime_error("hydro conservative flux correction targeted a non-authoritative local ghost cell");
      }
      if (cell_index >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro flux correction dense row is outside Cartesian geometry map");
      }
      const std::uint32_t geometry_row = m_geometry_row_by_dense_row[cell_index];
      hydro::HydroConservedState owned = m_conserved.loadCell(geometry_row);
      owned.mass_density_comoving += correction.delta_mass_density_comoving;
      owned.momentum_density_x_comoving += correction.delta_momentum_density_x_comoving;
      owned.momentum_density_y_comoving += correction.delta_momentum_density_y_comoving;
      owned.momentum_density_z_comoving += correction.delta_momentum_density_z_comoving;
      owned.total_energy_density_comoving += correction.delta_total_energy_density_comoving;
      owned.metal_mass_density_comoving += correction.delta_metal_mass_density_comoving;
      m_conserved.storeCell(geometry_row, owned);
      ++applied_flux_corrections;
      applied_flux_delta_l1 += std::abs(correction.delta_mass_density_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_momentum_density_x_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_momentum_density_y_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_momentum_density_z_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_total_energy_density_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_metal_mass_density_comoving);
    }
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.conservative_ghost_sync",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.godunov",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "restored imported hydro ghosts and applied conservative owner-side boundary corrections",
          .payload = {{"restored_ghost_cells", std::to_string(ghost_sync_report.restored_ghost_cells)},
                      {"local_correction_records", std::to_string(ghost_sync_report.correction_records.size())},
                      {"global_correction_records", std::to_string(global_flux_corrections.size())},
                      {"applied_flux_corrections", std::to_string(applied_flux_corrections)},
                      {"rejected_remote_delta_l1", std::to_string(ghost_sync_report.rejected_remote_delta_l1)},
                      {"applied_flux_delta_l1", std::to_string(applied_flux_delta_l1)},
                      {"imported_mpi_ghosts", std::to_string(m_last_remote_ghost_boundary_report.imported_ghosts)},
                      {"requested_remote_cells", std::to_string(m_last_remote_ghost_boundary_report.requested_cells)},
                      {"received_remote_cells", std::to_string(m_last_remote_ghost_boundary_report.received_cells)},
                      {"remote_interface_faces", std::to_string(m_last_remote_ghost_boundary_report.interface_faces)},
                      {"remote_stale_or_invalid_payloads", std::to_string(m_last_remote_ghost_boundary_report.stale_or_invalid_payloads)},
                      {"remote_request_bytes", std::to_string(m_last_remote_ghost_boundary_report.request_bytes)},
                      {"remote_payload_bytes", std::to_string(m_last_remote_ghost_boundary_report.payload_bytes)}}});
    }

    // Part 1 stage transaction. m_conserved is solver-owned staging storage;
    // SimulationState gas truth is written only inside commit() below, after
    // the staged state passed admissibility and the (communicator-wide) vote.
    // A rejection therefore leaves density, momentum, energy, metals,
    // identity maps and time metadata exactly as they were at step start.
    const internal::MpiStageCollectiveDecision hydro_collective_decision(m_mpi_context);
    core::StageTransaction hydro_transaction(
        "hydro.fixed_grid.update", m_memory_governor, &hydro_collective_decision);
    if (active_batch_reservation.has_value()) {
      hydro_transaction.adoptReservation(std::move(*active_batch_reservation));
      active_batch_reservation.reset();
    }
    if (composition_transport_reservation.pending()) {
      hydro_transaction.adoptReservation(std::move(composition_transport_reservation));
    }
    std::string staged_failure;
    bool staged_admissible = source_retry_failure.empty() &&
        stagedConservedStateIsAdmissible(staged_failure);
    if (!source_retry_failure.empty()) {
      staged_failure = "hydro source interval incomplete: " + source_retry_failure;
    }
    if (staged_admissible && composition_flux != nullptr) {
      // Composition deltas use the step-start sidecar and need the
      // post-hydro mass only for reconciliation inside commit().
      try {
        composition_flux->prepare(context.state);
        admitCompositionSparseActivation(
            context.state, *composition_flux, m_memory_governor,
            "hydro.fixed_grid.adaptive_chemistry_rows");
      } catch (const core::MemoryAdmissionError& error) {
        staged_admissible = false;
        staged_failure = std::string("composition sparse capacity: ") + error.what();
        hydro_transaction.reject(core::StageRejectionKind::kResourceRefused, staged_failure);
      } catch (const std::exception& error) {
        staged_admissible = false;
        staged_failure = std::string("composition transport: ") + error.what();
      }
    }
    if (!hydro_transaction.validate(
            staged_admissible, core::StageRejectionKind::kNumericalRetry, staged_failure)) {
      throw core::StageTransactionRejected(
          hydro_transaction.outcome().rejection,
          "fixed-grid hydro update rejected before commit (" +
              std::string(core::stageRejectionKindLabel(hydro_transaction.outcome().rejection)) +
              "): " + hydro_transaction.outcome().reason);
    }
    (void)hydro_transaction.commit([&]() {
    for (std::size_t cell_index = 0; cell_index < context.state.cells.size(); ++cell_index) {
      const std::uint32_t gas_cell_owner_rank = internal::gasCellOwnerRankForLocalRow(
          context.state,
          static_cast<std::uint32_t>(cell_index),
          particle_row_by_id,
          "hydro callback primitive store");
      if (gas_cell_owner_rank != world_rank) {
        continue;
      }
      if (cell_index >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro primitive store dense row is outside Cartesian geometry map");
      }
      const std::uint32_t geometry_row = m_geometry_row_by_dense_row[cell_index];
      const hydro::HydroConservedState stored_conserved = m_conserved.loadCell(geometry_row);
      hydro::HydroPrimitiveState primitive =
          hydro::HydroCoreSolver::primitiveFromConserved(stored_conserved, k_gamma_adiabatic);
      if (m_conserved.auxiliaryInternalEnergyEnabled()) {
        // Dual energy: the solver wrote the SELECTED thermal energy (E - K or
        // the auxiliary value) into the lane; the non-conservative difference
        // is ledgered in HydroConservationReport::dual_energy_sync_delta.
        const double rho_e = m_conserved.auxiliaryInternalEnergyDensityComoving()[geometry_row];
        if (std::isfinite(rho_e) && rho_e > 0.0) {
          primitive.specific_internal_energy_code = rho_e / primitive.rho_comoving;
          primitive.pressure_comoving = (k_gamma_adiabatic - 1.0) * rho_e;
          primitive.signal_speed_squared_code =
              k_gamma_adiabatic * primitive.pressure_comoving / primitive.rho_comoving;
        }
      }
      if (m_effective_ism_closure != nullptr) {
        const auto closure = m_effective_ism_closure->evaluate(
            cell_index,
            stored_conserved,
            primitive,
            update.scale_factor,
            source_context.redshift);
        if (closure.valid && closure.pressure_comoving > 0.0) {
          primitive.pressure_comoving = closure.pressure_comoving;
          primitive.signal_speed_squared_code = closure.signal_speed_squared_code;
          primitive.uses_effective_ism = closure.uses_effective_ism;
        }
      }
      context.state.gas_cells.density_code[cell_index] = primitive.rho_comoving;
      context.state.gas_cells.pressure_code[cell_index] = primitive.pressure_comoving;
      context.state.gas_cells.internal_energy_code[cell_index] =
          primitive.specific_internal_energy_code;
      context.state.gas_cells.sound_speed_code[cell_index] =
          std::sqrt(std::max(primitive.signal_speed_squared_code, 0.0));
      double temperature_factor = m_temperature_k_per_internal_energy_code;
      if (context.state.gas_composition.enabled()) {
        // CHUI Part 2: composition-dependent mu (amu) from thermochemistry.
        const double mu_amu = context.state.gas_composition.mean_molecular_weight[cell_index];
        if (std::isfinite(mu_amu) && mu_amu > 0.0) {
          temperature_factor *= mu_amu * (1.66053906660e-24 / 1.67262192369e-24) /
              hydro::k_part1_reference_mean_molecular_weight;
        }
      }
      context.state.gas_cells.temperature_code[cell_index] =
          temperature_factor * primitive.specific_internal_energy_code;
      context.state.cells.mass_code[cell_index] = primitive.rho_comoving * m_geometry.cell_volume_comoving;
      context.state.gas_cells.metal_mass_code[cell_index] =
          std::clamp(primitive.metallicity_mass_fraction, 0.0, 1.0) *
          context.state.cells.mass_code[cell_index];
      context.state.gas_cells.velocity_x_peculiar[cell_index] = primitive.vel_x_peculiar;
      context.state.gas_cells.velocity_y_peculiar[cell_index] = primitive.vel_y_peculiar;
      context.state.gas_cells.velocity_z_peculiar[cell_index] = primitive.vel_z_peculiar;
    }
    if (composition_flux != nullptr) {
      m_last_composition_transport = composition_flux->commit(context.state);
      // Z is a derived mirror of element (+dust) inventory when composition
      // is authoritative: the hydro passive metal lane is overwritten.
      (void)physics::synchronizeMetalMirrorFromComposition(context.state, std::span<const std::uint32_t>{});
    }
    });

    internal::synchronizeParentParticleCompatibilityMirrors(
        context.state, world_rank, "hydro callback parent compatibility mirror");
    persistEffectiveIsmDiagnostics(context.state);
    context.state.bumpGravitySourceGeneration();
    recordHydroFloorDiagnostics(context);
  }

 private:
  [[nodiscard]] amr::AmrThermodynamicCachePolicy amrThermodynamicCachePolicy(
      double scale_factor, double redshift) const {
    return amr::AmrThermodynamicCachePolicy{
        .temperature_k_per_internal_energy_code = m_temperature_k_per_internal_energy_code,
        .closure = m_effective_ism_closure.get(),
        .scale_factor = scale_factor,
        .redshift = redshift,
        .dual_energy = m_solver.numericsPolicy().dual_energy.mode != hydro::HydroDualEnergyMode::kOff,
        .dual_energy_eta = m_solver.numericsPolicy().dual_energy.thermal_fraction_switch};
  }

  // Validation predicate of the fixed-grid hydro transaction: every staged
  // real-cell conserved lane is finite, mass is positive and the thermal
  // energy E - |p|^2/(2 rho) is non-negative. O(N_cells) streaming check over
  // the solver-owned staging copy; no allocation.
  [[nodiscard]] bool stagedConservedStateIsAdmissible(std::string& failure) const {
    const std::size_t real_cells = m_geometry.cellCount();
    for (std::size_t row = 0; row < real_cells; ++row) {
      const hydro::HydroConservedState cell = m_conserved.loadCell(row);
      const bool finite = std::isfinite(cell.mass_density_comoving) &&
          std::isfinite(cell.momentum_density_x_comoving) &&
          std::isfinite(cell.momentum_density_y_comoving) &&
          std::isfinite(cell.momentum_density_z_comoving) &&
          std::isfinite(cell.total_energy_density_comoving) &&
          std::isfinite(cell.metal_mass_density_comoving);
      if (!finite || !(cell.mass_density_comoving > 0.0)) {
        failure = "staged conserved state is non-finite or has non-positive mass at solver row " +
            std::to_string(row);
        return false;
      }
      const double kinetic = 0.5 *
          (cell.momentum_density_x_comoving * cell.momentum_density_x_comoving +
           cell.momentum_density_y_comoving * cell.momentum_density_y_comoving +
           cell.momentum_density_z_comoving * cell.momentum_density_z_comoving) /
          cell.mass_density_comoving;
      if (cell.total_energy_density_comoving < kinetic) {
        failure = "staged conserved state has negative thermal energy at solver row " +
            std::to_string(row);
        return false;
      }
    }
    return true;
  }

  // Per-active-cell transient budget. The optional CTU / dual-energy lanes
  // scale with touched cells (<= 1 + 2 * 6 per active cell).
  [[nodiscard]] std::uint64_t hydroBatchScratchBytesPerCell() const {
    constexpr std::uint64_t k_touched_per_active_cell = 13U;
    const hydro::HydroActiveBatchPolicy features =
        m_solver.requiredWorkspaceFeatures(hydro::HydroActiveBatchPolicy{});
    std::uint64_t bytes = hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell;
    if (features.ctu_transverse_workspace) {
      bytes += k_touched_per_active_cell * sizeof(std::array<hydro::HydroConservedState, 3>);
    }
    if (features.dual_energy_workspace) {
      bytes += k_touched_per_active_cell * 2U * sizeof(double);
    }
    return bytes;
  }

  [[nodiscard]] hydro::HydroActiveBatchPolicy selectHydroActiveBatchPolicy() const {
    std::uint64_t requested = m_config.numerics.hydro_active_batch_max_cells == 0U
        ? static_cast<std::uint64_t>(hydro::k_hydro_automatic_active_batch_max_cells)
        : m_config.numerics.hydro_active_batch_max_cells;
    if (m_memory_governor != nullptr) {
      const core::DeterministicBatchSizingResult sizing =
          core::selectDeterministicBatchSize(
              m_memory_governor->snapshot(),
              core::DeterministicBatchSizingPolicy{
                  .requested_max_items = requested,
                  .bytes_per_item = hydroBatchScratchBytesPerCell(),
                  .fixed_reserve_bytes = 0U,
                  .minimum_items = 1U,
                  .alignment_items = hydro::k_hydro_active_batch_alignment_cells,
                  .headroom_use_basis_points = 10000U,
              });
      if (sizing.selected_items == 0U) {
        throw std::runtime_error(
            "hydro active-batch headroom cannot admit one active cell");
      }
      requested = sizing.selected_items;
    }
    requested = std::min<std::uint64_t>(
        requested, static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()));
    return m_solver.requiredWorkspaceFeatures(hydro::HydroActiveBatchPolicy{
        .max_active_cells = static_cast<std::size_t>(std::max<std::uint64_t>(1U, requested))});
  }

  [[nodiscard]] std::optional<core::MemoryReservation> reserveHydroActiveBatch(
      const hydro::HydroActiveBatchPolicy& batch_policy) const {
    if (m_memory_governor == nullptr) {
      return std::nullopt;
    }
    const std::uint64_t cells = static_cast<std::uint64_t>(batch_policy.max_active_cells);
    const std::uint64_t bytes_per_cell = hydroBatchScratchBytesPerCell();
    if (cells > std::numeric_limits<std::uint64_t>::max() / bytes_per_cell) {
      throw std::overflow_error("hydro active-batch reservation byte count overflow");
    }
    auto reservation = m_memory_governor->reserve(
        core::MemoryClass::kScratchArena,
        cells * bytes_per_cell,
        "hydro.active_batch");
    reservation.commit();
    return std::optional<core::MemoryReservation>(std::move(reservation));
  }

  [[nodiscard]] HydroFloorApplicationDiagnostics auditConfiguredHydroFloors(
      const core::SimulationState& state) const {
    HydroFloorApplicationDiagnostics diagnostics;
    for (std::size_t row = 0; row < state.gas_cells.size(); ++row) {
      const double density = state.gas_cells.density_code[row];
      const double pressure = state.gas_cells.pressure_code[row];
      if (density < m_config.numerics.hydro_density_floor_code) {
        ++diagnostics.density_floor_cell_count;
        diagnostics.density_floor_correction_l1_code +=
            m_config.numerics.hydro_density_floor_code - density;
      }
      if (pressure < m_config.numerics.hydro_pressure_floor_code) {
        ++diagnostics.pressure_floor_cell_count;
        diagnostics.pressure_floor_correction_l1_code +=
            m_config.numerics.hydro_pressure_floor_code - pressure;
      }
    }
    return diagnostics;
  }

  void recordHydroFloorDiagnostics(const core::StepContext& context) const {
    if (context.profiler_session == nullptr) return;
    context.profiler_session->recordEvent(core::RuntimeEvent{
        .event_kind = "hydro.floor_policy",
        .severity = (m_last_floor_diagnostics.density_floor_cell_count != 0U ||
                     m_last_floor_diagnostics.pressure_floor_cell_count != 0U)
            ? core::RuntimeEventSeverity::kWarning
            : core::RuntimeEventSeverity::kInfo,
        .subsystem = "hydro.floors",
        .step_index = context.integrator_state.step_index,
        .simulation_time_code = context.integrator_state.current_time_code,
        .scale_factor = context.integrator_state.current_scale_factor,
        .message = "configured hydro floor policy audit before the hydro update",
        .payload = {
            {"density_floor_code", std::to_string(m_config.numerics.hydro_density_floor_code)},
            {"pressure_floor_code", std::to_string(m_config.numerics.hydro_pressure_floor_code)},
            {"density_floor_cell_count", std::to_string(m_last_floor_diagnostics.density_floor_cell_count)},
            {"pressure_floor_cell_count", std::to_string(m_last_floor_diagnostics.pressure_floor_cell_count)},
            {"density_floor_correction_l1_code", std::to_string(m_last_floor_diagnostics.density_floor_correction_l1_code)},
            {"pressure_floor_correction_l1_code", std::to_string(m_last_floor_diagnostics.pressure_floor_correction_l1_code)},
        }});
  }

  void persistEffectiveIsmDiagnostics(core::SimulationState& state) const {
    if (m_effective_ism_energy_source == nullptr || m_effective_ism_closure == nullptr) return;
    const physics::EffectiveIsmEnergyLedger latest = m_effective_ism_energy_source->ledger();
    double cumulative_added = 0.0;
    double cumulative_removed = 0.0;
    double cumulative_net = 0.0;
    std::uint64_t cumulative_cells = 0U;
    bool legacy_density_ledger_discarded = false;
    if (const core::ModuleSidecarBlock* prior = state.sidecars.find("effective_multiphase_ism");
        prior != nullptr) {
      std::string text;
      text.reserve(prior->payload.size());
      for (const std::byte value : prior->payload) text.push_back(static_cast<char>(value));
      std::istringstream lines(text);
      std::string line;
      while (std::getline(lines, line)) {
        const auto split = line.find('=');
        if (split == std::string::npos) continue;
        const std::string key = line.substr(0, split);
        const std::string value = line.substr(split + 1U);
        try {
          if (key == "cumulative.energy_added_code") cumulative_added = std::stod(value);
          else if (key == "cumulative.energy_removed_code") cumulative_removed = std::stod(value);
          else if (key == "cumulative.net_energy_adjustment_code") cumulative_net = std::stod(value);
          else if (key == "cumulative.adjusted_cell_count") cumulative_cells = std::stoull(value);
        } catch (const std::exception&) {
          cumulative_added = cumulative_removed = cumulative_net = 0.0;
          cumulative_cells = 0U;
          break;
        }
      }
    }
    if (const core::ModuleSidecarBlock* prior = state.sidecars.find("effective_multiphase_ism");
        prior != nullptr) {
      std::string text;
      text.reserve(prior->payload.size());
      for (const std::byte value : prior->payload) text.push_back(static_cast<char>(value));
      const std::string version_key = "energy_ledger_version=" +
          std::to_string(physics::kEffectiveIsmEnergyLedgerVersion);
      const bool carried_legacy_flag =
          text.find("legacy_density_ledger_discarded=1") != std::string::npos;
      if (text.find(version_key) == std::string::npos) {
        // Pre-F7 sidecar: cumulative values are energy-density sums without
        // cell volume. They cannot be converted after the fact, so discard
        // them explicitly and record that the cumulative ledger restarted.
        cumulative_added = 0.0;
        cumulative_removed = 0.0;
        cumulative_net = 0.0;
        cumulative_cells = 0U;
        legacy_density_ledger_discarded = true;
      } else if (carried_legacy_flag) {
        legacy_density_ledger_discarded = true;
      }
    }
    cumulative_added += latest.energy_added_code;
    cumulative_removed += latest.energy_removed_code;
    cumulative_net += latest.net_energy_adjustment_code;
    cumulative_cells += latest.adjusted_cell_count;
    std::ostringstream stream;
    stream << "module_name=effective_multiphase_ism\n";
    stream << "schema_version=" << physics::kEffectiveMultiphaseEosSchemaVersion << "\n";
    stream << "parameter_set=" << m_effective_ism_closure->table().config().parameter_set << "\n";
    stream << "table_hash=" << m_effective_ism_closure->table().tableHashHex() << "\n";
    stream << "cooling_reference=" << m_effective_ism_closure->table().coolingReferenceDescription() << "\n";
    stream << "energy_ledger_version=" << physics::kEffectiveIsmEnergyLedgerVersion << "\n";
    stream << "energy_ledger_units=integrated_code_energy\n";
    stream << "legacy_density_ledger_discarded=" << (legacy_density_ledger_discarded ? 1 : 0) << "\n";
    stream << "latest.energy_added_code=" << latest.energy_added_code << "\n";
    stream << "latest.energy_removed_code=" << latest.energy_removed_code << "\n";
    stream << "latest.net_energy_adjustment_code=" << latest.net_energy_adjustment_code << "\n";
    stream << "latest.adjusted_cell_count=" << latest.adjusted_cell_count << "\n";
    stream << "cumulative.energy_added_code=" << cumulative_added << "\n";
    stream << "cumulative.energy_removed_code=" << cumulative_removed << "\n";
    stream << "cumulative.net_energy_adjustment_code=" << cumulative_net << "\n";
    stream << "cumulative.adjusted_cell_count=" << cumulative_cells << "\n";
    const std::string text = stream.str();
    core::ModuleSidecarBlock block;
    block.module_name = "effective_multiphase_ism";
    block.schema_version = physics::kEffectiveMultiphaseEosSchemaVersion;
    block.payload.resize(text.size());
    for (std::size_t i = 0; i < text.size(); ++i) block.payload[i] = static_cast<std::byte>(text[i]);
    state.sidecars.upsert(std::move(block));
  }

  [[nodiscard]] static hydro::HydroBoundaryKind hydroBoundaryKindFromModePolicy(
      core::BoundaryCondition boundary_condition) {
    switch (boundary_condition) {
      case core::BoundaryCondition::kPeriodic:
        return hydro::HydroBoundaryKind::kPeriodic;
      case core::BoundaryCondition::kOpen:
        return hydro::HydroBoundaryKind::kOpen;
      case core::BoundaryCondition::kReflective:
        return hydro::HydroBoundaryKind::kReflective;
    }
    return hydro::HydroBoundaryKind::kOpen;
  }

  [[nodiscard]] std::array<double, 3> ghostCenterForBoundaryFace(
      const core::SimulationState& state,
      const hydro::HydroGhostCell& ghost) const {
    const std::uint32_t dense_row = m_dense_row_by_geometry_row.at(ghost.owner_real_cell);
    std::array<double, 3> center{
        state.cells.center_x_comoving[dense_row],
        state.cells.center_y_comoving[dense_row],
        state.cells.center_z_comoving[dense_row]};
    const double sign = ghost.side == hydro::HydroFaceSide::kLower ? -1.0 : 1.0;
    switch (ghost.axis) {
      case hydro::HydroFaceAxis::kX:
        center[0] += sign * m_geometry.cell_width_x_comoving;
        break;
      case hydro::HydroFaceAxis::kY:
        center[1] += sign * m_geometry.cell_width_y_comoving;
        break;
      case hydro::HydroFaceAxis::kZ:
        center[2] += sign * m_geometry.cell_width_z_comoving;
        break;
    }
    return center;
  }

  [[nodiscard]] bool isBoundaryGeometryRow(std::size_t geometry_row) const {
    const auto ijk = m_geometry.cellIjk(geometry_row);
    return ijk[0] == 0U || ijk[0] + 1U == m_geometry.nx ||
        ijk[1] == 0U || ijk[1] + 1U == m_geometry.ny ||
        ijk[2] == 0U || ijk[2] + 1U == m_geometry.nz;
  }

  [[nodiscard]] std::vector<HydroBoundaryCellAdvertisement> buildLocalHydroBoundaryAdvertisements(
      const core::SimulationState& state,
      const std::unordered_map<std::uint64_t, std::uint32_t>& particle_row_by_id,
      std::uint64_t hydro_sync_epoch,
      std::uint64_t decomposition_epoch) const {
    std::vector<HydroBoundaryCellAdvertisement> records;
    records.reserve(state.cells.size());
    const std::uint32_t world_rank = static_cast<std::uint32_t>(m_mpi_context.worldRank());
    for (std::uint32_t dense_row = 0; dense_row < state.cells.size(); ++dense_row) {
      if (internal::gasCellOwnerRankForLocalRow(state, dense_row, particle_row_by_id, "hydro boundary advertisement") !=
          world_rank) {
        continue;
      }
      if (dense_row >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro boundary advertisement dense row is outside geometry map");
      }
      const std::size_t geometry_row = m_geometry_row_by_dense_row[dense_row];
      if (!isBoundaryGeometryRow(geometry_row)) {
        continue;
      }
      const core::GasCellIdentityRecord& identity = internal::gasCellIdentityRecordForLocalRow(
          state,
          dense_row,
          "hydro boundary advertisement");
      records.push_back(HydroBoundaryCellAdvertisement{
          .gas_cell_id = identity.gas_cell_id,
          .owner_rank = static_cast<int>(world_rank),
          .hydro_sync_epoch = hydro_sync_epoch,
          .decomposition_epoch = decomposition_epoch,
          .center_x_comoving = state.cells.center_x_comoving[dense_row],
          .center_y_comoving = state.cells.center_y_comoving[dense_row],
          .center_z_comoving = state.cells.center_z_comoving[dense_row]});
    }
    return records;
  }

  [[nodiscard]] std::optional<HydroBoundaryCellAdvertisement> findRemoteBoundaryAdvertisement(
      std::span<const HydroBoundaryCellAdvertisement> advertisements,
      const std::array<double, 3>& owner_center,
      hydro::HydroFaceAxis axis,
      hydro::HydroFaceSide side,
      std::uint64_t hydro_sync_epoch) const {
    const double tolerance = 1.0e-8 * std::max({
        1.0,
        std::abs(m_geometry.cell_width_x_comoving),
        std::abs(m_geometry.cell_width_y_comoving),
        std::abs(m_geometry.cell_width_z_comoving)});
    std::optional<HydroBoundaryCellAdvertisement> match;
    double best_normal_distance = std::numeric_limits<double>::infinity();
    for (const HydroBoundaryCellAdvertisement& record : advertisements) {
      if (record.owner_rank == m_mpi_context.worldRank()) {
        continue;
      }
      if (record.hydro_sync_epoch != hydro_sync_epoch) {
        continue;
      }
      double normal_distance = 0.0;
      bool transverse_match = false;
      switch (axis) {
        case hydro::HydroFaceAxis::kX:
          normal_distance = record.center_x_comoving - owner_center[0];
          transverse_match = nearlySameCoordinate(record.center_y_comoving, owner_center[1], tolerance) &&
              nearlySameCoordinate(record.center_z_comoving, owner_center[2], tolerance);
          break;
        case hydro::HydroFaceAxis::kY:
          normal_distance = record.center_y_comoving - owner_center[1];
          transverse_match = nearlySameCoordinate(record.center_x_comoving, owner_center[0], tolerance) &&
              nearlySameCoordinate(record.center_z_comoving, owner_center[2], tolerance);
          break;
        case hydro::HydroFaceAxis::kZ:
          normal_distance = record.center_z_comoving - owner_center[2];
          transverse_match = nearlySameCoordinate(record.center_x_comoving, owner_center[0], tolerance) &&
              nearlySameCoordinate(record.center_y_comoving, owner_center[1], tolerance);
          break;
      }
      const bool correct_side = side == hydro::HydroFaceSide::kUpper
          ? normal_distance > tolerance
          : normal_distance < -tolerance;
      if (!transverse_match || !correct_side) {
        continue;
      }
      const double abs_distance = std::abs(normal_distance);
      if (abs_distance + tolerance < best_normal_distance) {
        best_normal_distance = abs_distance;
        match = record;
      } else if (std::abs(abs_distance - best_normal_distance) <= tolerance) {
        throw std::runtime_error("hydro remote boundary discovery found duplicate remote cells for one face");
      }
    }
    return match;
  }

  [[nodiscard]] HydroRemoteGhostBoundaryReport refreshRemoteHydroBoundaryGhosts(
      core::StepContext& context,
      const std::unordered_map<std::uint64_t, std::uint32_t>& particle_row_by_id) {
    HydroRemoteGhostBoundaryReport report;
    if (!m_mpi_context.isEnabled() || m_mpi_context.worldSize() <= 1) {
      return report;
    }
    const std::uint64_t hydro_sync_epoch =
        context.integrator_state.step_index * core::integrationStageCount() +
        core::integrationStageIndex(context.stage) + 1U;
    const std::uint64_t decomposition_epoch = context.state.gasCellIdentityGeneration();
    const auto local_advertisements = buildLocalHydroBoundaryAdvertisements(
        context.state,
        particle_row_by_id,
        hydro_sync_epoch,
        decomposition_epoch);
    const auto advertisements = exchangeHydroBoundaryAdvertisements(m_mpi_context, local_advertisements);
    report.boundary_metadata_sent = local_advertisements.size();
    report.boundary_metadata_received = advertisements.size();

    std::vector<parallel::HydroGhostCellRequest> requests;
    std::unordered_map<std::uint64_t, std::size_t> ghost_slot_by_face_key;
    m_remote_non_authority_faces.clear();
    for (std::size_t face_index = 0; face_index < m_geometry.faces.size(); ++face_index) {
      hydro::HydroFace& face = m_geometry.faces[face_index];
      if (face.ghost_cell_slot != hydro::k_invalid_ghost_cell_slot ||
          face.neighbor_cell == hydro::k_invalid_cell_index ||
          face.owner_cell >= m_dense_row_by_geometry_row.size() ||
          face.neighbor_cell >= m_dense_row_by_geometry_row.size()) {
        continue;
      }
      const std::uint32_t owner_dense_row = m_dense_row_by_geometry_row[face.owner_cell];
      const std::uint32_t neighbor_dense_row = m_dense_row_by_geometry_row[face.neighbor_cell];
      const std::uint32_t owner_rank = internal::gasCellOwnerRankForLocalRow(
          context.state,
          owner_dense_row,
          particle_row_by_id,
          "hydro internal remote boundary discovery");
      const std::uint32_t neighbor_rank = internal::gasCellOwnerRankForLocalRow(
          context.state,
          neighbor_dense_row,
          particle_row_by_id,
          "hydro internal remote boundary discovery");
      if (owner_rank == neighbor_rank) {
        continue;
      }
      if (owner_rank != static_cast<std::uint32_t>(m_mpi_context.worldRank())) {
        m_remote_non_authority_faces.insert(face_index);
        continue;
      }
      const core::GasCellIdentityRecord& owner_identity = internal::gasCellIdentityRecordForLocalRow(
          context.state,
          owner_dense_row,
          "hydro internal remote boundary discovery");
      const core::GasCellIdentityRecord& remote_identity = internal::gasCellIdentityRecordForLocalRow(
          context.state,
          neighbor_dense_row,
          "hydro internal remote boundary discovery");
      const std::uint64_t face_key =
          hydroInterfaceFaceKey(owner_identity.gas_cell_id, remote_identity.gas_cell_id, face.axis);
      const std::size_t ghost_slot = m_geometry.ghost_cells.size();
      const std::size_t ghost_cell = m_geometry.cellCount() + ghost_slot;
      m_geometry.ghost_cells.push_back(hydro::HydroGhostCell{
          .owner_real_cell = face.owner_cell,
          .source_real_cell = face.owner_cell,
          .ghost_cell = ghost_cell,
          .ghost_slot = ghost_slot,
          .boundary_kind = hydro::HydroBoundaryKind::kImportedMpi,
          .axis = face.axis,
          .side = hydro::HydroFaceSide::kUpper,
          .mutation_rights = hydro::HydroGhostMutationRights::kReadOnlyImported,
          .origin_gas_cell_id = remote_identity.gas_cell_id,
          .origin_rank = static_cast<int>(neighbor_rank),
          .hydro_sync_epoch = hydro_sync_epoch,
          .decomposition_epoch = decomposition_epoch});
      face.neighbor_cell = ghost_cell;
      face.neighbor_plus_cell = hydro::k_invalid_cell_index;
      face.neighbor_plus2_cell = hydro::k_invalid_cell_index;
      face.ghost_cell_slot = ghost_slot;
      if (!ghost_slot_by_face_key.emplace(face_key, ghost_slot).second) {
        throw std::runtime_error("hydro internal remote boundary discovery produced duplicate face keys");
      }
      requests.push_back(parallel::HydroGhostCellRequest{
          .descriptor = parallel::HydroGhostCellDescriptor{
              .gas_cell_id = remote_identity.gas_cell_id,
              .owner_rank = static_cast<int>(neighbor_rank),
              .consumer_rank = m_mpi_context.worldRank(),
              .hydro_sync_epoch = hydro_sync_epoch,
              .decomposition_epoch = decomposition_epoch,
              .boundary_state_only = true},
          .face_key = face_key,
          .axis = hydroAxisCode(face.axis),
          .side = 1U});
    }
    if (m_conserved.size() < m_geometry.totalCellStorageCount()) {
      m_conserved.resize(m_geometry.totalCellStorageCount());
    }
    for (hydro::HydroGhostCell& ghost : m_geometry.ghost_cells) {
      if (ghost.boundary_kind == hydro::HydroBoundaryKind::kImportedMpi &&
          ghost.hydro_sync_epoch != hydro_sync_epoch) {
        ghost.boundary_kind = hydroBoundaryKindFromModePolicy(m_mode_policy.hydro_boundary);
        ghost.mutation_rights = hydro::HydroGhostMutationRights::kWritablePhysicalBoundaryScratch;
        ghost.origin_gas_cell_id = 0U;
        ghost.origin_rank = -1;
        ghost.hydro_sync_epoch = 0U;
        ghost.decomposition_epoch = 0U;
      }
      const std::uint32_t owner_dense_row = m_dense_row_by_geometry_row.at(ghost.owner_real_cell);
      if (internal::gasCellOwnerRankForLocalRow(
              context.state,
              owner_dense_row,
              particle_row_by_id,
              "hydro remote boundary discovery") != static_cast<std::uint32_t>(m_mpi_context.worldRank())) {
        continue;
      }
      const std::array<double, 3> owner_center{
          context.state.cells.center_x_comoving[owner_dense_row],
          context.state.cells.center_y_comoving[owner_dense_row],
          context.state.cells.center_z_comoving[owner_dense_row]};
      const auto remote = findRemoteBoundaryAdvertisement(
          advertisements,
          owner_center,
          ghost.axis,
          ghost.side,
          hydro_sync_epoch);
      if (!remote.has_value()) {
        continue;
      }
      const core::GasCellIdentityRecord& local_identity = internal::gasCellIdentityRecordForLocalRow(
          context.state,
          owner_dense_row,
          "hydro remote boundary discovery");
      const std::uint64_t face_key = hydroInterfaceFaceKey(local_identity.gas_cell_id, remote->gas_cell_id, ghost.axis);
      if (!ghost_slot_by_face_key.emplace(face_key, ghost.ghost_slot).second) {
        throw std::runtime_error("hydro remote boundary discovery produced duplicate face keys");
      }
      ghost.boundary_kind = hydro::HydroBoundaryKind::kImportedMpi;
      ghost.mutation_rights = hydro::HydroGhostMutationRights::kReadOnlyImported;
      ghost.source_real_cell = ghost.owner_real_cell;
      ghost.origin_gas_cell_id = remote->gas_cell_id;
      ghost.origin_rank = remote->owner_rank;
      ghost.hydro_sync_epoch = hydro_sync_epoch;
      ghost.decomposition_epoch = decomposition_epoch;
      requests.push_back(parallel::HydroGhostCellRequest{
          .descriptor = parallel::HydroGhostCellDescriptor{
              .gas_cell_id = remote->gas_cell_id,
              .owner_rank = remote->owner_rank,
              .consumer_rank = m_mpi_context.worldRank(),
              .hydro_sync_epoch = hydro_sync_epoch,
              .decomposition_epoch = decomposition_epoch,
              .boundary_state_only = true},
          .face_key = face_key,
          .axis = hydroAxisCode(ghost.axis),
          .side = hydroSideCode(ghost.side)});
    }

    const auto global_requests = parallel::executeBlockingHydroGhostCellRequestExchange(
        m_mpi_context,
        requests,
        hydro_sync_epoch);
    std::vector<parallel::HydroGhostCellPayloadRecord> local_payloads;
    for (const parallel::HydroGhostCellRequest& request : global_requests) {
      if (request.descriptor.owner_rank != m_mpi_context.worldRank()) {
        continue;
      }
      const auto row = context.state.rowForGasCellId(request.descriptor.gas_cell_id);
      if (!row.has_value()) {
        throw std::runtime_error("hydro ghost request targeted an unknown owner gas_cell_id");
      }
      if (internal::gasCellOwnerRankForLocalRow(
              context.state,
              *row,
              particle_row_by_id,
              "hydro ghost payload response") != static_cast<std::uint32_t>(m_mpi_context.worldRank())) {
        throw std::runtime_error("hydro ghost request targeted a local non-authoritative gas cell");
      }
      if (*row >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro ghost payload response dense row is outside geometry map");
      }
      const hydro::HydroConservedState state = m_conserved.loadCell(m_geometry_row_by_dense_row[*row]);
      local_payloads.push_back(parallel::HydroGhostCellPayloadRecord{
          .descriptor = request.descriptor,
          .face_key = request.face_key,
          .mass_density_comoving = state.mass_density_comoving,
          .momentum_density_x_comoving = state.momentum_density_x_comoving,
          .momentum_density_y_comoving = state.momentum_density_y_comoving,
          .momentum_density_z_comoving = state.momentum_density_z_comoving,
          .total_energy_density_comoving = state.total_energy_density_comoving,
          .metal_mass_density_comoving = state.metal_mass_density_comoving});
    }
    const auto global_payloads = parallel::executeBlockingHydroGhostCellPayloadExchange(
        m_mpi_context,
        local_payloads,
        hydro_sync_epoch);

    std::unordered_set<std::uint64_t> received_face_keys;
    for (const parallel::HydroGhostCellPayloadRecord& payload : global_payloads) {
      if (payload.descriptor.consumer_rank != m_mpi_context.worldRank()) {
        continue;
      }
      const auto slot_it = ghost_slot_by_face_key.find(payload.face_key);
      if (slot_it == ghost_slot_by_face_key.end()) {
        ++report.stale_or_invalid_payloads;
        continue;
      }
      if (!received_face_keys.insert(payload.face_key).second) {
        throw std::runtime_error("hydro ghost payload exchange returned duplicate payload for one face");
      }
      hydro::HydroGhostCell& ghost = m_geometry.ghost_cells.at(slot_it->second);
      if (payload.descriptor.gas_cell_id != ghost.origin_gas_cell_id ||
          payload.descriptor.owner_rank != ghost.origin_rank ||
          payload.descriptor.hydro_sync_epoch != ghost.hydro_sync_epoch ||
          payload.descriptor.decomposition_epoch != ghost.decomposition_epoch) {
        ++report.stale_or_invalid_payloads;
        continue;
      }
      m_conserved.storeCell(ghost.ghost_cell, hydro::HydroConservedState{
          .mass_density_comoving = payload.mass_density_comoving,
          .momentum_density_x_comoving = payload.momentum_density_x_comoving,
          .momentum_density_y_comoving = payload.momentum_density_y_comoving,
          .momentum_density_z_comoving = payload.momentum_density_z_comoving,
          .total_energy_density_comoving = payload.total_energy_density_comoving,
          .metal_mass_density_comoving = payload.metal_mass_density_comoving});
    }
    if (received_face_keys.size() != requests.size()) {
      throw std::runtime_error("hydro remote boundary ghost refresh did not receive every requested ghost cell");
    }
    report.requested_cells = requests.size();
    report.received_cells = received_face_keys.size();
    report.imported_ghosts = received_face_keys.size();
    report.interface_faces = requests.size();
    report.request_bytes = static_cast<std::uint64_t>(requests.size() * sizeof(parallel::HydroGhostCellRequest));
    report.payload_bytes =
        static_cast<std::uint64_t>(local_payloads.size() * sizeof(parallel::HydroGhostCellPayloadRecord));
    return report;
  }

  // Directed AMR boundary-cell exchange (sparse remote ghost sources) shared
  // by the synchronized distributed sweep and every level clock of the
  // distributed recursive driver. level_filter empty = all local patches.
  // Collective; admission failures are agreed before any payload moves.
  struct DirectedHydroGhostExchange {
    std::vector<parallel::AmrPatchPayloadRecord> local_patch_records;
    std::vector<amr::DistributedAmrRemotePatch> remote_patches;
    std::unordered_map<std::uint64_t, int> owner_rank_by_patch_id;
    parallel::DirectedAmrExchangeDiagnostics diagnostics{};
    core::MemoryReservation transport_reservation;
    std::vector<core::MemoryReservation> sparse_ghost_reservations;
  };

  [[nodiscard]] DirectedHydroGhostExchange exchangeDirectedHydroGhosts(
      core::StepContext& context,
      const amr::ProductionAmrHydroOptions& amr_options,
      std::span<const int> level_filter,
      std::uint64_t sequence) {
    DirectedHydroGhostExchange result;
    result.local_patch_records =
        internal::buildMigrationAmrPatchPayloadRecords(
            context.state, m_mpi_context.worldRank());
    if (!level_filter.empty()) {
      // G1.1: only patches whose state is current at this level clock take
      // part (every rank applies the same rule, so requests stay matched).
      std::erase_if(result.local_patch_records, [level_filter](const parallel::AmrPatchPayloadRecord& record) {
        return std::find(level_filter.begin(), level_filter.end(), static_cast<int>(record.level)) ==
            level_filter.end();
      });
    }
    const auto& local_patch_records = result.local_patch_records;
    parallel::DirectedAmrExchangeDiagnostics& directed_amr_diagnostics = result.diagnostics;
    const std::size_t wire_record_bytes = sizeof(parallel::AmrPatchCellPayloadRecord);
    std::size_t transport_round_limit_bytes = parallel::mpiTransportRoundLimitBytes();
    core::MemoryReservation& transport_reservation = result.transport_reservation;
    std::exception_ptr transport_admission_failure;
    try {
      if (m_memory_governor != nullptr) {
        const std::uint64_t simultaneous_record_bytes =
            core::checkedMemoryBytesAdd(
                static_cast<std::uint64_t>(wire_record_bytes),
                static_cast<std::uint64_t>(wire_record_bytes),
                "hydro remote ghost simultaneous wire-record bytes");
        const std::uint64_t requested_records = static_cast<std::uint64_t>(
            transport_round_limit_bytes / wire_record_bytes);
        const core::DeterministicBatchSizingResult sizing =
            core::selectDeterministicBatchSize(
                m_memory_governor->snapshot(),
                core::DeterministicBatchSizingPolicy{
                    .requested_max_items = std::max<std::uint64_t>(requested_records, 1U),
                    .bytes_per_item = simultaneous_record_bytes,
                    .fixed_reserve_bytes = 0U,
                    .minimum_items = 1U,
                    .alignment_items = 1U,
                    .headroom_use_basis_points = 10000U,
                });
        if (sizing.selected_items == 0U) {
          throw std::runtime_error(
              "hydro remote ghost transport cannot admit one send/receive wire record per rank");
        }
        transport_round_limit_bytes = core::checkedIntegralNarrow<std::size_t>(
            sizing.selected_items * static_cast<std::uint64_t>(wire_record_bytes),
            "hydro remote ghost selected transport round bytes");
        const std::uint64_t transport_reservation_bytes = sizing.selected_bytes;
        transport_reservation = m_memory_governor->reserve(
            core::MemoryClass::kCommunication,
            transport_reservation_bytes,
            "hydro.remote_ghost.transport");
        transport_reservation.commit();
      }
    } catch (...) {
      transport_admission_failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        transport_admission_failure,
        "hydro remote ghost transport memory admission");

    std::unordered_map<std::uint64_t, amr::DistributedAmrRemotePatch> remote_patch_store;
    std::vector<core::MemoryReservation>& sparse_ghost_reservations = result.sparse_ghost_reservations;
    std::uint64_t sparse_ghost_capacity_bytes = 0U;
    std::uint64_t sparse_ghost_metadata_bytes = 0U;

    const auto sparse_admission =
        [&](int peer_rank,
            std::span<const parallel::AmrPatchPayloadRecord> remote_records,
            std::span<const parallel::AmrPatchBoundaryCellRequest> requests,
            std::uint64_t logical_receive_count) {
          std::unordered_map<std::uint64_t, const parallel::AmrPatchPayloadRecord*> record_by_id;
          record_by_id.reserve(remote_records.size());
          for (const auto& record : remote_records) {
            record_by_id.emplace(record.patch_id, &record);
          }
          std::uint64_t expected_total = 0U;
          std::vector<std::pair<const parallel::AmrPatchPayloadRecord*, std::size_t>> planned;
          planned.reserve(requests.size());
          for (const auto& request : requests) {
            const auto record_it = record_by_id.find(request.patch_id);
            if (record_it == record_by_id.end()) {
              throw std::runtime_error(
                  "directed AMR sparse admission references missing remote patch metadata");
            }
            const auto* patch_record = record_it->second;
            if (patch_record->owner_rank != peer_rank) {
              throw std::runtime_error(
                  "directed AMR sparse admission patch owner does not match peer rank");
            }
            const std::size_t expected = parallel::directedAmrPatchBoundaryCellCount(
                *patch_record, request);
            expected_total = core::checkedMemoryBytesAdd(
                expected_total,
                static_cast<std::uint64_t>(expected),
                "directed AMR sparse admission expected record count");
            planned.emplace_back(patch_record, expected);
          }
          if (expected_total != logical_receive_count) {
            throw std::runtime_error(
                "directed AMR sparse admission logical receive count mismatch");
          }
          const std::uint64_t requested_sparse_bytes =
              core::checkedIntegralNarrow<std::uint64_t>(
                  core::checkedSizeMultiply(
                      core::checkedIntegralNarrow<std::size_t>(
                          logical_receive_count,
                          "hydro sparse ghost receive count"),
                      sizeof(amr::AmrHydroSparseRemoteCell),
                      "hydro sparse ghost requested bytes"),
                  "hydro sparse ghost requested bytes uint64");
          core::MemoryReservation sparse_reservation;
          if (m_memory_governor != nullptr) {
            sparse_reservation = m_memory_governor->reserve(
                core::MemoryClass::kCommunication,
                requested_sparse_bytes,
                "hydro.remote_ghost.sparse_store");
          }
          std::uint64_t actual_peer_capacity_bytes = 0U;
          for (const auto& [patch_record, expected] : planned) {
            amr::DistributedAmrRemotePatch remote;
            remote.patch = amr::PatchDescriptor{
                .patch_id = patch_record->patch_id,
                .parent_patch_id = patch_record->parent_patch_id,
                .level = static_cast<std::uint8_t>(patch_record->level),
                .morton_key = patch_record->morton_key,
                .origin_comov = {patch_record->origin_x_comoving, patch_record->origin_y_comoving, patch_record->origin_z_comoving},
                .extent_comov = {patch_record->extent_x_comoving, patch_record->extent_y_comoving, patch_record->extent_z_comoving},
                .cell_dims = {patch_record->cell_dim_x, patch_record->cell_dim_y, patch_record->cell_dim_z}};
            remote.owner_rank = patch_record->owner_rank;
            remote.ghost_hydro_epoch = patch_record->decomposition_epoch;
            remote.expected_ghost_hydro_epoch = patch_record->decomposition_epoch;
            remote.boundary_cells.reserve(expected);
            const std::uint64_t capacity_bytes =
                core::checkedIntegralNarrow<std::uint64_t>(
                    core::checkedSizeMultiply(
                        remote.boundary_cells.capacity(),
                        sizeof(amr::AmrHydroSparseRemoteCell),
                        "hydro sparse remote patch capacity bytes"),
                    "hydro sparse remote patch capacity bytes uint64");
            sparse_ghost_capacity_bytes = core::checkedMemoryBytesAdd(
                sparse_ghost_capacity_bytes,
                capacity_bytes,
                "hydro accumulated sparse ghost capacity bytes");
            actual_peer_capacity_bytes = core::checkedMemoryBytesAdd(
                actual_peer_capacity_bytes,
                capacity_bytes,
                "hydro peer sparse ghost capacity bytes");
            const auto [it, inserted] = remote_patch_store.emplace(
                remote.patch.patch_id, std::move(remote));
            if (!inserted) {
              throw std::runtime_error(
                  "directed AMR sparse admission found duplicate remote patch metadata");
            }
          }
          if (m_memory_governor != nullptr) {
            if (actual_peer_capacity_bytes > requested_sparse_bytes) {
              auto top_up = m_memory_governor->reserve(
                  core::MemoryClass::kCommunication,
                  actual_peer_capacity_bytes - requested_sparse_bytes,
                  "hydro.remote_ghost.sparse_store_capacity_top_up");
              top_up.commit();
              sparse_ghost_reservations.push_back(std::move(top_up));
            }
            sparse_reservation.commit();
            sparse_ghost_reservations.push_back(std::move(sparse_reservation));
          }
        };

    const auto sparse_consumer =
        [&](int peer_rank, std::span<const parallel::AmrPatchCellPayloadRecord> records) {
          for (const auto& cell_record : records) {
            auto remote_it = remote_patch_store.find(cell_record.patch_id);
            if (remote_it == remote_patch_store.end()) {
              throw std::runtime_error(
                  "directed AMR sparse consumer received a cell before patch admission");
            }
            amr::DistributedAmrRemotePatch& remote = remote_it->second;
            if (remote.owner_rank != peer_rank ||
                cell_record.owner_rank != peer_rank) {
              throw std::runtime_error(
                  "directed AMR sparse consumer received stale remote ownership metadata");
            }
            if (remote.boundary_cells.size() >= remote.boundary_cells.capacity()) {
              throw std::runtime_error(
                  "directed AMR sparse consumer exceeded its governor-admitted boundary capacity");
            }
            const double metallicity_mass_fraction = cell_record.mass_code > 0.0
                ? std::clamp(
                      cell_record.metal_mass_code / cell_record.mass_code,
                      0.0,
                      1.0)
                : 0.0;
            remote.boundary_cells.push_back(amr::AmrHydroSparseRemoteCell{
                .patch_local_cell = cell_record.local_cell_offset,
                .gas_cell_id = cell_record.gas_cell_id,
                .conserved = hydro::HydroCoreSolver::conservedFromPrimitive(
                    hydro::HydroPrimitiveState{
                        .rho_comoving = cell_record.density_code,
                        .vel_x_peculiar = cell_record.velocity_x_peculiar,
                        .vel_y_peculiar = cell_record.velocity_y_peculiar,
                        .vel_z_peculiar = cell_record.velocity_z_peculiar,
                        .pressure_comoving = cell_record.pressure_code,
                        .metallicity_mass_fraction = metallicity_mass_fraction},
                    k_gamma_adiabatic),
                .internal_energy_code = cell_record.internal_energy_code > 0.0
                    ? cell_record.internal_energy_code
                    : std::numeric_limits<double>::quiet_NaN()});
          }
        };

    const auto directed_amr_exchange = parallel::executeBlockingDirectedAmrPatchPayloadExchange(
        m_mpi_context,
        local_patch_records,
        [&context, this](
            std::span<const parallel::AmrPatchBoundaryCellRequest> requests,
            std::uint64_t first_record,
            std::size_t max_records,
            std::vector<parallel::AmrPatchCellPayloadRecord>& output) {
          internal::fillMigrationAmrPatchBoundaryCellPayloadChunk(
              context.state,
              m_mpi_context.worldRank(),
              requests,
              first_record,
              max_records,
              output);
        },
        sparse_consumer,
        sparse_admission,
        transport_round_limit_bytes,
        sequence,
        parallel::AmrInterfaceTopologyOptions{
            .periodic_box_extent_comoving =
                amr_options.physical_boundary_kind == hydro::HydroBoundaryKind::kPeriodic
                    ? amr_options.periodic_box_extent_comoving
                    : std::array<double, 3>{0.0, 0.0, 0.0},
            .ghost_layers = amr::productionAmrHydroGhostLayers(amr_options)});
    directed_amr_diagnostics = directed_amr_exchange.diagnostics;

    std::unordered_map<std::uint64_t, int>& owner_rank_by_patch_id = result.owner_rank_by_patch_id;
    owner_rank_by_patch_id.reserve(local_patch_records.size() + directed_amr_exchange.patch_payloads_received.size());
    for (const auto& record : local_patch_records) {
      const auto [it, inserted] = owner_rank_by_patch_id.emplace(record.patch_id, record.owner_rank);
      if (!inserted && it->second != record.owner_rank) {
        throw std::runtime_error("directed AMR hydro found conflicting local patch ownership metadata");
      }
    }

    std::vector<amr::DistributedAmrRemotePatch>& remote_patches = result.remote_patches;
    remote_patches.reserve(directed_amr_exchange.patch_payloads_received.size());
    for (const auto& patch_record : directed_amr_exchange.patch_payloads_received) {
      if (patch_record.owner_rank == m_mpi_context.worldRank()) {
        throw std::runtime_error("directed AMR hydro imported a local patch as a remote ghost");
      }
      owner_rank_by_patch_id.emplace(patch_record.patch_id, patch_record.owner_rank);
      auto remote_it = remote_patch_store.find(patch_record.patch_id);
      if (remote_it == remote_patch_store.end()) {
        throw std::runtime_error(
            "directed AMR hydro remote interface patch has no sparse boundary store");
      }
      auto& cells = remote_it->second.boundary_cells;
      std::sort(
          cells.begin(), cells.end(),
          [](const auto& lhs, const auto& rhs) {
            return lhs.patch_local_cell < rhs.patch_local_cell;
          });
      for (std::size_t i = 1; i < cells.size(); ++i) {
        if (cells[i - 1U].patch_local_cell == cells[i].patch_local_cell) {
          throw std::runtime_error(
              "directed AMR hydro received duplicate sparse remote boundary cell");
        }
      }
      remote_patches.push_back(std::move(remote_it->second));
    }
    remote_patch_store.clear();
    remote_patch_store.rehash(0U);
    sparse_ghost_metadata_bytes = core::checkedIntegralNarrow<std::uint64_t>(
        core::checkedSizeMultiply(
            remote_patches.capacity(),
            sizeof(amr::DistributedAmrRemotePatch),
            "hydro sparse remote patch metadata bytes"),
        "hydro sparse remote patch metadata bytes uint64");
    const std::uint64_t remote_ghost_resident_high_water =
        core::checkedMemoryBytesAdd(
            core::checkedMemoryBytesAdd(
                directed_amr_diagnostics.communication_workspace_high_water_bytes,
                sparse_ghost_capacity_bytes,
                "hydro remote ghost transport plus sparse cell capacity"),
            sparse_ghost_metadata_bytes,
            "hydro remote ghost total resident high-water");
    m_remote_ghost_sparse_cell_high_water = std::max(
        m_remote_ghost_sparse_cell_high_water,
        directed_amr_diagnostics.directed_patch_cell_records_received);
    m_remote_ghost_capacity_high_water_bytes = std::max(
        m_remote_ghost_capacity_high_water_bytes,
        sparse_ghost_capacity_bytes);
    m_remote_ghost_metadata_high_water_bytes = std::max(
        m_remote_ghost_metadata_high_water_bytes,
        sparse_ghost_metadata_bytes);
    m_communication_high_water_bytes = std::max(
        m_communication_high_water_bytes,
        remote_ghost_resident_high_water);
    m_last_communication_traffic_bytes =
        directed_amr_diagnostics.patch_cell_payload_bytes;
    return result;
  }

  // G1.1: workflow-side collectives of the distributed recursive AMR driver.
  // Retained remote coarse intervals: for each level with an open interval,
  // the remote boundary cells (stable patch_local_cell order) at the
  // interval start and end, exchanged at the level clock. Their storage is
  // admitted through the MemoryGovernor ("hydro.recursive.remote_interval").
  class DistributedRecursiveHooks final : public amr::RecursiveAmrDistributedHooks {
   public:
    DistributedRecursiveHooks(HydroAmrRuntimeImpl& runtime,
                              core::StepContext& context,
                              const amr::ProductionAmrHydroOptions& options,
                              std::vector<amr::PatchDescriptor>& regrid_neighbours,
                              std::vector<std::pair<std::uint64_t, int>>* remote_cell_owner)
        : m_runtime(runtime),
          m_context(context),
          m_options(options),
          m_regrid_neighbours(regrid_neighbours),
          m_remote_cell_owner(remote_cell_owner),
          m_sequence(context.integrator_state.step_index) {}

    void reduceFailure(std::exception_ptr local_failure, const char* stage) override {
      m_runtime.m_mpi_context.rethrowCollectivePreparationFailure(local_failure, stage);
    }

    void globalLevelUnion(std::array<std::uint64_t, 4>& level_bits) override {
      std::array<std::uint64_t, 256> counts{};
      for (std::size_t level = 0; level < counts.size(); ++level) {
        counts[level] = (level_bits[level / 64U] >> (level % 64U)) & 1ULL;
      }
      m_runtime.m_mpi_context.allreduceSumUint64sInPlace(counts);
      level_bits = {};
      for (std::size_t level = 0; level < counts.size(); ++level) {
        if (counts[level] != 0U) level_bits[level / 64U] |= 1ULL << (level % 64U);
      }
    }

    [[nodiscard]] std::uint64_t globalMax(std::uint64_t value) override {
      return m_runtime.m_mpi_context.allreduceMaxUint64(value);
    }

    [[nodiscard]] std::span<const amr::PatchDescriptor> prepare(double /*t_begin*/) override {
      DirectedHydroGhostExchange exchange =
          m_runtime.exchangeDirectedHydroGhosts(m_context, m_options, {}, nextSequence());
      std::exception_ptr failure;
      try {
        absorb(exchange);
        m_neighbours.clear();
        m_neighbours.reserve(exchange.remote_patches.size());
        for (const amr::DistributedAmrRemotePatch& remote : exchange.remote_patches) {
          m_neighbours.push_back(remote.patch);
          m_regrid_neighbours.push_back(remote.patch);
        }
        m_remote_patch_count = exchange.remote_patches.size();
      } catch (...) {
        failure = std::current_exception();
      }
      reduceFailure(failure, "distributed recursive AMR preparation");
      return m_neighbours;
    }

    [[nodiscard]] std::span<const amr::DistributedAmrRemotePatch> beginLevelSweep(
        int level, int coarser_level, int finer_level, double t_begin, double t_end, bool capture_start) override {
      std::array<int, 2> filter{level, finer_level};
      DirectedHydroGhostExchange exchange = m_runtime.exchangeDirectedHydroGhosts(
          m_context, m_options, std::span<const int>(filter.data(), finer_level >= 0 ? 2U : 1U), nextSequence());
      std::exception_ptr failure;
      try {
        absorb(exchange);
        m_sweep_sources.clear();
        m_sweep_sources.reserve(exchange.remote_patches.size() +
                                (coarser_level >= 0 ? intervalStore(coarser_level).size() : 0U));
        for (amr::DistributedAmrRemotePatch& remote : exchange.remote_patches) {
          remote.state_time_code = t_begin;  // levels {level, finer}: current at t_begin
          if (capture_start && static_cast<int>(remote.patch.level) == level) {
            RetainedInterval interval;
            interval.patch = remote.patch;
            interval.owner_rank = remote.owner_rank;
            interval.ghost_hydro_epoch = remote.ghost_hydro_epoch;
            interval.t_start = t_begin;
            interval.t_end = t_end;
            interval.start = remote.boundary_cells;
            admitRetained(interval.start.capacity());
            intervalStore(level).push_back(std::move(interval));
          }
          m_sweep_sources.push_back(std::move(remote));
        }
        if (coarser_level >= 0) {
          // Interpolated copies live for this sweep only; admitted per sweep.
          std::size_t interpolated_cells = 0U;
          for (const RetainedInterval& interval : intervalStore(coarser_level)) {
            interpolated_cells = core::checkedSizeAdd(interpolated_cells, interval.start.size(),
                                                      "recursive interpolated remote cells");
          }
          m_sweep_reservation.release();
          if (m_runtime.m_memory_governor != nullptr && interpolated_cells != 0U) {
            m_sweep_reservation = m_runtime.m_memory_governor->reserve(
                core::MemoryClass::kCommunication,
                core::checkedIntegralNarrow<std::uint64_t>(
                    core::checkedSizeMultiply(interpolated_cells, sizeof(amr::AmrHydroSparseRemoteCell),
                                              "recursive interpolated remote bytes"),
                    "recursive interpolated remote bytes"),
                "hydro.recursive.remote_interpolated");
            m_sweep_reservation.commit();
          }
          for (const RetainedInterval& interval : intervalStore(coarser_level)) {
            m_sweep_sources.push_back(interpolate(interval, t_begin));
          }
        }
      } catch (...) {
        failure = std::current_exception();
      }
      reduceFailure(failure, "distributed recursive AMR level ghost exchange");
      return m_sweep_sources;
    }

    void exchangeLevelFluxRegisters(core::SimulationState& state,
                                    std::span<const amr::FluxRegisterEntry> outbound,
                                    const amr::ProductionAmrHydroOptions& level_options) override {
      std::vector<parallel::AmrFluxRegisterPayloadRecord> payloads;
      std::exception_ptr failure;
      try {
        payloads.reserve(outbound.size());
        for (const amr::FluxRegisterEntry& entry : outbound) {
          const auto owner = m_owner_rank_by_patch_id.find(entry.coarse_patch_id);
          if (owner == m_owner_rank_by_patch_id.end()) {
            throw std::runtime_error("distributed recursive AMR register half has an unknown coarse owner");
          }
          payloads.push_back(fluxRegisterPayloadFromEntry(
              entry, m_runtime.m_mpi_context.worldRank(), owner->second, state));
        }
      } catch (...) {
        failure = std::current_exception();
      }
      reduceFailure(failure, "distributed recursive AMR register packing");
      const auto global = parallel::executeBlockingAmrFluxRegisterPayloadExchange(
          m_runtime.m_mpi_context, payloads, nextSequence(), m_runtime.m_memory_governor);
      try {
        std::vector<std::pair<int, std::size_t>> mine;
        for (std::size_t g = 0; g < global.size(); ++g) {
          if (global[g].owner_rank == m_runtime.m_mpi_context.worldRank()) mine.emplace_back(global[g].source_rank, g);
        }
        // Deterministic merge order: source rank, then wire order.
        std::stable_sort(mine.begin(), mine.end(),
                         [](const auto& l, const auto& r) { return l.first < r.first; });
        std::vector<std::pair<std::uint64_t, int>> received_keys;
        received_keys.reserve(mine.size());
        for (const auto& [rank, g] : mine) {
          received_keys.emplace_back(global[g].register_key, rank);
        }
        std::sort(received_keys.begin(), received_keys.end());
        for (std::size_t key = 1U; key < received_keys.size(); ++key) {
          if (received_keys[key - 1U] == received_keys[key]) {
            throw std::runtime_error(
                "distributed recursive AMR received a duplicate reflux register key from one source rank");
          }
        }
        std::vector<amr::FluxRegisterEntry> inbound;
        inbound.reserve(mine.size());
        for (const auto& [rank, g] : mine) inbound.push_back(fluxRegisterEntryFromPayload(global[g]));
        if (!inbound.empty()) (void)amr::mergeFluxRegistersIntoPendingStore(state, inbound, level_options);
        m_flux_records_sent += payloads.size();
        m_flux_records_received += inbound.size();
      } catch (...) {
        failure = std::current_exception();
      }
      reduceFailure(failure, "distributed recursive AMR register merge");
    }

    void captureLevelEnd(int level, int finer_level, double t_end) override {
      std::array<int, 2> filter{level, finer_level};
      DirectedHydroGhostExchange exchange = m_runtime.exchangeDirectedHydroGhosts(
          m_context, m_options, std::span<const int>(filter.data(), finer_level >= 0 ? 2U : 1U), nextSequence());
      std::exception_ptr failure;
      try {
        absorb(exchange);
        auto& store = intervalStore(level);
        std::size_t matched = 0U;
        for (amr::DistributedAmrRemotePatch& remote : exchange.remote_patches) {
          if (static_cast<int>(remote.patch.level) != level) continue;
          const auto it = std::find_if(store.begin(), store.end(), [&remote](const RetainedInterval& interval) {
            return interval.patch.patch_id == remote.patch.patch_id;
          });
          if (it == store.end()) {
            throw std::runtime_error(
                "distributed recursive AMR: remote coarse patch appeared at the interval end without a start state");
          }
          if (std::abs(it->t_end - t_end) > 1.0e-12 * std::max({1.0, std::abs(t_end), std::abs(it->t_end)})) {
            throw std::runtime_error("distributed recursive AMR: remote coarse interval end time mismatch");
          }
          it->end = std::move(remote.boundary_cells);
          admitRetained(it->end.capacity());
          ++matched;
        }
        if (matched != store.size()) {
          throw std::runtime_error(
              "distributed recursive AMR: a remote coarse patch lost its interval end state (fail closed)");
        }
      } catch (...) {
        failure = std::current_exception();
      }
      reduceFailure(failure, "distributed recursive AMR level end exchange");
    }

    void retireLevel(int level) noexcept override {
      if (level >= 0 && static_cast<std::size_t>(level) < m_intervals.size()) {
        m_intervals[static_cast<std::size_t>(level)].clear();
      }
    }

    [[nodiscard]] const parallel::DirectedAmrExchangeDiagnostics& diagnostics() const noexcept {
      return m_diagnostics;
    }
    [[nodiscard]] std::size_t remotePatchCount() const noexcept { return m_remote_patch_count; }
    [[nodiscard]] std::uint64_t fluxRecordsSent() const noexcept { return m_flux_records_sent; }
    [[nodiscard]] std::uint64_t fluxRecordsReceived() const noexcept { return m_flux_records_received; }

   private:
    struct RetainedInterval {
      amr::PatchDescriptor patch{};
      int owner_rank = 0;
      std::uint64_t ghost_hydro_epoch = 0U;
      double t_start = 0.0;
      double t_end = 0.0;
      std::vector<amr::AmrHydroSparseRemoteCell> start;
      std::vector<amr::AmrHydroSparseRemoteCell> end;
    };

    [[nodiscard]] std::uint64_t nextSequence() noexcept { return m_sequence++; }

    [[nodiscard]] std::vector<RetainedInterval>& intervalStore(int level) {
      if (level < 0 || level >= 256) throw std::out_of_range("distributed recursive AMR level outside 0..255");
      if (m_intervals.size() <= static_cast<std::size_t>(level)) m_intervals.resize(static_cast<std::size_t>(level) + 1U);
      return m_intervals[static_cast<std::size_t>(level)];
    }

    void admitRetained(std::size_t cell_capacity) {
      if (m_runtime.m_memory_governor == nullptr || cell_capacity == 0U) return;
      core::MemoryReservation reservation = m_runtime.m_memory_governor->reserve(
          core::MemoryClass::kCommunication,
          core::checkedIntegralNarrow<std::uint64_t>(
              core::checkedSizeMultiply(cell_capacity, sizeof(amr::AmrHydroSparseRemoteCell),
                                        "recursive remote interval bytes"),
              "recursive remote interval bytes"),
          "hydro.recursive.remote_interval");
      reservation.commit();
      m_retained_reservations.push_back(std::move(reservation));
    }

    // Linear interpolation in conserved variables with the ghost filler's
    // endpoint snapping; the result is labelled with the fill time.
    [[nodiscard]] amr::DistributedAmrRemotePatch interpolate(const RetainedInterval& interval, double t) const {
      if (interval.end.size() != interval.start.size()) {
        throw std::runtime_error("distributed recursive AMR: remote coarse interval start/end cell sets differ");
      }
      const double duration = interval.t_end - interval.t_start;
      if (!(duration > 0.0) || !std::isfinite(duration)) {
        throw std::runtime_error("distributed recursive AMR: invalid remote coarse interval");
      }
      double alpha = (t - interval.t_start) / duration;
      const double tol = 1.0e-12 * std::max({1.0, std::abs(interval.t_start), std::abs(interval.t_end), std::abs(t)});
      if (alpha < -tol || alpha > 1.0 + tol) {
        throw std::runtime_error("distributed recursive AMR: fill time outside the remote coarse interval");
      }
      if (alpha <= tol) alpha = 0.0;
      if (alpha >= 1.0 - tol) alpha = 1.0;
      amr::DistributedAmrRemotePatch out;
      out.patch = interval.patch;
      out.owner_rank = interval.owner_rank;
      out.ghost_hydro_epoch = interval.ghost_hydro_epoch;
      out.expected_ghost_hydro_epoch = interval.ghost_hydro_epoch;
      out.state_time_code = t;
      out.boundary_cells.reserve(interval.start.size());
      for (std::size_t c = 0; c < interval.start.size(); ++c) {
        const amr::AmrHydroSparseRemoteCell& a = interval.start[c];
        const amr::AmrHydroSparseRemoteCell& b = interval.end[c];
        if (a.patch_local_cell != b.patch_local_cell || a.gas_cell_id != b.gas_cell_id) {
          throw std::runtime_error("distributed recursive AMR: remote coarse identity changed inside the interval");
        }
        amr::AmrHydroSparseRemoteCell cell = a;
        cell.conserved = (1.0 - alpha) * a.conserved + alpha * b.conserved;
        // Specific internal energy interpolated with the same weights (NaN
        // stays NaN: the dual-energy lane then falls back to E - K, counted).
        cell.internal_energy_code = (1.0 - alpha) * a.internal_energy_code + alpha * b.internal_energy_code;
        if (!(cell.conserved.mass_density_comoving > 0.0) || !std::isfinite(cell.conserved.total_energy_density_comoving)) {
          throw std::runtime_error("distributed recursive AMR: interpolated remote coarse state is unphysical");
        }
        out.boundary_cells.push_back(cell);
      }
      return out;
    }

    void absorb(DirectedHydroGhostExchange& exchange) {
      for (const auto& [id, owner] : exchange.owner_rank_by_patch_id) m_owner_rank_by_patch_id.emplace(id, owner);
      if (m_remote_cell_owner != nullptr) {
        for (const amr::DistributedAmrRemotePatch& remote : exchange.remote_patches) {
          for (const amr::AmrHydroSparseRemoteCell& cell : remote.boundary_cells) {
            m_remote_cell_owner->emplace_back(cell.gas_cell_id, remote.owner_rank);
          }
        }
      }
      const auto& d = exchange.diagnostics;
      m_diagnostics.candidate_peer_count = std::max(m_diagnostics.candidate_peer_count, d.candidate_peer_count);
      m_diagnostics.neighbor_peer_count = std::max(m_diagnostics.neighbor_peer_count, d.neighbor_peer_count);
      m_diagnostics.remote_interface_count = std::max(m_diagnostics.remote_interface_count, d.remote_interface_count);
      m_diagnostics.directed_patch_descriptor_records_sent += d.directed_patch_descriptor_records_sent;
      m_diagnostics.directed_patch_descriptor_records_received += d.directed_patch_descriptor_records_received;
      m_diagnostics.directed_patch_cell_records_sent += d.directed_patch_cell_records_sent;
      m_diagnostics.directed_patch_cell_records_received += d.directed_patch_cell_records_received;
      m_diagnostics.control_plane_bytes += d.control_plane_bytes;
      m_diagnostics.patch_descriptor_bytes += d.patch_descriptor_bytes;
      m_diagnostics.patch_cell_payload_bytes += d.patch_cell_payload_bytes;
      m_diagnostics.communication_workspace_high_water_bytes = std::max(
          m_diagnostics.communication_workspace_high_water_bytes, d.communication_workspace_high_water_bytes);
      m_diagnostics.patch_cell_send_capacity_high_water_bytes = std::max(
          m_diagnostics.patch_cell_send_capacity_high_water_bytes, d.patch_cell_send_capacity_high_water_bytes);
      m_diagnostics.patch_cell_receive_capacity_high_water_bytes = std::max(
          m_diagnostics.patch_cell_receive_capacity_high_water_bytes, d.patch_cell_receive_capacity_high_water_bytes);
      m_diagnostics.patch_cell_transport_round_count += d.patch_cell_transport_round_count;
      m_diagnostics.remote_patch_ghost_count = std::max(m_diagnostics.remote_patch_ghost_count, d.remote_patch_ghost_count);
    }

    HydroAmrRuntimeImpl& m_runtime;
    core::StepContext& m_context;
    const amr::ProductionAmrHydroOptions& m_options;
    std::vector<amr::PatchDescriptor>& m_regrid_neighbours;
    std::vector<std::pair<std::uint64_t, int>>* m_remote_cell_owner = nullptr;
    std::uint64_t m_sequence = 0U;
    // Leases are declared before (destroyed after) the buffers they cover.
    std::vector<core::MemoryReservation> m_retained_reservations;
    core::MemoryReservation m_sweep_reservation;
    std::vector<amr::PatchDescriptor> m_neighbours;
    std::vector<amr::DistributedAmrRemotePatch> m_sweep_sources;
    std::vector<std::vector<RetainedInterval>> m_intervals;
    std::unordered_map<std::uint64_t, int> m_owner_rank_by_patch_id;
    parallel::DirectedAmrExchangeDiagnostics m_diagnostics{};
    std::size_t m_remote_patch_count = 0U;
    std::uint64_t m_flux_records_sent = 0U;
    std::uint64_t m_flux_records_received = 0U;
  };

  // G1.2: accepted-step guard for recursive subcycling, the same rule as the
  // coordinator's level-relaxed candidate: a level-L cell must be stable for
  // ITS substep dt / r^(L - L_min,global). Collective verdict (a violation on
  // any rank rejects on every rank before any exchange or mutation).
  void verifyAcceptedRecursiveAmrHydroCfl(core::StepContext& context, std::uint32_t refinement_ratio) {
    std::uint64_t local_min_level = std::numeric_limits<std::uint64_t>::max();
    for (std::size_t patch = 0; patch < context.state.patches.size(); ++patch) {
      if (context.state.patches.cell_count[patch] == 0U || context.state.patches.level[patch] < 0) continue;
      local_min_level = std::min<std::uint64_t>(local_min_level,
                                                static_cast<std::uint64_t>(context.state.patches.level[patch]));
    }
    const std::uint64_t global_min_level = m_mpi_context.allreduceMinUint64(local_min_level);
    std::exception_ptr failure;
    try {
      const double dt = context.integrator_state.dt_time_code;
      if (!(dt > 0.0)) throw std::invalid_argument("recursive AMR CFL guard requires dt_time_code > 0");
      const double cfl_scale_factor = m_config.units.coordinate_frame == core::CoordinateFrame::kComoving
          ? std::max(context.timeline_step.scale_factor_begin, 1.0e-12) : 1.0;
      core::HydroCflDiagnostics worst;
      bool have_worst = false;
      for (std::uint32_t row = 0; row < context.state.cells.size(); ++row) {
        const std::uint32_t patch = context.state.cells.patch_index[row];
        if (patch >= context.state.patches.size()) {
          throw std::runtime_error("recursive AMR CFL guard: gas cell without patch ownership");
        }
        const auto& p = context.state.patches;
        const std::array<double, 3> widths{
            p.extent_x_comoving[patch] / static_cast<double>(std::max<std::uint32_t>(p.cell_dim_x[patch], 1U)),
            p.extent_y_comoving[patch] / static_cast<double>(std::max<std::uint32_t>(p.cell_dim_y[patch], 1U)),
            p.extent_z_comoving[patch] / static_cast<double>(std::max<std::uint32_t>(p.cell_dim_z[patch], 1U))};
        const int offset = p.level[patch] - static_cast<int>(global_min_level);
        if (offset < 0) throw std::runtime_error("recursive AMR CFL guard: level below the global minimum");
        const double level_dt = dt / std::pow(static_cast<double>(refinement_ratio), static_cast<double>(offset));
        const core::DirectionalCflTimeStepInput input{
            .cell_width_axis_code = widths,
            .velocity_axis_code = {context.state.gas_cells.velocity_x_peculiar[row],
                                   context.state.gas_cells.velocity_y_peculiar[row],
                                   context.state.gas_cells.velocity_z_peculiar[row]},
            .sound_speed_code = std::max(context.state.gas_cells.sound_speed_code[row], 0.0),
            .coordinate_frame = m_config.units.coordinate_frame,
            .scale_factor = cfl_scale_factor,
            .stability_rule = core::hydroCflStabilityRuleFor(m_config),
            .active_axis = core::activeHydroAxes(widths, core::hydroDomainExtentCode(m_config)),
        };
        const core::HydroCflDiagnostics diagnostics = core::makeHydroCflDiagnostics(
            row, input, 0.4, level_dt, context.state.gas_cells.gas_cell_id[row], p.patch_id[patch]);
        if (!have_worst || diagnostics.safety_factor < worst.safety_factor) {
          worst = diagnostics;
          have_worst = true;
        }
      }
      if (have_worst) {
        m_last_hydro_cfl_diagnostics = worst;
        core::assertHydroCflStable(worst);
      }
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(failure, "recursive AMR accepted-step CFL guard");
  }

  void runProductionAmrHydroPath(core::StepContext& context) {
    // Hydro transport/source terms own the full accepted interval. Evaluate
    // background terms at its midpoint; AMR fine substeps derive their own
    // interval/evaluation epochs instead of inheriting the coarse end state.
    const double hydro_scale_factor = context.cosmology_background != nullptr
        ? context.timeline_step.scale_factor_midpoint
        : context.integrator_state.current_scale_factor;
    const bool cosmological_comoving =
        m_mode_policy.cosmological_comoving_frame &&
        m_config.units.coordinate_frame == core::CoordinateFrame::kComoving;
    const double hubble_rate_code = cosmological_comoving && context.cosmology_background != nullptr
        ? context.cosmology_background->hubbleSi(hydro_scale_factor) *
              context.timeline_step.time_si_per_code
        : 0.0;

    const hydro::HydroUpdateContext update{
        .dt_code = context.integrator_state.dt_time_code,
        .time_begin_code = context.timeline_step.time_begin_code,
        .time_end_code = context.timeline_step.time_end_code,
        .scale_factor_begin = context.timeline_step.scale_factor_begin,
        .scale_factor = std::max(1.0e-12, hydro_scale_factor),
        .scale_factor_end = context.timeline_step.scale_factor_end,
        .hubble_rate_begin_code = context.timeline_step.hubble_begin_code,
        .hubble_rate_code = hubble_rate_code,
        .hubble_rate_end_code = context.timeline_step.hubble_end_code,
        .comoving_coordinates =
            m_config.units.coordinate_frame == core::CoordinateFrame::kComoving,
    };
    const double density_frame_factor =
        m_config.units.coordinate_frame == core::CoordinateFrame::kComoving
            ? 1.0 / std::pow(update.scale_factor, 3.0)
            : 1.0;
    const double physical_density_cgs_per_density_code =
        m_units.densityCodeToCgs(density_frame_factor);
    hydro::HydroSourceContext source_context{
        .update = update,
        .gravity_accel_x_peculiar = {},
        .gravity_accel_y_peculiar = {},
        .gravity_accel_z_peculiar = {},
        .mass_density_physical_cgs = {},
        .hydrogen_number_density_cgs = {},
        .metallicity_mass_fraction = {},
        .temperature_k = {},
        .thermodynamic_closure = m_effective_ism_closure.get(),
        .redshift = std::max(0.0, (update.scale_factor > 0.0 ? (1.0 / update.scale_factor - 1.0) : 0.0)),
    };

    // Gravity and homogeneous expansion are integrated by the gravity-owned
    // KDK half-kicks around hydro/source evolution. Hydro owns transport plus
    // local thermochemical source terms only, preventing double application.
    std::vector<const hydro::HydroSourceTerm*> sources;
    if (m_cooling_source != nullptr) {
      sources.push_back(m_cooling_source.get());
    }
    if (m_effective_ism_energy_source != nullptr) {
      m_effective_ism_energy_source->resetLedger();
      sources.push_back(m_effective_ism_energy_source.get());
    }
    const hydro::HydroActiveBatchPolicy active_batch_policy = selectHydroActiveBatchPolicy();
    std::optional<core::MemoryReservation> active_batch_reservation;
    std::exception_ptr batch_admission_failure;
    try {
      active_batch_reservation = reserveHydroActiveBatch(active_batch_policy);
    } catch (...) {
      batch_admission_failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        batch_admission_failure, "hydro active-batch memory admission");

    if (m_config.numerics.amr_time_stepping == core::AmrTimeSteppingConfigMode::kRecursiveSubcycling) {
      verifyAcceptedRecursiveAmrHydroCfl(context, amr::ProductionAmrHydroOptions{}.refinement_ratio);
    }
    // T02: the adapter lives on this stack frame for the whole recursive
    // advance; the AMR driver holds only a non-owning pointer to it.
    std::optional<RecursiveAmrLevelSourceAdapter> recursive_source_adapter;
    if (m_config.numerics.amr_time_stepping ==
            core::AmrTimeSteppingConfigMode::kRecursiveSubcycling &&
        m_source_runtime != nullptr) {
      recursive_source_adapter.emplace(*m_source_runtime, context);
    }
    amr::ProductionAmrHydroOptions amr_options{
        .physical_boundary_kind = hydroBoundaryKindFromModePolicy(m_mode_policy.hydro_boundary),
        .adiabatic_index = k_gamma_adiabatic,
        .density_floor = m_config.numerics.hydro_density_floor_code,
        .pressure_floor = m_config.numerics.hydro_pressure_floor_code,
        .state_time_code = context.integrator_state.current_time_code,
        .ghost_fill_time_code = context.integrator_state.current_time_code,
        .active_batch_policy = active_batch_policy,
        .physical_density_cgs_per_density_code = physical_density_cgs_per_density_code,
        .regrid_memory_governor = m_memory_governor,
        .recursive_level_source_operator = recursive_source_adapter.has_value()
            ? &*recursive_source_adapter
            : nullptr,
        .reconstruction_scheme = hydroReconstructionSchemeFromConfig(m_config),
        .slope_limiter = hydroSlopeLimiterFromConfig(m_config),
        // Patch coordinates share the comoving coordinate frame of the
        // configured box (same convention as the fixed-grid layout).
        .periodic_box_extent_comoving = {
            m_config.cosmology.box_size_x_mpc_comoving,
            m_config.cosmology.box_size_y_mpc_comoving,
            m_config.cosmology.box_size_z_mpc_comoving},
        .thermodynamic_cache = amrThermodynamicCachePolicy(update.scale_factor, source_context.redshift)};
    amr::ProductionAmrHydroDiagnostics amr_diagnostics;
    std::size_t remote_patch_count = 0;
    std::size_t remote_flux_register_count = 0;
    std::size_t inbound_flux_register_count = 0;
    parallel::DirectedAmrExchangeDiagnostics directed_amr_diagnostics{};
    // Interface metadata of remote neighbour patches (CHUI Part 3 regrid
    // balance checks); filled by the distributed exchange below.
    std::vector<amr::PatchDescriptor> regrid_remote_neighbors;
    // CHUI Part 2 composition transport on the production AMR path. The
    // distributed path would need remote composition ghosts and composition
    // reflux payloads, which are Part 3 deliverables: refuse collectively
    // BEFORE any patch advances instead of silently dropping species.
    std::unique_ptr<physics::CompositionFluxAccumulator> composition_flux;
    // CHUI Part 3: distributed composition transport is supported through the
    // donor-fraction / remote-correction exchanges (Part 2 refused it here).
    bool distributed_composition = false;
    std::vector<std::pair<std::uint64_t, int>> remote_cell_owner;
    {
      const std::uint64_t composition_ranks = m_mpi_context.allreduceSumUint64(
          context.state.gas_composition.enabled() ? 1ULL : 0ULL);
      distributed_composition =
          composition_ranks > 0U && m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1;
      if (context.state.gas_composition.enabled()) {
        composition_flux = std::make_unique<physics::CompositionFluxAccumulator>(context.state);
        source_context.face_mass_flux_sink = composition_flux.get();
      }
    }
    if (m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1 &&
        m_config.numerics.amr_time_stepping == core::AmrTimeSteppingConfigMode::kRecursiveSubcycling) {
      // G1.1: distributed recursive Berger-Oliger subcycling. One directed
      // ghost exchange per level clock (plus the interval-end exchange of
      // levels with children), remote coarse ghosts temporally interpolated,
      // remote register halves merged into the coarse owner's pending store
      // after every sweep. Ranks without cells participate in every level.
      DistributedRecursiveHooks hooks(
          *this, context, amr_options, regrid_remote_neighbors,
          distributed_composition ? &remote_cell_owner : nullptr);
      amr_diagnostics = amr::advanceDistributedProductionAmrHydroRecursive(
          context.state,
          context.active_set.cell_indices,
          update,
          source_context,
          m_solver,
          *m_riemann_solver,
          sources,
          amr_options,
          m_mpi_context.worldRank(),
          context.state.gasCellIdentityGeneration(),
          hooks);
      directed_amr_diagnostics = hooks.diagnostics();
      directed_amr_diagnostics.outbound_reflux_count = hooks.fluxRecordsSent();
      directed_amr_diagnostics.inbound_reflux_count = hooks.fluxRecordsReceived();
      directed_amr_diagnostics.directed_flux_records_sent = hooks.fluxRecordsSent();
      directed_amr_diagnostics.directed_flux_records_received = hooks.fluxRecordsReceived();
      remote_patch_count = hooks.remotePatchCount();
      remote_flux_register_count = hooks.fluxRecordsSent();
      inbound_flux_register_count = hooks.fluxRecordsReceived();
      std::sort(remote_cell_owner.begin(), remote_cell_owner.end());
      remote_cell_owner.erase(std::unique(remote_cell_owner.begin(), remote_cell_owner.end()), remote_cell_owner.end());
      std::sort(regrid_remote_neighbors.begin(), regrid_remote_neighbors.end(),
                [](const amr::PatchDescriptor& l, const amr::PatchDescriptor& r) { return l.patch_id < r.patch_id; });
      regrid_remote_neighbors.erase(
          std::unique(regrid_remote_neighbors.begin(), regrid_remote_neighbors.end(),
                      [](const amr::PatchDescriptor& l, const amr::PatchDescriptor& r) { return l.patch_id == r.patch_id; }),
          regrid_remote_neighbors.end());
      m_communication_high_water_bytes = std::max(
          m_communication_high_water_bytes, directed_amr_diagnostics.communication_workspace_high_water_bytes);
      m_last_communication_traffic_bytes = directed_amr_diagnostics.patch_cell_payload_bytes;
    } else if (m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1) {
      DirectedHydroGhostExchange ghost_exchange =
          exchangeDirectedHydroGhosts(context, amr_options, {}, context.integrator_state.step_index);
      directed_amr_diagnostics = ghost_exchange.diagnostics;
      auto& owner_rank_by_patch_id = ghost_exchange.owner_rank_by_patch_id;
      std::vector<amr::DistributedAmrRemotePatch>& remote_patches = ghost_exchange.remote_patches;
      remote_patch_count = remote_patches.size();
      regrid_remote_neighbors.reserve(remote_patches.size());
      for (const amr::DistributedAmrRemotePatch& remote : remote_patches) {
        regrid_remote_neighbors.push_back(remote.patch);
        if (distributed_composition) {
          for (const amr::AmrHydroSparseRemoteCell& cell : remote.boundary_cells) {
            remote_cell_owner.emplace_back(cell.gas_cell_id, remote.owner_rank);
          }
        }
      }
      std::sort(remote_cell_owner.begin(), remote_cell_owner.end());
      std::vector<amr::FluxRegisterEntry> outbound_remote_entries;
      std::vector<amr::FluxRegisterEntry> deferred_local_entries;
      if (amr::hasProductionAmrHydroCoverage(context.state)) {
        amr_diagnostics = amr::advanceDistributedProductionAmrHydro(
            context.state,
            context.active_set.cell_indices,
            update,
            source_context,
            m_solver,
            *m_riemann_solver,
            sources,
            amr_options,
            amr::DistributedAmrHydroExchange{
                .local_rank = m_mpi_context.worldRank(),
                .ghost_hydro_epoch = context.state.gasCellIdentityGeneration(),
                .expected_ghost_hydro_epoch = context.state.gasCellIdentityGeneration(),
                .remote_patches = remote_patches,
                .outbound_remote_flux_registers = &outbound_remote_entries,
                .deferred_local_flux_registers = &deferred_local_entries});
      } else if (context.state.cells.size() != 0U || context.state.patches.size() != 0U) {
        throw std::runtime_error(
            "distributed AMR hydro rank has partial local patch coverage; only an explicitly empty rank may skip advancement");
      }
      std::vector<parallel::AmrFluxRegisterPayloadRecord> local_flux_payloads;
      local_flux_payloads.reserve(outbound_remote_entries.size());
      for (const amr::FluxRegisterEntry& entry : outbound_remote_entries) {
        const auto owner_it = owner_rank_by_patch_id.find(entry.coarse_patch_id);
        if (owner_it == owner_rank_by_patch_id.end()) {
          throw std::runtime_error("distributed AMR hydro remote reflux entry has unknown coarse patch owner");
        }
        local_flux_payloads.push_back(fluxRegisterPayloadFromEntry(
            entry, m_mpi_context.worldRank(), owner_it->second, context.state));
      }
      const auto global_flux_payloads = parallel::executeBlockingAmrFluxRegisterPayloadExchange(
          m_mpi_context, local_flux_payloads, context.integrator_state.step_index,
          m_memory_governor);
      directed_amr_diagnostics.outbound_reflux_count = static_cast<std::uint64_t>(local_flux_payloads.size());
      directed_amr_diagnostics.inbound_reflux_count = static_cast<std::uint64_t>(global_flux_payloads.size());
      directed_amr_diagnostics.directed_flux_records_sent = static_cast<std::uint64_t>(local_flux_payloads.size());
      directed_amr_diagnostics.directed_flux_records_received = static_cast<std::uint64_t>(global_flux_payloads.size());
      directed_amr_diagnostics.flux_payload_bytes =
          (static_cast<std::uint64_t>(local_flux_payloads.size()) + static_cast<std::uint64_t>(global_flux_payloads.size())) *
          sizeof(parallel::AmrFluxRegisterPayloadRecord);
      std::vector<amr::FluxRegisterEntry> inbound_entries;
      std::vector<int> inbound_source_ranks;
      inbound_entries.reserve(global_flux_payloads.size());
      inbound_source_ranks.reserve(global_flux_payloads.size());
      for (const auto& payload : global_flux_payloads) {
        if (payload.owner_rank != m_mpi_context.worldRank()) {
          continue;
        }
        inbound_entries.push_back(fluxRegisterEntryFromPayload(payload));
        inbound_source_ranks.push_back(payload.source_rank);
      }
      // Validate duplicate keys via sorted adjacent comparison instead of
      // a population-scale hash set.
      // One key may legitimately arrive from several source ranks (fine faces
      // of one coarse face on different ranks); the same key twice from the
      // same source rank is a protocol error.
      if (inbound_entries.size() > 1U) {
        std::vector<std::pair<std::uint64_t, int>> sort_keys;
        sort_keys.reserve(inbound_entries.size());
        for (std::size_t e = 0; e < inbound_entries.size(); ++e) {
          sort_keys.emplace_back(inbound_entries[e].register_key, inbound_source_ranks[e]);
        }
        std::sort(sort_keys.begin(), sort_keys.end());
        for (std::size_t i = 1; i < sort_keys.size(); ++i) {
          if (sort_keys[i - 1U] == sort_keys[i]) {
            throw std::runtime_error("distributed AMR hydro received duplicate reflux register key for local owner");
          }
        }
      }
      inbound_flux_register_count = inbound_entries.size();
      remote_flux_register_count = local_flux_payloads.size();
      // G2.3: complete every split interface once: local halves deferred by
      // the sweep first, then inbound halves in source-rank order, merged by
      // key. Previously each half was applied alone and skipped as incomplete,
      // so a coarse/fine interface split across ranks was never refluxed.
      std::vector<amr::FluxRegisterEntry> split_interface_entries;
      split_interface_entries.reserve(deferred_local_entries.size() + inbound_entries.size());
      split_interface_entries.insert(split_interface_entries.end(), deferred_local_entries.begin(),
                                     deferred_local_entries.end());
      {
        std::vector<std::size_t> order(inbound_entries.size());
        for (std::size_t e = 0; e < order.size(); ++e) order[e] = e;
        std::stable_sort(order.begin(), order.end(), [&](std::size_t l, std::size_t r) {
          return inbound_source_ranks[l] < inbound_source_ranks[r];
        });
        for (const std::size_t e : order) split_interface_entries.push_back(inbound_entries[e]);
      }
      const std::vector<amr::FluxRegisterEntry> merged_split_entries =
          amr::mergeFluxRegisterEntriesByKey(split_interface_entries);
      if (!merged_split_entries.empty()) {
        const auto local_descriptors = amr::buildProductionAmrPatchDescriptors(context.state);
        const amr::RefluxDiagnostics remote_reflux = amr::applyFluxRegistersToSimulationState(
            context.state, merged_split_entries, local_descriptors, k_gamma_adiabatic,
            amrThermodynamicCachePolicy(update.scale_factor, source_context.redshift));
        amr_diagnostics.reflux.complete_register_count += remote_reflux.complete_register_count;
        amr_diagnostics.reflux.skipped_incomplete_register_count += remote_reflux.skipped_incomplete_register_count;
        amr_diagnostics.reflux.skipped_area_mismatch_count += remote_reflux.skipped_area_mismatch_count;
        amr_diagnostics.reflux.skipped_missing_target_count += remote_reflux.skipped_missing_target_count;
        amr_diagnostics.reflux.corrected_cells += remote_reflux.corrected_cells;
        amr_diagnostics.reflux.corrected_mass_code += remote_reflux.corrected_mass_code;
        amr_diagnostics.reflux.corrected_momentum_x_code += remote_reflux.corrected_momentum_x_code;
        amr_diagnostics.reflux.corrected_momentum_y_code += remote_reflux.corrected_momentum_y_code;
        amr_diagnostics.reflux.corrected_momentum_z_code += remote_reflux.corrected_momentum_z_code;
        amr_diagnostics.reflux.corrected_total_energy_code += remote_reflux.corrected_total_energy_code;
        amr_diagnostics.reflux.corrected_metal_mass_code += remote_reflux.corrected_metal_mass_code;
        amr_diagnostics.reflux.corrected_energy_code += remote_reflux.corrected_energy_code;
        amr_diagnostics.reflux.corrected_internal_energy_code += remote_reflux.corrected_internal_energy_code;
        amr_diagnostics.reflux.dual_energy_reflux_cells += remote_reflux.dual_energy_reflux_cells;
        amr_diagnostics.reflux.dual_energy_reflux_thermal_correction_code +=
            remote_reflux.dual_energy_reflux_thermal_correction_code;
      }
    } else if (m_config.numerics.amr_time_stepping == core::AmrTimeSteppingConfigMode::kRecursiveSubcycling) {
      // CHUI Part 3 Milestone I: recursive multilevel subcycling of the gas
      // hierarchy inside the coordinator's accepted interval (the interval
      // belongs to the coarsest level; finer levels advance dt / r^dL).
      amr_diagnostics = amr::advanceProductionAmrHydroRecursive(
          context.state,
          context.active_set.cell_indices,
          update,
          source_context,
          m_solver,
          *m_riemann_solver,
          sources,
          amr_options);
    } else {
      amr_diagnostics = amr::advanceProductionAmrHydro(
          context.state,
          context.active_set.cell_indices,
          update,
          source_context,
          m_solver,
          *m_riemann_solver,
          sources,
          amr_options);
    }
    m_last_amr_hydro_diagnostics = amr_diagnostics;
    if (distributed_composition) {
      completeDistributedCompositionTransport(context, composition_flux.get(), remote_cell_owner);
    } else if (composition_flux != nullptr) {
      // Single rank: the composition update is prepared (may throw = hard
      // stage failure) and committed immediately after the hydro sweep.
      composition_flux->prepare(context.state);
      admitCompositionSparseActivation(
          context.state, *composition_flux, m_memory_governor,
          "hydro.amr.adaptive_chemistry_rows");
      m_last_composition_transport = composition_flux->commit(context.state);
      (void)physics::synchronizeMetalMirrorFromComposition(context.state, std::span<const std::uint32_t>{});
    }
    m_communication_high_water_bytes = std::max(
        m_communication_high_water_bytes,
        directed_amr_diagnostics.communication_workspace_high_water_bytes);

    // Effective-EOS pressure / signal speed and the temperature cache are
    // written by the AMR scatter and reflux themselves (amr_thermodynamic_cache
    // policy), for exactly the rows they mutate; the former whole-rank
    // closure pass here (which also read the stale pressure cache through a
    // max()) is no longer needed.
    internal::synchronizeParentParticleCompatibilityMirrors(
        context.state,
        m_current_world_rank,
        "production AMR hydro parent compatibility mirror");
    m_last_directed_amr_diagnostics = directed_amr_diagnostics;
    runAutomaticAmrRegrid(context, amr_options, regrid_remote_neighbors);
    // H01 telemetry: CTU interface outcomes summed over ranks (the per-rank
    // values above are summed over patches and recursive level clocks). The
    // reductions are rank-uniform and independent of profiler availability.
    const std::uint64_t global_ctu_corrected =
        m_mpi_context.allreduceSumUint64(amr_diagnostics.ctu_transverse_corrected_states);
    const std::uint64_t global_ctu_degraded =
        m_mpi_context.allreduceSumUint64(amr_diagnostics.ctu_transverse_degraded_states);
    const std::uint64_t global_ctu_ghost_corrected =
        m_mpi_context.allreduceSumUint64(amr_diagnostics.ctu_transverse_ghost_corrected_states);
    const std::uint64_t global_ctu_ghost_axes =
        m_mpi_context.allreduceSumUint64(amr_diagnostics.ctu_transverse_ghost_axes);
    const std::uint64_t global_ctu_boundary_proxy =
        m_mpi_context.allreduceSumUint64(amr_diagnostics.ctu_transverse_boundary_proxy_states);
    const std::uint64_t global_ctu_inadmissible =
        m_mpi_context.allreduceSumUint64(amr_diagnostics.ctu_transverse_inadmissible_states);
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.amr_production_stage",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.amr",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "advanced production AMR hydro patches from SimulationState-authoritative gas rows",
          .payload = {{"patch_count", std::to_string(amr_diagnostics.patch_count)},
                      {"advanced_patch_count", std::to_string(amr_diagnostics.advanced_patch_count)},
                      {"riemann_solve_count", std::to_string(amr_diagnostics.riemann_solve_count)},
                      {"positivity_fallback_count", std::to_string(amr_diagnostics.positivity_fallback_count)},
                      {"stencil_degradation_count", std::to_string(amr_diagnostics.stencil_degradation_count)},
                      {"ctu_transverse_degraded_states", std::to_string(amr_diagnostics.ctu_transverse_degraded_states)},
                      {"ctu_transverse_ghost_corrected_states", std::to_string(amr_diagnostics.ctu_transverse_ghost_corrected_states)},
                      {"ctu_transverse_ghost_axes", std::to_string(amr_diagnostics.ctu_transverse_ghost_axes)},
                      {"ctu_transverse_corrected_states", std::to_string(amr_diagnostics.ctu_transverse_corrected_states)},
                      {"ctu_transverse_boundary_proxy_states",
                       std::to_string(amr_diagnostics.ctu_transverse_boundary_proxy_states)},
                      {"ctu_transverse_inadmissible_states",
                       std::to_string(amr_diagnostics.ctu_transverse_inadmissible_states)},
                      {"global_ctu_transverse_corrected_states", std::to_string(global_ctu_corrected)},
                      {"global_ctu_transverse_degraded_states", std::to_string(global_ctu_degraded)},
                      {"global_ctu_transverse_ghost_corrected_states", std::to_string(global_ctu_ghost_corrected)},
                      {"global_ctu_transverse_ghost_axes", std::to_string(global_ctu_ghost_axes)},
                      {"global_ctu_transverse_boundary_proxy_states", std::to_string(global_ctu_boundary_proxy)},
                      {"global_ctu_transverse_inadmissible_states", std::to_string(global_ctu_inadmissible)},
                      {"dual_energy_auxiliary_cells", std::to_string(amr_diagnostics.dual_energy_auxiliary_cells)},
                      {"dual_energy_unresolved_ghosts", std::to_string(amr_diagnostics.dual_energy_unresolved_ghosts)},
                      {"dual_energy_thermal_correction_code",
                       std::to_string(amr_diagnostics.dual_energy_thermal_correction_code)},
                      {"dual_energy_reflux_cells", std::to_string(amr_diagnostics.reflux.dual_energy_reflux_cells)},
                      {"dual_energy_reflux_thermal_correction_code",
                       std::to_string(amr_diagnostics.reflux.dual_energy_reflux_thermal_correction_code)},
                      {"active_cell_count", std::to_string(amr_diagnostics.active_cell_count)},
                      {"active_face_count", std::to_string(amr_diagnostics.active_face_count)},
                      {"flux_register_entry_count", std::to_string(amr_diagnostics.flux_register_entry_count)},
                      {"reflux_corrected_cells", std::to_string(amr_diagnostics.reflux.corrected_cells)},
                      {"reflux_corrected_mass", std::to_string(amr_diagnostics.reflux.corrected_mass_code)},
                      {"reflux_corrected_momentum_x", std::to_string(amr_diagnostics.reflux.corrected_momentum_x_code)},
                      {"reflux_corrected_momentum_y", std::to_string(amr_diagnostics.reflux.corrected_momentum_y_code)},
                      {"reflux_corrected_momentum_z", std::to_string(amr_diagnostics.reflux.corrected_momentum_z_code)},
                      {"reflux_corrected_total_energy", std::to_string(amr_diagnostics.reflux.corrected_total_energy_code)},
                      {"reflux_corrected_metal_mass", std::to_string(amr_diagnostics.reflux.corrected_metal_mass_code)},
                      {"reflux_corrected_energy", std::to_string(amr_diagnostics.reflux.corrected_energy_code)},
                      {"reflux_corrected_internal_energy", std::to_string(amr_diagnostics.reflux.corrected_internal_energy_code)},
                      {"reflux_complete_register_count", std::to_string(amr_diagnostics.reflux.complete_register_count)},
                      {"reflux_skipped_incomplete_register_count", std::to_string(amr_diagnostics.reflux.skipped_incomplete_register_count)},
                      {"reflux_skipped_area_mismatch_count", std::to_string(amr_diagnostics.reflux.skipped_area_mismatch_count)},
                      {"reflux_skipped_missing_target_count", std::to_string(amr_diagnostics.reflux.skipped_missing_target_count)},
                      {"distributed_remote_patch_count", std::to_string(remote_patch_count)},
                      {"distributed_remote_flux_register_count", std::to_string(remote_flux_register_count)},
                      {"distributed_inbound_flux_register_count", std::to_string(inbound_flux_register_count)},
                      {"hydro_active_batch_capacity_cells", std::to_string(amr_diagnostics.active_batch_capacity_cells)},
                      {"hydro_face_batch_capacity", std::to_string(amr_diagnostics.face_batch_capacity)},
                      {"hydro_scratch_high_water_bytes", std::to_string(amr_diagnostics.scratch_high_water_bytes)},
                      {"hydro_geometry_capacity_bytes", std::to_string(amr_diagnostics.geometry_capacity_bytes)},
                      {"hydro_prepared_ghost_capacity_bytes", std::to_string(amr_diagnostics.prepared_ghost_capacity_bytes)},
                      {"hydro_max_patch_conserved_bytes", std::to_string(amr_diagnostics.max_patch_conserved_bytes)},
                      {"hydro_batch_scratch_bytes_per_cell", std::to_string(hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell)},
                      {"directed_amr_candidate_peer_count", std::to_string(directed_amr_diagnostics.candidate_peer_count)},
                      {"directed_amr_neighbor_peer_count", std::to_string(directed_amr_diagnostics.neighbor_peer_count)},
                      {"directed_amr_patch_descriptor_records_sent", std::to_string(directed_amr_diagnostics.directed_patch_descriptor_records_sent)},
                      {"directed_amr_patch_descriptor_records_received", std::to_string(directed_amr_diagnostics.directed_patch_descriptor_records_received)},
                      {"directed_amr_patch_cell_records_sent", std::to_string(directed_amr_diagnostics.directed_patch_cell_records_sent)},
                      {"directed_amr_patch_cell_records_received", std::to_string(directed_amr_diagnostics.directed_patch_cell_records_received)},
                      {"directed_amr_flux_records_sent", std::to_string(directed_amr_diagnostics.directed_flux_records_sent)},
                      {"directed_amr_flux_records_received", std::to_string(directed_amr_diagnostics.directed_flux_records_received)},
                      {"directed_amr_control_plane_bytes", std::to_string(directed_amr_diagnostics.control_plane_bytes)},
                      {"directed_amr_patch_descriptor_bytes", std::to_string(directed_amr_diagnostics.patch_descriptor_bytes)},
                      {"directed_amr_patch_cell_payload_bytes", std::to_string(directed_amr_diagnostics.patch_cell_payload_bytes)},
                      {"directed_amr_patch_cell_send_capacity_high_water_bytes", std::to_string(directed_amr_diagnostics.patch_cell_send_capacity_high_water_bytes)},
                      {"directed_amr_patch_cell_receive_capacity_high_water_bytes", std::to_string(directed_amr_diagnostics.patch_cell_receive_capacity_high_water_bytes)},
                      {"directed_amr_communication_workspace_high_water_bytes", std::to_string(directed_amr_diagnostics.communication_workspace_high_water_bytes)},
                      {"directed_amr_remote_ghost_sparse_cell_high_water", std::to_string(m_remote_ghost_sparse_cell_high_water)},
                      {"directed_amr_remote_ghost_capacity_high_water_bytes", std::to_string(m_remote_ghost_capacity_high_water_bytes)},
                      {"directed_amr_remote_ghost_metadata_high_water_bytes", std::to_string(m_remote_ghost_metadata_high_water_bytes)},
                      {"directed_amr_total_resident_high_water_bytes", std::to_string(m_communication_high_water_bytes)},
                      {"directed_amr_recent_communication_traffic_bytes", std::to_string(m_last_communication_traffic_bytes)},
                      {"directed_amr_patch_cell_transport_round_count", std::to_string(directed_amr_diagnostics.patch_cell_transport_round_count)},
                      {"directed_amr_flux_payload_bytes", std::to_string(directed_amr_diagnostics.flux_payload_bytes)},
                      {"directed_amr_remote_interface_count", std::to_string(directed_amr_diagnostics.remote_interface_count)},
                      {"directed_amr_inbound_reflux_count", std::to_string(directed_amr_diagnostics.inbound_reflux_count)},
                      {"directed_amr_outbound_reflux_count", std::to_string(directed_amr_diagnostics.outbound_reflux_count)}}});
    }
  }

  // CHUI Part 3 Checkpoint B: distributed composition transport. Every rank
  // participates in the three collective exchanges (possibly with empty
  // streams). Order: donor requests -> donor fractions (owner step-start
  // sidecar; no rank has committed composition yet) -> deferred prepare ->
  // remote corrections to owners -> validation -> commit. Any failure is
  // made collective before the commit, so no rank mutates composition alone.
  void completeDistributedCompositionTransport(
      core::StepContext& context,
      physics::CompositionFluxAccumulator* composition_flux,
      std::span<const std::pair<std::uint64_t, int>> remote_cell_owner) {
    const int rank = m_mpi_context.worldRank();
    const auto ownerOf = [&](std::uint64_t id) {
      const auto it = std::lower_bound(remote_cell_owner.begin(), remote_cell_owner.end(),
                                       std::pair<std::uint64_t, int>{id, std::numeric_limits<int>::min()});
      if (it == remote_cell_owner.end() || it->first != id) {
        throw std::runtime_error("distributed composition: remote cell " + std::to_string(id) +
                                 " has no known owner rank (fail closed)");
      }
      return it->second;
    };
    const std::uint64_t sequence = context.integrator_state.step_index;
    std::exception_ptr failure;
    std::vector<parallel::CompositionDonorRequestRecord> requests;
    try {
      if (composition_flux != nullptr) {
        for (const std::uint64_t id : composition_flux->requiredRemoteDonorIds()) {
          requests.push_back(parallel::CompositionDonorRequestRecord{
              .gas_cell_id = id, .requester_rank = rank, .owner_rank = ownerOf(id)});
        }
      }
    } catch (...) { failure = std::current_exception(); }
    m_mpi_context.rethrowCollectivePreparationFailure(failure, "distributed composition donor requests");
    const auto inbound_requests = parallel::exchangeCompositionDonorRequests(m_mpi_context, requests, sequence);
    std::vector<parallel::CompositionDonorFractionRecord> responses;
    try {
      const auto& sidecar = context.state.gas_composition;
      if (!inbound_requests.empty() && !sidecar.enabled()) {
        throw std::runtime_error("distributed composition: donor requested from a rank without composition");
      }
      const std::size_t lanes = sidecar.enabled() ? sidecar.laneCount() : 0U;
      if (lanes > parallel::k_max_composition_lanes) {
        throw std::runtime_error("distributed composition: lane count exceeds the wire record capacity");
      }
      const auto& d = sidecar.descriptor();
      for (const auto& request : inbound_requests) {
        const auto row = context.state.rowForGasCellId(request.gas_cell_id);
        if (!row.has_value()) {
          throw std::runtime_error("distributed composition: donor request for a cell this rank does not own");
        }
        double gas_mass = 0.0;
        for (std::uint32_t l = 0; l < d.element_lane_count; ++l) gas_mass += sidecar.value(d.element_lane_begin + l, *row);
        for (std::uint32_t l = 0; l < d.dust_lane_count; ++l) gas_mass += sidecar.value(d.dust_lane_begin + l, *row);
        if (!(gas_mass > 0.0)) {
          throw std::runtime_error("distributed composition: requested donor has no composition mass");
        }
        parallel::CompositionDonorFractionRecord response{
            .gas_cell_id = request.gas_cell_id, .owner_rank = rank, .requester_rank = request.requester_rank,
            .lane_count = static_cast<std::uint32_t>(lanes)};
        for (std::size_t l = 0; l < lanes; ++l) response.fraction[l] = sidecar.value(l, *row) / gas_mass;
        responses.push_back(response);
      }
    } catch (...) { failure = std::current_exception(); }
    m_mpi_context.rethrowCollectivePreparationFailure(failure, "distributed composition donor responses");
    auto fractions = parallel::exchangeCompositionDonorFractions(m_mpi_context, responses, sequence);
    std::vector<parallel::CompositionCorrectionRecord> corrections;
    try {
      if (composition_flux != nullptr) {
        std::sort(fractions.begin(), fractions.end(),
                  [](const auto& l, const auto& r) { return l.gas_cell_id < r.gas_cell_id; });
        const std::size_t lanes = context.state.gas_composition.laneCount();
        std::vector<std::uint64_t> ids;
        std::vector<double> values;
        for (const auto& record : fractions) {
          if (record.lane_count != lanes) throw std::runtime_error("distributed composition: donor lane mismatch");
          ids.push_back(record.gas_cell_id);
          values.insert(values.end(), record.fraction.begin(), record.fraction.begin() + lanes);
        }
        composition_flux->setRemoteDonorFractions(ids, values, lanes);
        composition_flux->prepareDeferred(context.state);
        for (const auto& outbound : composition_flux->outboundRemoteCorrections()) {
          parallel::CompositionCorrectionRecord record{
              .gas_cell_id = outbound.gas_cell_id, .source_rank = rank, .owner_rank = ownerOf(outbound.gas_cell_id),
              .lane_count = static_cast<std::uint32_t>(outbound.lane_delta_code.size())};
          std::copy(outbound.lane_delta_code.begin(), outbound.lane_delta_code.end(), record.lane_delta_code.begin());
          corrections.push_back(record);
        }
      }
    } catch (...) { failure = std::current_exception(); }
    m_mpi_context.rethrowCollectivePreparationFailure(failure, "distributed composition deferred prepare");
    auto inbound_corrections = parallel::exchangeCompositionCorrections(m_mpi_context, corrections, sequence);
    try {
      // Deterministic application order independent of arrival order.
      std::stable_sort(inbound_corrections.begin(), inbound_corrections.end(), [](const auto& l, const auto& r) {
        return l.gas_cell_id != r.gas_cell_id ? l.gas_cell_id < r.gas_cell_id : l.source_rank < r.source_rank;
      });
      if (!inbound_corrections.empty() && composition_flux == nullptr) {
        throw std::runtime_error("distributed composition: correction delivered to a rank without composition");
      }
      for (const auto& record : inbound_corrections) {
        composition_flux->addInboundCorrection(
            context.state, record.gas_cell_id,
            std::span<const double>(record.lane_delta_code.data(), record.lane_count));
      }
      if (composition_flux != nullptr) {
        composition_flux->validatePrepared(context.state);
        admitCompositionSparseActivation(
            context.state, *composition_flux, m_memory_governor,
            "hydro.distributed_amr.adaptive_chemistry_rows");
      }
    } catch (...) { failure = std::current_exception(); }
    m_mpi_context.rethrowCollectivePreparationFailure(failure, "distributed composition validation");
    if (composition_flux != nullptr) {
      m_last_composition_transport = composition_flux->commit(context.state);
      (void)physics::synchronizeMetalMirrorFromComposition(context.state, std::span<const std::uint32_t>{});
    }
  }

  // CHUI Part 3 Checkpoint B: automatic production refine/derefine (F6).
  // Runs after the AMR hydro sweep, its reflux and composition commit, i.e.
  // at a synchronized hierarchy boundary of the accepted interval; the
  // planner itself refuses (empty plan) when temporal histories or pending
  // flux registers exist. Every rank calls the collective id-window
  // reductions whether or not it has work, so MPI never deadlocks; a local
  // failure is turned into a collective error before any later stage.
  void runAutomaticAmrRegrid(
      core::StepContext& context,
      const amr::ProductionAmrHydroOptions& amr_options,
      std::span<const amr::PatchDescriptor> remote_neighbors) {
    if (m_config.numerics.amr_regrid_mode != core::AmrRegridConfigMode::kAutomatic) {
      return;
    }
    const std::uint64_t interval = static_cast<std::uint64_t>(m_config.numerics.amr_regrid_interval_steps);
    const std::uint64_t step = context.integrator_state.step_index + 1U;
    if (step % interval != 0U) {
      return;
    }
    const std::uint64_t epoch = step / interval;

    // All rank-local preparation that can allocate or throw is separated by a
    // collective failure gate before the next payload collective.  This is
    // deliberately more conservative than the old implementation, where
    // cooling evaluation / predictive-event preparation could throw on one
    // rank before its peers entered regrid planning.
    amr::ProductionAmrRegridPolicy policy;
    core::MemoryReservation cooling_time_reservation;      // outlives covered vector
    core::MemoryReservation predictive_local_reservation; // outlives local points
    core::MemoryReservation predictive_global_reservation;// outlives gathered/index scratch
    std::vector<double> cooling_time_code;
    std::vector<PredictiveRefinementPoint> local_predictive_points;
    std::vector<PredictiveRefinementPoint> global_predictive_points;
    std::vector<std::uint32_t> forced_rows;
    std::uint64_t local_forecast_events = 0U;
    std::exception_ptr failure;

    try {
      policy.criteria.use_mass_threshold = m_config.numerics.amr_refine_cell_mass_code > 0.0;
      policy.criteria.mass_threshold_code = m_config.numerics.amr_refine_cell_mass_code;
      policy.criteria.use_gradient_indicator = m_config.numerics.amr_refine_jump_threshold > 0.0;
      policy.criteria.gradient_threshold = m_config.numerics.amr_refine_jump_threshold;
      policy.criteria.use_jeans_condition = m_config.numerics.amr_refine_jeans_cells > 0.0;
      policy.criteria.jeans_resolution_cells = m_config.numerics.amr_refine_jeans_cells;
      policy.criteria.use_particle_count = false;
      policy.criteria.derefine_hysteresis = m_config.numerics.amr_derefine_hysteresis;
      policy.criteria.gravitational_constant_code = core::newtonGravitationalConstantCode(core::makeUnitSystem(
          m_config.units.length_unit, m_config.units.mass_unit, m_config.units.velocity_unit));

      policy.physics.cooling_length_cells = m_config.numerics.amr_refine_cooling_length_cells;
      policy.physics.metallicity_jump = m_config.numerics.amr_refine_metallicity_jump;
      policy.physics.species_jump = m_config.numerics.amr_refine_species_jump;
      policy.physics.vorticity = m_config.numerics.amr_refine_vorticity;
      policy.max_level = static_cast<std::uint8_t>(std::clamp(m_config.numerics.amr_max_level, 0, 255));
      policy.max_operations = memoryAwareOptionalRegridBudget(
          context.state, m_memory_governor,
          static_cast<std::size_t>(m_config.numerics.amr_regrid_max_operations));
      m_last_regrid_optional_budget = policy.max_operations;
      policy.pass = (epoch % 2U == 1U)
          ? amr::AmrRegridPassKind::kRefineOnly
          : amr::AmrRegridPassKind::kDerefineOnly;
      policy.scale_factor = m_mode_policy.cosmological_comoving_frame
          ? std::max(context.timeline_step.scale_factor_end, 1.0e-12)
          : 1.0;
      policy.periodic = amr_options.physical_boundary_kind == hydro::HydroBoundaryKind::kPeriodic;
      policy.periodic_box_extent_comoving = amr_options.periodic_box_extent_comoving;
      const double softening_kpc = m_config.numerics.gravity_softening_gas_kpc_comoving > 0.0
          ? m_config.numerics.gravity_softening_gas_kpc_comoving
          : m_config.numerics.gravity_softening_kpc_comoving;
      const core::UnitSystem units = core::makeUnitSystem(
          m_config.units.length_unit, m_config.units.mass_unit, m_config.units.velocity_unit);
      constexpr double k_kpc_m = 3.0856775814913673e19;
      policy.gravity_softening_code = softening_kpc > 0.0
          ? softening_kpc * k_kpc_m / units.length_si_per_code : 0.0;
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        failure, "automatic AMR regrid policy preparation");

    // Cooling-length refinement owns one dense-row cooling-time vector.  The
    // reservation is collective-before-allocation; every later evaluation
    // failure is reduced before predictive-AMR communication begins.
    if (policy.physics.cooling_length_cells > 0.0 && m_memory_governor != nullptr) {
      failure = nullptr;
      try {
        cooling_time_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            static_cast<std::uint64_t>(core::checkedSizeMultiply(
                context.state.cells.size(), sizeof(double), "regrid cooling times")),
            "hydro.amr_regrid.cooling_times");
        cooling_time_reservation.commit();
      } catch (...) {
        failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(
          failure, "regrid cooling-time admission");
    }
    if (policy.physics.cooling_length_cells > 0.0) {
      std::exception_ptr cooling_time_failure;
      try {
        const bool chemistry_history = context.state.gas_composition.enabled() &&
            context.state.gas_composition.chemistry_last_cooling_time_seconds.size() ==
                context.state.cells.size();
        // t_cool = e / (Lambda_net n_H^2 / rho) with physical density; net
        // heating (Lambda_net <= 0) imposes no cooling-length constraint (+inf).
        constexpr double k_proton_mass_g = 1.67262192369e-24;
        const double a = std::max(context.integrator_state.current_scale_factor, 1.0e-12);
        const bool comoving = m_config.units.coordinate_frame == core::CoordinateFrame::kComoving;
        const double redshift = m_mode_policy.cosmological_comoving_frame ? std::max(0.0, 1.0 / a - 1.0) : 0.0;
        cooling_time_code.assign(context.state.cells.size(), std::numeric_limits<double>::infinity());
        // A03: Part 2 thermochemistry without an accepted per-row history
        // (first interval, restart of a pre-history schema, regridded rows)
        // evaluates the SAME thermochemistry model's cooling time through the
        // source owner instead of refusing the mandatory criterion.
        const bool source_cooling_times = !chemistry_history &&
            !m_cooling_rate_provider.has_value() && m_source_runtime != nullptr &&
            m_source_runtime->fillAmrCoolingTimes(context, cooling_time_code);
        for (std::uint32_t row = 0;
             !source_cooling_times && row < context.state.cells.size(); ++row) {
          if (chemistry_history) {
            const double seconds =
                context.state.gas_composition.chemistry_last_cooling_time_seconds[row];
            if (seconds > 0.0 && std::isfinite(seconds)) {
              cooling_time_code[row] = seconds / m_units.timeSiPerCode();
            }
            continue;
          }
          if (!m_cooling_rate_provider.has_value()) {
            throw std::logic_error(
                "cooling-length refinement has neither accepted Part 2 cooling-time history "
                "nor a legacy cooling-rate provider");
          }
          const double e = context.state.gas_cells.internal_energy_code[row];
          const double rho_code = context.state.gas_cells.density_code[row] * (comoving ? 1.0 / (a * a * a) : 1.0);
          const double rho_cgs = m_units.densityCodeToCgs(rho_code);
          if (!(e > 0.0) || !(rho_cgs > 0.0)) continue;
          const double mass = context.state.cells.mass_code[row];
          const double z = mass > 0.0 ? std::clamp(context.state.gas_cells.metal_mass_code[row] / mass, 0.0, 1.0) : 0.0;
          const double n_h = (1.0 - m_cooling_helium_mass_fraction - z) * rho_cgs / k_proton_mass_g;
          const physics::CoolingRateResult rates = m_cooling_rate_provider->lookupRates(physics::CoolingRateQuery{
              .temperature_k = e * m_cooling_temperature_per_energy_k,
              .hydrogen_number_density_cgs = n_h,
              .metallicity_mass_fraction = z,
              .redshift = redshift});
          const double net = rates.netCoolingRateErgCm3S();
          if (!(net > 0.0)) continue;
          const double rate_code = net * n_h * n_h / rho_cgs * m_cooling_specific_rate_cgs_to_code;
          if (rate_code > 0.0) cooling_time_code[row] = e / rate_code;
        }
        policy.physics.cooling_time_code = cooling_time_code;
      } catch (...) {
        cooling_time_failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(
          cooling_time_failure, "regrid cooling-time evaluation");
    }
    policy.max_level = static_cast<std::uint8_t>(std::clamp(m_config.numerics.amr_max_level, 0, 255));
    policy.max_operations = static_cast<std::size_t>(m_config.numerics.amr_regrid_max_operations);
    // Alternating passes keep the 2:1 balance without exchanging decisions.
    policy.pass = (epoch % 2U == 1U) ? amr::AmrRegridPassKind::kRefineOnly : amr::AmrRegridPassKind::kDerefineOnly;
    policy.scale_factor = context.integrator_state.current_scale_factor > 0.0
        ? context.integrator_state.current_scale_factor : 1.0;
    policy.periodic = amr_options.physical_boundary_kind == hydro::HydroBoundaryKind::kPeriodic;
    policy.periodic_box_extent_comoving = amr_options.periodic_box_extent_comoving;
    {
      const double softening_kpc = m_config.numerics.gravity_softening_gas_kpc_comoving > 0.0
          ? m_config.numerics.gravity_softening_gas_kpc_comoving
          : m_config.numerics.gravity_softening_kpc_comoving;
      const core::UnitSystem units = core::makeUnitSystem(
          m_config.units.length_unit, m_config.units.mass_unit, m_config.units.velocity_unit);
      constexpr double k_kpc_m = 3.0856775814913673e19;
      policy.gravity_softening_code = softening_kpc > 0.0 ? softening_kpc * k_kpc_m / units.length_si_per_code : 0.0;
    }
    policy.force_resolution_action =
        m_config.numerics.amr_force_resolution_policy ==
                core::AmrForceResolutionPolicy::kAllowCountedReference
        ? amr::AmrForceResolutionAction::kAllowCountedReference
        : amr::AmrForceResolutionAction::kLimitOptionalRejectMandatory;

    core::V02ControlState v02_state;
    std::uint64_t v02_chemistry_reserve = 0U;
    std::uint64_t v02_optional_memory_budget = 0U;
    std::uint64_t regrid_baseline_bytes = 0U;
    std::uint64_t regrid_estimated_bytes_per_cell = 0U;
    const bool joint_v02 =
        m_config.numerics.v02_control_mode == core::V02ControlMode::kJointErrorResourceV1;
    std::exception_ptr v02_preparation_failure;
    try {
      regrid_baseline_bytes = core::memoryReportBaselineOwnedBytes(
          core::collectSimulationMemoryReport(context.state));
      regrid_estimated_bytes_per_cell =
          context.state.cells.size() == 0U
          ? 0U
          : regrid_baseline_bytes / context.state.cells.size() +
                (regrid_baseline_bytes % context.state.cells.size() != 0U
                     ? 1U
                     : 0U);
      if (joint_v02) {
        v02_state = core::loadV02ControlState(context.state);
        const auto& composition = context.state.gas_composition;
        const std::uint64_t inactive_rows = composition.sparseStorageEnabled()
            ? static_cast<std::uint64_t>(composition.rowCount() - composition.sparseActiveRowCount())
            : 0U;
        const std::uint64_t sparse_row_bytes = sizeof(std::uint64_t) +
            static_cast<std::uint64_t>(composition.sparseLaneCount()) * sizeof(double);
        const std::uint64_t ranks = static_cast<std::uint64_t>(
            std::max(1, m_mpi_context.worldSize()));
        const std::uint64_t observed_growth_rows = std::min(
            inactive_rows,
            (v02_state.chemistry_promoted_last_interval + ranks - 1U) / ranks);
        const std::uint64_t observed_growth_bytes =
            observed_growth_rows > 0U && sparse_row_bytes >
                    std::numeric_limits<std::uint64_t>::max() / observed_growth_rows
            ? std::numeric_limits<std::uint64_t>::max()
            : observed_growth_rows * sparse_row_bytes;
        const std::uint64_t maximum_growth_bytes =
            inactive_rows > 0U && sparse_row_bytes >
                    std::numeric_limits<std::uint64_t>::max() / inactive_rows
            ? std::numeric_limits<std::uint64_t>::max()
            : inactive_rows * sparse_row_bytes;
        v02_chemistry_reserve = std::min(
            maximum_growth_bytes,
            std::max(v02_state.chemistry_growth_reserve_bytes -
                         v02_state.chemistry_growth_reserve_bytes / 10U,
                     observed_growth_bytes));
        const core::MemoryGovernorSnapshot memory = m_memory_governor != nullptr
            ? m_memory_governor->snapshot()
            : core::MemoryGovernorSnapshot{
                  .headroom_bytes = std::numeric_limits<std::uint64_t>::max()};
        const std::uint64_t fixed_reserve = regrid_baseline_bytes >
                std::numeric_limits<std::uint64_t>::max() - v02_chemistry_reserve
            ? std::numeric_limits<std::uint64_t>::max()
            : regrid_baseline_bytes + v02_chemistry_reserve;
        const std::uint64_t available = memory.headroom_bytes ==
                std::numeric_limits<std::uint64_t>::max()
            ? memory.headroom_bytes
            : memory.headroom_bytes > fixed_reserve
                ? memory.headroom_bytes - fixed_reserve
                : 0U;
        if (available == std::numeric_limits<std::uint64_t>::max()) {
          v02_optional_memory_budget = available;
        } else {
          v02_optional_memory_budget = static_cast<std::uint64_t>(
              static_cast<long double>(available) *
              m_config.numerics.v02_optional_memory_headroom_fraction);
        }
        policy.optional_control = amr::OptionalRefinementControl{
            .enabled = true,
            .benefit_floor = std::max(
                1.0, 1.0 + 0.9 * (v02_state.optional_benefit_floor - 1.0)),
            .memory_budget_bytes = v02_optional_memory_budget,
            .estimated_bytes_per_cell = regrid_estimated_bytes_per_cell,
            .cell_step_budget = std::numeric_limits<std::uint64_t>::max(),
            .mpi_face_cell_budget = std::numeric_limits<std::uint64_t>::max(),
            .churn_cell_budget = std::numeric_limits<std::uint64_t>::max()};
      }
      if (!joint_v02) {
        policy.optional_control.estimated_bytes_per_cell =
            regrid_estimated_bytes_per_cell;
      }
    } catch (...) {
      v02_preparation_failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        v02_preparation_failure, "V02 AMR controller preparation");

    // A04 predictive event-aware refinement.  Local event chronology is
    // evaluated first; the resulting current/drifted positions are then
    // all-gathered so every destination rank can force its own host cell.
    // This removes the previous decomposition-dependent "drifted host on
    // another rank is not forced" limitation.
    const double predictive_age = m_config.numerics.amr_predictive_feedback_age_years;
    const double lookahead_years = m_config.numerics.amr_predictive_event_lookahead_years;
    if (predictive_age > 0.0 || lookahead_years > 0.0) {
      std::vector<amr::PatchDescriptor> descriptors;
      std::size_t maximum_forecast_positions = 0U;
      core::MemoryReservation forecast_positions_reservation;
      std::exception_ptr forecast_positions_admission_failure;
      std::vector<std::array<double, 3>> forecast_positions;
      try {
        descriptors = amr::buildProductionAmrPatchDescriptors(context.state);
        maximum_forecast_positions = core::checkedSizeMultiply(
            context.state.star_particles.size(),
            lookahead_years > 0.0 ? 2U : 1U,
            "predictive refinement local position capacity");
        if (m_memory_governor != nullptr && maximum_forecast_positions > 0U) {
          forecast_positions_reservation = m_memory_governor->reserve(
              core::MemoryClass::kPhaseResident,
              static_cast<std::uint64_t>(core::checkedSizeMultiply(
                  maximum_forecast_positions, sizeof(std::array<double, 3>),
                  "predictive refinement local positions")),
              "hydro.amr_regrid.local_predictive_positions");
          forecast_positions_reservation.commit();
        }
        forecast_positions.reserve(maximum_forecast_positions);
        if (forecast_positions.capacity() > maximum_forecast_positions) {
          throw std::length_error(
              "predictive refinement position allocator exceeded admitted capacity");
        }
      } catch (...) {
        forecast_positions_admission_failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(
          forecast_positions_admission_failure,
          "predictive refinement local position admission");
      const auto host_row = [&](const std::array<double, 3>& position) {
        for (const amr::PatchDescriptor& patch : descriptors) {
          bool inside = true;
          for (std::size_t c = 0; c < 3U && inside; ++c) {
            inside = position[c] >= patch.origin_comov[c] && position[c] < patch.origin_comov[c] + patch.extent_comov[c];
          }
          if (!inside) continue;
          const auto cell = amr::patchLocalCellForPoint(patch, position, "predictive refinement host cell");
          const auto rows = amr::rowsForPatchLocalCellOrder(context.state, patch, "predictive refinement");
          if (cell.linear_index < rows.size()) forced_rows.push_back(rows[cell.linear_index]);
          return;
        }
      };
      // R2 (A04) event-aware forecast. For every star the canonical
      // population record (restart/rollback state) gives its deterministic
      // event chronology; a COPY is advanced to the expected cumulative count
      // at age + lookahead, so the forecast consumes no event ordinal and
      // changes no budget. Expectation mode (no stochastic events) predicts
      // an "event" when >= 0.5 expected SNe fall in the window (HEURISTIC).
      // The host cell now and at the drifted position x + v dt / a are forced.
      std::uint64_t forecast_events = 0U;
      std::exception_ptr forecast_failure;
      try {
        std::optional<physics::SspModel> ssp;
        std::optional<physics::StellarPopulationStore> store;
        if (lookahead_years > 0.0) {
          physics::SspConfig ssp_config;
          ssp_config.imf_mass_min_msun = m_config.physics.ssp_imf_mass_min_msun;
          ssp_config.imf_mass_max_msun = m_config.physics.ssp_imf_mass_max_msun;
          ssp_config.cc_mass_min_msun = m_config.physics.ssp_cc_mass_min_msun;
          ssp_config.cc_mass_max_msun = m_config.physics.ssp_cc_mass_max_msun;
          ssp_config.snia_per_msun = m_config.physics.ssp_snia_per_msun;
          ssp_config.snia_delay_min_years = m_config.physics.ssp_snia_delay_min_years;
          ssp_config.snia_delay_max_years = m_config.physics.ssp_snia_delay_max_years;
          ssp_config.stochastic_events = m_config.physics.ssp_stochastic_sn_events;
          ssp.emplace(ssp_config);  // event counts only (lifetimes + IMF); no yields used
          store.emplace(context.state);
        }
        constexpr double k_solar_mass_kg = 1.98847e30;
        constexpr double k_year_s = 3.15576e7;
        const double msun_code = k_solar_mass_kg / m_units.mass_si_per_code;
        const double a = std::max(context.integrator_state.current_scale_factor, 1.0e-12);
        const double drift_scale = m_config.units.coordinate_frame == core::CoordinateFrame::kComoving
            ? 1.0 / a : 1.0;
        const double lookahead_code = lookahead_years * k_year_s / m_units.timeSiPerCode();
        for (std::size_t s = 0; s < context.state.star_particles.size(); ++s) {
          const double age = context.state.star_particles.stellar_age_years_last[s];
          const std::uint32_t p = context.state.star_particles.particle_index[s];
          if (p >= context.state.particles.size()) {
            throw std::runtime_error("predictive AMR stellar particle index is outside particle storage");
          }
          bool forced = predictive_age > 0.0 && age < predictive_age;
          if (!forced && ssp.has_value()) {
            const std::uint64_t pid = context.state.particle_sidecar.particle_id[p];
            if (store->has(pid)) {
              const physics::StellarPopulationRecord record = store->get(pid);
              forced = physics::forecastStellarEvent(
                  *ssp, record, age, lookahead_years,
                  context.state.star_particles.birth_mass_code[s] / msun_code,
                  m_config.physics.ssp_stochastic_sn_events);
              if (forced) ++local_forecast_events;
            }
          }
          if (!forced) continue;
          PredictiveRefinementPoint position{
              .x_comoving = context.state.particles.position_x_comoving[p],
              .y_comoving = context.state.particles.position_y_comoving[p],
              .z_comoving = context.state.particles.position_z_comoving[p]};
          local_predictive_points.push_back(position);
          if (lookahead_code > 0.0) {
            PredictiveRefinementPoint drifted{
                .x_comoving = position.x_comoving +
                    context.state.particles.velocity_x_peculiar[p] * lookahead_code * drift_scale,
                .y_comoving = position.y_comoving +
                    context.state.particles.velocity_y_peculiar[p] * lookahead_code * drift_scale,
                .z_comoving = position.z_comoving +
                    context.state.particles.velocity_z_peculiar[p] * lookahead_code * drift_scale};
            if (amr_options.physical_boundary_kind == hydro::HydroBoundaryKind::kPeriodic) {
              double* coordinates[3] = {&drifted.x_comoving, &drifted.y_comoving, &drifted.z_comoving};
              for (std::size_t c = 0; c < 3U; ++c) {
                const double box = amr_options.periodic_box_extent_comoving[c];
                if (box > 0.0) {
                  *coordinates[c] -= box * std::floor(*coordinates[c] / box);
                }
              }
            }
            local_predictive_points.push_back(drifted);
          }
        }
      }
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        failure, "predictive AMR local forecast preparation");

    if (predictive_enabled) {
      std::uint64_t local_point_count = 0U;
      failure = nullptr;
      try {
        local_point_count = core::checkedIntegralNarrow<std::uint64_t>(
            local_predictive_points.size(), "predictive AMR local point count");
      } catch (...) {
        failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(
          failure, "predictive AMR local point-count preparation");
      const std::uint64_t global_point_count =
          m_mpi_context.allreduceSumUint64(local_point_count);
      if (m_memory_governor != nullptr) {
        failure = nullptr;
        try {
          const std::size_t global_points = core::checkedIntegralNarrow<std::size_t>(
              global_point_count, "predictive AMR global point count");
          const std::size_t point_bytes = core::checkedSizeMultiply(
              global_points, sizeof(PredictiveRefinementPoint),
              "predictive AMR global point bytes");
          // allgatherPredictiveRefinementPoints briefly owns the gathered byte
          // buffer and the typed point vector at the same time.
          const std::size_t point_exchange_bytes = core::checkedSizeMultiply(
              point_bytes, 2U, "predictive AMR point exchange coexistence bytes");
          const std::size_t patch_count = context.state.patches.size();
          const std::size_t descriptor_bytes = core::checkedSizeMultiply(
              patch_count, sizeof(amr::PatchDescriptor),
              "predictive AMR descriptor bytes");
          struct PredictiveLookupBudgetEntry {
            std::size_t descriptor_index;
            double lower_x;
            double upper_x;
          };
          const std::size_t lookup_bytes = core::checkedSizeMultiply(
              patch_count, sizeof(PredictiveLookupBudgetEntry),
              "predictive AMR lookup bytes");
          const std::size_t prefix_bytes = core::checkedSizeMultiply(
              patch_count, sizeof(double), "predictive AMR prefix bytes");
          const std::size_t patch_id_index_bytes = core::checkedSizeMultiply(
              patch_count, sizeof(std::pair<std::uint64_t, std::size_t>),
              "predictive AMR patch-id index bytes");
          const std::size_t offset_bytes = core::checkedSizeMultiply(
              core::checkedSizeAdd(patch_count, 1U, "predictive AMR offset count"),
              sizeof(std::size_t), "predictive AMR offset bytes");
          const std::size_t row_cache_bytes = core::checkedSizeMultiply(
              context.state.cells.size(), sizeof(std::uint32_t),
              "predictive AMR row-cache bytes");
          const std::size_t forced_bytes = core::checkedSizeMultiply(
              std::min(global_points, context.state.cells.size()),
              sizeof(std::uint32_t), "predictive AMR forced-row bytes");
          std::size_t total = point_exchange_bytes;
          for (const std::size_t bytes : {descriptor_bytes, lookup_bytes, prefix_bytes,
                                         patch_id_index_bytes, offset_bytes,
                                         row_cache_bytes, forced_bytes}) {
            total = core::checkedSizeAdd(total, bytes,
                                         "predictive AMR global scratch bytes");
          }
          // Bounded allowance for vector allocator bookkeeping and the small
          // MPI count/displacement arrays used by the bounded all-gather.
          total = core::checkedSizeAdd(
              total,
              core::checkedSizeMultiply(
                  core::checkedSizeAdd(patch_count,
                      static_cast<std::size_t>(std::max(1, m_mpi_context.worldSize())),
                      "predictive AMR overhead item count"),
                  64U, "predictive AMR allocator overhead bytes"),
              "predictive AMR scratch plus allocator overhead");
          const std::uint64_t total_u64 = core::checkedIntegralNarrow<std::uint64_t>(
              total, "predictive AMR global scratch reservation");
          if (total_u64 > 0U) {
            predictive_global_reservation = m_memory_governor->reserve(
                core::MemoryClass::kScratchArena, total_u64,
                "hydro.amr_regrid.predictive_global_lookup");
            predictive_global_reservation.commit();
          }
        } catch (...) {
          failure = std::current_exception();
        }
        m_mpi_context.rethrowCollectivePreparationFailure(
            failure, "predictive AMR global lookup admission");
      }

      global_predictive_points = allgatherPredictiveRefinementPoints(
          m_mpi_context, local_predictive_points);
      failure = nullptr;
      try {
        const std::vector<amr::PatchDescriptor> descriptors =
            amr::buildProductionAmrPatchDescriptors(context.state);
        struct PatchLookupEntry {
          std::size_t descriptor_index = 0U;
          double lower_x = 0.0;
          double upper_x = 0.0;
        };
        std::vector<PatchLookupEntry> lookup;
        lookup.reserve(descriptors.size());
        for (std::size_t i = 0; i < descriptors.size(); ++i) {
          lookup.push_back(PatchLookupEntry{
              .descriptor_index = i,
              .lower_x = descriptors[i].origin_comov[0],
              .upper_x = descriptors[i].origin_comov[0] + descriptors[i].extent_comov[0]});
        }
        std::sort(lookup.begin(), lookup.end(), [&](const PatchLookupEntry& lhs, const PatchLookupEntry& rhs) {
          if (lhs.lower_x != rhs.lower_x) return lhs.lower_x < rhs.lower_x;
          return descriptors[lhs.descriptor_index].patch_id < descriptors[rhs.descriptor_index].patch_id;
        });
        std::vector<double> prefix_max_upper;
        prefix_max_upper.reserve(lookup.size());
        double max_upper = -std::numeric_limits<double>::infinity();
        for (const PatchLookupEntry& entry : lookup) {
          max_upper = std::max(max_upper, entry.upper_x);
          prefix_max_upper.push_back(max_upper);
        }
        // Flat patch-local row cache: no unordered_map/node allocations and
        // therefore an exact O(Npatch + Ncell) admitted workspace.
        std::vector<std::pair<std::uint64_t, std::size_t>> descriptor_by_patch_id;
        descriptor_by_patch_id.reserve(descriptors.size());
        std::vector<std::size_t> row_offsets(descriptors.size() + 1U, 0U);
        for (std::size_t i = 0; i < descriptors.size(); ++i) {
          descriptor_by_patch_id.emplace_back(descriptors[i].patch_id, i);
          row_offsets[i + 1U] = core::checkedSizeAdd(
              row_offsets[i], amr::checkedPatchCellCount(
                  descriptors[i], "predictive refinement row-cache size"),
              "predictive refinement row-cache prefix");
        }
        if (row_offsets.back() != context.state.cells.size()) {
          throw std::runtime_error(
              "predictive AMR patch descriptors do not cover the local gas-cell extent exactly");
        }
        std::sort(descriptor_by_patch_id.begin(), descriptor_by_patch_id.end());
        std::vector<std::uint32_t> row_cache(
            row_offsets.back(), core::kInvalidGasCellRow);
        for (std::uint32_t row = 0U; row < context.state.cells.size(); ++row) {
          const auto* identity = context.state.gas_cell_identity.findByLocalRow(row);
          if (identity == nullptr) {
            throw std::runtime_error(
                "predictive AMR row-cache construction found a gas row without identity");
          }
          const auto by_id = std::lower_bound(
              descriptor_by_patch_id.begin(), descriptor_by_patch_id.end(),
              std::pair<std::uint64_t, std::size_t>{identity->owning_patch_id, 0U},
              [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
          if (by_id == descriptor_by_patch_id.end() ||
              by_id->first != identity->owning_patch_id) {
            throw std::runtime_error(
                "predictive AMR row-cache construction found an unknown owning patch");
          }
          const std::size_t descriptor_index = by_id->second;
          const auto cell = amr::patchLocalCellForRow(
              context.state, descriptors[descriptor_index], row,
              "predictive refinement row-cache mapping");
          const std::size_t flat = core::checkedSizeAdd(
              row_offsets[descriptor_index], cell.linear_index,
              "predictive refinement flat row index");
          if (flat >= row_cache.size() || row_cache[flat] != core::kInvalidGasCellRow) {
            throw std::runtime_error(
                "predictive AMR row-cache construction found a duplicate local cell");
          }
          row_cache[flat] = row;
        }
        if (std::find(row_cache.begin(), row_cache.end(), core::kInvalidGasCellRow) !=
            row_cache.end()) {
          throw std::runtime_error(
              "predictive AMR row-cache construction found incomplete patch coverage");
        }
        forced_rows.reserve(std::min<std::size_t>(
            global_predictive_points.size(), context.state.cells.size()));

        const auto host_row = [&](const PredictiveRefinementPoint& point) {
          const auto upper = std::upper_bound(
              lookup.begin(), lookup.end(), point.x_comoving,
              [](double x, const PatchLookupEntry& entry) { return x < entry.lower_x; });
          std::size_t pos = static_cast<std::size_t>(std::distance(lookup.begin(), upper));
          while (pos > 0U) {
            if (prefix_max_upper[pos - 1U] <= point.x_comoving) break;
            const PatchLookupEntry& entry = lookup[pos - 1U];
            const amr::PatchDescriptor& patch = descriptors[entry.descriptor_index];
            const bool inside_x = point.x_comoving >= patch.origin_comov[0] &&
                point.x_comoving < patch.origin_comov[0] + patch.extent_comov[0];
            const bool inside_y = point.y_comoving >= patch.origin_comov[1] &&
                point.y_comoving < patch.origin_comov[1] + patch.extent_comov[1];
            const bool inside_z = point.z_comoving >= patch.origin_comov[2] &&
                point.z_comoving < patch.origin_comov[2] + patch.extent_comov[2];
            if (inside_x && inside_y && inside_z) {
              const auto cell = amr::patchLocalCellForPoint(
                  patch,
                  std::array<double, 3>{point.x_comoving, point.y_comoving, point.z_comoving},
                  "predictive refinement host cell");
              const std::size_t flat = core::checkedSizeAdd(
                  row_offsets[entry.descriptor_index], cell.linear_index,
                  "predictive refinement host row");
              if (flat >= row_cache.size() || row_cache[flat] == core::kInvalidGasCellRow) {
                throw std::runtime_error("predictive AMR host cell exceeds cached patch rows");
              }
              forced_rows.push_back(row_cache[flat]);
              return;
            }
            --pos;
          }
        };
        for (const PredictiveRefinementPoint& point : global_predictive_points) {
          host_row(point);
        }
        std::sort(forced_rows.begin(), forced_rows.end());
        forced_rows.erase(std::unique(forced_rows.begin(), forced_rows.end()), forced_rows.end());
        policy.forced_refine_rows = forced_rows;
      } catch (...) {
        failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(
          failure, "predictive AMR distributed host resolution");
      m_last_predictive_forecast_events = m_mpi_context.allreduceSumUint64(local_forecast_events);
    } else {
      m_last_predictive_forecast_events = 0U;
    }

    amr::ProductionAmrRegridPlan plan;
    failure = nullptr;
    try {
      plan = amr::planProductionAmrRegrid(context.state, remote_neighbors, policy);
      if (plan.mandatory_truncated_by_budget != 0U) {
        throw std::logic_error(
            "automatic AMR planner truncated a mandatory refinement request; hard criteria must bypass the optional budget");
      }
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        failure, "automatic AMR regrid planning");

    // A hard criterion that has reached max_level cannot be repaired by a
    // later balancing pass.  Fail collectively instead of silently continuing
    // below the configured scientific resolution requirement.
    std::uint64_t local_mandatory_max_level = 0U;
    std::uint64_t local_mandatory_balance = 0U;
    failure = nullptr;
    try {
      local_mandatory_max_level = core::checkedIntegralNarrow<std::uint64_t>(
          plan.mandatory_blocked_by_max_level,
          "mandatory AMR max-level block count");
      local_mandatory_balance = core::checkedIntegralNarrow<std::uint64_t>(
          plan.mandatory_blocked_by_balance,
          "mandatory AMR balance block count");
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        failure, "mandatory AMR block-count preparation");
    const std::uint64_t global_mandatory_max_level =
        m_mpi_context.allreduceSumUint64(local_mandatory_max_level);
    const std::uint64_t global_mandatory_balance =
        m_mpi_context.allreduceSumUint64(local_mandatory_balance);
    if (global_mandatory_max_level != 0U) {
      throw std::runtime_error(
          "automatic AMR cannot satisfy a mandatory refinement criterion because amr_max_level was reached");
    }
    if (global_mandatory_balance != 0U) {
      throw std::runtime_error(
          "automatic AMR cannot satisfy a mandatory refinement criterion without first restoring 2:1 balance; "
          "refusing the step rather than advancing an under-resolved hierarchy");
    }

    std::uint64_t local_max_patch = 0U;
    for (const std::uint64_t id : context.state.patches.patch_id) {
      local_max_patch = std::max(local_max_patch, id);
    }
    std::uint64_t local_max_cell = 0U;
    for (const std::uint64_t id : context.state.gas_cells.gas_cell_id) {
      local_max_cell = std::max(local_max_cell, id);
    }
    const std::uint64_t global_max_patch = m_mpi_context.allreduceMaxUint64(local_max_patch);
    const std::uint64_t global_max_cell = m_mpi_context.allreduceMaxUint64(local_max_cell);
    const std::uint64_t patch_window = m_mpi_context.allreduceMaxUint64(plan.required_patch_ids);
    const std::uint64_t cell_window = m_mpi_context.allreduceMaxUint64(plan.required_gas_cell_ids);
    if (patch_window == 0U && cell_window == 0U) {
      // Persist rejected/no-op decisions as well: otherwise the controller's
      // hysteresis and resource history would only advance when topology did.
      core::MemoryReservation v02_state_reservation;
      failure = nullptr;
      try {
        if (joint_v02 && m_memory_governor != nullptr) {
          v02_state_reservation = m_memory_governor->reserve(
              core::MemoryClass::kPhaseResident, 4096U,
              "hydro.amr_regrid.v02_controller_state");
          v02_state_reservation.commit();
        }
        persist_v02_regrid_state(context.state);
      } catch (...) {
        failure = std::current_exception();
      }
      m_mpi_context.rethrowCollectivePreparationFailure(
          failure, "automatic AMR V02 no-op state publication");
      return;
    }
    amr::ProductionAmrRegridIdWindow window;
    failure = nullptr;
    try {
      window = amr::rankStripedAmrRegridIdWindow(
          global_max_patch, global_max_cell, m_mpi_context.worldRank(),
          std::max(1, m_mpi_context.worldSize()), patch_window, cell_window);
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        failure, "automatic AMR identity-window preparation");
    amr::ProductionAmrRegridSummary summary;
    failure = nullptr;
    try {
      summary = amr::applyProductionAmrRegridPlan(
          context.state, plan, window, amr_options);
    } catch (...) {
      failure = std::current_exception();
    }
    m_mpi_context.rethrowCollectivePreparationFailure(
        failure, "automatic AMR regrid application");
    m_last_regrid_summary = summary;
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.amr_automatic_regrid",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.amr",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "automatic conservative AMR regrid at synchronized hierarchy boundary",
          .payload = {{"pass", policy.pass == amr::AmrRegridPassKind::kRefineOnly ? "refine" : "derefine"},
                      {"refinements", std::to_string(summary.applied_refinements)},
                      {"derefinements", std::to_string(summary.applied_derefinements)},
                      {"created_gas_cells", std::to_string(summary.created_gas_cells)},
                      {"retired_gas_cells", std::to_string(summary.retired_gas_cells)},
                      {"blocked_by_balance", std::to_string(plan.blocked_by_balance)},
                      {"blocked_by_max_level", std::to_string(plan.blocked_by_max_level)},
                      {"skipped_remote_sibling", std::to_string(plan.skipped_remote_sibling)},
                      {"optional_operation_budget", std::to_string(m_last_regrid_optional_budget)},
                      {"truncated_by_budget", std::to_string(plan.truncated_by_budget)},
                      {"mandatory_refine_requests", std::to_string(plan.mandatory_refine_requests)},
                      {"predictive_forecast_event_stars", std::to_string(m_last_predictive_forecast_events)},
                      {"optional_refine_requests", std::to_string(plan.optional_refine_requests)},
                      {"mandatory_blocked_by_balance", std::to_string(plan.mandatory_blocked_by_balance)},
                      {"mandatory_blocked_by_max_level", std::to_string(plan.mandatory_blocked_by_max_level)},
                      {"mandatory_truncated_by_budget", std::to_string(plan.mandatory_truncated_by_budget)},
                      {"optional_truncated_by_budget", std::to_string(plan.optional_truncated_by_budget)},
                      {"optional_blocked_by_force_resolution",
                       std::to_string(plan.optional_blocked_by_force_resolution)},
                      {"optional_admitted_cell_steps",
                       std::to_string(plan.optional_admitted_cell_steps)},
                      {"optional_admitted_mpi_face_cells",
                       std::to_string(plan.optional_admitted_mpi_face_cells)},
                      {"optional_admitted_churn_cells",
                       std::to_string(plan.optional_admitted_churn_cells)},
                      {"optional_admitted_memory_bytes",
                       std::to_string(plan.optional_admitted_memory_bytes)},
                      {"v02_chemistry_growth_reserve_bytes",
                       std::to_string(v02_chemistry_reserve)},
                      {"v02_optional_memory_budget_bytes",
                       std::to_string(v02_optional_memory_budget)},
                      {"refined_below_gravity_softening_unqualified",
                       std::to_string(plan.refined_below_gravity_softening)},
                      {"skipped_not_quiescent", plan.skipped_not_quiescent ? "true" : "false"},
                      {"mass_change_code", std::to_string(summary.mass_change_code)},
                      {"total_energy_change_code", std::to_string(summary.total_energy_change_code)}}});
    }
  }

  struct HydroGeometryCacheKey {
    std::size_t cell_count = 0;
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;
    double origin_x_comoving = 0.0;
    double origin_y_comoving = 0.0;
    double origin_z_comoving = 0.0;
    double cell_width_x_comoving = 0.0;
    double cell_width_y_comoving = 0.0;
    double cell_width_z_comoving = 0.0;
    double dt_time_code = 0.0;
    core::BoundaryCondition hydro_boundary = core::BoundaryCondition::kPeriodic;
    std::uint64_t patch_signature = 0;
    std::uint64_t gas_cell_identity_generation = 0;
    std::uint64_t cell_index_generation = 0;
    std::uint64_t row_mapping_signature = 0;

    [[nodiscard]] bool operator==(const HydroGeometryCacheKey& rhs) const noexcept {
      return cell_count == rhs.cell_count &&
          nx == rhs.nx &&
          ny == rhs.ny &&
          nz == rhs.nz &&
          origin_x_comoving == rhs.origin_x_comoving &&
          origin_y_comoving == rhs.origin_y_comoving &&
          origin_z_comoving == rhs.origin_z_comoving &&
          cell_width_x_comoving == rhs.cell_width_x_comoving &&
          cell_width_y_comoving == rhs.cell_width_y_comoving &&
          cell_width_z_comoving == rhs.cell_width_z_comoving &&
          dt_time_code == rhs.dt_time_code &&
          hydro_boundary == rhs.hydro_boundary &&
          patch_signature == rhs.patch_signature &&
          gas_cell_identity_generation == rhs.gas_cell_identity_generation &&
          cell_index_generation == rhs.cell_index_generation &&
          row_mapping_signature == rhs.row_mapping_signature;
    }
  };

  [[nodiscard]] static std::uint64_t patchSignature(const core::PatchSoa& patches) {
    std::uint64_t signature = 1469598103934665603ULL;
    const auto mix = [&signature](std::uint64_t value) {
      signature ^= value;
      signature *= 1099511628211ULL;
    };
    mix(static_cast<std::uint64_t>(patches.size()));
    for (std::size_t patch = 0; patch < patches.size(); ++patch) {
      mix(patches.patch_id[patch]);
      mix(static_cast<std::uint64_t>(static_cast<std::int64_t>(patches.level[patch])));
      mix(patches.first_cell[patch]);
      mix(patches.cell_count[patch]);
      mix(patches.owning_rank[patch]);
    }
    return signature;
  }

  void rebuildGeometryIfNeeded(const core::SimulationState& state, double dt_time_code) {
    const std::size_t cell_count = state.cells.size();
    if (cell_count == 0) {
      m_cached_geometry_key.reset();
      m_geometry = {};
      return;
    }

    const internal::CartesianGasCellRowLayout layout = requireCartesianGasCellRowLayout(
        state, m_config, "reference hydro workflow");
    const hydro::HydroCartesianPatchSpec& spec = layout.spec;
    const HydroGeometryCacheKey key{
        .cell_count = cell_count,
        .nx = spec.nx,
        .ny = spec.ny,
        .nz = spec.nz,
        .origin_x_comoving = spec.origin_x_comoving,
        .origin_y_comoving = spec.origin_y_comoving,
        .origin_z_comoving = spec.origin_z_comoving,
        .cell_width_x_comoving = spec.cell_width_x_comoving,
        .cell_width_y_comoving = spec.cell_width_y_comoving,
        .cell_width_z_comoving = spec.cell_width_z_comoving,
        .dt_time_code = dt_time_code,
        .hydro_boundary = m_mode_policy.hydro_boundary,
        .patch_signature = patchSignature(state.patches),
        .gas_cell_identity_generation = state.gasCellIdentityGeneration(),
        .cell_index_generation = state.cellIndexGeneration(),
        .row_mapping_signature = layout.mapping_signature,
    };
    if (m_cached_geometry_key.has_value() && *m_cached_geometry_key == key) {
      return;
    }

    m_cached_geometry_key = key;
    m_dense_row_by_geometry_row = layout.dense_row_by_geometry_row;
    m_geometry_row_by_dense_row = layout.geometry_row_by_dense_row;
    m_geometry = hydro::makeCartesianPatchGeometry(spec);
    hydro::appendCartesianBoundaryGhostFaces(
        m_geometry,
        hydroBoundaryKindFromModePolicy(m_mode_policy.hydro_boundary));
    // Complete along-normal stencils through the boundary ghost shell:
    // periodic boxes wrap to real cells (full order everywhere), open and
    // reflective boundaries use the ghost row as the stencil cell.
    (void)hydro::completeReconstructionStencils(m_geometry);
    // The reconstruction (and its dt/(a dx) Courant factors) is rebuilt per
    // step once the flux-update epoch is known; see the hydro callback.
  }

  [[nodiscard]] hydro::HydroActiveSetView buildActiveFaceView(std::span<const std::uint32_t> active_cells) {
    m_active_cells.clear();
    m_active_cells.reserve(active_cells.size());
    for (const std::uint32_t dense_row : active_cells) {
      if (dense_row >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("active gas-cell row is outside Cartesian row map");
      }
      m_active_cells.push_back(m_geometry_row_by_dense_row[dense_row]);
    }
    std::unordered_set<std::size_t> active_cell_lookup(m_active_cells.begin(), m_active_cells.end());
    m_active_faces.clear();
    for (std::size_t face_index = 0; face_index < m_geometry.faces.size(); ++face_index) {
      const hydro::HydroFace& face = m_geometry.faces[face_index];
      if (m_remote_non_authority_faces.contains(face_index)) {
        continue;
      }
      std::size_t neighbor_for_activity = face.neighbor_cell;
      if (face.ghost_cell_slot != hydro::k_invalid_ghost_cell_slot) {
        const hydro::HydroGhostCell& ghost = m_geometry.ghost_cells.at(face.ghost_cell_slot);
        if (ghost.boundary_kind == hydro::HydroBoundaryKind::kImportedMpi &&
            m_current_world_rank > static_cast<std::uint32_t>(ghost.origin_rank)) {
          continue;
        }
        if (ghost.boundary_kind == hydro::HydroBoundaryKind::kPeriodic) {
          // A periodic box face appears twice (lower-side ghost face of cell
          // 0 and upper-side ghost face of cell n-1); both deposit on the
          // same two real cells. Evaluate only the upper-side copy so the
          // boundary flux is applied exactly once, and treat its wrapped
          // source cell as the neighbor for activity.
          if (ghost.side == hydro::HydroFaceSide::kLower) {
            continue;
          }
          neighbor_for_activity = ghost.source_real_cell;
        }
      }
      if (active_cell_lookup.contains(face.owner_cell) ||
          (neighbor_for_activity != hydro::k_invalid_cell_index &&
           active_cell_lookup.contains(neighbor_for_activity))) {
        m_active_faces.push_back(face_index);
      }
    }
    return hydro::HydroActiveSetView{.active_cells = m_active_cells, .active_faces = m_active_faces};
  }

  [[nodiscard]] std::optional<std::pair<std::uint64_t, std::uint32_t>> patchIdentityForCellRow(
      const core::SimulationState& state,
      std::uint32_t cell_index) const {
    const auto* identity = state.gas_cell_identity.findByLocalRow(cell_index);
    if (identity == nullptr || identity->gas_cell_id == 0U) {
      throw std::runtime_error("hydro CFL diagnostics rejected missing gas-cell identity");
    }
    if (state.patches.size() == 0U) {
      return std::nullopt;
    }
    if (identity->owning_patch_id == 0U || cell_index >= m_geometry_row_by_dense_row.size()) {
      throw std::runtime_error("hydro CFL diagnostics rejected incomplete fixed-patch identity metadata");
    }
    const auto patch_it = std::find(
        state.patches.patch_id.begin(), state.patches.patch_id.end(), identity->owning_patch_id);
    if (patch_it == state.patches.patch_id.end()) {
      throw std::runtime_error("hydro CFL diagnostics rejected an identity patch absent from PatchSoa");
    }
    const std::size_t patch_index = static_cast<std::size_t>(std::distance(state.patches.patch_id.begin(), patch_it));
    if (state.cells.patch_index[cell_index] != patch_index) {
      throw std::runtime_error("hydro CFL diagnostics rejected stale dense patch-index metadata");
    }
    return std::pair<std::uint64_t, std::uint32_t>{
        identity->owning_patch_id,
        m_geometry_row_by_dense_row[cell_index]};
  }

  [[nodiscard]] core::HydroCflDiagnostics hydroCflDiagnosticsForCell(
      const core::SimulationState& state,
      std::uint32_t cell_index,
      double accepted_dt_time_code,
      double cfl_scale_factor) const {
    const core::DirectionalCflTimeStepInput input{
        .cell_width_axis_code = {
            m_geometry.cell_width_x_comoving,
            m_geometry.cell_width_y_comoving,
            m_geometry.cell_width_z_comoving},
        .velocity_axis_code = {
            state.gas_cells.velocity_x_peculiar[cell_index],
            state.gas_cells.velocity_y_peculiar[cell_index],
            state.gas_cells.velocity_z_peculiar[cell_index]},
        .sound_speed_code = std::max(state.gas_cells.sound_speed_code[cell_index], 0.0),
        .coordinate_frame = m_config.units.coordinate_frame,
        .scale_factor = cfl_scale_factor,
        .stability_rule = core::hydroCflStabilityRuleFor(m_config),
        .active_axis = core::activeHydroAxes(
            {m_geometry.cell_width_x_comoving, m_geometry.cell_width_y_comoving, m_geometry.cell_width_z_comoving},
            core::hydroDomainExtentCode(m_config)),
    };
    const auto patch_identity = patchIdentityForCellRow(state, cell_index);
    return core::makeHydroCflDiagnostics(
        cell_index,
        input,
        0.4,
        accepted_dt_time_code,
        state.gas_cell_identity.findByLocalRow(cell_index) != nullptr
            ? state.gas_cell_identity.findByLocalRow(cell_index)->gas_cell_id
            : 0U,
        patch_identity.has_value() ? std::optional<std::uint64_t>(patch_identity->first) : std::nullopt,
        patch_identity.has_value() ? std::optional<std::uint32_t>(patch_identity->second) : std::nullopt);
  }

  void verifyAcceptedHydroCfl(
      core::StepContext& context,
      std::span<const std::uint32_t> active_cells) {
    if (context.integrator_state.dt_time_code <= 0.0) {
      throw std::invalid_argument("hydro CFL guard requires dt_time_code > 0");
    }
    const std::size_t cell_count = context.state.cells.size();
    if (cell_count == 0) {
      return;
    }

    std::vector<std::uint32_t> all_cells;
    std::span<const std::uint32_t> cells_to_check = active_cells;
    if (cells_to_check.empty()) {
      all_cells.resize(cell_count);
      std::iota(all_cells.begin(), all_cells.end(), 0U);
      cells_to_check = all_cells;
    }

    core::HydroCflDiagnostics worst;
    bool have_worst = false;
    for (const std::size_t cell : cells_to_check) {
      if (cell >= cell_count) {
        throw std::out_of_range("hydro CFL guard active cell index out of range");
      }
      const double cfl_scale_factor = m_config.units.coordinate_frame == core::CoordinateFrame::kComoving
          ? std::max(context.timeline_step.scale_factor_begin, 1.0e-12)
          : 1.0;
      const core::HydroCflDiagnostics diagnostics = hydroCflDiagnosticsForCell(
          context.state,
          static_cast<std::uint32_t>(cell),
          context.integrator_state.dt_time_code,
          cfl_scale_factor);
      if (!have_worst || diagnostics.safety_factor < worst.safety_factor) {
        worst = diagnostics;
        have_worst = true;
      }
    }
    if (!have_worst) {
      return;
    }
    m_last_hydro_cfl_diagnostics = worst;
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.cfl_guard",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.godunov",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "checked accepted hydro timestep against local directional CFL bound",
          .payload = {{"local_row", std::to_string(worst.local_row)},
                      {"gas_cell_id", worst.has_gas_cell_id ? std::to_string(worst.gas_cell_id) : "n/a"},
                      {"patch_id", worst.has_patch_id ? std::to_string(worst.patch_id) : "n/a"},
                      {"patch_row", worst.has_patch_row ? std::to_string(worst.patch_row) : "n/a"},
                      {"proposed_dt_time_code", std::to_string(worst.proposed_dt_time_code)},
                      {"accepted_dt_time_code", std::to_string(worst.accepted_dt_time_code)},
                      {"cfl_number", std::to_string(worst.cfl_number)},
                      {"safety_factor", std::to_string(worst.safety_factor)},
                      {"velocity_x_code", std::to_string(worst.velocity_axis_code[0])},
                      {"velocity_y_code", std::to_string(worst.velocity_axis_code[1])},
                      {"velocity_z_code", std::to_string(worst.velocity_axis_code[2])},
                      {"sound_speed_code", std::to_string(worst.sound_speed_code)}}});
    }
    core::assertHydroCflStable(worst);
  }

  const core::SimulationConfig& m_config;
  const core::ModePolicy& m_mode_policy;
  const parallel::MpiContext& m_mpi_context;
  core::MemoryGovernor* m_memory_governor = nullptr;
  core::UnitSystem m_units;
  std::uint32_t m_current_world_rank = 0;
  double m_temperature_k_per_internal_energy_code = 0.0;
  hydro::HydroCoreSolver m_solver;
  std::unique_ptr<hydro::HydroReconstruction> m_reconstruction;
  std::unique_ptr<hydro::HydroRiemannSolver> m_riemann_solver;
  std::unique_ptr<physics::CoolingHeatingSource> m_cooling_source;
  std::optional<physics::CoolingRateProvider> m_cooling_rate_provider;
  std::uint64_t m_last_predictive_forecast_events = 0U;
  std::size_t m_last_regrid_optional_budget = 0U;
  double m_cooling_specific_rate_cgs_to_code = 0.0;
  double m_cooling_temperature_per_energy_k = 0.0;
  double m_cooling_helium_mass_fraction = 0.24;
  physics::CompositionTransportCommitReport m_last_composition_transport{};
  SourceRuntime* m_source_runtime = nullptr;  // non-owning, composition lifetime
  std::unique_ptr<physics::EffectiveIsmThermodynamicClosure> m_effective_ism_closure;
  std::unique_ptr<physics::EffectiveIsmEnergyRelaxationSource> m_effective_ism_energy_source;
  hydro::HydroConservedStateSoa m_conserved;
  hydro::HydroScratchBuffers m_scratch;
  hydro::HydroPatchGeometry m_geometry;
  // Dense CellSoa row <-> physical Cartesian solver-row maps.
  std::vector<std::uint32_t> m_dense_row_by_geometry_row;
  std::vector<std::uint32_t> m_geometry_row_by_dense_row;
  hydro::HydroProfileEvent m_last_hydro_profile{};
  core::HydroCflDiagnostics m_last_hydro_cfl_diagnostics{};
  internal::SolverGhostRefreshReport m_last_ghost_refresh{};
  HydroRemoteGhostBoundaryReport m_last_remote_ghost_boundary_report{};
  parallel::DirectedAmrExchangeDiagnostics m_last_directed_amr_diagnostics{};
  amr::ProductionAmrHydroDiagnostics m_last_amr_hydro_diagnostics{};
  amr::ProductionAmrRegridSummary m_last_regrid_summary{};
  std::uint64_t m_communication_high_water_bytes = 0U;
  std::uint64_t m_last_communication_traffic_bytes = 0U;
  std::uint64_t m_remote_ghost_sparse_cell_high_water = 0U;
  std::uint64_t m_remote_ghost_capacity_high_water_bytes = 0U;
  std::uint64_t m_remote_ghost_metadata_high_water_bytes = 0U;
  HydroFloorApplicationDiagnostics m_last_floor_diagnostics{};
  parallel::GhostCacheLifecycle m_ghost_cache_lifecycle{};
  std::vector<std::size_t> m_active_cells;
  std::vector<std::size_t> m_active_faces;
  std::unordered_set<std::size_t> m_remote_non_authority_faces;
  std::optional<HydroGeometryCacheKey> m_cached_geometry_key;
};
}  // namespace

std::unique_ptr<HydroAmrRuntime> makeHydroAmrRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const RuntimeServices& services,
    std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table) {
  return std::make_unique<HydroAmrRuntimeImpl>(
      config, mode_policy, services, std::move(effective_eos_table));
}

}  // namespace cosmosim::workflows
