#include "cosmosim/workflows/gravity_runtime.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/cuda_runtime.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/openmp_runtime.hpp"
#include "cosmosim/workflows/runtime_module_registry.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/gravity/gravity_memory.hpp"
#include "cosmosim/gravity/gravity_source_snapshot.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/particle_ghost_runtime.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/gravity_force_endpoint.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::workflows {
namespace {

[[nodiscard]] std::string formatRuntimeDouble(double value) {
  std::ostringstream stream;
  stream << std::scientific
         << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
  return stream.str();
}

constexpr std::array<core::IntegrationStage, 1> k_hydro_stage = {core::IntegrationStage::kHydroUpdate};

// Rank-agreement checks gathered into one collective (2026-10-25; hardened
// 2026-10-26). Each equality entry v occupies lanes v and ~v of a single
// MPI_MIN allreduce; ~ is an order-reversing bijection on uint64, so
// min(~v) = ~max(v) and min(v) == max(v) iff every rank holds v.
// Min-only entries return the global minimum.
//
// The lane layout has a fixed capacity, so the MPI element count never
// depends on how many entries a rank appended. Structural lanes make a
// mismatched entry count / order (name fingerprint), a capacity overflow and
// rank-local failures captured by recordLocalFailure() visible on every rank:
// verify() then throws identically everywhere (the failing rank rethrows its
// own exception) instead of deadlocking or reducing misaligned lanes.
class RankConsensusBatch {
 public:
  static constexpr std::size_t k_capacity = 32U;

  void requireEqual(std::string_view name, std::uint64_t value) { append(name, value, true); }
  // Returns an index into the verified results for a global-minimum entry.
  std::size_t globalMinimum(std::string_view name, std::uint64_t value) {
    append(name, value, false);
    return m_count - 1U;
  }
  // Captures a rank-local failure raised while the batched facts were being
  // assembled; reported collectively by verify().
  void recordLocalFailure(std::exception_ptr failure) noexcept {
    if (!m_local_failure) m_local_failure = std::move(failure);
  }

  void verify(const parallel::MpiContext& mpi_context, std::uint64_t world_size) {
    constexpr std::size_t c = k_capacity;
    const std::size_t used = std::min(m_count, c);
    // [values | complements | count, ~count | fingerprint, ~fingerprint |
    //  ~overflow | ~failure | failing rank (UINT64_MAX when none)]
    std::array<std::uint64_t, 2U * c + 7U> lanes{};
    for (std::size_t i = 0; i < c; ++i) {
      const std::uint64_t value = i < used ? m_entries[i].value : 0U;
      lanes[i] = value;
      lanes[c + i] = ~value;
    }
    const std::uint64_t rank = static_cast<std::uint64_t>(std::max(mpi_context.worldRank(), 0));
    lanes[2U * c + 0U] = m_count;
    lanes[2U * c + 1U] = ~static_cast<std::uint64_t>(m_count);
    lanes[2U * c + 2U] = m_fingerprint;
    lanes[2U * c + 3U] = ~m_fingerprint;
    lanes[2U * c + 4U] = ~static_cast<std::uint64_t>(m_count > c ? 1U : 0U);
    lanes[2U * c + 5U] = ~static_cast<std::uint64_t>(m_local_failure ? 1U : 0U);
    lanes[2U * c + 6U] = m_local_failure ? rank : std::numeric_limits<std::uint64_t>::max();
    mpi_context.allreduceMinUint64sInPlace(lanes);
    if (~lanes[2U * c + 5U] != 0U) {
      if (m_local_failure) std::rethrow_exception(m_local_failure);
      throw std::runtime_error(
          "TreePM rank-consensus batch: a peer rank failed while assembling consensus facts "
          "(first failing rank=" + std::to_string(lanes[2U * c + 6U]) +
          ", world_size=" + std::to_string(world_size) + ")");
    }
    if (lanes[2U * c + 0U] != ~lanes[2U * c + 1U] || lanes[2U * c + 2U] != ~lanes[2U * c + 3U]) {
      throw std::logic_error(
          "TreePM rank-consensus batch structure differs across ranks (entry count or order); local entries=" +
          std::to_string(m_count));
    }
    if (~lanes[2U * c + 4U] != 0U) {
      throw std::logic_error(
          "TreePM rank-consensus batch exceeds its fixed capacity of " + std::to_string(c) + " entries");
    }
    for (std::size_t i = 0; i < used; ++i) {
      const std::uint64_t global_min = lanes[i];
      const std::uint64_t global_max = ~lanes[c + i];
      if (m_entries[i].equality && global_min != global_max) {
        throw std::runtime_error(
            "TreePM cadence rank-consensus failure for " + std::string(m_entries[i].name) +
            ": local=" + std::to_string(m_entries[i].value) +
            ", global_min=" + std::to_string(global_min) +
            ", global_max=" + std::to_string(global_max) +
            ", world_size=" + std::to_string(world_size));
      }
      m_minimum[i] = global_min;
      m_maximum[i] = global_max;
    }
    m_verified = true;
  }
  [[nodiscard]] std::uint64_t minimum(std::size_t index) const {
    requireVerified(index);
    return m_minimum[index];
  }
  [[nodiscard]] std::uint64_t maximum(std::size_t index) const {
    requireVerified(index);
    return m_maximum[index];
  }

 private:
  struct Entry {
    std::string_view name;
    std::uint64_t value = 0U;
    bool equality = true;
  };
  void append(std::string_view name, std::uint64_t value, bool equality) noexcept {
    // FNV-1a over (name, kind) in append order: detects reordered or
    // differently named entry sequences across ranks.
    for (const char ch : name) {
      m_fingerprint = (m_fingerprint ^ static_cast<unsigned char>(ch)) * 1099511628211ULL;
    }
    m_fingerprint = (m_fingerprint ^ (equality ? 0x45U : 0x4DU)) * 1099511628211ULL;
    if (m_count < k_capacity) m_entries[m_count] = Entry{name, value, equality};
    ++m_count;  // overflow is reported collectively by verify()
  }
  void requireVerified(std::size_t index) const {
    if (!m_verified || index >= std::min(m_count, k_capacity)) {
      throw std::logic_error("rank consensus result requested before verification");
    }
  }
  std::array<Entry, k_capacity> m_entries{};
  std::array<std::uint64_t, k_capacity> m_minimum{};
  std::array<std::uint64_t, k_capacity> m_maximum{};
  std::size_t m_count = 0U;
  std::uint64_t m_fingerprint = 14695981039346656037ULL;
  std::exception_ptr m_local_failure;
  bool m_verified = false;
};

[[nodiscard]] gravity::PmAssignmentScheme toPmAssignmentScheme(
    core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic:
      return gravity::PmAssignmentScheme::kCic;
    case core::TreePmAssignmentScheme::kTsc:
      return gravity::PmAssignmentScheme::kTsc;
  }
  throw std::runtime_error("unhandled TreePm assignment scheme enum value");
}

[[nodiscard]] std::string treePmAssignmentSchemeName(core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic:
      return "cic";
    case core::TreePmAssignmentScheme::kTsc:
      return "tsc";
  }
  throw std::runtime_error("unhandled TreePm assignment scheme enum value");
}

[[nodiscard]] gravity::TreeOpeningCriterion toTreeOpeningCriterion(
    core::TreePmOpeningCriterion opening_criterion) {
  switch (opening_criterion) {
    case core::TreePmOpeningCriterion::kGeometric:
      return gravity::TreeOpeningCriterion::kBarnesHutGeometric;
    case core::TreePmOpeningCriterion::kComDistance:
      return gravity::TreeOpeningCriterion::kBarnesHutComDistance;
    case core::TreePmOpeningCriterion::kRelativeForceError:
      return gravity::TreeOpeningCriterion::kRelativeForceError;
  }
  throw std::runtime_error("unhandled TreePM opening criterion enum value");
}

[[nodiscard]] std::string treePmOpeningCriterionName(
    core::TreePmOpeningCriterion opening_criterion) {
  switch (opening_criterion) {
    case core::TreePmOpeningCriterion::kGeometric:
      return "geometric";
    case core::TreePmOpeningCriterion::kComDistance:
      return "com_distance";
    case core::TreePmOpeningCriterion::kRelativeForceError:
      return "relative_force_error";
  }
  throw std::runtime_error("unhandled TreePM opening criterion enum value");
}



gravity::TreePmCoordinator makeRuntimeAwareTreePmCoordinator(
    const core::SimulationConfig& config,
    const gravity::PmGridShape& pm_grid_shape,
    const RuntimeServices& services) {
  const parallel::MpiContext& mpi_context = services.mpi_context;
  mpi_context.validateExpectedWorldSizeOrThrow(config.parallel.mpi_ranks_expected);
  const auto layout = parallel::makePmSlabLayout(
      pm_grid_shape.nx,
      pm_grid_shape.ny,
      pm_grid_shape.nz,
      mpi_context.worldSize(),
      mpi_context.worldRank());
  return gravity::TreePmCoordinator(
      pm_grid_shape, layout, mpi_context, services.memory_governor);
}


[[nodiscard]] double newtonGCodeFromUnits(const core::UnitSystem& units) {
  return core::newtonGravitationalConstantCode(units);
}

[[nodiscard]] std::string pmDecompositionModeName(core::PmDecompositionMode mode) {
  switch (mode) {
    case core::PmDecompositionMode::kSlab:
      return "slab";
    case core::PmDecompositionMode::kPencil:
      return "pencil";
  }
  throw std::runtime_error("unhandled PM decomposition mode enum value");
}

[[nodiscard]] gravity::TreeSofteningSpeciesPolicy speciesSofteningByTag(
    const core::SimulationConfig& config) {
  gravity::TreeSofteningSpeciesPolicy values{};
  values.enabled = config.numerics.gravity_softening_gas_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_star_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_tracer_kpc_comoving > 0.0;
  values.epsilon_comoving_by_species.fill(config.numerics.gravity_softening_kpc_comoving * 1.0e-3);
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)] =
      (config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_dark_matter_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)] =
      (config.numerics.gravity_softening_gas_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_gas_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)] =
      (config.numerics.gravity_softening_star_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_star_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)] =
      (config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_black_hole_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kTracer)] =
      (config.numerics.gravity_softening_tracer_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_tracer_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kTracer)];
  return values;
}

[[nodiscard]] bool hasSpeciesSpecificSoftening(const core::SimulationConfig& config) {
  return config.numerics.gravity_softening_gas_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_star_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_tracer_kpc_comoving > 0.0;
}

[[nodiscard]] std::string describeSofteningPolicy(
    const core::SimulationConfig& config,
    const core::SimulationState* state = nullptr) {
  const bool per_species = hasSpeciesSpecificSoftening(config);
  const bool per_particle = state != nullptr && std::any_of(
      state->particle_sidecar.has_gravity_softening_override.begin(),
      state->particle_sidecar.has_gravity_softening_override.end(),
      [](std::uint8_t flag) { return flag != 0U; });
  if (per_particle && per_species) {
    return "comoving_species_plus_particle_override";
  }
  if (per_particle) {
    return "comoving_particle_override";
  }
  if (per_species) {
    return "comoving_species";
  }
  return "comoving_fixed";
}

void maybeInitializeParticleSofteningFromSpeciesPolicy(
    core::SimulationState& state,
    const core::SimulationConfig& config) {
  if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
    return;
  }
  if (!hasSpeciesSpecificSoftening(config)) {
    return;
  }
  const auto by_species = speciesSofteningByTag(config);
  state.particle_sidecar.gravity_softening_comoving.resize(state.particles.size(), 0.0);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const std::size_t species_tag = static_cast<std::size_t>(state.particleSpeciesTag(i));
    if (species_tag < by_species.epsilon_comoving_by_species.size()) {
      state.particle_sidecar.gravity_softening_comoving[i] = by_species.epsilon_comoving_by_species[species_tag];
    } else {
      state.particle_sidecar.gravity_softening_comoving[i] = config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
    }
  }
}

[[nodiscard]] double wrapPeriodicPosition(double position_comoving, double box_size_mpc_comoving) {
  if (box_size_mpc_comoving <= 0.0) {
    return position_comoving;
  }
  double wrapped = std::fmod(position_comoving, box_size_mpc_comoving);
  if (wrapped < 0.0) {
    wrapped += box_size_mpc_comoving;
  }
  if (wrapped >= box_size_mpc_comoving) {
    wrapped -= box_size_mpc_comoving;
  }
  return wrapped;
}

[[nodiscard]] bool hasHdf5Extension(const std::filesystem::path& path) {
  const std::string ext = path.extension().string();
  return ext == ".h5" || ext == ".hdf5";
}

[[nodiscard]] std::unordered_set<std::uint64_t> loadZoomParticleIdsFromText(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open zoom region file: " + path.string());
  }
  std::unordered_set<std::uint64_t> ids;
  std::string line;
  while (std::getline(in, line)) {
    const auto comment = line.find('#');
    if (comment != std::string::npos) {
      line.erase(comment);
    }
    std::istringstream stream(line);
    std::uint64_t id = 0;
    while (stream >> id) {
      ids.insert(id);
    }
  }
  return ids;
}

#if COSMOSIM_ENABLE_HDF5
[[nodiscard]] std::unordered_set<std::uint64_t> loadZoomParticleIdsFromHdf5(const std::filesystem::path& path) {
  hid_t file = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    throw std::runtime_error("failed to open zoom region HDF5 file: " + path.string());
  }
  auto close_file = [&]() { H5Fclose(file); };
  const char* candidates[] = {"/Region/ParticleIDs", "/ParticleIDs", "Region/ParticleIDs", "ParticleIDs"};
  hid_t dataset = -1;
  for (const char* candidate : candidates) {
    if (H5Lexists(file, candidate, H5P_DEFAULT) > 0) {
      dataset = H5Dopen2(file, candidate, H5P_DEFAULT);
      if (dataset >= 0) {
        break;
      }
    }
  }
  if (dataset < 0) {
    close_file();
    throw std::runtime_error("zoom region HDF5 file lacks ParticleIDs dataset: " + path.string());
  }
  hid_t space = H5Dget_space(dataset);
  if (space < 0) {
    H5Dclose(dataset); close_file();
    throw std::runtime_error("failed to query zoom region HDF5 dataspace");
  }
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_ndims(space) != 1 || H5Sget_simple_extent_dims(space, dims, nullptr) != 1) {
    H5Sclose(space); H5Dclose(dataset); close_file();
    throw std::runtime_error("zoom region ParticleIDs dataset must be 1D");
  }
  std::vector<std::uint64_t> values(static_cast<std::size_t>(dims[0]), 0ULL);
  const hid_t dtype = H5Dget_type(dataset);
  const H5T_class_t cls = H5Tget_class(dtype);
  if (cls != H5T_INTEGER) {
    H5Tclose(dtype); H5Sclose(space); H5Dclose(dataset); close_file();
    throw std::runtime_error("zoom region ParticleIDs dataset must be integer typed");
  }
  if (!values.empty() && H5Dread(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    H5Tclose(dtype); H5Sclose(space); H5Dclose(dataset); close_file();
    throw std::runtime_error("failed to read zoom region ParticleIDs dataset");
  }
  H5Tclose(dtype); H5Sclose(space); H5Dclose(dataset); close_file();
  return std::unordered_set<std::uint64_t>(values.begin(), values.end());
}
#endif


class GravityRuntimeImpl final : public GravityRuntime {
 public:
  struct GravityHealthSummary {
    std::uint64_t cheap_checks_executed = 0;
    std::uint64_t heavy_checks_executed = 0;
    std::uint64_t warning_count = 0;
    std::uint64_t fatal_count = 0;
    std::uint64_t pm_field_non_finite_count = 0;
    std::uint64_t force_non_finite_count = 0;
    std::uint64_t force_abnormal_count = 0;
    std::uint64_t illegal_sync_state_count = 0;
    std::uint64_t zoom_sanity_failure_count = 0;
    std::uint64_t decomposition_sanity_failure_count = 0;
  };

  enum class PmSyncSurface : std::uint8_t {
    kKickPre = 0,
    kForceRefresh = 1,
    kKickPost = 2,
  };

  struct PmLongRangeKickDecision {
    PmSyncSurface sync_surface = PmSyncSurface::kKickPre;
    std::uint64_t gravity_kick_opportunity = 0;
    bool refresh_long_range_field = false;
    std::uint64_t field_version = 0;
    std::uint64_t last_refresh_opportunity = 0;
    std::uint64_t field_built_step_index = 0;
    double field_built_scale_factor = 1.0;
  };

  GravityRuntimeImpl(
      const core::SimulationConfig& config,
      const core::ModePolicy& mode_policy,
      const RuntimeServices& services,
      std::filesystem::path zoom_region_path = {})
      : m_services(services),
        m_config(config),
        m_mode_policy(mode_policy),
        m_pm_update_cadence_steps(static_cast<std::uint64_t>(config.numerics.treepm_update_cadence_steps)),
        m_pm_grid_shape(gravity::PmGridShape{
            static_cast<std::size_t>(config.numerics.treepm_pm_grid_nx),
            static_cast<std::size_t>(config.numerics.treepm_pm_grid_ny),
            static_cast<std::size_t>(config.numerics.treepm_pm_grid_nz)}),
        m_mesh_spacing_x_mpc_comoving(
            config.cosmology.box_size_x_mpc_comoving / static_cast<double>(m_pm_grid_shape.nx)),
        m_mesh_spacing_y_mpc_comoving(
            config.cosmology.box_size_y_mpc_comoving / static_cast<double>(m_pm_grid_shape.ny)),
        m_mesh_spacing_z_mpc_comoving(
            config.cosmology.box_size_z_mpc_comoving / static_cast<double>(m_pm_grid_shape.nz)),
        m_tree_pm_coordinator(makeRuntimeAwareTreePmCoordinator(
            config, m_pm_grid_shape, services)),
        m_zoom_region_path(std::move(zoom_region_path)) {
    if (!m_pm_grid_shape.isValid()) {
      throw std::runtime_error("numerics.treepm_pm_grid_n{xyz} must all be > 0");
    }
    m_tree_pm_options.pm_options.box_size_x_mpc_comoving = config.cosmology.box_size_x_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_y_mpc_comoving = config.cosmology.box_size_y_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_z_mpc_comoving = config.cosmology.box_size_z_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_mpc_comoving = config.cosmology.box_size_x_mpc_comoving;
    m_tree_pm_options.pm_options.scale_factor = 1.0;
    const core::UnitSystem runtime_units = core::makeUnitSystem(
        config.units.length_unit,
        config.units.mass_unit,
        config.units.velocity_unit);
    m_gravitational_constant_code =
        core::newtonGravitationalConstantCode(runtime_units);
    m_tree_pm_options.pm_options.gravitational_constant_code =
        m_gravitational_constant_code;
    m_tree_pm_options.pm_options.assignment_scheme =
        toPmAssignmentScheme(config.numerics.treepm_assignment_scheme);
    m_tree_pm_options.pm_options.force_gradient_backend =
        config.numerics.treepm_pm_force_gradient == core::TreePmForceGradient::kFiniteDifference4
        ? gravity::PmForceGradientBackend::kFiniteDifference4
        : gravity::PmForceGradientBackend::kSpectral;
    m_tree_pm_options.pm_options.enable_window_deconvolution =
        config.numerics.treepm_enable_window_deconvolution;
    m_tree_pm_options.pm_options.decomposition_mode = config.numerics.treepm_pm_decomposition_mode;
    m_tree_pm_options.pm_options.boundary_condition = mode_policy.gravity_boundary ==
            core::GravityBoundaryModel::kPeriodicPoisson
        ? gravity::PmBoundaryCondition::kPeriodic
        : gravity::PmBoundaryCondition::kIsolatedOpen;
    m_tree_pm_options.pm_options.isolated_open_root_workspace_limit_bytes =
        config.parallel.isolated_pm_root_workspace_limit_bytes;
    const parallel::MpiContext& mpi_context = m_services.mpi_context;
    const core::CudaRuntimeInfo cuda_runtime = core::queryCudaRuntime();
    m_runtime_topology = parallel::buildDistributedExecutionTopology(
        m_pm_grid_shape.nx,
        m_pm_grid_shape.ny,
        m_pm_grid_shape.nz,
        mpi_context,
        config.parallel.mpi_ranks_expected,
        config.parallel.gpu_devices,
        cuda_runtime.runtime_available,
        cuda_runtime.visible_device_count,
        pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode));
    if (m_runtime_topology.usesCuda()) {
      core::setCudaDeviceOrThrow(m_runtime_topology.device_assignment.assigned_device_index);
      m_tree_pm_options.pm_options.execution_policy = core::ExecutionPolicy::kCuda;
      m_tree_pm_options.pm_options.data_residency = gravity::PmDataResidencyPolicy::kPreferDevice;
    }

    m_tree_pm_options.tree_options.opening_criterion =
        toTreeOpeningCriterion(config.numerics.treepm_tree_opening_criterion);
    m_tree_pm_options.tree_options.opening_theta = config.numerics.treepm_tree_opening_theta;
    m_tree_pm_options.tree_options.relative_force_tolerance =
        config.numerics.treepm_tree_relative_force_tolerance;
    m_tree_pm_options.tree_options.relative_force_acceleration_floor_code =
        config.numerics.treepm_tree_relative_force_acceleration_floor;
    m_tree_pm_options.spatial_work_history_enabled = config.parallel.decomposition_spatial_work_enabled;
    m_tree_pm_options.topology_refit_enabled = config.numerics.treepm_topology_refit_enabled;
    m_tree_pm_options.identical_source_tree_reuse_enabled =
        config.numerics.treepm_identical_source_tree_reuse_enabled;
    m_tree_pm_options.acceptance_policy = config.numerics.treepm_adaptive_acceptance_enabled
        ? gravity::TreePmAcceptancePolicy::kAdaptiveRelative
        : gravity::TreePmAcceptancePolicy::kStrictReference;
    m_tree_pm_options.gaussian_pair_lookup_enabled = config.numerics.treepm_gaussian_pair_lookup_enabled;
    m_tree_pm_options.full_residual_lookup_enabled = config.numerics.treepm_full_residual_lookup_enabled;
    m_tree_pm_options.leaf_simd_enabled = config.numerics.treepm_leaf_simd_enabled;
    m_tree_pm_options.radial_coefficients_enabled = config.numerics.treepm_radial_coefficients_enabled;
    m_tree_pm_options.residual_block_size = static_cast<std::size_t>(config.numerics.treepm_residual_block_size);
    m_tree_pm_options.tree_options.max_leaf_size = static_cast<std::size_t>(config.numerics.treepm_tree_max_leaf_size);
    m_tree_pm_options.full_mac_diagnostics = config.numerics.treepm_full_mac_diagnostics;
    m_tree_pm_options.adaptive_maximum_opening_angle =
        config.numerics.treepm_adaptive_maximum_opening_angle;
    m_tree_pm_options.adaptive_error_estimator =
        config.numerics.treepm_adaptive_error_estimator == core::TreePmAdaptiveErrorEstimator::kRelativeMonopole
        ? gravity::TreePmAdaptiveErrorEstimator::kRelativeMonopole
        : gravity::TreePmAdaptiveErrorEstimator::kResidualDerivative;
    m_tree_pm_options.cutoff_node_policy =
        config.numerics.treepm_cutoff_node_policy == core::TreePmCutoffNodePolicy::kComDistance
        ? gravity::TreePmCutoffNodePolicy::kComDistance
        : gravity::TreePmCutoffNodePolicy::kContainNode;
    m_tree_pm_options.group_walk_enabled = config.numerics.treepm_group_walk_enabled;
    m_tree_pm_options.allow_stale_long_range_field =
        config.numerics.treepm_pm_reuse_policy == core::TreePmPmReusePolicy::kStaleWithinCadence;
    m_tree_pm_options.group_walk_target_group_size =
        static_cast<std::size_t>(config.numerics.treepm_group_walk_target_size);
    m_tree_pm_options.tree_options.gravitational_constant_code =
        m_gravitational_constant_code;
    m_tree_pm_options.tree_options.softening.kernel = gravity::TreeSofteningKernel::kPlummer;
    m_tree_pm_options.tree_options.softening.epsilon_comoving =
        config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
    m_tree_pm_species_softening = speciesSofteningByTag(config);
    m_tree_pm_options.split_policy = gravity::makeTreePmSplitPolicyFromMeshSpacing(
        config.numerics.treepm_asmth_cells,
        config.numerics.treepm_rcut_cells,
        std::cbrt(
            m_mesh_spacing_x_mpc_comoving * m_mesh_spacing_y_mpc_comoving *
            m_mesh_spacing_z_mpc_comoving));
    m_tree_pm_options.enable_zoom_long_range_correction =
        config.mode.zoom_long_range_strategy ==
        core::ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection;
    m_tree_pm_options.zoom_focused_pm_shape = gravity::PmGridShape{
        static_cast<std::size_t>(std::max(config.mode.zoom_focused_pm_grid_nx, 0)),
        static_cast<std::size_t>(std::max(config.mode.zoom_focused_pm_grid_ny, 0)),
        static_cast<std::size_t>(std::max(config.mode.zoom_focused_pm_grid_nz, 0))};
    m_tree_pm_options.zoom_region_center_x_comoving = config.mode.zoom_region_center_x_mpc_comoving;
    m_tree_pm_options.zoom_region_center_y_comoving = config.mode.zoom_region_center_y_mpc_comoving;
    m_tree_pm_options.zoom_region_center_z_comoving = config.mode.zoom_region_center_z_mpc_comoving;
    m_tree_pm_options.zoom_region_radius_comoving = config.mode.zoom_region_radius_mpc_comoving;
    m_tree_pm_options.zoom_contamination_radius_comoving =
        config.mode.zoom_contamination_radius_mpc_comoving > 0.0
        ? config.mode.zoom_contamination_radius_mpc_comoving
        : config.mode.zoom_region_radius_mpc_comoving;
    m_tree_pm_options.tree_exchange_batch_bytes = config.numerics.treepm_tree_exchange_batch_bytes;
    m_tree_pm_options.pm_options.routing_exchange_batch_bytes =
        config.numerics.treepm_pm_exchange_batch_bytes;
    m_tree_pm_options.zoom_high_res_allgather_limit_bytes =
        config.parallel.zoom_high_res_allgather_limit_bytes;
    m_pm_assignment_scheme = treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme);
    m_pm_backend = gravity::PmSolver::fftBackendName();
    if (m_runtime_topology.usesCuda()) {
      m_pm_backend += "+cuda_cic";
    }
  }

  [[nodiscard]] std::size_t pmGridSize() const noexcept { return m_pm_grid_shape.nx; }
  [[nodiscard]] const gravity::PmGridShape& pmGridShape() const noexcept { return m_pm_grid_shape; }
  [[nodiscard]] std::uint64_t longRangeRefreshCount() const noexcept { return m_long_range_refresh_count; }
  [[nodiscard]] std::uint64_t longRangeReuseCount() const noexcept { return m_long_range_reuse_count; }
  [[nodiscard]] int pmCadenceSteps() const noexcept { return static_cast<int>(m_pm_update_cadence_steps); }
  [[nodiscard]] double gravitationalConstantCode() const noexcept {
    return m_gravitational_constant_code;
  }
  [[nodiscard]] std::span<const ReferenceWorkflowReport::TreePmCadenceRecord> cadenceRecords() const noexcept {
    return m_cadence_records;
  }
  void appendAcceptedCadenceRecord(ReferenceWorkflowReport::TreePmCadenceRecord record) {
    if (m_cadence_records.size() == k_recent_cadence_limit) {
      m_cadence_records.erase(m_cadence_records.begin());
    }
    m_cadence_records.push_back(std::move(record));
  }
  void appendCadenceRecord(ReferenceWorkflowReport::TreePmCadenceRecord record) {
    if (!m_cycle_transaction_open) {
      appendAcceptedCadenceRecord(std::move(record));
      return;
    }
    if (m_pending_cadence_records.size() == k_recent_cadence_limit) {
      m_pending_cadence_records.erase(m_pending_cadence_records.begin());
      ++m_pending_cadence_records_dropped;
    }
    m_pending_cadence_records.push_back(std::move(record));
  }
  void beginCycleTransaction() override {
    if (m_cycle_transaction_open) {
      throw std::logic_error("multirate PM cycle transaction is already open");
    }
    // Reserve both bounded windows here so commit (noexcept) never allocates.
    m_cadence_records.reserve(k_recent_cadence_limit);
    m_pending_cadence_records.reserve(k_recent_cadence_limit);
    m_cycle_transaction_open = true;
    m_pending_long_range_refreshes = 0U;
    m_pending_long_range_reuses = 0U;
    m_pending_cadence_records.clear();
    m_pending_cadence_records_dropped = 0U;
  }
  void commitCycleTransaction() noexcept override {
    if (!m_cycle_transaction_open) return;
    // Publish the journal in order through the same bounded window the
    // non-transactional path uses (oldest accepted records are evicted only
    // by newer accepted records).
    for (auto& record : m_pending_cadence_records) appendAcceptedCadenceRecord(std::move(record));
    m_pending_cadence_records.clear();
    m_pending_long_range_refreshes = 0U;
    m_pending_long_range_reuses = 0U;
    m_pending_cadence_records_dropped = 0U;
    m_cycle_transaction_open = false;
  }
  void rollbackCycleTransaction(
      const TransactionalForceEndpoint& opening_endpoint,
      const core::SimulationState& restored_state,
      std::uint64_t layout_certificate) override {
    if (!m_cycle_transaction_open) {
      throw std::logic_error("multirate PM cycle rollback without an open transaction");
    }
    // 1. Forces first: the row-aligned restore validates every precondition
    //    before its first write, so a refusal leaves the runtime (journal
    //    included) exactly as it was; the caller fails collectively.
    restoreTransactionalForceEndpoint(opening_endpoint, restored_state, layout_certificate);
    // 2. Diagnostics: the attempt's trajectory counts move to execution cost.
    m_long_range_refresh_count -= m_pending_long_range_refreshes;
    m_long_range_reuse_count -= m_pending_long_range_reuses;
    m_rejected_long_range_refresh_count += m_pending_long_range_refreshes;
    m_rejected_long_range_reuse_count += m_pending_long_range_reuses;
    m_rejected_cadence_record_count +=
        static_cast<std::uint64_t>(m_pending_cadence_records.size()) + m_pending_cadence_records_dropped;
    m_pending_cadence_records.clear();
    m_pending_long_range_refreshes = 0U;
    m_pending_long_range_reuses = 0U;
    m_pending_cadence_records_dropped = 0U;
    m_cycle_transaction_open = false;
    // 3. Derived state. Tree build/refit candidates and topology maps depend
    //    on the positions the rejected attempt visited: force a full build.
    m_tree_pm_coordinator.invalidateSourceHistoryForRollback();
    // Particle ghost layer: next solve refreshes it from restored positions.
    m_ghost_cache_lifecycle.valid = false;
    ++m_ghost_cache_lifecycle.invalidation_count;
    // Borrowed-source certificates are keyed by generations that the caller
    // has advanced; clear them explicitly as well (no consumer may rely on
    // the generation check alone).
    m_active_source_indices = {};
    m_borrowed_particle_index_generation = 0U;
    m_borrowed_gravity_source_generation = 0U;
    m_borrowed_eligibility_certificate_valid = false;
    m_source_predicted_inactive_count = 0U;
  }
  [[nodiscard]] bool cycleTransactionOpen() const noexcept override { return m_cycle_transaction_open; }
  [[nodiscard]] std::uint64_t rejectedCadenceRecordCount() const noexcept override {
    return m_rejected_cadence_record_count;
  }
  [[nodiscard]] std::uint64_t rejectedLongRangeRefreshCount() const noexcept override {
    return m_rejected_long_range_refresh_count;
  }
  [[nodiscard]] std::uint64_t rejectedLongRangeReuseCount() const noexcept override {
    return m_rejected_long_range_reuse_count;
  }
  [[nodiscard]] const parallel::DistributedExecutionTopology& runtimeTopology() const noexcept { return m_runtime_topology; }
  [[nodiscard]] core::MemoryReport memoryReport() const {
    core::MemoryReportBuilder builder;
    const auto add = [&builder](
                         core::MemorySubsystem subsystem,
                         core::MemoryLifetime lifetime,
                         std::string label,
                         const auto& container) {
      const std::uint64_t capacity_bytes = core::ownedCapacityBytesForContainer(container);
      builder.addEntry(core::MemoryEntry{
          .subsystem = subsystem,
          .lifetime = lifetime,
          .memory_class = lifetime == core::MemoryLifetime::kPersistent
                              ? core::MemoryClass::kPersistentCache
                              : core::MemoryClass::kPhaseResident,
          .label = std::move(label),
          .current_size_bytes = core::currentSizeBytesForContainer(container),
          .owned_capacity_bytes = capacity_bytes,
          .high_water_bytes = capacity_bytes,
          .estimated_next_step_bytes = 0U,
          .uncertainty_note = "retained vector capacity is the observed allocation high-water; next-step requirement is not predicted",
      });
    };
    const auto active = [&](std::string label, const auto& container) {
      add(core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
          std::move(label), container);
    };
    const auto persistent = [&](std::string label, const auto& container) {
      add(core::MemorySubsystem::kSidecars, core::MemoryLifetime::kPersistent,
          std::move(label), container);
    };

    active("gravity_runtime.active_accel_x", m_active_accel_x);
    active("gravity_runtime.active_accel_y", m_active_accel_y);
    active("gravity_runtime.active_accel_z", m_active_accel_z);
    active("gravity_runtime.active_previous_acceleration", m_active_previous_acceleration_magnitude);
    persistent("gravity_runtime.recent_cadence_records", m_cadence_records);
    persistent("gravity_runtime.multirate_cadence_journal", m_pending_cadence_records);
    std::uint64_t cadence_string_bytes = 0U;
    for (const auto* window : {&m_cadence_records, &m_pending_cadence_records}) for (const auto& record : *window) {
      for (const auto* value : {&record.stage_name, &record.pm_sync_surface, &record.pm_refresh_reason}) {
        cadence_string_bytes = core::checkedMemoryBytesAdd(cadence_string_bytes,
            core::checkedSizeAdd(value->capacity(), 1U, "cadence string terminator"), "cadence retained strings");
      }
    }
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kSidecars,
        .lifetime = core::MemoryLifetime::kPersistent, .memory_class = core::MemoryClass::kDiagnostic,
        .label = "gravity_runtime.recent_cadence_strings", .current_size_bytes = cadence_string_bytes,
        .owned_capacity_bytes = cadence_string_bytes, .high_water_bytes = cadence_string_bytes,
        .uncertainty_note = "conservative string capacity including SSO; allocator overhead uses process reserve"});
    const std::uint64_t bounded_feedback_bytes = sizeof(m_last_tree_pm_diagnostics) +
        sizeof(m_last_decomposition_measurements.spatial_tree_work_per_target);
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kTree,
        .lifetime = core::MemoryLifetime::kPersistent,
        .memory_class = core::MemoryClass::kDiagnostic,
        .label = "gravity_runtime.bounded_force_feedback",
        .current_size_bytes = bounded_feedback_bytes,
        .owned_capacity_bytes = bounded_feedback_bytes,
        .high_water_bytes = bounded_feedback_bytes});
    persistent("gravity_runtime.particle_accel_x", m_particle_accel_x);
    persistent("gravity_runtime.particle_accel_y", m_particle_accel_y);
    persistent("gravity_runtime.particle_accel_z", m_particle_accel_z);
    persistent("gravity_runtime.cell_accel_x", m_cell_accel_x);
    persistent("gravity_runtime.cell_accel_y", m_cell_accel_y);
    persistent("gravity_runtime.cell_accel_z", m_cell_accel_z);
    active("gravity_runtime.active_slot_by_particle", m_active_slot_by_particle);
    active("gravity_runtime.active_slot_by_cell", m_active_slot_by_cell);
    active("gravity_runtime.owned_local_index_by_particle", m_owned_local_index_by_global);
    active("gravity_runtime.owned_local_index_by_cell", m_owned_local_index_by_cell);
    active("gravity_runtime.owned_leaf_cell_mask", m_owned_leaf_cell_mask);
    active("gravity_runtime.source_x", m_local_source_x);
    active("gravity_runtime.source_y", m_local_source_y);
    active("gravity_runtime.source_z", m_local_source_z);
    active("gravity_runtime.source_mass", m_local_source_mass);
    active("gravity_runtime.source_species", m_local_source_species_tag);
    active("gravity_runtime.source_softening", m_local_source_softening_comoving);
    active("gravity_runtime.source_softening_mask", m_local_source_softening_override_mask);
    active("gravity_runtime.source_particle_row", m_local_source_particle_row);
    active("gravity_runtime.source_cell_row", m_local_source_cell_row);
    active("gravity_runtime.target_particle_row", m_active_target_particle_row);
    active("gravity_runtime.target_cell_row", m_active_target_cell_row);
    active("gravity_runtime.source_high_res_mask", m_source_is_high_res);
    active("gravity_runtime.target_high_res_mask", m_active_is_high_res);
    active("gravity_runtime.local_active_source_index", m_local_active_indices);
    active("gravity_runtime.local_active_particle_row", m_local_active_global_indices);
    active("gravity_runtime.force_refresh_particle_rows", m_force_refresh_particle_indices);
    active("gravity_runtime.authoritative_top_domain_leaves", m_authoritative_top_domain_leaves);
    persistent("gravity_runtime.authoritative_top_domain_seed_leaves", m_authoritative_top_domain_seed_leaves);
    persistent("gravity_runtime.particle_force_cache_valid", m_particle_force_cache_valid);
    persistent("gravity_runtime.hierarchical_force_generation", m_hierarchical_force_generation);
    persistent("gravity_runtime.hierarchical_pm_x", m_hierarchical_pm_x);
    persistent("gravity_runtime.hierarchical_pm_y", m_hierarchical_pm_y);
    persistent("gravity_runtime.hierarchical_pm_z", m_hierarchical_pm_z);
    persistent("gravity_runtime.hierarchical_tree_x", m_hierarchical_tree_x);
    persistent("gravity_runtime.hierarchical_tree_y", m_hierarchical_tree_y);
    persistent("gravity_runtime.hierarchical_tree_z", m_hierarchical_tree_z);
    persistent("gravity_runtime.hierarchical_nearest_pair_r2", m_hierarchical_nearest_r2);
    persistent("gravity_runtime.cell_force_cache_valid", m_cell_force_cache_valid);

    core::MemoryReport owned_report = std::move(builder).finish();
    core::MemoryReport coordinator_report = m_tree_pm_coordinator.memoryReport();
    const std::array<core::MemoryReport, 2> reports{
        std::move(owned_report), std::move(coordinator_report)};
    core::MemoryReport merged = core::mergeMemoryReports(reports);
    if (m_has_pre_run_memory_estimate) {
      merged.notes.push_back(
          "gravity_pre_run_known_peak_bytes=" +
          std::to_string(m_last_pre_run_memory_estimate.known_peak_bytes));
      merged.notes.push_back(
          "gravity_pre_run_external_backend_unknown_bytes=" +
          std::to_string(m_last_pre_run_memory_estimate.external_backend_unknown_bytes));
      merged.notes.push_back(
          "gravity_pre_run_budget_required_bytes=" +
          std::to_string(m_last_pre_run_memory_estimate.budget_required_bytes));
    }
    if (m_has_process_memory_estimate) {
      merged.notes.push_back(
          "dmo_process_pre_run_known_owned_peak_bytes=" +
          std::to_string(m_last_process_memory_estimate.known_owned_peak_bytes));
      merged.notes.push_back(
          "dmo_process_pre_run_external_unknown_reserve_bytes=" +
          std::to_string(m_last_process_memory_estimate.external_unknown_reserve_bytes));
      merged.notes.push_back(
          "dmo_process_pre_run_safety_margin_bytes=" +
          std::to_string(m_last_process_memory_estimate.safety_margin_bytes));
      merged.notes.push_back(
          "dmo_process_pre_run_budget_required_bytes=" +
          std::to_string(m_last_process_memory_estimate.budget_required_bytes));
      merged.notes.push_back(
          "dmo_process_pre_run_aggregate_required_bytes=" +
          std::to_string(m_last_process_memory_estimate.aggregate_required_bytes));
    }
    merged.notes.push_back(
        std::string("gravity_source_representation=") +
        (m_source_representation == gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
             ? "borrowed_homogeneous_dmo"
             : "materialized_generic"));
    return merged;
  }
  [[nodiscard]] const parallel::DecompositionRuntimeMeasurements& lastRuntimeDecompositionMeasurements() const noexcept {
    return m_last_decomposition_measurements;
  }
  void shutdownMpiResources() override {
    m_tree_pm_coordinator.shutdownMpiResources();
  }
  [[nodiscard]] std::uint64_t decompositionEpoch() const noexcept { return m_decomposition_epoch.value; }

  void restoreDecompositionEpoch(std::uint64_t decomposition_epoch) {
    if (m_has_long_range_field || m_last_committed_field_version != 0U) {
      throw std::logic_error("TreePM decomposition epoch must be restored before solver use");
    }
    m_decomposition_epoch = gravity::DecompositionEpoch{decomposition_epoch};
  }

  void publishAuthoritativeTopDomainLeavesSpan() {
    if (!m_authoritative_top_domain_source_generation_valid) {
      // Freshness was invalidated (ownership change before reinstall) or no
      // geometry was ever installed. Retained stale leaves must not be
      // consumed by TreePM; publish an empty routing view so selection falls
      // back conservatively.
      m_tree_pm_options.authoritative_domain_leaves = {};
      m_tree_pm_options.authoritative_geometry_source_generation =
          gravity::GravitySourceGeneration{m_authoritative_top_domain_source_generation};
      return;
    }
    m_tree_pm_options.authoritative_domain_leaves =
        std::span<const parallel::TopDomainLeaf>(
            m_authoritative_top_domain_leaves.data(),
            m_authoritative_top_domain_leaves.size());
    m_tree_pm_options.authoritative_geometry_source_generation =
        gravity::GravitySourceGeneration{m_authoritative_top_domain_source_generation};
  }

  // O(N_local) leaf refit for the force-solve path. Cached-kick paths return
  // before this call, so the scan runs only when a solve will execute. The
  // refit always consumes the stable seed leaves (decomposition-local leaf
  // identities and SFC intervals), never the previously published refit
  // result, so currently-empty seed groups keep their routing partition
  // identity across ordinary drift/source-generation refreshes. Seed leaves
  // retain ownership/epoch; only derived bounds and the source-generation
  // stamp of the published set advance.
  void refreshAuthoritativeTopDomainGeometryForSolve(
      const core::StepContext& context,
      const parallel::MpiContext& mpi_context) {
    m_domain_geometry_refreshed_this_solve = false;
    m_domain_geometry_refreshed_leaf_count = 0U;
    m_domain_geometry_out_of_seed_range_source_count = 0U;
    const std::uint64_t current_source_generation = context.state.gravitySourceGeneration();
    if ((m_authoritative_top_domain_leaves.empty() &&
         m_authoritative_top_domain_seed_leaves.empty()) &&
        !m_authoritative_top_domain_source_generation_valid) {
      publishAuthoritativeTopDomainLeavesSpan();
      return;
    }
    if (!m_authoritative_top_domain_source_generation_valid) {
      publishAuthoritativeTopDomainLeavesSpan();
      return;
    }
    if (m_authoritative_top_domain_source_generation_valid &&
        m_authoritative_top_domain_source_generation == current_source_generation &&
        m_source_predicted_inactive_count == 0U) {
      publishAuthoritativeTopDomainLeavesSpan();
      return;
    }
    parallel::DecompositionConfig refit_config;
    refit_config.world_size = std::max(mpi_context.worldSize(), 1);
    refit_config.domain_x_min_comov = 0.0;
    refit_config.domain_x_max_comov = m_config.cosmology.box_size_x_mpc_comoving;
    refit_config.domain_y_min_comov = 0.0;
    refit_config.domain_y_max_comov = m_config.cosmology.box_size_y_mpc_comoving;
    refit_config.domain_z_min_comov = 0.0;
    refit_config.domain_z_max_comov = m_config.cosmology.box_size_z_mpc_comoving;
    parallel::TopDomainGeometryRefitDiagnostics refit_diagnostics;
    std::exception_ptr refit_failure;
    std::vector<parallel::TopDomainLeaf> refreshed;
    // Ordinary refit seeds from the stable decomposition-local leaf set. The
    // previously published (possibly empty-omitted) result is never the seed;
    // fall back to it only if a legacy state populated current leaves without
    // a recorded seed.
    const std::span<const parallel::TopDomainLeaf> refit_seed =
        !m_authoritative_top_domain_seed_leaves.empty()
            ? std::span<const parallel::TopDomainLeaf>(
                  m_authoritative_top_domain_seed_leaves)
            : std::span<const parallel::TopDomainLeaf>(
                  m_authoritative_top_domain_leaves);
    try {
      refreshed = parallel::refitAuthoritativeTopDomainLeaves(
          refit_seed,
          sourcePosX(context.state),
          sourcePosY(context.state),
          sourcePosZ(context.state),
          refit_config,
          std::max(mpi_context.worldRank(), 0),
          m_decomposition_epoch.value,
          &refit_diagnostics);
    } catch (...) {
      refit_failure = std::current_exception();
    }
    FailureCoordinator(m_services).rethrowCollectiveFailure(
        refit_failure, "gravity authoritative top-domain geometry refit");
    m_authoritative_top_domain_leaves = std::move(refreshed);
    m_authoritative_top_domain_source_generation = current_source_generation;
    m_authoritative_top_domain_source_generation_valid = true;
    m_domain_geometry_refreshed_this_solve = true;
    m_domain_geometry_refreshed_leaf_count = refit_diagnostics.refreshed_leaf_count;
    m_domain_geometry_out_of_seed_range_source_count =
        refit_diagnostics.out_of_seed_range_source_count;
    publishAuthoritativeTopDomainLeavesSpan();
  }

  void commitParticleDecompositionChange() {
    m_decomposition_epoch = gravity::nextGravityIdentity(
        m_decomposition_epoch, "TreePM particle-decomposition epoch overflow");
    // A committed ownership change immediately invalidates routing geometry
    // produced under the previous epoch. Stale seed/current leaves are
    // retained only for diagnostics and are marked unusable (freshness false;
    // publish exposes an empty routing view) until the workflow installs
    // fresh decomposition leaves. Ownership data, migration state, and the
    // decomposition epoch itself are not cleared here.
    m_authoritative_top_domain_source_generation_valid = false;
    // Acceleration lanes are dense-row mirrors and are not part of the
    // migration payload. Invalidate immediately so an output boundary after
    // migration serializes an honest cache.valid=false state instead of
    // rejecting or persisting stale row ownership.
    m_force_cache_valid = false;
    m_force_cache_particle_index_generation = 0U;
    m_force_cache_cell_index_generation = 0U;
    m_force_cache_gas_identity_generation = 0U;
    m_force_cache_source_generation = 0U;
    m_particle_force_cache_valid.clear();
    m_cell_force_cache_valid.clear();
    m_hierarchical_force_generation.clear();
    m_hierarchical_nearest_r2.clear();
    m_hierarchical_pm_source_generation = 0U;
    m_hierarchical_pm_field_version = 0U;
    m_active_source_indices = {};
    m_borrowed_particle_index_generation = 0U;
    m_borrowed_gravity_source_generation = 0U;
    if (m_source_representation ==
        gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      m_source_representation =
          gravity::GravitySourceRepresentation::kMaterializedGeneric;
    }
  }

  void installAuthoritativeTopDomainLeaves(
      std::span<const parallel::TopDomainLeaf> leaves,
      std::uint64_t source_generation) override {
    for (const parallel::TopDomainLeaf& leaf : leaves) {
      if (leaf.owner_rank != m_services.mpi_context.worldRank()) {
        throw std::invalid_argument(
            "gravity authoritative top-domain install received a leaf owned by another rank");
      }
      if (leaf.decomposition_epoch != m_decomposition_epoch.value) {
        throw std::invalid_argument(
            "gravity authoritative top-domain install received a stale decomposition epoch");
      }
    }
    // True ownership/decomposition events (segment start, restart, migration
    // commit) replace both the stable seed set and the current published set
    // with the same fresh decomposition leaves, stamped with the source
    // generation they cover. Ordinary drift/source refreshes never write the
    // seed set.
    m_authoritative_top_domain_seed_leaves.assign(leaves.begin(), leaves.end());
    m_authoritative_top_domain_leaves.assign(leaves.begin(), leaves.end());
    m_authoritative_top_domain_source_generation = source_generation;
    m_authoritative_top_domain_source_generation_valid = true;
    publishAuthoritativeTopDomainLeavesSpan();
  }

  [[nodiscard]] bool authoritativeDomainGeometryMatches(
      std::uint64_t source_generation) const noexcept override {
    return m_authoritative_top_domain_source_generation_valid &&
        m_authoritative_top_domain_source_generation == source_generation;
  }

  [[nodiscard]] std::span<const double> cellAccelX() const noexcept { return m_cell_accel_x; }
  [[nodiscard]] std::span<const double> cellAccelY() const noexcept { return m_cell_accel_y; }
  [[nodiscard]] std::span<const double> cellAccelZ() const noexcept { return m_cell_accel_z; }
  [[nodiscard]] std::span<const double> particleAccelX() const noexcept { return m_particle_accel_x; }
  [[nodiscard]] std::span<const double> particleAccelY() const noexcept { return m_particle_accel_y; }
  [[nodiscard]] std::span<const double> particleAccelZ() const noexcept { return m_particle_accel_z; }

  [[nodiscard]] io::GravityForceCachePersistentView restartForceCacheView(
      const core::SimulationState& state) const {
    io::GravityForceCachePersistentView cache;
    cache.valid = m_force_cache_valid;
    if (!cache.valid) {
      return cache;
    }
    if (m_force_cache_particle_index_generation != state.particleIndexGeneration() ||
        m_force_cache_cell_index_generation != state.cellIndexGeneration() ||
        m_force_cache_gas_identity_generation != state.gasCellIdentityGeneration() ||
        m_force_cache_source_generation != state.gravitySourceGeneration()) {
      throw std::runtime_error(
          "gravity force cache cannot be checkpointed after particle/cell ownership or gas identity changed without refresh");
    }
    if (m_particle_accel_x.size() != state.particles.size() ||
        m_particle_accel_y.size() != state.particles.size() ||
        m_particle_accel_z.size() != state.particles.size() ||
        m_cell_accel_x.size() != state.cells.size() ||
        m_cell_accel_y.size() != state.cells.size() ||
        m_cell_accel_z.size() != state.cells.size()) {
      throw std::runtime_error(
          "gravity force cache cannot be checkpointed because its dense lanes do not match SimulationState");
    }
    cache.particle_id = state.particle_sidecar.particle_id;
    cache.gas_cell_id = state.gas_cells.gas_cell_id;
    cache.particle_accel_x_comoving = m_particle_accel_x;
    cache.particle_accel_y_comoving = m_particle_accel_y;
    cache.particle_accel_z_comoving = m_particle_accel_z;
    cache.cell_accel_x_comoving = m_cell_accel_x;
    cache.cell_accel_y_comoving = m_cell_accel_y;
    cache.cell_accel_z_comoving = m_cell_accel_z;
    // v24: persist the synchronized hierarchical split endpoint whenever it
    // is complete at the checkpoint positions, so a resumed run reuses the
    // exact force the uninterrupted run reuses (no re-evaluation whose
    // relative MAC would see a different force history).
    if (hierarchicalEndpointCompleteAt(state)) {
      cache.hierarchical_endpoint_valid = true;
      cache.hierarchical_pm_field_version = m_hierarchical_pm_field_version;
      cache.hierarchical_tree_accel_x_comoving = m_hierarchical_tree_x;
      cache.hierarchical_tree_accel_y_comoving = m_hierarchical_tree_y;
      cache.hierarchical_tree_accel_z_comoving = m_hierarchical_tree_z;
      cache.hierarchical_pm_accel_x_comoving = m_hierarchical_pm_x;
      cache.hierarchical_pm_accel_y_comoving = m_hierarchical_pm_y;
      cache.hierarchical_pm_accel_z_comoving = m_hierarchical_pm_z;
      cache.hierarchical_nearest_pair_r2_comoving = m_hierarchical_nearest_r2;
    }
    return cache;
  }

  // Local completeness of the synchronized hierarchical endpoint: every row
  // carries Tree and PM components (and the nearest-pair lane) evaluated at
  // the current source positions, with the PM component from one committed
  // field. Only an all-active synchronization endpoint establishes this.
  [[nodiscard]] bool hierarchicalEndpointCompleteAt(const core::SimulationState& state) const noexcept {
    const std::uint64_t generation = state.gravitySourceGeneration();
    return m_hierarchical_pm_source_generation == generation && m_hierarchical_pm_field_version != 0U &&
        hierarchicalTreeEndpointCompleteAt(state);
  }

  // Tree component and total history of every row at the current positions;
  // the PM lanes must exist (cycle-opening field) but need not be current.
  [[nodiscard]] bool hierarchicalTreeEndpointCompleteAt(const core::SimulationState& state) const noexcept {
    const std::size_t n = state.particles.size();
    const std::uint64_t generation = state.gravitySourceGeneration();
    if (!m_force_cache_valid || state.cells.size() != 0U ||
        m_force_cache_particle_index_generation != state.particleIndexGeneration() ||
        m_force_cache_source_generation != generation || m_hierarchical_pm_field_version == 0U ||
        m_hierarchical_force_generation.size() != n || m_particle_force_cache_valid.size() != n ||
        m_particle_accel_x.size() != n || m_particle_accel_y.size() != n || m_particle_accel_z.size() != n ||
        m_hierarchical_tree_x.size() != n || m_hierarchical_tree_y.size() != n ||
        m_hierarchical_tree_z.size() != n || m_hierarchical_pm_x.size() != n ||
        m_hierarchical_pm_y.size() != n || m_hierarchical_pm_z.size() != n ||
        m_hierarchical_nearest_r2.size() != n) {
      return false;
    }
    for (std::size_t row = 0; row < n; ++row) {
      if (m_hierarchical_force_generation[row] != generation || m_particle_force_cache_valid[row] == 0U) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] bool hierarchicalEndpointForceReusable(
      const core::SimulationState& state, const core::IntegratorState& integrator) const override {
    return integrator.pm_long_range_field_valid &&
        m_hierarchical_pm_field_version == integrator.pm_sync_state.fieldVersion() &&
        hierarchicalEndpointCompleteAt(state);
  }

  [[nodiscard]] bool hierarchicalShortRangeEndpointReusable(
      const core::SimulationState& state) const override {
    return hierarchicalTreeEndpointCompleteAt(state);
  }

  [[nodiscard]] std::span<const double> hierarchicalNearestPairR2() const noexcept override {
    return m_hierarchical_nearest_r2;
  }

  // Ordinary restart import: the checkpoint may come from a different row
  // order, so particle rows are matched by stable ID
  // (internal::remapForceEndpointById). Transactional rollback inside a
  // frozen-layout PM cycle uses restoreTransactionalForceEndpoint instead.
  void importRestartForceCache(
      const io::GravityForceCachePersistentState& cache,
      const core::SimulationState& state) {
    if (!cache.valid) {
      m_particle_accel_x.clear();
      m_particle_accel_y.clear();
      m_particle_accel_z.clear();
      m_cell_accel_x.clear();
      m_cell_accel_y.clear();
      m_cell_accel_z.clear();
      m_particle_force_cache_valid.clear();
      m_cell_force_cache_valid.clear();
      m_force_cache_valid = false;
      m_force_cache_particle_index_generation = 0U;
      m_force_cache_cell_index_generation = 0U;
      m_force_cache_gas_identity_generation = 0U;
      m_force_cache_source_generation = 0U;
      return;
    }
    const std::size_t n = state.particles.size();
    const std::size_t cell_count = state.cells.size();
    if (cache.cell_accel_x_comoving.size() != cell_count || cache.cell_accel_y_comoving.size() != cell_count ||
        cache.cell_accel_z_comoving.size() != cell_count) {
      throw std::runtime_error(
          "restart gravity force cache gas-cell extent does not match authoritative SimulationState");
    }
    for (const auto* lane : {&cache.cell_accel_x_comoving, &cache.cell_accel_y_comoving, &cache.cell_accel_z_comoving}) {
      for (const double value : *lane) {
        if (!std::isfinite(value)) {
          throw std::runtime_error("restart gravity force cache contains non-finite values");
        }
      }
    }
    if (cache.gas_cell_id.size() != cell_count) {
      throw std::runtime_error(
          "restart gravity force cache gas-cell identity extent does not match acceleration lanes");
    }
    std::unordered_map<std::uint64_t, std::size_t> gas_index_by_id;
    gas_index_by_id.reserve(cell_count);
    for (std::size_t index = 0; index < cell_count; ++index) {
      if (cache.gas_cell_id[index] == 0U || !gas_index_by_id.emplace(cache.gas_cell_id[index], index).second) {
        throw std::runtime_error(
            "restart gravity force cache gas-cell identity lane contains a zero or duplicate stable ID");
      }
    }
    // v24 hierarchical endpoint. Absent (pre-v24 or not at a synchronized
    // hierarchical endpoint) leaves the split lanes unmaterialized, so the
    // hierarchical segment performs its synchronized bootstrap evaluation.
    const bool endpoint = cache.hierarchical_endpoint_valid;
    m_particle_accel_x.resize(n);
    m_particle_accel_y.resize(n);
    m_particle_accel_z.resize(n);
    if (endpoint) {
      for (auto* lane : {&m_hierarchical_tree_x, &m_hierarchical_tree_y, &m_hierarchical_tree_z,
                         &m_hierarchical_pm_x, &m_hierarchical_pm_y, &m_hierarchical_pm_z,
                         &m_hierarchical_nearest_r2}) {
        lane->resize(n);
      }
    }
    internal::remapForceEndpointById(cache, state.particle_sidecar.particle_id,
        endpoint ? forceEndpointLanes() : internal::ForceEndpointLanes{
            .total_x = m_particle_accel_x, .total_y = m_particle_accel_y, .total_z = m_particle_accel_z},
        endpoint);
    m_cell_accel_x.resize(cell_count);
    m_cell_accel_y.resize(cell_count);
    m_cell_accel_z.resize(cell_count);
    for (std::size_t row = 0; row < cell_count; ++row) {
      const auto it = gas_index_by_id.find(state.gas_cells.gas_cell_id[row]);
      if (it == gas_index_by_id.end()) {
        throw std::runtime_error("restart gravity force cache is missing a persistent gas-cell ID");
      }
      const std::size_t cached_row = it->second;
      m_cell_accel_x[row] = cache.cell_accel_x_comoving[cached_row];
      m_cell_accel_y[row] = cache.cell_accel_y_comoving[cached_row];
      m_cell_accel_z[row] = cache.cell_accel_z_comoving[cached_row];
    }
    m_force_cache_valid = true;
    m_force_cache_particle_index_generation = state.particleIndexGeneration();
    m_force_cache_cell_index_generation = state.cellIndexGeneration();
    m_force_cache_gas_identity_generation = state.gasCellIdentityGeneration();
    m_force_cache_source_generation = state.gravitySourceGeneration();
    m_particle_force_cache_valid.assign(n, 1U);
    m_cell_force_cache_valid.assign(cell_count, 1U);
    m_hierarchical_force_generation.clear();
    m_hierarchical_pm_source_generation = 0U;
    m_hierarchical_pm_field_version = 0U;
    if (!endpoint) {
      m_hierarchical_nearest_r2.clear();
      return;
    }
    m_hierarchical_force_generation.assign(n, state.gravitySourceGeneration());
    m_hierarchical_pm_source_generation = state.gravitySourceGeneration();
    m_hierarchical_pm_field_version = cache.hierarchical_pm_field_version;
  }

  // Runtime-owned dense endpoint lanes (current extents).
  [[nodiscard]] internal::ForceEndpointLanes forceEndpointLanes() noexcept {
    return internal::ForceEndpointLanes{
        .total_x = m_particle_accel_x, .total_y = m_particle_accel_y, .total_z = m_particle_accel_z,
        .tree_x = m_hierarchical_tree_x, .tree_y = m_hierarchical_tree_y, .tree_z = m_hierarchical_tree_z,
        .pm_x = m_hierarchical_pm_x, .pm_y = m_hierarchical_pm_y, .pm_z = m_hierarchical_pm_z,
        .nearest_pair_r2 = m_hierarchical_nearest_r2};
  }

  void captureTransactionalForceEndpoint(
      const core::SimulationState& state, std::uint64_t layout_certificate,
      TransactionalForceEndpoint& endpoint) const override {
    if (!hierarchicalEndpointCompleteAt(state)) {
      endpoint.valid = false;
      throw std::logic_error("multirate PM cycle opening lacks a complete synchronized force endpoint");
    }
    internal::captureForceEndpointRowAligned(
        internal::ConstForceEndpointLanes{
            .total_x = m_particle_accel_x, .total_y = m_particle_accel_y, .total_z = m_particle_accel_z,
            .tree_x = m_hierarchical_tree_x, .tree_y = m_hierarchical_tree_y, .tree_z = m_hierarchical_tree_z,
            .pm_x = m_hierarchical_pm_x, .pm_y = m_hierarchical_pm_y, .pm_z = m_hierarchical_pm_z,
            .nearest_pair_r2 = m_hierarchical_nearest_r2},
        internal::RowAlignedEndpointStamp{
            .particle_count = static_cast<std::uint64_t>(state.particles.size()),
            .particle_index_generation = state.particleIndexGeneration(),
            .layout_certificate = layout_certificate},
        m_hierarchical_pm_field_version, endpoint);
  }

  // Transactional rollback restore (frozen layout): row-aligned copy of the
  // cycle-opening endpoint, then the force-history metadata re-stamped at
  // the restored (advanced) source generation. Every precondition, including
  // the extents of the per-row metadata lanes, is checked before the first
  // write, so a refusal changes nothing. No ID map, no allocation.
  void validateTransactionalForceEndpoint(
      const TransactionalForceEndpoint& endpoint, const core::SimulationState& state,
      std::uint64_t layout_certificate) const override {
    const std::size_t n = state.particles.size();
    if (state.cells.size() != 0U || !m_cell_accel_x.empty() || !m_cell_force_cache_valid.empty()) {
      throw std::logic_error("transactional force endpoint restore refused: gas cells present (hierarchical KDK is DMO-only)");
    }
    if (m_particle_force_cache_valid.size() != n || m_hierarchical_force_generation.size() != n) {
      throw std::logic_error("transactional force endpoint restore refused: per-row force metadata does not cover every particle row");
    }
    internal::validateRowAlignedForceEndpoint(
        endpoint,
        internal::RowAlignedEndpointStamp{
            .particle_count = static_cast<std::uint64_t>(n),
            .particle_index_generation = state.particleIndexGeneration(),
            .layout_certificate = layout_certificate},
        internal::ConstForceEndpointLanes{
            .total_x = m_particle_accel_x, .total_y = m_particle_accel_y, .total_z = m_particle_accel_z,
            .tree_x = m_hierarchical_tree_x, .tree_y = m_hierarchical_tree_y, .tree_z = m_hierarchical_tree_z,
            .pm_x = m_hierarchical_pm_x, .pm_y = m_hierarchical_pm_y, .pm_z = m_hierarchical_pm_z,
            .nearest_pair_r2 = m_hierarchical_nearest_r2});
  }

  void restoreTransactionalForceEndpoint(
      const TransactionalForceEndpoint& endpoint, const core::SimulationState& state,
      std::uint64_t layout_certificate) {
    validateTransactionalForceEndpoint(endpoint, state, layout_certificate);
    internal::restoreForceEndpointRowAligned(
        endpoint,
        internal::RowAlignedEndpointStamp{
            .particle_count = static_cast<std::uint64_t>(state.particles.size()),
            .particle_index_generation = state.particleIndexGeneration(),
            .layout_certificate = layout_certificate},
        forceEndpointLanes());
    const std::uint64_t generation = state.gravitySourceGeneration();
    std::fill(m_particle_force_cache_valid.begin(), m_particle_force_cache_valid.end(), std::uint8_t{1U});
    std::fill(m_hierarchical_force_generation.begin(), m_hierarchical_force_generation.end(), generation);
    m_force_cache_valid = true;
    m_force_cache_particle_index_generation = state.particleIndexGeneration();
    m_force_cache_cell_index_generation = state.cellIndexGeneration();
    m_force_cache_gas_identity_generation = state.gasCellIdentityGeneration();
    m_force_cache_source_generation = generation;
    m_hierarchical_pm_source_generation = generation;
    m_hierarchical_pm_field_version = endpoint.pm_field_version;
  }

  [[nodiscard]] static PmSyncSurface toPmSyncSurface(core::IntegrationStage stage) {
    if (stage == core::IntegrationStage::kGravityKickPre) {
      return PmSyncSurface::kKickPre;
    }
    if (stage == core::IntegrationStage::kForceRefresh) {
      return PmSyncSurface::kForceRefresh;
    }
    if (stage == core::IntegrationStage::kGravityKickPost) {
      return PmSyncSurface::kKickPost;
    }
    throw std::invalid_argument("TreePM sync surface is only defined on legal gravity force-refresh/kick stages");
  }

  [[nodiscard]] static std::string_view pmSyncSurfaceName(PmSyncSurface surface) {
    switch (surface) {
      case PmSyncSurface::kKickPre:
        return "kick_pre";
      case PmSyncSurface::kForceRefresh:
        return "force_refresh";
      case PmSyncSurface::kKickPost:
        return "kick_post";
    }
    return "unknown";
  }
  [[nodiscard]] static std::string_view pmRefreshReasonName(core::PmRefreshDirective::Reason reason) {
    switch (reason) {
      case core::PmRefreshDirective::Reason::kNone:
        return "none";
      case core::PmRefreshDirective::Reason::kInitialForceBootstrap:
        return "initial_force_bootstrap";
      case core::PmRefreshDirective::Reason::kScheduledForceRefreshStage:
        return "scheduled_force_refresh_stage";
      case core::PmRefreshDirective::Reason::kSourceMutationForceRefresh:
        return "source_mutation_force_refresh";
    }
    return "unknown";
  }

  [[nodiscard]] RuntimeTaskMemoryEstimate estimateMemory(
      const core::SimulationState& state) const override {
    const auto bytes = estimateIncrementalMemory(state);
    // This estimate runs at the serial task boundary, before dispatcher and
    // owner admission. Optional committed caches must not block mandatory RAM.
    m_tree_pm_coordinator.trimOptionalAccelerationCaches(bytes);
    return {bytes, true, {}};
  }

  [[nodiscard]] std::uint64_t estimateIncrementalMemory(
      const core::SimulationState& state) const {
    const auto representation = borrowedHomogeneousDmoBaseEligible(state)
        ? gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
        : gravity::GravitySourceRepresentation::kMaterializedGeneric;
    return estimateIncrementalMemoryForRepresentation(state, representation);
  }

  [[nodiscard]] std::uint64_t estimateIncrementalMemoryForRepresentation(
      const core::SimulationState& state,
      gravity::GravitySourceRepresentation representation) const {
    const std::size_t particle_count = state.particles.size();
    const std::size_t cell_count = state.cells.size();
    const parallel::MpiContext& mpi_context = m_services.mpi_context;
    const std::uint64_t gravity_baseline_before =
        core::memoryReportBaselineOwnedBytes(memoryReport());
    const std::size_t conservative_count = core::checkedSizeAdd(
            particle_count, cell_count,
            "gravity governed conservative source/target count");
    const gravity::GravityMemoryEstimate governed_peak =
            gravity::estimateGravityMemory(gravity::GravityMemoryEstimateInput{
                .local_source_count = static_cast<std::uint64_t>(conservative_count),
                .local_target_count = static_cast<std::uint64_t>(conservative_count),
                .local_particle_count = static_cast<std::uint64_t>(particle_count),
                .local_cell_count = static_cast<std::uint64_t>(cell_count),
                .source_representation = representation,
                .source_softening_uniform =
                    representation == gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo,
                .periodic_fft_backed_density =
                    m_tree_pm_options.pm_options.boundary_condition ==
                        gravity::PmBoundaryCondition::kPeriodic &&
                    !m_runtime_topology.usesCuda(),
                .pm_force_gradient_backend = m_tree_pm_options.pm_options.force_gradient_backend,
                .relative_force_mac_enabled =
                    (m_tree_pm_options.acceptance_policy == gravity::TreePmAcceptancePolicy::kAdaptiveRelative ||
                 m_tree_pm_options.tree_options.opening_criterion ==
                    gravity::TreeOpeningCriterion::kRelativeForceError),
                .hierarchical_kdk_enabled = m_config.numerics.hierarchical_max_rung > 0,
                .gaussian_pair_lookup_enabled = m_config.numerics.treepm_gaussian_pair_lookup_enabled,
                .residual_block_size = m_tree_pm_options.residual_block_size,
                .tree_leaf_size = m_tree_pm_options.tree_options.max_leaf_size,
                .multipole_order = m_tree_pm_options.tree_options.multipole_order,
                .pm_shape = m_pm_grid_shape,
                .assignment_scheme = m_tree_pm_options.pm_options.assignment_scheme,
                .decomposition_mode = m_tree_pm_options.pm_options.decomposition_mode,
                .mpi_rank_count = static_cast<std::uint32_t>(
                    std::max(m_runtime_topology.world_size, 1)),
                .mpi_world_rank = std::max(mpi_context.worldRank(), 0),
                .zoom_enabled = m_tree_pm_options.enable_zoom_long_range_correction,
                .zoom_pm_shape = m_tree_pm_options.zoom_focused_pm_shape,
                .periodic_tree_coordinates =
                    m_tree_pm_options.pm_options.boundary_condition ==
                    gravity::PmBoundaryCondition::kPeriodic,
                .cuda_resident = m_runtime_topology.usesCuda(),
                .tree_exchange_batch_bytes =
                    m_tree_pm_options.tree_exchange_batch_bytes,
                .pm_exchange_batch_bytes =
                    m_tree_pm_options.pm_options.routing_exchange_batch_bytes,
                .backend_unknown_reserve_bytes = 0U,
                .safety_margin_fraction = 0.0,
            });
    // The production communication arena is a separately governed physical
    // owner. Exclude it from the phase-resident reservation so the same bytes
    // are not reserved twice by the MemoryGovernor.
    if (governed_peak.communication_arena_bytes > governed_peak.known_peak_bytes) {
      throw std::logic_error(
          "gravity communication arena estimate exceeds the total governed peak");
    }
    std::uint64_t phase_resident_peak =
        governed_peak.known_peak_bytes - governed_peak.communication_arena_bytes;
    if (m_config.numerics.treepm_gaussian_pair_lookup_enabled) {
      // The coefficient table has its own persistent commitment. Keep it in
      // whole-runtime preflight, exclude it from this overlapping phase lease.
      const auto table_bytes = sizeof(gravity::TreePmGaussianCoefficientTable);
      if (table_bytes > phase_resident_peak) throw std::logic_error("Gaussian table exceeds modeled gravity peak");
      phase_resident_peak -= table_bytes;
    }
    if (representation == gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      // The sole uint32 identity lane is physically owned and admitted by the
      // governed TransientStepWorkspace scratch arena. The gravity estimate
      // keeps it in process/preflight modeling, but the phase lease must not
      // reserve the same physical bytes a second time.
      const std::uint64_t scratch_index_bytes = static_cast<std::uint64_t>(
          core::checkedSizeMultiply(
              particle_count, sizeof(std::uint32_t),
              "borrowed DMO target/source scratch byte count"));
      if (scratch_index_bytes > phase_resident_peak) {
        throw std::logic_error(
            "borrowed DMO scratch index estimate exceeds gravity phase peak");
      }
      phase_resident_peak -= scratch_index_bytes;
    }
    // Recent cadence has bounded control-path strings and records. Model its
    // full warm capacity before any force event grows it; observed capacities
    // are reconciled through memoryReport(). Opaque allocator costs remain in
    // the existing process allocator reserve.
    phase_resident_peak = core::checkedMemoryBytesAdd(phase_resident_peak,
        256U * (sizeof(ReferenceWorkflowReport::TreePmCadenceRecord) + 256U),
        "gravity bounded cadence phase model");
    return incrementalMemoryBeyondRetained(
        phase_resident_peak, gravity_baseline_before);
  }

  void execute(GravityStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::gravityContext(
        view,
        {{RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleGravitySource, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroPrimitiveState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kSchedulerTruth, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kGravityAcceleration, RuntimeResourceAccessMode::kWrite},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kReadWrite}});
    const bool hierarchical = context.hierarchical_kdk.enabled;
    const bool is_kick_stage = context.stage == core::IntegrationStage::kGravityKickPre ||
        context.stage == core::IntegrationStage::kGravityKickPost;
    const parallel::MpiContext& mpi_context = m_services.mpi_context;
    const std::uint64_t world_size = static_cast<std::uint64_t>(mpi_context.worldSize());
    const std::size_t particle_count = context.state.particles.size();
    const std::size_t cell_count = context.state.cells.size();

    if (hierarchical && (!borrowedHomogeneousDmoBaseEligible(context.state) || cell_count != 0U ||
        context.particle_scheduler == nullptr)) {
      throw std::logic_error("hierarchical gravity requires borrowed homogeneous DMO and scheduler authority");
    }
    if (hierarchical && is_kick_stage && !context.pm_refresh_directive.force_refresh_surface) {
      applyHierarchicalKick(context);
      return;
    }

    // Admit the conservative non-arena gravity high-water before any O(N)
    // cache/vector growth in this phase. The shared communication arena owns
    // and commits its separate Communication reservation before first use; this
    // phase reservation covers the remaining incremental TreePM footprint over
    // the retained gravity baseline. Every rank coordinates admission before
    // later solver collectives, so a low-headroom rank cannot strand peers.
    const std::uint64_t gravity_baseline_before =
        core::memoryReportBaselineOwnedBytes(memoryReport());
    std::uint64_t process_baseline_before = 0U;
    core::MemoryReservation gravity_phase_reservation;
    core::MemoryReservation gravity_fallback_reservation;
    std::exception_ptr gravity_admission_failure;
    try {
      if (m_services.memory_governor != nullptr) {
        process_baseline_before =
            m_services.memory_governor->snapshot().baseline_owned_bytes;
        if (gravity_baseline_before > process_baseline_before) {
          throw std::logic_error(
              "gravity governor baseline is stale: retained gravity capacity exceeds process baseline");
        }
        const std::uint64_t incremental_bytes = estimateIncrementalMemory(context.state);
        gravity_phase_reservation = m_services.memory_governor->reserve(
            core::MemoryClass::kPhaseResident, incremental_bytes,
            "gravity.treepm.phase_peak");
        gravity_phase_reservation.commit();
      }
    } catch (...) {
      gravity_admission_failure = std::current_exception();
    }
    FailureCoordinator(m_services).rethrowCollectiveFailure(
        gravity_admission_failure, "gravity TreePM memory admission");

    const auto admitAllParticleIndexScratch = [&](bool force_all_particles = false) {
      if (context.workspace == nullptr ||
          (!force_all_particles &&
           (context.stage != core::IntegrationStage::kGravityKickPre ||
            context.active_set.hasParticleSubset(particle_count)))) {
        return;
      }
      std::exception_ptr scratch_admission_failure;
      try {
        if (context.workspace->gravity_particle_index_scratch.size() != particle_count) {
          context.workspace->prepareGravityParticleIndexScratch(particle_count);
        }
      } catch (...) {
        scratch_admission_failure = std::current_exception();
      }
      FailureCoordinator(m_services).rethrowCollectiveFailure(
          scratch_admission_failure,
          "gravity all-particle scratch admission");
      auto scratch = context.workspace->gravity_particle_index_scratch;
      if (scratch.size() != particle_count) {
        throw std::logic_error(
            "gravity all-particle scratch admission returned an unexpected extent");
      }
      for (std::size_t particle_index = 0; particle_index < particle_count;
           ++particle_index) {
        if (particle_index > std::numeric_limits<std::uint32_t>::max()) {
          throw std::overflow_error(
              "gravity all-particle scratch exceeds the uint32 local-index contract");
        }
        scratch[particle_index] = static_cast<std::uint32_t>(particle_index);
      }
    };

    const bool particle_cache_generation_changed =
        m_force_cache_particle_index_generation != context.state.particleIndexGeneration() ||
        m_particle_force_cache_valid.size() != particle_count;
    const bool cell_cache_generation_changed =
        m_force_cache_cell_index_generation != context.state.cellIndexGeneration() ||
        m_force_cache_gas_identity_generation != context.state.gasCellIdentityGeneration() ||
        m_cell_force_cache_valid.size() != cell_count;
    const bool source_cache_generation_changed =
        m_force_cache_source_generation != context.state.gravitySourceGeneration();
    // A previous force is an accuracy scale, never a valid kick cache at the
    // new epoch. Admit the known single drift transition only on compact DMO
    // at force_refresh; population/ownership/gas changes and skipped source
    // generations fall back geometrically. Preserve no second force vector.
    const bool adaptive_drift_history_compatible =
        m_tree_pm_options.acceptance_policy == gravity::TreePmAcceptancePolicy::kAdaptiveRelative &&
        m_force_cache_valid && !particle_cache_generation_changed &&
        !cell_cache_generation_changed && context.state.hasHomogeneousDmoMetadata() &&
        m_source_representation == gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo &&
        context.stage == core::IntegrationStage::kForceRefresh &&
        m_force_cache_source_generation > 0U &&
        m_force_cache_source_generation < std::numeric_limits<std::uint64_t>::max() &&
        context.state.gravitySourceGeneration() == m_force_cache_source_generation + 1U;

    if (particle_cache_generation_changed) {
      // Dense-row identity is not a migration payload. A row is compatible only
      // after its authoritative source was refreshed in the current generation.
      m_particle_force_cache_valid.assign(particle_count, 0U);
      m_force_cache_particle_index_generation =
          context.state.particleIndexGeneration();
      m_force_cache_valid = false;
    }
    if (cell_cache_generation_changed) {
      // AMR refinement/coarsening and gas-cell migration change both dense row
      // indices and stable gas ownership. Invalidate together so a post-AMR
      // kick cannot reuse accelerations tied to a superseded leaf population.
      m_cell_force_cache_valid.assign(cell_count, 0U);
      m_force_cache_cell_index_generation = context.state.cellIndexGeneration();
      m_force_cache_gas_identity_generation = context.state.gasCellIdentityGeneration();
      m_force_cache_valid = false;
    }
    if (source_cache_generation_changed && !hierarchical) {
      if (!adaptive_drift_history_compatible) {
        m_particle_force_cache_valid.assign(particle_count, 0U);
        m_cell_force_cache_valid.assign(cell_count, 0U);
      }
      m_force_cache_valid = false;
    }
    const bool local_force_cache_incompatible = !m_force_cache_valid ||
        m_force_cache_particle_index_generation != context.state.particleIndexGeneration() ||
        m_force_cache_cell_index_generation != context.state.cellIndexGeneration() ||
        m_force_cache_gas_identity_generation != context.state.gasCellIdentityGeneration() ||
        m_force_cache_source_generation != context.state.gravitySourceGeneration();
    // One collective carries the cache-compatibility vote (min over ranks of
    // "compatible") together with the stage-entry agreement checks that the
    // force-refresh path previously issued as four separate min/max pairs.
    RankConsensusBatch entry_consensus;
    const std::size_t cache_compatible_lane = entry_consensus.globalMinimum(
        "force_cache_compatible", local_force_cache_incompatible ? 0ULL : 1ULL);
    entry_consensus.requireEqual(
        "pm_sync_event_presence", context.pm_refresh_directive.has_sync_event ? 1ULL : 0ULL);
    entry_consensus.requireEqual("gravity_source_generation", context.state.gravitySourceGeneration());
    entry_consensus.verify(mpi_context, world_size);
    const bool needs_force_cache_rebuild = entry_consensus.minimum(cache_compatible_lane) == 0U;
    const bool needs_initial_force_bootstrap =
        context.stage == core::IntegrationStage::kGravityKickPre &&
        !context.integrator_state.pm_long_range_field_valid;
    if (needs_initial_force_bootstrap && !context.pm_refresh_directive.initial_cache_bootstrap_allowed) {
      throw std::runtime_error(
          "TreePM initial force bootstrap requested outside an integrator-authorized global boundary");
    }
    const bool authoritative_pm_refresh =
        context.pm_refresh_directive.has_sync_event &&
        context.pm_refresh_directive.refresh_long_range_field;
    const bool is_force_refresh_stage = context.pm_refresh_directive.force_refresh_surface ||
        needs_initial_force_bootstrap || authoritative_pm_refresh ||
        (is_kick_stage && needs_force_cache_rebuild);
    if (!is_kick_stage && !is_force_refresh_stage) {
      throw std::logic_error("gravity handler received an unregistered stage");
    }
    const bool initial_bootstrap_boundary_exception =
        context.pm_refresh_directive.reason ==
            core::PmRefreshDirective::Reason::kInitialForceBootstrap &&
        context.pm_refresh_directive.initial_cache_bootstrap_allowed;
    if (context.pm_refresh_directive.has_sync_event &&
        context.pm_refresh_directive.refresh_long_range_field &&
        !context.boundary.pm_refresh_allowed &&
        !initial_bootstrap_boundary_exception) {
      throw std::runtime_error("TreePM long-range PM refresh reached an illegal integration boundary");
    }
    if (context.stage == core::IntegrationStage::kForceRefresh && !context.pm_refresh_directive.force_refresh_surface) {
      throw std::runtime_error("TreePM force-refresh stage lacks an integrator-issued PM refresh directive");
    }
    if (is_kick_stage && !is_force_refresh_stage) {
      // The drift view is constructed only after the pre-kick handler returns.
      // Admit its one governed all-particle index owner here even on a pure
      // cache hit, where the full force-refresh process preflight is skipped.
      admitAllParticleIndexScratch();
      applyCachedKick(context);
      return;
    }

    if (m_runtime_topology.world_size != mpi_context.worldSize()) {
      throw std::runtime_error("reference workflow runtime topology world_size drifted from MPI context");
    }

    // pm_sync_event_presence and gravity_source_generation were verified by
    // entry_consensus above; the PM decision lanes are verified in one batch
    // once the local decision is complete (pm_decision_consensus below).
    RankConsensusBatch pm_decision_consensus;

    const internal::SolverGhostRefreshReport gravity_ghost_refresh =
        internal::refreshParticleGhostsForSolver(
        context, mpi_context, "gravity.treepm", &m_ghost_cache_lifecycle);

    std::span<const std::uint32_t> force_target_particles =
        context.active_set.particle_indices;
    const bool borrowed_dmo_base_eligible =
        borrowedHomogeneousDmoBaseEligible(context.state);
    bool force_targets_use_identity_scratch =
        context.workspace != nullptr &&
        borrowed_dmo_base_eligible &&
        isIdentityParticleSpan(force_target_particles, particle_count);
    SourcePredictionEpoch source_prediction_epoch =
        context.pm_refresh_directive.requires_predicted_inactive_sources
        ? SourcePredictionEpoch::kStepEnd
        : SourcePredictionEpoch::kNone;
    if (is_kick_stage && needs_force_cache_rebuild) {
      // A dense-row ownership transition invalidates acceleration rows by
      // identity, not merely by extent. Rebuild every locally owned row in one
      // coordinated solve so later local pre-kicks never mix stale cache rows
      // with current ones. Sources whose persistent coordinates lag the
      // current scheduler time are predicted to the pre-kick epoch.
      if (borrowed_dmo_base_eligible && context.workspace != nullptr) {
        releaseVectorCapacity(m_force_refresh_particle_indices);
        force_target_particles = {};
        force_targets_use_identity_scratch = true;
      } else {
        m_force_refresh_particle_indices.clear();
        m_force_refresh_particle_indices.reserve(particle_count);
        const std::uint32_t local_rank =
            static_cast<std::uint32_t>(std::max(mpi_context.worldRank(), 0));
        for (std::size_t particle_index = 0; particle_index < particle_count;
             ++particle_index) {
          if (context.state.particleOwningRank(particle_index) ==
              local_rank) {
            m_force_refresh_particle_indices.push_back(
                static_cast<std::uint32_t>(particle_index));
          }
        }
        force_target_particles = m_force_refresh_particle_indices;
      }
      // A pre-kick rebuild evaluates the force cache at the step-begin
      // epoch. A post-kick rebuild evaluates at step end, where source-term
      // births already exist. Treating both kick stages as step-begin would
      // attempt to back-predict newborn sources to an epoch before creation.
      source_prediction_epoch =
          context.stage == core::IntegrationStage::kGravityKickPre
          ? SourcePredictionEpoch::kStepBegin
          : SourcePredictionEpoch::kStepEnd;
    }

    // Budget against a conservative per-rank upper bound before materializing
    // O(N) gravity staging/tree/communication buffers. Canonical state is
    // deliberately outside this gravity-only budget and is reported by core.
    const std::size_t conservative_local_source_count = core::checkedSizeAdd(
        particle_count, cell_count, "gravity pre-run local source count");
    const std::size_t force_target_count = force_targets_use_identity_scratch
        ? particle_count
        : force_target_particles.size();
    const bool prediction_requires_inactive_coordinates =
        source_prediction_epoch != SourcePredictionEpoch::kNone &&
        !isIdentityParticleSpan(context.active_set.particle_indices, particle_count);
    const bool borrowed_homogeneous_dmo_preflight =
        borrowed_dmo_base_eligible && (force_targets_use_identity_scratch || hierarchical) &&
        !prediction_requires_inactive_coordinates;

    // The early phase reservation is allowed to credit the state-only DMO
    // opportunity. If the stage shape proves that the borrowed identity view
    // cannot be used, collectively re-admit the generic increment before any
    // materialized O(N) source/map allocation occurs.
    if (borrowed_dmo_base_eligible &&
        !borrowed_homogeneous_dmo_preflight &&
        m_services.memory_governor != nullptr) {
      std::exception_ptr fallback_admission_failure;
      try {
        const std::uint64_t borrowed_increment =
            estimateIncrementalMemoryForRepresentation(
                context.state,
                gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo);
        const std::uint64_t generic_increment =
            estimateIncrementalMemoryForRepresentation(
                context.state,
                gravity::GravitySourceRepresentation::kMaterializedGeneric);
        const std::uint64_t additional_bytes = generic_increment > borrowed_increment
            ? generic_increment - borrowed_increment
            : 0U;
        gravity_fallback_reservation = m_services.memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            additional_bytes,
            "gravity.treepm.generic_fallback");
        gravity_fallback_reservation.commit();
      } catch (...) {
        fallback_admission_failure = std::current_exception();
      }
      FailureCoordinator(m_services).rethrowCollectiveFailure(
          fallback_admission_failure,
          "gravity TreePM generic fallback memory admission");
    }

    const std::size_t conservative_local_target_count = core::checkedSizeAdd(
        force_target_count, cell_count, "gravity pre-run local target count");
    m_last_pre_run_memory_estimate = gravity::estimateGravityMemory(
        gravity::GravityMemoryEstimateInput{
            .local_source_count = static_cast<std::uint64_t>(conservative_local_source_count),
            .local_target_count = static_cast<std::uint64_t>(conservative_local_target_count),
            .local_particle_count = static_cast<std::uint64_t>(particle_count),
            .local_cell_count = static_cast<std::uint64_t>(cell_count),
            .source_representation = borrowed_homogeneous_dmo_preflight
                ? gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
                : gravity::GravitySourceRepresentation::kMaterializedGeneric,
            .source_softening_uniform = borrowed_homogeneous_dmo_preflight,
            .periodic_fft_backed_density =
                m_tree_pm_options.pm_options.boundary_condition ==
                    gravity::PmBoundaryCondition::kPeriodic &&
                !m_runtime_topology.usesCuda(),
            .pm_force_gradient_backend = m_tree_pm_options.pm_options.force_gradient_backend,
            .relative_force_mac_enabled =
                (m_tree_pm_options.acceptance_policy == gravity::TreePmAcceptancePolicy::kAdaptiveRelative ||
             m_tree_pm_options.tree_options.opening_criterion ==
                gravity::TreeOpeningCriterion::kRelativeForceError),
            .hierarchical_kdk_enabled = m_config.numerics.hierarchical_max_rung > 0,
            .gaussian_pair_lookup_enabled = m_config.numerics.treepm_gaussian_pair_lookup_enabled,
            .residual_block_size = m_tree_pm_options.residual_block_size,
            .tree_leaf_size = m_tree_pm_options.tree_options.max_leaf_size,
            .multipole_order = m_tree_pm_options.tree_options.multipole_order,
            .pm_shape = m_pm_grid_shape,
            .assignment_scheme = m_tree_pm_options.pm_options.assignment_scheme,
            .decomposition_mode = m_tree_pm_options.pm_options.decomposition_mode,
            .mpi_rank_count = static_cast<std::uint32_t>(std::max(m_runtime_topology.world_size, 1)),
            .mpi_world_rank = std::max(mpi_context.worldRank(), 0),
            .zoom_enabled = m_tree_pm_options.enable_zoom_long_range_correction,
            .zoom_pm_shape = m_tree_pm_options.zoom_focused_pm_shape,
            .periodic_tree_coordinates =
                m_tree_pm_options.pm_options.boundary_condition == gravity::PmBoundaryCondition::kPeriodic,
            .cuda_resident = m_runtime_topology.usesCuda(),
            .tree_exchange_batch_bytes = m_tree_pm_options.tree_exchange_batch_bytes,
            .pm_exchange_batch_bytes = m_tree_pm_options.pm_options.routing_exchange_batch_bytes,
            .backend_unknown_reserve_bytes =
                m_config.parallel.gravity_backend_unknown_reserve_bytes,
            .safety_margin_fraction =
                m_config.parallel.gravity_memory_safety_margin_fraction,
        });
    m_has_pre_run_memory_estimate = true;
    gravity::enforceGravityMemoryBudget(
        m_last_pre_run_memory_estimate,
        m_config.parallel.gravity_memory_budget_bytes);

    if (m_config.parallel.process_memory_budget_bytes != 0U &&
        context.particle_scheduler == nullptr) {
      throw std::runtime_error(
          "DMO process memory budget is enabled but the gravity stage has no authoritative particle scheduler for owned-memory preflight");
    }
    const std::uint64_t scheduler_current_size_bytes = context.particle_scheduler != nullptr
        ? context.particle_scheduler->logicalSizeBytes()
        : 0U;
    const std::uint64_t scheduler_owned_bytes = context.particle_scheduler != nullptr
        ? context.particle_scheduler->ownedCapacityBytes()
        : 0U;
    const std::uint64_t scheduler_high_water_bytes = context.particle_scheduler != nullptr
        ? context.particle_scheduler->ownedCapacityHighWaterBytes()
        : 0U;
    const core::MemoryReport canonical_runtime_report =
        core::collectSimulationMemoryReport(context.state, context.workspace);
    m_last_process_memory_estimate = gravity::estimateDmoProcessMemory(
        canonical_runtime_report,
        m_last_pre_run_memory_estimate,
        gravity::DmoProcessMemoryPolicy{
            .mpi_rank_count = static_cast<std::uint32_t>(
                std::max(m_runtime_topology.world_size, 1)),
            .scheduler_current_size_bytes = scheduler_current_size_bytes,
            .scheduler_owned_bytes = scheduler_owned_bytes,
            .scheduler_high_water_bytes = scheduler_high_water_bytes,
            .output_restart_overlap_bytes =
                m_config.parallel.process_output_restart_overlap_bytes,
            .mpi_external_reserve_bytes =
                m_config.parallel.process_mpi_unknown_reserve_bytes,
            .fftw_external_reserve_bytes =
                m_config.parallel.process_fftw_unknown_reserve_bytes,
            .hdf5_external_reserve_bytes =
                m_config.parallel.process_hdf5_unknown_reserve_bytes,
            .allocator_external_reserve_bytes =
                m_config.parallel.process_allocator_unknown_reserve_bytes,
            .safety_margin_fraction =
                m_config.parallel.process_memory_safety_margin_fraction,
        });

    // Replace the estimate helper's equal-rank aggregate with the exact
    // configured-topology reduction.  Particle decomposition can be uneven,
    // so production diagnostics must report the sum/max of the rank-local
    // preflights rather than local_required * world_size.  All ranks enter
    // these reductions before any rank is allowed to reject its process
    // budget, preventing a low-memory rank from throwing while peers continue
    // into later collectives.
    // Two collectives (one SUM vector, one MIN vector carrying ~value for the
    // maxima) replace the former five scalar reductions; the budget failure
    // vote is folded into the SUM vector because it depends only on the local
    // estimate.
    const std::uint64_t process_budget_bytes =
        m_config.parallel.process_memory_budget_bytes;
    const bool local_process_budget_failure =
        process_budget_bytes != 0U &&
        m_last_process_memory_estimate.budget_required_bytes >
            process_budget_bytes;
    std::array<std::uint64_t, 3U> process_memory_sums{
        m_last_process_memory_estimate.modeled_subtotal_bytes,
        m_last_process_memory_estimate.budget_required_bytes,
        local_process_budget_failure ? 1ULL : 0ULL};
    std::array<std::uint64_t, 2U> process_memory_maxima{
        ~m_last_process_memory_estimate.modeled_subtotal_bytes,
        ~m_last_process_memory_estimate.budget_required_bytes};
    mpi_context.allreduceSumUint64sInPlace(process_memory_sums);
    mpi_context.allreduceMinUint64sInPlace(process_memory_maxima);
    const std::uint64_t global_modeled_subtotal_bytes = process_memory_sums[0];
    const std::uint64_t rank_max_modeled_subtotal_bytes = ~process_memory_maxima[0];
    const std::uint64_t global_required_bytes = process_memory_sums[1];
    const std::uint64_t rank_max_required_bytes = ~process_memory_maxima[1];
    const std::uint64_t failed_process_budget_rank_count = process_memory_sums[2];

    m_last_process_memory_estimate.aggregate_required_bytes =
        global_required_bytes;
    m_last_process_memory_estimate.report.notes.push_back(
        "dmo_process_runtime_exact_aggregate_required_bytes=" +
        std::to_string(global_required_bytes));
    m_last_process_memory_estimate.report.notes.push_back(
        "dmo_process_runtime_rank_max_required_bytes=" +
        std::to_string(rank_max_required_bytes));
    auto& distributed = m_last_process_memory_estimate.report.distributed;
    distributed.valid = true;
    distributed.rank_count = std::max(m_runtime_topology.world_size, 1);
    distributed.local_owned_bytes =
        m_last_process_memory_estimate.modeled_subtotal_bytes;
    distributed.global_sum_owned_bytes = global_modeled_subtotal_bytes;
    distributed.rank_max_owned_bytes = rank_max_modeled_subtotal_bytes;
    if (global_modeled_subtotal_bytes != 0U) {
      const long double mean_bytes =
          static_cast<long double>(global_modeled_subtotal_bytes) /
          static_cast<long double>(distributed.rank_count);
      distributed.max_to_mean_imbalance_ratio = static_cast<double>(
          static_cast<long double>(rank_max_modeled_subtotal_bytes) /
          mean_bytes);
    }
    m_has_process_memory_estimate = true;

    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "memory.dmo_process_preflight",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "core.memory",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = hierarchical ? context.timeline_step.time_end_code : context.integrator_state.current_time_code,
          .scale_factor = hierarchical ? context.timeline_step.scale_factor_end : context.integrator_state.current_scale_factor,
          .message = "authoritative DMO per-rank process memory preflight",
          .payload = {
              {"known_owned_peak_bytes",
               std::to_string(m_last_process_memory_estimate.known_owned_peak_bytes)},
              {"external_unknown_reserve_bytes",
               std::to_string(m_last_process_memory_estimate.external_unknown_reserve_bytes)},
              {"modeled_subtotal_bytes",
               std::to_string(m_last_process_memory_estimate.modeled_subtotal_bytes)},
              {"safety_margin_bytes",
               std::to_string(m_last_process_memory_estimate.safety_margin_bytes)},
              {"required_bytes",
               std::to_string(m_last_process_memory_estimate.budget_required_bytes)},
              {"configured_process_budget_bytes",
               std::to_string(m_config.parallel.process_memory_budget_bytes)},
              {"rank_max_required_bytes", std::to_string(rank_max_required_bytes)},
              {"aggregate_required_bytes", std::to_string(global_required_bytes)},
              {"rank_count", std::to_string(distributed.rank_count)},
          },
      });
    }

    if (failed_process_budget_rank_count != 0U) {
      if (local_process_budget_failure) {
        gravity::enforceDmoProcessMemoryBudget(
            m_last_process_memory_estimate, process_budget_bytes);
      }
      throw std::runtime_error(
          "DMO process memory preflight failed on " +
          std::to_string(failed_process_budget_rank_count) +
          " rank(s); local_required_bytes=" +
          std::to_string(m_last_process_memory_estimate.budget_required_bytes) +
          ", configured_process_budget_bytes=" +
          std::to_string(process_budget_bytes) +
          ", rank_max_required_bytes=" +
          std::to_string(rank_max_required_bytes) +
          ", aggregate_required_bytes=" +
          std::to_string(global_required_bytes));
    }

    // After every rank has completed the authoritative process-memory
    // preflight, admit the one physical all-particle uint32 index lane needed
    // by the subsequent direct drift view. The pre-kick no longer constructs
    // an N-sized compatibility vector before this governed owner exists.
    admitAllParticleIndexScratch(force_targets_use_identity_scratch);

    if (force_targets_use_identity_scratch) {
      if (context.workspace == nullptr) {
        throw std::logic_error(
            "borrowed homogeneous DMO force targets require the governed workspace index lane");
      }
      force_target_particles = context.workspace->gravity_particle_index_scratch;
    }

    rebuildOwnedParticleCompactView(
        context,
        force_target_particles,
        source_prediction_epoch,
        context.active_set.particle_indices);
    // Refresh derived top-domain leaf bounds against the source snapshot this
    // solve will evaluate (post-prediction compact sources). Cached-kick
    // paths returned earlier, so this O(N) scan never runs on a pure cache hit.
    refreshAuthoritativeTopDomainGeometryForSolve(context, mpi_context);
    m_local_kick_particle_count = 0U;
    if (m_source_representation ==
        gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      m_local_kick_particle_count = context.active_set.particle_indices.size();
    } else {
      for (const std::uint32_t particle_index :
           context.active_set.particle_indices) {
        if (particle_index < m_owned_local_index_by_global.size() &&
            m_owned_local_index_by_global[particle_index] >= 0) {
          ++m_local_kick_particle_count;
        }
      }
    }
    const auto active_source_indices = activeSourceIndices();
    m_active_accel_x.assign(active_source_indices.size(), 0.0);
    m_active_accel_y.assign(active_source_indices.size(), 0.0);
    m_active_accel_z.assign(active_source_indices.size(), 0.0);
    if (m_source_representation ==
        gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      releaseVectorCapacity(m_active_slot_by_particle);
      releaseVectorCapacity(m_active_slot_by_cell);
    } else {
      m_active_slot_by_particle.assign(particle_count, -1);
      m_active_slot_by_cell.assign(cell_count, -1);
    }
    constexpr std::uint32_t no_target_row =
        std::numeric_limits<std::uint32_t>::max();
    for (std::size_t target_slot = 0;
         m_source_representation == gravity::GravitySourceRepresentation::kMaterializedGeneric &&
         target_slot < active_source_indices.size();
         ++target_slot) {
      if (targetParticleRow(target_slot) != no_target_row) {
        m_active_slot_by_particle[targetParticleRow(target_slot)] =
            static_cast<int>(target_slot);
      }
      if (targetCellRow(target_slot) != no_target_row) {
        m_active_slot_by_cell[targetCellRow(target_slot)] =
            static_cast<int>(target_slot);
      }
    }
    if (m_tree_pm_options.enable_zoom_long_range_correction) {
      ensureZoomMembershipLoaded(context.state);
      const double box_size_x = m_config.cosmology.box_size_x_mpc_comoving;
      const double box_size_y = m_config.cosmology.box_size_y_mpc_comoving;
      const double box_size_z = m_config.cosmology.box_size_z_mpc_comoving;
      m_source_is_high_res.assign(m_local_source_x.size(), 0U);
      m_active_is_high_res.assign(active_source_indices.size(), 0U);
      for (std::size_t i = 0; i < m_local_source_x.size(); ++i) {
        if (!m_zoom_high_res_particle_ids.empty()) {
          const std::uint32_t particle_row = m_local_source_particle_row[i];
          if (particle_row != no_target_row) {
            const std::uint64_t particle_id =
                context.state.particle_sidecar.particle_id[particle_row];
            m_source_is_high_res[i] =
                m_zoom_high_res_particle_ids.contains(particle_id) ? 1U : 0U;
            continue;
          }
          const std::uint32_t cell_row = m_local_source_cell_row[i];
          if (cell_row != no_target_row) {
            const auto parent_id =
                context.state.parentParticleIdForGasCellRow(cell_row);
            if (parent_id.has_value()) {
              m_source_is_high_res[i] =
                  m_zoom_high_res_particle_ids.contains(*parent_id) ? 1U : 0U;
              continue;
            }
          }
          // A parentless refined cell has no particle lineage authority. Fall
          // back to the explicit spatial zoom region rather than inventing one.
        }
        const double dx = m_local_source_x[i] - m_tree_pm_options.zoom_region_center_x_comoving;
        const double dy = m_local_source_y[i] - m_tree_pm_options.zoom_region_center_y_comoving;
        const double dz = m_local_source_z[i] - m_tree_pm_options.zoom_region_center_z_comoving;
        const double wrapped_dx = dx - box_size_x * std::nearbyint(dx / box_size_x);
        const double wrapped_dy = dy - box_size_y * std::nearbyint(dy / box_size_y);
        const double wrapped_dz = dz - box_size_z * std::nearbyint(dz / box_size_z);
        const double r = std::sqrt(
            wrapped_dx * wrapped_dx + wrapped_dy * wrapped_dy + wrapped_dz * wrapped_dz);
        m_source_is_high_res[i] =
            (r <= m_tree_pm_options.zoom_region_radius_comoving) ? 1U : 0U;
      }
      for (std::size_t i = 0; i < active_source_indices.size(); ++i) {
        m_active_is_high_res[i] = m_source_is_high_res[active_source_indices[i]];
      }
      m_tree_pm_options.source_is_high_res = m_source_is_high_res;
      m_tree_pm_options.active_is_high_res = m_active_is_high_res;
    } else {
      m_source_is_high_res.clear();
      m_active_is_high_res.clear();
      m_tree_pm_options.source_is_high_res = {};
      m_tree_pm_options.active_is_high_res = {};
    }

    if (hierarchical) {
      std::exception_ptr split_cache_failure;
      try {
        if (m_hierarchical_force_generation.size() != particle_count ||
            m_hierarchical_nearest_r2.size() != particle_count) {
          m_hierarchical_force_generation.assign(particle_count, m_force_cache_valid &&
              m_force_cache_source_generation == context.state.gravitySourceGeneration()
              ? context.state.gravitySourceGeneration() : 0U);
          m_hierarchical_pm_x.assign(particle_count, 0.0);
          m_hierarchical_pm_y.assign(particle_count, 0.0);
          m_hierarchical_pm_z.assign(particle_count, 0.0);
          m_hierarchical_tree_x.assign(particle_count, 0.0);
          m_hierarchical_tree_y.assign(particle_count, 0.0);
          m_hierarchical_tree_z.assign(particle_count, 0.0);
          m_hierarchical_nearest_r2.assign(particle_count, std::numeric_limits<double>::infinity());
        }
      } catch (...) { split_cache_failure = std::current_exception(); }
      FailureCoordinator(m_services).rethrowCollectiveFailure(split_cache_failure,
          "hierarchical split force cache materialization");
    }
    m_active_previous_acceleration_magnitude.clear();
    const bool relative_mac_cache_compatible =
        (m_tree_pm_options.tree_options.opening_criterion == gravity::TreeOpeningCriterion::kRelativeForceError ||
         m_tree_pm_options.acceptance_policy == gravity::TreePmAcceptancePolicy::kAdaptiveRelative) &&
        (m_force_cache_valid || adaptive_drift_history_compatible || hierarchical) &&
        m_force_cache_particle_index_generation == context.state.particleIndexGeneration() &&
        m_force_cache_cell_index_generation == context.state.cellIndexGeneration() &&
        m_force_cache_gas_identity_generation == context.state.gasCellIdentityGeneration() &&
        m_particle_accel_x.size() == particle_count &&
        m_particle_accel_y.size() == particle_count &&
        m_particle_accel_z.size() == particle_count &&
        m_cell_accel_x.size() == cell_count &&
        m_cell_accel_y.size() == cell_count &&
        m_cell_accel_z.size() == cell_count;
    if (relative_mac_cache_compatible) {
      m_active_previous_acceleration_magnitude.resize(
          active_source_indices.size(),
          std::numeric_limits<double>::quiet_NaN());
      for (std::size_t target_slot = 0; target_slot < active_source_indices.size(); ++target_slot) {
        const std::uint32_t particle_row = targetParticleRow(target_slot);
        if (particle_row != no_target_row &&
            particle_row < m_particle_force_cache_valid.size() &&
            (m_particle_force_cache_valid[particle_row] != 0U || adaptive_drift_history_compatible) &&
            (!hierarchical || (m_hierarchical_force_generation.size() == particle_count &&
             (m_hierarchical_force_generation[particle_row] >= context.hierarchical_kdk.coarse_source_generation ||
              context.hierarchical_kdk.force_only_synchronization)))) {
          m_active_previous_acceleration_magnitude[target_slot] = std::sqrt(
              m_particle_accel_x[particle_row] * m_particle_accel_x[particle_row] +
              m_particle_accel_y[particle_row] * m_particle_accel_y[particle_row] +
              m_particle_accel_z[particle_row] * m_particle_accel_z[particle_row]);
          continue;
        }
        const std::uint32_t cell_row = targetCellRow(target_slot);
        if (cell_row != no_target_row &&
            cell_row < m_cell_force_cache_valid.size() &&
            m_cell_force_cache_valid[cell_row] != 0U) {
          m_active_previous_acceleration_magnitude[target_slot] = std::sqrt(
              m_cell_accel_x[cell_row] * m_cell_accel_x[cell_row] +
              m_cell_accel_y[cell_row] * m_cell_accel_y[cell_row] +
              m_cell_accel_z[cell_row] * m_cell_accel_z[cell_row]);
        }
      }
    }

    const bool capture_pm = hierarchical && context.hierarchical_kdk.include_long_range_force;
    if (capture_pm && !isIdentityParticleSpan(active_source_indices, particle_count)) {
      throw std::logic_error("hierarchical PM endpoint must capture all owned target rows");
    }
    gravity::TreePmForceAccumulatorView accumulator{
        .active_particle_index = active_source_indices,
        .accel_x_comoving = m_active_accel_x,
        .accel_y_comoving = m_active_accel_y,
        .accel_z_comoving = m_active_accel_z,
        .previous_acceleration_magnitude_code = m_active_previous_acceleration_magnitude,
        .long_range_accel_x_comoving = capture_pm ? std::span<double>(m_hierarchical_pm_x) : std::span<double>{},
        .long_range_accel_y_comoving = capture_pm ? std::span<double>(m_hierarchical_pm_y) : std::span<double>{},
        .long_range_accel_z_comoving = capture_pm ? std::span<double>(m_hierarchical_pm_z) : std::span<double>{},
        .short_range_accel_x_comoving = hierarchical ? std::span<double>(m_hierarchical_tree_x) : std::span<double>{},
        .short_range_accel_y_comoving = hierarchical ? std::span<double>(m_hierarchical_tree_y) : std::span<double>{},
        .short_range_accel_z_comoving = hierarchical ? std::span<double>(m_hierarchical_tree_z) : std::span<double>{},
        .short_range_nearest_pair_r2_comoving =
            hierarchical ? std::span<double>(m_hierarchical_nearest_r2) : std::span<double>{},
    };

    // The integrator owns PM cadence state. At every rank-coordinated force
    // evaluation it issues a concrete sync event. The active particle set may
    // be local, but all ranks participate and inactive sources are predicted,
    // so production cadence one advances PM truth at that same force epoch.
    PmLongRangeKickDecision decision{
        .sync_surface = toPmSyncSurface(context.stage),
        .gravity_kick_opportunity = context.integrator_state.pm_sync_state.gravityKickOpportunity(),
        .refresh_long_range_field = false,
        .field_version = context.integrator_state.pm_sync_state.fieldVersion(),
        .last_refresh_opportunity = context.integrator_state.pm_sync_state.lastRefreshOpportunity(),
        .field_built_step_index = context.integrator_state.pm_sync_state.lastRefreshStepIndex(),
        .field_built_scale_factor = context.integrator_state.pm_sync_state.lastRefreshScaleFactor(),
    };

    double force_evaluation_scale_factor = 0.0;
    try {
      if (context.pm_refresh_directive.has_sync_event) {
        if (!context.pm_refresh_directive.cadence_opportunity_allowed) {
          throw std::runtime_error(
              "TreePM received a PM sync event outside an integrator-authorized refresh opportunity");
        }
        const core::PmSyncEvent sync_event{
            .gravity_kick_opportunity = context.pm_refresh_directive.gravity_kick_opportunity,
            .refresh_long_range_field = context.pm_refresh_directive.refresh_long_range_field,
            .field_version = context.pm_refresh_directive.field_version,
            .last_refresh_opportunity = context.pm_refresh_directive.last_refresh_opportunity,
            .field_built_step_index = context.pm_refresh_directive.field_built_step_index,
            .field_built_scale_factor = context.pm_refresh_directive.field_built_scale_factor,
        };
        pm_decision_consensus.requireEqual(
            "pm_sync_stage", static_cast<std::uint64_t>(context.pm_refresh_directive.sync_stage));
        pm_decision_consensus.requireEqual(
            "pm_refresh_reason", static_cast<std::uint64_t>(context.pm_refresh_directive.reason));
        pm_decision_consensus.requireEqual("gravity_kick_opportunity", sync_event.gravity_kick_opportunity);
        // Refresh vote equality across ranks replaces the former sum vote: the
        // decision is either unanimous refresh or unanimous reuse.
        pm_decision_consensus.requireEqual(
            "pm_refresh_vote", sync_event.refresh_long_range_field ? 1ULL : 0ULL);

        decision.gravity_kick_opportunity = sync_event.gravity_kick_opportunity;
        decision.refresh_long_range_field = sync_event.refresh_long_range_field;
        decision.field_version = sync_event.field_version;
        decision.last_refresh_opportunity = sync_event.last_refresh_opportunity;
        decision.field_built_step_index = sync_event.field_built_step_index;
        decision.field_built_scale_factor = sync_event.field_built_scale_factor;
      } else {
        if (!context.integrator_state.pm_long_range_field_valid || decision.field_version == 0U) {
          throw std::runtime_error(
              "TreePM local force refresh attempted before a valid long-range PM field existed");
        }
        if (context.pm_refresh_directive.cadence_opportunity_allowed) {
          throw std::runtime_error(
              "TreePM cadence opportunity was declared without an accompanying integrator-owned PM sync event");
        }
      }

      pm_decision_consensus.requireEqual("long_range_field_version", decision.field_version);
      pm_decision_consensus.requireEqual(
          "last_long_range_refresh_opportunity", decision.last_refresh_opportunity);
      pm_decision_consensus.requireEqual(
          "long_range_field_built_step_index", decision.field_built_step_index);
      pm_decision_consensus.requireEqual(
          "long_range_field_built_scale_factor_bits",
          std::bit_cast<std::uint64_t>(decision.field_built_scale_factor));
      force_evaluation_scale_factor = hierarchical
          ? context.pm_refresh_directive.force_evaluation_scale_factor : decision.field_built_scale_factor;
      if (decision.refresh_long_range_field) {
        force_evaluation_scale_factor =
            context.pm_refresh_directive.force_evaluation_scale_factor;
        if (!std::isfinite(force_evaluation_scale_factor) ||
            force_evaluation_scale_factor <= 0.0 ||
            force_evaluation_scale_factor != decision.field_built_scale_factor) {
          throw std::runtime_error(
              "TreePM refresh force-evaluation scale factor disagrees with integrator field-build metadata");
        }
      }
      if (!std::isfinite(force_evaluation_scale_factor) ||
          force_evaluation_scale_factor <= 0.0) {
        throw std::runtime_error(
            "TreePM force evaluation requires a finite positive scale factor");
      }
      pm_decision_consensus.requireEqual(
          "force_evaluation_scale_factor_bits",
          std::bit_cast<std::uint64_t>(force_evaluation_scale_factor));
    } catch (...) {
      // Rank-local failure: reported on every rank by verify() below, so no
      // rank throws while its peers wait in the consensus collective.
      pm_decision_consensus.recordLocalFailure(std::current_exception());
    }
    // One collective for every PM decision lane (formerly up to eighteen
    // min/max allreduces plus a refresh-vote sum per force solve).
    pm_decision_consensus.verify(mpi_context, world_size);
    // At a post-drift global refresh this is timeline_step.scale_factor_end,
    // even though IntegratorState.current_scale_factor remains at the step
    // beginning until commitStep(). A reuse retains the integrator-owned field
    // build scale as validity metadata; the returned Newtonian kernel itself is
    // scale-free and the KDK integrator owns the cosmological kick factor.
    m_tree_pm_options.pm_options.scale_factor =
        m_mode_policy.cosmological_comoving_frame
        ? force_evaluation_scale_factor
        : 1.0;
    // Dense-row generations are rank-local mirrors and can legitimately differ
    // after uneven compaction. The workflow advances this shared epoch exactly
    // once, on every rank, after a committed particle-ownership transition.
    m_tree_pm_options.decomposition_epoch = m_decomposition_epoch;
    // Gravity source identity is owned by authoritative SimulationState and is
    // intentionally independent from PM cadence sequencing. PM field version
    // answers "which mesh refresh?"; source generation answers "which source
    // state?"; force epoch answers "which evaluation opportunity?".
    m_tree_pm_options.source_generation = gravity::GravitySourceGeneration{
        context.state.gravitySourceGeneration()};
    m_tree_pm_options.source_layout_generation = m_source_representation ==
            gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
        ? context.state.particleIndexGeneration() : 0U;
    // Homogeneous DMO has no mass-changing source physics: row masses change
    // only by reorder/migration, which advance the dense particle index
    // generation. Offset by one so generation 0 still certifies. Materialized
    // generic sources (gas, star formation, feedback) stay uncertified (0).
    m_tree_pm_options.source_mass_generation = m_source_representation ==
            gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
        ? context.state.particleIndexGeneration() + 1U : 0U;
    // Geometry freshness metadata was published by the pre-solve refit; keep
    // the stamp aligned with the authoritative source generation this solve.
    publishAuthoritativeTopDomainLeavesSpan();
    m_tree_pm_options.pm_field_version = gravity::PmFieldVersion{decision.field_version};
    m_tree_pm_options.short_range_only = hierarchical && !capture_pm;
    // Fine-tick solves skip telemetry-only collectives; coarse PM endpoints and
    // global KDK keep them (option must be rank-uniform: it derives from the
    // collectively agreed hierarchical directive).
    m_tree_pm_options.reduce_global_diagnostics = !hierarchical || capture_pm;
    m_tree_pm_options.force_epoch = gravity::ForceEvaluationEpoch{
        .sequence = hierarchical ? context.state.gravitySourceGeneration() : decision.gravity_kick_opportunity,
        .scale_factor = force_evaluation_scale_factor,
    };
    if (!decision.refresh_long_range_field && !hierarchical &&
        context.integrator_state.pm_source_generation !=
            context.state.gravitySourceGeneration()) {
      throw std::runtime_error(
          "TreePM PM-field reuse rejected because authoritative gravity sources changed since the committed mesh refresh");
    }

    if (m_source_representation ==
        gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      if (m_borrowed_particle_index_generation !=
              context.state.particleIndexGeneration() ||
          m_borrowed_gravity_source_generation !=
              context.state.gravitySourceGeneration()) {
        throw std::runtime_error(
            "borrowed homogeneous DMO gravity view became stale before TreePM solve");
      }
      m_tree_pm_options.tree_options.softening.epsilon_comoving =
          m_tree_pm_species_softening.enabled
          ? m_tree_pm_species_softening.epsilon_comoving_by_species[
                static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)]
          : m_config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
    } else {
      m_tree_pm_options.tree_options.softening.epsilon_comoving =
          m_config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
    }

    const gravity::GravitySourceSnapshot source_snapshot{
        .pos_x_comoving = sourcePosX(context.state),
        .pos_y_comoving = sourcePosY(context.state),
        .pos_z_comoving = sourcePosZ(context.state),
        .mass_code = sourceMass(context.state),
        .species_tag = sourceSpecies(),
        .generation = m_tree_pm_options.source_generation,
        .evaluation_epoch = m_tree_pm_options.force_epoch,
        .particle_index_generation = context.state.particleIndexGeneration(),
        .representation = m_source_representation,
        .contains_predicted_coordinates = m_source_predicted_inactive_count != 0U,
    };
    source_snapshot.validate();

    if (m_particle_accel_x.size() != particle_count) {
      m_particle_accel_x.assign(particle_count, 0.0);
      m_particle_accel_y.assign(particle_count, 0.0);
      m_particle_accel_z.assign(particle_count, 0.0);
      m_particle_force_cache_valid.assign(particle_count, 0U);
    }
    if (m_cell_accel_x.size() != context.state.cells.size()) {
      m_cell_accel_x.assign(context.state.cells.size(), 0.0);
      m_cell_accel_y.assign(context.state.cells.size(), 0.0);
      m_cell_accel_z.assign(context.state.cells.size(), 0.0);
      m_cell_force_cache_valid.assign(context.state.cells.size(), 0U);
    } else if (m_cell_force_cache_valid.size() != context.state.cells.size()) {
      m_cell_force_cache_valid.assign(context.state.cells.size(), 0U);
    }
    const gravity::TreeSofteningView softening_view{
        .source_species_tag = sourceSpecies(),
        .source_particle_epsilon_comoving = m_source_representation ==
                gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
            ? std::span<const double>{}
            : std::span<const double>(m_local_source_softening_comoving.data(),
                                      m_local_source_softening_comoving.size()),
        .source_particle_epsilon_override_mask = m_source_representation ==
                gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
            ? std::span<const std::uint8_t>{}
            : std::span<const std::uint8_t>(m_local_source_softening_override_mask.data(),
                                            m_local_source_softening_override_mask.size()),
         .target_species_tag = {},
         .target_particle_epsilon_comoving = {},
         .target_particle_epsilon_override_mask = {},

        .species_policy = m_tree_pm_species_softening,
    };
    gravity::TreePmProfileEvent tree_pm_profile;
    m_tree_pm_coordinator.solveActiveSetWithPmCadence(
        source_snapshot.pos_x_comoving,
        source_snapshot.pos_y_comoving,
        source_snapshot.pos_z_comoving,
        source_snapshot.mass_code,
        accumulator,
        m_tree_pm_options,
        decision.refresh_long_range_field,
        &tree_pm_profile,
        &m_last_tree_pm_diagnostics,
        softening_view);
    const std::uint64_t expected_pm_solve_count =
        decision.refresh_long_range_field ? 1U : 0U;
    if (m_last_tree_pm_diagnostics.pm_solve_count != expected_pm_solve_count ||
        m_last_tree_pm_diagnostics.pm_reuse_count !=
            (m_tree_pm_options.short_range_only ? 0U : 1U - expected_pm_solve_count)) {
      throw std::runtime_error(
          "TreePM solver PM solve/reuse outcome disagrees with the integrator-owned cadence directive");
    }

    m_last_decomposition_measurements = parallel::DecompositionRuntimeMeasurements{
        .tree_pair_evaluations_recent = m_last_tree_pm_diagnostics.residual_pair_evaluations,
        .tree_remote_request_bytes_recent = m_last_tree_pm_diagnostics.residual_remote_request_bytes +
            m_last_tree_pm_diagnostics.residual_remote_response_bytes,
        .pm_mesh_cells_touched_recent = m_tree_pm_options.short_range_only ? 0U :
            static_cast<std::uint64_t>(m_tree_pm_coordinator.slabLayout().localCellCount()),
        .pm_fft_transpose_bytes_recent = tree_pm_profile.pm_profile.fft_transpose_bytes,
        .amr_patch_cells_updated_recent = static_cast<std::uint64_t>(context.state.cells.size()),
        .hydro_face_fluxes_recent = 0,
        .ghost_exchange_bytes_recent = gravity_ghost_refresh.sent_bytes + gravity_ghost_refresh.received_bytes,
        // tree_short_range_ms is the inclusive residual-phase wall time
        // (source preprocess + build + multipoles + traversal + LET). The
        // build/multipole/traversal subphase fields are its children and must
        // not be added again.
        .tree_wall_ms_recent = tree_pm_profile.tree_short_range_ms,
        // total_ms is the single inclusive PM-phase measurement; assign/fft/
        // poisson/gradient/interpolate subphases overlap and must not be summed.
        .pm_wall_ms_recent = tree_pm_profile.pm_profile.total_ms,
        .gpu_kernel_ms_recent = tree_pm_profile.pm_profile.device_kernel_ms,
        .accelerator_occupancy_fraction_recent = (tree_pm_profile.pm_profile.device_kernel_ms > 0.0) ? 1.0 : 0.0,
        .has_measurements = true,
        .has_spatial_tree_work = m_config.parallel.decomposition_spatial_work_enabled &&
            m_last_tree_pm_diagnostics.spatial_work_history_solves != 0U,
        .spatial_tree_work_per_target = m_last_tree_pm_diagnostics.spatial_work_per_target,
        .incoming_tree_pair_evaluations_recent = m_last_tree_pm_diagnostics.incoming_remote_pair_evaluations,
    };

    const bool allow_heavy_reference_checks =
        m_config.analysis.diagnostics_execution_policy ==
        core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
    GravityHealthSummary gravity_health = validateGravityHealth(
        context,
        decision,
        allow_heavy_reference_checks);

    if (decision.refresh_long_range_field) {
      if (hierarchical) {
        m_hierarchical_pm_source_generation = context.state.gravitySourceGeneration();
        m_hierarchical_pm_field_version = decision.field_version;
      }
      m_has_long_range_field = true;
      ++m_long_range_refresh_count;
      if (m_cycle_transaction_open) ++m_pending_long_range_refreshes;
    } else if (!m_tree_pm_options.short_range_only) {
      ++m_long_range_reuse_count;
      if (m_cycle_transaction_open) ++m_pending_long_range_reuses;
    }
    context.pm_refresh_directive.solver_executed = true;
    m_last_committed_field_version = decision.field_version;

    const std::uint64_t inactive_particles_skipped = static_cast<std::uint64_t>(
        context.state.particles.size() - m_local_kick_particle_count);
    // Recent detail is bounded; refresh/reuse counters above remain lifetime
    // aggregates for the final report. No run-length-sized terminal copy.
    // Inside a multirate PM cycle transaction the record is journaled and
    // published only when the cycle is accepted.
    appendCadenceRecord(ReferenceWorkflowReport::TreePmCadenceRecord{
        .step_index = context.integrator_state.step_index,
        .stage_name = std::string(core::integrationStageName(context.stage)),
        .pm_sync_surface = std::string(pmSyncSurfaceName(decision.sync_surface)),
        .pm_refresh_reason = std::string(pmRefreshReasonName(context.pm_refresh_directive.reason)),
        .gravity_kick_opportunity = decision.gravity_kick_opportunity,
        .field_version = decision.field_version,
        .last_refresh_opportunity = decision.last_refresh_opportunity,
        .field_built_step_index = decision.field_built_step_index,
        .field_built_scale_factor = decision.field_built_scale_factor,
        .field_age_in_kick_opportunities =
            decision.gravity_kick_opportunity - decision.last_refresh_opportunity,
        .active_particles_kicked = hierarchical ? 0U : static_cast<std::uint64_t>(m_local_kick_particle_count),
        .inactive_particles_skipped = inactive_particles_skipped,
        .refreshed_long_range_field = decision.refresh_long_range_field,
        .short_range_only = m_tree_pm_options.short_range_only,
        .active_force_targets = static_cast<std::uint64_t>(m_local_kick_particle_count),
    });
    if (context.profiler_session != nullptr) {
      context.profiler_session->counters().addCount("treepm.local.targets", m_last_tree_pm_diagnostics.local_traversal.targets);
      context.profiler_session->counters().addCount("treepm.local.visited_nodes", m_last_tree_pm_diagnostics.local_traversal.visited_nodes);
      context.profiler_session->counters().addCount("treepm.local.opened_nodes", m_last_tree_pm_diagnostics.local_traversal.opened_nodes);
      context.profiler_session->counters().addCount("treepm.local.accepted_internal_multipoles", m_last_tree_pm_diagnostics.local_traversal.accepted_internal_multipoles);
      context.profiler_session->counters().addCount("treepm.local.accepted_leaves", m_last_tree_pm_diagnostics.local_traversal.accepted_leaves);
      context.profiler_session->counters().addCount("treepm.local.direct_pair_evaluations", m_last_tree_pm_diagnostics.local_traversal.direct_pair_evaluations);
      context.profiler_session->counters().addCount("treepm.local.cutoff_pruned_nodes", m_last_tree_pm_diagnostics.local_traversal.cutoff_pruned_nodes);
      context.profiler_session->counters().addCount("treepm.local.cutoff_skipped_pairs", m_last_tree_pm_diagnostics.local_traversal.cutoff_skipped_pairs);
      context.profiler_session->counters().addCount("treepm.local.selected_mac_rejections", m_last_tree_pm_diagnostics.local_traversal.selected_mac_rejections);
      context.profiler_session->counters().addCount("treepm.local.relative_mac_rejections", m_last_tree_pm_diagnostics.local_traversal.relative_mac_rejections);
      context.profiler_session->counters().addCount("treepm.local.skipped_mac_evaluations", m_last_tree_pm_diagnostics.local_traversal.skipped_mac_evaluations);
      context.profiler_session->counters().addCount("treepm.local.maximum_angle_rejections", m_last_tree_pm_diagnostics.local_traversal.maximum_angle_rejections);
      context.profiler_session->counters().addCount("treepm.local.strict_envelope_rejections", m_last_tree_pm_diagnostics.local_traversal.strict_envelope_rejections);
      context.profiler_session->counters().addCount("treepm.local.softening_rejections", m_last_tree_pm_diagnostics.local_traversal.softening_rejections);
      context.profiler_session->counters().addCount("treepm.local.near_node_rejections", m_last_tree_pm_diagnostics.local_traversal.near_node_rejections);
      context.profiler_session->counters().addCount("treepm.local.cutoff_containment_rejections", m_last_tree_pm_diagnostics.local_traversal.cutoff_containment_rejections);
      context.profiler_session->counters().addCount("treepm.local.geometric_history_fallbacks", m_last_tree_pm_diagnostics.local_traversal.geometric_history_fallbacks);
      context.profiler_session->counters().addCount("treepm.incoming.targets", m_last_tree_pm_diagnostics.incoming_traversal.targets);
      context.profiler_session->counters().addCount("treepm.incoming.visited_nodes", m_last_tree_pm_diagnostics.incoming_traversal.visited_nodes);
      context.profiler_session->counters().addCount("treepm.incoming.opened_nodes", m_last_tree_pm_diagnostics.incoming_traversal.opened_nodes);
      context.profiler_session->counters().addCount("treepm.incoming.accepted_internal_multipoles", m_last_tree_pm_diagnostics.incoming_traversal.accepted_internal_multipoles);
      context.profiler_session->counters().addCount("treepm.incoming.accepted_leaves", m_last_tree_pm_diagnostics.incoming_traversal.accepted_leaves);
      context.profiler_session->counters().addCount("treepm.incoming.direct_pair_evaluations", m_last_tree_pm_diagnostics.incoming_traversal.direct_pair_evaluations);
      context.profiler_session->counters().addCount("treepm.incoming.cutoff_pruned_nodes", m_last_tree_pm_diagnostics.incoming_traversal.cutoff_pruned_nodes);
      context.profiler_session->counters().addCount("treepm.incoming.cutoff_skipped_pairs", m_last_tree_pm_diagnostics.incoming_traversal.cutoff_skipped_pairs);
      context.profiler_session->counters().addCount("treepm.incoming.selected_mac_rejections", m_last_tree_pm_diagnostics.incoming_traversal.selected_mac_rejections);
      context.profiler_session->counters().addCount("treepm.incoming.relative_mac_rejections", m_last_tree_pm_diagnostics.incoming_traversal.relative_mac_rejections);
      context.profiler_session->counters().addCount("treepm.incoming.skipped_mac_evaluations", m_last_tree_pm_diagnostics.incoming_traversal.skipped_mac_evaluations);
      context.profiler_session->counters().addCount("treepm.incoming.maximum_angle_rejections", m_last_tree_pm_diagnostics.incoming_traversal.maximum_angle_rejections);
      context.profiler_session->counters().addCount("treepm.incoming.strict_envelope_rejections", m_last_tree_pm_diagnostics.incoming_traversal.strict_envelope_rejections);
      context.profiler_session->counters().addCount("treepm.incoming.softening_rejections", m_last_tree_pm_diagnostics.incoming_traversal.softening_rejections);
      context.profiler_session->counters().addCount("treepm.incoming.near_node_rejections", m_last_tree_pm_diagnostics.incoming_traversal.near_node_rejections);
      context.profiler_session->counters().addCount("treepm.incoming.cutoff_containment_rejections", m_last_tree_pm_diagnostics.incoming_traversal.cutoff_containment_rejections);
      context.profiler_session->counters().addCount("treepm.incoming.geometric_history_fallbacks", m_last_tree_pm_diagnostics.incoming_traversal.geometric_history_fallbacks);
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = m_tree_pm_options.short_range_only ? "gravity.short_range_force" : "gravity.pm_long_range_field",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = hierarchical ? context.timeline_step.time_end_code : context.integrator_state.current_time_code,
          .scale_factor = hierarchical ? context.timeline_step.scale_factor_end : context.integrator_state.current_scale_factor,
          .message = m_tree_pm_options.short_range_only
              ? "active short-range force evaluated; PM mesh solve and interpolation omitted by the coarse split operator"
              : (decision.refresh_long_range_field
                  ? "PM long-range field refreshed for gravity force evaluation"
                  : "PM long-range field reused for gravity force evaluation"),
          .payload = {{"short_range_only", m_tree_pm_options.short_range_only ? "true" : "false"},
                      {"hierarchical_kdk", hierarchical ? "true" : "false"},
                      {"active_force_targets", std::to_string(m_local_kick_particle_count)},
                      {"stage", std::string(core::integrationStageName(context.stage))},
                      {"pm_sync_surface", std::string(pmSyncSurfaceName(decision.sync_surface))},
                      {"gravity_kick_opportunity", std::to_string(decision.gravity_kick_opportunity)},
                      {"field_version", std::to_string(decision.field_version)},
                      {"field_built_step_index", std::to_string(decision.field_built_step_index)},
                      {"field_built_scale_factor", formatRuntimeDouble(decision.field_built_scale_factor)},
                      {"pm_update_cadence_steps", std::to_string(context.integrator_state.pm_sync_state.cadenceSteps())},
                      {"pm_grid", std::to_string(m_pm_grid_shape.nx) + "x" + std::to_string(m_pm_grid_shape.ny) +
                              "x" + std::to_string(m_pm_grid_shape.nz)},
                      {"pm_assignment_scheme", m_pm_assignment_scheme},
                      {"pm_window_deconvolution", m_config.numerics.treepm_enable_window_deconvolution ? "true" : "false"},
                      {"pm_force_gradient", m_config.numerics.treepm_pm_force_gradient ==
                              core::TreePmForceGradient::kFiniteDifference4 ? "finite_difference_4" : "spectral"},
                      {"asmth_cells", formatRuntimeDouble(m_config.numerics.treepm_asmth_cells)},
                      {"rcut_cells", formatRuntimeDouble(m_config.numerics.treepm_rcut_cells)},
                      {"tree_opening_criterion", treePmOpeningCriterionName(
                          m_config.numerics.treepm_tree_opening_criterion)},
                      {"tree_opening_theta", formatRuntimeDouble(
                          m_config.numerics.treepm_tree_opening_theta)},
                      {"tree_relative_force_tolerance", formatRuntimeDouble(
                          m_config.numerics.treepm_tree_relative_force_tolerance)},
                      {"tree_relative_force_acceleration_floor", formatRuntimeDouble(
                          m_config.numerics.treepm_tree_relative_force_acceleration_floor)},
                      {"relative_mac_previous_force_available",
                          relative_mac_cache_compatible ? "true" : "false"},
                      {"mesh_spacing_x_mpc_comoving", formatRuntimeDouble(m_mesh_spacing_x_mpc_comoving)},
                      {"mesh_spacing_y_mpc_comoving", formatRuntimeDouble(m_mesh_spacing_y_mpc_comoving)},
                      {"mesh_spacing_z_mpc_comoving", formatRuntimeDouble(m_mesh_spacing_z_mpc_comoving)},
                      {"split_scale_mpc_comoving", formatRuntimeDouble(m_tree_pm_options.split_policy.split_scale_comoving)},
                      {"cutoff_radius_mpc_comoving", formatRuntimeDouble(m_tree_pm_options.split_policy.cutoff_radius_comoving)},
                      {"softening_policy", describeSofteningPolicy(m_config, &context.state)},
                      {"softening_kernel", "plummer"},
                      {"softening_epsilon_kpc_comoving", formatRuntimeDouble(m_config.numerics.gravity_softening_kpc_comoving)},
                      {"pm_fft_backend", m_pm_backend},
                      {"pm_backend_capability", std::string(gravity::pmBackendCapabilityName(gravity::pmBackendCapability()))},
                      {"gravity_acceptance_profile_id", "unverified_current_source"},
                      {"active_particles_kicked", std::to_string(hierarchical ? 0U : m_local_kick_particle_count)},
                      {"gravity_source_representation",
                       m_source_representation == gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
                           ? "borrowed_homogeneous_dmo"
                           : "materialized_generic"},
                      {"borrowed_source_count",
                       std::to_string(m_source_representation == gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
                           ? source_snapshot.size()
                           : 0U)},
                      {"materialized_source_count",
                       std::to_string(m_source_representation == gravity::GravitySourceRepresentation::kMaterializedGeneric
                           ? source_snapshot.size()
                           : 0U)},
                      {"inactive_particles_skipped", std::to_string(inactive_particles_skipped)},
                      {"ghost_refresh_sent_bytes", std::to_string(gravity_ghost_refresh.sent_bytes)},
                      {"ghost_refresh_received_bytes", std::to_string(gravity_ghost_refresh.received_bytes)},
                      {"ghost_refresh_committed_slots", std::to_string(gravity_ghost_refresh.committed_slots)},
                      {"predicted_inactive_source_particles", std::to_string(m_source_predicted_inactive_count)},
                      {"predicted_inactive_source_count", std::to_string(m_source_predicted_inactive_count)},
                      {"source_prediction_drift_integrals", std::to_string(m_source_prediction_drift_integrals)},
                      {"predicted_inactive_sources_required",
                          context.pm_refresh_directive.requires_predicted_inactive_sources ? "true" : "false"},
                      {"pm_refresh_reason", std::string(pmRefreshReasonName(context.pm_refresh_directive.reason))},
                      {"pm_refresh_force_eval_scale_factor", formatRuntimeDouble(context.pm_refresh_directive.force_evaluation_scale_factor)},
                      {"pm_reuse_policy", m_config.numerics.treepm_pm_reuse_policy == core::TreePmPmReusePolicy::kMultirateOuterKick
                           ? "multirate_outer_kick"
                           : (m_tree_pm_options.allow_stale_long_range_field ? "stale_within_cadence" : "exact")},
                      {"pm_field_age_steps", std::to_string(
                          context.integrator_state.step_index >= decision.field_built_step_index
                              ? context.integrator_state.step_index - decision.field_built_step_index : 0U)},
                      {"refreshed_long_range_field", decision.refresh_long_range_field ? "true" : "false"}},
      });
      const auto decomposition_descriptor = gravity::describePmDecomposition(
          m_pm_grid_shape,
          m_tree_pm_coordinator.slabLayout(),
          m_tree_pm_options.pm_options.decomposition_mode);
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.treepm_let",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = hierarchical ? context.timeline_step.time_end_code : context.integrator_state.current_time_code,
          .scale_factor = hierarchical ? context.timeline_step.scale_factor_end : context.integrator_state.current_scale_factor,
          .message = "locality-driven TreePM short-range communication metrics",
          .payload = {
              {"tree_acceptance_policy", m_config.numerics.treepm_adaptive_acceptance_enabled ? "adaptive_relative" : "strict_reference"},
              {"tree_adaptive_error_estimator",
               m_config.numerics.treepm_adaptive_error_estimator == core::TreePmAdaptiveErrorEstimator::kRelativeMonopole
                   ? "relative_monopole" : "residual_derivative"},
              {"tree_cutoff_node_policy",
               m_config.numerics.treepm_cutoff_node_policy == core::TreePmCutoffNodePolicy::kComDistance
                   ? "com_distance" : "contain"},
              {"full_residual_lookup_requested", m_tree_pm_options.full_residual_lookup_enabled ? "true" : "false"},
              {"radial_coefficients_enabled", m_tree_pm_options.radial_coefficients_enabled ? "true" : "false"},
              {"leaf_simd_requested", m_tree_pm_options.leaf_simd_enabled ? "true" : "false"},
              {"residual_lookup_targets", std::to_string(m_last_tree_pm_diagnostics.local_traversal.residual_lookup_targets)},
              {"simd_leaf_targets", std::to_string(m_last_tree_pm_diagnostics.local_traversal.simd_leaf_targets)},
              {"group_walk_requested", m_last_tree_pm_diagnostics.group_walk_requested != 0U ? "true" : "false"},
              {"group_walk_status", gravity::treePmGroupWalkStatusName(
                  static_cast<gravity::TreePmGroupWalkStatus>(m_last_tree_pm_diagnostics.group_walk_status))},
              {"group_walk_target_size", std::to_string(m_tree_pm_options.group_walk_target_group_size)},
              {"group_walk_groups", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_groups)},
              {"group_walk_grouped_targets", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_grouped_targets)},
              {"group_walk_target_group_size_max", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_target_group_size_max)},
              {"group_walk_fallback_targets", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_fallback_targets)},
              {"group_walk_fallback_no_leaf", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_fallback_no_leaf)},
              {"group_walk_fallback_softening", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_fallback_softening)},
              {"group_walk_fallback_history", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_fallback_history)},
              {"group_walk_fallback_geometry", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_fallback_geometry)},
              {"group_walk_candidate_leaves", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_candidate_leaves)},
              {"group_walk_culled_leaves", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_culled_leaves)},
              {"group_walk_candidate_particles", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_candidate_particles)},
              {"group_walk_list_entries_max", std::to_string(m_last_tree_pm_diagnostics.local_traversal.group_walk_list_entries_max)},
              {"incoming_group_walk_groups", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.group_walk_groups)},
              {"incoming_group_walk_grouped_targets", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.group_walk_grouped_targets)},
              {"incoming_group_walk_fallback_targets", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.group_walk_fallback_targets)},
              {"incoming_group_walk_culled_leaves", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.group_walk_culled_leaves)},
              {"topology_cache_reuses", std::to_string(m_last_tree_pm_diagnostics.topology_cache_reuses)},
              {"topology_cache_rebuilds", std::to_string(m_last_tree_pm_diagnostics.topology_cache_rebuilds)},
              {"mass_uniformity_scans", std::to_string(m_last_tree_pm_diagnostics.mass_uniformity_scans)},
              {"mass_uniformity_reuses", std::to_string(m_last_tree_pm_diagnostics.mass_uniformity_reuses)},
              {"uniform_source_mass_lane_released", std::to_string(m_last_tree_pm_diagnostics.uniform_source_mass_lane_released)},
              {"skipped_geometry_evaluations", std::to_string(m_last_tree_pm_diagnostics.local_traversal.skipped_geometry_evaluations)},
              {"residual_block_size", std::to_string(m_tree_pm_options.residual_block_size)},
              {"tree_max_leaf_size", std::to_string(m_tree_pm_options.tree_options.max_leaf_size)},
              {"gaussian_pair_kernel", m_tree_pm_options.gaussian_pair_lookup_enabled ? "lookup_experimental" : "analytic"},
              {"mac_accounting", m_tree_pm_options.full_mac_diagnostics ? "full_overlapping" : "fast_guard_short_circuit"},
              {"adaptive_maximum_opening_angle", formatRuntimeDouble(m_tree_pm_options.adaptive_maximum_opening_angle)},
              {"tree_full_rebuild_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.full_rebuild_ms)},
              {"tree_build_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.build_ms)},
              {"tree_source_identity_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.source_identity_ms)},
              {"tree_source_softening_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.source_softening_ms)},
              {"tree_morton_ordering_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.morton_ordering_ms)},
              {"tree_topology_build_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.topology_build_ms)},
              {"tree_multipole_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.multipole_ms)},
              {"tree_source_preprocess_ms", formatRuntimeDouble(tree_pm_profile.source_preprocess_ms)},
              {"tree_actual_node_count", std::to_string(tree_pm_profile.tree_profile.actual_node_count)},
              {"tree_maximum_depth", std::to_string(tree_pm_profile.tree_profile.maximum_depth)},
              {"tree_source_identity_used_generation", std::to_string(tree_pm_profile.tree_profile.source_identity_used_generation)},
              {"tree_build_count", std::to_string(tree_pm_profile.tree_profile.build_count)},
              {"tree_reuse_attempted", std::to_string(tree_pm_profile.tree_profile.reuse_attempted)},
              {"tree_reuse_accepted", std::to_string(tree_pm_profile.tree_profile.reuse_accepted)},
              {"tree_rebuild_reason", std::to_string(tree_pm_profile.tree_profile.rebuild_reason)},
              {"tree_topology_validity_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.topology_validity_ms)},
              {"tree_refit_attempted", std::to_string(tree_pm_profile.tree_profile.refit_attempted)},
              {"tree_refit_accepted", std::to_string(tree_pm_profile.tree_profile.refit_accepted)},
              {"tree_refit_rejected", std::to_string(tree_pm_profile.tree_profile.refit_rejected)},
              {"tree_geometry_preparation_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.geometry_preparation_ms)},
              {"tree_refit_ms", formatRuntimeDouble(tree_pm_profile.tree_profile.refit_ms)},
              {"spatial_work_history_solves", std::to_string(m_last_tree_pm_diagnostics.spatial_work_history_solves)},
              {"local_targets", std::to_string(m_last_tree_pm_diagnostics.local_traversal.targets)},
              {"local_visited_nodes", std::to_string(m_last_tree_pm_diagnostics.local_traversal.visited_nodes)},
              {"local_opened_nodes", std::to_string(m_last_tree_pm_diagnostics.local_traversal.opened_nodes)},
              {"local_accepted_internal_multipoles", std::to_string(m_last_tree_pm_diagnostics.local_traversal.accepted_internal_multipoles)},
              {"local_accepted_leaves", std::to_string(m_last_tree_pm_diagnostics.local_traversal.accepted_leaves)},
              {"local_direct_pair_evaluations", std::to_string(m_last_tree_pm_diagnostics.local_traversal.direct_pair_evaluations)},
              {"local_cutoff_pruned_nodes", std::to_string(m_last_tree_pm_diagnostics.local_traversal.cutoff_pruned_nodes)},
              {"local_cutoff_skipped_pairs", std::to_string(m_last_tree_pm_diagnostics.local_traversal.cutoff_skipped_pairs)},
              {"local_selected_mac_rejections", std::to_string(m_last_tree_pm_diagnostics.local_traversal.selected_mac_rejections)},
              {"local_relative_mac_rejections", std::to_string(m_last_tree_pm_diagnostics.local_traversal.relative_mac_rejections)},
              {"local_skipped_mac_evaluations", std::to_string(m_last_tree_pm_diagnostics.local_traversal.skipped_mac_evaluations)},
              {"local_maximum_angle_rejections", std::to_string(m_last_tree_pm_diagnostics.local_traversal.maximum_angle_rejections)},
              {"local_strict_envelope_rejections", std::to_string(m_last_tree_pm_diagnostics.local_traversal.strict_envelope_rejections)},
              {"local_softening_rejections", std::to_string(m_last_tree_pm_diagnostics.local_traversal.softening_rejections)},
              {"local_near_node_rejections", std::to_string(m_last_tree_pm_diagnostics.local_traversal.near_node_rejections)},
              {"local_cutoff_containment_rejections", std::to_string(m_last_tree_pm_diagnostics.local_traversal.cutoff_containment_rejections)},
              {"local_geometric_history_fallbacks", std::to_string(m_last_tree_pm_diagnostics.local_traversal.geometric_history_fallbacks)},
              {"incoming_targets", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.targets)},
              {"incoming_visited_nodes", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.visited_nodes)},
              {"incoming_opened_nodes", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.opened_nodes)},
              {"incoming_accepted_internal_multipoles", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.accepted_internal_multipoles)},
              {"incoming_accepted_leaves", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.accepted_leaves)},
              {"incoming_direct_pair_evaluations", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.direct_pair_evaluations)},
              {"incoming_cutoff_pruned_nodes", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.cutoff_pruned_nodes)},
              {"incoming_cutoff_skipped_pairs", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.cutoff_skipped_pairs)},
              {"incoming_selected_mac_rejections", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.selected_mac_rejections)},
              {"incoming_relative_mac_rejections", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.relative_mac_rejections)},
              {"incoming_skipped_mac_evaluations", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.skipped_mac_evaluations)},
              {"incoming_maximum_angle_rejections", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.maximum_angle_rejections)},
              {"incoming_strict_envelope_rejections", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.strict_envelope_rejections)},
              {"incoming_softening_rejections", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.softening_rejections)},
              {"incoming_near_node_rejections", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.near_node_rejections)},
              {"incoming_cutoff_containment_rejections", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.cutoff_containment_rejections)},
              {"incoming_geometric_history_fallbacks", std::to_string(m_last_tree_pm_diagnostics.incoming_traversal.geometric_history_fallbacks)},
              {"worker_targets_min", formatRuntimeDouble(m_last_tree_pm_diagnostics.worker_targets_min)},
              {"worker_targets_max", formatRuntimeDouble(m_last_tree_pm_diagnostics.worker_targets_max)},
              {"worker_visits_max", formatRuntimeDouble(m_last_tree_pm_diagnostics.worker_visits_max)},
              {"worker_pairs_max", formatRuntimeDouble(m_last_tree_pm_diagnostics.worker_pairs_max)},
              {"worker_multipoles_max", formatRuntimeDouble(m_last_tree_pm_diagnostics.worker_multipoles_max)},
              {"worker_work_ms_min", formatRuntimeDouble(m_last_tree_pm_diagnostics.worker_work_ms_min)},
              {"worker_work_ms_max", formatRuntimeDouble(m_last_tree_pm_diagnostics.worker_work_ms_max)},
              {"worker_work_ms_sum", formatRuntimeDouble(m_last_tree_pm_diagnostics.worker_work_ms_sum)},
              {"worker_work_ms_mean", formatRuntimeDouble(m_last_tree_pm_diagnostics.worker_region_count == 0U ? 0.0 : m_last_tree_pm_diagnostics.worker_work_ms_sum / static_cast<double>(m_last_tree_pm_diagnostics.worker_region_count))},
              {"worker_region_count", std::to_string(m_last_tree_pm_diagnostics.worker_region_count)},
              {"pm_assign_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.assign_ms)},
              {"pm_density_routing_wait_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.density_routing_wait_ms)},
              {"pm_fft_forward_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.fft_forward_ms)},
              {"pm_poisson_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.poisson_ms)},
              {"pm_gradient_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.gradient_ms)},
              {"pm_fft_inverse_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.fft_inverse_ms)},
              {"pm_interpolate_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.interpolate_ms)},
              {"pm_halo_exchange_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.halo_exchange_ms)},
              {"pm_total_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.total_ms)},
              {"pm_gradient_factor_ms_x", formatRuntimeDouble(tree_pm_profile.pm_profile.gradient_factor_ms[0])},
              {"pm_inverse_axis_ms_x", formatRuntimeDouble(tree_pm_profile.pm_profile.inverse_axis_ms[0])},
              {"pm_normalization_axis_ms_x", formatRuntimeDouble(tree_pm_profile.pm_profile.normalization_axis_ms[0])},
              {"pm_gradient_factor_ms_y", formatRuntimeDouble(tree_pm_profile.pm_profile.gradient_factor_ms[1])},
              {"pm_inverse_axis_ms_y", formatRuntimeDouble(tree_pm_profile.pm_profile.inverse_axis_ms[1])},
              {"pm_normalization_axis_ms_y", formatRuntimeDouble(tree_pm_profile.pm_profile.normalization_axis_ms[1])},
              {"pm_gradient_factor_ms_z", formatRuntimeDouble(tree_pm_profile.pm_profile.gradient_factor_ms[2])},
              {"pm_inverse_axis_ms_z", formatRuntimeDouble(tree_pm_profile.pm_profile.inverse_axis_ms[2])},
              {"pm_normalization_axis_ms_z", formatRuntimeDouble(tree_pm_profile.pm_profile.normalization_axis_ms[2])},
              {"pm_decomposition_architecture", std::string(
                  gravity::pmDecompositionTopologyName(decomposition_descriptor.topology))},
              {"top_level_domain_leaf_count", std::to_string(m_last_tree_pm_diagnostics.top_level_domain_leaf_count)},
              {"authoritative_domain_leaf_count", std::to_string(m_last_tree_pm_diagnostics.authoritative_domain_leaf_count)},
              {"domain_geometry_source_generation", std::to_string(
                  m_last_tree_pm_diagnostics.domain_geometry_source_generation)},
              {"current_gravity_source_generation", std::to_string(
                  m_last_tree_pm_diagnostics.current_gravity_source_generation)},
              {"domain_geometry_fresh", m_last_tree_pm_diagnostics.domain_geometry_fresh != 0U ? "true" : "false"},
              {"domain_geometry_refreshed", m_domain_geometry_refreshed_this_solve ? "true" : "false"},
              {"domain_geometry_fallback_used",
               m_last_tree_pm_diagnostics.domain_geometry_fallback_used != 0U ? "true" : "false"},
              {"domain_geometry_fallback_reason", gravity::treePmDomainGeometryFallbackReasonName(
                  static_cast<gravity::TreePmDomainGeometryFallbackReason>(
                      m_last_tree_pm_diagnostics.domain_geometry_fallback_reason))},
              {"domain_geometry_refreshed_leaf_count",
               std::to_string(m_domain_geometry_refreshed_leaf_count)},
              {"domain_geometry_out_of_seed_range_source_count",
               std::to_string(m_domain_geometry_out_of_seed_range_source_count)},
              {"domain_geometry_uncovered_source_count", std::to_string(
                  m_last_tree_pm_diagnostics.domain_geometry_uncovered_source_count)},
              {"openmp_compiled", m_last_tree_pm_diagnostics.openmp_compiled != 0U ? "true" : "false"},
               {"openmp_configured_workers", std::to_string(m_last_tree_pm_diagnostics.openmp_configured_workers)},

              {"openmp_observed_workers", std::to_string(m_last_tree_pm_diagnostics.openmp_observed_workers)},
              {"residual_local_target_count", std::to_string(m_last_tree_pm_diagnostics.residual_local_target_count)},
              {"residual_incoming_target_count",
               std::to_string(m_last_tree_pm_diagnostics.residual_incoming_target_count)},
              {"residual_worker_scratch_high_water_bytes",
               std::to_string(m_last_tree_pm_diagnostics.residual_worker_scratch_high_water_bytes)},
              {"domain_hierarchy_node_count", std::to_string(m_last_tree_pm_diagnostics.domain_hierarchy_node_count)},
              {"domain_cache_hits", std::to_string(m_last_tree_pm_diagnostics.domain_cache_hit_count)},
              {"domain_cache_misses", std::to_string(m_last_tree_pm_diagnostics.domain_cache_miss_count)},
              {"graph_cache_hits", std::to_string(m_last_tree_pm_diagnostics.graph_cache_hit_count)},
              {"graph_cache_misses", std::to_string(m_last_tree_pm_diagnostics.graph_cache_miss_count)},
              {"graph_cache_superset_hits", std::to_string(m_last_tree_pm_diagnostics.graph_cache_superset_hit_count)},
              {"let_candidate_peer_count", std::to_string(m_last_tree_pm_diagnostics.let_candidate_peer_count)},
              {"communicating_peer_count", std::to_string(m_last_tree_pm_diagnostics.communicating_peer_count)},
              {"let_exported_target_count", std::to_string(m_last_tree_pm_diagnostics.let_exported_target_count)},
              {"let_imported_target_count", std::to_string(m_last_tree_pm_diagnostics.let_imported_target_count)},
              {"let_wire_bytes_sent", std::to_string(m_last_tree_pm_diagnostics.let_wire_bytes_sent)},
              {"let_wire_bytes_received", std::to_string(m_last_tree_pm_diagnostics.let_wire_bytes_received)},
              {"let_high_water_bytes", std::to_string(m_last_tree_pm_diagnostics.let_high_water_bytes)},
              {"let_wire_buffer_high_water_bytes", std::to_string(
                  m_last_tree_pm_diagnostics.let_wire_buffer_high_water_bytes)},
              {"let_known_workspace_high_water_bytes", std::to_string(
                  m_last_tree_pm_diagnostics.let_known_workspace_high_water_bytes)},
              {"communication_arena_capacity_bytes", std::to_string(m_last_tree_pm_diagnostics.communication_arena_capacity_bytes)},
              {"communication_arena_logical_high_water_bytes", std::to_string(m_last_tree_pm_diagnostics.communication_arena_logical_high_water_bytes)},
              {"communication_arena_pm_density_high_water_bytes", std::to_string(m_last_tree_pm_diagnostics.communication_arena_phase_high_water_bytes[0])},
              {"communication_arena_pm_halo_high_water_bytes", std::to_string(m_last_tree_pm_diagnostics.communication_arena_phase_high_water_bytes[1])},
              {"communication_arena_pm_interpolation_high_water_bytes", std::to_string(m_last_tree_pm_diagnostics.communication_arena_phase_high_water_bytes[2])},
              {"communication_arena_tree_exchange_high_water_bytes", std::to_string(m_last_tree_pm_diagnostics.communication_arena_phase_high_water_bytes[3])},
              {"communication_arena_growth_events", std::to_string(m_last_tree_pm_diagnostics.communication_arena_growth_events)},
              {"communication_arena_unused_bytes", std::to_string(m_last_tree_pm_diagnostics.communication_arena_unused_bytes)},
              {"pm_density_routing_planned_bytes", std::to_string(m_last_tree_pm_diagnostics.pm_density_routing_planned_bytes)},
              {"pm_interpolation_routing_planned_bytes", std::to_string(m_last_tree_pm_diagnostics.pm_interpolation_routing_planned_bytes)},
              {"pm_density_routing_ceiling_bytes", std::to_string(m_last_tree_pm_diagnostics.pm_density_routing_ceiling_bytes)},
              {"pm_interpolation_routing_ceiling_bytes", std::to_string(m_last_tree_pm_diagnostics.pm_interpolation_routing_ceiling_bytes)},
              {"exported_targets_per_requested_target", formatRuntimeDouble(m_last_tree_pm_diagnostics.exported_targets_per_requested_target)},
              {"let_discovery_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_discovery_ms)},
              {"let_graph_setup_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_graph_setup_ms)},
              {"let_communication_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_communication_ms)},
              {"let_overlap_local_work_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_overlap_local_work_ms)},
              {"let_communication_wait_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_communication_wait_ms)},
              {"let_overlap_efficiency", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_overlap_efficiency)},
              {"let_remote_traversal_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_remote_traversal_ms)},
              {"incoming_request_decode_validation_ms", formatRuntimeDouble(
                  m_last_tree_pm_diagnostics.incoming_request_decode_validation_ms)},
              {"incoming_remote_target_compute_ms", formatRuntimeDouble(
                  m_last_tree_pm_diagnostics.incoming_remote_target_compute_ms)},
              {"incoming_response_encode_pack_ms", formatRuntimeDouble(
                  m_last_tree_pm_diagnostics.incoming_response_encode_pack_ms)},
              {"protocol_validation_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.protocol_validation_ms)},
              {"protocol_consensus_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.protocol_consensus_ms)},
              {"response_exchange_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.response_exchange_ms)},
              {"local_pair_evaluations", std::to_string(m_last_tree_pm_diagnostics.local_pair_evaluations)},
              {"incoming_remote_pair_evaluations",
               std::to_string(m_last_tree_pm_diagnostics.incoming_remote_pair_evaluations)},
              {"remote_pair_evaluations",
               std::to_string(m_last_tree_pm_diagnostics.incoming_remote_pair_evaluations)},
              {"total_pair_evaluations",
               std::to_string(m_last_tree_pm_diagnostics.residual_pair_evaluations)},
              {"remote_request_packet_imbalance_ratio", formatRuntimeDouble(
                  m_last_tree_pm_diagnostics.residual_remote_request_packet_imbalance_ratio)},
              {"pm_routed_density_records", std::to_string(tree_pm_profile.pm_profile.routed_density_records)},
              {"pm_routed_force_requests", std::to_string(tree_pm_profile.pm_profile.routed_force_requests)},
              {"pm_routed_potential_requests", std::to_string(tree_pm_profile.pm_profile.routed_potential_requests)},
              {"pm_routed_density_peer_count", std::to_string(tree_pm_profile.pm_profile.routed_density_peer_count)},
              {"pm_routed_force_peer_count", std::to_string(tree_pm_profile.pm_profile.routed_force_peer_count)},
              {"pm_routed_potential_peer_count", std::to_string(tree_pm_profile.pm_profile.routed_potential_peer_count)},
              {"pm_routed_mpi_bytes_sent", std::to_string(tree_pm_profile.pm_profile.routed_mpi_bytes_sent)},
              {"pm_routed_mpi_bytes_received", std::to_string(tree_pm_profile.pm_profile.routed_mpi_bytes_received)},
              {"pm_routed_send_buffer_high_water_bytes", std::to_string(
                  tree_pm_profile.pm_profile.routed_send_buffer_high_water_bytes)},
              {"pm_routed_receive_buffer_high_water_bytes", std::to_string(
                  tree_pm_profile.pm_profile.routed_receive_buffer_high_water_bytes)},
              {"pm_routed_combined_buffer_high_water_bytes", std::to_string(
                  tree_pm_profile.pm_profile.routed_combined_buffer_high_water_bytes)},
              {"pm_routed_workspace_high_water_bytes", std::to_string(
                  tree_pm_profile.pm_profile.routed_workspace_high_water_bytes)},
              {"pm_routed_mpi_wait_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.routed_mpi_wait_ms)},
          },
      });
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.zoom_force_diagnostics",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = hierarchical ? context.timeline_step.time_end_code : context.integrator_state.current_time_code,
          .scale_factor = hierarchical ? context.timeline_step.scale_factor_end : context.integrator_state.current_scale_factor,
          .message = "zoom force decomposition and contamination diagnostics",
          .payload = {
              {"force_l2_pm_global", formatRuntimeDouble(m_last_tree_pm_diagnostics.force_l2_pm_global)},
              {"force_l2_pm_zoom_correction", formatRuntimeDouble(m_last_tree_pm_diagnostics.force_l2_pm_zoom_correction)},
              {"force_l2_tree_short_range", formatRuntimeDouble(m_last_tree_pm_diagnostics.force_l2_tree_short_range)},
              {"force_l2_total", formatRuntimeDouble(m_last_tree_pm_diagnostics.force_l2_total)},
              {"zoom_high_res_source_count", std::to_string(m_last_tree_pm_diagnostics.zoom_high_res_source_count)},
              {"zoom_low_res_source_count", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_source_count)},
              {"zoom_low_res_contamination_count", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_contamination_count)},
              {"zoom_low_res_contamination_mass_code", formatRuntimeDouble(m_last_tree_pm_diagnostics.zoom_low_res_contamination_mass_code)},
              {"zoom_membership_source", m_zoom_membership_source},
          },
      });
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.health_summary",
          .severity = (gravity_health.fatal_count > 0U)
              ? core::RuntimeEventSeverity::kFatal
              : (gravity_health.warning_count > 0U)
                    ? core::RuntimeEventSeverity::kWarning
                    : core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = hierarchical ? context.timeline_step.time_end_code : context.integrator_state.current_time_code,
          .scale_factor = hierarchical ? context.timeline_step.scale_factor_end : context.integrator_state.current_scale_factor,
          .message = "gravity health summary across PM field, force, sync, zoom, and decomposition checks",
          .payload = {
              {"cheap_checks_executed", std::to_string(gravity_health.cheap_checks_executed)},
              {"heavy_checks_executed", std::to_string(gravity_health.heavy_checks_executed)},
              {"warning_count", std::to_string(gravity_health.warning_count)},
              {"fatal_count", std::to_string(gravity_health.fatal_count)},
              {"pm_field_non_finite_count", std::to_string(gravity_health.pm_field_non_finite_count)},
              {"force_non_finite_count", std::to_string(gravity_health.force_non_finite_count)},
              {"force_abnormal_count", std::to_string(gravity_health.force_abnormal_count)},
              {"illegal_sync_state_count", std::to_string(gravity_health.illegal_sync_state_count)},
              {"zoom_sanity_failure_count", std::to_string(gravity_health.zoom_sanity_failure_count)},
              {"decomposition_sanity_failure_count", std::to_string(gravity_health.decomposition_sanity_failure_count)},
              {"heavy_reference_checks_opt_in", allow_heavy_reference_checks ? "true" : "false"},
          },
      });
    }

    if (is_kick_stage && !hierarchical) {
      applyActiveKickFromFreshForce(context);
    }

    constexpr std::uint32_t no_scatter_row =
        std::numeric_limits<std::uint32_t>::max();
    const auto solved_active_indices = activeSourceIndices();
    for (std::size_t target_slot = 0; target_slot < solved_active_indices.size(); ++target_slot) {
      const std::uint32_t particle_row = targetParticleRow(target_slot);
      if (particle_row != no_scatter_row) {
        m_particle_accel_x[particle_row] = m_active_accel_x[target_slot];
        m_particle_accel_y[particle_row] = m_active_accel_y[target_slot];
        m_particle_accel_z[particle_row] = m_active_accel_z[target_slot];
        if (hierarchical) {
          // Keep TOTAL history for relative accuracy/timestep proposals. Fine
          // solves return only Tree, so add the coarse endpoint PM history;
          // this value is never mistaken for a fresh long-range force.
          if (m_tree_pm_options.short_range_only) {
            m_particle_accel_x[particle_row] += m_hierarchical_pm_x[particle_row];
            m_particle_accel_y[particle_row] += m_hierarchical_pm_y[particle_row];
            m_particle_accel_z[particle_row] += m_hierarchical_pm_z[particle_row];
          }
          m_hierarchical_force_generation[particle_row] = context.state.gravitySourceGeneration();
        }
        m_particle_force_cache_valid[particle_row] = 1U;
      }
      const std::uint32_t cell_row = targetCellRow(target_slot);
      if (cell_row != no_scatter_row) {
        m_cell_accel_x[cell_row] = m_active_accel_x[target_slot];
        m_cell_accel_y[cell_row] = m_active_accel_y[target_slot];
        m_cell_accel_z[cell_row] = m_active_accel_z[target_slot];
        m_cell_force_cache_valid[cell_row] = 1U;
      }
    }

    m_force_cache_particle_index_generation = context.state.particleIndexGeneration();
    m_force_cache_cell_index_generation = context.state.cellIndexGeneration();
    m_force_cache_gas_identity_generation = context.state.gasCellIdentityGeneration();
    m_force_cache_source_generation = context.state.gravitySourceGeneration();
    m_force_cache_valid = true;
    const std::uint32_t local_rank =
        static_cast<std::uint32_t>(std::max(mpi_context.worldRank(), 0));
    const std::uint32_t gas_species_tag =
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    for (std::size_t particle_index = 0; !hierarchical && particle_index < particle_count;
         ++particle_index) {
      if (context.state.particleSpeciesTag(particle_index) ==
          gas_species_tag) {
        continue;
      }
      if (context.state.particleOwningRank(particle_index) ==
              local_rank &&
          m_particle_force_cache_valid[particle_index] == 0U) {
        m_force_cache_valid = false;
        break;
      }
    }
    if (m_force_cache_valid &&
        m_source_representation == gravity::GravitySourceRepresentation::kMaterializedGeneric) {
      for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
        if (m_owned_leaf_cell_mask[cell_index] != 0U &&
            m_cell_force_cache_valid[cell_index] == 0U) {
          m_force_cache_valid = false;
          break;
        }
      }
    }

    if (is_kick_stage && hierarchical) applyHierarchicalKick(context);

    if (gravity_fallback_reservation.valid()) {
      gravity_fallback_reservation.release();
    }
    if (gravity_phase_reservation.valid()) {
      const std::uint64_t gravity_baseline_after =
          core::memoryReportBaselineOwnedBytes(memoryReport());
      const std::uint64_t non_gravity_baseline =
          process_baseline_before - gravity_baseline_before;
      const std::uint64_t reconciled_process_baseline =
          core::checkedMemoryBytesAdd(
              non_gravity_baseline, gravity_baseline_after,
              "gravity retained-capacity baseline reconciliation");
      gravity_phase_reservation.reconcileBaselineOwnedAndRelease(
          reconciled_process_baseline);
    }
    // Borrowed spans may alias the step workspace. Keep only the explicit
    // representation/generation policy after the solve; rebuild the span at the
    // next force evaluation rather than retaining a pointer across phase reset.
    if (m_source_representation ==
        gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      m_active_source_indices = {};
    }

  }

 private:
  [[nodiscard]] static double kickFactorForStage(const core::StepContext& context) {
    if (context.stage == core::IntegrationStage::kGravityKickPre) {
      return context.timeline_step.first_kick_factor_code;
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost) {
      return context.timeline_step.second_kick_factor_code;
    }
    throw std::invalid_argument("gravity kick factor requested outside KDK kick stage");
  }

  // Drag-free source kick for magnetized gas (MHD owns the expansion source).
  [[nodiscard]] static double sourceKickFactorForStage(const core::StepContext& context) {
    if (context.stage == core::IntegrationStage::kGravityKickPre) {
      return context.timeline_step.first_source_kick_factor_code;
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost) {
      return context.timeline_step.second_source_kick_factor_code;
    }
    throw std::invalid_argument("source kick factor requested outside KDK kick stage");
  }

  [[nodiscard]] static double hubbleDragFactorForStage(const core::StepContext& context) {
    if (context.stage == core::IntegrationStage::kGravityKickPre) {
      return context.timeline_step.first_hubble_drag_factor;
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost) {
      return context.timeline_step.second_hubble_drag_factor;
    }
    throw std::invalid_argument("Hubble drag factor requested outside KDK kick stage");
  }

  static void applyPeculiarVelocityKick(
      core::SimulationState& state,
      std::uint32_t particle_index,
      double accel_x_code,
      double accel_y_code,
      double accel_z_code,
      double kick_factor_code,
      double hubble_drag_factor) {
    state.particles.velocity_x_peculiar[particle_index] =
        hubble_drag_factor * state.particles.velocity_x_peculiar[particle_index] + kick_factor_code * accel_x_code;
    state.particles.velocity_y_peculiar[particle_index] =
        hubble_drag_factor * state.particles.velocity_y_peculiar[particle_index] + kick_factor_code * accel_y_code;
    state.particles.velocity_z_peculiar[particle_index] =
        hubble_drag_factor * state.particles.velocity_z_peculiar[particle_index] + kick_factor_code * accel_z_code;
  }

  static void applyGasCellKdkKick(
      core::SimulationState& state,
      std::uint32_t cell_index,
      double accel_x_code,
      double accel_y_code,
      double accel_z_code,
      double kick_factor_code,
      double hubble_drag_factor) {
    if (cell_index >= state.gas_cells.size()) {
      throw std::out_of_range("gas KDK kick cell index is outside GasCellSidecar");
    }
    if (!std::isfinite(hubble_drag_factor) || hubble_drag_factor <= 0.0) {
      throw std::invalid_argument("gas KDK kick requires finite positive Hubble drag factor");
    }
    state.gas_cells.velocity_x_peculiar[cell_index] =
        hubble_drag_factor * state.gas_cells.velocity_x_peculiar[cell_index] +
        kick_factor_code * accel_x_code;
    state.gas_cells.velocity_y_peculiar[cell_index] =
        hubble_drag_factor * state.gas_cells.velocity_y_peculiar[cell_index] +
        kick_factor_code * accel_y_code;
    state.gas_cells.velocity_z_peculiar[cell_index] =
        hubble_drag_factor * state.gas_cells.velocity_z_peculiar[cell_index] +
        kick_factor_code * accel_z_code;

    // For gamma=5/3, homogeneous adiabatic expansion gives
    // specific internal energy u_int propto a^-2. The timeline drag factor is
    // a_begin/a_end for this half-kick, so thermal state scales as drag^2.
    // Density/mass are unchanged by a source kick; the finite-volume transport
    // stage owns those conserved quantities.
    constexpr double gamma_minus_one = 2.0 / 3.0;
    const double thermal_drag =
        std::pow(hubble_drag_factor, 3.0 * gamma_minus_one);
    state.gas_cells.internal_energy_code[cell_index] *= thermal_drag;
    state.gas_cells.pressure_code[cell_index] *= thermal_drag;
    state.gas_cells.temperature_code[cell_index] *= thermal_drag;
    state.gas_cells.sound_speed_code[cell_index] *= std::sqrt(thermal_drag);
  }

  // Magnetized-gas KDK source kick (CHUI-MHD-REF 4.3/8.5 expansion ownership,
  // pass-02 reconciliation R-EXP): the MHD ideal stage applies -H m_c and
  // -H (2K_c + 3p_c + M_c) at its SSPRK2 stage endpoints, so this kick adds
  // only the gravitational acceleration du = A * integral(dt/a^2). No Hubble
  // drag and no adiabatic thermal scaling are applied here (that would be a
  // second expansion owner). Gravitational work enters through the kinetic
  // energy of the updated velocity when the MHD total energy is rebuilt from
  // the gas lanes at the next kHydroUpdate.
  static void applyMagnetizedGasCellSourceKick(
      core::SimulationState& state,
      std::uint32_t cell_index,
      double accel_x_code,
      double accel_y_code,
      double accel_z_code,
      double source_kick_factor_code) {
    if (cell_index >= state.gas_cells.size()) {
      throw std::out_of_range("magnetized gas kick cell index is outside GasCellSidecar");
    }
    if (!std::isfinite(source_kick_factor_code) || source_kick_factor_code < 0.0) {
      throw std::invalid_argument("magnetized gas kick requires a finite source kick factor");
    }
    state.gas_cells.velocity_x_peculiar[cell_index] += source_kick_factor_code * accel_x_code;
    state.gas_cells.velocity_y_peculiar[cell_index] += source_kick_factor_code * accel_y_code;
    state.gas_cells.velocity_z_peculiar[cell_index] += source_kick_factor_code * accel_z_code;
  }

  void applyGasKick(core::StepContext& context, std::uint32_t cell_index, double ax, double ay,
                    double az, double kick_factor, double hubble_drag) const {
    if (m_config.mhd.enabled) {
      applyMagnetizedGasCellSourceKick(context.state, cell_index, ax, ay, az,
                                       sourceKickFactorForStage(context));
      return;
    }
    applyGasCellKdkKick(context.state, cell_index, ax, ay, az, kick_factor, hubble_drag);
  }

  void applyActiveKickFromFreshForce(core::StepContext& context) {
    const double kick_factor = kickFactorForStage(context);
    const double hubble_drag = hubbleDragFactorForStage(context);
    for (const std::uint32_t particle_index :
         context.active_set.particle_indices) {
      const std::optional<std::size_t> active_slot =
          activeSlotForParticle(particle_index);
      if (!active_slot.has_value()) {
        continue;
      }
      applyPeculiarVelocityKick(
          context.state,
          particle_index,
          m_active_accel_x[*active_slot],
          m_active_accel_y[*active_slot],
          m_active_accel_z[*active_slot],
          kick_factor,
          hubble_drag);
    }
    if (m_source_representation ==
        gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      if (!context.active_set.cell_indices.empty()) {
        throw std::logic_error(
            "borrowed homogeneous DMO fresh-force path unexpectedly received gas-cell targets");
      }
      return;
    }
    for (const std::uint32_t cell_index : context.active_set.cell_indices) {
      if (cell_index >= m_active_slot_by_cell.size()) {
        throw std::out_of_range("gas KDK fresh-force cell index is outside the active slot map");
      }
      const int active_slot = m_active_slot_by_cell[cell_index];
      if (active_slot < 0) {
        continue;
      }
      applyGasKick(
          context,
          cell_index,
          m_active_accel_x[static_cast<std::size_t>(active_slot)],
          m_active_accel_y[static_cast<std::size_t>(active_slot)],
          m_active_accel_z[static_cast<std::size_t>(active_slot)],
          kick_factor,
          hubble_drag);
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost &&
        !context.active_set.cell_indices.empty()) {
      internal::synchronizeParentParticleCompatibilityMirrors(
          context.state,
          static_cast<std::uint32_t>(std::max(m_services.mpi_context.worldRank(), 0)),
          "gravity post-kick gas compatibility mirror");
    }
  }

  void applyHierarchicalKick(core::StepContext& context) {
    if (!context.hierarchical_kdk.enabled || context.particle_scheduler == nullptr ||
        m_force_cache_particle_index_generation != context.state.particleIndexGeneration() ||
        m_hierarchical_force_generation.size() != context.state.particles.size() ||
        m_hierarchical_tree_x.size() != context.state.particles.size() ||
        m_hierarchical_tree_y.size() != context.state.particles.size() ||
        m_hierarchical_tree_z.size() != context.state.particles.size() ||
        m_hierarchical_pm_x.size() != context.state.particles.size() ||
        m_hierarchical_pm_y.size() != context.state.particles.size() ||
        m_hierarchical_pm_z.size() != context.state.particles.size() ||
        context.hierarchical_kdk.evaluation_tick != context.particle_scheduler->currentTick()) {
      throw std::logic_error("hierarchical kick rejected incompatible row/scheduler epoch");
    }
    for (const auto row : context.active_set.particle_indices) {
      if (row >= m_hierarchical_force_generation.size() ||
          m_hierarchical_force_generation[row] != context.state.gravitySourceGeneration() ||
          m_particle_force_cache_valid[row] == 0U) {
        throw std::logic_error("hierarchical kick rejected a Tree force outside its exact source epoch");
      }
      const auto bin = context.particle_scheduler->binIndex(row);
      const double tree_factor = context.hierarchical_kdk.tree_kick_factor_code.at(bin);
      const double pm_factor = context.hierarchical_kdk.pm_kick_factor_code;
      if (pm_factor != 0.0 &&
          (m_hierarchical_pm_source_generation != context.state.gravitySourceGeneration() ||
           m_hierarchical_pm_field_version != context.integrator_state.pm_sync_state.fieldVersion())) {
        throw std::logic_error("hierarchical PM kick rejected a component outside its synchronized field epoch");
      }
      if (!std::isfinite(tree_factor) || tree_factor <= 0.0 ||
          !std::isfinite(pm_factor) || pm_factor < 0.0) {
        throw std::logic_error("hierarchical kick factor is not finite with legal sign");
      }
      // PM enters only synchronized coarse endpoint half-kicks. Tree was
      // accumulated independently; never recover it by PM/total subtraction.
      applyPeculiarVelocityKick(context.state, row,
          tree_factor * m_hierarchical_tree_x[row] + pm_factor * m_hierarchical_pm_x[row],
          tree_factor * m_hierarchical_tree_y[row] + pm_factor * m_hierarchical_pm_y[row],
          tree_factor * m_hierarchical_tree_z[row] + pm_factor * m_hierarchical_pm_z[row],
          1.0, 1.0);
    }
  }

  void applyCachedKick(core::StepContext& context) {
    if (!m_force_cache_valid) {
      throw std::runtime_error("gravity kick requested before a coherent force-refresh boundary");
    }
    if (m_force_cache_particle_index_generation != context.state.particleIndexGeneration() ||
        m_force_cache_cell_index_generation != context.state.cellIndexGeneration() ||
        m_force_cache_gas_identity_generation != context.state.gasCellIdentityGeneration() ||
        m_force_cache_source_generation != context.state.gravitySourceGeneration()) {
      throw std::runtime_error(
          "gravity cached kick rejected because particle/cell source identity changed since force refresh");
    }
    const std::size_t particle_count = context.state.particles.size();
    if (m_particle_accel_x.size() != particle_count ||
        m_particle_accel_y.size() != particle_count ||
        m_particle_accel_z.size() != particle_count ||
        m_particle_force_cache_valid.size() != particle_count) {
      throw std::runtime_error(
          "gravity cached kick rejected an incomplete particle acceleration cache");
    }
    const std::uint32_t local_rank =
        static_cast<std::uint32_t>(std::max(m_services.mpi_context.worldRank(), 0));
    const std::uint32_t gas_species_tag =
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    const double kick_factor = kickFactorForStage(context);
    const double hubble_drag = hubbleDragFactorForStage(context);
    for (const std::uint32_t particle_index : context.active_set.particle_indices) {
      if (particle_index >= particle_count) {
        throw std::out_of_range("gravity cached kick particle index out of range");
      }
      if (context.state.particleOwningRank(particle_index) != local_rank ||
          context.state.particleSpeciesTag(particle_index) == gas_species_tag) {
        continue;
      }
      if (m_particle_force_cache_valid[particle_index] == 0U) {
        throw std::runtime_error(
            "gravity cached kick rejected an invalid per-particle force row");
      }
      applyPeculiarVelocityKick(
          context.state,
          particle_index,
          m_particle_accel_x[particle_index],
          m_particle_accel_y[particle_index],
          m_particle_accel_z[particle_index],
          kick_factor,
          hubble_drag);
    }
    if (m_cell_accel_x.size() != context.state.cells.size() ||
        m_cell_accel_y.size() != context.state.cells.size() ||
        m_cell_accel_z.size() != context.state.cells.size() ||
        m_cell_force_cache_valid.size() != context.state.cells.size()) {
      throw std::runtime_error(
          "gravity cached kick rejected an incomplete gas-cell acceleration cache");
    }
    const auto authoritative_rows = internal::selectAuthoritativeGravitySourceRows(
        context.state,
        local_rank,
        static_cast<std::uint32_t>(std::max(m_services.mpi_context.worldSize(), 1)),
        "gas KDK cached kick");
    std::vector<std::uint8_t> is_owned_leaf(context.state.cells.size(), 0U);
    for (const std::uint32_t cell_row : authoritative_rows.gas_cell_rows) {
      is_owned_leaf[cell_row] = 1U;
    }
    for (const std::uint32_t cell_index : context.active_set.cell_indices) {
      if (cell_index >= context.state.cells.size()) {
        throw std::out_of_range("gas KDK cached-force cell index is outside CellSoa");
      }
      if (is_owned_leaf[cell_index] == 0U) {
        continue;
      }
      if (m_cell_force_cache_valid[cell_index] == 0U) {
        throw std::runtime_error(
            "gas KDK cached kick rejected an invalid per-cell force row");
      }
      applyGasKick(
          context,
          cell_index,
          m_cell_accel_x[cell_index],
          m_cell_accel_y[cell_index],
          m_cell_accel_z[cell_index],
          kick_factor,
          hubble_drag);
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost &&
        !context.active_set.cell_indices.empty()) {
      internal::synchronizeParentParticleCompatibilityMirrors(
          context.state, local_rank, "gravity post-kick gas compatibility mirror");
    }
  }

  void emitGravityEvent(
      core::StepContext& context,
      core::RuntimeEventSeverity severity,
      std::string message,
      std::unordered_map<std::string, std::string> payload = {}) const {
    if (context.profiler_session == nullptr) {
      return;
    }
    context.profiler_session->recordEvent(core::RuntimeEvent{
        .event_kind = "gravity.health_check",
        .severity = severity,
        .subsystem = "gravity.treepm",
        .step_index = context.integrator_state.step_index,
        .simulation_time_code = context.integrator_state.current_time_code,
        .scale_factor = context.integrator_state.current_scale_factor,
        .message = std::move(message),
        .payload = std::move(payload),
    });
  }

  [[nodiscard]] GravityHealthSummary validateGravityHealth(
      core::StepContext& context,
      const PmLongRangeKickDecision& decision,
      bool allow_heavy_reference_checks) {
    GravityHealthSummary summary;
    auto failFatal = [&](std::string message, std::unordered_map<std::string, std::string> payload = {}) {
      ++summary.fatal_count;
      emitGravityEvent(context, core::RuntimeEventSeverity::kFatal, message, payload);
      throw std::runtime_error("fatal gravity-state check failed: " + message);
    };
    auto addWarning = [&](std::string message, std::unordered_map<std::string, std::string> payload = {}) {
      ++summary.warning_count;
      emitGravityEvent(context, core::RuntimeEventSeverity::kWarning, message, payload);
    };

    ++summary.cheap_checks_executed;
    const auto health_active_indices = activeSourceIndices();
    const std::size_t health_source_count = sourceMass(context.state).size();
    if (m_source_representation == gravity::GravitySourceRepresentation::kMaterializedGeneric &&
        (health_active_indices.size() != m_active_target_particle_row.size() ||
         health_active_indices.size() != m_active_target_cell_row.size())) {
      ++summary.decomposition_sanity_failure_count;
      failFatal("local gravity target vectors diverged after compact-view rebuild");
    }
    if (health_active_indices.size() > health_source_count) {
      ++summary.decomposition_sanity_failure_count;
      failFatal("active local subset exceeds local source population");
    }

    ++summary.cheap_checks_executed;
    if (!std::isfinite(m_last_tree_pm_diagnostics.force_l2_pm_global) ||
        !std::isfinite(m_last_tree_pm_diagnostics.force_l2_pm_zoom_correction) ||
        !std::isfinite(m_last_tree_pm_diagnostics.force_l2_tree_short_range) ||
        !std::isfinite(m_last_tree_pm_diagnostics.force_l2_total)) {
      ++summary.pm_field_non_finite_count;
      failFatal("PM/tree force norm diagnostics contain NaN or Inf");
    }

    ++summary.cheap_checks_executed;
    for (std::size_t i = 0; i < m_active_accel_x.size(); ++i) {
      if (!std::isfinite(m_active_accel_x[i]) ||
          !std::isfinite(m_active_accel_y[i]) ||
          !std::isfinite(m_active_accel_z[i])) {
        ++summary.force_non_finite_count;
        failFatal(
            "active-set acceleration contains NaN or Inf",
            {{"active_slot", std::to_string(i)}});
      }
    }

    ++summary.cheap_checks_executed;
    if (context.integrator_state.pm_sync_state.lastRefreshOpportunity() > context.integrator_state.pm_sync_state.gravityKickOpportunity() ||
        context.integrator_state.pm_sync_state.fieldVersion() < m_last_committed_field_version) {
      ++summary.illegal_sync_state_count;
      failFatal("gravity sync-state regressed across kick opportunities");
    }
    if (decision.refresh_long_range_field) {
      if (decision.field_version != context.integrator_state.pm_sync_state.fieldVersion() + 1U) {
        ++summary.illegal_sync_state_count;
        failFatal("refresh decision must increment field_version by exactly one");
      }
    } else if (decision.field_version != context.integrator_state.pm_sync_state.fieldVersion()) {
      ++summary.illegal_sync_state_count;
      failFatal("reuse decision must not mutate field_version");
    }

    ++summary.cheap_checks_executed;
    if (m_mode_policy.zoom_region_expected &&
        m_tree_pm_options.zoom_region_radius_comoving <= 0.0) {
      ++summary.zoom_sanity_failure_count;
      failFatal("zoom mode requires a strictly positive zoom region radius");
    }
    if (!m_mode_policy.zoom_region_expected &&
        m_tree_pm_options.enable_zoom_long_range_correction &&
        m_last_tree_pm_diagnostics.zoom_high_res_source_count == 0U &&
        m_last_tree_pm_diagnostics.zoom_low_res_source_count == 0U) {
      ++summary.zoom_sanity_failure_count;
      addWarning("zoom long-range correction enabled without detected zoom population");
    }

    if (allow_heavy_reference_checks) {
      ++summary.heavy_checks_executed;
      const double total_force = m_last_tree_pm_diagnostics.force_l2_total;
      const double pm_component = m_last_tree_pm_diagnostics.force_l2_pm_global;
      const double tree_component = m_last_tree_pm_diagnostics.force_l2_tree_short_range;
      if (total_force > 0.0 && std::isfinite(total_force)) {
        const double component_sum = pm_component + tree_component;
        const double relative_gap = std::abs(component_sum - total_force) / total_force;
        if (relative_gap > 1.0) {
          ++summary.force_abnormal_count;
          addWarning(
              "force decomposition appears abnormal under heavy reference check",
              {{"relative_gap", std::to_string(relative_gap)}});
        }
      }
      ++summary.heavy_checks_executed;
      if (m_mode_policy.zoom_region_expected &&
          m_last_tree_pm_diagnostics.zoom_low_res_contamination_count > 0U) {
        ++summary.force_abnormal_count;
        addWarning(
            "zoom contamination detected in low-resolution source particles",
            {{"zoom_low_res_contamination_count",
              std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_contamination_count)}});
      }
    }

    return summary;
  }

  void ensureZoomMembershipLoaded(const core::SimulationState& state) {
    if (m_zoom_membership_loaded) {
      return;
    }
    m_zoom_membership_loaded = true;
    if (!m_mode_policy.zoom_region_expected || m_config.mode.zoom_region_file.empty()) {
      m_zoom_membership_source = m_mode_policy.zoom_region_expected ? "configured_spherical_region" : "disabled";
      return;
    }
    const std::filesystem::path path = m_zoom_region_path.empty()
        ? std::filesystem::path(m_config.mode.zoom_region_file)
        : m_zoom_region_path;
    if (hasHdf5Extension(path)) {
#if COSMOSIM_ENABLE_HDF5
      m_zoom_high_res_particle_ids = loadZoomParticleIdsFromHdf5(path);
      m_zoom_membership_source = "particle_id_file_hdf5";
#else
      throw std::runtime_error("zoom_region_file requires HDF5 support in this build");
#endif
    } else {
      m_zoom_high_res_particle_ids = loadZoomParticleIdsFromText(path);
      m_zoom_membership_source = "particle_id_file_text";
    }
    for (const std::uint64_t id : state.particle_sidecar.particle_id) {
      (void)id;
    }
  }

  template <typename T>
  static void releaseVectorCapacity(std::vector<T>& values) {
    std::vector<T>().swap(values);
  }

  [[nodiscard]] static bool isIdentityParticleSpan(
      std::span<const std::uint32_t> rows,
      std::size_t particle_count) noexcept {
    if (rows.size() != particle_count) {
      return false;
    }
    for (std::size_t i = 0; i < particle_count; ++i) {
      if (rows[i] != i) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] bool borrowedHomogeneousDmoBaseEligible(
      const core::SimulationState& state) const {
    const bool cacheable_compact_state =
        state.hasHomogeneousDmoMetadata() &&
        state.particle_sidecar.has_gravity_softening_override.empty();
    if (cacheable_compact_state && m_borrowed_eligibility_certificate_valid &&
        m_borrowed_eligibility_particle_index_generation == state.particleIndexGeneration() &&
        m_borrowed_eligibility_gravity_source_generation == state.gravitySourceGeneration()) {
      return m_borrowed_eligibility_result;
    }
    const auto publish_certificate = [&](bool eligible) {
      if (cacheable_compact_state) {
        m_borrowed_eligibility_certificate_valid = true;
        m_borrowed_eligibility_particle_index_generation = state.particleIndexGeneration();
        m_borrowed_eligibility_gravity_source_generation = state.gravitySourceGeneration();
        m_borrowed_eligibility_result = eligible;
      }
      return eligible;
    };
    if (m_tree_pm_options.enable_zoom_long_range_correction ||
        state.cells.size() != 0U) {
      return publish_certificate(false);
    }
    const std::size_t particle_count = state.particles.size();
    if (particle_count > static_cast<std::size_t>(
            std::numeric_limits<std::uint32_t>::max())) {
      return publish_certificate(false);
    }
    if (!state.particles.isConsistent() || !state.particle_sidecar.isConsistent()) {
      return publish_certificate(false);
    }
    if (!state.hasHomogeneousDmoMetadata() &&
        (state.particle_sidecar.species_tag.size() != particle_count ||
         state.particle_sidecar.owning_rank.size() != particle_count)) {
      return publish_certificate(false);
    }
    if (!state.particle_sidecar.has_gravity_softening_override.empty() &&
        state.particle_sidecar.has_gravity_softening_override.size() != particle_count) {
      return publish_certificate(false);
    }

    const std::uint32_t local_rank =
        static_cast<std::uint32_t>(std::max(m_runtime_topology.world_rank, 0));
    const std::uint32_t dm_species_tag =
        static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
    const bool periodic_sources =
        m_tree_pm_options.pm_options.boundary_condition ==
        gravity::PmBoundaryCondition::kPeriodic;
    const double box_size_x = m_config.cosmology.box_size_x_mpc_comoving;
    const double box_size_y = m_config.cosmology.box_size_y_mpc_comoving;
    const double box_size_z = m_config.cosmology.box_size_z_mpc_comoving;
    const bool compact_policy_proves_species_and_owner = state.hasHomogeneousDmoMetadata();
    for (std::size_t row = 0; row < particle_count; ++row) {
      if (!compact_policy_proves_species_and_owner &&
          (state.particleOwningRank(row) != local_rank ||
           state.particleSpeciesTag(row) != dm_species_tag)) {
        return publish_certificate(false);
      }
      if (!state.particle_sidecar.has_gravity_softening_override.empty() &&
          state.particle_sidecar.has_gravity_softening_override[row] != 0U) {
        return publish_certificate(false);
      }
      const double x = state.particles.position_x_comoving[row];
      const double y = state.particles.position_y_comoving[row];
      const double z = state.particles.position_z_comoving[row];
      const double mass = state.particles.mass_code[row];
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
          !std::isfinite(mass) || mass < 0.0) {
        throw std::runtime_error(
            "borrowed homogeneous DMO gravity source has non-finite coordinates or negative/non-finite mass");
      }
      if (periodic_sources &&
          (x < 0.0 || x >= box_size_x || y < 0.0 || y >= box_size_y ||
           z < 0.0 || z >= box_size_z)) {
        return publish_certificate(false);
      }
    }

    const double dm_softening = m_tree_pm_species_softening.enabled
        ? m_tree_pm_species_softening.epsilon_comoving_by_species[
              static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)]
        : m_tree_pm_options.tree_options.softening.epsilon_comoving;
    constexpr double max_certified_softening_to_split_ratio = 0.20;
    if (!std::isfinite(dm_softening) || dm_softening < 0.0 ||
        dm_softening > max_certified_softening_to_split_ratio *
            m_tree_pm_options.split_policy.split_scale_comoving) {
      return publish_certificate(false);
    }
    return publish_certificate(true);
  }

  [[nodiscard]] bool borrowedHomogeneousDmoEligible(
      const core::StepContext& context,
      std::span<const std::uint32_t> active_particles,
      bool prediction_requires_inactive_coordinates) const {
    return !prediction_requires_inactive_coordinates &&
        borrowedHomogeneousDmoBaseEligible(context.state) &&
        (context.hierarchical_kdk.enabled ||
         isIdentityParticleSpan(active_particles, context.state.particles.size()));
  }

  void releaseGenericGravityRepresentationCapacity() {
    releaseVectorCapacity(m_active_slot_by_particle);
    releaseVectorCapacity(m_active_slot_by_cell);
    releaseVectorCapacity(m_owned_local_index_by_global);
    releaseVectorCapacity(m_owned_local_index_by_cell);
    releaseVectorCapacity(m_owned_leaf_cell_mask);
    releaseVectorCapacity(m_local_source_x);
    releaseVectorCapacity(m_local_source_y);
    releaseVectorCapacity(m_local_source_z);
    releaseVectorCapacity(m_local_source_mass);
    releaseVectorCapacity(m_local_source_species_tag);
    releaseVectorCapacity(m_local_source_softening_comoving);
    releaseVectorCapacity(m_local_source_softening_override_mask);
    releaseVectorCapacity(m_local_source_particle_row);
    releaseVectorCapacity(m_local_source_cell_row);
    releaseVectorCapacity(m_active_target_particle_row);
    releaseVectorCapacity(m_active_target_cell_row);
    releaseVectorCapacity(m_local_active_indices);
    releaseVectorCapacity(m_local_active_global_indices);
    releaseVectorCapacity(m_force_refresh_particle_indices);
    releaseVectorCapacity(m_source_is_high_res);
    releaseVectorCapacity(m_active_is_high_res);
  }

  [[nodiscard]] std::span<const double> sourcePosX(
      const core::SimulationState& state) const noexcept {
    return m_source_representation ==
            gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
        ? std::span<const double>(state.particles.position_x_comoving.data(),
                                  state.particles.position_x_comoving.size())
        : std::span<const double>(m_local_source_x.data(), m_local_source_x.size());
  }
  [[nodiscard]] std::span<const double> sourcePosY(
      const core::SimulationState& state) const noexcept {
    return m_source_representation ==
            gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
        ? std::span<const double>(state.particles.position_y_comoving.data(),
                                  state.particles.position_y_comoving.size())
        : std::span<const double>(m_local_source_y.data(), m_local_source_y.size());
  }
  [[nodiscard]] std::span<const double> sourcePosZ(
      const core::SimulationState& state) const noexcept {
    return m_source_representation ==
            gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
        ? std::span<const double>(state.particles.position_z_comoving.data(),
                                  state.particles.position_z_comoving.size())
        : std::span<const double>(m_local_source_z.data(), m_local_source_z.size());
  }
  [[nodiscard]] std::span<const double> sourceMass(
      const core::SimulationState& state) const noexcept {
    return m_source_representation ==
            gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
        ? std::span<const double>(state.particles.mass_code.data(),
                                  state.particles.mass_code.size())
        : std::span<const double>(m_local_source_mass.data(), m_local_source_mass.size());
  }
  [[nodiscard]] std::span<const std::uint32_t> sourceSpecies() const noexcept {
    return m_source_representation ==
            gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo
        ? std::span<const std::uint32_t>{}
        : std::span<const std::uint32_t>(m_local_source_species_tag.data(),
                                         m_local_source_species_tag.size());
  }
  [[nodiscard]] std::span<const std::uint32_t> activeSourceIndices() const noexcept {
    return m_active_source_indices;
  }
  [[nodiscard]] std::uint32_t targetParticleRow(std::size_t target_slot) const {
    if (m_source_representation ==
        gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      return m_active_source_indices[target_slot];
    }
    return m_active_target_particle_row.at(target_slot);
  }
  [[nodiscard]] std::uint32_t targetCellRow(std::size_t target_slot) const {
    if (m_source_representation ==
        gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      return std::numeric_limits<std::uint32_t>::max();
    }
    return m_active_target_cell_row.at(target_slot);
  }

  [[nodiscard]] std::optional<std::size_t> activeSlotForParticle(
      std::uint32_t particle_row) const {
    if (m_source_representation ==
        gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      if (particle_row >= m_active_accel_x.size()) {
        throw std::out_of_range(
            "borrowed homogeneous DMO particle row exceeds active acceleration extent");
      }
      return static_cast<std::size_t>(particle_row);
    }
    if (particle_row >= m_active_slot_by_particle.size()) {
      throw std::out_of_range(
          "gravity kick particle index is outside the fresh-force slot map");
    }
    const int active_slot = m_active_slot_by_particle[particle_row];
    if (active_slot < 0) {
      return std::nullopt;
    }
    return static_cast<std::size_t>(active_slot);
  }

  enum class SourcePredictionEpoch : std::uint8_t {
    kNone = 0,
    kStepBegin = 1,
    kStepEnd = 2,
  };

  void rebuildOwnedParticleCompactView(
      const core::StepContext& context,
      std::span<const std::uint32_t> active_particles,
      SourcePredictionEpoch prediction_epoch = SourcePredictionEpoch::kNone,
      std::span<const std::uint32_t> prediction_current_particles = {}) {
    const core::SimulationState& state = context.state;
    const std::size_t particle_count = state.particles.size();
    const std::size_t cell_count = state.cells.size();
    const bool prediction_requires_inactive_coordinates =
        prediction_epoch != SourcePredictionEpoch::kNone &&
        !isIdentityParticleSpan(prediction_current_particles, particle_count);
    m_active_source_indices = {};
    if (borrowedHomogeneousDmoEligible(
            context, active_particles, prediction_requires_inactive_coordinates)) {
      if (context.workspace == nullptr ||
          context.workspace->gravity_particle_index_scratch.size() != particle_count) {
        throw std::logic_error(
            "borrowed homogeneous DMO gravity requires the admitted all-particle index lane");
      }
      releaseGenericGravityRepresentationCapacity();
      m_source_representation =
          gravity::GravitySourceRepresentation::kBorrowedHomogeneousDmo;
      m_borrowed_particle_index_generation = state.particleIndexGeneration();
      m_borrowed_gravity_source_generation = state.gravitySourceGeneration();
      m_active_source_indices = context.hierarchical_kdk.enabled
          ? active_particles : std::span<const std::uint32_t>(context.workspace->gravity_particle_index_scratch);
      m_source_predicted_inactive_count = 0U;
      return;
    }

    m_source_representation =
        gravity::GravitySourceRepresentation::kMaterializedGeneric;
    m_borrowed_particle_index_generation = 0U;
    m_borrowed_gravity_source_generation = 0U;
    state.requireGasCellIdentityMapCoversDenseRows(
        "gravity authoritative gas source rebuild");
    const std::uint32_t world_rank =
        static_cast<std::uint32_t>(m_runtime_topology.world_rank);
    constexpr std::uint32_t no_row = std::numeric_limits<std::uint32_t>::max();

    m_owned_local_index_by_global.assign(particle_count, -1);
    m_owned_local_index_by_cell.assign(cell_count, -1);
    m_owned_leaf_cell_mask.assign(cell_count, 0U);
    m_local_source_x.clear();
    m_local_source_y.clear();
    m_local_source_z.clear();
    m_local_source_mass.clear();
    m_local_source_species_tag.clear();
    m_local_source_softening_comoving.clear();
    m_local_source_softening_override_mask.clear();
    m_local_source_particle_row.clear();
    m_local_source_cell_row.clear();
    m_source_predicted_inactive_count = 0;

    std::vector<std::uint8_t> active_mask;
    if (prediction_epoch != SourcePredictionEpoch::kNone) {
      if (!state.hasHomogeneousDmoMetadata() &&
          (state.particle_sidecar.last_drift_time_code.size() != particle_count ||
           state.particle_sidecar.last_drift_scale_factor.size() != particle_count)) {
        throw std::runtime_error(
            "PM source prediction requires logical particle drift epochs");
      }
      active_mask.assign(particle_count, 0U);
      for (const std::uint32_t global_index : prediction_current_particles) {
        if (global_index >= particle_count) {
          throw std::out_of_range(
              "gravity callback active particle index out of range");
        }
        active_mask[global_index] = 1U;
      }
    }

    const bool periodic_sources =
        m_tree_pm_options.pm_options.boundary_condition ==
        gravity::PmBoundaryCondition::kPeriodic;
    const double box_size_x = m_config.cosmology.box_size_x_mpc_comoving;
    const double box_size_y = m_config.cosmology.box_size_y_mpc_comoving;
    const double box_size_z = m_config.cosmology.box_size_z_mpc_comoving;
    const bool has_softening_values =
        !state.particle_sidecar.gravity_softening_comoving.empty();
    const bool has_softening_masks =
        !state.particle_sidecar.has_gravity_softening_override.empty();
    const std::uint32_t gas_species_tag =
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    const auto authoritative_source_rows =
        internal::selectAuthoritativeGravitySourceRows(
            state,
            world_rank,
            static_cast<std::uint32_t>(std::max(m_runtime_topology.world_size, 1)),
            "gravity authoritative source selection");
    const std::size_t source_capacity = core::checkedSizeAdd(
        authoritative_source_rows.particle_rows.size(),
        authoritative_source_rows.gas_cell_rows.size(),
        "gravity authoritative source capacity");
    m_local_source_x.reserve(source_capacity);
    m_local_source_y.reserve(source_capacity);
    m_local_source_z.reserve(source_capacity);
    m_local_source_mass.reserve(source_capacity);
    m_local_source_species_tag.reserve(source_capacity);
    m_local_source_particle_row.reserve(source_capacity);
    m_local_source_cell_row.reserve(source_capacity);
    if (has_softening_values) {
      m_local_source_softening_comoving.reserve(source_capacity);
    }
    if (has_softening_masks) {
      m_local_source_softening_override_mask.reserve(source_capacity);
    }
    const double gas_softening_comoving = m_tree_pm_species_softening.enabled
        ? m_tree_pm_species_softening.epsilon_comoving_by_species[
              static_cast<std::size_t>(core::ParticleSpecies::kGas)]
        : m_tree_pm_options.tree_options.softening.epsilon_comoving;
    constexpr double max_certified_softening_to_split_ratio = 0.20;
    const double split_scale_comoving =
        m_tree_pm_options.split_policy.split_scale_comoving;

    const auto validateCertifiedTreePmSoftening = [&] (
        std::uint32_t species_tag,
        double softening_value,
        std::uint8_t softening_override) {
      double effective_softening =
          m_tree_pm_options.tree_options.softening.epsilon_comoving;
      if (softening_override != 0U) {
        effective_softening = softening_value;
      } else if (m_tree_pm_species_softening.enabled) {
        if (!core::isValidParticleSpeciesTag(species_tag)) {
          throw std::invalid_argument(
              "gravity source species tag is outside the canonical range");
        }
        effective_softening =
            m_tree_pm_species_softening.epsilon_comoving_by_species[species_tag];
      }
      if (!std::isfinite(effective_softening) || effective_softening < 0.0) {
        throw std::invalid_argument(
            "gravity source effective softening must be finite and non-negative");
      }
      if (effective_softening >
          max_certified_softening_to_split_ratio * split_scale_comoving) {
        throw std::invalid_argument(
            "TreePM softening exceeds the certified epsilon/r_s <= 0.20 envelope; "
            "reduce gravitational softening or increase the TreePM split scale");
      }
    };

    const auto appendSource = [&](double source_x,
                                  double source_y,
                                  double source_z,
                                  double source_mass,
                                  std::uint32_t species_tag,
                                  double softening_value,
                                  std::uint8_t softening_override,
                                  std::uint32_t particle_row,
                                  std::uint32_t cell_row) {
      if (!std::isfinite(source_x) || !std::isfinite(source_y) ||
          !std::isfinite(source_z) || !std::isfinite(source_mass) ||
          source_mass < 0.0) {
        throw std::runtime_error(
            "gravity source rebuild encountered non-finite coordinates or a negative/non-finite mass");
      }
      validateCertifiedTreePmSoftening(
          species_tag, softening_value, softening_override);
      if (m_local_source_x.size() >=
          static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
        throw std::overflow_error(
            "gravity local source count exceeds TreePM 32-bit source-index contract");
      }
      if (periodic_sources) {
        source_x = wrapPeriodicPosition(source_x, box_size_x);
        source_y = wrapPeriodicPosition(source_y, box_size_y);
        source_z = wrapPeriodicPosition(source_z, box_size_z);
      }
      const auto source_index =
          static_cast<std::uint32_t>(m_local_source_x.size());
      m_local_source_x.push_back(source_x);
      m_local_source_y.push_back(source_y);
      m_local_source_z.push_back(source_z);
      m_local_source_mass.push_back(source_mass);
      m_local_source_species_tag.push_back(species_tag);
      m_local_source_particle_row.push_back(particle_row);
      m_local_source_cell_row.push_back(cell_row);
      if (has_softening_values) {
        m_local_source_softening_comoving.push_back(softening_value);
      }
      if (has_softening_masks) {
        m_local_source_softening_override_mask.push_back(softening_override);
      }
      return source_index;
    };

    // Exact memoization of the cosmological drift integral (2026-10-26). The
    // factor depends only on (source scale, evaluation scale); inactive
    // sources share a handful of last-drift epochs (block time bins, common
    // DMO epoch), so the 64-point quadrature runs once per distinct epoch
    // instead of once per source. Identical inputs reproduce bitwise
    // identical factors, so predicted coordinates are unchanged.
    struct DriftFactorMemoEntry {
      std::uint64_t source_scale_bits = 0U;
      std::uint64_t evaluation_scale_bits = 0U;
      double factor = 0.0;
    };
    std::array<DriftFactorMemoEntry, 16U> drift_factor_memo{};
    std::size_t drift_factor_memo_size = 0U;
    std::size_t drift_factor_memo_next = 0U;
    m_source_prediction_drift_integrals = 0U;
    const auto memoized_drift_factor = [&](double source_scale, double evaluation_scale) {
      const auto source_bits = std::bit_cast<std::uint64_t>(source_scale);
      const auto evaluation_bits = std::bit_cast<std::uint64_t>(evaluation_scale);
      for (std::size_t i = 0; i < drift_factor_memo_size; ++i) {
        const DriftFactorMemoEntry& entry = drift_factor_memo[i];
        if (entry.source_scale_bits == source_bits && entry.evaluation_scale_bits == evaluation_bits) {
          return entry.factor;
        }
      }
      const double factor = core::computeComovingDriftFactor(
          *context.cosmology_background, source_scale, evaluation_scale, 64) /
          context.timeline_step.time_si_per_code;
      drift_factor_memo[drift_factor_memo_next] = DriftFactorMemoEntry{source_bits, evaluation_bits, factor};
      drift_factor_memo_next = (drift_factor_memo_next + 1U) % drift_factor_memo.size();
      drift_factor_memo_size = std::min(drift_factor_memo_size + 1U, drift_factor_memo.size());
      ++m_source_prediction_drift_integrals;
      return factor;
    };

    // Collisionless and compact-object particle state remains particle-owned.
    // Gas-tagged generic particles are compatibility/lineage mirrors only and
    // must not duplicate the authoritative gas-cell mass in the gravity source
    // population. Remote ghost particle copies are also excluded: TreePM owns
    // its explicit remote-source exchange from each rank's authoritative
    // owned sources.
    for (const std::uint32_t global_index :
         authoritative_source_rows.particle_rows) {
      double source_x = state.particles.position_x_comoving[global_index];
      double source_y = state.particles.position_y_comoving[global_index];
      double source_z = state.particles.position_z_comoving[global_index];
      if (prediction_epoch != SourcePredictionEpoch::kNone &&
          active_mask[global_index] == 0U) {
        const double source_time =
            state.particleLastDriftTimeCode(global_index);
        const double source_scale =
            state.particleLastDriftScaleFactor(global_index);
        const double evaluation_time =
            prediction_epoch == SourcePredictionEpoch::kStepBegin
            ? context.timeline_step.time_begin_code
            : context.timeline_step.time_end_code;
        const double evaluation_scale =
            prediction_epoch == SourcePredictionEpoch::kStepBegin
            ? context.timeline_step.scale_factor_begin
            : context.timeline_step.scale_factor_end;
        if (source_time > evaluation_time + 1.0e-12 || source_scale <= 0.0) {
          std::ostringstream message;
          message << "inactive PM source has an invalid or future drift epoch"
                  << " source_index=" << global_index
                  << " source_time=" << source_time
                  << " evaluation_time=" << evaluation_time
                  << " source_scale=" << source_scale
                  << " evaluation_scale=" << evaluation_scale;
          throw std::runtime_error(message.str());
        }
        // Prediction contract: this predictor applies the global-KDK drift
        // operator x += u * integral(dt/a) (computeComovingDriftFactor, the
        // same factor as CosmologicalTimeline::drift_factor_code) to the
        // stored peculiar-velocity lane. Hierarchical KDK instead stores
        // u = p/a at the common position epoch with an explicit Hubble drag
        // and drift a_begin * integral(dt/a^2), and it drifts every source
        // canonically each fine tick; a hierarchical source must therefore
        // never need a non-zero predicted drift here. Fail closed if it would.
        if (context.hierarchical_kdk.enabled &&
            (source_scale != evaluation_scale || source_time != evaluation_time)) {
          throw std::logic_error(
              "hierarchical KDK sources must be canonically drifted to the force epoch; "
              "the global-KDK momentum predictor does not apply to u = p/a lanes");
        }
        double inactive_drift_factor = evaluation_time - source_time;
        if (context.cosmology_background != nullptr) {
          inactive_drift_factor = memoized_drift_factor(source_scale, evaluation_scale);
        }
        if (!std::isfinite(inactive_drift_factor) ||
            inactive_drift_factor < -1.0e-14) {
          throw std::runtime_error(
              "inactive PM source prediction produced an invalid drift factor");
        }
        source_x += state.particles.velocity_x_peculiar[global_index] *
            inactive_drift_factor;
        source_y += state.particles.velocity_y_peculiar[global_index] *
            inactive_drift_factor;
        source_z += state.particles.velocity_z_peculiar[global_index] *
            inactive_drift_factor;
        ++m_source_predicted_inactive_count;
      }
      const double softening_value = has_softening_values
          ? state.particle_sidecar.gravity_softening_comoving[global_index]
          : 0.0;
      const std::uint8_t softening_override = has_softening_masks &&
              state.particle_sidecar.hasGravitySofteningOverride(global_index)
          ? 1U
          : 0U;
      const std::uint32_t source_index = appendSource(
          source_x,
          source_y,
          source_z,
          state.particles.mass_code[global_index],
          state.particleSpeciesTag(global_index),
          softening_value,
          softening_override,
          static_cast<std::uint32_t>(global_index),
          no_row);
      m_owned_local_index_by_global[global_index] =
          static_cast<int>(source_index);
    }

    // Gas-cell state is authoritative. The shared selector excludes covered
    // coarse cells and never relies on optional parent-particle lineage.
    for (const std::uint32_t cell_index :
         authoritative_source_rows.gas_cell_rows) {
      m_owned_leaf_cell_mask[cell_index] = 1U;
      const auto gas_source = internal::authoritativeGasGravitySource(
          state, cell_index, "gravity authoritative gas source rebuild");
      const std::uint32_t source_index = appendSource(
          gas_source.x_comoving,
          gas_source.y_comoving,
          gas_source.z_comoving,
          gas_source.mass_code,
          gas_species_tag,
          gas_softening_comoving,
          0U,
          no_row,
          cell_index);
      m_owned_local_index_by_cell[cell_index] = static_cast<int>(source_index);
    }

    m_local_active_indices.clear();
    m_local_active_global_indices.clear();
    m_active_target_particle_row.clear();
    m_active_target_cell_row.clear();
    const std::size_t target_capacity = core::checkedSizeAdd(
        active_particles.size(), authoritative_source_rows.gas_cell_rows.size(),
        "gravity authoritative target capacity");
    m_local_active_indices.reserve(target_capacity);
    m_local_active_global_indices.reserve(active_particles.size());
    m_active_target_particle_row.reserve(target_capacity);
    m_active_target_cell_row.reserve(target_capacity);

    const auto appendTarget = [&](std::uint32_t source_index,
                                  std::uint32_t particle_row,
                                  std::uint32_t cell_row) {
      m_local_active_indices.push_back(source_index);
      m_active_target_particle_row.push_back(particle_row);
      m_active_target_cell_row.push_back(cell_row);
    };

    for (const std::uint32_t global_index : active_particles) {
      if (global_index >= particle_count) {
        throw std::out_of_range(
            "gravity callback active particle index out of range");
      }
      if (state.particleSpeciesTag(global_index) == gas_species_tag) {
        continue;
      }
      const int local_index = m_owned_local_index_by_global[global_index];
      if (local_index < 0) {
        continue;
      }
      const auto source_index = static_cast<std::uint32_t>(local_index);
      appendTarget(source_index, global_index, no_row);
      m_local_active_global_indices.push_back(global_index);
    }

    // Rung-zero production currently owns globally synchronized force
    // refreshes. Compute every owned leaf gas-cell target at those boundaries
    // so hydro and restart caches never depend on an optional lineage parent.
    // Future local gas time bins can narrow this list only after carrying a
    // per-cell force epoch through the scheduler contract.
    for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      const int local_index = m_owned_local_index_by_cell[cell_index];
      if (local_index < 0) {
        continue;
      }
      appendTarget(
          static_cast<std::uint32_t>(local_index),
          no_row,
          static_cast<std::uint32_t>(cell_index));
    }
    m_active_source_indices = m_local_active_indices;
  }

  const RuntimeServices& m_services;
  const core::SimulationConfig& m_config;
  const core::ModePolicy& m_mode_policy;
  parallel::DistributedExecutionTopology m_runtime_topology{};
  std::uint64_t m_pm_update_cadence_steps = 1;
  gravity::PmGridShape m_pm_grid_shape{};
  double m_mesh_spacing_x_mpc_comoving = 0.0;
  double m_mesh_spacing_y_mpc_comoving = 0.0;
  double m_mesh_spacing_z_mpc_comoving = 0.0;
  double m_gravitational_constant_code = 0.0;
  std::string m_pm_assignment_scheme = "unknown";
  std::string m_pm_backend = "unknown";
  mutable gravity::TreePmCoordinator m_tree_pm_coordinator;
  gravity::TreePmOptions m_tree_pm_options;
  gravity::TreeSofteningSpeciesPolicy m_tree_pm_species_softening{};
  gravity::GravitySourceRepresentation m_source_representation =
      gravity::GravitySourceRepresentation::kMaterializedGeneric;
  std::uint64_t m_borrowed_particle_index_generation = 0U;
  std::uint64_t m_borrowed_gravity_source_generation = 0U;
  mutable bool m_borrowed_eligibility_certificate_valid = false;
  mutable bool m_borrowed_eligibility_result = false;
  mutable std::uint64_t m_borrowed_eligibility_particle_index_generation = 0U;
  mutable std::uint64_t m_borrowed_eligibility_gravity_source_generation = 0U;
  std::span<const std::uint32_t> m_active_source_indices{};
  std::vector<double> m_active_accel_x;
  std::vector<double> m_active_accel_y;
  std::vector<double> m_active_accel_z;
  std::vector<double> m_active_previous_acceleration_magnitude;
  std::vector<double> m_particle_accel_x;
  std::vector<double> m_particle_accel_y;
  std::vector<double> m_particle_accel_z;
  std::vector<double> m_cell_accel_x;
  std::vector<double> m_cell_accel_y;
  std::vector<double> m_cell_accel_z;
  std::vector<int> m_active_slot_by_particle;
  std::vector<int> m_active_slot_by_cell;
  std::vector<int> m_owned_local_index_by_global;
  std::vector<int> m_owned_local_index_by_cell;
  std::vector<std::uint8_t> m_owned_leaf_cell_mask;
  std::vector<double> m_local_source_x;
  std::vector<double> m_local_source_y;
  std::vector<double> m_local_source_z;
  std::vector<double> m_local_source_mass;
  std::vector<std::uint32_t> m_local_source_species_tag;
  std::vector<double> m_local_source_softening_comoving;
  std::vector<std::uint8_t> m_local_source_softening_override_mask;
  std::vector<std::uint32_t> m_local_source_particle_row;
  std::vector<std::uint32_t> m_local_source_cell_row;
  std::vector<std::uint32_t> m_active_target_particle_row;
  std::vector<std::uint32_t> m_active_target_cell_row;
  std::vector<std::uint8_t> m_source_is_high_res;

  std::vector<std::uint8_t> m_active_is_high_res;
  std::filesystem::path m_zoom_region_path;
  bool m_zoom_membership_loaded = false;
  std::unordered_set<std::uint64_t> m_zoom_high_res_particle_ids;
  std::string m_zoom_membership_source = "disabled";
  std::vector<std::uint32_t> m_local_active_indices;
  std::vector<std::uint32_t> m_local_active_global_indices;
  std::vector<std::uint32_t> m_force_refresh_particle_indices;
  // Stable decomposition-local leaf identities/SFC intervals replaced only at
  // segment/restart/ownership-change install; survives ordinary drift refits
  // (may include currently-empty leaves).
  std::vector<parallel::TopDomainLeaf> m_authoritative_top_domain_seed_leaves;
  // Current published routing leaves for the stamped source generation; may
  // omit currently-empty seed groups; refitted from the seed set.
  std::vector<parallel::TopDomainLeaf> m_authoritative_top_domain_leaves;
  std::uint64_t m_authoritative_top_domain_source_generation = 0U;
  bool m_authoritative_top_domain_source_generation_valid = false;
  bool m_domain_geometry_refreshed_this_solve = false;
  std::uint64_t m_domain_geometry_refreshed_leaf_count = 0U;
  std::uint64_t m_domain_geometry_out_of_seed_range_source_count = 0U;
  std::size_t m_local_kick_particle_count = 0U;
  gravity::TreePmDiagnostics m_last_tree_pm_diagnostics{};
  parallel::DecompositionRuntimeMeasurements m_last_decomposition_measurements{};
  bool m_has_long_range_field = false;
  std::uint64_t m_last_committed_field_version = 0;
  std::uint64_t m_long_range_refresh_count = 0;
  std::uint64_t m_long_range_reuse_count = 0;
  bool m_force_cache_valid = false;
  std::uint64_t m_force_cache_particle_index_generation = 0U;
  std::uint64_t m_force_cache_cell_index_generation = 0U;
  std::uint64_t m_force_cache_gas_identity_generation = 0U;
  std::uint64_t m_force_cache_source_generation = 0U;
  std::vector<std::uint8_t> m_particle_force_cache_valid;
  // Optional hierarchy only: current per-target Tree generation and coarse PM
  // history. No cold particle metadata, duplicate source view, or second tree.
  std::uint64_t m_hierarchical_pm_source_generation = 0U;
  std::uint64_t m_hierarchical_pm_field_version = 0U;
  std::vector<std::uint64_t> m_hierarchical_force_generation;
  std::vector<double> m_hierarchical_pm_x;
  std::vector<double> m_hierarchical_pm_y;
  std::vector<double> m_hierarchical_pm_z;
  std::vector<double> m_hierarchical_tree_x;
  std::vector<double> m_hierarchical_tree_y;
  std::vector<double> m_hierarchical_tree_z;
  // Nearest distinct direct-pair separation squared from the latest Tree
  // solve of each row (+infinity: none inside r_cut). Same freshness as the
  // Tree lanes (m_hierarchical_force_generation).
  std::vector<double> m_hierarchical_nearest_r2;
  std::vector<std::uint8_t> m_cell_force_cache_valid;
  gravity::DecompositionEpoch m_decomposition_epoch{};
  gravity::GravityMemoryEstimate m_last_pre_run_memory_estimate{};
  bool m_has_pre_run_memory_estimate = false;
  gravity::DmoProcessMemoryEstimate m_last_process_memory_estimate{};
  bool m_has_process_memory_estimate = false;
  parallel::GhostCacheLifecycle m_ghost_cache_lifecycle{};
  std::uint64_t m_source_predicted_inactive_count = 0;
  // Distinct cosmological drift integrals evaluated by the last prediction
  // (memoized per (source epoch, evaluation epoch); telemetry only).
  std::uint64_t m_source_prediction_drift_integrals = 0;
  std::vector<ReferenceWorkflowReport::TreePmCadenceRecord> m_cadence_records;
  // Bounded recent cadence window (accepted trajectory) and the multirate
  // cycle transaction journal (same bound; older journal entries would be
  // evicted from the window by newer ones on commit anyway).
  static constexpr std::size_t k_recent_cadence_limit = 256U;
  bool m_cycle_transaction_open = false;
  std::uint64_t m_pending_long_range_refreshes = 0U;
  std::uint64_t m_pending_long_range_reuses = 0U;
  std::vector<ReferenceWorkflowReport::TreePmCadenceRecord> m_pending_cadence_records;
  std::uint64_t m_pending_cadence_records_dropped = 0U;
  // Execution-cost totals of rejected multirate attempts.
  std::uint64_t m_rejected_long_range_refresh_count = 0U;
  std::uint64_t m_rejected_long_range_reuse_count = 0U;
  std::uint64_t m_rejected_cadence_record_count = 0U;
};

}  // namespace

void initializeGravityState(
    core::SimulationState& state,
    const core::SimulationConfig& config) {
  maybeInitializeParticleSofteningFromSpeciesPolicy(state, config);
}

std::string describeGravitySofteningPolicy(
    const core::SimulationConfig& config,
    const core::SimulationState& state) {
  return describeSofteningPolicy(config, &state);
}

std::unique_ptr<GravityRuntime> makeGravityRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const RuntimeServices& services,
    std::filesystem::path zoom_region_path) {
  return std::make_unique<GravityRuntimeImpl>(
      config, mode_policy, services, std::move(zoom_region_path));
}

}  // namespace cosmosim::workflows
