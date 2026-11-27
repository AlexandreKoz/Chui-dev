#pragma once

#include <complex>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/gravity/gravity_communication_arena.hpp"
#include "cosmosim/gravity/pm_spectral_hook.hpp"
#include "cosmosim/core/execution_policy.hpp"
#include "cosmosim/gravity/tree_index.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::gravity {

enum class PmBackendCapability : std::uint8_t {
  kUnavailable = 0,
  kDiagnosticNaiveDft = 1,
  kProductionFftw = 2,
  // Reserved for the future fully device-resident cuFFT path. The present
  // CUDA staging implementation intentionally does not claim this tier.
  kProductionCudaFft = 3,
};

struct PmBackendArchitecture {
  PmBackendCapability capability = PmBackendCapability::kUnavailable;
  bool host_fft = false;
  bool device_assignment = false;
  bool device_fft = false;
  bool device_green_function = false;
  bool device_gradient = false;
  bool device_interpolation = false;
  bool persistent_device_buffers = false;
  bool persistent_residency = false;
  bool distributed_device_fft = false;
};

// The solver-facing ownership contract is intentionally topology-neutral. The
// current production implementation supplies serial/x-slab and transposed
// spectral slab extents; a future 2-D pencil backend can populate the same
// descriptor without changing TreePM force ownership.
enum class PmDecompositionTopology : std::uint8_t {
  kSerial = 0,
  kXSlab = 1,
  kXSlabTransposedSpectral = 2,
  kPencil2D = 3,
};

struct PmOwnershipExtent3D {
  std::size_t x_begin = 0;
  std::size_t x_count = 0;
  std::size_t y_begin = 0;
  std::size_t y_count = 0;
  std::size_t z_begin = 0;
  std::size_t z_count = 0;
};

struct PmDecompositionDescriptor {
  PmDecompositionTopology topology = PmDecompositionTopology::kSerial;
  std::size_t process_grid_x = 1;
  std::size_t process_grid_y = 1;
  PmOwnershipExtent3D real_extent{};
  PmOwnershipExtent3D spectral_extent{};
};

enum class PmAssignmentScheme {
  kCic,
  kTsc,
};

enum class PmBoundaryCondition {
  kPeriodic,
  kIsolatedOpen,
};

// Real-space PM force from the periodic potential spectrum phi_k.
// kSpectral (reference, default): F_a = IFFT(-i k_a phi_k) per axis, three
//   inverse FFTs; keeps a pristine phi_k copy (PlanResources::potential_k)
//   because complex-to-real FFTW transforms destroy their input.
// kFiniteDifference4 (opt-in, unqualified): one inverse FFT of phi_k, then the
//   fourth-order centered difference
//     F_a(i) = -[8(phi_{i+1} - phi_{i-1}) - (phi_{i+2} - phi_{i-2})] / (12 h_a)
//   on the periodic mesh (h_a = L_a/n_a; two x halo planes exchanged with
//   the neighboring slab ranks). Never allocates potential_k. Its transfer
//   function is sin-like below the spectral i k_a, so forces near the mesh
//   scale are weaker than the reference; requires owner force validation.
//   Periodic host solves only; every distributed rank must own >= 2 x planes.
enum class PmForceGradientBackend : std::uint8_t {
  kSpectral = 0,
  kFiniteDifference4 = 1,
};

struct PmGridShape {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;

  [[nodiscard]] std::size_t cellCount() const;
  [[nodiscard]] bool isValid() const;
};

// Deterministic CHUI-owned storage required by one cached FFT/Poisson plan.
// Backend-internal FFTW/cuFFT allocations are intentionally excluded and must
// be carried as an explicit external reserve by the production memory policy.
struct PmPlanResourcesMemoryEstimate {
  // Zero since the 2026-10-25 in-place layout: the padded real view aliases
  // the first complex spectral array (real_view_aliases_complex).
  std::uint64_t real_array_bytes = 0U;
  std::uint64_t complex_spectral_array_bytes = 0U;
  // Full scalar spectral operators are not retained after M48-09. This field
  // remains for source compatibility and is always zero.
  std::uint64_t scalar_spectral_array_bytes = 0U;
  std::uint64_t spectral_axis_metadata_bytes = 0U;
  // kFiniteDifference4 on distributed slabs: four retained x halo planes
  // plus two send-staging planes (2026-10-26; no per-solve allocation).
  std::uint64_t finite_difference_halo_bytes = 0U;
  std::uint64_t total_owned_bytes = 0U;
  std::uint64_t logical_local_complex_cells = 0U;
  std::uint64_t allocated_local_complex_cells = 0U;
  bool used_backend_allocation_query = false;
  bool real_view_aliases_complex = false;
};

[[nodiscard]] PmPlanResourcesMemoryEstimate estimatePmPlanResourcesMemory(
    PmGridShape shape,
    const parallel::PmSlabLayout& layout,
    core::PmDecompositionMode decomposition_mode,
    PmForceGradientBackend force_gradient_backend = PmForceGradientBackend::kSpectral);

struct PmGlobalCell {
  std::size_t x = 0;
  std::size_t y = 0;
  std::size_t z = 0;
};

// Concrete numerical ownership view used by PM assignment/interpolation and
// backend remap code. It is a small value type: no virtual dispatch occurs in
// particle/cell loops. The current implementation adapts serial/x-slab FFTW
// ownership; a future pencil backend can supply the same operations without
// rewriting routing kernels.
class PmDecompositionView {
 public:
  PmDecompositionView(
      PmGridShape shape,
      parallel::PmSlabLayout layout,
      core::PmDecompositionMode mode);
  PmDecompositionView(
      PmGridShape shape,
      PmDecompositionDescriptor descriptor,
      int world_rank,
      int world_size);

  [[nodiscard]] const PmDecompositionDescriptor& descriptor() const noexcept;
  [[nodiscard]] PmOwnershipExtent3D realExtentForRank(int rank) const;
  [[nodiscard]] int realOwnerRank(PmGlobalCell global_cell) const;
  [[nodiscard]] int spectralOwnerRank(PmGlobalCell global_cell) const;
  [[nodiscard]] bool ownsRealCell(
      std::size_t global_x, std::size_t global_y, std::size_t global_z) const noexcept;
  [[nodiscard]] std::size_t globalToLocalRealIndex(
      std::size_t global_x, std::size_t global_y, std::size_t global_z) const;
  [[nodiscard]] bool ownsSpectralCell(
      std::size_t global_x, std::size_t global_y, std::size_t global_z) const noexcept;
  [[nodiscard]] std::size_t globalToLocalSpectralIndex(
      std::size_t global_x, std::size_t global_y, std::size_t global_z) const;

 private:
  PmGridShape m_shape{};
  parallel::PmSlabLayout m_layout{};
  core::PmDecompositionMode m_mode = core::PmDecompositionMode::kSlab;
  PmDecompositionDescriptor m_descriptor{};
  int m_world_rank = 0;
  int m_world_size = 1;
};

enum class PmDataResidencyPolicy {
  kHostOnly,
  kPreferDevice,
};

// M1A structural routing ceiling for CHUI-owned PM communication workspace.
// It bounds the two simultaneously retained wire buffers plus the solver-owned
// fixed routing metadata. Persistent PM/FFT arrays and MPI/FFTW internals are
// outside this contract.
inline constexpr std::uint64_t k_pm_routing_workspace_target_bytes =
    128ULL * 1024ULL * 1024ULL;
// Keep deterministic fixed slack below the acceptance ceiling so minor
// rank-scale metadata growth cannot consume the final bytes of the contract.
// All known CHUI-owned routing buffers/cursors are still accounted explicitly.
inline constexpr std::uint64_t k_pm_routing_workspace_headroom_bytes =
    64ULL * 1024ULL;
inline constexpr std::size_t k_pm_routing_max_wire_record_bytes = 96U;
inline constexpr std::uint64_t k_pm_routing_modeled_workspace_limit_bytes =
    k_pm_routing_workspace_target_bytes - k_pm_routing_workspace_headroom_bytes;
static_assert(k_pm_routing_workspace_headroom_bytes < k_pm_routing_workspace_target_bytes);

struct PmRoutingCapacityModel {
  std::uint64_t configured_per_peer_max_bytes = 0U;
  std::uint64_t effective_per_peer_payload_bytes = 0U;
  std::uint64_t max_send_payload_bytes = 0U;
  std::uint64_t max_receive_payload_bytes = 0U;
  std::uint64_t fixed_workspace_bytes = 0U;
  std::uint64_t max_simultaneous_workspace_bytes = 0U;
};

// Deterministic capacity model used by density and reverse-interpolation
// routing. The effective peer payload is record-aligned, never exceeds the
// configured per-peer maximum, and is reduced when necessary so send+receive
// payload plus fixed solver-owned routing metadata fits the aggregate limit.
[[nodiscard]] PmRoutingCapacityModel modelPmRoutingCapacity(
    int world_size,
    std::uint64_t configured_per_peer_max_bytes,
    std::uint64_t aggregate_workspace_limit_bytes,
    std::uint64_t fixed_workspace_bytes,
    std::size_t wire_record_bytes);

struct PmSolveOptions {
  double box_size_x_mpc_comoving = 0.0;
  double box_size_y_mpc_comoving = 0.0;
  double box_size_z_mpc_comoving = 0.0;
  // Legacy scalar compatibility lane; if axis-aware lengths are zero, this value is used for all axes.
  double box_size_mpc_comoving = 0.0;
  double scale_factor = 1.0;
  double gravitational_constant_code = 1.0;
  PmAssignmentScheme assignment_scheme = PmAssignmentScheme::kCic;
  bool enable_window_deconvolution = false;
  core::ExecutionPolicy execution_policy = core::ExecutionPolicy::kHostSerial;
  PmDataResidencyPolicy data_residency = PmDataResidencyPolicy::kHostOnly;
  core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;
  PmBoundaryCondition boundary_condition = PmBoundaryCondition::kPeriodic;
  // Optional TreePM long-range Gaussian split scale. <=0 disables filtering.
  double tree_pm_split_scale_comoving = 0.0;
  // Configured per-peer wire-payload maximum for distributed PM routing. The
  // solver may derive a smaller effective per-peer value so simultaneous
  // send+receive workspace remains inside k_pm_routing_workspace_target_bytes.
  // The configured value therefore remains a compatibility-preserving upper
  // bound, not a promise to allocate that much for every peer at once.
  std::uint64_t routing_exchange_batch_bytes = 16ULL * 1024ULL * 1024ULL;
  // Root-owned transient workspace budget for the isolated/open gathered PM path.
  // The implementation is correct inside this envelope but is not a scalable
  // distributed isolated PM solver.
  std::uint64_t isolated_open_root_workspace_limit_bytes = 256ULL * 1024ULL * 1024ULL;
  PmForceGradientBackend force_gradient_backend = PmForceGradientBackend::kSpectral;
  // CHUI-EXT-PHYS-REF-v1 N03: optional isotropic source multiplier M(|k|) of
  // the periodic Poisson operator (pm_spectral_hook.hpp). Not owned. Null keeps
  // the predecessor operator bit-for-bit. Refused by isolated/open and CUDA solves.
  PmSpectralSourceHook* spectral_source_hook = nullptr;
};

struct PmProfileEvent {
  std::uint64_t bytes_moved = 0;
  std::uint64_t routed_density_records = 0;
  std::uint64_t routed_force_requests = 0;
  std::uint64_t routed_potential_requests = 0;
  std::uint64_t routed_density_peer_count = 0;
  std::uint64_t routed_force_peer_count = 0;
  std::uint64_t routed_potential_peer_count = 0;
  std::uint64_t routed_mpi_bytes_sent = 0;
  std::uint64_t routed_mpi_bytes_received = 0;
  std::uint64_t routed_send_buffer_high_water_bytes = 0;
  std::uint64_t routed_receive_buffer_high_water_bytes = 0;
  std::uint64_t routed_combined_buffer_high_water_bytes = 0;
  std::uint64_t routed_workspace_high_water_bytes = 0;
  std::uint64_t force_halo_cache_hits = 0;
  // Counts rebuilds of the small O(nx+ny+nz) wave-number/window metadata
  // used to derive spectral operators. No full spectral scalar kernel is cached.
  std::uint64_t spectral_operator_rebuilds = 0;
  std::uint64_t isolated_open_root_workspace_estimate_bytes = 0;
  std::uint64_t isolated_open_root_workspace_limit_bytes = 0;
  std::uint64_t isolated_open_gather_bytes = 0;
  double assign_ms = 0.0;
  double fft_forward_ms = 0.0;
  double poisson_ms = 0.0;
  double gradient_ms = 0.0;
  double fft_inverse_ms = 0.0;
  double fft_transpose_ms = 0.0;
  double interpolate_ms = 0.0;
  // Inclusive wall time of the PM operation measured once around the long-
  // range refresh + force interpolation + optional zoom-correction path.
  // Subphase fields above overlap (FFT transpose sits inside forward/inverse
  // passes, etc.) and must not be summed for load-balance feedback; total_ms
  // is the single non-overlapping value for that purpose.
  double total_ms = 0.0;
  double routed_mpi_wait_ms = 0.0;
  // Nested inside assign_ms (deposition + routing + normalization inclusive).
  double density_routing_wait_ms = 0.0;
  // Disjoint per-axis stages within gradient_ms; inverse excludes real copy.
  std::array<double, 3U> gradient_factor_ms{};
  std::array<double, 3U> inverse_axis_ms{};
  std::array<double, 3U> normalization_axis_ms{};
  double halo_exchange_ms = 0.0;  // inclusive exchange + cache preparation/commit
  double transfer_h2d_ms = 0.0;
  double transfer_d2h_ms = 0.0;
  double device_kernel_ms = 0.0;
  std::uint64_t fft_transpose_bytes = 0;
};

class PmProfiler {
 public:
  void reset();
  void append(const PmProfileEvent& event);
  [[nodiscard]] const PmProfileEvent& totals() const;

 private:
  PmProfileEvent m_totals{};
};

class PmGridStorage {
 public:
  explicit PmGridStorage(PmGridShape shape);
  PmGridStorage(PmGridShape shape, parallel::PmSlabLayout layout);

  [[nodiscard]] const PmGridShape& shape() const;
  [[nodiscard]] const parallel::PmSlabLayout& slabLayout() const;
  [[nodiscard]] bool ownsFullDomain() const noexcept;
  [[nodiscard]] std::size_t localCellCount() const;
  [[nodiscard]] std::uint64_t storageGeneration() const noexcept;

  // Non-const density access is the explicit compatibility materialization
  // boundary for compact physical-cell storage. Const access never allocates
  // and returns an empty span until compact density has been materialized.
  [[nodiscard]] std::span<double> density();
  [[nodiscard]] std::span<const double> density() const;
  [[nodiscard]] bool hasDensityStorage() const noexcept;
  void ensureDensityStorage();

  // Non-const access is the explicit compatibility materialization boundary for
  // real-space potential storage. Const access never allocates and returns an
  // empty span until potential has been materialized by a potential-producing
  // solve or ensurePotentialStorage().
  [[nodiscard]] std::span<double> potential();
  [[nodiscard]] std::span<const double> potential() const;
  [[nodiscard]] bool hasPotentialStorage() const noexcept;
  void ensurePotentialStorage();

  [[nodiscard]] std::span<double> force_x();
  [[nodiscard]] std::span<const double> force_x() const;

  [[nodiscard]] std::span<double> force_y();
  [[nodiscard]] std::span<const double> force_y() const;

  [[nodiscard]] std::span<double> force_z();
  [[nodiscard]] std::span<const double> force_z() const;

  struct ForceHaloWriteView {
    std::span<double> left;
    std::span<double> right;
  };

  void clearForceHaloCache();
  // Transactional distributed refresh: allocate/retain the six final cache
  // lanes up front, keep valid=false while X/Y/Z are filled, and publish only
  // after every component exchange succeeds and reports identical metadata.
  void beginForceHaloCacheRefresh(std::size_t halo_depth_x);
  [[nodiscard]] ForceHaloWriteView forceHaloWriteView(std::size_t component);
  void commitForceHaloCacheRefresh(
      std::size_t halo_depth_x,
      int left_peer_rank,
      int right_peer_rank,
      std::uint64_t exchange_sequence);
  [[nodiscard]] bool hasForceHaloCache() const noexcept;
  [[nodiscard]] bool tryLoadForceFromHalo(
      std::size_t global_x,
      std::size_t global_y,
      std::size_t global_z,
      double& force_x_value,
      double& force_y_value,
      double& force_z_value) const;

  [[nodiscard]] std::size_t linearIndex(std::size_t ix, std::size_t iy, std::size_t iz) const;
  void clear();
  void appendMemoryReport(core::MemoryReportBuilder& builder) const;

 private:
  friend class PmSolver;

  // M48-09 specialization boundary: periodic FFT-backed production must not
  // retain an earlier compatibility density allocation.
  void releaseDensityStorage();

  struct ForceHaloCache {
    std::vector<double> left_force_x;
    std::vector<double> left_force_y;
    std::vector<double> left_force_z;
    std::vector<double> right_force_x;
    std::vector<double> right_force_y;
    std::vector<double> right_force_z;
    std::size_t halo_depth_x = 0;
    int left_peer_rank = -1;
    int right_peer_rank = -1;
    std::uint64_t exchange_sequence = 0;
    bool valid = false;
  };

  PmGridShape m_shape;
  parallel::PmSlabLayout m_layout;
  std::vector<double> m_density;
  std::vector<double> m_potential;
  std::vector<double> m_force_x;
  std::vector<double> m_force_y;
  std::vector<double> m_force_z;
  ForceHaloCache m_force_halo_cache;
  std::uint64_t m_storage_generation = 1U;
};

class PmSolver {
 public:
  struct PmMassSourceView {
    std::span<const double> pos_x_comoving;
    std::span<const double> pos_y_comoving;
    std::span<const double> pos_z_comoving;
    std::span<const double> mass_code;
  };

  enum class PmForceCoordinateLayout : std::uint8_t {
    // Coordinate spans are compact and have exactly one row per active target.
    kCompactActive = 0,
    // Coordinate spans are source storage; coordinate_source_index explicitly
    // maps each active target to its source row. This mode remains explicit
    // even when the active set is empty.
    kIndexedSource = 1,
  };

  enum class PmForceOutputLayout : std::uint8_t {
    // Output spans are compact and have the same length/order as active_particle_index.
    kCompactActive = 0,
    // Output spans are global particle lanes; active_particle_index scatters compact PM
    // interpolation results into the selected rows only.
    kIndexedGlobal = 1,
  };

  struct PmForceTargetView {
    std::span<const std::uint32_t> active_particle_index;
    std::span<const double> pos_x_comoving;
    std::span<const double> pos_y_comoving;
    std::span<const double> pos_z_comoving;
    std::span<double> accel_x_comoving;
    std::span<double> accel_y_comoving;
    std::span<double> accel_z_comoving;
    PmForceCoordinateLayout coordinate_layout = PmForceCoordinateLayout::kCompactActive;
    PmForceOutputLayout output_layout = PmForceOutputLayout::kCompactActive;
    // Used only with kIndexedSource. Its extent equals the active-target count,
    // while the coordinate spans retain the independent source-storage extent.
    // An empty index span therefore unambiguously means zero active targets, not
    // a change of coordinate representation.
    std::span<const TreeLocalIndex> coordinate_source_index;
  };

  explicit PmSolver(PmGridShape shape);
  ~PmSolver();
  PmSolver(PmSolver&&) noexcept;
  PmSolver& operator=(PmSolver&&) noexcept;
  PmSolver(const PmSolver&) = delete;
  PmSolver& operator=(const PmSolver&) = delete;

  [[nodiscard]] const PmGridShape& shape() const;
  void appendMemoryReport(core::MemoryReportBuilder& builder) const;
  // Destroy cached FFT/FFTW backend plans while the owning MPI session is
  // guaranteed active. Idempotent; normal distributed runtime calls this
  // before MPI_Finalize.
  void shutdownBackendResources();

  void assignDensity(
      PmGridStorage& grid,
      const PmMassSourceView& source_view,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr) const;

  void assignDensity(
      PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<const double> mass,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr) const;

  // Periodic scale-free comoving force-kernel contract on the mesh. Density is
  // deposited per comoving cell volume and Ψ is the kernel potential whose
  // gradient is multiplied by the cosmological KDK factor outside this solver:
  //   ∇²Ψ(x) = 4π G [ρ_com(x) - ρ̄_com]
  // Fourier convention (for k != 0):
  //   Ψ_k = -4π G δρ_com,k / k²
  //   A_i(k) = -i k_i Ψ_k
  // Periodic zero mode policy:
  //   Ψ_{k=0} = 0 and therefore A_{k=0} = 0.
  //
  // After return, grid.potential() and grid.force_{x,y,z}() are populated and
  // available for direct mesh inspection and interpolation.
  void solvePoissonPeriodic(PmGridStorage& grid, const PmSolveOptions& options, PmProfileEvent* profile = nullptr);
  // Production TreePM force-only contract. It computes the same periodic force
  // fields without materializing PmGridStorage's real-space potential lane.
  // Existing potential-producing PM APIs retain solvePoissonPeriodic().
  void solvePoissonPeriodicForcesOnly(
      PmGridStorage& grid,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr);
  void solvePoissonIsolatedOpen(PmGridStorage& grid, const PmSolveOptions& options, PmProfileEvent* profile = nullptr);

  void interpolateForces(
      const PmGridStorage& grid,
      const PmForceTargetView& target_view,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr) const;

  // Workload-aware distributed routing (2026-10-25; admission contract
  // repaired 2026-10-26). Collective on MPI_COMM_WORLD for distributed slabs:
  // every rank performs its local validation (span extents, index ranges,
  // box lengths, options) inside one transaction and then exactly one
  // MPI_Allreduce(MAX); a local failure on any rank makes every rank throw
  // (the failing rank rethrows its own exception, peers report the lowest
  // failing rank). Serial meshes validate locally and return zeros.
  // The returned plan sizes the communication arena before leasing. The
  // cached agreement is withdrawn by the next routing call of the same kind,
  // which re-measures its own workload and fails collectively if the plan
  // does not bound it (changed coordinates, decomposition or options), so a
  // stale plan can never undersize a reservation. Must be called by every
  // rank or by none before the matching routing call.
  struct PmRoutingWorkloadPlan {
    std::uint64_t phase_arena_bytes = 0U;
    std::uint64_t send_reservation_bytes = 0U;
    std::uint64_t receive_reservation_bytes = 0U;
    std::uint64_t configured_ceiling_bytes = 0U;  // former (P-1) x per-peer send+receive
    std::uint64_t global_max_peer_records = 0U;
    bool distributed = false;
  };
  [[nodiscard]] PmRoutingWorkloadPlan planDensityRouting(
      const PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<const double> mass,
      const PmSolveOptions& options) const;
  [[nodiscard]] PmRoutingWorkloadPlan planForceInterpolationRouting(
      const PmGridStorage& grid,
      const PmForceTargetView& target_view,
      const PmSolveOptions& options) const;

  void interpolateForces(
      const PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<double> accel_x,
      std::span<double> accel_y,
      std::span<double> accel_z,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr) const;
  // For slab-distributed meshes, interpolation uses reverse communication:
  // particle owners send weighted stencil requests to slab owners and receive
  // weighted contributions back, then accumulate locally in particle order.

  // Assignment-scheme transpose gather of mesh potential values to particles.
  // This uses the same geometric stencil/convention as interpolateForces.
  void interpolatePotential(
      const PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<double> potential,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr) const;

  void solveForParticles(
      PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<const double> mass,
      std::span<double> accel_x,
      std::span<double> accel_y,
      std::span<double> accel_z,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr);

  [[nodiscard]] static bool fftBackendAvailable();
  [[nodiscard]] static bool cudaBackendAvailable();
  [[nodiscard]] static std::string fftBackendName();
  [[nodiscard]] std::size_t cachedPlanCount() const;
  [[nodiscard]] std::size_t planBuildCount() const;

 private:
  friend class TreePmCoordinator;

  // TreePmCoordinator attaches its currently active shared communication
  // lease only for the lexical duration of one PM communication call.
  // Standalone PM calls use a lazily allocated bounded internal arena.
  void attachCommunicationArena(GravityCommunicationArena* arena) const;
  void detachCommunicationArena() const noexcept;

  void assignDensityPeriodicFftBacked(
      PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<const double> mass,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr) const;
  void assignDensityImpl(
      PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<const double> mass,
      const PmSolveOptions& options,
      bool periodic_fft_backed,
      PmProfileEvent* profile) const;
  void solvePoissonPeriodicImpl(
      PmGridStorage& grid,
      const PmSolveOptions& options,
      bool materialize_potential,
      PmProfileEvent* profile);
  void interpolateForcesImpl(
      const PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<const TreeLocalIndex> coordinate_source_index,
      PmForceCoordinateLayout coordinate_layout,
      std::span<double> accel_x,
      std::span<double> accel_y,
      std::span<double> accel_z,
      const PmSolveOptions& options,
      PmProfileEvent* profile) const;

  class Impl;
  PmGridShape m_shape;
  std::unique_ptr<Impl> m_impl;
};

[[nodiscard]] bool treePmSupportedByBuild();
[[nodiscard]] PmBackendCapability pmBackendCapability() noexcept;
[[nodiscard]] std::string_view pmBackendCapabilityName(PmBackendCapability capability) noexcept;
[[nodiscard]] bool pmBackendProductionReady() noexcept;
[[nodiscard]] PmBackendArchitecture pmBackendArchitecture() noexcept;
[[nodiscard]] std::string_view pmDecompositionTopologyName(PmDecompositionTopology topology) noexcept;
[[nodiscard]] PmDecompositionDescriptor describePmDecomposition(
    PmGridShape shape,
    const parallel::PmSlabLayout& layout,
    core::PmDecompositionMode mode);
void requireTreePmSupportOrThrow(
    core::GravitySolver gravity_solver,
    bool allow_diagnostic_naive_dft = false);

}  // namespace cosmosim::gravity
