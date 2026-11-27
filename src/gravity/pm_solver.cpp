#include "cosmosim/gravity/pm_solver.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <limits>
#include <memory_resource>
#include <new>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/device_buffer.hpp"
#include "cosmosim/gravity/tree_pm_split_kernel.hpp"
#if COSMOSIM_HAVE_OPENMP
#include <omp.h>
#endif
#if COSMOSIM_ENABLE_CUDA
#include <cuda_runtime.h>
#include "cosmosim/gravity/pm_cuda_kernels.hpp"
#endif

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif
#if COSMOSIM_ENABLE_FFTW
#include <fftw3.h>
#if COSMOSIM_ENABLE_MPI
#include <fftw3-mpi.h>
#endif
#endif

namespace cosmosim::gravity {
namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

#if COSMOSIM_ENABLE_FFTW
// Threaded FFTW (libfftw3_omp) when linked: each plan uses the OpenMP team
// size current at plan creation. fftw_init_threads must precede every other
// FFTW call, including fftw_mpi_init, and is process-global and idempotent.
// Transforms stay deterministic for a fixed plan rigor and thread count; a
// different thread count may change results at round-off level.
void configureFftwPlannerThreads() {
#if COSMOSIM_FFTW_THREADS && COSMOSIM_HAVE_OPENMP
  static const bool k_threads_ready = fftw_init_threads() != 0;
  if (k_threads_ready) fftw_plan_with_nthreads(std::max(1, omp_get_max_threads()));
#endif
}
#endif

#if COSMOSIM_ENABLE_MPI
[[nodiscard]] std::string pmMpiErrorText(int error_code) {
  char message[MPI_MAX_ERROR_STRING]{};
  int length = 0;
  if (MPI_Error_string(error_code, message, &length) != MPI_SUCCESS || length <= 0) {
    return "MPI error code " + std::to_string(error_code);
  }
  return std::string(message, static_cast<std::size_t>(length));
}

void requirePmMpiSuccess(int error_code, std::string_view context) {
  if (error_code != MPI_SUCCESS) {
    throw std::runtime_error(
        std::string(context) + " failed: " + pmMpiErrorText(error_code));
  }
}

bool queryActiveMpiWorld(int& world_size, int& world_rank) noexcept {
  world_size = 1;
  world_rank = 0;
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized == 0) {
    return false;
  }
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized != 0) {
    return false;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  return true;
}
#endif

[[nodiscard]] std::size_t wrapIndex(std::ptrdiff_t i, std::size_t n) {
  const std::ptrdiff_t n_signed = static_cast<std::ptrdiff_t>(n);
  std::ptrdiff_t wrapped = i % n_signed;
  if (wrapped < 0) {
    wrapped += n_signed;
  }
  return static_cast<std::size_t>(wrapped);
}

[[nodiscard]] double wrapPosition(double x, double box_size) {
  const double wrapped = std::fmod(x, box_size);
  return wrapped < 0.0 ? wrapped + box_size : wrapped;
}

[[nodiscard]] bool positionInsideOpenDomain(double x, double box_size) {
  return x >= 0.0 && x < box_size;
}

[[nodiscard]] double sinc(double x) {
  if (std::abs(x) < 1.0e-12) {
    return 1.0;
  }
  return std::sin(x) / x;
}

struct PmAxisStencil1d {
  std::array<std::ptrdiff_t, 3> offsets{};
  std::array<double, 3> weights{};
  std::size_t count = 0;
};

struct PmDensityContributionRecord {
  std::uint32_t origin_rank = 0;
  std::uint32_t destination_rank = 0;
  std::uint32_t record_sequence = 0;
  std::uint32_t global_ix = 0;
  std::uint32_t global_iy = 0;
  std::uint32_t global_iz = 0;
  std::uint64_t exchange_epoch = 0;
  double mass_contribution = 0.0;
};

struct PmForceContributionRecord {
  std::uint32_t source_rank = 0;
  std::uint32_t origin_rank = 0;
  std::uint32_t request_sequence = 0;
  std::uint32_t origin_particle_index = 0;
  std::uint64_t exchange_epoch = 0;
  double accel_x = 0.0;
  double accel_y = 0.0;
  double accel_z = 0.0;
};

struct PmPotentialContributionRecord {
  std::uint32_t source_rank = 0;
  std::uint32_t origin_rank = 0;
  std::uint32_t request_sequence = 0;
  std::uint32_t origin_particle_index = 0;
  std::uint64_t exchange_epoch = 0;
  double potential = 0.0;
};

struct PmDensityPlaneRecord {
  std::uint32_t origin_rank = 0;
  std::uint32_t destination_rank = 0;
  std::uint32_t record_sequence = 0;
  std::uint32_t plane_count = 0;
  std::array<std::uint32_t, 3> global_ix{};
  std::uint64_t exchange_epoch = 0;
  double y_grid = 0.0;
  double z_grid = 0.0;
  double mass_code = 0.0;
  std::array<double, 3> x_weight{};
};

struct PmPlaneInterpolationRequestRecord {
  std::uint32_t origin_rank = 0;
  std::uint32_t destination_rank = 0;
  std::uint32_t request_sequence = 0;
  std::uint32_t origin_particle_index = 0;
  std::uint32_t plane_count = 0;
  std::array<std::uint32_t, 3> global_ix{};
  std::uint64_t exchange_epoch = 0;
  double y_grid = 0.0;
  double z_grid = 0.0;
  std::array<double, 3> x_weight{};
};

[[nodiscard, maybe_unused]] std::string pmRoutingDiagnostic(
    std::string_view stage,
    int receiver_rank,
    int sender_rank,
    std::uint64_t exchange_epoch,
    std::uint32_t request_sequence,
    std::uint32_t origin_particle_index,
    std::string_view detail) {
  return std::string(stage) + " routing failure local_receiver_rank=" + std::to_string(receiver_rank) +
      " sender_rank=" + std::to_string(sender_rank) + " exchange_epoch=" + std::to_string(exchange_epoch) +
      " request_sequence=" + std::to_string(request_sequence) +
      " origin_particle_index=" + std::to_string(origin_particle_index) + ": " + std::string(detail);
}

constexpr std::uint32_t k_pm_wire_magic = 0x31574d50U;  // "PMW1" in little-endian byte order.
constexpr std::uint32_t k_pm_wire_version = 1U;

enum class PmWireRecordKind : std::uint32_t {
  kDensityContribution = 1U,
  kForceResponse = 3U,
  kPotentialResponse = 5U,
  kDensityPlane = 6U,
  kForcePlaneRequest = 7U,
  kPotentialPlaneRequest = 8U,
};

constexpr std::size_t k_pm_density_wire_bytes = 56U;
constexpr std::size_t k_pm_density_plane_wire_bytes = 96U;
constexpr std::size_t k_pm_plane_interpolation_request_wire_bytes = 96U;
constexpr std::size_t k_pm_force_response_wire_bytes = 64U;
constexpr std::size_t k_pm_potential_response_wire_bytes = 48U;
static_assert(sizeof(double) == sizeof(std::uint64_t), "PM wire protocol requires binary64-sized doubles");
static_assert(
    std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53,
    "PM wire protocol requires IEEE-754 binary64 doubles");

void writePmWireU32(std::span<std::uint8_t> bytes, std::size_t& offset, std::uint32_t value) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(value)) {
    throw std::overflow_error("PM wire encoder exceeded its fixed record size");
  }
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    bytes[offset++] = static_cast<std::uint8_t>((value >> shift) & 0xffU);
  }
}

void writePmWireU64(std::span<std::uint8_t> bytes, std::size_t& offset, std::uint64_t value) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(value)) {
    throw std::overflow_error("PM wire encoder exceeded its fixed record size");
  }
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    bytes[offset++] = static_cast<std::uint8_t>((value >> shift) & 0xffU);
  }
}

void writePmWireDouble(std::span<std::uint8_t> bytes, std::size_t& offset, double value) {
  writePmWireU64(bytes, offset, std::bit_cast<std::uint64_t>(value));
}

[[nodiscard]] std::uint32_t readPmWireU32(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset,
    std::string_view context) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(std::uint32_t)) {
    throw std::runtime_error(std::string(context) + " PM wire record is truncated");
  }
  std::uint32_t value = 0U;
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] std::uint64_t readPmWireU64(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset,
    std::string_view context) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(std::uint64_t)) {
    throw std::runtime_error(std::string(context) + " PM wire record is truncated");
  }
  std::uint64_t value = 0U;
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] double readPmWireDouble(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset,
    std::string_view context) {
  return std::bit_cast<double>(readPmWireU64(bytes, offset, context));
}

void writePmWireHeader(
    std::span<std::uint8_t> bytes,
    std::size_t& offset,
    PmWireRecordKind kind) {
  writePmWireU32(bytes, offset, k_pm_wire_magic);
  writePmWireU32(bytes, offset, k_pm_wire_version);
  writePmWireU32(bytes, offset, static_cast<std::uint32_t>(kind));
}

void readAndValidatePmWireHeader(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset,
    PmWireRecordKind expected_kind,
    std::string_view context) {
  const std::uint32_t magic = readPmWireU32(bytes, offset, context);
  const std::uint32_t version = readPmWireU32(bytes, offset, context);
  const std::uint32_t kind = readPmWireU32(bytes, offset, context);
  if (magic != k_pm_wire_magic) {
    throw std::runtime_error(std::string(context) + " PM wire magic is invalid");
  }
  if (version != k_pm_wire_version) {
    throw std::runtime_error(
        std::string(context) + " PM wire version is unsupported: received=" + std::to_string(version) +
        " expected=" + std::to_string(k_pm_wire_version));
  }
  if (kind != static_cast<std::uint32_t>(expected_kind)) {
    throw std::runtime_error(
        std::string(context) + " PM wire record kind is invalid: received=" + std::to_string(kind) +
        " expected=" + std::to_string(static_cast<std::uint32_t>(expected_kind)));
  }
}

[[maybe_unused]] void encodePmDensityRecord(
    const PmDensityContributionRecord& record,
    std::span<std::uint8_t> bytes) {
  if (bytes.size() != k_pm_density_wire_bytes) {
    throw std::invalid_argument("PM density wire encoder requires an exact 56-byte record");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, PmWireRecordKind::kDensityContribution);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.destination_rank);
  writePmWireU32(bytes, offset, record.record_sequence);
  writePmWireU32(bytes, offset, record.global_ix);
  writePmWireU32(bytes, offset, record.global_iy);
  writePmWireU32(bytes, offset, record.global_iz);
  writePmWireU32(bytes, offset, 0U);
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.mass_contribution);
  if (offset != bytes.size()) {
    throw std::logic_error("PM density wire encoder size contract drifted");
  }
}

[[nodiscard, maybe_unused]] PmDensityContributionRecord decodePmDensityRecord(
    std::span<const std::uint8_t> bytes) {
  constexpr std::string_view context = "PM density contribution";
  if (bytes.size() != k_pm_density_wire_bytes) {
    throw std::runtime_error("PM density wire decoder requires an exact 56-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, PmWireRecordKind::kDensityContribution, context);
  PmDensityContributionRecord record;
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.destination_rank = readPmWireU32(bytes, offset, context);
  record.record_sequence = readPmWireU32(bytes, offset, context);
  record.global_ix = readPmWireU32(bytes, offset, context);
  record.global_iy = readPmWireU32(bytes, offset, context);
  record.global_iz = readPmWireU32(bytes, offset, context);
  const std::uint32_t reserved = readPmWireU32(bytes, offset, context);
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.mass_contribution = readPmWireDouble(bytes, offset, context);
  if (reserved != 0U || offset != bytes.size()) {
    throw std::runtime_error("PM density wire record has nonzero reserved flags or trailing data");
  }
  return record;
}

void encodePmDensityPlaneRecord(
    const PmDensityPlaneRecord& record,
    std::span<std::uint8_t> bytes) {
  if (bytes.size() != k_pm_density_plane_wire_bytes) {
    throw std::invalid_argument("PM density-plane wire encoder requires an exact 96-byte record");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, PmWireRecordKind::kDensityPlane);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.destination_rank);
  writePmWireU32(bytes, offset, record.record_sequence);
  writePmWireU32(bytes, offset, record.plane_count);
  for (const std::uint32_t ix : record.global_ix) {
    writePmWireU32(bytes, offset, ix);
  }
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.y_grid);
  writePmWireDouble(bytes, offset, record.z_grid);
  writePmWireDouble(bytes, offset, record.mass_code);
  for (const double weight : record.x_weight) {
    writePmWireDouble(bytes, offset, weight);
  }
  if (offset != bytes.size()) {
    throw std::logic_error("PM density-plane wire encoder size contract drifted");
  }
}

[[nodiscard]] PmDensityPlaneRecord decodePmDensityPlaneRecord(
    std::span<const std::uint8_t> bytes) {
  constexpr std::string_view context = "PM density-plane request";
  if (bytes.size() != k_pm_density_plane_wire_bytes) {
    throw std::runtime_error("PM density-plane wire decoder requires an exact 96-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, PmWireRecordKind::kDensityPlane, context);
  PmDensityPlaneRecord record;
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.destination_rank = readPmWireU32(bytes, offset, context);
  record.record_sequence = readPmWireU32(bytes, offset, context);
  record.plane_count = readPmWireU32(bytes, offset, context);
  for (std::uint32_t& ix : record.global_ix) {
    ix = readPmWireU32(bytes, offset, context);
  }
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.y_grid = readPmWireDouble(bytes, offset, context);
  record.z_grid = readPmWireDouble(bytes, offset, context);
  record.mass_code = readPmWireDouble(bytes, offset, context);
  for (double& weight : record.x_weight) {
    weight = readPmWireDouble(bytes, offset, context);
  }
  if (offset != bytes.size()) {
    throw std::runtime_error("PM density-plane wire record has trailing data");
  }
  return record;
}

void encodePmPlaneInterpolationRequest(
    const PmPlaneInterpolationRequestRecord& record,
    PmWireRecordKind kind,
    std::span<std::uint8_t> bytes) {
  if ((kind != PmWireRecordKind::kForcePlaneRequest &&
       kind != PmWireRecordKind::kPotentialPlaneRequest) ||
      bytes.size() != k_pm_plane_interpolation_request_wire_bytes) {
    throw std::invalid_argument(
        "PM plane-interpolation wire encoder received an invalid kind or record size");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, kind);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.destination_rank);
  writePmWireU32(bytes, offset, record.request_sequence);
  writePmWireU32(bytes, offset, record.origin_particle_index);
  writePmWireU32(bytes, offset, record.plane_count);
  for (const std::uint32_t ix : record.global_ix) {
    writePmWireU32(bytes, offset, ix);
  }
  writePmWireU32(bytes, offset, 0U);
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.y_grid);
  writePmWireDouble(bytes, offset, record.z_grid);
  for (const double weight : record.x_weight) {
    writePmWireDouble(bytes, offset, weight);
  }
  if (offset != bytes.size()) {
    throw std::logic_error("PM plane-interpolation wire encoder size contract drifted");
  }
}

[[nodiscard]] PmPlaneInterpolationRequestRecord decodePmPlaneInterpolationRequest(
    std::span<const std::uint8_t> bytes,
    PmWireRecordKind expected_kind,
    std::string_view context) {
  if (expected_kind != PmWireRecordKind::kForcePlaneRequest &&
      expected_kind != PmWireRecordKind::kPotentialPlaneRequest) {
    throw std::invalid_argument("PM plane-interpolation decoder received an invalid record kind");
  }
  if (bytes.size() != k_pm_plane_interpolation_request_wire_bytes) {
    throw std::runtime_error("PM plane-interpolation wire decoder requires an exact 96-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, expected_kind, context);
  PmPlaneInterpolationRequestRecord record;
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.destination_rank = readPmWireU32(bytes, offset, context);
  record.request_sequence = readPmWireU32(bytes, offset, context);
  record.origin_particle_index = readPmWireU32(bytes, offset, context);
  record.plane_count = readPmWireU32(bytes, offset, context);
  for (std::uint32_t& ix : record.global_ix) {
    ix = readPmWireU32(bytes, offset, context);
  }
  const std::uint32_t reserved = readPmWireU32(bytes, offset, context);
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.y_grid = readPmWireDouble(bytes, offset, context);
  record.z_grid = readPmWireDouble(bytes, offset, context);
  for (double& weight : record.x_weight) {
    weight = readPmWireDouble(bytes, offset, context);
  }
  if (reserved != 0U || offset != bytes.size()) {
    throw std::runtime_error(
        "PM plane-interpolation wire record has nonzero reserved flags or trailing data");
  }
  return record;
}

[[maybe_unused]] void encodePmForceResponse(
    const PmForceContributionRecord& record,
    std::span<std::uint8_t> bytes) {
  if (bytes.size() != k_pm_force_response_wire_bytes) {
    throw std::invalid_argument("PM force response wire encoder requires an exact 64-byte record");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, PmWireRecordKind::kForceResponse);
  writePmWireU32(bytes, offset, record.source_rank);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.request_sequence);
  writePmWireU32(bytes, offset, record.origin_particle_index);
  writePmWireU32(bytes, offset, 0U);
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.accel_x);
  writePmWireDouble(bytes, offset, record.accel_y);
  writePmWireDouble(bytes, offset, record.accel_z);
  if (offset != bytes.size()) {
    throw std::logic_error("PM force response wire encoder size contract drifted");
  }
}

[[nodiscard, maybe_unused]] PmForceContributionRecord decodePmForceResponse(
    std::span<const std::uint8_t> bytes) {
  constexpr std::string_view context = "PM force response";
  if (bytes.size() != k_pm_force_response_wire_bytes) {
    throw std::runtime_error("PM force response wire decoder requires an exact 64-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, PmWireRecordKind::kForceResponse, context);
  PmForceContributionRecord record;
  record.source_rank = readPmWireU32(bytes, offset, context);
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.request_sequence = readPmWireU32(bytes, offset, context);
  record.origin_particle_index = readPmWireU32(bytes, offset, context);
  const std::uint32_t reserved = readPmWireU32(bytes, offset, context);
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.accel_x = readPmWireDouble(bytes, offset, context);
  record.accel_y = readPmWireDouble(bytes, offset, context);
  record.accel_z = readPmWireDouble(bytes, offset, context);
  if (reserved != 0U || offset != bytes.size()) {
    throw std::runtime_error("PM force response wire record has nonzero reserved flags or trailing data");
  }
  return record;
}

[[maybe_unused]] void encodePmPotentialResponse(
    const PmPotentialContributionRecord& record,
    std::span<std::uint8_t> bytes) {
  if (bytes.size() != k_pm_potential_response_wire_bytes) {
    throw std::invalid_argument("PM potential response wire encoder requires an exact 48-byte record");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, PmWireRecordKind::kPotentialResponse);
  writePmWireU32(bytes, offset, record.source_rank);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.request_sequence);
  writePmWireU32(bytes, offset, record.origin_particle_index);
  writePmWireU32(bytes, offset, 0U);
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.potential);
  if (offset != bytes.size()) {
    throw std::logic_error("PM potential response wire encoder size contract drifted");
  }
}

[[nodiscard, maybe_unused]] PmPotentialContributionRecord decodePmPotentialResponse(
    std::span<const std::uint8_t> bytes) {
  constexpr std::string_view context = "PM potential response";
  if (bytes.size() != k_pm_potential_response_wire_bytes) {
    throw std::runtime_error("PM potential response wire decoder requires an exact 48-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, PmWireRecordKind::kPotentialResponse, context);
  PmPotentialContributionRecord record;
  record.source_rank = readPmWireU32(bytes, offset, context);
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.request_sequence = readPmWireU32(bytes, offset, context);
  record.origin_particle_index = readPmWireU32(bytes, offset, context);
  const std::uint32_t reserved = readPmWireU32(bytes, offset, context);
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.potential = readPmWireDouble(bytes, offset, context);
  if (reserved != 0U || offset != bytes.size()) {
    throw std::runtime_error("PM potential response wire record has nonzero reserved flags or trailing data");
  }
  return record;
}

[[nodiscard, maybe_unused]] std::size_t checkedPmWireByteCount(
    std::size_t record_count,
    std::size_t record_bytes,
    std::string_view context) {
  if (record_bytes == 0U || record_count > std::numeric_limits<std::size_t>::max() / record_bytes) {
    throw std::overflow_error(std::string(context) + " PM wire byte count overflows size_t");
  }
  return record_count * record_bytes;
}

template <typename Record, typename Encoder>
void encodePmWireRecords(
    std::span<const Record> records,
    std::size_t record_bytes,
    std::vector<std::uint8_t>& wire,
    Encoder&& encoder,
    std::string_view context) {
  wire.resize(checkedPmWireByteCount(records.size(), record_bytes, context));
  for (std::size_t i = 0; i < records.size(); ++i) {
    encoder(records[i], std::span<std::uint8_t>(wire).subspan(i * record_bytes, record_bytes));
  }
}

template <typename Record, typename Decoder>
void decodePmWireRecords(
    std::span<const std::uint8_t> wire,
    std::size_t record_bytes,
    std::vector<Record>& records,
    Decoder&& decoder,
    std::string_view context) {
  if (record_bytes == 0U || wire.size() % record_bytes != 0U) {
    throw std::runtime_error(std::string(context) + " PM wire payload is not record-aligned");
  }
  const std::size_t record_count = wire.size() / record_bytes;
  records.resize(record_count);
  for (std::size_t i = 0; i < record_count; ++i) {
    records[i] = decoder(wire.subspan(i * record_bytes, record_bytes));
  }
}

#if COSMOSIM_ENABLE_MPI
template <typename ByteVector>
[[nodiscard]] const std::uint8_t* nonNullPmWireData(const ByteVector& bytes) {
  static constexpr std::uint8_t empty_payload = 0U;
  return bytes.empty() ? &empty_payload : bytes.data();
}

template <typename ByteVector>
[[nodiscard]] std::uint8_t* nonNullPmWireData(ByteVector& bytes) {
  static std::uint8_t empty_payload = 0U;
  return bytes.empty() ? &empty_payload : bytes.data();
}

void validatePmExchangeEpochConsensus(std::uint64_t exchange_epoch, std::string_view context) {
  std::uint64_t minimum_epoch = 0U;
  std::uint64_t maximum_epoch = 0U;
  requirePmMpiSuccess(
      MPI_Allreduce(&exchange_epoch, &minimum_epoch, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD),
      "PM exchange-epoch minimum MPI_Allreduce");
  requirePmMpiSuccess(
      MPI_Allreduce(&exchange_epoch, &maximum_epoch, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD),
      "PM exchange-epoch maximum MPI_Allreduce");
  if (minimum_epoch != maximum_epoch) {
    throw std::runtime_error(
        std::string(context) + " ranks disagree on PM exchange epoch: minimum=" +
        std::to_string(minimum_epoch) + " maximum=" + std::to_string(maximum_epoch));
  }
}

template <typename Callable>
void measurePmMpiWait(double& accumulated_ms, Callable&& call) {
  const auto start = std::chrono::steady_clock::now();
  call();
  const auto stop = std::chrono::steady_clock::now();
  accumulated_ms += std::chrono::duration<double, std::milli>(stop - start).count();
}

enum class PmCollectiveEntryKind : int {
  kAssignDensity = 1,
  kSolvePoissonPeriodic = 2,
  kInterpolateForces = 3,
  kInterpolatePotential = 4,
  kSolvePoissonIsolatedOpen = 5,
  kInterpolateForceTargetView = 6,
  kSolveForParticles = 7,
  kSolvePoissonPeriodicForcesOnly = 8,
};

void validatePmCollectiveEntryConsensus(
    PmCollectiveEntryKind entry_kind,
    bool rank_local_serial_layout,
    const PmGridShape& shape,
    const PmSolveOptions& options,
    double& accumulated_mpi_wait_ms) {
  // Public PM entry is collective over MPI_COMM_WORLD. Bind every
  // communicator-global mesh/solver control that can select a different
  // exchange, FFT plan, or kernel before any rank enters a layout-specific
  // phase. Rank-local target coordinate/output layouts are deliberately not
  // fingerprinted: empty/uneven ranks may legitimately use indexed-source or
  // compact-active target views while participating in the same PM protocol.
  // Those representation contracts are validated in the coordinated local
  // preflight below rather than promoted to false communicator invariants.
  const std::array<std::uint64_t, 20> local_fingerprint{
      static_cast<std::uint64_t>(entry_kind),
      rank_local_serial_layout ? 1U : 0U,
      static_cast<std::uint64_t>(shape.nx),
      static_cast<std::uint64_t>(shape.ny),
      static_cast<std::uint64_t>(shape.nz),
      std::bit_cast<std::uint64_t>(options.box_size_mpc_comoving),
      std::bit_cast<std::uint64_t>(options.box_size_x_mpc_comoving),
      std::bit_cast<std::uint64_t>(options.box_size_y_mpc_comoving),
      std::bit_cast<std::uint64_t>(options.box_size_z_mpc_comoving),
      std::bit_cast<std::uint64_t>(options.scale_factor),
      std::bit_cast<std::uint64_t>(options.gravitational_constant_code),
      static_cast<std::uint64_t>(options.assignment_scheme),
      options.enable_window_deconvolution ? 1U : 0U,
      // The gradient backend shares this lane so the fingerprint width stays fixed.
      static_cast<std::uint64_t>(options.execution_policy) |
          (static_cast<std::uint64_t>(options.force_gradient_backend) << 8U),
      static_cast<std::uint64_t>(options.data_residency),
      static_cast<std::uint64_t>(options.decomposition_mode),
      static_cast<std::uint64_t>(options.boundary_condition),
      // The N03 spectral hook identity shares the split-scale lane so the
      // fingerprint width stays fixed; any rank-divergent hook state differs.
      std::bit_cast<std::uint64_t>(options.tree_pm_split_scale_comoving) ^
          (options.spectral_source_hook == nullptr ? 0U : options.spectral_source_hook->fingerprint()),
      options.routing_exchange_batch_bytes,
      options.isolated_open_root_workspace_limit_bytes,
  };
  std::array<std::uint64_t, 20> minimum_fingerprint{};
  std::array<std::uint64_t, 20> maximum_fingerprint{};
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    requirePmMpiSuccess(
        MPI_Allreduce(
            local_fingerprint.data(),
            minimum_fingerprint.data(),
            static_cast<int>(local_fingerprint.size()),
            MPI_UINT64_T,
            MPI_MIN,
            MPI_COMM_WORLD),
        "PM collective fingerprint minimum MPI_Allreduce");
  });
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    requirePmMpiSuccess(
        MPI_Allreduce(
            local_fingerprint.data(),
            maximum_fingerprint.data(),
            static_cast<int>(local_fingerprint.size()),
            MPI_UINT64_T,
            MPI_MAX,
            MPI_COMM_WORLD),
        "PM collective fingerprint maximum MPI_Allreduce");
  });
  if (minimum_fingerprint != maximum_fingerprint) {
    throw std::runtime_error(
        "PM collective entry disagreement: api_kind_min=" + std::to_string(minimum_fingerprint[0]) +
        " api_kind_max=" + std::to_string(maximum_fingerprint[0]) +
        " rank_local_serial_min=" + std::to_string(minimum_fingerprint[1]) +
        " rank_local_serial_max=" + std::to_string(maximum_fingerprint[1]) +
        " mesh_min=(" + std::to_string(minimum_fingerprint[2]) + "," +
        std::to_string(minimum_fingerprint[3]) + "," +
        std::to_string(minimum_fingerprint[4]) + ") mesh_max=(" +
        std::to_string(maximum_fingerprint[2]) + "," +
        std::to_string(maximum_fingerprint[3]) + "," +
        std::to_string(maximum_fingerprint[4]) + ")");
  }
}

void throwIfPmPayloadValidationFailed(
    std::string_view stage,
    int world_rank,
    int world_size,
    const std::string& local_error,
    double& accumulated_mpi_wait_ms) {
  const int local_failure = local_error.empty() ? 0 : 1;
  int global_failure_count = 0;
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    requirePmMpiSuccess(
        MPI_Allreduce(&local_failure, &global_failure_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD),
        "PM payload-validation failure-count MPI_Allreduce");
  });
  if (global_failure_count == 0) {
    return;
  }

  const int local_failure_rank = local_failure != 0 ? world_rank : world_size;
  int first_failure_rank = world_size;
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    requirePmMpiSuccess(
        MPI_Allreduce(&local_failure_rank, &first_failure_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD),
        "PM payload-validation first-rank MPI_Allreduce");
  });
  const std::string detail = local_failure != 0
      ? " local_detail=" + local_error
      : " local_detail=peer_rank_rejected_payload";
  throw std::runtime_error(
      std::string(stage) + " coordinated PM payload validation failed: failure_count=" +
      std::to_string(global_failure_count) + " first_failure_rank=" + std::to_string(first_failure_rank) +
      detail);
}

template <typename Callable>
void runPmCoordinatedPhase(
    std::string_view stage,
    int world_rank,
    int world_size,
    double& accumulated_mpi_wait_ms,
    Callable&& phase) {
  std::string local_error;
  try {
    phase();
  } catch (const std::exception& error) {
    local_error = error.what();
  } catch (...) {
    local_error = "non-standard exception";
  }
  throwIfPmPayloadValidationFailed(
      stage,
      world_rank,
      world_size,
      local_error,
      accumulated_mpi_wait_ms);
}

[[nodiscard]] int checkedMpiRecordLimitForByteTransport(
    std::size_t record_size,
    std::string_view context) {
  const std::uint64_t max_mpi_int = static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  if (record_size == 0U || record_size > max_mpi_int) {
    throw std::invalid_argument(
        std::string(context) + " record size=" + std::to_string(record_size) +
        " cannot be represented by MPI_BYTE int arguments");
  }
  return static_cast<int>(max_mpi_int / static_cast<std::uint64_t>(record_size));
}

[[nodiscard]] int checkedMpiByteCount(
    std::uint64_t byte_value,
    std::string_view context) {
  const std::uint64_t max_mpi_int = static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  if (byte_value > max_mpi_int) {
    throw std::invalid_argument(
        std::string(context) + " exceeds MPI int byte-count/displacement capacity: bytes=" +
        std::to_string(byte_value));
  }
  return static_cast<int>(byte_value);
}

[[nodiscard]] int checkedMpiRecordCountOrDisplacement(
    std::size_t record_value,
    std::size_t record_size,
    std::string_view context,
    std::string_view quantity,
    int rank) {
  const int record_limit = checkedMpiRecordLimitForByteTransport(record_size, context);
  if (record_value > static_cast<std::size_t>(record_limit)) {
    throw std::invalid_argument(
        std::string(context) + " " + std::string(quantity) + " exceeds MPI_BYTE record limit on rank " +
        std::to_string(rank) + ": records=" + std::to_string(record_value) +
        " record_size=" + std::to_string(record_size) +
        " max_records=" + std::to_string(record_limit));
  }
  return static_cast<int>(record_value);
}

[[nodiscard]] std::size_t checkedMpiReceivedRecordCount(
    int record_count,
    std::size_t record_size,
    std::string_view context,
    int rank) {
  if (record_count < 0) {
    throw std::invalid_argument(
        std::string(context) + " received negative record count from rank " + std::to_string(rank) +
        ": count=" + std::to_string(record_count));
  }
  const std::size_t count = static_cast<std::size_t>(record_count);
  static_cast<void>(checkedMpiRecordCountOrDisplacement(
      count, record_size, context, "received record count", rank));
  return count;
}

[[nodiscard]] std::size_t checkedMpiRecordTotal(
    std::size_t current_total,
    std::size_t record_count,
    std::size_t record_size,
    std::string_view context,
    int rank) {
  if (record_count > std::numeric_limits<std::size_t>::max() - current_total) {
    throw std::invalid_argument(
        std::string(context) + " record total overflows size_t while processing rank " + std::to_string(rank));
  }
  const std::size_t next_total = current_total + record_count;
  static_cast<void>(checkedMpiRecordCountOrDisplacement(
      next_total, record_size, context, "cumulative record total", rank));
  return next_total;
}

void checkedMpiRecordLayoutToByteLayout(
    std::span<const int> record_counts,
    std::span<const int> record_displacements,
    std::size_t record_size,
    std::span<int> byte_counts,
    std::span<int> byte_displacements,
    std::string_view context) {
  if (record_counts.size() != record_displacements.size() ||
      record_counts.size() != byte_counts.size() ||
      record_counts.size() != byte_displacements.size()) {
    throw std::invalid_argument(
        std::string(context) + " MPI record/byte count-displacement span sizes do not match");
  }

  const int record_limit = checkedMpiRecordLimitForByteTransport(record_size, context);
  const std::uint64_t max_mpi_int = static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  const std::uint64_t record_size_u64 = static_cast<std::uint64_t>(record_size);
  for (std::size_t rank = 0; rank < record_counts.size(); ++rank) {
    const int record_count = record_counts[rank];
    const int record_displacement = record_displacements[rank];
    if (record_count < 0 || record_displacement < 0) {
      throw std::invalid_argument(
          std::string(context) + " contains negative record count/displacement for rank " +
          std::to_string(rank) + ": count=" + std::to_string(record_count) +
          " displacement=" + std::to_string(record_displacement));
    }

    const std::uint64_t count_u64 = static_cast<std::uint64_t>(record_count);
    const std::uint64_t displacement_u64 = static_cast<std::uint64_t>(record_displacement);
    if (count_u64 > static_cast<std::uint64_t>(record_limit) ||
        displacement_u64 > static_cast<std::uint64_t>(record_limit)) {
      throw std::invalid_argument(
          std::string(context) + " record count/displacement exceeds MPI_BYTE record limit for rank " +
          std::to_string(rank) + ": count=" + std::to_string(record_count) +
          " displacement=" + std::to_string(record_displacement) +
          " record_size=" + std::to_string(record_size) +
          " max_records=" + std::to_string(record_limit));
    }

    const std::uint64_t byte_count = count_u64 * record_size_u64;
    const std::uint64_t byte_displacement = displacement_u64 * record_size_u64;
    if (byte_count > max_mpi_int || byte_displacement > max_mpi_int ||
        byte_count > max_mpi_int - byte_displacement) {
      throw std::invalid_argument(
          std::string(context) + " byte count/displacement exceeds MPI int range for rank " +
          std::to_string(rank) + ": byte_count=" + std::to_string(byte_count) +
          " byte_displacement=" + std::to_string(byte_displacement));
    }
    byte_counts[rank] = static_cast<int>(byte_count);
    byte_displacements[rank] = static_cast<int>(byte_displacement);
  }
}
#endif

[[nodiscard]] PmAxisStencil1d makeAxisStencil(double grid_position, PmAssignmentScheme scheme) {
  PmAxisStencil1d stencil{};
  switch (scheme) {
    case PmAssignmentScheme::kCic: {
      const double base = std::floor(grid_position);
      const std::ptrdiff_t i0 = static_cast<std::ptrdiff_t>(base);
      const double t = grid_position - base;
      stencil.offsets = {i0, i0 + 1, 0};
      stencil.weights = {1.0 - t, t, 0.0};
      stencil.count = 2;
      return stencil;
    }
    case PmAssignmentScheme::kTsc: {
      const double center = std::floor(grid_position + 0.5);
      const std::ptrdiff_t ic = static_cast<std::ptrdiff_t>(center);
      const double delta = grid_position - center;
      const double w_m1 = 0.5 * std::pow(0.5 - delta, 2.0);
      const double w_0 = 0.75 - delta * delta;
      const double w_p1 = 0.5 * std::pow(0.5 + delta, 2.0);
      stencil.offsets = {ic - 1, ic, ic + 1};
      stencil.weights = {w_m1, w_0, w_p1};
      stencil.count = 3;
      return stencil;
    }
  }
  throw std::invalid_argument("Unknown PM assignment scheme in makeAxisStencil");
}

struct PmXPlaneGroup {
  int destination_rank = -1;
  std::uint32_t plane_count = 0U;
  std::array<std::uint32_t, 3> global_ix{};
  std::array<double, 3> x_weight{};
};

[[nodiscard]] std::size_t makePmXPlaneGroups(
    const PmAxisStencil1d& sx,
    bool periodic,
    const PmGridShape& shape,
    const PmDecompositionView& decomposition_view,
    std::array<PmXPlaneGroup, 3>& groups) {
  groups = {};
  std::size_t group_count = 0U;
  for (std::size_t dx = 0; dx < sx.count; ++dx) {
    if (!periodic &&
        (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<std::ptrdiff_t>(shape.nx))) {
      continue;
    }
    const std::size_t ix = periodic
        ? wrapIndex(sx.offsets[dx], shape.nx)
        : static_cast<std::size_t>(sx.offsets[dx]);
    const int destination_rank = decomposition_view.realOwnerRank(PmGlobalCell{ix, 0U, 0U});
    std::size_t group_index = group_count;
    for (std::size_t i = 0; i < group_count; ++i) {
      if (groups[i].destination_rank == destination_rank) {
        group_index = i;
        break;
      }
    }
    if (group_index == group_count) {
      if (group_count >= groups.size()) {
        throw std::logic_error("PM x-stencil generated more destination groups than x planes");
      }
      groups[group_index].destination_rank = destination_rank;
      ++group_count;
    }
    auto& group = groups[group_index];
    if (group.plane_count >= group.global_ix.size()) {
      throw std::logic_error("PM x-plane destination group exceeded fixed stencil capacity");
    }
    const std::size_t slot = static_cast<std::size_t>(group.plane_count);
    group.global_ix[slot] = static_cast<std::uint32_t>(ix);
    group.x_weight[slot] = sx.weights[dx];
    ++group.plane_count;
  }
  return group_count;
}

template <typename ByteVector>
void resizePmWireBufferBounded(
    ByteVector& buffer,
    std::size_t required_bytes,
    std::size_t policy_limit_bytes,
    std::string_view context) {
  if (required_bytes > policy_limit_bytes) {
    throw std::runtime_error(
        std::string(context) + " requires " + std::to_string(required_bytes) +
        " bytes, exceeding aggregate PM routing buffer policy " +
        std::to_string(policy_limit_bytes));
  }
  if (required_bytes > buffer.capacity() || buffer.capacity() > policy_limit_bytes) {
    // Grow (or shrink after a policy reduction) via a fresh exact-size vector
    // rather than resize() growth. This avoids implementation growth heuristics
    // retaining up to ~2x the logical payload and makes owned-capacity high water
    // track the configured routing policy.
    ByteVector replacement(buffer.get_allocator());
    replacement.resize(required_bytes);
    buffer.swap(replacement);
  } else {
    buffer.resize(required_bytes);
  }
  if (buffer.capacity() > policy_limit_bytes) {
    throw std::runtime_error(
        std::string(context) + " retained vector capacity exceeds the PM routing policy");
  }
}

[[nodiscard]] std::size_t pmRoutingParticlesPerRound(
    std::uint64_t per_peer_batch_bytes,
    std::size_t wire_record_bytes) {
  if (wire_record_bytes == 0U || per_peer_batch_bytes < wire_record_bytes) {
    throw std::invalid_argument(
        "PM routing exchange batch must hold at least one fixed wire record");
  }
  const std::uint64_t records = per_peer_batch_bytes / wire_record_bytes;
  return static_cast<std::size_t>(std::min<std::uint64_t>(
      records,
      static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())));
}

[[nodiscard]] int assignmentWindowExponent(PmAssignmentScheme scheme) {
  switch (scheme) {
    case PmAssignmentScheme::kCic:
      return 2;
    case PmAssignmentScheme::kTsc:
      return 3;
  }
  throw std::invalid_argument("Unknown PM assignment scheme in assignmentWindowExponent");
}

[[nodiscard]] std::uint64_t bytesForGridSweep(std::size_t cell_count) {
  const std::size_t bytes = core::checkedSizeMultiply(
      cell_count, sizeof(double), "PM grid sweep bytes");
  return core::checkedIntegralNarrow<std::uint64_t>(bytes, "PM grid sweep bytes");
}

[[nodiscard]] std::uint64_t bytesForParticles(std::size_t particle_count) {
  const std::size_t bytes_per_particle = core::checkedSizeMultiply(
      sizeof(double), 4U, "PM particle profile bytes per particle");
  const std::size_t bytes = core::checkedSizeMultiply(
      particle_count, bytes_per_particle, "PM particle profile bytes");
  return core::checkedIntegralNarrow<std::uint64_t>(bytes, "PM particle profile bytes");
}

[[nodiscard]] std::size_t checkedProduct(std::size_t a, std::size_t b, std::string_view context) {
  return core::checkedSizeMultiply(a, b, context);
}

struct PmPlanLocalStorageLayout {
  std::size_t logical_local_complex_cells = 0U;
  std::size_t allocated_local_complex_cells = 0U;
  std::size_t real_element_count = 0U;
  std::size_t real_z_stride = 0U;
  std::size_t transposed_local_ny = 0U;
  std::size_t transposed_begin_y = 0U;
  bool is_distributed = false;
  bool spectral_transposed = false;
  bool used_backend_allocation_query = false;
};

[[nodiscard]] PmPlanLocalStorageLayout determinePmPlanLocalStorageLayout(
    PmGridShape shape,
    const parallel::PmSlabLayout& layout,
    core::PmDecompositionMode decomposition_mode,
    bool require_active_distributed_backend) {
  if (!shape.isValid()) {
    throw std::invalid_argument("PM plan storage layout requires a valid grid shape");
  }
  if (!layout.isValid() ||
      layout.global_nx != shape.nx ||
      layout.global_ny != shape.ny ||
      layout.global_nz != shape.nz) {
    throw std::invalid_argument(
        "PM plan storage layout requires a slab layout matching the PM grid shape");
  }

  PmPlanLocalStorageLayout storage;
  const std::size_t nz_complex = shape.nz / 2U + 1U;
  storage.logical_local_complex_cells = checkedProduct(
      checkedProduct(layout.local_nx(), shape.ny, "PM plan local Fourier extent"),
      nz_complex,
      "PM plan local Fourier extent");
  storage.allocated_local_complex_cells = storage.logical_local_complex_cells;
  // In-place r2c storage: every real row is padded to 2 (nz/2 + 1) doubles
  // and the real view aliases the complex allocation (2 doubles per cell).
  storage.real_z_stride = checkedProduct(2U, nz_complex, "PM plan real stride");

  if (layout.world_size <= 1) {
    storage.real_element_count = checkedProduct(
        2U, storage.allocated_local_complex_cells, "PM plan serial padded real extent");
    return storage;
  }

  storage.is_distributed = true;
  storage.real_z_stride = checkedProduct(2U, nz_complex, "PM plan distributed real stride");

#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  if (queryActiveMpiWorld(mpi_world_size, mpi_world_rank)) {
    if (mpi_world_size != layout.world_size || mpi_world_rank != layout.world_rank) {
      throw std::invalid_argument(
          "PM slab layout world metadata must match active MPI_COMM_WORLD for distributed FFT storage");
    }

    configureFftwPlannerThreads();
    fftw_mpi_init();
    ptrdiff_t backend_local_nx = 0;
    ptrdiff_t backend_begin_x = 0;
    ptrdiff_t backend_local_ny = 0;
    ptrdiff_t backend_begin_y = 0;
    ptrdiff_t backend_alloc_local = 0;
    if (decomposition_mode == core::PmDecompositionMode::kPencil) {
      backend_alloc_local = fftw_mpi_local_size_3d_transposed(
          static_cast<ptrdiff_t>(shape.nx),
          static_cast<ptrdiff_t>(shape.ny),
          static_cast<ptrdiff_t>(nz_complex),
          MPI_COMM_WORLD,
          &backend_local_nx,
          &backend_begin_x,
          &backend_local_ny,
          &backend_begin_y);
    } else {
      backend_alloc_local = fftw_mpi_local_size_3d(
          static_cast<ptrdiff_t>(shape.nx),
          static_cast<ptrdiff_t>(shape.ny),
          static_cast<ptrdiff_t>(nz_complex),
          MPI_COMM_WORLD,
          &backend_local_nx,
          &backend_begin_x);
    }
    if (backend_alloc_local < 0) {
      throw std::runtime_error(
          "FFTW MPI reported a negative local allocation size for distributed PM storage");
    }
    const bool local_extent_mismatch =
        backend_local_nx != static_cast<ptrdiff_t>(layout.local_nx());
    // FFTW leaves local_0_start unspecified for a zero-width input slab. The
    // origin carries no ownership meaning in that legal case.
    const bool nonempty_origin_mismatch =
        backend_local_nx > 0 &&
        backend_begin_x != static_cast<ptrdiff_t>(layout.owned_x.begin_x);
    if (local_extent_mismatch || nonempty_origin_mismatch) {
      std::ostringstream message;
      message << "PM slab layout is incompatible with FFTW MPI ownership for this communicator: rank="
              << layout.world_rank << ", configured=[" << layout.owned_x.begin_x << ','
              << layout.owned_x.end_x << "), backend=[" << backend_begin_x << ','
              << (backend_begin_x + backend_local_nx) << ')';
      throw std::invalid_argument(message.str());
    }

    storage.allocated_local_complex_cells = std::max<std::size_t>(
        1U, static_cast<std::size_t>(backend_alloc_local));
    if (storage.allocated_local_complex_cells < storage.logical_local_complex_cells) {
      throw std::runtime_error(
          "FFTW MPI local allocation is smaller than CHUI's logical PM Fourier extent");
    }
    if (decomposition_mode == core::PmDecompositionMode::kPencil) {
      if (backend_local_ny < 0 || backend_begin_y < 0) {
        throw std::runtime_error(
            "FFTW MPI reported negative transposed PM ownership metadata");
      }
      storage.spectral_transposed = true;
      storage.transposed_local_ny = static_cast<std::size_t>(backend_local_ny);
      storage.transposed_begin_y = static_cast<std::size_t>(backend_begin_y);
      storage.logical_local_complex_cells = checkedProduct(
          checkedProduct(
              storage.transposed_local_ny, shape.nx,
              "PM plan transposed logical spectral extent"),
          nz_complex,
          "PM plan transposed logical spectral extent");
      if (storage.allocated_local_complex_cells <
          storage.logical_local_complex_cells) {
        throw std::runtime_error(
            "FFTW MPI local allocation is smaller than CHUI's transposed logical PM Fourier extent");
      }
    }
    storage.used_backend_allocation_query = true;
  }
#endif

  if (!storage.used_backend_allocation_query) {
    if (require_active_distributed_backend) {
      throw std::runtime_error(
          "distributed PM plan storage requires an active FFTW-MPI session");
    }
    // Conservative backend-independent fallback for preflight/compile-only
    // environments. A transposed FFT can require the larger of input x-slab
    // or output y-slab spectral storage.
    std::size_t conservative_complex = storage.logical_local_complex_cells;
    if (decomposition_mode == core::PmDecompositionMode::kPencil) {
      const parallel::PmSlabRange y_range = parallel::pmOwnedXRangeForRank(
          shape.ny, layout.world_size, layout.world_rank);
      const std::size_t transposed_complex = checkedProduct(
          checkedProduct(shape.nx, y_range.extentX(),
                         "PM plan transposed Fourier extent"),
          nz_complex,
          "PM plan transposed Fourier extent");
      conservative_complex = std::max(conservative_complex, transposed_complex);
      storage.spectral_transposed = true;
      storage.transposed_local_ny = y_range.extentX();
      storage.transposed_begin_y = y_range.begin_x;
      storage.logical_local_complex_cells = transposed_complex;
    }
    storage.allocated_local_complex_cells =
        std::max<std::size_t>(1U, conservative_complex);
  }

  storage.real_element_count = checkedProduct(
      2U,
      storage.allocated_local_complex_cells,
      "PM plan distributed padded real extent");
  return storage;
}

[[nodiscard]] std::uint64_t checkedBytesForCount(
    std::size_t count,
    std::size_t element_size,
    std::string_view context) {
  if (element_size != 0U && count > std::numeric_limits<std::uint64_t>::max() / element_size) {
    throw std::overflow_error(std::string(context) + " byte estimate overflows uint64_t");
  }
  return static_cast<std::uint64_t>(count) * static_cast<std::uint64_t>(element_size);
}

[[nodiscard]] std::uint64_t checkedAddBytes(
    std::uint64_t a,
    std::uint64_t b,
    std::string_view context) {
  if (b > std::numeric_limits<std::uint64_t>::max() - a) {
    throw std::overflow_error(std::string(context) + " byte estimate overflows uint64_t");
  }
  return a + b;
}

// Workload-aware routing (2026-10-25). Round sizing derives from a
// rank-independent model (configured batch bytes, world size and the exact
// per-kind routing metadata formula), so every rank partitions its local
// population into identically sized rounds. The live wire buffers and the
// communication arena are then sized from the measured remote-record
// workload instead of (P-1) x effective_per_peer_payload_bytes.
enum class PmRoutingPhaseKind : std::uint8_t { kDensity, kForceInterpolation };

[[nodiscard]] std::size_t pmRoutingWireRecordBytes(PmRoutingPhaseKind kind) noexcept {
  return kind == PmRoutingPhaseKind::kDensity ? k_pm_density_plane_wire_bytes
                                              : k_pm_plane_interpolation_request_wire_bytes;
}

// Exactly the int count/displacement lanes plus the size_t cursor lane that
// DensityExchangeBuffers (8 int lanes) and PlaneInterpolationExchangeBuffers
// (12 int lanes) allocate per rank.
[[nodiscard]] std::uint64_t pmRoutingMetadataBytes(PmRoutingPhaseKind kind, int world_size) {
  const std::uint64_t int_lanes = kind == PmRoutingPhaseKind::kDensity ? 8U : 12U;
  return checkedBytesForCount(
      static_cast<std::size_t>(std::max(world_size, 1)),
      static_cast<std::size_t>(int_lanes * sizeof(int) + sizeof(std::size_t)),
      "PM routing metadata formula");
}

[[nodiscard]] PmRoutingCapacityModel pmPlannedRoutingModel(
    PmRoutingPhaseKind kind, int world_size, std::uint64_t configured_per_peer_bytes) {
  return modelPmRoutingCapacity(
      world_size, configured_per_peer_bytes, k_pm_routing_modeled_workspace_limit_bytes,
      pmRoutingMetadataBytes(kind, world_size), pmRoutingWireRecordBytes(kind));
}

struct PmRemotePlaneWorkload {
  // Largest record count sent to any single remote peer in any round.
  std::uint64_t max_peer_records = 0U;
  // Largest total remote record count of any round (this rank's send side).
  std::uint64_t max_round_records = 0U;
};

// Counts remote x-plane groups exactly as pass 1 of the density and force
// routing loops does. grid_x(p) returns the target's x coordinate in grid
// units, or std::nullopt for entries the routing loop skips.
template <typename GridX>
[[nodiscard]] PmRemotePlaneWorkload measurePmRemotePlaneWorkload(
    std::size_t count,
    std::size_t particles_per_round,
    int world_rank,
    int world_size,
    bool periodic,
    PmAssignmentScheme scheme,
    const PmGridShape& shape,
    const PmDecompositionView& decomposition_view,
    const GridX& grid_x) {
  if (particles_per_round == 0U) {
    throw std::invalid_argument("PM routing workload requires a positive round capacity");
  }
  PmRemotePlaneWorkload workload;
  std::vector<std::uint64_t> per_peer(static_cast<std::size_t>(std::max(world_size, 1)), 0U);
  for (std::size_t begin = 0; begin < count; begin += particles_per_round) {
    const std::size_t end = std::min(count, begin + particles_per_round);
    std::fill(per_peer.begin(), per_peer.end(), std::uint64_t{0U});
    std::uint64_t round_records = 0U;
    for (std::size_t p = begin; p < end; ++p) {
      const std::optional<double> x = grid_x(p);
      if (!x.has_value()) continue;
      const PmAxisStencil1d sx = makeAxisStencil(*x, scheme);
      std::array<PmXPlaneGroup, 3> groups{};
      const std::size_t group_count = makePmXPlaneGroups(sx, periodic, shape, decomposition_view, groups);
      for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
        const int destination_rank = groups[group_i].destination_rank;
        if (destination_rank == world_rank) continue;
        ++per_peer[static_cast<std::size_t>(destination_rank)];
        ++round_records;
      }
    }
    workload.max_round_records = std::max(workload.max_round_records, round_records);
    for (const std::uint64_t records : per_peer) {
      workload.max_peer_records = std::max(workload.max_peer_records, records);
    }
  }
  return workload;
}

// Bounded wire reservations and phase arena bytes from a collectively agreed
// workload: send <= this rank's largest round, receive <= (P-1) x the global
// largest per-peer round payload; both never exceed the configured model.
struct PmRoutingReservation {
  std::uint64_t send_bytes = 0U;
  std::uint64_t receive_bytes = 0U;
  std::uint64_t phase_arena_bytes = 0U;
};

[[nodiscard]] PmRoutingReservation pmRoutingReservation(
    PmRoutingPhaseKind kind,
    const PmRoutingCapacityModel& model,
    int world_size,
    std::uint64_t local_max_round_records,
    std::uint64_t global_max_peer_records) {
  const std::size_t record_bytes = pmRoutingWireRecordBytes(kind);
  const std::uint64_t remote_peers = static_cast<std::uint64_t>(std::max(world_size - 1, 0));
  PmRoutingReservation reservation;
  reservation.send_bytes = std::min<std::uint64_t>(
      model.max_send_payload_bytes,
      checkedBytesForCount(static_cast<std::size_t>(local_max_round_records), record_bytes,
                           "PM routing send reservation"));
  const std::uint64_t receive_records = global_max_peer_records > 0U &&
          remote_peers > std::numeric_limits<std::uint64_t>::max() / global_max_peer_records
      ? std::numeric_limits<std::uint64_t>::max()
      : remote_peers * global_max_peer_records;
  reservation.receive_bytes = receive_records >
          model.max_receive_payload_bytes / static_cast<std::uint64_t>(record_bytes)
      ? model.max_receive_payload_bytes
      : receive_records * static_cast<std::uint64_t>(record_bytes);
  reservation.phase_arena_bytes = checkedAddBytes(
      checkedAddBytes(pmRoutingMetadataBytes(kind, world_size), k_pm_routing_workspace_headroom_bytes,
                      "PM routing phase arena bytes"),
      checkedAddBytes(reservation.send_bytes, reservation.receive_bytes, "PM routing phase arena bytes"),
      "PM routing phase arena bytes");
  return reservation;
}

[[nodiscard]] std::uint64_t estimateIsolatedRootWorkspaceBytes(
    const PmGridShape& shape,
    bool distributed_slabs) {
  const std::size_t physical_cells =
      checkedProduct(checkedProduct(shape.nx, shape.ny, "isolated PM physical grid"),
                     shape.nz,
                     "isolated PM physical grid");
  const std::size_t pad_nx = checkedProduct(2U, shape.nx, "isolated PM padded nx");
  const std::size_t pad_ny = checkedProduct(2U, shape.ny, "isolated PM padded ny");
  const std::size_t pad_nz = checkedProduct(2U, shape.nz, "isolated PM padded nz");
  const std::size_t padded_cells =
      checkedProduct(checkedProduct(pad_nx, pad_ny, "isolated PM padded grid"),
                     pad_nz,
                     "isolated PM padded grid");
  std::uint64_t bytes = 0;
  if (distributed_slabs) {
    bytes = checkedAddBytes(
        bytes,
        checkedBytesForCount(physical_cells, sizeof(double), "isolated PM root density gather"),
        "isolated PM root workspace");
  }
  bytes = checkedAddBytes(
      bytes,
      checkedBytesForCount(padded_cells, sizeof(std::complex<double>) * 3U, "isolated PM padded complex workspace"),
      "isolated PM root workspace");
  bytes = checkedAddBytes(
      bytes,
      checkedBytesForCount(padded_cells, sizeof(double), "isolated PM padded real workspace"),
      "isolated PM root workspace");
  bytes = checkedAddBytes(
      bytes,
      checkedBytesForCount(physical_cells, sizeof(double) * 4U, "isolated PM extracted fields"),
      "isolated PM root workspace");
  return bytes;
}

[[nodiscard]] std::string isolatedPmGuardMessage(
    const PmGridShape& shape,
    int world_size,
    int world_rank,
    std::uint64_t estimated_bytes,
    std::uint64_t limit_bytes) {
  return "isolated/open PM root-gather guard exceeded on rank " + std::to_string(world_rank) +
      ": grid_shape=(" + std::to_string(shape.nx) + "," + std::to_string(shape.ny) + "," +
      std::to_string(shape.nz) + ") ranks=" + std::to_string(world_size) +
      " estimated_root_bytes=" + std::to_string(estimated_bytes) +
      " configured_limit_bytes=" + std::to_string(limit_bytes) +
      " route=root_gather_scatter policy=bounded_small_grid_only";
}

struct BoxLengths {
  double lx = 0.0;
  double ly = 0.0;
  double lz = 0.0;
};

[[nodiscard]] BoxLengths effectiveBoxLengths(const PmSolveOptions& options) {
  const double scalar = options.box_size_mpc_comoving;
  const double lx = options.box_size_x_mpc_comoving > 0.0 ? options.box_size_x_mpc_comoving : scalar;
  const double ly = options.box_size_y_mpc_comoving > 0.0 ? options.box_size_y_mpc_comoving : scalar;
  const double lz = options.box_size_z_mpc_comoving > 0.0 ? options.box_size_z_mpc_comoving : scalar;
  return BoxLengths{.lx = lx, .ly = ly, .lz = lz};
}

// Semantic identity of a routing workload plan: everything except the
// particle coordinates that determines round partitioning, x-plane
// destinations and wire record sizes. Coordinates are deliberately not keyed:
// consumers re-measure the current workload and require it to be bounded by
// the plan (see agreePmRoutingWorkload), which is strictly stronger than any
// pointer/count/generation key and needs no source-generation plumbing.
struct PmRoutingPlanKey {
  int world_size = 0;
  int world_rank = 0;
  std::size_t nx = 0U;
  std::size_t ny = 0U;
  std::size_t nz = 0U;
  std::size_t owned_x_begin = 0U;
  std::size_t owned_x_end = 0U;
  std::uint64_t configured_batch_bytes = 0U;
  PmAssignmentScheme assignment_scheme = PmAssignmentScheme::kCic;
  PmBoundaryCondition boundary_condition = PmBoundaryCondition::kPeriodic;
  core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;
  std::uint64_t box_x_bits = 0U;
  std::uint64_t box_y_bits = 0U;
  std::uint64_t box_z_bits = 0U;
  [[nodiscard]] bool operator==(const PmRoutingPlanKey&) const = default;
};

[[nodiscard]] PmRoutingPlanKey makePmRoutingPlanKey(
    const PmGridShape& shape, const parallel::PmSlabLayout& layout, const PmSolveOptions& options) {
  const BoxLengths lengths = effectiveBoxLengths(options);
  return PmRoutingPlanKey{
      .world_size = layout.world_size,
      .world_rank = layout.world_rank,
      .nx = shape.nx,
      .ny = shape.ny,
      .nz = shape.nz,
      .owned_x_begin = layout.owned_x.begin_x,
      .owned_x_end = layout.owned_x.end_x,
      .configured_batch_bytes = options.routing_exchange_batch_bytes,
      .assignment_scheme = options.assignment_scheme,
      .boundary_condition = options.boundary_condition,
      .decomposition_mode = options.decomposition_mode,
      .box_x_bits = std::bit_cast<std::uint64_t>(lengths.lx),
      .box_y_bits = std::bit_cast<std::uint64_t>(lengths.ly),
      .box_z_bits = std::bit_cast<std::uint64_t>(lengths.lz),
  };
}

// Non-throwing plan check run inside a routing call's workload phase. The
// plan bounds the current call iff the key matches and the freshly measured
// local workload does not exceed the planned values (count, per-peer maximum,
// largest local round). A violated plan is reported through the agreement
// reduction so every rank fails together.
[[nodiscard]] bool pmRoutingPlanBoundsWorkload(
    const PmRoutingPlanKey& planned_key,
    std::uint64_t planned_global_max_count,
    std::uint64_t planned_global_max_peer_records,
    std::uint64_t planned_local_max_round_records,
    const PmRoutingPlanKey& current_key,
    std::uint64_t current_count,
    const PmRemotePlaneWorkload& current) noexcept {
  return planned_key == current_key && current_count <= planned_global_max_count &&
      current.max_peer_records <= planned_global_max_peer_records &&
      current.max_round_records <= planned_local_max_round_records;
}

#if COSMOSIM_ENABLE_MPI
struct PmRoutingAgreement {
  std::uint64_t global_max_count = 0U;
  std::uint64_t global_max_peer_records = 0U;
  bool planned = false;
};

// Routing-call agreement (collective on MPI_COMM_WORLD; exactly one
// MPI_Allreduce(MAX) on every rank, reached only after the preceding
// coordinated workload phase succeeded everywhere). Lanes:
//   [count, M_local, plan_present, plan_absent, plan_violated, world_size - violating_rank]
// Plan presence must be rank-uniform (planning is collective); a mixed
// presence or any violated plan throws identically on every rank before any
// routing buffer is leased. The agreed count/M come from the *current*
// measurement, so reservations are always computed from present demand; a
// valid plan only certifies that the communication arena was configured
// large enough before the lease.
[[nodiscard]] PmRoutingAgreement agreePmRoutingWorkload(
    std::string_view stage,
    int world_rank,
    int world_size,
    std::uint64_t local_count,
    std::uint64_t local_max_peer_records,
    bool plan_present,
    bool plan_violated,
    double& accumulated_mpi_wait_ms) {
  const std::array<std::uint64_t, 6U> local{
      local_count,
      local_max_peer_records,
      plan_present ? 1U : 0U,
      plan_present ? 0U : 1U,
      plan_violated ? 1U : 0U,
      plan_violated ? static_cast<std::uint64_t>(world_size - world_rank) : 0U};
  std::array<std::uint64_t, 6U> global{};
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    requirePmMpiSuccess(
        MPI_Allreduce(local.data(), global.data(), 6, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD),
        std::string(stage) + " routing agreement MPI_Allreduce");
  });
  if (global[2] != 0U && global[3] != 0U) {
    throw std::runtime_error(
        std::string(stage) + " routing workload plan presence differs across ranks; plan*Routing() "
        "must be called collectively before the routing call");
  }
  if (global[4] != 0U) {
    const std::uint64_t first_rank = static_cast<std::uint64_t>(world_size) - global[5];
    throw std::runtime_error(
        std::string(stage) + " routing workload plan does not bound the current workload (coordinates, "
        "decomposition or routing options changed after planning); first violating rank=" +
        std::to_string(first_rank) + (plan_violated ? " (this rank)" : " (peer rank)"));
  }
  return PmRoutingAgreement{
      .global_max_count = global[0], .global_max_peer_records = global[1], .planned = global[2] != 0U};
}
#endif

void validateOptions(const PmGridShape& shape, const PmSolveOptions& options) {
  switch (options.assignment_scheme) {
    case PmAssignmentScheme::kCic:
    case PmAssignmentScheme::kTsc:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid assignment_scheme");
  }
  switch (options.boundary_condition) {
    case PmBoundaryCondition::kPeriodic:
    case PmBoundaryCondition::kIsolatedOpen:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid boundary_condition");
  }
  switch (options.data_residency) {
    case PmDataResidencyPolicy::kHostOnly:
    case PmDataResidencyPolicy::kPreferDevice:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid data_residency");
  }
  switch (options.execution_policy) {
    case core::ExecutionPolicy::kHostSerial:
    case core::ExecutionPolicy::kCuda:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid execution_policy");
  }
  switch (options.decomposition_mode) {
    case core::PmDecompositionMode::kSlab:
    case core::PmDecompositionMode::kPencil:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid decomposition_mode");
  }
  if (!shape.isValid()) {
    throw std::invalid_argument("PM grid shape must be non-zero in all dimensions");
  }
  if (!std::isfinite(options.box_size_mpc_comoving) ||
      !std::isfinite(options.box_size_x_mpc_comoving) ||
      !std::isfinite(options.box_size_y_mpc_comoving) ||
      !std::isfinite(options.box_size_z_mpc_comoving)) {
    throw std::invalid_argument("PM solve requires finite scalar and axis-aware box lengths");
  }
  const BoxLengths lengths = effectiveBoxLengths(options);
  if (!std::isfinite(lengths.lx) || !std::isfinite(lengths.ly) || !std::isfinite(lengths.lz) ||
      lengths.lx <= 0.0 || lengths.ly <= 0.0 || lengths.lz <= 0.0) {
    throw std::invalid_argument("PM solve requires finite positive effective box lengths on all axes");
  }
  if (!std::isfinite(options.scale_factor) || options.scale_factor <= 0.0) {
    throw std::invalid_argument("PM solve requires finite scale_factor > 0");
  }
  if (!std::isfinite(options.gravitational_constant_code) || options.gravitational_constant_code <= 0.0) {
    throw std::invalid_argument("PM solve requires finite gravitational_constant_code > 0");
  }
  if (!std::isfinite(options.tree_pm_split_scale_comoving)) {
    throw std::invalid_argument("PM solve requires a finite TreePM split scale");
  }
  if (options.routing_exchange_batch_bytes <
      std::max(k_pm_density_plane_wire_bytes, k_pm_plane_interpolation_request_wire_bytes)) {
    throw std::invalid_argument(
        "PM solve routing_exchange_batch_bytes is too small for one routing record");
  }
  if (options.execution_policy == core::ExecutionPolicy::kCuda && options.data_residency == PmDataResidencyPolicy::kHostOnly) {
    throw std::invalid_argument(
        "execution_policy=cuda requires data_residency=kPreferDevice for explicit host/device ownership");
  }
  if (options.execution_policy == core::ExecutionPolicy::kCuda && options.assignment_scheme != PmAssignmentScheme::kCic) {
    throw std::invalid_argument(
        "execution_policy=cuda currently supports only assignment_scheme=cic in this build");
  }
  if (options.spectral_source_hook != nullptr &&
      (options.execution_policy == core::ExecutionPolicy::kCuda ||
       options.boundary_condition != PmBoundaryCondition::kPeriodic)) {
    throw std::invalid_argument(
        "PM spectral_source_hook (N03 neutrino source) is supported only by the periodic host Poisson operator");
  }
  if (options.execution_policy == core::ExecutionPolicy::kCuda &&
      options.boundary_condition != PmBoundaryCondition::kPeriodic) {
    throw std::invalid_argument(
        "execution_policy=cuda currently supports only boundary_condition=periodic in this build");
  }
  switch (options.force_gradient_backend) {
    case PmForceGradientBackend::kSpectral:
      break;
    case PmForceGradientBackend::kFiniteDifference4:
      // The isolated/open and CUDA solvers keep their own gradient contracts.
      if (options.boundary_condition != PmBoundaryCondition::kPeriodic ||
          options.execution_policy != core::ExecutionPolicy::kHostSerial) {
        throw std::invalid_argument(
            "PM finite-difference force gradient supports only periodic host solves");
      }
      if (shape.nx < 4U || shape.ny < 4U || shape.nz < 4U) {
        throw std::invalid_argument(
            "PM finite-difference force gradient requires at least four cells per axis");
      }
      break;
    default:
      throw std::invalid_argument("PM solve has invalid force_gradient_backend");
  }
  if (options.boundary_condition == PmBoundaryCondition::kIsolatedOpen &&
      options.enable_window_deconvolution) {
    throw std::invalid_argument("isolated PM currently does not support window deconvolution");
  }
}

void validateSingleRankFullDomainGridContract(const PmGridStorage& grid, std::string_view callsite) {
  if (!grid.ownsFullDomain()) {
    throw std::invalid_argument(
        std::string(callsite) +
        " requires full-domain PM slab ownership in single-rank mode; use a valid distributed slab layout for multi-rank PM execution");
  }
}

}  // namespace

PmRoutingCapacityModel modelPmRoutingCapacity(
    int world_size,
    std::uint64_t configured_per_peer_max_bytes,
    std::uint64_t aggregate_workspace_limit_bytes,
    std::uint64_t fixed_workspace_bytes,
    std::size_t wire_record_bytes) {
  if (world_size < 1) {
    throw std::invalid_argument("PM routing capacity model requires world_size >= 1");
  }
  if (wire_record_bytes == 0U) {
    throw std::invalid_argument("PM routing capacity model requires a nonzero wire record size");
  }
  const std::uint64_t record_bytes = static_cast<std::uint64_t>(wire_record_bytes);
  if (fixed_workspace_bytes > aggregate_workspace_limit_bytes) {
    throw std::invalid_argument(
        "PM routing fixed workspace exceeds aggregate workspace limit");
  }

  PmRoutingCapacityModel model{
      .configured_per_peer_max_bytes = configured_per_peer_max_bytes,
      .fixed_workspace_bytes = fixed_workspace_bytes,
  };
  if (world_size == 1) {
    model.max_simultaneous_workspace_bytes = fixed_workspace_bytes;
    return model;
  }
  if (configured_per_peer_max_bytes < record_bytes) {
    throw std::invalid_argument(
        "PM routing configured per-peer maximum cannot hold one wire record");
  }

  const std::uint64_t remote_peers = static_cast<std::uint64_t>(world_size - 1);
  if (remote_peers > std::numeric_limits<std::uint64_t>::max() / 2U) {
    throw std::overflow_error("PM routing peer-slot count overflows uint64 capacity");
  }
  const std::uint64_t simultaneous_peer_slots = 2U * remote_peers;
  const std::uint64_t payload_budget =
      aggregate_workspace_limit_bytes - fixed_workspace_bytes;
  const std::uint64_t aggregate_limited_peer_bytes =
      payload_budget / simultaneous_peer_slots;
  std::uint64_t effective_peer_bytes =
      std::min(configured_per_peer_max_bytes, aggregate_limited_peer_bytes);
  effective_peer_bytes -= effective_peer_bytes % record_bytes;
  if (effective_peer_bytes < record_bytes) {
    throw std::invalid_argument(
        "PM routing aggregate workspace cannot hold one record per active peer direction");
  }
  if (remote_peers > std::numeric_limits<std::uint64_t>::max() / effective_peer_bytes) {
    throw std::overflow_error("PM routing directional payload capacity overflows uint64");
  }
  const std::uint64_t directional_payload = remote_peers * effective_peer_bytes;
  if (directional_payload > std::numeric_limits<std::uint64_t>::max() / 2U) {
    throw std::overflow_error("PM routing combined payload capacity overflows uint64");
  }
  const std::uint64_t combined_payload = 2U * directional_payload;
  if (combined_payload > aggregate_workspace_limit_bytes - fixed_workspace_bytes) {
    throw std::logic_error("PM routing capacity model violated aggregate workspace limit");
  }

  model.effective_per_peer_payload_bytes = effective_peer_bytes;
  model.max_send_payload_bytes = directional_payload;
  model.max_receive_payload_bytes = directional_payload;
  model.max_simultaneous_workspace_bytes = fixed_workspace_bytes + combined_payload;
  return model;
}

class PmSolver::Impl {
 public:
  struct PlanKey {
    int world_size = 1;
    int world_rank = 0;
    std::size_t owned_begin_x = 0;
    std::size_t owned_end_x = 0;
    core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;

    [[nodiscard]] bool operator==(const PlanKey& other) const noexcept {
      return world_size == other.world_size && world_rank == other.world_rank && owned_begin_x == other.owned_begin_x &&
          owned_end_x == other.owned_end_x && decomposition_mode == other.decomposition_mode;
    }
  };

  struct PlanKeyHasher {
    [[nodiscard]] std::size_t operator()(const PlanKey& key) const noexcept {
      constexpr std::size_t k_hash_magic = static_cast<std::size_t>(0x9e3779b97f4a7c15ULL);
      std::size_t seed = static_cast<std::size_t>(key.world_size) * static_cast<std::size_t>(1315423911U) +
          static_cast<std::size_t>(key.world_rank);
      seed ^= key.owned_begin_x + k_hash_magic + (seed << 6U) + (seed >> 2U);
      seed ^= key.owned_end_x + k_hash_magic + (seed << 6U) + (seed >> 2U);
      seed ^= static_cast<std::size_t>(key.decomposition_mode) + k_hash_magic + (seed << 6U) + (seed >> 2U);
      return seed;
    }
  };

  struct PlanResources {
    parallel::PmSlabLayout layout{};
    // In-place FFT storage (2026-10-25): `fourier` is the single owner of the
    // FFT working set. The real-space view aliases the same bytes as 2 *
    // fourier.size() doubles in FFTW's padded r2c layout (rows of
    // real_z_stride = 2 (nz/2 + 1) doubles; serial and distributed alike).
    // std::complex<double> is layout-compatible with double[2] ([complex.numbers]).
    // Forward r2c and inverse c2r transforms run in place; the real view and
    // the spectrum are never simultaneously meaningful.
    std::vector<std::complex<double>> fourier;
    [[nodiscard]] std::span<double> realView() noexcept {
      return {reinterpret_cast<double*>(fourier.data()), 2U * fourier.size()};
    }
    // Offset of the first physical z cell of (local_ix, iy) in the padded real
    // view. Physical cells are [base, base + nz); [base + nz, base +
    // real_z_stride) is r2c padding and is never read as mesh data.
    [[nodiscard]] std::size_t realRowBase(std::size_t local_ix, std::size_t iy) const noexcept {
      return (local_ix * layout.global_ny + iy) * real_z_stride;
    }
    // Pristine phi_k copy for the spectral gradient backend only (complex-to-
    // real transforms destroy their input). Allocated on the first spectral
    // solve and retained; never allocated by kFiniteDifference4.
    std::vector<std::complex<double>> potential_k;
    // kFiniteDifference4 distributed x halo planes: [left -2, left -1,
    // right +1, right +2] followed by two send-staging planes; each ny*nz doubles.
    std::vector<double> fd_halo_planes;
    // M48-09 derives the Poisson operator from small per-axis metadata rather
    // than retaining one scalar per allocated spectral coefficient.
    std::vector<double> kx_by_ix;
    std::vector<double> ky_by_iy;
    std::vector<double> kz_by_iz;
    std::vector<double> window_x_by_ix;
    std::vector<double> window_y_by_iy;
    std::vector<double> window_z_by_iz;
#if COSMOSIM_ENABLE_FFTW
    fftw_plan forward_plan = nullptr;
    fftw_plan inverse_plan = nullptr;
#endif
    bool is_distributed = false;
    bool spectral_transposed = false;
    std::size_t real_z_stride = 0;
    std::size_t transposed_local_ny = 0;
    std::size_t transposed_begin_y = 0;
    std::size_t logical_local_complex_cells = 0U;
    std::size_t allocated_local_complex_cells = 0U;
    bool spectral_axis_metadata_ready = false;
    double cached_lx = 0.0;
    double cached_ly = 0.0;
    double cached_lz = 0.0;
    bool cached_window_deconvolution = false;
    PmAssignmentScheme cached_assignment_scheme = PmAssignmentScheme::kCic;

    // Direct periodic-density lifecycle. The pointer is identity only; the PM
    // plan owns the bytes. Forward FFT invalidates this authority.
    const PmGridStorage* periodic_density_grid = nullptr;
    std::uint64_t periodic_density_generation = 0U;
    std::uint64_t periodic_density_grid_generation = 0U;
    bool periodic_density_ready = false;
  };

  struct DensityExchangeBuffers {
    explicit DensityExchangeBuffers(std::pmr::memory_resource* resource, std::size_t ranks)
        : send_wire(resource), recv_wire(resource), send_counts(resource), send_displs(resource),
          recv_counts(resource), recv_displs(resource), send_counts_bytes(resource),
          send_displs_bytes(resource), recv_counts_bytes(resource), recv_displs_bytes(resource),
          cursor(resource) {
      send_counts.assign(ranks, 0);
      send_displs.assign(ranks, 0);
      recv_counts.assign(ranks, 0);
      recv_displs.assign(ranks, 0);
      send_counts_bytes.assign(ranks, 0);
      send_displs_bytes.assign(ranks, 0);
      recv_counts_bytes.assign(ranks, 0);
      recv_displs_bytes.assign(ranks, 0);
      cursor.assign(ranks, 0U);
    }

    std::uint64_t workspace_high_water_bytes = 0U;
    std::pmr::vector<std::uint8_t> send_wire;
    std::pmr::vector<std::uint8_t> recv_wire;
    std::pmr::vector<int> send_counts;
    std::pmr::vector<int> send_displs;
    std::pmr::vector<int> recv_counts;
    std::pmr::vector<int> recv_displs;
    std::pmr::vector<int> send_counts_bytes;
    std::pmr::vector<int> send_displs_bytes;
    std::pmr::vector<int> recv_counts_bytes;
    std::pmr::vector<int> recv_displs_bytes;
    std::pmr::vector<std::size_t> cursor;
  };

  struct PlaneInterpolationExchangeBuffers {
    explicit PlaneInterpolationExchangeBuffers(std::pmr::memory_resource* resource, std::size_t ranks)
        : send_wire(resource), recv_wire(resource), send_counts(resource), send_displs(resource),
          recv_counts(resource), recv_displs(resource), send_counts_bytes(resource),
          send_displs_bytes(resource), recv_counts_bytes(resource), recv_displs_bytes(resource),
          send_response_counts_bytes(resource), send_response_displs_bytes(resource),
          recv_response_counts_bytes(resource), recv_response_displs_bytes(resource), cursor(resource) {
      send_counts.assign(ranks, 0);
      send_displs.assign(ranks, 0);
      recv_counts.assign(ranks, 0);
      recv_displs.assign(ranks, 0);
      send_counts_bytes.assign(ranks, 0);
      send_displs_bytes.assign(ranks, 0);
      recv_counts_bytes.assign(ranks, 0);
      recv_displs_bytes.assign(ranks, 0);
      send_response_counts_bytes.assign(ranks, 0);
      send_response_displs_bytes.assign(ranks, 0);
      recv_response_counts_bytes.assign(ranks, 0);
      recv_response_displs_bytes.assign(ranks, 0);
      cursor.assign(ranks, 0U);
    }

    std::uint64_t workspace_high_water_bytes = 0U;
    // Two physical wire buffers are reused across request/response phases.
    std::pmr::vector<std::uint8_t> send_wire;
    std::pmr::vector<std::uint8_t> recv_wire;
    std::pmr::vector<int> send_counts;
    std::pmr::vector<int> send_displs;
    std::pmr::vector<int> recv_counts;
    std::pmr::vector<int> recv_displs;
    std::pmr::vector<int> send_counts_bytes;
    std::pmr::vector<int> send_displs_bytes;
    std::pmr::vector<int> recv_counts_bytes;
    std::pmr::vector<int> recv_displs_bytes;
    std::pmr::vector<int> send_response_counts_bytes;
    std::pmr::vector<int> send_response_displs_bytes;
    std::pmr::vector<int> recv_response_counts_bytes;
    std::pmr::vector<int> recv_response_displs_bytes;
    std::pmr::vector<std::size_t> cursor;
  };


  struct IndexedTargetWorkspace {
    std::vector<double> compact_ax;
    std::vector<double> compact_ay;
    std::vector<double> compact_az;
  };

#if COSMOSIM_ENABLE_CUDA
  // Persistent owner for the CUDA PM staging backend. The currently validated
  // path still performs FFT/Poisson/gradient on the host, but stream and device
  // storage no longer have per-call ownership. Reserved FFT/filter workspaces
  // and an opaque plan slot make the resident backend boundary explicit without
  // claiming cuFFT execution in builds where it has not been implemented.
  struct CudaBackendWorkspace {
    cudaStream_t stream = nullptr;
    core::DeviceBufferDouble pos_x;
    core::DeviceBufferDouble pos_y;
    core::DeviceBufferDouble pos_z;
    core::DeviceBufferDouble mass;
    core::DeviceBufferDouble density;
    core::DeviceBufferDouble force_x;
    core::DeviceBufferDouble force_y;
    core::DeviceBufferDouble force_z;
    core::DeviceBufferDouble accel_x;
    core::DeviceBufferDouble accel_y;
    core::DeviceBufferDouble accel_z;
    core::DeviceBufferDouble fft_workspace;
    core::DeviceBufferDouble green_filter_workspace;
    std::uintptr_t fft_plan_handle = 0U;

    ~CudaBackendWorkspace() {
      if (stream != nullptr) {
        (void)cudaStreamDestroy(stream);
      }
    }

    [[nodiscard]] core::CudaStreamView prepare(
        std::size_t particle_count, std::size_t grid_cell_count) {
      if (stream == nullptr &&
          cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) != cudaSuccess) {
        throw std::runtime_error("Failed to create persistent CUDA PM stream");
      }
      pos_x.resize(particle_count);
      pos_y.resize(particle_count);
      pos_z.resize(particle_count);
      mass.resize(particle_count);
      density.resize(grid_cell_count);
      force_x.resize(grid_cell_count);
      force_y.resize(grid_cell_count);
      force_z.resize(grid_cell_count);
      accel_x.resize(particle_count);
      accel_y.resize(particle_count);
      accel_z.resize(particle_count);
      return core::CudaStreamView{stream};
    }
  };
#endif

  explicit Impl(PmGridShape shape) : m_shape(shape) {}

  void shutdownBackendResources() {
#if COSMOSIM_ENABLE_FFTW
#if COSMOSIM_ENABLE_MPI
    bool has_distributed_plan = false;
    for (const auto& [_, plan] : m_plan_cache) {
      has_distributed_plan = has_distributed_plan ||
          (plan.is_distributed &&
           (plan.forward_plan != nullptr || plan.inverse_plan != nullptr));
    }
    if (has_distributed_plan) {
      int world_size = 1;
      int world_rank = 0;
      if (!queryActiveMpiWorld(world_size, world_rank)) {
        throw std::runtime_error(
            "PM FFTW-MPI plan shutdown requires an active MPI session");
      }
    }
#endif
    for (auto& [_, plan] : m_plan_cache) {
      if (plan.forward_plan != nullptr) {
        fftw_destroy_plan(plan.forward_plan);
        plan.forward_plan = nullptr;
      }
      if (plan.inverse_plan != nullptr) {
        fftw_destroy_plan(plan.inverse_plan);
        plan.inverse_plan = nullptr;
      }
    }
    if (m_isolated_workspace.forward_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.forward_plan);
      m_isolated_workspace.forward_plan = nullptr;
    }
    if (m_isolated_workspace.inverse_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.inverse_plan);
      m_isolated_workspace.inverse_plan = nullptr;
    }
#endif
    m_plan_cache.clear();
    m_active_key.reset();
  }

  void shutdownBackendResourcesNoexcept() noexcept {
#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
    bool has_distributed_plan = false;
    for (const auto& [_, plan] : m_plan_cache) {
      has_distributed_plan = has_distributed_plan ||
          (plan.is_distributed &&
           (plan.forward_plan != nullptr || plan.inverse_plan != nullptr));
    }
    if (has_distributed_plan) {
      int world_size = 1;
      int world_rank = 0;
      if (!queryActiveMpiWorld(world_size, world_rank)) {
        std::fprintf(stderr,
            "CHUI MPI LIFETIME ERROR: PM FFTW-MPI plans reached destructor outside an active MPI session; explicit backend shutdown ordering was violated\n");
        // The process is already outside legal MPI teardown. Do not call FFTW
        // destroy routines for MPI plans here; process teardown reclaims them.
        for (auto& [_, plan] : m_plan_cache) {
          if (plan.is_distributed) {
            plan.forward_plan = nullptr;
            plan.inverse_plan = nullptr;
          }
        }
      }
    }
#endif
    try {
      shutdownBackendResources();
    } catch (...) {
      // Destructors must stay noexcept. Any late distributed plan was surfaced
      // above; for other backend teardown failures process exit remains the
      // final containment boundary.
      m_plan_cache.clear();
      m_active_key.reset();
    }
  }

  ~Impl() { shutdownBackendResourcesNoexcept(); }

  struct IsolatedOpenWorkspace {
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;
    std::vector<std::complex<double>> rho_k;
    std::vector<std::complex<double>> kernel_k;
    std::vector<std::complex<double>> scratch;
    std::vector<double> potential_real;
    bool kernel_ready = false;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double split_scale = -1.0;
#if COSMOSIM_ENABLE_FFTW
    fftw_plan forward_plan = nullptr;
    fftw_plan inverse_plan = nullptr;
#endif
  };

  void ensureIsolatedWorkspace(std::size_t nx, std::size_t ny, std::size_t nz) {
    const bool shape_changed = nx != m_isolated_workspace.nx || ny != m_isolated_workspace.ny || nz != m_isolated_workspace.nz;
    if (!shape_changed) {
      return;
    }
#if COSMOSIM_ENABLE_FFTW
    if (m_isolated_workspace.forward_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.forward_plan);
      m_isolated_workspace.forward_plan = nullptr;
    }
    if (m_isolated_workspace.inverse_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.inverse_plan);
      m_isolated_workspace.inverse_plan = nullptr;
    }
#endif
    m_isolated_workspace.nx = nx;
    m_isolated_workspace.ny = ny;
    m_isolated_workspace.nz = nz;
    const std::size_t n_tot = nx * ny * nz;
    m_isolated_workspace.rho_k.assign(n_tot, {0.0, 0.0});
    m_isolated_workspace.kernel_k.assign(n_tot, {0.0, 0.0});
    m_isolated_workspace.scratch.assign(n_tot, {0.0, 0.0});
    m_isolated_workspace.potential_real.assign(n_tot, 0.0);
    m_isolated_workspace.kernel_ready = false;
    m_isolated_workspace.dx = 0.0;
    m_isolated_workspace.dy = 0.0;
    m_isolated_workspace.dz = 0.0;
    m_isolated_workspace.split_scale = -1.0;
#if COSMOSIM_ENABLE_FFTW
    configureFftwPlannerThreads();
    m_isolated_workspace.forward_plan = fftw_plan_dft_3d(
        static_cast<int>(nx),
        static_cast<int>(ny),
        static_cast<int>(nz),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE);
    m_isolated_workspace.inverse_plan = fftw_plan_dft_3d(
        static_cast<int>(nx),
        static_cast<int>(ny),
        static_cast<int>(nz),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE);
    if (m_isolated_workspace.forward_plan == nullptr || m_isolated_workspace.inverse_plan == nullptr) {
      throw std::runtime_error("failed to create isolated PM FFTW plans");
    }
#endif
  }

  void isolatedForward(std::span<std::complex<double>> field) {
    std::copy(field.begin(), field.end(), m_isolated_workspace.scratch.begin());
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(m_isolated_workspace.forward_plan);
    std::copy(m_isolated_workspace.scratch.begin(), m_isolated_workspace.scratch.end(), field.begin());
#else
    naiveComplexDft(field, /*forward=*/true);
#endif
  }

  void isolatedInverse(std::span<std::complex<double>> field) {
    std::copy(field.begin(), field.end(), m_isolated_workspace.scratch.begin());
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(m_isolated_workspace.inverse_plan);
    std::copy(m_isolated_workspace.scratch.begin(), m_isolated_workspace.scratch.end(), field.begin());
#else
    naiveComplexDft(field, /*forward=*/false);
#endif
  }

  [[nodiscard]] IsolatedOpenWorkspace& isolatedWorkspace() { return m_isolated_workspace; }

  [[nodiscard]] PlanResources& planForLayout(const parallel::PmSlabLayout& layout, core::PmDecompositionMode decomposition_mode) {
    const PlanKey key{
        .world_size = layout.world_size,
        .world_rank = layout.world_rank,
        .owned_begin_x = layout.owned_x.begin_x,
        .owned_end_x = layout.owned_x.end_x,
        .decomposition_mode = decomposition_mode,
    };
    auto it = m_plan_cache.find(key);
#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      const int local_cache_hit = it != m_plan_cache.end() ? 1 : 0;
      int global_cache_hit_count = 0;
      requirePmMpiSuccess(
          MPI_Allreduce(
              &local_cache_hit,
              &global_cache_hit_count,
              1,
              MPI_INT,
              MPI_SUM,
              MPI_COMM_WORLD),
          "PmSolver distributed FFT plan-cache agreement");
      if (global_cache_hit_count != 0 && global_cache_hit_count != layout.world_size) {
        throw std::runtime_error(
            "Distributed PM FFT plan cache state diverged across ranks; all ranks must build or reuse together");
      }
      if (global_cache_hit_count == layout.world_size) {
        m_active_key = key;
        return it->second;
      }
    } else
#endif
    if (it != m_plan_cache.end()) {
      m_active_key = key;
      return it->second;
    }

    PlanResources plan{};
    plan.layout = layout;
    const std::size_t nz_complex = m_shape.nz / 2U + 1U;
    PmPlanLocalStorageLayout storage_layout;
#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
    double plan_mpi_wait_ms = 0.0;
    if (layout.world_size > 1) {
      runPmCoordinatedPhase(
          "PmSolver distributed FFT storage-layout preparation",
          layout.world_rank,
          layout.world_size,
          plan_mpi_wait_ms,
          [&]() {
            storage_layout = determinePmPlanLocalStorageLayout(
                m_shape, layout, decomposition_mode, true);
          });
    } else
#endif
    {
      storage_layout = determinePmPlanLocalStorageLayout(
          m_shape, layout, decomposition_mode, false);
    }
    const std::size_t expected_local_complex_size =
        storage_layout.logical_local_complex_cells;
    const std::size_t allocated_local_complex_size =
        storage_layout.allocated_local_complex_cells;
    plan.real_z_stride = storage_layout.real_z_stride;
    plan.allocated_local_complex_cells = allocated_local_complex_size;

#if COSMOSIM_ENABLE_FFTW
#if COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      plan.is_distributed = storage_layout.is_distributed;
      plan.spectral_transposed = storage_layout.spectral_transposed;
      plan.transposed_local_ny = storage_layout.transposed_local_ny;
      plan.transposed_begin_y = storage_layout.transposed_begin_y;
      plan.fourier.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
    }
#endif
#if COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      configureFftwPlannerThreads();
      plan.forward_plan = fftw_mpi_plan_dft_r2c_3d(
          static_cast<ptrdiff_t>(m_shape.nx),
          static_cast<ptrdiff_t>(m_shape.ny),
          static_cast<ptrdiff_t>(m_shape.nz),
          plan.realView().data(),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          MPI_COMM_WORLD,
          decomposition_mode == core::PmDecompositionMode::kPencil
              ? (COSMOSIM_FFTW_PLANNER_FLAGS | FFTW_MPI_TRANSPOSED_OUT)
              : COSMOSIM_FFTW_PLANNER_FLAGS);
      plan.inverse_plan = fftw_mpi_plan_dft_c2r_3d(
          static_cast<ptrdiff_t>(m_shape.nx),
          static_cast<ptrdiff_t>(m_shape.ny),
          static_cast<ptrdiff_t>(m_shape.nz),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          plan.realView().data(),
          MPI_COMM_WORLD,
          decomposition_mode == core::PmDecompositionMode::kPencil
              ? (COSMOSIM_FFTW_PLANNER_FLAGS | FFTW_MPI_TRANSPOSED_IN)
              : COSMOSIM_FFTW_PLANNER_FLAGS);
      runPmCoordinatedPhase(
          "PmSolver distributed FFT plan creation",
          layout.world_rank,
          layout.world_size,
          plan_mpi_wait_ms,
          [&]() {
            if (plan.forward_plan == nullptr || plan.inverse_plan == nullptr) {
              throw std::runtime_error("Failed to create distributed FFTW plans for PM solver");
            }
          });
    } else
#endif
    {
      plan.real_z_stride = storage_layout.real_z_stride;
      plan.fourier.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
      configureFftwPlannerThreads();
      plan.forward_plan = fftw_plan_dft_r2c_3d(
          static_cast<int>(m_shape.nx),
          static_cast<int>(m_shape.ny),
          static_cast<int>(m_shape.nz),
          plan.realView().data(),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          COSMOSIM_FFTW_PLANNER_FLAGS);
      plan.inverse_plan = fftw_plan_dft_c2r_3d(
          static_cast<int>(m_shape.nx),
          static_cast<int>(m_shape.ny),
          static_cast<int>(m_shape.nz),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          plan.realView().data(),
          COSMOSIM_FFTW_PLANNER_FLAGS);
    }
    if (plan.forward_plan == nullptr || plan.inverse_plan == nullptr) {
      throw std::runtime_error("Failed to create FFTW plans for PM solver");
    }
#else
    if (!layout.ownsFullDomain()) {
      throw std::invalid_argument(
          "PM solver naive DFT fallback requires full-domain slab ownership; distributed PM requires COSMOSIM_ENABLE_FFTW=ON and COSMOSIM_ENABLE_MPI=ON");
    }
    plan.fourier.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
#endif

    plan.logical_local_complex_cells = plan.spectral_transposed
        ? checkedProduct(
              checkedProduct(
                  plan.transposed_local_ny, m_shape.nx,
                  "PM transposed logical spectral extent"),
              nz_complex,
              "PM transposed logical spectral extent")
        : expected_local_complex_size;
    if (plan.logical_local_complex_cells > plan.allocated_local_complex_cells) {
      throw std::runtime_error(
          "PM logical spectral extent exceeds the backend allocation");
    }
    // In-place layout invariant (every backend, including the naive DFT): the
    // aliased real view must cover every physical row of this rank's x slab
    // at the padded stride 2 (nz/2 + 1), so realRowBase() + nz never leaves
    // the complex allocation.
    if (plan.real_z_stride != checkedProduct(2U, nz_complex, "PM padded real stride") ||
        checkedProduct(checkedProduct(layout.local_nx(), m_shape.ny, "PM padded real coverage"),
                       plan.real_z_stride, "PM padded real coverage") >
            checkedProduct(2U, plan.fourier.size(), "PM padded real coverage")) {
      throw std::logic_error(
          "PM in-place real view does not cover the owned slab at the padded r2c stride");
    }

#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      PlanResources* committed_plan = nullptr;
      runPmCoordinatedPhase(
          "PmSolver distributed FFT plan cache commit",
          layout.world_rank,
          layout.world_size,
          plan_mpi_wait_ms,
          [&]() {
            auto [insert_it, inserted] = m_plan_cache.emplace(key, std::move(plan));
            if (!inserted) {
              throw std::logic_error("Distributed PM FFT plan cache changed during coordinated build");
            }
            ++m_plan_build_count;
            m_active_key = key;
            committed_plan = &insert_it->second;
          });
      return *committed_plan;
    }
#endif
    auto [insert_it, inserted] = m_plan_cache.emplace(key, std::move(plan));
    (void)inserted;
    ++m_plan_build_count;
    m_active_key = key;
    return insert_it->second;
  }

  [[nodiscard]] std::span<double> realGrid() { return activePlan().realView(); }
  [[nodiscard]] std::span<std::complex<double>> fourierGrid() { return activePlan().fourier; }
  // Spectral-gradient phi_k copy, sized to the backend Fourier allocation on
  // first use and retained for later spectral solves.
  [[nodiscard]] std::span<std::complex<double>> potentialScratch() {
    auto& plan = activePlan();
    if (plan.potential_k.size() != plan.fourier.size()) {
      plan.potential_k.assign(plan.fourier.size(), std::complex<double>(0.0, 0.0));
    }
    return plan.potential_k;
  }
  [[nodiscard]] std::span<double> finiteDifferenceHaloPlanes(std::size_t plane_values) {
    auto& plan = activePlan();
    // [left -2, left -1, right +1, right +2, send staging x2]
    const std::size_t required = checkedProduct(6U, plane_values, "PM finite-difference halo planes");
    if (plan.fd_halo_planes.size() != required) plan.fd_halo_planes.assign(required, 0.0);
    return plan.fd_halo_planes;
  }

  void invalidatePeriodicDensityForGrid(const PmGridStorage* grid) noexcept {
    for (auto& [_, plan] : m_plan_cache) {
      if (plan.periodic_density_grid == grid) {
        plan.periodic_density_ready = false;
        plan.periodic_density_grid = nullptr;
        plan.periodic_density_generation = 0U;
        plan.periodic_density_grid_generation = 0U;
      }
    }
  }

  void publishPeriodicDensity(
      PlanResources& plan,
      const PmGridStorage& grid) {
    if (m_periodic_density_generation ==
        std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error(
          "PM periodic density generation would overflow; recreate the solver");
    }
    ++m_periodic_density_generation;
    plan.periodic_density_grid = &grid;
    plan.periodic_density_generation = m_periodic_density_generation;
    plan.periodic_density_grid_generation = grid.storageGeneration();
    plan.periodic_density_ready = true;
  }

  [[nodiscard]] bool periodicDensityReadyFor(
      const PlanResources& plan,
      const PmGridStorage& grid) const noexcept {
    return plan.periodic_density_ready &&
        plan.periodic_density_grid == &grid &&
        plan.periodic_density_generation != 0U &&
        plan.periodic_density_grid_generation == grid.storageGeneration();
  }

  [[nodiscard]] bool periodicDensityTokenTargetsGrid(
      const PlanResources& plan,
      const PmGridStorage& grid) const noexcept {
    return plan.periodic_density_ready &&
        plan.periodic_density_grid == &grid &&
        plan.periodic_density_generation != 0U;
  }

  void consumePeriodicDensity(PlanResources& plan) noexcept {
    plan.periodic_density_ready = false;
    plan.periodic_density_grid = nullptr;
    plan.periodic_density_generation = 0U;
    plan.periodic_density_grid_generation = 0U;
  }

  [[nodiscard]] bool ensureSpectralAxisMetadata(
      PlanResources& plan,
      const BoxLengths& lengths,
      const PmSolveOptions& options,
      const PmGridShape& shape) {
    if (plan.spectral_axis_metadata_ready &&
        plan.cached_lx == lengths.lx &&
        plan.cached_ly == lengths.ly &&
        plan.cached_lz == lengths.lz &&
        plan.cached_window_deconvolution == options.enable_window_deconvolution &&
        plan.cached_assignment_scheme == options.assignment_scheme) {
      return false;
    }

    const std::size_t nz_complex = shape.nz / 2U + 1U;
    plan.kx_by_ix.resize(shape.nx);
    plan.ky_by_iy.resize(shape.ny);
    plan.kz_by_iz.resize(nz_complex);
    plan.window_x_by_ix.assign(shape.nx, 1.0);
    plan.window_y_by_iy.assign(shape.ny, 1.0);
    plan.window_z_by_iz.assign(nz_complex, 1.0);

    const double dkx = 2.0 * k_pi / lengths.lx;
    const double dky = 2.0 * k_pi / lengths.ly;
    const double dkz = 2.0 * k_pi / lengths.lz;
    const int window_exponent = assignmentWindowExponent(options.assignment_scheme);

    for (std::size_t ix = 0; ix < shape.nx; ++ix) {
      const std::ptrdiff_t mode = ix <= shape.nx / 2U
          ? static_cast<std::ptrdiff_t>(ix)
          : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(shape.nx);
      const double kx = dkx * static_cast<double>(mode);
      plan.kx_by_ix[ix] = kx;
      if (options.enable_window_deconvolution) {
        plan.window_x_by_ix[ix] = std::pow(
            sinc(0.5 * kx * lengths.lx / static_cast<double>(shape.nx)),
            static_cast<double>(window_exponent));
      }
    }
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      const std::ptrdiff_t mode = iy <= shape.ny / 2U
          ? static_cast<std::ptrdiff_t>(iy)
          : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(shape.ny);
      const double ky = dky * static_cast<double>(mode);
      plan.ky_by_iy[iy] = ky;
      if (options.enable_window_deconvolution) {
        plan.window_y_by_iy[iy] = std::pow(
            sinc(0.5 * ky * lengths.ly / static_cast<double>(shape.ny)),
            static_cast<double>(window_exponent));
      }
    }
    for (std::size_t iz = 0; iz < nz_complex; ++iz) {
      const double kz = dkz * static_cast<double>(iz);
      plan.kz_by_iz[iz] = kz;
      if (options.enable_window_deconvolution) {
        plan.window_z_by_iz[iz] = std::pow(
            sinc(0.5 * kz * lengths.lz / static_cast<double>(shape.nz)),
            static_cast<double>(window_exponent));
      }
    }

    plan.spectral_axis_metadata_ready = true;
    plan.cached_lx = lengths.lx;
    plan.cached_ly = lengths.ly;
    plan.cached_lz = lengths.lz;
    plan.cached_window_deconvolution = options.enable_window_deconvolution;
    plan.cached_assignment_scheme = options.assignment_scheme;
    return true;
  }

  template <typename Callback>
  void forEachLogicalSpectralMode(
      const PlanResources& plan,
      const PmGridShape& shape,
      Callback&& callback) const {
    const std::size_t nz_complex = shape.nz / 2U + 1U;
    std::size_t visited = 0U;
    if (plan.spectral_transposed) {
#if COSMOSIM_HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+ : visited)
#endif
      for (std::size_t local_iy = 0; local_iy < plan.transposed_local_ny; ++local_iy) {
        const std::size_t iy = plan.transposed_begin_y + local_iy;
        for (std::size_t ix = 0; ix < shape.nx; ++ix) {
          for (std::size_t iz = 0; iz < nz_complex; ++iz) {
            const std::size_t index = (local_iy * shape.nx + ix) * nz_complex + iz;
            callback(
                index, ix, iy, iz,
                plan.kx_by_ix[ix], plan.ky_by_iy[iy], plan.kz_by_iz[iz]);
            ++visited;
          }
        }
      }
    } else {
      const std::size_t global_x_begin = plan.layout.owned_x.begin_x;
#if COSMOSIM_HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+ : visited)
#endif
      for (std::size_t local_ix = 0; local_ix < plan.layout.local_nx(); ++local_ix) {
        const std::size_t ix = global_x_begin + local_ix;
        for (std::size_t iy = 0; iy < shape.ny; ++iy) {
          for (std::size_t iz = 0; iz < nz_complex; ++iz) {
            const std::size_t index = (local_ix * shape.ny + iy) * nz_complex + iz;
            callback(
                index, ix, iy, iz,
                plan.kx_by_ix[ix], plan.ky_by_iy[iy], plan.kz_by_iz[iz]);
            ++visited;
          }
        }
      }
    }
    if (visited != plan.logical_local_complex_cells) {
      throw std::logic_error(
          "PM logical spectral traversal disagrees with the plan's logical extent");
    }
  }

  // Deterministic serial variant (N03 observation: storage-order sums).
  template <typename Callback>
  void forEachLogicalSpectralModeSerial(
      const PlanResources& plan,
      const PmGridShape& shape,
      Callback&& callback) const {
    const std::size_t nz_complex = shape.nz / 2U + 1U;
    std::size_t visited = 0U;
    if (plan.spectral_transposed) {
      for (std::size_t local_iy = 0; local_iy < plan.transposed_local_ny; ++local_iy) {
        const std::size_t iy = plan.transposed_begin_y + local_iy;
        for (std::size_t ix = 0; ix < shape.nx; ++ix) {
          for (std::size_t iz = 0; iz < nz_complex; ++iz) {
            const std::size_t index = (local_iy * shape.nx + ix) * nz_complex + iz;
            callback(
                index, ix, iy, iz,
                plan.kx_by_ix[ix], plan.ky_by_iy[iy], plan.kz_by_iz[iz]);
            ++visited;
          }
        }
      }
    } else {
      const std::size_t global_x_begin = plan.layout.owned_x.begin_x;
      for (std::size_t local_ix = 0; local_ix < plan.layout.local_nx(); ++local_ix) {
        const std::size_t ix = global_x_begin + local_ix;
        for (std::size_t iy = 0; iy < shape.ny; ++iy) {
          for (std::size_t iz = 0; iz < nz_complex; ++iz) {
            const std::size_t index = (local_ix * shape.ny + iy) * nz_complex + iz;
            callback(
                index, ix, iy, iz,
                plan.kx_by_ix[ix], plan.ky_by_iy[iy], plan.kz_by_iz[iz]);
            ++visited;
          }
        }
      }
    }
    if (visited != plan.logical_local_complex_cells) {
      throw std::logic_error(
          "PM logical spectral traversal disagrees with the plan's logical extent");
    }
  }

  [[nodiscard]] bool applyDerivedPoissonOperator(
      PlanResources& plan,
      const BoxLengths& lengths,
      const PmSolveOptions& options,
      const PmGridShape& shape) {
    const bool rebuilt = ensureSpectralAxisMetadata(plan, lengths, options, shape);
    auto fourier = fourierGrid();
    if (fourier.size() != plan.allocated_local_complex_cells ||
        plan.logical_local_complex_cells > fourier.size()) {
      throw std::logic_error("PM Fourier storage extent disagrees with plan metadata");
    }

    const double prefactor = -4.0 * k_pi * options.gravitational_constant_code;
    // N03: optional isotropic source hook. Observation (binned mode power of
    // the deposited cold density) precedes any modification of the modes; the
    // collective finalize runs on every rank, including empty spectral slabs.
    PmSpectralSourceHook* const hook = options.spectral_source_hook;
    if (hook != nullptr) {
      double k_axis_min = std::numeric_limits<double>::infinity();
      double kx_max = 0.0;
      double ky_max = 0.0;
      double kz_max = 0.0;
      const auto scan = [&k_axis_min](const auto& axis, double& axis_max) {
        for (const double value : axis) {
          const double magnitude = std::abs(value);
          axis_max = std::max(axis_max, magnitude);
          if (magnitude > 0.0) k_axis_min = std::min(k_axis_min, magnitude);
        }
      };
      scan(plan.kx_by_ix, kx_max);
      scan(plan.ky_by_iy, ky_max);
      scan(plan.kz_by_iz, kz_max);
      hook->requireWavenumberDomain(
          k_axis_min, std::sqrt(kx_max * kx_max + ky_max * ky_max + kz_max * kz_max));
      hook->beginSolve(options.scale_factor);
      std::vector<double> bin_power;
      std::vector<double> bin_modes;
      if (hook->needsObservation()) {
        bin_power.assign(hook->binCount(), 0.0);
        bin_modes.assign(hook->binCount(), 0.0);
        const std::size_t nz_complex = shape.nz / 2U + 1U;
        const bool nz_even = (shape.nz % 2U) == 0U;
        forEachLogicalSpectralModeSerial(
            plan, shape,
            [&](std::size_t index, std::size_t, std::size_t, std::size_t iz, double kx, double ky, double kz) {
              const double k2 = kx * kx + ky * ky + kz * kz;
              if (k2 == 0.0) return;
              const std::size_t bin = hook->binIndex(std::sqrt(k2));
              if (bin == PmSpectralSourceHook::k_no_bin) return;
              const double multiplicity = (iz == 0U || (nz_even && iz == nz_complex - 1U)) ? 1.0 : 2.0;
              bin_power[bin] += multiplicity * std::norm(fourier[index]);
              bin_modes[bin] += multiplicity;
            });
      }
      hook->finalizeObservation(bin_power, bin_modes);
    }
    forEachLogicalSpectralMode(
        plan, shape,
        [&](std::size_t index,
            std::size_t ix,
            std::size_t iy,
            std::size_t iz,
            double kx,
            double ky,
            double kz) {
          const double k2 = kx * kx + ky * ky + kz * kz;
          if (k2 == 0.0) {
            fourier[index] = std::complex<double>(0.0, 0.0);
            return;
          }
          double window_correction = 1.0;
          if (options.enable_window_deconvolution) {
            const double transfer_window =
                plan.window_x_by_ix[ix] *
                plan.window_y_by_iy[iy] *
                plan.window_z_by_iz[iz];
            window_correction =
                1.0 / std::max(transfer_window * transfer_window, 1.0e-12);
          }
          double split_filter = 1.0;
          if (options.tree_pm_split_scale_comoving > 0.0) {
            split_filter = treePmGaussianFourierLongRangeFilterUnchecked(
                std::sqrt(k2), options.tree_pm_split_scale_comoving);
          }
          const double source_multiplier =
              hook == nullptr ? 1.0 : hook->multiplier(std::sqrt(k2));
          const double poisson_operator =
              prefactor * window_correction * split_filter * source_multiplier / k2;
          fourier[index] *= poisson_operator;
        });

    // FFTW MPI may over-allocate beyond the logical spectral extent. The old
    // zero-initialized full scalar operator implicitly zeroed that tail;
    // preserve the same fail-closed behavior explicitly.
    std::fill(
        fourier.begin() + static_cast<std::ptrdiff_t>(plan.logical_local_complex_cells),
        fourier.end(),
        std::complex<double>(0.0, 0.0));
    return rebuilt;
  }

  void fillGradientSpectrum(
      PlanResources& plan,
      const BoxLengths& lengths,
      const PmGridShape& shape,
      int axis) {
    if (axis < 0 || axis > 2) {
      throw std::invalid_argument("PM spectral gradient axis must be 0, 1, or 2");
    }
    if (!plan.spectral_axis_metadata_ready ||
        plan.cached_lx != lengths.lx ||
        plan.cached_ly != lengths.ly ||
        plan.cached_lz != lengths.lz) {
      throw std::logic_error(
          "PM gradient requested without current spectral axis metadata");
    }

    auto dst = fourierGrid();
    const auto potential_k = std::span<const std::complex<double>>(
        plan.potential_k.data(), plan.potential_k.size());
    std::fill(dst.begin(), dst.end(), std::complex<double>(0.0, 0.0));

    const auto first_derivative_mode = [](std::size_t mode_index,
                                          std::size_t mode_count,
                                          double wave_number) {
      // Preserve the existing real-DFT first-derivative convention: only the
      // Nyquist coefficient of the differentiated even-length axis is zero.
      return mode_count % 2U == 0U && mode_index == mode_count / 2U
          ? 0.0
          : wave_number;
    };

    forEachLogicalSpectralMode(
        plan, shape,
        [&](std::size_t index,
            std::size_t ix,
            std::size_t iy,
            std::size_t iz,
            double kx,
            double ky,
            double kz) {
          const double gradient = axis == 0
              ? first_derivative_mode(ix, shape.nx, kx)
              : (axis == 1
                  ? first_derivative_mode(iy, shape.ny, ky)
                  : first_derivative_mode(iz, shape.nz, kz));
          dst[index] =
              std::complex<double>(0.0, -gradient) * potential_k[index];
        });
  }

  double forwardFft() {
    const auto start = std::chrono::steady_clock::now();
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(activePlan().forward_plan);
#else
    naiveForwardDft();
#endif
    const auto stop = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
  }

  double inverseFft() {
    const auto start = std::chrono::steady_clock::now();
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(activePlan().inverse_plan);
#else
    naiveInverseDft();
#endif
    const auto stop = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
  }

  [[nodiscard]] std::size_t planCount() const noexcept { return m_plan_cache.size(); }
  [[nodiscard]] std::size_t planBuildCount() const noexcept { return m_plan_build_count; }

  [[nodiscard]] std::uint64_t nextDistributedExchangeEpoch() {
    if (m_next_distributed_exchange_epoch == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error(
          "PmSolver distributed exchange epoch would overflow; recreate the solver before another exchange");
    }
    const std::uint64_t epoch = m_next_distributed_exchange_epoch;
    ++m_next_distributed_exchange_epoch;
    return epoch;
  }

  // required_bytes is the workload-derived phase requirement. An external
  // (coordinator-owned) arena must already hold it — the coordinator grows it
  // at the idle boundary before leasing; the standalone fallback grows itself
  // here (it is idle between calls) instead of retaining the 128 MiB ceiling.
  [[nodiscard]] std::pmr::memory_resource* communicationResource(
      GravityCommunicationArena::Phase phase,
      std::optional<GravityCommunicationArena::Lease>& fallback_lease,
      std::uint64_t required_bytes = k_pm_routing_workspace_target_bytes) {
    if (m_external_communication_arena != nullptr) {
      if (m_external_communication_arena->activePhase() != phase) {
        throw std::logic_error("PmSolver external gravity communication arena phase mismatch");
      }
      if (m_external_communication_arena->capacityBytes() < required_bytes) {
        throw std::logic_error(
            "PmSolver external gravity communication arena is smaller than the planned routing workload; "
            "the coordinator must grow it before leasing");
      }
      return m_external_communication_resource;
    }
    if (m_standalone_communication_arena == nullptr) {
      m_standalone_communication_arena = std::make_unique<GravityCommunicationArena>(
          nullptr, "pm_solver.communication_arena");
    }
    m_standalone_communication_arena->configure(std::max<std::uint64_t>(required_bytes, 1U));
    fallback_lease.emplace(m_standalone_communication_arena->begin(phase));
    return fallback_lease->resource();
  }

  // One collective workload plan per routing phase, produced by
  // PmSolver::plan*Routing() and withdrawn by the next routing call of the
  // same kind (success or failure). Lifetime is operation-scoped: the plan is
  // only a certificate that the communication arena was configured large
  // enough; the routing call re-measures its own workload and requires the
  // plan to bound it (key equality plus count / per-peer / round maxima),
  // agreed collectively, so a stale plan can never undersize a reservation.
  struct RoutingWorkloadPlanCache {
    bool valid = false;
    PmRoutingPlanKey key{};
    std::uint64_t global_max_count = 0U;
    std::uint64_t global_max_peer_records = 0U;
    std::uint64_t local_max_round_records = 0U;
  };
  RoutingWorkloadPlanCache m_density_routing_plan;
  RoutingWorkloadPlanCache m_force_routing_plan;

  void attachCommunicationArena(GravityCommunicationArena* arena) {
    if (arena == nullptr || arena->activePhase() == GravityCommunicationArena::Phase::kIdle) {
      throw std::invalid_argument("PmSolver requires an active external communication arena lease");
    }
    // External production ownership and standalone fallback ownership are
    // mutually exclusive.  Do not retain a full fallback after attachment.
    m_standalone_communication_arena.reset();
    m_external_communication_arena = arena;
    // The resource is phase-valid only while the coordinator's lease lives.
    switch (arena->activePhase()) {
      case GravityCommunicationArena::Phase::kPmDensity:
      case GravityCommunicationArena::Phase::kPmInterpolation:
        break;
      default:
        throw std::invalid_argument("PmSolver external arena phase is not a PM routing phase");
    }
    m_external_communication_resource = arena->activeResource();
  }

  void detachCommunicationArena() noexcept {
    m_external_communication_resource = nullptr;
    m_external_communication_arena = nullptr;
  }

  void recordCommunicationHighWater(
      GravityCommunicationArena::Phase phase,
      std::uint64_t bytes,
      std::optional<GravityCommunicationArena::Lease>& fallback_lease) {
    if (phase == GravityCommunicationArena::Phase::kPmDensity) {
      m_density_exchange_logical_high_water_bytes =
          std::max(m_density_exchange_logical_high_water_bytes, bytes);
    } else if (phase == GravityCommunicationArena::Phase::kPmInterpolation) {
      m_plane_interpolation_logical_high_water_bytes =
          std::max(m_plane_interpolation_logical_high_water_bytes, bytes);
    } else {
      throw std::logic_error("PmSolver attempted to record a non-PM communication phase");
    }
    if (m_external_communication_arena != nullptr) {
      m_external_communication_arena->recordActiveLogicalHighWater(bytes);
    } else if (fallback_lease.has_value()) {
      fallback_lease->recordLogicalHighWater(bytes);
    }
  }

  [[nodiscard]] IndexedTargetWorkspace& indexedTargetWorkspace() noexcept {
    return m_indexed_target_workspace;
  }
#if COSMOSIM_ENABLE_CUDA
  [[nodiscard]] CudaBackendWorkspace& cudaWorkspace() noexcept {
    return m_cuda_workspace;
  }
#endif
  void appendOwnedMemoryReport(core::MemoryReportBuilder& builder) const {
    appendMemoryReport(builder);
  }

#if COSMOSIM_ENABLE_CUDA
  // Current CUDA staging strategy. All transfer/deposit/host-field-solve/
  // interpolation choreography lives behind Impl so a future resident strategy
  // can replace these stages without changing solveForParticles().
  void solveCudaStagingBackend(
      PmSolver& solver,
      PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<const double> mass,
      std::span<double> accel_x,
      std::span<double> accel_y,
      std::span<double> accel_z,
      const PmSolveOptions& options,
      PmProfileEvent* profile) {
    auto& cuda_workspace = m_cuda_workspace;
    const core::CudaStreamView stream_view = cuda_workspace.prepare(
        pos_x.size(), m_shape.cellCount());
    const cudaStream_t stream = cuda_workspace.stream;
    auto& pos_x_device = cuda_workspace.pos_x;
    auto& pos_y_device = cuda_workspace.pos_y;
    auto& pos_z_device = cuda_workspace.pos_z;
    auto& mass_device = cuda_workspace.mass;
    auto& density_device = cuda_workspace.density;
    auto& force_x_device = cuda_workspace.force_x;
    auto& force_y_device = cuda_workspace.force_y;
    auto& force_z_device = cuda_workspace.force_z;
    auto& accel_x_device = cuda_workspace.accel_x;
    auto& accel_y_device = cuda_workspace.accel_y;
    auto& accel_z_device = cuda_workspace.accel_z;

    {
      const auto copy_h2d_start = std::chrono::steady_clock::now();
      pos_x_device.copyFromHostAsync(pos_x, stream_view);
      pos_y_device.copyFromHostAsync(pos_y, stream_view);
      pos_z_device.copyFromHostAsync(pos_z, stream_view);
      mass_device.copyFromHostAsync(mass, stream_view);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing H2D particle copy");
      }
      const auto copy_h2d_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_h2d_ms += std::chrono::duration<double, std::milli>(copy_h2d_stop - copy_h2d_start).count();
      }

      const BoxLengths lengths = effectiveBoxLengths(options);

      const auto kernel_assign_start = std::chrono::steady_clock::now();
      pmCudaAssignDensityCic(
          PmCudaAssignLaunch{
              pos_x.size(),
              m_shape.nx,
              m_shape.ny,
              m_shape.nz,
              lengths.lx,
              lengths.ly,
              lengths.lz},
          pos_x_device.data(),
          pos_y_device.data(),
          pos_z_device.data(),
          mass_device.data(),
          density_device.data(),
          stream_view);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing PM assignment kernel");
      }
      const auto kernel_assign_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->device_kernel_ms +=
            std::chrono::duration<double, std::milli>(kernel_assign_stop - kernel_assign_start).count();
      }

      const auto copy_density_start = std::chrono::steady_clock::now();
      density_device.copyToHostAsync(grid.density(), stream_view);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing D2H density copy");
      }
      const auto copy_density_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_d2h_ms += std::chrono::duration<double, std::milli>(copy_density_stop - copy_density_start).count();
      }
      const double cell_volume =
          (lengths.lx * lengths.ly * lengths.lz) / static_cast<double>(m_shape.cellCount());
      for (double& density_cell : grid.density()) {
        density_cell /= cell_volume;
      }

      solver.solvePoissonPeriodic(grid, options, profile);

      const auto copy_forces_start = std::chrono::steady_clock::now();
      force_x_device.copyFromHostAsync(grid.force_x(), stream_view);
      force_y_device.copyFromHostAsync(grid.force_y(), stream_view);
      force_z_device.copyFromHostAsync(grid.force_z(), stream_view);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing H2D force copy");
      }
      const auto copy_forces_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_h2d_ms += std::chrono::duration<double, std::milli>(copy_forces_stop - copy_forces_start).count();
      }

      const auto kernel_interp_start = std::chrono::steady_clock::now();
      pmCudaInterpolateForcesCic(
          PmCudaInterpLaunch{
              pos_x.size(),
              m_shape.nx,
              m_shape.ny,
              m_shape.nz,
              lengths.lx,
              lengths.ly,
              lengths.lz},
          pos_x_device.data(),
          pos_y_device.data(),
          pos_z_device.data(),
          force_x_device.data(),
          force_y_device.data(),
          force_z_device.data(),
          accel_x_device.data(),
          accel_y_device.data(),
          accel_z_device.data(),
          stream_view);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing PM interpolation kernel");
      }
      const auto kernel_interp_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->device_kernel_ms +=
            std::chrono::duration<double, std::milli>(kernel_interp_stop - kernel_interp_start).count();
      }

      const auto copy_accel_start = std::chrono::steady_clock::now();
      accel_x_device.copyToHostAsync(accel_x, stream_view);
      accel_y_device.copyToHostAsync(accel_y, stream_view);
      accel_z_device.copyToHostAsync(accel_z, stream_view);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing D2H acceleration copy");
      }
      const auto copy_accel_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_d2h_ms += std::chrono::duration<double, std::milli>(copy_accel_stop - copy_accel_start).count();
      }
    }

    if (profile != nullptr) {
      profile->bytes_moved += bytesForParticles(pos_x.size()) * 2U;
      profile->bytes_moved += bytesForGridSweep(m_shape.cellCount()) * 4U;
    }
    return;
  }
#endif

 private:
  [[nodiscard]] PlanResources& activePlan() {
    if (!m_active_key.has_value()) {
      throw std::logic_error("PM solver plan has not been initialized for the active slab layout");
    }
    return m_plan_cache.at(*m_active_key);
  }

#if !COSMOSIM_ENABLE_FFTW
  void naiveComplexDft(std::span<std::complex<double>> field, bool forward) {
    const std::size_t nx = m_isolated_workspace.nx;
    const std::size_t ny = m_isolated_workspace.ny;
    const std::size_t nz = m_isolated_workspace.nz;
    const std::size_t total = nx * ny * nz;
    std::vector<std::complex<double>> out(total, {0.0, 0.0});
    const double sign = forward ? -1.0 : 1.0;
    for (std::size_t kx = 0; kx < nx; ++kx) {
      for (std::size_t ky = 0; ky < ny; ++ky) {
        for (std::size_t kz = 0; kz < nz; ++kz) {
          std::complex<double> acc(0.0, 0.0);
          for (std::size_t x = 0; x < nx; ++x) {
            for (std::size_t y = 0; y < ny; ++y) {
              for (std::size_t z = 0; z < nz; ++z) {
                const double phase = sign * 2.0 * k_pi *
                    (static_cast<double>(kx * x) / static_cast<double>(nx) +
                     static_cast<double>(ky * y) / static_cast<double>(ny) +
                     static_cast<double>(kz * z) / static_cast<double>(nz));
                const std::complex<double> euler(std::cos(phase), std::sin(phase));
                acc += field[(x * ny + y) * nz + z] * euler;
              }
            }
          }
          out[(kx * ny + ky) * nz + kz] = acc;
        }
      }
    }
    if (!forward) {
      const double inv_total = 1.0 / static_cast<double>(total);
      for (auto& v : out) {
        v *= inv_total;
      }
    }
    std::copy(out.begin(), out.end(), field.begin());
  }

  // Diagnostic naive DFTs. The FFT working set is in place (the real view
  // aliases `fourier`), so each transform first stages its input in a
  // transient copy; this fallback exists only for small non-FFTW test meshes.
  void naiveForwardDft() {
    if (!activePlan().layout.ownsFullDomain()) {
      throw std::logic_error(
          "PM solver naiveForwardDft requires full-domain ownership; distributed PM is unavailable without FFTW/MPI");
    }
    auto& plan = activePlan();
    const std::size_t stride = plan.real_z_stride;
    std::vector<double> real(m_shape.cellCount());
    {
      const auto view = plan.realView();
      for (std::size_t row = 0; row < m_shape.nx * m_shape.ny; ++row) {
        std::copy_n(view.begin() + static_cast<std::ptrdiff_t>(row * stride), m_shape.nz,
                    real.begin() + static_cast<std::ptrdiff_t>(row * m_shape.nz));
      }
    }
    auto& fourier = plan.fourier;
    const std::size_t nz_complex = m_shape.nz / 2U + 1U;
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < nz_complex; ++iz) {
          std::complex<double> acc(0.0, 0.0);
          for (std::size_t x = 0; x < m_shape.nx; ++x) {
            for (std::size_t y = 0; y < m_shape.ny; ++y) {
              for (std::size_t z = 0; z < m_shape.nz; ++z) {
                const double phase = -2.0 * k_pi *
                    (static_cast<double>(ix * x) / static_cast<double>(m_shape.nx) +
                     static_cast<double>(iy * y) / static_cast<double>(m_shape.ny) +
                     static_cast<double>(iz * z) / static_cast<double>(m_shape.nz));
                const std::complex<double> euler(std::cos(phase), std::sin(phase));
                const std::size_t rindex = (x * m_shape.ny + y) * m_shape.nz + z;
                acc += real[rindex] * euler;
              }
            }
          }
          const std::size_t cindex = (ix * m_shape.ny + iy) * nz_complex + iz;
          fourier[cindex] = acc;
        }
      }
    }
  }

  void naiveInverseDft() {
    if (!activePlan().layout.ownsFullDomain()) {
      throw std::logic_error(
          "PM solver naiveInverseDft requires full-domain ownership; distributed PM is unavailable without FFTW/MPI");
    }
    auto& plan = activePlan();
    const std::vector<std::complex<double>> fourier(plan.fourier.begin(), plan.fourier.end());
    auto real = plan.realView();
    const std::size_t stride = plan.real_z_stride;
    const std::size_t total = m_shape.cellCount();
    std::fill(real.begin(), real.end(), 0.0);
    const std::size_t nz_complex = m_shape.nz / 2U + 1U;

    for (std::size_t x = 0; x < m_shape.nx; ++x) {
      for (std::size_t y = 0; y < m_shape.ny; ++y) {
        for (std::size_t z = 0; z < m_shape.nz; ++z) {
          std::complex<double> acc(0.0, 0.0);
          for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
            for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
              for (std::size_t iz = 0; iz < nz_complex; ++iz) {
                const double phase = 2.0 * k_pi *
                    (static_cast<double>(ix * x) / static_cast<double>(m_shape.nx) +
                     static_cast<double>(iy * y) / static_cast<double>(m_shape.ny) +
                     static_cast<double>(iz * z) / static_cast<double>(m_shape.nz));
                const std::complex<double> euler(std::cos(phase), std::sin(phase));
                const std::size_t cindex = (ix * m_shape.ny + iy) * nz_complex + iz;
                if (iz == 0 || (m_shape.nz % 2U == 0U && iz == m_shape.nz / 2U)) {
                  acc += fourier[cindex] * euler;
                } else {
                  acc += fourier[cindex] * euler + std::conj(fourier[cindex]) * std::conj(euler);
                }
              }
            }
          }
          real[(x * m_shape.ny + y) * stride + z] = acc.real() / static_cast<double>(total);
        }
      }
    }
  }
#endif

  void appendMemoryReport(core::MemoryReportBuilder& builder) const {
    const auto checked_add = [](std::uint64_t lhs, std::uint64_t rhs) {
      if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
        throw std::overflow_error("PM memory report byte sum overflow");
      }
      return lhs + rhs;
    };
    const auto accumulate = [&checked_add](
                                std::uint64_t& current,
                                std::uint64_t& capacity,
                                const auto& container) {
      current = checked_add(current, core::currentSizeBytesForContainer(container));
      capacity = checked_add(capacity, core::ownedCapacityBytesForContainer(container));
    };
    const auto emit = [&builder](
                          core::MemorySubsystem subsystem,
                          core::MemoryLifetime lifetime,
                          std::string label,
                          std::uint64_t current,
                          std::uint64_t capacity,
                          std::uint64_t high_water,
                          std::string note) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = subsystem,
          .lifetime = lifetime,
          .label = std::move(label),
          .current_size_bytes = current,
          .owned_capacity_bytes = capacity,
          .high_water_bytes = std::max(capacity, high_water),
          .estimated_next_step_bytes = 0U,
          .uncertainty_note = std::move(note),
      });
    };

    {
      std::uint64_t current = 0U;
      std::uint64_t capacity = 0U;
      accumulate(current, capacity, m_indexed_target_workspace.compact_ax);
      accumulate(current, capacity, m_indexed_target_workspace.compact_ay);
      accumulate(current, capacity, m_indexed_target_workspace.compact_az);
      emit(core::MemorySubsystem::kScratch, core::MemoryLifetime::kTransient,
           "pm_solver.indexed_target_workspace", current, capacity, capacity,
           "retained indexed-global output scatter scratch only; indexed coordinates borrow canonical source lanes directly");
    }

    {
      std::uint64_t current = 0U;
      std::uint64_t capacity = 0U;
      for (const auto& [_, resources] : m_plan_cache) {
        accumulate(current, capacity, resources.fourier);
        accumulate(current, capacity, resources.potential_k);
        accumulate(current, capacity, resources.fd_halo_planes);
      }
      emit(core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kPersistent,
           "pm_solver.plan_cache_owned_arrays", current, capacity, capacity,
           "owned FFT real + Fourier working + preserved potential spectrum (spectral gradient backend only) + finite-difference x halo planes (finite-difference backend, distributed only); no full Poisson-kernel lane");
    }

    {
      std::uint64_t current = 0U;
      std::uint64_t capacity = 0U;
      for (const auto& [_, resources] : m_plan_cache) {
        accumulate(current, capacity, resources.kx_by_ix);
        accumulate(current, capacity, resources.ky_by_iy);
        accumulate(current, capacity, resources.kz_by_iz);
        accumulate(current, capacity, resources.window_x_by_ix);
        accumulate(current, capacity, resources.window_y_by_iy);
        accumulate(current, capacity, resources.window_z_by_iz);
      }
      emit(core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kPersistent,
           "pm_solver.plan_spectral_axis_metadata", current, capacity, capacity,
           "O(nx+ny+nz) wave-number and assignment-window tables used to derive the Poisson operator");
    }

    {
      std::uint64_t current = 0U;
      std::uint64_t capacity = 0U;
      accumulate(current, capacity, m_isolated_workspace.rho_k);
      accumulate(current, capacity, m_isolated_workspace.kernel_k);
      accumulate(current, capacity, m_isolated_workspace.scratch);
      accumulate(current, capacity, m_isolated_workspace.potential_real);
      emit(core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kTransient,
           "pm_solver.isolated_open_workspace", current, capacity, capacity,
           "retained owned arrays; distributed root/backend allocations are reported separately");
    }

    {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kMpiBuffers,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kCommunication,
          .label = "pm_solver.density_exchange_logical_high_water",
          .current_size_bytes = 0U,
          .owned_capacity_bytes = 0U,
          .high_water_bytes = m_density_exchange_logical_high_water_bytes,
          .estimated_next_step_bytes = 0U,
          .estimate_only = true,
          .uncertainty_note = "logical non-owning PM density usage; physical storage is the attached TreePM communication arena or the standalone fallback arena",
      });
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kMpiBuffers,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = core::MemoryClass::kCommunication,
          .label = "pm_solver.plane_interpolation_exchange_logical_high_water",
          .current_size_bytes = 0U,
          .owned_capacity_bytes = 0U,
          .high_water_bytes = m_plane_interpolation_logical_high_water_bytes,
          .estimated_next_step_bytes = 0U,
          .estimate_only = true,
          .uncertainty_note = "logical non-owning PM interpolation usage; request/response reuse the same attached physical arena",
      });
      if (m_standalone_communication_arena != nullptr &&
          m_standalone_communication_arena->configured()) {
        const std::uint64_t capacity =
            m_standalone_communication_arena->capacityBytes();
        builder.addEntry(core::MemoryEntry{
            .subsystem = core::MemorySubsystem::kMpiBuffers,
            .lifetime = core::MemoryLifetime::kTransient,
            .memory_class = core::MemoryClass::kCommunication,
            .label = "pm_solver.communication_arena_fallback",
            .current_size_bytes = 0U,
            .owned_capacity_bytes = capacity,
            .high_water_bytes = capacity,
            .estimated_next_step_bytes = capacity,
            .uncertainty_note = "standalone PM-only bounded fallback owner; absent from externally arena-backed TreePM production",
        });
      }
    }

#if COSMOSIM_ENABLE_CUDA
    const auto emit_device = [&builder](std::string label, const core::DeviceBufferDouble& buffer) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kPmMesh,
          .lifetime = core::MemoryLifetime::kPersistent,
          .label = std::move(label),
          .current_size_bytes = static_cast<std::uint64_t>(buffer.sizeBytes()),
          .owned_capacity_bytes = static_cast<std::uint64_t>(buffer.sizeBytes()),
          .high_water_bytes = static_cast<std::uint64_t>(buffer.highWaterBytes()),
          .estimated_next_step_bytes = 0U,
          .uncertainty_note = "persistent CUDA allocation owned by PmSolver; CUDA/cuFFT runtime internals excluded",
      });
    };
    emit_device("pm_solver.cuda.pos_x", m_cuda_workspace.pos_x);
    emit_device("pm_solver.cuda.pos_y", m_cuda_workspace.pos_y);
    emit_device("pm_solver.cuda.pos_z", m_cuda_workspace.pos_z);
    emit_device("pm_solver.cuda.mass", m_cuda_workspace.mass);
    emit_device("pm_solver.cuda.density", m_cuda_workspace.density);
    emit_device("pm_solver.cuda.force_x", m_cuda_workspace.force_x);
    emit_device("pm_solver.cuda.force_y", m_cuda_workspace.force_y);
    emit_device("pm_solver.cuda.force_z", m_cuda_workspace.force_z);
    emit_device("pm_solver.cuda.accel_x", m_cuda_workspace.accel_x);
    emit_device("pm_solver.cuda.accel_y", m_cuda_workspace.accel_y);
    emit_device("pm_solver.cuda.accel_z", m_cuda_workspace.accel_z);
    emit_device("pm_solver.cuda.fft_workspace", m_cuda_workspace.fft_workspace);
    emit_device("pm_solver.cuda.green_filter_workspace", m_cuda_workspace.green_filter_workspace);
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kPmMesh,
        .lifetime = core::MemoryLifetime::kUnknown,
        .label = "pm_solver.cuda.plan_internal_workspace",
        .estimate_only = true,
        .uncertainty_note = "CUDA/cuFFT plan/runtime allocations are backend-owned and not measurable through the current API",
    });
#endif
  }

  PmGridShape m_shape;
  std::unordered_map<PlanKey, PlanResources, PlanKeyHasher> m_plan_cache;
  std::optional<PlanKey> m_active_key;
  std::size_t m_plan_build_count = 0;
  std::uint64_t m_next_distributed_exchange_epoch = 1;
  mutable GravityCommunicationArena* m_external_communication_arena = nullptr;
  mutable std::pmr::memory_resource* m_external_communication_resource = nullptr;
  mutable std::unique_ptr<GravityCommunicationArena> m_standalone_communication_arena;
  mutable std::uint64_t m_periodic_density_generation = 0U;
  mutable std::uint64_t m_density_exchange_logical_high_water_bytes = 0U;
  mutable std::uint64_t m_plane_interpolation_logical_high_water_bytes = 0U;
  IndexedTargetWorkspace m_indexed_target_workspace{};
  IsolatedOpenWorkspace m_isolated_workspace{};
#if COSMOSIM_ENABLE_CUDA
  CudaBackendWorkspace m_cuda_workspace{};
#endif
};

std::size_t PmGridShape::cellCount() const {
  return checkedProduct(checkedProduct(nx, ny, "PmGridShape::cellCount"), nz, "PmGridShape::cellCount");
}

bool PmGridShape::isValid() const {
  return nx > 0 && ny > 0 && nz > 0;
}

PmPlanResourcesMemoryEstimate estimatePmPlanResourcesMemory(
    PmGridShape shape,
    const parallel::PmSlabLayout& layout,
    core::PmDecompositionMode decomposition_mode,
    PmForceGradientBackend force_gradient_backend) {
  const PmPlanLocalStorageLayout storage_layout =
      determinePmPlanLocalStorageLayout(
          shape, layout, decomposition_mode, false);

  const auto checked_bytes = [](std::size_t count, std::size_t element_bytes,
                                std::string_view context) -> std::uint64_t {
    const std::size_t bytes = checkedProduct(count, element_bytes, context);
    return static_cast<std::uint64_t>(bytes);
  };

  PmPlanResourcesMemoryEstimate estimate;
  estimate.logical_local_complex_cells = static_cast<std::uint64_t>(storage_layout.logical_local_complex_cells);
  estimate.allocated_local_complex_cells = static_cast<std::uint64_t>(storage_layout.allocated_local_complex_cells);
  estimate.used_backend_allocation_query = storage_layout.used_backend_allocation_query;

  // In-place FFT: the padded real view aliases the complex allocation, so it
  // owns no bytes of its own (real_element_count = 2 x allocated complex).
  (void)storage_layout.real_element_count;
  estimate.real_array_bytes = 0U;
  estimate.real_view_aliases_complex = true;

  estimate.complex_spectral_array_bytes = checked_bytes(
      checkedProduct(
          force_gradient_backend == PmForceGradientBackend::kSpectral ? 2U : 1U,
          storage_layout.allocated_local_complex_cells,
                     "PM plan-memory complex spectral arrays"),
      sizeof(std::complex<double>),
      "PM plan-memory complex spectral bytes");
  estimate.scalar_spectral_array_bytes = 0U;
  const std::size_t nz_complex = shape.nz / 2U + 1U;
  const std::size_t axis_value_count = checkedProduct(
      2U,
      core::checkedSizeAdd(
          core::checkedSizeAdd(shape.nx, shape.ny, "PM plan-memory axis metadata"),
          nz_complex,
          "PM plan-memory axis metadata"),
      "PM plan-memory wave-number/window axis metadata");
  estimate.spectral_axis_metadata_bytes = checked_bytes(
      axis_value_count,
      sizeof(double),
      "PM plan-memory spectral axis metadata bytes");
  const auto checked_add_bytes = [](std::uint64_t lhs, std::uint64_t rhs,
                                    std::string_view context) {
    if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
      throw std::overflow_error(std::string(context) + " overflow");
    }
    return lhs + rhs;
  };
  const std::uint64_t spectral_bytes = checked_add_bytes(
      estimate.complex_spectral_array_bytes,
      estimate.spectral_axis_metadata_bytes,
      "PM plan-memory spectral byte sum");
  estimate.total_owned_bytes = checked_add_bytes(
      estimate.real_array_bytes,
      spectral_bytes,
      "PM plan-memory total byte sum");
  // Finite-difference gradients exchange two x planes per side with the
  // neighboring slab ranks; four halo planes plus two send-staging planes
  // are retained per plan (no per-solve allocation).
  if (force_gradient_backend == PmForceGradientBackend::kFiniteDifference4 && layout.world_size > 1) {
    estimate.finite_difference_halo_bytes = checked_bytes(
        checkedProduct(6U, checkedProduct(shape.ny, shape.nz, "PM finite-difference halo plane"),
                       "PM finite-difference halo planes"),
        sizeof(double), "PM finite-difference halo bytes");
    estimate.total_owned_bytes = checked_add_bytes(
        estimate.total_owned_bytes, estimate.finite_difference_halo_bytes,
        "PM plan-memory total byte sum");
  }
  return estimate;
}

void PmProfiler::reset() {
  m_totals = {};
}

void PmProfiler::append(const PmProfileEvent& event) {
  m_totals.bytes_moved += event.bytes_moved;
  m_totals.routed_density_records += event.routed_density_records;
  m_totals.routed_force_requests += event.routed_force_requests;
  m_totals.routed_potential_requests += event.routed_potential_requests;
  m_totals.routed_density_peer_count += event.routed_density_peer_count;
  m_totals.routed_force_peer_count += event.routed_force_peer_count;
  m_totals.routed_potential_peer_count += event.routed_potential_peer_count;
  m_totals.routed_mpi_bytes_sent += event.routed_mpi_bytes_sent;
  m_totals.routed_mpi_bytes_received += event.routed_mpi_bytes_received;
  m_totals.routed_send_buffer_high_water_bytes = std::max(
      m_totals.routed_send_buffer_high_water_bytes, event.routed_send_buffer_high_water_bytes);
  m_totals.routed_receive_buffer_high_water_bytes = std::max(
      m_totals.routed_receive_buffer_high_water_bytes, event.routed_receive_buffer_high_water_bytes);
  m_totals.routed_combined_buffer_high_water_bytes = std::max(
      m_totals.routed_combined_buffer_high_water_bytes, event.routed_combined_buffer_high_water_bytes);
  m_totals.routed_workspace_high_water_bytes = std::max(
      m_totals.routed_workspace_high_water_bytes, event.routed_workspace_high_water_bytes);
  m_totals.force_halo_cache_hits += event.force_halo_cache_hits;
  m_totals.spectral_operator_rebuilds += event.spectral_operator_rebuilds;
  m_totals.isolated_open_root_workspace_estimate_bytes =
      std::max(m_totals.isolated_open_root_workspace_estimate_bytes,
               event.isolated_open_root_workspace_estimate_bytes);
  m_totals.isolated_open_root_workspace_limit_bytes =
      std::max(m_totals.isolated_open_root_workspace_limit_bytes,
               event.isolated_open_root_workspace_limit_bytes);
  m_totals.isolated_open_gather_bytes += event.isolated_open_gather_bytes;
  m_totals.density_routing_wait_ms += event.density_routing_wait_ms;
  m_totals.halo_exchange_ms += event.halo_exchange_ms;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    m_totals.gradient_factor_ms[axis] += event.gradient_factor_ms[axis];
    m_totals.inverse_axis_ms[axis] += event.inverse_axis_ms[axis];
    m_totals.normalization_axis_ms[axis] += event.normalization_axis_ms[axis];
  }
  m_totals.assign_ms += event.assign_ms;
  m_totals.fft_forward_ms += event.fft_forward_ms;
  m_totals.poisson_ms += event.poisson_ms;
  m_totals.gradient_ms += event.gradient_ms;
  m_totals.fft_inverse_ms += event.fft_inverse_ms;
  m_totals.fft_transpose_ms += event.fft_transpose_ms;
  m_totals.interpolate_ms += event.interpolate_ms;
  m_totals.total_ms += event.total_ms;
  m_totals.routed_mpi_wait_ms += event.routed_mpi_wait_ms;
  m_totals.transfer_h2d_ms += event.transfer_h2d_ms;
  m_totals.transfer_d2h_ms += event.transfer_d2h_ms;
  m_totals.device_kernel_ms += event.device_kernel_ms;
  m_totals.fft_transpose_bytes += event.fft_transpose_bytes;
}

const PmProfileEvent& PmProfiler::totals() const {
  return m_totals;
}

PmGridStorage::PmGridStorage(PmGridShape shape)
    : PmGridStorage(
          shape,
          parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, /*world_size=*/1, /*world_rank=*/0)) {}

PmGridStorage::PmGridStorage(PmGridShape shape, parallel::PmSlabLayout layout)
    : m_shape(shape),
      m_layout(std::move(layout)),
      m_force_x(m_layout.localCellCount(), 0.0),
      m_force_y(m_layout.localCellCount(), 0.0),
      m_force_z(m_layout.localCellCount(), 0.0) {
  if (!m_shape.isValid()) {
    throw std::invalid_argument("PM grid shape must be valid");
  }
  if (!m_layout.isValid()) {
    throw std::invalid_argument("PM slab layout must be valid");
  }
  if (m_layout.global_nx != m_shape.nx || m_layout.global_ny != m_shape.ny || m_layout.global_nz != m_shape.nz) {
    throw std::invalid_argument("PM slab layout global shape must match PM grid shape");
  }
}

const PmGridShape& PmGridStorage::shape() const {
  return m_shape;
}

const parallel::PmSlabLayout& PmGridStorage::slabLayout() const {
  return m_layout;
}

bool PmGridStorage::ownsFullDomain() const noexcept {
  return m_layout.ownsFullDomain();
}

std::size_t PmGridStorage::localCellCount() const {
  return m_layout.localCellCount();
}

std::uint64_t PmGridStorage::storageGeneration() const noexcept {
  return m_storage_generation;
}

std::span<double> PmGridStorage::density() {
  ensureDensityStorage();
  return m_density;
}

std::span<const double> PmGridStorage::density() const {
  return m_density;
}

bool PmGridStorage::hasDensityStorage() const noexcept {
  return m_density.size() == localCellCount();
}

void PmGridStorage::ensureDensityStorage() {
  if (hasDensityStorage()) {
    return;
  }
  if (m_storage_generation == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("PM grid storage generation overflow");
  }
  m_density.assign(localCellCount(), 0.0);
  ++m_storage_generation;
}

void PmGridStorage::releaseDensityStorage() {
  if (m_density.empty()) {
    return;
  }
  if (m_storage_generation == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("PM grid storage generation overflow");
  }
  std::vector<double>{}.swap(m_density);
  ++m_storage_generation;
}

std::span<double> PmGridStorage::potential() {
  ensurePotentialStorage();
  return m_potential;
}

std::span<const double> PmGridStorage::potential() const {
  return m_potential;
}

bool PmGridStorage::hasPotentialStorage() const noexcept {
  return m_potential.size() == localCellCount();
}

void PmGridStorage::ensurePotentialStorage() {
  if (hasPotentialStorage()) {
    return;
  }
  m_potential.assign(localCellCount(), 0.0);
}

std::span<double> PmGridStorage::force_x() {
  return m_force_x;
}

std::span<const double> PmGridStorage::force_x() const {
  return m_force_x;
}

std::span<double> PmGridStorage::force_y() {
  return m_force_y;
}

std::span<const double> PmGridStorage::force_y() const {
  return m_force_y;
}

std::span<double> PmGridStorage::force_z() {
  return m_force_z;
}

std::span<const double> PmGridStorage::force_z() const {
  return m_force_z;
}

void PmGridStorage::clearForceHaloCache() {
  // Invalidation is intentionally capacity-preserving.  The six final force
  // halo lanes are legitimate phase-retained state and are reused across PM
  // refreshes; only validity/publication metadata is cleared here.
  m_force_halo_cache.valid = false;
  m_force_halo_cache.halo_depth_x = 0U;
  m_force_halo_cache.left_peer_rank = -1;
  m_force_halo_cache.right_peer_rank = -1;
  m_force_halo_cache.exchange_sequence = 0U;
}

void PmGridStorage::beginForceHaloCacheRefresh(std::size_t halo_depth_x) {
  m_force_halo_cache.valid = false;
  m_force_halo_cache.halo_depth_x = 0U;
  m_force_halo_cache.left_peer_rank = -1;
  m_force_halo_cache.right_peer_rank = -1;
  m_force_halo_cache.exchange_sequence = 0U;
  const std::size_t plane_size = core::checkedSizeMultiply(
      m_layout.global_ny, m_layout.global_nz,
      "PM force halo cache plane size");
  const std::size_t value_count = core::checkedSizeMultiply(
      halo_depth_x, plane_size,
      "PM force halo cache value count");
  m_force_halo_cache.left_force_x.resize(value_count);
  m_force_halo_cache.left_force_y.resize(value_count);
  m_force_halo_cache.left_force_z.resize(value_count);
  m_force_halo_cache.right_force_x.resize(value_count);
  m_force_halo_cache.right_force_y.resize(value_count);
  m_force_halo_cache.right_force_z.resize(value_count);
}

PmGridStorage::ForceHaloWriteView PmGridStorage::forceHaloWriteView(
    std::size_t component) {
  if (m_force_halo_cache.valid) {
    throw std::logic_error("PM force halo cache cannot be mutated after publication");
  }
  switch (component) {
    case 0U:
      return ForceHaloWriteView{
          std::span<double>(m_force_halo_cache.left_force_x),
          std::span<double>(m_force_halo_cache.right_force_x)};
    case 1U:
      return ForceHaloWriteView{
          std::span<double>(m_force_halo_cache.left_force_y),
          std::span<double>(m_force_halo_cache.right_force_y)};
    case 2U:
      return ForceHaloWriteView{
          std::span<double>(m_force_halo_cache.left_force_z),
          std::span<double>(m_force_halo_cache.right_force_z)};
    default:
      throw std::out_of_range("PM force halo cache component must be 0, 1, or 2");
  }
}

void PmGridStorage::commitForceHaloCacheRefresh(
    std::size_t halo_depth_x,
    int left_peer_rank,
    int right_peer_rank,
    std::uint64_t exchange_sequence) {
  const std::size_t plane_size = core::checkedSizeMultiply(
      m_layout.global_ny, m_layout.global_nz,
      "PM force halo cache commit plane size");
  const std::size_t expected_values = core::checkedSizeMultiply(
      halo_depth_x, plane_size,
      "PM force halo cache commit value count");
  const auto require_extent = [&](const std::vector<double>& values) {
    if (values.size() != expected_values) {
      throw std::invalid_argument(
          "PM force halo cache component extent changed during refresh");
    }
  };
  require_extent(m_force_halo_cache.left_force_x);
  require_extent(m_force_halo_cache.left_force_y);
  require_extent(m_force_halo_cache.left_force_z);
  require_extent(m_force_halo_cache.right_force_x);
  require_extent(m_force_halo_cache.right_force_y);
  require_extent(m_force_halo_cache.right_force_z);
  if (halo_depth_x == 0U || plane_size == 0U) {
    clearForceHaloCache();
    return;
  }
  m_force_halo_cache.halo_depth_x = halo_depth_x;
  m_force_halo_cache.left_peer_rank = left_peer_rank;
  m_force_halo_cache.right_peer_rank = right_peer_rank;
  m_force_halo_cache.exchange_sequence = exchange_sequence;
  m_force_halo_cache.valid = true;
}

bool PmGridStorage::hasForceHaloCache() const noexcept {
  return m_force_halo_cache.valid;
}

bool PmGridStorage::tryLoadForceFromHalo(
    std::size_t global_x,
    std::size_t global_y,
    std::size_t global_z,
    double& force_x_value,
    double& force_y_value,
    double& force_z_value) const {
  if (!m_force_halo_cache.valid || global_y >= m_layout.global_ny || global_z >= m_layout.global_nz ||
      m_force_halo_cache.halo_depth_x == 0 || m_layout.global_nx == 0) {
    return false;
  }
  const std::size_t depth = m_force_halo_cache.halo_depth_x;
  const std::size_t plane_size = m_layout.global_ny * m_layout.global_nz;
  const auto plane_offset = [&](std::size_t halo_x) {
    return halo_x * plane_size + global_y * m_layout.global_nz + global_z;
  };

  if (m_force_halo_cache.left_peer_rank >= 0) {
    for (std::size_t halo_x = 0; halo_x < depth; ++halo_x) {
      const std::size_t halo_global_x = (m_layout.owned_x.begin_x + m_layout.global_nx - depth + halo_x) %
          m_layout.global_nx;
      if (halo_global_x == global_x) {
        const std::size_t index = plane_offset(halo_x);
        force_x_value = m_force_halo_cache.left_force_x[index];
        force_y_value = m_force_halo_cache.left_force_y[index];
        force_z_value = m_force_halo_cache.left_force_z[index];
        return true;
      }
    }
  }

  if (m_force_halo_cache.right_peer_rank >= 0) {
    for (std::size_t halo_x = 0; halo_x < depth; ++halo_x) {
      const std::size_t halo_global_x = (m_layout.owned_x.end_x + halo_x) % m_layout.global_nx;
      if (halo_global_x == global_x) {
        const std::size_t index = plane_offset(halo_x);
        force_x_value = m_force_halo_cache.right_force_x[index];
        force_y_value = m_force_halo_cache.right_force_y[index];
        force_z_value = m_force_halo_cache.right_force_z[index];
        return true;
      }
    }
  }
  return false;
}

std::size_t PmGridStorage::linearIndex(std::size_t ix, std::size_t iy, std::size_t iz) const {
  return m_layout.localLinearIndex(ix, iy, iz);
}

void PmGridStorage::clear() {
  if (m_storage_generation == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("PM grid storage generation overflow");
  }
  ++m_storage_generation;
  if (!m_density.empty()) {
    std::fill(m_density.begin(), m_density.end(), 0.0);
  }
  if (!m_potential.empty()) {
    std::fill(m_potential.begin(), m_potential.end(), 0.0);
  }
  std::fill(m_force_x.begin(), m_force_x.end(), 0.0);
  std::fill(m_force_y.begin(), m_force_y.end(), 0.0);
  std::fill(m_force_z.begin(), m_force_z.end(), 0.0);
  clearForceHaloCache();
}

void PmGridStorage::appendMemoryReport(core::MemoryReportBuilder& builder) const {
  const auto add = [&builder](std::string label, const auto& container) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kPmMesh,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .current_size_bytes = core::currentSizeBytesForContainer(container),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = bytes,
                                       .estimated_next_step_bytes = 0U,
                                       .uncertainty_note = "retained PM-grid capacity is an observed allocation high-water; next-step requirement is not predicted"});
  };
  add("pm_mesh.density", m_density);
  add("pm_mesh.potential", m_potential);
  add("pm_mesh.force_x", m_force_x);
  add("pm_mesh.force_y", m_force_y);
  add("pm_mesh.force_z", m_force_z);
  const auto add_force_halo_cache = [&builder](
      std::string label, const auto& container) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kPmMesh,
        .lifetime = core::MemoryLifetime::kPersistent,
        .memory_class = core::MemoryClass::kPersistentCache,
        .label = std::move(label),
        .current_size_bytes = core::currentSizeBytesForContainer(container),
        .owned_capacity_bytes = bytes,
        .high_water_bytes = bytes,
        .estimated_next_step_bytes = 0U,
        .uncertainty_note =
            "retained final force-halo cache; survives communication-arena reset until PM interpolation completes",
    });
  };
  add_force_halo_cache("pm_mesh.force_halo_left_x", m_force_halo_cache.left_force_x);
  add_force_halo_cache("pm_mesh.force_halo_left_y", m_force_halo_cache.left_force_y);
  add_force_halo_cache("pm_mesh.force_halo_left_z", m_force_halo_cache.left_force_z);
  add_force_halo_cache("pm_mesh.force_halo_right_x", m_force_halo_cache.right_force_x);
  add_force_halo_cache("pm_mesh.force_halo_right_y", m_force_halo_cache.right_force_y);
  add_force_halo_cache("pm_mesh.force_halo_right_z", m_force_halo_cache.right_force_z);
}

PmSolver::PmSolver(PmGridShape shape) : m_shape(shape), m_impl(std::make_unique<Impl>(shape)) {
  if (!shape.isValid()) {
    throw std::invalid_argument("PM solver requires valid shape");
  }
}

PmSolver::~PmSolver() = default;
PmSolver::PmSolver(PmSolver&&) noexcept = default;
PmSolver& PmSolver::operator=(PmSolver&&) noexcept = default;

const PmGridShape& PmSolver::shape() const {
  return m_shape;
}

void PmSolver::appendMemoryReport(core::MemoryReportBuilder& builder) const {
  m_impl->appendOwnedMemoryReport(builder);
}

void PmSolver::attachCommunicationArena(GravityCommunicationArena* arena) const {
  m_impl->attachCommunicationArena(arena);
}

void PmSolver::detachCommunicationArena() const noexcept {
  m_impl->detachCommunicationArena();
}

void PmSolver::shutdownBackendResources() {
  m_impl->shutdownBackendResources();
}

void PmSolver::assignDensity(
    PmGridStorage& grid,
    const PmMassSourceView& source_view,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
  assignDensity(
      grid,
      source_view.pos_x_comoving,
      source_view.pos_y_comoving,
      source_view.pos_z_comoving,
      source_view.mass_code,
      options,
      profile);
}

void PmSolver::assignDensity(
    PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
  assignDensityImpl(
      grid, pos_x, pos_y, pos_z, mass, options,
      /*periodic_fft_backed=*/false, profile);
}

void PmSolver::assignDensityPeriodicFftBacked(
    PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
  assignDensityImpl(
      grid, pos_x, pos_y, pos_z, mass, options,
      /*periodic_fft_backed=*/true, profile);
}

void PmSolver::assignDensityImpl(
    PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const PmSolveOptions& options,
    bool periodic_fft_backed,
    PmProfileEvent* profile) const {
  const PmDecompositionView decomposition_view(
      m_shape, grid.slabLayout(), options.decomposition_mode);
#if COSMOSIM_ENABLE_MPI
  double distributed_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  // Public PM entry is collective whenever the hosting MPI world has peers.
  // The vote precedes any layout-based branch, so divergent rank-local layout
  // metadata cannot let one rank throw while another enters an exchange.
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kAssignDensity,
        rank_local_serial_layout,
        m_shape,
        options,
        distributed_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::assignDensity distributed API preflight",
        mpi_world_rank,
        mpi_world_size,
        distributed_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in assignDensity");
          }
          if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != mass.size()) {
            throw std::invalid_argument("Particle coordinate/mass spans must match in assignDensity");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::assignDensity requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::assignDensity slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
            throw std::invalid_argument("PmSolver::assignDensity mesh dimensions exceed fixed-width routing indices");
          }
        });
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in assignDensity");
  }
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != mass.size()) {
    throw std::invalid_argument("Particle coordinate/mass spans must match in assignDensity");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::assignDensity requires a valid PM slab layout");
  }
  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if !COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    throw std::invalid_argument(
        "PmSolver::assignDensity distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
#else
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument(
          "PmSolver::assignDensity slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::assignDensity");
  }
#endif

  if (periodic_fft_backed &&
      options.boundary_condition != PmBoundaryCondition::kPeriodic) {
    throw std::invalid_argument(
        "FFT-backed PM density storage is valid only for periodic solves");
  }

  const auto start = std::chrono::steady_clock::now();
  Impl::PlanResources* density_plan = nullptr;
  std::span<double> density_storage;
  std::size_t density_z_stride = m_shape.nz;
  if (periodic_fft_backed) {
    // The production specialization has one physical density owner. Release
    // any compact compatibility capacity retained by an earlier generic use
    // before publishing a new FFT-backed density generation.
    m_impl->invalidatePeriodicDensityForGrid(&grid);
    grid.releaseDensityStorage();
    Impl::PlanResources& plan =
        m_impl->planForLayout(grid.slabLayout(), options.decomposition_mode);
    density_plan = &plan;
    plan.periodic_density_ready = false;
    plan.periodic_density_grid = nullptr;
    plan.periodic_density_generation = 0U;
    plan.periodic_density_grid_generation = 0U;
    density_storage = plan.realView();
    density_z_stride = plan.real_z_stride;
    std::fill(density_storage.begin(), density_storage.end(), 0.0);
  } else {
    m_impl->invalidatePeriodicDensityForGrid(&grid);
    grid.ensureDensityStorage();
    density_storage = grid.density();
    std::fill(density_storage.begin(), density_storage.end(), 0.0);
  }

  const auto density_index = [&](
      std::size_t global_ix,
      std::size_t global_iy,
      std::size_t global_iz) -> std::size_t {
    if (!grid.slabLayout().ownsGlobalX(global_ix) ||
        global_iy >= m_shape.ny || global_iz >= m_shape.nz) {
      throw std::out_of_range(
          "PM density destination index is outside the owned physical mesh");
    }
    const std::size_t local_ix = global_ix - grid.slabLayout().owned_x.begin_x;
    return (local_ix * m_shape.ny + global_iy) * density_z_stride + global_iz;
  };

  const BoxLengths lengths = effectiveBoxLengths(options);
  const double inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  const double inv_dy = static_cast<double>(m_shape.ny) / lengths.ly;
  const double inv_dz = static_cast<double>(m_shape.nz) / lengths.lz;
  if (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
      m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
      m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw std::invalid_argument("PmSolver::assignDensity mesh dimensions exceed fixed-width routing indices");
  }
  std::optional<std::size_t> first_non_finite_source_index;
  for (std::size_t p = 0; p < mass.size(); ++p) {
    if (!std::isfinite(pos_x[p]) || !std::isfinite(pos_y[p]) || !std::isfinite(pos_z[p]) ||
        !std::isfinite(mass[p]) || mass[p] < 0.0) {
      first_non_finite_source_index = p;
      break;
    }
  }
  if (!distributed_slabs && first_non_finite_source_index.has_value()) {
    throw std::invalid_argument(
        "PmSolver::assignDensity requires finite particle coordinates and finite non-negative masses; particle_index=" +
        std::to_string(*first_non_finite_source_index));
  }

  const auto accumulate_owned = [&](const PmDensityContributionRecord& record) {
    if (record.global_ix >= m_shape.nx || record.global_iy >= m_shape.ny || record.global_iz >= m_shape.nz) {
      throw std::invalid_argument("PmSolver::assignDensity received out-of-range PM contribution record");
    }
    if (!grid.slabLayout().ownsGlobalX(record.global_ix)) {
      throw std::invalid_argument("PmSolver::assignDensity received contribution for non-owned PM slab x-index");
    }
    if (!std::isfinite(record.mass_contribution) || record.mass_contribution < 0.0) {
      throw std::invalid_argument(
          "PmSolver::assignDensity received a non-finite or negative mass contribution");
    }
    density_storage[density_index(
        record.global_ix, record.global_iy, record.global_iz)] +=
        record.mass_contribution;
  };

  if (!distributed_slabs) {
    for (std::size_t p = 0; p < mass.size(); ++p) {
      const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
      if (!periodic &&
          (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
        continue;
      }
      const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
      const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
      const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;

      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        if (!periodic && (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<int>(m_shape.nx))) {
          continue;
        }
        const std::size_t ix = periodic ? wrapIndex(sx.offsets[dx], m_shape.nx) : static_cast<std::size_t>(sx.offsets[dx]);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic && (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<int>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic ? wrapIndex(sy.offsets[dy], m_shape.ny) : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic && (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<int>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic ? wrapIndex(sz.offsets[dz], m_shape.nz) : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            accumulate_owned(PmDensityContributionRecord{
                .global_ix = static_cast<std::uint32_t>(ix),
                .global_iy = static_cast<std::uint32_t>(iy),
                .global_iz = static_cast<std::uint32_t>(iz),
                .mass_contribution = mass[p] * weight,
            });
          }
        }
      }
    }
  } else {
#if COSMOSIM_ENABLE_MPI
    const int world_size = grid.slabLayout().world_size;
    const int world_rank = grid.slabLayout().world_rank;
    double routed_mpi_wait_ms = distributed_preflight_mpi_wait_ms;
    const int local_sources_valid = first_non_finite_source_index.has_value() ? 0 : 1;
    int all_sources_valid = 0;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(&local_sources_valid, &all_sources_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD),
          "PmSolver::assignDensity source validation MPI_Allreduce");
    });
    if (all_sources_valid == 0) {
      const std::string local_detail = first_non_finite_source_index.has_value()
          ? " local_first_invalid_particle_index=" + std::to_string(*first_non_finite_source_index)
          : " local_sources_valid=true";
      throw std::invalid_argument(
          "PmSolver::assignDensity rejected non-finite coordinates or masses on at least one MPI rank;" +
          local_detail);
    }

    std::uint64_t exchange_epoch = 0U;
    runPmCoordinatedPhase(
        "PmSolver::assignDensity exchange epoch allocation",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() { exchange_epoch = m_impl->nextDistributedExchangeEpoch(); });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      validatePmExchangeEpochConsensus(exchange_epoch, "PmSolver::assignDensity");
    });

    std::optional<GravityCommunicationArena::Lease> standalone_communication_lease;
    std::optional<Impl::DensityExchangeBuffers> exchange_storage;
    Impl::DensityExchangeBuffers* exchange_ptr = nullptr;
    PmRoutingCapacityModel routing_capacity{};
    std::size_t routing_buffer_limit = 0U;
    std::size_t particles_per_round = 0U;
    const bool routing_periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    const auto density_grid_x = [&](std::size_t p) -> std::optional<double> {
      if (!routing_periodic &&
          (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
        return std::nullopt;
      }
      return (routing_periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
    };
    // Routing state machine (collective, every rank, fixed order):
    //   (1) coordinated workload phase: withdraw any plan, build the
    //       rank-independent round model, measure the current local remote
    //       workload and check it against the plan without throwing;
    //   (2) agreePmRoutingWorkload: one MAX allreduce over current maxima and
    //       plan presence/violation lanes (uniform failure on mismatch);
    //   (3) coordinated preflight phase: lease and reserve from the agreed
    //       current workload;
    //   (4) bounded data rounds.
    bool local_plan_present = false;
    bool local_plan_violated = false;
    PmRemotePlaneWorkload local_workload{};
    runPmCoordinatedPhase(
        "PmSolver::assignDensity routing workload",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
          const Impl::RoutingWorkloadPlanCache plan = m_impl->m_density_routing_plan;
          m_impl->m_density_routing_plan = {};  // withdrawn whatever the outcome
          local_plan_present = plan.valid;
          routing_capacity = pmPlannedRoutingModel(
              PmRoutingPhaseKind::kDensity, world_size, options.routing_exchange_batch_bytes);
          if (routing_capacity.max_send_payload_bytes >
              static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            throw std::overflow_error("PM density routing buffer limit exceeds size_t capacity");
          }
          routing_buffer_limit =
              static_cast<std::size_t>(routing_capacity.max_send_payload_bytes);
          particles_per_round = pmRoutingParticlesPerRound(
              routing_capacity.effective_per_peer_payload_bytes,
              k_pm_density_plane_wire_bytes);
          local_workload = measurePmRemotePlaneWorkload(
              mass.size(), particles_per_round, world_rank, world_size, routing_periodic,
              options.assignment_scheme, m_shape, decomposition_view, density_grid_x);
          if (local_plan_present) {
            local_plan_violated = !pmRoutingPlanBoundsWorkload(
                plan.key, plan.global_max_count, plan.global_max_peer_records, plan.local_max_round_records,
                makePmRoutingPlanKey(m_shape, grid.slabLayout(), options),
                static_cast<std::uint64_t>(mass.size()), local_workload);
          }
        });
    // Replaces the former max-source reduction (one MAX over agreed lanes).
    const PmRoutingAgreement routing_agreement = agreePmRoutingWorkload(
        "PmSolver::assignDensity", world_rank, world_size,
        static_cast<std::uint64_t>(mass.size()), local_workload.max_peer_records,
        local_plan_present, local_plan_violated, routed_mpi_wait_ms);
    const std::uint64_t global_max_source_count = routing_agreement.global_max_count;
    const std::uint64_t global_max_peer_records = routing_agreement.global_max_peer_records;
    const PmRoutingReservation wire_reservation = pmRoutingReservation(
        PmRoutingPhaseKind::kDensity, routing_capacity, world_size,
        local_workload.max_round_records, global_max_peer_records);
    runPmCoordinatedPhase(
        "PmSolver::assignDensity routing preflight",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
          std::pmr::memory_resource* communication_resource = m_impl->communicationResource(
              GravityCommunicationArena::Phase::kPmDensity,
              standalone_communication_lease,
              wire_reservation.phase_arena_bytes);
          exchange_storage.emplace(
              communication_resource, static_cast<std::size_t>(world_size));
          exchange_ptr = &*exchange_storage;
          // The phase arena is monotonic/resettable. Reserve both bounded wire
          // maxima once so later exact-size resizes never accumulate abandoned
          // growth blocks inside the same logical lease. The reservations are
          // workload bounds (see pmRoutingReservation), not the configured
          // (P-1) x per-peer ceiling.
          exchange_ptr->send_wire.reserve(static_cast<std::size_t>(wire_reservation.send_bytes));
          exchange_ptr->recv_wire.reserve(static_cast<std::size_t>(wire_reservation.receive_bytes));
        });
    auto& exchange = *exchange_ptr;
    // Fixed routing metadata of this phase (exact lane formula; the planned
    // round model is built on the same value on every rank).
    const std::uint64_t routing_metadata_capacity =
        pmRoutingMetadataBytes(PmRoutingPhaseKind::kDensity, world_size);
    const std::uint64_t local_source_count = static_cast<std::uint64_t>(mass.size());
    const std::uint64_t round_capacity = static_cast<std::uint64_t>(particles_per_round);
    const std::uint64_t round_count = global_max_source_count == 0U
        ? 0U
        : 1U + (global_max_source_count - 1U) / round_capacity;
    const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    std::uint32_t next_record_sequence = 0U;

    const auto accumulate_plane_record = [&](const PmDensityPlaneRecord& record) {
      if (record.plane_count == 0U || record.plane_count > 3U ||
          record.destination_rank != static_cast<std::uint32_t>(world_rank) ||
          record.exchange_epoch != exchange_epoch ||
          !std::isfinite(record.y_grid) || !std::isfinite(record.z_grid) ||
          !std::isfinite(record.mass_code) || record.mass_code < 0.0) {
        throw std::invalid_argument("PmSolver::assignDensity received invalid bounded density-plane record");
      }
      const PmAxisStencil1d sy = makeAxisStencil(record.y_grid, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(record.z_grid, options.assignment_scheme);
      for (std::size_t plane = 0; plane < static_cast<std::size_t>(record.plane_count); ++plane) {
        const std::size_t ix = static_cast<std::size_t>(record.global_ix[plane]);
        if (ix >= m_shape.nx || !grid.slabLayout().ownsGlobalX(ix) ||
            !std::isfinite(record.x_weight[plane]) || record.x_weight[plane] < 0.0) {
          throw std::invalid_argument("PmSolver::assignDensity density-plane x ownership/weight is invalid");
        }
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic &&
              (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic &&
                (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = record.x_weight[plane] * sy.weights[dy] * sz.weights[dz];
            density_storage[density_index(ix, iy, iz)] +=
                record.mass_code * weight;
          }
        }
      }
    };

    std::uint64_t logical_records_sent = 0U;
    std::uint64_t peers_touched = 0U;
    std::uint64_t total_wire_sent = 0U;
    std::uint64_t total_wire_received = 0U;
    std::uint64_t send_high_water = 0U;
    std::uint64_t receive_high_water = 0U;
    std::uint64_t combined_high_water = 0U;
    std::uint64_t workspace_high_water = routing_metadata_capacity;

    for (std::uint64_t round = 0U; round < round_count; ++round) {
      const std::uint64_t begin_u64 = round * round_capacity;
      const std::uint64_t end_u64 = std::min(
          local_source_count,
          begin_u64 > std::numeric_limits<std::uint64_t>::max() - round_capacity
              ? local_source_count
              : begin_u64 + round_capacity);
      const std::size_t begin = static_cast<std::size_t>(begin_u64);
      const std::size_t end = static_cast<std::size_t>(end_u64);

      std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
      std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
      std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
      std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
      std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
      std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
      std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
      std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
      exchange.send_wire.clear();
      exchange.recv_wire.clear();

      runPmCoordinatedPhase(
          "PmSolver::assignDensity bounded density send preparation",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            // Pass 1 counts one compact record per particle/destination x-plane group.
            for (std::size_t p = begin; p < end; ++p) {
              if (!periodic &&
                  (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
                   !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
                   !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
                continue;
              }
              const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
              const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const int destination_rank = groups[group_i].destination_rank;
                if (destination_rank == world_rank) {
                  continue;
                }
                int& count = exchange.send_counts[static_cast<std::size_t>(destination_rank)];
                if (count == std::numeric_limits<int>::max()) {
                  throw std::overflow_error("PM density per-peer record count exceeds MPI int capacity");
                }
                ++count;
              }
            }

            std::size_t total_send_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = static_cast<std::size_t>(
                  exchange.send_counts[static_cast<std::size_t>(rank)]);
              const std::uint64_t peer_bytes = checkedBytesForCount(
                  count, k_pm_density_plane_wire_bytes, "PM density per-peer routing batch");
              if (peer_bytes > routing_capacity.effective_per_peer_payload_bytes) {
                throw std::runtime_error("PM density routing exceeded effective per-peer batch high-water");
              }
              exchange.send_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_send_records,
                      k_pm_density_plane_wire_bytes,
                      "PmSolver::assignDensity bounded density exchange",
                      "send record displacement",
                      rank);
              total_send_records = checkedMpiRecordTotal(
                  total_send_records,
                  count,
                  k_pm_density_plane_wire_bytes,
                  "PmSolver::assignDensity bounded density exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.send_wire,
                checkedPmWireByteCount(
                    total_send_records,
                    k_pm_density_plane_wire_bytes,
                    "PmSolver::assignDensity bounded density exchange"),
                routing_buffer_limit,
                "PmSolver::assignDensity send buffer");

            for (int rank = 0; rank < world_size; ++rank) {
              exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
                  exchange.send_displs[static_cast<std::size_t>(rank)]);
            }

            // Pass 2 deposits owner-local x planes immediately and encodes only
            // remote x-plane groups directly into the bounded wire buffer.
            for (std::size_t p = begin; p < end; ++p) {
              if (!periodic &&
                  (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
                   !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
                   !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
                continue;
              }
              const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
              const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
              const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;
              const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const auto& group = groups[group_i];
                PmDensityPlaneRecord record{
                    .origin_rank = static_cast<std::uint32_t>(world_rank),
                    .destination_rank = static_cast<std::uint32_t>(group.destination_rank),
                    .record_sequence = next_record_sequence,
                    .plane_count = group.plane_count,
                    .global_ix = group.global_ix,
                    .exchange_epoch = exchange_epoch,
                    .y_grid = y,
                    .z_grid = z,
                    .mass_code = mass[p],
                    .x_weight = group.x_weight,
                };
                if (group.destination_rank == world_rank) {
                  accumulate_plane_record(record);
                  continue;
                }
                if (next_record_sequence == std::numeric_limits<std::uint32_t>::max()) {
                  throw std::overflow_error("PmSolver::assignDensity routing sequence exceeds uint32 capacity");
                }
                auto& cursor = exchange.cursor[static_cast<std::size_t>(group.destination_rank)];
                const std::size_t byte_offset = cursor * k_pm_density_plane_wire_bytes;
                encodePmDensityPlaneRecord(
                    record,
                    std::span<std::uint8_t>(exchange.send_wire).subspan(
                        byte_offset, k_pm_density_plane_wire_bytes));
                ++cursor;
                ++next_record_sequence;
              }
            }
          });

      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoall(
                exchange.send_counts.data(), 1, MPI_INT,
                exchange.recv_counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD),
            "PmSolver::assignDensity count MPI_Alltoall");
      });

      runPmCoordinatedPhase(
          "PmSolver::assignDensity bounded density receive layout",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            std::size_t total_recv_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = checkedMpiReceivedRecordCount(
                  exchange.recv_counts[static_cast<std::size_t>(rank)],
                  k_pm_density_plane_wire_bytes,
                  "PmSolver::assignDensity bounded density exchange",
                  rank);
              exchange.recv_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_recv_records,
                      k_pm_density_plane_wire_bytes,
                      "PmSolver::assignDensity bounded density exchange",
                      "received record displacement",
                      rank);
              total_recv_records = checkedMpiRecordTotal(
                  total_recv_records,
                  count,
                  k_pm_density_plane_wire_bytes,
                  "PmSolver::assignDensity bounded density exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.recv_wire,
                checkedPmWireByteCount(
                    total_recv_records,
                    k_pm_density_plane_wire_bytes,
                    "PmSolver::assignDensity bounded density exchange"),
                routing_buffer_limit,
                "PmSolver::assignDensity receive buffer");
            checkedMpiRecordLayoutToByteLayout(
                exchange.send_counts,
                exchange.send_displs,
                k_pm_density_plane_wire_bytes,
                exchange.send_counts_bytes,
                exchange.send_displs_bytes,
                "PmSolver::assignDensity bounded density send MPI_Alltoallv");
            checkedMpiRecordLayoutToByteLayout(
                exchange.recv_counts,
                exchange.recv_displs,
                k_pm_density_plane_wire_bytes,
                exchange.recv_counts_bytes,
                exchange.recv_displs_bytes,
                "PmSolver::assignDensity bounded density receive MPI_Alltoallv");
          });

      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.send_wire),
                exchange.send_counts_bytes.data(),
                exchange.send_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.recv_wire),
                exchange.recv_counts_bytes.data(),
                exchange.recv_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::assignDensity payload MPI_Alltoallv");
      });

      std::string local_density_validation_error;
      try {
        for (int source_rank = 0; source_rank < world_size; ++source_rank) {
          const int begin_record = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
          const int count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
          std::optional<std::uint32_t> previous_sequence;
          for (int i = 0; i < count; ++i) {
            const std::size_t record_index = static_cast<std::size_t>(begin_record + i);
            const PmDensityPlaneRecord record = decodePmDensityPlaneRecord(
                std::span<const std::uint8_t>(exchange.recv_wire).subspan(
                    record_index * k_pm_density_plane_wire_bytes,
                    k_pm_density_plane_wire_bytes));
            if (record.origin_rank != static_cast<std::uint32_t>(source_rank) ||
                record.destination_rank != static_cast<std::uint32_t>(world_rank) ||
                record.exchange_epoch != exchange_epoch) {
              throw std::invalid_argument("PM density-plane routing metadata does not match MPI segment");
            }
            if (previous_sequence.has_value() && record.record_sequence <= *previous_sequence) {
              throw std::invalid_argument("PM density-plane sequence is non-monotonic within sender segment");
            }
            previous_sequence = record.record_sequence;
            accumulate_plane_record(record);
          }
        }
      } catch (const std::exception& error) {
        local_density_validation_error = error.what();
      } catch (...) {
        local_density_validation_error = "non-standard exception while consuming bounded density payload";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::assignDensity bounded density receive",
          world_rank,
          world_size,
          local_density_validation_error,
          routed_mpi_wait_ms);

      const std::uint64_t round_send_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      const std::uint64_t round_recv_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());
      const std::uint64_t round_send_capacity = static_cast<std::uint64_t>(exchange.send_wire.capacity());
      const std::uint64_t round_recv_capacity = static_cast<std::uint64_t>(exchange.recv_wire.capacity());
      const std::uint64_t round_combined_capacity = checkedAddBytes(
          round_send_capacity,
          round_recv_capacity,
          "PM density combined routing buffer capacity");
      const std::uint64_t round_workspace_capacity = checkedAddBytes(
          routing_metadata_capacity,
          round_combined_capacity,
          "PM density routing workspace capacity");
      if (round_workspace_capacity > k_pm_routing_workspace_target_bytes) {
        throw std::logic_error("PM density routing workspace exceeded the M1A aggregate target");
      }
      logical_records_sent += static_cast<std::uint64_t>(
          exchange.send_wire.size() / k_pm_density_plane_wire_bytes);
      peers_touched += static_cast<std::uint64_t>(std::count_if(
          exchange.send_counts.begin(), exchange.send_counts.end(), [](int count) { return count > 0; }));
      total_wire_sent += round_send_bytes;
      total_wire_received += round_recv_bytes;
      send_high_water = std::max(send_high_water, round_send_capacity);
      receive_high_water = std::max(receive_high_water, round_recv_capacity);
      combined_high_water = std::max(combined_high_water, round_combined_capacity);
      workspace_high_water = std::max(workspace_high_water, round_workspace_capacity);
      exchange.workspace_high_water_bytes = std::max(
          exchange.workspace_high_water_bytes, round_workspace_capacity);
      m_impl->recordCommunicationHighWater(
          GravityCommunicationArena::Phase::kPmDensity,
          round_workspace_capacity,
          standalone_communication_lease);
    }

    if (profile != nullptr) {
      profile->routed_density_records += logical_records_sent;
      profile->routed_density_peer_count += peers_touched;
      profile->routed_mpi_bytes_sent += total_wire_sent;
      profile->routed_mpi_bytes_received += total_wire_received;
      profile->routed_send_buffer_high_water_bytes = std::max(
          profile->routed_send_buffer_high_water_bytes, send_high_water);
      profile->routed_receive_buffer_high_water_bytes = std::max(
          profile->routed_receive_buffer_high_water_bytes, receive_high_water);
      profile->routed_combined_buffer_high_water_bytes = std::max(
          profile->routed_combined_buffer_high_water_bytes, combined_high_water);
      profile->routed_workspace_high_water_bytes = std::max(
          profile->routed_workspace_high_water_bytes, workspace_high_water);
      profile->routed_mpi_wait_ms += routed_mpi_wait_ms;
      profile->density_routing_wait_ms += routed_mpi_wait_ms;
      profile->bytes_moved += total_wire_sent + total_wire_received;
    }
#endif
  }

  const double cell_volume =
      (lengths.lx * lengths.ly * lengths.lz) / static_cast<double>(m_shape.cellCount());
  for (std::size_t local_ix = 0; local_ix < grid.slabLayout().local_nx(); ++local_ix) {
    for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
      const std::size_t row_base =
          (local_ix * m_shape.ny + iy) * density_z_stride;
      for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
        density_storage[row_base + iz] /= cell_volume;
      }
    }
  }
  if (density_plan != nullptr) {
    m_impl->publishPeriodicDensity(*density_plan, grid);
  }

  const auto stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->assign_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    profile->bytes_moved += bytesForGridSweep(m_shape.cellCount());
    profile->bytes_moved += bytesForParticles(pos_x.size());
  }
}

void PmSolver::solvePoissonPeriodic(
    PmGridStorage& grid,
    const PmSolveOptions& options,
    PmProfileEvent* profile) {
  solvePoissonPeriodicImpl(grid, options, /*materialize_potential=*/true, profile);
}

void PmSolver::solvePoissonPeriodicForcesOnly(
    PmGridStorage& grid,
    const PmSolveOptions& options,
    PmProfileEvent* profile) {
  solvePoissonPeriodicImpl(grid, options, /*materialize_potential=*/false, profile);
}

void PmSolver::solvePoissonPeriodicImpl(
    PmGridStorage& grid,
    const PmSolveOptions& options,
    bool materialize_potential,
    PmProfileEvent* profile) {
#if COSMOSIM_ENABLE_MPI
  double distributed_plan_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        materialize_potential
            ? PmCollectiveEntryKind::kSolvePoissonPeriodic
            : PmCollectiveEntryKind::kSolvePoissonPeriodicForcesOnly,
        rank_local_serial_layout,
        m_shape,
        options,
        distributed_plan_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::solvePoissonPeriodic distributed API preflight",
        mpi_world_rank,
        mpi_world_size,
        distributed_plan_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonPeriodic");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::solvePoissonPeriodic requires valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::solvePoissonPeriodic slab layout world metadata must match MPI_COMM_WORLD");
          }
        });
    if (profile != nullptr) {
      profile->routed_mpi_wait_ms += distributed_plan_preflight_mpi_wait_ms;
    }
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonPeriodic");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::solvePoissonPeriodic requires valid PM slab layout");
  }
  if (grid.slabLayout().world_size > 1) {
#if !(COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI)
    throw std::invalid_argument(
        "PmSolver::solvePoissonPeriodic distributed slabs require COSMOSIM_ENABLE_FFTW=ON and COSMOSIM_ENABLE_MPI=ON");
#endif
  } else if (!grid.ownsFullDomain()) {
    throw std::invalid_argument("PmSolver::solvePoissonPeriodic single-rank path requires full-domain slab ownership");
  }

  auto& plan = m_impl->planForLayout(grid.slabLayout(), options.decomposition_mode);
  auto real = m_impl->realGrid();
  auto fourier = m_impl->fourierGrid();

  const bool fft_backed_density_ready =
      m_impl->periodicDensityReadyFor(plan, grid);
  if (!fft_backed_density_ready &&
      m_impl->periodicDensityTokenTargetsGrid(plan, grid)) {
    throw std::logic_error(
        "Periodic FFT-backed density token is stale for the current PM grid generation");
  }
  if (!fft_backed_density_ready) {
    const std::span<const double> compact_density =
        static_cast<const PmGridStorage&>(grid).density();
    if (!grid.hasDensityStorage() || compact_density.size() != grid.localCellCount()) {
      throw std::logic_error(
          "Periodic PM solve has neither current FFT-backed density nor materialized compact density");
    }
    std::fill(real.begin(), real.end(), 0.0);
    for (std::size_t local_ix = 0; local_ix < grid.slabLayout().local_nx(); ++local_ix) {
      const std::size_t global_ix = grid.slabLayout().globalXFromLocal(local_ix);
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        const std::size_t compact_base = grid.linearIndex(global_ix, iy, 0);
        const std::size_t fft_base =
            plan.realRowBase(local_ix, iy);
        std::copy_n(
            compact_density.begin() + static_cast<std::ptrdiff_t>(compact_base),
            m_shape.nz,
            real.begin() + static_cast<std::ptrdiff_t>(fft_base));
      }
    }
  }

  // Physical reductions deliberately skip FFTW row padding. The direct-density
  // path deposited and normalized only these locations; padding remains scratch.
  double local_density_sum = 0.0;
  for (std::size_t local_ix = 0; local_ix < grid.slabLayout().local_nx(); ++local_ix) {
    for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
      const std::size_t fft_base =
          plan.realRowBase(local_ix, iy);
      for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
        local_density_sum += real[fft_base + iz];
      }
    }
  }
  double global_density_sum = local_density_sum;
#if COSMOSIM_ENABLE_MPI
  if (grid.slabLayout().world_size > 1) {
    MPI_Allreduce(&local_density_sum, &global_density_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  const double mean_density =
      global_density_sum / static_cast<double>(m_shape.cellCount());
#if COSMOSIM_HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t local_ix = 0; local_ix < grid.slabLayout().local_nx(); ++local_ix) {
    for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
      const std::size_t fft_base =
          plan.realRowBase(local_ix, iy);
      for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
        real[fft_base + iz] -= mean_density;
      }
    }
  }

  // FFT execution destroys authoritative real-space density contents. Clear
  // the generation/identity token before handing the allocation to FFTW.
  m_impl->consumePeriodicDensity(plan);
  const double forward_fft_ms = m_impl->forwardFft();
  if (profile != nullptr) {
    profile->fft_forward_ms += forward_fft_ms;
    if (plan.spectral_transposed) {
      profile->fft_transpose_ms += forward_fft_ms;
    }
  }

  const auto poisson_start = std::chrono::steady_clock::now();
  const BoxLengths lengths = effectiveBoxLengths(options);
  const bool spectral_operator_rebuilt =
      m_impl->applyDerivedPoissonOperator(plan, lengths, options, m_shape);
  if (profile != nullptr && spectral_operator_rebuilt) {
    ++profile->spectral_operator_rebuilds;
  }
  const auto poisson_stop = std::chrono::steady_clock::now();

  if (profile != nullptr) {
    profile->poisson_ms += std::chrono::duration<double, std::milli>(poisson_stop - poisson_start).count();
  }

  if (options.force_gradient_backend == PmForceGradientBackend::kFiniteDifference4) {
    // One inverse FFT of phi_k into the plan real array, then a fourth-order
    // centered difference F = -grad(phi) on the periodic mesh. The phi_k copy
    // (potential_k) is never needed. See PmForceGradientBackend.
    const auto grad_start = std::chrono::steady_clock::now();
    const double fft_time = m_impl->inverseFft();
    if (profile != nullptr) {
      profile->fft_inverse_ms += fft_time;
      if (plan.spectral_transposed) profile->fft_transpose_ms += fft_time;
    }
    const auto phi = std::span<const double>(m_impl->realGrid());
#if COSMOSIM_ENABLE_FFTW
    const double normalization = 1.0 / static_cast<double>(m_shape.cellCount());
#else
    const double normalization = 1.0;  // the naive inverse DFT is normalized
#endif
    const std::size_t nx = m_shape.nx, ny = m_shape.ny, nz = m_shape.nz;
    const std::size_t local_nx = grid.slabLayout().local_nx();
    const std::size_t stride = plan.real_z_stride;
    const std::size_t plane_values = ny * nz;
    const bool full_domain = !plan.is_distributed;
    std::span<double> halo;
    if (!full_domain) {
#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
      // Every rank must own two x planes so the +-2 stencil reaches only the
      // immediate slab neighbours. The check is collective and uniform.
      unsigned long long local_planes = static_cast<unsigned long long>(local_nx);
      unsigned long long minimum_planes = 0ULL;
      requirePmMpiSuccess(MPI_Allreduce(&local_planes, &minimum_planes, 1, MPI_UNSIGNED_LONG_LONG,
                                        MPI_MIN, MPI_COMM_WORLD),
                          "PM finite-difference slab-depth consensus");
      if (minimum_planes < 2ULL) {
        throw std::invalid_argument(
            "PM finite-difference force gradient requires every MPI rank to own at least two x planes");
      }
      if (plane_values > static_cast<std::size_t>(std::numeric_limits<int>::max() / 2)) {
        throw std::overflow_error("PM finite-difference halo exchange exceeds MPI int count");
      }
      halo = m_impl->finiteDifferenceHaloPlanes(plane_values);
      const std::span<double> send_planes = halo.subspan(4U * plane_values, 2U * plane_values);
      const int world_size = grid.slabLayout().world_size;
      const int world_rank = grid.slabLayout().world_rank;
      const int left = (world_rank + world_size - 1) % world_size;
      const int right = (world_rank + 1) % world_size;
      const auto pack = [&](std::size_t first_local_plane) {
        for (std::size_t p = 0; p < 2U; ++p) {
          for (std::size_t iy = 0; iy < ny; ++iy) {
            const std::size_t src = ((first_local_plane + p) * ny + iy) * stride;
            std::copy_n(phi.begin() + static_cast<std::ptrdiff_t>(src), nz,
                        send_planes.begin() + static_cast<std::ptrdiff_t>(p * plane_values + iy * nz));
          }
        }
      };
      const int count = static_cast<int>(2U * plane_values);
      // Last two planes go right and arrive as the right rank's left halo.
      pack(local_nx - 2U);
      requirePmMpiSuccess(MPI_Sendrecv(send_planes.data(), count, MPI_DOUBLE, right, 7101,
                                       halo.data(), count, MPI_DOUBLE, left, 7101,
                                       MPI_COMM_WORLD, MPI_STATUS_IGNORE),
                          "PM finite-difference left halo exchange");
      // First two planes go left and arrive as the left rank's right halo.
      pack(0U);
      requirePmMpiSuccess(MPI_Sendrecv(send_planes.data(), count, MPI_DOUBLE, left, 7102,
                                       halo.data() + 2U * plane_values, count, MPI_DOUBLE, right, 7102,
                                       MPI_COMM_WORLD, MPI_STATUS_IGNORE),
                          "PM finite-difference right halo exchange");
      if (profile != nullptr) {
        profile->bytes_moved += 4U * plane_values * sizeof(double);
      }
#else
      throw std::logic_error("distributed PM finite-difference gradient requires FFTW-MPI");
#endif
    }
    // Plane pointer for local x offset t in [-2, local_nx + 1]; rows have
    // stride `stride` inside the FFT array and nz inside the halo planes.
    const auto plane_row = [&](std::ptrdiff_t t, std::size_t iy) -> const double* {
      const auto n_local = static_cast<std::ptrdiff_t>(local_nx);
      if (t >= 0 && t < n_local) {
        return phi.data() + (static_cast<std::size_t>(t) * ny + iy) * stride;
      }
      if (full_domain) {
        const auto wrapped = static_cast<std::size_t>((t % n_local + n_local) % n_local);
        return phi.data() + (wrapped * ny + iy) * stride;
      }
      const std::size_t halo_plane = t < 0 ? static_cast<std::size_t>(t + 2) : static_cast<std::size_t>(t - n_local) + 2U;
      return halo.data() + halo_plane * plane_values + iy * nz;
    };
    const BoxLengths fd_lengths = effectiveBoxLengths(options);
    const double scale_x = -normalization * static_cast<double>(nx) / (12.0 * fd_lengths.lx);
    const double scale_y = -normalization * static_cast<double>(ny) / (12.0 * fd_lengths.ly);
    const double scale_z = -normalization * static_cast<double>(nz) / (12.0 * fd_lengths.lz);
    if (materialize_potential) grid.ensurePotentialStorage();
    const std::span<double> potential_lane = materialize_potential ? grid.potential() : std::span<double>{};
    const std::array<std::span<double>, 3U> force_lanes{grid.force_x(), grid.force_y(), grid.force_z()};
#if COSMOSIM_HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t local_ix = 0; local_ix < local_nx; ++local_ix) {
      const std::size_t global_ix = grid.slabLayout().globalXFromLocal(local_ix);
      const auto lx = static_cast<std::ptrdiff_t>(local_ix);
      for (std::size_t iy = 0; iy < ny; ++iy) {
        const double* row = plane_row(lx, iy);
        const double* xm2 = plane_row(lx - 2, iy);
        const double* xm1 = plane_row(lx - 1, iy);
        const double* xp1 = plane_row(lx + 1, iy);
        const double* xp2 = plane_row(lx + 2, iy);
        const double* ym2 = plane_row(lx, (iy + 2U * ny - 2U) % ny);
        const double* ym1 = plane_row(lx, (iy + ny - 1U) % ny);
        const double* yp1 = plane_row(lx, (iy + 1U) % ny);
        const double* yp2 = plane_row(lx, (iy + 2U) % ny);
        const std::size_t compact_base = grid.linearIndex(global_ix, iy, 0);
        for (std::size_t iz = 0; iz < nz; ++iz) {
          const std::size_t zm2 = (iz + 2U * nz - 2U) % nz, zm1 = (iz + nz - 1U) % nz;
          const std::size_t zp1 = (iz + 1U) % nz, zp2 = (iz + 2U) % nz;
          force_lanes[0][compact_base + iz] =
              scale_x * (8.0 * (xp1[iz] - xm1[iz]) - (xp2[iz] - xm2[iz]));
          force_lanes[1][compact_base + iz] =
              scale_y * (8.0 * (yp1[iz] - ym1[iz]) - (yp2[iz] - ym2[iz]));
          force_lanes[2][compact_base + iz] =
              scale_z * (8.0 * (row[zp1] - row[zm1]) - (row[zp2] - row[zm2]));
          if (!potential_lane.empty()) potential_lane[compact_base + iz] = row[iz] * normalization;
        }
      }
    }
    if (profile != nullptr) {
      profile->gradient_ms += std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - grad_start).count();
      profile->bytes_moved += bytesForGridSweep(m_shape.cellCount()) * (materialize_potential ? 5U : 4U);
    }
    return;
  }

  auto potential_k = m_impl->potentialScratch();
  std::copy(fourier.begin(), fourier.end(), potential_k.begin());

  auto inverse_current_spectrum_into = [this, profile, &plan, &grid](std::span<double> dst, int axis = -1) {
    const double fft_time = m_impl->inverseFft();
    const auto normalization_start = std::chrono::steady_clock::now();
    auto real_values = m_impl->realGrid();
#if COSMOSIM_ENABLE_FFTW
    const double normalization = 1.0 / static_cast<double>(m_shape.cellCount());
#else
    const double normalization = 1.0;  // the naive inverse DFT is normalized
#endif
    // Padded in-place rows (real_z_stride) -> compact force/potential lanes,
    // identical for serial and distributed slabs.
#if COSMOSIM_HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t local_ix = 0; local_ix < plan.layout.local_nx(); ++local_ix) {
      const std::size_t global_ix = plan.layout.globalXFromLocal(local_ix);
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        const std::size_t compact_base = grid.linearIndex(global_ix, iy, 0);
        const std::size_t fftw_base = plan.realRowBase(local_ix, iy);
        for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
          dst[compact_base + iz] = real_values[fftw_base + iz] * normalization;
        }
      }
    }
    if (profile != nullptr) {
      if (axis >= 0) {
        profile->inverse_axis_ms[static_cast<std::size_t>(axis)] += fft_time;
        profile->normalization_axis_ms[static_cast<std::size_t>(axis)] +=
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - normalization_start).count();
      }
      profile->fft_inverse_ms += fft_time;
      if (plan.spectral_transposed) {
        profile->fft_transpose_ms += fft_time;
      }
    }
  };

  const auto grad_start = std::chrono::steady_clock::now();
  if (materialize_potential) {
    grid.ensurePotentialStorage();
    std::copy(potential_k.begin(), potential_k.end(), fourier.begin());
    inverse_current_spectrum_into(grid.potential());
  }

  const std::array<std::span<double>, 3U> force_lanes{
      grid.force_x(), grid.force_y(), grid.force_z()};
  for (int axis = 0; axis < 3; ++axis) {
    const auto factor_start = std::chrono::steady_clock::now();
    m_impl->fillGradientSpectrum(plan, lengths, m_shape, axis);
    if (profile != nullptr) {
      profile->gradient_factor_ms[static_cast<std::size_t>(axis)] +=
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - factor_start).count();
    }
    inverse_current_spectrum_into(force_lanes[static_cast<std::size_t>(axis)], axis);
  }
  const auto grad_stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->gradient_ms += std::chrono::duration<double, std::milli>(grad_stop - grad_start).count();
  }

  if (profile != nullptr) {
    profile->bytes_moved += bytesForGridSweep(m_shape.cellCount()) *
        (materialize_potential ? 6U : 5U);
    if (plan.spectral_transposed) {
      profile->fft_transpose_bytes +=
          static_cast<std::uint64_t>(sizeof(std::complex<double>)) *
          static_cast<std::uint64_t>(fourier.size()) *
          (materialize_potential ? 8ULL : 6ULL);
    }
  }
}

void PmSolver::solvePoissonIsolatedOpen(PmGridStorage& grid, const PmSolveOptions& options, PmProfileEvent* profile) {
#if COSMOSIM_ENABLE_MPI
  double isolated_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kSolvePoissonIsolatedOpen,
        rank_local_serial_layout,
        m_shape,
        options,
        isolated_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::solvePoissonIsolatedOpen API preflight",
        mpi_world_rank,
        mpi_world_size,
        isolated_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny ||
              grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonIsolatedOpen");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument(
                "PmSolver::solvePoissonIsolatedOpen requires valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::solvePoissonIsolatedOpen slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (m_shape.nx < 3 || m_shape.ny < 3 || m_shape.nz < 3) {
            throw std::invalid_argument(
                "PmSolver::solvePoissonIsolatedOpen requires nx,ny,nz >= 3 for one-sided boundary gradients");
          }
          const bool distributed_layout = !rank_local_serial_layout;
          const std::uint64_t workspace_estimate =
              estimateIsolatedRootWorkspaceBytes(m_shape, distributed_layout);
          if (options.isolated_open_root_workspace_limit_bytes == 0U ||
              workspace_estimate > options.isolated_open_root_workspace_limit_bytes) {
            throw std::runtime_error(isolatedPmGuardMessage(
                m_shape,
                distributed_layout ? grid.slabLayout().world_size : 1,
                distributed_layout ? grid.slabLayout().world_rank : 0,
                workspace_estimate,
                options.isolated_open_root_workspace_limit_bytes));
          }
          if (grid.localCellCount() >
              static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::overflow_error("isolated PM local slab count exceeds MPI int limit");
          }
          if (distributed_layout) {
            std::size_t displacement = 0U;
            for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
              const auto owned = parallel::pmOwnedXRangeForRank(
                  grid.slabLayout().global_nx, grid.slabLayout().world_size, rank);
              const std::size_t count = checkedProduct(
                  checkedProduct(
                      owned.extentX(),
                      grid.slabLayout().global_ny,
                      "isolated PM preflight slab count"),
                  grid.slabLayout().global_nz,
                  "isolated PM preflight slab count");
              if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
                  displacement > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::overflow_error(
                    "isolated PM gather/scatter count or displacement exceeds MPI int limit");
              }
              if (count > std::numeric_limits<std::size_t>::max() - displacement) {
                throw std::overflow_error(
                    "isolated PM preflight gather/scatter displacement overflows size_t");
              }
              displacement += count;
            }
          }
        });
    if (profile != nullptr) {
      profile->routed_mpi_wait_ms += isolated_preflight_mpi_wait_ms;
    }
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonIsolatedOpen");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen requires valid PM slab layout");
  }
#if !COSMOSIM_ENABLE_MPI
  if (grid.slabLayout().world_size > 1) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
#endif

  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else
#endif
  if (!grid.ownsFullDomain()) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen single-rank path requires full-domain ownership");
  }

  const BoxLengths lengths = effectiveBoxLengths(options);
  const double dx = lengths.lx / static_cast<double>(m_shape.nx);
  const double dy = lengths.ly / static_cast<double>(m_shape.ny);
  const double dz = lengths.lz / static_cast<double>(m_shape.nz);
  if (m_shape.nx < 3 || m_shape.ny < 3 || m_shape.nz < 3) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen requires nx,ny,nz >= 3 for one-sided boundary gradients");
  }

  constexpr int k_root = 0;
  const bool is_root = !distributed_slabs || grid.slabLayout().world_rank == k_root;
  const auto poisson_start = std::chrono::steady_clock::now();

  const int guard_world_size = distributed_slabs ? grid.slabLayout().world_size : 1;
  const int guard_world_rank = distributed_slabs ? grid.slabLayout().world_rank : 0;
  const std::uint64_t isolated_root_workspace_estimate_bytes =
      estimateIsolatedRootWorkspaceBytes(m_shape, distributed_slabs);
  if (options.isolated_open_root_workspace_limit_bytes == 0U ||
      isolated_root_workspace_estimate_bytes > options.isolated_open_root_workspace_limit_bytes) {
    throw std::runtime_error(isolatedPmGuardMessage(
        m_shape,
        guard_world_size,
        guard_world_rank,
        isolated_root_workspace_estimate_bytes,
        options.isolated_open_root_workspace_limit_bytes));
  }

  std::vector<double> global_density;
#if COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    const std::size_t global_size = checkedProduct(
        checkedProduct(grid.slabLayout().global_nx, grid.slabLayout().global_ny, "isolated PM distributed gather"),
        grid.slabLayout().global_nz,
        "isolated PM distributed gather");
    std::vector<int> recv_counts(static_cast<std::size_t>(grid.slabLayout().world_size), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(grid.slabLayout().world_size), 0);
    std::size_t disp = 0;
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const auto owned = parallel::pmOwnedXRangeForRank(grid.slabLayout().global_nx, grid.slabLayout().world_size, rank);
      const std::size_t count = owned.extentX() * grid.slabLayout().global_ny * grid.slabLayout().global_nz;
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("isolated PM distributed slab count exceeds MPI int limit");
      }
      if (disp > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("isolated PM distributed slab displacement exceeds MPI int limit");
      }
      recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
      recv_displs[static_cast<std::size_t>(rank)] = static_cast<int>(disp);
      disp += count;
    }
    if (is_root) {
      global_density.assign(global_size, 0.0);
    }
    MPI_Gatherv(
        grid.density().data(),
        static_cast<int>(grid.density().size()),
        MPI_DOUBLE,
        is_root ? global_density.data() : nullptr,
        recv_counts.data(),
        recv_displs.data(),
        MPI_DOUBLE,
        k_root,
        MPI_COMM_WORLD);
  }
#endif

  const std::size_t pad_nx = 2U * m_shape.nx;
  const std::size_t pad_ny = 2U * m_shape.ny;
  const std::size_t pad_nz = 2U * m_shape.nz;
  const std::size_t pad_total =
      checkedProduct(checkedProduct(pad_nx, pad_ny, "isolated PM padded workspace"),
                     pad_nz,
                     "isolated PM padded workspace");
  std::vector<double> full_potential;
  std::vector<double> full_force_x;
  std::vector<double> full_force_y;
  std::vector<double> full_force_z;

  if (is_root) {
    m_impl->ensureIsolatedWorkspace(pad_nx, pad_ny, pad_nz);
    auto& ws = m_impl->isolatedWorkspace();
    std::fill(ws.rho_k.begin(), ws.rho_k.end(), std::complex<double>(0.0, 0.0));
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
          const std::size_t src = distributed_slabs ? (ix * m_shape.ny + iy) * m_shape.nz + iz : grid.linearIndex(ix, iy, iz);
          const std::size_t dst = (ix * pad_ny + iy) * pad_nz + iz;
          ws.rho_k[dst] = {(distributed_slabs ? global_density[src] : grid.density()[src]), 0.0};
        }
      }
    }
    m_impl->isolatedForward(ws.rho_k);

    if (!ws.kernel_ready || ws.dx != dx || ws.dy != dy || ws.dz != dz || ws.split_scale != options.tree_pm_split_scale_comoving) {
      ws.dx = dx;
      ws.dy = dy;
      ws.dz = dz;
      ws.split_scale = options.tree_pm_split_scale_comoving;
      std::fill(ws.kernel_k.begin(), ws.kernel_k.end(), std::complex<double>(0.0, 0.0));
      for (std::size_t ix = 0; ix < pad_nx; ++ix) {
        const std::ptrdiff_t sx = ix <= m_shape.nx ? static_cast<std::ptrdiff_t>(ix)
                                                   : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(pad_nx);
        const double rx = static_cast<double>(sx) * dx;
        for (std::size_t iy = 0; iy < pad_ny; ++iy) {
          const std::ptrdiff_t sy = iy <= m_shape.ny ? static_cast<std::ptrdiff_t>(iy)
                                                     : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(pad_ny);
          const double ry = static_cast<double>(sy) * dy;
          for (std::size_t iz = 0; iz < pad_nz; ++iz) {
            const std::ptrdiff_t sz = iz <= m_shape.nz ? static_cast<std::ptrdiff_t>(iz)
                                                       : static_cast<std::ptrdiff_t>(iz) - static_cast<std::ptrdiff_t>(pad_nz);
            const double rz = static_cast<double>(sz) * dz;
            const double r = std::sqrt(rx * rx + ry * ry + rz * rz);
            double kernel = 0.0;
            if (r > 0.0) {
              kernel = -1.0 / r;
              if (options.tree_pm_split_scale_comoving > 0.0) {
                kernel = -std::erf(0.5 * r / options.tree_pm_split_scale_comoving) / r;
              }
            } else if (options.tree_pm_split_scale_comoving > 0.0) {
              // lim_{r->0} -erf(r/(2*r_s))/r = -1/(sqrt(pi)*r_s).
              // Leaving the origin at zero instead reverses the centered
              // force in the cell immediately adjacent to a point source.
              kernel = -1.0 / (std::sqrt(std::acos(-1.0)) * options.tree_pm_split_scale_comoving);
            } else {
              // The unsplit 1/r Green function has a finite *cell average*.
              // Approximate the integral over the source voxel by symmetric
              // midpoint cubature (even subdivision count avoids r=0).
              // This preserves attraction near a cell-centered point source
              // without introducing an arbitrary particle softening length.
              constexpr int k_subdivisions = 16;
              double inverse_radius_sum = 0.0;
              for (int qx = 0; qx < k_subdivisions; ++qx) {
                const double cx = (static_cast<double>(qx) + 0.5) / k_subdivisions - 0.5;
                for (int qy = 0; qy < k_subdivisions; ++qy) {
                  const double cy = (static_cast<double>(qy) + 0.5) / k_subdivisions - 0.5;
                  for (int qz = 0; qz < k_subdivisions; ++qz) {
                    const double cz = (static_cast<double>(qz) + 0.5) / k_subdivisions - 0.5;
                    inverse_radius_sum += 1.0 / std::sqrt(
                        (cx * dx) * (cx * dx) +
                        (cy * dy) * (cy * dy) +
                        (cz * dz) * (cz * dz));
                  }
                }
              }
              kernel = -inverse_radius_sum / static_cast<double>(
                  k_subdivisions * k_subdivisions * k_subdivisions);
            }
            ws.kernel_k[(ix * pad_ny + iy) * pad_nz + iz] = {kernel, 0.0};
          }
        }
      }
      m_impl->isolatedForward(ws.kernel_k);
      ws.kernel_ready = true;
    }

    for (std::size_t i = 0; i < pad_total; ++i) {
      ws.rho_k[i] *= ws.kernel_k[i];
    }
    m_impl->isolatedInverse(ws.rho_k);

    const double cell_volume = dx * dy * dz;
    const double prefactor = options.gravitational_constant_code * cell_volume;
    full_potential.assign(m_shape.cellCount(), 0.0);
    full_force_x.assign(m_shape.cellCount(), 0.0);
    full_force_y.assign(m_shape.cellCount(), 0.0);
    full_force_z.assign(m_shape.cellCount(), 0.0);
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
          const std::size_t physical_index = (ix * m_shape.ny + iy) * m_shape.nz + iz;
          const std::size_t padded_index = (ix * pad_ny + iy) * pad_nz + iz;
          double value = ws.rho_k[padded_index].real();
#if COSMOSIM_ENABLE_FFTW
          value /= static_cast<double>(pad_total);
#endif
          full_potential[physical_index] = prefactor * value;
        }
      }
    }

    const auto fieldPotential = [&](std::size_t ix, std::size_t iy, std::size_t iz) -> double {
      return full_potential[(ix * m_shape.ny + iy) * m_shape.nz + iz];
    };
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
          const std::size_t center = (ix * m_shape.ny + iy) * m_shape.nz + iz;
          const auto grad_x = [&]() {
            if (ix == 0) return (-3.0 * fieldPotential(0, iy, iz) + 4.0 * fieldPotential(1, iy, iz) - fieldPotential(2, iy, iz)) / (2.0 * dx);
            if (ix + 1U == m_shape.nx) return (3.0 * fieldPotential(m_shape.nx - 1U, iy, iz) - 4.0 * fieldPotential(m_shape.nx - 2U, iy, iz) + fieldPotential(m_shape.nx - 3U, iy, iz)) / (2.0 * dx);
            return (fieldPotential(ix + 1U, iy, iz) - fieldPotential(ix - 1U, iy, iz)) / (2.0 * dx);
          }();
          const auto grad_y = [&]() {
            if (iy == 0) return (-3.0 * fieldPotential(ix, 0, iz) + 4.0 * fieldPotential(ix, 1, iz) - fieldPotential(ix, 2, iz)) / (2.0 * dy);
            if (iy + 1U == m_shape.ny) return (3.0 * fieldPotential(ix, m_shape.ny - 1U, iz) - 4.0 * fieldPotential(ix, m_shape.ny - 2U, iz) + fieldPotential(ix, m_shape.ny - 3U, iz)) / (2.0 * dy);
            return (fieldPotential(ix, iy + 1U, iz) - fieldPotential(ix, iy - 1U, iz)) / (2.0 * dy);
          }();
          const auto grad_z = [&]() {
            if (iz == 0) return (-3.0 * fieldPotential(ix, iy, 0) + 4.0 * fieldPotential(ix, iy, 1) - fieldPotential(ix, iy, 2)) / (2.0 * dz);
            if (iz + 1U == m_shape.nz) return (3.0 * fieldPotential(ix, iy, m_shape.nz - 1U) - 4.0 * fieldPotential(ix, iy, m_shape.nz - 2U) + fieldPotential(ix, iy, m_shape.nz - 3U)) / (2.0 * dz);
            return (fieldPotential(ix, iy, iz + 1U) - fieldPotential(ix, iy, iz - 1U)) / (2.0 * dz);
          }();
          full_force_x[center] = -grad_x;
          full_force_y[center] = -grad_y;
          full_force_z[center] = -grad_z;
        }
      }
    }
  }

#if COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    std::vector<int> recv_counts(static_cast<std::size_t>(grid.slabLayout().world_size), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(grid.slabLayout().world_size), 0);
    std::size_t disp = 0;
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const auto owned = parallel::pmOwnedXRangeForRank(grid.slabLayout().global_nx, grid.slabLayout().world_size, rank);
      const std::size_t count = owned.extentX() * grid.slabLayout().global_ny * grid.slabLayout().global_nz;
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("isolated PM scatter slab count exceeds MPI int limit");
      }
      if (disp > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("isolated PM scatter slab displacement exceeds MPI int limit");
      }
      recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
      recv_displs[static_cast<std::size_t>(rank)] = static_cast<int>(disp);
      disp += count;
    }
    MPI_Scatterv(is_root ? full_potential.data() : nullptr, recv_counts.data(), recv_displs.data(), MPI_DOUBLE, grid.potential().data(), static_cast<int>(grid.potential().size()), MPI_DOUBLE, k_root, MPI_COMM_WORLD);
    MPI_Scatterv(is_root ? full_force_x.data() : nullptr, recv_counts.data(), recv_displs.data(), MPI_DOUBLE, grid.force_x().data(), static_cast<int>(grid.force_x().size()), MPI_DOUBLE, k_root, MPI_COMM_WORLD);
    MPI_Scatterv(is_root ? full_force_y.data() : nullptr, recv_counts.data(), recv_displs.data(), MPI_DOUBLE, grid.force_y().data(), static_cast<int>(grid.force_y().size()), MPI_DOUBLE, k_root, MPI_COMM_WORLD);
    MPI_Scatterv(is_root ? full_force_z.data() : nullptr, recv_counts.data(), recv_displs.data(), MPI_DOUBLE, grid.force_z().data(), static_cast<int>(grid.force_z().size()), MPI_DOUBLE, k_root, MPI_COMM_WORLD);
  } else
#endif
  {
    std::copy(full_potential.begin(), full_potential.end(), grid.potential().begin());
    std::copy(full_force_x.begin(), full_force_x.end(), grid.force_x().begin());
    std::copy(full_force_y.begin(), full_force_y.end(), grid.force_y().begin());
    std::copy(full_force_z.begin(), full_force_z.end(), grid.force_z().begin());
  }

  if (profile != nullptr) {
    const auto stop = std::chrono::steady_clock::now();
    profile->poisson_ms += std::chrono::duration<double, std::milli>(stop - poisson_start).count();
    profile->isolated_open_root_workspace_estimate_bytes =
        std::max(profile->isolated_open_root_workspace_estimate_bytes, isolated_root_workspace_estimate_bytes);
    profile->isolated_open_root_workspace_limit_bytes = options.isolated_open_root_workspace_limit_bytes;
    if (distributed_slabs) {
      const std::uint64_t local_gather_scatter_bytes =
          checkedBytesForCount(
              checkedProduct(5U, grid.localCellCount(), "isolated PM gather/scatter profile"),
              sizeof(double),
              "isolated PM gather/scatter profile");
      profile->isolated_open_gather_bytes += local_gather_scatter_bytes;
      profile->bytes_moved += local_gather_scatter_bytes;
    }
  }
}

namespace {

// Shared collective body of PmSolver::plan*Routing() — a rank-consistent
// planning transaction (2026-10-26 repair of the round-2 planner, which
// validated span extents and box lengths *before* the collective so one rank
// could throw while its peers entered MPI_Allreduce).
//
// Contract (distributed slabs, MPI_COMM_WORLD):
//   1. Every rank executes validate_local(), key construction and the local
//      workload measurement inside one try block. No rank-local exception can
//      escape before step 2.
//   2. Every rank executes exactly one MPI_Allreduce(MAX) over
//      [count, M_local, failure_kind, world_size - failing_rank].
//   3. If any rank failed, every rank throws: the failing rank rethrows its
//      original exception unchanged; peers throw a runtime_error naming the
//      lowest failing rank and the failure class. No cache is published and no
//      rank proceeds to routing data exchange.
// failure_kind: 1 = recoverable input/contract validation failure, 2 =
// allocation failure (std::bad_alloc; reported collectively on a best-effort
// basis — a rank that cannot allocate the collective itself cannot recover).
// MPI errors raised by the reduction itself are not recoverable: they throw
// from requirePmMpiSuccess with no further collective.
// Serial meshes (one rank, or a rank owning the full domain) run the same
// validation locally and perform no collective.
template <typename ValidateLocal, typename GridX>
PmSolver::PmRoutingWorkloadPlan planPmRoutingCollective(
    PmRoutingPhaseKind kind,
    const PmGridShape& shape,
    const PmGridStorage& grid,
    const PmSolveOptions& options,
    std::size_t count,
    const ValidateLocal& validate_local,
    const GridX& grid_x,
    PmRoutingPlanKey& key,
    std::uint64_t& global_max_count,
    std::uint64_t& global_max_peer_records,
    std::uint64_t& local_max_round_records) {
  PmSolver::PmRoutingWorkloadPlan plan;
#if COSMOSIM_ENABLE_MPI
  const int world_size = grid.slabLayout().world_size;
  const int world_rank = grid.slabLayout().world_rank;
  if (world_size <= 1 || grid.ownsFullDomain()) {
    validate_local();
    return plan;
  }
  std::exception_ptr local_failure;
  std::uint64_t failure_kind = 0U;
  PmRoutingCapacityModel model{};
  PmRemotePlaneWorkload local{};
  try {
    validate_local();
    validateOptions(shape, options);
    key = makePmRoutingPlanKey(shape, grid.slabLayout(), options);
    model = pmPlannedRoutingModel(kind, world_size, options.routing_exchange_batch_bytes);
    const std::size_t particles_per_round = pmRoutingParticlesPerRound(
        model.effective_per_peer_payload_bytes, pmRoutingWireRecordBytes(kind));
    const PmDecompositionView decomposition_view(shape, grid.slabLayout(), options.decomposition_mode);
    local = measurePmRemotePlaneWorkload(
        count, particles_per_round, world_rank, world_size,
        options.boundary_condition == PmBoundaryCondition::kPeriodic,
        options.assignment_scheme, shape, decomposition_view, grid_x);
  } catch (const std::bad_alloc&) {
    local_failure = std::current_exception();
    failure_kind = 2U;
  } catch (...) {
    local_failure = std::current_exception();
    failure_kind = 1U;
  }
  const std::uint64_t failing_rank_code = local_failure
      ? static_cast<std::uint64_t>(world_size - world_rank) : 0U;
  std::array<std::uint64_t, 4U> local_maxima{
      static_cast<std::uint64_t>(count), local.max_peer_records, failure_kind, failing_rank_code};
  std::array<std::uint64_t, 4U> global_maxima{};
  requirePmMpiSuccess(
      MPI_Allreduce(local_maxima.data(), global_maxima.data(), 4, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD),
      "PmSolver routing workload plan MPI_Allreduce");
  if (global_maxima[2] != 0U) {
    if (local_failure) std::rethrow_exception(local_failure);
    const std::uint64_t first_rank = static_cast<std::uint64_t>(world_size) - global_maxima[3];
    throw std::runtime_error(
        std::string("PM ") + (kind == PmRoutingPhaseKind::kDensity ? "density" : "force-interpolation") +
        " routing workload planning failed on peer rank " + std::to_string(first_rank) +
        (global_maxima[2] >= 2U ? " (allocation failure)" : " (input validation failure)") +
        "; no rank entered routing data exchange");
  }
  global_max_count = global_maxima[0];
  global_max_peer_records = global_maxima[1];
  local_max_round_records = local.max_round_records;
  const PmRoutingReservation reservation = pmRoutingReservation(
      kind, model, world_size, local.max_round_records, global_max_peer_records);
  plan.distributed = true;
  plan.phase_arena_bytes = reservation.phase_arena_bytes;
  plan.send_reservation_bytes = reservation.send_bytes;
  plan.receive_reservation_bytes = reservation.receive_bytes;
  plan.configured_ceiling_bytes = checkedAddBytes(
      model.max_send_payload_bytes, model.max_receive_payload_bytes, "PM routing configured ceiling");
  plan.global_max_peer_records = global_max_peer_records;
#else
  static_cast<void>(kind); static_cast<void>(shape); static_cast<void>(grid);
  static_cast<void>(options); static_cast<void>(count); static_cast<void>(grid_x);
  static_cast<void>(key); static_cast<void>(global_max_count);
  static_cast<void>(global_max_peer_records); static_cast<void>(local_max_round_records);
  validate_local();
#endif
  return plan;
}

}  // namespace

PmSolver::PmRoutingWorkloadPlan PmSolver::planDensityRouting(
    const PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const PmSolveOptions& options) const {
  // Lengths are resolved inside the planning transaction (validate_local) so
  // an invalid box cannot make this rank leave before its peers' collective.
  BoxLengths lengths{};
  double inv_dx = 0.0;
  const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
  const auto validate_local = [&]() {
    if (pos_y.size() != pos_x.size() || pos_z.size() != pos_x.size() || mass.size() != pos_x.size()) {
      throw std::invalid_argument("PmSolver::planDensityRouting source extents must match");
    }
    lengths = effectiveBoxLengths(options);
    if (!(std::isfinite(lengths.lx) && lengths.lx > 0.0 && std::isfinite(lengths.ly) && lengths.ly > 0.0 &&
          std::isfinite(lengths.lz) && lengths.lz > 0.0)) {
      throw std::invalid_argument("PmSolver::planDensityRouting requires finite positive box lengths");
    }
    inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  };
  const auto grid_x = [&](std::size_t p) -> std::optional<double> {
    // Non-finite sources are rejected collectively by assignDensity itself.
    if (!std::isfinite(pos_x[p]) || !std::isfinite(pos_y[p]) || !std::isfinite(pos_z[p])) {
      return std::nullopt;
    }
    if (!periodic &&
        (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
         !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
         !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
      return std::nullopt;
    }
    return (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
  };
  // Any previous plan of this kind is withdrawn before the transaction, so a
  // failed planning call never leaves an older plan behind.
  auto& cache = m_impl->m_density_routing_plan;
  cache = {};
  PmRoutingPlanKey key{};
  std::uint64_t global_max_count = 0U, global_max_peer_records = 0U, local_max_round_records = 0U;
  const PmRoutingWorkloadPlan plan = planPmRoutingCollective(
      PmRoutingPhaseKind::kDensity, m_shape, grid, options, pos_x.size(), validate_local, grid_x,
      key, global_max_count, global_max_peer_records, local_max_round_records);
  if (plan.distributed) {
    cache = Impl::RoutingWorkloadPlanCache{
        .valid = true,
        .key = key,
        .global_max_count = global_max_count,
        .global_max_peer_records = global_max_peer_records,
        .local_max_round_records = local_max_round_records,
    };
  }
  return plan;
}

PmSolver::PmRoutingWorkloadPlan PmSolver::planForceInterpolationRouting(
    const PmGridStorage& grid,
    const PmForceTargetView& target_view,
    const PmSolveOptions& options) const {
  const bool indexed = target_view.coordinate_layout == PmForceCoordinateLayout::kIndexedSource;
  const std::size_t count = indexed ? target_view.coordinate_source_index.size()
                                    : target_view.pos_x_comoving.size();
  const auto pos_x = target_view.pos_x_comoving;
  const auto pos_y = target_view.pos_y_comoving;
  const auto pos_z = target_view.pos_z_comoving;
  BoxLengths lengths{};
  double inv_dx = 0.0;
  const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
  const auto validate_local = [&]() {
    if (pos_y.size() != pos_x.size() || pos_z.size() != pos_x.size()) {
      throw std::invalid_argument("PmSolver::planForceInterpolationRouting coordinate extents must match");
    }
    if (indexed) {
      for (const auto row : target_view.coordinate_source_index) {
        if (static_cast<std::size_t>(row) >= pos_x.size()) {
          throw std::out_of_range(
              "PmSolver::planForceInterpolationRouting coordinate_source_index row " +
              std::to_string(static_cast<std::size_t>(row)) + " exceeds coordinate extent " +
              std::to_string(pos_x.size()));
        }
      }
    }
    lengths = effectiveBoxLengths(options);
    if (!(std::isfinite(lengths.lx) && lengths.lx > 0.0 && std::isfinite(lengths.ly) && lengths.ly > 0.0 &&
          std::isfinite(lengths.lz) && lengths.lz > 0.0)) {
      throw std::invalid_argument("PmSolver::planForceInterpolationRouting requires finite positive box lengths");
    }
    inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  };
  const auto grid_x = [&](std::size_t p) -> std::optional<double> {
    const std::size_t row = indexed ? static_cast<std::size_t>(target_view.coordinate_source_index[p]) : p;
    if (!std::isfinite(pos_x[row]) || !std::isfinite(pos_y[row]) || !std::isfinite(pos_z[row])) {
      return std::nullopt;  // rejected collectively by interpolateForces itself
    }
    if (!periodic &&
        (!positionInsideOpenDomain(pos_x[row], lengths.lx) ||
         !positionInsideOpenDomain(pos_y[row], lengths.ly) ||
         !positionInsideOpenDomain(pos_z[row], lengths.lz))) {
      return std::nullopt;
    }
    return (periodic ? wrapPosition(pos_x[row], lengths.lx) : pos_x[row]) * inv_dx;
  };
  auto& cache = m_impl->m_force_routing_plan;
  cache = {};
  PmRoutingPlanKey key{};
  std::uint64_t global_max_count = 0U, global_max_peer_records = 0U, local_max_round_records = 0U;
  const PmRoutingWorkloadPlan plan = planPmRoutingCollective(
      PmRoutingPhaseKind::kForceInterpolation, m_shape, grid, options, count, validate_local, grid_x,
      key, global_max_count, global_max_peer_records, local_max_round_records);
  if (plan.distributed) {
    cache = Impl::RoutingWorkloadPlanCache{
        .valid = true,
        .key = key,
        .global_max_count = global_max_count,
        .global_max_peer_records = global_max_peer_records,
        .local_max_round_records = local_max_round_records,
    };
  }
  return plan;
}

void PmSolver::interpolateForces(
    const PmGridStorage& grid,
    const PmForceTargetView& target_view,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
#if COSMOSIM_ENABLE_MPI
  double target_view_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kInterpolateForceTargetView,
        rank_local_serial_layout,
        m_shape,
        options,
        target_view_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces target-view API preflight",
        mpi_world_rank,
        mpi_world_size,
        target_view_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny ||
              grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in interpolateForces target view");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces target view requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces target-view slab metadata must match MPI_COMM_WORLD");
          }
          const std::size_t local_active_count = target_view.active_particle_index.size();
          const bool indexed_coordinates =
              target_view.coordinate_layout == PmForceCoordinateLayout::kIndexedSource;
          if (target_view.pos_x_comoving.size() != target_view.pos_y_comoving.size() ||
              target_view.pos_x_comoving.size() != target_view.pos_z_comoving.size()) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces coordinate source extents must match");
          }
          if (indexed_coordinates) {
            if (target_view.coordinate_source_index.size() != local_active_count) {
              throw std::invalid_argument(
                  "PmSolver::interpolateForces indexed coordinate map extent must match active count");
            }
          } else {
            if (!target_view.coordinate_source_index.empty()) {
              throw std::invalid_argument(
                  "PmSolver::interpolateForces compact coordinates must not provide source indices");
            }
            if (local_active_count != target_view.pos_x_comoving.size()) {
              throw std::invalid_argument(
                  "PmSolver::interpolateForces compact coordinate extent must match active count");
            }
          }
          switch (target_view.output_layout) {
            case PmForceOutputLayout::kCompactActive:
              if (local_active_count != target_view.accel_x_comoving.size() ||
                  local_active_count != target_view.accel_y_comoving.size() ||
                  local_active_count != target_view.accel_z_comoving.size()) {
                throw std::invalid_argument(
                    "PmSolver::interpolateForces compact active output extents must match active count");
              }
              break;
            case PmForceOutputLayout::kIndexedGlobal:
              for (const std::uint32_t particle_index : target_view.active_particle_index) {
                if (particle_index >= target_view.accel_x_comoving.size() ||
                    particle_index >= target_view.accel_y_comoving.size() ||
                    particle_index >= target_view.accel_z_comoving.size()) {
                  throw std::out_of_range(
                      "PmSolver::interpolateForces indexed active particle output index out of range");
                }
              }
              break;
            default:
              throw std::invalid_argument("PmSolver::interpolateForces unknown PM force output layout");
          }
        });
    if (profile != nullptr) {
      profile->routed_mpi_wait_ms += target_view_preflight_mpi_wait_ms;
    }
  }
#endif
  const std::size_t active_count = target_view.active_particle_index.size();
  const bool indexed_coordinates =
      target_view.coordinate_layout == PmForceCoordinateLayout::kIndexedSource;
  if (target_view.pos_x_comoving.size() != target_view.pos_y_comoving.size() ||
      target_view.pos_x_comoving.size() != target_view.pos_z_comoving.size()) {
    throw std::invalid_argument(
        "PmSolver::interpolateForces coordinate source extents must match");
  }
  if (indexed_coordinates) {
    if (target_view.coordinate_source_index.size() != active_count) {
      throw std::invalid_argument(
          "PmSolver::interpolateForces indexed coordinate map extent must match active count");
    }
  } else {
    if (!target_view.coordinate_source_index.empty()) {
      throw std::invalid_argument(
          "PmSolver::interpolateForces compact coordinates must not provide source indices");
    }
    if (active_count != target_view.pos_x_comoving.size()) {
      throw std::invalid_argument(
          "PmSolver::interpolateForces compact coordinate extent must match active count");
    }
  }

  auto& indexed_workspace = m_impl->indexedTargetWorkspace();
  switch (target_view.output_layout) {
    case PmForceOutputLayout::kCompactActive: {
      if (active_count != target_view.accel_x_comoving.size() ||
          active_count != target_view.accel_y_comoving.size() ||
          active_count != target_view.accel_z_comoving.size()) {
        throw std::invalid_argument(
            "PmSolver::interpolateForces compact active output extents must match active count");
      }
      interpolateForcesImpl(
          grid,
          target_view.pos_x_comoving,
          target_view.pos_y_comoving,
          target_view.pos_z_comoving,
          target_view.coordinate_source_index,
          target_view.coordinate_layout,
          target_view.accel_x_comoving,
          target_view.accel_y_comoving,
          target_view.accel_z_comoving,
          options,
          profile);
      return;
    }
    case PmForceOutputLayout::kIndexedGlobal: {
      for (const std::uint32_t particle_index : target_view.active_particle_index) {
        if (particle_index >= target_view.accel_x_comoving.size() ||
            particle_index >= target_view.accel_y_comoving.size() ||
            particle_index >= target_view.accel_z_comoving.size()) {
          throw std::out_of_range(
              "PmSolver::interpolateForces indexed active particle output index out of range");
        }
      }
      indexed_workspace.compact_ax.assign(active_count, 0.0);
      indexed_workspace.compact_ay.assign(active_count, 0.0);
      indexed_workspace.compact_az.assign(active_count, 0.0);
      interpolateForcesImpl(
          grid,
          target_view.pos_x_comoving,
          target_view.pos_y_comoving,
          target_view.pos_z_comoving,
          target_view.coordinate_source_index,
          target_view.coordinate_layout,
          indexed_workspace.compact_ax,
          indexed_workspace.compact_ay,
          indexed_workspace.compact_az,
          options,
          profile);
      for (std::size_t active_i = 0; active_i < active_count; ++active_i) {
        const std::uint32_t particle_index = target_view.active_particle_index[active_i];
        target_view.accel_x_comoving[particle_index] = indexed_workspace.compact_ax[active_i];
        target_view.accel_y_comoving[particle_index] = indexed_workspace.compact_ay[active_i];
        target_view.accel_z_comoving[particle_index] = indexed_workspace.compact_az[active_i];
      }
      return;
    }
  }
  throw std::invalid_argument("PmSolver::interpolateForces unknown PM force output layout");
}

void PmSolver::interpolateForces(
    const PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<double> accel_x,
    std::span<double> accel_y,
    std::span<double> accel_z,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
  interpolateForcesImpl(
      grid, pos_x, pos_y, pos_z, {}, PmForceCoordinateLayout::kCompactActive,
      accel_x, accel_y, accel_z, options, profile);
}

void PmSolver::interpolateForcesImpl(
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
    PmProfileEvent* profile) const {
  const bool indexed_coordinates =
      coordinate_layout == PmForceCoordinateLayout::kIndexedSource;
  const std::size_t target_count =
      indexed_coordinates ? coordinate_source_index.size() : pos_x.size();

  const auto validate_coordinate_contract = [&](bool validate_source_indices) {
    if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size()) {
      throw std::invalid_argument(
          "PmSolver::interpolateForces coordinate source extents must match");
    }
    switch (coordinate_layout) {
      case PmForceCoordinateLayout::kCompactActive:
        if (!coordinate_source_index.empty()) {
          throw std::invalid_argument(
              "PmSolver::interpolateForces compact coordinates must not provide source indices");
        }
        break;
      case PmForceCoordinateLayout::kIndexedSource:
        if (validate_source_indices) {
          for (std::size_t i = 0; i < coordinate_source_index.size(); ++i) {
            if (static_cast<std::size_t>(coordinate_source_index[i]) >= pos_x.size()) {
              throw std::out_of_range(
                  "PmSolver::interpolateForces coordinate source index is out of range; active_index=" +
                  std::to_string(i));
            }
          }
        }
        break;
      default:
        throw std::invalid_argument(
            "PmSolver::interpolateForces received an unsupported coordinate layout");
    }
    if (target_count != accel_x.size() || target_count != accel_y.size() ||
        target_count != accel_z.size()) {
      throw std::invalid_argument(
          "PmSolver::interpolateForces target/output extents must match");
    }
  };

  const PmDecompositionView decomposition_view(
      m_shape, grid.slabLayout(), options.decomposition_mode);
  bool coordinate_indices_collectively_validated = false;
#if COSMOSIM_ENABLE_MPI
  double distributed_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kInterpolateForces,
        rank_local_serial_layout,
        m_shape,
        options,
        distributed_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces distributed API preflight",
        mpi_world_rank,
        mpi_world_size,
        distributed_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in interpolateForces");
          }
          validate_coordinate_contract(/*validate_source_indices=*/true);
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::interpolateForces requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces mesh dimensions exceed fixed-width routing indices");
          }
        });
    coordinate_indices_collectively_validated = true;
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in interpolateForces");
  }
  validate_coordinate_contract(!coordinate_indices_collectively_validated);
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::interpolateForces requires a valid PM slab layout");
  }

  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if !COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    throw std::invalid_argument("PmSolver::interpolateForces distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
  validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolateForces");
#else
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument("PmSolver::interpolateForces slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolateForces");
  }
#endif

  // The coordinate source is non-owning. The layout decision is invariant
  // for the call, so this predictable row indirection replaces the former
  // O(active_count) XYZ gather without branching inside the PM stencil loops.
  const auto coordinate_row = [&](std::size_t target_row) noexcept -> std::size_t {
    return indexed_coordinates
        ? static_cast<std::size_t>(coordinate_source_index[target_row])
        : target_row;
  };

  const auto start = std::chrono::steady_clock::now();

  const BoxLengths lengths = effectiveBoxLengths(options);
  const double inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  const double inv_dy = static_cast<double>(m_shape.ny) / lengths.ly;
  const double inv_dz = static_cast<double>(m_shape.nz) / lengths.lz;

  if (distributed_slabs &&
      (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
       m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
       m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))) {
    throw std::invalid_argument("PmSolver::interpolateForces mesh dimensions exceed fixed-width routing indices");
  }
  std::optional<std::size_t> first_non_finite_target_index;
  for (std::size_t p = 0; p < target_count; ++p) {
    const std::size_t source_row = coordinate_row(p);
    if (!std::isfinite(pos_x[source_row]) || !std::isfinite(pos_y[source_row]) ||
        !std::isfinite(pos_z[source_row])) {
      first_non_finite_target_index = p;
      break;
    }
  }
  if (!distributed_slabs && first_non_finite_target_index.has_value()) {
    throw std::invalid_argument(
        "PmSolver::interpolateForces requires finite particle coordinates; particle_index=" +
        std::to_string(*first_non_finite_target_index));
  }

  if (!distributed_slabs) {
#if COSMOSIM_HAVE_OPENMP
#pragma omp parallel for schedule(static) if(target_count >= 1024U)
#endif
    for (std::size_t p = 0; p < target_count; ++p) {
      const std::size_t source_row = coordinate_row(p);
      const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
      if (!periodic &&
          (!positionInsideOpenDomain(pos_x[source_row], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[source_row], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[source_row], lengths.lz))) {
        accel_x[p] = 0.0;
        accel_y[p] = 0.0;
        accel_z[p] = 0.0;
        continue;
      }
      const double x = (periodic ? wrapPosition(pos_x[source_row], lengths.lx) : pos_x[source_row]) * inv_dx;
      const double y = (periodic ? wrapPosition(pos_y[source_row], lengths.ly) : pos_y[source_row]) * inv_dy;
      const double z = (periodic ? wrapPosition(pos_z[source_row], lengths.lz) : pos_z[source_row]) * inv_dz;

      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      double gx = 0.0;
      double gy = 0.0;
      double gz = 0.0;

      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        if (!periodic && (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<int>(m_shape.nx))) {
          continue;
        }
        const std::size_t ix = periodic ? wrapIndex(sx.offsets[dx], m_shape.nx) : static_cast<std::size_t>(sx.offsets[dx]);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic && (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<int>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic ? wrapIndex(sy.offsets[dy], m_shape.ny) : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic && (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<int>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic ? wrapIndex(sz.offsets[dz], m_shape.nz) : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            const std::size_t index = grid.linearIndex(ix, iy, iz);
            gx += weight * grid.force_x()[index];
            gy += weight * grid.force_y()[index];
            gz += weight * grid.force_z()[index];
          }
        }
      }

      accel_x[p] = gx;
      accel_y[p] = gy;
      accel_z[p] = gz;
    }
  } else {
#if COSMOSIM_ENABLE_MPI
    const int world_size = grid.slabLayout().world_size;
    const int world_rank = grid.slabLayout().world_rank;
    double routed_mpi_wait_ms = distributed_preflight_mpi_wait_ms;
    const int local_targets_valid = first_non_finite_target_index.has_value() ? 0 : 1;
    int all_targets_valid = 0;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(&local_targets_valid, &all_targets_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD),
          "PmSolver::interpolateForces target validation MPI_Allreduce");
    });
    if (all_targets_valid == 0) {
      const std::string local_detail = first_non_finite_target_index.has_value()
          ? " local_first_invalid_particle_index=" + std::to_string(*first_non_finite_target_index)
          : " local_targets_valid=true";
      throw std::invalid_argument(
          "PmSolver::interpolateForces rejected non-finite target coordinates on at least one MPI rank;" +
          local_detail);
    }

    std::uint64_t exchange_epoch = 0U;
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces exchange epoch allocation",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() { exchange_epoch = m_impl->nextDistributedExchangeEpoch(); });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      validatePmExchangeEpochConsensus(exchange_epoch, "PmSolver::interpolateForces");
    });

    std::optional<GravityCommunicationArena::Lease> standalone_communication_lease;
    std::optional<Impl::PlaneInterpolationExchangeBuffers> exchange_storage;
    Impl::PlaneInterpolationExchangeBuffers* exchange_ptr = nullptr;
    PmRoutingCapacityModel routing_capacity{};
    std::size_t routing_buffer_limit = 0U;
    std::size_t particles_per_round = 0U;
    const bool routing_periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    const auto target_grid_x = [&](std::size_t p) -> std::optional<double> {
      const std::size_t source_row = coordinate_row(p);
      if (!routing_periodic &&
          (!positionInsideOpenDomain(pos_x[source_row], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[source_row], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[source_row], lengths.lz))) {
        return std::nullopt;
      }
      return (routing_periodic ? wrapPosition(pos_x[source_row], lengths.lx) : pos_x[source_row]) * inv_dx;
    };
    // Same four-step routing state machine as assignDensity (workload phase,
    // agreement reduction, preflight phase, bounded rounds).
    bool local_plan_present = false;
    bool local_plan_violated = false;
    PmRemotePlaneWorkload local_workload{};
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces routing workload",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
          const Impl::RoutingWorkloadPlanCache plan = m_impl->m_force_routing_plan;
          m_impl->m_force_routing_plan = {};  // withdrawn whatever the outcome
          local_plan_present = plan.valid;
          routing_capacity = pmPlannedRoutingModel(
              PmRoutingPhaseKind::kForceInterpolation, world_size, options.routing_exchange_batch_bytes);
          if (routing_capacity.max_send_payload_bytes >
              static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            throw std::overflow_error("PM force routing buffer limit exceeds size_t capacity");
          }
          routing_buffer_limit =
              static_cast<std::size_t>(routing_capacity.max_send_payload_bytes);
          particles_per_round = pmRoutingParticlesPerRound(
              routing_capacity.effective_per_peer_payload_bytes,
              k_pm_plane_interpolation_request_wire_bytes);
          local_workload = measurePmRemotePlaneWorkload(
              target_count, particles_per_round, world_rank, world_size, routing_periodic,
              options.assignment_scheme, m_shape, decomposition_view, target_grid_x);
          if (local_plan_present) {
            local_plan_violated = !pmRoutingPlanBoundsWorkload(
                plan.key, plan.global_max_count, plan.global_max_peer_records, plan.local_max_round_records,
                makePmRoutingPlanKey(m_shape, grid.slabLayout(), options),
                static_cast<std::uint64_t>(target_count), local_workload);
          }
        });
    const std::uint64_t local_target_count = static_cast<std::uint64_t>(target_count);
    // Replaces the former max-target reduction (one MAX over agreed lanes).
    const PmRoutingAgreement routing_agreement = agreePmRoutingWorkload(
        "PmSolver::interpolateForces", world_rank, world_size,
        local_target_count, local_workload.max_peer_records,
        local_plan_present, local_plan_violated, routed_mpi_wait_ms);
    const std::uint64_t global_max_target_count = routing_agreement.global_max_count;
    const std::uint64_t global_max_peer_records = routing_agreement.global_max_peer_records;
    // Responses (64 B) are compacted in place inside the request buffers
    // (96 B), so the request bounds also bound the response traffic.
    const PmRoutingReservation wire_reservation = pmRoutingReservation(
        PmRoutingPhaseKind::kForceInterpolation, routing_capacity, world_size,
        local_workload.max_round_records, global_max_peer_records);
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces routing preflight",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
          std::pmr::memory_resource* communication_resource = m_impl->communicationResource(
              GravityCommunicationArena::Phase::kPmInterpolation,
              standalone_communication_lease,
              wire_reservation.phase_arena_bytes);
          exchange_storage.emplace(
              communication_resource, static_cast<std::size_t>(world_size));
          exchange_ptr = &*exchange_storage;
          exchange_ptr->send_wire.reserve(static_cast<std::size_t>(wire_reservation.send_bytes));
          exchange_ptr->recv_wire.reserve(static_cast<std::size_t>(wire_reservation.receive_bytes));
        });
    auto& exchange = *exchange_ptr;
    const std::uint64_t routing_metadata_capacity =
        pmRoutingMetadataBytes(PmRoutingPhaseKind::kForceInterpolation, world_size);
    const std::uint64_t round_capacity = static_cast<std::uint64_t>(particles_per_round);
    const std::uint64_t round_count = global_max_target_count == 0U
        ? 0U
        : 1U + (global_max_target_count - 1U) / round_capacity;
    const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    std::uint32_t next_request_sequence = 0U;

    std::fill(accel_x.begin(), accel_x.end(), 0.0);
    std::fill(accel_y.begin(), accel_y.end(), 0.0);
    std::fill(accel_z.begin(), accel_z.end(), 0.0);

    const auto evaluate_plane_request = [&](const PmPlaneInterpolationRequestRecord& request) {
      if (request.plane_count == 0U || request.plane_count > 3U ||
          request.destination_rank != static_cast<std::uint32_t>(world_rank) ||
          request.exchange_epoch != exchange_epoch ||
          !std::isfinite(request.y_grid) || !std::isfinite(request.z_grid)) {
        throw std::invalid_argument(
            "PmSolver::interpolateForces received invalid bounded force-plane request");
      }
      const PmAxisStencil1d sy = makeAxisStencil(request.y_grid, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(request.z_grid, options.assignment_scheme);
      PmForceContributionRecord response{
          .source_rank = static_cast<std::uint32_t>(world_rank),
          .origin_rank = request.origin_rank,
          .request_sequence = request.request_sequence,
          .origin_particle_index = request.origin_particle_index,
          .exchange_epoch = request.exchange_epoch,
          .accel_x = 0.0,
          .accel_y = 0.0,
          .accel_z = 0.0,
      };
      for (std::size_t plane = 0; plane < static_cast<std::size_t>(request.plane_count); ++plane) {
        const std::size_t ix = static_cast<std::size_t>(request.global_ix[plane]);
        if (ix >= m_shape.nx || !grid.slabLayout().ownsGlobalX(ix) ||
            !std::isfinite(request.x_weight[plane])) {
          throw std::invalid_argument(
              "PmSolver::interpolateForces force-plane x ownership/weight is invalid");
        }
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic &&
              (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic &&
                (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const std::size_t local_index = decomposition_view.globalToLocalRealIndex(ix, iy, iz);
            if (!std::isfinite(grid.force_x()[local_index]) ||
                !std::isfinite(grid.force_y()[local_index]) ||
                !std::isfinite(grid.force_z()[local_index])) {
              throw std::invalid_argument(
                  "PmSolver::interpolateForces encountered non-finite owner-local mesh force data");
            }
            const double weight = request.x_weight[plane] * sy.weights[dy] * sz.weights[dz];
            response.accel_x += weight * grid.force_x()[local_index];
            response.accel_y += weight * grid.force_y()[local_index];
            response.accel_z += weight * grid.force_z()[local_index];
          }
        }
      }
      return response;
    };

    std::uint64_t logical_requests_sent = 0U;
    std::uint64_t peers_touched = 0U;
    std::uint64_t total_wire_sent = 0U;
    std::uint64_t total_wire_received = 0U;
    std::uint64_t send_high_water = 0U;
    std::uint64_t receive_high_water = 0U;
    std::uint64_t combined_high_water = 0U;
    std::uint64_t workspace_high_water = routing_metadata_capacity;

    for (std::uint64_t round = 0U; round < round_count; ++round) {
      const std::uint64_t begin_u64 = round * round_capacity;
      const std::uint64_t end_u64 = std::min(
          local_target_count,
          begin_u64 > std::numeric_limits<std::uint64_t>::max() - round_capacity
              ? local_target_count
              : begin_u64 + round_capacity);
      const std::size_t begin = static_cast<std::size_t>(begin_u64);
      const std::size_t end = static_cast<std::size_t>(end_u64);
      const std::uint32_t round_sequence_begin = next_request_sequence;

      std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
      std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
      std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
      std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
      std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
      std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
      std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
      std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
      std::fill(exchange.send_response_counts_bytes.begin(), exchange.send_response_counts_bytes.end(), 0);
      std::fill(exchange.send_response_displs_bytes.begin(), exchange.send_response_displs_bytes.end(), 0);
      std::fill(exchange.recv_response_counts_bytes.begin(), exchange.recv_response_counts_bytes.end(), 0);
      std::fill(exchange.recv_response_displs_bytes.begin(), exchange.recv_response_displs_bytes.end(), 0);
      exchange.send_wire.clear();
      exchange.recv_wire.clear();

      runPmCoordinatedPhase(
          "PmSolver::interpolateForces bounded request preparation",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            // Pass 1 counts at most one request per particle/destination x-plane group.
            for (std::size_t p = begin; p < end; ++p) {
              const std::size_t source_row = coordinate_row(p);
              if (p > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
                throw std::invalid_argument(
                    "PmSolver::interpolateForces origin particle index exceeds routing token limit");
              }
              if (!periodic &&
                  (!positionInsideOpenDomain(pos_x[source_row], lengths.lx) ||
                   !positionInsideOpenDomain(pos_y[source_row], lengths.ly) ||
                   !positionInsideOpenDomain(pos_z[source_row], lengths.lz))) {
                continue;
              }
              const double x = (periodic ? wrapPosition(pos_x[source_row], lengths.lx) : pos_x[source_row]) * inv_dx;
              const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const int destination_rank = groups[group_i].destination_rank;
                if (destination_rank == world_rank) {
                  continue;
                }
                int& count = exchange.send_counts[static_cast<std::size_t>(destination_rank)];
                if (count == std::numeric_limits<int>::max()) {
                  throw std::overflow_error("PM force per-peer request count exceeds MPI int capacity");
                }
                ++count;
              }
            }

            std::size_t total_send_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = static_cast<std::size_t>(
                  exchange.send_counts[static_cast<std::size_t>(rank)]);
              const std::uint64_t peer_bytes = checkedBytesForCount(
                  count, k_pm_plane_interpolation_request_wire_bytes, "PM force per-peer routing batch");
              if (peer_bytes > routing_capacity.effective_per_peer_payload_bytes) {
                throw std::runtime_error(
                    "PM force routing exceeded effective per-peer batch high-water");
              }
              exchange.send_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_send_records,
                      k_pm_plane_interpolation_request_wire_bytes,
                      "PmSolver::interpolateForces bounded request exchange",
                      "send record displacement",
                      rank);
              total_send_records = checkedMpiRecordTotal(
                  total_send_records,
                  count,
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolateForces bounded request exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.send_wire,
                checkedPmWireByteCount(
                    total_send_records,
                    k_pm_plane_interpolation_request_wire_bytes,
                    "PmSolver::interpolateForces bounded request exchange"),
                routing_buffer_limit,
                "PmSolver::interpolateForces request send buffer");

            for (int rank = 0; rank < world_size; ++rank) {
              exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
                  exchange.send_displs[static_cast<std::size_t>(rank)]);
            }

            // Pass 2 evaluates local x planes immediately and encodes only remote
            // x-plane groups. The receiver expands the y/z stencil and returns one
            // accumulated force response per request.
            for (std::size_t p = begin; p < end; ++p) {
              const std::size_t source_row = coordinate_row(p);
              if (!periodic &&
                  (!positionInsideOpenDomain(pos_x[source_row], lengths.lx) ||
                   !positionInsideOpenDomain(pos_y[source_row], lengths.ly) ||
                   !positionInsideOpenDomain(pos_z[source_row], lengths.lz))) {
                continue;
              }
              const double x = (periodic ? wrapPosition(pos_x[source_row], lengths.lx) : pos_x[source_row]) * inv_dx;
              const double y = (periodic ? wrapPosition(pos_y[source_row], lengths.ly) : pos_y[source_row]) * inv_dy;
              const double z = (periodic ? wrapPosition(pos_z[source_row], lengths.lz) : pos_z[source_row]) * inv_dz;
              const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const auto& group = groups[group_i];
                PmPlaneInterpolationRequestRecord request{
                    .origin_rank = static_cast<std::uint32_t>(world_rank),
                    .destination_rank = static_cast<std::uint32_t>(group.destination_rank),
                    .request_sequence = next_request_sequence,
                    .origin_particle_index = static_cast<std::uint32_t>(p),
                    .plane_count = group.plane_count,
                    .global_ix = group.global_ix,
                    .exchange_epoch = exchange_epoch,
                    .y_grid = y,
                    .z_grid = z,
                    .x_weight = group.x_weight,
                };
                if (group.destination_rank == world_rank) {
                  const PmForceContributionRecord local = evaluate_plane_request(request);
                  accel_x[p] += local.accel_x;
                  accel_y[p] += local.accel_y;
                  accel_z[p] += local.accel_z;
                  continue;
                }
                if (next_request_sequence == std::numeric_limits<std::uint32_t>::max()) {
                  throw std::overflow_error(
                      "PmSolver::interpolateForces routing sequence exceeds uint32 capacity");
                }
                auto& cursor = exchange.cursor[static_cast<std::size_t>(group.destination_rank)];
                encodePmPlaneInterpolationRequest(
                    request,
                    PmWireRecordKind::kForcePlaneRequest,
                    std::span<std::uint8_t>(exchange.send_wire).subspan(
                        cursor * k_pm_plane_interpolation_request_wire_bytes,
                        k_pm_plane_interpolation_request_wire_bytes));
                ++cursor;
                ++next_request_sequence;
              }
            }
          });

      const std::uint64_t request_send_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoall(
                exchange.send_counts.data(), 1, MPI_INT,
                exchange.recv_counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD),
            "PmSolver::interpolateForces request-count MPI_Alltoall");
      });

      runPmCoordinatedPhase(
          "PmSolver::interpolateForces bounded request receive layout",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            std::size_t total_recv_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = checkedMpiReceivedRecordCount(
                  exchange.recv_counts[static_cast<std::size_t>(rank)],
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolateForces bounded request exchange",
                  rank);
              exchange.recv_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_recv_records,
                      k_pm_plane_interpolation_request_wire_bytes,
                      "PmSolver::interpolateForces bounded request exchange",
                      "received record displacement",
                      rank);
              total_recv_records = checkedMpiRecordTotal(
                  total_recv_records,
                  count,
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolateForces bounded request exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.recv_wire,
                checkedPmWireByteCount(
                    total_recv_records,
                    k_pm_plane_interpolation_request_wire_bytes,
                    "PmSolver::interpolateForces bounded request exchange"),
                routing_buffer_limit,
                "PmSolver::interpolateForces request receive buffer");
            checkedMpiRecordLayoutToByteLayout(
                exchange.send_counts,
                exchange.send_displs,
                k_pm_plane_interpolation_request_wire_bytes,
                exchange.send_counts_bytes,
                exchange.send_displs_bytes,
                "PmSolver::interpolateForces bounded request send MPI_Alltoallv");
            checkedMpiRecordLayoutToByteLayout(
                exchange.recv_counts,
                exchange.recv_displs,
                k_pm_plane_interpolation_request_wire_bytes,
                exchange.recv_counts_bytes,
                exchange.recv_displs_bytes,
                "PmSolver::interpolateForces bounded request receive MPI_Alltoallv");
          });
      const std::uint64_t request_recv_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());

      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.send_wire),
                exchange.send_counts_bytes.data(),
                exchange.send_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.recv_wire),
                exchange.recv_counts_bytes.data(),
                exchange.recv_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::interpolateForces request MPI_Alltoallv");
      });

      std::string local_force_request_validation_error;
      try {
        // Compact responses in-place into the receive-request buffer. Since a
        // response (64 B) is smaller than a request (96 B), forward compaction
        // cannot overwrite a request that has not yet been decoded.
        for (int source_rank = 0; source_rank < world_size; ++source_rank) {
          const int request_begin = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
          const int request_count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
          exchange.send_response_counts_bytes[static_cast<std::size_t>(source_rank)] =
              checkedMpiByteCount(
                  checkedBytesForCount(
                      static_cast<std::size_t>(request_count),
                      k_pm_force_response_wire_bytes,
                      "PM force response sender count"),
                  "PM force response sender count");
          exchange.send_response_displs_bytes[static_cast<std::size_t>(source_rank)] =
              checkedMpiByteCount(
                  checkedBytesForCount(
                      static_cast<std::size_t>(request_begin),
                      k_pm_force_response_wire_bytes,
                      "PM force response sender displacement"),
                  "PM force response sender displacement");
          std::optional<std::uint32_t> previous_sequence;
          for (int i = 0; i < request_count; ++i) {
            const std::size_t request_index = static_cast<std::size_t>(request_begin + i);
            const PmPlaneInterpolationRequestRecord request = decodePmPlaneInterpolationRequest(
                std::span<const std::uint8_t>(exchange.recv_wire).subspan(
                    request_index * k_pm_plane_interpolation_request_wire_bytes,
                    k_pm_plane_interpolation_request_wire_bytes),
                PmWireRecordKind::kForcePlaneRequest,
                "PmSolver::interpolateForces force-plane request");
            if (request.origin_rank != static_cast<std::uint32_t>(source_rank) ||
                request.destination_rank != static_cast<std::uint32_t>(world_rank) ||
                request.exchange_epoch != exchange_epoch) {
              throw std::invalid_argument(
                  "PM force-plane routing metadata does not match MPI sender/receiver segment");
            }
            if (previous_sequence.has_value() && request.request_sequence <= *previous_sequence) {
              throw std::invalid_argument(
                  "PM force-plane request sequence is non-monotonic within sender segment");
            }
            previous_sequence = request.request_sequence;
            const PmForceContributionRecord response = evaluate_plane_request(request);
            const std::size_t response_index =
                static_cast<std::size_t>(exchange.send_response_displs_bytes[static_cast<std::size_t>(source_rank)]) /
                    k_pm_force_response_wire_bytes +
                static_cast<std::size_t>(i);
            encodePmForceResponse(
                response,
                std::span<std::uint8_t>(exchange.recv_wire).subspan(
                    response_index * k_pm_force_response_wire_bytes,
                    k_pm_force_response_wire_bytes));
          }
        }
      } catch (const std::exception& error) {
        local_force_request_validation_error = error.what();
      } catch (...) {
        local_force_request_validation_error =
            "non-standard exception while consuming bounded force-plane requests";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::interpolateForces bounded request receive",
          world_rank,
          world_size,
          local_force_request_validation_error,
          routed_mpi_wait_ms);

      const std::size_t total_received_requests = exchange.recv_wire.size() /
          k_pm_plane_interpolation_request_wire_bytes;
      resizePmWireBufferBounded(
          exchange.recv_wire,
          checkedPmWireByteCount(
              total_received_requests,
              k_pm_force_response_wire_bytes,
              "PmSolver::interpolateForces compact response send"),
          routing_buffer_limit,
          "PmSolver::interpolateForces response send buffer");

      const std::size_t total_sent_requests = exchange.send_wire.size() /
          k_pm_plane_interpolation_request_wire_bytes;
      for (int rank = 0; rank < world_size; ++rank) {
        exchange.recv_response_counts_bytes[static_cast<std::size_t>(rank)] =
            checkedMpiByteCount(
                checkedBytesForCount(
                    static_cast<std::size_t>(exchange.send_counts[static_cast<std::size_t>(rank)]),
                    k_pm_force_response_wire_bytes,
                    "PM force response receive count"),
                "PM force response receive count");
        exchange.recv_response_displs_bytes[static_cast<std::size_t>(rank)] =
            checkedMpiByteCount(
                checkedBytesForCount(
                    static_cast<std::size_t>(exchange.send_displs[static_cast<std::size_t>(rank)]),
                    k_pm_force_response_wire_bytes,
                    "PM force response receive displacement"),
                "PM force response receive displacement");
      }
      resizePmWireBufferBounded(
          exchange.send_wire,
          checkedPmWireByteCount(
              total_sent_requests,
              k_pm_force_response_wire_bytes,
              "PmSolver::interpolateForces compact response receive"),
          routing_buffer_limit,
          "PmSolver::interpolateForces response receive buffer");

      const std::uint64_t response_send_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());
      const std::uint64_t response_recv_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.recv_wire),
                exchange.send_response_counts_bytes.data(),
                exchange.send_response_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.send_wire),
                exchange.recv_response_counts_bytes.data(),
                exchange.recv_response_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::interpolateForces response MPI_Alltoallv");
      });

      std::string local_force_response_validation_error;
      try {
        // Re-generate the compact request stream instead of retaining an O(N)
        // request registry. This provides exact response identity validation with
        // only O(world_size) cursor state.
        for (int rank = 0; rank < world_size; ++rank) {
          exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
              exchange.send_displs[static_cast<std::size_t>(rank)]);
        }
        std::uint32_t expected_sequence = round_sequence_begin;
        for (std::size_t p = begin; p < end; ++p) {
          const std::size_t source_row = coordinate_row(p);
          if (!periodic &&
              (!positionInsideOpenDomain(pos_x[source_row], lengths.lx) ||
               !positionInsideOpenDomain(pos_y[source_row], lengths.ly) ||
               !positionInsideOpenDomain(pos_z[source_row], lengths.lz))) {
            continue;
          }
          const double x = (periodic ? wrapPosition(pos_x[source_row], lengths.lx) : pos_x[source_row]) * inv_dx;
          const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
          std::array<PmXPlaneGroup, 3> groups{};
          const std::size_t group_count = makePmXPlaneGroups(
              sx, periodic, m_shape, decomposition_view, groups);
          for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
            const int sender_rank = groups[group_i].destination_rank;
            if (sender_rank == world_rank) {
              continue;
            }
            auto& cursor = exchange.cursor[static_cast<std::size_t>(sender_rank)];
            const std::size_t segment_end = static_cast<std::size_t>(
                exchange.send_displs[static_cast<std::size_t>(sender_rank)] +
                exchange.send_counts[static_cast<std::size_t>(sender_rank)]);
            if (cursor >= segment_end) {
              throw std::invalid_argument("PM force response segment ended before issued requests");
            }
            const PmForceContributionRecord response = decodePmForceResponse(
                std::span<const std::uint8_t>(exchange.send_wire).subspan(
                    cursor * k_pm_force_response_wire_bytes,
                    k_pm_force_response_wire_bytes));
            if (response.source_rank != static_cast<std::uint32_t>(sender_rank) ||
                response.origin_rank != static_cast<std::uint32_t>(world_rank) ||
                response.exchange_epoch != exchange_epoch ||
                response.request_sequence != expected_sequence ||
                response.origin_particle_index != static_cast<std::uint32_t>(p) ||
                !std::isfinite(response.accel_x) ||
                !std::isfinite(response.accel_y) ||
                !std::isfinite(response.accel_z)) {
              throw std::invalid_argument(
                  "PM force response does not match the regenerated bounded request stream");
            }
            accel_x[p] += response.accel_x;
            accel_y[p] += response.accel_y;
            accel_z[p] += response.accel_z;
            ++cursor;
            ++expected_sequence;
          }
        }
        if (expected_sequence != next_request_sequence) {
          throw std::logic_error("PM force response validation sequence accounting drifted");
        }
        for (int rank = 0; rank < world_size; ++rank) {
          const std::size_t expected_end = static_cast<std::size_t>(
              exchange.send_displs[static_cast<std::size_t>(rank)] +
              exchange.send_counts[static_cast<std::size_t>(rank)]);
          if (exchange.cursor[static_cast<std::size_t>(rank)] != expected_end) {
            throw std::invalid_argument("PM force response segment contains unconsumed records");
          }
        }
      } catch (const std::exception& error) {
        local_force_response_validation_error = error.what();
      } catch (...) {
        local_force_response_validation_error =
            "non-standard exception while validating bounded force responses";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::interpolateForces bounded response receive",
          world_rank,
          world_size,
          local_force_response_validation_error,
          routed_mpi_wait_ms);

      logical_requests_sent += static_cast<std::uint64_t>(total_sent_requests);
      peers_touched += static_cast<std::uint64_t>(std::count_if(
          exchange.send_counts.begin(), exchange.send_counts.end(), [](int count) { return count > 0; }));
      total_wire_sent += request_send_bytes + response_send_bytes;
      total_wire_received += request_recv_bytes + response_recv_bytes;
      const std::uint64_t send_wire_capacity =
          static_cast<std::uint64_t>(exchange.send_wire.capacity());
      const std::uint64_t recv_wire_capacity =
          static_cast<std::uint64_t>(exchange.recv_wire.capacity());
      // The two physical buffers swap direction for the response phase, so
      // each directional high-water must consider both retained capacities.
      const std::uint64_t directional_capacity =
          std::max(send_wire_capacity, recv_wire_capacity);
      const std::uint64_t combined_capacity = checkedAddBytes(
          send_wire_capacity,
          recv_wire_capacity,
          "PM force combined routing buffer capacity");
      const std::uint64_t workspace_capacity = checkedAddBytes(
          routing_metadata_capacity,
          combined_capacity,
          "PM force routing workspace capacity");
      if (workspace_capacity > k_pm_routing_workspace_target_bytes) {
        throw std::logic_error("PM force routing workspace exceeded the M1A aggregate target");
      }
      send_high_water = std::max(send_high_water, directional_capacity);
      receive_high_water = std::max(receive_high_water, directional_capacity);
      combined_high_water = std::max(combined_high_water, combined_capacity);
      workspace_high_water = std::max(workspace_high_water, workspace_capacity);
      exchange.workspace_high_water_bytes = std::max(
          exchange.workspace_high_water_bytes, workspace_capacity);
      m_impl->recordCommunicationHighWater(
          GravityCommunicationArena::Phase::kPmInterpolation,
          workspace_capacity,
          standalone_communication_lease);
    }

    if (profile != nullptr) {
      profile->routed_force_requests += logical_requests_sent;
      profile->routed_force_peer_count += peers_touched;
      profile->routed_mpi_bytes_sent += total_wire_sent;
      profile->routed_mpi_bytes_received += total_wire_received;
      profile->routed_send_buffer_high_water_bytes = std::max(
          profile->routed_send_buffer_high_water_bytes, send_high_water);
      profile->routed_receive_buffer_high_water_bytes = std::max(
          profile->routed_receive_buffer_high_water_bytes, receive_high_water);
      profile->routed_combined_buffer_high_water_bytes = std::max(
          profile->routed_combined_buffer_high_water_bytes, combined_high_water);
      profile->routed_workspace_high_water_bytes = std::max(
          profile->routed_workspace_high_water_bytes, workspace_high_water);
      profile->routed_mpi_wait_ms += routed_mpi_wait_ms;
      profile->bytes_moved += total_wire_sent + total_wire_received;
    }
#endif
  }

  const auto stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->interpolate_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    profile->bytes_moved += bytesForParticles(target_count);
  }
}

void PmSolver::interpolatePotential(
    const PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<double> potential,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
  const PmDecompositionView decomposition_view(
      m_shape, grid.slabLayout(), options.decomposition_mode);
#if COSMOSIM_ENABLE_MPI
  double distributed_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kInterpolatePotential,
        rank_local_serial_layout,
        m_shape,
        options,
        distributed_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential distributed API preflight",
        mpi_world_rank,
        mpi_world_size,
        distributed_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in interpolatePotential");
          }
          if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != potential.size()) {
            throw std::invalid_argument("Particle coordinate/potential spans must match in interpolatePotential");
          }
          if (!grid.hasPotentialStorage()) {
            throw std::logic_error(
                "PmSolver::interpolatePotential requires materialized real-space PM potential storage");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::interpolatePotential requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::interpolatePotential slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
            throw std::invalid_argument(
                "PmSolver::interpolatePotential mesh dimensions exceed fixed-width routing indices");
          }
        });
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in interpolatePotential");
  }
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != potential.size()) {
    throw std::invalid_argument("Particle coordinate/potential spans must match in interpolatePotential");
  }
  if (!grid.hasPotentialStorage()) {
    throw std::logic_error(
        "PmSolver::interpolatePotential requires materialized real-space PM potential storage");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::interpolatePotential requires a valid PM slab layout");
  }

  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if !COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    throw std::invalid_argument("PmSolver::interpolatePotential distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
  validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolatePotential");
#else
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument(
          "PmSolver::interpolatePotential slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolatePotential");
  }
#endif

  const auto start = std::chrono::steady_clock::now();
  const BoxLengths lengths = effectiveBoxLengths(options);
  const double inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  const double inv_dy = static_cast<double>(m_shape.ny) / lengths.ly;
  const double inv_dz = static_cast<double>(m_shape.nz) / lengths.lz;
  if (distributed_slabs &&
      (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
       m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
       m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))) {
    throw std::invalid_argument("PmSolver::interpolatePotential mesh dimensions exceed fixed-width routing indices");
  }
  std::optional<std::size_t> first_non_finite_target_index;
  for (std::size_t p = 0; p < pos_x.size(); ++p) {
    if (!std::isfinite(pos_x[p]) || !std::isfinite(pos_y[p]) || !std::isfinite(pos_z[p])) {
      first_non_finite_target_index = p;
      break;
    }
  }
  if (!distributed_slabs && first_non_finite_target_index.has_value()) {
    throw std::invalid_argument(
        "PmSolver::interpolatePotential requires finite particle coordinates; particle_index=" +
        std::to_string(*first_non_finite_target_index));
  }

  if (!distributed_slabs) {
    for (std::size_t p = 0; p < pos_x.size(); ++p) {
      const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
      if (!periodic &&
          (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
        potential[p] = 0.0;
        continue;
      }
      const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
      const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
      const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;

      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      double phi = 0.0;
      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        if (!periodic && (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<std::ptrdiff_t>(m_shape.nx))) {
          continue;
        }
        const std::size_t ix = periodic
            ? wrapIndex(sx.offsets[dx], m_shape.nx)
            : static_cast<std::size_t>(sx.offsets[dx]);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic && (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic && (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            phi += weight * grid.potential()[grid.linearIndex(ix, iy, iz)];
          }
        }
      }
      potential[p] = phi;
    }
  } else {
#if COSMOSIM_ENABLE_MPI
    const int world_size = grid.slabLayout().world_size;
    const int world_rank = grid.slabLayout().world_rank;
    double routed_mpi_wait_ms = distributed_preflight_mpi_wait_ms;
    const int local_targets_valid = first_non_finite_target_index.has_value() ? 0 : 1;
    int all_targets_valid = 0;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(&local_targets_valid, &all_targets_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD),
          "PmSolver::interpolatePotential target validation MPI_Allreduce");
    });
    if (all_targets_valid == 0) {
      const std::string local_detail = first_non_finite_target_index.has_value()
          ? " local_first_invalid_particle_index=" + std::to_string(*first_non_finite_target_index)
          : " local_targets_valid=true";
      throw std::invalid_argument(
          "PmSolver::interpolatePotential rejected non-finite target coordinates on at least one MPI rank;" +
          local_detail);
    }

    std::uint64_t exchange_epoch = 0U;
    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential exchange epoch allocation",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() { exchange_epoch = m_impl->nextDistributedExchangeEpoch(); });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      validatePmExchangeEpochConsensus(exchange_epoch, "PmSolver::interpolatePotential");
    });

    std::optional<GravityCommunicationArena::Lease> standalone_communication_lease;
    std::optional<Impl::PlaneInterpolationExchangeBuffers> exchange_storage;
    Impl::PlaneInterpolationExchangeBuffers* exchange_ptr = nullptr;
    std::uint64_t routing_metadata_capacity = 0U;
    PmRoutingCapacityModel routing_capacity{};
    std::size_t routing_buffer_limit = 0U;
    std::size_t particles_per_round = 0U;
    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential routing preflight",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
          std::pmr::memory_resource* communication_resource = m_impl->communicationResource(
              GravityCommunicationArena::Phase::kPmInterpolation,
              standalone_communication_lease);
          exchange_storage.emplace(
              communication_resource, static_cast<std::size_t>(world_size));
          exchange_ptr = &*exchange_storage;
          auto& prepared_exchange = *exchange_ptr;
          const auto add_routing_metadata = [&](const auto& values) {
            routing_metadata_capacity = checkedAddBytes(
                routing_metadata_capacity,
                core::ownedCapacityBytesForContainer(values),
                "PM potential routing metadata capacity");
          };
          add_routing_metadata(prepared_exchange.send_counts);
          add_routing_metadata(prepared_exchange.send_displs);
          add_routing_metadata(prepared_exchange.recv_counts);
          add_routing_metadata(prepared_exchange.recv_displs);
          add_routing_metadata(prepared_exchange.send_counts_bytes);
          add_routing_metadata(prepared_exchange.send_displs_bytes);
          add_routing_metadata(prepared_exchange.recv_counts_bytes);
          add_routing_metadata(prepared_exchange.recv_displs_bytes);
          add_routing_metadata(prepared_exchange.send_response_counts_bytes);
          add_routing_metadata(prepared_exchange.send_response_displs_bytes);
          add_routing_metadata(prepared_exchange.recv_response_counts_bytes);
          add_routing_metadata(prepared_exchange.recv_response_displs_bytes);
          add_routing_metadata(prepared_exchange.cursor);
          routing_capacity = modelPmRoutingCapacity(
              world_size,
              options.routing_exchange_batch_bytes,
              k_pm_routing_modeled_workspace_limit_bytes,
              routing_metadata_capacity,
              k_pm_plane_interpolation_request_wire_bytes);
          if (routing_capacity.max_send_payload_bytes >
              static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            throw std::overflow_error("PM potential routing buffer limit exceeds size_t capacity");
          }
          routing_buffer_limit =
              static_cast<std::size_t>(routing_capacity.max_send_payload_bytes);
          // The phase arena is monotonic/resettable. Reserve both bounded wire
          // maxima once so later exact-size swaps/resizes never accumulate
          // abandoned growth blocks inside the same logical lease.
          prepared_exchange.send_wire.reserve(routing_buffer_limit);
          prepared_exchange.recv_wire.reserve(routing_buffer_limit);
          particles_per_round = pmRoutingParticlesPerRound(
              routing_capacity.effective_per_peer_payload_bytes,
              k_pm_plane_interpolation_request_wire_bytes);
        });
    auto& exchange = *exchange_ptr;
    const std::uint64_t local_target_count = static_cast<std::uint64_t>(pos_x.size());
    std::uint64_t global_max_target_count = 0U;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(
              &local_target_count,
              &global_max_target_count,
              1,
              MPI_UINT64_T,
              MPI_MAX,
              MPI_COMM_WORLD),
          "PmSolver::interpolatePotential max-target MPI_Allreduce");
    });
    const std::uint64_t round_capacity = static_cast<std::uint64_t>(particles_per_round);
    const std::uint64_t round_count = global_max_target_count == 0U
        ? 0U
        : 1U + (global_max_target_count - 1U) / round_capacity;
    const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    std::uint32_t next_request_sequence = 0U;

    std::fill(potential.begin(), potential.end(), 0.0);

    const auto evaluate_plane_request = [&](const PmPlaneInterpolationRequestRecord& request) {
      if (request.plane_count == 0U || request.plane_count > 3U ||
          request.destination_rank != static_cast<std::uint32_t>(world_rank) ||
          request.exchange_epoch != exchange_epoch ||
          !std::isfinite(request.y_grid) || !std::isfinite(request.z_grid)) {
        throw std::invalid_argument(
            "PmSolver::interpolatePotential received invalid bounded potential-plane request");
      }
      const PmAxisStencil1d sy = makeAxisStencil(request.y_grid, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(request.z_grid, options.assignment_scheme);
      PmPotentialContributionRecord response{
          .source_rank = static_cast<std::uint32_t>(world_rank),
          .origin_rank = request.origin_rank,
          .request_sequence = request.request_sequence,
          .origin_particle_index = request.origin_particle_index,
          .exchange_epoch = request.exchange_epoch,
          .potential = 0.0,
      };
      for (std::size_t plane = 0; plane < static_cast<std::size_t>(request.plane_count); ++plane) {
        const std::size_t ix = static_cast<std::size_t>(request.global_ix[plane]);
        if (ix >= m_shape.nx || !grid.slabLayout().ownsGlobalX(ix) ||
            !std::isfinite(request.x_weight[plane])) {
          throw std::invalid_argument(
              "PmSolver::interpolatePotential potential-plane x ownership/weight is invalid");
        }
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic &&
              (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic &&
                (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const std::size_t local_index = decomposition_view.globalToLocalRealIndex(ix, iy, iz);
            if (!std::isfinite(grid.potential()[local_index])) {
              throw std::invalid_argument(
                  "PmSolver::interpolatePotential encountered non-finite owner-local mesh potential data");
            }
            const double weight = request.x_weight[plane] * sy.weights[dy] * sz.weights[dz];
            response.potential += weight * grid.potential()[local_index];
          }
        }
      }
      return response;
    };

    std::uint64_t logical_requests_sent = 0U;
    std::uint64_t peers_touched = 0U;
    std::uint64_t total_wire_sent = 0U;
    std::uint64_t total_wire_received = 0U;
    std::uint64_t send_high_water = 0U;
    std::uint64_t receive_high_water = 0U;
    std::uint64_t combined_high_water = 0U;
    std::uint64_t workspace_high_water = routing_metadata_capacity;

    for (std::uint64_t round = 0U; round < round_count; ++round) {
      const std::uint64_t begin_u64 = round * round_capacity;
      const std::uint64_t end_u64 = std::min(
          local_target_count,
          begin_u64 > std::numeric_limits<std::uint64_t>::max() - round_capacity
              ? local_target_count
              : begin_u64 + round_capacity);
      const std::size_t begin = static_cast<std::size_t>(begin_u64);
      const std::size_t end = static_cast<std::size_t>(end_u64);
      const std::uint32_t round_sequence_begin = next_request_sequence;

      std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
      std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
      std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
      std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
      std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
      std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
      std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
      std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
      std::fill(exchange.send_response_counts_bytes.begin(), exchange.send_response_counts_bytes.end(), 0);
      std::fill(exchange.send_response_displs_bytes.begin(), exchange.send_response_displs_bytes.end(), 0);
      std::fill(exchange.recv_response_counts_bytes.begin(), exchange.recv_response_counts_bytes.end(), 0);
      std::fill(exchange.recv_response_displs_bytes.begin(), exchange.recv_response_displs_bytes.end(), 0);
      exchange.send_wire.clear();
      exchange.recv_wire.clear();

      runPmCoordinatedPhase(
          "PmSolver::interpolatePotential bounded request preparation",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            // Pass 1 counts at most one request per particle/destination x-plane group.
            for (std::size_t p = begin; p < end; ++p) {
              if (p > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
                throw std::invalid_argument(
                    "PmSolver::interpolatePotential origin particle index exceeds routing token limit");
              }
              if (!periodic &&
                  (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
                   !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
                   !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
                continue;
              }
              const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
              const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const int destination_rank = groups[group_i].destination_rank;
                if (destination_rank == world_rank) {
                  continue;
                }
                int& count = exchange.send_counts[static_cast<std::size_t>(destination_rank)];
                if (count == std::numeric_limits<int>::max()) {
                  throw std::overflow_error("PM potential per-peer request count exceeds MPI int capacity");
                }
                ++count;
              }
            }

            std::size_t total_send_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = static_cast<std::size_t>(
                  exchange.send_counts[static_cast<std::size_t>(rank)]);
              const std::uint64_t peer_bytes = checkedBytesForCount(
                  count, k_pm_plane_interpolation_request_wire_bytes, "PM potential per-peer routing batch");
              if (peer_bytes > routing_capacity.effective_per_peer_payload_bytes) {
                throw std::runtime_error(
                    "PM potential routing exceeded effective per-peer batch high-water");
              }
              exchange.send_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_send_records,
                      k_pm_plane_interpolation_request_wire_bytes,
                      "PmSolver::interpolatePotential bounded request exchange",
                      "send record displacement",
                      rank);
              total_send_records = checkedMpiRecordTotal(
                  total_send_records,
                  count,
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolatePotential bounded request exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.send_wire,
                checkedPmWireByteCount(
                    total_send_records,
                    k_pm_plane_interpolation_request_wire_bytes,
                    "PmSolver::interpolatePotential bounded request exchange"),
                routing_buffer_limit,
                "PmSolver::interpolatePotential request send buffer");

            for (int rank = 0; rank < world_size; ++rank) {
              exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
                  exchange.send_displs[static_cast<std::size_t>(rank)]);
            }

            // Pass 2 evaluates local x planes immediately and encodes only remote
            // x-plane groups. The receiver expands the y/z stencil and returns one
            // accumulated potential response per request.
            for (std::size_t p = begin; p < end; ++p) {
              if (!periodic &&
                  (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
                   !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
                   !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
                continue;
              }
              const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
              const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
              const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;
              const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const auto& group = groups[group_i];
                PmPlaneInterpolationRequestRecord request{
                    .origin_rank = static_cast<std::uint32_t>(world_rank),
                    .destination_rank = static_cast<std::uint32_t>(group.destination_rank),
                    .request_sequence = next_request_sequence,
                    .origin_particle_index = static_cast<std::uint32_t>(p),
                    .plane_count = group.plane_count,
                    .global_ix = group.global_ix,
                    .exchange_epoch = exchange_epoch,
                    .y_grid = y,
                    .z_grid = z,
                    .x_weight = group.x_weight,
                };
                if (group.destination_rank == world_rank) {
                  const PmPotentialContributionRecord local = evaluate_plane_request(request);
                  potential[p] += local.potential;
                  continue;
                }
                if (next_request_sequence == std::numeric_limits<std::uint32_t>::max()) {
                  throw std::overflow_error(
                      "PmSolver::interpolatePotential routing sequence exceeds uint32 capacity");
                }
                auto& cursor = exchange.cursor[static_cast<std::size_t>(group.destination_rank)];
                encodePmPlaneInterpolationRequest(
                    request,
                    PmWireRecordKind::kPotentialPlaneRequest,
                    std::span<std::uint8_t>(exchange.send_wire).subspan(
                        cursor * k_pm_plane_interpolation_request_wire_bytes,
                        k_pm_plane_interpolation_request_wire_bytes));
                ++cursor;
                ++next_request_sequence;
              }
            }
          });

      const std::uint64_t request_send_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoall(
                exchange.send_counts.data(), 1, MPI_INT,
                exchange.recv_counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD),
            "PmSolver::interpolatePotential request-count MPI_Alltoall");
      });

      runPmCoordinatedPhase(
          "PmSolver::interpolatePotential bounded request receive layout",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            std::size_t total_recv_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = checkedMpiReceivedRecordCount(
                  exchange.recv_counts[static_cast<std::size_t>(rank)],
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolatePotential bounded request exchange",
                  rank);
              exchange.recv_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_recv_records,
                      k_pm_plane_interpolation_request_wire_bytes,
                      "PmSolver::interpolatePotential bounded request exchange",
                      "received record displacement",
                      rank);
              total_recv_records = checkedMpiRecordTotal(
                  total_recv_records,
                  count,
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolatePotential bounded request exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.recv_wire,
                checkedPmWireByteCount(
                    total_recv_records,
                    k_pm_plane_interpolation_request_wire_bytes,
                    "PmSolver::interpolatePotential bounded request exchange"),
                routing_buffer_limit,
                "PmSolver::interpolatePotential request receive buffer");
            checkedMpiRecordLayoutToByteLayout(
                exchange.send_counts,
                exchange.send_displs,
                k_pm_plane_interpolation_request_wire_bytes,
                exchange.send_counts_bytes,
                exchange.send_displs_bytes,
                "PmSolver::interpolatePotential bounded request send MPI_Alltoallv");
            checkedMpiRecordLayoutToByteLayout(
                exchange.recv_counts,
                exchange.recv_displs,
                k_pm_plane_interpolation_request_wire_bytes,
                exchange.recv_counts_bytes,
                exchange.recv_displs_bytes,
                "PmSolver::interpolatePotential bounded request receive MPI_Alltoallv");
          });
      const std::uint64_t request_recv_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());

      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.send_wire),
                exchange.send_counts_bytes.data(),
                exchange.send_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.recv_wire),
                exchange.recv_counts_bytes.data(),
                exchange.recv_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::interpolatePotential request MPI_Alltoallv");
      });

      std::string local_potential_request_validation_error;
      try {
        // Compact responses in-place into the receive-request buffer. Since a
        // response (48 B) is smaller than a request (96 B), forward compaction
        // cannot overwrite a request that has not yet been decoded.
        for (int source_rank = 0; source_rank < world_size; ++source_rank) {
          const int request_begin = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
          const int request_count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
          exchange.send_response_counts_bytes[static_cast<std::size_t>(source_rank)] =
              checkedMpiByteCount(
                  checkedBytesForCount(
                      static_cast<std::size_t>(request_count),
                      k_pm_potential_response_wire_bytes,
                      "PM potential response sender count"),
                  "PM potential response sender count");
          exchange.send_response_displs_bytes[static_cast<std::size_t>(source_rank)] =
              checkedMpiByteCount(
                  checkedBytesForCount(
                      static_cast<std::size_t>(request_begin),
                      k_pm_potential_response_wire_bytes,
                      "PM potential response sender displacement"),
                  "PM potential response sender displacement");
          std::optional<std::uint32_t> previous_sequence;
          for (int i = 0; i < request_count; ++i) {
            const std::size_t request_index = static_cast<std::size_t>(request_begin + i);
            const PmPlaneInterpolationRequestRecord request = decodePmPlaneInterpolationRequest(
                std::span<const std::uint8_t>(exchange.recv_wire).subspan(
                    request_index * k_pm_plane_interpolation_request_wire_bytes,
                    k_pm_plane_interpolation_request_wire_bytes),
                PmWireRecordKind::kPotentialPlaneRequest,
                "PmSolver::interpolatePotential potential-plane request");
            if (request.origin_rank != static_cast<std::uint32_t>(source_rank) ||
                request.destination_rank != static_cast<std::uint32_t>(world_rank) ||
                request.exchange_epoch != exchange_epoch) {
              throw std::invalid_argument(
                  "PM potential-plane routing metadata does not match MPI sender/receiver segment");
            }
            if (previous_sequence.has_value() && request.request_sequence <= *previous_sequence) {
              throw std::invalid_argument(
                  "PM potential-plane request sequence is non-monotonic within sender segment");
            }
            previous_sequence = request.request_sequence;
            const PmPotentialContributionRecord response = evaluate_plane_request(request);
            const std::size_t response_index =
                static_cast<std::size_t>(exchange.send_response_displs_bytes[static_cast<std::size_t>(source_rank)]) /
                    k_pm_potential_response_wire_bytes +
                static_cast<std::size_t>(i);
            encodePmPotentialResponse(
                response,
                std::span<std::uint8_t>(exchange.recv_wire).subspan(
                    response_index * k_pm_potential_response_wire_bytes,
                    k_pm_potential_response_wire_bytes));
          }
        }
      } catch (const std::exception& error) {
        local_potential_request_validation_error = error.what();
      } catch (...) {
        local_potential_request_validation_error =
            "non-standard exception while consuming bounded potential-plane requests";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::interpolatePotential bounded request receive",
          world_rank,
          world_size,
          local_potential_request_validation_error,
          routed_mpi_wait_ms);

      const std::size_t total_received_requests = exchange.recv_wire.size() /
          k_pm_plane_interpolation_request_wire_bytes;
      resizePmWireBufferBounded(
          exchange.recv_wire,
          checkedPmWireByteCount(
              total_received_requests,
              k_pm_potential_response_wire_bytes,
              "PmSolver::interpolatePotential compact response send"),
          routing_buffer_limit,
          "PmSolver::interpolatePotential response send buffer");

      const std::size_t total_sent_requests = exchange.send_wire.size() /
          k_pm_plane_interpolation_request_wire_bytes;
      for (int rank = 0; rank < world_size; ++rank) {
        exchange.recv_response_counts_bytes[static_cast<std::size_t>(rank)] =
            checkedMpiByteCount(
                checkedBytesForCount(
                    static_cast<std::size_t>(exchange.send_counts[static_cast<std::size_t>(rank)]),
                    k_pm_potential_response_wire_bytes,
                    "PM potential response receive count"),
                "PM potential response receive count");
        exchange.recv_response_displs_bytes[static_cast<std::size_t>(rank)] =
            checkedMpiByteCount(
                checkedBytesForCount(
                    static_cast<std::size_t>(exchange.send_displs[static_cast<std::size_t>(rank)]),
                    k_pm_potential_response_wire_bytes,
                    "PM potential response receive displacement"),
                "PM potential response receive displacement");
      }
      resizePmWireBufferBounded(
          exchange.send_wire,
          checkedPmWireByteCount(
              total_sent_requests,
              k_pm_potential_response_wire_bytes,
              "PmSolver::interpolatePotential compact response receive"),
          routing_buffer_limit,
          "PmSolver::interpolatePotential response receive buffer");

      const std::uint64_t response_send_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());
      const std::uint64_t response_recv_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.recv_wire),
                exchange.send_response_counts_bytes.data(),
                exchange.send_response_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.send_wire),
                exchange.recv_response_counts_bytes.data(),
                exchange.recv_response_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::interpolatePotential response MPI_Alltoallv");
      });

      std::string local_potential_response_validation_error;
      try {
        // Re-generate the compact request stream instead of retaining an O(N)
        // request registry. This provides exact response identity validation with
        // only O(world_size) cursor state.
        for (int rank = 0; rank < world_size; ++rank) {
          exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
              exchange.send_displs[static_cast<std::size_t>(rank)]);
        }
        std::uint32_t expected_sequence = round_sequence_begin;
        for (std::size_t p = begin; p < end; ++p) {
          if (!periodic &&
              (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
               !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
               !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
            continue;
          }
          const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
          const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
          std::array<PmXPlaneGroup, 3> groups{};
          const std::size_t group_count = makePmXPlaneGroups(
              sx, periodic, m_shape, decomposition_view, groups);
          for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
            const int sender_rank = groups[group_i].destination_rank;
            if (sender_rank == world_rank) {
              continue;
            }
            auto& cursor = exchange.cursor[static_cast<std::size_t>(sender_rank)];
            const std::size_t segment_end = static_cast<std::size_t>(
                exchange.send_displs[static_cast<std::size_t>(sender_rank)] +
                exchange.send_counts[static_cast<std::size_t>(sender_rank)]);
            if (cursor >= segment_end) {
              throw std::invalid_argument("PM potential response segment ended before issued requests");
            }
            const PmPotentialContributionRecord response = decodePmPotentialResponse(
                std::span<const std::uint8_t>(exchange.send_wire).subspan(
                    cursor * k_pm_potential_response_wire_bytes,
                    k_pm_potential_response_wire_bytes));
            if (response.source_rank != static_cast<std::uint32_t>(sender_rank) ||
                response.origin_rank != static_cast<std::uint32_t>(world_rank) ||
                response.exchange_epoch != exchange_epoch ||
                response.request_sequence != expected_sequence ||
                response.origin_particle_index != static_cast<std::uint32_t>(p) ||
                !std::isfinite(response.potential)) {
              throw std::invalid_argument(
                  "PM potential response does not match the regenerated bounded request stream");
            }
            potential[p] += response.potential;
            ++cursor;
            ++expected_sequence;
          }
        }
        if (expected_sequence != next_request_sequence) {
          throw std::logic_error("PM potential response validation sequence accounting drifted");
        }
        for (int rank = 0; rank < world_size; ++rank) {
          const std::size_t expected_end = static_cast<std::size_t>(
              exchange.send_displs[static_cast<std::size_t>(rank)] +
              exchange.send_counts[static_cast<std::size_t>(rank)]);
          if (exchange.cursor[static_cast<std::size_t>(rank)] != expected_end) {
            throw std::invalid_argument("PM potential response segment contains unconsumed records");
          }
        }
      } catch (const std::exception& error) {
        local_potential_response_validation_error = error.what();
      } catch (...) {
        local_potential_response_validation_error =
            "non-standard exception while validating bounded potential responses";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::interpolatePotential bounded response receive",
          world_rank,
          world_size,
          local_potential_response_validation_error,
          routed_mpi_wait_ms);

      logical_requests_sent += static_cast<std::uint64_t>(total_sent_requests);
      peers_touched += static_cast<std::uint64_t>(std::count_if(
          exchange.send_counts.begin(), exchange.send_counts.end(), [](int count) { return count > 0; }));
      total_wire_sent += request_send_bytes + response_send_bytes;
      total_wire_received += request_recv_bytes + response_recv_bytes;
      const std::uint64_t send_wire_capacity =
          static_cast<std::uint64_t>(exchange.send_wire.capacity());
      const std::uint64_t recv_wire_capacity =
          static_cast<std::uint64_t>(exchange.recv_wire.capacity());
      // The two physical buffers swap direction for the response phase, so
      // each directional high-water must consider both retained capacities.
      const std::uint64_t directional_capacity =
          std::max(send_wire_capacity, recv_wire_capacity);
      const std::uint64_t combined_capacity = checkedAddBytes(
          send_wire_capacity,
          recv_wire_capacity,
          "PM potential combined routing buffer capacity");
      const std::uint64_t workspace_capacity = checkedAddBytes(
          routing_metadata_capacity,
          combined_capacity,
          "PM potential routing workspace capacity");
      if (workspace_capacity > k_pm_routing_workspace_target_bytes) {
        throw std::logic_error("PM potential routing workspace exceeded the M1A aggregate target");
      }
      send_high_water = std::max(send_high_water, directional_capacity);
      receive_high_water = std::max(receive_high_water, directional_capacity);
      combined_high_water = std::max(combined_high_water, combined_capacity);
      workspace_high_water = std::max(workspace_high_water, workspace_capacity);
      exchange.workspace_high_water_bytes = std::max(
          exchange.workspace_high_water_bytes, workspace_capacity);
      m_impl->recordCommunicationHighWater(
          GravityCommunicationArena::Phase::kPmInterpolation,
          workspace_capacity,
          standalone_communication_lease);
    }

    if (profile != nullptr) {
      profile->routed_potential_requests += logical_requests_sent;
      profile->routed_potential_peer_count += peers_touched;
      profile->routed_mpi_bytes_sent += total_wire_sent;
      profile->routed_mpi_bytes_received += total_wire_received;
      profile->routed_send_buffer_high_water_bytes = std::max(
          profile->routed_send_buffer_high_water_bytes, send_high_water);
      profile->routed_receive_buffer_high_water_bytes = std::max(
          profile->routed_receive_buffer_high_water_bytes, receive_high_water);
      profile->routed_combined_buffer_high_water_bytes = std::max(
          profile->routed_combined_buffer_high_water_bytes, combined_high_water);
      profile->routed_workspace_high_water_bytes = std::max(
          profile->routed_workspace_high_water_bytes, workspace_high_water);
      profile->routed_mpi_wait_ms += routed_mpi_wait_ms;
      profile->bytes_moved += total_wire_sent + total_wire_received;
    }
#endif
  }

  const auto stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->interpolate_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    profile->bytes_moved += bytesForParticles(pos_x.size());
  }
}

void PmSolver::solveForParticles(
    PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::span<double> accel_x,
    std::span<double> accel_y,
    std::span<double> accel_z,
    const PmSolveOptions& options,
    PmProfileEvent* profile) {
#if COSMOSIM_ENABLE_MPI
  double solve_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kSolveForParticles,
        rank_local_serial_layout,
        m_shape,
        options,
        solve_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::solveForParticles API preflight",
        mpi_world_rank,
        mpi_world_size,
        solve_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny ||
              grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in solveForParticles");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::solveForParticles requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::solveForParticles slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() ||
              pos_x.size() != mass.size() || pos_x.size() != accel_x.size() ||
              pos_x.size() != accel_y.size() || pos_x.size() != accel_z.size()) {
            throw std::invalid_argument(
                "Particle coordinate/mass/acceleration spans must match in solveForParticles");
          }
#if !COSMOSIM_ENABLE_CUDA
          if (options.execution_policy == core::ExecutionPolicy::kCuda) {
            throw std::runtime_error(
                "PM solve requested execution_policy=cuda, but this build has COSMOSIM_ENABLE_CUDA=OFF");
          }
#else
          if (options.execution_policy == core::ExecutionPolicy::kCuda &&
              !rank_local_serial_layout) {
            throw std::invalid_argument(
                "Distributed PM solve does not support execution_policy=cuda");
          }
#endif
        });
    if (profile != nullptr) {
      profile->routed_mpi_wait_ms += solve_preflight_mpi_wait_ms;
    }
  }
#endif
  // Backend-independent API contract. Keep this outside the MPI branch so a
  // single-rank CUDA solve cannot bypass the source/target shape checks that
  // the host and distributed paths rely on. Multi-rank callers have already
  // passed the coordinated equivalent above, making these checks deterministic.
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny ||
      grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in solveForParticles");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::solveForParticles requires a valid PM slab layout");
  }
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() ||
      pos_x.size() != mass.size() || pos_x.size() != accel_x.size() ||
      pos_x.size() != accel_y.size() || pos_x.size() != accel_z.size()) {
    throw std::invalid_argument(
        "Particle coordinate/mass/acceleration spans must match in solveForParticles");
  }
#if !COSMOSIM_ENABLE_MPI
  if (grid.slabLayout().world_size != 1) {
    throw std::invalid_argument(
        "PmSolver::solveForParticles distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
#endif
  if (grid.slabLayout().world_size == 1) {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::solveForParticles");
  }
#if !COSMOSIM_ENABLE_CUDA
  if (options.execution_policy == core::ExecutionPolicy::kCuda) {
    throw std::runtime_error(
        "PM solve requested execution_policy=cuda, but this build has COSMOSIM_ENABLE_CUDA=OFF");
  }
#endif

  if (options.execution_policy == core::ExecutionPolicy::kCuda) {
#if COSMOSIM_ENABLE_CUDA
    m_impl->solveCudaStagingBackend(
        *this, grid, pos_x, pos_y, pos_z, mass, accel_x, accel_y, accel_z,
        options, profile);
    return;
#else
    throw std::runtime_error("PM solve requested execution_policy=cuda, but this build has COSMOSIM_ENABLE_CUDA=OFF");
#endif
  }

  assignDensity(grid, pos_x, pos_y, pos_z, mass, options, profile);
  if (options.boundary_condition == PmBoundaryCondition::kPeriodic) {
    solvePoissonPeriodic(grid, options, profile);
  } else {
    solvePoissonIsolatedOpen(grid, options, profile);
  }
  interpolateForces(grid, pos_x, pos_y, pos_z, accel_x, accel_y, accel_z, options, profile);
}

bool PmSolver::fftBackendAvailable() {
#if COSMOSIM_ENABLE_FFTW
  return true;
#else
  return false;
#endif
}

bool PmSolver::cudaBackendAvailable() {
#if COSMOSIM_ENABLE_CUDA
  int device_count = 0;
  return cudaGetDeviceCount(&device_count) == cudaSuccess && device_count > 0;
#else
  return false;
#endif
}

std::string PmSolver::fftBackendName() {
#if COSMOSIM_ENABLE_FFTW
  return "fftw";
#else
  return "naive_dft";
#endif
}

std::size_t PmSolver::cachedPlanCount() const {
  return m_impl->planCount();
}

std::size_t PmSolver::planBuildCount() const {
  return m_impl->planBuildCount();
}

bool treePmSupportedByBuild() {
  return pmBackendProductionReady();
}

PmBackendCapability pmBackendCapability() noexcept {
#if COSMOSIM_ENABLE_FFTW
  return PmBackendCapability::kProductionFftw;
#else
  return PmBackendCapability::kDiagnosticNaiveDft;
#endif
}

std::string_view pmBackendCapabilityName(PmBackendCapability capability) noexcept {
  switch (capability) {
    case PmBackendCapability::kUnavailable:
      return "unavailable";
    case PmBackendCapability::kDiagnosticNaiveDft:
      return "diagnostic_naive_dft";
    case PmBackendCapability::kProductionFftw:
      return "production_fftw";
    case PmBackendCapability::kProductionCudaFft:
      return "production_cuda_fft";
  }
  return "unknown";
}

bool pmBackendProductionReady() noexcept {
  const PmBackendCapability capability = pmBackendCapability();
  return capability == PmBackendCapability::kProductionFftw ||
      capability == PmBackendCapability::kProductionCudaFft;
}

PmBackendArchitecture pmBackendArchitecture() noexcept {
  PmBackendArchitecture architecture;
  architecture.capability = pmBackendCapability();
#if COSMOSIM_ENABLE_FFTW
  architecture.host_fft = true;
#endif
#if COSMOSIM_ENABLE_CUDA
  // The present CUDA path performs assignment/interpolation on device but
  // returns through the host FFT solver. Keep the missing resident stages
  // explicit so provenance/capability code cannot over-promote it.
  architecture.device_assignment = true;
  architecture.device_interpolation = true;
  architecture.persistent_device_buffers = true;
  // The staging backend now retains its stream and device workspaces across
  // solves. Full persistent residency remains false because density/force still
  // cross the host boundary for FFT/Poisson/gradient.
  architecture.persistent_residency = false;
#endif
  return architecture;
}

std::string_view pmDecompositionTopologyName(PmDecompositionTopology topology) noexcept {
  switch (topology) {
    case PmDecompositionTopology::kSerial:
      return "serial";
    case PmDecompositionTopology::kXSlab:
      return "x_slab";
    case PmDecompositionTopology::kXSlabTransposedSpectral:
      return "x_slab_transposed_spectral";
    case PmDecompositionTopology::kPencil2D:
      return "pencil_2d";
  }
  return "unknown";
}

PmDecompositionView::PmDecompositionView(
    PmGridShape shape,
    parallel::PmSlabLayout layout,
    core::PmDecompositionMode mode)
    : m_shape(shape),
      m_layout(std::move(layout)),
      m_mode(mode),
      m_descriptor(describePmDecomposition(m_shape, m_layout, m_mode)),
      m_world_rank(m_layout.world_rank),
      m_world_size(m_layout.world_size) {}

PmDecompositionView::PmDecompositionView(
    PmGridShape shape,
    PmDecompositionDescriptor descriptor,
    int world_rank,
    int world_size)
    : m_shape(shape),
      m_descriptor(std::move(descriptor)),
      m_world_rank(world_rank),
      m_world_size(world_size) {
  if (!m_shape.isValid() || m_world_size <= 0 || m_world_rank < 0 ||
      m_world_rank >= m_world_size) {
    throw std::invalid_argument("PM decomposition view requires valid mesh/rank metadata");
  }
  if (m_descriptor.process_grid_x == 0U || m_descriptor.process_grid_y == 0U ||
      m_descriptor.process_grid_x * m_descriptor.process_grid_y !=
          static_cast<std::size_t>(m_world_size)) {
    throw std::invalid_argument("PM decomposition process grid does not match world size");
  }
}

const PmDecompositionDescriptor& PmDecompositionView::descriptor() const noexcept {
  return m_descriptor;
}

PmOwnershipExtent3D PmDecompositionView::realExtentForRank(int rank) const {
  if (rank < 0 || rank >= m_world_size) {
    throw std::out_of_range("PM real-space extent rank is outside decomposition world");
  }
  const auto block_range = [](std::size_t count, std::size_t parts, std::size_t coord) {
    const std::size_t base = count / parts;
    const std::size_t remainder = count % parts;
    const std::size_t begin = coord * base + std::min(coord, remainder);
    const std::size_t extent = base + (coord < remainder ? 1U : 0U);
    return std::pair<std::size_t, std::size_t>{begin, extent};
  };
  if (m_descriptor.topology == PmDecompositionTopology::kPencil2D) {
    const std::size_t rank_u = static_cast<std::size_t>(rank);
    const std::size_t px = rank_u % m_descriptor.process_grid_x;
    const std::size_t py = rank_u / m_descriptor.process_grid_x;
    const auto [x_begin, x_count] = block_range(
        m_shape.nx, m_descriptor.process_grid_x, px);
    const auto [y_begin, y_count] = block_range(
        m_shape.ny, m_descriptor.process_grid_y, py);
    return PmOwnershipExtent3D{
        .x_begin = x_begin,
        .x_count = x_count,
        .y_begin = y_begin,
        .y_count = y_count,
        .z_begin = 0U,
        .z_count = m_shape.nz,
    };
  }
  const parallel::PmSlabRange owned =
      parallel::pmOwnedXRangeForRank(m_shape.nx, m_world_size, rank);
  return PmOwnershipExtent3D{
      .x_begin = owned.begin_x,
      .x_count = owned.extentX(),
      .y_begin = 0U,
      .y_count = m_shape.ny,
      .z_begin = 0U,
      .z_count = m_shape.nz,
  };
}

int PmDecompositionView::realOwnerRank(PmGlobalCell global_cell) const {
  if (global_cell.x >= m_shape.nx || global_cell.y >= m_shape.ny ||
      global_cell.z >= m_shape.nz) {
    throw std::out_of_range("PM real-space global cell is outside the mesh");
  }
  if (m_descriptor.topology == PmDecompositionTopology::kPencil2D) {
    const auto owner_coord = [](std::size_t index, std::size_t count, std::size_t parts) {
      const std::size_t base = count / parts;
      const std::size_t remainder = count % parts;
      const std::size_t wide = (base + 1U) * remainder;
      if (index < wide) {
        return index / (base + 1U);
      }
      return remainder + (index - wide) / std::max<std::size_t>(base, 1U);
    };
    const std::size_t px = owner_coord(
        global_cell.x, m_shape.nx, m_descriptor.process_grid_x);
    const std::size_t py = owner_coord(
        global_cell.y, m_shape.ny, m_descriptor.process_grid_y);
    return static_cast<int>(py * m_descriptor.process_grid_x + px);
  }
  return parallel::pmOwnerRankForGlobalX(
      m_shape.nx, m_world_size, global_cell.x);
}

int PmDecompositionView::spectralOwnerRank(PmGlobalCell global_cell) const {
  if (global_cell.x >= m_shape.nx || global_cell.y >= m_shape.ny ||
      global_cell.z >= m_shape.nz) {
    throw std::out_of_range("PM spectral global cell is outside the mesh");
  }
  if (m_descriptor.topology == PmDecompositionTopology::kXSlabTransposedSpectral) {
    return parallel::pmOwnerRankForGlobalX(
        m_shape.ny, m_world_size, global_cell.y);
  }
  return realOwnerRank(global_cell);
}

bool PmDecompositionView::ownsRealCell(
    std::size_t global_x, std::size_t global_y, std::size_t global_z) const noexcept {
  if (global_x >= m_shape.nx || global_y >= m_shape.ny || global_z >= m_shape.nz) {
    return false;
  }
  try {
    return realOwnerRank(PmGlobalCell{global_x, global_y, global_z}) == m_world_rank;
  } catch (...) {
    return false;
  }
}

std::size_t PmDecompositionView::globalToLocalRealIndex(
    std::size_t global_x, std::size_t global_y, std::size_t global_z) const {
  if (!ownsRealCell(global_x, global_y, global_z)) {
    throw std::out_of_range("PM real-space global index is not owned by this backend view");
  }
  const PmOwnershipExtent3D extent = realExtentForRank(m_world_rank);
  return core::checkedSizeAdd(
      core::checkedSizeProduct3(
          global_x - extent.x_begin, extent.y_count, extent.z_count,
          "PM real local index"),
      core::checkedSizeAdd(
          core::checkedSizeProduct3(
              global_y - extent.y_begin, extent.z_count, 1U,
              "PM real local index"),
          global_z - extent.z_begin,
          "PM real local index"),
      "PM real local index");
}

bool PmDecompositionView::ownsSpectralCell(
    std::size_t global_x, std::size_t global_y, std::size_t global_z) const noexcept {
  const auto& e = m_descriptor.spectral_extent;
  return e.x_count != 0U && e.y_count != 0U && e.z_count != 0U &&
      global_x >= e.x_begin && global_x - e.x_begin < e.x_count &&
      global_y >= e.y_begin && global_y - e.y_begin < e.y_count &&
      global_z >= e.z_begin && global_z - e.z_begin < e.z_count;
}

std::size_t PmDecompositionView::globalToLocalSpectralIndex(
    std::size_t global_x, std::size_t global_y, std::size_t global_z) const {
  if (!ownsSpectralCell(global_x, global_y, global_z)) {
    throw std::out_of_range(
        "PM spectral global index is not owned by this backend view");
  }
  const auto& e = m_descriptor.spectral_extent;
  const std::size_t lx = global_x - e.x_begin;
  const std::size_t ly = global_y - e.y_begin;
  const std::size_t lz = global_z - e.z_begin;
  return core::checkedSizeAdd(
      core::checkedSizeProduct3(lx, e.y_count, e.z_count,
                                "PM spectral local index"),
      core::checkedSizeAdd(
          core::checkedSizeProduct3(ly, e.z_count, 1U, "PM spectral local index"),
          lz, "PM spectral local index"),
      "PM spectral local index");
}

PmDecompositionDescriptor describePmDecomposition(
    PmGridShape shape,
    const parallel::PmSlabLayout& layout,
    core::PmDecompositionMode mode) {
  if (!shape.isValid() || !layout.isValid() ||
      layout.global_nx != shape.nx || layout.global_ny != shape.ny || layout.global_nz != shape.nz) {
    throw std::invalid_argument("PM decomposition descriptor requires a valid layout matching the mesh shape");
  }
  PmDecompositionDescriptor descriptor;
  descriptor.process_grid_x = static_cast<std::size_t>(std::max(layout.world_size, 1));
  descriptor.process_grid_y = 1U;
  descriptor.real_extent = PmOwnershipExtent3D{
      .x_begin = layout.owned_x.begin_x,
      .x_count = layout.local_nx(),
      .y_begin = 0U,
      .y_count = shape.ny,
      .z_begin = 0U,
      .z_count = shape.nz,
  };
  // Spectral ownership is backend-specific. The current transposed-slab FFTW
  // implementation fills it internally; zero counts here deliberately mean
  // "owned by the selected backend contract", not a fabricated pencil range.
  if (layout.ownsFullDomain()) {
    descriptor.topology = PmDecompositionTopology::kSerial;
    descriptor.spectral_extent = descriptor.real_extent;
  } else if (mode == core::PmDecompositionMode::kPencil) {
    descriptor.topology = PmDecompositionTopology::kXSlabTransposedSpectral;
  } else {
    descriptor.topology = PmDecompositionTopology::kXSlab;
  }
  return descriptor;
}

void requireTreePmSupportOrThrow(
    core::GravitySolver gravity_solver,
    bool allow_diagnostic_naive_dft) {
  if (gravity_solver != core::GravitySolver::kTreePm) {
    return;
  }
  const PmBackendCapability capability = pmBackendCapability();
  if (pmBackendProductionReady()) {
    return;
  }
  if (capability == PmBackendCapability::kDiagnosticNaiveDft && allow_diagnostic_naive_dft) {
    return;
  }
  throw std::runtime_error(
      "TreePM production run requires a production FFT backend; this build exposes only " +
      std::string(pmBackendCapabilityName(capability)) +
      ". Enable FFTW or use the explicit diagnostic naive-DFT override for tiny validation cases.");
}

}  // namespace cosmosim::gravity
