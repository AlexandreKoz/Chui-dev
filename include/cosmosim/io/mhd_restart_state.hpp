#pragma once

#include <cstdint>
#include <string>
#include <string_view>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::io {

// Exact-restart magnetic authority group of the MHD checkpoint schema
// cosmosim_restart_v26 (CHUI-MHD-REF v1.0 section 13.2 with compatibility
// amendment R1: v25 is the existing nonmagnetic dual-energy schema and is
// never reused for magnetic fields).
//
// Group /state/mhd (group attribute schema_version = 1):
//   patch_id                 uint64[P]
//   geometry_generation      uint64[P]
//   refinement_level         uint8[P]
//   face_dims_xyz            uint32[P,3]   (cell counts nx, ny, nz)
//   patch_origin_comoving    float64[P,3]
//   cell_width_comoving      float64[P,3]
//   face_offset_x/_y/_z      uint64[P+1]   (exclusive prefix sums, [0] = 0)
//   face_flux_x/_y/_z_comoving float64[sum]
// Attributes: reference_contract (string), coefficient_table_sha256 (string),
// divergence_max_norm (float64), scale_factor (float64),
// magnetic_units = "code_comoving_lorentz_heaviside", face_flux_contract.
//
// STATUS: codec only. writeRestartCheckpointHdf5 still fails closed for a
// magnetized state; wiring this group into the v26 writer/reader pipeline
// (collective MPI publication, integrity readback, config-hash admission) is
// recorded as open work in docs/architecture/chui_mhd_implementation_pass01_status.md.

inline constexpr std::string_view k_mhd_restart_schema_name = "cosmosim_restart_v26";
inline constexpr std::uint32_t k_mhd_restart_schema_version = 26U;
inline constexpr std::uint32_t k_mhd_restart_group_schema_version = 1U;
inline constexpr std::string_view k_mhd_restart_magnetic_units = "code_comoving_lorentz_heaviside";
// Bounded HDF5 staging per rank (CHUI-MHD-PERF v1.0 section 12, MP08).
inline constexpr std::uint64_t k_mhd_restart_chunk_bytes = 8ULL * 1024ULL * 1024ULL;

struct MhdRestartGroupMetadata {
  std::string reference_contract;
  std::string coefficient_table_sha256;
  double divergence_max_norm = 0.0;
  double scale_factor = 1.0;
};

// What a reader requires before admitting the magnetic state.
struct MhdRestartExpectation {
  std::string reference_contract;
  std::string coefficient_table_sha256;  // empty = provider without a table
  double scale_factor = 1.0;
  double scale_factor_relative_tolerance = 1.0e-12;
  double divergence_tolerance = 1.0e-11;
};

// Write /state/mhd into an open HDF5 file (hid_t passed as int64 to keep this
// header free of hdf5.h). Validates the storage first (MagneticFaceSoa::
// isConsistent, finite values); writes flux lanes in <= 8 MiB hyperslabs
// directly from the canonical lanes (no full mirror copy).
void writeMhdRestartGroup(std::int64_t hdf5_file,
                          const core::MagneticFaceSoa& faces,
                          const MhdRestartGroupMetadata& metadata);

// Read and FULLY validate /state/mhd into a fresh MagneticFaceSoa (version,
// contract, units, table SHA, scale factor, shapes, offsets P+1 monotone from
// 0, checked counts, finite values, unique patch ids, recomputed per-patch
// discrete divergence <= tolerance). Throws std::runtime_error with the prefix
// MHD_RESTART_SCHEMA_MISMATCH; the caller swaps the result into
// SimulationState only after every other restart group has also validated.
[[nodiscard]] core::MagneticFaceSoa readMhdRestartGroup(std::int64_t hdf5_file,
                                                        const MhdRestartExpectation& expected,
                                                        MhdRestartGroupMetadata* metadata_out);

}  // namespace cosmosim::io
