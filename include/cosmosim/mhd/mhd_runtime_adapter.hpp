#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_state.hpp"
#include "cosmosim/mhd/mhd_step_operator.hpp"
#include "cosmosim/mhd/mhd_transport.hpp"

namespace cosmosim::mhd {

// Adapter between core-owned plain state (core::SimulationState,
// core::MagneticFaceSoa) and the MHD solver types. core never includes mhd;
// the workflow layer composes them through this adapter. Nothing here runs a
// physics step; production dispatch is still rejected at capability admission
// (runtime capability mhd_production_dispatch).

// One resolved patch binding between core storage and the MHD solver.
struct MhdAdapterPatchBinding {
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::int32_t level = 0;
  MhdPatchGeometry geometry{};
};

// Translate the typed config into the numerics policy. The only place core
// config enums are mapped; unsupported selections throw
// MHD_UNSUPPORTED_PHYSICS instead of degrading silently.
[[nodiscard]] MhdNumericsPolicy mhdNumericsPolicyFromConfig(const core::MhdConfig& config);

[[nodiscard]] MhdPatchBoundary mhdBoundaryFromCondition(core::BoundaryCondition condition);

// Build a complete, freshly allocated MagneticFaceSoa for an ordered set of
// patch bindings (deterministic order: as supplied; callers pass the AMR
// owner's stable order). Checked arithmetic on every count and offset;
// duplicate patch ids, invalid geometry or count overflow throw. All faces
// are zero; the caller seeds them (initial conditions) or transfers them
// (regrid/restart) and then swaps the result into SimulationState, so old and
// new layouts coexist only for the duration of the transaction.
[[nodiscard]] core::MagneticFaceSoa buildMagneticFaceLayout(
    std::span<const MhdAdapterPatchBinding> bindings,
    std::string coefficient_table_sha256);

// Deep validation of a core::MagneticFaceSoa: record/geometry validity,
// unique (patch_id), offset tiling with checked arithmetic, exact lane sizes
// and finite face values. Throws MHD_INVALID_FIELD / MHD_STALE_GEOMETRY.
void validateMagneticFaceSoa(const core::MagneticFaceSoa& storage);

// Locate a patch record by (patch_id, geometry_generation); returns
// storage.patches.size() when absent.
[[nodiscard]] std::size_t findMagneticFacePatch(const core::MagneticFaceSoa& storage,
                                                std::uint64_t patch_id,
                                                std::uint64_t geometry_generation);

[[nodiscard]] MhdPatchGeometry magneticPatchGeometry(const core::MagneticFacePatchRecord& record);

[[nodiscard]] MhdFaceFluxView makeMagneticFaceView(const core::MagneticFaceSoa& storage,
                                                   const MhdAdapterPatchBinding& binding,
                                                   std::size_t patch_index);

[[nodiscard]] MhdFaceFluxMutableView makeMutableMagneticFaceView(
    core::MagneticFaceSoa& storage,
    const MhdAdapterPatchBinding& binding,
    std::size_t patch_index);

// Gather one patch of canonical gas rows into the MHD comoving conserved
// state. `dense_row_by_patch_cell[cell]` maps the patch-local cell index
// (i + nx (j + ny k)) to a dense SimulationState gas row (from the stable
// gas-cell identity map, never a persistent dense-row assumption). The total
// energy is the single MHD authority
//   E_c = rho_c e + rho_c |u|^2 / 2 + |B_c|^2 / (2 a),
// with B_c derived from the authoritative face fluxes at scale factor a.
void loadMhdGasFromState(const core::SimulationState& state,
                         const MhdPatchGeometry& geometry,
                         std::span<const std::uint32_t> dense_row_by_patch_cell,
                         const MhdFaceFluxView& faces,
                         double scale_factor,
                         MhdGasStateSoa& gas);

// Scatter an accepted MHD patch state back to canonical gas rows. Specific
// internal energy is recovered as e = (E_c - K_c - |B_c|^2/(2a)) / rho_c with
// B_c from the accepted faces and the same EOS gamma as the solver. There is
// NO floor: a non-positive density or thermal energy throws
// MHD_NONPOSITIVE_STATE before any row is written. The pressure cache is
// refreshed; temperature/sound-speed caches must be refreshed by the existing
// gas thermodynamic refresh owner.
void storeMhdGasToState(core::SimulationState& state,
                        const MhdPatchGeometry& geometry,
                        std::span<const std::uint32_t> dense_row_by_patch_cell,
                        const MhdGasStateSoa& gas,
                        const MhdFaceFluxView& faces,
                        double scale_factor,
                        double gamma_adiabatic);

// Comoving magnetic energy of one stored patch: sum |B_c|^2/(2a) V.
[[nodiscard]] double patchMagneticEnergy(const core::MagneticFaceSoa& storage,
                                         std::size_t patch_index,
                                         double scale_factor);

// Construct and admit the configured coefficient provider:
//   constant -> ConstantMhdCoefficientProvider (cgs, signed Hall);
//   table    -> TableMhdCoefficientProvider with the configured path,
//               composition/grain ids and optional pinned SHA-256;
//   chemistry/none -> MHD_UNSUPPORTED_PHYSICS / MHD_BAD_COEFFICIENT.
// validateAdmission() has been called on the returned provider.
[[nodiscard]] std::unique_ptr<MhdCoefficientProvider> makeConfiguredCoefficientProvider(
    const core::MhdConfig& config);

}  // namespace cosmosim::mhd
