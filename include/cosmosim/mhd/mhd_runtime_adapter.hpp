#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_state.hpp"
#include "cosmosim/mhd/mhd_step_operator.hpp"
#include "cosmosim/mhd/mhd_transport.hpp"

namespace cosmosim::core {
class UnitSystem;
}  // namespace cosmosim::core

namespace cosmosim::mhd {

// One resolved patch binding between core storage and the MHD solver. It is
// rebuilt whenever the (patch_id, geometry_generation, dims) key changes.
struct MhdAdapterPatchBinding {
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::int32_t level = 0;
  MhdPatchGeometry geometry{};
};

// Build the numerics policy from the typed config. This is the only place the
// core config enums are translated; the solver itself has no config dependency.
[[nodiscard]] MhdNumericsPolicy mhdNumericsPolicyFromConfig(const core::MhdConfig& config);

// Map a core boundary condition to the MHD physical boundary class.
[[nodiscard]] MhdPatchBoundary mhdBoundaryFromCondition(core::BoundaryCondition condition);

// Ensure a persistent face-flux patch exists for the binding. The first call
// for a key allocates the concatenated lanes; a geometry/dimension change
// reallocates from zero. Returns the patch index.
[[nodiscard]] std::size_t ensureMagneticFacePatch(core::MagneticFaceSoa& storage,
                                                  const MhdAdapterPatchBinding& binding);

// Locate an existing patch record; returns size when absent.
[[nodiscard]] std::size_t findMagneticFacePatch(const core::MagneticFaceSoa& storage,
                                                std::uint64_t patch_id,
                                                std::uint64_t geometry_generation);

[[nodiscard]] MhdFaceFluxView makeMagneticFaceView(
    const core::MagneticFaceSoa& storage,
    const MhdAdapterPatchBinding& binding,
    std::size_t patch_index);

[[nodiscard]] MhdFaceFluxMutableView makeMutableMagneticFaceView(
    core::MagneticFaceSoa& storage,
    const MhdAdapterPatchBinding& binding,
    std::size_t patch_index);

// Convert between canonical core gas rows and the MHD conserved state.
// `dense_row_by_geometry_row` maps geometry-order cells to dense CellSoa rows.
// Generalized to wrapping index type is omitted deliberately: patches are
// contiguous geometry-order slices of the fixed layout.
void loadMhdGasFromState(const core::SimulationState& state,
                         const MhdPatchGeometry& geometry,
                         std::size_t first_geometry_row,
                         std::span<const std::uint32_t> dense_row_by_geometry_row,
                         MhdGasStateSoa& gas);

void storeMhdGasToState(core::SimulationState& state,
                        const MhdPatchGeometry& geometry,
                        std::size_t first_geometry_row,
                        std::span<const std::uint32_t> geometry_row_by_dense_row,
                        const MhdGasStateSoa& gas,
                        double density_floor_code,
                        double pressure_floor_code);

// Initialise an all-zero persistent flux patch from a uniform comoving field
// (Bx, By, Bz). Used only when a run has no explicit face-flux input. Returns
// false when the patch is already non-empty or the geometry is invalid.
[[nodiscard]] bool seedUniformComovingField(core::MagneticFaceSoa& storage,
                                            std::size_t patch_index,
                                            double bx_comoving,
                                            double by_comoving,
                                            double bz_comoving);

// Sum the magnetic energy of one patch (comoving code units), for diagnostics.
[[nodiscard]] double patchMagneticEnergy(const core::MagneticFaceSoa& storage,
                                         std::size_t patch_index);

// Persistent backing for a configured coefficient provider. Table and chemistry
// providers fail closed until their audited data contracts exist; only the
// constant provider is currently instantiable.
class MhdConfiguredCoefficientProvider final : public MhdCoefficientProvider {
 public:
  explicit MhdConfiguredCoefficientProvider(const core::MhdConfig& config);
  ~MhdConfiguredCoefficientProvider() override = default;

  [[nodiscard]] MhdCoefficientProviderKind kind() const noexcept override {
    return MhdCoefficientProviderKind::kConstant;
  }
  [[nodiscard]] std::string provenanceId() const override;
  [[nodiscard]] std::optional<std::string> tableSha256() const override {
    return std::nullopt;
  }
  void validateAdmission() const override;
  [[nodiscard]] MhdTransportCoefficients evaluate(
      std::size_t cell_index,
      const MhdCellPhysicalProperties& properties) const override;

 private:
  ConstantMhdCoefficientProvider m_constant;
};

}  // namespace cosmosim::mhd
