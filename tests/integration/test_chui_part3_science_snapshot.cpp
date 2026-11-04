// CHUI Part 3 Checkpoint G — composition and stellar-population science
// snapshot export (CHUI_CompositionSchemaVersion = 1,
// CHUI_StellarPopulationSchemaVersion = 1).
// STATUS: written as source only. NOT RUN — USER-RESERVED. Skips (77) without
// an HDF5 build.
//
// Invariants:
//  * a gas snapshot of a run with a composition sidecar carries
//    /PartType0/CHUI_ElementMassFractions [N, 11] (+ species, dust, SGS, mu,
//    tags, dwell) whose row sums (elements + dust) equal 1 to 1e-12 and whose
//    values equal lane mass / composition gas mass of the same stable
//    GasCellIDs row;
//  * canonical GADGET/AREPO datasets are unchanged (Masses, Density,
//    InternalEnergy, Metallicity present);
//  * a run without composition writes NO composition datasets (absent means
//    disabled; never zero-filled).
#include "../support/test_require.hpp"

#include <cmath>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "../support/test_temp_workspace.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace {

#if COSMOSIM_ENABLE_HDF5
namespace core = cosmosim::core;
namespace ph = cosmosim::physics;

core::SimulationState makeGasState(bool composition) {
  core::SimulationState state;
  state.resizeParticles(0U);
  state.resizeCells(3U);
  for (std::size_t row = 0; row < 3U; ++row) {
    state.cells.center_x_comoving[row] = 0.1 + 0.2 * static_cast<double>(row);
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.mass_code[row] = 2.0 + static_cast<double>(row);
    state.cells.patch_index[row] = 0U;
    state.gas_cells.density_code[row] = 4.0 + static_cast<double>(row);
    state.gas_cells.internal_energy_code[row] = 5.0 + static_cast<double>(row);
    state.gas_cells.metal_mass_code[row] = (0.01 + 0.01 * static_cast<double>(row)) * state.cells.mass_code[row];
  }
  state.replaceGasCellIdentityRecords({
      {.gas_cell_id = 7001U, .parent_particle_id = std::nullopt, .owning_patch_id = 42U, .local_cell_row = 0U},
      {.gas_cell_id = 7002U, .parent_particle_id = std::nullopt, .owning_patch_id = 42U, .local_cell_row = 1U},
      {.gas_cell_id = 7003U, .parent_particle_id = std::nullopt, .owning_patch_id = 42U, .local_cell_row = 2U}});
  state.rebuildSpeciesIndex();
  if (composition) {
    ph::enableGasComposition(state, ph::CompositionLayout::make({
        .network = ph::CompositionNetworkKind::kHybridLight, .dust = ph::CompositionDustKind::kTwoReservoirV1,
        .sgs_energy = true}), {});
  }
  return state;
}

std::filesystem::path writeSnapshot(const core::SimulationState& state, const char* name) {
  static auto config = core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = name;
  const std::filesystem::path path =
      cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath(std::string(name) + ".hdf5");
  cosmosim::io::SnapshotWritePayload payload;
  payload.state = &state;
  payload.config = &config;
  payload.normalized_config_text = "schema_version = 1\n";
  payload.provenance = core::makeProvenanceRecord(name, name, 0, payload.normalized_config_text);
  cosmosim::io::writeGadgetArepoSnapshotHdf5(path, payload);
  return path;
}

bool exists(hid_t file, const char* path) { return H5Lexists(file, path, H5P_DEFAULT) > 0; }

std::vector<double> read2d(hid_t file, const char* path, hsize_t& rows, hsize_t& cols) {
  const hid_t ds = H5Dopen2(file, path, H5P_DEFAULT);
  CHUI_REQUIRE(ds >= 0);
  const hid_t space = H5Dget_space(ds);
  hsize_t dims[2] = {0, 0};
  CHUI_REQUIRE(H5Sget_simple_extent_ndims(space) == 2);
  H5Sget_simple_extent_dims(space, dims, nullptr);
  rows = dims[0];
  cols = dims[1];
  std::vector<double> values(static_cast<std::size_t>(rows * cols));
  CHUI_REQUIRE(H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0);
  H5Sclose(space);
  H5Dclose(ds);
  return values;
}

void testCompositionExported() {
  const core::SimulationState state = makeGasState(true);
  const auto path = writeSnapshot(state, "chui_part3_composition_snapshot");
  const hid_t file = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  CHUI_REQUIRE(file >= 0);
  for (const char* canonical : {"/PartType0/Masses", "/PartType0/Density", "/PartType0/InternalEnergy",
                                "/PartType0/Metallicity", "/PartType0/GasCellIDs"}) {
    CHUI_REQUIRE(exists(file, canonical));
  }
  hsize_t rows = 0, cols = 0, dust_rows = 0, dust_cols = 0;
  const auto elements = read2d(file, "/PartType0/CHUI_ElementMassFractions", rows, cols);
  const auto dust = read2d(file, "/PartType0/CHUI_DustMassFractions", dust_rows, dust_cols);
  CHUI_REQUIRE(rows == 3U && cols == 11U && dust_rows == 3U);
  for (const char* optional : {"/PartType0/CHUI_SpeciesMassFractions", "/PartType0/CHUI_SgsSpecificEnergy",
                               "/PartType0/CHUI_MeanMolecularWeight", "/PartType0/CHUI_ChemistryModelTag",
                               "/PartType0/CHUI_IsmRegimeTag"}) {
    CHUI_REQUIRE(exists(file, optional));
  }
  const auto& sidecar = state.gas_composition;
  for (hsize_t r = 0; r < rows; ++r) {
    double sum = 0.0;
    for (hsize_t c = 0; c < cols; ++c) sum += elements[r * cols + c];
    for (hsize_t c = 0; c < dust_cols; ++c) sum += dust[r * dust_cols + c];
    CHUI_REQUIRE(std::abs(sum - 1.0) < 1.0e-12);
    const double gas = ph::compositionGasMass(sidecar, r);
    CHUI_REQUIRE(std::abs(elements[r * cols + 4U] -
                          sidecar.value(sidecar.descriptor().element_lane_begin + 4U, r) / gas) < 1.0e-14);
  }
  H5Fclose(file);
}

void testDisabledCompositionWritesNothing() {
  const core::SimulationState state = makeGasState(false);
  const auto path = writeSnapshot(state, "chui_part3_no_composition_snapshot");
  const hid_t file = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  CHUI_REQUIRE(file >= 0);
  CHUI_REQUIRE(!exists(file, "/PartType0/CHUI_ElementMassFractions"));
  CHUI_REQUIRE(!exists(file, "/PartType0/CHUI_DustMassFractions"));
  CHUI_REQUIRE(exists(file, "/PartType0/Metallicity"));
  H5Fclose(file);
}
#endif

}  // namespace

int main() {
#if COSMOSIM_ENABLE_HDF5
  testCompositionExported();
  testDisabledCompositionWritesNothing();
  return 0;
#else
  return 77;  // legitimate skip: HDF5 not built
#endif
}
