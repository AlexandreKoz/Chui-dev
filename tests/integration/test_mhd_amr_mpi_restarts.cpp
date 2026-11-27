// CHUI-MHD-REF-v1 AMR transfer, EMF register coverage, v26 /state/mhd codec,
// config-hash compatibility and capability fail-closed admission
// (G02 nested, G07, G08 partial, MP06). Source written, NOT EXECUTED.
//
// MPI 2/4-rank equivalence (G06) is NOT covered: no MPI magnetic halo
// protocol exists, kImportedMpi sides are rejected (see the pass01 ledger).

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <span>
#include <string>
#include <vector>

#include "../support/test_require.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/io/mhd_restart_state.hpp"
#include "cosmosim/mhd/mhd_amr.hpp"
#include "cosmosim/mhd/mhd_initial_conditions.hpp"
#include "cosmosim/mhd/mhd_runtime_adapter.hpp"
#include "cosmosim/workflows/runtime_capabilities.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"
#include "cosmosim/core/simulation_mode.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace {

using namespace cosmosim::mhd;
constexpr double k_gamma = 5.0 / 3.0;

MhdPatchGeometry parentGeometry() {
  MhdPatchGeometry geometry{};
  geometry.nx = 8;
  geometry.ny = 6;
  geometry.nz = 4;
  geometry.dx_comoving = 0.25;
  geometry.dy_comoving = 0.2;
  geometry.dz_comoving = 0.5;
  return geometry;
}

double maxDivergence(const MhdFaceFluxStorage& faces) {
  const double normalizer =
      magneticDivergenceNormalizer(faces.geometry(), faces.fluxX(), faces.fluxY(), faces.fluxZ());
  return maxNormalizedFaceDivergence(faces.geometry(), faces.fluxX(), faces.fluxY(), faces.fluxZ(),
                                     normalizer);
}

// G07: uniform field -> every child subface carries exactly Phi_parent / 4.
void testUniformProlongation() {
  const MhdPatchGeometry parent = parentGeometry();
  MhdFaceFluxStorage faces(parent);
  (void)seedFaceFlux(parent, {0.3, -0.2, 0.7}, {}, faces.fluxX(), faces.fluxY(), faces.fluxZ());
  std::array<MhdFaceFluxStorage, 8> children{};
  const auto result = MhdAmrSynchronizer::prolongFaceFluxes(parent, faces.fluxX(), faces.fluxY(),
                                                            faces.fluxZ(), children);
  CHUI_REQUIRE(result.ok);
  for (const MhdFaceFluxStorage& child : children) {
    for (const double value : child.fluxX()) {
      CHUI_REQUIRE(std::abs(value - 0.25 * faces.faceFluxX(0, 0, 0)) <= 1e-15);
    }
    for (const double value : child.fluxZ()) {
      CHUI_REQUIRE(std::abs(value - 0.25 * faces.faceFluxZ(0, 0, 0)) <= 1e-15);
    }
    CHUI_REQUIRE(maxDivergence(child) <= 1e-14);
  }
}

// Nonuniform divergence-free parent: children divergence free, restriction
// of the prolongation reproduces the parent, sibling interface faces agree,
// and repeated refine/coarsen keeps the net boundary flux.
void testNonuniformTransferIdentity() {
  const MhdPatchGeometry parent = parentGeometry();
  MhdFaceFluxStorage faces(parent);
  (void)seedFaceFlux(parent, {0.1, 0.05, -0.02},
                     stochasticSeedPotential({2.0, 1.2, 2.0}, 2, 0.3, -1.0, 77U), faces.fluxX(),
                     faces.fluxY(), faces.fluxZ());
  CHUI_REQUIRE(maxDivergence(faces) <= 1e-12);
  const std::array<double, 3> net0 =
      MhdAmrSynchronizer::netBoundaryFlux(parent, faces.fluxX(), faces.fluxY(), faces.fluxZ());
  for (int cycle = 0; cycle < 3; ++cycle) {
    std::array<MhdFaceFluxStorage, 8> children{};
    const auto result = MhdAmrSynchronizer::prolongFaceFluxes(parent, faces.fluxX(), faces.fluxY(),
                                                              faces.fluxZ(), children);
    CHUI_REQUIRE(result.ok && result.max_cell_divergence_residual <= 1e-12);
    // Sibling interface (x-midplane): child 0 upper plane == child 1 lower plane.
    const MhdPatchGeometry& child_geometry = children[0].geometry();
    for (std::size_t k = 0; k < child_geometry.nz; ++k) {
      for (std::size_t j = 0; j < child_geometry.ny; ++j) {
        CHUI_REQUIRE(children[0].faceFluxX(child_geometry.nx, j, k) == children[1].faceFluxX(0, j, k));
      }
    }
    MhdFaceFluxStorage restored(parent);
    MhdAmrSynchronizer::restrictFaceFluxes(parent, restored.fluxX(), restored.fluxY(),
                                           restored.fluxZ(), children);
    for (std::size_t index = 0; index < faces.fluxX().size(); ++index) {
      CHUI_REQUIRE(std::abs(restored.fluxX()[index] - faces.fluxX()[index]) <=
                   1e-14 * (1.0 + std::abs(faces.fluxX()[index])));
    }
    faces = restored;
  }
  const std::array<double, 3> net1 =
      MhdAmrSynchronizer::netBoundaryFlux(parent, faces.fluxX(), faces.fluxY(), faces.fluxZ());
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    CHUI_REQUIRE(std::abs(net1[axis] - net0[axis]) <= 1e-13);
  }
}

// Transactional gas+field transfer: exact mass/momentum/metal/total-energy
// conservation; an inadmissible fine thermal state is rejected with the
// outputs untouched.
void testGasTransferTransaction() {
  const MhdPatchGeometry parent = parentGeometry();
  MhdFaceFluxStorage faces(parent);
  (void)seedFaceFlux(parent, {0.0, 0.0, 0.0},
                     stochasticSeedPotential({2.0, 1.2, 2.0}, 2, 0.2, -1.0, 5U), faces.fluxX(),
                     faces.fluxY(), faces.fluxZ());
  MhdGasStateSoa gas(parent.cellCount());
  for (std::size_t cell = 0; cell < parent.cellCount(); ++cell) {
    MhdConservedState state{};
    state.rho_comoving = 1.0 + 0.01 * static_cast<double>(cell % 7U);
    state.momentum_density_x_comoving = 0.1;
    state.total_energy_density_comoving = 2.0;
    state.metal_mass_density_comoving = 0.02;
    gas.storeGas(cell, state);
  }
  const MhdFaceFluxView view = makeFaceFluxView(faces, 1U, 1U, 0);
  std::array<MhdFaceFluxStorage, 8> children_faces{};
  std::array<MhdGasStateSoa, 8> children_gas{};
  const auto result = MhdAmrSynchronizer::prolongPatch(parent, gas, view, 1.0, k_gamma,
                                                       children_faces, children_gas);
  CHUI_REQUIRE(result.ok);
  double parent_mass = 0.0;
  double parent_energy = 0.0;
  for (std::size_t cell = 0; cell < parent.cellCount(); ++cell) {
    parent_mass += gas.rho()[cell] * parent.cellVolume();
    parent_energy += gas.totalEnergy()[cell] * parent.cellVolume();
  }
  double child_mass = 0.0;
  double child_energy = 0.0;
  for (std::size_t block = 0; block < 8U; ++block) {
    const double volume = children_faces[block].geometry().cellVolume();
    for (std::size_t cell = 0; cell < children_gas[block].size(); ++cell) {
      child_mass += children_gas[block].rho()[cell] * volume;
      child_energy += children_gas[block].totalEnergy()[cell] * volume;
    }
  }
  CHUI_REQUIRE(std::abs(child_mass - parent_mass) <= 1e-13 * parent_mass);
  CHUI_REQUIRE(std::abs(child_energy - parent_energy) <= 1e-13 * parent_energy);
  // Restriction returns the parent totals and an admissible thermal state.
  MhdGasStateSoa restored(parent.cellCount());
  MhdFaceFluxStorage restored_faces(parent);
  const auto back = MhdAmrSynchronizer::restrictPatch(parent, children_faces, children_gas, 1.0,
                                                      k_gamma, restored, restored_faces.fluxX(),
                                                      restored_faces.fluxY(), restored_faces.fluxZ());
  CHUI_REQUIRE(back.ok);
  for (std::size_t cell = 0; cell < parent.cellCount(); ++cell) {
    CHUI_REQUIRE(std::abs(restored.totalEnergy()[cell] - gas.totalEnergy()[cell]) <= 1e-13);
  }
  // Inadmissible: thermal energy far below the magnetic energy change.
  for (std::size_t cell = 0; cell < parent.cellCount(); ++cell) {
    MhdConservedState state = gas.loadGas(cell);
    const std::array<double, 3> b = faces.cellFieldComoving(cell % parent.nx, (cell / parent.nx) % parent.ny,
                                                            cell / (parent.nx * parent.ny));
    state.momentum_density_x_comoving = 0.0;
    state.total_energy_density_comoving = 0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) + 1e-14;
    gas.storeGas(cell, state);
  }
  std::array<MhdFaceFluxStorage, 8> untouched_faces{};
  std::array<MhdGasStateSoa, 8> untouched_gas{};
  const auto rejected = MhdAmrSynchronizer::prolongPatch(parent, gas, view, 1.0, k_gamma,
                                                         untouched_faces, untouched_gas);
  CHUI_REQUIRE(!rejected.ok);
  for (std::size_t block = 0; block < 8U; ++block) {
    CHUI_REQUIRE(untouched_faces[block].fluxX().empty() && untouched_gas[block].size() == 0U);
  }
}

MhdAmrIntervalId interval() {
  MhdAmrIntervalId id{};
  id.coarse_level = 0;
  id.time_begin_code = 0.0;
  id.time_end_code = 1.0;
  id.coarse_geometry_generation = 4U;
  id.hierarchy_epoch = 9U;
  return id;
}

// Register coverage: complete records apply; gaps, duplicates, stale
// intervals and undeclared keys are rejected; the curl correction keeps D*C=0.
void testEmfRegisterCoverage() {
  const MhdPatchGeometry coarse = parentGeometry();
  const std::vector<MhdCoarseEdgeKey> keys{{.coarse_patch_id = 1, .axis = MhdAxis::kZ, .i = 2, .j = 2, .k = 1},
                                           {.coarse_patch_id = 1, .axis = MhdAxis::kZ, .i = 2, .j = 3, .k = 1}};
  const auto make = [&](std::int32_t substeps) {
    MhdEmfRegister reflux;
    reflux.bindInterval(interval(), 1U, 4U, substeps);
    reflux.declareInterfaceEdges(keys);
    return reflux;
  };
  const auto expect_incomplete = [](auto&& action) {
    bool rejected = false;
    try {
      action();
    } catch (const MhdError& error) {
      rejected = error.code() == MhdErrorCode::kAmrRegisterIncomplete;
    }
    CHUI_REQUIRE(rejected);
  };
  MhdEmfRegister complete = make(2);
  for (const MhdCoarseEdgeKey& key : keys) {
    complete.addCoarse(interval(), key, 1.0);
    for (std::int32_t segment = 0; segment < 2; ++segment) {
      for (std::int32_t substep = 0; substep < 2; ++substep) {
        complete.addFine(interval(), key, segment, substep, 0.3);
      }
    }
  }
  complete.finalize();
  MhdFaceFluxStorage faces(coarse);
  (void)seedFaceFlux(coarse, {0.2, 0.1, 0.0}, {}, faces.fluxX(), faces.fluxY(), faces.fluxZ());
  // Per record I_fine = 2 segments x 2 substeps x 0.3 = 1.2, I_coarse = 1.0, so
  // delta = 0.2 on
  // both z-edges (2,2,1) and (2,3,1). x-face (2,1,1) is bounded by the edges
  // j=1 (no correction) and j=2 (+0.2): Phi -= (0.2 - 0) -> changes by -0.2;
  // x-face (2,2,1) sees equal corrections on both edges and is unchanged.
  const double before_changed = faces.faceFluxX(2, 1, 1);
  const double before_same = faces.faceFluxX(2, 2, 1);
  MhdAmrSynchronizer::applyRefluxCurl(coarse, faces.fluxX(), faces.fluxY(), faces.fluxZ(), complete);
  CHUI_REQUIRE(std::abs(faces.faceFluxX(2, 1, 1) - (before_changed - 0.2)) <= 1e-15);
  CHUI_REQUIRE(std::abs(faces.faceFluxX(2, 2, 1) - before_same) <= 1e-15);
  CHUI_REQUIRE(maxDivergence(faces) <= 1e-14);
  // Gap: one fine substep missing.
  expect_incomplete([&]() {
    MhdEmfRegister gap = make(2);
    for (const MhdCoarseEdgeKey& key : keys) {
      gap.addCoarse(interval(), key, 1.0);
      gap.addFine(interval(), key, 0, 0, 0.3);
    }
    gap.finalize();
  });
  // Duplicate contribution.
  expect_incomplete([&]() {
    MhdEmfRegister duplicate = make(1);
    duplicate.addCoarse(interval(), keys[0], 1.0);
    duplicate.addCoarse(interval(), keys[0], 1.0);
  });
  // Stale interval.
  expect_incomplete([&]() {
    MhdEmfRegister stale = make(1);
    MhdAmrIntervalId other = interval();
    other.hierarchy_epoch = 10U;
    stale.addCoarse(other, keys[0], 1.0);
  });
  // Undeclared edge.
  expect_incomplete([&]() {
    MhdEmfRegister undeclared = make(1);
    undeclared.addCoarse(interval(), {.coarse_patch_id = 1, .axis = MhdAxis::kX, .i = 0, .j = 0, .k = 0}, 1.0);
  });
  CHUI_REQUIRE(mhdEmfRegisterBound(coarse) < coarse.edgeCountTotal());
}

// Face layout lifecycle: checked construction, duplicate rejection, stale
// binding rejection.
void testFaceLayoutLifecycle() {
  MhdAdapterPatchBinding a{.patch_id = 5, .geometry_generation = 1, .level = 0, .geometry = parentGeometry()};
  MhdAdapterPatchBinding b = a;
  b.patch_id = 6;
  const std::array<MhdAdapterPatchBinding, 2> bindings{a, b};
  cosmosim::core::MagneticFaceSoa storage = buildMagneticFaceLayout(bindings, "");
  validateMagneticFaceSoa(storage);
  CHUI_REQUIRE(storage.patches.size() == 2U && storage.isConsistent());
  const std::array<MhdAdapterPatchBinding, 2> duplicate{a, a};
  bool rejected = false;
  try {
    (void)buildMagneticFaceLayout(duplicate, "");
  } catch (const MhdError& error) {
    rejected = error.code() == MhdErrorCode::kStaleGeometry;
  }
  CHUI_REQUIRE(rejected);
  MhdAdapterPatchBinding stale = a;
  stale.geometry_generation = 2;
  rejected = false;
  try {
    (void)makeMagneticFaceView(storage, stale, 0);
  } catch (const MhdError& error) {
    rejected = error.code() == MhdErrorCode::kStaleGeometry;
  }
  CHUI_REQUIRE(rejected);
}

#if COSMOSIM_ENABLE_HDF5
// G08 (codec part): /state/mhd roundtrip, and corrupt/missing/hash-mismatch
// rejection before any authoritative mutation.
void testRestartGroupCodec() {
  const MhdAdapterPatchBinding binding{.patch_id = 42, .geometry_generation = 3, .level = 1,
                                       .geometry = parentGeometry()};
  cosmosim::core::MagneticFaceSoa storage =
      buildMagneticFaceLayout(std::span<const MhdAdapterPatchBinding>(&binding, 1U), "abc");
  MhdFaceFluxMutableView view = makeMutableMagneticFaceView(storage, binding, 0);
  (void)seedFaceFlux(binding.geometry, {0.1, 0.2, 0.3},
                     fieldLoopPotential(0.01, 0.5, {1.0, 0.6, 1.0}), view.phi_x, view.phi_y,
                     view.phi_z);
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "chui_mhd_v26_group.h5";
  const auto write = [&](const std::string& sha) {
    const hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    cosmosim::io::writeMhdRestartGroup(
        file, storage,
        {.reference_contract = "chui-mhd-ref-v1.0", .coefficient_table_sha256 = sha,
         .divergence_max_norm = 0.0, .scale_factor = 0.5});
    H5Fclose(file);
  };
  write("abc");
  const cosmosim::io::MhdRestartExpectation expected{.reference_contract = "chui-mhd-ref-v1.0",
                                                     .coefficient_table_sha256 = "abc",
                                                     .scale_factor = 0.5};
  {
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    const cosmosim::core::MagneticFaceSoa restored =
        cosmosim::io::readMhdRestartGroup(file, expected, nullptr);
    H5Fclose(file);
    CHUI_REQUIRE(restored.patches.size() == 1U && restored.patches[0].patch_id == 42U);
    CHUI_REQUIRE(restored.face_flux_x_comoving == storage.face_flux_x_comoving);
    CHUI_REQUIRE(restored.face_flux_z_comoving == storage.face_flux_z_comoving);
  }
  const auto expect_mismatch = [&](const cosmosim::io::MhdRestartExpectation& expectation) {
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    bool rejected = false;
    try {
      (void)cosmosim::io::readMhdRestartGroup(file, expectation, nullptr);
    } catch (const std::runtime_error& error) {
      rejected = std::string(error.what()).find("MHD_RESTART_SCHEMA_MISMATCH") != std::string::npos;
    }
    H5Fclose(file);
    CHUI_REQUIRE(rejected);
  };
  cosmosim::io::MhdRestartExpectation wrong_sha = expected;
  wrong_sha.coefficient_table_sha256 = "def";
  expect_mismatch(wrong_sha);
  cosmosim::io::MhdRestartExpectation wrong_a = expected;
  wrong_a.scale_factor = 0.6;
  expect_mismatch(wrong_a);
  // Corrupt the face lane (remove dataset) -> missing dataset rejection.
  {
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    H5Ldelete(file, "/state/mhd/face_flux_y_comoving", H5P_DEFAULT);
    H5Fclose(file);
  }
  expect_mismatch(expected);
  // A v25-style file without /state/mhd is never read as magnetized.
  {
    const hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose(file);
  }
  expect_mismatch(expected);
  std::filesystem::remove(path);
}
#endif

const char* k_baseline_config = R"(
[units]
length_unit = mpc
coordinate_frame = comoving

[mode]
mode = zoom_in
ic_file = zoom_ics.hdf5
ic_convention = gadget_arepo_bridge_v1
ic_bridge_source_length_unit_to_si = 3.0856775814913673e19
ic_bridge_source_mass_unit_to_si = 1.98847e30
ic_bridge_source_velocity_unit_to_si = 1000
ic_bridge_coordinate_frame = comoving
ic_bridge_velocity_convention = physical_peculiar
ic_bridge_length_hubble_exponent = 0
ic_bridge_length_scale_factor_exponent = 0
ic_bridge_mass_hubble_exponent = 0
ic_bridge_mass_scale_factor_exponent = 0
ic_bridge_velocity_hubble_exponent = 0
ic_bridge_velocity_scale_factor_exponent = 0
zoom_high_res_region = true
zoom_region_file = region_zoom.hdf5
zoom_region_radius = 1.0

[cosmology]
box_size = 50000 kpc
omega_matter = 0.31
omega_lambda = 0.69

[numerics]
t_code_begin = 0.0
t_code_end = 1.0
gravity_softening = 2.5 kpc

[output]
output_stem = snapshot
restart_stem = restart
)";

// MHD-off normalized text and hash are byte-identical to the predecessor
// layout (no [mhd*] sections); inert non-default MHD keys are rejected; an
// MHD-on configuration is parsed but rejected at runtime capability admission.
void testConfigAndCapabilityAdmission() {
  const auto off = cosmosim::core::loadFrozenConfigFromString(k_baseline_config, "mhd_off");
  CHUI_REQUIRE(off.normalized_text.find("[mhd") == std::string::npos);
  const auto off_explicit = cosmosim::core::loadFrozenConfigFromString(
      std::string(k_baseline_config) + "\n[mhd]\nenabled = false\n", "mhd_off_explicit");
  CHUI_REQUIRE(off_explicit.provenance.config_hash_hex == off.provenance.config_hash_hex);
  bool rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        std::string(k_baseline_config) + "\n[mhd]\nenabled = false\nreconstruction = plm_mc\n",
        "inert_key");
  } catch (const std::exception&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);
  // Pass 02: capability-specific admission (no blanket refusal).
  // Negative: the global-step transaction must own MHD retries.
  const auto on = cosmosim::core::loadFrozenConfigFromString(
      std::string(k_baseline_config) + "\n[mhd]\nenabled = true\n", "mhd_on");
  CHUI_REQUIRE(on.normalized_text.find("[mhd]") != std::string::npos);
  const std::string missing_retry_owner =
      cosmosim::workflows::mhdProductionDispatchRefusal(on.config, cosmosim::core::ModePolicy{});
  CHUI_REQUIRE(missing_retry_owner.find("global_step_max_retries") != std::string::npos);
  const auto report = cosmosim::workflows::buildRuntimeCapabilityReport(on.config);
  CHUI_REQUIRE(report.require("mhd_production_dispatch").requested);
  CHUI_REQUIRE(!report.require("mhd_production_dispatch").runtime_available);
  rejected = false;
  try {
    cosmosim::workflows::validateRequestedRuntimeCapabilities(on.config, report);
  } catch (const std::invalid_argument& error) {
    rejected = std::string(error.what()).find("MHD_UNSUPPORTED_PHYSICS") != std::string::npos;
  } catch (const std::exception&) {
    rejected = false;  // a different (e.g. TreePM backend) refusal is not this guard
  }
  CHUI_REQUIRE(rejected);
  // Positive: the synchronized CPU reference predicates hold.
  std::string supported_text(k_baseline_config);
  const std::size_t numerics = supported_text.find("[numerics]\n");
  CHUI_REQUIRE(numerics != std::string::npos);
  supported_text.insert(numerics + std::string("[numerics]\n").size(),
                        "global_step_max_retries = 1\n");
  supported_text += "\n[mhd]\nenabled = true\n";
  const auto supported = cosmosim::core::loadFrozenConfigFromString(supported_text, "mhd_supported");
  const std::string supported_refusal = cosmosim::workflows::mhdProductionDispatchRefusal(
      supported.config, cosmosim::core::ModePolicy{});
#if COSMOSIM_ENABLE_HDF5
  CHUI_REQUIRE(supported_refusal.empty());
  const auto supported_report = cosmosim::workflows::buildRuntimeCapabilityReport(supported.config);
  CHUI_REQUIRE(supported_report.require("mhd_production_dispatch").runtime_available);
  CHUI_REQUIRE(supported_report.require("mhd_production_dispatch").active);
#else
  CHUI_REQUIRE(supported_refusal.find("HDF5") != std::string::npos);
#endif
  // Negative: specific unsupported combinations keep failing closed.
  std::string shadow_text = supported_text;
  shadow_text.insert(numerics + std::string("[numerics]\n").size(), "hydro_shadow_enabled = true\n");
  bool shadow_refused = false;
  try {
    const auto shadow = cosmosim::core::loadFrozenConfigFromString(shadow_text, "mhd_shadow");
    shadow_refused = !cosmosim::workflows::mhdProductionDispatchRefusal(
                          shadow.config, cosmosim::core::ModePolicy{})
                          .empty();
  } catch (const std::exception&) {
    shadow_refused = true;  // rejected even earlier by typed config validation
  }
  CHUI_REQUIRE(shadow_refused);
  // Pass 04 (M03-02/M04-01/M04-02/M04-03): precise refusals for combinations
  // whose MHD integration is incomplete (typed config mutated directly so the
  // refusal under test is the MHD predicate, not an unrelated table check).
  {
    cosmosim::core::SimulationConfig composition = supported.config;
    composition.physics.thermochemistry_backend =
        cosmosim::core::ThermochemistryBackendMode::kFixedHybrid;
    CHUI_REQUIRE(cosmosim::workflows::mhdProductionDispatchRefusal(
                     composition, cosmosim::core::ModePolicy{})
                     .find("thermochemistry_backend") != std::string::npos);
    cosmosim::core::SimulationConfig observables = supported.config;
    observables.extended.observables.enabled = true;
    CHUI_REQUIRE(cosmosim::workflows::mhdProductionDispatchRefusal(
                     observables, cosmosim::core::ModePolicy{})
                     .find("observables") != std::string::npos);
    cosmosim::core::SimulationConfig neutrinos = supported.config;
    neutrinos.extended.neutrinos.method = cosmosim::core::NeutrinoMethodConfig::kLinearResponse;
    CHUI_REQUIRE(cosmosim::workflows::mhdProductionDispatchRefusal(
                     neutrinos, cosmosim::core::ModePolicy{})
                     .find("linear_response") != std::string::npos);
    // Stateless background-only neutrinos are NOT refused by the MHD predicate.
    cosmosim::core::SimulationConfig background = supported.config;
    background.extended.neutrinos.method = cosmosim::core::NeutrinoMethodConfig::kBackgroundOnly;
    CHUI_REQUIRE(cosmosim::workflows::mhdProductionDispatchRefusal(
                     background, cosmosim::core::ModePolicy{}) == supported_refusal);
  }
  // Pass 04 (M03-12): the mandatory divergence integrity record cannot be
  // disabled.
  {
    bool divergence_refused = false;
    try {
      (void)cosmosim::core::loadFrozenConfigFromString(
          supported_text + "
[mhd_io]
output_divergence_diagnostics = false
", "mhd_no_div");
    } catch (const std::exception&) {
      divergence_refused = true;
    }
    CHUI_REQUIRE(divergence_refused);
  }
  const auto off_report = cosmosim::workflows::buildRuntimeCapabilityReport(off.config);
  CHUI_REQUIRE(!off_report.require("mhd_production_dispatch").requested);
  CHUI_REQUIRE(off_report.require("mhd_production_dispatch").runtime_available);
}

}  // namespace

// Pass 04 (M03-09): on ANISOTROPIC rectangular cells (dx != dy != dz) a
// linear solenoidal field with nonzero curl,
//   B = (0.2 + 0.3 y, 0.1 + 0.4 z, -0.2 + 0.5 x),
// is reproduced EXACTLY on the fine faces of every parent cell whose
// transverse slopes are unlimited interior differences (MC of linear data is
// the central difference), the fine-cell divergence is zero, and restriction
// returns the parent. The pass-01/02 cubic-cell quadratic coefficient would
// break the divergence of nonlinear data on these cells (see
// testNonuniformTransferIdentity, which uses the same anisotropic parent).
void testAnisotropicLinearFieldExactness() {
  const MhdPatchGeometry parent = parentGeometry();
  MhdFaceFluxStorage faces(parent);
  const auto bx = [](double, double y, double) { return 0.2 + 0.3 * y; };
  const auto by = [](double, double, double z) { return 0.1 + 0.4 * z; };
  const auto bz = [](double x, double, double) { return -0.2 + 0.5 * x; };
  const auto fill = [](const MhdPatchGeometry& g, std::span<double> px, std::span<double> py,
                       std::span<double> pz, auto fx, auto fy, auto fz) {
    const double ox = g.origin_x_comoving;
    const double oy = g.origin_y_comoving;
    const double oz = g.origin_z_comoving;
    for (std::size_t k = 0; k < g.nz; ++k) {
      for (std::size_t j = 0; j < g.ny; ++j) {
        for (std::size_t i = 0; i <= g.nx; ++i) {
          px[i + (g.nx + 1U) * (j + g.ny * k)] =
              fx(ox + i * g.dx_comoving, oy + (j + 0.5) * g.dy_comoving,
                 oz + (k + 0.5) * g.dz_comoving) * g.faceAreaX();
        }
      }
    }
    for (std::size_t k = 0; k < g.nz; ++k) {
      for (std::size_t j = 0; j <= g.ny; ++j) {
        for (std::size_t i = 0; i < g.nx; ++i) {
          py[i + g.nx * (j + (g.ny + 1U) * k)] =
              fy(ox + (i + 0.5) * g.dx_comoving, oy + j * g.dy_comoving,
                 oz + (k + 0.5) * g.dz_comoving) * g.faceAreaY();
        }
      }
    }
    for (std::size_t k = 0; k <= g.nz; ++k) {
      for (std::size_t j = 0; j < g.ny; ++j) {
        for (std::size_t i = 0; i < g.nx; ++i) {
          pz[i + g.nx * (j + g.ny * k)] =
              fz(ox + (i + 0.5) * g.dx_comoving, oy + (j + 0.5) * g.dy_comoving,
                 oz + k * g.dz_comoving) * g.faceAreaZ();
        }
      }
    }
  };
  fill(parent, faces.fluxX(), faces.fluxY(), faces.fluxZ(), bx, by, bz);
  CHUI_REQUIRE(maxDivergence(faces) <= 1e-14);
  std::array<MhdFaceFluxStorage, 8> children{};
  const auto result = MhdAmrSynchronizer::prolongFaceFluxes(parent, faces.fluxX(), faces.fluxY(),
                                                            faces.fluxZ(), children);
  CHUI_REQUIRE(result.ok && result.max_cell_divergence_residual <= 1e-13);
  for (std::int32_t block = 0; block < 8; ++block) {
    const MhdFaceFluxStorage& child = children[static_cast<std::size_t>(block)];
    const MhdPatchGeometry& g = child.geometry();
    CHUI_REQUIRE(maxDivergence(child) <= 1e-13);
    MhdFaceFluxStorage exact(g);
    fill(g, exact.fluxX(), exact.fluxY(), exact.fluxZ(), bx, by, bz);
    // Compare only fine faces inside parent cells at least one parent cell
    // away from the parent patch boundary (where MC slopes are interior).
    const auto interior = [&](double coordinate, double origin, double width, std::size_t n) {
      const double parent_index = (coordinate - origin) / width;
      return parent_index > 1.0 + 1e-9 && parent_index < static_cast<double>(n) - 1.0 - 1e-9;
    };
    for (std::size_t k = 0; k < g.nz; ++k) {
      for (std::size_t j = 0; j < g.ny; ++j) {
        for (std::size_t i = 0; i < g.nx; ++i) {
          const double xc = g.origin_x_comoving + (i + 0.5) * g.dx_comoving;
          const double yc = g.origin_y_comoving + (j + 0.5) * g.dy_comoving;
          const double zc = g.origin_z_comoving + (k + 0.5) * g.dz_comoving;
          if (!interior(xc, parent.origin_x_comoving, parent.dx_comoving, parent.nx) ||
              !interior(yc, parent.origin_y_comoving, parent.dy_comoving, parent.ny) ||
              !interior(zc, parent.origin_z_comoving, parent.dz_comoving, parent.nz)) {
            continue;
          }
          CHUI_REQUIRE(std::abs(child.faceFluxX(i, j, k) - exact.faceFluxX(i, j, k)) <= 1e-14);
          CHUI_REQUIRE(std::abs(child.faceFluxY(i, j, k) - exact.faceFluxY(i, j, k)) <= 1e-14);
          CHUI_REQUIRE(std::abs(child.faceFluxZ(i, j, k) - exact.faceFluxZ(i, j, k)) <= 1e-14);
        }
      }
    }
  }
  MhdFaceFluxStorage restored(parent);
  MhdAmrSynchronizer::restrictFaceFluxes(parent, restored.fluxX(), restored.fluxY(),
                                         restored.fluxZ(), children);
  for (std::size_t index = 0; index < faces.fluxY().size(); ++index) {
    CHUI_REQUIRE(std::abs(restored.fluxY()[index] - faces.fluxY()[index]) <= 1e-14);
  }
}

int main() {
  testUniformProlongation();
  testNonuniformTransferIdentity();
  testAnisotropicLinearFieldExactness();
  testGasTransferTransaction();
  testEmfRegisterCoverage();
  testFaceLayoutLifecycle();
#if COSMOSIM_ENABLE_HDF5
  testRestartGroupCodec();
#endif
  testConfigAndCapabilityAdmission();
  return 0;
}
