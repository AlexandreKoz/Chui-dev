#include "cosmosim/io/snapshot_hdf5.hpp"

#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/stellar_population.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "cosmosim/core/periodic_domain.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/version.hpp"
#include "cosmosim/io/io_contract.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "io/internal/snapshot_conversion.hpp"
#include "io/internal/snapshot_field_contract.hpp"
#include "io/internal/snapshot_readiness.hpp"
#include "io/internal/snapshot_set_internal.hpp"
#include "io/internal/transactional_file.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/governed_scratch_arena.hpp"
#include "io/internal/sidecar_row_lookup.hpp"
#include "core/internal/sha256.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif
#if COSMOSIM_ENABLE_MPI && COSMOSIM_HDF5_PARALLEL
#include <mpi.h>
#endif

namespace cosmosim::io {
namespace {

constexpr std::uint32_t k_species_dark_matter =
    static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
constexpr std::uint32_t k_species_gas = static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
constexpr std::uint32_t k_species_star = static_cast<std::uint32_t>(core::ParticleSpecies::kStar);
constexpr std::uint32_t k_species_black_hole =
    static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);
constexpr std::uint32_t k_species_tracer = static_cast<std::uint32_t>(core::ParticleSpecies::kTracer);

[[nodiscard]] std::size_t mapSpeciesTagToPartType(std::uint32_t species_tag) {
  if (species_tag == k_species_gas) {
    return 0;
  }
  if (species_tag == k_species_dark_matter) {
    return 1;
  }
  if (species_tag == k_species_tracer) {
    return 3;
  }
  if (species_tag == k_species_star) {
    return 4;
  }
  if (species_tag == k_species_black_hole) {
    return 5;
  }
  throw std::invalid_argument("snapshot HDF5: invalid particle species tag " + std::to_string(species_tag));
}

[[nodiscard]] std::uint32_t speciesEnumToTag(core::ParticleSpecies species) {
  return static_cast<std::uint32_t>(species);
}

[[nodiscard]] std::uint32_t mapPartTypeToSpeciesTag(
    std::size_t part_type,
    const SnapshotReadOptions& options,
    bool chui_authored) {
  if (part_type >= options.part_type_species_map.size()) {
    throw std::invalid_argument("snapshot HDF5: invalid PartType index " + std::to_string(part_type));
  }
  if (options.part_type_species_map[part_type].has_value()) {
    return speciesEnumToTag(*options.part_type_species_map[part_type]);
  }
  if (part_type == 0U) return k_species_gas;
  if (part_type == 1U) return k_species_dark_matter;
  if (part_type == 4U) return k_species_star;
  if (part_type == 5U) return k_species_black_hole;
  if (chui_authored && part_type == 3U) return k_species_tracer;
  throw std::runtime_error(
      "snapshot import: populated PartType" + std::to_string(part_type) +
      " is scientifically ambiguous for the selected dialect; provide an explicit part_type_species_map");
}

[[nodiscard]] std::string toTypeAliasPath(std::size_t type_index) {
  return "/ParticleType" + std::to_string(type_index);
}

#if COSMOSIM_ENABLE_HDF5

thread_local std::uint64_t g_snapshot_max_dataset_bytes =
    std::numeric_limits<std::uint64_t>::max();
thread_local std::uint64_t g_snapshot_max_attribute_bytes =
    std::numeric_limits<std::uint64_t>::max();

class SnapshotReadBudgetScope {
 public:
  explicit SnapshotReadBudgetScope(const SnapshotReadBudget& budget) noexcept
      : m_previous_dataset(g_snapshot_max_dataset_bytes),
        m_previous_attribute(g_snapshot_max_attribute_bytes) {
    g_snapshot_max_dataset_bytes = budget.max_dataset_bytes;
    g_snapshot_max_attribute_bytes = budget.max_attribute_bytes;
  }
  SnapshotReadBudgetScope(const SnapshotReadBudgetScope&) = delete;
  SnapshotReadBudgetScope& operator=(const SnapshotReadBudgetScope&) = delete;
  ~SnapshotReadBudgetScope() {
    g_snapshot_max_dataset_bytes = m_previous_dataset;
    g_snapshot_max_attribute_bytes = m_previous_attribute;
  }

 private:
  std::uint64_t m_previous_dataset;
  std::uint64_t m_previous_attribute;
};

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t handle = -1) : m_handle(handle) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept : m_handle(other.m_handle) { other.m_handle = -1; }
  Hdf5Handle& operator=(Hdf5Handle&& other) noexcept {
    if (this != &other) {
      close();
      m_handle = other.m_handle;
      other.m_handle = -1;
    }
    return *this;
  }
  ~Hdf5Handle() { close(); }

  [[nodiscard]] hid_t get() const { return m_handle; }
  [[nodiscard]] bool valid() const { return m_handle >= 0; }

 private:
  void close() {
    if (m_handle >= 0) {
      const H5I_type_t type = H5Iget_type(m_handle);
      if (type == H5I_FILE) {
        H5Fclose(m_handle);
      } else if (type == H5I_GROUP) {
        H5Gclose(m_handle);
      } else if (type == H5I_DATASET) {
        H5Dclose(m_handle);
      } else if (type == H5I_DATASPACE) {
        H5Sclose(m_handle);
      } else if (type == H5I_ATTR) {
        H5Aclose(m_handle);
      } else if (type == H5I_DATATYPE) {
        H5Tclose(m_handle);
      } else if (type == H5I_GENPROP_LST) {
        H5Pclose(m_handle);
      }
      m_handle = -1;
    }
  }

  hid_t m_handle = -1;
};

[[nodiscard]] bool hdf5PathExists(hid_t parent, const std::string& path) {
  return H5Lexists(parent, path.c_str(), H5P_DEFAULT) > 0;
}

[[nodiscard]] bool hdf5AttributeExists(hid_t location, const std::string& key) {
  const htri_t exists = H5Aexists(location, key.c_str());
  return exists > 0;
}

void writeScalarStringAttribute(hid_t location, const std::string& key, std::string_view value) {
  if (static_cast<std::uint64_t>(value.size()) > g_snapshot_max_attribute_bytes) {
    throw std::length_error("snapshot string attribute exceeds max_attribute_bytes: " + key);
  }
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  if (!scalar_space.valid()) {
    throw std::runtime_error("failed to allocate scalar dataspace for attribute: " + key);
  }

  Hdf5Handle string_type(H5Tcopy(H5T_C_S1));
  if (!string_type.valid()) {
    throw std::runtime_error("failed to allocate HDF5 string datatype");
  }
  // HDF5 1.10 stores large fixed-width string attributes in the object header
  // and can fail with "object header message is too large" for runtime
  // provenance blobs.  Variable-length scalar strings keep the attribute
  // metadata small while preserving the same logical schema.  The reader
  // already accepts both fixed-width and variable-length HDF5 strings.
  constexpr std::size_t k_variable_string_threshold = 4096U;
  const bool variable_length = value.size() >= k_variable_string_threshold;
  if (variable_length) {
    if (H5Tset_size(string_type.get(), H5T_VARIABLE) < 0) {
      throw std::runtime_error("failed to configure variable-length HDF5 string datatype");
    }
  } else if (H5Tset_size(
                 string_type.get(),
                 std::max<std::size_t>(std::size_t{1}, value.size() + 1)) < 0 ||
             H5Tset_strpad(string_type.get(), H5T_STR_NULLTERM) < 0) {
    throw std::runtime_error("failed to configure HDF5 string datatype");
  }

  Hdf5Handle attr(H5Acreate2(location, key.c_str(), string_type.get(), scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid()) {
    throw std::runtime_error("failed to create attribute: " + key);
  }

  std::string local(value);
  if (variable_length) {
    const char* raw = local.c_str();
    if (H5Awrite(attr.get(), string_type.get(), &raw) < 0) {
      throw std::runtime_error("failed to write variable string attribute: " + key);
    }
  } else if (H5Awrite(attr.get(), string_type.get(), local.c_str()) < 0) {
    throw std::runtime_error("failed to write attribute: " + key);
  }
}

void writeScalarDoubleAttribute(hid_t location, const std::string& key, double value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, key.c_str(), H5T_IEEE_F64LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error("failed to write double attribute: " + key);
  }
}

void writeScalarUint32Attribute(hid_t location, const std::string& key, std::uint32_t value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, key.c_str(), H5T_STD_U32LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error("failed to write uint32 attribute: " + key);
  }
}

void writeScalarUint64Attribute(hid_t location, const std::string& key, std::uint64_t value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, key.c_str(), H5T_STD_U64LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_UINT64, &value) < 0) {
    throw std::runtime_error("failed to write uint64 attribute: " + key);
  }
}

[[nodiscard]] const char* metalSpeciesModeLabel(core::MetalSpeciesMode mode) {
  switch (mode) {
    case core::MetalSpeciesMode::kTotalOnly:
      return "total_only";
    case core::MetalSpeciesMode::kCoreElements:
      return "core_elements";
  }
  throw std::logic_error("unhandled metal species mode");
}

[[nodiscard]] const char* metalDiffusionModelLabel(core::MetalDiffusionModel model) {
  switch (model) {
    case core::MetalDiffusionModel::kNone:
      return "none";
    case core::MetalDiffusionModel::kSmagorinsky:
      return "smagorinsky";
  }
  throw std::logic_error("unhandled metal diffusion model");
}

[[nodiscard]] const char* metalDiffusionIntegratorLabel(
    core::MetalDiffusionTimeIntegrator integrator) {
  switch (integrator) {
    case core::MetalDiffusionTimeIntegrator::kExplicitSubcycling:
      return "explicit_subcycling";
    case core::MetalDiffusionTimeIntegrator::kRkl2:
      return "rkl2";
  }
  throw std::logic_error("unhandled metal diffusion integrator");
}

void writeHeaderArrays(
    hid_t header_group,
    const std::array<std::uint64_t, 6>& local_part_count,
    const std::array<std::uint64_t, 6>& global_part_count,
    const std::array<double, 6>& mass_table,
    const core::SimulationConfig& config,
    const core::SimulationState& state,
    const SnapshotSetMemberInfo& member,
    const internal::SnapshotConversionContext& conversion) {
  std::array<std::uint32_t, 6> local_low{};
  std::array<std::uint32_t, 6> global_low{};
  std::array<std::uint32_t, 6> global_high{};
  for (std::size_t i = 0; i < local_part_count.size(); ++i) {
    if (local_part_count[i] > static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max())) {
      throw std::runtime_error(
          "snapshot writer: one file cannot contain more than UINT32_MAX rows per PartType; "
          "increase NumFilesPerSnapshot / repartition the output set");
    }
    local_low[i] = static_cast<std::uint32_t>(local_part_count[i]);
    global_low[i] = static_cast<std::uint32_t>(global_part_count[i] & 0xffffffffULL);
    global_high[i] = static_cast<std::uint32_t>(global_part_count[i] >> 32U);
  }

  hsize_t vector_dims[1] = {6};
  Hdf5Handle vector_space(H5Screate_simple(1, vector_dims, nullptr));
  if (!vector_space.valid()) {
    throw std::runtime_error("failed to create vector dataspace for Header attributes");
  }

  Hdf5Handle attr_num_this(
      H5Acreate2(header_group, "NumPart_ThisFile", H5T_STD_U32LE, vector_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle attr_num_total(
      H5Acreate2(header_group, "NumPart_Total", H5T_STD_U32LE, vector_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle attr_num_total_high(H5Acreate2(
      header_group, "NumPart_Total_HighWord", H5T_STD_U32LE, vector_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle attr_mass(
      H5Acreate2(header_group, "MassTable", H5T_IEEE_F64LE, vector_space.get(), H5P_DEFAULT, H5P_DEFAULT));

  if (!attr_num_this.valid() || H5Awrite(attr_num_this.get(), H5T_NATIVE_UINT32, local_low.data()) < 0 ||
      !attr_num_total.valid() || H5Awrite(attr_num_total.get(), H5T_NATIVE_UINT32, global_low.data()) < 0 ||
      !attr_num_total_high.valid() || H5Awrite(attr_num_total_high.get(), H5T_NATIVE_UINT32, global_high.data()) < 0 ||
      !attr_mass.valid() || H5Awrite(attr_mass.get(), H5T_NATIVE_DOUBLE, mass_table.data()) < 0) {
    throw std::runtime_error("failed to write required Header vector attributes");
  }

  writeScalarDoubleAttribute(header_group, "Time", state.metadata.scale_factor);
  writeScalarDoubleAttribute(header_group, "Redshift", (1.0 / state.metadata.scale_factor) - 1.0);
  writeScalarDoubleAttribute(
      header_group, "BoxSize",
      conversion.boxSizeMpcToStored(config.cosmology.box_size_mpc_comoving));

  // CHUI extensions preserve the exact axis-aware box in physical Mpc/h-independent
  // semantic metadata while the canonical BoxSize follows the selected export dialect.
  writeScalarDoubleAttribute(header_group, "CHUIBoxSizeX_MpcComoving", config.cosmology.box_size_x_mpc_comoving);
  writeScalarDoubleAttribute(header_group, "CHUIBoxSizeY_MpcComoving", config.cosmology.box_size_y_mpc_comoving);
  writeScalarDoubleAttribute(header_group, "CHUIBoxSizeZ_MpcComoving", config.cosmology.box_size_z_mpc_comoving);
  {
    const std::array<double, 3> box_size_vec = {
        config.cosmology.box_size_x_mpc_comoving,
        config.cosmology.box_size_y_mpc_comoving,
        config.cosmology.box_size_z_mpc_comoving};
    hsize_t dims[1] = {3};
    Hdf5Handle vec_space(H5Screate_simple(1, dims, nullptr));
    Hdf5Handle vec_attr(H5Acreate2(
        header_group, "CHUIBoxSizeVec_MpcComoving", H5T_IEEE_F64LE,
        vec_space.get(), H5P_DEFAULT, H5P_DEFAULT));
    if (!vec_space.valid() || !vec_attr.valid() ||
        H5Awrite(vec_attr.get(), H5T_NATIVE_DOUBLE, box_size_vec.data()) < 0) {
      throw std::runtime_error("failed to write axis-aware BoxSize header attributes");
    }
  }
  writeScalarDoubleAttribute(header_group, "Omega0", config.cosmology.omega_matter);
  writeScalarDoubleAttribute(header_group, "OmegaLambda", config.cosmology.omega_lambda);
  writeScalarDoubleAttribute(header_group, "OmegaBaryon", config.cosmology.omega_baryon);
  writeScalarDoubleAttribute(header_group, "HubbleParam", config.cosmology.hubble_param);

  writeScalarUint32Attribute(header_group, "NumFilesPerSnapshot", member.num_files_per_snapshot);
  writeScalarUint32Attribute(header_group, "CHUISnapshotMemberIndex", member.member_index);
  if (!member.generation_id.empty()) {
    writeScalarStringAttribute(header_group, "CHUISnapshotGenerationID", member.generation_id);
  }
  writeScalarUint32Attribute(header_group, "Flag_Sfr", config.physics.enable_star_formation ? 1u : 0u);
  writeScalarUint32Attribute(header_group, "Flag_Cooling", config.physics.enable_cooling ? 1u : 0u);
  writeScalarUint32Attribute(header_group, "Flag_StellarAge", 1u);
  writeScalarUint32Attribute(header_group, "Flag_Metals", 1u);
  writeScalarUint32Attribute(header_group, "Flag_Feedback", config.physics.enable_feedback ? 1u : 0u);
}

Hdf5Handle createDatasetProperties(std::size_t element_count, std::size_t element_width, const SnapshotIoPolicy& policy) {
  Hdf5Handle properties(H5Pcreate(H5P_DATASET_CREATE));
  if (!properties.valid()) {
    throw std::runtime_error("failed to create HDF5 dataset creation properties");
  }

  const hsize_t chunk_rows = std::max<hsize_t>(1, std::min<hsize_t>(static_cast<hsize_t>(policy.chunk_particle_count), static_cast<hsize_t>(element_count)));
  if (element_width == 1) {
    hsize_t chunk_dims[1] = {chunk_rows};
    if (H5Pset_chunk(properties.get(), 1, chunk_dims) < 0) {
      throw std::runtime_error("failed to set chunking for 1D dataset");
    }
  } else {
    hsize_t chunk_dims[2] = {chunk_rows, static_cast<hsize_t>(element_width)};
    if (H5Pset_chunk(properties.get(), 2, chunk_dims) < 0) {
      throw std::runtime_error("failed to set chunking for 2D dataset");
    }
  }

  if (policy.enable_compression) {
    if (H5Pset_deflate(properties.get(), static_cast<unsigned int>(policy.compression_level)) < 0) {
      throw std::runtime_error("failed to set deflate compression level");
    }
  }

  return properties;
}

struct StreamingDataset {
  Hdf5Handle handle;
  std::size_t width = 1;
};

[[nodiscard]] StreamingDataset createStreamingDataset(
    hid_t group,
    std::string_view name,
    hid_t file_type,
    std::size_t rows,
    std::size_t width,
    const SnapshotIoPolicy& policy) {
  // CHUI Part 3: any positive width (composition lanes are [N, lanes]).
  if (rows == 0 || width == 0) {
    throw std::invalid_argument("invalid streaming dataset shape");
  }
  Hdf5Handle properties = createDatasetProperties(rows, width, policy);
  Hdf5Handle dataspace;
  if (width == 1) {
    hsize_t dims[1] = {static_cast<hsize_t>(rows)};
    dataspace = Hdf5Handle(H5Screate_simple(1, dims, nullptr));
  } else {
    hsize_t dims[2] = {static_cast<hsize_t>(rows), static_cast<hsize_t>(width)};
    dataspace = Hdf5Handle(H5Screate_simple(2, dims, nullptr));
  }
  if (!dataspace.valid()) {
    throw std::runtime_error("failed to create streaming dataspace for " + std::string(name));
  }
  StreamingDataset result;
  result.width = width;
  result.handle = Hdf5Handle(H5Dcreate2(
      group, std::string(name).c_str(), file_type, dataspace.get(), H5P_DEFAULT,
      properties.get(), H5P_DEFAULT));
  if (!result.handle.valid()) {
    throw std::runtime_error("failed to create streaming dataset: " + std::string(name));
  }
  return result;
}

void writeStreamingChunk(
    const StreamingDataset& dataset,
    hid_t native_type,
    const void* values,
    std::size_t row_offset,
    std::size_t row_count) {
  if (row_count == 0) {
    return;
  }
  Hdf5Handle file_space(H5Dget_space(dataset.handle.get()));
  if (!file_space.valid()) {
    throw std::runtime_error("failed to open streaming dataset file space");
  }
  Hdf5Handle mem_space;
  if (dataset.width == 1) {
    hsize_t start[1] = {static_cast<hsize_t>(row_offset)};
    hsize_t count[1] = {static_cast<hsize_t>(row_count)};
    if (H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, start, nullptr, count, nullptr) < 0) {
      throw std::runtime_error("failed to select 1D output hyperslab");
    }
    mem_space = Hdf5Handle(H5Screate_simple(1, count, nullptr));
  } else {
    hsize_t start[2] = {static_cast<hsize_t>(row_offset), 0};
    hsize_t count[2] = {static_cast<hsize_t>(row_count), static_cast<hsize_t>(dataset.width)};
    if (H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, start, nullptr, count, nullptr) < 0) {
      throw std::runtime_error("failed to select 2D output hyperslab");
    }
    mem_space = Hdf5Handle(H5Screate_simple(2, count, nullptr));
  }
  if (!mem_space.valid() ||
      H5Dwrite(dataset.handle.get(), native_type, mem_space.get(), file_space.get(), H5P_DEFAULT, values) < 0) {
    throw std::runtime_error("failed to write streaming dataset chunk");
  }
}

struct SnapshotWriteMetrics {
  std::uint64_t peak_staging_bytes = 0;
  std::uint64_t logical_bytes_written = 0;
  std::uint64_t chunk_write_count = 0;

  void observe(std::uint64_t staging_bytes, std::uint64_t logical_bytes) {
    peak_staging_bytes = std::max(peak_staging_bytes, staging_bytes);
    logical_bytes_written += logical_bytes;
    ++chunk_write_count;
  }
};

void validateDatasetShape(
    hid_t dataset,
    int expected_rank,
    std::size_t start,
    std::size_t count,
    std::optional<hsize_t> expected_second_dimension,
    std::string_view dataset_name) {
  if (dataset < 0) {
    throw std::runtime_error("failed to open dataset: " + std::string(dataset_name));
  }
  Hdf5Handle space(H5Dget_space(dataset));
  if (!space.valid()) {
    throw std::runtime_error("failed to inspect dataset shape: " + std::string(dataset_name));
  }
  const int rank = H5Sget_simple_extent_ndims(space.get());
  if (rank != expected_rank) {
    throw std::runtime_error(
        "dataset rank mismatch for " + std::string(dataset_name) +
        ": expected " + std::to_string(expected_rank) +
        ", found " + std::to_string(rank));
  }
  std::array<hsize_t, 2> dims{};
  if (H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr) < 0) {
    throw std::runtime_error("failed to read dataset dimensions: " + std::string(dataset_name));
  }
  if (start > static_cast<std::size_t>(dims[0]) ||
      count > static_cast<std::size_t>(dims[0]) - start) {
    throw std::runtime_error("dataset row extent is smaller than Header count: " + std::string(dataset_name));
  }
  if (expected_second_dimension.has_value() && dims[1] != *expected_second_dimension) {
    throw std::runtime_error(
        "dataset component extent mismatch for " + std::string(dataset_name));
  }
  Hdf5Handle datatype(H5Dget_type(dataset));
  if (!datatype.valid()) {
    throw std::runtime_error(
        "failed to inspect dataset datatype: " + std::string(dataset_name));
  }
  const std::size_t scalar_bytes = H5Tget_size(datatype.get());
  if (scalar_bytes == 0U) {
    throw std::runtime_error(
        "dataset has zero-size datatype: " + std::string(dataset_name));
  }
  std::uint64_t element_count = 1U;
  for (int axis = 0; axis < rank; ++axis) {
    const std::uint64_t extent = static_cast<std::uint64_t>(dims[static_cast<std::size_t>(axis)]);
    if (extent != 0U &&
        element_count > std::numeric_limits<std::uint64_t>::max() / extent) {
      throw std::length_error(
          "dataset element count overflows uint64: " + std::string(dataset_name));
    }
    element_count *= extent;
  }
  const std::uint64_t scalar_bytes_u64 = static_cast<std::uint64_t>(scalar_bytes);
  if (element_count != 0U &&
      scalar_bytes_u64 > g_snapshot_max_dataset_bytes / element_count) {
    throw std::length_error(
        "snapshot dataset exceeds max_dataset_bytes read budget: " +
        std::string(dataset_name));
  }
}

void readDatasetChunk1d(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<double>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  validateDatasetShape(dataset.get(), 1, start, count, std::nullopt, dataset_name);
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  hsize_t file_offset[1] = {static_cast<hsize_t>(start)};
  hsize_t file_count[1] = {static_cast<hsize_t>(count)};
  H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);
  hsize_t mem_dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle mem_space(H5Screate_simple(1, mem_dims, nullptr));
  out.resize(count);
  if (!dataset.valid() || H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, mem_space.get(), file_space.get(), H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed to read 1D dataset: " + dataset_name);
  }
}


void readDatasetChunk1dU8(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<std::uint8_t>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  validateDatasetShape(dataset.get(), 1, start, count, std::nullopt, dataset_name);
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  hsize_t file_offset[1] = {static_cast<hsize_t>(start)};
  hsize_t file_count[1] = {static_cast<hsize_t>(count)};
  H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);
  hsize_t mem_dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle mem_space(H5Screate_simple(1, mem_dims, nullptr));
  out.resize(count);
  if (!dataset.valid() || H5Dread(dataset.get(), H5T_NATIVE_UINT8, mem_space.get(), file_space.get(), H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed to read 1D uint8 dataset: " + dataset_name);
  }
}
void readDatasetChunk2d(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<double>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  validateDatasetShape(dataset.get(), 2, start, count, 3, dataset_name);
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  hsize_t file_offset[2] = {static_cast<hsize_t>(start), 0};
  hsize_t file_count[2] = {static_cast<hsize_t>(count), 3};
  H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);
  hsize_t mem_dims[2] = {static_cast<hsize_t>(count), 3};
  Hdf5Handle mem_space(H5Screate_simple(2, mem_dims, nullptr));
  out.resize(count * 3);
  if (!dataset.valid() || H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, mem_space.get(), file_space.get(), H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed to read 2D dataset: " + dataset_name);
  }
}

void readDatasetChunkIds(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<std::uint64_t>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  validateDatasetShape(dataset.get(), 1, start, count, std::nullopt, dataset_name);
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  hsize_t file_offset[1] = {static_cast<hsize_t>(start)};
  hsize_t file_count[1] = {static_cast<hsize_t>(count)};
  H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);
  hsize_t mem_dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle mem_space(H5Screate_simple(1, mem_dims, nullptr));
  out.resize(count);
  if (!dataset.valid() || H5Dread(dataset.get(), H5T_NATIVE_UINT64, mem_space.get(), file_space.get(), H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed to read ids dataset: " + dataset_name);
  }
}

void readDatasetChunkU32(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<std::uint32_t>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  validateDatasetShape(dataset.get(), 1, start, count, std::nullopt, dataset_name);
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  hsize_t file_offset[1] = {static_cast<hsize_t>(start)};
  hsize_t file_count[1] = {static_cast<hsize_t>(count)};
  H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);
  hsize_t mem_dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle mem_space(H5Screate_simple(1, mem_dims, nullptr));
  out.resize(count);
  if (!dataset.valid() || H5Dread(dataset.get(), H5T_NATIVE_UINT32, mem_space.get(), file_space.get(), H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed to read uint32 dataset: " + dataset_name);
  }
}


void observeDatasetStoragePolicy(
    hid_t group,
    const std::string& dataset_name,
    SnapshotIoReport& report) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    return;
  }
  Hdf5Handle properties(H5Dget_create_plist(dataset.get()));
  if (!properties.valid()) {
    return;
  }
  std::size_t chunk_rows = 0U;
  const H5D_layout_t layout = H5Pget_layout(properties.get());
  if (layout == H5D_CHUNKED) {
    Hdf5Handle space(H5Dget_space(dataset.get()));
    const int rank = space.valid() ? H5Sget_simple_extent_ndims(space.get()) : -1;
    if (rank == 1 || rank == 2) {
      std::array<hsize_t, 2> chunk{};
      if (H5Pget_chunk(properties.get(), rank, chunk.data()) >= 0) {
        chunk_rows = static_cast<std::size_t>(chunk[0]);
      }
    }
  }
  bool compressed = false;
  int compression_level = 0;
  const int filter_count = H5Pget_nfilters(properties.get());
  for (int filter_index = 0; filter_index < filter_count; ++filter_index) {
    unsigned int flags = 0U;
    std::array<unsigned int, 16> client_data{};
    std::size_t client_data_count = client_data.size();
    const H5Z_filter_t filter = H5Pget_filter2(
        properties.get(), static_cast<unsigned int>(filter_index), &flags,
        &client_data_count, client_data.data(), 0, nullptr, nullptr);
    if (filter == H5Z_FILTER_DEFLATE) {
      compressed = true;
      if (client_data_count > 0U) {
        compression_level = static_cast<int>(client_data[0]);
      }
    }
  }
  if (!report.storage_policy_known) {
    report.storage_policy_known = true;
    report.chunk_particle_count = chunk_rows;
    report.compression_enabled = compressed;
    report.compression_level = compression_level;
    return;
  }
  if (report.chunk_particle_count != chunk_rows ||
      report.compression_enabled != compressed ||
      report.compression_level != compression_level) {
    report.storage_policy_mixed = true;
  }
}

[[nodiscard]] std::string pickAlias(
    hid_t group,
    const GadgetArepoFieldAliases& aliases,
    SnapshotIoReport& report,
    const std::string& context_path,
    bool required) {
  for (std::string_view alias : aliases.read_aliases) {
    if (hdf5PathExists(group, std::string(alias))) {
      report.present_aliases.push_back(context_path + "=" + std::string(alias));
      return std::string(alias);
    }
  }
  if (required) {
    throw std::runtime_error("required dataset missing under: " + context_path);
  }
  return {};
}

void readScalarStringAttribute(
    hid_t location,
    const std::string& key,
    std::string& out_value) {
  if (H5Aexists(location, key.c_str()) <= 0) return;
  Hdf5Handle attr(H5Aopen(location, key.c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    return;
  }
  Hdf5Handle type(H5Aget_type(attr.get()));
  if (!type.valid() || H5Tget_class(type.get()) != H5T_STRING) {
    throw std::runtime_error("snapshot attribute is not a string: " + key);
  }
  if (H5Tis_variable_str(type.get()) > 0) {
    char* raw = nullptr;
    if (H5Aread(attr.get(), type.get(), &raw) < 0) {
      throw std::runtime_error("failed reading variable string attribute: " + key);
    }
    const std::size_t length = raw != nullptr ? std::char_traits<char>::length(raw) : 0U;
    if (length > g_snapshot_max_attribute_bytes) {
      if (raw != nullptr) H5free_memory(raw);
      throw std::length_error(
          "snapshot string attribute exceeds max_attribute_bytes: " + key);
    }
    out_value.assign(raw != nullptr ? raw : "", length);
    if (raw != nullptr) H5free_memory(raw);
    return;
  }
  const std::size_t length = H5Tget_size(type.get());
  if (length == 0U ||
      static_cast<std::uint64_t>(length) > g_snapshot_max_attribute_bytes) {
    throw std::length_error(
        "snapshot string attribute exceeds max_attribute_bytes: " + key);
  }
  std::string buffer(length, '\0');
  if (H5Aread(attr.get(), type.get(), buffer.data()) < 0) {
    throw std::runtime_error("failed reading attribute: " + key);
  }
  const std::size_t terminator = buffer.find('\0');
  if (terminator != std::string::npos) {
    buffer.resize(terminator);
  }
  out_value = std::move(buffer);
}

void requireLogicalScalarAttributeShape(hid_t attribute, const std::string& key) {
  Hdf5Handle space(H5Aget_space(attribute));
  if (!space.valid()) throw std::runtime_error("failed reading attribute dataspace: " + key);
  const int rank = H5Sget_simple_extent_ndims(space.get());
  if (rank == 0) return;
  if (rank == 1) {
    hsize_t extent = 0;
    if (H5Sget_simple_extent_dims(space.get(), &extent, nullptr) < 0 || extent != 1U) {
      throw std::runtime_error("logical scalar attribute must have rank 0 or shape [1]: " + key);
    }
    return;
  }
  throw std::runtime_error("logical scalar attribute must have rank 0 or shape [1]: " + key);
}

[[nodiscard]] bool readScalarUint32Attribute(hid_t location, const std::string& key, std::uint32_t& out_value) {
  if (H5Aexists(location, key.c_str()) <= 0) return false;
  Hdf5Handle attr(H5Aopen(location, key.c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    return false;
  }
  requireLogicalScalarAttributeShape(attr.get(), key);
  if (H5Aread(attr.get(), H5T_NATIVE_UINT32, &out_value) < 0) {
    throw std::runtime_error("failed reading attribute: " + key);
  }
  return true;
}

[[nodiscard]] bool readScalarUint64Attribute(hid_t location, const std::string& key, std::uint64_t& out_value) {
  if (H5Aexists(location, key.c_str()) <= 0) return false;
  Hdf5Handle attr(H5Aopen(location, key.c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    return false;
  }
  requireLogicalScalarAttributeShape(attr.get(), key);
  if (H5Aread(attr.get(), H5T_NATIVE_UINT64, &out_value) < 0) {
    throw std::runtime_error("failed reading attribute: " + key);
  }
  return true;
}

[[nodiscard]] bool readScalarDoubleAttribute(hid_t location, const std::string& key, double& out_value) {
  if (H5Aexists(location, key.c_str()) <= 0) return false;
  Hdf5Handle attr(H5Aopen(location, key.c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    return false;
  }
  requireLogicalScalarAttributeShape(attr.get(), key);
  if (H5Aread(attr.get(), H5T_NATIVE_DOUBLE, &out_value) < 0) {
    throw std::runtime_error("failed reading attribute: " + key);
  }
  return true;
}

void readHeaderCountArray(
    hid_t header_group,
    const char* attribute_name,
    std::array<std::uint64_t, 6>& out_counts) {
  if (H5Aexists(header_group, attribute_name) <= 0) {
    throw std::runtime_error(
        std::string("Header/") + attribute_name + " missing");
  }
  Hdf5Handle attr(H5Aopen(header_group, attribute_name, H5P_DEFAULT));
  Hdf5Handle type(attr.valid() ? H5Aget_type(attr.get()) : -1);
  if (!attr.valid() || !type.valid()) {
    throw std::runtime_error(
        std::string("Header/") + attribute_name + " unreadable");
  }
  const std::size_t element_bytes = H5Tget_size(type.get());
  if (element_bytes >= sizeof(std::uint64_t)) {
    if (H5Aread(attr.get(), H5T_NATIVE_UINT64, out_counts.data()) < 0) {
      throw std::runtime_error(
          std::string("Header/") + attribute_name + " unreadable as uint64");
    }
    return;
  }
  std::array<std::uint32_t, 6> raw{};
  if (H5Aread(attr.get(), H5T_NATIVE_UINT32, raw.data()) < 0) {
    throw std::runtime_error(
        std::string("Header/") + attribute_name + " unreadable as uint32");
  }
  for (std::size_t i = 0; i < out_counts.size(); ++i) {
    out_counts[i] = raw[i];
  }
}

void readHeaderCounts(
    hid_t header_group,
    std::array<std::uint64_t, 6>& out_counts) {
  readHeaderCountArray(header_group, "NumPart_ThisFile", out_counts);
}

void readHeaderGlobalCounts(
    hid_t header_group,
    const std::array<std::uint64_t, 6>& fallback_local,
    std::array<std::uint64_t, 6>& out_counts) {
  if (H5Aexists(header_group, "NumPart_Total") <= 0) {
    out_counts = fallback_local;
    return;
  }
  Hdf5Handle low_attr(H5Aopen(header_group, "NumPart_Total", H5P_DEFAULT));
  Hdf5Handle low_type(low_attr.valid() ? H5Aget_type(low_attr.get()) : -1);
  if (!low_attr.valid() || !low_type.valid()) {
    throw std::runtime_error("Header/NumPart_Total unreadable");
  }
  if (H5Tget_size(low_type.get()) >= sizeof(std::uint64_t)) {
    if (H5Aread(low_attr.get(), H5T_NATIVE_UINT64, out_counts.data()) < 0) {
      throw std::runtime_error("Header/NumPart_Total unreadable as uint64");
    }
    return;
  }
  std::array<std::uint32_t, 6> low{};
  std::array<std::uint32_t, 6> high{};
  if (H5Aread(low_attr.get(), H5T_NATIVE_UINT32, low.data()) < 0) {
    throw std::runtime_error("Header/NumPart_Total unreadable as uint32");
  }
  if (H5Aexists(header_group, "NumPart_Total_HighWord") > 0) {
    Hdf5Handle high_attr(
        H5Aopen(header_group, "NumPart_Total_HighWord", H5P_DEFAULT));
    if (!high_attr.valid() ||
        H5Aread(high_attr.get(), H5T_NATIVE_UINT32, high.data()) < 0) {
      throw std::runtime_error("Header/NumPart_Total_HighWord unreadable");
    }
  }
  for (std::size_t i = 0; i < out_counts.size(); ++i) {
    out_counts[i] = static_cast<std::uint64_t>(low[i]) |
        (static_cast<std::uint64_t>(high[i]) << 32U);
  }
}

[[nodiscard]] std::size_t checkedSnapshotCountToSize(
    std::uint64_t value,
    std::string_view context) {
  if (value > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw std::length_error(std::string(context) + ": count exceeds addressable size_t range");
  }
  return static_cast<std::size_t>(value);
}

void readOptionalHeaderMassTable(
    hid_t header_group,
    std::array<double, 6>& out_mass_table) {
  out_mass_table.fill(0.0);
  if (H5Aexists(header_group, "MassTable") <= 0) return;
  Hdf5Handle attr_mass(H5Aopen(header_group, "MassTable", H5P_DEFAULT));
  if (attr_mass.valid()) {
    if (H5Aread(attr_mass.get(), H5T_NATIVE_DOUBLE, out_mass_table.data()) < 0) {
      throw std::runtime_error("Header/MassTable unreadable");
    }
  }
}

void readOptionalHeaderDouble(hid_t header_group, const std::string& key, double default_value, double& out_value) {
  out_value = default_value;
  if (H5Aexists(header_group, key.c_str()) <= 0) return;
  Hdf5Handle attr(H5Aopen(header_group, key.c_str(), H5P_DEFAULT));
  if (attr.valid()) {
    requireLogicalScalarAttributeShape(attr.get(), key);
    if (H5Aread(attr.get(), H5T_NATIVE_DOUBLE, &out_value) < 0) {
      throw std::runtime_error("failed reading optional Header attribute: " + key);
    }
  }
}

#endif

}  // namespace

const ScienceSnapshotSchemaMap& scienceSnapshotSchemaMap() {
  static const ScienceSnapshotSchemaMap k_schema{};
  return k_schema;
}

const GadgetArepoSchemaMap& gadgetArepoSchemaMap() {
  return scienceSnapshotSchemaMap();
}

void SnapshotReadResult::requireEvolutionReady() const {
  if (!report.evolution_ready) {
    std::string message =
        "snapshot read result is analysis-only/partial and is not safe to evolve";
    if (!report.evolution_readiness_reasons.empty()) {
      message += ": ";
      for (std::size_t i = 0; i < report.evolution_readiness_reasons.size(); ++i) {
        if (i != 0U) message += "; ";
        message += report.evolution_readiness_reasons[i];
      }
    }
    throw std::runtime_error(message);
  }
}

void writeScienceSnapshotHdf5(
    const std::filesystem::path& output_path,
    const SnapshotWritePayload& payload,
    const SnapshotIoPolicy& policy) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(output_path);
  static_cast<void>(payload);
  static_cast<void>(policy);
  throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: snapshot writer unavailable");
#else
  if (payload.state == nullptr || payload.config == nullptr) {
    throw std::runtime_error("snapshot writer requires non-null state and config");
  }
  if (policy.chunk_particle_count == 0) {
    throw std::runtime_error("chunk_particle_count must be > 0");
  }
  if (policy.enable_compression && (policy.compression_level < 0 || policy.compression_level > 9)) {
    throw std::runtime_error("compression_level must be in [0, 9]");
  }
  if (payload.normalized_config_text.empty()) {
    throw std::runtime_error(
        "snapshot writer requires the exact normalized configuration text for provenance");
  }
  if (payload.provenance.config_hash_hex.empty()) {
    throw std::runtime_error(
        "snapshot writer requires a non-empty provenance configuration hash");
  }
  if (payload.collective_single_file) {
#if !(COSMOSIM_ENABLE_MPI && COSMOSIM_HDF5_PARALLEL)
    throw std::runtime_error(
        "collective single-file science snapshots require MPI-enabled Parallel HDF5");
#else
    if (!policy.caller_managed_publication) {
      throw std::runtime_error(
          "collective single-file science snapshots require caller-managed transactional publication");
    }
    if (policy.enable_compression) {
      throw std::runtime_error(
          "collective single-file science snapshots do not enable filtered compression in the qualified backend");
    }
#endif
  }

  const core::SimulationState& state = *payload.state;
  const core::SimulationConfig& config = *payload.config;
  std::optional<physics::EffectiveMultiphaseEosTable> effective_eos_table;
  if (config.physics.enable_star_formation &&
      config.physics.star_formation_model ==
          core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
    const core::UnitSystem units = core::makeUnitSystem(
        config.units.length_unit, config.units.mass_unit, config.units.velocity_unit);
    effective_eos_table.emplace(
        physics::makeEffectiveMultiphaseEosConfig(config.physics),
        units,
        physics::makeEffectiveIsmReferenceCoolingProvider(config.physics));
  }
  state.requireGasCellIdentityMapCoversDenseRows("snapshot writer");
  if (state.cells.size() != state.gas_cells.size()) {
    throw std::runtime_error(
        "snapshot writer: CellSoa and GasCellSidecar row counts must match");
  }

  if (!state.hasCertifiedPersistentParticleIds()) {
    core::MemoryReservation particle_id_validation_reservation;
    if (payload.memory_governor != nullptr) {
      particle_id_validation_reservation = payload.memory_governor->reserve(
          core::MemoryClass::kDiagnostic,
          core::OwnershipValidationWorkspace::k_id_validation_workspace_limit_bytes,
          "io.snapshot.particle_id_validation");
      particle_id_validation_reservation.commit();
    }
    {
      core::OwnershipValidationWorkspace particle_id_validation_scratch;
      if (!state.validatePersistentParticleIds(particle_id_validation_scratch)) {
        throw std::runtime_error(
            "snapshot writer: persistent particle IDs must be nonzero and unique");
      }
    }
  } else if (!state.validatePersistentParticleIds()) {
    throw std::runtime_error(
        "snapshot writer: persistent particle IDs must be nonzero and unique");
  }
  const std::uint64_t lookup_bytes = internal::sidecarLookupWorkspaceBytes(
      state.star_particles.size(), state.tracers.size(), state.black_holes.size());
  core::GovernedScratchArena sidecar_lookup_arena(
      payload.memory_governor, core::MemoryClass::kDiagnostic,
      lookup_bytes, "io.snapshot.sidecar_lookup");
  internal::SidecarRowLookup star_row_by_particle(
      state.star_particles.particle_index, sidecar_lookup_arena.resource(), "stellar");
  internal::SidecarRowLookup tracer_row_by_particle(
      state.tracers.particle_index, sidecar_lookup_arena.resource(), "tracer");
  internal::SidecarRowLookup black_hole_row_by_particle(
      state.black_holes.particle_index, sidecar_lookup_arena.resource(), "black-hole");

  std::array<std::uint64_t, 6> count_by_type{};
  std::array<double, 6> mass_table{};

  count_by_type[0] = static_cast<std::uint64_t>(state.cells.size());
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const std::size_t type_index =
        mapSpeciesTagToPartType(state.particleSpeciesTag(i));
    if (type_index != 0U) {
      ++count_by_type[type_index];
    }
  }

  const SnapshotDialect write_dialect =
      internal::resolveSnapshotWriteDialect(policy.dialect, config);
  const internal::SnapshotConversionContext conversion =
      internal::makeSnapshotConversionContext(
          write_dialect, config, state.metadata.scale_factor);
  SnapshotSetMemberInfo member = payload.set_member;
  if (member.num_files_per_snapshot == 0U) {
    throw std::invalid_argument("snapshot writer: NumFilesPerSnapshot must be positive");
  }
  if (member.member_index >= member.num_files_per_snapshot) {
    throw std::invalid_argument("snapshot writer: member index is outside snapshot set");
  }
  if (member.generation_id.empty()) {
    if (member.num_files_per_snapshot != 1U) {
      throw std::invalid_argument(
          "snapshot writer: multifile output requires one shared non-empty generation id");
    }
    std::string generation_material = output_path.filename().string();
    generation_material += "\n" + payload.provenance.config_hash_hex;
    generation_material += "\n" + std::to_string(state.metadata.scale_factor);
    const std::string generation_hash = core::internal::sha256Hex(generation_material);
    member.generation_id = "snapshot_" + generation_hash.substr(0U, 16U);
  }
  const std::array<std::uint64_t, 6> global_count_by_type =
      member.has_global_part_count ? member.global_part_count : count_by_type;
  if (member.num_files_per_snapshot > 1U && !member.has_global_part_count) {
    throw std::invalid_argument(
        "snapshot writer: multifile output requires explicit global particle counts");
  }
  if (payload.collective_single_file) {
    if (member.member_index != 0U || member.num_files_per_snapshot != 1U ||
        !member.has_global_part_count) {
      throw std::invalid_argument(
          "collective single-file snapshot requires member_index=0, NumFilesPerSnapshot=1, and global counts");
    }
    for (std::size_t type = 0; type < 6U; ++type) {
      if (payload.collective_file_row_offset[type] > global_count_by_type[type] ||
          count_by_type[type] > global_count_by_type[type] - payload.collective_file_row_offset[type]) {
        throw std::runtime_error("collective snapshot rank hyperslab exceeds global dataset extent");
      }
    }
  }

  std::unique_ptr<internal::TransactionalFileTarget> transaction;
  std::filesystem::path write_path = output_path;
  if (!policy.caller_managed_publication) {
    transaction = std::make_unique<internal::TransactionalFileTarget>(
        output_path, ".part",
        policy.durable_publication
            ? internal::FileDurability::kDurablePublication
            : internal::FileDurability::kAtomicVisibility);
    write_path = transaction->temporaryPath();
  }
  SnapshotWriteMetrics write_metrics;
  {
  Hdf5Handle file_access;
#if COSMOSIM_ENABLE_MPI && COSMOSIM_HDF5_PARALLEL
  if (payload.collective_single_file) {
    file_access = Hdf5Handle(H5Pcreate(H5P_FILE_ACCESS));
    if (!file_access.valid() ||
        H5Pset_fapl_mpio(file_access.get(), MPI_COMM_WORLD, MPI_INFO_NULL) < 0 ||
        H5Pset_all_coll_metadata_ops(file_access.get(), 1) < 0 ||
        H5Pset_coll_metadata_write(file_access.get(), 1) < 0) {
      throw std::runtime_error(
          "failed configuring Parallel-HDF5 MPI-IO/collective-metadata file access");
    }
  }
#endif
  Hdf5Handle file(H5Fcreate(
      write_path.string().c_str(), H5F_ACC_TRUNC,
      H5P_DEFAULT, file_access.valid() ? file_access.get() : H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error(
        "failed creating snapshot file: " + write_path.string());
  }

  const auto& schema = scienceSnapshotSchemaMap();
  const auto& shared_names = sharedIoContractNames();
  writeScalarStringAttribute(
      file.get(), std::string(shared_names.file_kind_attribute), shared_names.science_snapshot_file_kind);
  Hdf5Handle header_group(H5Gcreate2(file.get(), std::string(schema.header_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!header_group.valid()) {
    throw std::runtime_error("failed creating /Header group");
  }

  writeHeaderArrays(
      header_group.get(),
      payload.collective_single_file ? global_count_by_type : count_by_type,
      global_count_by_type, mass_table, config, state, member, conversion);
  writeScalarStringAttribute(header_group.get(), "CosmoSimSchemaName", schema.schema_name);
  writeScalarUint32Attribute(header_group.get(), "CosmoSimSchemaVersion", schema.schema_version);
  writeScalarStringAttribute(
      header_group.get(), "CHUISnapshotDialect",
      internal::snapshotDialectLabel(write_dialect));
  writeScalarStringAttribute(header_group.get(), "CHUIUnitLength", config.units.length_unit);
  writeScalarStringAttribute(header_group.get(), "CHUIUnitMass", config.units.mass_unit);
  writeScalarStringAttribute(header_group.get(), "CHUIUnitVelocity", config.units.velocity_unit);
  writeScalarStringAttribute(
      header_group.get(), "CHUIVelocityStorageConvention",
      write_dialect == SnapshotDialect::kArepoFormat3 ||
              write_dialect == SnapshotDialect::kGadget4Hdf5
          ? "peculiar_velocity_div_sqrt_a"
          : "chui_internal_peculiar_velocity");
  writeScalarStringAttribute(header_group.get(), "NamingRulesVersion", "1.0");
  writeScalarStringAttribute(header_group.get(), "FileNamingRulesVersion", "1.0");
  writeScalarUint32Attribute(header_group.get(), "CHUILocalIndexWidthBits", 32U);
  writeScalarStringAttribute(
      header_group.get(), "CHUILocalIndexPolicy",
      "uint32_dense_local_rows_per_snapshot_member;global_counts_are_uint64");
  writeScalarStringAttribute(header_group.get(), "CosmoSimBuild", core::buildProvenance());
  writeScalarStringAttribute(
      header_group.get(), "MetalSpeciesMode",
      metalSpeciesModeLabel(config.physics.metal_species_mode));
  writeScalarStringAttribute(
      header_group.get(), "MetalDiffusionModel",
      metalDiffusionModelLabel(config.physics.metal_diffusion_model));
  writeScalarStringAttribute(
      header_group.get(), "MetalDiffusionTimeIntegrator",
      metalDiffusionIntegratorLabel(config.physics.metal_diffusion_time_integrator));
  writeScalarDoubleAttribute(
      header_group.get(), "MetalDiffusionCoefficient",
      config.physics.metal_diffusion_coefficient);
  writeScalarStringAttribute(
      header_group.get(), "StellarEvolutionTablePath",
      config.physics.stellar_evolution_table_path.empty()
          ? "builtin_zero_yield"
          : config.physics.stellar_evolution_table_path);

  Hdf5Handle config_group(
      H5Gcreate2(file.get(), std::string(schema.config_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeScalarStringAttribute(config_group.get(), std::string(schema.config_normalized_attribute), payload.normalized_config_text);
  Hdf5Handle parameters_group(
      H5Gcreate2(file.get(), std::string(schema.parameters_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!parameters_group.valid()) {
    throw std::runtime_error("failed creating /Parameters group");
  }
  writeScalarStringAttribute(parameters_group.get(), "normalized_param_txt", payload.normalized_config_text);
  Hdf5Handle units_group(
      H5Gcreate2(file.get(), std::string(schema.units_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!units_group.valid()) {
    throw std::runtime_error("failed creating /Units group");
  }
  writeScalarStringAttribute(units_group.get(), "LengthUnit", config.units.length_unit);
  writeScalarStringAttribute(units_group.get(), "MassUnit", config.units.mass_unit);
  writeScalarStringAttribute(units_group.get(), "VelocityUnit", config.units.velocity_unit);
  writeScalarStringAttribute(units_group.get(), "CoordinateFrame",
      config.units.coordinate_frame == core::CoordinateFrame::kComoving ? "comoving" : "physical");
  writeScalarStringAttribute(units_group.get(), "ExternalDialect", internal::snapshotDialectLabel(write_dialect));

  Hdf5Handle provenance_group(
      H5Gcreate2(file.get(), std::string(schema.provenance_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeScalarStringAttribute(provenance_group.get(), "schema_version", payload.provenance.schema_version);
  writeScalarStringAttribute(provenance_group.get(), "config_schema_name", payload.provenance.config_schema_name);
  writeScalarStringAttribute(
      provenance_group.get(), "config_schema_version", payload.provenance.config_schema_version);
  writeScalarStringAttribute(provenance_group.get(), "git_sha", payload.provenance.git_sha.empty() ? payload.git_sha : payload.provenance.git_sha);
  writeScalarStringAttribute(provenance_group.get(), "compiler_id", payload.provenance.compiler_id);
  writeScalarStringAttribute(provenance_group.get(), "compiler_version", payload.provenance.compiler_version);
  writeScalarStringAttribute(provenance_group.get(), "build_preset", payload.provenance.build_preset);
  writeScalarStringAttribute(provenance_group.get(), "enabled_features", payload.provenance.enabled_features);
  writeScalarStringAttribute(provenance_group.get(), "config_hash_hex", payload.provenance.config_hash_hex);
  writeScalarStringAttribute(
      provenance_group.get(),
      "normalized_config_hash_hex",
      payload.provenance.normalized_config_hash_hex.empty() ? payload.provenance.config_hash_hex
                                                            : payload.provenance.normalized_config_hash_hex);
  writeScalarStringAttribute(provenance_group.get(), "raw_input_config", payload.provenance.raw_input_config);
  writeScalarStringAttribute(provenance_group.get(), "normalized_config", payload.provenance.normalized_config);
  writeScalarStringAttribute(
      provenance_group.get(), "derived_runtime_state", payload.provenance.derived_runtime_state);
  writeScalarStringAttribute(provenance_group.get(), "timestamp_utc", payload.provenance.timestamp_utc);
  writeScalarStringAttribute(provenance_group.get(), "hardware_summary", payload.provenance.hardware_summary);
  if (payload.provenance.schema_version == "provenance_v7") {
    writeScalarStringAttribute(provenance_group.get(), "integrity_digest_algorithm", payload.provenance.integrity_digest_algorithm);
    writeScalarStringAttribute(provenance_group.get(), "normalized_config_sha256_hex", payload.provenance.normalized_config_sha256_hex);
    writeScalarStringAttribute(provenance_group.get(), "compiler_flags", payload.provenance.compiler_flags);
    writeScalarStringAttribute(provenance_group.get(), "cpu_model", payload.provenance.cpu_model);
    writeScalarUint32Attribute(provenance_group.get(), "logical_thread_count", payload.provenance.logical_thread_count);
    writeScalarUint32Attribute(provenance_group.get(), "physical_core_count", payload.provenance.physical_core_count);
    writeScalarUint64Attribute(provenance_group.get(), "system_ram_bytes", payload.provenance.system_ram_bytes);
    writeScalarStringAttribute(provenance_group.get(), "host_name", payload.provenance.host_name);
    writeScalarStringAttribute(provenance_group.get(), "gpu_summary", payload.provenance.gpu_summary);
    writeScalarStringAttribute(provenance_group.get(), "cuda_runtime_version", payload.provenance.cuda_runtime_version);
    writeScalarStringAttribute(provenance_group.get(), "cuda_driver_version", payload.provenance.cuda_driver_version);
    writeScalarStringAttribute(provenance_group.get(), "mpi_summary", payload.provenance.mpi_summary);
    writeScalarUint32Attribute(provenance_group.get(), "mpi_world_size", cosmosim::core::checkedIntegralNarrow<std::uint32_t>(payload.provenance.mpi_world_size, "snapshot provenance mpi_world_size"));
    writeScalarUint32Attribute(provenance_group.get(), "mpi_node_local_rank", cosmosim::core::checkedIntegralNarrow<std::uint32_t>(payload.provenance.mpi_node_local_rank, "snapshot provenance mpi_node_local_rank"));
    writeScalarStringAttribute(provenance_group.get(), "deterministic_mode", payload.provenance.deterministic_mode);
  }
  writeScalarUint32Attribute(
      provenance_group.get(),
      "gravity_treepm_pm_grid",
      static_cast<std::uint32_t>(std::max(payload.provenance.gravity_treepm_pm_grid, 0)));
  writeScalarUint32Attribute(
      provenance_group.get(),
      "gravity_treepm_pm_grid_nx",
      static_cast<std::uint32_t>(std::max(payload.provenance.gravity_treepm_pm_grid_nx, 0)));
  writeScalarUint32Attribute(
      provenance_group.get(),
      "gravity_treepm_pm_grid_ny",
      static_cast<std::uint32_t>(std::max(payload.provenance.gravity_treepm_pm_grid_ny, 0)));
  writeScalarUint32Attribute(
      provenance_group.get(),
      "gravity_treepm_pm_grid_nz",
      static_cast<std::uint32_t>(std::max(payload.provenance.gravity_treepm_pm_grid_nz, 0)));
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_treepm_assignment_scheme",
      payload.provenance.gravity_treepm_assignment_scheme);
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_treepm_window_deconvolution",
      payload.provenance.gravity_treepm_window_deconvolution ? "true" : "false");
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_asmth_cells",
      payload.provenance.gravity_treepm_asmth_cells);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_rcut_cells",
      payload.provenance.gravity_treepm_rcut_cells);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_mesh_spacing_mpc_comoving",
      payload.provenance.gravity_treepm_mesh_spacing_mpc_comoving);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_mesh_spacing_x_mpc_comoving",
      payload.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_mesh_spacing_y_mpc_comoving",
      payload.provenance.gravity_treepm_mesh_spacing_y_mpc_comoving);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_mesh_spacing_z_mpc_comoving",
      payload.provenance.gravity_treepm_mesh_spacing_z_mpc_comoving);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_split_scale_mpc_comoving",
      payload.provenance.gravity_treepm_split_scale_mpc_comoving);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_cutoff_radius_mpc_comoving",
      payload.provenance.gravity_treepm_cutoff_radius_mpc_comoving);
  writeScalarUint32Attribute(
      provenance_group.get(),
      "gravity_treepm_update_cadence_steps",
      static_cast<std::uint32_t>(std::max(payload.provenance.gravity_treepm_update_cadence_steps, 0)));
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_treepm_tree_opening_criterion",
      payload.provenance.gravity_treepm_tree_opening_criterion);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_tree_opening_theta",
      payload.provenance.gravity_treepm_tree_opening_theta);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_tree_relative_force_tolerance",
      payload.provenance.gravity_treepm_tree_relative_force_tolerance);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_tree_relative_force_acceleration_floor",
      payload.provenance.gravity_treepm_tree_relative_force_acceleration_floor);
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_softening_policy",
      payload.provenance.gravity_softening_policy);
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_softening_kernel",
      payload.provenance.gravity_softening_kernel);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_softening_epsilon_kpc_comoving",
      payload.provenance.gravity_softening_epsilon_kpc_comoving);
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_pm_fft_backend",
      payload.provenance.gravity_pm_fft_backend);
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_pm_backend_capability",
      payload.provenance.gravity_pm_backend_capability);
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_acceptance_profile_id",
      payload.provenance.gravity_acceptance_profile_id);

  for (std::size_t type_index = 0; type_index < schema.part_type_group.size(); ++type_index) {
    const std::size_t row_count = static_cast<std::size_t>(count_by_type[type_index]);
    const std::size_t file_row_count = static_cast<std::size_t>(
        payload.collective_single_file ? global_count_by_type[type_index] : count_by_type[type_index]);
    const std::size_t file_row_base = static_cast<std::size_t>(
        payload.collective_single_file ? payload.collective_file_row_offset[type_index] : 0U);
    if (file_row_count == 0U) {
      continue;
    }

    Hdf5Handle type_group(H5Gcreate2(
        file.get(), std::string(schema.part_type_group[type_index]).c_str(),
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    if (!type_group.valid()) {
      throw std::runtime_error("failed creating part type group");
    }

    auto coords_ds = createStreamingDataset(
        type_group.get(), schema.coordinates.canonical_name, H5T_IEEE_F64LE,
        file_row_count, 3, policy);
    auto velocities_ds = createStreamingDataset(
        type_group.get(), schema.velocities.canonical_name, H5T_IEEE_F64LE,
        file_row_count, 3, policy);
    auto masses_ds = createStreamingDataset(
        type_group.get(), schema.masses.canonical_name, H5T_IEEE_F64LE,
        file_row_count, 1, policy);
    auto ids_ds = createStreamingDataset(
        type_group.get(), schema.particle_ids.canonical_name, H5T_STD_U64LE,
        file_row_count, 1, policy);

    std::optional<StreamingDataset> softening_ds;
    std::optional<StreamingDataset> softening_override_ds;
    const bool has_particle_softening = type_index != 0U &&
        (payload.collective_single_file
             ? payload.collective_write_particle_softening
             : !state.particle_sidecar.gravity_softening_comoving.empty());
    if (has_particle_softening) {
      softening_ds.emplace(createStreamingDataset(
          type_group.get(), "GravitySofteningComoving", H5T_IEEE_F64LE,
          file_row_count, 1, policy));
      if (payload.collective_single_file
              ? payload.collective_write_particle_softening_override
              : !state.particle_sidecar.has_gravity_softening_override.empty()) {
        softening_override_ds.emplace(createStreamingDataset(
            type_group.get(), "GravitySofteningOverrideMask", H5T_STD_U8LE,
            file_row_count, 1, policy));
      }
    }

    std::unordered_map<std::string, StreamingDataset> extra;
    auto add_double = [&](std::string name) {
      extra.emplace(name, createStreamingDataset(
          type_group.get(), name, H5T_IEEE_F64LE, file_row_count, 1, policy));
    };
    auto add_u64 = [&](std::string name) {
      extra.emplace(name, createStreamingDataset(
          type_group.get(), name, H5T_STD_U64LE, file_row_count, 1, policy));
    };
    auto add_u32 = [&](std::string name) {
      extra.emplace(name, createStreamingDataset(
          type_group.get(), name, H5T_STD_U32LE, file_row_count, 1, policy));
    };
    auto add_u8 = [&](std::string name) {
      extra.emplace(name, createStreamingDataset(
          type_group.get(), name, H5T_STD_U8LE, file_row_count, 1, policy));
    };

    if (type_index == 0U) {
      writeScalarUint32Attribute(type_group.get(), "GasIdentitySchemaVersion", 1U);
      writeScalarStringAttribute(type_group.get(), "GasParticleIdsSemantics", "stable_gas_cell_id");
      writeScalarStringAttribute(
          type_group.get(), "GasParentIdentitySemantics",
          "ParentParticleIDs valid where HasParentParticle=1");
    }
    // The required native field list is shared with the independent validator,
    // so schema additions cannot silently update only one side of the contract.
    internal::forEachRequiredChuiSnapshotField(type_index, [&](const auto& field) {
      switch (field.storage) {
        case internal::SnapshotFieldStorage::kFloat64: add_double(field.name); break;
        case internal::SnapshotFieldStorage::kUint64: add_u64(field.name); break;
        case internal::SnapshotFieldStorage::kUint32: add_u32(field.name); break;
        case internal::SnapshotFieldStorage::kUint8: add_u8(field.name); break;
      }
    });
    if (type_index == 0U && policy.write_optional_pressure) {
      add_double("Pressure");
    }
    // CHUI Part 3 Checkpoint G: science-state export of the composition
    // sidecar (schema CHUI_CompositionSchemaVersion = 3, see
    // docs/output_schema.md). Present iff the run carries a composition
    // sidecar (identical decision on every rank: the layout is part of the
    // configuration); absent datasets mean "physics disabled", never zero.
    // Values are mass fractions of the composition gas mass (elements + dust)
    // so sum(elements) + sum(dust) == 1 per cell and they are unit-free.
    std::optional<StreamingDataset> composition_elements_ds;
    std::optional<StreamingDataset> composition_species_ds;
    std::optional<StreamingDataset> composition_dust_ds;
    const bool composition_out = type_index == 0U && state.gas_composition.enabled();
    const core::GasCompositionLayoutDescriptor composition_descriptor =
        composition_out ? state.gas_composition.descriptor() : core::GasCompositionLayoutDescriptor{};
    const bool composition_sgs = composition_out &&
        composition_descriptor.sgs_energy_lane != core::GasCompositionLayoutDescriptor::k_no_lane;
    if (composition_out) {
      const physics::CompositionLayout layout = physics::CompositionLayout::fromDescriptor(composition_descriptor);
      const auto names = [&](std::uint32_t begin, std::uint32_t count) {
        std::string joined;
        for (std::uint32_t l = begin; l < begin + count; ++l) {
          if (!joined.empty()) joined += ",";
          joined += layout.lanes().at(l).name;
        }
        return joined;
      };
      writeScalarUint32Attribute(type_group.get(), "CHUI_CompositionSchemaVersion",
                                 physics::k_gas_composition_schema);
      writeScalarStringAttribute(type_group.get(), "CHUI_CompositionLayoutName", composition_descriptor.layout_name);
      writeScalarStringAttribute(type_group.get(), "CHUI_CompositionFractionBasis",
                                 "lane_mass / (sum(element lanes) + sum(dust lanes))");
      writeScalarStringAttribute(
          type_group.get(), "CHUI_DetailedChemistryStorage",
          state.gas_composition.sparseStorageEnabled()
              ? "adaptive_sparse_by_gas_cell_id"
              : "dense_reference");
      writeScalarStringAttribute(type_group.get(), "CHUI_ElementNames",
                                 names(composition_descriptor.element_lane_begin, composition_descriptor.element_lane_count));
      composition_elements_ds.emplace(createStreamingDataset(
          type_group.get(), "CHUI_ElementMassFractions", H5T_IEEE_F64LE, file_row_count,
          composition_descriptor.element_lane_count, policy));
      if (composition_descriptor.species_lane_count != 0U) {
        writeScalarStringAttribute(type_group.get(), "CHUI_SpeciesNames",
                                   names(composition_descriptor.species_lane_begin, composition_descriptor.species_lane_count));
        composition_species_ds.emplace(createStreamingDataset(
            type_group.get(), "CHUI_SpeciesMassFractions", H5T_IEEE_F64LE, file_row_count,
            composition_descriptor.species_lane_count, policy));
      }
      if (composition_descriptor.dust_lane_count != 0U) {
        writeScalarStringAttribute(type_group.get(), "CHUI_DustNames",
                                   names(composition_descriptor.dust_lane_begin, composition_descriptor.dust_lane_count));
        composition_dust_ds.emplace(createStreamingDataset(
            type_group.get(), "CHUI_DustMassFractions", H5T_IEEE_F64LE, file_row_count,
            composition_descriptor.dust_lane_count, policy));
      }
      if (composition_sgs) add_double("CHUI_SgsSpecificEnergy");
      add_double("CHUI_MeanMolecularWeight");
      add_u8("CHUI_ChemistryModelTag");
      add_u8("CHUI_DetailedChemistryStateMaterialized");
      add_u8("CHUI_IsmRegimeTag");
      add_double("CHUI_ChemistryModelDwellYears");
      add_double("CHUI_IsmRegimeDwellYears");
      add_double("CHUI_ChemistryLastTemperatureK");
      add_double("CHUI_ChemistryLastCoolingTimeSeconds");
      add_double("CHUI_ChemistryLastErrorEstimate");
      add_double("CHUI_ChemistryLastElementFractionResidual");
      add_double("CHUI_ChemistryLastChargeResidual");
      add_double("CHUI_ChemistryCumulativeRadiativeLossErg");
      add_double("CHUI_ChemistryCumulativeHeatingGainErg");
    }
    // Stellar population science state (PartType4), present iff a Part 2
    // stellar population model is configured (rank-uniform decision).
    std::optional<StreamingDataset> star_channel_mass_ds;
    std::optional<StreamingDataset> star_channel_energy_ds;
    std::optional<StreamingDataset> star_birth_elements_ds;
    const bool population_out = type_index == 4U &&
        config.physics.stellar_population_model != core::StellarPopulationModelConfig::kLegacyTable;
    std::optional<physics::StellarPopulationStore> population_store;
    if (population_out) {
      writeScalarUint32Attribute(
          type_group.get(), "CHUI_StellarPopulationSchemaVersion",
          physics::k_stellar_population_schema);
      writeScalarStringAttribute(type_group.get(), "CHUI_StellarChannelOrder", "AGB,CC,Ia");
      star_channel_mass_ds.emplace(createStreamingDataset(
          type_group.get(), "CHUI_StellarReturnedMassByChannel", H5T_IEEE_F64LE, file_row_count, 3U, policy));
      star_channel_energy_ds.emplace(createStreamingDataset(
          type_group.get(), "CHUI_StellarFeedbackEnergyByChannelErg", H5T_IEEE_F64LE, file_row_count, 3U, policy));
      star_birth_elements_ds.emplace(createStreamingDataset(
          type_group.get(), "CHUI_BirthElementMassFractions", H5T_IEEE_F64LE, file_row_count,
          physics::k_chemical_element_count, policy));
      add_double("CHUI_BirthTimeCode");
      add_u32("CHUI_SnCcRealized");
      add_u32("CHUI_SnIaRealized");
      add_u64("CHUI_EventStreamKey");
      add_u8("CHUI_PopulationRecordInitialized");
      add_double("CHUI_ChronologyTimeCode");
      add_u64("CHUI_ChronologyStepIndex");
      add_u8("CHUI_ChronologyInitialized");
      population_store.emplace(state);
    }

    const std::size_t chunk_rows = std::min(policy.chunk_particle_count, row_count);
    std::vector<double> coords(chunk_rows * 3U);
    std::vector<double> velocities(chunk_rows * 3U);
    std::vector<double> masses(chunk_rows);
    std::vector<std::uint64_t> ids(chunk_rows);
    std::vector<double> softening(chunk_rows);
    std::vector<std::uint8_t> override_mask(chunk_rows);
    std::vector<double> d0(chunk_rows), d1(chunk_rows), d2(chunk_rows), d3(chunk_rows), d4(chunk_rows), d5(chunk_rows), d6(chunk_rows), d7(chunk_rows), d8(chunk_rows), d9(chunk_rows);
    std::vector<std::uint64_t> u640(chunk_rows), u641(chunk_rows), u642(chunk_rows);
    std::vector<std::uint32_t> u320(chunk_rows), u321(chunk_rows);
    std::vector<std::uint8_t> u80(chunk_rows), u81(chunk_rows);
    std::vector<std::uint32_t> particle_indices;
    particle_indices.reserve(chunk_rows);

    const auto write_double_extra = [&](std::string_view name, const double* values, std::size_t offset, std::size_t count) {
      writeStreamingChunk(extra.at(std::string(name)), H5T_NATIVE_DOUBLE, values, file_row_base + offset, count);
    };
    const auto write_u64_extra = [&](std::string_view name, const std::uint64_t* values, std::size_t offset, std::size_t count) {
      writeStreamingChunk(extra.at(std::string(name)), H5T_NATIVE_UINT64, values, file_row_base + offset, count);
    };
    const auto write_u32_extra = [&](std::string_view name, const std::uint32_t* values, std::size_t offset, std::size_t count) {
      writeStreamingChunk(extra.at(std::string(name)), H5T_NATIVE_UINT32, values, file_row_base + offset, count);
    };
    const auto write_u8_extra = [&](std::string_view name, const std::uint8_t* values, std::size_t offset, std::size_t count) {
      writeStreamingChunk(extra.at(std::string(name)), H5T_NATIVE_UINT8, values, file_row_base + offset, count);
    };

    auto write_base = [&](std::size_t output_offset, std::size_t count) {
      writeStreamingChunk(coords_ds, H5T_NATIVE_DOUBLE, coords.data(), file_row_base + output_offset, count);
      writeStreamingChunk(velocities_ds, H5T_NATIVE_DOUBLE, velocities.data(), file_row_base + output_offset, count);
      writeStreamingChunk(masses_ds, H5T_NATIVE_DOUBLE, masses.data(), file_row_base + output_offset, count);
      writeStreamingChunk(ids_ds, H5T_NATIVE_UINT64, ids.data(), file_row_base + output_offset, count);
      if (softening_ds.has_value()) {
        writeStreamingChunk(*softening_ds, H5T_NATIVE_DOUBLE, softening.data(), file_row_base + output_offset, count);
      }
      if (softening_override_ds.has_value()) {
        writeStreamingChunk(*softening_override_ds, H5T_NATIVE_UINT8, override_mask.data(), file_row_base + output_offset, count);
      }
    };

    if (type_index == 0U) {
      const core::UnitSystem units = core::makeUnitSystem(
          config.units.length_unit, config.units.mass_unit, config.units.velocity_unit);
      for (std::size_t offset = 0; offset < row_count; offset += chunk_rows) {
        const std::size_t count = std::min(chunk_rows, row_count - offset);
        std::fill_n(u80.begin(), count, 0U);
        std::fill_n(d3.begin(), count, 0.0);
        std::fill_n(d4.begin(), count, 0.0);
        std::fill_n(d5.begin(), count, 0.0);
        std::fill_n(d6.begin(), count, 0.0);
        for (std::size_t j = 0; j < count; ++j) {
          const std::size_t gas_row = offset + j;
          const auto* identity = state.gas_cell_identity.findByLocalRow(
              core::checkedLocalCellRow(gas_row, "snapshot gas row"));
          if (identity == nullptr || identity->gas_cell_id == 0U) {
            throw std::runtime_error("snapshot writer: gas-cell identity map does not cover dense rows");
          }
          coords[j * 3U + 0U] = conversion.positionToStored(state.cells.center_x_comoving[gas_row]);
          coords[j * 3U + 1U] = conversion.positionToStored(state.cells.center_y_comoving[gas_row]);
          coords[j * 3U + 2U] = conversion.positionToStored(state.cells.center_z_comoving[gas_row]);
          velocities[j * 3U + 0U] = conversion.velocityToStored(state.gas_cells.velocity_x_peculiar[gas_row]);
          velocities[j * 3U + 1U] = conversion.velocityToStored(state.gas_cells.velocity_y_peculiar[gas_row]);
          velocities[j * 3U + 2U] = conversion.velocityToStored(state.gas_cells.velocity_z_peculiar[gas_row]);
          masses[j] = conversion.massToStored(state.cells.mass_code[gas_row]);
          ids[j] = identity->gas_cell_id;
          d0[j] = conversion.internalEnergyToStored(state.gas_cells.internal_energy_code[gas_row]);
          d1[j] = conversion.densityComovingToStored(state.gas_cells.density_code[gas_row]);
          const double gas_mass = state.cells.mass_code[gas_row];
          d2[j] = gas_mass > 0.0
              ? std::clamp(state.gas_cells.metal_mass_code[gas_row] / gas_mass, 0.0, 1.0)
              : 0.0;
          d7[j] = conversion.pressureComovingToStored(
              state.gas_cells.pressure_code[gas_row]);
          d8[j] = state.gas_cells.temperature_code[gas_row];
          d9[j] = state.gas_cells.sound_speed_code[gas_row];
          u640[j] = identity->gas_cell_id;
          u641[j] = identity->parent_particle_id.value_or(0U);
          u81[j] = identity->parent_particle_id.has_value() ? 1U : 0U;
          u642[j] = identity->owning_patch_id;

          double sfr_code = 0.0;
          if (effective_eos_table.has_value()) {
            const double a = std::max(state.metadata.scale_factor, 1.0e-12);
            const double density_phys = config.units.coordinate_frame == core::CoordinateFrame::kComoving
                ? state.gas_cells.density_code[gas_row] / (a * a * a)
                : state.gas_cells.density_code[gas_row];
            const auto equilibrium = effective_eos_table->lookup(density_phys);
            if (equilibrium.above_threshold && equilibrium.valid &&
                state.gas_cells.internal_energy_code[gas_row] <=
                    equilibrium.entry.specific_internal_energy_eff_code *
                    (1.0 + config.physics.sf_effective_hot_excess_tolerance)) {
              const double long_lived_factor =
                  config.physics.sf_effective_birth_mass_convention ==
                          core::EffectiveIsmBirthMassConvention::kLongLivedMass
                      ? (1.0 - config.physics.sf_effective_massive_star_fraction)
                      : 1.0;
              sfr_code = gas_mass * long_lived_factor * equilibrium.entry.cold_mass_fraction /
                  std::max(equilibrium.entry.star_formation_timescale_code, 1.0e-30);
              d4[j] = equilibrium.entry.cold_mass_fraction;
              d5[j] = conversion.pressureComovingToStored(
                  equilibrium.entry.pressure_phys_code *
                  (config.units.coordinate_frame == core::CoordinateFrame::kComoving ? a * a * a : 1.0));
              d6[j] = conversion.internalEnergyToStored(
                  equilibrium.entry.specific_internal_energy_eff_code);
              u80[j] = 1U;
            }
          } else if (
              config.physics.enable_star_formation &&
              config.physics.star_formation_model == core::StarFormationModelKind::kLegacySchmidtThreshold &&
              state.gas_cells.density_code[gas_row] >= config.physics.sf_density_threshold_code &&
              state.gas_cells.temperature_code[gas_row] <= config.physics.sf_temperature_threshold_k) {
            const double g_code = core::newtonGravitationalConstantCode(units);
            const double t_ff = std::sqrt(
                3.0 * 3.14159265358979323846 /
                (32.0 * g_code * std::max(state.gas_cells.density_code[gas_row], 1.0e-30)));
            sfr_code = config.physics.sf_epsilon_ff * gas_mass / std::max(t_ff, 1.0e-30);
          }
          d3[j] = conversion.starFormationRateCodeToStored(sfr_code);
        }
        write_base(offset, count);
        write_double_extra("InternalEnergy", d0.data(), offset, count);
        write_double_extra("Density", d1.data(), offset, count);
        write_double_extra("Metallicity", d2.data(), offset, count);
        if (policy.write_optional_pressure) {
          write_double_extra("Pressure", d7.data(), offset, count);
        }
        write_double_extra("CHUI_TemperatureCode", d8.data(), offset, count);
        write_double_extra("CHUI_SoundSpeedCode", d9.data(), offset, count);
        write_double_extra("StarFormationRate", d3.data(), offset, count);
        write_double_extra("ColdCloudMassFraction", d4.data(), offset, count);
        write_double_extra("EffectivePressure", d5.data(), offset, count);
        write_double_extra("EffectiveInternalEnergy", d6.data(), offset, count);
        write_u8_extra("IsOnEffectiveEos", u80.data(), offset, count);
        write_u64_extra("GasCellIDs", u640.data(), offset, count);
        write_u64_extra("ParentParticleIDs", u641.data(), offset, count);
        write_u8_extra("HasParentParticle", u81.data(), offset, count);
        write_u64_extra("OwningPatchIDs", u642.data(), offset, count);
        if (composition_out) {
          const auto& sidecar = state.gas_composition;
          const auto& cd = composition_descriptor;
          std::vector<double> lanes_buffer;
          const auto write_lanes = [&](const StreamingDataset& dataset, std::uint32_t begin, std::uint32_t width) {
            lanes_buffer.assign(static_cast<std::size_t>(width) * count, 0.0);
            for (std::size_t j = 0; j < count; ++j) {
              const std::size_t row = offset + j;
              const double gas = physics::compositionGasMass(sidecar, row);
              for (std::uint32_t l = 0; l < width; ++l) {
                lanes_buffer[j * width + l] = gas > 0.0 ? sidecar.value(begin + l, row) / gas : 0.0;
              }
            }
            writeStreamingChunk(dataset, H5T_NATIVE_DOUBLE, lanes_buffer.data(), file_row_base + offset, count);
          };
          write_lanes(*composition_elements_ds, cd.element_lane_begin, cd.element_lane_count);
          if (composition_species_ds.has_value()) write_lanes(*composition_species_ds, cd.species_lane_begin, cd.species_lane_count);
          if (composition_dust_ds.has_value()) write_lanes(*composition_dust_ds, cd.dust_lane_begin, cd.dust_lane_count);
          for (std::size_t j = 0; j < count; ++j) {
            const std::size_t row = offset + j;
            const double gas = physics::compositionGasMass(sidecar, row);
            d0[j] = composition_sgs && gas > 0.0 ? sidecar.value(cd.sgs_energy_lane, row) / gas : 0.0;
            d1[j] = sidecar.mean_molecular_weight[row];
            d2[j] = sidecar.model_dwell_years[row];
            d3[j] = sidecar.ism_regime_dwell_years[row];
            u80[j] = sidecar.chemistry_model_tag[row];
            u81[j] = sidecar.ism_regime_tag[row];
          }
          if (composition_sgs) write_double_extra("CHUI_SgsSpecificEnergy", d0.data(), offset, count);
          write_double_extra("CHUI_MeanMolecularWeight", d1.data(), offset, count);
          write_double_extra("CHUI_ChemistryModelDwellYears", d2.data(), offset, count);
          write_double_extra("CHUI_IsmRegimeDwellYears", d3.data(), offset, count);
          write_u8_extra("CHUI_ChemistryModelTag", u80.data(), offset, count);
          write_u8_extra("CHUI_IsmRegimeTag", u81.data(), offset, count);
          for (std::size_t j = 0; j < count; ++j) {
            u80[j] = sidecar.hasExplicitSparseRow(offset + j) ? 1U : 0U;
          }
          write_u8_extra(
              "CHUI_DetailedChemistryStateMaterialized", u80.data(), offset, count);
          for (std::size_t j = 0; j < count; ++j) {
            const std::size_t row = offset + j;
            d0[j] = sidecar.chemistry_last_temperature_k[row];
            d1[j] = sidecar.chemistry_last_cooling_time_seconds[row];
            d2[j] = sidecar.chemistry_last_error_estimate[row];
            d3[j] = sidecar.chemistry_last_element_fraction_residual[row];
            d4[j] = sidecar.chemistry_last_charge_residual[row];
            d5[j] = sidecar.chemistry_cumulative_radiative_loss_erg[row];
            d6[j] = sidecar.chemistry_cumulative_heating_gain_erg[row];
          }
          write_double_extra("CHUI_ChemistryLastTemperatureK", d0.data(), offset, count);
          write_double_extra("CHUI_ChemistryLastCoolingTimeSeconds", d1.data(), offset, count);
          write_double_extra("CHUI_ChemistryLastErrorEstimate", d2.data(), offset, count);
          write_double_extra("CHUI_ChemistryLastElementFractionResidual", d3.data(), offset, count);
          write_double_extra("CHUI_ChemistryLastChargeResidual", d4.data(), offset, count);
          write_double_extra("CHUI_ChemistryCumulativeRadiativeLossErg", d5.data(), offset, count);
          write_double_extra("CHUI_ChemistryCumulativeHeatingGainErg", d6.data(), offset, count);
        }
        const std::uint64_t staging_bytes = static_cast<std::uint64_t>(
            count * (6U * sizeof(double) + sizeof(double) + sizeof(std::uint64_t) +
                     17U * sizeof(double) + 3U * sizeof(std::uint64_t) + 3U * sizeof(std::uint8_t)));
        const std::uint64_t logical_bytes = staging_bytes;
        write_metrics.observe(staging_bytes, logical_bytes);
      }
    } else {
      std::size_t output_offset = 0U;
      for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
        if (mapSpeciesTagToPartType(state.particleSpeciesTag(particle_index)) != type_index) {
          continue;
        }
        particle_indices.push_back(core::checkedLocalParticleRow(particle_index, "snapshot particle row"));
        if (particle_indices.size() < chunk_rows && output_offset + particle_indices.size() < row_count) {
          continue;
        }
        const std::size_t count = particle_indices.size();
        std::fill_n(override_mask.begin(), count, 0U);
        for (std::size_t j = 0; j < count; ++j) {
          const std::size_t idx = particle_indices[j];
          coords[j * 3U + 0U] = conversion.positionToStored(state.particles.position_x_comoving[idx]);
          coords[j * 3U + 1U] = conversion.positionToStored(state.particles.position_y_comoving[idx]);
          coords[j * 3U + 2U] = conversion.positionToStored(state.particles.position_z_comoving[idx]);
          velocities[j * 3U + 0U] = conversion.velocityToStored(state.particles.velocity_x_peculiar[idx]);
          velocities[j * 3U + 1U] = conversion.velocityToStored(state.particles.velocity_y_peculiar[idx]);
          velocities[j * 3U + 2U] = conversion.velocityToStored(state.particles.velocity_z_peculiar[idx]);
          masses[j] = conversion.massToStored(state.particles.mass_code[idx]);
          ids[j] = state.particle_sidecar.particle_id[idx];
          if (softening_ds.has_value()) {
            softening[j] = conversion.softeningComovingToStored(
                state.particle_sidecar.gravity_softening_comoving[idx]);
            if (softening_override_ds.has_value()) {
              override_mask[j] = state.particle_sidecar.has_gravity_softening_override[idx];
            }
          }

          if (type_index == 4U) {
            const std::size_t row = star_row_by_particle.rowFor(
                static_cast<std::uint32_t>(idx));
            d0[j] = state.star_particles.metallicity_mass_fraction[row];
            d1[j] = state.star_particles.formation_scale_factor[row];
            d2[j] = conversion.massToStored(state.star_particles.birth_mass_code[row]);
            u640[j] = state.star_particles.birth_key[row];
            u641[j] = state.star_particles.parent_gas_cell_id[row];
            u642[j] = state.star_particles.birth_tick[row];
            u320[j] = state.star_particles.birth_ordinal[row];
            d3[j] = state.star_particles.stellar_age_years_last[row];
            d4[j] = conversion.massToStored(state.star_particles.stellar_returned_mass_cumulative_code[row]);
            d5[j] = conversion.massToStored(state.star_particles.stellar_returned_metals_cumulative_code[row]);
            d6[j] = conversion.massToStored(state.star_particles.stellar_newly_synthesized_metals_cumulative_code[row]);
            d7[j] = state.star_particles.stellar_feedback_energy_cumulative_erg[row];
          } else if (type_index == 3U) {
            const std::size_t row = tracer_row_by_particle.rowFor(
                static_cast<std::uint32_t>(idx));
            u640[j] = state.tracers.parent_particle_id[row];
            u641[j] = state.tracers.injection_step[row];
            u320[j] = state.tracers.host_cell_index[row];
            d0[j] = state.tracers.mass_fraction_of_host[row];
            d1[j] = conversion.massToStored(state.tracers.last_host_mass_code[row]);
            d2[j] = conversion.massToStored(state.tracers.cumulative_exchanged_mass_code[row]);
          } else if (type_index == 5U) {
            const std::size_t row = black_hole_row_by_particle.rowFor(
                static_cast<std::uint32_t>(idx));
            d0[j] = conversion.massToStored(state.black_holes.subgrid_mass_code[row]);
            d1[j] = conversion.starFormationRateCodeToStored(state.black_holes.accretion_rate_code[row]);
            d2[j] = state.black_holes.feedback_energy_code[row];
            d3[j] = state.black_holes.eddington_ratio[row];
            d4[j] = conversion.massToStored(state.black_holes.cumulative_accreted_mass_code[row]);
            d5[j] = state.black_holes.cumulative_feedback_energy_code[row];
            d6[j] = state.black_holes.duty_cycle_active_time_code[row];
            d7[j] = state.black_holes.duty_cycle_total_time_code[row];
            u320[j] = state.black_holes.host_cell_index[row];
          }
        }
        write_base(output_offset, count);
        if (type_index == 4U) {
          write_double_extra("Metallicity", d0.data(), output_offset, count);
          write_double_extra("StellarFormationTime", d1.data(), output_offset, count);
          write_double_extra("BirthMass", d2.data(), output_offset, count);
          write_u64_extra("StarFormationBirthKey", u640.data(), output_offset, count);
          write_u64_extra("ParentGasCellID", u641.data(), output_offset, count);
          write_u64_extra("BirthIntegrationTick", u642.data(), output_offset, count);
          write_u32_extra("BirthOrdinal", u320.data(), output_offset, count);
          write_double_extra("CHUI_StellarAgeYearsLast", d3.data(), output_offset, count);
          write_double_extra("CHUI_StellarReturnedMassCumulative", d4.data(), output_offset, count);
          write_double_extra("CHUI_StellarReturnedMetalsCumulative", d5.data(), output_offset, count);
          write_double_extra("CHUI_StellarNewlySynthesizedMetalsCumulative", d6.data(), output_offset, count);
          write_double_extra("CHUI_StellarFeedbackEnergyCumulativeErg", d7.data(), output_offset, count);
          for (std::size_t j = 0; j < count; ++j) {
            const std::size_t row = star_row_by_particle.rowFor(particle_indices[j]);
            d0[j] = conversion.massToStored(state.star_particles.stellar_deposited_mass_cumulative_code[row]);
            d1[j] = conversion.massToStored(state.star_particles.stellar_deposited_metals_cumulative_code[row]);
            d2[j] = state.star_particles.stellar_deposited_feedback_energy_cumulative_erg[row];
          }
          write_double_extra("CHUI_StellarDepositedMassCumulative", d0.data(), output_offset, count);
          write_double_extra("CHUI_StellarDepositedMetalsCumulative", d1.data(), output_offset, count);
          write_double_extra("CHUI_StellarDepositedFeedbackEnergyCumulativeErg", d2.data(), output_offset, count);
          if (population_out) {
            std::vector<double> channel(count * 3U);
            std::vector<double> energy(count * 3U);
            std::vector<double> birth(count * physics::k_chemical_element_count, 0.0);
            for (std::size_t j = 0; j < count; ++j) {
              const std::size_t row = star_row_by_particle.rowFor(particle_indices[j]);
              for (std::size_t c = 0; c < 3U; ++c) {
                channel[j * 3U + c] = conversion.massToStored(
                    state.star_particles.stellar_returned_mass_channel_cumulative_code[c][row]);
                energy[j * 3U + c] = state.star_particles.stellar_feedback_energy_channel_cumulative_erg[c][row];
              }
              const std::uint64_t pid = state.particle_sidecar.particle_id[particle_indices[j]];
              const bool has = population_store->has(pid);
              const physics::StellarPopulationRecord record =
                  has ? population_store->get(pid) : physics::StellarPopulationRecord{};
              for (std::size_t e = 0; e < physics::k_chemical_element_count; ++e) {
                birth[j * physics::k_chemical_element_count + e] = record.birth_element_mass_fraction[e];
              }
              d0[j] = record.birth_time_code;
              u320[j] = record.cc_realized;
              u321[j] = record.ia_realized;
              u640[j] = record.event_stream_key;
              u80[j] = has && record.initialized != 0U ? 1U : 0U;
              d1[j] = record.chronology_time_code;
              u641[j] = record.chronology_step_index;
              u81[j] = has && record.chronology_initialized != 0U ? 1U : 0U;
            }
            writeStreamingChunk(*star_channel_mass_ds, H5T_NATIVE_DOUBLE, channel.data(), file_row_base + output_offset, count);
            writeStreamingChunk(*star_channel_energy_ds, H5T_NATIVE_DOUBLE, energy.data(), file_row_base + output_offset, count);
            writeStreamingChunk(*star_birth_elements_ds, H5T_NATIVE_DOUBLE, birth.data(), file_row_base + output_offset, count);
            write_double_extra("CHUI_BirthTimeCode", d0.data(), output_offset, count);
            write_u32_extra("CHUI_SnCcRealized", u320.data(), output_offset, count);
            write_u32_extra("CHUI_SnIaRealized", u321.data(), output_offset, count);
            write_u64_extra("CHUI_EventStreamKey", u640.data(), output_offset, count);
            write_u8_extra("CHUI_PopulationRecordInitialized", u80.data(), output_offset, count);
            write_double_extra("CHUI_ChronologyTimeCode", d1.data(), output_offset, count);
            write_u64_extra("CHUI_ChronologyStepIndex", u641.data(), output_offset, count);
            write_u8_extra("CHUI_ChronologyInitialized", u81.data(), output_offset, count);
          }
        } else if (type_index == 3U) {
          write_u64_extra("TracerParentParticleID", u640.data(), output_offset, count);
          write_u64_extra("TracerInjectionStep", u641.data(), output_offset, count);
          write_u32_extra("TracerHostCellIndex", u320.data(), output_offset, count);
          write_double_extra("TracerMassFractionOfHost", d0.data(), output_offset, count);
          write_double_extra("TracerLastHostMassCode", d1.data(), output_offset, count);
          write_double_extra("TracerCumulativeExchangedMassCode", d2.data(), output_offset, count);
        } else if (type_index == 5U) {
          write_double_extra("CHUI_BHSubgridMass", d0.data(), output_offset, count);
          write_double_extra("CHUI_BHAccretionRateMsunPerYr", d1.data(), output_offset, count);
          write_double_extra("CHUI_BHFeedbackEnergyCode", d2.data(), output_offset, count);
          write_double_extra("CHUI_BHEddingtonRatio", d3.data(), output_offset, count);
          write_double_extra("CHUI_BHCumulativeAccretedMass", d4.data(), output_offset, count);
          write_double_extra("CHUI_BHCumulativeFeedbackEnergyCode", d5.data(), output_offset, count);
          write_double_extra("CHUI_BHDutyCycleActiveTimeCode", d6.data(), output_offset, count);
          write_double_extra("CHUI_BHDutyCycleTotalTimeCode", d7.data(), output_offset, count);
          write_u32_extra("CHUI_BHHostCellIndex", u320.data(), output_offset, count);
        }
        const std::uint64_t staging_bytes = static_cast<std::uint64_t>(
            count * (7U * sizeof(double) + sizeof(std::uint64_t) +
                     8U * sizeof(double) + 3U * sizeof(std::uint64_t) +
                     2U * sizeof(std::uint32_t) + 2U * sizeof(std::uint8_t)));
        write_metrics.observe(staging_bytes, staging_bytes);
        output_offset += count;
        particle_indices.clear();
      }
      if (!particle_indices.empty()) {
        throw std::logic_error("snapshot writer internal chunk flush invariant failed");
      }
      if (output_offset != row_count) {
        throw std::logic_error("snapshot writer emitted unexpected PartType row count");
      }
    }

    if (policy.write_particle_type_alias_groups) {
      const std::string alias_group_path = toTypeAliasPath(type_index);
      if (H5Lcreate_hard(
              file.get(), std::string(schema.part_type_group[type_index]).c_str(),
              file.get(), alias_group_path.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0) {
        throw std::runtime_error("failed to create ParticleType compatibility alias");
      }
    }
  }

  writeScalarUint64Attribute(header_group.get(), "CHUIPeakStagingBytes", write_metrics.peak_staging_bytes);
  writeScalarUint64Attribute(header_group.get(), "CHUILogicalBytesWritten", write_metrics.logical_bytes_written);
  writeScalarUint64Attribute(header_group.get(), "CHUIChunkWriteCount", write_metrics.chunk_write_count);
  if (H5Fflush(file.get(), H5F_SCOPE_GLOBAL) < 0) {
    throw std::runtime_error("failed to flush snapshot file");
  }
  }  // close all HDF5 identifiers before filesystem publication
  if (transaction) {
    transaction->publish();
    internal::writeSnapshotMemberIntegritySidecar(
        output_path, member, policy.durable_publication);
  }
#endif
}

void verifySingleFileScienceSnapshotPartitionHdf5(
    const std::filesystem::path& input_path,
    const SnapshotWritePayload& payload,
    const SnapshotIoPolicy& policy) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(input_path);
  static_cast<void>(payload);
  static_cast<void>(policy);
  throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: snapshot verifier unavailable");
#else
  if (!payload.collective_single_file || payload.state == nullptr || payload.config == nullptr) {
    throw std::invalid_argument("single-file partition verification requires a collective payload with state/config");
  }
  const core::SimulationState& state = *payload.state;
  const core::SimulationConfig& config = *payload.config;
  const SnapshotDialect dialect = internal::resolveSnapshotWriteDialect(policy.dialect, config);
  const internal::SnapshotConversionContext conversion = internal::makeSnapshotConversionContext(
      dialect, config, state.metadata.scale_factor);
  Hdf5Handle file(H5Fopen(input_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed opening collective science snapshot for partition readback: " + input_path.string());
  }
  const auto& schema = scienceSnapshotSchemaMap();
  constexpr std::size_t k_verify_chunk = 4096U;
  const auto require_double_equal = [](double actual, double expected, std::string_view what) {
    if (actual != expected) {
      throw std::runtime_error("collective science snapshot partition readback mismatch in " + std::string(what));
    }
  };

  // Gas rows preserve dense local CellSoa ordering within each rank's global hyperslab.
  if (state.cells.size() != 0U) {
    Hdf5Handle group(H5Gopen2(file.get(), std::string(schema.part_type_group[0]).c_str(), H5P_DEFAULT));
    if (!group.valid()) throw std::runtime_error("collective science snapshot missing /PartType0");
    const std::size_t base = static_cast<std::size_t>(payload.collective_file_row_offset[0]);
    for (std::size_t offset = 0; offset < state.cells.size(); offset += k_verify_chunk) {
      const std::size_t count = std::min(k_verify_chunk, state.cells.size() - offset);
      std::vector<double> coords, velocities, masses;
      std::vector<std::uint64_t> ids;
      readDatasetChunk2d(group.get(), std::string(schema.coordinates.canonical_name), base + offset, count, coords);
      readDatasetChunk2d(group.get(), std::string(schema.velocities.canonical_name), base + offset, count, velocities);
      readDatasetChunk1d(group.get(), std::string(schema.masses.canonical_name), base + offset, count, masses);
      readDatasetChunkIds(group.get(), std::string(schema.particle_ids.canonical_name), base + offset, count, ids);
      for (std::size_t j = 0; j < count; ++j) {
        const std::size_t row = offset + j;
        const auto* identity = state.gas_cell_identity.findByLocalRow(
            core::checkedLocalCellRow(row, "collective snapshot readback gas row"));
        if (identity == nullptr || ids[j] != identity->gas_cell_id) {
          throw std::runtime_error("collective science snapshot gas ParticleIDs readback mismatch");
        }
        require_double_equal(coords[3*j+0], conversion.positionToStored(state.cells.center_x_comoving[row]), "gas Coordinates");
        require_double_equal(coords[3*j+1], conversion.positionToStored(state.cells.center_y_comoving[row]), "gas Coordinates");
        require_double_equal(coords[3*j+2], conversion.positionToStored(state.cells.center_z_comoving[row]), "gas Coordinates");
        require_double_equal(velocities[3*j+0], conversion.velocityToStored(state.gas_cells.velocity_x_peculiar[row]), "gas Velocities");
        require_double_equal(velocities[3*j+1], conversion.velocityToStored(state.gas_cells.velocity_y_peculiar[row]), "gas Velocities");
        require_double_equal(velocities[3*j+2], conversion.velocityToStored(state.gas_cells.velocity_z_peculiar[row]), "gas Velocities");
        require_double_equal(masses[j], conversion.massToStored(state.cells.mass_code[row]), "gas Masses");
      }
    }
  }

  for (std::size_t type_index = 1U; type_index < 6U; ++type_index) {
    std::vector<std::uint32_t> rows;
    rows.reserve(std::min<std::size_t>(k_verify_chunk, state.particles.size()));
    std::size_t output_offset = 0U;
    Hdf5Handle group;
    const std::size_t base = static_cast<std::size_t>(payload.collective_file_row_offset[type_index]);
    auto flush = [&]() {
      if (rows.empty()) return;
      if (!group.valid()) {
        group = Hdf5Handle(H5Gopen2(file.get(), std::string(schema.part_type_group[type_index]).c_str(), H5P_DEFAULT));
        if (!group.valid()) throw std::runtime_error("collective science snapshot missing populated PartType group");
      }
      std::vector<double> coords, velocities, masses;
      std::vector<std::uint64_t> ids;
      readDatasetChunk2d(group.get(), std::string(schema.coordinates.canonical_name), base + output_offset, rows.size(), coords);
      readDatasetChunk2d(group.get(), std::string(schema.velocities.canonical_name), base + output_offset, rows.size(), velocities);
      readDatasetChunk1d(group.get(), std::string(schema.masses.canonical_name), base + output_offset, rows.size(), masses);
      readDatasetChunkIds(group.get(), std::string(schema.particle_ids.canonical_name), base + output_offset, rows.size(), ids);
      for (std::size_t j = 0; j < rows.size(); ++j) {
        const std::size_t row = rows[j];
        if (ids[j] != state.particle_sidecar.particle_id[row]) {
          throw std::runtime_error("collective science snapshot particle ID readback mismatch");
        }
        require_double_equal(coords[3*j+0], conversion.positionToStored(state.particles.position_x_comoving[row]), "particle Coordinates");
        require_double_equal(coords[3*j+1], conversion.positionToStored(state.particles.position_y_comoving[row]), "particle Coordinates");
        require_double_equal(coords[3*j+2], conversion.positionToStored(state.particles.position_z_comoving[row]), "particle Coordinates");
        require_double_equal(velocities[3*j+0], conversion.velocityToStored(state.particles.velocity_x_peculiar[row]), "particle Velocities");
        require_double_equal(velocities[3*j+1], conversion.velocityToStored(state.particles.velocity_y_peculiar[row]), "particle Velocities");
        require_double_equal(velocities[3*j+2], conversion.velocityToStored(state.particles.velocity_z_peculiar[row]), "particle Velocities");
        require_double_equal(masses[j], conversion.massToStored(state.particles.mass_code[row]), "particle Masses");
      }
      output_offset += rows.size();
      rows.clear();
    };
    for (std::size_t row = 0; row < state.particles.size(); ++row) {
      if (mapSpeciesTagToPartType(state.particleSpeciesTag(row)) != type_index) continue;
      rows.push_back(core::checkedLocalParticleRow(row, "collective snapshot readback particle row"));
      if (rows.size() == k_verify_chunk) flush();
    }
    flush();
  }
#endif
}

void publishSingleFileScienceSnapshot(
    const std::filesystem::path& temporary_path,
    const std::filesystem::path& final_path,
    const SnapshotSetMemberInfo& member,
    bool durable_publication) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(temporary_path);
  static_cast<void>(final_path);
  static_cast<void>(member);
  static_cast<void>(durable_publication);
  throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: snapshot publication unavailable");
#else
  if (member.member_index != 0U || member.num_files_per_snapshot != 1U ||
      !member.has_global_part_count || member.generation_id.empty()) {
    throw std::invalid_argument("single-file snapshot publication received invalid logical member metadata");
  }
  internal::atomicReplaceFile(
      temporary_path, final_path,
      durable_publication ? internal::FileDurability::kDurablePublication
                          : internal::FileDurability::kAtomicVisibility);
#endif
}

void writeGadgetArepoSnapshotHdf5(
    const std::filesystem::path& output_path,
    const SnapshotWritePayload& payload,
    const SnapshotIoPolicy& policy) {
  writeScienceSnapshotHdf5(output_path, payload, policy);
}

SnapshotReadResult readGadgetArepoSnapshotHdf5(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    const SnapshotReadOptions& options) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(input_path);
  static_cast<void>(config);
  static_cast<void>(options);
  throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: snapshot reader unavailable");
#else
  SnapshotReadResult result;
  SnapshotReadBudgetScope read_budget_scope(options.budget);
  Hdf5Handle file(H5Fopen(input_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed opening snapshot file: " + input_path.string());
  }

  const auto& schema = gadgetArepoSchemaMap();
  const auto& shared_names = sharedIoContractNames();
  if (hdf5AttributeExists(file.get(), std::string(shared_names.file_kind_attribute))) {
    readScalarStringAttribute(file.get(), std::string(shared_names.file_kind_attribute), result.report.file_kind);
    if (result.report.file_kind == shared_names.restart_checkpoint_file_kind) {
      throw std::runtime_error(
          "snapshot reader rejected HDF5 file-kind '" + result.report.file_kind +
          "'; expected science_snapshot or legacy GADGET/AREPO particle snapshot");
    }
  } else {
    result.report.file_kind = "legacy_or_external_snapshot";
  }
  result.report.restart_compatible = false;
  const bool chui_authored =
      result.report.file_kind == shared_names.science_snapshot_file_kind;

  Hdf5Handle header_group(H5Gopen2(file.get(), std::string(schema.header_group).c_str(), H5P_DEFAULT));
  if (!header_group.valid()) {
    throw std::runtime_error("snapshot missing " + std::string(schema.header_group));
  }

  readOptionalHeaderDouble(header_group.get(), "Time", 1.0, result.state.metadata.scale_factor);
  result.report.header_time = result.state.metadata.scale_factor;
  readOptionalHeaderDouble(
      header_group.get(),
      "Redshift",
      result.report.header_time > 0.0 ? 1.0 / result.report.header_time - 1.0 : 0.0,
      result.report.header_redshift);
  double scalar_box_size_stored = config.cosmology.box_size_mpc_comoving;
  readOptionalHeaderDouble(
      header_group.get(), "BoxSize", config.cosmology.box_size_mpc_comoving,
      scalar_box_size_stored);
  readOptionalHeaderDouble(
      header_group.get(), "Omega0", config.cosmology.omega_matter, result.report.header_omega_matter);
  readOptionalHeaderDouble(
      header_group.get(), "OmegaLambda", config.cosmology.omega_lambda, result.report.header_omega_lambda);
  readOptionalHeaderDouble(
      header_group.get(), "OmegaBaryon", config.cosmology.omega_baryon, result.report.header_omega_baryon);
  readOptionalHeaderDouble(
      header_group.get(), "HubbleParam", config.cosmology.hubble_param, result.report.header_hubble_param);
  result.report.schema_name = std::string(schema.schema_name);
  readScalarStringAttribute(header_group.get(), "CosmoSimSchemaName", result.report.schema_name);

  std::uint32_t schema_version = schema.schema_version;
  {
    Hdf5Handle attr(H5Aopen(header_group.get(), "CosmoSimSchemaVersion", H5P_DEFAULT));
    if (attr.valid()) {
      H5Aread(attr.get(), H5T_NATIVE_UINT32, &schema_version);
    }
  }
  result.report.schema_version = schema_version;

  SnapshotDialect read_dialect = options.dialect;
  std::string stored_dialect_label;
  readScalarStringAttribute(header_group.get(), "CHUISnapshotDialect", stored_dialect_label);
  if (!stored_dialect_label.empty()) {
    if (stored_dialect_label == "chui_native") read_dialect = SnapshotDialect::kChuiNative;
    else if (stored_dialect_label == "arepo_format3") read_dialect = SnapshotDialect::kArepoFormat3;
    else if (stored_dialect_label == "gadget4_hdf5") read_dialect = SnapshotDialect::kGadget4Hdf5;
    else throw std::runtime_error("snapshot reader: unknown CHUISnapshotDialect '" + stored_dialect_label + "'");
  } else if (read_dialect == SnapshotDialect::kAuto) {
    // Historical CHUI gadget_arepo_v<=5 files stored internal values directly.
    if (chui_authored || result.report.schema_name.rfind("gadget_arepo_v", 0) == 0) {
      read_dialect = SnapshotDialect::kChuiNative;
    } else {
      throw std::runtime_error(
          "external snapshot import requires an explicit SnapshotDialect; PartType names alone do not define units/velocity semantics");
    }
  }
  result.report.dialect = read_dialect;
  internal::SnapshotConversionContext conversion =
      internal::makeSnapshotConversionContext(read_dialect, config, result.state.metadata.scale_factor);
  conversion.hubble_param = result.report.header_hubble_param;

  double axis_x = 0.0;
  double axis_y = 0.0;
  double axis_z = 0.0;
  if (readScalarDoubleAttribute(header_group.get(), "CHUIBoxSizeX_MpcComoving", axis_x) &&
      readScalarDoubleAttribute(header_group.get(), "CHUIBoxSizeY_MpcComoving", axis_y) &&
      readScalarDoubleAttribute(header_group.get(), "CHUIBoxSizeZ_MpcComoving", axis_z)) {
    result.report.header_box_size_x = axis_x;
    result.report.header_box_size_y = axis_y;
    result.report.header_box_size_z = axis_z;
  } else {
    readOptionalHeaderDouble(header_group.get(), "CosmoSimBoxSizeX", conversion.boxSizeStoredToMpc(scalar_box_size_stored), result.report.header_box_size_x);
    readOptionalHeaderDouble(header_group.get(), "CosmoSimBoxSizeY", conversion.boxSizeStoredToMpc(scalar_box_size_stored), result.report.header_box_size_y);
    readOptionalHeaderDouble(header_group.get(), "CosmoSimBoxSizeZ", conversion.boxSizeStoredToMpc(scalar_box_size_stored), result.report.header_box_size_z);
  }
  // SnapshotReadReport preserves the public/header contract in comoving Mpc,
  // while decoded CHUI-native positions are expressed in the configured code
  // length unit. Bounds validation must therefore compare like units. Routing
  // the Mpc header values through the dialect-aware box conversion gives the
  // corresponding code-space box lengths without changing report/schema units.
  const double header_box_size_x_code = conversion.positionFromStored(
      conversion.boxSizeMpcToStored(result.report.header_box_size_x));
  const double header_box_size_y_code = conversion.positionFromStored(
      conversion.boxSizeMpcToStored(result.report.header_box_size_y));
  const double header_box_size_z_code = conversion.positionFromStored(
      conversion.boxSizeMpcToStored(result.report.header_box_size_z));

  std::array<std::uint64_t, 6> header_counts{};
  std::array<double, 6> header_mass_table{};
  readHeaderCounts(header_group.get(), header_counts);
  std::array<std::uint64_t, 6> global_header_counts{};
  readHeaderGlobalCounts(header_group.get(), header_counts, global_header_counts);
  readOptionalHeaderMassTable(header_group.get(), header_mass_table);
  result.report.local_part_count = header_counts;
  result.report.global_part_count = global_header_counts;
  static_cast<void>(readScalarUint32Attribute(
      header_group.get(), "NumFilesPerSnapshot", result.report.num_files_per_snapshot));
  static_cast<void>(readScalarUint32Attribute(
      header_group.get(), "CHUISnapshotMemberIndex", result.report.member_index));
  if (H5Aexists(header_group.get(), "CHUISnapshotGenerationID") > 0) {
    readScalarStringAttribute(
        header_group.get(), "CHUISnapshotGenerationID", result.report.generation_id);
  }
  if (result.report.num_files_per_snapshot == 0U) {
    throw std::runtime_error("snapshot reader: NumFilesPerSnapshot must be positive");
  }
  std::uint64_t total_count_u64 = 0U;
  for (std::uint64_t count : header_counts) {
    if (count > options.budget.max_particles - total_count_u64) {
      throw std::length_error("snapshot reader: particle-count read budget exceeded");
    }
    total_count_u64 += count;
  }
  if (header_counts[0] > options.budget.max_gas_cells) {
    throw std::length_error("snapshot reader: gas-cell read budget exceeded");
  }
  const std::uint64_t sidecar_rows = header_counts[0] + header_counts[3] +
      header_counts[4] + header_counts[5];
  if (sidecar_rows > options.budget.max_sidecar_rows) {
    throw std::length_error("snapshot reader: sidecar-row read budget exceeded");
  }
  const std::uint64_t base_bytes_per_particle =
      7U * sizeof(double) + sizeof(std::uint64_t) + 4U * sizeof(std::uint32_t);
  if (total_count_u64 != 0U &&
      total_count_u64 > options.budget.max_materialized_bytes / base_bytes_per_particle) {
    throw std::length_error("snapshot reader: materialized-byte read budget exceeded");
  }
  const std::size_t total_count = checkedSnapshotCountToSize(total_count_u64, "snapshot reader");
  result.state.resizeParticles(total_count);
  const core::ModePolicy mode_policy = core::buildModePolicy(config.mode);
  const bool periodic_external_positions = !chui_authored &&
      mode_policy.gravity_boundary == core::GravityBoundaryModel::kPeriodicPoisson;
  bool source_id_has_zero = false;
  bool source_id_has_uint64_max = false;
  bool source_id_affine_plus_one = false;
  std::uint64_t periodically_wrapped_components = 0U;
  std::vector<std::uint32_t> tracer_particle_index;
  std::vector<std::uint64_t> tracer_parent_particle_id;
  std::vector<std::uint64_t> tracer_injection_step;
  std::vector<std::uint32_t> tracer_host_cell_index;
  std::vector<double> tracer_mass_fraction_of_host;
  std::vector<double> tracer_last_host_mass_code;
  std::vector<double> tracer_cumulative_exchanged_mass_code;
  std::vector<std::uint32_t> gas_particle_index;
  std::vector<double> gas_internal_energy_code;
  std::vector<double> gas_density_code;
  std::vector<double> gas_pressure_code;
  std::vector<double> gas_temperature_code;
  std::vector<double> gas_sound_speed_code;
  std::vector<double> gas_metallicity_mass_fraction;
  std::vector<std::uint64_t> gas_cell_id;
  std::vector<std::uint64_t> gas_parent_particle_id;
  std::vector<std::uint8_t> gas_has_parent_particle;
  std::vector<std::uint64_t> gas_owning_patch_id;
  std::vector<std::uint32_t> star_particle_index;
  std::vector<double> star_metallicity_mass_fraction;
  std::vector<double> star_formation_scale_factor;
  std::vector<double> star_birth_mass_code;
  std::vector<std::uint64_t> star_birth_key;
  std::vector<std::uint64_t> star_parent_gas_cell_id;
  std::vector<std::uint64_t> star_birth_tick;
  std::vector<std::uint32_t> star_birth_ordinal;
  std::vector<double> star_age_years_last;
  std::vector<double> star_returned_mass_cumulative_code;
  std::vector<double> star_returned_metals_cumulative_code;
  std::vector<double> star_newly_synthesized_metals_cumulative_code;
  std::vector<double> star_feedback_energy_cumulative_erg;
  std::vector<double> star_deposited_mass_cumulative_code;
  std::vector<double> star_deposited_metals_cumulative_code;
  std::vector<double> star_deposited_feedback_energy_cumulative_erg;
  std::vector<std::uint32_t> black_hole_particle_index;
  std::vector<std::uint32_t> black_hole_host_cell_index;
  std::vector<double> black_hole_subgrid_mass_code;
  std::vector<double> black_hole_accretion_rate_code;
  std::vector<double> black_hole_feedback_energy_code;
  std::vector<double> black_hole_eddington_ratio;
  std::vector<double> black_hole_cumulative_accreted_mass_code;
  std::vector<double> black_hole_cumulative_feedback_energy_code;
  std::vector<double> black_hole_duty_cycle_active_time_code;
  std::vector<double> black_hole_duty_cycle_total_time_code;

  auto missing_double_field = [&](std::string_view field, std::size_t count, std::vector<double>& target, double documented_default) {
    if (options.missing_field_policy == SnapshotMissingFieldPolicy::kReject) {
      throw std::runtime_error("snapshot import: required scientific field is missing: " + std::string(field));
    }
    target.assign(count, documented_default);
    if (options.missing_field_policy == SnapshotMissingFieldPolicy::kMarkUnavailable) {
      result.report.unavailable_fields.push_back(std::string(field));
    } else {
      result.report.defaulted_fields.push_back(std::string(field) + "=documented_default");
    }
  };

  auto missing_u64_field = [&](std::string_view field, std::size_t count, std::vector<std::uint64_t>& target, std::uint64_t documented_default) {
    if (options.missing_field_policy == SnapshotMissingFieldPolicy::kReject) {
      throw std::runtime_error("snapshot import: required scientific field is missing: " + std::string(field));
    }
    target.assign(count, documented_default);
    if (options.missing_field_policy == SnapshotMissingFieldPolicy::kMarkUnavailable) {
      result.report.unavailable_fields.push_back(std::string(field));
    } else {
      result.report.defaulted_fields.push_back(std::string(field) + "=documented_default");
    }
  };

  auto missing_u32_field = [&](std::string_view field, std::size_t count, std::vector<std::uint32_t>& target, std::uint32_t documented_default) {
    if (options.missing_field_policy == SnapshotMissingFieldPolicy::kReject) {
      throw std::runtime_error("snapshot import: required scientific field is missing: " + std::string(field));
    }
    target.assign(count, documented_default);
    if (options.missing_field_policy == SnapshotMissingFieldPolicy::kMarkUnavailable) {
      result.report.unavailable_fields.push_back(std::string(field));
    } else {
      result.report.defaulted_fields.push_back(std::string(field) + "=documented_default");
    }
  };

  std::size_t global_offset = 0;
  for (std::size_t type_index = 0; type_index < header_counts.size(); ++type_index) {
    const std::size_t local_count = static_cast<std::size_t>(header_counts[type_index]);
    if (local_count == 0) {
      continue;
    }

    Hdf5Handle group(H5Gopen2(file.get(), std::string(schema.part_type_group[type_index]).c_str(), H5P_DEFAULT));
    if (!group.valid()) {
      const std::string alias_group = toTypeAliasPath(type_index);
      group = Hdf5Handle(H5Gopen2(file.get(), alias_group.c_str(), H5P_DEFAULT));
      if (!group.valid()) {
        throw std::runtime_error("missing particle group for type " + std::to_string(type_index));
      }
      result.report.present_aliases.push_back(std::string(schema.part_type_group[type_index]) + "=" + alias_group);
    }

    const std::string coordinates_name = pickAlias(
        group.get(),
        schema.coordinates,
        result.report,
        std::string(schema.part_type_group[type_index]) + "/Coordinates",
        true);
    observeDatasetStoragePolicy(group.get(), coordinates_name, result.report);
    const std::string velocities_name = pickAlias(
        group.get(),
        schema.velocities,
        result.report,
        std::string(schema.part_type_group[type_index]) + "/Velocities",
        options.require_velocities);
    const std::string ids_name = pickAlias(
        group.get(),
        schema.particle_ids,
        result.report,
        std::string(schema.part_type_group[type_index]) + "/ParticleIDs",
        options.require_ids);
    const std::string masses_name = pickAlias(
        group.get(),
        schema.masses,
        result.report,
        std::string(schema.part_type_group[type_index]) + "/Masses",
        false);

    std::vector<double> coords_chunk;
    std::vector<double> vel_chunk;
    std::vector<double> mass_chunk;
    std::vector<std::uint64_t> ids_chunk;
    std::vector<std::uint64_t> tracer_parent_chunk;
    std::vector<std::uint64_t> tracer_step_chunk;
    std::vector<std::uint32_t> tracer_host_chunk;
    std::vector<double> tracer_fraction_chunk;
    std::vector<double> tracer_last_host_mass_chunk;
    std::vector<double> tracer_exchange_chunk;
    std::vector<double> softening_chunk;
    std::vector<std::uint8_t> softening_override_mask_chunk;
    std::vector<double> gas_internal_energy_chunk;
    std::vector<double> gas_density_chunk;
    std::vector<double> gas_pressure_chunk;
    std::vector<double> gas_temperature_chunk;
    std::vector<double> gas_sound_speed_chunk;
    std::vector<double> gas_metallicity_chunk;
    std::vector<std::uint64_t> gas_cell_id_chunk;
    std::vector<std::uint64_t> gas_parent_particle_id_chunk;
    std::vector<std::uint8_t> gas_has_parent_particle_chunk;
    std::vector<std::uint64_t> gas_owning_patch_id_chunk;
    std::vector<double> star_metallicity_chunk;
    std::vector<double> star_formation_time_chunk;
    std::vector<double> star_birth_mass_chunk;
    std::vector<std::uint64_t> star_birth_key_chunk;
    std::vector<std::uint64_t> star_parent_gas_cell_id_chunk;
    std::vector<std::uint64_t> star_birth_tick_chunk;
    std::vector<std::uint32_t> star_birth_ordinal_chunk;
    std::vector<double> star_age_chunk;
    std::vector<double> star_returned_mass_chunk;
    std::vector<double> star_returned_metals_chunk;
    std::vector<double> star_new_metals_chunk;
    std::vector<double> star_feedback_energy_chunk;
    std::vector<double> star_deposited_mass_chunk;
    std::vector<double> star_deposited_metals_chunk;
    std::vector<double> star_deposited_energy_chunk;
    std::vector<double> bh_subgrid_mass_chunk;
    std::vector<double> bh_accretion_rate_chunk;
    std::vector<double> bh_feedback_energy_chunk;
    std::vector<double> bh_eddington_chunk;
    std::vector<double> bh_cumulative_accreted_chunk;
    std::vector<double> bh_cumulative_feedback_chunk;
    std::vector<double> bh_duty_active_chunk;
    std::vector<double> bh_duty_total_chunk;
    std::vector<std::uint32_t> bh_host_cell_chunk;

    readDatasetChunk2d(group.get(), coordinates_name, 0, local_count, coords_chunk);
    if (!velocities_name.empty()) {
      readDatasetChunk2d(group.get(), velocities_name, 0, local_count, vel_chunk);
    } else {
      missing_double_field(
          std::string(schema.part_type_group[type_index]) + "/Velocities",
          local_count * 3U, vel_chunk, 0.0);
    }

    if (!ids_name.empty()) {
      readDatasetChunkIds(group.get(), ids_name, 0, local_count, ids_chunk);
    } else {
      if (!options.allow_generated_particle_ids ||
          options.missing_field_policy !=
              SnapshotMissingFieldPolicy::kFillDocumentedDefault) {
        throw std::runtime_error(
            "snapshot import: ParticleIDs are missing; generated IDs require "
            "allow_generated_particle_ids=true and FillDocumentedDefault");
      }
      ids_chunk.resize(local_count);
      for (std::size_t i = 0; i < local_count; ++i) {
        ids_chunk[i] = static_cast<std::uint64_t>(global_offset + i + 1);
      }
      result.report.defaulted_fields.push_back(
          std::string(schema.part_type_group[type_index]) +
          "/ParticleIDs=generated_explicitly");
    }

    if (!masses_name.empty()) {
      readDatasetChunk1d(group.get(), masses_name, 0, local_count, mass_chunk);
    } else if (options.allow_mass_table_fallback) {
      const double constant_mass = header_mass_table[type_index];
      if (!(constant_mass > 0.0)) {
        throw std::runtime_error(
            "Masses missing and Header/MassTable fallback is unavailable for type " +
            std::to_string(type_index));
      }
      mass_chunk.assign(local_count, constant_mass);
      result.report.defaulted_fields.push_back(
          std::string(schema.part_type_group[type_index]) + "/Masses=MassTable");
    } else {
      throw std::runtime_error("Masses missing and fallback disabled");
    }
    if (hdf5PathExists(group.get(), "GravitySofteningComoving")) {
      readDatasetChunk1d(group.get(), "GravitySofteningComoving", 0, local_count, softening_chunk);
    } else {
      softening_chunk.clear();
    }
    if (hdf5PathExists(group.get(), "GravitySofteningOverrideMask")) {
      readDatasetChunk1dU8(group.get(), "GravitySofteningOverrideMask", 0, local_count, softening_override_mask_chunk);
    } else {
      softening_override_mask_chunk.clear();
    }
    if (type_index == 0U) {
      if (hdf5PathExists(group.get(), "InternalEnergy")) {
        readDatasetChunk1d(
            group.get(), "InternalEnergy", 0, local_count,
            gas_internal_energy_chunk);
      } else {
        missing_double_field(
            "/PartType0/InternalEnergy", local_count, gas_internal_energy_chunk, 0.0);
      }
      if (hdf5PathExists(group.get(), "Density")) {
        readDatasetChunk1d(
            group.get(), "Density", 0, local_count, gas_density_chunk);
      } else {
        missing_double_field(
            "/PartType0/Density", local_count, gas_density_chunk, 0.0);
      }
      if (hdf5PathExists(group.get(), "Pressure")) {
        readDatasetChunk1d(
            group.get(), "Pressure", 0, local_count, gas_pressure_chunk);
      } else {
        missing_double_field("/PartType0/Pressure", local_count, gas_pressure_chunk, 0.0);
      }
      if (hdf5PathExists(group.get(), "CHUI_TemperatureCode")) {
        readDatasetChunk1d(
            group.get(), "CHUI_TemperatureCode", 0, local_count,
            gas_temperature_chunk);
      } else {
        missing_double_field(
            "/PartType0/CHUI_TemperatureCode", local_count,
            gas_temperature_chunk, 0.0);
      }
      if (hdf5PathExists(group.get(), "CHUI_SoundSpeedCode")) {
        readDatasetChunk1d(
            group.get(), "CHUI_SoundSpeedCode", 0, local_count,
            gas_sound_speed_chunk);
      } else {
        missing_double_field(
            "/PartType0/CHUI_SoundSpeedCode", local_count,
            gas_sound_speed_chunk, 0.0);
      }
      if (hdf5PathExists(group.get(), "Metallicity")) {
        readDatasetChunk1d(
            group.get(), "Metallicity", 0, local_count,
            gas_metallicity_chunk);
      } else {
        missing_double_field(
            "/PartType0/Metallicity", local_count, gas_metallicity_chunk, 0.0);
      }
      if (hdf5PathExists(group.get(), "GasCellIDs")) {
        readDatasetChunkIds(
            group.get(), "GasCellIDs", 0, local_count, gas_cell_id_chunk);
      } else {
        gas_cell_id_chunk = ids_chunk;
        result.report.defaulted_fields.push_back(
            "/PartType0/GasCellIDs=ParticleIDs");
      }

      if (hdf5PathExists(group.get(), "ParentParticleIDs")) {
        readDatasetChunkIds(
            group.get(), "ParentParticleIDs", 0, local_count,
            gas_parent_particle_id_chunk);
        if (hdf5PathExists(group.get(), "HasParentParticle")) {
          readDatasetChunk1dU8(
              group.get(), "HasParentParticle", 0, local_count,
              gas_has_parent_particle_chunk);
        } else {
          gas_has_parent_particle_chunk.resize(local_count, 0U);
          for (std::size_t i = 0; i < local_count; ++i) {
            gas_has_parent_particle_chunk[i] =
                gas_parent_particle_id_chunk[i] != 0U ? 1U : 0U;
          }
          result.report.defaulted_fields.push_back(
              "/PartType0/HasParentParticle=ParentParticleIDs!=0");
        }
      } else if (chui_authored && schema_version < 6U) {
        // Historical CHUI science snapshots were particle-bound: the PartType0
        // ParticleID was the parent particle identity. This reconstruction is
        // versioned legacy behavior, not a general AREPO/GADGET assumption.
        gas_parent_particle_id_chunk = ids_chunk;
        gas_has_parent_particle_chunk.assign(local_count, 1U);
        result.report.defaulted_fields.push_back(
            "/PartType0/ParentParticleIDs=ParticleIDs(legacy_chui_schema)");
      } else {
        missing_u64_field(
            "/PartType0/ParentParticleIDs", local_count,
            gas_parent_particle_id_chunk, 0U);
        gas_has_parent_particle_chunk.assign(local_count, 0U);
      }

      if (hdf5PathExists(group.get(), "OwningPatchIDs")) {
        readDatasetChunkIds(
            group.get(), "OwningPatchIDs", 0, local_count,
            gas_owning_patch_id_chunk);
      } else {
        missing_u64_field(
            "/PartType0/OwningPatchIDs", local_count,
            gas_owning_patch_id_chunk, 0U);
      }
    }
    if (type_index == 4) {
      if (hdf5PathExists(group.get(), "Metallicity")) {
        readDatasetChunk1d(group.get(), "Metallicity", 0, local_count, star_metallicity_chunk);
      } else {
        missing_double_field("/PartType4/Metallicity", local_count, star_metallicity_chunk, 0.0);
      }
      if (hdf5PathExists(group.get(), "StellarFormationTime")) {
        readDatasetChunk1d(group.get(), "StellarFormationTime", 0, local_count, star_formation_time_chunk);
      } else {
        missing_double_field(
            "/PartType4/StellarFormationTime", local_count, star_formation_time_chunk,
            result.state.metadata.scale_factor);
      }
      if (hdf5PathExists(group.get(), "BirthMass")) {
        readDatasetChunk1d(group.get(), "BirthMass", 0, local_count, star_birth_mass_chunk);
      } else {
        if (options.missing_field_policy == SnapshotMissingFieldPolicy::kReject) {
          throw std::runtime_error("snapshot import: required scientific field is missing: /PartType4/BirthMass");
        }
        star_birth_mass_chunk = mass_chunk;
        if (options.missing_field_policy == SnapshotMissingFieldPolicy::kMarkUnavailable) {
          result.report.unavailable_fields.push_back("/PartType4/BirthMass");
        } else {
          result.report.defaulted_fields.push_back("/PartType4/BirthMass=current_mass_explicit_default");
        }
      }
      if (hdf5PathExists(group.get(), "StarFormationBirthKey")) {
        readDatasetChunkIds(group.get(), "StarFormationBirthKey", 0, local_count, star_birth_key_chunk);
      } else {
        missing_u64_field("/PartType4/StarFormationBirthKey", local_count, star_birth_key_chunk, 0U);
      }
      if (hdf5PathExists(group.get(), "ParentGasCellID")) {
        readDatasetChunkIds(group.get(), "ParentGasCellID", 0, local_count, star_parent_gas_cell_id_chunk);
      } else {
        missing_u64_field("/PartType4/ParentGasCellID", local_count, star_parent_gas_cell_id_chunk, 0U);
      }
      if (hdf5PathExists(group.get(), "BirthIntegrationTick")) {
        readDatasetChunkIds(group.get(), "BirthIntegrationTick", 0, local_count, star_birth_tick_chunk);
      } else {
        missing_u64_field("/PartType4/BirthIntegrationTick", local_count, star_birth_tick_chunk, 0U);
      }
      if (hdf5PathExists(group.get(), "BirthOrdinal")) {
        readDatasetChunkU32(group.get(), "BirthOrdinal", 0, local_count, star_birth_ordinal_chunk);
      } else {
        missing_u32_field("/PartType4/BirthOrdinal", local_count, star_birth_ordinal_chunk, 0U);
      }
      auto read_star_extension = [&](std::string_view name, std::vector<double>& values) {
        if (hdf5PathExists(group.get(), std::string(name))) {
          readDatasetChunk1d(group.get(), std::string(name), 0, local_count, values);
        } else if (chui_authored && schema_version >= 6U) {
          missing_double_field(std::string("/PartType4/") + std::string(name), local_count, values, 0.0);
        } else {
          values.assign(local_count, 0.0);
          result.report.unavailable_fields.push_back(std::string("/PartType4/") + std::string(name));
        }
      };
      read_star_extension("CHUI_StellarAgeYearsLast", star_age_chunk);
      read_star_extension("CHUI_StellarReturnedMassCumulative", star_returned_mass_chunk);
      read_star_extension("CHUI_StellarReturnedMetalsCumulative", star_returned_metals_chunk);
      read_star_extension("CHUI_StellarNewlySynthesizedMetalsCumulative", star_new_metals_chunk);
      read_star_extension("CHUI_StellarFeedbackEnergyCumulativeErg", star_feedback_energy_chunk);
      read_star_extension("CHUI_StellarDepositedMassCumulative", star_deposited_mass_chunk);
      read_star_extension("CHUI_StellarDepositedMetalsCumulative", star_deposited_metals_chunk);
      read_star_extension("CHUI_StellarDepositedFeedbackEnergyCumulativeErg", star_deposited_energy_chunk);
    }
    if (type_index == 3) {
      if (hdf5PathExists(group.get(), "TracerParentParticleID")) {
        readDatasetChunkIds(group.get(), "TracerParentParticleID", 0, local_count, tracer_parent_chunk);
      } else {
        missing_u64_field("/PartType3/TracerParentParticleID", local_count, tracer_parent_chunk, 0U);
      }
      if (hdf5PathExists(group.get(), "TracerInjectionStep")) {
        readDatasetChunkIds(group.get(), "TracerInjectionStep", 0, local_count, tracer_step_chunk);
      } else {
        missing_u64_field("/PartType3/TracerInjectionStep", local_count, tracer_step_chunk, 0U);
      }
      if (hdf5PathExists(group.get(), "TracerHostCellIndex")) {
        readDatasetChunkU32(group.get(), "TracerHostCellIndex", 0, local_count, tracer_host_chunk);
      } else {
        missing_u32_field(
            "/PartType3/TracerHostCellIndex", local_count, tracer_host_chunk,
            core::kInvalidGasCellRow);
      }
      if (hdf5PathExists(group.get(), "TracerMassFractionOfHost")) {
        readDatasetChunk1d(group.get(), "TracerMassFractionOfHost", 0, local_count, tracer_fraction_chunk);
      } else {
        missing_double_field("/PartType3/TracerMassFractionOfHost", local_count, tracer_fraction_chunk, 0.0);
      }
      if (hdf5PathExists(group.get(), "TracerLastHostMassCode")) {
        readDatasetChunk1d(group.get(), "TracerLastHostMassCode", 0, local_count, tracer_last_host_mass_chunk);
      } else {
        missing_double_field("/PartType3/TracerLastHostMassCode", local_count, tracer_last_host_mass_chunk, 0.0);
      }
      if (hdf5PathExists(group.get(), "TracerCumulativeExchangedMassCode")) {
        readDatasetChunk1d(group.get(), "TracerCumulativeExchangedMassCode", 0, local_count, tracer_exchange_chunk);
      } else {
        missing_double_field("/PartType3/TracerCumulativeExchangedMassCode", local_count, tracer_exchange_chunk, 0.0);
      }
    }
    if (type_index == 5) {
      auto read_bh_double = [&](std::string_view name, std::vector<double>& values) {
        if (hdf5PathExists(group.get(), std::string(name))) {
          readDatasetChunk1d(group.get(), std::string(name), 0, local_count, values);
        } else if (chui_authored && schema_version >= 6U) {
          missing_double_field(std::string("/PartType5/") + std::string(name), local_count, values, 0.0);
        } else {
          values.assign(local_count, 0.0);
          result.report.unavailable_fields.push_back(std::string("/PartType5/") + std::string(name));
        }
      };
      read_bh_double("CHUI_BHSubgridMass", bh_subgrid_mass_chunk);
      read_bh_double("CHUI_BHAccretionRateMsunPerYr", bh_accretion_rate_chunk);
      read_bh_double("CHUI_BHFeedbackEnergyCode", bh_feedback_energy_chunk);
      read_bh_double("CHUI_BHEddingtonRatio", bh_eddington_chunk);
      read_bh_double("CHUI_BHCumulativeAccretedMass", bh_cumulative_accreted_chunk);
      read_bh_double("CHUI_BHCumulativeFeedbackEnergyCode", bh_cumulative_feedback_chunk);
      read_bh_double("CHUI_BHDutyCycleActiveTimeCode", bh_duty_active_chunk);
      read_bh_double("CHUI_BHDutyCycleTotalTimeCode", bh_duty_total_chunk);
      if (hdf5PathExists(group.get(), "CHUI_BHHostCellIndex")) {
        readDatasetChunkU32(group.get(), "CHUI_BHHostCellIndex", 0, local_count, bh_host_cell_chunk);
      } else {
        missing_u32_field(
            "/PartType5/CHUI_BHHostCellIndex", local_count, bh_host_cell_chunk,
            core::kInvalidGasCellRow);
      }
    }

    for (std::size_t i = 0; i < local_count; ++i) {
      const std::size_t global_i = global_offset + i;
      double position_x = conversion.positionFromStored(coords_chunk[i * 3 + 0]);
      double position_y = conversion.positionFromStored(coords_chunk[i * 3 + 1]);
      double position_z = conversion.positionFromStored(coords_chunk[i * 3 + 2]);
      const double velocity_x = conversion.velocityFromStored(vel_chunk[i * 3 + 0]);
      const double velocity_y = conversion.velocityFromStored(vel_chunk[i * 3 + 1]);
      const double velocity_z = conversion.velocityFromStored(vel_chunk[i * 3 + 2]);
      const double mass_code = conversion.massFromStored(mass_chunk[i]);
      if (!std::isfinite(position_x) || !std::isfinite(position_y) ||
          !std::isfinite(position_z) || !std::isfinite(velocity_x) ||
          !std::isfinite(velocity_y) || !std::isfinite(velocity_z) ||
          !std::isfinite(mass_code) || mass_code < 0.0) {
        throw std::runtime_error(
            "snapshot reader: non-finite phase-space value or negative mass");
      }
      if (periodic_external_positions) {
        const double wrapped_x = core::canonicalPeriodicPosition(position_x, header_box_size_x_code);
        const double wrapped_y = core::canonicalPeriodicPosition(position_y, header_box_size_y_code);
        const double wrapped_z = core::canonicalPeriodicPosition(position_z, header_box_size_z_code);
        periodically_wrapped_components += static_cast<std::uint64_t>(wrapped_x != position_x);
        periodically_wrapped_components += static_cast<std::uint64_t>(wrapped_y != position_y);
        periodically_wrapped_components += static_cast<std::uint64_t>(wrapped_z != position_z);
        position_x = wrapped_x; position_y = wrapped_y; position_z = wrapped_z;
      }
      std::uint64_t internal_particle_id = ids_chunk[i];
      if (!chui_authored) {
        source_id_has_zero = source_id_has_zero || internal_particle_id == 0U;
        source_id_has_uint64_max = source_id_has_uint64_max ||
            internal_particle_id == std::numeric_limits<std::uint64_t>::max();
      }
      if (chui_authored) {
        const double tolerance = 1.0e-10 * std::max(
            {1.0, header_box_size_x_code,
             header_box_size_y_code, header_box_size_z_code});
        if (position_x < -tolerance ||
            position_x > header_box_size_x_code + tolerance ||
            position_y < -tolerance ||
            position_y > header_box_size_y_code + tolerance ||
            position_z < -tolerance ||
            position_z > header_box_size_z_code + tolerance) {
          throw std::runtime_error(
              "snapshot reader: CHUI-authored coordinate lies outside declared box bounds");
        }
      }
      result.state.particles.position_x_comoving[global_i] = position_x;
      result.state.particles.position_y_comoving[global_i] = position_y;
      result.state.particles.position_z_comoving[global_i] = position_z;
      result.state.particles.velocity_x_peculiar[global_i] = velocity_x;
      result.state.particles.velocity_y_peculiar[global_i] = velocity_y;
      result.state.particles.velocity_z_peculiar[global_i] = velocity_z;
      result.state.particles.mass_code[global_i] = mass_code;
      result.state.particles.time_bin[global_i] = 0;
      result.state.particle_sidecar.particle_id[global_i] = internal_particle_id;
      result.state.particle_sidecar.species_tag[global_i] = mapPartTypeToSpeciesTag(type_index, options, chui_authored);
      result.state.particle_sidecar.owning_rank[global_i] = 0;
      if (!softening_chunk.empty()) {
        if (result.state.particle_sidecar.gravity_softening_comoving.empty()) {
          result.state.particle_sidecar.gravity_softening_comoving.resize(result.state.particles.size(), 0.0);
        }
        result.state.particle_sidecar.gravity_softening_comoving[global_i] = conversion.softeningComovingFromStored(softening_chunk[i]);
        if (!softening_override_mask_chunk.empty()) {
          if (result.state.particle_sidecar.has_gravity_softening_override.empty()) {
            result.state.particle_sidecar.has_gravity_softening_override.resize(result.state.particles.size(), 0U);
          }
          result.state.particle_sidecar.has_gravity_softening_override[global_i] = softening_override_mask_chunk[i];
        }
      }
      result.state.species.count_by_species[result.state.particle_sidecar.species_tag[global_i]] += 1;
      if (type_index == 0) {
        gas_particle_index.push_back(static_cast<std::uint32_t>(global_i));
        const double gas_u = conversion.internalEnergyFromStored(
            gas_internal_energy_chunk[i]);
        const double gas_rho = conversion.densityComovingFromStored(
            gas_density_chunk[i]);
        const double gas_pressure = conversion.pressureComovingFromStored(
            gas_pressure_chunk[i]);
        const double gas_metallicity = gas_metallicity_chunk[i];
        if (!std::isfinite(gas_u) || gas_u < 0.0 ||
            !std::isfinite(gas_rho) || gas_rho < 0.0 ||
            !std::isfinite(gas_pressure) || gas_pressure < 0.0 ||
            !std::isfinite(gas_temperature_chunk[i]) ||
            gas_temperature_chunk[i] < 0.0 ||
            !std::isfinite(gas_sound_speed_chunk[i]) ||
            gas_sound_speed_chunk[i] < 0.0 ||
            !std::isfinite(gas_metallicity) || gas_metallicity < 0.0 ||
            gas_metallicity > 1.0) {
          throw std::runtime_error(
              "snapshot reader: invalid gas thermodynamic/metallicity value");
        }
        gas_internal_energy_code.push_back(gas_u);
        gas_density_code.push_back(gas_rho);
        gas_pressure_code.push_back(gas_pressure);
        gas_temperature_code.push_back(gas_temperature_chunk[i]);
        gas_sound_speed_code.push_back(gas_sound_speed_chunk[i]);
        gas_metallicity_mass_fraction.push_back(gas_metallicity);
        gas_cell_id.push_back(gas_cell_id_chunk[i]);
        gas_parent_particle_id.push_back(gas_parent_particle_id_chunk[i]);
        gas_has_parent_particle.push_back(gas_has_parent_particle_chunk[i]);
        gas_owning_patch_id.push_back(gas_owning_patch_id_chunk[i]);
      }
      if (type_index == 4) {
        star_particle_index.push_back(static_cast<std::uint32_t>(global_i));
        star_metallicity_mass_fraction.push_back(star_metallicity_chunk[i]);
        star_formation_scale_factor.push_back(star_formation_time_chunk[i]);
        star_birth_mass_code.push_back(conversion.massFromStored(star_birth_mass_chunk[i]));
        star_birth_key.push_back(star_birth_key_chunk[i]);
        star_parent_gas_cell_id.push_back(star_parent_gas_cell_id_chunk[i]);
        star_birth_tick.push_back(star_birth_tick_chunk[i]);
        star_birth_ordinal.push_back(star_birth_ordinal_chunk[i]);
        star_age_years_last.push_back(star_age_chunk[i]);
        star_returned_mass_cumulative_code.push_back(conversion.massFromStored(star_returned_mass_chunk[i]));
        star_returned_metals_cumulative_code.push_back(conversion.massFromStored(star_returned_metals_chunk[i]));
        star_newly_synthesized_metals_cumulative_code.push_back(conversion.massFromStored(star_new_metals_chunk[i]));
        star_feedback_energy_cumulative_erg.push_back(star_feedback_energy_chunk[i]);
        star_deposited_mass_cumulative_code.push_back(conversion.massFromStored(star_deposited_mass_chunk[i]));
        star_deposited_metals_cumulative_code.push_back(conversion.massFromStored(star_deposited_metals_chunk[i]));
        star_deposited_feedback_energy_cumulative_erg.push_back(star_deposited_energy_chunk[i]);
      }
      if (type_index == 3) {
        tracer_particle_index.push_back(static_cast<std::uint32_t>(global_i));
        tracer_parent_particle_id.push_back(tracer_parent_chunk[i]);
        tracer_injection_step.push_back(tracer_step_chunk[i]);
        tracer_host_cell_index.push_back(tracer_host_chunk[i]);
        tracer_mass_fraction_of_host.push_back(tracer_fraction_chunk[i]);
        tracer_last_host_mass_code.push_back(conversion.massFromStored(tracer_last_host_mass_chunk[i]));
        tracer_cumulative_exchanged_mass_code.push_back(conversion.massFromStored(tracer_exchange_chunk[i]));
      }
      if (type_index == 5) {
        black_hole_particle_index.push_back(static_cast<std::uint32_t>(global_i));
        black_hole_host_cell_index.push_back(bh_host_cell_chunk[i]);
        black_hole_subgrid_mass_code.push_back(conversion.massFromStored(bh_subgrid_mass_chunk[i]));
        black_hole_accretion_rate_code.push_back(conversion.starFormationRateStoredToCode(bh_accretion_rate_chunk[i]));
        black_hole_feedback_energy_code.push_back(bh_feedback_energy_chunk[i]);
        black_hole_eddington_ratio.push_back(bh_eddington_chunk[i]);
        black_hole_cumulative_accreted_mass_code.push_back(conversion.massFromStored(bh_cumulative_accreted_chunk[i]));
        black_hole_cumulative_feedback_energy_code.push_back(bh_cumulative_feedback_chunk[i]);
        black_hole_duty_cycle_active_time_code.push_back(bh_duty_active_chunk[i]);
        black_hole_duty_cycle_total_time_code.push_back(bh_duty_total_chunk[i]);
      }
    }

    global_offset += local_count;
  }

  source_id_affine_plus_one = !chui_authored && source_id_has_zero;
  if (source_id_affine_plus_one) {
    if (source_id_has_uint64_max) {
      throw std::runtime_error(
          "snapshot reader: external ParticleIDs contain both 0 and UINT64_MAX; +1 normalization would overflow");
    }
    for (std::uint64_t& particle_id : result.state.particle_sidecar.particle_id) {
      ++particle_id;
    }
    for (std::size_t i = 0; i < gas_parent_particle_id.size(); ++i) {
      if (gas_has_parent_particle[i] != 0U) {
        if (gas_parent_particle_id[i] == std::numeric_limits<std::uint64_t>::max()) {
          throw std::runtime_error(
              "snapshot reader: gas ParentParticleIDs overflow +1 external-ID normalization");
        }
        ++gas_parent_particle_id[i];
      }
    }
    for (std::uint64_t& parent_id : tracer_parent_particle_id) {
      if (parent_id == std::numeric_limits<std::uint64_t>::max()) {
        throw std::runtime_error(
            "snapshot reader: tracer parent ParticleID overflows +1 external-ID normalization");
      }
      ++parent_id;
    }
  }

  result.state.metadata.run_name = config.output.run_name;
  result.state.resizeCells(gas_particle_index.size());
  for (std::size_t i = 0; i < gas_particle_index.size(); ++i) {
    const std::size_t particle_index = gas_particle_index[i];
    result.state.cells.center_x_comoving[i] = result.state.particles.position_x_comoving[particle_index];
    result.state.cells.center_y_comoving[i] = result.state.particles.position_y_comoving[particle_index];
    result.state.cells.center_z_comoving[i] = result.state.particles.position_z_comoving[particle_index];
    result.state.cells.mass_code[i] = result.state.particles.mass_code[particle_index];
    result.state.cells.time_bin[i] = 0U;
    result.state.cells.patch_index[i] = 0U;
    result.state.gas_cells.gas_cell_id[i] = gas_cell_id[i] != 0U
        ? gas_cell_id[i]
        : result.state.particle_sidecar.particle_id[particle_index];
    if (gas_has_parent_particle[i] > 1U) {
      throw std::runtime_error(
          "snapshot reader: HasParentParticle values must be 0 or 1");
    }
    if (gas_has_parent_particle[i] != 0U &&
        gas_parent_particle_id[i] == 0U) {
      throw std::runtime_error(
          "snapshot reader: present gas parent identity must be nonzero");
    }
    result.state.gas_cells.parent_particle_id[i] =
        gas_has_parent_particle[i] != 0U ? gas_parent_particle_id[i] : 0U;
    result.state.gas_cells.velocity_x_peculiar[i] =
        result.state.particles.velocity_x_peculiar[particle_index];
    result.state.gas_cells.velocity_y_peculiar[i] =
        result.state.particles.velocity_y_peculiar[particle_index];
    result.state.gas_cells.velocity_z_peculiar[i] =
        result.state.particles.velocity_z_peculiar[particle_index];
    result.state.gas_cells.density_code[i] = gas_density_code[i];
    result.state.gas_cells.pressure_code[i] = gas_pressure_code[i];
    result.state.gas_cells.internal_energy_code[i] = gas_internal_energy_code[i];
    result.state.gas_cells.temperature_code[i] = gas_temperature_code[i];
    result.state.gas_cells.sound_speed_code[i] = gas_sound_speed_code[i];
    result.state.gas_cells.metal_mass_code[i] =
        gas_metallicity_mass_fraction[i] * result.state.cells.mass_code[i];
  }
  result.state.star_particles.resize(star_particle_index.size());
  for (std::size_t i = 0; i < star_particle_index.size(); ++i) {
    result.state.star_particles.particle_index[i] = star_particle_index[i];
    result.state.star_particles.formation_scale_factor[i] = star_formation_scale_factor[i];
    result.state.star_particles.birth_mass_code[i] = star_birth_mass_code[i];
    result.state.star_particles.metallicity_mass_fraction[i] =
        std::clamp(star_metallicity_mass_fraction[i], 0.0, 1.0);
    result.state.star_particles.birth_key[i] = star_birth_key[i];
    result.state.star_particles.parent_gas_cell_id[i] = star_parent_gas_cell_id[i];
    result.state.star_particles.birth_tick[i] = star_birth_tick[i];
    result.state.star_particles.birth_ordinal[i] = star_birth_ordinal[i];
    result.state.star_particles.stellar_age_years_last[i] = star_age_years_last[i];
    result.state.star_particles.stellar_returned_mass_cumulative_code[i] = star_returned_mass_cumulative_code[i];
    result.state.star_particles.stellar_returned_metals_cumulative_code[i] = star_returned_metals_cumulative_code[i];
    result.state.star_particles.stellar_newly_synthesized_metals_cumulative_code[i] = star_newly_synthesized_metals_cumulative_code[i];
    result.state.star_particles.stellar_feedback_energy_cumulative_erg[i] = star_feedback_energy_cumulative_erg[i];
    result.state.star_particles.stellar_deposited_mass_cumulative_code[i] = star_deposited_mass_cumulative_code[i];
    result.state.star_particles.stellar_deposited_metals_cumulative_code[i] = star_deposited_metals_cumulative_code[i];
    result.state.star_particles.stellar_deposited_feedback_energy_cumulative_erg[i] = star_deposited_feedback_energy_cumulative_erg[i];
  }
  result.state.tracers.resize(tracer_particle_index.size());
  for (std::size_t i = 0; i < tracer_particle_index.size(); ++i) {
    result.state.tracers.particle_index[i] = tracer_particle_index[i];
    result.state.tracers.parent_particle_id[i] = tracer_parent_particle_id[i];
    result.state.tracers.injection_step[i] = tracer_injection_step[i];
    result.state.tracers.host_cell_index[i] = tracer_host_cell_index[i];
    result.state.tracers.mass_fraction_of_host[i] = tracer_mass_fraction_of_host[i];
    result.state.tracers.last_host_mass_code[i] = tracer_last_host_mass_code[i];
    result.state.tracers.cumulative_exchanged_mass_code[i] = tracer_cumulative_exchanged_mass_code[i];
  }
  result.state.black_holes.resize(black_hole_particle_index.size());
  for (std::size_t i = 0; i < black_hole_particle_index.size(); ++i) {
    result.state.black_holes.particle_index[i] = black_hole_particle_index[i];
    result.state.black_holes.host_cell_index[i] = black_hole_host_cell_index[i];
    result.state.black_holes.subgrid_mass_code[i] = black_hole_subgrid_mass_code[i];
    result.state.black_holes.accretion_rate_code[i] = black_hole_accretion_rate_code[i];
    result.state.black_holes.feedback_energy_code[i] = black_hole_feedback_energy_code[i];
    result.state.black_holes.eddington_ratio[i] = black_hole_eddington_ratio[i];
    result.state.black_holes.cumulative_accreted_mass_code[i] = black_hole_cumulative_accreted_mass_code[i];
    result.state.black_holes.cumulative_feedback_energy_code[i] = black_hole_cumulative_feedback_energy_code[i];
    result.state.black_holes.duty_cycle_active_time_code[i] = black_hole_duty_cycle_active_time_code[i];
    result.state.black_holes.duty_cycle_total_time_code[i] = black_hole_duty_cycle_total_time_code[i];
  }
  result.state.rebuildSpeciesIndex();
  if (!gas_particle_index.empty()) {
    std::vector<core::GasCellIdentityRecord> identity_records;
    identity_records.reserve(gas_particle_index.size());
    for (std::size_t i = 0; i < gas_particle_index.size(); ++i) {
      identity_records.push_back(core::GasCellIdentityRecord{
          .gas_cell_id = result.state.gas_cells.gas_cell_id[i],
          .parent_particle_id = gas_has_parent_particle[i] != 0U
              ? std::optional<std::uint64_t>(gas_parent_particle_id[i])
              : std::nullopt,
          .owning_patch_id = gas_owning_patch_id[i],
          .local_cell_row = static_cast<std::uint32_t>(i),
      });
    }
    result.state.restoreGasCellIdentityRecords(std::move(identity_records), 1U);
  }

  Hdf5Handle config_group(H5Gopen2(file.get(), std::string(schema.config_group).c_str(), H5P_DEFAULT));
  if (config_group.valid()) {
    readScalarStringAttribute(
        config_group.get(),
        std::string(schema.config_normalized_attribute),
        result.normalized_config_text);
  }

  Hdf5Handle provenance_group(H5Gopen2(file.get(), std::string(schema.provenance_group).c_str(), H5P_DEFAULT));
  if (provenance_group.valid()) {
    readScalarStringAttribute(provenance_group.get(), "schema_version", result.provenance.schema_version);
    readScalarStringAttribute(provenance_group.get(), "config_schema_name", result.provenance.config_schema_name);
    readScalarStringAttribute(
        provenance_group.get(), "config_schema_version", result.provenance.config_schema_version);
    readScalarStringAttribute(provenance_group.get(), "git_sha", result.provenance.git_sha);
    readScalarStringAttribute(provenance_group.get(), "compiler_id", result.provenance.compiler_id);
    readScalarStringAttribute(provenance_group.get(), "compiler_version", result.provenance.compiler_version);
    readScalarStringAttribute(provenance_group.get(), "build_preset", result.provenance.build_preset);
    readScalarStringAttribute(provenance_group.get(), "enabled_features", result.provenance.enabled_features);
    readScalarStringAttribute(provenance_group.get(), "config_hash_hex", result.provenance.config_hash_hex);
    readScalarStringAttribute(
        provenance_group.get(), "normalized_config_hash_hex", result.provenance.normalized_config_hash_hex);
    readScalarStringAttribute(provenance_group.get(), "raw_input_config", result.provenance.raw_input_config);
    readScalarStringAttribute(provenance_group.get(), "normalized_config", result.provenance.normalized_config);
    readScalarStringAttribute(
        provenance_group.get(), "derived_runtime_state", result.provenance.derived_runtime_state);
    readScalarStringAttribute(provenance_group.get(), "timestamp_utc", result.provenance.timestamp_utc);
    readScalarStringAttribute(provenance_group.get(), "hardware_summary", result.provenance.hardware_summary);
    if (result.provenance.schema_version == "provenance_v7") {
      readScalarStringAttribute(provenance_group.get(), "integrity_digest_algorithm", result.provenance.integrity_digest_algorithm);
      readScalarStringAttribute(provenance_group.get(), "normalized_config_sha256_hex", result.provenance.normalized_config_sha256_hex);
      readScalarStringAttribute(provenance_group.get(), "compiler_flags", result.provenance.compiler_flags);
      readScalarStringAttribute(provenance_group.get(), "cpu_model", result.provenance.cpu_model);
      static_cast<void>(readScalarUint32Attribute(provenance_group.get(), "logical_thread_count", result.provenance.logical_thread_count));
      static_cast<void>(readScalarUint32Attribute(provenance_group.get(), "physical_core_count", result.provenance.physical_core_count));
      static_cast<void>(readScalarUint64Attribute(provenance_group.get(), "system_ram_bytes", result.provenance.system_ram_bytes));
      readScalarStringAttribute(provenance_group.get(), "host_name", result.provenance.host_name);
      readScalarStringAttribute(provenance_group.get(), "gpu_summary", result.provenance.gpu_summary);
      readScalarStringAttribute(provenance_group.get(), "cuda_runtime_version", result.provenance.cuda_runtime_version);
      readScalarStringAttribute(provenance_group.get(), "cuda_driver_version", result.provenance.cuda_driver_version);
      readScalarStringAttribute(provenance_group.get(), "mpi_summary", result.provenance.mpi_summary);
      std::uint32_t mpi_world_size = 0U;
      std::uint32_t mpi_node_local_rank = 0U;
      if (readScalarUint32Attribute(provenance_group.get(), "mpi_world_size", mpi_world_size)) {
        result.provenance.mpi_world_size = cosmosim::core::checkedIntegralNarrow<int>(mpi_world_size, "snapshot provenance mpi_world_size");
      }
      if (readScalarUint32Attribute(provenance_group.get(), "mpi_node_local_rank", mpi_node_local_rank)) {
        result.provenance.mpi_node_local_rank = cosmosim::core::checkedIntegralNarrow<int>(mpi_node_local_rank, "snapshot provenance mpi_node_local_rank");
      }
      readScalarStringAttribute(provenance_group.get(), "deterministic_mode", result.provenance.deterministic_mode);
    }
    std::uint32_t gravity_pm_grid = 0;
    if (readScalarUint32Attribute(provenance_group.get(), "gravity_treepm_pm_grid", gravity_pm_grid)) {
      result.provenance.gravity_treepm_pm_grid = static_cast<int>(gravity_pm_grid);
    }
    std::uint32_t gravity_pm_grid_nx = 0;
    std::uint32_t gravity_pm_grid_ny = 0;
    std::uint32_t gravity_pm_grid_nz = 0;
    if (readScalarUint32Attribute(provenance_group.get(), "gravity_treepm_pm_grid_nx", gravity_pm_grid_nx)) {
      result.provenance.gravity_treepm_pm_grid_nx = static_cast<int>(gravity_pm_grid_nx);
    }
    if (readScalarUint32Attribute(provenance_group.get(), "gravity_treepm_pm_grid_ny", gravity_pm_grid_ny)) {
      result.provenance.gravity_treepm_pm_grid_ny = static_cast<int>(gravity_pm_grid_ny);
    }
    if (readScalarUint32Attribute(provenance_group.get(), "gravity_treepm_pm_grid_nz", gravity_pm_grid_nz)) {
      result.provenance.gravity_treepm_pm_grid_nz = static_cast<int>(gravity_pm_grid_nz);
    }
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_treepm_assignment_scheme",
        result.provenance.gravity_treepm_assignment_scheme);
    std::string gravity_deconvolution;
    readScalarStringAttribute(
        provenance_group.get(), "gravity_treepm_window_deconvolution", gravity_deconvolution);
    if (!gravity_deconvolution.empty()) {
      result.provenance.gravity_treepm_window_deconvolution = gravity_deconvolution == "true";
    }
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_asmth_cells",
        result.provenance.gravity_treepm_asmth_cells));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_rcut_cells",
        result.provenance.gravity_treepm_rcut_cells));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_mesh_spacing_mpc_comoving",
        result.provenance.gravity_treepm_mesh_spacing_mpc_comoving));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_mesh_spacing_x_mpc_comoving",
        result.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_mesh_spacing_y_mpc_comoving",
        result.provenance.gravity_treepm_mesh_spacing_y_mpc_comoving));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_mesh_spacing_z_mpc_comoving",
        result.provenance.gravity_treepm_mesh_spacing_z_mpc_comoving));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_split_scale_mpc_comoving",
        result.provenance.gravity_treepm_split_scale_mpc_comoving));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_cutoff_radius_mpc_comoving",
        result.provenance.gravity_treepm_cutoff_radius_mpc_comoving));
    std::uint32_t gravity_pm_cadence = 0;
    if (readScalarUint32Attribute(
            provenance_group.get(), "gravity_treepm_update_cadence_steps", gravity_pm_cadence)) {
      result.provenance.gravity_treepm_update_cadence_steps = static_cast<int>(gravity_pm_cadence);
    }
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_treepm_tree_opening_criterion",
        result.provenance.gravity_treepm_tree_opening_criterion);
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_tree_opening_theta",
        result.provenance.gravity_treepm_tree_opening_theta));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_tree_relative_force_tolerance",
        result.provenance.gravity_treepm_tree_relative_force_tolerance));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_tree_relative_force_acceleration_floor",
        result.provenance.gravity_treepm_tree_relative_force_acceleration_floor));
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_softening_policy",
        result.provenance.gravity_softening_policy);
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_softening_kernel",
        result.provenance.gravity_softening_kernel);
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_softening_epsilon_kpc_comoving",
        result.provenance.gravity_softening_epsilon_kpc_comoving));
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_pm_fft_backend",
        result.provenance.gravity_pm_fft_backend);
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_pm_backend_capability",
        result.provenance.gravity_pm_backend_capability);
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_acceptance_profile_id",
        result.provenance.gravity_acceptance_profile_id);
    if (result.provenance.gravity_treepm_pm_grid_nx == 0 && result.provenance.gravity_treepm_pm_grid > 0) {
      result.provenance.gravity_treepm_pm_grid_nx = result.provenance.gravity_treepm_pm_grid;
      result.provenance.gravity_treepm_pm_grid_ny = result.provenance.gravity_treepm_pm_grid;
      result.provenance.gravity_treepm_pm_grid_nz = result.provenance.gravity_treepm_pm_grid;
    }
    if (result.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving == 0.0 &&
        result.provenance.gravity_treepm_mesh_spacing_mpc_comoving > 0.0) {
      result.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving =
          result.provenance.gravity_treepm_mesh_spacing_mpc_comoving;
      result.provenance.gravity_treepm_mesh_spacing_y_mpc_comoving =
          result.provenance.gravity_treepm_mesh_spacing_mpc_comoving;
      result.provenance.gravity_treepm_mesh_spacing_z_mpc_comoving =
          result.provenance.gravity_treepm_mesh_spacing_mpc_comoving;
    }
  }

  const core::MemoryReport memory_report = core::collectSimulationMemoryReport(result.state);
  result.report.materialized_state_bytes = memory_report.totals.persistent_total_bytes;
  if (result.report.materialized_state_bytes > options.budget.max_materialized_bytes) {
    throw std::length_error(
        "snapshot reader: actual materialized SimulationState exceeds max_materialized_bytes");
  }
  if (source_id_affine_plus_one) {
    result.report.present_aliases.push_back("external_id_mapping=zero_present_plus_one_v1");
  }
  if (periodically_wrapped_components > 0U) {
    result.report.present_aliases.push_back(
        "periodic_position_components_wrapped=" + std::to_string(periodically_wrapped_components));
  }
  internal::updateSnapshotReadiness(result.state, &result.report);
  return result;
#endif
}

}  // namespace cosmosim::io
