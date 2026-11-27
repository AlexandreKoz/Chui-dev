#include "cosmosim/io/mhd_restart_state.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/build_config.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::io {
namespace {

[[noreturn]] void fail(const std::string& detail) {
  throw std::runtime_error("MHD_RESTART_SCHEMA_MISMATCH: " + detail);
}

#if COSMOSIM_ENABLE_HDF5

class Handle {
 public:
  using Closer = herr_t (*)(hid_t);
  Handle(hid_t id, Closer closer) : m_id(id), m_closer(closer) {}
  ~Handle() {
    if (m_id >= 0) {
      m_closer(m_id);
    }
  }
  Handle(const Handle&) = delete;
  Handle& operator=(const Handle&) = delete;
  [[nodiscard]] hid_t get() const noexcept { return m_id; }
  [[nodiscard]] bool valid() const noexcept { return m_id >= 0; }

 private:
  hid_t m_id;
  Closer m_closer;
};

template <typename T>
[[nodiscard]] hid_t nativeType();
template <>
hid_t nativeType<std::uint64_t>() { return H5T_NATIVE_UINT64; }
template <>
hid_t nativeType<std::uint32_t>() { return H5T_NATIVE_UINT32; }
template <>
hid_t nativeType<std::uint8_t>() { return H5T_NATIVE_UINT8; }
template <>
hid_t nativeType<double>() { return H5T_NATIVE_DOUBLE; }

template <typename T>
void writeDataset(hid_t group, const char* name, const T* data, std::uint64_t rows,
                  std::uint64_t columns) {
  const std::array<hsize_t, 2> dims{static_cast<hsize_t>(rows), static_cast<hsize_t>(columns)};
  const int rank = columns == 1U ? 1 : 2;
  Handle space(H5Screate_simple(rank, dims.data(), nullptr), H5Sclose);
  Handle dataset(H5Dcreate2(group, name, nativeType<T>(), space.get(), H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT),
                 H5Dclose);
  if (!space.valid() || !dataset.valid()) {
    fail(std::string("cannot create dataset ") + name);
  }
  if (rows * columns == 0U) {
    return;
  }
  // Bounded staging: hyperslabs of at most k_mhd_restart_chunk_bytes rows.
  const std::uint64_t row_bytes = columns * sizeof(T);
  const std::uint64_t rows_per_chunk =
      std::max<std::uint64_t>(1U, k_mhd_restart_chunk_bytes / std::max<std::uint64_t>(1U, row_bytes));
  for (std::uint64_t begin = 0; begin < rows; begin += rows_per_chunk) {
    const std::uint64_t count = std::min(rows_per_chunk, rows - begin);
    const std::array<hsize_t, 2> start{static_cast<hsize_t>(begin), 0};
    const std::array<hsize_t, 2> extent{static_cast<hsize_t>(count), static_cast<hsize_t>(columns)};
    Handle file_space(H5Dget_space(dataset.get()), H5Sclose);
    if (H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, start.data(), nullptr,
                            extent.data(), nullptr) < 0) {
      fail(std::string("cannot select hyperslab for ") + name);
    }
    Handle memory_space(H5Screate_simple(rank, extent.data(), nullptr), H5Sclose);
    if (H5Dwrite(dataset.get(), nativeType<T>(), memory_space.get(), file_space.get(),
                 H5P_DEFAULT, data + begin * columns) < 0) {
      fail(std::string("cannot write dataset ") + name);
    }
  }
}

template <typename T>
[[nodiscard]] std::vector<T> readDataset(hid_t group, const char* name, std::uint64_t rows,
                                         std::uint64_t columns) {
  if (H5Lexists(group, name, H5P_DEFAULT) <= 0) {
    fail(std::string("missing dataset /state/mhd/") + name);
  }
  Handle dataset(H5Dopen2(group, name, H5P_DEFAULT), H5Dclose);
  Handle space(H5Dget_space(dataset.get()), H5Sclose);
  const int rank = H5Sget_simple_extent_ndims(space.get());
  std::array<hsize_t, 2> dims{0, 1};
  if (rank < 1 || rank > 2) {
    fail(std::string("dataset ") + name + " has an unsupported rank");
  }
  H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr);
  const std::uint64_t found_columns = rank == 1 ? 1U : dims[1];
  if (dims[0] != rows || found_columns != columns) {
    fail(std::string("dataset ") + name + " has a shape inconsistent with the patch table");
  }
  Handle type(H5Dget_type(dataset.get()), H5Tclose);
  if (H5Tget_size(type.get()) != sizeof(T) ||
      H5Tget_class(type.get()) != H5Tget_class(nativeType<T>())) {
    fail(std::string("dataset ") + name + " has the wrong element type");
  }
  std::vector<T> values(static_cast<std::size_t>(rows * columns));
  if (!values.empty() && H5Dread(dataset.get(), nativeType<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                 values.data()) < 0) {
    fail(std::string("cannot read dataset ") + name + " (truncated or corrupt)");
  }
  return values;
}

// Stream one face-flux lane directly into the candidate (final) storage in
// hyperslabs of at most k_mhd_restart_chunk_bytes (MP08): shape and element
// type are validated BEFORE the destination is sized; there is no
// intermediate full-size copy.
void readFluxLaneInto(hid_t group, const char* name, std::uint64_t expected_count,
                      core::AlignedVector<double>& lane) {
  if (H5Lexists(group, name, H5P_DEFAULT) <= 0) {
    fail(std::string("missing dataset /state/mhd/") + name);
  }
  Handle dataset(H5Dopen2(group, name, H5P_DEFAULT), H5Dclose);
  Handle space(H5Dget_space(dataset.get()), H5Sclose);
  if (H5Sget_simple_extent_ndims(space.get()) != 1) {
    fail(std::string("dataset ") + name + " must be one-dimensional");
  }
  hsize_t dims = 0;
  H5Sget_simple_extent_dims(space.get(), &dims, nullptr);
  if (dims != expected_count) {
    fail(std::string("dataset ") + name + " length disagrees with the patch offset table");
  }
  Handle type(H5Dget_type(dataset.get()), H5Tclose);
  if (H5Tget_size(type.get()) != sizeof(double) || H5Tget_class(type.get()) != H5T_FLOAT) {
    fail(std::string("dataset ") + name + " has the wrong element type");
  }
  if (expected_count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    fail(std::string("dataset ") + name + " exceeds the addressable size");
  }
  lane.resize(static_cast<std::size_t>(expected_count));
  const std::uint64_t chunk = std::max<std::uint64_t>(1U, k_mhd_restart_chunk_bytes / sizeof(double));
  for (std::uint64_t begin = 0; begin < expected_count; begin += chunk) {
    const hsize_t start = begin;
    const hsize_t count = std::min<std::uint64_t>(chunk, expected_count - begin);
    Handle file_space(H5Dget_space(dataset.get()), H5Sclose);
    if (H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, &start, nullptr, &count, nullptr) < 0) {
      fail(std::string("cannot select a hyperslab of ") + name);
    }
    Handle memory_space(H5Screate_simple(1, &count, nullptr), H5Sclose);
    if (H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, memory_space.get(), file_space.get(), H5P_DEFAULT,
                lane.data() + begin) < 0) {
      fail(std::string("cannot read dataset ") + name + " (truncated or corrupt)");
    }
  }
}

void writeStringAttribute(hid_t object, const char* name, const std::string& value) {
  Handle type(H5Tcopy(H5T_C_S1), H5Tclose);
  H5Tset_size(type.get(), std::max<std::size_t>(1U, value.size()));
  H5Tset_strpad(type.get(), H5T_STR_NULLPAD);
  Handle space(H5Screate(H5S_SCALAR), H5Sclose);
  Handle attribute(H5Acreate2(object, name, type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT),
                   H5Aclose);
  const std::string padded = value.empty() ? std::string(1, '\0') : value;
  if (!attribute.valid() || H5Awrite(attribute.get(), type.get(), padded.data()) < 0) {
    fail(std::string("cannot write attribute ") + name);
  }
}

[[nodiscard]] std::string readStringAttribute(hid_t object, const char* name) {
  if (H5Aexists(object, name) <= 0) {
    fail(std::string("missing attribute ") + name);
  }
  Handle attribute(H5Aopen(object, name, H5P_DEFAULT), H5Aclose);
  Handle type(H5Aget_type(attribute.get()), H5Tclose);
  if (H5Tget_class(type.get()) != H5T_STRING || H5Tis_variable_str(type.get()) > 0) {
    fail(std::string("attribute ") + name + " must be a fixed-length string");
  }
  const std::size_t size = H5Tget_size(type.get());
  std::vector<char> buffer(size + 1U, '\0');
  if (H5Aread(attribute.get(), type.get(), buffer.data()) < 0) {
    fail(std::string("cannot read attribute ") + name);
  }
  return std::string(buffer.data(), strnlen(buffer.data(), size));
}

template <typename T>
void writeScalarAttribute(hid_t object, const char* name, T value) {
  Handle space(H5Screate(H5S_SCALAR), H5Sclose);
  Handle attribute(H5Acreate2(object, name, nativeType<T>(), space.get(), H5P_DEFAULT, H5P_DEFAULT),
                   H5Aclose);
  if (!attribute.valid() || H5Awrite(attribute.get(), nativeType<T>(), &value) < 0) {
    fail(std::string("cannot write attribute ") + name);
  }
}

template <typename T>
[[nodiscard]] T readScalarAttribute(hid_t object, const char* name) {
  if (H5Aexists(object, name) <= 0) {
    fail(std::string("missing attribute ") + name);
  }
  Handle attribute(H5Aopen(object, name, H5P_DEFAULT), H5Aclose);
  T value{};
  if (H5Aread(attribute.get(), nativeType<T>(), &value) < 0) {
    fail(std::string("cannot read attribute ") + name);
  }
  return value;
}

[[nodiscard]] std::uint64_t checkedMul(std::uint64_t a, std::uint64_t b) {
  if (a != 0U && b > std::numeric_limits<std::uint64_t>::max() / a) {
    fail("face count overflow");
  }
  return a * b;
}

[[nodiscard]] std::uint64_t checkedAdd(std::uint64_t a, std::uint64_t b) {
  if (b > std::numeric_limits<std::uint64_t>::max() - a) {
    fail("face offset overflow");
  }
  return a + b;
}

#endif  // COSMOSIM_ENABLE_HDF5

}  // namespace

void writeMhdRestartGroup(std::int64_t hdf5_file, const core::MagneticFaceSoa& faces,
                          const MhdRestartGroupMetadata& metadata) {
#if !COSMOSIM_ENABLE_HDF5
  (void)hdf5_file;
  (void)faces;
  (void)metadata;
  fail("MHD restart group requires COSMOSIM_ENABLE_HDF5=ON");
#else
  if (!faces.isConsistent()) {
    fail("refusing to write an inconsistent magnetic face state");
  }
  const bool empty_rank = faces.patches.empty();
  if (empty_rank && (!metadata.allow_empty_rank || !faces.face_flux_x_comoving.empty() ||
                     !faces.face_flux_y_comoving.empty() || !faces.face_flux_z_comoving.empty())) {
    fail("refusing to write an empty magnetic face state outside a run-wide magnetized "
         "checkpoint (or with orphan flux values)");
  }
  const auto finite = [](const auto& lane) {
    return std::all_of(lane.begin(), lane.end(), [](double v) { return std::isfinite(v); });
  };
  if (!finite(faces.face_flux_x_comoving) || !finite(faces.face_flux_y_comoving) ||
      !finite(faces.face_flux_z_comoving) || !std::isfinite(metadata.divergence_max_norm) ||
      !(metadata.scale_factor > 0.0)) {
    fail("refusing to write non-finite magnetic authority or metadata");
  }
  const hid_t file = static_cast<hid_t>(hdf5_file);
  if (H5Lexists(file, "/state", H5P_DEFAULT) <= 0) {
    Handle state(H5Gcreate2(file, "/state", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    if (!state.valid()) {
      fail("cannot create /state");
    }
  }
  Handle group(H5Gcreate2(file, "/state/mhd", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
  if (!group.valid()) {
    fail("cannot create /state/mhd (already present?)");
  }
  const std::uint64_t count = faces.patches.size();
  std::vector<std::uint64_t> ids(count);
  std::vector<std::uint64_t> generations(count);
  std::vector<std::uint8_t> levels(count);
  std::vector<std::uint32_t> dims(3U * count);
  std::vector<double> origins(3U * count);
  std::vector<double> widths(3U * count);
  std::array<std::vector<std::uint64_t>, 3> offsets{};
  for (auto& lane : offsets) {
    lane.assign(count + 1U, 0U);
  }
  for (std::size_t index = 0; index < count; ++index) {
    const core::MagneticFacePatchRecord& record = faces.patches[index];
    if (record.level < 0 || record.level > 255) {
      fail("refinement level outside the uint8 schema range");
    }
    ids[index] = record.patch_id;
    generations[index] = record.geometry_generation;
    levels[index] = static_cast<std::uint8_t>(record.level);
    dims[3U * index + 0U] = record.nx;
    dims[3U * index + 1U] = record.ny;
    dims[3U * index + 2U] = record.nz;
    origins[3U * index + 0U] = record.origin_x_comoving;
    origins[3U * index + 1U] = record.origin_y_comoving;
    origins[3U * index + 2U] = record.origin_z_comoving;
    widths[3U * index + 0U] = record.cell_width_x_comoving;
    widths[3U * index + 1U] = record.cell_width_y_comoving;
    widths[3U * index + 2U] = record.cell_width_z_comoving;
    offsets[0][index] = record.flux_x_offset;
    offsets[1][index] = record.flux_y_offset;
    offsets[2][index] = record.flux_z_offset;
  }
  offsets[0][count] = faces.face_flux_x_comoving.size();
  offsets[1][count] = faces.face_flux_y_comoving.size();
  offsets[2][count] = faces.face_flux_z_comoving.size();
  writeScalarAttribute<std::uint32_t>(group.get(), "schema_version",
                                      k_mhd_restart_group_schema_version);
  writeStringAttribute(group.get(), "reference_contract", metadata.reference_contract);
  writeStringAttribute(group.get(), "coefficient_table_sha256", metadata.coefficient_table_sha256);
  writeStringAttribute(group.get(), "magnetic_units", std::string(k_mhd_restart_magnetic_units));
  writeStringAttribute(group.get(), "face_flux_contract", faces.contract);
  writeScalarAttribute<double>(group.get(), "divergence_max_norm", metadata.divergence_max_norm);
  writeScalarAttribute<double>(group.get(), "scale_factor", metadata.scale_factor);
  if (empty_rank) {
    writeScalarAttribute<std::uint32_t>(group.get(), "empty_rank", 1U);
  }
  writeDataset(group.get(), "patch_id", ids.data(), count, 1U);
  writeDataset(group.get(), "geometry_generation", generations.data(), count, 1U);
  writeDataset(group.get(), "refinement_level", levels.data(), count, 1U);
  writeDataset(group.get(), "face_dims_xyz", dims.data(), count, 3U);
  writeDataset(group.get(), "patch_origin_comoving", origins.data(), count, 3U);
  writeDataset(group.get(), "cell_width_comoving", widths.data(), count, 3U);
  writeDataset(group.get(), "face_offset_x", offsets[0].data(), count + 1U, 1U);
  writeDataset(group.get(), "face_offset_y", offsets[1].data(), count + 1U, 1U);
  writeDataset(group.get(), "face_offset_z", offsets[2].data(), count + 1U, 1U);
  writeDataset(group.get(), "face_flux_x_comoving", faces.face_flux_x_comoving.data(),
               faces.face_flux_x_comoving.size(), 1U);
  writeDataset(group.get(), "face_flux_y_comoving", faces.face_flux_y_comoving.data(),
               faces.face_flux_y_comoving.size(), 1U);
  writeDataset(group.get(), "face_flux_z_comoving", faces.face_flux_z_comoving.data(),
               faces.face_flux_z_comoving.size(), 1U);
#endif
}

core::MagneticFaceSoa readMhdRestartGroup(std::int64_t hdf5_file,
                                          const MhdRestartExpectation& expected,
                                          MhdRestartGroupMetadata* metadata_out) {
#if !COSMOSIM_ENABLE_HDF5
  (void)hdf5_file;
  (void)expected;
  (void)metadata_out;
  fail("MHD restart group requires COSMOSIM_ENABLE_HDF5=ON");
#else
  const hid_t file = static_cast<hid_t>(hdf5_file);
  if (H5Lexists(file, "/state", H5P_DEFAULT) <= 0 ||
      H5Lexists(file, "/state/mhd", H5P_DEFAULT) <= 0) {
    fail("v26 checkpoint lacks the authoritative /state/mhd group");
  }
  Handle group(H5Gopen2(file, "/state/mhd", H5P_DEFAULT), H5Gclose);
  if (readScalarAttribute<std::uint32_t>(group.get(), "schema_version") !=
      k_mhd_restart_group_schema_version) {
    fail("/state/mhd group schema_version mismatch");
  }
  MhdRestartGroupMetadata metadata{};
  metadata.reference_contract = readStringAttribute(group.get(), "reference_contract");
  metadata.coefficient_table_sha256 = readStringAttribute(group.get(), "coefficient_table_sha256");
  metadata.divergence_max_norm = readScalarAttribute<double>(group.get(), "divergence_max_norm");
  metadata.scale_factor = readScalarAttribute<double>(group.get(), "scale_factor");
  if (readStringAttribute(group.get(), "magnetic_units") != k_mhd_restart_magnetic_units) {
    fail("magnetic units/normalization mismatch");
  }
  const std::string face_flux_contract = readStringAttribute(group.get(), "face_flux_contract");
  if (metadata.reference_contract != expected.reference_contract) {
    fail("reference contract mismatch (found " + metadata.reference_contract + ")");
  }
  if (metadata.coefficient_table_sha256 != expected.coefficient_table_sha256) {
    fail("coefficient table SHA-256 mismatch");
  }
  if (!(metadata.scale_factor > 0.0) ||
      std::abs(metadata.scale_factor - expected.scale_factor) >
          expected.scale_factor_relative_tolerance * expected.scale_factor) {
    fail("scale factor of the checkpoint does not match the restart timeline");
  }
  if (H5Lexists(group.get(), "patch_id", H5P_DEFAULT) <= 0) {
    fail("missing dataset /state/mhd/patch_id");
  }
  std::uint64_t count = 0U;
  {
    Handle dataset(H5Dopen2(group.get(), "patch_id", H5P_DEFAULT), H5Dclose);
    Handle space(H5Dget_space(dataset.get()), H5Sclose);
    hsize_t dims = 0;
    if (H5Sget_simple_extent_ndims(space.get()) != 1) {
      fail("patch_id must be one-dimensional");
    }
    H5Sget_simple_extent_dims(space.get(), &dims, nullptr);
    count = dims;
  }
  const bool empty_rank_marker = H5Aexists(group.get(), "empty_rank") > 0 &&
      readScalarAttribute<std::uint32_t>(group.get(), "empty_rank") == 1U;
  if (H5Aexists(group.get(), "empty_rank") > 0 && !empty_rank_marker) {
    fail("invalid empty_rank marker value");
  }
  if (count == 0U && !empty_rank_marker) {
    fail("magnetized checkpoint carries no patches and no empty_rank marker (truncated or "
         "corrupt group)");
  }
  if (count != 0U && empty_rank_marker) {
    fail("empty_rank marker on a non-empty patch table (corrupt group)");
  }
  // Bound the patch table before allocating it (malicious/corrupt extents).
  if (count > k_mhd_restart_max_patches_per_rank) {
    fail("patch table exceeds the per-rank bound of the v26 schema");
  }
  // Canonical dataset spelling (pass-02 schema decision R-SPELL): the
  // repository naming rulebook's _comoving suffix; REF section 13.2's
  // abbreviated _comov spelling is not an accepted alias.
  for (const char* superseded :
       {"face_flux_x_comov", "face_flux_y_comov", "face_flux_z_comov"}) {
    if (H5Lexists(group.get(), superseded, H5P_DEFAULT) > 0) {
      fail(std::string("non-canonical dataset spelling /state/mhd/") + superseded +
           " (canonical v26 names are face_flux_{x,y,z}_comoving)");
    }
  }
  const auto ids = readDataset<std::uint64_t>(group.get(), "patch_id", count, 1U);
  const auto generations = readDataset<std::uint64_t>(group.get(), "geometry_generation", count, 1U);
  const auto levels = readDataset<std::uint8_t>(group.get(), "refinement_level", count, 1U);
  const auto dims = readDataset<std::uint32_t>(group.get(), "face_dims_xyz", count, 3U);
  const auto origins = readDataset<double>(group.get(), "patch_origin_comoving", count, 3U);
  const auto widths = readDataset<double>(group.get(), "cell_width_comoving", count, 3U);
  const auto offset_x = readDataset<std::uint64_t>(group.get(), "face_offset_x", count + 1U, 1U);
  const auto offset_y = readDataset<std::uint64_t>(group.get(), "face_offset_y", count + 1U, 1U);
  const auto offset_z = readDataset<std::uint64_t>(group.get(), "face_offset_z", count + 1U, 1U);

  core::MagneticFaceSoa faces{};
  faces.contract = face_flux_contract;
  faces.coefficient_table_sha256 = metadata.coefficient_table_sha256;
  faces.divergence_max_normalized = metadata.divergence_max_norm;
  faces.patches.resize(static_cast<std::size_t>(count));
  std::set<std::uint64_t> seen;
  std::array<std::uint64_t, 3> cursor{};
  if (offset_x[0] != 0U || offset_y[0] != 0U || offset_z[0] != 0U) {
    fail("face offsets must start at zero");
  }
  for (std::size_t index = 0; index < count; ++index) {
    core::MagneticFacePatchRecord& record = faces.patches[index];
    record.patch_id = ids[index];
    record.geometry_generation = generations[index];
    record.level = levels[index];
    record.nx = dims[3U * index + 0U];
    record.ny = dims[3U * index + 1U];
    record.nz = dims[3U * index + 2U];
    record.origin_x_comoving = origins[3U * index + 0U];
    record.origin_y_comoving = origins[3U * index + 1U];
    record.origin_z_comoving = origins[3U * index + 2U];
    record.cell_width_x_comoving = widths[3U * index + 0U];
    record.cell_width_y_comoving = widths[3U * index + 1U];
    record.cell_width_z_comoving = widths[3U * index + 2U];
    if (!seen.insert(record.patch_id).second) {
      fail("duplicate patch_id in /state/mhd");
    }
    const std::uint64_t nx = record.nx;
    const std::uint64_t ny = record.ny;
    const std::uint64_t nz = record.nz;
    const std::array<std::uint64_t, 3> counts{checkedMul(checkedMul(nx + 1U, ny), nz),
                                              checkedMul(checkedMul(nx, ny + 1U), nz),
                                              checkedMul(checkedMul(nx, ny), nz + 1U)};
    const std::array<const std::vector<std::uint64_t>*, 3> offsets{&offset_x, &offset_y, &offset_z};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if ((*offsets[axis])[index] != cursor[axis]) {
        fail("face offsets are not a monotone exclusive prefix sum");
      }
      cursor[axis] = checkedAdd(cursor[axis], counts[axis]);
      if ((*offsets[axis])[index + 1U] != cursor[axis]) {
        fail("face offsets disagree with the patch face counts");
      }
    }
    record.flux_x_offset = offset_x[index];
    record.flux_y_offset = offset_y[index];
    record.flux_z_offset = offset_z[index];
  }
  readFluxLaneInto(group.get(), "face_flux_x_comoving", cursor[0], faces.face_flux_x_comoving);
  readFluxLaneInto(group.get(), "face_flux_y_comoving", cursor[1], faces.face_flux_y_comoving);
  readFluxLaneInto(group.get(), "face_flux_z_comoving", cursor[2], faces.face_flux_z_comoving);
  if (!faces.isConsistent()) {
    fail("restored magnetic face state is inconsistent (geometry/offsets)");
  }
  const auto finite = [](const auto& lane) {
    return std::all_of(lane.begin(), lane.end(), [](double v) { return std::isfinite(v); });
  };
  if (!finite(faces.face_flux_x_comoving) || !finite(faces.face_flux_y_comoving) ||
      !finite(faces.face_flux_z_comoving)) {
    fail("non-finite face flux in the checkpoint");
  }
  // Recompute the discrete divergence of every patch (G03 normalization).
  for (const core::MagneticFacePatchRecord& record : faces.patches) {
    const std::uint64_t nx = record.nx;
    const std::uint64_t ny = record.ny;
    const std::uint64_t nz = record.nz;
    const double* px = faces.face_flux_x_comoving.data() + record.flux_x_offset;
    const double* py = faces.face_flux_y_comoving.data() + record.flux_y_offset;
    const double* pz = faces.face_flux_z_comoving.data() + record.flux_z_offset;
    const std::uint64_t fx = (nx + 1U) * ny * nz;
    const std::uint64_t fy = nx * (ny + 1U) * nz;
    const std::uint64_t fz = nx * ny * (nz + 1U);
    double maximum_flux = 0.0;
    for (std::uint64_t i = 0; i < fx; ++i) maximum_flux = std::max(maximum_flux, std::abs(px[i]));
    for (std::uint64_t i = 0; i < fy; ++i) maximum_flux = std::max(maximum_flux, std::abs(py[i]));
    for (std::uint64_t i = 0; i < fz; ++i) maximum_flux = std::max(maximum_flux, std::abs(pz[i]));
    const double area = std::max({record.cell_width_y_comoving * record.cell_width_z_comoving,
                                  record.cell_width_x_comoving * record.cell_width_z_comoving,
                                  record.cell_width_x_comoving * record.cell_width_y_comoving});
    const double normalizer = std::max(maximum_flux, 1.0e-30 * area);
    double divergence = 0.0;
    for (std::uint64_t k = 0; k < nz; ++k) {
      for (std::uint64_t j = 0; j < ny; ++j) {
        for (std::uint64_t i = 0; i < nx; ++i) {
          const double value = px[i + 1U + (nx + 1U) * (j + ny * k)] - px[i + (nx + 1U) * (j + ny * k)] +
              py[i + nx * (j + 1U + (ny + 1U) * k)] - py[i + nx * (j + (ny + 1U) * k)] +
              pz[i + nx * (j + ny * (k + 1U))] - pz[i + nx * (j + ny * k)];
          divergence = std::max(divergence, std::abs(value));
        }
      }
    }
    if (divergence / normalizer > expected.divergence_tolerance) {
      fail("restored face flux violates the discrete divergence contract");
    }
  }
  if (metadata_out != nullptr) {
    *metadata_out = metadata;
  }
  return faces;
#endif
}

}  // namespace cosmosim::io
