#include "cosmosim/mhd/mhd_transport.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <set>
#include <sstream>
#include <stdexcept>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/units.hpp"
#include "internal/mhd_sha256.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::mhd {
namespace {

// CODATA 2018 values (cgs). Duplicated locally because cosmosim_mhd depends on
// cosmosim_core only; the values match physics/baryon_composition.hpp.
constexpr double k_proton_mass_g = 1.67262192369e-24;
constexpr double k_boltzmann_erg_per_k = 1.380649e-16;
// Field regularization B_reg = 1e-30 * B_unit_code (section 4.1).
constexpr double k_b_regularization_code = k_mhd_b_zero_regularization_fraction;

[[noreturn]] void throwCoefficient(std::string owner, std::string detail) {
  MhdFailureContext context{};
  context.owner = std::move(owner);
  context.detail = std::move(detail);
  throw MhdError(MhdErrorCode::kBadCoefficient, std::move(context));
}

[[noreturn]] void throwInvalid(std::string owner, std::string detail) {
  MhdFailureContext context{};
  context.owner = std::move(owner);
  context.detail = std::move(detail);
  throw MhdError(MhdErrorCode::kInvalidField, std::move(context));
}

[[nodiscard]] double spacingAlong(const MhdPatchGeometry& geometry, std::size_t axis) noexcept {
  return axis == 0U ? geometry.dx_comoving
                    : (axis == 1U ? geometry.dy_comoving : geometry.dz_comoving);
}

[[nodiscard]] std::size_t edgeIndex(const MhdPatchGeometry& geometry, std::size_t axis,
                                    const MhdIndex3& e) noexcept {
  const auto i = static_cast<std::size_t>(e[0]);
  const auto j = static_cast<std::size_t>(e[1]);
  const auto k = static_cast<std::size_t>(e[2]);
  return axis == 0U ? geometry.edgeIndexX(i, j, k)
                    : (axis == 1U ? geometry.edgeIndexY(i, j, k) : geometry.edgeIndexZ(i, j, k));
}

[[nodiscard]] std::span<const double> emfLane(const MhdEdgeEmfStorage& emf,
                                              std::size_t axis) noexcept {
  return axis == 0U ? emf.emfX() : (axis == 1U ? emf.emfY() : emf.emfZ());
}

[[nodiscard]] std::span<double> emfLane(MhdEdgeEmfStorage& emf, std::size_t axis) noexcept {
  return axis == 0U ? emf.emfX() : (axis == 1U ? emf.emfY() : emf.emfZ());
}

}  // namespace

std::string_view mhdCoefficientProviderName(MhdCoefficientProviderKind kind) noexcept {
  switch (kind) {
    case MhdCoefficientProviderKind::kNone: return "none";
    case MhdCoefficientProviderKind::kConstant: return "constant";
    case MhdCoefficientProviderKind::kTable: return "table";
    case MhdCoefficientProviderKind::kChemistry: return "chemistry";
  }
  return "unknown";
}

double physicalDiffusivityToComovingCode(double eta_phys_cgs,
                                         double scale_factor,
                                         const core::UnitSystem& units) {
  if (!std::isfinite(eta_phys_cgs) || !(scale_factor > 0.0) || !std::isfinite(scale_factor)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double length_cm_per_code = units.length_si_per_code * 100.0;
  const double time_s_per_code = units.timeSiPerCode();
  if (!(length_cm_per_code > 0.0) || !(time_s_per_code > 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Physical diffusivity in code units, then the comoving operator factor
  // 1/a^2 (section 4.3). This is the ONLY place either factor is applied.
  const double code_per_cgs = time_s_per_code / (length_cm_per_code * length_cm_per_code);
  return eta_phys_cgs * code_per_cgs / (scale_factor * scale_factor);
}

// ---------------------------------------------------------------------------
// Constant provider
// ---------------------------------------------------------------------------

ConstantMhdCoefficientProvider::ConstantMhdCoefficientProvider(double eta_ohmic_cgs,
                                                               double eta_hall_cgs,
                                                               double eta_ambipolar_cgs)
    : m_eta_ohmic_cgs(eta_ohmic_cgs),
      m_eta_hall_cgs(eta_hall_cgs),
      m_eta_ambipolar_cgs(eta_ambipolar_cgs) {}

std::string ConstantMhdCoefficientProvider::provenanceId() const {
  std::ostringstream stream;
  stream.precision(17);
  stream << "constant_eta_cgs[O=" << m_eta_ohmic_cgs << ",H=" << m_eta_hall_cgs
         << ",A=" << m_eta_ambipolar_cgs << ']';
  return stream.str();
}

void ConstantMhdCoefficientProvider::validateAdmission() const {
  const auto check = [](double value, const char* name) {
    if (!std::isfinite(value)) {
      throwCoefficient("mhd_coefficient_provider",
                       std::string("non-finite constant coefficient ") + name);
    }
  };
  check(m_eta_ohmic_cgs, "eta_ohmic");
  check(m_eta_hall_cgs, "eta_hall");
  check(m_eta_ambipolar_cgs, "eta_ambipolar");
  if (m_eta_ohmic_cgs < 0.0 || m_eta_ambipolar_cgs < 0.0) {
    throwCoefficient("mhd_coefficient_provider",
                     "Ohmic and ambipolar diffusivities must be nonnegative");
  }
}

MhdTransportCoefficients ConstantMhdCoefficientProvider::evaluate(
    std::size_t cell_index,
    const MhdCellPhysicalProperties& properties) const {
  (void)cell_index;
  MhdTransportCoefficients coefficients{};
  coefficients.eta_ohmic_phys_cgs = m_eta_ohmic_cgs;
  coefficients.eta_hall_phys_cgs = m_eta_hall_cgs;
  coefficients.eta_ambipolar_phys_cgs = m_eta_ambipolar_cgs;
  coefficients.valid = std::isfinite(m_eta_ohmic_cgs) && std::isfinite(m_eta_hall_cgs) &&
      std::isfinite(m_eta_ambipolar_cgs) && m_eta_ohmic_cgs >= 0.0 &&
      m_eta_ambipolar_cgs >= 0.0;
  coefficients.closure_valid = properties.closure_valid;
  return coefficients;
}

// ---------------------------------------------------------------------------
// Chemistry provider: fail closed until a charged-species contract exists.
// ---------------------------------------------------------------------------

std::string ChemistryMhdCoefficientProvider::provenanceId() const {
  return "chemistry_transport_unavailable";
}

void ChemistryMhdCoefficientProvider::validateAdmission() const {
  MhdFailureContext context{};
  context.owner = "mhd_coefficient_provider";
  context.detail =
      "chemistry-sourced Ohmic/Hall/ambipolar transport requires a live charged-species, "
      "collision-rate, charge-neutrality and grain model; the repository does not "
      "establish that validated contract (amendment R5)";
  throw MhdError(MhdErrorCode::kUnsupportedPhysics, std::move(context));
}

MhdTransportCoefficients ChemistryMhdCoefficientProvider::evaluate(
    std::size_t cell_index,
    const MhdCellPhysicalProperties& properties) const {
  (void)cell_index;
  (void)properties;
  MhdTransportCoefficients coefficients{};
  coefficients.valid = false;
  coefficients.closure_valid = false;
  return coefficients;
}

// ---------------------------------------------------------------------------
// Versioned table provider
//
// HDF5 layout (contract chui_mhd_transport_table_v1; root attribute
// ContractId = "chui_mhd_transport_table_v1"):
//   /Axes/LogDensityCgs          float64[nrho]   attr Units="log10(g/cm^3)"
//   /Axes/LogTemperatureK        float64[nT]     attr Units="log10(K)"
//   /Axes/LogMagneticFieldGauss  float64[nB]     attr Units="log10(G)"
//   /Axes/CompositionIds         string[ncomp]
//   /Axes/GrainModelIds          string[ngrain]
//   /EtaOhmicCgs, /EtaHallCgs, /EtaAmbipolarCgs
//                                float64[nrho][nT][nB][ncomp][ngrain]
//                                attr Units="cm^2/s"
//   /BZeroLimit/EtaOhmicCgs, /BZeroLimit/EtaHallCgs, /BZeroLimit/EtaAmbipolarCgs
//                                float64[nrho][nT][ncomp][ngrain]
//                                attr Units="cm^2/s"
// The BZeroLimit group is the declared exactly-zero-field limit (section 8.1);
// it is used only for b_phys_gauss == 0 and never extrapolated from log B.
// ---------------------------------------------------------------------------

struct TableMhdCoefficientProvider::Impl {
  std::vector<double> log_rho;
  std::vector<double> log_temperature;
  std::vector<double> log_b;
  std::vector<std::string> composition_ids;
  std::vector<std::string> grain_model_ids;
  std::vector<double> eta_ohmic;
  std::vector<double> eta_hall;
  std::vector<double> eta_ambipolar;
  std::vector<double> zero_ohmic;
  std::vector<double> zero_hall;
  std::vector<double> zero_ambipolar;
  std::size_t composition_index = 0;
  std::size_t grain_index = 0;
  bool admitted = false;
  std::string admission_failure;
};

namespace {

[[nodiscard]] bool strictlyIncreasingFinite(const std::vector<double>& axis) {
  if (axis.size() < 2U) {
    return false;
  }
  for (std::size_t index = 0; index < axis.size(); ++index) {
    if (!std::isfinite(axis[index])) {
      return false;
    }
    if (index > 0U && !(axis[index] > axis[index - 1U])) {
      return false;
    }
  }
  return true;
}

// Index of the lower bracketing node for value in [axis.front(), axis.back()].
[[nodiscard]] std::size_t bracketIndex(const std::vector<double>& axis, double value) {
  const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
  std::size_t index = static_cast<std::size_t>(upper - axis.begin());
  index = index == 0U ? 0U : index - 1U;
  return std::min(index, axis.size() - 2U);
}

[[nodiscard]] double bracketWeight(const std::vector<double>& axis, std::size_t lower,
                                   double value) {
  return (value - axis[lower]) / (axis[lower + 1U] - axis[lower]);
}

#if COSMOSIM_ENABLE_HDF5

[[nodiscard]] std::vector<hsize_t> datasetDims(hid_t dataset) {
  const hid_t space = H5Dget_space(dataset);
  const int rank = H5Sget_simple_extent_ndims(space);
  std::vector<hsize_t> dims(static_cast<std::size_t>(std::max(rank, 0)));
  if (rank > 0) {
    H5Sget_simple_extent_dims(space, dims.data(), nullptr);
  }
  H5Sclose(space);
  return dims;
}

[[nodiscard]] std::string readStringAttribute(hid_t object, const char* name) {
  if (H5Aexists(object, name) <= 0) {
    return {};
  }
  const hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
  if (attribute < 0) {
    return {};
  }
  const hid_t type = H5Aget_type(attribute);
  std::string value;
  if (H5Tget_class(type) == H5T_STRING) {
    if (H5Tis_variable_str(type) > 0) {
      char* buffer = nullptr;
      const hid_t memory_type = H5Tcopy(H5T_C_S1);
      H5Tset_size(memory_type, H5T_VARIABLE);
      if (H5Aread(attribute, memory_type, &buffer) >= 0 && buffer != nullptr) {
        value = buffer;
        H5free_memory(buffer);
      }
      H5Tclose(memory_type);
    } else {
      const std::size_t size = H5Tget_size(type);
      std::vector<char> buffer(size + 1U, '\0');
      if (H5Aread(attribute, type, buffer.data()) >= 0) {
        value.assign(buffer.data(), strnlen(buffer.data(), size));
      }
    }
  }
  H5Tclose(type);
  H5Aclose(attribute);
  return value;
}

[[nodiscard]] std::vector<double> readDoubleDataset(hid_t file, const char* name,
                                                    const std::vector<hsize_t>& expected_dims,
                                                    const char* expected_units) {
  if (H5Lexists(file, name, H5P_DEFAULT) <= 0) {
    throwCoefficient("mhd_table_provider", std::string("missing dataset ") + name);
  }
  const hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  if (dataset < 0) {
    throwCoefficient("mhd_table_provider", std::string("cannot open dataset ") + name);
  }
  const std::vector<hsize_t> dims = datasetDims(dataset);
  if (!expected_dims.empty() && dims != expected_dims) {
    H5Dclose(dataset);
    throwCoefficient("mhd_table_provider",
                     std::string("dataset ") + name + " has a shape inconsistent with its axes");
  }
  if (expected_units != nullptr) {
    const std::string units = readStringAttribute(dataset, "Units");
    if (units != expected_units) {
      H5Dclose(dataset);
      throwCoefficient("mhd_table_provider", std::string("dataset ") + name +
                                                 " must declare Units=\"" + expected_units +
                                                 "\" (found \"" + units + "\")");
    }
  }
  std::size_t total = 1U;
  for (const hsize_t dim : dims) {
    total *= static_cast<std::size_t>(dim);
  }
  std::vector<double> values(total, 0.0);
  if (total > 0U && H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            values.data()) < 0) {
    H5Dclose(dataset);
    throwCoefficient("mhd_table_provider", std::string("cannot read dataset ") + name);
  }
  H5Dclose(dataset);
  return values;
}

[[nodiscard]] std::vector<std::string> readStringDataset(hid_t file, const char* name) {
  if (H5Lexists(file, name, H5P_DEFAULT) <= 0) {
    throwCoefficient("mhd_table_provider", std::string("missing dataset ") + name);
  }
  const hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  const hid_t type = H5Dget_type(dataset);
  const hid_t space = H5Dget_space(dataset);
  const hssize_t count = H5Sget_simple_extent_npoints(space);
  std::vector<std::string> values;
  if (H5Tget_class(type) != H5T_STRING || count <= 0) {
    H5Sclose(space);
    H5Tclose(type);
    H5Dclose(dataset);
    throwCoefficient("mhd_table_provider",
                     std::string("dataset ") + name + " must be a non-empty string array");
  }
  if (H5Tis_variable_str(type) > 0) {
    std::vector<char*> buffer(static_cast<std::size_t>(count), nullptr);
    const hid_t memory_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(memory_type, H5T_VARIABLE);
    const herr_t status = H5Dread(dataset, memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                  buffer.data());
    for (char* entry : buffer) {
      values.emplace_back(entry != nullptr ? entry : "");
    }
    if (status >= 0) {
#if H5_VERSION_GE(1, 12, 0)
      H5Treclaim(memory_type, space, H5P_DEFAULT, buffer.data());
#else
      H5Dvlen_reclaim(memory_type, space, H5P_DEFAULT, buffer.data());
#endif
    }
    H5Tclose(memory_type);
    if (status < 0) {
      values.clear();
    }
  } else {
    const std::size_t size = H5Tget_size(type);
    std::vector<char> buffer(static_cast<std::size_t>(count) * size, '\0');
    if (H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) >= 0) {
      for (hssize_t index = 0; index < count; ++index) {
        const char* begin = buffer.data() + static_cast<std::size_t>(index) * size;
        values.emplace_back(begin, strnlen(begin, size));
      }
    }
  }
  H5Sclose(space);
  H5Tclose(type);
  H5Dclose(dataset);
  if (values.size() != static_cast<std::size_t>(count)) {
    throwCoefficient("mhd_table_provider", std::string("cannot read dataset ") + name);
  }
  return values;
}

#endif  // COSMOSIM_ENABLE_HDF5

}  // namespace

TableMhdCoefficientProvider::TableMhdCoefficientProvider(std::filesystem::path table_path,
                                                         Options options)
    : m_impl(std::make_unique<Impl>()), m_path(std::move(table_path)), m_options(std::move(options)) {
#if COSMOSIM_ENABLE_HDF5
  if (!std::filesystem::is_regular_file(m_path)) {
    throwCoefficient("mhd_table_provider",
                     "coefficient table does not exist: " + m_path.string());
  }
  // Digest the exact bytes that are then parsed (provenance and restart pin).
  m_sha256 = internal::sha256FileHex(m_path);
  const hid_t file = H5Fopen(m_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    throwCoefficient("mhd_table_provider", "cannot open coefficient table " + m_path.string());
  }
  try {
    const hid_t root = H5Gopen2(file, "/", H5P_DEFAULT);
    const std::string contract = readStringAttribute(root, "ContractId");
    H5Gclose(root);
    if (contract != "chui_mhd_transport_table_v1") {
      throwCoefficient("mhd_table_provider",
                       "table root attribute ContractId must be chui_mhd_transport_table_v1 "
                       "(found \"" + contract + "\")");
    }
    m_impl->log_rho = readDoubleDataset(file, "/Axes/LogDensityCgs", {}, "log10(g/cm^3)");
    m_impl->log_temperature = readDoubleDataset(file, "/Axes/LogTemperatureK", {}, "log10(K)");
    m_impl->log_b = readDoubleDataset(file, "/Axes/LogMagneticFieldGauss", {}, "log10(G)");
    m_impl->composition_ids = readStringDataset(file, "/Axes/CompositionIds");
    m_impl->grain_model_ids = readStringDataset(file, "/Axes/GrainModelIds");
    const std::vector<hsize_t> full_dims{
        m_impl->log_rho.size(), m_impl->log_temperature.size(), m_impl->log_b.size(),
        m_impl->composition_ids.size(), m_impl->grain_model_ids.size()};
    const std::vector<hsize_t> zero_dims{m_impl->log_rho.size(), m_impl->log_temperature.size(),
                                         m_impl->composition_ids.size(),
                                         m_impl->grain_model_ids.size()};
    m_impl->eta_ohmic = readDoubleDataset(file, "/EtaOhmicCgs", full_dims, "cm^2/s");
    m_impl->eta_hall = readDoubleDataset(file, "/EtaHallCgs", full_dims, "cm^2/s");
    m_impl->eta_ambipolar = readDoubleDataset(file, "/EtaAmbipolarCgs", full_dims, "cm^2/s");
    m_impl->zero_ohmic = readDoubleDataset(file, "/BZeroLimit/EtaOhmicCgs", zero_dims, "cm^2/s");
    m_impl->zero_hall = readDoubleDataset(file, "/BZeroLimit/EtaHallCgs", zero_dims, "cm^2/s");
    m_impl->zero_ambipolar =
        readDoubleDataset(file, "/BZeroLimit/EtaAmbipolarCgs", zero_dims, "cm^2/s");
  } catch (...) {
    H5Fclose(file);
    throw;
  }
  H5Fclose(file);
#else
  throwCoefficient("mhd_table_provider",
                   "table coefficient provider requires a COSMOSIM_ENABLE_HDF5=ON build");
#endif
}

TableMhdCoefficientProvider::~TableMhdCoefficientProvider() = default;

std::string TableMhdCoefficientProvider::provenanceId() const {
  std::ostringstream stream;
  stream << "table[sha256=" << m_sha256 << ";composition=" << m_options.composition_id
         << ";grain=" << m_options.grain_model_id << ']';
  return stream.str();
}

void TableMhdCoefficientProvider::validateAdmission() const {
#if COSMOSIM_ENABLE_HDF5
  Impl& impl = *m_impl;
  if (!strictlyIncreasingFinite(impl.log_rho) || !strictlyIncreasingFinite(impl.log_temperature) ||
      !strictlyIncreasingFinite(impl.log_b)) {
    throwCoefficient("mhd_table_provider",
                     "continuous table axes must be finite, strictly increasing and have at "
                     "least two nodes");
  }
  const auto unique = [](const std::vector<std::string>& ids) {
    return std::set<std::string>(ids.begin(), ids.end()).size() == ids.size();
  };
  if (!unique(impl.composition_ids) || !unique(impl.grain_model_ids)) {
    throwCoefficient("mhd_table_provider", "categorical table ids must be unique");
  }
  const auto locate = [](const std::vector<std::string>& ids, const std::string& id,
                         const char* what) {
    const auto found = std::find(ids.begin(), ids.end(), id);
    if (id.empty() || found == ids.end()) {
      throwCoefficient("mhd_table_provider",
                       std::string(what) + " id \"" + id + "\" is not present in the table");
    }
    return static_cast<std::size_t>(found - ids.begin());
  };
  impl.composition_index = locate(impl.composition_ids, m_options.composition_id, "composition");
  impl.grain_index = locate(impl.grain_model_ids, m_options.grain_model_id, "grain model");
  const auto nonnegative = [](const std::vector<double>& data) {
    return std::all_of(data.begin(), data.end(),
                       [](double v) { return std::isfinite(v) && v >= 0.0; });
  };
  const auto finite = [](const std::vector<double>& data) {
    return std::all_of(data.begin(), data.end(), [](double v) { return std::isfinite(v); });
  };
  if (!nonnegative(impl.eta_ohmic) || !nonnegative(impl.eta_ambipolar) ||
      !nonnegative(impl.zero_ohmic) || !nonnegative(impl.zero_ambipolar)) {
    throwCoefficient("mhd_table_provider",
                     "Ohmic/ambipolar table entries must be finite and nonnegative");
  }
  if (!finite(impl.eta_hall) || !finite(impl.zero_hall)) {
    throwCoefficient("mhd_table_provider", "Hall table entries must be finite");
  }
  if (!m_options.expected_sha256.empty() && m_options.expected_sha256 != m_sha256) {
    throwCoefficient("mhd_table_provider",
                     "coefficient table SHA-256 " + m_sha256 +
                         " does not match the pinned digest " + m_options.expected_sha256);
  }
  impl.admitted = true;
#else
  throwCoefficient("mhd_table_provider",
                   "table provider requires a COSMOSIM_ENABLE_HDF5=ON build");
#endif
}

MhdTransportCoefficients TableMhdCoefficientProvider::evaluate(
    std::size_t cell_index,
    const MhdCellPhysicalProperties& properties) const {
  (void)cell_index;
  MhdTransportCoefficients coefficients{};
  coefficients.closure_valid = properties.closure_valid;
  const Impl& impl = *m_impl;
  if (!impl.admitted) {
    coefficients.valid = false;
    return coefficients;
  }
  if (!(properties.rho_phys_cgs > 0.0) || !(properties.temperature_k > 0.0) ||
      !std::isfinite(properties.rho_phys_cgs) || !std::isfinite(properties.temperature_k) ||
      !std::isfinite(properties.b_phys_gauss) || properties.b_phys_gauss < 0.0) {
    coefficients.valid = false;
    return coefficients;
  }
  const double log_rho = std::log10(properties.rho_phys_cgs);
  const double log_temperature = std::log10(properties.temperature_k);
  const bool zero_field = properties.b_phys_gauss == 0.0;
  const double log_b = zero_field ? 0.0 : std::log10(properties.b_phys_gauss);
  const auto inside = [](const std::vector<double>& axis, double value) {
    return value >= axis.front() && value <= axis.back();
  };
  // Strict coverage: outside the table is fatal for the caller; no
  // extrapolation, clamping or replacement.
  if (!inside(impl.log_rho, log_rho) || !inside(impl.log_temperature, log_temperature) ||
      (!zero_field && !inside(impl.log_b, log_b))) {
    coefficients.valid = false;
    return coefficients;
  }
  const std::size_t i = bracketIndex(impl.log_rho, log_rho);
  const std::size_t j = bracketIndex(impl.log_temperature, log_temperature);
  const double ti = bracketWeight(impl.log_rho, i, log_rho);
  const double tj = bracketWeight(impl.log_temperature, j, log_temperature);
  const std::size_t n_t = impl.log_temperature.size();
  const std::size_t n_b = impl.log_b.size();
  const std::size_t n_comp = impl.composition_ids.size();
  const std::size_t n_grain = impl.grain_model_ids.size();
  const std::size_t category = impl.composition_index * n_grain + impl.grain_index;

  std::array<std::size_t, 8> corner{};
  std::array<double, 8> weight{};
  std::size_t corners = 0;
  if (zero_field) {
    for (std::size_t di = 0; di < 2U; ++di) {
      for (std::size_t dj = 0; dj < 2U; ++dj) {
        corner[corners] = ((i + di) * n_t + (j + dj)) * n_comp * n_grain + category;
        weight[corners] = (di == 0U ? 1.0 - ti : ti) * (dj == 0U ? 1.0 - tj : tj);
        ++corners;
      }
    }
  } else {
    const std::size_t k = bracketIndex(impl.log_b, log_b);
    const double tk = bracketWeight(impl.log_b, k, log_b);
    for (std::size_t di = 0; di < 2U; ++di) {
      for (std::size_t dj = 0; dj < 2U; ++dj) {
        for (std::size_t dk = 0; dk < 2U; ++dk) {
          corner[corners] =
              (((i + di) * n_t + (j + dj)) * n_b + (k + dk)) * n_comp * n_grain + category;
          weight[corners] = (di == 0U ? 1.0 - ti : ti) * (dj == 0U ? 1.0 - tj : tj) *
              (dk == 0U ? 1.0 - tk : tk);
          ++corners;
        }
      }
    }
  }
  // Ohmic/ambipolar: monotone log-space interpolation when every bracketing
  // node is strictly positive; linear (still monotone, nonnegative) when a
  // node is exactly zero. Hall: sign-aware linear interpolation, so a zero
  // crossing is reproduced exactly where the linear interpolant vanishes.
  const auto positive = [&](const std::vector<double>& data) {
    bool all_positive = true;
    for (std::size_t c = 0; c < corners; ++c) {
      all_positive = all_positive && data[corner[c]] > 0.0;
    }
    double value = 0.0;
    if (all_positive) {
      for (std::size_t c = 0; c < corners; ++c) {
        value += weight[c] * std::log(data[corner[c]]);
      }
      return std::exp(value);
    }
    for (std::size_t c = 0; c < corners; ++c) {
      value += weight[c] * data[corner[c]];
    }
    return value;
  };
  const auto signed_linear = [&](const std::vector<double>& data) {
    double value = 0.0;
    for (std::size_t c = 0; c < corners; ++c) {
      value += weight[c] * data[corner[c]];
    }
    return value;
  };
  coefficients.eta_ohmic_phys_cgs = positive(zero_field ? impl.zero_ohmic : impl.eta_ohmic);
  coefficients.eta_ambipolar_phys_cgs =
      positive(zero_field ? impl.zero_ambipolar : impl.eta_ambipolar);
  coefficients.eta_hall_phys_cgs = signed_linear(zero_field ? impl.zero_hall : impl.eta_hall);
  coefficients.valid = std::isfinite(coefficients.eta_ohmic_phys_cgs) &&
      std::isfinite(coefficients.eta_hall_phys_cgs) &&
      std::isfinite(coefficients.eta_ambipolar_phys_cgs) &&
      coefficients.eta_ohmic_phys_cgs >= 0.0 && coefficients.eta_ambipolar_phys_cgs >= 0.0;
  return coefficients;
}

// ---------------------------------------------------------------------------
// RKL2 coefficients
// ---------------------------------------------------------------------------

MhdRkl2Coefficients makeRkl2Coefficients(std::int32_t stages) {
  if (stages < 2) {
    throwInvalid("mhd_rkl2", "RKL2 requires at least two stages");
  }
  MhdRkl2Coefficients coefficients{};
  const auto count = static_cast<std::size_t>(stages) + 1U;
  const double s = static_cast<double>(stages);
  coefficients.stages = stages;
  coefficients.w1 = 4.0 / (s * s + s - 2.0);
  coefficients.b.assign(count, 0.0);
  coefficients.mu.assign(count, 0.0);
  coefficients.nu.assign(count, 0.0);
  coefficients.mu_tilde.assign(count, 0.0);
  coefficients.gamma_tilde.assign(count, 0.0);
  coefficients.b[0] = 1.0 / 3.0;
  coefficients.b[1] = 1.0 / 3.0;
  for (std::size_t j = 2; j < count; ++j) {
    const double jd = static_cast<double>(j);
    coefficients.b[j] = (jd * jd + jd - 2.0) / (2.0 * jd * (jd + 1.0));
  }
  coefficients.mu_tilde[1] = coefficients.b[1] * coefficients.w1;
  for (std::size_t j = 2; j < count; ++j) {
    const double jd = static_cast<double>(j);
    coefficients.mu[j] = (2.0 * jd - 1.0) / jd * coefficients.b[j] / coefficients.b[j - 1U];
    coefficients.nu[j] = -(jd - 1.0) / jd * coefficients.b[j] / coefficients.b[j - 2U];
    coefficients.mu_tilde[j] = coefficients.mu[j] * coefficients.w1;
    coefficients.gamma_tilde[j] = -(1.0 - coefficients.b[j - 1U]) * coefficients.mu_tilde[j];
  }
  return coefficients;
}

double rkl2StabilityFactor(std::int32_t stages) noexcept {
  if (stages <= 1) {
    return 1.0;
  }
  const double s = static_cast<double>(stages);
  return (s * s + s - 2.0) / 4.0;
}

std::int32_t selectRkl2StageCount(double ratio, std::int32_t max_stages) noexcept {
  if (!(ratio > 1.0)) {
    return 0;
  }
  if (!std::isfinite(ratio) || max_stages < 2) {
    return max_stages + 1;
  }
  // Least s with (s^2 + s - 2)/4 >= ratio: s >= (sqrt(9 + 16 ratio) - 1)/2.
  const double estimate = std::ceil((std::sqrt(9.0 + 16.0 * ratio) - 1.0) / 2.0);
  if (!(estimate <= static_cast<double>(max_stages))) {
    return max_stages + 1;
  }
  auto stages = static_cast<std::int32_t>(std::max(2.0, estimate));
  while (stages > 2 && rkl2StabilityFactor(stages - 1) >= ratio) {
    --stages;
  }
  while (rkl2StabilityFactor(stages) < ratio) {
    ++stages;
  }
  return stages <= max_stages ? stages : max_stages + 1;
}

MhdRkl2Plan planRkl2(double ratio, std::int32_t max_stages_per_cycle) {
  if (!std::isfinite(ratio) || ratio < 0.0 || max_stages_per_cycle < 2) {
    throwInvalid("mhd_rkl2", "invalid RKL2 planning request");
  }
  MhdRkl2Plan plan{};
  // A single cycle with s = max covers ratio_max; split into the least number
  // of equal cycles so each cycle's ratio is covered by <= max stages.
  const double ratio_max = rkl2StabilityFactor(max_stages_per_cycle);
  const double cycles = std::max(1.0, std::ceil(ratio / ratio_max));
  if (cycles > static_cast<double>(std::numeric_limits<std::int32_t>::max())) {
    MhdFailureContext context{};
    context.owner = "mhd_rkl2";
    context.detail = "diffusion interval requires an unrepresentable number of RKL2 cycles";
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(context));
  }
  plan.cycles = static_cast<std::int32_t>(cycles);
  const double per_cycle = ratio / cycles;
  plan.stages_per_cycle = std::max<std::int32_t>(2, selectRkl2StageCount(per_cycle,
                                                                          max_stages_per_cycle));
  return plan;
}

MhdTransportStageDiagnostics& MhdTransportStageDiagnostics::operator+=(
    const MhdTransportStageDiagnostics& rhs) noexcept {
  explicit_diffusion_substeps += rhs.explicit_diffusion_substeps;
  rkl_cycle_count += rhs.rkl_cycle_count;
  rkl_stage_count += rhs.rkl_stage_count;
  hall_subcycle_count += rhs.hall_subcycle_count;
  hall_rejected_subcycle_count += rhs.hall_rejected_subcycle_count;
  coefficient_evaluations += rhs.coefficient_evaluations;
  joule_heating_integral += rhs.joule_heating_integral;
  hall_work_integral += rhs.hall_work_integral;
  max_current_density = std::max(max_current_density, rhs.max_current_density);
  const auto merge_min = [](double lhs, double value) {
    if (!(value > 0.0)) {
      return lhs;
    }
    return lhs > 0.0 ? std::min(lhs, value) : value;
  };
  min_diffusion_dt_code = merge_min(min_diffusion_dt_code, rhs.min_diffusion_dt_code);
  min_hall_dt_code = merge_min(min_hall_dt_code, rhs.min_hall_dt_code);
  return *this;
}

// ---------------------------------------------------------------------------
// Physical property model
// ---------------------------------------------------------------------------

void MhdPhysicalPropertyModel::validate() const {
  if (units == nullptr) {
    throwCoefficient("mhd_property_model", "no unit system bound to the MHD property model");
  }
  units->validate();
  if (!(mean_molecular_weight > 0.0) || !std::isfinite(mean_molecular_weight)) {
    throwCoefficient("mhd_property_model",
                     "mean molecular weight must be bound from the composition closure");
  }
}

double MhdPhysicalPropertyModel::pressureUnitCgs() const {
  // P_unit = M_unit / (L_unit T_unit^2) = M_unit V_unit^2 / L_unit^3.
  const double mass_g = units->mass_si_per_code * 1.0e3;
  const double length_cm = units->length_si_per_code * 1.0e2;
  const double velocity_cm_per_s = units->velocity_si_per_code * 1.0e2;
  return mass_g * velocity_cm_per_s * velocity_cm_per_s / (length_cm * length_cm * length_cm);
}

double MhdPhysicalPropertyModel::fieldCodeToGauss(double b_physical_code) const {
  // Normalized Lorentz-Heaviside code field: B_code^2 is an energy density in
  // code pressure units, so B_LH,cgs = B_code sqrt(P_unit_cgs); Gaussian
  // B_G = sqrt(4 pi) B_LH.
  return b_physical_code * std::sqrt(4.0 * std::numbers::pi * pressureUnitCgs());
}

MhdCellPhysicalProperties MhdPhysicalPropertyModel::evaluate(double rho_comoving_code,
                                                             double pressure_comoving_code,
                                                             double b_comoving_magnitude_code,
                                                             double scale_factor) const {
  MhdCellPhysicalProperties properties{};
  const double a3 = scale_factor * scale_factor * scale_factor;
  const double rho_phys_code = rho_comoving_code / a3;
  const double pressure_phys_code = pressure_comoving_code / a3;
  const double b_phys_code = b_comoving_magnitude_code / (scale_factor * scale_factor);
  properties.rho_phys_cgs = units->densityCodeToCgs(rho_phys_code);
  const double pressure_cgs = pressure_phys_code * pressureUnitCgs();
  properties.temperature_k = mean_molecular_weight * k_proton_mass_g * pressure_cgs /
      (k_boltzmann_erg_per_k * properties.rho_phys_cgs);
  properties.b_phys_gauss = fieldCodeToGauss(b_phys_code);
  properties.closure_valid = std::isfinite(properties.rho_phys_cgs) &&
      properties.rho_phys_cgs > 0.0 && std::isfinite(properties.temperature_k) &&
      properties.temperature_k > 0.0 && std::isfinite(properties.b_phys_gauss);
  return properties;
}

// ---------------------------------------------------------------------------
// Workspaces
// ---------------------------------------------------------------------------

void MhdBoundaryFaceEnergyFluxes::resize(const MhdPatchGeometry& geometry) {
  nx = geometry.nx;
  ny = geometry.ny;
  nz = geometry.nz;
  const std::array<std::size_t, 3> n{nx, ny, nz};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::size_t count = n[(axis + 1U) % 3U] * n[(axis + 2U) % 3U];
    for (std::size_t side = 0; side < 2U; ++side) {
      if (values[axis][side].size() != count) {
        values[axis][side].resize(count, 0.0);
      }
    }
  }
}

double& MhdBoundaryFaceEnergyFluxes::at(std::size_t axis, std::size_t side, std::size_t t1,
                                        std::size_t t2) noexcept {
  const std::array<std::size_t, 3> n{nx, ny, nz};
  return values[axis][side][t1 + n[(axis + 1U) % 3U] * t2];
}

double MhdBoundaryFaceEnergyFluxes::at(std::size_t axis, std::size_t side, std::size_t t1,
                                       std::size_t t2) const noexcept {
  const std::array<std::size_t, 3> n{nx, ny, nz};
  return values[axis][side][t1 + n[(axis + 1U) % 3U] * t2];
}

std::uint64_t MhdBoundaryFaceEnergyFluxes::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (const auto& sides : values) {
    for (const MhdDoubleArray& lane : sides) {
      bytes += lane.ownedCapacityBytes();
    }
  }
  return bytes;
}

void MhdNonIdealWorkspace::resize(const MhdPatchGeometry& geometry, bool super_stepping_lanes) {
  nx = geometry.nx;
  ny = geometry.ny;
  nz = geometry.nz;
  const std::size_t padded_cells = (nx + 2U) * (ny + 2U) * (nz + 2U);
  for (MhdDoubleArray* lane : {&eta_ohmic_code, &eta_hall_code, &eta_ambipolar_code}) {
    if (lane->size() != padded_cells) {
      lane->resize(padded_cells, 0.0);
    }
  }
  boundary_energy.resize(geometry);
  const std::size_t cells = geometry.cellCount();
  if (energy_rate_stage.size() != cells) {
    energy_rate_stage.resize(cells, 0.0);
  }
  if (super_stepping_lanes) {
    const std::array<std::size_t, 3> faces{geometry.faceCountX(), geometry.faceCountY(),
                                           geometry.faceCountZ()};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (yjm2_faces[axis].size() != faces[axis]) {
        yjm2_faces[axis].resize(faces[axis], 0.0);
      }
    }
    for (MhdDoubleArray* lane : {&yjm2_energy, &energy_rate_y0}) {
      if (lane->size() != cells) {
        lane->resize(cells, 0.0);
      }
    }
    if (emf_y0.geometry().nx != nx || emf_y0.geometry().ny != ny || emf_y0.geometry().nz != nz) {
      emf_y0.resize(geometry);
    }
  }
}

std::uint64_t MhdNonIdealWorkspace::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = eta_ohmic_code.ownedCapacityBytes() + eta_hall_code.ownedCapacityBytes() +
      eta_ambipolar_code.ownedCapacityBytes() + boundary_energy.ownedCapacityBytes() +
      yjm2_energy.ownedCapacityBytes() + emf_y0.ownedCapacityBytes() +
      energy_rate_y0.ownedCapacityBytes() + energy_rate_stage.ownedCapacityBytes();
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    bytes += yjm2_faces[axis].ownedCapacityBytes();
  }
  return bytes;
}

// ---------------------------------------------------------------------------
// Non-ideal rate operator
// ---------------------------------------------------------------------------

MhdNonIdealRateOperator::MhdNonIdealRateOperator(const MhdCoefficientProvider* provider,
                                                 MhdPhysicalPropertyModel properties,
                                                 MhdStiffnessLimits limits,
                                                 bool enable_ohmic,
                                                 bool enable_hall,
                                                 bool enable_ambipolar)
    : m_provider(provider),
      m_properties(properties),
      m_limits(limits),
      m_enable_ohmic(enable_ohmic),
      m_enable_hall(enable_hall),
      m_enable_ambipolar(enable_ambipolar) {
  if (m_provider == nullptr) {
    throwCoefficient("mhd_nonideal", "non-ideal rate operator requires a coefficient provider");
  }
  m_properties.validate();
  if (!(m_limits.cfl_diffusion > 0.0) || !(m_limits.cfl_hall > 0.0) ||
      m_limits.rkl_max_stages < 2 || m_limits.rkl_switch_ratio < 1 ||
      m_limits.hall_subcycle_cost_cap < 1 || m_limits.rkl_stage_cost_cap < 2) {
    throwInvalid("mhd_nonideal", "invalid non-ideal stiffness limits");
  }
}

void MhdNonIdealRateOperator::evaluateCoefficients(const MhdPatchGeometry& geometry,
                                                   const MhdPrimitivePatch& primitives,
                                                   double scale_factor,
                                                   MhdNonIdealWorkspace& workspace,
                                                   MhdTransportStageDiagnostics& diagnostics) const {
  if (primitives.ghost < 1U || primitives.nx != geometry.nx || primitives.ny != geometry.ny ||
      primitives.nz != geometry.nz) {
    throwInvalid("mhd_nonideal", "primitive cache does not match the patch");
  }
  workspace.resize(geometry, false);
  const double sqrt_a = std::sqrt(scale_factor);
  const auto nx = static_cast<std::ptrdiff_t>(geometry.nx);
  const auto ny = static_cast<std::ptrdiff_t>(geometry.ny);
  const auto nz = static_cast<std::ptrdiff_t>(geometry.nz);
  for (std::ptrdiff_t k = -1; k <= nz; ++k) {
    for (std::ptrdiff_t j = -1; j <= ny; ++j) {
      for (std::ptrdiff_t i = -1; i <= nx; ++i) {
        const MhdPrimitiveState state = primitives.load({i, j, k});
        const double b_comoving = sqrt_a *
            std::sqrt(state.bx_comoving * state.bx_comoving +
                      state.by_comoving * state.by_comoving +
                      state.bz_comoving * state.bz_comoving);
        const MhdCellPhysicalProperties properties = m_properties.evaluate(
            state.rho_comoving, state.pressure_comoving, b_comoving, scale_factor);
        const std::size_t padded = workspace.coefficientIndex(i, j, k);
        const MhdTransportCoefficients physical = m_provider->evaluate(padded, properties);
        ++diagnostics.coefficient_evaluations;
        if (!physical.valid || !physical.closure_valid || !properties.closure_valid) {
          std::ostringstream stream;
          stream << "coefficient provider rejected cell (" << i << ',' << j << ',' << k
                 << ") rho_cgs=" << properties.rho_phys_cgs << " T=" << properties.temperature_k
                 << " B_G=" << properties.b_phys_gauss
                 << " (invalid, closure-invalid or outside table coverage)";
          throwCoefficient("mhd_nonideal", stream.str());
        }
        if (physical.eta_ohmic_phys_cgs < 0.0 || physical.eta_ambipolar_phys_cgs < 0.0) {
          throwCoefficient("mhd_nonideal", "negative Ohmic or ambipolar diffusivity");
        }
        const auto convert = [&](double value, bool enabled) {
          if (!enabled) {
            return 0.0;
          }
          const double code = physicalDiffusivityToComovingCode(value, scale_factor,
                                                                *m_properties.units);
          if (!std::isfinite(code)) {
            throwCoefficient("mhd_nonideal",
                             "diffusivity unit conversion produced a non-finite value");
          }
          return code;
        };
        workspace.eta_ohmic_code[padded] = convert(physical.eta_ohmic_phys_cgs, m_enable_ohmic);
        workspace.eta_hall_code[padded] = convert(physical.eta_hall_phys_cgs, m_enable_hall);
        workspace.eta_ambipolar_code[padded] =
            convert(physical.eta_ambipolar_phys_cgs, m_enable_ambipolar);
      }
    }
  }
}

double MhdNonIdealRateOperator::diffusionStableTimeStep(
    const MhdPatchGeometry& geometry, const MhdNonIdealWorkspace& workspace) const {
  const double inverse_square_sum = 1.0 / (geometry.dx_comoving * geometry.dx_comoving) +
      1.0 / (geometry.dy_comoving * geometry.dy_comoving) +
      1.0 / (geometry.dz_comoving * geometry.dz_comoving);
  double maximum = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::size_t padded = workspace.coefficientIndex(
            static_cast<std::ptrdiff_t>(i), static_cast<std::ptrdiff_t>(j),
            static_cast<std::ptrdiff_t>(k));
        maximum = std::max(maximum, workspace.eta_ohmic_code[padded] +
                                        workspace.eta_ambipolar_code[padded]);
      }
    }
  }
  if (!(maximum > 0.0)) {
    return std::numeric_limits<double>::infinity();
  }
  const double dt = m_limits.cfl_diffusion / (maximum * inverse_square_sum);
  if (!std::isfinite(dt) || !(dt > 0.0)) {
    MhdFailureContext context{};
    context.owner = "mhd_diffusion";
    context.detail = "diffusion stable timestep underflowed or is non-finite";
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(context));
  }
  return dt;
}

double MhdNonIdealRateOperator::hallStableTimeStep(const MhdPatchGeometry& geometry,
                                                   const MhdNonIdealWorkspace& workspace) const {
  const double inverse_square_sum = 1.0 / (geometry.dx_comoving * geometry.dx_comoving) +
      1.0 / (geometry.dy_comoving * geometry.dy_comoving) +
      1.0 / (geometry.dz_comoving * geometry.dz_comoving);
  double maximum = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::size_t padded = workspace.coefficientIndex(
            static_cast<std::ptrdiff_t>(i), static_cast<std::ptrdiff_t>(j),
            static_cast<std::ptrdiff_t>(k));
        maximum = std::max(maximum, std::abs(workspace.eta_hall_code[padded]));
      }
    }
  }
  if (!(maximum > 0.0)) {
    return std::numeric_limits<double>::infinity();
  }
  const double diffusive_form = m_limits.cfl_hall / (maximum * inverse_square_sum);
  const double whistler_fan = m_limits.cfl_hall / (std::numbers::pi * maximum * inverse_square_sum);
  const double dt = std::min(diffusive_form, whistler_fan);
  if (!std::isfinite(dt) || !(dt > 0.0)) {
    MhdFailureContext context{};
    context.owner = "mhd_hall";
    context.detail = "Hall stable timestep underflowed or is non-finite";
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(context));
  }
  return dt;
}

void MhdNonIdealRateOperator::computeEmf(MhdNonIdealOperatorKind kind,
                                         const MhdPatchGeometry& geometry,
                                         const MhdPatchBoundary& boundary,
                                         const MhdStaggeredFieldPatch& staggered,
                                         double scale_factor,
                                         const MhdNonIdealWorkspace& workspace,
                                         MhdEdgeEmfStorage& emf,
                                         MhdTransportStageDiagnostics& diagnostics) const {
  if (staggered.tangential_ghost < 2U || staggered.nx != geometry.nx ||
      staggered.ny != geometry.ny || staggered.nz != geometry.nz) {
    throwInvalid("mhd_nonideal", "staggered field cache does not match the patch");
  }
  if (emf.geometry().nx != geometry.nx || emf.geometry().ny != geometry.ny ||
      emf.geometry().nz != geometry.nz) {
    emf.resize(geometry);
  }
  const double sqrt_a = std::sqrt(scale_factor);
  const std::array<std::ptrdiff_t, 3> n{static_cast<std::ptrdiff_t>(geometry.nx),
                                        static_cast<std::ptrdiff_t>(geometry.ny),
                                        static_cast<std::ptrdiff_t>(geometry.nz)};
  const auto field = [&](std::size_t normal, const MhdIndex3& c) {
    return sqrt_a * staggered.at(static_cast<MhdAxis>(normal), c);  // B_c
  };
  // (1) Compatible edge current J_c = d_a B_b - d_b B_a at any edge whose own
  // (c) cell index lies in -1..n_c (dual-mesh curl of the face field).
  const auto edgeCurrent = [&](std::size_t c, const MhdIndex3& e) {
    const std::size_t a = (c + 1U) % 3U;
    const std::size_t b = (c + 2U) % 3U;
    MhdIndex3 face_b_hi = e;  // b-face: [b]=node, [a]=cell, [c]=cell
    MhdIndex3 face_b_lo = e;
    face_b_lo[a] -= 1;
    MhdIndex3 face_a_hi = e;  // a-face: [a]=node, [b]=cell, [c]=cell
    MhdIndex3 face_a_lo = e;
    face_a_lo[b] -= 1;
    return (field(b, face_b_hi) - field(b, face_b_lo)) / spacingAlong(geometry, a) -
        (field(a, face_a_hi) - field(a, face_a_lo)) / spacingAlong(geometry, b);
  };
  double max_current = 0.0;

  // (2) Edge EMFs.
  for (std::size_t c = 0; c < 3U; ++c) {
    const std::size_t a = (c + 1U) % 3U;
    const std::size_t b = (c + 2U) % 3U;
    const double edge_length = spacingAlong(geometry, c);
    const double da = spacingAlong(geometry, a);
    const double db = spacingAlong(geometry, b);
    const double cell_volume = geometry.cellVolume();
    std::span<double> lane = emfLane(emf, c);
    for (std::ptrdiff_t kc = 0; kc < n[c]; ++kc) {
      for (std::ptrdiff_t nb = 0; nb <= n[b]; ++nb) {
        for (std::ptrdiff_t na = 0; na <= n[a]; ++na) {
          MhdIndex3 e{};
          e[c] = kc;
          e[a] = na;
          e[b] = nb;
          std::array<double, 3> current{};
          current[c] = edgeCurrent(c, e);
          max_current = std::max(max_current, std::abs(current[c]));
          // J_a: four a-edges (a-cell na-1|na, b-node nb, c-node kc|kc+1).
          double sum_a = 0.0;
          double sum_b = 0.0;
          for (int da_offset = -1; da_offset <= 0; ++da_offset) {
            for (int dc_offset = 0; dc_offset <= 1; ++dc_offset) {
              MhdIndex3 edge_a{};
              edge_a[a] = na + da_offset;
              edge_a[b] = nb;
              edge_a[c] = kc + dc_offset;
              sum_a += edgeCurrent(a, edge_a);
            }
          }
          for (int db_offset = -1; db_offset <= 0; ++db_offset) {
            for (int dc_offset = 0; dc_offset <= 1; ++dc_offset) {
              MhdIndex3 edge_b{};
              edge_b[b] = nb + db_offset;
              edge_b[a] = na;
              edge_b[c] = kc + dc_offset;
              sum_b += edgeCurrent(b, edge_b);
            }
          }
          current[a] = 0.25 * sum_a;
          current[b] = 0.25 * sum_b;
          // Edge-averaged field and coefficients.
          MhdIndex3 face_a_lo{};
          face_a_lo[c] = kc;
          face_a_lo[a] = na;
          face_a_lo[b] = nb - 1;
          MhdIndex3 face_a_hi = face_a_lo;
          face_a_hi[b] = nb;
          MhdIndex3 face_b_lo{};
          face_b_lo[c] = kc;
          face_b_lo[b] = nb;
          face_b_lo[a] = na - 1;
          MhdIndex3 face_b_hi = face_b_lo;
          face_b_hi[a] = na;
          std::array<double, 3> b_edge{};
          b_edge[a] = 0.5 * (field(a, face_a_lo) + field(a, face_a_hi));
          b_edge[b] = 0.5 * (field(b, face_b_lo) + field(b, face_b_hi));
          double sum_c = 0.0;
          double eta_o = 0.0;
          double eta_h = 0.0;
          double eta_a = 0.0;
          for (int qa = -1; qa <= 0; ++qa) {
            for (int qb = -1; qb <= 0; ++qb) {
              MhdIndex3 face_c{};
              face_c[a] = na + qa;
              face_c[b] = nb + qb;
              face_c[c] = kc;
              sum_c += field(c, face_c);
              face_c[c] = kc + 1;
              sum_c += field(c, face_c);
              MhdIndex3 cell{};
              cell[a] = na + qa;
              cell[b] = nb + qb;
              cell[c] = kc;
              const std::size_t padded = workspace.coefficientIndex(cell[0], cell[1], cell[2]);
              eta_o += 0.25 * workspace.eta_ohmic_code[padded];
              eta_h += 0.25 * workspace.eta_hall_code[padded];
              eta_a += 0.25 * workspace.eta_ambipolar_code[padded];
            }
          }
          b_edge[c] = 0.125 * sum_c;
          double electric = 0.0;
          double edge_weight = 1.0;  // share of this edge's dual volume on the patch
          for (const std::size_t node_axis : {a, b}) {
            if (e[node_axis] == 0 || e[node_axis] == n[node_axis]) {
              edge_weight *= 0.5;
            }
          }
          if (kind == MhdNonIdealOperatorKind::kDiffusion) {
            const std::array<double, 3> field_vector = nonIdealElectricField(
                current, b_edge, eta_o, 0.0, eta_a, k_b_regularization_code);
            electric = field_vector[c];
            diagnostics.joule_heating_integral +=
                edge_weight * electric * current[c] * cell_volume / scale_factor;
          } else {
            const std::array<double, 3> field_vector = nonIdealElectricField(
                current, b_edge, 0.0, eta_h, 0.0, k_b_regularization_code);
            electric = field_vector[c];
            diagnostics.hall_work_integral +=
                edge_weight * electric * current[c] * cell_volume / scale_factor;
            // Whistler-HLL upwinding with reconstructed staggered edge states.
            const auto toward = [&](std::size_t normal, MhdIndex3 face, std::size_t along,
                                    double direction) {
              MhdIndex3 minus = face;
              MhdIndex3 plus = face;
              minus[along] -= 1;
              plus[along] += 1;
              const double center = field(normal, face);
              return center + direction * 0.5 *
                  mcLimitedSlope(field(normal, minus), center, field(normal, plus));
            };
            const double whistler_a = std::numbers::pi * std::abs(eta_h) / da;
            const double whistler_b = std::numbers::pi * std::abs(eta_h) / db;
            const double jump_b_across_a =
                toward(b, face_b_hi, a, -1.0) - toward(b, face_b_lo, a, +1.0);
            const double jump_a_across_b =
                toward(a, face_a_hi, b, -1.0) - toward(a, face_a_lo, b, +1.0);
            electric += 0.5 * whistler_a * jump_b_across_a - 0.5 * whistler_b * jump_a_across_b;
          }
          if (!std::isfinite(electric)) {
            throwInvalid("mhd_nonideal", "non-finite non-ideal edge electric field");
          }
          lane[edgeIndex(geometry, c, e)] = electric * edge_length;
        }
      }
    }
  }
  diagnostics.max_current_density = std::max(diagnostics.max_current_density, max_current);
  (void)MhdEmfBuilder::applyEdgeBoundaryRules(geometry, boundary, emf);
}

void MhdNonIdealRateOperator::applyEnergyFlux(const MhdPatchGeometry& geometry,
                                              const MhdPrimitivePatch& primitives,
                                              const MhdEdgeEmfStorage& emf,
                                              double scale_factor,
                                              double weight,
                                              MhdNonIdealWorkspace& workspace,
                                              std::span<double> energy_target) const {
  if (energy_target.size() != geometry.cellCount()) {
    throwInvalid("mhd_nonideal", "energy target does not match the patch");
  }
  workspace.boundary_energy.resize(geometry);
  const double sqrt_a = std::sqrt(scale_factor);
  const std::array<std::ptrdiff_t, 3> n{static_cast<std::ptrdiff_t>(geometry.nx),
                                        static_cast<std::ptrdiff_t>(geometry.ny),
                                        static_cast<std::ptrdiff_t>(geometry.nz)};
  const auto cellIndex = [&](const MhdIndex3& c) {
    return geometry.cellIndex(static_cast<std::size_t>(c[0]), static_cast<std::size_t>(c[1]),
                              static_cast<std::size_t>(c[2]));
  };
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::size_t t1 = (axis + 1U) % 3U;
    const std::size_t t2 = (axis + 2U) % 3U;
    const double inverse_dn = 1.0 / spacingAlong(geometry, axis);
    const double inverse_dt1 = 1.0 / spacingAlong(geometry, t1);
    const double inverse_dt2 = 1.0 / spacingAlong(geometry, t2);
    const std::span<const double> emf_t1 = emfLane(emf, t1);
    const std::span<const double> emf_t2 = emfLane(emf, t2);
    const MhdDoubleArray& b_t1 = primitives.field(static_cast<MhdAxis>(t1));
    const MhdDoubleArray& b_t2 = primitives.field(static_cast<MhdAxis>(t2));
    for (std::ptrdiff_t c2 = 0; c2 < n[t2]; ++c2) {
      for (std::ptrdiff_t c1 = 0; c1 < n[t1]; ++c1) {
        for (std::ptrdiff_t f = 0; f <= n[axis]; ++f) {
          // E_t1 on the two t1-edges bounding the face (t2-nodes c2, c2+1).
          MhdIndex3 edge{};
          edge[axis] = f;
          edge[t1] = c1;
          edge[t2] = c2;
          double e_t1 = emf_t1[edgeIndex(geometry, t1, edge)];
          edge[t2] = c2 + 1;
          e_t1 = 0.5 * (e_t1 + emf_t1[edgeIndex(geometry, t1, edge)]) * inverse_dt1;
          // E_t2 on the two t2-edges bounding the face (t1-nodes c1, c1+1).
          edge[t2] = c2;
          edge[t1] = c1;
          double e_t2 = emf_t2[edgeIndex(geometry, t2, edge)];
          edge[t1] = c1 + 1;
          e_t2 = 0.5 * (e_t2 + emf_t2[edgeIndex(geometry, t2, edge)]) * inverse_dt2;
          // Face field: mean of the two adjacent cell-centred B_c.
          MhdIndex3 lower_cell{};
          lower_cell[axis] = f - 1;
          lower_cell[t1] = c1;
          lower_cell[t2] = c2;
          MhdIndex3 upper_cell = lower_cell;
          upper_cell[axis] = f;
          const double face_b_t1 =
              0.5 * sqrt_a * (b_t1[primitives.index(lower_cell)] + b_t1[primitives.index(upper_cell)]);
          const double face_b_t2 =
              0.5 * sqrt_a * (b_t2[primitives.index(lower_cell)] + b_t2[primitives.index(upper_cell)]);
          // (E x B)_n for the right-handed basis (n, t1, t2), comoving 1/a.
          const double flux = (e_t1 * face_b_t2 - e_t2 * face_b_t1) / scale_factor;
          if (f == 0 || f == n[axis]) {
            workspace.boundary_energy.at(axis, f == 0 ? 0U : 1U, static_cast<std::size_t>(c1),
                                         static_cast<std::size_t>(c2)) = flux;
          } else {
            energy_target[cellIndex(lower_cell)] -= weight * flux * inverse_dn;
            energy_target[cellIndex(upper_cell)] += weight * flux * inverse_dn;
          }
        }
      }
    }
  }
}

std::array<double, 3> nonIdealElectricField(const std::array<double, 3>& current_density_code,
                                            const std::array<double, 3>& magnetic_field_comoving,
                                            double eta_ohmic_code,
                                            double eta_hall_code,
                                            double eta_ambipolar_code,
                                            double b_regularization) noexcept {
  const std::array<double, 3>& j = current_density_code;
  const std::array<double, 3>& b = magnetic_field_comoving;
  const double magnitude = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
  const double regularized = std::max(magnitude, b_regularization);
  std::array<double, 3> b_hat{0.0, 0.0, 0.0};
  if (magnitude > 0.0) {
    b_hat = {b[0] / regularized, b[1] / regularized, b[2] / regularized};
  }
  const double j_dot_b = j[0] * b_hat[0] + j[1] * b_hat[1] + j[2] * b_hat[2];
  const std::array<double, 3> j_cross_b{j[1] * b_hat[2] - j[2] * b_hat[1],
                                        j[2] * b_hat[0] - j[0] * b_hat[2],
                                        j[0] * b_hat[1] - j[1] * b_hat[0]};
  std::array<double, 3> electric{};
  for (std::size_t component = 0; component < 3U; ++component) {
    electric[component] = eta_ohmic_code * j[component] +
        eta_hall_code * j_cross_b[component] +
        eta_ambipolar_code * (j[component] - j_dot_b * b_hat[component]);
  }
  return electric;
}

void applyBoundaryEnergyFluxes(const MhdPatchGeometry& geometry,
                               const MhdBoundaryFaceEnergyFluxes& fluxes,
                               double weight,
                               std::span<double> energy_target) {
  const std::array<std::size_t, 3> n{geometry.nx, geometry.ny, geometry.nz};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::size_t t1 = (axis + 1U) % 3U;
    const std::size_t t2 = (axis + 2U) % 3U;
    const double inverse_dn = 1.0 / spacingAlong(geometry, axis);
    for (std::size_t side = 0; side < 2U; ++side) {
      const double sign = side == 0U ? 1.0 : -1.0;
      for (std::size_t c2 = 0; c2 < n[t2]; ++c2) {
        for (std::size_t c1 = 0; c1 < n[t1]; ++c1) {
          std::array<std::size_t, 3> c{};
          c[axis] = side == 0U ? 0U : n[axis] - 1U;
          c[t1] = c1;
          c[t2] = c2;
          energy_target[geometry.cellIndex(c[0], c[1], c[2])] +=
              sign * weight * fluxes.at(axis, side, c1, c2) * inverse_dn;
        }
      }
    }
  }
}

}  // namespace cosmosim::mhd
