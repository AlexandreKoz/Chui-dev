#include "cosmosim/mhd/mhd_transport.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "cosmosim/core/units.hpp"
#include "internal/mhd_sha256.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::mhd {
namespace {

constexpr double k_tiny_positive = 1.0e-300;

// L(Phi) = -curl(E) for one edge-EMF field. Applying the constrained-transport
// operator to a zero face array yields exactly the negative curl, so the same
// operator defines both the accept step and the super-time-stepping increment.
void inductionRateFromEmf(const MhdPatchGeometry& geometry,
                          const MhdEdgeEmfStorage& emf_rate,
                          std::span<double> inc_x,
                          std::span<double> inc_y,
                          std::span<double> inc_z) {
  std::fill(inc_x.begin(), inc_x.end(), 0.0);
  std::fill(inc_y.begin(), inc_y.end(), 0.0);
  std::fill(inc_z.begin(), inc_z.end(), 0.0);
  applyConstrainedTransportUpdate(geometry, inc_x, inc_y, inc_z,
                                  emf_rate.emfX(), emf_rate.emfY(), emf_rate.emfZ(), 1.0);
}

// Cell-centred curl of the derived cell field. Indices clamp to the patch; the
// caller owns halo correctness, this is a patch-local operator.
void cellCentredCurrent(const MhdPatchGeometry& geometry,
                        const MhdDoubleArray& bx,
                        const MhdDoubleArray& by,
                        const MhdDoubleArray& bz,
                        MhdDoubleArray& jx,
                        MhdDoubleArray& jy,
                        MhdDoubleArray& jz) {
  const std::size_t nx = geometry.nx;
  const std::size_t ny = geometry.ny;
  const std::size_t nz = geometry.nz;
  const std::size_t cells = geometry.cellCount();
  jx.resize(cells, 0.0);
  jy.resize(cells, 0.0);
  jz.resize(cells, 0.0);
  if (bx.size() != cells || by.size() != cells || bz.size() != cells) {
    return;
  }
  const auto at = [&](std::size_t i, std::size_t j, std::size_t k) {
    return geometry.cellIndex(std::min(i, nx - 1U), std::min(j, ny - 1U), std::min(k, nz - 1U));
  };
  const auto minus = [](std::size_t index) { return index > 0U ? index - 1U : 0U; };
  for (std::size_t k = 0; k < nz; ++k) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t i = 0; i < nx; ++i) {
        const std::size_t ip = at(i + 1U, j, k);
        const std::size_t im = at(minus(i), j, k);
        const std::size_t jp = at(i, j + 1U, k);
        const std::size_t jm = at(i, minus(j), k);
        const std::size_t kp = at(i, j, k + 1U);
        const std::size_t km = at(i, j, minus(k));
        const double di = static_cast<double>((i + 1U < nx ? i + 1U : i) -
                                              (i > 0U ? i - 1U : i));
        const double dj = static_cast<double>((j + 1U < ny ? j + 1U : j) -
                                              (j > 0U ? j - 1U : j));
        const double dk = static_cast<double>((k + 1U < nz ? k + 1U : k) -
                                              (k > 0U ? k - 1U : k));
        const double dx = geometry.dx_comoving * std::max(1.0, di);
        const double dy = geometry.dy_comoving * std::max(1.0, dj);
        const double dz = geometry.dz_comoving * std::max(1.0, dk);
        const std::size_t cell = geometry.cellIndex(i, j, k);
        jx[cell] = (bz[jp] - bz[jm]) / dy - (by[kp] - by[km]) / dz;
        jy[cell] = (bx[kp] - bx[km]) / dz - (bz[ip] - bz[im]) / dx;
        jz[cell] = (by[ip] - by[im]) / dx - (bx[jp] - bx[jm]) / dy;
      }
    }
  }
}

struct CellFieldLanes {
  MhdDoubleArray bx;
  MhdDoubleArray by;
  MhdDoubleArray bz;
};

bool refreshCellField(const MhdPatchGeometry& geometry,
                      std::span<const double> phi_x,
                      std::span<const double> phi_y,
                      std::span<const double> phi_z,
                      CellFieldLanes& lanes) {
  if (phi_x.size() != geometry.faceCountX() ||
      phi_y.size() != geometry.faceCountY() ||
      phi_z.size() != geometry.faceCountZ()) {
    return false;
  }
  const std::size_t cells = geometry.cellCount();
  lanes.bx.resize(cells, 0.0);
  lanes.by.resize(cells, 0.0);
  lanes.bz.resize(cells, 0.0);
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::size_t cell = geometry.cellIndex(i, j, k);
        lanes.bx[cell] = 0.5 * (phi_x[geometry.faceIndexX(i, j, k)] +
                               phi_x[geometry.faceIndexX(i + 1U, j, k)]) /
            geometry.faceAreaX();
        lanes.by[cell] = 0.5 * (phi_y[geometry.faceIndexY(i, j, k)] +
                               phi_y[geometry.faceIndexY(i, j + 1U, k)]) /
            geometry.faceAreaY();
        lanes.bz[cell] = 0.5 * (phi_z[geometry.faceIndexZ(i, j, k)] +
                               phi_z[geometry.faceIndexZ(i, j, k + 1U)]) /
            geometry.faceAreaZ();
      }
    }
  }
  return true;
}

// Arithmetic cell-to-edge average with patch clamping.
void averageCellToEdges(const MhdPatchGeometry& geometry,
                        MhdAxis axis,
                        const MhdDoubleArray& cell_values,
                        MhdEdgeEmfStorage& emf_rate_field) {
  const std::size_t nx = geometry.nx;
  const std::size_t ny = geometry.ny;
  const std::size_t nz = geometry.nz;
  const auto cell_index_clamped = [&](std::size_t i, std::size_t j, std::size_t k) {
    const std::size_t ci = i < nx ? i : nx - 1U;
    const std::size_t cj = j < ny ? j : ny - 1U;
    const std::size_t ck = k < nz ? k : nz - 1U;
    return geometry.cellIndex(ci, cj, ck);
  };
  if (axis == MhdAxis::kX) {
    for (std::size_t k = 0; k <= nz; ++k) {
      for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
          const std::size_t i0 = i > 0 ? i - 1U : 0U;
          const std::size_t j0 = j > 0 ? j - 1U : 0U;
          const std::size_t k0 = k > 0 ? k - 1U : 0U;
          double sum = 0.0;
          sum += cell_values[cell_index_clamped(i0, j0, k0)];
          sum += cell_values[cell_index_clamped(i0, j, k0)];
          sum += cell_values[cell_index_clamped(i0, j0, k)];
          sum += cell_values[cell_index_clamped(i0, j, k)];
          emf_rate_field.emfX()[geometry.edgeIndexX(i, j, k)] = 0.25 * sum;
        }
      }
    }
  } else if (axis == MhdAxis::kY) {
    for (std::size_t k = 0; k <= nz; ++k) {
      for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
          const std::size_t i0 = i > 0 ? i - 1U : 0U;
          const std::size_t j0 = j > 0 ? j - 1U : 0U;
          const std::size_t k0 = k > 0 ? k - 1U : 0U;
          double sum = 0.0;
          sum += cell_values[cell_index_clamped(i0, j0, k0)];
          sum += cell_values[cell_index_clamped(i0, j0, k)];
          sum += cell_values[cell_index_clamped(i, j0, k0)];
          sum += cell_values[cell_index_clamped(i, j0, k)];
          emf_rate_field.emfY()[geometry.edgeIndexY(i, j, k)] = 0.25 * sum;
        }
      }
    }
  } else {
    for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
          const std::size_t i0 = i > 0 ? i - 1U : 0U;
          const std::size_t j0 = j > 0 ? j - 1U : 0U;
          const std::size_t k0 = k > 0 ? k - 1U : 0U;
          double sum = 0.0;
          sum += cell_values[cell_index_clamped(i0, j0, k0)];
          sum += cell_values[cell_index_clamped(i, j0, k0)];
          sum += cell_values[cell_index_clamped(i0, j, k0)];
          sum += cell_values[cell_index_clamped(i, j, k0)];
          emf_rate_field.emfZ()[geometry.edgeIndexZ(i, j, k)] = 0.25 * sum;
        }
      }
    }
  }
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
  if (!std::isfinite(eta_phys_cgs) || !(scale_factor > 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double length_cm_per_code = units.length_si_per_code * 100.0;
  if (!(length_cm_per_code > 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double time_s_per_code = units.timeSiPerCode();
  const double conversion = time_s_per_code / (length_cm_per_code * length_cm_per_code);
  return eta_phys_cgs * conversion / (scale_factor * scale_factor);
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
  stream << "constant_eta_cgs[O=" << m_eta_ohmic_cgs << ",H=" << m_eta_hall_cgs
         << ",A=" << m_eta_ambipolar_cgs << ']';
  return stream.str();
}

void ConstantMhdCoefficientProvider::validateAdmission() const {
  const auto check = [](double value, const char* name) {
    if (!std::isfinite(value)) {
      MhdFailureContext context{};
      context.owner = "mhd_coefficient_provider";
      context.detail = std::string("non-finite constant coefficient ") + name;
      throw MhdError(MhdErrorCode::kBadCoefficient, std::move(context));
    }
  };
  check(m_eta_ohmic_cgs, "eta_ohmic");
  check(m_eta_hall_cgs, "eta_hall");
  check(m_eta_ambipolar_cgs, "eta_ambipolar");
  if (m_eta_ohmic_cgs < 0.0 || m_eta_ambipolar_cgs < 0.0) {
    MhdFailureContext context{};
    context.owner = "mhd_coefficient_provider";
    context.detail = "Ohmic and ambipolar diffusivities must be nonnegative";
    throw MhdError(MhdErrorCode::kBadCoefficient, std::move(context));
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
  coefficients.valid = true;
  coefficients.closure_valid = properties.closure_valid;
  coefficients.provenance_id = provenanceId();
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
      "collision-rate and grain model; the current snapshot does not establish that contract";
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
// ---------------------------------------------------------------------------

struct TableMhdCoefficientProvider::Impl {
  std::vector<double> log_rho;
  std::vector<double> log_temperature;
  std::vector<double> log_b;
  std::vector<std::string> composition_ids;
  std::vector<std::string> grain_model_ids;
  // Flattened [nrho][nT][nB][ncomp][ngrain]
  std::vector<double> eta_ohmic;
  std::vector<double> eta_hall;
  std::vector<double> eta_ambipolar;
  std::vector<double> b_zero_limit;
};

#if COSMOSIM_ENABLE_HDF5
namespace {

std::vector<double> readDoubleDataset(hid_t file, const char* name) {
  const hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  if (dataset < 0) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_table_provider",
                                     .detail = std::string("missing dataset ") + name});
  }
  const hid_t space = H5Dget_space(dataset);
  const int rank = H5Sget_simple_extent_ndims(space);
  std::vector<hsize_t> dims(static_cast<std::size_t>(std::max(rank, 0)));
  if (rank > 0) {
    H5Sget_simple_extent_dims(space, dims.data(), nullptr);
  }
  std::size_t total = 1U;
  for (const hsize_t dim : dims) {
    total *= static_cast<std::size_t>(dim);
  }
  std::vector<double> values(total, 0.0);
  if (total > 0U) {
    const herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                  H5P_DEFAULT, values.data());
    if (status < 0) {
      H5Sclose(space);
      H5Dclose(dataset);
      throw MhdError(MhdErrorCode::kBadCoefficient,
                     MhdFailureContext{.owner = "mhd_table_provider",
                                       .detail = std::string("cannot read dataset ") + name});
    }
  }
  H5Sclose(space);
  H5Dclose(dataset);
  return values;
}

std::vector<std::string> readStringDataset(hid_t file, const char* name) {
  const hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  if (dataset < 0) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_table_provider",
                                     .detail = std::string("missing dataset ") + name});
  }
  const hid_t type = H5Dget_type(dataset);
  const std::size_t size = H5Tget_size(type);
  const hid_t space = H5Dget_space(dataset);
  const hsize_t count = H5Sget_simple_extent_npoints(space);
  std::vector<char> buffer(static_cast<std::size_t>(count) * size, '\0');
  if (count > 0U) {
    H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
  }
  std::vector<std::string> values;
  values.reserve(static_cast<std::size_t>(count));
  for (hsize_t index = 0; index < count; ++index) {
    const char* begin = buffer.data() + static_cast<std::size_t>(index) * size;
    std::string value(begin, strnlen(begin, size));
    values.push_back(std::move(value));
  }
  H5Sclose(space);
  H5Tclose(type);
  H5Dclose(dataset);
  return values;
}

std::string sha256OfFile(const std::filesystem::path& path) {
  // The shared internal SHA-256 helper implements the textual digest contract
  // used by every CHUI provenance pin.
  return internal::sha256FileHex(path);
}

std::size_t bracketIndex(const std::vector<double>& axis, double value) {
  if (axis.size() < 2U) {
    return 0U;
  }
  if (value <= axis.front()) {
    return 0U;
  }
  if (value >= axis.back()) {
    return axis.size() - 2U;
  }
  const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
  const std::size_t index = static_cast<std::size_t>(upper - axis.begin());
  return index == 0U ? 0U : std::min(index - 1U, axis.size() - 2U);
}

double interpolationWeight(const std::vector<double>& axis, std::size_t lower, double value) {
  if (axis.size() < 2U) {
    return 0.0;
  }
  const double x0 = axis[lower];
  const double x1 = axis[lower + 1U];
  if (!(x1 > x0)) {
    return 0.0;
  }
  return std::clamp((value - x0) / (x1 - x0), 0.0, 1.0);
}

}  // namespace
#endif  // COSMOSIM_ENABLE_HDF5

TableMhdCoefficientProvider::TableMhdCoefficientProvider(std::filesystem::path table_path,
                                                         Options options)
    : m_impl(std::make_unique<Impl>()), m_path(std::move(table_path)), m_options(std::move(options)) {
#if COSMOSIM_ENABLE_HDF5
  if (!std::filesystem::exists(m_path)) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_table_provider",
                                     .detail = "coefficient table does not exist: " +
                                         m_path.string()});
  }
  const hid_t file = H5Fopen(m_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_table_provider",
                                     .detail = "cannot open coefficient table " + m_path.string()});
  }
  const auto close_on_error = [file]() { H5Fclose(file); };
  try {
    m_impl->log_rho = readDoubleDataset(file, "/Axes/LogDensityCgs");
    m_impl->log_temperature = readDoubleDataset(file, "/Axes/LogTemperatureK");
    m_impl->log_b = readDoubleDataset(file, "/Axes/LogMagneticFieldGauss");
    m_impl->composition_ids = readStringDataset(file, "/Axes/CompositionIds");
    m_impl->grain_model_ids = readStringDataset(file, "/Axes/GrainModelIds");
    m_impl->eta_ohmic = readDoubleDataset(file, "/EtaOhmicCgs");
    m_impl->eta_hall = readDoubleDataset(file, "/EtaHallCgs");
    m_impl->eta_ambipolar = readDoubleDataset(file, "/EtaAmbipolarCgs");
    m_impl->b_zero_limit = readDoubleDataset(file, "/BZeroLimit");
  } catch (...) {
    close_on_error();
    throw;
  }
  H5Fclose(file);
  m_sha256 = sha256OfFile(m_path);
#else
  throw MhdError(MhdErrorCode::kBadCoefficient,
                 MhdFailureContext{.owner = "mhd_table_provider",
                                   .detail =
                                       "table coefficient provider requires "
                                       "COSMOSIM_ENABLE_HDF5=ON"});
#endif
}

TableMhdCoefficientProvider::~TableMhdCoefficientProvider() = default;

std::string TableMhdCoefficientProvider::provenanceId() const {
  std::ostringstream stream;
  stream << "table[" << m_path.string() << ";composition=" << m_options.composition_id
         << ";grain=" << m_options.grain_model_id << ']';
  return stream.str();
}

void TableMhdCoefficientProvider::validateAdmission() const {
#if COSMOSIM_ENABLE_HDF5
  const std::size_t expected =
      m_impl->log_rho.size() * m_impl->log_temperature.size() * m_impl->log_b.size() *
      m_impl->composition_ids.size() * m_impl->grain_model_ids.size();
  if (m_impl->eta_ohmic.size() != expected || m_impl->eta_hall.size() != expected ||
      m_impl->eta_ambipolar.size() != expected) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_table_provider",
                                     .detail = "coefficient table shape does not match its axes"});
  }
  const auto strictly_increasing = [](const std::vector<double>& axis) {
    for (std::size_t index = 1U; index < axis.size(); ++index) {
      if (!(axis[index] > axis[index - 1U])) {
        return false;
      }
    }
    return true;
  };
  if (!strictly_increasing(m_impl->log_rho) ||
      !strictly_increasing(m_impl->log_temperature) ||
      !strictly_increasing(m_impl->log_b)) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_table_provider",
                                     .detail = "coefficient table axes must be sorted and unique"});
  }
  if (std::find(m_impl->composition_ids.begin(), m_impl->composition_ids.end(),
                m_options.composition_id) == m_impl->composition_ids.end()) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_table_provider",
                                     .detail = "composition id not present in table: " +
                                         m_options.composition_id});
  }
  if (std::find(m_impl->grain_model_ids.begin(), m_impl->grain_model_ids.end(),
                m_options.grain_model_id) == m_impl->grain_model_ids.end()) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_table_provider",
                                     .detail = "grain model id not present in table: " +
                                         m_options.grain_model_id});
  }
  if (!m_options.expected_sha256.empty() && m_options.expected_sha256 != m_sha256) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_table_provider",
                                     .detail = "coefficient table SHA-256 does not match pin"});
  }
#else
  throw MhdError(MhdErrorCode::kBadCoefficient,
                 MhdFailureContext{.owner = "mhd_table_provider",
                                   .detail = "table provider requires COSMOSIM_ENABLE_HDF5=ON"});
#endif
}

MhdTransportCoefficients TableMhdCoefficientProvider::evaluate(
    std::size_t cell_index,
    const MhdCellPhysicalProperties& properties) const {
  MhdTransportCoefficients coefficients{};
  coefficients.closure_valid = properties.closure_valid;
  coefficients.provenance_id = provenanceId();
#if COSMOSIM_ENABLE_HDF5
  if (!(properties.rho_phys_cgs > 0.0) || !(properties.temperature_k > 0.0) ||
      !(properties.b_phys_gauss > 0.0)) {
    coefficients.valid = false;
    coefficients.closure_valid = false;
    return coefficients;
  }
  (void)cell_index;
  const double log_rho = std::log10(properties.rho_phys_cgs);
  const double log_temperature = std::log10(properties.temperature_k);
  const double log_b = std::log10(properties.b_phys_gauss);
  const auto& rho_axis = m_impl->log_rho;
  const auto& temperature_axis = m_impl->log_temperature;
  const auto& b_axis = m_impl->log_b;
  if (log_rho < rho_axis.front() || log_rho > rho_axis.back() ||
      log_temperature < temperature_axis.front() || log_temperature > temperature_axis.back() ||
      log_b < b_axis.front() || log_b > b_axis.back()) {
    // No extrapolation: outside coverage is fatal for the caller.
    coefficients.valid = false;
    coefficients.closure_valid = false;
    return coefficients;
  }
  const std::size_t i = bracketIndex(rho_axis, log_rho);
  const std::size_t j = bracketIndex(temperature_axis, log_temperature);
  const std::size_t k = bracketIndex(b_axis, log_b);
  const double ti = interpolationWeight(rho_axis, i, log_rho);
  const double tj = interpolationWeight(temperature_axis, j, log_temperature);
  const double tk = interpolationWeight(b_axis, k, log_b);
  const std::size_t composition_index = static_cast<std::size_t>(
      std::find(m_impl->composition_ids.begin(), m_impl->composition_ids.end(),
                m_options.composition_id) - m_impl->composition_ids.begin());
  const std::size_t grain_index = static_cast<std::size_t>(
      std::find(m_impl->grain_model_ids.begin(), m_impl->grain_model_ids.end(),
                m_options.grain_model_id) - m_impl->grain_model_ids.begin());
  const std::size_t n_b = b_axis.size();
  const std::size_t n_t = temperature_axis.size();
  const std::size_t n_rho = rho_axis.size();
  const std::size_t n_comp = m_impl->composition_ids.size();
  const std::size_t n_grain = m_impl->grain_model_ids.size();
  const auto flatten = [&](std::size_t a, std::size_t b, std::size_t c) {
    return ((a * n_t + b) * n_b + c) * n_comp * n_grain +
        composition_index * n_grain + grain_index;
  };
  const auto sample = [&](const std::vector<double>& data) {
    const double c000 = data[flatten(i, j, k)];
    const double c100 = data[flatten(i + 1U, j, k)];
    const double c010 = data[flatten(i, j + 1U, k)];
    const double c110 = data[flatten(i + 1U, j + 1U, k)];
    const double c001 = data[flatten(i, j, k + 1U)];
    const double c101 = data[flatten(i + 1U, j, k + 1U)];
    const double c011 = data[flatten(i, j + 1U, k + 1U)];
    const double c111 = data[flatten(i + 1U, j + 1U, k + 1U)];
    const double c00 = c000 * (1.0 - ti) + c100 * ti;
    const double c10 = c010 * (1.0 - ti) + c110 * ti;
    const double c01 = c001 * (1.0 - ti) + c101 * ti;
    const double c11 = c011 * (1.0 - ti) + c111 * ti;
    const double c0 = c00 * (1.0 - tj) + c10 * tj;
    const double c1 = c01 * (1.0 - tj) + c11 * tj;
    return c0 * (1.0 - tk) + c1 * tk;
  };
  (void)n_rho;
  coefficients.eta_ohmic_phys_cgs = sample(m_impl->eta_ohmic);
  coefficients.eta_hall_phys_cgs = sample(m_impl->eta_hall);
  coefficients.eta_ambipolar_phys_cgs = sample(m_impl->eta_ambipolar);
  coefficients.valid = std::isfinite(coefficients.eta_ohmic_phys_cgs) &&
      std::isfinite(coefficients.eta_hall_phys_cgs) &&
      std::isfinite(coefficients.eta_ambipolar_phys_cgs);
#endif
  return coefficients;
}

// ---------------------------------------------------------------------------
// RKL2 coefficients
// ---------------------------------------------------------------------------

MhdRkl2Coefficients makeRkl2Coefficients(std::int32_t stages) {
  MhdRkl2Coefficients coefficients{};
  if (stages <= 1) {
    coefficients.stages = 1;
    coefficients.w1 = 1.0;
    coefficients.mu = {1.0};
    coefficients.nu = {0.0};
    coefficients.mu_tilde = {1.0};
    coefficients.gamma_tilde = {1.0};
    return coefficients;
  }
  const std::int32_t s = stages;
  const double s_double = static_cast<double>(s);
  // Meyer-Balsara-Aslam (2012) RKL2 recurrence, with the first-stage
  // initialization b_1 = 1 and gamma_tilde_j = b_{j-1} * mu_tilde_j. For s = 2
  // this reduces exactly to SSPRK2 (verified in the unit tests).
  coefficients.stages = s;
  coefficients.w1 = 4.0 / (s_double * s_double + s_double - 2.0);
  coefficients.mu.assign(static_cast<std::size_t>(s) + 1U, 0.0);
  coefficients.nu.assign(static_cast<std::size_t>(s) + 1U, 0.0);
  coefficients.mu_tilde.assign(static_cast<std::size_t>(s) + 1U, 0.0);
  coefficients.gamma_tilde.assign(static_cast<std::size_t>(s) + 1U, 0.0);
  coefficients.mu[1] = 1.0;
  coefficients.mu_tilde[1] = 1.0;
  coefficients.gamma_tilde[1] = 1.0;
  std::vector<double> b(static_cast<std::size_t>(s) + 1U, 0.0);
  b[1] = 1.0;
  for (std::int32_t j = 2; j <= s; ++j) {
    const double jd = static_cast<double>(j);
    b[static_cast<std::size_t>(j)] = (jd * jd + jd - 2.0) / (2.0 * jd * (jd + 1.0));
  }
  for (std::int32_t j = 2; j <= s; ++j) {
    const double jd = static_cast<double>(j);
    const double bj = b[static_cast<std::size_t>(j)];
    coefficients.mu[static_cast<std::size_t>(j)] =
        (2.0 * jd - 1.0) / jd * bj / b[static_cast<std::size_t>(j - 1)];
    if (j >= 3) {
      coefficients.nu[static_cast<std::size_t>(j)] =
          -(jd - 1.0) / jd * bj / b[static_cast<std::size_t>(j - 2)];
    }
    coefficients.mu_tilde[static_cast<std::size_t>(j)] =
        coefficients.mu[static_cast<std::size_t>(j)] * coefficients.w1;
    coefficients.gamma_tilde[static_cast<std::size_t>(j)] =
        b[static_cast<std::size_t>(j - 1)] *
        coefficients.mu_tilde[static_cast<std::size_t>(j)];
  }
  return coefficients;
}

double rkl2StabilityFactor(std::int32_t stages) noexcept {
  if (stages <= 1) {
    return 1.0;
  }
  const double s = static_cast<double>(stages);
  return 2.0 * s * (s + 1.0) / (s * s + s - 2.0);
}

std::int32_t selectRkl2StageCount(double ratio, std::int32_t max_stages) noexcept {
  if (!(ratio > 1.0) || max_stages <= 1) {
    return 0;
  }
  for (std::int32_t s = 2; s <= max_stages; ++s) {
    if (rkl2StabilityFactor(s) >= ratio) {
      return s;
    }
  }
  return max_stages + 1;
}

// ---------------------------------------------------------------------------
// Non-ideal workspace and operators
// ---------------------------------------------------------------------------

void MhdNonIdealWorkspace::resize(const MhdPatchGeometry& geometry) {
  emf_interval.resize(geometry);
  emf_stage.resize(geometry);
  const std::size_t cells = geometry.cellCount();
  current_x.resize(cells, 0.0);
  current_y.resize(cells, 0.0);
  current_z.resize(cells, 0.0);
  energy_stage.resize(cells, 0.0);
  energy_rhs.resize(cells, 0.0);
  phi_x_stage.resize(geometry.faceCountX(), 0.0);
  phi_y_stage.resize(geometry.faceCountY(), 0.0);
  phi_z_stage.resize(geometry.faceCountZ(), 0.0);
}

std::uint64_t MhdNonIdealWorkspace::ownedCapacityBytes() const noexcept {
  return emf_interval.ownedCapacityBytes() + emf_stage.ownedCapacityBytes() +
      current_x.ownedCapacityBytes() + current_y.ownedCapacityBytes() +
      current_z.ownedCapacityBytes() + energy_stage.ownedCapacityBytes() +
      energy_rhs.ownedCapacityBytes() + phi_x_stage.ownedCapacityBytes() +
      phi_y_stage.ownedCapacityBytes() + phi_z_stage.ownedCapacityBytes();
}

MhdNonIdealRateOperator::MhdNonIdealRateOperator(const MhdCoefficientProvider* provider,
                                                 const core::UnitSystem* units,
                                                 MhdStiffnessLimits limits,
                                                 bool enable_ohmic,
                                                 bool enable_hall,
                                                 bool enable_ambipolar)
    : m_provider(provider),
      m_units(units),
      m_limits(limits),
      m_enable_ohmic(enable_ohmic),
      m_enable_hall(enable_hall),
      m_enable_ambipolar(enable_ambipolar) {}

void MhdNonIdealRateOperator::evaluateCoefficients(
    std::span<const MhdCellPhysicalProperties> properties,
    double scale_factor,
    std::vector<MhdTransportCoefficients>& out_code_units) const {
  if (m_provider == nullptr || m_units == nullptr) {
    throw MhdError(MhdErrorCode::kBadCoefficient,
                   MhdFailureContext{.owner = "mhd_nonideal",
                                     .detail = "non-ideal rate operator has no coefficient "
                                               "provider or unit system"});
  }
  out_code_units.resize(properties.size());
  for (std::size_t cell = 0; cell < properties.size(); ++cell) {
    MhdTransportCoefficients physical = m_provider->evaluate(cell, properties[cell]);
    if (!physical.valid || !physical.closure_valid) {
      MhdFailureContext context{};
      context.owner = "mhd_nonideal";
      context.detail = "coefficient provider rejected cell " + std::to_string(cell);
      throw MhdError(MhdErrorCode::kBadCoefficient, std::move(context));
    }
    MhdTransportCoefficients code_units = physical;
    code_units.eta_ohmic_phys_cgs =
        physicalDiffusivityToComovingCode(physical.eta_ohmic_phys_cgs, scale_factor, *m_units);
    code_units.eta_hall_phys_cgs =
        physicalDiffusivityToComovingCode(physical.eta_hall_phys_cgs, scale_factor, *m_units);
    code_units.eta_ambipolar_phys_cgs =
        physicalDiffusivityToComovingCode(physical.eta_ambipolar_phys_cgs, scale_factor, *m_units);
    if (!std::isfinite(code_units.eta_ohmic_phys_cgs) ||
        !std::isfinite(code_units.eta_hall_phys_cgs) ||
        !std::isfinite(code_units.eta_ambipolar_phys_cgs)) {
      MhdFailureContext context{};
      context.owner = "mhd_nonideal";
      context.detail = "diffusivity unit conversion produced a non-finite value";
      throw MhdError(MhdErrorCode::kBadCoefficient, std::move(context));
    }
    out_code_units[cell] = code_units;
  }
}

void MhdNonIdealRateOperator::computeEdgeCurrent(const MhdPatchGeometry& geometry,
                                                 std::span<const double> phi_x,
                                                 std::span<const double> phi_y,
                                                 std::span<const double> phi_z,
                                                 MhdDoubleArray& jx,
                                                 MhdDoubleArray& jy,
                                                 MhdDoubleArray& jz) const {
  CellFieldLanes lanes{};
  if (!refreshCellField(geometry, phi_x, phi_y, phi_z, lanes)) {
    return;
  }
  cellCentredCurrent(geometry, lanes.bx, lanes.by, lanes.bz, jx, jy, jz);
}

void MhdNonIdealRateOperator::computeDiffusionRate(
    const MhdPatchGeometry& geometry,
    std::span<const double> phi_x,
    std::span<const double> phi_y,
    std::span<const double> phi_z,
    std::span<const MhdCellPhysicalProperties> properties,
    double scale_factor,
    MhdEdgeEmfStorage& emf_rate,
    std::span<double> energy_rate_density) const {
  emf_rate.resize(geometry);
  const std::size_t cells = geometry.cellCount();
  if (energy_rate_density.size() != cells) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_nonideal",
                                     .detail = "energy rate lane has the wrong size"});
  }
  CellFieldLanes lanes{};
  if (!refreshCellField(geometry, phi_x, phi_y, phi_z, lanes)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_nonideal",
                                     .detail = "face flux lanes do not match patch geometry"});
  }
  MhdDoubleArray jx;
  MhdDoubleArray jy;
  MhdDoubleArray jz;
  cellCentredCurrent(geometry, lanes.bx, lanes.by, lanes.bz, jx, jy, jz);
  std::vector<MhdTransportCoefficients> coefficients;
  evaluateCoefficients(properties, scale_factor, coefficients);

  MhdDoubleArray ex(cells, 0.0);
  MhdDoubleArray ey(cells, 0.0);
  MhdDoubleArray ez(cells, 0.0);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    MhdTransportCoefficients effective{};
    effective.eta_ohmic_phys_cgs = m_enable_ohmic ? coefficients[cell].eta_ohmic_phys_cgs : 0.0;
    effective.eta_hall_phys_cgs = 0.0;  // diffusion operator excludes Hall
    effective.eta_ambipolar_phys_cgs =
        m_enable_ambipolar ? coefficients[cell].eta_ambipolar_phys_cgs : 0.0;
    effective.valid = true;
    const std::array<double, 3> current{jx[cell], jy[cell], jz[cell]};
    const std::array<double, 3> field{lanes.bx[cell], lanes.by[cell], lanes.bz[cell]};
    const std::array<double, 3> electric =
        nonIdealElectricField(current, field, effective, scale_factor);
    ex[cell] = electric[0];
    ey[cell] = electric[1];
    ez[cell] = electric[2];
    // Conservative Poynting heating: -J . E (Ohmic + ambipolar).
    energy_rate_density[cell] = -(current[0] * electric[0] + current[1] * electric[1] +
                                  current[2] * electric[2]);
  }
  averageCellToEdges(geometry, MhdAxis::kX, ex, emf_rate);
  averageCellToEdges(geometry, MhdAxis::kY, ey, emf_rate);
  averageCellToEdges(geometry, MhdAxis::kZ, ez, emf_rate);
  for (double& value : emf_rate.emfX()) {
    value *= geometry.dx_comoving;
  }
  for (double& value : emf_rate.emfY()) {
    value *= geometry.dy_comoving;
  }
  for (double& value : emf_rate.emfZ()) {
    value *= geometry.dz_comoving;
  }
}

void MhdNonIdealRateOperator::computeHallRate(
    const MhdPatchGeometry& geometry,
    std::span<const double> phi_x,
    std::span<const double> phi_y,
    std::span<const double> phi_z,
    std::span<const MhdCellPhysicalProperties> properties,
    double scale_factor,
    MhdEdgeEmfStorage& emf_rate,
    std::span<double> energy_rate_density) const {
  emf_rate.resize(geometry);
  const std::size_t cells = geometry.cellCount();
  if (energy_rate_density.size() != cells) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_nonideal",
                                     .detail = "energy rate lane has the wrong size"});
  }
  CellFieldLanes lanes{};
  if (!refreshCellField(geometry, phi_x, phi_y, phi_z, lanes)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_nonideal",
                                     .detail = "face flux lanes do not match patch geometry"});
  }
  MhdDoubleArray jx;
  MhdDoubleArray jy;
  MhdDoubleArray jz;
  cellCentredCurrent(geometry, lanes.bx, lanes.by, lanes.bz, jx, jy, jz);
  std::vector<MhdTransportCoefficients> coefficients;
  evaluateCoefficients(properties, scale_factor, coefficients);

  MhdDoubleArray ex(cells, 0.0);
  MhdDoubleArray ey(cells, 0.0);
  MhdDoubleArray ez(cells, 0.0);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    MhdTransportCoefficients effective{};
    effective.eta_hall_phys_cgs = m_enable_hall ? coefficients[cell].eta_hall_phys_cgs : 0.0;
    effective.valid = true;
    const std::array<double, 3> current{jx[cell], jy[cell], jz[cell]};
    const std::array<double, 3> field{lanes.bx[cell], lanes.by[cell], lanes.bz[cell]};
    const std::array<double, 3> electric =
        nonIdealElectricField(current, field, effective, scale_factor);
    ex[cell] = electric[0];
    ey[cell] = electric[1];
    ez[cell] = electric[2];
    // The ideal Hall contribution does no work; J.E_H is identically zero, and
    // the discrete value is retained only as a numerical residual diagnostic.
    energy_rate_density[cell] = -(current[0] * electric[0] + current[1] * electric[1] +
                                  current[2] * electric[2]);
  }
  averageCellToEdges(geometry, MhdAxis::kX, ex, emf_rate);
  averageCellToEdges(geometry, MhdAxis::kY, ey, emf_rate);
  averageCellToEdges(geometry, MhdAxis::kZ, ez, emf_rate);
  for (double& value : emf_rate.emfX()) {
    value *= geometry.dx_comoving;
  }
  for (double& value : emf_rate.emfY()) {
    value *= geometry.dy_comoving;
  }
  for (double& value : emf_rate.emfZ()) {
    value *= geometry.dz_comoving;
  }
}

double MhdNonIdealRateOperator::diffusionStableTimeStep(
    const MhdPatchGeometry& geometry,
    std::span<const MhdCellPhysicalProperties> properties,
    double scale_factor,
    const core::UnitSystem& units) const {
  (void)units;
  std::vector<MhdTransportCoefficients> coefficients;
  evaluateCoefficients(properties, scale_factor, coefficients);
  double maximum = 0.0;
  for (const MhdTransportCoefficients& cell : coefficients) {
    maximum = std::max(maximum, std::abs(cell.eta_ohmic_phys_cgs) +
                                    std::abs(cell.eta_ambipolar_phys_cgs));
  }
  const double spacing = std::min({geometry.dx_comoving, geometry.dy_comoving,
                                   geometry.dz_comoving});
  if (maximum <= k_tiny_positive) {
    return std::numeric_limits<double>::infinity();
  }
  return m_limits.cfl_diffusion * spacing * spacing / maximum;
}

double MhdNonIdealRateOperator::hallStableTimeStep(
    const MhdPatchGeometry& geometry,
    std::span<const MhdCellPhysicalProperties> properties,
    double scale_factor,
    const core::UnitSystem& units) const {
  (void)units;
  std::vector<MhdTransportCoefficients> coefficients;
  evaluateCoefficients(properties, scale_factor, coefficients);
  double maximum = 0.0;
  double maximum_field = 0.0;
  for (std::size_t cell = 0; cell < coefficients.size(); ++cell) {
    maximum = std::max(maximum, std::abs(coefficients[cell].eta_hall_phys_cgs));
    maximum_field = std::max(maximum_field, std::abs(properties[cell].b_phys_gauss));
  }
  const double spacing = std::min({geometry.dx_comoving, geometry.dy_comoving,
                                   geometry.dz_comoving});
  const double whistler = maximum * std::max(maximum_field, k_tiny_positive);
  if (whistler <= k_tiny_positive) {
    return std::numeric_limits<double>::infinity();
  }
  return m_limits.cfl_hall * spacing * spacing / whistler;
}

// ---------------------------------------------------------------------------
// Integrators
// ---------------------------------------------------------------------------

MhdDiffusionIntegrator::MhdDiffusionIntegrator(MhdStiffnessLimits limits)
    : m_limits(limits) {}

void MhdDiffusionIntegrator::integrate(
    const MhdNonIdealRateOperator& rate_operator,
    const MhdPatchGeometry& geometry,
    MhdNonIdealWorkspace& workspace,
    std::span<double> phi_x,
    std::span<double> phi_y,
    std::span<double> phi_z,
    std::span<double> total_energy_density,
    std::span<const MhdCellPhysicalProperties> properties,
    double scale_factor,
    double dt_code,
    MhdTransportStageDiagnostics& diagnostics) const {
  if (!rate_operator.nonidealEnabled()) {
    return;
  }
  workspace.resize(geometry);
  const std::size_t cells = geometry.cellCount();
  if (phi_x.size() != geometry.faceCountX() || phi_y.size() != geometry.faceCountY() ||
      phi_z.size() != geometry.faceCountZ() || total_energy_density.size() != cells ||
      properties.size() != cells) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_diffusion",
                                     .detail = "diffusion integrator operand sizes disagree"});
  }
  if (!(dt_code > 0.0)) {
    return;
  }
  const double stable_dt =
      rate_operator.diffusionStableTimeStep(geometry, properties, scale_factor,
                                            core::UnitSystem{});
  const double ratio = std::isfinite(stable_dt) ? dt_code / stable_dt : 0.0;
  const std::int32_t stages = selectRkl2StageCount(ratio, m_limits.rkl_max_stages);
  if (stages > m_limits.rkl_max_stages) {
    MhdFailureContext context{};
    context.owner = "mhd_diffusion";
    context.detail = "diffusion stiffness requires more RKL2 stages than the configured cap";
    context.predicted_stages = static_cast<std::uint64_t>(stages);
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(context));
  }
  const std::int32_t effective_stages = stages <= 0 ? 1 : stages;
  if (static_cast<std::int64_t>(effective_stages) * static_cast<std::int64_t>(cells) >
      static_cast<std::int64_t>(m_limits.rkl_stage_cost_cap)) {
    diagnostics.rkl_rejected_subcycle_count += 1U;
  }

  const auto apply_rate = [&](std::span<const double> state_x,
                              std::span<const double> state_y,
                              std::span<const double> state_z,
                              std::span<const double> state_energy) {
    std::vector<MhdCellPhysicalProperties> stage_properties(properties.begin(), properties.end());
    rate_operator.computeDiffusionRate(geometry, state_x, state_y, state_z,
                                       stage_properties, scale_factor,
                                       workspace.emf_interval, workspace.energy_rhs.span());
    std::copy(state_energy.begin(), state_energy.end(), workspace.energy_stage.begin());
  };

  apply_rate(phi_x, phi_y, phi_z, total_energy_density);

  std::vector<double> increment_x(phi_x.size(), 0.0);
  std::vector<double> increment_y(phi_y.size(), 0.0);
  std::vector<double> increment_z(phi_z.size(), 0.0);

  if (effective_stages == 1) {
    // Single explicit step: Y_1 = Y_0 + dt * (-curl(E(Y_0))).
    inductionRateFromEmf(geometry, workspace.emf_interval, increment_x, increment_y,
                         increment_z);
    for (std::size_t index = 0; index < phi_x.size(); ++index) {
      phi_x[index] += dt_code * increment_x[index];
    }
    for (std::size_t index = 0; index < phi_y.size(); ++index) {
      phi_y[index] += dt_code * increment_y[index];
    }
    for (std::size_t index = 0; index < phi_z.size(); ++index) {
      phi_z[index] += dt_code * increment_z[index];
    }
    for (std::size_t cell = 0; cell < cells; ++cell) {
      total_energy_density[cell] += dt_code * workspace.energy_rhs[cell];
    }
    diagnostics.rkl_cycle_count += 1U;
    diagnostics.rkl_stage_count += 1U;
    return;
  }

  // RKL2 super-time-stepping over the same D operator.
  const MhdRkl2Coefficients coefficients = makeRkl2Coefficients(effective_stages);
  std::vector<double> y0_x(phi_x.begin(), phi_x.end());
  std::vector<double> y0_y(phi_y.begin(), phi_y.end());
  std::vector<double> y0_z(phi_z.begin(), phi_z.end());
  std::vector<double> y0_e(total_energy_density.begin(), total_energy_density.end());
  // Y_1 = Y_0 + dt L(Y_0)
  inductionRateFromEmf(geometry, workspace.emf_interval, increment_x, increment_y,
                       increment_z);
  for (std::size_t index = 0; index < phi_x.size(); ++index) {
    phi_x[index] = y0_x[index] + dt_code * increment_x[index];
  }
  for (std::size_t index = 0; index < phi_y.size(); ++index) {
    phi_y[index] = y0_y[index] + dt_code * increment_y[index];
  }
  for (std::size_t index = 0; index < phi_z.size(); ++index) {
    phi_z[index] = y0_z[index] + dt_code * increment_z[index];
  }
  for (std::size_t cell = 0; cell < cells; ++cell) {
    total_energy_density[cell] = y0_e[cell] + dt_code * workspace.energy_rhs[cell];
  }
  std::vector<double> previous_x = y0_x;
  std::vector<double> previous_y = y0_y;
  std::vector<double> previous_z = y0_z;
  std::vector<double> previous_e = y0_e;
  std::vector<double> current_x = phi_x.begin(), current_y = phi_y.begin(),
                      current_z = phi_z.begin(), current_e = total_energy_density.begin();
  diagnostics.rkl_stage_count += 1U;
  for (std::int32_t stage = 2; stage <= effective_stages; ++stage) {
    apply_rate(phi_x, phi_y, phi_z, total_energy_density);
    inductionRateFromEmf(geometry, workspace.emf_interval, increment_x, increment_y,
                         increment_z);
    const double mu_tilde = coefficients.mu_tilde[static_cast<std::size_t>(stage)];
    const double nu_tilde = stage >= 3
        ? coefficients.nu[static_cast<std::size_t>(stage)] * coefficients.w1
        : 0.0;
    const double gamma_tilde = coefficients.gamma_tilde[static_cast<std::size_t>(stage)];
    const double y0_weight = 1.0 - mu_tilde - nu_tilde;
    for (std::size_t index = 0; index < current_x.size(); ++index) {
      const double value = mu_tilde * current_x[index] + nu_tilde * previous_x[index] +
          gamma_tilde * dt_code * increment_x[index] + y0_weight * y0_x[index];
      previous_x[index] = current_x[index];
      current_x[index] = value;
      phi_x[index] = value;
    }
    for (std::size_t index = 0; index < current_y.size(); ++index) {
      const double value = mu_tilde * current_y[index] + nu_tilde * previous_y[index] +
          gamma_tilde * dt_code * increment_y[index] + y0_weight * y0_y[index];
      previous_y[index] = current_y[index];
      current_y[index] = value;
      phi_y[index] = value;
    }
    for (std::size_t index = 0; index < current_z.size(); ++index) {
      const double value = mu_tilde * current_z[index] + nu_tilde * previous_z[index] +
          gamma_tilde * dt_code * increment_z[index] + y0_weight * y0_z[index];
      previous_z[index] = current_z[index];
      current_z[index] = value;
      phi_z[index] = value;
    }
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const double value = mu_tilde * current_e[cell] + nu_tilde * previous_e[cell] +
          gamma_tilde * dt_code * workspace.energy_rhs[cell] + y0_weight * y0_e[cell];
      previous_e[cell] = current_e[cell];
      current_e[cell] = value;
      total_energy_density[cell] = value;
    }
    diagnostics.rkl_stage_count += 1U;
  }
  diagnostics.rkl_cycle_count += 1U;
}

MhdHallIntegrator::MhdHallIntegrator(MhdStiffnessLimits limits)
    : m_limits(limits) {}

void MhdHallIntegrator::integrate(
    const MhdNonIdealRateOperator& rate_operator,
    const MhdPatchGeometry& geometry,
    MhdNonIdealWorkspace& workspace,
    std::span<double> phi_x,
    std::span<double> phi_y,
    std::span<double> phi_z,
    std::span<double> total_energy_density,
    std::span<const MhdCellPhysicalProperties> properties,
    double scale_factor,
    double dt_code,
    MhdTransportStageDiagnostics& diagnostics) const {
  if (!rate_operator.nonidealEnabled() || !(dt_code > 0.0)) {
    return;
  }
  workspace.resize(geometry);
  const std::size_t cells = geometry.cellCount();
  if (phi_x.size() != geometry.faceCountX() || phi_y.size() != geometry.faceCountY() ||
      phi_z.size() != geometry.faceCountZ() || total_energy_density.size() != cells ||
      properties.size() != cells) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   MhdFailureContext{.owner = "mhd_hall",
                                     .detail = "Hall integrator operand sizes disagree"});
  }
  const double stable_dt = rate_operator.hallStableTimeStep(
      geometry, properties, scale_factor, core::UnitSystem{});
  std::int64_t subcycles = 1;
  if (std::isfinite(stable_dt) && stable_dt > 0.0) {
    subcycles = static_cast<std::int64_t>(std::ceil(dt_code / stable_dt));
  }
  if (subcycles > static_cast<std::int64_t>(m_limits.hall_subcycle_cost_cap)) {
    MhdFailureContext context{};
    context.owner = "mhd_hall";
    context.detail = "Hall subcycle count exceeds the configured cost cap";
    context.predicted_stages = static_cast<std::uint64_t>(subcycles);
    throw MhdError(MhdErrorCode::kStiffnessBudgetExceeded, std::move(context));
  }
  const double subcycle_dt = dt_code / static_cast<double>(subcycles);
  std::vector<double> increment_x(phi_x.size(), 0.0);
  std::vector<double> increment_y(phi_y.size(), 0.0);
  std::vector<double> increment_z(phi_z.size(), 0.0);
  for (std::int64_t subcycle = 0; subcycle < subcycles; ++subcycle) {
    std::vector<MhdCellPhysicalProperties> stage_properties(properties.begin(), properties.end());
    rate_operator.computeHallRate(geometry, phi_x, phi_y, phi_z, stage_properties,
                                  scale_factor, workspace.emf_stage,
                                  workspace.energy_rhs.span());
    inductionRateFromEmf(geometry, workspace.emf_stage, increment_x, increment_y,
                         increment_z);
    // Forward-Euler predictor, then SSPRK2-style averaging with the corrector.
    std::vector<double> predicted_x(phi_x.size(), 0.0);
    std::vector<double> predicted_y(phi_y.size(), 0.0);
    std::vector<double> predicted_z(phi_z.size(), 0.0);
    std::vector<double> predicted_e(total_energy_density.begin(), total_energy_density.end());
    for (std::size_t index = 0; index < phi_x.size(); ++index) {
      predicted_x[index] = phi_x[index] + subcycle_dt * increment_x[index];
    }
    for (std::size_t index = 0; index < phi_y.size(); ++index) {
      predicted_y[index] = phi_y[index] + subcycle_dt * increment_y[index];
    }
    for (std::size_t index = 0; index < phi_z.size(); ++index) {
      predicted_z[index] = phi_z[index] + subcycle_dt * increment_z[index];
    }
    for (std::size_t cell = 0; cell < cells; ++cell) {
      predicted_e[cell] += subcycle_dt * workspace.energy_rhs[cell];
    }
    rate_operator.computeHallRate(geometry, predicted_x, predicted_y, predicted_z,
                                  stage_properties, scale_factor, workspace.emf_interval,
                                  workspace.energy_rhs.span());
    inductionRateFromEmf(geometry, workspace.emf_interval, increment_x, increment_y,
                         increment_z);
    for (std::size_t index = 0; index < phi_x.size(); ++index) {
      phi_x[index] = 0.5 * (phi_x[index] + predicted_x[index] +
                            subcycle_dt * increment_x[index]);
    }
    for (std::size_t index = 0; index < phi_y.size(); ++index) {
      phi_y[index] = 0.5 * (phi_y[index] + predicted_y[index] +
                            subcycle_dt * increment_y[index]);
    }
    for (std::size_t index = 0; index < phi_z.size(); ++index) {
      phi_z[index] = 0.5 * (phi_z[index] + predicted_z[index] +
                            subcycle_dt * increment_z[index]);
    }
    for (std::size_t cell = 0; cell < cells; ++cell) {
      total_energy_density[cell] =
          0.5 * (total_energy_density[cell] + predicted_e[cell] +
                 subcycle_dt * workspace.energy_rhs[cell]);
    }
    diagnostics.hall_subcycle_count += 1U;
  }
}

std::array<double, 3> nonIdealElectricField(
    const std::array<double, 3>& current_density_code,
    const std::array<double, 3>& magnetic_field_comoving,
    const MhdTransportCoefficients& coefficients_code_units,
    double scale_factor) noexcept {
  (void)scale_factor;
  const std::array<double, 3> j = current_density_code;
  const std::array<double, 3> b = magnetic_field_comoving;
  const double b_magnitude = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
  std::array<double, 3> b_hat{0.0, 0.0, 0.0};
  if (b_magnitude > 0.0) {
    b_hat = {b[0] / b_magnitude, b[1] / b_magnitude, b[2] / b_magnitude};
  }
  const double j_dot_b = j[0] * b_hat[0] + j[1] * b_hat[1] + j[2] * b_hat[2];
  std::array<double, 3> field{};
  if (coefficients_code_units.eta_hall_phys_cgs != 0.0) {
    // J x b_hat
    const std::array<double, 3> j_cross_b{
        j[1] * b_hat[2] - j[2] * b_hat[1],
        j[2] * b_hat[0] - j[0] * b_hat[2],
        j[0] * b_hat[1] - j[1] * b_hat[0]};
    field[0] += coefficients_code_units.eta_hall_phys_cgs * j_cross_b[0];
    field[1] += coefficients_code_units.eta_hall_phys_cgs * j_cross_b[1];
    field[2] += coefficients_code_units.eta_hall_phys_cgs * j_cross_b[2];
  }
  const double ohmic = coefficients_code_units.eta_ohmic_phys_cgs;
  const double ambipolar = coefficients_code_units.eta_ambipolar_phys_cgs;
  for (int component = 0; component < 3; ++component) {
    field[static_cast<std::size_t>(component)] += ohmic * j[static_cast<std::size_t>(component)];
    field[static_cast<std::size_t>(component)] +=
        ambipolar * (j[static_cast<std::size_t>(component)] - j_dot_b * b_hat[static_cast<std::size_t>(component)]);
  }
  return field;
}

}  // namespace cosmosim::mhd
