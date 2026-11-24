#include "cosmosim/mhd/mhd_validation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

namespace cosmosim::mhd {

std::string_view mhdErrorCodeName(MhdErrorCode code) noexcept {
  switch (code) {
    case MhdErrorCode::kInvalidField: return "MHD_INVALID_FIELD";
    case MhdErrorCode::kDivergenceContract: return "MHD_DIVERGENCE_CONTRACT";
    case MhdErrorCode::kStaleGeometry: return "MHD_STALE_GEOMETRY";
    case MhdErrorCode::kBadCoefficient: return "MHD_BAD_COEFFICIENT";
    case MhdErrorCode::kUnsupportedPhysics: return "MHD_UNSUPPORTED_PHYSICS";
    case MhdErrorCode::kStiffnessBudgetExceeded: return "MHD_STIFFNESS_BUDGET_EXCEEDED";
    case MhdErrorCode::kNonPositiveState: return "MHD_NONPOSITIVE_STATE";
    case MhdErrorCode::kAmrRegisterIncomplete: return "MHD_AMR_REGISTER_INCOMPLETE";
    case MhdErrorCode::kRestartSchemaMismatch: return "MHD_RESTART_SCHEMA_MISMATCH";
    case MhdErrorCode::kMemoryBudgetExceeded: return "MHD_MEMORY_BUDGET_EXCEEDED";
    case MhdErrorCode::kMpiEpochMismatch: return "MHD_MPI_EPOCH_MISMATCH";
  }
  return "MHD_UNKNOWN";
}

MhdError::MhdError(MhdErrorCode code, MhdFailureContext context)
    : std::runtime_error(buildMessage(code, context)),
      m_code(code),
      m_context(std::move(context)) {}

std::string MhdError::buildMessage(MhdErrorCode code, const MhdFailureContext& context) {
  std::ostringstream stream;
  stream << mhdErrorCodeName(code);
  if (!context.owner.empty()) stream << " [" << context.owner << ']';
  if (!context.detail.empty()) stream << ": " << context.detail;
  stream << " (rank=" << context.rank << ", patch=" << context.patch_id
         << ", geometry_generation=" << context.geometry_generation;
  if (context.face_key != 0U) stream << ", face_key=" << context.face_key;
  if (context.edge_key != 0U) stream << ", edge_key=" << context.edge_key;
  stream << ", t=" << context.time_code << ", a=" << context.scale_factor;
  if (context.requested_bytes != 0U) {
    stream << ", requested_bytes=" << context.requested_bytes
           << ", admitted_bytes=" << context.admitted_bytes
           << ", safe_limit_bytes=" << context.safe_limit_bytes;
  }
  if (context.predicted_stages != 0U) {
    stream << ", predicted_stages=" << context.predicted_stages;
  }
  stream << ')';
  return stream.str();
}

std::string MhdProvenance::describe() const {
  std::ostringstream stream;
  stream << "MhdSchemeVersion=" << scheme_version
         << ";profile=" << reference_profile
         << ";representation=" << magnetic_representation
         << ";reconstruction=" << reconstruction
         << ";riemann=" << riemann_solver << '/' << riemann_fallback
         << ";emf=" << emf_scheme
         << ";time_integrator=" << time_integrator
         << ";divergence=" << divergence_policy
         << ";nonideal_provider=" << nonideal_provider
         << ";diffusivity_units=" << dimensional_diffusivity_convention;
  if (!coefficient_table_sha256.empty()) {
    stream << ";coefficient_table_sha256=" << coefficient_table_sha256;
  }
  return stream.str();
}

double normalizedDivergence(std::span<const double> signed_cell_flux_imbalance,
                            double normalizer) {
  double max_abs = 0.0;
  for (const double value : signed_cell_flux_imbalance) {
    max_abs = std::max(max_abs, std::abs(value));
  }
  if (!(normalizer > 0.0) || !std::isfinite(normalizer)) {
    return max_abs > 0.0 ? std::numeric_limits<double>::infinity() : 0.0;
  }
  return max_abs / normalizer;
}

void requireDivergenceContract(double max_normalized_divergence,
                               double tolerance,
                               const MhdFailureContext& context) {
  if (!std::isfinite(max_normalized_divergence) || max_normalized_divergence > tolerance) {
    MhdFailureContext detail = context;
    std::ostringstream stream;
    stream << "max normalized discrete divergence " << max_normalized_divergence
           << " exceeds tolerance " << tolerance;
    detail.detail = stream.str();
    throw MhdError(MhdErrorCode::kDivergenceContract, std::move(detail));
  }
}

void requirePositiveFinite(double value,
                           std::string_view what,
                           const MhdFailureContext& context) {
  if (!std::isfinite(value) || !(value > 0.0)) {
    MhdFailureContext detail = context;
    detail.detail = "non-positive or non-finite " + std::string(what) + " (" +
        std::to_string(value) + ")";
    throw MhdError(MhdErrorCode::kNonPositiveState, std::move(detail));
  }
}

}  // namespace cosmosim::mhd
