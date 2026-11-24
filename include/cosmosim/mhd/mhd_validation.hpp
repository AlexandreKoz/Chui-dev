#pragma once

#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cosmosim::mhd {

// Fatal MHD diagnostic family (CHUI-MHD-REF v1.0, section 15). Every failure
// path in this subsystem carries one of these codes; no MHD failure is reported
// through a generic exception without an identified code.
enum class MhdErrorCode : std::uint8_t {
  kInvalidField = 0,
  kDivergenceContract,
  kStaleGeometry,
  kBadCoefficient,
  kUnsupportedPhysics,
  kStiffnessBudgetExceeded,
  kNonPositiveState,
  kAmrRegisterIncomplete,
  kRestartSchemaMismatch,
  kMemoryBudgetExceeded,
  kMpiEpochMismatch,
};

[[nodiscard]] std::string_view mhdErrorCodeName(MhdErrorCode code) noexcept;

// Stable reference-profile identifiers. A profile name is a promise about the
// numerical contract, not evidence that it has been validated.
inline constexpr std::string_view k_mhd_reference_ideal_profile = "mhd_ref_ideal_v01";
inline constexpr std::string_view k_mhd_reference_nonideal_profile = "mhd_ref_nonideal_v01";
inline constexpr std::string_view k_mhd_reference_contract_id = "chui-mhd-ref-v1.0";
inline constexpr std::uint32_t k_mhd_scheme_version = 1U;

struct MhdFailureContext {
  std::string owner;                 // subsystem/owner label
  std::string detail;                // human-readable specifics
  std::int64_t rank = -1;
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::uint64_t face_key = 0U;
  std::uint64_t edge_key = 0U;
  double time_code = 0.0;
  double scale_factor = 1.0;
  std::uint64_t requested_bytes = 0U;
  std::uint64_t admitted_bytes = 0U;
  std::uint64_t safe_limit_bytes = 0U;
  std::uint64_t predicted_stages = 0U;
};

// Exception type used by every MHD code path that must fail closed.
class MhdError : public std::runtime_error {
 public:
  MhdError(MhdErrorCode code, MhdFailureContext context);

  [[nodiscard]] MhdErrorCode code() const noexcept { return m_code; }
  [[nodiscard]] const MhdFailureContext& context() const noexcept { return m_context; }

 private:
  [[nodiscard]] static std::string buildMessage(MhdErrorCode code,
                                                const MhdFailureContext& context);
  MhdErrorCode m_code;
  MhdFailureContext m_context;
};

// Provenance fields that must travel with a magnetized run. They are emitted
// into normalized provenance and compared on restart admission.
struct MhdProvenance {
  std::string reference_profile = std::string(k_mhd_reference_ideal_profile);
  std::string magnetic_representation = "face_flux_ct";
  std::string reconstruction = "ppm";
  std::string riemann_solver = "hlld";
  std::string riemann_fallback = "hlle";
  std::string emf_scheme = "uct_hll";
  std::string time_integrator = "ssprk2";
  std::string divergence_policy = "strict_ct";
  std::string nonideal_provider = "none";
  std::string coefficient_table_sha256;
  std::string dimensional_diffusivity_convention = "cm2_per_s";
  std::uint32_t scheme_version = k_mhd_scheme_version;

  [[nodiscard]] std::string describe() const;
};

// Discrete-divergence acceptance helper shared by the runtime and tests.
// Returns the maximum normalized cell flux imbalance. `normalizer` is supplied
// by the caller (maximum oriented flux magnitude, reference field*area, or the
// specification epsilon scale) so a zero field cannot report artificial
// excellence.
[[nodiscard]] double normalizedDivergence(
    std::span<const double> signed_cell_flux_imbalance,
    double normalizer);

// Throws MHD_DIVERGENCE_CONTRACT when max_norm exceeds the declared bound.
void requireDivergenceContract(double max_normalized_divergence,
                               double tolerance,
                               const MhdFailureContext& context);

// Validate a finite, strictly positive amplitude used by reconstruction /
// Riemann code. Throws MHD_NONPOSITIVE_STATE on violation.
void requirePositiveFinite(double value,
                           std::string_view what,
                           const MhdFailureContext& context);

}  // namespace cosmosim::mhd
