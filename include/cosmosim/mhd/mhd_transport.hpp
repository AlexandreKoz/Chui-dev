#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/mhd/mhd_induction.hpp"
#include "cosmosim/mhd/mhd_riemann.hpp"

namespace cosmosim::core {
struct UnitSystem;
}

namespace cosmosim::mhd {

// Physical (cgs) magnetic diffusivities. Ohmic and ambipolar are nonnegative;
// Hall is signed. Length^2/time in cm^2/s.
struct MhdTransportCoefficients {
  double eta_ohmic_phys_cgs = 0.0;
  double eta_hall_phys_cgs = 0.0;
  double eta_ambipolar_phys_cgs = 0.0;
  bool valid = false;
  bool closure_valid = true;
  std::string provenance_id;
  double max_stable_substep_hint_code = 0.0;
};

// Convert a physical cgs diffusivity to comoving code units exactly once:
// eta_comov_code = eta_cgs * (T_unit_s / L_unit_cm^2) / a^2.
[[nodiscard]] double physicalDiffusivityToComovingCode(double eta_phys_cgs,
                                                       double scale_factor,
                                                       const core::UnitSystem& units);

struct MhdCellPhysicalProperties {
  double rho_phys_cgs = 0.0;
  double temperature_k = 0.0;
  double b_phys_gauss = 0.0;
  bool closure_valid = true;
};

enum class MhdCoefficientProviderKind : std::uint8_t {
  kNone = 0,
  kConstant = 1,
  kTable = 2,
  kChemistry = 3,
};

[[nodiscard]] std::string_view mhdCoefficientProviderName(MhdCoefficientProviderKind kind) noexcept;

// Read-only provider. Implementations may cache a rank-scope table but must
// never mutate gas, field, chemistry or scheduler authority.
class MhdCoefficientProvider {
 public:
  virtual ~MhdCoefficientProvider() = default;

  [[nodiscard]] virtual MhdCoefficientProviderKind kind() const noexcept = 0;
  [[nodiscard]] virtual std::string provenanceId() const = 0;
  [[nodiscard]] virtual std::optional<std::string> tableSha256() const { return std::nullopt; }

  // Validated admission. Throws MHD_BAD_COEFFICIENT / MHD_UNSUPPORTED_PHYSICS.
  virtual void validateAdmission() const = 0;

  // Batch-free per-cell query used only outside the hot edge loop; callers may
  // hoist it to a patch-local coefficient cache at stage granularity.
  [[nodiscard]] virtual MhdTransportCoefficients evaluate(
      std::size_t cell_index,
      const MhdCellPhysicalProperties& properties) const = 0;
};

class ConstantMhdCoefficientProvider final : public MhdCoefficientProvider {
 public:
  ConstantMhdCoefficientProvider(double eta_ohmic_cgs,
                                 double eta_hall_cgs,
                                 double eta_ambipolar_cgs);

  [[nodiscard]] MhdCoefficientProviderKind kind() const noexcept override {
    return MhdCoefficientProviderKind::kConstant;
  }
  [[nodiscard]] std::string provenanceId() const override;
  void validateAdmission() const override;
  [[nodiscard]] MhdTransportCoefficients evaluate(
      std::size_t cell_index,
      const MhdCellPhysicalProperties& properties) const override;

 private:
  double m_eta_ohmic_cgs;
  double m_eta_hall_cgs;
  double m_eta_ambipolar_cgs;
};

// Monotone versioned HDF5 table. Axes:
//   (log_rho_phys_cgs, log_temperature_k, log_b_phys_gauss,
//    composition_id, grain_model_id)
// Datasets: EtaOhmicCgs, EtaHallCgs, EtaAmbipolarCgs, BZeroLimit.
// Continuous axes must be sorted and unique; categorical axes are exact-match
// and never interpolated. Outside coverage is fatal; no extrapolation.
class TableMhdCoefficientProvider final : public MhdCoefficientProvider {
 public:
  struct Options {
    std::string composition_id;
    std::string grain_model_id;
    std::string expected_sha256;  // optional pinned digest
  };

  TableMhdCoefficientProvider(std::filesystem::path table_path, Options options);
  ~TableMhdCoefficientProvider() override;

  TableMhdCoefficientProvider(const TableMhdCoefficientProvider&) = delete;
  TableMhdCoefficientProvider& operator=(const TableMhdCoefficientProvider&) = delete;

  [[nodiscard]] MhdCoefficientProviderKind kind() const noexcept override {
    return MhdCoefficientProviderKind::kTable;
  }
  [[nodiscard]] std::string provenanceId() const override;
  [[nodiscard]] std::optional<std::string> tableSha256() const override { return m_sha256; }
  void validateAdmission() const override;
  [[nodiscard]] MhdTransportCoefficients evaluate(
      std::size_t cell_index,
      const MhdCellPhysicalProperties& properties) const override;

  struct Impl;

 private:
  std::unique_ptr<Impl> m_impl;
  std::filesystem::path m_path;
  Options m_options;
  std::string m_sha256;
};

// Detailed charged-carrier/grain transport. REF-v1 requires a live chemistry
// system that provides every charged species, collision rate, charge-neutrality
// relation and grain property. That contract is not established in this
// snapshot, so requested chemistry transport fails closed at admission rather
// than inventing an ionization fraction.
class ChemistryMhdCoefficientProvider final : public MhdCoefficientProvider {
 public:
  ChemistryMhdCoefficientProvider() = default;

  [[nodiscard]] MhdCoefficientProviderKind kind() const noexcept override {
    return MhdCoefficientProviderKind::kChemistry;
  }
  [[nodiscard]] std::string provenanceId() const override;
  void validateAdmission() const override;
  [[nodiscard]] MhdTransportCoefficients evaluate(
      std::size_t cell_index,
      const MhdCellPhysicalProperties& properties) const override;
};

// ---------------------------------------------------------------------------
// RKL2 (second-order Runge-Kutta-Legendre) super-time-stepping.
//
// Coefficients follow the Meyer-Balsara-Aslam recurrence:
//   b_j = (j^2 + j - 2) / (2 j (j + 1)),  a_j = 1 - b_j,
//   w_1 = 4 / (s^2 + s - 2),
//   mu_j   = (2j - 1)/j * b_j / b_{j-1},
//   nu_j   = -(j - 1)/j * b_j / b_{j-2},
//   mu_tilde_j = mu_j * w_1,   gamma_tilde_j = -a_{j-1} * mu_tilde_j.
// The formal explicit stability interval of the s-stage scheme is
//   dt <= 2 s (s + 1) / (s^2 + s - 2) * dt_explicit_single_stage,
// and the implementation verifies the recurrence itself in the unit tests.
// ---------------------------------------------------------------------------
struct MhdRkl2Coefficients {
  std::int32_t stages = 1;
  double w1 = 1.0;
  std::vector<double> mu;
  std::vector<double> nu;
  std::vector<double> mu_tilde;
  std::vector<double> gamma_tilde;
};

[[nodiscard]] MhdRkl2Coefficients makeRkl2Coefficients(std::int32_t stages);
[[nodiscard]] double rkl2StabilityFactor(std::int32_t stages) noexcept;
// Least deterministic stage count whose stability interval covers `ratio`,
// bounded by `max_stages`. Returns 0 when ratio <= 1 (single explicit step).
[[nodiscard]] std::int32_t selectRkl2StageCount(double ratio, std::int32_t max_stages) noexcept;

// ---------------------------------------------------------------------------
// Non-ideal operators.
// ---------------------------------------------------------------------------

struct MhdStiffnessLimits {
  double cfl_diffusion = 0.15;
  double cfl_hall = 0.10;
  std::int32_t rkl_switch_ratio = 4;
  std::int32_t rkl_max_stages = 128;
  std::int32_t hall_subcycle_cost_cap = 4096;
  std::int32_t rkl_stage_cost_cap = 4096;
};

struct MhdTransportStageDiagnostics {
  std::uint64_t rkl_cycle_count = 0;
  std::uint64_t rkl_stage_count = 0;
  std::uint64_t hall_subcycle_count = 0;
  std::uint64_t hall_rejected_subcycle_count = 0;
  std::uint64_t coefficient_evaluations = 0;
  double jacobi_heating_integral = 0.0;  // diagnostic only
  double max_current_density = 0.0;
};

// Stage-local workspace for the non-ideal operators. All arrays are patch
// sized; RKL/Hall subcycles reuse them and never multiply residency.
struct MhdNonIdealWorkspace {
  MhdEdgeEmfStorage emf_interval{};
  MhdEdgeEmfStorage emf_stage{};
  MhdDoubleArray current_x{};
  MhdDoubleArray current_y{};
  MhdDoubleArray current_z{};
  MhdDoubleArray energy_stage{};
  MhdDoubleArray energy_rhs{};
  MhdDoubleArray phi_x_stage{};
  MhdDoubleArray phi_y_stage{};
  MhdDoubleArray phi_z_stage{};

  void resize(const MhdPatchGeometry& geometry);
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

// Non-ideal induction rate provider. Computes d(Phi)/dt and the conservative
// non-ideal Poynting energy rate at a stage state. Coefficients are reevaluated
// at every call from the supplied physical cell properties.
class MhdNonIdealRateOperator {
 public:
  MhdNonIdealRateOperator(const MhdCoefficientProvider* provider,
                          const core::UnitSystem* units,
                          MhdStiffnessLimits limits,
                          bool enable_ohmic,
                          bool enable_hall,
                          bool enable_ambipolar);

  // Ohmic + ambipolar (the D operator). Fills `emf_integral` with the
  // *rate* edge EMF (per unit time) so the caller multiplies by its stage dt.
  void computeDiffusionRate(const MhdPatchGeometry& geometry,
                            std::span<const double> phi_x,
                            std::span<const double> phi_y,
                            std::span<const double> phi_z,
                            std::span<const MhdCellPhysicalProperties> properties,
                            double scale_factor,
                            MhdEdgeEmfStorage& emf_rate,
                            std::span<double> energy_rate_density) const;

  // Hall (the H operator), whistler-expanded HLL induction EMF.
  void computeHallRate(const MhdPatchGeometry& geometry,
                       std::span<const double> phi_x,
                       std::span<const double> phi_y,
                       std::span<const double> phi_z,
                       std::span<const MhdCellPhysicalProperties> properties,
                       double scale_factor,
                       MhdEdgeEmfStorage& emf_rate,
                       std::span<double> energy_rate_density) const;

  [[nodiscard]] double diffusionStableTimeStep(
      const MhdPatchGeometry& geometry,
      std::span<const MhdCellPhysicalProperties> properties,
      double scale_factor,
      const core::UnitSystem& units) const;

  [[nodiscard]] double hallStableTimeStep(
      const MhdPatchGeometry& geometry,
      std::span<const MhdCellPhysicalProperties> properties,
      double scale_factor,
      const core::UnitSystem& units) const;

  [[nodiscard]] MhdStiffnessLimits limits() const noexcept { return m_limits; }
  [[nodiscard]] bool nonidealEnabled() const noexcept {
    return m_enable_ohmic || m_enable_hall || m_enable_ambipolar;
  }

 private:
  void evaluateCoefficients(std::span<const MhdCellPhysicalProperties> properties,
                            double scale_factor,
                            std::vector<MhdTransportCoefficients>& out_code_units) const;
  void computeEdgeCurrent(const MhdPatchGeometry& geometry,
                          std::span<const double> phi_x,
                          std::span<const double> phi_y,
                          std::span<const double> phi_z,
                          MhdDoubleArray& jx,
                          MhdDoubleArray& jy,
                          MhdDoubleArray& jz) const;

  const MhdCoefficientProvider* m_provider = nullptr;
  const core::UnitSystem* m_units = nullptr;
  MhdStiffnessLimits m_limits{};
  bool m_enable_ohmic = false;
  bool m_enable_hall = false;
  bool m_enable_ambipolar = false;
};

// SSPRK2 / RKL2 driver for the D operator. Updates face fluxes and the
// conservative total-energy lane together; no artificial conductivity and no
// hidden energy rescale.
class MhdDiffusionIntegrator {
 public:
  explicit MhdDiffusionIntegrator(MhdStiffnessLimits limits);

  void integrate(const MhdNonIdealRateOperator& rate_operator,
                 const MhdPatchGeometry& geometry,
                 MhdNonIdealWorkspace& workspace,
                 std::span<double> phi_x,
                 std::span<double> phi_y,
                 std::span<double> phi_z,
                 std::span<double> total_energy_density,
                 std::span<const MhdCellPhysicalProperties> properties,
                 double scale_factor,
                 double dt_code,
                 MhdTransportStageDiagnostics& diagnostics) const;

  [[nodiscard]] MhdStiffnessLimits limits() const noexcept { return m_limits; }

 private:
  MhdStiffnessLimits m_limits;
};

// Explicit SSPRK2 Hall driver with bounded subcycling.
class MhdHallIntegrator {
 public:
  explicit MhdHallIntegrator(MhdStiffnessLimits limits);

  void integrate(const MhdNonIdealRateOperator& rate_operator,
                 const MhdPatchGeometry& geometry,
                 MhdNonIdealWorkspace& workspace,
                 std::span<double> phi_x,
                 std::span<double> phi_y,
                 std::span<double> phi_z,
                 std::span<double> total_energy_density,
                 std::span<const MhdCellPhysicalProperties> properties,
                 double scale_factor,
                 double dt_code,
                 MhdTransportStageDiagnostics& diagnostics) const;

 private:
  MhdStiffnessLimits m_limits;
};

// Non-ideal edge EMF from the constitutive law:
//   E_O = eta_O J, E_H = eta_H (J x b_hat), E_A = eta_A (J - (J.b_hat)b_hat).
// `b_hat` uses the regularized field direction; at exactly B=0 the Ohmic term
// is well defined (J), Hall/ambipolar vanish in their physical limit and an
// undefined coefficient table query is rejected before this function.
[[nodiscard]] std::array<double, 3> nonIdealElectricField(
    const std::array<double, 3>& current_density_code,
    const std::array<double, 3>& magnetic_field_comoving,
    const MhdTransportCoefficients& coefficients_code_units,
    double scale_factor) noexcept;

}  // namespace cosmosim::mhd
