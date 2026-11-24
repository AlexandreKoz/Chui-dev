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

// Versioned HDF5 table (contract chui_mhd_transport_table_v1; full layout in
// src/mhd/mhd_transport.cpp and docs/configuration.md). Axes:
//   (log_rho_phys_cgs, log_temperature_k, log_b_phys_gauss,
//    composition_id, grain_model_id)
// Datasets EtaOhmicCgs, EtaHallCgs, EtaAmbipolarCgs and the exact zero-field
// group BZeroLimit/{EtaOhmicCgs,EtaHallCgs,EtaAmbipolarCgs}. Continuous axes
// must be finite, sorted and unique; categorical axes are exact-match and never
// interpolated. Ohmic/ambipolar use monotone log interpolation, Hall signed
// linear interpolation. Outside coverage returns valid=false (fatal for the
// caller); no extrapolation or clamping. The SHA-256 of the file bytes is the
// provenance/restart pin. evaluate() returns valid=false until
// validateAdmission() has succeeded.
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
// RKL2 (second-order Runge-Kutta-Legendre) super-time-stepping,
// Meyer, Balsara & Aslam (2014, J. Comput. Phys. 257, 594), eqs. (16)-(17):
//   b_0 = b_1 = b_2 = 1/3,  b_j = (j^2 + j - 2) / (2 j (j + 1)) for j >= 2,
//   a_j = 1 - b_j,  w_1 = 4 / (s^2 + s - 2),
//   Y_1 = Y_0 + mu~_1 tau L(Y_0),  mu~_1 = b_1 w_1,
//   Y_j = mu_j Y_{j-1} + nu_j Y_{j-2} + (1 - mu_j - nu_j) Y_0
//         + mu~_j tau L(Y_{j-1}) + gamma~_j tau L(Y_0),         j = 2..s,
//   mu_j = (2j-1)/j * b_j/b_{j-1},  nu_j = -(j-1)/j * b_j/b_{j-2},
//   mu~_j = mu_j w_1,  gamma~_j = -a_{j-1} mu~_j.
// Linear stability: tau <= dt_explicit * (s^2 + s - 2) / 4, where
// dt_explicit is the forward-Euler (equivalently SSPRK2) diffusive bound.
// Vectors are indexed 0..s; index 0 is unused except b_0.
// ---------------------------------------------------------------------------
struct MhdRkl2Coefficients {
  std::int32_t stages = 0;
  double w1 = 0.0;
  std::vector<double> b;
  std::vector<double> mu;
  std::vector<double> nu;
  std::vector<double> mu_tilde;
  std::vector<double> gamma_tilde;
};

// Throws MHD_INVALID_FIELD for stages < 2.
[[nodiscard]] MhdRkl2Coefficients makeRkl2Coefficients(std::int32_t stages);
// (s^2 + s - 2) / 4 for s >= 2; 1 for s <= 1.
[[nodiscard]] double rkl2StabilityFactor(std::int32_t stages) noexcept;
// Least s >= 2 with rkl2StabilityFactor(s) >= ratio. Returns 0 when
// ratio <= 1 (no super-stepping needed) and max_stages + 1 when even
// max_stages does not cover ratio (caller partitions into cycles).
[[nodiscard]] std::int32_t selectRkl2StageCount(double ratio, std::int32_t max_stages) noexcept;

// Deterministic partition of a diffusion interval of `ratio` explicit steps:
// the least number of equal RKL2 cycles whose per-cycle stage count does not
// exceed `max_stages_per_cycle`.
struct MhdRkl2Plan {
  std::int32_t cycles = 0;
  std::int32_t stages_per_cycle = 0;
  [[nodiscard]] std::int64_t totalStages() const noexcept {
    return static_cast<std::int64_t>(cycles) * static_cast<std::int64_t>(stages_per_cycle);
  }
};
[[nodiscard]] MhdRkl2Plan planRkl2(double ratio, std::int32_t max_stages_per_cycle);

// ---------------------------------------------------------------------------
// Non-ideal operators (CHUI-MHD-REF v1.0 section 8).
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
  std::uint64_t explicit_diffusion_substeps = 0;
  std::uint64_t rkl_cycle_count = 0;
  std::uint64_t rkl_stage_count = 0;
  std::uint64_t hall_subcycle_count = 0;
  std::uint64_t hall_rejected_subcycle_count = 0;
  std::uint64_t coefficient_evaluations = 0;
  // Discrete edge-integrated sum of J.E_ni / a (comoving volume), times the
  // stage weight: a diagnostic of the physical Ohmic+ambipolar heating rate.
  // Never added to the energy equation (the conservative Poynting flux carries
  // it). The Hall part is reported separately and must vanish analytically.
  double joule_heating_integral = 0.0;
  double hall_work_integral = 0.0;
  double max_current_density = 0.0;
  double min_diffusion_dt_code = 0.0;
  double min_hall_dt_code = 0.0;

  MhdTransportStageDiagnostics& operator+=(const MhdTransportStageDiagnostics& rhs) noexcept;
};

// The single owner of MHD unit transforms from comoving code state to the
// physical cgs properties requested by coefficient providers:
//   rho_phys = rho_c / a^3,  p_phys = p_c / a^3,  B_phys = B_c / a^2,
//   T = mu m_p p_phys / (k_B rho_phys)  (calorically perfect gas),
//   B_gauss = sqrt(4 pi) * B_LH,cgs with B_LH,cgs = B_code * sqrt(P_unit_cgs)
// (normalized Lorentz-Heaviside code field; section 4.1). No other MHD module
// introduces 4 pi or mu_0 factors.
struct MhdPhysicalPropertyModel {
  const core::UnitSystem* units = nullptr;
  // Mean molecular weight of the gas in units of the proton mass. There is no
  // default: the workflow must bind it from the authoritative composition
  // closure; a non-positive value is rejected.
  double mean_molecular_weight = 0.0;

  void validate() const;
  [[nodiscard]] MhdCellPhysicalProperties evaluate(double rho_comoving_code,
                                                   double pressure_comoving_code,
                                                   double b_comoving_magnitude_code,
                                                   double scale_factor) const;
  [[nodiscard]] double pressureUnitCgs() const;
  [[nodiscard]] double fieldCodeToGauss(double b_physical_code) const;
};

enum class MhdNonIdealOperatorKind : std::uint8_t {
  kDiffusion = 0,  // Ohmic + ambipolar (the D operator)
  kHall = 1,       // Hall with whistler-HLL upwinding (the H operator)
};

// Non-ideal Poynting energy fluxes on the six patch-boundary planes (one
// value per boundary face; tangential coordinates are the cyclic successors
// of the normal axis). Reconciled from the deterministic owner like the gas
// boundary fluxes so both sides of a shared face see the same flux.
struct MhdBoundaryFaceEnergyFluxes {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  std::array<std::array<MhdDoubleArray, 2>, 3> values{};

  void resize(const MhdPatchGeometry& geometry);
  [[nodiscard]] double& at(std::size_t axis, std::size_t side, std::size_t t1,
                           std::size_t t2) noexcept;
  [[nodiscard]] double at(std::size_t axis, std::size_t side, std::size_t t1,
                          std::size_t t2) const noexcept;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

// Patch-local non-ideal scratch (kScratchArena class). Stage counts never
// multiply residency: the D/H integrators reuse these lanes for every stage.
// To respect the PERF 512 B/active-cell phase+scratch cap, the non-ideal
// operator does NOT own an EMF array or a Y_0 copy: it writes its EMF into the
// level stepper's MhdStageWorkspace::emf and uses the stepper's stage-0
// snapshot as Y_0 (the ideal and non-ideal phases are never co-live). Edge
// currents are recomputed from the staggered field on demand (4 face reads)
// instead of being stored.
struct MhdNonIdealWorkspace {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  // Cell coefficients in comoving code units on cells -1..n (one ghost layer).
  MhdDoubleArray eta_ohmic_code;
  MhdDoubleArray eta_hall_code;
  MhdDoubleArray eta_ambipolar_code;
  // Face Poynting energy fluxes on patch-boundary planes.
  MhdBoundaryFaceEnergyFluxes boundary_energy{};
  // Super-stepping lanes: Y_{j-2} faces+energy, L(Y_0) as EMF + energy rate.
  std::array<MhdDoubleArray, 3> yjm2_faces{};
  MhdDoubleArray yjm2_energy;
  MhdEdgeEmfStorage emf_y0{};
  MhdDoubleArray energy_rate_y0;
  // Energy rate of the current stage (also used by explicit SSPRK2 / Hall).
  MhdDoubleArray energy_rate_stage;

  void resize(const MhdPatchGeometry& geometry, bool super_stepping_lanes);
  [[nodiscard]] std::size_t coefficientIndex(std::ptrdiff_t i, std::ptrdiff_t j,
                                             std::ptrdiff_t k) const noexcept {
    return static_cast<std::size_t>(i + 1) +
        (nx + 2U) * (static_cast<std::size_t>(j + 1) + (ny + 2U) * static_cast<std::size_t>(k + 1));
  }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

// Non-ideal rate kernels on one patch. Coefficients are re-evaluated from the
// SUPPLIED stage state (the padded primitive cache with ghosts filled) at
// every call; nothing is frozen across stages.
class MhdNonIdealRateOperator {
 public:
  MhdNonIdealRateOperator(const MhdCoefficientProvider* provider,
                          MhdPhysicalPropertyModel properties,
                          MhdStiffnessLimits limits,
                          bool enable_ohmic,
                          bool enable_hall,
                          bool enable_ambipolar);

  [[nodiscard]] bool diffusionEnabled() const noexcept {
    return m_enable_ohmic || m_enable_ambipolar;
  }
  [[nodiscard]] bool hallEnabled() const noexcept { return m_enable_hall; }
  [[nodiscard]] const MhdStiffnessLimits& limits() const noexcept { return m_limits; }

  // Evaluate physical properties and provider coefficients on cells -1..n of
  // the stage state and convert them exactly once to comoving code units
  // (eta_code = eta_cgs * T_unit / L_unit^2 / a^2). Disabled terms are stored
  // as exact zeros. Throws MHD_BAD_COEFFICIENT for an invalid/out-of-coverage
  // coefficient, a closure-invalid cell, eta_O < 0 or eta_A < 0.
  void evaluateCoefficients(const MhdPatchGeometry& geometry,
                            const MhdPrimitivePatch& primitives,
                            double scale_factor,
                            MhdNonIdealWorkspace& workspace,
                            MhdTransportStageDiagnostics& diagnostics) const;

  // Explicit stability bounds from the evaluated coefficients (comoving code
  // units, comoving spacings), CHUI-MHD-REF v1.0 sections 8.3-8.4:
  //   dt_D = cfl_diffusion / [max(eta_O + eta_A) * sum_i dx_i^-2]
  //   dt_H = cfl_hall / [pi * max|eta_H| * sum_i dx_i^-2]
  // The Hall bound is the stricter of the diffusive-form bound
  // 0.10/[max|eta_H| sum dx^-2] and the whistler-fan bound with grid Nyquist
  // whistler speed c_w,i = pi |eta_H| / dx_i (sum_i c_w,i/dx_i), which is
  // always smaller by pi. The ambipolar projection has unit operator norm, so
  // eta_O + eta_A bounds the anisotropic operator spectrum. +infinity when
  // the corresponding coefficients are exactly zero.
  [[nodiscard]] double diffusionStableTimeStep(const MhdPatchGeometry& geometry,
                                               const MhdNonIdealWorkspace& workspace) const;
  [[nodiscard]] double hallStableTimeStep(const MhdPatchGeometry& geometry,
                                          const MhdNonIdealWorkspace& workspace) const;

  // Compatible edge current J = curl_x B_c from the padded staggered field
  // (B_c = sqrt(a) B_eff), then the edge EMF of the requested operator:
  //   D: E = eta_O J + eta_A (J - (J.b)b)
  //   H: E = eta_H (J x b) + whistler-HLL upwinding
  //        + (c_w,a/2)(B_b(+a) - B_b(-a)) - (c_w,b/2)(B_a(+b) - B_a(-b))
  // with edge coefficients averaged from the four adjacent cells, b from the
  // edge-averaged field regularized by B_reg = 1e-30 (code units), and
  // reconstructed (PLM-MC) staggered fields at the edge. Output EMFs are edge
  // line integrals per unit time (written to `emf`) with conducting-wall and
  // periodic edge rules applied.
  void computeEmf(MhdNonIdealOperatorKind kind,
                  const MhdPatchGeometry& geometry,
                  const MhdPatchBoundary& boundary,
                  const MhdStaggeredFieldPatch& staggered,
                  double scale_factor,
                  const MhdNonIdealWorkspace& workspace,
                  MhdEdgeEmfStorage& emf,
                  MhdTransportStageDiagnostics& diagnostics) const;

  // Conservative non-ideal Poynting energy flux F_n = (E_c x B_c)_n / a on
  // every face, from the (reconciled) edge EMF in `emf` and the face-averaged
  // cell field B_c. Interior faces apply weight * (-div_x F) directly into
  // `energy_target`; boundary-plane fluxes are recorded in
  // workspace.boundary_energy for owner reconciliation.
  void applyEnergyFlux(const MhdPatchGeometry& geometry,
                       const MhdPrimitivePatch& primitives,
                       const MhdEdgeEmfStorage& emf,
                       double scale_factor,
                       double weight,
                       MhdNonIdealWorkspace& workspace,
                       std::span<double> energy_target) const;

 private:
  const MhdCoefficientProvider* m_provider = nullptr;
  MhdPhysicalPropertyModel m_properties{};
  MhdStiffnessLimits m_limits{};
  bool m_enable_ohmic = false;
  bool m_enable_hall = false;
  bool m_enable_ambipolar = false;
};

// Apply reconciled boundary-plane Poynting energy fluxes:
//   lower plane -> cell 0 gains +w F/dx_n; upper plane -> cell n-1 loses.
void applyBoundaryEnergyFluxes(const MhdPatchGeometry& geometry,
                               const MhdBoundaryFaceEnergyFluxes& fluxes,
                               double weight,
                               std::span<double> energy_target);

// Non-ideal edge electric field from the constitutive law:
//   E_O = eta_O J, E_H = eta_H (J x b_hat), E_A = eta_A (J - (J.b_hat) b_hat),
// b_hat = B / max(|B|, b_regularization). At B = 0 b_hat = 0: the Hall term
// vanishes and the ambipolar term reduces to eta_A J, so a provider must
// supply its physical zero-field limit of eta_A (no direction is invented).
[[nodiscard]] std::array<double, 3> nonIdealElectricField(
    const std::array<double, 3>& current_density_code,
    const std::array<double, 3>& magnetic_field_comoving,
    double eta_ohmic_code,
    double eta_hall_code,
    double eta_ambipolar_code,
    double b_regularization) noexcept;

}  // namespace cosmosim::mhd
