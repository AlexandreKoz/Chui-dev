#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/mhd/mhd_amr.hpp"
#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_induction.hpp"
#include "cosmosim/mhd/mhd_riemann.hpp"
#include "cosmosim/mhd/mhd_state.hpp"
#include "cosmosim/mhd/mhd_transport.hpp"
#include "cosmosim/mhd/mhd_validation.hpp"

namespace cosmosim::mhd {

enum class MhdReconstructionMode : std::uint8_t { kPlmMc = 0, kPpm = 1 };
enum class MhdTimeIntegrator : std::uint8_t { kSsprk2 = 0 };
enum class MhdDivergencePolicy : std::uint8_t { kStrictCt = 0 };
enum class MhdEmfScheme : std::uint8_t { kUctHll = 0 };
enum class MhdSpatialMethod : std::uint8_t { kFiniteVolumeUnsplit = 0 };
enum class MhdMagneticRepresentation : std::uint8_t { kFaceFluxCt = 0 };

[[nodiscard]] std::string_view mhdReconstructionName(MhdReconstructionMode mode) noexcept;

struct MhdNumericsPolicy {
  MhdReconstructionMode reconstruction = MhdReconstructionMode::kPpm;
  MhdRiemann::Solver riemann = MhdRiemann::Solver::kHlld;
  MhdEmfScheme emf_scheme = MhdEmfScheme::kUctHll;
  MhdTimeIntegrator time_integrator = MhdTimeIntegrator::kSsprk2;
  MhdDivergencePolicy divergence_policy = MhdDivergencePolicy::kStrictCt;
  double gamma_adiabatic = k_mhd_gamma_default;
  double cfl_ideal = 0.30;
  std::int32_t max_step_retries = 8;
  double fallback_fraction_limit = 0.001;
  double floor_energy_fraction_limit = 1.0e-8;
  // Audited legacy hydro floors, reused only for flux evaluation. Floor-added
  // energy into the conserved state is counted separately.
  double density_floor_code = 1.0e-10;
  double pressure_floor_code = 1.0e-10;
  bool deterministic = true;
  MhdStiffnessLimits stiffness{};

  [[nodiscard]] std::size_t ghostCells() const noexcept {
    return reconstruction == MhdReconstructionMode::kPpm ? 3U : 2U;
  }
};

// Gas patch state owned by the caller (ultimately core::SimulationState). The
// total-energy lane is thermal + kinetic + magnetic in MHD mode.
class MhdGasStateSoa final {
 public:
  MhdGasStateSoa() = default;
  explicit MhdGasStateSoa(std::size_t cell_count);
  void resize(std::size_t cell_count);
  [[nodiscard]] std::size_t size() const noexcept { return m_rho.size(); }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;

  [[nodiscard]] std::span<double> rho() noexcept { return m_rho.span(); }
  [[nodiscard]] std::span<const double> rho() const noexcept { return m_rho.span(); }
  [[nodiscard]] std::span<double> momentumX() noexcept { return m_mx.span(); }
  [[nodiscard]] std::span<const double> momentumX() const noexcept { return m_mx.span(); }
  [[nodiscard]] std::span<double> momentumY() noexcept { return m_my.span(); }
  [[nodiscard]] std::span<const double> momentumY() const noexcept { return m_my.span(); }
  [[nodiscard]] std::span<double> momentumZ() noexcept { return m_mz.span(); }
  [[nodiscard]] std::span<const double> momentumZ() const noexcept { return m_mz.span(); }
  [[nodiscard]] std::span<double> totalEnergy() noexcept { return m_energy.span(); }
  [[nodiscard]] std::span<const double> totalEnergy() const noexcept { return m_energy.span(); }
  [[nodiscard]] std::span<double> metalMass() noexcept { return m_metal.span(); }
  [[nodiscard]] std::span<const double> metalMass() const noexcept { return m_metal.span(); }

  [[nodiscard]] MhdConservedState load(std::size_t cell) const noexcept;
  void store(std::size_t cell, const MhdConservedState& state) noexcept;
  void copyFrom(const MhdGasStateSoa& other);
  void combineSsprk2(const MhdGasStateSoa& previous,
                     const MhdGasStateSoa& first_stage,
                     const MhdGasStateSoa& first_stage_rhs,
                     double dt_code);

 private:
  MhdDoubleArray m_rho;
  MhdDoubleArray m_mx;
  MhdDoubleArray m_my;
  MhdDoubleArray m_mz;
  MhdDoubleArray m_energy;
  MhdDoubleArray m_metal;
};

// Padded primitive cache. Ghost layers are filled from physical boundary kinds
// (periodic wrap, open zero-gradient, reflecting conducting wall). It is
// phase-local scratch and never an authority.
struct MhdPrimitivePatch {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  std::size_t ghost = 0;
  std::size_t px = 0;
  std::size_t py = 0;
  std::size_t pz = 0;

  MhdDoubleArray rho;
  MhdDoubleArray vx;
  MhdDoubleArray vy;
  MhdDoubleArray vz;
  MhdDoubleArray pressure;
  MhdDoubleArray bx;
  MhdDoubleArray by;
  MhdDoubleArray bz;
  MhdDoubleArray metal;

  void resize(std::size_t patch_nx, std::size_t patch_ny, std::size_t patch_nz,
              std::size_t ghost_cells);
  [[nodiscard]] std::size_t index(std::ptrdiff_t i, std::ptrdiff_t j,
                                  std::ptrdiff_t k) const noexcept {
    const auto pi = static_cast<std::size_t>(i + static_cast<std::ptrdiff_t>(ghost));
    const auto pj = static_cast<std::size_t>(j + static_cast<std::ptrdiff_t>(ghost));
    const auto pk = static_cast<std::size_t>(k + static_cast<std::ptrdiff_t>(ghost));
    return pi + px * (pj + py * pk);
  }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

// Patch-local phase workspace. One instance per active patch/task; edges and
// stage buffers are reused across SSPRK2/RKL/Hall subcycles.
struct MhdStageWorkspace {
  MhdPrimitivePatch primitives{};
  MhdFaceSignalPatch signals{};
  MhdEdgeEmfStorage emf_rate{};
  MhdGasStateSoa gas_rhs{};
  MhdDoubleArray phi_rhs_x{};
  MhdDoubleArray phi_rhs_y{};
  MhdDoubleArray phi_rhs_z{};
  std::vector<MhdCellPhysicalProperties> cell_properties{};
  std::vector<MhdTransportCoefficients> coefficient_cache{};
  MhdNonIdealWorkspace nonideal{};

  void resize(const MhdPatchGeometry& geometry,
              MhdReconstructionMode reconstruction,
              bool nonideal_enabled);
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

struct MhdAdvanceResult {
  bool accepted = false;
  double dt_used_code = 0.0;
  std::int32_t retry_count = 0;
  std::uint64_t face_count = 0U;
  std::uint64_t hlle_fallback_faces = 0U;
  std::uint64_t first_order_faces = 0U;
  double fallback_fraction = 0.0;
  double floor_added_energy = 0.0;
  double divergence_max_normalized = 0.0;
  double divergence_normalizer = 0.0;
  MhdConservationBudget budget{};
  MhdTransportStageDiagnostics transport{};
  std::uint64_t ideal_stage_evaluations = 0U;
  std::uint32_t accepted_stage_ordinal = 0U;
};

// The reference ideal + non-ideal stepped operator. All algorithms are selected
// before the face loop (no per-face virtual dispatch or allocation).
class MhdStepOperator {
 public:
  explicit MhdStepOperator(MhdNumericsPolicy policy = {});

  [[nodiscard]] const MhdNumericsPolicy& policy() const noexcept { return m_policy; }

  // Advance one accepted interval. On failure the caller's gas/face state is
  // restored to its pre-call value (the operator snapshots internally and never
  // leaves a partially committed interval). Throws MhdError on an unrecoverable
  // physics failure after the retry budget is exhausted.
  [[nodiscard]] MhdAdvanceResult advance(const MhdPatchGeometry& geometry,
                                         const MhdPatchBoundary& boundary,
                                         MhdGasStateSoa& gas,
                                         MhdFaceFluxMutableView faces,
                                         const MhdStageContext& context,
                                         const MhdCoefficientProvider* coefficient_provider,
                                         MhdDiagnosticSink* diagnostics = nullptr);

  // Ideal fast-magnetosonic CFL candidate for one patch.
  [[nodiscard]] double estimateIdealTimeStep(const MhdPatchGeometry& geometry,
                                             const MhdGasStateSoa& gas,
                                             const MhdFaceFluxView& faces,
                                             double cfl) const;

  // Primitive recovery used by the operator and by diagnostics.
  [[nodiscard]] bool recoverPatchPrimitives(const MhdPatchGeometry& geometry,
                                            const MhdGasStateSoa& gas,
                                            const MhdFaceFluxView& faces,
                                            double scale_factor,
                                            MhdStageWorkspace& workspace,
                                            std::size_t* nonpositive_cells = nullptr) const;

  // Sum the conservative budget of a patch (thermal/kinetic/magnetic split).
  [[nodiscard]] MhdConservationBudget patchBudget(const MhdPatchGeometry& geometry,
                                                  const MhdGasStateSoa& gas,
                                                  const MhdFaceFluxView& faces,
                                                  double scale_factor) const;

  // Fill the ghost primitive layers for a patch boundary. Exposed for tests.
  static void fillGhostPrimitives(const MhdPatchGeometry& geometry,
                                  const MhdPatchBoundary& boundary,
                                  MhdPrimitivePatch& primitives);

  // One-dimensional reconstruction helpers used by the face sweep and tested
  // directly. `c` is the cell value; the stencil runs along the face normal.
  [[nodiscard]] static std::pair<double, double> reconstructPlmMc(
      double c_minus2, double c_minus, double c, double c_plus, double c_plus2,
      bool face_is_upper) noexcept;
  [[nodiscard]] static std::pair<double, double> reconstructPpm(
      double c_minus2, double c_minus, double c, double c_plus, double c_plus2,
      bool face_is_upper) noexcept;

 private:
  void computeIdealRate(const MhdPatchGeometry& geometry,
                        const MhdPatchBoundary& boundary,
                        const MhdGasStateSoa& gas,
                        std::span<const double> phi_x,
                        std::span<const double> phi_y,
                        std::span<const double> phi_z,
                        const MhdStageContext& context,
                        MhdStageWorkspace& workspace,
                        std::span<double> phi_rhs_x,
                        std::span<double> phi_rhs_y,
                        std::span<double> phi_rhs_z,
                        std::uint64_t* hlle_faces,
                        std::uint64_t* first_order_faces) const;

  void applyExpansionSourceHalf(std::span<double> rho,
                                std::span<double> mx,
                                std::span<double> my,
                                std::span<double> mz,
                                std::span<double> energy,
                                const MhdPatchGeometry& geometry,
                                std::span<const double> phi_x,
                                std::span<const double> phi_y,
                                std::span<const double> phi_z,
                                double scale_factor,
                                double hubble_rate_code,
                                double dt_half_code) const;

  MhdNumericsPolicy m_policy;
};

}  // namespace cosmosim::mhd
