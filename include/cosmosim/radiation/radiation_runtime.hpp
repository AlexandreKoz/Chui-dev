#pragma once

// CHUI-RT-BH-AGN-REF-v1 radiation runtime (RT01-RT07 composition):
// canonical state allocation, two-stage PD-ARS IMEX subcycling over one
// accepted source interval, time-integrated emission, cosmological edge
// transfer, staged atomic radiation+matter commit and ledgers.
//
// IMEX (He, Wibking & Krumholz 2024, MNRAS 531, 1228, eqs. 28-29; first
// stage evaluated at the END time t_{n+1}, both stages implicit in S):
//   U*      = U^n + dt T(U^n, t_n)     + dt   S(U*, t_{n+1})
//   U^{n+1} = U^n + dt/2 [T(U^n,t_n) + T(U*,t_{n+1})] + dt/2 [S(U*) + S(U^{n+1})]
// S(U*) is recovered exactly as (U* - U_rhs)/dt from the stage-1 solve (no
// re-evaluation). Emission j is part of T, so the PD-ARS weights integrate a
// constant luminosity exactly once over the substep (never full luminosity
// in both stages). Streaming is second order; asymptotic diffusion is
// first-order in time — no global second-order claim. The explicit
// reference integrator (kExplicitReference) uses U^{n+1} = U^n + dt T(U^n)
// followed by the implicit S(dt) solve (split comparator only).
// Transport CFL always applies (dt_rt = cfl / sum_d c_hat/dx_phys,d); the
// last substep is shortened to land exactly on the interval end; rejected
// substeps roll back and retry with dt/2 up to max_retries; exhaustion
// throws with cell/regime information (no partial commit).
// Commit: the radiation lanes and the matter exchange (gas total energy and
// momentum through SourceTransactionBatch, channel kRtMatterExchange) are
// staged together and applied only after the whole interval succeeded.

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "cosmosim/core/config_rt_bh_agn.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/rt_bh_agn_frame.hpp"
#include "cosmosim/radiation/radiation_matter_coupling.hpp"
#include "cosmosim/radiation/radiation_spectrum.hpp"
#include "cosmosim/radiation/radiation_transport.hpp"

namespace cosmosim::radiation {

enum class RadiationCouplingModel : std::uint8_t { kNone = 0, kGreyThermalReference = 1, kPhotoionization = 2 };

struct RadiationRuntimeSettings {
  core::RadiationTransportConfig config{};
  RadiationSpectralDefinition spectrum{};
  RadiationTransportSettings transport{};
  RadiationCouplingConstants coupling{};
  RadiationCouplingModel coupling_model = RadiationCouplingModel::kNone;
};

// Time-integrated source deposit: comoving energy-density rate added to one
// record cell and group over the whole interval (isotropic: no net flux).
struct RadiationSourceDeposit {
  std::uint32_t record = 0U;
  std::uint32_t cell = 0U;
  std::uint32_t group = 0U;
  double comoving_energy_rate = 0.0;
};

struct RadiationIntervalInputs {
  double dt_proper_code = 0.0;
  double scale_factor_begin = 1.0;
  double scale_factor_end = 1.0;
  double hubble_begin_code = 0.0;
  double hubble_end_code = 0.0;
  std::uint64_t interval_index = 0U;
  physics::ProperFrameContext frame_end{};
  std::vector<RadiationSourceDeposit> sources;
};

struct RadiationIntervalReport {
  std::uint32_t substeps = 0U;
  std::uint32_t retries = 0U;
  std::uint64_t newton_iterations = 0U;
  std::uint64_t first_order_fallbacks = 0U;
  std::uint64_t hll_fallback_faces = 0U;
  double emitted_comoving_energy = 0.0;
  double boundary_outflow_comoving_energy = 0.0;
  double redshift_leaked_comoving_energy = 0.0;
  double gas_energy_exchange_code = 0.0;          // proper energy given to gas (sum)
  double radiation_energy_change_comoving = 0.0;  // sum E_comov V_comov change
  double rsla_invariant_residual = 0.0;           // |d(E_gas + (c/c_hat) E_rad) - sources + outflow| / scale
  std::vector<std::uint32_t> touched_gas_rows;
};

// Allocate zero (vacuum) moments on the owned leaf patches of `state`.
void allocateRadiationState(core::SimulationState& state, const RadiationSpectralDefinition& spectrum,
                            double c_hat_over_c, double scale_factor, std::uint32_t world_rank);
// Verify the radiation records still match the leaf patch tiling (epoch).
[[nodiscard]] bool radiationRecordsMatchPatches(const core::SimulationState& state, std::uint32_t world_rank,
                                                std::string* reason = nullptr);
// Predicted canonical bytes: 4 lanes * G * 8 * N_cells (PERF01 scale law).
[[nodiscard]] std::uint64_t radiationCanonicalBytes(std::uint64_t cells, std::uint32_t groups);

class RadiationRuntime {
 public:
  explicit RadiationRuntime(RadiationRuntimeSettings settings);
  [[nodiscard]] const RadiationRuntimeSettings& settings() const noexcept { return m_settings; }
  [[nodiscard]] RadiationIntervalReport advanceInterval(core::SimulationState& state,
                                                        const RadiationIntervalInputs& inputs) const;

 private:
  RadiationRuntimeSettings m_settings;
};

}  // namespace cosmosim::radiation
