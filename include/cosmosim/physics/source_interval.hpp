#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/stage_transaction.hpp"

namespace cosmosim::physics {

// ---------------------------------------------------------------------------
// CHUI Part 2 single source-term contract shared by thermochemistry, dust,
// ISM/SGS, star formation, stellar evolution and feedback.
//
// Equations of the time mapping
// -----------------------------
// Source rates are physical (proper) quantities. The accepted integrator
// interval [t_begin, t_end] in code time maps to the physical elapsed
// interval
//   isolated frames (a == 1):       dt_phys = (t_end - t_begin) * s_per_code
//   cosmological frames (comoving): dt_phys = t_cosmic(a_end) - t_cosmic(a_begin)
// using the background's exact cosmic-time integral when available (the
// SourceRuntime already owns that conversion for stellar ages). Proper
// densities are rho_phys = rho_code / a^3 (comoving) or rho_code (isolated)
// evaluated at the source epoch, converted with UnitSystem::densityCodeToCgs.
// Proper specific internal energy equals the code value in both frames: the
// code's e is the peculiar/thermal specific energy, so u_phys = e * v_unit^2
// (no factor of a; Hubble cooling of e is the hydro/gravity owner's term and
// is NEVER re-applied here).
// ---------------------------------------------------------------------------
struct SourceInterval {
  double time_begin_code = 0.0;
  double time_end_code = 0.0;
  double scale_factor_begin = 1.0;
  double scale_factor_end = 1.0;
  // Physical elapsed time of the interval in seconds (> 0 for a real step).
  double physical_seconds = 0.0;
  // Epoch at which densities/radiation are evaluated (Lie split: step end).
  double evaluation_scale_factor = 1.0;
  double evaluation_redshift = 0.0;
  bool cosmological = false;
  std::uint64_t accepted_step_index = 0U;

  [[nodiscard]] bool valid() const noexcept;
};

enum class SourceAdvanceStatus : std::uint8_t {
  // The ENTIRE accepted interval was integrated with bounded error.
  kCompleted = 0,
  // A local solve could not finish: the stage must be retried with the
  // suggested shorter interval. Nothing was committed.
  kRetryShorterInterval = 1,
  // The requested model cannot be integrated (missing data, forbidden
  // extrapolation, invalid state). Nothing was committed.
  kRejectedInvalid = 2,
};

[[nodiscard]] std::string_view sourceAdvanceStatusLabel(SourceAdvanceStatus status) noexcept;

// Per-operation ledger. Energies in erg (physical), masses in code mass.
struct SourceLedger {
  double radiative_loss_erg = 0.0;          // energy removed from gas by radiation
  double heating_gain_erg = 0.0;            // photo/photoelectric/other heating
  double chemical_binding_delta_erg = 0.0;  // change of chemical potential energy
  double floor_injection_erg = 0.0;         // energy added by temperature floors
  double thermal_energy_delta_erg = 0.0;    // net change of thermal energy
  double element_mass_residual_code = 0.0;  // max |sum element after - before|
  double charge_residual = 0.0;             // max |sum charge| per H nucleus
  double dust_to_gas_transfer_code = 0.0;   // + = solids returned to gas
  double gas_to_dust_transfer_code = 0.0;   // + = gas locked into solids
  std::uint64_t cells = 0U;
  std::uint64_t accepted_substeps = 0U;
  std::uint64_t rejected_substeps = 0U;
  std::uint64_t newton_iterations = 0U;
  std::uint64_t floor_activations = 0U;
  // Accepted substeps whose metal-table evaluation was clamped in n_H or z
  // (only possible under metal_cooling_table_range = clamp_density_redshift).
  std::uint64_t metal_table_clamps = 0U;
  double max_error_estimate = 0.0;
  double max_newton_residual = 0.0;

  SourceLedger& operator+=(const SourceLedger& rhs);
};

struct SourceAdvanceResult {
  SourceAdvanceStatus status = SourceAdvanceStatus::kCompleted;
  // Physical seconds actually integrated. Equals the request on success.
  double completed_physical_seconds = 0.0;
  double requested_physical_seconds = 0.0;
  // Largest interval (physical seconds) the failing cell could integrate;
  // a legal shorter-step request for the time coordinator.
  double suggested_retry_physical_seconds = 0.0;
  std::string diagnostic;
  SourceLedger ledger;
  std::vector<std::uint32_t> touched_rows;

  [[nodiscard]] bool completed() const noexcept { return status == SourceAdvanceStatus::kCompleted; }
};

// Exact interval closure test. The integrator accumulates time in long
// double and the remainder must vanish to a scale-aware tolerance that cannot
// hide a finite unprocessed remainder: |t_done - t_req| <= 4 eps t_req.
[[nodiscard]] bool sourceIntervalClosed(long double completed_seconds,
                                        long double requested_seconds) noexcept;

// Thrown inside StageTransaction::prepare when a source cannot complete the
// interval (maps to StageRejectionKind::kNumericalRetry). Carries the legal
// shorter-interval request.
class SourceIntervalIncompleteError : public core::StageNumericalRetryRequest {
 public:
  SourceIntervalIncompleteError(const std::string& message, double suggested_retry_physical_seconds)
      : core::StageNumericalRetryRequest(message),
        m_suggested_retry_physical_seconds(suggested_retry_physical_seconds) {}
  [[nodiscard]] double suggestedRetryPhysicalSeconds() const noexcept {
    return m_suggested_retry_physical_seconds;
  }

 private:
  double m_suggested_retry_physical_seconds = 0.0;
};

// Thrown when a selected production model cannot run (missing data file,
// failed checksum, forbidden extrapolation, unsupported topology). It is
// NEVER converted into a cheaper model silently.
class SourceModelUnavailableError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

}  // namespace cosmosim::physics
