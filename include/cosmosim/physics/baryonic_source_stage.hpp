#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "cosmosim/physics/ism_dust_sgs.hpp"
#include "cosmosim/physics/source_interval.hpp"
#include "cosmosim/physics/thermochemistry_operator.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

// ===========================================================================
// CHUI Part 2 local source-stage atomicity (corrective gate).
//
// One stage owns, for the active rows of one accepted interval:
//   ISM regime decision (+ dwell)  -> chemistry-fidelity decision (+ dwell)
//   -> thermochemistry (skipped for effective-regime rows)
//   -> live dust (accretion / sputtering)  -> SGS dissipation.
// prepare(): reads ONLY the step-start state, evaluates every decision and
//   update on a private one-row scratch copy of the composition, and stages
//   per row: all composition lanes, e, mu, both tags, both dwell counters and
//   the pending diagnostics. Throws SourceIntervalIncompleteError (or another
//   exception) without touching SimulationState.
// commit(): writes every staged row, re-derives the Z mirror, publishes the
//   pending diagnostics/ledgers and returns the rows whose P/T/c_s must be
//   refreshed. It is only called after the StageTransaction collective vote.
// discard(): drops staging and pending diagnostics (rejected interval).
// Memory: (lanes + 16) doubles per active row, admitted by the caller through
// the StageTransaction before prepare(); no full-state copy.
// Scope: this is LOCAL source-stage atomicity. A rejection after other
// stages of the same step have committed (hydro, AMR sweep) is the global
// distributed step-rollback problem owned by Part 3.
// ===========================================================================
struct BaryonicSourceStageConfig {
  std::optional<IsmRegimeConfig> regime;
  std::optional<DustConfig> dust;
  std::optional<SgsConfig> sgs;
};

struct BaryonicSourceStageReport {
  SourceLedger chemistry{};
  ChemistryShadowReport shadow{};
  ChemistryAdaptiveReport adaptive{};
  IsmRegimeReport regime{};
  DustUpdateLedger dust{};
  SgsLedger sgs{};
  std::uint64_t skipped_effective_rows = 0U;
};

class BaryonicSourceStage {
 public:
  BaryonicSourceStage(const ThermochemistryOperator& chemistry, BaryonicSourceStageConfig config);

  [[nodiscard]] static std::uint64_t stagingBytes(std::size_t rows, std::size_t lanes) noexcept {
    return static_cast<std::uint64_t>(rows) *
        (static_cast<std::uint64_t>(lanes) + 16U) * sizeof(double);
  }

  // rows empty => all gas rows.
  void prepare(const core::SimulationState& state, const SourceInterval& interval,
               std::span<const std::uint32_t> rows);
  std::vector<std::uint32_t> commit(core::SimulationState& state);
  void discard() noexcept;

  [[nodiscard]] bool prepared() const noexcept { return m_prepared; }
  [[nodiscard]] std::size_t sparseRowsToActivate() const noexcept {
    return m_sparse_rows_to_activate;
  }
  // Published only by commit(); unchanged by a rejected interval.
  [[nodiscard]] const BaryonicSourceStageReport& lastCommittedReport() const noexcept { return m_committed; }

 private:
  struct StagedMeta {
    std::uint32_t row = 0;
    std::uint8_t ism_tag = 0;
    std::uint8_t chemistry_tag = 0;
    bool chemistry_ran = false;
    bool refresh = false;
    double ism_dwell_years = 0.0;
    double chemistry_dwell_years = 0.0;
    double internal_energy_code = 0.0;
    double mean_molecular_weight = 0.0;
    double last_temperature_k = 0.0;
    double last_cooling_time_seconds = 0.0;
    double last_error_estimate = 0.0;
    double last_element_fraction_residual = 0.0;
    double last_charge_residual = 0.0;
    double cumulative_radiative_loss_erg = 0.0;
    double cumulative_heating_gain_erg = 0.0;
  };

  const ThermochemistryOperator& m_chemistry;
  BaryonicSourceStageConfig m_config;
  std::vector<StagedMeta> m_meta;
  std::vector<double> m_lanes;  // row-major: staged row x lane
  std::size_t m_lane_count = 0;
  std::size_t m_sparse_rows_to_activate = 0U;
  BaryonicSourceStageReport m_pending{};
  BaryonicSourceStageReport m_committed{};
  bool m_prepared = false;
};

}  // namespace cosmosim::physics
