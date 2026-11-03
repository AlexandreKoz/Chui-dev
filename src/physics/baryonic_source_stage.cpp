#include "cosmosim/physics/baryonic_source_stage.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

BaryonicSourceStage::BaryonicSourceStage(const ThermochemistryOperator& chemistry, BaryonicSourceStageConfig config)
    : m_chemistry(chemistry), m_config(config) {}

void BaryonicSourceStage::discard() noexcept {
  m_meta.clear();
  m_lanes.clear();
  m_pending = BaryonicSourceStageReport{};
  m_prepared = false;
}

void BaryonicSourceStage::prepare(const core::SimulationState& state, const SourceInterval& interval,
                                  std::span<const std::uint32_t> rows) {
  discard();
  if (!interval.valid()) throw std::invalid_argument("baryonic source stage: invalid source interval");
  const auto& sidecar = state.gas_composition;
  const CompositionLayout& layout = m_chemistry.layout();
  if (!sidecar.enabled() || !(sidecar.descriptor() == layout.descriptor())) {
    throw SourceModelUnavailableError("baryonic source stage requires the gas composition sidecar with layout '" +
                                      layout.descriptor().layout_name + "'");
  }
  sidecar.requireAlignedWith(state.gas_cells.gas_cell_id, "BaryonicSourceStage::prepare");
  const ThermochemistryUnitBridge& units = m_chemistry.config().units;
  m_lane_count = sidecar.laneCount();
  const std::size_t row_count = rows.empty() ? state.cells.size() : rows.size();
  m_meta.reserve(row_count);
  m_lanes.reserve(row_count * m_lane_count);
  core::GasCompositionSidecar scratch = sidecar.emptyLike(1U);  // one row of private scratch
  const double a = interval.evaluation_scale_factor;
  const double time_code_s = units.length_code_to_cm / units.velocity_code_to_cm_s;
  const double dt_code_time = interval.physical_seconds / time_code_s;
  double min_retry = std::numeric_limits<double>::infinity();
  std::string first_failure;
  std::size_t failures = 0;

  for (std::size_t i = 0; i < row_count; ++i) {
    const std::uint32_t row = rows.empty() ? static_cast<std::uint32_t>(i) : rows[i];
    if (row >= state.cells.size()) throw std::out_of_range("baryonic source stage: active row outside gas state");
    StagedMeta meta;
    meta.row = row;
    meta.internal_energy_code = state.gas_cells.internal_energy_code[row];
    meta.mean_molecular_weight = sidecar.mean_molecular_weight[row];
    // 1) ISM regime (provisional).
    meta.ism_tag = sidecar.ism_regime_tag[row];
    meta.ism_dwell_years = sidecar.ism_regime_dwell_years[row];
    if (m_config.regime.has_value()) {
      const IsmRegimeDecision d = decideIsmRegime(state, *m_config.regime, units, interval, row);
      meta.ism_tag = d.tag;
      meta.ism_dwell_years = d.dwell_years;
      if (d.held_by_dwell) ++m_pending.regime.held_by_dwell;
      if (d.changed) {
        m_pending.regime.changed_rows.push_back(row);
        meta.refresh = true;  // closure ownership changed: re-derive P/c_s
        if (d.tag == core::k_gas_ism_regime_effective) ++m_pending.regime.to_effective;
        else ++m_pending.regime.to_resolved;
      }
      if (d.tag == core::k_gas_ism_regime_effective) ++m_pending.regime.effective_cells;
      else ++m_pending.regime.resolved_cells;
    }
    // 2) Chemistry fidelity (provisional).
    const ChemistryFidelityDecision fidelity = m_chemistry.decideChemistryFidelity(state, interval, row);
    meta.chemistry_tag = fidelity.tag;
    meta.chemistry_dwell_years = fidelity.dwell_years;
    if (fidelity.promoted) ++m_pending.adaptive.promoted;
    if (fidelity.demoted) ++m_pending.adaptive.demoted;
    if (fidelity.held_by_dwell) ++m_pending.adaptive.held_by_dwell;
    if (fidelity.tag == k_chemistry_tag_detailed) ++m_pending.adaptive.detailed_cells;
    else ++m_pending.adaptive.hybrid_cells;
    // Private copy of the row's composition.
    for (std::size_t l = 0; l < m_lane_count; ++l) scratch.setValue(l, 0, sidecar.value(l, row));
    double temperature_k = state.gas_cells.temperature_code[row];
    // 3) Thermochemistry (resolved regime only — no double cooling).
    if (meta.ism_tag == core::k_gas_ism_regime_effective) {
      ++m_pending.skipped_effective_rows;
    } else {
      const ThermochemistryRowResult r =
          m_chemistry.advanceRow(state, interval, row, meta.chemistry_tag, m_pending.shadow);
      if (r.status != SourceAdvanceStatus::kCompleted) {
        ++failures;
        if (first_failure.empty()) first_failure = r.diagnostic;
        min_retry = std::min(min_retry, r.suggested_retry_seconds);
        continue;
      }
      m_pending.chemistry += r.ledger;
      meta.chemistry_ran = true;
      meta.refresh = true;
      meta.internal_energy_code = r.staged.internal_energy_code;
      meta.mean_molecular_weight = r.staged.mean_molecular_weight;
      temperature_k = r.staged.temperature_k;
      for (std::size_t s = 0; s < m_chemistry.solver().speciesCount(); ++s) {
        scratch.setValue(layout.speciesLane(s), 0, r.staged.species_mass_code[s]);
      }
    }
    // 4) Live dust on the private copy.
    if (m_config.dust.has_value() && scratch.descriptor().dust_lane_count > 0U) {
      const double rho_cgs = state.gas_cells.density_code[row] * (units.comoving ? 1.0 / (a * a * a) : 1.0) *
          units.density_code_to_cgs;
      const double mass = compositionGasMass(scratch, 0);
      const double x_h = mass > 0.0 ? scratch.value(layout.elementLane(ChemicalElement::kH), 0) / mass : 0.0;
      const double n_h = x_h * rho_cgs / (chemicalElementInfo(ChemicalElement::kH).atomic_mass_amu * k_amu_g);
      advanceDustCell(scratch, layout, 0, n_h, temperature_k, interval.physical_seconds, *m_config.dust,
                      m_pending.dust);
    }
    // 5) SGS dissipation (exact transfer SGS -> thermal).
    if (m_config.sgs.has_value() && layout.hasSgs()) {
      const std::uint32_t lane = scratch.descriptor().sgs_energy_lane;
      const double mass = state.cells.mass_code[row];
      const double volume = mass / std::max(state.gas_cells.density_code[row], 1.0e-300);
      const double dx = std::cbrt(volume) * (units.comoving ? a : 1.0);
      const SgsDissipation d = sgsDissipate(scratch.value(lane, 0), mass, dt_code_time, dx, *m_config.sgs);
      scratch.setValue(lane, 0, d.sgs_energy);
      if (d.dissipated > 0.0 && mass > 0.0) {
        meta.internal_energy_code += d.dissipated / mass;
        m_pending.sgs.dissipated_to_heat_code += d.dissipated;
        meta.refresh = true;
      }
    }
    m_meta.push_back(meta);
    for (std::size_t l = 0; l < m_lane_count; ++l) m_lanes.push_back(scratch.value(l, 0));
  }
  if (failures > 0U) {
    discard();
    throw SourceIntervalIncompleteError(
        "baryonic source stage could not integrate the accepted interval in " + std::to_string(failures) +
            " cell(s); first: " + first_failure,
        std::isfinite(min_retry) ? min_retry : 0.0);
  }
  m_prepared = true;
}

std::vector<std::uint32_t> BaryonicSourceStage::commit(core::SimulationState& state) {
  if (!m_prepared) throw std::logic_error("BaryonicSourceStage::commit without a successful prepare");
  auto& sidecar = state.gas_composition;
  std::vector<std::uint32_t> refresh;
  std::vector<std::uint32_t> rows;
  refresh.reserve(m_meta.size());
  rows.reserve(m_meta.size());
  for (std::size_t i = 0; i < m_meta.size(); ++i) {
    const StagedMeta& meta = m_meta[i];
    for (std::size_t l = 0; l < m_lane_count; ++l) sidecar.setValue(l, meta.row, m_lanes[i * m_lane_count + l]);
    state.gas_cells.internal_energy_code[meta.row] = meta.internal_energy_code;
    sidecar.mean_molecular_weight[meta.row] = meta.mean_molecular_weight;
    sidecar.ism_regime_tag[meta.row] = meta.ism_tag;
    sidecar.ism_regime_dwell_years[meta.row] = meta.ism_dwell_years;
    sidecar.chemistry_model_tag[meta.row] = meta.chemistry_tag;
    sidecar.model_dwell_years[meta.row] = meta.chemistry_dwell_years;
    rows.push_back(meta.row);
    if (meta.refresh) refresh.push_back(meta.row);
  }
  if (m_config.dust.has_value()) (void)synchronizeMetalMirrorFromComposition(state, rows);
  m_committed = std::move(m_pending);
  m_meta.clear();
  m_lanes.clear();
  m_pending = BaryonicSourceStageReport{};
  m_prepared = false;
  return refresh;
}

}  // namespace cosmosim::physics
