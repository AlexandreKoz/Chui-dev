#include "cosmosim/physics/source_transaction.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

const char* sourceChannelLabel(SourceChannel channel) noexcept {
  switch (channel) {
    case SourceChannel::kBhAccretionWithdrawal: return "bh_accretion_withdrawal";
    case SourceChannel::kBhSeedConversion: return "bh_seed_conversion";
    case SourceChannel::kBhMerger: return "bh_merger";
    case SourceChannel::kBhDynamicalFrictionReaction: return "bh_dynamical_friction_reaction";
    case SourceChannel::kAgnThermal: return "agn_thermal";
    case SourceChannel::kAgnKinetic: return "agn_kinetic";
    case SourceChannel::kAgnWind: return "agn_wind";
    case SourceChannel::kAgnJet: return "agn_jet";
    case SourceChannel::kAgnMagnetic: return "agn_magnetic";
    case SourceChannel::kAgnNozzle: return "agn_nozzle";
    case SourceChannel::kRtMatterExchange: return "rt_matter_exchange";
  }
  return "unknown";
}

ConservedIncrement& ConservedIncrement::operator+=(const ConservedIncrement& rhs) {
  mass_code += rhs.mass_code;
  for (std::size_t a = 0; a < 3; ++a) {
    momentum_code[a] += rhs.momentum_code[a];
    angular_momentum_code[a] += rhs.angular_momentum_code[a];
  }
  total_energy_code += rhs.total_energy_code;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) element_mass_code[e] += rhs.element_mass_code[e];
  metal_mass_code += rhs.metal_mass_code;
  return *this;
}

ConservedIncrement ConservedIncrement::negated() const {
  ConservedIncrement out;
  out.mass_code = -mass_code;
  for (std::size_t a = 0; a < 3; ++a) {
    out.momentum_code[a] = -momentum_code[a];
    out.angular_momentum_code[a] = -angular_momentum_code[a];
  }
  out.total_energy_code = -total_energy_code;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) out.element_mass_code[e] = -element_mass_code[e];
  out.metal_mass_code = -metal_mass_code;
  return out;
}

namespace {

[[nodiscard]] bool finiteIncrement(const ConservedIncrement& d) noexcept {
  if (!std::isfinite(d.mass_code) || !std::isfinite(d.total_energy_code) || !std::isfinite(d.metal_mass_code)) {
    return false;
  }
  for (std::size_t a = 0; a < 3; ++a) {
    if (!std::isfinite(d.momentum_code[a]) || !std::isfinite(d.angular_momentum_code[a])) return false;
  }
  for (double v : d.element_mass_code) {
    if (!std::isfinite(v)) return false;
  }
  return true;
}

// Relative closure residual of a scalar ledger: |sum| / max(sum |terms|, floor).
[[nodiscard]] double closureResidual(double a, double b, double c) noexcept {
  const double scale = std::max(std::abs(a) + std::abs(b) + std::abs(c), 1.0e-300);
  return std::abs(a + b - c) / scale;
}

}  // namespace

SourceTransactionBatch::SourceTransactionBatch(TransactionFloors floors) : m_floors(floors) {}

void SourceTransactionBatch::add(SourceTransactionPlan plan) {
  if (m_phase != TransactionPhase::kPlanned) {
    throw std::logic_error("SourceTransactionBatch: plans may only be added before prepare()");
  }
  m_plans.push_back(std::move(plan));
}

void SourceTransactionBatch::reject(std::string reason) {
  m_phase = TransactionPhase::kRejected;
  m_report.phase = m_phase;
  m_report.failure = std::move(reason);
}

bool SourceTransactionBatch::prepare(const core::SimulationState& state) {
  if (m_phase != TransactionPhase::kPlanned) {
    throw std::logic_error("SourceTransactionBatch: prepare() called twice");
  }
  // Deterministic plan order: by transaction identity.
  std::sort(m_plans.begin(), m_plans.end(),
            [](const SourceTransactionPlan& a, const SourceTransactionPlan& b) { return a.id < b.id; });
  for (std::size_t i = 1; i < m_plans.size(); ++i) {
    if (m_plans[i].id == m_plans[i - 1].id) {
      reject("duplicate source transaction identity (event replay)");
      return false;
    }
  }
  const std::uint64_t epoch = state.cellIndexGeneration();
  std::vector<RowAggregate> rows;
  for (const SourceTransactionPlan& plan : m_plans) {
    if (plan.topology_epoch != epoch) {
      reject("source transaction topology epoch is stale (planned against a different cell index generation)");
      return false;
    }
    if (!finiteIncrement(plan.source_increment) || !finiteIncrement(plan.external_increment)) {
      reject("source transaction carries a non-finite source/external increment");
      return false;
    }
    for (const GasCellIncrement& cell : plan.cells) {
      if (cell.local_row >= state.cells.size() || cell.local_row >= state.gas_cells.size()) {
        reject("source transaction recipient row is out of range");
        return false;
      }
      if (state.gas_cells.gas_cell_id[cell.local_row] != cell.gas_cell_id || cell.gas_cell_id == 0U) {
        reject("source transaction recipient identity does not match the stable gas_cell_id");
        return false;
      }
      if (!finiteIncrement(cell.delta) || !(cell.withdrawal_fraction >= 0.0) || cell.withdrawal_fraction > 1.0) {
        reject("source transaction carries a non-finite delta or invalid withdrawal fraction");
        return false;
      }
      rows.push_back(RowAggregate{cell.gas_cell_id, cell.local_row, cell.delta, 1.0 - cell.withdrawal_fraction});
    }
  }
  std::sort(rows.begin(), rows.end(),
            [](const RowAggregate& a, const RowAggregate& b) { return a.gas_cell_id < b.gas_cell_id; });
  m_rows.clear();
  for (const RowAggregate& r : rows) {
    if (!m_rows.empty() && m_rows.back().gas_cell_id == r.gas_cell_id) {
      m_rows.back().delta += r.delta;
      // Withdrawal fractions are all evaluated against the same pre-state,
      // so they add: retained = 1 - sum f.
      m_rows.back().retained_fraction -= (1.0 - r.retained_fraction);
    } else {
      m_rows.push_back(r);
    }
  }
  // Admissibility on the aggregated image (no mutation).
  for (const RowAggregate& r : m_rows) {
    const std::uint32_t row = r.row;
    const double m = state.cells.mass_code[row];
    const double vx = state.gas_cells.velocity_x_peculiar[row];
    const double vy = state.gas_cells.velocity_y_peculiar[row];
    const double vz = state.gas_cells.velocity_z_peculiar[row];
    const double u = state.gas_cells.internal_energy_code[row];
    if (!(r.retained_fraction >= 0.0)) {
      std::ostringstream msg;
      msg << "overlapping withdrawals exceed the available mass of gas_cell_id " << r.gas_cell_id;
      reject(msg.str());
      return false;
    }
    const double m_new = m + r.delta.mass_code;
    if (!(m_new >= m_floors.mass_floor_code) || !std::isfinite(m_new)) {
      std::ostringstream msg;
      msg << "transaction would violate the mass floor of gas_cell_id " << r.gas_cell_id;
      reject(msg.str());
      return false;
    }
    const double px = m * vx + r.delta.momentum_code[0];
    const double py = m * vy + r.delta.momentum_code[1];
    const double pz = m * vz + r.delta.momentum_code[2];
    const double e_total = m * (u + 0.5 * (vx * vx + vy * vy + vz * vz)) + r.delta.total_energy_code;
    const double u_new = (e_total - (px * px + py * py + pz * pz) / (2.0 * m_new)) / m_new;
    if (!std::isfinite(u_new) || u_new < m_floors.specific_internal_energy_floor_code) {
      std::ostringstream msg;
      msg << "transaction would produce an inadmissible specific internal energy (" << u_new
          << ") in gas_cell_id " << r.gas_cell_id;
      reject(msg.str());
      return false;
    }
  }
  m_phase = TransactionPhase::kPrepared;
  m_report.phase = m_phase;
  m_report.plans = m_plans.size();
  return true;
}

bool SourceTransactionBatch::validate() {
  if (m_phase != TransactionPhase::kPrepared) {
    throw std::logic_error("SourceTransactionBatch: validate() requires a prepared batch");
  }
  double worst = 0.0;
  for (const SourceTransactionPlan& plan : m_plans) {
    ConservedIncrement gas{};
    for (const GasCellIncrement& c : plan.cells) gas += c.delta;
    double r = closureResidual(gas.mass_code, plan.source_increment.mass_code, plan.external_increment.mass_code);
    for (std::size_t a = 0; a < 3; ++a) {
      r = std::max(r, closureResidual(gas.momentum_code[a], plan.source_increment.momentum_code[a],
                                      plan.external_increment.momentum_code[a]));
    }
    r = std::max(r, closureResidual(gas.total_energy_code, plan.source_increment.total_energy_code,
                                    plan.external_increment.total_energy_code));
    worst = std::max(worst, r);
    if (r > plan.relative_tolerance) {
      std::ostringstream msg;
      msg << "source transaction ledger does not close (channel " << sourceChannelLabel(plan.id.channel)
          << ", source " << plan.id.source_id << ", residual " << r << ")";
      reject(msg.str());
      m_report.max_relative_closure_residual = worst;
      return false;
    }
  }
  m_report.max_relative_closure_residual = worst;
  m_phase = TransactionPhase::kValidated;
  m_report.phase = m_phase;
  return true;
}

bool SourceTransactionBatch::authorize(bool participants_agree) {
  if (m_phase != TransactionPhase::kValidated && m_phase != TransactionPhase::kRejected) {
    throw std::logic_error("SourceTransactionBatch: authorize() requires a validated batch");
  }
  if (m_phase == TransactionPhase::kRejected) return false;
  if (!participants_agree) {
    reject("a participant rank rejected the distributed source transaction");
    return false;
  }
  m_phase = TransactionPhase::kAuthorized;
  m_report.phase = m_phase;
  return true;
}

std::vector<std::uint32_t> SourceTransactionBatch::commit(core::SimulationState& state) {
  if (m_phase != TransactionPhase::kAuthorized) {
    throw std::logic_error("SourceTransactionBatch: commit() requires an authorized, not yet committed batch");
  }
  // Mark committed first: a second commit can never replay this batch even if
  // an internal error escapes below (the workflow then aborts to checkpoint).
  m_phase = TransactionPhase::kCommitted;
  m_report.phase = m_phase;
  const bool composition = state.gas_composition.enabled();
  std::optional<CompositionLayout> layout;
  if (composition) layout = CompositionLayout::fromDescriptor(state.gas_composition.descriptor());
  std::vector<std::uint32_t> touched;
  touched.reserve(m_rows.size());
  for (const RowAggregate& r : m_rows) {
    const std::uint32_t row = r.row;
    const double m = state.cells.mass_code[row];
    const double vx = state.gas_cells.velocity_x_peculiar[row];
    const double vy = state.gas_cells.velocity_y_peculiar[row];
    const double vz = state.gas_cells.velocity_z_peculiar[row];
    const double u = state.gas_cells.internal_energy_code[row];
    const double m_new = m + r.delta.mass_code;
    const double px = m * vx + r.delta.momentum_code[0];
    const double py = m * vy + r.delta.momentum_code[1];
    const double pz = m * vz + r.delta.momentum_code[2];
    const double e_total = m * (u + 0.5 * (vx * vx + vy * vy + vz * vz)) + r.delta.total_energy_code;
    const double u_new = (e_total - (px * px + py * py + pz * pz) / (2.0 * m_new)) / m_new;
    state.gas_cells.velocity_x_peculiar[row] = px / m_new;
    state.gas_cells.velocity_y_peculiar[row] = py / m_new;
    state.gas_cells.velocity_z_peculiar[row] = pz / m_new;
    state.gas_cells.internal_energy_code[row] = u_new;
    state.gas_cells.density_code[row] *= m_new / m;
    state.cells.mass_code[row] = m_new;
    if (composition) {
      if (r.retained_fraction < 1.0) state.gas_composition.scaleRow(row, std::clamp(r.retained_fraction, 0.0, 1.0));
      for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
        if (r.delta.element_mass_code[e] > 0.0) {
          state.gas_composition.addValue(layout->elementLane(static_cast<ChemicalElement>(e)), row,
                                         r.delta.element_mass_code[e]);
        }
      }
      // Injected material enters as neutral atoms (convention of the stellar
      // mechanical feedback commit); thermochemistry re-solves the species.
      const std::size_t species = state.gas_composition.descriptor().species_lane_count;
      if (species >= k_hybrid_species_count) {
        if (r.delta.element_mass_code[0] > 0.0) {
          state.gas_composition.addValue(layout->speciesLane(static_cast<std::size_t>(ChemicalSpecies::kHI)), row,
                                         r.delta.element_mass_code[0]);
        }
        if (r.delta.element_mass_code[1] > 0.0) {
          state.gas_composition.addValue(layout->speciesLane(static_cast<std::size_t>(ChemicalSpecies::kHeI)), row,
                                         r.delta.element_mass_code[1]);
        }
      }
    } else {
      const double metal_before = state.gas_cells.metal_mass_code[row];
      double metal = metal_before * std::clamp(r.retained_fraction, 0.0, 1.0);
      if (r.delta.metal_mass_code > 0.0) metal += r.delta.metal_mass_code;
      state.gas_cells.metal_mass_code[row] = std::max(metal, 0.0);
    }
    m_report.committed_to_gas += r.delta;
    touched.push_back(row);
  }
  for (const SourceTransactionPlan& plan : m_plans) m_report.committed_to_sources += plan.source_increment;
  if (composition && !touched.empty()) (void)synchronizeMetalMirrorFromComposition(state, touched);
  if (!touched.empty()) state.notifyGasThermodynamicMutation();
  std::sort(touched.begin(), touched.end());
  m_report.touched_rows = touched.size();
  return touched;
}

}  // namespace cosmosim::physics
