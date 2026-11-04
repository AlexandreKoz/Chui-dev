#include "cosmosim/physics/supernova_feedback.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/ism_dust_sgs.hpp"

namespace cosmosim::physics {

double kimOstriker2015TerminalMomentumMsunKms(double energy_per_sn_erg, double n_h_cm3) {
  const double e51 = std::max(energy_per_sn_erg, 0.0) / 1.0e51;
  const double n = std::max(n_h_cm3, 1.0e-6);  // fit validity floor (KO15 tested n >= 0.001)
  return 2.17e5 * std::pow(e51, 0.93) * std::pow(n, -0.13);
}

MechanicalFeedbackDeposition::MechanicalFeedbackDeposition(MechanicalFeedbackConfig config, FeedbackUnits units)
    : m_config(config), m_units(units) {
  if (m_config.neighbor_count == 0U || !(m_units.mass_code_to_g > 0.0) || !(m_units.velocity_code_to_cm_s > 0.0) ||
      !(m_units.length_code_to_cm > 0.0) || !(m_units.density_code_to_cgs > 0.0) ||
      m_config.sgs_energy_fraction < 0.0 || m_config.sgs_energy_fraction > 1.0) {
    throw std::invalid_argument("MechanicalFeedbackDeposition: invalid configuration");
  }
}

void MechanicalFeedbackDeposition::setCandidates(const core::SimulationState& state,
                                                 std::vector<std::uint32_t> owned_leaf_rows) {
  std::sort(owned_leaf_rows.begin(), owned_leaf_rows.end(), [&](std::uint32_t a, std::uint32_t b) {
    const double xa = state.cells.center_x_comoving[a], xb = state.cells.center_x_comoving[b];
    return xa < xb || (xa == xb && a < b);
  });
  m_sorted_rows = std::move(owned_leaf_rows);
  m_sorted_x.resize(m_sorted_rows.size());
  for (std::size_t i = 0; i < m_sorted_rows.size(); ++i) m_sorted_x[i] = state.cells.center_x_comoving[m_sorted_rows[i]];
  m_local_count = static_cast<std::uint32_t>(state.cells.size());
  m_remote.clear();
  m_row_slot.assign(state.cells.size(), 0U);
  m_rows.clear();
  m_ledger = FeedbackLedger{};
}

void MechanicalFeedbackDeposition::setRemoteCandidates(const core::SimulationState& state,
                                                       std::vector<FeedbackRemoteCandidate> remote) {
  if (!m_rows.empty()) {
    throw std::logic_error("feedback remote candidates must be set before any source is staged");
  }
  std::sort(remote.begin(), remote.end(),
            [](const FeedbackRemoteCandidate& l, const FeedbackRemoteCandidate& r) { return l.gas_cell_id < r.gas_cell_id; });
  remote.erase(std::unique(remote.begin(), remote.end(),
                           [](const FeedbackRemoteCandidate& l, const FeedbackRemoteCandidate& r) {
                             return l.gas_cell_id == r.gas_cell_id;
                           }),
               remote.end());
  m_remote = std::move(remote);
  for (std::uint32_t i = 0; i < m_remote.size(); ++i) {
    if (m_remote[i].gas_cell_id == 0U || state.rowForGasCellId(m_remote[i].gas_cell_id).has_value()) {
      throw std::invalid_argument("feedback remote candidate is zero or already local");
    }
    m_sorted_rows.push_back(m_local_count + i);
  }
  std::sort(m_sorted_rows.begin(), m_sorted_rows.end(), [&](std::uint32_t a, std::uint32_t b) {
    const CandidateProps pa = props(state, a);
    const CandidateProps pb = props(state, b);
    return pa.center[0] < pb.center[0] || (pa.center[0] == pb.center[0] && pa.stable_id < pb.stable_id);
  });
  m_sorted_x.resize(m_sorted_rows.size());
  for (std::size_t i = 0; i < m_sorted_rows.size(); ++i) m_sorted_x[i] = props(state, m_sorted_rows[i]).center[0];
  m_row_slot.assign(static_cast<std::size_t>(m_local_count) + m_remote.size(), 0U);
}

MechanicalFeedbackDeposition::CandidateProps MechanicalFeedbackDeposition::props(
    const core::SimulationState& state, std::uint32_t candidate) const {
  CandidateProps result;
  if (candidate >= m_local_count) {
    const FeedbackRemoteCandidate& r = m_remote.at(candidate - m_local_count);
    result.center = r.center_comoving;
    result.mass = r.mass_code;
    result.momentum = r.momentum_code;
    result.density = r.density_code;
    result.x_h = r.hydrogen_mass_fraction;
    result.stable_id = r.gas_cell_id;
    return result;
  }
  const std::uint32_t row = candidate;
  result.center = {state.cells.center_x_comoving[row], state.cells.center_y_comoving[row],
                   state.cells.center_z_comoving[row]};
  result.mass = state.cells.mass_code[row];
  result.momentum = {result.mass * state.gas_cells.velocity_x_peculiar[row],
                     result.mass * state.gas_cells.velocity_y_peculiar[row],
                     result.mass * state.gas_cells.velocity_z_peculiar[row]};
  result.density = state.gas_cells.density_code[row];
  result.x_h = state.gas_composition.enabled()
      ? state.gas_composition.value(state.gas_composition.descriptor().element_lane_begin, row) /
            std::max(compositionGasMass(state.gas_composition, row), 1.0e-300)
      : 0.76;
  result.stable_id = state.gas_cells.gas_cell_id[row];
  return result;
}

double MechanicalFeedbackDeposition::localStencilRadius2(const core::SimulationState& state,
                                                         const std::array<double, 3>& position) const {
  const auto best = nearest(state, position);
  if (best.size() < m_config.neighbor_count) return std::numeric_limits<double>::infinity();
  return best.back().first;
}

std::vector<std::uint32_t> MechanicalFeedbackDeposition::localCandidatesWithin(
    const core::SimulationState& state, const std::array<double, 3>& position, double r2_max) const {
  std::vector<std::uint32_t> rows;
  for (const auto& [r2, candidate] : nearest(state, position)) {
    if (candidate < m_local_count && r2 <= r2_max) rows.push_back(candidate);
  }
  return rows;
}

std::vector<FeedbackCellDelta> MechanicalFeedbackDeposition::remoteDeltas(int source_rank) const {
  std::vector<FeedbackCellDelta> out;
  for (const RowDelta& d : m_rows) {
    if (d.row < m_local_count) continue;
    const FeedbackRemoteCandidate& r = m_remote.at(d.row - m_local_count);
    FeedbackCellDelta w;
    w.gas_cell_id = r.gas_cell_id;
    w.source_rank = source_rank;
    w.owner_rank = r.owner_rank;
    w.mass = d.mass;
    w.momentum = d.momentum;
    w.thermal_energy = d.thermal_energy;
    w.sgs_energy = d.sgs_energy;
    w.intended_kinetic = d.intended_kinetic;
    w.elements = d.elements;
    w.dust = d.dust;
    w.swept_mass_for_dust = d.swept_mass_for_dust;
    out.push_back(w);
  }
  std::sort(out.begin(), out.end(),
            [](const FeedbackCellDelta& l, const FeedbackCellDelta& r) { return l.gas_cell_id < r.gas_cell_id; });
  return out;
}

void MechanicalFeedbackDeposition::addImportedDelta(const core::SimulationState& state, const FeedbackCellDelta& in) {
  const auto row = state.rowForGasCellId(in.gas_cell_id);
  if (!row.has_value() || *row >= m_local_count) {
    throw std::runtime_error("feedback: imported delta targets a cell this rank does not own (fail closed)");
  }
  RowDelta& d = delta(*row);
  d.mass += in.mass;
  for (int c = 0; c < 3; ++c) d.momentum[c] += in.momentum[c];
  d.thermal_energy += in.thermal_energy;
  d.sgs_energy += in.sgs_energy;
  d.intended_kinetic += in.intended_kinetic;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) d.elements[e] += in.elements[e];
  for (std::size_t r = 0; r < k_dust_reservoir_count; ++r) d.dust[r] += in.dust[r];
  d.swept_mass_for_dust += in.swept_mass_for_dust;
}

std::vector<std::pair<double, std::uint32_t>> MechanicalFeedbackDeposition::nearest(
    const core::SimulationState& state, const std::array<double, 3>& p) const {
  std::vector<std::pair<double, std::uint32_t>> best;  // (r^2, row), max-heap by r^2
  const std::size_t n = m_sorted_rows.size();
  if (n == 0U) return best;
  const std::size_t k = std::min<std::size_t>(m_config.neighbor_count, n);
  const auto image = [&](double d, int axis) {
    if (m_units.periodic && m_units.box_extent_comoving[axis] > 0.0) {
      d -= m_units.box_extent_comoving[axis] * std::round(d / m_units.box_extent_comoving[axis]);
    }
    return d;
  };
  // Ties are broken by the stable gas-cell id (not the rank-local row), so
  // the stencil is identical under any MPI decomposition.
  const auto key_less = [&](const std::pair<double, std::uint32_t>& l, const std::pair<double, std::uint32_t>& r) {
    return l.first < r.first ||
        (l.first == r.first && props(state, l.second).stable_id < props(state, r.second).stable_id);
  };
  const auto consider = [&](std::size_t idx) {
    const std::uint32_t row = m_sorted_rows[idx];
    const std::array<double, 3> center = props(state, row).center;
    const double dx = image(center[0] - p[0], 0);
    const double dy = image(center[1] - p[1], 1);
    const double dz = image(center[2] - p[2], 2);
    const double r2 = dx * dx + dy * dy + dz * dz;
    const std::pair<double, std::uint32_t> item{r2, row};
    if (best.size() < k) {
      best.push_back(item);
      std::push_heap(best.begin(), best.end(), key_less);
    } else if (key_less(item, best.front())) {
      std::pop_heap(best.begin(), best.end(), key_less);
      best.back() = item;
      std::push_heap(best.begin(), best.end(), key_less);
    }
  };
  // Expand left/right along the x-sorted list (circularly when periodic)
  // until the x separation alone exceeds the current k-th distance.
  const std::size_t start = static_cast<std::size_t>(
      std::lower_bound(m_sorted_x.begin(), m_sorted_x.end(), p[0]) - m_sorted_x.begin());
  std::size_t visited = 0;
  std::ptrdiff_t left = static_cast<std::ptrdiff_t>(start) - 1;
  std::size_t right = start;
  bool left_open = true, right_open = true;
  while ((left_open || right_open) && visited < n) {
    const double bound = best.size() < k ? std::numeric_limits<double>::infinity() : best.front().first;
    if (right_open) {
      std::size_t idx = right;
      if (idx >= n) {
        if (m_units.periodic) idx -= n; else { right_open = false; }
      }
      if (right_open) {
        const double dx = image(m_sorted_x[idx] - p[0], 0);
        if (best.size() >= k && dx * dx > bound) right_open = false;
        else { consider(idx); ++visited; ++right; }
      }
    }
    if (left_open && visited < n) {
      std::ptrdiff_t idx = left;
      if (idx < 0) {
        if (m_units.periodic) idx += static_cast<std::ptrdiff_t>(n); else { left_open = false; }
      }
      if (left_open) {
        const double dx = image(m_sorted_x[static_cast<std::size_t>(idx)] - p[0], 0);
        if (best.size() >= k && dx * dx > bound) left_open = false;
        else { consider(static_cast<std::size_t>(idx)); ++visited; --left; }
      }
    }
  }
  std::sort(best.begin(), best.end(), key_less);
  return best;
}

MechanicalFeedbackDeposition::RowDelta& MechanicalFeedbackDeposition::delta(std::uint32_t row) {
  if (m_row_slot[row] == 0U) {
    m_rows.push_back(RowDelta{.row = row});
    m_row_slot[row] = static_cast<std::uint32_t>(m_rows.size());
  }
  return m_rows[m_row_slot[row] - 1U];
}

void MechanicalFeedbackDeposition::stage(const core::SimulationState& state, const FeedbackSourceBudget& source) {
  if (!(source.ejecta_mass_code >= 0.0) || !(source.supernova_energy_erg >= 0.0) ||
      !(source.supernova_count >= 0.0) || !(source.presn_momentum_code >= 0.0) ||
      !(source.wind_energy_erg >= 0.0) || !(source.wind_momentum_code >= 0.0) ||
      !(source.association_dispersion_energy_code >= 0.0)) {
    throw std::invalid_argument("feedback source budget must be non-negative");
  }
  if (source.ejecta_mass_code == 0.0 && source.supernova_energy_erg == 0.0 && source.presn_momentum_code == 0.0 &&
      source.wind_energy_erg == 0.0 && source.wind_momentum_code == 0.0 &&
      source.association_dispersion_energy_code == 0.0) return;
  const auto neighbors = nearest(state, source.position_comoving);
  if (neighbors.empty()) {
    throw std::runtime_error("feedback source " + std::to_string(source.source_id) +
                             " has no owned leaf gas cell to receive its budget (fail closed)");
  }
  const double a = m_units.scale_factor;
  const double erg_per_code = m_units.mass_code_to_g * m_units.velocity_code_to_cm_s * m_units.velocity_code_to_cm_s;
  const double e_code = source.supernova_energy_erg / erg_per_code;
  const double msun_kms_code = (k_solar_mass_g / m_units.mass_code_to_g) * (1.0e5 / m_units.velocity_code_to_cm_s);
  const std::size_t k = neighbors.size();

  struct Target {
    std::uint32_t row;
    std::array<double, 3> dir;
    double w;
    double mass;
    std::array<double, 3> momentum;  // current (state + staged) momentum
    double n_h;
  };
  std::vector<Target> targets(k);
  double w_sum = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    const std::uint32_t row = neighbors[i].second;
    const CandidateProps cell = props(state, row);
    const double rho = std::max(cell.density, 1.0e-300);
    double mass = cell.mass;
    std::array<double, 3> momentum = cell.momentum;
    if (!m_config.step_start_reference && m_row_slot[row] != 0U) {
      const RowDelta& d = m_rows[m_row_slot[row] - 1U];
      mass += d.mass;
      for (int c = 0; c < 3; ++c) momentum[c] += d.momentum[c];
    }
    const double volume = cell.mass / rho;
    const double h = std::cbrt(volume);
    std::array<double, 3> r{};
    const std::array<double, 3> center = cell.center;
    double r2 = 0.0;
    for (int c = 0; c < 3; ++c) {
      double d = center[c] - source.position_comoving[c];
      if (m_units.periodic && m_units.box_extent_comoving[c] > 0.0) {
        d -= m_units.box_extent_comoving[c] * std::round(d / m_units.box_extent_comoving[c]);
      }
      r[c] = d;
      r2 += d * d;
    }
    const double rr = std::sqrt(r2);
    std::array<double, 3> dir{};
    if (rr > 1.0e-12 * h) for (int c = 0; c < 3; ++c) dir[c] = r[c] / rr;
    const double w = std::pow(volume, 2.0 / 3.0) / std::max(r2, 0.25 * h * h);
    const double x_h = cell.x_h;
    const double rho_cgs = rho * (m_units.comoving ? 1.0 / (a * a * a) : 1.0) * m_units.density_code_to_cgs;
    const double n_h = x_h * rho_cgs / (chemicalElementInfo(ChemicalElement::kH).atomic_mass_amu * k_amu_g);
    targets[i] = Target{row, dir, w, mass, momentum, n_h};
    w_sum += w;
  }
  for (Target& t : targets) t.w /= w_sum;

  // Momentum magnitudes (star frame).
  const double m_ej = source.ejecta_mass_code;
  const double n_sn = source.supernova_count;
  const double e_per_sn = n_sn > 0.0 ? source.supernova_energy_erg / n_sn : 0.0;
  std::vector<std::array<double, 3>> dp_star(k);
  double radial_requested = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    const Target& t = targets[i];
    double p = 0.0;
    if (e_code > 0.0) {
      const double p_energy = std::sqrt(2.0 * e_code * t.w * (t.w * m_ej + t.mass));
      const double p_term = n_sn > 0.0
          ? t.w * n_sn * kimOstriker2015TerminalMomentumMsunKms(e_per_sn, t.n_h) * msun_kms_code
          : p_energy;
      if (p_energy <= p_term) { p = p_energy; ++m_ledger.resolved_cells; }
      else { p = p_term; ++m_ledger.unresolved_cells; }
    }
    p += t.w * (source.presn_momentum_code + source.wind_momentum_code);
    radial_requested += p;
    for (int c = 0; c < 3; ++c) dp_star[i][c] = p * t.dir[c];
  }
  // Remove the net vector of the asymmetric stencil (mass weighted).
  std::array<double, 3> net{};
  double mass_after_sum = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    for (int c = 0; c < 3; ++c) net[c] += dp_star[i][c];
    mass_after_sum += targets[i].mass + targets[i].w * m_ej;
  }
  for (std::size_t i = 0; i < k; ++i) {
    const double share = (targets[i].mass + targets[i].w * m_ej) / mass_after_sum;
    for (int c = 0; c < 3; ++c) dp_star[i][c] -= net[c] * share;
  }
  // Energy closure with a common momentum scale s in [0, 1].
  const double v2_star = source.velocity_peculiar[0] * source.velocity_peculiar[0] +
      source.velocity_peculiar[1] * source.velocity_peculiar[1] +
      source.velocity_peculiar[2] * source.velocity_peculiar[2];
  // Lab-frame ejecta kinetic energy = 1/2 M V^2 + internal association
  // dispersion (zero for single stars), so merged sources conserve energy.
  const double e_available = e_code + source.wind_energy_erg / erg_per_code + 0.5 * m_ej * v2_star +
      source.association_dispersion_energy_code;
  const auto kinetic_change = [&](double s) {
    double dke = 0.0;
    for (std::size_t i = 0; i < k; ++i) {
      const Target& t = targets[i];
      const double dm = t.w * m_ej;
      double p_new2 = 0.0, p_old2 = 0.0;
      for (int c = 0; c < 3; ++c) {
        const double pn = t.momentum[c] + s * dp_star[i][c] + dm * source.velocity_peculiar[c];
        p_new2 += pn * pn;
        p_old2 += t.momentum[c] * t.momentum[c];
      }
      dke += p_new2 / (2.0 * (t.mass + dm)) - (t.mass > 0.0 ? p_old2 / (2.0 * t.mass) : 0.0);
    }
    return dke;
  };
  double s = 1.0;
  double dke = kinetic_change(1.0);
  if (dke > e_available) {
    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < 80; ++it) {
      const double mid = 0.5 * (lo + hi);
      if (kinetic_change(mid) > e_available) hi = mid; else lo = mid;
    }
    s = lo;
    dke = kinetic_change(s);
    m_ledger.momentum_scaled_down_sources += 1.0;
  }
  const double thermal_total = std::max(e_available - dke, 0.0);
  // Stage deltas.
  std::array<double, 3> injected{};
  double radial_delivered = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    const Target& t = targets[i];
    RowDelta& d = delta(t.row);
    const double dm = t.w * m_ej;
    {
      // Intended kinetic change of this target for this source (used by the
      // per-cell closure at commit in step_start_reference mode).
      double p_new2 = 0.0, p_old2 = 0.0;
      for (int c = 0; c < 3; ++c) {
        const double pn = t.momentum[c] + s * dp_star[i][c] + dm * source.velocity_peculiar[c];
        p_new2 += pn * pn;
        p_old2 += t.momentum[c] * t.momentum[c];
      }
      d.intended_kinetic += p_new2 / (2.0 * (t.mass + dm)) - (t.mass > 0.0 ? p_old2 / (2.0 * t.mass) : 0.0);
    }
    if (t.row >= m_local_count) ++m_ledger.remote_target_cells;
    d.mass += dm;
    for (int c = 0; c < 3; ++c) {
      const double dp = s * dp_star[i][c] + dm * source.velocity_peculiar[c];
      d.momentum[c] += dp;
      injected[c] += dp;
      radial_delivered += s * dp_star[i][c] * t.dir[c];
    }
    const double thermal = t.w * thermal_total;
    d.sgs_energy += m_config.sgs_energy_fraction * thermal;
    d.thermal_energy += (1.0 - m_config.sgs_energy_fraction) * thermal;
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) d.elements[e] += t.w * source.ejecta_elements_code[e];
    for (std::size_t r = 0; r < k_dust_reservoir_count; ++r) d.dust[r] += t.w * source.ejecta_dust_code[r];
    if (n_sn > 0.0) {
      d.swept_mass_for_dust += m_config.dust_sn_destruction_efficiency * t.w * n_sn *
          (m_config.dust_sn_swept_mass_per_1e51_erg_g / m_units.mass_code_to_g) * (e_per_sn / 1.0e51);
    }
  }
  double residual2 = 0.0;
  for (int c = 0; c < 3; ++c) {
    const double r = injected[c] - m_ej * source.velocity_peculiar[c];
    residual2 += r * r;
  }
  m_ledger.net_vector_momentum_residual_code = std::max(m_ledger.net_vector_momentum_residual_code, std::sqrt(residual2));
  m_ledger.radial_momentum_requested_code += radial_requested;
  m_ledger.radial_momentum_delivered_code += radial_delivered;
  m_ledger.mass_injected_code += m_ej;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) m_ledger.elements_injected_code[e] += source.ejecta_elements_code[e];
  for (std::size_t r = 0; r < k_dust_reservoir_count; ++r) m_ledger.dust_injected_code += source.ejecta_dust_code[r];
  m_ledger.energy_budget_code += e_available;
  m_ledger.kinetic_energy_injected_code += dke;
  m_ledger.thermal_energy_injected_code += (1.0 - m_config.sgs_energy_fraction) * thermal_total;
  m_ledger.sgs_energy_injected_code += m_config.sgs_energy_fraction * thermal_total;
  ++m_ledger.sources;
}

std::vector<std::uint32_t> MechanicalFeedbackDeposition::commit(core::SimulationState& state) {
  std::sort(m_rows.begin(), m_rows.end(), [](const RowDelta& a, const RowDelta& b) { return a.row < b.row; });
  std::vector<std::uint32_t> touched;
  touched.reserve(m_rows.size());
  const bool composition = state.gas_composition.enabled();
  std::optional<CompositionLayout> layout;
  if (composition) layout = CompositionLayout::fromDescriptor(state.gas_composition.descriptor());
  for (const RowDelta& d : m_rows) {
    if (d.row >= m_local_count) {
      continue;  // remote target: delivered to its owner via remoteDeltas()
    }
    const std::uint32_t row = d.row;
    const double m_old = state.cells.mass_code[row];
    const double m_new = m_old + d.mass;
    if (!(m_new > 0.0)) throw std::runtime_error("feedback commit produced a non-positive gas mass");
    const double vx = state.gas_cells.velocity_x_peculiar[row];
    const double vy = state.gas_cells.velocity_y_peculiar[row];
    const double vz = state.gas_cells.velocity_z_peculiar[row];
    double closure_thermal = 0.0;
    if (m_config.step_start_reference) {
      // Per-cell energy closure: the sources were evaluated independently
      // against the step-start state, so the actual kinetic change of the
      // summed momentum differs from the sum of intended changes by the
      // cross terms; that difference is moved into thermal energy (exact
      // total energy), and any part the thermal floor cannot absorb is
      // ledgered as unclosed.
      const double px = m_old * vx + d.momentum[0];
      const double py = m_old * vy + d.momentum[1];
      const double pz = m_old * vz + d.momentum[2];
      const double ke_new = (px * px + py * py + pz * pz) / (2.0 * m_new);
      const double ke_old = 0.5 * m_old * (vx * vx + vy * vy + vz * vz);
      closure_thermal = d.intended_kinetic - (ke_new - ke_old);
      const double thermal_available = m_old * state.gas_cells.internal_energy_code[row] + d.thermal_energy;
      if (closure_thermal < -thermal_available) {
        m_ledger.kinetic_cross_term_unclosed_code += -thermal_available - closure_thermal;
        closure_thermal = -thermal_available;
      }
      m_ledger.kinetic_cross_term_to_thermal_code += closure_thermal;
    }
    state.gas_cells.velocity_x_peculiar[row] = (m_old * vx + d.momentum[0]) / m_new;
    state.gas_cells.velocity_y_peculiar[row] = (m_old * vy + d.momentum[1]) / m_new;
    state.gas_cells.velocity_z_peculiar[row] = (m_old * vz + d.momentum[2]) / m_new;
    const bool sgs_lane = composition && state.gas_composition.descriptor().sgs_energy_lane !=
        core::GasCompositionLayoutDescriptor::k_no_lane;
    // Without an SGS reservoir the SGS share is thermal (never lost).
    const double thermal = d.thermal_energy + (sgs_lane ? 0.0 : d.sgs_energy) + closure_thermal;
    state.gas_cells.internal_energy_code[row] =
        (m_old * state.gas_cells.internal_energy_code[row] + thermal) / m_new;
    state.gas_cells.density_code[row] *= m_new / std::max(m_old, 1.0e-300);
    state.cells.mass_code[row] = m_new;
    if (composition) {
      auto& sidecar = state.gas_composition;
      for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
        sidecar.addValue(layout->elementLane(static_cast<ChemicalElement>(e)), row, d.elements[e]);
      }
      // Ejecta enter as neutral atoms (H -> HI, He -> HeI, C -> C+ per the
      // detailed-network convention); the next thermochemistry step evolves
      // them from the post-feedback temperature (no stale electron fraction is
      // reused: species are re-solved before any rate uses them).
      const std::size_t species = sidecar.descriptor().species_lane_count;
      if (species >= k_hybrid_species_count) {
        sidecar.addValue(layout->speciesLane(static_cast<std::size_t>(ChemicalSpecies::kHI)), row, d.elements[0]);
        sidecar.addValue(layout->speciesLane(static_cast<std::size_t>(ChemicalSpecies::kHeI)), row, d.elements[1]);
      }
      if (species >= k_detailed_species_count) {
        sidecar.addValue(layout->speciesLane(static_cast<std::size_t>(ChemicalSpecies::kCII)), row,
                         d.elements[static_cast<std::size_t>(ChemicalElement::kC)]);
      }
      for (std::size_t r = 0; r < sidecar.descriptor().dust_lane_count; ++r) {
        sidecar.addValue(layout->dustLane(r), row, d.dust[r]);
      }
      if (sgs_lane) sidecar.addValue(sidecar.descriptor().sgs_energy_lane, row, d.sgs_energy);
      if (d.swept_mass_for_dust > 0.0) {
        m_ledger.dust_destroyed_code +=
            destroyDustBySupernova(sidecar, *layout, row, d.swept_mass_for_dust, compositionGasMass(sidecar, row));
      }
      sidecar.mean_molecular_weight[row] = std::numeric_limits<double>::quiet_NaN();  // re-derived by chemistry
    } else {
      double metals = 0.0;
      for (std::size_t e = 2; e < k_chemical_element_count; ++e) metals += d.elements[e];
      state.gas_cells.metal_mass_code[row] += metals;
    }
    touched.push_back(row);
  }
  if (composition) (void)synchronizeMetalMirrorFromComposition(state, touched);
  for (const RowDelta& d : m_rows) m_row_slot[d.row] = 0U;
  m_rows.clear();
  return touched;
}

namespace {

[[nodiscard]] double minimumImage(double d, bool periodic, double extent) noexcept {
  if (periodic && extent > 0.0) d -= extent * std::round(d / extent);
  return d;
}

}  // namespace

std::vector<std::uint64_t> discoverAssociationLeaders(
    std::span<const AssociationSourceDescriptor> sources, double radius_code, double window_years,
    bool periodic, const std::array<double, 3>& box_extent_comoving) {
  if (!(radius_code > 0.0) || !std::isfinite(radius_code) || !(window_years > 0.0) ||
      !std::isfinite(window_years)) {
    throw std::invalid_argument("discoverAssociationLeaders: radius and window must be finite and positive");
  }
  const std::size_t n = sources.size();
  {
    std::vector<std::uint64_t> ids(n);
    for (std::size_t i = 0; i < n; ++i) ids[i] = sources[i].source_id;
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) {
      throw std::invalid_argument("discoverAssociationLeaders: duplicate source id");
    }
  }
  // Uniform hash grid with cell width >= radius: linked pairs are always in
  // the same or an adjacent cell. Periodic axes use an integer number of
  // cells over the box so the wrap is exact.
  std::array<std::int64_t, 3> cells_per_axis{0, 0, 0};
  std::array<double, 3> cell_width{radius_code, radius_code, radius_code};
  for (int c = 0; c < 3; ++c) {
    if (periodic && box_extent_comoving[c] > 0.0) {
      cells_per_axis[c] = std::max<std::int64_t>(
          1, static_cast<std::int64_t>(std::floor(box_extent_comoving[c] / radius_code)));
      cell_width[c] = box_extent_comoving[c] / static_cast<double>(cells_per_axis[c]);
    }
  }
  struct Entry {
    std::array<std::int64_t, 3> key{};
    std::uint64_t id = 0U;
    std::size_t index = 0U;
  };
  std::vector<Entry> entries(n);
  for (std::size_t i = 0; i < n; ++i) {
    Entry e;
    e.id = sources[i].source_id;
    e.index = i;
    for (int c = 0; c < 3; ++c) {
      double x = sources[i].position_comoving[c];
      if (!std::isfinite(x) || !std::isfinite(sources[i].birth_years)) {
        throw std::invalid_argument("discoverAssociationLeaders: non-finite source position or birth time");
      }
      if (cells_per_axis[c] > 0) {
        x -= box_extent_comoving[c] * std::floor(x / box_extent_comoving[c]);
        e.key[c] = std::clamp<std::int64_t>(static_cast<std::int64_t>(std::floor(x / cell_width[c])), 0,
                                             cells_per_axis[c] - 1);
      } else {
        e.key[c] = static_cast<std::int64_t>(std::floor(x / cell_width[c]));
      }
    }
    entries[i] = e;
  }
  std::sort(entries.begin(), entries.end(), [](const Entry& l, const Entry& r) {
    return l.key != r.key ? l.key < r.key : l.id < r.id;
  });
  // Union-find over sorted positions; the root of every set is its minimum id
  // (union attaches the larger-id root below the smaller), so the result is
  // independent of the visit order.
  std::vector<std::size_t> parent(n);
  for (std::size_t i = 0; i < n; ++i) parent[i] = i;
  const auto find = [&parent](std::size_t x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  const auto unite = [&](std::size_t a, std::size_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return;
    if (entries[a].id < entries[b].id) parent[b] = a; else parent[a] = b;
  };
  const double r2 = radius_code * radius_code;
  std::vector<std::array<std::int64_t, 3>> neighbor_keys;
  neighbor_keys.reserve(27U);
  for (std::size_t i = 0; i < n; ++i) {
    neighbor_keys.clear();
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          std::array<std::int64_t, 3> k{entries[i].key[0] + dx, entries[i].key[1] + dy, entries[i].key[2] + dz};
          for (int c = 0; c < 3; ++c) {
            if (cells_per_axis[c] > 0) k[c] = ((k[c] % cells_per_axis[c]) + cells_per_axis[c]) % cells_per_axis[c];
          }
          neighbor_keys.push_back(k);
        }
      }
    }
    std::sort(neighbor_keys.begin(), neighbor_keys.end());
    neighbor_keys.erase(std::unique(neighbor_keys.begin(), neighbor_keys.end()), neighbor_keys.end());
    const AssociationSourceDescriptor& a = sources[entries[i].index];
    for (const auto& k : neighbor_keys) {
      auto it = std::lower_bound(entries.begin(), entries.end(), k,
                                 [](const Entry& e, const std::array<std::int64_t, 3>& key) { return e.key < key; });
      for (; it != entries.end() && it->key == k; ++it) {
        const std::size_t j = static_cast<std::size_t>(it - entries.begin());
        if (j <= i) continue;  // each unordered pair once
        const AssociationSourceDescriptor& b = sources[it->index];
        if (std::abs(a.birth_years - b.birth_years) > window_years) continue;
        double d2 = 0.0;
        for (int c = 0; c < 3; ++c) {
          const double d = minimumImage(b.position_comoving[c] - a.position_comoving[c], periodic,
                                        box_extent_comoving[c]);
          d2 += d * d;
        }
        if (d2 <= r2) unite(i, j);
      }
    }
  }
  std::vector<std::uint64_t> leaders(n);
  for (std::size_t i = 0; i < n; ++i) leaders[entries[i].index] = entries[find(i)].id;
  return leaders;
}

FeedbackSourceBudget mergeAssociationMembers(
    std::span<const FeedbackSourceBudget> members, bool periodic, const std::array<double, 3>& box_extent_comoving) {
  if (members.empty()) throw std::invalid_argument("mergeAssociationMembers: empty association");
  std::vector<std::size_t> order(members.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](std::size_t l, std::size_t r) { return members[l].source_id < members[r].source_id; });
  for (std::size_t i = 1; i < order.size(); ++i) {
    if (members[order[i]].source_id == members[order[i - 1U]].source_id) {
      throw std::invalid_argument("mergeAssociationMembers: duplicate member source id");
    }
  }
  const FeedbackSourceBudget& leader = members[order.front()];
  if (order.size() == 1U) return leader;
  FeedbackSourceBudget group = leader;
  group.member_count = 0U;
  group.ejecta_mass_code = 0.0;
  group.ejecta_elements_code = {};
  group.ejecta_dust_code = {};
  group.supernova_count = 0.0;
  group.supernova_energy_erg = 0.0;
  group.presn_momentum_code = 0.0;
  group.wind_energy_erg = 0.0;
  group.wind_momentum_code = 0.0;
  group.association_dispersion_energy_code = 0.0;
  std::array<double, 3> weighted_offset{};
  std::array<double, 3> plain_offset{};
  std::array<double, 3> momentum{};
  std::array<double, 3> velocity_sum{};
  for (const std::size_t idx : order) {
    const FeedbackSourceBudget& m = members[idx];
    group.member_count += m.member_count;
    group.ejecta_mass_code += m.ejecta_mass_code;
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) group.ejecta_elements_code[e] += m.ejecta_elements_code[e];
    for (std::size_t r = 0; r < k_dust_reservoir_count; ++r) group.ejecta_dust_code[r] += m.ejecta_dust_code[r];
    group.supernova_count += m.supernova_count;
    group.supernova_energy_erg += m.supernova_energy_erg;
    group.presn_momentum_code += m.presn_momentum_code;
    group.wind_energy_erg += m.wind_energy_erg;
    group.wind_momentum_code += m.wind_momentum_code;
    group.association_dispersion_energy_code += m.association_dispersion_energy_code;
    for (int c = 0; c < 3; ++c) {
      const double d = minimumImage(m.position_comoving[c] - leader.position_comoving[c], periodic,
                                    box_extent_comoving[c]);
      weighted_offset[c] += m.ejecta_mass_code * d;
      plain_offset[c] += d;
      momentum[c] += m.ejecta_mass_code * m.velocity_peculiar[c];
      velocity_sum[c] += m.velocity_peculiar[c];
    }
  }
  const double n = static_cast<double>(order.size());
  const double mass = group.ejecta_mass_code;
  for (int c = 0; c < 3; ++c) {
    double x = leader.position_comoving[c] + (mass > 0.0 ? weighted_offset[c] / mass : plain_offset[c] / n);
    if (periodic && box_extent_comoving[c] > 0.0) {
      x -= box_extent_comoving[c] * std::floor(x / box_extent_comoving[c]);
    }
    group.position_comoving[c] = x;
    group.velocity_peculiar[c] = mass > 0.0 ? momentum[c] / mass : velocity_sum[c] / n;
  }
  // Internal kinetic energy in the centre-of-momentum frame (non-negative by
  // construction; summed in id order).
  if (mass > 0.0) {
    for (const std::size_t idx : order) {
      const FeedbackSourceBudget& m = members[idx];
      double dv2 = 0.0;
      for (int c = 0; c < 3; ++c) {
        const double dv = m.velocity_peculiar[c] - group.velocity_peculiar[c];
        dv2 += dv * dv;
      }
      group.association_dispersion_energy_code += 0.5 * m.ejecta_mass_code * dv2;
    }
  }
  group.source_id = leader.source_id;
  return group;
}

std::vector<FeedbackSourceBudget> mergeAssociationSources(
    std::span<const FeedbackSourceBudget> sources, std::span<const double> birth_years,
    double radius_code, double window_years, bool periodic, const std::array<double, 3>& box_extent_comoving) {
  if (sources.size() != birth_years.size()) {
    throw std::invalid_argument("mergeAssociationSources: birth time list does not match sources");
  }
  std::vector<AssociationSourceDescriptor> descriptors(sources.size());
  for (std::size_t i = 0; i < sources.size(); ++i) {
    descriptors[i].source_id = sources[i].source_id;
    descriptors[i].position_comoving = sources[i].position_comoving;
    descriptors[i].birth_years = birth_years[i];
  }
  const std::vector<std::uint64_t> leaders =
      discoverAssociationLeaders(descriptors, radius_code, window_years, periodic, box_extent_comoving);
  std::vector<std::size_t> order(sources.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::size_t l, std::size_t r) {
    return leaders[l] != leaders[r] ? leaders[l] < leaders[r] : sources[l].source_id < sources[r].source_id;
  });
  std::vector<FeedbackSourceBudget> merged;
  std::vector<FeedbackSourceBudget> members;
  std::size_t begin = 0U;
  while (begin < order.size()) {
    std::size_t end = begin + 1U;
    while (end < order.size() && leaders[order[end]] == leaders[order[begin]]) ++end;
    members.clear();
    for (std::size_t k = begin; k < end; ++k) members.push_back(sources[order[k]]);
    merged.push_back(mergeAssociationMembers(members, periodic, box_extent_comoving));
    begin = end;
  }
  return merged;
}

}  // namespace cosmosim::physics
