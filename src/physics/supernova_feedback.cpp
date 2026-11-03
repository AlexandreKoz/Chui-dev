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
  m_row_slot.assign(state.cells.size(), 0U);
  m_rows.clear();
  m_ledger = FeedbackLedger{};
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
  const auto consider = [&](std::size_t idx) {
    const std::uint32_t row = m_sorted_rows[idx];
    const double dx = image(state.cells.center_x_comoving[row] - p[0], 0);
    const double dy = image(state.cells.center_y_comoving[row] - p[1], 1);
    const double dz = image(state.cells.center_z_comoving[row] - p[2], 2);
    const double r2 = dx * dx + dy * dy + dz * dz;
    if (best.size() < k) {
      best.emplace_back(r2, row);
      std::push_heap(best.begin(), best.end());
    } else if (r2 < best.front().first || (r2 == best.front().first && row < best.front().second)) {
      std::pop_heap(best.begin(), best.end());
      best.back() = {r2, row};
      std::push_heap(best.begin(), best.end());
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
  std::sort(best.begin(), best.end());
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
      !(source.supernova_count >= 0.0) || !(source.presn_momentum_code >= 0.0)) {
    throw std::invalid_argument("feedback source budget must be non-negative");
  }
  if (source.ejecta_mass_code == 0.0 && source.supernova_energy_erg == 0.0 && source.presn_momentum_code == 0.0) return;
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
    const double rho = std::max(state.gas_cells.density_code[row], 1.0e-300);
    double mass = state.cells.mass_code[row];
    std::array<double, 3> momentum{mass * state.gas_cells.velocity_x_peculiar[row],
                                   mass * state.gas_cells.velocity_y_peculiar[row],
                                   mass * state.gas_cells.velocity_z_peculiar[row]};
    if (m_row_slot[row] != 0U) {
      const RowDelta& d = m_rows[m_row_slot[row] - 1U];
      mass += d.mass;
      for (int c = 0; c < 3; ++c) momentum[c] += d.momentum[c];
    }
    const double volume = state.cells.mass_code[row] / rho;
    const double h = std::cbrt(volume);
    std::array<double, 3> r{};
    const std::array<double, 3> center{state.cells.center_x_comoving[row], state.cells.center_y_comoving[row],
                                       state.cells.center_z_comoving[row]};
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
    const double x_h = state.gas_composition.enabled()
        ? state.gas_composition.value(state.gas_composition.descriptor().element_lane_begin, row) /
              std::max(compositionGasMass(state.gas_composition, row), 1.0e-300)
        : 0.76;
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
    p += t.w * source.presn_momentum_code;
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
  const double e_available = e_code + 0.5 * m_ej * v2_star;
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
    const std::uint32_t row = d.row;
    const double m_old = state.cells.mass_code[row];
    const double m_new = m_old + d.mass;
    if (!(m_new > 0.0)) throw std::runtime_error("feedback commit produced a non-positive gas mass");
    const double vx = state.gas_cells.velocity_x_peculiar[row];
    const double vy = state.gas_cells.velocity_y_peculiar[row];
    const double vz = state.gas_cells.velocity_z_peculiar[row];
    state.gas_cells.velocity_x_peculiar[row] = (m_old * vx + d.momentum[0]) / m_new;
    state.gas_cells.velocity_y_peculiar[row] = (m_old * vy + d.momentum[1]) / m_new;
    state.gas_cells.velocity_z_peculiar[row] = (m_old * vz + d.momentum[2]) / m_new;
    const bool sgs_lane = composition && state.gas_composition.descriptor().sgs_energy_lane !=
        core::GasCompositionLayoutDescriptor::k_no_lane;
    // Without an SGS reservoir the SGS share is thermal (never lost).
    const double thermal = d.thermal_energy + (sgs_lane ? 0.0 : d.sgs_energy);
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

}  // namespace cosmosim::physics
