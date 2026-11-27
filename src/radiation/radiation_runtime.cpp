#include "cosmosim/radiation/radiation_runtime.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/source_transaction.hpp"
#include "cosmosim/radiation/radiation_closure.hpp"
#include "cosmosim/radiation/radiation_cosmology.hpp"

namespace cosmosim::radiation {
namespace {

[[nodiscard]] std::vector<std::uint32_t> ownedLeafPatches(const core::SimulationState& state, std::uint32_t rank) {
  std::vector<std::uint64_t> parents;
  parents.reserve(state.patches.size());
  for (std::size_t p = 0; p < state.patches.size(); ++p) {
    if (state.patches.parent_patch_id[p] != 0U) parents.push_back(state.patches.parent_patch_id[p]);
  }
  std::sort(parents.begin(), parents.end());
  std::vector<std::uint32_t> leaves;
  for (std::size_t p = 0; p < state.patches.size(); ++p) {
    if (std::binary_search(parents.begin(), parents.end(), state.patches.patch_id[p])) continue;  // covered parent
    if (state.patches.owning_rank.size() == state.patches.size() && state.patches.owning_rank[p] != rank) continue;
    leaves.push_back(static_cast<std::uint32_t>(p));
  }
  return leaves;
}

[[nodiscard]] core::RadiationPatchRecord recordFor(const core::SimulationState& state, std::uint32_t p) {
  core::RadiationPatchRecord r;
  r.patch_id = state.patches.patch_id[p];
  r.geometry_generation = state.patches.morton_key[p];
  r.level = state.patches.level[p];
  r.nx = state.patches.cell_dim_x[p];
  r.ny = state.patches.cell_dim_y[p];
  r.nz = state.patches.cell_dim_z[p];
  r.origin_x_comoving = state.patches.origin_x_comoving[p];
  r.origin_y_comoving = state.patches.origin_y_comoving[p];
  r.origin_z_comoving = state.patches.origin_z_comoving[p];
  if (r.nx == 0U || r.ny == 0U || r.nz == 0U) {
    throw std::runtime_error("radiation: leaf patch without explicit geometry (production RT needs patch lanes)");
  }
  r.cell_width_x_comoving = state.patches.extent_x_comoving[p] / r.nx;
  r.cell_width_y_comoving = state.patches.extent_y_comoving[p] / r.ny;
  r.cell_width_z_comoving = state.patches.extent_z_comoving[p] / r.nz;
  if (static_cast<std::uint64_t>(r.nx) * r.ny * r.nz != state.patches.cell_count[p]) {
    throw std::runtime_error("radiation: patch cell_count disagrees with its dimensions");
  }
  return r;
}

struct GasWork {
  std::vector<std::uint32_t> row;   // per record cell (flattened like one group of the lanes)
  std::vector<double> density;      // proper code
  std::vector<double> tpu;          // K per code specific energy
  std::vector<double> etot;         // proper code energy density
  std::vector<double> px, py, pz;   // proper code momentum density
};

}  // namespace

std::uint64_t radiationCanonicalBytes(std::uint64_t cells, std::uint32_t groups) {
  return 4ULL * groups * 8ULL * cells;
}

void allocateRadiationState(core::SimulationState& state, const RadiationSpectralDefinition& sp, double c_hat_over_c,
                            double a, std::uint32_t rank) {
  core::RadiationMomentSoa rad;
  rad.spectral_profile = sp.profile;
  rad.spectral_definition_sha256 = sp.sha256;
  rad.group_count = static_cast<std::uint32_t>(sp.groupCount());
  rad.group_edges_ev = sp.edges_ev;
  rad.c_hat_over_c = c_hat_over_c;
  rad.scale_factor = a;
  std::uint64_t offset = 0U;
  for (const std::uint32_t p : ownedLeafPatches(state, rank)) {
    core::RadiationPatchRecord r = recordFor(state, p);
    r.value_offset = offset;
    offset += static_cast<std::uint64_t>(r.nx) * r.ny * r.nz * rad.group_count;
    rad.patches.push_back(r);
  }
  rad.energy_density_comoving.assign(offset, 0.0);
  rad.flux_x_comoving.assign(offset, 0.0);
  rad.flux_y_comoving.assign(offset, 0.0);
  rad.flux_z_comoving.assign(offset, 0.0);
  if (!rad.isConsistent()) throw std::logic_error("allocateRadiationState produced an inconsistent state");
  state.radiation = std::move(rad);
}

bool radiationRecordsMatchPatches(const core::SimulationState& state, std::uint32_t rank, std::string* reason) {
  const std::vector<std::uint32_t> leaves = ownedLeafPatches(state, rank);
  if (leaves.size() != state.radiation.patches.size()) {
    if (reason != nullptr) *reason = "leaf patch count changed (regrid/migration without radiation transfer)";
    return false;
  }
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    const core::RadiationPatchRecord expected = recordFor(state, leaves[i]);
    const core::RadiationPatchRecord& got = state.radiation.patches[i];
    if (expected.patch_id != got.patch_id || expected.nx != got.nx || expected.ny != got.ny || expected.nz != got.nz ||
        expected.level != got.level || expected.origin_x_comoving != got.origin_x_comoving ||
        expected.origin_y_comoving != got.origin_y_comoving || expected.origin_z_comoving != got.origin_z_comoving) {
      if (reason != nullptr) *reason = "radiation record does not match leaf patch " + std::to_string(expected.patch_id);
      return false;
    }
  }
  return true;
}

RadiationRuntime::RadiationRuntime(RadiationRuntimeSettings settings) : m_settings(std::move(settings)) {
  if (m_settings.spectrum.groupCount() == 0U) throw std::invalid_argument("RadiationRuntime: no radiation groups");
  if (!(m_settings.coupling.c_code > 0.0) || !(m_settings.coupling.c_hat_code > 0.0) ||
      m_settings.coupling.c_hat_code > m_settings.coupling.c_code * (1.0 + 1.0e-15)) {
    throw std::invalid_argument("RadiationRuntime: requires 0 < c_hat <= c");
  }
  if (m_settings.coupling_model == RadiationCouplingModel::kGreyThermalReference && !m_settings.spectrum.grey) {
    throw std::invalid_argument("RadiationRuntime: grey thermal coupling requires the grey profile");
  }
}

RadiationIntervalReport RadiationRuntime::advanceInterval(core::SimulationState& state,
                                                          const RadiationIntervalInputs& in) const {
  RadiationIntervalReport report;
  core::RadiationMomentSoa& rad = state.radiation;
  if (!rad.isConsistent() || rad.group_count != m_settings.spectrum.groupCount()) {
    throw std::runtime_error("RadiationRuntime: canonical radiation state is inconsistent with the spectrum");
  }
  if (rad.spectral_definition_sha256 != m_settings.spectrum.sha256) {
    throw std::runtime_error("RadiationRuntime: stored spectral definition hash differs from the configured one");
  }
  if (!(in.dt_proper_code > 0.0)) throw std::invalid_argument("RadiationRuntime: dt must be > 0");
  const double c = m_settings.coupling.c_code;
  const double ch = m_settings.coupling.c_hat_code;
  const std::uint32_t ng = rad.group_count;
  const std::size_t nvalues = rad.energy_density_comoving.size();
  const RadiationTopology topology(rad, m_settings.transport);
  const RadiationLevelSolver solver(rad, topology, m_settings.transport);
  // F01: reserve the interval work copies before allocating them.
  core::MemoryReservation workspace_reservation;
  if (m_settings.governor != nullptr) {
    std::size_t record_cells = 0;
    for (const core::RadiationPatchRecord& record : rad.patches) {
      record_cells += static_cast<std::size_t>(record.nx) * record.ny * record.nz;
    }
    workspace_reservation = m_settings.governor->reserve(
        core::MemoryClass::kPhaseResident,
        radiationIntervalWorkspaceBytes(nvalues, record_cells, m_settings.coupling_model != RadiationCouplingModel::kNone),
        "radiation_interval_workspace");
  }

  // Working state (E, G = F c_hat/c), comoving.
  RadiationWorkState u;
  u.resize(nvalues);
  for (std::size_t i = 0; i < nvalues; ++i) {
    u.e[i] = rad.energy_density_comoving[i];
    u.gx[i] = rad.flux_x_comoving[i] * ch / c;
    u.gy[i] = rad.flux_y_comoving[i] * ch / c;
    u.gz[i] = rad.flux_z_comoving[i] * ch / c;
  }
  // Source rates (comoving energy density per proper time).
  std::vector<double> source(nvalues, 0.0);
  for (const RadiationSourceDeposit& d : in.sources) {
    if (d.record >= rad.patches.size() || d.group >= ng) throw std::out_of_range("RadiationRuntime: bad source deposit");
    const core::RadiationPatchRecord& r = rad.patches[d.record];
    const std::size_t nc = static_cast<std::size_t>(r.nx) * r.ny * r.nz;
    if (d.cell >= nc || !(d.comoving_energy_rate >= 0.0)) throw std::out_of_range("RadiationRuntime: bad source cell");
    source[r.value_offset + static_cast<std::size_t>(d.group) * nc + d.cell] += d.comoving_energy_rate;
  }
  // Gas work arrays (one entry per record cell).
  std::unordered_map<std::uint64_t, std::uint32_t> patch_by_id;
  for (std::size_t p = 0; p < state.patches.size(); ++p) patch_by_id[state.patches.patch_id[p]] = static_cast<std::uint32_t>(p);
  std::vector<std::size_t> cell_base(rad.patches.size() + 1U, 0U);
  for (std::size_t r = 0; r < rad.patches.size(); ++r) {
    cell_base[r + 1U] = cell_base[r] + static_cast<std::size_t>(rad.patches[r].nx) * rad.patches[r].ny * rad.patches[r].nz;
  }
  const std::size_t ncells = cell_base.back();
  GasWork gas;
  const bool coupled = m_settings.coupling_model != RadiationCouplingModel::kNone;
  if (coupled) {
    gas.row.resize(ncells);
    gas.density.resize(ncells);
    gas.tpu.resize(ncells);
    gas.etot.resize(ncells);
    gas.px.resize(ncells);
    gas.py.resize(ncells);
    gas.pz.resize(ncells);
    for (std::size_t r = 0; r < rad.patches.size(); ++r) {
      const auto it = patch_by_id.find(rad.patches[r].patch_id);
      if (it == patch_by_id.end()) throw std::runtime_error("RadiationRuntime: record without hydro patch");
      const std::uint32_t first = state.patches.first_cell[it->second];
      for (std::size_t c0 = 0; c0 < cell_base[r + 1U] - cell_base[r]; ++c0) {
        const std::size_t k = cell_base[r] + c0;
        const std::uint32_t row = first + static_cast<std::uint32_t>(c0);
        gas.row[k] = row;
        const double rho = in.frame_end.properDensity(state.gas_cells.density_code[row]);
        const double uu = state.gas_cells.internal_energy_code[row];
        const double tt = state.gas_cells.temperature_code[row];
        if (!(rho > 0.0) || !(uu > 0.0) || !(tt > 0.0)) {
          throw std::runtime_error("RadiationRuntime: coupled cell has non-positive density/energy/temperature cache");
        }
        gas.density[k] = rho;
        gas.tpu[k] = tt / uu;
        const double vx = state.gas_cells.velocity_x_peculiar[row];
        const double vy = state.gas_cells.velocity_y_peculiar[row];
        const double vz = state.gas_cells.velocity_z_peculiar[row];
        gas.etot[k] = rho * (uu + 0.5 * (vx * vx + vy * vy + vz * vz));
        gas.px[k] = rho * vx;
        gas.py[k] = rho * vy;
        gas.pz[k] = rho * vz;
      }
    }
  }
  const GasWork gas_initial = gas;
  const double gas_energy_initial = [&] {
    double s = 0.0;
    for (std::size_t r = 0; r < rad.patches.size() && coupled; ++r) {
      const double v = rad.patches[r].cell_width_x_comoving * rad.patches[r].cell_width_y_comoving *
          rad.patches[r].cell_width_z_comoving;
      for (std::size_t k = cell_base[r]; k < cell_base[r + 1U]; ++k) s += gas.etot[k] * in.frame_end.properVolume(v);
    }
    return s;
  }();

  auto value_index = [&](std::size_t r, std::uint32_t g, std::size_t c0) {
    const std::size_t nc = cell_base[r + 1U] - cell_base[r];
    return static_cast<std::size_t>(rad.patches[r].value_offset) + static_cast<std::size_t>(g) * nc + c0;
  };
  // Implicit local solve U - alpha S(U) = U_rhs at scale factor a.
  auto solveLocal = [&](RadiationWorkState& w, GasWork& gw, double alpha, double a) -> bool {
    if (!coupled) return true;
    const double a4 = a * a * a * a;
    for (std::size_t r = 0; r < rad.patches.size(); ++r) {
      for (std::size_t c0 = 0; c0 < cell_base[r + 1U] - cell_base[r]; ++c0) {
        const std::size_t k = cell_base[r] + c0;
        RadiationMatterCell cell;
        cell.density = gw.density[k];
        cell.temperature_per_specific_energy = gw.tpu[k];
        cell.velocity = {gw.px[k] / gw.density[k], gw.py[k] / gw.density[k], gw.pz[k] / gw.density[k]};
        const double ke = 0.5 * (gw.px[k] * gw.px[k] + gw.py[k] * gw.py[k] + gw.pz[k] * gw.pz[k]) / gw.density[k];
        cell.specific_internal_energy = (gw.etot[k] - ke) / gw.density[k];
        if (!(cell.specific_internal_energy > 0.0)) return false;
        if (m_settings.coupling_model == RadiationCouplingModel::kGreyThermalReference) {
          const std::size_t vi = value_index(r, 0U, c0);
          double e = w.e[vi] / a4;
          std::array<double, 3> gv{w.gx[vi] / a4, w.gy[vi] / a4, w.gz[vi] / a4};
          const RadiationLocalResult res = solveGreyThermalImplicit(cell, e, gv, alpha, m_settings.coupling);
          report.newton_iterations += res.iterations;
          if (!res.converged) return false;
          w.e[vi] = e * a4;
          w.gx[vi] = gv[0] * a4;
          w.gy[vi] = gv[1] * a4;
          w.gz[vi] = gv[2] * a4;
          gw.etot[k] += res.gas_energy_change;
          gw.px[k] += res.gas_momentum_change[0];
          gw.py[k] += res.gas_momentum_change[1];
          gw.pz[k] += res.gas_momentum_change[2];
        } else {
          std::vector<double> e(ng);
          std::vector<std::array<double, 3>> gv(ng);
          for (std::uint32_t g = 0; g < ng; ++g) {
            const std::size_t vi = value_index(r, g, c0);
            e[g] = w.e[vi] / a4;
            gv[g] = {w.gx[vi] / a4, w.gy[vi] / a4, w.gz[vi] / a4};
          }
          const auto row = gw.row[k];
          if (state.gas_composition.enabled()) {
            const auto layout = physics::CompositionLayout::fromDescriptor(state.gas_composition.descriptor());
            const double mass = state.cells.mass_code[row];
            const double rho_cgs = gw.density[k] * m_settings.coupling.density_cgs_per_code;
            constexpr double k_amu_g = 1.66053906660e-24;
            const physics::ChemicalSpecies species[3] = {physics::ChemicalSpecies::kHI, physics::ChemicalSpecies::kHeI,
                                                         physics::ChemicalSpecies::kHeII};
            for (int s = 0; s < 3; ++s) {
              const double lane = state.gas_composition.value(layout.speciesLane(static_cast<std::size_t>(species[s])), row);
              const double amu = physics::chemicalSpeciesInfo(species[s]).mass_amu;
              cell.absorber_number_density_cgs[static_cast<std::size_t>(s)] = lane / mass * rho_cgs / (amu * k_amu_g);
            }
          }
          const RadiationLocalResult res = solvePhotoionizationImplicit(cell, m_settings.spectrum, e, gv, alpha,
                                                                        m_settings.coupling);
          for (std::uint32_t g = 0; g < ng; ++g) {
            const std::size_t vi = value_index(r, g, c0);
            w.e[vi] = e[g] * a4;
            w.gx[vi] = gv[g][0] * a4;
            w.gy[vi] = gv[g][1] * a4;
            w.gz[vi] = gv[g][2] * a4;
          }
          gw.px[k] += res.gas_momentum_change[0];
          gw.py[k] += res.gas_momentum_change[1];
          gw.pz[k] += res.gas_momentum_change[2];
        }
      }
    }
    return true;
  };
  auto add = [](RadiationWorkState& out, const RadiationWorkState& a, double wa, const RadiationWorkState& b, double wb) {
    for (std::size_t i = 0; i < out.size(); ++i) {
      out.e[i] = a.e[i] + wa * b.e[i];
      out.gx[i] = a.gx[i] + wa * b.gx[i];
      out.gy[i] = a.gy[i] + wa * b.gy[i];
      out.gz[i] = a.gz[i] + wa * b.gz[i];
    }
    (void)wb;
  };
  auto rateWithSources = [&](const RadiationWorkState& w, double inv_a, bool first_order, RadiationWorkState& rate,
                             RadiationStageDiagnostics& diag) {
    if (!solver.transportRate(w, inv_a, first_order, rate, diag)) return false;
    for (std::size_t i = 0; i < nvalues; ++i) rate.e[i] += source[i];
    return true;
  };

  const double total = in.dt_proper_code;
  double t = 0.0;
  std::uint32_t retries_left = m_settings.config.max_retries;
  double h_next = 0.0;
  RadiationWorkState t0, t1, tmp, ustar, unew;
  t0.resize(nvalues);
  t1.resize(nvalues);
  tmp.resize(nvalues);
  ustar.resize(nvalues);
  unew.resize(nvalues);
  if (workspace_reservation.pending()) workspace_reservation.commit();
  while (t < total * (1.0 - 1.0e-14)) {
    if (report.substeps >= m_settings.config.max_subcycles) {
      throw std::runtime_error("RadiationRuntime: physics.rt_max_subcycles exhausted before the interval end");
    }
    const double frac0 = t / total;
    const double a0 = in.scale_factor_begin + frac0 * (in.scale_factor_end - in.scale_factor_begin);
    double h = std::min(solver.stableTimeStep(m_settings.config.cfl, a0), total - t);
    if (h_next > 0.0) h = std::min(h, h_next);
    const double frac1 = (t + h) / total;
    const double a1 = in.scale_factor_begin + frac1 * (in.scale_factor_end - in.scale_factor_begin);
    const GasWork gas_saved = gas;
    RadiationStageDiagnostics diag;
    bool ok = rateWithSources(u, 1.0 / a0, false, t0, diag);
    add(tmp, u, h, t0, 0.0);
    GasWork gas_stage = gas;
    if (ok && m_settings.config.integrator == core::RtIntegratorConfig::kExplicitReference) {
      unew = tmp;
      ok = solveLocal(unew, gas_stage, h, a1);
    } else if (ok) {
      ustar = tmp;
      ok = solveLocal(ustar, gas_stage, h, a1);
      if (ok) {
        // S(U*) from the stage-1 solve: (U* - U_rhs)/h (radiation and gas).
        GasWork gas_s1 = gas_stage;
        ok = rateWithSources(ustar, 1.0 / a1, false, t1, diag);
        for (std::size_t i = 0; i < nvalues && ok; ++i) {
          const double s_e = (ustar.e[i] - tmp.e[i]) / h;
          const double s_x = (ustar.gx[i] - tmp.gx[i]) / h;
          const double s_y = (ustar.gy[i] - tmp.gy[i]) / h;
          const double s_z = (ustar.gz[i] - tmp.gz[i]) / h;
          unew.e[i] = u.e[i] + 0.5 * h * (t0.e[i] + t1.e[i]) + 0.5 * h * s_e;
          unew.gx[i] = u.gx[i] + 0.5 * h * (t0.gx[i] + t1.gx[i]) + 0.5 * h * s_x;
          unew.gy[i] = u.gy[i] + 0.5 * h * (t0.gy[i] + t1.gy[i]) + 0.5 * h * s_y;
          unew.gz[i] = u.gz[i] + 0.5 * h * (t0.gz[i] + t1.gz[i]) + 0.5 * h * s_z;
        }
        if (ok && coupled) {
          for (std::size_t k = 0; k < gas.etot.size(); ++k) {
            gas_stage.etot[k] = gas_saved.etot[k] + 0.5 * (gas_s1.etot[k] - gas_saved.etot[k]);
            gas_stage.px[k] = gas_saved.px[k] + 0.5 * (gas_s1.px[k] - gas_saved.px[k]);
            gas_stage.py[k] = gas_saved.py[k] + 0.5 * (gas_s1.py[k] - gas_saved.py[k]);
            gas_stage.pz[k] = gas_saved.pz[k] + 0.5 * (gas_s1.pz[k] - gas_saved.pz[k]);
          }
        }
        if (ok) ok = solveLocal(unew, gas_stage, 0.5 * h, a1);
      }
    }
    std::uint64_t bad = 0U;
    if (ok && !solver.admissible(unew, bad)) {
      // Whole-stage first-order GLF fallback (explicit part), then re-solve.
      ++report.first_order_fallbacks;
      gas_stage = gas_saved;
      RadiationStageDiagnostics fdiag;
      ok = rateWithSources(u, 1.0 / a0, true, t0, fdiag);
      if (ok) {
        add(unew, u, h, t0, 0.0);
        ok = solveLocal(unew, gas_stage, h, a1) && solver.admissible(unew, bad);
      }
    }
    if (!ok) {
      if (retries_left == 0U) {
        std::ostringstream msg;
        msg << "RadiationRuntime: substep failed after " << m_settings.config.max_retries
            << " retries (inadmissible cells " << bad << ", local implicit solve or transport failure); nothing "
            << "was committed";
        throw std::runtime_error(msg.str());
      }
      --retries_left;
      ++report.retries;
      gas = gas_saved;
      h_next = 0.5 * h;
      continue;
    }
    // Cosmological edge transfer (exact exponential, H at the substep midpoint).
    const double hubble = in.hubble_begin_code + 0.5 * (frac0 + frac1) * (in.hubble_end_code - in.hubble_begin_code);
    if (hubble > 0.0 && !m_settings.spectrum.grey) {
      RadiationRedshiftLedger ledger;
      std::vector<double> e(ng), gx(ng), gy(ng), gz(ng);
      for (std::size_t r = 0; r < rad.patches.size(); ++r) {
        const double vol = rad.patches[r].cell_width_x_comoving * rad.patches[r].cell_width_y_comoving *
            rad.patches[r].cell_width_z_comoving;
        for (std::size_t c0 = 0; c0 < cell_base[r + 1U] - cell_base[r]; ++c0) {
          for (std::uint32_t g = 0; g < ng; ++g) {
            const std::size_t vi = value_index(r, g, c0);
            e[g] = unew.e[vi];
            gx[g] = unew.gx[vi];
            gy[g] = unew.gy[vi];
            gz[g] = unew.gz[vi];
          }
          applyFrequencyRedshiftCell(m_settings.spectrum, hubble, h, e, gx, gy, gz, vol, ledger);
          for (std::uint32_t g = 0; g < ng; ++g) {
            const std::size_t vi = value_index(r, g, c0);
            unew.e[vi] = e[g];
            unew.gx[vi] = gx[g];
            unew.gy[vi] = gy[g];
            unew.gz[vi] = gz[g];
          }
        }
      }
      report.redshift_leaked_comoving_energy += ledger.leaked_comoving_energy;
    }
    for (std::size_t r = 0; r < rad.patches.size(); ++r) {
      const std::size_t nc = cell_base[r + 1U] - cell_base[r];
      const double vol = rad.patches[r].cell_width_x_comoving * rad.patches[r].cell_width_y_comoving *
          rad.patches[r].cell_width_z_comoving;
      for (std::uint32_t g = 0; g < ng; ++g) {
        for (std::size_t c0 = 0; c0 < nc; ++c0) report.emitted_comoving_energy += source[value_index(r, g, c0)] * vol * h;
      }
    }
    report.boundary_outflow_comoving_energy += diag.boundary_energy_outflow_code * h;
    report.hll_fallback_faces += diag.hll_fallback_faces;
    u = unew;
    gas = gas_stage;
    t += h;
    h_next = 0.0;
    ++report.substeps;
  }

  // Ledger: radiation energy change and RSLA invariant (isolated frames).
  double e_before = 0.0;
  double e_after = 0.0;
  for (std::size_t r = 0; r < rad.patches.size(); ++r) {
    const double vol = rad.patches[r].cell_width_x_comoving * rad.patches[r].cell_width_y_comoving *
        rad.patches[r].cell_width_z_comoving;
    const std::size_t nc = cell_base[r + 1U] - cell_base[r];
    for (std::uint32_t g = 0; g < ng; ++g) {
      for (std::size_t c0 = 0; c0 < nc; ++c0) {
        const std::size_t vi = value_index(r, g, c0);
        e_before += rad.energy_density_comoving[vi] * vol;
        e_after += u.e[vi] * vol;
      }
    }
  }
  report.radiation_energy_change_comoving = e_after - e_before;

  // Stage the atomic commit: gas exchange through a transaction, then lanes.
  physics::SourceTransactionBatch batch(physics::TransactionFloors{1.0e-300, 0.0});
  if (coupled) {
    physics::SourceTransactionPlan plan;
    plan.id = physics::SourceTransactionId{0U, in.interval_index, physics::SourceChannel::kRtMatterExchange, 0U};
    plan.topology_epoch = state.cellIndexGeneration();
    double gas_energy_final = 0.0;
    for (std::size_t r = 0; r < rad.patches.size(); ++r) {
      const double vol = in.frame_end.properVolume(rad.patches[r].cell_width_x_comoving *
                                                   rad.patches[r].cell_width_y_comoving *
                                                   rad.patches[r].cell_width_z_comoving);
      for (std::size_t k = cell_base[r]; k < cell_base[r + 1U]; ++k) {
        physics::GasCellIncrement inc;
        inc.local_row = gas.row[k];
        inc.gas_cell_id = state.gas_cells.gas_cell_id[gas.row[k]];
        inc.delta.total_energy_code = (gas.etot[k] - gas_initial.etot[k]) * vol;
        inc.delta.momentum_code = {(gas.px[k] - gas_initial.px[k]) * vol, (gas.py[k] - gas_initial.py[k]) * vol,
                                   (gas.pz[k] - gas_initial.pz[k]) * vol};
        gas_energy_final += gas.etot[k] * vol;
        if (inc.delta.total_energy_code != 0.0 || inc.delta.momentum_code[0] != 0.0 ||
            inc.delta.momentum_code[1] != 0.0 || inc.delta.momentum_code[2] != 0.0) {
          plan.cells.push_back(inc);
        }
      }
    }
    physics::ConservedIncrement total_gas{};
    for (const auto& cinc : plan.cells) total_gas += cinc.delta;
    plan.source_increment = total_gas.negated();  // supplied by the radiation field
    report.gas_energy_exchange_code = gas_energy_final - gas_energy_initial;
    if (in.scale_factor_begin == in.scale_factor_end && in.hubble_begin_code == 0.0 &&
        in.hubble_end_code == 0.0) {
      const double a = in.scale_factor_end;
      const double vol_scale = a * a * a / (a * a * a * a);  // comoving E V_comov -> physical E V_phys
      const double invariant_change = report.gas_energy_exchange_code +
          (c / ch) * vol_scale * (report.radiation_energy_change_comoving - report.emitted_comoving_energy +
                                  report.boundary_outflow_comoving_energy);
      const double scale = std::max({std::abs(report.gas_energy_exchange_code),
                                     (c / ch) * vol_scale * std::abs(e_after) + 1.0e-300});
      report.rsla_invariant_residual = std::abs(invariant_change) / scale;
    }
    if (!plan.cells.empty()) batch.add(std::move(plan));
  }
  if (batch.size() > 0U) {
    if (!batch.prepare(state) || !batch.validate() || !batch.authorize(true)) {
      throw std::runtime_error("RadiationRuntime: matter exchange commit rejected: " + batch.report().failure);
    }
  }
  // Both validated: commit radiation lanes and gas together.
  for (std::size_t i = 0; i < nvalues; ++i) {
    rad.energy_density_comoving[i] = u.e[i];
    rad.flux_x_comoving[i] = u.gx[i] * c / ch;
    rad.flux_y_comoving[i] = u.gy[i] * c / ch;
    rad.flux_z_comoving[i] = u.gz[i] * c / ch;
  }
  rad.scale_factor = in.scale_factor_end;
  if (batch.size() > 0U) report.touched_gas_rows = batch.commit(state);
  return report;
}

std::uint64_t radiationIntervalWorkspaceBytes(std::size_t nvalues, std::size_t ncells, bool coupled) {
  const std::uint64_t work = (6U * 4U + 1U) * static_cast<std::uint64_t>(nvalues) * sizeof(double);
  const std::uint64_t gas = coupled ? 3U * static_cast<std::uint64_t>(ncells) * (sizeof(std::uint32_t) + 6U * sizeof(double))
                                    : 0U;
  return work + gas;
}

}  // namespace cosmosim::radiation
