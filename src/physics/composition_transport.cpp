#include "cosmosim/physics/composition_transport.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/baryon_composition.hpp"

#include <algorithm>
#include <map>

namespace cosmosim::physics {

CompositionFluxAccumulator::CompositionFluxAccumulator(const core::SimulationState& state)
    : m_state(state) {}

void CompositionFluxAccumulator::beginAdvance(
    const hydro::HydroPatchGeometry& geometry,
    std::span<const std::uint32_t> dense_row_by_real_cell,
    std::span<const std::uint64_t> ghost_origin_gas_cell_id_by_slot) {
  if (dense_row_by_real_cell.size() < geometry.cellCount()) {
    throw std::invalid_argument(
        "CompositionFluxAccumulator: dense-row map does not cover every real geometry row");
  }
  m_geometry = &geometry;
  m_dense_row_by_real_cell = dense_row_by_real_cell;
  m_ghost_origin_by_slot = ghost_origin_gas_cell_id_by_slot;
}

void CompositionFluxAccumulator::endAdvance() noexcept {
  m_geometry = nullptr;
  m_dense_row_by_real_cell = {};
  m_ghost_origin_by_slot = {};
}

std::uint32_t CompositionFluxAccumulator::resolveRow(std::size_t storage_row, bool& is_real,
                                                     std::uint64_t& remote_id) const {
  is_real = false;
  remote_id = 0U;
  if (storage_row == hydro::k_invalid_cell_index) return k_unresolved;
  if (storage_row < m_geometry->cellCount()) {
    is_real = true;
    return m_dense_row_by_real_cell[storage_row];
  }
  for (std::size_t slot = 0; slot < m_geometry->ghost_cells.size(); ++slot) {
    const hydro::HydroGhostCell& ghost = m_geometry->ghost_cells[slot];
    if (ghost.ghost_cell != storage_row) continue;
    std::uint64_t origin = slot < m_ghost_origin_by_slot.size() ? m_ghost_origin_by_slot[slot] : 0U;
    if (origin == 0U) origin = ghost.origin_gas_cell_id;
    if (origin != 0U) {
      const auto row = m_state.rowForGasCellId(origin);
      if (!row.has_value()) {
        // CHUI Part 3: remote ghost source, resolved by stable id through the
        // donor-fraction exchange before prepare (fail closed there if absent).
        remote_id = origin;
        return k_unresolved;
      }
      return *row;
    }
    if (ghost.boundary_kind == hydro::HydroBoundaryKind::kImportedMpi) {
      throw std::runtime_error(
          "composition transport: imported MPI ghost has no stable source gas_cell_id (fail closed)");
    }
    // Open/reflective physical boundary ghosts copy the owner state: the
    // inflow composition is the owner's (zero-gradient boundary).
    if (ghost.owner_real_cell < m_geometry->cellCount()) {
      return m_dense_row_by_real_cell[ghost.owner_real_cell];
    }
    return k_unresolved;
  }
  return k_unresolved;
}

void CompositionFluxAccumulator::recordFaceMassTransfer(
    std::size_t face_index, std::size_t owner_cell, std::size_t neighbor_cell,
    double mass_transfer_code) {
  if (m_geometry == nullptr) {
    throw std::logic_error("CompositionFluxAccumulator: face recorded outside beginAdvance/endAdvance");
  }
  if (mass_transfer_code == 0.0) return;
  if (!std::isfinite(mass_transfer_code)) {
    throw std::runtime_error("composition transport: non-finite face mass transfer");
  }
  const hydro::HydroFace& face = m_geometry->faces.at(face_index);
  hydro::HydroFluxRegisterFaceRole role = hydro::HydroFluxRegisterFaceRole::kNone;
  std::uint64_t coarse_gas_cell_id = 0U;
  if (face_index < m_geometry->flux_register_faces.size()) {
    role = m_geometry->flux_register_faces[face_index].role;
    coarse_gas_cell_id = m_geometry->flux_register_faces[face_index].coarse_gas_cell_id;
  }
  bool owner_real = false;
  bool neighbor_real = false;
  Record record;
  record.mass = mass_transfer_code;
  record.owner_row = resolveRow(owner_cell, owner_real, record.owner_remote_id);
  // neighbor_cell is the conservative TARGET (periodic wrap already
  // resolved by the solver). When it is invalid, the donor for inflow is the
  // face's stencil neighbour (ghost) — resolve it for the donor fraction.
  const std::size_t donor_neighbor =
      neighbor_cell != hydro::k_invalid_cell_index ? neighbor_cell : face.neighbor_cell;
  record.neighbor_row = resolveRow(donor_neighbor, neighbor_real, record.neighbor_remote_id);
  if (neighbor_cell == hydro::k_invalid_cell_index) neighbor_real = false;

  record.update_owner = owner_real ? 1U : 0U;
  record.update_neighbor = neighbor_real ? 1U : 0U;
  if (role == hydro::HydroFluxRegisterFaceRole::kCoarse) {
    // Coarse side of a coarse/fine interface: fine faces own the update.
    record.update_owner = 0U;
    record.update_neighbor = 0U;
  } else if (role == hydro::HydroFluxRegisterFaceRole::kFine) {
    // Fine side: mirror the transfer onto the coarse cell.
    const auto coarse_row = m_state.rowForGasCellId(coarse_gas_cell_id);
    if (!coarse_row.has_value() && coarse_gas_cell_id == 0U) {
      throw std::runtime_error(
          "composition transport: coarse register cell has no stable identity (fail closed)");
    }
    // A remote coarse cell (CHUI Part 3) receives its mirrored update as an
    // owner-targeted correction; the donor fraction for inflow from it comes
    // from the donor exchange (same stable id).
    const std::uint32_t coarse = coarse_row.has_value() ? *coarse_row : k_unresolved;
    const std::uint64_t coarse_remote = coarse_row.has_value() ? 0U : coarse_gas_cell_id;
    if (owner_real) {
      record.neighbor_row = coarse;
      record.neighbor_remote_id = coarse_remote;
      record.update_neighbor = 1U;
    } else {
      record.owner_row = coarse;
      record.owner_remote_id = coarse_remote;
      record.update_owner = 1U;
    }
  }
  if (record.update_owner == 0U && record.update_neighbor == 0U) return;
  const std::uint32_t donor = mass_transfer_code >= 0.0 ? record.owner_row : record.neighbor_row;
  const std::uint64_t donor_remote = mass_transfer_code >= 0.0 ? record.owner_remote_id : record.neighbor_remote_id;
  if (donor == k_unresolved && donor_remote == 0U) {
    throw std::runtime_error(
        "composition transport: donor cell of a face could not be resolved to local stable identity");
  }
  m_records.push_back(record);
}

std::vector<std::uint64_t> CompositionFluxAccumulator::requiredRemoteDonorIds() const {
  std::vector<std::uint64_t> ids;
  for (const Record& record : m_records) {
    if (record.owner_remote_id != 0U) ids.push_back(record.owner_remote_id);
    if (record.neighbor_remote_id != 0U) ids.push_back(record.neighbor_remote_id);
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  return ids;
}

void CompositionFluxAccumulator::setRemoteDonorFractions(
    std::span<const std::uint64_t> sorted_ids, std::span<const double> fractions_row_major, std::size_t lanes) {
  if (!std::is_sorted(sorted_ids.begin(), sorted_ids.end()) ||
      std::adjacent_find(sorted_ids.begin(), sorted_ids.end()) != sorted_ids.end() ||
      fractions_row_major.size() != sorted_ids.size() * lanes) {
    throw std::invalid_argument("composition transport: remote donor table is unsorted, duplicated or mis-sized");
  }
  m_remote_ids.assign(sorted_ids.begin(), sorted_ids.end());
  m_remote_fractions.assign(fractions_row_major.begin(), fractions_row_major.end());
  m_remote_lanes = lanes;
}

const double* CompositionFluxAccumulator::remoteFractions(std::uint64_t gas_cell_id) const {
  const auto it = std::lower_bound(m_remote_ids.begin(), m_remote_ids.end(), gas_cell_id);
  if (it == m_remote_ids.end() || *it != gas_cell_id) {
    return nullptr;
  }
  return m_remote_fractions.data() + static_cast<std::size_t>(it - m_remote_ids.begin()) * m_remote_lanes;
}

void CompositionFluxAccumulator::accumulateDeltas(const core::SimulationState& state) {
  const auto& sidecar = state.gas_composition;
  m_prepared = false;
  m_delta.clear();
  m_outbound.clear();
  if (!sidecar.enabled()) {
    m_records.clear();
    return;
  }
  sidecar.requireAlignedWith(state.gas_cells.gas_cell_id, "CompositionFluxAccumulator::prepare");
  const auto& d = sidecar.descriptor();
  const std::size_t lanes = sidecar.laneCount();
  const std::size_t rows = sidecar.rowCount();
  // Deltas are accumulated against the UNMODIFIED step-start sidecar, so the
  // result is independent of patch / record order.
  m_delta.assign(lanes * rows, 0.0);
  const auto gas_mass_of = [&](std::uint32_t row) {
    double m = 0.0;
    for (std::uint32_t l = 0; l < d.element_lane_count; ++l) m += sidecar.value(d.element_lane_begin + l, row);
    for (std::uint32_t l = 0; l < d.dust_lane_count; ++l) m += sidecar.value(d.dust_lane_begin + l, row);
    return m;
  };
  // Remote corrections accumulated by stable id (sorted map => deterministic).
  std::map<std::uint64_t, std::vector<double>> outbound;
  std::vector<double> fraction(lanes, 0.0);
  for (const Record& record : m_records) {
    const bool donor_is_owner = record.mass >= 0.0;
    const std::uint32_t donor = donor_is_owner ? record.owner_row : record.neighbor_row;
    const std::uint64_t donor_remote = donor_is_owner ? record.owner_remote_id : record.neighbor_remote_id;
    if (donor_remote != 0U) {
      const double* remote = remoteFractions(donor_remote);
      if (remote == nullptr || m_remote_lanes != lanes) {
        throw std::runtime_error(
            "composition transport: remote donor " + std::to_string(donor_remote) +
            " has no exchanged composition fractions (fail closed; no species may be dropped)");
      }
      std::copy(remote, remote + lanes, fraction.begin());
    } else {
      if (donor >= rows) {
        throw std::runtime_error("composition transport: record row outside the composition sidecar");
      }
      const double donor_mass = gas_mass_of(donor);
      if (!(donor_mass > 0.0)) {
        throw std::runtime_error(
            "composition transport: donor cell has no composition mass but hydro moved mass out of it");
      }
      for (std::size_t l = 0; l < lanes; ++l) fraction[l] = sidecar.value(l, donor) / donor_mass;
    }
    const auto apply = [&](bool update, std::uint32_t row, std::uint64_t remote_id, double sign) {
      if (!update) return;
      if (remote_id != 0U) {
        auto& delta = outbound[remote_id];
        delta.resize(lanes, 0.0);
        for (std::size_t l = 0; l < lanes; ++l) delta[l] += sign * record.mass * fraction[l];
        return;
      }
      if (row >= rows) {
        throw std::runtime_error("composition transport: record row outside the composition sidecar");
      }
      for (std::size_t l = 0; l < lanes; ++l) m_delta[l * rows + row] += sign * record.mass * fraction[l];
    };
    apply(record.update_owner != 0U, record.owner_row, record.owner_remote_id, -1.0);
    apply(record.update_neighbor != 0U, record.neighbor_row, record.neighbor_remote_id, +1.0);
  }
  m_outbound.reserve(outbound.size());
  for (auto& [id, delta] : outbound) {
    m_outbound.push_back(RemoteCorrection{.gas_cell_id = id, .lane_delta_code = std::move(delta)});
  }
  m_records.clear();
}

void CompositionFluxAccumulator::addInboundCorrection(
    const core::SimulationState& state, std::uint64_t gas_cell_id, std::span<const double> lane_delta_code) {
  const auto& sidecar = state.gas_composition;
  if (!sidecar.enabled()) {
    throw std::logic_error("composition transport: inbound correction without composition sidecar");
  }
  const auto row = state.rowForGasCellId(gas_cell_id);
  if (!row.has_value() || lane_delta_code.size() != sidecar.laneCount()) {
    throw std::runtime_error(
        "composition transport: inbound remote correction targets a non-local cell or has a lane mismatch (fail closed)");
  }
  const std::size_t rows = sidecar.rowCount();
  for (std::size_t l = 0; l < lane_delta_code.size(); ++l) {
    if (!std::isfinite(lane_delta_code[l])) {
      throw std::runtime_error("composition transport: non-finite inbound remote correction");
    }
    m_delta[l * rows + *row] += lane_delta_code[l];
  }
}

void CompositionFluxAccumulator::validatePrepared(const core::SimulationState& state) {
  const auto& sidecar = state.gas_composition;
  if (!sidecar.enabled()) {
    m_prepared = true;
    return;
  }
  const std::size_t lanes = sidecar.laneCount();
  const std::size_t rows = sidecar.rowCount();
  for (std::size_t l = 0; l < lanes; ++l) {
    const auto values = sidecar.lane(l);
    for (std::size_t row = 0; row < rows; ++row) {
      const double updated = values[row] + m_delta[l * rows + row];
      const double scale = std::max(state.cells.mass_code[row], 1.0e-300);
      if (updated < -1.0e-12 * scale || !std::isfinite(updated)) {
        throw std::runtime_error(
            "composition transport produced a negative/non-finite lane " + std::to_string(l) +
            " at row " + std::to_string(row) + "; hydro CFL/positivity contract violated");
      }
    }
  }
  m_prepared = true;
}

void CompositionFluxAccumulator::prepareDeferred(const core::SimulationState& state) {
  accumulateDeltas(state);
}

void CompositionFluxAccumulator::prepare(const core::SimulationState& state) {
  accumulateDeltas(state);
  if (!m_outbound.empty()) {
    throw std::runtime_error(
        "composition transport: remote corrections pending but the single-rank prepare() was used (fail closed)");
  }
  validatePrepared(state);
}

CompositionTransportCommitReport CompositionFluxAccumulator::commit(core::SimulationState& state) noexcept {
  CompositionTransportCommitReport report;
  auto& sidecar = state.gas_composition;
  if (!m_prepared || !sidecar.enabled() || m_delta.empty()) {
    m_prepared = false;
    return report;
  }
  const auto& d = sidecar.descriptor();
  const std::size_t rows = sidecar.rowCount();
  for (std::size_t l = 0; l < sidecar.laneCount(); ++l) {
    auto values = sidecar.lane(l);
    for (std::size_t row = 0; row < rows; ++row) {
      values[row] = std::max(values[row] + m_delta[l * rows + row], 0.0);
    }
  }
  // Reconcile with the hydro cell mass. Differences arise only from the
  // hydro density floor (mass created) or round-off; the correction adds or
  // removes mass WITH THE CELL'S OWN COMPOSITION and is ledgered.
  for (std::size_t row = 0; row < rows; ++row) {
    double comp = 0.0;
    for (std::uint32_t l = 0; l < d.element_lane_count; ++l) comp += sidecar.value(d.element_lane_begin + l, row);
    for (std::uint32_t l = 0; l < d.dust_lane_count; ++l) comp += sidecar.value(d.dust_lane_begin + l, row);
    const double mass = state.cells.mass_code[row];
    if (!(comp > 0.0) || !(mass > 0.0)) continue;
    const double rel = std::abs(comp - mass) / mass;
    report.max_relative_mismatch = std::max(report.max_relative_mismatch, rel);
    if (rel > 0.0) {
      const double factor = mass / comp;
      for (std::size_t l = 0; l < sidecar.laneCount(); ++l) sidecar.lane(l)[row] *= factor;
      report.reconciled_mass_code += mass - comp;
      if (rel > 1.0e-10) ++report.reconciled_rows;
    }
  }
  m_delta.clear();
  m_prepared = false;
  return report;
}

}  // namespace cosmosim::physics
