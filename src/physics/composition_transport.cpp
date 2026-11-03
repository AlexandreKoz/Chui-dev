#include "cosmosim/physics/composition_transport.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/baryon_composition.hpp"

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

std::uint32_t CompositionFluxAccumulator::resolveRow(std::size_t storage_row, bool& is_real) const {
  is_real = false;
  if (storage_row == hydro::k_invalid_cell_index) return k_unresolved;
  if (storage_row < m_geometry->cellCount()) {
    is_real = true;
    return m_dense_row_by_real_cell[storage_row];
  }
  for (std::size_t slot = 0; slot < m_geometry->ghost_cells.size(); ++slot) {
    const hydro::HydroGhostCell& ghost = m_geometry->ghost_cells[slot];
    if (ghost.ghost_cell != storage_row) continue;
    if (ghost.boundary_kind == hydro::HydroBoundaryKind::kImportedMpi) {
      throw std::runtime_error(
          "composition transport: face touches an imported MPI ghost; distributed composition ghost "
          "exchange is a Part 3 deliverable and the stage is refused (no species may be dropped)");
    }
    std::uint64_t origin = slot < m_ghost_origin_by_slot.size() ? m_ghost_origin_by_slot[slot] : 0U;
    if (origin == 0U) origin = ghost.origin_gas_cell_id;
    if (origin != 0U) {
      const auto row = m_state.rowForGasCellId(origin);
      if (!row.has_value()) {
        throw std::runtime_error(
            "composition transport: ghost source gas_cell_id " + std::to_string(origin) +
            " is not local; remote AMR composition ghosts are Part 3 scope (fail closed)");
      }
      return *row;
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
  record.owner_row = resolveRow(owner_cell, owner_real);
  // neighbor_cell is the conservative TARGET (periodic wrap already
  // resolved by the solver). When it is invalid, the donor for inflow is the
  // face's stencil neighbour (ghost) — resolve it for the donor fraction.
  const std::size_t donor_neighbor =
      neighbor_cell != hydro::k_invalid_cell_index ? neighbor_cell : face.neighbor_cell;
  record.neighbor_row = resolveRow(donor_neighbor, neighbor_real);
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
    if (!coarse_row.has_value()) {
      throw std::runtime_error(
          "composition transport: coarse register cell is not local; distributed composition reflux "
          "is Part 3 scope (fail closed)");
    }
    if (owner_real) {
      record.neighbor_row = *coarse_row;
      record.update_neighbor = 1U;
    } else {
      record.owner_row = *coarse_row;
      record.update_owner = 1U;
    }
  }
  if (record.update_owner == 0U && record.update_neighbor == 0U) return;
  const std::uint32_t donor = mass_transfer_code >= 0.0 ? record.owner_row : record.neighbor_row;
  if (donor == k_unresolved) {
    throw std::runtime_error(
        "composition transport: donor cell of a face could not be resolved to local stable identity");
  }
  m_records.push_back(record);
}

void CompositionFluxAccumulator::prepare(const core::SimulationState& state) {
  const auto& sidecar = state.gas_composition;
  m_prepared = false;
  m_delta.clear();
  if (!sidecar.enabled()) {
    m_records.clear();
    m_prepared = true;
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
  for (const Record& record : m_records) {
    const std::uint32_t donor = record.mass >= 0.0 ? record.owner_row : record.neighbor_row;
    if (donor >= rows ||
        (record.update_owner != 0U && record.owner_row >= rows) ||
        (record.update_neighbor != 0U && record.neighbor_row >= rows)) {
      throw std::runtime_error("composition transport: record row outside the composition sidecar");
    }
    const double donor_mass = gas_mass_of(donor);
    if (!(donor_mass > 0.0)) {
      throw std::runtime_error(
          "composition transport: donor cell has no composition mass but hydro moved mass out of it");
    }
    for (std::size_t l = 0; l < lanes; ++l) {
      const double moved = record.mass * sidecar.value(l, donor) / donor_mass;
      if (record.update_owner != 0U) m_delta[l * rows + record.owner_row] -= moved;
      if (record.update_neighbor != 0U) m_delta[l * rows + record.neighbor_row] += moved;
    }
  }
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
  m_records.clear();
  m_prepared = true;
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
