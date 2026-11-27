// CHUI-EXT-PHYS-REF-v1 read-only science accessor (source only, NOT compiled
// or executed in pass 01).

#include "cosmosim/observables/gas_science_view.hpp"

#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::observables {
namespace {

namespace k = core::constants;

// Persistent species slot order (physics::ChemicalSpecies, layout v1) and the
// particle masses in atomic mass units of the slots used here.
constexpr std::size_t k_slot_hi = 0;
constexpr std::size_t k_slot_hii = 1;
constexpr std::size_t k_slot_hm = 2;
constexpr std::size_t k_slot_h2p = 4;
constexpr std::size_t k_slot_heii = 6;
constexpr std::size_t k_slot_heiii = 7;
constexpr double k_amu_h = 1.00782503223;
constexpr double k_amu_he = 4.00260325413;
constexpr double k_amu_e = 5.48579909065e-4;

[[nodiscard]] double slotMassAmu(std::size_t slot) {
  switch (slot) {
    case 0: return k_amu_h;                    // HI
    case 1: return k_amu_h - k_amu_e;          // HII
    case 2: return k_amu_h + k_amu_e;          // H-
    case 3: return 2.0 * k_amu_h;              // H2
    case 4: return 2.0 * k_amu_h - k_amu_e;    // H2+
    case 5: return k_amu_he;                   // HeI
    case 6: return k_amu_he - k_amu_e;         // HeII
    case 7: return k_amu_he - 2.0 * k_amu_e;   // HeIII
    default: break;
  }
  throw std::out_of_range("species slot outside the hydrogen/helium network");
}

}  // namespace

GasScienceAccessor::GasScienceAccessor(const core::SimulationState& state, const core::UnitSystem& units,
                                       double scale_factor)
    : m_state(state), m_units(units), m_a(scale_factor), m_rows(state.cells.size()) {
  m_units.validate();
  if (!(scale_factor > 0.0)) throw std::invalid_argument("science accessor needs a > 0");
  if (state.gas_cells.size() != m_rows) throw std::invalid_argument("gas sidecar rows differ from cell rows");
  const auto& d = state.gas_composition.descriptor();
  if (state.gas_composition.enabled() && d.species_lane_count >= 8U && state.gas_composition.rowCount() == m_rows) {
    m_species_begin = d.species_lane_begin;
  }
  m_row_volume_comoving_code.assign(m_rows, 0.0);
  const core::PatchSoa& patches = state.patches;
  for (std::size_t p = 0; p < patches.size(); ++p) {
    if (patches.cell_count[p] == 0U) continue;
    if (patches.cell_dim_x[p] == 0U || patches.extent_x_comoving[p] == 0.0) {
      throw std::invalid_argument("science accessor: patch with legacy undeclared geometry (no cell volume)");
    }
    const double volume = patches.extent_x_comoving[p] / patches.cell_dim_x[p] * patches.extent_y_comoving[p] /
                          patches.cell_dim_y[p] * patches.extent_z_comoving[p] / patches.cell_dim_z[p];
    for (std::uint32_t c = 0; c < patches.cell_count[p]; ++c) {
      const std::size_t row = static_cast<std::size_t>(patches.first_cell[p]) + c;
      if (row >= m_rows) throw std::out_of_range("patch cell range exceeds the cell rows");
      m_row_volume_comoving_code[row] = volume;
    }
  }
}

double GasScienceAccessor::properVolume(std::size_t row) const {
  const double v = m_row_volume_comoving_code.at(row);
  if (!(v > 0.0)) throw std::runtime_error("cell row has no leaf-patch volume");
  const double l = m_units.length_si_per_code;
  return v * l * l * l * m_a * m_a * m_a;
}

spectral::ProviderValue GasScienceAccessor::speciesDensity(std::size_t row, std::size_t slot) const {
  if (m_species_begin == k_no_species) {
    return spectral::ProviderValue::unavailable("no species-resolved composition layout (>= 8 species lanes) in state");
  }
  const double mass_code = m_state.gas_composition.value(m_species_begin + slot, row);
  if (!(mass_code >= 0.0)) return spectral::ProviderValue::unavailable("negative species mass lane");
  if (mass_code == 0.0) return spectral::ProviderValue::physicalZero();
  const double mass_kg = mass_code * m_units.mass_si_per_code;
  return spectral::ProviderValue::available(mass_kg / (slotMassAmu(slot) * k::k_atomic_mass_unit_si) / properVolume(row));
}

spectral::ProviderValue GasScienceAccessor::neutralHydrogenDensity(std::size_t row) const {
  return speciesDensity(row, k_slot_hi);
}

spectral::ProviderValue GasScienceAccessor::protonDensity(std::size_t row) const {
  return speciesDensity(row, k_slot_hii);
}

spectral::ProviderValue GasScienceAccessor::electronDensity(std::size_t row) const {
  const auto hii = speciesDensity(row, k_slot_hii);
  if (!hii.known()) return hii;
  const auto hm = speciesDensity(row, k_slot_hm);
  const auto h2p = speciesDensity(row, k_slot_h2p);
  const auto heii = speciesDensity(row, k_slot_heii);
  const auto heiii = speciesDensity(row, k_slot_heiii);
  const double ne = hii.value + h2p.value + heii.value + 2.0 * heiii.value - hm.value;
  if (ne < 0.0) return spectral::ProviderValue::unavailable("charge neutrality gives n_e < 0 (inconsistent species)");
  return ne == 0.0 ? spectral::ProviderValue::physicalZero() : spectral::ProviderValue::available(ne);
}

spectral::ProviderValue GasScienceAccessor::ionZ2Density(std::size_t row) const {
  const auto hii = speciesDensity(row, k_slot_hii);
  if (!hii.known()) return hii;
  const double value = hii.value + speciesDensity(row, k_slot_h2p).value + speciesDensity(row, k_slot_heii).value +
                       4.0 * speciesDensity(row, k_slot_heiii).value;
  return value == 0.0 ? spectral::ProviderValue::physicalZero() : spectral::ProviderValue::available(value);
}

spectral::ProviderValue GasScienceAccessor::temperature(std::size_t row) const {
  const double t = m_state.gas_cells.temperature_code.at(row);
  if (!(t > 0.0) || !std::isfinite(t)) return spectral::ProviderValue::unavailable("temperature cache not positive");
  return spectral::ProviderValue::available(t);
}

std::array<double, 3> GasScienceAccessor::peculiarVelocity(std::size_t row) const {
  const double v = m_units.velocity_si_per_code;
  return {m_state.gas_cells.velocity_x_peculiar.at(row) * v, m_state.gas_cells.velocity_y_peculiar.at(row) * v,
          m_state.gas_cells.velocity_z_peculiar.at(row) * v};
}

void GasScienceAccessor::attachMagneticField(std::span<const double> bx, std::span<const double> by,
                                             std::span<const double> bz) {
  if (bx.size() != m_rows || by.size() != m_rows || bz.size() != m_rows) {
    throw std::invalid_argument("magnetic field samples must match the cell rows");
  }
  m_bx = bx;
  m_by = by;
  m_bz = bz;
}

std::array<double, 3> GasScienceAccessor::magneticField(std::size_t row) const {
  if (m_bx.empty()) throw std::runtime_error("magnetic field unavailable (no CT authority sample attached)");
  return {m_bx[row], m_by[row], m_bz[row]};
}

std::string GasScienceAccessor::provenance() const {
  std::ostringstream text;
  text << "gas_science_accessor_v1;a=" << m_a << ";species=" << (hasSpecies() ? "composition_lanes" : "unavailable")
       << ";magnetic=" << (hasMagneticField() ? "attached" : "unavailable") << ";electrons=h_he_charge_neutrality";
  return text.str();
}

}  // namespace cosmosim::observables
