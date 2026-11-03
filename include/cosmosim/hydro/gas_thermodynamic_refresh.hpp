#pragma once

#include <cstdint>
#include <span>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {

// Single thermodynamic invalidation/refresh contract (audit F4).
//
// Authority. For a gas cell the authoritative thermodynamic inputs are
//   SimulationState::gas_cells.density_code        (comoving mass density)
//   SimulationState::gas_cells.internal_energy_code (specific internal energy u)
//   SimulationState::gas_cells.velocity_*_peculiar  (peculiar velocity)
//   SimulationState::gas_cells.metal_mass_code       (composition)
// together with the selected closure. pressure_code, sound_speed_code and
// temperature_code are DERIVED caches. Any stage that changes an
// authoritative input outside the hydro primitive store (feedback,
// enrichment, accretion, regrid, cooling outside hydro) must call
// SimulationState::notifyGasThermodynamicMutation() and then this refresh
// before the next CFL decision, output or reconstruction.
//
// Closure. Ideal gas: P = (gamma-1) rho u, c^2 = gamma P / rho. If a
// HydroThermodynamicClosure is supplied (effective multiphase ISM), it is
// evaluated at (scale_factor, redshift) from the comoving conserved state and
// overrides P and c^2 when valid (including branch transitions caused by the
// mutation, e.g. heated gas leaving the effective EOS).
//
// Floors. Density and pressure floors are applied to the DERIVED pressure
// only; the authoritative density/internal energy are never modified here,
// so this refresh cannot create or destroy mass or energy. Floor activations
// are counted.
//
// Temperature. With temperature_k_per_internal_energy_code > 0 the cache is
// set to T = factor * u (the runtime-wide convention, see
// temperatureKelvinPerInternalEnergyCode). With factor 0 (legacy/focused
// callers) it is rescaled by u_new / u_before; cells without a positive
// previous u or temperature then keep their value and are counted in
// temperature_unscaled_cells.
//
// Units. density_code is the comoving mass density and u a specific energy
// in code units; P and c_s are produced in the same comoving code frame the
// hydro solver uses, and the closure receives the evaluation-epoch scale
// factor so it can form physical densities.

// T[K] per code specific internal energy for fully ionized primordial gas
// (mu = 0.59), T = (gamma-1) mu m_p u / k_B with u in (code velocity)^2.
// Matches CoolingModelConfig::temperature_per_internal_energy_k.
[[nodiscard]] double temperatureKelvinPerInternalEnergyCode(
    double adiabatic_index, double velocity_si_per_code);

struct GasThermodynamicRefreshPolicy {
  double adiabatic_index = 5.0 / 3.0;
  double density_floor_code = 0.0;
  double pressure_floor_code = 0.0;
  // Epoch at which the closure is evaluated (the epoch of the mutation).
  double scale_factor = 1.0;
  double redshift = 0.0;
  const HydroThermodynamicClosure* closure = nullptr;
  double temperature_k_per_internal_energy_code = 0.0;
};

struct GasThermodynamicTouch {
  std::uint32_t cell_row = 0;
  // Specific internal energy before the FIRST mutation of this cell in the
  // stage; used only for the temperature-cache rescaling.
  double internal_energy_before_code = 0.0;
};

struct GasThermodynamicRefreshReport {
  std::uint64_t refreshed_cells = 0;
  std::uint64_t effective_closure_cells = 0;
  std::uint64_t pressure_floor_cells = 0;
  std::uint64_t temperature_unscaled_cells = 0;
  // max over refreshed cells of c_new / c_old (c_old > 0), the factor by
  // which the local hydro CFL bound tightens.
  double max_sound_speed_increase_ratio = 1.0;
  double max_sound_speed_code = 0.0;
};

// Refresh the derived thermodynamic caches of the touched cells and mark the
// state's gas thermodynamic caches fresh. Touches may arrive in any order;
// duplicate rows must have been collapsed by the caller keeping the first
// internal_energy_before_code (see FeedbackThermodynamicTouchLog). Rows out
// of range throw std::out_of_range before any cache is written.
GasThermodynamicRefreshReport refreshGasThermodynamics(
    core::SimulationState& state,
    std::span<const GasThermodynamicTouch> touched_cells,
    const GasThermodynamicRefreshPolicy& policy);

// Refresh every gas cell (used after bulk mutations such as regrid). The
// temperature cache is left untouched because no per-cell previous energy is
// available; callers that need temperature consistency use the touched form.
GasThermodynamicRefreshReport refreshAllGasThermodynamics(
    core::SimulationState& state,
    const GasThermodynamicRefreshPolicy& policy);

}  // namespace cosmosim::hydro
