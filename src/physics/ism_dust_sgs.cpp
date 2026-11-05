#include "cosmosim/physics/ism_dust_sgs.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {
namespace {
constexpr double k_newton_g_cgs = 6.67430e-8;
constexpr double k_pi = 3.14159265358979323846;

[[nodiscard]] double harmonicMeanPositive(double left, double right) noexcept {
  return left > 0.0 && right > 0.0 ? 2.0 * left * right / (left + right) : 0.0;
}

struct SgsFaceTransfer {
  std::array<double, 3> momentum_to_lower{};
  double energy_to_lower = 0.0;
  double stable_dt = std::numeric_limits<double>::infinity();
  bool active = false;
};

// Traceless viscous traction on a face normal to `axis`, oriented from the
// lower-coordinate cell to the upper-coordinate cell. Only the normal
// derivative is reconstructed at this face; the complete cell-centred
// contraction used for production is supplied separately by the AMR topology.
[[nodiscard]] SgsFaceTransfer sgsFaceTransfer(
    const std::array<double, 3>& lower_velocity,
    const std::array<double, 3>& upper_velocity,
    double lower_density,
    double upper_density,
    std::uint8_t axis,
    double area,
    double distance,
    double dt,
    double coefficient) noexcept {
  SgsFaceTransfer result;
  if (axis >= 3U || !(area > 0.0) || !(distance > 0.0) || !(dt > 0.0) ||
      !(coefficient > 0.0)) {
    return result;
  }
  const double density = harmonicMeanPositive(lower_density, upper_density);
  if (!(density > 0.0)) return result;

  double strain[3][3]{};
  for (std::size_t component = 0U; component < 3U; ++component) {
    const double gradient =
        (upper_velocity[component] - lower_velocity[component]) / distance;
    strain[component][axis] += 0.5 * gradient;
    strain[axis][component] += 0.5 * gradient;
  }
  const double trace_third =
      (strain[0][0] + strain[1][1] + strain[2][2]) / 3.0;
  double contraction = 0.0;
  for (std::size_t i = 0U; i < 3U; ++i) {
    strain[i][i] -= trace_third;
    for (std::size_t j = 0U; j < 3U; ++j) {
      contraction += strain[i][j] * strain[i][j];
    }
  }
  const double magnitude = std::sqrt(std::max(2.0 * contraction, 0.0));
  const double viscosity = coefficient * distance * distance * magnitude;
  if (!(viscosity > 0.0) || !std::isfinite(viscosity)) return result;
  result.stable_dt = 0.4 / (6.0 * coefficient * magnitude);

  const std::array<double, 3> face_velocity{
      0.5 * (lower_velocity[0] + upper_velocity[0]),
      0.5 * (lower_velocity[1] + upper_velocity[1]),
      0.5 * (lower_velocity[2] + upper_velocity[2])};
  const double scale = area * dt;
  for (std::size_t component = 0U; component < 3U; ++component) {
    const double traction = 2.0 * density * viscosity * strain[component][axis];
    result.momentum_to_lower[component] = scale * traction;
    result.energy_to_lower += scale * face_velocity[component] * traction;
  }
  result.active = std::isfinite(result.energy_to_lower) &&
      std::all_of(result.momentum_to_lower.begin(), result.momentum_to_lower.end(),
                  [](double value) { return std::isfinite(value); });
  return result;
}
}  // namespace

// --------------------------------------------------------------------------
// ISM regimes
// --------------------------------------------------------------------------
double jeansLengthCm(double temperature_k, double mean_molecular_weight, double rho_cgs, double adiabatic_index,
                     double extra_velocity_dispersion_sq_cm2_s2) {
  if (!(rho_cgs > 0.0) || !(mean_molecular_weight > 0.0)) return std::numeric_limits<double>::infinity();
  const double cs2 = adiabatic_index * k_boltzmann_erg_per_k * std::max(temperature_k, 0.0) /
          (mean_molecular_weight * k_amu_g) + std::max(extra_velocity_dispersion_sq_cm2_s2, 0.0);
  return std::sqrt(k_pi * cs2 / (k_newton_g_cgs * rho_cgs));
}

IsmRegimeDecision decideIsmRegime(const core::SimulationState& state, const IsmRegimeConfig& config,
                                  const ThermochemistryUnitBridge& units, const SourceInterval& interval,
                                  std::uint32_t row) {
  const auto& sidecar = state.gas_composition;
  IsmRegimeDecision decision;
  decision.tag = sidecar.ism_regime_tag[row];
  decision.dwell_years = sidecar.ism_regime_dwell_years[row] + interval.physical_seconds / k_seconds_per_year;
  const double a = interval.evaluation_scale_factor;
  const double frame_density = units.comoving ? 1.0 / (a * a * a) : 1.0;
  const double frame_length = units.comoving ? a : 1.0;
  const double h = config.hysteresis_fraction;
  const auto& descriptor = sidecar.descriptor();
  const double rho_code = state.gas_cells.density_code[row];
  const double rho_cgs = rho_code * frame_density * units.density_code_to_cgs;
  double mass = 0.0;
  for (std::uint32_t l = 0; l < descriptor.element_lane_count; ++l) mass += sidecar.value(descriptor.element_lane_begin + l, row);
  for (std::uint32_t l = 0; l < descriptor.dust_lane_count; ++l) mass += sidecar.value(descriptor.dust_lane_begin + l, row);
  const double x_h = mass > 0.0 ? sidecar.value(descriptor.element_lane_begin, row) / mass : 0.75;
  const double n_h = x_h * rho_cgs / (chemicalElementInfo(ChemicalElement::kH).atomic_mass_amu * k_amu_g);
  double mu = sidecar.mean_molecular_weight[row];
  if (!(mu > 0.0) || !std::isfinite(mu)) mu = 1.22;  // neutral reference until first chemistry commit
  const double u_cgs = state.gas_cells.internal_energy_code[row] * units.velocity_code_to_cm_s * units.velocity_code_to_cm_s;
  const double t_k = (config.adiabatic_index - 1.0) * mu * k_amu_g * u_cgs / k_boltzmann_erg_per_k;
  const double volume_code = state.cells.mass_code[row] / std::max(rho_code, 1.0e-300);
  const double dx_cm = std::cbrt(std::max(volume_code, 0.0)) * frame_length * units.length_code_to_cm;
  const double lambda = jeansLengthCm(t_k, mu, rho_cgs, config.adiabatic_index);
  const double required = config.jeans_resolution_cells * dx_cm;
  std::uint8_t desired = decision.tag;
  if (decision.tag == core::k_gas_ism_regime_effective) {
    if (lambda > required * (1.0 + h) || n_h < config.effective_threshold_n_h_cm3 * (1.0 - h)) {
      desired = core::k_gas_ism_regime_resolved;
    }
  } else if (lambda < required * (1.0 - h) && n_h >= config.effective_threshold_n_h_cm3 * (1.0 + h)) {
    desired = core::k_gas_ism_regime_effective;
  }
  if (desired != decision.tag) {
    if (decision.dwell_years < config.minimum_dwell_years) {
      decision.held_by_dwell = true;
    } else {
      decision.tag = desired;
      decision.dwell_years = 0.0;
      decision.changed = true;
    }
  }
  return decision;
}

IsmRegimeReport updateIsmRegimes(core::SimulationState& state, const IsmRegimeConfig& config,
                                 const ThermochemistryUnitBridge& units, const SourceInterval& interval,
                                 std::span<const std::uint32_t> rows) {
  IsmRegimeReport report;
  auto& sidecar = state.gas_composition;
  if (!sidecar.enabled()) {
    throw std::runtime_error("ISM regime manager requires the gas composition sidecar");
  }
  sidecar.requireAlignedWith(state.gas_cells.gas_cell_id, "updateIsmRegimes");
  const std::size_t count = rows.empty() ? state.cells.size() : rows.size();
  // Decide for all rows from the unmodified state first, then apply.
  std::vector<std::pair<std::uint32_t, IsmRegimeDecision>> decisions;
  decisions.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    const std::uint32_t row = rows.empty() ? static_cast<std::uint32_t>(i) : rows[i];
    decisions.emplace_back(row, decideIsmRegime(state, config, units, interval, row));
  }
  for (const auto& [row, d] : decisions) {
    sidecar.ism_regime_tag[row] = d.tag;
    sidecar.ism_regime_dwell_years[row] = d.dwell_years;
    if (d.held_by_dwell) ++report.held_by_dwell;
    if (d.changed) {
      report.changed_rows.push_back(row);
      if (d.tag == core::k_gas_ism_regime_effective) ++report.to_effective; else ++report.to_resolved;
    }
    if (d.tag == core::k_gas_ism_regime_effective) ++report.effective_cells; else ++report.resolved_cells;
  }
  return report;
}

// --------------------------------------------------------------------------
// Dust
// --------------------------------------------------------------------------
namespace {
// Largest dust mass of reservoir d that the gas element lanes can still
// supply, limited by the key element of the stoichiometry.
double condensableDustMass(const core::GasCompositionSidecar& sidecar, const CompositionLayout& layout,
                           std::size_t row, std::size_t d) {
  const auto& fractions = dustElementMassFractions(static_cast<DustReservoir>(d));
  const bool detailed = sidecar.descriptor().species_lane_count >= k_detailed_species_count;
  double co_mass = 0.0, cii_mass = 0.0;
  if (detailed) {
    co_mass = sidecar.value(layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kCO)), row);
    cii_mass = sidecar.value(layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kCII)), row);
  }
  const double a_c = chemicalElementInfo(ChemicalElement::kC).atomic_mass_amu;
  const double a_o = chemicalElementInfo(ChemicalElement::kO).atomic_mass_amu;
  double limit = std::numeric_limits<double>::infinity();
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    if (fractions[e] <= 0.0) continue;
    double available = sidecar.value(layout.elementLane(static_cast<ChemicalElement>(e)), row);
    // Molecule-bound nuclei cannot condense: CO keeps its C and O.
    if (detailed && e == static_cast<std::size_t>(ChemicalElement::kC)) available = cii_mass;
    if (detailed && e == static_cast<std::size_t>(ChemicalElement::kO)) available -= co_mass * a_o / (a_c + a_o);
    limit = std::min(limit, std::max(available, 0.0) / fractions[e]);
  }
  return std::isfinite(limit) ? std::max(limit, 0.0) : 0.0;
}

void transferGasToDust(core::GasCompositionSidecar& sidecar, const CompositionLayout& layout, std::size_t row,
                       std::size_t d, double dust_mass_delta) {
  // dust_mass_delta > 0: gas -> dust; < 0: dust -> gas. Stoichiometric.
  const bool detailed_row_was_explicit =
      !sidecar.sparseStorageEnabled() || sidecar.hasExplicitSparseRow(row);
  const auto& fractions = dustElementMassFractions(static_cast<DustReservoir>(d));
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    if (fractions[e] <= 0.0) continue;
    const std::uint32_t lane = layout.elementLane(static_cast<ChemicalElement>(e));
    sidecar.setValue(lane, row, std::max(sidecar.value(lane, row) - fractions[e] * dust_mass_delta, 0.0));
  }
  const std::uint32_t dust_lane = layout.dustLane(d);
  sidecar.setValue(dust_lane, row, std::max(sidecar.value(dust_lane, row) + dust_mass_delta, 0.0));
  // Species lanes of the networks (H/He/C/O) are untouched by refractory
  // transfers except carbon: C+ (and CO) track gas-phase carbon in the
  // detailed network, so the carbon removed/returned is applied to C+.
  if (fractions[static_cast<std::size_t>(ChemicalElement::kC)] > 0.0 &&
      sidecar.descriptor().species_lane_count >= k_detailed_species_count &&
      detailed_row_was_explicit) {
    const std::uint32_t cii = layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kCII));
    const double delta_c = fractions[static_cast<std::size_t>(ChemicalElement::kC)] * dust_mass_delta;
    sidecar.setValue(cii, row, std::max(sidecar.value(cii, row) - delta_c, 0.0));
  }
}
}  // namespace

void advanceDustCell(core::GasCompositionSidecar& sidecar, const CompositionLayout& layout, std::size_t row,
                     double n_h_cm3, double temperature_k, double dt_seconds, const DustConfig& config,
                     DustUpdateLedger& ledger) {
  if (sidecar.descriptor().dust_lane_count == 0U || !(dt_seconds > 0.0)) return;
  double metals = 0.0, mass = 0.0;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    const double m = sidecar.value(layout.elementLane(static_cast<ChemicalElement>(e)), row);
    mass += m;
    if (e >= 2) metals += m;
  }
  for (std::size_t d = 0; d < sidecar.descriptor().dust_lane_count; ++d) {
    const double dust = sidecar.value(layout.dustLane(d), row);
    mass += dust;
    metals += dust;
  }
  if (!(mass > 0.0)) return;
  const double z_rel = (metals / mass) / solarMetalMassFraction();
  const double t = std::max(temperature_k, 1.0);
  const double n_h = std::max(n_h_cm3, 0.0);
  const double tau_sp = config.sputtering_timescale_ref_s / std::max(n_h / 1.0e-3, 1.0e-300) /
      (std::pow(2.0e6 / t, 2.5) + 1.0);
  const double tau_acc = (z_rel > 0.0 && n_h > 0.0)
      ? config.accretion_timescale_ref_s / ((n_h / 100.0) * std::sqrt(t / 50.0) * z_rel)
      : std::numeric_limits<double>::infinity();
  for (std::size_t d = 0; d < sidecar.descriptor().dust_lane_count; ++d) {
    const double d0 = sidecar.value(layout.dustLane(d), row);
    if (!(d0 > 0.0)) continue;  // growth needs seed grains; sputtering needs dust
    // Operator split within the cell: sputtering (exact) then accretion (exact).
    const double after_sputter = d0 * std::exp(-3.0 * dt_seconds / tau_sp);
    transferGasToDust(sidecar, layout, row, d, after_sputter - d0);
    ledger.dust_to_gas_code += d0 - after_sputter;
    if (std::isfinite(tau_acc) && after_sputter > 0.0) {
      const double d_max = after_sputter + condensableDustMass(sidecar, layout, row, d);
      const double grown = d_max / (1.0 + (d_max / after_sputter - 1.0) * std::exp(-dt_seconds / tau_acc));
      const double delta = std::clamp(grown - after_sputter, 0.0, d_max - after_sputter);
      transferGasToDust(sidecar, layout, row, d, delta);
      ledger.gas_to_dust_code += delta;
    }
  }
  ++ledger.cells;
}

double destroyDustBySupernova(core::GasCompositionSidecar& sidecar, const CompositionLayout& layout,
                              std::size_t row, double swept_gas_mass_code, double cell_gas_mass_code) {
  if (sidecar.descriptor().dust_lane_count == 0U || !(cell_gas_mass_code > 0.0)) return 0.0;
  const double fraction = std::clamp(swept_gas_mass_code / cell_gas_mass_code, 0.0, 1.0);
  double destroyed_total = 0.0;
  for (std::size_t d = 0; d < sidecar.descriptor().dust_lane_count; ++d) {
    const double destroyed = fraction * sidecar.value(layout.dustLane(d), row);
    transferGasToDust(sidecar, layout, row, d, -destroyed);
    destroyed_total += destroyed;
  }
  return destroyed_total;
}

EjectaPartition partitionEjectaIntoDust(const std::array<double, k_chemical_element_count>& ejecta,
                                        double condensation_efficiency) {
  EjectaPartition result;
  result.gas_elements = ejecta;
  const double delta = std::clamp(condensation_efficiency, 0.0, 1.0);
  if (delta == 0.0) return result;
  for (std::size_t d = 0; d < k_dust_reservoir_count; ++d) {
    const auto& fractions = dustElementMassFractions(static_cast<DustReservoir>(d));
    double limit = std::numeric_limits<double>::infinity();
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
      if (fractions[e] > 0.0) limit = std::min(limit, result.gas_elements[e] / fractions[e]);
    }
    const double dust = std::isfinite(limit) ? delta * std::max(limit, 0.0) : 0.0;
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
      result.gas_elements[e] = std::max(result.gas_elements[e] - fractions[e] * dust, 0.0);
    }
    result.dust[d] = dust;
  }
  return result;
}

// --------------------------------------------------------------------------
// SGS
// --------------------------------------------------------------------------
std::uint64_t SgsShearTransferWorkspace::ownedCapacityBytes() const noexcept {
  const std::uint64_t doubles = static_cast<std::uint64_t>(
      velocity_x.capacity() + velocity_y.capacity() + velocity_z.capacity() +
      internal_energy.capacity() + sgs_energy.capacity());
  return doubles * sizeof(double) +
      static_cast<std::uint64_t>(touched.capacity()) * sizeof(std::uint8_t);
}

std::uint64_t SgsShearTransferWorkspace::requiredBytes(std::size_t rows) noexcept {
  constexpr std::uint64_t width = 5U * sizeof(double) + sizeof(std::uint8_t);
  if (rows > std::numeric_limits<std::uint64_t>::max() / width) {
    return std::numeric_limits<std::uint64_t>::max();
  }
  return static_cast<std::uint64_t>(rows) * width;
}

void SgsShearTransferWorkspace::clear() noexcept {
  velocity_x.clear();
  velocity_y.clear();
  velocity_z.clear();
  internal_energy.clear();
  sgs_energy.clear();
  touched.clear();
  prepared = false;
}

SgsLedger prepareSgsShearTransfer(
    const core::SimulationState& state,
    std::span<const std::uint8_t> active_rows,
    std::span<const double> trace_free_strain_magnitude_code,
    std::span<const MetalDiffusionFace> local_faces,
    std::span<const SgsRemoteShearFace> remote_faces,
    double physical_length_per_stored_length,
    double dt_code_time,
    const SgsConfig& config,
    SgsShearTransferWorkspace& workspace) {
  workspace.clear();
  const std::size_t rows = state.cells.size();
  if (active_rows.size() != rows || trace_free_strain_magnitude_code.size() != rows ||
      state.gas_cells.size() != rows || !(physical_length_per_stored_length > 0.0) ||
      !std::isfinite(physical_length_per_stored_length) || !(dt_code_time > 0.0) ||
      !std::isfinite(dt_code_time) || !(config.eddy_viscosity_coefficient > 0.0) ||
      !std::isfinite(config.eddy_viscosity_coefficient)) {
    throw std::invalid_argument("SGS shear transfer received inconsistent state, interval or closure coefficient");
  }
  const auto& sidecar = state.gas_composition;
  if (!sidecar.enabled() ||
      sidecar.descriptor().sgs_energy_lane == core::GasCompositionLayoutDescriptor::k_no_lane) {
    throw std::invalid_argument("SGS shear transfer requires the authoritative SGS energy lane");
  }
  sidecar.requireAlignedWith(state.gas_cells.gas_cell_id, "prepareSgsShearTransfer");
  const std::uint32_t sgs_lane = sidecar.descriptor().sgs_energy_lane;

  workspace.velocity_x.assign(rows, 0.0);  // staged momentum delta x
  workspace.velocity_y.assign(rows, 0.0);  // staged momentum delta y
  workspace.velocity_z.assign(rows, 0.0);  // staged momentum delta z
  workspace.internal_energy.assign(rows, 0.0);  // staged resolved-total energy delta
  workspace.sgs_energy.assign(rows, 0.0);  // staged SGS production
  workspace.touched.assign(rows, 0U);
  SgsLedger ledger;

  const auto velocity = [&state](std::uint32_t row) {
    return std::array<double, 3>{
        state.gas_cells.velocity_x_peculiar[row],
        state.gas_cells.velocity_y_peculiar[row],
        state.gas_cells.velocity_z_peculiar[row]};
  };
  const double volume_factor = physical_length_per_stored_length *
      physical_length_per_stored_length * physical_length_per_stored_length;
  const auto physicalDensity = [&state, volume_factor](std::uint32_t row) {
    const double stored_density = state.gas_cells.density_code[row];
    return stored_density > 0.0 ? stored_density / volume_factor : 0.0;
  };
  const auto accumulate = [&](std::uint32_t row, const SgsFaceTransfer& transfer,
                              double sign) {
    if (row >= rows || active_rows[row] == 0U || !transfer.active) return;
    workspace.velocity_x[row] += sign * transfer.momentum_to_lower[0];
    workspace.velocity_y[row] += sign * transfer.momentum_to_lower[1];
    workspace.velocity_z[row] += sign * transfer.momentum_to_lower[2];
    workspace.internal_energy[row] += sign * transfer.energy_to_lower;
    workspace.touched[row] = 1U;
    ledger.resolved_stress_work_code += sign * transfer.energy_to_lower;
    ledger.momentum_impulse_l1_code +=
        std::abs(transfer.momentum_to_lower[0]) +
        std::abs(transfer.momentum_to_lower[1]) +
        std::abs(transfer.momentum_to_lower[2]);
  };

  for (const MetalDiffusionFace& face : local_faces) {
    if (face.left_cell >= rows || face.right_cell >= rows || face.axis >= 3U ||
        active_rows[face.left_cell] == 0U || active_rows[face.right_cell] == 0U) {
      continue;
    }
    const SgsFaceTransfer transfer = sgsFaceTransfer(
        velocity(face.left_cell), velocity(face.right_cell),
        physicalDensity(face.left_cell), physicalDensity(face.right_cell),
        face.axis, face.area_code, face.center_distance_code, dt_code_time,
        config.eddy_viscosity_coefficient);
    if (!transfer.active) continue;
    if (dt_code_time > transfer.stable_dt *
            (1.0 + 16.0 * std::numeric_limits<double>::epsilon())) {
      throw core::StageNumericalRetryRequest(
          "SGS face stress exceeds its explicit viscous timestep (requested " +
          std::to_string(dt_code_time) + ", stable " +
          std::to_string(transfer.stable_dt) + ")");
    }
    accumulate(face.left_cell, transfer, +1.0);
    accumulate(face.right_cell, transfer, -1.0);
    ++ledger.shear_faces;
  }
  for (const SgsRemoteShearFace& face : remote_faces) {
    if (face.local_row >= rows || face.axis >= 3U ||
        active_rows[face.local_row] == 0U ||
        (face.local_side != -1 && face.local_side != 1)) {
      continue;
    }
    const std::array<double, 3> local_velocity = velocity(face.local_row);
    const bool local_is_lower = face.local_side > 0;
    const SgsFaceTransfer transfer = sgsFaceTransfer(
        local_is_lower ? local_velocity : face.remote_velocity_peculiar,
        local_is_lower ? face.remote_velocity_peculiar : local_velocity,
        local_is_lower ? physicalDensity(face.local_row) : face.remote_density_physical_code,
        local_is_lower ? face.remote_density_physical_code : physicalDensity(face.local_row),
        face.axis, face.area_code, face.center_distance_code, dt_code_time,
        config.eddy_viscosity_coefficient);
    if (!transfer.active) continue;
    if (dt_code_time > transfer.stable_dt *
            (1.0 + 16.0 * std::numeric_limits<double>::epsilon())) {
      throw core::StageNumericalRetryRequest(
          "SGS remote-face stress exceeds its explicit viscous timestep (requested " +
          std::to_string(dt_code_time) + ", stable " +
          std::to_string(transfer.stable_dt) + ")");
    }
    accumulate(face.local_row, transfer, local_is_lower ? +1.0 : -1.0);
    ++ledger.shear_faces;
  }

  for (std::uint32_t row = 0U; row < rows; ++row) {
    if (active_rows[row] == 0U) continue;
    const double mass = state.cells.mass_code[row];
    const double density = state.gas_cells.density_code[row];
    const double strain = trace_free_strain_magnitude_code[row];
    if (!(mass > 0.0) || !(density > 0.0) || !(strain >= 0.0) ||
        !std::isfinite(strain)) {
      if (!std::isfinite(strain)) {
        throw std::runtime_error("SGS shear transfer encountered a non-finite AMR strain estimate");
      }
      continue;
    }
    const double volume_physical = (mass / density) * volume_factor;
    const double filter_length = std::cbrt(std::max(volume_physical, 0.0));
    const double viscosity = config.eddy_viscosity_coefficient *
        filter_length * filter_length * strain;
    // Three-dimensional explicit viscous bound. The 0.4 safety factor is the
    // same parabolic margin used by the repository's explicit diffusion
    // reference; it is a numerical stability limit, not a fitted SGS rate.
    constexpr double k_viscous_cfl = 0.4;
    const double stable_dt = strain > 0.0
        ? k_viscous_cfl /
            (6.0 * config.eddy_viscosity_coefficient * strain)
        : std::numeric_limits<double>::infinity();
    if (dt_code_time > stable_dt *
            (1.0 + 16.0 * std::numeric_limits<double>::epsilon())) {
      throw core::StageNumericalRetryRequest(
          "SGS shear transfer exceeds its explicit viscous timestep (requested " +
          std::to_string(dt_code_time) + ", stable " +
          std::to_string(stable_dt) + ")");
    }
    const double produced = mass * viscosity * strain * strain * dt_code_time;
    if (!(produced >= 0.0) || !std::isfinite(produced)) {
      throw std::runtime_error("SGS shear production is non-finite");
    }
    if (produced > 0.0) {
      workspace.internal_energy[row] -= produced;
      workspace.sgs_energy[row] += produced;
      workspace.touched[row] = 1U;
      ledger.produced_from_resolved_code += produced;
      ++ledger.shear_cells;
    }
  }

  for (std::uint32_t row = 0U; row < rows; ++row) {
    if (workspace.touched[row] == 0U) continue;
    const double mass = state.cells.mass_code[row];
    if (!(mass > 0.0)) {
      throw std::runtime_error("SGS shear transfer touched a non-positive-mass gas cell");
    }
    const std::array<double, 3> old_velocity = velocity(row);
    const std::array<double, 3> new_momentum{
        mass * old_velocity[0] + workspace.velocity_x[row],
        mass * old_velocity[1] + workspace.velocity_y[row],
        mass * old_velocity[2] + workspace.velocity_z[row]};
    const double old_kinetic = 0.5 * mass *
        (old_velocity[0] * old_velocity[0] +
         old_velocity[1] * old_velocity[1] +
         old_velocity[2] * old_velocity[2]);
    const double new_kinetic = 0.5 *
        (new_momentum[0] * new_momentum[0] +
         new_momentum[1] * new_momentum[1] +
         new_momentum[2] * new_momentum[2]) / mass;
    const double new_resolved_total = mass * state.gas_cells.internal_energy_code[row] +
        old_kinetic + workspace.internal_energy[row];
    const double new_internal = (new_resolved_total - new_kinetic) / mass;
    const double new_sgs = sidecar.value(sgs_lane, row) + workspace.sgs_energy[row];
    if (!(new_internal > 0.0) || !std::isfinite(new_internal) ||
        !(new_sgs >= 0.0) || !std::isfinite(new_sgs)) {
      throw core::StageNumericalRetryRequest(
          "SGS shear transfer would violate positive resolved/SGS energy; "
          "retry with a shorter interval");
    }
    workspace.velocity_x[row] = new_momentum[0] / mass;
    workspace.velocity_y[row] = new_momentum[1] / mass;
    workspace.velocity_z[row] = new_momentum[2] / mass;
    workspace.internal_energy[row] = new_internal;
    workspace.sgs_energy[row] = new_sgs;
  }
  workspace.prepared = true;
  return ledger;
}

void commitSgsShearTransfer(
    core::SimulationState& state,
    SgsShearTransferWorkspace& workspace) noexcept {
  if (!workspace.prepared) return;
  const std::uint32_t lane = state.gas_composition.descriptor().sgs_energy_lane;
  for (std::uint32_t row = 0U; row < workspace.touched.size(); ++row) {
    if (workspace.touched[row] == 0U) continue;
    state.gas_cells.velocity_x_peculiar[row] = workspace.velocity_x[row];
    state.gas_cells.velocity_y_peculiar[row] = workspace.velocity_y[row];
    state.gas_cells.velocity_z_peculiar[row] = workspace.velocity_z[row];
    state.gas_cells.internal_energy_code[row] = workspace.internal_energy[row];
    state.gas_composition.setValue(lane, row, workspace.sgs_energy[row]);
  }
  workspace.prepared = false;
}

double sgsSpecificEnergyAfter(double k0, double dt, double dx, double c_eps) noexcept {
  if (!(k0 > 0.0) || !(dt > 0.0) || !(dx > 0.0)) return std::max(k0, 0.0);
  const double denom = 1.0 + c_eps * std::sqrt(k0) * dt / (2.0 * dx);
  return k0 / (denom * denom);
}

SgsDissipation sgsDissipate(double sgs_energy, double mass, double dt_code_time, double cell_size_code,
                            const SgsConfig& config) noexcept {
  SgsDissipation out;
  out.sgs_energy = std::max(sgs_energy, 0.0);
  if (!(mass > 0.0)) return out;
  const double k0 = out.sgs_energy / mass;
  const double k1 = sgsSpecificEnergyAfter(k0, dt_code_time, cell_size_code, config.dissipation_coefficient);
  out.sgs_energy = mass * k1;
  out.dissipated = mass * (k0 - k1);
  return out;
}

double advanceSgsCell(core::SimulationState& state, const CompositionLayout& layout, std::size_t row,
                      double dt_code_time, double cell_size_code, const SgsConfig& config, SgsLedger& ledger) {
  auto& sidecar = state.gas_composition;
  if (!layout.hasSgs()) return 0.0;
  const std::uint32_t lane = sidecar.descriptor().sgs_energy_lane;
  const double mass = state.cells.mass_code[row];
  if (!(mass > 0.0)) return 0.0;
  const double e_sgs = sidecar.value(lane, row);
  const double k0 = e_sgs / mass;
  const double k1 = sgsSpecificEnergyAfter(k0, dt_code_time, cell_size_code, config.dissipation_coefficient);
  const double dissipated = mass * (k0 - k1);
  sidecar.setValue(lane, row, mass * k1);
  state.gas_cells.internal_energy_code[row] += dissipated / mass;  // heat, exact transfer
  ledger.dissipated_to_heat_code += dissipated;
  return dissipated;
}

}  // namespace cosmosim::physics
