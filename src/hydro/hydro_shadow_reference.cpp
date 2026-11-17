#include "cosmosim/hydro/hydro_shadow_reference.hpp"

#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace cosmosim::hydro {
namespace {

[[nodiscard]] std::uint64_t mix64(std::uint64_t x) noexcept {
  x ^= x >> 30U;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27U;
  x *= 0x94d049bb133111ebULL;
  return x ^ (x >> 31U);
}

[[nodiscard]] std::uint64_t candidateScore(
    const HydroShadowCandidate& candidate, std::uint64_t step, std::uint64_t seed) noexcept {
  return mix64(mix64(hydroShadowEffectiveKey(candidate)) ^
               mix64(static_cast<std::uint32_t>(candidate.level)) ^ mix64(step) ^ mix64(seed));
}

[[nodiscard]] bool sameCandidate(
    const HydroShadowCandidate& a, const HydroShadowCandidate& b) noexcept {
  return hydroShadowEffectiveKey(a) == hydroShadowEffectiveKey(b) && a.level == b.level;
}

[[nodiscard]] std::uint64_t checkedMul(std::uint64_t a, std::uint64_t b) {
  if (b != 0U && a > std::numeric_limits<std::uint64_t>::max() / b) {
    throw std::overflow_error("hydro shadow byte multiplication overflow");
  }
  return a * b;
}

[[nodiscard]] std::uint64_t checkedAdd(std::uint64_t a, std::uint64_t b) {
  if (b > std::numeric_limits<std::uint64_t>::max() - a) {
    throw std::overflow_error("hydro shadow byte addition overflow");
  }
  return a + b;
}

[[nodiscard]] std::uint64_t vectorBytes(std::size_t count, std::size_t width) {
  return checkedMul(static_cast<std::uint64_t>(count), static_cast<std::uint64_t>(width));
}

[[nodiscard]] double kineticDensity(const HydroConservedState& state) noexcept {
  const double rho = std::max(state.mass_density_comoving, 1.0e-30);
  const double momentum_sq = state.momentum_density_x_comoving * state.momentum_density_x_comoving +
      state.momentum_density_y_comoving * state.momentum_density_y_comoving +
      state.momentum_density_z_comoving * state.momentum_density_z_comoving;
  return 0.5 * momentum_sq / rho;
}

[[nodiscard]] double specificThermalEnergy(const HydroConservedState& state) noexcept {
  const double rho = std::max(state.mass_density_comoving, 1.0e-30);
  return std::max((state.total_energy_density_comoving - kineticDensity(state)) / rho, 0.0);
}

[[nodiscard]] double selectedThermalEnergy(
    const HydroConservedStateSoa& state, std::size_t row) {
  const auto conserved = state.loadCell(row);
  if (state.auxiliaryInternalEnergyEnabled()) {
    const auto auxiliary = state.auxiliaryInternalEnergyDensityComoving()[row];
    if (std::isfinite(auxiliary) && auxiliary > 0.0 && conserved.mass_density_comoving > 0.0) {
      return auxiliary / conserved.mass_density_comoving;
    }
  }
  return specificThermalEnergy(conserved);
}

[[nodiscard]] double compositionGasMassOf(const core::GasCompositionSidecar& sidecar, std::size_t row) {
  const auto& d = sidecar.descriptor();
  double mass = 0.0;
  for (std::uint32_t l = 0; l < d.element_lane_count; ++l) mass += sidecar.value(d.element_lane_begin + l, row);
  for (std::uint32_t l = 0; l < d.dust_lane_count; ++l) mass += sidecar.value(d.dust_lane_begin + l, row);
  return mass;
}

// Composition-dependent temperature factor, identical to the production
// thermodynamic refresh (rowTemperatureFactor): base * mu * (amu/m_p) / 0.59.
[[nodiscard]] double temperatureFactorForMu(double base, double mu_amu) noexcept {
  constexpr double k_amu_over_proton_mass = 1.66053906660e-24 / 1.67262192369e-24;
  constexpr double k_reference_mu = 0.59;
  if (base > 0.0 && std::isfinite(mu_amu) && mu_amu > 0.0) {
    return base * mu_amu * k_amu_over_proton_mass / k_reference_mu;
  }
  return base;
}

// The provider never reads SimulationState or an authoritative mutable cache.
// Density and temperature are reconstructed from each PRIVATE intermediate
// state. With a composition-aware closure (temperature factor and private mu
// available) T = factor(mu) * u, exactly the production refresh; otherwise the
// captured temperature is rescaled by u/u0 (frozen-composition fallback,
// counted by the caller's metadata as the legacy closure).
class ShadowSourceProperties final : public HydroCellSourcePropertyProvider {
 public:
  ShadowSourceProperties(const HydroConservedStateSoa& reference,
                         const HydroConservedStateSoa& initial,
                         std::span<const HydroCellSourceProperties> baseline,
                         double initial_scale, double substep_scale, bool comoving,
                         double temperature_factor, double density_cgs_per_code,
                         bool comoving_density,
                         const core::GasCompositionSidecar* composition)
      : m_reference(reference), m_initial(initial), m_baseline(baseline),
        m_initial_scale(initial_scale), m_substep_scale(substep_scale), m_comoving(comoving),
        m_temperature_factor(temperature_factor), m_density_cgs_per_code(density_cgs_per_code),
        m_comoving_density(comoving_density), m_composition(composition) {}

  [[nodiscard]] HydroCellSourceProperties sourcePropertiesForCell(std::size_t cell_index) const override {
    const auto initial = m_initial.loadCell(cell_index);
    const auto current = m_reference.loadCell(cell_index);
    const auto baseline = m_baseline.at(cell_index);
    HydroCellSourceProperties result = baseline;
    // density_cgs_per_code is expressed at a = 1 when the density frame is
    // comoving: the physical density is recomputed at this substep's epoch.
    const double frame = m_comoving_density && m_substep_scale > 0.0
        ? 1.0 / (m_substep_scale * m_substep_scale * m_substep_scale) : 1.0;
    if (m_density_cgs_per_code > 0.0) {
      const double scale = current.mass_density_comoving * frame * m_density_cgs_per_code /
          std::max(baseline.mass_density_physical_cgs, 1.0e-300);
      result.mass_density_physical_cgs = current.mass_density_comoving * frame * m_density_cgs_per_code;
      result.hydrogen_number_density_cgs *= baseline.mass_density_physical_cgs > 0.0 ? scale : 0.0;
    } else {
      const double ratio = current.mass_density_comoving /
          std::max(initial.mass_density_comoving, 1.0e-30);
      const double dilution = m_comoving
          ? std::pow(m_initial_scale / m_substep_scale, 3.0) : 1.0;
      result.mass_density_physical_cgs *= ratio * dilution;
      result.hydrogen_number_density_cgs *= ratio * dilution;
    }
    const double u_current = selectedThermalEnergy(m_reference, cell_index);
    double mu = std::numeric_limits<double>::quiet_NaN();
    if (m_composition != nullptr && m_composition->enabled() && cell_index < m_composition->rowCount()) {
      mu = m_composition->mean_molecular_weight[cell_index];
    }
    if (m_temperature_factor > 0.0 && std::isfinite(mu) && mu > 0.0) {
      result.temperature_k = temperatureFactorForMu(m_temperature_factor, mu) * std::max(u_current, 0.0);
    } else {
      const double u_initial = selectedThermalEnergy(m_initial, cell_index);
      result.temperature_k *= u_current / std::max(u_initial, 1.0e-30);
    }
    result.metallicity_mass_fraction = current.mass_density_comoving > 0.0
        ? std::clamp(current.metal_mass_density_comoving / current.mass_density_comoving, 0.0, 1.0)
        : baseline.metallicity_mass_fraction;
    return result;
  }

 private:
  const HydroConservedStateSoa& m_reference;
  const HydroConservedStateSoa& m_initial;
  std::span<const HydroCellSourceProperties> m_baseline;
  double m_initial_scale;
  double m_substep_scale;
  bool m_comoving;
  double m_temperature_factor;
  double m_density_cgs_per_code;
  bool m_comoving_density;
  const core::GasCompositionSidecar* m_composition;
};

// PRIVATE composition transport: the same donor-cell rule as the production
// CompositionFluxAccumulator (species flux = F_rho * X_donor with donor
// fractions taken from the halfstep-start composition), applied only to the
// sample's private real rows. The private reference has no reflux, so every
// face (including coarse-role faces) transports composition with exactly the
// private mass flux; Sum(lanes) therefore follows the private mass.
class ShadowCompositionSink final : public HydroFaceMassFluxSink {
 public:
  ShadowCompositionSink(core::GasCompositionSidecar& composition,
                        const HydroPatchGeometry& geometry,
                        std::span<const double> ghost_donor_fractions,
                        std::vector<std::size_t>& slot_by_storage_row)
      : m_composition(composition), m_geometry(geometry),
        m_ghost_fractions(ghost_donor_fractions), m_slot_by_row(slot_by_storage_row) {}

  void recordFaceMassTransfer(std::size_t face_index, std::size_t owner_cell,
                              std::size_t neighbor_cell, double mass_transfer_code) override {
    if (mass_transfer_code == 0.0) return;
    if (!std::isfinite(mass_transfer_code)) {
      throw std::runtime_error("shadow composition: non-finite private face mass transfer");
    }
    std::size_t donor_row = owner_cell;
    if (mass_transfer_code < 0.0) {
      donor_row = neighbor_cell != k_invalid_cell_index
          ? neighbor_cell : m_geometry.faces.at(face_index).neighbor_cell;
    }
    m_records.push_back(Record{owner_cell, neighbor_cell, donor_row, mass_transfer_code});
  }

  // Applies all recorded transfers against the halfstep-start composition.
  // Returns the number of faces whose ghost donor composition was unknown
  // (the receiving cell's own fractions are used and the cell is marked).
  std::uint64_t apply(std::vector<std::uint8_t>& composition_unknown_cell) {
    const std::size_t lanes = m_composition.laneCount();
    const std::size_t rows = m_composition.rowCount();
    const std::size_t real = m_geometry.cellCount();
    std::vector<double> delta(lanes * rows, 0.0);
    std::vector<double> fraction(lanes, 0.0);
    std::uint64_t unknown_faces = 0U;
    const auto fractions_of_real = [&](std::size_t row, std::vector<double>& out) {
      const double mass = compositionGasMassOf(m_composition, row);
      if (!(mass > 0.0)) {
        throw std::runtime_error("shadow composition: private donor has no composition mass");
      }
      for (std::size_t l = 0; l < lanes; ++l) out[l] = m_composition.value(l, row) / mass;
    };
    for (const Record& record : m_records) {
      bool unknown = false;
      if (record.donor < real) {
        fractions_of_real(record.donor, fraction);
      } else {
        const std::size_t slot = record.donor < m_slot_by_row.size()
            ? m_slot_by_row[record.donor] : k_invalid_ghost_cell_slot;
        const bool have = slot != k_invalid_ghost_cell_slot &&
            (slot + 1U) * lanes <= m_ghost_fractions.size() &&
            std::isfinite(m_ghost_fractions[slot * lanes]);
        if (have) {
          for (std::size_t l = 0; l < lanes; ++l) fraction[l] = m_ghost_fractions[slot * lanes + l];
        } else {
          unknown = true;
          const std::size_t receiver = record.owner < real ? record.owner : record.neighbor;
          if (receiver >= real) continue;
          fractions_of_real(receiver, fraction);
        }
      }
      if (unknown) ++unknown_faces;
      const auto deposit = [&](std::size_t row, double sign) {
        if (row >= real) return;
        if (unknown) composition_unknown_cell[row] = 1U;
        for (std::size_t l = 0; l < lanes; ++l) delta[l * rows + row] += sign * record.mass * fraction[l];
      };
      deposit(record.owner, -1.0);
      if (record.neighbor != k_invalid_cell_index) deposit(record.neighbor, +1.0);
    }
    for (std::size_t l = 0; l < lanes; ++l) {
      for (std::size_t row = 0; row < rows; ++row) {
        const double d = delta[l * rows + row];
        if (d != 0.0) m_composition.addValue(l, row, d);
      }
    }
    m_records.clear();
    return unknown_faces;
  }

  void reserve(std::size_t faces) { m_records.reserve(faces); }

 private:
  struct Record {
    std::size_t owner = k_invalid_cell_index;
    std::size_t neighbor = k_invalid_cell_index;
    std::size_t donor = k_invalid_cell_index;
    double mass = 0.0;
  };
  core::GasCompositionSidecar& m_composition;
  const HydroPatchGeometry& m_geometry;
  std::span<const double> m_ghost_fractions;
  std::vector<std::size_t>& m_slot_by_row;
  std::vector<Record> m_records;
};

[[nodiscard]] std::array<double, 7> lanes(const HydroConservedState& state, double auxiliary) {
  return {state.mass_density_comoving, state.momentum_density_x_comoving,
          state.momentum_density_y_comoving, state.momentum_density_z_comoving,
          state.total_energy_density_comoving, state.metal_mass_density_comoving, auxiliary};
}

[[nodiscard]] std::array<double, 6> physicalLanes(const HydroConservedState& state) {
  const auto l = lanes(state, 0.0);
  return {l[0], l[1], l[2], l[3], l[4], l[5]};
}

[[nodiscard]] bool finiteState(const HydroConservedState& state) {
  for (double value : physicalLanes(state)) if (!std::isfinite(value)) return false;
  return state.mass_density_comoving > 0.0;
}

[[nodiscard]] bool selfAuthority(HydroShadowGhostProvenance provenance) noexcept {
  return provenance == HydroShadowGhostProvenance::kSelfPhysical ||
         provenance == HydroShadowGhostProvenance::kSelfPeriodicWrap;
}

[[nodiscard]] HydroShadowComparisonClass hydroClassFor(bool temporal_refined,
                                                       std::uint64_t cone_excluded) noexcept {
  if (cone_excluded != 0U) return HydroShadowComparisonClass::kBoundedProjection;
  return temporal_refined ? HydroShadowComparisonClass::kTimeOnlyCtuPpm
                          : HydroShadowComparisonClass::kMixedSpatialTemporal;
}

void addMetric(HydroShadowMetric& metric, double produced, double expected, double volume) {
  if (!std::isfinite(produced) || !std::isfinite(expected)) {
    throw std::runtime_error("non-finite lane in shadow comparison");
  }
  const double error = std::abs(produced - expected);
  metric.weighted_abs_difference += error * volume;
  metric.weighted_reference_magnitude += std::abs(expected) * volume;
  metric.linf_difference = std::max(metric.linf_difference, error);
}

// Region <-> SoA conversion. The authoritative gas state between operators
// is (rho, momentum, e, Z); E is derived as rho (e + |v|^2/2).
void recomputeRegionEnergy(HydroShadowPrivateRegion& region) {
  for (std::size_t i = 0; i < region.cell_count; ++i) {
    const double rho = std::max(region.mass_density[i], 1.0e-300);
    const double momentum_sq = region.momentum_x[i] * region.momentum_x[i] +
        region.momentum_y[i] * region.momentum_y[i] + region.momentum_z[i] * region.momentum_z[i];
    region.total_energy[i] = rho * std::max(region.specific_internal_energy[i], 0.0) + 0.5 * momentum_sq / rho;
  }
}

void regionToSoa(const HydroShadowPrivateRegion& region, HydroConservedStateSoa& soa) {
  for (std::size_t i = 0; i < region.cell_count; ++i) {
    soa.storeCell(i, HydroConservedState{
        .mass_density_comoving = region.mass_density[i],
        .momentum_density_x_comoving = region.momentum_x[i],
        .momentum_density_y_comoving = region.momentum_y[i],
        .momentum_density_z_comoving = region.momentum_z[i],
        .total_energy_density_comoving = region.total_energy[i],
        .metal_mass_density_comoving = region.metal_density[i]});
    if (soa.auxiliaryInternalEnergyEnabled()) {
      soa.auxiliaryInternalEnergyDensityComoving()[i] =
          region.mass_density[i] * region.specific_internal_energy[i];
    }
  }
}

void soaToRegion(const HydroConservedStateSoa& soa, HydroShadowPrivateRegion& region) {
  for (std::size_t i = 0; i < region.cell_count; ++i) {
    const auto cell = soa.loadCell(i);
    region.mass_density[i] = cell.mass_density_comoving;
    region.momentum_x[i] = cell.momentum_density_x_comoving;
    region.momentum_y[i] = cell.momentum_density_y_comoving;
    region.momentum_z[i] = cell.momentum_density_z_comoving;
    region.metal_density[i] = cell.metal_mass_density_comoving;
    region.specific_internal_energy[i] = selectedThermalEnergy(soa, i);
  }
  recomputeRegionEnergy(region);
}

// target += weight * (b - a) on the authoritative lanes and composition.
void addRegionDifference(HydroShadowPrivateRegion& target, const HydroShadowPrivateRegion& a,
                         const HydroShadowPrivateRegion& b, double weight) {
  if (a.cell_count != target.cell_count || b.cell_count != target.cell_count) {
    throw std::runtime_error("shadow chronology forcing has inconsistent region extents");
  }
  for (std::size_t i = 0; i < target.cell_count; ++i) {
    target.mass_density[i] += weight * (b.mass_density[i] - a.mass_density[i]);
    target.momentum_x[i] += weight * (b.momentum_x[i] - a.momentum_x[i]);
    target.momentum_y[i] += weight * (b.momentum_y[i] - a.momentum_y[i]);
    target.momentum_z[i] += weight * (b.momentum_z[i] - a.momentum_z[i]);
    target.metal_density[i] += weight * (b.metal_density[i] - a.metal_density[i]);
    target.specific_internal_energy[i] += weight * (b.specific_internal_energy[i] - a.specific_internal_energy[i]);
  }
  if (target.composition.enabled()) {
    if (!a.composition.enabled() || !b.composition.enabled() ||
        a.composition.laneCount() != target.composition.laneCount() ||
        b.composition.laneCount() != target.composition.laneCount()) {
      throw std::runtime_error("shadow chronology forcing has inconsistent composition layouts");
    }
    for (std::size_t l = 0; l < target.composition.laneCount(); ++l) {
      for (std::size_t i = 0; i < target.cell_count; ++i) {
        const double d = b.composition.value(l, i) - a.composition.value(l, i);
        if (d != 0.0) target.composition.addValue(l, i, weight * d);
      }
    }
  }
  recomputeRegionEnergy(target);
}

[[nodiscard]] std::size_t storageSlotMap(const HydroPatchGeometry& geometry, std::vector<std::size_t>& map) {
  map.assign(geometry.totalCellStorageCount(), k_invalid_ghost_cell_slot);
  for (std::size_t slot = 0; slot < geometry.ghost_cells.size(); ++slot) {
    const std::size_t row = geometry.ghost_cells[slot].ghost_cell;
    if (row < map.size()) map[row] = slot;
  }
  return map.size();
}

// Sub-region of the listed rows (production smooth sources advance only the
// scheduler-active rows; inactive rows are untouched by the private replay).
[[nodiscard]] HydroShadowPrivateRegion extractRows(const HydroShadowPrivateRegion& region,
                                                   std::span<const std::size_t> rows) {
  HydroShadowPrivateRegion sub;
  sub.resize(rows.size());
  sub.cell_volume_comoving = region.cell_volume_comoving;
  if (region.composition.enabled()) sub.composition = region.composition.emptyLike(rows.size(), false);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const std::size_t r = rows[i];
    sub.gas_cell_ids[i] = region.gas_cell_ids[r];
    sub.mass_density[i] = region.mass_density[r];
    sub.momentum_x[i] = region.momentum_x[r];
    sub.momentum_y[i] = region.momentum_y[r];
    sub.momentum_z[i] = region.momentum_z[r];
    sub.total_energy[i] = region.total_energy[r];
    sub.metal_density[i] = region.metal_density[r];
    sub.specific_internal_energy[i] = region.specific_internal_energy[r];
    if (region.composition.enabled()) {
      sub.composition.unpackRow(i, region.composition.packRow(r), region.gas_cell_ids[r]);
      sub.composition.mean_molecular_weight[i] = region.composition.mean_molecular_weight[r];
    }
  }
  return sub;
}

void scatterRows(const HydroShadowPrivateRegion& sub, std::span<const std::size_t> rows,
                 HydroShadowPrivateRegion& region) {
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const std::size_t r = rows[i];
    region.mass_density[r] = sub.mass_density[i];
    region.momentum_x[r] = sub.momentum_x[i];
    region.momentum_y[r] = sub.momentum_y[i];
    region.momentum_z[r] = sub.momentum_z[i];
    region.total_energy[r] = sub.total_energy[i];
    region.metal_density[r] = sub.metal_density[i];
    region.specific_internal_energy[r] = sub.specific_internal_energy[i];
    if (region.composition.enabled() && sub.composition.enabled()) {
      region.composition.unpackRow(r, sub.composition.packRow(i), region.gas_cell_ids[r]);
      region.composition.mean_molecular_weight[r] = sub.composition.mean_molecular_weight[i];
    }
  }
}

// Read-only view of the sample data the private engine needs. Built by
// HydroShadowReference members (the only code with access to Sample).
struct SampleView {
  std::array<double, 5> epoch_scale_factor{};
  std::array<double, 5> epoch_hubble_code{};
  const HydroThermodynamicClosure* closure = nullptr;
  std::span<const std::uint8_t> effective_closure_allowed{};
  double source_redshift = 0.0;
  std::span<const HydroCellSourceProperties> captured_properties{};
  bool source_properties_captured = false;
  double temperature_k_per_internal_energy_code = 0.0;
  double density_cgs_per_code = 0.0;
  bool comoving_density = false;
};

// Private halfstep engine: one HydroCoreSolver/Riemann/scratch built from the
// sample's production policies. Every call touches only private storage.
class PrivateHydroEngine {
 public:
  explicit PrivateHydroEngine(const HydroConservedStateSoa& initial,
                              const HydroPatchGeometry& geometry,
                              std::span<const std::size_t> active_cells,
                              std::span<const std::size_t> active_faces,
                              const HydroNumericsPolicy& numerics,
                              const HydroReconstructionPolicy& reconstruction_policy,
                              HydroRiemannSolverKind riemann_kind)
      : m_initial(initial), m_geometry(geometry),
        m_active{.active_cells = active_cells, .active_faces = active_faces},
        m_reconstruction_policy(reconstruction_policy),
        m_solver(reconstruction_policy.adiabatic_index, [&numerics] {
          HydroNumericsPolicy policy = numerics;
          policy.integration_scheme = HydroIntegrationScheme::kCornerTransportUpwind;
          return policy;
        }()),
        m_riemann(makeHydroRiemannSolver(riemann_kind)) {
    m_batch = m_solver.requiredWorkspaceFeatures(HydroActiveBatchPolicy{.max_active_cells = 4096U});
    (void)m_scratch.prepareForActiveSet(
        active_cells.size(), active_faces.size(), geometry.cellCount(),
        geometry.totalCellStorageCount(), geometry.faceGhostSlotExtent(), m_batch);
  }

  // Writes the ghost rows for substep `substep` of `substeps`. Substep 0 uses
  // the exact production step-start ghosts for EVERY slot (identical boundary
  // data to the production advance). Later substeps refresh self-authority
  // slots from the private state and interpolate endpoint-authority slots to
  // the substep start time; step-start-only slots keep their value.
  void setGhosts(HydroConservedStateSoa& state,
                 std::span<const HydroShadowGhostTimeAuthority> authority,
                 std::span<const HydroConservedState> endpoint,
                 std::span<const double> endpoint_auxiliary,
                 std::uint32_t substep, std::uint32_t substeps) const {
    const double fraction = static_cast<double>(substep) / static_cast<double>(substeps);
    if (substep != 0U) {
      fillHydroBoundaryGhostCells(state, m_geometry, m_reconstruction_policy.adiabatic_index);
      if (state.auxiliaryInternalEnergyEnabled()) {
        fillHydroBoundaryGhostAuxiliaryInternalEnergy(state, m_geometry);
      }
    }
    for (std::size_t slot = 0; slot < m_geometry.ghost_cells.size(); ++slot) {
      const std::size_t row = m_geometry.ghost_cells[slot].ghost_cell;
      const bool self = slot < authority.size() &&
          authority[slot] == HydroShadowGhostTimeAuthority::kSelfRefreshed;
      if (substep != 0U && self) continue;
      HydroConservedState value = m_initial.loadCell(row);
      double auxiliary = m_initial.auxiliaryInternalEnergyEnabled()
          ? m_initial.auxiliaryInternalEnergyDensityComoving()[row]
          : std::numeric_limits<double>::quiet_NaN();
      if (substep != 0U && slot < authority.size() &&
          authority[slot] == HydroShadowGhostTimeAuthority::kEndpointInterpolated &&
          slot < endpoint.size()) {
        value = (1.0 - fraction) * value + fraction * endpoint[slot];
        const double end_aux = slot < endpoint_auxiliary.size()
            ? endpoint_auxiliary[slot] : std::numeric_limits<double>::quiet_NaN();
        auxiliary = std::isfinite(auxiliary) && std::isfinite(end_aux)
            ? (1.0 - fraction) * auxiliary + fraction * end_aux
            : std::numeric_limits<double>::quiet_NaN();
      }
      state.storeCell(row, value);
      if (state.auxiliaryInternalEnergyEnabled()) {
        state.auxiliaryInternalEnergyDensityComoving()[row] = auxiliary;
      }
    }
  }

  // One private CTU/PPM substep over [f_begin, f_end] of `update`.
  void advance(HydroConservedStateSoa& state, const SampleView& sample,
               const HydroUpdateContext& update, std::uint32_t substep, std::uint32_t substeps,
               std::span<const HydroSourceTerm* const> sources,
               const HydroConservedStateSoa& property_initial,
               core::GasCompositionSidecar* composition,
               HydroFaceMassFluxSink* sink, std::uint64_t& fallbacks) {
    HydroUpdateContext sub = update;
    if (substeps > 1U) {
      const double f0 = static_cast<double>(substep) / static_cast<double>(substeps);
      const double f1 = static_cast<double>(substep + 1U) / static_cast<double>(substeps);
      const double fm = 0.5 * (f0 + f1);
      const std::size_t i0 = static_cast<std::size_t>(std::lround(4.0 * f0));
      const std::size_t i1 = static_cast<std::size_t>(std::lround(4.0 * f1));
      const std::size_t im = static_cast<std::size_t>(std::lround(4.0 * fm));
      sub.dt_code = update.dt_code * (f1 - f0);
      sub.time_begin_code = std::lerp(update.time_begin_code, update.time_end_code, f0);
      sub.time_end_code = std::lerp(update.time_begin_code, update.time_end_code, f1);
      sub.scale_factor_begin = sample.epoch_scale_factor[i0];
      sub.scale_factor_end = sample.epoch_scale_factor[i1];
      sub.scale_factor = sample.epoch_scale_factor[im];
      sub.hubble_rate_begin_code = sample.epoch_hubble_code[i0];
      sub.hubble_rate_end_code = sample.epoch_hubble_code[i1];
      sub.hubble_rate_code = sample.epoch_hubble_code[im];
    }
    HydroReconstructionPolicy policy = m_reconstruction_policy;
    policy.dt_over_cell_width_code = hydroTransportCourantFactors(
        sub.dt_code, sub.scale_factor, m_geometry.cell_width_x_comoving,
        m_geometry.cell_width_y_comoving, m_geometry.cell_width_z_comoving);
    policy.dt_over_dx_code = policy.dt_over_cell_width_code[0];
    PpmReconstruction reconstruction(policy);
    HydroSourceContext source{};
    source.update = sub;
    source.thermodynamic_closure = sample.closure;
    source.effective_closure_allowed = sample.effective_closure_allowed;
    source.redshift = sub.comoving_coordinates && sub.scale_factor > 0.0
        ? std::max(0.0, 1.0 / sub.scale_factor - 1.0) : sample.source_redshift;
    source.face_mass_flux_sink = sink;
    ShadowSourceProperties provider(
        state, property_initial, sample.captured_properties,
        update.scale_factor, sub.scale_factor, sub.comoving_coordinates,
        sample.temperature_k_per_internal_energy_code, sample.density_cgs_per_code,
        sample.comoving_density, composition);
    if (sample.source_properties_captured) source.source_property_provider = &provider;
    HydroProfileEvent profile;
    m_solver.advancePatchActiveSetBatchedWithScratch(
        state, m_geometry, m_active, sub, reconstruction, *m_riemann, sources,
        source, m_scratch, nullptr, m_batch, &profile, nullptr);
    fallbacks += profile.riemann_fallback_count + profile.positivity_fallback_count +
        profile.stencil_degradation_count;
  }

  [[nodiscard]] std::uint64_t scratchCapacityBytes() const { return m_scratch.ownedCapacityBytes(); }

 private:
  const HydroConservedStateSoa& m_initial;
  const HydroPatchGeometry& m_geometry;
  HydroActiveSetView m_active;
  HydroReconstructionPolicy m_reconstruction_policy;
  HydroCoreSolver m_solver;
  std::unique_ptr<HydroRiemannSolver> m_riemann;
  HydroActiveBatchPolicy m_batch{};
  HydroScratchBuffers m_scratch;
};

template <typename T>
void putWire(std::vector<std::uint8_t>& wire, T value) {
  static_assert(sizeof(T) == 8U);
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  for (std::size_t byte = 0; byte < 8U; ++byte) wire.push_back(static_cast<std::uint8_t>(bits >> (8U * byte)));
}

template <typename T>
[[nodiscard]] T getWire(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  static_assert(sizeof(T) == 8U);
  if (offset + 8U > bytes.size()) throw std::runtime_error("truncated shadow region result wire");
  std::uint64_t bits = 0U;
  for (std::size_t byte = 0; byte < 8U; ++byte) {
    bits |= static_cast<std::uint64_t>(bytes[offset + byte]) << (8U * byte);
  }
  offset += 8U;
  T value{};
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

}  // namespace

// ---------------------------------------------------------------------------
// Identity and ranking (C04-04)
// ---------------------------------------------------------------------------

std::uint64_t hydroShadowPhysicalRegionKey(
    int level, const std::array<std::int64_t, 3>& lattice_origin,
    const std::array<std::uint32_t, 3>& extent) noexcept {
  std::uint64_t key = mix64(0x43303421ULL ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(level)));
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    key = mix64(key ^ static_cast<std::uint64_t>(lattice_origin[axis]));
    key = mix64(key ^ (static_cast<std::uint64_t>(extent[axis]) << 1U));
  }
  return key == 0U ? 1U : key;
}

std::optional<std::array<std::int64_t, 3>> hydroShadowLatticeOrigin(
    const std::array<double, 3>& origin_comoving,
    const std::array<double, 3>& domain_origin_comoving,
    const std::array<double, 3>& cell_width_comoving) noexcept {
  std::array<std::int64_t, 3> result{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const double width = cell_width_comoving[axis];
    if (!(std::isfinite(width) && width > 0.0)) return std::nullopt;
    const double cells = (origin_comoving[axis] - domain_origin_comoving[axis]) / width;
    if (!std::isfinite(cells) || std::abs(cells) > 9.0e15) return std::nullopt;
    const double rounded = std::nearbyint(cells);
    if (std::abs(cells - rounded) > 1.0e-6) return std::nullopt;
    result[axis] = static_cast<std::int64_t>(rounded);
  }
  return result;
}

std::uint64_t hydroShadowEffectiveKey(const HydroShadowCandidate& candidate) noexcept {
  if (candidate.region_key != 0U) return candidate.region_key;
  const std::uint64_t legacy = mix64(candidate.patch_id ^ 0x6c656761637900ULL);
  return legacy == 0U ? 1U : legacy;
}

bool hydroShadowCandidateLess(
    const HydroShadowCandidate& lhs, const HydroShadowCandidate& rhs,
    std::uint64_t step_index, std::uint64_t seed) noexcept {
  return std::tuple{candidateScore(lhs, step_index, seed), hydroShadowEffectiveKey(lhs), lhs.level} <
         std::tuple{candidateScore(rhs, step_index, seed), hydroShadowEffectiveKey(rhs), rhs.level};
}

std::vector<HydroShadowCandidate> rankHydroShadowCandidates(
    std::span<const HydroShadowCandidate> candidates, std::size_t k,
    std::uint64_t step_index, std::uint64_t seed) {
  if (k == 0U || candidates.empty()) return {};
  std::vector<HydroShadowCandidate> ranked(candidates.begin(), candidates.end());
  std::sort(ranked.begin(), ranked.end(), [step_index, seed](const auto& lhs, const auto& rhs) {
    return hydroShadowCandidateLess(lhs, rhs, step_index, seed);
  });
  ranked.erase(std::unique(ranked.begin(), ranked.end(), sameCandidate), ranked.end());
  if (ranked.size() > k) ranked.resize(k);
  return ranked;
}

double hydroShadowScaleFactorAtPhysicalFraction(
    const core::LambdaCdmBackground* background, double begin, double end, double fraction) {
  if (!(fraction >= 0.0 && fraction <= 1.0)) {
    throw std::invalid_argument("hydro shadow physical-time fraction must lie in [0,1]");
  }
  if (fraction == 0.0 || begin == end) return begin;
  if (fraction == 1.0) return end;
  if (background == nullptr || !(begin > 0.0) || !(end > begin)) return begin + fraction * (end - begin);
  const double total = background->cosmicTimeIntervalSi(begin, end);
  if (!(std::isfinite(total) && total >= 0.0)) {
    throw std::runtime_error("hydro shadow interval has an invalid cosmological elapsed time");
  }
  const double target = fraction * total;
  double lower = begin;
  double upper = end;
  for (std::uint32_t iteration = 0U; iteration < 64U; ++iteration) {
    const double midpoint = lower + 0.5 * (upper - lower);
    if (background->cosmicTimeIntervalSi(begin, midpoint) < target) {
      lower = midpoint;
    } else {
      upper = midpoint;
    }
  }
  return lower + 0.5 * (upper - lower);
}

std::string_view hydroShadowComparisonClassToString(HydroShadowComparisonClass value) noexcept {
  switch (value) {
    case HydroShadowComparisonClass::kTimeOnlyCtuPpm: return "time_only_ctu_ppm_subdivision";
    case HydroShadowComparisonClass::kMixedSpatialTemporal: return "mixed_spatial_temporal_scheme";
    case HydroShadowComparisonClass::kCoupledFullPhysicsStrang: return "coupled_full_physics_strang_interval";
    case HydroShadowComparisonClass::kBoundedProjection: return "bounded_projection";
  }
  return "unknown";
}

std::string_view hydroShadowGhostProvenanceToString(HydroShadowGhostProvenance value) noexcept {
  switch (value) {
    case HydroShadowGhostProvenance::kSelfPhysical: return "self_physical";
    case HydroShadowGhostProvenance::kSelfPeriodicWrap: return "self_periodic_wrap";
    case HydroShadowGhostProvenance::kImportedSameLevel: return "imported_same_level";
    case HydroShadowGhostProvenance::kImportedSameLevelRemote: return "imported_same_level_mpi";
    case HydroShadowGhostProvenance::kImportedEdgeDiagonal: return "imported_edge_diagonal";
    case HydroShadowGhostProvenance::kCoarseToFine: return "coarse_to_fine";
    case HydroShadowGhostProvenance::kCoarseToFineTemporal: return "coarse_to_fine_temporal";
    case HydroShadowGhostProvenance::kFineToCoarse: return "fine_to_coarse";
    case HydroShadowGhostProvenance::kPhysicalClosureMapped: return "physical_closure_mapped";
    case HydroShadowGhostProvenance::kUnknownImported: return "unknown_imported";
  }
  return "unknown";
}

HydroShadowGhostProvenance hydroShadowDefaultGhostProvenance(
    const HydroPatchGeometry& geometry, std::size_t ghost_slot) noexcept {
  if (ghost_slot >= geometry.ghost_cells.size()) return HydroShadowGhostProvenance::kUnknownImported;
  const HydroGhostCell& ghost = geometry.ghost_cells[ghost_slot];
  if (ghost.transverse_axis < 3U) return HydroShadowGhostProvenance::kImportedEdgeDiagonal;
  if (ghost.mutation_rights == HydroGhostMutationRights::kWritablePhysicalBoundaryScratch &&
      ghost.boundary_kind != HydroBoundaryKind::kImportedMpi &&
      ghost.boundary_kind != HydroBoundaryKind::kInterior &&
      ghost.owner_real_cell < geometry.cellCount() &&
      ghost.source_real_cell < geometry.cellCount()) {
    return ghost.boundary_kind == HydroBoundaryKind::kPeriodic
        ? HydroShadowGhostProvenance::kSelfPeriodicWrap
        : HydroShadowGhostProvenance::kSelfPhysical;
  }
  return HydroShadowGhostProvenance::kUnknownImported;
}

std::vector<std::uint8_t> hydroShadowDependencyConeMask(
    const HydroPatchGeometry& geometry, std::span<const std::uint8_t> stale_slot, int radius) {
  const std::size_t count = geometry.cellCount();
  std::vector<std::uint8_t> mask(count, 0U);
  if (radius < 0) return mask;
  const std::array<std::size_t, 3> dims{geometry.nx, geometry.ny, geometry.nz};
  std::vector<std::size_t> owners;
  for (std::size_t slot = 0; slot < geometry.ghost_cells.size() && slot < stale_slot.size(); ++slot) {
    if (stale_slot[slot] != 0U && geometry.ghost_cells[slot].owner_real_cell < count) {
      owners.push_back(geometry.ghost_cells[slot].owner_real_cell);
    }
  }
  std::sort(owners.begin(), owners.end());
  owners.erase(std::unique(owners.begin(), owners.end()), owners.end());
  const auto r = static_cast<std::int64_t>(radius);
  for (const std::size_t owner : owners) {
    const auto o = geometry.cellIjk(owner);
    std::array<std::int64_t, 3> lo{};
    std::array<std::int64_t, 3> hi{};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const auto extent = static_cast<std::int64_t>(dims[axis]);
      const auto center = static_cast<std::int64_t>(o[axis]);
      lo[axis] = dims[axis] <= 1U ? 0 : std::max<std::int64_t>(0, center - r);
      hi[axis] = dims[axis] <= 1U ? 0 : std::min<std::int64_t>(extent - 1, center + r);
    }
    for (std::int64_t k = lo[2]; k <= hi[2]; ++k) {
      for (std::int64_t j = lo[1]; j <= hi[1]; ++j) {
        for (std::int64_t i = lo[0]; i <= hi[0]; ++i) {
          mask[geometry.linearCellIndex(static_cast<std::size_t>(i), static_cast<std::size_t>(j),
                                        static_cast<std::size_t>(k))] = 1U;
        }
      }
    }
  }
  return mask;
}

void HydroShadowPrivateRegion::resize(std::size_t count) {
  cell_count = count;
  gas_cell_ids.assign(count, 0U);
  for (auto* lane : {&mass_density, &momentum_x, &momentum_y, &momentum_z, &total_energy,
                     &metal_density, &specific_internal_energy}) {
    lane->assign(count, 0.0);
  }
}

std::uint64_t HydroShadowPrivateRegion::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = static_cast<std::uint64_t>(gas_cell_ids.capacity()) * sizeof(std::uint64_t);
  for (const auto* lane : {&mass_density, &momentum_x, &momentum_y, &momentum_z, &total_energy,
                           &metal_density, &specific_internal_energy}) {
    bytes += static_cast<std::uint64_t>(lane->capacity()) * sizeof(double);
  }
  return bytes + (composition.enabled() ? composition.ownedCapacityBytes() : 0U);
}

std::uint64_t HydroShadowPrivateRegion::estimatedBytes(std::size_t count, std::size_t composition_lanes) {
  // 7 double lanes + id lane + composition lanes + 13 composition metadata
  // doubles/bytes per row + id, plus a fixed allocator/object allowance.
  const std::uint64_t per_row = checkedAdd(checkedMul(8U + 1U, sizeof(double)),
      checkedMul(static_cast<std::uint64_t>(composition_lanes) + 16U, sizeof(double)));
  return checkedAdd(checkedMul(static_cast<std::uint64_t>(count), per_row), 4096U);
}

void appendHydroShadowRegionResultWire(const HydroShadowRegionResult& result, std::vector<std::uint8_t>& wire) {
  const std::size_t start = wire.size();
  putWire<std::uint64_t>(wire, result.region_key);
  putWire<std::int64_t>(wire, static_cast<std::int64_t>(result.level));
  putWire<std::uint64_t>(wire, static_cast<std::uint64_t>(result.comparison_class));
  putWire<std::uint64_t>(wire, result.compared_cells);
  putWire<std::uint64_t>(wire, result.compared_auxiliary_cells);
  putWire<std::uint64_t>(wire, result.cone_excluded_cells);
  putWire<std::uint64_t>(wire, result.compared_composition_cells);
  putWire<double>(wire, result.compared_volume);
  for (const auto& metric : result.components) {
    putWire<double>(wire, metric.weighted_abs_difference);
    putWire<double>(wire, metric.weighted_reference_magnitude);
    putWire<double>(wire, metric.linf_difference);
  }
  putWire<double>(wire, result.composition.weighted_abs_difference);
  putWire<double>(wire, result.composition.weighted_reference_magnitude);
  putWire<double>(wire, result.composition.linf_difference);
  for (double value : result.reference_delta) putWire<double>(wire, value);
  for (double value : result.production_delta) putWire<double>(wire, value);
  while (wire.size() - start < k_hydro_shadow_region_result_wire_bytes) wire.push_back(0U);
}

HydroShadowRegionResult readHydroShadowRegionResultWire(std::span<const std::uint8_t> bytes) {
  if (bytes.size() != k_hydro_shadow_region_result_wire_bytes) {
    throw std::runtime_error("hydro shadow region result wire has the wrong width");
  }
  std::size_t offset = 0U;
  HydroShadowRegionResult result;
  result.region_key = getWire<std::uint64_t>(bytes, offset);
  const std::int64_t level = getWire<std::int64_t>(bytes, offset);
  if (level < 0 || level > std::numeric_limits<int>::max()) {
    throw std::runtime_error("hydro shadow region result has an invalid level");
  }
  result.level = static_cast<int>(level);
  const std::uint64_t comparison_class = getWire<std::uint64_t>(bytes, offset);
  if (comparison_class >= k_hydro_shadow_comparison_class_count) {
    throw std::runtime_error("hydro shadow region result has an invalid comparison class");
  }
  result.comparison_class = static_cast<HydroShadowComparisonClass>(comparison_class);
  result.compared_cells = getWire<std::uint64_t>(bytes, offset);
  result.compared_auxiliary_cells = getWire<std::uint64_t>(bytes, offset);
  result.cone_excluded_cells = getWire<std::uint64_t>(bytes, offset);
  result.compared_composition_cells = getWire<std::uint64_t>(bytes, offset);
  result.compared_volume = getWire<double>(bytes, offset);
  for (auto& metric : result.components) {
    metric.weighted_abs_difference = getWire<double>(bytes, offset);
    metric.weighted_reference_magnitude = getWire<double>(bytes, offset);
    metric.linf_difference = getWire<double>(bytes, offset);
  }
  result.composition.weighted_abs_difference = getWire<double>(bytes, offset);
  result.composition.weighted_reference_magnitude = getWire<double>(bytes, offset);
  result.composition.linf_difference = getWire<double>(bytes, offset);
  for (double& value : result.reference_delta) value = getWire<double>(bytes, offset);
  for (double& value : result.production_delta) value = getWire<double>(bytes, offset);
  return result;
}

// ---------------------------------------------------------------------------
// Epoch, selection and capture
// ---------------------------------------------------------------------------

HydroShadowReference::HydroShadowReference(HydroShadowOptions options, core::MemoryGovernor* governor)
    : m_options(options), m_governor(governor) {
  if (options.cadence == 0U) {
    throw std::invalid_argument("hydro shadow requires a nonzero cadence");
  }
}

void HydroShadowReference::beginEpoch(std::uint64_t step_index) {
  discardPending();
  m_step_index = step_index;
  m_storage_cells_charged = 0U;
  m_workspace_bytes_charged = 0U;
  m_report = {};
  m_global_candidates.clear();
  m_selected_once.clear();
  m_region_results.clear();
}

void HydroShadowReference::setGlobalCandidates(std::vector<HydroShadowCandidate> ranked) {
  m_global_candidates = rankHydroShadowCandidates(ranked, m_options.max_patches_per_rank,
                                                 m_step_index, m_options.seed);
  m_report.globally_selected_regions = m_global_candidates.size();
  m_selected_once.clear();
}

bool HydroShadowReference::isGloballySelected(std::uint64_t region_key, int level) const noexcept {
  const HydroShadowCandidate probe{.patch_id = 0U, .level = level, .region_key = region_key};
  return std::any_of(m_global_candidates.begin(), m_global_candidates.end(),
                     [&](const auto& candidate) { return sameCandidate(candidate, probe); });
}

HydroShadowReference::Sample HydroShadowReference::capture(
    std::uint64_t patch_id, int level,
    const HydroConservedStateSoa& before, const HydroPatchGeometry& geometry,
    const HydroActiveSetView& active, const HydroUpdateContext& update,
    const HydroNumericsPolicy& production_numerics,
    HydroReconstructionScheme production_reconstruction,
    const HydroReconstructionPolicy& production_policy,
    HydroRiemannSolverKind production_riemann,
    std::span<const HydroSourceTerm* const> source_terms,
    const HydroSourceContext& source_context) {
  HydroShadowPatchCapture spec;
  spec.identity = HydroShadowCandidate{.patch_id = patch_id, .level = level, .region_key = 0U};
  spec.before = &before;
  spec.geometry = &geometry;
  spec.active = active;
  spec.update = update;
  spec.production_numerics = production_numerics;
  spec.production_reconstruction = production_reconstruction;
  spec.production_policy = production_policy;
  spec.production_riemann = production_riemann;
  spec.source_terms = source_terms;
  spec.source_context = &source_context;
  return capture(spec);
}

HydroShadowReference::Sample HydroShadowReference::capture(const HydroShadowPatchCapture& spec) {
  Sample sample;
  if (!m_options.enabled) return sample;
  ++m_report.eligible;
  if ((m_step_index % m_options.cadence) != 0U) {
    ++m_report.skipped_cadence;
    return sample;
  }
  if (spec.before == nullptr || spec.geometry == nullptr || spec.source_context == nullptr) {
    ++m_report.skipped_geometry;
    return sample;
  }
  const HydroConservedStateSoa& before = *spec.before;
  const HydroPatchGeometry& geometry = *spec.geometry;
  const HydroSourceContext& source_context = *spec.source_context;
  // Global rank-independent selection, prepared once before conditional
  // patch loops. A selected physical region may be sampled only once per
  // epoch, also when recursive AMR invokes the same patch several times.
  const HydroShadowCandidate identity = spec.identity;
  const bool domain_selected = std::any_of(
      m_global_candidates.begin(), m_global_candidates.end(),
      [&](const auto& candidate) { return sameCandidate(candidate, identity); });
  const bool region_selected = std::any_of(spec.regions.begin(), spec.regions.end(), [&](const auto& region) {
    return std::any_of(m_global_candidates.begin(), m_global_candidates.end(),
                       [&](const auto& candidate) { return sameCandidate(candidate, region.identity); });
  });
  if ((!domain_selected && !region_selected) ||
      std::any_of(m_selected_once.begin(), m_selected_once.end(),
                  [&](const auto& candidate) { return sameCandidate(candidate, identity); })) {
    ++m_report.skipped_selection;
    return sample;
  }
  ++m_report.attempted;
  // Source terms are cloned and evaluated on the independent reference.
  // Discrete or unknown source types opt out instead of mutating live ledgers.
  if (std::any_of(spec.source_terms.begin(), spec.source_terms.end(),
                  [](const auto* source) { return source == nullptr; })) {
    ++m_report.skipped_unsupported_source;
    ++m_report.skipped_sources_or_time;
    return sample;
  }
  // The production composition face sink is replaced by a PRIVATE sink over
  // private composition rows. Without the row maps the composition cannot be
  // replayed and the comparison would silently mislabel passives: refuse.
  const bool composition = source_context.face_mass_flux_sink != nullptr;
  if (composition && (spec.composition == nullptr || !spec.composition->enabled() ||
                      spec.composition_real_rows.size() < geometry.cellCount())) {
    ++m_report.skipped_external_passives;
    return sample;
  }
  const std::size_t count = geometry.cellCount();
  if (count == 0U || before.size() != geometry.totalCellStorageCount() ||
      spec.active.active_cells.empty()) {
    ++m_report.skipped_geometry;
    return sample;
  }
  // The production active set is authoritative: every active cell/face index
  // must be valid. Partial activity and suppressed duplicate periodic faces
  // are reproduced exactly (the private solver uses the same view).
  if (std::any_of(spec.active.active_cells.begin(), spec.active.active_cells.end(),
                  [count](std::size_t row) { return row >= count; }) ||
      std::any_of(spec.active.active_faces.begin(), spec.active.active_faces.end(),
                  [&](std::size_t face) { return face >= geometry.faces.size(); })) {
    ++m_report.skipped_geometry;
    return sample;
  }
  const bool temporal_refined =
      spec.production_numerics.integration_scheme == HydroIntegrationScheme::kCornerTransportUpwind &&
      spec.production_reconstruction == HydroReconstructionScheme::kPiecewiseParabolic;
  // Ghost provenance and time authority.
  const std::size_t slots = geometry.ghost_cells.size();
  std::vector<HydroShadowGhostProvenance> provenance(slots);
  for (std::size_t slot = 0; slot < slots; ++slot) {
    provenance[slot] = slot < spec.ghost_provenance.size()
        ? spec.ghost_provenance[slot] : hydroShadowDefaultGhostProvenance(geometry, slot);
  }
  bool any_imported = false;
  for (const auto value : provenance) any_imported = any_imported || !selfAuthority(value);
  // A two-halfstep reference needs a midpoint value for every non-self ghost.
  // Without a promised endpoint it is still legal: those slots become
  // step-start-only and their dependency cones are excluded. If that leaves
  // no comparable active cell, refuse up front with a specific reason.
  if (temporal_refined && any_imported && !spec.endpoint_ghosts_will_be_provided) {
    std::vector<std::uint8_t> stale(slots, 0U);
    for (std::size_t slot = 0; slot < slots; ++slot) stale[slot] = selfAuthority(provenance[slot]) ? 0U : 1U;
    const auto cone = hydroShadowDependencyConeMask(geometry, stale, k_hydro_shadow_halfstep_dependency_radius);
    if (std::none_of(spec.active.active_cells.begin(), spec.active.active_cells.end(),
                     [&](std::size_t row) { return cone[row] == 0U; })) {
      ++m_report.skipped_temporal_ghost;
      return sample;
    }
  }
  if (before.size() > m_options.max_storage_cells_per_rank -
          std::min(m_storage_cells_charged, m_options.max_storage_cells_per_rank)) {
    ++m_report.skipped_resource;
    return sample;
  }
  const std::size_t lanes = composition ? spec.composition->laneCount() : 0U;
  const HydroActiveBatchPolicy batch{.max_active_cells = 4096U,
                                    .ctu_transverse_workspace = true,
                                    .dual_energy_workspace =
                                        spec.production_numerics.dual_energy.mode != HydroDualEnergyMode::kOff};
  std::uint64_t requested = 0U;
  try {
    // Admission covers every retained or simultaneously live private buffer:
    // initial copy + working copy (+ chronology working copy) of all storage
    // rows (7 lanes), endpoint ghosts, production result, geometry/active
    // copies, source snapshots, composition rows/ghost fractions, slot maps,
    // two solver scratch envelopes (vector growth/allocator slack).
    const std::uint64_t lane_count = 6U + (before.auxiliaryInternalEnergyEnabled() ? 1U : 0U);
    const std::uint64_t state_bytes = checkedMul(checkedMul(static_cast<std::uint64_t>(before.size()),
                                                            lane_count), 3U * sizeof(double));
    std::uint64_t retained = checkedAdd(state_bytes,
        checkedMul(vectorBytes(slots, sizeof(HydroConservedState) + sizeof(double) + 3U), 1U));
    retained = checkedAdd(retained, vectorBytes(count, sizeof(HydroConservedState) + sizeof(double)));
    retained = checkedAdd(retained, vectorBytes(count, sizeof(HydroConservedState)));  // baseline
    retained = checkedAdd(retained, vectorBytes(geometry.faces.size(), sizeof(HydroFace)));
    retained = checkedAdd(retained, vectorBytes(geometry.ghost_cells.size(), sizeof(HydroGhostCell)));
    retained = checkedAdd(retained, vectorBytes(geometry.flux_register_faces.size(), sizeof(HydroFluxRegisterFace)));
    retained = checkedAdd(retained, vectorBytes(spec.active.active_cells.size() + spec.active.active_faces.size(),
                                                sizeof(std::size_t)));
    retained = checkedAdd(retained, vectorBytes(before.size(), sizeof(std::size_t)));  // slot map
    retained = checkedAdd(retained, vectorBytes(count, sizeof(HydroCellSourceProperties) + 2U));
    retained = checkedAdd(retained, checkedMul(spec.source_terms.size(),
                                               sizeof(std::unique_ptr<HydroSourceTerm>) + 512U));
    if (composition) {
      retained = checkedAdd(retained, checkedMul(2U, HydroShadowPrivateRegion::estimatedBytes(count, lanes)));
      retained = checkedAdd(retained, vectorBytes(checkedMul(slots, lanes + 1U), sizeof(double)));
      retained = checkedAdd(retained, vectorBytes(checkedMul(count, lanes), sizeof(double)));  // deltas
      retained = checkedAdd(retained, vectorBytes(spec.active.active_faces.size(), 4U * sizeof(std::size_t)));
    }
    std::size_t region_cells = 0U;
    for (const auto& region : spec.regions) region_cells += region.cells.size();
    retained = checkedAdd(retained, vectorBytes(region_cells, sizeof(std::size_t)));
    const std::uint64_t scratch_bytes = HydroScratchBuffers::workspaceBytes(
        count, spec.active.active_faces.size(), count, before.size(), geometry.faceGhostSlotExtent(), batch);
    requested = checkedAdd(retained, checkedMul(2U, scratch_bytes));
    if (requested > m_options.max_workspace_bytes_per_rank -
                        std::min(m_workspace_bytes_charged, m_options.max_workspace_bytes_per_rank)) {
      ++m_report.skipped_resource;
      return sample;
    }
    if (m_governor != nullptr) {
      sample.lease = m_governor->reserve(core::MemoryClass::kDiagnostic, requested,
                                         "hydro.shadow_reference.private_state_and_scratch");
      sample.lease.commit();
    }
    sample.identity = identity;
    if (spec.regions.empty()) {
      sample.regions.push_back(HydroShadowComparisonRegion{.identity = identity, .cells = {}});
    } else {
      for (const auto& region : spec.regions) {
        if (std::any_of(m_global_candidates.begin(), m_global_candidates.end(),
                        [&](const auto& candidate) { return sameCandidate(candidate, region.identity); })) {
          sample.regions.push_back(region);
        }
      }
    }
    sample.initial = before;
    sample.geometry = geometry;
    sample.active_cells.assign(spec.active.active_cells.begin(), spec.active.active_cells.end());
    sample.active_faces.assign(spec.active.active_faces.begin(), spec.active.active_faces.end());
    sample.update = spec.update;
    sample.numerics = spec.production_numerics;
    sample.production_reconstruction = spec.production_reconstruction;
    sample.shadow_sources.reserve(spec.source_terms.size());
    for (const HydroSourceTerm* source : spec.source_terms) {
      auto copy = source->cloneForShadow();
      if (!copy) {
        ++m_report.skipped_unsupported_source;
        ++m_report.skipped_sources_or_time;
        return {};
      }
      sample.shadow_sources.push_back(std::move(copy));
    }
    sample.closure = source_context.thermodynamic_closure;
    sample.source_redshift = source_context.redshift;
    if (!source_context.effective_closure_allowed.empty()) {
      // Reject malformed source authority rather than silently changing
      // the effective/resolved regime of any sampled real cell.
      if (source_context.effective_closure_allowed.size() < count) {
        ++m_report.skipped_unsupported_source;
        return {};
      }
      sample.effective_closure_allowed.assign(
          source_context.effective_closure_allowed.begin(),
          source_context.effective_closure_allowed.begin() + static_cast<std::ptrdiff_t>(count));
    }
    if (!spec.source_terms.empty()) {
      sample.captured_properties.reserve(count);
      for (std::size_t row = 0; row < count; ++row) {
        HydroCellSourceProperties properties{};
        if (source_context.source_property_provider != nullptr) {
          properties = source_context.source_property_provider->sourcePropertiesForCell(row);
        } else {
          if (row < source_context.mass_density_physical_cgs.size()) {
            properties.mass_density_physical_cgs = source_context.mass_density_physical_cgs[row];
          }
          if (row < source_context.hydrogen_number_density_cgs.size()) {
            properties.hydrogen_number_density_cgs = source_context.hydrogen_number_density_cgs[row];
          }
          if (row < source_context.metallicity_mass_fraction.size()) {
            properties.metallicity_mass_fraction = source_context.metallicity_mass_fraction[row];
          }
          if (row < source_context.temperature_k.size()) {
            properties.temperature_k = source_context.temperature_k[row];
          }
        }
        sample.captured_properties.push_back(properties);
      }
      sample.source_properties_captured = true;
    }
    sample.reconstruction_policy = spec.production_policy;
    sample.riemann_kind = spec.production_riemann;
    sample.provenance = std::move(provenance);
    sample.authority.resize(slots);
    for (std::size_t slot = 0; slot < slots; ++slot) {
      sample.authority[slot] = selfAuthority(sample.provenance[slot])
          ? HydroShadowGhostTimeAuthority::kSelfRefreshed
          : (spec.endpoint_ghosts_will_be_provided && temporal_refined
                 ? HydroShadowGhostTimeAuthority::kEndpointInterpolated
                 : HydroShadowGhostTimeAuthority::kStepStartOnly);
      ++m_report.ghost_slots_by_provenance[static_cast<std::size_t>(sample.provenance[slot])];
    }
    sample.closed_physical_boundary = !any_imported;
    sample.temporal_refined = temporal_refined;
    sample.needs_endpoint = temporal_refined && any_imported && spec.endpoint_ghosts_will_be_provided;
    sample.partial_active = spec.active.active_cells.size() != count;
    sample.source_active = !spec.source_terms.empty() || sample.closure != nullptr;
    sample.evolving_comoving = spec.update.comoving_coordinates &&
        (spec.update.scale_factor_begin != spec.update.scale_factor_end ||
         spec.update.hubble_rate_begin_code != spec.update.hubble_rate_end_code);
    sample.temperature_k_per_internal_energy_code = spec.temperature_k_per_internal_energy_code;
    sample.density_cgs_per_code = spec.density_cgs_per_code;
    sample.comoving_density = spec.comoving_density;
    sample.retain_for_chronology = spec.retain_for_chronology || composition ||
        chronologyOpen(hydroShadowEffectiveKey(identity));
    // Exact private halfstep epochs between the production endpoints.
    for (std::size_t q = 0; q < 5U; ++q) {
      const double f = static_cast<double>(q) / 4.0;
      const double a = spec.update.comoving_coordinates
          ? hydroShadowScaleFactorAtPhysicalFraction(spec.cosmology, spec.update.scale_factor_begin,
                                                     spec.update.scale_factor_end, f)
          : spec.update.scale_factor;
      sample.epoch_scale_factor[q] = std::max(a, 1.0e-12);
      sample.epoch_hubble_code[q] = spec.cosmology != nullptr && spec.time_si_per_code > 0.0 &&
              spec.update.comoving_coordinates
          ? spec.cosmology->hubbleSi(sample.epoch_scale_factor[q]) * spec.time_si_per_code
          : std::lerp(spec.update.hubble_rate_begin_code, spec.update.hubble_rate_end_code, f);
    }
    if (!spec.update.comoving_coordinates) {
      for (std::size_t q = 0; q < 5U; ++q) {
        sample.epoch_scale_factor[q] = spec.update.scale_factor;
        sample.epoch_hubble_code[q] = std::lerp(spec.update.hubble_rate_begin_code,
                                                spec.update.hubble_rate_end_code,
                                                static_cast<double>(q) / 4.0);
      }
    }
    if (composition) {
      const core::GasCompositionSidecar& production = *spec.composition;
      sample.composition = production.emptyLike(count, false);
      sample.gas_cell_ids.assign(count, 0U);
      for (std::size_t row = 0; row < count; ++row) {
        const std::uint32_t dense = spec.composition_real_rows[row];
        if (dense == core::kInvalidGasCellRow || dense >= production.rowCount()) {
          throw std::runtime_error("shadow composition: sampled real cell has no composition row");
        }
        const std::vector<double> packed = production.packRow(dense);
        sample.composition.unpackRow(row, packed, production.gas_cell_id[dense]);
        sample.composition.mean_molecular_weight[row] = production.mean_molecular_weight[dense];
        sample.gas_cell_ids[row] = production.gas_cell_id[dense];
      }
      sample.ghost_donor_fractions.assign(slots * lanes, std::numeric_limits<double>::quiet_NaN());
      for (std::size_t slot = 0; slot < slots; ++slot) {
        std::uint32_t dense = slot < spec.composition_ghost_rows.size()
            ? spec.composition_ghost_rows[slot] : core::kInvalidGasCellRow;
        if (selfAuthority(sample.provenance[slot])) {
          const std::size_t source = geometry.ghost_cells[slot].source_real_cell;
          dense = source < count ? spec.composition_real_rows[source] : core::kInvalidGasCellRow;
        }
        if (dense == core::kInvalidGasCellRow || dense >= production.rowCount()) continue;
        const double mass = compositionGasMassOf(production, dense);
        if (!(mass > 0.0)) continue;
        for (std::size_t l = 0; l < lanes; ++l) {
          sample.ghost_donor_fractions[slot * lanes + l] = production.value(l, dense) / mass;
        }
      }
    }
  } catch (const core::MemoryAdmissionError&) {
    ++m_report.skipped_resource;
    return {};
  } catch (const std::bad_alloc&) {
    ++m_report.skipped_resource;
    return {};
  } catch (const std::overflow_error&) {
    ++m_report.skipped_resource;
    return {};
  } catch (const std::exception&) {
    // A malformed source snapshot/closure is not a memory refusal.
    // Production is separately responsible for its own state validation.
    ++m_report.skipped_numerical;
    return {};
  }
  sample.admitted_bytes = requested;
  sample.valid = true;
  m_storage_cells_charged += before.size();
  m_workspace_bytes_charged += requested;
  m_selected_once.push_back(identity);
  ++m_report.selected;
  m_report.workspace_high_water_bytes = std::max(m_report.workspace_high_water_bytes, requested);
  return sample;
}

void HydroShadowReference::observeProduction(Sample&& sample, const HydroConservedStateSoa& production_after) {
  if (!sample) return;
  try {
    const std::size_t count = sample.geometry.cellCount();
    if (production_after.size() < count) throw std::runtime_error("shadow production result is truncated");
    sample.production_after.resize(count);
    sample.production_after_auxiliary.assign(count, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t row = 0; row < count; ++row) {
      sample.production_after[row] = production_after.loadCell(row);
      if (production_after.auxiliaryInternalEnergyEnabled()) {
        sample.production_after_auxiliary[row] = production_after.auxiliaryInternalEnergyDensityComoving()[row];
      }
    }
    sample.production_observed = true;
  } catch (const std::exception&) {
    ++m_report.skipped_numerical;
    return;
  }
  if (sample.needs_endpoint) {
    m_pending.push_back(std::move(sample));
    return;
  }
  runHydroComparison(sample);
  if (sample.retain_for_chronology && sample.valid) m_pending.push_back(std::move(sample));
}

void HydroShadowReference::compare(
    Sample&& sample, const HydroConservedStateSoa& production_after,
    const HydroPatchGeometry&, const HydroActiveSetView&, const HydroUpdateContext&,
    const HydroNumericsPolicy&, HydroReconstructionScheme) {
  observeProduction(std::move(sample), production_after);
  finalizePendingEndpoints();
}

std::vector<HydroShadowReference::EndpointRequest> HydroShadowReference::pendingEndpointRequests() const {
  std::vector<EndpointRequest> requests;
  for (const Sample& sample : m_pending) {
    if (sample.needs_endpoint && !sample.endpoint_received && !sample.hydro_compared) {
      requests.push_back(EndpointRequest{.identity = sample.identity,
                                         .time_end_code = sample.update.time_end_code,
                                         .ghost_slot_count = sample.geometry.ghost_cells.size()});
    }
  }
  return requests;
}

std::size_t HydroShadowReference::pendingEndpointCount() const noexcept {
  return static_cast<std::size_t>(std::count_if(m_pending.begin(), m_pending.end(), [](const Sample& sample) {
    return sample.needs_endpoint && !sample.endpoint_received && !sample.hydro_compared;
  }));
}

void HydroShadowReference::provideEndpointGhosts(
    std::uint64_t region_key, double time_end_code,
    std::span<const HydroConservedState> ghost_states,
    std::span<const double> ghost_auxiliary,
    std::span<const std::uint8_t> slot_valid) noexcept {
  for (Sample& sample : m_pending) {
    if (hydroShadowEffectiveKey(sample.identity) != region_key || sample.endpoint_received ||
        sample.hydro_compared) {
      continue;
    }
    const double tolerance = 1.0e-12 * std::max({1.0, std::abs(time_end_code),
                                                 std::abs(sample.update.time_end_code)});
    const std::size_t slots = sample.geometry.ghost_cells.size();
    if (std::abs(time_end_code - sample.update.time_end_code) > tolerance ||
        ghost_states.size() != slots || slot_valid.size() != slots ||
        (!ghost_auxiliary.empty() && ghost_auxiliary.size() != slots)) {
      ++m_report.skipped_endpoint_geometry_mismatch;
      sample.valid = false;
      return;
    }
    try {
      sample.endpoint_ghosts.assign(ghost_states.begin(), ghost_states.end());
      sample.endpoint_ghost_auxiliary.assign(slots, std::numeric_limits<double>::quiet_NaN());
      if (!ghost_auxiliary.empty()) {
        std::copy(ghost_auxiliary.begin(), ghost_auxiliary.end(), sample.endpoint_ghost_auxiliary.begin());
      }
      sample.endpoint_slot_valid.assign(slot_valid.begin(), slot_valid.end());
      for (std::size_t slot = 0; slot < slots; ++slot) {
        if (sample.authority[slot] == HydroShadowGhostTimeAuthority::kEndpointInterpolated &&
            (slot_valid[slot] == 0U || !finiteState(ghost_states[slot]))) {
          sample.authority[slot] = HydroShadowGhostTimeAuthority::kStepStartOnly;
        }
      }
      sample.endpoint_received = true;
    } catch (...) {
      ++m_report.skipped_resource;
      sample.valid = false;
    }
    return;
  }
}

void HydroShadowReference::finalizePendingEndpoints() noexcept {
  for (Sample& sample : m_pending) {
    if (!sample.valid || sample.hydro_compared || !sample.production_observed) continue;
    if (sample.needs_endpoint && !sample.endpoint_received) {
      ++m_report.skipped_endpoint_unavailable;
      for (auto& authority : sample.authority) {
        if (authority == HydroShadowGhostTimeAuthority::kEndpointInterpolated) {
          authority = HydroShadowGhostTimeAuthority::kStepStartOnly;
        }
      }
    }
    runHydroComparison(sample);
  }
  std::erase_if(m_pending, [](const Sample& sample) {
    return !sample.valid || (sample.hydro_compared && !sample.retain_for_chronology);
  });
}

void HydroShadowReference::accumulateRegionResult(const HydroShadowRegionResult& result) {
  m_region_results.push_back(result);
  for (std::size_t k = 0; k < 7U; ++k) {
    m_report.components[k].weighted_abs_difference += result.components[k].weighted_abs_difference;
    m_report.components[k].weighted_reference_magnitude += result.components[k].weighted_reference_magnitude;
    m_report.components[k].linf_difference =
        std::max(m_report.components[k].linf_difference, result.components[k].linf_difference);
  }
  for (std::size_t k = 0; k < 6U; ++k) {
    m_report.reference_delta[k] += result.reference_delta[k];
    m_report.production_delta[k] += result.production_delta[k];
  }
  m_report.compared_cells += result.compared_cells;
  m_report.compared_auxiliary_cells += result.compared_auxiliary_cells;
  m_report.compared_composition_cells += result.compared_composition_cells;
  m_report.cone_excluded_cells += result.cone_excluded_cells;
  m_report.compared_volume += result.compared_volume;
  ++m_report.completed_by_class[static_cast<std::size_t>(result.comparison_class)];
}

void HydroShadowReference::runHydroComparison(Sample& sample) noexcept {
  if (!sample.valid || sample.hydro_compared) return;
  sample.hydro_compared = true;
  // All sample work is private and non-authoritative, including Riemann
  // fallback counters; no production source terms, registers or flux sinks.
  try {
    const HydroPatchGeometry& geometry = sample.geometry;
    const std::size_t count = geometry.cellCount();
    HydroConservedStateSoa reference = sample.initial;
    SampleView view;
    view.epoch_scale_factor = sample.epoch_scale_factor;
    view.epoch_hubble_code = sample.epoch_hubble_code;
    view.closure = sample.closure;
    view.effective_closure_allowed = sample.effective_closure_allowed;
    view.source_redshift = sample.source_redshift;
    view.captured_properties = sample.captured_properties;
    view.source_properties_captured = sample.source_properties_captured;
    view.temperature_k_per_internal_energy_code = sample.temperature_k_per_internal_energy_code;
    view.density_cgs_per_code = sample.density_cgs_per_code;
    view.comoving_density = sample.comoving_density;
    PrivateHydroEngine engine(sample.initial, geometry, sample.active_cells, sample.active_faces,
                              sample.numerics, sample.reconstruction_policy, sample.riemann_kind);
    const std::uint32_t substeps = sample.temporal_refined ? 2U : 1U;
    std::vector<const HydroSourceTerm*> shadow_sources;
    shadow_sources.reserve(sample.shadow_sources.size());
    for (const auto& source : sample.shadow_sources) shadow_sources.push_back(source.get());
    // Private composition copy for the hydro-only interval.
    core::GasCompositionSidecar composition = sample.composition;
    std::vector<std::size_t> slot_by_row;
    (void)storageSlotMap(geometry, slot_by_row);
    std::vector<std::uint8_t> composition_unknown(count, 0U);
    std::optional<ShadowCompositionSink> sink;
    if (composition.enabled()) {
      sink.emplace(composition, geometry, sample.ghost_donor_fractions, slot_by_row);
      sink->reserve(sample.active_faces.size());
    }
    for (std::uint32_t substep = 0; substep < substeps; ++substep) {
      engine.setGhosts(reference, sample.authority, sample.endpoint_ghosts,
                       sample.endpoint_ghost_auxiliary, substep, substeps);
      engine.advance(reference, view, sample.update, substep, substeps, shadow_sources,
                     sample.initial, composition.enabled() ? &composition : nullptr,
                     sink.has_value() ? &*sink : nullptr, m_report.fallback_count);
      if (sink.has_value()) m_report.composition_remote_donor_faces += sink->apply(composition_unknown);
    }
    // Only the halfstep-2 boundary data can be stale: cone exclusion for
    // step-start-only slots of a two-halfstep reference.
    std::vector<std::uint8_t> stale(sample.authority.size(), 0U);
    std::uint64_t stale_slots = 0U;
    for (std::size_t slot = 0; slot < sample.authority.size(); ++slot) {
      switch (sample.authority[slot]) {
        case HydroShadowGhostTimeAuthority::kSelfRefreshed: ++m_report.ghost_slots_self_refreshed; break;
        case HydroShadowGhostTimeAuthority::kEndpointInterpolated:
          ++m_report.ghost_slots_endpoint_interpolated;
          break;
        case HydroShadowGhostTimeAuthority::kStepStartOnly:
          ++m_report.ghost_slots_step_start_only;
          if (substeps > 1U) {
            stale[slot] = 1U;
            ++stale_slots;
          }
          break;
      }
    }
    const std::vector<std::uint8_t> cone = stale_slots == 0U
        ? std::vector<std::uint8_t>(count, 0U)
        : hydroShadowDependencyConeMask(geometry, stale, k_hydro_shadow_halfstep_dependency_radius);
    std::uint64_t observed_capacity_bytes = checkedAdd(
        sample.initial.ownedCapacityBytes(), reference.ownedCapacityBytes());
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes, engine.scratchCapacityBytes());
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes,
        vectorBytes(sample.production_after.capacity(), sizeof(HydroConservedState)));
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes,
        vectorBytes(sample.endpoint_ghosts.capacity(), sizeof(HydroConservedState)));
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes,
        vectorBytes(sample.captured_properties.capacity(), sizeof(HydroCellSourceProperties)));
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes,
        vectorBytes(slot_by_row.capacity(), sizeof(std::size_t)));
    if (composition.enabled()) {
      observed_capacity_bytes = checkedAdd(observed_capacity_bytes,
          checkedAdd(composition.ownedCapacityBytes(), sample.composition.ownedCapacityBytes()));
    }
    m_report.observed_capacity_high_water_bytes = std::max(
        m_report.observed_capacity_high_water_bytes, observed_capacity_bytes);
    std::vector<std::uint8_t> active_flag(count, 0U);
    for (std::size_t row : sample.active_cells) active_flag[row] = 1U;
    const double volume = geometry.cell_volume_comoving;
    const bool dual_energy = reference.auxiliaryInternalEnergyEnabled();
    std::uint64_t completed_regions = 0U;
    for (const HydroShadowComparisonRegion& region : sample.regions) {
      HydroShadowRegionResult result;
      result.region_key = hydroShadowEffectiveKey(region.identity);
      result.level = region.identity.level;
      const auto visit = [&](std::size_t row) {
        if (row >= count || active_flag[row] == 0U) return;
        if (cone[row] != 0U) {
          ++result.cone_excluded_cells;
          return;
        }
        const auto& produced = sample.production_after[row];
        const auto expected = reference.loadCell(row);
        if (!finiteState(produced) || !finiteState(expected)) {
          throw std::runtime_error("non-finite state in shadow hydro comparison");
        }
        const double aux_ref = dual_energy ? reference.auxiliaryInternalEnergyDensityComoving()[row] : 0.0;
        const double aux_prod = dual_energy ? sample.production_after_auxiliary[row] : 0.0;
        const bool aux = dual_energy && std::isfinite(aux_ref) && std::isfinite(aux_prod);
        const auto left = lanes(produced, aux_prod);
        const auto right = lanes(expected, aux_ref);
        for (std::size_t k = 0; k < (aux ? 7U : 6U); ++k) addMetric(result.components[k], left[k], right[k], volume);
        if (aux) ++result.compared_auxiliary_cells;
        const auto baseline = physicalLanes(sample.initial.loadCell(row));
        const auto ref_lanes = physicalLanes(expected);
        const auto prod_lanes = physicalLanes(produced);
        for (std::size_t k = 0; k < 6U; ++k) {
          result.reference_delta[k] += (ref_lanes[k] - baseline[k]) * volume;
          result.production_delta[k] += (prod_lanes[k] - baseline[k]) * volume;
        }
        ++result.compared_cells;
        result.compared_volume += volume;
      };
      if (region.cells.empty()) {
        for (std::size_t row : sample.active_cells) visit(row);
      } else {
        for (std::size_t row : region.cells) visit(row);
      }
      result.comparison_class = hydroClassFor(sample.temporal_refined, result.cone_excluded_cells);
      if (result.compared_cells == 0U) {
        ++m_report.skipped_geometry;
        continue;
      }
      accumulateRegionResult(result);
      ++completed_regions;
    }
    // Retain the private hydro-only composition for the post-commit check
    // (unknown-donor cells lie inside the composition cone used there).
    if (composition.enabled()) sample.composition = std::move(composition);
    if (completed_regions == 0U) {
      sample.valid = false;
      return;
    }
    ++m_report.completed;
    if (sample.closed_physical_boundary) ++m_report.completed_closed_boundary_patches;
    if (sample.source_active) ++m_report.completed_source_active_patches;
    if (sample.evolving_comoving) ++m_report.completed_evolving_comoving_patches;
    if (sample.endpoint_received) ++m_report.completed_endpoint_ghost_patches;
    if (sample.partial_active) ++m_report.completed_partial_active_patches;
    m_report.reference_substeps += substeps;
  } catch (const std::exception&) {
    // Shadow is a diagnostic. A failed reference must not change a step,
    // acceptance decision, source history or physical model. Count and skip.
    ++m_report.skipped_numerical;
    sample.valid = false;
  }
}

// ---------------------------------------------------------------------------
// Full-physics chronology (C04-02)
// ---------------------------------------------------------------------------

bool HydroShadowReference::openChronology(const HydroShadowCandidate& domain,
                                          std::span<const std::uint64_t> gas_cell_ids,
                                          std::size_t composition_lanes, bool level_clock) {
  if (!m_options.enabled || gas_cell_ids.empty()) return false;
  const std::uint64_t key = hydroShadowEffectiveKey(domain);
  if (chronologyOpen(key)) return true;
  try {
    // Six retained point snapshots + working region + forcing copies.
    const std::uint64_t bytes = checkedMul(
        k_hydro_shadow_chronology_point_count + 4U,
        HydroShadowPrivateRegion::estimatedBytes(gas_cell_ids.size(), composition_lanes));
    if (bytes > m_options.max_workspace_bytes_per_rank -
                    std::min(m_workspace_bytes_charged, m_options.max_workspace_bytes_per_rank)) {
      ++m_report.chronology_skipped_resource;
      return false;
    }
    ChronologyRecord record;
    if (m_governor != nullptr) {
      record.lease = m_governor->reserve(core::MemoryClass::kDiagnostic, bytes,
                                         "hydro.shadow_reference.chronology_snapshots");
      record.lease.commit();
    }
    record.domain = domain;
    record.gas_cell_ids.assign(gas_cell_ids.begin(), gas_cell_ids.end());
    record.composition_lanes = composition_lanes;
    record.level_clock = level_clock;
    m_workspace_bytes_charged += bytes;
    m_chronologies.push_back(std::move(record));
    ++m_report.chronology_started;
    return true;
  } catch (const std::exception&) {
    ++m_report.chronology_skipped_resource;
    return false;
  }
}

bool HydroShadowReference::chronologyOpen(std::uint64_t domain_key) const noexcept {
  return std::any_of(m_chronologies.begin(), m_chronologies.end(), [&](const ChronologyRecord& record) {
    return hydroShadowEffectiveKey(record.domain) == domain_key;
  });
}

std::vector<HydroShadowCandidate> HydroShadowReference::openChronologyDomains() const {
  std::vector<HydroShadowCandidate> domains;
  domains.reserve(m_chronologies.size());
  for (const ChronologyRecord& record : m_chronologies) domains.push_back(record.domain);
  return domains;
}

std::span<const std::uint64_t> HydroShadowReference::chronologyGasCellIds(std::uint64_t domain_key) const noexcept {
  for (const ChronologyRecord& record : m_chronologies) {
    if (hydroShadowEffectiveKey(record.domain) == domain_key) return record.gas_cell_ids;
  }
  return {};
}

void HydroShadowReference::recordChronologyPoint(std::uint64_t domain_key, HydroShadowChronologyPoint point,
                                                 const HydroShadowPrivateRegion& snapshot) noexcept {
  for (ChronologyRecord& record : m_chronologies) {
    if (hydroShadowEffectiveKey(record.domain) != domain_key) continue;
    if (snapshot.cell_count != record.gas_cell_ids.size() || snapshot.gas_cell_ids != record.gas_cell_ids) {
      dropChronology(domain_key, &HydroShadowReport::chronology_skipped_region_changed);
      return;
    }
    try {
      record.points[static_cast<std::size_t>(point)] = snapshot;
    } catch (...) {
      dropChronology(domain_key, &HydroShadowReport::chronology_skipped_resource);
      return;
    }
    return;
  }
}

std::vector<HydroShadowCandidate> HydroShadowReference::compositionCheckDomains() const {
  std::vector<HydroShadowCandidate> domains;
  for (const Sample& sample : m_pending) {
    if (sample.valid && sample.hydro_compared && sample.composition.enabled() && !sample.composition_checked) {
      domains.push_back(sample.identity);
    }
  }
  return domains;
}

std::span<const std::uint64_t> HydroShadowReference::sampleGasCellIds(std::uint64_t domain_key) const noexcept {
  for (const Sample& sample : m_pending) {
    if (hydroShadowEffectiveKey(sample.identity) == domain_key) return sample.gas_cell_ids;
  }
  return {};
}

void HydroShadowReference::compareCompositionAfterCommit(std::uint64_t domain_key,
                                                         const HydroShadowPrivateRegion& snapshot,
                                                         double endpoint_time_code) noexcept {
  // Hydro-only composition check: the private composition transported with
  // the private mass fluxes vs the committed production composition at the
  // same physical endpoint (post-hydro, pre-regrid). Ghost donor fractions
  // are step-start values, so cells within the donor cone of any imported
  // slot are excluded (production uses step-start donors for the whole step;
  // the private second halfstep would need midpoint donors).
  for (Sample& sample : m_pending) {
    if (hydroShadowEffectiveKey(sample.identity) != domain_key || !sample.valid || !sample.hydro_compared ||
        !sample.composition.enabled() || sample.composition_checked || !snapshot.composition.enabled()) {
      continue;
    }
    sample.composition_checked = true;
    // A recursive fine-level sample ends inside the coarse interval: the
    // committed post-hydro composition is a later endpoint. Refuse.
    const double tolerance = 1.0e-9 * std::max({1.0, std::abs(endpoint_time_code),
                                                std::abs(sample.update.time_end_code)});
    if (!(std::abs(sample.update.time_end_code - endpoint_time_code) <= tolerance)) {
      ++m_report.skipped_composition_endpoint_mismatch;
      return;
    }
    try {
      if (snapshot.cell_count != sample.geometry.cellCount() ||
          snapshot.composition.laneCount() != sample.composition.laneCount()) {
        ++m_report.skipped_endpoint_geometry_mismatch;
        return;
      }
      HydroShadowRegionResult result;
      result.region_key = domain_key;
      result.level = sample.identity.level;
      result.comparison_class = sample.temporal_refined ? HydroShadowComparisonClass::kTimeOnlyCtuPpm
                                                        : HydroShadowComparisonClass::kMixedSpatialTemporal;
      std::vector<std::uint8_t> imported(sample.authority.size(), 0U);
      for (std::size_t slot = 0; slot < imported.size(); ++slot) {
        imported[slot] = sample.authority[slot] == HydroShadowGhostTimeAuthority::kSelfRefreshed ? 0U : 1U;
      }
      const auto cone = hydroShadowDependencyConeMask(sample.geometry, imported,
                                                      k_hydro_shadow_composition_dependency_radius);
      const std::size_t lanes = sample.composition.laneCount();
      for (std::size_t row : sample.active_cells) {
        if (row >= snapshot.cell_count || cone[row] != 0U) continue;
        for (std::size_t l = 0; l < lanes; ++l) {
          addMetric(result.composition, snapshot.composition.value(l, row), sample.composition.value(l, row), 1.0);
        }
        ++result.compared_composition_cells;
      }
      if (result.compared_composition_cells != 0U) {
        accumulateRegionResult(result);
        ++m_report.completed_composition_patches;
      }
    } catch (...) {
      ++m_report.skipped_numerical;
    }
    return;
  }
}

void HydroShadowReference::dropChronology(std::uint64_t domain_key, std::uint64_t HydroShadowReport::*reason) noexcept {
  const auto before = m_chronologies.size();
  std::erase_if(m_chronologies, [&](const ChronologyRecord& record) {
    return hydroShadowEffectiveKey(record.domain) == domain_key;
  });
  if (m_chronologies.size() != before && reason != nullptr) ++(m_report.*reason);
}

void HydroShadowReference::noteChronologyUnsupported(std::uint64_t HydroShadowReport::*reason,
                                                     std::uint64_t count) noexcept {
  if (reason != nullptr) m_report.*reason += count;
}

void HydroShadowReference::finalizeChronologies(const HydroShadowSmoothSourceReplay* replay,
                                                std::optional<int> level_filter) noexcept {
  std::vector<std::uint64_t> finished;
  for (ChronologyRecord& record : m_chronologies) {
    const std::uint64_t key = hydroShadowEffectiveKey(record.domain);
    if (level_filter.has_value() && (!record.level_clock || record.domain.level != *level_filter)) continue;
    finished.push_back(key);
    if (!level_filter.has_value() && record.level_clock) {
      // A level window that never reached its post half end in this attempt.
      ++m_report.chronology_skipped_missing_point;
      continue;
    }
    if (replay == nullptr) {
      ++m_report.chronology_skipped_no_replay_operator;
      continue;
    }
    Sample* sample = nullptr;
    for (Sample& candidate : m_pending) {
      if (hydroShadowEffectiveKey(candidate.identity) == key && candidate.hydro_compared &&
          candidate.valid) {
        sample = &candidate;
      }
    }
    if (sample == nullptr) {
      ++m_report.chronology_skipped_hydro_sample_missing;
      continue;
    }
    runChronology(record, *sample, *replay);
  }
  std::erase_if(m_chronologies, [&](const ChronologyRecord& record) {
    return std::find(finished.begin(), finished.end(), hydroShadowEffectiveKey(record.domain)) != finished.end();
  });
  if (!level_filter.has_value()) {
    std::erase_if(m_pending, [](const Sample& sample) { return sample.hydro_compared || !sample.valid; });
  }
}

void HydroShadowReference::runChronology(ChronologyRecord& record, Sample& sample,
                                         const HydroShadowSmoothSourceReplay& replay) noexcept {
  using P = HydroShadowChronologyPoint;
  const auto point = [&record](P p) -> const std::optional<HydroShadowPrivateRegion>& {
    return record.points[static_cast<std::size_t>(p)];
  };
  for (const P required : {P::kPreHalfBegin, P::kPreHalfAfterTransport, P::kPreHalfEnd,
                           P::kPostHalfAfterThermochemistry, P::kPostHalfAfterTransport, P::kPostHalfEnd}) {
    if (!point(required).has_value()) {
      ++m_report.chronology_skipped_missing_point;
      return;
    }
  }
  const std::size_t count = sample.geometry.cellCount();
  if (record.gas_cell_ids.size() != count) {
    ++m_report.chronology_skipped_region_changed;
    return;
  }
  try {
    const auto status_ok = [this](HydroShadowReplayStatus status) {
      switch (status) {
        case HydroShadowReplayStatus::kCompleted: return true;
        case HydroShadowReplayStatus::kModelUnavailable: ++m_report.chronology_skipped_replay_unavailable; return false;
        case HydroShadowReplayStatus::kIntervalIncomplete: ++m_report.chronology_skipped_replay_incomplete; return false;
        case HydroShadowReplayStatus::kInvalidState: ++m_report.skipped_numerical; return false;
      }
      return false;
    };
    // Production smooth sources advance exactly the scheduler-active rows of
    // the interval (the same active set as the sampled hydro advance); the
    // private replay is restricted to those rows, inactive rows are untouched.
    const auto smooth = [&](HydroShadowPrivateRegion& region, double f0, double f1) {
      if (!sample.partial_active) return replay.advanceSmooth(region, f0, f1);
      HydroShadowPrivateRegion sub = extractRows(region, sample.active_cells);
      const HydroShadowReplayStatus status = replay.advanceSmooth(sub, f0, f1);
      if (status == HydroShadowReplayStatus::kCompleted) scatterRows(sub, sample.active_cells, region);
      return status;
    };
    // Replay identity: the private S1 from the production post-transport
    // state must reproduce the production hydro input (same operator, same
    // data). A nonzero value measures replay infidelity, not physics error.
    {
      HydroShadowPrivateRegion probe = *point(P::kPreHalfAfterTransport);
      if (!status_ok(smooth(probe, 0.0, 0.5))) return;
      const HydroShadowPrivateRegion& target = *point(P::kPreHalfEnd);
      double identity = 0.0;
      for (std::size_t i = 0; i < count; ++i) {
        identity = std::max(identity, std::abs(probe.specific_internal_energy[i] - target.specific_internal_energy[i]) /
                                          std::max(std::abs(target.specific_internal_energy[i]), 1.0e-300));
      }
      m_report.chronology_replay_identity_linf = std::max(m_report.chronology_replay_identity_linf, identity);
    }
    const HydroShadowPrivateRegion& u0 = *point(P::kPreHalfBegin);
    HydroShadowPrivateRegion x = u0;
    HydroConservedStateSoa soa = sample.initial;
    SampleView view;
    view.epoch_scale_factor = sample.epoch_scale_factor;
    view.epoch_hubble_code = sample.epoch_hubble_code;
    view.closure = sample.closure;
    view.effective_closure_allowed = sample.effective_closure_allowed;
    view.source_redshift = sample.source_redshift;
    view.captured_properties = sample.captured_properties;
    view.source_properties_captured = sample.source_properties_captured;
    view.temperature_k_per_internal_energy_code = sample.temperature_k_per_internal_energy_code;
    view.density_cgs_per_code = sample.density_cgs_per_code;
    view.comoving_density = sample.comoving_density;
    PrivateHydroEngine engine(sample.initial, sample.geometry, sample.active_cells, sample.active_faces,
                              sample.numerics, sample.reconstruction_policy, sample.riemann_kind);
    std::vector<const HydroSourceTerm*> shadow_sources;
    for (const auto& source : sample.shadow_sources) shadow_sources.push_back(source.get());
    std::vector<std::size_t> slot_by_row;
    (void)storageSlotMap(sample.geometry, slot_by_row);
    std::vector<std::uint8_t> composition_unknown(count, 0U);
    const auto hydro_halfstep = [&](std::uint32_t substep) {
      regionToSoa(x, soa);
      engine.setGhosts(soa, sample.authority, sample.endpoint_ghosts,
                       sample.endpoint_ghost_auxiliary, substep, 2U);
      // Level windows: production records the composition face transfers
      // during the level sweep but commits them only at the coarse end, so
      // inside the window the production composition is NOT advected; the
      // private window mirrors that deferral exactly.
      std::optional<ShadowCompositionSink> sink;
      if (x.composition.enabled() && !record.level_clock) {
        sink.emplace(x.composition, sample.geometry, sample.ghost_donor_fractions, slot_by_row);
        sink->reserve(sample.active_faces.size());
      }
      engine.advance(soa, view, sample.update, substep, 2U, shadow_sources, sample.initial,
                     x.composition.enabled() ? &x.composition : nullptr,
                     sink.has_value() ? &*sink : nullptr, m_report.fallback_count);
      if (sink.has_value()) m_report.composition_remote_donor_faces += sink->apply(composition_unknown);
      soaToRegion(soa, x);
      // Production commits composition transport and then re-derives the
      // metal mirror from the inventory; the private chronology does too.
      if (x.composition.enabled() && !record.level_clock) {
        replay.synchronizeDerivedMirrors(x);
        recomputeRegionEnergy(x);
      }
    };
    const HydroShadowPrivateRegion& p1 = *point(P::kPreHalfAfterTransport);
    const HydroShadowPrivateRegion& p4 = *point(P::kPostHalfAfterThermochemistry);
    const HydroShadowPrivateRegion& p5 = *point(P::kPostHalfAfterTransport);
    const HydroShadowPrivateRegion& w = *point(P::kPostHalfEnd);
    // Refined T01 chronology over the same accepted interval [t0, t1]:
    //   [S(q0,q1) H(q0,q2) S(q1,q2)] [S(q2,q3) H(q2,q4) S(q3,q4)] + events.
    // Smooth thermochemistry/ISM/dust/SGS is recomputed privately with the
    // production operator; element diffusion/SGS shear (neighbour-coupled)
    // enter as the captured production increments apportioned linearly in
    // time over the matching quarter intervals; discrete events enter as the
    // captured production increment at the interval end (immutable forcing).
    addRegionDifference(x, u0, p1, 0.5);
    if (!status_ok(smooth(x, 0.0, 0.25))) return;
    hydro_halfstep(0U);
    if (!status_ok(smooth(x, 0.25, 0.5))) return;
    addRegionDifference(x, u0, p1, 0.5);
    addRegionDifference(x, p4, p5, 0.5);
    if (!status_ok(smooth(x, 0.5, 0.75))) return;
    hydro_halfstep(1U);
    if (!status_ok(smooth(x, 0.75, 1.0))) return;
    addRegionDifference(x, p4, p5, 0.5);
    addRegionDifference(x, p5, w, 1.0);
    // Comparable cells: outside the hydro dependency cone of stale slots;
    // composition additionally outside the donor cone of every imported slot
    // and away from faces with unknown ghost donors.
    std::vector<std::uint8_t> stale(sample.authority.size(), 0U);
    std::vector<std::uint8_t> imported(sample.authority.size(), 0U);
    for (std::size_t slot = 0; slot < stale.size(); ++slot) {
      stale[slot] = sample.authority[slot] == HydroShadowGhostTimeAuthority::kStepStartOnly ? 1U : 0U;
      imported[slot] = sample.authority[slot] == HydroShadowGhostTimeAuthority::kSelfRefreshed ? 0U : 1U;
    }
    const auto cone = hydroShadowDependencyConeMask(sample.geometry, stale, k_hydro_shadow_halfstep_dependency_radius);
    const auto composition_cone = hydroShadowDependencyConeMask(
        sample.geometry, imported, k_hydro_shadow_composition_dependency_radius);
    const double volume = sample.geometry.cell_volume_comoving;
    std::uint64_t completed_regions = 0U;
    for (const HydroShadowComparisonRegion& region : sample.regions) {
      HydroShadowRegionResult result;
      result.region_key = hydroShadowEffectiveKey(region.identity);
      result.level = region.identity.level;
      result.comparison_class = HydroShadowComparisonClass::kCoupledFullPhysicsStrang;
      const auto visit = [&](std::size_t i) {
        if (i >= count) return;
        if (cone[i] != 0U) {
          ++result.cone_excluded_cells;
          return;
        }
        const std::array<double, 7> produced{w.mass_density[i], w.momentum_x[i], w.momentum_y[i], w.momentum_z[i],
                                             w.total_energy[i], w.metal_density[i],
                                             w.mass_density[i] * w.specific_internal_energy[i]};
        const std::array<double, 7> expected{x.mass_density[i], x.momentum_x[i], x.momentum_y[i], x.momentum_z[i],
                                             x.total_energy[i], x.metal_density[i],
                                             x.mass_density[i] * x.specific_internal_energy[i]};
        for (std::size_t k = 0; k < 7U; ++k) addMetric(result.components[k], produced[k], expected[k], volume);
        ++result.compared_auxiliary_cells;
        const std::array<double, 6> base{u0.mass_density[i], u0.momentum_x[i], u0.momentum_y[i], u0.momentum_z[i],
                                         u0.total_energy[i], u0.metal_density[i]};
        for (std::size_t k = 0; k < 6U; ++k) {
          result.reference_delta[k] += (expected[k] - base[k]) * volume;
          result.production_delta[k] += (produced[k] - base[k]) * volume;
        }
        ++result.compared_cells;
        result.compared_volume += volume;
        if (x.composition.enabled() && w.composition.enabled() && composition_cone[i] == 0U &&
            composition_unknown[i] == 0U) {
          for (std::size_t l = 0; l < x.composition.laneCount(); ++l) {
            addMetric(result.composition, w.composition.value(l, i), x.composition.value(l, i), 1.0);
          }
          ++result.compared_composition_cells;
        }
      };
      if (region.cells.empty()) {
        for (std::size_t i = 0; i < count; ++i) visit(i);
      } else {
        for (std::size_t i : region.cells) visit(i);
      }
      if (result.compared_cells == 0U) {
        ++m_report.skipped_geometry;
        continue;
      }
      accumulateRegionResult(result);
      ++completed_regions;
    }
    if (completed_regions == 0U) return;
    ++m_report.completed_chronology_regions;
  } catch (const std::exception&) {
    ++m_report.skipped_numerical;
  }
}

void HydroShadowReference::discardPending() noexcept {
  for (const Sample& sample : m_pending) {
    if (sample.valid && !sample.hydro_compared) ++m_report.discarded_rejected_attempt;
  }
  m_report.discarded_rejected_attempt += m_chronologies.size();
  m_pending.clear();
  m_chronologies.clear();
}

}  // namespace cosmosim::hydro
