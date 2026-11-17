#include "cosmosim/hydro/hydro_shadow_reference.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <tuple>
#include <new>
#include <utility>
#include <numeric>
#include <stdexcept>
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
  return mix64(mix64(candidate.patch_id) ^ mix64(static_cast<std::uint32_t>(candidate.level)) ^
               mix64(step) ^ mix64(seed));
}

[[nodiscard]] bool sameCandidate(
    const HydroShadowCandidate& a, const HydroShadowCandidate& b) noexcept {
  return a.patch_id == b.patch_id && a.level == b.level;
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

[[nodiscard]] double specificThermalEnergy(const HydroConservedState& state) {
  const double rho = std::max(state.mass_density_comoving, 1.0e-30);
  const double momentum_sq = state.momentum_density_x_comoving * state.momentum_density_x_comoving +
      state.momentum_density_y_comoving * state.momentum_density_y_comoving +
      state.momentum_density_z_comoving * state.momentum_density_z_comoving;
  return std::max((state.total_energy_density_comoving - 0.5 * momentum_sq / rho) / rho, 0.0);
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

// The provider never reads SimulationState or an authoritative mutable cache.
// Composition fractions and thermodynamic conversion are frozen over the
// sampled HYDRO transport interval; density, temperature and cosmological
// dilution are reconstructed from each private intermediate state.
class ShadowSourceProperties final : public HydroCellSourcePropertyProvider {
 public:
  ShadowSourceProperties(const HydroConservedStateSoa& reference,
                         const HydroConservedStateSoa& initial,
                         std::span<const HydroCellSourceProperties> baseline,
                         double initial_scale, double substep_scale, bool comoving)
      : m_reference(reference), m_initial(initial), m_baseline(baseline),
        m_initial_scale(initial_scale), m_substep_scale(substep_scale), m_comoving(comoving) {}

  [[nodiscard]] HydroCellSourceProperties sourcePropertiesForCell(std::size_t cell_index) const override {
    const auto initial = m_initial.loadCell(cell_index);
    const auto current = m_reference.loadCell(cell_index);
    const auto baseline = m_baseline.at(cell_index);
    const double ratio = current.mass_density_comoving /
        std::max(initial.mass_density_comoving, 1.0e-30);
    const double dilution = m_comoving
        ? std::pow(m_initial_scale / m_substep_scale, 3.0) : 1.0;
    HydroCellSourceProperties result = baseline;
    result.mass_density_physical_cgs *= ratio * dilution;
    result.hydrogen_number_density_cgs *= ratio * dilution;
    const double u_initial = selectedThermalEnergy(m_initial, cell_index);
    const double u_current = selectedThermalEnergy(m_reference, cell_index);
    result.temperature_k *= u_current / std::max(u_initial, 1.0e-30);
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
};

[[nodiscard]] std::array<double, 7> lanes(const HydroConservedState& state, double auxiliary) {
  return {state.mass_density_comoving, state.momentum_density_x_comoving,
          state.momentum_density_y_comoving, state.momentum_density_z_comoving,
          state.total_energy_density_comoving, state.metal_mass_density_comoving, auxiliary};
}

// Exact self-contained boundary evolution: periodic/open/reflective ghost
// values are recomputed from this PRIVATE state at each substep. Imported MPI,
// same-level and coarse-fine ghosts are not self-contained and may NOT be
// reused as if they had a midpoint AMR temporal boundary history.
[[nodiscard]] bool closedPhysicalBoundary(const HydroPatchGeometry& geometry) {
  if (geometry.ghost_cells.empty()) return false;
  return std::all_of(geometry.ghost_cells.begin(), geometry.ghost_cells.end(),
      [&](const HydroGhostCell& ghost) {
        return ghost.mutation_rights == HydroGhostMutationRights::kWritablePhysicalBoundaryScratch &&
               (ghost.boundary_kind == HydroBoundaryKind::kPeriodic ||
                ghost.boundary_kind == HydroBoundaryKind::kReflective ||
                ghost.boundary_kind == HydroBoundaryKind::kOpen) &&
               ghost.owner_real_cell < geometry.cellCount() &&
               ghost.source_real_cell < geometry.cellCount() &&
               ghost.ghost_cell < geometry.totalCellStorageCount();
      });
}

[[nodiscard]] bool interiorCell(const HydroPatchGeometry& geometry, std::size_t row, int guard) {
  if (row >= geometry.cellCount()) return false;
  const auto ijk = geometry.cellIjk(row);
  const std::array<std::size_t, 3> dims{geometry.nx, geometry.ny, geometry.nz};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (dims[axis] <= 1U) continue;  // suppressed dimension
    const std::size_t g = static_cast<std::size_t>(guard);
    if (ijk[axis] < g || ijk[axis] + g >= dims[axis]) return false;
  }
  return true;
}

[[nodiscard]] std::array<double, 6> physicalLanes(const HydroConservedState& state) {
  const auto l = lanes(state, 0.0);
  return {l[0], l[1], l[2], l[3], l[4], l[5]};
}

[[nodiscard]] bool finiteState(const HydroConservedState& state) {
  for (double value : physicalLanes(state)) if (!std::isfinite(value)) return false;
  return state.mass_density_comoving > 0.0;
}

}  // namespace

bool hydroShadowCandidateLess(
    const HydroShadowCandidate& lhs, const HydroShadowCandidate& rhs,
    std::uint64_t step_index, std::uint64_t seed) noexcept {
  return std::tuple{candidateScore(lhs, step_index, seed), lhs.patch_id, lhs.level} <
         std::tuple{candidateScore(rhs, step_index, seed), rhs.patch_id, rhs.level};
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

HydroShadowReference::HydroShadowReference(HydroShadowOptions options, core::MemoryGovernor* governor)
    : m_options(options), m_governor(governor) {
  if (options.cadence == 0U) {
    throw std::invalid_argument("hydro shadow requires a nonzero cadence");
  }
}

void HydroShadowReference::beginEpoch(std::uint64_t step_index) {
  m_step_index = step_index;
  m_storage_cells_charged = 0U;
  m_workspace_bytes_charged = 0U;
  m_report = {};
  m_global_candidates.clear();
  m_selected_once.clear();
}

void HydroShadowReference::setGlobalCandidates(std::vector<HydroShadowCandidate> ranked) {
  m_global_candidates = rankHydroShadowCandidates(ranked, m_options.max_patches_per_rank,
                                                 m_step_index, m_options.seed);
  m_selected_once.clear();
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
  Sample sample;
  if (!m_options.enabled) return sample;
  ++m_report.eligible;
  if ((m_step_index % m_options.cadence) != 0U) {
    ++m_report.skipped_cadence;
    return sample;
  }
  // Global rank-independent selection, prepared once before conditional
  // patch loops. A selected physical patch may be sampled only once per epoch,
  // also when recursive AMR invokes the same patch in several subintervals.
  const HydroShadowCandidate identity{patch_id, level};
  if (std::none_of(m_global_candidates.begin(), m_global_candidates.end(),
                   [&](const auto& candidate) { return sameCandidate(candidate, identity); }) ||
      std::any_of(m_selected_once.begin(), m_selected_once.end(),
                  [&](const auto& candidate) { return sameCandidate(candidate, identity); })) {
    ++m_report.skipped_selection;
    return sample;
  }
  // Source terms are cloned and evaluated on the independent reference.
  // Discrete or unknown source types opt out instead of mutating live ledgers.
  if (std::any_of(source_terms.begin(), source_terms.end(),
                  [](const auto* source) { return source == nullptr; })) {
    ++m_report.skipped_unsupported_source;
    return sample;
  }
  // Composition sidecars are not part of HydroConservedStateSoa: do not
  // claim source consistency with an independently transported element/ion
  // population while production uses a separate face-flux sink.
  if (source_context.face_mass_flux_sink != nullptr) {
    ++m_report.skipped_external_passives;
    return sample;
  }
  const std::size_t count = geometry.cellCount();
  const bool temporal_refined =
      production_numerics.integration_scheme == HydroIntegrationScheme::kCornerTransportUpwind &&
      production_reconstruction == HydroReconstructionScheme::kPiecewiseParabolic;
  const bool closed_boundary = closedPhysicalBoundary(geometry);
  const int guard = closed_boundary ? 0 : (temporal_refined ? 7 : 4);
  if (count == 0U || before.size() != geometry.totalCellStorageCount() ||
      active.active_cells.size() != count ||
      active.active_faces.size() != geometry.faces.size() ||
      !std::any_of(active.active_cells.begin(), active.active_cells.end(),
                   [&](std::size_t row) { return interiorCell(geometry, row, guard); })) {
    if (count != 0U && !closed_boundary &&
        std::any_of(geometry.ghost_cells.begin(), geometry.ghost_cells.end(),
                    [](const HydroGhostCell& ghost) {
                      return ghost.mutation_rights == HydroGhostMutationRights::kReadOnlyImported ||
                             ghost.boundary_kind == HydroBoundaryKind::kImportedMpi;
                    })) {
      ++m_report.skipped_temporal_ghost;
    } else {
      ++m_report.skipped_geometry;
    }
    return sample;
  }
  if (before.size() > m_options.max_storage_cells_per_rank -
          std::min(m_storage_cells_charged, m_options.max_storage_cells_per_rank)) {
    ++m_report.skipped_resource;
    return sample;
  }
  // One private conserved SoA copy and one full CTU/PPM scratch. Double the
  // scratch budget for vector growth/allocator slack. This lease coexists
  // with the authoritative solver's already admitted scratch.
  const HydroActiveBatchPolicy batch{.max_active_cells = 4096U,
                                    .ctu_transverse_workspace = true,
                                    .dual_energy_workspace = production_numerics.dual_energy.mode != HydroDualEnergyMode::kOff};
  std::uint64_t requested = 0U;
  try {
    const std::uint64_t state_bytes = checkedMul(
        checkedMul(static_cast<std::uint64_t>(before.size()),
                   18U + (before.auxiliaryInternalEnergyEnabled() ? 2U : 0U)),
        sizeof(double));
    const std::uint64_t scratch_bytes = HydroScratchBuffers::workspaceBytes(
        count, active.active_faces.size(), count, before.size(), geometry.faceGhostSlotExtent(), batch);
    requested = checkedAdd(state_bytes, checkedMul(2U, scratch_bytes));
    requested = checkedAdd(requested, checkedMul(count, sizeof(HydroCellSourceProperties)));
    requested = checkedAdd(requested, checkedMul(count, sizeof(std::uint8_t)));
    requested = checkedAdd(requested, checkedMul(source_terms.size(),
                                               sizeof(std::unique_ptr<HydroSourceTerm>) + 512U));
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
    sample.initial = before;
    sample.shadow_sources.reserve(source_terms.size());
    for (const HydroSourceTerm* source : source_terms) {
      auto copy = source->cloneForShadow();
      if (!copy) {
        ++m_report.skipped_unsupported_source;
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
          source_context.effective_closure_allowed.begin() + count);
    }
    if (!source_terms.empty()) {
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
    sample.reconstruction_policy = production_policy;
    sample.riemann_kind = production_riemann;
    sample.closed_physical_boundary = closed_boundary;
    sample.source_active = !source_terms.empty() || sample.closure != nullptr;
    sample.evolving_comoving = update.comoving_coordinates &&
        (update.scale_factor_begin != update.scale_factor_end ||
         update.hubble_rate_begin_code != update.hubble_rate_end_code);
  } catch (const core::MemoryAdmissionError&) {
    ++m_report.skipped_resource;
    return sample;
  } catch (const std::bad_alloc&) {
    ++m_report.skipped_resource;
    return sample;
  } catch (const std::overflow_error&) {
    ++m_report.skipped_resource;
    return sample;
  } catch (const std::exception&) {
    // A malformed source snapshot/closure is not a memory refusal.
    // Production is separately responsible for its own state validation.
    ++m_report.skipped_numerical;
    return sample;
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

void HydroShadowReference::compare(
    Sample&& sample, const HydroConservedStateSoa& production_after,
    const HydroPatchGeometry& geometry, const HydroActiveSetView& active,
    const HydroUpdateContext& update, const HydroNumericsPolicy& production_numerics,
    HydroReconstructionScheme production_reconstruction) {
  if (!sample) return;
  // All sample work is private and non-authoritative, including Riemann
  // fallback counters; no production source terms, registers or flux sinks.
  try {
    // Save just the initial real interior for independent trajectory deltas.
    // This is covered by the second six-lane state allowance in capture().
    std::vector<HydroConservedState> baseline;
    baseline.reserve(geometry.cellCount());
    for (std::size_t row = 0; row < geometry.cellCount(); ++row) {
      baseline.push_back(sample.initial.loadCell(row));
    }
    HydroConservedStateSoa reference = sample.initial;
    HydroNumericsPolicy reference_policy = production_numerics;
    reference_policy.integration_scheme = HydroIntegrationScheme::kCornerTransportUpwind;
    HydroCoreSolver solver(sample.reconstruction_policy.adiabatic_index, reference_policy);
    auto riemann = makeHydroRiemannSolver(sample.riemann_kind);
    HydroScratchBuffers scratch;
    const HydroActiveBatchPolicy batch = solver.requiredWorkspaceFeatures(
        HydroActiveBatchPolicy{.max_active_cells = 4096U});
    (void)scratch.prepareForActiveSet(
        active.active_cells.size(), active.active_faces.size(), geometry.cellCount(),
        reference.size(), geometry.faceGhostSlotExtent(), batch);
    const bool temporal_refined =
        production_numerics.integration_scheme == HydroIntegrationScheme::kCornerTransportUpwind &&
        production_reconstruction == HydroReconstructionScheme::kPiecewiseParabolic;
    const std::uint32_t substeps = temporal_refined ? 2U : 1U;
    std::vector<const HydroSourceTerm*> shadow_sources;
    shadow_sources.reserve(sample.shadow_sources.size());
    for (const auto& source : sample.shadow_sources) shadow_sources.push_back(source.get());
    for (std::uint32_t substep = 0; substep < substeps; ++substep) {
      HydroUpdateContext sub_update = update;
      if (temporal_refined) {
        const double fraction_begin = static_cast<double>(substep) / 2.0;
        const double fraction_end = static_cast<double>(substep + 1U) / 2.0;
        sub_update.dt_code = update.dt_code / 2.0;
        sub_update.time_begin_code = std::lerp(update.time_begin_code, update.time_end_code, fraction_begin);
        sub_update.time_end_code = std::lerp(update.time_begin_code, update.time_end_code, fraction_end);
        sub_update.scale_factor_begin = std::lerp(update.scale_factor_begin, update.scale_factor_end, fraction_begin);
        sub_update.scale_factor_end = std::lerp(update.scale_factor_begin, update.scale_factor_end, fraction_end);
        sub_update.scale_factor = 0.5 * (sub_update.scale_factor_begin + sub_update.scale_factor_end);
        sub_update.hubble_rate_begin_code = std::lerp(update.hubble_rate_begin_code, update.hubble_rate_end_code, fraction_begin);
        sub_update.hubble_rate_end_code = std::lerp(update.hubble_rate_begin_code, update.hubble_rate_end_code, fraction_end);
        sub_update.hubble_rate_code = 0.5 * (sub_update.hubble_rate_begin_code + sub_update.hubble_rate_end_code);
      }
      HydroReconstructionPolicy reconstruction_policy = sample.reconstruction_policy;
      reconstruction_policy.dt_over_cell_width_code = hydroTransportCourantFactors(
          sub_update.dt_code, sub_update.scale_factor,
          geometry.cell_width_x_comoving, geometry.cell_width_y_comoving,
          geometry.cell_width_z_comoving);
      reconstruction_policy.dt_over_dx_code = reconstruction_policy.dt_over_cell_width_code[0];
      PpmReconstruction reconstruction(reconstruction_policy);
      HydroSourceContext sub_source{};
      sub_source.update = sub_update;
      sub_source.thermodynamic_closure = sample.closure;
      sub_source.effective_closure_allowed = sample.effective_closure_allowed;
      sub_source.redshift = sub_update.comoving_coordinates && sub_update.scale_factor > 0.0
          ? std::max(0.0, 1.0 / sub_update.scale_factor - 1.0) : sample.source_redshift;
      // This view is valid for the whole solver invocation; it recomputes
      // source physical properties from this substep's independent state.
      ShadowSourceProperties property_provider(
          reference, sample.initial, sample.captured_properties,
          update.scale_factor, sub_update.scale_factor, sub_update.comoving_coordinates);
      if (sample.source_properties_captured) {
        sub_source.source_property_provider = &property_provider;
      }
      if (sample.closed_physical_boundary) {
        fillHydroBoundaryGhostCells(reference, geometry, sample.reconstruction_policy.adiabatic_index);
        fillHydroBoundaryGhostAuxiliaryInternalEnergy(reference, geometry);
      }
      HydroProfileEvent profile;
      solver.advancePatchActiveSetBatchedWithScratch(
          reference, geometry, active, sub_update, reconstruction, *riemann, shadow_sources,
          sub_source, scratch, nullptr, batch, &profile, nullptr);
      m_report.fallback_count += profile.riemann_fallback_count +
          profile.positivity_fallback_count + profile.stencil_degradation_count;
    }
    // Account for BOTH independently resident copies and the baseline, plus
    // admitted source/closure snapshots. This capacity metric does not claim
    // to measure allocator metadata or opaque physics-table storage.
    std::uint64_t observed_capacity_bytes = checkedAdd(
        sample.initial.ownedCapacityBytes(), reference.ownedCapacityBytes());
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes, scratch.ownedCapacityBytes());
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes,
        checkedMul(baseline.capacity(), sizeof(HydroConservedState)));
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes,
        checkedMul(sample.captured_properties.capacity(), sizeof(HydroCellSourceProperties)));
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes,
        checkedMul(sample.effective_closure_allowed.capacity(), sizeof(std::uint8_t)));
    observed_capacity_bytes = checkedAdd(observed_capacity_bytes,
        checkedMul(sample.shadow_sources.capacity(), sizeof(std::unique_ptr<HydroSourceTerm>)));
    m_report.observed_capacity_high_water_bytes = std::max(
        m_report.observed_capacity_high_water_bytes, observed_capacity_bytes);
    // Restrict comparison to the domain unaffected by the copied ghost shell:
    // PPM/CTU stencil support is contained within a four-cell ring per full
    // step and a seven-cell ring for two halfsteps. On narrow AMR patches we
    // report an unsupported sample, not a spurious small error.
    const int guard = temporal_refined ? 7 : 4;
    std::uint64_t compared = 0U;
    std::uint64_t auxiliary_compared = 0U;
    const double volume = geometry.cell_volume_comoving;
    std::array<HydroShadowMetric, 7> metrics{};
    std::array<double, 6> reference_delta{};
    std::array<double, 6> production_delta{};
    for (std::size_t row : active.active_cells) {
      if (!sample.closed_physical_boundary && !interiorCell(geometry, row, guard)) continue;
      const auto produced = production_after.loadCell(row);
      const auto expected = reference.loadCell(row);
      if (!finiteState(produced) || !finiteState(expected)) {
        throw std::runtime_error("non-finite state in shadow hydro comparison");
      }
      const bool dual_energy = reference.auxiliaryInternalEnergyEnabled() &&
          production_after.auxiliaryInternalEnergyEnabled();
      const double aux_ref = dual_energy ? reference.auxiliaryInternalEnergyDensityComoving()[row] : 0.0;
      const double aux_prod = dual_energy ? production_after.auxiliaryInternalEnergyDensityComoving()[row] : 0.0;
      const auto left = lanes(produced, aux_prod);
      const auto right = lanes(expected, aux_ref);
      if (dual_energy) ++auxiliary_compared;
      for (std::size_t k = 0; k < (dual_energy ? 7U : 6U); ++k) {
        if (!std::isfinite(left[k]) || !std::isfinite(right[k])) {
          throw std::runtime_error("non-finite auxiliary state in shadow hydro comparison");
        }
        const double error = std::abs(left[k] - right[k]);
        metrics[k].weighted_abs_difference += error * volume;
        metrics[k].weighted_reference_magnitude += std::abs(right[k]) * volume;
        metrics[k].linf_difference = std::max(metrics[k].linf_difference, error);
      }
      // A local-cell conservation difference is not a global conservation
      // residual; report both trajectory deltas separately to prevent misuse.
      const auto baseline_lanes = physicalLanes(baseline[row]);
      const auto ref_lanes = physicalLanes(expected);
      const auto prod_lanes = physicalLanes(produced);
      for (std::size_t k = 0; k < 6U; ++k) {
        reference_delta[k] += (ref_lanes[k] - baseline_lanes[k]) * volume;
        production_delta[k] += (prod_lanes[k] - baseline_lanes[k]) * volume;
      }
      ++compared;
    }
    if (compared == 0U) {
      ++m_report.skipped_geometry;
      return;
    }
    for (std::size_t k = 0; k < 7U; ++k) {
      m_report.components[k].weighted_abs_difference += metrics[k].weighted_abs_difference;
      m_report.components[k].weighted_reference_magnitude += metrics[k].weighted_reference_magnitude;
      m_report.components[k].linf_difference = std::max(m_report.components[k].linf_difference,
                                                       metrics[k].linf_difference);
    }
    for (std::size_t k = 0; k < 6U; ++k) {
      m_report.reference_delta[k] += reference_delta[k];
      m_report.production_delta[k] += production_delta[k];
    }
    ++m_report.completed;
    if (sample.closed_physical_boundary) ++m_report.completed_closed_boundary_patches;
    if (sample.source_active) ++m_report.completed_source_active_patches;
    if (sample.evolving_comoving) ++m_report.completed_evolving_comoving_patches;
    m_report.compared_cells += compared;
    m_report.compared_auxiliary_cells += auxiliary_compared;
    m_report.compared_volume += volume * static_cast<double>(compared);
    m_report.reference_substeps += substeps;
  } catch (const std::exception&) {
    // Shadow is a diagnostic. A failed reference must not change a step,
    // acceptance decision, source history or physical model. Count and skip.
    ++m_report.skipped_numerical;
  }
}

} // namespace cosmosim::hydro
