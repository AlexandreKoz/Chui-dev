#include "cosmosim/hydro/hydro_shadow_reference.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
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

[[nodiscard]] std::array<double, 7> lanes(const HydroConservedState& state, double auxiliary) {
  return {state.mass_density_comoving, state.momentum_density_x_comoving,
          state.momentum_density_y_comoving, state.momentum_density_z_comoving,
          state.total_energy_density_comoving, state.metal_mass_density_comoving, auxiliary};
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

HydroShadowReference::HydroShadowReference(HydroShadowOptions options, core::MemoryGovernor* governor)
    : m_options(options), m_governor(governor) {
  if (options.cadence == 0U) {
    throw std::invalid_argument("hydro shadow requires a nonzero cadence and patch budget");
  }
}

void HydroShadowReference::beginEpoch(std::uint64_t step_index) {
  m_step_index = step_index;
  m_storage_cells_charged = 0U;
  m_workspace_bytes_charged = 0U;
  m_report = {};
}

HydroShadowReference::Sample HydroShadowReference::capture(
    std::uint64_t patch_id, int level,
    const HydroConservedStateSoa& before, const HydroPatchGeometry& geometry,
    const HydroActiveSetView& active, const HydroUpdateContext& update,
    const HydroNumericsPolicy& production_numerics,
    HydroReconstructionScheme production_reconstruction,
    std::size_t source_count, bool effective_closure_active) {
  Sample sample;
  if (!m_options.enabled) return sample;
  ++m_report.eligible;
  if ((m_step_index % m_options.cadence) != 0U) {
    ++m_report.skipped_cadence;
    return sample;
  }
  // Stable ID/level/step/seed selection, independent of vector iteration order.
  if ((mix64(patch_id ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(level)) << 32U) ^
             mix64(m_step_index ^ m_options.seed)) & 3ULL) != 0ULL ||
      m_report.selected >= m_options.max_patches_per_rank) {
    ++m_report.skipped_selection;
    return sample;
  }
  // Keep the same physical model: source operators, effective pressure and
  // evolving cosmological background require a separately time-centred shadow
  // source/ghost history and cannot be compared to a transport-only reference.
  if (source_count != 0U || effective_closure_active ||
      (update.comoving_coordinates &&
       (update.scale_factor_begin != update.scale_factor_end ||
        update.hubble_rate_begin_code != update.hubble_rate_end_code))) {
    ++m_report.skipped_sources_or_time;
    return sample;
  }
  const std::size_t count = geometry.cellCount();
  const bool temporal_refined =
      production_numerics.integration_scheme == HydroIntegrationScheme::kCornerTransportUpwind &&
      production_reconstruction == HydroReconstructionScheme::kPiecewiseParabolic;
  const int guard = temporal_refined ? 7 : 4; // two PPM/CTU halfsteps vs one full PPM/CTU step
  if (count == 0U || before.size() != geometry.totalCellStorageCount() ||
      active.active_cells.size() != count ||
      active.active_faces.size() != geometry.faces.size() ||
      !std::any_of(active.active_cells.begin(), active.active_cells.end(),
                   [&](std::size_t row) { return interiorCell(geometry, row, guard); })) {
    ++m_report.skipped_geometry;
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
    const std::uint64_t state_bytes = static_cast<std::uint64_t>(before.size()) *
        (12U + (before.auxiliaryInternalEnergyEnabled() ? 1U : 0U)) * sizeof(double);
    const std::uint64_t scratch_bytes = HydroScratchBuffers::workspaceBytes(
        count, active.active_faces.size(), count, before.size(), geometry.faceGhostSlotExtent(), batch);
    if (scratch_bytes > (std::numeric_limits<std::uint64_t>::max() - state_bytes) / 2U) {
      ++m_report.skipped_resource;
      return sample;
    }
    requested = state_bytes + 2U * scratch_bytes;
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
  } catch (const core::MemoryAdmissionError&) {
    ++m_report.skipped_resource;
    return sample;
  } catch (const std::bad_alloc&) {
    ++m_report.skipped_resource;
    return sample;
  } catch (const std::exception&) {
    // Optional diagnostic admission must not cause a physical step rejection.
    ++m_report.skipped_resource;
    return sample;
  }
  sample.admitted_bytes = requested;
  sample.valid = true;
  m_storage_cells_charged += before.size();
  m_workspace_bytes_charged += requested;
  ++m_report.selected;
  m_report.workspace_high_water_bytes = std::max(m_report.workspace_high_water_bytes, requested);
  return sample;
}

void HydroShadowReference::compare(
    Sample&& sample, const HydroConservedStateSoa& production_after,
    const HydroPatchGeometry& geometry, const HydroActiveSetView& active,
    const HydroUpdateContext& update, const HydroNumericsPolicy& production_numerics,
    HydroReconstructionScheme production_reconstruction,
    double adiabatic_index, double density_floor, double pressure_floor) {
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
    HydroConservedStateSoa reference = std::move(sample.initial);
    HydroNumericsPolicy reference_policy = production_numerics;
    reference_policy.integration_scheme = HydroIntegrationScheme::kCornerTransportUpwind;
    HydroCoreSolver solver(adiabatic_index, reference_policy);
    HllcRiemannSolver riemann;
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
    HydroSourceContext source_context{};
    source_context.update = update;
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
      HydroReconstructionPolicy reconstruction_policy{};
      reconstruction_policy.adiabatic_index = adiabatic_index;
      reconstruction_policy.rho_floor = std::max(density_floor, 1.0e-30);
      reconstruction_policy.pressure_floor = std::max(pressure_floor, 1.0e-30);
      reconstruction_policy.dt_over_cell_width_code = hydroTransportCourantFactors(
          sub_update.dt_code, sub_update.scale_factor,
          geometry.cell_width_x_comoving, geometry.cell_width_y_comoving,
          geometry.cell_width_z_comoving);
      reconstruction_policy.dt_over_dx_code = reconstruction_policy.dt_over_cell_width_code[0];
      PpmReconstruction reconstruction(reconstruction_policy);
      source_context.update = sub_update;
      HydroProfileEvent profile;
      solver.advancePatchActiveSetBatchedWithScratch(
          reference, geometry, active, sub_update, reconstruction, riemann, {},
          source_context, scratch, nullptr, batch, &profile, nullptr);
      m_report.fallback_count += profile.riemann_fallback_count +
          profile.positivity_fallback_count + profile.stencil_degradation_count;
    }
    m_report.observed_capacity_high_water_bytes = std::max(
        m_report.observed_capacity_high_water_bytes,
        reference.ownedCapacityBytes() + scratch.ownedCapacityBytes() +
            static_cast<std::uint64_t>(baseline.capacity()) * sizeof(HydroConservedState));
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
      if (!interiorCell(geometry, row, guard)) continue;
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
