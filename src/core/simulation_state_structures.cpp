#include "cosmosim/core/simulation_state.hpp"

#include <cmath>
#include <limits>

namespace cosmosim::core {

void ParticleSoa::resize(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalParticleCount, "particle", "ParticleSoa::resize");
  position_x_comoving.resize(count);
  position_y_comoving.resize(count);
  position_z_comoving.resize(count);
  velocity_x_peculiar.resize(count);
  velocity_y_peculiar.resize(count);
  velocity_z_peculiar.resize(count);
  mass_code.resize(count);
  time_bin.resize(count);
}

std::size_t ParticleSoa::size() const noexcept { return position_x_comoving.size(); }

bool ParticleSoa::isConsistent() const noexcept {
  const std::size_t expected = position_x_comoving.size();
  if (position_y_comoving.size() != expected || position_z_comoving.size() != expected ||
      velocity_x_peculiar.size() != expected || velocity_y_peculiar.size() != expected ||
      velocity_z_peculiar.size() != expected || mass_code.size() != expected ||
      (!time_bin.empty() && time_bin.size() != expected)) {
    return false;
  }
  for (std::size_t i = 0; i < expected; ++i) {
    if (!std::isfinite(position_x_comoving[i]) || !std::isfinite(position_y_comoving[i]) ||
        !std::isfinite(position_z_comoving[i]) || !std::isfinite(velocity_x_peculiar[i]) ||
        !std::isfinite(velocity_y_peculiar[i]) || !std::isfinite(velocity_z_peculiar[i]) ||
        !std::isfinite(mass_code[i]) || mass_code[i] < 0.0) {
      return false;
    }
  }
  return true;
}

void ParticleSidecar::resize(std::size_t count) {
  particle_id.resize(count);
  sfc_key.resize(count);
  species_tag.resize(count);
  particle_flags.resize(count);
  owning_rank.resize(count);
  const std::size_t old_drift_size = last_drift_time_code.size();
  last_drift_time_code.resize(count, 0.0);
  last_drift_scale_factor.resize(count, 1.0);
  for (std::size_t i = old_drift_size; i < count; ++i) {
    last_drift_scale_factor[i] = 1.0;
  }
  if (!gravity_softening_comoving.empty()) {
    gravity_softening_comoving.resize(count);
  }
  if (!has_gravity_softening_override.empty()) {
    has_gravity_softening_override.resize(count, 0U);
    if (gravity_softening_comoving.empty()) {
      gravity_softening_comoving.resize(count, 0.0);
    }
  }
}

std::size_t ParticleSidecar::size() const noexcept { return particle_id.size(); }

bool ParticleSidecar::isConsistent() const noexcept {
  const std::size_t expected = particle_id.size();
  if ((!sfc_key.empty() && sfc_key.size() != expected) ||
      (!species_tag.empty() && species_tag.size() != expected) ||
      (!particle_flags.empty() && particle_flags.size() != expected) ||
      (!owning_rank.empty() && owning_rank.size() != expected) ||
      (!last_drift_time_code.empty() && last_drift_time_code.size() != expected) ||
      (!last_drift_scale_factor.empty() && last_drift_scale_factor.size() != expected) ||
      (!gravity_softening_comoving.empty() && gravity_softening_comoving.size() != expected) ||
      (!has_gravity_softening_override.empty() && has_gravity_softening_override.size() != expected) ||
      (!has_gravity_softening_override.empty() && gravity_softening_comoving.empty())) {
    return false;
  }
  if (last_drift_time_code.size() != last_drift_scale_factor.size()) return false;
  for (std::size_t i = 0; i < last_drift_time_code.size(); ++i) {
    if (!std::isfinite(last_drift_time_code[i]) || !std::isfinite(last_drift_scale_factor[i]) ||
        last_drift_scale_factor[i] <= 0.0) {
      return false;
    }
  }
  for (std::size_t i = 0; i < gravity_softening_comoving.size(); ++i) {
    const bool has_softening_override =
        !has_gravity_softening_override.empty() && has_gravity_softening_override[i] != 0U;
    if (has_softening_override &&
        (!std::isfinite(gravity_softening_comoving[i]) || !(gravity_softening_comoving[i] > 0.0))) {
      return false;
    }
  }
  for (const auto flag : has_gravity_softening_override) {
    if (flag > 1U) {
      return false;
    }
  }
  return true;
}

bool ParticleSidecar::hasGravitySofteningOverride(std::size_t particle_index) const {
  if (particle_index >= particle_id.size()) {
    throw std::out_of_range("ParticleSidecar.hasGravitySofteningOverride: particle index out of range");
  }
  return !has_gravity_softening_override.empty() && has_gravity_softening_override[particle_index] != 0U;
}

double ParticleSidecar::gravitySofteningOverride(std::size_t particle_index) const {
  if (!hasGravitySofteningOverride(particle_index)) {
    throw std::runtime_error("ParticleSidecar.gravitySofteningOverride: particle has no override");
  }
  if (gravity_softening_comoving.empty() || particle_index >= gravity_softening_comoving.size()) {
    throw std::runtime_error("ParticleSidecar.gravitySofteningOverride: override value sidecar missing");
  }
  return gravity_softening_comoving[particle_index];
}

void ParticleSidecar::setGravitySofteningOverride(std::size_t particle_index, double epsilon_comoving) {
  if (particle_index >= particle_id.size()) {
    throw std::out_of_range("ParticleSidecar.setGravitySofteningOverride: particle index out of range");
  }
  if (gravity_softening_comoving.empty()) {
    gravity_softening_comoving.resize(particle_id.size(), 0.0);
  }
  if (has_gravity_softening_override.empty()) {
    has_gravity_softening_override.resize(particle_id.size(), 0U);
  }
  gravity_softening_comoving[particle_index] = epsilon_comoving;
  has_gravity_softening_override[particle_index] = 1U;
}

void ParticleSidecar::clearGravitySofteningOverride(std::size_t particle_index) {
  if (particle_index >= particle_id.size()) {
    throw std::out_of_range("ParticleSidecar.clearGravitySofteningOverride: particle index out of range");
  }
  if (!has_gravity_softening_override.empty()) {
    has_gravity_softening_override[particle_index] = 0U;
  }
}

void CellSoa::resize(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalCellCount, "gas-cell", "CellSoa::resize");
  center_x_comoving.resize(count);
  center_y_comoving.resize(count);
  center_z_comoving.resize(count);
  mass_code.resize(count);
  time_bin.resize(count);
  patch_index.resize(count);
}

std::size_t CellSoa::size() const noexcept { return center_x_comoving.size(); }

bool CellSoa::isConsistent() const noexcept {
  const std::size_t expected = center_x_comoving.size();
  if (center_y_comoving.size() != expected || center_z_comoving.size() != expected ||
      mass_code.size() != expected || time_bin.size() != expected || patch_index.size() != expected) {
    return false;
  }
  for (std::size_t i = 0; i < expected; ++i) {
    if (!std::isfinite(center_x_comoving[i]) || !std::isfinite(center_y_comoving[i]) ||
        !std::isfinite(center_z_comoving[i]) || !std::isfinite(mass_code[i]) || mass_code[i] < 0.0) {
      return false;
    }
  }
  return true;
}

void GasCellSidecar::resize(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalCellCount, "gas-cell", "GasCellSidecar::resize");
  gas_cell_id.resize(count);
  parent_particle_id.resize(count);
  velocity_x_peculiar.resize(count);
  velocity_y_peculiar.resize(count);
  velocity_z_peculiar.resize(count);
  density_code.resize(count);
  pressure_code.resize(count);
  internal_energy_code.resize(count);
  metal_mass_code.resize(count);
  temperature_code.resize(count);
  sound_speed_code.resize(count);
}

std::size_t GasCellSidecar::size() const noexcept { return density_code.size(); }

bool GasCellSidecar::isConsistent() const noexcept {
  const std::size_t expected = density_code.size();
  if (gas_cell_id.size() != expected || parent_particle_id.size() != expected ||
      velocity_x_peculiar.size() != expected || velocity_y_peculiar.size() != expected ||
      velocity_z_peculiar.size() != expected || pressure_code.size() != expected ||
      internal_energy_code.size() != expected || metal_mass_code.size() != expected ||
      temperature_code.size() != expected || sound_speed_code.size() != expected) {
    return false;
  }
  for (std::size_t i = 0; i < expected; ++i) {
    if (!std::isfinite(velocity_x_peculiar[i]) || !std::isfinite(velocity_y_peculiar[i]) ||
        !std::isfinite(velocity_z_peculiar[i]) || !std::isfinite(density_code[i]) || density_code[i] < 0.0 ||
        !std::isfinite(pressure_code[i]) || pressure_code[i] < 0.0 ||
        !std::isfinite(internal_energy_code[i]) || internal_energy_code[i] < 0.0 ||
        !std::isfinite(metal_mass_code[i]) || metal_mass_code[i] < 0.0 ||
        !std::isfinite(temperature_code[i]) || temperature_code[i] < 0.0 ||
        !std::isfinite(sound_speed_code[i]) || sound_speed_code[i] < 0.0) {
      return false;
    }
  }
  return true;
}

void StarParticleSidecar::resize(std::size_t count) {
  particle_index.resize(count);
  formation_scale_factor.resize(count);
  birth_mass_code.resize(count);
  metallicity_mass_fraction.resize(count);
  birth_key.resize(count);
  parent_gas_cell_id.resize(count);
  birth_tick.resize(count);
  birth_ordinal.resize(count);
  stellar_age_years_last.resize(count);
  stellar_returned_mass_cumulative_code.resize(count);
  stellar_returned_metals_cumulative_code.resize(count);
  stellar_newly_synthesized_metals_cumulative_code.resize(count);
  stellar_feedback_energy_cumulative_erg.resize(count);
  enrichment_carry_mass_code.resize(count);
  enrichment_carry_metals_code.resize(count);
  enrichment_carry_feedback_energy_erg.resize(count);
  enrichment_carry_momentum_code.resize(count);
  stellar_deposited_mass_cumulative_code.resize(count);
  stellar_deposited_metals_cumulative_code.resize(count);
  stellar_deposited_feedback_energy_cumulative_erg.resize(count);
  for (std::size_t channel = 0; channel < stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    stellar_returned_mass_channel_cumulative_code[channel].resize(count);
    stellar_returned_metals_channel_cumulative_code[channel].resize(count);
    stellar_feedback_energy_channel_cumulative_erg[channel].resize(count);
  }
}

std::size_t StarParticleSidecar::size() const noexcept { return particle_index.size(); }

bool StarParticleSidecar::isConsistent() const noexcept {
  const std::size_t expected = particle_index.size();
  if (formation_scale_factor.size() != expected || birth_mass_code.size() != expected ||
      metallicity_mass_fraction.size() != expected || birth_key.size() != expected ||
      parent_gas_cell_id.size() != expected || birth_tick.size() != expected ||
      birth_ordinal.size() != expected || stellar_age_years_last.size() != expected ||
      stellar_returned_mass_cumulative_code.size() != expected ||
      stellar_returned_metals_cumulative_code.size() != expected ||
      stellar_newly_synthesized_metals_cumulative_code.size() != expected ||
      stellar_feedback_energy_cumulative_erg.size() != expected ||
      enrichment_carry_mass_code.size() != expected ||
      enrichment_carry_metals_code.size() != expected ||
      enrichment_carry_feedback_energy_erg.size() != expected ||
      enrichment_carry_momentum_code.size() != expected ||
      stellar_deposited_mass_cumulative_code.size() != expected ||
      stellar_deposited_metals_cumulative_code.size() != expected ||
      stellar_deposited_feedback_energy_cumulative_erg.size() != expected) {
    return false;
  }

  for (std::size_t channel = 0; channel < stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    if (stellar_returned_mass_channel_cumulative_code[channel].size() != expected ||
        stellar_returned_metals_channel_cumulative_code[channel].size() != expected ||
        stellar_feedback_energy_channel_cumulative_erg[channel].size() != expected) {
      return false;
    }
  }
  for (std::size_t i = 0; i < expected; ++i) {
    if (!std::isfinite(formation_scale_factor[i]) || formation_scale_factor[i] < 0.0 ||
        !std::isfinite(birth_mass_code[i]) || birth_mass_code[i] < 0.0 ||
        !std::isfinite(metallicity_mass_fraction[i]) || metallicity_mass_fraction[i] < 0.0 ||
        metallicity_mass_fraction[i] > 1.0 || !std::isfinite(stellar_age_years_last[i]) ||
        stellar_age_years_last[i] < 0.0 || !std::isfinite(stellar_returned_mass_cumulative_code[i]) ||
        stellar_returned_mass_cumulative_code[i] < 0.0 ||
        !std::isfinite(stellar_returned_metals_cumulative_code[i]) ||
        stellar_returned_metals_cumulative_code[i] < 0.0 ||
        !std::isfinite(stellar_feedback_energy_cumulative_erg[i]) ||
        stellar_feedback_energy_cumulative_erg[i] < 0.0) {
      return false;
    }
  }
  return true;
}

void BlackHoleParticleSidecar::resize(std::size_t count) {
  particle_index.resize(count);
  host_cell_index.resize(count, kInvalidGasCellRow);
  subgrid_mass_code.resize(count);
  accretion_rate_code.resize(count);
  feedback_energy_code.resize(count);
  eddington_ratio.resize(count);
  cumulative_accreted_mass_code.resize(count);
  cumulative_feedback_energy_code.resize(count);
  duty_cycle_active_time_code.resize(count);
  duty_cycle_total_time_code.resize(count);
}

std::size_t BlackHoleParticleSidecar::size() const noexcept { return particle_index.size(); }

bool BlackHoleParticleSidecar::isConsistent() const noexcept {
  const std::size_t expected = particle_index.size();
  if (host_cell_index.size() != expected || subgrid_mass_code.size() != expected ||
      accretion_rate_code.size() != expected || feedback_energy_code.size() != expected ||
      eddington_ratio.size() != expected || cumulative_accreted_mass_code.size() != expected ||
      cumulative_feedback_energy_code.size() != expected || duty_cycle_active_time_code.size() != expected ||
      duty_cycle_total_time_code.size() != expected) {
    return false;
  }
  for (std::size_t i = 0; i < expected; ++i) {
    if (!std::isfinite(subgrid_mass_code[i]) || subgrid_mass_code[i] < 0.0 ||
        !std::isfinite(accretion_rate_code[i]) || accretion_rate_code[i] < 0.0 ||
        !std::isfinite(feedback_energy_code[i]) || feedback_energy_code[i] < 0.0 ||
        !std::isfinite(eddington_ratio[i]) || eddington_ratio[i] < 0.0 ||
        !std::isfinite(cumulative_accreted_mass_code[i]) || cumulative_accreted_mass_code[i] < 0.0 ||
        !std::isfinite(cumulative_feedback_energy_code[i]) || cumulative_feedback_energy_code[i] < 0.0 ||
        !std::isfinite(duty_cycle_active_time_code[i]) || duty_cycle_active_time_code[i] < 0.0 ||
        !std::isfinite(duty_cycle_total_time_code[i]) || duty_cycle_total_time_code[i] < 0.0) {
      return false;
    }
  }
  return true;
}

void TracerParticleSidecar::resize(std::size_t count) {
  particle_index.resize(count);
  parent_particle_id.resize(count);
  injection_step.resize(count);
  host_cell_index.resize(count, kInvalidGasCellRow);
  mass_fraction_of_host.resize(count);
  last_host_mass_code.resize(count);
  cumulative_exchanged_mass_code.resize(count);
}

std::size_t TracerParticleSidecar::size() const noexcept { return particle_index.size(); }

bool TracerParticleSidecar::isConsistent() const noexcept {
  const std::size_t expected = particle_index.size();
  if (parent_particle_id.size() != expected || injection_step.size() != expected ||
      host_cell_index.size() != expected || mass_fraction_of_host.size() != expected ||
      last_host_mass_code.size() != expected || cumulative_exchanged_mass_code.size() != expected) {
    return false;
  }
  for (std::size_t i = 0; i < expected; ++i) {
    if (!std::isfinite(mass_fraction_of_host[i]) || mass_fraction_of_host[i] < 0.0 ||
        mass_fraction_of_host[i] > 1.0 || !std::isfinite(last_host_mass_code[i]) ||
        last_host_mass_code[i] < 0.0 || !std::isfinite(cumulative_exchanged_mass_code[i]) ||
        cumulative_exchanged_mass_code[i] < 0.0) {
      return false;
    }
  }
  return true;
}

bool MagneticFaceSoa::isConsistent() const noexcept {
  if (patches.empty()) {
    return face_flux_x_comoving.empty() && face_flux_y_comoving.empty() &&
        face_flux_z_comoving.empty();
  }
  // Offsets are concatenation cursors; each record's slice must lie entirely
  // inside its lane and the records must tile the lanes in order.
  std::uint64_t expected_x = 0U;
  std::uint64_t expected_y = 0U;
  std::uint64_t expected_z = 0U;
  for (const MagneticFacePatchRecord& patch : patches) {
    if (patch.nx == 0U || patch.ny == 0U || patch.nz == 0U) {
      return false;
    }
    if (!(patch.cell_width_x_comoving > 0.0) || !(patch.cell_width_y_comoving > 0.0) ||
        !(patch.cell_width_z_comoving > 0.0)) {
      return false;
    }
    if (!std::isfinite(patch.origin_x_comoving) || !std::isfinite(patch.origin_y_comoving) ||
        !std::isfinite(patch.origin_z_comoving)) {
      return false;
    }
    // nx, ny, nz are uint32, so (n+1) fits in uint64 and every pairwise
    // product of two factors < 2^33 fits; the triple product is checked.
    const std::uint64_t nx = patch.nx;
    const std::uint64_t ny = patch.ny;
    const std::uint64_t nz = patch.nz;
    const auto checked_triple = [](std::uint64_t a, std::uint64_t b, std::uint64_t c,
                                   std::uint64_t& out) {
      const std::uint64_t ab = a * b;
      if (c != 0U && ab > std::numeric_limits<std::uint64_t>::max() / c) {
        return false;
      }
      out = ab * c;
      return true;
    };
    std::uint64_t count_x = 0U;
    std::uint64_t count_y = 0U;
    std::uint64_t count_z = 0U;
    if (!checked_triple(nx + 1U, ny, nz, count_x) || !checked_triple(nx, ny + 1U, nz, count_y) ||
        !checked_triple(nx, ny, nz + 1U, count_z)) {
      return false;
    }
    if (!std::isfinite(patch.cell_width_x_comoving) || !std::isfinite(patch.cell_width_y_comoving) ||
        !std::isfinite(patch.cell_width_z_comoving) || patch.level < 0) {
      return false;
    }
    if (patch.flux_x_offset != expected_x || patch.flux_y_offset != expected_y ||
        patch.flux_z_offset != expected_z) {
      return false;
    }
    if (count_x > std::numeric_limits<std::uint64_t>::max() - expected_x ||
        count_y > std::numeric_limits<std::uint64_t>::max() - expected_y ||
        count_z > std::numeric_limits<std::uint64_t>::max() - expected_z) {
      return false;
    }
    expected_x += count_x;
    expected_y += count_y;
    expected_z += count_z;
  }
  return face_flux_x_comoving.size() == expected_x &&
      face_flux_y_comoving.size() == expected_y &&
      face_flux_z_comoving.size() == expected_z;
}

void PatchSoa::resize(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalPatchCount, "AMR patch", "PatchSoa::resize");
  patch_id.resize(count);
  level.resize(count);
  first_cell.resize(count);
  cell_count.resize(count);
  parent_patch_id.resize(count);
  morton_key.resize(count);
  origin_x_comoving.resize(count);
  origin_y_comoving.resize(count);
  origin_z_comoving.resize(count);
  extent_x_comoving.resize(count);
  extent_y_comoving.resize(count);
  extent_z_comoving.resize(count);
  cell_dim_x.resize(count);
  cell_dim_y.resize(count);
  cell_dim_z.resize(count);
  owning_rank.resize(count);
}

std::size_t PatchSoa::size() const noexcept { return patch_id.size(); }

bool PatchSoa::isConsistent() const noexcept {
  const std::size_t expected = patch_id.size();
  if (level.size() != expected || first_cell.size() != expected || cell_count.size() != expected ||
      parent_patch_id.size() != expected || morton_key.size() != expected ||
      origin_x_comoving.size() != expected || origin_y_comoving.size() != expected ||
      origin_z_comoving.size() != expected || extent_x_comoving.size() != expected ||
      extent_y_comoving.size() != expected || extent_z_comoving.size() != expected ||
      cell_dim_x.size() != expected || cell_dim_y.size() != expected || cell_dim_z.size() != expected ||
      owning_rank.size() != expected) {
    return false;
  }
  for (std::size_t i = 0; i < expected; ++i) {
    const bool has_any_geometry_value = extent_x_comoving[i] != 0.0 ||
        extent_y_comoving[i] != 0.0 || extent_z_comoving[i] != 0.0 ||
        cell_dim_x[i] != 0U || cell_dim_y[i] != 0U || cell_dim_z[i] != 0U;
    if (cell_count[i] == 0U) {
      continue;
    }
    // A fully zero geometry tuple is the documented fixed-patch compatibility
    // sentinel. It is structurally valid; workflow code must derive and verify
    // the physical Cartesian layout from cell centers before using it. Partial
    // or non-positive geometry, by contrast, is never a valid PatchSoa
    // descriptor. Even the sentinel must carry finite origins.
    if (!has_any_geometry_value) {
      if (!std::isfinite(origin_x_comoving[i]) || !std::isfinite(origin_y_comoving[i]) ||
          !std::isfinite(origin_z_comoving[i])) {
        return false;
      }
      continue;
    }
    if (!(extent_x_comoving[i] > 0.0 && extent_y_comoving[i] > 0.0 && extent_z_comoving[i] > 0.0) ||
        cell_dim_x[i] == 0U || cell_dim_y[i] == 0U || cell_dim_z[i] == 0U) {
      return false;
    }
    const std::size_t dim_x = static_cast<std::size_t>(cell_dim_x[i]);
    const std::size_t dim_y = static_cast<std::size_t>(cell_dim_y[i]);
    const std::size_t dim_z = static_cast<std::size_t>(cell_dim_z[i]);
    if (dim_x > std::numeric_limits<std::size_t>::max() / dim_y ||
        dim_x * dim_y > std::numeric_limits<std::size_t>::max() / dim_z) {
      return false;
    }
    const std::size_t dims_product = dim_x * dim_y * dim_z;
    if (dims_product != cell_count[i]) {
      return false;
    }
    if (!std::isfinite(origin_x_comoving[i]) || !std::isfinite(origin_y_comoving[i]) ||
        !std::isfinite(origin_z_comoving[i]) || !std::isfinite(extent_x_comoving[i]) ||
        !std::isfinite(extent_y_comoving[i]) || !std::isfinite(extent_z_comoving[i])) {
      return false;
    }
  }
  return true;
}

}  // namespace cosmosim::core
