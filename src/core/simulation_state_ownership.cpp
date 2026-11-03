#include "cosmosim/core/simulation_state.hpp"

#include <algorithm>
#include <limits>
#include "cosmosim/core/memory_governor.hpp"
#include <stdexcept>

namespace cosmosim::core {

std::uint64_t OwnershipValidationWorkspace::requiredBytes(
    std::size_t particle_count, std::size_t cell_count,
    bool need_star_rows, bool need_bh_rows, bool need_tracer_rows,
    bool need_cell_owner, bool include_id_validation_workspace) {
  const auto mul = [](std::uint64_t count, std::uint64_t width) {
    if (count != 0U && width > std::numeric_limits<std::uint64_t>::max() / count) {
      throw std::overflow_error("ownership validation scratch byte overflow");
    }
    return count * width;
  };
  std::uint64_t bytes = include_id_validation_workspace
      ? k_id_validation_workspace_limit_bytes : 0U;
  const std::uint64_t particles_u64 = static_cast<std::uint64_t>(particle_count);
  if (need_star_rows) {
    bytes = checkedMemoryBytesAdd(bytes, mul(particles_u64, sizeof(std::uint8_t)),
                                  "ownership validation star markers");
  }
  if (need_bh_rows) {
    bytes = checkedMemoryBytesAdd(bytes, mul(particles_u64, sizeof(std::uint8_t)),
                                  "ownership validation BH markers");
  }
  if (need_tracer_rows) {
    bytes = checkedMemoryBytesAdd(bytes, mul(particles_u64, sizeof(std::uint8_t)),
                                  "ownership validation tracer markers");
  }
  if (need_cell_owner) {
    bytes = checkedMemoryBytesAdd(
        bytes, mul(static_cast<std::uint64_t>(cell_count), sizeof(std::uint32_t)),
        "ownership validation cell ownership");
  }
  return bytes;
}

void OwnershipValidationWorkspace::resize(
    std::size_t particle_count, std::size_t cell_count,
    bool need_star_rows, bool need_bh_rows, bool need_tracer_rows,
    bool need_cell_owner) {
  (void)requiredBytes(
      particle_count, cell_count, need_star_rows, need_bh_rows,
      need_tracer_rows, need_cell_owner, false);
  const auto resize_marker = [particle_count](std::vector<std::uint8_t>& lane, bool needed) {
    if (needed) {
      lane.assign(particle_count, 0U);
    } else {
      std::vector<std::uint8_t>{}.swap(lane);
    }
  };
  resize_marker(star_rows, need_star_rows);
  resize_marker(bh_rows, need_bh_rows);
  resize_marker(tracer_rows, need_tracer_rows);
  if (need_cell_owner) {
    cell_owner.resize(cell_count);
  } else {
    std::vector<std::uint32_t>{}.swap(cell_owner);
  }
}

std::uint64_t OwnershipValidationWorkspace::ownedCapacityBytes() const {
  std::uint64_t bytes = 0U;
  const auto add = [&bytes](std::uint64_t count, std::uint64_t width) {
    if (count != 0U && width > std::numeric_limits<std::uint64_t>::max() / count) {
      throw std::overflow_error("ownership validation retained capacity overflow");
    }
    bytes = checkedMemoryBytesAdd(bytes, count * width,
                                  "ownership validation retained capacity");
  };
  add(static_cast<std::uint64_t>(id_bucket.capacity()), sizeof(std::uint64_t));
  add(static_cast<std::uint64_t>(star_rows.capacity()), sizeof(std::uint8_t));
  add(static_cast<std::uint64_t>(bh_rows.capacity()), sizeof(std::uint8_t));
  add(static_cast<std::uint64_t>(tracer_rows.capacity()), sizeof(std::uint8_t));
  add(static_cast<std::uint64_t>(cell_owner.capacity()), sizeof(std::uint32_t));
  return bytes;
}

void SimulationState::resizeParticles(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalParticleCount, "particle", "SimulationState::resizeParticles");
  if (hasHomogeneousDmoMetadata()) {
    if (count != particles.size()) {
      throw std::logic_error(
          "SimulationState::resizeParticles requires explicit governed metadata materialization before changing a compact-DMO population");
    }
    return;
  }
  particles.resize(count);
  particle_sidecar.resize(count);
  bumpParticleIdentityGeneration();
  bumpParticleIndexGeneration();
}

void SimulationState::resizeCells(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalCellCount, "gas-cell", "SimulationState::resizeCells");
  cells.resize(count);
  gas_cells.resize(count);
  gas_cell_identity.clear();
  bumpCellIndexGeneration();
}

void SimulationState::resizePatches(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalPatchCount, "AMR patch", "SimulationState::resizePatches");
  patches.resize(count);
}

bool SimulationState::validateOwnershipInvariants() const {
  OwnershipValidationWorkspace scratch;
  return validateOwnershipInvariantsImpl(scratch, false);
}

bool SimulationState::validateOwnershipInvariants(
    OwnershipValidationWorkspace& scratch) const {
  return validateOwnershipInvariantsImpl(scratch, true);
}

bool SimulationState::validateOwnershipInvariantsImpl(
    OwnershipValidationWorkspace& scratch, bool use_bounded_id_scratch) const {
  if (!particles.isConsistent() || !particle_sidecar.isConsistent() || !cells.isConsistent() ||
      !gas_cells.isConsistent() || !patches.isConsistent() || !star_particles.isConsistent() ||
      !black_holes.isConsistent() || !tracers.isConsistent()) {
    return false;
  }

  if (particles.size() != particle_sidecar.size()) {
    return false;
  }
  if (hasHomogeneousDmoMetadata()) {
    if ((!particles.time_bin.empty() && particles.time_bin.size() != particles.size()) || !particle_sidecar.species_tag.empty() ||
        !particle_sidecar.particle_flags.empty() || !particle_sidecar.owning_rank.empty() ||
        !particle_sidecar.last_drift_time_code.empty() || !particle_sidecar.last_drift_scale_factor.empty()) {
      return false;
    }
    if (species.count_by_species[particleSpeciesIndex(ParticleSpecies::kDarkMatter)] != particles.size()) return false;
    for (std::size_t i = 0; i < k_particle_species_count; ++i) {
      if (i != particleSpeciesIndex(ParticleSpecies::kDarkMatter) && species.count_by_species[i] != 0U) return false;
    }
  } else {
    const std::size_t n = particles.size();
    if (particles.time_bin.size() != n || particle_sidecar.species_tag.size() != n ||
        particle_sidecar.particle_flags.size() != n || particle_sidecar.owning_rank.size() != n ||
        particle_sidecar.last_drift_time_code.size() != n || particle_sidecar.last_drift_scale_factor.size() != n) {
      return false;
    }
  }

  const bool need_star_rows = star_particles.size() != 0U;
  const bool need_bh_rows = black_holes.size() != 0U;
  const bool need_tracer_rows = tracers.size() != 0U;
  const bool need_cell_owner = patches.size() != 0U && cells.size() != 0U;
  scratch.resize(
      particles.size(), cells.size(), need_star_rows, need_bh_rows,
      need_tracer_rows, need_cell_owner);
  if (!(use_bounded_id_scratch
            ? validateUniqueParticleIds(scratch) : validateUniqueParticleIds())) {
    return false;
  }

  if (cells.size() != gas_cells.size()) {
    return false;
  }

  if (!gas_cell_identity.isConsistent() || !gas_cell_identity.coversDenseLocalRows(cells.size()) ||
      !gasCellIdentityMapMatchesSidecarLanes()) {
    return false;
  }

  if (!hasHomogeneousDmoMetadata() && !species.isConsistentWith(particle_sidecar)) {
    return false;
  }

  auto& star_rows_by_particle = scratch.star_rows;
  auto& bh_rows_by_particle = scratch.bh_rows;
  auto& tracer_rows_by_particle = scratch.tracer_rows;

  if (patches.size() != 0U) {
    auto& cell_owner = scratch.cell_owner;
    std::fill(cell_owner.begin(), cell_owner.end(),
              std::numeric_limits<std::uint32_t>::max());
    for (std::size_t patch = 0; patch < patches.size(); ++patch) {
      const std::uint64_t begin = patches.first_cell[patch];
      const std::uint64_t count = patches.cell_count[patch];
      if (begin > cells.size() || count > static_cast<std::uint64_t>(cells.size()) - begin) {
        return false;
      }
      const std::uint64_t end = begin + count;
      for (std::uint64_t cell = begin; cell < end; ++cell) {
        const auto cell_index = static_cast<std::size_t>(cell);
        if (cell_owner[cell_index] != std::numeric_limits<std::uint32_t>::max()) {
          return false;
        }
        cell_owner[cell_index] = checkedLocalPatchRow(patch, "SimulationState::validateOwnershipInvariants patch row");
      }
    }
    for (std::size_t cell = 0; cell < cells.size(); ++cell) {
      if (cell_owner[cell] == std::numeric_limits<std::uint32_t>::max() ||
          cells.patch_index[cell] != cell_owner[cell]) {
        return false;
      }
    }
  } else if (cells.size() != 0U) {
    // Non-AMR hydro states are allowed to carry cells without PatchSoa ownership.
    // In that mode patch_index is not authoritative and must remain the neutral zero value.
    for (const std::uint32_t patch_index : cells.patch_index) {
      if (patch_index != 0U) {
        return false;
      }
    }
  }

  for (std::size_t i = 0; i < star_particles.size(); ++i) {
    const auto index = star_particles.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (particleSpeciesTag(index) != static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
      return false;
    }
    if (++star_rows_by_particle[index] != 1) {
      return false;
    }
  }

  for (std::size_t i = 0; i < black_holes.size(); ++i) {
    const auto index = black_holes.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (black_holes.host_cell_index[i] != kInvalidGasCellRow &&
        black_holes.host_cell_index[i] >= cells.size()) {
      return false;
    }
    if (particleSpeciesTag(index) != static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
      return false;
    }
    if (++bh_rows_by_particle[index] != 1) {
      return false;
    }
  }

  for (std::size_t i = 0; i < tracers.size(); ++i) {
    const auto index = tracers.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (tracers.host_cell_index[i] != kInvalidGasCellRow &&
        tracers.host_cell_index[i] >= cells.size()) {
      return false;
    }
    if (tracers.mass_fraction_of_host[i] < 0.0 || tracers.last_host_mass_code[i] < 0.0) {
      return false;
    }
    if (particleSpeciesTag(index) != static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
      return false;
    }
    if (++tracer_rows_by_particle[index] != 1) {
      return false;
    }
  }

  for (std::size_t particle_index = 0; particle_index < particles.size(); ++particle_index) {
    const auto species_tag = particleSpeciesTag(particle_index);
    const bool has_star_row = !star_rows_by_particle.empty() && star_rows_by_particle[particle_index] == 1;
    const bool has_bh_row = !bh_rows_by_particle.empty() && bh_rows_by_particle[particle_index] == 1;
    const bool has_tracer_row = !tracer_rows_by_particle.empty() && tracer_rows_by_particle[particle_index] == 1;

    if (species_tag == static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
      if (!has_star_row || has_bh_row || has_tracer_row) {
        return false;
      }
    } else if (species_tag == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
      if (has_star_row || !has_bh_row || has_tracer_row) {
        return false;
      }
    } else if (species_tag == static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
      if (has_star_row || has_bh_row || !has_tracer_row) {
        return false;
      }
    } else if (has_star_row || has_bh_row || has_tracer_row) {
      return false;
    }
  }

  return true;
}

std::uint64_t SimulationState::particleIndexGeneration() const noexcept {
  return m_particle_index_generation;
}

std::uint64_t SimulationState::particleMassGeneration() const noexcept {
  return m_particle_mass_generation;
}

void SimulationState::notifyParticleMassMutation() noexcept {
  ++m_particle_mass_generation;
}

std::uint64_t SimulationState::cellIndexGeneration() const noexcept {
  return m_cell_index_generation;
}

std::uint64_t SimulationState::gravitySourceGeneration() const noexcept {
  return m_gravity_source_generation;
}

std::uint64_t SimulationState::gasThermodynamicMutationEpoch() const noexcept {
  return m_gas_thermodynamic_mutation_epoch;
}

std::uint64_t SimulationState::gasThermodynamicRefreshEpoch() const noexcept {
  return m_gas_thermodynamic_refresh_epoch;
}

void SimulationState::notifyGasThermodynamicMutation() noexcept {
  if (m_gas_thermodynamic_mutation_epoch != std::numeric_limits<std::uint64_t>::max()) {
    ++m_gas_thermodynamic_mutation_epoch;
  }
}

void SimulationState::markGasThermodynamicsRefreshed() noexcept {
  m_gas_thermodynamic_refresh_epoch = m_gas_thermodynamic_mutation_epoch;
}

bool SimulationState::gasThermodynamicCachesFresh() const noexcept {
  return m_gas_thermodynamic_refresh_epoch == m_gas_thermodynamic_mutation_epoch;
}

void SimulationState::bumpGravitySourceGeneration() noexcept {
  if (m_gravity_source_generation != std::numeric_limits<std::uint64_t>::max()) {
    ++m_gravity_source_generation;
  }
}

void SimulationState::bumpParticleIndexGeneration() noexcept {
  ++m_particle_mass_generation;  // any dense-row change may move masses between rows
  ++m_particle_index_generation;
}

void SimulationState::bumpParticleIdentityGeneration() noexcept {
  if (m_particle_identity_generation != std::numeric_limits<std::uint64_t>::max()) {
    ++m_particle_identity_generation;
  }
  m_particle_identity_validation_certificate = {};
}

void SimulationState::bumpCellIndexGeneration() noexcept {
  ++m_cell_index_generation;
}

}  // namespace cosmosim::core
