#include "cosmosim/physics/black_hole/bh_engine_state.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <type_traits>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

static_assert(std::is_trivially_copyable_v<BlackHoleEngineRecord>,
              "BlackHoleEngineRecord bytes are the persistent chui_bh_engine_v1 schema");

BlackHoleEngineStore::BlackHoleEngineStore(const core::SimulationState& state) {
  const core::ModuleSidecarBlock* block = state.sidecars.find(k_bh_engine_module);
  if (block == nullptr) return;
  if (!block->particle_indexed || block->schema_version != k_bh_engine_schema ||
      block->row_stride_bytes != sizeof(BlackHoleEngineRecord)) {
    throw std::runtime_error("chui_bh_engine_v1 sidecar has an incompatible schema/stride");
  }
  std::vector<std::pair<std::uint64_t, std::size_t>> order;
  order.reserve(block->rowCount());
  for (std::size_t r = 0; r < block->rowCount(); ++r) order.emplace_back(block->particle_id_by_row[r], r);
  std::sort(order.begin(), order.end());
  m_ids.reserve(order.size());
  m_records.reserve(order.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    if (i > 0U && order[i].first == order[i - 1].first) {
      throw std::runtime_error("chui_bh_engine_v1 sidecar has a duplicated particle id");
    }
    BlackHoleEngineRecord record{};
    const auto bytes = block->rowPayload(order[i].second);
    std::memcpy(&record, bytes.data(), sizeof(record));
    m_ids.push_back(order[i].first);
    m_records.push_back(record);
  }
}

bool BlackHoleEngineStore::has(std::uint64_t particle_id) const {
  return std::binary_search(m_ids.begin(), m_ids.end(), particle_id);
}

BlackHoleEngineRecord BlackHoleEngineStore::get(std::uint64_t particle_id) const {
  const auto it = std::lower_bound(m_ids.begin(), m_ids.end(), particle_id);
  if (it == m_ids.end() || *it != particle_id) {
    throw std::out_of_range("BlackHoleEngineStore: no record for BH particle id");
  }
  return m_records[static_cast<std::size_t>(it - m_ids.begin())];
}

void BlackHoleEngineStore::put(std::uint64_t particle_id, const BlackHoleEngineRecord& record) {
  const auto it = std::lower_bound(m_ids.begin(), m_ids.end(), particle_id);
  const auto index = static_cast<std::size_t>(it - m_ids.begin());
  if (it != m_ids.end() && *it == particle_id) {
    m_records[index] = record;
    return;
  }
  m_ids.insert(it, particle_id);
  m_records.insert(m_records.begin() + static_cast<std::ptrdiff_t>(index), record);
}

void BlackHoleEngineStore::erase(std::uint64_t particle_id) {
  const auto it = std::lower_bound(m_ids.begin(), m_ids.end(), particle_id);
  if (it == m_ids.end() || *it != particle_id) return;
  const auto index = static_cast<std::ptrdiff_t>(it - m_ids.begin());
  m_ids.erase(it);
  m_records.erase(m_records.begin() + index);
}

std::size_t BlackHoleEngineStore::adoptLegacyBlackHoles(const core::SimulationState& state) {
  std::size_t adopted = 0U;
  for (std::size_t b = 0; b < state.black_holes.size(); ++b) {
    const std::uint32_t p = state.black_holes.particle_index[b];
    if (p >= state.particle_sidecar.particle_id.size()) {
      throw std::runtime_error("BlackHoleEngineStore: BH sidecar row references a missing particle");
    }
    const std::uint64_t id = state.particle_sidecar.particle_id[p];
    if (has(id)) continue;
    BlackHoleEngineRecord record{};
    // Adoption policy: the legacy subgrid mass is the horizon mass; any
    // difference between particle and subgrid mass (legacy runs keep them
    // equal) is NOT reinterpreted as a disk; it is refused below.
    record.horizon_mass_code = state.black_holes.subgrid_mass_code[b];
    const double particle_mass = state.particles.mass_code[p];
    const double tol = 1.0e-12 * std::max(std::abs(particle_mass), 1.0e-300);
    if (std::abs(particle_mass - record.horizon_mass_code) > tol) {
      throw std::runtime_error(
          "BlackHoleEngineStore: cannot adopt a legacy BH whose particle mass differs from its subgrid mass "
          "(no rule exists to reinterpret the difference as disk or ejecta mass)");
    }
    record.flags = kBhFlagInitialized;
    put(id, record);
    ++adopted;
  }
  return adopted;
}

void BlackHoleEngineStore::store(core::SimulationState& state) const {
  core::ModuleSidecarBlock block;
  block.module_name = k_bh_engine_module;
  block.schema_version = k_bh_engine_schema;
  block.particle_indexed = true;
  block.row_stride_bytes = static_cast<std::uint32_t>(sizeof(BlackHoleEngineRecord));
  block.requirement.kind = core::ModuleSidecarRequirementKind::kSpeciesMask;
  block.requirement.species_mask = 1U << static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);
  std::vector<std::uint64_t> live;
  live.reserve(state.black_holes.size());
  for (std::size_t b = 0; b < state.black_holes.size(); ++b) {
    const std::uint32_t p = state.black_holes.particle_index[b];
    if (p < state.particle_sidecar.particle_id.size()) live.push_back(state.particle_sidecar.particle_id[p]);
  }
  std::sort(live.begin(), live.end());
  for (const std::uint64_t id : live) {
    if (!has(id)) {
      throw std::runtime_error("BlackHoleEngineStore: a live BH has no chui_bh_engine_v1 record");
    }
  }
  block.particle_id_by_row.reserve(live.size());
  block.payload.reserve(live.size() * sizeof(BlackHoleEngineRecord));
  for (std::size_t i = 0; i < m_ids.size(); ++i) {
    if (!std::binary_search(live.begin(), live.end(), m_ids[i])) continue;
    block.particle_id_by_row.push_back(m_ids[i]);
    const auto* bytes = reinterpret_cast<const std::byte*>(&m_records[i]);
    block.payload.insert(block.payload.end(), bytes, bytes + sizeof(BlackHoleEngineRecord));
  }
  state.sidecars.upsert(std::move(block));
}

}  // namespace cosmosim::physics
