#include "cosmosim/physics/baryon_composition.hpp"

#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <cmath>
#include <stdexcept>
#include <string>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/retained_capacity_transaction.hpp"
#include "cosmosim/physics/metal_ion_network.hpp"

namespace cosmosim::physics {
namespace {

constexpr std::array<ChemicalElementInfo, k_chemical_element_count> k_elements{{
    {"H", 1.00794, 12.00},
    {"He", 4.002602, 10.93},
    {"C", 12.0107, 8.43},
    {"N", 14.0067, 7.83},
    {"O", 15.9994, 8.69},
    {"Ne", 20.1797, 7.93},
    {"Mg", 24.3050, 7.60},
    {"Si", 28.0855, 7.51},
    {"S", 32.065, 7.12},
    {"Ca", 40.078, 6.34},
    {"Fe", 55.845, 7.50},
}};

constexpr double A_H = 1.00794;
constexpr double A_He = 4.002602;
constexpr double A_C = 12.0107;
constexpr double A_O = 15.9994;

// Species mass is the neutral-atomic-weight sum of the nuclear content so the
// species lanes of one element sum EXACTLY to that element lane. The electron
// mass (5.5e-4 amu per free electron) is not tracked as a separate reservoir;
// this is a documented O(1e-4) mass-accounting convention, not a leak:
// element nuclei are conserved exactly.
constexpr std::array<ChemicalSpeciesInfo, k_chemical_species_count> k_species{{
    {"HI", 1, 0, 0, 0, 0, A_H},
    {"HII", 1, 0, 0, 0, +1, A_H},
    {"H-", 1, 0, 0, 0, -1, A_H},
    {"H2", 2, 0, 0, 0, 0, 2.0 * A_H},
    {"H2+", 2, 0, 0, 0, +1, 2.0 * A_H},
    {"HeI", 0, 1, 0, 0, 0, A_He},
    {"HeII", 0, 1, 0, 0, +1, A_He},
    {"HeIII", 0, 1, 0, 0, +2, A_He},
    {"CII", 0, 0, 1, 0, +1, A_C},
    {"CO", 0, 0, 1, 1, 0, A_C + A_O},
}};

[[nodiscard]] std::array<double, k_chemical_element_count> computeSolarMassFractions() {
  std::array<double, k_chemical_element_count> fractions{};
  double total = 0.0;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    const double number_per_h = std::pow(10.0, k_elements[e].solar_log_epsilon - 12.0);
    fractions[e] = number_per_h * k_elements[e].atomic_mass_amu;
    total += fractions[e];
  }
  for (double& f : fractions) f /= total;
  return fractions;
}

[[nodiscard]] std::array<double, k_chemical_element_count> stoichiometricFractions(
    std::initializer_list<std::pair<ChemicalElement, double>> atoms) {
  std::array<double, k_chemical_element_count> fractions{};
  double total = 0.0;
  for (const auto& [element, count] : atoms) {
    const double mass = count * k_elements[static_cast<std::size_t>(element)].atomic_mass_amu;
    fractions[static_cast<std::size_t>(element)] += mass;
    total += mass;
  }
  for (double& f : fractions) f /= total;
  return fractions;
}

}  // namespace

const std::array<ChemicalElementInfo, k_chemical_element_count>& chemicalElementTable() {
  return k_elements;
}

const ChemicalElementInfo& chemicalElementInfo(ChemicalElement element) {
  return k_elements.at(static_cast<std::size_t>(element));
}

const std::array<double, k_chemical_element_count>& solarElementMassFractions() {
  static const std::array<double, k_chemical_element_count> fractions = computeSolarMassFractions();
  return fractions;
}

double solarMetalMassFraction() {
  const auto& f = solarElementMassFractions();
  double z = 0.0;
  for (std::size_t e = 2; e < k_chemical_element_count; ++e) z += f[e];
  return z;
}

const ChemicalSpeciesInfo& chemicalSpeciesInfo(ChemicalSpecies species) {
  return k_species.at(static_cast<std::size_t>(species));
}

CompositionLayout CompositionLayout::make(const CompositionLayoutRequest& request) {
  if (request.adaptive_sparse_detailed &&
      request.network != CompositionNetworkKind::kDetailedCarbonOxygen) {
    throw std::invalid_argument(
        "adaptive sparse chemistry storage requires the detailed C/O network");
  }
  CompositionLayout layout;
  layout.m_request = request;
  for (const auto& element : k_elements) {
    layout.m_lanes.push_back({"element_" + std::string(element.symbol), "mass_code"});
  }
  std::size_t species_count = 0U;
  std::string name = "chui_elements11";
  switch (request.network) {
    case CompositionNetworkKind::kElementsOnly: species_count = 0U; break;
    case CompositionNetworkKind::kHybridLight:
      species_count = k_hybrid_species_count;
      name += "_hybrid_light_v1";
      break;
    case CompositionNetworkKind::kDetailedCarbonOxygen:
      species_count = k_detailed_species_count;
      name += "_detailed_co_v1";
      break;
  }
  const auto species_begin = static_cast<std::uint32_t>(layout.m_lanes.size());
  for (std::size_t s = 0; s < species_count; ++s) {
    layout.m_lanes.push_back({"species_" + std::string(k_species[s].name), "mass_code"});
  }
  const auto ion_begin = static_cast<std::uint32_t>(layout.m_lanes.size());
  const MetalIonElementSet ion_set = MetalIonElementSet::fromMask(request.metal_ion_element_mask);
  if (!ion_set.empty()) {
    if (request.network != CompositionNetworkKind::kDetailedCarbonOxygen) {
      throw std::invalid_argument(
          "metal-ion ladders require the detailed C/O composition network (fixed_detailed)");
    }
    for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
      const auto element = static_cast<MetalIonElement>(e);
      if (!ion_set.contains(element)) continue;
      for (std::size_t k = 0; k < metalIonStageCount(element); ++k) {
        layout.m_lanes.push_back({"ion_" + metalIonStageName(element, k), "mass_code"});
      }
    }
    name += "_ions_" + ion_set.label() + "_v1";
  }
  const auto ion_count = static_cast<std::uint32_t>(ion_set.stageCount());
  const auto dust_begin = static_cast<std::uint32_t>(layout.m_lanes.size());
  std::size_t dust_count = 0U;
  if (request.dust == CompositionDustKind::kTwoReservoirV1) {
    layout.m_lanes.push_back({"dust_silicate_MgFeSiO4", "mass_code"});
    layout.m_lanes.push_back({"dust_carbonaceous_C", "mass_code"});
    dust_count = 2U;
    name += "_dust2_v1";
  }
  std::uint32_t sgs_lane = core::GasCompositionLayoutDescriptor::k_no_lane;
  if (request.sgs_energy) {
    sgs_lane = static_cast<std::uint32_t>(layout.m_lanes.size());
    layout.m_lanes.push_back({"sgs_turbulent_energy", "energy_code"});
    name += "_sgs_v1";
  }
  static_assert(k_metal_ion_stage_capacity == 7U + 9U + 15U + 27U);
  if (layout.m_lanes.size() > k_max_constructible_composition_lane_count) {
    throw std::logic_error("CompositionLayout::make exceeded the constructible composition lane bound");
  }
  std::string signature_text = name;
  for (const auto& lane : layout.m_lanes) {
    signature_text += "|" + lane.name + ":" + lane.unit;
  }
  auto& d = layout.m_descriptor;
  d.layout_version = k_layout_version;
  d.layout_signature = core::fnv1a64(signature_text);
  d.layout_name = name;
  d.lane_count = static_cast<std::uint32_t>(layout.m_lanes.size());
  d.element_lane_begin = 0U;
  d.element_lane_count = static_cast<std::uint32_t>(k_chemical_element_count);
  d.species_lane_begin = species_begin;
  d.species_lane_count = static_cast<std::uint32_t>(species_count);
  d.dust_lane_begin = dust_begin;
  d.dust_lane_count = static_cast<std::uint32_t>(dust_count);
  d.ion_lane_begin = ion_count == 0U ? 0U : ion_begin;
  d.ion_lane_count = ion_count;
  d.sgs_energy_lane = sgs_lane;
  if (request.adaptive_sparse_detailed) {
    layout.m_sparse_default_source_lanes = {
        d.element_lane_begin + static_cast<std::uint32_t>(ChemicalElement::kC),
        core::GasCompositionLayoutDescriptor::k_no_lane};
    // Ion lanes follow C+/CO contiguously; implicit value zero.
    layout.m_sparse_default_source_lanes.insert(
        layout.m_sparse_default_source_lanes.end(), ion_count,
        core::GasCompositionLayoutDescriptor::k_no_lane);
  }
  if (!d.isConsistent()) {
    throw std::logic_error("CompositionLayout::make produced an inconsistent descriptor");
  }
  return layout;
}

CompositionLayout CompositionLayout::fromDescriptor(const core::GasCompositionLayoutDescriptor& descriptor) {
  CompositionLayoutRequest request;
  request.network = descriptor.species_lane_count == 0U ? CompositionNetworkKind::kElementsOnly
      : descriptor.species_lane_count == k_hybrid_species_count ? CompositionNetworkKind::kHybridLight
                                                                : CompositionNetworkKind::kDetailedCarbonOxygen;
  request.dust = descriptor.dust_lane_count == 0U ? CompositionDustKind::kNone : CompositionDustKind::kTwoReservoirV1;
  request.sgs_energy = descriptor.sgs_energy_lane != core::GasCompositionLayoutDescriptor::k_no_lane;
  request.metal_ion_element_mask =
      MetalIonElementSet::fromStageCount(descriptor.ion_lane_count).mask();
  CompositionLayout layout = make(request);
  if (!(layout.descriptor() == descriptor)) {
    throw std::runtime_error("composition descriptor '" + descriptor.layout_name +
                             "' is not a known CHUI layout (signature mismatch)");
  }
  return layout;
}

std::uint32_t CompositionLayout::speciesLane(std::size_t species_slot) const {
  if (species_slot >= m_descriptor.species_lane_count) {
    throw std::out_of_range("CompositionLayout: species slot is not part of the active network");
  }
  return m_descriptor.species_lane_begin + static_cast<std::uint32_t>(species_slot);
}

std::uint32_t CompositionLayout::ionStageLane(std::size_t packed_stage) const {
  if (packed_stage >= m_descriptor.ion_lane_count) {
    throw std::out_of_range("CompositionLayout: metal-ion stage is not part of the active layout");
  }
  return m_descriptor.ion_lane_begin + static_cast<std::uint32_t>(packed_stage);
}

std::uint32_t CompositionLayout::dustLane(std::size_t dust_slot) const {
  if (dust_slot >= m_descriptor.dust_lane_count) {
    throw std::out_of_range("CompositionLayout: dust slot is not part of the active layout");
  }
  return m_descriptor.dust_lane_begin + static_cast<std::uint32_t>(dust_slot);
}

const std::array<double, k_chemical_element_count>& dustElementMassFractions(DustReservoir reservoir) {
  static const auto silicate = stoichiometricFractions({{ChemicalElement::kMg, 1.0},
                                                        {ChemicalElement::kFe, 1.0},
                                                        {ChemicalElement::kSi, 1.0},
                                                        {ChemicalElement::kO, 4.0}});
  static const auto carbon = stoichiometricFractions({{ChemicalElement::kC, 1.0}});
  return reservoir == DustReservoir::kSilicate ? silicate : carbon;
}

std::array<double, k_chemical_element_count> initialElementMasses(
    double gas_mass_code, double metal_mass_code, const PrimordialCompositionPolicy& policy) {
  if (!(gas_mass_code >= 0.0) || !(metal_mass_code >= 0.0) || metal_mass_code > gas_mass_code) {
    throw std::invalid_argument("initialElementMasses requires 0 <= metal mass <= gas mass");
  }
  if (!(policy.helium_mass_fraction > 0.0) || policy.helium_mass_fraction >= 1.0) {
    throw std::invalid_argument("initialElementMasses requires 0 < Y_p < 1");
  }
  std::array<double, k_chemical_element_count> masses{};
  if (gas_mass_code == 0.0) return masses;
  const auto& solar = solarElementMassFractions();
  const double z_sun = solarMetalMassFraction();
  const double z = metal_mass_code / gas_mass_code;
  for (std::size_t e = 2; e < k_chemical_element_count; ++e) {
    masses[e] = metal_mass_code * solar[e] / z_sun;
  }
  double y = policy.helium_mass_fraction +
      (solar[static_cast<std::size_t>(ChemicalElement::kHe)] - policy.helium_mass_fraction) * z / z_sun;
  y = std::clamp(y, 0.0, std::max(0.0, 1.0 - z));
  masses[static_cast<std::size_t>(ChemicalElement::kHe)] = y * gas_mass_code;
  masses[static_cast<std::size_t>(ChemicalElement::kH)] =
      std::max(0.0, gas_mass_code - metal_mass_code - y * gas_mass_code);
  return masses;
}

void prolongComposition(std::span<const double> parent_lanes,
                        std::span<const double> child_mass_share,
                        std::span<double> child_lanes_row_major) {
  const std::size_t lanes = parent_lanes.size();
  if (child_lanes_row_major.size() != lanes * child_mass_share.size()) {
    throw std::invalid_argument("prolongComposition: child buffer size mismatch");
  }
  double share_sum = 0.0;
  for (const double share : child_mass_share) {
    if (!(share >= 0.0)) throw std::invalid_argument("prolongComposition: negative mass share");
    share_sum += share;
  }
  if (std::abs(share_sum - 1.0) > 1.0e-12) {
    throw std::invalid_argument("prolongComposition: child mass shares must sum to one");
  }
  for (std::size_t c = 0; c < child_mass_share.size(); ++c) {
    for (std::size_t l = 0; l < lanes; ++l) {
      child_lanes_row_major[c * lanes + l] = parent_lanes[l] * child_mass_share[c];
    }
  }
}

void restrictComposition(std::span<const double> child_lanes_row_major,
                         std::size_t lane_count,
                         std::span<double> parent_lanes) {
  if (lane_count == 0U || parent_lanes.size() != lane_count ||
      child_lanes_row_major.size() % lane_count != 0U) {
    throw std::invalid_argument("restrictComposition: buffer size mismatch");
  }
  std::fill(parent_lanes.begin(), parent_lanes.end(), 0.0);
  const std::size_t children = child_lanes_row_major.size() / lane_count;
  for (std::size_t c = 0; c < children; ++c) {
    for (std::size_t l = 0; l < lane_count; ++l) {
      parent_lanes[l] += child_lanes_row_major[c * lane_count + l];
    }
  }
}

double compositionGasMass(const core::GasCompositionSidecar& sidecar, std::size_t row) {
  const auto& d = sidecar.descriptor();
  double mass = sidecar.elementMassSum(row);
  for (std::uint32_t l = 0; l < d.dust_lane_count; ++l) {
    mass += sidecar.value(d.dust_lane_begin + l, row);
  }
  return mass;
}

double synchronizeMetalMirrorFromComposition(core::SimulationState& state,
                                             const CompositionLayout& layout,
                                             std::span<const std::uint32_t> rows) {
  if (state.gas_composition.enabled() && !(state.gas_composition.descriptor() == layout.descriptor())) {
    throw std::runtime_error("synchronizeMetalMirrorFromComposition: layout mismatch");
  }
  return synchronizeMetalMirrorFromComposition(state, rows);
}

double synchronizeMetalMirrorFromComposition(core::SimulationState& state,
                                             std::span<const std::uint32_t> rows) {
  auto& sidecar = state.gas_composition;
  if (!sidecar.enabled()) return 0.0;
  const auto& descriptor = sidecar.descriptor();
  if (descriptor.element_lane_count != k_chemical_element_count ||
      descriptor.dust_lane_count > k_dust_reservoir_count) {
    throw std::runtime_error("synchronizeMetalMirrorFromComposition: unsupported layout");
  }
  double max_drift = 0.0;
  const auto sync_row = [&](std::size_t row) {
    double metals = 0.0;
    for (std::size_t e = 2; e < k_chemical_element_count; ++e) {
      metals += sidecar.value(descriptor.element_lane_begin + e, row);
    }
    for (std::size_t d = 0; d < descriptor.dust_lane_count; ++d) {
      const auto& frac = dustElementMassFractions(static_cast<DustReservoir>(d));
      double metal_fraction = 0.0;
      for (std::size_t e = 2; e < k_chemical_element_count; ++e) metal_fraction += frac[e];
      metals += metal_fraction * sidecar.value(descriptor.dust_lane_begin + d, row);
    }
    const double old_value = state.gas_cells.metal_mass_code[row];
    const double scale = std::max({std::abs(old_value), std::abs(metals), 1.0e-300});
    max_drift = std::max(max_drift, std::abs(old_value - metals) / scale);
    state.gas_cells.metal_mass_code[row] = metals;
  };
  if (rows.empty()) {
    for (std::size_t row = 0; row < sidecar.rowCount(); ++row) sync_row(row);
  } else {
    for (const std::uint32_t row : rows) sync_row(row);
  }
  return max_drift;
}

void enableGasComposition(core::SimulationState& state, const CompositionLayout& layout,
                          const PrimordialCompositionPolicy& policy) {
  auto& sidecar = state.gas_composition;
  if (sidecar.enabled()) {
    if (!(sidecar.descriptor() == layout.descriptor()) ||
        sidecar.sparseStorageEnabled() != layout.sparseDetailedStorage()) {
      throw std::runtime_error(
          "enableGasComposition: state already carries composition layout '" +
          sidecar.descriptor().layout_name + "' but '" + layout.descriptor().layout_name +
          "' with a different adaptive storage policy was requested; an explicit, "
          "versioned conversion is required");
    }
    return;
  }
  const std::size_t rows = state.cells.size();
  sidecar.configure(layout.descriptor(), rows, layout.sparseLaneBegin(),
                    layout.sparseDefaultSourceLanes());
  for (std::size_t row = 0; row < rows; ++row) {
    sidecar.gas_cell_id[row] = state.gas_cells.gas_cell_id[row];
    if (sidecar.sparseStorageEnabled()) {
      sidecar.chemistry_model_tag[row] = k_chemistry_tag_hybrid;
    }
    const double mass = std::max(state.cells.mass_code[row], 0.0);
    const double metals = std::clamp(state.gas_cells.metal_mass_code[row], 0.0, mass);
    const auto elements = initialElementMasses(mass, metals, policy);
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
      sidecar.setValue(layout.elementLane(static_cast<ChemicalElement>(e)), row, elements[e]);
    }
    // Species start in the neutral atomic state; the first thermochemistry
    // solve relaxes them (an explicit, documented initialization; equilibrium
    // initialization is available through the equilibrium backend).
    const std::size_t species = sidecar.descriptor().species_lane_count;
    if (species >= k_hybrid_species_count) {
      sidecar.setValue(layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kHI)), row,
                       elements[static_cast<std::size_t>(ChemicalElement::kH)]);
      sidecar.setValue(layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kHeI)), row,
                       elements[static_cast<std::size_t>(ChemicalElement::kHe)]);
    }
    if (species >= k_detailed_species_count && !sidecar.sparseStorageEnabled()) {
      sidecar.setValue(layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kCII)), row,
                       elements[static_cast<std::size_t>(ChemicalElement::kC)]);
    }
  }
}

std::array<double, k_chemical_element_count> removeGasFractionFromComposition(
    core::GasCompositionSidecar& sidecar, const CompositionLayout& layout, std::size_t row,
    double removed_fraction) {
  if (!(removed_fraction >= 0.0) || removed_fraction > 1.0) {
    throw std::invalid_argument("removeGasFractionFromComposition: fraction must be in [0, 1]");
  }
  std::array<double, k_chemical_element_count> removed{};
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    removed[e] = removed_fraction *
        sidecar.value(layout.elementLane(static_cast<ChemicalElement>(e)), row);
  }
  for (std::size_t d = 0; d < sidecar.descriptor().dust_lane_count; ++d) {
    const double dust = removed_fraction * sidecar.value(layout.dustLane(d), row);
    const auto& frac = dustElementMassFractions(static_cast<DustReservoir>(d));
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) removed[e] += dust * frac[e];
  }
  sidecar.scaleRow(row, 1.0 - removed_fraction);
  return removed;
}

namespace {
template <class T>
void appendPod(std::vector<std::byte>& out, const T& value) {
  const auto* bytes = reinterpret_cast<const std::byte*>(&value);
  out.insert(out.end(), bytes, bytes + sizeof(T));
}
template <class T>
T readPod(const std::vector<std::byte>& in, std::size_t& offset) {
  if (offset + sizeof(T) > in.size()) throw std::runtime_error("chui_gas_composition_v1 block is truncated");
  T value;
  std::memcpy(&value, in.data() + offset, sizeof(T));
  offset += sizeof(T);
  return value;
}
}  // namespace

void packGasCompositionIntoModuleSidecar(core::SimulationState& state) {
  const auto& sidecar = state.gas_composition;
  if (!sidecar.enabled()) return;
  sidecar.requireAlignedWith(state.gas_cells.gas_cell_id, "packGasCompositionIntoModuleSidecar");
  core::ModuleSidecarBlock block;
  block.module_name = k_gas_composition_module;
  block.schema_version = k_gas_composition_schema;
  const auto& d = sidecar.descriptor();
  appendPod(block.payload, d.layout_version);
  appendPod(block.payload, d.layout_signature);
  appendPod(block.payload, d.lane_count);
  appendPod(block.payload, static_cast<std::uint64_t>(sidecar.rowCount()));
  for (std::size_t row = 0; row < sidecar.rowCount(); ++row) {
    appendPod(block.payload, sidecar.gas_cell_id[row]);
    for (const double v : sidecar.packRow(row)) appendPod(block.payload, v);
  }
  state.sidecars.upsert(std::move(block));
}

bool restoreGasCompositionFromModuleSidecar(core::SimulationState& state, const CompositionLayout& layout) {
  const core::ModuleSidecarBlock* block = state.sidecars.find(k_gas_composition_module);
  if (block == nullptr) return false;
  if (block->schema_version != 1U && block->schema_version != 2U &&
      block->schema_version != k_gas_composition_schema) {
    throw std::runtime_error("unsupported chui_gas_composition schema version");
  }
  std::size_t offset = 0;
  const auto version = readPod<std::uint32_t>(block->payload, offset);
  const auto signature = readPod<std::uint64_t>(block->payload, offset);
  const auto lanes = readPod<std::uint32_t>(block->payload, offset);
  const auto rows = readPod<std::uint64_t>(block->payload, offset);
  const auto& d = layout.descriptor();
  if (version != d.layout_version || signature != d.layout_signature || lanes != d.lane_count) {
    throw std::runtime_error(
        "checkpointed gas composition layout differs from the configured layout '" + d.layout_name +
        "'; changing the composition model on resume requires an explicit, versioned conversion");
  }
  if (rows != state.cells.size()) {
    throw std::runtime_error("checkpointed gas composition row count differs from the restored gas cells");
  }
  core::GasCompositionSidecar restored;
  restored.configure(d, state.cells.size(), layout.sparseLaneBegin(),
                     layout.sparseDefaultSourceLanes());
  std::unordered_map<std::uint64_t, std::uint32_t> row_by_id;
  row_by_id.reserve(state.cells.size());
  for (std::uint32_t r = 0; r < state.cells.size(); ++r) row_by_id.emplace(state.gas_cells.gas_cell_id[r], r);
  const std::size_t stored_width = block->schema_version == 1U
      ? restored.laneCount() + 4U
      : block->schema_version == 2U
          ? restored.laneCount() + 11U
          : restored.packedRowWidth();
  if (restored.sparseStorageEnabled()) {
    std::size_t scan = offset;
    std::size_t active_rows = 0U;
    for (std::uint64_t i = 0; i < rows; ++i) {
      (void)readPod<std::uint64_t>(block->payload, scan);
      double carbon = 0.0;
      double cii = 0.0;
      double co = 0.0;
      double chemistry_tag = static_cast<double>(k_chemistry_tag_hybrid);
      double explicit_marker = 0.0;
      for (std::size_t value_index = 0; value_index < stored_width; ++value_index) {
        const double value = readPod<double>(block->payload, scan);
        if (value_index == layout.elementLane(ChemicalElement::kC)) carbon = value;
        if (value_index == restored.sparseLaneBegin()) cii = value;
        if (value_index == restored.sparseLaneBegin() + 1U) co = value;
        if (value_index == restored.laneCount()) chemistry_tag = value;
        if (block->schema_version == 3U && value_index + 1U == stored_width) {
          explicit_marker = value;
        }
      }
      const bool active = explicit_marker != 0.0 || cii != carbon || co != 0.0 ||
          (block->schema_version < 3U &&
           chemistry_tag == static_cast<double>(k_chemistry_tag_detailed));
      if (active) ++active_rows;
    }
    core::RetainedCapacityTransaction sparse_capacity(
        [&restored] { return restored.ownedCapacityBytes(); });
    restored.planSparseRowCapacity(sparse_capacity, active_rows);
    sparse_capacity.execute(nullptr, core::MemoryClass::kCanonicalPersistent,
                            "restart.adaptive_chemistry_rows");
  }
  std::vector<double> stored_payload(stored_width);
  std::vector<double> payload(restored.packedRowWidth(), 0.0);
  std::vector<std::uint8_t> filled(state.cells.size(), 0U);
  for (std::uint64_t i = 0; i < rows; ++i) {
    const auto id = readPod<std::uint64_t>(block->payload, offset);
    for (double& v : stored_payload) v = readPod<double>(block->payload, offset);
    std::fill(payload.begin(), payload.end(), 0.0);
    std::copy(stored_payload.begin(), stored_payload.end(), payload.begin());
    const auto it = row_by_id.find(id);
    if (it == row_by_id.end() || filled[it->second] != 0U) {
      throw std::runtime_error("checkpointed gas composition row " + std::to_string(id) +
                               " has no unique restored gas cell");
    }
    // Schema 1/2 stored every lane densely. Preserve a detailed row whenever
    // its values differ from the sparse projection; schema 3 carries the
    // explicit-history marker even when the current values equal defaults.
    if (block->schema_version < 3U && restored.sparseStorageEnabled()) {
      const std::size_t cii = restored.sparseLaneBegin();
      const double implicit_cii = payload[layout.elementLane(ChemicalElement::kC)];
      if (payload[cii] != implicit_cii || payload[cii + 1U] != 0.0 ||
          payload[restored.laneCount()] == static_cast<double>(k_chemistry_tag_detailed)) {
        payload.back() = 1.0;
      }
    }
    restored.unpackRow(it->second, payload, id);
    filled[it->second] = 1U;
  }
  if (offset != block->payload.size()) {
    throw std::runtime_error(
        "checkpointed gas composition block has trailing bytes for its schema");
  }
  state.gas_composition = std::move(restored);
  return true;
}

}  // namespace cosmosim::physics
