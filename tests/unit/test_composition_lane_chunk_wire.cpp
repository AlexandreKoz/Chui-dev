// CHUI patch 50 — P0-49-01 chunked composition wire (NOT RUN — USER-RESERVED).
//
// Invariants:
//  * every constructible composition layout (up to the C/O/Si/Fe + dust + SGS
//    82-lane layout) is representable: split -> exchange (single rank) ->
//    reassemble returns the exact lane vector and per-cell scalar;
//  * reassembly is order independent and refuses a missing chunk, a duplicate
//    chunk, a foreign layout signature, a stale sequence and a misaddressed
//    record (exactly-once application per (cell, source));
//  * the static wire bound covers the largest constructible layout;
//  * the C02 dataset energy contract parses ionization potentials and the
//    energy convention, and a dataset without them reports no contract.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <vector>

#include "cosmosim/parallel/composition_exchange.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/metal_ion_network.hpp"
#include "../support/test_require.hpp"

namespace {

namespace par = cosmosim::parallel;
namespace ph = cosmosim::physics;

template <class F>
bool throws(F&& f) {
  try {
    f();
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

void testStaticBound() {
  CHUI_REQUIRE(par::k_max_composition_layout_lanes >= ph::k_max_constructible_composition_lane_count);
  ph::CompositionLayoutRequest request;
  request.network = ph::CompositionNetworkKind::kDetailedCarbonOxygen;
  request.metal_ion_element_mask = 0x0FU;  // C, O, Si, Fe
  request.dust = ph::CompositionDustKind::kTwoReservoirV1;
  request.sgs_energy = true;
  const ph::CompositionLayout layout = ph::CompositionLayout::make(request);
  CHUI_REQUIRE(layout.descriptor().lane_count == ph::k_max_constructible_composition_lane_count);
  CHUI_REQUIRE(par::compositionLaneChunkCount(layout.descriptor().lane_count) == 3U);
}

void testRoundTripAndRefusals() {
  const par::MpiContext serial(false, 1, 0);
  for (const std::uint32_t lanes : {1U, 31U, 32U, 33U, 37U, 64U, 79U, 82U}) {
    const par::CompositionLaneStreamIdentity identity{
        .stream = par::CompositionLaneStream::kCorrection,
        .layout_signature = 0x1234ULL,
        .lane_count = lanes,
        .exchange_sequence = 7U};
    std::vector<par::CompositionLaneChunkRecord> records;
    std::vector<int> destinations;
    for (std::uint64_t id : {5ULL, 3ULL}) {
      std::vector<double> values(lanes);
      for (std::uint32_t l = 0; l < lanes; ++l) values[l] = static_cast<double>(id * 1000U + l);
      par::appendCompositionLaneChunks(identity, id, 0, 0, static_cast<double>(id), values, records, destinations);
    }
    std::reverse(records.begin(), records.end());  // arrival order must not matter
    std::reverse(destinations.begin(), destinations.end());
    const auto received = par::exchangeCompositionLaneChunks(serial, records, destinations, 7U).records;
    const auto assembled = par::assembleCompositionLaneChunks(received, identity, 0);
    CHUI_REQUIRE(assembled.size() == 2U);
    CHUI_REQUIRE(assembled.gas_cell_id[0] == 3U && assembled.gas_cell_id[1] == 5U);
    for (std::size_t e = 0; e < 2U; ++e) {
      CHUI_REQUIRE(assembled.scalar[e] == static_cast<double>(assembled.gas_cell_id[e]));
      const auto lv = assembled.lanesOf(e);
      for (std::uint32_t l = 0; l < lanes; ++l) {
        CHUI_REQUIRE(lv[l] == static_cast<double>(assembled.gas_cell_id[e] * 1000U + l));
      }
    }
    if (records.size() > 2U) {
      auto missing = received;
      missing.pop_back();
      CHUI_REQUIRE(throws([&] { (void)par::assembleCompositionLaneChunks(missing, identity, 0); }));
      auto duplicate = received;
      duplicate.push_back(received.front());
      duplicate.push_back(received.front());
      CHUI_REQUIRE(throws([&] { (void)par::assembleCompositionLaneChunks(duplicate, identity, 0); }));
    }
    auto foreign = identity;
    foreign.layout_signature = 0x9999ULL;
    CHUI_REQUIRE(throws([&] { (void)par::assembleCompositionLaneChunks(received, foreign, 0); }));
    auto stale = identity;
    stale.exchange_sequence = 8U;
    CHUI_REQUIRE(throws([&] { (void)par::assembleCompositionLaneChunks(received, stale, 0); }));
    CHUI_REQUIRE(throws([&] { (void)par::assembleCompositionLaneChunks(received, identity, 1); }));
  }
  CHUI_REQUIRE(throws([] { (void)par::compositionLaneChunkCount(par::k_max_composition_layout_lanes + 1U); }));
}

void testIonEnergyContract() {
  const char* base =
      "schema chui_metal_ion_rates_v1\n"
      "provenance synthetic unit-test fixture (not physical data)\n"
      "license test-only\n"
      "verification unverified\n"
      "temperature_grid_log10_k 2 4 5\n";
  std::string text = base;
  for (int k = 0; k < 7; ++k) {
    const std::string stage = std::to_string(k);
    text += "rate C " + stage + " collisional_ionization " + (k == 6 ? "-inf -inf" : "-10 -9") + "\n";
    text += "rate C " + stage + " recombination " + (k == 0 ? "-inf -inf" : "-12 -12") + "\n";
    text += "rate C " + stage + " cooling -23 -22\n";
  }
  const auto no_contract = ph::MetalIonRateTable::parse(text, "fixture", ph::MetalIonElementSet::parse("C"));
  CHUI_REQUIRE(!no_contract.hasIonizationEnergyContract(ph::MetalIonElement::kC));
  std::string with = text + "energy_convention collisional_ionization_work_in_cooling\n";
  // Fixture values (arbitrary positive test numbers, not atomic data).
  for (int k = 0; k < 6; ++k) with += "ionization_potential C " + std::to_string(k) + " " + std::to_string(k + 1) + "\n";
  const auto contract = ph::MetalIonRateTable::parse(with, "fixture", ph::MetalIonElementSet::parse("C"));
  CHUI_REQUIRE(contract.hasIonizationEnergyContract(ph::MetalIonElement::kC));
  std::vector<double> all_bare(7, 0.0);
  all_bare[6] = 1.0;
  const double stored = contract.storedIonizationEnergyPerNucleusErg(ph::MetalIonElement::kC, all_bare);
  CHUI_REQUIRE(std::abs(stored - 21.0 * 1.602176634e-12) <= 1.0e-24);
  CHUI_REQUIRE(throws([&] {
    (void)ph::MetalIonRateTable::parse(text + "ionization_potential C 6 1.0\n", "fixture",
                                       ph::MetalIonElementSet::parse("C"));
  }));
}

}  // namespace

int main() {
  testStaticBound();
  testRoundTripAndRefusals();
  testIonEnergyContract();
  return 0;
}
