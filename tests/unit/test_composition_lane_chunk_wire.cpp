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
//  * patch 52 (P1-51-02): the governed typed exchange returns a lease that
//    covers exactly the received bytes, survives moves, is released with the
//    storage (releaseStorage / destruction) and a refused admission throws
//    without leaving a committed reservation behind;
//  * the C02 dataset energy contract parses ionization potentials and the
//    energy convention, and a dataset without them reports no contract.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
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

void testGovernedLeaseLifetime() {
  namespace core = cosmosim::core;
  const par::MpiContext serial(false, 1, 0);
  struct Rec {
    std::uint64_t id = 0U;
    double value = 0.0;
  };
  const std::vector<Rec> outbound{{1U, 1.0}, {2U, 2.0}, {3U, 3.0}};
  const std::vector<int> destinations{0, 0, 0};
  core::MemoryGovernor governor;
  {
    auto received = par::exchangeGovernedRecords<Rec>(serial, outbound, destinations, &governor, "test.lease");
    CHUI_REQUIRE(received.records.size() == 3U);
    CHUI_REQUIRE(received.lease.committed());
    CHUI_REQUIRE(received.lease.bytes() == 3U * sizeof(Rec));
    CHUI_REQUIRE(governor.snapshot().committed_bytes == 3U * sizeof(Rec));
    par::GovernedRecords<Rec> moved;
    moved = std::move(received);  // ownership transfer keeps one lease
    CHUI_REQUIRE(governor.snapshot().committed_bytes == 3U * sizeof(Rec));
    moved.releaseStorage();
    CHUI_REQUIRE(moved.records.empty() && moved.records.capacity() == 0U);
    CHUI_REQUIRE(governor.snapshot().committed_bytes == 0U);
  }
  CHUI_REQUIRE(governor.snapshot().committed_bytes == 0U);
  core::MemoryGovernorPolicy tight;
  tight.hard_limit_bytes = 8U;  // below one record: admission must refuse
  core::MemoryGovernor refusing(tight);
  CHUI_REQUIRE(throws([&] {
    (void)par::exchangeGovernedRecords<Rec>(serial, outbound, destinations, &refusing, "test.refusal");
  }));
  CHUI_REQUIRE(refusing.snapshot().committed_bytes == 0U && refusing.snapshot().reserved_bytes == 0U);
  // Empty stream on every rank: no admission, no allocation.
  const auto empty = par::exchangeGovernedRecords<Rec>(serial, std::vector<Rec>{}, std::vector<int>{}, &governor,
                                                       "test.empty");
  CHUI_REQUIRE(empty.records.empty() && !empty.lease.valid());
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
  testGovernedLeaseLifetime();
  testIonEnergyContract();
  return 0;
}
