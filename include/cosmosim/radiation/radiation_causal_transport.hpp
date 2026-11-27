#pragma once

// CHUI-EXT-PHYS-REF-v1 R04/R05: causal ray transport and anisotropic-scattering
// Monte Carlo on the shared read-only leaf geometry (lightcone::LeafGeometry).
// Reference backends, not dispatched by RadiationRuntime in pass 01.
//
// R04 causal rays (time dependent, distinct from the steady-state
// infinite-speed traceDirectRays): a packet emitted at t_e from x_e along n
// is at x_e + c (t - t_e) n at time t. advance(t_end) moves every packet the
// remaining c (t_end - t_now) through the leaf cells, absorbing
// E (1 - exp(-chi l)) per segment into the cell ledger (exact exponential),
// and ledgers escape at open faces. A packet whose energy falls below
// min_energy_fraction of its emitted energy is terminated into the explicit
// `unresolved` ledger (reported, never dropped). Exactly-once closure:
//   emitted = absorbed + escaped + unresolved + in_flight   (to round-off).
// The packet queue is bounded (max_packets): an emission that would exceed it
// is REFUSED before any state changes. Packet state (position history,
// energy, emission time, ids) serializes for restart.
//
// R05 Monte Carlo: analog packets with absorption and Henyey-Greenstein
// anisotropic elastic scattering (asymmetry g in (-1, 1); g = 0 isotropic),
//   mu = (1 + g^2 - ((1 - g^2)/(1 - g + 2 g xi))^2) / (2 g),
// sampled optical depths tau = -ln xi, counter-based reproducible RNG keyed
// by (seed, packet id, event), deposition of absorbed packets in their cell,
// escape at open faces, truncation after max_events into a reported ledger,
// and batch means with standard errors. Physical scattering kernels (HG) are
// distinct from the isotropic elastic model of radiation_angular.

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/lightcone/ray_geometry.hpp"

namespace cosmosim::radiation {

struct CausalPacket {
  std::uint64_t id = 0;
  std::array<double, 3> origin{};
  std::array<double, 3> direction{};
  double emission_time = 0.0;
  double emitted_energy = 0.0;
  double energy = 0.0;
  double travelled = 0.0;  // path length from origin (comoving/proper as the caller's frame)
};

struct CausalLedger {
  double emitted = 0.0;
  double absorbed = 0.0;
  double escaped = 0.0;
  double unresolved = 0.0;
  [[nodiscard]] double closureResidual(double in_flight) const noexcept {
    return emitted - absorbed - escaped - unresolved - in_flight;
  }
};

class CausalRayTransport {
 public:
  CausalRayTransport(const lightcone::LeafGeometry& geometry, double signal_speed, std::size_t max_packets,
                     double min_energy_fraction, std::size_t max_segments_per_step);
  // Isotropic emission of `energy` at `time` over `directions` equal-energy
  // packets (Fibonacci sphere). Throws (no partial emission) above capacity.
  void emitIsotropic(const std::array<double, 3>& position, double energy, double time, std::size_t directions,
                     std::uint64_t source_id);
  // Moves all packets to t_end; chi indexed by geometry cell row; absorbed
  // energy is added to absorbed_per_cell (same indexing).
  void advance(double t_end, std::span<const double> chi_per_cell, std::span<double> absorbed_per_cell);
  [[nodiscard]] double currentTime() const noexcept { return m_time; }
  [[nodiscard]] const std::vector<CausalPacket>& packets() const noexcept { return m_packets; }
  [[nodiscard]] const CausalLedger& ledger() const noexcept { return m_ledger; }
  [[nodiscard]] double inFlightEnergy() const noexcept;
  [[nodiscard]] std::vector<std::byte> serialize() const;
  void restore(std::span<const std::byte> bytes);

 private:
  const lightcone::LeafGeometry& m_geometry;
  double m_c = 1.0;
  std::size_t m_max_packets = 0;
  double m_min_fraction = 0.0;
  std::size_t m_max_segments = 0;
  double m_time = 0.0;
  std::vector<CausalPacket> m_packets;
  CausalLedger m_ledger;
};

struct HgMonteCarloSettings {
  std::uint64_t seed = 1U;
  std::uint64_t packets = 100000U;
  std::uint64_t batches = 10U;
  std::uint32_t max_events = 100000U;
  double asymmetry_g = 0.0;
  std::size_t max_segments_per_flight = 1U << 16U;
};

struct HgMonteCarloResult {
  std::vector<double> absorbed_per_cell;
  double emitted = 0.0;
  double escaped = 0.0;
  double absorbed = 0.0;
  double truncated = 0.0;
  double escaped_fraction_mean = 0.0;
  double escaped_fraction_standard_error = 0.0;
  std::uint64_t scatterings = 0U;
};

[[nodiscard]] double sampleHenyeyGreensteinCosine(double g, double xi);
// Rotates `direction` by polar cosine mu and azimuth phi about itself.
[[nodiscard]] std::array<double, 3> scatterDirection(const std::array<double, 3>& direction, double mu, double phi);

[[nodiscard]] HgMonteCarloResult runHenyeyGreensteinMonteCarlo(const lightcone::LeafGeometry& geometry,
                                                               std::span<const double> chi_absorption,
                                                               std::span<const double> chi_scattering,
                                                               const std::array<double, 3>& source_position,
                                                               double source_energy,
                                                               const HgMonteCarloSettings& settings);

}  // namespace cosmosim::radiation
