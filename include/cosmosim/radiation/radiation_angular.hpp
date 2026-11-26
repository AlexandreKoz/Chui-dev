#pragma once

// CHUI-RT-BH-AGN-REF-v1, P1-RT01: angular backends over the radiation record
// tiling. STATUS: P1 source implementation, NOT qualified; production
// admission refuses physics.rt_transport != m1. Available to explicitly named
// reference API calls and tests.
//
// 1. Direct adaptive rays (point sources). Initial directions: deterministic
//    Fibonacci-sphere quadrature of N0 equal solid-angle rays; a ray splits
//    into 4 children when Omega r^2 > A_cell / q_split (Abel & Wandelt 2002
//    criterion), children offset in the local tangent plane and carrying 1/4
//    of the parent budget. Traversal: exact cell-boundary (Amanatides-Woo DDA)
//    marching inside a record, point location to cross into the next record
//    (AMR leaves traversed once per segment). Each segment absorbs
//    1 - exp(-tau), tau = chi ds, of the ray's remaining energy, deposited in
//    the cell's absorbed ledger; energy reaching the domain boundary is
//    ledgered as escaped. Explicitly a STEADY-STATE (infinite propagation
//    speed) limit; deterministic ray IDs (parent id * 4 + child); bounded
//    queue (max_rays) with backpressure (overflow terminates the ray into an
//    explicit "unresolved" ledger, never dropped silently). Direct photons are
//    absorbed or escape exactly once and are NEVER re-emitted into M1.
// 2. VET closure: normalized Eddington tensor f = P/E from an angular solution
//    (sum_k w_k I_k n_k n_k / sum_k w_k I_k), symmetrized, projected onto the
//    PSD cone (Jacobi eigen-decomposition, negative eigenvalues set to 0) and
//    renormalized to trace 1; the projection size is reported. The M1 HLL
//    Riemann solver is NOT reused with a VET tensor (a distinct transport
//    backend is required; not implemented).
// 3. Monte Carlo reference / post-processing: packets with absorption and
//    isotropic elastic scattering, counter-based SplitMix64 RNG keyed by
//    (seed, packet id, event) for reproducibility, absorbed-energy estimator
//    per cell with per-cell sample variance, escaped fraction; batched.

#include <array>
#include <cstdint>
#include <functional>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/rt_bh_agn_frame.hpp"

namespace cosmosim::radiation {

// Opacity per proper code length of (record, cell); grey or one group.
using RadiationOpacityFunction = std::function<double(std::uint32_t record, std::uint32_t cell)>;

struct RecordCellRef {
  std::uint32_t record = 0U;
  std::uint32_t cell = 0U;
};

class RadiationRecordLocator {
 public:
  RadiationRecordLocator(const core::RadiationMomentSoa& moments, const physics::PeriodicBox& box,
                         std::array<double, 3> domain_origin);
  [[nodiscard]] bool locate(const std::array<double, 3>& position_stored, RecordCellRef& out) const;
  // Distance (stored units) to the exit of the cell along unit direction d.
  [[nodiscard]] double exitDistance(const RecordCellRef& ref, const std::array<double, 3>& position,
                                    const std::array<double, 3>& direction) const;
  [[nodiscard]] double cellArea(const RecordCellRef& ref) const;
  [[nodiscard]] bool insideDomain(const std::array<double, 3>& p) const;
  [[nodiscard]] std::array<double, 3> wrap(const std::array<double, 3>& p) const { return m_box.wrap(p); }

 private:
  const core::RadiationMomentSoa& m_moments;
  physics::PeriodicBox m_box;
  std::array<double, 3> m_origin{};
};

struct RayTracingSettings {
  std::uint32_t initial_rays = 48U;
  double split_quality = 3.0;      // q_split
  std::uint32_t max_rays = 1U << 20U;
  double min_energy_fraction = 1.0e-12;  // termination threshold (energy fraction, ledgered)
};

struct RayTracingResult {
  std::vector<double> absorbed_energy;  // per (record, cell) flattened by record cell order
  double escaped_energy = 0.0;
  double terminated_energy = 0.0;        // below threshold or queue overflow (unresolved ledger)
  std::uint64_t rays_traced = 0U;
  std::uint64_t splits = 0U;
  std::uint64_t queue_overflows = 0U;
};

[[nodiscard]] RayTracingResult traceDirectRays(const RadiationRecordLocator& locator,
                                               const core::RadiationMomentSoa& moments,
                                               const std::array<double, 3>& source_position_stored,
                                               double source_energy, const RadiationOpacityFunction& opacity,
                                               const physics::ProperFrameContext& frame,
                                               const RayTracingSettings& settings);

struct VetTensor {
  std::array<double, 9> f{};      // normalized Eddington tensor (row-major), trace 1, PSD
  double projection_size = 0.0;   // Frobenius norm of the PSD/trace correction
  bool defined = false;
};

// Angular samples: directions (unit), intensities I_k >= 0, weights w_k > 0.
[[nodiscard]] VetTensor vetFromAngularSamples(const std::vector<std::array<double, 3>>& directions,
                                              const std::vector<double>& intensities,
                                              const std::vector<double>& weights);

struct MonteCarloSettings {
  std::uint64_t seed = 0x5eed5eedULL;
  std::uint64_t packets = 100000U;
  std::uint64_t batch = 4096U;
  std::uint32_t max_events_per_packet = 100000U;
  double scattering_albedo = 0.0;  // fraction of interactions that scatter (isotropic, elastic)
};

struct MonteCarloResult {
  std::vector<double> absorbed_energy;
  std::vector<double> absorbed_energy_variance;  // per-cell sample variance of the estimator
  double escaped_energy = 0.0;
  double escaped_fraction_standard_error = 0.0;
  std::uint64_t packets = 0U;
  std::uint64_t truncated_packets = 0U;
};

[[nodiscard]] std::uint64_t radiationCounterRandom(std::uint64_t seed, std::uint64_t stream,
                                                   std::uint64_t counter) noexcept;

[[nodiscard]] MonteCarloResult runMonteCarloReference(const RadiationRecordLocator& locator,
                                                      const core::RadiationMomentSoa& moments,
                                                      const std::array<double, 3>& source_position_stored,
                                                      double source_energy, const RadiationOpacityFunction& opacity,
                                                      const physics::ProperFrameContext& frame,
                                                      const MonteCarloSettings& settings);

}  // namespace cosmosim::radiation
