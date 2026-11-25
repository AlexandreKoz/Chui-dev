// CHUI-MHD-REF-v1 pass-02 contracts of the synchronized leaf-hierarchy
// stepper. NOT RUN in the source-only campaign (UNVALIDATED); see
// docs/architecture/chui_mhd_implementation_pass02_status.md.
//
//  * non-uniform patch extents, nonzero domain origin, batch-size invariance
//  * two-level leaf tiling: exact conservation, strict CT, interface-face
//    consistency, uniform-state invariance across the coarse/fine interface
//  * transactional atomicity: entry mirror mismatch (read-only refusal),
//    external retry (nothing committed), foreign failure in a late stage,
//    memory refusal before allocation, non-mutating estimator, attempted vs
//    accepted counters
//  * two-rank equivalence through an in-process collective transport (the
//    exact exchange protocol of MPI runs), including a zero-patch rank

#include <algorithm>
#include <array>
#include <barrier>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

#include "../support/test_require.hpp"
#include "../unit/mhd/mhd_test_fixtures.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/mhd/mhd_initial_conditions.hpp"
#include "cosmosim/mhd/mhd_level_stepper.hpp"
#include "cosmosim/mhd/mhd_transport.hpp"

namespace {

using cosmosim::mhd::MhdBoundaryKind;
using cosmosim::mhd::MhdConservationBudget;
using cosmosim::mhd::MhdConservedState;
using cosmosim::mhd::MhdError;
using cosmosim::mhd::MhdErrorCode;
using cosmosim::mhd::MhdFaceFluxStorage;
using cosmosim::mhd::MhdGasStateSoa;
using cosmosim::mhd::MhdLevelAdvanceOptions;
using cosmosim::mhd::MhdLevelAdvanceResult;
using cosmosim::mhd::MhdLevelDomain;
using cosmosim::mhd::MhdLevelPatchBinding;
using cosmosim::mhd::MhdLevelStepper;
using cosmosim::mhd::MhdLevelStepperConfig;
using cosmosim::mhd::MhdNumericsPolicy;
using cosmosim::mhd::MhdPatchGeometry;
using cosmosim::mhd::MhdRemotePatchDescriptor;
using cosmosim::mhd::MhdRetryMode;
using cosmosim::mhd::MhdStageContext;
using cosmosim::mhd::test::GasPrimitive;
using cosmosim::mhd::test::GasProfile;

constexpr double k_gamma = 5.0 / 3.0;
constexpr double k_two_pi = 6.283185307179586476925286766559;

struct PatchSpec {
  std::uint64_t id = 0U;
  std::int32_t level = 0;
  std::array<std::int64_t, 3> origin{};
  std::array<std::int64_t, 3> cells{};
  std::int32_t rank = 0;
};

// Explicit leaf hierarchy with caller-owned storage (stable addresses).
class Hierarchy {
 public:
  Hierarchy(MhdLevelDomain domain, std::array<double, 3> base_spacing, std::vector<PatchSpec> specs)
      : m_domain(domain), m_base_spacing(base_spacing), m_specs(std::move(specs)) {
    m_gas.reserve(m_specs.size());
    m_faces.reserve(m_specs.size());
    for (const PatchSpec& spec : m_specs) {
      MhdPatchGeometry geometry{};
      const double factor = std::ldexp(1.0, -(spec.level - domain.base_level));
      geometry.nx = static_cast<std::size_t>(spec.cells[0]);
      geometry.ny = static_cast<std::size_t>(spec.cells[1]);
      geometry.nz = static_cast<std::size_t>(spec.cells[2]);
      geometry.dx_comoving = base_spacing[0] * factor;
      geometry.dy_comoving = base_spacing[1] * factor;
      geometry.dz_comoving = base_spacing[2] * factor;
      geometry.origin_x_comoving =
          domain.origin_comoving[0] + static_cast<double>(spec.origin[0]) * geometry.dx_comoving;
      geometry.origin_y_comoving =
          domain.origin_comoving[1] + static_cast<double>(spec.origin[1]) * geometry.dy_comoving;
      geometry.origin_z_comoving =
          domain.origin_comoving[2] + static_cast<double>(spec.origin[2]) * geometry.dz_comoving;
      m_gas.emplace_back(geometry.cellCount());
      m_faces.emplace_back(geometry);
    }
  }

  void initialize(const GasProfile& profile, std::array<double, 3> uniform_b,
                  const cosmosim::mhd::MhdEdgePotential& potential) {
    for (std::size_t index = 0; index < m_specs.size(); ++index) {
      MhdFaceFluxStorage& faces = m_faces[index];
      const double divergence = cosmosim::mhd::seedFaceFlux(
          faces.geometry(), uniform_b, potential, faces.fluxX(), faces.fluxY(), faces.fluxZ());
      CHUI_REQUIRE(divergence < 1.0e-12);
      const MhdPatchGeometry& g = faces.geometry();
      for (std::size_t k = 0; k < g.nz; ++k) {
        for (std::size_t j = 0; j < g.ny; ++j) {
          for (std::size_t i = 0; i < g.nx; ++i) {
            const std::array<double, 3> x{
                g.origin_x_comoving + (static_cast<double>(i) + 0.5) * g.dx_comoving,
                g.origin_y_comoving + (static_cast<double>(j) + 0.5) * g.dy_comoving,
                g.origin_z_comoving + (static_cast<double>(k) + 0.5) * g.dz_comoving};
            const GasPrimitive w = profile(x);
            const std::array<double, 3> b = faces.cellFieldComoving(i, j, k);
            MhdConservedState state{};
            state.rho_comoving = w.rho;
            state.momentum_density_x_comoving = w.rho * w.velocity[0];
            state.momentum_density_y_comoving = w.rho * w.velocity[1];
            state.momentum_density_z_comoving = w.rho * w.velocity[2];
            state.total_energy_density_comoving = w.pressure / (k_gamma - 1.0) +
                0.5 * w.rho *
                    (w.velocity[0] * w.velocity[0] + w.velocity[1] * w.velocity[1] +
                     w.velocity[2] * w.velocity[2]) +
                0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
            state.metal_mass_density_comoving = w.rho * w.metal_fraction;
            m_gas[index].storeGas(g.cellIndex(i, j, k), state);
          }
        }
      }
    }
  }

  [[nodiscard]] std::vector<MhdLevelPatchBinding> bindings(std::int32_t rank) {
    std::vector<MhdLevelPatchBinding> out;
    for (std::size_t index = 0; index < m_specs.size(); ++index) {
      const PatchSpec& spec = m_specs[index];
      if (spec.rank != rank) {
        continue;
      }
      MhdLevelPatchBinding binding{};
      binding.patch_id = spec.id;
      binding.geometry_generation = 1U;
      binding.level = spec.level;
      binding.origin_cell = spec.origin;
      binding.geometry = m_faces[index].geometry();
      binding.gas = &m_gas[index];
      binding.faces = cosmosim::mhd::makeMutableFaceFluxView(m_faces[index], spec.id, 1U, spec.level);
      binding.owner_rank = rank;
      out.push_back(binding);
    }
    return out;
  }

  [[nodiscard]] std::vector<MhdRemotePatchDescriptor> remote(std::int32_t rank) const {
    std::vector<MhdRemotePatchDescriptor> out;
    for (const PatchSpec& spec : m_specs) {
      if (spec.rank == rank) {
        continue;
      }
      out.push_back({spec.id, 1U, spec.level, spec.origin, spec.cells, spec.rank});
    }
    return out;
  }

  // Exact sums over all leaves (independent of the stepper).
  [[nodiscard]] std::array<double, 5> totals() const {
    std::array<double, 5> sum{};
    for (std::size_t index = 0; index < m_specs.size(); ++index) {
      const MhdPatchGeometry& g = m_faces[index].geometry();
      const double volume = g.cellVolume();
      for (std::size_t cell = 0; cell < g.cellCount(); ++cell) {
        const MhdConservedState state = m_gas[index].loadGas(cell);
        sum[0] += state.rho_comoving * volume;
        sum[1] += state.momentum_density_x_comoving * volume;
        sum[2] += state.momentum_density_y_comoving * volume;
        sum[3] += state.momentum_density_z_comoving * volume;
        sum[4] += state.total_energy_density_comoving * volume;
      }
    }
    return sum;
  }

  [[nodiscard]] std::size_t size() const noexcept { return m_specs.size(); }
  [[nodiscard]] MhdGasStateSoa& gas(std::size_t index) { return m_gas[index]; }
  [[nodiscard]] MhdFaceFluxStorage& faces(std::size_t index) { return m_faces[index]; }
  [[nodiscard]] const MhdLevelDomain& domain() const noexcept { return m_domain; }

  // Bitwise snapshot of every lane (for atomicity checks).
  [[nodiscard]] std::vector<double> snapshot() const {
    std::vector<double> out;
    for (std::size_t index = 0; index < m_specs.size(); ++index) {
      for (std::size_t lane = 0; lane < MhdGasStateSoa::k_lane_count; ++lane) {
        const auto values = m_gas[index].lane(lane);
        out.insert(out.end(), values.begin(), values.end());
      }
      for (const auto values : {m_faces[index].fluxX(), m_faces[index].fluxY(), m_faces[index].fluxZ()}) {
        out.insert(out.end(), values.begin(), values.end());
      }
    }
    return out;
  }

 private:
  MhdLevelDomain m_domain;
  std::array<double, 3> m_base_spacing;
  std::vector<PatchSpec> m_specs;
  std::vector<MhdGasStateSoa> m_gas;
  std::vector<MhdFaceFluxStorage> m_faces;
};

GasPrimitive smoothProfile(std::array<double, 3> x) {
  GasPrimitive w{};
  w.rho = 1.0 + 0.2 * std::sin(k_two_pi * x[0]) * std::cos(k_two_pi * x[1]);
  w.velocity = {0.3, -0.1 * std::sin(k_two_pi * x[2]), 0.05};
  w.pressure = 0.6 + 0.1 * std::cos(k_two_pi * x[0]);
  w.metal_fraction = 0.01;
  return w;
}

MhdStageContext isolated(double t, double dt) {
  return cosmosim::mhd::test::isolatedInterval(t, dt);
}

double maxDifference(const std::vector<double>& lhs, const std::vector<double>& rhs) {
  CHUI_REQUIRE(lhs.size() == rhs.size());
  double maximum = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    maximum = std::max(maximum, std::abs(lhs[i] - rhs[i]));
  }
  return maximum;
}

void runSteps(MhdLevelStepper& stepper, int steps, double dt) {
  double t = 0.0;
  for (int step = 0; step < steps; ++step) {
    const MhdLevelAdvanceResult result = stepper.advance(isolated(t, dt), {});
    CHUI_REQUIRE(result.accepted);
    CHUI_REQUIRE(result.divergence_max_normalized <= 1.0e-11);
    t += dt;
  }
}

// --- Non-uniform extents, nonzero origin, batch invariance -----------------

void testNonUniformExtentsAndOrigin() {
  MhdLevelDomain domain = cosmosim::mhd::test::periodicDomain();
  domain.cells = {16, 8, 8};
  const std::array<double, 3> spacing{1.0 / 16, 1.0 / 8, 1.0 / 8};
  Hierarchy single(domain, spacing, {{7, 0, {0, 0, 0}, {16, 8, 8}}});
  Hierarchy split(domain, spacing,
                  {{31, 0, {0, 0, 0}, {6, 8, 8}}, {5, 0, {6, 0, 0}, {4, 8, 8}},
                   {19, 0, {10, 0, 0}, {6, 8, 8}}});
  MhdLevelDomain shifted_domain = domain;
  shifted_domain.origin_comoving = {10.0, -3.0, 5.0};
  Hierarchy shifted(shifted_domain, spacing,
                    {{31, 0, {0, 0, 0}, {6, 8, 8}}, {5, 0, {6, 0, 0}, {4, 8, 8}},
                     {19, 0, {10, 0, 0}, {6, 8, 8}}});
  const auto potential = cosmosim::mhd::fieldLoopPotential(1.0e-3, 0.3, {0.5, 0.5, 0.5});
  single.initialize(smoothProfile, {0.05, 0.02, 0.0}, potential);
  split.initialize(smoothProfile, {0.05, 0.02, 0.0}, potential);
  // The shifted hierarchy carries identical cell data (profile evaluated at the
  // unshifted positions) so only the origin bookkeeping differs.
  shifted.initialize([](std::array<double, 3> x) {
                       return smoothProfile({x[0] - 10.0, x[1] + 3.0, x[2] - 5.0});
                     },
                     {0.05, 0.02, 0.0},
                     cosmosim::mhd::fieldLoopPotential(1.0e-3, 0.3, {10.5, -2.5, 5.5}));
  MhdLevelStepper a(MhdNumericsPolicy{}, domain, single.bindings(0));
  MhdLevelStepper b(MhdNumericsPolicy{}, domain, split.bindings(0));
  MhdLevelStepperConfig one_slot{};
  one_slot.batch_patches = 1U;
  MhdLevelStepper c(MhdNumericsPolicy{}, shifted_domain, shifted.bindings(0), {}, one_slot);
  const double dt = 0.5 * a.estimateStableTimeStep(isolated(0.0, 1.0), {});
  runSteps(a, 3, dt);
  runSteps(b, 3, dt);
  runSteps(c, 3, dt);
  const std::array<double, 5> ta = single.totals();
  const std::array<double, 5> tb = split.totals();
  for (std::size_t q = 0; q < 5U; ++q) {
    CHUI_REQUIRE(std::abs(ta[q] - tb[q]) <= 1.0e-12 * std::max(1.0, std::abs(ta[q])));
  }
  // Split vs shifted-origin split with one batch slot: identical arithmetic.
  CHUI_REQUIRE(maxDifference(split.snapshot(), shifted.snapshot()) <= 1.0e-12);
}

// --- Two-level leaf tiling ---------------------------------------------------

Hierarchy makeTwoLevel(std::int32_t fine_rank = 0, std::int32_t coarse_rank = 0) {
  MhdLevelDomain domain = cosmosim::mhd::test::periodicDomain();
  domain.cells = {16, 8, 8};
  // Level 0 leaves x in [0,4) and [12,16); level 1 leaf covers base x in [4,12).
  return Hierarchy(domain, {1.0 / 16, 1.0 / 8, 1.0 / 8},
                   {{40, 0, {0, 0, 0}, {4, 8, 8}, coarse_rank},
                    {11, 0, {12, 0, 0}, {4, 8, 8}, coarse_rank},
                    {23, 1, {8, 0, 0}, {16, 16, 16}, fine_rank}});
}

void testTwoLevelConservationAndCt() {
  Hierarchy h = makeTwoLevel();
  // Uniform field (exact interface consistency) plus a gas pulse crossing
  // the coarse/fine interfaces.
  h.initialize(smoothProfile, {0.1, 0.05, 0.02}, {});
  MhdLevelStepper stepper(MhdNumericsPolicy{}, h.domain(), h.bindings(0));
  CHUI_REQUIRE(stepper.levelCount() == 2);
  CHUI_REQUIRE(stepper.coarseFineInterfaceFaceCount() > 0U);
  CHUI_REQUIRE(stepper.validateCanonicalConsistency() <= 1.0e-12);
  const std::array<double, 5> begin = h.totals();
  const double dt = 0.5 * stepper.estimateStableTimeStep(isolated(0.0, 1.0), {});
  runSteps(stepper, 4, dt);
  const std::array<double, 5> end = h.totals();
  for (std::size_t q = 0; q < 5U; ++q) {
    // G04: closed periodic conservation (no sources, no floors).
    CHUI_REQUIRE(std::abs(end[q] - begin[q]) <= 1.0e-12 * std::max(1.0, std::abs(begin[q])));
  }
  // Coarse interface faces equal the sum of the fine faces after every step.
  CHUI_REQUIRE(stepper.validateCanonicalConsistency() <= 1.0e-12);
}

void testTwoLevelUniformInvariance() {
  Hierarchy h = makeTwoLevel();
  h.initialize(
      [](std::array<double, 3>) {
        GasPrimitive w{};
        w.rho = 1.3;
        w.velocity = {0.2, -0.1, 0.05};
        w.pressure = 0.7;
        w.metal_fraction = 0.02;
        return w;
      },
      {0.1, -0.2, 0.3}, {});
  const std::vector<double> before = h.snapshot();
  MhdLevelStepper stepper(MhdNumericsPolicy{}, h.domain(), h.bindings(0));
  const double dt = 0.5 * stepper.estimateStableTimeStep(isolated(0.0, 1.0), {});
  runSteps(stepper, 3, dt);
  // Restriction, prolongation and the interface overrides reproduce a
  // uniform state exactly up to roundoff.
  CHUI_REQUIRE(maxDifference(before, h.snapshot()) <= 1.0e-12);
}

// --- Transactions -------------------------------------------------------------

void testEntryMirrorMismatchIsReadOnly() {
  MhdLevelDomain domain = cosmosim::mhd::test::periodicDomain();
  domain.cells = {16, 8, 8};
  Hierarchy h(domain, {1.0 / 16, 1.0 / 8, 1.0 / 8},
              {{3, 0, {0, 0, 0}, {8, 8, 8}}, {9, 0, {8, 0, 0}, {8, 8, 8}}});
  h.initialize(smoothProfile, {0.1, 0.0, 0.0}, {});
  // Corrupt one shared x-face mirror (lower plane of patch 9 = upper of 3).
  h.faces(1).fluxX()[0] += 1.0e-6;
  const std::vector<double> before = h.snapshot();
  MhdLevelStepper stepper(MhdNumericsPolicy{}, h.domain(), h.bindings(0));
  bool refused = false;
  try {
    (void)stepper.advance(isolated(0.0, 1.0e-4), {});
  } catch (const MhdError& error) {
    refused = error.code() == MhdErrorCode::kInvalidField;
  }
  CHUI_REQUIRE(refused);
  // Validation is read-only: the mismatch was neither overwritten nor rolled.
  CHUI_REQUIRE(maxDifference(before, h.snapshot()) == 0.0);
  // The estimator does not normalize mirrors either.
  (void)stepper.estimateStableTimeStep(isolated(0.0, 1.0), {});
  CHUI_REQUIRE(maxDifference(before, h.snapshot()) == 0.0);
}

void testExternalRetryCommitsNothing() {
  Hierarchy h = makeTwoLevel();
  h.initialize(smoothProfile, {0.1, 0.05, 0.02}, {});
  const std::vector<double> before = h.snapshot();
  MhdLevelStepper stepper(MhdNumericsPolicy{}, h.domain(), h.bindings(0));
  const double cfl_dt = stepper.estimateStableTimeStep(isolated(0.0, 1.0), {});
  MhdLevelAdvanceResult result{};
  bool threw = false;
  try {
    result = stepper.advance(isolated(0.0, 200.0 * cfl_dt), {});
  } catch (const MhdError& error) {
    threw = true;
    CHUI_REQUIRE(error.code() == MhdErrorCode::kNonPositiveState);
  }
  if (!threw) {
    CHUI_REQUIRE(!result.accepted);
    CHUI_REQUIRE(result.dt_suggested_code == 0.5 * result.dt_requested_code);
    CHUI_REQUIRE(result.rejected_attempts == 1U);
    // Accepted-only physics counters stay empty; attempted counters record cost.
    CHUI_REQUIRE(result.ideal.faces_solved == 0U);
    CHUI_REQUIRE(result.ideal_attempted.faces_solved > 0U ||
                 result.ideal_stage_evaluations_attempted > 0U);
  }
  CHUI_REQUIRE(maxDifference(before, h.snapshot()) == 0.0);
}

// Provider that throws a foreign (non-MHD) exception after `limit` queries:
// a late failure after successful earlier stages.
class FailingProvider final : public cosmosim::mhd::MhdCoefficientProvider {
 public:
  explicit FailingProvider(std::size_t limit) : m_limit(limit) {}
  [[nodiscard]] cosmosim::mhd::MhdCoefficientProviderKind kind() const noexcept override {
    return cosmosim::mhd::MhdCoefficientProviderKind::kConstant;
  }
  [[nodiscard]] std::string provenanceId() const override { return "test_failing_provider"; }
  void validateAdmission() const override {}
  [[nodiscard]] cosmosim::mhd::MhdTransportCoefficients evaluate(
      std::size_t, const cosmosim::mhd::MhdCellPhysicalProperties&) const override {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (++m_calls > m_limit) {
      throw std::runtime_error("injected provider failure");
    }
    cosmosim::mhd::MhdTransportCoefficients c{};
    c.eta_ohmic_phys_cgs = 1.0e-3;
    c.valid = true;
    return c;
  }

 private:
  std::size_t m_limit;
  mutable std::size_t m_calls = 0U;
  mutable std::mutex m_mutex;
};

cosmosim::core::UnitSystem cgsUnits() {
  cosmosim::core::UnitSystem units{};
  units.length_unit = "cm";
  units.mass_unit = "g";
  units.velocity_unit = "cm/s";
  units.length_si_per_code = 1.0e-2;
  units.mass_si_per_code = 1.0e-3;
  units.velocity_si_per_code = 1.0e-2;
  return units;
}

void testForeignLateFailureRollsBack() {
  MhdLevelDomain domain = cosmosim::mhd::test::periodicDomain();
  domain.cells = {16, 8, 8};
  Hierarchy h(domain, {1.0 / 16, 1.0 / 8, 1.0 / 8},
              {{3, 0, {0, 0, 0}, {8, 8, 8}}, {9, 0, {8, 0, 0}, {8, 8, 8}}});
  h.initialize(smoothProfile, {0.1, 0.05, 0.0},
               cosmosim::mhd::circularAlfvenPotential(0.01, 1.0));
  const std::vector<double> before = h.snapshot();
  const cosmosim::core::UnitSystem units = cgsUnits();
  cosmosim::mhd::MhdPhysicalPropertyModel properties{};
  properties.units = &units;
  properties.mean_molecular_weight = 0.6;
  // Enough queries for the forecast and the first D half, not for the rest.
  FailingProvider provider(3U * 10U * 10U * 10U * 2U);
  const cosmosim::mhd::MhdNonIdealRateOperator rate(&provider, properties,
                                                    cosmosim::mhd::MhdStiffnessLimits{}, true,
                                                    false, false);
  MhdLevelAdvanceOptions options{};
  options.nonideal = &rate;
  MhdLevelStepper stepper(MhdNumericsPolicy{}, h.domain(), h.bindings(0));
  bool refused = false;
  try {
    (void)stepper.advance(isolated(0.0, 1.0e-4), options);
  } catch (const MhdError& error) {
    // A foreign exception is agreed as a non-retryable MHD failure (or the
    // provider adapter maps it to MHD_BAD_COEFFICIENT); never retried.
    refused = error.code() == MhdErrorCode::kInvalidField ||
        error.code() == MhdErrorCode::kBadCoefficient;
  }
  CHUI_REQUIRE(refused);
  CHUI_REQUIRE(maxDifference(before, h.snapshot()) == 0.0);
}

void testMemoryRefusalBeforeAllocation() {
  Hierarchy h = makeTwoLevel();
  h.initialize(smoothProfile, {0.1, 0.0, 0.0}, {});
  MhdLevelStepperConfig config{};
  config.require_memory_admission = true;
  MhdLevelStepper stepper(MhdNumericsPolicy{}, h.domain(), h.bindings(0), {}, config);
  cosmosim::core::MemoryGovernorPolicy policy{};
  policy.hard_limit_bytes = 4096U;
  cosmosim::core::MemoryGovernor governor(policy);
  bool refused = false;
  try {
    stepper.admitMemory(governor, {});
  } catch (const MhdError& error) {
    refused = error.code() == MhdErrorCode::kMemoryBudgetExceeded;
  }
  CHUI_REQUIRE(refused);
  CHUI_REQUIRE(stepper.ownedCapacityBytes() == 0U);
  // Without admission a production-configured stepper refuses to allocate.
  refused = false;
  try {
    (void)stepper.advance(isolated(0.0, 1.0e-4), {});
  } catch (const MhdError& error) {
    refused = error.code() == MhdErrorCode::kMemoryBudgetExceeded;
  }
  CHUI_REQUIRE(refused);
  CHUI_REQUIRE(stepper.ownedCapacityBytes() == 0U);
  // The honest plan separates level-resident state from batch scratch.
  const cosmosim::mhd::MhdLevelMemoryPlan plan = stepper.planMemory({});
  CHUI_REQUIRE(plan.level_resident_bytes > 0U && plan.batch_scratch_bytes > 0U);
  CHUI_REQUIRE(plan.active_cells == 2U * 4U * 8U * 8U + 16U * 16U * 16U);
}

// --- Two-rank equivalence through an in-process collective transport -------

class InProcessHub {
 public:
  explicit InProcessHub(int size) : m_size(size), m_barrier(size), m_outbox(size), m_u64(size),
                                    m_double(size), m_sums(size) {}
  int m_size;
  std::barrier<> m_barrier;
  std::vector<std::vector<std::vector<std::uint8_t>>> m_outbox;
  std::vector<std::uint64_t> m_u64;
  std::vector<double> m_double;
  std::vector<std::vector<double>> m_sums;
};

class InProcessTransport final : public cosmosim::mhd::MhdHaloTransport {
 public:
  InProcessTransport(InProcessHub& hub, int rank) : m_hub(hub), m_rank(rank) {}
  [[nodiscard]] int rank() const noexcept override { return m_rank; }
  [[nodiscard]] int size() const noexcept override { return m_hub.m_size; }
  [[nodiscard]] std::vector<std::vector<std::uint8_t>> exchange(
      const std::vector<std::vector<std::uint8_t>>& send) override {
    m_hub.m_outbox[static_cast<std::size_t>(m_rank)] = send;
    m_hub.m_barrier.arrive_and_wait();
    std::vector<std::vector<std::uint8_t>> received(static_cast<std::size_t>(m_hub.m_size));
    for (int r = 0; r < m_hub.m_size; ++r) {
      if (r != m_rank) {
        received[static_cast<std::size_t>(r)] =
            m_hub.m_outbox[static_cast<std::size_t>(r)][static_cast<std::size_t>(m_rank)];
      }
    }
    m_hub.m_barrier.arrive_and_wait();
    return received;
  }
  [[nodiscard]] std::uint64_t allreduceMaxUint64(std::uint64_t value) override {
    m_hub.m_u64[static_cast<std::size_t>(m_rank)] = value;
    m_hub.m_barrier.arrive_and_wait();
    const std::uint64_t result = *std::max_element(m_hub.m_u64.begin(), m_hub.m_u64.end());
    m_hub.m_barrier.arrive_and_wait();
    return result;
  }
  [[nodiscard]] double allreduceMinDouble(double value) override {
    m_hub.m_double[static_cast<std::size_t>(m_rank)] = value;
    m_hub.m_barrier.arrive_and_wait();
    const double result = *std::min_element(m_hub.m_double.begin(), m_hub.m_double.end());
    m_hub.m_barrier.arrive_and_wait();
    return result;
  }
  [[nodiscard]] double allreduceMaxDouble(double value) override {
    m_hub.m_double[static_cast<std::size_t>(m_rank)] = value;
    m_hub.m_barrier.arrive_and_wait();
    const double result = *std::max_element(m_hub.m_double.begin(), m_hub.m_double.end());
    m_hub.m_barrier.arrive_and_wait();
    return result;
  }
  void allreduceSumDoubles(std::span<double> values) override {
    m_hub.m_sums[static_cast<std::size_t>(m_rank)].assign(values.begin(), values.end());
    m_hub.m_barrier.arrive_and_wait();
    for (std::size_t i = 0; i < values.size(); ++i) {
      double sum = 0.0;
      for (int r = 0; r < m_hub.m_size; ++r) {
        sum += m_hub.m_sums[static_cast<std::size_t>(r)][i];  // rank order: deterministic
      }
      values[i] = sum;
    }
    m_hub.m_barrier.arrive_and_wait();
  }

 private:
  InProcessHub& m_hub;
  int m_rank;
};

void runTwoRanks(Hierarchy& h, double dt, int steps) {
  InProcessHub hub(2);
  std::array<std::exception_ptr, 2> failures{};
  std::array<std::thread, 2> threads;
  for (int rank = 0; rank < 2; ++rank) {
    threads[static_cast<std::size_t>(rank)] = std::thread([&, rank]() {
      try {
        InProcessTransport transport(hub, rank);
        MhdLevelStepperConfig config{};
        config.transport = &transport;
        MhdLevelStepper stepper(MhdNumericsPolicy{}, h.domain(), h.bindings(rank), h.remote(rank),
                                config);
        double t = 0.0;
        for (int step = 0; step < steps; ++step) {
          const MhdLevelAdvanceResult result = stepper.advance(isolated(t, dt), {});
          CHUI_REQUIRE(result.accepted);
          t += dt;
        }
      } catch (...) {
        failures[static_cast<std::size_t>(rank)] = std::current_exception();
      }
    });
  }
  for (std::thread& thread : threads) {
    thread.join();
  }
  for (const std::exception_ptr& failure : failures) {
    if (failure) {
      std::rethrow_exception(failure);
    }
  }
}

void testTwoRankEquivalence() {
  Hierarchy serial = makeTwoLevel();
  Hierarchy split = makeTwoLevel(/*fine_rank=*/1, /*coarse_rank=*/0);
  Hierarchy idle = makeTwoLevel(/*fine_rank=*/0, /*coarse_rank=*/0);  // rank 1 owns nothing
  for (Hierarchy* h : {&serial, &split, &idle}) {
    h->initialize(smoothProfile, {0.1, 0.05, 0.02}, {});
  }
  MhdLevelStepper reference(MhdNumericsPolicy{}, serial.domain(), serial.bindings(0));
  const double dt = 0.5 * reference.estimateStableTimeStep(isolated(0.0, 1.0), {});
  runSteps(reference, 2, dt);
  runTwoRanks(split, dt, 2);
  runTwoRanks(idle, dt, 2);
  // Same plan, same owners, same arithmetic order: bitwise agreement expected;
  // the gate tolerance is roundoff.
  CHUI_REQUIRE(maxDifference(serial.snapshot(), split.snapshot()) <= 1.0e-13);
  CHUI_REQUIRE(maxDifference(serial.snapshot(), idle.snapshot()) <= 1.0e-13);
}

}  // namespace

int main() {
  testNonUniformExtentsAndOrigin();
  testTwoLevelConservationAndCt();
  testTwoLevelUniformInvariance();
  testEntryMirrorMismatchIsReadOnly();
  testExternalRetryCommitsNothing();
  testForeignLateFailureRollsBack();
  testMemoryRefusalBeforeAllocation();
  testTwoRankEquivalence();
  return 0;
}
