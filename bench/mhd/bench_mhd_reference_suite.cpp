// CHUI-MHD-PERF-v1 benchmark runner (BP01-BP10) for the patch/level kernels.
//
// Emits one JSON document on stdout. Every value is produced by this process
// at run time; nothing is hard-coded. Items that require the production MHD
// dispatch (not wired), external codes (Athena++/AthenaK, IDEFIX, AREPO) or
// hardware counters are reported as "not_evaluated" with a reason — never as
// a pass. Timing protocol (PERF section 14.2): one untimed warm-up, then
// `repeats` independent samples, reporting median and median absolute
// deviation. Source written, NOT EXECUTED in the source-only campaign.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include "../../tests/unit/mhd/mhd_test_fixtures.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/mhd/mhd_initial_conditions.hpp"
#include "cosmosim/mhd/mhd_level_stepper.hpp"
#include "cosmosim/mhd/mhd_memory.hpp"

namespace {

using namespace cosmosim::mhd;
using cosmosim::mhd::test::GasPrimitive;
using cosmosim::mhd::test::TestLevel;
constexpr double k_two_pi = 2.0 * std::numbers::pi;

struct Statistics {
  double median = 0.0;
  double mad = 0.0;
  std::size_t samples = 0;
};

Statistics summarize(std::vector<double> values) {
  Statistics stats{};
  stats.samples = values.size();
  if (values.empty()) return stats;
  std::sort(values.begin(), values.end());
  stats.median = values[values.size() / 2U];
  std::vector<double> deviations;
  for (const double value : values) deviations.push_back(std::abs(value - stats.median));
  std::sort(deviations.begin(), deviations.end());
  stats.mad = deviations[deviations.size() / 2U];
  return stats;
}

// Times `body` (which builds and runs one complete problem) and returns the
// wall-time statistics plus the scalar accuracy metric of the last sample.
Statistics timeRepeated(int repeats, const std::function<double()>& body, double* metric) {
  (void)body();  // warm-up
  std::vector<double> seconds;
  for (int sample = 0; sample < repeats; ++sample) {
    const auto begin = std::chrono::steady_clock::now();
    *metric = body();
    const auto end = std::chrono::steady_clock::now();
    seconds.push_back(std::chrono::duration<double>(end - begin).count());
  }
  return summarize(seconds);
}

std::string entry(const std::string& id, const std::string& status, const std::string& body) {
  return "    {\"id\": \"" + id + "\", \"status\": \"" + status + "\"" +
      (body.empty() ? std::string() : ", " + body) + "}";
}

std::string number(const char* name, double value) {
  std::ostringstream stream;
  stream.precision(10);
  stream << '"' << name << "\": " << value;
  return stream.str();
}

double alfvenProblem(std::int64_t cells, MhdReconstructionMode mode) {
  const double amplitude = 0.1;
  TestLevel level({cells, 4, 4}, {cells, 4, 4}, {1.0 / cells, 0.25, 0.25},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize(
      [&](std::array<double, 3> x) {
        GasPrimitive w{};
        w.pressure = 0.1;
        w.velocity = {0.0, -amplitude * std::sin(k_two_pi * x[0]), -amplitude * std::cos(k_two_pi * x[0])};
        return w;
      },
      {1.0, 0.0, 0.0}, circularAlfvenPotential(amplitude, 1.0), 1.0, k_mhd_gamma_default);
  MhdNumericsPolicy policy{};
  policy.reconstruction = mode;
  MhdLevelStepper stepper = level.makeStepper(policy);
  double time = 0.0;
  while (time < 1.0 - 1e-14) {
    MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(time, 1.0);
    interval.dt_code = std::min(0.8 * stepper.estimateStableTimeStep(interval, {}), 1.0 - time);
    const MhdLevelAdvanceResult result = stepper.advance(interval, {});
    if (!result.accepted) return -1.0;
    time += result.dt_accepted_code;
  }
  double error = 0.0;
  for (std::size_t i = 0; i < static_cast<std::size_t>(cells); ++i) {
    const double x = (static_cast<double>(i) + 0.5) / static_cast<double>(cells);
    error += std::abs(level.faces(0).cellFieldComoving(i, 1, 1)[1] - amplitude * std::sin(k_two_pi * x)) /
        static_cast<double>(cells);
  }
  return error;
}

double tiledStageProblem(std::int64_t patch_cells, std::int64_t tiles, std::uint64_t* scratch) {
  const std::int64_t cells = patch_cells * tiles;
  TestLevel level({cells, cells, patch_cells}, {patch_cells, patch_cells, patch_cells},
                  {1.0 / cells, 1.0 / cells, 1.0 / cells}, cosmosim::mhd::test::periodicDomain());
  level.initialize([](std::array<double, 3> x) {
                     GasPrimitive w{};
                     w.velocity = {std::sin(k_two_pi * x[1]), std::sin(k_two_pi * x[0]), 0.0};
                     return w;
                   },
                   {0.0, 0.0, 0.0}, fieldLoopPotential(1e-3, 0.3, {0.5, 0.5, 0.5}), 1.0,
                   k_mhd_gamma_default);
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(0.0, 1.0);
  interval.dt_code = stepper.estimateStableTimeStep(interval, {});
  const MhdLevelAdvanceResult result = stepper.advance(interval, {});
  *scratch = result.scratch_capacity_bytes;
  return result.divergence_max_normalized;
}

double diffusionProblem(double eta, MhdTransportStageDiagnostics* diagnostics) {
  cosmosim::core::UnitSystem units{};
  units.length_si_per_code = 1.0e-2;
  units.mass_si_per_code = 1.0e-3;
  units.velocity_si_per_code = 1.0e-2;
  TestLevel level({64, 4, 4}, {64, 4, 4}, {1.0 / 64, 0.25, 0.25}, cosmosim::mhd::test::periodicDomain());
  level.initialize([](std::array<double, 3>) { return GasPrimitive{}; }, {0.0, 0.0, 0.0},
                   [](MhdAxis axis, std::array<double, 3> start, double length) {
                     return axis == MhdAxis::kZ ? 1e-3 * std::cos(k_two_pi * start[0]) / k_two_pi * length : 0.0;
                   },
                   1.0, k_mhd_gamma_default);
  ConstantMhdCoefficientProvider provider(eta, 0.0, 0.0);
  MhdPhysicalPropertyModel model{};
  model.units = &units;
  model.mean_molecular_weight = 0.6;
  const MhdNonIdealRateOperator rate(&provider, model, MhdStiffnessLimits{}, true, false, false);
  MhdLevelAdvanceOptions options{};
  options.nonideal = &rate;
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  const MhdLevelAdvanceResult result =
      stepper.advance(cosmosim::mhd::test::isolatedInterval(0.0, 2.0e-3), options);
  *diagnostics = result.transport;
  return result.accepted ? 1.0 : 0.0;
}

}  // namespace

int main(int argc, char** argv) {
  const int repeats = argc > 1 ? std::max(1, std::atoi(argv[1])) : 7;
  std::vector<std::string> entries;

  // BP01: Alfven wave error vs wall time at 32/64/128 cells per wavelength.
  for (const MhdReconstructionMode mode : {MhdReconstructionMode::kPlmMc, MhdReconstructionMode::kPpm}) {
    std::string body = "\"reconstruction\": \"" + std::string(mhdReconstructionName(mode)) + "\", \"runs\": [";
    for (const std::int64_t cells : {32, 64, 128}) {
      double error = 0.0;
      const Statistics stats = timeRepeated(repeats, [&]() { return alfvenProblem(cells, mode); }, &error);
      body += "{" + number("cells_per_wavelength", static_cast<double>(cells)) + ", " +
          number("l1_error_by", error) + ", " + number("wall_seconds_median", stats.median) + ", " +
          number("wall_seconds_mad", stats.mad) + "}" + (cells == 128 ? "" : ", ");
    }
    body += "], \"external_reference\": \"not_evaluated: Athena-family run on matched hardware unavailable\"";
    entries.push_back(entry("BP01", "measured_local_only", body));
  }
  // BP02/BP03: robustness and dissipation problems need production output and
  // matched external references for the time-to-accuracy comparison.
  entries.push_back(entry("BP02", "not_evaluated",
                          "\"reason\": \"Brio-Wu/blast time-to-accuracy requires frozen high-resolution "
                          "reference data and matched external runs\""));
  entries.push_back(entry("BP03", "not_evaluated",
                          "\"reason\": \"Orszag-Tang/field-loop Pareto study requires external codes\""));
  // BP04: tiled 32^3 patches, one interval; divergence and retained scratch.
  {
    std::uint64_t scratch = 0U;
    double divergence = 0.0;
    const Statistics stats =
        timeRepeated(repeats, [&]() { return tiledStageProblem(32, 2, &scratch); }, &divergence);
    entries.push_back(entry("BP04", "measured_local_only",
                            number("patch_cells", 32.0) + ", " + number("tiles_per_axis", 2.0) + ", " +
                                number("wall_seconds_median", stats.median) + ", " +
                                number("divergence_max_normalized", divergence) + ", " +
                                number("scratch_capacity_bytes", static_cast<double>(scratch)) +
                                ", \"amr_reflux\": \"not_evaluated: no production AMR MHD dispatch\""));
  }
  // BP05: explicit SSPRK2 vs RKL2 crossover (stage counts are measured).
  for (const double eta : {1.0e-3, 1.0e-1, 1.0}) {
    MhdTransportStageDiagnostics diagnostics{};
    double accepted = 0.0;
    const Statistics stats = timeRepeated(repeats, [&]() { return diffusionProblem(eta, &diagnostics); }, &accepted);
    entries.push_back(entry(
        "BP05", accepted > 0.0 ? "measured_local_only" : "rejected",
        number("eta_code", eta) + ", " + number("wall_seconds_median", stats.median) + ", " +
            number("explicit_substeps", static_cast<double>(diagnostics.explicit_diffusion_substeps)) + ", " +
            number("rkl_cycles", static_cast<double>(diagnostics.rkl_cycle_count)) + ", " +
            number("rkl_stages", static_cast<double>(diagnostics.rkl_stage_count)) +
            ", \"external_reference\": \"not_evaluated: IDEFIX matched run unavailable\""));
  }
  entries.push_back(entry("BP06", "not_evaluated",
                          "\"reason\": \"Hall cost curve requires a calibrated stiffness sweep campaign; "
                          "kernel exists, scheduler refusal path is unit-tested only\""));
  entries.push_back(entry("BP07", "not_evaluated",
                          "\"reason\": \"no physically validated transport table is distributed\""));
  entries.push_back(entry("BP08", "not_evaluated",
                          "\"reason\": \"cloud collapse needs gravity/source coupling via the production "
                          "dispatch, which is fail-closed\""));
  entries.push_back(entry("BP09", "not_evaluated",
                          "\"reason\": \"cosmological workflow requires the production MHD dispatch\""));
  // BP10: allocation model vs retained capacity at 16/32/64 (exact formulas).
  {
    std::string body = "\"patches\": [";
    for (const std::size_t n : {16U, 32U, 64U}) {
      MhdMemoryRequest request{};
      request.geometry.nx = request.geometry.ny = request.geometry.nz = n;
      request.geometry.dx_comoving = request.geometry.dy_comoving = request.geometry.dz_comoving = 1.0;
      const MhdPatchMemoryEstimate estimate = estimatePatchMemory(request);
      body += "{" + number("n", static_cast<double>(n)) + ", " +
          number("canonical_face_bytes", static_cast<double>(estimate.canonical_face_flux_bytes)) + ", " +
          number("phase_plus_scratch_bytes", static_cast<double>(estimate.phasePlusScratchBytes())) + ", " +
          number("bytes_per_active_cell", scratchBytesPerActiveCell(request, 1U)) + "}" + (n == 64U ? "" : ", ");
    }
    body += "], \"rss\": \"not_evaluated: process RSS sampling belongs to the production profiler\"";
    entries.push_back(entry("BP10", "model_only", body));
  }

  std::cout << "{\n  \"suite\": \"chui_mhd_perf_v1\",\n  \"repeats\": " << repeats
            << ",\n  \"certification\": \"none\",\n  \"benchmarks\": [\n";
  for (std::size_t index = 0; index < entries.size(); ++index) {
    std::cout << entries[index] << (index + 1U == entries.size() ? "\n" : ",\n");
  }
  std::cout << "  ]\n}\n";
  return 0;
}
