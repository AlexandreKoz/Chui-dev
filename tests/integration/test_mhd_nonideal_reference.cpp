// CHUI-MHD-REF-v1 non-ideal reference integration tests (N01, N03, N04, N05,
// N06/N07 smoke, RKL2 recurrence, coefficient units, stiffness caps).
// Source written, NOT EXECUTED in the source-only campaign.

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <numbers>
#include <string>
#include <vector>

#include "../unit/mhd/mhd_test_fixtures.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/mhd/mhd_initial_conditions.hpp"
#include "cosmosim/mhd/mhd_level_stepper.hpp"
#include "cosmosim/mhd/mhd_transport.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace {

using namespace cosmosim::mhd;
using cosmosim::mhd::test::GasPrimitive;
using cosmosim::mhd::test::TestLevel;
constexpr double k_gamma = 5.0 / 3.0;
constexpr double k_two_pi = 2.0 * std::numbers::pi;

// Code units = cgs (L = 1 cm, M = 1 g, V = 1 cm/s => T = 1 s), so
// eta_code = eta_cgs / a^2 exactly.
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

MhdPhysicalPropertyModel properties(const cosmosim::core::UnitSystem& units) {
  MhdPhysicalPropertyModel model{};
  model.units = &units;
  model.mean_molecular_weight = 0.6;
  return model;
}

// RKL2 applied to the scalar test problem y' = -lambda y. Stability at the
// edge of the interval tau <= dt_euler (s^2+s-2)/4 and second order accuracy.
double rkl2Scalar(std::int32_t stages, double lambda_tau) {
  const MhdRkl2Coefficients c = makeRkl2Coefficients(stages);
  const auto L = [lambda_tau](double y) { return -lambda_tau * y; };  // tau L(y)
  const double y0 = 1.0;
  const double ly0 = L(y0);
  double yjm2 = y0;
  double yjm1 = y0 + c.mu_tilde[1] * ly0;
  for (std::int32_t j = 2; j <= stages; ++j) {
    const auto n = static_cast<std::size_t>(j);
    const double yj = c.mu[n] * yjm1 + c.nu[n] * yjm2 + (1.0 - c.mu[n] - c.nu[n]) * y0 +
        c.mu_tilde[n] * L(yjm1) + c.gamma_tilde[n] * ly0;
    yjm2 = yjm1;
    yjm1 = yj;
  }
  return yjm1;
}

void testRkl2Recurrence() {
  for (const std::int32_t s : {2, 3, 5, 9, 17, 64, 128}) {
    // Forward-Euler limit for y' = -lambda y is lambda dt = 2.
    const double edge = 2.0 * rkl2StabilityFactor(s) * (1.0 - 1.0e-9);
    for (int sample = 1; sample <= 200; ++sample) {
      const double lambda_tau = edge * sample / 200.0;
      CHUI_REQUIRE(std::abs(rkl2Scalar(s, lambda_tau)) <= 1.0 + 1.0e-12);
    }
    // Second order: local error O((lambda tau)^3).
    const double e1 = std::abs(rkl2Scalar(s, 1.0e-2) - std::exp(-1.0e-2));
    const double e2 = std::abs(rkl2Scalar(s, 5.0e-3) - std::exp(-5.0e-3));
    CHUI_REQUIRE(e1 / e2 > 7.0);  // ~8 for a third-order local error
  }
  CHUI_REQUIRE(selectRkl2StageCount(0.5, 128) == 0);
  for (const double ratio : {1.5, 4.0, 37.0, 1000.0}) {
    const std::int32_t s = selectRkl2StageCount(ratio, 128);
    CHUI_REQUIRE(s >= 2 && rkl2StabilityFactor(s) >= ratio && rkl2StabilityFactor(s - 1) < ratio);
  }
  const MhdRkl2Plan plan = planRkl2(1.0e5, 128);
  CHUI_REQUIRE(plan.stages_per_cycle <= 128 && plan.cycles >= 25);
  CHUI_REQUIRE(static_cast<double>(plan.cycles) * rkl2StabilityFactor(plan.stages_per_cycle) >= 1.0e5);
}

void testCoefficientUnitConversion() {
  cosmosim::core::UnitSystem units = cgsUnits();
  CHUI_REQUIRE(std::abs(physicalDiffusivityToComovingCode(3.0, 1.0, units) - 3.0) <= 1e-15);
  CHUI_REQUIRE(std::abs(physicalDiffusivityToComovingCode(3.0, 0.5, units) - 12.0) <= 1e-14);
  // kpc / km/s units: T/L^2 = 1 / (L V).
  units.length_si_per_code = 3.0856775814913673e19;
  units.velocity_si_per_code = 1.0e3;
  const double expected = 1.0e26 / (3.0856775814913673e21 * 1.0e5);
  CHUI_REQUIRE(std::abs(physicalDiffusivityToComovingCode(1.0e26, 1.0, units) - expected) <=
               1e-12 * expected);
}

struct Mode {
  std::complex<double> by;
  std::complex<double> bz;
};

Mode fieldMode(TestLevel& level, std::int64_t cells) {
  Mode mode{};
  const MhdFaceFluxStorage& faces = level.faces(0);
  for (std::int64_t i = 0; i < cells; ++i) {
    const double x = (static_cast<double>(i) + 0.5) / static_cast<double>(cells);
    const std::array<double, 3> b = faces.cellFieldComoving(static_cast<std::size_t>(i), 1, 1);
    const std::complex<double> phase = std::exp(std::complex<double>(0.0, -k_two_pi * x));
    mode.by += b[1] * phase / static_cast<double>(cells);
    mode.bz += b[2] * phase / static_cast<double>(cells);
  }
  return mode;
}

// N01: Ohmic decay of B_y = B0 sin(k x) at 128 cells per wavelength; decay
// rate within 1%; closed periodic total energy conserved, thermal gain equal
// to magnetic loss.
void testOhmicDecay() {
  const std::int64_t cells = 128;
  const cosmosim::core::UnitSystem units = cgsUnits();
  const double eta = 0.05;
  const double b0 = 1.0e-3;
  TestLevel level({cells, 4, 4}, {cells, 4, 4}, {1.0 / cells, 0.25, 0.25},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize([](std::array<double, 3>) { return GasPrimitive{}; }, {0.0, 0.0, 0.0},
                   [b0](MhdAxis axis, std::array<double, 3> start, double length) {
                     // A_z = B0 cos(kx)/k  =>  B_y = -dA_z/dx = B0 sin(kx).
                     return axis == MhdAxis::kZ ? b0 * std::cos(k_two_pi * start[0]) / k_two_pi * length
                                                : 0.0;
                   },
                   1.0, k_gamma);
  ConstantMhdCoefficientProvider provider(eta, 0.0, 0.0);
  provider.validateAdmission();
  const MhdNonIdealRateOperator rate(&provider, properties(units), MhdStiffnessLimits{}, true,
                                     false, false);
  MhdLevelAdvanceOptions options{};
  options.nonideal = &rate;
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  const MhdConservationBudget begin = stepper.levelBudget(1.0);
  const double initial = std::abs(fieldMode(level, cells).by);
  const double end_time = 0.5;
  double time = 0.0;
  while (time < end_time * (1.0 - 1.0e-14)) {
    MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(time, 1.0);
    // Ideal CFL: 0.3 / (c_s (1/dx + 1/dy + 1/dz)) = 0.3 / (1.29 * 136) ~ 1.7e-3.
    interval.dt_code = std::min(1.5e-3, end_time - time);
    const MhdLevelAdvanceResult result = stepper.advance(interval, options);
    CHUI_REQUIRE(result.accepted);
    time += result.dt_accepted_code;
  }
  const double measured_rate = -std::log(std::abs(fieldMode(level, cells).by) / initial) / end_time;
  const double expected_rate = eta * k_two_pi * k_two_pi;
  CHUI_REQUIRE(std::abs(measured_rate - expected_rate) <= 0.01 * expected_rate);
  const MhdConservationBudget end = stepper.levelBudget(1.0);
  CHUI_REQUIRE(std::abs(end.energy_total - begin.energy_total) <= 1e-12 * begin.energy_total);
  const double magnetic_loss = begin.energy_magnetic - end.energy_magnetic;
  const double thermal_gain = end.energy_thermal - begin.energy_thermal;
  CHUI_REQUIRE(magnetic_loss > 0.0);
  CHUI_REQUIRE(std::abs(thermal_gain + (end.energy_kinetic - begin.energy_kinetic) - magnetic_loss) <=
               1e-9 * begin.energy_total);
}

// N03/N04: Hall whistler on a weak guide field (Alfven frequency negligible):
// B_y + i B_z rotates at omega = eta_H k^2 with a sense fixed by sign(eta_H);
// total energy is conserved and the discrete J.E_H work vanishes.
double whistlerRotation(double eta_hall, double* hall_work, double* energy_drift,
                        std::uint64_t* rate_evaluations) {
  const std::int64_t cells = 128;
  const cosmosim::core::UnitSystem units = cgsUnits();
  const double delta = 1.0e-5;
  TestLevel level({cells, 4, 4}, {cells, 4, 4}, {1.0 / cells, 0.25, 0.25},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize([](std::array<double, 3>) { return GasPrimitive{}; }, {1.0e-4, 0.0, 0.0},
                   circularAlfvenPotential(delta, 1.0), 1.0, k_gamma);
  ConstantMhdCoefficientProvider provider(0.0, eta_hall, 0.0);
  const MhdNonIdealRateOperator rate(&provider, properties(units), MhdStiffnessLimits{}, false,
                                     true, false);
  MhdLevelAdvanceOptions options{};
  options.nonideal = &rate;
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  const MhdConservationBudget begin = stepper.levelBudget(1.0);
  const Mode initial = fieldMode(level, cells);
  const double end_time = 0.25;
  double time = 0.0;
  double work = 0.0;
  std::uint64_t evaluations = 0U;
  while (time < end_time * (1.0 - 1.0e-14)) {
    MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(time, 1.0);
    interval.dt_code = std::min(1.5e-3, end_time - time);
    const MhdLevelAdvanceResult result = stepper.advance(interval, options);
    CHUI_REQUIRE(result.accepted);
    work += result.transport.hall_work_integral;
    evaluations += 2U * result.transport.hall_subcycle_count;  // two SSPRK2 stages each
    time += result.dt_accepted_code;
  }
  *rate_evaluations = evaluations;
  const Mode final_mode = fieldMode(level, cells);
  const MhdConservationBudget end = stepper.levelBudget(1.0);
  *hall_work = work;
  *energy_drift = std::abs(end.energy_total - begin.energy_total) / begin.energy_total;
  // Rotation of the transverse field vector at the first Fourier mode.
  const std::complex<double> circular_initial = initial.by + std::complex<double>(0.0, 1.0) * initial.bz;
  const std::complex<double> circular_final = final_mode.by + std::complex<double>(0.0, 1.0) * final_mode.bz;
  return std::arg(circular_final / circular_initial) / end_time;
}

void testHallWhistler() {
  const double eta_hall = 0.05;
  double work = 0.0;
  double drift = 0.0;
  std::uint64_t evaluations = 0U;
  const double omega = whistlerRotation(eta_hall, &work, &drift, &evaluations);
  const double expected = eta_hall * k_two_pi * k_two_pi;
  CHUI_REQUIRE(std::abs(std::abs(omega) - expected) <= 0.02 * expected);  // N03 phase speed
  CHUI_REQUIRE(drift <= 1e-12);                                            // N04 energy
  // J . E_H vanishes pointwise; the discrete diagnostic sums E_c J_c over the
  // three edge families (different locations), so it vanishes to O((k dx)^2).
  // Bound it against the Joule rate an equal Ohmic diffusivity would produce,
  // eta k^2 |delta B|^2 V per evaluation (V = 1, |delta B| = 1e-5).
  const double joule_scale = eta_hall * k_two_pi * k_two_pi * 1.0e-10 *
      static_cast<double>(std::max<std::uint64_t>(1U, evaluations));
  CHUI_REQUIRE(std::isfinite(work) && std::abs(work) <= 0.05 * joule_scale);
  double work_reversed = 0.0;
  double drift_reversed = 0.0;
  std::uint64_t evaluations_reversed = 0U;
  const double omega_reversed =
      whistlerRotation(-eta_hall, &work_reversed, &drift_reversed, &evaluations_reversed);
  CHUI_REQUIRE(omega * omega_reversed < 0.0);  // polarity reversal with sign(eta_H)
  CHUI_REQUIRE(std::abs(std::abs(omega_reversed) - std::abs(omega)) <= 1e-6 * expected);
}

// N06/N07 smoke: all three terms with a stiff Ohmic coefficient (RKL2 path),
// closed-box energy and strict CT.
void testCombinedStrang() {
  const cosmosim::core::UnitSystem units = cgsUnits();
  TestLevel level({32, 32, 4}, {16, 16, 4}, {1.0 / 32, 1.0 / 32, 0.25},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize(
      [](std::array<double, 3> x) {
        GasPrimitive w{};
        w.velocity = {0.1 * std::sin(k_two_pi * x[1]), 0.1 * std::sin(k_two_pi * x[0]), 0.0};
        return w;
      },
      {0.05, 0.0, 0.0}, fieldLoopPotential(1.0e-2, 0.3, {0.5, 0.5, 0.5}), 1.0, k_gamma);
  // eta_O + eta_A = 0.51: dt_D = 0.15/(0.51 * 2064) = 1.4e-4, so the half
  // interval dt/2 = 1e-3 has R ~ 7 > rkl_switch_ratio (RKL2 branch); the ideal
  // CFL at dt = 2e-3 is ~0.2.
  ConstantMhdCoefficientProvider provider(0.5, 0.003, 0.01);
  const MhdNonIdealRateOperator rate(&provider, properties(units), MhdStiffnessLimits{}, true, true,
                                     true);
  MhdLevelAdvanceOptions options{};
  options.nonideal = &rate;
  options.produce_edge_integrals = true;
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  const MhdConservationBudget begin = stepper.levelBudget(1.0);
  const std::vector<double> phi0(level.faces(0).fluxX().begin(), level.faces(0).fluxX().end());
  MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(0.0, 2.0e-3);
  const MhdLevelAdvanceResult result = stepper.advance(interval, options);
  CHUI_REQUIRE(result.accepted);
  CHUI_REQUIRE(result.transport.rkl_cycle_count > 0U);
  CHUI_REQUIRE(result.transport.hall_subcycle_count > 0U);
  CHUI_REQUIRE(result.divergence_max_normalized <= 1e-11);
  const MhdConservationBudget end = stepper.levelBudget(1.0);
  CHUI_REQUIRE(std::abs(end.energy_total - begin.energy_total) <= 1e-12 * begin.energy_total);
  // Time-integrated EMF reproduces the accepted face change: Phi1 = Phi0 - C(I).
  const MhdPatchGeometry& geometry = level.binding(0).geometry;
  MhdFaceFluxStorage check(geometry);
  std::copy(phi0.begin(), phi0.end(), check.fluxX().begin());
  const MhdEdgeEmfStorage& integral = stepper.edgeIntegral(0);
  applyConstrainedTransportUpdate(geometry, check.fluxX(), check.fluxY(), check.fluxZ(),
                                  integral.emfX(), integral.emfY(), integral.emfZ(), 1.0);
  for (std::size_t index = 0; index < phi0.size(); ++index) {
    CHUI_REQUIRE(std::abs(check.fluxX()[index] - level.faces(0).fluxX()[index]) <=
                 1e-12 * (std::abs(phi0[index]) + 1e-12));
  }
}

// Stiffness: a Hall cap of 1 subcycle cannot cover the interval; the stepper
// halves (no mutation before the forecast) and either succeeds on a shorter
// interval or fails with MHD_STIFFNESS_BUDGET_EXCEEDED leaving state intact.
void testStiffnessCap() {
  const cosmosim::core::UnitSystem units = cgsUnits();
  TestLevel level({16, 4, 4}, {16, 4, 4}, {1.0 / 16, 0.25, 0.25}, cosmosim::mhd::test::periodicDomain());
  level.initialize([](std::array<double, 3>) { return GasPrimitive{}; }, {0.1, 0.0, 0.0},
                   circularAlfvenPotential(1e-3, 1.0), 1.0, k_gamma);
  ConstantMhdCoefficientProvider provider(0.0, 10.0, 0.0);
  MhdStiffnessLimits limits{};
  limits.hall_subcycle_cost_cap = 1;
  const MhdNonIdealRateOperator rate(&provider, properties(units), limits, false, true, false);
  MhdLevelAdvanceOptions options{};
  options.nonideal = &rate;
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  const MhdConservationBudget begin = stepper.levelBudget(1.0);
  MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(0.0, 5.0e-3);
  try {
    const MhdLevelAdvanceResult result = stepper.advance(interval, options);
    CHUI_REQUIRE(!result.accepted || result.retry_count > 0);
  } catch (const MhdError& error) {
    CHUI_REQUIRE(error.code() == MhdErrorCode::kStiffnessBudgetExceeded);
    const MhdConservationBudget after = stepper.levelBudget(1.0);
    CHUI_REQUIRE(after.energy_total == begin.energy_total && after.mass == begin.mass);
  }
}

// B = 0 with ambipolar and Hall: well-defined zero-field limit (no NaN, no
// invented Hall direction), field stays zero.
void testZeroFieldNonIdeal() {
  const cosmosim::core::UnitSystem units = cgsUnits();
  TestLevel level({8, 8, 4}, {8, 8, 4}, {0.125, 0.125, 0.25}, cosmosim::mhd::test::periodicDomain());
  level.initialize([](std::array<double, 3>) { return GasPrimitive{}; }, {0.0, 0.0, 0.0}, {}, 1.0,
                   k_gamma);
  ConstantMhdCoefficientProvider provider(0.0, 0.01, 0.01);
  const MhdNonIdealRateOperator rate(&provider, properties(units), MhdStiffnessLimits{}, false, true,
                                     true);
  MhdLevelAdvanceOptions options{};
  options.nonideal = &rate;
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  CHUI_REQUIRE(stepper.advance(cosmosim::mhd::test::isolatedInterval(0.0, 1e-3), options).accepted);
  for (const double value : level.faces(0).fluxY()) CHUI_REQUIRE(value == 0.0);
}

void testProviderAdmission() {
  ConstantMhdCoefficientProvider negative(-1.0, 0.0, 0.0);
  bool rejected = false;
  try {
    negative.validateAdmission();
  } catch (const MhdError& error) {
    rejected = error.code() == MhdErrorCode::kBadCoefficient;
  }
  CHUI_REQUIRE(rejected);
  ChemistryMhdCoefficientProvider chemistry;
  rejected = false;
  try {
    chemistry.validateAdmission();
  } catch (const MhdError& error) {
    rejected = error.code() == MhdErrorCode::kUnsupportedPhysics;
  }
  CHUI_REQUIRE(rejected);
}

#if COSMOSIM_ENABLE_HDF5
// N05 (table): a SYNTHETIC VERIFICATION FIXTURE written by the test (power-law
// Ohmic/ambipolar, sign-changing Hall, explicit BZeroLimit). It is not a
// physical transport table and carries no microphysical claim.
void writeAttribute(hid_t object, const char* name, const std::string& value) {
  const hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, value.size());
  const hid_t space = H5Screate(H5S_SCALAR);
  const hid_t attribute = H5Acreate2(object, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, type, value.data());
  H5Aclose(attribute);
  H5Sclose(space);
  H5Tclose(type);
}

void writeDoubles(hid_t file, const char* name, const std::vector<hsize_t>& dims,
                  const std::vector<double>& values, const char* units) {
  const hid_t space = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
  const hid_t dataset = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
  writeAttribute(dataset, "Units", units);
  H5Dclose(dataset);
  H5Sclose(space);
}

void writeStrings(hid_t file, const char* name, const std::vector<std::string>& values) {
  const hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, 16);
  std::vector<char> buffer(16U * values.size(), '\0');
  for (std::size_t index = 0; index < values.size(); ++index) {
    values[index].copy(buffer.data() + 16U * index, 16U);
  }
  const hsize_t count = values.size();
  const hid_t space = H5Screate_simple(1, &count, nullptr);
  const hid_t dataset = H5Dcreate2(file, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
  H5Dclose(dataset);
  H5Sclose(space);
  H5Tclose(type);
}

void testTableProvider() {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "chui_mhd_transport_fixture_v1.h5";
  const std::vector<double> log_rho{-24.0, -22.0, -20.0};
  const std::vector<double> log_t{1.0, 3.0};
  const std::vector<double> log_b{-6.0, -3.0};
  std::vector<double> ohmic;
  std::vector<double> hall;
  std::vector<double> ambipolar;
  for (const double r : log_rho) {
    for (const double t : log_t) {
      for (const double b : log_b) {
        ohmic.push_back(std::pow(10.0, 18.0 + 0.5 * r + 0.0 * t + 0.0 * b + 11.0));
        hall.push_back(1.0e17 * (r + 22.0));  // sign change at log rho = -22
        ambipolar.push_back(std::pow(10.0, 20.0 + 2.0 * b));
      }
    }
  }
  std::vector<double> zero_ohmic(log_rho.size() * log_t.size(), 7.0e16);
  std::vector<double> zero_hall(log_rho.size() * log_t.size(), 0.0);
  std::vector<double> zero_ambi(log_rho.size() * log_t.size(), 0.0);
  {
    const hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    const hid_t root = H5Gopen2(file, "/", H5P_DEFAULT);
    writeAttribute(root, "ContractId", "chui_mhd_transport_table_v1");
    H5Gclose(root);
    H5Gclose(H5Gcreate2(file, "/Axes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    H5Gclose(H5Gcreate2(file, "/BZeroLimit", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    writeDoubles(file, "/Axes/LogDensityCgs", {3}, log_rho, "log10(g/cm^3)");
    writeDoubles(file, "/Axes/LogTemperatureK", {2}, log_t, "log10(K)");
    writeDoubles(file, "/Axes/LogMagneticFieldGauss", {2}, log_b, "log10(G)");
    writeStrings(file, "/Axes/CompositionIds", {"solar"});
    writeStrings(file, "/Axes/GrainModelIds", {"none"});
    writeDoubles(file, "/EtaOhmicCgs", {3, 2, 2, 1, 1}, ohmic, "cm^2/s");
    writeDoubles(file, "/EtaHallCgs", {3, 2, 2, 1, 1}, hall, "cm^2/s");
    writeDoubles(file, "/EtaAmbipolarCgs", {3, 2, 2, 1, 1}, ambipolar, "cm^2/s");
    writeDoubles(file, "/BZeroLimit/EtaOhmicCgs", {3, 2, 1, 1}, zero_ohmic, "cm^2/s");
    writeDoubles(file, "/BZeroLimit/EtaHallCgs", {3, 2, 1, 1}, zero_hall, "cm^2/s");
    writeDoubles(file, "/BZeroLimit/EtaAmbipolarCgs", {3, 2, 1, 1}, zero_ambi, "cm^2/s");
    H5Fclose(file);
  }
  TableMhdCoefficientProvider provider(path, {.composition_id = "solar", .grain_model_id = "none",
                                              .expected_sha256 = {}});
  provider.validateAdmission();
  CHUI_REQUIRE(provider.tableSha256().has_value() && provider.tableSha256()->size() == 64U);
  MhdCellPhysicalProperties inside{.rho_phys_cgs = 1.0e-21, .temperature_k = 100.0,
                                   .b_phys_gauss = 1.0e-4, .closure_valid = true};
  const MhdTransportCoefficients value = provider.evaluate(0, inside);
  CHUI_REQUIRE(value.valid);
  // Log interpolation is exact for the power-law Ohmic and ambipolar entries.
  CHUI_REQUIRE(std::abs(value.eta_ohmic_phys_cgs - std::pow(10.0, 29.0 - 10.5)) <=
               1e-10 * value.eta_ohmic_phys_cgs);
  CHUI_REQUIRE(std::abs(value.eta_ambipolar_phys_cgs - std::pow(10.0, 12.0)) <=
               1e-10 * value.eta_ambipolar_phys_cgs);
  CHUI_REQUIRE(value.eta_hall_phys_cgs > 0.0);  // log rho = -21 > -22
  inside.rho_phys_cgs = 1.0e-22;
  CHUI_REQUIRE(std::abs(provider.evaluate(0, inside).eta_hall_phys_cgs) <= 1e-3);  // exact zero crossing
  inside.b_phys_gauss = 0.0;  // exact zero field -> BZeroLimit
  CHUI_REQUIRE(std::abs(provider.evaluate(0, inside).eta_ohmic_phys_cgs - 7.0e16) <= 1e3);
  inside.b_phys_gauss = 1.0;  // log B = 0 > -3: outside coverage, no extrapolation
  CHUI_REQUIRE(!provider.evaluate(0, inside).valid);
  // A wrong SHA-256 pin is rejected at admission.
  TableMhdCoefficientProvider pinned(path, {.composition_id = "solar", .grain_model_id = "none",
                                            .expected_sha256 = std::string(64, '0')});
  bool rejected = false;
  try {
    pinned.validateAdmission();
  } catch (const MhdError& error) {
    rejected = error.code() == MhdErrorCode::kBadCoefficient;
  }
  CHUI_REQUIRE(rejected);
  std::filesystem::remove(path);
}
#endif

}  // namespace

int main() {
  testRkl2Recurrence();
  testCoefficientUnitConversion();
  testProviderAdmission();
  testZeroFieldNonIdeal();
  testStiffnessCap();
  testCombinedStrang();
  testOhmicDecay();
  testHallWhistler();
#if COSMOSIM_ENABLE_HDF5
  testTableProvider();
#endif
  return 0;
}
