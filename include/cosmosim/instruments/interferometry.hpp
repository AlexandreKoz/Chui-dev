#pragma once

// CHUI-EXT-PHYS-REF-v1 O01-O05: standalone virtual interferometer. Consumes an
// immutable sky (point components or a pixelized cube); never touches the
// simulation state.
//
// Geometry (Thompson, Moran & Swenson, "Interferometry and Synthesis in Radio
// Astronomy", eqs. 4.1 and 4.15; recalled lead): station ENU (E, N, U) at
// latitude phi -> equatorial (X, Y, Z) = (-sin phi N + cos phi U, E,
// cos phi N + sin phi U); baseline L = X_q - X_p; for hour angle H and
// declination delta of the phase centre
//   u = (sin H Lx + cos H Ly) / lambda,
//   v = (-sin d cos H Lx + sin d sin H Ly + cos d Lz) / lambda,
//   w = (cos d cos H Lx - cos d sin H Ly + sin d Lz) / lambda,
// H(t) = H0 + Omega_sid (t - t0), Omega_sid = 2 pi / 86164.0905 s (declared).
// Measurement equation (direct reference, full polarization, linear feeds):
//   V_pq = sum_s J_p(s) B_s J_q(s)^H exp(-2 pi i (u l + v m + w (n - 1))),
//   n = sqrt(1 - l^2 - m^2), B = [[I + Q, U + i V], [U - i V, I - Q]]
//   (Stokes normalization: XX + YY = 2 I for unit Jones; point components
//   carry flux densities, so no 1/n factor; declared).
// Jones chain per station J = G(nu) D Ionosphere(nu) E(s, nu):
//   G = diag(g_x, g_y) exp(-2 pi i nu tau) (complex gain, bandpass delay),
//   D = [[1, d_x], [-d_y, 1]] (polarization leakage),
//   Ionosphere: scalar phase e^{i phi}, phi = e^2 TEC / (4 pi eps0 m_e c nu)
//     (cold plasma, first order; sign convention declared) times Faraday
//     rotation R(chi), chi = RM_ion lambda^2,
//   E: scalar analytic voltage beam (none | Gaussian exp(-theta^2 / (2 s^2)),
//     s = FWHM_power / (2 sqrt(2 ln 2)) * sqrt(2) | Airy 2 J1(x)/x,
//     x = pi D sin(theta) / lambda). Measured/EM beams are not bundled: a
//     missing beam is an explicit `none` choice, never a silent default.
// Averaging: every visibility is the mean over `time_subsamples` x the
// channel's Gauss-Legendre frequency nodes (explicit time and bandwidth
// smearing). Noise: complex Gaussian with sigma_re = sigma_im =
// SEFD / sqrt(2 dnu tau) per correlation (radiometer equation, unit
// efficiency, declared), counter-based deterministic per (seed, baseline,
// time, channel, correlation). Flags: flagged samples are not computed and
// are written with flag = 1 and zero data. Output is streamed per time step
// through a sink (no resident full visibility cube): N_bl = N_a (N_a - 1)/2.
// Accelerated backends (GPU direct sum, NUFFT, antenna matrix): FAIL_CLOSED.

#include <array>
#include <complex>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "cosmosim/spectral/frequency_axis.hpp"

namespace cosmosim::instruments {

using Complex = std::complex<double>;
using Jones = std::array<Complex, 4>;  // row-major 2 x 2

struct Station {
  std::string name;
  std::array<double, 3> enu_m{};
  Complex gain_x{1.0, 0.0};
  Complex gain_y{1.0, 0.0};
  double delay_s = 0.0;
  Complex leakage_x{0.0, 0.0};
  Complex leakage_y{0.0, 0.0};
  double ionosphere_tec_m2 = 0.0;      // electrons per m^2 (1 TECU = 1e16)
  double ionosphere_rm_rad_m2 = 0.0;
};

enum class BeamModel : std::uint8_t { kNone = 0, kGaussian = 1, kAiry = 2 };
struct BeamSettings {
  BeamModel model = BeamModel::kNone;
  double fwhm_rad = 0.0;      // power-beam FWHM at reference frequency (Gaussian)
  double reference_hz = 0.0;  // Gaussian width scales as lambda
  double dish_diameter_m = 0.0;  // Airy
};

struct SkyComponent {
  double l = 0.0;  // direction cosines relative to the phase centre
  double m = 0.0;
  std::vector<std::array<double, 4>> stokes_per_channel;  // Jy, (I, Q, U, V)
};

struct ObservationSettings {
  double latitude_rad = 0.0;
  double phase_centre_dec_rad = 0.0;
  double hour_angle_start_rad = 0.0;
  double integration_s = 1.0;
  std::size_t time_steps = 1;
  std::uint32_t time_subsamples = 1;
  spectral::FrequencyAxis frequencies;
  double sefd_jy = 0.0;        // 0 = noiseless
  std::uint64_t noise_seed = 1;
  BeamSettings beam;
};

struct Baseline {
  std::uint32_t p = 0;
  std::uint32_t q = 0;
};

struct VisibilityBatch {
  std::size_t time_index = 0;
  double hour_angle_rad = 0.0;
  std::vector<Baseline> baselines;
  std::vector<std::array<double, 3>> uvw_m;    // per baseline (metres, at the batch mid-time)
  std::vector<std::array<Complex, 4>> data;    // [baseline * channels + channel] (XX, XY, YX, YY)
  std::vector<std::uint8_t> flags;             // same indexing
};

using VisibilitySink = std::function<void(const VisibilityBatch&)>;
// flag(baseline_index, time_index, channel) -> true when flagged
using FlagFunction = std::function<bool(std::size_t, std::size_t, std::size_t)>;

[[nodiscard]] std::array<double, 3> enuToEquatorial(const std::array<double, 3>& enu, double latitude_rad);
[[nodiscard]] std::array<double, 3> baselineUvwMetres(const std::array<double, 3>& baseline_xyz, double hour_angle,
                                                      double declination);
[[nodiscard]] double beamVoltage(const BeamSettings& beam, double l, double m, double frequency_hz);
[[nodiscard]] Jones stationJones(const Station& station, const BeamSettings& beam, double l, double m, double nu);
[[nodiscard]] Jones brightnessMatrix(const std::array<double, 4>& stokes);
[[nodiscard]] double thermalNoiseSigmaJy(double sefd_jy, double bandwidth_hz, double integration_s);
[[nodiscard]] std::uint64_t baselineCount(std::uint64_t stations) noexcept;

struct RimeReport {
  std::uint64_t visibilities = 0;
  std::uint64_t flagged = 0;
  std::uint64_t component_evaluations = 0;
};

// Direct RIME (CPU reference). Streams one batch per time step to `sink`.
RimeReport simulateVisibilitiesDirect(const std::vector<Station>& stations, const std::vector<SkyComponent>& sky,
                                      const ObservationSettings& settings, const VisibilitySink& sink,
                                      const FlagFunction& flags = {});

// Native self-describing product writer (schema chui_visibilities_v1).
class NativeVisibilityWriter {
 public:
  NativeVisibilityWriter(std::string path, const std::vector<Station>& stations, const ObservationSettings& settings,
                         std::string provenance);
  void consume(const VisibilityBatch& batch);
  void finish();

 private:
  std::string m_path;
  std::string m_partial;
  bool m_open = false;
  std::uint64_t m_batches = 0;
};

// OSKAR interoperability export (inputs only; OSKAR is not executed):
// sky model text 'RA_deg Dec_deg I Q U V ref_freq_hz spectral_index' for the
// first channel and a telescope directory with layout.txt (E N U metres) and
// position.txt (longitude latitude altitude). Format: recalled lead of the
// OSKAR documentation, to verify in the interoperability test campaign.
void exportOskarInputs(const std::string& directory, const std::vector<Station>& stations,
                       const std::vector<SkyComponent>& sky, double phase_centre_ra_rad, double phase_centre_dec_rad,
                       double longitude_rad, double latitude_rad, double reference_frequency_hz);

}  // namespace cosmosim::instruments
