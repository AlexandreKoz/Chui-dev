// CHUI-EXT-PHYS-REF-v1 O01-O05 virtual interferometer (source only, NOT
// compiled or executed in pass 01).

#include "cosmosim/instruments/interferometry.hpp"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::instruments {
namespace {

namespace k = core::constants;
constexpr double k_sidereal_rate = 2.0 * k::k_pi / 86164.0905;  // rad s^-1 (sidereal day, declared)

[[nodiscard]] Jones multiply(const Jones& a, const Jones& b) {
  return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3], a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]};
}

[[nodiscard]] Jones hermitian(const Jones& a) { return {std::conj(a[0]), std::conj(a[2]), std::conj(a[1]), std::conj(a[3])}; }

// Direction-independent station terms A = G D Ion (2 x 2) at frequency nu.
[[nodiscard]] Jones stationTerms(const Station& s, double nu) {
  const Complex delay = std::polar(1.0, -2.0 * k::k_pi * nu * s.delay_s);
  const Jones g{s.gain_x * delay, 0.0, 0.0, s.gain_y * delay};
  const Jones d{1.0, s.leakage_x, -s.leakage_y, 1.0};
  const double e = k::k_elementary_charge_si;
  const double phase = e * e * s.ionosphere_tec_m2 /
                       (4.0 * k::k_pi * k::k_vacuum_permittivity_si * k::k_electron_mass_si * k::k_speed_of_light_si * nu);
  const double lambda = k::k_speed_of_light_si / nu;
  const double chi = s.ionosphere_rm_rad_m2 * lambda * lambda;
  const Complex ph = std::polar(1.0, phase);
  const Jones ion{ph * std::cos(chi), -ph * std::sin(chi), ph * std::sin(chi), ph * std::cos(chi)};
  return multiply(multiply(g, d), ion);
}

[[nodiscard]] double gaussian(std::uint64_t seed, std::uint64_t a, std::uint64_t b) {
  // SplitMix64 counter streams -> Box-Muller.
  const auto mix = [](std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31U);
  };
  const std::uint64_t h = mix(mix(seed) ^ mix(a) ^ (mix(b) << 1U));
  const double u1 = (static_cast<double>(mix(h) >> 11U) + 0.5) / 9007199254740992.0;
  const double u2 = (static_cast<double>(mix(h + 1U) >> 11U) + 0.5) / 9007199254740992.0;
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * k::k_pi * u2);
}

// J1 from the periodic Bessel integral J1(x) = (1/2pi) int_{-pi}^{pi} cos(t - x sin t) dt,
// trapezoid on 256 nodes (exponentially convergent for a periodic integrand,
// accurate to ~1e-14 for x < 100); large x uses the leading asymptotic form.
// Avoids std::cyl_bessel_j, which libc++ does not provide.
[[nodiscard]] double besselJ1(double x) {
  if (std::abs(x) > 100.0) {
    return std::sqrt(2.0 / (k::k_pi * std::abs(x))) * std::cos(std::abs(x) - 0.75 * k::k_pi) * (x < 0.0 ? -1.0 : 1.0);
  }
  constexpr int n = 256;
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    const double t = -k::k_pi + 2.0 * k::k_pi * i / n;
    sum += std::cos(t - x * std::sin(t));
  }
  return sum / n;
}

void putBytes(std::ofstream& out, const void* data, std::size_t n) {
  out.write(static_cast<const char*>(data), static_cast<std::streamsize>(n));
}

}  // namespace

std::array<double, 3> enuToEquatorial(const std::array<double, 3>& enu, double lat) {
  const double e = enu[0];
  const double n = enu[1];
  const double u = enu[2];
  return {-std::sin(lat) * n + std::cos(lat) * u, e, std::cos(lat) * n + std::sin(lat) * u};
}

std::array<double, 3> baselineUvwMetres(const std::array<double, 3>& l, double h, double d) {
  return {std::sin(h) * l[0] + std::cos(h) * l[1],
          -std::sin(d) * std::cos(h) * l[0] + std::sin(d) * std::sin(h) * l[1] + std::cos(d) * l[2],
          std::cos(d) * std::cos(h) * l[0] - std::cos(d) * std::sin(h) * l[1] + std::sin(d) * l[2]};
}

double beamVoltage(const BeamSettings& beam, double l, double m, double nu) {
  const double r2 = l * l + m * m;
  if (r2 > 1.0) return 0.0;  // below the horizon of the tangent plane
  const double sin_theta = std::sqrt(r2);
  switch (beam.model) {
    case BeamModel::kNone:
      return 1.0;
    case BeamModel::kGaussian: {
      if (!(beam.fwhm_rad > 0.0) || !(beam.reference_hz > 0.0)) throw std::invalid_argument("Gaussian beam needs FWHM and reference frequency");
      const double fwhm = beam.fwhm_rad * beam.reference_hz / nu;
      const double sigma_v = fwhm / (2.0 * std::sqrt(std::log(2.0)));
      const double theta = std::asin(sin_theta);
      return std::exp(-theta * theta / (2.0 * sigma_v * sigma_v));
    }
    case BeamModel::kAiry: {
      if (!(beam.dish_diameter_m > 0.0)) throw std::invalid_argument("Airy beam needs a dish diameter");
      const double x = k::k_pi * beam.dish_diameter_m * sin_theta * nu / k::k_speed_of_light_si;
      return x < 1.0e-8 ? 1.0 : 2.0 * besselJ1(x) / x;
    }
  }
  return 1.0;
}

Jones stationJones(const Station& station, const BeamSettings& beam, double l, double m, double nu) {
  const Jones a = stationTerms(station, nu);
  const double e = beamVoltage(beam, l, m, nu);
  return {a[0] * e, a[1] * e, a[2] * e, a[3] * e};
}

Jones brightnessMatrix(const std::array<double, 4>& s) {
  return {Complex(s[0] + s[1], 0.0), Complex(s[2], s[3]), Complex(s[2], -s[3]), Complex(s[0] - s[1], 0.0)};
}

double thermalNoiseSigmaJy(double sefd, double bandwidth, double tau) {
  if (!(bandwidth > 0.0) || !(tau > 0.0)) throw std::invalid_argument("noise needs bandwidth and integration > 0");
  return sefd / std::sqrt(2.0 * bandwidth * tau);
}

std::uint64_t baselineCount(std::uint64_t stations) noexcept { return stations * (stations - 1U) / 2U; }

RimeReport simulateVisibilitiesDirect(const std::vector<Station>& stations, const std::vector<SkyComponent>& sky,
                                      const ObservationSettings& settings, const VisibilitySink& sink,
                                      const FlagFunction& flags) {
  if (stations.size() < 2U) throw std::invalid_argument("interferometer needs >= 2 stations");
  const std::size_t channels = settings.frequencies.channelCount();
  if (channels == 0U || settings.time_steps == 0U || settings.time_subsamples == 0U || !(settings.integration_s > 0.0)) {
    throw std::invalid_argument("observation needs channels, time steps, subsamples and integration > 0");
  }
  for (const SkyComponent& c : sky) {
    if (c.stokes_per_channel.size() != channels) throw std::invalid_argument("sky component channel count mismatch");
    if (c.l * c.l + c.m * c.m >= 1.0) throw std::invalid_argument("sky component outside the tangent plane");
  }
  RimeReport report;
  std::vector<std::array<double, 3>> xyz(stations.size());
  for (std::size_t s = 0; s < stations.size(); ++s) xyz[s] = enuToEquatorial(stations[s].enu_m, settings.latitude_rad);
  std::vector<Baseline> baselines;
  for (std::uint32_t p = 0; p < stations.size(); ++p) {
    for (std::uint32_t q = p + 1U; q < stations.size(); ++q) baselines.push_back({p, q});
  }
  const std::vector<double>& weights = settings.frequencies.nodeWeights();
  for (std::size_t t = 0; t < settings.time_steps; ++t) {
    VisibilityBatch batch;
    batch.time_index = t;
    const double t_mid = (static_cast<double>(t) + 0.5) * settings.integration_s;
    batch.hour_angle_rad = settings.hour_angle_start_rad + k_sidereal_rate * t_mid;
    batch.baselines = baselines;
    batch.uvw_m.resize(baselines.size());
    batch.data.assign(baselines.size() * channels, {Complex{}, Complex{}, Complex{}, Complex{}});
    batch.flags.assign(baselines.size() * channels, 0U);
    for (std::size_t b = 0; b < baselines.size(); ++b) {
      std::array<double, 3> l{};
      for (std::size_t a = 0; a < 3U; ++a) l[a] = xyz[baselines[b].q][a] - xyz[baselines[b].p][a];
      batch.uvw_m[b] = baselineUvwMetres(l, batch.hour_angle_rad, settings.phase_centre_dec_rad);
      for (std::size_t c = 0; c < channels; ++c) {
        const std::size_t idx = b * channels + c;
        if (flags && flags(b, t, c)) {
          batch.flags[idx] = 1U;
          ++report.flagged;
          continue;
        }
        const std::vector<double> nodes = settings.frequencies.nodeFrequenciesHz(c);
        std::array<Complex, 4> sum{};
        for (std::uint32_t ts = 0; ts < settings.time_subsamples; ++ts) {
          const double t_sub = (static_cast<double>(t) + (ts + 0.5) / settings.time_subsamples) * settings.integration_s;
          const double h = settings.hour_angle_start_rad + k_sidereal_rate * t_sub;
          const auto uvw = baselineUvwMetres(l, h, settings.phase_centre_dec_rad);
          for (std::size_t f = 0; f < nodes.size(); ++f) {
            const double nu = nodes[f];
            const double inv_lambda = nu / k::k_speed_of_light_si;
            Jones coherency{};
            for (const SkyComponent& comp : sky) {
              const double n = std::sqrt(1.0 - comp.l * comp.l - comp.m * comp.m);
              const double phase = -2.0 * k::k_pi * inv_lambda * (uvw[0] * comp.l + uvw[1] * comp.m + uvw[2] * (n - 1.0));
              const double e = beamVoltage(settings.beam, comp.l, comp.m, nu);
              const Jones bm = brightnessMatrix(comp.stokes_per_channel[c]);
              const Complex factor = e * e * std::polar(1.0, phase);
              for (std::size_t q = 0; q < 4U; ++q) coherency[q] += factor * bm[q];
              ++report.component_evaluations;
            }
            const Jones v = multiply(multiply(stationTerms(stations[baselines[b].p], nu), coherency),
                                     hermitian(stationTerms(stations[baselines[b].q], nu)));
            const double w = weights[f] / settings.time_subsamples;
            for (std::size_t q = 0; q < 4U; ++q) sum[q] += w * v[q];
          }
        }
        if (settings.sefd_jy > 0.0) {
          const double sigma = thermalNoiseSigmaJy(settings.sefd_jy, settings.frequencies.widthHz(c), settings.integration_s);
          for (std::size_t q = 0; q < 4U; ++q) {
            const std::uint64_t key = ((static_cast<std::uint64_t>(b) * settings.time_steps + t) * channels + c) * 4U + q;
            sum[q] += Complex(sigma * gaussian(settings.noise_seed, key, 0U), sigma * gaussian(settings.noise_seed, key, 1U));
          }
        }
        batch.data[idx] = sum;
        ++report.visibilities;
      }
    }
    sink(batch);
  }
  return report;
}

NativeVisibilityWriter::NativeVisibilityWriter(std::string path, const std::vector<Station>& stations,
                                               const ObservationSettings& settings, std::string provenance)
    : m_path(std::move(path)), m_partial(m_path + ".partial") {
  std::ofstream out(m_partial, std::ios::binary | std::ios::trunc);
  if (!out) throw std::runtime_error("cannot open visibility output '" + m_partial + "'");
  std::ostringstream header;
  header.precision(17);
  header << "schema=chui_visibilities_v1\nbyte_order=native\ncorrelations=XX,XY,YX,YY\nunits=Jy\n"
         << "stokes_normalization=XX+YY=2I\nfourier_sign=exp(-2 pi i (u l + v m + w (n - 1)))\n"
         << "frequencies=" << settings.frequencies.describe() << "\nintegration_s=" << settings.integration_s
         << "\ntime_steps=" << settings.time_steps << "\ntime_subsamples=" << settings.time_subsamples
         << "\nlatitude_rad=" << settings.latitude_rad << "\nphase_centre_dec_rad=" << settings.phase_centre_dec_rad
         << "\nhour_angle_start_rad=" << settings.hour_angle_start_rad << "\nsefd_jy=" << settings.sefd_jy
         << "\nnoise_seed=" << settings.noise_seed << "\nbeam_model=" << static_cast<int>(settings.beam.model)
         << "\nstations=" << stations.size() << '\n';
  for (const Station& s : stations) header << "station=" << s.name << ',' << s.enu_m[0] << ',' << s.enu_m[1] << ',' << s.enu_m[2] << '\n';
  header << "provenance=" << provenance << "\nrecord=u64 time_index, f64 hour_angle, u64 n_baselines, n_baselines x (u32 p, u32 q, 3 f64 uvw_m), "
         << "n_baselines x channels x (8 f64 re/im XX XY YX YY, u8 flag)\n";
  const std::string text = header.str();
  const char magic[8] = {'C', 'H', 'U', 'I', 'V', 'I', 'S', '1'};
  const std::uint64_t length = text.size();
  putBytes(out, magic, 8U);
  putBytes(out, &length, sizeof(length));
  putBytes(out, text.data(), text.size());
  if (!out) throw std::runtime_error("visibility header write failed");
  m_open = true;
}

void NativeVisibilityWriter::consume(const VisibilityBatch& batch) {
  if (!m_open) throw std::logic_error("visibility writer already finished");
  std::ofstream out(m_partial, std::ios::binary | std::ios::app);
  const std::uint64_t t = batch.time_index;
  const std::uint64_t nb = batch.baselines.size();
  putBytes(out, &t, sizeof(t));
  putBytes(out, &batch.hour_angle_rad, sizeof(double));
  putBytes(out, &nb, sizeof(nb));
  for (std::size_t b = 0; b < batch.baselines.size(); ++b) {
    putBytes(out, &batch.baselines[b].p, sizeof(std::uint32_t));
    putBytes(out, &batch.baselines[b].q, sizeof(std::uint32_t));
    putBytes(out, batch.uvw_m[b].data(), 3U * sizeof(double));
  }
  for (std::size_t i = 0; i < batch.data.size(); ++i) {
    for (const Complex& c : batch.data[i]) {
      const double re = c.real();
      const double im = c.imag();
      putBytes(out, &re, sizeof(double));
      putBytes(out, &im, sizeof(double));
    }
    putBytes(out, &batch.flags[i], 1U);
  }
  if (!out) throw std::runtime_error("visibility batch write failed");
  ++m_batches;
}

void NativeVisibilityWriter::finish() {
  if (!m_open) return;
  m_open = false;
  std::filesystem::rename(m_partial, m_path);
}

void exportOskarInputs(const std::string& directory, const std::vector<Station>& stations,
                       const std::vector<SkyComponent>& sky, double ra0, double dec0, double lon, double lat,
                       double ref_hz) {
  namespace fs = std::filesystem;
  const fs::path dir(directory);
  fs::create_directories(dir / "telescope");
  {
    std::ofstream layout(dir / "telescope" / "layout.txt", std::ios::trunc);
    layout.precision(17);
    for (const Station& s : stations) layout << s.enu_m[0] << ' ' << s.enu_m[1] << ' ' << s.enu_m[2] << '\n';
    std::ofstream position(dir / "telescope" / "position.txt", std::ios::trunc);
    position.precision(17);
    position << lon * 180.0 / k::k_pi << ' ' << lat * 180.0 / k::k_pi << " 0\n";
    if (!layout || !position) throw std::runtime_error("OSKAR telescope export failed");
  }
  std::ofstream model(dir / "sky_model.osm", std::ios::trunc);
  model.precision(17);
  model << "# CHUI-EXT-PHYS-REF-v1 export: RA_deg Dec_deg I Q U V ref_freq_hz spectral_index (first channel)\n";
  for (const SkyComponent& c : sky) {
    const double n = std::sqrt(1.0 - c.l * c.l - c.m * c.m);
    const double dec = std::asin(c.m * std::cos(dec0) + n * std::sin(dec0));
    const double ra = ra0 + std::atan2(c.l, n * std::cos(dec0) - c.m * std::sin(dec0));
    const auto& s = c.stokes_per_channel.front();
    model << ra * 180.0 / k::k_pi << ' ' << dec * 180.0 / k::k_pi << ' ' << s[0] << ' ' << s[1] << ' ' << s[2] << ' '
          << s[3] << ' ' << ref_hz << " 0\n";
  }
  if (!model) throw std::runtime_error("OSKAR sky model export failed");
}

}  // namespace cosmosim::instruments
