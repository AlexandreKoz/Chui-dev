// CHUI-EXT-PHYS-REF-v1 I05 inline observable pipeline (source only, NOT
// compiled or executed in pass 01).

#include "cosmosim/observables/observable_pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "core/internal/sha256.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/lightcone/ray_geometry.hpp"
#include "cosmosim/observables/cmb_secondaries.hpp"
#include "cosmosim/observables/fits_writer.hpp"
#include "cosmosim/observables/gas_science_view.hpp"
#include "cosmosim/observables/radio_transfer.hpp"

namespace cosmosim::observables {
namespace {

namespace k = core::constants;

[[nodiscard]] std::string fileName(const core::SimulationConfig& config, const std::string& product, std::uint64_t step) {
  std::ostringstream name;
  name << config.extended.observables.output_stem << '_' << product << "_step" << step;
  return name.str();
}

}  // namespace

ObservablePipeline::ObservablePipeline(const core::SimulationConfig& config)
    : m_config(config),
      m_units(core::makeUnitSystem(config.units.length_unit, config.units.mass_unit, config.units.velocity_unit)),
      m_background(core::buildModePolicy(config.mode).cosmological_comoving_frame
                       ? core::LambdaCdmBackground(core::cosmologyBackgroundConfigFromSimulation(config))
                       : core::LambdaCdmBackground()) {
  const core::ObservablesConfig& o = config.extended.observables;
  if (!o.enabled) throw std::invalid_argument("ObservablePipeline constructed with observables.enabled=false");
  if (o.hi_cube && o.hi_spin_model == core::HiSpinModelConfig::kCollisionalTable) {
    std::ifstream input(o.hi_collision_table_path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open 21 cm collision table '" + o.hi_collision_table_path + "'");
    const std::string bytes((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    if (core::internal::sha256Hex(bytes) != o.hi_collision_table_sha256) {
      throw std::runtime_error("21 cm collision table SHA-256 mismatch");
    }
    m_collision_table =
        std::make_unique<CollisionalCouplingTable>(parseCollisionalCouplingTable(bytes, o.hi_collision_table_path));
  }
}

bool ObservablePipeline::due(std::uint64_t step) const noexcept {
  const auto interval = m_config.extended.observables.interval_steps;
  return interval > 0U && step % interval == 0U;
}

ObservablePipelineReport ObservablePipeline::execute(const core::SimulationState& state, std::uint64_t step,
                                                     double a) {
  ObservablePipelineReport report;
  report.step = step;
  const core::ObservablesConfig& o = m_config.extended.observables;
  const bool cosmological = core::buildModePolicy(m_config.mode).cosmological_comoving_frame;
  // Domain = bounding box of the populated leaves (periodic box for cosmological runs).
  std::array<double, 3> lo{1e300, 1e300, 1e300};
  std::array<double, 3> hi{-1e300, -1e300, -1e300};
  for (std::size_t p = 0; p < state.patches.size(); ++p) {
    if (state.patches.cell_count[p] == 0U) continue;
    const std::array<double, 3> origin{state.patches.origin_x_comoving[p], state.patches.origin_y_comoving[p],
                                       state.patches.origin_z_comoving[p]};
    const std::array<double, 3> extent{state.patches.extent_x_comoving[p], state.patches.extent_y_comoving[p],
                                       state.patches.extent_z_comoving[p]};
    for (std::size_t ax = 0; ax < 3U; ++ax) {
      lo[ax] = std::min(lo[ax], origin[ax]);
      hi[ax] = std::max(hi[ax], origin[ax] + extent[ax]);
    }
  }
  if (!(hi[0] > lo[0])) throw std::runtime_error("observables: no populated leaf patch (no gas geometry)");
  const std::array<double, 3> extent{hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
  const auto geometry = lightcone::LeafGeometry::fromPatches(
      state.patches, lo, extent, cosmological ? lightcone::RayBoundary::kPeriodic : lightcone::RayBoundary::kOpen);
  const GasScienceAccessor gas(state, m_units, a);
  const double length_m = m_units.length_si_per_code;
  const std::size_t axis = static_cast<std::size_t>(o.los_axis);
  const std::size_t t1 = (axis + 1U) % 3U;
  const std::size_t t2 = (axis + 2U) % 3U;
  const std::size_t n = o.sightlines_per_side;
  std::array<double, 3> dir{0.0, 0.0, 0.0};
  dir[axis] = 1.0;
  const double hubble = cosmological ? m_background.hubbleSi(a) : 0.0;
  const double redshift = 1.0 / a - 1.0;
  const double los_length_code = extent[axis];
  const double t_cmb = cosmological && m_config.cosmology.radiation_species.enabled ? m_background.cmbTemperatureK(a) : 0.0;

  std::vector<double> forest_flux;
  std::vector<double> mean_flux_per_sightline;
  std::vector<double> mean_power;
  std::vector<double> mean_pdf;
  std::vector<double> power_k;
  std::vector<double> hi_cube;
  std::vector<double> y_map, tau_map, ksz_map, dm_map;
  VelocityAxis forest_axis;
  if (o.lyman_alpha_forest) {
    forest_axis.u_min = 0.0;
    forest_axis.u_max = a * hubble * los_length_code * length_m;
    forest_axis.pixels = o.forest_pixels;
    forest_axis.periodic = true;
    forest_flux.assign(static_cast<std::size_t>(n) * n * o.forest_pixels, 0.0);
  }
  VelocityAxis hi_axis;
  if (o.hi_cube) {
    hi_axis.u_min = -o.hi_velocity_half_width_km_s * 1.0e3;
    hi_axis.u_max = o.hi_velocity_half_width_km_s * 1.0e3;
    hi_axis.pixels = o.hi_channels;
    hi_axis.periodic = false;
    hi_cube.assign(static_cast<std::size_t>(o.hi_channels) * n * n, 0.0);
  }
  if (o.sz_maps) {
    y_map.assign(static_cast<std::size_t>(n) * n, 0.0);
    tau_map = y_map;
    ksz_map = y_map;
  }
  if (o.dispersion_map) dm_map.assign(static_cast<std::size_t>(n) * n, 0.0);
  const LineProfileModel forest_model =
      o.forest_profile == core::ForestProfileConfig::kVoigt ? LineProfileModel::kVoigt : LineProfileModel::kDopplerOnly;

  for (std::size_t j = 0; j < n; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      std::array<double, 3> start = lo;
      start[t1] = lo[t1] + (static_cast<double>(i) + 0.5) / static_cast<double>(n) * extent[t1];
      start[t2] = lo[t2] + (static_cast<double>(j) + 0.5) / static_cast<double>(n) * extent[t2];
      const auto segments = geometry.traceStraightRay(start, dir, los_length_code, std::size_t{1} << 22U);
      report.segments += segments.size();
      ++report.sightlines;
      std::vector<LineSegment> forest_segments;
      std::vector<LineSegment> hi_segments;
      std::vector<SzSegment> sz_segments;
      std::vector<DispersionSegment> dm_segments;
      for (const auto& seg : segments) {
        const std::size_t row = seg.cell_row;
        const double path = a * seg.length() * length_m;
        const auto v = gas.peculiarVelocity(row);
        const double v_los = v[0] * dir[0] + v[1] * dir[1] + v[2] * dir[2];
        if (o.lyman_alpha_forest || o.hi_cube) {
          LineSegment line;
          line.path_proper_m = path;
          line.absorber_density_m3 = gas.neutralHydrogenDensity(row).require("n_HI");
          line.temperature_k = gas.temperature(row).require("T");
          if (o.lyman_alpha_forest) {
            LineSegment f = line;
            f.velocity_begin_m_s = a * hubble * seg.s_begin * length_m + v_los;
            f.velocity_end_m_s = a * hubble * seg.s_end * length_m + v_los;
            forest_segments.push_back(f);
          }
          if (o.hi_cube) {
            LineSegment h = line;
            const double centre = 0.5 * los_length_code;
            h.velocity_begin_m_s = a * hubble * (seg.s_begin - centre) * length_m + v_los;
            h.velocity_end_m_s = a * hubble * (seg.s_end - centre) * length_m + v_los;
            if (o.hi_spin_model == core::HiSpinModelConfig::kSaturatedKinetic) {
              h.spin_temperature_k = line.temperature_k;
            } else {
              SpinTemperatureInputs sin;
              sin.kinetic_temperature_k = line.temperature_k;
              sin.radiation_temperature_k = t_cmb;
              sin.collisional_coupling = m_collision_table->coupling(
                  line.temperature_k, t_cmb, line.absorber_density_m3, gas.electronDensity(row).require("n_e"),
                  gas.protonDensity(row).require("n_HII"));
              sin.lyman_alpha_coupling = 0.0;  // no admitted Ly-a coupling provider (L03 spin coupling FAIL_CLOSED)
              h.spin_temperature_k = spinTemperature(sin);
            }
            hi_segments.push_back(h);
          }
        }
        if (o.sz_maps) {
          SzSegment s;
          s.path_proper_m = path;
          s.electron_density_m3 = gas.electronDensity(row).require("n_e");
          s.electron_temperature_k = gas.temperature(row).require("T");
          s.radial_velocity_m_s = v_los;  // observer at -infinity along -dir: +dir recedes
          sz_segments.push_back(s);
        }
        if (o.dispersion_map) {
          DispersionSegment d;
          d.path_proper_m = path;
          d.electron_density_m3 = gas.electronDensity(row).require("n_e");
          d.redshift = cosmological ? redshift : 0.0;
          dm_segments.push_back(d);
        }
      }
      const std::size_t pixel = j * n + i;
      if (o.lyman_alpha_forest) {
        const auto tau = lymanAlphaOpticalDepth(forest_segments, forest_axis, forest_model);
        const ForestStatistics stats = forestStatistics(tau, forest_axis, o.forest_pdf_bins, redshift, forest_segments);
        for (std::size_t p = 0; p < tau.size(); ++p) forest_flux[pixel * o.forest_pixels + p] = std::exp(-tau[p]);
        mean_flux_per_sightline.push_back(stats.mean_flux);
        if (mean_power.empty()) {
          mean_power.assign(stats.power_m_per_s.size(), 0.0);
          power_k = stats.power_k_s_per_m;
          mean_pdf.assign(stats.flux_pdf.size(), 0.0);
        }
        for (std::size_t q = 0; q < mean_power.size(); ++q) mean_power[q] += stats.power_m_per_s[q] / (n * n);
        for (std::size_t q = 0; q < mean_pdf.size(); ++q) mean_pdf[q] += stats.flux_pdf[q] / (n * n);
      }
      if (o.hi_cube) {
        const double background = cosmological ? t_cmb : 0.0;
        const HiTransferModel model = o.hi_transfer == core::HiTransferConfig::kFullOpticalDepth
                                          ? HiTransferModel::kFullOpticalDepth
                                          : HiTransferModel::kOpticallyThin;
        const auto tb = hiBrightnessTemperature(hi_segments, hi_axis, model, background);
        for (std::size_t c = 0; c < tb.size(); ++c) {
          hi_cube[(c * n + j) * n + i] = cosmological ? (tb[c] - background) / (1.0 + redshift) : tb[c];
        }
      }
      if (o.sz_maps) {
        const SzIntegrals sz = integrateSz(sz_segments);
        y_map[pixel] = sz.compton_y;
        tau_map[pixel] = sz.optical_depth;
        ksz_map[pixel] = sz.ksz_delta_t_over_t;
      }
      if (o.dispersion_map) dm_map[pixel] = integrateDispersion(dm_segments).dm_m2;
    }
  }

  const std::filesystem::path dir_out(m_config.output.output_directory);
  std::filesystem::create_directories(dir_out);
  const std::string provenance = gas.provenance() + ";a=" + std::to_string(a) + ";los_axis=" +
                                 core::losAxisConfigToString(o.los_axis) + ";sightlines_per_side=" + std::to_string(n) +
                                 ";geometry=leaf_patches;finite_volume=periodic_box_no_replication";
  const auto transverse_axes = [&](FitsImage& image) {
    for (const std::size_t t : {t1, t2}) {
      FitsAxis ax;
      ax.length = n;
      ax.ctype = "COMOV_" + std::string(1, static_cast<char>('X' + t));
      ax.cunit = "code_length";
      ax.crpix = 1.0;
      ax.crval = lo[t] + 0.5 * extent[t] / static_cast<double>(n);
      ax.cdelt = extent[t] / static_cast<double>(n);
      image.axes.push_back(ax);
    }
  };
  const auto write = [&](const std::string& product, FitsImage image) {
    image.history.push_back("CHUI-EXT-PHYS-REF-v1 " + product + " (source-only pass 01 implementation, NOT_QUALIFIED)");
    image.history.push_back(provenance);
    const std::string path = (dir_out / (fileName(m_config, product, step) + ".fits")).string();
    writeFitsImage(path, image);
    report.written_files.push_back(path);
  };
  if (o.lyman_alpha_forest) {
    FitsImage image;
    FitsAxis vel;
    vel.length = o.forest_pixels;
    vel.ctype = "VELO";
    vel.cunit = "m/s";
    vel.crval = forest_axis.center(0);
    vel.cdelt = forest_axis.width();
    image.axes.push_back(vel);
    FitsAxis los;
    los.length = static_cast<std::size_t>(n) * n;
    los.ctype = "SIGHTLIN";
    image.axes.push_back(los);
    image.data = forest_flux;
    image.bunit = "transmitted_flux";
    image.history.push_back(std::string("profile=") + core::forestProfileConfigToString(o.forest_profile));
    write("lya_flux", std::move(image));
  }
  if (o.hi_cube) {
    FitsImage image;
    transverse_axes(image);
    FitsAxis vel;
    vel.length = o.hi_channels;
    vel.ctype = "VRAD";
    vel.cunit = "m/s";
    vel.crval = hi_axis.center(0);
    vel.cdelt = hi_axis.width();
    image.axes.push_back(vel);
    image.data = hi_cube;
    image.bunit = "K";
    image.history.push_back(std::string("spin=") + core::hiSpinModelConfigToString(o.hi_spin_model) +
                            ";transfer=" + core::hiTransferConfigToString(o.hi_transfer));
    write("hi_cube", std::move(image));
  }
  if (o.sz_maps) {
    for (const auto& [name, data] : {std::pair<std::string, std::vector<double>*>{"compton_y", &y_map},
                                     std::pair<std::string, std::vector<double>*>{"tau_e", &tau_map},
                                     std::pair<std::string, std::vector<double>*>{"ksz_dt_over_t", &ksz_map}}) {
      FitsImage image;
      transverse_axes(image);
      image.data = *data;
      image.bunit = "dimensionless";
      write(name, std::move(image));
    }
  }
  if (o.dispersion_map) {
    FitsImage image;
    transverse_axes(image);
    image.data = dm_map;
    image.bunit = "m-2";
    image.history.push_back("RM not produced: no magnetic-field sample attached to the science accessor");
    write("dispersion_measure", std::move(image));
  }
  // Text summary (statistics and provenance).
  std::ostringstream summary;
  summary.precision(17);
  summary << "schema = chui_observables_summary_v1\nstep = " << step << "\nscale_factor = " << a
          << "\nprovenance = " << provenance << "\nqualification = NOT_QUALIFIED\n";
  if (o.lyman_alpha_forest) {
    double mean = 0.0;
    for (const double f : mean_flux_per_sightline) mean += f / static_cast<double>(mean_flux_per_sightline.size());
    summary << "lya_mean_flux = " << mean << "\nlya_flux_pdf =";
    for (const double v : mean_pdf) summary << ' ' << v;
    summary << "\nlya_power_k_s_per_m =";
    for (const double v : power_k) summary << ' ' << v;
    summary << "\nlya_power_m_per_s =";
    for (const double v : mean_power) summary << ' ' << v;
    summary << '\n';
  }
  const std::string summary_path = (dir_out / (fileName(m_config, "summary", step) + ".txt")).string();
  {
    std::ofstream out(summary_path + ".partial", std::ios::trunc);
    out << summary.str();
    if (!out) throw std::runtime_error("cannot write observables summary '" + summary_path + "'");
  }
  std::filesystem::rename(summary_path + ".partial", summary_path);
  report.written_files.push_back(summary_path);
  return report;
}

}  // namespace cosmosim::observables
