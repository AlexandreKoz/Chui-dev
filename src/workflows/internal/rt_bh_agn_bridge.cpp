#include "workflows/internal/rt_bh_agn_bridge.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

#include "cosmosim/analysis/halo_workflow.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/physics/black_hole/bh_seed_provider.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/radiation/radiation_sources.hpp"

namespace cosmosim::workflows::internal {
namespace {

template <typename T>
void appendRecord(std::vector<std::uint8_t>& out, const T& record) {
  const auto* bytes = reinterpret_cast<const std::uint8_t*>(&record);
  out.insert(out.end(), bytes, bytes + sizeof(T));
}

template <typename T>
std::vector<T> decodeRecords(const std::vector<std::uint8_t>& bytes) {
  if (bytes.size() % sizeof(T) != 0U) throw std::runtime_error("BH support exchange: truncated wire payload");
  std::vector<T> out(bytes.size() / sizeof(T));
  if (!out.empty()) std::memcpy(out.data(), bytes.data(), bytes.size());
  return out;
}

[[nodiscard]] double kpcToCode(double kpc, const core::UnitSystem& units) {
  return kpc * core::constants::k_kiloparsec_si / units.length_si_per_code;
}

}  // namespace

MpiBhSupportExchange::MpiBhSupportExchange(const parallel::MpiContext& mpi, const core::SimulationState& state,
                                           std::span<const std::uint8_t> mask, const physics::PeriodicBox& box)
    : m_mpi(mpi), m_box(box) {
  std::array<double, 7> local{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                              std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
                              -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), 0.0};
  const physics::BhOwnedLeafView leaves{mask};
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    if (!leaves.owned(row)) continue;
    std::array<double, 3> hw{};
    if (!physics::bhCellHalfWidthStored(state, static_cast<std::uint32_t>(row), hw)) continue;
    const double x[3] = {state.cells.center_x_comoving[row], state.cells.center_y_comoving[row],
                         state.cells.center_z_comoving[row]};
    for (int a = 0; a < 3; ++a) {
      local[a] = std::min(local[a], x[a] - hw[a]);
      local[3 + a] = std::max(local[3 + a], x[a] + hw[a]);
    }
    local[6] = 1.0;
  }
  // Stars contribute to the torque decomposition: include them in the box.
  for (std::size_t s = 0; s < state.star_particles.size(); ++s) {
    const std::uint32_t p = state.star_particles.particle_index[s];
    if (p >= state.particles.size()) continue;
    const double x[3] = {state.particles.position_x_comoving[p], state.particles.position_y_comoving[p],
                         state.particles.position_z_comoving[p]};
    for (int a = 0; a < 3; ++a) {
      local[a] = std::min(local[a], x[a]);
      local[3 + a] = std::max(local[3 + a], x[a]);
    }
    local[6] = 1.0;
  }
  std::vector<std::uint8_t> bytes;
  appendRecord(bytes, local);
  const std::vector<std::uint8_t> all = m_mpi.allgatherBytesBounded(bytes);
  m_rank_boxes = decodeRecords<std::array<double, 7>>(all);
  if (m_rank_boxes.size() != static_cast<std::size_t>(m_mpi.worldSize())) {
    throw std::runtime_error("BH support exchange: rank box allgather size mismatch");
  }
}

bool MpiBhSupportExchange::sphereIntersectsRank(const physics::BhSupportQuery& q, int rank) const {
  const auto& b = m_rank_boxes[static_cast<std::size_t>(rank)];
  if (b[6] == 0.0) return false;
  double d2 = 0.0;
  for (int a = 0; a < 3; ++a) {
    const double lo = b[a];
    const double hi = b[3 + a];
    double best = std::numeric_limits<double>::max();
    const int images = (m_box.periodic && m_box.extent[a] > 0.0) ? 1 : 0;
    for (int img = -images; img <= images; ++img) {
      const double x = q.position_stored[a] + img * m_box.extent[a];
      const double d = x < lo ? lo - x : (x > hi ? x - hi : 0.0);
      best = std::min(best, d);
    }
    d2 += best * best;
  }
  return d2 <= q.support_radius_stored * q.support_radius_stored;
}

std::vector<physics::BhSupportQuery> MpiBhSupportExchange::routeQueries(std::span<const physics::BhSupportQuery> q) {
  const int world = m_mpi.worldSize();
  std::vector<std::vector<std::uint8_t>> send(static_cast<std::size_t>(world));
  for (const physics::BhSupportQuery& query : q) {
    for (int r = 0; r < world; ++r) {
      if (sphereIntersectsRank(query, r)) appendRecord(send[static_cast<std::size_t>(r)], query);
    }
  }
  const auto received = parallel::exchangeBoundedAlltoallBytes(m_mpi, send);
  m_last_routed.clear();
  for (const auto& bytes : received) {
    const auto records = decodeRecords<physics::BhSupportQuery>(bytes);
    m_last_routed.insert(m_last_routed.end(), records.begin(), records.end());
  }
  std::sort(m_last_routed.begin(), m_last_routed.end(), [](const auto& a, const auto& b) {
    return std::tie(a.origin_rank, a.query_index) < std::tie(b.origin_rank, b.query_index);
  });
  return m_last_routed;
}

std::vector<physics::BhSupportReply> MpiBhSupportExchange::returnReplies(std::span<const physics::BhSupportReply> r) {
  if (r.size() != m_last_routed.size()) {
    throw std::logic_error("BH support exchange: reply count differs from the routed query count");
  }
  std::vector<std::vector<std::uint8_t>> send(static_cast<std::size_t>(m_mpi.worldSize()));
  for (std::size_t i = 0; i < r.size(); ++i) {
    appendRecord(send[m_last_routed[i].origin_rank], r[i]);
  }
  const auto received = parallel::exchangeBoundedAlltoallBytes(m_mpi, send);
  std::vector<physics::BhSupportReply> out;
  for (const auto& bytes : received) {
    const auto records = decodeRecords<physics::BhSupportReply>(bytes);
    out.insert(out.end(), records.begin(), records.end());
  }
  std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) {
    return std::tie(a.query_index, a.responder_rank) < std::tie(b.query_index, b.responder_rank);
  });
  m_last_routed.clear();
  return out;
}

bool MpiBhSupportExchange::allAgree(bool local_ok) { return m_mpi.allreduceMinUint64(local_ok ? 1U : 0U) == 1U; }

std::uint32_t MpiBhSupportExchange::worldRank() const noexcept {
  return static_cast<std::uint32_t>(m_mpi.worldRank());
}

RtBhAgnSourceBridge::RtBhAgnSourceBridge(const core::SimulationConfig& config, const core::UnitSystem& units,
                                         const parallel::MpiContext& mpi, std::array<double, 3> box, bool periodic)
    : m_config(config), m_units(units), m_mpi(mpi) {
  m_box.extent = box;
  m_box.periodic = periodic;
  const core::PhysicsConfig& p = config.physics;
  if (p.bh_engine.model == core::BhEngineModelConfig::kReferenceV1) {
    const double soft_kpc = config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0
        ? config.numerics.gravity_softening_black_hole_kpc_comoving
        : config.numerics.gravity_softening_kpc_comoving;
    m_engine.emplace(physics::makeBlackHoleEngineParameters(p, units, kpcToCode(soft_kpc, units), m_box, 1.0e-20));
  }
  if (p.agn.model == core::AgnFeedbackModelConfig::kChannelsV1) m_agn.emplace(p.agn);
  if (p.rt.enabled) {
    radiation::RadiationRuntimeSettings rs;
    rs.config = p.rt;
    radiation::SpectralShape shape;
    switch (p.rt.source_shape) {
      case core::RtSourceShapeConfig::kSedTable:
        if (p.rt.agn_sources) {
          m_agn_sed = std::make_unique<radiation::TabulatedSpectrum>(radiation::loadTabulatedSpectrum(
              p.rt.agn_sed_table_path, p.rt.agn_sed_table_sha256, "chui_agn_sed_table_v1"));
          shape.kind = radiation::SpectralShapeKind::kTable;
          shape.table = m_agn_sed.get();
        } else {
          shape.kind = radiation::SpectralShapeKind::kFlat;  // no sources: shape only sets quadrature weights
        }
        break;
      case core::RtSourceShapeConfig::kBlackbodyReference:
        shape.kind = radiation::SpectralShapeKind::kBlackbody;
        shape.blackbody_temperature_k = p.rt.source_blackbody_temperature_k;
        break;
      case core::RtSourceShapeConfig::kPowerLawReference:
        shape.kind = radiation::SpectralShapeKind::kPowerLaw;
        shape.power_law_index = p.rt.source_power_law_index;
        break;
    }
    const std::string profile = core::rtSpectralProfileConfigToString(p.rt.spectral_profile);
    rs.spectrum = radiation::buildRadiationSpectralDefinition(profile, p.rt.spectrum_upper_edge_ev, shape);
    m_agn_group_fractions = radiation::radiationGroupLuminosityFractions(rs.spectrum, shape);
    const double c_code = core::constants::k_speed_of_light_si / units.velocity_si_per_code;
    rs.coupling.c_code = c_code;
    rs.coupling.c_hat_code = p.rt.c_fraction * c_code;
    const double mass_g = units.mass_si_per_code * 1.0e3;
    const double length_cm = units.length_si_per_code * 1.0e2;
    rs.coupling.kappa_absorption_code = p.rt.constant_absorption_cm2_g * mass_g / (length_cm * length_cm);
    rs.coupling.kappa_scattering_code = p.rt.constant_scattering_cm2_g * mass_g / (length_cm * length_cm);
    const double energy_density_si_per_code = units.mass_si_per_code * units.velocity_si_per_code *
        units.velocity_si_per_code / std::pow(units.length_si_per_code, 3.0);
    rs.coupling.radiation_constant_code = physics::k_radiation_constant_si / energy_density_si_per_code;
    rs.coupling.length_cm_per_code = length_cm;
    rs.coupling.energy_density_cgs_per_code = energy_density_si_per_code * 10.0;
    rs.coupling.density_cgs_per_code = units.densityCodeToCgs(1.0);
    rs.coupling.relative_tolerance = p.rt.newton_relative_tolerance;
    rs.coupling.max_iterations = p.rt.newton_max_iterations;
    rs.transport.riemann = p.rt.riemann_solver;
    rs.transport.reconstruction = p.rt.reconstruction;
    rs.transport.boundary = p.rt.boundary;
    rs.transport.domain_extent_comoving = box;
    rs.transport.c_hat_code = rs.coupling.c_hat_code;
    if (p.rt.opacity_model == core::RtOpacityModelConfig::kConstantReference && rs.spectrum.grey) {
      rs.coupling_model = radiation::RadiationCouplingModel::kGreyThermalReference;
    } else if (p.rt.opacity_model == core::RtOpacityModelConfig::kPhotoionization) {
      rs.coupling_model = radiation::RadiationCouplingModel::kPhotoionization;
    }
    m_radiation.emplace(std::move(rs));
  }
}

std::vector<physics::BhHaloSeedCandidate> RtBhAgnSourceBridge::serialHaloCandidates(const core::SimulationState& state,
                                                                                   const RtBhAgnEpochInputs& in) const {
  if (m_mpi.isEnabled() && m_mpi.worldSize() > 1) {
    throw std::runtime_error("halo_fof_serial seeding is a single-rank provider (distributed seeding FAIL_CLOSED)");
  }
  analysis::FofConfig fof;
  fof.linking_length_factor_mean_interparticle = m_config.analysis.halo_fof_linking_length_factor;
  fof.min_group_size = static_cast<std::uint64_t>(std::max(m_config.analysis.halo_fof_min_group_size, 1));
  fof.include_gas = m_config.analysis.halo_include_gas;
  fof.include_stars = m_config.analysis.halo_include_stars;
  fof.include_black_holes = true;  // occupancy needs BH membership
  const analysis::FofHaloFinder finder(fof);
  const analysis::HaloCatalog catalog = finder.buildCatalog(state, m_config, in.step_index, in.scale_factor_end);
  std::unordered_map<std::uint64_t, std::size_t> index;
  std::vector<physics::BhHaloSeedCandidate> out;
  std::vector<double> radius2;
  for (const analysis::HaloCatalogEntry& h : catalog.halos) {
    physics::BhHaloSeedCandidate c;
    c.halo_key = h.min_particle_id;
    c.halo_mass_code = h.total_mass_code;
    c.center_stored = h.center_of_mass_comov;
    index[h.halo_id] = out.size();
    out.push_back(c);
    radius2.push_back(0.0);
  }
  for (std::size_t pi = 0; pi < catalog.halo_id_by_particle.size() && pi < state.particles.size(); ++pi) {
    const auto it = index.find(catalog.halo_id_by_particle[pi]);
    if (it == index.end()) continue;
    physics::BhHaloSeedCandidate& c = out[it->second];
    if (state.particleSpeciesTag(pi) == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) c.occupied = true;
    const auto d = m_box.delta(c.center_stored, {state.particles.position_x_comoving[pi],
                                                  state.particles.position_y_comoving[pi],
                                                  state.particles.position_z_comoving[pi]});
    radius2[it->second] = std::max(radius2[it->second], d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  }
  // Donor: densest owned leaf gas cell within the member radius (ties: lower id).
  std::vector<physics::BhSupportQuery> queries(out.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    queries[i].position_stored = out[i].center_stored;
    queries[i].support_radius_stored = std::sqrt(radius2[i]);
    queries[i].query_index = i;
  }
  std::vector<double> best(out.size(), -1.0);
  if (!queries.empty()) {
    const physics::BhProximityFilter filter(queries, m_box);
    const physics::BhOwnedLeafView leaves{in.owned_leaf_mask};
    std::vector<std::uint32_t> hits;
    for (std::size_t row = 0; row < state.cells.size(); ++row) {
      if (!leaves.owned(row)) continue;
      filter.matches({state.cells.center_x_comoving[row], state.cells.center_y_comoving[row],
                      state.cells.center_z_comoving[row]}, hits);
      for (const std::uint32_t qi : hits) {
        const double rho = state.gas_cells.density_code[row];
        const std::uint64_t id = state.gas_cells.gas_cell_id[row];
        if (rho > best[qi] || (rho == best[qi] && id < out[qi].donor_gas_cell_id)) {
          best[qi] = rho;
          out[qi].donor_gas_cell_id = id;
        }
      }
    }
  }
  return out;
}

RtBhAgnEpochReport RtBhAgnSourceBridge::executeDiscreteEpoch(core::SimulationState& state, const RtBhAgnEpochInputs& in) {
  RtBhAgnEpochReport report;
  std::unique_ptr<physics::BhSupportExchange> exchange;
  if (m_mpi.isEnabled() && m_mpi.worldSize() > 1) {
    exchange = std::make_unique<MpiBhSupportExchange>(m_mpi, state, in.owned_leaf_mask, m_box);
  } else {
    exchange = std::make_unique<physics::LocalBhSupportExchange>();
  }
  const physics::ProperFrameContext frame = in.cosmological
      ? physics::ProperFrameContext::cosmological(m_units, in.scale_factor_end, in.hubble_end_si)
      : physics::ProperFrameContext::isolated(m_units);
  const double dt_proper = in.physical_seconds / m_units.timeSiPerCode();
  if (!(dt_proper > 0.0)) throw std::runtime_error("RT/BH/AGN epoch requires a positive physical interval");
  const physics::BhOwnedLeafView leaves{in.owned_leaf_mask};
  auto merge_rows = [&](const std::vector<std::uint32_t>& rows) {
    report.touched_gas_rows.insert(report.touched_gas_rows.end(), rows.begin(), rows.end());
  };
  if (m_engine.has_value()) {
    physics::BhIntervalContext ctx;
    ctx.frame = frame;
    ctx.time_begin_code = in.time_begin_code;
    ctx.time_end_code = in.time_end_code;
    ctx.dt_proper_code = dt_proper;
    ctx.interval_index = in.interval_index;
    ctx.leaves = leaves;
    ctx.particle_scheduler = in.particle_scheduler;
    report.engine = m_engine->advanceInterval(state, ctx, *exchange);
    report.gas_mutated = report.gas_mutated || report.engine.gas_mutated;
    report.particle_mass_mutated = report.particle_mass_mutated || report.engine.particle_mass_mutated;
    merge_rows(report.engine.touched_gas_rows);
    report.removed_particle_ids = report.engine.removed_particle_ids;
    const core::BlackHoleEngineConfig& bh = m_config.physics.bh_engine;
    if (bh.seed_model == core::BhSeedModelConfig::kHaloFofSerial && bh.seed_cadence_steps > 0U &&
        in.step_index % bh.seed_cadence_steps == 0U) {
      const auto candidates = serialHaloCandidates(state, in);
      physics::BhSeedingParameters sp;
      sp.halo_mass_threshold_code = m_config.physics.bh_seed_halo_mass_threshold_code;
      sp.seed_mass_code = m_config.physics.bh_seed_mass_code;
      physics::LocalParticleIdRegistry local_registry;
      physics::ParticleIdPrecommit& ids = in.id_precommit != nullptr ? *in.id_precommit : local_registry;
      const physics::BhSeedingReport seeds =
          physics::applyHaloSeeding(state, candidates, sp, in.step_index, in.interval_index, in.time_end_code, ids);
      report.seeded = seeds.seeded;
      report.seeded_particle_ids = seeds.seeded_particle_ids;
      merge_rows(seeds.touched_gas_rows);
      if (seeds.seeded > 0U) {
        report.gas_mutated = true;
        report.particle_mass_mutated = true;
      }
    }
  }
  if (m_agn.has_value()) {
    physics::AgnInjectionContext actx;
    actx.frame = frame;
    actx.time_end_code = in.time_end_code;
    actx.interval_index = in.interval_index;
    actx.leaves = leaves;
    actx.box = m_box;
    report.agn = m_agn->apply(state, actx, *exchange);
    report.gas_mutated = report.gas_mutated || report.agn.gas_mutated;
    report.particle_mass_mutated = true;
    merge_rows(report.agn.touched_gas_rows);
  }
  if (m_radiation.has_value()) {
    const std::uint32_t rank = static_cast<std::uint32_t>(m_mpi.worldRank());
    if (state.radiation.empty()) {
      radiation::allocateRadiationState(state, m_radiation->settings().spectrum, m_config.physics.rt.c_fraction,
                                        in.scale_factor_begin, rank);
    }
    std::string reason;
    if (!radiation::radiationRecordsMatchPatches(state, rank, &reason)) {
      throw std::runtime_error("RT records no longer match the leaf patch tiling: " + reason +
                               " (regrid/migration of radiation moments is not implemented; FAIL_CLOSED)");
    }
    radiation::RadiationIntervalInputs rin;
    rin.dt_proper_code = dt_proper;
    rin.scale_factor_begin = in.scale_factor_begin;
    rin.scale_factor_end = in.scale_factor_end;
    rin.hubble_begin_code = in.hubble_begin_si * m_units.timeSiPerCode();
    rin.hubble_end_code = in.hubble_end_si * m_units.timeSiPerCode();
    rin.interval_index = in.interval_index;
    rin.frame_end = frame;
    if (m_config.physics.rt.agn_sources) {
      (void)radiation::buildAgnRadiationSources(state, m_agn_group_fractions, frame, m_box, dt_proper, rin.sources);
    }
    report.radiation = m_radiation->advanceInterval(state, rin);
    merge_rows(report.radiation.touched_gas_rows);
    report.gas_mutated = report.gas_mutated || !report.radiation.touched_gas_rows.empty();
  }
  std::sort(report.touched_gas_rows.begin(), report.touched_gas_rows.end());
  report.touched_gas_rows.erase(std::unique(report.touched_gas_rows.begin(), report.touched_gas_rows.end()),
                                report.touched_gas_rows.end());
  return report;
}

}  // namespace cosmosim::workflows::internal
