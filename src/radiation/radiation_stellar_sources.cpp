// CHUI-EXT-PHYS-REF-v1 F02 stellar luminosity provider (source only, NOT
// compiled or executed in pass 01).

#include "cosmosim/radiation/radiation_stellar_sources.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>

#include "core/internal/sha256.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/lightcone/ray_geometry.hpp"

namespace cosmosim::radiation {
namespace {

[[noreturn]] void fail(const std::string& origin, const std::string& message) {
  throw std::runtime_error("stellar luminosity table '" + origin + "': " + message);
}

[[nodiscard]] std::string trim(const std::string& s) {
  std::size_t b = 0;
  std::size_t e = s.size();
  while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])) ) --e;
  return s.substr(b, e - b);
}

[[nodiscard]] std::vector<double> parseList(const std::string& text, const std::string& origin) {
  std::vector<double> out;
  std::stringstream stream(text);
  std::string item;
  while (std::getline(stream, item, ',')) {
    const std::string t = trim(item);
    char* end = nullptr;
    const double v = std::strtod(t.c_str(), &end);
    if (t.empty() || end != t.c_str() + static_cast<std::ptrdiff_t>(t.size()) || !std::isfinite(v)) {
      fail(origin, "malformed list value '" + item + "'");
    }
    out.push_back(v);
  }
  return out;
}

[[nodiscard]] std::size_t bracket(const std::vector<double>& grid, double x, const char* what, const std::string& origin) {
  if (!(x >= grid.front()) || !(x <= grid.back())) {
    fail(origin, std::string(what) + " " + std::to_string(x) + " outside the tabulated range [" +
                     std::to_string(grid.front()) + ", " + std::to_string(grid.back()) + "] (no extrapolation)");
  }
  const auto it = std::upper_bound(grid.begin(), grid.end(), x);
  std::size_t hi = static_cast<std::size_t>(std::distance(grid.begin(), it));
  hi = std::clamp<std::size_t>(hi, 1U, grid.size() - 1U);
  return hi - 1U;
}

}  // namespace

namespace {

[[nodiscard]] double interpolateLane(const StellarLuminosityTable& t, const std::vector<double>& lane, std::size_t g,
                                     double age, double z) {
  const std::size_t nz = t.metallicity.size();
  const std::size_t ng = t.groupCount();
  if (g >= ng) throw std::out_of_range("stellar luminosity group index");
  const std::size_t ia = bracket(t.age_yr, age, "age_yr", t.sha256);
  const double ta = (age - t.age_yr[ia]) / (t.age_yr[ia + 1U] - t.age_yr[ia]);
  std::size_t iz = 0;
  double tz = 0.0;
  if (nz > 1U) {
    iz = bracket(t.metallicity, z, "metallicity", t.sha256);
    tz = (std::log10(z) - std::log10(t.metallicity[iz])) /
         (std::log10(t.metallicity[iz + 1U]) - std::log10(t.metallicity[iz]));
  } else if (std::abs(z - t.metallicity.front()) > 1.0e-12 * t.metallicity.front()) {
    fail(t.sha256, "single-metallicity table queried at a different Z (no extrapolation)");
  }
  const auto at = [&](std::size_t a, std::size_t zz) { return lane[(a * nz + zz) * ng + g]; };
  const std::size_t iz1 = nz > 1U ? iz + 1U : iz;
  const double lo = (1.0 - tz) * at(ia, iz) + tz * at(ia, iz1);
  const double hi = (1.0 - tz) * at(ia + 1U, iz) + tz * at(ia + 1U, iz1);
  return (1.0 - ta) * lo + ta * hi;
}

}  // namespace

double StellarLuminosityTable::luminosity(std::size_t g, double age, double z) const {
  return interpolateLane(*this, l_erg_s_per_msun, g, age, z);
}

double StellarLuminosityTable::photonRate(std::size_t g, double age, double z) const {
  return interpolateLane(*this, q_photons_s_per_msun, g, age, z);
}

double StellarLuminosityTable::integratedEnergy(std::size_t g, double age0, double age1, double z) const {
  if (!(age1 >= age0)) throw std::invalid_argument("stellar emission interval must have age1 >= age0");
  if (age1 == age0) return 0.0;
  // Piecewise-linear interpolant in age: trapezoid over the union of nodes is exact.
  std::vector<double> nodes{age0};
  for (const double a : age_yr) {
    if (a > age0 && a < age1) nodes.push_back(a);
  }
  nodes.push_back(age1);
  double sum = 0.0;
  for (std::size_t i = 0; i + 1U < nodes.size(); ++i) {
    sum += 0.5 * (nodes[i + 1U] - nodes[i]) * (luminosity(g, nodes[i], z) + luminosity(g, nodes[i + 1U], z));
  }
  return sum * core::constants::k_seconds_per_megayear * 1.0e-6;  // yr -> s: erg per Msun
}

StellarLuminosityTable parseStellarLuminosityTable(const std::string& bytes, const std::string& origin) {
  StellarLuminosityTable table;
  std::map<std::string, std::string> meta;
  struct Row {
    double age;
    double z;
    std::size_t g;
    double l;
    double q;
  };
  std::vector<Row> rows;
  std::istringstream stream(bytes);
  std::string line;
  while (std::getline(stream, line)) {
    const std::string t = trim(line);
    if (t.empty()) continue;
    if (t.front() == '#') {
      const auto colon = t.find(':');
      if (colon != std::string::npos) meta[trim(t.substr(1, colon - 1))] = trim(t.substr(colon + 1));
      continue;
    }
    std::istringstream row_stream(t);
    Row row{};
    double g = 0.0;
    row_stream >> row.age >> row.z >> g >> row.l >> row.q;
    if (row_stream.fail() || !std::isfinite(row.age) || !std::isfinite(row.l) || !std::isfinite(row.q) ||
        row.l < 0.0 || row.q < 0.0 || g < 0.0 || g != std::floor(g)) {
      fail(origin, "malformed data row '" + t + "'");
    }
    row.g = static_cast<std::size_t>(g);
    rows.push_back(row);
  }
  const auto need = [&](const char* key) {
    const auto it = meta.find(key);
    if (it == meta.end() || it->second.empty()) fail(origin, std::string("metadata '") + key + "' missing");
    return it->second;
  };
  if (need("schema") != k_ssp_group_luminosity_schema) fail(origin, "schema is not chui_ssp_group_luminosity_v1");
  table.source = need("source");
  table.imf = need("imf");
  table.group_edges_ev = parseList(need("group_edges_ev"), origin);
  table.metallicity = parseList(need("metallicities"), origin);
  if (need("columns") != "age_yr metallicity group L_erg_s_per_msun Q_photons_s_per_msun") {
    fail(origin, "columns must be 'age_yr metallicity group L_erg_s_per_msun Q_photons_s_per_msun'");
  }
  if (table.group_edges_ev.size() < 2U) fail(origin, "need at least one group");
  for (std::size_t i = 1; i < table.group_edges_ev.size(); ++i) {
    if (!(table.group_edges_ev[i] > table.group_edges_ev[i - 1U])) fail(origin, "group edges must increase");
  }
  for (std::size_t i = 0; i < table.metallicity.size(); ++i) {
    if (!(table.metallicity[i] > 0.0) || (i > 0 && !(table.metallicity[i] > table.metallicity[i - 1U]))) {
      fail(origin, "metallicities must be > 0 and strictly increasing");
    }
  }
  for (const Row& row : rows) {
    if (table.age_yr.empty() || row.age > table.age_yr.back()) {
      table.age_yr.push_back(row.age);
    } else if (row.age < table.age_yr.back()) {
      fail(origin, "rows must be ordered by ascending age");
    }
  }
  if (table.age_yr.size() < 2U || table.age_yr.front() != 0.0) {
    fail(origin, "ages must start at 0 and contain at least two nodes (newborn populations without extrapolation)");
  }
  const std::size_t na = table.age_yr.size();
  const std::size_t nz = table.metallicity.size();
  const std::size_t ng = table.groupCount();
  if (rows.size() != na * nz * ng) fail(origin, "the (age, Z, group) tensor grid is incomplete");
  table.l_erg_s_per_msun.assign(na * nz * ng, -1.0);
  table.q_photons_s_per_msun.assign(na * nz * ng, -1.0);
  for (const Row& row : rows) {
    const std::size_t ia = static_cast<std::size_t>(
        std::distance(table.age_yr.begin(), std::lower_bound(table.age_yr.begin(), table.age_yr.end(), row.age)));
    const auto zit = std::find_if(table.metallicity.begin(), table.metallicity.end(),
                                  [&](double z) { return std::abs(z - row.z) <= 1.0e-12 * z; });
    if (zit == table.metallicity.end() || row.g >= ng) fail(origin, "row metallicity/group not on the declared grid");
    const std::size_t iz = static_cast<std::size_t>(std::distance(table.metallicity.begin(), zit));
    const std::size_t index = (ia * nz + iz) * ng + row.g;
    if (table.l_erg_s_per_msun[index] >= 0.0) fail(origin, "duplicated (age, Z, group) row");
    table.l_erg_s_per_msun[index] = row.l;
    table.q_photons_s_per_msun[index] = row.q;
  }
  table.sha256 = core::internal::sha256Hex(bytes);
  return table;
}

StellarLuminosityTable loadStellarLuminosityTable(const std::string& path, const std::string& expected_sha256) {
  if (path.empty()) fail(path, "empty path; no stellar population synthesis table is bundled or synthesized");
  std::ifstream input(path, std::ios::binary);
  if (!input) fail(path, "cannot open file");
  const std::string bytes((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  const std::string sha = core::internal::sha256Hex(bytes);
  if (sha != expected_sha256) fail(path, "SHA-256 mismatch (expected " + expected_sha256 + ", got " + sha + ")");
  return parseStellarLuminosityTable(bytes, path);
}

void requireMatchingGroupEdges(const StellarLuminosityTable& table, const RadiationSpectralDefinition& definition) {
  if (definition.grey || table.group_edges_ev.size() != definition.edges_ev.size()) {
    fail(table.sha256, "group count differs from the RT spectral definition (" + definition.profile + ")");
  }
  for (std::size_t i = 0; i < table.group_edges_ev.size(); ++i) {
    if (std::abs(table.group_edges_ev[i] - definition.edges_ev[i]) > 1.0e-9 * definition.edges_ev[i]) {
      fail(table.sha256, "group edge " + std::to_string(i) + " differs from the RT spectral definition");
    }
  }
}

StellarRadiationSourceReport buildStellarRadiationSources(const core::SimulationState& state,
                                                          const StellarLuminosityTable& table,
                                                          const RadiationSpectralDefinition& definition,
                                                          const StellarSourceInputs& in,
                                                          const physics::ProperFrameContext& frame,
                                                          const physics::PeriodicBox& box,
                                                          std::vector<RadiationSourceDeposit>& deposits) {
  StellarRadiationSourceReport report;
  if (in.background == nullptr) throw std::invalid_argument("stellar sources need the cosmology background authority");
  if (!(in.dt_proper_code > 0.0)) return report;
  requireMatchingGroupEdges(table, definition);
  const core::RadiationMomentSoa& rad = state.radiation;
  if (rad.patches.empty() || rad.group_count != table.groupCount()) {
    throw std::runtime_error("stellar sources: radiation records are absent or have a different group count");
  }
  // Bucket-indexed leaf locator over the radiation records (record-local rows).
  std::vector<lightcone::LeafPatch> leaves;
  leaves.reserve(rad.patches.size());
  for (std::uint32_t r = 0; r < rad.patches.size(); ++r) {
    const core::RadiationPatchRecord& rr = rad.patches[r];
    lightcone::LeafPatch leaf;
    leaf.patch_id = rr.patch_id;
    leaf.level = rr.level;
    leaf.origin = {rr.origin_x_comoving, rr.origin_y_comoving, rr.origin_z_comoving};
    leaf.extent = {rr.cell_width_x_comoving * rr.nx, rr.cell_width_y_comoving * rr.ny, rr.cell_width_z_comoving * rr.nz};
    leaf.dims = {rr.nx, rr.ny, rr.nz};
    leaf.first_cell = 0U;
    leaves.push_back(leaf);
  }
  std::array<double, 3> domain_lo{1e300, 1e300, 1e300};
  std::array<double, 3> domain_hi{-1e300, -1e300, -1e300};
  for (const auto& leaf : leaves) {
    for (std::size_t a = 0; a < 3U; ++a) {
      domain_lo[a] = std::min(domain_lo[a], leaf.origin[a]);
      domain_hi[a] = std::max(domain_hi[a], leaf.origin[a] + leaf.extent[a]);
    }
  }
  const lightcone::LeafGeometry locator(std::move(leaves), domain_lo,
                                        {domain_hi[0] - domain_lo[0], domain_hi[1] - domain_lo[1],
                                         domain_hi[2] - domain_lo[2]},
                                        lightcone::RayBoundary::kOpen);
  const double erg_to_code = 1.0e-7 / (in.units.mass_si_per_code * in.units.velocity_si_per_code *
                                       in.units.velocity_si_per_code);
  const double msun_per_code_mass = in.units.mass_si_per_code / core::constants::k_solar_mass_si;
  const double a4 = std::pow(in.scale_factor_end, 4);
  const double year_s = core::constants::k_seconds_per_megayear * 1.0e-6;
  const std::size_t ng = table.groupCount();
  for (std::size_t s = 0; s < state.star_particles.particle_index.size(); ++s) {
    const double a_form = state.star_particles.formation_scale_factor[s];
    if (!(a_form > 0.0) || a_form > in.scale_factor_end) continue;  // not yet formed in this interval
    const double a_start = std::max(a_form, in.scale_factor_begin);
    const double age0 = in.background->cosmicTimeIntervalSi(a_form, a_start) / year_s;
    const double age1 = in.background->cosmicTimeIntervalSi(a_form, in.scale_factor_end) / year_s;
    const double z = state.star_particles.metallicity_mass_fraction[s];
    const double mass_msun = state.star_particles.birth_mass_code[s] * msun_per_code_mass;
    const std::uint32_t p = state.star_particles.particle_index[s];
    const std::array<double, 3> x = box.wrap({state.particles.position_x_comoving[p], state.particles.position_y_comoving[p],
                                              state.particles.position_z_comoving[p]});
    lightcone::LeafCellRef ref;
    if (!locator.locate(x, ref)) {
      throw std::runtime_error("stellar source outside every radiation record (records must tile the gas domain)");
    }
    const core::RadiationPatchRecord& rr = rad.patches[ref.leaf];
    const double v_phys = frame.properVolume(rr.cell_width_x_comoving * rr.cell_width_y_comoving * rr.cell_width_z_comoving);
    ++report.sources;
    for (std::size_t g = 0; g < ng; ++g) {
      const double energy_erg = mass_msun * table.integratedEnergy(g, age0, age1, z);
      if (!(energy_erg > 0.0)) continue;
      const double energy_code = energy_erg * erg_to_code;
      const double rate_phys = energy_code / in.dt_proper_code / v_phys;
      deposits.push_back(RadiationSourceDeposit{ref.leaf, ref.cell_row, static_cast<std::uint32_t>(g), rate_phys * a4});
      report.emitted_energy_code += energy_code;
      const double mean_photon_erg = definition.groups[g].mean_photon_energy_ev * k_ev_erg;
      report.transported_photon_equivalent += energy_erg / mean_photon_erg;
      // Diagnostic photon integral of the table (trapezoid at the interval ends).
      report.table_photons += 0.5 * (age1 - age0) * year_s * mass_msun *
                              (table.photonRate(g, age0, z) + table.photonRate(g, age1, z));
    }
  }
  return report;
}

}  // namespace cosmosim::radiation
