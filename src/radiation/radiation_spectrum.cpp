#include "cosmosim/radiation/radiation_spectrum.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "core/internal/sha256.hpp"

namespace cosmosim::radiation {
namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

// Gauss-Legendre nodes/weights on [-1, 1] via Newton on P_n.
void gaussLegendre(std::uint32_t n, std::vector<double>& x, std::vector<double>& w) {
  x.assign(n, 0.0);
  w.assign(n, 0.0);
  for (std::uint32_t i = 0; i < n; ++i) {
    double z = std::cos(k_pi * (i + 0.75) / (n + 0.5));
    double pp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p1 = 1.0;
      double p2 = 0.0;
      for (std::uint32_t j = 0; j < n; ++j) {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      const double dz = p1 / pp;
      z -= dz;
      if (std::abs(dz) < 1.0e-15) break;
    }
    x[i] = -z;
    w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
  }
}

[[nodiscard]] double hydrogenic(double e_ev, double threshold_ev, double sigma0) noexcept {
  if (e_ev < threshold_ev) return 0.0;
  const double ratio = e_ev / threshold_ev;
  if (ratio - 1.0 < 1.0e-12) return sigma0;
  const double eps = std::sqrt(ratio - 1.0);
  return sigma0 * std::pow(1.0 / ratio, 4.0) * std::exp(4.0 - 4.0 * std::atan(eps) / eps) /
      (1.0 - std::exp(-2.0 * k_pi / eps));
}

[[nodiscard]] double vernerHeI(double e_ev) noexcept {
  if (e_ev < k_ionization_potential_hei_ev) return 0.0;
  constexpr double e0 = 13.61, sigma0_mb = 949.2, ya = 1.469, p = 3.188, yw = 2.039, y0 = 0.4434, y1 = 2.136;
  const double x = e_ev / e0 - y0;
  const double y = std::sqrt(x * x + y1 * y1);
  const double f = ((x - 1.0) * (x - 1.0) + yw * yw) * std::pow(y, 0.5 * p - 5.5) * std::pow(1.0 + std::sqrt(y / ya), -p);
  return sigma0_mb * 1.0e-18 * f;
}

}  // namespace

double ionizationPotentialEv(RadiationAbsorber a) noexcept {
  switch (a) {
    case RadiationAbsorber::kHI: return k_ionization_potential_hi_ev;
    case RadiationAbsorber::kHeI: return k_ionization_potential_hei_ev;
    case RadiationAbsorber::kHeII: return k_ionization_potential_heii_ev;
  }
  return 0.0;
}

double photoionizationCrossSectionCm2(RadiationAbsorber a, double e) noexcept {
  switch (a) {
    case RadiationAbsorber::kHI: return hydrogenic(e, k_ionization_potential_hi_ev, 6.30e-18);
    case RadiationAbsorber::kHeII: return hydrogenic(e, k_ionization_potential_heii_ev, 6.30e-18 / 4.0);
    case RadiationAbsorber::kHeI: return vernerHeI(e);
  }
  return 0.0;
}

double TabulatedSpectrum::evaluate(double e) const {
  if (energy_ev.size() < 2U || !(e >= energy_ev.front()) || !(e <= energy_ev.back())) {
    throw std::domain_error("TabulatedSpectrum: photon energy outside the tabulated domain (no extrapolation)");
  }
  const auto it = std::upper_bound(energy_ev.begin(), energy_ev.end(), e);
  const std::size_t i = std::min<std::size_t>(static_cast<std::size_t>(it - energy_ev.begin()), energy_ev.size() - 1U);
  const std::size_t lo = i - 1U;
  const double y0 = l_nu[lo];
  const double y1 = l_nu[i];
  if (!(y0 > 0.0) || !(y1 > 0.0)) {
    const double t = (e - energy_ev[lo]) / (energy_ev[i] - energy_ev[lo]);
    return y0 + t * (y1 - y0);
  }
  const double t = std::log(e / energy_ev[lo]) / std::log(energy_ev[i] / energy_ev[lo]);
  return std::exp(std::log(y0) + t * (std::log(y1) - std::log(y0)));
}

TabulatedSpectrum loadTabulatedSpectrum(const std::string& path, const std::string& expected_sha256,
                                        const std::string& expected_schema) {
  if (path.empty()) throw std::runtime_error("spectral table path is empty; no SED is bundled or synthesized");
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open spectral table '" + path + "'");
  const std::string bytes((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  const std::string sha = core::internal::sha256Hex(bytes);
  if (sha != expected_sha256) {
    throw std::runtime_error("spectral table '" + path + "' SHA-256 mismatch (expected " + expected_sha256 + ", got " +
                             sha + ")");
  }
  TabulatedSpectrum table;
  table.sha256 = sha;
  std::istringstream stream(bytes);
  std::string line;
  while (std::getline(stream, line)) {
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos) continue;
    if (line[first] == '#') {
      const auto pos = line.find("schema:");
      if (pos != std::string::npos) {
        std::istringstream head(line.substr(pos + 7U));
        head >> table.schema;
      }
      continue;
    }
    std::istringstream row(line);
    double e = 0.0;
    double l = 0.0;
    row >> e >> l;
    if (row.fail() || !std::isfinite(e) || !std::isfinite(l) || !(e > 0.0) || l < 0.0) {
      throw std::runtime_error("malformed spectral table row in '" + path + "'");
    }
    if (!table.energy_ev.empty() && !(e > table.energy_ev.back())) {
      throw std::runtime_error("spectral table '" + path + "' energies are not strictly increasing");
    }
    table.energy_ev.push_back(e);
    table.l_nu.push_back(l);
  }
  if (table.schema != expected_schema) {
    throw std::runtime_error("spectral table '" + path + "' declares schema '" + table.schema + "', expected '" +
                             expected_schema + "'");
  }
  if (table.energy_ev.size() < 2U) throw std::runtime_error("spectral table '" + path + "' has fewer than 2 rows");
  return table;
}

double SpectralShape::evaluate(double e) const {
  switch (kind) {
    case SpectralShapeKind::kTable:
      if (table == nullptr) throw std::logic_error("SpectralShape: table shape without a table");
      return table->evaluate(e);
    case SpectralShapeKind::kBlackbody: {
      // B_nu ~ nu^3 / (exp(h nu / k T) - 1)
      const double x = e / (k_boltzmann_ev_per_k * blackbody_temperature_k);
      return x > 700.0 ? 0.0 : e * e * e / std::expm1(x);
    }
    case SpectralShapeKind::kPowerLaw:
      return std::pow(e, -power_law_index);
    case SpectralShapeKind::kFlat:
      return 1.0;
  }
  return 0.0;
}

std::string SpectralShape::describe() const {
  std::ostringstream s;
  s.precision(17);
  switch (kind) {
    case SpectralShapeKind::kTable: s << "table:" << (table != nullptr ? table->schema + ":" + table->sha256 : "none"); break;
    case SpectralShapeKind::kBlackbody: s << "blackbody_reference:T_K=" << blackbody_temperature_k; break;
    case SpectralShapeKind::kPowerLaw: s << "power_law_reference:index=" << power_law_index; break;
    case SpectralShapeKind::kFlat: s << "flat"; break;
  }
  return s.str();
}

std::vector<double> radiationProfileEdgesEv(const std::string& profile, double upper) {
  if (profile == "grey") return {};
  if (profile == "ionizing3") {
    if (!(upper > k_ionization_potential_heii_ev)) {
      throw std::invalid_argument("ionizing3 requires an explicit upper edge above 54.417765 eV");
    }
    return {k_ionization_potential_hi_ev, k_ionization_potential_hei_ev, k_ionization_potential_heii_ev, upper};
  }
  if (profile == "broad_uv_ir") return {0.01, 1.0, 5.8, 11.2, 13.6, 24.6, 54.4, 100.0};
  throw std::invalid_argument("unknown radiation spectral profile '" + profile + "'");
}

RadiationSpectralDefinition buildRadiationSpectralDefinition(const std::string& profile, double upper,
                                                            const SpectralShape& shape, std::uint32_t nq) {
  RadiationSpectralDefinition def;
  def.profile = profile;
  def.quadrature_points_per_group = nq;
  def.shape_description = shape.describe();
  def.edges_ev = radiationProfileEdgesEv(profile, upper);
  std::ostringstream canon;
  canon.precision(17);
  canon << "chui_rt_spectral_definition_v1\nprofile=" << profile << "\nshape=" << def.shape_description
        << "\nquadrature=gauss_legendre_log_nu:" << nq << "\nsigma=HI:hydrogenic_OF06,HeI:verner96,HeII:hydrogenic_OF06\n";
  if (def.edges_ev.empty()) {
    def.grey = true;
    RadiationGroup g;
    def.groups.push_back(g);
    canon << "edges=grey\n";
    def.sha256 = core::internal::sha256Hex(canon.str());
    return def;
  }
  std::vector<double> x;
  std::vector<double> w;
  gaussLegendre(nq, x, w);
  canon << "edges=";
  for (double e : def.edges_ev) canon << e << ",";
  canon << "\n";
  for (std::size_t gi = 0; gi + 1U < def.edges_ev.size(); ++gi) {
    RadiationGroup g;
    g.lower_ev = def.edges_ev[gi];
    g.upper_ev = def.edges_ev[gi + 1U];
    const double l0 = std::log(g.lower_ev);
    const double l1 = std::log(g.upper_ev);
    double energy = 0.0;
    double number = 0.0;
    std::array<double, k_radiation_absorber_count> se{};
    std::array<double, k_radiation_absorber_count> sn{};
    for (std::uint32_t q = 0; q < nq; ++q) {
      const double le = 0.5 * (l1 - l0) * x[q] + 0.5 * (l1 + l0);
      const double e = std::exp(le);
      const double jac = 0.5 * (l1 - l0) * w[q] * e;  // dE = E dlnE
      const double s = shape.evaluate(e);
      energy += jac * s;
      number += jac * s / e;
      for (std::size_t a = 0; a < k_radiation_absorber_count; ++a) {
        const double sigma = photoionizationCrossSectionCm2(static_cast<RadiationAbsorber>(a), e);
        se[a] += jac * s * sigma;
        sn[a] += jac * s * sigma / e;
      }
    }
    if (!(energy > 0.0) || !(number > 0.0)) {
      throw std::domain_error("radiation group [" + std::to_string(g.lower_ev) + ", " + std::to_string(g.upper_ev) +
                              "] eV has no spectral support under the declared shape");
    }
    g.shape_energy_integral = energy;
    g.mean_photon_energy_ev = energy / number;
    g.lower_edge_density_per_ev = shape.evaluate(g.lower_ev) / energy;
    for (std::size_t a = 0; a < k_radiation_absorber_count; ++a) {
      g.sigma_energy_cm2[a] = se[a] / energy;
      g.sigma_number_cm2[a] = sn[a] / number;
      g.absorbed_energy_per_ionization_ev[a] = sn[a] > 0.0 ? se[a] / sn[a] : 0.0;
    }
    def.groups.push_back(g);
  }
  def.sha256 = core::internal::sha256Hex(canon.str());
  return def;
}

std::vector<double> radiationGroupLuminosityFractions(const RadiationSpectralDefinition& def, const SpectralShape&) {
  std::vector<double> fractions(def.groups.size(), 0.0);
  if (def.grey) {
    fractions[0] = 1.0;
    return fractions;
  }
  double total = 0.0;
  for (const RadiationGroup& g : def.groups) total += g.shape_energy_integral;
  if (!(total > 0.0)) throw std::domain_error("radiation source shape has no support in the modeled groups");
  for (std::size_t g = 0; g < def.groups.size(); ++g) fractions[g] = def.groups[g].shape_energy_integral / total;
  return fractions;
}

}  // namespace cosmosim::radiation
