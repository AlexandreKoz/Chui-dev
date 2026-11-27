// CHUI-EXT-PHYS-REF-v1 CR04-CR07 spectral cosmic rays (source only, NOT
// compiled or executed in pass 01).

#include "cosmosim/cosmic_rays/cosmic_ray_spectrum.hpp"

#include <algorithm>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::cosmic_rays {
namespace {

namespace k = core::constants;
constexpr double k_four_pi = 4.0 * k::k_pi;

struct GaussRule {
  std::array<double, 24> x{};
  std::array<double, 24> w{};
};

[[nodiscard]] GaussRule buildRule() {
  GaussRule rule;
  constexpr int n = 24;
  for (int i = 0; i < n; ++i) {
    double x = std::cos(k::k_pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
    double dp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0;
      double p1 = x;
      for (int j = 2; j <= n; ++j) {
        const double p2 = ((2.0 * j - 1.0) * x * p1 - (j - 1.0) * p0) / static_cast<double>(j);
        p0 = p1;
        p1 = p2;
      }
      dp = static_cast<double>(n) * (x * p1 - p0) / (x * x - 1.0);
      const double step = p1 / dp;
      x -= step;
      if (std::abs(step) < 1.0e-16) break;
    }
    rule.x[static_cast<std::size_t>(i)] = x;
    rule.w[static_cast<std::size_t>(i)] = 2.0 / ((1.0 - x * x) * dp * dp);
  }
  return rule;
}

[[nodiscard]] const GaussRule& rule24() {
  static const GaussRule rule = buildRule();
  return rule;
}

// int_{lo}^{hi} g(p) dp in ln p with 24 nodes.
template <typename G>
[[nodiscard]] double integrateLog(double lo, double hi, const G& g) {
  const double l0 = std::log(lo);
  const double l1 = std::log(hi);
  const double half = 0.5 * (l1 - l0);
  const double mid = 0.5 * (l1 + l0);
  const GaussRule& r = rule24();
  double s = 0.0;
  for (std::size_t i = 0; i < r.x.size(); ++i) {
    const double p = std::exp(mid + half * r.x[i]);
    s += r.w[i] * g(p) * p;  // dp = p d ln p
  }
  return half * s;
}

[[nodiscard]] double kinetic(double p, double rest) { return (std::sqrt(1.0 + p * p) - 1.0) * rest; }

[[nodiscard]] double expm1OverX(double s) { return std::abs(s) < 1.0e-8 ? 1.0 + 0.5 * s : std::expm1(s) / s; }

}  // namespace

double crSpeciesRestEnergyJ(CrSpecies species) noexcept {
  const double c2 = k::k_speed_of_light_si * k::k_speed_of_light_si;
  return species == CrSpecies::kProton ? k::k_proton_mass_si * c2 : k::k_electron_mass_si * c2;
}

CrMomentumGrid CrMomentumGrid::logarithmic(double p_min, double p_max, std::size_t bins) {
  if (!(p_min > 0.0) || !(p_max > p_min) || bins == 0U) throw std::invalid_argument("invalid CR momentum grid");
  CrMomentumGrid g;
  g.edges.resize(bins + 1U);
  for (std::size_t i = 0; i <= bins; ++i) {
    g.edges[i] = p_min * std::pow(p_max / p_min, static_cast<double>(i) / static_cast<double>(bins));
  }
  return g;
}

void CrMomentumGrid::validate() const {
  if (edges.size() < 2U) throw std::invalid_argument("CR momentum grid needs >= 1 bin");
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (!(edges[i] > 0.0) || (i > 0 && !(edges[i] > edges[i - 1U]))) {
      throw std::invalid_argument("CR momentum edges must be > 0 and strictly increasing");
    }
  }
}

double powerLawNumber(double lo, double hi, double q, double f0) {
  const double ln_r = std::log(hi / lo);
  return k_four_pi * f0 * lo * lo * lo * ln_r * expm1OverX((3.0 - q) * ln_r);
}

double powerLawEnergy(double lo, double hi, double q, double f0, double rest) {
  return k_four_pi * f0 * integrateLog(lo, hi, [&](double p) { return p * p * std::pow(p / lo, -q) * kinetic(p, rest); });
}

double powerLawPressure(double lo, double hi, double q, double f0, double rest) {
  return k_four_pi / 3.0 * f0 *
         integrateLog(lo, hi, [&](double p) { return p * p * std::pow(p / lo, -q) * rest * p * p / std::sqrt(1.0 + p * p); });
}

CrPowerLawBin reconstructPowerLaw(double lo, double hi, double number, double energy, double rest, double q_min,
                                  double q_max) {
  if (number < 0.0 || energy < 0.0 || !std::isfinite(number) || !std::isfinite(energy)) {
    throw std::invalid_argument("CR bin moments must be finite and >= 0");
  }
  CrPowerLawBin bin;
  if (number == 0.0) return bin;  // empty bin
  const double target = energy / number;
  const auto mean_t = [&](double q) { return powerLawEnergy(lo, hi, q, 1.0, rest) / powerLawNumber(lo, hi, q, 1.0); };
  const double t_hard = mean_t(q_min);  // largest admissible <T>
  const double t_soft = mean_t(q_max);
  double q = 0.0;
  if (target >= t_hard) {
    q = q_min;
    bin.clamped = target > t_hard * (1.0 + 1.0e-12);
  } else if (target <= t_soft) {
    q = q_max;
    bin.clamped = target < t_soft * (1.0 - 1.0e-12);
  } else {
    double a = q_min;
    double b = q_max;
    for (int it = 0; it < 200 && b - a > 1.0e-13; ++it) {
      const double mid = 0.5 * (a + b);
      if (mean_t(mid) > target) {
        a = mid;
      } else {
        b = mid;
      }
    }
    q = 0.5 * (a + b);
  }
  bin.q = q;
  bin.f0 = number / powerLawNumber(lo, hi, q, 1.0);
  return bin;
}

CrSpectralState::CrSpectralState(CrMomentumGrid grid, CrSpectralSettings settings, std::size_t cells,
                                 core::MemoryGovernor* governor)
    : m_grid(std::move(grid)), m_settings(settings), m_cells(cells) {
  m_grid.validate();
  if (!(m_settings.q_max > m_settings.q_min)) throw std::invalid_argument("CR slope bounds must satisfy q_min < q_max");
  const std::uint64_t bytes = bytesFor(m_grid.bins(), cells, m_settings.representation);
  if (governor != nullptr) {
    m_reservation = governor->reserve(core::MemoryClass::kCanonicalPersistent, bytes, "cosmic_ray_spectral_state");
  }
  m_number.assign(m_grid.bins() * cells, 0.0);
  if (m_settings.representation == CrSpectralRepresentation::kPiecewisePowerLaw) {
    m_energy.assign(m_grid.bins() * cells, 0.0);
  }
  if (m_reservation.pending()) m_reservation.commit();
}

std::uint64_t CrSpectralState::bytesFor(std::size_t bins, std::size_t cells, CrSpectralRepresentation rep) {
  const std::uint64_t lanes = rep == CrSpectralRepresentation::kPiecewisePowerLaw ? 2U : 1U;
  return lanes * static_cast<std::uint64_t>(bins) * cells * sizeof(double);
}

std::span<double> CrSpectralState::number(std::size_t bin) { return {m_number.data() + bin * m_cells, m_cells}; }
std::span<double> CrSpectralState::energy(std::size_t bin) {
  if (m_energy.empty()) return {};
  return {m_energy.data() + bin * m_cells, m_cells};
}
std::span<const double> CrSpectralState::number(std::size_t bin) const { return {m_number.data() + bin * m_cells, m_cells}; }
std::span<const double> CrSpectralState::energy(std::size_t bin) const {
  if (m_energy.empty()) return {};
  return {m_energy.data() + bin * m_cells, m_cells};
}

CrPowerLawBin CrSpectralState::reconstruct(std::size_t bin, std::size_t cell) const {
  const double lo = m_grid.edges.at(bin);
  const double hi = m_grid.edges.at(bin + 1U);
  const double n = m_number[bin * m_cells + cell];
  const double rest = crSpeciesRestEnergyJ(m_settings.species);
  if (m_settings.representation == CrSpectralRepresentation::kFixedBins) {
    CrPowerLawBin b;
    b.q = 0.0;
    b.f0 = n > 0.0 ? n / powerLawNumber(lo, hi, 0.0, 1.0) : 0.0;
    return b;
  }
  return reconstructPowerLaw(lo, hi, n, m_energy[bin * m_cells + cell], rest, m_settings.q_min, m_settings.q_max);
}

double CrSpectralState::binEnergy(std::size_t bin, std::size_t cell) const {
  if (m_settings.representation == CrSpectralRepresentation::kPiecewisePowerLaw) return m_energy[bin * m_cells + cell];
  const CrPowerLawBin b = reconstruct(bin, cell);
  return powerLawEnergy(m_grid.edges[bin], m_grid.edges[bin + 1U], b.q, b.f0, crSpeciesRestEnergyJ(m_settings.species));
}

double CrSpectralState::totalEnergy(std::size_t cell) const {
  double sum = 0.0;
  for (std::size_t b = 0; b < m_grid.bins(); ++b) sum += binEnergy(b, cell);
  return sum;
}

double CrSpectralState::totalNumber(std::size_t cell) const {
  double sum = 0.0;
  for (std::size_t b = 0; b < m_grid.bins(); ++b) sum += m_number[b * m_cells + cell];
  return sum;
}

double CrSpectralState::pressure(std::size_t cell) const {
  double sum = 0.0;
  const double rest = crSpeciesRestEnergyJ(m_settings.species);
  for (std::size_t b = 0; b < m_grid.bins(); ++b) {
    const CrPowerLawBin bin = reconstruct(b, cell);
    if (bin.f0 > 0.0) sum += powerLawPressure(m_grid.edges[b], m_grid.edges[b + 1U], bin.q, bin.f0, rest);
  }
  return sum;
}

double CrSpectralState::distribution(std::size_t cell, double p) const {
  for (std::size_t b = 0; b < m_grid.bins(); ++b) {
    if (p >= m_grid.edges[b] && p <= m_grid.edges[b + 1U]) {
      const CrPowerLawBin bin = reconstruct(b, cell);
      return bin.f0 * std::pow(p / m_grid.edges[b], -bin.q);
    }
  }
  return 0.0;  // outside the grid: no particles are represented there
}

double momentumTimestep(const CrMomentumGrid& grid, const CrMomentumRate& rate, double cfl) {
  double dt = 1.0e300;
  for (std::size_t b = 0; b < grid.bins(); ++b) {
    const double dlnp = std::log(grid.edges[b + 1U] / grid.edges[b]);
    const double r = std::max(std::abs(rate(grid.edges[b])) / grid.edges[b],
                              std::abs(rate(grid.edges[b + 1U])) / grid.edges[b + 1U]);
    if (r > 0.0) dt = std::min(dt, cfl * dlnp / r);
  }
  return dt;
}

bool advanceMomentumSpace(CrSpectralState& state, std::size_t cell, const CrMomentumRate& rate, double dt,
                          CrSpectralLedger& ledger) {
  const CrMomentumGrid& grid = state.grid();
  const std::size_t nb = grid.bins();
  const bool power_law = state.settings().representation == CrSpectralRepresentation::kPiecewisePowerLaw;
  const double rest = crSpeciesRestEnergyJ(state.settings().species);
  std::vector<CrPowerLawBin> recon(nb);
  for (std::size_t b = 0; b < nb; ++b) {
    recon[b] = state.reconstruct(b, cell);
    if (recon[b].clamped) ++ledger.clamped_slopes;
  }
  const auto f_in_bin = [&](std::size_t b, double p) { return recon[b].f0 * std::pow(p / grid.edges[b], -recon[b].q); };
  std::vector<double> phi_n(nb + 1U, 0.0);
  std::vector<double> phi_e(nb + 1U, 0.0);
  for (std::size_t e = 0; e <= nb; ++e) {
    const double p = grid.edges[e];
    const double a = rate(p);
    double f_up = 0.0;
    if (a < 0.0 && e < nb) f_up = f_in_bin(e, p);           // flow downward from the bin above the edge
    if (a > 0.0 && e > 0U) f_up = f_in_bin(e - 1U, p);      // flow upward from the bin below the edge
    phi_n[e] = k_four_pi * p * p * f_up * a;
    phi_e[e] = kinetic(p, rest) * phi_n[e];
  }
  std::vector<double> new_n(nb);
  std::vector<double> new_e(nb, 0.0);
  double work_total = 0.0;
  for (std::size_t b = 0; b < nb; ++b) {
    new_n[b] = state.number(b)[cell] + dt * (phi_n[b] - phi_n[b + 1U]);
    if (power_law) {
      const double lo = grid.edges[b];
      const double hi = grid.edges[b + 1U];
      const double work = recon[b].f0 > 0.0
          ? k_four_pi * integrateLog(lo, hi, [&](double p) {
              return p * p * f_in_bin(b, p) * rest * p / std::sqrt(1.0 + p * p) * rate(p);
            })
          : 0.0;
      new_e[b] = state.energy(b)[cell] + dt * (phi_e[b] - phi_e[b + 1U] + work);
      work_total += dt * work;
      if (new_e[b] < 0.0) return false;
    }
    if (new_n[b] < 0.0) return false;
  }
  for (std::size_t b = 0; b < nb; ++b) {
    state.number(b)[cell] = new_n[b];
    if (power_law) state.energy(b)[cell] = new_e[b];
  }
  if (phi_n[0] < 0.0) {
    ledger.number_lost_low -= phi_n[0] * dt;
    ledger.energy_lost_low -= phi_e[0] * dt;
  }
  if (phi_n[nb] > 0.0) {
    ledger.number_lost_high += phi_n[nb] * dt;
    ledger.energy_lost_high += phi_e[nb] * dt;
  }
  ledger.work += work_total;
  return true;
}

}  // namespace cosmosim::cosmic_rays
