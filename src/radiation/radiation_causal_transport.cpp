// CHUI-EXT-PHYS-REF-v1 R04/R05 causal rays and HG Monte Carlo (source only,
// NOT compiled or executed in pass 01).

#include "cosmosim/radiation/radiation_causal_transport.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

#include "cosmosim/radiation/radiation_angular.hpp"
#include "cosmosim/radiation/radiation_reference_transport.hpp"

namespace cosmosim::radiation {
namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

[[nodiscard]] double uniform01(std::uint64_t seed, std::uint64_t stream, std::uint64_t counter) {
  const std::uint64_t bits = radiationCounterRandom(seed, stream, counter);
  return (static_cast<double>(bits >> 11U) + 0.5) * (1.0 / 9007199254740992.0);
}

[[nodiscard]] double domainDiagonal(const lightcone::LeafGeometry& g) {
  const auto e = g.domainExtent();
  return std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
}

}  // namespace

CausalRayTransport::CausalRayTransport(const lightcone::LeafGeometry& geometry, double signal_speed,
                                       std::size_t max_packets, double min_energy_fraction,
                                       std::size_t max_segments_per_step)
    : m_geometry(geometry),
      m_c(signal_speed),
      m_max_packets(max_packets),
      m_min_fraction(min_energy_fraction),
      m_max_segments(max_segments_per_step) {
  if (!(m_c > 0.0) || m_max_packets == 0U || !(m_min_fraction >= 0.0) || !(m_min_fraction < 1.0) ||
      m_max_segments == 0U) {
    throw std::invalid_argument("invalid causal ray transport settings");
  }
  m_packets.reserve(m_max_packets);
}

void CausalRayTransport::emitIsotropic(const std::array<double, 3>& position, double energy, double time,
                                       std::size_t directions, std::uint64_t source_id) {
  if (!(energy >= 0.0) || !std::isfinite(energy)) throw std::invalid_argument("emission energy must be finite >= 0");
  if (time < m_time) throw std::invalid_argument("causal emission before the current transport time");
  if (m_packets.size() + directions > m_max_packets) {
    throw std::runtime_error("causal ray queue capacity exceeded (emission refused; raise max_packets)");
  }
  const auto dirs = fibonacciDirections(directions);
  const double share = energy / static_cast<double>(directions);
  for (std::size_t d = 0; d < dirs.size(); ++d) {
    CausalPacket packet;
    packet.id = source_id * 1000003ULL + d;
    packet.origin = position;
    packet.direction = dirs[d];
    packet.emission_time = time;
    packet.emitted_energy = share;
    packet.energy = share;
    m_packets.push_back(packet);
  }
  m_ledger.emitted += energy;
}

void CausalRayTransport::advance(double t_end, std::span<const double> chi, std::span<double> absorbed) {
  if (t_end < m_time) throw std::invalid_argument("causal transport cannot run backwards");
  if (chi.size() < m_geometry.cellCount() || absorbed.size() < m_geometry.cellCount()) {
    throw std::invalid_argument("per-cell spans smaller than the geometry cell count");
  }
  std::vector<CausalPacket> survivors;
  survivors.reserve(m_packets.size());
  for (CausalPacket& packet : m_packets) {
    const double t_start = std::max(m_time, packet.emission_time);
    const double target = m_c * (t_end - packet.emission_time);  // total path length at t_end
    const double length = target - packet.travelled;
    if (!(length > 0.0) || t_end <= t_start) {
      survivors.push_back(packet);
      continue;
    }
    std::array<double, 3> position{};
    for (std::size_t a = 0; a < 3U; ++a) position[a] = packet.origin[a] + packet.travelled * packet.direction[a];
    const auto segments = m_geometry.traceStraightRay(position, packet.direction, length, m_max_segments);
    double covered = 0.0;
    bool terminated = false;
    for (const auto& seg : segments) {
      const double l = seg.length();
      const double attenuation = std::exp(-chi[seg.cell_row] * l);
      const double deposit = packet.energy * (1.0 - attenuation);
      absorbed[seg.cell_row] += deposit;
      m_ledger.absorbed += deposit;
      packet.energy *= attenuation;
      covered += l;
      if (packet.energy < m_min_fraction * packet.emitted_energy) {
        m_ledger.unresolved += packet.energy;
        packet.energy = 0.0;
        terminated = true;
        break;
      }
    }
    if (terminated) continue;
    packet.travelled += covered;
    if (covered < length * (1.0 - 1.0e-12)) {
      // The open boundary stopped the ray before the requested length.
      m_ledger.escaped += packet.energy;
      continue;
    }
    survivors.push_back(packet);
  }
  m_packets = std::move(survivors);
  m_time = t_end;
}

double CausalRayTransport::inFlightEnergy() const noexcept {
  double sum = 0.0;
  for (const CausalPacket& p : m_packets) sum += p.energy;
  return sum;
}

std::vector<std::byte> CausalRayTransport::serialize() const {
  std::vector<std::byte> out;
  const auto put = [&](const void* data, std::size_t n) {
    const auto* b = static_cast<const std::byte*>(data);
    out.insert(out.end(), b, b + n);
  };
  const std::uint64_t count = m_packets.size();
  put(&m_time, sizeof(m_time));
  put(&m_ledger, sizeof(m_ledger));
  put(&count, sizeof(count));
  for (const CausalPacket& p : m_packets) put(&p, sizeof(p));
  return out;
}

void CausalRayTransport::restore(std::span<const std::byte> bytes) {
  std::size_t offset = 0;
  const auto get = [&](void* data, std::size_t n) {
    if (offset + n > bytes.size()) throw std::runtime_error("causal ray restart payload truncated");
    std::memcpy(data, bytes.data() + offset, n);
    offset += n;
  };
  double time = 0.0;
  CausalLedger ledger;
  std::uint64_t count = 0;
  get(&time, sizeof(time));
  get(&ledger, sizeof(ledger));
  get(&count, sizeof(count));
  if (count > m_max_packets) throw std::runtime_error("causal ray restart exceeds the packet capacity");
  std::vector<CausalPacket> packets(static_cast<std::size_t>(count));
  for (CausalPacket& p : packets) get(&p, sizeof(p));
  if (offset != bytes.size()) throw std::runtime_error("causal ray restart payload has trailing bytes");
  m_time = time;
  m_ledger = ledger;
  m_packets = std::move(packets);
}

double sampleHenyeyGreensteinCosine(double g, double xi) {
  if (!(g > -1.0) || !(g < 1.0)) throw std::invalid_argument("HG asymmetry must lie in (-1, 1)");
  if (std::abs(g) < 1.0e-8) return 2.0 * xi - 1.0;
  const double t = (1.0 - g * g) / (1.0 - g + 2.0 * g * xi);
  return std::clamp((1.0 + g * g - t * t) / (2.0 * g), -1.0, 1.0);
}

std::array<double, 3> scatterDirection(const std::array<double, 3>& d, double mu, double phi) {
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  // Orthonormal basis (u, v, d).
  std::array<double, 3> u{};
  if (std::abs(d[2]) < 0.9) {
    u = {-d[1], d[0], 0.0};
  } else {
    u = {0.0, -d[2], d[1]};
  }
  const double un = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  for (double& x : u) x /= un;
  const std::array<double, 3> v{d[1] * u[2] - d[2] * u[1], d[2] * u[0] - d[0] * u[2], d[0] * u[1] - d[1] * u[0]};
  std::array<double, 3> out{};
  for (std::size_t a = 0; a < 3U; ++a) {
    out[a] = mu * d[a] + sin_theta * (std::cos(phi) * u[a] + std::sin(phi) * v[a]);
  }
  const double n = std::sqrt(out[0] * out[0] + out[1] * out[1] + out[2] * out[2]);
  for (double& x : out) x /= n;
  return out;
}

HgMonteCarloResult runHenyeyGreensteinMonteCarlo(const lightcone::LeafGeometry& geometry,
                                                 std::span<const double> chi_a, std::span<const double> chi_s,
                                                 const std::array<double, 3>& source, double source_energy,
                                                 const HgMonteCarloSettings& settings) {
  if (chi_a.size() < geometry.cellCount() || chi_s.size() < geometry.cellCount()) {
    throw std::invalid_argument("opacity spans smaller than the geometry");
  }
  if (settings.packets == 0U || settings.batches == 0U || settings.packets % settings.batches != 0U) {
    throw std::invalid_argument("packets must be a positive multiple of batches");
  }
  HgMonteCarloResult result;
  result.absorbed_per_cell.assign(geometry.cellCount(), 0.0);
  const double packet_energy = source_energy / static_cast<double>(settings.packets);
  const double chunk = std::max(domainDiagonal(geometry), geometry.finestCellWidth());
  const std::uint64_t per_batch = settings.packets / settings.batches;
  std::vector<double> batch_escape(settings.batches, 0.0);
  for (std::uint64_t id = 0; id < settings.packets; ++id) {
    std::uint64_t counter = 0;
    const double mu0 = 2.0 * uniform01(settings.seed, id, counter++) - 1.0;
    const double phi0 = 2.0 * k_pi * uniform01(settings.seed, id, counter++);
    const double s0 = std::sqrt(std::max(0.0, 1.0 - mu0 * mu0));
    std::array<double, 3> dir{s0 * std::cos(phi0), s0 * std::sin(phi0), mu0};
    std::array<double, 3> pos = source;
    bool alive = true;
    for (std::uint32_t event = 0; event < settings.max_events && alive; ++event) {
      double tau_target = -std::log(uniform01(settings.seed, id, counter++));
      bool interacted = false;
      while (!interacted) {
        const auto segments = geometry.traceStraightRay(pos, dir, chunk, settings.max_segments_per_flight);
        double covered = 0.0;
        for (const auto& seg : segments) {
          const double chi_t = chi_a[seg.cell_row] + chi_s[seg.cell_row];
          const double dtau = chi_t * seg.length();
          if (dtau >= tau_target && chi_t > 0.0) {
            const double l = tau_target / chi_t;
            // pos is the flight start; the interaction point lies s_begin + l along dir.
            for (std::size_t a = 0; a < 3U; ++a) pos[a] += (seg.s_begin + l) * dir[a];
            pos = geometry.wrap(pos);
            const double albedo = chi_s[seg.cell_row] / chi_t;
            if (uniform01(settings.seed, id, counter++) < albedo) {
              const double mu = sampleHenyeyGreensteinCosine(settings.asymmetry_g, uniform01(settings.seed, id, counter++));
              const double phi = 2.0 * k_pi * uniform01(settings.seed, id, counter++);
              dir = scatterDirection(dir, mu, phi);
              ++result.scatterings;
            } else {
              result.absorbed_per_cell[seg.cell_row] += packet_energy;
              result.absorbed += packet_energy;
              alive = false;
            }
            interacted = true;
            break;
          }
          tau_target -= dtau;
          covered += seg.length();
        }
        if (interacted) break;
        if (covered < chunk * (1.0 - 1.0e-12)) {  // left through an open face
          result.escaped += packet_energy;
          batch_escape[id / per_batch] += 1.0;
          alive = false;
          break;
        }
        for (std::size_t a = 0; a < 3U; ++a) pos[a] += covered * dir[a];
        pos = geometry.wrap(pos);
      }
    }
    if (alive) result.truncated += packet_energy;  // reached max_events: reported, not dropped
  }
  result.emitted = source_energy;
  double mean = 0.0;
  for (double& b : batch_escape) {
    b /= static_cast<double>(per_batch);
    mean += b;
  }
  mean /= static_cast<double>(settings.batches);
  double var = 0.0;
  for (const double b : batch_escape) var += (b - mean) * (b - mean);
  const double nb = static_cast<double>(settings.batches);
  result.escaped_fraction_mean = mean;
  result.escaped_fraction_standard_error = settings.batches > 1U ? std::sqrt(var / (nb - 1.0) / nb) : 0.0;
  return result;
}

}  // namespace cosmosim::radiation
