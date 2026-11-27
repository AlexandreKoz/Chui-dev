// CHUI-EXT-PHYS-REF-v1 S02 frequency axis (source only, NOT compiled in pass 01).

#include "cosmosim/spectral/frequency_axis.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::spectral {
namespace {

void gaussLegendreUnit(std::uint32_t n, std::vector<double>& nodes, std::vector<double>& weights) {
  nodes.assign(n, 0.0);
  weights.assign(n, 0.0);
  if (n == 1U) {
    nodes[0] = 0.5;
    weights[0] = 1.0;
    return;
  }
  for (std::uint32_t i = 0; i < n; ++i) {
    double x = std::cos(core::constants::k_pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
    double dp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0;
      double p1 = x;
      for (std::uint32_t k = 2; k <= n; ++k) {
        const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / static_cast<double>(k);
        p0 = p1;
        p1 = p2;
      }
      dp = static_cast<double>(n) * (x * p1 - p0) / (x * x - 1.0);
      const double step = p1 / dp;
      x -= step;
      if (std::abs(step) < 1.0e-16) break;
    }
    // Map [-1, 1] -> [0, 1]; weights normalized to sum 1.
    nodes[i] = 0.5 * (x + 1.0);
    weights[i] = 1.0 / ((1.0 - x * x) * dp * dp);
  }
}

}  // namespace

FrequencyAxis::FrequencyAxis(std::vector<double> edges_hz, std::uint32_t subsamples)
    : m_edges(std::move(edges_hz)), m_subsamples(subsamples) {
  if (m_edges.size() < 2U) throw std::invalid_argument("FrequencyAxis needs at least two edges");
  if (subsamples < 1U || subsamples > 64U) throw std::invalid_argument("FrequencyAxis subsamples must lie in [1, 64]");
  for (std::size_t i = 0; i < m_edges.size(); ++i) {
    if (!std::isfinite(m_edges[i]) || !(m_edges[i] > 0.0)) throw std::invalid_argument("frequency edges must be > 0");
    if (i > 0 && !(m_edges[i] > m_edges[i - 1])) throw std::invalid_argument("frequency edges must increase strictly");
  }
  gaussLegendreUnit(subsamples, m_unit_nodes, m_unit_weights);
}

FrequencyAxis FrequencyAxis::uniform(double nu_min_hz, double nu_max_hz, std::size_t channels,
                                     std::uint32_t subsamples) {
  if (channels == 0U || !(nu_max_hz > nu_min_hz)) throw std::invalid_argument("invalid uniform frequency axis");
  std::vector<double> edges(channels + 1U);
  for (std::size_t i = 0; i <= channels; ++i) {
    edges[i] = nu_min_hz + (nu_max_hz - nu_min_hz) * static_cast<double>(i) / static_cast<double>(channels);
  }
  return FrequencyAxis(std::move(edges), subsamples);
}

double FrequencyAxis::centerHz(std::size_t channel) const {
  return 0.5 * (m_edges.at(channel) + m_edges.at(channel + 1U));
}

double FrequencyAxis::widthHz(std::size_t channel) const { return m_edges.at(channel + 1U) - m_edges.at(channel); }

std::vector<double> FrequencyAxis::nodeFrequenciesHz(std::size_t channel) const {
  const double lo = m_edges.at(channel);
  const double width = widthHz(channel);
  std::vector<double> out(m_unit_nodes.size());
  for (std::size_t i = 0; i < out.size(); ++i) out[i] = lo + width * m_unit_nodes[i];
  return out;
}

std::string FrequencyAxis::describe() const {
  std::ostringstream text;
  text.precision(17);
  text << "edges_hz=";
  for (std::size_t i = 0; i < m_edges.size(); ++i) text << (i ? "," : "") << m_edges[i];
  text << ";subsamples=" << m_subsamples;
  return text.str();
}

}  // namespace cosmosim::spectral
