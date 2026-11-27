#pragma once

// CHUI-EXT-PHYS-REF-v1 S02/RA05/O01: observer-frame spectral axis shared by
// observables (Ly-a, HI, continuum/Stokes) and instruments. Channels are
// defined by explicit, strictly increasing EDGES in Hz (no implicit centers);
// a channel value is the channel average of the integrand computed with a
// fixed per-channel Gauss-Legendre rule (`subsamples` nodes), so bandwidth
// averaging (smearing / depolarization) is explicit and convergence-testable.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cosmosim::spectral {

class FrequencyAxis {
 public:
  FrequencyAxis() = default;
  // Throws unless edges are finite, > 0, strictly increasing and >= 2 values.
  FrequencyAxis(std::vector<double> edges_hz, std::uint32_t subsamples);
  // Uniform channels in frequency.
  [[nodiscard]] static FrequencyAxis uniform(double nu_min_hz, double nu_max_hz, std::size_t channels,
                                             std::uint32_t subsamples);

  [[nodiscard]] std::size_t channelCount() const noexcept { return m_edges.empty() ? 0U : m_edges.size() - 1U; }
  [[nodiscard]] const std::vector<double>& edgesHz() const noexcept { return m_edges; }
  [[nodiscard]] double centerHz(std::size_t channel) const;
  [[nodiscard]] double widthHz(std::size_t channel) const;
  [[nodiscard]] std::uint32_t subsamples() const noexcept { return m_subsamples; }
  // Gauss-Legendre node frequencies and weights (weights sum to 1) of a channel.
  [[nodiscard]] std::vector<double> nodeFrequenciesHz(std::size_t channel) const;
  [[nodiscard]] const std::vector<double>& nodeWeights() const noexcept { return m_unit_weights; }
  [[nodiscard]] std::string describe() const;

 private:
  std::vector<double> m_edges;
  std::uint32_t m_subsamples = 1U;
  std::vector<double> m_unit_nodes;    // on [0, 1]
  std::vector<double> m_unit_weights;  // sum to 1
};

}  // namespace cosmosim::spectral
