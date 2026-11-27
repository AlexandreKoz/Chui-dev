#pragma once

// CHUI-EXT-PHYS-REF-v1 I05: inline read-only observable pipeline.
//   .param.txt [observables] -> ObservablesConfig (validated) ->
//   AnalysisRuntime::executeDiagnostics (cadence observables.interval_steps)
//   -> ObservablePipeline::execute(state) -> lightcone::LeafGeometry over the
//   populated patches + GasScienceAccessor -> line/SZ/DM kernels -> FITS
//   products and a provenance summary in output.output_directory.
// Products (axis-aligned parallel sightlines, sightlines_per_side^2):
//   lyman_alpha_forest: flux image [sightline, pixel] + mean flux / PDF / P_F(k);
//   hi_cube: brightness-temperature cube [channel, y, x] (delta T_b observer
//            frame in cosmological runs, T_b above zero background otherwise);
//   sz_maps: Compton-y, tau_e, kSZ Delta T/T images;
//   dispersion_map: DM image (RM only when a magnetic sample is attached).
// The pipeline never mutates the state; any unavailable input refuses with
// the provider reason (no zero substitution). Single rank only (admission).

#include <cstdint>
#include <memory>
#include <string>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/observables/line_transfer.hpp"

namespace cosmosim::observables {

struct ObservablePipelineReport {
  std::uint64_t step = 0;
  std::uint64_t sightlines = 0;
  std::uint64_t segments = 0;
  std::vector<std::string> written_files;
};

class ObservablePipeline {
 public:
  explicit ObservablePipeline(const core::SimulationConfig& config);
  [[nodiscard]] bool due(std::uint64_t step) const noexcept;
  ObservablePipelineReport execute(const core::SimulationState& state, std::uint64_t step, double scale_factor);

 private:
  core::SimulationConfig m_config;
  core::UnitSystem m_units;
  core::LambdaCdmBackground m_background;
  std::unique_ptr<CollisionalCouplingTable> m_collision_table;
};

}  // namespace cosmosim::observables
