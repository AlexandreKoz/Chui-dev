#pragma once

// CHUI-EXT-PHYS-REF-v1, I02: typed options for the extended-physics campaign
// (massive-neutrino perturbation methods, cosmic rays, observables and
// instruments). Keys live in the authoritative `.param.txt` sections
// [neutrinos], [cosmic_rays], [observables] and [instruments]; the
// relativistic *background* species live in [cosmology] (CosmologyConfig).
//
// Disabled defaults reproduce the predecessor behaviour exactly: the
// normalized dump emits a section only when its model is enabled (config hash
// and restart admission of existing runs unchanged) and any non-default value
// on a disabled model is rejected as an inert key. Admission
// (validateExtendedPhysicsConfig, src/core/config_extended_physics.cpp)
// distinguishes: integrated unqualified paths (admitted), reference kernels
// whose runtime integration edge is missing (refused with the exact edge), and
// research exclusions (no key exists).

#include <cstdint>
#include <string>
#include <vector>

namespace cosmosim::core {

// --------------------------------------------------------------- neutrinos
enum class NeutrinoMethodConfig : std::uint8_t {
  kNone = 0,             // no massive-neutrino perturbation treatment
  kBackgroundOnly = 1,   // N01: homogeneous FD background only (requires cosmology.radiation_background)
  kTransferRatio = 2,    // N03a: instantaneous T_nu/T_cb ratio from versioned Boltzmann tables
  kLinearResponse = 3,   // N03b: time-dependent Newtonian linear response to the evolving potential
  kParticles = 4,        // N04: ordinary collisionless particles (kernels only; admission refused)
  kDeltaF = 5,           // N04: weighted delta-f particles (kernels only; admission refused)
};
enum class BoltzmannTableFormatConfig : std::uint8_t {
  kChuiV1 = 0,        // chui_boltzmann_transfer_v1 (self-describing metadata header)
  kClassTk = 1,       // CLASS *_tk.dat (column names from the header line)
  kCambTransfer = 2,  // CAMB transfer output, 13-column layout of CAMB >= 1.0
};
enum class BoltzmannGaugeConfig : std::uint8_t { kSynchronous = 0, kNewtonian = 1 };

struct NeutrinoPerturbationConfig {
  NeutrinoMethodConfig method = NeutrinoMethodConfig::kNone;
  BoltzmannTableFormatConfig table_format = BoltzmannTableFormatConfig::kChuiV1;
  BoltzmannGaugeConfig table_gauge = BoltzmannGaugeConfig::kSynchronous;
  std::vector<std::string> table_paths;     // one file per redshift (ascending a)
  std::vector<std::string> table_sha256;    // 64-hex per file, same order
  std::vector<double> table_redshifts;      // declared redshift per file (verified vs header when present)
  std::uint32_t response_k_bins = 64U;      // logarithmic |k| bins of the response history
  std::uint32_t response_history_capacity = 4096U;  // bounded history (refuse, never drop)
  double response_max_mean_velocity_c = 0.0;  // validity: <|p|>/(m c) at the start must be below (required)
  double response_k_min_h_mpc = 0.0;        // response table domain (required, > 0)
  double response_k_max_h_mpc = 0.0;        // must cover the PM Nyquist wavenumber
  bool operator==(const NeutrinoPerturbationConfig&) const = default;
};

// ------------------------------------------------------------- cosmic rays
enum class CosmicRayModelConfig : std::uint8_t {
  kNone = 0,
  kGrey = 1,        // CR01/CR02 kernels exist; runtime stage missing (refused)
  kTwoMoment = 2,   // CR03 kernels exist; runtime stage missing (refused)
  kSpectral = 3,    // CR04-CR06 kernels exist; runtime stage missing (refused)
};

struct CosmicRayConfig {
  CosmicRayModelConfig model = CosmicRayModelConfig::kNone;
  bool operator==(const CosmicRayConfig&) const = default;
};

struct ExtendedPhysicsConfig {
  NeutrinoPerturbationConfig neutrinos;
  CosmicRayConfig cosmic_rays;
  bool operator==(const ExtendedPhysicsConfig&) const = default;
};

[[nodiscard]] std::string neutrinoMethodConfigToString(NeutrinoMethodConfig value);
[[nodiscard]] std::string boltzmannTableFormatConfigToString(BoltzmannTableFormatConfig value);
[[nodiscard]] std::string boltzmannGaugeConfigToString(BoltzmannGaugeConfig value);
[[nodiscard]] std::string cosmicRayModelConfigToString(CosmicRayModelConfig value);

}  // namespace cosmosim::core
