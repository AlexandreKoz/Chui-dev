#pragma once

// CHUI-EXT-PHYS-REF-v1 S02 / N02 / N03: versioned Boltzmann transfer-table
// adapters. No native Einstein-Boltzmann engine exists in CHUI; CLASS or CAMB
// output is consumed read-only with explicit provenance, gauge, k units,
// redshift and coverage. Every query outside the table domain throws (no
// extrapolation in k or a).
//
// Canonical in-memory epoch (all formats are converted to it):
//   k_h_mpc         wavenumber [h/Mpc], strictly increasing
//   delta_cb        cold (CDM + baryon) density transfer, mass weighted with
//                   background weights Omega_cdm, Omega_b
//   delta_ncdm      total massive-neutrino density transfer, weighted by the
//                   per-eigenstate energy densities rho_i(a) of the SAME
//                   background authority (core::LambdaCdmBackground)
//   theta_cb, theta_ncdm (optional) velocity divergence with respect to
//                   CONFORMAL time in units of 1/Mpc (CLASS convention, c = 1).
// All transfer columns of one epoch share one (arbitrary) normalization, so
// only ratios and the epoch-to-epoch shape are used by CHUI.
//
// Formats:
//   chui_v1      self-describing: '# schema: chui_boltzmann_transfer_v1' plus
//                '# key: value' metadata (gauge, redshift, k_unit = h/Mpc,
//                source_code, source_version, ncdm_masses_ev, columns: ...)
//   class_tk     CLASS '*_tk.dat': column names from the 'N:name' header; uses
//                d_cdm, d_b, d_ncdm[i] (count must equal the configured
//                eigenstates) and t_cdm/t_b/t_ncdm[i] when present (Newtonian
//                gauge). In synchronous gauge t_cdm is absent by construction.
//   camb_transfer 13 columns of CAMB >= 1.0 transfer output
//                (k/h CDM baryon photon nu mass_nu total no_nu total_de Weyl
//                v_CDM v_b v_b-v_c); delta_cb = no_nu, delta_ncdm = mass_nu.
//                CAMB provides no massive-neutrino velocity, so theta_ncdm is
//                absent and linear-response initialization from CAMB is refused.

#include <cstddef>
#include <string>
#include <vector>

#include "cosmosim/core/config_extended_physics.hpp"
#include "cosmosim/core/cosmology.hpp"

namespace cosmosim::spectral {

struct BoltzmannTransferEpoch {
  std::string path;
  std::string sha256;
  std::string format;          // chui_v1 | class_tk | camb_transfer
  std::string source_code;     // CLASS | CAMB | declared in chui_v1
  std::string source_version;  // as recorded in the file, or "undeclared"
  core::BoltzmannGaugeConfig gauge = core::BoltzmannGaugeConfig::kSynchronous;
  double redshift = 0.0;
  bool redshift_from_header = false;
  std::vector<double> k_h_mpc;
  std::vector<double> delta_cb;
  std::vector<double> delta_ncdm;
  std::vector<double> theta_cb;    // empty when not provided
  std::vector<double> theta_ncdm;  // empty when not provided
  [[nodiscard]] bool hasVelocities() const noexcept { return !theta_cb.empty() && !theta_ncdm.empty(); }
  [[nodiscard]] double kMin() const { return k_h_mpc.front(); }
  [[nodiscard]] double kMax() const { return k_h_mpc.back(); }
};

struct BoltzmannLoadRequest {
  std::string path;
  std::string expected_sha256;
  core::BoltzmannTableFormatConfig format = core::BoltzmannTableFormatConfig::kChuiV1;
  core::BoltzmannGaugeConfig declared_gauge = core::BoltzmannGaugeConfig::kSynchronous;
  double declared_redshift = 0.0;
  double omega_baryon = 0.0;  // background weight for delta_cb (Omega_b)
};

// Loads one epoch; verifies SHA-256, schema/columns, monotone k, finite
// values, declared-vs-header gauge and redshift (|dz| <= 1e-6 (1+z)), and the
// massive-eigenstate count against `background`. Throws std::runtime_error
// with the exact defect.
[[nodiscard]] BoltzmannTransferEpoch loadBoltzmannTransferEpoch(const BoltzmannLoadRequest& request,
                                                                const core::LambdaCdmBackground& background);

// Parses already-read bytes (used by the loader and by tests; no hashing).
[[nodiscard]] BoltzmannTransferEpoch parseBoltzmannTransferEpoch(const std::string& bytes,
                                                                 const BoltzmannLoadRequest& request,
                                                                 const core::LambdaCdmBackground& background);

// Ordered set of epochs (ascending scale factor) with interpolation:
// in ln k (linear in the transferred ratio), and in ln a between bracketing
// epochs. All epochs must share gauge and cover [k_min, k_max].
class BoltzmannTransferSet {
 public:
  BoltzmannTransferSet() = default;
  explicit BoltzmannTransferSet(std::vector<BoltzmannTransferEpoch> epochs);

  [[nodiscard]] std::size_t epochCount() const noexcept { return m_epochs.size(); }
  [[nodiscard]] const BoltzmannTransferEpoch& epoch(std::size_t i) const { return m_epochs.at(i); }
  [[nodiscard]] double scaleFactorMin() const;
  [[nodiscard]] double scaleFactorMax() const;
  [[nodiscard]] double commonKMin() const noexcept { return m_k_min; }
  [[nodiscard]] double commonKMax() const noexcept { return m_k_max; }

  // delta_ncdm / delta_cb at (k, a); throws outside the k or a domain.
  [[nodiscard]] double ncdmToColdRatio(double k_h_mpc, double scale_factor) const;
  // theta_X / delta_cb in 1/Mpc (conformal-time CLASS units); requires velocities.
  [[nodiscard]] double thetaNcdmOverDeltaCold(double k_h_mpc, double scale_factor) const;
  [[nodiscard]] double thetaColdOverDeltaCold(double k_h_mpc, double scale_factor) const;
  // Canonical text digest of all epochs' provenance (paths excluded, hashes included).
  [[nodiscard]] std::string provenanceDigest() const;

 private:
  enum class Quantity : unsigned char { kNcdmRatio, kThetaNcdm, kThetaCold };
  [[nodiscard]] double epochValue(const BoltzmannTransferEpoch& epoch, Quantity quantity, double k_h_mpc) const;
  [[nodiscard]] double interpolate(Quantity quantity, double k_h_mpc, double scale_factor) const;
  std::vector<BoltzmannTransferEpoch> m_epochs;
  double m_k_min = 0.0;
  double m_k_max = 0.0;
};

}  // namespace cosmosim::spectral
