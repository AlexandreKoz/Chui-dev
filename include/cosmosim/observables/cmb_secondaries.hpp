#pragma once

// CHUI-EXT-PHYS-REF-v1 C01-C03: CMB primary inputs and late-time secondaries.
// No early-Universe Einstein-Boltzmann solver: primary spectra come from
// external CLASS/CAMB tables; the CMB background (T_cmb(a), u_cmb) comes only
// from the S01 authority.
//
// C01 tables: `chui_cmb_cl_v1` text (metadata: source_code, source_version,
//   units = muK2, normalization = l(l+1)Cl/2pi or Cl, lensed = yes/no,
//   cosmology_hash; columns ell TT EE TE [BB] [PP]). Flat-sky Gaussian
//   realization of a field with C_l on an N x N map of side theta [rad]:
//   a_k = sqrt(C_l(k) A / 2) (g1 + i g2) with Hermitian symmetry,
//   A = (N dtheta)^2, T(x) = (1/A) sum_k a_k e^{i k.x} (reference FFT),
//   counter-based deterministic Gaussian draws. Flat-sky small-angle
//   approximation declared; no full-sky claims.
// C03 ray integrals over ordered segments with PROPER path lengths:
//   y = (sigma_T / (m_e c^2)) int n_e k T_e dl,
//   tau = sigma_T int n_e dl,
//   kSZ: Delta T / T = -(sigma_T / c) int n_e v_r dl with v_r > 0 for gas
//   receding from the observer (convention stated; observer frame).
//   Nonrelativistic tSZ spectrum (Kompaneets, derived):
//   Delta T / T = y (x coth(x/2) - 4), x = h nu / (k T_cmb).
//   Relativistic tSZ: FAIL_CLOSED (no verified kernel in pass 01).
//   ISW / Rees-Sciama: Delta T / T = (2 / c^2) int dPhi/dt dt along the
//   light path, dt = a dchi / c, with dPhi/dt from two potential snapshots
//   (linear in a between them: controlled temporal interpolation; refuses
//   without two snapshots bracketing the segment epoch).
// C02 lensing on planes (flat sky, comoving transverse coordinates x = chi theta):
//   Poisson on plane j: lap_perp psi_j = (3 H0^2 Omega_m / c^2) Sigma_j,
//   Sigma_j = int_plane delta / a dchi  (comoving surface density contrast),
//   solved with the reference 2-D FFT (periodic plane), alpha_j = grad psi_j.
//   Born: kappa = (3 H0^2 Omega_m / (2 c^2)) sum_j chi_j (chi_s - chi_j)/chi_s Sigma_j(chi_j theta).
//   Multi-plane (ray shooting, lens-lens coupling):
//     theta_{j+1} = theta_j - alpha_j(x_j),  x_{j+1} = x_j + (chi_{j+1} - chi_j) theta_{j+1},
//   bilinear periodic interpolation of alpha_j, beta = x_s / chi_s, and
//   kappa = 1 - tr(d beta / d theta) / 2 by central finite differences.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "cosmosim/core/cosmology.hpp"

namespace cosmosim::observables {

struct CmbSpectraTable {
  std::string sha256;
  std::string source_code;
  std::string source_version;
  bool lensed = false;
  std::vector<double> ell;
  std::vector<double> cl_tt, cl_ee, cl_te;  // C_l in muK^2 (normalization removed)
  [[nodiscard]] double clTT(double ell_value) const;  // linear interpolation, throws outside
};

[[nodiscard]] CmbSpectraTable parseCmbSpectraTable(const std::string& bytes, const std::string& origin);

struct FlatSkyMap {
  std::size_t n = 0;
  double side_rad = 0.0;
  std::vector<double> values;  // row-major [iy][ix]
};

[[nodiscard]] FlatSkyMap gaussianRealizationFlatSky(const CmbSpectraTable& table, std::size_t n, double side_rad,
                                                    std::uint64_t seed);

struct SzSegment {
  double path_proper_m = 0.0;
  double electron_density_m3 = 0.0;
  double electron_temperature_k = 0.0;
  double radial_velocity_m_s = 0.0;  // > 0 receding
};

struct SzIntegrals {
  double compton_y = 0.0;
  double optical_depth = 0.0;
  double ksz_delta_t_over_t = 0.0;
};

[[nodiscard]] SzIntegrals integrateSz(const std::vector<SzSegment>& segments);
// Nonrelativistic thermal SZ spectral factor x coth(x/2) - 4 at frequency nu for T_cmb.
[[nodiscard]] double tszSpectralFactor(double frequency_hz, double t_cmb_k);

struct IswSegment {
  double comoving_length_m = 0.0;
  double scale_factor = 1.0;     // epoch of the segment (emission time)
  double potential_a0 = 0.0;     // Phi at snapshot a0 at the segment position [m^2/s^2]
  double potential_a1 = 0.0;     // Phi at snapshot a1
};
// dPhi/dt from the linear-in-a interpolation between snapshots a0 < a1:
// dPhi/dt = (Phi1 - Phi0)/(a1 - a0) * a H(a). Refuses segments outside [a0, a1].
[[nodiscard]] double integrateIsw(const std::vector<IswSegment>& segments, double a0, double a1,
                                  const core::LambdaCdmBackground& background);

struct LensPlane {
  double chi_m = 0.0;             // comoving distance of the plane
  double scale_factor = 1.0;
  std::size_t n = 0;              // n x n periodic plane
  double side_m = 0.0;            // comoving side length
  std::vector<double> sigma_m;    // Sigma_j = int delta / a dchi  [m], row-major
};

struct LensingResult {
  std::size_t n = 0;              // rays per side
  double field_rad = 0.0;
  std::vector<double> kappa_born;
  std::vector<double> kappa_multiplane;
  std::vector<double> deflection_x_rad, deflection_y_rad;  // theta - beta
};

// Rays on an n_rays x n_rays grid over field_rad centred on the optical axis
// (plane coordinates are periodic and centred on the axis); source at chi_s.
[[nodiscard]] LensingResult rayShootLensPlanes(const std::vector<LensPlane>& planes, double chi_source_m,
                                               double hubble0_si, double omega_matter, std::size_t n_rays,
                                               double field_rad);

}  // namespace cosmosim::observables
