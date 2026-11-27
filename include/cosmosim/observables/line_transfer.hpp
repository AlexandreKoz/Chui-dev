#pragma once

// CHUI-EXT-PHYS-REF-v1 L01 / RA01: deterministic line transfer along ordered
// leaf-ray segments (Ly-a forest absorption, HI 21 cm full-optical-depth
// emission/absorption). Shared by snapshot sightlines and lightcone rays.
//
// Velocity coordinate. Pixels and segments live on one velocity axis u [m/s]:
// snapshot sightlines use the proper Hubble velocity u = a H(a) s plus the
// LOS peculiar velocity (periodic with period a H L_box through a periodic
// box); lightcone rays use the log-wavelength velocity u = c ln[(1+z)(1+v/c)]
// (first-order Doppler, declared) so the same kernel applies. Inside one
// segment the gas state is constant and u varies linearly (Hubble flow), so
// the line profile is AVERAGED over the segment's velocity interval by
// Gauss-Legendre sub-quadrature (3 nodes per b/4 sub-interval): narrow lines
// in coarse cells are never point-sampled.
//
// Ly-a absorption (L01): tau(u_p) = sum_seg n_HI l sigma_int lambda_0 <phi_v>,
//   phi_v(dv) = H(a, dv/b) / (sqrt(pi) b),  b = sqrt(2 k T / m_H),
//   a = A lambda_0 / (4 pi b)  (Voigt; kDopplerOnly sets a = 0, a NAMED
//   comparator that is wrong for damped systems). F = exp(-tau).
// HI 21 cm (RA01): absorption coefficient (Einstein relations, g1/g0 = 3,
//   n1 = 3/4 n_HI, h nu << k T_s):
//   alpha = 3 c^3 A10 h n_HI / (32 pi nu0^2 k T_s) phi_v,
//   ordered formal solution in the Rayleigh-Jeans brightness temperature
//   from the FAR end toward the observer: T_b <- T_b e^-dtau + T_s (1 - e^-dtau),
//   starting from the declared background (T_cmb(z) or 0). Thin and
//   saturated (T_s = T_K) shortcuts are separately named options; the
//   homogeneous-slab formula is never applied to a whole nonuniform sightline.
// Spin temperature (Field 1958; recalled lead):
//   T_s^-1 = (T_g^-1 + x_c T_K^-1 + x_a T_c^-1) / (1 + x_c + x_a),
//   x_c = T_* / (A10 T_g) sum_i kappa_10^i(T_K) n_i, T_* = h nu0 / k,
//   x_a = 4 P_a T_* / (27 A10 T_g) (recalled lead; P_a from an explicit
//   provider, e.g. the L03 scattering-rate estimator), T_c = T_K declared.
//   kappa_10 only from an external versioned table (chui_21cm_collision_rates_v1).

#include <cstddef>
#include <string>
#include <vector>

#include "cosmosim/spectral/atomic_data.hpp"

namespace cosmosim::observables {

enum class LineProfileModel : unsigned char { kVoigt = 0, kDopplerOnly = 1 };

struct LineSegment {
  double path_proper_m = 0.0;
  double velocity_begin_m_s = 0.0;
  double velocity_end_m_s = 0.0;
  double absorber_density_m3 = 0.0;  // n_HI
  double temperature_k = 0.0;        // kinetic temperature (line width)
  double spin_temperature_k = 0.0;   // 21 cm only
};

struct VelocityAxis {
  double u_min = 0.0;
  double u_max = 0.0;
  std::size_t pixels = 0;
  bool periodic = false;  // snapshot sightlines through a periodic box
  [[nodiscard]] double width() const noexcept { return (u_max - u_min) / static_cast<double>(pixels); }
  [[nodiscard]] double center(std::size_t p) const noexcept { return u_min + (static_cast<double>(p) + 0.5) * width(); }
};

// Line-profile average over [u0, u1] at pixel velocity u_p (per m/s).
[[nodiscard]] double segmentAveragedProfile(const spectral::AtomicLine& line, LineProfileModel model, double temperature_k,
                                            double u0, double u1, double u_p, const VelocityAxis& axis);

[[nodiscard]] std::vector<double> lymanAlphaOpticalDepth(const std::vector<LineSegment>& segments,
                                                         const VelocityAxis& axis, LineProfileModel model);

struct ForestStatistics {
  double mean_flux = 0.0;
  std::vector<double> flux_pdf;          // normalized histogram on [0, 1]
  std::vector<double> power_k_s_per_m;   // k [s/m] of the 1D power
  std::vector<double> power_m_per_s;     // P_F(k) [m/s] of delta_F = F/<F> - 1
  double equivalent_width_m = 0.0;       // observed-frame W_lambda for rest lambda (1 + z)
  double neutral_column_m2 = 0.0;
};

// Pixels must be a power of two for the 1D power spectrum (reference FFT).
[[nodiscard]] ForestStatistics forestStatistics(const std::vector<double>& tau, const VelocityAxis& axis,
                                                std::size_t pdf_bins, double redshift,
                                                const std::vector<LineSegment>& segments);

struct SpinTemperatureInputs {
  double kinetic_temperature_k = 0.0;
  double radiation_temperature_k = 0.0;  // T_gamma (CMB at z, from S01)
  double collisional_coupling = 0.0;     // x_c
  double lyman_alpha_coupling = 0.0;     // x_a
  double color_temperature_k = 0.0;      // T_c (declared = T_K unless supplied)
};
[[nodiscard]] double spinTemperature(const SpinTemperatureInputs& in);

// kappa_10 table: rows 'T_K kappa_HH_cm3_s kappa_eH_cm3_s kappa_pH_cm3_s',
// '# schema: chui_21cm_collision_rates_v1' and '# reference:' required.
struct CollisionalCouplingTable {
  std::string sha256;
  std::string reference;
  std::vector<double> temperature_k;
  std::vector<double> kappa_hh, kappa_eh, kappa_ph;  // cm^3 s^-1
  // x_c for densities in m^-3; throws outside the temperature range.
  [[nodiscard]] double coupling(double t_k, double t_gamma, double n_h, double n_e, double n_p) const;
};
[[nodiscard]] CollisionalCouplingTable parseCollisionalCouplingTable(const std::string& bytes, const std::string& origin);
// x_a from a Ly-a scattering rate per atom P_a [1/s].
[[nodiscard]] double lymanAlphaCoupling(double scattering_rate_per_atom_s, double t_gamma);

enum class HiTransferModel : unsigned char {
  kFullOpticalDepth = 0,  // ordered formal solution (reference)
  kOpticallyThin = 1,     // named shortcut: T_b = T_bg + sum dtau (T_s - T_bg)
};

// Brightness temperature per channel after ordered transfer. Segments are
// ordered from the observer outward (s increasing); transfer runs from the
// far end. background_k is T_b entering the far end.
[[nodiscard]] std::vector<double> hiBrightnessTemperature(const std::vector<LineSegment>& segments,
                                                          const VelocityAxis& axis, HiTransferModel model,
                                                          double background_k);

}  // namespace cosmosim::observables
