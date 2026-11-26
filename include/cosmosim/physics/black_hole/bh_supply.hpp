#pragma once

// CHUI-RT-BH-AGN-REF-v1, BH companion section 4 (BH03): neighborhood
// sufficient statistics and gas-supply estimators.
//
// Supply rates are REQUESTED supply at the disk boundary, never accepted
// transfer. Hot (T >= T_split) and cold (T < T_split) gas form a sharp
// partition of unity, so no mass is offered to both estimators.
//
// Implemented closures (provenance in docs/black_hole_engine.md):
//  * Bondi-Hoyle-Lyttleton, unboosted alpha = 1 reference, optional named
//    constant boost alpha:  mdot = 4 pi alpha G^2 M^2 rho_h/(c_s,h^2+v_rel,h^2)^{3/2}.
//  * Angular-momentum-suppressed comparator, Rosas-Guevara et al. (2015,
//    MNRAS 454, 1038) as used by EAGLE (Schaye et al. 2015 eq. 9):
//    mdot = mdot_Bondi min(1, (c_s/V_phi)^3 / C_visc).
//  * Gravitational-torque supply, Hopkins & Quataert (2011) as used by SIMBA
//    (Dave et al. 2019 eqs. 1-2), with eps_T reinterpreted as a disk-boundary
//    supply normalization (no retention factor), M_d^gal the RESOLVED
//    kinematic disk (counter-rotating bulge decomposition, M_bulge = 2 M_counter,
//    applied separately to cold gas and stars) and R0 the physical kernel radius.
// Zero eligible disk or gas gives zero torque supply. Every input of a rate is
// recorded in BhSupplyEstimate together with adequacy flags.

#include <array>
#include <cstdint>

#include "cosmosim/physics/rt_bh_agn_frame.hpp"

namespace cosmosim::physics {

// Kernel radius candidates evaluated in one sufficient-statistic pass.
inline constexpr std::size_t k_bh_kernel_radius_bins = 16U;

// Compact Wendland C2 kernel shape W(q) = (1-q)^4 (1+4q), q = r/h in [0,1].
[[nodiscard]] double bhKernelShape(double q) noexcept;
// Fraction of a box cell [center +- half_width] (proper lengths, offsets
// relative to the kernel center) weighted by W: integral_cell W dV / V_cell,
// evaluated with 2x2x2 Gauss-Legendre quadrature (exact for the polynomial
// inside the support when the cell is fully inside; a quadrature estimate at
// the support boundary).
[[nodiscard]] double bhCellKernelWeight(const std::array<double, 3>& offset_proper,
                                        const std::array<double, 3>& half_width_proper,
                                        double kernel_radius_proper) noexcept;

// Mergeable (rank-reducible) sufficient statistics of one kernel support.
// Lengths/volumes proper code units; velocities peculiar; positions are
// minimum-image offsets from the BH (proper). Weighted sums use the cell
// kernel weight w in [0, 1].
struct BhPhaseSums {
  double weight_mass = 0.0;                 // sum w m
  double weight_volume = 0.0;               // sum w V
  double weight_mass_cs2 = 0.0;             // sum w m c_s^2
  std::array<double, 3> weight_momentum{};  // sum w m v
  std::array<double, 3> weight_angular_momentum{};  // sum w m (r x (v - v_bh))
  double weight_mass_radius = 0.0;          // sum w m |r|
  double mass = 0.0;                        // unweighted sum m (support mass)
  std::uint64_t cells = 0U;
  void merge(const BhPhaseSums& rhs);
};

struct BhKernelStatistics {
  BhPhaseSums hot{};
  BhPhaseSums cold{};
  // Stars inside the support (unweighted; decomposition input).
  double star_mass = 0.0;
  std::array<double, 3> star_mass_position{};   // sum m r
  std::array<double, 3> star_momentum{};        // sum m v
  std::array<double, 3> star_r_cross_p{};       // sum m r x v
  std::uint64_t stars = 0U;
  // Cold-gas decomposition sums (unweighted).
  std::array<double, 3> cold_mass_position{};
  std::array<double, 3> cold_momentum{};
  std::array<double, 3> cold_r_cross_p{};
  // Pass 2 (axis known): counter-rotating masses.
  double cold_counter_rotating_mass = 0.0;
  double star_counter_rotating_mass = 0.0;
  double min_cell_width_proper = 0.0;  // 0 = no cells
  void merge(const BhKernelStatistics& rhs);
};

struct BhRadiusBinSums {
  std::array<double, k_bh_kernel_radius_bins> enclosed_gas_mass{};
  std::array<std::uint64_t, k_bh_kernel_radius_bins> enclosed_cells{};
  void merge(const BhRadiusBinSums& rhs);
};

// Log-spaced candidate radii in [r_min, r_max] (proper).
[[nodiscard]] std::array<double, k_bh_kernel_radius_bins> bhKernelRadiusCandidates(double r_min, double r_max);
// Smallest candidate whose enclosed mass reaches the target (r_max if none).
[[nodiscard]] std::size_t bhSelectKernelRadius(const BhRadiusBinSums& sums, double target_mass_code);

struct BhSupplyInputs {
  double bh_mass_code = 0.0;
  std::array<double, 3> bh_velocity_peculiar{};
  double kernel_radius_proper = 0.0;
  double alpha = 1.0;                 // Bondi boost (1 = unboosted reference)
  double c_visc = 0.0;                // > 0 selects angular suppression
  bool use_bondi = true;
  bool use_torque = true;
  double torque_epsilon = 0.0;
  std::uint32_t kernel_min_cells = 1U;
  std::uint32_t torque_min_samples = 1U;
  double softening_proper = 0.0;
};

struct BhSupplyEstimate {
  double hot_rate_code = 0.0;         // requested hot supply [mass_code / time_code]
  double cold_rate_code = 0.0;        // requested cold supply
  double bondi_unsuppressed_rate_code = 0.0;
  double angular_suppression_factor = 1.0;
  double hot_density_proper = 0.0;
  double hot_sound_speed = 0.0;
  double hot_relative_speed = 0.0;
  double hot_rotation_speed = 0.0;
  double bondi_radius_proper = 0.0;
  double bondi_radius_over_dx = 0.0;
  double bondi_radius_over_softening = 0.0;
  double enclosed_baryon_mass_code = 0.0;
  double disk_mass_code = 0.0;        // resolved galactic disk M_d^gal (NOT the subgrid reservoir)
  double disk_fraction = 0.0;         // f_d = M_d / M_enc
  double disk_gas_fraction = 0.0;     // f_gas = M_d,gas / M_d
  bool kernel_adequate = false;
  bool torque_adequate = false;
  bool bondi_estimator_rejected = false;
  [[nodiscard]] double totalRateCode() const noexcept { return hot_rate_code + cold_rate_code; }
};

// Cold-disk axis: direction of the baryonic angular momentum about the BH
// relative to the support's bulk velocity (zero vector when undefined).
[[nodiscard]] std::array<double, 3> bhDiskAxis(const BhKernelStatistics& stats);

[[nodiscard]] double bondiHoyleRate(double alpha, double newton_g, double bh_mass, double density_proper,
                                    double sound_speed, double relative_speed);
[[nodiscard]] double rosasGuevaraSuppression(double sound_speed, double rotation_speed, double c_visc) noexcept;
// SIMBA/Hopkins-Quataert torque rate in Msun/yr from physical inputs.
[[nodiscard]] double hopkinsQuataertTorqueRateMsunPerYear(double epsilon_t, double disk_fraction,
                                                          double bh_mass_msun, double enclosed_mass_msun,
                                                          double radius_pc, double disk_mass_msun,
                                                          double disk_gas_fraction) noexcept;

[[nodiscard]] BhSupplyEstimate estimateBhSupply(const BhKernelStatistics& stats, const BhSupplyInputs& inputs,
                                                const ProperFrameContext& frame);

// P1-BH01 resolved supply: inward and net mass flux through a moving
// spherical control surface S_r of radius r0 estimated from the cells whose
// centres lie in the shell [r0 - dr/2, r0 + dr/2] (embedded volume quadrature,
// proper units): Mdot_in = sum_i m_i max(-(v_i - v_s) . n_i, 0) / dr,
// Mdot_net = -sum_i m_i (v_i - v_s) . n_i / dr; v_s includes the Hubble flow
// of the proper surface when cosmological. Measurement only: debiting uses
// the accepted withdrawal. P1: admission refuses bh_supply_model=resolved_flux.
struct BhResolvedFluxSample {
  double mass_code = 0.0;
  double radius_proper = 0.0;
  double radial_velocity_relative = 0.0;  // (v - v_surface) . n
};
struct BhResolvedFlux {
  double inward_rate_code = 0.0;
  double net_rate_code = 0.0;
  std::uint64_t samples = 0U;
};
[[nodiscard]] BhResolvedFlux resolvedControlSurfaceFlux(const BhResolvedFluxSample* samples, std::size_t count,
                                                        double shell_radius_proper, double shell_width_proper) noexcept;

}  // namespace cosmosim::physics
