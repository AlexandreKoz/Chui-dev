#pragma once

// CHUI-RT-BH-AGN-REF-v1, shared contract B03: unit and frame conversions.
//
// Single conversion authority for the RT, BH and AGN source operators
// (master design section 5). All physical source laws are evaluated in PROPER
// quantities at the declared evaluation scale factor and converted back to the
// canonical solver lanes exactly once:
//
//   rho_phys   = rho_stored / a^3        (comoving storage), rho_stored otherwise
//   dx_phys    = a dx_stored             (comoving storage)
//   V_phys     = a^3 V_stored,  A_phys = a^2 A_stored
//   v_pec      = stored velocity lane    (already physical peculiar velocity)
//   dt_phys    = integral da/(a H(a))    (owned by SourceInterval / cosmology)
//   isolated   : a = 1, H = 0
//
// Masses, momenta (m v_pec) and energies (m v^2) are frame invariant in the
// repository's convention (specific internal energy is the peculiar thermal
// energy; Hubble cooling of e is owned by hydro/gravity and is never
// re-applied by source operators). Radiation moments use a separate a^4
// convention (radiation/radiation_cosmology.hpp) and must never be converted
// with the a^3 matter rule.
//
// This header contains no physics model; it is the explicit transform table
// required by B03 so that every module states its conversions in one place.

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/units.hpp"

namespace cosmosim::physics {

struct ProperFrameContext {
  double scale_factor = 1.0;
  double hubble_si = 0.0;          // H(a) [s^-1]; 0 for isolated frames
  bool comoving_storage = false;   // positions/densities stored comoving
  double length_si_per_code = 1.0;
  double mass_si_per_code = 1.0;
  double velocity_si_per_code = 1.0;

  [[nodiscard]] static ProperFrameContext isolated(const core::UnitSystem& units) {
    ProperFrameContext frame;
    frame.length_si_per_code = units.length_si_per_code;
    frame.mass_si_per_code = units.mass_si_per_code;
    frame.velocity_si_per_code = units.velocity_si_per_code;
    return frame;
  }
  [[nodiscard]] static ProperFrameContext cosmological(const core::UnitSystem& units, double scale_factor,
                                                       double hubble_si) {
    if (!(scale_factor > 0.0) || !std::isfinite(scale_factor) || !std::isfinite(hubble_si)) {
      throw std::invalid_argument("ProperFrameContext: cosmological frame requires finite a > 0 and finite H");
    }
    ProperFrameContext frame = isolated(units);
    frame.scale_factor = scale_factor;
    frame.hubble_si = hubble_si;
    frame.comoving_storage = true;
    return frame;
  }

  void validate() const {
    if (!(scale_factor > 0.0) || !std::isfinite(scale_factor) || !(length_si_per_code > 0.0) ||
        !(mass_si_per_code > 0.0) || !(velocity_si_per_code > 0.0) || !std::isfinite(hubble_si)) {
      throw std::invalid_argument("ProperFrameContext: invalid frame or unit system");
    }
    if (!comoving_storage && (scale_factor != 1.0 || hubble_si != 0.0)) {
      throw std::invalid_argument("ProperFrameContext: isolated frames require a = 1 and H = 0");
    }
  }

  [[nodiscard]] double lengthScale() const noexcept { return comoving_storage ? scale_factor : 1.0; }
  [[nodiscard]] double properLength(double stored_length_code) const noexcept {
    return stored_length_code * lengthScale();
  }
  [[nodiscard]] double storedLength(double proper_length_code) const noexcept {
    return proper_length_code / lengthScale();
  }
  [[nodiscard]] double properArea(double stored_area_code) const noexcept {
    const double s = lengthScale();
    return stored_area_code * s * s;
  }
  [[nodiscard]] double properVolume(double stored_volume_code) const noexcept {
    const double s = lengthScale();
    return stored_volume_code * s * s * s;
  }
  [[nodiscard]] double properDensity(double stored_density_code) const noexcept {
    const double s = lengthScale();
    return stored_density_code / (s * s * s);
  }
  [[nodiscard]] double storedDensity(double proper_density_code) const noexcept {
    const double s = lengthScale();
    return proper_density_code * s * s * s;
  }

  [[nodiscard]] double timeSiPerCode() const noexcept { return length_si_per_code / velocity_si_per_code; }
  [[nodiscard]] double energySiPerCode() const noexcept {
    return mass_si_per_code * velocity_si_per_code * velocity_si_per_code;
  }
  [[nodiscard]] double powerSiPerCode() const noexcept { return energySiPerCode() / timeSiPerCode(); }
  [[nodiscard]] double newtonGCode() const noexcept {
    const double t = timeSiPerCode();
    return core::constants::k_newton_g_si * mass_si_per_code * t * t /
        (length_si_per_code * length_si_per_code * length_si_per_code);
  }
  [[nodiscard]] double speedOfLightCode() const noexcept {
    return core::constants::k_speed_of_light_si / velocity_si_per_code;
  }
  [[nodiscard]] double hubbleCode() const noexcept { return hubble_si * timeSiPerCode(); }
  [[nodiscard]] double secondsToCodeTime(double seconds) const noexcept { return seconds / timeSiPerCode(); }
  [[nodiscard]] double codeTimeToSeconds(double t_code) const noexcept { return t_code * timeSiPerCode(); }
};

// Mass-energy and accretion constants used by the BH/AGN ledgers, in SI.
// CODATA 2018 values; no calibration content.
inline constexpr double k_proton_mass_si = 1.67262192369e-27;
inline constexpr double k_thomson_cross_section_si = 6.6524587321e-29;
// Solar mass and parsec are taken from core::constants (single authority).
inline constexpr double k_year_si = 3.15576e7;
inline constexpr double k_electron_volt_si = 1.602176634e-19;
inline constexpr double k_planck_constant_si = 6.62607015e-34;
inline constexpr double k_boltzmann_constant_si = 1.380649e-23;
inline constexpr double k_radiation_constant_si = 7.565733250280007e-16;  // 4 sigma_SB / c [J m^-3 K^-4]

}  // namespace cosmosim::physics

namespace cosmosim::physics {

// Minimum-image displacement in STORED (comoving for cosmological runs)
// coordinates; proper offsets are frame.properLength(delta).
struct PeriodicBox {
  std::array<double, 3> extent{0.0, 0.0, 0.0};
  bool periodic = false;
  [[nodiscard]] std::array<double, 3> delta(const std::array<double, 3>& from,
                                            const std::array<double, 3>& to) const noexcept {
    std::array<double, 3> d{to[0] - from[0], to[1] - from[1], to[2] - from[2]};
    if (periodic) {
      for (std::size_t a = 0; a < 3; ++a) {
        if (extent[a] > 0.0) d[a] -= extent[a] * std::round(d[a] / extent[a]);
      }
    }
    return d;
  }
  [[nodiscard]] std::array<double, 3> wrap(std::array<double, 3> p) const noexcept {
    if (periodic) {
      for (std::size_t a = 0; a < 3; ++a) {
        if (extent[a] > 0.0) {
          p[a] -= extent[a] * std::floor(p[a] / extent[a]);
        }
      }
    }
    return p;
  }
};

}  // namespace cosmosim::physics
