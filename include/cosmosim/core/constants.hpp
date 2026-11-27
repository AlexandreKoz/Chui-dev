#pragma once

namespace cosmosim::core::constants {

// Mathematical constant pi.
constexpr double k_pi = 3.141592653589793238462643383279502884;
// Newtonian gravitational constant in SI units [m^3 kg^-1 s^-2].
constexpr double k_newton_g_si = 6.67430e-11;                // m^3 kg^-1 s^-2
// Speed of light in vacuum [m s^-1].
constexpr double k_speed_of_light_si = 2.99792458e8;         // m s^-1
// Parsec in SI units [m].
constexpr double k_parsec_si = 3.0856775814913673e16;        // m
// Kiloparsec in SI units [m].
constexpr double k_kiloparsec_si = 1.0e3 * k_parsec_si;      // m
// Megaparsec in SI units [m].
constexpr double k_megaparsec_si = 1.0e6 * k_parsec_si;      // m
// Solar mass in SI units [kg].
constexpr double k_solar_mass_si = 1.98847e30;               // kg
// Kilometer in SI units [m].
constexpr double k_kilometer_si = 1.0e3;                     // m
// Centimeter in SI units [m].
constexpr double k_centimeter_si = 1.0e-2;                   // m
// Gram in SI units [kg].
constexpr double k_gram_si = 1.0e-3;                         // kg
// Seconds per megayear [s].
constexpr double k_seconds_per_megayear = 3.15576e13;        // s
// Hubble scale conversion 100 km/s/Mpc expressed in SI [s^-1].
constexpr double k_hubble_100_km_s_mpc_si =
    (100.0 * k_kilometer_si) / k_megaparsec_si;              // s^-1

// ---------------------------------------------------------------------------
// Fundamental physical constants shared by the extended-physics modules
// (CHUI-EXT-PHYS-REF-v1). Values marked "exact" are SI-defining constants of
// the 2019 SI redefinition; the others are CODATA 2018 recommended values.
// Modules must use these rather than private copies when they are added.
// Planck constant [J s] (exact).
constexpr double k_planck_si = 6.62607015e-34;
// Reduced Planck constant hbar = h / (2 pi) [J s].
constexpr double k_reduced_planck_si = k_planck_si / (2.0 * k_pi);
// Boltzmann constant [J K^-1] (exact).
constexpr double k_boltzmann_si = 1.380649e-23;
// Elementary charge [C] (exact); also one electron volt in joules.
constexpr double k_elementary_charge_si = 1.602176634e-19;
constexpr double k_electron_volt_si = k_elementary_charge_si;  // J
// Electron mass [kg] (CODATA 2018).
constexpr double k_electron_mass_si = 9.1093837015e-31;
// Proton mass [kg] (CODATA 2018).
constexpr double k_proton_mass_si = 1.67262192369e-27;
// Atomic mass constant [kg] (CODATA 2018).
constexpr double k_atomic_mass_unit_si = 1.66053906660e-27;
// Thomson cross section [m^2] (CODATA 2018).
constexpr double k_thomson_cross_section_si = 6.6524587321e-29;
// Vacuum electric permittivity [F m^-1] (CODATA 2018).
constexpr double k_vacuum_permittivity_si = 8.8541878128e-12;
// Riemann zeta values used by Fermi-Dirac / Bose-Einstein moment identities
// (mathematical constants, 16 significant digits).
constexpr double k_zeta_3 = 1.202056903159594285;
constexpr double k_zeta_5 = 1.036927755143369926;
// Radiation constant a_rad = pi^2 k^4 / (15 hbar^3 c^3) [J m^-3 K^-4], derived.
constexpr double k_radiation_constant_si =
    (k_pi * k_pi * k_boltzmann_si * k_boltzmann_si * k_boltzmann_si * k_boltzmann_si) /
    (15.0 * k_reduced_planck_si * k_reduced_planck_si * k_reduced_planck_si *
     k_speed_of_light_si * k_speed_of_light_si * k_speed_of_light_si);

}  // namespace cosmosim::core::constants
