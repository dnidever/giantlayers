"""LTE H⁻ continuum absorption using the ATLAS9 cross-section tables."""

from bisect import bisect_right
import math

# HMINOP tabulates bound-free data from Mathisen (1984) and
# free-free data from Bell & Berrington (1987).
# https://wwwuser.oats.inaf.it/fiorella.castelli/sources/atlas9g/atlas9mem.for
BF_WAVELENGTH = (18.0, 19.6, 21.4, 23.6, 26.4, 29.8, 34.3, 40.4, 49.1, 62.6, 111.3, 112.1, 112.67,
 112.95, 113.05, 113.1, 113.2, 113.23, 113.5, 114.4, 121.0, 139.0, 164.0, 175.0, 200.0,
 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0,
 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0,
 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0, 1000.0, 1025.0, 1050.0, 1075.0,
 1100.0, 1125.0, 1150.0, 1175.0, 1200.0, 1225.0, 1250.0, 1275.0, 1300.0, 1325.0, 1350.0,
 1375.0, 1400.0, 1425.0, 1450.0, 1475.0, 1500.0, 1525.0, 1550.0, 1575.0, 1600.0, 1610.0,
 1620.0, 1630.0, 1643.91)
BF_VALUES = (0.067, 0.088, 0.117, 0.155, 0.206, 0.283, 0.414, 0.703, 1.24, 2.33, 11.6, 13.9, 24.3,
 66.7, 95.0, 56.6, 20.0, 14.6, 8.5, 7.1, 5.43, 5.91, 7.29, 7.918, 9.453, 11.08, 12.75,
 14.46, 16.19, 17.92, 19.65, 21.35, 23.02, 24.65, 26.24, 27.77, 29.23, 30.62, 31.94,
 33.17, 34.32, 35.37, 36.32, 37.17, 37.91, 38.54, 39.07, 39.48, 39.77, 39.95, 40.01,
 39.95, 39.77, 39.48, 39.06, 38.53, 37.89, 37.13, 36.25, 35.28, 34.19, 33.01, 31.72,
 30.34, 28.87, 27.33, 25.71, 24.02, 22.26, 20.46, 18.62, 16.74, 14.85, 12.95, 11.07,
 9.211, 7.407, 5.677, 4.052, 2.575, 1.302, 0.8697, 0.4974, 0.1989, 0.0)
FF_WAVENUMBERS = (0.5, 0.4, 0.35, 0.3, 0.25, 0.2, 0.18, 0.16, 0.14, 0.12, 0.1, 0.09, 0.08, 0.07, 0.06,
 0.05, 0.04, 0.03, 0.02, 0.01, 0.008, 0.006)
FF_THETA = (0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6)
FF_SHORT = (0.0178, 0.0222, 0.0308, 0.0402, 0.0498, 0.0596, 0.0695, 0.0795, 0.0896, 0.131, 0.172,
 0.0228, 0.028, 0.0388, 0.0499, 0.0614, 0.0732, 0.0851, 0.0972, 0.11, 0.16, 0.211,
 0.0277, 0.0342, 0.0476, 0.0615, 0.076, 0.0908, 0.105, 0.121, 0.136, 0.199, 0.262,
 0.0364, 0.0447, 0.0616, 0.0789, 0.0966, 0.114, 0.132, 0.15, 0.169, 0.243, 0.318, 0.052,
 0.0633, 0.0859, 0.108, 0.131, 0.154, 0.178, 0.201, 0.225, 0.321, 0.418, 0.0791, 0.0959,
 0.129, 0.161, 0.194, 0.227, 0.26, 0.293, 0.327, 0.463, 0.602, 0.0965, 0.117, 0.157,
 0.195, 0.234, 0.272, 0.311, 0.351, 0.39, 0.549, 0.711, 0.121, 0.146, 0.195, 0.241,
 0.288, 0.334, 0.381, 0.428, 0.475, 0.667, 0.861, 0.154, 0.188, 0.249, 0.309, 0.367,
 0.424, 0.482, 0.539, 0.597, 0.83, 1.07, 0.208, 0.25, 0.332, 0.409, 0.484, 0.557, 0.63,
 0.702, 0.774, 1.06, 1.36, 0.293, 0.354, 0.468, 0.576, 0.677, 0.777, 0.874, 0.969, 1.06,
 1.45, 1.83)
FF_LONG = (0.358, 0.432, 0.572, 0.702, 0.825, 0.943, 1.06, 1.17, 1.28, 1.73, 2.17, 0.448, 0.539,
 0.711, 0.871, 1.02, 1.16, 1.29, 1.43, 1.57, 2.09, 2.6, 0.579, 0.699, 0.924, 1.13, 1.33,
 1.51, 1.69, 1.86, 2.02, 2.67, 3.31, 0.781, 0.94, 1.24, 1.52, 1.78, 2.02, 2.26, 2.48,
 2.69, 3.52, 4.31, 1.11, 1.34, 1.77, 2.17, 2.53, 2.87, 3.2, 3.51, 3.8, 4.92, 5.97, 1.73,
 2.08, 2.74, 3.37, 3.9, 4.5, 5.01, 5.5, 5.95, 7.59, 9.06, 3.04, 3.65, 4.8, 5.86, 6.86,
 7.79, 8.67, 9.5, 10.3, 13.2, 15.6, 6.79, 8.16, 10.7, 13.1, 15.3, 17.4, 19.4, 21.2,
 23.0, 29.5, 35.0, 27.0, 32.4, 42.6, 51.9, 60.7, 68.9, 76.8, 84.2, 91.4, 117.0, 140.0,
 42.3, 50.6, 66.4, 80.8, 94.5, 107.0, 120.0, 131.0, 142.0, 183.0, 219.0, 75.1, 90.0,
 118.0, 144.0, 168.0, 191.0, 212.0, 234.0, 253.0, 325.0, 388.0)

K_B = 1.380649e-16  # erg/K
K_EV = 8.617333262e-5  # eV/K
M_H = 1.6735575e-24  # g
HC_OVER_K_NM = 1.438776877e7  # nm K
SAHA = 2.4148e15  # cm^-3 K^-3/2, before statistical weights
FF_LOG_WAVELENGTH = tuple(math.log(91.134 / value) for value in FF_WAVENUMBERS)
FF_LOG_VALUES = tuple(math.log(value * 1e-26) for value in FF_SHORT + FF_LONG)


def bracket(grid, value):
    """Return the interpolation bracket, clamping to the table endpoints."""
    value = max(grid[0], min(value, grid[-1]))
    i = min(bisect_right(grid, value) - 1, len(grid) - 2)
    return i, (value - grid[i]) / (grid[i + 1] - grid[i])


def bound_free_cross_section(wavelength_nm):
    """Photodetachment cross section in cm² per H⁻ ion."""
    if not BF_WAVELENGTH[0] <= wavelength_nm < BF_WAVELENGTH[-1]:
        return 0.0
    i, f = bracket(BF_WAVELENGTH, wavelength_nm)
    return 1e-18 * (BF_VALUES[i] + f * (BF_VALUES[i + 1] - BF_VALUES[i]))


def mass_opacity(wavelength_nm, temperature, electron_density, neutral_h_per_gram):
    """Return bound-free and free-free absorption in cm²/g."""
    pair = electron_density * neutral_h_per_gram
    ion_factor = math.exp(0.754209 / (K_EV * temperature)) / (
        4 * SAHA * temperature**1.5
    )
    stimulated = -math.expm1(-HC_OVER_K_NM / (wavelength_nm * temperature))
    bf = bound_free_cross_section(wavelength_nm) * stimulated * ion_factor * pair
    iw, fw = bracket(FF_LOG_WAVELENGTH, math.log(wavelength_nm))
    it, ft = bracket(FF_THETA, 5040 / temperature)
    per_pair = []
    for j in (it, it + 1):
        i = iw * len(FF_THETA) + j
        low, high = FF_LOG_VALUES[i], FF_LOG_VALUES[i + len(FF_THETA)]
        per_pair.append(math.exp(low + fw * (high - low)) * 5040 * K_B / FF_THETA[j])
    ff = pair * (per_pair[0] + ft * (per_pair[1] - per_pair[0]))
    return bf, ff


def estimate_gas_state(temperature, gas_pressure):
    """Provisional LTE hydrogen ionization and fixed trace electron donors."""
    hydrogen_fraction = 0.70
    helium_per_h = 0.28 / (4 * hydrogen_fraction)
    donors_per_h = 1e-4
    particles = gas_pressure / (K_B * temperature)
    saha = SAHA * temperature**1.5 * math.exp(-13.5984 / (K_EV * temperature))
    # For hydrogen ionization fraction x, n_e/n_H = x + donors_per_h.
    # The ideal-gas pressure includes both donor ions and their electrons.
    particles_per_h = 1 + helium_per_h + 2 * donors_per_h
    a = particles + saha
    b = particles * donors_per_h + saha * (particles_per_h - 1)
    c = saha * particles_per_h
    x = 2 * c / (b + math.sqrt(b * b + 4 * a * c))
    n_h = particles / (particles_per_h + x)
    return n_h * (x + donors_per_h), (1 - x) * hydrogen_fraction / M_H
