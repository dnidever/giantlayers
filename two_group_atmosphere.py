"""Semi-infinite, plane-parallel LTE benchmark with two gray opacity groups."""

import argparse
import csv
import math
import sys

from gray_atmosphere import MU, SIGMA


def solve(teff, logg, k1, k2, layers=81, tau_bottom=100.0):
    """Sample the radiative-equilibrium solution on a reference-depth grid.

    Both groups have half the integrated Planck function at every temperature.
    The boundary conditions are zero incident surface intensity and diffusion
    at infinite depth. tau_bottom sets the output range, not a boundary.
    """
    target_flux = SIGMA * teff**4
    h_total = target_flux / (4.0 * math.pi)
    gravity = 10.0**logg
    reference_opacity = 2.0 / (1.0 / k1 + 1.0 / k2)
    r1 = k1 / reference_opacity
    r2 = k2 / reference_opacity

    # In reference optical depth, J_i' = 3 r_i H_i and
    # H_i' = r_i (J_i - B/2). The linear diffusion mode carries H_total.
    # The decaying mode has J_2 perturbation = -(r2/r1) times that of J_1.
    slope = 3.0 * h_total
    decay_rate = math.sqrt(3.0 * r1 * r2)
    ratio = r2 / r1
    surface_factor1 = 1.0 + decay_rate / (3.0 * MU * r1)
    surface_factor2 = 1.0 + decay_rate / (3.0 * MU * r2)
    diffusion1 = slope / (6.0 * r1)
    diffusion2 = slope / (6.0 * r2)
    amplitude = ((diffusion1 - diffusion2) / MU
                 / (surface_factor1 + ratio * surface_factor2))
    intercept = 2.0 * (diffusion1 / MU - amplitude * surface_factor1)

    for i in range(layers):
        tau = math.expm1(math.log1p(tau_bottom) * i / (layers - 1))
        perturbation = amplitude * math.exp(-decay_rate * tau)
        diffusion_source = (intercept + slope * tau) / 2.0
        j1 = diffusion_source + perturbation
        j2 = diffusion_source - ratio * perturbation
        source = 2.0 * (r1 * j1 + r2 * j2) / (r1 + r2)
        h1 = diffusion1 - decay_rate * perturbation / (3.0 * r1)
        h2 = diffusion2 + decay_rate * ratio * perturbation / (3.0 * r2)
        column_mass = tau / reference_opacity
        yield (tau, column_mass, (math.pi * source / SIGMA)**0.25,
               gravity * column_mass, 4.0 * math.pi * h1, 4.0 * math.pi * h2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--teff", type=float, default=4300, help="Effective temperature in K")
    parser.add_argument("--logg", type=float, default=1.5, help="log10 surface gravity in cgs")
    parser.add_argument("--opacity1", type=float, default=0.1, help="Group 1 opacity in cm²/g")
    parser.add_argument("--opacity2", type=float, default=1.0, help="Group 2 opacity in cm²/g")
    parser.add_argument("--layers", type=int, default=81)
    parser.add_argument("--tau-bottom", type=float, default=100.0)
    args = parser.parse_args()
    positive = (args.teff, args.opacity1, args.opacity2, args.tau_bottom)
    if not all(math.isfinite(x) and x > 0 for x in positive) or args.layers < 2:
        parser.error("Use positive finite temperature, opacities and bottom depth, and at least two layers.")
    if not math.isfinite(args.logg):
        parser.error("Surface gravity must be finite in log units.")

    writer = csv.writer(sys.stdout)
    writer.writerow(("tau_reference", "column_mass_g_cm2", "temperature_K",
                     "pressure_dyn_cm2", "group1_flux_erg_cm2_s", "group2_flux_erg_cm2_s"))
    writer.writerows(solve(args.teff, args.logg, args.opacity1, args.opacity2,
                           args.layers, args.tau_bottom))


if __name__ == "__main__":
    main()
