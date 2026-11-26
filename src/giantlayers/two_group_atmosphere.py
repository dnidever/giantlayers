"""Semi-infinite, plane-parallel LTE benchmark with two gray opacity groups."""

import argparse
import csv
import math
import sys

from .gray_atmosphere import MU, SIGMA


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


def solve_numerical(teff, logg, k1, k2, layers=81, tau_bottom=100.0):
    """Solve the two-group moment equations on a finite depth grid."""
    total_h = SIGMA * teff**4 / (4.0 * math.pi)
    gravity = 10.0**logg
    reference_opacity = 2.0 / (1.0 / k1 + 1.0 / k2)
    r1, r2 = k1 / reference_opacity, k2 / reference_opacity
    r_sum = r1 + r2
    decay_squared = 3.0 * r1 * r2
    tau = [
        math.expm1(math.log1p(tau_bottom) * i / (layers - 1))
        for i in range(layers)
    ]

    # Eliminate the LTE source with local radiative equilibrium, and
    # eliminate H_2 with total flux conservation. For D = J_1 - J_2,
    # the resulting scalar equation is D'' = 3 r_1 r_2 D.
    #
    # At the surface, J_i = H_i / MU gives D' = a D + c.
    # Integrate D'' over the surface half-cell to impose this condition.
    a = 1.5 * r_sum * MU
    c = 1.5 * (r1 - r2) * total_h
    step = tau[1] - tau[0]
    diagonal = 1.0 / step + a + 0.5 * step * decay_squared
    upper = [-1.0 / step / diagonal]
    rhs = [-c / diagonal]

    # Assemble and eliminate the interior tridiagonal rows. The last
    # depth has D = 0, the deep diffusion boundary for this benchmark.
    for i in range(1, layers - 1):
        left = tau[i] - tau[i - 1]
        right = tau[i + 1] - tau[i]
        lower = -1.0 / left
        diagonal = (
            1.0 / left + 1.0 / right + 0.5 * (left + right) * decay_squared
            - lower * upper[-1]
        )
        upper.append(-1.0 / right / diagonal)
        rhs.append(-lower * rhs[-1] / diagonal)

    difference = [0.0] * layers
    for i in range(layers - 2, -1, -1):
        difference[i] = rhs[i] - upper[i] * difference[i + 1]

    surface_h1 = 0.5 * (total_h + MU * difference[0])
    weighted_j_surface = (surface_h1 / r1 + (total_h - surface_h1) / r2) / MU
    for i, depth in enumerate(tau):
        d = difference[i]
        # (J_1/r_1 + J_2/r_2)' = 3 H_total, exactly.
        weighted_j = weighted_j_surface + 3.0 * total_h * depth
        j2 = (weighted_j - d / r1) / (1.0 / r1 + 1.0 / r2)
        j1 = j2 + d
        source = 2.0 * (r1 * j1 + r2 * j2) / r_sum

        if i == 0:
            derivative = a * d + c
        else:
            step = depth - tau[i - 1]
            # Convert the face gradient to the node using D'' at the node.
            derivative = (d - difference[i - 1]) / step + 0.5 * step * decay_squared * d
        h1 = (derivative / 3.0 + r2 * total_h) / r_sum
        column_mass = depth / reference_opacity
        yield (depth, column_mass, (math.pi * source / SIGMA)**0.25,
               gravity * column_mass, 4.0 * math.pi * h1,
               4.0 * math.pi * (total_h - h1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--teff", type=float, default=4300, help="Effective temperature in K")
    parser.add_argument("--logg", type=float, default=1.5, help="log10 surface gravity in cgs")
    parser.add_argument("--opacity1", type=float, default=0.1, help="Group 1 opacity in cm²/g")
    parser.add_argument("--opacity2", type=float, default=1.0, help="Group 2 opacity in cm²/g")
    parser.add_argument("--layers", type=int, default=81)
    parser.add_argument("--tau-bottom", type=float, default=100.0)
    parser.add_argument("--numerical", action="store_true",
                        help="Use the finite-depth numerical solver")
    args = parser.parse_args()
    positive = (args.teff, args.opacity1, args.opacity2, args.tau_bottom)
    if not all(math.isfinite(x) and x > 0 for x in positive) or args.layers < 2:
        parser.error("Use positive finite temperature, opacities and bottom depth, and at least two layers.")
    if not math.isfinite(args.logg):
        parser.error("Surface gravity must be finite in log units.")

    writer = csv.writer(sys.stdout)
    writer.writerow(("tau_reference", "column_mass_g_cm2", "temperature_K",
                     "pressure_dyn_cm2", "group1_flux_erg_cm2_s", "group2_flux_erg_cm2_s"))
    solver = solve_numerical if args.numerical else solve
    writer.writerows(solver(args.teff, args.logg, args.opacity1, args.opacity2,
                            args.layers, args.tau_bottom))


if __name__ == "__main__":
    main()
