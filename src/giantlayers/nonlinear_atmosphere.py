"""Two-group LTE atmosphere with illustrative temperature-dependent opacity."""

import argparse
import csv
import math
import sys

from .gray_atmosphere import MU, SIGMA
from .opacity_table import OpacityTable


def radiation_field(depths, opacity1, opacity2):
    """Return B/(F/4π) for a frozen opacity profile and unit total H."""
    layers = len(depths)
    conductance = []
    diffusion_share = []
    for i in range(layers - 1):
        a = (opacity1[i] + opacity1[i + 1]) / 2
        b = (opacity2[i] + opacity2[i + 1]) / 2
        width = depths[i + 1] - depths[i]
        conductance.append(1 / (3 * (a + b) * width))
        diffusion_share.append(b / (a + b))

    # H_1' = r_1 r_2 (J_1 - J_2)/(r_1 + r_2) in reference optical depth.
    # Use a conservative half-cell at the surface and set J_1 = J_2
    # at the inner boundary, where the groups approach local diffusion.
    upper = []
    rhs = []
    for i in range(layers - 1):
        right = conductance[i]
        a, b = opacity1[i], opacity2[i]
        width = (depths[i + 1] - depths[max(0, i - 1)]) / 2
        absorption = width * a * b / (a + b)
        if i == 0:
            diagonal = right + MU / 2 + absorption
            forcing = diffusion_share[i] - 0.5
        else:
            left = conductance[i - 1]
            diagonal = left + right + absorption + left * upper[-1]
            forcing = diffusion_share[i] - diffusion_share[i - 1] + left * rhs[-1]
        upper.append(-right / diagonal)
        rhs.append(forcing / diagonal)

    difference = [0.0] * layers
    for i in range(layers - 2, -1, -1):
        difference[i] = rhs[i] - upper[i] * difference[i + 1]

    h1_surface = (1 + MU * difference[0]) / 2
    j_sum = 1 / MU
    source = []
    for i, d in enumerate(difference):
        a, b = opacity1[i], opacity2[i]
        if i:
            width = depths[i] - depths[i - 1]
            h1 = diffusion_share[i - 1] + conductance[i - 1] * (d - difference[i - 1])
            a_face = (a + opacity1[i - 1]) / 2
            b_face = (b + opacity2[i - 1]) / 2
            j_sum += 3 * width * (a_face * h1 + b_face * (1 - h1))
        source.append(j_sum + (a - b) / (a + b) * d)
    return source, h1_surface


def solve(teff, logg, layers=321, tau_bottom=100.0, initial_scale=1.0,
          strength=1.0, tolerance=1e-9, max_iterations=200, opacity_table=None):
    """Iterate the temperature, returning the profile and convergence history."""
    k1, k2 = 0.1, 1.0  # cm²/g in the constant-opacity limit
    reference = 2 / (1 / k1 + 1 / k2)
    depths = [
        math.expm1(math.log1p(tau_bottom) * i / (layers - 1))
        for i in range(layers)
    ]
    pressure = [10**logg * depth / reference for depth in depths]
    temperature = [
        initial_scale * teff * (0.75 * (depth + MU))**0.25 for depth in depths
    ]
    history = []
    for _ in range(max_iterations):
        if opacity_table is None:
            # Bounded opacity response for numerical experiments, with P in cgs.
            group_opacity = [
                (k1 * math.exp(strength * (
                    0.5 * math.tanh(math.log(t / teff)) + 0.25 * p / (p + 100)
                )), k2)
                for t, p in zip(temperature, pressure)
            ]
        else:
            group_opacity = [opacity_table(t, p) for t, p in zip(temperature, pressure)]
        opacity1 = [a / reference for a, _ in group_opacity]
        opacity2 = [b / reference for _, b in group_opacity]
        source, surface_share = radiation_field(depths, opacity1, opacity2)
        target = [teff * (b / 4)**0.25 for b in source]
        residual = max(abs(new / old - 1) for old, new in zip(temperature, target))
        history.append(residual)
        if residual <= tolerance:
            rows = [
                (depth, p, t, a, b)
                for depth, p, t, (a, b) in zip(depths, pressure, temperature, group_opacity)
            ]
            return rows, history, surface_share
        temperature = [(old + new) / 2 for old, new in zip(temperature, target)]

    raise RuntimeError(f"Temperature iteration failed. Final residual {history[-1]:.3e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--teff", type=float, default=4300)
    parser.add_argument("--logg", type=float, default=1.5)
    parser.add_argument("--layers", type=int, default=321)
    parser.add_argument("--tau-bottom", type=float, default=100)
    parser.add_argument("--initial-scale", type=float, default=1)
    parser.add_argument("--strength", type=float, default=1)
    parser.add_argument("--opacity-table", metavar="CSV",
                        help="Read both group opacities from a CSV table")
    args = parser.parse_args()
    positive = (args.teff, args.tau_bottom, args.initial_scale)
    if (not all(math.isfinite(x) and x > 0 for x in positive)
            or not math.isfinite(args.logg)
            or not 0 <= args.strength <= 2 or args.layers < 2):
        parser.error("Use finite inputs, positive temperature, depth and scale, "
                     "strength between 0 and 2, and at least two layers.")
    try:
        if args.opacity_table is not None:
            args.opacity_table = OpacityTable(args.opacity_table)
        rows, history, surface_share = solve(**vars(args))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("tau_reference", "pressure_dyn_cm2", "temperature_K",
                     "group1_opacity_cm2_g", "group2_opacity_cm2_g"))
    writer.writerows(rows)
    print(f"Iterations: {len(history)}, max relative temperature residual: {history[-1]:.3e}",
          file=sys.stderr)
    print(f"Group 1 surface flux: {surface_share * SIGMA * args.teff**4:.6e} erg cm⁻² s⁻¹",
          file=sys.stderr)


if __name__ == "__main__":
    main()
