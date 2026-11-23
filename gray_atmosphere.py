"""Plane-parallel gray LTE atmosphere with Eddington two-stream closure."""

import argparse
import csv
import math
import sys

SIGMA = 5.670374419e-5  # erg cm^-2 s^-1 K^-4
MU = 1.0 / math.sqrt(3.0)


def propagate(intensity, source_start, source_end, optical_distance):
    """Integrate a ray exactly across a segment with a linear source."""
    absorbed = -math.expm1(-optical_distance)
    # The integral of the source slope has a removable singularity at zero.
    if optical_distance < 1e-4:
        slope_weight = optical_distance * (
            0.5 - optical_distance / 6.0 + optical_distance**2 / 24.0
        )
    else:
        slope_weight = 1.0 - absorbed / optical_distance
    return (
        intensity * math.exp(-optical_distance)
        + source_start * absorbed
        + (source_end - source_start) * slope_weight
    )


def solve(teff, logg, opacity, layers, tau_bottom):
    """Return rows in cgs units, with optical depth increasing inward."""
    gravity = 10.0**logg
    target_flux = SIGMA * teff**4
    eddington_flux = target_flux / (4.0 * math.pi)
    tau = [
        math.expm1(math.log1p(tau_bottom) * i / (layers - 1))
        for i in range(layers)
    ]

    # Radiative equilibrium gives B = J and dH/dtau = 0.
    # The moment equation gives dJ/dtau = 3H. No incident surface
    # intensity implies J(0) = H / MU in the two-stream quadrature.
    source = [eddington_flux * (1.0 / MU + 3.0 * depth) for depth in tau]

    # Check the moment solution with a formal ray sweep. Use the diffusion
    # boundary for the outgoing ray at the bottom and vacuum at the top.
    outgoing = [0.0] * layers
    outgoing[-1] = source[-1] + 3.0 * MU * eddington_flux
    for i in range(layers - 2, -1, -1):
        distance = (tau[i + 1] - tau[i]) / MU
        outgoing[i] = propagate(outgoing[i + 1], source[i + 1], source[i], distance)

    incoming = 0.0
    for i, depth in enumerate(tau):
        if i:
            distance = (depth - tau[i - 1]) / MU
            incoming = propagate(incoming, source[i - 1], source[i], distance)
        flux = 2.0 * math.pi * MU * (outgoing[i] - incoming)
        column_mass = depth / opacity
        yield (
            depth,
            column_mass,
            (math.pi * source[i] / SIGMA)**0.25,
            gravity * column_mass,
            flux,
            (flux - target_flux) / target_flux,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--teff", type=float, default=4300, help="Effective temperature in K")
    parser.add_argument("--logg", type=float, default=1.5, help="log10 surface gravity in cgs")
    parser.add_argument("--opacity", type=float, default=0.1, help="Gray opacity in cm²/g")
    parser.add_argument("--layers", type=int, default=81)
    parser.add_argument("--tau-bottom", type=float, default=100)
    args = parser.parse_args()
    values = (args.teff, args.opacity, args.tau_bottom)
    if not all(math.isfinite(x) and x > 0 for x in values) or args.layers < 2:
        parser.error("Use positive finite temperature, opacity and bottom depth, and at least two layers.")
    if not math.isfinite(args.logg):
        parser.error("Surface gravity must be finite in log units.")

    rows = list(solve(args.teff, args.logg, args.opacity, args.layers, args.tau_bottom))
    writer = csv.writer(sys.stdout)
    writer.writerow(("tau", "column_mass_g_cm2", "temperature_K", "pressure_dyn_cm2",
                     "flux_erg_cm2_s", "relative_flux_residual"))
    writer.writerows(rows)
    gravity = 10.0**args.logg
    hydrostatic_error = max(
        abs((right[3] - left[3]) / (right[1] - left[1]) / gravity - 1.0)
        for left, right in zip(rows, rows[1:])
    )
    print(f"Maximum relative flux residual: {max(abs(row[-1]) for row in rows):.3e}", file=sys.stderr)
    print(f"Maximum relative hydrostatic residual: {hydrostatic_error:.3e}", file=sys.stderr)


if __name__ == "__main__":
    main()
