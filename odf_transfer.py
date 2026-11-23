"""Solve LTE absorption transfer on a fixed gray atmosphere with ATLAS9 ODFs."""

import argparse
import csv
import math
import sys

from atlas9_odf import Atlas9ODF, C_NM_PER_S, LOG_P, quadrature
from gray_atmosphere import MU, SIGMA, propagate, solve as gray_atmosphere

PLANCK = 6.62607015e-27  # erg s
BOLTZMANN = 1.380649e-16  # erg/K
LIGHT = C_NM_PER_S * 1e-7  # cm/s


def planck_nu(frequency, temperature):
    exponent = PLANCK * frequency / (BOLTZMANN * temperature)
    if exponent > 700:
        return 0.0
    return 2 * PLANCK * frequency**3 / (LIGHT**2 * math.expm1(exponent))


def ray_flux(column_mass, source, opacity):
    """Return net two-stream flux density at each depth."""
    distances = [
        (right - left) * (a + b) / (2 * MU)
        for left, right, a, b in zip(column_mass, column_mass[1:], opacity, opacity[1:])
    ]
    outgoing = [0.0] * len(source)
    # At the base use B + mu * dB/dtau, estimated from the last cell.
    outgoing[-1] = source[-1] + (source[-1] - source[-2]) / distances[-1]
    for i in range(len(source) - 2, -1, -1):
        outgoing[i] = propagate(outgoing[i + 1], source[i + 1], source[i], distances[i])

    incoming = 0.0
    flux = [2 * math.pi * MU * outgoing[0]]
    for i, distance in enumerate(distances):
        incoming = propagate(incoming, source[i], source[i + 1], distance)
        flux.append(2 * math.pi * MU * (outgoing[i + 1] - incoming))
    return flux


def solve(odf, teff, logg, continuum, layers=161, tau_bottom=100.0):
    """Return the band-integrated surface fluxes and total flux versus depth."""
    profile = list(gray_atmosphere(teff, logg, continuum, layers, tau_bottom))
    mass = [row[1] for row in profile]
    temperature = [row[2] for row in profile]
    pressure = [row[3] for row in profile]

    line_columns = None
    if odf is not None:
        pressure_floor = 10**LOG_P[0]
        line_columns = list(zip(*(
            odf.line_opacity(t, max(p, pressure_floor))
            for t, p in zip(temperature, pressure)
        )))

    bands = []
    total_flux = [0.0] * layers
    for band, (short, long, weights) in enumerate(quadrature()):
        frequency = 0.5 * C_NM_PER_S * (1 / short + 1 / long)
        source = [planck_nu(frequency, t) for t in temperature]
        if line_columns is None:
            ranked_opacities = ([continuum] * layers,)
            weights = (sum(weights),)
        else:
            ranked_opacities = (
                [continuum + opacity for opacity in rank]
                for rank in zip(*line_columns[band])
            )
        surface_flux = 0.0
        for opacity, weight in zip(ranked_opacities, weights):
            flux = ray_flux(mass, source, opacity)
            surface_flux += weight * flux[0]
            for i, value in enumerate(flux):
                total_flux[i] += weight * value
        bands.append((short, long, surface_flux))
    return bands, total_flux


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("odf", help="Path to an ATLAS9 BIG ODF binary")
    parser.add_argument("--teff", type=float, default=4300)
    parser.add_argument("--logg", type=float, default=1.5)
    parser.add_argument("--continuum", type=float, default=0.1,
                        help="Constant continuum absorption in cm²/g")
    parser.add_argument("--layers", type=int, default=161)
    parser.add_argument("--tau-bottom", type=float, default=100)
    parser.add_argument("--continuum-only", action="store_true",
                        help="Skip line opacity for comparison")
    args = parser.parse_args()
    positive = (args.teff, args.continuum, args.tau_bottom)
    if (not all(math.isfinite(x) and x > 0 for x in positive)
            or not math.isfinite(args.logg) or args.layers < 2):
        parser.error("Use positive finite temperature, continuum and depth, "
                     "finite logg, and at least two layers.")
    try:
        odf = None if args.continuum_only else Atlas9ODF(args.odf)
        bands, flux = solve(odf, args.teff, args.logg, args.continuum,
                            args.layers, args.tau_bottom)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("wavelength_min_nm", "wavelength_max_nm", "surface_flux_erg_cm2_s"))
    writer.writerows(bands)
    reference = SIGMA * args.teff**4
    print(f"Bolometric surface flux / (sigma Teff^4): {flux[0] / reference:.6f}",
          file=sys.stderr)
    print(f"Maximum relative flux residual through the fixed atmosphere: "
          f"{max(abs(value / reference - 1) for value in flux):.6f}", file=sys.stderr)


if __name__ == "__main__":
    main()
