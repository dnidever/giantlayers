"""Solve LTE absorption transfer on a fixed gray atmosphere with ATLAS9 ODFs."""

import argparse
import csv
import math
import sys

from atlas9_odf import Atlas9ODF, C_NM_PER_S, LOG_P, quadrature
from gray_atmosphere import MU, SIGMA, propagate, solve as gray_atmosphere
import hminus
from marcs_model import read_model
from marcs_opacity import MarcsOpacity

PLANCK = 6.62607015e-27  # erg s
BOLTZMANN = 1.380649e-16  # erg/K
LIGHT = C_NM_PER_S * 1e-7  # cm/s


def planck_nu(frequency, temperature):
    exponent = PLANCK * frequency / (BOLTZMANN * temperature)
    if exponent > 700:
        return 0.0
    return 2 * PLANCK * frequency**3 / (LIGHT**2 * math.expm1(exponent))


def ray_flux(column_mass, source, opacity, thermal_bottom=False):
    """Return net two-stream flux density at each depth."""
    distances = [
        (right - left) * (a + b) / (2 * MU)
        for left, right, a, b in zip(column_mass, column_mass[1:], opacity, opacity[1:])
    ]
    outgoing = [0.0] * len(source)
    # At the base use B + mu * dB/dtau, estimated from the last cell.
    outgoing[-1] = source[-1]
    if not thermal_bottom or (sum(distances) >= 10 and distances[-1] > 0):
        outgoing[-1] += (source[-1] - source[-2]) / distances[-1]
    for i in range(len(source) - 2, -1, -1):
        outgoing[i] = propagate(outgoing[i + 1], source[i + 1], source[i], distances[i])

    incoming = 0.0
    flux = [2 * math.pi * MU * outgoing[0]]
    for i, distance in enumerate(distances):
        incoming = propagate(incoming, source[i], source[i + 1], distance)
        flux.append(2 * math.pi * MU * (outgoing[i + 1] - incoming))
    return flux


def scattering_ray_flux(column_mass, source, absorption, scattering):
    """Solve two-stream LTE transfer with isotropic coherent scattering.

    The symmetric intensity u obeys u'' = epsilon * (u - B) in ray
    optical depth, where epsilon is the true-absorption fraction.
    """
    extinction = [a + s for a, s in zip(absorption, scattering)]
    epsilon = [a / total for a, total in zip(absorption, extinction)]
    steps = [
        (right - left) * (a + b) / (2 * MU)
        for left, right, a, b in zip(column_mass, column_mass[1:], extinction, extinction[1:])
    ]
    bottom_outgoing = source[-1]
    if sum(steps) >= 10:
        bottom_outgoing += (source[-1] - source[-2]) / steps[-1]

    upper = []
    rhs = []
    # Integrate over each cell, including the surface and bottom half-cells.
    # Apply u' = u at the surface and u' = I_out - u at the bottom.
    last = len(source) - 1
    for i, thermal in enumerate(source):
        left = steps[i - 1] if i else 0.0
        right = steps[i] if i < last else 0.0
        left_inverse = 1 / left if i else 0.0
        right_inverse = 1 / right if i < last else 0.0
        coupling = 0.5 * (left + right) * epsilon[i]
        diagonal = left_inverse + right_inverse + coupling
        forcing = coupling * thermal
        if i == 0:
            diagonal += 1
        else:
            diagonal += left_inverse * upper[-1]
            forcing += left_inverse * rhs[-1]
        if i == last:
            diagonal += 1
            forcing += bottom_outgoing
        upper.append(-right_inverse / diagonal)
        rhs.append(forcing / diagonal)

    mean_intensity = [0.0] * len(source)
    mean_intensity[-1] = rhs[-1]
    for i in range(last - 1, -1, -1):
        mean_intensity[i] = rhs[i] - upper[i] * mean_intensity[i + 1]

    flux_factor = 4 * math.pi * MU
    flux = [flux_factor * mean_intensity[0]]
    for i in range(1, last):
        left_gradient = (mean_intensity[i] - mean_intensity[i - 1]) / steps[i - 1]
        curvature = epsilon[i] * (mean_intensity[i] - source[i])
        flux.append(flux_factor * (left_gradient + 0.5 * steps[i - 1] * curvature))
    flux.append(flux_factor * (bottom_outgoing - mean_intensity[-1]))
    return flux


def solve(odf, teff, logg, continuum, layers=161, tau_bottom=100.0, include_hminus=False, model=None, tabulated_continuum=None, continuum_scattering=False):
    """Return the band-integrated surface fluxes and total flux versus depth."""
    if continuum_scattering and tabulated_continuum is None:
        raise ValueError("Continuum scattering requires a MARCS continuum table")
    if tabulated_continuum is not None:
        if model is None or include_hminus:
            raise ValueError("Tabulated continuum requires a MARCS model and already includes H⁻")
        tabulated_continuum.check_profile(model)
    gas = None
    if model is None:
        profile = list(gray_atmosphere(teff, logg, continuum, layers, tau_bottom))
        mass = [row[1] for row in profile]
        temperature = [row[2] for row in profile]
        pressure = [row[3] for row in profile]
        if include_hminus:
            gas = [hminus.estimate_gas_state(t, p) for t, p in zip(temperature, pressure)]
    else:
        mass = [layer.column_mass for layer in model.layers]
        temperature = [layer.temperature for layer in model.layers]
        pressure = [layer.gas_pressure for layer in model.layers]
        layers = len(mass)
        if include_hminus:
            gas = [(layer.electron_density, layer.neutral_h_per_gram) for layer in model.layers]

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
        continuum_opacity = [continuum] * layers
        if tabulated_continuum is not None:
            wavelength = C_NM_PER_S / frequency
            continuum_opacity, scatter = zip(*tabulated_continuum.at_wavelength(wavelength))
        elif gas is not None:
            wavelength = C_NM_PER_S / frequency
            continuum_opacity = [
                continuum + sum(hminus.mass_opacity(wavelength, t, ne, hydrogen))
                for t, (ne, hydrogen) in zip(temperature, gas)
            ]
        if line_columns is None:
            ranked_opacities = (continuum_opacity,)
            weights = (sum(weights),)
        else:
            ranked_opacities = (
                [background + opacity for background, opacity in zip(continuum_opacity, rank)]
                for rank in zip(*line_columns[band])
            )
        surface_flux = 0.0
        for opacity, weight in zip(ranked_opacities, weights):
            if continuum_scattering:
                flux = scattering_ray_flux(mass, source, opacity, scatter)
            else:
                flux = ray_flux(mass, source, opacity,
                                thermal_bottom=tabulated_continuum is not None)
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
    parser.add_argument("--hminus", action="store_true",
                        help="Add LTE H⁻ absorption to the constant continuum")
    parser.add_argument("--marcs-model", metavar="FILE",
                        help="Use a fixed MARCS .mod or .mod.gz profile with plane-parallel transfer")
    parser.add_argument("--marcs-continuum", metavar="FILE",
                        help="Use true absorption from a matching MARCS .opa or .opa.gz")
    parser.add_argument("--continuum-scattering", action="store_true",
                        help="Include isotropic coherent scattering from the MARCS continuum")
    args = parser.parse_args()
    positive = (args.teff, args.continuum, args.tau_bottom)
    if (not all(math.isfinite(x) and x > 0 for x in positive)
            or not math.isfinite(args.logg) or args.layers < 2):
        parser.error("Use positive finite temperature, continuum and depth, "
                     "finite logg, and at least two layers.")
    try:
        odf = None if args.continuum_only else Atlas9ODF(args.odf)
        model = read_model(args.marcs_model) if args.marcs_model else None
        tabulated = MarcsOpacity(args.marcs_continuum) if args.marcs_continuum else None
        bands, flux = solve(odf, args.teff, args.logg, args.continuum,
                            args.layers, args.tau_bottom, include_hminus=args.hminus, model=model,
                            tabulated_continuum=tabulated,
                            continuum_scattering=args.continuum_scattering)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("wavelength_min_nm", "wavelength_max_nm", "surface_flux_erg_cm2_s"))
    writer.writerows(bands)
    reference = SIGMA * (model.teff if model else args.teff)**4
    print(f"Bolometric surface flux / (sigma Teff^4): {flux[0] / reference:.6f}",
          file=sys.stderr)
    print(f"Maximum relative flux residual through the fixed atmosphere: "
          f"{max(abs(value / reference - 1) for value in flux):.6f}", file=sys.stderr)


if __name__ == "__main__":
    main()
