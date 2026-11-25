"""Compute a spherical LTE atmosphere from stellar parameters and opacity data."""

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
import sys

from benchmark_sampling import nested_wavelengths
from gray_atmosphere import SIGMA
from lte_chemistry import Thermochemistry
from lte_opacity import read_lines
from marcs_writer import MarcsDocument, write_document
from molecular_cross_section import MolecularCrossSection
from native_marcs import build_document
from native_transfer import read_wavelengths
from spherical_heating import balance_all_layers, evaluate, report_balance_step, write_balance_csv
from spherical_hydrostatic import HydrostaticBoundary, SOLAR_MASS


def initialize(teff, logg, mass, abundances, *, layers=56, surface_pressure=1.0,
               column_bottom=700.0, gray_tau_bottom=32.0, vmic=2.0, alpha=1.5,
               chemistry=None):
    """Build an analytic gray starting profile and independent mass grid.

    Mass is in solar masses, pressure and column mass are cgs, and
    abundances are absolute log epsilon values indexed by element.
    """
    chemistry = chemistry or Thermochemistry()
    if (layers < 4 or not all(math.isfinite(x) and x > 0 for x in
                             (teff, mass, surface_pressure, column_bottom))
            or not math.isfinite(logg) or not math.isfinite(gray_tau_bottom)
            or gray_tau_bottom <= 1
            or not all(math.isfinite(x) and x >= 0 for x in (vmic, alpha))):
        raise ValueError("Use positive stellar and grid parameters, at least four layers, and bottom gray depth > 1")
    if (not isinstance(abundances, dict) or abundances.get("H") != 12
            or not abundances.keys() <= chemistry.atoms.keys()
            or not all(isinstance(x, (float, int)) and math.isfinite(x) for x in abundances.values())):
        raise ValueError("Supply absolute log epsilon values for known elements, including H = 12")
    gravity = 10.0**logg
    if not math.isfinite(gravity) or gravity <= 0:
        raise ValueError("Surface gravity must be positive and finite")
    # Distribute depths in log(1 + m/m_scale) to cover the surface and
    # interior. A gray opacity increasing with m places the hot part of
    # the analytic starting profile deeper than a constant-opacity seed.
    scale = column_bottom/140
    bottom = math.log1p(column_bottom/scale)
    photosphere = math.log1p(column_bottom/scale * gray_tau_bottom**-0.25)
    anchor = max(1, min(layers - 2, round((layers - 1)*photosphere/bottom)))
    axis = [photosphere*(i/anchor-1) for i in range(anchor)]
    axis.extend((bottom-photosphere)*i/(layers-1-anchor) for i in range(layers-anchor))
    column = tuple(scale*math.expm1(x+photosphere) for x in axis)
    tau = [gray_tau_bottom*(m/column_bottom)**4 for m in column]
    temperature = [teff*(0.75*(x + 0.39))**0.25 for x in tau]
    if not chemistry.temperatures[0] <= temperature[0] <= temperature[-1] <= chemistry.temperatures[-1]:
        raise ValueError("Gray seed temperatures fall outside the LTE chemistry table")
    vector = [-99.0] * 92
    for element, value in abundances.items():
        vector[chemistry.atoms[element]["atomic_number"] - 1] = value
    # Reuse the transfer layer record shape. Hydrostatic pressure, density,
    # and depth are populated before any opacity or convection calculation.
    effective_abundances = {
        element: vector[atom["atomic_number"] - 1] for element, atom in chemistry.atoms.items()
    }
    boundary = HydrostaticBoundary(column, mass*SOLAR_MASS, gravity, surface_pressure,
                                   anchor, effective_abundances, chemistry)
    flux = SIGMA * teff**4
    header = (
        "GiantLayers gray seed",
        f"{teff:.12g} Teff [K]",
        f"{flux:.17g} Flux [erg/cm2/s]",
        f"{gravity:.17g} Surface gravity [cm/s2]",
        f"{vmic:.12g} Microturbulence [km/s]",
        f"{mass:.17g} Mass [Msun]",
        "",
        f"{boundary.reference_radius:.17g} Reference radius [cm]",
        "",
        f"{alpha:.12g} 8.00 0.076 0.00 Convection parameters",
        "",
    )
    structure = tuple((x, x, 0.0, t, 0.0, surface_pressure, 0.0, 0.0)
                      for x, t in zip(axis, temperature))
    thermo = tuple((x, column_bottom/tau[-1], 0.0, 0.0, 0.0, 0.0, m)
                   for x, m in zip(axis, column))
    return MarcsDocument(header, tuple(vector), structure, thermo, ()), boundary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wavelengths", type=Path, help="Vacuum wavelength grid in Å")
    parser.add_argument("abundances", type=Path, help="JSON element-to-log-epsilon mapping with H = 12")
    parser.add_argument("--teff", type=float, required=True)
    parser.add_argument("--logg", type=float, required=True)
    parser.add_argument("--mass", type=float, required=True, help="Mass in solar masses")
    parser.add_argument("--lines", type=Path, help="Atomic line-list CSV")
    parser.add_argument("--molecular-cross-section", action="append", default=[], type=Path)
    parser.add_argument("--layers", type=int, default=56)
    parser.add_argument("--surface-pressure", type=float, default=1.0, help="Gas pressure in dyn/cm²")
    parser.add_argument("--column-bottom", type=float, default=700.0, help="Bottom column mass in g/cm²")
    parser.add_argument("--gray-tau-bottom", type=float, default=32.0, help="Seed optical depth at bottom")
    parser.add_argument("--vmic", type=float, default=2.0, help="Microturbulence in km/s")
    parser.add_argument("--alpha", type=float, default=1.5, help="Mixing length in pressure scale heights")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--core-rays", type=int, default=24)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--upper-knots", type=int, default=12)
    parser.add_argument("--kernel", type=Path, help="Optional compiled spherical transfer library")
    parser.add_argument("--flux-tolerance", type=float, default=0.01,
                        help="Maximum fractional luminosity error for convergence")
    parser.add_argument("--output-model", type=Path, help="Write a MARCS-format .mod or .mod.gz record")
    parser.add_argument("--verify-full", action="store_true",
                        help="Check the final profile at all input wavelengths")
    args = parser.parse_args()
    try:
        abundances = json.loads(args.abundances.read_text(encoding="utf-8"))
        seed, boundary = initialize(
            args.teff, args.logg, args.mass, abundances, layers=args.layers,
            surface_pressure=args.surface_pressure, column_bottom=args.column_bottom,
            gray_tau_bottom=args.gray_tau_bottom, vmic=args.vmic, alpha=args.alpha,
        )
        reference = read_wavelengths(args.wavelengths)
        wavelengths = nested_wavelengths(reference, args.stride)
        lines = read_lines(args.lines) if args.lines else ()
        tables = [MolecularCrossSection(path) for path in args.molecular_cross_section]
        model, result, _, _, _ = balance_all_layers(
            seed, wavelengths, lines, {}, steps=args.steps, tables=tables,
            core_rays=args.core_rays, vmic=args.vmic, kernel_path=args.kernel,
            upper_knots=args.upper_knots, boundary=boundary, progress=report_balance_step,
            flux_tolerance=args.flux_tolerance,
        )
        if args.verify_full and wavelengths != reference:
            print(f"Verifying on {len(reference):,} wavelengths", file=sys.stderr, flush=True)
            result = evaluate(model, reference, lines, {}, tables, args.core_rays,
                              args.vmic, chemistry=boundary.chemistry,
                              kernel_path=args.kernel, native_convection=True)
        budget = write_balance_csv(sys.stdout, seed, model, result, native_convection=True)
        maximum = max(abs(row[-1]) for row in budget)
        status = "converged" if maximum <= args.flux_tolerance else "unconverged"
        print(f"Luminosity balance {status} at {args.flux_tolerance:.3%} tolerance "
              f"(maximum error {maximum:.3%})", file=sys.stderr)
        if args.output_model:
            document = build_document(model, result, lines, tables, boundary, args.vmic)
            header = (f"{document.header[0]} ({status}, max luminosity error {maximum:.3%})",
                      *document.header[1:])
            document = replace(document, header=header)
            write_document(document, args.output_model)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
