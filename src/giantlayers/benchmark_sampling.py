"""Measure sampled opacity convergence on one fixed MARCS atmosphere."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

from .lte_chemistry import Thermochemistry
from .lte_opacity import read_lines
from .marcs_writer import read_document
from .molecular_cross_section import MolecularCrossSection
from .native_transfer import evaluate, read_wavelengths, reference_flux_budget


def nested_wavelengths(axis, stride, supplemental=()):
    """Retain common endpoints on nested grids and add fixed samples."""
    if stride < 1:
        raise ValueError("Strides must be positive")
    return tuple(sorted(set(axis[::stride]) | {axis[-1]} | set(supplemental)))


def compare_sampling(model, lines, composition, tables, strides, supplemental=(),
                     vmic=2.0, reference_wavelengths=None, quadrature="frequency"):
    """Yield flux and heating metrics for each grid on the same atmosphere."""
    if not tables:
        raise ValueError("Supply at least one molecular cross-section table")
    if reference_wavelengths is None:
        axis = tables[0].wavelengths
        if any(table.wavelengths != axis for table in tables[1:]):
            raise ValueError("Supply a reference grid or tables on a shared wavelength grid")
    else:
        axis = reference_wavelengths

    chemistry = Thermochemistry()
    for stride in strides:
        wavelengths = nested_wavelengths(axis, stride, supplemental)
        energy = evaluate(model, wavelengths, lines, composition, vmic,
                          chemistry, cross_sections=tables, quadrature=quadrature)
        budget = list(reference_flux_budget(model, energy.flux))
        yield (stride, len(wavelengths), energy.flux[0], budget[0][3],
               max(abs(row[-1]) for row in budget),
               abs(energy.closure) / budget[0][2])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", help="Reference MARCS .mod or .mod.gz")
    parser.add_argument("lines", help="Transition CSV, with only a header to use cross sections alone")
    parser.add_argument("composition", help="JSON abundance pattern")
    parser.add_argument("tables", nargs="+", help="Molecular cross-section tables")
    parser.add_argument("--strides", nargs="+", type=int, default=[4, 2, 1])
    parser.add_argument("--supplemental-wavelengths", type=Path, metavar="FILE")
    parser.add_argument("--vmic", type=float, default=2)
    parser.add_argument("--reference-wavelengths", metavar="FILE",
                        help="Use a plain or gzipped reference grid instead of the table grid")
    parser.add_argument("--quadrature", choices=("frequency", "wavelength"), default="frequency")
    args = parser.parse_args()
    if not math.isfinite(args.vmic) or args.vmic < 0:
        parser.error("Microturbulence must be finite and nonnegative")
    try:
        model = read_document(args.model)
        lines = read_lines(args.lines)
        with open(args.composition, encoding="utf-8") as stream:
            composition = json.load(stream)
        tables = [MolecularCrossSection(path) for path in args.tables]
        wavelengths = ()
        if args.supplemental_wavelengths:
            wavelengths = read_wavelengths(args.supplemental_wavelengths)
        reference = read_wavelengths(args.reference_wavelengths) if args.reference_wavelengths else None
        rows = compare_sampling(model, lines, composition, tables, args.strides,
                                wavelengths, args.vmic, reference, args.quadrature)
        writer = csv.writer(sys.stdout)
        writer.writerow(("stride", "wavelengths", "surface_radiative_flux_erg_cm2_s",
                         "relative_surface_total_flux_error", "max_relative_total_flux_error",
                         "relative_heating_closure_error"))
        for row in rows:
            writer.writerow(row)
            sys.stdout.flush()
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
