"""Compute a spherical LTE atmosphere from stellar parameters and opacity data."""

import argparse
from bisect import bisect_right
import csv
from dataclasses import replace
import json
import math
from pathlib import Path
import sys
from tempfile import NamedTemporaryFile

from atlas_writer import write_atlas
from benchmark_sampling import nested_wavelengths
from collision_induced import HitranCIA
from gray_atmosphere import SIGMA
from lte_chemistry import Thermochemistry
from lte_opacity import read_lines
from marcs_writer import MarcsDocument, write_document
from molecular_cross_section import MolecularCrossSection
from native_marcs import build_document, photospheric_anchor, rosseland_depth
from native_transfer import read_wavelengths
from spherical_heating import (balance_all_layers, calibrate_heating, evaluate, evaluate_parallel,
                               report_balance_step, surface_equilibrium_temperature, total_flux_budget,
                               write_balance_csv)
from spherical_hydrostatic import HydrostaticBoundary, SOLAR_MASS, reconstruct_for_transfer

from transfer_cache import input_fingerprint, read_residuals, transfer_key, write_residuals


def initialize(teff, logg, mass, abundances, *, layers=56, surface_pressure=1.0,
               column_bottom=700.0, gray_tau_bottom=32.0, vmic=2.0, alpha=1.5,
               chemistry=None, deep_layers=None):
    """Build an analytic gray starting profile and independent mass grid.

    Mass is in solar masses, pressure and column mass are cgs, and
    abundances are absolute log epsilon values indexed by element.
    Deep layers are counted below the gray seed's photospheric anchor.
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
    try:
        gravity = 10.0**logg
    except OverflowError as exc:
        raise ValueError("Surface gravity must be positive and finite") from exc
    if not math.isfinite(gravity) or gravity <= 0:
        raise ValueError("Surface gravity must be positive and finite")
    # Distribute depths in log(1 + m/m_scale) to cover the surface and
    # interior. A gray opacity increasing with m places the hot part of
    # the analytic starting profile deeper than a constant-opacity seed.
    scale = column_bottom/140
    bottom = math.log1p(column_bottom/scale)
    photosphere = math.log1p(column_bottom/scale * gray_tau_bottom**-0.25)
    anchor = max(1, min(layers - 2, round((layers - 1)*photosphere/bottom)))
    if deep_layers is not None:
        if not isinstance(deep_layers, int) or not 1 <= deep_layers <= layers-2:
            raise ValueError("Deep layers must leave at least one layer above and below the anchor")
        anchor = layers-1-deep_layers
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


def restart_temperatures(seed, path, interpolate=False):
    """Resume from CSV temperatures on the same radial column-mass grid."""
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"layer", "column_mass_g_cm2", "temperature_K"}
        if not required <= set(reader.fieldnames or ()):
            raise ValueError("Restart CSV requires layer, column mass, and temperature columns")
        rows = list(reader)
    mass = [float(row["column_mass_g_cm2"]) for row in rows]
    temp = [float(row["temperature_K"]) for row in rows]
    if (len(rows) < 2 or any(int(row["layer"]) != i for i, row in enumerate(rows, 1))
            or any(not math.isfinite(m) or m < 0 for m in mass)
            or any(a >= b for a, b in zip(mass, mass[1:]))
            or any(not math.isfinite(t) or t <= 0 for t in temp)):
        raise ValueError("Restart CSV must give finite temperatures on an increasing radial column-mass grid")
    target = [row[-1] for row in seed.thermodynamics]
    if interpolate:
        if ((target[0] < mass[0] and not math.isclose(target[0], mass[0], abs_tol=1e-9))
                or (target[-1] > mass[-1] and not math.isclose(target[-1], mass[-1], rel_tol=1e-9))):
            raise ValueError("The restart must cover the entire requested radial column-mass grid")
        scale = target[-1]/140
        axis = [math.log1p(x/scale) for x in mass]
        log_temp = list(map(math.log, temp))
        temp = []
        for m in target:
            x = max(axis[0], min(axis[-1], math.log1p(m/scale)))
            j = min(bisect_right(axis, x), len(axis)-1)
            fraction = (x-axis[j-1])/(axis[j]-axis[j-1])
            temp.append(math.exp(log_temp[j-1] + fraction*(log_temp[j]-log_temp[j-1])))
    elif len(mass) != len(target) or any(
        not math.isclose(a, b, rel_tol=1e-9) for a, b in zip(mass, target)
    ):
        raise ValueError("Restart CSV must match the current radial column-mass grid")
    structure = tuple((*row[:3], t, *row[4:]) for row, t in zip(seed.structure, temp))
    return replace(seed, structure=structure)


def write_checkpoint(path, initial, model, result):
    """Replace the checkpoint after a complete CSV has been written."""
    temporary = None
    try:
        with NamedTemporaryFile(mode="w", dir=path.parent, encoding="utf-8",
                                newline="", delete=False) as stream:
            temporary = Path(stream.name)
            write_balance_csv(stream, initial, model, result, native_convection=True)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def update_boundary(model, diagnostics, boundary, at_photosphere, radiative_force):
    anchor = photospheric_anchor(model, boundary, diagnostics) if at_photosphere else boundary.anchor
    force = diagnostics[9] if radiative_force else ()
    previous = boundary.radiative_acceleration or (0.0,)*len(force)
    force_error = max((abs(a-b) for a, b in zip(force, previous)), default=0.0)/boundary.gravity
    settled = abs(anchor-boundary.anchor) <= 1e-3 and force_error <= 1e-4
    return replace(boundary, anchor=anchor, radiative_acceleration=force), settled


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wavelengths", type=Path, help="Vacuum wavelength grid in Å")
    parser.add_argument("abundances", type=Path, help="JSON element-to-log-epsilon mapping with H = 12")
    parser.add_argument("--teff", type=float, required=True)
    parser.add_argument("--logg", type=float, required=True)
    parser.add_argument("--mass", type=float, required=True, help="Mass in solar masses")
    parser.add_argument("--lines", type=Path, help="Atomic line-list CSV")
    parser.add_argument("--molecular-cross-section", action="append", default=[], type=Path)
    parser.add_argument("--test-omit-hot-opacity", action="append", default=[], type=Path,
                        metavar="TABLE", help="Stress test: omit this molecular table above its maximum temperature")
    parser.add_argument("--cia", action="append", default=[], type=Path,
                        help="HITRAN collision-induced absorption table, repeatable")
    parser.add_argument("--layers", type=int, default=56)
    parser.add_argument("--deep-layers", type=int,
                        help="Number of grid points below the gray photospheric anchor")
    parser.add_argument("--surface-pressure", type=float, default=1.0, help="Gas pressure in dyn/cm²")
    parser.add_argument("--column-bottom", type=float, default=700.0, help="Bottom column mass in g/cm²")
    parser.add_argument("--gray-tau-bottom", type=float, default=32.0, help="Seed optical depth at bottom")
    parser.add_argument("--vmic", type=float, default=2.0, help="Microturbulence in km/s")
    parser.add_argument("--alpha", type=float, default=1.5, help="Mixing length in pressure scale heights")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--refine-stride", type=int, action="append", default=[],
                        help="Continue on a finer wavelength stride after the initial stage")
    parser.add_argument("--opacity-ramp", type=float, nargs="+", default=[],
                        help="Start with increasing fractions of line opacity before switching to full opacity")
    parser.add_argument("--core-rays", type=int, default=24)
    parser.add_argument("--steps", type=int, default=1,
                        help="Temperature corrections per stage (0 to evaluate the initial profile)")
    parser.add_argument("--upper-knots", type=int, default=12)
    parser.add_argument("--adaptive-knots", type=int, default=0,
                        help="Add up to this many temperature controls at large residuals per stage")
    parser.add_argument("--peak-error-weight", type=float, default=0.25,
                        help="Weight of the maximum error in the temperature step acceptance criterion")
    parser.add_argument("--jacobian-interval", type=int, default=3,
                        help="Steps between finite-difference Jacobian refreshes")
    parser.add_argument("--kernel", type=Path, help="Optional compiled spherical transfer library")
    parser.add_argument("--workers", type=int, default=1,
                        help="Process count for Jacobian perturbations and fixed-profile transfer")
    parser.add_argument("--flux-tolerance", type=float, default=0.01,
                        help="Maximum fractional luminosity error for convergence")
    parser.add_argument("--require-flux-convergence", action="store_true",
                        help="Exit with status 2 if the final luminosity error exceeds the tolerance")
    parser.add_argument("--min-bottom-rosseland-tau", type=float,
                        help="Exit with status 2 if the computed bottom Rosseland depth is smaller")
    parser.add_argument("--outer-heating-weight", type=float, default=0.0,
                        help="Include normalized optically thin radiative heating in the temperature objective")
    parser.add_argument("--heating-transition-tau", type=float,
                        help="Experimental local radiative balance above this depth and flux balance below it")
    parser.add_argument("--bracket-bottom-flux", action="store_true",
                        help="Stabilize the deepest luminosity residual with a bounded temperature correction")
    parser.add_argument("--bracket-deep-layers", type=int, default=1,
                        help="Choose the worst residual in this many deepest layers when bracketing")
    parser.add_argument("--temperature-smoothness", type=float, default=0.0,
                        help="Penalize short-scale log-temperature curvature during the fit")
    parser.add_argument("--output-atlas", type=Path, help="Write an ATLAS DECK6 .dat or .dat.gz record")
    parser.add_argument("--atlas-format", choices=("atlas9", "atlas12"), default="atlas9",
                        help="Format for --output-atlas (default: atlas9)")
    parser.add_argument("--output-model", type=Path, help="Write a MARCS-format .mod or .mod.gz record")
    parser.add_argument("--output-spectrum", type=Path, help="Write the emergent spectrum on the final grid")
    parser.add_argument("--restart-csv", type=Path, help="Start from temperatures on the same mass grid")
    parser.add_argument("--interpolate-restart", action="store_true",
                        help="Interpolate the restart temperatures onto a new mass grid")
    parser.add_argument("--checkpoint-csv", type=Path, help="Save a restart CSV after each accepted correction")
    parser.add_argument("--transfer-cache", type=Path,
                        help="Reuse a matching full-grid calibration from a verified restart")
    parser.add_argument("--adaptive-jacobian", action="store_true",
                        help="Refresh the Jacobian based on measured linearization error")
    parser.add_argument("--verify-full", action="store_true",
                        help="Check the final profile at all input wavelengths")
    parser.add_argument("--full-grid-correction", action="store_true",
                        help="Calibrate the final sampled-grid fit with a full-grid transfer")
    parser.add_argument("--equilibrate-surface", type=int, nargs="?", const=1, default=0, metavar="PASSES",
                        help="Refine the outermost temperature with the full-grid radiation field")
    parser.add_argument("--photospheric-anchor", action="store_true",
                        help="Set the reference gravity at the computed Rosseland photosphere")
    parser.add_argument("--radiative-force", action="store_true",
                        help="Include radiative acceleration in the hydrostatic pressure balance")
    args = parser.parse_args()
    try:
        if not math.isfinite(args.temperature_smoothness) or args.temperature_smoothness < 0:
            raise ValueError("Temperature smoothness weight must be finite and nonnegative")
        if args.workers < 1 or not 0 < args.flux_tolerance < 1:
            raise ValueError("Use at least one worker and a fractional flux tolerance between zero and one")
        if args.equilibrate_surface < 0:
            raise ValueError("Use a nonnegative number of surface correction passes")
        if args.heating_transition_tau is not None and (
            not math.isfinite(args.heating_transition_tau)
            or args.heating_transition_tau <= 0 or args.outer_heating_weight
        ):
            raise ValueError("--heating-transition-tau requires a positive depth and no --outer-heating-weight")
        if args.atlas_format == "atlas12":
            if not args.output_atlas:
                raise ValueError("--atlas-format atlas12 requires --output-atlas")
            if args.layers > 72:
                raise ValueError("The reference ATLAS12 reader supports at most 72 depth records")
        if args.output_atlas:
            reserved = (args.wavelengths, args.abundances, args.lines, args.kernel, args.restart_csv,
                        args.checkpoint_csv, args.output_model, args.output_spectrum,
                        *args.molecular_cross_section, *args.cia)
            if args.output_atlas.resolve() in {path.resolve() for path in reserved if path is not None}:
                raise ValueError("The ATLAS export must use a separate file")
        omitted = {path.resolve() for path in args.test_omit_hot_opacity}
        if omitted:
            if (args.steps != 0 or args.equilibrate_surface or not args.restart_csv
                    or args.photospheric_anchor or args.radiative_force):
                raise ValueError("--test-omit-hot-opacity requires --steps 0 and --restart-csv with structure and temperature updates disabled")
            if args.output_model or args.output_atlas:
                raise ValueError("The hot-opacity stress test supports diagnostic output only")
            if not omitted <= {path.resolve() for path in args.molecular_cross_section}:
                raise ValueError("Each --test-omit-hot-opacity table must also be supplied with --molecular-cross-section")
        minimum_tau = args.min_bottom_rosseland_tau
        if minimum_tau is not None and (not math.isfinite(minimum_tau) or minimum_tau <= 0):
            raise ValueError("Minimum bottom Rosseland depth must be positive and finite")
        abundances = json.loads(args.abundances.read_text(encoding="utf-8"))
        seed, boundary = initialize(
            args.teff, args.logg, args.mass, abundances, layers=args.layers,
            surface_pressure=args.surface_pressure, column_bottom=args.column_bottom,
            gray_tau_bottom=args.gray_tau_bottom, vmic=args.vmic, alpha=args.alpha,
            deep_layers=args.deep_layers,
        )
        reference = read_wavelengths(args.wavelengths)
        strides = [args.stride, *args.refine_stride]
        if any(a <= b for a, b in zip(strides, strides[1:])):
            raise ValueError("Refinement strides must decrease at each stage")
        opacity_scales = [*args.opacity_ramp, 1.0]
        if any(not 0 <= a < b for a, b in zip(opacity_scales, opacity_scales[1:])):
            raise ValueError("Use strictly increasing opacity ramp fractions from zero up to one")
        stages = [(args.stride, scale) for scale in args.opacity_ramp]
        stages.extend((stride, 1.0) for stride in strides)
        fingerprint = None
        if args.transfer_cache:
            opacity_paths = [args.lines, args.kernel, *args.molecular_cross_section, None, *args.cia]
            reserved = [args.wavelengths, args.abundances, args.restart_csv,
                        args.checkpoint_csv, args.output_model, args.output_atlas, args.output_spectrum,
                        *opacity_paths]
            if args.transfer_cache.resolve() in {path.resolve() for path in reserved if path is not None}:
                raise ValueError("The transfer cache must be separate from input and output files")
            try:
                with args.transfer_cache.open("rb") as stream:
                    recognized = stream.read(19).startswith(b'{"payload":{"key":')
            except FileNotFoundError:
                recognized = True
            except OSError:
                recognized = False
            if recognized:
                fingerprint = input_fingerprint(opacity_paths)
            else:
                print("Ignoring transfer cache path without a recognized cache file", file=sys.stderr)
                args.transfer_cache = None

        lines = read_lines(args.lines) if args.lines else ()
        molecular = [
            MolecularCrossSection(path, omit_hot=path.resolve() in omitted)
            for path in args.molecular_cross_section
        ]
        for table in molecular:
            if table.omit_hot:
                print(f"Stress test: omitting {table.isotopologue} cross sections above "
                      f"{table.temperatures[-1]:g} K", file=sys.stderr)
        cia = [HitranCIA(path) for path in args.cia]
        known_cia = boundary.chemistry.atoms.keys() | boundary.chemistry.molecules.keys()
        if any(not set(table.species) <= known_cia for table in cia):
            raise ValueError("CIA partners must be neutral atoms or molecules in the LTE chemistry table")
        if len({tuple(sorted(table.species)) for table in cia}) != len(cia):
            raise ValueError("Provide one CIA table per collision pair")
        tables = [*molecular, *cia]
        if args.interpolate_restart and not args.restart_csv:
            raise ValueError("--interpolate-restart requires --restart-csv")
        if args.restart_csv:
            seed = restart_temperatures(seed, args.restart_csv, args.interpolate_restart)
        model = seed
        def cache_key():
            return transfer_key(fingerprint, model, reference, args.core_rays, args.vmic,
                                args.workers, boundary.radiative_acceleration)

        def checkpoint(current, diagnostics):
            write_checkpoint(args.checkpoint_csv, seed, current, diagnostics)

        last_full_transfer = None

        def transfer_on_grid(surface_sample=None, line_scale=1.0, grid=None, surface_mean_sample=None):
            nonlocal last_full_transfer
            grid = wavelengths if grid is None else grid
            if args.workers > 1:
                transfer = evaluate_parallel(
                    model, grid, lines, {}, tables, args.core_rays, args.vmic, boundary.chemistry,
                    args.kernel, args.workers, boundary.radiative_acceleration, line_scale,
                    surface_mean_sample, surface_sample,
                )
            else:
                transfer = evaluate(
                    model, grid, lines, {}, tables, args.core_rays, args.vmic,
                    chemistry=boundary.chemistry, kernel_path=args.kernel, native_convection=True,
                    hydrostatic_force=boundary.radiative_acceleration, surface_sample=surface_sample,
                    line_scale=line_scale, surface_mean_sample=surface_mean_sample,
                )
            if args.transfer_cache and grid == reference and line_scale == 1.0:
                last_full_transfer = transfer
            return transfer

        update_structure = args.photospheric_anchor or args.radiative_force
        for index, (stride, line_scale) in enumerate(stages, 1):
            wavelengths = nested_wavelengths(reference, stride)
            if args.steps == 0 and args.equilibrate_surface and not update_structure:
                model, _ = reconstruct_for_transfer(model, boundary=boundary)
                break
            print(f"Solving on {len(wavelengths):,} wavelengths (stride {stride}, "
                  f"line opacity {line_scale:.0%})",
                  file=sys.stderr, flush=True)
            flux_offset = heating_offset = heating_scale = ()
            if args.full_grid_correction and index == len(stages) and args.steps > 0 and wavelengths != reference:
                print(f"Calibrating the residuals on {len(reference):,} wavelengths",
                      file=sys.stderr, flush=True)
                model, _ = reconstruct_for_transfer(model, boundary=boundary)
                sampled_transfer = transfer_on_grid()
                cached = (read_residuals(args.transfer_cache, cache_key(), len(model.structure))
                          if args.transfer_cache else None)
                if cached is None:
                    exact_transfer = transfer_on_grid(grid=reference)
                    exact_flux = tuple(row[-1] for row in total_flux_budget(model, exact_transfer, True))
                    exact_heating = exact_transfer[8]
                else:
                    exact_flux, exact_heating = cached
                    print("Reused verified full-grid calibration", file=sys.stderr, flush=True)
                sampled = total_flux_budget(model, sampled_transfer, True)
                flux_offset = tuple(exact-row[-1] for row, exact in zip(sampled, exact_flux))
                if args.outer_heating_weight or args.heating_transition_tau is not None:
                    heating_scale, heating_offset = calibrate_heating(sampled_transfer[8], exact_heating)
                    print(f"Full-grid starting surface net heating / LTE emission "
                          f"{exact_heating[0]:+.3%}", file=sys.stderr, flush=True)
                print("Iteration errors use the calibrated full-grid residuals",
                      file=sys.stderr, flush=True)
            if args.steps == 0 and args.workers > 1:
                model, _ = reconstruct_for_transfer(model, boundary=boundary)
                result = transfer_on_grid(line_scale=line_scale)
                errors = [row[-1] for row in total_flux_budget(model, result, True)]
                rms = math.sqrt(math.fsum(x*x for x in errors)/len(errors))
                report_balance_step(0, (rms, max(map(abs, errors)), errors[0], errors[-1]))
            else:
                model, result, _, _, _ = balance_all_layers(
                    model, wavelengths, lines, {}, steps=args.steps, tables=tables,
                    core_rays=args.core_rays, vmic=args.vmic, kernel_path=args.kernel,
                    upper_knots=args.upper_knots, boundary=boundary, progress=report_balance_step,
                    flux_tolerance=args.flux_tolerance, jacobian_interval=args.jacobian_interval,
                    outer_heating_weight=args.outer_heating_weight,
                    heating_transition_tau=args.heating_transition_tau,
                    on_accept=checkpoint if args.checkpoint_csv else None,
                    adaptive_knots=args.adaptive_knots, peak_weight=args.peak_error_weight,
                    line_scale=line_scale, jacobian_workers=args.workers,
                    flux_offset=flux_offset, heating_offset=heating_offset, heating_scale=heating_scale,
                    temperature_smoothness=args.temperature_smoothness,
                    bracket_bottom_flux=args.bracket_bottom_flux,
                    bracket_deep_layers=args.bracket_deep_layers,
                    adaptive_jacobian=args.adaptive_jacobian,
                    jacobian_progress=lambda count: print(
                        f"Refreshing Jacobian ({count} atmosphere perturbations)", file=sys.stderr, flush=True
                    ),
                )
            if update_structure and index < len(stages):
                boundary, _ = update_boundary(model, result, boundary, args.photospheric_anchor,
                                              args.radiative_force)
        if ((args.verify_full or args.full_grid_correction) and wavelengths != reference
                and (not args.equilibrate_surface or update_structure)):
            print(f"Verifying on {len(reference):,} wavelengths", file=sys.stderr, flush=True)
            wavelengths = reference

            result = transfer_on_grid()
        if update_structure:
            for _ in range(4):
                updated, settled = update_boundary(model, result, boundary, args.photospheric_anchor,
                                                   args.radiative_force)
                if settled:
                    break
                boundary = updated
                model, _ = reconstruct_for_transfer(model, boundary=boundary)
                result = transfer_on_grid()
            _, settled = update_boundary(model, result, boundary, args.photospheric_anchor,
                                         args.radiative_force)
            print(f"Boundary updates {'settled' if settled else 'reached the iteration limit'}",
                  file=sys.stderr)
        if args.equilibrate_surface:
            wavelengths = reference
            field = []
            print(f"Solving the outermost thermal balance on {len(reference):,} wavelengths",
                  file=sys.stderr, flush=True)
            result = transfer_on_grid(surface_mean_sample=lambda wave, nu, weight, mean: field.append(mean))
            for _ in range(args.equilibrate_surface):
                temperature = surface_equilibrium_temperature(model, wavelengths, field, lines, tables,
                                                              args.vmic, boundary.chemistry)
                if temperature is None:
                    print("Surface balance was not bracketed within the opacity temperature range", file=sys.stderr)
                    break
                if temperature == model.structure[0][3]:
                    break
                original = model
                top = model.structure[0]
                profile = replace(model, structure=((*top[:3], temperature, *top[4:]), *model.structure[1:]))
                model, _ = reconstruct_for_transfer(profile, boundary=boundary)
                field = []
                trial = transfer_on_grid(surface_mean_sample=lambda wave, nu, weight, mean: field.append(mean))
                old_max = max(abs(row[-1]) for row in total_flux_budget(original, result, True))
                new_max = max(abs(row[-1]) for row in total_flux_budget(model, trial, True))
                if abs(trial[8][0]) < abs(result[8][0]) and new_max <= max(args.flux_tolerance, old_max):
                    print(f"Accepted surface temperature {temperature:.2f} K "
                          f"(net heating {result[8][0]:+.3%} to {trial[8][0]:+.3%})", file=sys.stderr)
                    result = trial
                else:
                    print("Surface trial failed full-grid heating or luminosity acceptance", file=sys.stderr)
                    model = original
                    break
        if args.photospheric_anchor:
            measured = photospheric_anchor(model, boundary, result)
            print(f"Rosseland photosphere reference index {boundary.anchor+1:.4f}, "
                  f"measured {measured+1:.4f}", file=sys.stderr)
        if args.radiative_force:
            print(f"Peak radiative acceleration / reference gravity {max(result[9])/boundary.gravity:.3%}",
                  file=sys.stderr)
        if args.output_spectrum:
            with args.output_spectrum.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(("wavelength_vacuum_A", "frequency_Hz", "frequency_weight_Hz",
                                 "reference_flux_per_Hz", "reference_flux_per_A"))

                def sample(wavelength, frequency, weight, flux):
                    writer.writerow((wavelength, frequency, weight, flux, flux*frequency/wavelength))

                result = transfer_on_grid(sample)
        budget = write_balance_csv(sys.stdout, seed, model, result, native_convection=True)
        if (args.transfer_cache and wavelengths == reference and not args.test_omit_hot_opacity
                and (args.workers == 1 or result is last_full_transfer)):
            try:
                if fingerprint == input_fingerprint(opacity_paths):
                    write_residuals(args.transfer_cache, cache_key(), tuple(row[-1] for row in budget), result[8])
            except (OSError, ValueError) as exc:
                print(f"Could not save full-grid calibration: {exc}", file=sys.stderr)
        maximum = max(abs(row[-1]) for row in budget)
        status = "converged" if maximum <= args.flux_tolerance else "unconverged"
        print(f"Luminosity balance {status} at {args.flux_tolerance:.3%} tolerance "
              f"(maximum error {maximum:.3%})", file=sys.stderr)
        print(f"Surface net heating / LTE emission {result[8][0]:+.3%}", file=sys.stderr)
        final_depth = rosseland_depth(model, boundary, result)
        bottom_tau = final_depth[-1]
        if args.heating_transition_tau is not None:
            thin_heating = [abs(q) for tau, q in zip(final_depth, result[8])
                            if tau < args.heating_transition_tau]
            if thin_heating:
                print(f"Peak optically thin net heating / LTE emission {max(thin_heating):.3%} "
                      f"({len(thin_heating)} layers)", file=sys.stderr)
        print(f"Bottom Rosseland optical depth {bottom_tau:.6g}", file=sys.stderr)
        shallow = minimum_tau is not None and bottom_tau < minimum_tau
        if shallow:
            print(f"Bottom optical depth is below the required minimum {minimum_tau:g}", file=sys.stderr)
        for table in molecular:
            cold, hot, low, high = table.clamped_layers(model.structure)
            if any((cold, hot, low, high)):
                hot_policy = "omitted" if table.omit_hot else "clamped"
                print(f"{table.isotopologue} layers outside tabulated T "
                      f"(cold {cold}, hot {hot} {hot_policy}) or P (low {low}, high {high})",
                      file=sys.stderr)
        for table in cia:
            outside = sum(not table.temperatures[0] <= row[3] <= table.temperatures[-1]
                          for row in model.structure)
            if outside:
                print(f"{'-'.join(table.species)} CIA temperature clamped in {outside} layers",
                      file=sys.stderr)
        if args.output_model:
            document = build_document(model, result, lines, tables, boundary, args.vmic)
            header = (f"{document.header[0]} ({status}, max luminosity error {maximum:.3%})",
                      *document.header[1:])
            document = replace(document, header=header)
            write_document(document, args.output_model)
        if args.output_atlas:
            write_atlas(model, result, boundary, args.vmic, args.output_atlas, status=status,
                        atlas12=args.atlas_format == "atlas12")
        return 2 if shallow or (args.require_flux_convergence and maximum > args.flux_tolerance) else 0
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    raise SystemExit(main())
