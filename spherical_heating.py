"""Coherent spherical continuum scattering and LTE heating on a fixed MARCS profile."""

import argparse
import csv
from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
import sys

from benchmark_sampling import nested_wavelengths
from gray_atmosphere import propagate
from lte_chemistry import Thermochemistry, _linear_solve
from lte_opacity import layer_opacities, prepare_layers, read_lines
from marcs_writer import read_document
from molecular_cross_section import MolecularCrossSection
from native_transfer import read_wavelengths, reference_flux_budget, spectral_quadrature
from odf_transfer import planck_nu
from spherical_transfer import disk_rays, impact_parameters


@dataclass(frozen=True)
class SphericalRay:
    core: bool
    chord_scale: tuple
    mu_bottom: float
    mean_weights: tuple
    flux_weights: tuple


def angular_geometry(radii, core_rays):
    """Tabulate shell-by-shell angular weights for the shared impact rays."""
    disk = disk_rays(radii, core_rays)
    squared, impacts = impact_parameters(radii, core_rays)
    j_weights = [[0.0] * (len(ray[1]) + 1) for ray in disk]
    f_weights = [[0.0] * (len(ray[1]) + 1) for ray in disk]
    for layer, r2 in enumerate(squared):
        hits = [(i, p2) for i, (p2, last, _) in enumerate(impacts) if layer <= last]
        for (i, a), (j, b) in zip(hits, hits[1:]):
            flux_weight = math.pi * (b - a) / (2 * r2)
            # J = 1/2 integral of the sum of inward and outward intensities.
            mean_weight = (math.sqrt(max(0, 1 - a/r2)) - math.sqrt(max(0, 1 - b/r2))) / 4
            for ray in (i, j):
                j_weights[ray][layer] += mean_weight
                f_weights[ray][layer] += flux_weight
    return tuple(
        SphericalRay(core, tuple(scale), mu_bottom, tuple(jw), tuple(fw))
        for (core, scale, _, mu_bottom, _), jw, fw in zip(disk, j_weights, f_weights)
    )


def _endpoint_weight(distance):
    if distance < 1e-4:
        return distance * (0.5 - distance / 6 + distance**2 / 24)
    return 1 + math.expm1(-distance) / distance


class FormalTransfer:
    """Linear spherical formal solution at one wavelength."""

    def __init__(self, rays, optical_depths, thermal):
        self.rays = []
        self.thermal = thermal
        diagonal = [0.0] * len(thermal)
        for ray in rays:
            distance = [tau * scale for tau, scale in zip(optical_depths, ray.chord_scale)]
            endpoint = [_endpoint_weight(step) for step in distance]
            bottom = thermal[-1] if ray.core else None
            if ray.core and sum(distance) >= 10 and optical_depths[-1] > 0:
                bottom += ray.mu_bottom * (thermal[-1] - thermal[-2]) / optical_depths[-1]
            self.rays.append((ray, distance, bottom))
            last = len(distance)
            distance_to_tangent = 0.0
            for i in range(last, -1, -1):
                inbound = endpoint[i - 1] if i else 0.0
                outbound = endpoint[i] if i < last else 0.0
                if not ray.core:
                    # A local source perturbation also travels around the tangent.
                    if i == last:
                        outbound = inbound
                    else:
                        step = distance[i]
                        start_weight = -math.expm1(-step) - endpoint[i]
                        outbound += start_weight * math.exp(-step - 2 * distance_to_tangent)
                        distance_to_tangent += step
                        outbound += inbound * math.exp(-2 * distance_to_tangent)
                diagonal[i] += ray.mean_weights[i] * (inbound + outbound)
        self.diagonal = diagonal

    def moments(self, source, include_flux=True):
        mean = [0.0] * len(source)
        flux = [0.0] * len(source)
        for ray, distance, bottom in self.rays:
            last = len(distance)
            inward = [0.0]
            for i, step in enumerate(distance):
                inward.append(propagate(inward[-1], source[i], source[i + 1], step))
            outward = inward[-1] if bottom is None else bottom
            for i in range(last, -1, -1):
                if i < last:
                    outward = propagate(outward, source[i + 1], source[i], distance[i])
                mean[i] += ray.mean_weights[i] * (outward + inward[i])
                if include_flux:
                    flux[i] += ray.flux_weights[i] * (outward - inward[i])
        return mean, flux


def solve_wavelength(rays, radial_widths, density, thermal, absorption, scattering,
                     tolerance=1e-7, max_iterations=200):
    """Solve isotropic coherent scattering and return local true-absorption heating."""
    total_opacity = [a + s for a, s in zip(absorption, scattering)]
    extinction = [k * rho for k, rho in zip(total_opacity, density)]
    steps = [width * (a + b) / 2 for width, a, b in zip(radial_widths, extinction, extinction[1:])]
    transfer = FormalTransfer(rays, steps, thermal)
    scatter_fraction = [s / k if k else 0.0 for s, k in zip(scattering, total_opacity)]
    source = list(thermal)
    iterations = 0
    if any(scatter_fraction):
        # Diagonal accelerated lambda iteration. The boundary intensities
        # remain fixed while updating the scattering source.
        floor = max(thermal) * 1e-12
        if floor > 0:
            for iterations in range(1, max_iterations + 1):
                mean, _ = transfer.moments(source, include_flux=False)
                trial = [
                    current + ((1 - w) * b + w * j - current) / (1 - w * diag)
                    for current, b, w, j, diag in
                    zip(source, thermal, scatter_fraction, mean, transfer.diagonal)
                ]
                if max(abs(new - old) / max(abs(new), floor) for old, new in zip(source, trial)) <= tolerance:
                    source = trial
                    break
                source = trial
            else:
                raise RuntimeError(f"Spherical scattering failed to converge in {max_iterations} iterations")
    mean, flux = transfer.moments(source)
    heating = [4 * math.pi * k * (j - b) for k, j, b in zip(absorption, mean, thermal)]
    return mean, flux, heating, iterations


def shell_energy_budget(radii, density, flux, pointwise_heating):
    """Return conservative shell heating and a pointwise quadrature diagnostic."""
    areas = [4 * math.pi * radius**2 for radius in radii]
    luminosity = [area * f for area, f in zip(areas, flux)]
    mass_per_radius = [area * rho for area, rho in zip(areas, density)]
    shells = []
    pointwise_integral = 0.0
    for i, (a, b) in enumerate(zip(radii, radii[1:])):
        left = 0.5 * (a - b) * mass_per_radius[i]
        right = 0.5 * (a - b) * mass_per_radius[i + 1]
        net_heating = luminosity[i + 1] - luminosity[i]
        shells.append((net_heating / (left + right), net_heating))
        pointwise_integral += left * pointwise_heating[i] + right * pointwise_heating[i + 1]
    scale = max(map(abs, luminosity))
    mismatch = pointwise_integral - (luminosity[-1] - luminosity[0])
    return luminosity, shells, mismatch / scale if scale else 0.0


def evaluate(model, wavelengths, lines, composition, tables=(), core_rays=12, vmic=2.0,
             include_scattering=True, chemistry=None):
    """Integrate spherical flux and heating on the supplied wavelength grid."""
    if float(model.header[5].split()[0]) <= 0:
        raise ValueError("A spherical MARCS atmosphere is required")
    reference_radius = float(model.header[7].split()[0])
    radii = [reference_radius - row[2] for row in model.structure]
    rays = angular_geometry(radii, core_rays)
    widths = [a - b for a, b in zip(radii, radii[1:])]
    frequencies, weights = spectral_quadrature(wavelengths, "wavelength")
    layers = list(prepare_layers(model, lines, composition, vmic, chemistry))
    temperatures = [layer[0].temperature for layer in layers]
    density = [layer[0].mass_density for layer in layers]
    flux = [0.0] * len(layers)
    heating = [0.0] * len(layers)
    max_iterations = 0
    for start in range(0, len(wavelengths), 128):
        end = start + 128
        sample = wavelengths[start:end]
        grids = [(table, table.on_grid(sample)) for table in tables]
        opacities = [list(layer_opacities(layer, sample, grids)) for layer in layers]
        for columns, nu, weight in zip(zip(*opacities), frequencies[start:end], weights[start:end]):
            absorption = [line + continuum for line, continuum, _ in columns]
            scattering = [scatter if include_scattering else 0.0 for _, _, scatter in columns]
            thermal = [planck_nu(nu, t) for t in temperatures]
            try:
                _, f, q, iterations = solve_wavelength(rays, widths, density, thermal, absorption, scattering)
            except RuntimeError as exc:
                raise RuntimeError(f"At {nu:.8g} Hz, {exc}") from exc
            max_iterations = max(max_iterations, iterations)
            for i, (local_flux, local_heating) in enumerate(zip(f, q)):
                flux[i] += weight * local_flux
                heating[i] += weight * local_heating
    luminosity, shells, mismatch = shell_energy_budget(radii, density, flux, heating)
    return flux, heating, luminosity, shells, max_iterations, mismatch


def correct_temperatures(model, wavelengths, lines, composition, count, steps=1,
                         tables=(), core_rays=12, vmic=2.0, include_scattering=True):
    """Reduce outer-shell total luminosity divergence at fixed gas pressure."""
    if not 1 <= count < len(model.structure) or steps < 1:
        raise ValueError("Select outer layers to adjust and leave at least one deeper layer fixed")
    chemistry = Thermochemistry()
    reference_radius = float(model.header[7].split()[0])
    target = 4 * math.pi * reference_radius**2 * float(model.header[2].split()[0])
    convection = [target * row[5] for row in model.thermodynamics]
    baseline = model.structure

    def trial(temperatures):
        structure = tuple((*row[:3], t, *row[4:])
                          for row, t in zip(baseline, temperatures)) + baseline[count:]
        current = replace(model, structure=structure)
        result = evaluate(current, wavelengths, lines, composition, tables, core_rays,
                          vmic, include_scattering, chemistry)
        luminosity = [r + c for r, c in zip(result[2], convection)]
        residual = [(b - a) / target for a, b in zip(luminosity[:count], luminosity[1:count+1])]
        surface = luminosity[0] / target - 1
        return current, result, residual, surface

    temperatures = [row[3] for row in baseline[:count]]
    current, result, residual, surface = trial(temperatures)
    history = [(max(map(abs, residual)), surface)]
    for _ in range(steps):
        if history[-1][0] < 1e-8:
            break
        columns = []
        for i, value in enumerate(temperatures):
            delta = 0.001 * value
            if value + delta > chemistry.temperatures[-1]:
                delta = -delta
            perturbed = temperatures.copy()
            perturbed[i] += delta
            shifted = trial(perturbed)[2]
            columns.append([(b - a) / delta for a, b in zip(residual, shifted)])
        direction = _linear_solve([list(row) for row in zip(*columns)], [-x for x in residual])
        # Keep the Newton direction and cap the largest fractional change.
        fractional = max(abs(d) / t for d, t in zip(direction, temperatures))
        scale = min(1.0, 0.03 / fractional) if fractional else 1.0
        for _ in range(8):
            candidate = [t + scale * d for t, d in zip(temperatures, direction)]
            if all(chemistry.temperatures[0] <= t <= chemistry.temperatures[-1] for t in candidate):
                new_model, new_result, new_residual, new_surface = trial(candidate)
                worst = max(map(abs, new_residual))
                if worst < history[-1][0]:
                    temperatures = candidate
                    current, result, residual = new_model, new_result, new_residual
                    history.append((worst, new_surface))
                    break
            scale *= 0.5
        else:
            break
    return current, result, history


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model")
    parser.add_argument("lines")
    parser.add_argument("composition")
    parser.add_argument("wavelengths")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--core-rays", type=int, default=12)
    parser.add_argument("--vmic", type=float, default=2)
    parser.add_argument("--absorption-only", action="store_true")
    parser.add_argument("--molecular-cross-section", action="append", default=[])
    parser.add_argument("--adjust-layers", type=int, default=0)
    parser.add_argument("--steps", type=int, default=1)
    args = parser.parse_args()
    if not math.isfinite(args.vmic) or args.vmic < 0 or args.adjust_layers < 0 or args.steps < 1:
        parser.error("Use nonnegative microturbulence and layer count, and at least one step")
    try:
        model = read_document(args.model)
        reference = read_wavelengths(args.wavelengths)
        wavelengths = nested_wavelengths(reference, args.stride)
        lines = read_lines(args.lines)
        with Path(args.composition).open(encoding="utf-8") as stream:
            composition = json.load(stream)
        tables = [MolecularCrossSection(path) for path in args.molecular_cross_section]
        initial = model
        if args.adjust_layers:
            model, result, history = correct_temperatures(
                model, wavelengths, lines, composition, args.adjust_layers, args.steps,
                tables, args.core_rays, args.vmic, include_scattering=not args.absorption_only,
            )
            for step, (residual, surface) in enumerate(history):
                print(f"Step {step}: max outer-shell luminosity residual {residual:.4e}, "
                      f"surface total-flux error {surface:+.3%}", file=sys.stderr)
        else:
            result = evaluate(model, wavelengths, lines, composition, tables,
                              args.core_rays, args.vmic, include_scattering=not args.absorption_only)
        flux, heating, luminosity, shells, iterations, mismatch = result
        budget = list(reference_flux_budget(model, flux))
    except (IndexError, OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("layer", "initial_temperature_K", "temperature_K", "column_mass_g_cm2", "radiative_flux_erg_cm2_s",
                     "pointwise_heating_erg_g_s", "radiative_luminosity_erg_s",
                     "shell_to_next_heating_erg_g_s", "shell_to_next_heating_erg_s",
                     "marcs_convective_flux_erg_cm2_s", "total_flux_erg_cm2_s",
                     "target_flux_erg_cm2_s", "relative_total_flux_error"))
    for i, (old, row, thermo, f, q, lum, shell, diagnostic) in enumerate(
        zip(initial.structure, model.structure, model.thermodynamics, flux, heating, luminosity,
            [*shells, ("", "")], budget), 1
    ):
        writer.writerow((i, old[3], row[3], thermo[-1], f, q, lum, *shell, *diagnostic))
    print(f"Surface total-flux error: {budget[0][-1]:+.3%}", file=sys.stderr)
    print(f"Maximum scattering iterations per wavelength: {iterations}", file=sys.stderr)
    print(f"Relative mismatch from trapezoidal integration of pointwise heating: {mismatch:+.3e}", file=sys.stderr)


if __name__ == "__main__":
    main()
