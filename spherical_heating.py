"""Coherent spherical continuum scattering and LTE heating on a fixed MARCS profile."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import dataclass, replace
import json
import math
from multiprocessing import get_context
from pathlib import Path
import sys

from benchmark_sampling import nested_wavelengths
from gray_atmosphere import propagate
from lte_chemistry import Thermochemistry, _linear_solve
from hminus import K_EV
from lte_opacity import H_EV_S, LIGHT, layer_opacities, prepare_layers, read_lines
from mixing_length import evaluate as evaluate_convection
from marcs_writer import read_document
from molecular_cross_section import MolecularCrossSection
from native_transfer import read_wavelengths, reference_flux_budget, spectral_quadrature
from odf_transfer import planck_nu
from spherical_hydrostatic import boundary_from_marcs, increments, integrate, reconstruct_for_transfer
from spherical_kernel_loader import CompiledKernel, ScatteringIterationLimit
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
            if ray.core and optical_depths[-1] > 0:
                # Blend the diffusion boundary over the optical-depth
                # transition so a temperature perturbation cannot switch it abruptly.
                x = max(0.0, min(1.0, (sum(distance) - 5)/10))
                increment = ray.mu_bottom * (thermal[-1] - thermal[-2]) / optical_depths[-1]
                bottom = max(0.0, bottom + x*x*(3-2*x)*increment)
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


def direct_scattering_source(transfer, thermal, fraction):
    """Solve the linear scattering equation when iterative transfer stalls."""
    n = len(thermal)
    empty = [0.0] * n
    background, _ = transfer.moments(empty, include_flux=False)
    rhs = [(1 - w)*b + w*j for w, b, j in zip(fraction, thermal, background)]
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        unit = empty.copy()
        unit[i] = 1.0
        response, _ = transfer.moments(unit, include_flux=False)
        for j, (w, base, lit) in enumerate(zip(fraction, background, response)):
            matrix[j][i] = (i == j) - w*(lit - base)
    return _linear_solve(matrix, rhs)


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
                source = direct_scattering_source(transfer, thermal, scatter_fraction)
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
             include_scattering=True, chemistry=None, kernel_path=None, native_convection=False,
             hydrostatic_force=(), surface_sample=None, line_scale=1.0):
    """Integrate spherical flux and heating on the supplied wavelength grid."""
    if float(model.header[5].split()[0]) <= 0:
        raise ValueError("A spherical MARCS atmosphere is required")
    reference_radius = float(model.header[7].split()[0])
    radii = [reference_radius - row[2] for row in model.structure]
    rays = angular_geometry(radii, core_rays)
    solver = CompiledKernel(kernel_path, rays) if kernel_path else solve_wavelength
    widths = [a - b for a, b in zip(radii, radii[1:])]
    frequencies, weights = spectral_quadrature(wavelengths, "wavelength")
    layers = list(prepare_layers(model, lines, composition, vmic, chemistry))
    temperatures = [layer[0].temperature for layer in layers]
    density = [layer[0].mass_density for layer in layers]
    flux = [0.0] * len(layers)
    heating = [0.0] * len(layers)
    max_iterations = 0
    rosseland_weight = [0.0] * len(layers)
    rosseland_inverse = [0.0] * len(layers)
    emission = [0.0] * len(layers)
    radiative_acceleration = [0.0] * len(layers)
    for start in range(0, len(wavelengths), 128):
        end = start + 128
        sample = wavelengths[start:end]
        grids = [(table, table.on_grid(sample)) for table in tables]
        opacities = [list(layer_opacities(layer, sample, grids)) for layer in layers]
        for wavelength, columns, nu, weight in zip(sample, zip(*opacities),
                                                    frequencies[start:end], weights[start:end]):
            absorption = [line_scale*line + continuum for line, continuum, _ in columns]
            scattering = [scatter if include_scattering else 0.0 for _, _, scatter in columns]
            thermal = [planck_nu(nu, t) for t in temperatures]
            if native_convection:
                for i, (temperature, source, absorb, (_, _, scatter)) in enumerate(
                    zip(temperatures, thermal, absorption, columns)
                ):
                    x = H_EV_S * nu / (K_EV * temperature)
                    weighted = weight * source * x / (-math.expm1(-x))
                    rosseland_weight[i] += weighted
                    rosseland_inverse[i] += weighted / (absorb + scatter)
                    emission[i] += weight * 4 * math.pi * absorb * source
            try:
                _, f, q, iterations = solver(rays, widths, density, thermal, absorption, scattering)
            except ScatteringIterationLimit:
                _, f, q, _ = solve_wavelength(rays, widths, density, thermal, absorption,
                                              scattering, max_iterations=0)
                iterations = 200
            except RuntimeError as exc:
                raise RuntimeError(f"At {nu:.8g} Hz, {exc}") from exc
            if surface_sample is not None:
                surface_sample(wavelength, nu, weight, f[0]*(radii[0]/reference_radius)**2)
            max_iterations = max(max_iterations, iterations)
            for i, (local_flux, local_heating) in enumerate(zip(f, q)):
                flux[i] += weight * local_flux
                heating[i] += weight * local_heating
                if native_convection:
                    radiative_acceleration[i] += weight*local_flux*(absorption[i]+scattering[i])/LIGHT
    luminosity, shells, mismatch = shell_energy_budget(radii, density, flux, heating)
    result = (flux, heating, luminosity, shells, max_iterations, mismatch)
    if native_convection:
        means = [a / b for a, b in zip(rosseland_weight, rosseland_inverse)]
        convection = tuple(evaluate_convection(model, composition, means, chemistry, hydrostatic_force))
        normalized_heating = tuple(q/e if e else 0.0 for q, e in zip(heating, emission))
        return (*result, tuple(row[11] for row in convection), convection, normalized_heating,
                tuple(radiative_acceleration), (rosseland_weight, rosseland_inverse, emission, density))
    return result


def evaluate_parallel(model, wavelengths, lines, mixture, tables, core_rays, vmic,
                      chemistry, kernel_path, workers, hydrostatic_force=(), line_scale=1.0):
    """Integrate a fixed native atmosphere over overlapping wavelength bands."""
    if len(wavelengths) < 2 or workers < 2:
        raise ValueError("Parallel transfer requires two or more wavelengths and workers")
    count = min(workers, len(wavelengths)-1)
    indices = [(len(wavelengths)-1)*i//count for i in range(count+1)]
    with ProcessPoolExecutor(max_workers=count, mp_context=get_context("spawn")) as executor:
        futures = [
            executor.submit(evaluate, model, wavelengths[a:b+1], lines, mixture, tables,
                            core_rays, vmic, True, chemistry, kernel_path, True,
                            hydrostatic_force, None, line_scale)
            for a, b in zip(indices, indices[1:])
        ]
        parts = [future.result() for future in futures]
    # The shared endpoints carry half-weights for the two neighboring
    # intervals. Their sum yields the original wavelength quadrature.
    flux = [math.fsum(values) for values in zip(*(part[0] for part in parts))]
    heating = [math.fsum(values) for values in zip(*(part[1] for part in parts))]
    acceleration = [math.fsum(values) for values in zip(*(part[9] for part in parts))]
    numerator, denominator, emission = (
        [math.fsum(values) for values in zip(*(part[10][i] for part in parts))]
        for i in range(3)
    )
    density = parts[0][10][3]
    radius = float(model.header[7].split()[0])
    radii = [radius-row[2] for row in model.structure]
    luminosity, shells, mismatch = shell_energy_budget(radii, density, flux, heating)
    means = [a/b for a, b in zip(numerator, denominator)]
    convection = tuple(evaluate_convection(model, mixture, means, chemistry, hydrostatic_force))
    normalized = tuple(q/e if e else 0.0 for q, e in zip(heating, emission))
    return (flux, heating, luminosity, shells, max(part[4] for part in parts), mismatch,
            tuple(row[11] for row in convection), convection, normalized, acceleration,
            (numerator, denominator, emission, density))


def correct_temperatures(model, wavelengths, lines, composition, count, steps=1,
                         tables=(), core_rays=12, vmic=2.0, include_scattering=True,
                         luminosity_target=False, kernel_path=None, hydrostatic=False,
                         native_convection=False, deep_layers=0):
    """Reduce outer-shell luminosity divergence, optionally updating hydrostatic structure."""
    if not 1 <= count < len(model.structure) or steps < 1 or deep_layers < 0:
        raise ValueError("Invalid outer layers, deep layers or iteration count")
    if count + deep_layers >= len(model.structure) or (deep_layers and not native_convection):
        raise ValueError("Leave an interior anchor and use native convection for deep adjustments")
    chemistry = Thermochemistry()
    boundary = boundary_from_marcs(model, chemistry, composition)[0] if hydrostatic else None
    reference_radius = boundary.reference_radius if boundary else float(model.header[7].split()[0])
    target = 4 * math.pi * reference_radius**2 * float(model.header[2].split()[0])
    convection = [target * row[5] for row in model.thermodynamics]
    baseline = model.structure
    # Adjust a photospheric temperature mode independently of the outer
    # layers. Its support tapers to zero before the convective interior.
    deep_start = len(baseline) - deep_layers
    mode = [0.0] * count + [
        max(0.0, 1.0 - abs(row[0]) / 1.5) for row in baseline[count:deep_start]
    ] + [0.0] * deep_layers
    pivot = count + deep_layers
    if luminosity_target and not any(mode):
        raise ValueError("Leave photospheric layers near Rosseland optical depth one free for the luminosity mode")

    def temperature_profile(temperatures):
        amplitude = temperatures[-1] if luminosity_target else 0.0
        return temperatures[:count] + [
            row[3] + amplitude * w
            for row, w in zip(baseline[count:deep_start], mode[count:deep_start])
        ] + temperatures[count:pivot]

    def trial(temperatures):
        structure = tuple((*row[:3], t, *row[4:])
                          for row, t in zip(baseline, temperature_profile(temperatures)))
        current = replace(model, structure=structure)
        if boundary:
            current, _ = reconstruct_for_transfer(current, boundary=boundary)
        result = evaluate(current, wavelengths, lines, composition, tables, core_rays,
                          vmic, include_scattering, chemistry, kernel_path, native_convection)
        active_convection = convection
        if native_convection:
            active_convection = [
                4 * math.pi * (reference_radius - row[2])**2 * flux
                for row, flux in zip(current.structure, result[6])
            ]
        luminosity = [r + c for r, c in zip(result[2], active_convection)]
        residual = [(b - a) / target for a, b in zip(luminosity[:count], luminosity[1:count+1])]
        surface = luminosity[0] / target - 1
        deep = [value / target - 1 for value in luminosity[deep_start:]]
        residual.extend(deep)
        if luminosity_target:
            residual.append(surface)
        diagnostics = (max(map(abs, residual[:count])), surface)
        if deep_layers:
            diagnostics = (*diagnostics, max(map(abs, deep)))
        return current, result, residual, diagnostics

    temperatures = [row[3] for row in (*baseline[:count], *baseline[deep_start:])]
    if luminosity_target:
        temperatures.append(0.0)
    current, result, residual, diagnostics = trial(temperatures)
    history = [diagnostics]
    for _ in range(steps):
        if max(map(abs, residual)) < 1e-8:
            break
        columns = []
        for i, value in enumerate(temperatures):
            delta = 0.001 * (value if i < pivot else float(model.header[1].split()[0]))
            if i < pivot and value + delta > chemistry.temperatures[-1]:
                delta = -delta
            perturbed = temperatures.copy()
            perturbed[i] += delta
            shifted = trial(perturbed)[2]
            columns.append([(b - a) / delta for a, b in zip(residual, shifted)])
        direction = _linear_solve([list(row) for row in zip(*columns)], [-x for x in residual])
        # Keep the Newton direction and cap the largest fractional change.
        fractional = max(abs(d) / t for d, t in zip(direction[:pivot], temperatures[:pivot]))
        if luminosity_target:
            fractional = max(fractional, max(abs(direction[-1] * w) / row[3]
                                             for row, w in zip(baseline, mode)))
        scale = min(1.0, 0.03 / fractional) if fractional else 1.0
        for _ in range(8):
            candidate = [t + scale * d for t, d in zip(temperatures, direction)]
            if all(chemistry.temperatures[0] <= t <= chemistry.temperatures[-1]
                   for t in temperature_profile(candidate)):
                new_model, new_result, new_residual, diagnostics = trial(candidate)
                prior = history[-1]
                if (diagnostics[0] <= prior[0] and max(map(abs, new_residual)) < max(map(abs, residual))
                        and (not luminosity_target or abs(diagnostics[1]) <= abs(prior[1]))
                        and (not deep_layers or diagnostics[2] <= prior[2])):
                    temperatures = candidate
                    current, result, residual = new_model, new_result, new_residual
                    history.append(diagnostics)
                    break
            scale *= 0.5
        else:
            break
    return current, result, history


def global_temperature_knots(convective_fractions, upper_knots=5):
    """Resolve each layer of the deep energy-carrying convection zone."""
    count = len(convective_fractions)
    start = count
    while start and convective_fractions[start-1] >= 1e-4:
        start -= 1
    if start == 0:
        return tuple(range(count)), start
    upper_count = min(start, upper_knots)
    if upper_count == 1:
        upper = [0]
    else:
        upper = [round(i*(start-1)/(upper_count-1)) for i in range(upper_count)]
    return (*upper, *range(start, count)), start


def interpolate_corrections(axis, knots, values):
    correction = [0.0] * len(axis)
    for lo, hi, a, b in zip(knots, knots[1:], values, values[1:]):
        for i in range(lo, hi):
            fraction = (axis[i]-axis[lo])/(axis[hi]-axis[lo])
            correction[i] = a + fraction*(b-a)
    correction[knots[-1]] = values[-1]
    return correction


def temperature_direction(columns, residual):
    """Regularize weak temperature controls to keep the trial step local."""
    maximum = max(map(abs, residual))
    weights = [1 + 3*(abs(x)/maximum)**4 for x in residual]
    weighted = [[w*x for w, x in zip(weights, column)] for column in columns]
    weighted_residual = [w*x for w, x in zip(weights, residual)]
    normal = [[math.fsum(a*b for a, b in zip(x, y)) for y in weighted] for x in weighted]
    rhs = [-math.fsum(a*b for a, b in zip(column, weighted_residual)) for column in weighted]
    damping = 1e-8 * max(normal[i][i] for i in range(len(columns)))
    if damping == 0:
        return rhs
    for i, row in enumerate(normal):
        row[i] += damping
    direction = _linear_solve([row.copy() for row in normal], rhs.copy())
    if max(map(abs, direction)) <= 0.03:
        return direction

    # Minimize the same linearized least-squares objective within the
    # per-control step bounds. A weak mode at its bound does not shrink
    # the step at all the other controls.
    direction = [0.0] * len(columns)
    for _ in range(100):
        largest_change = 0.0
        for i, row in enumerate(normal):
            correction = (rhs[i]-math.fsum(a*x for a, x in zip(row, direction)))/row[i]
            updated = max(-0.03, min(0.03, direction[i]+correction))
            largest_change = max(largest_change, abs(updated-direction[i]))
            direction[i] = updated
        if largest_change < 1e-7:
            break
    return direction


def merit_descent(columns, residual, layer_count, peak_weight):
    """Find a fallback direction using the line search's own objective."""
    norm = math.sqrt(math.fsum(x*x for x in residual)*layer_count)
    if not norm:
        return [0.0]*len(columns)
    peak = max(range(len(residual)), key=lambda i: abs(residual[i]))
    derivative = [x/norm for x in residual]
    derivative[peak] += math.copysign(peak_weight, residual[peak])
    gradient = [math.fsum(a*b for a, b in zip(column, derivative)) for column in columns]
    scale = max(map(abs, gradient))
    return [-0.01*x/scale for x in gradient] if scale else gradient


def update_columns(columns, step, previous, current):
    """Update the finite-difference Jacobian to satisfy the accepted secant."""
    square = math.fsum(x*x for x in step)
    remainder = [
        (b - a - math.fsum(dx*column[i] for dx, column in zip(step, columns)))/square
        for i, (a, b) in enumerate(zip(previous, current))
    ]
    return [[x + dx*r for x, r in zip(column, remainder)] for dx, column in zip(step, columns)]


@dataclass(frozen=True)
class TrialTransfer:
    boundary: object
    wavelengths: tuple
    lines: tuple
    mixture: dict
    tables: tuple
    core_rays: int
    vmic: float
    scattering: bool
    chemistry: Thermochemistry
    kernel_path: Path
    line_scale: float
    luminosity_target: float
    flux_offset: tuple = ()

    def __call__(self, profile):
        current, _ = reconstruct_for_transfer(profile, boundary=self.boundary)
        force = getattr(self.boundary, "radiative_acceleration", ())
        result = evaluate(current, self.wavelengths, self.lines, self.mixture, self.tables,
                          self.core_rays, self.vmic, self.scattering, self.chemistry,
                          self.kernel_path, True, force, None, self.line_scale)
        radius = self.boundary.reference_radius
        convection = [
            4*math.pi*(radius-row[2])**2*flux/self.luminosity_target
            for row, flux in zip(current.structure, result[6])
        ]
        offset = self.flux_offset or (0.0,)*len(convection)
        residual = [r/self.luminosity_target + c - 1 + adjustment
                    for r, c, adjustment in zip(result[2], convection, offset)]
        return current, result, residual, convection


def initialize_trial_worker(transfer):
    global _worker_transfer
    _worker_transfer = transfer


def run_trial(profile):
    return _worker_transfer(profile)


def balance_all_layers(model, wavelengths, lines, composition, steps=1,
                       tables=(), core_rays=12, vmic=2.0, include_scattering=True,
                       kernel_path=None, upper_knots=5, boundary=None, progress=None,
                       flux_tolerance=1e-8, jacobian_interval=1, outer_heating_weight=0.0,
                       on_accept=None, adaptive_knots=0, peak_weight=0.25, line_scale=1.0,
                       jacobian_workers=1, flux_offset=()):
    """Fit absolute total luminosity across the entire hydrostatic structure."""
    if steps < 0 or upper_knots < 2:
        raise ValueError("Use a nonnegative step count and at least two upper temperature knots")
    if not math.isfinite(flux_tolerance) or not 0 < flux_tolerance < 1:
        raise ValueError("Flux tolerance must be a fractional error between zero and one")
    if jacobian_interval < 1:
        raise ValueError("Jacobian interval must be at least one")
    if not math.isfinite(outer_heating_weight) or outer_heating_weight < 0:
        raise ValueError("Outer heating weight must be finite and nonnegative")
    if adaptive_knots < 0 or not math.isfinite(peak_weight) or peak_weight < 0:
        raise ValueError("Adaptive temperature controls and peak error weight must be nonnegative")
    if not math.isfinite(line_scale) or not 0 <= line_scale <= 1:
        raise ValueError("Line opacity scale must be in [0, 1]")
    if jacobian_workers < 1:
        raise ValueError("Use at least one worker for the temperature Jacobian")
    if boundary is None:
        chemistry = Thermochemistry()
        boundary, _ = boundary_from_marcs(model, chemistry, composition)
    else:
        chemistry = boundary.chemistry
    radius = boundary.reference_radius
    luminosity_target = 4 * math.pi * radius**2 * float(model.header[2].split()[0])
    baseline = model.structure
    flux_offset = tuple(flux_offset) or (0.0,)*len(baseline)
    if len(flux_offset) != len(baseline) or not all(map(math.isfinite, flux_offset)):
        raise ValueError("The luminosity offset must contain a finite value at every layer")
    axis = tuple(row[0] for row in baseline)
    if any(a >= b for a, b in zip(axis, axis[1:])):
        raise ValueError("Reference log Rosseland depth must increase inward")

    evaluate_trial = TrialTransfer(boundary, wavelengths, lines, composition, tables, core_rays,
                                   vmic, include_scattering, chemistry, kernel_path,
                                   line_scale, luminosity_target, flux_offset)

    def metrics(residual):
        return (math.sqrt(math.fsum(x*x for x in residual)/len(residual)),
                max(map(abs, residual)), residual[0], residual[-1])

    current, result, flux_residual, convection = evaluate_trial(model)
    heating_weights = []
    if outer_heating_weight:
        rosseland = [row[4] for row in result[7]]
        depths = integrate(increments(boundary.column_mass, rosseland), 0.0)
        heating_weights = [outer_heating_weight*math.exp(-tau/0.03) for tau in depths]

    def objective(evaluation, errors):
        return [*errors, *(w*q for w, q in zip(heating_weights, evaluation[8]))] if heating_weights else errors

    def merit(errors):
        return math.sqrt(math.fsum(x*x for x in errors)/len(baseline)) + peak_weight*max(map(abs, errors))

    residual = objective(result, flux_residual)
    knots, convective_start = global_temperature_knots(convection, upper_knots)
    parameters = [0.0] * len(knots)
    history = [metrics(flux_residual)]
    if progress is not None:
        progress(0, history[-1])

    def trial_profile(values):
        correction = interpolate_corrections(axis, knots, values)
        structure = tuple(
            (*row[:3], row[3]*(1+x), *row[4:]) for row, x in zip(baseline, correction)
        )
        if any(not chemistry.temperatures[0] <= row[3] <= chemistry.temperatures[-1]
               for row in structure):
            return None
        return replace(model, structure=structure)

    def finite_difference():
        profiles, deltas = [], []
        for i in range(len(knots)):
            perturb = parameters.copy()
            perturb[i] += 0.001
            profile = trial_profile(perturb)
            if profile is None:
                perturb[i] -= 0.002
                profile = trial_profile(perturb)
            if profile is None:
                raise RuntimeError("Cannot perturb a temperature knot within the LTE table")
            profiles.append(profile)
            deltas.append(perturb[i] - parameters[i])
        if jacobian_workers == 1:
            results = list(map(evaluate_trial, profiles))
        else:
            with ProcessPoolExecutor(max_workers=min(jacobian_workers, len(profiles)),
                                     mp_context=get_context("spawn"), initializer=initialize_trial_worker,
                                     initargs=(evaluate_trial,)) as executor:
                results = list(executor.map(run_trial, profiles))
        columns = []
        for shifted, delta in zip(results, deltas):
            trial_residual = objective(shifted[1], shifted[2])
            columns.append([(b-a)/delta for a, b in zip(residual, trial_residual)])
        return columns

    def line_search(direction):
        if not any(direction):
            return None
        previous_merit = merit(residual)
        scale = 1.0
        for _ in range(6):
            candidate = [p + scale*d for p, d in zip(parameters, direction)]
            profile = trial_profile(candidate)
            if profile is not None:
                proposed = evaluate_trial(profile)
                proposed_residual = objective(proposed[1], proposed[2])
                if merit(proposed_residual) < previous_merit:
                    return candidate, proposed, proposed_residual
            scale *= 0.5
        return None

    columns = None
    for step in range(steps):
        if max(map(abs, residual)) <= flux_tolerance:
            break
        fresh = columns is None or step % jacobian_interval == 0
        if fresh:
            columns = finite_difference()
        accepted = line_search(temperature_direction(columns, residual))
        if accepted is None and not fresh:
            columns = finite_difference()
            accepted = line_search(temperature_direction(columns, residual))
        if accepted is None:
            accepted = line_search(merit_descent(columns, residual, len(baseline), peak_weight))
        if accepted is None:
            break
        candidate, (current, result, flux_residual, _), next_residual = accepted
        if jacobian_interval > 1:
            delta = [b-a for a, b in zip(parameters, candidate)]
            columns = update_columns(columns, delta, residual, next_residual)
        parameters, residual = candidate, next_residual
        history.append(metrics(flux_residual))
        if progress is not None:
            progress(len(history)-1, history[-1])
        if on_accept is not None:
            on_accept(current, result)
        if adaptive_knots and step + 1 < steps:
            worst = max((i for i in range(len(axis)) if i not in knots),
                        key=lambda i: abs(flux_residual[i]), default=None)
            if worst is not None and abs(flux_residual[worst]) > max(flux_tolerance, 0.7*history[-1][1]):
                correction = interpolate_corrections(axis, knots, parameters)
                knots = tuple(sorted((*knots, worst)))
                parameters = [correction[i] for i in knots]
                columns = None
                adaptive_knots -= 1
    return current, result, history, knots, convective_start


def report_balance_step(step, metrics):
    rms, maximum, surface, bottom = metrics
    print(f"Step {step}: luminosity error RMS {rms:.3%}, maximum {maximum:.3%}, "
          f"surface {surface:+.3%}, bottom {bottom:+.3%}", file=sys.stderr, flush=True)


def total_flux_budget(model, result, native_convection):
    flux = result[0]
    budget = list(reference_flux_budget(model, flux))
    convection = result[6] if native_convection else [row[0] for row in budget]
    return [(c, f+c, target, (f+c-target)/target)
            for c, f, (_, _, target, _) in zip(convection, flux, budget)]


def write_balance_csv(stream, initial, model, result, native_convection):
    flux, heating, luminosity, shells, _, _ = result[:6]
    budget = total_flux_budget(model, result, native_convection)
    writer = csv.writer(stream)
    columns = ("layer", "initial_temperature_K", "temperature_K", "column_mass_g_cm2", "radiative_flux_erg_cm2_s",
                     "pointwise_heating_erg_g_s", "radiative_luminosity_erg_s",
                     "shell_to_next_heating_erg_g_s", "shell_to_next_heating_erg_s",
                     "convective_flux_erg_cm2_s" if native_convection else "marcs_convective_flux_erg_cm2_s",
                     "total_flux_erg_cm2_s",
                     "target_flux_erg_cm2_s", "relative_total_flux_error",
                     "gas_pressure_dyn_cm2", "electron_pressure_dyn_cm2", "density_g_cm3", "radius_cm")
    extra = [()] * len(flux)
    if native_convection:
        columns += ("rosseland_mean_cm2_g", "relative_net_heating", "radiative_acceleration_cm_s2")
        extra = [(row[4], heating, force) for row, heating, force in zip(*result[7:10])]
    writer.writerow(columns)
    reference_radius = float(model.header[7].split()[0])
    for i, (old, row, thermo, f, q, lum, shell, diagnostic, native) in enumerate(
        zip(initial.structure, model.structure, model.thermodynamics, flux, heating, luminosity,
            [*shells, ("", "")], budget, extra), 1
    ):
        writer.writerow((i, old[3], row[3], thermo[-1], f, q, lum, *shell, *diagnostic,
                         row[5], row[4], thermo[2], reference_radius - row[2], *native))
    return budget


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
    parser.add_argument("--luminosity-target", action="store_true",
                        help="Also fit the absolute MARCS luminosity using a photospheric temperature mode")
    parser.add_argument("--verify-full", action="store_true",
                        help="Reevaluate the accepted profile on the full input wavelength grid")
    parser.add_argument("--kernel", type=Path,
                        help="Optional shared library compiled from spherical_kernel.c")
    parser.add_argument("--hydrostatic", action="store_true",
                        help="Solve hydrostatic pressure and radius before transfer and every temperature trial")
    parser.add_argument("--native-convection", action="store_true",
                        help="Compute convection from native Rosseland opacity and LTE thermodynamics")
    parser.add_argument("--adjust-deep", type=int, default=0,
                        help="Also fit the absolute flux at this many bottom layers")
    parser.add_argument("--global-balance", action="store_true",
                        help="Fit total luminosity at every layer with native convection and hydrostatic balance")
    parser.add_argument("--upper-knots", type=int, default=5,
                        help="Number of temperature control points above the deep convective region")
    args = parser.parse_args()
    if not math.isfinite(args.vmic) or args.vmic < 0 or args.adjust_layers < 0 or args.steps < 1:
        parser.error("Use nonnegative microturbulence and layer count, and at least one step")
    if args.luminosity_target and not (args.adjust_layers or args.global_balance):
        parser.error("--luminosity-target requires --adjust-layers or --global-balance")
    if args.verify_full and not (args.adjust_layers or args.global_balance):
        parser.error("--verify-full requires --adjust-layers or --global-balance")
    if args.adjust_deep and (not args.native_convection or not args.adjust_layers):
        parser.error("--adjust-deep requires --native-convection and --adjust-layers")
    if args.global_balance and (args.adjust_layers or args.adjust_deep):
        parser.error("--global-balance selects its own temperature controls")
    if args.global_balance:
        args.hydrostatic = args.native_convection = True
    try:
        model = read_document(args.model)
        reference = read_wavelengths(args.wavelengths)
        wavelengths = nested_wavelengths(reference, args.stride)
        lines = read_lines(args.lines)
        with Path(args.composition).open(encoding="utf-8") as stream:
            composition = json.load(stream)
        tables = [MolecularCrossSection(path) for path in args.molecular_cross_section]
        initial = model
        if args.hydrostatic and not (args.adjust_layers or args.global_balance):
            model, solution = reconstruct_for_transfer(model, mixture=composition)
            print(f"Hydrostatic iterations: {solution.iterations}; "
                  f"relative cell pressure residual: {solution.pressure_residual:.3e}", file=sys.stderr)
        if args.global_balance:
            model, result, history, knots, start = balance_all_layers(
                model, wavelengths, lines, composition, args.steps, tables, args.core_rays,
                args.vmic, not args.absorption_only, args.kernel, args.upper_knots,
                progress=report_balance_step,
            )
            print(f"Temperature controls: {', '.join(str(i+1) for i in knots)}", file=sys.stderr)
            if start < len(model.structure):
                print(f"Deep convective region: layers {start+1}–{len(model.structure)}",
                      file=sys.stderr)
        elif args.adjust_layers:
            model, result, history = correct_temperatures(
                model, wavelengths, lines, composition, args.adjust_layers, args.steps,
                tables, args.core_rays, args.vmic, include_scattering=not args.absorption_only,
                luminosity_target=args.luminosity_target, kernel_path=args.kernel, hydrostatic=args.hydrostatic,
                native_convection=args.native_convection, deep_layers=args.adjust_deep,
            )
            for step, (residual, surface, *deep) in enumerate(history):
                deep_text = f", max deep luminosity error {deep[0]:.3%}" if deep else ""
                print(f"Step {step}: max outer-shell luminosity residual {residual:.4e}, "
                      f"surface total-flux error {surface:+.3%}{deep_text}", file=sys.stderr)
        else:
            result = evaluate(model, wavelengths, lines, composition, tables,
                              args.core_rays, args.vmic, include_scattering=not args.absorption_only,
                              kernel_path=args.kernel, native_convection=args.native_convection)
        if args.verify_full and wavelengths != reference:
            print(f"Verifying the corrected profile on {len(reference):,} wavelengths", file=sys.stderr)
            result = evaluate(model, reference, lines, composition, tables,
                              args.core_rays, args.vmic, include_scattering=not args.absorption_only,
                              kernel_path=args.kernel, native_convection=args.native_convection)
        budget = write_balance_csv(sys.stdout, initial, model, result, args.native_convection)
    except (IndexError, OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    iterations, mismatch = result[4:6]
    print(f"Surface total-flux error: {budget[0][-1]:+.3%}", file=sys.stderr)
    print(f"Maximum scattering iterations per wavelength: {iterations}", file=sys.stderr)
    print(f"Relative mismatch from trapezoidal integration of pointwise heating: {mismatch:+.3e}", file=sys.stderr)


if __name__ == "__main__":
    main()
