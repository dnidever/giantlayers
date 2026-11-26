"""Compare plane-parallel and spherical LTE surface transfer on a MARCS structure."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

from .gray_atmosphere import propagate
from .lte_opacity import layer_opacities, prepare_layers, read_lines
from .marcs_writer import read_document
from .molecular_cross_section import MolecularCrossSection
from .native_transfer import read_wavelengths, spectral_quadrature
from .odf_transfer import planck_nu


def impact_parameters(radii, core_rays):
    """Return squared shell radii and ordered ray intersections."""
    squared = [radius**2 for radius in radii]
    impacts = [(squared[-1] * i / core_rays, len(radii) - 1, True)
               for i in range(core_rays + 1)]
    # Include both limits at the core's edge. A ray just inside receives
    # bottom-boundary emission, while a ray just outside traverses the shell.
    impacts.extend((squared[i], i, False) for i in range(len(radii) - 1, -1, -1))
    return squared, impacts

def disk_rays(radii, core_rays):
    """Build rays and disk weights on a decreasing radial grid.

    Rays in the apparent core are uniform in impact parameter squared.
    Shell-intersecting rays turn at the atmospheric grid radii.
    """
    if (len(radii) < 2 or not all(math.isfinite(r) and r > 0 for r in radii)
            or any(a <= b for a, b in zip(radii, radii[1:])) or core_rays < 2):
        raise ValueError("Use decreasing positive radii and at least two core intervals")

    squared, impacts = impact_parameters(radii, core_rays)
    rays = []
    for i, (impact, last, core) in enumerate(impacts):
        left = impact - impacts[i - 1][0] if i else 0
        right = impacts[i + 1][0] - impact if i + 1 < len(impacts) else 0
        weight = math.pi * (left + right) / (2 * squared[0])
        chord = [math.sqrt(max(0, r2 - impact)) for r2 in squared[:last + 1]]
        # Rationalized chord differences avoid cancellation far from the limb.
        scale = [(outer + inner) / (a + b)
                 for outer, inner, a, b in zip(radii, radii[1:], chord, chord[1:])]
        mu_surface = chord[0] / radii[0]
        mu_bottom = chord[-1] / radii[last]
        rays.append((core, scale, mu_surface, mu_bottom, weight))
    return rays


def _surface_intensity(source, optical_steps, core, radial_bottom_step, mu_bottom):
    """Integrate the formal solution with zero incident intensity."""
    last = len(optical_steps)
    if core:
        intensity = source[last]
        if sum(optical_steps) >= 10 and radial_bottom_step > 0:
            intensity += mu_bottom * (source[last] - source[last - 1]) / radial_bottom_step
    else:
        intensity = 0.0
        for i, step in enumerate(optical_steps):
            intensity = propagate(intensity, source[i], source[i + 1], step)

    for i in range(last - 1, -1, -1):
        intensity = propagate(intensity, source[i + 1], source[i], optical_steps[i])
    return intensity


def surface_fluxes(source, radial_steps, rays):
    """Return spherical and planar surface flux densities with matched disk weights."""
    spherical = plane = 0.0
    for core, chord_scale, mu, mu_bottom, weight in rays:
        if not mu:
            continue
        ray_depth = [step * scale for step, scale in zip(radial_steps, chord_scale)]
        spherical += weight * _surface_intensity(source, ray_depth, core, radial_steps[-1], mu_bottom)
        plane += weight * _surface_intensity(source, [step / mu for step in radial_steps],
                                          True, radial_steps[-1], mu)
    # The outer tangent has zero intensity in spherical geometry. In the
    # planar limit it sees the surface source function for positive opacity.
    plane += rays[-1][-1] * source[0]
    return spherical, plane


def sample_fluxes(model, wavelengths, lines, composition, tables=(), core_rays=12, vmic=2.0):
    """Yield absorption-only surface flux densities and per-sample weights."""
    if float(model.header[5].split()[0]) <= 0:
        raise ValueError("A spherical MARCS structure with nonzero stellar mass is required")
    reference_radius = float(model.header[7].split()[0])
    radii = [reference_radius - row[2] for row in model.structure]
    rays = disk_rays(radii, core_rays)
    radial_widths = [a - b for a, b in zip(radii, radii[1:])]
    frequency, weights = spectral_quadrature(wavelengths, "wavelength")
    layers = list(prepare_layers(model, lines, composition, vmic))
    for start in range(0, len(wavelengths), 256):
        end = start + 256
        sample = wavelengths[start:end]
        grids = [(table, table.on_grid(sample)) for table in tables]
        opacities = [list(layer_opacities(layer, sample, grids)) for layer in layers]
        for columns, wavelength, nu, weight in zip(
            zip(*opacities), sample, frequency[start:end], weights[start:end]
        ):
            extinction = [(line + continuum) * layer[0].mass_density
                          for (line, continuum, _), layer in zip(columns, layers)]
            steps = [width * (a + b) / 2 for width, a, b in
                     zip(radial_widths, extinction, extinction[1:])]
            source = [planck_nu(nu, state[0].temperature) for state in layers]
            sphere, plane = surface_fluxes(source, steps, rays)
            yield wavelength, sphere, plane, weight


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model")
    parser.add_argument("lines")
    parser.add_argument("composition")
    parser.add_argument("wavelengths")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--core-rays", type=int, default=12)
    parser.add_argument("--vmic", type=float, default=2)
    parser.add_argument("--molecular-cross-section", action="append", default=[])
    args = parser.parse_args()
    if args.stride < 1 or not math.isfinite(args.vmic) or args.vmic < 0:
        parser.error("Use a positive stride and finite nonnegative microturbulence")
    try:
        model = read_document(args.model)
        lines = read_lines(args.lines)
        with Path(args.composition).open(encoding="utf-8") as stream:
            composition = json.load(stream)
        axis = read_wavelengths(args.wavelengths)
        wavelengths = tuple(sorted(set(axis[::args.stride]) | {axis[-1]}))
        tables = [MolecularCrossSection(path) for path in args.molecular_cross_section]
        rows = list(sample_fluxes(model, wavelengths, lines, composition, tables,
                                 args.core_rays, args.vmic))
    except (IndexError, OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("wavelength_vacuum_A", "spherical_flux_per_Hz", "planar_flux_per_Hz",
                     "frequency_weight_Hz"))
    writer.writerows(rows)
    sphere = math.fsum(row[1] * row[3] for row in rows)
    plane = math.fsum(row[2] * row[3] for row in rows)
    reference_radius = float(model.header[7].split()[0])
    outer_radius = reference_radius - model.structure[0][2]
    target = float(model.header[2].split()[0])
    spherical_ratio = sphere * (outer_radius / reference_radius)**2 / target
    local_plane_ratio = plane / (target * (reference_radius / outer_radius)**2)
    print(f"Integrated spherical/reference and planar/local-target flux ratios: "
          f"{spherical_ratio:.8f}, {local_plane_ratio:.8f}", file=sys.stderr)


if __name__ == "__main__":
    main()
