"""Solve spherical hydrostatic structure at a prescribed temperature profile."""

import argparse
import csv
from dataclasses import dataclass, replace
import math
import sys

from hminus import K_B
from lte_chemistry import Thermochemistry
from lte_opacity import composition
from marcs_writer import read_document

G = 6.67430e-8
SOLAR_MASS = 1.98847e33


@dataclass(frozen=True)
class HydrostaticProfile:
    pressure: tuple
    radius: tuple
    density: tuple
    electron_density: tuple
    iterations: int
    pressure_residual: float
    radius_residual: float


def increments(axis, values):
    """Integrate over successive cells using the trapezoidal rule."""
    return [(b - a) * (x + y) / 2 for a, b, x, y in
            zip(axis, axis[1:], values, values[1:])]


def integrate(cell_values, initial):
    profile = [initial]
    for change in cell_values:
        profile.append(profile[-1] + change)
    return profile


def solve(column_mass, temperature, abundances, stellar_mass, gravity, surface_pressure,
          reference_index, chemistry=None, tolerance=1e-8, max_iterations=60):
    """Solve dP/dm = GM/r² and dr/dm = -1/rho on a radial column-mass grid.

    Use cgs units, with stellar mass in grams. The column mass increases
    inward and is the radial integral of density. The layer at
    reference_index is anchored to sqrt(GM/gravity).
    """
    n = len(column_mass)
    if (n < 2 or len(temperature) != n or not 0 <= reference_index < n
            or not all(math.isfinite(x) and x >= 0 for x in column_mass)
            or any(a >= b for a, b in zip(column_mass, column_mass[1:]))
            or not all(math.isfinite(x) and x > 0 for x in
                       (stellar_mass, gravity, surface_pressure, tolerance))
            or max_iterations < 1):
        raise ValueError("Use increasing nonnegative column mass, a valid reference layer, and positive boundaries")
    chemistry = chemistry or Thermochemistry()
    gm = G * stellar_mass
    reference_radius = math.sqrt(gm / gravity)
    radii = [reference_radius] * n
    for iteration in range(1, max_iterations + 1):
        local_g = [gm / r**2 for r in radii]
        pressure = integrate(increments(column_mass, local_g), surface_pressure)
        gas = [chemistry.solve(t, p, abundances) for t, p in zip(temperature, pressure)]
        density = [state.mass_density for state in gas]
        thickness = integrate(increments(column_mass, [1 / rho for rho in density]), 0.0)
        anchor = thickness[reference_index]
        updated_radius = [reference_radius + (anchor - depth) for depth in thickness]
        if any(r <= 0 for r in updated_radius):
            raise RuntimeError("Hydrostatic update reached a nonpositive radius")
        new_gravity = [gm / r**2 for r in updated_radius]
        expected = increments(column_mass, new_gravity)
        actual = [b-a for a, b in zip(pressure, pressure[1:])]
        pressure_error = max(abs(a-b)/b for a, b in zip(actual, expected))
        radius_error = max(abs(a-b)/reference_radius for a,b in zip(radii, updated_radius))
        if max(pressure_error, radius_error) <= tolerance:
            electrons = tuple(state.electron_density for state in gas)
            return HydrostaticProfile(tuple(pressure), tuple(updated_radius), tuple(density),
                                      electrons, iteration, pressure_error, radius_error)
        radii = updated_radius
    raise RuntimeError(f"Hydrostatic solution did not converge in {max_iterations} iterations")


def marcs_radial_grid(model):
    """Build benchmark radial column mass from MARCS radii and densities."""
    reference_radius = float(model.header[7].split()[0])
    radii = [reference_radius - row[2] for row in model.structure]
    density = [row[2] for row in model.thermodynamics]

    def logarithmic_mean(a, b):
        return a if a == b else (b-a) / (math.log(b)-math.log(a))

    # Integrate density assuming log rho is linear between published layers.
    # This supplies a comparison grid independently of the solver's radial
    # integration, which uses a trapezoid in inverse density.
    cells = [(a-b) * logarithmic_mean(x, y)
             for a, b, x, y in zip(radii, radii[1:], density, density[1:])]
    return integrate(cells, 0.0), radii


def solve_marcs_benchmark(model, chemistry=None, mixture=None):
    """Reconstruct hydrostatic conditions on the reference radial mass grid."""
    mass = float(model.header[5].split()[0]) * SOLAR_MASS
    gravity = float(model.header[3].split()[0])
    anchor = min(range(len(model.structure)), key=lambda i: abs(model.structure[i][0]))
    if mass <= 0 or not math.isclose(model.structure[anchor][0], 0.0, abs_tol=1e-8):
        raise ValueError("Use a spherical MARCS model with a layer at log(tau_Ross) = 0")
    column_mass, original_radius = marcs_radial_grid(model)
    chemistry = chemistry or Thermochemistry()
    abundances = {el: model.abundances[atom["atomic_number"] - 1]
                  for el, atom in chemistry.atoms.items()}
    if mixture is not None:
        abundances = composition(abundances, mixture)
    solution = solve(column_mass, [row[3] for row in model.structure], abundances,
                     mass, gravity, model.structure[0][5], anchor, chemistry)
    return solution, column_mass, original_radius, anchor


def reconstruct_for_transfer(model, chemistry=None, mixture=None):
    """Build a transfer input with a fixed temperature and new hydrostatic structure.

    Rosseland depth labels and the convective fraction are reference data.
    This in-memory record is intended for transfer, since the chemistry
    and convection tables in a complete model would also require updating.
    """
    solution, mass, _, anchor = solve_marcs_benchmark(model, chemistry, mixture)
    reference_radius = solution.radius[anchor]
    header = list(model.header)
    _, label = header[7].split(maxsplit=1)
    header[7] = f"{reference_radius:.17e} {label}"
    structure = tuple(
        (*row[:2], reference_radius - r, row[3], K_B * row[3] * ne, p, *row[6:])
        for row, r, ne, p in zip(model.structure, solution.radius, solution.electron_density,
                                 solution.pressure)
    )
    thermo = tuple(
        (*row[:2], rho, *row[3:6], column)
        for row, rho, column in zip(model.thermodynamics, solution.density, mass)
    )
    return replace(model, header=tuple(header), structure=structure, thermodynamics=thermo), solution


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", help="Spherical MARCS reference .mod or .mod.gz")
    args = parser.parse_args()
    try:
        model = read_document(args.model)
        result, column_mass, reference_radius, _ = solve_marcs_benchmark(model)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("layer", "log_tau_ross", "radial_column_mass_g_cm2", "temperature_K",
                     "gas_pressure_dyn_cm2", "marcs_gas_pressure_dyn_cm2",
                     "radius_cm", "marcs_radius_cm", "density_g_cm3", "marcs_density_g_cm3"))
    for i, (m, state, p, r, r_ref, rho, thermo) in enumerate(zip(
        column_mass, model.structure, result.pressure, result.radius, reference_radius,
        result.density, model.thermodynamics,
    ), 1):
        writer.writerow((i, state[0], m, state[3], p, state[5], r, r_ref, rho, thermo[2]))
    print(f"Iterations: {result.iterations}, max relative cell hydrostatic residual: "
          f"{result.pressure_residual:.3e}", file=sys.stderr)


if __name__ == "__main__":
    main()
