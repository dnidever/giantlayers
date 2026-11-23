"""Sampled LTE transfer and outer-layer temperature corrections with native opacities."""

import argparse
import csv
from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
import sys

from lte_chemistry import Thermochemistry
from lte_opacity import C_ANGSTROM_S, H_EV_S, opacities, read_lines
from marcs_writer import read_document
from molecular_cross_section import MolecularCrossSection
from odf_transfer import planck_nu, scattering_ray_flux
from hminus import K_EV


@dataclass(frozen=True)
class EnergyBalance:
    flux: tuple
    heating: tuple
    cell_heating: tuple
    emission: tuple
    emission_derivative: tuple

    def residual(self, count):
        return max(abs(q) / e for q, e in zip(self.heating[:count], self.emission[:count]))

    @property
    def closure(self):
        return math.fsum(self.cell_heating) - (self.flux[-1] - self.flux[0])


def reference_flux_budget(model, radiative_flux):
    """Compare against MARCS total flux with its tabulated convection held fixed.

    The reference target scales with radius in spherical MARCS models.
    GiantLayers' transfer remains plane-parallel in this diagnostic.
    """
    surface_flux = float(model.header[2].split()[0])
    spherical = float(model.header[5].split()[0]) > 0
    reference_radius = float(model.header[7].split()[0])
    if len(radiative_flux) != len(model.structure):
        raise ValueError("Radiative flux must have one value per MARCS layer")
    for flux, state, thermo in zip(radiative_flux, model.structure, model.thermodynamics):
        radius = reference_radius - state[2] if spherical else reference_radius
        target = surface_flux * (reference_radius / radius)**2
        convection = target * thermo[5]
        total = flux + convection
        yield convection, total, target, (total - target) / target


def spectral_quadrature(wavelengths):
    if (len(wavelengths) < 2 or not all(math.isfinite(w) and w > 0 for w in wavelengths)
            or any(a >= b for a, b in zip(wavelengths, wavelengths[1:]))):
        raise ValueError("Supply at least two increasing positive vacuum wavelengths")
    frequencies = tuple(C_ANGSTROM_S / w for w in wavelengths)
    intervals = [a - b for a, b in zip(frequencies, frequencies[1:])]
    weights = (intervals[0] / 2,
               *((a + b) / 2 for a, b in zip(intervals, intervals[1:])),
               intervals[-1] / 2)
    return frequencies, weights


def evaluate(model, wavelengths, lines, composition, vmic=2.0, chemistry=None,
             cross_sections=()):
    """Integrate sampled flux and net radiative heating over the supplied grid."""
    frequency, weights = spectral_quadrature(wavelengths)
    mass = [row[-1] for row in model.thermodynamics]
    temperature = [row[3] for row in model.structure]
    if any(a >= b for a, b in zip(mass, mass[1:])):
        raise ValueError("Column mass must increase inward")
    n_layers = len(mass)
    rows = iter(opacities(model, wavelengths, lines, composition, vmic, chemistry, cross_sections))
    # Native opacity yields each depth's wavelengths consecutively.
    thermal_layers, scattering_layers = [], []
    for _ in mass:
        layer = [next(rows) for _ in wavelengths]
        thermal_layers.append([row[2] + row[3] for row in layer])
        scattering_layers.append([row[4] for row in layer])

    totals = [[0.0] * n_layers for _ in range(5)]
    for j, (nu, weight) in enumerate(zip(frequency, weights)):
        absorption = [row[j] for row in thermal_layers]
        scattering = [row[j] for row in scattering_layers]
        source = [planck_nu(nu, t) for t in temperature]
        flux, heating, cells = scattering_ray_flux(
            mass, source, absorption, scattering, diagnostics=True
        )
        for i, (f, q, cell, t, b, k) in enumerate(
            zip(flux, heating, cells, temperature, source, absorption)
        ):
            emission = 4 * math.pi * k * b
            x = H_EV_S * nu / (K_EV * t)
            derivative = emission * x / (-t * math.expm1(-x))
            for profile, value in zip(totals, (f, q, cell, emission, derivative)):
                profile[i] += weight * value
    return EnergyBalance(*(tuple(profile) for profile in totals))


def correct_temperatures(model, wavelengths, lines, composition, count, steps=1,
                         vmic=2.0, flux_guard=False, cross_sections=()):
    """Attempt damped corrections at fixed pressure and column mass.

    The deeper layers act as a fixed boundary. Recompute chemistry, opacity and
    transfer for each trial. Keep only steps that reduce the outer-layer residual.
    """
    if not 1 <= count < len(model.structure) or steps < 1:
        raise ValueError("Select a positive number of outer layers and keep deeper layers fixed")
    chemistry = Thermochemistry()
    balance = evaluate(model, wavelengths, lines, composition, vmic, chemistry, cross_sections)
    history = [balance.residual(count)]
    surface_error = abs(next(reference_flux_budget(model, balance.flux))[3])
    for _ in range(steps):
        if history[-1] < 1e-8:
            break
        changes = []
        for row, heating, derivative in zip(
            model.structure[:count], balance.heating, balance.emission_derivative
        ):
            limit = 0.05 * row[3]
            changes.append(max(-limit, min(limit, heating / derivative)))
        for backtrack in range(8):
            scale = 2 ** -backtrack
            structure = tuple(
                (*row[:3], row[3] + scale * delta, *row[4:])
                for row, delta in zip(model.structure[:count], changes)
            ) + model.structure[count:]
            trial_model = replace(model, structure=structure)
            if not all(chemistry.temperatures[0] <= row[3] <= chemistry.temperatures[-1]
                       for row in structure[:count]):
                continue
            trial = evaluate(trial_model, wavelengths, lines, composition, vmic, chemistry, cross_sections)
            residual = trial.residual(count)
            trial_surface_error = abs(next(reference_flux_budget(trial_model, trial.flux))[3])
            if residual < history[-1] and (not flux_guard or trial_surface_error <= surface_error):
                model, balance = trial_model, trial
                history.append(residual)
                surface_error = trial_surface_error
                break
        else:
            break
    return model, balance, history


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", help="Reference MARCS .mod or .mod.gz")
    parser.add_argument("lines", help="Native transition CSV")
    parser.add_argument("composition", help="JSON abundance pattern")
    parser.add_argument("wavelengths", help="Text file with ascending vacuum wavelengths in Å")
    parser.add_argument("output", type=Path, help="New energy-balance CSV")
    parser.add_argument("--vmic", type=float, default=2)
    parser.add_argument("--adjust-layers", type=int, default=0)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--flux-guard", action="store_true",
                        help="Reject temperature trials that worsen the surface total-flux error")
    parser.add_argument("--molecular-cross-section", action="append", default=[], metavar="FILE",
                        help="Add a molecular cross-section JSON or JSON.gz table")
    args = parser.parse_args()
    if not math.isfinite(args.vmic) or args.vmic < 0 or args.adjust_layers < 0 or args.steps < 1:
        parser.error("Use nonnegative microturbulence and layer count, and at least one step")
    try:
        model = read_document(args.model)
        lines = read_lines(args.lines)
        with open(args.composition, encoding="utf-8") as stream:
            composition = json.load(stream)
        wavelengths = tuple(map(float, Path(args.wavelengths).read_text().split()))
        tables = [MolecularCrossSection(path) for path in args.molecular_cross_section]
        if args.adjust_layers:
            result, balance, history = correct_temperatures(
                model, wavelengths, lines, composition, args.adjust_layers, args.steps, args.vmic,
                flux_guard=args.flux_guard, cross_sections=tables,
            )
            print("Outer-layer radiative residual: " + " -> ".join(f"{x:.6g}" for x in history),
                  file=sys.stderr)
            print(f"Accepted temperature corrections: {len(history) - 1}/{args.steps}", file=sys.stderr)
        else:
            result = model
            balance = evaluate(model, wavelengths, lines, composition, args.vmic, cross_sections=tables)
        budget = list(reference_flux_budget(model, balance.flux))
        with args.output.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("layer", "column_mass_g_cm2", "initial_temperature_K", "temperature_K",
                             "radiative_flux_erg_cm2_s", "heating_erg_g_s", "cell_heating_erg_cm2_s",
                             "marcs_convective_flux_erg_cm2_s", "total_flux_erg_cm2_s",
                             "target_flux_erg_cm2_s", "relative_total_flux_error"))
            for i, (thermo, initial, updated, flux, heating, cell, diagnostic) in enumerate(
                zip(model.thermodynamics, model.structure, result.structure,
                    balance.flux, balance.heating, balance.cell_heating, budget), 1
            ):
                writer.writerow((i, thermo[-1], initial[3], updated[3], flux, heating, cell, *diagnostic))
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    print(f"Surface total-flux error: {budget[0][3]:+.3%}", file=sys.stderr)
    scale = max(map(abs, balance.flux))
    print(f"Relative heating/flux closure error: {abs(balance.closure) / scale:.3e}", file=sys.stderr)


if __name__ == "__main__":
    main()
