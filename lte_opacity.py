"""Native LTE line and continuum mass opacity on a fixed MARCS structure."""

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path

from hminus import K_EV, mass_opacity as hminus_opacity
from lte_chemistry import AMU, Thermochemistry
from marcs_writer import read_document

C_ANGSTROM_S = 2.99792458e18
LIGHT = C_ANGSTROM_S * 1e-8
H_EV_S = 4.135667696e-15
LINE_STRENGTH = 0.02654009  # pi e²/(m_e c), cm² Hz in cgs
THOMSON = 6.6524587321e-25  # cm²


@dataclass(frozen=True)
class Line:
    wavelength: float
    species: str
    lower_energy: float
    log_gf: float
    damping_s: float = 0.0


def read_lines(path):
    """Read vacuum Å, lower excitation in eV and optional damping in s⁻¹."""
    with open(path, newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"wavelength_vacuum_A", "species", "lower_energy_eV", "log_gf"}
        if not required <= set(reader.fieldnames or ()):
            raise ValueError(f"Line list requires: {', '.join(sorted(required))}")
        result = []
        for row in reader:
            line = Line(float(row["wavelength_vacuum_A"]), row["species"].strip(),
                        float(row["lower_energy_eV"]), float(row["log_gf"]),
                        float(row.get("damping_s") or 0))
            if (not all(math.isfinite(x) for x in (line.wavelength, line.lower_energy,
                                                   line.log_gf, line.damping_s))
                    or line.wavelength <= 0 or min(line.lower_energy, line.damping_s) < 0):
                raise ValueError(f"Invalid transition on line {reader.line_num}")
            result.append(line)
    return result


def composition(reference, config):
    """Apply a metal offset and absolute log epsilon overrides."""
    if not config.keys() <= {"metallicity_offset", "log_epsilon"}:
        raise ValueError("Unknown composition setting")
    offset = config.get("metallicity_offset", 0.0)
    changes = config.get("log_epsilon", {})
    if (not isinstance(changes, dict) or not changes.keys() <= reference.keys()
            or not math.isfinite(offset) or changes.get("H", 12) != 12):
        raise ValueError("Use known element symbols, a finite metal offset, and A(H) = 12")
    return {e: changes.get(e, value + (offset if e not in ("H", "He") else 0))
            for e, value in reference.items()}


def continuum_opacity(wavelength, gas, ion_density):
    """Return LTE absorption and coherent electron scattering in cm²/g."""
    frequency = C_ANGSTROM_S / wavelength
    kt_ev = K_EV * gas.temperature
    photon_ev = H_EV_S * frequency
    stimulated = -math.expm1(-photon_ev / kt_ev)
    neutral_h = gas.number_densities["H"] / gas.mass_density
    bound_free, free_free = hminus_opacity(
        wavelength / 10, gas.temperature, gas.electron_density, neutral_h
    )
    # Hydrogenic H I bound-free estimate, with unit Gaunt factor.
    hydrogen_bf = 0.0
    rydberg = 13.5984
    for n in range(1, 9):
        if photon_ev >= rydberg / n**2:
            population = neutral_h * 2 * n**2 / gas.partition_functions["H"]
            population *= math.exp(-rydberg * (1 - 1 / n**2) / kt_ev)
            hydrogen_bf += population * (6.30e-18 / n**5) * (rydberg / photon_ev)**3
    # Thermal electron–ion free-free, singly charged ions and unit Gaunt factor.
    ion_ff = (3.692e8 * gas.electron_density * ion_density
              / (gas.mass_density * math.sqrt(gas.temperature) * frequency**3))
    scattering = THOMSON * gas.electron_density / gas.mass_density
    return bound_free + free_free + stimulated * (hydrogen_bf + ion_ff), scattering


def profile_parameters(line, gas, vmic):
    center = C_ANGSTROM_S / line.wavelength
    thermal = 2 * K_EV * gas.temperature * 1.602176634e-12 / (gas.masses_amu[line.species] * AMU)
    doppler = center / LIGHT * math.sqrt(thermal + (vmic * 1e5)**2)
    gaussian_width = 2 * math.sqrt(math.log(2)) * doppler
    lorentz_width = line.damping_s / (2 * math.pi)
    # Thompson–Cox–Hastings normalized pseudo-Voigt approximation.
    g, l = gaussian_width, lorentz_width
    width = (g**5 + 2.69269*g**4*l + 2.42843*g**3*l*l + 4.47163*g*g*l**3
             + 0.07842*g*l**4 + l**5)**0.2
    ratio = l / width
    mixing = min(1.0, 1.36603*ratio - 0.47719*ratio**2 + 0.11116*ratio**3)
    kt = K_EV * gas.temperature
    strength = (LINE_STRENGTH * 10**line.log_gf * gas.number_densities[line.species]
                / (gas.mass_density * gas.partition_functions[line.species])
                * math.exp(-line.lower_energy / kt) * -math.expm1(-H_EV_S * center / kt))
    return center, width, mixing, strength


def line_opacity(frequency, transitions):
    """Sum the sampled line absorption in cm²/g."""
    total = 0.0
    for center, width, mixing, strength in transitions:
        displacement = 2 * (frequency - center) / width
        gaussian = math.sqrt(math.log(2) / math.pi) * math.exp(-math.log(2) * displacement**2)
        lorentz = 1 / (math.pi * (1 + displacement**2))
        total += 2 * strength / width * ((1 - mixing) * gaussian + mixing * lorentz)
    return total


def opacities(model, wavelengths, lines, config, vmic=2.0, chemistry=None):
    """Yield opacity rows, computing a fresh gas state at each model depth."""
    chemistry = chemistry or Thermochemistry()
    reference = {element: model.abundances[atom["atomic_number"] - 1]
                 for element, atom in chemistry.atoms.items()}
    abundances = composition(reference, config)
    for index, row in enumerate(model.structure, start=1):
        gas = chemistry.solve(row[3], row[5], abundances)
        unknown = {line.species for line in lines} - gas.partition_functions.keys()
        if unknown:
            raise ValueError(f"Unsupported line species: {', '.join(sorted(unknown))}")
        ions = math.fsum(n for species, n in gas.number_densities.items() if species.endswith("+"))
        transitions = [profile_parameters(line, gas, vmic) for line in lines]
        for wavelength in wavelengths:
            absorption, scattering = continuum_opacity(wavelength, gas, ions)
            line = line_opacity(C_ANGSTROM_S / wavelength, transitions)
            yield (index, wavelength, line, absorption, scattering,
                   gas.electron_density, gas.mass_density)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", help="Reference MARCS .mod or .mod.gz")
    parser.add_argument("lines", help="Transition CSV")
    parser.add_argument("composition", help="JSON abundance pattern")
    parser.add_argument("start", type=float, help="First vacuum wavelength in Å")
    parser.add_argument("stop", type=float, help="Upper wavelength bound in Å")
    parser.add_argument("step", type=float, help="Wavelength spacing in Å")
    parser.add_argument("output", type=Path, help="New opacity CSV")
    parser.add_argument("--vmic", type=float, default=2, help="Microturbulence in km/s")
    args = parser.parse_args()
    start, stop, step = args.start, args.stop, args.step
    if (not all(math.isfinite(v) for v in (start, stop, step, args.vmic))
            or not (0 < start < stop and 0 < step <= stop-start and args.vmic >= 0)):
        parser.error("Use a positive increasing wavelength grid and nonnegative microturbulence")
    wavelengths = [start + i*step for i in range(int((stop-start)/step + 1e-9) + 1)]
    try:
        model = read_document(args.model)
        lines = read_lines(args.lines)
        with open(args.composition, encoding="utf-8") as stream:
            config = json.load(stream)
        # Calculate before opening the output so a failed model cannot leave a partial CSV.
        rows = list(opacities(model, wavelengths, lines, config, args.vmic))
        with args.output.open("x", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("layer", "wavelength_vacuum_A", "line_opacity_cm2_g",
                             "continuum_absorption_cm2_g", "continuum_scattering_cm2_g",
                             "electron_density_cm3", "density_g_cm3"))
            writer.writerows(rows)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
