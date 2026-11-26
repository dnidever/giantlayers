"""Convert the revised MARCS VALD oscillator strengths to a transition CSV."""

import argparse
import csv
from pathlib import Path
import sys

from .lte_chemistry import Thermochemistry


def air_to_vacuum(air_angstrom):
    """Convert optical standard-air wavelengths using the Edlén refractivity."""
    wavenumber_squared = (10000 / air_angstrom)**2
    refractivity = (8342.13 + 2406030 / (130 - wavenumber_squared)
                    + 15997 / (38.9 - wavenumber_squared)) * 1e-8
    return air_angstrom * (1 + refractivity)


def read_tables(paths):
    elements = {data["atomic_number"]: element
                for element, data in Thermochemistry().atoms.items()}
    transitions = {}
    for path in paths:
        with Path(path).open(encoding="ascii") as stream:
            for line_number, line in enumerate(stream, 1):
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                try:
                    z, ion, air, energy, _, log_gf = line.split()[:6]
                    species = elements[int(z)]
                    stage = int(ion)
                    air, energy, log_gf = map(float, (air, energy, log_gf))
                    if stage not in (1, 2) or air < 2000 or energy < 0:
                        raise ValueError("Unsupported ion stage, wavelength, or excitation")
                except (KeyError, ValueError) as exc:
                    raise ValueError(f"Invalid VALD row {line_number} in {path}") from exc
                if log_gf == -9.99:
                    continue
                species += "+" if stage == 2 else ""
                transitions[(air, species, energy)] = log_gf

    for (air, species, energy), log_gf in sorted(transitions.items()):
        yield air_to_vacuum(air), species, energy, log_gf


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", nargs="+", help="Published astroph_gfVALD*.dat files, oldest first")
    args = parser.parse_args()
    try:
        rows = list(read_tables(args.tables))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("wavelength_vacuum_A", "species", "lower_energy_eV", "log_gf"))
    writer.writerows(rows)


if __name__ == "__main__":
    main()
