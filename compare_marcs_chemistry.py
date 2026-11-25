"""Compare the GiantLayers equation of state at fixed MARCS temperatures and pressures."""

import argparse
import csv
import math
from pathlib import Path
from statistics import median
import sys

from hminus import K_B
from lte_chemistry import Thermochemistry
from marcs_writer import PARTIAL_PRESSURES, read_document


SPECIES = ("electron", "CO", "CN", "H2O", "TiO")
PRESSURE_COLUMNS = {
    species: (group, column)
    for group, header in enumerate(PARTIAL_PRESSURES)
    for column, species in enumerate(header.replace("H I", "H").split()[1:])
}


def compare_chemistry(model, chemistry, minimum_depth=-4.0, maximum_depth=1.0):
    """Summarize log partial-pressure offsets in the selected Rosseland interval."""
    abundances = {element: model.abundances[record["atomic_number"]-1]
                  for element, record in chemistry.atoms.items()}
    offsets = {species: [] for species in SPECIES}
    for index, row in enumerate(model.structure):
        depth, _, _, temperature, electron_pressure, gas_pressure = row[:6]
        if not minimum_depth <= depth <= maximum_depth:
            continue
        gas = chemistry.solve(temperature, gas_pressure, abundances)
        thermal = K_B * temperature
        if electron_pressure > 0 and gas.electron_density > 0:
            offsets["electron"].append(math.log10(gas.electron_density*thermal/electron_pressure))
        for species in SPECIES[1:]:
            group, column = PRESSURE_COLUMNS[species]
            reference = model.partial_pressures[group][index][column]
            density = gas.number_densities[species]
            if reference > -90 and density > 0:
                offsets[species].append(math.log10(density*thermal)-reference)
    if not any(offsets.values()):
        raise ValueError("No comparable partial pressures in the selected Rosseland interval")
    for species, values in offsets.items():
        if values:
            rms = math.sqrt(math.fsum(x*x for x in values)/len(values))
            yield species, len(values), median(values), rms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("models", nargs="+", type=Path, help="Complete MARCS .mod or .mod.gz records")
    parser.add_argument("--min-log-tau", type=float, default=-4.0)
    parser.add_argument("--max-log-tau", type=float, default=1.0)
    args = parser.parse_args()
    if (not all(math.isfinite(depth) for depth in (args.min_log_tau, args.max_log_tau))
            or args.min_log_tau > args.max_log_tau):
        parser.error("Supply a finite, ordered optical depth interval")
    chemistry = Thermochemistry()
    try:
        rows = [(path.name, *row) for path in args.models for row in compare_chemistry(
            read_document(path), chemistry, args.min_log_tau, args.max_log_tau)]
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("model", "species", "layers", "median_log_pressure_difference_dex",
                     "rms_log_pressure_difference_dex"))
    writer.writerows(rows)


if __name__ == "__main__":
    main()
