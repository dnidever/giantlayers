"""Compare two MARCS-format structures over their shared Rosseland depth range."""

import argparse
from bisect import bisect_right
import csv
import math
from pathlib import Path
import sys

from .marcs_writer import read_document


def compare(reference, model):
    """Interpolate temperature and log pressures onto the reference depth grid."""
    # The format rounds log optical depth. Keep the innermost value of any
    # runs of equal labels so that each interpolation cell has nonzero width.
    rows = list({row[0]: row for row in model.structure}.values())
    depths = [row[0] for row in rows]
    if len(depths) < 2 or any(a >= b for a, b in zip(depths, depths[1:])):
        raise ValueError("The model must have increasing Rosseland depth")
    if any(row[4] <= 0 or row[5] <= 0 for row in rows):
        raise ValueError("The model must have positive electron and gas pressures")
    for index, row in enumerate(reference.structure, 1):
        depth = row[0]
        if not depths[0] <= depth <= depths[-1]:
            continue
        high = min(bisect_right(depths, depth), len(depths) - 1)
        a, b = rows[high-1:high+1]
        weight = (depth-a[0])/(b[0]-a[0])
        temperature = (1-weight)*a[3] + weight*b[3]
        pe = (1-weight)*math.log10(a[4]) + weight*math.log10(b[4])
        pg = (1-weight)*math.log10(a[5]) + weight*math.log10(b[5])
        if min(row[3:6]) <= 0:
            raise ValueError("The reference requires positive temperature and pressures")
        yield index, depth, row[3], temperature, temperature-row[3], pe-math.log10(row[4]), pg-math.log10(row[5])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("model", type=Path)
    args = parser.parse_args()
    try:
        rows = list(compare(read_document(args.reference), read_document(args.model)))
        if not rows:
            raise ValueError("The models have no shared Rosseland depth points")
        writer = csv.writer(sys.stdout)
        writer.writerow(("reference_layer", "log_tau_ross", "reference_temperature_K", "model_temperature_K",
                         "temperature_difference_K", "electron_pressure_difference_dex", "gas_pressure_difference_dex"))
        writer.writerows(rows)
        rms = math.sqrt(math.fsum(row[4]**2 for row in rows)/len(rows))
        print(f"Compared {len(rows)} layers. Temperature RMS difference {rms:.1f} K.", file=sys.stderr)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
