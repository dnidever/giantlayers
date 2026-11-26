"""Convert Kurucz GFALL fixed-width atomic records to GiantLayers transitions."""

import argparse
import csv
import math
from pathlib import Path
import sys

from lte_chemistry import Thermochemistry
from marcs_atomic_lines import air_to_vacuum


COLUMNS = (
    "wavelength_vacuum_A", "species", "lower_energy_eV", "log_gf", "damping_s",
    "electron_damping_cm3_s", "neutral_h_damping_cm3_s",
    "electron_damping_exponent", "neutral_h_damping_exponent",
)


def read_kurucz(paths, *, minimum_nm=0, maximum_nm=math.inf, vacuum=False, doppler_only=False):
    """Read neutral and singly ionized lines from standard wavelength-format records.

    Bounds are in the wavelength medium of the input. Standard Kurucz lists
    tabulate vacuum wavelengths below 200 nm and air wavelengths above it.
    Supply vacuum=True for lists already tabulated in vacuum throughout.
    """
    elements = {atom["atomic_number"]: symbol for symbol, atom in Thermochemistry().atoms.items()}
    for path in paths:
        with Path(path).open(encoding="ascii") as stream:
            for row, record in enumerate(stream, 1):
                if not record.strip():
                    continue
                try:
                    wavelength = float(record[:11])
                    if not math.isfinite(wavelength) or wavelength <= 0:
                        raise ValueError("invalid wavelength")
                    if not minimum_nm <= wavelength <= maximum_nm:
                        continue
                    species = float(record[18:24])
                    atomic_number = int(species)
                    ion = round(100*(species-atomic_number))
                    if ion not in (0, 1) or atomic_number not in elements:
                        continue
                    log_gf = float(record[11:18])
                    low_energy = min(abs(float(record[24:36])), abs(float(record[52:64])))
                    rates = tuple(float(record[i:i+6].strip() or 0) for i in (80, 86, 92))
                    if not all(map(math.isfinite, (log_gf, low_energy, *rates))):
                        raise ValueError("invalid strength, energy, or damping coefficient")
                    damping = tuple(10.0**value if value and not doppler_only else 0.0 for value in rates)
                except (ValueError, OverflowError) as exc:
                    raise ValueError(f"Invalid Kurucz record {row} in {path}") from exc
                wave_a = 10*wavelength
                if not vacuum and wavelength >= 200:
                    wave_a = air_to_vacuum(wave_a)
                yield (wave_a, elements[atomic_number]+("+" if ion else ""),
                       low_energy/8065.544005, log_gf, *damping, 0.0, 0.3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", nargs="+", type=Path)
    parser.add_argument("--wavelength-nm", type=float, nargs=2, metavar=("MIN", "MAX"))
    parser.add_argument("--vacuum", action="store_true", help="All input wavelengths are in vacuum")
    parser.add_argument("--doppler-only", action="store_true", help="Omit radiative and collisional damping")
    args = parser.parse_args()
    lo, hi = args.wavelength_nm or (0, math.inf)
    if not 0 <= lo < hi:
        parser.error("Supply increasing nonnegative wavelength bounds")
    writer = csv.writer(sys.stdout)
    writer.writerow(COLUMNS)
    try:
        writer.writerows(read_kurucz(args.tables, minimum_nm=lo, maximum_nm=hi,
                                     vacuum=args.vacuum, doppler_only=args.doppler_only))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
