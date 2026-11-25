"""Compare GiantLayers band-integrated flux with published MARCS samples."""

import argparse
from bisect import bisect_right
import csv
import gzip
import math
from pathlib import Path
import sys


DEFAULT_BANDS_NM = ((400, 700), (700, 1000), (1000, 2500), (2500, 5000), (5000, 10000))


def read_samples(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as stream:
        return [float(value) for value in stream.read().split()]


def compare(bands, wavelength_angstrom, flux_per_angstrom):
    """Yield comparisons for bands fully covered by the sampled reference."""
    if (len(wavelength_angstrom) != len(flux_per_angstrom) or len(wavelength_angstrom) < 2
            or any(a >= b for a, b in zip(wavelength_angstrom, wavelength_angstrom[1:]))
            or not all(math.isfinite(v) for v in wavelength_angstrom + flux_per_angstrom)):
        raise ValueError("MARCS flux and wavelength arrays must be finite and share a sorted grid")

    integral = [0.0]
    for left, right, a, b in zip(
        wavelength_angstrom, wavelength_angstrom[1:], flux_per_angstrom, flux_per_angstrom[1:]
    ):
        integral.append(integral[-1] + (right - left) * (a + b) / 2)

    def cumulative(wavelength):
        i = min(bisect_right(wavelength_angstrom, wavelength) - 1, len(integral) - 2)
        delta = wavelength - wavelength_angstrom[i]
        slope = ((flux_per_angstrom[i + 1] - flux_per_angstrom[i])
                 / (wavelength_angstrom[i + 1] - wavelength_angstrom[i]))
        return integral[i] + delta * (flux_per_angstrom[i] + 0.5 * delta * slope)

    for short, long, flux in bands:
        if wavelength_angstrom[0] <= short * 10 < long * 10 <= wavelength_angstrom[-1]:
            reference = cumulative(long * 10) - cumulative(short * 10)
            yield short, long, flux, reference, flux / reference if reference > 0 else ""


def sampled_bands(path, wavelengths, bands):
    """Integrate a native spectrum on exactly the supplied MARCS grid."""
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as stream:
        samples = [(float(row["wavelength_vacuum_A"]), float(row["reference_flux_per_A"]))
                   for row in csv.DictReader(stream)]
    if [row[0] for row in samples] != wavelengths:
        raise ValueError("GiantLayers samples must use the complete, ordered MARCS wavelength grid")
    fluxes = [row[1] for row in samples]
    requested = [(short, long, 0.0) for short, long in bands]
    result = [(short, long, integral) for short, long, _, integral, _
              in compare(requested, wavelengths, fluxes)]
    if len(result) != len(requested):
        raise ValueError("All requested bands must lie inside the shared wavelength grid")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bands", help="Band CSV, or native wavelength CSV with --sampled")
    parser.add_argument("wavelengths", help="MARCS flx_wavelengths.vac, optionally gzipped")
    parser.add_argument("flux", help="MARCS model .flx, optionally gzipped")
    parser.add_argument("--sampled", action="store_true",
                        help="Integrate the native --output-spectrum CSV on the full MARCS grid")
    parser.add_argument("--band", nargs=2, type=float, action="append", metavar=("SHORT_NM", "LONG_NM"),
                        help="Band edges for --sampled, repeatable")
    args = parser.parse_args()
    try:
        wavelengths = read_samples(args.wavelengths)
        if args.sampled:
            bands = sampled_bands(args.bands, wavelengths, args.band or DEFAULT_BANDS_NM)
        else:
            if args.band:
                raise ValueError("--band requires --sampled")
            with open(args.bands, newline="", encoding="utf-8") as stream:
                bands = [
                    tuple(float(row[key]) for key in (
                        "wavelength_min_nm", "wavelength_max_nm", "surface_flux_erg_cm2_s"
                    ))
                    for row in csv.DictReader(stream)
                ]
        comparison = list(compare(bands, wavelengths, read_samples(args.flux)))
        if not comparison:
            raise ValueError("No bands are fully covered by the MARCS wavelength grid")
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("wavelength_min_nm", "wavelength_max_nm", "giantlayers_flux_erg_cm2_s",
                     "marcs_flux_erg_cm2_s", "giantlayers_to_marcs"))
    writer.writerows(comparison)


if __name__ == "__main__":
    main()
