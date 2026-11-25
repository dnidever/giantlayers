"""Convert a local ExoMolOP TauREx opacity table to GiantLayers JSON.gz."""

import argparse
import gzip
import json
import math
from pathlib import Path
import os
from tempfile import NamedTemporaryFile


def select(axis, limits, scale, name, points=None):
    values = [float(x)*scale for x in axis]
    if (not all(math.isfinite(x) and x > 0 for x in values)
            or any(a >= b for a, b in zip(values, values[1:]))):
        raise ValueError(f"{name} values must be positive and increasing")
    if points is None:
        indices = [i for i, x in enumerate(values) if limits is None or limits[0] <= x <= limits[1]]
    else:
        if limits is not None:
            raise ValueError(f"Choose either a {name.lower()} range or tabulated points")
        selected = set()
        for point in points:
            matches = [i for i, x in enumerate(values) if math.isclose(x, point, rel_tol=1e-8)]
            if len(matches) != 1:
                raise ValueError(f"{name} point {point:g} does not match one tabulated sample")
            selected.add(matches[0])
        indices = sorted(selected)
    if len(indices) < 2:
        raise ValueError(f"Select at least two {name.lower()} samples")
    return indices, [values[i] for i in indices]


def convert(source, destination, species, isotope_fraction=1.0,
            temperature_range=None, pressure_range=None, source_url=None, *,
            temperature_points=None, pressure_points=None):
    """Read T in K, pressure in bar, wavenumber in cm⁻¹ and cross section in cm²."""
    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError("Install h5py to convert TauREx input; the solver reads the resulting JSON without h5py") from exc
    if not 0 < isotope_fraction <= 1:
        raise ValueError("Use an isotopic fraction in (0, 1]")
    if not species:
        raise ValueError("Provide a chemical species name")
    destination = Path(destination)
    output = None
    try:
        with h5py.File(source, "r") as table:
            box = table["xsecarr"]
            ti, temperatures = select(table["t"][:], temperature_range, 1.0, "Temperature", temperature_points)
            pi, pressures = select(table["p"][:], pressure_range, 1e6, "Pressure", pressure_points)
            _, wavenumbers = select(table["bin_edges"][:], None, 1.0, "Wavenumber")
            units = box.attrs.get("units")
            if isinstance(units, bytes):
                units = units.decode("ascii")
            if units != "cm^2/molecule" or box.shape != (len(table["p"]), len(table["t"]), len(wavenumbers)):
                raise ValueError("Expected TauREx opacity in cm^2/molecule ordered pressure × temperature × wavenumber")
            isotope = table["key_iso_ll"][0]
            if isinstance(isotope, bytes):
                isotope = isotope.decode("ascii")
            meta = dict(
                species=species, isotopologue=str(isotope), isotope_fraction=isotope_fraction,
                source=source_url or Path(source).name, temperature_K=temperatures,
                pressure_dyn_cm2=pressures,
                wavelength_vacuum_A=[1e8/wn for wn in reversed(wavenumbers)],
            )
            with NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
                output = Path(handle.name)
                with gzip.open(handle, "wt", encoding="utf-8") as stream:
                    stream.write(json.dumps(meta, separators=(",", ":"))[:-1] + ',"cross_section_cm2":[')
                    for i, t in enumerate(ti):
                        stream.write(",[" if i else "[")
                        for j, p in enumerate(pi):
                            if j:
                                stream.write(",")
                            row = box[p, t, :].tolist()[::-1]
                            if not all(math.isfinite(v) and v >= 0 for v in row):
                                raise ValueError(f"Invalid cross section at T={temperatures[i]:g} K and P={pressures[j]:g} dyn/cm²")
                            json.dump(row, stream, separators=(",", ":"), allow_nan=False)
                        stream.write("]")
                    stream.write("]}\n")
        os.replace(output, destination)
        output = None
    finally:
        if output is not None:
            output.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path, help="Output path ending in .json.gz")
    parser.add_argument("--species", required=True)
    parser.add_argument("--isotope-fraction", type=float, default=1.0)
    temperature = parser.add_mutually_exclusive_group()
    temperature.add_argument("--temperature-range", type=float, nargs=2, metavar=("MIN_K", "MAX_K"))
    temperature.add_argument("--temperature-points", type=float, nargs="+", metavar="K")
    pressure = parser.add_mutually_exclusive_group()
    pressure.add_argument("--pressure-range", type=float, nargs=2, metavar=("MIN_CGS", "MAX_CGS"))
    pressure.add_argument("--pressure-points", type=float, nargs="+", metavar="CGS")
    parser.add_argument("--source-url")
    args = parser.parse_args()
    if args.destination.suffix != ".gz":
        parser.error("Output path must end in .gz")
    try:
        convert(args.source, args.destination, args.species, args.isotope_fraction,
                args.temperature_range, args.pressure_range, args.source_url,
                temperature_points=args.temperature_points, pressure_points=args.pressure_points)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
