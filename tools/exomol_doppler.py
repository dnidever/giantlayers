"""Generate experimental pressure-independent LTE tables from ExoMol lines."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import numpy as np


def wavenumber_edges(wavelengths):
    """Return ascending geometric bin boundaries in inverse centimetres."""
    centers = 1e8 / np.asarray(wavelengths, dtype=float)[::-1]
    if (len(centers) < 2 or not np.all(np.isfinite(centers))
            or centers[0] <= 0 or not np.all(np.diff(centers) > 0)):
        raise ValueError("Template wavelengths must be positive and strictly increasing")
    inner = np.sqrt(centers[1:] * centers[:-1])
    return np.r_[centers[0]**2 / inner[0], inner, centers[-1]**2 / inner[-1]]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template", type=Path, help="Portable GiantLayers table (.json or .json.gz)")
    parser.add_argument("states", type=Path, help="ExoMol .states.bz2 file")
    parser.add_argument("transitions", type=Path, help="ExoMol .trans.bz2 file")
    parser.add_argument("partition", type=Path, help="ExoMol .pf file")
    parser.add_argument("output", type=Path, help="Output compressed JSON file")
    parser.add_argument("--kernel", type=Path, required=True, help="Compiled exomol_doppler.c")
    parser.add_argument("--mass-amu", type=float, required=True)
    parser.add_argument("--vmic", type=float, default=2.0, help="Microturbulent speed in km/s")
    parser.add_argument("--extend-temperatures", nargs="+", type=float, default=[])
    parser.add_argument("--capped-control", type=Path, help="Write the same calculation on the template temperatures")
    args = parser.parse_args()

    if (not np.isfinite(args.mass_amu) or args.mass_amu <= 0
            or not np.isfinite(args.vmic) or args.vmic < 0):
        parser.error("Use positive isotopologue mass and nonnegative microturbulence")
    opener = gzip.open if args.template.suffix == ".gz" else open
    with opener(args.template, "rt", encoding="utf-8") as stream:
        template = json.load(stream)
    edges = wavenumber_edges(template["wavelength_vacuum_A"])
    original = template["temperature_K"]
    temperatures = [*original, *args.extend_temperatures]
    if not np.all(np.isfinite(temperatures)) or not np.all(np.diff(temperatures) > 0):
        parser.error("Temperatures must be finite and strictly increasing")
    states = np.loadtxt(args.states, usecols=(0, 1, 2), ndmin=2)
    partition = np.loadtxt(args.partition, usecols=(0, 1), ndmin=2)
    if (temperatures[0] < partition[0, 0] or temperatures[-1] > partition[-1, 0]
            or not np.all(np.diff(partition[:, 0]) > 0)):
        parser.error("Partition-function temperatures must be sorted and cover the requested range")
    partitions = np.interp(temperatures, partition[:, 0], partition[:, 1])
    bins = len(edges) - 1
    with TemporaryDirectory() as directory:
        config = Path(directory) / "input"
        with config.open("w", encoding="ascii") as stream:
            print(bins, len(temperatures), int(states[:, 0].max()), args.mass_amu, args.vmic*1e5,
                  file=stream)
            np.savetxt(stream, edges, fmt="%.16g")
            np.savetxt(stream, np.column_stack((temperatures, partitions)), fmt="%.16g")
            np.savetxt(stream, states, fmt=["%d", "%.12g", "%d"])
        with (Path(directory) / "cross-sections").open("w+b") as stream:
            with subprocess.Popen(["bzcat", args.transitions], stdout=subprocess.PIPE) as source:
                result = subprocess.run([args.kernel.resolve(), config], stdin=source.stdout, stdout=stream)
                source.stdout.close()
                if source.wait() or result.returncode:
                    parser.error("ExoMol line integration failed")
            stream.seek(0)
            sigma = np.fromfile(stream).reshape(len(temperatures), bins)
    pressures = template["pressure_dyn_cm2"]
    layers = [[row[::-1].tolist()] * len(pressures) for row in sigma]
    header = {key: template[key] for key in ("species", "pressure_dyn_cm2", "wavelength_vacuum_A")}
    header["isotopologue"] = template.get("isotopologue", template["species"])
    header["isotope_fraction"] = template.get("isotope_fraction", 1.0)
    header["provenance"] = {
        "profile": "Thermal and microturbulent Doppler, geometric wavenumber-bin averages",
        "mass_amu": args.mass_amu,
        "vmic_km_s": args.vmic,
        "input_sha256": {path.name: sha256(path) for path in (args.states, args.transitions, args.partition)},
    }
    outputs = [(args.output, len(temperatures))]
    if args.capped_control:
        outputs.append((args.capped_control, len(original)))
    for path, count in outputs:
        with gzip.open(path, "wt", encoding="utf-8", compresslevel=3) as stream:
            json.dump(dict(header, temperature_K=temperatures[:count], cross_section_cm2=layers[:count]),
                      stream, separators=(",", ":"), allow_nan=False)


if __name__ == "__main__":
    main()
