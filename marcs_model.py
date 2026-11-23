"""Read the layers needed for fixed-profile transfer from a MARCS .mod file."""

from dataclasses import dataclass
import gzip
import math
from pathlib import Path

from hminus import K_B


@dataclass(frozen=True)
class Layer:
    column_mass: float
    temperature: float
    gas_pressure: float
    electron_density: float
    neutral_h_per_gram: float


@dataclass(frozen=True)
class MarcsModel:
    identifier: str
    teff: float
    layers: tuple


def _table(lines, header, count):
    start = next((i + 1 for i, line in enumerate(lines) if header(line.split())), None)
    if start is None:
        raise ValueError("Missing required MARCS depth table")
    rows = []
    for index, line in enumerate(lines[start:start + count], start=1):
        values = line.split()
        if not values or values[0] != str(index):
            raise ValueError(f"Invalid MARCS depth row {index}")
        rows.append([float(value) for value in values[1:]])
    if len(rows) != count:
        raise ValueError("Truncated MARCS depth table")
    return rows


def read_model(path):
    """Read .mod or .mod.gz, including tabulated Pe and H I partial pressure."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as stream:
        lines = stream.readlines()
    try:
        identifier = lines[0].strip()
        teff = float(next(line.split()[0] for line in lines if "Teff [K]" in line))
        count = int(next(line.split()[0] for line in lines if "Number of depth points" in line))
        structure = _table(lines, lambda h: h[:7] == [
            "k", "lgTauR", "lgTau5", "Depth", "T", "Pe", "Pg"
        ], count)
        density = _table(lines, lambda h: h[:4] == ["k", "lgTauR", "KappaRoss", "Density"], count)
        partial = _table(lines, lambda h: h[:5] == ["k", "lgPgas", "H", "I", "H-"], count)
        layers = []
        for state, mass, species in zip(structure, density, partial):
            t, pe, pg = state[3:6]
            rho, column_mass = mass[2], mass[6]
            if not all(math.isfinite(v) and v > 0 for v in (t, pe, pg, rho, column_mass)):
                raise ValueError("Invalid thermodynamic quantity in MARCS model")
            kt = K_B * t
            layers.append(Layer(column_mass, t, pg, pe / kt, 10**species[1] / (kt * rho)))
    except (StopIteration, IndexError, OverflowError) as exc:
        raise ValueError("Incomplete or malformed MARCS model") from exc
    if (not math.isfinite(teff) or teff <= 0 or len(layers) < 2
            or any(left.column_mass >= right.column_mass for left, right in zip(layers, layers[1:]))):
        raise ValueError("MARCS model requires positive Teff and increasing column mass")
    return MarcsModel(identifier, teff, tuple(layers))
