"""Compare luminosity closure and structure across GiantLayers depth grids."""

import argparse
from bisect import bisect_right
import csv
import json
import math
from pathlib import Path


MASS = "column_mass_g_cm2"
TEMPERATURE = "temperature_K"
PRESSURE = "gas_pressure_dyn_cm2"
RESIDUAL = "relative_total_flux_error"


def read_profile(path):
    """Read the solver's CSV on its radial column-mass coordinate."""
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fields = (MASS, TEMPERATURE, PRESSURE, RESIDUAL)
        if not set(fields) <= set(reader.fieldnames or ()):
            raise ValueError(f"{path} lacks required solver CSV columns")
        try:
            rows = [tuple(float(row[field]) for field in fields) for row in reader]
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{path} contains a missing or invalid numeric value") from exc
    if (
        len(rows) < 2
        or any(not all(map(math.isfinite, row)) or row[0] < 0 or min(row[1:3]) <= 0 for row in rows)
        or any(a[0] >= b[0] for a, b in zip(rows, rows[1:]))
    ):
        raise ValueError(f"{path} requires increasing nonnegative mass and positive finite structure")
    return rows


def interpolate(axis, rows, mass):
    hi = min(bisect_right(axis, mass), len(axis)-1)
    a, b = rows[hi-1:hi+1]
    weight = (mass-a[0])/(b[0]-a[0])
    return tuple(x + weight*(y-x) for x, y in zip(a[1:], b[1:]))


def summarize(rows):
    index, worst = max(enumerate(rows, 1), key=lambda pair: abs(pair[1][3]))
    return {
        "layers": len(rows),
        "rms_luminosity_error": math.sqrt(math.fsum(r[3]**2 for r in rows)/len(rows)),
        "max_absolute_luminosity_error": abs(worst[3]),
        "worst_layer": index,
        "worst_column_mass_g_cm2": worst[0],
        "surface_luminosity_error": rows[0][3],
        "bottom_luminosity_error": rows[-1][3],
    }


def audit(coarse, fine, flux_tolerance=0.01, temperature_tolerance=None, pressure_tolerance=None):
    """Measure structure changes on the union of mass nodes in the overlap."""
    axes = [[row[0] for row in profile] for profile in (coarse, fine)]
    lo, hi = max(axis[0] for axis in axes), min(axis[-1] for axis in axes)
    if lo >= hi:
        raise ValueError("The two profiles must share a nonzero column-mass interval")
    nodes = sorted({m for axis in axes for m in axis if lo <= m <= hi})
    a, b = [[(m, t, math.log10(p)) for m, t, p, _ in profile] for profile in (coarse, fine)]
    differences = []
    for m in nodes:
        t1, p1 = interpolate(axes[0], a, m)
        t2, p2 = interpolate(axes[1], b, m)
        differences.append((m, abs(t2-t1), abs(p2-p1)))
    worst_t = max(differences, key=lambda row: row[1])
    worst_p = max(differences, key=lambda row: row[2])
    summaries = [summarize(rows) for rows in (coarse, fine)]
    checks = {
        "same_mass_extent": all(math.isclose(axes[0][i], axes[1][i], rel_tol=1e-9)
                                for i in (0, -1)),
        "both_in_luminosity_balance": all(
            s["max_absolute_luminosity_error"] <= flux_tolerance for s in summaries
        ),
    }
    if temperature_tolerance is not None:
        checks["temperature_change_within_tolerance"] = worst_t[1] <= temperature_tolerance
    if pressure_tolerance is not None:
        checks["gas_pressure_change_within_tolerance"] = worst_p[2] <= pressure_tolerance
    return {
        "coarse": summaries[0], "fine": summaries[1],
        "overlap_column_mass_g_cm2": [lo, hi], "comparison_nodes": len(nodes),
        "max_temperature_difference_K": worst_t[1], "temperature_difference_mass_g_cm2": worst_t[0],
        "max_gas_pressure_difference_dex": worst_p[2], "gas_pressure_difference_mass_g_cm2": worst_p[0],
        "tolerances": {"luminosity_fraction": flux_tolerance, "temperature_K": temperature_tolerance,
                       "gas_pressure_dex": pressure_tolerance},
        "checks": checks,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coarse", type=Path)
    parser.add_argument("fine", type=Path)
    parser.add_argument("--flux-tolerance", type=float, default=0.01)
    parser.add_argument("--temperature-tolerance", type=float, help="Maximum temperature change in K")
    parser.add_argument("--pressure-tolerance", type=float, help="Maximum log gas-pressure change in dex")
    parser.add_argument("--require-pass", action="store_true", help="Return 2 if any configured check fails")
    args = parser.parse_args()
    tolerances = (args.flux_tolerance, args.temperature_tolerance, args.pressure_tolerance)
    if any(x is not None and (not math.isfinite(x) or x < 0) for x in tolerances):
        parser.error("Tolerances must be finite and nonnegative")
    try:
        result = audit(read_profile(args.coarse), read_profile(args.fine), *tolerances)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, indent=2, allow_nan=False))
    return 2 if args.require_pass and not all(result["checks"].values()) else 0


if __name__ == "__main__":
    raise SystemExit(main())
