"""Interpolate portable molecular cross-section tables without external packages."""

from array import array
from bisect import bisect_right
import gzip
import json
import math
from pathlib import Path


def _axis(values, name):
    if (len(values) < 2 or not all(math.isfinite(x) and x > 0 for x in values)
            or any(a >= b for a, b in zip(values, values[1:]))):
        raise ValueError(f"{name} must contain increasing positive values")
    return tuple(values)


def _bracket(axis, value):
    value = max(axis[0], min(axis[-1], value))
    i = min(bisect_right(axis, value) - 1, len(axis) - 2)
    return i, (value - axis[i]) / (axis[i + 1] - axis[i])


class MolecularCrossSection:
    """Cross sections in cm² per isotopologue on a T, pressure, wavelength grid."""

    def __init__(self, path, *, omit_hot=False):
        self.omit_hot = omit_hot
        path = Path(path)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
        self.species = data["species"]
        self.isotopologue = data.get("isotopologue", self.species)
        self.isotope_fraction = data.get("isotope_fraction", 1.0)
        if not 0 < self.isotope_fraction <= 1:
            raise ValueError("Isotopologue fraction must lie in (0, 1]")
        self.temperatures = _axis(data["temperature_K"], "Temperature axis")
        self.pressures = _axis(data["pressure_dyn_cm2"], "Pressure axis")
        self.log_pressures = tuple(map(math.log, self.pressures))
        self.wavelengths = _axis(data["wavelength_vacuum_A"], "Wavelength axis")
        sections = data["cross_section_cm2"]
        if (len(sections) != len(self.temperatures)
                or any(len(layer) != len(self.pressures) for layer in sections)
                or any(len(row) != len(self.wavelengths)
                       or not all(math.isfinite(value) and value >= 0 for value in row)
                       for layer in sections for row in layer)):
            raise ValueError("Cross sections must be a complete nonnegative T × P × wavelength grid")
        self.cross_sections = tuple(tuple(array("d", row) for row in layer) for layer in sections)

    def on_grid(self, wavelengths):
        """Precompute wavelength interpolation for a transfer grid."""
        return tuple(
            _bracket(self.wavelengths, wavelength)
            if self.wavelengths[0] <= wavelength <= self.wavelengths[-1] else None
            for wavelength in wavelengths
        )

    def clamped_layers(self, structure):
        """Count layers where cross-section interpolation clamps to an axis edge."""
        cold = hot = low = high = 0
        for row in structure:
            cold += row[3] < self.temperatures[0]
            hot += row[3] > self.temperatures[-1]
            low += row[5] < self.pressures[0]
            high += row[5] > self.pressures[-1]
        return cold, hot, low, high

    def mass_opacity(self, gas, pressure, grid):
        """Return molecular absorption on the requested grid in cm²/g."""
        density = gas.number_densities.get(self.species)
        if density is None:
            raise ValueError(f"Missing chemistry data for {self.species}")
        if density == 0 or (self.omit_hot and gas.temperature > self.temperatures[-1]):
            return [0.0] * len(grid)
        it, ft = _bracket(self.temperatures, gas.temperature)
        ip, fp = _bracket(self.log_pressures, math.log(pressure))
        table = self.cross_sections
        a, b, c, d = table[it][ip], table[it][ip + 1], table[it + 1][ip], table[it + 1][ip + 1]
        wa, wb, wc, wd = (1-ft)*(1-fp), (1-ft)*fp, ft*(1-fp), ft*fp
        abundance = density * self.isotope_fraction / gas.mass_density
        result = []
        for point in grid:
            if point is None:
                result.append(0.0)
                continue
            i, fw = point
            left = wa*a[i] + wb*b[i] + wc*c[i] + wd*d[i]
            j = i + 1
            right = wa*a[j] + wb*b[j] + wc*c[j] + wd*d[j]
            result.append(abundance * ((1 - fw) * left + fw * right))
        return result
