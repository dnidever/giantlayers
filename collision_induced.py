"""HITRAN collision-induced absorption on a shared wavenumber grid."""

from array import array
from bisect import bisect_right
import gzip
import math
from pathlib import Path


class HitranCIA:
    """Collision-pair coefficients in cm⁵ per particle pair."""

    is_continuum = True

    def __init__(self, path):
        path = Path(path)
        opener = gzip.open if path.suffix == ".gz" else open
        temperatures = array("d")
        coefficients = []
        wavenumbers = array("d")
        pair = None
        with opener(path, "rt", encoding="ascii") as stream:
            for header in stream:
                if not header.strip():
                    continue
                species = tuple(header[:20].strip().split("-"))
                count = int(header[40:47])
                temperature = float(header[47:54])
                if (len(species) != 2 or not all(species) or count < 2
                        or not math.isfinite(temperature) or temperature <= 0
                        or (temperatures and temperature <= temperatures[-1])):
                    raise ValueError(f"Invalid or unordered CIA header in {path}")
                if pair is None:
                    pair = species
                if species != pair:
                    raise ValueError("Use one collision pair per CIA file")
                axis, values = array("d"), array("d")
                for _ in range(count):
                    fields = stream.readline().split()
                    if len(fields) != 2:
                        raise ValueError(f"Truncated CIA block in {path}")
                    wave, coefficient = map(float, fields)
                    if not all(math.isfinite(x) and x >= 0 for x in (wave, coefficient)):
                        raise ValueError(f"Invalid CIA sample in {path}")
                    axis.append(wave)
                    values.append(coefficient)
                if not wavenumbers:
                    if any(a >= b for a, b in zip(axis, axis[1:])):
                        raise ValueError("CIA wavenumbers must increase")
                    wavenumbers = axis
                if axis != wavenumbers:
                    raise ValueError("Use a CIA file with matching wavenumber grids at every temperature")
                temperatures.append(temperature)
                coefficients.append(values)
        if not temperatures:
            raise ValueError(f"Empty CIA file {path}")
        self.species = pair
        self.temperatures = temperatures
        self.wavenumbers = wavenumbers
        self.coefficients = coefficients

    def on_grid(self, wavelengths):
        grid = []
        for wavelength in wavelengths:
            wave = 1e8/wavelength
            if not self.wavenumbers[0] <= wave <= self.wavenumbers[-1]:
                grid.append(None)
                continue
            i = min(bisect_right(self.wavenumbers, wave), len(self.wavenumbers)-1) - 1
            lo, hi = self.wavenumbers[i:i+2]
            grid.append((i, (wave-lo)/(hi-lo)))
        return tuple(grid)

    def mass_opacity(self, gas, pressure, grid):
        """Return absorption in cm²/g from the two LTE number densities."""
        density_a, density_b = (gas.number_densities.get(atom, 0.0) for atom in self.species)
        factor = density_a*density_b/gas.mass_density
        if not factor:
            return [0.0] * len(grid)
        t = max(self.temperatures[0], min(self.temperatures[-1], gas.temperature))
        hi = min(bisect_right(self.temperatures, t), len(self.temperatures)-1)
        lo = max(0, hi-1)
        weight = (t-self.temperatures[lo])/(self.temperatures[hi]-self.temperatures[lo]) if hi != lo else 0.0
        a, b = self.coefficients[lo], self.coefficients[hi]
        output = []
        for point in grid:
            if point is None:
                output.append(0.0)
                continue
            i, frequency_weight = point
            left = (1-frequency_weight)*a[i] + frequency_weight*a[i+1]
            right = (1-frequency_weight)*b[i] + frequency_weight*b[i+1]
            output.append(factor*((1-weight)*left + weight*right))
        return output
