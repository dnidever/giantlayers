"""Read and interpolate a rectangular two-group mass-opacity table."""

from bisect import bisect_right
import csv
import math


class OpacityTable:
    """Bilinear interpolation in log temperature, pressure, and opacity."""

    def __init__(self, path):
        columns = (
            "temperature_K", "pressure_dyn_cm2",
            "group1_opacity_cm2_g", "group2_opacity_cm2_g",
        )
        values = {}
        temperatures = set()
        pressures = set()
        with open(path, newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or not set(columns) <= set(reader.fieldnames):
                raise ValueError(f"Opacity table requires columns: {', '.join(columns)}")
            for line, row in enumerate(reader, start=2):
                try:
                    t, p, k1, k2 = (float(row[column]) for column in columns)
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Invalid opacity table value on line {line}") from exc
                if not all(math.isfinite(x) and x > 0 for x in (t, p, k1, k2)):
                    raise ValueError(f"Opacity table line {line} requires positive finite values")
                if (t, p) in values:
                    raise ValueError(f"Duplicate temperature and pressure on line {line}")
                temperatures.add(t)
                pressures.add(p)
                values[t, p] = (math.log(k1), math.log(k2))

        self.temperatures = sorted(temperatures)
        self.pressures = sorted(pressures)
        if (min(len(temperatures), len(pressures)) < 2
                or len(values) != len(temperatures) * len(pressures)):
            raise ValueError("Opacity table must be a complete rectangular grid of at least 2 × 2")
        self.values = values

    @staticmethod
    def _bracket(axis, value, name):
        if not math.isfinite(value) or not axis[0] <= value <= axis[-1]:
            raise ValueError(f"{name} {value:g} outside opacity table range [{axis[0]:g}, {axis[-1]:g}]")
        i = min(bisect_right(axis, value) - 1, len(axis) - 2)
        low, high = axis[i:i + 2]
        fraction = math.log(value / low) / math.log(high / low)
        return low, high, fraction

    def __call__(self, temperature, pressure):
        t0, t1, ft = self._bracket(self.temperatures, temperature, "Temperature (K)")
        if not math.isfinite(pressure) or pressure < 0:
            raise ValueError("Pressure must be finite and nonnegative")
        # At the vacuum boundary and below the lowest tabulated pressure,
        # hold the mass opacity at the lowest-pressure row.
        p0, p1, fp = self._bracket(
            self.pressures, max(pressure, self.pressures[0]), "Pressure (dyn/cm²)"
        )
        result = []
        for group in range(2):
            low = (1 - ft) * self.values[t0, p0][group] + ft * self.values[t1, p0][group]
            high = (1 - ft) * self.values[t0, p1][group] + ft * self.values[t1, p1][group]
            result.append(math.exp((1 - fp) * low + fp * high))
        return tuple(result)
