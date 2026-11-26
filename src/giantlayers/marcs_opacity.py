"""Read depth-dependent continuum absorption and scattering from MARCS."""

from bisect import bisect_right
import gzip
import math
from pathlib import Path


def take(tokens, count):
    return [next(tokens) for _ in range(count)]


class MarcsOpacity:
    """Wavelength-dependent continuum opacity for a matching MARCS model."""

    def __init__(self, path):
        path = Path(path)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="ascii") as stream:
            tokens = iter(stream.read().split())
        try:
            label = next(tokens)
            count = int(next(tokens))
            reference_wavelength = float(next(tokens))
            wavelengths = int(next(tokens))
            if label != "MRXF" or count < 2 or wavelengths < 2:
                raise ValueError("Invalid MARCS continuum header")
            self.wavelengths = tuple(float(value) / 10 for value in take(tokens, wavelengths))
            if (not all(math.isfinite(w) and w > 0 for w in self.wavelengths)
                    or any(a >= b for a, b in zip(self.wavelengths, self.wavelengths[1:]))):
                raise ValueError("MARCS continuum wavelengths must be positive and increasing")
            self.layers = []
            for _ in range(count):
                geometry, tau, t, pe, pg, rho, velocity, reference = (
                    float(value) for value in take(tokens, 8)
                )
                ratios = tuple(map(float, take(tokens, 2 * wavelengths)))
                if (not math.isfinite(reference) or reference <= 0
                        or not all(math.isfinite(x) and x >= 0 for x in ratios)):
                    raise ValueError("Invalid MARCS continuum opacity in depth table")
                self.layers.append((t, pg, reference, ratios))
        except StopIteration as exc:
            raise ValueError("Truncated MARCS continuum opacity table") from exc
        self.reference_wavelength_nm = reference_wavelength / 10

    def check_profile(self, model):
        """Check that the continuum and structure files describe the same layers."""
        if len(model.layers) != len(self.layers) or any(
            not (math.isclose(layer.temperature, t, rel_tol=5e-5)
                 and math.isclose(layer.gas_pressure, pg, rel_tol=5e-4))
            for layer, (t, pg, _, _) in zip(model.layers, self.layers)
        ):
            raise ValueError("MARCS continuum temperature and pressure do not match the model")

    def at_wavelength(self, wavelength_nm):
        """Return absorption and scattering in cm²/g for each layer."""
        axis = self.wavelengths
        wavelength_nm = min(max(wavelength_nm, axis[0]), axis[-1])
        i = min(bisect_right(axis, wavelength_nm) - 1, len(axis) - 2)
        fraction = (wavelength_nm - axis[i]) / (axis[i + 1] - axis[i])
        for _, _, scale, ratios in self.layers:
            low = 2 * i
            yield tuple(
                scale * (ratios[j] + fraction * (ratios[j + 2] - ratios[j]))
                for j in (low, low + 1)
            )
