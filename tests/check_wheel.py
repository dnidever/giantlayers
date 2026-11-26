"""Check a binary wheel from outside its source checkout."""

from importlib.resources import files
import math
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from giantlayers.spherical_kernel_loader import CompiledKernel, bundled_kernel_path
from giantlayers.spherical_heating import angular_geometry, solve_wavelength


def main():
    library = bundled_kernel_path()
    assert library is not None, "The wheel is missing the compiled kernel"
    rays = angular_geometry((10.0, 9.95, 9.9, 9.85, 9.8), 12)
    inputs = (rays, (0.05,)*4, (1,)*5, (1, 1.1, 1.3, 1.6, 1.7), (1e-5,)*5, (500,)*5)
    expected = solve_wavelength(*inputs, max_iterations=0)
    native = CompiledKernel(library, rays)(*inputs, max_iterations=0)
    for expected_values, actual_values in zip(expected[:3], native[:3]):
        for a, b in zip(expected_values, actual_values):
            assert math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12), (a, b)

    example = files("giantlayers.examples")
    with TemporaryDirectory() as directory:
        model = Path(directory) / "trial.mod"
        command = [
            sys.executable, "-m", "giantlayers",
            str(example.joinpath("continuum-wavelengths.dat")),
            str(example.joinpath("solar-gs98.json")),
            "--teff", "3500", "--logg", "1", "--mass", "1",
            "--layers", "12", "--upper-knots", "4", "--steps", "1",
            "--output-model", str(model),
        ]
        for options, backend in [([], "compiled"), (["--python-kernel"], "Python")]:
            result = subprocess.run(command + options, capture_output=True, text=True)
            assert result.returncode == 0, result.stderr
            assert f"Transfer kernel: {backend}" in result.stderr, result.stderr
            assert model.stat().st_size > 0 and result.stdout
            model.unlink()


if __name__ == "__main__":
    main()
