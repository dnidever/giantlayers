# GiantLayers

GiantLayers is an experimental spherical LTE stellar atmosphere solver for cool giants. It starts from effective temperature, surface gravity, stellar mass, and an elemental abundance mixture. The solver couples hydrostatic structure, LTE chemistry, wavelength-dependent radiative transfer with continuum scattering, and mixing-length convection. It writes a CSV with the luminosity residual at each layer and can export a MARCS-format model.

## Install and run

From a repository checkout, install with Python 3.9 or newer.

```sh
python3 -m pip install .
giantlayers examples/continuum-wavelengths.dat examples/solar-gs98.json \
  --teff 3500 --logg 1 --mass 1 --layers 12 --upper-knots 4 \
  --steps 1 --output-model trial.mod > trial.csv
```

This continuum-only example exercises the solver and exporter. For a line-blanketed calculation, provide an atomic transition CSV with `--lines` and molecular cross sections with repeatable `--molecular-cross-section` options. The output identifies unconverged profiles. Check the maximum luminosity error and sensitivity to wavelength and depth resolution before using a model.

See the [quickstart](docs/quickstart.md) for input formats, opacity sources, and the optional compiled transfer kernel. The gray Eddington benchmark is available through `python3 gray_atmosphere.py`. Run the regression tests with `python3 -m unittest discover`.

## Documentation

```sh
python3 -m pip install -r docs/requirements.txt
python3 -m mkdocs serve
```
