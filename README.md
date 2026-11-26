# GiantLayers

GiantLayers is an experimental spherical LTE stellar atmosphere solver for cool giants. It starts from effective temperature, surface gravity, stellar mass, and an elemental abundance mixture. The solver couples hydrostatic structure, LTE chemistry, wavelength-dependent radiative transfer with continuum scattering, and mixing-length convection. It writes a CSV with the luminosity residual at each layer and can export MARCS, ATLAS9, and ATLAS12 atmosphere files.

## Install and run

GiantLayers requires Python 3.9 or later. Install the alpha release from PyPI with `python -m pip install --pre giantlayers`. To install a repository checkout, use `python -m pip install .`.

Binary wheels use the included C transfer kernel automatically, with no compiler required to install or run the wheel. Use `--python-kernel` to select Python transfer or `--kernel PATH` for a custom shared library. On platforms without a matching binary wheel, pip can install the universal wheel and GiantLayers uses Python transfer.

The package includes a small solar-composition example. Run it from any directory after installing.

```sh
EXAMPLE_DIR=$(python -c "from importlib.resources import files; print(files('giantlayers.examples'))")
giantlayers "$EXAMPLE_DIR/continuum-wavelengths.dat" "$EXAMPLE_DIR/solar-gs98.json" \
  --teff 3500 --logg 1 --mass 1 --layers 12 --upper-knots 4 \
  --steps 1 --output-model trial.mod > trial.csv
```

This continuum-only example exercises the solver and exporter. For a line-blanketed calculation, provide an atomic transition CSV with `--lines` and molecular cross sections with repeatable `--molecular-cross-section` options. The output identifies unconverged profiles. Check the maximum luminosity error and sensitivity to wavelength and depth resolution before using a model.

See the [Quickstart](docs/quickstart.md) for installation and export, and [Advanced usage](docs/advanced.md) for opacity preparation and solver diagnostics. The gray Eddington benchmark is available through `python -m giantlayers.gray_atmosphere`. Run the regression tests with `python -m unittest discover`.

## Documentation

```sh
python3 -m pip install -r docs/requirements.txt
python3 -m mkdocs serve
```
