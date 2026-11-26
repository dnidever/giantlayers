# Quickstart

## Install from PyPI

GiantLayers requires Python 3.9 or later. The current release is an alpha.

```sh
python -m pip install --pre giantlayers
giantlayers --help
```

On supported platforms, pip installs a wheel with the C transfer kernel, and GiantLayers uses it automatically. The universal wheel uses Python transfer and also installs without a compiler.

## Run the bundled example

Copy the included continuum wavelength grid and GS98 solar abundance mixture into an empty working directory.

```sh
python -c "from importlib.resources import files; from pathlib import Path; examples = files('giantlayers.examples'); names = ('continuum-wavelengths.dat', 'solar-gs98.json'); [Path(name).write_bytes(examples.joinpath(name).read_bytes()) for name in names]"
giantlayers continuum-wavelengths.dat solar-gs98.json --teff 3500 --logg 1 --mass 1 --layers 12 --upper-knots 4 --steps 1 --output-model trial.mod > trial.csv
```

This starts a spherical atmosphere at 3,500 K, log g = 1 in cgs units, and one solar mass. It runs one temperature correction with continuum opacity and writes the resulting atmosphere and diagnostics. Expect an unconverged status from this installation example.

## Supply your composition and opacities

Edit a copy of `solar-gs98.json` to set your composition. The file maps each element symbol to an absolute logarithmic abundance, log₁₀(Nₓ/Nₕ) + 12. Keep hydrogen at `"H": 12`. Supply the full mixture you intend to model because omitted elements are assigned a trace abundance.

For a line-blanketed calculation, provide a wavelength grid in vacuum Å, an atomic transition CSV via `--lines`, and one or more molecular tables via repeatable `--molecular-cross-section` options. The following command assumes you have prepared those input files. [Prepare atomic and molecular opacity inputs](advanced.md)

```sh
giantlayers flx_wavelengths.vac.gz mixture.json --teff 3500 --logg 1 --mass 1 --lines atomic-lines.csv --molecular-cross-section tio.json.gz --molecular-cross-section water.json.gz --stride 120 --steps 8 --verify-full --output-model atmosphere.mod > balance.csv
```

`--stride 120` fits on every 120th input wavelength, including the final endpoint. `--verify-full` reevaluates the final profile on the full supplied grid. The fitting stride and iteration budget are starting settings to adjust for your inputs.

## Check and export the result

- Inspect the maximum luminosity residual and the residuals by layer in the CSV. The default luminosity tolerance is 1%. Use `--require-flux-convergence` in scripts to return status 2 when the final profile exceeds that tolerance.
- Inspect molecular table coverage and compare the result at finer depth and wavelength resolution before using it quantitatively. The full-grid check measures luminosity closure for the supplied opacities. [Validation benchmarks](validation.md)
- Use `--output-model atmosphere.mod` for MARCS, `--output-atlas atmosphere.dat` for ATLAS9, or add `--atlas-format atlas12` for ATLAS12. You can request MARCS and ATLAS output in the same run. [Export options and format details](advanced.md#control-sampling-restarts-and-export)

## Plot the atmosphere

The development version can save the final temperature, gas and electron pressure, density, and luminosity residual profiles versus Rosseland optical depth. Install the plotting extra from a repository checkout with `python -m pip install '.[plot]'`, then add `--output-plot profiles.png` to your run. Use a `.pdf` or `.svg` filename for vector output. The plot shows the convergence status and the luminosity tolerance band.

[How GiantLayers works](how-it-works.md) · [Advanced usage](advanced.md)
