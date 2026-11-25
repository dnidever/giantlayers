# Validation

## 3,500 K giant

The main reference is the spherical MARCS GS98 giant at 3,500 K, log g = 1, one solar mass, and 2 km/s microturbulence. GiantLayers starts from an analytic temperature profile and independent radial column-mass grid. The reference supplies the comparison structure, abundance pattern, and wavelength axis. [MARCS model downloads](https://marcs.astro.uu.se/data.html) · [MARCS atmosphere methods](https://arxiv.org/abs/0805.0554)

The trial uses a revised MARCS atomic transition subset and ExoMolOP cross sections for TiO, H₂O, VO, and CO. Run the atomic transition converter and the ExoMolOP converter described in the [Quickstart](quickstart.md). The example abundance file supplies the GS98 mixture.

```sh
cc -O3 -std=c99 -shared -fPIC spherical_kernel.c -o giantlayers-transfer.so -lm
giantlayers flx_wavelengths.vac.gz examples/solar-gs98.json \
  --teff 3500 --logg 1 --mass 1 --layers 32 --core-rays 24 \
  --stride 400 --refine-stride 160 --steps 8 --upper-knots 10 \
  --kernel ./giantlayers-transfer.so --lines revised-atomic-lines.csv \
  --molecular-cross-section tio.json.gz --molecular-cross-section water.json.gz \
  --molecular-cross-section vo.json.gz --molecular-cross-section co.json.gz \
  --checkpoint-csv checkpoint.csv --output-model trial.mod.gz > trial.csv
```

The startup and refinement grids contain 253 and 631 points, respectively. Reuse the same opacities and structure arguments with `--restart-csv trial.csv --stride 1 --steps 0 --workers 4` to evaluate the profile at all 100,724 MARCS wavelengths. Use the per-layer luminosity errors to distinguish convergence on the startup samples from convergence on the full grid.

In the first independent 631-point trial, the RMS luminosity error was 5.00% and the largest layer error was 14.63%. Restarting the 631-point calculation with `--upper-knots 13 --adaptive-knots 4 --peak-error-weight 1` reduced these to 2.26% and 6.83% after three accepted steps. Evaluating that refined temperature profile on the complete grid gave 2.84% RMS, 6.95% maximum, +1.28% at the surface, and +0.17% at the bottom. Compared at common Rosseland depths, temperatures differed from MARCS by 120 K RMS over 55 reference layers. This run used 32 radial layers and the four molecular tables above.

A subsequent fit used 1,261 wavelength samples, all 32 temperature controls, and `--full-grid-correction`. The final transfer on all 100,724 samples gave 0.52% RMS luminosity error, 0.863% maximum, −0.25% at the surface, and +0.06% at the innermost layer. The temperature RMS difference from MARCS was 125 K over 55 shared reference layers. This establishes luminosity convergence for that 32-layer calculation. Ten layers used the upper temperature edge of the molecular opacity tables. The 56-layer refinement and the carbon-rich trials require further convergence work.

Compare temperature and pressures over the common optical depths with `python3 -m compare_atmospheres marcs-reference.mod.gz trial.mod.gz`. Confirm that each dataset covers the line-forming layers and check the results on finer depth and angle grids.

## Analytic and software checks

Run `python3 -m unittest discover -v`. The tests cover the constant-density hydrostatic sphere, the effect of radiative acceleration, scattering transfer, temperature corrections, chemical pressure and charge conservation, opacity interpolation, restart interpolation, and MARCS-format export.

## Depth-resolution audit

Hold the stellar inputs, opacity data, wavelength grid, angular quadrature, and boundary options fixed. Solve at two depth resolutions and compare the solver CSVs:

```sh
python3 -m audit_resolution coarse.csv fine.csv --require-pass > depth-audit.json
```

The JSON reports luminosity closure for each calculation and identifies the layer and column mass with its largest error. It also measures the maximum temperature and log gas-pressure differences on the union of the two depth grids over their shared radial column-mass interval. Temperature and log pressure are interpolated linearly in column mass. Choose structure tolerances for the intended science with `--temperature-tolerance` in kelvin and `--pressure-tolerance` in dex. With `--require-pass`, status 2 identifies a failed configured check. Use the CLI and input files to establish that the two CSVs share the intended physical and numerical settings.

For the 3,500 K fiducial, a restarted 56-layer trial was fitted on 631 wavelengths and evaluated at fixed temperature on all 100,724 wavelengths. Both full-grid evaluations used the same radial mass extent, optical inputs, angular quadrature, and original boundary options. The 56-layer trial gave 1.89% RMS and 7.83% maximum luminosity error. Its worst layer was at a radial column mass of 540.2 g/cm², just above the layers carrying appreciable convective flux. The maximum temperature difference between the profiles was 775 K at 589.0 g/cm². The two calculations have yet to establish depth convergence. The 56-layer profile used the upper temperature edge of the molecular opacity tables in 18 layers.

A fixed-temperature remapping check used the same 1,261 wavelength samples on both grids, with the 32-layer temperature profile interpolated onto 56 layers and no temperature corrections. The maximum luminosity error rose from 1.96% to 94.54%. At the worst layer (540.2 g/cm²), the mixing-length convective flux was zero, so this discrepancy appeared directly in the radiative flux. The grid change includes both the new transfer cells and the hydrostatic reconstruction of the interpolated profile. This test measures sensitivity to that combined change.

Use `--deep-layers 18` with `--layers 56` to place 18 grid points below the gray photospheric anchor while keeping the anchor column and both mass endpoints fixed. Use `--interpolate-restart` when continuing from a differently distributed grid. Compare the resulting temperature structure after separately checking luminosity closure on each grid.
