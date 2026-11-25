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

Using all 56 independent temperature controls and calibrating the 1,261-point fit against the full wavelength grid reduced the 56-layer trial's maximum error from 7.83% to 0.896% in one accepted correction. The final full-grid RMS error was 0.263%, with −0.080% at the surface. The 32- and 56-layer calculations both met the 1% luminosity criterion and still differed by up to 775 K at common radial column mass. The 56-layer temperature RMS difference from MARCS was 104 K over 54 shared Rosseland depths. Luminosity closure and structural resolution are separate validation criteria.

On the shared optical-depth range, the two flux-converged solutions differ by 145 K RMS. From log Rosseland depth −3 to 0, the difference is 126 K RMS. A separate emergent-spectrum check on the same 1,261 wavelength samples found 13.6% more integrated 400–700 nm flux and 1.0% less integrated 1–2.5 μm flux in the 56-layer model. These are rectangular-band integration diagnostics on the sampled spectra. They show why bolometric flux alone is insufficient to validate the predicted spectrum.

Holding the 56-layer temperatures fixed on the 1,261-point wavelength grid, increasing the core ray count from 24 to 48 changed layer luminosity errors by at most 0.076 percentage points. Increasing from 48 to 96 changed them by at most 0.023 percentage points. These angular sensitivity checks are smaller than the changes produced by remapping the temperature structure.

To repeat the full-grid correction of the 56-layer trial, use the same stellar, opacity, and boundary inputs with `--layers 56 --restart-csv depth-audit.csv --stride 80 --steps 12 --upper-knots 56 --jacobian-interval 4 --peak-error-weight 0.8 --outer-heating-weight 0.005 --workers 6 --full-grid-correction`. The temperature correction itself uses the sampled grid and the final check uses the full grid.

A second starting profile was made by interpolating the 32-layer solution to the original 56-layer grid. After 16 corrections, that trial had 1.085% maximum luminosity error on the full wavelength grid. Its bottom reached Rosseland optical depth 16.6 and carried 5.3% of the target luminosity by convection. The 56-layer calculation with 0.896% maximum error reached Rosseland optical depth 104.7 and carried 78.8% by convection at the bottom. Their bottom temperatures differed by 1,504 K at the same radial column mass. These runs establish strong sensitivity to the starting profile and the inner thermal structure.

GiantLayers reports the final bottom Rosseland optical depth. Set `--min-bottom-rosseland-tau` to impose a minimum useful for a specific calculation and return status 2 when it is missed. Choose and verify the inner mass boundary for the intended model. Monochromatic optical depth in transparent windows requires a separate check.
