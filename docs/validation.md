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

A monochromatic check used radial extinction integrated with the same shell widths and trapezoid rule as the transfer kernel. On the 1,261-point spectra, 24.3% of the emergent radiative flux from the shallower 56-layer trial came from samples with bottom radial extinction optical depth below 10. The deeper 56-layer solution had no samples contributing flux below that threshold. Over the same wavelength samples, the shallower trial had 26.1% less integrated 400–700 nm flux than the deeper solution. These diagnostics quantify the sensitivity of the trial spectrum to the inner structure. The illustrative threshold of 10 is a radial extinction diagnostic and does not by itself certify a diffusion boundary.

The outermost layer of the initial flux-converged 56-layer solution had net radiative cooling equal to 31.7% of its LTE thermal emission. The solver CSV reports this diagnostic as `relative_net_heating`. With `--outer-heating-weight` enabled, `--full-grid-correction` now calibrates this term alongside the luminosity residual. As with the luminosity correction, the calibration is fixed during the sampled fit and the final full-grid transfer independently measures the actual residuals.

Recalibrating the second starting profile gave a full-grid maximum residual of 1.081% and left a 1,498 K bottom-temperature difference relative to the flux-converged 56-layer solution. Moving 18 nodes below the gray photospheric anchor and recalibrating separately gave a 1.48% maximum full-grid luminosity error and a bottom Rosseland optical depth of 81.3. These trials leave measurable sensitivity to both starting profile and depth placement.

Giving the outer-heating objective weight 0.25 in the full-grid-calibrated 56-layer fit reduced the measured maximum luminosity error from 0.896% to 0.367% and the outermost cooling fraction from 31.7% to 20.4%. The final bottom Rosseland depth was 106.4. Over the 30 layers with log Rosseland depth at or below log₁₀(0.03), the RMS net-heating fraction fell from 5.98% to 3.87%. One cold surface layer and 17 hot layers used molecular table boundary values. On the 51 MARCS depths covered by both exports, their temperature RMS differences were 104 K before and 120 K after the heating correction.

A second heating calibration at stride 40 (2,520 samples) and outer heating weight 0.35 gave 0.164% RMS and 0.608% maximum luminosity error on all 100,724 wavelengths. The final bottom Rosseland depth was 106.6. The outermost cooling fraction was 18.5% and the RMS net-heating fraction over the 30 layers at Rosseland depth at or below 0.03 was 3.47%. These full-grid measurements show that the flux-converged solution retains a local thermal imbalance. The original MARCS Fortran reader opened the resulting .mod file.

## Cool-surface opacity and heating

The cool-surface tests extend the same TiO, H₂O, VO, and CO ExoMolOP tables with the published 1,800, 1,900, and 2,000 K slices. The previous 2,200–3,400 K rows, 10–10,000 dyn/cm² pressure points, and native wavelengths are preserved. A local copy of each original TauREx HDF5 file can be converted with:

```sh
python3 -m exomol_import input.h5 output.json.gz --species TiO --isotope-fraction 0.7372 \
  --temperature-points 1800 1900 2000 2200 2400 2600 2800 3000 3200 3400 \
  --pressure-points 10 100 1000 10000
```

Use the appropriate species and isotopic fraction for each table. Explicit points match tabulated values, and unavailable choices are errors.

The heating fit calibrates the absorbed-to-emitted power ratio at the starting profile. If the sampled and full-grid normalized heating values are q_s and q_f, the calibration multiplies 1 + q_s by (1 + q_f)/(1 + q_s) evaluated at that profile. The multiplier stays fixed during the sampled fit. The final full-grid transfer reports the actual surface heating fraction and luminosity residuals.

The denser 5,038-point fit with the original opacity tables had a maximum full-grid luminosity error of 0.730% and a surface cooling fraction of 14.4%. At the same starting temperatures, adding the colder opacity data changed the measured surface cooling from 18.5% to 19.2%. After fitting with the colder data and the absorbed-to-emitted ratio calibration, the full-grid maximum luminosity error was 0.926% and the surface cooling fraction was 9.17%. The surface temperature was within the tabulated range and the hot layers still used the upper table edge.

Use `--equilibrate-surface` for one additional correction to the outermost temperature, or pass a number to repeat it. GiantLayers measures the full-grid mean intensity, solves local thermal balance while holding that radiation field and surface gas pressure fixed, then reconstructs and evaluates the atmosphere on the complete wavelength grid. It accepts the correction when the measured local imbalance decreases and the maximum luminosity error stays within tolerance or improves from an already unconverged profile. Other layer temperatures and hydrostatic boundary parameters stay fixed during this step. A first application took the fiducial surface from 1,875 K to 1,850 K, reducing measured cooling to 5.36% and maximum luminosity error to 0.886%.

Two more accepted surface-only corrections brought the surface to 1,828 K and the measured surface cooling to 1.61% of its LTE emission. Across the 30 layers at Rosseland depth at or below 0.03, the RMS normalized heating residual was 0.81% and the largest was 2.84%. The final full-grid luminosity errors were 0.223% RMS and 0.851% maximum. Bottom Rosseland optical depth was 106.7. The MARCS Fortran reader opened the export.

Relative to the model with the same extended opacity data before the surface-only corrections, the sampled 400–700 nm band flux changed by +0.077% and the 1–2.5 μm band flux by −0.135% on 2,520 wavelengths. At the same 47 MARCS Rosseland-depth comparison points, the temperature RMS difference increased from 144 to 171 K. These checks establish a better thermal balance for the supplied inputs. The molecular opacities still use the hot table edge in 17 layers and the low-pressure edge in two layers.
