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

## Hot molecular opacity sensitivity

Repeat a zero-step restart with `--test-omit-hot-opacity TABLE` to remove a selected molecular table's absorption only where the local temperature exceeds its highest tabulated temperature. Supply the same path with `--molecular-cross-section`. Repeat the diagnostic flag to test several tables together. Use the same model and boundary inputs for the control and stress test. The option requires `--steps 0 --restart-csv PROFILE` and disables model export for the test. The control uses the usual clamped cross sections, with molecular abundances still computed at each layer's actual temperature.

For the 56-layer fiducial on 2,520 wavelength samples, omitting all four tables above 3,400 K shifted the layer luminosity residuals by as much as 4.50 percentage points relative to the sampled control. Tested individually, the largest shifts were 2.89 points for TiO, 1.72 for CO, 0.079 for H₂O, and 0.053 for VO. The control and joint stress test had maximum absolute residuals of 0.667% and 4.31%, respectively. These tests measure the leverage of the clamped opacity and prioritize TiO and CO for hotter tabulations. They do not estimate the error of the clamped cross sections against true high-temperature opacity.

On all 100,724 MARCS wavelengths, the joint omission increased the maximum absolute luminosity residual from 0.851% to 4.167%. The largest shift was 4.374 percentage points in layer 41 at 3,551 K. The control and stress test used identical temperatures and column masses. The surface cooling fraction changed from 1.615% to 1.304% of local LTE emission.

The [published TiO ToTo line list](https://doi.org/10.1093/mnras/stz1818) is designed for temperatures below 5,000 K. It can supply additional TiO cross-section slices beyond the 3,400 K limit of the imported grid.

## Recompute hot molecular opacity

The offline Doppler converter calculates LTE absorption directly from ExoMol states, Einstein A coefficients, and the tabulated partition function. It integrates each thermal and microturbulent Doppler profile over geometric wavenumber bins defined by the wavelength grid of a portable GiantLayers table. The output has the same species, isotope fraction, pressure axis, and wavelengths as the template. It recomputes all temperatures and repeats the same pressure-independent cross sections at each pressure. The generated table records hashes of all three ExoMol inputs.

The converter uses NumPy and a C compiler during table preparation. The atmosphere solver reads the resulting JSON file through `--molecular-cross-section` with its usual dependencies. A matched capped table isolates the effect of temperature coverage from differences between the Doppler calculation and the ExoMolOP reference cross sections.

For example, from a repository checkout with the TiO files and portable ExoMolOP template in the current directory:

```sh
cc -O3 -o exomol_doppler tools/exomol_doppler.c -lm
python tools/exomol_doppler.py tio-exomolop.json.gz \
  48Ti-16O__Toto.states.bz2 48Ti-16O__Toto.trans.bz2 48Ti-16O__Toto.pf \
  tio-extended.json.gz --kernel ./exomol_doppler --mass-amu 63.943846 \
  --extend-temperatures 3500 3600 3800 4000 4250 4500 4750 4900 \
  --capped-control tio-capped.json.gz
```

The TiO example uses the published [ExoMol ToTo files](https://www.exomol.com/data/molecules/TiO/48Ti-16O/Toto/) and a template ending at 3,400 K. The offline converter omits collisional and natural broadening. Use its tables to test the sensitivity to thermal populations and Doppler widths before adopting them in new atmosphere calculations.

## Temperature smoothness and layer structure

Use `--temperature-smoothness WEIGHT` to add a penalty for short-scale oscillations to a temperature fit. At every interior layer, the solver compares log temperature with a linear interpolation from the neighboring layers in the reference log Rosseland coordinate. It includes these weighted differences in the least-squares objective and in trial-step acceptance. Luminosity convergence is still measured from the unpenalized flux errors. The default weight of zero preserves the original objective. Test sensitivity to the weight and to the layer grid before interpreting a feature in the temperature structure.

## Compare native and MARCS flux samples

The MARCS GS98 archive contains the `.flx.gz` record for the fiducial atmosphere. Pair it with the published `flx_wavelengths.vac.gz` grid. Run GiantLayers on the complete grid with `--output-spectrum native.csv` and then compare broad wavelength integrals:

```sh
python -m compare_marcs_flux native.csv flx_wavelengths.vac.gz fiducial.flx.gz --sampled
```

Use repeatable `--band SHORT_NM LONG_NM` arguments for other ranges. The comparison requires an exact match between the native and MARCS wavelength arrays. These broad-band trapezoidal integrals compare statistical samples of flux per wavelength normalized to each model's reference radius. Line strengths and narrow-band integrals require dedicated spectrum synthesis. [MARCS flux format](https://marcs.astro.uu.se/docs.html)

## Full-grid results for the high-temperature opacity trial

The matched Doppler extensions use TiO ToTo through 4,900 K and CO Li2015 through 8,500 K, with the same recalculated samples up to 3,400 K in the capped controls. The published fiducial contains 100,724 wavelength samples. At the original GiantLayers temperature and hydrostatic structure, the extension moved individual luminosity residuals by up to 0.759 percentage points with CO alone, and 0.758 points with both molecules. Both largest shifts were inward at layer 50 (4,312 K). The two extension-induced changes differ by at most 0.007 percentage points across the depth grid. Four layers remain above the extended TiO grid, and all layers are covered in CO.

The joint capped and extended Doppler models evaluated at the original temperatures had peak luminosity errors of 5.231% and 5.484%, respectively. The model using the original ExoMolOP inputs had a 0.851% peak error at those same temperatures. Changing the opacity preparation and line inputs moves the luminosity more than adding hotter slices in this benchmark. Refit the temperatures before using the Doppler inputs in a model. [TiO ToTo line list](https://www.exomol.com/data/molecules/TiO/48Ti-16O/Toto/) · [CO Li2015 line list](https://exomol.com/data/molecules/CO/12C-16O/Li2015/)

With the original opacity inputs, applying temperature smoothness weight 0.25 reduced adjacent inward temperature decreases from 12 to 3, and the largest decrease fell from 203 to 138 K. The full-grid peak luminosity error was 1.756%, up from 0.851%. Surface cooling reached 18.86%, up from 1.61% of the local LTE emission. Over the same 47 MARCS Rosseland depths, the temperature RMS difference fell from 171 to 102 K. The original MARCS Fortran reader accepted the smoother export.

The original GiantLayers fiducial's emergent flux integrated over the full grid was 0.33% below the published MARCS flux. On the same wavelengths it was 9.13% lower at 400–700 nm, 12.51% lower at 700–1,000 nm, and 1.45% higher at 1,000–2,500 nm. These band integrals show the wavelength redistribution that remains even with good total luminosity closure. [MARCS GS98 flux archive](https://marcs.astro.uu.se/data.html)

Starting from the original 56-layer fiducial, one accepted correction with the extended Doppler inputs reduced the full-grid maximum luminosity error from 5.484% to 0.803%. The final RMS luminosity error was 0.320%, the surface cooling was 3.24%, and the bottom Rosseland depth was 110.7. Adjacent inward temperature decreases fell from 12 to seven. On the same 47 MARCS depths, the temperature RMS difference fell from 171 to 148 K. The MARCS Fortran reader accepted the refitted export.

For that refit, use the extended TiO and CO tables, the original ExoMolOP water and VO tables, the same revised atomic list and stellar parameters, and the following settings:

```sh
--restart-csv fiducial.csv --steps 7 --stride 120 --workers 6 --upper-knots 56 \
--jacobian-interval 4 --peak-error-weight 0.65 --outer-heating-weight 0.14 \
--temperature-smoothness 0.05 --full-grid-correction
```

The refitted Doppler model's integrated flux was 0.29% below MARCS over the full sampled grid. It was 14.22% lower at 400–700 nm, 10.92% lower at 700–1,000 nm, and 1.82% higher at 1,000–2,500 nm. The improved luminosity closure and temperature RMS therefore coexist with substantial optical-band errors. These measurements refer to the temperature refit before any surface-only correction.

Two surface-only corrections to the Doppler refit lowered its surface temperature from 1,883.3 to 1,868.3 K and reduced surface cooling from 3.24% to 1.11%. The full-grid maximum luminosity error remained 0.803%. A separate, weaker-regularization refit using the original opacities reached 0.857% maximum luminosity error, with nine adjacent inward temperature decreases and 8.96% surface cooling.

## Diatomic equilibrium pressure units

The bundled Barklem–Collet diatomic dissociation constants now convert the published pascal values to dyn/cm². This shifts their tabulated `log10(Kp)` by +1. The NASA-based polyatomic entries already use cgs units. All atmosphere trials described above used the earlier diatomic table and should be recomputed with the current table. [Barklem and Collet data](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A%2BA/588/A96)

At the official MARCS temperatures and gas pressures, median TiO differences for three compositions went from +0.93–0.94 dex to between −0.010 and −0.003 dex after the correction. Across 46 layers per model from log Rosseland depth −4 to +1, TiO RMS differences are 0.026 dex and CO RMS differences are 0.010–0.011 dex. Electron-pressure RMS differences are below 0.0004 dex. CN still differs by about −0.08 dex in the median. [Official MARCS composition grids](https://marcs.astro.uu.se/data.html)

Run the fixed-structure check against any complete MARCS atmosphere with:

```sh
python -m compare_marcs_chemistry reference.mod.gz --min-log-tau -4 --max-log-tau 1
```

The command uses the official element-by-element abundances, temperatures, and gas pressures. It reports the differences in electron, CO, CN, water, and TiO partial pressures before any atmospheric temperature or density fitting.

## Composition benchmarks after the unit correction

The four benchmarks use official spherical MARCS atmospheres at 3,500 K, log g = 1, one solar mass, and 2 km/s microturbulence. The alpha pair shares the MARCS metallicity label −0.50 and uses alpha labels +0.40 and +0.00. The heavy-CN mixture uses the official C-depleted, N-enhanced solar-metallicity record. The solver takes each reference's complete abundance vector and independently reconstructs hydrostatic structure on the same 56-layer GiantLayers mass grid. Its starting temperatures come from the earlier GiantLayers GS98 model. It uses revised atomic lines, extended TiO and CO Doppler tables, and ExoMolOP water and VO. CN participates in the equation of state. The radiative transfer in these trials includes no CN line opacity. [Official MARCS grids](https://marcs.astro.uu.se/data.html)

| Composition | C/O by number | Peak luminosity error | Temperature RMS at log Rosseland depth −2 to +1 |
| --- | --- | --- | --- |
| GS98 | 0.490 | 0.867% | 66 K |
| Heavy CN, after recalibration | 0.224 | 0.943% | 107 K |
| Alpha-enhanced | 0.214 | 17.685% | 186 K |
| Alpha-poor | 0.537 | 44.830% | 393 K |

All peak errors come from the final independent evaluation on the 100,724-wavelength MARCS grid. The heavy-CN run initially reached 1.176%. Recalibrating the 631-sample fit at that profile and allowing three more corrections brought it under 1%, with surface heating still 8.41% of local LTE emission. The two main alpha fits had their largest luminosity errors at layer 55. A separate alpha-poor run initialized from the official MARCS temperature profile interpolated in log gas pressure reached 12.030% peak error. Its Rosseland-depth overlap extended only to log depth +0.2. The MARCS Fortran reader accepted all six exported trial files.
