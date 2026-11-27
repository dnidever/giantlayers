# Advanced usage

The [Quickstart](quickstart.md) shows how to build a spherical atmosphere from stellar parameters. Use these workflows to prepare opacity inputs, compare against reference atmospheres, or isolate individual parts of the solve. Commands that take a MARCS model as input are diagnostics or comparisons on a supplied structure.

## Run with a two-group opacity table

The nonlinear prototype can read both mass opacities from a CSV file.

```sh
python -m giantlayers.nonlinear_atmosphere --opacity-table opacities.csv > atmosphere.csv
```

Use the following header and provide one row for each temperature–pressure pair on a complete rectangular grid with at least two values on each axis. All input values must be positive and finite. Opacities are in cm²/g.

```csv
temperature_K,pressure_dyn_cm2,group1_opacity_cm2_g,group2_opacity_cm2_g
```

The solver interpolates in log temperature, log pressure, and log opacity. Pressures below the lowest grid point, including the vacuum surface, use the lowest-pressure row. A temperature outside the table or a pressure above its upper limit stops the run. Make sure the table covers the temperatures encountered during iteration as well as the final profile.

The two groups retain equal Planck weights. Tables supplied here should contain the group means appropriate to that simplified model.

## Compare with a MARCS atmosphere

Use a matching MARCS model structure, continuum table, and sampled surface flux from the same grid point. Keep the published MARCS wavelength grid alongside the flux file. [MARCS downloads](https://marcs.astro.uu.se/data.html) · [File formats and wavelength grid](https://marcs.astro.uu.se/docs.html)

```sh
python -m giantlayers.odf_transfer p00big2.bdf --marcs-model model.mod.gz \
  --marcs-continuum model.opa.gz > bands.csv
python -m giantlayers.compare_marcs_flux bands.csv flx_wavelengths.vac.gz model.flx.gz > comparison.csv
```

GiantLayers uses the MARCS true continuum absorption with ATLAS9 line opacity on the fixed temperature and pressure structure. The transfer uses plane-parallel two-stream geometry. Add `--continuum-scattering` to include isotropic coherent scattering from the MARCS table. The comparison integrates published flux samples only over bands fully covered by the MARCS wavelength grid.

### Inspect radiative heating

Add `--heating-output heating.csv` to write the heating profile while keeping the band flux CSV on standard output. The profile includes column mass, local net heating per gram, integrated heating in each depth cell, and net outward radiative flux. A positive heating rate means the matter absorbs more energy than it emits. The diagnostic uses the conservative two-stream solver and prints the closure error between the summed cell heating and the change in radiative flux across the atmosphere.

## Write a MARCS atmosphere file

The exporter writes a complete MARCS model as plain text or gzip. The command below reads an existing model and re-exports its header, 92 elemental abundances, structure, thermodynamics, and partial pressures.

```sh
python -m giantlayers.marcs_writer input.mod.gz output.mod
```

From Python, call `marcs_writer.write_document(model, path)` with a complete `MarcsDocument`. Supply every table and the abundance pattern for the same physical model. The gray and two-group CSV outputs contain only a subset of the fields required for export.

## Calculate opacities for a custom composition

Run the native LTE opacity solver with Python 3.9 or later. It uses the Python standard library and the thermochemistry data included in the checkout.

Save a composition as JSON. The optional metal offset shifts every element heavier than helium in dex relative to the input MARCS atmosphere. Individual values override the shifted abundances on the absolute scale where hydrogen is 12.

```json
{
  "metallicity_offset": -0.3,
  "log_epsilon": {"C": 8.72, "N": 7.92, "O": 8.83}
}
```

Supply a transition CSV with vacuum wavelength in Å, species, lower excitation energy in eV, and log gf. Use element symbols for neutral atoms, append `+` for singly ionized atoms, and use formulas such as `CO` for molecules. An optional `damping_s` column supplies the total Lorentz damping rate in s⁻¹. With no damping rate, the profile is Doppler broadened.

```csv
wavelength_vacuum_A,species,lower_energy_eV,log_gf
15570.299575,CO,6.996,-6.129123
```

```sh
python -m giantlayers.lte_opacity model.mod.gz transitions.csv mixture.json 15570 15585 0.05 opacity.csv --vmic 2
```

The reference MARCS structure supplies the temperature and gas pressure. The solver recomputes the mass density, electron density and species populations at each layer. The output contains sampled line absorption, continuum absorption, and continuum scattering in cm²/g. For isotopologue lines, include the isotopic fraction in the oscillator strength used in the transition CSV.

The bundled data cover 1,000–10,000 K and include neutral atoms, first ions, H⁻, 20 neutral diatomic molecules, and H₂O. Water equilibrium uses the NASA Glenn data through 6,000 K. Continuum terms include H⁻ bound-free and free-free absorption, hydrogenic H I bound-free absorption, electron–ion free-free absorption, and Thomson scattering. The hydrogenic continuum terms use unit Gaunt factors. Additional polyatomic and charged molecules, higher ionization stages, metal photoionization, and calculated collisional line damping are future extensions.

[Atomic and diatomic thermochemistry data](https://doi.org/10.1051/0004-6361/201526961) · [H⁻ continuum data in the ATLAS9 HMINOP routine](https://wwwuser.oats.inaf.it/fiorella.castelli/sources/atlas9g/atlas9mem.for) · [APOGEE DR17 line data reference](https://ui.adsabs.harvard.edu/abs/2021AJ....161..254S/abstract)

### Evaluate radiative heating with native opacities

Provide a text file containing increasing vacuum wavelengths in Å. The transfer uses the same grid for all layers and integrates in frequency using trapezoidal weights. Include broad wavelength coverage and enough samples to resolve the lines when interpreting the result as a bolometric energy balance.

```sh
python -m giantlayers.native_transfer model.mod.gz transitions.csv mixture.json wavelengths.txt heating.csv
```

Add `--adjust-layers 6 --steps 3` to try up to three damped temperature corrections in the outermost six layers. Gas pressure and column mass stay fixed. Each trial recomputes LTE equilibrium, opacity, and transfer, and is accepted only if it reduces the maximum outer-layer net heating relative to thermal emission. Deeper layers supply a fixed lower structure.

The CSV reports initial and corrected temperature, radiative flux, local heating, and integrated cell heating. These corrections are a radiative test on a fixed pressure structure. A converged atmosphere also requires adequate opacity coverage, a flux constraint, and treatment of convective energy transport.

### Use the fiducial cool giant

Use the spherical MARCS GS98 model at 3,500 K, log g 1.0, one solar mass, and 2 km/s microturbulence for the cool-giant benchmark. The model header reports [Fe/H] = +0.05 and [α/Fe] = +0.11 on the MARCS reference scale.

Extract this member from the MARCS GS98 structure archive. Use an empty JSON object for the published abundance pattern. [MARCS GS98 archive](https://marcs.astro.uu.se/data.html)

```text
s3500_g+1.0_m1.0_t02_gs_z+0.05_a+0.11_c+0.08_n+0.09_o+0.12_r+0.00_s+0.00.mod.gz
```

The energy-balance CSV includes the radiative flux, published MARCS convective contribution, and total flux error at each depth. The reference flux is scaled by the inverse square of the model radius for spherical models. The transfer still uses plane-parallel geometry. Add `--flux-guard` to a temperature trial to reject steps that reduce the local heating imbalance while worsening the surface flux error. The diagnostic keeps the MARCS convective flux fixed.

### Add tabulated TiO and water absorption

Both native opacity and transfer commands accept repeatable `--molecular-cross-section FILE` options. A portable JSON or gzip-compressed JSON table has the following fields.

| Field | Contents |
| --- | --- |
| `species` | Chemistry species such as `TiO` or `H2O` |
| `isotopologue` | Optional data identifier |
| `isotope_fraction` | Fraction of the parent species represented by the table. Defaults to 1. |
| `temperature_K`, `pressure_dyn_cm2`, `wavelength_vacuum_A` | Increasing grid axes |
| `cross_section_cm2` | Three-dimensional array ordered by temperature, pressure, and wavelength. Values are in cm² per isotopologue. |

The code scales each cross section by the molecule population and isotope fraction. It interpolates in temperature, log pressure, and wavelength, using endpoint temperature or pressure values outside those axes. Wavelengths outside the table contribute zero absorption. Represent each isotopologue once by choosing either its cross section or its transitions on overlapping wavelengths.

The first fiducial experiment sampled ExoMol ToTo and POKAZATEL cross sections over 0.3–50 μm on a 2,200–3,400 K grid. Deeper, hotter layers use the hottest tabulated cross section for each molecule that remains in the chemistry calculation. Use the same wavelength grid when comparing runs. [TiO ToTo](https://exomol.com/data/data-types/opacity/TiO/48Ti-16O/Toto/) · [H₂O POKAZATEL](https://zenodo.org/records/5722861) · [NASA Glenn thermochemistry](https://ntrs.nasa.gov/citations/20020085330)

### Check wavelength sampling

Run the sampling benchmark with cross-section tables on a common wavelength axis. It keeps the reference atmosphere fixed and reports surface flux, the maximum total-flux error with depth, and the heating/flux closure error on each nested grid. Supply any additional wavelengths outside the tabulated range or in line windows through `--supplemental-wavelengths`.

```sh
python -m giantlayers.benchmark_sampling model.mod.gz transitions.csv mixture.json tio.json.gz h2o.json.gz --strides 4 2 1 > sampling.csv
```

For expanded optical and infrared coverage, include ExoMol HyVO and Li2015 CO tables on the same grid. If the CO table replaces the explicit CO transitions, remove those transitions from the input CSV. A header-only transition CSV is valid when using cross sections alone. These files add molecular opacity over the common sampled range and inherit the temperature and pressure coverage of the source tables. [VO HyVO](https://exomol.com/data/molecules/VO/51V-16O/HyVO/) · [CO Li2015](https://www.exomol.com/data/molecules/CO/12C-16O/Li2015/)

### Run on the published MARCS wavelength grid

Pass the published vacuum grid directly to native transfer with `--quadrature wavelength`. The reader accepts the original gzip file, and the transfer processes the opacity in wavelength blocks while retaining quadrature weights across block boundaries.

```sh
python -m giantlayers.native_transfer model.mod.gz transitions.csv mixture.json flx_wavelengths.vac.gz full-grid.csv --quadrature wavelength --molecular-cross-section tio.json.gz
```

To compare coarse startup sampling with the full published grid, pass the same input file to the sampling benchmark. Stride 10 retains every tenth sample and both endpoints. Stride 1 uses all 100,724 wavelengths.

```sh
python -m giantlayers.benchmark_sampling model.mod.gz transitions.csv mixture.json tio.json.gz --reference-wavelengths flx_wavelengths.vac.gz --strides 10 1 --quadrature wavelength > marcs-sampling.csv
```

Use the native-resolution opacity tables when comparing the grids. The transfer evaluates the available opacities at the MARCS sample positions. This matches the published wavelength axis and the wavelength integration convention. The geometry, opacity inventory, and convection model also affect how closely the calculated atmosphere reproduces MARCS. [MARCS wavelength download](https://marcs.astro.uu.se/docs.html) · [MARCS methods](https://arxiv.org/abs/0805.0554)

### Compare spherical and plane-parallel surface transfer

Convert the published MARCS oscillator-strength revisions from VALD air wavelengths into the native vacuum-wavelength transition format. These files provide a subset of the atomic lines used by MARCS. The converted transitions use Doppler broadening. [MARCS auxiliary atomic data](https://marcs.astro.uu.se/docs.html#Auxiliary)

```sh
python -m giantlayers.marcs_atomic_lines astroph_gfVALD1.dat astroph_gfVALD2.dat > revised-atomic-lines.csv
```

An atomic transition CSV can supply `damping_s` for a density-independent Lorentz damping rate in s⁻¹. Add `electron_damping_cm3_s` and `neutral_h_damping_cm3_s` for rate coefficients at 10,000 K. The optional `electron_damping_exponent` and `neutral_h_damping_exponent` multiply those coefficients by a power of `T / 10,000 K`. GiantLayers multiplies each coefficient by its local perturber number density and adds the three rates before computing the line profile. Missing coefficients and exponents default to zero. Convert logarithmic coefficients to linear units when preparing the CSV.

Convert a Kurucz wavelength-format GFALL or per-ion `.pos` list with `python -m giantlayers.kurucz_atomic_lines gf2600.pos --wavelength-nm 400 700 > iron.csv`. The converter accepts neutral and singly ionized atoms and uses the rates supplied in the line records. It assumes temperature-independent electron coefficients and a neutral H exponent of 0.3 for this conversion. Use `--vacuum` for a list supplied entirely in vacuum or `--doppler-only` for a coverage-only comparison. When combining line lists, replace overlapping transitions to avoid counting them twice.

Use a spherical model and the same opacity inputs to compare surface transfer in both geometries. Pass several `--molecular-cross-section` options to include multiple absorbers. The spherical solver traces rays through atmospheric shells and integrates the emergent intensity over the apparent disk. The plane-parallel case uses the same surface angles and quadrature weights. Both modes use LTE true absorption and fixed temperatures and gas pressures for this benchmark.

```sh
python -m giantlayers.spherical_transfer model.mod.gz revised-atomic-lines.csv mixture.json flx_wavelengths.vac.gz --stride 50 --core-rays 24 --molecular-cross-section tio.json.gz > spherical-comparison.csv
```

The CSV contains monochromatic surface fluxes per Hz and frequency weights. The printed spherical ratio normalizes the flux to the MARCS reference radius. The planar ratio uses the corresponding target flux at the outer boundary. Reduce `--stride` and increase `--core-rays` to test sensitivity to the wavelength and ray grids.

### Compute spherical heating with scattering

Use the same spherical model and opacity inputs with the depth-resolved solver. It iterates the isotropic coherent continuum scattering source on the shared impact-parameter rays. Net heating is calculated from true absorption and the difference between the mean intensity and the LTE Planck function. Use `--absorption-only` to isolate the effect of scattering.

```sh
python -m giantlayers.spherical_heating model.mod.gz revised-atomic-lines.csv mixture.json flx_wavelengths.vac.gz --stride 100 --core-rays 24 --molecular-cross-section tio.json.gz > spherical-heating.csv
```

The CSV reports integrated radiative flux, luminosity, and pointwise true-absorption heating at each layer. The shell columns give the mean net heating and energy gain between the current layer and the next deeper layer, calculated from the luminosity difference and shell mass. Positive heating means the gas gains energy. The shell values balance the luminosity change by construction. Standard error reports how far a separate trapezoidal integration of the pointwise rates departs from that budget. Use that difference to assess whether pointwise heating is suitable for a temperature update. The MARCS temperature, pressure, and convection profiles are held fixed.

### Try spherical temperature corrections

Add `--adjust-layers 3 --steps 2` to the spherical heating command to try two damped temperature corrections in the outer three layers. Each trial recomputes LTE chemistry, opacities, and spherical scattering. The correction uses finite differences to reduce the change in radiative plus convective luminosity across the selected shells. The deeper temperature profile, radii, gas pressure and published MARCS convective contribution stay fixed.

The solver reports the maximum outer-shell luminosity residual and the surface total-flux error separately. The correction targets local shell balance. Fitting the absolute emergent flux requires a separate constraint. The CSV includes the initial and corrected temperature.

### Fit the absolute luminosity and check the full wavelength grid

Add `--luminosity-target` to fit the MARCS surface luminosity alongside the outer-shell luminosity differences. The correction adds one smooth temperature mode centered at Rosseland optical depth one. It tapers to zero at log optical depths −1.5 and +1.5. Gas pressure, radii, and the reference convective flux stay fixed.

With a coarse startup stride, `--verify-full` reevaluates the accepted profile on every input wavelength and writes the full-grid diagnostics. The verification holds the accepted temperatures fixed. For larger runs, the optional C99 kernel accelerates the ray solve while chemistry and opacity remain in Python.

The main `giantlayers` command automatically uses the C kernel when installed from a binary wheel. The commands below compile a custom library from a source checkout for the heating diagnostic, which takes an explicit `--kernel` path.

```sh
cc -O3 -std=c99 -shared -fPIC src/giantlayers/spherical_kernel.c -o spherical_kernel.so -lm
python -m giantlayers.spherical_heating model.mod.gz revised-atomic-lines.csv mixture.json flx_wavelengths.vac.gz --stride 400 --adjust-layers 3 --steps 5 --luminosity-target --verify-full --kernel ./spherical_kernel.so --molecular-cross-section tio.json.gz > luminosity-check.csv
```

### Solve pressure and radius on a fixed temperature profile

Run the spherical hydrostatic benchmark on a MARCS reference model.

```sh
python -m giantlayers.spherical_hydrostatic model.mod.gz > hydrostatic.csv
```

The benchmark takes the temperature, composition, and outer gas pressure from MARCS. It integrates the published density over radial distance to define a comparison column-mass grid. GiantLayers then computes its own gas pressure, density, and radius at each depth using spherical gravity and its native LTE equation of state. The reference radius is set by the input mass and gravity at the layer with Rosseland optical depth one. This version solves gas-pressure balance. Radiation and turbulent pressure are outside that equation.

The Python `solve` function also accepts a radial column-mass grid and boundary values supplied independently of MARCS.

### Compare transfer after hydrostatic reconstruction

Add `--hydrostatic` to the spherical heating command to rebuild pressure, electron pressure, density, and radii before evaluating transfer. The reconstruction and transfer use the same input abundance pattern. The temperatures and reference convective fractions stay fixed. Run the command with and without the option using identical wavelength, line, and cross-section inputs to isolate the effect of hydrostatic reconstruction.

The reconstructed CSV reports radial column mass measured from the upper boundary. The original comparison reports the model's published column mass. Match rows by layer index. The reconstructed radius at Rosseland depth one follows the input mass and gravity, so its reference luminosity normalization uses that radius.

### Couple the hydrostatic and temperature updates

Combine `--hydrostatic` with `--adjust-layers` to solve hydrostatic balance before each temperature trial, including the finite-difference perturbations and the line-search evaluations. The original radial mass grid, surface gas pressure, stellar mass, and gravity define the boundary throughout the run. Each trial recomputes pressure, electron pressure, density, radius, opacities, and scattering. With `--luminosity-target`, it also fits the target at the reconstructed reference radius.

```sh
python -m giantlayers.spherical_heating model.mod.gz revised-atomic-lines.csv mixture.json flx_wavelengths.vac.gz --hydrostatic --adjust-layers 3 --steps 2 --luminosity-target --stride 400 --molecular-cross-section tio.json.gz > coupled-trial.csv
```

The output includes the updated temperatures and hydrostatic structure. The Rosseland layer labels and convective luminosity fractions are carried over from MARCS for this benchmark.

### Calculate mixing-length convection

Run the Henyey mixing-length calculation on the fiducial pressure and temperature profile. GiantLayers computes the logarithmic temperature gradient, estimates the heat capacity and thermal expansion from its native LTE chemistry, and reports the convective velocity and flux beside the MARCS values. The default benchmark uses the published MARCS Rosseland opacities to isolate the convection calculation. [MARCS mixing-length equations](https://arxiv.org/abs/0805.0554)

```sh
python -m giantlayers.mixing_length model.mod.gz > convection.csv
```

Pass `--wavelengths flx_wavelengths.vac.gz` to compute the Rosseland mean from GiantLayers opacities. Add the atomic transition list with `--lines` and molecular tables with repeatable `--molecular-cross-section` arguments. Use `--hydrostatic` to evaluate convection on the reconstructed pressure and radius profile. These commands diagnose convection on a fixed temperature structure.

In the development version, use `--lte-radiation-pressure` with the standalone `giantlayers.mixing_length` diagnostic to test a total-pressure gradient at Rosseland optical depths of one and greater. In those layers, it evaluates the thermal expansion and enthalpy derivatives at constant gas plus LTE radiation pressure. The diagnostic uses the isotropic LTE approximation for radiation pressure. This option is for fixed-structure comparisons. The coupled atmosphere command uses the gas-pressure formulation.

### Recompute convection during the coupled solve

Add `--native-convection` to use the Henyey flux in place of the benchmark's published convection. Each trial computes a Rosseland mean from the same line and continuum opacities used by spherical transfer. It then recomputes the convective flux from the trial temperature and gas pressure gradients and the native LTE equation of state.

```sh
python -m giantlayers.spherical_heating model.mod.gz revised-atomic-lines.csv mixture.json flx_wavelengths.vac.gz --hydrostatic --native-convection --adjust-layers 3 --adjust-deep 2 --steps 1 --luminosity-target --stride 400 --molecular-cross-section tio.json.gz > coupled-convection.csv
```

The outer temperature variables target adjacent-shell luminosity balance. The photospheric mode targets absolute surface luminosity, and `--adjust-deep` varies the deepest temperatures individually to target total luminosity at those layers. The solver reports these errors separately. This trial adjusts selected depths and retains the original radial mass grid and outer pressure.


### Fit the total luminosity at every layer

Use `--global-balance` to solve for the temperature profile with native spherical transfer, hydrostatic balance, and Henyey mixing-length convection. The target is the same total luminosity at every radius. Each trial recomputes the gas structure, chemistry, wavelength-dependent opacities, and radiative and convective fluxes on the input model's radial mass grid. The input mass, gravity, and surface gas pressure set the boundaries.

GiantLayers identifies the contiguous region extending inward to the model's deepest layer where convection carries at least 0.01% of the target luminosity. Each layer in that region gets its own temperature variable. Above it, fractional temperature corrections are interpolated in the input model's log Rosseland depth between five control points by default. Change that number with `--upper-knots`. The solver evaluates the total luminosity error at every layer. Its line search accepts a step when the RMS error plus one quarter of the maximum absolute error decreases.

```sh
python -m giantlayers.spherical_heating model.mod.gz revised-atomic-lines.csv mixture.json flx_wavelengths.vac.gz --global-balance --stride 400 --core-rays 24 --steps 1 --kernel ./spherical_kernel.so --molecular-cross-section tio.json.gz > global-luminosity.csv
```

The CSV contains the corrected structure and the relative total-flux error at every layer. Increase `--steps` to attempt additional corrections. For a wavelength-resolution check, add `--verify-full` to reevaluate the resulting profile across the full input wavelength grid.

## Control sampling, restarts, and export

Use `--stride` to subsample the wavelength grid and repeat `--refine-stride` to solve on denser subsets. `--steps` sets the correction budget per grid. For long runs, `--checkpoint-csv` saves the accepted profile at each step. Use `--restart-csv` with the same mass-grid parameters to resume and `--steps 0` to evaluate the profile without temperature corrections.

Use the experimental `--heating-transition-tau 0.01` option to fit local radiative balance at shallower Rosseland depths and total luminosity deeper down. This option replaces `--outer-heating-weight` and retains a surface luminosity constraint. The partition is set at the start of each fitting stage. The final report gives both the luminosity error and the peak normalized heating over the final optically thin layers. With a 100,724-point input grid, `--stride 120 --refine-stride 12 --full-grid-correction` uses about 840 and 8,400 points for successive temperature fits and verifies the result on the full grid. Add a final `--refine-stride 1` to fit on the full grid when that compute cost is practical. [MARCS energy-balance methods](https://www.aanda.org/articles/aa/pdf/2008/30/aa09724-08.pdf)

The MARCS-format export contains native thermochemical partial pressures, Rosseland and 5000 Å optical depth, convection, and the calculated gas structure. By default the radius is anchored at the fixed mass coordinate selected from the seed. Use `--photospheric-anchor` to update the reference radius coordinate to the computed Rosseland photosphere. The export uses an isotropic LTE estimate for radiation pressure and zero turbulent pressure. The H₂⁺ column uses the absent-species sentinel, and the Other column sums the unlisted supported molecules. Files with large residuals are labeled unconverged.

The ATLAS export includes the input abundance mixture and all native depth points. The title identifies the spherical solution and its luminosity convergence status. The deck stores column mass and omits shell radii. Consumers that require a plane-parallel atmosphere or exactly 72 layers need to account for those assumptions. The `PRADK` field estimates surface radiation pressure from the emergent radiative flux assuming isotropic intensity over the outgoing hemisphere. [Published ATLAS9 grid files](https://wwwuser.oats.inaf.it/fiorella.castelli/grids/gridp00k2odfnew/ap00k2tab.html)

Use `--output-atlas trial.dat --atlas-format atlas12` for the Castelli ATLAS12 punch layout. It uses the same `READ DECK6` depth convention with nine columns and includes a separate `ABUNDANCE TABLE` with the model's uniform mixture. The reference ATLAS12 reader supports at most 72 depth points. Give ATLAS12 the uncompressed file with `READ PUNCH`. The export describes the GiantLayers atmosphere and retains its native depth grid. [ATLAS12 source and input examples](https://wwwuser.oats.inaf.it/fiorella.castelli/sources/atlas12.html)

## Diagnose radiative equilibrium in the outer layers

The native CSV includes the Rosseland mean and the ratio of net radiative heating to LTE thermal emission at each layer. An optically thin layer can have substantial local heating while barely changing the total luminosity. Set `--outer-heating-weight 0.1` to include that local balance in the trial objective. The weight falls exponentially with the initial Rosseland depth on a scale of 0.03. The luminosity status printed by the command is based on total flux. Check the heating column as well when fitting both terms.

The continuum includes coherent electron scattering and Rayleigh scattering by H I, H₂, and He I. The Rayleigh terms use long-wavelength expansions on the red side of the first resonance. The transfer solves those opacities with an isotropic coherent source. When accelerated scattering iterations reach their budget, GiantLayers solves the linear source equation directly. [H I cross section](https://adsabs.harvard.edu/pdf/1964ApJ...140.1381S) · [H₂ cross section](https://adsabs.harvard.edu/pdf/1962apj...136..690d) · [He I cross section](https://arxiv.org/abs/1601.01005)

## Convert ExoMolOP opacity tables

Download a TauREx HDF5 cross-section table for the desired isotopologue and convert it into the JSON format used by the solver. The optional h5py dependency is used by the converter. The transfer command reads the resulting table using Python's standard library.

```sh
python -m pip install --pre 'giantlayers[exomol]'
python -m giantlayers.exomol_import TiO.xsec.TauREx.h5 tio.json.gz --species TiO \
  --isotope-fraction 0.7372 --temperature-range 2200 3400 --pressure-range 10 10000
```

The output keeps the native wavelength samples in increasing vacuum Å and converts pressure from bar to dyn/cm². The temperature and pressure ranges select tabulated points inclusively. Each selected axis needs at least two points. The solver clamps values outside the supplied temperature and pressure range. Provide tables that cover the line-forming layers and include the contributing isotopologues for your mixture. [ExoMolOP data and formats](https://arxiv.org/abs/2009.00687)


## Update the spherical hydrostatic boundary

With `--photospheric-anchor`, GiantLayers integrates the native Rosseland opacity through the radial column mass and locates optical depth one between layers. The solver sets the gravity specified by `--logg` at that interpolated radius. The grid must bracket the photosphere.

The native CSV also contains the radiative acceleration from the frequency integral of flux times mass extinction divided by the speed of light. With `--radiative-force`, GiantLayers subtracts that acceleration from gravity when it updates gas pressure and the mixing-length pressure scale. The solver requires positive effective gravity on the mass grid. By default, the acceleration is reported and hydrostatic balance uses gravity alone.

The optional boundary updates use the latest solved profile between wavelength stages. After the final wavelength stage, GiantLayers iterates the boundary and hydrostatic structure at fixed temperature and reports the resulting luminosity error.

## Check depth resolution and compare with a reference

In the development version, combine `--deep-layers` with `--deep-grid-power 2` to concentrate the deep points toward the inner boundary while retaining the same number of points and the same depth endpoints. Use `--restart-csv` with `--interpolate-restart` to move a solved temperature profile to a different number of layers. The previous mass grid must cover the new one. GiantLayers recomputes the gas structure and transfer on the new grid. The solver reports how many final layers used molecular cross sections clamped to the edge of their temperature or pressure axes.

Compare the exported model with a published model at the same stellar parameters. The comparison interpolates only over the shared Rosseland depth range and reports the temperature difference in kelvin and the pressure differences in dex.

```sh
python -m giantlayers.compare_atmospheres published.mod.gz trial.mod.gz > structure-comparison.csv
```

## Carbon-rich mixtures

Change the carbon abundance in the input JSON to set C/O. The native LTE equation of state includes CO, CN, C₂, HCN, and acetylene. HCN and acetylene use NASA Glenn Gibbs energies for the dissociation equilibrium and ExoMol partition functions for explicit line strengths. The tabulated dissociation constants cover 1,000–6,000 K. Above the available partition-function temperatures (4,000 K for HCN and 5,000 K for acetylene), the line-strength partition function uses the closest tabulated value. [NASA Glenn thermochemistry](https://ntrs.nasa.gov/citations/20020085330) · [HCN partition function](https://exomol.com/data/data-types/opacity/HCN/1H-12C-14N/Harris/) · [Acetylene line list](https://academic.oup.com/mnras/article/493/2/1531/5716685)

To include carbon-bearing molecular absorption, convert the published CN and C₂ cross sections and add them with `--molecular-cross-section`. The abundance inputs independently set the molecular number densities used to convert the cross sections to mass opacity. [CN MoLLIST](https://www.exomol.com/data/data-types/linelist/CN/12C-14N/MoLLIST-CN/) · [C₂ 8states](https://exomol.com/data/data-types/opacity/C2/12C2/8states/)

Use `--adaptive-knots 4` to allow up to four additional temperature controls at layers with large luminosity residuals during each stage. `--peak-error-weight` changes the weight assigned to the maximum residual in the acceptance criterion. The final CSV records the luminosity error of every layer.

For mixtures that strongly change the line blanketing, start with reduced line opacity using `--opacity-ramp 0.1 0.5`. GiantLayers solves on the startup wavelength grid at each supplied fraction, then at full line opacity. The continuum and LTE chemistry are evaluated throughout. The optional finer wavelength stages use full opacity.

Use `--workers 4` to distribute finite-difference temperature perturbations among processes. For a fixed-profile `--steps 0` or `--verify-full` evaluation, the processes integrate separate wavelength bands and combine the flux, Rosseland moments, and radiative acceleration. Spectrum output uses a serial sweep.

Use `--full-grid-correction` to calibrate the last sampled stage against a transfer on the complete input grid. The solver holds the initial difference between the sampled and full-grid luminosity residuals fixed while it computes temperature corrections on the sampled grid. The command then recomputes the final luminosity on the full grid and reports the measured residual.

## Export the emergent spectrum

Pass `--output-spectrum spectrum.csv` to evaluate the final atmosphere and write the outgoing surface flux on the final wavelength grid. The spectrum reports flux per Hz and per Å at the model's reference radius, plus the integration weight of each sample. Summing flux per Hz times its weight recovers the integrated radiative surface flux at the reference radius.

## Add hydrogen and helium collision-induced absorption

Download the rectangular H₂–H₂, H₂–He, H₂–H, and He–H files from HITRAN and pass each selected pair with `--cia FILE.cia`. The reader also accepts gzip-compressed tables. It interpolates in temperature and wavenumber, then multiplies the coefficient by the two LTE species number densities and divides by mass density. CIA is continuum absorption and retains its full weight during an opacity ramp. Temperature outside a table's axis uses the nearest tabulated value. Wavenumbers outside the tabulated range contribute zero. [HITRAN CIA files](https://hitran.org/cia/) · [Format and reference data](https://hitran.org/data/CIA/CIA_Readme_2024.pdf)


## Use the solver from Python

```python
import json
from importlib.resources import files
from giantlayers import initialize
from giantlayers.native_transfer import read_wavelengths
from giantlayers.spherical_heating import balance_all_layers

examples = files("giantlayers.examples")
with examples.joinpath("solar-gs98.json").open(encoding="utf-8") as stream:
    abundances = json.load(stream)
seed, boundary = initialize(3500, 1, 1, abundances, layers=12)
wavelengths = read_wavelengths(examples.joinpath("continuum-wavelengths.dat"))
model, diagnostics, history, knots, convective_start = balance_all_layers(
    seed, wavelengths, (), {}, boundary=boundary, core_rays=12, steps=1
)
```

This example makes one correction with continuum extinction. Pass loaded atomic transitions through the third argument and loaded molecular tables through `tables`. When supplying your own initial seed and boundary, use the same abundance mixture in the structure and opacity calculations.
