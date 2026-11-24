# Quickstart

Clone the repository and run the prototype with Python 3. The solver uses the Python standard library.

```sh
git clone https://github.com/dnidever/giantlayers.git
cd giantlayers
python3 gray_atmosphere.py --teff 4300 --logg 1.5 > atmosphere.csv
```

This sets the effective temperature to 4300 K and the base-10 logarithm of surface gravity to 1.5 in cgs units. The default grid has 81 layers, extends to gray optical depth 100, and uses a constant opacity of 0.1 cm²/g.

## Model output

| CSV column | Meaning |
| --- | --- |
| `tau` | Gray optical depth, increasing inward from zero |
| `column_mass_g_cm2` | Column mass in g/cm² |
| `temperature_K` | Temperature in K |
| `pressure_dyn_cm2` | Gas pressure in dyn/cm² |
| `flux_erg_cm2_s` | Net outward radiative flux in erg/(cm² s) |
| `relative_flux_residual` | Fractional difference between the ray-sweep flux and the prescribed flux |

The program prints maximum relative flux and hydrostatic residuals to standard error, so they stay visible when the CSV goes to a file. For grid and opacity options, run:

```sh
python3 gray_atmosphere.py --help
```

## Run with a two-group opacity table

The nonlinear prototype can read both mass opacities from a CSV file.

```sh
python3 nonlinear_atmosphere.py --opacity-table opacities.csv > atmosphere.csv
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
python3 odf_transfer.py p00big2.bdf --marcs-model model.mod.gz \
  --marcs-continuum model.opa.gz > bands.csv
python3 compare_marcs_flux.py bands.csv flx_wavelengths.vac.gz model.flx.gz > comparison.csv
```

GiantLayers uses the MARCS true continuum absorption with ATLAS9 line opacity on the fixed temperature and pressure structure. The transfer uses plane-parallel two-stream geometry. Add `--continuum-scattering` to include isotropic coherent scattering from the MARCS table. The comparison integrates published flux samples only over bands fully covered by the MARCS wavelength grid.

### Inspect radiative heating

Add `--heating-output heating.csv` to write the heating profile while keeping the band flux CSV on standard output. The profile includes column mass, local net heating per gram, integrated heating in each depth cell, and net outward radiative flux. A positive heating rate means the matter absorbs more energy than it emits. The diagnostic uses the conservative two-stream solver and prints the closure error between the summed cell heating and the change in radiative flux across the atmosphere.

## Write a MARCS atmosphere file

The exporter writes a complete MARCS model as plain text or gzip. The command below reads an existing model and re-exports its header, 92 elemental abundances, structure, thermodynamics, and partial pressures.

```sh
python3 marcs_writer.py input.mod.gz output.mod
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
python3 lte_opacity.py model.mod.gz transitions.csv mixture.json 15570 15585 0.05 opacity.csv --vmic 2
```

The reference MARCS structure supplies the temperature and gas pressure. The solver recomputes the mass density, electron density and species populations at each layer. The output contains sampled line absorption, continuum absorption, and continuum scattering in cm²/g. For isotopologue lines, include the isotopic fraction in the oscillator strength used in the transition CSV.

The bundled data cover 1,000–10,000 K and include neutral atoms, first ions, H⁻, 20 neutral diatomic molecules, and H₂O. Water equilibrium uses the NASA Glenn data through 6,000 K. Continuum terms include H⁻ bound-free and free-free absorption, hydrogenic H I bound-free absorption, electron–ion free-free absorption, and Thomson scattering. The hydrogenic continuum terms use unit Gaunt factors. Additional polyatomic and charged molecules, higher ionization stages, metal photoionization, and calculated collisional line damping are future extensions.

[Atomic and diatomic thermochemistry data](https://doi.org/10.1051/0004-6361/201526961) · [H⁻ continuum data in the ATLAS9 HMINOP routine](https://wwwuser.oats.inaf.it/fiorella.castelli/sources/atlas9g/atlas9mem.for) · [APOGEE DR17 line data reference](https://ui.adsabs.harvard.edu/abs/2021AJ....161..254S/abstract)

### Evaluate radiative heating with native opacities

Provide a text file containing increasing vacuum wavelengths in Å. The transfer uses the same grid for all layers and integrates in frequency using trapezoidal weights. Include broad wavelength coverage and enough samples to resolve the lines when interpreting the result as a bolometric energy balance.

```sh
python3 native_transfer.py model.mod.gz transitions.csv mixture.json wavelengths.txt heating.csv
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
python3 benchmark_sampling.py model.mod.gz transitions.csv mixture.json tio.json.gz h2o.json.gz --strides 4 2 1 > sampling.csv
```

For expanded optical and infrared coverage, include ExoMol HyVO and Li2015 CO tables on the same grid. If the CO table replaces the explicit CO transitions, remove those transitions from the input CSV. A header-only transition CSV is valid when using cross sections alone. These files add molecular opacity over the common sampled range and inherit the temperature and pressure coverage of the source tables. [VO HyVO](https://exomol.com/data/molecules/VO/51V-16O/HyVO/) · [CO Li2015](https://www.exomol.com/data/molecules/CO/12C-16O/Li2015/)

### Run on the published MARCS wavelength grid

Pass the published vacuum grid directly to native transfer with `--quadrature wavelength`. The reader accepts the original gzip file, and the transfer processes the opacity in wavelength blocks while retaining quadrature weights across block boundaries.

```sh
python3 native_transfer.py model.mod.gz transitions.csv mixture.json flx_wavelengths.vac.gz full-grid.csv --quadrature wavelength --molecular-cross-section tio.json.gz
```

To compare coarse startup sampling with the full published grid, pass the same input file to the sampling benchmark. Stride 10 retains every tenth sample and both endpoints. Stride 1 uses all 100,724 wavelengths.

```sh
python3 benchmark_sampling.py model.mod.gz transitions.csv mixture.json tio.json.gz --reference-wavelengths flx_wavelengths.vac.gz --strides 10 1 --quadrature wavelength > marcs-sampling.csv
```

Use the native-resolution opacity tables when comparing the grids. The transfer evaluates the available opacities at the MARCS sample positions. This matches the published wavelength axis and the wavelength integration convention. The geometry, opacity inventory, and convection model also affect how closely the calculated atmosphere reproduces MARCS. [MARCS wavelength download](https://marcs.astro.uu.se/docs.html) · [MARCS methods](https://arxiv.org/abs/0805.0554)

### Compare spherical and plane-parallel surface transfer

Convert the published MARCS oscillator-strength revisions from VALD air wavelengths into the native vacuum-wavelength transition format. These files provide a subset of the atomic lines used by MARCS. The converted transitions use Doppler broadening. [MARCS auxiliary atomic data](https://marcs.astro.uu.se/docs.html#Auxiliary)

```sh
python3 marcs_atomic_lines.py astroph_gfVALD1.dat astroph_gfVALD2.dat > revised-atomic-lines.csv
```

Use a spherical model and the same opacity inputs to compare surface transfer in both geometries. Pass several `--molecular-cross-section` options to include multiple absorbers. The spherical solver traces rays through atmospheric shells and integrates the emergent intensity over the apparent disk. The plane-parallel case uses the same surface angles and quadrature weights. Both modes use LTE true absorption and fixed temperatures and gas pressures for this benchmark.

```sh
python3 spherical_transfer.py model.mod.gz revised-atomic-lines.csv mixture.json flx_wavelengths.vac.gz --stride 50 --core-rays 24 --molecular-cross-section tio.json.gz > spherical-comparison.csv
```

The CSV contains monochromatic surface fluxes per Hz and frequency weights. The printed spherical ratio normalizes the flux to the MARCS reference radius. The planar ratio uses the corresponding target flux at the outer boundary. Reduce `--stride` and increase `--core-rays` to test sensitivity to the wavelength and ray grids.

### Compute spherical heating with scattering

Use the same spherical model and opacity inputs with the depth-resolved solver. It iterates the isotropic coherent continuum scattering source on the shared impact-parameter rays. Net heating is calculated from true absorption and the difference between the mean intensity and the LTE Planck function. Use `--absorption-only` to isolate the effect of scattering.

```sh
python3 spherical_heating.py model.mod.gz revised-atomic-lines.csv mixture.json flx_wavelengths.vac.gz --stride 100 --core-rays 24 --molecular-cross-section tio.json.gz > spherical-heating.csv
```

The CSV reports integrated radiative flux, luminosity, and pointwise true-absorption heating at each layer. The shell columns give the mean net heating and energy gain between the current layer and the next deeper layer, calculated from the luminosity difference and shell mass. Positive heating means the gas gains energy. The shell values balance the luminosity change by construction. Standard error reports how far a separate trapezoidal integration of the pointwise rates departs from that budget. Use that difference to assess whether pointwise heating is suitable for a temperature update. The MARCS temperature, pressure, and convection profiles are held fixed.

### Try spherical temperature corrections

Add `--adjust-layers 3 --steps 2` to the spherical heating command to try two damped temperature corrections in the outer three layers. Each trial recomputes LTE chemistry, opacities, and spherical scattering. The correction uses finite differences to reduce the change in radiative plus convective luminosity across the selected shells. The deeper temperature profile, radii, gas pressure and published MARCS convective contribution stay fixed.

The solver reports the maximum outer-shell luminosity residual and the surface total-flux error separately. The correction targets local shell balance. Fitting the absolute emergent flux requires a separate constraint. The CSV includes the initial and corrected temperature.
