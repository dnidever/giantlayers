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

The bundled data cover 1,000–10,000 K and include neutral atoms, first ions, H⁻, and 20 neutral diatomic molecules. Continuum terms include H⁻ bound-free and free-free absorption, hydrogenic H I bound-free absorption, electron–ion free-free absorption, and Thomson scattering. The hydrogenic continuum terms use unit Gaunt factors. Polyatomic and charged molecules other than H⁻, higher ionization stages, metal photoionization, and calculated collisional line damping are future extensions.

[Atomic and diatomic thermochemistry data](https://doi.org/10.1051/0004-6361/201526961) · [H⁻ continuum data in the ATLAS9 HMINOP routine](https://wwwuser.oats.inaf.it/fiorella.castelli/sources/atlas9g/atlas9mem.for) · [APOGEE DR17 line data reference](https://ui.adsabs.harvard.edu/abs/2021AJ....161..254S/abstract)
