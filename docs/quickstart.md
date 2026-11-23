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
