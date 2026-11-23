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
