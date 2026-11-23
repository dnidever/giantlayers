# GiantLayers

GiantLayers is an early-stage stellar atmosphere solver for cool giants.

The current prototype computes a gray, plane-parallel LTE atmosphere in hydrostatic and radiative equilibrium. It uses constant opacity and an Eddington two-stream closure. A formal ray sweep checks the radiative flux of the moment solution.

## Run the prototype

Use Python 3. No third-party packages are required.

```sh
python3 gray_atmosphere.py --teff 4300 --logg 1.5 > atmosphere.csv
```

The CSV contains optical depth, column mass, temperature, gas pressure, radiative flux, and relative flux residual. The program prints the maximum relative flux and hydrostatic residuals to standard error. Run with `--help` for the grid and opacity options.

The prototype is a numerical benchmark for future frequency-dependent opacity, convection, and spherical geometry.
