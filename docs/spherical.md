# Spherical LTE atmosphere model

GiantLayers begins with an analytic gray temperature profile on a radial column-mass grid. Its inputs are effective temperature, surface gravity, stellar mass, a chemical abundance mixture, and opacity data. A reference radius follows from R = √(GM/g). The optional photospheric update locates this radius at Rosseland optical depth one.

## Hydrostatic structure and chemistry

The radial column mass increases inward. GiantLayers integrates the gas-pressure gradient dP/dm = GM/r² − g_rad and the depth gradient dr/dm = −1/ρ. The outer gas pressure is a boundary input. Set `--radiative-force` to update g_rad from the calculated radiation field. The pressure solution requires positive effective gravity.

LTE chemistry enforces each element's abundance, charge neutrality, and the ideal-gas particle pressure. The network uses neutral and singly ionized atoms, H⁻, and supported molecules. Partition functions and equilibrium constants determine the species populations. Those populations set the electron pressure, mass density, molecular absorption, and relevant thermodynamic derivatives. [Atomic and diatomic data](https://doi.org/10.1051/0004-6361/201526961) · [NASA Glenn polyatomic equilibrium](https://ntrs.nasa.gov/citations/20020085330)

## Transfer, convection, and luminosity balance

The transfer sweeps incoming and outgoing rays on their spherical impact parameters. The outer boundary has zero incident radiation. A smooth transition from a thermal inner boundary to diffusion uses the optical thickness of each ray. The solver uses LTE emissivity, isotropic coherent scattering, and a linear source function in each ray segment.

Atomic transitions, molecular cross sections, and continuum extinction are evaluated at each supplied wavelength. The continuum includes hydrogen absorption, electron–ion free-free absorption, electron and Rayleigh scattering, and optional collision-induced absorption. The same transfer calculates the Rosseland mean for Henyey mixing-length convection. [MARCS methods and convection](https://arxiv.org/abs/0805.0554)

The temperature solver targets L = 4πr²(F_rad + F_conv) = 4πR²σT_eff⁴ at each layer. It perturbs the temperatures at selected depths, recomputes the hydrostatic and radiation solutions, and updates the temperatures from the response of the luminosity residuals. An optional normalized heating constraint supplies information about optically thin layers where luminosity varies little with local temperature. Temperature corrections are bounded and tested before acceptance.

## Interpreting a result

Use the CSV maximum luminosity residual to check convergence at every layer. Surface flux, agreement with a published temperature profile, and a low RMS residual each provide different information. Compare the same model on finer wavelength and depth grids and inspect the molecular table coverage report.

The scientific scope is static, dust-free atmospheres in LTE. The result depends on the supplied opacity inventory and on the pressure and mass endpoints. The MARCS-format export includes an isotropic LTE estimate of radiation pressure and uses zero turbulent pressure. See [Validation](validation.md) for the benchmarks and the [Quickstart](quickstart.md) for export commands.
