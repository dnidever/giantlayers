# Gray atmosphere model

The prototype assumes a static, plane-parallel atmosphere in local thermodynamic equilibrium. Gravity and gray absorption opacity are constant. Radiation carries the entire flux.

## Structure and transfer

With column mass *m* increasing inward from the surface, gray optical depth is τ = κ*m* and hydrostatic equilibrium gives gas pressure *P* = *gm*. Here κ is the absorption opacity and *g* is surface gravity.

The two-stream quadrature uses outward and inward rays with direction cosines ±1/√3. Radiative equilibrium fixes the bolometric flux to σ times the fourth power of the effective temperature. The solver derives the temperature profile from the gray moment equations with zero incident intensity at the surface.

## Numerical checks

The formal ray sweep integrates the transfer equation across each depth interval using a linear source function. It uses the surface vacuum boundary for the incoming ray and the diffusion boundary for the outgoing ray at the bottom. The flux residual compares the resulting flux with the prescribed bolometric flux at every depth. The hydrostatic residual checks the pressure gradient with respect to column mass against the input gravity.

This model provides a benchmark for future frequency-dependent opacity, convection, and spherical geometry.
