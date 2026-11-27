"""Henyey mixing-length convection on a specified spherical atmosphere."""

import argparse
from bisect import bisect_right
import csv
import json
import math
from pathlib import Path
import sys

from .benchmark_sampling import nested_wavelengths
from .gray_atmosphere import SIGMA
from .hminus import K_B, K_EV
from .lte_chemistry import Thermochemistry
from .lte_opacity import H_EV_S, LIGHT, composition, layer_opacities, prepare_layers, read_lines
from .marcs_writer import read_document
from .molecular_cross_section import MolecularCrossSection
from .native_transfer import read_wavelengths, spectral_quadrature
from .odf_transfer import planck_nu
from .spherical_hydrostatic import G, SOLAR_MASS, reconstruct_for_transfer


def lte_radiation_pressure(temperature):
    """Isotropic LTE radiation pressure in dyn/cm²."""
    return 4 * SIGMA * temperature**4 / (3 * LIGHT)


def slopes(axis, values, quadratic_endpoints=False):
    """Estimate knot derivatives with one-sided secants at the endpoints."""
    slopes_at_edges = [(y-x)/(b-a) for a, b, x, y in
                       zip(axis, axis[1:], values, values[1:])]
    result = [slopes_at_edges[0]]
    for i in range(1, len(axis)-1):
        left, right = axis[i]-axis[i-1], axis[i+1]-axis[i]
        result.append((right*slopes_at_edges[i-1]+left*slopes_at_edges[i])/(left+right))
    result.append(slopes_at_edges[-1])
    if quadratic_endpoints and len(axis) > 2:
        result[0] -= ((axis[1]-axis[0])/(axis[2]-axis[0])
                      * (slopes_at_edges[1]-slopes_at_edges[0]))
        result[-1] += ((axis[-1]-axis[-2])/(axis[-1]-axis[-3])
                       * (slopes_at_edges[-1]-slopes_at_edges[-2]))
    return result


def interpolate(axis, values, x):
    i = min(max(0, bisect_right(axis, x)-1), len(axis)-2)
    fraction = (x-axis[i])/(axis[i+1]-axis[i])
    return (1-fraction)*values[i] + fraction*values[i+1]


class Enthalpy:
    """Estimate LTE enthalpy from populations and tabulated temperature derivatives."""

    def __init__(self, chemistry):
        self.chemistry = chemistry
        self.log_t = [math.log(t) for t in chemistry.temperatures]
        self.atomic_slopes = {}
        for atom, data in chemistry.atoms.items():
            for species, key in ((atom, "neutral_q"), (atom+"+", "ion_q")):
                self.atomic_slopes[species] = slopes(self.log_t, [math.log(q) for q in data[key]])
        self.molecular_slopes = {}
        for species, data in chemistry.molecules.items():
            count = sum(t <= data.get("maximum_formation_temperature", math.inf)
                        for t in chemistry.temperatures)
            axis = self.log_t[:count]
            values = [math.log(10)*p for p in data["log10_kp_cgs"][:count]]
            self.molecular_slopes[species] = (axis, slopes(axis, values))

    def __call__(self, gas):
        kt = K_B * gas.temperature
        log_t = math.log(gas.temperature)
        energy = {"e-": 2.5*kt, "H-": 2.5*kt - 0.754209*K_B/K_EV}
        for atom, data in self.chemistry.atoms.items():
            if atom not in gas.number_densities:
                continue
            ion = atom + "+"
            energy[atom] = kt * (2.5 + interpolate(self.log_t, self.atomic_slopes[atom], log_t))
            energy[ion] = kt * (2.5 + interpolate(self.log_t, self.atomic_slopes[ion], log_t))
            energy[ion] += data["ionization_ev"] * K_B / K_EV
        for name, data in self.chemistry.molecules.items():
            if not gas.number_densities.get(name, 0):
                continue
            axis, derivative = self.molecular_slopes[name]
            # van 't Hoff: the dissociation enthalpy is kT d(log Kp)/d(log T).
            energy[name] = (sum(energy[atom] for atom in data["atoms"])
                            - kt*interpolate(axis, derivative, log_t))
        return math.fsum(n*energy[name] for name, n in gas.number_densities.items() if n) / gas.mass_density

    def derivatives(self, temperature, pressure, abundances, *, include_radiation=False):
        """Return Cp, the density derivative Q, and the adiabatic gradient."""
        lower = max(self.chemistry.temperatures[0], temperature*0.998)
        upper = min(self.chemistry.temperatures[-1], temperature*1.002)
        gas = self.chemistry.solve(temperature, pressure, abundances)
        total_pressure = pressure
        if include_radiation:
            total_pressure += lte_radiation_pressure(temperature)
            # Keep the pressure perturbation safely inside the positive gas-pressure domain.
            upper = min(upper, temperature * (1 + 0.05 * pressure /
                                               lte_radiation_pressure(temperature)))
            a, b = (self.chemistry.solve(t, total_pressure-lte_radiation_pressure(t),
                                         abundances, initial=gas) for t in (lower, upper))
            ha, hb = (self(state) + 4*lte_radiation_pressure(state.temperature)/state.mass_density
                      for state in (a, b))
        else:
            a, b = (self.chemistry.solve(t, pressure, abundances, initial=gas) for t in (lower, upper))
            ha, hb = self(a), self(b)
        cp = (hb-ha)/(upper-lower)
        q = math.log(a.mass_density/b.mass_density)/math.log(upper/lower)
        if cp <= 0 or q <= 0:
            raise ValueError(f"Invalid thermodynamic derivative at T={temperature:g} K and P={pressure:g}")
        adiabatic = total_pressure*q/(gas.mass_density*temperature*cp)
        return gas.mass_density, cp, q, adiabatic


def henyey_flux(temperature, pressure, density, gravity, rosseland, gradient, cp, q,
                adiabatic, alpha=1.5, y=0.076, nu=8.0):
    """Evaluate the Henyey convective flux and velocity in cgs units."""
    excess = gradient - adiabatic
    if excess <= 0 or alpha == 0 or rosseland == 0:
        return 0.0, 0.0
    hp = pressure/(gravity*density)
    tau = density*rosseland*alpha*hp
    inverse_efficiency = 8*SIGMA*temperature**3*tau/(density*cp*(1+y*tau*tau))
    driving = alpha**2*gravity*hp*q*excess/nu
    # Solve v² + v/B = A² (nabla - nabla_ad), where Gamma = Bv.
    velocity = 2*driving/(math.hypot(inverse_efficiency, 2*math.sqrt(driving))+inverse_efficiency)
    contrast = velocity**2*nu/(alpha**2*gravity*hp*q)
    return 0.5*density*cp*temperature*velocity*alpha*contrast, velocity


def native_rosseland(model, wavelengths, lines, mixture, tables, chemistry, vmic):
    """Compute native Rosseland mean total extinction on the sampled grid."""
    frequencies, weights = spectral_quadrature(wavelengths, "wavelength")
    grids = [(table, table.on_grid(wavelengths)) for table in tables]
    for layer in prepare_layers(model, lines, mixture, vmic, chemistry):
        temperature = layer[0].temperature
        numerator = denominator = 0.0
        for frequency, weight, (line, absorb, scatter) in zip(
            frequencies, weights, layer_opacities(layer, wavelengths, grids)
        ):
            x = H_EV_S*frequency/(K_EV*temperature)
            weighted = weight*planck_nu(frequency, temperature)*x/(-math.expm1(-x))
            numerator += weighted
            denominator += weighted/(line+absorb+scatter)
        yield numerator/denominator


def evaluate(model, mixture, rosseland=None, chemistry=None, radiative_acceleration=(), *,
             include_radiation=False):
    """Compute the mixing-length structure on the supplied pressure and radius grid."""
    chemistry = chemistry or Thermochemistry()
    enthalpy = Enthalpy(chemistry)
    reference = {e: model.abundances[row["atomic_number"]-1] for e, row in chemistry.atoms.items()}
    abundances = composition(reference, mixture)
    alpha, nu, y = map(float, model.header[9].split()[:3])
    if not all(math.isfinite(x) for x in (alpha, nu, y)) or alpha < 0 or y < 0 or nu <= 0:
        raise ValueError("Invalid mixing-length parameters in the MARCS header")
    radius = float(model.header[7].split()[0])
    stellar_mass = float(model.header[5].split()[0])*SOLAR_MASS
    if stellar_mass <= 0:
        raise ValueError("Use a spherical MARCS reference model")
    p = [row[5] for row in model.structure]
    t = [row[3] for row in model.structure]
    if any(a >= b for a, b in zip(p, p[1:])):
        raise ValueError("Gas pressure must increase inward")
    log_t = list(map(math.log, t))
    gradients = slopes(list(map(math.log, p)), log_t, quadratic_endpoints=True)
    if include_radiation:
        pressures = [gas + lte_radiation_pressure(temp) for gas, temp in zip(p, t)]
        first_thick = next((i for i, row in enumerate(model.structure) if row[0] >= 0), len(p))
        pivot = max(0, first_thick - 1)
        if first_thick < len(p) and len(p) - pivot >= 2:
            thick_pressures = pressures[pivot:]
            if any(a >= b for a, b in zip(thick_pressures, thick_pressures[1:])):
                raise ValueError("Gas plus LTE radiation pressure must increase in optically thick layers")
            thick_gradients = slopes(list(map(math.log, thick_pressures)), log_t[pivot:],
                                     quadratic_endpoints=True)
            gradients[first_thick:] = thick_gradients[first_thick-pivot:]
    ref_flux = float(model.header[2].split()[0])
    if rosseland is None:
        rosseland = [row[1] for row in model.thermodynamics]
    force = radiative_acceleration if len(radiative_acceleration) else (0.0,)*len(p)
    if len(force) != len(p):
        raise ValueError("Radiative acceleration must give a value at every layer")
    for i, (row, pressure, temperature, gradient, opacity, thermo) in enumerate(
        zip(model.structure, p, t, gradients, rosseland, model.thermodynamics), 1
    ):
        r = radius - row[2]
        thick_radiation = include_radiation and row[0] >= 0
        density, cp, q, adiabatic = enthalpy.derivatives(temperature, pressure, abundances,
                                                        include_radiation=thick_radiation)
        # Radiation pressure contributes directly to the total-pressure gradient.
        gravity = G*stellar_mass/r**2 - (0 if thick_radiation else force[i-1])
        mixing_pressure = pressures[i-1] if thick_radiation else pressure
        if gravity <= 0:
            raise ValueError("Mixing-length convection requires positive effective gravity")
        flux, velocity = henyey_flux(temperature, mixing_pressure, density, gravity,
                                     opacity, gradient, cp, q, adiabatic, alpha, y, nu)
        target = ref_flux*(radius/r)**2
        yield (i, row[0], temperature, pressure, opacity, gradient, adiabatic, cp, q, velocity,
               thermo[4], flux, flux/target, thermo[5])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model")
    parser.add_argument("--composition", type=Path)
    parser.add_argument("--hydrostatic", action="store_true")
    parser.add_argument("--wavelengths", type=Path, help="Compute native Rosseland opacities on this grid")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--lines", type=Path)
    parser.add_argument("--vmic", type=float, default=2)
    parser.add_argument("--molecular-cross-section", action="append", default=[])
    parser.add_argument("--lte-radiation-pressure", action="store_true",
                        help="Include LTE radiation in convection at Rosseland depth >= 1")
    args = parser.parse_args()
    try:
        chemistry = Thermochemistry()
        model = read_document(args.model)
        mixture = json.loads(args.composition.read_text()) if args.composition else {}
        if args.hydrostatic:
            model, _ = reconstruct_for_transfer(model, chemistry, mixture)
        rosseland = None
        if args.wavelengths:
            wave = nested_wavelengths(read_wavelengths(args.wavelengths), args.stride)
            lines = read_lines(args.lines) if args.lines else []
            tables = [MolecularCrossSection(path) for path in args.molecular_cross_section]
            rosseland = list(native_rosseland(model, wave, lines, mixture, tables, chemistry, args.vmic))
        rows = list(evaluate(model, mixture, rosseland, chemistry,
                             include_radiation=args.lte_radiation_pressure))
    except (IndexError, OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    writer = csv.writer(sys.stdout)
    writer.writerow(("layer", "log_tau_ross", "temperature_K", "gas_pressure_dyn_cm2",
                     "rosseland_opacity_cm2_g", "temperature_gradient", "adiabatic_gradient",
                     "cp_erg_g_K", "thermal_expansion_q", "convective_velocity_cm_s",
                     "marcs_convective_velocity_cm_s", "convective_flux_erg_cm2_s",
                     "convective_flux_fraction", "marcs_convective_flux_fraction"))
    writer.writerows(rows)


if __name__ == "__main__":
    main()
