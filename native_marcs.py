"""Populate MARCS-format fields from a computed GiantLayers atmosphere."""

from bisect import bisect_left
from dataclasses import replace
import math

from gray_atmosphere import SIGMA
from hminus import K_B
from lte_chemistry import AMU
from lte_opacity import LIGHT, layer_opacities, prepare_layers
from spherical_hydrostatic import G, SOLAR_MASS, increments, integrate


_PRESSURE_SPECIES = (
    ("H", "H-", "H2", "H2+", "H2O", "OH", "CH", "CO", "CN", "C2"),
    ("N2", "O2", "NO", "NH", "TiO", "C2H2", "HCN", "C2H", "HS", "SiH", "C3H"),
    ("C3", "CS", "SiC", "SiC2", "NS", "SiN", "SiO", "SO", "S2", "SiS"),
)
_LISTED_SPECIES = frozenset(name for group in _PRESSURE_SPECIES for name in group)
# Reference values used only for the approximate composition labels.
_GS98_ALPHA = {"O": 8.83, "Mg": 7.58, "Si": 7.55, "S": 7.33, "Ca": 6.36, "Ti": 5.02}
_SOLAR_LUMINOSITY = 3.828e33


def _log(value):
    return math.log10(value) if value > 1e-99 else -99.0


def optical_depth(column, opacity, exterior_column):
    return integrate(increments(column, opacity), exterior_column*opacity[0])


def exterior_column(model, boundary):
    outer_radius = boundary.reference_radius - model.structure[0][2]
    return boundary.surface_pressure*outer_radius**2/(G*boundary.stellar_mass)


def photospheric_anchor(model, boundary, result):
    """Interpolate the radial layer index at Rosseland optical depth one."""
    rosseland = [row[4] for row in result[7]]
    depth = optical_depth(boundary.column_mass, rosseland, exterior_column(model, boundary))
    if not depth[0] < 1 < depth[-1]:
        raise ValueError("The column-mass grid must bracket Rosseland optical depth one")
    hi = bisect_left(depth, 1.0)
    lo = hi - 1
    return lo + (1 - depth[lo])/(depth[hi] - depth[lo])


def _partial_pressures(gas, pressure, chemistry):
    factor = K_B*gas.temperature
    densities = gas.number_densities
    rows = [tuple(_log(factor*densities.get(name, 0.0)) for name in group)
            for group in _PRESSURE_SPECIES]
    other = math.fsum(densities.get(name, 0.0) for name in chemistry.molecules
                      if name not in _LISTED_SPECIES)
    return (_log(pressure), *rows[0]), rows[1], (*rows[2], _log(factor*other))


def build_document(model, result, lines, tables, boundary, vmic):
    """Fill the export tables from the accepted LTE structure and opacity grid."""
    chemistry = boundary.chemistry
    wavelength = (5000.0,)
    grids = [(table, table.on_grid(wavelength)) for table in tables]
    layers = list(prepare_layers(model, lines, {}, vmic, chemistry))
    k5000 = [sum(next(layer_opacities(layer, wavelength, grids))) for layer in layers]
    convection = result[7]
    rosseland = [row[4] for row in convection]
    radius = boundary.reference_radius
    exterior = exterior_column(model, boundary)
    ross = list(map(_log, optical_depth(boundary.column_mass, rosseland, exterior)))
    tau5000 = list(map(_log, optical_depth(boundary.column_mass, k5000, exterior)))
    structure, thermo, partial = [], [], []
    for row, layer, lr, l5, conv, column in zip(
        model.structure, layers, ross, tau5000, convection, boundary.column_mass
    ):
        gas, pressure = layer[:2]
        temperature = gas.temperature
        radiation_pressure = 4*SIGMA*temperature**4/(3*LIGHT)
        structure.append((lr, l5, row[2], temperature, row[4], pressure, radiation_pressure, 0.0))
        total_particles = math.fsum(gas.number_densities.values())
        molecular_weight = gas.mass_density/(total_particles*AMU)
        thermo.append((lr, conv[4], gas.mass_density, molecular_weight, conv[9], conv[12],
                       column+exterior))
        partial.append(_partial_pressures(gas, pressure, chemistry))

    abundances = boundary.abundances
    masses = {element: 10**(value-12)*chemistry.atoms[element]["mass_amu"]
              for element, value in abundances.items()}
    total = math.fsum(masses.values())
    x, y = masses["H"]/total, masses.get("He", 0.0)/total
    z = math.fsum(v for e, v in masses.items() if e not in ("H", "He"))/total
    iron = abundances.get("Fe", -99.0)
    fe = iron - 7.50 if iron > -90 else -99.0
    alpha = [abundances[e]-a-fe for e, a in _GS98_ALPHA.items() if abundances.get(e, -99) > -90]
    alpha_label = math.fsum(alpha)/len(alpha) if alpha and fe > -90 else -99.0
    flux = float(model.header[2].split()[0])
    luminosity = 4*math.pi*radius**2*flux/_SOLAR_LUMINOSITY
    header = list(model.header)
    header[0] = "GiantLayers spherical LTE model"
    teff = float(header[1].split()[0])
    vmic, mass = float(header[4].split()[0]), boundary.stellar_mass/SOLAR_MASS
    alpha_mlt, nu, mixing_y, beta = map(float, header[9].split()[:4])
    header[1] = f"{teff:7.0f} Teff [K]"
    header[2] = f"{flux:12.4E} Flux [erg/cm2/s]"
    header[3] = f"{boundary.gravity:12.4E} Surface gravity [cm/s2]"
    header[4] = f"{vmic:5.1f} Microturbulence [km/s]"
    header[5] = f"{mass:5.1f} Mass [Msun]"
    header[6] = f"{fe:+6.2f}{alpha_label:+6.2f} Approximate [Fe/H] and mean [alpha/Fe], GS98"
    header[7] = f"{radius:12.4E} Reference radius [cm] at fixed mass anchor"
    header[8] = f"{luminosity:12.4E} Luminosity [Lsun]"
    header[9] = f"{alpha_mlt:6.2f}{nu:5.2f}{mixing_y:6.3f}{beta:5.2f} Convection parameters"
    header[10] = f"{x:.6f} {y:.6f} {z:.6e} are X, Y and Z"
    return replace(model, header=tuple(header), structure=tuple(structure),
                   thermodynamics=tuple(thermo), partial_pressures=tuple(zip(*partial)))
