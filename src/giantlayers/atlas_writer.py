"""Export native atmospheres as Kurucz ATLAS READ DECK6 records."""

import gzip
import math
from pathlib import Path

from .hminus import K_B
from .lte_opacity import LIGHT
from .native_marcs import exterior_column


DECK_COLUMNS = "RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB"
# Conventional ATLAS9 deck flags. The originating opacities are described in TITLE.
OPACITY_FLAGS = (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0)


def abundance_records(abundances, *, atlas12=False):
    """Write ATLAS number fractions for H and He and log fractions for metals."""
    nuclei = [10**(value-12) if value > -99 else 0.0 for value in abundances]
    total = math.fsum(nuclei)
    records = [
        "ABUNDANCE SCALE   1.00000 ABUNDANCE CHANGE"
        + "".join(f"{i:2d}{fraction/total:8.5f}" for i, fraction in enumerate(nuclei[:2], 1))
    ]
    metals = [math.log10(value/total) if value else -99.0 for value in nuclei[2:]]
    metals.extend([-99.0]*(99-len(nuclei)))
    for start in range(0, len(metals), 6):
        records.append(" ABUNDANCE CHANGE" + "".join(
            f"{i:3d}{value:7.2f}" for i, value in enumerate(metals[start:start+6], start+3)
        ))
    if atlas12:
        symbols = (
            "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
            "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te "
            "I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os "
            "Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es"
        ).split()
        records.extend((" ABUNDANCE TABLE",
                        f"    1H {nuclei[0]/total:10.6f}       2He{nuclei[1]/total:10.6f}"))
        for start in range(0, len(metals), 5):
            records.append("".join(
                f"{i:5d}{symbols[i-1]:2s}{max(value, -99):7.3f}{0.0:6.3f}"
                for i, value in enumerate(metals[start:start+5], start+3)
            ))
    return records


def sound_speed(temperature, pressure, density, convection, boundary):
    """Evaluate the equilibrium adiabatic sound speed from the native EOS."""
    adiabatic, q = convection[6], convection[8]
    step = 0.002
    low, high = (boundary.chemistry.solve(temperature, pressure*math.exp(sign*step),
                                         boundary.abundances).mass_density for sign in (-1, 1))
    derivative = math.log(high/low)/(2*step) - q*adiabatic
    if derivative <= 0 or not math.isfinite(derivative):
        raise ValueError("Invalid adiabatic compressibility for ATLAS sound speed")
    return math.sqrt(pressure/(density*derivative))


def write_atlas(model, result, boundary, vmic, path, *, status, atlas12=False):
    """Write a spherical GiantLayers solution as a depth-only ATLAS deck."""
    path = Path(path)
    layers = len(model.structure)
    if atlas12 and layers > 72:
        raise ValueError("The reference ATLAS12 reader supports at most 72 depth records")
    if not 1 <= layers <= 999:
        raise ValueError("ATLAS READ DECK6 permits 1–999 depth records")
    if any(len(column) != layers for column in (
        boundary.column_mass, model.thermodynamics, result[0], result[7], result[9]
    )):
        raise ValueError("ATLAS export requires complete transfer diagnostics")
    teff = float(model.header[1].split()[0])
    logg = math.log10(boundary.gravity)
    alpha = float(model.header[9].split()[0])
    teff_field, gravity_field = f"{teff:.0f}.", f"{logg:8.5f}"
    alpha_field = f"{alpha:6.2f}"
    if len(teff_field) > 7 or len(gravity_field) > 8 or len(alpha_field) > 6:
        raise ValueError("Stellar parameters exceed ATLAS header field widths")
    title = f"GiantLayers spherical LTE export, luminosity {status}, native opacity"
    # The ATLAS12 card parser needs a delimiter before a full-width negative value.
    gravity_separator = " " if atlas12 and gravity_field.startswith("-") else ""
    records = [
        f"TEFF {teff_field:>7}  GRAVITY{gravity_separator}{gravity_field} LTE",
        f"TITLE {title[:74]:74s}",
        " OPACITY IFOP" + "".join(f"{flag:2d}" for flag in OPACITY_FLAGS),
        f" CONVECTION {'ON' if alpha else 'OFF':<3}{alpha_field} TURBULENCE OFF"
        + f"{0.:6.2f}"*4,
        *abundance_records(model.abundances, atlas12=atlas12),
        f"READ DECK6{layers:3d} {DECK_COLUMNS}" + ("" if atlas12 else ", FLXCNV,VCONV,VELSND"),
    ]
    exterior = exterior_column(model, boundary)
    previous = 0.0
    for mass, row, thermo, conv, acceleration in zip(
        boundary.column_mass, model.structure, model.thermodynamics, result[7], result[9]
    ):
        temperature, electron_pressure, pressure = row[3:6]
        rhox = mass + exterior
        values = (pressure, electron_pressure/(K_B*temperature), conv[4], acceleration,
                  vmic*1e5, conv[11], conv[9])
        if not atlas12:
            values += (sound_speed(temperature, pressure, thermo[2], conv, boundary),)
        if rhox <= previous or not all(map(math.isfinite, (rhox, temperature, *values))):
            raise ValueError("ATLAS export requires finite fields and increasing positive column mass")
        previous = rhox
        fields = [f"{rhox:15.8E}", f"{temperature:9.1f}", *(f"{value:10.3E}" for value in values)]
        if any(len(field) != width for field, width in zip(fields, (15, 9, *(10,)*len(values)))):
            raise ValueError("Atmosphere values exceed ATLAS depth field widths")
        records.append("".join(fields))
    if not math.isfinite(result[0][0]) or result[0][0] < 0:
        raise ValueError("ATLAS export requires a finite, nonnegative emergent radiative flux")
    # Outgoing isotropic hemisphere, P_rad = (2/3) F_rad / c.
    records.append(f"PRADK{2*result[0][0]/(3*LIGHT):11.4E}")
    if atlas12:
        # READ PUNCH closes the model stream when it reaches this card.
        records.append("BEGIN")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="ascii", newline="\n") as stream:
        stream.write("\n".join(records) + "\n")
