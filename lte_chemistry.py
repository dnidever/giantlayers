"""Native ideal-gas LTE chemical equilibrium at fixed temperature and gas pressure."""

from bisect import bisect_right
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
import json
import math
from pathlib import Path

from hminus import K_B, K_EV, SAHA

AMU = 1.66053906660e-24
DEFAULT_TABLE = Path(__file__).with_name("data") / "barklem_collet.json"


@dataclass(frozen=True)
class GasState:
    temperature: float
    electron_density: float
    mass_density: float
    number_densities: dict
    partition_functions: dict
    masses_amu: dict
    residual: float


class Thermochemistry:
    """Tabulated ionization and neutral molecular equilibrium data."""

    def __init__(self, path=None):
        if path is None:
            path = DEFAULT_TABLE if DEFAULT_TABLE.is_file() else files("giantlayers_data").joinpath(DEFAULT_TABLE.name)
        elif not isinstance(path, Path):
            path = Path(path)
        with path.open(encoding="utf-8") as stream:
            values = json.load(stream)
        self.temperatures = values["temperatures"]
        self.atoms = values["atoms"]
        self.molecules = values["molecules"]

    def solve(self, temperature, gas_pressure, abundances, tolerance=1e-10, initial=None):
        """Solve element conservation, charge neutrality and the ideal-gas law.

        Abundances are log epsilon values keyed by element symbol, with H=12.
        Number densities are in cm⁻³, and gas pressure is in dyn/cm².
        """
        try:
            return self._solve_equilibrium(temperature, gas_pressure, abundances, tolerance, initial)
        except RuntimeError:
            if initial is not None:
                return self.solve(temperature, gas_pressure, abundances, tolerance)
            start = min(2400.0, self.temperatures[-1])
            if temperature >= start:
                raise
            # A warmer solution supplies a starting point for molecule-rich
            # states in which an atomic gas is a poor initial approximation.
            state = self._solve_equilibrium(start, gas_pressure, abundances, tolerance, None)
            steps = math.ceil((start-temperature)/75.0)
            for i in range(1, steps+1):
                t = start + (temperature-start)*(i/steps)
                state = self._solve_equilibrium(t, gas_pressure, abundances, tolerance, state)
            return state

    @lru_cache(maxsize=512)
    def _equilibrium_terms(self, temperature, abundance_items):
        """Reuse exact mass-action terms across pressure solves at the same temperature."""
        abundances = dict(abundance_items)
        molecules = {
            name: data for name, data in self.molecules.items()
            if all(e in abundances for e in data["atoms"])
        }
        coupled = {"H", *(e for data in molecules.values() for e in data["atoms"]
                          if temperature <= data.get("maximum_formation_temperature", math.inf))}
        elements = tuple(e for e in abundances if e in coupled)
        indices = {element: i for i, element in enumerate(elements)}
        abundance_ratios = {e: 10 ** (value - 12) for e, value in abundances.items()}
        ratios = [abundance_ratios[e] for e in elements]
        e_index = len(elements)
        kt = K_B * temperature
        index = min(bisect_right(self.temperatures, temperature) - 1, len(self.temperatures) - 2)
        low, high = self.temperatures[index:index + 2]
        fraction = math.log(temperature / low) / math.log(high / low)
        inverse_fraction = (1 / temperature - 1 / low) / (1 / high - 1 / low)

        def partition(values):
            return math.exp((1 - fraction) * math.log(values[index])
                            + fraction * math.log(values[index + 1]))

        # Each species stores its log mass-action prefactor and exponents in
        # free neutral atom and electron densities.
        species = []
        elemental_terms = [[] for _ in elements]
        positive = []
        negative = []
        atomic_only = []
        q = {}
        masses = {}

        def add(name, nuclei, powers, prefactor=0.0, charge=0):
            k = len(species)
            species.append((name, powers, prefactor))
            for element, coefficient in Counter(nuclei).items():
                elemental_terms[indices[element]].append((k, coefficient))
            if charge > 0:
                positive.append((k, charge))
            elif charge < 0:
                negative.append((k, -charge))

        for element in abundances:
            atom = self.atoms[element]
            neutral, ion = element, element + "+"
            q[neutral] = partition(atom["neutral_q"])
            q[ion] = partition(atom["ion_q"])
            masses[neutral] = masses[ion] = atom["mass_amu"]
            saha = (math.log(2 * SAHA * temperature**1.5 * q[ion] / q[neutral])
                    - atom["ionization_ev"] / (K_EV * temperature))
            if element in indices:
                i = indices[element]
                add(neutral, (element,), ((i, 1),))
                add(ion, (element,), ((i, 1), (e_index, -1)), saha, charge=1)
            else:
                atomic_only.append((neutral, ion, abundance_ratios[element], saha))

        for name, data in molecules.items():
            nuclei = data["atoms"]
            q[name] = partition(data["q"])
            masses[name] = sum(masses[e] for e in nuclei)
            if temperature > data.get("maximum_formation_temperature", math.inf):
                continue
            kp = data["log10_kp_cgs"]
            log_kp = (1 - inverse_fraction) * kp[index] + inverse_fraction * kp[index + 1]
            powers = tuple((indices[e], n) for e, n in Counter(nuclei).items())
            add(name, nuclei, powers, (len(nuclei) - 1) * math.log(kt) - math.log(10) * log_kp)

        h = indices["H"]
        # LTE attachment of H + e, with ground-state H and H⁻ weights 2 and 1.
        log_attachment = (0.754209 / (K_EV * temperature)
                          - math.log(4 * SAHA * temperature**1.5))
        add("H-", ("H",), ((h, 1), (e_index, 1)), log_attachment, charge=-1)
        add("e-", (), ((e_index, 1),), charge=-1)

        mass_per_hydrogen = AMU * math.fsum(ratio*masses[e] for e, ratio in abundance_ratios.items())
        return (elements, e_index, ratios, abundance_ratios, kt, species, elemental_terms,
                positive, negative, atomic_only, q, masses, mass_per_hydrogen)

    def _solve_equilibrium(self, temperature, gas_pressure, abundances, tolerance, initial):
        if (not math.isfinite(temperature) or not self.temperatures[0] <= temperature <= self.temperatures[-1]
                or not math.isfinite(gas_pressure) or gas_pressure <= 0):
            raise ValueError("Temperature must lie in the chemistry grid and pressure must be positive")
        if (abundances.get("H") != 12 or not abundances.keys() <= self.atoms.keys()
                or not all(math.isfinite(value) for value in abundances.values())):
            raise ValueError("Supply finite log epsilon abundances for known elements, with H=12")

        (elements, e_index, ratios, abundance_ratios, kt, species, elemental_terms,
         positive, negative, atomic_only, q, masses, mass_per_hydrogen) = self._equilibrium_terms(
            temperature, tuple(abundances.items())
        )
        particles = gas_pressure / kt
        size = e_index + 2
        all_particles = tuple((i, 1) for i in range(len(species)))
        log_ratios = [math.log(ratio) for ratio in ratios]
        log_particles = math.log(particles)

        def system(x):
            populations = [math.exp(c + sum(x[i] * p for i, p in powers))
                           for _, powers, c in species]

            # Atoms outside the molecular network have an analytic Saha
            # solution at the trial electron and hydrogen nuclei densities.
            neutral_and_ion = []
            atomic_number = net_positive = positive_derivative = 0.0
            nuclei_density = math.exp(x[-1])
            for _, _, ratio, saha in atomic_only:
                atom_density = nuclei_density*ratio
                odds = math.exp(-abs(saha-x[e_index]))
                minority = odds/(1+odds)
                majority = 1/(1+odds)
                neutral, ionized = (minority, majority) if saha > x[e_index] else (majority, minority)
                neutral_and_ion.append((atom_density*neutral, atom_density*ionized))
                atomic_number += atom_density
                net_positive += atom_density*ionized
                positive_derivative -= atom_density*neutral*ionized

            def total(terms, amount=0.0, electron_derivative=0.0):
                derivative = [0.0] * size
                derivative[e_index] = electron_derivative
                derivative[-1] = amount
                for j, multiplicity in terms:
                    contribution = multiplicity * populations[j]
                    amount += contribution
                    for i, power in species[j][1]:
                        derivative[i] += power * contribution
                return math.log(amount), [d / amount for d in derivative]

            residual, jacobian = [], []
            for log_ratio, terms in zip(log_ratios, elemental_terms):
                amount, row = total(terms)
                residual.append(amount - log_ratio - x[-1])
                row[-1] = -1.0
                jacobian.append(row)
            pos, dpos = total(positive, net_positive, positive_derivative)
            neg, dneg = total(negative)
            residual.append(pos - neg)
            jacobian.append([a - b for a, b in zip(dpos, dneg)])
            number, row = total(all_particles, atomic_number)
            residual.append(number - log_particles)
            jacobian.append(row)
            return residual, jacobian, populations, neutral_and_ion

        nuclei_scale = particles / math.fsum(abundance_ratios.values())
        x = [math.log(r * nuclei_scale) for r in ratios]
        x.extend((math.log(particles * 1e-4), math.log(nuclei_scale)))
        if initial is not None:
            x = [math.log(initial.number_densities[e]) for e in elements]
            x.extend((math.log(initial.electron_density),
                      math.log(initial.mass_density/mass_per_hydrogen)))
        residual, jacobian, populations, atomic_populations = system(x)
        for _ in range(60):
            error = max(map(abs, residual))
            if error < tolerance:
                numbers = dict.fromkeys(q, 0.0)
                numbers.update((entry[0], n) for entry, n in zip(species, populations))
                for (neutral, ion, _, _), (n, c) in zip(atomic_only, atomic_populations):
                    numbers[neutral], numbers[ion] = n, c
                density = math.exp(x[-1]) * mass_per_hydrogen
                return GasState(temperature, numbers["e-"], density, numbers, q.copy(), masses.copy(), error)
            step = _linear_solve(jacobian, [-r for r in residual])
            scale = min(1.0, 4.0 / max(map(abs, step)))
            for _ in range(20):
                trial = [a + scale * b for a, b in zip(x, step)]
                next_residual, next_jacobian, next_populations, next_atomic = system(trial)
                if max(map(abs, next_residual)) < error:
                    x = trial
                    residual, jacobian, populations = next_residual, next_jacobian, next_populations
                    atomic_populations = next_atomic
                    break
                scale *= 0.5
            else:
                raise RuntimeError(f"LTE chemistry step stalled at residual {error:.3g}")
        raise RuntimeError("LTE chemistry did not converge")


def _linear_solve(matrix, rhs):
    """Pivoted elimination for the small equilibrium Jacobian."""
    size = len(rhs)
    for i in range(size):
        pivot = max(range(i, size), key=lambda j: abs(matrix[j][i]))
        matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
        rhs[i], rhs[pivot] = rhs[pivot], rhs[i]
        diagonal = matrix[i][i]
        if abs(diagonal) < 1e-15:
            raise RuntimeError("Singular chemical-equilibrium Jacobian")
        columns = [k for k in range(i + 1, size) if matrix[i][k]]
        for j in range(i + 1, size):
            if matrix[j][i] == 0:
                continue
            factor = matrix[j][i] / diagonal
            for k in columns:
                matrix[j][k] -= factor * matrix[i][k]
            rhs[j] -= factor * rhs[i]
    for i in range(size - 1, -1, -1):
        rhs[i] = (rhs[i] - sum(matrix[i][j] * rhs[j] for j in range(i + 1, size))) / matrix[i][i]
    return rhs
