"""Independent initialization, equilibrium, and export checks."""

from dataclasses import replace
import math
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from giantlayers.collision_induced import HitranCIA
from giantlayers import initialize, restart_temperatures
from giantlayers.hminus import K_B, K_EV
from giantlayers.lte_chemistry import AMU, Thermochemistry
from giantlayers.lte_opacity import C_ANGSTROM_S, H_EV_S, profile_parameters, rayleigh_cross_sections, read_lines
from giantlayers.marcs_writer import read_document, write_document
from giantlayers.native_marcs import build_document
from giantlayers.spherical_heating import (
    evaluate, evaluate_parallel, log_temperature_curvature, surface_equilibrium_temperature,
)
from giantlayers.native_transfer import spectral_quadrature
from giantlayers.odf_transfer import planck_nu
from giantlayers.spherical_hydrostatic import reconstruct_for_transfer


class NativeAtmosphereChecks(unittest.TestCase):
    def test_collisional_damping_matches_lorentz_rate_at_layer_conditions(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/"lines.csv"
            path.write_text(
                "wavelength_vacuum_A,species,lower_energy_eV,log_gf,damping_s,"
                "electron_damping_cm3_s,neutral_h_damping_cm3_s,neutral_h_damping_exponent\n"
                "5000,Fe,1,-1,1e7,2e-6,3e-8,0.3\n"
                "5000,Fe,1,-1\n"
            )
            broadened, doppler = read_lines(path)
        gas = Thermochemistry().solve(3500, 100, {"H": 12, "Fe": 7.5})
        rate = 1e7 + 2e-6*gas.electron_density + 3e-8*gas.number_densities["H"]*.35**.3
        self.assertEqual(profile_parameters(broadened, gas, 2),
                         profile_parameters(replace(doppler, damping_s=rate), gas, 2))
        self.assertGreater(profile_parameters(broadened, gas, 2)[2], 0)

    def test_temperature_curvature_accounts_for_nonuniform_layer_spacing(self):
        axis = (0, 0.1, 1)
        linear = [(0, 0, 0, math.exp(8 + 0.2*x)) for x in axis]
        self.assertAlmostEqual(next(log_temperature_curvature(axis, linear)), 0)
        linear[1] = (0, 0, 0, math.exp(8.05))
        self.assertAlmostEqual(next(log_temperature_curvature(axis, linear)), 0.03)

    def test_deep_grid_power_concentrates_points_and_preserves_endpoints(self):
        options = dict(layers=80, deep_layers=32, column_bottom=1750, surface_pressure=0.5)
        abundance = {"H": 12}
        _, plain = initialize(3500, 0, 1, abundance, **options)
        _, concentrated = initialize(3500, 0, 1, abundance, deep_grid_power=2, **options)
        anchor = 80 - 1 - 32
        self.assertEqual(concentrated.column_mass[:anchor+1], plain.column_mass[:anchor+1])
        self.assertEqual(concentrated.column_mass[-1], plain.column_mass[-1])
        self.assertTrue(all(a < b for a,b in zip(concentrated.column_mass, concentrated.column_mass[1:])))
        self.assertGreater(sum(x > 1695 for x in concentrated.column_mass),
                           sum(x > 1695 for x in plain.column_mass))
        with self.assertRaisesRegex(ValueError, "coincident layers"):
            initialize(3500, 0, 1, abundance, deep_grid_power=1e6, **options)

    def test_deep_grid_refines_below_same_reference_column(self):
        _, original = initialize(3500, 1, 1, {"H": 12}, layers=56)
        _, refined = initialize(3500, 1, 1, {"H": 12}, layers=56, deep_layers=18)
        self.assertEqual(len(refined.column_mass)-1-refined.anchor, 18)
        self.assertEqual(refined.column_mass[refined.anchor], original.column_mass[original.anchor])
        self.assertEqual(refined.column_mass[0], original.column_mass[0])
        self.assertEqual(refined.column_mass[-1], original.column_mass[-1])
        self.assertTrue(all(a < b for a, b in zip(refined.column_mass, refined.column_mass[1:])))

    def test_parallel_wavelength_bands_preserve_flux_and_rosseland_mean(self):
        seed, boundary = initialize(3500, 1, 1, {"H": 12, "He": 10.93}, layers=7)
        model, _ = reconstruct_for_transfer(seed, boundary=boundary)
        wavelength = (1000., 3500., 7000., 15000., 40000., 120000.)
        serial_samples, parallel_samples = [], []
        serial_flux, parallel_flux = [], []
        serial = evaluate(model, wavelength, (), {}, core_rays=4, chemistry=boundary.chemistry,
                          native_convection=True, surface_mean_sample=lambda *x: serial_samples.append(x),
                          surface_sample=lambda *x: serial_flux.append(x))
        with patch("giantlayers.spherical_heating.MIN_WAVELENGTHS_PER_BAND", 2):
            parallel = evaluate_parallel(model, wavelength, (), {}, (), 4, 2, boundary.chemistry, None, 2,
                                         surface_mean_sample=lambda *x: parallel_samples.append(x),
                                         surface_sample=lambda *x: parallel_flux.append(x))
        self.assertEqual(serial_samples, parallel_samples)
        self.assertEqual(serial_flux, parallel_flux)
        for i in (0, 1, 6, 8, 9):
            for a, b in zip(serial[i], parallel[i]):
                self.assertTrue(math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12))
        for a, b in zip(serial[7], parallel[7]):
            self.assertTrue(math.isclose(a[4], b[4], rel_tol=1e-10))

    def test_surface_balance_recovers_blackbody_radiation_temperature(self):
        seed, boundary = initialize(3500, 1, 1, {"H": 12, "He": 10.93}, layers=7)
        model, _ = reconstruct_for_transfer(seed, boundary=boundary)
        wavelengths = (1000., 3500., 7000., 15000., 40000., 120000.)
        frequencies, _ = spectral_quadrature(wavelengths, "wavelength")
        field = tuple(planck_nu(nu, 2400) for nu in frequencies)
        temperature = surface_equilibrium_temperature(model, wavelengths, field, (), (), 2,
                                                      boundary.chemistry)
        self.assertAlmostEqual(temperature, 2400, delta=.01)

    def test_cold_dense_chemistry_preserves_pressure_and_charge(self):
        chemistry = Thermochemistry()
        abundances = {"H": 12, "He": 10.93, "C": 9.0, "N": 7.92, "O": 8.83,
                      "Fe": 7.50, "Ar": 6.40}
        gas = chemistry.solve(1000, 1e5, abundances)
        particles = math.fsum(gas.number_densities.values())
        positive = math.fsum(n for name, n in gas.number_densities.items() if name.endswith("+"))
        negative = gas.electron_density + gas.number_densities["H-"]
        self.assertAlmostEqual(particles*K_B*1000/1e5, 1.0, places=9)
        self.assertAlmostEqual(positive/negative, 1.0, places=9)

    def test_reused_chemistry_terms_keep_pressure_and_composition_independent(self):
        chemistry = Thermochemistry()
        reference = Thermochemistry()
        abundances = {"H": 12, "He": 10.93, "C": 8.4, "O": 8.7}
        first = chemistry.solve(2700, 10, abundances)
        first.partition_functions["H"] = -1
        first.masses_amu["H"] = -1
        for pressure, carbon in ((100, 8.4), (100, 9.1), (10, 8.4)):
            mixture = dict(abundances, C=carbon)
            self.assertEqual(chemistry.solve(2700, pressure, mixture),
                             reference.solve(2700, pressure, mixture))

    def test_cia_interpolates_wavenumber_temperature_and_density(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/"pair.cia"
            records = []
            for t, a, b in ((2000, 1, 3), (3000, 2, 6)):
                records.extend((f"{'H2-H2':>20}{100:10.3f}{200:10.3f}{2:7d}{t:7.1f}",
                                f"100 {a:e}", f"200 {b:e}"))
            path.write_text("\n".join(records))
            table = HitranCIA(path)
        gas = SimpleNamespace(number_densities={"H2": 10.0}, mass_density=2.0, temperature=2500.0)
        opacity = table.mass_opacity(gas, 1e3, table.on_grid((1e8/150, 1e8/250)))
        self.assertEqual(opacity, [150.0, 0.0])

    def test_tio_dissociation_constant_is_in_cgs_pressure_units(self):
        chemistry = Thermochemistry()
        temperature = 3000
        index = chemistry.temperatures.index(temperature)
        titanium, oxygen = (chemistry.atoms[name] for name in ("Ti", "O"))
        ti_mass, o_mass = titanium["mass_amu"], oxygen["mass_amu"]
        reduced_mass = AMU * ti_mass * o_mass / (ti_mass + o_mass)
        partition_ratio = (titanium["neutral_q"][index] * oxygen["neutral_q"][index]
                           / chemistry.molecules["TiO"]["q"][index])
        thermal = K_B * temperature
        expected = (2*math.pi*reduced_mass*thermal/(6.62607015e-27**2))**1.5
        expected *= thermal * partition_ratio * math.exp(-6.869/(K_EV*temperature))
        gas = chemistry.solve(temperature, 1000, {"H": 12, "Ti": 5.02, "O": 8.83})
        densities = gas.number_densities
        actual = thermal * densities["Ti"] * densities["O"] / densities["TiO"]
        self.assertAlmostEqual(actual/expected, 1.0, delta=0.003)

    def test_polyatomic_carbon_species_follow_mass_action(self):
        chemistry = Thermochemistry()
        values = {"H": 12, "C": 9.0, "O": 8.83, "N": 7.92}
        t = 3000
        gas = chemistry.solve(t, 1000, values)
        pressure = {name: n*K_B*t for name, n in gas.number_densities.items()}
        for species in ("HCN", "C2H2"):
            data = chemistry.molecules[species]
            ratio = math.prod(pressure[a] for a in data["atoms"])/pressure[species]
            self.assertAlmostEqual(math.log10(ratio), data["log10_kp_cgs"][3])

    def test_restart_depth_refinement_preserves_endpoints(self):
        seed, boundary = initialize(3500, 1, 1, {"H": 12}, layers=7)
        with TemporaryDirectory() as directory:
            path = Path(directory)/"restart.csv"
            path.write_text(f"layer,column_mass_g_cm2,temperature_K\n1,0,3000\n2,{boundary.column_mass[-1]},6000\n")
            with self.assertRaisesRegex(ValueError, "match the current"):
                restart_temperatures(seed, path)
            restarted = restart_temperatures(seed, path, interpolate=True)
        temperatures = [row[3] for row in restarted.structure]
        self.assertAlmostEqual(temperatures[0], 3000)
        self.assertAlmostEqual(temperatures[-1], 6000)
        self.assertTrue(all(a < b for a, b in zip(temperatures, temperatures[1:])))

    def test_rayleigh_cross_sections(self):
        def cross_sections(wavelength):
            return rayleigh_cross_sections(wavelength, H_EV_S*C_ANGSTROM_S/wavelength)

        h, h2, he = cross_sections(5000)
        self.assertAlmostEqual(h/1e-27, 1.026, places=3)
        self.assertAlmostEqual(h2/1e-27, 1.3884, places=4)
        self.assertGreater(he, 0)
        for red, near in zip(cross_sections(100000), cross_sections(50000)):
            self.assertAlmostEqual(near/red, 16, delta=0.02)
        self.assertEqual(cross_sections(500), (0, 0, 0))

    def test_warm_equilibrium_matches_cold_solution(self):
        chemistry = Thermochemistry()
        abundance = {"H": 12, "He": 10.93, "C": 9.0, "O": 8.83, "Fe": 7.50}
        seed = chemistry.solve(3499, 250.1, abundance)
        warm = chemistry.solve(3500, 250, abundance, initial=seed)
        cold = chemistry.solve(3500, 250, abundance)
        for a, b in zip(warm.number_densities.values(), cold.number_densities.values()):
            self.assertTrue(math.isclose(a, b, rel_tol=2e-9))
        self.assertTrue(math.isclose(warm.mass_density, cold.mass_density, rel_tol=2e-9))

    def test_native_structure_exports_calculated_depth_and_abundances(self):
        seed, boundary = initialize(3500, 1, 1, {"H": 12, "He": 10.93, "C": 9.0}, layers=5)
        model, _ = reconstruct_for_transfer(seed, boundary=boundary)
        # Native transfer normally provides these Rosseland and convection fields.
        convection = [[0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0]] * 5
        document = build_document(model, (None,) * 7 + (convection,), (), (), boundary, 2)
        with TemporaryDirectory() as directory:
            path = Path(directory) / "test.mod.gz"
            write_document(document, path)
            reread = read_document(path)
        self.assertEqual(reread.abundances[5], 9.0)
        self.assertEqual(len(reread.structure), 5)
        self.assertEqual(len(reread.partial_pressures), 3)
        self.assertEqual(reread.header[4][:5], "  2.0")
        self.assertEqual(reread.header[6][:12], "-99.00-99.00")
        self.assertGreater(float(reread.header[7][:12]), 0)
        self.assertTrue(all(a[0] < b[0] for a, b in zip(reread.structure, reread.structure[1:])))
        self.assertTrue(all(row[5] > 0 for row in reread.structure))


if __name__ == "__main__":
    unittest.main()
