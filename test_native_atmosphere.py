"""Independent initialization, equilibrium, and export checks."""

import math
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from collision_induced import HitranCIA
from giantlayers import initialize, restart_temperatures
from hminus import K_B
from lte_chemistry import Thermochemistry
from lte_opacity import C_ANGSTROM_S, H_EV_S, rayleigh_cross_sections
from marcs_writer import read_document, write_document
from native_marcs import build_document
from spherical_heating import evaluate, evaluate_parallel
from spherical_hydrostatic import reconstruct_for_transfer


class NativeAtmosphereChecks(unittest.TestCase):
    def test_parallel_wavelength_bands_preserve_flux_and_rosseland_mean(self):
        seed, boundary = initialize(3500, 1, 1, {"H": 12, "He": 10.93}, layers=7)
        model, _ = reconstruct_for_transfer(seed, boundary=boundary)
        wavelength = (1000., 3500., 7000., 15000., 40000., 120000.)
        serial = evaluate(model, wavelength, (), {}, core_rays=4, chemistry=boundary.chemistry,
                          native_convection=True)
        parallel = evaluate_parallel(model, wavelength, (), {}, (), 4, 2, boundary.chemistry, None, 2)
        for i in (0, 1, 6, 8, 9):
            for a, b in zip(serial[i], parallel[i]):
                self.assertTrue(math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12))
        for a, b in zip(serial[7], parallel[7]):
            self.assertTrue(math.isclose(a[4], b[4], rel_tol=1e-10))

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
