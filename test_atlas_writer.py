import gzip
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class AtlasExportChecks(unittest.TestCase):
    def test_cli_exports_native_composition_and_fixed_width_deck(self):
        root = Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            composition = path / "abundances.json"
            composition.write_text('{"H": 12, "He": 11, "C": 8.7, "O": 8.3}')
            atlas = path / "trial.dat.gz"
            marcs = path / "trial.mod"
            completed = subprocess.run(
                [sys.executable, "-m", "giantlayers",
                 str(root / "examples/continuum-wavelengths.dat"), str(composition),
                 "--teff", "3500", "--logg", "-.5", "--mass", "1", "--layers", "8", "--steps", "0",
                 "--vmic", "2.5", "--output-atlas", str(atlas), "--output-model", str(marcs)],
                cwd=root, capture_output=True, text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(marcs.is_file())
            with gzip.open(atlas, "rt", encoding="ascii") as stream:
                records = stream.read().splitlines()

        self.assertEqual(records[0], "TEFF   3500.  GRAVITY-0.50000 LTE")
        self.assertIn("spherical LTE", records[1])
        header = next(i for i, line in enumerate(records) if line.startswith("READ DECK6"))
        self.assertEqual(int(records[header][10:13]), 8)
        hydrogen = float(records[4][44:52])
        carbon = float(records[5][50:57])
        oxygen = float(records[5][70:77])
        self.assertAlmostEqual(carbon - math.log10(hydrogen) + 12, 8.7, delta=0.0051)
        self.assertAlmostEqual(oxygen - math.log10(hydrogen) + 12, 8.3, delta=0.0051)
        previous_mass = 0.0
        for record in records[header+1:header+9]:
            self.assertEqual(len(record), 104)
            mass = float(record[:15])
            self.assertGreater(mass, previous_mass)
            previous_mass = mass
            columns = [float(record[24+i:34+i]) for i in range(0, 80, 10)]
            self.assertGreater(float(record[15:24]), 0)
            self.assertGreater(columns[1], 0)  # Electron number density.
            self.assertEqual(columns[4], 2.5e5)
            self.assertGreater(columns[7], 0)  # Sound speed.
        self.assertEqual(records[header+9][:5], "PRADK")
        self.assertGreater(float(records[header+9][5:]), 0)
