"""Keep cached calibration data exact and tied to the transfer inputs."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from giantlayers.transfer_cache import input_fingerprint, read_residuals, transfer_key, write_residuals


class TransferCacheTests(unittest.TestCase):
    def test_mismatched_or_modified_cache_is_a_miss(self):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "transfer.json"
            write_residuals(path, "valid", [0.1, -0.25], [0.25, 0.5])
            self.assertEqual(read_residuals(path, "valid", 2), ((0.1, -0.25), (0.25, 0.5)))
            self.assertIsNone(read_residuals(path, "different", 2))
            self.assertIsNone(read_residuals(path, "valid", 3))
            altered = json.loads(path.read_text())
            altered["payload"]["flux"][0] = 0
            path.write_text(json.dumps(altered))
            self.assertIsNone(read_residuals(path, "valid", 2))

    def test_fingerprint_tracks_input_content_and_transfer_grid(self):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "lines.dat"
            path.write_text("A")
            original = input_fingerprint([path])
            path.write_text("B")
            self.assertNotEqual(original, input_fingerprint([path]))
        model = SimpleNamespace(header=(), abundances=(), structure=((3500,),), thermodynamics=())
        original_key = transfer_key(original, model, (1000, 2000), 24, 2, 4, ())
        self.assertNotEqual(original_key, transfer_key(original, model, (1000, 2001), 24, 2, 4, ()))
        self.assertNotEqual(original_key, transfer_key(original, model, (1000, 2000), 24, 2, 1, ()))
