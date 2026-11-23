"""Read the Castelli–Kurucz ATLAS9 BIG line-opacity distribution tables."""

from bisect import bisect_right
import math
from pathlib import Path
import struct

# Grid and weights from LINOP and BLOCKBIG in Castelli's ATLAS9 source:
# https://wwwuser.oats.inaf.it/fiorella.castelli/sources/atlas9g/atlas9mem.for
# LINOP overrides the last BLOCKBIG edge to 1,000,000 nm.
WAVELENGTH_EDGES_NM = (
    8.97666, 9.2, 9.5, 9.7173, 9.8159, 10.109, 10.3, 10.4645, 10.65, 10.8855, 11.2, 11.6, 11.9556,
    12.3, 12.6656, 12.7468, 12.935, 13.1603, 13.3172, 13.467, 13.8625, 14.1046, 14.5478, 14.9,
    15.3, 15.7487, 16.0035, 16.4, 16.8, 17.2509, 17.4326, 17.9366, 18.1011, 18.2444, 18.978,
    19.2242, 19.551, 20.1199, 20.2359, 20.8, 21.3, 21.9478, 22.0198, 22.5807, 22.7838, 23.2,
    23.6522, 24.3, 25.0074, 25.4661, 25.8924, 26.1432, 26.6, 27.1, 27.6, 28.0926, 28.6587, 29.3,
    29.9585, 30.2626, 31.3, 32.3, 33.3, 34.3, 35.2793, 36.141, 37, 38.016, 38.9664, 40, 41,
    41.8799, 42.5, 43.5, 44.7326, 45.5, 46.5, 47.5, 48.5018, 49.5, 50.4259, 50.8763, 52, 54, 56,
    57.4061, 59.5, 61.5, 63.5, 65.5379, 67, 69, 71, 72.2401, 74, 76, 78, 80, 82, 84, 86, 88, 90,
    91.1753, 94, 98, 102, 106, 110.031, 113, 116, 120, 123.929, 128, 132, 136, 140, 144.434, 148,
    151.435, 157, 162.151, 167.403, 171, 175, 180, 184, 188, 193, 197.47, 202, 207.132, 210, 215,
    220, 225, 230, 235, 240, 245, 251.126, 255, 260, 265, 270, 275, 280, 285, 290, 300, 310, 320,
    330, 340, 350, 360, 364.702, 370, 380, 390, 400, 410, 420, 430, 450, 460, 470, 480, 490, 500,
    510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690,
    700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 810, 820.583, 830, 840, 850, 860, 870,
    880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000, 1025, 1050, 1075, 1100,
    1125, 1150, 1175, 1200, 1225, 1250, 1275, 1300, 1325, 1350, 1375, 1400, 1425, 1458.82, 1475,
    1500, 1525, 1550, 1575, 1600, 1640, 1680, 1720, 1760, 1800, 1840, 1880, 1920, 1960, 2000,
    2050, 2100, 2150, 2200, 2250, 2279.4, 2300, 2350, 2400, 2450, 2500, 2550, 2600, 2650, 2700,
    2750, 2800, 2850, 2900, 2950, 3000, 3050, 3100, 3150, 3200, 3282.34, 3400, 3500, 3600, 3700,
    3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200,
    5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6600, 6800, 7000,
    7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 1e+06,
)
LOG_T = tuple(
    x / 100
    for x in (*range(330, 371, 2), *range(373, 401, 3), *range(405, 531, 5))
)
LOG_P = tuple(x / 2 for x in range(-8, 17))

# The binary values are stored in ascending opacity order. LINOP visits
# them in reverse order and starts at the smallest weight.
RANK_WEIGHTS = (0.1,) * 9 + (3 / 60, 2 / 60, 1 / 60)
C_NM_PER_S = 2.99792458e17
STEPS = struct.Struct("<12h")
MARKER = struct.Struct("<I")
PAYLOAD_BYTES = len(LOG_P) * STEPS.size
RECORD_BYTES = PAYLOAD_BYTES + 2 * MARKER.size


def quadrature():
    """Yield wavelength bounds in nm and the 12 frequency weights in Hz."""
    for short, long in zip(WAVELENGTH_EDGES_NM, WAVELENGTH_EDGES_NM[1:]):
        bandwidth = C_NM_PER_S * (1 / short - 1 / long)
        yield short, long, tuple(bandwidth * weight for weight in RANK_WEIGHTS)


def _bracket(grid, value, name):
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    log_value = math.log10(value)
    if not grid[0] <= log_value <= grid[-1]:
        raise ValueError(f"log10({name}) outside the ATLAS9 grid: {log_value:g}")
    i = min(bisect_right(grid, log_value) - 1, len(grid) - 2)
    return i, (log_value - grid[i]) / (grid[i + 1] - grid[i])


class Atlas9ODF:
    """Little-endian BIG ODF with unformatted Fortran record markers."""

    def __init__(self, path):
        self._data = Path(path).read_bytes()
        records = (len(WAVELENGTH_EDGES_NM) - 1) * len(LOG_T)
        if len(self._data) != records * RECORD_BYTES:
            raise ValueError("Unexpected size for a 328-band ATLAS9 BIG ODF")
        for start in range(0, len(self._data), RECORD_BYTES):
            head = MARKER.unpack_from(self._data, start)[0]
            tail = MARKER.unpack_from(self._data, start + MARKER.size + PAYLOAD_BYTES)[0]
            if head != PAYLOAD_BYTES or tail != PAYLOAD_BYTES:
                raise ValueError(f"Invalid ATLAS9 record markers at byte {start}")

    def _packed(self, band, t, p):
        record = band * len(LOG_T) + t
        start = record * RECORD_BYTES + MARKER.size + p * STEPS.size
        return STEPS.unpack_from(self._data, start)

    def line_opacity(self, temperature, pressure):
        """Yield one tuple of 12 line mass opacities in cm²/g per band."""
        it, ft = _bracket(LOG_T, temperature, "temperature")
        ip, fp = _bracket(LOG_P, pressure, "gas pressure")
        for band in range(len(WAVELENGTH_EDGES_NM) - 1):
            corners = (
                self._packed(band, it, ip),
                self._packed(band, it, ip + 1),
                self._packed(band, it + 1, ip),
                self._packed(band, it + 1, ip + 1),
            )
            yield tuple(
                10 ** (0.001 * ((1 - ft) * ((1 - fp) * a + fp * b)
                                + ft * ((1 - fp) * c + fp * d)))
                for a, b, c, d in zip(*corners)
            )
