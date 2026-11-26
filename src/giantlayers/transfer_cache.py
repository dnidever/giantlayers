"""Persist full-grid calibration residuals for an exactly matching restart."""

from hashlib import sha256
import json
import math
from pathlib import Path
import platform
import sys
from tempfile import NamedTemporaryFile


from .lte_chemistry import DEFAULT_TABLE


TRANSFER_MODULES = (
    "collision_induced", "gray_atmosphere", "hminus", "lte_chemistry", "lte_opacity",
    "marcs_writer", "mixing_length", "molecular_cross_section", "native_transfer",
    "odf_transfer", "spherical_heating", "spherical_hydrostatic",
    "spherical_kernel_loader", "spherical_transfer",
)


def encode(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False).encode("utf-8")


def input_fingerprint(opacity_paths):
    """Identify the runtime, transfer implementation, and input file contents."""
    directory = Path(__file__).resolve().parent
    paths = [Path(__file__), DEFAULT_TABLE, *(directory / f"{name}.py" for name in TRANSFER_MODULES),
             *opacity_paths]
    hashes = [sha256(path.read_bytes()).hexdigest() if path is not None else None for path in paths]
    return sha256(encode([sys.version, platform.platform(), platform.machine(), hashes])).hexdigest()


def transfer_key(fingerprint, model, wavelengths, core_rays, vmic, workers, force):
    return sha256(encode([fingerprint, model.header, model.abundances, model.structure,
                          model.thermodynamics, wavelengths, core_rays, vmic, workers, force])).hexdigest()


def read_residuals(path, key, layers):
    try:
        data = json.loads(path.read_bytes())
        payload = data["payload"]
        if payload["key"] != key or data["checksum"] != sha256(encode(payload)).hexdigest():
            return None
        flux, heating = payload["flux"], payload["heating"]
        if len(flux) == len(heating) == layers and all(map(math.isfinite, flux + heating)):
            return tuple(flux), tuple(heating)
    except (OSError, ValueError, TypeError, KeyError, OverflowError):
        pass
    return None


def write_residuals(path, key, flux, heating):
    payload = dict(key=key, flux=flux, heating=heating)
    data = encode(dict(payload=payload, checksum=sha256(encode(payload)).hexdigest()))
    temporary = None
    try:
        with NamedTemporaryFile(dir=path.parent, mode="wb", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(data)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
