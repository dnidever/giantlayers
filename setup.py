"""Build the optional spherical transfer kernel."""

import os
import sys

from setuptools import Extension, setup


kernel = Extension(
    "giantlayers._spherical_kernel",
    [
        "src/giantlayers/kernel_module.c",
        "src/giantlayers/spherical_kernel.c",
    ],
    define_macros=[
        ("Py_LIMITED_API", "0x03090000"),
        ("GIANTLAYERS_EXTENSION", "1"),
    ],
    libraries=[] if sys.platform == "win32" else ["m"],
    py_limited_api=True,
    optional=True,
)

setup(ext_modules=[] if os.environ.get("GIANTLAYERS_PURE_PYTHON") == "1" else [kernel])
