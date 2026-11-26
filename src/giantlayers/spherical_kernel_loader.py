"""Optional ctypes interface to the standard-C spherical wavelength kernel."""

import ctypes as ct
from pathlib import Path


class ScatteringIterationLimit(RuntimeError):
    """The compiled iteration exhausted its budget at one wavelength."""


class CompiledKernel:
    def __init__(self, path, rays):
        library = ct.CDLL(str(Path(path).resolve()))
        self.solve = library.spherical_solve
        integers = ct.POINTER(ct.c_int)
        doubles = ct.POINTER(ct.c_double)
        self.solve.argtypes = [ct.c_int, ct.c_int, integers, integers, *([doubles] * 8),
                               ct.c_double, ct.c_int, doubles, doubles, doubles]
        self.solve.restype = ct.c_int
        self.ray_count = len(rays)
        offsets = [0]
        scale = []
        for ray in rays:
            scale.extend((*ray.chord_scale, 0))
            offsets.append(len(scale))
        self.offsets = (ct.c_int * len(offsets))(*offsets)
        self.core = (ct.c_int * len(rays))(*(int(ray.core) for ray in rays))
        self.scale = (ct.c_double * len(scale))(*scale)
        self.mu = (ct.c_double * len(rays))(*(ray.mu_bottom for ray in rays))
        self.mean_weights = (ct.c_double * len(scale))(*(x for ray in rays for x in ray.mean_weights))
        self.flux_weights = (ct.c_double * len(scale))(*(x for ray in rays for x in ray.flux_weights))
        self.count = max(len(ray.mean_weights) for ray in rays)
        block = ct.c_double * self.count
        self.mean, self.flux, self.heating = block(), block(), block()

    def __call__(self, rays, widths, density, thermal, absorption, scattering,
                 tolerance=1e-7, max_iterations=200):
        extinction = [(a + s) * rho for a, s, rho in zip(absorption, scattering, density)]
        radial = [w * (a + b) / 2 for w, a, b in zip(widths, extinction, extinction[1:])]
        vector = ct.c_double * self.count
        iterations = self.solve(
            self.count, self.ray_count, self.offsets, self.core, self.scale, self.mu,
            self.mean_weights, self.flux_weights, (ct.c_double * len(radial))(*radial),
            vector(*thermal), vector(*absorption), vector(*scattering), tolerance, max_iterations,
            self.mean, self.flux, self.heating,
        )
        if iterations == -2:
            raise RuntimeError("The compiled spherical solver could not allocate ray workspace")
        if iterations < 0:
            raise ScatteringIterationLimit(f"Spherical scattering failed at limit {max_iterations}")
        return self.mean, self.flux, self.heating, iterations
