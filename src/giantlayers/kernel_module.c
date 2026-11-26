/* Stable-ABI module marker. Python calls the numerical kernel through ctypes. */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_spherical_kernel",
    "Bundled GiantLayers spherical transfer kernel.",
    -1,
    NULL
};

PyMODINIT_FUNC PyInit__spherical_kernel(void) {
    return PyModule_Create(&module);
}
