"""Visualisation routines for mHM analysis outputs.

The package exposes plotting entry points lazily so importing the package does
not immediately import plotting backends.
"""

import importlib

_MODULE_EXPORTS = {
    "taylor_diagram": "taylor_diagram",
}

_ATTR_EXPORTS = {
    "calc_tim_mean": ("taylor_diagram", "calc_tim_mean"),
    "generate_taylor_diagram": ("taylor_diagram", "generate_taylor_diagram"),
    "mask_nan": ("taylor_diagram", "mask_nan"),
    "plot_2d_map": ("2d_map", "run"),
    "prepare_da": ("taylor_diagram", "prepare_da"),
}

__all__ = [
    "calc_tim_mean",
    "generate_taylor_diagram",
    "mask_nan",
    "plot_2d_map",
    "prepare_da",
    "taylor_diagram",
]


def __getattr__(name):
    if name in _ATTR_EXPORTS:
        module_name, attr_name = _ATTR_EXPORTS[name]
        module = importlib.import_module(f".{module_name}", __name__)
        attr = getattr(module, attr_name)
        globals()[name] = attr
        return attr
    if name in _MODULE_EXPORTS:
        module = importlib.import_module(f".{_MODULE_EXPORTS[name]}", __name__)
        globals()[name] = module
        return module
    error_msg = f"module '{__name__}' has no attribute '{name}'"
    raise AttributeError(error_msg)


def __dir__():
    return sorted(set(globals().keys()) | set(__all__))
