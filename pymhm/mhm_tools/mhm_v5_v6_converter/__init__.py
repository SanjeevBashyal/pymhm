"""Converters for migrating mHM v5 setup inputs to mHM v6 structures.

Imports are resolved lazily to keep importing :mod:`mhm_tools` lightweight.
"""

import importlib

_MODULE_EXPORTS = {
    "landcover_ascii_to_nc": "landcover_ascii_to_nc",
}

_ATTR_EXPORTS = {
    "add_time_bounds_cf": ("landcover_ascii_to_nc", "add_time_bounds_cf"),
    "convert_lc_ascii_to_nc": ("landcover_ascii_to_nc", "convert_lc_ascii_to_nc"),
    "parse_nml_for_landcover": (
        "landcover_ascii_to_nc",
        "parse_nml_for_landcover",
    ),
}

__all__ = [
    "add_time_bounds_cf",
    "convert_lc_ascii_to_nc",
    "landcover_ascii_to_nc",
    "parse_nml_for_landcover",
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
