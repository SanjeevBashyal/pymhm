"""Legacy tools retained for backwards-compatible workflows.

This package contains tools that are still available but no longer belong to
the main public package groups. Imports are resolved lazily following the old
``pre`` and ``post`` package pattern.
"""

import importlib

_MODULE_EXPORTS = {
    "create_mhm_restart_file": "create_mhm_restart_file",
    "subdomain_masks": "subdomain_masks",
}

_ATTR_EXPORTS = {
    "CreateSubdomainMasks": ("subdomain_masks", "CreateSubdomainMasks"),
    "Grid": ("create_mhm_restart_file", "Grid"),
    "LatLon": ("create_mhm_restart_file", "LatLon"),
    "MHMRestartFile": ("create_mhm_restart_file", "MHMRestartFile"),
    "MPRRunner": ("create_mhm_restart_file", "MPRRunner"),
    "MorphFiles": ("create_mhm_restart_file", "MorphFiles"),
    "create_id_gauges": ("create_id_gauges", "create_id_gauges"),
    "create_subdomain_masks": ("subdomain_masks", "create_subdomain_masks"),
    "write_gauge_id": ("create_id_gauges", "write_gauge_id"),
}

__all__ = [
    "CreateSubdomainMasks",
    "Grid",
    "LatLon",
    "MHMRestartFile",
    "MPRRunner",
    "MorphFiles",
    "create_id_gauges",
    "create_mhm_restart_file",
    "create_subdomain_masks",
    "subdomain_masks",
    "write_gauge_id",
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
