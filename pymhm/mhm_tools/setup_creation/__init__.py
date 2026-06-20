"""Setup creation routines for mHM and mRM domains.

This package contains catchment delineation, setup cropping, PET calculation,
and forcing preparation utilities. Submodules and public entry points are
imported lazily to avoid loading optional geospatial dependencies at package
import time.
"""

import importlib

_MODULE_EXPORTS = {
    "catchment": "catchment",
    "latlon": "latlon",
    "pet_calc": "pet_calc",
    "prepare_mhm_forcings": "prepare_mhm_forcings",
}

_ATTR_EXPORTS = {
    "Catchment": ("catchment", "Catchment"),
    "Gauge": ("catchment", "Gauge"),
    "LatlonFiles": ("crop_mhm_setup", "LatlonFiles"),
    "calculate_pet": ("pet_calc", "calculate_pet"),
    "convert_units": ("prepare_mhm_forcings", "convert_units"),
    "create_catchment": ("catchment", "create_catchment"),
    "create_cell_area": ("catchment", "create_cell_area"),
    "create_latlon": ("latlon", "create_latlon"),
    "crop_file": ("crop_mhm_setup", "crop_file"),
    "crop_file_with_header": ("crop_mhm_setup", "crop_file_with_header"),
    "crop_mhm_setup": ("crop_mhm_setup", "crop_mhm_setup"),
    "e_rad_calculator": ("pet_calc", "e_rad_calculator"),
    "is_data_global": ("catchment", "is_data_global"),
    "merge_catchment": ("catchment", "merge_catchment"),
    "pet_calculator": ("pet_calc", "pet_calculator"),
    "prepare_forcings": ("prepare_mhm_forcings", "prepare_forcings"),
    "regrid_mask": ("crop_mhm_setup", "regrid_mask"),
    "validate_tmin_tmax": ("pet_calc", "validate_tmin_tmax"),
    "write_gauges_out": ("catchment", "write_gauges_out"),
    "xy_to_latlon": ("latlon", "xy_to_latlon"),
}

__all__ = [
    "Catchment",
    "Gauge",
    "LatlonFiles",
    "calculate_pet",
    "catchment",
    "convert_units",
    "create_catchment",
    "create_cell_area",
    "create_latlon",
    "crop_file",
    "crop_file_with_header",
    "crop_mhm_setup",
    "e_rad_calculator",
    "is_data_global",
    "latlon",
    "merge_catchment",
    "pet_calc",
    "pet_calculator",
    "prepare_forcings",
    "prepare_mhm_forcings",
    "regrid_mask",
    "validate_tmin_tmax",
    "write_gauges_out",
    "xy_to_latlon",
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
