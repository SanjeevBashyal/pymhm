"""Data processing routines for mHM.

This package contains tools for converting, merging, regridding, and deriving
gridded data products. Submodules and public entry points are imported lazily to
keep import time low and avoid loading optional scientific dependencies unless
they are needed.
"""

import importlib

_MODULE_EXPORTS = {
    "bankfull": "bankfull",
    "difference": "difference",
    "long_term_mean": "long_term_mean",
    "merge": "merge",
    "ratio": "ratio",
    "relative_difference": "relative_difference",
}

_ATTR_EXPORTS = {
    "PlotOptions": ("difference", "PlotOptions"),
    "OutputOptions": ("difference", "OutputOptions"),
    "bankfull_discharge": ("bankfull", "bankfull_discharge"),
    "calc_diff": ("difference", "calc_diff"),
    "calc_ratio": ("ratio", "calc_ratio"),
    "calc_rel_diff": ("relative_difference", "calc_rel_diff"),
    "cal_long_term_mean": ("long_term_mean", "cal_long_term_mean"),
    "fill_dataarray_with_nearest": ("fill_nearest", "fill_dataarray_with_nearest"),
    "fill_nearest": ("fill_nearest", "fill_nearest"),
    "merge_files": ("merge", "merge_files"),
    "merge_files_from_folder": ("merge", "merge_files_from_folder"),
    "nearest_indices": ("fill_nearest", "nearest_indices"),
    "read_mask": ("fill_nearest", "read_mask"),
    "regrid": ("regrid", "regrid"),
    "regrid_file": ("regrid", "regrid_file"),
    "regrid_xarray": ("regrid", "regrid_xarray"),
}

__all__ = [
    "OutputOptions",
    "PlotOptions",
    "bankfull",
    "bankfull_discharge",
    "cal_long_term_mean",
    "calc_diff",
    "calc_ratio",
    "calc_rel_diff",
    "difference",
    "fill_dataarray_with_nearest",
    "fill_nearest",
    "long_term_mean",
    "merge",
    "merge_files",
    "merge_files_from_folder",
    "nearest_indices",
    "ratio",
    "read_mask",
    "regrid",
    "regrid_file",
    "regrid_xarray",
    "relative_difference",
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
