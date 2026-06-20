"""Evaluation routines for mHM outputs and reference data.

This package contains discharge, hydrograph, gridded-data, and run-overview
evaluation helpers. Submodules and public entry points are imported lazily to
keep package imports lightweight.
"""

import importlib

_MODULE_EXPORTS = {
    "discharge_evaluation": "discharge_evaluation",
    "hydrograph": "hydrograph",
    "mhm_run_overview": "mhm_run_overview",
}

_ATTR_EXPORTS = {
    "Catchment": ("hydrograph", "Catchment"),
    "EvalDataset": ("gridded_data_evaluation", "EvalDataset"),
    "Hydrograph": ("hydrograph", "Hydrograph"),
    "NmlStringAssignment": ("mhm_run_overview", "NmlStringAssignment"),
    "Objectives": ("hydrograph", "Objectives"),
    "collect_input_netcdf_files": (
        "mhm_run_overview",
        "collect_input_netcdf_files",
    ),
    "collect_output_flux_files": (
        "mhm_run_overview",
        "collect_output_flux_files",
    ),
    "compare_input_with_ref": ("gridded_data_evaluation", "compare_input_with_ref"),
    "create_mhm_run_overview": ("mhm_run_overview", "create_mhm_run_overview"),
    "evaludate_discharge_data": (
        "discharge_evaluation",
        "evaludate_discharge_data",
    ),
    "get_hydrograph_from_path": ("hydrograph", "get_hydrograph_from_path"),
    "gridded_data_evaluation": (
        "gridded_data_evaluation",
        "gridded_data_evaluation",
    ),
    "parse_namelist_string_assignments": (
        "mhm_run_overview",
        "parse_namelist_string_assignments",
    ),
    "plot_cdf": ("discharge_evaluation", "plot_cdf"),
    "plot_kde": ("discharge_evaluation", "plot_kde"),
    "plot_map": ("discharge_evaluation", "plot_map"),
    "resolve_namelist_path": ("mhm_run_overview", "resolve_namelist_path"),
}

__all__ = [
    "Catchment",
    "EvalDataset",
    "Hydrograph",
    "NmlStringAssignment",
    "Objectives",
    "collect_input_netcdf_files",
    "collect_output_flux_files",
    "compare_input_with_ref",
    "create_mhm_run_overview",
    "discharge_evaluation",
    "evaludate_discharge_data",
    "get_hydrograph_from_path",
    "gridded_data_evaluation",
    "hydrograph",
    "mhm_run_overview",
    "parse_namelist_string_assignments",
    "plot_cdf",
    "plot_kde",
    "plot_map",
    "resolve_namelist_path",
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
