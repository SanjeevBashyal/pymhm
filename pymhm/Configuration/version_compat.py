# -*- coding: utf-8 -*-
"""Version-specific namelist value translations."""
from __future__ import annotations

import re
from datetime import datetime
from typing import Any

from ..project_layout import version_key as layout_version_key
from .domain_info import domain_count as current_domain_count
from .namelist import canonical_name
from .path_defaults import configuration_path_defaults


ConfigValues = dict[str, dict[str, Any]]


def render_values_for_version(
        version_text: str,
        kind: str,
        values: ConfigValues,
        dialog: Any | None = None) -> ConfigValues:
    """Return values in the block/variable layout expected by a template."""
    if kind == "mhm" and layout_version_key(version_text) == "v5.13":
        return v513_mhm_values(values, dialog)
    return values


def v513_mhm_values(values: ConfigValues, dialog: Any | None = None) -> ConfigValues:
    """Translate schema-driven mHM values to the mHM 5.13 namelist layout."""
    count = _configured_domain_count(values, dialog)
    merged = _merge_defaults(configuration_path_defaults(count), values)
    result: ConfigValues = {}

    _project_description(result, merged)
    _mainconfig(result, merged, count, dialog)
    _mainconfig_mhm_mrm(result, merged, count, dialog)
    _mainconfig_mrm(result, merged)
    _directories_general(result, merged, count)
    _directories_mhm(result, merged, count)
    _directories_mrm(result, merged, count)
    _process_selection(result, merged)
    _soil_lai_landcover(result, merged)
    _directories_mpr(result, merged, count)
    _time_periods(result, merged, count)
    _meteo_weights(result, merged)
    _optimization(result, merged)
    _finalize_indexed_values(result)
    return result


def _configured_domain_count(values: ConfigValues, dialog: Any | None) -> int:
    actual = 0
    if dialog is not None:
        try:
            actual = int(current_domain_count(dialog))
        except Exception:
            actual = 0
    configured = _value(values, "config_project", "n_domains", None)
    try:
        configured = int(configured)
    except (TypeError, ValueError):
        configured = 0
    return max(1, configured or actual or 1)


def _merge_defaults(defaults: ConfigValues, values: ConfigValues) -> ConfigValues:
    merged = {
        canonical_name(block): dict(block_values or {})
        for block, block_values in (defaults or {}).items()
    }
    for block, block_values in (values or {}).items():
        block_key = canonical_name(block)
        merged.setdefault(block_key, {})
        for key, value in (block_values or {}).items():
            if _is_blank(value) and canonical_name(key) in merged[block_key]:
                continue
            merged[block_key][canonical_name(key)] = value
    return merged


def _is_blank(value: Any) -> bool:
    return value is None or value == "" or value == []


def _block(values: ConfigValues, block: str) -> dict[str, Any]:
    block_key = canonical_name(block)
    if block_key in values:
        return values.get(block_key, {})
    for candidate, block_values in values.items():
        if canonical_name(candidate) == block_key and isinstance(block_values, dict):
            return {
                canonical_name(name): value
                for name, value in block_values.items()
            }
    return {}


def _value(
        values: ConfigValues,
        block: str,
        name: str,
        default: Any = None) -> Any:
    return _block(values, block).get(canonical_name(name), default)


def _domain_value(
        values: ConfigValues,
        block: str,
        name: str,
        index: int,
        default: Any = None) -> Any:
    block_values = _block(values, block)
    key = canonical_name(name)
    indexed_key = f"{key}__domain{index}"
    if indexed_key in block_values:
        return block_values[indexed_key]
    if key not in block_values:
        return default
    value = block_values[key]
    if isinstance(value, tuple):
        value = list(value)
    if isinstance(value, dict) and "__indexed__" in value:
        indexed = value.get("__indexed__", {})
        return indexed.get(str(index), indexed.get(index, default))
    if isinstance(value, list):
        if value and isinstance(value[0], list):
            return value[index - 1] if index <= len(value) else value[-1]
        return value[index - 1] if index <= len(value) else (
            value[-1] if value else default)
    return value


def _domain_value_any(
        values: ConfigValues,
        block: str,
        names: tuple[str, ...],
        index: int,
        default: Any = None) -> Any:
    for name in names:
        value = _domain_value(values, block, name, index, None)
        if not _is_blank(value):
            return value
    return default


def _first_domain(
        values: ConfigValues,
        block: str,
        name: str,
        default: Any = None) -> Any:
    return _domain_value(values, block, name, 1, _value(values, block, name, default))


def _put(result: ConfigValues, block: str, name: str, value: Any) -> None:
    if value is None:
        return
    result.setdefault(block, {})[name] = value


def _put_index(
        result: ConfigValues,
        block: str,
        name: str,
        index: int,
        value: Any) -> None:
    if value is None:
        return
    block_values = result.setdefault(block, {})
    indexed = block_values.setdefault(name, {"__indexed__": {}})
    indexed["__indexed__"][index] = value


def _finalize_indexed_values(result: ConfigValues) -> None:
    """Replace temporary indexed mappings with converter-ready arrays."""
    for block_values in result.values():
        for name, value in list(block_values.items()):
            if not isinstance(value, dict) or "__indexed__" not in value:
                continue
            indexed = value["__indexed__"]
            indices = sorted(indexed)
            if indices != list(range(1, len(indices) + 1)):
                raise ValueError(f"indexed values for '{name}' must be contiguous")
            block_values[name] = [indexed[index] for index in indices]


def _path_text(value: Any, fallback: str = "") -> str:
    if _is_blank(value):
        value = fallback
    if isinstance(value, (list, tuple)):
        value = value[0] if value else fallback
    return str(value or fallback).strip().replace("\\", "/")


def _dir_from_path(value: Any, fallback: str) -> str:
    path = _path_text(value, fallback)
    if not path:
        return ""
    if path.endswith("/"):
        return path
    if "/" not in path:
        return path + "/"
    name = path.rsplit("/", 1)[-1]
    if "." not in name:
        return path + "/"
    return path.rsplit("/", 1)[0] + "/"


def _basename(value: Any, fallback: str) -> str:
    path = _path_text(value, fallback)
    if not path:
        return fallback
    return path.rstrip("/").rsplit("/", 1)[-1] or fallback


def _dialog_resolution(dialog: Any | None, method_name: str) -> float | None:
    if dialog is None or not hasattr(dialog, method_name):
        return None
    try:
        value = getattr(dialog, method_name)()
    except Exception:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if value > 0 else None


def _resolution(
        values: ConfigValues,
        block: str,
        name: str,
        index: int,
        dialog_value: float | None) -> Any:
    if dialog_value:
        return dialog_value
    value = _domain_value(values, block, name, index, None)
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if numeric > 0 else None


def _resolution_any(
        values: ConfigValues,
        block: str,
        names: tuple[str, ...],
        index: int,
        dialog_value: float | None) -> Any:
    if dialog_value:
        return dialog_value
    for name in names:
        value = _domain_value(values, block, name, index, None)
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if numeric > 0:
            return numeric
    return None


def _any_domain_bool(
        values: ConfigValues,
        block: str,
        name: str,
        count: int,
        default: bool = False) -> bool:
    found = False
    for index in range(1, count + 1):
        value = _domain_value(values, block, name, index, None)
        if value is None:
            continue
        found = True
        if bool(value):
            return True
    return default if not found else False


def _project_description(result: ConfigValues, values: ConfigValues) -> None:
    for name in (
            "project_details",
            "setup_description",
            "simulation_type",
            "Conventions",
            "contact",
            "mHM_details",
            "history"):
        _put(
            result,
            "project_description",
            name,
            _value(values, "config_project", name, None),
        )


def _mainconfig(
        result: ConfigValues,
        values: ConfigValues,
        count: int,
        dialog: Any | None) -> None:
    l1_resolution = _dialog_resolution(dialog, "current_l1_resolution")
    coord_system = _value(values, "config_project", "iFlag_cordinate_sys", None)
    if coord_system is None:
        coord_system = _first_domain(
            values,
            "config_coupling",
            "morph_grid_coordsys",
            _first_domain(values, "config_coupling", "hydro_grid_coordsys", 0),
        )
    _put(result, "mainconfig", "iFlag_cordinate_sys", coord_system)
    _put(result, "mainconfig", "nDomains", count)
    _put(
        result,
        "mainconfig",
        "write_restart",
        _any_domain_bool(values, "config_mhm", "write_restart", count, False),
    )
    for index in range(1, count + 1):
        _put_index(
            result,
            "mainconfig",
            "resolution_Hydrology",
            index,
            _resolution_any(
                values,
                "config_mhm",
                ("resolution_Hydrology", "resolution"),
                index,
                l1_resolution,
            ),
        )
        _put_index(
            result,
            "mainconfig",
            "L0Domain",
            index,
            _domain_value(values, "config_mhm", "L0Domain", index, index),
        )
        _put_index(
            result,
            "mainconfig",
            "read_opt_domain_data",
            index,
            _domain_value(
                values,
                "config_mhm",
                "read_opt_domain_data",
                index,
                0,
            ),
        )


def _mainconfig_mhm_mrm(
        result: ConfigValues,
        values: ConfigValues,
        count: int,
        dialog: Any | None) -> None:
    l11_resolution = _dialog_resolution(dialog, "current_l11_resolution")
    for index in range(1, count + 1):
        _put_index(
            result,
            "mainconfig_mhm_mrm",
            "mhm_file_RestartIn",
            index,
            _path_text(_domain_value_any(
                values,
                "config_mhm",
                ("mhm_file_RestartIn", "restart_input_path"),
                index,
                "restart/mhm_restart_in.nc")),
        )
        _put_index(
            result,
            "mainconfig_mhm_mrm",
            "mrm_file_RestartIn",
            index,
            _path_text(
                _domain_value_any(
                    values,
                    "config_mhm",
                    ("mrm_file_RestartIn",),
                    index,
                    None,
                )
                or _domain_value(
                    values,
                    "config_mrm",
                    "restart_input_path",
                    index,
                    "restart/mrm_restart_in.nc",
                )
            ),
        )
        _put_index(
            result,
            "mainconfig_mhm_mrm",
            "resolution_Routing",
            index,
            _resolution_any(
                values,
                "config_mrm",
                ("resolution_Routing", "resolution"),
                index,
                l11_resolution,
            ),
        )
    _put(
        result,
        "mainconfig_mhm_mrm",
        "timestep",
        _first_domain(values, "config_time", "timestep", 1),
    )
    _put(
        result,
        "mainconfig_mhm_mrm",
        "read_restart",
        _any_domain_bool(values, "config_mhm", "read_restart", count, False),
    )
    for name in ("optimize", "optimize_restart", "opti_method", "opti_function"):
        _put(
            result,
            "mainconfig_mhm_mrm",
            name,
            _value(values, "config_optimize", name, None),
        )


def _mainconfig_mrm(result: ConfigValues, values: ConfigValues) -> None:
    _put(
        result,
        "mainconfig_mrm",
        "gw_coupling",
        _value(values, "config_mrm", "gw_coupling", None),
    )


def _directories_general(
        result: ConfigValues,
        values: ConfigValues,
        count: int) -> None:
    first_morph_dir = _dir_from_path(
        _domain_value(values, "config_input", "dem_path", 1, "data/static/morph/dem.asc"),
        "data/static/morph/",
    )
    _put(result, "directories_general", "dirConfigOut", "output/")
    _put(result, "directories_general", "dirCommonFiles", first_morph_dir)
    for index in range(1, count + 1):
        morph_dir = _dir_from_path(
            _domain_value(
                values,
                "config_input",
                "dem_path",
                index,
                "data/static/morph/dem.asc"),
            "data/static/morph/",
        )
        land_cover = _domain_value(
            values,
            "config_mpr",
            "land_cover_path",
            index,
            _domain_value(values, "config_input", "lai_class_path", index, "data/static/morph/lc.asc"),
        )
        output_dir = _dir_from_path(
            _domain_value(values, "config_mhm", "output_path", index, "output/mhm_output.nc"),
            "output/",
        )
        _put_index(result, "directories_general", "dir_Morpho", index, morph_dir)
        _put_index(
            result,
            "directories_general",
            "dir_LCover",
            index,
            _dir_from_path(land_cover, morph_dir),
        )
        _put_index(
            result,
            "directories_general",
            "mhm_file_RestartOut",
            index,
            _path_text(_domain_value_any(
                values,
                "config_mhm",
                ("mhm_file_RestartOut", "restart_output_path"),
                index,
                "restart/mhm_restart_out.nc")),
        )
        _put_index(
            result,
            "directories_general",
            "mrm_file_RestartOut",
            index,
            _path_text(
                _domain_value_any(
                    values,
                    "config_mhm",
                    ("mrm_file_RestartOut",),
                    index,
                    None,
                )
                or _domain_value(
                    values,
                    "config_mrm",
                    "restart_output_path",
                    index,
                    "restart/mrm_restart_out.nc",
                )
            ),
        )
        _put_index(result, "directories_general", "dir_Out", index, output_dir)
        _put_index(
            result,
            "directories_general",
            "file_LatLon",
            index,
            _path_text(_domain_value(
                values,
                "config_input",
                "latlon_path",
                index,
                "data/latlon.nc")),
        )


def _directories_mhm(
        result: ConfigValues,
        values: ConfigValues,
        count: int) -> None:
    _put(result, "directories_mHM", "inputFormat_meteo_forcings", "nc")
    for index in range(1, count + 1):
        meteo_defaults = {
            "pre_path": "data/meteo/pre/pre.nc",
            "temp_path": "data/meteo/tavg/tavg.nc",
            "pet_path": "data/meteo/pet/pet.nc",
            "tmin_path": "data/meteo/tmin/tmin.nc",
            "tmax_path": "data/meteo/tmax/tmax.nc",
            "netrad_path": "data/meteo/netrad/netrad.nc",
            "eabs_path": "data/meteo/eabs/eabs.nc",
            "wind_path": "data/meteo/windspeed/windspeed.nc",
            "ssrd_path": "data/meteo/ssrd/ssrd.nc",
        }
        path_map = {
            "dir_Precipitation": "pre_path",
            "dir_Temperature": "temp_path",
            "dir_ReferenceET": "pet_path",
            "dir_MinTemperature": "tmin_path",
            "dir_MaxTemperature": "tmax_path",
            "dir_NetRadiation": "netrad_path",
            "dir_absVapPressure": "eabs_path",
            "dir_windspeed": "wind_path",
            "dir_Radiation": "ssrd_path",
        }
        for target, source in path_map.items():
            _put_index(
                result,
                "directories_mHM",
                target,
                index,
                _dir_from_path(
                    _domain_value(
                        values,
                        "config_input",
                        source,
                        index,
                        meteo_defaults[source]),
                    _dir_from_path(meteo_defaults[source], "data/meteo/"),
                ),
            )
        _put_index(
            result,
            "directories_mHM",
            "time_step_model_inputs",
            index,
            _domain_value(values, "config_input", "chunking", index, 0),
        )


def _directories_mrm(
        result: ConfigValues,
        values: ConfigValues,
        count: int) -> None:
    for index in range(1, count + 1):
        output_dir = _dir_from_path(
            _domain_value(values, "config_mhm", "output_path", index, "output/mhm_output.nc"),
            "output/",
        )
        _put_index(
            result,
            "directories_mRM",
            "dir_Gauges",
            index,
            "data/observation/streamflow/",
        )
        _put_index(result, "directories_mRM", "dir_Total_Runoff", index, output_dir)
        _put_index(
            result,
            "directories_mRM",
            "dir_Bankfull_Runoff",
            index,
            "data/static/morph/",
        )


def _process_selection(result: ConfigValues, values: ConfigValues) -> None:
    mapping = (
        (1, "interception"),
        (2, "snow"),
        (3, "soil_moisture"),
        (4, "direct_runoff"),
        (5, "pet"),
        (6, "interflow"),
        (7, "percolation"),
        (8, "routing"),
        (9, "baseflow"),
        (10, "neutrons"),
        (11, "temperature_routing"),
    )
    defaults = (1, 1, 1, 1, 0, 1, 1, 3, 1, 0, 0)
    process_cases = []
    for index, name in mapping:
        value = _value(values, "config_processes", name, None)
        if value is None:
            value = defaults[index - 1]
        if name == "pet" and value == -2:
            value = 0
        process_cases.append(value)
    _put(result, "processSelection", "processCase", process_cases)


def _soil_lai_landcover(result: ConfigValues, values: ConfigValues) -> None:
    _put(
        result,
        "soildata",
        "iFlag_soilDB",
        _first_domain(values, "config_mpr", "soil_db_mode", 0),
    )
    _put(
        result,
        "soildata",
        "tillageDepth",
        _first_domain(values, "config_mpr", "tillage_depth", None),
    )
    _put(
        result,
        "soildata",
        "nSoilHorizons_mHM",
        _first_domain(values, "config_mpr", "n_horizons", None),
    )
    soil_depth = _first_domain(values, "config_mpr", "soil_depth", None)
    if isinstance(soil_depth, list):
        for index, depth in enumerate(soil_depth, start=1):
            _put_index(result, "soildata", "soil_Depth", index, depth)
    _put(
        result,
        "LAI_data_information",
        "timeStep_LAI_input",
        _first_domain(values, "config_mpr", "lai_time_step", None),
    )
    _put(result, "LAI_data_information", "inputFormat_gridded_LAI", "nc")
    _put(
        result,
        "LCover_MPR",
        "fracSealed_cityArea",
        _first_domain(values, "config_mpr", "fracSealed_cityArea", None),
    )
    land_cover_path = _first_domain(
        values,
        "config_mpr",
        "land_cover_path",
        _first_domain(values, "config_input", "lai_class_path", "data/static/morph/lc.asc"),
    )
    periods = _value(values, "config_time", "eval_Per", None)
    first_period = periods[0] if isinstance(periods, list) and periods else {}
    last_period = periods[-1] if isinstance(periods, list) and periods else {}
    start_year = _component_value(first_period, "yStart") or (
        (_date_parts(_first_domain(values, "config_time", "sim_start", None)) or [None])[0]
        or (_date_parts(_first_domain(values, "config_time", "eval_start", None)) or [None])[0]
        or 1900
    )
    end_year = _component_value(last_period, "yEnd") or (
        (_date_parts(_first_domain(values, "config_time", "sim_end", None)) or [None])[0]
        or 2100
    )
    _put(result, "LCover", "nLCoverScene", 1)
    _put_index(result, "LCover", "LCoverYearStart", 1, start_year)
    _put_index(result, "LCover", "LCoverYearEnd", 1, end_year)
    _put_index(
        result,
        "LCover",
        "LCoverfName",
        1,
        _basename(land_cover_path, "lc.asc"),
    )
    evap_coeff = _first_domain(values, "config_mhm", "evap_coeff", None)
    _put(result, "panEvapo", "evap_coeff", evap_coeff)


def _directories_mpr(
        result: ConfigValues,
        values: ConfigValues,
        count: int) -> None:
    for index in range(1, count + 1):
        _put_index(
            result,
            "directories_MPR",
            "dir_gridded_LAI",
            index,
            _dir_from_path(
                _domain_value(values, "config_mpr", "lai_path", index, "data/lai/lai.nc"),
                "data/lai/",
            ),
        )


def _date_parts(value: Any) -> tuple[int, int, int] | None:
    if _is_blank(value):
        return None
    text = str(value).strip().replace("T", " ")
    match = re.match(r"^(\d{4})-(\d{1,2})-(\d{1,2})", text)
    if not match:
        return None
    return tuple(int(part) for part in match.groups())


def _parse_date(value: Any) -> datetime | None:
    parts = _date_parts(value)
    if not parts:
        return None
    try:
        return datetime(parts[0], parts[1], parts[2])
    except ValueError:
        return None


def _component_value(value: Any, name: str, default: Any = None) -> Any:
    """Return a derived component using forgiving input-name matching."""
    if not isinstance(value, dict):
        return default
    key = canonical_name(name)
    for component, component_value in value.items():
        if canonical_name(component) == key:
            return component_value
    return default


def _evaluation_periods(value: Any) -> list[dict[str, Any]]:
    """Normalize configured evaluation-period objects to exact components."""
    if not isinstance(value, list):
        return []
    names = ("yStart", "mStart", "dStart", "yEnd", "mEnd", "dEnd")
    periods = []
    for item in value:
        if not isinstance(item, dict):
            return []
        period = {
            name: _component_value(item, name)
            for name in names
            if _component_value(item, name) is not None
        }
        periods.append(period)
    return periods


def _time_periods(
        result: ConfigValues,
        values: ConfigValues,
        count: int) -> None:
    periods = _evaluation_periods(
        _value(values, "config_time", "eval_Per", None))
    warming = _value(values, "config_time", "warming_Days", None)
    if periods:
        if warming is not None:
            _put(result, "time_periods", "warming_Days", warming)
        _put(result, "time_periods", "eval_Per", periods)
        return

    starts: list[tuple[int, int, int] | None] = []
    ends: list[tuple[int, int, int] | None] = []
    warming_days: list[int | None] = []
    for index in range(1, count + 1):
        sim_start = _domain_value(values, "config_time", "sim_start", index, None)
        eval_start = _domain_value(values, "config_time", "eval_start", index, sim_start)
        sim_end = _domain_value(values, "config_time", "sim_end", index, None)
        starts.append(_date_parts(eval_start or sim_start))
        ends.append(_date_parts(sim_end))
        sim_dt = _parse_date(sim_start)
        eval_dt = _parse_date(eval_start)
        warming_days.append(
            max(0, (eval_dt - sim_dt).days) if sim_dt and eval_dt else None)

    if all(day is not None for day in warming_days):
        _put(result, "time_periods", "warming_Days", warming_days)

    if len([item for item in starts if item is not None]) != count:
        return
    if len([item for item in ends if item is not None]) != count:
        return
    _put(result, "time_periods", "eval_Per", [
        {
            "yStart": starts[index][0],
            "mStart": starts[index][1],
            "dStart": starts[index][2],
            "yEnd": ends[index][0],
            "mEnd": ends[index][1],
            "dEnd": ends[index][2],
        }
        for index in range(count)
    ])


def _meteo_weights(result: ConfigValues, values: ConfigValues) -> None:
    _put(
        result,
        "nightDayRatio",
        "read_meteo_weights",
        _first_domain(values, "config_meteo", "read_meteo_weights", None),
    )
    for source, target in (
            ("frac_night_pre", "fnight_prec"),
            ("frac_night_pet", "fnight_pet"),
            ("frac_night_temp", "fnight_temp"),
            ("frac_night_ssrd", "fnight_ssrd"),
            ("frac_night_strd", "fnight_strd")):
        _put(
            result,
            "nightDayRatio",
            target,
            _first_domain(values, "config_meteo", source, None),
        )


def _optimization(result: ConfigValues, values: ConfigValues) -> None:
    for name in (
            "nIterations",
            "seed",
            "dds_r",
            "sa_temp",
            "sce_ngs",
            "sce_npg",
            "sce_nps",
            "mcmc_opti",
            "mcmc_error_params"):
        source = "iterations" if canonical_name(name) == "niterations" else name
        _put(result, "Optimization", name, _value(values, "config_optimize", source, None))
    _put(
        result,
        "baseflow_config",
        "BFI_calc",
        _value(values, "config_optimize", "BFI_calc", None),
    )
