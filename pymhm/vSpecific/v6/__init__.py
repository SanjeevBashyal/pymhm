"""mHM v6 namelist initial values."""

from __future__ import annotations

from typing import Any

from ..common import (
    domain_count,
    geology_count,
    geoparameter,
    namelist_settings,
    repeat,
    resolution,
)


def build_dimensions(dialog: Any) -> dict[str, int]:
    return {
        "max_layers": 10,
        "n_domains": domain_count(dialog),
        "n_geo_units": geology_count(dialog),
    }


def build_initial_values(dialog: Any) -> dict[str, Any]:
    count = domain_count(dialog)
    settings = namelist_settings(dialog)
    geo_count = geology_count(dialog)
    l1 = resolution(dialog, "current_l1_resolution")
    l11 = resolution(dialog, "current_l11_resolution")
    resolution_values: dict[str, Any] = {}
    if l1 is not None:
        resolution_values["hydro"] = repeat(l1, count)
    if l11 is not None:
        resolution_values["route"] = repeat(l11, count)

    values: dict[str, Any] = {
        "main": {
            "config_project": {
                "n_domains": count,
                "n_geo_units": geo_count,
                "max_layers": 10,
            },
            "config_resolution": resolution_values,
            "config_input": {
                "latlon_path": repeat("data/latlon.nc", count),
                "pre_path": repeat("data/meteo/pre/pre.nc", count),
                "pet_path": repeat("data/meteo/pet/pet.nc", count),
                "temp_path": repeat("data/meteo/tavg/tavg.nc", count),
                "tann_path": repeat("data/meteo/tann/tann.nc", count),
                "tmin_path": repeat("data/meteo/tmin/tmin.nc", count),
                "tmax_path": repeat("data/meteo/tmax/tmax.nc", count),
                "ssrd_path": repeat("data/meteo/ssrd/ssrd.nc", count),
                "strd_path": repeat("data/meteo/strd/strd.nc", count),
                "netrad_path": repeat("data/meteo/netrad/netrad.nc", count),
                "eabs_path": repeat("data/meteo/eabs/eabs.nc", count),
                "wind_path": repeat("data/meteo/windspeed/windspeed.nc", count),
                "hydro_mask_path": repeat("data/static/morph/dem.asc", count),
                "dem_path": repeat("data/static/morph/dem.asc", count),
                "slope_path": repeat("data/static/morph/slope.asc", count),
                "aspect_path": repeat("data/static/morph/aspect.asc", count),
                "fdir_path": repeat("data/static/morph/fdir.asc", count),
                "facc_path": repeat("data/static/morph/facc.asc", count),
                "geo_class_path": repeat("data/static/morph/geology_class.asc", count),
                "soil_class_path": repeat("data/static/morph/soil_class.asc", count),
                "lai_class_path": repeat("data/static/morph/lc.asc", count),
                "morph_mask_path": repeat("data/static/morph/dem.asc", count),
            },
            "config_mpr": {
                "land_cover_path": repeat("data/static/morph/lc.asc", count),
                "lai_path": repeat("data/lai/lai.nc", count),
                "soil_lut_path": repeat(
                    "data/static/morph/soil_classdefinition.txt", count
                ),
                "geo_lut_path": repeat(
                    "data/static/morph/geology_classdefinition.txt", count
                ),
                "lai_lut_path": repeat(
                    "data/static/morph/LAI_classdefinition.txt", count
                ),
                "restart_input_path": repeat("restart/mpr_restart_in.nc", count),
                "restart_output_path": repeat("restart/mpr_restart_out.nc", count),
            },
            "config_mhm": {
                "output_path": repeat("output/mhm_output.nc", count),
                "restart_input_path": repeat("restart/mhm_restart_in.nc", count),
                "restart_output_path": repeat("restart/mhm_restart_out.nc", count),
            },
            "config_mrm": {
                "scc_gauges_path": repeat("data/static/morph/idgauges.asc", count),
                "output_path": repeat("output/mrm_output.nc", count),
                "output_node_path": repeat("output/mrm_node_output.nc", count),
                "restart_input_path": repeat("restart/mrm_restart_in.nc", count),
                "restart_output_path": repeat("restart/mrm_restart_out.nc", count),
                "diagnostics_path": repeat("output/mrm_diagnostics.nc", count),
            },
            "config_time": {"time_step": repeat(1, count)},
        }
    }
    geo_values = geoparameter(dialog, geo_count)
    if geo_values is not None:
        values["parameter"] = {"geoparameter": {"GeoParam": geo_values}}
    _apply_land_cover(values["main"], settings, count)
    _apply_soil(values["main"], settings, count)
    return values


def _apply_land_cover(main, settings, count):
    land_cover = settings.get("land_cover", {})
    if not isinstance(land_cover, dict):
        return
    path = land_cover.get("output_path") or land_cover.get("path")
    if path:
        main["config_mpr"]["land_cover_path"] = repeat(str(path), count)
    variable = land_cover.get("variable")
    if variable:
        main["config_mpr"]["land_cover_var"] = repeat(str(variable), count)


def _apply_soil(main, settings, count):
    soil = settings.get("soil", {})
    horizons = soil.get("horizons", []) if isinstance(soil, dict) else []
    if not isinstance(soil, dict):
        return
    mode = int(soil.get("soil_db_mode", 1 if horizons else 0))
    path = soil.get("output_path") or soil.get("path")
    if path:
        field = "soil_horizon_class_path" if mode == 1 else "soil_class_path"
        main["config_input"][field] = repeat(str(path), count)
    variable = soil.get("variable")
    if variable:
        main["config_input"]["soil_class_var"] = repeat(str(variable), count)
    lookup = soil.get("classdefinition_path") or soil.get("lut_path")
    if lookup:
        main["config_mpr"]["soil_lut_path"] = repeat(str(lookup), count)
    if not horizons:
        main["config_mpr"]["soil_db_mode"] = repeat(mode, count)
        return
    lower_depths = []
    for horizon in horizons[:10]:
        if not isinstance(horizon, dict):
            continue
        try:
            lower_depths.append(int(round(float(horizon["lower_depth"]))))
        except (KeyError, TypeError, ValueError):
            continue
    if not lower_depths:
        return
    main["config_mpr"].update(
        {
            "soil_db_mode": repeat(mode, count),
            "tillage_depth": repeat(
                int(soil.get("tillage_depth") or lower_depths[0]), count
            ),
            "n_layers": repeat(len(lower_depths), count),
            "soil_depth": [
                repeat(depth, count) for depth in lower_depths
            ]
            + [repeat(0, count) for _ in range(10 - len(lower_depths))],
        }
    )


__all__ = ["build_dimensions", "build_initial_values"]
