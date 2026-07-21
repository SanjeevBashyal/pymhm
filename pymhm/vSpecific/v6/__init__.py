"""mHM v6 namelist initial values."""

from __future__ import annotations

from typing import Any

from ..common import domain_count, geology_count, geoparameter, repeat, resolution


def build_dimensions(dialog: Any) -> dict[str, int]:
    return {
        "max_layers": 10,
        "n_domains": domain_count(dialog),
        "n_geo_units": geology_count(dialog),
    }


def build_initial_values(dialog: Any) -> dict[str, Any]:
    count = domain_count(dialog)
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
    return values


__all__ = ["build_dimensions", "build_initial_values"]
