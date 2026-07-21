"""mHM 5.13 legacy namelist initial values."""

from __future__ import annotations

from typing import Any

from ..common import (
    coordinate_flag,
    domain_count,
    geology_count,
    geoparameter,
    repeat,
    resolution,
)


def build_dimensions(dialog: Any) -> dict[str, int]:
    return {
        "max_domains": domain_count(dialog),
        "max_geo_units": geology_count(dialog),
        "max_layers": 10,
    }


def build_initial_values(dialog: Any) -> dict[str, Any]:
    count = domain_count(dialog)
    l1 = resolution(dialog, "current_l1_resolution")
    l11 = resolution(dialog, "current_l11_resolution")
    mainconfig: dict[str, Any] = {
        "iFlag_cordinate_sys": coordinate_flag(dialog),
        "nDomains": count,
        "L0Domain": list(range(1, count + 1)),
        "read_opt_domain_data": repeat(0, count),
    }
    if l1 is not None:
        mainconfig["resolution_Hydrology"] = repeat(l1, count)

    shared: dict[str, Any] = {
        "mhm_file_RestartIn": repeat("restart/mhm_restart_in.nc", count),
        "mrm_file_RestartIn": repeat("restart/mrm_restart_in.nc", count),
        "timestep": 1,
    }
    if l11 is not None:
        shared["resolution_Routing"] = repeat(l11, count)

    values: dict[str, Any] = {
        "main": {
            "mainconfig": mainconfig,
            "mainconfig_mhm_mrm": shared,
            "directories_general": {
                "dirConfigOut": "output/",
                "dirCommonFiles": "data/static/morph/",
                "dir_Morpho": repeat("data/static/morph/", count),
                "dir_LCover": repeat("data/static/morph/", count),
                "mhm_file_RestartOut": repeat("restart/mhm_restart_out.nc", count),
                "mrm_file_RestartOut": repeat("restart/mrm_restart_out.nc", count),
                "dir_Out": repeat("output/", count),
                "file_LatLon": repeat("data/latlon.nc", count),
            },
            "directories_mHM": {
                "inputFormat_meteo_forcings": "nc",
                "dir_Precipitation": repeat("data/meteo/pre/", count),
                "dir_Temperature": repeat("data/meteo/tavg/", count),
                "dir_ReferenceET": repeat("data/meteo/pet/", count),
                "dir_MinTemperature": repeat("data/meteo/tmin/", count),
                "dir_MaxTemperature": repeat("data/meteo/tmax/", count),
                "dir_NetRadiation": repeat("data/meteo/netrad/", count),
                "dir_absVapPressure": repeat("data/meteo/eabs/", count),
                "dir_windspeed": repeat("data/meteo/windspeed/", count),
                "dir_Radiation": repeat("data/meteo/ssrd/", count),
                "time_step_model_inputs": repeat(0, count),
            },
            "directories_mRM": {
                "dir_Gauges": repeat("data/observation/streamflow/", count),
                "dir_Total_Runoff": repeat("output/", count),
                "dir_Bankfull_Runoff": repeat("data/static/morph/", count),
            },
            "processSelection": {
                "processCase": [1, 1, 1, 1, 0, 1, 1, 3, 1, 0, 0],
            },
            "LCover": {
                "nLCoverScene": 1,
                "LCoverYearStart": [1900],
                "LCoverYearEnd": [2100],
                "LCoverfName": ["lc.asc"],
            },
            "directories_MPR": {
                "dir_gridded_LAI": repeat("data/lai/", count),
            },
        }
    }
    geo_count = geology_count(dialog)
    geo_values = geoparameter(dialog, geo_count, component_major=False)
    if geo_values is not None:
        values["parameter"] = {"geoparameter": {"GeoParam": geo_values}}
    return values


__all__ = ["build_dimensions", "build_initial_values"]
