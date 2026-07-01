# -*- coding: utf-8 -*-
"""Generated defaults for schema path fields."""

from __future__ import annotations

from typing import Any

from .namelist import canonical_name


def repeat(value: Any, count: int) -> list[Any]:
    """Return a domain-length list."""
    return [value for _ in range(max(1, int(count)))]


def configuration_path_defaults(domain_count: int) -> dict[str, dict[str, Any]]:
    """Return generated configuration defaults keyed by canonical block/name."""
    count = max(1, int(domain_count or 1))

    static = "data/static/morph"
    meteo = "data/meteo"
    lai = "data/lai"
    output = "output"
    restart = "restart"

    defaults = {
        "config_project": {
            "n_domains": count,
            "nDomains": count,
            "iFlag_cordinate_sys": 0,
        },
        "config_input": {
            "latlon_path": repeat("data/latlon.nc", count),
            "pre_path": repeat(f"{meteo}/pre/pre.nc", count),
            "pet_path": repeat(f"{meteo}/pet/pet.nc", count),
            "temp_path": repeat(f"{meteo}/tavg/tavg.nc", count),
            "tann_path": repeat(f"{meteo}/tann/tann.nc", count),
            "tmin_path": repeat(f"{meteo}/tmin/tmin.nc", count),
            "tmax_path": repeat(f"{meteo}/tmax/tmax.nc", count),
            "ssrd_path": repeat(f"{meteo}/ssrd/ssrd.nc", count),
            "strd_path": repeat(f"{meteo}/strd/strd.nc", count),
            "netrad_path": repeat(f"{meteo}/netrad/netrad.nc", count),
            "eabs_path": repeat(f"{meteo}/eabs/eabs.nc", count),
            "wind_path": repeat(f"{meteo}/windspeed/windspeed.nc", count),
            "dem_path": repeat(f"{static}/dem.asc", count),
            "slope_path": repeat(f"{static}/slope.asc", count),
            "aspect_path": repeat(f"{static}/aspect.asc", count),
            "fdir_path": repeat(f"{static}/fdir.asc", count),
            "facc_path": repeat(f"{static}/facc.asc", count),
            "geo_class_path": repeat(f"{static}/geology_class.asc", count),
            "soil_class_path": repeat(f"{static}/soil_class.asc", count),
            "lai_class_path": repeat(f"{static}/lc.asc", count),
            "morph_mask_path": repeat(f"{static}/dem.asc", count),
            "hydro_mask_path": repeat(f"{static}/dem.asc", count),
        },
        "config_mpr": {
            "land_cover_path": repeat(f"{static}/lc.asc", count),
            "lai_path": repeat(f"{lai}/lai.nc", count),
            "soil_lut_path": repeat(
                f"{static}/soil_classdefinition.txt", count),
            "geo_lut_path": repeat(
                f"{static}/geology_classdefinition.txt", count),
            "lai_lut_path": repeat(
                f"{static}/LAI_classdefinition.txt", count),
            "restart_input_path": repeat(
                f"{restart}/mpr_restart_in.nc", count),
            "restart_output_path": repeat(
                f"{restart}/mpr_restart_out.nc", count),
        },
        "config_mhm": {
            "L0Domain": list(range(1, count + 1)),
            "read_opt_domain_data": repeat(0, count),
            "output_path": repeat(f"{output}/mhm_output.nc", count),
            "mhm_file_RestartIn": repeat(
                f"{restart}/mhm_restart_in.nc", count),
            "mrm_file_RestartIn": repeat(
                f"{restart}/mrm_restart_in.nc", count),
            "mhm_file_RestartOut": repeat(
                f"{restart}/mhm_restart_out.nc", count),
            "mrm_file_RestartOut": repeat(
                f"{restart}/mrm_restart_out.nc", count),
            "restart_input_path": repeat(
                f"{restart}/mhm_restart_in.nc", count),
            "restart_output_path": repeat(
                f"{restart}/mhm_restart_out.nc", count),
        },
        "config_mrm": {
            "scc_gauges_path": repeat(f"{static}/idgauges.asc", count),
            "output_path": repeat(f"{output}/mrm_output.nc", count),
            "output_node_path": repeat(f"{output}/mrm_node_output.nc", count),
            "restart_input_path": repeat(
                f"{restart}/mrm_restart_in.nc", count),
            "restart_output_path": repeat(
                f"{restart}/mrm_restart_out.nc", count),
            "diagnostics_path": repeat(
                f"{output}/mrm_diagnostics.nc", count),
        },
        "config_time": {
            "timestep": repeat(1, count),
            "time_step": repeat(1, count),
        },
    }

    return {
        canonical_name(block): {
            canonical_name(name): value
            for name, value in values.items()
        }
        for block, values in defaults.items()
    }
