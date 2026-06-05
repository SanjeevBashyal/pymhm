# -*- coding: utf-8 -*-
"""Generated defaults for schema path fields."""

from __future__ import annotations

from .namelist import canonical_name


def repeat(value, count):
    """Return a domain-length list."""
    return [value for _ in range(max(1, int(count)))]


def configuration_path_defaults(domain_count):
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
        },
        "config_input": {
            "latlon_path": repeat("data/latlon.nc", count),
            "pre_path": repeat(f"{meteo}/pre/pre.nc", count),
            "temp_path": repeat(f"{meteo}/tavg/tavg.nc", count),
            "tmin_path": repeat(f"{meteo}/tmin/tmin.nc", count),
            "tmax_path": repeat(f"{meteo}/tmax/tmax.nc", count),
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
            "output_path": repeat(f"{output}/mhm_output.nc", count),
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
    }

    return {
        canonical_name(block): {
            canonical_name(name): value
            for name, value in values.items()
        }
        for block, values in defaults.items()
    }
