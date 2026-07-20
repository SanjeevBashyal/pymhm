# -*- coding: utf-8 -*-
"""Geology setup adapters for :mod:`mhm_tools`."""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from .._bundled import ensure_bundled_mhm_tools
from ..logging import capture_messages
from ._categorical import format_categorical_file


def rasterize_geology_map(
    input_file: str | Path,
    dem_file: str | Path,
    output_file: str | Path,
    mapping_field: str,
    *,
    lookup_table: str | Path,
    lookup_mapping_field: str,
    input_crs=None,
    dem_crs=None,
    log: Callable[[str], None] | None = None,
) -> Path:
    """Rasterize and classify a vector geology map on the DEM grid."""
    ensure_bundled_mhm_tools()
    from mhm_tools.pre import rasterize_map_data

    with capture_messages(log):
        kwargs = {
            "input_file": input_file,
            "dem_file": dem_file,
            "output_file": output_file,
            "mapping_field": mapping_field,
            "lookup_table": lookup_table,
            "lookup_mapping_field": lookup_mapping_field,
            "lookup_value_field": "GEOLOGY_CLASS",
        }
        if dem_crs is not None:
            kwargs["dem_crs"] = dem_crs
        if input_crs is not None:
            kwargs["input_crs"] = input_crs
        return Path(rasterize_map_data(**kwargs))


def format_geology_file(
    input_file: str | Path,
    dem_file: str | Path,
    output_file: str | Path,
    lookup_table: str | Path,
    mapping_field: str,
    *,
    classdefinition_file: str | Path | None = None,
    input_crs=None,
    dem_crs=None,
    log: Callable[[str], None] | None = None,
) -> Path:
    """Reclassify a raster geology map on the DEM grid."""
    ensure_bundled_mhm_tools()
    from mhm_tools.pre import format_geology_data

    with capture_messages(log):
        return format_categorical_file(
            formatter=format_geology_data,
            kind="geology",
            input_file=input_file,
            dem_file=dem_file,
            output_file=output_file,
            lookup_table=lookup_table,
            mapping_field=mapping_field,
            classdefinition_file=classdefinition_file,
            input_crs=input_crs,
            dem_crs=dem_crs,
        )


def write_geology_classdefinition_file(
    lookup_table: str | Path,
    output_file: str | Path,
    *,
    log: Callable[[str], None] | None = None,
) -> Path:
    """Write an mHM geology classdefinition through :mod:`mhm_tools`."""
    ensure_bundled_mhm_tools()
    from mhm_tools.pre import write_geology_classdefinition

    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    with capture_messages(log):
        return Path(write_geology_classdefinition(lookup_table, output))


__all__ = [
    "format_geology_file",
    "rasterize_geology_map",
    "write_geology_classdefinition_file",
]
