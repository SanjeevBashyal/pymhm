# -*- coding: utf-8 -*-
"""Soil setup adapters for :mod:`mhm_tools`."""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from .._bundled import ensure_bundled_mhm_tools
from ..logging import capture_messages


def rasterize_soil_map(
    input_file: str | Path,
    dem_file: str | Path,
    output_file: str | Path,
    mapping_field: str,
    *,
    lookup_table: str | Path,
    lookup_mapping_field: str,
    log: Callable[[str], None] | None = None,
) -> Path:
    """Rasterize and classify a vector soil map on the DEM grid."""
    ensure_bundled_mhm_tools()
    from mhm_tools.pre import rasterize_map_data

    with capture_messages(log):
        return Path(
            rasterize_map_data(
                input_file=input_file,
                dem_file=dem_file,
                output_file=output_file,
                mapping_field=mapping_field,
                lookup_table=lookup_table,
                lookup_mapping_field=lookup_mapping_field,
                lookup_value_field="SOIL_CLASS",
            )
        )


def format_soil_file(
    input_file: str | Path,
    output_file: str | Path,
    lookup_table: str | Path,
    mapping_field: str,
    reference_file: str | Path,
    *,
    classdefinition_file: str | Path | None = None,
    log: Callable[[str], None] | None = None,
) -> Path:
    """Reclassify a raster soil map and align it to the DEM grid."""
    ensure_bundled_mhm_tools()
    from mhm_tools.pre import format_soil_data

    output = Path(output_file)
    with capture_messages(log):
        return Path(
            format_soil_data(
                input_file=input_file,
                output_path=output.parent,
                lookup_table=lookup_table,
                mapping_field=mapping_field,
                output_file=output,
                classdefinition_file=classdefinition_file,
                reference_file=reference_file,
            )
        )


def write_soil_classdefinition_file(
    lookup_table: str | Path,
    output_file: str | Path,
    *,
    log: Callable[[str], None] | None = None,
) -> Path:
    """Write an mHM soil classdefinition through :mod:`mhm_tools`."""
    ensure_bundled_mhm_tools()
    from mhm_tools.pre import write_soil_classdefinition

    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    with capture_messages(log):
        return Path(write_soil_classdefinition(lookup_table, output))


__all__ = [
    "format_soil_file",
    "rasterize_soil_map",
    "write_soil_classdefinition_file",
]
