# -*- coding: utf-8 -*-
"""Shared adapter for categorical mhm-tools raster formatters."""
from __future__ import annotations

import os
from pathlib import Path
from tempfile import TemporaryDirectory


def format_categorical_file(
    formatter,
    kind,
    input_file,
    dem_file,
    output_file,
    lookup_table,
    mapping_field,
    classdefinition_file=None,
    input_crs=None,
    dem_crs=None,
):
    """Run a CLI-shaped formatter and convert its result to GeoTIFF."""
    from mhm_tools.common.file_handler import get_raster_data, write_xarray_to_file

    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f"{kind}_", dir=output.parent) as temporary:
        temporary_path = Path(temporary)
        formatted = formatter(
            input_file=Path(input_file),
            dem_file=Path(dem_file),
            output_path=temporary_path,
            lookup_table=Path(lookup_table),
            mapping_field=mapping_field,
            output_type="nc",
            input_crs=input_crs,
            dem_crs=dem_crs,
        )
        raster = get_raster_data(formatted)
        try:
            write_xarray_to_file(raster, output, crs=raster.rio.crs)
        finally:
            raster.close()

        if classdefinition_file is not None:
            definition = temporary_path / f"{kind}_classdefinition.txt"
            destination = Path(classdefinition_file)
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(definition, destination)
    return output


__all__ = ["format_categorical_file"]
