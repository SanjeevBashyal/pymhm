# -*- coding: utf-8 -*-
"""Path-oriented adapter for categorical :mod:`mhm_tools` formatters."""
from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path
from tempfile import TemporaryDirectory

from .._bundled import ensure_bundled_mhm_tools
from ..logging import capture_messages

_OUTPUTS = {
    "soil": ("format_soil_data", "soil_class.tif", "soil_classdefinition.txt"),
    "geology": (
        "format_geology_data",
        "geology_class.tif",
        "geology_classdefinition.txt",
    ),
    "lc": ("format_lc_data", "lc.tif", None),
}


def read_categorical_lookup_table(lookup_table):
    """Read CSV lookups with mHM-tools and delimited TXT lookups with pandas."""
    lookup_path = Path(lookup_table)
    if lookup_path.suffix.lower() != ".txt":
        ensure_bundled_mhm_tools()
        from mhm_tools.common.format_data import read_lookup_table

        return read_lookup_table(lookup_path)

    import pandas as pd

    try:
        table = pd.read_csv(
            lookup_path,
            sep=None,
            engine="python",
            encoding="utf-8-sig",
        )
    except Exception as error:
        raise ValueError(
            f"Could not read lookup table {lookup_path}: {error}"
        ) from error
    if table.empty:
        raise ValueError(f"Lookup table {lookup_path} is empty.")
    return table


def prepare_categorical_file(
    kind: str,
    input_file: str | Path,
    dem_file: str | Path,
    output_file: str | Path,
    lookup_table: str | Path,
    mapping_field: str,
    class_field: str,
    *,
    is_vector: bool = False,
    classdefinition_file: str | Path | None = None,
    input_crs=None,
    dem_crs=None,
    log: Callable[[str], None] | None = None,
) -> Path:
    """Prepare an aligned soil, geology, or land-cover GeoTIFF."""
    if kind not in _OUTPUTS:
        raise ValueError("kind must be 'soil', 'geology', or 'lc'.")

    formatter_name, raster_name, definition_name = _OUTPUTS[kind]
    if classdefinition_file is not None and definition_name is None:
        raise ValueError("Land-cover formatting does not write a classdefinition.")

    input_path = Path(input_file)
    dem_path = Path(dem_file)
    lookup_path = Path(lookup_table)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    definition_path = (
        Path(classdefinition_file) if classdefinition_file is not None else None
    )
    if definition_path is not None:
        definition_path.parent.mkdir(parents=True, exist_ok=True)

    ensure_bundled_mhm_tools()
    from mhm_tools import pre

    formatter = getattr(pre, formatter_name)
    with capture_messages(log):
        with TemporaryDirectory(
            prefix=f"pymhm_{kind}_", dir=output_path.parent
        ) as temporary:
            temporary_path = Path(temporary)
            lookup_for_tools = lookup_path
            if lookup_path.suffix.lower() == ".txt":
                lookup_for_tools = temporary_path / "lookup.csv"
                read_categorical_lookup_table(lookup_path).to_csv(
                    lookup_for_tools,
                    index=False,
                )
            formatted_input = input_path
            formatter_mapping_field = mapping_field
            formatter_input_crs = input_crs

            if is_vector:
                formatted_input = temporary_path / f"{kind}_rasterized.tif"
                rasterize_kwargs = {
                    "input_file": input_path,
                    "dem_file": dem_path,
                    "output_file": formatted_input,
                    "mapping_field": mapping_field,
                    "lookup_table": lookup_for_tools,
                    "lookup_mapping_field": mapping_field,
                    "lookup_value_field": class_field,
                }
                if input_crs is not None:
                    rasterize_kwargs["input_crs"] = input_crs
                if dem_crs is not None:
                    rasterize_kwargs["dem_crs"] = dem_crs
                pre.rasterize_map_data(**rasterize_kwargs)
                formatter_mapping_field = class_field
                formatter_input_crs = None

            formatter_kwargs = {
                "input_file": formatted_input,
                "dem_file": dem_path,
                "output_path": temporary_path,
                "lookup_table": lookup_for_tools,
                "mapping_field": formatter_mapping_field,
                "class_field": class_field,
                "output_type": "tif",
            }
            if formatter_input_crs is not None:
                formatter_kwargs["input_crs"] = formatter_input_crs
            if dem_crs is not None:
                formatter_kwargs["dem_crs"] = dem_crs
            formatter(**formatter_kwargs)

            raster_output = temporary_path / raster_name
            if not raster_output.is_file():
                raise FileNotFoundError(
                    f"mhm-tools did not create the expected raster: {raster_output}"
                )

            definition_output = (
                temporary_path / definition_name
                if definition_path is not None and definition_name is not None
                else None
            )
            if definition_output is not None and not definition_output.is_file():
                raise FileNotFoundError(
                    "mhm-tools did not create the expected classdefinition: "
                    f"{definition_output}"
                )

            os.replace(raster_output, output_path)
            if definition_output is not None and definition_path is not None:
                os.replace(definition_output, definition_path)

    return output_path


__all__ = ["prepare_categorical_file", "read_categorical_lookup_table"]
