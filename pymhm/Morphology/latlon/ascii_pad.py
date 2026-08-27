# -*- coding: utf-8 -*-
"""Pad an Arc/Info ASCII grid onto a larger aligned header by streaming text.

An ASCII grid is row-major
text, so padding is a line-by-line rewrite with constant memory instead.
"""
from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any, Mapping


HEADER_KEYS = ("ncols", "nrows", "xllcorner", "yllcorner", "cellsize")


def read_ascii_header(path) -> dict:
    """Return the six header values of an Arc/Info ASCII grid."""
    header: dict[str, Any] = {}
    with open(path, "r", encoding="utf-8") as stream:
        for _ in range(8):
            position = stream.tell()
            line = stream.readline()
            if not line:
                break
            parts = line.split()
            if len(parts) != 2:
                stream.seek(position)
                break
            key = parts[0].lower()
            if key in ("ncols", "nrows"):
                header[key] = int(float(parts[1]))
            elif key in (
                    "xllcorner", "yllcorner", "xllcenter", "yllcenter",
                    "cellsize", "nodata_value"):
                header[key] = float(parts[1])
            else:
                stream.seek(position)
                break
    cellsize = header.get("cellsize")
    if "xllcenter" in header and cellsize:
        header["xllcorner"] = header["xllcenter"] - 0.5 * cellsize
    if "yllcenter" in header and cellsize:
        header["yllcorner"] = header["yllcenter"] - 0.5 * cellsize
    if not all(key in header for key in HEADER_KEYS):
        raise ValueError(f"Not a readable Arc/Info ASCII header: {path}")
    return header


def ascii_window_offsets(source: Mapping[str, Any],
                         target: Mapping[str, Any]) -> tuple[int, int]:
    """Return the (row, column) offset of the source inside the target grid.

    ``target index = source index + offset``. Raises when the grids do not share
    a cell size and cell boundaries, because padding would then shift the data.
    """
    cellsize = float(target["cellsize"])
    tolerance = abs(cellsize) * 1e-6
    if abs(float(source["cellsize"]) - cellsize) > tolerance:
        raise ValueError(
            f"Source cell size {source['cellsize']} does not match the target "
            f"cell size {cellsize}; it cannot be padded without resampling."
        )
    column = (float(source["xllcorner"]) - float(target["xllcorner"])) / cellsize
    source_top = float(source["yllcorner"]) + int(source["nrows"]) * cellsize
    target_top = float(target["yllcorner"]) + int(target["nrows"]) * cellsize
    row = (target_top - source_top) / cellsize
    for name, value in (("row", row), ("column", column)):
        if abs(value - round(value)) > 1e-6:
            raise ValueError(
                f"The source grid is not aligned to the target {name} grid.")
    return int(round(row)), int(round(column))


def _format_header(target: Mapping[str, Any], nodata) -> str:
    return (
        f"ncols         {int(target['ncols'])}\n"
        f"nrows         {int(target['nrows'])}\n"
        f"xllcorner     {float(target['xllcorner'])}\n"
        f"yllcorner     {float(target['yllcorner'])}\n"
        f"cellsize      {float(target['cellsize'])}\n"
        f"NODATA_value  {nodata}\n"
    )


def pad_ascii_grid(path, target: Mapping[str, Any], nodata=-9999, pad=None) -> Path:
    """Pad an ASCII grid onto ``target`` in place, streaming row by row.

    Returns the path unchanged when the grid already matches the target. Cells
    outside the source extent become ``pad``, or ``nodata`` when no pad value
    is given. ``nodata`` stays the declared NODATA_value either way, so a class
    layer can be padded with a valid class without redefining its nodata.
    """
    path = Path(path)
    source = read_ascii_header(path)
    tolerance = abs(float(target["cellsize"])) * 1e-6
    if (
            int(source["ncols"]) == int(target["ncols"])
            and int(source["nrows"]) == int(target["nrows"])
            and abs(float(source["xllcorner"]) - float(target["xllcorner"])) <= tolerance
            and abs(float(source["yllcorner"]) - float(target["yllcorner"])) <= tolerance
            and abs(float(source["cellsize"]) - float(target["cellsize"])) <= tolerance):
        return path

    row_offset, column_offset = ascii_window_offsets(source, target)
    source_rows, source_columns = int(source["nrows"]), int(source["ncols"])
    target_rows, target_columns = int(target["nrows"]), int(target["ncols"])
    if row_offset < 0 or column_offset < 0 or (
            row_offset + source_rows > target_rows
            or column_offset + source_columns > target_columns):
        raise ValueError(
            "The target grid does not contain the source grid; padding only "
            "expands, it never crops."
        )

    text = str(nodata if pad is None else pad)
    left = f"{text} " * column_offset
    right = f" {text}" * (target_columns - column_offset - source_columns)
    blank = " ".join([text] * target_columns)
    temporary = path.with_name(f".{path.name}.pad.tmp")
    temporary.unlink(missing_ok=True)
    try:
        with open(path, "r", encoding="utf-8") as source_stream, \
                open(temporary, "w", encoding="utf-8") as output_stream:
            for _ in range(6):
                source_stream.readline()
            output_stream.write(_format_header(target, nodata))
            for _ in range(row_offset):
                output_stream.write(blank + "\n")
            written = 0
            for line in source_stream:
                values = line.strip()
                if not values:
                    continue
                output_stream.write(f"{left}{values}{right}\n")
                written += 1
            if written != source_rows:
                raise ValueError(
                    f"{path.name} has {written} data row(s) but its header "
                    f"declares {source_rows}."
                )
            for _ in range(target_rows - row_offset - source_rows):
                output_stream.write(blank + "\n")
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    os.replace(temporary, path)
    return path


__all__ = [
    "ascii_window_offsets",
    "pad_ascii_grid",
    "read_ascii_header",
]
