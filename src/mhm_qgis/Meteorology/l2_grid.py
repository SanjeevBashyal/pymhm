# -*- coding: utf-8 -*-
"""Exact L2 grid placement and validation for meteorology forcing.

An mHM-ready input that already sits on the target L2 grid is placed by cell
lookup and nodata padding. Only a genuinely misaligned or reprojected input is
resampled, and then with nearest-neighbour only.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Mapping

from ..grid_resolution import (
    axes_match_header,
    read_header_file,
    reindex_to_header,
)


def slice_and_pad(da, x_dim: str, y_dim: str, header: Mapping[str, Any]):
    """Place an already aligned source on the exact L2 grid, padding with nodata."""
    return reindex_to_header(da, x_dim, y_dim, dict(header))


def assert_matches_header(dataset, variable: str, header: Mapping[str, Any]) -> None:
    """Raise when a prepared L2 dataset does not match the saved L2 header."""
    rows = int(header["nrows"])
    cols = int(header["ncols"])
    sizes = dataset[variable].sizes
    if sizes.get("lat") != rows or sizes.get("lon") != cols:
        raise ValueError(
            f"{variable}: prepared L2 grid is "
            f"{sizes.get('lat')}x{sizes.get('lon')}, expected {rows}x{cols}."
        )


def assert_header_file_matches(path: Path | str, header: Mapping[str, Any]) -> None:
    """Raise when a written L2 header file differs from the saved L2 header."""
    written = read_header_file(path)
    if written is None:
        raise ValueError(f"Could not read the written L2 header: {path}")
    if (
            int(written["ncols"]) != int(header["ncols"])
            or int(written["nrows"]) != int(header["nrows"])):
        raise ValueError(f"Written L2 header has the wrong dimensions: {path}")
    tolerance = max(abs(float(header["cellsize"])), 1.0) * 1e-9
    for key in ("xllcorner", "yllcorner", "cellsize"):
        if not math.isclose(
                float(written[key]), float(header[key]),
                rel_tol=0.0, abs_tol=tolerance):
            raise ValueError(f"Written L2 header {key} does not match: {path}")


def assert_netcdf_matches_header(path, variable, header) -> None:
    """Raise when a written L2 NetCDF does not sit on the saved L2 header."""
    import xarray as xr

    with xr.open_dataset(path) as dataset:
        if variable not in dataset:
            raise ValueError(f"{path} has no {variable} variable.")
        assert_matches_header(dataset, variable, header)


__all__ = [
    "assert_header_file_matches",
    "assert_matches_header",
    "assert_netcdf_matches_header",
    "axes_match_header",
    "slice_and_pad",
]
