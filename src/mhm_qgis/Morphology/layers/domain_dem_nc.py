# -*- coding: utf-8 -*-
"""Write one domain DEM as NetCDF for the v6 layout.

v6 expects `data/<domain>/dem.nc` instead of `dem.asc`. The content is the same
masked elevation on the common L0 grid, so mHM can use it as both the domain
elevation and the domain mask.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Mapping

from ..file_tasks import _aligned_l0_window, _dependencies, _rasterize_target_mask


NODATA = -9999.0
BLOCK_ROWS = 256


def write_domain_dem_netcdf(
        source_path,
        output_path,
        target_header: Mapping[str, Any],
        mask_vector,
        *,
        reference_path=None) -> Path:
    """Mask the cropped L0 DEM with a domain polygon and write `dem.nc`."""
    import numpy as np
    from netCDF4 import Dataset

    _np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    source = gdal.Open(str(source_path), gdal.GA_ReadOnly)
    if source is None:
        raise RuntimeError(f"Could not open the cropped L0 DEM: {source_path}")
    projection = source.GetProjection()
    window = _aligned_l0_window(source, source_path, target_header, reference_path)
    band = source.GetRasterBand(1)
    source_nodata = band.GetNoDataValue()
    keep = _rasterize_target_mask(mask_vector, target_header, projection)

    rows = int(target_header["nrows"])
    columns = int(target_header["ncols"])
    cellsize = float(target_header["cellsize"])
    xmin = float(target_header["xllcorner"])
    ymax = float(target_header["yllcorner"]) + rows * cellsize

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.unlink(missing_ok=True)
    chunks = (min(rows, 128), min(columns, 1024))
    try:
        with Dataset(temporary, "w", format="NETCDF4") as handle:
            handle.createDimension("x", columns)
            handle.createDimension("y", rows)
            x = handle.createVariable("x", "f8", ("x",))
            x.setncatts({
                "axis": "X",
                "standard_name": "projection_x_coordinate",
                "long_name": "x coordinate of projection",
                "units": "m",
            })
            x[:] = xmin + (np.arange(columns) + 0.5) * cellsize
            y = handle.createVariable("y", "f8", ("y",))
            y.setncatts({
                "axis": "Y",
                "standard_name": "projection_y_coordinate",
                "long_name": "y coordinate of projection",
                "units": "m",
            })
            y[:] = ymax - (np.arange(rows) + 0.5) * cellsize
            dem = handle.createVariable(
                "dem", "f8", ("y", "x"),
                zlib=True, complevel=1, chunksizes=chunks, fill_value=NODATA,
            )
            dem.setncatts({
                "long_name": "elevation",
                "standard_name": "height_above_mean_sea_level",
                "units": "m",
                "missing_value": NODATA,
            })

            blank = np.full((1, columns), NODATA, dtype="float64")
            for start in range(0, rows, BLOCK_ROWS):
                stop = min(start + BLOCK_ROWS, rows)
                block = np.repeat(blank, stop - start, axis=0)
                overlap_start = max(start, window["target_y"])
                overlap_stop = min(stop, window["target_y"] + window["copy_rows"])
                if window["copy_rows"] > 0 and overlap_stop > overlap_start:
                    read = band.ReadAsArray(
                        window["source_x"],
                        window["source_y"] + (overlap_start - window["target_y"]),
                        window["copy_cols"],
                        overlap_stop - overlap_start,
                    )
                    if read is None:
                        raise RuntimeError(
                            f"Could not read the L0 DEM: {source_path}")
                    block[
                        overlap_start - start: overlap_stop - start,
                        window["target_x"]:
                        window["target_x"] + window["copy_cols"],
                    ] = read
                if source_nodata is not None:
                    block = np.where(block == source_nodata, NODATA, block)
                dem[start:stop, :] = np.where(keep[start:stop, :], block, NODATA)
            handle.Conventions = "CF-1.8"
            handle.title = "Domain elevation prepared by mhm_qgis"
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    finally:
        band = None
        source = None
    if output.exists():
        output.unlink()
    os.replace(temporary, output)
    return output


__all__ = ["write_domain_dem_netcdf"]
