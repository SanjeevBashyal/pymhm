# -*- coding: utf-8 -*-
"""QGIS-free LAI source reading, resampling, and cropping.

This is the whole LAI pipeline expressed over paths and primitives, so it can
run inside a ``QgsTask`` worker. `Morphology/layers/lai.py` is the thin QGIS
boundary that resolves widgets and CRSs into these arguments.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import NamedTuple

from ...lai_temporal import prepare_lai_temporal
from .lai_l0 import (
    PAD_VALUE,
    ResampleSampler,
    WindowSampler,
    lai_window_offsets,
    output_byte_size,
    stream_lai_grid,
)


NODATA = -9999.0
LAI_MAX_BYTES_ENV = "MHM_QGIS_LAI_MAX_BYTES"
# Measured on a 100x upsample of real GIMMS LAI: bilinear output compresses only
# 1.80:1 on the DEM grid and 1.85:1 on the L0 grid, because every cell differs.
# Nearest-neighbour reached 70:1 on its long runs of repeated values, so do not
# reuse that figure here. The guard re-reads free space at each stage, so an
# accurate ratio lets stage 1 through and stops stage 2 once the disk is short.
LAI_ASSUMED_COMPRESSION = 1.8


class LaiTargetGrid(NamedTuple):
    """Target cell centres plus a lazy per-row-block WGS84 mesh generator."""

    x_centers: object
    y_centers: object
    row_lonlat: object
    geographic: bool


def is_geographic_crs_string(crs_string) -> bool:
    """Return True when a CRS string describes a lon/lat system."""
    text = str(crs_string or "").strip()
    if not text:
        return False
    if text.upper() in {"EPSG:4326", "OGC:CRS84", "CRS84"}:
        return True
    try:
        from pyproj import CRS

        return bool(CRS.from_user_input(text).is_geographic)
    except Exception:
        return False


def lai_grid_byte_size(steps: int, nrows: int, ncols: int) -> int:
    """Return the double-precision byte size of the placed LAI array."""
    return output_byte_size(steps, nrows, ncols)


def assert_lai_output_fits(steps: int, nrows: int, ncols: int, folder) -> int:
    """Reject an LAI request that cannot fit on disk, and return its size.

    Memory is bounded by the streaming writer, so the remaining limit is the
    output file. It is written compressed, but refuse outright when even a
    generous compression estimate cannot fit in the free space.
    """
    import shutil

    required = lai_grid_byte_size(steps, nrows, ncols)
    override = os.environ.get(LAI_MAX_BYTES_ENV, "").strip()
    if override:
        try:
            limit = int(override)
        except ValueError:
            limit = 0
        if limit > 0 and required > limit:
            raise MemoryError(
                f"LAI would write {required / 1024 ** 3:.1f} GiB: "
                f"{int(steps)} time step(s) on a {int(ncols)} x {int(nrows)} "
                f"L0 grid, over the {LAI_MAX_BYTES_ENV} limit of "
                f"{limit / 1024 ** 3:.1f} GiB."
            )
        return required

    try:
        free = shutil.disk_usage(str(folder)).free
    except OSError:
        return required
    if required / LAI_ASSUMED_COMPRESSION > free:
        raise MemoryError(
            f"LAI needs about {required / 1024 ** 3:.1f} GiB uncompressed "
            f"({int(steps)} time step(s) on a {int(ncols)} x {int(nrows)} L0 "
            f"grid) and only {free / 1024 ** 3:.1f} GiB is free on the output "
            "volume. Choose 'Long Term Mean Monthly Gridded Data' as the LAI "
            f"target timestep, use a smaller model extent, or set "
            f"{LAI_MAX_BYTES_ENV} to override this check."
        )
    return required


def lazy_target_grid(x_centers, y_centers, crs_string) -> LaiTargetGrid:
    """Return target cell centres plus a per-row-block WGS84 mesh generator.

    The full mesh is never built. On a 13320 x 6120 grid it would be two
    652 MiB arrays before any LAI value has been read.
    """
    import numpy as np

    x_centers = np.asarray(x_centers, dtype=np.float64)
    y_centers = np.asarray(y_centers, dtype=np.float64)
    geographic = is_geographic_crs_string(crs_string)

    def row_lonlat(start, stop):
        rows = y_centers[start:stop]
        if geographic:
            lon = np.broadcast_to(x_centers, (rows.size, x_centers.size))
            lat = np.repeat(rows[:, None], x_centers.size, axis=1)
            return lon, lat
        from pyproj import Transformer

        x_grid, y_grid = np.meshgrid(x_centers, rows)
        transform = Transformer.from_crs(crs_string, "EPSG:4326", always_xy=True)
        lon, lat = transform.transform(x_grid, y_grid)
        return np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)

    return LaiTargetGrid(x_centers, y_centers, row_lonlat, geographic)


def target_grid_from_header(header, crs_string) -> LaiTargetGrid:
    """Return the target grid described by an mHM-style header."""
    import numpy as np

    ncols = int(header["ncols"])
    nrows = int(header["nrows"])
    cellsize = float(header["cellsize"])
    if ncols <= 0 or nrows <= 0 or cellsize <= 0:
        raise ValueError("Grid header has invalid dimensions or cell size.")
    xmin = float(header["xllcorner"])
    ymax = float(header["yllcorner"]) + nrows * cellsize
    x_centers = xmin + (np.arange(ncols, dtype=np.float64) + 0.5) * cellsize
    y_centers = ymax - (np.arange(nrows, dtype=np.float64) + 0.5) * cellsize
    return lazy_target_grid(x_centers, y_centers, crs_string)


def target_grid_from_raster(raster_path, crs_string=None) -> tuple[LaiTargetGrid, str]:
    """Return the target grid of a raster file, read without QGIS."""
    import numpy as np
    from osgeo import gdal, osr

    dataset = gdal.Open(str(raster_path))
    if dataset is None:
        raise RuntimeError(f"Could not open the raster: {raster_path}")
    transform = dataset.GetGeoTransform()
    width, height = dataset.RasterXSize, dataset.RasterYSize
    if width <= 0 or height <= 0:
        raise ValueError(f"Raster has invalid dimensions: {raster_path}")
    if transform[2] or transform[4]:
        raise ValueError(f"Rotated grids are not supported: {raster_path}")
    projection = dataset.GetProjection()
    dataset = None

    x_centers = transform[0] + (np.arange(width, dtype=np.float64) + 0.5) * transform[1]
    y_centers = transform[3] + (np.arange(height, dtype=np.float64) + 0.5) * transform[5]

    resolved = str(crs_string or "")
    if projection:
        reference = osr.SpatialReference(wkt=projection)
        code = reference.GetAuthorityCode(None)
        authority = reference.GetAuthorityName(None)
        if authority and code:
            resolved = f"{authority}:{code}"
        elif not resolved:
            resolved = projection
    if not resolved:
        raise ValueError(f"Could not determine the CRS of {raster_path}")
    return lazy_target_grid(x_centers, y_centers, resolved), resolved


def find_coordinate(dataset, coordinate_type: str) -> str | None:
    """Find a 1D latitude or longitude coordinate in a NetCDF dataset."""
    if coordinate_type == "lat":
        candidates = ("lat", "latitude", "y")
        standard_name = "latitude"
        axis = "Y"
    else:
        candidates = ("lon", "longitude", "x")
        standard_name = "longitude"
        axis = "X"

    for name in candidates:
        if name in dataset.variables and dataset[name].ndim == 1:
            return name

    for name in dataset.variables:
        variable = dataset[name]
        if variable.ndim != 1:
            continue
        attrs = variable.attrs
        if str(attrs.get("standard_name", "")).lower() == standard_name:
            return name
        if str(attrs.get("axis", "")).upper() == axis:
            return name
    return None


def find_time_dimension(data_array, lat_dim, lon_dim):
    """Return the single non-spatial dimension, or None when it is ambiguous."""
    dimensions = [dim for dim in data_array.dims if dim not in (lat_dim, lon_dim)]
    if len(dimensions) != 1:
        return None
    return dimensions[0]


def find_variable(dataset, source_variable, lat_dim: str, lon_dim: str, log=None):
    """Find the LAI data variable to process."""
    candidates = []
    if source_variable:
        candidates.append(source_variable)
        basename = os.path.basename(source_variable)
        if basename not in candidates:
            candidates.append(basename)

    for candidate in candidates:
        if candidate in dataset.data_vars:
            data_array = dataset[candidate]
            if lat_dim in data_array.dims and lon_dim in data_array.dims:
                return data_array
            raise ValueError(
                f"Selected NetCDF variable '{candidate}' does not use the "
                "detected latitude/longitude dimensions."
            )

    for candidate in ("lai", "LAI", "leaf_area_index", "Leaf_Area_Index"):
        if candidate in dataset.data_vars:
            data_array = dataset[candidate]
            if lat_dim in data_array.dims and lon_dim in data_array.dims:
                return data_array

    for name, data_array in dataset.data_vars.items():
        if lat_dim not in data_array.dims or lon_dim not in data_array.dims:
            continue
        if find_time_dimension(data_array, lat_dim, lon_dim):
            if log:
                log(f"Using LAI variable '{name}'.")
            return data_array

    raise ValueError(
        "Could not find a LAI variable with latitude, longitude, and time "
        "dimensions."
    )


def use_coordinate_dimensions(data_array, lat_coord, lon_coord, lat_dim, lon_dim):
    """Use latitude and longitude coordinates as the interpolation dimensions."""
    swap = {}
    if lat_coord != lat_dim:
        swap[lat_dim] = lat_coord
    if lon_coord != lon_dim:
        swap[lon_dim] = lon_coord
    if swap:
        data_array = data_array.swap_dims(swap)
        lat_dim, lon_dim = lat_coord, lon_coord
    return data_array, lat_dim, lon_dim


def time_attrs_for_step(time_step: int) -> dict:
    """Return the time-axis attributes for a prepared LAI time step."""
    if time_step == 1:
        return {"long_name": "month", "units": "month"}
    return {"standard_name": "time", "axis": "T"}


def coordinate_dataset(
        x_centers,
        y_centers,
        crs_string,
        description,
        time_values,
        time_bounds,
        time_attrs):
    """Build the small coordinate-only dataset written before the LAI cube."""
    import numpy as np
    import xarray as xr

    dataset = xr.Dataset(
        data_vars={"time_bnds": (("time", "bnds"), time_bounds)},
        coords={
            "time": time_values,
            "yc": np.asarray(y_centers, dtype=np.float64),
            "xc": np.asarray(x_centers, dtype=np.float64),
            "bnds": np.arange(2, dtype=np.int8),
        },
        attrs={
            "description": description,
            "projection": str(crs_string or "").lower(),
        },
    )
    dataset["time"].attrs.update(dict(time_attrs or {}))
    dataset["time"].attrs["bounds"] = "time_bnds"
    units = "degrees" if is_geographic_crs_string(crs_string) else "m"
    dataset["yc"].attrs.update({"axis": "Y", "units": units})
    dataset["xc"].attrs.update({"axis": "X", "units": units})
    return dataset


def resample_lai_file_to_grid(
        source_path,
        output_path,
        *,
        target_grid,
        crs_string,
        source_variable=None,
        input_resolution="Monthly",
        target_timestep="Long Term Mean Monthly Gridded Data",
        description="Monthly LAI resampled to the target grid",
        method="bilinear",
        blank_fill=None,
        task=None,
        log=None) -> str:
    """Aggregate LAI in time, then stream it onto ``target_grid``.

    Only the small source cube is held; the target is filled and written one
    block of rows at a time for one time step at a time, so peak memory does
    not grow with the record length or the grid size.
    """
    import numpy as np
    import xarray as xr

    dataset = xr.open_dataset(source_path)
    try:
        lat_coord = find_coordinate(dataset, "lat")
        lon_coord = find_coordinate(dataset, "lon")
        if lat_coord is None or lon_coord is None:
            raise ValueError(
                "LAI NetCDF must contain 1D latitude and longitude coordinates.")

        promote = [
            name for name in (lat_coord, lon_coord) if name in dataset.data_vars
        ]
        if promote:
            dataset = dataset.set_coords(promote)

        lat_dim = dataset[lat_coord].dims[0]
        lon_dim = dataset[lon_coord].dims[0]
        lai_data = find_variable(dataset, source_variable, lat_dim, lon_dim, log)
        lai_data, lat_dim, lon_dim = use_coordinate_dimensions(
            lai_data, lat_coord, lon_coord, lat_dim, lon_dim)
        time_dim = find_time_dimension(lai_data, lat_dim, lon_dim)
        if time_dim is None:
            raise ValueError(
                "LAI variable must contain exactly one temporal dimension.")
        extra = [
            dim for dim in lai_data.dims
            if dim not in (time_dim, lat_dim, lon_dim)
        ]
        if extra:
            raise ValueError(
                f"LAI variable has unsupported extra dimension(s): {', '.join(extra)}.")

        temporal = prepare_lai_temporal(
            lai_data, input_resolution, target_timestep, time_dim=time_dim)
        lai_data = temporal.data.sortby(lat_coord).sortby(lon_coord)
        x_centers, y_centers = target_grid.x_centers, target_grid.y_centers
        steps = int(lai_data.sizes.get("time", 1))
        required = assert_lai_output_fits(
            steps, len(y_centers), len(x_centers),
            os.path.dirname(str(output_path)) or ".")
        if log:
            log(
                f"LAI: placing {steps} time step(s) on a "
                f"{len(x_centers)} x {len(y_centers)} grid "
                f"({required / 1024 ** 3:.1f} GiB uncompressed), streamed in "
                f"row blocks using {method} interpolation."
            )

        # The source stays small enough to hold; only the target is huge.
        source_values = np.asarray(
            lai_data.transpose("time", lat_coord, lon_coord).values)
        attrs = dict(lai_data.attrs)
        for reserved in (
                "_FillValue", "missing_value", "scale_factor",
                "add_offset", "coordinates", "grid_mapping"):
            attrs.pop(reserved, None)
        attrs.setdefault("long_name", "leaf area index")
        attrs.setdefault("units", "1")
        attrs["nodata_value"] = NODATA

        return stream_lai_grid(
            output_path,
            coordinate_dataset=coordinate_dataset(
                x_centers, y_centers, crs_string, description,
                lai_data["time"].values, temporal.time_bounds,
                time_attrs_for_step(temporal.time_step),
            ),
            sampler=ResampleSampler(
                source_values,
                np.asarray(lai_data[lat_coord].values, dtype=np.float64),
                np.asarray(lai_data[lon_coord].values, dtype=np.float64),
                method=method,
                blank_fill=blank_fill,
            ),
            x_centers=x_centers,
            y_centers=y_centers,
            row_lonlat=target_grid.row_lonlat,
            lai_attrs=attrs,
            task=task,
            log=log,
        )
    finally:
        dataset.close()


def window_copy_lai_file(
        source_path,
        output_path,
        target_header,
        crs_string,
        description,
        mask=None,
        pad_value=PAD_VALUE,
        task=None,
        log=None) -> str:
    """Copy an aligned LAI file onto ``target_header`` without resampling.

    The staged LAI already sits on the filled DEM cell grid, so the common
    extent is reached by an integer window copy, exactly like the raster
    layers. Cells beyond the staged extent take ``pad_value`` rather than
    nodata, so widening the grid leaves no holes inside the model domain. One
    time step of one row block is held at a time.
    """
    import numpy as np
    import xarray as xr
    from netCDF4 import Dataset

    target_grid = target_grid_from_header(target_header, crs_string)
    if mask is not None and tuple(np.shape(mask)) != (
            int(target_header["nrows"]), int(target_header["ncols"])):
        raise ValueError("The watershed mask does not match the L0 grid.")

    with xr.open_dataset(source_path) as staged:
        source_x = np.asarray(staged["xc"].values, dtype=np.float64)
        source_y = np.asarray(staged["yc"].values, dtype=np.float64)
        time_values = staged["time"].values
        time_bounds = np.asarray(staged["time_bnds"].values)
        time_attrs = dict(staged["time"].attrs)
        lai_attrs = dict(staged["lai"].attrs)
    row_offset, column_offset = lai_window_offsets(
        source_x, source_y, target_header)

    handle = Dataset(str(source_path), "r")
    try:
        sampler = WindowSampler(
            handle.variables["lai"],
            row_offset,
            column_offset,
            pad_value=pad_value,
        )
        # The expanded extent makes this output no smaller than the staged one,
        # so the disk guard applies here too.
        required = assert_lai_output_fits(
            sampler.steps,
            int(target_header["nrows"]),
            int(target_header["ncols"]),
            os.path.dirname(str(output_path)) or ".")
        if log:
            log(
                f"LAI: copying {sampler.steps} time step(s) onto a "
                f"{int(target_header['ncols'])} x {int(target_header['nrows'])} "
                f"grid at row offset {row_offset}, column offset "
                f"{column_offset} ({required / 1024 ** 3:.1f} GiB uncompressed)."
            )
        return stream_lai_grid(
            output_path,
            coordinate_dataset=coordinate_dataset(
                target_grid.x_centers, target_grid.y_centers, crs_string,
                description, time_values, time_bounds, time_attrs,
            ),
            sampler=sampler,
            x_centers=target_grid.x_centers,
            y_centers=target_grid.y_centers,
            row_lonlat=target_grid.row_lonlat,
            lai_attrs=lai_attrs,
            mask=mask,
            task=task,
            log=log,
        )
    finally:
        handle.close()


def run_lai_resample(options, task=None, log=None) -> str:
    """Resample LAI onto the filled DEM grid from a plain options mapping.

    This is the QgsTask entry point: ``options`` holds only paths, strings and
    numbers snapshotted on the main thread.
    """
    target_grid, crs_string = target_grid_from_raster(
        options["filled_dem"], options.get("crs_string") or None)
    return resample_lai_file_to_grid(
        options["source_path"],
        options["output_path"],
        target_grid=target_grid,
        crs_string=crs_string,
        source_variable=options.get("source_variable"),
        input_resolution=options.get("input_resolution") or "Monthly",
        target_timestep=(
            options.get("target_timestep")
            or "Long Term Mean Monthly Gridded Data"
        ),
        description="Monthly LAI resampled to the filled DEM grid",
        method=options.get("method") or "bilinear",
        blank_fill=options.get("blank_fill", 0.0),
        task=task,
        log=log,
    )


__all__ = [
    "LAI_ASSUMED_COMPRESSION",
    "LAI_MAX_BYTES_ENV",
    "LaiTargetGrid",
    "assert_lai_output_fits",
    "coordinate_dataset",
    "find_coordinate",
    "find_time_dimension",
    "find_variable",
    "is_geographic_crs_string",
    "lai_grid_byte_size",
    "lazy_target_grid",
    "resample_lai_file_to_grid",
    "run_lai_resample",
    "target_grid_from_header",
    "target_grid_from_raster",
    "time_attrs_for_step",
    "use_coordinate_dimensions",
    "window_copy_lai_file",
]
