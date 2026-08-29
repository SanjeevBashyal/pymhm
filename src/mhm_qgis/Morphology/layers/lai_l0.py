# -*- coding: utf-8 -*-
"""Place a LAI series on a target grid without holding the whole cube.

LAI follows the same two-stage path as every other layer. Execute All resamples
the source onto the filled DEM grid; Morphology Setup then window-copies that
staged file onto the common L0 extent and masks it. Both stages run through the
one streaming writer below, which fills the target one block of rows at a time
for one time step at a time, so peak memory does not grow with the record
length or the grid size.
"""
from __future__ import annotations

import os
from pathlib import Path


NODATA = -9999.0
# Cells the source cube does not reach are padded with a valid LAI instead of
# nodata, so widening the grid to the model extent leaves no holes for mHM.
PAD_VALUE = 0.0
DEFAULT_BLOCK_BYTES = 32 * 1024 ** 2
CHUNK_ROWS = 128
CHUNK_COLS = 1024
# Level 1 compresses an upsampled LAI grid about as well as level 4 (76:1 vs
# 78:1 on a 100x upsample) in roughly half the time, and the time dominates.
LAI_COMPRESS_LEVEL = 1
COORD_COMPRESS_LEVEL = 4


def block_row_count(ncols: int, block_bytes: int = DEFAULT_BLOCK_BYTES) -> int:
    """Return how many target rows one working block should cover."""
    per_row = max(1, int(ncols)) * 8
    rows = max(1, int(block_bytes) // per_row)
    # Keep blocks a whole number of chunk rows so no write is a partial chunk.
    rows = max(CHUNK_ROWS, (rows // CHUNK_ROWS) * CHUNK_ROWS)
    return min(rows, 4096)


def output_byte_size(steps: int, nrows: int, ncols: int) -> int:
    """Return the uncompressed float64 size of the placed LAI array."""
    return int(steps) * int(nrows) * int(ncols) * 8


def nearest_source_indices(values, axis):
    """Return nearest indices into an ascending 1-D axis, clamped at the ends."""
    import numpy as np

    axis = np.asarray(axis, dtype="float64")
    values = np.asarray(values, dtype="float64")
    if axis.size == 0:
        raise ValueError("The LAI source axis is empty.")
    if axis.size == 1:
        return np.zeros(values.shape, dtype="int64")
    right = np.clip(np.searchsorted(axis, values), 1, axis.size - 1)
    left = right - 1
    take_right = np.abs(values - axis[right]) < np.abs(values - axis[left])
    return np.where(take_right, right, left).astype("int64")


def bracketing_weights(values, axis):
    """Return lower/upper indices and the upper weight for linear interpolation.

    Values outside the source extent clamp to the edge cell, which extends the
    boundary value rather than inventing one.
    """
    import numpy as np

    axis = np.asarray(axis, dtype="float64")
    values = np.asarray(values, dtype="float64")
    if axis.size == 0:
        raise ValueError("The LAI source axis is empty.")
    if axis.size == 1:
        zeros = np.zeros(values.shape, dtype="int64")
        return zeros, zeros, np.zeros(values.shape, dtype="float64")
    upper = np.clip(np.searchsorted(axis, values), 1, axis.size - 1)
    lower = upper - 1
    span = axis[upper] - axis[lower]
    weight = np.where(span != 0, (values - axis[lower]) / np.where(span != 0, span, 1.0), 0.0)
    return lower, upper, np.clip(weight, 0.0, 1.0)


def longitude_to_source_convention(values, source_lon):
    """Shift target longitudes into the source longitude convention."""
    import numpy as np

    source_lon = np.asarray(source_lon, dtype="float64")
    if source_lon.size and np.nanmin(source_lon) >= 0 and np.nanmax(source_lon) > 180:
        return np.mod(values, 360.0)
    return values


def separable_axes(lon_block, lat_block):
    """Return 1-D lon/lat axes when the block varies along one dimension only.

    A geographic target grid has one longitude per column and one latitude per
    row, so the lookup collapses to 1-D index arrays and an ``ix_`` gather. That
    is several times faster than indexing with full 2-D index arrays.
    """
    import numpy as np

    if lon_block.ndim != 2 or lat_block.ndim != 2:
        return None
    if not np.array_equal(lon_block[0], lon_block[-1]):
        return None
    if not np.array_equal(lat_block[:, 0], lat_block[:, -1]):
        return None
    return lon_block[0], lat_block[:, 0]


def lai_window_offsets(source_x, source_y, target_header):
    """Return the integer (row, column) offset of the target grid into a source.

    ``source index = target index + offset`` on both axes. Raises when the two
    grids do not share a cell size and cell boundaries, because the copy would
    then silently shift the data.
    """
    import numpy as np

    source_x = np.asarray(source_x, dtype="float64")
    source_y = np.asarray(source_y, dtype="float64")
    if source_x.size < 2 or source_y.size < 2:
        raise ValueError("The staged LAI grid needs at least two cells per axis.")
    cellsize = float(target_header["cellsize"])
    tolerance = abs(cellsize) * 1e-6

    source_cell_x = float(abs(source_x[1] - source_x[0]))
    source_cell_y = float(abs(source_y[1] - source_y[0]))
    for name, value in (("x", source_cell_x), ("y", source_cell_y)):
        if abs(value - cellsize) > tolerance:
            raise ValueError(
                f"Staged LAI {name} cell size {value} does not match the target "
                f"cell size {cellsize}; it cannot be cropped without resampling."
            )

    # yc runs from the top down, so its first centre is half a cell below ymax.
    source_ymax = float(source_y[0]) + cellsize / 2.0
    source_xmin = float(source_x[0]) - cellsize / 2.0
    target_ymax = float(target_header["yllcorner"]) + int(
        target_header["nrows"]) * cellsize
    target_xmin = float(target_header["xllcorner"])

    row_offset = (source_ymax - target_ymax) / cellsize
    column_offset = (target_xmin - source_xmin) / cellsize
    for name, value in (("row", row_offset), ("column", column_offset)):
        if abs(value - round(value)) > 1e-6:
            raise ValueError(
                f"The staged LAI grid is not aligned to the target {name} grid."
            )
    return int(round(row_offset)), int(round(column_offset))


class ResampleSampler:
    """Sample an in-memory source cube onto target rows by lon/lat lookup."""

    def __init__(self, source_values, source_lat, source_lon,
                 method="bilinear", blank_fill=None):
        """``blank_fill`` replaces cells the source cannot supply, if given."""
        method = str(method or "bilinear").strip().lower()
        if method not in {"bilinear", "nearest"}:
            raise ValueError(f"Unsupported LAI resampling method: {method}")
        self.values = source_values
        self.lat = source_lat
        self.lon = source_lon
        self.method = method
        self.blank_fill = None if blank_fill is None else float(blank_fill)
        self._plan = None

    @property
    def steps(self) -> int:
        return int(self.values.shape[0])

    def prepare_block(self, _start, _stop, lon_block, lat_block):
        import numpy as np

        lon_block = longitude_to_source_convention(lon_block, self.lon)
        axes = separable_axes(lon_block, lat_block)
        separable = axes is not None
        if separable:
            lon_values, lat_values = axes
        else:
            lon_values, lat_values = lon_block, lat_block

        if self.method == "nearest":
            rows = nearest_source_indices(lat_values, self.lat)
            columns = nearest_source_indices(lon_values, self.lon)
            gather = np.ix_(rows, columns) if separable else (rows, columns)
            self._plan = ("nearest", gather)
            return

        row0, row1, row_weight = bracketing_weights(lat_values, self.lat)
        col0, col1, col_weight = bracketing_weights(lon_values, self.lon)
        if separable:
            corners = tuple(
                np.ix_(rows, columns)
                for rows in (row0, row1)
                for columns in (col0, col1)
            )
            row_weight = row_weight[:, None]
            col_weight = col_weight[None, :]
        else:
            corners = tuple(
                (rows, columns)
                for rows in (row0, row1)
                for columns in (col0, col1)
            )
        weights = (
            (1.0 - row_weight) * (1.0 - col_weight),
            (1.0 - row_weight) * col_weight,
            row_weight * (1.0 - col_weight),
            row_weight * col_weight,
        )
        self._plan = ("bilinear", corners, weights)

    def sample(self, step):
        import numpy as np

        values = np.asarray(self.values[step], dtype="float64")
        if self._plan[0] == "nearest":
            return self._filled(values[self._plan[1]])

        _kind, corners, weights = self._plan
        total = None
        weight_sum = None
        for gather, weight in zip(corners, weights):
            corner = values[gather]
            valid = np.isfinite(corner)
            contribution = np.where(valid, corner, 0.0) * weight
            share = np.where(valid, weight, 0.0)
            total = contribution if total is None else total + contribution
            weight_sum = share if weight_sum is None else weight_sum + share
        # Renormalise so a cell next to missing source data stays usable
        # instead of being poisoned by a single NaN corner.
        return self._filled(np.divide(
            total, weight_sum,
            out=np.full(total.shape, np.nan),
            where=weight_sum > 0,
        ))

    def _filled(self, block):
        """Replace cells the source could not supply with ``blank_fill``."""
        import numpy as np

        if self.blank_fill is None:
            return block
        return np.where(np.isfinite(block), block, self.blank_fill)


class WindowSampler:
    """Copy an aligned source variable by integer window, padding the rest.

    ``pad_value`` fills the cells the source window does not cover; it is a
    valid LAI rather than nodata so that expanding the extent never creates
    holes inside the model domain.
    """

    def __init__(self, variable, row_offset, column_offset, pad_value=PAD_VALUE):
        self.variable = variable
        self.row_offset = int(row_offset)
        self.column_offset = int(column_offset)
        self.pad_value = float(pad_value)
        self._source_rows = int(variable.shape[1])
        self._source_columns = int(variable.shape[2])
        self._window = None

    @property
    def steps(self) -> int:
        return int(self.variable.shape[0])

    def prepare_block(self, start, stop, _lon_block, _lat_block):
        rows = stop - start
        source_start = start + self.row_offset
        read_start = max(source_start, 0)
        read_stop = min(source_start + rows, self._source_rows)
        column_start = max(self.column_offset, 0)
        column_stop = min(self.column_offset + self._columns, self._source_columns)
        self._window = (
            rows,
            read_start,
            read_stop,
            read_start - source_start,
            column_start,
            column_stop,
            column_start - self.column_offset,
        )

    def bind(self, columns):
        """Record the target column count before the first block is prepared."""
        self._columns = int(columns)

    def sample(self, step):
        import numpy as np

        (rows, read_start, read_stop, target_row,
         column_start, column_stop, target_column) = self._window
        block = np.full((rows, self._columns), self.pad_value, dtype="float64")
        if read_stop > read_start and column_stop > column_start:
            values = np.asarray(
                self.variable[step, read_start:read_stop, column_start:column_stop],
                dtype="float64",
            )
            block[
                target_row:target_row + (read_stop - read_start),
                target_column:target_column + (column_stop - column_start),
            ] = values
        return block


def stream_lai_grid(
        output_path,
        *,
        coordinate_dataset,
        sampler,
        x_centers,
        y_centers,
        row_lonlat,
        lai_attrs=None,
        mask=None,
        block_bytes: int = DEFAULT_BLOCK_BYTES,
        task=None,
        log=None) -> str:
    """Write a LAI NetCDF one block of target rows and one time step at a time.

    ``coordinate_dataset`` is a small xarray dataset carrying time, time_bnds,
    yc, xc and the global attributes; it is written first so xarray handles the
    CF time encoding. ``row_lonlat(start, stop)`` returns the WGS84 longitude
    and latitude meshes for target rows ``[start, stop)``.
    """
    import numpy as np
    from netCDF4 import Dataset

    nrows = int(len(y_centers))
    ncols = int(len(x_centers))
    steps = int(sampler.steps)
    rows_per_block = block_row_count(ncols, block_bytes)
    if hasattr(sampler, "bind"):
        sampler.bind(ncols)

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.unlink(missing_ok=True)
    coordinate_dataset.to_netcdf(temporary)

    chunks = (1, min(nrows, CHUNK_ROWS), min(ncols, CHUNK_COLS))
    try:
        with Dataset(temporary, "a") as handle:
            lai = handle.createVariable(
                "lai", "f8", ("time", "yc", "xc"),
                zlib=True, complevel=LAI_COMPRESS_LEVEL,
                chunksizes=chunks, fill_value=NODATA,
            )
            lai.setncatts(dict(lai_attrs or {}))
            latitude = handle.createVariable(
                "lat", "f8", ("yc", "xc"),
                zlib=True, complevel=COORD_COMPRESS_LEVEL, chunksizes=chunks[1:],
            )
            latitude.setncatts(
                {"units": "degrees_north", "long_name": "latitude"})
            longitude = handle.createVariable(
                "lon", "f8", ("yc", "xc"),
                zlib=True, complevel=COORD_COMPRESS_LEVEL, chunksizes=chunks[1:],
            )
            longitude.setncatts(
                {"units": "degrees_east", "long_name": "longitude"})

            for start in range(0, nrows, rows_per_block):
                if task is not None and task.isCanceled():
                    raise RuntimeError("Task cancelled.")
                stop = min(start + rows_per_block, nrows)
                lon_block, lat_block = row_lonlat(start, stop)
                lon_block = np.asarray(lon_block, dtype="float64")
                lat_block = np.asarray(lat_block, dtype="float64")
                longitude[start:stop, :] = lon_block
                latitude[start:stop, :] = lat_block

                sampler.prepare_block(start, stop, lon_block, lat_block)
                block_mask = None if mask is None else np.asarray(
                    mask[start:stop, :], dtype=bool)

                for step in range(steps):
                    if task is not None and task.isCanceled():
                        raise RuntimeError("Task cancelled.")
                    placed = sampler.sample(step)
                    if block_mask is not None:
                        placed = np.where(block_mask, placed, NODATA)
                    lai[step, start:stop, :] = placed

                if task is not None:
                    task.setProgress(100.0 * stop / nrows)
                if log:
                    log(
                        f"LAI rows {stop}/{nrows} written "
                        f"({steps} time step(s))."
                    )
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise

    if output.exists():
        output.unlink()
    os.replace(temporary, output)
    return str(output)


__all__ = [
    "DEFAULT_BLOCK_BYTES",
    "NODATA",
    "PAD_VALUE",
    "ResampleSampler",
    "WindowSampler",
    "block_row_count",
    "lai_window_offsets",
    "bracketing_weights",
    "longitude_to_source_convention",
    "nearest_source_indices",
    "output_byte_size",
    "separable_axes",
    "stream_lai_grid",
]
