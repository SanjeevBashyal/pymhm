"""QGIS-free raster and watershed jobs suitable for ``QgsTask`` workers."""
from __future__ import annotations

import math
import os
import shutil
from pathlib import Path


def _dependencies():
    import numpy as np
    import pyflwdir
    from affine import Affine
    from osgeo import gdal, ogr, osr

    return np, pyflwdir, Affine, gdal, ogr, osr


def _cancelled(task):
    if task is not None and task.isCanceled():
        raise RuntimeError("Task cancelled.")


def _progress(task, value):
    if task is not None:
        task.setProgress(float(value))


def _read_raster(path, *, as_float=False):
    np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    dataset = gdal.Open(str(path))
    if dataset is None:
        raise RuntimeError(f"Could not open raster: {path}")
    band = dataset.GetRasterBand(1)
    array = band.ReadAsArray()
    if as_float:
        array = array.astype(np.float32)
    result = {
        "array": array,
        "nodata": band.GetNoDataValue(),
        "geotransform": dataset.GetGeoTransform(),
        "projection": dataset.GetProjection(),
        "rows": dataset.RasterYSize,
        "cols": dataset.RasterXSize,
    }
    band = None
    dataset = None
    return result


def _write_raster(path, array, reference, nodata, data_type):
    _np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    path = os.path.abspath(str(path))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        os.remove(path)
    dataset = gdal.GetDriverByName("GTiff").Create(
        path,
        int(reference["cols"]),
        int(reference["rows"]),
        1,
        data_type,
        options=["COMPRESS=LZW"],
    )
    if dataset is None:
        raise RuntimeError(f"Could not create raster: {path}")
    dataset.SetGeoTransform(reference["geotransform"])
    if reference.get("projection"):
        dataset.SetProjection(reference["projection"])
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(float(nodata))
    band.WriteArray(array)
    band.FlushCache()
    band = None
    dataset = None
    return path


def _same_grid(first, second):
    try:
        a = _read_raster(first)
        b = _read_raster(second)
    except RuntimeError:
        return False
    if (a["rows"], a["cols"]) != (b["rows"], b["cols"]):
        return False
    return all(
        math.isclose(float(x), float(y), rel_tol=1e-9, abs_tol=1e-9)
        for x, y in zip(a["geotransform"], b["geotransform"])
    ) and a.get("projection", "") == b.get("projection", "")


def _aligned_l0_window(source, source_path, target_header, reference_path=None):
    """Validate an L0 raster against the target header and return its copy window."""
    _np, _pfd, _affine, gdal, _ogr, osr = _dependencies()
    transform = source.GetGeoTransform()
    cellsize = float(target_header["cellsize"])
    tolerance = max(abs(cellsize), 1.0) * 1e-8
    if (
            not math.isclose(float(transform[1]), cellsize, rel_tol=0.0, abs_tol=tolerance)
            or not math.isclose(float(transform[5]), -cellsize, rel_tol=0.0, abs_tol=tolerance)
            or abs(float(transform[2])) > tolerance
            or abs(float(transform[4])) > tolerance):
        raise ValueError(
            f"L0 raster is not on the filled-DEM resolution: {source_path}"
        )

    source_xmin = float(transform[0])
    source_ymax = float(transform[3])
    source_ymin = source_ymax - source.RasterYSize * cellsize
    target_xmin = float(target_header["xllcorner"])
    target_ymin = float(target_header["yllcorner"])
    target_cols = int(target_header["ncols"])
    target_rows = int(target_header["nrows"])
    target_ymax = target_ymin + target_rows * cellsize

    def exact_offset(value, label):
        nearest = round(value)
        if not math.isclose(value, nearest, rel_tol=0.0, abs_tol=1e-7):
            raise ValueError(
                f"L0 raster {label} is not aligned to the filled-DEM cell grid: "
                f"{source_path}"
            )
        return int(nearest)

    destination_x = exact_offset(
        (source_xmin - target_xmin) / cellsize, "x origin"
    )
    destination_y = exact_offset(
        (target_ymax - source_ymax) / cellsize, "y origin"
    )

    if reference_path:
        reference = gdal.Open(str(reference_path), gdal.GA_ReadOnly)
        if reference is None:
            raise RuntimeError(f"Could not open L0 reference raster: {reference_path}")
        reference_transform = reference.GetGeoTransform()
        exact_offset((source_xmin - reference_transform[0]) / cellsize, "x origin")
        exact_offset((source_ymax - reference_transform[3]) / cellsize, "y origin")
        source_srs = osr.SpatialReference(wkt=source.GetProjection())
        reference_srs = osr.SpatialReference(wkt=reference.GetProjection())
        if source.GetProjection() and reference.GetProjection() and not source_srs.IsSame(reference_srs):
            raise ValueError(f"L0 raster CRS differs from the filled DEM: {source_path}")
        reference = None

    source_x = max(0, -destination_x)
    source_y = max(0, -destination_y)
    target_x = max(0, destination_x)
    target_y = max(0, destination_y)
    return {
        "cellsize": cellsize,
        "target_xmin": target_xmin,
        "target_ymax": target_ymax,
        "target_cols": target_cols,
        "target_rows": target_rows,
        "source_x": source_x,
        "source_y": source_y,
        "target_x": target_x,
        "target_y": target_y,
        "copy_cols": min(source.RasterXSize - source_x, target_cols - target_x),
        "copy_rows": min(source.RasterYSize - source_y, target_rows - target_y),
    }


def _copy_aligned_l0_raster(
        source_path,
        output_path,
        target_header,
        reference_path=None,
        mask=None,
        pad_value=None,
        task=None):
    """Window-copy an aligned L0 raster, padding it and optionally masking it.

    ``pad_value`` is written to the cells the source window does not reach and
    defaults to the band nodata. Categorical layers pass a valid class instead,
    so expanding a layer to the common extent never leaves nodata inside the
    model domain. Masked cells always keep the band nodata.
    """
    np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    source = gdal.Open(str(source_path), gdal.GA_ReadOnly)
    if source is None:
        raise RuntimeError(f"Could not open L0 raster: {source_path}")
    window = _aligned_l0_window(source, source_path, target_header, reference_path)
    cellsize = window["cellsize"]
    target_xmin = window["target_xmin"]
    target_ymax = window["target_ymax"]
    target_cols = window["target_cols"]
    target_rows = window["target_rows"]
    source_x = window["source_x"]
    source_y = window["source_y"]
    target_x = window["target_x"]
    target_y = window["target_y"]
    copy_cols = window["copy_cols"]
    copy_rows = window["copy_rows"]

    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.stem}.tmp{output.suffix or '.tif'}")
    temporary.unlink(missing_ok=True)
    data_type = source.GetRasterBand(1).DataType
    destination = gdal.GetDriverByName("GTiff").Create(
        str(temporary),
        target_cols,
        target_rows,
        source.RasterCount,
        data_type,
        options=["TILED=YES", "COMPRESS=LZW", "SPARSE_OK=TRUE"],
    )
    if destination is None:
        raise RuntimeError(f"Could not create cropped L0 raster: {output_path}")
    destination.SetGeoTransform(
        (target_xmin, cellsize, 0.0, target_ymax, 0.0, -cellsize)
    )
    destination.SetProjection(source.GetProjection())

    try:
        for band_index in range(1, source.RasterCount + 1):
            _cancelled(task)
            source_band = source.GetRasterBand(band_index)
            target_band = destination.GetRasterBand(band_index)
            nodata = source_band.GetNoDataValue()
            if nodata is None:
                nodata = float(target_header.get("nodata_value", -9999.0))
            target_band.SetNoDataValue(float(nodata))
            target_band.Fill(float(nodata if pad_value is None else pad_value))
            if copy_cols > 0 and copy_rows > 0:
                block_rows = min(512, copy_rows)
                for row in range(0, copy_rows, block_rows):
                    _cancelled(task)
                    height = min(block_rows, copy_rows - row)
                    values = source_band.ReadAsArray(
                        source_x, source_y + row, copy_cols, height
                    )
                    if values is None:
                        raise RuntimeError(f"Could not read L0 raster window: {source_path}")
                    if mask is not None:
                        keep = mask[
                            target_y + row: target_y + row + height,
                            target_x: target_x + copy_cols,
                        ]
                        values = np.where(keep, values, nodata).astype(values.dtype)
                    target_band.WriteArray(values, target_x, target_y + row)
            target_band.FlushCache()
            _progress(task, 100.0 * band_index / source.RasterCount)
        destination.FlushCache()
    except Exception:
        destination = None
        source = None
        temporary.unlink(missing_ok=True)
        raise
    destination = None
    source = None
    os.replace(temporary, output)
    return str(output)


def crop_aligned_l0_raster(
        source_path,
        output_path,
        target_header,
        *,
        reference_path=None,
        pad_value=None,
        task=None):
    """Window-copy an aligned L0 raster and pad outside its extent.

    Cells beyond the source extent take ``pad_value``, or the band nodata when
    no pad value is given.
    """
    return _copy_aligned_l0_raster(
        source_path,
        output_path,
        target_header,
        reference_path=reference_path,
        pad_value=pad_value,
        task=task,
    )


def _rasterize_target_mask(vector_path, target_header, projection):
    """Return a boolean mask of the target grid cells covered by a polygon layer."""
    _np, _pfd, _affine, gdal, ogr, _osr = _dependencies()
    cellsize = float(target_header["cellsize"])
    cols = int(target_header["ncols"])
    rows = int(target_header["nrows"])
    xmin = float(target_header["xllcorner"])
    ymax = float(target_header["yllcorner"]) + rows * cellsize
    dataset = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Byte)
    if dataset is None:
        raise RuntimeError("Could not allocate the L0 mask grid.")
    dataset.SetGeoTransform((xmin, cellsize, 0.0, ymax, 0.0, -cellsize))
    if projection:
        dataset.SetProjection(projection)
    source = ogr.Open(str(vector_path))
    if source is None:
        raise RuntimeError(f"Could not open the mask layer: {vector_path}")
    status = gdal.RasterizeLayer(dataset, [1], source.GetLayer(0), burn_values=[1])
    source = None
    if status != 0:
        dataset = None
        raise RuntimeError(f"Could not rasterize the mask layer: {vector_path}")
    values = dataset.GetRasterBand(1).ReadAsArray().astype(bool)
    dataset = None
    return values


def mask_aligned_l0_raster(
        source_path,
        output_path,
        target_header,
        mask_vector,
        *,
        reference_path=None,
        pad_value=None,
        task=None):
    """Apply a polygon mask to an aligned L0 raster without resampling it.

    Cells outside the polygon become nodata; cells beyond the source extent
    take ``pad_value`` as in :func:`crop_aligned_l0_raster`.
    """
    _np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    source = gdal.Open(str(source_path), gdal.GA_ReadOnly)
    if source is None:
        raise RuntimeError(f"Could not open L0 raster: {source_path}")
    projection = source.GetProjection()
    source = None
    return _copy_aligned_l0_raster(
        source_path,
        output_path,
        target_header,
        reference_path=reference_path,
        mask=_rasterize_target_mask(mask_vector, target_header, projection),
        pad_value=pad_value,
        task=task,
    )


def write_domain_dem_ascii(
        source_path,
        output_path,
        target_header,
        mask_vector,
        *,
        reference_path=None,
        nodata=-9999,
        task=None):
    """Write one domain DEM as Arc/Info ASCII on the common L0 grid.

    The source is the cropped L0 DEM, so every domain shares the model extent and
    matches the other inputs cell for cell; only the cells inside the domain
    polygon keep their value. Rows stream straight to text, so memory stays flat
    whatever the grid size.
    """
    np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    source = gdal.Open(str(source_path), gdal.GA_ReadOnly)
    if source is None:
        raise RuntimeError(f"Could not open the cropped L0 DEM: {source_path}")
    projection = source.GetProjection()
    window = _aligned_l0_window(source, source_path, target_header, reference_path)
    band = source.GetRasterBand(1)
    keep = _rasterize_target_mask(mask_vector, target_header, projection)

    cellsize = float(target_header["cellsize"])
    rows = int(target_header["nrows"])
    columns = int(target_header["ncols"])
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.unlink(missing_ok=True)

    blank = np.full(columns, float(nodata), dtype="float64")
    source_nodata = band.GetNoDataValue()
    try:
        with open(temporary, "w", encoding="utf-8") as stream:
            stream.write(
                f"ncols         {columns}\n"
                f"nrows         {rows}\n"
                f"xllcorner     {float(target_header['xllcorner'])}\n"
                f"yllcorner     {float(target_header['yllcorner'])}\n"
                f"cellsize      {cellsize}\n"
                f"NODATA_value  {nodata}\n"
            )
            for row in range(rows):
                _cancelled(task)
                values = blank
                source_row = row - window["target_y"] + window["source_y"]
                if (
                        window["copy_rows"] > 0
                        and window["target_y"] <= row
                        < window["target_y"] + window["copy_rows"]):
                    read = band.ReadAsArray(
                        window["source_x"], source_row, window["copy_cols"], 1)
                    if read is None:
                        raise RuntimeError(
                            f"Could not read the L0 DEM row {row}: {source_path}")
                    values = blank.copy()
                    values[
                        window["target_x"]:
                        window["target_x"] + window["copy_cols"]
                    ] = read[0]
                    if source_nodata is not None:
                        values = np.where(
                            values == source_nodata, float(nodata), values)
                values = np.where(keep[row], values, float(nodata))
                stream.write(" ".join(_ascii_value(value) for value in values))
                stream.write("\n")
                if task is not None and rows:
                    _progress(task, 100.0 * (row + 1) / rows)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    finally:
        band = None
        source = None
    os.replace(temporary, output)
    return str(output)


def _ascii_value(value) -> str:
    """Format one ASCII grid cell, keeping integers integral."""
    if value == int(value):
        return str(int(value))
    return repr(float(value))


def _prepared_dem_source(source, source_crs, target_crs, reprojected_path):
    _np, _pfd, _affine, gdal, _ogr, osr = _dependencies()
    source = str(source)
    if "|" in source and os.path.isfile(source.split("|", 1)[0]):
        source = source.split("|", 1)[0]
    if not target_crs or source_crs == target_crs:
        return source
    source_ref = osr.SpatialReference()
    target_ref = osr.SpatialReference()
    if source_ref.SetFromUserInput(str(source_crs or "")) != 0:
        return source
    if target_ref.SetFromUserInput(str(target_crs)) != 0 or source_ref.IsSame(target_ref):
        return source
    os.makedirs(os.path.dirname(reprojected_path), exist_ok=True)
    result = gdal.Warp(
        str(reprojected_path),
        source,
        dstSRS=str(target_crs),
        resampleAlg="near",
        multithread=True,
    )
    if result is None:
        raise RuntimeError("DEM reprojection failed.")
    result = None
    return str(reprojected_path)


def fill_dem_file(
    source,
    output_path,
    *,
    source_crs="",
    target_crs="",
    reprojected_path="",
    task=None,
):
    """Reproject when needed and fill DEM depressions using file-only inputs."""
    np, pyflwdir, _affine, gdal, _ogr, _osr = _dependencies()
    _cancelled(task)
    prepared = _prepared_dem_source(
        source, source_crs, target_crs, reprojected_path
    )
    _progress(task, 15)
    if os.path.exists(output_path) and _same_grid(output_path, prepared):
        _progress(task, 100)
        return {"filled_dem": str(output_path), "dem_source": prepared, "reused": True}
    reference = _read_raster(prepared, as_float=True)
    dem = reference["array"].copy()
    invalid = ~np.isfinite(dem)
    nodata = reference["nodata"]
    if nodata is not None and np.isfinite(nodata):
        invalid |= np.isclose(dem, nodata)
    dem[invalid] = -9999.0
    _cancelled(task)
    _progress(task, 30)
    filled, changed = pyflwdir.dem.fill_depressions(dem, nodata=-9999.0)
    filled = filled.astype(np.float32)
    filled[invalid] = -9999.0
    _cancelled(task)
    _progress(task, 85)
    path = _write_raster(output_path, filled, reference, -9999.0, gdal.GDT_Float32)
    _progress(task, 100)
    return {
        "filled_dem": path,
        "dem_source": prepared,
        "filled_cells": int(np.count_nonzero(changed)) if changed is not None else 0,
        "reused": False,
    }


#: mHM-tools derivative -> (result key, canonical geometry-folder file name).
DEM_DERIVATIVE_OUTPUTS = {
    "dem_filled": ("filled_dem", "1_dem_filled.tif"),
    "slope": ("slope", "1_dem_slope.tif"),
    "aspect": ("aspect", "1_dem_aspect.tif"),
    "facc": ("flow_accumulation", "2_flow_accumulation.tif"),
    "fdir": ("flow_direction", "2_flow_direction.tif"),
}


def dem_derivative_files(
    source,
    output_folder,
    *,
    source_crs="",
    target_crs="",
    reprojected_path="",
    compute=None,
    log=None,
    task=None,
):
    """Create the filled DEM, slope, aspect, facc and fdir in one mHM-tools pass.

    Reprojection and the same-grid reuse check stay here because they are
    project concerns; the derivatives themselves come from mHM-tools. `compute`
    lets the caller run that step out of process -- it holds every derivative in
    memory at once, which is several gigabytes on a country-scale DEM.
    """
    _cancelled(task)
    prepared = _prepared_dem_source(
        source, source_crs, target_crs, reprojected_path
    )
    _progress(task, 10)
    folder = Path(output_folder)
    targets = {
        key: folder / name for key, (_result, name) in DEM_DERIVATIVE_OUTPUTS.items()
    }
    result = {
        result_key: str(targets[key])
        for key, (result_key, _name) in DEM_DERIVATIVE_OUTPUTS.items()
    }
    result["dem_source"] = prepared
    if all(path.exists() for path in targets.values()) and _same_grid(
        str(targets["dem_filled"]), prepared
    ):
        _progress(task, 100)
        return {**result, "reused": True}

    staging = folder / "0_dem_derivatives"
    try:
        written = (compute or _compute_dem_derivatives)(prepared, staging, log)
        _cancelled(task)
        _progress(task, 85)
        for key, target in targets.items():
            target.parent.mkdir(parents=True, exist_ok=True)
            Path(written[key]).replace(target)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    _progress(task, 100)
    return {**result, "reused": False}


def _compute_dem_derivatives(prepared, staging, log):
    """Run the mHM-tools derivative pass in this process."""
    from ..mhm_tools_adapter import create_dem_derivative_files

    return create_dem_derivative_files(prepared, staging, "tif", log=log)


def _flow_context(filled_dem):
    np, pyflwdir, Affine, _gdal, _ogr, osr = _dependencies()
    reference = _read_raster(filled_dem, as_float=True)
    dem = reference["array"].copy()
    invalid = ~np.isfinite(dem)
    nodata = reference["nodata"]
    if nodata is not None and np.isfinite(nodata):
        invalid |= np.isclose(dem, nodata)
    dem[invalid] = -9999.0
    spatial_ref = osr.SpatialReference()
    geographic = bool(
        reference.get("projection")
        and spatial_ref.ImportFromWkt(reference["projection"]) == 0
        and spatial_ref.IsGeographic()
    )
    transform = Affine.from_gdal(*reference["geotransform"])
    return {
        "np": np,
        "reference": reference,
        "invalid": invalid,
        "transform": transform,
        "flwdir": pyflwdir.from_dem(
            dem, nodata=-9999.0, transform=transform, latlon=geographic
        ),
    }


def hydrology_files(
    filled_dem,
    *,
    accumulation_path="",
    area_path="",
    direction_path="",
    channel_path="",
    threshold_cells=None,
    task=None,
):
    """Create requested flow rasters and channel vector in one flow-grid pass."""
    np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    _cancelled(task)
    context = _flow_context(filled_dem)
    flwdir = context["flwdir"]
    invalid = context["invalid"]
    reference = context["reference"]
    result = {}
    _progress(task, 20)
    cells = None
    if accumulation_path or area_path or channel_path:
        cells = flwdir.upstream_area(unit="cell")
    if accumulation_path:
        values = np.rint(cells).astype(np.int32)
        values[invalid] = -9999
        result["flow_accumulation"] = _write_raster(
            accumulation_path, values, reference, -9999, gdal.GDT_Int32
        )
    if area_path:
        values = flwdir.upstream_area(unit="m2").astype(np.float32)
        values[invalid] = -9999.0
        result["flow_accumulation_area"] = _write_raster(
            area_path, values, reference, -9999.0, gdal.GDT_Float32
        )
    _cancelled(task)
    _progress(task, 55)
    if direction_path:
        values = flwdir.to_array().astype(np.uint8)
        values[invalid] = 247
        result["flow_direction"] = _write_raster(
            direction_path, values, reference, 247, gdal.GDT_Byte
        )
    if channel_path:
        result["channel_network"] = _write_channel_network(
            channel_path, cells, context, threshold_cells
        )
    _progress(task, 100)
    return result


def terrain_files(filled_dem, slope_path, aspect_path, scale=1.0, *, task=None):
    """Create slope and aspect rasters with GDAL's file API."""
    _np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    outputs = {}
    for index, (name, path, options) in enumerate(
        (
            (
                "slope",
                slope_path,
                gdal.DEMProcessingOptions(
                    format="GTiff",
                    slopeFormat="percent",
                    scale=float(scale),
                    computeEdges=False,
                    alg="Horn",
                    creationOptions=["COMPRESS=LZW"],
                ),
            ),
            (
                "aspect",
                aspect_path,
                gdal.DEMProcessingOptions(
                    format="GTiff",
                    zeroForFlat=False,
                    computeEdges=False,
                    alg="Horn",
                    creationOptions=["COMPRESS=LZW"],
                ),
            ),
        )
    ):
        _cancelled(task)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.exists(path):
            outputs[name] = str(path)
        else:
            result = gdal.DEMProcessing(str(path), str(filled_dem), name, options=options)
            if result is None:
                raise RuntimeError(f"Could not create {name} raster.")
            result = None
            outputs[name] = str(path)
        _progress(task, 50 * (index + 1))
    return outputs


def _cell_area(reference, osr):
    transform = reference["geotransform"]
    width, height = abs(float(transform[1])), abs(float(transform[5]))
    spatial_ref = osr.SpatialReference()
    if (
        reference.get("projection")
        and spatial_ref.ImportFromWkt(reference["projection"]) == 0
        and spatial_ref.IsGeographic()
    ):
        latitude = transform[3] + transform[5] * reference["rows"] / 2.0
        return max(width * 111320.0 * math.cos(math.radians(latitude)) * height * 110574.0, 1.0)
    return max(width * height, 1.0)


def _remove_vector(path, ogr):
    if os.path.exists(path):
        ogr.GetDriverByName("ESRI Shapefile").DeleteDataSource(str(path))


def _write_channel_network(path, accumulation, context, threshold_cells):
    np, _pfd, _affine, _gdal, ogr, osr = _dependencies()
    valid = accumulation[accumulation > 0]
    if valid.size == 0:
        raise RuntimeError("No valid flow-accumulation cells were found.")
    maximum = int(np.nanmax(valid))
    if threshold_cells is None:
        threshold = max(
            1,
            int(round(10_000_000.0 / _cell_area(context["reference"], osr))),
        )
        if threshold > maximum:
            threshold = max(1, int(np.nanpercentile(valid, 95)))
    else:
        threshold = int(threshold_cells)
        if threshold < 1 or threshold > maximum:
            raise ValueError(
                f"Channel threshold must be between 1 and {maximum} cells."
            )
    stream_mask = accumulation >= threshold
    order = context["flwdir"].stream_order("strahler", mask=stream_mask)
    streams = context["flwdir"].streams(
        mask=stream_mask, strord=order, uparea=accumulation
    )
    if not streams:
        raise RuntimeError("pyflwdir did not produce any stream features.")
    path = os.path.abspath(str(path))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    _remove_vector(path, ogr)
    dataset = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(path)
    spatial_ref = osr.SpatialReference()
    srs = None
    if context["reference"].get("projection"):
        spatial_ref.ImportFromWkt(context["reference"]["projection"])
        srs = spatial_ref
    layer = dataset.CreateLayer(Path(path).stem, srs, ogr.wkbLineString)
    for name, field_type in (
        ("Order", ogr.OFTInteger),
        ("UpArea", ogr.OFTReal),
        ("idx", ogr.OFTInteger),
        ("idx_ds", ogr.OFTInteger),
        ("pit", ogr.OFTInteger),
    ):
        layer.CreateField(ogr.FieldDefn(name, field_type))
    for item in streams:
        coordinates = item.get("geometry", {}).get("coordinates", [])
        if len(coordinates) < 2:
            continue
        geometry = ogr.Geometry(ogr.wkbLineString)
        for x, y in coordinates:
            geometry.AddPoint_2D(float(x), float(y))
        properties = item.get("properties", {})
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(geometry)
        feature.SetField("Order", int(properties.get("strord", 1)))
        feature.SetField("UpArea", float(properties.get("uparea", 0.0)))
        feature.SetField("idx", int(properties.get("idx", -1)))
        feature.SetField("idx_ds", int(properties.get("idx_ds", -1)))
        feature.SetField("pit", int(bool(properties.get("pit", False))))
        layer.CreateFeature(feature)
        feature = None
    dataset = None
    if not os.path.exists(path):
        raise RuntimeError("Could not write the channel network.")
    return path


def _polygonize(raster_path, vector_path):
    _np, _pfd, _affine, gdal, ogr, osr = _dependencies()
    vector_path = os.path.abspath(str(vector_path))
    os.makedirs(os.path.dirname(vector_path), exist_ok=True)
    _remove_vector(vector_path, ogr)
    raster = gdal.Open(str(raster_path))
    band = raster.GetRasterBand(1)
    spatial_ref = osr.SpatialReference()
    srs = None
    if raster.GetProjection():
        spatial_ref.ImportFromWkt(raster.GetProjection())
        srs = spatial_ref
    dataset = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(vector_path)
    layer = dataset.CreateLayer(Path(vector_path).stem, srs, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn("DN", ogr.OFTInteger))
    status = gdal.Polygonize(band, band.GetMaskBand(), layer, 0, [])
    dataset = None
    band = None
    raster = None
    if status != 0 or not os.path.exists(vector_path):
        raise RuntimeError("Could not polygonize the watershed raster.")
    return vector_path


def _delineate(context, x, y, raster_path, vector_path, basin_id):
    np, _pfd, _affine, gdal, _ogr, osr = _dependencies()
    inverse = ~context["transform"]
    col_value, row_value = inverse * (float(x), float(y))
    row, col = int(math.floor(row_value)), int(math.floor(col_value))
    reference = context["reference"]
    if row < 0 or col < 0 or row >= reference["rows"] or col >= reference["cols"]:
        raise ValueError("Outlet coordinates are outside the filled DEM grid.")
    if context["invalid"][row, col]:
        raise ValueError("Outlet coordinates are outside the valid filled DEM domain.")
    center_x, center_y = context["transform"] * (col + 0.5, row + 0.5)
    basin = context["flwdir"].basins(
        xy=([center_x], [center_y]), ids=[int(basin_id)]
    )
    mask = basin == int(basin_id)
    if not np.any(mask):
        raise RuntimeError("The delineated watershed is empty.")
    values = np.where(mask, int(basin_id), 0).astype(np.int32)
    raster_path = _write_raster(
        raster_path, values, reference, 0, gdal.GDT_Int32
    )
    saved_vector = _polygonize(raster_path, vector_path) if vector_path else None
    area = context.get("area")
    if area is None:
        area = context["flwdir"].upstream_area(unit="m2")
        context["area"] = area
    catchment = float(area[row, col])
    if not np.isfinite(catchment) or catchment <= 0:
        catchment = int(np.count_nonzero(mask)) * _cell_area(reference, osr)
    return {
        "raster_path": raster_path,
        "vector_path": saved_vector,
        "cell_center": (float(center_x), float(center_y)),
        "catchment_area_m2": catchment,
    }


def delineate_outlet_file(
    filled_dem, x, y, raster_path, vector_path, basin_id=1, *, task=None
):
    """Delineate and polygonize one outlet using no Qt or QGIS objects."""
    _cancelled(task)
    _progress(task, 10)
    result = _delineate(
        _flow_context(filled_dem), x, y, raster_path, vector_path, basin_id
    )
    _progress(task, 100)
    return result


def _dem_domain(context, basin_id, raster_path, vector_path):
    np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    values = np.where(~context["invalid"], int(basin_id), 0).astype(np.int32)
    raster = _write_raster(raster_path, values, context["reference"], 0, gdal.GDT_Int32)
    return raster, _polygonize(raster, vector_path)


def _merge_vectors(paths, output_path):
    _np, _pfd, _affine, _gdal, ogr, osr = _dependencies()
    output_path = os.path.abspath(str(output_path))
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    _remove_vector(output_path, ogr)
    inputs = [ogr.Open(str(path)) for path in paths]
    if not inputs or any(dataset is None for dataset in inputs):
        raise RuntimeError("A domain polygon could not be opened.")
    first = inputs[0].GetLayer(0)
    source_srs = first.GetSpatialRef()
    srs = source_srs.Clone() if source_srs is not None else None
    output = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(output_path)
    layer = output.CreateLayer(Path(output_path).stem, srs, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn("DN", ogr.OFTInteger))
    for dataset in inputs:
        source = dataset.GetLayer(0)
        source.ResetReading()
        for item in source:
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(item.GetGeometryRef().Clone())
            feature.SetField("DN", int(item.GetField("DN") or 0))
            layer.CreateFeature(feature)
            feature = None
    output = None
    inputs = []
    if not os.path.exists(output_path):
        raise RuntimeError("Could not merge the active domain polygons.")
    return output_path


def materialize_domain_dem_file(filled_dem, watershed_vector, output_path):
    """Mask and crop the L0 DEM to one watershed and write Arc/Info ASCII."""
    _np, _pfd, _affine, gdal, _ogr, _osr = _dependencies()
    source = gdal.Open(str(filled_dem))
    if source is None:
        raise RuntimeError(f"Could not open the filled DEM: {filled_dem}")
    transform = source.GetGeoTransform()
    if transform[2] or transform[4]:
        raise ValueError("Rotated DEM grids cannot be written as mHM ASCII files.")
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.asc")
    for path in (temporary, Path(f"{temporary}.aux.xml"), temporary.with_suffix(".prj")):
        path.unlink(missing_ok=True)
    options = gdal.WarpOptions(
        format="AAIGrid",
        cutlineDSName=str(watershed_vector),
        cropToCutline=True,
        dstNodata=-9999.0,
        xRes=abs(float(transform[1])),
        yRes=abs(float(transform[5])),
        targetAlignedPixels=True,
        resampleAlg="near",
        multithread=False,
        creationOptions=["DECIMAL_PRECISION=8"],
    )
    result = gdal.Warp(str(temporary), source, options=options)
    source = None
    if result is None:
        raise RuntimeError(f"Could not create domain DEM: {output}")
    result.FlushCache()
    result = None
    if not temporary.is_file():
        raise RuntimeError(f"Could not create domain DEM: {output}")
    os.replace(temporary, output)
    temporary_projection = temporary.with_suffix(".prj")
    if temporary_projection.is_file():
        os.replace(temporary_projection, output.with_suffix(".prj"))
    Path(f"{temporary}.aux.xml").unlink(missing_ok=True)
    return str(output)


def delineate_domains_file(
    filled_dem,
    delineations,
    *,
    dem_domain=None,
    merged_path="",
    task=None,
):
    """Delineate all selected domains and merge their polygons."""
    context = _flow_context(filled_dem)
    results = {}
    vectors = []
    total = max(1, len(delineations) + int(dem_domain is not None))
    for index, item in enumerate(delineations):
        _cancelled(task)
        outlet_id, x, y, raster, vector, domain_id, *domain_dem = item
        result = _delineate(context, x, y, raster, vector, domain_id)
        # Domain DEMs are written during Morphology Setup on the common L0 grid.
        if domain_dem:
            result["dem_path"] = str(domain_dem[0])
        results[str(outlet_id)] = result
        vectors.append(result["vector_path"])
        _progress(task, 80 * (index + 1) / total)
    if dem_domain is not None:
        domain_id, raster, vector, *domain_dem = dem_domain
        _raster, saved_vector = _dem_domain(
            context, domain_id, raster, vector
        )
        dem_result = str(domain_dem[0]) if domain_dem else None
        vectors.append(saved_vector)
    if not vectors:
        raise ValueError("Select at least one domain outlet.")
    merged = _merge_vectors(vectors, merged_path)
    _progress(task, 100)
    return {
        "outlets": results,
        "dem_domain_path": dem_result if dem_domain is not None else None,
        "merged_path": merged,
    }


__all__ = [
    "crop_aligned_l0_raster",
    "delineate_domains_file",
    "delineate_outlet_file",
    "fill_dem_file",
    "hydrology_files",
    "mask_aligned_l0_raster",
    "materialize_domain_dem_file",
    "write_domain_dem_ascii",
    "terrain_files",
]
