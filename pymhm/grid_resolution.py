# -*- coding: utf-8 -*-
"""Grid-resolution helpers for morphology, meteorology, and lat/lon files."""
from __future__ import annotations

import json
import math
import os
from decimal import Decimal, ROUND_CEILING
from pathlib import Path

from .project_layout import geometry_folder, meteo_folder


METEO_GRID_METADATA = "meteo_grid_metadata.json"
NODATA_VALUE = -9999.0
PROJECTED_CELLSIZE_PRECISION = 8
GEOGRAPHIC_CELLSIZE_PRECISION = 8


def is_geographic_unit(unit: str | None) -> bool:
    """Return True when a unit label describes geographic degrees."""
    text = (unit or "").strip().lower()
    return text in {"deg", "degree", "degrees"}


def cellsize_precision_for_unit(unit: str | None) -> int:
    """Return internal cellsize precision for projected or geographic grids."""
    if is_geographic_unit(unit):
        return GEOGRAPHIC_CELLSIZE_PRECISION
    return PROJECTED_CELLSIZE_PRECISION


def ceil_to_precision(value, precision: int) -> float:
    """Ceil a numeric value to a fixed number of decimal places."""
    numeric = float(value)
    if not math.isfinite(numeric):
        return numeric
    quantum = Decimal("1").scaleb(-int(precision))
    return float(
        Decimal(str(numeric)).quantize(
            quantum,
            rounding=ROUND_CEILING,
        )
    )


def ceil_cellsize(value, unit: str | None = None,
                  precision: int | None = None) -> float:
    """Ceil a grid cellsize using the configured internal CRS precision."""
    if precision is None:
        precision = cellsize_precision_for_unit(unit)
    return ceil_to_precision(value, precision)


def display_precision_for_unit(unit: str | None) -> int:
    """Return display precision for map/grid values in the given unit."""
    if is_geographic_unit(unit):
        return 6
    return 3


def format_resolution(value, unit: str | None = None, precision: int | None = None):
    """Return a compact display string for a grid resolution."""
    if precision is None:
        precision = display_precision_for_unit(unit)
    try:
        text = f"{float(value):.{precision}f}".rstrip("0").rstrip(".")
    except (TypeError, ValueError):
        return ""
    return text or "0"


def qgis_crs_unit_label(crs) -> str:
    """Return a short map-unit label for a QGIS CRS."""
    if crs is None or not crs.isValid():
        return ""
    if crs.isGeographic():
        return "deg"
    try:
        from qgis.core import QgsUnitTypes

        unit_label = QgsUnitTypes.toAbbreviatedString(crs.mapUnits())
        return unit_label or "m"
    except Exception:
        return "m"


def qgis_crs_to_authid(crs) -> str:
    """Return the best CRS authid available for storage and pyproj/QGIS reuse."""
    if crs is None or not crs.isValid():
        return ""
    if crs.authid():
        return crs.authid()
    if crs.postgisSrid():
        return f"EPSG:{crs.postgisSrid()}"
    return ""


def raster_resolution_info(layer) -> dict | None:
    """Return scalar resolution information from a raster layer."""
    if layer is None or not layer.isValid():
        return None

    extent = layer.extent()
    width = int(layer.width())
    height = int(layer.height())
    if width <= 0 or height <= 0:
        return None

    x_resolution = abs((extent.xMaximum() - extent.xMinimum()) / width)
    y_resolution = abs((extent.yMaximum() - extent.yMinimum()) / height)
    crs = layer.crs()
    unit = qgis_crs_unit_label(crs)
    x_resolution = ceil_cellsize(x_resolution, unit)
    y_resolution = ceil_cellsize(y_resolution, unit)
    resolution = ceil_cellsize((x_resolution + y_resolution) / 2.0, unit)
    return {
        "resolution": resolution,
        "x_resolution": x_resolution,
        "y_resolution": y_resolution,
        "unit": unit,
        "crs_authid": qgis_crs_to_authid(crs),
        "cellsize_precision": cellsize_precision_for_unit(unit),
        "header": {
            "ncols": width,
            "nrows": height,
            "xllcorner": float(extent.xMinimum()),
            "yllcorner": float(extent.yMinimum()),
            "cellsize": float(resolution),
            "nodata_value": NODATA_VALUE,
            "unit": unit,
            "cellsize_precision": cellsize_precision_for_unit(unit),
        },
    }


def read_header_file(path, unit: str | None = None) -> dict | None:
    """Read an mHM/ESRI-style header file into normalized lowercase keys."""
    path = Path(path)
    if not path.exists():
        return None

    header = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            key = parts[0].lower()
            value = parts[1]
            if key in ("ncols", "nrows"):
                header[key] = int(float(value))
            elif key in ("xllcorner", "yllcorner", "xllcenter", "yllcenter", "cellsize", "nodata_value"):
                header[key] = float(value)

    if "xllcenter" in header:
        header["xllcorner"] = header["xllcenter"] - 0.5 * header.get("cellsize", 1.0)
    if "yllcenter" in header:
        header["yllcorner"] = header["yllcenter"] - 0.5 * header.get("cellsize", 1.0)
    header.setdefault("nodata_value", NODATA_VALUE)
    if unit is not None and "cellsize" in header:
        header["cellsize"] = ceil_cellsize(header["cellsize"], unit)
        header["unit"] = unit
        header["cellsize_precision"] = cellsize_precision_for_unit(unit)

    required = {"ncols", "nrows", "xllcorner", "yllcorner", "cellsize"}
    if not required.issubset(header):
        return None
    return header


def write_header_file(path, header: dict) -> None:
    """Write a normalized grid header."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "ncols         {ncols}\n"
        "nrows         {nrows}\n"
        "xllcorner     {xllcorner}\n"
        "yllcorner     {yllcorner}\n"
        "cellsize      {cellsize}\n"
        "NODATA_value  {nodata_value}\n".format(
            ncols=int(header["ncols"]),
            nrows=int(header["nrows"]),
            xllcorner=header["xllcorner"],
            yllcorner=header["yllcorner"],
            cellsize=header["cellsize"],
            nodata_value=header.get("nodata_value", NODATA_VALUE),
        ),
        encoding="utf-8",
    )


def metadata_path(project_folder) -> Path:
    """Return the project-local meteo grid metadata path."""
    return Path(meteo_folder(project_folder)) / METEO_GRID_METADATA


def save_meteo_grid_metadata(project_folder, metadata: dict) -> None:
    """Persist selected L2 grid metadata for UI restore and latlon creation."""
    path = metadata_path(project_folder)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)


def load_meteo_grid_metadata(project_folder) -> dict | None:
    """Load persisted L2 grid metadata."""
    path = metadata_path(project_folder)
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        if isinstance(data, dict):
            return data
    except Exception:
        return None
    return None


def nearest_integer_multiple(
        raw_resolution: float,
        base_resolution: float,
        unit: str | None = None) -> tuple[float, int]:
    """Return the nearest resolution that is an integer multiple of the base."""
    if raw_resolution <= 0 or base_resolution <= 0:
        raise ValueError("Grid resolutions must be positive.")
    base_resolution = ceil_cellsize(base_resolution, unit)
    ratio = max(1, int(round(raw_resolution / base_resolution)))
    return ceil_cellsize(ratio * base_resolution, unit), ratio


def resolution_is_multiple(
        coarse_resolution: float,
        fine_resolution: float,
        tolerance=None,
        unit: str | None = None) -> bool:
    """Return True when coarse is an integer multiple of fine."""
    if coarse_resolution <= 0 or fine_resolution <= 0:
        return False
    coarse_resolution = ceil_cellsize(coarse_resolution, unit)
    fine_resolution = ceil_cellsize(fine_resolution, unit)
    if tolerance is None:
        tolerance = max(10 ** -cellsize_precision_for_unit(unit), 1e-9)
    ratio = coarse_resolution / fine_resolution
    return math.isclose(ratio, round(ratio), abs_tol=tolerance, rel_tol=0.0)


def possible_resolutions(
        fine_resolution: float,
        coarse_resolution: float,
        unit: str | None = None) -> list[float]:
    """Return all resolutions between fine and coarse compatible with both."""
    fine_resolution = ceil_cellsize(fine_resolution, unit)
    coarse_resolution = ceil_cellsize(coarse_resolution, unit)
    if not resolution_is_multiple(coarse_resolution, fine_resolution, unit=unit):
        return []
    ratio = int(round(coarse_resolution / fine_resolution))
    values = []
    for divisor in range(1, ratio + 1):
        if ratio % divisor == 0:
            quotient = ratio // divisor
            if quotient == 1:
                values.append(float(coarse_resolution))
            else:
                values.append(ceil_cellsize(coarse_resolution / quotient, unit))
    return values


def header_bounds(header: dict) -> tuple[float, float, float, float]:
    """Return xmin, xmax, ymin, ymax from a grid header."""
    xmin = float(header["xllcorner"])
    ymin = float(header["yllcorner"])
    xmax = xmin + int(header["ncols"]) * float(header["cellsize"])
    ymax = ymin + int(header["nrows"]) * float(header["cellsize"])
    return xmin, xmax, ymin, ymax


def header_for_bounds(bounds, cellsize: float, unit: str | None = None) -> dict:
    """Create a header snapped outward to cellsize boundaries."""
    cellsize = ceil_cellsize(cellsize, unit)
    xmin, xmax, ymin, ymax = bounds
    xll = math.floor(min(xmin, xmax) / cellsize) * cellsize
    yll = math.floor(min(ymin, ymax) / cellsize) * cellsize
    xur = math.ceil(max(xmin, xmax) / cellsize) * cellsize
    yur = math.ceil(max(ymin, ymax) / cellsize) * cellsize
    ncols = max(1, int(round((xur - xll) / cellsize)))
    nrows = max(1, int(round((yur - yll) / cellsize)))
    return {
        "ncols": ncols,
        "nrows": nrows,
        "xllcorner": float(xll),
        "yllcorner": float(yll),
        "cellsize": float(cellsize),
        "nodata_value": NODATA_VALUE,
        "unit": unit or "",
        "cellsize_precision": cellsize_precision_for_unit(unit),
    }


def header_for_existing_bounds(
        reference_header: dict,
        cellsize: float,
        unit: str | None = None) -> dict:
    """Create a compatible header on the same lower-left and upper-right bounds."""
    cellsize = ceil_cellsize(cellsize, unit)
    xmin, xmax, ymin, ymax = header_bounds(reference_header)
    ncols = max(1, int(round((xmax - xmin) / cellsize)))
    nrows = max(1, int(round((ymax - ymin) / cellsize)))
    return {
        "ncols": ncols,
        "nrows": nrows,
        "xllcorner": float(xmin),
        "yllcorner": float(ymin),
        "cellsize": float(cellsize),
        "nodata_value": NODATA_VALUE,
        "unit": unit or reference_header.get("unit", ""),
        "cellsize_precision": cellsize_precision_for_unit(
            unit or reference_header.get("unit", "")),
    }


def selected_dem_layer(dialog):
    """Return the currently selected DEM layer, if available."""
    combo = getattr(dialog, "mMapLayerComboBox_dem", None)
    if combo is None:
        return None
    return combo.currentLayer()


def selected_dem_resolution_info(dialog) -> dict | None:
    """Return L0 resolution info from the currently selected DEM."""
    return raster_resolution_info(selected_dem_layer(dialog))


def project_crs_from_dialog(dialog):
    """Return the processing CRS, falling back to selected DEM CRS."""
    crs = None
    if hasattr(dialog, "get_crs"):
        crs = dialog.get_crs()
    if crs is not None and crs.isValid():
        return crs

    layer = selected_dem_layer(dialog)
    if layer is not None:
        crs = layer.crs()
        if crs.isValid():
            return crs
    return crs


def watershed_extent_in_crs(dialog, target_crs):
    """Return the merged watershed extent in target CRS, falling back to DEM extent."""
    from qgis.core import QgsCoordinateTransform, QgsProject, QgsRasterLayer, QgsVectorLayer

    candidates = []
    if getattr(dialog, "project_folder", None):
        geom = Path(geometry_folder(dialog.project_folder))
        candidates.extend([
            geom / "Watersheds" / "4_watershed_merged_vector.shp",
            geom / "4_watershed_merged_vector.shp",
        ])

    for path in candidates:
        if not path.exists():
            continue
        layer = QgsVectorLayer(str(path), "Watershed_Merged", "ogr")
        if layer.isValid():
            extent = _transform_extent(layer.extent(), layer.crs(), target_crs)
            return extent, str(path)

    layer = selected_dem_layer(dialog)
    if layer is not None and layer.isValid():
        extent = _transform_extent(layer.extent(), layer.crs(), target_crs)
        return extent, layer.source()

    if getattr(dialog, "project_folder", None):
        masked_dem = Path(geometry_folder(dialog.project_folder)) / "1_dem_filled_masked.tif"
        if masked_dem.exists():
            layer = QgsRasterLayer(str(masked_dem), "DEM_Masked")
            if layer.isValid():
                extent = _transform_extent(layer.extent(), layer.crs(), target_crs)
                return extent, str(masked_dem)

    return None, ""


def _transform_extent(extent, source_crs, target_crs):
    """Transform an extent between QGIS CRSs when needed."""
    from qgis.core import QgsCoordinateTransform, QgsProject

    if source_crs is None or target_crs is None:
        return extent
    if not source_crs.isValid() or not target_crs.isValid():
        return extent
    if source_crs.authid() == target_crs.authid():
        return extent
    transform = QgsCoordinateTransform(source_crs, target_crs, QgsProject.instance())
    transform.setBallparkTransformsAreAppropriate(True)
    return transform.transformBoundingBox(extent)


def extent_to_bounds(extent) -> tuple[float, float, float, float]:
    """Return xmin, xmax, ymin, ymax from a QGIS extent."""
    return (
        float(extent.xMinimum()),
        float(extent.xMaximum()),
        float(extent.yMinimum()),
        float(extent.yMaximum()),
    )


def extent_center_wgs84(extent, source_crs) -> tuple[float, float]:
    """Return lon/lat center of an extent."""
    from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsPointXY, QgsProject

    center = QgsPointXY(
        (extent.xMinimum() + extent.xMaximum()) / 2.0,
        (extent.yMinimum() + extent.yMaximum()) / 2.0,
    )
    if source_crs is not None and source_crs.isValid() and not source_crs.isGeographic():
        wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
        transform = QgsCoordinateTransform(source_crs, wgs84, QgsProject.instance())
        transform.setBallparkTransformsAreAppropriate(True)
        center = transform.transform(center)
    return float(center.x()), float(center.y())


def raw_meteo_resolution(nc_folder, target_crs, center_lonlat) -> dict:
    """Estimate raw meteorology grid spacing in the target CRS units."""
    from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsPointXY, QgsProject
    import numpy as np
    import xarray as xr

    sample = next(iter(sorted(Path(nc_folder).glob("ERA5_Land_*.nc"))), None)
    if sample is None:
        raise FileNotFoundError(f"No ERA5_Land_*.nc files found in {nc_folder}")

    dataset = None
    last_error = None
    for engine in ("netcdf4", "h5netcdf", "scipy", None):
        try:
            kwargs = {} if engine is None else {"engine": engine}
            dataset = xr.open_dataset(sample, **kwargs)
            break
        except Exception as e:
            last_error = e
    if dataset is None:
        raise RuntimeError(f"Could not open meteo NetCDF {sample}: {last_error}")

    with dataset as ds:
        lat_name = "lat" if "lat" in ds.coords else "latitude"
        lon_name = "lon" if "lon" in ds.coords else "longitude"
        if lat_name not in ds.coords or lon_name not in ds.coords:
            raise KeyError("Could not find latitude/longitude coordinates in meteo NetCDF.")

        lat_values = np.asarray(ds[lat_name].values, dtype="float64")
        lon_values = np.asarray(ds[lon_name].values, dtype="float64")

    lat_diffs = np.abs(np.diff(np.sort(np.unique(lat_values))))
    lon_diffs = np.abs(np.diff(np.sort(np.unique(lon_values))))
    if lat_diffs.size == 0 or lon_diffs.size == 0:
        raise ValueError("Meteo grid must have at least two latitude and longitude coordinates.")

    lat_step = float(np.nanmedian(lat_diffs))
    lon_step = float(np.nanmedian(lon_diffs))
    degree_resolution = (lat_step + lon_step) / 2.0

    if target_crs is None or not target_crs.isValid() or target_crs.isGeographic():
        return {
            "resolution": degree_resolution,
            "x_resolution": lon_step,
            "y_resolution": lat_step,
            "unit": "deg",
            "source_file": str(sample),
        }

    lon_center, lat_center = center_lonlat
    wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
    transform = QgsCoordinateTransform(wgs84, target_crs, QgsProject.instance())
    transform.setBallparkTransformsAreAppropriate(True)

    p0 = transform.transform(QgsPointXY(lon_center, lat_center))
    px = transform.transform(QgsPointXY(lon_center + lon_step, lat_center))
    py = transform.transform(QgsPointXY(lon_center, lat_center + lat_step))
    dx = math.hypot(px.x() - p0.x(), px.y() - p0.y())
    dy = math.hypot(py.x() - p0.x(), py.y() - p0.y())
    return {
        "resolution": (dx + dy) / 2.0,
        "x_resolution": dx,
        "y_resolution": dy,
        "degree_resolution": degree_resolution,
        "degree_x_resolution": lon_step,
        "degree_y_resolution": lat_step,
        "unit": qgis_crs_unit_label(target_crs),
        "source_file": str(sample),
    }


def header_center_coordinates(header: dict) -> tuple[list[float], list[float]]:
    """Return x center coordinates and descending y center coordinates."""
    cellsize = float(header["cellsize"])
    xll = float(header["xllcorner"])
    yll = float(header["yllcorner"])
    ncols = int(header["ncols"])
    nrows = int(header["nrows"])
    x_values = [xll + (index + 0.5) * cellsize for index in range(ncols)]
    y_values = [yll + (index + 0.5) * cellsize for index in range(nrows)]
    y_values.reverse()
    return x_values, y_values


def target_lon_lat_from_header(header: dict, source_crs) -> tuple[list[float], list[float]]:
    """Return 1D lon/lat target axes for a projected grid header."""
    from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsPointXY, QgsProject

    x_values, y_values = header_center_coordinates(header)
    if source_crs is None or not source_crs.isValid() or source_crs.isGeographic():
        return x_values, y_values

    xmin, xmax, ymin, ymax = header_bounds(header)
    mid_x = (xmin + xmax) / 2.0
    mid_y = (ymin + ymax) / 2.0
    wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
    transform = QgsCoordinateTransform(source_crs, wgs84, QgsProject.instance())
    transform.setBallparkTransformsAreAppropriate(True)

    lon_values = [
        float(transform.transform(QgsPointXY(x_value, mid_y)).x())
        for x_value in x_values
    ]
    lat_values = [
        float(transform.transform(QgsPointXY(mid_x, y_value)).y())
        for y_value in y_values
    ]
    return lon_values, lat_values


def build_meteo_l2_grid(dialog, nc_folder) -> dict:
    """Build adjusted L2 grid metadata and target lat/lon axes for meteo forcing."""
    l0_info = getattr(dialog, "_grid_l0_info", None)
    if not l0_info and hasattr(dialog, "filled_dem_resolution_info"):
        l0_info = dialog.filled_dem_resolution_info()
    if not l0_info:
        raise ValueError("Filled DEM is required before preparing meteorology data.")

    target_crs = project_crs_from_dialog(dialog)
    if target_crs is None or not target_crs.isValid():
        raise ValueError("Please set a valid processing CRS before preparing meteorology data.")

    extent, extent_source = watershed_extent_in_crs(dialog, target_crs)
    if extent is None:
        raise ValueError("Could not determine watershed or DEM extent for meteo clipping.")

    center_lonlat = extent_center_wgs84(extent, target_crs)
    raw = raw_meteo_resolution(nc_folder, target_crs, center_lonlat)
    target_unit = qgis_crs_unit_label(target_crs)
    l0_resolution = ceil_cellsize(l0_info["resolution"], target_unit)
    adjusted_resolution, ratio = nearest_integer_multiple(
        raw["resolution"],
        l0_resolution,
        target_unit,
    )
    l2_header = header_for_bounds(extent_to_bounds(extent), adjusted_resolution, target_unit)
    lon_values, lat_values = target_lon_lat_from_header(l2_header, target_crs)
    wgs84_bounds = (
        float(min(lon_values)),
        float(max(lon_values)),
        float(min(lat_values)),
        float(max(lat_values)),
    )
    metadata = {
        "version": 1,
        "l0_resolution": l0_resolution,
        "l0_unit": target_unit,
        "l2_resolution": adjusted_resolution,
        "l2_unit": target_unit,
        "l2_ratio_to_l0": ratio,
        "raw_meteo_resolution": raw["resolution"],
        "raw_meteo_unit": raw["unit"],
        "raw_meteo_x_resolution": raw.get("x_resolution"),
        "raw_meteo_y_resolution": raw.get("y_resolution"),
        "raw_meteo_degree_resolution": raw.get("degree_resolution", raw["resolution"]),
        "source_file": raw.get("source_file", ""),
        "extent_source": os.path.basename(str(extent_source).split("|")[0]),
        "wgs84_bounds": wgs84_bounds,
        "l2_bounds": header_bounds(l2_header),
        "crs_authid": qgis_crs_to_authid(target_crs),
        "l2_header": l2_header,
    }
    return {
        "metadata": metadata,
        "header": l2_header,
        "lon": lon_values,
        "lat": lat_values,
        "bounds": wgs84_bounds,
    }
