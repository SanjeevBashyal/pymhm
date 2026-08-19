# -*- coding: utf-8 -*-
"""Grid-resolution helpers for morphology, meteorology, and lat/lon files."""
from __future__ import annotations

import json
import math
import os
from decimal import Decimal, ROUND_FLOOR
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


def floor_to_precision(value, precision: int) -> float:
    """Floor a numeric value to a fixed number of decimal places."""
    numeric = float(value)
    if not math.isfinite(numeric):
        return numeric
    quantum = Decimal("1").scaleb(-int(precision))
    return float(
        Decimal(str(numeric)).quantize(
            quantum,
            rounding=ROUND_FLOOR,
        )
    )


def floor_cellsize(value, unit: str | None = None,
                   precision: int | None = None) -> float:
    """Floor a grid cellsize using the configured internal CRS precision."""
    if precision is None:
        precision = cellsize_precision_for_unit(unit)
    return floor_to_precision(value, precision)


def ceil_cellsize(value, unit: str | None = None,
                  precision: int | None = None) -> float:
    """Compatibility wrapper; cell sizes are floor-normalized now."""
    return floor_cellsize(value, unit, precision)


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

    exact_x_resolution = abs((extent.xMaximum() - extent.xMinimum()) / width)
    exact_y_resolution = abs((extent.yMaximum() - extent.yMinimum()) / height)
    exact_resolution = (exact_x_resolution + exact_y_resolution) / 2.0
    crs = layer.crs()
    unit = qgis_crs_unit_label(crs)
    x_resolution = ceil_cellsize(exact_x_resolution, unit)
    y_resolution = ceil_cellsize(exact_y_resolution, unit)
    resolution = ceil_cellsize((x_resolution + y_resolution) / 2.0, unit)
    return {
        "resolution": resolution,
        "x_resolution": x_resolution,
        "y_resolution": y_resolution,
        # Unrounded values. Grid construction must use these: flooring a
        # repeating cell size such as 1/1200 degrees drifts by a fraction of a
        # cell per row, which breaks the L0/L2 alignment over tall rasters.
        "exact_resolution": exact_resolution,
        "exact_x_resolution": exact_x_resolution,
        "exact_y_resolution": exact_y_resolution,
        "unit": unit,
        "crs_authid": qgis_crs_to_authid(crs),
        "cellsize_precision": cellsize_precision_for_unit(unit),
        "header": {
            "ncols": width,
            "nrows": height,
            "xllcorner": float(extent.xMinimum()),
            "yllcorner": float(extent.yMinimum()),
            "cellsize": float(resolution),
            "exact_cellsize": float(exact_resolution),
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
        # Record the unit but keep the written cell size exactly as stored. A
        # grid-defining cell size must never be re-rounded on the way in.
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
    return integer_resolution_factor(
        coarse_resolution,
        fine_resolution,
        unit=unit,
        tolerance=tolerance,
    ) is not None


def integer_resolution_factor(
        coarse_resolution: float,
        fine_resolution: float,
        unit: str | None = None,
        tolerance=None) -> int | None:
    """Return integer coarse/fine factor when the ratio is close enough.

    The inputs are compared as given. Rounding them to a fixed number of
    decimals first turns an exact ratio such as 0.1 / (1/1200) into 119.99999
    and the factor is then rejected.
    """
    if coarse_resolution <= 0 or fine_resolution <= 0:
        return None
    coarse_resolution = float(coarse_resolution)
    fine_resolution = float(fine_resolution)
    if tolerance is None:
        tolerance = 1e-7
    ratio = coarse_resolution / fine_resolution
    nearest = int(round(ratio))
    if nearest < 1:
        return None
    if abs(ratio - nearest) <= tolerance:
        return nearest
    return None


def possible_resolutions(
        fine_resolution: float,
        coarse_resolution: float,
        unit: str | None = None) -> list[float]:
    """Return all resolutions between fine and coarse compatible with both."""
    fine_resolution = float(fine_resolution)
    coarse_resolution = float(coarse_resolution)
    ratio = integer_resolution_factor(
        coarse_resolution,
        fine_resolution,
        unit=unit,
    )
    if ratio is None:
        return []
    values = []
    for divisor in range(1, ratio + 1):
        if ratio % divisor == 0:
            quotient = ratio // divisor
            values.append(float(coarse_resolution) / quotient)
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
    cellsize = float(cellsize)
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


def header_for_aligned_bounds(
        bounds,
        cellsize: float,
        anchor_x: float,
        anchor_y: float,
        unit: str | None = None) -> dict:
    """Create the smallest outward-snapped header aligned to an anchor grid."""
    cellsize = float(cellsize)
    xmin, xmax, ymin, ymax = bounds
    xll = anchor_x + math.floor(
        (min(xmin, xmax) - anchor_x) / cellsize
    ) * cellsize
    yll = anchor_y + math.floor(
        (min(ymin, ymax) - anchor_y) / cellsize
    ) * cellsize
    xur = anchor_x + math.ceil(
        (max(xmin, xmax) - anchor_x) / cellsize
    ) * cellsize
    yur = anchor_y + math.ceil(
        (max(ymin, ymax) - anchor_y) / cellsize
    ) * cellsize
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
    cellsize = float(cellsize)
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


def aligned_l0_l2_headers(
        bounds,
        l0_header: dict,
        multiplier: int,
        unit: str | None = None) -> tuple[dict, dict]:
    """Return the smallest common extent with exactly n-by-n L0 cells per L2 cell."""
    try:
        ratio = int(multiplier)
    except (TypeError, ValueError) as error:
        raise ValueError("The L2 resolution multiplier must be an integer.") from error
    if ratio < 1:
        raise ValueError("The L2 resolution multiplier must be at least 1.")

    l0_cellsize = float(l0_header["cellsize"])
    if not math.isfinite(l0_cellsize) or l0_cellsize <= 0:
        raise ValueError("The L0 cell size must be positive.")
    anchor_x = float(l0_header["xllcorner"])
    anchor_y = float(l0_header["yllcorner"])
    l2_cellsize = ratio * l0_cellsize
    xmin, xmax, ymin, ymax = map(float, bounds)

    def grid_index(value, anchor, *, upper):
        quotient = (value - anchor) / l2_cellsize
        nearest = round(quotient)
        tolerance = 1e-9 * max(1.0, abs(quotient))
        if math.isclose(quotient, nearest, rel_tol=0.0, abs_tol=tolerance):
            return int(nearest)
        return int(math.ceil(quotient) if upper else math.floor(quotient))

    x0 = grid_index(min(xmin, xmax), anchor_x, upper=False)
    x1 = grid_index(max(xmin, xmax), anchor_x, upper=True)
    y0 = grid_index(min(ymin, ymax), anchor_y, upper=False)
    y1 = grid_index(max(ymin, ymax), anchor_y, upper=True)
    resolved_unit = unit or l0_header.get("unit", "")
    l2 = {
        "ncols": max(1, x1 - x0),
        "nrows": max(1, y1 - y0),
        "xllcorner": anchor_x + x0 * l2_cellsize,
        "yllcorner": anchor_y + y0 * l2_cellsize,
        "cellsize": l2_cellsize,
        "nodata_value": NODATA_VALUE,
        "unit": resolved_unit,
        "cellsize_precision": cellsize_precision_for_unit(resolved_unit),
    }
    l0 = {
        **l2,
        "ncols": int(l2["ncols"]) * ratio,
        "nrows": int(l2["nrows"]) * ratio,
        "cellsize": l0_cellsize,
    }
    validate_l0_l2_alignment(l0, l2, ratio)
    return l0, l2


def l0_header_from_l2(
        l2_header: dict,
        l0_cellsize: float,
        multiplier: int) -> dict:
    """Derive an exact L0 header from a saved L2 header and multiplier."""
    ratio = int(multiplier)
    l2_cellsize = float(l2_header["cellsize"])
    tolerance = max(abs(l2_cellsize), 1.0) * 1e-9
    if ratio < 1 or not math.isclose(
            l2_cellsize, ratio * float(l0_cellsize),
            rel_tol=0.0, abs_tol=tolerance):
        raise ValueError("The saved L2 resolution is not an exact multiple of L0.")
    l0 = {
        **dict(l2_header),
        "ncols": int(l2_header["ncols"]) * ratio,
        "nrows": int(l2_header["nrows"]) * ratio,
        "cellsize": float(l0_cellsize),
    }
    validate_l0_l2_alignment(l0, l2_header, ratio)
    return l0


def validate_l0_l2_alignment(l0_header: dict, l2_header: dict, multiplier: int) -> None:
    """Raise when two headers do not share an exact integer L0/L2 contract."""
    ratio = int(multiplier)
    if ratio < 1:
        raise ValueError("The L2 resolution multiplier must be at least 1.")
    l0_cellsize = float(l0_header["cellsize"])
    l2_cellsize = float(l2_header["cellsize"])
    tolerance = max(abs(l2_cellsize), abs(l0_cellsize), 1.0) * 1e-9
    if not math.isclose(
            l2_cellsize, ratio * l0_cellsize,
            rel_tol=0.0, abs_tol=tolerance):
        raise ValueError("L2 cell size must equal multiplier x L0 cell size.")
    if (
            int(l0_header["ncols"]) != int(l2_header["ncols"]) * ratio
            or int(l0_header["nrows"]) != int(l2_header["nrows"]) * ratio):
        raise ValueError("L0 dimensions must equal L2 dimensions x multiplier.")
    for key in ("xllcorner", "yllcorner"):
        if not math.isclose(
                float(l0_header[key]), float(l2_header[key]),
                rel_tol=0.0, abs_tol=tolerance):
            raise ValueError("L0 and L2 must share the same lower-left boundary.")
    if any(
            not math.isclose(left, right, rel_tol=0.0, abs_tol=tolerance)
            for left, right in zip(
                header_bounds(l0_header), header_bounds(l2_header))):
        raise ValueError("L0 and L2 must share the same outer extent.")


def axes_match_header(x_values, y_values, header: dict) -> bool:
    """Return True when source axes sit exactly on the header cell grid."""
    import numpy as np

    cellsize = float(header["cellsize"])
    tolerance = abs(cellsize) * 1e-6
    target_x, target_y = header_center_coordinates(header)
    pairs = ((x_values, target_x), (y_values, sorted(target_y)))
    for values, target in pairs:
        values = np.sort(np.asarray(values, dtype="float64").ravel())
        if values.size < 1:
            return False
        if values.size > 1:
            spacing = float(np.abs(np.diff(values)).max())
            if abs(spacing - cellsize) > tolerance:
                return False
        offset = (float(values[0]) - float(target[0])) / cellsize
        if abs(offset - round(offset)) > 1e-6:
            return False
    return True


def reindex_to_header(da, x_dim: str, y_dim: str, header: dict, fill_value=None):
    """Place an already aligned array on the exact header grid by cell lookup."""
    if fill_value is None:
        fill_value = float("nan")
    x_values, y_values = header_center_coordinates(header)
    return da.reindex(
        {x_dim: x_values, y_dim: y_values},
        method="nearest",
        tolerance=abs(float(header["cellsize"])) * 1e-6,
        fill_value=fill_value,
    )


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


def merged_watershed_extent_in_crs(dialog, target_crs):
    """Return the merged watershed extent in the requested CRS."""
    from qgis.core import QgsVectorLayer

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

    return None, ""


def watershed_extent_in_crs(dialog, target_crs):
    """Return a watershed extent, retaining the historic DEM fallbacks."""
    from qgis.core import QgsRasterLayer

    extent, source = merged_watershed_extent_in_crs(dialog, target_crs)
    if extent is not None:
        return extent, source

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


def target_lon_lat_mesh_from_header(header: dict, source_crs):
    """Return exact WGS84 coordinates for every target grid cell centre."""
    import numpy as np

    x_values, y_values = header_center_coordinates(header)
    x_mesh, y_mesh = np.meshgrid(
        np.asarray(x_values, dtype="float64"),
        np.asarray(y_values, dtype="float64"),
    )
    if source_crs is None or not source_crs.isValid() or source_crs.isGeographic():
        return x_mesh, y_mesh

    try:
        from pyproj import Transformer

        source = qgis_crs_to_authid(source_crs)
        if not source:
            to_wkt = getattr(source_crs, "toWkt", None)
            source = to_wkt() if callable(to_wkt) else ""
        transform = Transformer.from_crs(source, "EPSG:4326", always_xy=True)
        return transform.transform(x_mesh, y_mesh)
    except Exception:
        from qgis.core import (QgsCoordinateReferenceSystem,
                               QgsCoordinateTransform, QgsPointXY, QgsProject)

        transform = QgsCoordinateTransform(
            source_crs,
            QgsCoordinateReferenceSystem("EPSG:4326"),
            QgsProject.instance(),
        )
        transform.setBallparkTransformsAreAppropriate(True)
        lon = np.empty_like(x_mesh)
        lat = np.empty_like(y_mesh)
        for row in range(x_mesh.shape[0]):
            for column in range(x_mesh.shape[1]):
                point = transform.transform(QgsPointXY(
                    float(x_mesh[row, column]),
                    float(y_mesh[row, column]),
                ))
                lon[row, column] = point.x()
                lat[row, column] = point.y()
        return lon, lat


def build_meteo_l2_grid(
        dialog,
        multiplier: int,
        raw_metadata: dict | None = None) -> dict:
    """Build the requested integer-multiple L2 grid on the L0 anchor."""
    l0_info = getattr(dialog, "_grid_l0_info", None)
    if not l0_info and hasattr(dialog, "filled_dem_resolution_info"):
        l0_info = dialog.filled_dem_resolution_info()
    if not l0_info:
        raise ValueError("Filled DEM is required before preparing meteorology data.")
    x_resolution = float(
        l0_info.get("exact_x_resolution")
        or l0_info.get("x_resolution", l0_info["resolution"])
    )
    y_resolution = float(
        l0_info.get("exact_y_resolution")
        or l0_info.get("y_resolution", l0_info["resolution"])
    )
    if abs(x_resolution - y_resolution) > max(
            abs(x_resolution), abs(y_resolution), 1.0) * 1e-6:
        raise ValueError(
            "The filled DEM must have square L0 cells before preparing "
            f"meteorology data (x={x_resolution}, y={y_resolution})."
        )

    target_crs = project_crs_from_dialog(dialog)
    if target_crs is None or not target_crs.isValid():
        raise ValueError("Please set a valid processing CRS before preparing meteorology data.")

    extent, extent_source = merged_watershed_extent_in_crs(dialog, target_crs)
    if extent is None:
        raise ValueError(
            "The merged watershed is required before preparing meteorology data."
        )

    target_unit = qgis_crs_unit_label(target_crs)
    # Anchor the grid on the filled DEM's own cell size, unrounded. A floored
    # cell size drifts against the source grid and the L0 window copy then
    # rejects every raster as misaligned.
    l0_resolution = float(
        l0_info.get("exact_resolution")
        or ceil_cellsize(l0_info["resolution"], target_unit)
    )
    try:
        ratio = int(multiplier)
    except (TypeError, ValueError):
        raise ValueError("The L2 resolution multiplier must be an integer.")
    if ratio < 1:
        raise ValueError("The L2 resolution multiplier must be at least 1.")
    l0_header = l0_info.get("header") or {}
    if "xllcorner" not in l0_header or "yllcorner" not in l0_header:
        raise ValueError("Could not determine the L0 grid anchor from the filled DEM.")
    exact_l0_header, l2_header = aligned_l0_l2_headers(
        extent_to_bounds(extent),
        {
            **l0_header,
            "cellsize": l0_resolution,
            "unit": target_unit,
        },
        ratio,
        target_unit,
    )
    lon_mesh, lat_mesh = target_lon_lat_mesh_from_header(l2_header, target_crs)
    middle_row = lon_mesh.shape[0] // 2
    middle_column = lon_mesh.shape[1] // 2
    lon_values = lon_mesh[middle_row, :].tolist()
    lat_values = lat_mesh[:, middle_column].tolist()
    wgs84_bounds = (
        float(lon_mesh.min()),
        float(lon_mesh.max()),
        float(lat_mesh.min()),
        float(lat_mesh.max()),
    )
    metadata = {
        "version": 2,
        "l0_resolution": l0_resolution,
        "l0_unit": target_unit,
        "l2_resolution": float(l2_header["cellsize"]),
        "l2_unit": target_unit,
        "l2_ratio_to_l0": ratio,
        "extent_source": os.path.basename(str(extent_source).split("|")[0]),
        "wgs84_bounds": wgs84_bounds,
        "l2_bounds": header_bounds(l2_header),
        "crs_authid": qgis_crs_to_authid(target_crs),
        "l2_header": l2_header,
        "l0_header": exact_l0_header,
    }
    if raw_metadata:
        metadata.update({
            "raw_meteo_resolution": raw_metadata.get("resolution"),
            "raw_meteo_unit": raw_metadata.get("unit", ""),
            "raw_meteo_x_resolution": raw_metadata.get("x_resolution"),
            "raw_meteo_y_resolution": raw_metadata.get("y_resolution"),
            "source_file": raw_metadata.get("source_file", ""),
        })
    return {
        "metadata": metadata,
        "header": l2_header,
        "lon": lon_values,
        "lat": lat_values,
        "sample_lon": lon_mesh,
        "sample_lat": lat_mesh,
        "bounds": wgs84_bounds,
    }
