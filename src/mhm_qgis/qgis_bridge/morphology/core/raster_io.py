# -*- coding: utf-8 -*-
"""Raster array IO, cell-area helpers, and pyflwdir raster context builders."""
from __future__ import annotations

from typing import Any

from ..common import os
from .base import BaseProcessingMixin
from ....applications import pyflwdir_handler
from .dependencies import PythonDependencyMixin


class RasterIOMixin(PythonDependencyMixin, BaseProcessingMixin):
    """Raster array IO, cell-area helpers, and pyflwdir raster context builders."""

    def _read_raster_array(
            self,
            raster_path: str,
            as_float: bool = False) -> dict[str, Any] | None:
        """Read a single-band raster into a NumPy array with GDAL metadata."""
        deps = self._get_python_morphology_deps()
        if not deps:
            return None

        gdal = deps["gdal"]
        ds = gdal.Open(raster_path)
        if ds is None:
            self.log_message(f"ERROR: Could not open raster: {raster_path}")
            return None

        band = ds.GetRasterBand(1)
        array = band.ReadAsArray()
        if as_float:
            array = array.astype(deps["np"].float32)

        metadata = {
            "array": array,
            "nodata": band.GetNoDataValue(),
            "geotransform": ds.GetGeoTransform(),
            "projection": ds.GetProjection(),
            "rows": ds.RasterYSize,
            "cols": ds.RasterXSize,
        }

        band = None
        ds = None
        return metadata

    def _read_raster_grid_metadata(self, raster_path: str) -> dict[str, Any] | None:
        """Read raster grid metadata without loading band values."""
        deps = self._get_python_morphology_deps()
        if not deps:
            return None

        gdal = deps["gdal"]
        ds = gdal.Open(raster_path)
        if ds is None:
            self.log_message(f"ERROR: Could not open raster: {raster_path}")
            return None

        metadata = {
            "geotransform": ds.GetGeoTransform(),
            "projection": ds.GetProjection(),
            "rows": ds.RasterYSize,
            "cols": ds.RasterXSize,
        }
        ds = None
        return metadata

    def _projection_matches(
            self,
            projection_a: str | None,
            projection_b: str | None,
            deps: dict[str, Any]) -> bool:
        """Return True when two WKT projections describe the same CRS."""
        if not projection_a and not projection_b:
            return True
        if not projection_a or not projection_b:
            return False

        osr = deps["osr"]
        spatial_ref_a = osr.SpatialReference()
        spatial_ref_b = osr.SpatialReference()
        if spatial_ref_a.ImportFromWkt(projection_a) != 0:
            return projection_a == projection_b
        if spatial_ref_b.ImportFromWkt(projection_b) != 0:
            return projection_a == projection_b
        return bool(spatial_ref_a.IsSame(spatial_ref_b))

    def _geotransform_matches(
            self,
            transform_a: tuple,
            transform_b: tuple,
            tolerance: float = 1e-9) -> bool:
        """Return True when two affine transforms are equal within tolerance."""
        for value_a, value_b in zip(transform_a, transform_b):
            scale = max(abs(float(value_a)), abs(float(value_b)), 1.0)
            if abs(float(value_a) - float(value_b)) > tolerance * scale:
                return False
        return True

    def _raster_grid_matches(
            self,
            raster_path: str,
            reference_path: str) -> bool:
        """Return True when a raster has the same grid and CRS as a reference raster."""
        if not raster_path or not reference_path:
            return False
        if not os.path.exists(raster_path) or not os.path.exists(reference_path):
            return False

        deps = self._get_python_morphology_deps()
        if not deps:
            return False

        raster = self._read_raster_grid_metadata(raster_path)
        reference = self._read_raster_grid_metadata(reference_path)
        if not raster or not reference:
            return False

        if raster["rows"] != reference["rows"] or raster["cols"] != reference["cols"]:
            return False
        if not self._geotransform_matches(
                raster["geotransform"], reference["geotransform"]):
            return False
        return self._projection_matches(
            raster.get("projection"), reference.get("projection"), deps)

    def _write_raster_array(
            self,
            output_path: str,
            array: Any,
            reference: dict[str, Any],
            nodata: float | int | None = None,
            gdal_type: Any | None = None) -> bool:
        """Write a NumPy array as a single-band GeoTIFF using reference metadata."""
        deps = self._get_python_morphology_deps()
        if not deps:
            return False

        gdal = deps["gdal"]
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        if os.path.exists(output_path):
            try:
                os.remove(output_path)
            except Exception as e:
                self.log_message(f"ERROR: Could not remove existing raster {output_path}: {e}")
                return False

        if gdal_type is None:
            gdal_type = gdal.GDT_Float32

        driver = gdal.GetDriverByName("GTiff")
        ds = driver.Create(
            output_path,
            int(reference["cols"]),
            int(reference["rows"]),
            1,
            gdal_type,
            options=["COMPRESS=LZW"]
        )
        if ds is None:
            self.log_message(f"ERROR: Could not create raster: {output_path}")
            return False

        ds.SetGeoTransform(reference["geotransform"])
        if reference.get("projection"):
            ds.SetProjection(reference["projection"])

        band = ds.GetRasterBand(1)
        if nodata is not None:
            band.SetNoDataValue(float(nodata))
        band.WriteArray(array)
        band.FlushCache()
        band = None
        ds = None
        if os.path.exists(output_path):
            self.mark_output_prepared(
                output_path,
                name=os.path.basename(output_path),
                loaded=False
            )
            return True
        return False

    @staticmethod
    def _projection_is_geographic(projection: str | None, osr: Any = None) -> bool:
        """Report whether a WKT projection is a geographic (lat/lon) CRS."""
        return pyflwdir_handler.is_geographic(projection)

    @staticmethod
    def _reference_cell_area_m2(reference: dict[str, Any]) -> float:
        """Approximate raster cell area in square metres for channel thresholds."""
        return pyflwdir_handler.cell_area_m2(reference)

    @staticmethod
    def _normalise_dem_array(
            array: Any,
            nodata: float | int | None) -> tuple[Any, Any, float]:
        """Prepare DEM data for pyflwdir using a finite -9999 nodata value."""
        dem, invalid = pyflwdir_handler.normalise_dem(array, nodata)
        return dem, invalid, pyflwdir_handler.NODATA

    def _build_flwdir_from_filled_dem(self) -> dict[str, Any] | None:
        """Build the pyflwdir flow grid from the filled DEM."""
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            self.log_message("ERROR: Filled DEM is missing. Run Fill DEM first.")
            return None

        reference = self._read_raster_array(self.filled_dem_path, as_float=True)
        if not reference:
            return None
        try:
            return pyflwdir_handler.flwdir_context(reference)
        except pyflwdir_handler.MissingDependencies as error:
            self.log_message(f"ERROR: {error}")
            return None
