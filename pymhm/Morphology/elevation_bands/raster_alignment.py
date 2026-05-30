# -*- coding: utf-8 -*-
"""Reference-grid comparison and nearest-neighbour raster alignment helpers."""
from ..common import os


class RasterAlignmentMixin:
    """Reference-grid comparison and nearest-neighbour raster alignment helpers."""

    def _raster_matches_reference(self, raster_reference, target_reference):
        """Check whether two rasters share shape, geotransform, and projection."""
        if raster_reference["rows"] != target_reference["rows"]:
            return False
        if raster_reference["cols"] != target_reference["cols"]:
            return False

        for current, target in zip(
                raster_reference["geotransform"],
                target_reference["geotransform"]):
            if abs(float(current) - float(target)) > 1e-7:
                return False

        raster_projection = raster_reference.get("projection") or ""
        target_projection = target_reference.get("projection") or ""
        if raster_projection and target_projection and raster_projection != target_projection:
            return False

        return True

    def _source_path_for_gdal(self, layer):
        """Return a GDAL-openable source path from a QGIS layer source."""
        source = layer.source()
        source_path = source.split("|")[0]
        return source_path if os.path.exists(source_path) else source

    def _reference_bounds(self, reference):
        """Return GDAL output bounds from a reference raster."""
        gt = reference["geotransform"]
        cols = int(reference["cols"])
        rows = int(reference["rows"])

        corners = [
            (0, 0),
            (cols, 0),
            (0, rows),
            (cols, rows),
        ]
        xs = []
        ys = []
        for col, row in corners:
            xs.append(gt[0] + col * gt[1] + row * gt[2])
            ys.append(gt[3] + col * gt[4] + row * gt[5])

        return (min(xs), min(ys), max(xs), max(ys))

    def _align_raster_to_reference(self, input_path, output_path, reference, deps, nodata=-9999):
        """Nearest-neighbour warp a categorical raster to a reference grid."""
        gdal = deps["gdal"]
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        if os.path.exists(output_path):
            try:
                os.remove(output_path)
            except Exception as e:
                self.log_message(
                    f"ERROR: Could not remove existing aligned raster {output_path}: {e}")
                return False

        warp_kwargs = {
            "format": "GTiff",
            "outputBounds": self._reference_bounds(reference),
            "width": int(reference["cols"]),
            "height": int(reference["rows"]),
            "resampleAlg": gdal.GRA_NearestNeighbour,
            "dstNodata": nodata,
            "outputType": gdal.GDT_Float32,
            "creationOptions": ["COMPRESS=LZW"],
        }
        if reference.get("projection"):
            warp_kwargs["dstSRS"] = reference["projection"]

        try:
            ds = gdal.Warp(
                output_path,
                input_path,
                options=gdal.WarpOptions(**warp_kwargs)
            )
            if ds is None:
                self.log_message(
                    f"ERROR: Failed to align raster to elevation-band grid: {input_path}")
                return False
            ds = None
        except Exception as e:
            self.log_message(f"ERROR: GDAL warp failed while aligning land cover: {e}")
            return False

        if os.path.exists(output_path):
            self.mark_output_prepared(
                output_path,
                name=os.path.basename(output_path),
                loaded=False
            )
            return True
        return False
