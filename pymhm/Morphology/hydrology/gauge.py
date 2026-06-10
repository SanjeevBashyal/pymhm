# -*- coding: utf-8 -*-
"""Gauge-position raster preparation from snapped pour points."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    QMessageBox,
    QgsCoordinateTransform,
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
)
from ..core.predecessors import PredecessorMixin
from ..watershed.pour_point_workflow import PourPointWorkflowMixin
from .outlets import (
    OutletCountMixin,
    StationIdError,
    find_station_id_field,
    station_id_int,
    station_id_text,
)


class GaugePositionMixin(PourPointWorkflowMixin, OutletCountMixin, PredecessorMixin):
    """Gauge-position raster preparation from snapped pour points."""

    def process_gauge_position(self) -> None:
        """Prepare the mHM idgauges raster from snapped gauged outlets."""
        if not self._ensure_filled_dem(self.fill_dem):
            return

        if not self.snapped_points_path or not os.path.exists(self.snapped_points_path):
            if not self._ensure_snapped_points(
                    self.snap_points,
                    self.process_channel_network,
                    self.process_flow_accumulation,
                    self.fill_dem):
                return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        self.gauge_position_path = os.path.join(
            geometry_folder, "2_gauge_position.tif")

        snapped_layer = QgsVectorLayer(
            self.snapped_points_path,
            "2_pour_points_snapped",
            "ogr",
        )
        if not snapped_layer.isValid():
            QMessageBox.critical(
                self.dialog,
                "Invalid Snapped Points",
                f"Could not read snapped pour points:\n{self.snapped_points_path}")
            self.log_message(
                f"ERROR: Could not read snapped pour points: {self.snapped_points_path}")
            return

        try:
            station_field = find_station_id_field(snapped_layer)
            gauge_features = self._gauge_features(snapped_layer, station_field)
        except StationIdError as e:
            QMessageBox.critical(self.dialog, "Missing STATION_ID", str(e))
            self.log_message(f"ERROR: {e}")
            return

        if hasattr(self, "update_gauged_outlet_count"):
            self.update_gauged_outlet_count()

        deps = self._get_python_morphology_deps()
        if not deps:
            return

        if (
                self.gauge_position_path
                and os.path.exists(self.gauge_position_path)
                and self._existing_gauge_position_matches(
                    self.gauge_position_path,
                    [item["station_int"] for item in gauge_features],
                    deps)):
            self.log_message(
                "Gauge Position already exists with matching STATION_ID values. Loading existing file...")
            self.load_layer(self.gauge_position_path, "2_Gauge_Position")
            return

        self.log_message(
            "Processing Gauge Position with STATION_ID burn values...")
        self.log_message(
            f"Gauge STATION_ID values: {', '.join(item['station_text'] for item in gauge_features)}")

        context = self._build_flwdir_from_filled_dem()
        if not context:
            self.log_message("Gauge Position processing failed.")
            return

        np = context["deps"]["np"]
        gdal = context["deps"]["gdal"]
        reference = context["reference"]
        stream_mask = self._stream_mask_from_context(context)
        if stream_mask is None or not np.any(stream_mask):
            QMessageBox.critical(
                self.dialog,
                "No Stream Cells",
                "Could not identify stream cells from the prepared flow accumulation.")
            self.log_message(
                "ERROR: Could not identify stream cells from flow accumulation.")
            return

        gauge_array = np.full(
            (int(reference["rows"]), int(reference["cols"])),
            -9999,
            dtype=np.int32,
        )

        point_transform = self._snapped_to_raster_transform(snapped_layer)
        for item in gauge_features:
            point = item["geometry"].asPoint()
            if point_transform is not None:
                point = point_transform.transform(point)

            row, col = self._point_to_row_col(
                point.x(),
                point.y(),
                reference["geotransform"],
                context["deps"],
            )
            snapped_row, snapped_col, distance_cells = self._nearest_stream_cell(
                row,
                col,
                stream_mask,
                np,
            )

            existing = gauge_array[snapped_row, snapped_col]
            if existing != -9999 and existing != item["station_int"]:
                self.log_message(
                    "WARNING: Multiple gauges snapped to the same raster cell; "
                    f"overwriting {existing} with {item['station_int']}.")

            gauge_array[snapped_row, snapped_col] = item["station_int"]
            self.log_message(
                f"  STATION_ID {item['station_text']} burned at row={snapped_row}, "
                f"col={snapped_col}, stream snap distance={distance_cells:.2f} cell(s).")

        if self._write_raster_array(
                self.gauge_position_path,
                gauge_array,
                reference,
                nodata=-9999,
                gdal_type=gdal.GDT_Int32):
            self.load_layer(self.gauge_position_path, "2_Gauge_Position")
            self.log_message(
                "Gauge Position processing completed successfully.")
        else:
            self.log_message("Gauge Position processing failed.")

    def _gauge_features(self, layer, station_field: str) -> list[dict]:
        """Return valid snapped gauge features with station ids."""
        gauge_features = []
        seen_station_ids = set()

        for feature in layer.getFeatures():
            station_text = station_id_text(feature.attribute(station_field))
            station_value = station_id_int(feature.attribute(station_field))
            if station_text in seen_station_ids:
                raise StationIdError(
                    f"Duplicate STATION_ID '{station_text}' found in snapped pour points.")
            seen_station_ids.add(station_text)

            geometry = feature.geometry()
            if geometry is None or geometry.isEmpty():
                raise StationIdError(
                    f"STATION_ID '{station_text}' has an empty geometry.")

            gauge_features.append({
                "station_text": station_text,
                "station_int": station_value,
                "geometry": geometry,
            })

        if not gauge_features:
            raise StationIdError(
                "The selected pour points layer does not contain any gauged outlets.")

        return gauge_features

    def _stream_mask_from_context(self, context: dict) -> object | None:
        """Recreate the channel network stream mask on the filled DEM grid."""
        deps = context["deps"]
        np = deps["np"]
        flwdir = context["flwdir"]
        reference = context["reference"]
        invalid_mask = context["invalid_mask"]

        flow_accumulation = flwdir.upstream_area(unit="cell")
        flow_accumulation[invalid_mask] = 0

        cell_area_m2 = self._reference_cell_area_m2(reference, deps)
        threshold_cells = max(1, int(round(10000000.0 / cell_area_m2)))
        valid_accumulation = flow_accumulation[flow_accumulation > 0]
        if valid_accumulation.size == 0:
            return None

        max_accumulation = int(np.nanmax(valid_accumulation))
        if threshold_cells > max_accumulation:
            threshold_cells = max(
                1,
                int(np.nanpercentile(valid_accumulation, 95)),
            )
            self.log_message(
                "Gauge stream threshold exceeded basin accumulation. "
                f"Using 95th percentile threshold: {threshold_cells} cells.")

        return flow_accumulation >= threshold_cells

    def _snapped_to_raster_transform(self, snapped_layer: object) -> object | None:
        """Return a coordinate transform when snapped points and DEM CRS differ."""
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            return None

        source_crs = snapped_layer.crs()
        target_crs = filled_dem_layer.crs()
        if not source_crs.isValid() or not target_crs.isValid():
            return None
        if source_crs.authid() == target_crs.authid():
            return None

        transform = QgsCoordinateTransform(
            source_crs,
            target_crs,
            QgsProject.instance(),
        )
        transform.setBallparkTransformsAreAppropriate(True)
        return transform

    def _point_to_row_col(
            self,
            x: float,
            y: float,
            geotransform: tuple,
            deps: dict) -> tuple[int, int]:
        """Convert map coordinates to raster row/column indices."""
        transform = deps["Affine"].from_gdal(*geotransform)
        col_float, row_float = (~transform) * (float(x), float(y))
        np = deps["np"]
        return int(np.floor(row_float)), int(np.floor(col_float))

    def _nearest_stream_cell(
            self,
            row: int,
            col: int,
            stream_mask: object,
            np: object,
            max_radius: int = 10) -> tuple[int, int, float]:
        """Return the nearest stream cell to a candidate row/column."""
        rows, cols = stream_mask.shape
        row = min(max(int(row), 0), rows - 1)
        col = min(max(int(col), 0), cols - 1)

        for radius in range(max_radius + 1):
            r0 = max(0, row - radius)
            r1 = min(rows, row + radius + 1)
            c0 = max(0, col - radius)
            c1 = min(cols, col + radius + 1)
            candidates = np.argwhere(stream_mask[r0:r1, c0:c1])
            if candidates.size:
                candidates[:, 0] += r0
                candidates[:, 1] += c0
                return self._closest_candidate(row, col, candidates, np)

        candidates = np.argwhere(stream_mask)
        return self._closest_candidate(row, col, candidates, np)

    def _closest_candidate(
            self,
            row: int,
            col: int,
            candidates: object,
            np: object) -> tuple[int, int, float]:
        """Choose the closest row/column candidate and report cell distance."""
        deltas = candidates.astype(float)
        deltas[:, 0] -= row
        deltas[:, 1] -= col
        distances = np.hypot(deltas[:, 0], deltas[:, 1])
        best_index = int(np.argmin(distances))
        best_row = int(candidates[best_index, 0])
        best_col = int(candidates[best_index, 1])
        return best_row, best_col, float(distances[best_index])

    def _existing_gauge_position_matches(
            self,
            path: str,
            station_values: list[int],
            deps: dict) -> bool:
        """Return True when an existing idgauges raster already has these IDs."""
        reference = self._read_raster_array(path, as_float=False)
        if not reference:
            return False

        np = deps["np"]
        array = reference["array"]
        nodata = reference["nodata"]
        valid = np.isfinite(array)
        if nodata is not None and np.isfinite(nodata):
            valid &= ~np.isclose(array, nodata)

        values = {
            int(value)
            for value in np.asarray(array[valid]).ravel()
            if int(value) > 0
        }
        return values == set(station_values)
