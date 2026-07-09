# -*- coding: utf-8 -*-
"""Preparation guards for prerequisite morphology outputs."""
from __future__ import annotations

from collections.abc import Callable

from ..common import (
    os,
    QgsRasterLayer,
    QgsVectorLayer,
)
from .base import BaseProcessingMixin
from .naming import NamingAndRangeMixin


class PredecessorMixin(BaseProcessingMixin, NamingAndRangeMixin):
    """Preparation guards for prerequisite morphology outputs."""

    def _restore_existing_path(
            self,
            attr_name: str,
            *filenames: str) -> str | None:
        """Restore an output path attribute if a known file already exists."""
        current_path = getattr(self, attr_name, None)
        if current_path and self.is_output_prepared(current_path):
            return current_path

        for filename in filenames:
            path = self._geometry_path(filename)
            if self.is_output_prepared(path):
                setattr(self, attr_name, path)
                return path

        return None

    def _remove_stale_raster_output(self, path: str) -> None:
        """Remove a stale raster output and common GDAL sidecars."""
        for candidate in (path, f"{path}.aux.xml"):
            if candidate and os.path.exists(candidate):
                try:
                    os.remove(candidate)
                except Exception as e:
                    self.log_message(f"WARNING: Could not remove stale file {candidate}: {e}")

    def _remove_stale_vector_output(self, path: str) -> None:
        """Remove a stale vector output and sidecars."""
        if hasattr(self, "_remove_vector_dataset"):
            self._remove_vector_dataset(path)
            return
        base, _extension = os.path.splitext(path)
        for extension in (".shp", ".shx", ".dbf", ".prj", ".cpg", ".qpj", ".fix"):
            candidate = base + extension
            if os.path.exists(candidate):
                try:
                    os.remove(candidate)
                except Exception as e:
                    self.log_message(f"WARNING: Could not remove stale file {candidate}: {e}")

    def _filled_dem_matches_current_input(self, path: str) -> bool:
        """Return True when a restored filled DEM still matches the selected DEM."""
        if not path or not os.path.exists(path):
            return False
        if not hasattr(self, "_existing_filled_dem_matches_current_dem"):
            return True

        dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        if not dem_layer:
            return False
        if hasattr(self, "check_and_reproject_dem"):
            dem_layer = self.check_and_reproject_dem(dem_layer)
        return self._existing_filled_dem_matches_current_dem(path, dem_layer)

    def _raster_output_matches_filled_dem(self, path: str) -> bool:
        """Return True when a raster output is aligned with the filled DEM."""
        if not path or not self.filled_dem_path:
            return False
        if not hasattr(self, "_raster_grid_matches"):
            return True
        return self._raster_grid_matches(path, self.filled_dem_path)

    def _vector_crs_matches_raster(self, vector_path: str, raster_path: str) -> bool:
        """Return True when a vector layer CRS matches a raster layer CRS."""
        if not vector_path or not raster_path:
            return False
        if not os.path.exists(vector_path) or not os.path.exists(raster_path):
            return False

        vector_layer = QgsVectorLayer(vector_path, "Existing_Vector", "ogr")
        raster_layer = QgsRasterLayer(raster_path, "Reference_Raster")
        if not vector_layer.isValid() or not raster_layer.isValid():
            return False

        vector_crs = vector_layer.crs()
        raster_crs = raster_layer.crs()
        if vector_crs.isValid() != raster_crs.isValid():
            return False
        if vector_crs.isValid() and vector_crs.authid() != raster_crs.authid():
            return False
        return True

    def _snapped_points_match_channel_network(
            self,
            snapped_path: str,
            channel_network_path: str) -> bool:
        """Return True when snapped points are current and include snap status."""
        if not snapped_path or not channel_network_path:
            return False
        if not os.path.exists(snapped_path) or not os.path.exists(channel_network_path):
            return False

        snapped_layer = QgsVectorLayer(snapped_path, "Existing_Snapped_Points", "ogr")
        channel_layer = QgsVectorLayer(channel_network_path, "Existing_Channel_Network", "ogr")
        if not snapped_layer.isValid() or not channel_layer.isValid():
            return False

        snapped_crs = snapped_layer.crs()
        channel_crs = channel_layer.crs()
        if snapped_crs.isValid() != channel_crs.isValid():
            return False
        if snapped_crs.isValid() and snapped_crs.authid() != channel_crs.authid():
            return False

        field_names = snapped_layer.fields().names()
        status_field = None
        for candidate in ("snap_status", "snap_statu"):
            if candidate in field_names:
                status_field = candidate
                break
        if status_field is None:
            return False

        valid_statuses = {"high_order", "closest", "failed"}
        found_feature = False
        for feature in snapped_layer.getFeatures():
            found_feature = True
            status = feature.attribute(status_field)
            if status is None or str(status) not in valid_statuses:
                return False
        return found_feature

    def _ensure_filled_dem(self, fill_dem: Callable[[], object]) -> bool:
        """Create the filled DEM if it has not been generated yet."""
        existing_path = self._restore_existing_path(
            "filled_dem_path", "1_dem_filled.tif", "1_dem_filled.sdat")
        if existing_path:
            if self._filled_dem_matches_current_input(existing_path):
                return True
            self.log_message(
                "Existing filled DEM does not match the current DEM CRS/grid. "
                "Recreating Fill DEM and dependent outputs.")
            if hasattr(self, "_remove_stale_flow_outputs"):
                self._remove_stale_flow_outputs()
            else:
                self._remove_stale_raster_output(existing_path)
            self.filled_dem_path = None

        self.log_message("Filled DEM is missing. Running Fill DEM first...")
        self.without_layer_loading(fill_dem)
        return bool(self.filled_dem_path and os.path.exists(self.filled_dem_path))

    def _ensure_flow_accumulation(
            self,
            process_flow_accumulation: Callable[[], object],
            fill_dem: Callable[[], object]) -> bool:
        """Create flow accumulation if it has not been generated yet."""
        if not self._ensure_filled_dem(fill_dem):
            return False

        existing_path = self._restore_existing_path(
            "flow_accumulation_path", "2_flow_accumulation.tif")
        if existing_path:
            if self._raster_output_matches_filled_dem(existing_path):
                area_path = self._restore_existing_path(
                    "flow_accumulation_area_path",
                    "2_flow_accumulation_area.tif",
                    "2_flow_accumulation_area.sdat"
                )
                if area_path and not self._raster_output_matches_filled_dem(area_path):
                    self.log_message(
                        "Existing flow accumulation area raster is stale. Recreating flow accumulation.")
                    self._remove_stale_raster_output(area_path)
                    self.flow_accumulation_area_path = None
                else:
                    return True
            self.log_message(
                "Existing flow accumulation raster is stale. Recreating flow accumulation.")
            self._remove_stale_raster_output(existing_path)
            self.flow_accumulation_path = None

        self.log_message("Flow accumulation is missing. Running Flow Accumulation first...")
        self.without_layer_loading(process_flow_accumulation)
        return bool(self.flow_accumulation_path and os.path.exists(self.flow_accumulation_path))

    def _ensure_flow_direction(
            self,
            process_flow_direction: Callable[[], object],
            fill_dem: Callable[[], object]) -> bool:
        """Create flow direction if it has not been generated yet."""
        if not self._ensure_filled_dem(fill_dem):
            return False

        existing_path = self._restore_existing_path(
            "flow_direction_path", "2_flow_direction.tif")
        if existing_path:
            if self._raster_output_matches_filled_dem(existing_path):
                return True
            self.log_message(
                "Existing flow direction raster is stale. Recreating flow direction.")
            self._remove_stale_raster_output(existing_path)
            self.flow_direction_path = None

        self.log_message("Flow direction is missing. Running Flow Direction first...")
        self.without_layer_loading(process_flow_direction)
        return bool(self.flow_direction_path and os.path.exists(self.flow_direction_path))

    def _ensure_channel_network(
            self,
            process_channel_network: Callable[[], object],
            process_flow_accumulation: Callable[[], object],
            fill_dem: Callable[[], object]) -> bool:
        """Create the channel network if it has not been generated yet."""
        if not self._ensure_flow_accumulation(process_flow_accumulation, fill_dem):
            return False

        existing_path = self._restore_existing_path(
            "channel_network_vector_path", "2_channel_network.shp")
        if existing_path:
            if self._vector_crs_matches_raster(existing_path, self.filled_dem_path):
                return True
            self.log_message(
                "Existing channel network CRS does not match the filled DEM. "
                "Recreating channel network.")
            self._remove_stale_vector_output(existing_path)
            self.channel_network_vector_path = None

        self.log_message("Channel network is missing. Running Create Network first...")
        self.without_layer_loading(process_channel_network)
        return bool(self.channel_network_vector_path and os.path.exists(self.channel_network_vector_path))

    def _ensure_snapped_points(
            self,
            snap_points: Callable[[], object],
            process_channel_network: Callable[[], object],
            process_flow_accumulation: Callable[[], object],
            fill_dem: Callable[[], object]) -> bool:
        """Snap pour points if the snapped output is not available."""
        if not self.check_prerequisites(needs_pour_points=True):
            return False

        if not self._ensure_channel_network(
                process_channel_network,
                process_flow_accumulation,
                fill_dem):
            return False

        existing_path = self._restore_existing_path(
            "snapped_points_path", "2_pour_points_snapped.shp")
        if existing_path:
            if self._snapped_points_match_channel_network(
                    existing_path, self.channel_network_vector_path):
                return True
            self.log_message(
                "Existing snapped pour points are stale or incomplete. Recreating snapped points.")
            self._remove_stale_vector_output(existing_path)
            self.snapped_points_path = None

        self.log_message("Snapped pour points are missing. Running Snap Points first...")
        self.without_layer_loading(snap_points)
        return bool(self.snapped_points_path and os.path.exists(self.snapped_points_path))

    def _ensure_merged_watershed(
            self,
            delineate_watershed: Callable[[], object],
            snap_points: Callable[[], object],
            process_channel_network: Callable[[], object],
            process_flow_accumulation: Callable[[], object],
            fill_dem: Callable[[], object]) -> bool:
        """Create the merged watershed vector if it has not been generated yet."""
        if self._restore_existing_path(
                "merged_watershed_path",
                os.path.join("Watersheds", "4_watershed_merged_vector.shp"),
                "4_watershed_merged_vector.shp"):
            return True

        if not self._ensure_snapped_points(
                snap_points,
                process_channel_network,
                process_flow_accumulation,
                fill_dem):
            return False

        self.log_message("Merged watershed is missing. Running Watershed Delineation first...")
        self.without_layer_loading(delineate_watershed)
        return bool(self.merged_watershed_path and os.path.exists(self.merged_watershed_path))
