# -*- coding: utf-8 -*-
"""Preparation guards for prerequisite morphology outputs."""
from __future__ import annotations

from collections.abc import Callable

from ..common import os
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

    def _ensure_filled_dem(self, fill_dem: Callable[[], object]) -> bool:
        """Create the filled DEM if it has not been generated yet."""
        if self._restore_existing_path(
                "filled_dem_path", "1_dem_filled.tif", "1_dem_filled.sdat"):
            return True

        self.log_message("Filled DEM is missing. Running Fill DEM first...")
        self.without_layer_loading(fill_dem)
        return bool(self.filled_dem_path and os.path.exists(self.filled_dem_path))

    def _ensure_flow_accumulation(
            self,
            process_flow_accumulation: Callable[[], object],
            fill_dem: Callable[[], object]) -> bool:
        """Create flow accumulation if it has not been generated yet."""
        if self._restore_existing_path("flow_accumulation_path", "2_flow_accumulation.tif"):
            self._restore_existing_path(
                "flow_accumulation_area_path",
                "2_flow_accumulation_area.tif",
                "2_flow_accumulation_area.sdat"
            )
            return True

        if not self._ensure_filled_dem(fill_dem):
            return False

        self.log_message("Flow accumulation is missing. Running Flow Accumulation first...")
        self.without_layer_loading(process_flow_accumulation)
        return bool(self.flow_accumulation_path and os.path.exists(self.flow_accumulation_path))

    def _ensure_flow_direction(
            self,
            process_flow_direction: Callable[[], object],
            fill_dem: Callable[[], object]) -> bool:
        """Create flow direction if it has not been generated yet."""
        if self._restore_existing_path("flow_direction_path", "2_flow_direction.tif"):
            return True

        if not self._ensure_filled_dem(fill_dem):
            return False

        self.log_message("Flow direction is missing. Running Flow Direction first...")
        self.without_layer_loading(process_flow_direction)
        return bool(self.flow_direction_path and os.path.exists(self.flow_direction_path))

    def _ensure_channel_network(
            self,
            process_channel_network: Callable[[], object],
            process_flow_accumulation: Callable[[], object],
            fill_dem: Callable[[], object]) -> bool:
        """Create the channel network if it has not been generated yet."""
        if self._restore_existing_path("channel_network_vector_path", "2_channel_network.shp"):
            return True

        if not self._ensure_flow_accumulation(process_flow_accumulation, fill_dem):
            return False

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
        if self._restore_existing_path("snapped_points_path", "2_pour_points_snapped.shp"):
            return True

        if not self.check_prerequisites(needs_pour_points=True):
            return False

        if not self._ensure_channel_network(
                process_channel_network,
                process_flow_accumulation,
                fill_dem):
            return False

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
