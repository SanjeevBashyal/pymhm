# -*- coding: utf-8 -*-
"""Processor setup, QGIS layer loading policy, and algorithm dispatch."""
from __future__ import annotations

from collections.abc import Callable
from typing import Any

from ..common import (
    json,
    processing,
    DialogUtils,
)
from .processing_state import ProcessingStateMixin


class BaseProcessingMixin(ProcessingStateMixin):
    """Processor setup, QGIS layer loading policy, and algorithm dispatch."""

    def __init__(self, dialog: Any) -> None:
        """
        Initialize the morphology processor with reference to the dialog.

        Args:
            dialog: The pymhmDialog instance (needs DialogUtils methods)
        """
        self.dialog = dialog
        # Copy DialogUtils methods to self
        self.log_message = dialog.log_message
        self.check_prerequisites = dialog.check_prerequisites
        # Note: load_layer is overridden below to support skip_loading flag
        self.get_dem_extent_and_resolution = dialog.get_dem_extent_and_resolution

        # DEM layer reference (reprojected if needed)
        self.dem_layer = None  # Will store the reprojected DEM layer if CRS differs

        # Processed layer references
        # Final land use layer with reclassified values (1, 2, or 3)
        self.land_use_layer = None

        # Paths for geometry processing outputs
        self.filled_dem_path = None
        self.flow_dir_path = None
        self.channel_network_vector_path = None
        self.snapped_points_path = None
        self.watershed_raster_path = None
        self.watershed_vector_path = None
        self.merged_watershed_path = None

        # Hydrological processing paths
        self.aspect_path = None
        self.slope_path = None
        self.flow_direction_path = None
        self.flow_accumulation_path = None
        self.flow_accumulation_area_path = None
        self.gauge_position_path = None
        
        # Layer processing paths
        self.geology_path = None
        
        # Flag to skip loading layers (used in execute_all_processing)
        self.skip_loading = False
        self._layer_loading_suppressed = 0
        self.processing_state_filename = "pymhm_processing_state.json"
        self.processing_state = {"version": 1, "outputs": {}}
        
        # Lat/Lon header information (L0, L1, L11, L2)
        self.L0 = None  # Dictionary with ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value
        self.L1 = None
        self.L11 = None
        self.L2 = None

    def load_layer(
            self,
            path: str,
            name: str,
            is_raster: bool = True) -> Any | None:
        """
        Override load_layer to respect batch/intermediate loading suppression.
        """
        if self.skip_loading or self._layer_loading_suppressed > 0:
            self.log_message(f"Prepared output without loading layer: {name}")
            self.mark_output_prepared(path, name=name, loaded=False)
            return

        self.mark_output_prepared(path, name=name, loaded=True)
        # Call the parent load_layer method
        return self.dialog.load_layer(path, name, is_raster)

    def run_processing_algorithm(
            self,
            name: str,
            params: dict[str, Any]) -> Any:
        """Run a processing algorithm and record any file outputs it creates."""
        result = self.dialog.run_processing_algorithm(name, params)
        self.record_processing_outputs(name, params, result)
        return result

    def without_layer_loading(
            self,
            callback: Callable[..., Any],
            *args: Any,
            **kwargs: Any) -> Any:
        """Run a callable while preventing QGIS layer additions."""
        self._layer_loading_suppressed += 1
        try:
            return callback(*args, **kwargs)
        finally:
            self._layer_loading_suppressed -= 1
