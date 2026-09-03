# -*- coding: utf-8 -*-
"""Processor setup, QGIS layer loading policy, and algorithm dispatch."""
from __future__ import annotations

from collections.abc import Callable
from typing import Any

from ..common import DialogUtils, json, os
from ... import layers
from ....core.session import MorphologySession
from ....qt import reporting
from ....core.utils.report import CRITICAL, INFORMATION, WARNING
from .processing_state import ProcessingStateMixin


class BaseProcessingMixin(ProcessingStateMixin):
    """Processor setup, QGIS layer loading policy, and algorithm dispatch."""

    def __init__(self, dialog: Any) -> None:
        """
        Initialize the morphology processor with reference to the dialog.

        Args:
            dialog: The MhmQgisDialog instance (needs DialogUtils methods)
        """
        self.dialog = dialog
        # Copy DialogUtils methods to self
        self.log_message = dialog.log_message
        self.check_prerequisites = dialog.check_prerequisites
        # Note: load_layer is overridden below to support skip_loading flag
        self.get_dem_extent_and_resolution = dialog.get_dem_extent_and_resolution
        # The single seam through which anything user-facing reaches Qt. Every
        # module below used to call QMessageBox itself; now they call warn /
        # error / confirm, which become `report=` / `confirm=` parameters once
        # these mixins are plain functions.
        self._report = reporting.dialog_reporter(dialog, log=dialog.log_message)
        self._confirm = reporting.dialog_confirmer(dialog)

        # DEM layer reference (reprojected if needed)
        self.dem_layer = None  # Will store the reprojected DEM layer if CRS differs

        # Processed layer references
        # Final land use layer with reclassified values (1, 2, or 3)
        self.land_use_layer = None

        # The run state lives on a MorphologySession. The path attributes
        # are properties over it, so a step that has been converted to a
        # function taking the session and one that still says `self.x` see
        # the same values. The properties go when the last mixin does.
        self._session = MorphologySession(
            project_folder=getattr(dialog, "project_folder", "") or "",
            log=dialog.log_message,
            report=self._report,
            load=self.load_layer,
        )
        
        # Flag to skip loading layers (used in execute_all_processing)
        self.skip_loading = False
        self._layer_loading_suppressed = 0
        self.processing_state_filename = "mhm_qgis_processing_state.json"
        self.processing_state = {"version": 1, "outputs": {}}
        

    @property
    def filled_dem_path(self):
        return self._state().filled_dem_path

    @filled_dem_path.setter
    def filled_dem_path(self, value):
        self._state().filled_dem_path = value

    @property
    def flow_dir_path(self):
        return self._state().flow_dir_path

    @flow_dir_path.setter
    def flow_dir_path(self, value):
        self._state().flow_dir_path = value

    @property
    def flow_direction_path(self):
        return self._state().flow_direction_path

    @flow_direction_path.setter
    def flow_direction_path(self, value):
        self._state().flow_direction_path = value

    @property
    def flow_accumulation_path(self):
        return self._state().flow_accumulation_path

    @flow_accumulation_path.setter
    def flow_accumulation_path(self, value):
        self._state().flow_accumulation_path = value

    @property
    def flow_accumulation_area_path(self):
        return self._state().flow_accumulation_area_path

    @flow_accumulation_area_path.setter
    def flow_accumulation_area_path(self, value):
        self._state().flow_accumulation_area_path = value

    @property
    def channel_network_vector_path(self):
        return self._state().channel_network_vector_path

    @channel_network_vector_path.setter
    def channel_network_vector_path(self, value):
        self._state().channel_network_vector_path = value

    @property
    def snapped_points_path(self):
        return self._state().snapped_points_path

    @snapped_points_path.setter
    def snapped_points_path(self, value):
        self._state().snapped_points_path = value

    @property
    def watershed_raster_path(self):
        return self._state().watershed_raster_path

    @watershed_raster_path.setter
    def watershed_raster_path(self, value):
        self._state().watershed_raster_path = value

    @property
    def watershed_vector_path(self):
        return self._state().watershed_vector_path

    @watershed_vector_path.setter
    def watershed_vector_path(self, value):
        self._state().watershed_vector_path = value

    @property
    def merged_watershed_path(self):
        return self._state().merged_watershed_path

    @merged_watershed_path.setter
    def merged_watershed_path(self, value):
        self._state().merged_watershed_path = value

    @property
    def gauge_position_path(self):
        return self._state().gauge_position_path

    @gauge_position_path.setter
    def gauge_position_path(self, value):
        self._state().gauge_position_path = value

    @property
    def aspect_path(self):
        return self._state().aspect_path

    @aspect_path.setter
    def aspect_path(self, value):
        self._state().aspect_path = value

    @property
    def slope_path(self):
        return self._state().slope_path

    @slope_path.setter
    def slope_path(self, value):
        self._state().slope_path = value

    @property
    def geology_path(self):
        return self._state().geology_path

    @geology_path.setter
    def geology_path(self, value):
        self._state().geology_path = value

    @property
    def categorical_ready_outputs(self):
        return self._state().categorical_ready_outputs

    @categorical_ready_outputs.setter
    def categorical_ready_outputs(self, value):
        self._state().categorical_ready_outputs = value

    @property
    def L0(self):
        return self._state().L0

    @L0.setter
    def L0(self, value):
        self._state().L0 = value

    @property
    def L1(self):
        return self._state().L1

    @L1.setter
    def L1(self, value):
        self._state().L1 = value

    @property
    def L11(self):
        return self._state().L11

    @L11.setter
    def L11(self, value):
        self._state().L11 = value

    @property
    def L2(self):
        return self._state().L2

    @L2.setter
    def L2(self, value):
        self._state().L2 = value

    def _reporter(self):
        """Return the report callback, building one for uninitialised objects.

        Tests subclass a single mixin and never run __init__, so there may be
        no `_report` yet; a bare object with no dialog just logs.
        """
        found = self.__dict__.get("_report")
        if found is None:
            dialog = getattr(self, "dialog", None)
            log = getattr(self, "log_message", None)
            found = (reporting.dialog_reporter(dialog, log=log)
                     if dialog is not None else reporting.log_reporter(log))
            self.__dict__["_report"] = found
        return found

    def warn(self, title: str, message: str) -> None:
        """Tell the user about a recoverable problem."""
        self._reporter()(WARNING, title, message)

    def error(self, title: str, message: str) -> None:
        """Tell the user about a failure."""
        self._reporter()(CRITICAL, title, message)

    def inform(self, title: str, message: str) -> None:
        """Tell the user something succeeded."""
        self._reporter()(INFORMATION, title, message)

    def confirm(self, title: str, question: str) -> bool:
        """Ask the user a yes/no question; False when there is no one to ask."""
        return self._confirm(title, question)

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

        display_path = path
        display_name = name
        if is_raster and hasattr(self, "_preferred_display_raster_path"):
            display_path = self._preferred_display_raster_path(path)
            display_name = self._preferred_display_layer_name(name, display_path)
            if display_path != path:
                self.log_message(
                    f"Loading prepared display raster instead: {os.path.basename(display_path)}")

        self.mark_output_prepared(display_path, name=display_name, loaded=True)
        return layers.load(
            display_path, display_name, is_raster=is_raster, log=self.log_message)

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
