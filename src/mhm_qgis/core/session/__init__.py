# -*- coding: utf-8 -*-
"""The state one morphology run carries between its steps.

Every step needs the same handful of things -- where the project is, which CRS
was chosen, which DEM was selected, and the paths of what earlier steps produced
-- and they used to reach them through `self` on a mixin that also owned the
behaviour. Holding them here instead lets each step be a plain function.

Deliberately free of Qt and QGIS. The selected inputs are stored as **sources**
(a path or provider URI) plus their CRS, not as `QgsMapLayer` objects: the four
things the old code actually asked a layer were its source, CRS, extent and
pixel size, and anything needing a live layer can open one through
`qgis_bridge.layers`. Keeping objects out is what lets `core/morphology` be
imported and tested without QGIS installed.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any


@dataclass
class MorphologySession:
    """Everything a morphology step needs that is not one of its arguments."""

    project_folder: str = ""

    #: The CRS chosen in the dialog -- the value, not the widget.
    crs: Any = None

    #: What the user selected, as a path or provider URI. `dem_source` is
    #: replaced by the reprojected DEM once one has been produced.
    dem_source: str = ""
    pour_points_source: str = ""
    #: Source of every selected input, keyed by kind ("dem", "lai", "soil"...).
    #: Filled at the boundary from the combo boxes, so no step reads a widget.
    input_sources: dict = field(default_factory=dict)
    #: Lookup configuration per categorical kind, likewise resolved up front.
    lookup_configs: dict = field(default_factory=dict)

    # --- what earlier steps produced -------------------------------------
    filled_dem_path: str | None = None
    flow_dir_path: str | None = None
    flow_direction_path: str | None = None
    flow_accumulation_path: str | None = None
    flow_accumulation_area_path: str | None = None
    channel_network_vector_path: str | None = None
    snapped_points_path: str | None = None
    watershed_raster_path: str | None = None
    watershed_vector_path: str | None = None
    merged_watershed_path: str | None = None
    gauge_position_path: str | None = None
    aspect_path: str | None = None
    slope_path: str | None = None
    geology_path: str | None = None
    categorical_ready_outputs: dict = field(default_factory=dict)

    # --- grid headers, filled in once the L0/L2 extent is known -----------
    L0: dict | None = None
    L1: dict | None = None
    L11: dict | None = None
    L2: dict | None = None
    #: Headers for every level, as the dialog computes them.
    grid_headers: dict | None = None
    l0_resolution: Any = None
    dem_resolution_info: Any = None

    # --- the shared processing-state registry ----------------------------
    #: The morphology sections of `mhm_qgis_processing_state.json`, as loaded.
    #: `core/handlers/state/morphology` reads and writes it.
    processing_state: dict = field(
        default_factory=lambda: {"version": 1, "outputs": {}, "workflows": {}, "grid": {}})
    #: Suppress adding produced layers to the map: set for a batch run, and
    #: incremented by `without_layer_loading` around a nested call.
    skip_loading: bool = False
    loading_suppressed: int = 0

    # --- how a step talks back -------------------------------------------
    #: Running commentary for the plugin log.
    log: Callable[[str], None] | None = None
    #: `report(level, title, message)` -- see `core/utils/report` for levels.
    report: Callable[[str, str, str], None] | None = None
    #: `load(path, name, is_raster=True)` -- puts a result on the map canvas.
    load: Callable[..., Any] | None = None
    #: `run(algorithm, params)` -- dispatches one QGIS Processing algorithm.
    #: Supplied by `qgis_bridge/processing`, so no core step imports QGIS.
    run: Callable[..., Any] | None = None

    def say(self, message: str) -> None:
        """Log when there is somewhere to log to."""
        if self.log is not None:
            self.log(message)

    def tell(self, level: str, title: str, message: str) -> None:
        """Report to the user when there is anyone to report to."""
        if self.report is not None:
            self.report(level, title, message)

    def source_for(self, kind: str) -> str:
        """Return the selected source for one input kind, or ''."""
        return self.input_sources.get(kind, "") or ""

    def execute(self, algorithm: str, params: dict):
        """Run one Processing algorithm, when the caller supplied a runner."""
        if self.run is None:
            return None
        return self.run(algorithm, params)

    def show(self, path, name, is_raster: bool = True):
        """Put a produced file on the map, if the caller wants that."""
        if self.load is not None:
            return self.load(path, name, is_raster)
        return None


__all__ = ["MorphologySession"]
