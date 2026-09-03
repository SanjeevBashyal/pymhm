# -*- coding: utf-8 -*-
"""Adapter: terrain derivatives, now computed in `core/morphology`.

The two methods below resolve the selected DEM into the session and hand off.
They exist because `qt/bindings` calls `processor.process_aspect` by name; they
go when the processor does.
"""
from __future__ import annotations

from ....core.morphology.hydrology import terrain
from ....core.utils.report import WARNING
from ..core.dem_inputs import DemInputMixin


class TerrainAnalysisMixin(DemInputMixin):
    """Terrain derivative outputs: aspect and slope."""

    def _dem_into_session(self) -> bool:
        """Put the reprojected DEM's source and CRS on the session."""
        if not self.dem_layer:
            selected = self.dialog.input_combo("dem").currentLayer()
            if not selected:
                self.warn("Input Error", "Please select an input DEM layer.")
                return False
            self.check_and_reproject_dem(selected)

        layer = self.dem_layer or self.dialog.input_combo("dem").currentLayer()
        if layer is None:
            self.warn("Input Error", "Please select an input DEM layer.")
            return False
        session = self.session
        session.dem_source = layer.source()
        session.crs = layer.crs()
        return True

    def process_aspect(self) -> None:
        """Process Aspect from the input DEM."""
        if self._dem_into_session():
            terrain.produce(self.session, terrain.ASPECT,
                            self.run_processing_algorithm)

    def process_slope(self) -> None:
        """Process Slope from the input DEM."""
        if self._dem_into_session():
            terrain.produce(self.session, terrain.SLOPE,
                            self.run_processing_algorithm)


__all__ = ["TerrainAnalysisMixin"]
