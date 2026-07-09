# -*- coding: utf-8 -*-
"""Flow direction and flow accumulation outputs."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    processing,
)
from ..core.predecessors import PredecessorMixin
from ..watershed.dem_fill import DemFillMixin


class FlowAnalysisMixin(DemFillMixin, PredecessorMixin):
    """Flow direction and flow accumulation outputs."""

    def process_flow_direction(self) -> None:
        """Process D8 flow direction with pyflwdir."""
        if not self._ensure_filled_dem(self.fill_dem):
            return

        # Check if flow direction already exists
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        self.flow_direction_path = os.path.join(geometry_folder, "2_flow_direction.tif")
        if self.flow_direction_path and os.path.exists(self.flow_direction_path):
            if self._raster_output_matches_filled_dem(self.flow_direction_path):
                self.log_message("Flow Direction already exists. Loading existing file...")
                self.load_layer(self.flow_direction_path, "2_Flow_Direction")
                return
            self.log_message(
                "Existing Flow Direction does not match the current filled DEM. Recreating it.")
            self._remove_stale_raster_output(self.flow_direction_path)

        self.log_message("Processing D8 Flow Direction with pyflwdir...")

        context = self._build_flwdir_from_filled_dem()
        if not context:
            self.log_message("Flow Direction processing failed.")
            return

        deps = context["deps"]
        np = deps["np"]
        gdal = deps["gdal"]
        flow_direction = context["flwdir"].to_array().astype(np.uint8)
        flow_direction[context["invalid_mask"]] = 247

        if self._write_raster_array(
            self.flow_direction_path,
            flow_direction,
            context["reference"],
            nodata=247,
            gdal_type=gdal.GDT_Byte,
        ):
            self.load_layer(self.flow_direction_path, "2_Flow_Direction")
            self.log_message("Flow Direction processing completed successfully.")
        else:
            self.log_message("Flow Direction processing failed.")

    def process_flow_accumulation(self) -> None:
        """Process cell-based flow accumulation with pyflwdir."""
        if not self._ensure_filled_dem(self.fill_dem):
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        self.flow_accumulation_area_path = os.path.join(
            geometry_folder, "2_flow_accumulation_area.tif"
        )
        self.flow_accumulation_path = os.path.join(
            geometry_folder, "2_flow_accumulation.tif"
        )

        # Check if final flow accumulation already exists
        if self.flow_accumulation_path and os.path.exists(self.flow_accumulation_path):
            if self._raster_output_matches_filled_dem(self.flow_accumulation_path):
                self.log_message(
                    "Flow Accumulation (pixels, integer) already exists. Loading existing file..."
                )
                self.load_layer(self.flow_accumulation_path, "2_Flow_Accumulation")
                return
            self.log_message(
                "Existing Flow Accumulation does not match the current filled DEM. Recreating it.")
            self._remove_stale_raster_output(self.flow_accumulation_path)
            self._remove_stale_raster_output(self.flow_accumulation_area_path)

        self.log_message("Processing Flow Accumulation with pyflwdir...")

        context = self._build_flwdir_from_filled_dem()
        if not context:
            self.log_message("ERROR: Flow Accumulation processing failed.")
            return

        deps = context["deps"]
        np = deps["np"]
        gdal = deps["gdal"]

        flow_accumulation = context["flwdir"].upstream_area(unit="cell")
        flow_accumulation = np.rint(flow_accumulation).astype(np.int32)
        flow_accumulation[context["invalid_mask"]] = -9999

        flow_accumulation_area = context["flwdir"].upstream_area(unit="m2")
        flow_accumulation_area = flow_accumulation_area.astype(np.float32)
        flow_accumulation_area[context["invalid_mask"]] = -9999.0

        wrote_cells = self._write_raster_array(
            self.flow_accumulation_path,
            flow_accumulation,
            context["reference"],
            nodata=-9999,
            gdal_type=gdal.GDT_Int32,
        )
        wrote_area = self._write_raster_array(
            self.flow_accumulation_area_path,
            flow_accumulation_area,
            context["reference"],
            nodata=-9999.0,
            gdal_type=gdal.GDT_Float32,
        )

        if wrote_cells:
            self.load_layer(self.flow_accumulation_path, "2_Flow_Accumulation")
            self.log_message(
                "Flow Accumulation (pixels, integer) processing completed successfully."
            )
            if wrote_area:
                self.log_message(
                    "Flow Accumulation area raster written for channel thresholding."
                )
            else:
                self.log_message(
                    "WARNING: Flow Accumulation area raster could not be written."
                )
        else:
            self.log_message("ERROR: Failed to write flow accumulation raster.")
