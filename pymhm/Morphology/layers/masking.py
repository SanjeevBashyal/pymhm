# -*- coding: utf-8 -*-
"""Watershed mask application for prepared morphology rasters."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    QMessageBox,
    QgsRasterLayer,
)
from ..watershed.watershed_delineation import WatershedDelineationMixin
from ...grid_resolution import header_bounds, raster_resolution_info


class MaskingMixin(WatershedDelineationMixin):
    """Watershed mask application for prepared morphology rasters."""

    def mask_all_layers(self) -> None:
        """Mask all prepared rasters using derived L0 resolution and L2 extent."""
        if not self.check_prerequisites():
            return

        if not self._ensure_filled_dem(self.fill_dem):
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        merged_watershed_path = self._restore_existing_path(
            "merged_watershed_path",
            os.path.join("Watersheds", "4_watershed_merged_vector.shp"),
            "4_watershed_merged_vector.shp",
        )

        if not merged_watershed_path or not os.path.exists(merged_watershed_path):
            if not self._ensure_merged_watershed(
                    self.delineate_watershed,
                    self.snap_points,
                    self.process_channel_network,
                    self.process_flow_accumulation,
                    self.fill_dem):
                return
            merged_watershed_path = self.merged_watershed_path

        if not merged_watershed_path or not os.path.exists(merged_watershed_path):
            return

        self.log_message("\n--- Masking all layers ---")

        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog,
                "Error",
                "Cannot read filled DEM layer.")
            return

        dem_extent = filled_dem_layer.extent()
        dem_crs = filled_dem_layer.crs()
        input_crs = self.dialog.get_crs()

        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog,
                "CRS Error",
                "Please set a valid input CRS.")
            return

        self.log_message(f"Filled DEM CRS: {dem_crs.authid()}")
        self.log_message(
            f"Original DEM extent: {dem_extent.xMinimum():.2f}, "
            f"{dem_extent.xMaximum():.2f}, {dem_extent.yMinimum():.2f}, "
            f"{dem_extent.yMaximum():.2f}")

        x_resolution, y_resolution = self._masking_l0_resolution(
            filled_dem_layer,
            dem_extent,
        )
        extent_str = self._masking_target_extent(dem_extent, input_crs)

        layers_to_mask = self._collect_layers_to_mask(geometry_folder)
        if not layers_to_mask:
            QMessageBox.warning(
                self.dialog,
                "No Layers",
                "No raster layers found to mask. Please process at least the filled DEM first.")
            return

        self.log_message(f"Masking {len(layers_to_mask)} raster layer(s)...")
        for layer_info in layers_to_mask:
            input_path = layer_info["input"]
            output_path = layer_info["output"]
            layer_name = layer_info["name"]
            load_output = layer_info.get("load", True)

            self.log_message(f"Masking {layer_name}...")
            params_mask = {
                "INPUT": input_path,
                "MASK": merged_watershed_path,
                "SOURCE_CRS": input_crs,
                "TARGET_CRS": input_crs,
                "NODATA": None,
                "ALPHA_BAND": False,
                "CROP_TO_CUTLINE": False,
                "KEEP_RESOLUTION": False,
                "SET_RESOLUTION": True,
                "X_RESOLUTION": x_resolution,
                "Y_RESOLUTION": y_resolution,
                "MULTITHREADING": False,
                "OPTIONS": None,
                "DATA_TYPE": 0,
                "EXTRA": "",
                "OUTPUT": output_path,
            }
            if extent_str:
                params_mask["TARGET_EXTENT"] = extent_str

            result = self.run_processing_algorithm(
                "gdal:cliprasterbymasklayer",
                params_mask,
            )

            if result and os.path.exists(output_path):
                if load_output:
                    self.load_layer(output_path, layer_name)
                self.log_message(f"{layer_name} masked successfully.")
            else:
                self.log_message(f"ERROR: Failed to mask {layer_name}.")

        self.log_message("Masking process completed.")

    def _masking_l0_resolution(self, filled_dem_layer, dem_extent):
        """Return L0 x/y resolution for masking."""
        l0_info = getattr(self.dialog, "_grid_l0_info", None)
        if not l0_info:
            l0_info = raster_resolution_info(filled_dem_layer)

        if l0_info:
            resolution = float(l0_info["resolution"])
            self.log_message(
                f"L0 resolution: {resolution:.6f} {l0_info.get('unit', '')}")
            return resolution, resolution

        x_resolution = abs(
            (dem_extent.xMaximum() - dem_extent.xMinimum()) / filled_dem_layer.width())
        y_resolution = abs(
            (dem_extent.yMaximum() - dem_extent.yMinimum()) / filled_dem_layer.height())
        self.log_message(
            f"L0 resolution from filled DEM: x={x_resolution:.6f}, y={y_resolution:.6f}")
        return x_resolution, y_resolution

    def _masking_target_extent(self, dem_extent, input_crs):
        """Return target extent string from prepared L2 grid or filled DEM."""
        l2_header = getattr(self.dialog, "_grid_l2_header", None)
        if l2_header:
            xmin, xmax, ymin, ymax = header_bounds(l2_header)
            extent_str = f"{xmin},{xmax},{ymin},{ymax} [{input_crs.authid()}]"
            self.log_message(f"Masking extent from prepared L2 grid: {extent_str}")
            return extent_str

        extent_str = (
            f"{dem_extent.xMinimum()},{dem_extent.xMaximum()},"
            f"{dem_extent.yMinimum()},{dem_extent.yMaximum()} "
            f"[{input_crs.authid()}]"
        )
        self.log_message(f"Masking extent from filled DEM: {extent_str}")
        return extent_str

    def _collect_layers_to_mask(self, geometry_folder):
        """Return prepared raster layers that should be masked."""
        layers_to_mask = []
        if self.filled_dem_path and os.path.exists(self.filled_dem_path):
            layers_to_mask.append({
                "input": self.filled_dem_path,
                "output": os.path.join(geometry_folder, "1_dem_filled_masked.tif"),
                "name": "1_DEM_Filled_Masked",
                "load": False,
            })

        layer_mapping = [
            ("aspect_path", "1_dem_aspect_masked.tif", "1_DEM_Aspect_Masked"),
            ("slope_path", "1_dem_slope_masked.tif", "1_DEM_Slope_Masked"),
            ("flow_accumulation_path", "2_flow_accumulation_masked.tif", "2_Flow_Accumulation_Masked"),
            ("flow_direction_path", "2_flow_direction_masked.tif", "2_Flow_Direction_Masked"),
            ("gauge_position_path", "2_gauge_position_masked.tif", "2_Gauge_Position_Masked"),
        ]
        for attr_name, output_filename, layer_name in layer_mapping:
            layer_path = getattr(self, attr_name, None)
            if layer_path and os.path.exists(layer_path):
                layers_to_mask.append({
                    "input": layer_path,
                    "output": os.path.join(geometry_folder, output_filename),
                    "name": layer_name,
                })

        if self.land_use_layer and os.path.exists(self.land_use_layer):
            layers_to_mask.append({
                "input": self.land_use_layer,
                "output": os.path.join(geometry_folder, "3_land_use_masked.tif"),
                "name": "3_Land_Use_Masked",
            })

        soil_path = os.path.join(geometry_folder, "3_soil.tif")
        if os.path.exists(soil_path):
            layers_to_mask.append({
                "input": soil_path,
                "output": os.path.join(geometry_folder, "3_soil_masked.tif"),
                "name": "3_Soil_Masked",
            })

        if self.geology_path and os.path.exists(self.geology_path):
            layers_to_mask.append({
                "input": self.geology_path,
                "output": os.path.join(geometry_folder, "3_geology_processed_masked.tif"),
                "name": "3_Geology_Masked",
            })

        return layers_to_mask
