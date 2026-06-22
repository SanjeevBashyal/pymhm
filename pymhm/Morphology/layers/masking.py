# -*- coding: utf-8 -*-
"""Cropping and watershed masking for prepared morphology rasters."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    QMessageBox,
)
from ..watershed.watershed_delineation import WatershedDelineationMixin
from ...grid_resolution import header_bounds, header_for_existing_bounds


class MaskingMixin(WatershedDelineationMixin):
    """Cropping and watershed masking for prepared morphology rasters."""

    def crop_all_layers(self) -> None:
        """Crop prepared morphology rasters to the common L0/L2 model extent."""
        if not self.check_prerequisites():
            return

        if not self._ensure_filled_dem(self.fill_dem):
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        try:
            l0_header = self._target_l0_header()
        except Exception as e:
            self.log_message(f"ERROR: Cannot crop morphology rasters: {e}")
            QMessageBox.warning(self.dialog, "Grid Configuration Error", str(e))
            return

        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog,
                "CRS Error",
                "Please set a valid input CRS.")
            return

        layers_to_crop = self._collect_layers_to_crop(geometry_folder)
        if not layers_to_crop:
            lai_cropped = self._crop_lai_netcdf_if_available(l0_header, input_crs)
            if lai_cropped:
                self.log_message("Crop all process completed.")
                return
            QMessageBox.warning(
                self.dialog,
                "No Layers",
                "No raster layers found to crop. Please process at least the filled DEM first.")
            return

        extent_str = self._extent_from_header(l0_header, input_crs)
        resolution = float(l0_header["cellsize"])

        self.log_message("\n--- Cropping all morphology rasters ---")
        self.log_message(f"Crop extent: {extent_str}")
        self.log_message(f"Crop resolution: {resolution:.6f}")
        self.log_message(f"Cropping {len(layers_to_crop)} raster layer(s)...")

        for layer_info in layers_to_crop:
            input_path = layer_info["input"]
            output_path = layer_info["crop"]
            masked_path = layer_info["masked"]
            layer_name = layer_info["name"]

            self.log_message(f"Cropping {layer_name}...")
            self._remove_existing_raster(output_path)
            self._remove_existing_raster(masked_path)
            params_crop = {
                "INPUT": input_path,
                "SOURCE_CRS": None,
                "TARGET_CRS": input_crs,
                "RESAMPLING": 0,
                "NODATA": -9999,
                "TARGET_RESOLUTION": resolution,
                "OPTIONS": None,
                "DATA_TYPE": 0,
                "TARGET_EXTENT": extent_str,
                "TARGET_EXTENT_CRS": input_crs,
                "MULTITHREADING": False,
                "EXTRA": "",
                "OUTPUT": output_path,
            }

            result = self.run_processing_algorithm(
                "gdal:warpreproject",
                params_crop,
            )

            if result and os.path.exists(output_path):
                self.mark_output_prepared(
                    output_path,
                    name=os.path.basename(output_path),
                    loaded=False,
                    algorithm="gdal:warpreproject",
                )
                self.log_message(f"{layer_name} cropped successfully.")
            else:
                self.log_message(f"ERROR: Failed to crop {layer_name}.")

        self._crop_lai_netcdf_if_available(l0_header, input_crs)
        self.log_message("Crop all process completed.")

    def mask_all_layers(self) -> None:
        """Mask cropped morphology rasters by merged watershed, leaving -9999 outside."""
        if not self.check_prerequisites():
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        try:
            l0_header = self._target_l0_header()
        except Exception as e:
            self.log_message(f"ERROR: Cannot mask morphology rasters: {e}")
            QMessageBox.warning(self.dialog, "Grid Configuration Error", str(e))
            return

        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog,
                "CRS Error",
                "Please set a valid input CRS.")
            return

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

        layers_to_mask = self._collect_layers_to_mask(geometry_folder)
        if not layers_to_mask:
            lai_masked = self._mask_lai_netcdf_if_available(
                l0_header,
                input_crs,
                merged_watershed_path,
            )
            if lai_masked:
                self.log_message("Masking process completed.")
                return
            QMessageBox.warning(
                self.dialog,
                "No Cropped Layers",
                "No cropped raster layers found. Please run Crop All first.")
            return

        extent_str = self._extent_from_header(l0_header, input_crs)
        resolution = float(l0_header["cellsize"])

        self.log_message("\n--- Masking all cropped morphology rasters ---")
        self.log_message(f"Mask extent: {extent_str}")
        self.log_message(f"Mask resolution: {resolution:.6f}")
        self.log_message(f"Masking {len(layers_to_mask)} raster layer(s)...")

        for layer_info in layers_to_mask:
            input_path = layer_info["input"]
            output_path = layer_info["output"]
            layer_name = layer_info["name"]

            self.log_message(f"Masking {layer_name}...")
            self._remove_existing_raster(output_path)
            params_mask = {
                "INPUT": input_path,
                "MASK": merged_watershed_path,
                "SOURCE_CRS": input_crs,
                "TARGET_CRS": input_crs,
                "NODATA": -9999,
                "ALPHA_BAND": False,
                "CROP_TO_CUTLINE": False,
                "KEEP_RESOLUTION": False,
                "SET_RESOLUTION": True,
                "X_RESOLUTION": resolution,
                "Y_RESOLUTION": resolution,
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
                self.mark_output_prepared(
                    output_path,
                    name=os.path.basename(output_path),
                    loaded=False,
                    algorithm="gdal:cliprasterbymasklayer",
                )
                self.log_message(f"{layer_name} masked successfully.")
            else:
                self.log_message(f"ERROR: Failed to mask {layer_name}.")

        self._mask_lai_netcdf_if_available(
            l0_header,
            input_crs,
            merged_watershed_path,
        )
        self.log_message("Masking process completed.")

    def _target_l0_header(self):
        """Return the L0 grid header derived from the configured model extent."""
        l2_header = getattr(self.dialog, "_grid_l2_header", None)
        if not l2_header and hasattr(self.dialog, "update_l2_resolution_from_metadata"):
            self.dialog.update_l2_resolution_from_metadata()
            l2_header = getattr(self.dialog, "_grid_l2_header", None)
        if not l2_header:
            raise ValueError("L2 grid is not available. Run and save meteorology data first.")

        l0_resolution = None
        if hasattr(self.dialog, "current_l0_resolution"):
            l0_resolution = self.dialog.current_l0_resolution()
        if not l0_resolution:
            l0_info = getattr(self.dialog, "_grid_l0_info", None)
            if not l0_info and hasattr(self.dialog, "filled_dem_resolution_info"):
                l0_info = self.dialog.filled_dem_resolution_info()
            if l0_info:
                l0_resolution = float(l0_info["resolution"])
        if not l0_resolution:
            raise ValueError("L0 resolution is not available. Select or prepare the filled DEM first.")

        return header_for_existing_bounds(l2_header, float(l0_resolution))

    def _extent_from_header(self, header, input_crs):
        """Return a QGIS Processing extent string from a grid header."""
        xmin, xmax, ymin, ymax = header_bounds(header)
        return f"{xmin},{xmax},{ymin},{ymax} [{input_crs.authid()}]"

    def _crop_lai_netcdf_if_available(self, l0_header, input_crs):
        """Crop LAI NetCDF data when the long-term monthly LAI input is selected."""
        crop_lai = getattr(self, "crop_lai_netcdf_to_l0", None)
        if crop_lai is None:
            return False
        return bool(crop_lai(l0_header, input_crs))

    def _mask_lai_netcdf_if_available(
            self,
            l0_header,
            input_crs,
            merged_watershed_path):
        """Mask LAI NetCDF data when the long-term monthly LAI input is selected."""
        mask_lai = getattr(self, "mask_lai_netcdf_to_l0", None)
        if mask_lai is None:
            return False
        return bool(mask_lai(l0_header, input_crs, merged_watershed_path))

    def _collect_layers_to_crop(self, geometry_folder):
        """Return prepared base morphology rasters that should be cropped."""
        layer_specs = [
            ("filled_dem_path", "1_dem_filled.tif", "1_DEM_Filled"),
            ("aspect_path", "1_dem_aspect.tif", "1_DEM_Aspect"),
            ("slope_path", "1_dem_slope.tif", "1_DEM_Slope"),
            ("flow_accumulation_path", "2_flow_accumulation.tif", "2_Flow_Accumulation"),
            ("flow_direction_path", "2_flow_direction.tif", "2_Flow_Direction"),
            ("gauge_position_path", "2_gauge_position.tif", "2_Gauge_Position"),
            ("land_use_layer", "3_land_use.tif", "3_Land_Use"),
            (None, "3_soil.tif", "3_Soil"),
            ("geology_path", "3_geology_processed.tif", "3_Geology"),
        ]

        layers_to_crop = []
        seen = set()
        for attr_name, filename, layer_name in layer_specs:
            expected_path = os.path.join(geometry_folder, filename)
            layer_path = getattr(self, attr_name, None) if attr_name else None
            if not layer_path or not os.path.exists(layer_path):
                layer_path = expected_path
            if not layer_path or not os.path.exists(layer_path):
                continue

            normalized = os.path.abspath(layer_path)
            if normalized in seen:
                continue
            seen.add(normalized)

            layers_to_crop.append({
                "input": layer_path,
                "crop": self._raster_variant_path(layer_path, "_crop"),
                "masked": self._raster_variant_path(layer_path, "_masked"),
                "name": layer_name,
            })
        return layers_to_crop

    def _collect_layers_to_mask(self, geometry_folder):
        """Return cropped raster layers that should be watershed-masked."""
        layers_to_mask = []
        for layer_info in self._collect_layers_to_crop(geometry_folder):
            crop_path = layer_info["crop"]
            if not crop_path or not os.path.exists(crop_path):
                continue
            layers_to_mask.append({
                "input": crop_path,
                "output": layer_info["masked"],
                "name": f"{layer_info['name']}_Masked",
            })
        return layers_to_mask

    def _raster_variant_path(self, raster_path, suffix):
        """Return a sidecar raster path with the requested suffix."""
        folder, filename = os.path.split(raster_path)
        stem, ext = os.path.splitext(filename)
        for existing_suffix in ("_crop", "_masked"):
            if stem.endswith(existing_suffix):
                stem = stem[: -len(existing_suffix)]
        return os.path.join(folder, f"{stem}{suffix}{ext or '.tif'}")

    def _remove_existing_raster(self, raster_path):
        """Remove an old raster output before overwriting it."""
        if not raster_path or not os.path.exists(raster_path):
            return
        try:
            os.remove(raster_path)
        except Exception as e:
            self.log_message(f"WARNING: Could not remove existing raster {raster_path}: {e}")

    def _preferred_display_raster_path(self, raster_path):
        """Return masked/cropped raster variant if it exists for inspection."""
        if not raster_path:
            return raster_path
        stem = os.path.splitext(os.path.basename(raster_path))[0]
        if stem.endswith("_masked") or stem.endswith("_crop"):
            return raster_path
        masked_path = self._raster_variant_path(raster_path, "_masked")
        if os.path.exists(masked_path):
            return masked_path
        crop_path = self._raster_variant_path(raster_path, "_crop")
        if os.path.exists(crop_path):
            return crop_path
        return raster_path

    def _preferred_display_layer_name(self, layer_name, raster_path):
        """Return a layer name matching the preferred display raster variant."""
        stem = os.path.splitext(os.path.basename(raster_path or ""))[0]
        if stem.endswith("_masked") and not layer_name.endswith("_Masked"):
            return f"{layer_name}_Masked"
        if stem.endswith("_crop") and not layer_name.endswith("_Crop"):
            return f"{layer_name}_Crop"
        return layer_name
