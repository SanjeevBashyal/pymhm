# -*- coding: utf-8 -*-
"""Watershed mask application for prepared morphology rasters."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    math,
    QMessageBox,
    QgsRasterLayer,
    QgsGeometry,
    QgsProject,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsRectangle,
    QgsPointXY,
)
from ..watershed.watershed_delineation import WatershedDelineationMixin


class MaskingMixin(WatershedDelineationMixin):
    """Watershed mask application for prepared morphology rasters."""

    def mask_all_layers(self) -> None:
        """
        Mask all raster layers using the merged watershed vector.
        Reads L0, L1, L2 values from UI and calculates extent based on L2.
        """
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        if not self._ensure_filled_dem(self.fill_dem):
            return
        
        # Check for merged watershed mask
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        merged_watershed_path = self._restore_existing_path(
            "merged_watershed_path",
            os.path.join("Watersheds", "4_watershed_merged_vector.shp"),
            "4_watershed_merged_vector.shp"
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
        
        # Read L0, L1, L2 values and units from UI
        try:
            l0_value = float(self.dialog.lineEdit_L0.text())
            l0_unit = self.dialog.comboBox_L0.currentText()
            l1_value = float(self.dialog.lineEdit_L1.text())
            l1_unit = self.dialog.comboBox_L1.currentText()
            l2_value = float(self.dialog.lineEdit_L2.text())
            l2_unit = self.dialog.comboBox_L2.currentText()
        except (ValueError, AttributeError) as e:
            QMessageBox.warning(
                self.dialog, "Input Error",
                f"Please enter valid numeric values for L0, L1, and L2.\nError: {e}")
            return
        
        if not l0_unit or not l1_unit or not l2_unit:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select units (m or °) for L0, L1, and L2.")
            return
        
        self.log_message(f"\n--- Masking all layers ---")
        self.log_message(f"L0: {l0_value} {l0_unit}")
        self.log_message(f"L1: {l1_value} {l1_unit}")
        self.log_message(f"L2: {l2_value} {l2_unit}")
        
        # Get filled DEM layer
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM layer.")
            return
        
        dem_extent = filled_dem_layer.extent()
        dem_crs = filled_dem_layer.crs()
        input_crs = self.dialog.get_crs()
        
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog, "CRS Error",
                "Please set a valid input CRS.")
            return
        
        self.log_message(f"Filled DEM CRS: {dem_crs.authid()}")
        self.log_message(
            f"Original DEM extent: {dem_extent.xMinimum():.2f}, "
            f"{dem_extent.xMaximum():.2f}, {dem_extent.yMinimum():.2f}, "
            f"{dem_extent.yMaximum():.2f}")
        
        # Calculate extent based on L2
        # For ERA5 grid snapping: if L2 is in degrees, do all operations in WGS84 degrees
        if l2_unit == "°":
            # Transform extent to WGS84 to snap to ERA5 grid
            wgs84_crs = QgsCoordinateReferenceSystem('EPSG:4326')
            
            # Transform extent to WGS84
            transform_to_wgs84 = QgsCoordinateTransform(
                dem_crs, wgs84_crs, QgsProject.instance())
            transform_to_wgs84.setBallparkTransformsAreAppropriate(True)
            dem_extent_wgs84 = transform_to_wgs84.transform(dem_extent)
            
            self.log_message(
                f"DEM extent in WGS84: {dem_extent_wgs84.xMinimum():.6f}, "
                f"{dem_extent_wgs84.xMaximum():.6f}, {dem_extent_wgs84.yMinimum():.6f}, "
                f"{dem_extent_wgs84.yMaximum():.6f}")
            
            # Round to L2 grid in degrees (snap to ERA5 grid)
            xmin_rounded = math.floor(dem_extent_wgs84.xMinimum() / l2_value) * l2_value
            xmax_rounded = math.ceil(dem_extent_wgs84.xMaximum() / l2_value) * l2_value
            ymin_rounded = math.floor(dem_extent_wgs84.yMinimum() / l2_value) * l2_value
            ymax_rounded = math.ceil(dem_extent_wgs84.yMaximum() / l2_value) * l2_value
            
            # Shift L2/2 inside on all 4 directions (in degrees)
            half_l2 = l2_value / 2.0
            xmin_shifted = xmin_rounded + half_l2
            xmax_shifted = xmax_rounded - half_l2
            ymin_shifted = ymin_rounded + half_l2
            ymax_shifted = ymax_rounded - half_l2
            
            # Validate shifted extent
            if xmin_shifted >= xmax_shifted or ymin_shifted >= ymax_shifted:
                QMessageBox.warning(
                    self.dialog, "Extent Error",
                    f"Invalid extent after shifting. L2 ({l2_value}°) may be too large relative to DEM extent.")
                return
            
            # Create shifted extent rectangle in WGS84
            shifted_extent_wgs84 = QgsRectangle(
                xmin_shifted, ymin_shifted, xmax_shifted, ymax_shifted)
            
            self.log_message(
                f"Extent rounded to L2 grid ({l2_value}°) and shifted L2/2 in WGS84: "
                f"{xmin_shifted:.6f}, {xmax_shifted:.6f}, {ymin_shifted:.6f}, {ymax_shifted:.6f}")
            
            # Transform back to input CRS
            transform_back = QgsCoordinateTransform(
                wgs84_crs, input_crs, QgsProject.instance())
            transform_back.setBallparkTransformsAreAppropriate(True)
            
            # Transform the four corners individually to ensure valid transformation
            corners = [
                QgsPointXY(xmin_shifted, ymin_shifted),  # Bottom-left
                QgsPointXY(xmax_shifted, ymin_shifted),  # Bottom-right
                QgsPointXY(xmax_shifted, ymax_shifted),  # Top-right
                QgsPointXY(xmin_shifted, ymax_shifted)   # Top-left
            ]
            
            transformed_corners = []
            for corner in corners:
                try:
                    transformed_corner = transform_back.transform(corner)
                    transformed_corners.append(transformed_corner)
                except Exception as e:
                    self.log_message(f"WARNING: Corner transformation failed: {e}")
                    # Fallback: use geometry transform
                    geom = QgsGeometry.fromPointXY(corner)
                    geom.transform(transform_back)
                    transformed_corners.append(geom.asPoint())
            
            # Create extent from transformed corners
            x_coords = [p.x() for p in transformed_corners]
            y_coords = [p.y() for p in transformed_corners]
            xmin_final = min(x_coords)
            xmax_final = max(x_coords)
            ymin_final = min(y_coords)
            ymax_final = max(y_coords)
            
            extent_str = f"{xmin_final},{xmax_final},{ymin_final},{ymax_final} [{input_crs.authid()}]"
            
            self.log_message(
                f"Extent transformed back to input CRS ({input_crs.authid()}): {extent_str}")
        else:
            # L2 is in meters - round directly in input CRS
            # Round to L2 value
            xmin_rounded = math.floor(dem_extent.xMinimum() / l2_value) * l2_value
            xmax_rounded = math.ceil(dem_extent.xMaximum() / l2_value) * l2_value
            ymin_rounded = math.floor(dem_extent.yMinimum() / l2_value) * l2_value
            ymax_rounded = math.ceil(dem_extent.yMaximum() / l2_value) * l2_value
            
            # Shift L2/2 inside on all 4 directions
            half_l2 = l2_value / 2.0
            xmin_shifted = xmin_rounded + half_l2
            xmax_shifted = xmax_rounded - half_l2
            ymin_shifted = ymin_rounded + half_l2
            ymax_shifted = ymax_rounded - half_l2
            
            # Validate shifted extent
            if xmin_shifted >= xmax_shifted or ymin_shifted >= ymax_shifted:
                QMessageBox.warning(
                    self.dialog, "Extent Error",
                    f"Invalid extent after shifting. L2 ({l2_value} m) may be too large relative to DEM extent.")
                return
            
            extent_str = f"{xmin_shifted},{xmax_shifted},{ymin_shifted},{ymax_shifted} [{input_crs.authid()}]"
            
            self.log_message(
                f"Extent rounded to L2 ({l2_value} m) and shifted L2/2 "
                f"({half_l2:.2f} m) inside in input CRS: {extent_str}")
        
        # Convert L0 to meters if in degrees
        l0_meters = l0_value
        if l0_unit == "°":
            # Convert degrees to meters (approximate: 1 degree ≈ 111,000 meters)
            # For more accuracy, could use cos(lat) for longitude, but using average
            center_lat = (dem_extent.yMinimum() + dem_extent.yMaximum()) / 2.0
            l0_meters = l0_value * 111000
            self.log_message(
                f"L0 converted from degrees to meters: {l0_meters:.2f} m "
                f"(using center latitude: {center_lat:.4f})")
        
        x_resolution = l0_meters
        y_resolution = l0_meters
        
        self.log_message(f"X Resolution (L0): {x_resolution:.2f} m")
        self.log_message(f"Y Resolution (L0): {y_resolution:.2f} m")
        
        # Get list of raster layers to mask
        layers_to_mask = []
        
        # Check for filled DEM
        if self.filled_dem_path and os.path.exists(self.filled_dem_path):
            layers_to_mask.append({
                'input': self.filled_dem_path,
                'output': os.path.join(geometry_folder, "1_dem_filled_masked.tif"),
                'name': "1_DEM_Filled_Masked",
                'load': False
            })
        
        # Check for other processed layers
        layer_mapping = [
            ('aspect_path', '1_dem_aspect_masked.tif', '1_DEM_Aspect_Masked'),
            ('slope_path', '1_dem_slope_masked.tif', '1_DEM_Slope_Masked'),
            ('flow_accumulation_path', '2_flow_accumulation_masked.tif', '2_Flow_Accumulation_Masked'),
            ('flow_direction_path', '2_flow_direction_masked.tif', '2_Flow_Direction_Masked'),
            ('gauge_position_path', '2_gauge_position_masked.tif', '2_Gauge_Position_Masked'),
        ]
        
        for attr_name, output_filename, layer_name in layer_mapping:
            layer_path = getattr(self, attr_name, None)
            if layer_path and os.path.exists(layer_path):
                layers_to_mask.append({
                    'input': layer_path,
                    'output': os.path.join(geometry_folder, output_filename),
                    'name': layer_name
                })
        
        # Check for land cover layer
        if self.land_use_layer and os.path.exists(self.land_use_layer):
            layers_to_mask.append({
                'input': self.land_use_layer,
                'output': os.path.join(geometry_folder, "3_land_use_masked.tif"),
                'name': "3_Land_Use_Masked"
            })
        
        # Check for soil layer
        soil_path = os.path.join(geometry_folder, "3_soil.tif")
        if os.path.exists(soil_path):
            layers_to_mask.append({
                'input': soil_path,
                'output': os.path.join(geometry_folder, "3_soil_masked.tif"),
                'name': "3_Soil_Masked"
            })
        
        # Check for geology layer
        if self.geology_path and os.path.exists(self.geology_path):
            layers_to_mask.append({
                'input': self.geology_path,
                'output': os.path.join(geometry_folder, "3_geology_processed_masked.tif"),
                'name': "3_Geology_Masked"
            })
        
        if not layers_to_mask:
            QMessageBox.warning(
                self.dialog, "No Layers",
                "No raster layers found to mask. Please process at least the filled DEM first.")
            return
        
        self.log_message(f"Masking {len(layers_to_mask)} raster layer(s)...")
        
        # Mask each layer
        for layer_info in layers_to_mask:
            input_path = layer_info['input']
            output_path = layer_info['output']
            layer_name = layer_info['name']
            load_output = layer_info.get('load', True)
            
            self.log_message(f"Masking {layer_name}...")
            
            # Build params - only include TARGET_EXTENT if we have a valid extent string
            params_mask = {
                'INPUT': input_path,
                'MASK': merged_watershed_path,
                'SOURCE_CRS': input_crs,
                'TARGET_CRS': input_crs,
                'NODATA': None,
                'ALPHA_BAND': False,
                'CROP_TO_CUTLINE': False,
                'KEEP_RESOLUTION': False,
                'SET_RESOLUTION': True,
                'X_RESOLUTION': x_resolution,
                'Y_RESOLUTION': y_resolution,
                'MULTITHREADING': False,
                'OPTIONS': None,
                'DATA_TYPE': 0,
                'EXTRA': '',
                'OUTPUT': output_path
            }
            
            # Only add TARGET_EXTENT if we have a valid extent string
            # GDAL will use TARGET_CRS for the extent automatically
            if extent_str:
                params_mask['TARGET_EXTENT'] = extent_str
            
            result = self.run_processing_algorithm(
                "gdal:cliprasterbymasklayer", params_mask)
            
            if result and os.path.exists(output_path):
                if load_output:
                    self.load_layer(output_path, layer_name)
                self.log_message(f"{layer_name} masked successfully.")
            else:
                self.log_message(f"ERROR: Failed to mask {layer_name}.")
        
        self.log_message("Masking process completed.")
