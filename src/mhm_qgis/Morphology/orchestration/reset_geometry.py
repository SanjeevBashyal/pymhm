# -*- coding: utf-8 -*-
"""Geometry output cleanup and in-memory path reset."""
from __future__ import annotations

from ..common import (QgsProject, QMessageBox, os, processing,
                      project_geometry_folder)
from ..core.base import BaseProcessingMixin


class ResetGeometryMixin(BaseProcessingMixin):
    """Geometry output cleanup and in-memory path reset."""

    def resetGeometry(self) -> None:
        """
        Reset geometry processing by:
        1. Removing all layers from QGIS interface that were created by the plugin
        2. Deleting all files in the geometry folder
        3. Resetting all path attributes
        """
        from qgis.PyQt.QtWidgets import QMessageBox

        # Confirm with user
        reply = QMessageBox.question(
            self.dialog, "Reset Geometry",
            "This will remove all geometry layers from QGIS and delete all files in the Geometry folder.\n\n"
            "Are you sure you want to continue?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No)
        
        if reply != QMessageBox.Yes:
            self.log_message("Reset cancelled by user.")
            return
        
        self.log_message("\n--- Resetting Geometry Processing ---")
        
        # Step 1: Remove layers from QGIS interface
        project = QgsProject.instance()
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        
        if not os.path.exists(geometry_folder):
            self.log_message("Geometry folder does not exist. Nothing to reset.")
            os.makedirs(geometry_folder, exist_ok=True)
        
        # Get all layers and identify those from geometry folder
        layers_to_remove = []
        for layer_id, layer in project.mapLayers().items():
            source = layer.source()
            # Check if layer source is in geometry folder
            if geometry_folder in source.replace('\\', '/'):
                layers_to_remove.append(layer_id)
                self.log_message(f"Removing layer from QGIS: {layer.name()}")
        
        # Remove layers
        for layer_id in layers_to_remove:
            project.removeMapLayer(layer_id)
        
        if layers_to_remove:
            self.log_message(f"Removed {len(layers_to_remove)} layer(s) from QGIS interface.")
        else:
            self.log_message("No layers found in QGIS interface to remove.")
        
        # Step 2: Delete all files in geometry folder
        deleted_count = 0
        try:
            for root, dirs, files in os.walk(geometry_folder):
                for file in files:
                    file_path = os.path.join(root, file)
                    try:
                        os.remove(file_path)
                        deleted_count += 1
                        self.log_message(f"Deleted: {os.path.basename(file_path)}")
                    except Exception as e:
                        self.log_message(f"WARNING: Could not delete {file_path}: {e}")
            
            # Also try to remove empty subdirectories (optional, but clean)
            for root, dirs, files in os.walk(geometry_folder, topdown=False):
                for dir_name in dirs:
                    dir_path = os.path.join(root, dir_name)
                    try:
                        os.rmdir(dir_path)
                        self.log_message(f"Removed empty directory: {dir_name}")
                    except:
                        pass  # Directory not empty or error, skip
        
        except Exception as e:
            self.log_message(f"ERROR while deleting files: {e}")
            QMessageBox.warning(
                self.dialog, "Reset Error",
                f"Some files could not be deleted. Check the log for details.\nError: {e}")
        
        if deleted_count > 0:
            self.log_message(f"Deleted {deleted_count} file(s) from Geometry folder.")
        else:
            self.log_message("No files found to delete in Geometry folder.")

        self._invalidate_domain_delineation_outputs()
        
        # Step 3: Reset all path attributes
        self.filled_dem_path = None
        self.flow_dir_path = None
        self.channel_network_vector_path = None
        self.snapped_points_path = None
        self.watershed_raster_path = None
        self.watershed_vector_path = None
        self.merged_watershed_path = None
        self.aspect_path = None
        self.slope_path = None
        self.flow_direction_path = None
        self.flow_accumulation_path = None
        self.flow_accumulation_area_path = None
        self.gauge_position_path = None
        self.geology_path = None
        self.land_use_layer = None
        self.categorical_ready_outputs = {}
        self.dem_layer = None
        self.processing_state = {"version": 1, "outputs": {}, "workflows": {}}
        self.save_processing_state()
        
        self.log_message("All geometry processing attributes reset.")
        self.log_message("Geometry reset completed successfully.")

    def _invalidate_domain_delineation_outputs(self) -> None:
        """Clear generated domain references while preserving user choices."""
        project_folder = getattr(self.dialog, "project_folder", None)
        if not project_folder:
            return

        from ..watershed.domain_state import load_state, save_state, state_path

        if not state_path(project_folder).is_file():
            return

        state = load_state(project_folder)
        output_keys = (
            "picked",
            "picked_point",
            "picked_x",
            "picked_y",
            "picked_crs",
            "catchment_area",
            "catchment_area_m2",
            "mask_path",
            "vector_path",
        )
        changed = False
        for record in state.get("outlets", {}).values():
            for key in output_keys:
                if key in record:
                    record.pop(key)
                    changed = True

        if changed:
            save_state(project_folder, state)
            self.log_message(
                "Invalidated saved domain watershed outputs; gauge, discharge, "
                "and domain selections were preserved."
            )
