# -*- coding: utf-8 -*-
"""Discovery and restoration of prepared project morphology outputs."""
from ..common import (
    os,
    project_geometry_folder,
    QgsRasterLayer,
    processing,
)


class ProjectStateMixin:
    """Discovery and restoration of prepared project morphology outputs."""

    def load_project_state(self):
        """Checks for existing output files and sets instance attributes to resume work."""
        self.log_message("\n--- Checking for existing project files... ---")
        self.load_processing_state()

        files_to_check = {
            # Original processing files
            'filled_dem_path': ("1_dem_filled.tif", "1_dem_filled.sdat"),
            'flow_dir_path': ("2_flow_direction.tif", "2_flow_direction.sdat"),
            'channel_network_vector_path': "2_channel_network.shp",
            'snapped_points_path': "2_pour_points_snapped.shp",
            'watershed_raster_path': ("4_watershed_raster.tif", "1_watershed_raster.sdat"),
            'watershed_vector_path': "1_watershed_final.shp",
            'merged_watershed_path': (
                os.path.join("Watersheds", "4_watershed_merged_vector.shp"),
                "4_watershed_merged_vector.shp"
            ),

            # New hydrological processing files
            'aspect_path': ("1_dem_aspect.tif", "1_dem_aspect.sdat"),
            'slope_path': ("1_dem_slope.tif", "1_dem_slope.sdat"),
            'flow_direction_path': "2_flow_direction.tif",
            'flow_accumulation_path': "2_flow_accumulation.tif",
            'flow_accumulation_area_path': ("2_flow_accumulation_area.tif", "2_flow_accumulation_area.sdat"),
            'gauge_position_path': "2_gauge_position.tif"
        }

        found_any = False
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        for attr, filenames in files_to_check.items():
            if isinstance(filenames, str):
                filenames = (filenames,)

            matched_path = None
            matched_filename = None
            for filename in filenames:
                expected_path = os.path.join(geometry_folder, filename)
                if os.path.exists(expected_path):
                    matched_path = expected_path
                    matched_filename = filename
                    break

            if matched_path:
                setattr(self, attr, matched_path)
                self.mark_output_prepared(matched_path, name=matched_filename)
                self.log_message(f"Found existing file: {matched_filename}")
                found_any = True
            else:
                setattr(self, attr, None)

        # Check for final land use layer separately
        land_use_path = os.path.join(geometry_folder, "3_land_use.tif")
        if os.path.exists(land_use_path):
            self.land_use_layer = land_use_path
            self.mark_output_prepared(land_use_path, name="3_land_use.tif")
            self.log_message(f"Found existing file: 3_land_use.tif")
            found_any = True

        # Check for reprojected DEM and load it into self.dem_layer (but don't add to project)
        reprojected_dem_path = os.path.join(
            geometry_folder, "0_dem_reprojected.tif")
        if os.path.exists(reprojected_dem_path):
            reprojected_layer = QgsRasterLayer(
                reprojected_dem_path, "1_DEM_Reprojected")
            if reprojected_layer.isValid():
                self.dem_layer = reprojected_layer
                self.log_message(
                    f"Found existing reprojected DEM: 0_dem_reprojected.tif")
                found_any = True
            else:
                self.log_message(
                    f"WARNING: Reprojected DEM file exists but is not valid.")

        # Keep filled DEM as an internal intermediate; do not add it to QGIS.
        if self.filled_dem_path and os.path.exists(self.filled_dem_path):
            self.log_message(
                "Found existing filled DEM intermediate: 1_dem_filled.tif")

        if found_any:
            self.log_message(
                "Project state loaded. You may skip completed steps.")
        else:
            self.log_message(
                "No existing project files found in 'Geometry' folder.")

    # --- Geometry Processing Methods ---
