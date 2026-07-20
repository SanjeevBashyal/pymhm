# -*- coding: utf-8 -*-
"""Execute-all morphology workflow orchestration."""
from __future__ import annotations

from ..common import (
    QMessageBox,
    morph_folder,
    os,
    project_geometry_folder,
)
from ..classification_writers import ClassificationWritersMixin
from ..hydrology.aggregate import HydrologyMixin
from ..latlon import LatLonMixin
from ..layers import LayerProcessingMixin
from ..watershed import WatershedMixin


class ExecuteAllMixin(
    WatershedMixin,
    HydrologyMixin,
    LayerProcessingMixin,
    LatLonMixin,
    ClassificationWritersMixin,
):
    """Execute-all morphology workflow orchestration."""

    def execute_all_processing(self, show_error_dialog=True) -> bool:
        """
        Execute all processing steps sequentially in the following order:
        1. filled dem
        2. slope
        3. aspect
        4. land cover
        5. soil
        6. geology
        7. flow_accumulation
        8. flow_direction
        9. idgauges (will be processed after snap points due to dependency)
        10. channel network
        11. snap
        12. upslope area
        13. crop all layers
        14. mask all cropped layers
        15. process lat/lon headers
        16. verify geology class definition
        17. verify soil outputs
        18. write all layers (convert to ASCII)
        """
        self.log_message("\n=== Starting Execute All Processing ===")

        # Check prerequisites first
        if not self.check_prerequisites():
            message = "Prerequisites check failed. Aborting Execute All."
            self.log_message(f"ERROR: {message}")
            self.mark_workflow_status("execute_all", "failed", message)
            return False

        self.skip_loading = True
        self.mark_workflow_status("execute_all", "running")
        try:
            def fail(message):
                self.log_message(f"ERROR: {message}")
                self.mark_workflow_status("execute_all", "failed", message)
                return False

            # Step 1: Fill DEM
            self.log_message("\n--- Step 1/18: Fill DEM ---")
            self.without_layer_loading(self.fill_dem)
            if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
                return fail("Fill DEM failed. Aborting Execute All.")

            # Step 2: Slope
            self.log_message("\n--- Step 2/18: Process Slope ---")
            self.process_slope()

            # Step 3: Aspect
            self.log_message("\n--- Step 3/18: Process Aspect ---")
            self.process_aspect()

            # Step 4: Land Cover
            self.log_message("\n--- Step 4/18: Process Land Cover ---")
            self.process_land_use()

            # Step 5: Soil
            self.log_message("\n--- Step 5/18: Process Soil ---")
            soil_raster = os.path.join(
                project_geometry_folder(self.dialog.project_folder),
                "3_soil.tif",
            )
            soil_definition = os.path.join(
                morph_folder(self.dialog.project_folder),
                "soil_classdefinition.txt",
            )
            if not self.process_soil(write_classdefinition=True):
                return fail("Soil processing failed. Aborting Execute All.")
            if not os.path.exists(soil_raster) or not os.path.exists(
                    soil_definition):
                return fail(
                    "Soil processing did not create both required outputs. "
                    "Aborting Execute All.")

            # Step 6: Geology
            self.log_message("\n--- Step 6/18: Process Geology ---")
            geology_raster = os.path.join(
                project_geometry_folder(self.dialog.project_folder),
                "3_geology_processed.tif",
            )
            geology_definition = os.path.join(
                morph_folder(self.dialog.project_folder),
                "geology_classdefinition.txt",
            )
            geology_metadata = os.path.join(
                project_geometry_folder(self.dialog.project_folder),
                "geology_class_metadata.json",
            )
            if not self.process_geology(write_classdefinition=True):
                return fail("Geology processing failed. Aborting Execute All.")
            if not all(
                os.path.exists(path)
                for path in (
                    geology_raster,
                    geology_definition,
                    geology_metadata,
                )
            ):
                return fail(
                    "Geology processing did not create all required outputs. "
                    "Aborting Execute All."
                )

            # Step 7: Flow Accumulation
            self.log_message("\n--- Step 7/18: Process Flow Accumulation ---")
            self.process_flow_accumulation()
            if not self.flow_accumulation_path or not os.path.exists(
                self.flow_accumulation_path
            ):
                return fail("Flow Accumulation failed. Aborting Execute All.")

            # Step 8: Flow Direction
            self.log_message("\n--- Step 8/18: Process Flow Direction ---")
            self.process_flow_direction()
            if not self.flow_direction_path or not os.path.exists(
                self.flow_direction_path
            ):
                return fail("Flow Direction failed. Aborting Execute All.")

            # Step 9: ID Gauges (Gauge Position) - Note: requires snap points, will be processed after step 11
            self.log_message(
                "\n--- Step 9/18: Process ID Gauges (deferred until after snap points) ---"
            )
            # This will be processed after snap points in step 11

            # Step 10: Channel Network
            self.log_message("\n--- Step 10/18: Process Channel Network ---")
            self.process_channel_network()
            if not self.channel_network_vector_path or not os.path.exists(
                self.channel_network_vector_path
            ):
                return fail("Channel Network failed. Aborting Execute All.")

            # Step 11: Snap Points
            self.log_message("\n--- Step 11/18: Snap Points ---")
            if not self.check_prerequisites(needs_pour_points=True):
                self.log_message(
                    "WARNING: Pour points not available. Skipping Snap Points step."
                )
            else:
                self.snap_points()
                # Now process ID Gauges (Step 9) after snap points are created
                if self.snapped_points_path and os.path.exists(
                    self.snapped_points_path
                ):
                    self.log_message(
                        "\n--- Processing Step 9/18: ID Gauges (now that snap points are available) ---"
                    )
                    self.process_gauge_position()

            # Step 12: Upslope Area (Delineate Watershed)
            self.log_message("\n--- Step 12/18: Delineate Watershed (Upslope Area) ---")
            if not self.snapped_points_path or not os.path.exists(
                self.snapped_points_path
            ):
                self.log_message(
                    "WARNING: Snapped points not available. Skipping Upslope Area step."
                )
            else:
                self.delineate_watershed()

            # Step 13: Crop All Layers
            self.log_message("\n--- Step 13/18: Crop All Layers ---")
            if not self.crop_all_layers(show_error_dialog=show_error_dialog):
                return fail("Crop All Layers failed. Aborting Execute All.")

            # Step 14: Mask All Layers
            self.log_message("\n--- Step 14/18: Mask All Layers ---")
            if not self.mask_all_layers(show_error_dialog=show_error_dialog):
                return fail("Mask All Layers failed. Aborting Execute All.")

            # Step 15: Process Lat/Lon Headers
            self.log_message("\n--- Step 15/18: Process Lat/Lon Headers ---")
            self.process_lat_lon()

            # Step 16: Verify Geology Class Definition
            self.log_message("\n--- Step 16/18: Verify Geology Class Definition ---")
            if not all(
                os.path.exists(path)
                for path in (
                    geology_raster,
                    geology_definition,
                    geology_metadata,
                )
            ):
                return fail(
                    "Required geology outputs are missing. "
                    "Aborting Execute All."
                )

            # Step 17: Verify Soil Outputs
            self.log_message("\n--- Step 17/18: Verify Soil Outputs ---")
            if not os.path.exists(soil_raster) or not os.path.exists(
                    soil_definition):
                return fail(
                    "Required soil outputs are missing. Aborting Execute All.")

            # Step 18: Write All Layers (Convert to ASCII)
            self.log_message(
                "\n--- Step 18/18: Write All Layers (Convert to ASCII) ---"
            )
            if not self.write_all_layers(show_error_dialog=show_error_dialog):
                return fail("Write All Layers failed. Aborting Execute All.")

            self.log_message("\n=== Execute All Processing Completed Successfully ===")
            self.mark_workflow_status(
                "execute_all",
                "completed",
                "Execute All Processing completed successfully.",
            )
            return True

        except Exception as e:
            message = f"Execute All Processing failed with exception: {str(e)}"
            self.log_message(
                f"\nERROR: {message}"
            )
            import traceback

            self.log_message(f"Traceback: {traceback.format_exc()}")
            self.mark_workflow_status("execute_all", "failed", message)
            if show_error_dialog:
                QMessageBox.critical(
                    self.dialog, "Error", f"Execute All Processing failed:\n{str(e)}"
                )
            return False
        finally:
            self.skip_loading = False
            self.log_message(
                "Execute All finished. Prepared outputs were recorded in the project processing state."
            )

    def execute_morph_setup_processing(self, show_error_dialog=True) -> bool:
        """Run Crop All, Mask All, and Write All as one workflow."""
        self.log_message("\n=== Starting Morphology Setup ===")

        if not self.check_prerequisites():
            message = "Prerequisites check failed. Aborting Morphology Setup."
            self.log_message(f"ERROR: {message}")
            self.mark_workflow_status("morph_setup", "failed", message)
            return False

        self.skip_loading = True
        self.mark_workflow_status("morph_setup", "running")
        try:
            def fail(message):
                self.log_message(f"ERROR: {message}")
                self.mark_workflow_status("morph_setup", "failed", message)
                return False

            self.log_message("\n--- Morphology Setup Step 1/3: Crop All Layers ---")
            if not self.crop_all_layers(show_error_dialog=show_error_dialog):
                return fail("Crop All Layers failed. Aborting Morphology Setup.")

            self.log_message("\n--- Morphology Setup Step 2/3: Mask All Layers ---")
            if not self.mask_all_layers(show_error_dialog=show_error_dialog):
                return fail("Mask All Layers failed. Aborting Morphology Setup.")

            self.log_message(
                "\n--- Morphology Setup Step 3/3: Write All Layers ---"
            )
            if not self.write_all_layers(show_error_dialog=show_error_dialog):
                return fail("Write All Layers failed. Aborting Morphology Setup.")

            self.log_message("\n=== Morphology Setup Completed Successfully ===")
            self.mark_workflow_status(
                "morph_setup",
                "completed",
                "Morphology Setup completed successfully.",
            )
            return True

        except Exception as e:
            message = f"Morphology Setup failed with exception: {str(e)}"
            self.log_message(f"\nERROR: {message}")
            import traceback

            self.log_message(f"Traceback: {traceback.format_exc()}")
            self.mark_workflow_status("morph_setup", "failed", message)
            if show_error_dialog:
                QMessageBox.critical(
                    self.dialog, "Error", f"Morphology Setup failed:\n{str(e)}"
                )
            return False
        finally:
            self.skip_loading = False
            self.log_message(
                "Morphology Setup finished. Prepared outputs were recorded in the project processing state."
            )
