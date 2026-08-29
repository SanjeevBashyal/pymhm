# -*- coding: utf-8 -*-
"""Execute-all morphology workflow orchestration."""
from __future__ import annotations

from ..common import (QMessageBox, morph_staging_folder, os,
                      project_geometry_folder)
from ..hydrology.aggregate import HydrologyMixin
from ..hydrology.outlets import configured_gauged_outlet_ids
from ..latlon import LatLonMixin
from ..layers import LayerProcessingMixin
from ..watershed import WatershedMixin


class ExecuteAllMixin(
    WatershedMixin,
    HydrologyMixin,
    LayerProcessingMixin,
    LatLonMixin,
):
    """Execute-all morphology workflow orchestration."""

    def execute_all_processing(self, show_error_dialog=True) -> bool:
        """Run every Execute All step synchronously on the calling thread.

        The `pushButton_executeAllMorphology` action does NOT come here: it goes
        through `MorphologyTaskBridge.start_execute_all()`, which runs the same
        stages as background `QgsTask` jobs. Keep the two stage lists in step --
        a stage added only here will silently never run from the UI.

        The bridge covers steps 1-3, 8 and 9 with a single mHM-tools
        `create-dem-derivatives` pass; the steps below reach the same files one
        at a time.

        Steps, in order:
        1. filled dem
        2. slope
        3. aspect
        4. land cover
        5. soil
        6. geology
        7. LAI resampled to the filled DEM grid
        8. flow_accumulation
        9. flow_direction
        10. idgauges (will be processed after snap points due to dependency)
        11. channel network
        12. snap
        13. reuse an existing manually delineated domain mask
        14. verify geology class definition
        15. verify soil outputs
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
            self.log_message("\n--- Step 1/15: Fill DEM ---")
            self.without_layer_loading(self.fill_dem)
            if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
                return fail("Fill DEM failed. Aborting Execute All.")

            # Step 2: Slope
            self.log_message("\n--- Step 2/15: Process Slope ---")
            self.process_slope()

            # Step 3: Aspect
            self.log_message("\n--- Step 3/15: Process Aspect ---")
            self.process_aspect()

            # Step 4: Land Cover
            self.log_message("\n--- Step 4/15: Process Land Cover ---")
            if not self.process_land_use():
                return fail("Land Cover processing failed. Aborting Execute All.")
            land_cover_ready = self._categorical_mode("lc") == "mhm_ready"
            advanced_land_cover = bool(
                getattr(
                    self.dialog,
                    "uses_advanced_categorical_input",
                    lambda _kind: False,
                )("lc")
            )
            if land_cover_ready and not advanced_land_cover and not os.path.exists(
                self.categorical_ready_outputs.get("lc", "")
            ):
                return fail(
                    "mHM-ready Land Cover output is missing. Aborting Execute All."
                )

            # Step 5: Soil
            self.log_message("\n--- Step 5/15: Process Soil ---")
            soil_raster = os.path.join(
                project_geometry_folder(self.dialog.project_folder),
                "3_soil.tif",
            )
            soil_definition = os.path.join(
                morph_staging_folder(self.dialog.project_folder),
                "soil_classdefinition.txt",
            )
            if not self.process_soil(write_classdefinition=True):
                return fail("Soil processing failed. Aborting Execute All.")
            advanced_soil = bool(
                getattr(
                    self.dialog,
                    "uses_advanced_categorical_input",
                    lambda _kind: False,
                )("soil")
            )
            if advanced_soil:
                version = self.dialog.comboBox_mHMversion.currentText().strip()
                if version.startswith("5"):
                    soil_raster = os.path.join(
                        morph_staging_folder(self.dialog.project_folder),
                        "soil_class.asc",
                    )
                    soil_definition = os.path.join(
                        morph_staging_folder(self.dialog.project_folder),
                        "soil_classdefinition.txt",
                    )
                else:
                    soil_raster = os.path.join(
                        morph_staging_folder(self.dialog.project_folder),
                        "soil_horizon_class.nc",
                    )
                    soil_definition = os.path.join(
                        morph_staging_folder(self.dialog.project_folder),
                        "soil_classdefinition_iFlag_soilDB_1.txt",
                    )
            elif self._categorical_mode("soil") == "mhm_ready":
                soil_raster = self.categorical_ready_outputs.get("soil", "")
            if not os.path.exists(soil_raster) or not os.path.exists(
                    soil_definition):
                return fail(
                    "Soil processing did not create both required outputs. "
                    "Aborting Execute All.")

            # Step 6: Geology
            self.log_message("\n--- Step 6/15: Process Geology ---")
            geology_raster = os.path.join(
                project_geometry_folder(self.dialog.project_folder),
                "3_geology_processed.tif",
            )
            geology_definition = os.path.join(
                morph_staging_folder(self.dialog.project_folder),
                "geology_classdefinition.txt",
            )
            geology_metadata = os.path.join(
                project_geometry_folder(self.dialog.project_folder),
                "geology_class_metadata.json",
            )
            if not self.process_geology(write_classdefinition=True):
                return fail("Geology processing failed. Aborting Execute All.")
            if self._categorical_mode("geology") == "mhm_ready":
                geology_raster = self.categorical_ready_outputs.get("geology", "")
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

            # Step 7: LAI on the filled DEM grid
            self.log_message(
                "\n--- Step 7/15: Resample LAI to the filled DEM grid ---"
            )
            if not self.resample_lai_to_dem_grid():
                return fail("LAI resampling failed. Aborting Execute All.")

            # Step 8: Flow Accumulation
            self.log_message("\n--- Step 8/15: Process Flow Accumulation ---")
            self.process_flow_accumulation()
            if not self.flow_accumulation_path or not os.path.exists(
                self.flow_accumulation_path
            ):
                return fail("Flow Accumulation failed. Aborting Execute All.")

            # Step 9: Flow Direction
            self.log_message("\n--- Step 9/15: Process Flow Direction ---")
            self.process_flow_direction()
            if not self.flow_direction_path or not os.path.exists(
                self.flow_direction_path
            ):
                return fail("Flow Direction failed. Aborting Execute All.")

            # Step 10: ID Gauges (Gauge Position) - Note: requires snap points, will be processed after step 11
            self.log_message(
                "\n--- Step 10/15: Process ID Gauges (deferred until after snap points) ---"
            )
            # This will be processed after snap points in step 12

            # Step 11: Channel Network
            self.log_message("\n--- Step 11/15: Process Channel Network ---")
            self.process_channel_network()
            if not self.channel_network_vector_path or not os.path.exists(
                self.channel_network_vector_path
            ):
                return fail("Channel Network failed. Aborting Execute All.")

            # Step 12: Snap Points
            self.log_message("\n--- Step 12/15: Snap Points ---")
            if not self.check_prerequisites(needs_pour_points=True):
                self.log_message(
                    "WARNING: Pour points not available. Skipping Snap Points step."
                )
            else:
                self.snap_points()
                # Now process ID Gauges (Step 10) after snap points are created
                if self.snapped_points_path and os.path.exists(
                    self.snapped_points_path
                ):
                    configured_ids = configured_gauged_outlet_ids(
                        self.dialog.project_folder
                    )
                    if configured_ids == []:
                        self.log_message(
                            "No outlets are marked as gauged. Skipping ID Gauges."
                        )
                    else:
                        self.log_message(
                            "\n--- Processing Step 10/15: ID Gauges (now that snap points are available) ---"
                        )
                        self.process_gauge_position()

            # Step 13: domain masks are deliberately created in the interactive
            # Domain Delineator, never by Execute All.
            self.log_message("\n--- Step 13/15: Check Domain Watershed Mask ---")
            merged_watershed = self._restore_existing_path(
                "merged_watershed_path",
                os.path.join(
                    "Watersheds",
                    "4_watershed_merged_vector.shp",
                ),
                "4_watershed_merged_vector.shp",
            )
            if merged_watershed:
                self.log_message(
                    "Reusing the saved merged domain watershed mask."
                )
            else:
                self.log_message(
                    "Domain delineation is handled separately. Use Domain "
                    "Delineator before Meteo & Morph Setup."
                )

            # Step 14: Verify Geology Class Definition
            self.log_message("\n--- Step 14/15: Verify Geology Class Definition ---")
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

            # Step 15: Verify Soil Outputs
            self.log_message("\n--- Step 15/15: Verify Soil Outputs ---")
            if not os.path.exists(soil_raster) or not os.path.exists(
                    soil_definition):
                return fail(
                    "Required soil outputs are missing. Aborting Execute All.")

            self.log_message("\n=== Morphology Preparation Completed Successfully ===")
            self.mark_workflow_status(
                "execute_all",
                "completed",
                "Morphology preparation completed successfully.",
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

    def execute_morph_setup_processing(
            self,
            show_error_dialog=True,
            workflow_key="morph_setup") -> bool:
        """Run crop, mask, lat/lon, and ASCII export as one workflow."""
        self.log_message("\n=== Starting Morphology Setup ===")

        if not self.check_prerequisites():
            message = "Prerequisites check failed. Aborting Morphology Setup."
            self.log_message(f"ERROR: {message}")
            self.mark_workflow_status(workflow_key, "failed", message)
            return False

        self.skip_loading = True
        self.mark_workflow_status(workflow_key, "running")
        try:
            def fail(message):
                self.log_message(f"ERROR: {message}")
                self.mark_workflow_status(workflow_key, "failed", message)
                return False

            steps = (
                (
                    "crop", "Crop All Layers",
                    lambda: self.crop_all_layers(
                        show_error_dialog=show_error_dialog),
                ),
                (
                    "mask", "Mask All Layers",
                    lambda: self.mask_all_layers(
                        show_error_dialog=show_error_dialog),
                ),
                ("latlon", "Create latlon.nc", self.process_lat_lon),
                (
                    "write", "Write All Layers",
                    lambda: self.write_all_layers(
                        show_error_dialog=show_error_dialog),
                ),
                # Publish first: it moves the staged advanced/mHM-ready
                # rasters into data/master on the L0 grid, and Write Domain
                # DEMs copies that shared set into every domain folder.
                (
                    "publish", "Publish Model Inputs",
                    self.align_advanced_inputs_to_l0,
                ),
                (
                    "domain_dems", "Write Domain DEMs",
                    self.write_domain_dems_to_l0,
                ),
            )
            total = len(steps)
            for index, (step, label, run) in enumerate(steps, start=1):
                status_key = f"{workflow_key}_{step}"
                self.log_message(
                    f"\n--- Morphology Setup Step {index}/{total}: {label} ---")
                self.mark_workflow_status(status_key, "running")
                if not run():
                    message = f"{label} failed. Aborting Morphology Setup."
                    self.mark_workflow_status(status_key, "failed", message)
                    return fail(message)
                self.mark_workflow_status(status_key, "completed")

            self.log_message("\n=== Morphology Setup Completed Successfully ===")
            self.mark_workflow_status(
                workflow_key,
                "completed",
                "Morphology Setup completed successfully.",
            )
            return True

        except Exception as e:
            message = f"Morphology Setup failed with exception: {str(e)}"
            self.log_message(f"\nERROR: {message}")
            import traceback

            self.log_message(f"Traceback: {traceback.format_exc()}")
            self.mark_workflow_status(workflow_key, "failed", message)
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
