# -*- coding: utf-8 -*-
"""Execute-all morphology workflow orchestration."""
from __future__ import annotations

from ..common import (
    os,
    QMessageBox,
    processing,
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
        ClassificationWritersMixin):
    """Execute-all morphology workflow orchestration."""

    def execute_all_processing(self) -> None:
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
        13. mask all layers
        14. process lat/lon headers
        15. write geology class definition
        16. write soil class definition
        17. write all layers (convert to ASCII)
        """
        self.log_message("\n=== Starting Execute All Processing ===")
        
        # Check prerequisites first
        if not self.check_prerequisites():
            self.log_message("ERROR: Prerequisites check failed. Aborting Execute All.")
            return
        
        try:
            # Step 1: Fill DEM
            self.log_message("\n--- Step 1/17: Fill DEM ---")
            self.fill_dem()
            if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
                self.log_message("ERROR: Fill DEM failed. Aborting Execute All.")
                return
            
            # Step 2: Slope
            self.log_message("\n--- Step 2/17: Process Slope ---")
            self.process_slope()
            
            # Step 3: Aspect
            self.log_message("\n--- Step 3/17: Process Aspect ---")
            self.process_aspect()
            
            # Step 4: Land Cover
            self.log_message("\n--- Step 4/17: Process Land Cover ---")
            self.process_land_use()
            
            # Step 5: Soil
            self.log_message("\n--- Step 5/17: Process Soil ---")
            self.process_soil()
            
            # Step 6: Geology
            self.log_message("\n--- Step 6/17: Process Geology ---")
            self.process_geology()
            
            # Step 7: Flow Accumulation
            self.log_message("\n--- Step 7/17: Process Flow Accumulation ---")
            self.process_flow_accumulation()
            if not self.flow_accumulation_path or not os.path.exists(self.flow_accumulation_path):
                self.log_message("ERROR: Flow Accumulation failed. Aborting Execute All.")
                return
            
            # Step 8: Flow Direction
            self.log_message("\n--- Step 8/17: Process Flow Direction ---")
            self.process_flow_direction()
            if not self.flow_direction_path or not os.path.exists(self.flow_direction_path):
                self.log_message("ERROR: Flow Direction failed. Aborting Execute All.")
                return
            
            # Step 9: ID Gauges (Gauge Position) - Note: requires snap points, will be processed after step 11
            self.log_message("\n--- Step 9/17: Process ID Gauges (deferred until after snap points) ---")
            # This will be processed after snap points in step 11
            
            # Step 10: Channel Network
            self.log_message("\n--- Step 10/17: Process Channel Network ---")
            self.process_channel_network()
            if not self.channel_network_vector_path or not os.path.exists(self.channel_network_vector_path):
                self.log_message("ERROR: Channel Network failed. Aborting Execute All.")
                return
            
            # Step 11: Snap Points
            self.log_message("\n--- Step 11/17: Snap Points ---")
            if not self.check_prerequisites(needs_pour_points=True):
                self.log_message("WARNING: Pour points not available. Skipping Snap Points step.")
            else:
                self.snap_points()
                # Now process ID Gauges (Step 9) after snap points are created
                if self.snapped_points_path and os.path.exists(self.snapped_points_path):
                    self.log_message("\n--- Processing Step 9/17: ID Gauges (now that snap points are available) ---")
                    self.process_gauge_position()
            
            # Step 12: Upslope Area (Delineate Watershed)
            self.log_message("\n--- Step 12/17: Delineate Watershed (Upslope Area) ---")
            if not self.snapped_points_path or not os.path.exists(self.snapped_points_path):
                self.log_message("WARNING: Snapped points not available. Skipping Upslope Area step.")
            else:
                self.delineate_watershed()
            
            # Step 13: Mask All Layers
            self.log_message("\n--- Step 13/17: Mask All Layers ---")
            if not self.merged_watershed_path or not os.path.exists(self.merged_watershed_path):
                self.log_message("WARNING: Merged watershed not available. Skipping Mask All Layers step.")
            else:
                self.mask_all_layers()
            
            # Step 14: Process Lat/Lon Headers
            self.log_message("\n--- Step 14/17: Process Lat/Lon Headers ---")
            self.process_lat_lon()
            
            # Step 15: Write Geology Class Definition
            self.log_message("\n--- Step 15/17: Write Geology Class Definition ---")
            self.geology_classification_writer()

            # Step 16: Write Soil Class Definition
            self.log_message("\n--- Step 16/17: Write Soil Class Definition ---")
            self.soil_classdefinition_writer()

            # Step 17: Write All Layers (Convert to ASCII)
            self.log_message("\n--- Step 17/17: Write All Layers (Convert to ASCII) ---")
            self.write_all_layers()
            
            self.log_message("\n=== Execute All Processing Completed Successfully ===")
            
        except Exception as e:
            self.log_message(f"\nERROR: Execute All Processing failed with exception: {str(e)}")
            import traceback
            self.log_message(f"Traceback: {traceback.format_exc()}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Execute All Processing failed:\n{str(e)}")
        finally:
            self.skip_loading = False
            self.log_message("Execute All finished. Prepared outputs were recorded in the project processing state.")
