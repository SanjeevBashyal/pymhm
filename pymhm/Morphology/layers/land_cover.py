# -*- coding: utf-8 -*-
"""Land-cover clipping and mHM class raster preparation."""
from ..common import (
    os,
    QMessageBox,
)


class LandCoverProcessingMixin:
    """Land-cover clipping and mHM class raster preparation."""

    def reclassify_land_use_raster(self, input_raster_path, output_path, lookup_mapping):
        """
        Reclassify land use raster values based on lookup table.

        Args:
            input_raster_path: Path to clipped land use raster
            output_path: Path for reclassified output
            lookup_mapping: Dictionary mapping grid values to type_int (1, 2, or 3)

        Returns:
            bool: True if successful, False otherwise
        """
        self.log_message(
            "Reclassifying land use raster based on lookup table...")

        # Build reclassification table
        # Format: [min1, max1, value1, min2, max2, value2, ...]
        reclass_table = []
        for grid_value, type_int in sorted(lookup_mapping.items()):
            reclass_table.extend(
                [str(grid_value-0.5), str(grid_value+0.5), str(type_int)])

        self.log_message(f"Reclassification table: {reclass_table}")

        params_reclass = {
            'INPUT_RASTER': input_raster_path,
            'RASTER_BAND': 1,
            'TABLE': reclass_table,
            'NO_DATA': -9999,
            'RANGE_BOUNDARIES': 0,  # min <= value <= max
            'NODATA_FOR_MISSING': False,
            'DATA_TYPE': 3,  # Int16 (since values are 1, 2, or 3)
            'CREATE_OPTIONS': None,
            'OUTPUT': output_path
        }

        result = self.run_processing_algorithm(
            "native:reclassifybytable", params_reclass)

        if result and os.path.exists(output_path):
            self.log_message(
                "Land use reclassification completed successfully.")
            return True
        else:
            self.log_message("ERROR: Land use reclassification failed.")
            return False

    def process_land_use(self):
        """Process Land Use layer with clipping and reclassification"""
        if not self.check_prerequisites():
            return

        layer = self.dialog.mMapLayerComboBox_land_cover.currentLayer()
        if not layer:
            QMessageBox.warning(self.dialog, "Input Error",
                                "Please select a Land Cover layer.")
            return

        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        clipped_path = os.path.join(geometry_folder, "3_land_use_clipped.tif")
        final_path = os.path.join(geometry_folder, "3_land_use.tif")

        # Check if final land use layer already exists
        if os.path.exists(final_path):
            self.log_message(
                "Final land use layer already exists. Loading existing file...")
            self.land_use_layer = final_path
            self.load_layer(final_path, "3_Land_Use")
            return

        # Step 1: Clip land use layer to DEM extent (display=False for intermediate layer)
        self.log_message("Clipping land use layer to DEM extent...")
        result = self.process_layer_with_dem_clipping(
            layer, "3_Land_Use_Clip", "3_land_use_clipped.tif", display=False)

        if not result:
            self.log_message("ERROR: Land use clipping failed.")
            return

        # Step 2: Load lookup table
        lookup_mapping = self.load_land_cover_lookup_table()

        # Step 3: Reclassify clipped raster
        success = self.reclassify_land_use_raster(
            result, final_path, lookup_mapping)

        if success:
            self.land_use_layer = final_path
            self.load_layer(final_path, "3_Land_Use")
            self.log_message(
                "Land Use layer processed, clipped, and reclassified successfully.")

            # Clean up intermediate clipped file
            try:
                if os.path.exists(clipped_path):
                    os.remove(clipped_path)
            except:
                pass
        else:
            self.log_message("ERROR: Land use reclassification failed.")
