# -*- coding: utf-8 -*-
"""ASCII export of masked morphology rasters."""
from __future__ import annotations

from ..common import (
    os,
    QMessageBox,
    geometry_folder,
    morph_folder,
)
from ..layers.masking import MaskingMixin
from ...mhm_tools_to_integrate.data_processing import (
    MorphologyAsciiLayer,
    prepare_morphology_ascii_files,
)


class AsciiExportMixin(MaskingMixin):
    """ASCII export of masked morphology rasters."""

    def write_all_layers(self, show_error_dialog=True) -> bool:
        """
        Convert all masked layers to mHM ASCII format through mhm_tools.
        The output grid is cropped/aligned to the L0 header derived from the
        common L2 model extent, and L0/L1/L11/L2 sizes are checked first.
        """
        self.log_message("\n--- Converting all masked layers to ASCII format ---")

        # Check prerequisites
        if not self.check_prerequisites():
            return False

        try:
            headers = self.dialog.grid_level_headers()
        except Exception as e:
            self.log_message(f"ERROR: Cannot prepare ASCII grid headers: {e}")
            if show_error_dialog:
                QMessageBox.warning(
                    self.dialog,
                    "Grid Configuration Error",
                    str(e),
                )
            return False

        geom_folder = geometry_folder(self.dialog.project_folder)

        if not os.path.exists(geom_folder):
            if show_error_dialog:
                QMessageBox.warning(
                    self.dialog, "Error",
                    "Geometry folder does not exist. Please process layers first.")
            return False

        morph_output_folder = morph_folder(self.dialog.project_folder)
        os.makedirs(morph_output_folder, exist_ok=True)

        # Mapping from masked layer filenames to output ASCII settings.
        filename_mapping = {
            "1_dem_filled_masked.tif": ("dem.asc", -9999, False),
            "1_dem_slope_masked.tif": ("slope.asc", -9999, False),
            "1_dem_aspect_masked.tif": ("aspect.asc", -9999, False),
            "2_flow_accumulation_masked.tif": ("facc.asc", -9999, True),
            "2_flow_direction_masked.tif": ("fdir.asc", -9999, True),
            "2_gauge_position_masked.tif": ("idgauges.asc", -9999, True),
            "3_land_use_masked.tif": ("lc.asc", -9999, True),
            "3_soil_masked.tif": ("soil_class.asc", -9999, True),
            "3_geology_processed_masked.tif": ("geology_class.asc", -9999, True),
        }

        # Find all masked .tif files in the geometry folder
        masked_layers = []
        for filename in os.listdir(geom_folder):
            if filename.endswith("_masked.tif"):
                input_path = os.path.join(geom_folder, filename)

                # Get output settings from mapping, or use default if not found
                output_settings = filename_mapping.get(filename)
                if output_settings:
                    output_filename, nodata_value, integer = output_settings
                else:
                    # Default: remove "_masked" and change extension to .asc
                    base_name = filename.replace("_masked.tif", "")
                    output_filename = f"{base_name}.asc"
                    nodata_value = -9999
                    integer = False
                    self.log_message(f"WARNING: No mapping found for {filename}, using default name: {output_filename}")

                output_path = os.path.join(morph_output_folder, output_filename)
                layer_name = filename.replace("_masked.tif", "").replace("_", " ").title()
                masked_layers.append(
                    MorphologyAsciiLayer(
                        input_path=input_path,
                        output_path=output_path,
                        name=layer_name,
                        nodata_value=nodata_value,
                        integer=integer,
                    )
                )

        if not masked_layers:
            self.log_message("No masked layers found. Please run Mask All before Write All.")
            if show_error_dialog:
                QMessageBox.warning(
                    self.dialog,
                    "No Masked Layers",
                    "No masked raster layers found. Please run Mask All first.")
            return False

        self.log_message(f"Found {len(masked_layers)} masked layer(s) to convert...")

        try:
            result = prepare_morphology_ascii_files(
                layers=masked_layers,
                headers=headers,
                overwrite=True,
                log=self.log_message,
            )
        except Exception as e:
            self.log_message(f"ERROR: ASCII conversion failed: {e}")
            if show_error_dialog:
                QMessageBox.warning(
                    self.dialog, "Conversion Error",
                    f"Failed to prepare ASCII files:\n{e}")
            return False

        for output_path in result.outputs.values():
            self.mark_output_prepared(
                str(output_path),
                name=os.path.basename(str(output_path)),
                loaded=False,
                algorithm="mhm_tools.common.file_handler",
            )

        self.log_message(f"\n--- Conversion Summary ---")
        self.log_message(f"Successfully converted: {len(result.outputs)} layer(s)")
        self.log_message(f"ASCII files saved to: {morph_output_folder}")
        self.log_message("ASCII conversion process completed.")
        return True
