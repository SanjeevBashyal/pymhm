# -*- coding: utf-8 -*-
"""ASCII export of masked morphology rasters."""
from ..common import (
    os,
    QMessageBox,
)


class AsciiExportMixin:
    """ASCII export of masked morphology rasters."""

    def write_all_layers(self):
        """
        Convert all masked layers to ASCII format using gdal:translate.
        Finds all masked .tif files in the geometry folder and converts them to .asc files
        with standardized names in the input/morph folder.
        """
        self.log_message("\n--- Converting all masked layers to ASCII format ---")
        
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        
        if not os.path.exists(geometry_folder):
            QMessageBox.warning(
                self.dialog, "Error",
                "Geometry folder does not exist. Please process layers first.")
            return
        
        # Output folder - save to input/morph folder
        input_folder = os.path.join(self.dialog.project_folder, "input")
        morph_folder = os.path.join(input_folder, "morph")
        os.makedirs(morph_folder, exist_ok=True)
        
        # Mapping from masked layer filenames to output ASCII filenames
        filename_mapping = {
            "1_dem_filled_masked.tif": "dem.asc",
            "1_dem_slope_masked.tif": "slope.asc",
            "1_dem_aspect_masked.tif": "aspect.asc",
            "2_flow_accumulation_masked.tif": "facc.asc",
            "2_flow_direction_masked.tif": "fdir.asc",
            "2_gauge_position_masked.tif": "idgauges.asc",
            "3_land_use_masked.tif": "luse.asc",
            "3_soil_masked.tif": "soil_class.asc",
            "3_geology_processed_masked.tif": "geology_class.asc"
        }
        
        # Find all masked .tif files in the geometry folder
        masked_layers = []
        for filename in os.listdir(geometry_folder):
            if filename.endswith("_masked.tif"):
                input_path = os.path.join(geometry_folder, filename)
                
                # Get output filename from mapping, or use default if not found
                output_filename = filename_mapping.get(filename)
                if not output_filename:
                    # Default: remove "_masked" and change extension to .asc
                    base_name = filename.replace("_masked.tif", "")
                    output_filename = f"{base_name}.asc"
                    self.log_message(f"WARNING: No mapping found for {filename}, using default name: {output_filename}")
                
                output_path = os.path.join(morph_folder, output_filename)
                layer_name = filename.replace("_masked.tif", "").replace("_", " ").title()
                masked_layers.append({
                    'input': input_path,
                    'output': output_path,
                    'name': layer_name,
                    'output_filename': output_filename
                })
        
        if not masked_layers:
            self.log_message("No masked layers found. Running Mask All Layers first...")
            self.without_layer_loading(self.mask_all_layers)
            for filename in os.listdir(geometry_folder):
                if filename.endswith("_masked.tif"):
                    input_path = os.path.join(geometry_folder, filename)
                    output_filename = filename_mapping.get(filename)
                    if not output_filename:
                        base_name = filename.replace("_masked.tif", "")
                        output_filename = f"{base_name}.asc"
                        self.log_message(
                            f"WARNING: No mapping found for {filename}, using default name: {output_filename}")

                    output_path = os.path.join(morph_folder, output_filename)
                    layer_name = filename.replace("_masked.tif", "").replace("_", " ").title()
                    masked_layers.append({
                        'input': input_path,
                        'output': output_path,
                        'name': layer_name,
                        'output_filename': output_filename
                    })

            if not masked_layers:
                return
        
        self.log_message(f"Found {len(masked_layers)} masked layer(s) to convert...")
        
        # Convert each masked layer to ASCII
        success_count = 0
        failed_count = 0
        
        for layer_info in masked_layers:
            input_path = layer_info['input']
            output_path = layer_info['output']
            layer_name = layer_info['name']
            output_filename = layer_info['output_filename']
            
            # Skip if output already exists
            if os.path.exists(output_path):
                self.log_message(f"{layer_name} ASCII file ({output_filename}) already exists. Skipping...")
                success_count += 1
                continue
            
            self.log_message(f"Converting {layer_name} to {output_filename}...")
            nodata_value = 247 if output_filename == "fdir.asc" else -9999
            
            # Use gdal:translate to convert to ASCII
            params_translate = {
                'INPUT': input_path,
                'TARGET_CRS': None,
                'NODATA': nodata_value,
                'COPY_SUBDATASETS': False,
                'OPTIONS': None,
                'EXTRA': '',
                'DATA_TYPE': 0,  # Use input data type
                'OUTPUT': output_path
            }
            
            result = self.run_processing_algorithm("gdal:translate", params_translate)
            
            if result and os.path.exists(output_path):
                self.log_message(f"{layer_name} converted to {output_filename} successfully.")
                success_count += 1
            else:
                self.log_message(f"ERROR: Failed to convert {layer_name} to {output_filename}.")
                failed_count += 1
        
        # Summary
        self.log_message(f"\n--- Conversion Summary ---")
        self.log_message(f"Successfully converted: {success_count} layer(s)")
        if failed_count > 0:
            self.log_message(f"Failed to convert: {failed_count} layer(s)")
        
        if success_count > 0:
            self.log_message(f"ASCII files saved to: {morph_folder}")
            self.log_message("ASCII conversion process completed.")
        else:
            QMessageBox.warning(
                self.dialog, "Conversion Error",
                "No layers were successfully converted to ASCII format.")
