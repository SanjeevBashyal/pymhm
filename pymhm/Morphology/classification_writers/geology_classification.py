# -*- coding: utf-8 -*-
"""mHM geology classdefinition writer."""
from __future__ import annotations

from ..common import (
    os,
    QMessageBox,
    morph_folder,
    QgsVectorLayer,
)
from ..core.base import BaseProcessingMixin


class GeologyClassificationWriterMixin(BaseProcessingMixin):
    """mHM geology classdefinition writer."""

    def geology_classification_writer(self) -> str | None:
        """Read the selected geology lookup table layer and write ``geology_classdefinition.txt``.

        Required lookup fields::

            Geo-Class,Karstic,Default_Parameter_Value
            1,0,100
            2,0,100

        Generated output format::

            nGeo_Formations  <count>
            GeoParam(i)   ClassUnit     Karstic      Description
                     1          1           0      GeoUnit-1
                     2          2           0      GeoUnit-2
            ...
            !<-END
            !***********************************
            ! NOTES
            !***********************************
            1 = Karstic
            0 = Non-karstic

        The generated classes must follow the ``geoparameter`` ordering in
        ``mhm_parameter.nml``.
        """
        self.log_message("\n--- Writing Geology Classification Definition File ---")
        
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        lookup_layer_combo = getattr(
            self.dialog, "mMapLayerComboBox_geologyLookup", None)
        lookup_layer = (
            lookup_layer_combo.currentLayer()
            if lookup_layer_combo is not None else None
        )

        if not lookup_layer:
            self.log_message(
                "ERROR: Geology lookup table layer not selected.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Geology lookup table layer.")
            return
        
        morph_output_folder = morph_folder(self.dialog.project_folder)
        os.makedirs(morph_output_folder, exist_ok=True)
        output_file = os.path.join(
            morph_output_folder, "geology_classdefinition.txt")
        
        try:
            import pandas as pd
            
            df = self._lookup_layer_to_dataframe(
                lookup_layer, pd, "geology")
            
            if df is None or df.empty:
                self.log_message("ERROR: Geology lookup table layer is empty or could not be read.")
                QMessageBox.critical(
                    self.dialog, "Error",
                    "Geology lookup table layer is empty or could not be read.")
                return None
            
            # Log column names found
            self.log_message(f"Geology lookup fields found: {', '.join(df.columns.tolist())}")
            
            # Check if required columns exist (case-insensitive)
            column_names_lower = {col.lower(): col for col in df.columns}
            geo_class_col = column_names_lower.get('geo-class')
            karstic_col = column_names_lower.get('karstic')
            default_param_col = column_names_lower.get('default_parameter_value')
            
            if not geo_class_col:
                self.log_message("ERROR: 'Geo-Class' field not found in geology lookup table.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"Geology lookup table must contain a 'Geo-Class' field.\nFound fields: {', '.join(df.columns.tolist())}")
                return None
            
            if not karstic_col:
                self.log_message("ERROR: 'Karstic' field not found in geology lookup table.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"Geology lookup table must contain a 'Karstic' field.\nFound fields: {', '.join(df.columns.tolist())}")
                return None
            
            if not default_param_col:
                self.log_message("ERROR: 'Default_Parameter_Value' field not found in geology lookup table.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"Geology lookup table must contain a 'Default_Parameter_Value' field.\nFound fields: {', '.join(df.columns.tolist())}")
                return None
            
            self.log_message(f"Using columns: '{geo_class_col}', '{karstic_col}', '{default_param_col}'")
            
            # Clean and sort data
            df_clean = df[[geo_class_col, karstic_col, default_param_col]].dropna()
            
            # Sort by Geo-Class to ensure proper ordering
            df_clean = df_clean.sort_values(by=geo_class_col)
            
            # Get count of formations
            n_geo_formations = len(df_clean)
            
            self.log_message(f"Found {n_geo_formations} geological formations")
            
            # Write output file
            with open(output_file, 'w', encoding='utf-8') as f:
                # Write header
                f.write(f"nGeo_Formations  {n_geo_formations}\n")
                f.write("\n")
                f.write("GeoParam(i)   ClassUnit     Karstic      Description\n")
                
                # Write data rows
                for idx, row in df_clean.iterrows():
                    geo_class = int(float(row[geo_class_col])) if pd.notna(row[geo_class_col]) else 0
                    karstic = int(float(row[karstic_col])) if pd.notna(row[karstic_col]) else 0
                    default_param = row[default_param_col] if pd.notna(row[default_param_col]) else 100
                    description = f"GeoUnit-{geo_class}"
                    
                    # Format: GeoParam(i)   ClassUnit     Karstic      Description
                    # Match the exact format from the example
                    #         1	          1           0      GeoUnit-1
                    f.write(f"{geo_class:10d}\t{geo_class:10d}     {karstic:10d}      {description}\n")
                
                # Write footer
                f.write("!<-END\n")
                f.write("!***********************************\n")
                f.write("! NOTES\n")
                f.write("!***********************************\n")
                f.write("1 = Karstic\n")
                f.write("0 = Non-karstic\n")
                f.write("IMPORTANT ::\n")
                f.write("   Ordering has to be according to the ordering in mhm_parameter.nml\n")
                f.write("   (namelist: geoparameter)\n")
            
            self.log_message(f"Geology classification definition file written successfully: {output_file}")
            self.mark_output_prepared(
                output_file,
                name="geology_classdefinition.txt",
                loaded=False
            )
            QMessageBox.information(
                self.dialog, "Success",
                f"Geology classification definition file created successfully.\n{output_file}")
            return output_file
            
        except ImportError:
            self.log_message("ERROR: pandas library is not installed.")
            QMessageBox.critical(
                self.dialog, "Error",
                "pandas library is required but not installed.\nPlease install it using: pip install pandas")
            return None
        except Exception as e:
            import traceback
            error_details = traceback.format_exc()
            self.log_message(
                f"ERROR writing geology classification file: {e}\n{error_details}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error writing geology classification file:\n{str(e)}")
            return None

    def _lookup_layer_to_dataframe(self, lookup_layer, pd, layer_label):
        """Convert a selected QGIS lookup table layer to a pandas DataFrame."""
        if not isinstance(lookup_layer, QgsVectorLayer) or not lookup_layer.isValid():
            self.log_message(f"ERROR: Selected {layer_label} lookup layer is not valid.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                f"Please select a valid {layer_label} lookup table layer.")
            return None

        field_names = lookup_layer.fields().names()
        rows = []
        for feature in lookup_layer.getFeatures():
            rows.append({
                field_name: feature[field_name]
                for field_name in field_names
            })

        self.log_message(
            f"Read {len(rows)} rows from {layer_label} lookup layer '{lookup_layer.name()}'.")
        return pd.DataFrame(rows, columns=field_names)
