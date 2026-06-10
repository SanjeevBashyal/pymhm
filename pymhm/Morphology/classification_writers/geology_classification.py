# -*- coding: utf-8 -*-
"""mHM geology classdefinition writer."""
from __future__ import annotations

from ..common import (
    os,
    QMessageBox,
    morph_folder,
)
from ..core.base import BaseProcessingMixin


class GeologyClassificationWriterMixin(BaseProcessingMixin):
    """mHM geology classdefinition writer."""

    def geology_classification_writer(self) -> str | None:
        """
        Read geological lookup file and produce geology_classdefinition.txt.
        
        Input CSV format:
        Geo-Class,Karstic,Default_Parameter_Value
        1,0,100
        2,0,100
        ...
        
        Output format:
        nGeo_Formations  <count>
        GeoParam(i)   ClassUnit     Karstic      Description
             1	          1           0      GeoUnit-1
             2     	      2           0      GeoUnit-2
        ...
        !<-END
        !***********************************
        ! NOTES
        !***********************************
        1 = Karstic
        0 = Non-karstic
        IMPORTANT ::
           Ordering has to be according to the ordering in mhm_parameter.nml
           (namelist: geoparameter)
        """
        self.log_message("\n--- Writing Geology Classification Definition File ---")
        
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        # Get geology lookup file path
        lookup_file = self.dialog.lineEdit_geology_lookup.text()
        
        if not lookup_file or not os.path.exists(lookup_file):
            self.log_message(
                "ERROR: Geology lookup table file not found or not specified.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Geology lookup table file.")
            return
        
        morph_output_folder = morph_folder(self.dialog.project_folder)
        os.makedirs(morph_output_folder, exist_ok=True)
        output_file = os.path.join(
            morph_output_folder, "geology_classdefinition.txt")
        
        try:
            import pandas as pd
            
            # Try different encodings
            encodings = ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252']
            df = None
            used_encoding = None
            
            for encoding in encodings:
                try:
                    df = pd.read_csv(lookup_file, encoding=encoding)
                    used_encoding = encoding
                    self.log_message(f"Successfully read file with encoding: {encoding}")
                    break
                except (UnicodeDecodeError, pd.errors.EmptyDataError) as e:
                    continue
            
            if df is None or df.empty:
                self.log_message("ERROR: Could not read file with any encoding or file is empty.")
                QMessageBox.critical(
                    self.dialog, "Error",
                    "Could not read geology lookup table file. Encoding issue or file is empty.")
                return None
            
            # Log column names found
            self.log_message(f"CSV columns found: {', '.join(df.columns.tolist())}")
            
            # Check if required columns exist (case-insensitive)
            column_names_lower = {col.lower(): col for col in df.columns}
            geo_class_col = column_names_lower.get('geo-class') or column_names_lower.get('geo_class')
            karstic_col = column_names_lower.get('karstic')
            default_param_col = column_names_lower.get('default_parameter_value') or column_names_lower.get('default parameter value')
            
            if not geo_class_col:
                self.log_message("ERROR: 'Geo-Class' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'Geo-Class' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            if not karstic_col:
                self.log_message("ERROR: 'Karstic' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'Karstic' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            if not default_param_col:
                self.log_message("ERROR: 'Default_Parameter_Value' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'Default_Parameter_Value' column.\nFound columns: {', '.join(df.columns.tolist())}")
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
