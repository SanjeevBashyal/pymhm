# -*- coding: utf-8 -*-
"""Soil lookup table parsing."""
from ..common import (
    os,
    QMessageBox,
)


class SoilLookupMixin:
    """Soil lookup table parsing."""

    def load_soil_lookup_table(self):
        """
        Load soil lookup table from CSV file using pandas.
        Maps Dominant_S to CLASS (soil class code).

        Returns:
            dict: Mapping of Dominant_S to CLASS, or None if failed
        """
        lookup_file = self.dialog.lineEdit_soil_lookup.text()

        if not lookup_file or not os.path.exists(lookup_file):
            self.log_message(
                "ERROR: Soil lookup table file not found or not specified.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Soil lookup table file.")
            return None

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
                    "Could not read soil lookup table file. Encoding issue or file is empty.")
                return None
            
            # Log column names found
            self.log_message(f"CSV columns found: {', '.join(df.columns.tolist())}")
            
            # Check if required columns exist (case-insensitive)
            column_names_lower = {col.lower(): col for col in df.columns}
            dominant_s_col = column_names_lower.get('dominant_s')
            class_col = column_names_lower.get('class')
            
            if not dominant_s_col:
                self.log_message("ERROR: 'Dominant_S' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'Dominant_S' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            if not class_col:
                self.log_message("ERROR: 'CLASS' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'CLASS' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            self.log_message(f"Using column '{dominant_s_col}' for Dominant_S and '{class_col}' for CLASS")
            
            # Create lookup mapping
            lookup_mapping = {}
            
            # Remove rows with missing values in required columns
            df_clean = df[[dominant_s_col, class_col]].dropna()
            
            # Process each row
            for idx, row in df_clean.iterrows():
                try:
                    # Get Dominant_S value (strip whitespace and convert to string)
                    dominant_s = str(row[dominant_s_col]).strip() if pd.notna(row[dominant_s_col]) else ''
                    
                    # Get CLASS value
                    class_value = row[class_col]
                    
                    if not dominant_s:
                        self.log_message(f"WARNING: Row {idx + 2} (header + 1-based) has empty Dominant_S value. Skipping.")
                        continue
                    
                    if pd.isna(class_value):
                        self.log_message(f"WARNING: Row {idx + 2} has empty CLASS value. Skipping.")
                        continue
                    
                    # Convert CLASS to integer
                    try:
                        class_code = int(float(class_value))  # Use float first to handle numeric strings
                    except (ValueError, TypeError):
                        self.log_message(
                            f"WARNING: Row {idx + 2} has invalid CLASS value '{class_value}'. Skipping.")
                        continue
                    
                    if class_code > 0:
                        lookup_mapping[dominant_s] = class_code
                        if len(lookup_mapping) <= 10:  # Only log first 10 mappings
                            self.log_message(f"  Mapped: {dominant_s} -> {class_code}")
                    else:
                        self.log_message(
                            f"WARNING: Row {idx + 2} has CLASS <= 0 ({class_code}). Skipping.")
                except Exception as e:
                    self.log_message(
                        f"WARNING: Error parsing lookup table row {idx + 2}: {e}")
                    continue

            self.log_message(f"Processed {len(df)} rows from CSV file.")
            
            if lookup_mapping:
                self.log_message(
                    f"Loaded {len(lookup_mapping)} entries from soil lookup table.")
                return lookup_mapping
            else:
                self.log_message(
                    "ERROR: Soil lookup table is empty or could not be parsed.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"Soil lookup table is empty or could not be parsed.\nProcessed {len(df)} rows but found no valid mappings.")
                return None
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
                f"ERROR reading soil lookup table: {e}\n{error_details}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error reading soil lookup table:\n{str(e)}")
            return None
