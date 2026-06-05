# -*- coding: utf-8 -*-
"""mHM soil classdefinition writer."""
from ..common import (
    os,
    QMessageBox,
    morph_folder,
)


class SoilClassificationWriterMixin:
    """mHM soil classdefinition writer."""

    def soil_classdefinition_writer(self):
        """
        Read soil lookup file and produce soil_classdefinition.txt.
        
        Input CSV format:
        CLASS,Dominant_S,SNAM,NLAYERS,...,SOL_Z1,SOL_BD1,CLAY1,SAND1,...,SOL_Z2,SOL_BD2,CLAY2,SAND2,...
        
        Output format:
        nSoil_Types  <total_rows>
        MU_GLOBAL	HORIZON	UD[mm]	LD[mm]	CLAY[%]	SAND[%]	BD[gcm-3]
        1	1	0	300	50.0	25.0	1.75
        1	2	300	1000	50.0	25.0	1.75
        2	1	0	300	19.85	53.45	1.41
        2	2	300	1000	31.33	45.22	1.5
        ...
        """
        self.log_message("\n--- Writing Soil Classification Definition File ---")
        
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        # Get soil lookup file path
        lookup_file = self.dialog.lineEdit_soil_lookup.text()
        
        if not lookup_file or not os.path.exists(lookup_file):
            self.log_message(
                "ERROR: Soil lookup table file not found or not specified.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Soil lookup table file.")
            return
        
        morph_output_folder = morph_folder(self.dialog.project_folder)
        os.makedirs(morph_output_folder, exist_ok=True)
        output_file = os.path.join(
            morph_output_folder, "soil_classdefinition.txt")
        
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
            class_col = column_names_lower.get('class')
            nlayers_col = column_names_lower.get('nlayers')
            
            if not class_col:
                self.log_message("ERROR: 'CLASS' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'CLASS' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            if not nlayers_col:
                self.log_message("ERROR: 'NLAYERS' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'NLAYERS' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            # Check for horizon columns (SOL_Z1, CLAY1, SAND1, SOL_BD1, etc.)
            horizon_columns = {}
            for i in range(1, 11):  # Check horizons 1-10
                sol_z_col = column_names_lower.get(f'sol_z{i}')
                clay_col = column_names_lower.get(f'clay{i}')
                sand_col = column_names_lower.get(f'sand{i}')
                sol_bd_col = column_names_lower.get(f'sol_bd{i}')
                
                if sol_z_col and clay_col and sand_col and sol_bd_col:
                    horizon_columns[i] = {
                        'sol_z': sol_z_col,
                        'clay': clay_col,
                        'sand': sand_col,
                        'sol_bd': sol_bd_col
                    }
            
            if not horizon_columns:
                self.log_message("ERROR: No horizon columns found (SOL_Z1, CLAY1, SAND1, SOL_BD1, etc.).")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    "CSV file must contain horizon columns (SOL_Z1, CLAY1, SAND1, SOL_BD1, etc.).")
                return None
            
            self.log_message(f"Found horizon columns for horizons: {sorted(horizon_columns.keys())}")
            
            # Clean data
            df_clean = df[[class_col, nlayers_col] + 
                         [col for h in horizon_columns.values() for col in h.values()]].dropna(subset=[class_col, nlayers_col])
            
            # Sort by CLASS to ensure proper ordering
            df_clean = df_clean.sort_values(by=class_col)
            
            # Build output rows
            output_rows = []
            
            for idx, row in df_clean.iterrows():
                mu_global = int(float(row[class_col])) if pd.notna(row[class_col]) else 0
                nlayers = int(float(row[nlayers_col])) if pd.notna(row[nlayers_col]) else 1
                
                # Limit nlayers to available horizons
                nlayers = min(nlayers, max(horizon_columns.keys()))
                
                prev_depth = 0  # Upper depth starts at 0 for first horizon
                
                for horizon in range(1, nlayers + 1):
                    if horizon not in horizon_columns:
                        self.log_message(f"WARNING: Horizon {horizon} data not available for CLASS {mu_global}. Skipping.")
                        continue
                    
                    cols = horizon_columns[horizon]
                    
                    # Get values
                    sol_z = float(row[cols['sol_z']]) if pd.notna(row[cols['sol_z']]) else 0
                    clay = float(row[cols['clay']]) if pd.notna(row[cols['clay']]) else 0
                    sand = float(row[cols['sand']]) if pd.notna(row[cols['sand']]) else 0
                    sol_bd = float(row[cols['sol_bd']]) if pd.notna(row[cols['sol_bd']]) else 0
                    
                    # Upper depth is previous horizon's lower depth (or 0 for first horizon)
                    ud = prev_depth
                    # Lower depth is SOL_Z for this horizon
                    ld = sol_z
                    
                    output_rows.append({
                        'mu_global': mu_global,
                        'horizon': horizon,
                        'ud': ud,
                        'ld': ld,
                        'clay': clay,
                        'sand': sand,
                        'bd': sol_bd
                    })
                    
                    # Update previous depth for next horizon
                    prev_depth = ld
            
            if not output_rows:
                self.log_message("ERROR: No valid soil horizon data found.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    "No valid soil horizon data found in the lookup table.")
                return None
            
            n_soil_types = len(output_rows)
            self.log_message(f"Found {n_soil_types} soil horizon entries")
            
            # Write output file
            with open(output_file, 'w', encoding='utf-8') as f:
                # Write header
                f.write(f"nSoil_Types  {n_soil_types}\n")
                f.write("\n")
                f.write("MU_GLOBAL\tHORIZON\tUD[mm]\tLD[mm]\tCLAY[%]\tSAND[%]\tBD[gcm-3]\n")
                
                # Write data rows
                for row_data in output_rows:
                    f.write(f"{row_data['mu_global']}\t"
                           f"{row_data['horizon']}\t"
                           f"{row_data['ud']:.0f}\t"
                           f"{row_data['ld']:.0f}\t"
                           f"{row_data['clay']:.2f}\t"
                           f"{row_data['sand']:.2f}\t"
                           f"{row_data['bd']:.2f}\n")
            
            self.log_message(f"Soil classification definition file written successfully: {output_file}")
            self.mark_output_prepared(
                output_file,
                name="soil_classdefinition.txt",
                loaded=False
            )
            QMessageBox.information(
                self.dialog, "Success",
                f"Soil classification definition file created successfully.\n{output_file}")
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
                f"ERROR writing soil classification file: {e}\n{error_details}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error writing soil classification file:\n{str(e)}")
            return None
