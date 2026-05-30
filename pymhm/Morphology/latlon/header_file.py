# -*- coding: utf-8 -*-
"""mHM grid header file writing."""
from ..common import os


class HeaderFileMixin:
    """mHM grid header file writing."""

    def write_header_file(self, L_info, output_path):
        """
        Helper function to write header.txt file for a given L level.
        
        Args:
            L_info: Dictionary with keys: ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value
            output_path: Directory path where header.txt should be created
        """
        header_file = os.path.join(output_path, "header.txt")
        
        try:
            with open(header_file, 'w', encoding='utf-8') as f:
                f.write(f"ncols\t{L_info['ncols']}\n")
                f.write(f"nrows\t{L_info['nrows']}\n")
                f.write(f"xllcorner\t{L_info['xllcorner']}\n")
                f.write(f"yllcorner\t{L_info['yllcorner']}\n")
                f.write(f"cellsize\t{L_info['cellsize']}\n")
                f.write(f"NODATA_value\t{L_info['NODATA_value']}\n")
                f.write("\n")
            
            self.log_message(f"Header file written: {header_file}")
            self.mark_output_prepared(
                header_file,
                name="header.txt",
                loaded=False
            )
            return True
        except Exception as e:
            self.log_message(f"ERROR: Failed to write header file {header_file}: {e}")
            return False
