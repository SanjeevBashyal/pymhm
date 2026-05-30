# -*- coding: utf-8 -*-
"""Lat/lon NetCDF generation."""
from ..common import (
    os,
    QMessageBox,
)


class LatLonNetcdfMixin:
    """Lat/lon NetCDF generation."""

    def create_latlon_nc_file(self, latlon_folder):
        """
        Create latlon.nc NetCDF file with L0, L1, and L11 information.
        Includes xc/yc coordinates (1D arrays in UTM meters) and lat/lon values (2D arrays) for each level.
        """
        try:
            import numpy as np
            import xarray as xr
            from pyproj import Transformer
            import time
            
            # Get input CRS for coordinate transformation
            input_crs = self.dialog.get_crs()
            if not input_crs.isValid():
                self.log_message("ERROR: Invalid input CRS. Cannot create latlon.nc file.")
                return False
            
            # Get CRS string for pyproj
            crs_string = None
            if input_crs.postgisSrid():
                crs_string = f"EPSG:{input_crs.postgisSrid()}"
            else:
                crs_string = input_crs.authid()
            
            if not crs_string:
                self.log_message("ERROR: Could not determine CRS string. Cannot create latlon.nc file.")
                return False
            
            # Create transformer from input CRS to WGS84 (EPSG:4326)
            transformer = Transformer.from_crs(crs_string, "EPSG:4326", always_xy=True)
            
            # Helper function to calculate xc/yc arrays (1D) for a given level
            def calculate_coordinate_arrays(L_info):
                """
                Calculate xc and yc 1D arrays for grid cell centers in UTM meters.
                xc represents column centers, yc represents row centers.
                """
                ncols = L_info['ncols']
                nrows = L_info['nrows']
                xllcorner = L_info['xllcorner']
                yllcorner = L_info['yllcorner']
                cellsize = L_info['cellsize']
                
                # Calculate xc (column centers) - 1D array
                # xc = xllcorner + (column_index + 0.5) * cellsize
                xc = np.array([xllcorner + (i + 0.5) * cellsize for i in range(ncols)], dtype=np.float32)
                
                # Calculate yc (row centers) - 1D array
                # yc = yllcorner + (row_index + 0.5) * cellsize
                # Note: yllcorner is bottom-left, so we go from bottom to top
                yc = np.array([yllcorner + (i + 0.5) * cellsize for i in range(nrows)], dtype=np.float32)
                
                return xc, yc
            
            # Helper function to calculate lat/lon 2D arrays
            def calculate_latlon_arrays(L_info, xc, yc):
                """
                Calculate lat and lon 2D arrays from xc/yc coordinates.
                Creates a meshgrid and transforms UTM coordinates to lat/lon.
                """
                ncols = L_info['ncols']
                nrows = L_info['nrows']
                
                # Create meshgrid of xc and yc
                x_grid, y_grid = np.meshgrid(xc, yc)
                
                # Transform to lat/lon (transformer expects x, y and returns lon, lat)
                lon_array, lat_array = transformer.transform(x_grid.flatten(), y_grid.flatten())
                
                # Reshape to 2D (nrows, ncols)
                lon_2d = np.array(lon_array, dtype=np.float64).reshape(nrows, ncols)
                lat_2d = np.array(lat_array, dtype=np.float64).reshape(nrows, ncols)
                
                return lat_2d, lon_2d
            
            # Calculate coordinates for each level
            self.log_message("Calculating coordinate arrays for L0...")
            xc_l0, yc_l0 = calculate_coordinate_arrays(self.L0)
            lat_l0, lon_l0 = calculate_latlon_arrays(self.L0, xc_l0, yc_l0)
            
            self.log_message("Calculating coordinate arrays for L1...")
            xc_l1, yc_l1 = calculate_coordinate_arrays(self.L1)
            lat_l1, lon_l1 = calculate_latlon_arrays(self.L1, xc_l1, yc_l1)
            
            self.log_message("Calculating coordinate arrays for L11...")
            xc_l11, yc_l11 = calculate_coordinate_arrays(self.L11)
            lat_l11, lon_l11 = calculate_latlon_arrays(self.L11, xc_l11, yc_l11)
            
            # Create xarray Dataset
            self.log_message("Creating NetCDF dataset...")
            
            # Create dimensions and coordinate variables
            ds = xr.Dataset(
                {
                    # L0 variables
                    'xc_l0': (['xc_l0'], xc_l0),
                    'yc_l0': (['yc_l0'], yc_l0),
                    'lon_l0': (['yc_l0', 'xc_l0'], lon_l0),
                    'lat_l0': (['yc_l0', 'xc_l0'], lat_l0),
                    
                    # L1 variables
                    'xc': (['xc'], xc_l1),
                    'yc': (['yc'], yc_l1),
                    'lon': (['yc', 'xc'], lon_l1),
                    'lat': (['yc', 'xc'], lat_l1),
                    
                    # L11 variables (same as L1)
                    'xc_l11': (['xc_l11'], xc_l11),
                    'yc_l11': (['yc_l11'], yc_l11),
                    'lon_l11': (['yc_l11', 'xc_l11'], lon_l11),
                    'lat_l11': (['yc_l11', 'xc_l11'], lat_l11),
                },
                attrs={
                    'description': 'lat lon file',
                    'projection': crs_string.lower(),
                    'history': f'Created {time.ctime()}'
                }
            )
            
            # Set variable attributes
            # L0 attributes
            ds['xc_l0'].attrs['axis'] = 'X'
            ds['yc_l0'].attrs['axis'] = 'Y'
            ds['lon_l0'].attrs['units'] = 'degrees_east'
            ds['lon_l0'].attrs['long_name'] = 'longitude at level 0'
            ds['lat_l0'].attrs['units'] = 'degrees_north'
            ds['lat_l0'].attrs['long_name'] = 'latitude at level 0'
            
            # L1 attributes
            ds['xc'].attrs['axis'] = 'X'
            ds['yc'].attrs['axis'] = 'Y'
            ds['lon'].attrs['units'] = 'degrees_east'
            ds['lon'].attrs['long_name'] = 'longitude at level 1'
            ds['lat'].attrs['units'] = 'degrees_north'
            ds['lat'].attrs['long_name'] = 'latitude at level 1'
            
            # L11 attributes
            ds['xc_l11'].attrs['axis'] = 'X'
            ds['yc_l11'].attrs['axis'] = 'Y'
            ds['lon_l11'].attrs['units'] = 'degrees_east'
            ds['lon_l11'].attrs['long_name'] = 'longitude at level 11'
            ds['lat_l11'].attrs['units'] = 'degrees_north'
            ds['lat_l11'].attrs['long_name'] = 'latitude at level 11'
            
            # Set encoding for compression
            encoding = {
                'lon_l0': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L0['nrows'], self.L0['ncols'])
                },
                'lat_l0': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L0['nrows'], self.L0['ncols'])
                },
                'lon': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L1['nrows'], self.L1['ncols'])
                },
                'lat': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L1['nrows'], self.L1['ncols'])
                },
                'lon_l11': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L11['nrows'], self.L11['ncols'])
                },
                'lat_l11': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L11['nrows'], self.L11['ncols'])
                },
            }
            
            # Save to NetCDF file
            output_file = os.path.join(latlon_folder, "latlon_1.nc")
            ds.to_netcdf(output_file, encoding=encoding)
            
            self.log_message(f"NetCDF file created successfully: {output_file}")
            self.mark_output_prepared(
                output_file,
                name="latlon_1.nc",
                loaded=False
            )
            return True
            
        except ImportError as e:
            self.log_message(f"ERROR: Required library not available: {e}")
            QMessageBox.warning(
                self.dialog, "Missing Dependency",
                f"Required library not available: {e}\n"
                "Please install: pip install xarray netCDF4 pyproj numpy")
            return False
        except Exception as e:
            import traceback
            error_details = traceback.format_exc()
            self.log_message(f"ERROR creating latlon.nc file: {e}\n{error_details}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error creating latlon.nc file:\n{str(e)}")
            return False
