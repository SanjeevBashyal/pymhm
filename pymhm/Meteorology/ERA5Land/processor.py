from __future__ import annotations

from pathlib import Path
from typing import Tuple, Optional, List

import numpy as np
import pandas as pd
import xarray as xr

from .dem import DEMSampler
from .netcdf import NetCDFProcessor


def _get_time_periods_from_files(folder: Path, pattern: str) -> List[str]:
    """Extract unique time periods from NetCDF filenames."""
    files = sorted(Path(folder).glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching '{pattern}' found in {folder}")
    
    time_periods = []
    for file_path in files:
        # Extract time period from filename (e.g., "ERA5_Land_total_precipitation_2015_01.nc" -> "2015_01")
        filename = file_path.stem
        parts = filename.split('_')
        if len(parts) >= 4:
            # Look for year_month pattern at the end
            time_part = '_'.join(parts[-2:])  # Last two parts should be year_month
            if time_part not in time_periods:
                time_periods.append(time_part)
    
    return sorted(time_periods)


def _open_single_era_file(folder: Path, pattern: str, time_period: str) -> xr.Dataset:
    """Open a single NetCDF file for a specific time period."""
    # Ensure folder is absolute path
    folder = Path(folder).resolve()
    
    # Find file matching the time period
    files = list(folder.glob(pattern))
    matching_files = [f for f in files if time_period in f.name]
    
    if not matching_files:
        # Try to find files with similar time periods for debugging
        similar_files = [f.name for f in files if time_period.split('_')[0] in f.name][:5]
        raise FileNotFoundError(f"No file found for time period {time_period} with pattern {pattern} in {folder}. Available files for this year: {similar_files}")
    
    target_file = matching_files[0]  # Take the first match
    
    # Try to open the file
    for eng in ("netcdf4", "scipy", "h5netcdf"):
        try:
            ds = xr.open_dataset(str(target_file), engine=eng)
            return ds
        except Exception as e:
            print(f"Warning: Failed to open {target_file} with engine {eng}: {e}")
            continue
    
    raise RuntimeError(f"Failed to open file {target_file} with any engine")


def _open_era_nc_folder(folder: Path, pattern: str) -> xr.Dataset:
    """Legacy function for backward compatibility - opens all files at once."""
    files = sorted(Path(folder).glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching '{pattern}' found in {folder}")

    def try_open_one(path: Path) -> xr.Dataset | None:
        for eng in ("netcdf4", "scipy", "h5netcdf"):
            try:
                return xr.open_dataset(str(path), engine=eng)
            except Exception:
                continue
        return None

    opened = []
    for p in files:
        ds = try_open_one(p)
        if ds is not None:
            opened.append(ds)

    if not opened:
        raise RuntimeError(
            f"Failed to open any files for pattern '{pattern}'. Ensure they are valid NetCDF."
        )

    try:
        merged = xr.combine_by_coords(opened, combine_attrs="override")
    except Exception:
        # Fallback to mfdataset if combine fails
        paths = [str(p) for p in files]
        merged = xr.open_mfdataset(paths, combine="by_coords")
    return merged


def _detect_lat_lon_columns(df: pd.DataFrame) -> Tuple[str, str]:
    if "lat" in df.columns:
        lat_col = "lat"
    elif "latitude" in df.columns:
        lat_col = "latitude"
    else:
        raise KeyError("Latitude column not found in dataframe")

    if "lon" in df.columns:
        lon_col = "lon"
    elif "longitude" in df.columns:
        lon_col = "longitude"
    else:
        raise KeyError("Longitude column not found in dataframe")
    return lat_col, lon_col


def _c_to_k(x: xr.DataArray) -> xr.DataArray:
    return x + 273.15


def _k_to_c(x: xr.DataArray) -> xr.DataArray:
    return x - 273.15


def _compute_rh_percent_from_t_tdew_k(t_k: xr.DataArray, td_k: xr.DataArray) -> xr.DataArray:
    """
    Calculate relative humidity using Magnus formula.
    
    Formula: RH(%) = exp((17.62 * Td) / (243.12 + Td) - (17.62 * Ta) / (243.12 + Ta)) * 100
    
    Args:
        t_k: Air temperature in Kelvin
        td_k: Dewpoint temperature in Kelvin
        
    Returns:
        Relative humidity in percentage (0-100)
    """
    # Convert to Celsius
    t_c = _k_to_c(t_k)
    td_c = _k_to_c(td_k)
    
    # Magnus formula
    num = 17.62 * td_c / (243.12 + td_c)
    den = 17.62 * t_c / (243.12 + t_c)
    rh = np.exp(num - den) * 100.0
    
    # Limit RH between 0 and 100
    return rh.clip(min=0.0, max=100.0)


def _process_single_time_period(nc_folder: Path, vector_path: Path, dem_path: Path, output_dir: Path, 
                               time_period: str, station_data: dict, station_coords: Optional[pd.DataFrame]) -> None:
    """Process a single time period and update station files incrementally."""
    print(f"Processing time period: {time_period}")
    
    # Open datasets for this time period only
    try:
        print(f"  Opening precipitation file for {time_period}...")
        ds_pcp = _open_single_era_file(nc_folder, "ERA5_Land_total_precipitation_*.nc", time_period)
        
        print(f"  Opening temperature file for {time_period}...")
        ds_tmp = _open_single_era_file(nc_folder, "ERA5_Land_2m_temperature_*.nc", time_period)
        
        print(f"  Opening dewpoint file for {time_period}...")
        ds_dew = _open_single_era_file(nc_folder, "ERA5_Land_2m_dewpoint_temperature_*.nc", time_period)
        
        print(f"  Opening solar radiation file for {time_period}...")
        ds_rad = _open_single_era_file(nc_folder, "ERA5_Land_surface_solar_radiation_downwards_*.nc", time_period)
        
        print(f"  Opening u-wind file for {time_period}...")
        ds_u = _open_single_era_file(nc_folder, "ERA5_Land_10m_u_component_of_wind*.nc", time_period)
        
        print(f"  Opening v-wind file for {time_period}...")
        ds_v = _open_single_era_file(nc_folder, "ERA5_Land_10m_v_component_of_wind*.nc", time_period)
        
    except FileNotFoundError as e:
        print(f"  Error: {e}")
        print(f"  Skipping time period {time_period}")
        return
    
    # Standardize time coordinate names
    if 'valid_time' in ds_pcp.coords:
        ds_pcp = ds_pcp.rename({'valid_time': 'time'})
    if 'valid_time' in ds_tmp.coords:
        ds_tmp = ds_tmp.rename({'valid_time': 'time'})
    if 'valid_time' in ds_dew.coords:
        ds_dew = ds_dew.rename({'valid_time': 'time'})
    if 'valid_time' in ds_rad.coords:
        ds_rad = ds_rad.rename({'valid_time': 'time'})
    if 'valid_time' in ds_u.coords:
        ds_u = ds_u.rename({'valid_time': 'time'})
    if 'valid_time' in ds_v.coords:
        ds_v = ds_v.rename({'valid_time': 'time'})

    # FIX 1: Use simple bounding box subsetting instead of vector clipping to avoid NaN issues
    # Get bounding box from vector
    import geopandas as gpd
    gdf = gpd.read_file(vector_path)
    bounds = gdf.total_bounds  # minx, miny, maxx, maxy
    
    # Subset datasets to bounding box (avoid NaN issues from complex vector clipping)
    ds_pcp_clip = ds_pcp.sel(
        latitude=slice(max(bounds[3], bounds[1]), min(bounds[3], bounds[1])),
        longitude=slice(min(bounds[0], bounds[2]), max(bounds[0], bounds[2]))
    )
    ds_tmp_clip = ds_tmp.sel(
        latitude=slice(max(bounds[3], bounds[1]), min(bounds[3], bounds[1])),
        longitude=slice(min(bounds[0], bounds[2]), max(bounds[0], bounds[2]))
    )
    ds_dew_clip = ds_dew.sel(
        latitude=slice(max(bounds[3], bounds[1]), min(bounds[3], bounds[1])),
        longitude=slice(min(bounds[0], bounds[2]), max(bounds[0], bounds[2]))
    )
    ds_rad_clip = ds_rad.sel(
        latitude=slice(max(bounds[3], bounds[1]), min(bounds[3], bounds[1])),
        longitude=slice(min(bounds[0], bounds[2]), max(bounds[0], bounds[2]))
    )
    ds_u_clip = ds_u.sel(
        latitude=slice(max(bounds[3], bounds[1]), min(bounds[3], bounds[1])),
        longitude=slice(min(bounds[0], bounds[2]), max(bounds[0], bounds[2]))
    )
    ds_v_clip = ds_v.sel(
        latitude=slice(max(bounds[3], bounds[1]), min(bounds[3], bounds[1])),
        longitude=slice(min(bounds[0], bounds[2]), max(bounds[0], bounds[2]))
    )

    # FIX 4: Flexible variable name detection
    # Define possible variable names for each type
    pcp_var_names = ["tp", "Total_precipitation_surface_1_Hour_Accumulation", "total_precipitation"]
    tmp_var_names = ["t2m", "Temperature_height_above_ground", "2m_temperature"]
    dew_var_names = ["d2m", "Dewpoint_temperature_height_above_ground", "2m_dewpoint_temperature"]
    rad_var_names = ["ssrd", "surface_solar_radiation_downwards", "Downward_Short-Wave_Radiation_Flux_surface"]
    u_var_names = ["u10", "10m_u_component_of_wind", "u-component_of_wind_height_above_ground"]
    v_var_names = ["v10", "10m_v_component_of_wind", "v-component_of_wind_height_above_ground"]
    
    # Find actual variable names
    pcp_var = next((v for v in pcp_var_names if v in ds_pcp_clip.data_vars), None)
    tmp_var = next((v for v in tmp_var_names if v in ds_tmp_clip.data_vars), None)
    dew_var = next((v for v in dew_var_names if v in ds_dew_clip.data_vars), None)
    rad_var = next((v for v in rad_var_names if v in ds_rad_clip.data_vars), None)
    u_var = next((v for v in u_var_names if v in ds_u_clip.data_vars), None)
    v_var = next((v for v in v_var_names if v in ds_v_clip.data_vars), None)
    
    if not pcp_var:
        print(f"Warning: No precipitation variable found for {time_period}. Available: {list(ds_pcp_clip.data_vars)}")
        return
    if not tmp_var:
        print(f"Warning: No temperature variable found for {time_period}. Available: {list(ds_tmp_clip.data_vars)}")
        return
    if not dew_var:
        print(f"Warning: No dewpoint variable found for {time_period}. Available: {list(ds_dew_clip.data_vars)}")
        return
    if not rad_var:
        print(f"Warning: No radiation variable found for {time_period}. Available: {list(ds_rad_clip.data_vars)}")
        return
    if not u_var or not v_var:
        print(f"Warning: Wind components not found for {time_period}")
        return

    # Process data for this time period
    # Precipitation: hourly to daily sum (m/day) -> convert to mm/day
    pcp_daily = (ds_pcp_clip[pcp_var].resample(time="1D").max()) * 1000.0

    # Temperature: hourly to daily max/min; convert to Celsius from Kelvin
    tmax_daily_k = ds_tmp_clip[tmp_var].resample(time="1D").max()
    tmin_daily_k = ds_tmp_clip[tmp_var].resample(time="1D").min()
    tmax_daily_c = _k_to_c(tmax_daily_k)
    tmin_daily_c = _k_to_c(tmin_daily_k)

    # Solar radiation: hourly to daily maximum (J/m²) -> Convert to MJ/m²
    slr_daily = (ds_rad_clip[rad_var].resample(time="1D").max()) / 1000000.0

    # Dewpoint temperature: hourly to daily median for RH calculation
    dew_daily_k = ds_dew_clip[dew_var].resample(time="1D").median()
    
    # FIX 2 & 3: Ensure all datasets have same coordinate grid before computing humidity
    # Round coordinates in xarray before any operations to ensure consistency
    def _round_xarray_coords(da, decimals=1):
        """Round latitude and longitude coordinates in xarray DataArray"""
        if 'latitude' in da.coords:
            da = da.assign_coords(latitude=da.latitude.values.round(decimals))
        if 'longitude' in da.coords:
            da = da.assign_coords(longitude=da.longitude.values.round(decimals))
        return da
    
    # Round all daily datasets to same precision BEFORE computing derived variables
    pcp_daily = _round_xarray_coords(pcp_daily)
    tmax_daily_c = _round_xarray_coords(tmax_daily_c)
    tmin_daily_c = _round_xarray_coords(tmin_daily_c)
    slr_daily = _round_xarray_coords(slr_daily)
    dew_daily_k = _round_xarray_coords(dew_daily_k)
    
    # Now compute humidity with aligned coordinates
    # Ensure all inputs have rounded coordinates
    tmax_daily_k = _round_xarray_coords(tmax_daily_k)
    tmin_daily_k = _round_xarray_coords(tmin_daily_k)
    tmean_daily_k = (tmax_daily_k + tmin_daily_k) / 2
    
    # Force align dewpoint to temperature grid using nearest neighbor to handle slight mismatches
    dew_daily_k = _round_xarray_coords(dew_daily_k)
    dew_daily_k = dew_daily_k.reindex_like(tmean_daily_k, method='nearest', tolerance=1e-2)
    
    hmd_values = _compute_rh_percent_from_t_tdew_k(tmean_daily_k, dew_daily_k)
    
    # Wind components: calculate resultant wind velocity and get daily median
    wind_speed = np.sqrt(ds_u_clip[u_var]**2 + ds_v_clip[v_var]**2)
    wnd_daily = wind_speed.resample(time="1D").median()
    wnd_daily = _round_xarray_coords(wnd_daily)

    # Convert to DataFrames - now all should have aligned coordinates
    df_pcp = pcp_daily.to_dataframe(name="pcp").reset_index()
    df_tmax = tmax_daily_c.to_dataframe(name="tmax").reset_index()
    df_tmin = tmin_daily_c.to_dataframe(name="tmin").reset_index()
    df_slr = slr_daily.to_dataframe(name="slr").reset_index()
    
    # FIX 2: Create hmd DataFrame directly from aligned xarray
    # Ensure humidity has exactly the same coordinates as temperature
    hmd_values = hmd_values.reindex_like(tmax_daily_c)
    df_hmd = hmd_values.to_dataframe(name='hmd').reset_index()
    
    df_wnd = wnd_daily.to_dataframe(name="wnd").reset_index()

    # Detect lat/lon column names
    lat_col, lon_col = _detect_lat_lon_columns(df_pcp)

    # FIX 3: Ensure consistent data types AND precision for coordinates
    # Convert to same dtype (float32) and round to same precision
    def _normalize_coords(df, lat_col, lon_col, decimals=1):
        """Convert coordinates to float32 and round to consistent precision"""
        if lat_col in df.columns:
            df[lat_col] = df[lat_col].round(decimals).astype('float32')
        if lon_col in df.columns:
            df[lon_col] = df[lon_col].round(decimals).astype('float32')
        return df
    
    df_pcp = _normalize_coords(df_pcp, lat_col, lon_col)
    df_tmax = _normalize_coords(df_tmax, lat_col, lon_col)
    df_tmin = _normalize_coords(df_tmin, lat_col, lon_col)
    df_slr = _normalize_coords(df_slr, lat_col, lon_col)
    df_hmd = _normalize_coords(df_hmd, lat_col, lon_col)
    df_wnd = _normalize_coords(df_wnd, lat_col, lon_col)

    # Keep only essential columns for each dataframe to avoid duplicate coord columns
    def _prune(df, value_col):
        cols = [c for c in df.columns]
        # Drop common non-spatial coord columns if present
        drop_candidates = ["number", "expver", "spatial_ref"]
        for dc in drop_candidates:
            if dc in df.columns:
                df = df.drop(columns=[dc])
        # Ensure we only retain time, lat, lon, and the value column
        keep = [c for c in ["time", lat_col, lon_col, value_col] if c in df.columns]
        return df[keep]

    df_pcp = _prune(df_pcp, "pcp")
    df_tmax = _prune(df_tmax, "tmax")
    df_tmin = _prune(df_tmin, "tmin")
    df_slr = _prune(df_slr, "slr")
    df_hmd = _prune(df_hmd, "hmd")
    df_wnd = _prune(df_wnd, "wnd")

    # Debug: Print coordinate info before merging
    print(f"  Coordinate check:")
    print(f"    df_pcp: {len(df_pcp)} rows, lat dtype={df_pcp[lat_col].dtype}, lon dtype={df_pcp[lon_col].dtype}")
    print(f"    df_slr: {len(df_slr)} rows, lat dtype={df_slr[lat_col].dtype}, lon dtype={df_slr[lon_col].dtype}")
    print(f"    df_hmd: {len(df_hmd)} rows, lat dtype={df_hmd[lat_col].dtype}, lon dtype={df_hmd[lon_col].dtype}")
    print(f"    df_hmd has {df_hmd['hmd'].notna().sum()} non-null hmd values out of {len(df_hmd)}")
    print(f"    Sample pcp coords: lat={df_pcp[lat_col].iloc[0]}, lon={df_pcp[lon_col].iloc[0]}")
    print(f"    Sample slr coords: lat={df_slr[lat_col].iloc[0]}, lon={df_slr[lon_col].iloc[0]}")
    print(f"    Sample hmd value: {df_hmd['hmd'].iloc[0]:.2f}%")
    
    # Merge all variables on lat/lon/time (inner join to align dates)
    df_tmp = pd.merge(df_tmax, df_tmin, on=["time", lat_col, lon_col], how="inner")
    print(f"  After merging tmax+tmin: {len(df_tmp)} rows")
    
    df = pd.merge(df_pcp, df_tmp, on=["time", lat_col, lon_col], how="inner")
    print(f"  After merging pcp+tmp: {len(df)} rows")
    
    df = pd.merge(df, df_slr, on=["time", lat_col, lon_col], how="inner")
    print(f"  After merging +slr: {len(df)} rows")
    
    df = pd.merge(df, df_hmd, on=["time", lat_col, lon_col], how="inner")
    print(f"  After merging +hmd: {len(df)} rows")
    
    df = pd.merge(df, df_wnd, on=["time", lat_col, lon_col], how="inner")
    print(f"  After merging +wnd: {len(df)} rows (final)")

    # Normalize longitudes to -180..180 for readability
    if df[lon_col].max() > 180:
        df[lon_col] = df[lon_col].where(df[lon_col] <= 180, df[lon_col] - 360)

    # Update station files for this time period
    for station_id, ((lat, lon), group) in enumerate(df.groupby([lat_col, lon_col]), start=1):
        station_name = f"station_{station_id}"
        
        # Initialize station data if not exists
        if station_id not in station_data:
            station_data[station_id] = {
                'name': station_name,
                'lat': float(lat),
                'lon': float(lon),
                # Elevation will be populated later once station_coords is available
                'elev': -99.0,
                'pcp_data': [],
                'tmp_data': [],
                'slr_data': [],
                'hmd_data': [],
                'wnd_data': []
            }
        
        # Sort by time and append data
        group_sorted = group.sort_values(by="time")
        
        # Append data to station records
        station_data[station_id]['pcp_data'].extend(group_sorted["pcp"].fillna(-99.0).tolist())
        station_data[station_id]['tmp_data'].extend(list(zip(group_sorted["tmax"].fillna(-99.0).tolist(), group_sorted["tmin"].fillna(-99.0).tolist())))
        station_data[station_id]['slr_data'].extend(group_sorted["slr"].fillna(-99.0).tolist())
        station_data[station_id]['hmd_data'].extend(group_sorted["hmd"].fillna(-99.0).tolist())
        station_data[station_id]['wnd_data'].extend(group_sorted["wnd"].fillna(-99.0).tolist())
    
    # Clean up memory
    del ds_pcp, ds_tmp, ds_dew, ds_rad, ds_u, ds_v
    del ds_pcp_clip, ds_tmp_clip, ds_dew_clip, ds_rad_clip, ds_u_clip, ds_v_clip
    del df_pcp, df_tmax, df_tmin, df_slr, df_hmd, df_wnd, df


def main(nc_folder: Path, vector_path: Path, dem_path: Path, output_dir: Path, start_date: Optional[str] = None, end_date: Optional[str] = None) -> None:
    """Main function that processes ERA5-Land data incrementally by time period."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Starting incremental ERA5-Land data processing...")
    
    # Get all available time periods from precipitation files (as reference)
    try:
        print(f"Looking for files in: {Path(nc_folder).resolve()}")
        time_periods = _get_time_periods_from_files(nc_folder, "ERA5_Land_total_precipitation_*.nc")
        print(f"Found {len(time_periods)} time periods: {time_periods[:5]}...")
    except FileNotFoundError as e:
        print(f"No precipitation files found: {e}")
        print(f"Please check the folder path: {Path(nc_folder).resolve()}")
        return
    
    # Filter time periods based on date range if specified
    if start_date is not None or end_date is not None:
        filtered_periods = []
        for period in time_periods:
            # Convert period (e.g., "2015_01") to datetime for comparison
            try:
                year, month = period.split('_')
                period_date = pd.to_datetime(f"{year}-{month}-01")
                
                include_period = True
                if start_date is not None:
                    start_dt = pd.to_datetime(start_date)
                    if period_date < start_dt:
                        include_period = False
                
                if end_date is not None and include_period:
                    end_dt = pd.to_datetime(end_date)
                    if period_date > end_dt:
                        include_period = False
                
                if include_period:
                    filtered_periods.append(period)
            except ValueError:
                print(f"Warning: Could not parse time period {period}")
                continue
        
        time_periods = filtered_periods
        print(f"Filtered to {len(time_periods)} time periods based on date range")
    
    if not time_periods:
        print("No time periods found after filtering.")
        return
    
    # Initialize station data storage
    station_data = {}
    station_coords = None
    
    # Process each time period incrementally
    for i, time_period in enumerate(time_periods, 1):
        print(f"Processing {i}/{len(time_periods)}: {time_period}")
        
        try:
            _process_single_time_period(nc_folder, vector_path, dem_path, output_dir, 
                                     time_period, station_data, station_coords)
            
            # Get station coordinates from first successful processing
            if station_coords is None and station_data:
                # Create station coordinates DataFrame
                coords_data = []
                for station_id, data in station_data.items():
                    coords_data.append({
                        'station_id': station_id,
                        'lat': data['lat'],
                        'lon': data['lon']
                    })
                station_coords = pd.DataFrame(coords_data)
                
                # Sample elevations for all stations
                sampler = DEMSampler(str(dem_path))
                station_coords["elevation_m"] = sampler.sample_many(
                    station_coords['lat'].tolist(),
                    station_coords['lon'].tolist(),
                )
                
                # Update station data with elevations
                for station_id, data in station_data.items():
                    elev_val = station_coords.loc[station_coords['station_id'] == station_id, 'elevation_m'].iloc[0]
                    data['elev'] = float(elev_val) if pd.notna(elev_val) else -99.0
                    
        except Exception as e:
            print(f"Error processing {time_period}: {e}")
            continue
    
    if not station_data:
        print("No station data was processed successfully.")
        return
    
    print(f"Writing station files for {len(station_data)} stations...")
    
    # Write station files
    station_records = {
        "pcp": [],
        "tmp": [],
        "slr": [],
        "hmd": [],
        "wnd": []
    }
    
    def write_swatplus_station_file(filepath: Path, station_name: str, var_type: str, 
                                  lat: float, lon: float, elev: float, 
                                  data: list, start_date_str: str) -> None:
        """
        Write a station file in SWAT+ format.
        
        Args:
            filepath: Path to output file
            station_name: Name of the station
            var_type: Type of variable ('pcp', 'tmp', 'slr', 'hmd', 'wnd')
            lat: Latitude
            lon: Longitude
            elev: Elevation
            data: List of data values. For 'tmp', each item is a tuple (max, min).
            start_date_str: Start date string (YYYYMMDD) - used to calculate years
        """
        
        # Parse start date to track years
        start_dt = pd.to_datetime(start_date_str, format="%Y%m%d")
        current_year = start_dt.year
        
        # Calculate number of years (approximate based on data length)
        # This is just for the header, exact count might need full traversal if strictly required,
        # but usually it's the span of data.
        # For simplicity, we'll write the header after processing or just use a placeholder if needed.
        # SWAT+ seems to want nbyr. Let's calculate it properly.
        
        # We need to iterate through data to format lines and count years
        lines = []
        
        current_date = start_dt
        years_seen = set()
        
        for val in data:
            year = current_date.year
            doy = current_date.dayofyear
            years_seen.add(year)
            
            if var_type == 'tmp':
                # val is (max, min)
                tmax, tmin = val
                lines.append(f"{year:<4} {doy:<4} {float(tmax):10.5f} {float(tmin):10.5f}")
            else:
                # val is single value
                lines.append(f"{year:<4} {doy:<4} {float(val):10.5f}")
            
            current_date += pd.Timedelta(days=1)
            
        nbyr = len(years_seen)
        tstep = 0 # Daily
        
        with filepath.open("w", newline="\n") as f:
            # Line 1: Title
            # Format: [filename]: [Description]
            # We'll use a generic description
            desc_map = {
                'pcp': 'Precipitation data',
                'tmp': 'Temperature data',
                'slr': 'Solar radiation data',
                'hmd': 'Relative humidity data',
                'wnd': 'Wind speed data'
            }
            f.write(f"{filepath.name}: {desc_map.get(var_type, 'Data')} - generated by pymhm\n")
            
            # Line 2: Header
            f.write(f"nbyr     tstep       lat       lon      elev\n")
            
            # Line 3: Header values
            # Format matches sample:   24         0    11.600    37.360  1790.000
            f.write(f"{nbyr:4} {tstep:9} {lat:9.3f} {lon:9.3f} {elev:9.3f}\n")
            
            # Data lines
            for line in lines:
                f.write(f"{line}\n")

    for station_id, data in station_data.items():
        station_name = data['name']
        lat = data['lat']
        lon = data['lon']
        elev = data.get('elev', -99.0)
        
        # Define filenames
        pcp_filename = f"{station_name}.pcp"
        tmp_filename = f"{station_name}.tem" # Note .tem extension
        slr_filename = f"{station_name}.slr"
        hmd_filename = f"{station_name}.hmd"
        wnd_filename = f"{station_name}.wnd"
        
        # Get start date from first time period
        if data['pcp_data']:
            first_period = time_periods[0]
            year, month = first_period.split('_')
            start_date_str = f"{year}{month}01"
            
            # Write files
            write_swatplus_station_file(output_dir / pcp_filename, station_name, 'pcp', lat, lon, elev, data['pcp_data'], start_date_str)
            station_records['pcp'].append(pcp_filename)
            
            write_swatplus_station_file(output_dir / tmp_filename, station_name, 'tmp', lat, lon, elev, data['tmp_data'], start_date_str)
            station_records['tmp'].append(tmp_filename)
            
            write_swatplus_station_file(output_dir / slr_filename, station_name, 'slr', lat, lon, elev, data['slr_data'], start_date_str)
            station_records['slr'].append(slr_filename)
            
            write_swatplus_station_file(output_dir / hmd_filename, station_name, 'hmd', lat, lon, elev, data['hmd_data'], start_date_str)
            station_records['hmd'].append(hmd_filename)
            
            write_swatplus_station_file(output_dir / wnd_filename, station_name, 'wnd', lat, lon, elev, data['wnd_data'], start_date_str)
            station_records['wnd'].append(wnd_filename)
    
    # Write entry files (.cli)
    # Format:
    # [filename]: [description]
    # filename
    # [list of files]
    
    cli_configs = [
        ('pcp.cli', 'pcp', 'Precipitation file names'),
        ('tmp.cli', 'tmp', 'Temperature file names'),
        ('slr.cli', 'slr', 'Solar radiation file names'),
        ('hmd.cli', 'hmd', 'Relative humidity file names'),
        ('wnd.cli', 'wnd', 'Wind speed file names')
    ]
    
    for cli_file, key, desc in cli_configs:
        with (output_dir / cli_file).open("w", newline="\n") as f:
            f.write(f"{cli_file}: {desc}\n")
            f.write("filename\n")
            for fname in sorted(station_records[key]):
                f.write(f"{fname}\n")
    
    print(f"Processing complete! Generated files for {len(station_data)} stations.")
    print(f"Output directory: {output_dir}")



if __name__ == "__main__":
    # Get the script directory and project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parents[1]  # Go up two levels from scripts/ to project root
    
    # Defaults for direct execution - use absolute paths
    NC_FOLDER = project_root / "1 Data" / "ERA5"
    VECTOR_FILE = project_root / "1 Data" / "Watershed" / "watershed.shp"
    DEM_FILE = project_root / "1 Data" / "Watershed" / "SRTM_WGS_84.tif"
    OUTPUT_DIR = project_root / "outputs"
    
    print(f"Script directory: {script_dir}")
    print(f"Project root: {project_root}")
    print(f"NC folder: {NC_FOLDER}")
    print(f"NC folder exists: {NC_FOLDER.exists()}")
    
    # Example usage with date range (optional)
    START_DATE = "2000-01-01"  # Format: YYYY-MM-DD
    END_DATE = "2015-12-31"    # Format: YYYY-MM-DD
    main(NC_FOLDER, VECTOR_FILE, DEM_FILE, OUTPUT_DIR, START_DATE, END_DATE)
    
    # Default execution without date filtering
    # main(NC_FOLDER, VECTOR_FILE, DEM_FILE, OUTPUT_DIR)
