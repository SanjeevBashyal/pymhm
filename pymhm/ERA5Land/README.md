# ERA5-Land Module

This module provides tools for downloading and processing ERA5-Land reanalysis data for hydrological modeling.

## Features

- **Download ERA5-Land data** from Copernicus Climate Data Store (CDS)
- **Process NetCDF files** to SWAT model input format
- **Extract meteorological variables**:
  - Precipitation (pcp)
  - Temperature max/min (tmp)
  - Solar radiation (slr)
  - Relative humidity (hmd)
  - Wind speed (wnd)
- **Generate station files** with proper coordinates and elevation
- **Handle multiple NetCDF formats** with flexible variable name detection

## Module Structure

```
ERA5_Land/
├── __init__.py          # Module initialization and main API
├── downloader.py        # Download ERA5-Land data from CDS
├── processor.py         # Process NetCDF to SWAT format
└── README.md           # This file
```

## Usage

### 1. Download ERA5-Land Data

```python
from pyhydrology.ERA5_Land import download_era5

# The downloader will be configured with:
# - CDS API credentials
# - Time range (years and months)
# - Geographic extent
# - Variables to download
# - Destination folder

# Run the downloader script directly:
# python pyhydrology/ERA5_Land/downloader.py
```

**Configuration in downloader.py:**
```python
url = "https://cds.climate.copernicus.eu/api/"
api_key = "your-cds-api-key"
years = [str(i) for i in range(1995,2025)]
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
extents = [25.0, 79.0, 31.0, 89.0]  # North, West, South, East
destination_folder = 'path/to/ERA5/data'
variables = ['2m_temperature', 'total_precipitation', ...]
```

### 2. Process ERA5-Land to SWAT Format

```python
from pathlib import Path
from pyhydrology.ERA5_Land.processor import main

# Set paths
NC_FOLDER = Path("./1 Data/ERA5")
VECTOR_FILE = Path("./1 Data/Watershed/watershed.shp")
DEM_FILE = Path("./1 Data/Watershed/SRTM_WGS_84.tif")
OUTPUT_DIR = Path("./outputs")

# Process with date range
START_DATE = "2015-01-01"  # Format: YYYY-MM-DD
END_DATE = "2019-12-31"
main(NC_FOLDER, VECTOR_FILE, DEM_FILE, OUTPUT_DIR, START_DATE, END_DATE)

# Or process all available data
main(NC_FOLDER, VECTOR_FILE, DEM_FILE, OUTPUT_DIR)
```

## Output Format

The processor generates SWAT-compatible files:

### stations.cli
```
id    name         lat         lon       elev   pcp                  tmp                  slr                 hmd                 wnd
1     station_1    30.5       79.2      150.0  station_1_pcp.txt    station_1_tmp.txt    station_1_slr.txt   station_1_hmd.txt   station_1_wnd.txt
2     station_2    30.6       79.3      165.0  station_2_pcp.txt    station_2_tmp.txt    station_2_slr.txt   station_2_hmd.txt   station_2_wnd.txt
```

### Station Data Files

Each station gets 5 data files:

**Precipitation (station_X_pcp.txt):**
```
20150101
0.0
5.2
12.7
...
```

**Temperature (station_X_tmp.txt):**
```
20150101
25.3,15.2
26.1,14.8
...
```

**Solar Radiation (station_X_slr.txt):**
```
20150101
25847692.0
26123456.0
...
```

**Relative Humidity (station_X_hmd.txt):**
```
20150101
65.2
68.5
...
```

**Wind Speed (station_X_wnd.txt):**
```
20150101
2.45
3.12
...
```

## Key Features

### Incremental Processing
- Processes one time period (month) at a time
- Low memory usage for large datasets
- Progress tracking and error handling

### Flexible Variable Detection
Handles different NetCDF formats with multiple possible variable names:
- Precipitation: `tp`, `Total_precipitation_surface_1_Hour_Accumulation`
- Temperature: `t2m`, `Temperature_height_above_ground`
- And more...

### Coordinate Normalization
- Automatically handles `time` vs `valid_time` coordinate names
- Rounds coordinates to consistent precision
- Converts float64/float32 to consistent dtype

### Robust Error Handling
- Skips missing files gracefully
- Detailed error messages
- Debug output for troubleshooting

## Requirements

- xarray
- pandas
- numpy
- geopandas
- cdsapi
- pyhydrology (NetCDFProcessor, DEMSampler)

## Notes

- **CDS API Key**: You need a Copernicus Climate Data Store account and API key
- **Memory Efficient**: Processes data incrementally to handle large datasets
- **SWAT Compatible**: Output format follows SWAT meteorological input requirements
- **Customizable**: Easy to modify date ranges, variables, and spatial extent

## Author

Pyhydrology Development Team


