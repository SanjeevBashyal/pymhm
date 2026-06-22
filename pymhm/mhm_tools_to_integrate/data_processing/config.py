"""
Configuration file for ERA5-Land module

This file contains default configurations and settings for the ERA5-Land module.
Users can modify these settings or provide their own configuration.
"""

from pathlib import Path

# Default CDS API settings
DEFAULT_CDS_URL = "https://cds.climate.copernicus.eu/api/"

# Default geographic extent (North, West, South, East)
DEFAULT_EXTENT = [25.0, 79.0, 31.0, 89.0]

# Default variables to download
DEFAULT_VARIABLES = [
    "2m_temperature",
    "total_precipitation",
    "surface_solar_radiation_downwards",
    "2m_dewpoint_temperature",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
]

# Variable name mappings for different NetCDF formats
VARIABLE_NAME_MAPPINGS = {
    "precipitation": [
        "tp",
        "Total_precipitation_surface_1_Hour_Accumulation",
        "total_precipitation",
    ],
    "temperature": ["t2m", "Temperature_height_above_ground", "2m_temperature"],
    "dewpoint": [
        "d2m",
        "Dewpoint_temperature_height_above_ground",
        "2m_dewpoint_temperature",
    ],
    "solar_radiation": [
        "ssrd",
        "surface_solar_radiation_downwards",
        "Downward_Short-Wave_Radiation_Flux_surface",
    ],
    "u_wind": [
        "u10",
        "10m_u_component_of_wind",
        "u-component_of_wind_height_above_ground",
    ],
    "v_wind": [
        "v10",
        "10m_v_component_of_wind",
        "v-component_of_wind_height_above_ground",
    ],
}

# Default coordinate precision (decimal places)
COORDINATE_PRECISION = 1

# Default data type for coordinates
COORDINATE_DTYPE = "float32"

# Missing value indicator
MISSING_VALUE = -99.0

# Default time range
DEFAULT_START_YEAR = 1995
DEFAULT_END_YEAR = 2025

# Default months (all months)
DEFAULT_MONTHS = [
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
]

# Default days (all days)
DEFAULT_DAYS = [
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "23",
    "24",
    "25",
    "26",
    "27",
    "28",
    "29",
    "30",
    "31",
]

# Default hours (all hours)
DEFAULT_HOURS = [
    "00:00",
    "01:00",
    "02:00",
    "03:00",
    "04:00",
    "05:00",
    "06:00",
    "07:00",
    "08:00",
    "09:00",
    "10:00",
    "11:00",
    "12:00",
    "13:00",
    "14:00",
    "15:00",
    "16:00",
    "17:00",
    "18:00",
    "19:00",
    "20:00",
    "21:00",
    "22:00",
    "23:00",
]

# File naming patterns
NC_FILE_PATTERN = "ERA5_Land_{variable}_{year}_{month}.nc"
STATION_FILE_PATTERNS = {
    "pcp": "station_{id}_pcp.txt",
    "tmp": "station_{id}_tmp.txt",
    "slr": "station_{id}_slr.txt",
    "hmd": "station_{id}_hmd.txt",
    "wnd": "station_{id}_wnd.txt",
}

# Output file names
STATIONS_CLI_FILE = "stations.cli"

# Magnus formula constants for relative humidity calculation
MAGNUS_A = 17.62
MAGNUS_B = 243.12


def get_default_paths(project_root: Path = None) -> dict:
    """
    Get default paths for data directories.

    Args:
        project_root: Root directory of the project. If None, uses current working directory.

    Returns:
        Dictionary with default paths
    """
    if project_root is None:
        project_root = Path.cwd()

    return {
        "nc_folder": project_root / "data_raw" / "meteo" / "ERA5Land",
        "vector_file": project_root / "data_raw" / "morpho" / "watershed.shp",
        "dem_file": project_root / "data_raw" / "morpho" / "SRTM_WGS_84.tif",
        "output_dir": project_root / "outputs",
    }
