"""
ERA5-Land Module for pymhm

This module contains tools for downloading, processing, and converting ERA5-Land 
reanalysis data for hydrological modeling applications, particularly SWAT format.

Main Components:
- download: Download ERA5-Land data from Copernicus CDS
- processor: Process and convert ERA5-Land NetCDF files to SWAT input format
- extractor: Extract data from compressed NetCDF archives

Author: pymhm Development Team
"""

from pathlib import Path

from .mhm import MissingDependencyError

__version__ = "0.1.0"
__all__ = [
    "MissingDependencyError",
    "download_era5",
    "inspect_era5_folder",
    "process_era5_to_mhm",
    "process_era5_to_swat",
]

# Module metadata
MODULE_DIR = Path(__file__).parent
SCRIPTS_DIR = MODULE_DIR / "scripts"


# Import main functions when needed
def download_era5(*args, **kwargs):
    """Download ERA5-Land data from CDS API"""
    from .downloader import main as download_main
    return download_main(*args, **kwargs)


def process_era5_to_swat(*args, **kwargs):
    """Process ERA5-Land NetCDF files to SWAT format"""
    from .processor import main as processor_main
    return processor_main(*args, **kwargs)


def inspect_era5_folder(*args, **kwargs):
    """Inspect ERA5-Land files available in a folder."""
    from .mhm import inspect_era5_folder as inspect_main
    return inspect_main(*args, **kwargs)


def process_era5_to_mhm(*args, **kwargs):
    """Process ERA5-Land NetCDF files to mHM forcing format."""
    from .mhm import process_era5_to_mhm as processor_main
    return processor_main(*args, **kwargs)


