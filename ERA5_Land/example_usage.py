"""
Example usage of the ERA5-Land module

This script demonstrates how to use the ERA5-Land module to:
1. Download ERA5-Land data from CDS
2. Process NetCDF files to SWAT format
"""

from pathlib import Path
from pyhydrology.ERA5_Land.processor import main as process_era5

# Example 1: Process ERA5-Land data to SWAT format
# -------------------------------------------------

def example_process_era5_to_swat():
    """Example: Process ERA5-Land NetCDF files to SWAT format"""
    
    # Define paths
    project_root = Path(__file__).resolve().parents[2]
    
    NC_FOLDER = project_root / "1 Data" / "ERA5"
    VECTOR_FILE = project_root / "1 Data" / "Watershed" / "watershed.shp"
    DEM_FILE = project_root / "1 Data" / "Watershed" / "SRTM_WGS_84.tif"
    OUTPUT_DIR = project_root / "outputs"
    
    # Option 1: Process specific date range
    START_DATE = "2015-01-01"  # Format: YYYY-MM-DD
    END_DATE = "2019-12-31"
    
    print("Processing ERA5-Land data to SWAT format...")
    print(f"  NC Folder: {NC_FOLDER}")
    print(f"  Date Range: {START_DATE} to {END_DATE}")
    print(f"  Output: {OUTPUT_DIR}")
    
    process_era5(
        nc_folder=NC_FOLDER,
        vector_path=VECTOR_FILE,
        dem_path=DEM_FILE,
        output_dir=OUTPUT_DIR,
        start_date=START_DATE,
        end_date=END_DATE
    )
    
    print(f"\nProcessing complete! Check {OUTPUT_DIR} for results.")
    print("\nGenerated files:")
    print("  - stations.cli: Station metadata with file references")
    print("  - station_X_pcp.txt: Precipitation data")
    print("  - station_X_tmp.txt: Temperature (max, min)")
    print("  - station_X_slr.txt: Solar radiation")
    print("  - station_X_hmd.txt: Relative humidity")
    print("  - station_X_wnd.txt: Wind speed")


def example_process_all_data():
    """Example: Process all available ERA5-Land data"""
    
    project_root = Path(__file__).resolve().parents[2]
    
    NC_FOLDER = project_root / "1 Data" / "ERA5"
    VECTOR_FILE = project_root / "1 Data" / "Watershed" / "watershed.shp"
    DEM_FILE = project_root / "1 Data" / "Watershed" / "SRTM_WGS_84.tif"
    OUTPUT_DIR = project_root / "outputs"
    
    print("Processing ALL available ERA5-Land data to SWAT format...")
    
    # Process without date filtering - will process all available files
    process_era5(
        nc_folder=NC_FOLDER,
        vector_path=VECTOR_FILE,
        dem_path=DEM_FILE,
        output_dir=OUTPUT_DIR
    )
    
    print("\nProcessing complete!")


# Example 2: Download ERA5-Land data
# -----------------------------------
# Note: You need to edit downloader.py to set your CDS API key and preferences
# Then run: python pyhydrology/ERA5_Land/downloader.py


if __name__ == "__main__":
    # Run the example
    print("=" * 60)
    print("ERA5-Land Module - Example Usage")
    print("=" * 60)
    print()
    
    # Uncomment the example you want to run:
    
    # Example 1: Process specific date range
    example_process_era5_to_swat()
    
    # Example 2: Process all available data
    # example_process_all_data()

