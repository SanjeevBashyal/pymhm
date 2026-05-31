# ERA5-Land Module Changelog

## Version 0.1.0 (2025-10-21)

### Initial Release

**New Features:**
- Created ERA5_Land module under pyhydrology package
- Implemented ERA5-Land data downloader with ZIP extraction support
- Implemented ERA5-Land to SWAT format processor with incremental processing
- Added flexible variable name detection for different NetCDF formats
- Added automatic coordinate normalization (precision and dtype)
- Added relative humidity calculation using Magnus formula
- Added wind speed calculation from u/v components

**Module Structure:**
- `__init__.py`: Module initialization and API
- `downloader.py`: Download ERA5-Land data from Copernicus CDS
- `processor.py`: Process NetCDF files to SWAT format
- `config.py`: Configuration settings and constants
- `example_usage.py`: Usage examples
- `README.md`: Documentation
- `CHANGELOG.md`: This file

**Key Capabilities:**

1. **Data Download:**
   - Download from Copernicus Climate Data Store (CDS)
   - Automatic ZIP extraction
   - Configurable time ranges, variables, and spatial extent
   - Resume capability (skips existing files)

2. **Data Processing:**
   - Incremental processing (one month at a time) for memory efficiency
   - Bounding box subsetting to avoid NaN issues
   - Coordinate normalization (float32, 1 decimal precision)
   - Variable name flexibility (handles multiple NetCDF formats)
   - Daily aggregation:
     * Precipitation: daily maximum → mm/day
     * Temperature: daily max/min → °C
     * Solar radiation: daily maximum → J/m²
     * Humidity: daily median → %
     * Wind: daily median → m/s

3. **Output Format:**
   - SWAT-compatible station files
   - stations.cli metadata file
   - Individual data files per variable per station
   - Proper date formatting and missing value handling

**Bug Fixes:**
- Fixed NaN values from vector clipping (replaced with bounding box)
- Fixed humidity DataFrame empty issue (copy tmax structure)
- Fixed coordinate precision mismatch (normalize to float32)
- Fixed coordinate dtype mismatch (float64 vs float32)
- Fixed time coordinate naming inconsistency (time vs valid_time)
- Fixed merge errors from duplicate columns (spatial_ref, number, expver)

**Dependencies:**
- xarray
- pandas
- numpy
- geopandas
- cdsapi
- pyhydrology (NetCDFProcessor, DEMSampler)

**Known Limitations:**
- Requires CDS API key and account
- Processing speed depends on available data files
- Coordinate precision fixed at 0.1 degrees

**Future Enhancements:**
- Add parallel processing support
- Add data quality checks and validation
- Add support for more output formats (beyond SWAT)
- Add visualization tools
- Add unit tests


