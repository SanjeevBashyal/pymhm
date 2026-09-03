# ERA5-Land processing

This package processes ERA5-Land NetCDF files that the user already has.
It does not download data or manage Copernicus CDS credentials.

The mhm_qgis dialog accepts separate precipitation, temperature, and optional PET
folders. It validates each folder, combines compatible files along time,
resamples them to the configured L2 grid, crops them to the common model
extent, and writes mHM forcing files below `mhm-plugin/data/meteo`.

The legacy SWAT conversion helpers remain available through
`process_era5_to_swat`. The mHM conversion helpers remain available through
`inspect_era5_folder` and `process_era5_to_mhm`.
