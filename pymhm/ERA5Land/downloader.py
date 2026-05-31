import cdsapi
import os
import shutil
import zipfile
from pathlib import Path

###################################################################################################

url = "https://cds.climate.copernicus.eu/api/"
api_key = "5fa3ac2c-dba0-4cd4-b2b0-3cd3937f9c0b"

###################################################################################################

years = [str(i) for i in range(1995,2025)]                    ## Time range
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
extents = [25.0, 79.0, 31.0, 89.0]                        ## North, West, South, East. Default: global

destination_folder = 'E:/0 Python/pyhydrology/1 Data/ERA5'  # Update with your ERA5 data directory

variables = ['2m_temperature', 'total_precipitation', 'surface_solar_radiation_downwards', '2m_dewpoint_temperature', '10m_u_component_of_wind', '10m_v_component_of_wind']

###################################################################################################

# Initialize the CDS API client
c = cdsapi.Client(key=api_key, url=url)

###################################################################################################

for iv in variables:
    for iy in years:
        for iM in months:
            # Ensure destination directory exists
            Path(destination_folder).mkdir(parents=True, exist_ok=True)

            final_nc = Path(destination_folder) / f'ERA5_Land_{iv}_{iy}_{iM}.nc'
            if final_nc.exists():
                print(f'{final_nc.name} already exists, skipping.')
                continue

            # Download as ZIP (CDS now often serves ZIP with data_0.nc inside)
            zip_name = f'ERA5_Land_{iv}_{iy}_{iM}.nc'
            zip_path = Path(zip_name).resolve()

            print(f'Downloading ZIP for {iv} {iy}-{iM} ...')
            c.retrieve(
                'reanalysis-era5-land',
                {
                    'product_type': 'reanalysis',
                    'format': 'zip',
                    'variable': [iv],
                    'year': [iy],
                    'month': [iM],
                    'day': [
                        '01','02','03',
                        '04','05','06',
                        '07','08','09',
                        '10','11','12',
                        '13','14','15',
                        '16','17','18',
                        '19','20','21',
                        '22','23','24',
                        '25','26','27',
                        '28','29','30',
                        '31'
                    ],
                    'time': [
                        '00:00','01:00','02:00',
                        '03:00','04:00','05:00',
                        '06:00','07:00','08:00',
                        '09:00','10:00','11:00',
                        '12:00','13:00','14:00',
                        '15:00','16:00','17:00',
                        '18:00','19:00','20:00',
                        '21:00','22:00','23:00'
                    ],
                    'area': extents,
                },
                str(zip_path)
            )

            # Extract the file from ZIP (there will be exactly 1 file inside)
            try:
                if not zipfile.is_zipfile(zip_path):
                    raise RuntimeError(f'Downloaded file is not a valid ZIP: {zip_path}')

                with zipfile.ZipFile(zip_path, 'r') as zf:
                    members = zf.namelist()
                    if not members:
                        raise FileNotFoundError('ZIP archive is empty')
                    
                    # Grab the first (and typically only) file in the archive
                    internal_file = members[0]
                    print(f'  Extracting: {internal_file}')
                    
                    # Extract to a temp location next to ZIP
                    extract_dir = zip_path.parent
                    extracted_path = Path(zf.extract(internal_file, path=extract_dir))

                # Move/rename to destination with expected .nc filename
                shutil.move(str(extracted_path), str(final_nc))
                print(f'Created {final_nc}')
            finally:
                # Clean up ZIP after extraction
                if zip_path.exists():
                    try:
                        zip_path.unlink()
                    except Exception:
                        pass

