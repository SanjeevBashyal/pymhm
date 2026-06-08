from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

import cdsapi

try:
    from ...project_layout import raw_meteo_folder
except ImportError:  # pragma: no cover - direct script execution

    def raw_meteo_folder(project_folder):
        return Path(project_folder) / "data_raw" / "meteo" / "ERA5Land"


DEFAULT_URL = "https://cds.climate.copernicus.eu/api/"
DEFAULT_API_KEY = "5fa3ac2c-dba0-4cd4-b2b0-3cd3937f9c0b"
DEFAULT_YEARS = [str(i) for i in range(1995, 2025)]
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
DEFAULT_EXTENTS = [25.0, 79.0, 31.0, 89.0]  # North, West, South, East
DEFAULT_VARIABLES = [
    "2m_temperature",
    "total_precipitation",
    "surface_solar_radiation_downwards",
    "2m_dewpoint_temperature",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
]
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
DEFAULT_TIMES = [
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


def destination_folder_for_project(project_folder) -> Path:
    """Return the selected project's raw meteorology download folder."""
    return Path(raw_meteo_folder(project_folder))


def main(
    project_folder=None,
    *,
    destination_folder=None,
    years=None,
    months=None,
    extents=None,
    variables=None,
    api_key=DEFAULT_API_KEY,
    url=DEFAULT_URL,
) -> Path:
    """Download ERA5-Land monthly NetCDF files to project/data_raw/meteo."""
    if destination_folder is None:
        if project_folder is None:
            raise ValueError(
                "project_folder is required when destination_folder is not set."
            )
        destination = destination_folder_for_project(project_folder)
    else:
        destination = Path(destination_folder)

    destination.mkdir(parents=True, exist_ok=True)

    years = list(years or DEFAULT_YEARS)
    months = list(months or DEFAULT_MONTHS)
    extents = list(extents or DEFAULT_EXTENTS)
    variables = list(variables or DEFAULT_VARIABLES)

    client = cdsapi.Client(key=api_key, url=url)

    for variable in variables:
        for year in years:
            for month in months:
                final_nc = destination / f"ERA5_Land_{variable}_{year}_{month}.nc"
                if final_nc.exists():
                    print(f"{final_nc.name} already exists, skipping.")
                    continue

                zip_path = destination / f"ERA5_Land_{variable}_{year}_{month}.zip"
                print(f"Downloading ZIP for {variable} {year}-{month} ...")
                client.retrieve(
                    "reanalysis-era5-land",
                    {
                        "product_type": "reanalysis",
                        "format": "zip",
                        "variable": [variable],
                        "year": [year],
                        "month": [month],
                        "day": DEFAULT_DAYS,
                        "time": DEFAULT_TIMES,
                        "area": extents,
                    },
                    str(zip_path),
                )

                try:
                    if not zipfile.is_zipfile(zip_path):
                        raise RuntimeError(
                            f"Downloaded file is not a valid ZIP: {zip_path}"
                        )

                    with zipfile.ZipFile(zip_path, "r") as archive:
                        members = archive.namelist()
                        if not members:
                            raise FileNotFoundError("ZIP archive is empty")

                        internal_file = members[0]
                        print(f"  Extracting: {internal_file}")
                        extracted_path = Path(
                            archive.extract(internal_file, path=destination)
                        )

                    shutil.move(str(extracted_path), str(final_nc))
                    print(f"Created {final_nc}")
                finally:
                    if zip_path.exists():
                        zip_path.unlink()

    return destination


if __name__ == "__main__":  # pragma: no cover
    import argparse

    parser = argparse.ArgumentParser(
        description="Download ERA5-Land files into a pymhm project."
    )
    parser.add_argument(
        "project_folder",
        help="Selected pymhm project folder; downloads go to data_raw/meteo.",
    )
    args = parser.parse_args()
    main(project_folder=args.project_folder)
