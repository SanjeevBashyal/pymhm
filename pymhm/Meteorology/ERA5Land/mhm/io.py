"""NetCDF and mHM header writing."""
from __future__ import annotations

from pathlib import Path


def write_netcdf(ds, variable: str, output_file: Path) -> None:
    """Write one mHM forcing NetCDF with conservative CF/xarray encoding."""
    output_file.parent.mkdir(parents=True, exist_ok=True)
    encoding = {
        variable: {
            "dtype": "float64",
            "_FillValue": -9999.0,
            "missing_value": -9999.0,
            "zlib": True,
            "complevel": 4,
        },
        "time": {
            "dtype": "int32",
            "units": "days since 1900-01-01 00:00:00",
            "calendar": "proleptic_gregorian",
        },
        "time_bnds": {
            "dtype": "int32",
            "units": "days since 1900-01-01 00:00:00",
            "calendar": "proleptic_gregorian",
        },
        "lat": {
            "dtype": "float64",
            "_FillValue": None,
        },
        "lon": {
            "dtype": "float64",
            "_FillValue": None,
        },
    }

    last_error = None
    for engine in ("netcdf4", "h5netcdf", "scipy", None):
        try:
            engine_encoding = encoding
            if engine == "scipy":
                engine_encoding = _strip_compression_encoding(encoding)
            kwargs = {} if engine is None else {"engine": engine}
            ds_to_write = _strip_attrs_managed_by_encoding(
                ds, engine_encoding)
            ds_to_write.to_netcdf(
                output_file, encoding=engine_encoding, **kwargs)
            return
        except Exception as e:
            last_error = e

    raise RuntimeError(f"Could not write NetCDF file {output_file}: {last_error}")


def write_header(ds, variable: str, header_file: Path, header: dict | None = None) -> None:
    """Write the mHM ASCII grid header that accompanies a forcing NetCDF."""
    header_file.parent.mkdir(parents=True, exist_ok=True)
    if header is not None:
        header_text = (
            f"ncols         {int(header['ncols'])}\n"
            f"nrows         {int(header['nrows'])}\n"
            f"xllcorner     {header['xllcorner']}\n"
            f"yllcorner     {header['yllcorner']}\n"
            f"cellsize      {header['cellsize']}\n"
            f"NODATA_value  {header.get('nodata_value', -9999.0)}\n"
        )
        header_file.write_text(header_text, encoding="utf-8")
        return

    lon = ds["lon"].values
    lat = ds["lat"].values

    if len(lon) < 2 or len(lat) < 2:
        raise ValueError("Cannot write mHM header for a one-cell meteo grid.")

    dx = abs(float(lon[1] - lon[0]))
    dy = abs(float(lat[1] - lat[0]))
    cellsize = round((dx + dy) / 2.0, 12)
    xllcorner = round(float(min(lon)) - cellsize / 2.0, 12)
    yllcorner = round(float(min(lat)) - cellsize / 2.0, 12)

    header = (
        f"ncols         {len(lon)}\n"
        f"nrows         {len(lat)}\n"
        f"xllcorner     {xllcorner}\n"
        f"yllcorner     {yllcorner}\n"
        f"cellsize      {cellsize}\n"
        "NODATA_value  -9999.0\n"
    )
    header_file.write_text(header, encoding="utf-8")


def _strip_compression_encoding(encoding: dict) -> dict:
    stripped = {}
    for variable, options in encoding.items():
        stripped[variable] = {
            key: value
            for key, value in options.items()
            if key not in {"zlib", "complevel", "shuffle", "chunksizes"}
        }
    return stripped


def _strip_attrs_managed_by_encoding(ds, encoding: dict):
    """Remove attrs that xarray expects to own through encoding."""
    ds_to_write = ds.copy(deep=False)
    for variable, options in encoding.items():
        if variable not in ds_to_write.variables:
            continue
        for key in options:
            ds_to_write[variable].attrs.pop(key, None)
    return ds_to_write
