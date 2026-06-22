"""CF metadata helpers for mHM NetCDF output."""
from __future__ import annotations


def add_time_bounds(ds, np, pd, xr):
    """Add daily time bounds required by mHM-style forcing files."""
    time_index = pd.DatetimeIndex(ds["time"].values)
    starts = time_index.values.astype("datetime64[ns]")
    ends = (time_index + pd.Timedelta(days=1)).values.astype("datetime64[ns]")
    bounds = xr.DataArray(
        data=np.stack([starts, ends], axis=1),
        dims=("time", "bnds"),
        coords={"time": ds["time"], "bnds": [0, 1]},
        name="time_bnds",
    )
    ds["time_bnds"] = bounds
    ds["time"].attrs.update({
        "standard_name": "time",
        "axis": "T",
        "bounds": "time_bnds",
    })
    return ds


def finalize_geospatial_metadata(ds, variable: str):
    """Make output coordinates self-describing for CF readers and Panoply."""
    ds = ds.drop_vars(
        [name for name in ("number", "expver", "spatial_ref") if name in ds],
        errors="ignore",
    )

    ds["lat"].attrs.clear()
    ds["lat"].attrs.update({
        "standard_name": "latitude",
        "long_name": "latitude",
        "units": "degrees_north",
        "axis": "Y",
    })
    ds["lon"].attrs.clear()
    ds["lon"].attrs.update({
        "standard_name": "longitude",
        "long_name": "longitude",
        "units": "degrees_east",
        "axis": "X",
    })
    ds["time"].attrs.setdefault("long_name", "time")
    ds[variable].attrs.pop("coordinates", None)
    ds[variable].attrs.pop("grid_mapping", None)
    return ds
