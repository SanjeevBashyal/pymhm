"""Spatial axis normalization and ERA5-Land grid helpers."""
from __future__ import annotations


def snap_bounds_to_grid(ds, bounds, np) -> tuple[float, float, float, float]:
    """Snap WGS84 bounds outward to available ERA5-Land grid centers."""
    west, east, south, north = bounds
    lon_values = np.sort(np.asarray(ds["lon"].values, dtype=float))
    lat_values = np.sort(np.asarray(ds["lat"].values, dtype=float))

    return (
        _coord_floor(lon_values, min(west, east)),
        _coord_ceil(lon_values, max(west, east)),
        _coord_floor(lat_values, min(south, north)),
        _coord_ceil(lat_values, max(south, north)),
    )


def subset_dataset(ds, bounds):
    """Subset a dataset to snapped WGS84 bounds without changing axis order."""
    west, east, south, north = bounds
    lon_ascending = float(ds["lon"].values[0]) <= float(ds["lon"].values[-1])
    lat_ascending = float(ds["lat"].values[0]) <= float(ds["lat"].values[-1])

    lon_slice = slice(west, east) if lon_ascending else slice(east, west)
    lat_slice = slice(south, north) if lat_ascending else slice(north, south)
    subset = ds.sel(lon=lon_slice, lat=lat_slice)

    if subset.sizes.get("lat", 0) == 0 or subset.sizes.get("lon", 0) == 0:
        raise ValueError(
            "ERA5-Land crop produced an empty grid. "
            f"Bounds were west={west}, east={east}, south={south}, north={north}.")
    return subset


def ensure_latitude_descending(ds):
    """Keep mHM output orientation north-to-south."""
    if ds.sizes.get("lat", 0) > 1 and float(ds["lat"][0]) < float(ds["lat"][-1]):
        return ds.sortby("lat", ascending=False)
    return ds


def normalize_spatial_axes(obj, np, decimals: int = 6):
    """Round and deduplicate 1D lat/lon axes so output axes stay unique."""
    for axis in ("lat", "lon"):
        if axis not in obj.coords:
            continue
        values = np.asarray(obj[axis].values, dtype="float64")
        if values.ndim != 1:
            raise ValueError(f"Coordinate {axis} must be one-dimensional.")

        rounded = np.round(values, decimals)
        obj = obj.assign_coords({axis: rounded})
        _, first_indices = np.unique(rounded, return_index=True)
        if len(first_indices) != len(rounded):
            obj = obj.isel({axis: np.sort(first_indices)})
    return obj


def force_target_grid(da, target_lat, target_lon, np):
    """Force a monthly daily array onto the first month's canonical grid."""
    lat_values = np.asarray(da["lat"].values, dtype="float64")
    lon_values = np.asarray(da["lon"].values, dtype="float64")
    if np.array_equal(lat_values, target_lat) and np.array_equal(lon_values, target_lon):
        return da

    lat_tolerance = _axis_tolerance(target_lat, np)
    lon_tolerance = _axis_tolerance(target_lon, np)
    return da.reindex(
        lat=target_lat,
        lon=target_lon,
        method="nearest",
        tolerance=max(lat_tolerance, lon_tolerance),
    )


def resample_to_target_grid(da, target_lat, target_lon, np):
    """Nearest-neighbor resample a data array to an explicit target grid."""
    if target_lat is None or target_lon is None:
        return da

    target_lat = np.asarray(target_lat, dtype="float64")
    target_lon = np.asarray(target_lon, dtype="float64")
    source = da.sortby("lat").sortby("lon")
    resampled = source.reindex(
        lat=target_lat,
        lon=target_lon,
        method="nearest",
    )
    return resampled.transpose("time", "lat", "lon")


def resample_to_target_points(
        da,
        sample_lat,
        sample_lon,
        target_lat,
        target_lon,
        np):
    """Nearest-neighbour sample at exact 2-D WGS84 cell centres."""
    import xarray as xr
    from scipy.spatial import cKDTree

    source = da.sortby("lat").sortby("lon").transpose("time", "lat", "lon")
    source_lon, source_lat = np.meshgrid(
        np.asarray(source["lon"].values, dtype="float64"),
        np.asarray(source["lat"].values, dtype="float64"),
    )
    source_points = np.column_stack(
        [source_lon.reshape(-1), source_lat.reshape(-1)]
    )
    target_points = np.column_stack([
        np.asarray(sample_lon, dtype="float64").reshape(-1),
        np.asarray(sample_lat, dtype="float64").reshape(-1),
    ])
    _, indices = cKDTree(source_points).query(target_points)
    values = np.asarray(source.values)
    sampled = values.reshape(values.shape[0], -1)[:, indices]
    sampled = sampled.reshape(
        values.shape[0],
        len(target_lat),
        len(target_lon),
    )
    return xr.DataArray(
        sampled,
        dims=("time", "lat", "lon"),
        coords={
            "time": source["time"],
            "lat": np.asarray(target_lat, dtype="float64"),
            "lon": np.asarray(target_lon, dtype="float64"),
        },
        name=da.name,
        attrs=dict(da.attrs),
    )


def _coord_floor(values, target: float) -> float:
    candidates = values[values <= target]
    if len(candidates):
        return float(candidates[-1])
    return float(values[0])


def _coord_ceil(values, target: float) -> float:
    candidates = values[values >= target]
    if len(candidates):
        return float(candidates[0])
    return float(values[-1])


def _axis_tolerance(values, np) -> float:
    values = np.asarray(values, dtype="float64")
    if len(values) < 2:
        return 1e-8
    diffs = np.diff(np.sort(values))
    diffs = diffs[diffs > 0]
    if len(diffs) == 0:
        return 1e-8
    return float(np.min(diffs) / 4.0)
