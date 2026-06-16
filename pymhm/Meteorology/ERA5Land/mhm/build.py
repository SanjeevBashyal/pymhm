"""Build daily mHM forcing datasets from monthly ERA5-Land files."""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable

from ..config import VARIABLE_NAME_MAPPINGS
from .aggregation import aggregate_daily, drop_duplicate_times
from .dataset import find_source_variable, open_dataset, standardize_dataset
from .dependencies import import_dependencies
from .grid import (
    ensure_latitude_descending,
    force_target_grid,
    normalize_spatial_axes,
    resample_to_target_grid,
    snap_bounds_to_grid,
    subset_dataset,
)
from .logging import log_message
from .metadata import add_time_bounds, finalize_geospatial_metadata
from .types import ForcingSpec


def build_daily_dataset(
    files: Iterable[Path],
    spec: ForcingSpec,
    bounds: tuple[float, float, float, float] | None,
    target_lat=None,
    target_lon=None,
    log: Callable[[str], None] | None = None,
):
    """Build one complete daily forcing dataset from monthly ERA5-Land files."""
    np, pd, xr = import_dependencies()
    explicit_target_lat = None if target_lat is None else np.asarray(target_lat, dtype="float64")
    explicit_target_lon = None if target_lon is None else np.asarray(target_lon, dtype="float64")
    has_explicit_target = explicit_target_lat is not None and explicit_target_lon is not None
    daily_arrays = []
    target_lat = None
    target_lon = None
    snapped_bounds = None
    aliases = VARIABLE_NAME_MAPPINGS.get(spec.source_key, [])

    for path in files:
        with open_dataset(path, xr) as ds_raw:
            ds = standardize_dataset(ds_raw)
            ds = normalize_spatial_axes(ds, np)
            source_variable = find_source_variable(ds, aliases)
            if source_variable is None:
                raise KeyError(
                    f"No {spec.source_key} variable found in {path.name}. "
                    f"Available variables: {list(ds.data_vars)}")

            ds = ds[[source_variable]]
            if snapped_bounds is None and bounds is not None:
                snapped_bounds = snap_bounds_to_grid(ds, bounds, np)
                log_message(
                    log,
                    "ERA5-Land crop bounds snapped to grid centers: "
                    f"west={snapped_bounds[0]:.6f}, east={snapped_bounds[1]:.6f}, "
                    f"south={snapped_bounds[2]:.6f}, north={snapped_bounds[3]:.6f}",
                )

            if snapped_bounds is not None:
                ds = subset_dataset(ds, snapped_bounds)

            ds = normalize_spatial_axes(ds, np)
            ds = ensure_latitude_descending(ds)
            daily = aggregate_daily(ds[source_variable], spec, np)
            daily = daily.astype("float64").load()
            daily = normalize_spatial_axes(daily, np)
            daily = ensure_latitude_descending(daily)
            daily = daily.transpose("time", "lat", "lon")

            if has_explicit_target:
                daily = resample_to_target_grid(
                    daily,
                    explicit_target_lat,
                    explicit_target_lon,
                    np,
                )
                if explicit_target_lat[0] < explicit_target_lat[-1]:
                    daily = ensure_latitude_descending(daily)
                target_lat = daily["lat"].values
                target_lon = daily["lon"].values
            elif target_lat is None or target_lon is None:
                target_lat = daily["lat"].values
                target_lon = daily["lon"].values
            else:
                daily = force_target_grid(daily, target_lat, target_lon, np)

            daily.name = spec.output_variable
            daily.attrs.clear()
            daily.attrs.update({
                "units": spec.units,
                "long_name": spec.long_name,
            })
            daily_arrays.append(daily)

    if not daily_arrays:
        raise RuntimeError(f"No daily data produced for {spec.output_variable}.")

    combined = xr.concat(
        daily_arrays,
        dim="time",
        coords="minimal",
        compat="override",
        join="override",
    )
    combined = combined.assign_coords(lat=target_lat, lon=target_lon)
    combined = combined.sortby("time")
    combined = normalize_spatial_axes(combined, np)
    combined = ensure_latitude_descending(combined)
    combined = drop_duplicate_times(combined, np)
    ds_out = combined.to_dataset(name=spec.output_variable)
    ds_out = add_time_bounds(ds_out, np, pd, xr)
    ds_out = finalize_geospatial_metadata(ds_out, spec.output_variable)
    ds_out.attrs.update({
        "Conventions": "CF-1.6",
        "source": "ERA5-Land",
        "history": "Prepared for mHM by pymhm",
    })
    return ds_out
