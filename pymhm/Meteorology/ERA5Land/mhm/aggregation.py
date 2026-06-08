"""Daily aggregation rules for mHM forcing variables."""
from __future__ import annotations

from .types import ForcingSpec


def aggregate_daily(da, spec: ForcingSpec, np):
    """Aggregate hourly ERA5-Land values to daily mHM forcing values."""
    if spec.aggregation == "sum":
        daily = da.resample(time="1D").sum(skipna=True)
    elif spec.aggregation == "era5land_daily_total":
        daily = aggregate_era5land_daily_total(da, np)
    elif spec.aggregation == "mean":
        daily = da.resample(time="1D").mean(skipna=True)
    elif spec.aggregation == "min":
        daily = da.resample(time="1D").min(skipna=True)
    elif spec.aggregation == "max":
        daily = da.resample(time="1D").max(skipna=True)
    else:
        raise ValueError(f"Unsupported aggregation: {spec.aggregation}")

    return daily * spec.scale + spec.offset


def aggregate_era5land_daily_total(da, np):
    """Derive daily precipitation totals from ERA5-Land accumulated tp."""
    fallback = da.resample(time="1D").max(skipna=True)
    try:
        zero_utc = da.where(da["time"].dt.hour == 0, drop=True)
    except Exception:
        return fallback

    if zero_utc.sizes.get("time", 0) == 0:
        return fallback

    shifted = zero_utc.assign_coords(
        time=zero_utc["time"] - np.timedelta64(1, "D"))
    shifted = shifted.reindex(time=fallback["time"])
    return shifted.combine_first(fallback)


def drop_duplicate_times(da, np):
    """Keep the first value for duplicated daily timestamps."""
    values = da["time"].values
    _, unique_indices = np.unique(values, return_index=True)
    if len(unique_indices) == len(values):
        return da
    return da.isel(time=np.sort(unique_indices))
