"""QGIS-free LAI temporal aggregation and repetition."""
from __future__ import annotations

from dataclasses import dataclass


TARGET_TIMESTEPS = {
    "daily gridded data": (-1, "D"),
    "monthly gridded data": (-2, "MS"),
    "yearly gridded data": (-3, "YS"),
    "long term mean monthly gridded data": (1, None),
}
_INPUT_DAYS = {
    "daily": 1.0,
    "semi-weekly": 3.5,
    "weekly": 7.0,
    "biweekly": 14.0,
    "monthly": 30.4375,
    "annual": 365.25,
}
_TARGET_DAYS = {-1: 1.0, -2: 30.4375, -3: 365.25}


@dataclass(frozen=True)
class LaiTemporalResult:
    data: object
    time_bounds: object
    time_step: int
    target_timestep: str


def lai_time_step(target_timestep: str) -> int:
    normalized = str(target_timestep or "").strip().lower()
    try:
        return TARGET_TIMESTEPS[normalized][0]
    except KeyError as error:
        raise ValueError(f"Unsupported LAI target timestep: {target_timestep}") from error


def prepare_lai_temporal(data, input_resolution, target_timestep, time_dim="time"):
    """Convert one LAI data array to the requested mHM time convention."""
    import numpy as np
    import pandas as pd

    input_name = str(input_resolution or "").strip().lower()
    if input_name not in _INPUT_DAYS:
        raise ValueError(f"Unsupported LAI input resolution: {input_resolution}")
    target_name = str(target_timestep or "").strip().lower()
    if target_name not in TARGET_TIMESTEPS:
        raise ValueError(f"Unsupported LAI target timestep: {target_timestep}")
    if time_dim not in data.dims:
        raise ValueError("The selected LAI variable has no time dimension.")
    if data.sizes.get(time_dim, 0) == 0:
        raise ValueError("The selected LAI variable has an empty time dimension.")

    values = data[time_dim].values
    if not np.issubdtype(np.asarray(values).dtype, np.datetime64):
        if target_name == "long term mean monthly gridded data" and data.sizes[time_dim] == 12:
            monthly = data.rename({time_dim: "time"}).assign_coords(
                time=np.arange(1, 13, dtype=np.int32)
            )
            bounds = np.column_stack(
                (np.arange(0, 12, dtype=float), np.arange(1, 13, dtype=float))
            )
            return LaiTemporalResult(monthly, bounds, 1, target_timestep)
        raise ValueError("LAI time coordinates must be decoded calendar dates.")

    data = data.sortby(time_dim)
    if target_name == "long term mean monthly gridded data":
        monthly = data.groupby(f"{time_dim}.month").mean(time_dim, skipna=True)
        if monthly.sizes.get("month") != 12:
            raise ValueError("Long-term monthly LAI requires observations for all 12 months.")
        monthly = monthly.rename({"month": "time"}).assign_coords(
            time=np.arange(1, 13, dtype=np.int32)
        )
        bounds = np.column_stack(
            (np.arange(0, 12, dtype=float), np.arange(1, 13, dtype=float))
        )
        return LaiTemporalResult(monthly, bounds, 1, target_timestep)

    target_step, frequency = TARGET_TIMESTEPS[target_name]
    source_days = _INPUT_DAYS[input_name]
    target_days = _TARGET_DAYS[target_step]
    if source_days < target_days:
        converted = data.resample({time_dim: frequency}).mean(skipna=True)
    elif source_days > target_days:
        times = pd.DatetimeIndex(data[time_dim].values)
        end = times[-1] + _source_offset(input_name)
        target_index = pd.date_range(
            start=times[0], end=end, freq=frequency, inclusive="left"
        )
        converted = data.reindex({time_dim: target_index}, method="ffill")
    else:
        converted = data

    if time_dim != "time":
        converted = converted.rename({time_dim: "time"})
    converted = converted.astype("float64")
    starts = pd.DatetimeIndex(converted["time"].values)
    ends = starts + _target_offset(target_step)
    bounds = np.column_stack(
        (starts.values.astype("datetime64[ns]"), ends.values.astype("datetime64[ns]"))
    )
    return LaiTemporalResult(converted, bounds, target_step, target_timestep)


def _source_offset(name):
    import pandas as pd

    return {
        "daily": pd.offsets.Day(1),
        "semi-weekly": pd.Timedelta(days=3, hours=12),
        "weekly": pd.offsets.Week(1),
        "biweekly": pd.offsets.Week(2),
        "monthly": pd.offsets.MonthBegin(1),
        "annual": pd.offsets.YearBegin(1),
    }[name]


def _target_offset(step):
    import pandas as pd

    return {
        -1: pd.offsets.Day(1),
        -2: pd.offsets.MonthBegin(1),
        -3: pd.offsets.YearBegin(1),
    }[step]


__all__ = [
    "LaiTemporalResult",
    "TARGET_TIMESTEPS",
    "lai_time_step",
    "prepare_lai_temporal",
]
