"""ERA5-Land file naming and inventory helpers."""
from __future__ import annotations

from pathlib import Path


def inspect_era5_folder(nc_folder) -> dict:
    """Inspect ERA5-Land monthly files without opening the NetCDF payloads."""
    counts: dict[str, int] = {}
    first: dict[str, str] = {}
    last: dict[str, str] = {}
    periods: dict[str, list[str]] = {}

    for path in sorted(Path(nc_folder).glob("ERA5_Land_*.nc")):
        parsed = parse_era5_filename(path)
        if not parsed:
            continue
        variable, period = parsed
        counts[variable] = counts.get(variable, 0) + 1
        periods.setdefault(variable, []).append(period)
        first.setdefault(variable, path.name)
        last[variable] = path.name

    return {
        "counts": dict(sorted(counts.items())),
        "first": first,
        "last": last,
        "periods": {key: sorted(value) for key, value in periods.items()},
    }


def parse_era5_filename(path: Path) -> tuple[str, str] | None:
    """Return ``(<variable>, <year_month>)`` for a standard ERA5-Land file."""
    parts = path.stem.split("_")
    if len(parts) < 5 or parts[0] != "ERA5" or parts[1] != "Land":
        return None

    year = parts[-2]
    month = parts[-1]
    if not (year.isdigit() and month.isdigit()):
        return None

    variable = "_".join(parts[2:-2])
    return variable, f"{year}_{month}"


def files_for_variable(folder: Path, variable: str) -> list[Path]:
    """Return all standard ERA5-Land monthly files for a source variable."""
    return sorted(folder.glob(f"ERA5_Land_{variable}_*.nc"))
