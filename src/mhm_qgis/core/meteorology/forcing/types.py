# -*- coding: utf-8 -*-
"""Value objects shared by meteorology inspection and preparation."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


MHM_READY = "mHM ready"
ERA5LAND = "ERA5Land"
_KINDS = {
    "pre": "precipitation",
    "precipitation": "precipitation",
    "temp": "temperature",
    "temperature": "temperature",
    "pet": "pet",
}


def normalize_kind(kind: str) -> str:
    """Return the canonical meteorology input kind."""
    normalized = _KINDS.get(str(kind).strip().lower())
    if not normalized:
        raise ValueError(f"Unsupported meteorology data type: {kind}")
    return normalized


def normalize_source(source: str) -> str:
    """Normalize UI and legacy spellings of a meteorology data source."""
    compact = str(source).strip().lower().replace("_", "").replace("-", "")
    compact = compact.replace(" ", "")
    if compact == "mhmready":
        return MHM_READY
    if compact == "era5land":
        return ERA5LAND
    raise ValueError(f"Unsupported meteorology data source: {source}")


@dataclass(frozen=True)
class MeteoFolderSpec:
    """One user-selected meteorology input folder."""

    kind: str
    folder: Path | str
    source: str
    crs: str | None = None

    def normalized(self) -> "MeteoFolderSpec":
        return MeteoFolderSpec(
            kind=normalize_kind(self.kind),
            folder=Path(self.folder),
            source=normalize_source(self.source),
            crs=self.crs,
        )


@dataclass(frozen=True)
class SpatialMetadata:
    """Spatial metadata shared by all NetCDF files in an input folder."""

    kind: str
    source: str
    files: tuple[Path, ...]
    shape: tuple[int, int]
    x_coordinate: str
    y_coordinate: str
    x_resolution: float
    y_resolution: float
    resolution: float
    bounds: tuple[float, float, float, float]
    crs: str | None
    unit: str


__all__ = [
    "ERA5LAND",
    "MHM_READY",
    "MeteoFolderSpec",
    "SpatialMetadata",
    "normalize_kind",
    "normalize_source",
]
