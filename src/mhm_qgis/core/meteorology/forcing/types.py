# -*- coding: utf-8 -*-
"""Value objects shared by meteorology inspection and preparation."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence


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
class TargetGrid:
    """WGS84 cell-centre axes and the matching projected mHM header."""

    lon: Sequence[float]
    lat: Sequence[float]
    header: Mapping[str, float]
    crs: str | None = None
    sample_lon: object | None = None
    sample_lat: object | None = None

    def validate(self) -> None:
        if len(self.lon) == 0 or len(self.lat) == 0:
            raise ValueError("The target meteorology grid is empty.")
        if int(self.header["ncols"]) != len(self.lon):
            raise ValueError("Target longitude count does not match header ncols.")
        if int(self.header["nrows"]) != len(self.lat):
            raise ValueError("Target latitude count does not match header nrows.")
        if (self.sample_lon is None) != (self.sample_lat is None):
            raise ValueError(
                "Both target sample longitude and latitude grids are required."
            )
        if self.sample_lon is not None:
            import numpy as np

            shape = (int(self.header["nrows"]), int(self.header["ncols"]))
            if (
                    np.shape(self.sample_lon) != shape
                    or np.shape(self.sample_lat) != shape):
                raise ValueError(
                    "Target sample coordinates do not match the header shape."
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


@dataclass(frozen=True)
class SpatialResolution:
    """A grid resolution expressed in one CRS."""

    resolution: float
    x_resolution: float
    y_resolution: float
    unit: str
    crs: str | None


__all__ = [
    "ERA5LAND",
    "MHM_READY",
    "MeteoFolderSpec",
    "SpatialMetadata",
    "SpatialResolution",
    "TargetGrid",
    "normalize_kind",
    "normalize_source",
]
