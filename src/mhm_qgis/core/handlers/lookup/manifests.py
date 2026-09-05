"""Validated manifests for time-varying land cover and soil horizons."""
from __future__ import annotations

import csv
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


LAND_USE_HEADER = ("StartYear", "EndYear", "FilePath")
SOIL_HEADER = (
    "Horizon",
    "Upper Depth",
    "Lower Depth",
    "Clay Layer",
    "Sand Layer",
    "Silt Layer",
    "Bulk Density Layer",
    "Bulk Density Unit",
)
BULK_DENSITY_UNITS = ("g/cm3", "kg/m3", "cg/cm3")
MAX_LAND_USE_PERIODS = 50
MAX_SOIL_HORIZONS = 10


@dataclass(frozen=True)
class LandUsePeriod:
    start_year: int
    end_year: int
    file_path: Path

    def as_dict(self) -> dict:
        return {
            "start_year": self.start_year,
            "end_year": self.end_year,
            "file_path": str(self.file_path),
        }


@dataclass(frozen=True)
class LandUseInput:
    periods: tuple[LandUsePeriod, ...]
    lookup_table: Path
    mapping_field: str
    class_field: str

    def as_dict(self) -> dict:
        return {
            "periods": [period.as_dict() for period in self.periods],
            "lookup_table": str(self.lookup_table),
            "mapping_field": self.mapping_field,
            "class_field": self.class_field,
        }


@dataclass(frozen=True)
class SoilHorizon:
    horizon: int
    upper_depth: float
    lower_depth: float
    clay_layer: Path
    sand_layer: Path
    silt_layer: Path
    bulk_density_layer: Path

    def as_dict(self) -> dict:
        return {
            "horizon": self.horizon,
            "upper_depth": self.upper_depth,
            "lower_depth": self.lower_depth,
            "clay_layer": str(self.clay_layer),
            "sand_layer": str(self.sand_layer),
            "silt_layer": str(self.silt_layer),
            "bulk_density_layer": str(self.bulk_density_layer),
        }


@dataclass(frozen=True)
class SoilInput:
    horizons: tuple[SoilHorizon, ...]
    bulk_density_unit: str

    @property
    def lower_depths(self) -> tuple[float, ...]:
        return tuple(horizon.lower_depth for horizon in self.horizons)

    def as_dict(self) -> dict:
        return {
            "bulk_density_unit": self.bulk_density_unit,
            "lower_depths": list(self.lower_depths),
            "horizons": [horizon.as_dict() for horizon in self.horizons],
        }


def validate_land_use_input(value: LandUseInput) -> None:
    """Validate paths and continuous, non-overlapping calendar periods."""
    if not value.periods:
        raise ValueError("Add at least one land-use layer.")
    if len(value.periods) > MAX_LAND_USE_PERIODS:
        raise ValueError(f"Land use supports at most {MAX_LAND_USE_PERIODS} periods.")
    previous_end = None
    for number, period in enumerate(value.periods, 1):
        if period.end_year < period.start_year:
            raise ValueError(f"Land-use layer {number} ends before it starts.")
        if previous_end is not None and period.start_year != previous_end + 1:
            raise ValueError("Land-use periods must be ordered and gap-free.")
        _require_file(period.file_path, f"land-use layer {number}")
        previous_end = period.end_year
    _require_file(value.lookup_table, "land-use lookup table")
    if not value.mapping_field.strip() or not value.class_field.strip():
        raise ValueError("Select both the land-use mapping and class fields.")


def validate_soil_input(value: SoilInput) -> None:
    """Validate contiguous soil horizons and their four source layers."""
    if value.bulk_density_unit not in BULK_DENSITY_UNITS:
        raise ValueError("Select the bulk-density unit.")
    if not value.horizons:
        raise ValueError("Add at least one soil horizon.")
    if len(value.horizons) > MAX_SOIL_HORIZONS:
        raise ValueError(f"Soil supports at most {MAX_SOIL_HORIZONS} horizons.")
    previous_lower = None
    for number, horizon in enumerate(value.horizons, 1):
        if horizon.horizon != number:
            raise ValueError("Soil horizons must be numbered consecutively.")
        if horizon.upper_depth < 0 or horizon.lower_depth <= horizon.upper_depth:
            raise ValueError(f"Soil horizon {number} has invalid depth bounds.")
        if number == 1 and abs(horizon.upper_depth) > 1e-9:
            raise ValueError("The first soil horizon must start at depth 0.")
        if (
            previous_lower is not None
            and abs(horizon.upper_depth - previous_lower) > 1e-9
        ):
            raise ValueError("Soil horizons must be contiguous.")
        for label, path in (
            ("clay", horizon.clay_layer),
            ("sand", horizon.sand_layer),
            ("silt", horizon.silt_layer),
            ("bulk-density", horizon.bulk_density_layer),
        ):
            _require_file(path, f"{label} layer for soil horizon {number}")
        previous_lower = horizon.lower_depth


def write_land_use_manifest(value: LandUseInput, output_file) -> Path:
    """Write ``format-data.csv`` for historical land-use processing."""
    validate_land_use_input(value)
    output = Path(output_file)
    rows = (
        (period.start_year, period.end_year, _manifest_path(period.file_path, output))
        for period in value.periods
    )
    return _write_csv(output, LAND_USE_HEADER, rows)


def write_soil_manifest(value: SoilInput, output_file) -> Path:
    """Write soil ``format-data.csv`` for multi-horizon processing.

    mHM-tools reads the manifest with a plain ``read_csv``, so the bulk-density
    unit is repeated as a column on every row rather than carried in a preamble
    line: a preamble would be consumed as the header.
    """
    validate_soil_input(value)
    output = Path(output_file)
    rows = (
        (
            horizon.horizon,
            _number(horizon.upper_depth),
            _number(horizon.lower_depth),
            _manifest_path(horizon.clay_layer, output),
            _manifest_path(horizon.sand_layer, output),
            _manifest_path(horizon.silt_layer, output),
            _manifest_path(horizon.bulk_density_layer, output),
            value.bulk_density_unit,
        )
        for horizon in value.horizons
    )
    return _write_csv(output, SOIL_HEADER, rows)


def _require_file(path: Path, label: str) -> None:
    if not Path(path).is_file():
        raise ValueError(f"Select a readable {label}: {path}")


def _manifest_path(path: Path, manifest: Path) -> str:
    source = Path(path).expanduser().resolve()
    try:
        return source.relative_to(manifest.parent.resolve()).as_posix()
    except ValueError:
        return str(source)


def _number(value: float):
    return int(value) if float(value).is_integer() else value


def _write_csv(
    output: Path,
    header: Iterable[str],
    rows: Iterable[Iterable],
) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="",
            prefix=f".{output.name}.",
            suffix=".tmp",
            dir=output.parent,
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, output)
    except Exception:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise
    return output


__all__ = [
    "BULK_DENSITY_UNITS",
    "MAX_LAND_USE_PERIODS",
    "MAX_SOIL_HORIZONS",
    "LandUseInput",
    "LandUsePeriod",
    "SoilHorizon",
    "SoilInput",
    "validate_land_use_input",
    "validate_soil_input",
    "write_land_use_manifest",
    "write_soil_manifest",
]
