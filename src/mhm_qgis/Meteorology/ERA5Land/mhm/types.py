"""Shared data structures for mHM forcing preparation."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ForcingSpec:
    """Definition of one mHM forcing output generated from ERA5-Land data."""

    output_variable: str
    output_folder: str
    source_key: str
    file_variable: str
    aggregation: str
    units: str
    long_name: str
    scale: float = 1.0
    offset: float = 0.0


@dataclass
class MeteoForcingResult:
    """Final files prepared for mHM."""

    outputs: dict[str, Path]
    headers: dict[str, Path]
    files_used: dict[str, int]
