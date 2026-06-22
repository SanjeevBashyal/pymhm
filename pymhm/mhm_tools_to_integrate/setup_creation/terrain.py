# -*- coding: utf-8 -*-
"""Backend-neutral morphology terrain product orchestration."""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol


LogCallback = Callable[[str], None]


class TerrainBackend(Protocol):
    """Minimal processing backend required by terrain operations."""

    def run(self, algorithm: str, params: dict):
        """Run an algorithm and return backend-specific result."""


@dataclass(frozen=True)
class TerrainProducts:
    """Paths for region-level terrain derivative outputs."""

    slope: Path
    aspect: Path
    flow_direction: Path | None = None
    flow_accumulation: Path | None = None
    flow_accumulation_area: Path | None = None

    def as_dict(self) -> dict[str, Path]:
        """Return available products keyed by semantic name."""
        return {
            key: value
            for key, value in {
                "slope": self.slope,
                "aspect": self.aspect,
                "flow_direction": self.flow_direction,
                "flow_accumulation": self.flow_accumulation,
                "flow_accumulation_area": self.flow_accumulation_area,
            }.items()
            if value is not None
        }


def aspect_params(dem_path: str | Path, output_path: str | Path) -> dict:
    """Return GDAL aspect algorithm parameters."""
    return {
        "INPUT": str(dem_path),
        "BAND": 1,
        "TRIG_ANGLE": False,
        "ZERO_FLAT": False,
        "COMPUTE_EDGES": False,
        "ZEVENBERGEN": False,
        "OPTIONS": None,
        "EXTRA": "",
        "OUTPUT": str(output_path),
    }


def slope_params(
        dem_path: str | Path,
        output_path: str | Path,
        scale: float = 1.0) -> dict:
    """Return GDAL slope algorithm parameters."""
    return {
        "INPUT": str(dem_path),
        "BAND": 1,
        "SCALE": scale,
        "AS_PERCENT": True,
        "COMPUTE_EDGES": False,
        "ZEVENBERGEN": False,
        "OPTIONS": None,
        "EXTRA": "",
        "OUTPUT": str(output_path),
    }


def create_terrain_products(
        dem_path: str | Path,
        output_folder: str | Path,
        backend: TerrainBackend,
        slope_scale: float = 1.0,
        log: LogCallback | None = None) -> TerrainProducts:
    """Create slope and aspect products using a caller-provided backend."""
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    products = TerrainProducts(
        slope=output_folder / "1_dem_slope.tif",
        aspect=output_folder / "1_dem_aspect.tif",
    )

    if log:
        log("Creating region aspect.")
    backend.run("gdal:aspect", aspect_params(dem_path, products.aspect))

    if log:
        log("Creating region slope.")
    backend.run("gdal:slope", slope_params(dem_path, products.slope, slope_scale))

    return products
