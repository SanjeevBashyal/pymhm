"""Resolve prepared morphology outputs for the main display selector."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .nml_settings import load_settings
from .project_layout import geometry_folder, morph_folder, workspace_folder


DISPLAY_KEYS = (
    "dem",
    "slope",
    "aspect",
    "flow_accumulation",
    "flow_direction",
    "idgauges",
    "land_cover",
    "soil",
    "geology",
    "lai",
    "domain_mask",
)


@dataclass(frozen=True)
class DisplayOutput:
    path: Path
    name: str
    is_raster: bool = True
    band: int | None = None
    variable: str | None = None


def land_cover_periods(project_folder: str | Path) -> list[dict]:
    land_cover = load_settings(project_folder).get("land_cover", {})
    scenes = land_cover.get("scenes", []) if isinstance(land_cover, dict) else []
    return [scene for scene in scenes if isinstance(scene, dict)]


def resolve_display_output(
    project_folder: str | Path,
    key: str,
    *,
    year: int | None = None,
) -> DisplayOutput | None:
    project_folder = str(project_folder)
    geometry = Path(geometry_folder(project_folder))
    morph = Path(morph_folder(project_folder))
    workspace = Path(workspace_folder(project_folder))
    settings = load_settings(project_folder)

    if key == "land_cover":
        return _land_cover_output(settings, workspace, geometry, morph, year)
    if key == "soil":
        configured = _configured_path(settings.get("soil"), workspace)
        candidates = [configured] if configured else []
        candidates += [
            morph / "soil_horizon_class.nc",
            morph / "soil_class.asc",
            *_raster_variants(geometry / "3_soil.tif"),
        ]
        return _first(candidates, "Soil", variable="soil_class")
    if key == "domain_mask":
        merged = geometry / "Watersheds" / "4_watershed_merged_vector.shp"
        if merged.is_file():
            return DisplayOutput(merged, "Domain Mask", is_raster=False)
        return _first(
            [geometry / "Watersheds" / "4_watershed_DEM.tif"],
            "Domain Mask",
        )

    candidates = {
        "dem": [*_raster_variants(geometry / "1_dem_filled.tif"), morph / "dem.asc"],
        "slope": [*_raster_variants(geometry / "1_dem_slope.tif"), morph / "slope.asc"],
        "aspect": [*_raster_variants(geometry / "1_dem_aspect.tif"), morph / "aspect.asc"],
        "flow_accumulation": [
            *_raster_variants(geometry / "2_flow_accumulation.tif"),
            morph / "facc.asc",
        ],
        "flow_direction": [
            *_raster_variants(geometry / "2_flow_direction.tif"),
            morph / "fdir.asc",
        ],
        "idgauges": [
            *_raster_variants(geometry / "2_gauge_position.tif"),
            morph / "idgauges.asc",
        ],
        "geology": [
            morph / "geology_class.asc",
            *_raster_variants(geometry / "3_geology_processed.tif"),
        ],
        "lai": [
            workspace / "data" / "lai" / "lai_masked.nc",
            workspace / "data" / "lai" / "lai.nc",
        ],
    }.get(key, [])
    return _first(
        candidates,
        key.replace("_", " ").title(),
        variable="lai" if key == "lai" else None,
    )


def _land_cover_output(settings, workspace, geometry, morph, year):
    land_cover = settings.get("land_cover", {})
    if not isinstance(land_cover, dict):
        land_cover = {}
    scenes = [scene for scene in land_cover.get("scenes", []) if isinstance(scene, dict)]
    selected_index = 0
    if scenes and year is not None:
        for index, scene in enumerate(scenes):
            try:
                if int(scene["start_year"]) <= year <= int(scene["end_year"]):
                    selected_index = index
                    break
            except (KeyError, TypeError, ValueError):
                continue
    if scenes:
        configured = _configured_path(scenes[selected_index], workspace)
        if configured and configured.is_file():
            return DisplayOutput(
                configured, "Land Cover", band=1, variable="land_cover"
            )
    configured = _configured_path(land_cover, workspace)
    if configured and configured.is_file():
        band = selected_index + 1 if configured.suffix.lower() == ".nc" else 1
        return DisplayOutput(
            configured, "Land Cover", band=band, variable="land_cover"
        )
    return _first(
        [
            morph / "lc_periods.nc",
            morph / "lc.asc",
            *_raster_variants(geometry / "3_land_use.tif"),
        ],
        "Land Cover",
        band=selected_index + 1,
        variable="land_cover",
    )


def _configured_path(value, workspace):
    if not isinstance(value, dict):
        return None
    raw = value.get("output_path") or value.get("path") or value.get("filename")
    if not raw:
        return None
    path = Path(str(raw))
    return path if path.is_absolute() else workspace / path


def _raster_variants(path: Path) -> tuple[Path, Path, Path]:
    return (
        path.with_name(f"{path.stem}_masked{path.suffix}"),
        path.with_name(f"{path.stem}_crop{path.suffix}"),
        path,
    )


def _first(candidates, name, band=None, variable=None):
    for candidate in candidates:
        if candidate is not None and Path(candidate).is_file():
            return DisplayOutput(
                Path(candidate), name, band=band, variable=variable
            )
    return None


__all__ = [
    "DISPLAY_KEYS",
    "DisplayOutput",
    "land_cover_periods",
    "resolve_display_output",
]
