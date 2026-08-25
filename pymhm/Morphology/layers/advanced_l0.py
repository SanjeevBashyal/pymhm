# -*- coding: utf-8 -*-
"""Publish staged morphology inputs into `data/master` on the common L0 grid.

Execute All writes land cover, soil, and geology outputs into
`Z Temp/Morphology/morph`. Nothing is model-ready until the common L0 extent is
known, which only happens once meteorology has fixed the L2 grid, so publication
is the last step of Morphology Setup: each staged raster is padded onto the L0
header, moved into `data/master/static/morph`, and the paths recorded in
`nml-settings.json` are rewritten to point at the published copies.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any, Mapping

from ...grid_resolution import CATEGORICAL_PAD_VALUE, LAI_PAD_VALUE
from ...nml_settings import load_settings, relative_workspace_path, save_settings
from ...project_layout import morph_folder, morph_staging_folder, workspace_folder
from ..latlon.ascii_morphology import pad_l0_file_to_header


# Rasters are padded onto the L0 grid; everything else is carried across as is.
RASTER_SUFFIXES = {".asc", ".nc"}
# Class definitions and projection sidecars travel with their raster untouched.
CARRIED_SUFFIXES = {".txt", ".prj", ".xml", ".cpg", ".csv"}
PATH_KEYS = ("output_path", "classdefinition_path")
SECTIONS = ("land_cover", "soil", "geology")
# Staged rasters are named after their model input, so the file name is what
# tells us how to pad them: (name prefixes, pad value, integer output).
PAD_SPECS = (
    (
        ("lc", "land_cover", "land_use", "landuse", "soil", "geology"),
        CATEGORICAL_PAD_VALUE,
        True,
    ),
    (("lai",), LAI_PAD_VALUE, False),
)


def pad_spec_for(path: Path) -> tuple[float | int | None, bool]:
    """Return the pad value and integer flag for one staged raster.

    Class layers pad with a valid class and LAI with a valid leaf area index;
    padding them with nodata would leave cells mHM cannot read inside the
    model domain. Anything unrecognised keeps the historical nodata padding.
    """
    stem = path.stem.lower()
    for prefixes, pad_value, integer in PAD_SPECS:
        if stem.startswith(prefixes):
            return pad_value, integer
    return None, True


def staged_files(project_folder) -> list[Path]:
    """Return the files Execute All staged for publication, rasters first."""
    staging = Path(morph_staging_folder(project_folder))
    if not staging.is_dir():
        return []
    files = [path for path in sorted(staging.iterdir()) if path.is_file()]
    return (
        [path for path in files if path.suffix.lower() in RASTER_SUFFIXES]
        + [path for path in files if path.suffix.lower() not in RASTER_SUFFIXES]
    )


def publish_model_inputs(
        project_folder,
        l0_header: Mapping[str, Any],
        *,
        log=None) -> list[Path]:
    """Pad every staged raster onto L0 and move the staged set into data/master."""
    master = Path(morph_folder(project_folder))
    master.mkdir(parents=True, exist_ok=True)
    published: list[Path] = []

    for path in staged_files(project_folder):
        suffix = path.suffix.lower()
        if suffix in RASTER_SUFFIXES:
            pad_value, integer = pad_spec_for(path)
            pad_l0_file_to_header(
                path,
                l0_header,
                nodata_value=-9999,
                integer=integer,
                pad_value=pad_value,
            )
            if log:
                padding = "nodata" if pad_value is None else str(pad_value)
                log(
                    f"{path.name} placed on the common L0 grid "
                    f"(padded with {padding})."
                )
        elif suffix not in CARRIED_SUFFIXES:
            if log:
                log(f"Skipping unexpected staged file: {path.name}")
            continue
        target = master / path.name
        if target.exists():
            target.unlink()
        shutil.move(str(path), str(target))
        published.append(target)
        if log:
            log(f"Published {target.name} to {master}.")

    _repoint_settings(project_folder, published, log=log)
    return published


def _repoint_settings(project_folder, published, log=None) -> None:
    """Rewrite the namelist handoff so every path names the published copy.

    Anything already sitting in `data/master` is repointed too, not just what
    this run moved, so a reused stage whose output was published earlier still
    ends up with a resolvable path.
    """
    master = Path(morph_folder(project_folder))
    by_name = {path.name: path for path in published}
    if master.is_dir():
        for path in master.iterdir():
            if path.is_file():
                by_name.setdefault(path.name, path)
    settings = load_settings(project_folder)
    changed = False

    for section in SECTIONS:
        record = settings.get(section)
        if not isinstance(record, dict):
            continue
        for key in PATH_KEYS:
            value = str(record.get(key, "") or "")
            target = by_name.get(Path(value).name) if value else None
            if target is None:
                continue
            updated = relative_workspace_path(project_folder, target)
            if updated != value:
                record[key] = updated
                changed = True
        scenes = record.get("scenes")
        if isinstance(scenes, list):
            for scene in scenes:
                if not isinstance(scene, dict):
                    continue
                value = str(scene.get("output_path", "") or "")
                target = by_name.get(Path(value).name) if value else None
                if target is None:
                    continue
                updated = relative_workspace_path(project_folder, target)
                if updated != value:
                    scene["output_path"] = updated
                    changed = True

    if changed:
        save_settings(project_folder, settings)
        if log:
            log("nml-settings.json now points at the published model inputs.")


def missing_model_inputs(project_folder) -> list[str]:
    """Return namelist paths that do not resolve to an existing file."""
    workspace = Path(workspace_folder(project_folder))
    settings = load_settings(project_folder)
    missing: list[str] = []

    def check(value):
        text = str(value or "")
        if not text:
            return
        path = Path(text)
        if not path.is_absolute():
            path = workspace / path
        if not path.is_file():
            missing.append(text)

    for section in SECTIONS + ("lai",):
        record = settings.get(section)
        if not isinstance(record, dict):
            continue
        for key in PATH_KEYS:
            check(record.get(key))
        scenes = record.get("scenes")
        if isinstance(scenes, list):
            for scene in scenes:
                if isinstance(scene, dict):
                    check(scene.get("output_path"))
    return missing


__all__ = [
    "missing_model_inputs",
    "pad_spec_for",
    "publish_model_inputs",
    "staged_files",
]
