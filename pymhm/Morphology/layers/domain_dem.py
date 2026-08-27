# -*- coding: utf-8 -*-
"""Write every domain DEM on the common L0 grid.

Delineation records each domain's polygon but not its DEM, because the model
extent is not known until meteorology has fixed the L2 grid. Morphology Setup
then masks the cropped L0 DEM with each polygon, so every `data/<domain>/dem.asc`
has the same matrix size as the shared morphology and meteorology inputs.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any, Mapping

from ..file_tasks import write_domain_dem_ascii
from .domain_dem_nc import write_domain_dem_netcdf
from ..watershed.domain_state import (
    active_domain_records,
    load_state as load_domain_state,
    resolve_output_path,
)
from ...project_layout import (
    domain_data_folder,
    domain_dem_path,
    is_v6,
    morph_folder,
)


# mHM reads its morphology per domain from `dir_Morpho`, so each domain needs
# the shared rasters beside its own dem.asc. Files are listed rather than
# globbed so an unrelated file in data/master never leaks into a domain, and
# they are grouped so that a group with several spellings -- v5 writes a soil
# class raster where v6 writes a horizon cube -- only counts as missing when
# none of them is there. Land cover, LAI, and idgauges are absent on purpose:
# the namelist points every domain at the shared copy under data/master.
DOMAIN_INPUT_GROUPS = (
    ("slope", ("slope.asc",)),
    ("aspect", ("aspect.asc",)),
    ("flow accumulation", ("facc.asc",)),
    ("flow direction", ("fdir.asc",)),
    ("soil class", ("soil_class.asc", "soil_horizon_class.nc")),
    (
        "soil class definition",
        ("soil_classdefinition.txt", "soil_classdefinition_iFlag_soilDB_1.txt"),
    ),
    ("geology class", ("geology_class.asc",)),
    ("geology class definition", ("geology_classdefinition.txt",)),
)
DOMAIN_INPUT_NAMES = tuple(
    name for _label, names in DOMAIN_INPUT_GROUPS for name in names
)
SIDECAR_SUFFIXES = (".prj",)


def copy_master_inputs(project_folder, directory, log=None) -> list[Path]:
    """Copy the shared morphology inputs from data/master into one domain.

    Run this after `advanced_l0.publish_model_inputs`: an advanced or
    mHM-ready soil, geology, or land-cover raster only reaches data/master at
    publication, so copying earlier would silently find nothing. An unchanged
    destination is left alone, so re-running Morphology Setup does not recopy
    gigabytes.
    """
    master = Path(morph_folder(project_folder))
    destination = Path(directory)
    destination.mkdir(parents=True, exist_ok=True)
    copied: list[Path] = []
    absent: list[str] = []

    for label, names in DOMAIN_INPUT_GROUPS:
        present = [master / name for name in names if (master / name).is_file()]
        if not present:
            absent.append(label)
            continue
        candidates = list(present)
        candidates.extend(
            source.with_suffix(suffix)
            for source in present
            for suffix in SIDECAR_SUFFIXES
        )
        for source in candidates:
            if not source.is_file():
                continue
            target = destination / source.name
            if _already_copied(source, target):
                continue
            temporary = target.with_name(f".{target.name}.tmp")
            temporary.unlink(missing_ok=True)
            try:
                shutil.copyfile(source, temporary)
                os.replace(temporary, target)
            except BaseException:
                temporary.unlink(missing_ok=True)
                raise
            copied.append(target)
    if log:
        log(
            f"Copied {len(copied)} shared morphology file(s) into "
            f"{destination}."
        )
        if absent:
            # A silent skip once left every domain with only its dem.asc, so
            # name what data/master could not supply.
            log(
                f"WARNING: {master} has no "
                + ", ".join(absent)
                + f"; {destination.name} will run without it."
            )
    return copied


def _already_copied(source: Path, target: Path) -> bool:
    """Return True when the destination already holds this exact content."""
    if not target.is_file():
        return False
    source_stat, target_stat = source.stat(), target.stat()
    return (
        source_stat.st_size == target_stat.st_size
        and target_stat.st_mtime_ns >= source_stat.st_mtime_ns
    )


def domain_dem_plan(project_folder) -> list[dict[str, Any]]:
    """Return one entry per active domain: its polygon and its target DEM."""
    state = load_domain_state(project_folder)
    plan: list[dict[str, Any]] = []
    for record in active_domain_records(state):
        if record.get("domain_id") is None:
            continue
        outlet_id = str(record.get("outlet_id", ""))
        name = "dem_extent" if record.get("is_dem_domain") else outlet_id
        polygon = _domain_polygon(project_folder, record)
        plan.append({
            "domain_id": int(record["domain_id"]),
            "outlet_id": outlet_id,
            "name": name,
            "polygon": polygon,
            "directory": domain_data_folder(project_folder, name),
            "dem_path": domain_dem_path(project_folder, name),
        })
    return plan


def _domain_polygon(project_folder, record) -> str:
    """Return the delineated polygon path for one domain."""
    if record.get("is_dem_domain"):
        from ...project_layout import geometry_folder

        return os.path.join(
            geometry_folder(project_folder),
            "Watersheds",
            "DomainDelineation",
            "4_watershed_DEM.shp",
        )
    value = record.get("vector_path")
    if not value:
        return ""
    return str(resolve_output_path(project_folder, value))


def write_domain_dems(
        project_folder,
        cropped_dem: str,
        l0_header: Mapping[str, Any],
        *,
        reference_path=None,
        log=None) -> list[dict[str, Any]]:
    """Mask the cropped L0 DEM with each domain polygon and write its dem.asc."""
    plan = domain_dem_plan(project_folder)
    if not plan:
        return []
    if not Path(cropped_dem).is_file():
        raise FileNotFoundError(
            f"The cropped L0 DEM is required before writing domain DEMs: "
            f"{cropped_dem}"
        )

    written: list[dict[str, Any]] = []
    for entry in plan:
        polygon = entry["polygon"]
        if not polygon or not Path(polygon).is_file():
            raise FileNotFoundError(
                f"Domain {entry['name']} has no delineated polygon: {polygon}"
            )
        os.makedirs(entry["directory"], exist_ok=True)
        if log:
            log(
                f"Writing domain {entry['name']} DEM on the common L0 grid "
                f"({int(l0_header['ncols'])} x {int(l0_header['nrows'])} cells)."
            )
        writer = (
            write_domain_dem_netcdf
            if str(entry["dem_path"]).endswith(".nc")
            else write_domain_dem_ascii
        )
        writer(
            cropped_dem,
            entry["dem_path"],
            l0_header,
            polygon,
            reference_path=reference_path,
        )
        # v6 reads every shared layer from master's input.nc, so a domain folder
        # holds only its own DEM.
        entry["copied"] = [] if is_v6(project_folder) else [
            str(path)
            for path in copy_master_inputs(
                project_folder, entry["directory"], log=log)
        ]
        written.append(entry)
        if log:
            log(f"Domain {entry['name']} DEM written: {entry['dem_path']}")
    return written


__all__ = [
    "DOMAIN_INPUT_GROUPS",
    "DOMAIN_INPUT_NAMES",
    "copy_master_inputs",
    "domain_dem_plan",
    "write_domain_dems",
]
