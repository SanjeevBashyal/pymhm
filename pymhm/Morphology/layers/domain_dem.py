# -*- coding: utf-8 -*-
"""Write every domain DEM on the common L0 grid.

Delineation records each domain's polygon but not its DEM, because the model
extent is not known until meteorology has fixed the L2 grid. Morphology Setup
then masks the cropped L0 DEM with each polygon, so every `data/<domain>/dem.asc`
has the same matrix size as the shared morphology and meteorology inputs.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Mapping

from ..file_tasks import write_domain_dem_ascii
from ..watershed.domain_state import (
    active_domain_records,
    load_state as load_domain_state,
    resolve_output_path,
)
from ...project_layout import domain_data_folder, domain_dem_path


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
        write_domain_dem_ascii(
            cropped_dem,
            entry["dem_path"],
            l0_header,
            polygon,
            reference_path=reference_path,
        )
        written.append(entry)
        if log:
            log(f"Domain {entry['name']} DEM written: {entry['dem_path']}")
    return written


__all__ = ["domain_dem_plan", "write_domain_dems"]
