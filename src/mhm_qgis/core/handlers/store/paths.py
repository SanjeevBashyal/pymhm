# -*- coding: utf-8 -*-
"""Where things live: pure path construction, no I/O and no version branching.

Deliberately the leaf of the store. `core/handlers/state` imports only this, so
the arrow between state and store never reverses -- a version-aware path helper
here would have to read the input state, and the state modules need
`workspace_folder` to find that file in the first place.
"""
from __future__ import annotations

import os

#: The plugin owns exactly one directory inside a user's project.
WORKSPACE_FOLDER_NAME = "mhm-plugin"


def plugin_root() -> str:
    """Return the plugin package root.

    This module sits three levels inside the package, so the root is three
    directories up -- not this file's own directory.
    """
    here = os.path.dirname(os.path.abspath(__file__))       # .../core/handlers/store
    return os.path.dirname(os.path.dirname(os.path.dirname(here)))


def workspace_folder(project_folder) -> str:
    """Return the plugin-owned workspace inside the selected project folder."""
    return os.path.join(str(project_folder), WORKSPACE_FOLDER_NAME)


def version_key(version_text: str | None) -> str:
    """Return the normalized version folder name used by templates."""
    text = (version_text or "").strip().lower()
    if text.startswith("5"):
        return "v5.13"
    return "v6"


def z_temp_folder(project_folder) -> str:
    """Return the helper-file folder for temporary project artifacts."""
    return os.path.join(workspace_folder(project_folder), "Z Temp")


def geometry_folder(project_folder) -> str:
    """Return the temporary Geometry folder."""
    return os.path.join(z_temp_folder(project_folder), "Geometry")


def morphology_staging_folder(project_folder) -> str:
    """Return the staging folder for model-ready morphology inputs."""
    return os.path.join(z_temp_folder(project_folder), "Morphology")


def morph_staging_folder(project_folder) -> str:
    """Return the staging folder for land cover, soil, and geology outputs.

    Execute All writes here. Morphology Setup publishes into
    `data/master/static/morph` only after crop, mask, and write succeed.
    """
    return os.path.join(morphology_staging_folder(project_folder), "morph")


def lai_staging_folder(project_folder) -> str:
    """Return the staging folder for LAI crop and mask outputs."""
    return os.path.join(morphology_staging_folder(project_folder), "lai")


def lai_dem_staging_path(project_folder) -> str:
    """Return the staged LAI resampled onto the filled DEM grid."""
    return os.path.join(lai_staging_folder(project_folder), "lai_dem.nc")


def data_folder(project_folder) -> str:
    """Return the mHM data folder."""
    return os.path.join(workspace_folder(project_folder), "data")


def master_data_folder(project_folder) -> str:
    """Return the shared model-input folder used by every domain."""
    return os.path.join(data_folder(project_folder), "master")


def domain_data_folder(project_folder, outlet_id) -> str:
    """Return one outlet-named domain folder, rejecting unsafe IDs."""
    name = str(outlet_id).strip()
    unsafe = not name or name in {".", ".."} or any(
        char in name for char in ("/", "\\", "\0")
    )
    if unsafe:
        raise ValueError(f"Outlet ID cannot be used as a domain directory: {outlet_id!r}")
    return os.path.join(data_folder(project_folder), name)


def data_raw_folder(project_folder) -> str:
    """Return the raw input data folder."""
    return os.path.join(workspace_folder(project_folder), "data_raw")


def static_folder(project_folder) -> str:
    """Return the shared static data folder."""
    return os.path.join(master_data_folder(project_folder), "static")


def meteo_folder(project_folder) -> str:
    """Return the meteorology forcing data folder."""
    return os.path.join(master_data_folder(project_folder), "meteo")


def raw_meteo_folder(project_folder) -> str:
    """Return the raw meteorology input data folder."""
    return os.path.join(data_raw_folder(project_folder), "meteo")


def output_folder(project_folder) -> str:
    """Return the mHM output folder."""
    return os.path.join(workspace_folder(project_folder), "output")


def restart_folder(project_folder) -> str:
    """Return the mHM restart folder."""
    return os.path.join(workspace_folder(project_folder), "restart")


def relative_project_path(project_folder, path) -> str:
    """Return a forward-slash path relative to the selected project root."""
    return os.path.relpath(str(path), str(project_folder)).replace("\\", "/")
