# -*- coding: utf-8 -*-
"""Which layout a project uses, and the folders that follow from it.

v5.13 and v6 disagree about where morphology, LAI and observations live, so
these helpers read the version out of the saved input state rather than making
every caller thread it through.

That read is cached per project folder. `morph_folder` has eleven importers and
is called inside loops; uncached, each call opened and parsed
`mhm_qgis_input_state.json`.
"""
from __future__ import annotations

import os

from ..file import json as jsonio
from .paths import (data_folder, geometry_folder, master_data_folder,
                    meteo_folder, morph_staging_folder,
                    morphology_staging_folder, output_folder, plugin_root,
                    raw_meteo_folder, restart_folder, static_folder,
                    version_key, workspace_folder, z_temp_folder,
                    data_raw_folder, domain_data_folder)

INPUT_STATE_FILENAME = "mhm_qgis_input_state.json"

#: project folder -> version key, with the state file's mtime as the guard.
_VERSION_CACHE: dict[str, tuple[float, str]] = {}


def input_state_path(project_folder) -> str:
    """Return the saved input state path for one project."""
    return os.path.join(workspace_folder(project_folder), INPUT_STATE_FILENAME)


def project_version(project_folder) -> str:
    """Return the mHM version the project was set up with.

    Read from the saved input state so the layout helpers can branch on version
    without every caller having to thread it through. Defaults to v5.13, which
    is the layout every existing project already uses.
    """
    path = input_state_path(project_folder)
    try:
        stamp = os.path.getmtime(path)
    except OSError:
        return "v5.13"
    cached = _VERSION_CACHE.get(str(project_folder))
    if cached is not None and cached[0] == stamp:
        return cached[1]
    state = jsonio.read(path)
    version = version_key(state.get("mhm_version")) if isinstance(state, dict) else "v5.13"
    _VERSION_CACHE[str(project_folder)] = (stamp, version)
    return version


def forget_project_version(project_folder=None) -> None:
    """Drop the cached version, for when the input state has just been written."""
    if project_folder is None:
        _VERSION_CACHE.clear()
    else:
        _VERSION_CACHE.pop(str(project_folder), None)


def is_v6(project_folder) -> bool:
    """Return True when the project uses the v6 input layout."""
    return project_version(project_folder) == "v6"


def project_template_root() -> str:
    """Return the project-template directory.

    A packaged plugin carries the templates inside itself; a source checkout
    keeps them beside the package at the repository root.
    """
    packaged = os.path.join(plugin_root(), "project-template")
    if os.path.isdir(packaged):
        return packaged
    # A source checkout keeps the templates under .scripts/, not beside the
    # package -- without this the mirror phase silently did nothing.
    repository = os.path.dirname(os.path.dirname(plugin_root()))
    scripted = os.path.join(repository, ".scripts", "project-template")
    if os.path.isdir(scripted):
        return scripted
    repository = os.path.dirname(os.path.dirname(plugin_root()))
    return os.path.join(repository, "project-template")


def project_template_dir(version_text: str | None) -> str | None:
    """Return the best matching project-template directory."""
    root = project_template_root()
    preferred = os.path.join(root, version_key(version_text))
    if os.path.isdir(preferred):
        return preferred

    fallback = os.path.join(root, "v5.13")
    if os.path.isdir(fallback):
        return fallback
    return None


def domain_dem_path(project_folder, outlet_id) -> str:
    """Return the per-domain DEM: NetCDF for v6, Arc/Info ASCII for v5.13."""
    name = "dem.nc" if is_v6(project_folder) else "dem.asc"
    return os.path.join(domain_data_folder(project_folder, outlet_id), name)


def morph_folder(project_folder, version: str | None = None) -> str:
    """Return the static morphology data folder."""
    if (version or project_version(project_folder)) == "v6":
        return os.path.join(master_data_folder(project_folder), "morph")
    return os.path.join(static_folder(project_folder), "morph")


def lai_folder(project_folder, version: str | None = None) -> str:
    """Return the LAI data folder.

    v6 carries LAI inside the same morphology input, so the two coincide.
    """
    resolved = version or project_version(project_folder)
    if resolved == "v6":
        return morph_folder(project_folder, version=resolved)
    return os.path.join(master_data_folder(project_folder), "lai")


def streamflow_observation_folder(project_folder, version: str | None = None) -> str:
    """Return the streamflow observation folder."""
    if (version or project_version(project_folder)) == "v6":
        return os.path.join(
            master_data_folder(project_folder), "gauge", "streamflow"
        )
    return os.path.join(
        master_data_folder(project_folder), "observation", "streamflow"
    )


def ensure_project_structure(project_folder, version_text=None) -> list[str]:
    """
    Create the plugin workspace structure from project-template.

    Only directories are created. Placeholder files are intentionally skipped so
    processing code does not mistake empty template files for prepared outputs.
    """
    created = []
    os.makedirs(str(project_folder), exist_ok=True)
    workspace = workspace_folder(project_folder)
    if not os.path.isdir(workspace):
        os.makedirs(workspace, exist_ok=True)
        created.append(workspace)
    os.makedirs(z_temp_folder(project_folder), exist_ok=True)
    os.makedirs(geometry_folder(project_folder), exist_ok=True)
    os.makedirs(morphology_staging_folder(project_folder), exist_ok=True)
    os.makedirs(morph_staging_folder(project_folder), exist_ok=True)

    # The caller's version wins: on a new project the state file does not exist
    # yet, so reading the version back would always answer v5.13 and build the
    # wrong tree for a v6 project.
    version = version_key(version_text) if version_text else project_version(
        project_folder)

    template = project_template_dir(version_text)
    if template:
        for root, dirs, _ in os.walk(template):
            dirs[:] = [
                name
                for name in dirs
                if not (name.startswith("<") and name.endswith(">"))
            ]
            for dirname in dirs:
                src = os.path.join(root, dirname)
                rel = os.path.relpath(src, template)
                dst = os.path.join(workspace, rel)
                if not os.path.isdir(dst):
                    os.makedirs(dst, exist_ok=True)
                    created.append(dst)

    for path in (
            data_folder(project_folder),
            master_data_folder(project_folder),
            data_raw_folder(project_folder),
            static_folder(project_folder),
            morph_folder(project_folder, version=version),
            meteo_folder(project_folder),
            raw_meteo_folder(project_folder),
            lai_folder(project_folder, version=version),
            streamflow_observation_folder(project_folder, version=version),
            output_folder(project_folder),
            restart_folder(project_folder)):
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
            created.append(path)

    return created
