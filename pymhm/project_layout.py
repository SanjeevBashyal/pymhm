# -*- coding: utf-8 -*-
"""Project folder layout helpers for pymhm."""

from __future__ import annotations

import os


WORKSPACE_FOLDER_NAME = "mhm-plugin"


def plugin_root() -> str:
    """Return the plugin package root."""
    return os.path.dirname(os.path.abspath(__file__))


def workspace_folder(project_folder) -> str:
    """Return the plugin-owned workspace inside the selected project folder."""
    return os.path.join(str(project_folder), WORKSPACE_FOLDER_NAME)


def version_key(version_text: str | None) -> str:
    """Return the normalized version folder name used by templates."""
    text = (version_text or "").strip().lower()
    if text.startswith("5"):
        return "v5.13"
    return "v6"


def project_template_dir(version_text: str | None) -> str | None:
    """Return the best matching project-template directory."""
    root = os.path.join(plugin_root(), "project-template")
    preferred = os.path.join(root, version_key(version_text))
    if os.path.isdir(preferred):
        return preferred

    fallback = os.path.join(root, "v5.13")
    if os.path.isdir(fallback):
        return fallback
    return None


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


def domain_dem_path(project_folder, outlet_id) -> str:
    return os.path.join(domain_data_folder(project_folder, outlet_id), "dem.asc")


def data_raw_folder(project_folder) -> str:
    """Return the raw input data folder."""
    return os.path.join(workspace_folder(project_folder), "data_raw")


def static_folder(project_folder) -> str:
    """Return the shared static data folder."""
    return os.path.join(master_data_folder(project_folder), "static")


def morph_folder(project_folder) -> str:
    """Return the static morphology data folder."""
    return os.path.join(static_folder(project_folder), "morph")


def meteo_folder(project_folder) -> str:
    """Return the meteorology forcing data folder."""
    return os.path.join(master_data_folder(project_folder), "meteo")


def raw_meteo_folder(project_folder) -> str:
    """Return the raw meteorology input data folder."""
    return os.path.join(data_raw_folder(project_folder), "meteo")


def lai_folder(project_folder) -> str:
    """Return the LAI data folder."""
    return os.path.join(master_data_folder(project_folder), "lai")


def streamflow_observation_folder(project_folder) -> str:
    """Return the streamflow observation folder."""
    return os.path.join(
        master_data_folder(project_folder), "observation", "streamflow"
    )


def output_folder(project_folder) -> str:
    """Return the mHM output folder."""
    return os.path.join(workspace_folder(project_folder), "output")


def restart_folder(project_folder) -> str:
    """Return the mHM restart folder."""
    return os.path.join(workspace_folder(project_folder), "restart")


def relative_project_path(project_folder, path) -> str:
    """Return a forward-slash path relative to the selected project root."""
    return os.path.relpath(str(path), str(project_folder)).replace("\\", "/")


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
            morph_folder(project_folder),
            meteo_folder(project_folder),
            raw_meteo_folder(project_folder),
            lai_folder(project_folder),
            streamflow_observation_folder(project_folder),
            output_folder(project_folder),
            restart_folder(project_folder)):
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
            created.append(path)

    return created
