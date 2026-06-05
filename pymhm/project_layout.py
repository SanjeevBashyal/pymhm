# -*- coding: utf-8 -*-
"""Project folder layout helpers for pymhm."""

from __future__ import annotations

import os


def plugin_root() -> str:
    """Return the plugin package root."""
    return os.path.dirname(os.path.abspath(__file__))


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
    return os.path.join(str(project_folder), "Z Temp")


def geometry_folder(project_folder) -> str:
    """Return the temporary Geometry folder."""
    return os.path.join(z_temp_folder(project_folder), "Geometry")


def data_folder(project_folder) -> str:
    """Return the mHM data folder."""
    return os.path.join(str(project_folder), "data")


def static_folder(project_folder) -> str:
    """Return the static data folder."""
    return os.path.join(data_folder(project_folder), "static")


def morph_folder(project_folder) -> str:
    """Return the static morphology data folder."""
    return os.path.join(static_folder(project_folder), "morph")


def meteo_folder(project_folder) -> str:
    """Return the meteorology forcing data folder."""
    return os.path.join(data_folder(project_folder), "meteo")


def lai_folder(project_folder) -> str:
    """Return the LAI data folder."""
    return os.path.join(data_folder(project_folder), "lai")


def streamflow_observation_folder(project_folder) -> str:
    """Return the streamflow observation folder."""
    return os.path.join(data_folder(project_folder), "observation", "streamflow")


def output_folder(project_folder) -> str:
    """Return the mHM output folder."""
    return os.path.join(str(project_folder), "output")


def restart_folder(project_folder) -> str:
    """Return the mHM restart folder."""
    return os.path.join(str(project_folder), "restart")


def relative_project_path(project_folder, path) -> str:
    """Return a forward-slash path relative to the project root."""
    return os.path.relpath(str(path), str(project_folder)).replace("\\", "/")


def ensure_project_structure(project_folder, version_text=None) -> list[str]:
    """
    Create the project directory structure from project-template.

    Only directories are created. Placeholder files are intentionally skipped so
    processing code does not mistake empty template files for prepared outputs.
    """
    created = []
    os.makedirs(str(project_folder), exist_ok=True)
    os.makedirs(z_temp_folder(project_folder), exist_ok=True)
    os.makedirs(geometry_folder(project_folder), exist_ok=True)

    template = project_template_dir(version_text)
    if template:
        for root, dirs, _ in os.walk(template):
            for dirname in dirs:
                src = os.path.join(root, dirname)
                rel = os.path.relpath(src, template)
                dst = os.path.join(str(project_folder), rel)
                if not os.path.isdir(dst):
                    os.makedirs(dst, exist_ok=True)
                    created.append(dst)

    for path in (
            data_folder(project_folder),
            static_folder(project_folder),
            morph_folder(project_folder),
            meteo_folder(project_folder),
            lai_folder(project_folder),
            streamflow_observation_folder(project_folder),
            output_folder(project_folder),
            restart_folder(project_folder)):
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
            created.append(path)

    return created
