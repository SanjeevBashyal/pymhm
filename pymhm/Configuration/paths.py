# -*- coding: utf-8 -*-
"""Path helpers for mHM configuration schemas and templates."""
from __future__ import annotations

import os
from typing import Union

from .constants import OUTPUT_FILE_OVERRIDES, OUTPUT_FILES, TEMPLATE_NAMES
from ..project_layout import version_key as layout_version_key

PathInput = Union[str, os.PathLike[str]]


def package_root() -> str:
    """Return the plugin package root directory."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def schemas_dir() -> str:
    """Return the directory containing namelist YAML schemas."""
    return os.path.join(package_root(), "nml-schemas")


def templates_root() -> str:
    """Return the directory containing namelist templates."""
    return os.path.join(package_root(), "nml-templates")


def version_key(version_text: str | None) -> str:
    """Map the version combo-box text to a template folder key."""
    return layout_version_key(version_text)


def template_dir(version_text: str) -> str:
    """Return the template directory for the selected mHM version."""
    return os.path.join(templates_root(), version_key(version_text))


def template_path(version_text: str, kind: str) -> str:
    """Return the first existing template path for a namelist kind."""
    directory = template_dir(version_text)
    for name in TEMPLATE_NAMES[kind]:
        path = os.path.join(directory, name)
        if os.path.exists(path):
            return path
    expected = ", ".join(TEMPLATE_NAMES[kind])
    raise FileNotFoundError(
        f"No template found for {kind!r} in {directory}. Expected: {expected}")


def namelist_output_dir(project_folder: PathInput) -> str:
    """Return the directory where root namelist files are written."""
    return str(project_folder)


def output_filename(kind: str, version_text: str | None = None) -> str:
    """Return the version-aware output namelist filename for a kind."""
    overrides = OUTPUT_FILE_OVERRIDES.get(version_key(version_text), {})
    return overrides.get(kind, OUTPUT_FILES[kind])


def output_path(
        project_folder: PathInput,
        kind: str,
        version_text: str | None = None) -> str:
    """Return the project-local output namelist path for a kind."""
    return os.path.join(
        namelist_output_dir(project_folder),
        output_filename(kind, version_text),
    )
