# -*- coding: utf-8 -*-
"""Path helpers for mHM configuration schemas and templates."""

import os

from .constants import OUTPUT_FILES, TEMPLATE_NAMES


def package_root():
    """Return the plugin package root directory."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def schemas_dir():
    """Return the directory containing namelist YAML schemas."""
    return os.path.join(package_root(), "nml-schemas")


def templates_root():
    """Return the directory containing namelist templates."""
    return os.path.join(package_root(), "nml-templates")


def version_key(version_text):
    """Map the version combo-box text to a template folder key."""
    text = (version_text or "").strip().lower()
    if text.startswith("5"):
        return "v5.13"
    return "v6"


def template_dir(version_text):
    """Return the template directory for the selected mHM version."""
    return os.path.join(templates_root(), version_key(version_text))


def template_path(version_text, kind):
    """Return the first existing template path for a namelist kind."""
    directory = template_dir(version_text)
    for name in TEMPLATE_NAMES[kind]:
        path = os.path.join(directory, name)
        if os.path.exists(path):
            return path
    expected = ", ".join(TEMPLATE_NAMES[kind])
    raise FileNotFoundError(
        f"No template found for {kind!r} in {directory}. Expected: {expected}")


def model_inputs_dir(project_folder):
    """Return the project-local model input directory."""
    return os.path.join(project_folder, "model_inputs")


def output_path(project_folder, kind):
    """Return the project-local output namelist path for a kind."""
    return os.path.join(model_inputs_dir(project_folder), OUTPUT_FILES[kind])
