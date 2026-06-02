# -*- coding: utf-8 -*-
"""Hydrology observation output paths."""
from __future__ import annotations

from pathlib import Path


def streamflow_observation_folder(project_folder) -> Path:
    """Return the project-local mHM streamflow observation folder."""
    return Path(project_folder) / "input" / "observations" / "streamflow"
