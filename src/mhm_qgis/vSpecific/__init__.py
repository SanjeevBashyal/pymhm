"""Version-specific nml-tools initial values."""

from __future__ import annotations

from typing import Any


def _adapter(version_text: str):
    if str(version_text or "").strip().startswith("5"):
        from . import v5_13 as adapter
    else:
        from . import v6 as adapter
    return adapter


def build_dimensions(version_text: str, dialog: Any) -> dict[str, int]:
    """Return nml-tools runtime dimensions for the selected mHM version."""
    return _adapter(version_text).build_dimensions(dialog)


def build_initial_values(version_text: str, dialog: Any) -> dict[str, Any]:
    """Return the profile/namelist/field overlay accepted by ``launch_gui``."""
    return _adapter(version_text).build_initial_values(dialog)


__all__ = [
    "build_dimensions",
    "build_initial_values",
]
