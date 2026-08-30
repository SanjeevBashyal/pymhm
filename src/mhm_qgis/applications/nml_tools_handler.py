# -*- coding: utf-8 -*-
"""Everything the plugin asks of nml-tools.

One call crosses this boundary: opening the namelist editor. The editor is
launched in-process and runs modally on the QGIS event loop, so it returns only
once the user closes it.
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any, Mapping

LogCallback = Callable[[str], None]


def launch_profile_editor(
    schemas_dir: str | Path,
    output_dir: str | Path,
    initial_values: Mapping[str, Any] | None = None,
    initial_dimensions: Mapping[str, int] | None = None,
    *,
    log: LogCallback | None = None,
) -> int:
    """Open the nml-tools namelist editor over one project.

    `configure_qtpy_api()` must already have run: nml-tools loads Qt on the
    first call, and the API choice has to be settled before that.
    """
    from .mhm_tools_handler import capture_messages

    # Looked up on the module rather than bound at import time, so the editor
    # stays patchable and nml-tools keeps deferring its Qt import until called.
    from nml_tools import gui

    with capture_messages(log, logger_names=("nml_tools",)):
        return gui.launch_gui(
            schemas_dir=schemas_dir,
            output_dir=output_dir,
            initial_values=initial_values,
            initial_dimensions=initial_dimensions,
        )


__all__ = ["launch_profile_editor"]
