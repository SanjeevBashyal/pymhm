# -*- coding: utf-8 -*-
"""Result recording helpers for meteorology forcing outputs."""
from __future__ import annotations

from pathlib import Path
from typing import Callable, TYPE_CHECKING

from .state import MeteorologyOutputState

if TYPE_CHECKING:
    from .ERA5Land.mhm.types import MeteoForcingResult


def record_forcing_outputs(
        state: MeteorologyOutputState,
        result: MeteoForcingResult) -> None:
    """Persist prepared NetCDF and header paths in the project state file."""
    for variable, output_path in result.outputs.items():
        state.mark_prepared(output_path, Path(output_path).name)
        header_path = result.headers.get(variable)
        if header_path:
            state.mark_prepared(header_path, Path(header_path).name)


def log_result_summary(
        result: MeteoForcingResult,
        log_message: Callable[[str], None]) -> None:
    """Log a compact summary for the QGIS console."""
    log_message("Meteorology forcing preparation completed.")
    log_message(f"Prepared variables: {', '.join(sorted(result.outputs))}")
