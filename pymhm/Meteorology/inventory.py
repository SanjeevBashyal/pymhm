# -*- coding: utf-8 -*-
"""ERA5-Land inventory reporting for the meteorology workflow."""
from __future__ import annotations

from typing import Callable, Any


def log_inventory(
        inventory: dict[str, Any],
        log_message: Callable[[str], None]) -> None:
    """Log a concise inventory of ERA5-Land monthly files."""
    log_message("ERA5-Land file inventory:")
    for variable, count in inventory.get("counts", {}).items():
        first = inventory.get("first", {}).get(variable, "")
        last = inventory.get("last", {}).get(variable, "")
        period = f" | {first} -> {last}" if first and last else ""
        log_message(f"  {variable}: {count} file(s){period}")
