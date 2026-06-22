# -*- coding: utf-8 -*-
"""Lat/lon setup creation adapter for mhm_tools."""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from .._bundled import ensure_bundled_mhm_tools
from ..logging import capture_messages


def create_latlon_file(
        out_file: str | Path,
        level0: dict | str | Path,
        level1,
        level11=None,
        level2=None,
        write_header_l0: str | Path | None = None,
        write_header_l1: str | Path | None = None,
        write_header_l11: str | Path | None = None,
        write_header_l2: str | Path | None = None,
        crs: str | None = None,
        dtype: str = "f4",
        compression: int = 9,
        add_bounds: bool = False,
        log: Callable[[str], None] | None = None) -> Path:
    """Create latlon.nc through mhm_tools and return the output path."""
    ensure_bundled_mhm_tools()
    from mhm_tools.setup_creation import create_latlon

    output = Path(out_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    with capture_messages(log):
        create_latlon(
            out_file=output,
            level0=level0,
            level1=level1,
            level11=level11,
            level2=level2,
            write_header_l0=write_header_l0,
            write_header_l1=write_header_l1,
            write_header_l11=write_header_l11,
            write_header_l2=write_header_l2,
            crs=crs,
            dtype=dtype,
            compression=compression,
            add_bounds=add_bounds,
        )
    return output


__all__ = ["create_latlon_file"]
