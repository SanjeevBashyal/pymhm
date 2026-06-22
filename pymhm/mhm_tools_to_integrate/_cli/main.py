# -*- coding: utf-8 -*-
"""CLI entry points for mhm_tools_to_integrate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ..data_processing.era5 import process_era5_to_mhm
from ..setup_creation.latlon import create_latlon_file


def _level_value(text: str):
    """Parse a level argument as path, JSON dictionary, or cell size."""
    path = Path(text)
    if path.exists():
        return path
    try:
        return float(text)
    except ValueError:
        pass
    try:
        value = json.loads(text)
    except json.JSONDecodeError as e:
        raise argparse.ArgumentTypeError(
            f"Expected an existing path, number, or JSON dictionary: {text}"
        ) from e
    if not isinstance(value, dict):
        raise argparse.ArgumentTypeError(
            f"Expected JSON dictionary for grid header: {text}")
    return value


def _bounds_value(text: str) -> tuple[float, float, float, float]:
    """Parse west,east,south,north bounds."""
    parts = [float(part.strip()) for part in text.split(",")]
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            "Bounds must be west,east,south,north.")
    return parts[0], parts[1], parts[2], parts[3]


def _header_value(text: str) -> dict:
    """Parse target header from JSON string or JSON file."""
    path = Path(text)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    value = json.loads(text)
    if not isinstance(value, dict):
        raise argparse.ArgumentTypeError("Header must be a JSON dictionary.")
    return value


def _add_latlon_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "latlon",
        help="Create latlon.nc through mhm_tools without using the mhm_tools CLI.",
    )
    parser.add_argument("-D", "--level0", required=True, type=_level_value)
    parser.add_argument("-H", "--level1", required=True, type=_level_value)
    parser.add_argument("-R", "--level11", type=_level_value)
    parser.add_argument("-M", "--level2", type=_level_value)
    parser.add_argument("-c", "--crs")
    parser.add_argument("-o", "--out-file", default="latlon.nc")
    parser.add_argument("-d", "--dtype", default="f4")
    parser.add_argument("-x", "--compression", type=int, default=9)
    parser.add_argument("-b", "--add-bounds", action="store_true")
    parser.set_defaults(func=_run_latlon)


def _add_era5_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "era5-to-mhm",
        help="Prepare mHM forcing NetCDF files from ERA5-Land monthly NetCDFs.",
    )
    parser.add_argument("--nc-folder", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--bounds", type=_bounds_value)
    parser.add_argument("--target-header", type=_header_value)
    parser.add_argument("--skip-existing", action="store_true")
    parser.set_defaults(func=_run_era5_to_mhm)


def _run_latlon(args) -> int:
    create_latlon_file(
        out_file=args.out_file,
        level0=args.level0,
        level1=args.level1,
        level11=args.level11,
        level2=args.level2,
        crs=args.crs,
        dtype=args.dtype,
        compression=args.compression,
        add_bounds=args.add_bounds,
        log=print,
    )
    return 0


def _run_era5_to_mhm(args) -> int:
    process_era5_to_mhm(
        nc_folder=args.nc_folder,
        output_root=args.output_root,
        bounds=args.bounds,
        target_header=args.target_header,
        skip_existing=args.skip_existing,
        log=print,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser."""
    parser = argparse.ArgumentParser(
        prog="mhm-tools-to-integrate",
        description="UI-free pymhm core workflows.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_latlon_parser(subparsers)
    _add_era5_parser(subparsers)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the command line interface."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
