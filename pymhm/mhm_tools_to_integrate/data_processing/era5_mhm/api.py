"""Public mHM forcing conversion API."""
from __future__ import annotations

from pathlib import Path
from typing import Callable

from .build import build_daily_dataset
from .dependencies import import_dependencies
from .inventory import files_for_variable
from .io import write_header, write_netcdf
from .logging import log_message
from .specs import FORCING_SPECS
from .types import MeteoForcingResult
from .validation import existing_output_is_valid


def process_era5_to_mhm(
    nc_folder,
    output_root,
    bounds: tuple[float, float, float, float] | None = None,
    target_lat=None,
    target_lon=None,
    target_header: dict | None = None,
    skip_existing: bool = True,
    log: Callable[[str], None] | None = None,
) -> MeteoForcingResult:
    """
    Convert ERA5-Land monthly NetCDF files to mHM meteorology forcing files.

    Parameters
    ----------
    nc_folder:
        Folder containing files named ``ERA5_Land_<variable>_<year>_<month>.nc``.
    output_root:
        Folder that will receive ``pre/pre.nc``, ``tavg/tavg.nc``,
        ``tmin/tmin.nc``, and ``tmax/tmax.nc``.
    bounds:
        Optional WGS84 bounds as ``(west, east, south, north)``. Bounds are
        snapped outward to the available ERA5-Land grid centers.
    target_lat, target_lon:
        Optional target coordinate axes used to nearest-neighbor resample the
        final daily fields.
    target_header:
        Optional projected mHM grid header to write beside each forcing file.
    skip_existing:
        Reuse an output when both ``<var>.nc`` and ``header.txt`` already exist.
    log:
        Optional logger callback.
    """
    np, _, xr = import_dependencies()

    nc_folder = Path(nc_folder)
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    outputs: dict[str, Path] = {}
    headers: dict[str, Path] = {}
    files_used: dict[str, int] = {}

    for spec in FORCING_SPECS:
        output_dir = output_root / spec.output_folder
        output_file = output_dir / f"{spec.output_variable}.nc"
        header_file = output_dir / "header.txt"

        if skip_existing and output_file.exists() and header_file.exists():
            if existing_output_is_valid(
                    output_file, spec.output_variable, np, xr, log):
                log_message(
                    log,
                    f"{spec.output_variable}: existing output found, reusing {output_file}",
                )
                outputs[spec.output_variable] = output_file
                headers[spec.output_variable] = header_file
                files_used[spec.output_variable] = 0
                continue

            log_message(
                log,
                f"{spec.output_variable}: existing output is invalid or from an older converter; regenerating {output_file}",
            )

        files = files_for_variable(nc_folder, spec.file_variable)
        if not files:
            raise FileNotFoundError(
                f"No ERA5-Land files found for {spec.file_variable} in {nc_folder}")

        output_dir.mkdir(parents=True, exist_ok=True)
        log_message(
            log,
            f"{spec.output_variable}: processing {len(files)} file(s) "
            f"from {files[0].name} to {files[-1].name}",
        )
        ds_out = build_daily_dataset(
            files,
            spec,
            bounds=bounds,
            target_lat=target_lat,
            target_lon=target_lon,
            log=log,
        )
        write_netcdf(ds_out, spec.output_variable, output_file)
        write_header(
            ds_out,
            spec.output_variable,
            header_file,
            header=target_header,
        )

        outputs[spec.output_variable] = output_file
        headers[spec.output_variable] = header_file
        files_used[spec.output_variable] = len(files)
        log_message(log, f"{spec.output_variable}: written {output_file}")

    return MeteoForcingResult(
        outputs=outputs,
        headers=headers,
        files_used=files_used,
    )
