# -*- coding: utf-8 -*-
"""Region-level data processing services independent of UI classes."""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

from .ascii_morphology import prepare_morphology_ascii_files
from .era5 import process_era5_to_mhm
from .._bundled import ensure_bundled_mhm_tools
from ..setup_creation.latlon import create_latlon_file


class RegionDataHandler:
    """Provide morphology, meteorology, ASCII, and NetCDF data to a region."""

    def __init__(
        self, region: Any | None = None, log: Callable[[str], None] | None = None
    ) -> None:
        self.region = region
        self.log = log

    def create_latlon(
        self,
        out_file: str | Path,
        level0: dict | str | Path,
        level1,
        level11=None,
        level2=None,
        crs: str | None = None,
    ) -> Path:
        """Create latlon.nc through mhm_tools."""
        if crs is None and self.region is not None:
            crs = getattr(self.region, "crs", None)
        result = create_latlon_file(
            out_file=out_file,
            level0=level0,
            level1=level1,
            level11=level11,
            level2=level2,
            crs=crs,
            log=self.log,
        )
        if self.region is not None:
            self.region.outputs["latlon"] = result
        return result

    def prepare_meteo_forcing(
        self,
        nc_folder: str | Path,
        output_root: str | Path,
        bounds: tuple[float, float, float, float] | None = None,
        target_lat=None,
        target_lon=None,
        target_header: dict | None = None,
        skip_existing: bool = True,
    ):
        """Prepare mHM meteorology forcing files from ERA5-Land inputs."""
        return process_era5_to_mhm(
            nc_folder=nc_folder,
            output_root=output_root,
            bounds=bounds,
            target_lat=target_lat,
            target_lon=target_lon,
            target_header=target_header,
            skip_existing=skip_existing,
            log=self.log,
        )

    def read_dataset(self, path: str | Path, **kwargs):
        """Read a NetCDF/raster-compatible dataset through mhm_tools."""
        ensure_bundled_mhm_tools()
        from mhm_tools.common.file_handler import get_xarray_ds_from_file

        return get_xarray_ds_from_file(path, **kwargs)

    def write_dataset(self, dataset, path: str | Path, **kwargs) -> Path:
        """Write an xarray dataset through mhm_tools."""
        ensure_bundled_mhm_tools()
        from mhm_tools.common.file_handler import write_xarray_to_file

        output = Path(path)
        write_xarray_to_file(dataset, output, **kwargs)
        if self.log:
            self.log(f"Dataset written: {output}")
        return output

    def write_ascii(self, dataset, path: str | Path, **kwargs) -> Path:
        """Write an xarray object to ASCII through mhm_tools."""
        ensure_bundled_mhm_tools()
        from mhm_tools.common.file_handler import write_xarray_to_ascii

        output = Path(path)
        write_xarray_to_ascii(dataset, output, **kwargs)
        if self.log:
            self.log(f"ASCII written: {output}")
        return output

    def write_header(self, header: dict, path: str | Path, **kwargs) -> Path:
        """Write an ESRI/mHM header through mhm_tools."""
        ensure_bundled_mhm_tools()
        from mhm_tools.common.esri_grid import write_header

        output = Path(path)
        write_header(output, header, **kwargs)
        if self.log:
            self.log(f"Header written: {output}")
        return output

    def prepare_morphology_ascii(
        self,
        layers,
        headers: dict,
        overwrite: bool = True,
    ):
        """Prepare static morphology ASCII files through mhm_tools."""
        result = prepare_morphology_ascii_files(
            layers=layers,
            headers=headers,
            overwrite=overwrite,
            log=self.log,
        )
        if self.region is not None:
            self.region.outputs.update(result.outputs)
            self.region.grid_headers.update(result.headers)
        return result
