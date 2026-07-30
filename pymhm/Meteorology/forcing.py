# -*- coding: utf-8 -*-
"""QGIS-free inspection and preparation of mHM meteorology forcing."""
from __future__ import annotations

from dataclasses import dataclass
import math
import os
from pathlib import Path
import shutil
import tempfile
from typing import Callable, Mapping, Sequence

from .ERA5Land.mhm.build import (
    build_daily_dataset,
)
from .ERA5Land.mhm.io import (
    write_header,
    write_netcdf,
)
from .ERA5Land.mhm.specs import (
    FORCING_SPECS,
)
from .ERA5Land.mhm.types import (
    ForcingSpec,
    MeteoForcingResult,
)


MHM_READY = "mHM ready"
ERA5LAND = "ERA5Land"
_KINDS = {"pre": "precipitation", "precipitation": "precipitation",
          "temp": "temperature", "temperature": "temperature", "pet": "pet"}
_ERA_ALIASES = {
    "precipitation": ("tp", "total_precipitation",
                      "Total_precipitation_surface_1_Hour_Accumulation"),
    "temperature": ("t2m", "2m_temperature",
                    "Temperature_height_above_ground"),
    "pet": ("pev", "potential_evaporation"),
}
_READY_FILES = {
    "temperature": ("tavg.nc", "tmin.nc", "tmax.nc"),
    "pet": ("pet.nc",),
}
_PET_SPEC = ForcingSpec(
    output_variable="pet",
    output_folder="pet",
    source_key="potential_evaporation",
    file_variable="potential_evaporation",
    aggregation="era5land_daily_total_negative",
    units="mm",
    long_name="daily potential evapotranspiration",
    scale=-1000.0,
)


@dataclass(frozen=True)
class MeteoFolderSpec:
    """One user-selected meteorology input folder."""

    kind: str
    folder: Path | str
    source: str
    crs: str | None = None

    def normalized(self) -> "MeteoFolderSpec":
        return MeteoFolderSpec(
            kind=normalize_kind(self.kind),
            folder=Path(self.folder),
            source=normalize_source(self.source),
            crs=self.crs,
        )


@dataclass(frozen=True)
class TargetGrid:
    """WGS84 cell-centre axes and the matching projected mHM header."""

    lon: Sequence[float]
    lat: Sequence[float]
    header: Mapping[str, float]
    crs: str | None = None
    sample_lon: object | None = None
    sample_lat: object | None = None

    def validate(self) -> None:
        if len(self.lon) == 0 or len(self.lat) == 0:
            raise ValueError("The target meteorology grid is empty.")
        if int(self.header["ncols"]) != len(self.lon):
            raise ValueError("Target longitude count does not match header ncols.")
        if int(self.header["nrows"]) != len(self.lat):
            raise ValueError("Target latitude count does not match header nrows.")
        if (self.sample_lon is None) != (self.sample_lat is None):
            raise ValueError(
                "Both target sample longitude and latitude grids are required."
            )
        if self.sample_lon is not None:
            import numpy as np

            shape = (int(self.header["nrows"]), int(self.header["ncols"]))
            if (
                    np.shape(self.sample_lon) != shape
                    or np.shape(self.sample_lat) != shape):
                raise ValueError(
                    "Target sample coordinates do not match the header shape."
                )


@dataclass(frozen=True)
class SpatialMetadata:
    """Spatial metadata shared by all NetCDF files in an input folder."""

    kind: str
    source: str
    files: tuple[Path, ...]
    shape: tuple[int, int]
    x_coordinate: str
    y_coordinate: str
    x_resolution: float
    y_resolution: float
    resolution: float
    bounds: tuple[float, float, float, float]
    crs: str | None
    unit: str


@dataclass(frozen=True)
class SpatialResolution:
    """A grid resolution expressed in one CRS."""

    resolution: float
    x_resolution: float
    y_resolution: float
    unit: str
    crs: str | None


@dataclass(frozen=True)
class _SpatialSignature:
    metadata: SpatialMetadata
    x_values: object
    y_values: object


def normalize_kind(kind: str) -> str:
    """Return the canonical meteorology input kind."""
    normalized = _KINDS.get(str(kind).strip().lower())
    if not normalized:
        raise ValueError(f"Unsupported meteorology data type: {kind}")
    return normalized


def normalize_source(source: str) -> str:
    """Normalize UI and legacy spellings of a meteorology data source."""
    compact = str(source).strip().lower().replace("_", "").replace("-", "")
    compact = compact.replace(" ", "")
    if compact == "mhmready":
        return MHM_READY
    if compact == "era5land":
        return ERA5LAND
    raise ValueError(f"Unsupported meteorology data source: {source}")


def inspect_meteo_folder(
        folder: Path | str,
        kind: str,
        source: str,
        crs_fallback: str | None = None) -> SpatialMetadata:
    """Validate one input folder and return its common spatial metadata."""
    spec = MeteoFolderSpec(kind, folder, source, crs_fallback).normalized()
    files = _files_for_spec(spec)
    signatures = [
        _inspect_file(path, spec, _ready_variable(path, spec))
        for path in files
    ]
    _validate_spatial_signatures(signatures)
    _validate_times(files, spec)
    return SpatialMetadata(
        **{
            **signatures[0].metadata.__dict__,
            "files": tuple(files),
        }
    )


def inspect_meteo_inputs(
        precipitation: MeteoFolderSpec,
        temperature: MeteoFolderSpec,
        pet: MeteoFolderSpec | None = None) -> dict[str, SpatialMetadata]:
    """Validate all selected meteorology inputs."""
    specs = [precipitation.normalized(), temperature.normalized()]
    if specs[0].kind != "precipitation":
        raise ValueError("The precipitation input must have kind='precipitation'.")
    if specs[1].kind != "temperature":
        raise ValueError("The temperature input must have kind='temperature'.")
    if pet is not None:
        pet = pet.normalized()
        if pet.kind != "pet":
            raise ValueError("The PET input must have kind='pet'.")
        specs.append(pet)

    return {
        spec.kind: inspect_meteo_folder(
            spec.folder, spec.kind, spec.source, spec.crs)
        for spec in specs
    }


def resolution_in_crs(
        metadata: SpatialMetadata,
        target_crs: str | None) -> SpatialResolution:
    """Express an inspected grid resolution in ``target_crs`` units."""
    source_crs = metadata.crs
    if not target_crs or not source_crs or _same_crs(source_crs, target_crs):
        return SpatialResolution(
            metadata.resolution,
            metadata.x_resolution,
            metadata.y_resolution,
            metadata.unit,
            target_crs or source_crs,
        )

    try:
        from pyproj import CRS, Transformer
    except Exception as exc:
        raise RuntimeError(
            "pyproj is required to transform meteorology resolution.") from exc

    west, east, south, north = metadata.bounds
    x0 = (west + east) / 2.0
    y0 = (south + north) / 2.0
    transform = Transformer.from_crs(source_crs, target_crs, always_xy=True)
    p0 = transform.transform(x0, y0)
    px = transform.transform(x0 + metadata.x_resolution, y0)
    py = transform.transform(x0, y0 + metadata.y_resolution)
    dx = math.hypot(px[0] - p0[0], px[1] - p0[1])
    dy = math.hypot(py[0] - p0[0], py[1] - p0[1])
    return SpatialResolution(
        resolution=(dx + dy) / 2.0,
        x_resolution=dx,
        y_resolution=dy,
        unit=_crs_unit(CRS.from_user_input(target_crs)),
        crs=target_crs,
    )


def process_meteo_inputs(
        precipitation: MeteoFolderSpec,
        temperature: MeteoFolderSpec,
        output_root: Path | str,
        target_grid: TargetGrid,
        pet: MeteoFolderSpec | None = None,
        log: Callable[[str], None] | None = None) -> MeteoForcingResult:
    """Compile, resample, crop, and write all selected meteorology inputs."""
    target_grid.validate()
    inspected = inspect_meteo_inputs(precipitation, temperature, pet)
    specs = {
        "precipitation": precipitation.normalized(),
        "temperature": temperature.normalized(),
    }
    if pet is not None:
        specs["pet"] = pet.normalized()

    final_output_root = Path(output_root)
    final_output_root.mkdir(parents=True, exist_ok=True)
    output_root = Path(tempfile.mkdtemp(
        prefix=".pymhm-meteo-",
        dir=final_output_root,
    ))
    outputs: dict[str, Path] = {}
    headers: dict[str, Path] = {}
    files_used: dict[str, int] = {}
    forcing_specs = {spec.output_variable: spec for spec in FORCING_SPECS}
    forcing_specs["pet"] = _PET_SPEC

    variables = {
        "precipitation": ("pre",),
        "temperature": ("tavg", "tmin", "tmax"),
        "pet": ("pet",),
    }
    try:
        for kind, folder_spec in specs.items():
            metadata = inspected[kind]
            _validate_target_coverage(metadata, target_grid)
            for variable in variables[kind]:
                if log:
                    log(f"{variable}: processing {len(metadata.files)} file(s)")
                if folder_spec.source == ERA5LAND:
                    ds_out = build_daily_dataset(
                        metadata.files,
                        forcing_specs[variable],
                        bounds=None,
                        target_lat=target_grid.lat,
                        target_lon=target_grid.lon,
                        target_sample_lat=target_grid.sample_lat,
                        target_sample_lon=target_grid.sample_lon,
                        log=log,
                    )
                else:
                    path = _ready_path(metadata.files, variable)
                    ds_out = _prepare_ready_dataset(
                        path, variable, folder_spec.crs, target_grid)
                if variable == "pre":
                    ds_out = _normalize_precipitation(ds_out)
                ds_out = _attach_target_coordinates(
                    ds_out,
                    variable,
                    target_grid,
                )

                output_dir = output_root / variable
                output_file = output_dir / f"{variable}.nc"
                header_file = output_dir / "header.txt"
                write_netcdf(ds_out, variable, output_file)
                write_header(
                    ds_out,
                    variable,
                    header_file,
                    header=dict(target_grid.header),
                )
                outputs[variable] = output_file
                headers[variable] = header_file
                files_used[variable] = (
                    len(metadata.files)
                    if folder_spec.source == ERA5LAND else 1
                )

        result = _publish_meteo_result(
            MeteoForcingResult(outputs, headers, files_used),
            final_output_root,
        )
        if log:
            for variable, output in result.outputs.items():
                log(f"{variable}: written {output}")
        return result
    finally:
        shutil.rmtree(output_root, ignore_errors=True)


def _publish_meteo_result(
        staged: MeteoForcingResult,
        output_root: Path) -> MeteoForcingResult:
    """Publish a complete staged result and roll back replacement failures."""
    staging_root = next(iter(staged.outputs.values())).parents[1]
    backup_root = staging_root / "backups"
    managed_variables = ("pre", "tavg", "tmin", "tmax", "pet")
    backups: list[tuple[Path, Path]] = []
    published: list[Path] = []

    try:
        for variable in managed_variables:
            for filename in (f"{variable}.nc", "header.txt"):
                final = output_root / variable / filename
                if not final.exists():
                    continue
                backup = backup_root / variable / filename
                backup.parent.mkdir(parents=True, exist_ok=True)
                os.replace(final, backup)
                backups.append((backup, final))

        final_outputs = {}
        final_headers = {}
        for variable, staged_output in staged.outputs.items():
            final_dir = output_root / variable
            final_dir.mkdir(parents=True, exist_ok=True)
            final_output = final_dir / f"{variable}.nc"
            final_header = final_dir / "header.txt"
            os.replace(staged_output, final_output)
            published.append(final_output)
            os.replace(staged.headers[variable], final_header)
            published.append(final_header)
            final_outputs[variable] = final_output
            final_headers[variable] = final_header

        return MeteoForcingResult(
            final_outputs,
            final_headers,
            dict(staged.files_used),
        )
    except Exception:
        for path in published:
            if path.exists():
                path.unlink()
        for backup, final in reversed(backups):
            if backup.exists():
                final.parent.mkdir(parents=True, exist_ok=True)
                os.replace(backup, final)
        raise


def _direct_netcdf_files(folder: Path) -> list[Path]:
    if not folder.exists() or not folder.is_dir():
        raise FileNotFoundError(f"Meteorology folder does not exist: {folder}")
    return sorted(
        path for path in folder.iterdir()
        if path.is_file() and path.suffix.lower() == ".nc")


def _files_for_spec(spec: MeteoFolderSpec) -> list[Path]:
    folder = Path(spec.folder)
    files = _direct_netcdf_files(folder)
    if spec.source == ERA5LAND:
        if not files:
            raise ValueError(
                f"{spec.kind}: ERA5Land folder must contain at least one .nc file.")
        return files

    if spec.kind == "precipitation":
        if len(files) != 1:
            raise ValueError(
                "mHM-ready precipitation requires exactly one direct .nc file; "
                f"found {len(files)}.")
        return files

    expected = set(_READY_FILES[spec.kind])
    actual = {path.name for path in files}
    if actual != expected:
        raise ValueError(
            f"mHM-ready {spec.kind} folder must contain exactly "
            f"{', '.join(_READY_FILES[spec.kind])}.")
    return [folder / name for name in _READY_FILES[spec.kind]]


def _ready_variable(path: Path, spec: MeteoFolderSpec) -> str | None:
    if spec.source == ERA5LAND:
        return None
    if spec.kind == "precipitation":
        return "pre"
    return path.stem


def _open_dataset(path: Path):
    import xarray as xr

    last_error = None
    for engine in ("netcdf4", "h5netcdf", "scipy", None):
        try:
            kwargs = {} if engine is None else {"engine": engine}
            return xr.open_dataset(
                path, decode_times=True, mask_and_scale=True, **kwargs)
        except Exception as exc:
            last_error = exc
    raise RuntimeError(f"Could not open NetCDF file {path}: {last_error}")


def _inspect_file(
        path: Path,
        spec: MeteoFolderSpec,
        expected_variable: str | None) -> _SpatialSignature:
    import numpy as np

    with _open_dataset(path) as ds:
        variable = _find_variable(ds, spec.kind, expected_variable)
        da = ds[variable]
        y_dim, x_dim = _spatial_dims(da)
        x_name, x_values = _axis_values(ds, x_dim, "x")
        y_name, y_values = _axis_values(ds, y_dim, "y")

        if x_values is None or y_values is None:
            lon_name, lon = _geographic_coordinate(ds, "lon", (y_dim, x_dim))
            lat_name, lat = _geographic_coordinate(ds, "lat", (y_dim, x_dim))
            if lon is None or lat is None:
                raise ValueError(
                    f"{path.name}: could not identify spatial coordinates.")
            x_name, y_name = lon_name, lat_name
            x_values, y_values = lon, lat

        x_values = np.asarray(x_values, dtype="float64")
        y_values = np.asarray(y_values, dtype="float64")
        if x_values.ndim == 1 and y_values.ndim == 1:
            dx = _axis_resolution(x_values, path, x_name)
            dy = _axis_resolution(y_values, path, y_name)
            bounds = (
                float(np.nanmin(x_values)), float(np.nanmax(x_values)),
                float(np.nanmin(y_values)), float(np.nanmax(y_values)))
            shape = (len(y_values), len(x_values))
        elif x_values.shape == y_values.shape == (
                da.sizes[y_dim], da.sizes[x_dim]):
            dx = _mesh_resolution(x_values, axis=1, path=path, name=x_name)
            dy = _mesh_resolution(y_values, axis=0, path=path, name=y_name)
            bounds = (
                float(np.nanmin(x_values)), float(np.nanmax(x_values)),
                float(np.nanmin(y_values)), float(np.nanmax(y_values)))
            shape = x_values.shape
        else:
            raise ValueError(f"{path.name}: spatial coordinates are not a 2-D grid.")

        crs = _dataset_crs(ds, da, x_name, y_name) or spec.crs
        unit = _coordinate_unit(ds, x_name, crs)
        metadata = SpatialMetadata(
            kind=spec.kind,
            source=spec.source,
            files=(path,),
            shape=shape,
            x_coordinate=x_name,
            y_coordinate=y_name,
            x_resolution=dx,
            y_resolution=dy,
            resolution=(dx + dy) / 2.0,
            bounds=bounds,
            crs=crs,
            unit=unit,
        )
        return _SpatialSignature(metadata, x_values, y_values)


def _find_variable(ds, kind: str, expected: str | None = None) -> str:
    if expected:
        if expected in ds.data_vars:
            return expected
        candidates = [
            name for name, variable in ds.data_vars.items()
            if "time" in variable.dims or "valid_time" in variable.dims
        ]
        if kind == "precipitation" and len(candidates) == 1:
            return candidates[0]
        raise ValueError(
            f"Expected variable '{expected}'; available variables: "
            f"{', '.join(ds.data_vars)}.")

    for alias in _ERA_ALIASES[kind]:
        if alias in ds.data_vars:
            return alias
    raise ValueError(
        f"No ERA5Land {kind} variable found; available variables: "
        f"{', '.join(ds.data_vars)}.")


def _spatial_dims(da) -> tuple[str, str]:
    dims = [dim for dim in da.dims if dim not in {"time", "valid_time", "bnds"}]
    extra = [dim for dim in dims if da.sizes[dim] > 1]
    if len(extra) != 2:
        raise ValueError(
            f"Variable {da.name!r} must have exactly two spatial dimensions.")
    return extra[-2], extra[-1]


def _axis_values(ds, dim: str, axis: str):
    import numpy as np

    candidates = (
        ("lon", "longitude", "xc", "x", "easting")
        if axis == "x" else
        ("lat", "latitude", "yc", "y", "northing")
    )
    for name in (dim,) + candidates:
        if name not in ds.variables:
            continue
        values = np.asarray(ds[name].values)
        if values.ndim == 1 and ds[name].dims == (dim,):
            return name, values
    return dim, None


def _geographic_coordinate(ds, axis: str, dims: tuple[str, str]):
    aliases = ("lon", "longitude") if axis == "lon" else ("lat", "latitude")
    for name in aliases:
        if name in ds.variables and ds[name].dims == dims:
            return name, ds[name].values
    return axis, None


def _axis_resolution(values, path: Path, name: str) -> float:
    import numpy as np

    unique = np.unique(values[np.isfinite(values)])
    if len(unique) < 2:
        raise ValueError(
            f"{path.name}: coordinate {name} needs at least two values.")
    if len(unique) != len(values):
        raise ValueError(f"{path.name}: coordinate {name} contains duplicates.")
    diffs = np.abs(np.diff(np.sort(unique)))
    resolution = float(np.nanmedian(diffs))
    if not np.allclose(diffs, resolution, rtol=1e-6, atol=1e-9):
        raise ValueError(f"{path.name}: coordinate {name} is not regularly spaced.")
    return resolution


def _mesh_resolution(values, axis: int, path: Path, name: str) -> float:
    import numpy as np

    diffs = np.abs(np.diff(values, axis=axis))
    diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
    if not diffs.size:
        raise ValueError(
            f"{path.name}: coordinate {name} needs at least two values.")
    return float(np.nanmedian(diffs))


def _dataset_crs(ds, da, x_name: str, y_name: str) -> str | None:
    if x_name.lower() in {"lon", "longitude"} and y_name.lower() in {
            "lat", "latitude"}:
        return "EPSG:4326"

    mapping_name = da.attrs.get("grid_mapping")
    candidates = [ds.attrs]
    if mapping_name and mapping_name in ds.variables:
        candidates.append(ds[mapping_name].attrs)
    for name in ("spatial_ref", "crs"):
        if name in ds.variables:
            candidates.append(ds[name].attrs)
    for attrs in candidates:
        for key in ("spatial_ref", "crs_wkt", "crs"):
            value = attrs.get(key)
            if value:
                return str(value)
    return None


def _coordinate_unit(ds, x_name: str, crs: str | None) -> str:
    units = str(ds[x_name].attrs.get("units", "")) if x_name in ds else ""
    if "degree" in units.lower() or x_name.lower() in {"lon", "longitude"}:
        return "deg"
    if crs:
        try:
            from pyproj import CRS

            return _crs_unit(CRS.from_user_input(crs))
        except Exception:
            pass
    return units or ""


def _crs_unit(crs) -> str:
    if getattr(crs, "is_geographic", False):
        return "deg"
    name = crs.axis_info[0].unit_name.lower() if crs.axis_info else ""
    if name in {"metre", "meter"}:
        return "m"
    if "foot" in name or "feet" in name:
        return "ft"
    return name


def _same_crs(left: str, right: str) -> bool:
    try:
        from pyproj import CRS

        return CRS.from_user_input(left) == CRS.from_user_input(right)
    except Exception:
        return str(left).strip().lower() == str(right).strip().lower()


def _validate_spatial_signatures(signatures: list[_SpatialSignature]) -> None:
    import numpy as np

    first = signatures[0]
    for current in signatures[1:]:
        if current.metadata.shape != first.metadata.shape:
            raise ValueError(
                "NetCDF files do not share the same spatial dimensions: "
                f"{first.metadata.files[0].name}={first.metadata.shape}, "
                f"{current.metadata.files[0].name}={current.metadata.shape}.")
        if (
                np.shape(current.x_values) != np.shape(first.x_values)
                or np.shape(current.y_values) != np.shape(first.y_values)
                or not np.allclose(current.x_values, first.x_values,
                                   rtol=1e-9, atol=1e-9, equal_nan=True)
                or not np.allclose(current.y_values, first.y_values,
                                   rtol=1e-9, atol=1e-9, equal_nan=True)):
            raise ValueError(
                "NetCDF files do not share identical spatial coordinates: "
                f"{first.metadata.files[0].name} and "
                f"{current.metadata.files[0].name}.")
        left, right = first.metadata.crs, current.metadata.crs
        if (left is None) != (right is None) or (
                left and right and not _same_crs(left, right)):
            raise ValueError(
                "NetCDF files do not share the same CRS: "
                f"{first.metadata.files[0].name} and "
                f"{current.metadata.files[0].name}.")


def _validate_times(files: list[Path], spec: MeteoFolderSpec) -> None:
    import numpy as np

    groups = [files]
    if spec.source == MHM_READY and spec.kind == "temperature":
        groups = [[path] for path in files]
        reference = _time_values(files[0])
        for path in files[1:]:
            if not np.array_equal(reference, _time_values(path)):
                raise ValueError(
                    "mHM-ready tavg, tmin, and tmax files must share time values.")

    for group in groups:
        values = [
            value
            for path in group
            for value in _time_values(path)
        ]
        if len(np.unique(values)) != len(values):
            raise ValueError(
                f"{spec.kind}: NetCDF files contain duplicate time values.")


def _time_values(path: Path):
    import numpy as np

    with _open_dataset(path) as ds:
        name = "time" if "time" in ds.coords else "valid_time"
        if name not in ds.coords:
            raise ValueError(f"{path.name}: missing time coordinate.")
        return np.asarray(ds[name].values)


def _ready_path(files: tuple[Path, ...], variable: str) -> Path:
    if len(files) == 1:
        return files[0]
    return next(path for path in files if path.name == f"{variable}.nc")


def _validate_target_coverage(
        metadata: SpatialMetadata,
        target: TargetGrid) -> None:
    """Reject a target whose cell centres fall outside the source grid."""
    import numpy as np

    target_x = np.asarray(
        target.sample_lon
        if target.sample_lon is not None else target.lon,
        dtype="float64",
    )
    target_y = np.asarray(
        target.sample_lat
        if target.sample_lat is not None else target.lat,
        dtype="float64",
    )
    if (
            target.crs
            and metadata.crs
            and _same_crs(metadata.crs, target.crs)):
        target_x, target_y = _header_center_axes(target.header, np)
    elif metadata.crs and not _same_crs(metadata.crs, "EPSG:4326"):
        try:
            from pyproj import Transformer
        except Exception as exc:
            raise RuntimeError(
                "pyproj is required to validate projected meteorology data.") from exc
        transform = Transformer.from_crs(
            "EPSG:4326", metadata.crs, always_xy=True)
        target_x, target_y = transform.transform(target_x, target_y)
    elif metadata.crs is None and metadata.unit not in {"deg", "degree"}:
        return

    west, east, south, north = metadata.bounds
    x_tolerance = metadata.x_resolution / 2.0 + 1e-9
    y_tolerance = metadata.y_resolution / 2.0 + 1e-9
    if (
            float(target_x.min()) < west - x_tolerance
            or float(target_x.max()) > east + x_tolerance
            or float(target_y.min()) < south - y_tolerance
            or float(target_y.max()) > north + y_tolerance):
        raise ValueError(
            f"The target grid is outside the {metadata.kind} source extent.")


def _prepare_ready_dataset(
        path: Path,
        variable: str,
        crs_fallback: str | None,
        target: TargetGrid):
    import numpy as np
    import pandas as pd
    import xarray as xr

    spec = MeteoFolderSpec(
        "precipitation" if variable == "pre" else
        "temperature" if variable.startswith("t") else "pet",
        path.parent,
        MHM_READY,
        crs_fallback,
    ).normalized()
    with _open_dataset(path) as ds:
        source_name = _find_variable(ds, spec.kind, variable)
        da = ds[source_name]
        time_name = "time" if "time" in da.dims else "valid_time"
        if time_name not in da.dims:
            raise ValueError(f"{path.name}: variable {source_name} has no time axis.")
        if time_name != "time":
            da = da.rename({time_name: "time"})
        spatial_dims = set(_spatial_dims(da))
        for dim in tuple(da.dims):
            if dim != "time" and dim not in spatial_dims and da.sizes[dim] == 1:
                da = da.isel({dim: 0}, drop=True)
        da = _resample_ready(da, ds, spec.crs, target, np, xr)
        da = da.astype("float64").load().sortby("time")

    da.name = variable
    attrs = {
        "pre": ("mm", "daily total precipitation"),
        "tavg": ("Celsius", "daily mean air temperature"),
        "tmin": ("Celsius", "daily minimum air temperature"),
        "tmax": ("Celsius", "daily maximum air temperature"),
        "pet": ("mm", "daily potential evapotranspiration"),
    }
    da.attrs = {"units": attrs[variable][0], "long_name": attrs[variable][1]}
    ds_out = da.to_dataset()
    starts = pd.DatetimeIndex(ds_out["time"].values)
    ds_out["time_bnds"] = xr.DataArray(
        np.stack([
            starts.values.astype("datetime64[ns]"),
            (starts + pd.Timedelta(days=1)).values.astype("datetime64[ns]"),
        ], axis=1),
        dims=("time", "bnds"),
        coords={"time": ds_out["time"], "bnds": [0, 1]},
    )
    ds_out["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "bounds": "time_bnds"})
    ds_out["lat"].attrs.update(
        {"standard_name": "latitude", "units": "degrees_north", "axis": "Y"})
    ds_out["lon"].attrs.update(
        {"standard_name": "longitude", "units": "degrees_east", "axis": "X"})
    ds_out.attrs.update({
        "Conventions": "CF-1.6",
        "source": "mHM ready",
        "history": "Prepared for mHM by pymhm",
    })
    return ds_out


def _resample_ready(da, ds, source_crs, target, np, xr):
    y_dim, x_dim = _spatial_dims(da)
    x_name, x_values = _axis_values(ds, x_dim, "x")
    y_name, y_values = _axis_values(ds, y_dim, "y")
    target_lon = np.asarray(target.lon, dtype="float64")
    target_lat = np.asarray(target.lat, dtype="float64")
    sample_lon = (
        np.asarray(target.sample_lon, dtype="float64")
        if target.sample_lon is not None else None
    )
    sample_lat = (
        np.asarray(target.sample_lat, dtype="float64")
        if target.sample_lat is not None else None
    )

    if x_values is not None and y_values is not None:
        source_crs = _dataset_crs(ds, da, x_name, y_name) or source_crs
        lon_name, lon = _geographic_coordinate(ds, "lon", (y_dim, x_dim))
        lat_name, lat = _geographic_coordinate(ds, "lat", (y_dim, x_dim))
        if source_crs is None and lon is not None and lat is not None:
            return _nearest_curvilinear(
                da, np.asarray(lon), np.asarray(lat),
                target_lon, target_lat, xr, sample_lon, sample_lat)
        x_target, y_target = target_lon, target_lat
        if (
                source_crs
                and target.crs
                and _same_crs(source_crs, target.crs)):
            x_target, y_target = _header_center_axes(target.header, np)
        elif sample_lon is not None:
            x_target, y_target = sample_lon, sample_lat
            if source_crs and not _same_crs(source_crs, "EPSG:4326"):
                try:
                    from pyproj import Transformer
                except Exception as exc:
                    raise RuntimeError(
                        "pyproj is required to resample projected "
                        "mHM-ready data."
                    ) from exc
                transform = Transformer.from_crs(
                    "EPSG:4326",
                    source_crs,
                    always_xy=True,
                )
                x_target, y_target = transform.transform(
                    x_target,
                    y_target,
                )
            return _nearest_rectilinear(
                da,
                np.asarray(x_values),
                np.asarray(y_values),
                x_target,
                y_target,
                target_lon,
                target_lat,
                xr,
            )
        elif source_crs and not _same_crs(source_crs, "EPSG:4326"):
            try:
                from pyproj import Transformer
            except Exception as exc:
                raise RuntimeError(
                    "pyproj is required to resample projected mHM-ready data.") from exc
            transform = Transformer.from_crs(
                "EPSG:4326", source_crs, always_xy=True)
            center_lat = float((target_lat.min() + target_lat.max()) / 2.0)
            center_lon = float((target_lon.min() + target_lon.max()) / 2.0)
            x_target = np.asarray([
                transform.transform(lon, center_lat)[0]
                for lon in target_lon])
            y_target = np.asarray([
                transform.transform(center_lon, lat)[1]
                for lat in target_lat])

        source = da.assign_coords(
            {x_dim: np.asarray(x_values), y_dim: np.asarray(y_values)})
        source = source.sortby(x_dim).sortby(y_dim)
        result = source.interp(
            {x_dim: x_target, y_dim: y_target}, method="nearest")
        result = result.rename({x_dim: "lon", y_dim: "lat"})
        return result.assign_coords(lon=target_lon, lat=target_lat).transpose(
            "time", "lat", "lon")

    lon_name, lon = _geographic_coordinate(ds, "lon", (y_dim, x_dim))
    lat_name, lat = _geographic_coordinate(ds, "lat", (y_dim, x_dim))
    if lon is None or lat is None:
        raise ValueError("mHM-ready data has no usable spatial coordinates.")
    return _nearest_curvilinear(
        da,
        np.asarray(lon),
        np.asarray(lat),
        target_lon,
        target_lat,
        xr,
        sample_lon,
        sample_lat,
    )


def _header_center_axes(header, np):
    cellsize = float(header["cellsize"])
    xll = float(header["xllcorner"])
    yll = float(header["yllcorner"])
    x_values = xll + (np.arange(int(header["ncols"])) + 0.5) * cellsize
    y_values = yll + (np.arange(int(header["nrows"])) + 0.5) * cellsize
    return x_values, y_values[::-1]


def _nearest_rectilinear(
        da,
        source_x,
        source_y,
        target_x,
        target_y,
        output_lon,
        output_lat,
        xr):
    """Sample a rectilinear source at exact 2-D target points."""
    import numpy as np

    source_x_mesh, source_y_mesh = np.meshgrid(source_x, source_y)
    return _nearest_curvilinear(
        da,
        source_x_mesh,
        source_y_mesh,
        output_lon,
        output_lat,
        xr,
        target_x,
        target_y,
    )


def _nearest_curvilinear(
        da,
        source_lon,
        source_lat,
        target_lon,
        target_lat,
        xr,
        sample_lon=None,
        sample_lat=None):
    import numpy as np

    try:
        from scipy.spatial import cKDTree
    except Exception as exc:
        raise RuntimeError(
            "scipy is required to resample curvilinear mHM-ready data.") from exc

    y_dim, x_dim = _spatial_dims(da)
    check_lon = sample_lon if sample_lon is not None else target_lon
    check_lat = sample_lat if sample_lat is not None else target_lat
    if (
            np.min(check_lon) < np.nanmin(source_lon)
            or np.max(check_lon) > np.nanmax(source_lon)
            or np.min(check_lat) < np.nanmin(source_lat)
            or np.max(check_lat) > np.nanmax(source_lat)):
        raise ValueError("The target grid is outside the source extent.")
    source_points = np.column_stack(
        [source_lon.reshape(-1), source_lat.reshape(-1)])
    if sample_lon is None or sample_lat is None:
        lon_mesh, lat_mesh = np.meshgrid(target_lon, target_lat)
    else:
        lon_mesh = np.asarray(sample_lon, dtype="float64")
        lat_mesh = np.asarray(sample_lat, dtype="float64")
    target_points = np.column_stack(
        [lon_mesh.reshape(-1), lat_mesh.reshape(-1)])
    _, indices = cKDTree(source_points).query(target_points)
    values = np.asarray(da.transpose("time", y_dim, x_dim).values)
    sampled = values.reshape(values.shape[0], -1)[:, indices]
    sampled = sampled.reshape(values.shape[0], len(target_lat), len(target_lon))
    return xr.DataArray(
        sampled,
        dims=("time", "lat", "lon"),
        coords={"time": da["time"], "lat": target_lat, "lon": target_lon},
    )


def _attach_target_coordinates(ds, variable: str, target: TargetGrid):
    """Attach exact 2-D cell-centre coordinates for projected target grids."""
    if target.sample_lon is None or target.sample_lat is None:
        return ds

    ds = ds.assign_coords({
        "lon2d": (
            ("lat", "lon"),
            target.sample_lon,
            {
                "standard_name": "longitude",
                "units": "degrees_east",
            },
        ),
        "lat2d": (
            ("lat", "lon"),
            target.sample_lat,
            {
                "standard_name": "latitude",
                "units": "degrees_north",
            },
        ),
    })
    ds[variable].attrs["coordinates"] = "lat2d lon2d"
    if target.crs:
        try:
            from pyproj import CRS

            crs = CRS.from_user_input(target.crs)
            crs_attrs = crs.to_cf()
            crs_attrs["spatial_ref"] = crs.to_wkt()
            ds["crs"] = 0
            ds["crs"].attrs.update(crs_attrs)
            ds[variable].attrs["grid_mapping"] = "crs"
        except Exception:
            pass
    return ds


def _normalize_precipitation(ds):
    """Constrain precipitation to its physical lower bound."""
    attrs = dict(ds["pre"].attrs)
    ds = ds.copy()
    ds["pre"] = ds["pre"].clip(min=0.0)
    ds["pre"].attrs = {**attrs, "valid_min": 0.0}
    history = ds.attrs.get("history", "")
    note = "Precipitation values below 0 mm clipped to 0."
    ds.attrs["history"] = f"{history}; {note}" if history else note
    return ds


__all__ = [
    "ERA5LAND",
    "MHM_READY",
    "MeteoFolderSpec",
    "MeteoForcingResult",
    "SpatialMetadata",
    "SpatialResolution",
    "TargetGrid",
    "inspect_meteo_folder",
    "inspect_meteo_inputs",
    "normalize_kind",
    "normalize_source",
    "process_meteo_inputs",
    "resolution_in_crs",
]
