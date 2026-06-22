# -*- coding: utf-8 -*-
"""mHM morphology ASCII preparation through mhm_tools file handlers."""
from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .._bundled import ensure_bundled_mhm_tools


Header = dict[str, float | int]
LogCallback = Callable[[str], None]


@dataclass(frozen=True)
class MorphologyAsciiLayer:
    """One raster-to-ASCII export request."""

    input_path: Path | str
    output_path: Path | str
    name: str = ""
    nodata_value: float | int = -9999
    integer: bool = False
    data_var: str | None = None


@dataclass(frozen=True)
class MorphologyAsciiResult:
    """Prepared morphology ASCII outputs and grid compatibility details."""

    outputs: dict[str, Path]
    headers: dict[str, Header]
    ratios: dict[str, int]


def prepare_morphology_ascii_files(
        layers: Iterable[MorphologyAsciiLayer | Mapping[str, Any]],
        headers: Mapping[str, Mapping[str, Any]],
        overwrite: bool = True,
        log: LogCallback | None = None) -> MorphologyAsciiResult:
    """Crop/align prepared morphology rasters to L0 and write mHM ASCII files."""
    standardized_headers, ratios = validate_grid_headers(headers)
    l0_header = standardized_headers["L0"]
    outputs: dict[str, Path] = {}

    _log(log, "Grid compatibility verified for ASCII export.")
    _log(
        log,
        "Grid ratios: "
        f"L0/L1={ratios['L0_to_L1']}, "
        f"L1/L2={ratios['L1_to_L2']}, "
        f"L11/L2={ratios['L11_to_L2']}.",
    )

    for layer in layers:
        spec = _layer_spec(layer)
        input_path = Path(spec.input_path)
        output_path = Path(spec.output_path)
        display_name = spec.name or output_path.name
        if not input_path.exists():
            raise FileNotFoundError(f"Input raster does not exist: {input_path}")
        if output_path.exists() and not overwrite:
            _log(log, f"{display_name} already exists. Skipping: {output_path}")
            outputs[display_name] = output_path
            continue

        _log(log, f"Preparing {display_name}: {input_path.name} -> {output_path.name}")
        dataset = _read_raster(input_path, spec.data_var)
        aligned = align_dataset_to_header(
            dataset,
            l0_header,
            data_var=spec.data_var,
            nodata_value=spec.nodata_value,
            integer=spec.integer,
        )
        _write_ascii(
            aligned,
            output_path,
            l0_header,
            nodata_value=spec.nodata_value,
        )
        outputs[display_name] = output_path
        _log(log, f"Written {display_name}: {output_path}")

    return MorphologyAsciiResult(
        outputs=outputs,
        headers=standardized_headers,
        ratios=ratios,
    )


def validate_grid_headers(
        headers: Mapping[str, Mapping[str, Any]]) -> tuple[dict[str, Header], dict[str, int]]:
    """Validate L0, L1, L11, and L2 headers for mHM-compatible dimensions."""
    required = ("L0", "L1", "L11", "L2")
    missing = [level for level in required if level not in headers]
    if missing:
        raise ValueError(f"Missing grid header(s): {', '.join(missing)}")

    standardized = {
        level: _standardize_header(headers[level])
        for level in required
    }
    _assert_same_extent(standardized)

    ratios = {
        "L0_to_L1": _resolution_ratio(
            standardized["L0"]["cellsize"],
            standardized["L1"]["cellsize"],
            "L0",
            "L1",
        ),
        "L1_to_L2": _resolution_ratio(
            standardized["L1"]["cellsize"],
            standardized["L2"]["cellsize"],
            "L1",
            "L2",
        ),
        "L11_to_L2": _resolution_ratio(
            standardized["L11"]["cellsize"],
            standardized["L2"]["cellsize"],
            "L11",
            "L2",
        ),
    }
    _assert_matrix_multiple(
        standardized["L0"],
        standardized["L1"],
        ratios["L0_to_L1"],
        "L0",
        "L1",
    )
    _assert_matrix_multiple(
        standardized["L1"],
        standardized["L2"],
        ratios["L1_to_L2"],
        "L1",
        "L2",
    )
    _assert_matrix_multiple(
        standardized["L11"],
        standardized["L2"],
        ratios["L11_to_L2"],
        "L11",
        "L2",
    )
    return standardized, ratios


def align_dataset_to_header(
        dataset,
        header: Mapping[str, Any],
        data_var: str | None = None,
        nodata_value: float | int = -9999,
        integer: bool = False):
    """Return a dataset sampled on the exact grid described by ``header``."""
    import numpy as np
    import xarray as xr
    ensure_bundled_mhm_tools()
    from mhm_tools.common.xarray_utils import get_coord_key, get_single_data_var

    header = _normalise_header(header)
    var_name = data_var
    if var_name is None:
        var_name = get_single_data_var(dataset)
    if var_name is None or var_name not in dataset:
        raise ValueError("Cannot determine a single raster data variable for ASCII export.")

    x_key = get_coord_key(dataset, lon=True)
    y_key = get_coord_key(dataset, lat=True)
    source = dataset.sortby(x_key).sortby(y_key)
    _assert_source_covers_header(source, header, x_key, y_key)

    target_x, target_y = _target_coordinates(header)
    da = source[var_name]
    da = _drop_singleton_non_spatial_dims(da, x_key, y_key)
    aligned = da.interp(
        {x_key: target_x, y_key: target_y},
        method="nearest",
    )
    aligned = aligned.assign_coords({
        x_key: xr.DataArray(
            target_x,
            dims=(x_key,),
            attrs={"axis": "X"},
        ),
        y_key: xr.DataArray(
            target_y,
            dims=(y_key,),
            attrs={"axis": "Y"},
        ),
    })
    aligned = _apply_nodata(aligned, nodata_value, integer)
    if aligned.sizes.get(y_key) != int(header["nrows"]) or aligned.sizes.get(x_key) != int(header["ncols"]):
        raise ValueError(
            "Aligned raster size does not match L0 header: "
            f"got {aligned.sizes.get(y_key)}x{aligned.sizes.get(x_key)}, "
            f"expected {int(header['nrows'])}x{int(header['ncols'])}."
        )
    return aligned.to_dataset(name=var_name)


def _read_raster(path: Path, data_var: str | None):
    ensure_bundled_mhm_tools()
    from mhm_tools.common.file_handler import get_xarray_ds_from_file

    return get_xarray_ds_from_file(
        path,
        var_name=data_var,
        force_decending_y=True,
    )


def _write_ascii(dataset, output_path: Path, header: Mapping[str, Any],
                 nodata_value: float | int) -> None:
    ensure_bundled_mhm_tools()
    from mhm_tools.common.file_handler import write_xarray_to_ascii
    from mhm_tools.common.xarray_utils import get_single_data_var

    output_path.parent.mkdir(parents=True, exist_ok=True)
    data_var = get_single_data_var(dataset)
    write_xarray_to_ascii(
        dataset,
        output_path,
        data_var=data_var,
        nodata_value=nodata_value,
        resolution=float(header["cellsize"]),
    )
    written = _read_ascii_header(output_path)
    _assert_header_matches(written, _normalise_header(header), output_path.name)


def _normalise_header(header: Mapping[str, Any]) -> Header:
    normalized = {str(key).lower(): value for key, value in dict(header).items()}
    if "nodata_value" not in normalized and "nodata" in normalized:
        normalized["nodata_value"] = normalized["nodata"]
    if "nodata_value" not in normalized:
        normalized["nodata_value"] = -9999.0
    return normalized


def _standardize_header(header: Mapping[str, Any]) -> Header:
    normalized = _normalise_header(header)
    if "xllcenter" in normalized and "xllcorner" not in normalized:
        normalized["xllcorner"] = (
            float(normalized["xllcenter"]) - 0.5 * float(normalized.get("cellsize", 1.0))
        )
    if "yllcenter" in normalized and "yllcorner" not in normalized:
        normalized["yllcorner"] = (
            float(normalized["yllcenter"]) - 0.5 * float(normalized.get("cellsize", 1.0))
        )
    required = {"ncols", "nrows", "xllcorner", "yllcorner", "cellsize"}
    missing = required - set(normalized)
    if missing:
        raise ValueError(f"Grid header is missing required key(s): {', '.join(sorted(missing))}")
    return {
        "ncols": int(float(normalized["ncols"])),
        "nrows": int(float(normalized["nrows"])),
        "xllcorner": float(normalized["xllcorner"]),
        "yllcorner": float(normalized["yllcorner"]),
        "cellsize": float(normalized["cellsize"]),
        "nodata_value": float(normalized.get("nodata_value", -9999.0)),
    }


def _resolution_ratio(
        fine_cellsize: float,
        coarse_cellsize: float,
        fine_name: str,
        coarse_name: str,
        tolerance: float = 1e-7) -> int:
    fine_cellsize = float(fine_cellsize)
    coarse_cellsize = float(coarse_cellsize)
    if fine_cellsize <= 0 or coarse_cellsize <= 0:
        raise ValueError("Grid cell sizes must be positive.")
    if fine_cellsize > coarse_cellsize:
        raise ValueError(
            f"{fine_name} cellsize ({fine_cellsize}) must be finer than or equal to "
            f"{coarse_name} cellsize ({coarse_cellsize})."
        )
    ratio = coarse_cellsize / fine_cellsize
    rounded = int(round(ratio))
    if abs(ratio - rounded) > tolerance:
        raise ValueError(
            f"{coarse_name}/{fine_name} resolution ratio must be an integer. "
            f"Got {ratio}."
        )
    return max(1, rounded)


def _layer_spec(layer: MorphologyAsciiLayer | Mapping[str, Any]) -> MorphologyAsciiLayer:
    if isinstance(layer, MorphologyAsciiLayer):
        return layer
    return MorphologyAsciiLayer(
        input_path=layer["input_path"],
        output_path=layer["output_path"],
        name=layer.get("name", ""),
        nodata_value=layer.get("nodata_value", -9999),
        integer=bool(layer.get("integer", False)),
        data_var=layer.get("data_var"),
    )


def _target_coordinates(header: Mapping[str, Any]):
    import numpy as np

    cellsize = float(header["cellsize"])
    xll = float(header["xllcorner"])
    yll = float(header["yllcorner"])
    ncols = int(header["ncols"])
    nrows = int(header["nrows"])
    x_values = xll + (np.arange(ncols, dtype="float64") + 0.5) * cellsize
    y_values = yll + (np.arange(nrows, dtype="float64") + 0.5) * cellsize
    return x_values, y_values[::-1]


def _header_bounds(header: Mapping[str, Any]) -> tuple[float, float, float, float]:
    cellsize = float(header["cellsize"])
    xmin = float(header["xllcorner"])
    ymin = float(header["yllcorner"])
    xmax = xmin + int(header["ncols"]) * cellsize
    ymax = ymin + int(header["nrows"]) * cellsize
    return xmin, xmax, ymin, ymax


def _assert_same_extent(headers: Mapping[str, Header]) -> None:
    reference = _header_bounds(headers["L0"])
    tolerance = max(float(h["cellsize"]) for h in headers.values()) * 1e-6
    for level, header in headers.items():
        bounds = _header_bounds(header)
        if any(abs(left - right) > tolerance for left, right in zip(reference, bounds)):
            raise ValueError(
                f"{level} extent {bounds} does not match L0 extent {reference}."
            )


def _assert_matrix_multiple(
        fine: Header,
        coarse: Header,
        ratio: int,
        fine_name: str,
        coarse_name: str) -> None:
    expected_cols = int(coarse["ncols"]) * ratio
    expected_rows = int(coarse["nrows"]) * ratio
    if int(fine["ncols"]) != expected_cols or int(fine["nrows"]) != expected_rows:
        raise ValueError(
            f"{fine_name} matrix must be {ratio} times {coarse_name}. "
            f"Got {fine_name}=({fine['nrows']} rows, {fine['ncols']} cols), "
            f"expected ({expected_rows} rows, {expected_cols} cols) from {coarse_name}."
        )


def _assert_header_matches(written: Mapping[str, Any], expected: Header, name: str) -> None:
    tolerance = float(expected["cellsize"]) * 1e-6
    for key in ("ncols", "nrows"):
        if int(written[key]) != int(expected[key]):
            raise ValueError(
                f"Written ASCII {name} has {key}={written[key]}, expected {expected[key]}."
            )
    for key in ("xllcorner", "yllcorner", "cellsize"):
        if abs(float(written[key]) - float(expected[key])) > tolerance:
            raise ValueError(
                f"Written ASCII {name} has {key}={written[key]}, expected {expected[key]}."
            )


def _read_ascii_header(path: Path) -> Header:
    header: dict[str, float | int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            key = parts[0].lower()
            if key == "nodata_value":
                header[key] = float(parts[1])
            elif key in {"ncols", "nrows"}:
                header[key] = int(float(parts[1]))
            elif key in {"xllcorner", "yllcorner", "xllcenter", "yllcenter", "cellsize"}:
                header[key] = float(parts[1])
            if len(header) >= 6:
                break
    return _standardize_header(header)


def _assert_source_covers_header(dataset, header: Header, x_key: str, y_key: str) -> None:
    import numpy as np

    x_values = np.asarray(dataset[x_key].values, dtype="float64")
    y_values = np.asarray(dataset[y_key].values, dtype="float64")
    x_res = _coord_resolution(x_values, float(header["cellsize"]))
    y_res = _coord_resolution(y_values, float(header["cellsize"]))
    source_bounds = (
        float(np.nanmin(x_values)) - 0.5 * x_res,
        float(np.nanmax(x_values)) + 0.5 * x_res,
        float(np.nanmin(y_values)) - 0.5 * y_res,
        float(np.nanmax(y_values)) + 0.5 * y_res,
    )
    target_bounds = _header_bounds(header)
    tolerance = max(float(header["cellsize"]), x_res, y_res) * 1e-4
    outside = (
        target_bounds[0] < source_bounds[0] - tolerance
        or target_bounds[1] > source_bounds[1] + tolerance
        or target_bounds[2] < source_bounds[2] - tolerance
        or target_bounds[3] > source_bounds[3] + tolerance
    )
    if outside:
        raise ValueError(
            "Input raster extent does not cover the target L0 extent. "
            f"Source={source_bounds}, target={target_bounds}."
        )


def _coord_resolution(values, fallback: float) -> float:
    import numpy as np

    unique_values = np.unique(values[np.isfinite(values)])
    if unique_values.size < 2:
        return fallback
    return float(np.nanmedian(np.abs(np.diff(np.sort(unique_values)))))


def _drop_singleton_non_spatial_dims(da, x_key: str, y_key: str):
    extra_dims = [dim for dim in da.dims if dim not in {x_key, y_key}]
    for dim in extra_dims:
        if da.sizes.get(dim) != 1:
            raise ValueError(
                f"ASCII export requires a 2D raster; dimension {dim!r} has "
                f"size {da.sizes.get(dim)}."
            )
        da = da.isel({dim: 0}, drop=True)
    return da


def _apply_nodata(da, nodata_value: float | int, integer: bool):
    source_nodata = da.attrs.get("_FillValue", da.attrs.get("nodata_value"))
    if source_nodata is not None:
        try:
            da = da.where(da != source_nodata, nodata_value)
        except Exception:
            pass
    da = da.fillna(nodata_value)
    if integer:
        da = da.round().astype("int32")
    da.attrs["nodata_value"] = nodata_value
    da.attrs["_FillValue"] = nodata_value
    return da


def _log(log: LogCallback | None, message: str) -> None:
    if log:
        log(message)


__all__ = [
    "MorphologyAsciiLayer",
    "MorphologyAsciiResult",
    "align_dataset_to_header",
    "prepare_morphology_ascii_files",
    "validate_grid_headers",
]
