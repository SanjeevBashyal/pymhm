# -*- coding: utf-8 -*-
"""mHM morphology ASCII preparation through mhm_tools file handlers."""
from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .....grid_resolution import (
    axes_match_header,
    ceil_cellsize,
    cellsize_precision_for_unit,
    reindex_to_header,
)


Header = dict[str, Any]
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
    unit = _headers_unit(standardized)

    ratios = {
        "L0_to_L1": _resolution_ratio(
            standardized["L0"]["cellsize"],
            standardized["L1"]["cellsize"],
            "L0",
            "L1",
            unit,
        ),
        "L1_to_L2": _resolution_ratio(
            standardized["L1"]["cellsize"],
            standardized["L2"]["cellsize"],
            "L1",
            "L2",
            unit,
        ),
        "L11_to_L2": _resolution_ratio(
            standardized["L11"]["cellsize"],
            standardized["L2"]["cellsize"],
            "L11",
            "L2",
            unit,
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
    import xarray as xr
    from .....applications.mhm_tools_handler import (xarray_coord_key,
                                                   xarray_single_data_var)

    header = _normalise_header(header)
    var_name = data_var
    if var_name is None:
        var_name = xarray_single_data_var(dataset)
    if var_name is None or var_name not in dataset:
        raise ValueError("Cannot determine a single raster data variable for ASCII export.")

    x_key = xarray_coord_key(dataset, lon=True)
    y_key = xarray_coord_key(dataset, lat=True)
    source = dataset.sortby(x_key).sortby(y_key)
    _assert_source_covers_header(source, header, x_key, y_key)

    target_x, target_y = _target_coordinates(header)
    da = source[var_name]
    da = _drop_singleton_non_spatial_dims(da, x_key, y_key)
    # Sample by nearest label, not by interpolation. The source y-centres come
    # from the GDAL geotransform (top-down) while _target_coordinates builds
    # them bottom-up: algebraically equal, ~3.6e-15 apart in float64. interp()
    # has no extrapolation domain, so that one ulp put the bottom row outside it
    # and _apply_nodata turned the resulting NaN into nodata -- an all-nodata
    # last row in every written layer. interp() also framed a genuinely coarser
    # source in nodata. sel() snaps to the nearest label and has no domain.
    # float64 is the input contract _apply_nodata already gets from
    # pad_dataset_to_header; interp() used to provide it by upcasting.
    aligned = da.sel(
        {x_key: target_x, y_key: target_y},
        method="nearest",
    ).astype("float64")
    # sel() returns the source labels, so this restores the exact header
    # coordinates that the written ASCII header is derived from. Load-bearing.
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


def pad_dataset_to_header(
        dataset,
        header: Mapping[str, Any],
        data_var: str | None = None,
        nodata_value: float | int = -9999,
        integer: bool = False,
        pad_value: float | int | None = None):
    """Return a dataset placed on ``header`` by exact cell lookup, padded out.

    Unlike :func:`align_dataset_to_header` this never resamples. The source must
    already sit on the header cell grid; cells outside it receive ``pad_value``,
    or ``nodata_value`` when no pad value is given. Non-spatial dimensions such
    as time or soil horizon are preserved.
    """
    from .....applications.mhm_tools_handler import (xarray_coord_key,
                                                   xarray_single_data_var)

    header = _standardize_header(header)
    var_name = data_var
    if var_name is None:
        var_name = xarray_single_data_var(dataset)
    if var_name is None or var_name not in dataset:
        raise ValueError("Cannot determine a single raster data variable for L0 padding.")

    x_key = xarray_coord_key(dataset, lon=True)
    y_key = xarray_coord_key(dataset, lat=True)
    source = dataset.sortby(x_key).sortby(y_key)
    _assert_cells_align(source, x_key, y_key, header)

    padded = reindex_to_header(
        source[var_name].astype("float64"),
        x_key,
        y_key,
        header,
        fill_value=float(nodata_value if pad_value is None else pad_value),
    )
    padded = _apply_nodata(padded, nodata_value, integer)
    if (
            padded.sizes.get(y_key) != int(header["nrows"])
            or padded.sizes.get(x_key) != int(header["ncols"])):
        raise ValueError(
            "Padded raster size does not match L0 header: "
            f"got {padded.sizes.get(y_key)}x{padded.sizes.get(x_key)}, "
            f"expected {int(header['nrows'])}x{int(header['ncols'])}."
        )
    return padded.to_dataset(name=var_name)


def pad_l0_file_to_header(
        path: Path | str,
        header: Mapping[str, Any],
        *,
        data_var: str | None = None,
        nodata_value: float | int = -9999,
        integer: bool = False,
        pad_value: float | int | None = None) -> Path:
    """Rewrite a model-ready raster onto the exact L0 header without resampling.

    Cells the source does not cover take ``pad_value``; the declared nodata
    stays ``nodata_value`` regardless.
    """
    path = Path(path)
    header = _standardize_header(header)
    if path.suffix.lower() == ".asc":
        # ASCII grids are padded by streaming text, which keeps memory flat.
        # Loading one costs gigabytes: 3.2 GiB for a 13201 x 6001 soil raster.
        from .pad import pad_ascii_grid

        return pad_ascii_grid(path, header, nodata=nodata_value, pad=pad_value)
    dataset = _read_raster(path, data_var)
    if _grid_matches_header(dataset, header):
        return path

    padded = pad_dataset_to_header(
        dataset,
        header,
        data_var=data_var,
        nodata_value=nodata_value,
        integer=integer,
        pad_value=pad_value,
    )
    temporary = path.with_name(f".{path.name}.tmp{path.suffix or '.nc'}")
    temporary.unlink(missing_ok=True)
    padded.to_netcdf(temporary)
    dataset.close()
    temporary.replace(path)
    return path


def _grid_matches_header(dataset, header: Mapping[str, Any]) -> bool:
    """Return True when a dataset already sits exactly on the header grid."""
    import numpy as np
    from .....applications.mhm_tools_handler import xarray_coord_key

    x_key = xarray_coord_key(dataset, lon=True)
    y_key = xarray_coord_key(dataset, lat=True)
    target_x, target_y = _target_coordinates(header)
    if (
            dataset.sizes.get(x_key) != len(target_x)
            or dataset.sizes.get(y_key) != len(target_y)):
        return False
    tolerance = float(header["cellsize"]) * 1e-6
    return bool(
        np.allclose(np.sort(dataset[x_key].values), target_x, rtol=0.0, atol=tolerance)
        and np.allclose(
            np.sort(dataset[y_key].values), np.sort(target_y),
            rtol=0.0, atol=tolerance)
    )


def _assert_cells_align(dataset, x_key, y_key, header) -> None:
    """Reject a source whose cells are not on the header cell grid."""
    if not axes_match_header(
            dataset[x_key].values, dataset[y_key].values, dict(header)):
        raise ValueError(
            "Source raster cells are not aligned to the L0 grid "
            f"(cell size {float(header['cellsize'])}); it cannot be placed on "
            "the common extent without resampling."
        )


def _read_raster(path: Path, data_var: str | None):
    from .....applications.mhm_tools_handler import read_xarray_dataset

    return read_xarray_dataset(
        path,
        var_name=data_var,
        force_decending_y=True,
    )


def _write_ascii(dataset, output_path: Path, header: Mapping[str, Any],
                 nodata_value: float | int) -> None:
    from .....applications.mhm_tools_handler import (write_xarray_ascii,
                                                   xarray_single_data_var)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    data_var = xarray_single_data_var(dataset)
    write_xarray_ascii(
        dataset,
        output_path,
        data_var=data_var,
        nodata_value=nodata_value,
        resolution=float(header["cellsize"]),
    )
    expected = _normalise_header(header)
    written = _read_ascii_header(output_path, unit=expected.get("unit", ""))
    _assert_header_matches(written, expected, output_path.name)


def _normalise_header(header: Mapping[str, Any]) -> Header:
    normalized = {str(key).lower(): value for key, value in dict(header).items()}
    if "nodata_value" not in normalized and "nodata" in normalized:
        normalized["nodata_value"] = normalized["nodata"]
    if "nodata_value" not in normalized:
        normalized["nodata_value"] = -9999.0
    if "unit" not in normalized and "crs_unit" in normalized:
        normalized["unit"] = normalized["crs_unit"]
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
    unit = str(normalized.get("unit", "") or "")
    # The cell size is grid-defining and is kept exactly as supplied; rounding
    # it here turns exact level ratios such as 120 into 119.99999.
    cellsize = float(normalized["cellsize"])
    return {
        "ncols": int(float(normalized["ncols"])),
        "nrows": int(float(normalized["nrows"])),
        "xllcorner": float(normalized["xllcorner"]),
        "yllcorner": float(normalized["yllcorner"]),
        "cellsize": cellsize,
        "nodata_value": float(normalized.get("nodata_value", -9999.0)),
        "unit": unit,
        "cellsize_precision": cellsize_precision_for_unit(unit),
    }


def _headers_unit(headers: Mapping[str, Header]) -> str:
    """Return the first non-empty grid unit recorded in the headers."""
    for header in headers.values():
        unit = str(header.get("unit", "") or "")
        if unit:
            return unit
    return ""


def _resolution_ratio(
        fine_cellsize: float,
        coarse_cellsize: float,
        fine_name: str,
        coarse_name: str,
        unit: str | None = None,
        tolerance: float | None = None) -> int:
    # Compare the cell sizes as given. Rounding them to a fixed number of
    # decimals turns an exact ratio such as 0.1 / (1/1200) into 119.99999.
    fine_cellsize = float(fine_cellsize)
    coarse_cellsize = float(coarse_cellsize)
    precision = cellsize_precision_for_unit(unit)
    if tolerance is None:
        tolerance = max(10 ** -precision, 1e-9)
    if fine_cellsize <= 0 or coarse_cellsize <= 0:
        raise ValueError("Grid cell sizes must be positive.")
    if fine_cellsize - coarse_cellsize > tolerance:
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


def _read_ascii_header(path: Path, unit: str | None = None) -> Header:
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
    if unit:
        header["unit"] = unit
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
    "pad_dataset_to_header",
    "pad_l0_file_to_header",
    "prepare_morphology_ascii_files",
    "validate_grid_headers",
]
