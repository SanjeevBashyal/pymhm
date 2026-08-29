# -*- coding: utf-8 -*-
"""Small path-based adapters for the installed :mod:`mhm_tools` package."""
from __future__ import annotations

import contextlib
import io
import logging
import os
from collections.abc import Callable, Generator
from pathlib import Path
from tempfile import TemporaryDirectory


LogCallback = Callable[[str], None]

_OUTPUTS = {
    "soil": ("format_soil_data", "soil_class.tif", "soil_classdefinition.txt"),
    "geology": (
        "format_geology_data",
        "geology_class.tif",
        "geology_classdefinition.txt",
    ),
    "lc": ("format_lc_data", "lc.tif", None),
}


class _StreamToCallback(io.TextIOBase):
    """Forward completed stream lines to a callback."""

    def __init__(self, callback: LogCallback, prefix: str = "") -> None:
        self._callback = callback
        self._prefix = prefix
        self._buffer = ""

    def writable(self) -> bool:
        return True

    def write(self, text: str) -> int:
        if not text:
            return 0
        self._buffer += str(text)
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            self._emit(line)
        return len(text)

    def flush(self) -> None:
        if self._buffer:
            self._emit(self._buffer)
            self._buffer = ""

    def _emit(self, line: str) -> None:
        line = line.rstrip("\r")
        if line:
            self._callback(f"{self._prefix}{line}")


class _CallbackLoggingHandler(logging.Handler):
    """Forward logging records to a callback."""

    def __init__(self, callback: LogCallback) -> None:
        super().__init__()
        self._callback = callback
        self.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self._callback(self.format(record))
        except Exception:
            self.handleError(record)


@contextlib.contextmanager
def capture_messages(
    callback: LogCallback | None,
    logger_names: tuple[str, ...] = ("mhm_tools",),
) -> Generator[None, None, None]:
    """Capture tool stdout, stderr, and logging into ``callback``."""
    if callback is None:
        yield
        return

    stdout = _StreamToCallback(callback)
    stderr = _StreamToCallback(callback, prefix="ERROR: ")
    handlers: list[tuple[logging.Logger, logging.Handler, int, bool]] = []
    try:
        for logger_name in logger_names:
            logger = logging.getLogger(logger_name)
            handler = _CallbackLoggingHandler(callback)
            handlers.append((logger, handler, logger.level, logger.propagate))
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
            logger.propagate = False

        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            yield
    finally:
        stdout.flush()
        stderr.flush()
        for logger, handler, old_level, old_propagate in handlers:
            logger.removeHandler(handler)
            logger.setLevel(old_level)
            logger.propagate = old_propagate


def read_categorical_lookup_table(lookup_table):
    """Read CSV lookups with mHM-tools and delimited TXT lookups with pandas."""
    lookup_path = Path(lookup_table)
    if lookup_path.suffix.lower() != ".txt":
        from mhm_tools.common.lookup_handler import read_lookup_table

        return read_lookup_table(lookup_path)

    import pandas as pd

    try:
        table = pd.read_csv(
            lookup_path,
            sep=None,
            engine="python",
            encoding="utf-8-sig",
        )
    except Exception as error:
        raise ValueError(
            f"Could not read lookup table {lookup_path}: {error}"
        ) from error
    if table.empty:
        raise ValueError(f"Lookup table {lookup_path} is empty.")
    return table


def _category_key(value):
    """Normalize a vector/lookup category like mHM-tools rasterization."""
    import pandas as pd

    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return text
    return str(int(number)) if number.is_integer() else text


def _field_key(value):
    text = str(value).strip().lstrip("*").split("[", 1)[0]
    return "".join(character.lower() for character in text if character.isalnum())


def resolve_vector_mapping_field(input_file, lookup_table, mapping_field):
    """Resolve a vector field by name or unambiguous lookup-value coverage."""
    import geopandas as gpd

    frame = gpd.read_file(input_file, ignore_geometry=True)
    normalized = _field_key(mapping_field)
    matches = [
        column
        for column in frame.columns
        if _field_key(column) == normalized
    ]
    if len(matches) == 1:
        return str(matches[0])

    lookup = read_categorical_lookup_table(lookup_table)
    lookup_column = next(
        (
            column
            for column in lookup.columns
            if _field_key(column) == normalized
        ),
        None,
    )
    if lookup_column is None:
        raise ValueError(f"Lookup mapping field {mapping_field!r} was not found.")
    lookup_values = {
        key for key in map(_category_key, lookup[lookup_column]) if key is not None
    }
    candidates = []
    for column in frame.columns:
        values = {key for key in map(_category_key, frame[column]) if key is not None}
        if values and values.issubset(lookup_values):
            candidates.append(str(column))
    if len(candidates) == 1:
        return candidates[0]
    available = ", ".join(str(column) for column in frame.columns)
    raise ValueError(
        f"Could not identify the vector field corresponding to lookup field "
        f"{mapping_field!r}. Available vector fields: {available or '<none>'}."
    )


def prepare_categorical_file(
    kind: str,
    input_file: str | Path,
    dem_file: str | Path,
    output_file: str | Path,
    lookup_table: str | Path,
    mapping_field: str,
    class_field: str,
    *,
    is_vector: bool = False,
    classdefinition_file: str | Path | None = None,
    input_crs=None,
    dem_crs=None,
    log: LogCallback | None = None,
) -> Path:
    """Prepare an aligned soil, geology, or land-cover GeoTIFF."""
    if kind not in _OUTPUTS:
        raise ValueError("kind must be 'soil', 'geology', or 'lc'.")

    formatter_name, raster_name, definition_name = _OUTPUTS[kind]
    if classdefinition_file is not None and definition_name is None:
        raise ValueError("Land-cover formatting does not write a classdefinition.")

    input_path = Path(input_file)
    dem_path = Path(dem_file)
    lookup_path = Path(lookup_table)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    definition_path = (
        Path(classdefinition_file) if classdefinition_file is not None else None
    )
    if definition_path is not None:
        definition_path.parent.mkdir(parents=True, exist_ok=True)

    from mhm_tools import pre

    formatter = getattr(pre, formatter_name)
    with capture_messages(log):
        with TemporaryDirectory(
            prefix=f"mhm_qgis_{kind}_", dir=output_path.parent
        ) as temporary:
            temporary_path = Path(temporary)
            lookup_for_tools = lookup_path
            if lookup_path.suffix.lower() == ".txt":
                lookup_for_tools = temporary_path / "lookup.csv"
                read_categorical_lookup_table(lookup_path).to_csv(
                    lookup_for_tools,
                    index=False,
                )
            formatted_input = input_path
            formatter_mapping_field = mapping_field
            formatter_input_crs = input_crs

            if is_vector:
                formatted_input = temporary_path / f"{kind}_rasterized.tif"
                vector_mapping_field = resolve_vector_mapping_field(
                    input_path, lookup_for_tools, mapping_field
                )
                rasterize_kwargs = {
                    "input_file": input_path,
                    "dem_file": dem_path,
                    "output_file": formatted_input,
                    "mapping_field": vector_mapping_field,
                    "lookup_table": lookup_for_tools,
                    "lookup_mapping_field": mapping_field,
                    "lookup_value_field": class_field,
                }
                if input_crs is not None:
                    rasterize_kwargs["input_crs"] = input_crs
                if dem_crs is not None:
                    rasterize_kwargs["dem_crs"] = dem_crs
                pre.rasterize_map_data(**rasterize_kwargs)
                formatter_mapping_field = class_field
                formatter_input_crs = None

            formatter_kwargs = {
                "input_file": formatted_input,
                "dem_file": dem_path,
                "output_path": temporary_path,
                "lookup_table": lookup_for_tools,
                "mapping_field": formatter_mapping_field,
                "class_field": class_field,
                "output_type": "tif",
            }
            if formatter_input_crs is not None:
                formatter_kwargs["input_crs"] = formatter_input_crs
            if dem_crs is not None:
                formatter_kwargs["dem_crs"] = dem_crs
            formatter(**formatter_kwargs)

            raster_output = temporary_path / raster_name
            if not raster_output.is_file():
                raise FileNotFoundError(
                    f"mhm-tools did not create the expected raster: {raster_output}"
                )

            definition_output = (
                temporary_path / definition_name
                if definition_path is not None and definition_name is not None
                else None
            )
            if definition_output is not None and not definition_output.is_file():
                raise FileNotFoundError(
                    "mhm-tools did not create the expected classdefinition: "
                    f"{definition_output}"
                )

            os.replace(raster_output, output_path)
            if definition_output is not None and definition_path is not None:
                os.replace(definition_output, definition_path)

    return output_path


def prepare_land_cover_periods(
    input_path: str | Path,
    dem_file: str | Path,
    output_path: str | Path,
    lookup_table: str | Path,
    mapping_field: str,
    class_field: str,
    *,
    output_type: str,
    input_crs=None,
    dem_crs=None,
    resampling: str = "auto",
    log: LogCallback | None = None,
) -> tuple[Path, ...]:
    """Format a land-cover period manifest through the public mHM-tools API."""
    from mhm_tools import pre

    output = Path(output_path)
    output.mkdir(parents=True, exist_ok=True)
    lookup_path = Path(lookup_table)
    kwargs = {
        "input_path": Path(input_path),
        "dem_file": Path(dem_file),
        "output_path": output,
        "lookup_table": lookup_path,
        "mapping_field": mapping_field,
        "class_field": class_field,
        "output_type": output_type,
        "resampling": resampling,
    }
    if input_crs is not None:
        kwargs["input_crs"] = input_crs
    if dem_crs is not None:
        kwargs["dem_crs"] = dem_crs
    if lookup_path.suffix.lower() == ".txt":
        with TemporaryDirectory(prefix="mhm_qgis_lc_lookup_", dir=output) as temporary:
            lookup_csv = Path(temporary) / "lookup.csv"
            read_categorical_lookup_table(lookup_path).to_csv(lookup_csv, index=False)
            kwargs["lookup_table"] = lookup_csv
            with capture_messages(log):
                outputs = tuple(Path(path) for path in pre.format_lc_periods(**kwargs))
    else:
        with capture_messages(log):
            outputs = tuple(Path(path) for path in pre.format_lc_periods(**kwargs))
    return _require_outputs(outputs)


def prepare_soil_horizons(
    input_path: str | Path,
    dem_file: str | Path,
    output_path: str | Path,
    *,
    output_type: str,
    input_crs=None,
    dem_crs=None,
    resampling: str = "auto",
    composition_step: float = 5.0,
    bulk_density_step: float = 0.1,
    log: LogCallback | None = None,
) -> tuple[Path, Path]:
    """Format a soil-horizon manifest through the public mHM-tools API."""
    from mhm_tools import pre

    output = Path(output_path)
    output.mkdir(parents=True, exist_ok=True)
    kwargs = {
        "input_path": Path(input_path),
        "dem_file": Path(dem_file),
        "output_path": output,
        "output_type": output_type,
        "resampling": resampling,
        "composition_step": composition_step,
        "bulk_density_step": bulk_density_step,
    }
    if input_crs is not None:
        kwargs["input_crs"] = input_crs
    if dem_crs is not None:
        kwargs["dem_crs"] = dem_crs
    with capture_messages(log):
        outputs = tuple(Path(path) for path in pre.format_soil_horizons(**kwargs))
    validated = _require_outputs(outputs)
    if len(validated) != 2:
        raise RuntimeError(
            "mHM-tools did not return soil data and classdefinition outputs."
        )
    return validated[0], validated[1]


def _require_outputs(outputs: tuple[Path, ...]) -> tuple[Path, ...]:
    if not outputs:
        raise RuntimeError("mHM-tools did not return any formatted outputs.")
    missing = [str(path) for path in outputs if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "mHM-tools did not create the expected output(s): " + ", ".join(missing)
        )
    return outputs


def create_dem_derivative_files(
    input_file: str | Path,
    output_folder: str | Path,
    output_extension: str = "tif",
    crs: str | None = None,
    log: LogCallback | None = None,
) -> dict[str, Path]:
    """Create the DEM derivatives through mHM-tools, keyed by derivative name."""
    from mhm_tools.pre.dem_derivatives import create_dem_derivatives

    folder = Path(output_folder)
    folder.mkdir(parents=True, exist_ok=True)
    with capture_messages(log):
        written = create_dem_derivatives(
            input_file=str(input_file),
            output_path=folder,
            output_extension=output_extension,
            crs=crs,
        )
    paths = tuple(Path(path) for path in written)
    _require_outputs(paths)
    return {path.stem: path for path in paths}


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
    log: LogCallback | None = None,
) -> Path:
    """Create ``latlon.nc`` through mHM-tools and return its path."""
    from mhm_tools.pre.latlon import create_latlon

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


__all__ = [
    "capture_messages",
    "create_dem_derivative_files",
    "create_latlon_file",
    "prepare_categorical_file",
    "prepare_land_cover_periods",
    "prepare_soil_horizons",
    "read_categorical_lookup_table",
    "resolve_vector_mapping_field",
]
