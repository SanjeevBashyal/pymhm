"""The LAI L0 writer must place values exactly while bounding its memory."""
from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from mhm_qgis import standalone  # noqa: E402

standalone.install(force=True)

from mhm_qgis.Morphology.layers.lai_l0 import (  # noqa: E402
    CHUNK_ROWS,
    NODATA,
    PAD_VALUE,
    ResampleSampler,
    WindowSampler,
    block_row_count,
    bracketing_weights,
    lai_window_offsets,
    longitude_to_source_convention,
    nearest_source_indices,
    output_byte_size,
    separable_axes,
    stream_lai_grid,
)

xr = pytest.importorskip("xarray")


SOURCE_LAT = np.array([10.0, 20.0, 30.0])
SOURCE_LON = np.array([100.0, 110.0])
STEPS = 3


def _source_values():
    """Return a 3-step cube whose values encode their own source cell."""
    base = np.array([[11.0, 12.0], [21.0, 22.0], [31.0, 32.0]])
    return np.stack([base + 100 * step for step in range(STEPS)])


def _target(nrows=7, ncols=5):
    # Target cell centres deliberately offset from the source points.
    y_centers = np.linspace(31.0, 9.0, nrows)
    x_centers = np.linspace(99.0, 111.0, ncols)
    return x_centers, y_centers


def _coordinates(x_centers, y_centers):
    return xr.Dataset(
        data_vars={
            "time_bnds": (
                ("time", "bnds"),
                np.column_stack(
                    [np.arange(STEPS), np.arange(1, STEPS + 1)]
                ).astype(float),
            )
        },
        coords={
            "time": np.arange(1, STEPS + 1, dtype="int32"),
            "yc": y_centers,
            "xc": x_centers,
            "bnds": np.arange(2, dtype="int8"),
        },
        attrs={"description": "test"},
    )


def _write(tmp_path, mask=None, nrows=7, ncols=5, block_bytes=None,
           method="nearest"):
    x_centers, y_centers = _target(nrows, ncols)

    def row_lonlat(start, stop):
        rows = y_centers[start:stop]
        return (
            np.broadcast_to(x_centers, (rows.size, ncols)),
            np.repeat(rows[:, None], ncols, axis=1),
        )

    kwargs = {} if block_bytes is None else {"block_bytes": block_bytes}
    path = stream_lai_grid(
        tmp_path / "lai.nc",
        coordinate_dataset=_coordinates(x_centers, y_centers),
        sampler=ResampleSampler(
            _source_values(), SOURCE_LAT, SOURCE_LON, method=method
        ),
        x_centers=x_centers,
        y_centers=y_centers,
        row_lonlat=row_lonlat,
        lai_attrs={"long_name": "leaf area index", "units": "1"},
        mask=mask,
        **kwargs,
    )
    return path, x_centers, y_centers


def _expected(x_centers, y_centers):
    values = _source_values()
    rows = np.abs(y_centers[:, None] - SOURCE_LAT[None, :]).argmin(axis=1)
    cols = np.abs(x_centers[:, None] - SOURCE_LON[None, :]).argmin(axis=1)
    return np.stack([values[step][np.ix_(rows, cols)] for step in range(STEPS)])


def test_nearest_indices_pick_the_closest_source_cell():
    assert list(nearest_source_indices([9.0, 14.0, 16.0, 99.0], SOURCE_LAT)) == [
        0, 0, 1, 2,
    ]
    # A single-point axis always resolves to its only cell.
    assert list(nearest_source_indices([5.0, 50.0], [7.0])) == [0, 0]


def test_longitudes_follow_the_source_convention():
    positive_only = np.array([0.0, 359.0])
    np.testing.assert_allclose(
        longitude_to_source_convention(np.array([-1.0]), positive_only), [359.0]
    )
    signed = np.array([-179.0, 179.0])
    np.testing.assert_allclose(
        longitude_to_source_convention(np.array([-1.0]), signed), [-1.0]
    )


def test_separable_axes_collapse_a_geographic_block_and_reject_a_curved_one():
    lon = np.broadcast_to(np.array([1.0, 2.0, 3.0]), (2, 3))
    lat = np.repeat(np.array([[10.0], [20.0]]), 3, axis=1)
    axes = separable_axes(lon, lat)
    assert axes is not None
    np.testing.assert_array_equal(axes[0], [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(axes[1], [10.0, 20.0])
    # A projected grid curves, so it cannot be reduced to two axes.
    assert separable_axes(lon + np.array([[0.0], [0.5]]), lat) is None


def test_streamed_values_match_a_direct_nearest_neighbour_placement(tmp_path):
    path, x_centers, y_centers = _write(tmp_path)

    with xr.open_dataset(path) as result:
        np.testing.assert_array_equal(
            result["lai"].values, _expected(x_centers, y_centers)
        )
        assert result["lai"].attrs["long_name"] == "leaf area index"
        assert result.attrs["description"] == "test"
        np.testing.assert_allclose(result["xc"].values, x_centers)
        # The 2-D lat/lon variables are written block by block.
        np.testing.assert_allclose(result["lon"].values[0], x_centers)
        np.testing.assert_allclose(result["lat"].values[:, 0], y_centers)


def test_many_small_blocks_produce_the_same_file_as_one_block(tmp_path):
    """Row blocking must not change a single value at a block boundary."""
    single, x_centers, y_centers = _write(tmp_path / "one", nrows=40)
    with xr.open_dataset(single) as result:
        expected = result["lai"].values.copy()

    blocked, _x, _y = _write(tmp_path / "many", nrows=40, block_bytes=1)
    with xr.open_dataset(blocked) as result:
        np.testing.assert_array_equal(result["lai"].values, expected)


def test_a_mask_leaves_nodata_outside_the_watershed(tmp_path):
    mask = np.zeros((7, 5), dtype=bool)
    mask[2:5, 1:4] = True
    path, x_centers, y_centers = _write(tmp_path, mask=mask)

    with xr.open_dataset(path, mask_and_scale=False) as result:
        values = result["lai"].values
    expected = _expected(x_centers, y_centers)
    np.testing.assert_array_equal(values[:, mask], expected[:, mask])
    assert np.all(values[:, ~mask] == NODATA)


def test_block_sizing_stays_on_whole_chunks_and_is_bounded():
    assert block_row_count(13320) % CHUNK_ROWS == 0
    assert block_row_count(13320, block_bytes=1) == CHUNK_ROWS
    assert block_row_count(1, block_bytes=10 ** 12) == 4096
    # The A26 case, stated so the size is visible in the test suite.
    assert output_byte_size(468, 6120, 13320) == 305204889600


def test_a_failed_write_leaves_no_partial_output(tmp_path):
    x_centers, y_centers = _target()

    def boom(_start, _stop):
        raise RuntimeError("row generation failed")

    with pytest.raises(RuntimeError, match="row generation failed"):
        stream_lai_grid(
            tmp_path / "lai.nc",
            coordinate_dataset=_coordinates(x_centers, y_centers),
            sampler=ResampleSampler(_source_values(), SOURCE_LAT, SOURCE_LON),
            x_centers=x_centers,
            y_centers=y_centers,
            row_lonlat=boom,
        )
    assert list(tmp_path.iterdir()) == []


def test_bilinear_weights_bracket_and_clamp_at_the_edges():
    lower, upper, weight = bracketing_weights(
        [5.0, 10.0, 15.0, 25.0, 99.0], SOURCE_LAT
    )
    assert list(lower) == [0, 0, 0, 1, 1]
    assert list(upper) == [1, 1, 1, 2, 2]
    # Below and above the source extent clamp to the edge value.
    np.testing.assert_allclose(weight, [0.0, 0.0, 0.5, 0.5, 1.0])


def test_bilinear_interpolates_between_source_cells(tmp_path):
    """A target point midway between four source cells gets their mean."""
    x_centers = np.array([105.0])          # midway between lon 100 and 110
    y_centers = np.array([15.0])           # midway between lat 10 and 20

    def row_lonlat(start, stop):
        rows = y_centers[start:stop]
        return (
            np.broadcast_to(x_centers, (rows.size, 1)),
            np.repeat(rows[:, None], 1, axis=1),
        )

    path = stream_lai_grid(
        tmp_path / "lai.nc",
        coordinate_dataset=_coordinates(x_centers, y_centers),
        sampler=ResampleSampler(
            _source_values(), SOURCE_LAT, SOURCE_LON, method="bilinear"
        ),
        x_centers=x_centers,
        y_centers=y_centers,
        row_lonlat=row_lonlat,
    )
    with xr.open_dataset(path) as result:
        # Corners of step 0 are 11, 12, 21, 22 -> mean 16.5.
        assert float(result["lai"][0, 0, 0].values) == pytest.approx(16.5)
        assert float(result["lai"][1, 0, 0].values) == pytest.approx(116.5)


def test_bilinear_ignores_missing_source_corners(tmp_path):
    """One NaN corner must not poison the cell; the rest are renormalised."""
    values = _source_values()
    values[0, 0, 0] = np.nan          # lat 10, lon 100 for the first step
    x_centers = np.array([105.0])
    y_centers = np.array([15.0])

    def row_lonlat(start, stop):
        rows = y_centers[start:stop]
        return (
            np.broadcast_to(x_centers, (rows.size, 1)),
            np.repeat(rows[:, None], 1, axis=1),
        )

    path = stream_lai_grid(
        tmp_path / "lai.nc",
        coordinate_dataset=_coordinates(x_centers, y_centers),
        sampler=ResampleSampler(values, SOURCE_LAT, SOURCE_LON, method="bilinear"),
        x_centers=x_centers,
        y_centers=y_centers,
        row_lonlat=row_lonlat,
    )
    with xr.open_dataset(path) as result:
        # Remaining corners 12, 21, 22 each carry weight 0.25 -> mean 18.333.
        assert float(result["lai"][0, 0, 0].values) == pytest.approx(
            (12.0 + 21.0 + 22.0) / 3.0
        )


def test_an_unsupported_resampling_method_is_rejected():
    with pytest.raises(ValueError, match="Unsupported LAI resampling method"):
        ResampleSampler(_source_values(), SOURCE_LAT, SOURCE_LON, method="cubic")


def test_window_offsets_locate_an_expanded_target_in_its_source():
    cell = 0.0008333333333333333
    # Staged DEM grid: 13201 x 6001 starting at the A26 DEM origin.
    source_x = 78.99958333333333 + (np.arange(13201) + 0.5) * cell
    source_y = 31.000416666666666 - (np.arange(6001) + 0.5) * cell
    target = {
        "ncols": 13320, "nrows": 6120,
        "xllcorner": 78.99958333333333,
        "yllcorner": 25.999583333333334,
        "cellsize": cell,
    }
    # The target extends 119 rows above the DEM, so the DEM starts at row 119.
    assert lai_window_offsets(source_x, source_y, target) == (-119, 0)


def test_window_offsets_reject_a_grid_that_would_need_resampling():
    cell = 0.001
    source_x = 10.0 + (np.arange(5) + 0.5) * cell
    source_y = 20.0 - (np.arange(5) + 0.5) * cell
    base = {"ncols": 5, "nrows": 5, "xllcorner": 10.0, "yllcorner": 19.995}

    lai_window_offsets(source_x, source_y, {**base, "cellsize": cell})
    with pytest.raises(ValueError, match="does not match the target cell size"):
        lai_window_offsets(source_x, source_y, {**base, "cellsize": 0.002})
    with pytest.raises(ValueError, match="not aligned to the target"):
        lai_window_offsets(
            source_x + cell / 3.0, source_y, {**base, "cellsize": cell}
        )


def test_window_copy_pads_outside_the_source_and_keeps_values(tmp_path):
    """The crop stage copies cells verbatim and pads the expansion with 0 LAI."""
    from netCDF4 import Dataset

    staged = tmp_path / "lai_dem.nc"
    values = np.arange(2 * 3 * 4, dtype="float64").reshape(2, 3, 4)
    with Dataset(staged, "w") as handle:
        handle.createDimension("time", 2)
        handle.createDimension("yc", 3)
        handle.createDimension("xc", 4)
        handle.createVariable("lai", "f8", ("time", "yc", "xc"))[:] = values

    x_centers = np.arange(6, dtype="float64")   # one extra column each side
    y_centers = np.arange(4, -1, -1, dtype="float64")

    def row_lonlat(start, stop):
        rows = y_centers[start:stop]
        return (
            np.broadcast_to(x_centers, (rows.size, 6)),
            np.repeat(rows[:, None], 6, axis=1),
        )

    with Dataset(staged, "r") as handle:
        path = stream_lai_grid(
            tmp_path / "lai_crop.nc",
            coordinate_dataset=xr.Dataset(
                data_vars={"time_bnds": (("time", "bnds"), np.zeros((2, 2)))},
                coords={
                    "time": np.arange(1, 3, dtype="int32"),
                    "yc": y_centers, "xc": x_centers,
                    "bnds": np.arange(2, dtype="int8"),
                },
            ),
            sampler=WindowSampler(
                handle.variables["lai"], row_offset=-1, column_offset=-1
            ),
            x_centers=x_centers,
            y_centers=y_centers,
            row_lonlat=row_lonlat,
        )

    with xr.open_dataset(path, mask_and_scale=False) as result:
        placed = result["lai"].values
    assert placed.shape == (2, 5, 6)
    np.testing.assert_array_equal(placed[:, 1:4, 1:5], values)
    # Widening the extent must not create nodata holes inside the domain.
    assert PAD_VALUE != NODATA
    assert np.all(placed[:, 0, :] == PAD_VALUE)
    assert np.all(placed[:, 4, :] == PAD_VALUE)
    assert np.all(placed[:, :, 0] == PAD_VALUE)
    assert np.all(placed[:, :, 5] == PAD_VALUE)


def _dem_raster(path, ncols=6, nrows=4, cell=0.5, xmin=99.0, ymax=31.0):
    """Write a small geographic DEM defining the stage-1 target grid."""
    from osgeo import gdal, osr

    dataset = gdal.GetDriverByName("GTiff").Create(
        str(path), ncols, nrows, 1, gdal.GDT_Float32
    )
    dataset.SetGeoTransform((xmin, cell, 0, ymax, 0, -cell))
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(4326)
    dataset.SetProjection(crs.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(np.zeros((nrows, ncols), dtype="float32"))
    dataset = None


def _lai_source(path):
    """Write a monthly LAI source covering the DEM extent."""
    times = np.array(
        [f"2001-{month:02d}-01" for month in range(1, 13)], dtype="datetime64[ns]"
    )
    values = np.arange(12 * 3 * 2, dtype="float32").reshape(12, 3, 2)
    xr.DataArray(
        values,
        dims=("time", "lat", "lon"),
        coords={"time": times, "lat": [30.0, 30.5, 31.0], "lon": [99.5, 101.5]},
        name="lai",
        attrs={"units": "1"},
    ).to_dataset().to_netcdf(path)


def test_run_lai_resample_stages_on_the_dem_grid_from_primitives(tmp_path):
    """The QgsTask entry point must work from paths and strings alone."""
    from mhm_qgis.Morphology.layers.lai_source import run_lai_resample

    dem = tmp_path / "dem.tif"
    source = tmp_path / "lai_src.nc"
    _dem_raster(dem)
    _lai_source(source)

    output = run_lai_resample({
        "source_path": str(source),
        "source_variable": None,
        "output_path": str(tmp_path / "lai_dem.nc"),
        "filled_dem": str(dem),
        "crs_string": "EPSG:4326",
        "input_resolution": "Monthly",
        "target_timestep": "Long Term Mean Monthly Gridded Data",
    })

    with xr.open_dataset(output) as result:
        # The staged matrix matches the DEM exactly, with 12 climatology months.
        assert result["lai"].shape == (12, 4, 6)
        np.testing.assert_allclose(
            result["xc"].values, 99.0 + (np.arange(6) + 0.5) * 0.5
        )
        np.testing.assert_allclose(
            result["yc"].values, 31.0 - (np.arange(4) + 0.5) * 0.5
        )
        values = result["lai"].values
        assert np.all(np.isfinite(values))
        # Bilinear must produce values the source does not contain.
        with xr.open_dataset(source) as raw:
            assert not set(np.unique(values).tolist()).issubset(
                set(np.unique(raw["lai"].values).tolist())
            )


def test_a_cancelled_task_stops_the_stream_and_leaves_no_output(tmp_path):
    from mhm_qgis.Morphology.layers.lai_source import run_lai_resample

    dem = tmp_path / "dem.tif"
    source = tmp_path / "lai_src.nc"
    _dem_raster(dem)
    _lai_source(source)

    class _Cancelled:
        def isCanceled(self):
            return True

        def setProgress(self, _value):
            pass

    output = tmp_path / "lai_dem.nc"
    with pytest.raises(RuntimeError, match="cancelled"):
        run_lai_resample({
            "source_path": str(source),
            "source_variable": None,
            "output_path": str(output),
            "filled_dem": str(dem),
            "crs_string": "EPSG:4326",
            "input_resolution": "Monthly",
            "target_timestep": "Long Term Mean Monthly Gridded Data",
        }, task=_Cancelled())
    assert not output.exists()


def test_blank_cells_are_filled_with_zero_when_asked(tmp_path):
    """Execute All fills gaps with 0: no leaf area is zero, not missing."""
    values = _source_values()
    values[:, 0, 0] = np.nan          # lat 10 / lon 100 missing in every step
    x_centers = np.array([100.0])     # exactly on the missing source cell
    y_centers = np.array([10.0])

    def row_lonlat(start, stop):
        rows = y_centers[start:stop]
        return (
            np.broadcast_to(x_centers, (rows.size, 1)),
            np.repeat(rows[:, None], 1, axis=1),
        )

    def write(blank_fill):
        return stream_lai_grid(
            tmp_path / f"lai_{blank_fill}.nc",
            coordinate_dataset=_coordinates(x_centers, y_centers),
            sampler=ResampleSampler(
                values, SOURCE_LAT, SOURCE_LON,
                method="bilinear", blank_fill=blank_fill,
            ),
            x_centers=x_centers,
            y_centers=y_centers,
            row_lonlat=row_lonlat,
        )

    with xr.open_dataset(write(0.0)) as result:
        assert np.all(result["lai"].values == 0.0)
    # Without the fill the same cell stays missing, so the fill is what changed it.
    with xr.open_dataset(write(None), mask_and_scale=False) as result:
        assert np.all(np.isnan(result["lai"].values))


def test_the_fill_leaves_real_values_untouched(tmp_path):
    values = _source_values()
    values[0, 0, 0] = np.nan
    x_centers = np.array([100.0, 110.0])
    y_centers = np.array([10.0])

    def row_lonlat(start, stop):
        rows = y_centers[start:stop]
        return (
            np.broadcast_to(x_centers, (rows.size, 2)),
            np.repeat(rows[:, None], 2, axis=1),
        )

    path = stream_lai_grid(
        tmp_path / "lai.nc",
        coordinate_dataset=_coordinates(x_centers, y_centers),
        sampler=ResampleSampler(
            values, SOURCE_LAT, SOURCE_LON, method="nearest", blank_fill=0.0),
        x_centers=x_centers,
        y_centers=y_centers,
        row_lonlat=row_lonlat,
    )
    with xr.open_dataset(path) as result:
        first = result["lai"].values[0, 0]
    # Only the missing cell became 0; its neighbour kept its value.
    assert first[0] == 0.0
    assert first[1] == 12.0


def test_the_resample_stage_fills_blanks_with_zero_by_default():
    """`run_lai_resample` is the Execute All entry point, so 0 is its default."""
    import inspect

    from mhm_qgis.Morphology.layers import lai_source

    source = inspect.getsource(lai_source.run_lai_resample)
    assert 'blank_fill=options.get("blank_fill", 0.0)' in source
