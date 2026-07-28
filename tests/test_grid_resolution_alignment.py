"""Focused checks for the shared L0/L2 extent contract."""

from pymhm.grid_resolution import (
    header_for_aligned_bounds,
    header_for_existing_bounds,
)


def test_aligned_extent_is_minimal_and_has_integer_l0_l2_dimensions():
    l2 = header_for_aligned_bounds(
        bounds=(115.0, 985.0, 218.0, 812.0),
        cellsize=30.0,
        anchor_x=100.0,
        anchor_y=200.0,
        unit="m",
    )
    l0 = header_for_existing_bounds(l2, cellsize=10.0, unit="m")

    assert l2 == {
        "ncols": 30,
        "nrows": 21,
        "xllcorner": 100.0,
        "yllcorner": 200.0,
        "cellsize": 30.0,
        "nodata_value": -9999.0,
        "unit": "m",
        "cellsize_precision": 8,
    }
    assert l0["ncols"] == 3 * l2["ncols"]
    assert l0["nrows"] == 3 * l2["nrows"]
