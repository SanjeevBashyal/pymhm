from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from mhm_qgis import standalone

standalone.install(force=True)

from mhm_qgis.lai_temporal import lai_time_step, prepare_lai_temporal  # noqa: E402


def _series(values, start, frequency):
    return xr.DataArray(
        np.asarray(values, dtype="float64"),
        dims=("time",),
        coords={"time": pd.date_range(start, periods=len(values), freq=frequency)},
    )


def test_finer_lai_is_aggregated_by_mean():
    daily = _series(np.arange(31), "2001-01-01", "D")

    result = prepare_lai_temporal(
        daily, "Daily", "Monthly Gridded Data"
    )

    assert result.time_step == -2
    assert result.data.values.tolist() == [15.0]
    assert result.time_bounds.shape == (1, 2)


def test_coarser_lai_is_repeated_without_interpolation():
    annual = _series([2.0, 4.0], "2000-01-01", "YS")

    result = prepare_lai_temporal(
        annual, "Annual", "Monthly Gridded Data"
    )

    assert result.data.sizes["time"] == 24
    assert np.all(result.data.values[:12] == 2.0)
    assert np.all(result.data.values[12:] == 4.0)


def test_long_term_monthly_lai_has_twelve_climatology_values():
    monthly = _series(np.arange(24), "2000-01-01", "MS")

    result = prepare_lai_temporal(
        monthly, "Monthly", "Long Term Mean Monthly Gridded Data"
    )

    assert result.time_step == 1
    assert result.data.sizes["time"] == 12
    assert result.data.values.tolist() == [float(value) for value in range(6, 18)]
    assert result.time_bounds.shape == (12, 2)
    assert lai_time_step("Daily Gridded Data") == -1



def test_lai_output_guard_is_about_disk_not_memory(monkeypatch, tmp_path):
    """Memory is bounded by streaming, so only the output size can still block."""
    import shutil as shutil_module

    from mhm_qgis.Morphology.layers import lai as lai_module
    from mhm_qgis.Morphology.layers.lai import (
        LAI_MAX_BYTES_ENV,
        assert_lai_output_fits,
        lai_grid_byte_size,
    )

    # The A26 case: 468 monthly steps on a 13320 x 6120 L0 grid.
    required = lai_grid_byte_size(468, 6120, 13320)
    assert required > 280 * 1024 ** 3

    # It is allowed when the volume has room for a compressed copy. Free space
    # is stubbed so the result does not depend on the host's disk.
    monkeypatch.delenv(LAI_MAX_BYTES_ENV, raising=False)
    monkeypatch.setattr(
        lai_module.shutil, "disk_usage",
        lambda _p: shutil_module._ntuple_diskusage(0, 0, required),
    )
    assert assert_lai_output_fits(468, 6120, 13320, tmp_path) == required
    assert assert_lai_output_fits(12, 500, 500, tmp_path)

    # An explicit budget still overrides.
    monkeypatch.setenv(LAI_MAX_BYTES_ENV, str(16 * 1024 ** 3))
    with pytest.raises(MemoryError, match="over the"):
        assert_lai_output_fits(468, 6120, 13320, tmp_path)
    assert assert_lai_output_fits(12, 6120, 13320, tmp_path)


def test_lai_output_guard_refuses_when_the_volume_cannot_hold_it(monkeypatch, tmp_path):
    import shutil as shutil_module

    from mhm_qgis.Morphology.layers import lai as lai_module

    monkeypatch.delenv(lai_module.LAI_MAX_BYTES_ENV, raising=False)
    monkeypatch.setattr(
        lai_module.shutil,
        "disk_usage",
        lambda _path: shutil_module._ntuple_diskusage(0, 0, 1024 ** 3),
    )
    with pytest.raises(MemoryError, match="free on the output volume"):
        lai_module.assert_lai_output_fits(468, 6120, 13320, tmp_path)


def test_the_disk_guard_uses_the_measured_bilinear_compression(monkeypatch, tmp_path):
    """Bilinear output compresses ~1.8:1, so the guard must not assume better."""
    import shutil as shutil_module

    from mhm_qgis.Morphology.layers import lai as lai_module

    monkeypatch.delenv(lai_module.LAI_MAX_BYTES_ENV, raising=False)
    # The measured bilinear ratio, not the 70:1 that nearest-neighbour reached.
    assert lai_module.LAI_ASSUMED_COMPRESSION == 1.8

    required = lai_module.lai_grid_byte_size(468, 6120, 13320)
    needed = required / lai_module.LAI_ASSUMED_COMPRESSION
    monkeypatch.setattr(
        lai_module.shutil, "disk_usage",
        lambda _p: shutil_module._ntuple_diskusage(0, 0, int(needed) - 1),
    )
    with pytest.raises(MemoryError, match="free on the output volume"):
        lai_module.assert_lai_output_fits(468, 6120, 13320, tmp_path)

    monkeypatch.setattr(
        lai_module.shutil, "disk_usage",
        lambda _p: shutil_module._ntuple_diskusage(0, 0, int(needed) + 1),
    )
    assert lai_module.assert_lai_output_fits(468, 6120, 13320, tmp_path) == required
