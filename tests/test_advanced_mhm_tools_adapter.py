"""Tests for the period and horizon mHM-tools adapter calls."""
from __future__ import annotations

import sys
from types import ModuleType

from mhm_qgis.mhm_tools_adapter import (
    prepare_land_cover_periods,
    prepare_soil_horizons,
)


def _install_fake(monkeypatch, **functions):
    package = ModuleType("mhm_tools")
    pre = ModuleType("mhm_tools.pre")
    for name, function in functions.items():
        setattr(pre, name, function)
    package.pre = pre
    monkeypatch.setitem(sys.modules, "mhm_tools", package)
    monkeypatch.setitem(sys.modules, "mhm_tools.pre", pre)


def test_land_cover_period_adapter_uses_public_api(tmp_path, monkeypatch):
    calls = []

    def format_lc_periods(**kwargs):
        calls.append(kwargs)
        output = kwargs["output_path"] / "lc_periods.nc"
        output.write_bytes(b"nc")
        return (output,)

    _install_fake(monkeypatch, format_lc_periods=format_lc_periods)
    output = prepare_land_cover_periods(
        tmp_path / "manifest",
        tmp_path / "dem.tif",
        tmp_path / "output",
        tmp_path / "lookup.csv",
        "source",
        "class",
        output_type="nc",
    )

    assert [path.name for path in output] == ["lc_periods.nc"]
    assert calls[0]["input_path"] == tmp_path / "manifest"
    assert calls[0]["resampling"] == "auto"


def test_soil_horizon_adapter_returns_data_and_definition(tmp_path, monkeypatch):
    calls = []

    def format_soil_horizons(**kwargs):
        calls.append(kwargs)
        data = kwargs["output_path"] / "soil_class.asc"
        definition = kwargs["output_path"] / "soil_classdefinition.txt"
        data.write_bytes(b"asc")
        definition.write_text("definition", encoding="utf-8")
        return data, definition

    _install_fake(monkeypatch, format_soil_horizons=format_soil_horizons)
    data, definition = prepare_soil_horizons(
        tmp_path / "manifest",
        tmp_path / "dem.tif",
        tmp_path / "output",
        output_type="asc",
    )

    assert data.name == "soil_class.asc"
    assert definition.name == "soil_classdefinition.txt"
    assert calls[0]["composition_step"] == 5.0
    assert calls[0]["bulk_density_step"] == 0.1
