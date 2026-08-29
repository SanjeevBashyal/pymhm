"""Tests for historical land-use and soil-horizon manifest records."""
from __future__ import annotations

import csv

import pytest

from mhm_qgis.advanced_input_manifests import (
    LandUseInput,
    LandUsePeriod,
    SoilHorizon,
    SoilInput,
    write_land_use_manifest,
    write_soil_manifest,
)


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"raster")
    return path


def test_writes_land_use_period_manifest(tmp_path):
    lookup = _touch(tmp_path / "lookup.csv")
    first = _touch(tmp_path / "lc_1990.tif")
    second = _touch(tmp_path / "lc_2000.tif")
    value = LandUseInput(
        (
            LandUsePeriod(1990, 1999, first),
            LandUsePeriod(2000, 2010, second),
        ),
        lookup,
        "source",
        "class",
    )

    output = write_land_use_manifest(value, tmp_path / "format-data.csv")

    with output.open(newline="", encoding="utf-8") as handle:
        assert list(csv.reader(handle)) == [
            ["StartYear", "EndYear", "FilePath"],
            ["1990", "1999", "lc_1990.tif"],
            ["2000", "2010", "lc_2000.tif"],
        ]


def test_land_use_periods_must_be_continuous(tmp_path):
    value = LandUseInput(
        (
            LandUsePeriod(1990, 1999, _touch(tmp_path / "a.tif")),
            LandUsePeriod(2001, 2010, _touch(tmp_path / "b.tif")),
        ),
        _touch(tmp_path / "lookup.csv"),
        "source",
        "class",
    )

    with pytest.raises(ValueError, match="gap-free"):
        write_land_use_manifest(value, tmp_path / "format-data.csv")


def test_writes_soil_unit_preamble_and_horizons(tmp_path):
    files = [_touch(tmp_path / f"layer_{index}.tif") for index in range(8)]
    value = SoilInput(
        (
            SoilHorizon(1, 0, 50, *files[:4]),
            SoilHorizon(2, 50, 200, *files[4:]),
        ),
        "kg/m3",
    )

    output = write_soil_manifest(value, tmp_path / "format-data.csv")
    lines = output.read_text(encoding="utf-8").splitlines()

    assert lines[0] == "Bulk Density Unit = kg/m3"
    assert lines[1] == (
        "Horizon,Upper Depth,Lower Depth,Clay Layer,Sand Layer,"
        "Silt Layer,Bulk Density Layer"
    )
    assert lines[2] == "1,0,50,layer_0.tif,layer_1.tif,layer_2.tif,layer_3.tif"
    assert lines[3] == "2,50,200,layer_4.tif,layer_5.tif,layer_6.tif,layer_7.tif"


def test_soil_horizons_must_start_at_zero_and_be_contiguous(tmp_path):
    files = [_touch(tmp_path / f"layer_{index}.tif") for index in range(4)]

    with pytest.raises(ValueError, match="start at depth 0"):
        write_soil_manifest(
            SoilInput((SoilHorizon(1, 10, 50, *files),), "g/cm3"),
            tmp_path / "format-data.csv",
        )


def test_soil_state_exposes_lower_depths(tmp_path):
    files = [_touch(tmp_path / f"layer_{index}.tif") for index in range(4)]
    value = SoilInput((SoilHorizon(1, 0, 300, *files),), "g/cm3")

    assert value.lower_depths == (300,)
    assert value.as_dict()["lower_depths"] == [300]
