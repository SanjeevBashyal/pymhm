"""Tests for QGIS-free domain-delineation state handling."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

# isort: off
from pymhm.Morphology.watershed.domain_state import (
    DEM_DOMAIN_ID,
    active_domain_records,
    default_state,
    domain_count,
    gauged_outlet_ids,
    load_state,
    resolve_input_path,
    resolve_output_path,
    save_state,
    serialize_input_path,
    serialize_output_path,
    state_path,
)
# isort: on


def test_missing_and_malformed_state_return_fresh_defaults(tmp_path: Path) -> None:
    project = tmp_path / "project"

    assert load_state(project) == default_state()

    path = state_path(project)
    path.parent.mkdir(parents=True)
    path.write_text("{not json", encoding="utf-8")
    assert load_state(project) == default_state()

    path.write_text('["not", "an", "object"]', encoding="utf-8")
    assert load_state(project) == default_state()


def test_path_serialization_uses_project_and_workspace_roots(
    tmp_path: Path,
) -> None:
    project = tmp_path / "project"
    discharge = project / "inputs" / "gauge.csv"
    external = tmp_path / "external" / "gauge.csv"
    mask = project / "mhm-plugin" / "Z Temp" / "Geometry" / "mask.tif"

    assert serialize_input_path(project, discharge) == "inputs/gauge.csv"
    assert resolve_input_path(project, "inputs/gauge.csv") == discharge.resolve()
    assert serialize_input_path(project, external) == str(external.resolve())
    assert resolve_input_path(project, external) == external.resolve()

    assert (
        serialize_output_path(project, mask)
        == "Z Temp/Geometry/mask.tif"
    )
    assert (
        resolve_output_path(project, "Z Temp/Geometry/mask.tif")
        == mask.resolve()
    )
    with pytest.raises(ValueError, match="plugin workspace"):
        serialize_output_path(project, external)
    with pytest.raises(ValueError, match="plugin workspace"):
        resolve_output_path(project, "../outside.tif")


def test_save_is_atomic_normalized_and_does_not_mutate_input(
    tmp_path: Path,
) -> None:
    project = tmp_path / "project"
    discharge = project / "inputs" / "gauge.csv"
    mask = project / "mhm-plugin" / "Z Temp" / "Geometry" / "mask.tif"
    vector = mask.with_suffix(".gpkg")
    state = {
        "version": 99,
        "pour_points_source": "inputs/outlets.gpkg",
        "outlet_id_field": "station",
        "dem_domain": True,
        "outlets": {
            7: {
                "is_gauged": True,
                "is_domain": True,
                "discharge_file": str(discharge),
                "mask_path": str(mask),
                "vector_path": str(vector),
            }
        },
    }

    path = save_state(project, state)

    assert path == project / "mhm-plugin" / path.name
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["version"] == 1
    assert saved["outlets"]["7"]["discharge_file"] == "inputs/gauge.csv"
    assert saved["outlets"]["7"]["mask_path"] == "Z Temp/Geometry/mask.tif"
    assert saved["outlets"]["7"]["vector_path"] == "Z Temp/Geometry/mask.gpkg"
    assert state["outlets"][7]["mask_path"] == str(mask)
    assert not list(path.parent.glob(f".{path.name}.*.tmp"))
    assert load_state(project) == saved


def test_active_domains_include_optional_dem_domain_and_gauges() -> None:
    state = {
        "dem_domain": True,
        "outlets": {
            "A": {"is_domain": True, "is_gauged": True},
            "B": {"is_domain": False, "is_gauged": True},
            "C": {"domain": True, "gauged": False},
        },
    }

    assert gauged_outlet_ids(state) == ["A", "B"]
    records = active_domain_records(state)
    assert [record["outlet_id"] for record in records] == [
        "A",
        "C",
        DEM_DOMAIN_ID,
    ]
    assert records[-1]["gauged_outlet_ids"] == ["A", "B"]
    assert records[0]["is_dem_domain"] is False
    assert records[-1]["is_dem_domain"] is True
    assert domain_count(state) == 3
