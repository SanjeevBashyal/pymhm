"""Focused integration tests for saved domain workflow choices."""
from __future__ import annotations

from pathlib import Path

import pytest

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

# isort: off
from pymhm.Morphology.core.predecessors import PredecessorMixin  # noqa: E402
from pymhm.Morphology.hydrology.gauge import GaugePositionMixin  # noqa: E402
from pymhm.Morphology.hydrology.outlets import OutletCountMixin  # noqa: E402
from pymhm.Morphology.orchestration.reset_geometry import (  # noqa: E402
    ResetGeometryMixin,
)
from pymhm.Morphology.watershed.domain_state import (  # noqa: E402
    load_state,
    save_state,
)
from pymhm.project_layout import geometry_folder  # noqa: E402
from pymhm.vSpecific.common import domain_count  # noqa: E402
# isort: on


class _Label:
    def __init__(self):
        self.value = ""

    def setText(self, value):
        self.value = value


class _Layer:
    def __init__(self, features=None, count=0):
        self._features = list(features or [])
        self._count = count

    def isValid(self):
        return True

    def getFeatures(self):
        return iter(self._features)

    def featureCount(self):
        return self._count


class _Feature:
    def __init__(self, outlet_id):
        self.outlet_id = outlet_id

    def attribute(self, _field):
        return self.outlet_id

    def geometry(self):
        return _Geometry()


class _Geometry:
    def isEmpty(self):
        return False


def test_saved_state_controls_gauge_and_domain_counts(tmp_path: Path) -> None:
    save_state(
        tmp_path,
        {
            "pour_points_source": "inputs/outlets.gpkg",
            "outlet_id_field": "outlet_id",
            "dem_domain": True,
            "outlets": {
                "1": {"is_gauged": True, "is_domain": True},
                "2": {"is_gauged": False, "is_domain": True},
            },
        },
    )
    label = _Label()
    selector = type(
        "Selector",
        (),
        {"currentLayer": lambda self: _Layer(count=20)},
    )()
    dialog = type(
        "Dialog",
        (),
        {
            "project_folder": str(tmp_path),
            "mMapLayerComboBox_pour_points": selector,
            "label_numberOfGaugedOutletsValue": label,
        },
    )()
    harness = type("Harness", (), {"dialog": dialog})()

    assert OutletCountMixin.update_gauged_outlet_count(harness) == "1"
    assert label.value == "1"
    assert domain_count(dialog) == 3


def test_configured_project_keeps_one_domain_minimum(tmp_path: Path) -> None:
    save_state(
        tmp_path,
        {
            "outlet_id_field": "outlet_id",
            "outlets": {"1": {"is_domain": False}},
        },
    )
    dialog = type(
        "Dialog", (), {"project_folder": str(tmp_path)}
    )()

    assert domain_count(dialog) == 1


def test_gauge_features_only_include_saved_gauged_ids() -> None:
    layer = _Layer([_Feature("1"), _Feature("2")])

    features = GaugePositionMixin._gauge_features(
        object(),
        layer,
        "outlet_id",
        ["2"],
    )

    assert [feature["station_text"] for feature in features] == ["2"]
    with pytest.raises(ValueError, match="No outlets are marked as gauged"):
        GaugePositionMixin._gauge_features(
            object(),
            layer,
            "outlet_id",
            [],
        )


def test_configured_state_prevents_legacy_automatic_delineation(
    tmp_path: Path,
) -> None:
    save_state(
        tmp_path,
        {
            "outlet_id_field": "outlet_id",
            "outlets": {"1": {"is_domain": True}},
        },
    )
    messages = []
    harness = type(
        "Harness",
        (),
        {
            "dialog": type(
                "Dialog", (), {"project_folder": str(tmp_path)}
            )(),
            "_restore_existing_path": lambda self, *args: None,
            "log_message": messages.append,
        },
    )()

    assert not PredecessorMixin._ensure_merged_watershed(
        harness,
        lambda: None,
        lambda: None,
        lambda: None,
        lambda: None,
        lambda: None,
    )
    assert any("Domain Delineator" in message for message in messages)


def test_legacy_delineation_remains_without_domain_state(
    tmp_path: Path,
) -> None:
    output = Path(geometry_folder(tmp_path)) / "legacy.shp"

    class Harness:
        def __init__(self):
            self.dialog = type(
                "Dialog", (), {"project_folder": str(tmp_path)}
            )()
            self.merged_watershed_path = None

        def _restore_existing_path(self, *args):
            return None

        def _ensure_snapped_points(self, *args):
            return True

        def without_layer_loading(self, callback):
            callback()

        def log_message(self, _message):
            pass

    harness = Harness()

    def delineate():
        output.parent.mkdir(parents=True, exist_ok=True)
        output.touch()
        harness.merged_watershed_path = str(output)

    assert PredecessorMixin._ensure_merged_watershed(
        harness,
        delineate,
        lambda: None,
        lambda: None,
        lambda: None,
        lambda: None,
    )


def test_geometry_reset_invalidates_only_domain_outputs(
    tmp_path: Path,
) -> None:
    mask = Path(geometry_folder(tmp_path)) / "Watersheds" / "mask.tif"
    vector = mask.with_suffix(".gpkg")
    save_state(
        tmp_path,
        {
            "pour_points_source": "inputs/outlets.gpkg",
            "outlet_id_field": "outlet_id",
            "dem_domain": True,
            "outlets": {
                "7": {
                    "is_gauged": True,
                    "is_domain": True,
                    "discharge_file": "inputs/gauge.csv",
                    "threshold_cells": 42,
                    "picked": {"x": 1.0, "y": 2.0},
                    "catchment_area_m2": 123.0,
                    "mask_path": str(mask),
                    "vector_path": str(vector),
                }
            },
        },
    )
    messages = []
    harness = type(
        "Harness",
        (),
        {
            "dialog": type(
                "Dialog", (), {"project_folder": str(tmp_path)}
            )(),
            "log_message": messages.append,
        },
    )()

    ResetGeometryMixin._invalidate_domain_delineation_outputs(harness)

    state = load_state(tmp_path)
    record = state["outlets"]["7"]
    assert record["is_gauged"] is True
    assert record["is_domain"] is True
    assert record["discharge_file"] == "inputs/gauge.csv"
    assert record["threshold_cells"] == 42
    assert "picked" not in record
    assert "catchment_area_m2" not in record
    assert "mask_path" not in record
    assert "vector_path" not in record
    assert state["dem_domain"] is True
