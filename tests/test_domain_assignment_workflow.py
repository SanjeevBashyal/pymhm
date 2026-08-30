"""QGIS-light tests for shared domain/gauge assignment validation."""
from __future__ import annotations

from pathlib import Path

import pytest

from mhm_qgis import standalone

standalone.install(force=True)

# isort: off
from qgis.core import QgsVectorLayer  # noqa: E402

from mhm_qgis.Morphology.hydrology.discharge_dialog import (  # noqa: E402
    OutletAssignment,
)
from mhm_qgis.Morphology.watershed.domain_state import (  # noqa: E402
    DOMAIN_MODE_DEM_EXTENT,
    DOMAIN_MODE_SNAPPED,
    gauge_records,
)
from mhm_qgis.Morphology.watershed.domain_workflow import (  # noqa: E402
    DomainWorkflow,
)
# isort: on


class _MainDialog:
    def __init__(self, project_folder: Path):
        self.project_folder = str(project_folder)
        self.messages = []

    def log_message(self, message):
        self.messages.append(message)


class _Processor:
    def __init__(self):
        self.processing_state = {"outputs": {}}


class _SnapProcessor(_Processor):
    def __init__(self, snapped_path: Path):
        super().__init__()
        self._snapped_path = str(snapped_path)
        self.snapped_points_path = str(snapped_path)
        self.removed = []
        self.channel_prepared = False

    def process_channel_network(self):
        pass

    def process_flow_accumulation(self):
        pass

    def fill_dem(self):
        pass

    def snap_points(self):
        pass

    def _ensure_channel_network(self, *_callbacks):
        self.channel_prepared = True
        return True

    def _remove_stale_vector_output(self, path):
        self.removed.append(path)
        Path(path).unlink()

    def _ensure_snapped_points(self, *_callbacks):
        self.snapped_points_path = self._snapped_path
        Path(self.snapped_points_path).write_text("fresh", encoding="utf-8")
        return True


class _PourPoints:
    def __init__(self, source: Path):
        self._source = str(source)

    def source(self):
        return self._source


class _UriTable:
    def __init__(self, path: Path):
        self.path = path

    def source(self):
        return self.path.as_uri() + "?type=csv&delimiter=,"

    def isValid(self):
        return True

    def name(self):
        return self.path.name


def _table(path: Path) -> QgsVectorLayer:
    path.write_text(
        "date,discharge\n2020-01-01,1.5\n2020-01-02,2.5\n",
        encoding="utf-8",
    )
    return QgsVectorLayer(str(path), path.name, "ogr")


def _workflow(tmp_path: Path, outlet_ids=("001", "2")) -> DomainWorkflow:
    return DomainWorkflow(
        _MainDialog(tmp_path),
        _Processor(),
        _PourPoints(tmp_path / "outlets.gpkg"),
        "outlet_id",
        outlet_ids,
    )


def test_gauges_are_all_validated_before_any_observation_is_written(tmp_path):
    workflow = _workflow(tmp_path)
    valid = _table(tmp_path / "valid.csv")
    assignments = [
        OutletAssignment("001", True, False, valid),
        OutletAssignment("2", True, False, None),
    ]

    with pytest.raises(ValueError, match="outlet 2"):
        workflow.validate_gauge_assignments(assignments)

    assert not list(tmp_path.rglob("*.txt"))


def test_assignment_state_separates_filename_text_from_numeric_gauge_id(tmp_path):
    workflow = _workflow(tmp_path)
    layer = _table(tmp_path / "valid.csv")
    assignments = [
        OutletAssignment("001", True, False, layer),
        OutletAssignment("2", False, False, None),
    ]
    prepared = workflow.validate_gauge_assignments(assignments)
    state = workflow.load_synced_state(DOMAIN_MODE_DEM_EXTENT, True)

    workflow.apply_assignment_records(state, assignments, prepared)
    workflow.validate_unique_state_gauge_ids(state)

    assert gauge_records(state) == [
        {
            "outlet_id": "001",
            "gauge_id": 1,
            "gauge_filename": "001.txt",
            "gauge_path": str(
                tmp_path
                / "mhm-plugin"
                / "data"
                / "master"
                / "observation"
                / "streamflow"
                / "001.txt"
            ),
            "domain_ids": [],
        }
    ]


def test_file_provider_uri_is_persisted_as_a_reloadable_local_path(tmp_path):
    workflow = _workflow(tmp_path, ("1",))
    path = tmp_path / "valid.csv"
    _table(path)
    assignments = [OutletAssignment("1", True, False, _UriTable(path))]

    prepared = workflow.validate_gauge_assignments(assignments)
    state = workflow.load_synced_state(DOMAIN_MODE_DEM_EXTENT, True)
    workflow.apply_assignment_records(state, assignments, prepared)

    assert state["outlets"]["1"]["discharge_file"] == str(path.resolve())


def test_numeric_gauge_ids_must_be_unique_even_when_filenames_differ(tmp_path):
    workflow = _workflow(tmp_path, ("001", "1"))
    layer = _table(tmp_path / "valid.csv")

    with pytest.raises(ValueError, match="same numeric gauge ID 1"):
        workflow.validate_gauge_assignments(
            [
                OutletAssignment("001", True, False, layer),
                OutletAssignment("1", True, False, layer),
            ]
        )


def test_assignment_rows_must_still_match_selected_feature_order(tmp_path):
    workflow = _workflow(tmp_path)

    with pytest.raises(ValueError, match="no longer match"):
        workflow._validated_assignment_list(
            [
                OutletAssignment("2", False),
                OutletAssignment("001", False),
            ]
        )


def test_snapped_mode_requires_at_least_one_selected_domain(tmp_path):
    workflow = _workflow(tmp_path)
    state = workflow.load_synced_state(DOMAIN_MODE_SNAPPED, False)

    with pytest.raises(ValueError, match="at least one domain"):
        workflow.require_active_domain(state)

    state["outlets"]["2"]["is_domain"] = True
    workflow.require_active_domain(state)


def test_snapped_points_are_regenerated_after_inputs_may_have_moved(tmp_path):
    snapped = tmp_path / "2_pour_points_snapped.shp"
    snapped.write_text("stale", encoding="utf-8")
    processor = _SnapProcessor(snapped)
    workflow = DomainWorkflow(
        _MainDialog(tmp_path),
        processor,
        _PourPoints(tmp_path / "outlets.gpkg"),
        "outlet_id",
        ["1"],
    )

    result = workflow.regenerate_snapped_points()

    assert processor.channel_prepared
    assert processor.removed == [str(snapped)]
    assert result == str(snapped)
    assert snapped.read_text(encoding="utf-8") == "fresh"

