"""Tests for project-local state files in the plugin workspace."""
from __future__ import annotations

from pathlib import Path

from mhm_qgis import standalone

standalone.install(force=True)

# isort: off
from mhm_qgis.core.handlers.state.meteo_outputs import MeteorologyOutputState  # noqa: E402
from mhm_qgis.qgis_bridge.morphology.core.processing_state import ProcessingStateMixin  # noqa: E402
from mhm_qgis.core.handlers.store.layout import workspace_folder  # noqa: E402
# isort: on


class _Dialog:
    def __init__(self, project_folder: Path) -> None:
        self.project_folder = str(project_folder)
        self.messages = []

    def log_message(self, message: str) -> None:
        self.messages.append(message)


class _MorphologyState(ProcessingStateMixin):
    def __init__(self, dialog: _Dialog) -> None:
        self.dialog = dialog
        self.processing_state_filename = "mhm_qgis_processing_state.json"
        self.processing_state = {"version": 1, "outputs": {}, "workflows": {}}
        self.log_message = dialog.log_message


def test_morphology_state_is_saved_in_workspace(tmp_path: Path) -> None:
    dialog = _Dialog(tmp_path)
    state = _MorphologyState(dialog)
    expected = Path(workspace_folder(tmp_path)) / "mhm_qgis_processing_state.json"

    assert Path(state.processing_state_path()) == expected
    state.save_processing_state()

    assert expected.is_file()


def test_meteorology_state_and_output_keys_use_workspace(tmp_path: Path) -> None:
    dialog = _Dialog(tmp_path)
    state = MeteorologyOutputState(dialog.project_folder)
    workspace = Path(workspace_folder(tmp_path))
    output = workspace / "data" / "meteo" / "pre" / "pre.nc"

    assert state.state_path() == workspace / "mhm_qgis_processing_state.json"
    assert state.output_key(output) == "data/meteo/pre/pre.nc"
    state.save({"version": 1, "outputs": {}})

    assert state.state_path().is_file()


def test_grid_contract_is_persisted_and_revalidated_on_resume(tmp_path: Path) -> None:
    from mhm_qgis.grid_resolution import aligned_l0_l2_headers

    dialog = _Dialog(tmp_path)
    state = _MorphologyState(dialog)
    l0_header, l2_header = aligned_l0_l2_headers(
        (99.0, 191.0, 199.0, 291.0),
        {"xllcorner": 100.0, "yllcorner": 200.0, "cellsize": 10.0, "unit": "m"},
        multiplier=3,
        unit="m",
    )

    saved = state.save_grid_contract(l0_header, l2_header, 3)
    assert saved["l2_ratio_to_l0"] == 3

    reloaded = _MorphologyState(_Dialog(tmp_path))
    reloaded.load_processing_state()
    contract = reloaded.saved_grid_contract()
    assert contract["l0_header"] == l0_header
    assert contract["l2_header"] == l2_header
    assert contract["l2_ratio_to_l0"] == 3


def test_an_inconsistent_saved_grid_contract_is_discarded(tmp_path: Path) -> None:
    state = _MorphologyState(_Dialog(tmp_path))
    state.processing_state["grid"] = {
        "l0_header": {
            "ncols": 9, "nrows": 9, "xllcorner": 70.0, "yllcorner": 170.0,
            "cellsize": 10.0,
        },
        # Dimensions no longer equal L2 dimensions times the multiplier.
        "l2_header": {
            "ncols": 5, "nrows": 3, "xllcorner": 70.0, "yllcorner": 170.0,
            "cellsize": 30.0,
        },
        "l2_ratio_to_l0": 3,
    }

    assert state.saved_grid_contract() is None
    assert any("grid contract" in message for message in state.dialog.messages)


def test_saving_the_registry_preserves_sections_other_writers_own(tmp_path: Path) -> None:
    """Wholesale writes used to erase the reuse fingerprints, disabling reuse."""
    from mhm_qgis.core.handlers.state.cache import cached_payload, load_state, store_payload

    state = _MorphologyState(_Dialog(tmp_path))
    state.load_processing_state()
    state.mark_workflow_status("execute_all", "running")

    # Another writer records a stage fingerprint after the registry was loaded.
    store_payload(tmp_path, "stages", "geology", "d1", {"outputs": ["g.tif"]})
    store_payload(tmp_path, "meteo_inspection", "pre", "d2", {"files": []})

    # The registry writes again from its stale in-memory copy.
    state.mark_workflow_status("execute_all", "completed")

    assert cached_payload(tmp_path, "stages", "geology", "d1") == {
        "outputs": ["g.tif"]
    }
    assert cached_payload(tmp_path, "meteo_inspection", "pre", "d2") == {"files": []}
    stored = load_state(tmp_path)
    assert stored["workflows"]["execute_all"]["status"] == "completed"


def test_a_corrupt_registry_file_does_not_lose_the_new_write(tmp_path: Path) -> None:
    from mhm_qgis.core.handlers.state.cache import load_state

    state = _MorphologyState(_Dialog(tmp_path))
    state.load_processing_state()
    path = Path(state.processing_state_path())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{ not json", encoding="utf-8")

    state.mark_workflow_status("execute_all", "completed")

    assert load_state(tmp_path)["workflows"]["execute_all"]["status"] == "completed"
