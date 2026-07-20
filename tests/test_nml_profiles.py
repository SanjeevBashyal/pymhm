"""Pure tests for namelist JSON profiles and version adapters."""
from __future__ import annotations

import importlib
import json
import sys
import types
from pathlib import Path

from nml_tools import json_to_namelist
from pymhm.dependency_bootstrap import dependency_from_requirement_line


def configuration_module(name: str):
    """Import a Configuration module without loading the QGIS processor."""
    root = Path(__file__).resolve().parents[1] / "pymhm" / "Configuration"
    package = types.ModuleType("pymhm.Configuration")
    package.__path__ = [str(root)]
    sys.modules["pymhm.Configuration"] = package
    return importlib.import_module(f"pymhm.Configuration.{name}")


class Dialog:
    def __init__(self, project_folder: Path):
        self.project_folder = str(project_folder)
        self.messages = []

    def log_message(self, message):
        self.messages.append(message)


def test_plugin_dependency_name_maps_to_nml_tools_import() -> None:
    dependency = dependency_from_requirement_line("nml-tools")
    assert dependency is not None
    assert dependency.import_name == "nml_tools"


def test_v513_derived_periods_are_converter_ready() -> None:
    profiles = configuration_module("profile_values")
    values = profiles.namelist_values(
        "5.13",
        "mhm",
        {
            "configproject": {"ndomains": 2},
            "configtime": {
                "warmingdays": [0, 180],
                "evalper": [
                    {
                        "yStart": 1990, "mStart": 1, "dStart": 1,
                        "yEnd": 1993, "mEnd": 12, "dEnd": 31,
                    },
                    {
                        "yStart": 1993, "mStart": 1, "dStart": 1,
                        "yEnd": 1996, "mEnd": 12, "dEnd": 31,
                    },
                ],
                "timestep": [24, 24],
            },
        },
        [],
    )

    assert values["time_periods"]["warming_Days"] == [0, 180]
    assert values["time_periods"]["eval_Per"][1]["yEnd"] == 1996
    rendered = json_to_namelist({"profile": "mhm", "values": values})
    assert "eval_Per(1)%yStart = 1990" in rendered
    assert "eval_Per(2)%dEnd = 31" in rendered
    assert "processCase(11) = 0" in rendered


def test_schema_adapter_restores_names_and_array_axes() -> None:
    profiles = configuration_module("profile_values")
    pages = [{
        "block": "config_input",
        "schema": {
            "properties": {
                "weights": {
                    "type": "array",
                    "x-fortran-shape": [2, "max_domains"],
                }
            }
        },
    }]

    values = profiles.schema_values(
        {
            "configinput": {
                "weights__domain1": [1, 2],
                "weights__domain2": [3, 4],
            }
        },
        pages,
        "6.0",
        "mhm",
    )

    assert values == {"config_input": {"weights": [[1, 3], [2, 4]]}}


def test_v513_parameter_and_output_layouts() -> None:
    profiles = configuration_module("profile_values")
    parameter_pages = [{
        "block": "petm2",
        "schema": {"properties": {"factor": {"type": "array"}}},
    }]
    parameters = profiles.namelist_values(
        "5.13",
        "parameters",
        {"petm2": {"factor": [0.0, 1.0, 0.5, 1, 1]}},
        parameter_pages,
    )
    assert list(parameters) == ["PET0"]

    output_pages = [{
        "block": "output_mhm",
        "schema": {
            "properties": {
                "output_frequency": {"type": "integer"},
                "out_interception": {"type": "boolean"},
                "out_Qsm": {"type": "boolean"},
            }
        },
    }]
    outputs = profiles.namelist_values(
        "5.13",
        "outputs",
        {"outputmhm": {
            "outputfrequency": -1,
            "outinterception": True,
            "outqsm": True,
        }},
        output_pages,
    )
    block = outputs["NLoutputResults"]
    assert block["timeStep_model_outputs"] == -1
    assert block["outputFlxState"][0] is True
    assert block["outputFlxState"][20] is True


def test_nmljson_profiles_and_application_state_are_separate(tmp_path: Path) -> None:
    nmljson = configuration_module("nmljson")
    state_module = configuration_module("state")
    dialog = Dialog(tmp_path)

    store = nmljson.NmlJsonStore(dialog)
    document = store.empty("5.13")
    store.set_profile(
        document,
        "mhm",
        {"time_periods": {"eval_Per": [{"yStart": 2000}]}},
        "5.13",
        editor_values={"configtime": {"evalper": [{"yStart": 2000}]}},
        dimensions={"max_domains": 1},
    )
    store.save(document)

    saved = json.loads((tmp_path / "nmljson.json").read_text())
    assert saved["file_profiles"]["parameters"]["values"] == {}
    assert saved["file_profiles"]["parameters"]["output_file"] == (
        "mhm_parameter.nml")
    profile = saved["file_profiles"]["mhm"]
    assert "eval_Per(1)%yStart = 2000" in json_to_namelist(profile)

    state_path = tmp_path / state_module.CONFIGURATION_STATE_FILENAME
    state_path.write_text(json.dumps({
        "version": 1,
        "mhm_version": "5.13",
        "settings": {"mhm": {"old": 1}},
        "window": {"tab": 2},
    }))
    state_store = state_module.ConfigurationStateStore(dialog)
    state = state_store.load()
    assert state["_legacy_namelist_settings"]["mhm"] == {"old": 1}
    state_store.save(state)
    persisted = json.loads(state_path.read_text())
    assert "settings" not in persisted
    assert "_legacy_namelist_settings" not in persisted
    assert persisted["window"] == {"tab": 2}


def test_processor_writes_profile_and_namelist_through_nml_tools(
        tmp_path: Path) -> None:
    sys.modules.setdefault("processing", types.ModuleType("processing"))
    processor_module = configuration_module("processor")
    dialog = Dialog(tmp_path)
    dialog.check_prerequisites = lambda: True
    processor = processor_module.ConfigurationProcessor(dialog)
    pages = [{
        "block": "simulation",
        "schema": {
            "properties": {
                "periods": {
                    "type": "array",
                    "x-fortran-shape": "max_domains",
                    "items": {
                        "type": "object",
                        "x-fortran-type": "period_t",
                        "properties": {"year": {"type": "integer"}},
                    },
                }
            }
        },
    }]

    processor.write_kind(
        "mhm",
        {"simulation": {"periods": [{"year": 2001}]}},
        pages=pages,
    )

    assert "periods(1)%year = 2001" in (tmp_path / "mhm.nml").read_text()
    saved = json.loads((tmp_path / "nmljson.json").read_text())
    profile = saved["file_profiles"]["mhm"]
    assert profile["values"] == {
        "simulation": {"periods": [{"year": 2001}]}}
    state = json.loads((tmp_path / "pymhm_configuration_state.json").read_text())
    assert "settings" not in state


def test_v513_parameter_profile_preserves_full_template(tmp_path: Path) -> None:
    sys.modules.setdefault("processing", types.ModuleType("processing"))
    processor_module = configuration_module("processor")
    dialog = Dialog(tmp_path)
    dialog.check_prerequisites = lambda: True
    dialog.comboBox_mHMversion = types.SimpleNamespace(
        currentText=lambda: "5.13")
    processor = processor_module.ConfigurationProcessor(dialog)

    processor.write_kind(
        "parameters",
        {"petm2": {"mincorrectionfactorpet": [0.7, 1.3, 0.95, 1, 1]}},
        pages=[],
    )

    profile = json.loads((tmp_path / "nmljson.json").read_text())[
        "file_profiles"]["parameters"]
    assert "PETminus1" in profile["values"]
    assert "PET0" in profile["values"]
    assert "soilmoisture4" in profile["values"]
    assert profile["values"]["PET0"]["minCorrectionFactorPET"][2] == 0.95


def test_cancelled_editor_does_not_write_configuration(
        tmp_path: Path, monkeypatch) -> None:
    sys.modules.setdefault("processing", types.ModuleType("processing"))
    processor_module = configuration_module("processor")
    dialog = Dialog(tmp_path)
    dialog.check_prerequisites = lambda: True
    processor = processor_module.ConfigurationProcessor(dialog)

    class CancelledEditor:
        Accepted = 1

        def __init__(self, *args, **kwargs):
            pass

        def exec_(self):
            return 0

    monkeypatch.setattr(
        processor_module, "NamelistEditorDialog", CancelledEditor)
    saved = processor.open_editor(
        "mhm",
        "Configuration",
        [{"block": "run", "schema": {"properties": {}}, "defaults": {}}],
        False,
    )

    assert saved is False
    assert not (tmp_path / "mhm.nml").exists()
    assert not (tmp_path / "nmljson.json").exists()
