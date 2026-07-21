"""Version-contract tests for the embedded nml-tools GUI projects."""

from __future__ import annotations

import json
import os
from pathlib import Path

from nml_tools import json_to_namelist
from nml_tools.gui.model import (
    empty_document,
    load_project,
    merge_initial_values,
    profile_values,
    render_profile,
)

from pymhm.vSpecific import build_dimensions, build_initial_values


SCHEMAS = Path(__file__).resolve().parents[1] / "pymhm" / "nml-schemas"


class Dialog:
    def __init__(self, project_folder: Path):
        self.project_folder = str(project_folder)

    def current_l1_resolution(self):
        return 1000.0

    def current_l11_resolution(self):
        return 2000.0

    def current_grid_unit(self):
        return "m"


class Field:
    def __init__(self, name: str):
        self._name = name

    def name(self):
        return self._name


class Feature:
    def __init__(self, is_domain: bool):
        self.is_domain = is_domain

    def attribute(self, name: str):
        return self.is_domain if name == "IS_DOMAIN" else None


class Layer:
    def isValid(self):
        return True

    def fields(self):
        return [Field("IS_DOMAIN")]

    def getFeatures(self):
        return [Feature(True), Feature(False), Feature(True)]


class LayerWidget:
    def currentLayer(self):
        return Layer()


def _document(version: str, dialog: Dialog):
    key = "v5.13" if version.startswith("5") else "v6"
    project = load_project(SCHEMAS / key, dialog.project_folder)
    dimensions = build_dimensions(version, dialog)
    document = empty_document(project)
    document["dimensions"] = dimensions
    document = merge_initial_values(
        document,
        build_initial_values(version, dialog),
        project,
    )
    return project, dimensions, document


def test_version_projects_have_model_filenames_and_dimensions(tmp_path: Path) -> None:
    v5 = load_project(SCHEMAS / "v5.13", tmp_path)
    v6 = load_project(SCHEMAS / "v6", tmp_path)

    assert v5.default_dimensions == {
        "max_domains": 100,
        "max_geo_units": 25,
        "max_layers": 10,
    }
    assert {profile.name: profile.default_file for profile in v5.profiles} == {
        "main": "mhm.nml",
        "output": "mhm_outputs.nml",
        "parameter": "mhm_parameter.nml",
    }
    assert {profile.name: profile.default_file for profile in v6.profiles} == {
        "main": "mhm.nml",
        "output": "mhm_outputs.nml",
        "parameter": "mhm_parameters.nml",
    }


def test_initial_values_render_version_specific_main_groups(tmp_path: Path) -> None:
    dialog = Dialog(tmp_path)

    v5, v5_dimensions, v5_document = _document("5.13", dialog)
    v5_main = v5.profile("main")
    rendered_v5 = render_profile(
        v5,
        v5_main,
        profile_values(v5_document, v5_main),
        v5_dimensions,
    )
    assert "&mainconfig\n" in rendered_v5
    assert "&directories_mHM\n" in rendered_v5
    assert "resolution_Hydrology(1) = 1000.0" in rendered_v5
    assert 'dir_Morpho(1) = "data/static/morph/"' in rendered_v5
    assert "&config_input" not in rendered_v5

    v6, v6_dimensions, v6_document = _document("6.0", dialog)
    v6_main = v6.profile("main")
    rendered_v6 = render_profile(
        v6,
        v6_main,
        profile_values(v6_document, v6_main),
        v6_dimensions,
    )
    assert "&config_input\n" in rendered_v6
    assert "&config_resolution\n" in rendered_v6
    assert "hydro(1) = 1000.0" in rendered_v6
    assert 'dem_path(1) = "data/static/morph/dem.asc"' in rendered_v6


def test_domain_dimensions_and_paths_follow_plugin_domains(tmp_path: Path) -> None:
    dialog = Dialog(tmp_path)
    dialog.mMapLayerComboBox_pour_points = LayerWidget()

    assert build_dimensions("5.13", dialog)["max_domains"] == 2
    assert build_dimensions("6.0", dialog)["n_domains"] == 2
    v5_values = build_initial_values("5.13", dialog)
    v6_values = build_initial_values("6.0", dialog)
    assert len(v5_values["main"]["directories_general"]["dir_Morpho"]) == 2
    assert len(v6_values["main"]["config_input"]["dem_path"]) == 2


def test_geoparameter_orientation_is_version_specific(tmp_path: Path) -> None:
    metadata = tmp_path / "Z Temp" / "Geometry" / "geology_class_metadata.json"
    metadata.parent.mkdir(parents=True)
    metadata.write_text(
        json.dumps(
            {
                "classes": [
                    {"geo_param": 1, "parameter_value": 42.0},
                    {"geo_param": 2, "parameter_value": 84.0},
                ]
            }
        ),
        encoding="utf-8",
    )
    dialog = Dialog(tmp_path)

    _, _, v5_document = _document("5.13", dialog)
    v5_geo = v5_document["file_profiles"]["parameter"]["values"][
        "geoparameter"
    ]["GeoParam"]
    assert v5_geo == [
        [1.0, 1000.0, 42.0, 1.0, 1.0],
        [1.0, 1000.0, 84.0, 1.0, 1.0],
    ]
    rendered_v5 = json_to_namelist(
        {"values": {"geoparameter": {"GeoParam": v5_geo}}}
    )
    assert "GeoParam(1,3) = 42.0" in rendered_v5
    assert "GeoParam(2,3) = 84.0" in rendered_v5

    _, _, v6_document = _document("6.0", dialog)
    v6_geo = v6_document["file_profiles"]["parameter"]["values"][
        "geoparameter"
    ]["GeoParam"]
    assert v6_geo[2] == [42.0, 84.0]
    rendered_v6 = json_to_namelist(
        {"values": {"geoparameter": {"GeoParam": v6_geo}}}
    )
    assert "GeoParam(3,1) = 42.0" in rendered_v6
    assert "GeoParam(3,2) = 84.0" in rendered_v6


def test_every_profile_form_constructs_and_renders(tmp_path: Path) -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from nml_tools.gui.app import ProfileDialog
    from qtpy.QtWidgets import QApplication

    application = QApplication.instance() or QApplication([])
    dialog = Dialog(tmp_path)
    for version in ("5.13", "6.0"):
        project, dimensions, document = _document(version, dialog)
        for profile in project.profiles:
            editor = ProfileDialog(
                project,
                profile,
                profile_values(document, profile),
                dimensions,
            )
            values = {
                page.name: editor.forms[page.name].values()
                for page in profile.pages
            }
            rendered = render_profile(project, profile, values, dimensions)
            assert rendered
            if version == "5.13" and profile.name == "output":
                assert "&NLoutputResults\n" in rendered
                assert "outputFlxState(1)" in rendered
            if version == "5.13" and profile.name == "parameter":
                assert "&PET0\n" in rendered
                assert "&PETminus1\n" in rendered
            editor.close()
    assert application is not None
