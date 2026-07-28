"""Smoke tests for the plain input-combo integration."""

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

# isort: off
from qgis.PyQt.QtWidgets import QApplication  # noqa: E402

from pymhm.input_selection import InputComboAdapter  # noqa: E402
from pymhm.pymhm_dialog import pymhmDialog  # noqa: E402
# isort: on

_APPLICATION = None


def _app():
    global _APPLICATION
    _APPLICATION = QApplication.instance() or _APPLICATION or QApplication([])
    return _APPLICATION


def test_dialog_maps_new_combos_and_discovers_relative_project_files(tmp_path):
    _app()
    raster = tmp_path / "0 GIS" / "dem.tif"
    raster.parent.mkdir()
    raster.touch()
    excluded = tmp_path / "Z Temp" / "ignored.tif"
    excluded.parent.mkdir()
    excluded.touch()

    dialog = pymhmDialog()
    dialog.project_folder = str(tmp_path)
    dialog.checkBox_enableFolderSearch.setChecked(True)
    dialog.refresh_input_sources()

    assert isinstance(dialog.mMapLayerComboBox_dem, InputComboAdapter)
    assert dialog.mMapLayerComboBox_dem.combo_box is dialog.comboBox_demInput
    labels = [
        dialog.comboBox_demInput.itemText(index)
        for index in range(dialog.comboBox_demInput.count())
    ]
    assert "0 GIS/dem.tif" in labels
    assert all("Z Temp" not in label for label in labels)
    dialog.close()


def test_project_without_saved_state_clears_previous_project_inputs(tmp_path):
    _app()
    first_project = tmp_path / "first"
    second_project = tmp_path / "second"
    first_project.mkdir()
    second_project.mkdir()
    raster = first_project / "soil.tif"
    raster.touch()

    dialog = pymhmDialog()
    dialog.project_folder = str(first_project)
    dialog.comboBox_soilInput.addItem(
        "soil.tif",
        {
            "origin": "file",
            "kind": "soil",
            "path": str(raster),
            "manual": True,
        },
    )
    dialog.comboBox_soilInput.setCurrentIndex(0)
    dialog.comboBox_soil_inputType.setCurrentIndex(1)
    dialog._categorical_lookup_configs["soil"] = {
        "lookup_table": str(first_project / "soil.csv"),
        "mapping_field": "source",
        "class_field": "class",
    }
    dialog.checkBox_enableFolderSearch.setChecked(True)

    dialog.project_folder = str(second_project)
    dialog.load_input_state()

    assert dialog.comboBox_soilInput.currentIndex() == -1
    assert dialog.comboBox_soil_inputType.currentIndex() == -1
    assert dialog.categorical_lookup_config("soil") is None
    assert not dialog.checkBox_enableFolderSearch.isChecked()
    dialog.close()
