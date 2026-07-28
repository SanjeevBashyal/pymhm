"""Smoke tests for the plain input-combo integration."""

import json
import os
from pathlib import Path

import numpy as np
import xarray as xr

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

# isort: off
from qgis.core import QgsCoordinateReferenceSystem  # noqa: E402
from qgis.PyQt.QtWidgets import QApplication  # noqa: E402

from pymhm.input_selection import InputComboAdapter  # noqa: E402
from pymhm.project_layout import workspace_folder  # noqa: E402
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
    excluded = tmp_path / "mhm-plugin" / "ignored.tif"
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
    assert all("mhm-plugin" not in label for label in labels)
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


def test_meteo_folder_source_and_multiplier_state_round_trip(tmp_path):
    _app()
    precipitation = tmp_path / "forcing" / "pre"
    temperature = tmp_path / "forcing" / "temp"
    precipitation.mkdir(parents=True)
    temperature.mkdir()

    dialog = pymhmDialog()
    assert dialog.comboBox_precipitationFile.count() == 0
    dialog.project_folder = str(tmp_path)
    dialog.refresh_meteo_folder_sources()
    dialog._loading_input_state = True
    try:
        dialog.comboBox_precipitationFile.setCurrentIndex(
            dialog._folder_combo_index(
                dialog.comboBox_precipitationFile,
                str(precipitation),
            )
        )
        dialog.comboBox_temperatureFile.setCurrentIndex(
            dialog._folder_combo_index(
                dialog.comboBox_temperatureFile,
                str(temperature),
            )
        )
        dialog.comboBox_precipitationDataSource.setCurrentIndex(0)
        dialog.comboBox_temperatureDataSource.setCurrentIndex(1)
        dialog.spinBox_L2ResolutionMultiplier.setValue(4)
    finally:
        dialog._loading_input_state = False
    dialog.save_input_state()
    dialog.close()

    state_path = Path(workspace_folder(tmp_path)) / "pymhm_input_state.json"
    state = json.loads(state_path.read_text(encoding="utf-8"))
    assert state["meteo_inputs"]["precipitation"] == {
        "folder": "forcing/pre",
        "source": "mHM ready",
    }
    assert state["meteo_inputs"]["l2_multiplier"] == 4

    restored = pymhmDialog()
    restored.project_folder = str(tmp_path)
    restored.load_input_state()

    assert restored.selected_meteo_folder("precipitation") == str(
        precipitation.resolve()
    )
    assert restored.selected_meteo_source("precipitation") == "mhm_ready"
    assert restored.selected_meteo_source("temperature") == "era5land"
    assert restored.spinBox_L2ResolutionMultiplier.value() == 4
    restored.close()


def test_combined_workflow_green_state_restores_from_workspace(tmp_path):
    _app()
    dialog = pymhmDialog()
    dialog.project_folder = str(tmp_path)
    dialog.morphology_processor.load_processing_state()
    dialog.morphology_processor.mark_workflow_status(
        "meteo_morph_setup",
        "completed",
    )
    dialog.close()

    restored = pymhmDialog()
    restored.project_folder = str(tmp_path)
    restored.morphology_processor.load_processing_state()
    restored.refresh_morphology_workflow_button_states()

    assert "#2e7d32" in restored.pushButton_executeMeteoMorphSetup.styleSheet()
    restored.close()


def test_precipitation_selection_displays_resolution_and_l0_multiplier(tmp_path):
    _app()
    folder = tmp_path / "forcing" / "pre"
    folder.mkdir(parents=True)
    dataset = xr.Dataset(
        {
            "pre": (
                ("time", "lat", "lon"),
                np.ones((1, 2, 2), dtype="float64"),
            )
        },
        coords={
            "time": np.asarray(["2001-01-01"], dtype="datetime64[ns]"),
            "lat": [51.0, 50.0],
            "lon": [10.0, 11.0],
        },
    )
    dataset.to_netcdf(folder / "pre.nc", engine="scipy")

    dialog = pymhmDialog()
    dialog.project_folder = str(tmp_path)
    dialog.mProjectionSelectionWidget_crs.setCrs(
        QgsCoordinateReferenceSystem("EPSG:4326")
    )
    dialog._grid_l0_info = {"resolution": 0.5, "unit": "deg"}
    dialog.refresh_meteo_folder_sources()
    dialog.comboBox_precipitationFile.setCurrentIndex(
        dialog._folder_combo_index(
            dialog.comboBox_precipitationFile,
            str(folder),
        )
    )
    dialog.comboBox_precipitationDataSource.setCurrentIndex(0)

    assert dialog.label_precipitationResolutionValue.text() == "1"
    assert dialog.label_precipitationResolutionUnit.text() == "deg"
    assert dialog.label_precipitationResolutionMultiplier.text() == "2.0"
    dialog.close()
