from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

from qgis.PyQt.QtWidgets import QApplication, QPushButton  # noqa: E402

from pymhm.pymhm_dialog import pymhmDialog  # noqa: E402
from pymhm.morphology_task_bridge import _saved_categorical_outputs  # noqa: E402
from pymhm.nml_settings import update_section  # noqa: E402
from pymhm.project_layout import ensure_project_structure  # noqa: E402


_APPLICATION = None


def _app():
    global _APPLICATION
    _APPLICATION = QApplication.instance() or _APPLICATION or QApplication([])
    return _APPLICATION


def test_land_cover_action_uses_path_task_bridge():
    _app()
    dialog = pymhmDialog()
    button = QPushButton()
    calls = []
    dialog.morphology_tasks.start_categorical = (
        lambda kind, **options: calls.append((kind, options)) or True
    )

    dialog.run_background_processor_action(
        "Land Use", lambda: (_ for _ in ()).throw(AssertionError()), button
    )

    assert calls[0][0] == "lc"
    assert calls[0][1]["controls"] == (button,)
    dialog.close()


def test_execute_all_uses_managed_pipeline():
    _app()
    dialog = pymhmDialog()
    calls = []
    dialog.morphology_tasks.start_execute_all = lambda: calls.append(True) or True

    dialog.start_execute_all_processing()

    assert calls == [True]
    dialog.close()


def test_saved_land_cover_and_soil_outputs_are_reusable(tmp_path):
    ensure_project_structure(tmp_path, "5.13")
    morph = tmp_path / "mhm-plugin/data/master/static/morph"
    land_cover = morph / "lc_2000_2010.asc"
    soil = morph / "soil_class.asc"
    definition = morph / "soil_classdefinition.txt"
    for path in (land_cover, soil, definition):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    update_section(
        tmp_path,
        "land_cover",
        {"scenes": [{"output_path": "data/master/static/morph/lc_2000_2010.asc"}]},
    )
    update_section(
        tmp_path,
        "soil",
        {
            "output_path": "data/master/static/morph/soil_class.asc",
            "classdefinition_path": (
                "data/master/static/morph/soil_classdefinition.txt"
            ),
        },
    )

    assert _saved_categorical_outputs(tmp_path, "lc") == (land_cover,)
    assert _saved_categorical_outputs(tmp_path, "soil") == (soil, definition)


def test_execute_all_stage_order_includes_lai_before_hydrology():
    """The button runs the task bridge, so the LAI stage must live there."""
    import inspect

    from pymhm.morphology_task_bridge import MorphologyTaskBridge

    source = inspect.getsource(MorphologyTaskBridge.start_execute_all)
    order = [
        line.split("self.start_", 1)[1].split("(", 1)[0]
        for line in source.splitlines()
        if "lambda next_step: self.start_" in line
    ]
    assert order == [
        "dem_derivatives", "categorical", "categorical", "categorical",
        "lai", "hydrology",
    ]


def test_lai_stage_is_skipped_without_a_netcdf_selection():
    _app()
    dialog = pymhmDialog()
    bridge = dialog.morphology_tasks
    submitted = []
    bridge.coordinator.submit = lambda *a, **k: submitted.append(a) or True
    dialog.morphology_processor._is_lai_long_term_monthly_netcdf_selected = (
        lambda: False
    )
    done = []

    assert bridge.start_lai(done=done.append) is True
    assert submitted == []          # nothing queued
    assert done == [None]           # the pipeline still advances
    dialog.close()


def test_lai_stage_submits_a_background_task_when_selected(tmp_path):
    _app()
    dialog = pymhmDialog()
    dialog.project_folder = str(tmp_path)
    bridge = dialog.morphology_tasks
    submitted = []
    bridge.coordinator.submit = (
        lambda key, label, fn, **k: submitted.append((key, label, k)) or True
    )
    processor = dialog.morphology_processor
    processor._is_lai_long_term_monthly_netcdf_selected = lambda: True
    processor.lai_task_options = lambda output_path: {
        "source_path": "/tmp/lai.nc",
        "source_variable": None,
        "output_path": str(tmp_path / "staged" / "lai_dem.nc"),
        "filled_dem": "/tmp/dem.tif",
        "crs_string": "EPSG:4326",
        "input_resolution": "Biweekly",
        "target_timestep": "Monthly Gridded Data",
    }

    assert bridge.start_lai() is True
    key, label, options = submitted[0]
    assert key == "execute-all-lai"
    assert "LAI" in label
    # It must run on the file resource and be cancellable like the other stages.
    assert options["resource"] == "morphology-files"
    assert options["task_aware"] is True
    dialog.close()


def test_unchanged_lai_inputs_reuse_the_staged_file(tmp_path):
    """The 153 GiB staged LAI must not be regenerated when nothing changed."""
    from pymhm.state_cache import fingerprint, store_payload

    _app()
    dialog = pymhmDialog()
    dialog.project_folder = str(tmp_path)
    bridge = dialog.morphology_tasks
    submitted = []
    bridge.coordinator.submit = lambda *a, **k: submitted.append(a) or True

    source = tmp_path / "lai_src.nc"
    dem = tmp_path / "dem.tif"
    staged = tmp_path / "staged" / "lai_dem.nc"
    for path in (source, dem):
        path.write_text("x", encoding="utf-8")
    staged.parent.mkdir(parents=True, exist_ok=True)
    staged.write_text("staged", encoding="utf-8")

    options = {
        "source_path": str(source),
        "source_variable": None,
        "output_path": str(staged),
        "filled_dem": str(dem),
        "crs_string": "EPSG:4326",
        "input_resolution": "Biweekly",
        "target_timestep": "Monthly Gridded Data",
    }
    processor = dialog.morphology_processor
    processor._is_lai_long_term_monthly_netcdf_selected = lambda: True
    processor.lai_task_options = lambda _output: dict(options)

    # Record the fingerprint the stage will compute for these inputs.
    digest = fingerprint(
        (options["source_path"], options["filled_dem"]),
        {
            "source_variable": None,
            "input_resolution": "Biweekly",
            "target_timestep": "Monthly Gridded Data",
            "method": "bilinear",
            "output_path": str(staged),
        },
    )
    store_payload(tmp_path, "stages", "lai", digest, {"outputs": [str(staged)]})

    done = []
    assert bridge.start_lai(done=done.append) is True
    assert submitted == []                       # no resample queued

    # Editing an input makes it newer than the output, which invalidates both
    # the fingerprint and the adoption heuristic, so the work is queued again.
    staged_mtime = staged.stat().st_mtime
    os.utime(source, (staged_mtime + 10, staged_mtime + 10))
    assert bridge.start_lai() is True
    assert len(submitted) == 1
    dialog.close()


def test_a_missing_staged_file_forces_the_lai_resample(tmp_path):
    from pymhm.state_cache import fingerprint, store_payload

    _app()
    dialog = pymhmDialog()
    dialog.project_folder = str(tmp_path)
    bridge = dialog.morphology_tasks
    submitted = []
    bridge.coordinator.submit = lambda *a, **k: submitted.append(a) or True

    source = tmp_path / "lai_src.nc"
    dem = tmp_path / "dem.tif"
    for path in (source, dem):
        path.write_text("x", encoding="utf-8")
    staged = tmp_path / "staged" / "lai_dem.nc"      # recorded but deleted

    options = {
        "source_path": str(source), "source_variable": None,
        "output_path": str(staged), "filled_dem": str(dem),
        "crs_string": "EPSG:4326", "input_resolution": "Monthly",
        "target_timestep": "Long Term Mean Monthly Gridded Data",
    }
    processor = dialog.morphology_processor
    processor._is_lai_long_term_monthly_netcdf_selected = lambda: True
    processor.lai_task_options = lambda _output: dict(options)
    digest = fingerprint(
        (options["source_path"], options["filled_dem"]),
        {
            "source_variable": None,
            "input_resolution": "Monthly",
            "target_timestep": "Long Term Mean Monthly Gridded Data",
            "method": "bilinear",
            "output_path": str(staged),
        },
    )
    store_payload(tmp_path, "stages", "lai", digest, {"outputs": [str(staged)]})

    assert bridge.start_lai() is True
    assert len(submitted) == 1                   # deleted output -> rework
    dialog.close()


def test_geology_outputs_already_on_disk_are_adopted_not_rebuilt(tmp_path):
    """A project prepared before the cache existed must not rebuild geology."""
    _app()
    dialog = pymhmDialog()
    dialog.project_folder = str(tmp_path)
    bridge = dialog.morphology_tasks
    submitted = []
    bridge.coordinator.submit = lambda *a, **k: submitted.append(a) or True
    messages = []
    dialog.log_message = messages.append

    inputs = {}
    for name in ("geology_input.shp", "1_dem_filled.tif", "lookup.csv"):
        path = tmp_path / name
        path.write_text("x", encoding="utf-8")
        inputs[name] = path
    outputs = {}
    for name in ("3_geology_processed.tif", "geology_classdefinition.txt",
                 "geology_class_metadata.json"):
        path = tmp_path / name
        path.write_text("y", encoding="utf-8")
        outputs[name] = path
    # Outputs postdate every input.
    newest = max(path.stat().st_mtime for path in inputs.values())
    for path in outputs.values():
        os.utime(path, (newest + 10, newest + 10))

    job = {
        "kind": "geology",
        "input_file": str(inputs["geology_input.shp"]),
        "filled_dem": str(inputs["1_dem_filled.tif"]),
        "lookup_table": str(inputs["lookup.csv"]),
        "output": str(outputs["3_geology_processed.tif"]),
        "definition": str(outputs["geology_classdefinition.txt"]),
        "metadata": str(outputs["geology_class_metadata.json"]),
        "mapping_field": "ID",
        "class_field": "CLASS",
        "is_vector": True,
        "input_crs": "EPSG:4326",
        "dem_crs": "EPSG:4326",
    }
    processor = dialog.morphology_processor
    processor.filled_dem_path = str(inputs["1_dem_filled.tif"])
    processor._categorical_mode = lambda _kind: "lookup_table"
    processor._record_categorical_nml = lambda *_a, **_k: None
    dialog.uses_advanced_categorical_input = lambda _kind: False
    bridge._lookup_job = lambda _kind: dict(job)

    assert bridge.start_categorical("geology", reuse_existing=True) is True
    assert submitted == []                                  # nothing rebuilt
    assert any("Adopted the existing geology output" in m for m in messages)

    # The adoption is recorded, so the next run reuses without re-adopting.
    from pymhm.state_cache import load_state

    assert "geology" in load_state(tmp_path)["stages"]
    messages.clear()
    assert bridge.start_categorical("geology", reuse_existing=True) is True
    assert submitted == []
    assert any("inputs are unchanged" in m for m in messages)
    dialog.close()
