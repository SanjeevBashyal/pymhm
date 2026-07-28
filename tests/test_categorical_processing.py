"""Focused tests for pymhm categorical workflow orchestration."""

from pathlib import Path

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

# isort: off
from qgis.core import (  # noqa: E402
    QgsCoordinateReferenceSystem,
    QgsRasterLayer,
    QgsVectorLayer,
)

from pymhm.Morphology.layers import categorical  # noqa: E402
# isort: on


class _Combo:
    def __init__(self, layer):
        self.layer = layer

    def currentLayer(self):
        return self.layer

    def source_path(self):
        return self.layer.source()


class _Dialog:
    def __init__(self, project, kind, layer, mode, lookup=None):
        self.project_folder = str(project)
        setattr(self, categorical._SPECS[kind]["combo"], _Combo(layer))
        self.mode = mode
        self.lookup = lookup

    def categorical_input_mode(self, _kind):
        return self.mode

    def categorical_lookup_config(self, _kind):
        return self.lookup

    def get_crs(self):
        return QgsCoordinateReferenceSystem("EPSG:32645")


class _Processor(categorical.CategoricalProcessingMixin):
    def __init__(self, dialog, dem):
        self.dialog = dialog
        self.filled_dem_path = str(dem)
        self.messages = []
        self.loaded = []
        self.prepared = []
        self.land_use_layer = None
        self.geology_path = None
        self.categorical_ready_outputs = {}
        self.prerequisite_calls = 0

    def check_prerequisites(self):
        self.prerequisite_calls += 1
        return True

    def _ensure_filled_dem(self, _callback):
        return True

    def log_message(self, message):
        self.messages.append(message)

    def load_layer(self, path, name):
        self.loaded.append((path, name))

    def mark_output_prepared(self, path, **_kwargs):
        self.prepared.append(path)


def _touch(path, content=b""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


def _lookup(path):
    return {
        "lookup_table": str(_touch(path)),
        "mapping_field": "source_code",
        "class_field": "class_id",
    }


def test_raster_land_cover_delegates_to_unified_adapter(tmp_path, monkeypatch):
    source = _touch(tmp_path / "input" / "land_cover.tif")
    geometry = tmp_path / "project" / "Z Temp" / "Geometry"
    dem = _touch(geometry / "1_dem_filled.tif")
    stale_crop = _touch(geometry / "3_land_use_crop.tif", b"old")
    stale_mask = _touch(geometry / "3_land_use_masked.tif", b"old")
    layer = QgsRasterLayer(str(source), "land cover")
    config = _lookup(tmp_path / "input" / "lookup.csv")
    processor = _Processor(
        _Dialog(tmp_path / "project", "lc", layer, "Lookup Table", config), dem
    )
    calls = {}

    def prepare(*args, **kwargs):
        calls["args"] = args
        calls["kwargs"] = kwargs
        _touch(Path(args[3]), b"lc")

    monkeypatch.setattr(categorical, "prepare_categorical_file", prepare)

    assert processor.process_land_use()
    assert calls["args"][0] == "lc"
    assert calls["args"][1] == str(source)
    assert calls["args"][4:7] == (
        config["lookup_table"],
        "source_code",
        "class_id",
    )
    assert calls["kwargs"]["is_vector"] is False
    assert Path(processor.land_use_layer).name == "3_land_use.tif"
    assert not stale_crop.exists()
    assert not stale_mask.exists()


def test_vector_soil_passes_vector_mode_and_definition_target(tmp_path, monkeypatch):
    source = _touch(tmp_path / "input" / "soil.shp")
    dem = _touch(tmp_path / "project" / "Z Temp" / "Geometry" / "1_dem_filled.tif")
    layer = QgsVectorLayer(str(source), "soil", "ogr")
    config = _lookup(tmp_path / "input" / "lookup.csv")
    processor = _Processor(
        _Dialog(tmp_path / "project", "soil", layer, "Lookup Table", config), dem
    )
    calls = {}

    def prepare(*args, **kwargs):
        calls["args"] = args
        calls["kwargs"] = kwargs
        _touch(Path(args[3]), b"soil")
        _touch(Path(kwargs["classdefinition_file"]), b"definition")

    monkeypatch.setattr(categorical, "prepare_categorical_file", prepare)

    assert processor.process_soil()
    assert calls["kwargs"]["is_vector"] is True
    assert Path(calls["kwargs"]["classdefinition_file"]).name == (
        "soil_classdefinition.txt"
    )
    assert Path(calls["args"][3]).name == "3_soil.tif"


def test_ready_soil_copies_raster_and_removes_stale_intermediates(
    tmp_path, monkeypatch
):
    project = tmp_path / "project"
    source = _touch(tmp_path / "input" / "soil.asc", b"ready")
    dem = _touch(project / "Z Temp" / "Geometry" / "1_dem_filled.tif")
    stale = _touch(project / "Z Temp" / "Geometry" / "3_soil.tif", b"stale")
    definition = _touch(
        project / "data" / "static" / "morph" / "soil_classdefinition.txt",
        b"definition",
    )
    layer = QgsRasterLayer(str(source), "soil")
    processor = _Processor(_Dialog(project, "soil", layer, "mHM_ready"), dem)
    monkeypatch.setattr(
        categorical,
        "prepare_categorical_file",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("ready data must not be formatted")
        ),
    )

    assert processor.process_soil()
    target = project / "data" / "static" / "morph" / "soil_class.asc"
    assert target.read_bytes() == b"ready"
    assert definition.is_file()
    assert not stale.exists()
    assert processor.categorical_ready_outputs["soil"] == str(target)
    assert processor.prerequisite_calls == 0


def test_ready_mode_rejects_vector_input(tmp_path, monkeypatch):
    source = _touch(tmp_path / "input" / "soil.shp")
    dem = _touch(tmp_path / "project" / "Z Temp" / "Geometry" / "dem.tif")
    processor = _Processor(
        _Dialog(
            tmp_path / "project",
            "soil",
            QgsVectorLayer(str(source), "soil", "ogr"),
            "mHM_ready",
        ),
        dem,
    )
    monkeypatch.setattr(categorical.QMessageBox, "warning", lambda *_args: None)

    assert not processor.process_soil()


def test_ready_land_cover_is_available_to_elevation_bands(tmp_path):
    from pymhm.Morphology.elevation_bands.band_landcover_helpers import \
        BandLandCoverHelperMixin

    ready = _touch(tmp_path / "project" / "data" / "static" / "morph" / "lc.asc")
    processor = _Processor.__new__(_Processor)
    processor.dialog = type("Dialog", (), {"project_folder": str(tmp_path / "project")})()
    processor.land_use_layer = None
    processor.categorical_ready_outputs = {"lc": str(ready)}

    result = BandLandCoverHelperMixin._ensure_land_use_raster(
        processor,
        lambda: (_ for _ in ()).throw(
            AssertionError("ready land cover must not be processed again")
        ),
    )

    assert result == str(ready)


def test_lookup_replaces_ready_state_and_files(tmp_path, monkeypatch):
    project = tmp_path / "project"
    source = _touch(tmp_path / "input" / "land_cover.tif")
    dem = _touch(project / "Z Temp" / "Geometry" / "1_dem_filled.tif")
    config = _lookup(tmp_path / "input" / "lookup.csv")
    processor = _Processor(
        _Dialog(
            project,
            "lc",
            QgsRasterLayer(str(source), "land cover"),
            "Lookup Table",
            config,
        ),
        dem,
    )
    morph = project / "data" / "static" / "morph"
    stale = [_touch(morph / f"lc{suffix}", b"stale") for suffix in (".asc", ".nc")]
    processor.categorical_ready_outputs["lc"] = str(stale[0])

    def prepare(*args, **_kwargs):
        _touch(Path(args[3]), b"lookup")

    monkeypatch.setattr(categorical, "prepare_categorical_file", prepare)

    assert processor.process_land_use()
    assert all(not path.exists() for path in stale)
    assert "lc" not in processor.categorical_ready_outputs
    assert (
        project / "Z Temp" / "Geometry" / "3_land_use.tif"
    ).read_bytes() == b"lookup"


def test_failed_geology_publish_restores_previous_outputs(tmp_path, monkeypatch):
    project = tmp_path / "project"
    source = _touch(tmp_path / "input" / "geology.tif")
    dem = _touch(project / "Z Temp" / "Geometry" / "1_dem_filled.tif")
    config = _lookup(tmp_path / "input" / "lookup.csv")
    processor = _Processor(
        _Dialog(
            project,
            "geology",
            QgsRasterLayer(str(source), "geology"),
            "Lookup Table",
            config,
        ),
        dem,
    )
    output = _touch(
        project / "Z Temp" / "Geometry" / "3_geology_processed.tif", b"old raster"
    )
    definition = _touch(
        project / "data" / "static" / "morph" / "geology_classdefinition.txt",
        b"old definition",
    )
    metadata = _touch(
        project / "Z Temp" / "Geometry" / "geology_class_metadata.json",
        b"old metadata",
    )

    def prepare(*args, **kwargs):
        _touch(Path(args[3]), b"new raster")
        _touch(Path(kwargs["classdefinition_file"]), b"new definition")

    def write_metadata(_lookup, _field, path):
        return _touch(Path(path), b"new metadata")

    real_replace = categorical.os.replace
    failed = False

    def fail_definition_publish(source_path, destination_path):
        nonlocal failed
        if Path(destination_path) == definition and not failed:
            failed = True
            raise OSError("simulated publish failure")
        return real_replace(source_path, destination_path)

    monkeypatch.setattr(categorical, "prepare_categorical_file", prepare)
    monkeypatch.setattr(categorical, "write_geology_metadata", write_metadata)
    monkeypatch.setattr(categorical.os, "replace", fail_definition_publish)
    monkeypatch.setattr(categorical.QMessageBox, "critical", lambda *_args: None)

    assert not processor.process_geology()
    assert output.read_bytes() == b"old raster"
    assert definition.read_bytes() == b"old definition"
    assert metadata.read_bytes() == b"old metadata"


def test_failed_backup_does_not_delete_existing_output(tmp_path, monkeypatch):
    project = tmp_path / "project"
    source = _touch(tmp_path / "input" / "land_cover.tif")
    dem = _touch(project / "Z Temp" / "Geometry" / "1_dem_filled.tif")
    config = _lookup(tmp_path / "input" / "lookup.csv")
    processor = _Processor(
        _Dialog(
            project,
            "lc",
            QgsRasterLayer(str(source), "land cover"),
            "Lookup Table",
            config,
        ),
        dem,
    )
    output = _touch(
        project / "Z Temp" / "Geometry" / "3_land_use.tif",
        b"existing",
    )

    def prepare(*args, **_kwargs):
        _touch(Path(args[3]), b"new")

    real_replace = categorical.os.replace

    def fail_backup(source_path, destination_path):
        if Path(source_path) == output:
            raise OSError("simulated backup failure")
        return real_replace(source_path, destination_path)

    monkeypatch.setattr(categorical, "prepare_categorical_file", prepare)
    monkeypatch.setattr(categorical.os, "replace", fail_backup)
    monkeypatch.setattr(categorical.QMessageBox, "critical", lambda *_args: None)

    assert not processor.process_land_use()
    assert output.read_bytes() == b"existing"
