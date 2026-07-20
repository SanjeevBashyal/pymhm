"""Focused tests for the path-based mhm-tools geology adapter."""

from pathlib import Path

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

from qgis.core import QgsRasterLayer, QgsVectorLayer  # noqa: E402

from pymhm.Morphology.layers import geology  # noqa: E402


class _Combo:
    def __init__(self, layer):
        self._layer = layer

    def currentLayer(self):
        return self._layer


class _Dialog:
    def __init__(self, project_folder, geology_layer, lookup_layer):
        self.project_folder = str(project_folder)
        self.mMapLayerComboBox_geology = _Combo(geology_layer)
        self.mMapLayerComboBox_geologyLookup = _Combo(lookup_layer)


class _Processor(geology.GeologyProcessingMixin):
    def __init__(self, dialog, filled_dem):
        self.dialog = dialog
        self.filled_dem_path = str(filled_dem)
        self.messages = []
        self.loaded = []
        self.prepared = []

    def check_prerequisites(self):
        return True

    def _ensure_filled_dem(self, _callback):
        return True

    def _geology_lookup_input(self):
        return self.dialog.mMapLayerComboBox_geologyLookup.currentLayer(), "map_code"

    def _required_lookup_field(self, _field_names, required_name):
        return required_name

    def log_message(self, message):
        self.messages.append(message)

    def load_layer(self, path, name):
        self.loaded.append((path, name))

    def mark_output_prepared(self, path, **_kwargs):
        self.prepared.append(path)

    def geology_class_metadata_writer(self):
        output = (
            Path(self.dialog.project_folder)
            / "Z Temp"
            / "Geometry"
            / "geology_class_metadata.json"
        )
        output.write_text("{}\n", encoding="utf-8")
        return str(output)


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return path


def test_vector_geology_uses_rasterizer_and_definition(tmp_path, monkeypatch):
    geology_file = _touch(tmp_path / "input" / "geology.gpkg")
    lookup_file = _touch(tmp_path / "input" / "lookup.gpkg")
    dem_file = _touch(tmp_path / "project" / "Z Temp" / "Geometry" / "dem.tif")
    processor = _Processor(
        _Dialog(
            tmp_path / "project",
            QgsVectorLayer(str(geology_file), "geology", "ogr"),
            QgsVectorLayer(str(lookup_file), "lookup", "ogr"),
        ),
        dem_file,
    )
    calls = {}

    def fake_rasterize(**kwargs):
        calls["rasterize"] = kwargs
        output = Path(kwargs["output_file"])
        output.write_bytes(b"geology raster")
        return output

    def fake_definition(**kwargs):
        calls["definition"] = kwargs
        output = Path(kwargs["output_file"])
        output.write_text("nGeo_Formations 1\n", encoding="utf-8")
        return output

    monkeypatch.setattr(geology, "rasterize_geology_map", fake_rasterize)
    monkeypatch.setattr(geology, "write_geology_classdefinition_file", fake_definition)

    assert processor.process_geology()
    assert calls["rasterize"]["dem_file"] == str(dem_file)
    assert calls["rasterize"]["lookup_table"] == str(lookup_file)
    definition = (
        tmp_path
        / "project"
        / "data"
        / "static"
        / "morph"
        / "geology_classdefinition.txt"
    )
    assert definition.read_text(encoding="utf-8") == "nGeo_Formations 1\n"


def test_raster_geology_passes_dem_to_formatter(tmp_path, monkeypatch):
    geology_file = _touch(tmp_path / "input" / "geology.tif")
    lookup_file = _touch(tmp_path / "input" / "lookup.gpkg")
    dem_file = _touch(tmp_path / "project" / "Z Temp" / "Geometry" / "dem.tif")
    processor = _Processor(
        _Dialog(
            tmp_path / "project",
            QgsRasterLayer(str(geology_file), "geology"),
            QgsVectorLayer(str(lookup_file), "lookup", "ogr"),
        ),
        dem_file,
    )
    calls = {}

    def fake_format(**kwargs):
        calls.update(kwargs)
        output = Path(kwargs["output_file"])
        output.write_bytes(b"geology raster")
        Path(kwargs["classdefinition_file"]).write_text(
            "nGeo_Formations 1\n", encoding="utf-8"
        )
        return output

    monkeypatch.setattr(geology, "format_geology_file", fake_format)

    assert processor.process_geology(write_classdefinition=False)
    assert calls["input_file"] == str(geology_file)
    assert calls["dem_file"] == str(dem_file)
    assert not (
        tmp_path
        / "project"
        / "data"
        / "static"
        / "morph"
        / "geology_classdefinition.txt"
    ).exists()


def test_geology_fails_when_required_metadata_is_not_written(tmp_path, monkeypatch):
    """A raster and definition alone are not a complete geology result."""
    geology_file = _touch(tmp_path / "input" / "geology.gpkg")
    lookup_file = _touch(tmp_path / "input" / "lookup.gpkg")
    dem_file = _touch(tmp_path / "project" / "Z Temp" / "Geometry" / "dem.tif")
    processor = _Processor(
        _Dialog(
            tmp_path / "project",
            QgsVectorLayer(str(geology_file), "geology", "ogr"),
            QgsVectorLayer(str(lookup_file), "lookup", "ogr"),
        ),
        dem_file,
    )
    old_definition = _touch(
        tmp_path
        / "project"
        / "data"
        / "static"
        / "morph"
        / "geology_classdefinition.txt"
    )
    old_definition.write_text("old definition\n", encoding="utf-8")
    old_metadata = _touch(
        tmp_path / "project" / "Z Temp" / "Geometry" / "geology_class_metadata.json"
    )
    old_metadata.write_text('{"old": true}\n', encoding="utf-8")

    def fake_rasterize(**kwargs):
        output = Path(kwargs["output_file"])
        output.write_bytes(b"geology raster")
        return output

    def fake_definition(**kwargs):
        output = Path(kwargs["output_file"])
        output.write_text("nGeo_Formations 1\n", encoding="utf-8")
        return output

    monkeypatch.setattr(geology, "rasterize_geology_map", fake_rasterize)
    monkeypatch.setattr(geology, "write_geology_classdefinition_file", fake_definition)
    monkeypatch.setattr(geology.QMessageBox, "critical", lambda *_args: None)
    processor.geology_class_metadata_writer = lambda: None

    assert not processor.process_geology()
    assert not (
        tmp_path / "project" / "Z Temp" / "Geometry" / "3_geology_processed.tif"
    ).exists()
    assert old_definition.read_text(encoding="utf-8") == "old definition\n"
    assert old_metadata.read_text(encoding="utf-8") == '{"old": true}\n'
