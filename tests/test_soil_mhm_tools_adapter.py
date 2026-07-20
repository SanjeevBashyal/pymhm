"""Focused tests for the path-based mhm-tools soil adapter."""

from pathlib import Path

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

from qgis.core import QgsRasterLayer, QgsVectorLayer  # noqa: E402

from pymhm.Morphology.layers import soil_raster  # noqa: E402
from pymhm.Morphology.core.soil_sources import local_layer_source  # noqa: E402


class _Combo:
    def __init__(self, layer):
        self._layer = layer

    def currentLayer(self):
        return self._layer


class _Dialog:
    def __init__(self, project_folder, soil_layer, lookup_layer):
        self.project_folder = str(project_folder)
        self.mMapLayerComboBox_soil = _Combo(soil_layer)
        self.mMapLayerComboBox_soilLookup = _Combo(lookup_layer)


class _Processor(soil_raster.SoilRasterMixin):
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

    def load_soil_lookup_table(self):
        self.soil_lookup_key_field = "map_code"
        return {"1": 1}

    def _required_lookup_field(self, _field_names, required_name):
        return required_name

    def log_message(self, message):
        self.messages.append(message)

    def load_layer(self, path, name):
        self.loaded.append((path, name))

    def mark_output_prepared(self, path, **_kwargs):
        self.prepared.append(path)


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return path


def test_local_layer_source_excludes_sublayer_uris(tmp_path):
    """Plain paths are reused, while selected provider sublayers are copied."""
    source = _touch(tmp_path / "soil.gpkg")
    plain_layer = QgsVectorLayer(str(source), "soil", "ogr")

    class _UriLayer:
        def source(self):
            return f"{source}|layername=selected_soil"

    assert local_layer_source(plain_layer) == str(source)
    assert local_layer_source(_UriLayer()) is None


def test_vector_process_uses_direct_paths_and_writes_both_outputs(
    tmp_path, monkeypatch
):
    """The vector UI path delegates classification and definition to mhm-tools."""
    soil_file = _touch(tmp_path / "input" / "soil.gpkg")
    lookup_file = _touch(tmp_path / "input" / "lookup.gpkg")
    dem_file = _touch(tmp_path / "project" / "Z Temp" / "Geometry" / "dem.tif")
    processor = _Processor(
        _Dialog(
            tmp_path / "project",
            QgsVectorLayer(str(soil_file), "soil", "ogr"),
            QgsVectorLayer(str(lookup_file), "lookup", "ogr"),
        ),
        dem_file,
    )
    calls = {}

    def fake_rasterize_soil_map(**kwargs):
        calls["rasterize"] = kwargs
        output = Path(kwargs["output_file"])
        output.write_bytes(b"soil raster")
        return output

    def fake_write_definition(**kwargs):
        calls["definition"] = kwargs
        output = Path(kwargs["output_file"])
        output.write_text("nSoil_Types 1\n", encoding="utf-8")
        return output

    monkeypatch.setattr(soil_raster, "rasterize_soil_map", fake_rasterize_soil_map)
    monkeypatch.setattr(
        soil_raster, "write_soil_classdefinition_file", fake_write_definition
    )
    monkeypatch.setattr(
        soil_raster,
        "materialize_vector_layer",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("plain local layers should not be materialized")
        ),
    )

    assert processor.process_soil()
    assert calls["rasterize"]["input_file"] == str(soil_file)
    assert calls["rasterize"]["dem_file"] == str(dem_file)
    assert calls["rasterize"]["lookup_table"] == str(lookup_file)
    assert Path(calls["rasterize"]["output_file"]).name == "3_soil.tif"
    definition = (
        tmp_path / "project" / "data" / "static" / "morph" / "soil_classdefinition.txt"
    )
    assert definition.read_text(encoding="utf-8") == "nSoil_Types 1\n"


def test_raster_process_discards_definition_when_legacy_flag_is_false(
    tmp_path, monkeypatch
):
    """The mandatory core definition output stays temporary when requested."""
    soil_file = _touch(tmp_path / "input" / "soil.tif")
    lookup_file = _touch(tmp_path / "input" / "lookup.gpkg")
    dem_file = _touch(tmp_path / "project" / "Z Temp" / "Geometry" / "dem.tif")
    processor = _Processor(
        _Dialog(
            tmp_path / "project",
            QgsRasterLayer(str(soil_file), "soil"),
            QgsVectorLayer(str(lookup_file), "lookup", "ogr"),
        ),
        dem_file,
    )
    calls = {}

    def fake_format_soil_file(**kwargs):
        calls.update(kwargs)
        output = Path(kwargs["output_file"])
        output.write_bytes(b"soil raster")
        Path(kwargs["classdefinition_file"]).write_text(
            "nSoil_Types 1\n", encoding="utf-8"
        )
        return output

    monkeypatch.setattr(soil_raster, "format_soil_file", fake_format_soil_file)
    monkeypatch.setattr(
        soil_raster,
        "materialize_vector_layer",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("plain local layers should not be materialized")
        ),
    )

    assert processor.process_soil(write_classdefinition=False)
    assert calls["input_file"] == str(soil_file)
    assert calls["dem_file"] == str(dem_file)
    assert calls["lookup_table"] == str(lookup_file)
    final_definition = (
        tmp_path / "project" / "data" / "static" / "morph" / "soil_classdefinition.txt"
    )
    assert not final_definition.exists()
    assert not Path(calls["classdefinition_file"]).exists()
