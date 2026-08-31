"""Soil, geology, and land-cover preparation through :mod:`mhm_tools`."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

from ...applications.mhm_tools_handler import prepare_categorical_file
from ..common import (QgsRasterLayer, QgsVectorLayer, QMessageBox,
                      morph_folder, morph_staging_folder,
                      project_geometry_folder)
from ..core.layer_preparation import LayerPreparationMixin
from ..core.predecessors import PredecessorMixin
from ..core.soil_sources import local_layer_source, materialize_vector_layer
from ..watershed.dem_fill import DemFillMixin
from ...geology_metadata import write_geology_metadata

_SPECS = {
    "lc": {
        "kind": "land_cover",
        "label": "Land Cover",
        "geometry": "3_land_use.tif",
        "layer_name": "3_Land_Use",
        "ready_name": "lc",
    },
    "soil": {
        "kind": "soil",
        "label": "Soil",
        "geometry": "3_soil.tif",
        "layer_name": "3_Soil",
        "ready_name": "soil_class",
        "definition": "soil_classdefinition.txt",
    },
    "geology": {
        "kind": "geology",
        "label": "Geology",
        "geometry": "3_geology_processed.tif",
        "layer_name": "3_Geology",
        "ready_name": "geology_class",
        "definition": "geology_classdefinition.txt",
    },
}


class CategoricalProcessingMixin(
    LayerPreparationMixin,
    DemFillMixin,
    PredecessorMixin,
):
    """Prepare categorical morphology inputs with one shared workflow."""

    def process_land_use(self) -> bool:
        uses_advanced = getattr(self.dialog, "uses_advanced_categorical_input", None)
        advanced = getattr(self.dialog, "process_advanced_categorical_input", None)
        if uses_advanced is not None and advanced is not None and uses_advanced("lc"):
            return bool(advanced("lc"))
        return self._process_categorical("lc")

    def process_soil(self, write_classdefinition=True) -> bool:
        uses_advanced = getattr(self.dialog, "uses_advanced_categorical_input", None)
        advanced = getattr(self.dialog, "process_advanced_categorical_input", None)
        if uses_advanced is not None and advanced is not None and uses_advanced("soil"):
            return bool(advanced("soil"))
        return self._process_categorical(
            "soil", write_classdefinition=write_classdefinition
        )

    def process_geology(self, write_classdefinition=True) -> bool:
        return self._process_categorical(
            "geology", write_classdefinition=write_classdefinition
        )

    def _process_categorical(self, kind, write_classdefinition=True):
        spec = _SPECS[kind]
        if not getattr(self.dialog, "project_folder", None):
            QMessageBox.critical(
                self.dialog,
                "Missing Input",
                "Please select a project folder before proceeding.",
            )
            return False

        layer, source, is_vector = self._categorical_input(kind, spec)
        if not source:
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                f"Please select a {spec['label']} input.",
            )
            return False

        mode = self._categorical_mode(kind)
        if mode == "mhm_ready":
            try:
                return self._copy_ready_categorical(kind, source, is_vector)
            except Exception as error:
                return self._categorical_error(spec, error)
        if mode != "lookup_table":
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                f"Please select an input type for {spec['label']}.",
            )
            return False
        if not self.check_prerequisites():
            return False
        if not self._ensure_filled_dem(self.fill_dem):
            return False

        config = self._categorical_lookup(kind)
        required = ("lookup_table", "mapping_field", "class_field")
        if not config or any(not config.get(key) for key in required):
            QMessageBox.warning(
                self.dialog,
                "Lookup Table Required",
                f"Configure the {spec['label']} lookup table and fields first.",
            )
            return False

        geometry = Path(project_geometry_folder(self.dialog.project_folder))
        staging = Path(morph_staging_folder(self.dialog.project_folder))
        geometry.mkdir(parents=True, exist_ok=True)
        staging.mkdir(parents=True, exist_ok=True)
        output = geometry / spec["geometry"]
        definition = (
            staging / spec["definition"]
            if write_classdefinition and spec.get("definition")
            else None
        )
        metadata = geometry / "geology_class_metadata.json"

        try:
            with TemporaryDirectory(prefix=f"mhm_qgis_{kind}_", dir=geometry) as temp:
                temporary = Path(temp)
                input_file = self._categorical_input_file(
                    layer, source, is_vector, temporary
                )
                prepared_output = temporary / spec["geometry"]
                prepared_definition = (
                    temporary / definition.name if definition is not None else None
                )
                prepared_metadata = (
                    temporary / metadata.name
                    if kind == "geology" and write_classdefinition
                    else None
                )
                if prepared_metadata is not None:
                    write_geology_metadata(
                        config["lookup_table"],
                        config["class_field"],
                        prepared_metadata,
                    )
                prepare_categorical_file(
                    kind,
                    input_file,
                    self.filled_dem_path,
                    prepared_output,
                    config["lookup_table"],
                    config["mapping_field"],
                    config["class_field"],
                    is_vector=is_vector,
                    classdefinition_file=prepared_definition,
                    input_crs=self._input_crs(layer),
                    dem_crs=self._dem_crs(),
                    log=self.log_message,
                )
                replacements = [(prepared_output, output)]
                if prepared_definition is not None:
                    replacements.append((prepared_definition, definition))
                if prepared_metadata is not None:
                    replacements.append((prepared_metadata, metadata))
                removals = self._ready_output_paths(kind)
                removals.extend(self._categorical_intermediate_paths(kind)[1:])
                self._publish_categorical_outputs(
                    replacements,
                    removals,
                    temporary,
                )
                self.categorical_ready_outputs.pop(kind, None)
            if prepared_metadata is not None:
                self.mark_output_prepared(
                    str(metadata),
                    name=metadata.name,
                    loaded=False,
                )
        except Exception as error:
            return self._categorical_error(spec, error)

        if kind == "lc":
            self.land_use_layer = str(output)
        elif kind == "geology":
            self.geology_path = str(output)
        self.load_layer(str(output), spec["layer_name"])
        self.mark_output_prepared(str(output), name=output.name)
        if definition:
            self.mark_output_prepared(
                str(definition), name=definition.name, loaded=False
            )
        if kind == "soil":
            record = getattr(self.dialog, "record_standard_soil_output", None)
            if record is not None:
                record(None, str(definition) if definition else None)
        self._record_categorical_nml(kind, output, definition)
        self.log_message(f"{spec['label']} data prepared successfully.")
        return True

    def _categorical_input(self, kind, spec):
        combo = self.dialog.input_combo(spec["kind"])
        layer = combo.currentLayer() if combo is not None else None
        source = ""
        if combo is not None and hasattr(combo, "source_path"):
            source = combo.source_path() or ""
        if not source and layer is not None:
            source = str(layer.source() or "")
        if not source:
            source_config = getattr(
                self.dialog, "categorical_source_config", lambda _kind: None
            )(kind)
            if source_config:
                source = str(source_config.get("input_path", "") or "")
        suffix = Path(source.split("|", 1)[0]).suffix.lower()
        return layer, source, isinstance(layer, QgsVectorLayer) or suffix == ".shp"

    def _categorical_input_file(self, layer, source, is_vector, temporary):
        local = local_layer_source(layer) if layer is not None else None
        if local:
            return local
        if is_vector and layer is not None:
            return materialize_vector_layer(layer, str(temporary / "input.gpkg"))
        if layer is not None:
            return self._materialize_categorical_raster(
                layer, str(temporary / "input.tif")
            )
        path = Path(source.split("|", 1)[0])
        if path.is_file():
            return str(path)
        raise ValueError(f"Input is not a readable local file: {source}")

    def _materialize_categorical_raster(self, layer, output):
        result = self.run_processing_algorithm(
            "gdal:translate",
            {
                "INPUT": layer,
                "TARGET_CRS": None,
                "NODATA": None,
                "COPY_SUBDATASETS": False,
                "OPTIONS": None,
                "EXTRA": "",
                "DATA_TYPE": 0,
                "OUTPUT": output,
            },
        )
        if not result or not Path(output).is_file():
            raise RuntimeError(
                f"Could not materialize raster layer '{layer.name()}'."
            )
        return output

    def _copy_ready_categorical(self, kind, source, is_vector):
        spec = _SPECS[kind]
        source_path = Path(source.split("|", 1)[0])
        if is_vector or source_path.suffix.lower() not in {".asc", ".nc", ".tif"}:
            QMessageBox.warning(
                self.dialog,
                "Invalid mHM-ready Input",
                f"{spec['label']} mHM-ready input must be an ASC, NetCDF, "
                "or TIFF raster.",
            )
            return False
        if "|" in source or not source_path.is_file():
            QMessageBox.warning(
                self.dialog,
                "Invalid mHM-ready Input",
                "mHM-ready input must be a directly accessible local raster file.",
            )
            return False

        target = (
            Path(morph_folder(self.dialog.project_folder))
            / f"{spec['ready_name']}{source_path.suffix.lower()}"
        )
        target.parent.mkdir(parents=True, exist_ok=True)
        definition = spec.get("definition")
        source_config = getattr(
            self.dialog, "categorical_source_config", lambda _kind: None
        )(kind) or {}
        definition_source = Path(
            str(source_config.get("classdefinition_path", "") or "")
        )
        if definition and not definition_source.is_file() and not (
            target.parent / definition
        ).is_file():
            QMessageBox.warning(
                self.dialog,
                "Missing Class Definition",
                f"mHM-ready {spec['label']} requires a class-definition file.",
            )
            return False

        with TemporaryDirectory(
            prefix=f"mhm_qgis_{kind}_ready_", dir=target.parent
        ) as temp:
            temporary = Path(temp)
            prepared = temporary / target.name
            shutil.copy2(source_path, prepared)
            replacements = [(prepared, target)]
            if definition and definition_source.is_file():
                prepared_definition = temporary / definition
                shutil.copy2(definition_source, prepared_definition)
                replacements.append((prepared_definition, target.parent / definition))
            removals = self._ready_output_paths(kind)
            removals.extend(self._categorical_intermediate_paths(kind))
            if source_path.suffix.lower() == ".asc":
                for suffix in (".prj", ".PRJ"):
                    sidecar = source_path.with_suffix(suffix)
                    if sidecar.is_file():
                        prepared_sidecar = prepared.with_suffix(suffix)
                        shutil.copy2(sidecar, prepared_sidecar)
                        replacements.append(
                            (prepared_sidecar, target.with_suffix(suffix))
                        )
                        break
            self._publish_categorical_outputs(replacements, removals, temporary)

        if kind == "lc":
            self.land_use_layer = None
        elif kind == "geology":
            self.geology_path = None
        ready = getattr(self, "categorical_ready_outputs", {})
        ready[kind] = str(target)
        self.categorical_ready_outputs = ready
        self.load_layer(str(target), spec["layer_name"])
        self.mark_output_prepared(str(target), name=target.name)
        if kind == "soil":
            record = getattr(self.dialog, "record_standard_soil_output", None)
            if record is not None:
                definition_path = (
                    target.parent / definition if definition else None
                )
                record(
                    str(target),
                    str(definition_path) if definition_path else None,
                )
        self.log_message(f"Copied mHM-ready {spec['label']} data to {target}.")
        self._record_categorical_nml(kind, target, target.parent / definition if definition else None)
        return True

    def _record_categorical_nml(self, kind, output, definition=None):
        """Record a prepared categorical input for the namelist handoff."""
        from ...core.handlers.state.nml_settings import relative_workspace_path, update_section

        section = "land_cover" if kind == "lc" else kind
        value = {
            "mode": self._categorical_mode(kind),
            "output_path": relative_workspace_path(self.dialog.project_folder, output),
            "variable": {
                "lc": "land_cover",
                "soil": "soil_class",
                "geology": "geology_class",
            }[kind],
        }
        if definition is not None and Path(definition).is_file():
            value["classdefinition_path"] = relative_workspace_path(
                self.dialog.project_folder, definition
            )
        config = self._categorical_lookup(kind) or {}
        for key in ("lookup_table", "mapping_field", "class_field"):
            if config.get(key):
                value[key] = config[key]
        update_section(self.dialog.project_folder, section, value)

    def _categorical_intermediate_paths(self, kind):
        base = (
            Path(project_geometry_folder(self.dialog.project_folder))
            / _SPECS[kind]["geometry"]
        )
        return [
            base,
            base.with_name(f"{base.stem}_crop{base.suffix}"),
            base.with_name(f"{base.stem}_masked{base.suffix}"),
        ]

    def _ready_output_paths(self, kind):
        spec = _SPECS[kind]
        folder = Path(morph_folder(self.dialog.project_folder))
        paths = [
            folder / f"{spec['ready_name']}{suffix}"
            for suffix in (".asc", ".nc", ".tif", ".tiff")
        ]
        paths.extend(
            folder / f"{spec['ready_name']}{suffix}"
            for suffix in (".prj", ".PRJ")
        )
        return paths

    @staticmethod
    def _publish_categorical_outputs(replacements, removals, temporary):
        destinations = [Path(destination) for _, destination in replacements]
        tracked = []
        for path in destinations + [Path(path) for path in removals]:
            if path not in tracked:
                tracked.append(path)

        backups = []
        published = []
        try:
            for index, path in enumerate(tracked):
                if path.is_file():
                    backup = temporary / f".backup_{index}_{path.name}"
                    os.replace(path, backup)
                    backups.append((backup, path))
            for source, destination in replacements:
                os.replace(source, destination)
                published.append(Path(destination))
        except Exception:
            for path in published:
                if path.is_file():
                    path.unlink()
            for backup, path in reversed(backups):
                if backup.is_file():
                    os.replace(backup, path)
            raise

    def _categorical_error(self, spec, error):
        self.log_message(f"ERROR preparing {spec['label']} data: {error}")
        QMessageBox.critical(
            self.dialog,
            f"{spec['label']} Processing Error",
            f"Could not prepare the {spec['label']} data:\n{error}",
        )
        return False

    def _categorical_mode(self, kind):
        method = getattr(self.dialog, "categorical_input_mode", None)
        if method is not None:
            value = (
                str(method(kind) or "")
                .strip()
                .lower()
                .replace("-", "_")
                .replace(" ", "_")
            )
            if "lookup_table" in value:
                return "lookup_table"
            if value == "mhm_ready":
                return value
            return value
        return ""

    def _categorical_lookup(self, kind):
        method = getattr(self.dialog, "categorical_lookup_config", None)
        return method(kind) if method is not None else None

    def _input_crs(self, layer):
        value = self._layer_crs_text(layer)
        return value or self._dialog_crs_text()

    def _dem_crs(self):
        return self._layer_crs_text(
            QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        ) or self._dialog_crs_text()

    def _dialog_crs_text(self):
        method = getattr(self.dialog, "get_crs", None)
        return self._crs_text(method() if method is not None else None)

    @classmethod
    def _layer_crs_text(cls, layer):
        return cls._crs_text(layer.crs() if layer is not None else None)

    @staticmethod
    def _crs_text(crs):
        if crs is None or not crs.isValid():
            return None
        return crs.authid() or crs.toWkt()


__all__ = ["CategoricalProcessingMixin"]
