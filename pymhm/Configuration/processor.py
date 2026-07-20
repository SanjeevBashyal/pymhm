# -*- coding: utf-8 -*-
"""Processor for schema-driven mHM configuration namelists."""
from __future__ import annotations

import copy
import os
from typing import Any

from nml_tools import json_to_namelist
from qgis.PyQt.QtWidgets import QFileDialog, QMessageBox

from ..grid_resolution import ceil_cellsize
from ..simulation_processor import SimulationProcessor
from ..project_layout import ensure_project_structure
from .constants import (
    PARAMETER_ALWAYS_BLOCKS,
    PROCESS_PARAMETER_MAP,
    STATUS_LABELS,
    STATUS_TITLES,
)
from .dialogs import NamelistEditorDialog
from .domain_info import (
    domain_count,
    domain_infos,
    geology_class_count,
    geology_class_rows,
    geology_parameter_values,
)
from .namelist import (
    canonical_name,
    template_block_order,
    template_values,
)
from .nmljson import NmlJsonStore, atomic_write_text
from .path_defaults import configuration_path_defaults
from .paths import namelist_output_dir, output_path, template_path, version_key
from .schema_loader import config_schemas, output_schemas, parameter_schema_lookup
from .profile_values import namelist_values
from .state import ConfigurationStateStore


ConfigPage = dict[str, Any]
ConfigValues = dict[str, dict[str, Any]]


class ConfigurationProcessor(SimulationProcessor):
    """Handle mHM configuration dialogs, namelist writing, and status labels."""

    def __init__(self, dialog: Any) -> None:
        super(ConfigurationProcessor, self).__init__(dialog)
        self.state_store = ConfigurationStateStore(dialog)
        self.nml_store = NmlJsonStore(dialog)
        self.configuration_state = self.state_store.empty()
        self.nml_document = self.nml_store.empty()
        self._last_mhm_values: ConfigValues | None = None

    def selected_version(self) -> str:
        """Return the selected mHM version, defaulting to schema-compatible v6."""
        combo = getattr(self.dialog, "comboBox_mHMversion", None)
        if combo is None:
            return "6.0"
        return combo.currentText().strip() or "6.0"

    def load_project_state(self) -> None:
        """Load saved Configuration-tab state for the selected project."""
        self.configuration_state = self.state_store.load()
        self.nml_document = self.nml_store.load(self.selected_version())
        version = (
            self.configuration_state.get("mhm_version")
            or self.nml_document.get("mhm_version")
        )
        if version and hasattr(self.dialog, "comboBox_mHMversion"):
            combo = self.dialog.comboBox_mHMversion
            index = combo.findText(version)
            if index >= 0:
                combo.setCurrentIndex(index)

        legacy = self.configuration_state.pop(
            "_legacy_namelist_settings", None)
        if isinstance(legacy, dict):
            self.migrate_legacy_settings(legacy)

        if hasattr(self.dialog, "lineEdit_loadConfiguration"):
            self.dialog.lineEdit_loadConfiguration.setText(
                self.nml_store.path() or "")

        self._last_mhm_values = self.kind_state("mhm") or None
        self.refresh_status_indicators()
        path = self.nml_store.path()
        if path and os.path.exists(path):
            self.log_message(f"Namelist configuration loaded: {path}")

    def browse_configuration_file(self) -> bool:
        """Import a namelist-profile JSON file selected by the user."""
        if not self.ensure_project_folder():
            return False

        file_path, _ = QFileDialog.getOpenFileName(
            self.dialog,
            "Select Configuration Settings File",
            "",
            "JSON (*.json);;All Files (*)",
        )
        if not file_path:
            return False
        try:
            import json
            with open(file_path, "r", encoding="utf-8") as profile_file:
                document = json.load(profile_file)
            if not isinstance(document, dict):
                raise ValueError("Namelist configuration is not a JSON object.")
            if not isinstance(document.get("file_profiles"), dict):
                raise ValueError("Namelist configuration has no file_profiles object.")
            self.nml_document = document
            self.nml_store.save(self.nml_document)
            if hasattr(self.dialog, "lineEdit_loadConfiguration"):
                self.dialog.lineEdit_loadConfiguration.setText(
                    self.nml_store.path() or file_path)
            version = document.get("mhm_version")
            if version and hasattr(self.dialog, "comboBox_mHMversion"):
                index = self.dialog.comboBox_mHMversion.findText(version)
                if index >= 0:
                    self.dialog.comboBox_mHMversion.setCurrentIndex(index)
            self.save_state()
            self._last_mhm_values = self.kind_state("mhm") or None
            self.refresh_status_indicators()
            self.log_message(f"Namelist configuration imported: {file_path}")
            return True
        except Exception as exc:
            self.report_error("Load Configuration Settings", exc)
            return False

    def save_state(self) -> None:
        """Persist current configuration state."""
        self.configuration_state["mhm_version"] = self.selected_version()
        self.state_store.save(self.configuration_state)
        if hasattr(self.dialog, "lineEdit_loadConfiguration"):
            self.dialog.lineEdit_loadConfiguration.setText(
                self.nml_store.path() or "")

    def handle_version_changed(self) -> None:
        """Persist the selected version and refresh project paths/status."""
        if self.dialog.project_folder:
            ensure_project_structure(
                self.dialog.project_folder,
                self.selected_version())
            self.save_state()
        self.refresh_status_indicators()

    def kind_state(self, kind: str) -> dict[str, Any]:
        """Return saved values for a namelist kind."""
        profile = self.nml_store.profile(self.nml_document, kind)
        values = profile.get("editor_values")
        return values if isinstance(values, dict) else {}

    def migrate_legacy_settings(self, settings: dict[str, Any]) -> None:
        """Move legacy state-file namelist settings into nmljson profiles."""
        mhm_values = settings.get("mhm") or {}
        pages_by_kind = {
            "mhm": self.build_config_pages(),
            "parameters": self.build_parameter_pages(mhm_values),
            "outputs": self.build_output_pages(),
        }
        for kind, pages in pages_by_kind.items():
            editor_values = settings.get(kind)
            profile = self.nml_store.profile(self.nml_document, kind)
            if not isinstance(editor_values, dict) or profile.get("values"):
                continue
            converted = namelist_values(
                self.selected_version(),
                kind,
                editor_values,
                pages,
                self.dialog,
            )
            self.nml_store.set_profile(
                self.nml_document,
                kind,
                converted,
                self.selected_version(),
                editor_values=editor_values,
                dimensions=self.profile_dimensions(kind),
            )
        self.nml_store.save(self.nml_document)
        self.state_store.save(self.configuration_state)

    def configure_mhm(self) -> bool:
        """Open the mHM configuration dialog and save mhm.nml on OK."""
        return self.open_editor(
            kind="mhm",
            title="mHM Configuration",
            pages=self.build_config_pages(),
            parameter_mode=False,
        )

    def configure_parameters(self) -> bool:
        """Open the parameter configuration dialog and save parameter namelist."""
        config_values = self.current_mhm_values()
        return self.open_editor(
            kind="parameters",
            title="mHM Parameters",
            pages=self.build_parameter_pages(config_values),
            parameter_mode=True,
        )

    def configure_outputs(self) -> bool:
        """Open the output configuration dialog and save mhm_outputs.nml."""
        return self.open_editor(
            kind="outputs",
            title="mHM Outputs",
            pages=self.build_output_pages(),
            parameter_mode=False,
        )

    def open_editor(
            self,
            kind: str,
            title: str,
            pages: list[ConfigPage],
            parameter_mode: bool) -> bool:
        """Open a namelist editor and write its output when accepted."""
        if not self.ensure_project_folder():
            return False
        if not pages:
            QMessageBox.warning(
                self.dialog,
                "Configuration",
                f"No schema pages are available for {title}.")
            return False

        try:
            editor = NamelistEditorDialog(
                self.dialog,
                title,
                pages,
                parameter_mode=parameter_mode,
                domain_infos=domain_infos(self.dialog),
            )
            if editor.exec_() != editor.Accepted:
                return False
            values = editor.collect_values()
            self.write_kind(kind, values, pages=pages)
            if kind == "mhm":
                self._last_mhm_values = values
            return True
        except Exception as exc:
            self.report_error(title, exc)
            return False

    def create_nml_files(self) -> bool:
        """Create all schema-driven namelist files using current/default values."""
        if not self.check_prerequisites():
            return False
        try:
            self.log_message("\n--- Creating schema-driven namelist files ---")
            config_pages = self.build_config_pages()
            config_values = self.kind_state("mhm") or self.default_values_from_pages(config_pages)
            config_values = self.apply_current_resolution_values(config_values)
            self.write_kind("mhm", config_values, pages=config_pages)
            self._last_mhm_values = config_values

            parameter_pages = self.build_parameter_pages(config_values)
            parameter_values = self.default_values_from_pages(
                parameter_pages,
                parameter_mode=True,
            )
            self.write_kind(
                "parameters",
                parameter_values,
                pages=parameter_pages,
            )

            output_pages = self.build_output_pages()
            output_values = (
                self.kind_state("outputs")
                or self.default_values_from_pages(output_pages)
            )
            self.write_kind("outputs", output_values, pages=output_pages)

            self.log_message("Schema-driven namelist files created successfully.")
            return True
        except Exception as exc:
            self.report_error("Create Namelists", exc)
            return False

    def write_kind(
            self,
            kind: str,
            values: dict[str, Any],
            pages: list[ConfigPage] | None = None) -> str:
        """Save one file profile and render it through nml-tools."""
        project_folder = self.dialog.project_folder
        ensure_project_structure(project_folder, self.selected_version())
        os.makedirs(namelist_output_dir(project_folder), exist_ok=True)
        destination = output_path(project_folder, kind, self.selected_version())
        if pages is None:
            pages = self.pages_for_kind(kind)
        source_values = (
            self.apply_current_resolution_values(values)
            if kind == "mhm" else values
        )
        if kind == "parameters" and version_key(self.selected_version()) == "v5.13":
            pages = self.build_parameter_pages({}, include_all=True)
            source_values = self.merge_generated_defaults(
                self.default_values_from_pages(pages, parameter_mode=True),
                values,
            )
        render_values = namelist_values(
            self.selected_version(),
            kind,
            source_values,
            pages,
            self.dialog,
        )
        self.nml_store.set_profile(
            self.nml_document,
            kind,
            render_values,
            self.selected_version(),
            editor_values=values,
            dimensions=self.profile_dimensions(kind),
        )
        profile = self.nml_store.profile(self.nml_document, kind)
        atomic_write_text(destination, json_to_namelist(profile))
        self.nml_store.save(self.nml_document)
        self.save_state()
        self.log_message(f"Written: {destination}")
        self.refresh_status_indicators()
        return destination

    def pages_for_kind(
            self,
            kind: str) -> list[ConfigPage]:
        """Return schema pages when a caller did not already build them."""
        if kind == "mhm":
            return self.build_config_pages()
        if kind == "parameters":
            return self.build_parameter_pages(self.current_mhm_values())
        if kind == "outputs":
            return self.build_output_pages()
        raise ValueError(f"Unknown namelist profile kind: {kind}")

    def profile_dimensions(self, kind: str) -> dict[str, int]:
        """Return dimensions needed to interpret the saved profile arrays."""
        if kind == "parameters":
            return {"max_geo_units": geology_class_count(self.dialog)}
        return {"max_domains": domain_count(self.dialog)}

    def build_config_pages(self) -> list[ConfigPage]:
        """Build page data for the mHM configuration dialog."""
        return self.pages_from_schemas(config_schemas(self.selected_version()), "mhm")

    def build_output_pages(self) -> list[ConfigPage]:
        """Build page data for the output configuration dialog."""
        return self.pages_from_schemas(output_schemas(self.selected_version()), "outputs")

    def build_parameter_pages(
            self,
            config_values: dict[str, Any],
            include_all: bool = False) -> list[ConfigPage]:
        """Build parameter pages required by selected mHM process cases."""
        schema_lookup = parameter_schema_lookup(self.selected_version())
        ordered_blocks = self.ordered_template_blocks("parameters")
        schema_blocks = [
            "petm2" if block == "pet0" else block
            for block in ordered_blocks
        ]
        selected_blocks = (
            schema_blocks
            if include_all
            else self.parameter_blocks_for_config(config_values)
        )
        selected_keys = {canonical_name(block) for block in selected_blocks}

        ordered_selected = [
            block for block in schema_blocks
            if canonical_name(block) in selected_keys
        ]
        for block in selected_blocks:
            key = canonical_name(block)
            if key not in [canonical_name(name) for name in ordered_selected]:
                ordered_selected.append(block)

        schemas = []
        for block in ordered_selected:
            schema = schema_lookup.get(canonical_name(block))
            if schema:
                schemas.append(schema)
        return self.pages_from_schemas(schemas, "parameters")

    def parameter_blocks_for_config(self, config_values: dict[str, Any]) -> list[str]:
        """Return parameter namelist blocks required by process selections."""
        process_values = self.process_values_from_config(config_values)
        selected_blocks = []
        for process_name, case_map in PROCESS_PARAMETER_MAP.items():
            case_value = process_values.get(canonical_name(process_name), 0)
            try:
                case_value = int(case_value)
            except (TypeError, ValueError):
                case_value = 0
            block = case_map.get(case_value)
            if block:
                selected_blocks.append(block)
        selected_blocks.extend(PARAMETER_ALWAYS_BLOCKS)
        return selected_blocks

    def process_values_from_config(self, config_values: dict[str, Any]) -> dict[str, Any]:
        """Extract config_processes values from current mHM config values."""
        processes = config_values.get(canonical_name("config_processes"), {})
        if processes:
            return processes
        defaults = self.default_values_from_pages(self.build_config_pages())
        return defaults.get(canonical_name("config_processes"), {})

    def current_mhm_values(self) -> dict[str, Any]:
        """Return current mHM config values from memory, file, or defaults."""
        if self._last_mhm_values:
            return self.apply_current_resolution_values(self._last_mhm_values)
        project_folder = self.dialog.project_folder
        if project_folder:
            values = self.kind_state("mhm")
            if values:
                return self.apply_current_resolution_values(values)
            path = output_path(project_folder, "mhm", self.selected_version())
            if os.path.exists(path):
                return self.apply_current_resolution_values(template_values(path))
        return self.apply_current_resolution_values(
            self.default_values_from_pages(self.build_config_pages()))

    def pages_from_schemas(
            self,
            schemas: list[dict[str, Any]],
            kind: str) -> list[ConfigPage]:
        """Create dialog page dictionaries from schemas and template defaults."""
        defaults = template_values(template_path(self.selected_version(), kind))
        saved = self.kind_state(kind)
        if kind == "mhm":
            generated = configuration_path_defaults(domain_count(self.dialog))
            generated = self.merge_generated_defaults(
                generated,
                self.current_resolution_defaults(),
            )
        else:
            generated = {}
        dynamic = self.current_resolution_defaults() if kind == "mhm" else {}
        pages = []
        for schema in schemas:
            block = schema.get("x-fortran-namelist")
            if not block:
                continue
            block_key = canonical_name(block)
            block_defaults = defaults.get(canonical_name(block), {})
            block_generated = generated.get(block_key, {})
            block_saved = saved.get(block_key, {})
            block_dynamic = dynamic.get(block_key, {})
            geology_rows = []
            if kind == "parameters" and block_key == "geoparameter":
                geology_rows = geology_class_rows(self.dialog)
                block_generated = self.generated_geoparameter_defaults(
                    block_defaults,
                    schema.get("properties", {}).get("GeoParam", {}),
                )
            page_defaults = {}
            for name in schema.get("properties", {}):
                key = canonical_name(name)
                if key in block_dynamic:
                    page_defaults[name] = block_dynamic[key]
                else:
                    page_defaults[name] = self.default_for_property(
                        block_defaults,
                        block_generated,
                        block_saved,
                        name,
                        schema["properties"][name],
                    )
            pages.append({
                "block": block,
                "title": schema.get("title") or block,
                "schema": schema,
                "defaults": page_defaults,
            })
            if block_key == "geoparameter":
                if not geology_rows:
                    geology_rows = geology_class_rows(self.dialog)
                pages[-1]["geo_classes"] = geology_rows
                pages[-1]["geo_class_count"] = (
                    len(geology_rows) or geology_class_count(self.dialog))
        return pages

    def merge_generated_defaults(
            self,
            base: ConfigValues,
            overlay: ConfigValues) -> ConfigValues:
        """Merge generated defaults without mutating the source dictionaries."""
        merged = {
            canonical_name(block): dict(values or {})
            for block, values in (base or {}).items()
        }
        for block, values in (overlay or {}).items():
            block_key = canonical_name(block)
            merged.setdefault(block_key, {})
            merged[block_key].update({
                canonical_name(name): value
                for name, value in (values or {}).items()
            })
        return merged

    def current_resolution_defaults(self) -> ConfigValues:
        """Return current L1/L11 resolutions as config defaults for all domains."""
        unit = self.current_grid_unit()
        domain_total = max(1, int(domain_count(self.dialog) or 1))
        version = version_key(self.selected_version())
        defaults: ConfigValues = {}

        l1_resolution = self.current_resolution_value("current_l1_resolution", unit)
        if l1_resolution is not None:
            name = "resolution_Hydrology" if version == "v5.13" else "resolution"
            defaults.setdefault(canonical_name("config_mhm"), {})[
                canonical_name(name)
            ] = [l1_resolution for _ in range(domain_total)]

        l11_resolution = self.current_resolution_value("current_l11_resolution", unit)
        if l11_resolution is not None:
            name = "resolution_Routing" if version == "v5.13" else "resolution"
            defaults.setdefault(canonical_name("config_mrm"), {})[
                canonical_name(name)
            ] = [l11_resolution for _ in range(domain_total)]

        return defaults

    def apply_current_resolution_values(self, values: ConfigValues) -> ConfigValues:
        """Overlay current L1/L11 resolution values onto mHM config values."""
        return self.merge_generated_defaults(values, self.current_resolution_defaults())

    def current_resolution_value(self, method_name: str, unit: str | None) -> float | None:
        """Read and round a current dialog resolution for configuration use."""
        if not hasattr(self.dialog, method_name):
            return None
        try:
            value = getattr(self.dialog, method_name)()
        except Exception:
            return None
        try:
            value = float(value)
        except (TypeError, ValueError):
            return None
        if value <= 0:
            return None
        return ceil_cellsize(value, unit)

    def current_grid_unit(self) -> str:
        """Return the active grid unit from the dialog, if available."""
        if hasattr(self.dialog, "current_grid_unit"):
            try:
                return self.dialog.current_grid_unit() or ""
            except Exception:
                return ""
        return ""

    def default_for_property(
            self,
            template_defaults: dict[str, Any],
            generated_defaults: dict[str, Any],
            saved_defaults: dict[str, Any],
            name: str,
            prop_schema: dict[str, Any]) -> Any:
        """Return defaults for a schema property from saved/generated/template values."""
        base_key = canonical_name(name)
        for source in (saved_defaults, generated_defaults, template_defaults):
            if base_key in source:
                return source[base_key]

        if canonical_name(name) == "geoparam":
            rows = self.indexed_defaults(
                template_defaults, base_key, "class")
            rows.update(self.indexed_defaults(
                generated_defaults, base_key, "class"))
            rows.update(self.indexed_defaults(
                saved_defaults, base_key, "class"))
            return rows

        indexed = self.indexed_defaults(template_defaults, base_key, "domain")
        indexed.update(self.indexed_defaults(
            generated_defaults, base_key, "domain"))
        indexed.update(self.indexed_defaults(
            saved_defaults, base_key, "domain"))
        if indexed:
            return {"__indexed__": indexed}
        return None

    def indexed_defaults(
            self,
            values: dict[str, Any],
            base_key: str,
            suffix: str) -> dict[str, Any]:
        """Return indexed defaults for base_key as {index: value}."""
        prefix = f"{base_key}__{suffix}"
        found = {}
        for key, value in values.items():
            if not key.startswith(prefix):
                continue
            index = key[len(prefix):]
            if index.isdigit():
                found[index] = value
        return found

    def generated_geoparameter_defaults(
            self,
            template_defaults: dict[str, Any],
            prop_schema: dict[str, Any]) -> dict[str, Any]:
        """Return GeoParam defaults generated from geology lookup metadata."""
        parameter_values = geology_parameter_values(self.dialog)
        if not parameter_values:
            return {}

        generated = {}
        template_rows = self.indexed_defaults(
            template_defaults,
            canonical_name("GeoParam"),
            "class",
        )
        for class_index, parameter_value in parameter_values.items():
            default = self.parameter_default(
                template_rows.get(str(class_index)),
                prop_schema,
            )
            default[2] = parameter_value
            generated[f"geoparam__class{class_index}"] = default
        return generated

    def default_values_from_pages(
            self,
            pages: list[ConfigPage],
            parameter_mode: bool = False) -> dict[str, Any]:
        """Collect default values from page data without opening an editor."""
        values = {}
        for page in pages:
            block_key = canonical_name(page["block"])
            values.setdefault(block_key, {})
            if parameter_mode and block_key == "geoparameter":
                prop_schema = page["schema"].get("properties", {}).get(
                    "GeoParam", {})
                default_rows = page["defaults"].get("GeoParam") or {}
                for class_index in self.geoparam_class_indices(page):
                    values[block_key][f"geoparam__class{class_index}"] = (
                        self.parameter_default(
                            default_rows.get(str(class_index)), prop_schema))
                continue

            for name, prop_schema in page["schema"].get("properties", {}).items():
                default = page["defaults"].get(name)
                if parameter_mode:
                    values[block_key][canonical_name(name)] = (
                        self.parameter_default(default, prop_schema))
                else:
                    values[block_key].update(
                        self.default_entries_for_property(
                            canonical_name(name), default, prop_schema))
        return values

    def default_entries_for_property(
            self,
            base_key: str,
            default: Any,
            prop_schema: dict[str, Any]) -> dict[str, Any]:
        """Return namelist-ready entries for a non-parameter property."""
        if prop_schema.get("type") != "array":
            return {base_key: self.general_default(default, prop_schema)}

        if isinstance(default, dict) and "__indexed__" in default:
            return {
                f"{base_key}__domain{index}": value
                for index, value in default["__indexed__"].items()
            }
        return {base_key: self.general_default(default, prop_schema)}

    def geoparam_class_indices(self, page: ConfigPage) -> list[int]:
        """Return GeoParam row indices from geology metadata or class count."""
        indices = []
        for row in page.get("geo_classes", []) or []:
            try:
                indices.append(int(row.get("geo_param")))
            except (AttributeError, TypeError, ValueError):
                continue
        if indices:
            return indices

        class_count = int(page.get("geo_class_count") or 16)
        return list(range(1, class_count + 1))

    def general_default(self, value: Any, prop_schema: dict[str, Any]) -> Any:
        """Return a non-parameter default from template, schema, or type."""
        if value is not None:
            return value
        prop_type = prop_schema.get("type")
        if prop_type == "array":
            examples = prop_schema.get("examples") or []
            if examples:
                return copy.deepcopy(
                    examples[0] if isinstance(examples[0], list) else examples)
            item = prop_schema.get("items", {})
            item_value = self.general_default(None, item)
            shape = prop_schema.get("x-fortran-shape")
            dimensions = shape if isinstance(shape, list) else [shape]

            def filled(index: int) -> Any:
                if index >= len(dimensions):
                    return copy.deepcopy(item_value)
                dimension = dimensions[index]
                count = (
                    domain_count(self.dialog)
                    if str(dimension).lower() == "max_domains"
                    else dimension if isinstance(dimension, int) else 1
                )
                return [filled(index + 1) for _ in range(max(1, count))]

            return filled(0)
        if prop_type == "object":
            return {
                name: self.general_default(None, component)
                for name, component in prop_schema.get("properties", {}).items()
            }
        if "default" in prop_schema:
            return prop_schema["default"]
        examples = prop_schema.get("examples") or []
        if examples:
            return examples[0]
        if prop_type == "boolean":
            return False
        if prop_type == "integer":
            return 0
        if prop_type == "number":
            return 0.0
        return ""

    def parameter_default(self, value: Any, prop_schema: dict[str, Any]) -> list[Any]:
        """Return [lower, upper, default, flag, scaling] for a parameter."""
        if isinstance(value, (list, tuple)) and len(value) >= 5:
            default = list(value[:5])
        else:
            examples = prop_schema.get("examples") or []
            example = examples[0] if examples else [0.0, 0.0, 0.0, 0, 1]
            default = list(example[:5])
        while len(default) < 5:
            default.append(0 if len(default) != 4 else 1)
        return default

    def ordered_template_blocks(self, kind: str) -> list[str]:
        """Return namelist block order for the selected version template."""
        return template_block_order(template_path(self.selected_version(), kind))

    def ensure_project_folder(self) -> bool:
        """Ensure a project folder is selected."""
        if self.dialog.project_folder:
            return True
        QMessageBox.warning(
            self.dialog,
            "Project Folder Required",
            "Select a project folder before preparing namelist files.")
        return False

    def refresh_status_indicators(self) -> None:
        """Update red/green status labels for saved namelist files."""
        project_folder = self.dialog.project_folder
        for kind, label_name in STATUS_LABELS.items():
            label = getattr(self.dialog, label_name, None)
            if label is None:
                continue
            path = (
                output_path(project_folder, kind, self.selected_version())
                if project_folder else ""
            )
            saved = bool(path and os.path.exists(path))
            self.set_status_label(label, kind, path, saved)

    def set_status_label(
            self,
            label: Any,
            kind: str,
            path: str,
            saved: bool) -> None:
        """Style one status label."""
        color = "#1f8f4d" if saved else "#b42318"
        label.setText("Saved" if saved else "Not saved")
        title = STATUS_TITLES.get(kind, kind)
        if saved:
            tooltip = f"{title} namelist saved:\n{path}"
        elif path:
            tooltip = f"{title} namelist has not been saved yet:\n{path}"
        else:
            tooltip = f"{title} namelist has not been saved yet."
        label.setToolTip(tooltip)
        label.setStyleSheet(
            "QLabel {"
            f"background-color: {color};"
            "color: white;"
            "border-radius: 8px;"
            "padding: 2px 8px;"
            "font-weight: 600;"
            "}"
        )

    def run_mhm(self) -> bool:
        """Run mHM in the selected project folder through the project terminal."""
        if not self.ensure_project_folder():
            return False

        self.log_message("\n--- Running mHM ---")
        project_folder = self.dialog.project_folder
        ensure_project_structure(project_folder, self.selected_version())
        self.log_message(f"Running mHM in project directory: {project_folder}")
        terminal = self.dialog.open_project_terminal()
        if terminal is None:
            return False
        self.log_message("Executing command in project terminal: mhm")
        return terminal.run_command("mhm", cwd=project_folder, show=True)

    def report_error(self, title: str, exc: Exception) -> None:
        """Log and show a configuration error."""
        self.log_message(f"ERROR: {title} failed. Details: {exc}")
        QMessageBox.critical(self.dialog, "Configuration Error", str(exc))
