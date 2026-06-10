# -*- coding: utf-8 -*-
"""Processor for schema-driven mHM configuration namelists."""
from __future__ import annotations

import os
import subprocess
from typing import Any

from qgis.PyQt.QtWidgets import QFileDialog, QMessageBox

from ..simulation_processor import SimulationProcessor
from ..project_layout import ensure_project_structure
from .constants import (
    PARAMETER_ALWAYS_BLOCKS,
    PROCESS_PARAMETER_MAP,
    STATUS_LABELS,
    STATUS_TITLES,
)
from .dialogs import NamelistEditorDialog
from .domain_info import domain_count, domain_infos, geology_class_count
from .namelist import (
    canonical_name,
    template_block_order,
    template_values,
    write_rendered_namelist,
)
from .path_defaults import configuration_path_defaults
from .paths import namelist_output_dir, output_path, template_path
from .schema_loader import config_schemas, output_schemas, parameter_schema_lookup
from .state import ConfigurationStateStore


ConfigPage = dict[str, Any]
ConfigValues = dict[str, dict[str, Any]]


class ConfigurationProcessor(SimulationProcessor):
    """Handle mHM configuration dialogs, namelist writing, and status labels."""

    def __init__(self, dialog: Any) -> None:
        super(ConfigurationProcessor, self).__init__(dialog)
        self.state_store = ConfigurationStateStore(dialog)
        self.configuration_state = self.state_store.empty()
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
        if hasattr(self.dialog, "lineEdit_loadConfiguration"):
            path = self.state_store.path()
            self.dialog.lineEdit_loadConfiguration.setText(path or "")

        version = self.configuration_state.get("mhm_version")
        if version and hasattr(self.dialog, "comboBox_mHMversion"):
            combo = self.dialog.comboBox_mHMversion
            index = combo.findText(version)
            if index >= 0:
                combo.setCurrentIndex(index)

        self._last_mhm_values = self.kind_state("mhm") or None
        self.refresh_status_indicators()
        path = self.state_store.path()
        if path and os.path.exists(path):
            self.log_message(f"Configuration state loaded: {path}")

    def browse_configuration_file(self) -> bool:
        """Load a configuration settings JSON file selected by the user."""
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
            with open(file_path, "r", encoding="utf-8") as state_file:
                state = json.load(state_file)
            if not isinstance(state, dict):
                raise ValueError("Configuration settings file is not a JSON object.")
            state.setdefault("settings", {})
            self.configuration_state = state
            self.save_state()
            if hasattr(self.dialog, "lineEdit_loadConfiguration"):
                self.dialog.lineEdit_loadConfiguration.setText(
                    self.state_store.path() or file_path)
            version = state.get("mhm_version")
            if version and hasattr(self.dialog, "comboBox_mHMversion"):
                index = self.dialog.comboBox_mHMversion.findText(version)
                if index >= 0:
                    self.dialog.comboBox_mHMversion.setCurrentIndex(index)
            self._last_mhm_values = self.kind_state("mhm") or None
            self.refresh_status_indicators()
            self.log_message(f"Configuration settings imported: {file_path}")
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
                self.state_store.path() or "")

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
        return (
            self.configuration_state
            .setdefault("settings", {})
            .setdefault(kind, {})
        )

    def set_kind_state(self, kind: str, values: dict[str, Any]) -> None:
        """Save values for a namelist kind into the state object."""
        self.configuration_state.setdefault("settings", {})[kind] = values
        self.save_state()

    def configure_mhm(self) -> bool:
        """Open the mHM configuration dialog and save mhm.nml on OK."""
        return self.open_editor(
            kind="mhm",
            title="mHM Configuration",
            pages=self.build_config_pages(),
            parameter_mode=False,
        )

    def configure_parameters(self) -> bool:
        """Open the parameter configuration dialog and save mhm_parameters.nml."""
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
            include_blocks = [page["block"] for page in pages] if parameter_mode else None
            self.set_kind_state(kind, values)
            self.write_kind(kind, values, include_blocks=include_blocks)
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
            self.set_kind_state("mhm", config_values)
            self.write_kind("mhm", config_values)
            self._last_mhm_values = config_values

            parameter_pages = self.build_parameter_pages(config_values)
            parameter_values = (
                self.kind_state("parameters")
                or self.default_values_from_pages(
                    parameter_pages, parameter_mode=True)
            )
            self.set_kind_state("parameters", parameter_values)
            self.write_kind(
                "parameters",
                parameter_values,
                include_blocks=[page["block"] for page in parameter_pages],
            )

            output_pages = self.build_output_pages()
            output_values = self.kind_state("outputs") or self.default_values_from_pages(output_pages)
            self.set_kind_state("outputs", output_values)
            self.write_kind("outputs", output_values)

            self.log_message("Schema-driven namelist files created successfully.")
            return True
        except Exception as exc:
            self.report_error("Create Namelists", exc)
            return False

    def write_kind(
            self,
            kind: str,
            values: dict[str, Any],
            include_blocks: list[str] | None = None) -> str:
        """Render and save one namelist kind."""
        project_folder = self.dialog.project_folder
        ensure_project_structure(project_folder, self.selected_version())
        os.makedirs(namelist_output_dir(project_folder), exist_ok=True)
        template = template_path(self.selected_version(), kind)
        destination = output_path(project_folder, kind)
        write_rendered_namelist(
            template, destination, values, include_blocks=include_blocks)
        self.log_message(f"Written: {destination}")
        self.refresh_status_indicators()
        return destination

    def build_config_pages(self) -> list[ConfigPage]:
        """Build page data for the mHM configuration dialog."""
        return self.pages_from_schemas(config_schemas(), "mhm")

    def build_output_pages(self) -> list[ConfigPage]:
        """Build page data for the output configuration dialog."""
        return self.pages_from_schemas(output_schemas(), "outputs")

    def build_parameter_pages(self, config_values: dict[str, Any]) -> list[ConfigPage]:
        """Build parameter pages required by selected mHM process cases."""
        schema_lookup = parameter_schema_lookup()
        selected_blocks = self.parameter_blocks_for_config(config_values)
        ordered_blocks = self.ordered_template_blocks("parameters")
        selected_keys = {canonical_name(block) for block in selected_blocks}

        ordered_selected = [
            block for block in ordered_blocks
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
            return self._last_mhm_values
        project_folder = self.dialog.project_folder
        if project_folder:
            values = self.kind_state("mhm")
            if values:
                return values
            path = output_path(project_folder, "mhm")
            if os.path.exists(path):
                return template_values(path)
        return self.default_values_from_pages(self.build_config_pages())

    def pages_from_schemas(
            self,
            schemas: list[dict[str, Any]],
            kind: str) -> list[ConfigPage]:
        """Create dialog page dictionaries from schemas and template defaults."""
        defaults = template_values(template_path(self.selected_version(), kind))
        saved = self.kind_state(kind)
        if kind == "mhm":
            generated = configuration_path_defaults(domain_count(self.dialog))
        else:
            generated = {}
        pages = []
        for schema in schemas:
            block = schema.get("x-fortran-namelist")
            if not block:
                continue
            block_defaults = defaults.get(canonical_name(block), {})
            block_generated = generated.get(canonical_name(block), {})
            block_saved = saved.get(canonical_name(block), {})
            page_defaults = {}
            for name in schema.get("properties", {}):
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
            if canonical_name(block) == "geoparameter":
                pages[-1]["geo_class_count"] = geology_class_count(self.dialog)
        return pages

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
                class_count = int(page.get("geo_class_count") or 16)
                for class_index in range(1, class_count + 1):
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

    def general_default(self, value: Any, prop_schema: dict[str, Any]) -> Any:
        """Return a non-parameter default from template, schema, or type."""
        if value is not None:
            return value
        if "default" in prop_schema:
            return prop_schema["default"]
        item_default = prop_schema.get("items", {}).get("default")
        if item_default is not None:
            return item_default
        examples = prop_schema.get("examples") or []
        if examples:
            return examples[0]
        prop_type = prop_schema.get("type")
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
            path = output_path(project_folder, kind) if project_folder else ""
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
        """Create root-level namelists and run mHM against the project folder."""
        if not self.check_prerequisites():
            return False

        self.log_message("\n--- Running mHM ---")
        self.log_message("Creating NML files...")
        if not self.create_nml_files():
            self.log_message("ERROR: Failed to create NML files. Cannot run mHM.")
            return False

        project_folder = self.dialog.project_folder
        ensure_project_structure(project_folder, self.selected_version())
        self.log_message(f"Running mHM in project directory: {project_folder}")
        command = ["mHM", "."]

        try:
            return self.run_mhm_command(command, project_folder)
        except FileNotFoundError:
            self.log_message("mHM not found in PATH. Trying conda environments...")

        for env_name in ("mhm", "mHM", "base", "conda"):
            command = ["conda", "run", "-n", env_name, "mHM", "."]
            try:
                if self.run_mhm_command(command, project_folder):
                    return True
            except FileNotFoundError:
                self.log_message("conda not found in PATH.")
                break

        QMessageBox.critical(
            self.dialog,
            "mHM Error",
            "Could not run mHM. Ensure mHM is available in PATH or a conda environment.")
        return False

    def run_mhm_command(self, command: list[str], cwd: str) -> bool:
        """Run one mHM command and stream its output to the plugin console."""
        self.log_message(f"Executing command: {' '.join(command)}")
        process = subprocess.Popen(
            command,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )

        self.log_message("mHM execution started...")
        for line in process.stdout:
            line = line.strip()
            if line:
                self.log_message(line)
        process.wait()

        if process.returncode == 0:
            self.log_message("mHM execution completed successfully.")
            QMessageBox.information(
                self.dialog,
                "mHM",
                "mHM execution completed successfully.")
            return True

        self.log_message(
            f"mHM command failed with return code {process.returncode}.")
        return False

    def report_error(self, title: str, exc: Exception) -> None:
        """Log and show a configuration error."""
        self.log_message(f"ERROR: {title} failed. Details: {exc}")
        QMessageBox.critical(self.dialog, "Configuration Error", str(exc))
