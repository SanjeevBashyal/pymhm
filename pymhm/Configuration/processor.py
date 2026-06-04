# -*- coding: utf-8 -*-
"""Processor for schema-driven mHM configuration namelists."""

import os

from qgis.PyQt.QtWidgets import QMessageBox

from ..simulation_processor import SimulationProcessor
from .constants import (
    PARAMETER_ALWAYS_BLOCKS,
    PROCESS_PARAMETER_MAP,
    STATUS_LABELS,
    STATUS_TITLES,
)
from .dialogs import NamelistEditorDialog
from .namelist import (
    canonical_name,
    template_block_order,
    template_values,
    write_rendered_namelist,
)
from .paths import model_inputs_dir, output_path, template_path
from .schema_loader import config_schemas, output_schemas, parameter_schema_lookup


class ConfigurationProcessor(SimulationProcessor):
    """Handle mHM configuration dialogs, namelist writing, and status labels."""

    def __init__(self, dialog):
        super(ConfigurationProcessor, self).__init__(dialog)
        self._last_mhm_values = None

    def selected_version(self):
        """Return the selected mHM version, defaulting to schema-compatible v6."""
        combo = getattr(self.dialog, "comboBox_mHMversion", None)
        if combo is None:
            return "6.0"
        return combo.currentText().strip() or "6.0"

    def configure_mhm(self):
        """Open the mHM configuration dialog and save mhm.nml on OK."""
        return self.open_editor(
            kind="mhm",
            title="mHM Configuration",
            pages=self.build_config_pages(),
            parameter_mode=False,
        )

    def configure_parameters(self):
        """Open the parameter configuration dialog and save mhm_parameters.nml."""
        config_values = self.current_mhm_values()
        return self.open_editor(
            kind="parameters",
            title="mHM Parameters",
            pages=self.build_parameter_pages(config_values),
            parameter_mode=True,
        )

    def configure_outputs(self):
        """Open the output configuration dialog and save mhm_outputs.nml."""
        return self.open_editor(
            kind="outputs",
            title="mHM Outputs",
            pages=self.build_output_pages(),
            parameter_mode=False,
        )

    def open_editor(self, kind, title, pages, parameter_mode):
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
                self.dialog, title, pages, parameter_mode=parameter_mode)
            if editor.exec_() != editor.Accepted:
                return False
            values = editor.collect_values()
            include_blocks = [page["block"] for page in pages] if parameter_mode else None
            self.write_kind(kind, values, include_blocks=include_blocks)
            if kind == "mhm":
                self._last_mhm_values = values
            return True
        except Exception as exc:
            self.report_error(title, exc)
            return False

    def create_nml_files(self):
        """Create all schema-driven namelist files using current/default values."""
        if not self.check_prerequisites():
            return False
        try:
            self.log_message("\n--- Creating schema-driven namelist files ---")
            config_pages = self.build_config_pages()
            config_values = self.default_values_from_pages(config_pages)
            self.write_kind("mhm", config_values)
            self._last_mhm_values = config_values

            parameter_pages = self.build_parameter_pages(config_values)
            parameter_values = self.default_values_from_pages(
                parameter_pages, parameter_mode=True)
            self.write_kind(
                "parameters",
                parameter_values,
                include_blocks=[page["block"] for page in parameter_pages],
            )

            output_pages = self.build_output_pages()
            output_values = self.default_values_from_pages(output_pages)
            self.write_kind("outputs", output_values)

            self.log_message("Schema-driven namelist files created successfully.")
            return True
        except Exception as exc:
            self.report_error("Create Namelists", exc)
            return False

    def write_kind(self, kind, values, include_blocks=None):
        """Render and save one namelist kind."""
        project_folder = self.dialog.project_folder
        os.makedirs(model_inputs_dir(project_folder), exist_ok=True)
        template = template_path(self.selected_version(), kind)
        destination = output_path(project_folder, kind)
        write_rendered_namelist(
            template, destination, values, include_blocks=include_blocks)
        self.log_message(f"Written: {destination}")
        self.refresh_status_indicators()
        return destination

    def build_config_pages(self):
        """Build page data for the mHM configuration dialog."""
        return self.pages_from_schemas(config_schemas(), "mhm")

    def build_output_pages(self):
        """Build page data for the output configuration dialog."""
        return self.pages_from_schemas(output_schemas(), "outputs")

    def build_parameter_pages(self, config_values):
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

    def parameter_blocks_for_config(self, config_values):
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

    def process_values_from_config(self, config_values):
        """Extract config_processes values from current mHM config values."""
        processes = config_values.get(canonical_name("config_processes"), {})
        if processes:
            return processes
        defaults = self.default_values_from_pages(self.build_config_pages())
        return defaults.get(canonical_name("config_processes"), {})

    def current_mhm_values(self):
        """Return current mHM config values from memory, file, or defaults."""
        if self._last_mhm_values:
            return self._last_mhm_values
        project_folder = self.dialog.project_folder
        if project_folder:
            path = output_path(project_folder, "mhm")
            if os.path.exists(path):
                return template_values(path)
        return self.default_values_from_pages(self.build_config_pages())

    def pages_from_schemas(self, schemas, kind):
        """Create dialog page dictionaries from schemas and template defaults."""
        defaults = template_values(template_path(self.selected_version(), kind))
        pages = []
        for schema in schemas:
            block = schema.get("x-fortran-namelist")
            if not block:
                continue
            block_defaults = defaults.get(canonical_name(block), {})
            page_defaults = {}
            for name in schema.get("properties", {}):
                page_defaults[name] = block_defaults.get(canonical_name(name))
            pages.append({
                "block": block,
                "title": schema.get("title") or block,
                "schema": schema,
                "defaults": page_defaults,
            })
        return pages

    def default_values_from_pages(self, pages, parameter_mode=False):
        """Collect default values from page data without opening an editor."""
        values = {}
        for page in pages:
            block_key = canonical_name(page["block"])
            values.setdefault(block_key, {})
            for name, prop_schema in page["schema"].get("properties", {}).items():
                default = page["defaults"].get(name)
                if parameter_mode:
                    values[block_key][canonical_name(name)] = (
                        self.parameter_default(default, prop_schema))
                else:
                    values[block_key][canonical_name(name)] = (
                        self.general_default(default, prop_schema))
        return values

    def general_default(self, value, prop_schema):
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

    def parameter_default(self, value, prop_schema):
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

    def ordered_template_blocks(self, kind):
        """Return namelist block order for the selected version template."""
        return template_block_order(template_path(self.selected_version(), kind))

    def ensure_project_folder(self):
        """Ensure a project folder is selected."""
        if self.dialog.project_folder:
            return True
        QMessageBox.warning(
            self.dialog,
            "Project Folder Required",
            "Select a project folder before preparing namelist files.")
        return False

    def refresh_status_indicators(self):
        """Update red/green status labels for saved namelist files."""
        project_folder = self.dialog.project_folder
        for kind, label_name in STATUS_LABELS.items():
            label = getattr(self.dialog, label_name, None)
            if label is None:
                continue
            path = output_path(project_folder, kind) if project_folder else ""
            saved = bool(path and os.path.exists(path))
            self.set_status_label(label, kind, path, saved)

    def set_status_label(self, label, kind, path, saved):
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

    def report_error(self, title, exc):
        """Log and show a configuration error."""
        self.log_message(f"ERROR: {title} failed. Details: {exc}")
        QMessageBox.critical(self.dialog, "Configuration Error", str(exc))
