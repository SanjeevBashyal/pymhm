# -*- coding: utf-8 -*-
"""Array editor widgets for schema-driven configuration dialogs."""

from __future__ import annotations

from qgis.PyQt import QtCore, QtWidgets


MONTH_LABELS = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December",
]
NIGHT_HOUR_LABELS = [
    "18:00", "19:00", "20:00", "21:00", "22:00", "23:00",
    "00:00", "01:00", "02:00", "03:00", "04:00", "05:00",
]


def shape_parts(prop_schema):
    """Return x-fortran-shape as a list."""
    shape = prop_schema.get("x-fortran-shape", [])
    if isinstance(shape, (list, tuple)):
        return list(shape)
    if shape:
        return [shape]
    return []


def is_max_domains(value):
    """Return True when a shape dimension means max_domains."""
    return str(value).lower() == "max_domains"


def is_boolean_item(item_schema):
    """Return True when array items are booleans."""
    return item_schema.get("type") == "boolean"


def is_string_item(item_schema):
    """Return True when array items are strings."""
    return item_schema.get("type") == "string"


def parse_scalar_text(text, item_schema):
    """Parse editor text according to an item schema."""
    prop_type = item_schema.get("type")
    text = str(text).strip()
    if prop_type == "boolean":
        return text.lower() in (".true.", "true", "1", "yes", "y")
    if prop_type == "integer":
        return int(float(text or 0))
    if prop_type == "number":
        return float(text or 0.0)
    return text


def value_text(value):
    """Format a value for an editor field."""
    if value is True:
        return ".true."
    if value is False:
        return ".false."
    return "" if value is None else str(value)


def scalar_editor(value, item_schema, parent):
    """Create an editor for one scalar array item."""
    if is_boolean_item(item_schema):
        checkbox = QtWidgets.QCheckBox(parent)
        checkbox.setChecked(bool(value))
        return checkbox
    return QtWidgets.QLineEdit(value_text(value), parent)


def scalar_editor_value(editor, item_schema):
    """Read one scalar array item editor."""
    if isinstance(editor, QtWidgets.QCheckBox):
        return editor.isChecked()
    return parse_scalar_text(editor.text(), item_schema)


def set_scalar_editor_value(editor, value):
    """Set one scalar array item editor."""
    if isinstance(editor, QtWidgets.QCheckBox):
        editor.setChecked(bool(value))
    else:
        editor.setText(value_text(value))


def derived_default(schema, value=None):
    """Return one simple derived value with intrinsic scalar components."""
    source = value if isinstance(value, dict) else {}
    result = {}
    for name, component in schema.get("properties", {}).items():
        if component.get("type") in ("object", "array"):
            raise ValueError("derived components must be intrinsic scalars")
        source_name = next(
            (candidate for candidate in source
             if str(candidate).lower() == name.lower()),
            None,
        )
        if source_name is not None:
            result[name] = source[source_name]
        elif "default" in component:
            result[name] = component["default"]
        elif component.get("type") == "boolean":
            result[name] = False
        elif component.get("type") in ("integer", "number"):
            result[name] = component.get("minimum", 0)
        else:
            result[name] = ""
    return result


class DerivedValueWidget(QtWidgets.QWidget):
    """Inline editor for one simple Fortran derived value."""

    def __init__(self, schema, value=None, parent=None):
        super(DerivedValueWidget, self).__init__(parent)
        self.schema = schema
        self.editors = {}
        layout = QtWidgets.QFormLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        values = derived_default(schema, value)
        for name, component in schema.get("properties", {}).items():
            editor = scalar_editor(values[name], component, self)
            self.editors[name] = editor
            layout.addRow(component.get("title") or name, editor)

    def value(self):
        return {
            name: scalar_editor_value(editor, self.schema["properties"][name])
            for name, editor in self.editors.items()
        }

    def value_map(self, base_key):
        return {base_key: self.value()}

    def reset_value(self, value):
        values = derived_default(self.schema, value)
        for name, editor in self.editors.items():
            set_scalar_editor_value(editor, values[name])


class DerivedArrayWidget(QtWidgets.QTableWidget):
    """Table editor for a one-dimensional array of derived values."""

    def __init__(self, schema, value, domain_infos, parent=None):
        self.schema = schema
        self.item_schema = schema.get("items", {})
        self.component_names = list(self.item_schema.get("properties", {}))
        self.domain_infos = domain_infos or [
            {"station_id": "1", "label": "1", "index": 1}]
        count = len(self.domain_infos) if (
            shape_parts(schema) and is_max_domains(shape_parts(schema)[0])
        ) else max(1, len(value) if isinstance(value, list) else 1)
        super(DerivedArrayWidget, self).__init__(
            count, len(self.component_names) + 1, parent)
        self.setHorizontalHeaderLabels(["Domain", *self.component_names])
        self.editors = []
        self.reset_value(value)
        self.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.verticalHeader().setVisible(False)
        self.setMinimumHeight(min(260, 72 + 32 * count))

    def reset_value(self, value):
        source = value if isinstance(value, list) else []
        self.editors = []
        for row in range(self.rowCount()):
            label = self.domain_infos[row]["label"] if (
                row < len(self.domain_infos)) else str(row + 1)
            self.setItem(row, 0, QtWidgets.QTableWidgetItem(str(label)))
            values = derived_default(
                self.item_schema,
                source[row] if row < len(source) else None,
            )
            row_editors = {}
            for column, name in enumerate(self.component_names, start=1):
                component = self.item_schema["properties"][name]
                editor = scalar_editor(values[name], component, self)
                self.setCellWidget(row, column, editor)
                row_editors[name] = editor
            self.editors.append(row_editors)

    def value(self):
        return [
            {
                name: scalar_editor_value(
                    editor, self.item_schema["properties"][name])
                for name, editor in row.items()
            }
            for row in self.editors
        ]

    def value_map(self, base_key):
        return {base_key: self.value()}


def normalize_list_value(value, count, default):
    """Normalize one list to a requested length."""
    if isinstance(value, tuple):
        value = list(value)
    if isinstance(value, dict) and "__indexed__" in value:
        indexed = value.get("__indexed__", {})
        value = indexed.get("1") or indexed.get(1) or default
    if not isinstance(value, list):
        value = [] if value in (None, "") else [value]
    if value and isinstance(value[0], list):
        value = value[0]
    values = list(value)
    fill = values[-1] if values else default
    while len(values) < count:
        values.append(fill)
    return values[:count]


class InlineArrayWidget(QtWidgets.QWidget):
    """Inline editor for one-value arrays that still serialize as arrays."""

    def __init__(self, prop_schema, value, parent=None):
        super(InlineArrayWidget, self).__init__(parent)
        self.prop_schema = prop_schema
        self.item_schema = prop_schema.get("items", {})
        self._value = normalize_list_value(value, 1, self.item_default())

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.editor = scalar_editor(self._value[0], self.item_schema, self)
        layout.addWidget(self.editor, 1)

    def item_default(self):
        """Return the scalar item default."""
        if "default" in self.item_schema:
            return self.item_schema["default"]
        if "default" in self.prop_schema:
            return self.prop_schema["default"]
        examples = self.prop_schema.get("examples") or []
        if examples:
            example = examples[0]
            if isinstance(example, list) and example:
                return example[0]
            return example
        if is_boolean_item(self.item_schema):
            return False
        if self.item_schema.get("type") == "integer":
            return 0
        if self.item_schema.get("type") == "number":
            return 0.0
        return ""

    def value(self):
        """Return this one-value array."""
        return [scalar_editor_value(self.editor, self.item_schema)]

    def value_map(self, base_key):
        """Return namelist values keyed for rendering."""
        return {base_key: self.value()}

    def reset_value(self, value):
        """Reset the inline editor to a new default."""
        self._value = normalize_list_value(value, 1, self.item_default())
        set_scalar_editor_value(self.editor, self._value[0])


class ArrayValueWidget(QtWidgets.QWidget):
    """Compact form widget that opens a dedicated array editor dialog."""

    def __init__(self, field_name, prop_schema, value, domain_infos, parent=None):
        super(ArrayValueWidget, self).__init__(parent)
        self.field_name = field_name
        self.prop_schema = prop_schema
        self.domain_infos = domain_infos or [
            {"station_id": "1", "label": "1", "index": 1}]
        self.item_schema = prop_schema.get("items", {})
        self.shape = shape_parts(prop_schema)
        self._multi_domain = (
            len(self.shape) >= 2 and is_max_domains(self.shape[-1]))
        self._domain_1d = (
            len(self.shape) == 1 and is_max_domains(self.shape[0]))
        self._value = self.normalized_value(value)

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.summary_label = QtWidgets.QLabel(self)
        self.edit_button = QtWidgets.QPushButton("Edit", self)
        self.edit_button.clicked.connect(self.open_editor)
        layout.addWidget(self.edit_button, 0)
        layout.addWidget(self.summary_label, 1)
        self.update_summary()

    def normalized_value(self, value):
        """Normalize incoming defaults to the editor value shape."""
        if self._multi_domain:
            return self.normalize_multi_domain(value)
        if self._domain_1d:
            return self.normalize_list(
                value, len(self.domain_infos), self.domain_default())
        labels = self.first_dimension_labels(value)
        return self.normalize_list(value, len(labels), self.item_default())

    def normalize_multi_domain(self, value):
        """Normalize a multidimensional array into one list per domain."""
        domain_count = len(self.domain_infos)
        row_count = len(self.first_dimension_labels(value))
        default_row = self.item_default()

        if isinstance(value, dict):
            indexed = value.get("__indexed__", {})
            rows = []
            for index in range(1, domain_count + 1):
                rows.append(self.normalize_list(
                    indexed.get(str(index)) or indexed.get(index),
                    row_count,
                    default_row,
                ))
            return rows

        if isinstance(value, list) and value and isinstance(value[0], list):
            rows = []
            for index in range(domain_count):
                source = value[index] if index < len(value) else value[0]
                rows.append(self.normalize_list(source, row_count, default_row))
            return rows

        base = self.normalize_list(value, row_count, default_row)
        return [list(base) for _ in range(domain_count)]

    def normalize_list(self, value, count, default):
        """Normalize one list to a requested length."""
        return normalize_list_value(value, count, default)

    def domain_default(self):
        """Return default for domain-indexed arrays."""
        return self.item_default()

    def item_default(self):
        """Return the scalar item default."""
        if "default" in self.item_schema:
            return self.item_schema["default"]
        if "default" in self.prop_schema:
            return self.prop_schema["default"]
        examples = self.prop_schema.get("examples") or []
        if examples:
            example = examples[0]
            if isinstance(example, list) and example:
                return example[0]
            return example
        prop_type = self.item_schema.get("type")
        if prop_type == "boolean":
            return False
        if prop_type == "integer":
            return 0
        if prop_type == "number":
            return 0.0
        return ""

    def first_dimension_labels(self, value=None):
        """Return labels for the first array dimension."""
        first_dim = self.shape[0] if self.shape else None
        if self.field_name.startswith("frac_night_"):
            return NIGHT_HOUR_LABELS
        if first_dim == 12:
            return MONTH_LABELS
        if first_dim == 24:
            return [f"{hour:02d}:00" for hour in range(24)]
        if str(first_dim).lower() == "max_layers":
            source = value[0] if (
                isinstance(value, list) and value and isinstance(value[0], list)
            ) else value
            count = len(source) if isinstance(source, list) else 3
            return [f"Layer {index}" for index in range(1, max(1, count) + 1)]
        if isinstance(first_dim, int):
            return [str(index) for index in range(1, first_dim + 1)]
        if self._domain_1d:
            return [item["label"] for item in self.domain_infos]
        if isinstance(value, list) and value:
            return [str(index) for index in range(1, len(value) + 1)]
        return ["1"]

    def open_editor(self):
        """Open the appropriate array editor."""
        if self._multi_domain:
            dialog = MultiDomainArrayDialog(
                self.field_name,
                self.first_dimension_labels(self._value),
                self.domain_infos,
                self._value,
                self.item_schema,
                self,
            )
        else:
            labels = (
                [item["label"] for item in self.domain_infos]
                if self._domain_1d
                else self.first_dimension_labels(self._value)
            )
            dialog = SingleArrayDialog(
                self.field_name,
                labels,
                self._value,
                self.item_schema,
                self,
            )

        if dialog.exec_() == dialog.Accepted:
            self._value = dialog.values()
            self.update_summary()

    def update_summary(self):
        """Refresh the compact value summary."""
        if self._multi_domain:
            self.summary_label.setText(
                f"{len(self._value)} domain page(s), "
                f"{len(self._value[0]) if self._value else 0} value(s) each")
        elif self._domain_1d and "resolution" in self.field_name.lower():
            self.summary_label.setText(
                ", ".join(value_text(value) for value in self._value))
        else:
            self.summary_label.setText(f"{len(self._value)} value(s)")

    def value_map(self, base_key):
        """Return namelist values keyed for rendering."""
        if self._multi_domain:
            return {
                f"{base_key}__domain{index + 1}": values
                for index, values in enumerate(self._value)
            }
        return {base_key: self._value}

    def value(self):
        """Return raw editor value for state persistence."""
        return self._value


class SingleArrayDialog(QtWidgets.QDialog):
    """Single-page array editor."""

    def __init__(self, title, labels, values, item_schema, parent=None):
        super(SingleArrayDialog, self).__init__(parent)
        self.item_schema = item_schema
        self.editors = []
        self.setWindowTitle(f"Edit {title}")
        self.resize(420, 520)

        layout = QtWidgets.QVBoxLayout(self)
        scroll = QtWidgets.QScrollArea(self)
        scroll.setWidgetResizable(True)
        content = QtWidgets.QWidget(scroll)
        grid = QtWidgets.QGridLayout(content)
        for row, label in enumerate(labels):
            grid.addWidget(QtWidgets.QLabel(str(label), content), row, 0)
            editor = scalar_editor(values[row], item_schema, content)
            self.editors.append(editor)
            grid.addWidget(editor, row, 1)
        scroll.setWidget(content)
        layout.addWidget(scroll)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel,
            self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def values(self):
        """Return edited values."""
        return [
            scalar_editor_value(editor, self.item_schema)
            for editor in self.editors
        ]


class MultiDomainArrayDialog(QtWidgets.QDialog):
    """Multi-page array editor with one page per domain/STATION_ID."""

    def __init__(self, title, row_labels, domain_infos, values, item_schema,
                 parent=None):
        super(MultiDomainArrayDialog, self).__init__(parent)
        self.row_labels = row_labels
        self.domain_infos = domain_infos
        self.item_schema = item_schema
        self.page_editors = []
        self.setWindowTitle(f"Edit {title}")
        self.resize(520, 560)

        layout = QtWidgets.QVBoxLayout(self)
        selector = QtWidgets.QHBoxLayout()
        selector.addWidget(QtWidgets.QLabel("STATION_ID:", self))
        self.page_combo = QtWidgets.QComboBox(self)
        for domain in domain_infos:
            self.page_combo.addItem(domain["label"])
        self.page_combo.currentIndexChanged.connect(self.goto_page)
        selector.addWidget(self.page_combo, 1)
        layout.addLayout(selector)

        self.stack = QtWidgets.QStackedWidget(self)
        for domain_index, domain in enumerate(domain_infos):
            page = self.build_domain_page(
                domain, values[domain_index], item_schema)
            self.stack.addWidget(page)
        layout.addWidget(self.stack, 1)

        controls = QtWidgets.QHBoxLayout()
        copy_button = QtWidgets.QPushButton("Copy to All", self)
        copy_button.clicked.connect(self.copy_current_to_all)
        self.back_button = QtWidgets.QPushButton("Back", self)
        self.back_button.clicked.connect(self.previous_page)
        self.next_button = QtWidgets.QPushButton("Next", self)
        self.next_button.clicked.connect(self.next_page)
        controls.addWidget(copy_button)
        controls.addStretch(1)
        controls.addWidget(self.back_button)
        controls.addWidget(self.next_button)
        layout.addLayout(controls)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel,
            self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.update_navigation()

    def build_domain_page(self, domain, values, item_schema):
        """Build one domain page."""
        scroll = QtWidgets.QScrollArea(self)
        scroll.setWidgetResizable(True)
        content = QtWidgets.QWidget(scroll)
        grid = QtWidgets.QGridLayout(content)
        editors = []
        heading = QtWidgets.QLabel(
            f"STATION_ID {domain['label']}", content)
        heading.setAlignment(QtCore.Qt.AlignCenter)
        grid.addWidget(heading, 0, 0, 1, 2)
        for offset, label in enumerate(self.row_labels, start=1):
            grid.addWidget(QtWidgets.QLabel(str(label), content), offset, 0)
            editor = scalar_editor(
                values[offset - 1], item_schema, content)
            editors.append(editor)
            grid.addWidget(editor, offset, 1)
        self.page_editors.append(editors)
        scroll.setWidget(content)
        return scroll

    def goto_page(self, index):
        self.stack.setCurrentIndex(index)
        self.update_navigation()

    def previous_page(self):
        self.page_combo.setCurrentIndex(max(0, self.stack.currentIndex() - 1))

    def next_page(self):
        self.page_combo.setCurrentIndex(
            min(self.stack.count() - 1, self.stack.currentIndex() + 1))

    def update_navigation(self):
        index = self.stack.currentIndex()
        self.back_button.setEnabled(index > 0)
        self.next_button.setEnabled(index < self.stack.count() - 1)

    def copy_current_to_all(self):
        """Copy values from the active domain page to every domain."""
        current = [
            scalar_editor_value(editor, self.item_schema)
            for editor in self.page_editors[self.stack.currentIndex()]
        ]
        for editors in self.page_editors:
            for editor, value in zip(editors, current):
                set_scalar_editor_value(editor, value)

    def values(self):
        """Return edited values per domain."""
        return [
            [
                scalar_editor_value(editor, self.item_schema)
                for editor in editors
            ]
            for editors in self.page_editors
        ]
