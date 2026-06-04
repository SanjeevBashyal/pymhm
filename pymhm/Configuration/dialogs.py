# -*- coding: utf-8 -*-
"""Qt dialogs for schema-driven namelist editing."""

from qgis.PyQt import QtCore, QtWidgets

from .namelist import canonical_name, parse_scalar, split_csv


class NamelistEditorDialog(QtWidgets.QDialog):
    """Multi-page editor for one namelist family."""

    def __init__(self, parent, title, pages, parameter_mode=False):
        super(NamelistEditorDialog, self).__init__(parent)
        self.pages = pages
        self.parameter_mode = parameter_mode
        self.field_widgets = {}
        self.parameter_widgets = {}

        self.setWindowTitle(title)
        self.resize(780, 620)
        self.build_ui()
        self.update_navigation()

    def build_ui(self):
        layout = QtWidgets.QVBoxLayout(self)

        page_bar = QtWidgets.QHBoxLayout()
        page_bar.addWidget(QtWidgets.QLabel("Page:"))
        self.page_combo = QtWidgets.QComboBox(self)
        for page in self.pages:
            self.page_combo.addItem(page["title"], page["block"])
        self.page_combo.currentIndexChanged.connect(self.goto_page)
        page_bar.addWidget(self.page_combo, 1)
        layout.addLayout(page_bar)

        self.stack = QtWidgets.QStackedWidget(self)
        for page in self.pages:
            if self.parameter_mode:
                self.stack.addWidget(self.build_parameter_page(page))
            else:
                self.stack.addWidget(self.build_general_page(page))
        layout.addWidget(self.stack, 1)

        button_bar = QtWidgets.QHBoxLayout()
        self.restore_button = QtWidgets.QPushButton("Restore Defaults", self)
        self.restore_button.clicked.connect(self.restore_current_page_defaults)
        self.back_button = QtWidgets.QPushButton("Back", self)
        self.back_button.clicked.connect(self.previous_page)
        self.next_button = QtWidgets.QPushButton("Next", self)
        self.next_button.clicked.connect(self.next_page)
        self.cancel_button = QtWidgets.QPushButton("Cancel", self)
        self.cancel_button.clicked.connect(self.reject)
        self.save_button = QtWidgets.QPushButton("OK", self)
        self.save_button.clicked.connect(self.accept)
        self.save_button.setDefault(True)

        button_bar.addWidget(self.restore_button)
        button_bar.addStretch(1)
        button_bar.addWidget(self.back_button)
        button_bar.addWidget(self.next_button)
        button_bar.addWidget(self.cancel_button)
        button_bar.addWidget(self.save_button)
        layout.addLayout(button_bar)

    def build_general_page(self, page):
        scroll = QtWidgets.QScrollArea(self)
        scroll.setWidgetResizable(True)
        content = QtWidgets.QWidget(scroll)
        form = QtWidgets.QFormLayout(content)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        form.setLabelAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)

        schema = page["schema"]
        defaults = page["defaults"]
        for name, prop_schema in schema.get("properties", {}).items():
            label = QtWidgets.QLabel(self.field_label(name, prop_schema), content)
            label.setToolTip(prop_schema.get("description", ""))
            widget = self.create_general_widget(prop_schema, defaults.get(name))
            widget.setToolTip(prop_schema.get("description", ""))
            self.field_widgets[(page["block"], name)] = {
                "widget": widget,
                "schema": prop_schema,
                "default": defaults.get(name),
            }
            form.addRow(label, widget)

        scroll.setWidget(content)
        return scroll

    def build_parameter_page(self, page):
        table = QtWidgets.QTableWidget(self)
        table.setColumnCount(6)
        table.setHorizontalHeaderLabels([
            "Parameter",
            "Lower bound",
            "Upper bound",
            "Value",
            "Flag",
            "Scaling",
        ])
        table.verticalHeader().setVisible(False)
        table.setAlternatingRowColors(True)
        table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)

        schema = page["schema"]
        defaults = page["defaults"]
        properties = list(schema.get("properties", {}).items())
        table.setRowCount(len(properties))
        for row, (name, prop_schema) in enumerate(properties):
            default = self.parameter_default(defaults.get(name), prop_schema)
            lower, upper, value, flag, scaling = default

            title = self.field_label(name, prop_schema)
            item = QtWidgets.QTableWidgetItem(title)
            item.setToolTip(prop_schema.get("description", ""))
            item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEditable)
            table.setItem(row, 0, item)

            lower_item = QtWidgets.QTableWidgetItem(str(lower))
            lower_item.setFlags(lower_item.flags() & ~QtCore.Qt.ItemIsEditable)
            table.setItem(row, 1, lower_item)

            upper_item = QtWidgets.QTableWidgetItem(str(upper))
            upper_item.setFlags(upper_item.flags() & ~QtCore.Qt.ItemIsEditable)
            table.setItem(row, 2, upper_item)

            value_edit = QtWidgets.QLineEdit(str(value), table)
            flag_combo = self.binary_combo(flag, table)
            scaling_combo = self.binary_combo(scaling, table)
            table.setCellWidget(row, 3, value_edit)
            table.setCellWidget(row, 4, flag_combo)
            table.setCellWidget(row, 5, scaling_combo)
            self.parameter_widgets[(page["block"], name)] = {
                "lower": lower,
                "upper": upper,
                "value": value_edit,
                "flag": flag_combo,
                "scaling": scaling_combo,
                "default": default,
            }

        table.resizeColumnsToContents()
        table.horizontalHeader().setStretchLastSection(True)
        return table

    def field_label(self, name, prop_schema):
        title = prop_schema.get("title") or name
        return f"{title} ({name})"

    def create_general_widget(self, prop_schema, value):
        value = self.general_default(value, prop_schema)
        if "enum" in prop_schema:
            combo = QtWidgets.QComboBox(self)
            for option in prop_schema.get("enum", []):
                combo.addItem(str(option), option)
            self.set_combo_value(combo, value)
            return combo

        prop_type = prop_schema.get("type")
        if prop_type == "boolean":
            checkbox = QtWidgets.QCheckBox(self)
            checkbox.setChecked(bool(value))
            return checkbox
        if prop_type == "integer":
            spin = QtWidgets.QSpinBox(self)
            spin.setRange(
                int(prop_schema.get("minimum", -2147483648)),
                int(prop_schema.get("maximum", 2147483647)))
            spin.setValue(int(value or 0))
            return spin
        if prop_type == "number":
            spin = QtWidgets.QDoubleSpinBox(self)
            spin.setDecimals(6)
            spin.setRange(
                float(prop_schema.get("minimum", -1.0e12)),
                float(prop_schema.get("maximum", 1.0e12)))
            spin.setValue(float(value or 0.0))
            return spin

        line_edit = QtWidgets.QLineEdit(self)
        line_edit.setText(self.value_to_text(value))
        return line_edit

    def binary_combo(self, value, parent):
        combo = QtWidgets.QComboBox(parent)
        combo.addItem("0", 0)
        combo.addItem("1", 1)
        self.set_combo_value(combo, int(float(value or 0)))
        return combo

    def set_combo_value(self, combo, value):
        for index in range(combo.count()):
            if str(combo.itemData(index)) == str(value):
                combo.setCurrentIndex(index)
                return
        if combo.count():
            combo.setCurrentIndex(0)

    def goto_page(self, index):
        self.stack.setCurrentIndex(index)
        self.update_navigation()

    def previous_page(self):
        index = max(0, self.stack.currentIndex() - 1)
        self.page_combo.setCurrentIndex(index)

    def next_page(self):
        index = min(self.stack.count() - 1, self.stack.currentIndex() + 1)
        self.page_combo.setCurrentIndex(index)

    def update_navigation(self):
        index = self.stack.currentIndex()
        self.back_button.setEnabled(index > 0)
        self.next_button.setEnabled(index < self.stack.count() - 1)

    def restore_current_page_defaults(self):
        page = self.pages[self.stack.currentIndex()]
        for name, prop_schema in page["schema"].get("properties", {}).items():
            if self.parameter_mode:
                self.restore_parameter_default(page["block"], name)
            else:
                self.restore_general_default(page["block"], name, prop_schema)

    def restore_general_default(self, block, name, prop_schema):
        info = self.field_widgets[(block, name)]
        widget = info["widget"]
        value = self.general_default(info["default"], prop_schema)
        if isinstance(widget, QtWidgets.QCheckBox):
            widget.setChecked(bool(value))
        elif isinstance(widget, QtWidgets.QSpinBox):
            widget.setValue(int(value or 0))
        elif isinstance(widget, QtWidgets.QDoubleSpinBox):
            widget.setValue(float(value or 0.0))
        elif isinstance(widget, QtWidgets.QComboBox):
            self.set_combo_value(widget, value)
        elif isinstance(widget, QtWidgets.QLineEdit):
            widget.setText(self.value_to_text(value))

    def restore_parameter_default(self, block, name):
        info = self.parameter_widgets[(block, name)]
        lower, upper, value, flag, scaling = info["default"]
        info["value"].setText(str(value))
        self.set_combo_value(info["flag"], int(float(flag or 0)))
        self.set_combo_value(info["scaling"], int(float(scaling or 0)))
        info["lower"] = lower
        info["upper"] = upper

    def collect_values(self):
        if self.parameter_mode:
            return self.collect_parameter_values()
        return self.collect_general_values()

    def collect_general_values(self):
        values = {}
        for (block, name), info in self.field_widgets.items():
            block_key = canonical_name(block)
            values.setdefault(block_key, {})
            values[block_key][canonical_name(name)] = self.read_general_widget(
                info["widget"], info["schema"])
        return values

    def collect_parameter_values(self):
        values = {}
        for (block, name), info in self.parameter_widgets.items():
            block_key = canonical_name(block)
            values.setdefault(block_key, {})
            values[block_key][canonical_name(name)] = [
                info["lower"],
                info["upper"],
                parse_scalar(info["value"].text()),
                int(info["flag"].currentData()),
                int(info["scaling"].currentData()),
            ]
        return values

    def read_general_widget(self, widget, prop_schema):
        if isinstance(widget, QtWidgets.QCheckBox):
            return widget.isChecked()
        if isinstance(widget, QtWidgets.QSpinBox):
            return widget.value()
        if isinstance(widget, QtWidgets.QDoubleSpinBox):
            return widget.value()
        if isinstance(widget, QtWidgets.QComboBox):
            return widget.currentData()
        if prop_schema.get("type") == "array":
            return self.parse_array_text(widget.text(), prop_schema)
        return self.parse_text_value(widget.text(), prop_schema)

    def parse_array_text(self, text, prop_schema):
        item_schema = prop_schema.get("items", {})
        if not text.strip():
            return ""
        return [
            self.parse_text_value(part.strip(), item_schema)
            for part in split_csv(text)
        ]

    def parse_text_value(self, text, prop_schema):
        prop_type = prop_schema.get("type")
        if prop_type == "boolean":
            return str(text).strip().lower() in (".true.", "true", "1", "yes")
        if prop_type == "integer":
            return int(float(str(text).strip() or 0))
        if prop_type == "number":
            return float(str(text).strip() or 0.0)
        return str(text)

    def general_default(self, value, prop_schema):
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
        if isinstance(value, (list, tuple)) and len(value) >= 5:
            default = list(value[:5])
        else:
            examples = prop_schema.get("examples") or []
            example = examples[0] if examples else [0.0, 0.0, 0.0, 0, 1]
            default = list(example[:5])
        while len(default) < 5:
            default.append(0 if len(default) != 4 else 1)
        return default

    def value_to_text(self, value):
        if isinstance(value, (list, tuple)):
            return ", ".join(str(item) for item in value)
        if value is True:
            return ".true."
        if value is False:
            return ".false."
        return "" if value is None else str(value)
