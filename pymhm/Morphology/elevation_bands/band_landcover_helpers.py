# -*- coding: utf-8 -*-
"""Land-cover helpers for elevation-band detail tables."""
from ..common import (
    os,
    project_geometry_folder,
    math,
    QMessageBox,
    processing,
)


class BandLandCoverHelperMixin:
    """Land-cover helpers for elevation-band detail tables."""

    def _ensure_land_use_raster(self):
        """Return a processed land-use raster, creating it if required."""
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        expected_path = os.path.join(geometry_folder, "3_land_use.tif")

        if self.land_use_layer and os.path.exists(self.land_use_layer):
            return self.land_use_layer
        if os.path.exists(expected_path):
            self.land_use_layer = expected_path
            return expected_path

        land_cover_layer = self.dialog.mMapLayerComboBox_land_cover.currentLayer()
        if not land_cover_layer:
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Please select a Land Cover raster/vector layer.")
            return None

        self.log_message("Processed land-use raster is missing. Running Land Use first...")
        self.without_layer_loading(self.process_land_use)

        if self.land_use_layer and os.path.exists(self.land_use_layer):
            return self.land_use_layer
        if os.path.exists(expected_path):
            self.land_use_layer = expected_path
            return expected_path

        self.log_message("ERROR: Land Use processing did not create 3_land_use.tif.")
        return None

    def _canonical_land_cover_class(self, value):
        """Canonical numeric class value for stable CSV headers/comparisons."""
        numeric_value = float(value)
        rounded_value = round(numeric_value)
        if math.isclose(numeric_value, rounded_value, abs_tol=1e-9):
            return int(rounded_value)
        return numeric_value

    def _land_cover_class_header(self, class_value, class_names=None):
        """Return a CSV-safe raw land-cover class area column name."""
        if isinstance(class_value, int):
            suffix = str(class_value)
        else:
            suffix = f"{class_value:g}".replace("-", "minus_").replace(".", "_")
        value_suffix = self._clean_output_name(suffix, "unknown")

        class_name = None
        if class_names:
            class_name = class_names.get(class_value)
            if class_name is None:
                try:
                    class_name = class_names.get(int(float(class_value)))
                except (TypeError, ValueError):
                    class_name = None

        if class_name:
            name_suffix = self._clean_output_name(class_name, "unknown")
            return f"land_cover_class_{value_suffix}_{name_suffix}_area"

        return f"land_cover_class_{value_suffix}_area"

    def _format_optional_float(self, value):
        """Format optional float values for CSV output."""
        if value is None:
            return ""
        return f"{float(value):.6f}"
