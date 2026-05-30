# -*- coding: utf-8 -*-
"""Vector dataset cleanup and filtered polygon-copy helpers."""
from ..common import (
    os,
    QgsVectorLayer,
    QgsVectorFileWriter,
    NULL,
)


class VectorIOMixin:
    """Vector dataset cleanup and filtered polygon-copy helpers."""

    def _remove_vector_dataset(self, path):
        """Remove a shapefile and its sidecar files if they exist."""
        base, _ = os.path.splitext(path)
        for extension in (".shp", ".shx", ".dbf", ".prj", ".cpg", ".qpj", ".fix"):
            sidecar = base + extension
            if os.path.exists(sidecar):
                try:
                    os.remove(sidecar)
                except Exception as e:
                    self.log_message(f"WARNING: Could not remove {sidecar}: {e}")

    def _copy_nonzero_polygons(self, input_path, output_path, field_name="DN"):
        """Copy polygonized basin features with non-zero IDs to a clean shapefile."""
        source_layer = QgsVectorLayer(input_path, "Watershed_Raw", "ogr")
        if not source_layer.isValid():
            self.log_message(f"ERROR: Could not load polygonized watershed: {input_path}")
            return False

        self._remove_vector_dataset(output_path)
        writer = QgsVectorFileWriter(
            output_path,
            "UTF-8",
            source_layer.fields(),
            source_layer.wkbType(),
            source_layer.crs(),
            "ESRI Shapefile"
        )
        if writer.hasError() != QgsVectorFileWriter.NoError:
            self.log_message(f"ERROR creating watershed vector: {writer.errorMessage()}")
            return False

        copied = 0
        for feature in source_layer.getFeatures():
            value = feature.attribute(field_name)
            try:
                keep_feature = value not in (None, NULL) and float(value) != 0.0
            except (TypeError, ValueError):
                keep_feature = False

            if keep_feature:
                writer.addFeature(feature)
                copied += 1

        del writer
        if copied == 0:
            self.log_message("WARNING: Polygonized watershed had no non-zero basin polygons.")
            return False

        if os.path.exists(output_path):
            self.mark_output_prepared(
                output_path,
                name=os.path.basename(output_path),
                loaded=False
            )
            return True
        return False
