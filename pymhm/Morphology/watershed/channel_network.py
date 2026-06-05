# -*- coding: utf-8 -*-
"""Channel network extraction from flow accumulation."""
from ..common import (
    os,
    project_geometry_folder,
    QgsRasterLayer,
    QgsFeature,
    QgsGeometry,
    QgsVectorFileWriter,
    QgsField,
    QgsFields,
    QgsWkbTypes,
    QgsPointXY,
    QVariant,
    processing,
)


class ChannelNetworkMixin:
    """Channel network extraction from flow accumulation."""

    def process_channel_network(self):
        """Create a stream network vector from pyflwdir flow accumulation."""
        if not self._ensure_filled_dem():
            return

        # Keep the UI workflow consistent with the notebook: filled DEM,
        # flow accumulation, then stream extraction.
        if not self._ensure_flow_accumulation():
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        self.channel_network_vector_path = os.path.join(
            geometry_folder, "2_channel_network.shp")

        # Check if Channel Network already exists, otherwise process it
        if self.channel_network_vector_path and os.path.exists(self.channel_network_vector_path):
            self.log_message(
                "Channel Network already exists. Loading existing file...")
            self.load_layer(self.channel_network_vector_path,
                            "2_Channel_Network", is_raster=False)
        else:
            self.log_message("Processing Channel Network with pyflwdir...")

            context = self._build_flwdir_from_filled_dem()
            if not context:
                self.log_message("Channel Network processing failed.")
                self.channel_network_vector_path = None
                return

            deps = context["deps"]
            np = deps["np"]
            flwdir = context["flwdir"]
            reference = context["reference"]
            invalid_mask = context["invalid_mask"]

            flow_accumulation = flwdir.upstream_area(unit="cell")
            flow_accumulation[invalid_mask] = 0

            cell_area_m2 = self._reference_cell_area_m2(reference, deps)
            channel_area_threshold_m2 = 10000000.0
            threshold_cells = max(
                1, int(round(channel_area_threshold_m2 / cell_area_m2)))

            valid_accumulation = flow_accumulation[flow_accumulation > 0]
            if valid_accumulation.size == 0:
                self.log_message("ERROR: No valid flow accumulation cells found.")
                self.channel_network_vector_path = None
                return

            max_accumulation = int(np.nanmax(valid_accumulation))
            if threshold_cells > max_accumulation:
                threshold_cells = max(1, int(np.nanpercentile(valid_accumulation, 95)))
                self.log_message(
                    "Channel threshold exceeded the basin accumulation. "
                    f"Using 95th percentile threshold: {threshold_cells} cells.")

            stream_mask = flow_accumulation >= threshold_cells
            stream_order = flwdir.stream_order("strahler", mask=stream_mask)
            features = flwdir.streams(
                mask=stream_mask,
                strord=stream_order,
                uparea=flow_accumulation
            )

            if not features:
                self.log_message("ERROR: pyflwdir did not produce any stream features.")
                self.channel_network_vector_path = None
                return

            fields = QgsFields()
            fields.append(QgsField("Order", QVariant.Int))
            fields.append(QgsField("UpArea", QVariant.Double))
            fields.append(QgsField("idx", QVariant.Int))
            fields.append(QgsField("idx_ds", QVariant.Int))
            fields.append(QgsField("pit", QVariant.Int))

            filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
            output_crs = filled_dem_layer.crs() if filled_dem_layer.isValid() else self.dialog.get_crs()

            self._remove_vector_dataset(self.channel_network_vector_path)
            writer = QgsVectorFileWriter(
                self.channel_network_vector_path,
                "UTF-8",
                fields,
                QgsWkbTypes.LineString,
                output_crs,
                "ESRI Shapefile"
            )
            if writer.hasError() != QgsVectorFileWriter.NoError:
                self.log_message(f"ERROR creating channel network: {writer.errorMessage()}")
                self.channel_network_vector_path = None
                return

            written = 0
            for pyflwdir_feature in features:
                coordinates = pyflwdir_feature.get("geometry", {}).get("coordinates", [])
                if len(coordinates) < 2:
                    continue

                properties = pyflwdir_feature.get("properties", {})
                qgs_feature = QgsFeature(fields)
                qgs_feature.setGeometry(QgsGeometry.fromPolylineXY([
                    QgsPointXY(float(x), float(y)) for x, y in coordinates
                ]))
                qgs_feature.setAttribute("Order", int(properties.get("strord", 1)))
                qgs_feature.setAttribute("UpArea", float(properties.get("uparea", 0.0)))
                qgs_feature.setAttribute("idx", int(properties.get("idx", -1)))
                qgs_feature.setAttribute("idx_ds", int(properties.get("idx_ds", -1)))
                qgs_feature.setAttribute("pit", int(bool(properties.get("pit", False))))
                writer.addFeature(qgs_feature)
                written += 1

            del writer

            if written > 0 and os.path.exists(self.channel_network_vector_path):
                self.mark_output_prepared(
                    self.channel_network_vector_path,
                    name="2_Channel_Network",
                    loaded=False
                )
                self.load_layer(self.channel_network_vector_path,
                                "2_Channel_Network", is_raster=False)
                self.log_message(
                    f"Channel Network processing completed with {written} stream segments "
                    f"(threshold: {threshold_cells} cells).")
            else:
                self.log_message("Channel Network processing failed.")
                self.channel_network_vector_path = None
