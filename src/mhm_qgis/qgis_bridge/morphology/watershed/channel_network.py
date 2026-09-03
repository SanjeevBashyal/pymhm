# -*- coding: utf-8 -*-
"""Channel network extraction from flow accumulation."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    QgsFeature,
    QgsGeometry,
    QgsFields,
    QgsWkbTypes,
    QgsPointXY,
)
from ...layers.compat import create_vector_file_writer, qgs_field
from ..core.vector_io import VectorIOMixin
from ..hydrology.flow import FlowAnalysisMixin
from ... import layers


class ChannelNetworkMixin(FlowAnalysisMixin, VectorIOMixin):
    """Channel network extraction from flow accumulation."""

    def process_channel_network(self) -> None:
        """Create a stream network vector from pyflwdir flow accumulation."""
        geometry_folder = project_geometry_folder(self.session.project_folder)
        output_path = os.path.join(geometry_folder, "2_channel_network.shp")
        self.channel_network_vector_path = self.create_channel_network(
            threshold_cells=None,
            output_path=output_path,
            force=False,
            load=True,
        )

    def create_channel_network(
            self,
            threshold_cells: int | None,
            output_path: str,
            force: bool = False,
            load: bool = True) -> str | None:
        """Create a channel network for an explicit cell threshold."""
        if not self._ensure_filled_dem(self.fill_dem):
            return None

        # Keep the UI workflow consistent with the notebook: filled DEM,
        # flow accumulation, then stream extraction.
        if not self._ensure_flow_accumulation(
                self.process_flow_accumulation,
                self.fill_dem):
            return None
        if not output_path:
            self.log_message("ERROR: Channel network output path is required.")
            return None
        output_path = os.path.abspath(output_path)

        # Check if Channel Network already exists, otherwise process it
        if os.path.exists(output_path) and not force:
            if self._vector_crs_matches_raster(
                    output_path, self.filled_dem_path):
                self.log_message(
                    "Channel Network already exists. Loading existing file...")
                if load:
                    self.load_layer(
                        output_path, "2_Channel_Network", is_raster=False)
                return output_path
            self.log_message(
                "Existing Channel Network CRS does not match the filled DEM. Recreating it.")
        if os.path.exists(output_path):
            self._remove_stale_vector_output(output_path)

        self.log_message("Processing Channel Network with pyflwdir...")

        context = self._build_flwdir_from_filled_dem()
        if not context:
            self.log_message("Channel Network processing failed.")
            return None

        deps = context["deps"]
        np = deps["np"]
        flwdir = context["flwdir"]
        reference = context["reference"]
        invalid_mask = context["invalid"]

        flow_accumulation = flwdir.upstream_area(unit="cell")
        flow_accumulation[invalid_mask] = 0

        valid_accumulation = flow_accumulation[flow_accumulation > 0]
        if valid_accumulation.size == 0:
            self.log_message("ERROR: No valid flow accumulation cells found.")
            return None

        max_accumulation = int(np.nanmax(valid_accumulation))
        if threshold_cells is None:
            cell_area_m2 = self._reference_cell_area_m2(reference)
            threshold = max(1, int(round(10000000.0 / cell_area_m2)))
            if threshold > max_accumulation:
                threshold = max(
                    1, int(np.nanpercentile(valid_accumulation, 95)))
                self.log_message(
                    "Channel threshold exceeded the basin accumulation. "
                    f"Using 95th percentile threshold: {threshold} cells.")
        else:
            try:
                threshold = int(threshold_cells)
            except (TypeError, ValueError):
                self.log_message(
                    "ERROR: Channel threshold must be a positive integer.")
                return None
            if threshold < 1:
                self.log_message(
                    "ERROR: Channel threshold must be a positive integer.")
                return None
            if threshold > max_accumulation:
                self.log_message(
                    "ERROR: Channel threshold exceeds the maximum flow "
                    f"accumulation ({max_accumulation} cells).")
                return None

        stream_mask = flow_accumulation >= threshold
        stream_order = flwdir.stream_order("strahler", mask=stream_mask)
        features = flwdir.streams(
            mask=stream_mask,
            strord=stream_order,
            uparea=flow_accumulation
        )

        if not features:
            self.log_message("ERROR: pyflwdir did not produce any stream features.")
            return None

        fields = QgsFields()
        fields.append(qgs_field("Order", "Int"))
        fields.append(qgs_field("UpArea", "Double"))
        fields.append(qgs_field("idx", "Int"))
        fields.append(qgs_field("idx_ds", "Int"))
        fields.append(qgs_field("pit", "Int"))

        output_crs = layers.crs_of(self.filled_dem_path) or self.session.crs

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        self._remove_vector_dataset(output_path)
        writer = create_vector_file_writer(
            output_path,
            fields,
            QgsWkbTypes.LineString,
            output_crs,
        )
        if writer.hasError():
            self.log_message(f"ERROR creating channel network: {writer.errorMessage()}")
            return None

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

        if written > 0 and os.path.exists(output_path):
            self.mark_output_prepared(
                output_path,
                name="2_Channel_Network",
                loaded=False
            )
            if load:
                self.load_layer(
                    output_path, "2_Channel_Network", is_raster=False)
            self.log_message(
                f"Channel Network processing completed with {written} stream segments "
                f"(threshold: {threshold} cells).")
            return output_path

        self.log_message("Channel Network processing failed.")
        return None
