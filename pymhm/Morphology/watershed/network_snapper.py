# -*- coding: utf-8 -*-
"""Custom high-order stream-network snapping implementation."""
from __future__ import annotations

from typing import Any

from ..common import (
    os,
    QMessageBox,
    QgsVectorLayer,
    QgsFeature,
    QgsGeometry,
    QgsSpatialIndex,
    QgsFields,
    QgsWkbTypes,
    create_vector_file_writer,
    qgs_field,
    processing,
)
from ..core.vector_io import VectorIOMixin


class NetworkSnapperMixin(VectorIOMixin):
    """Custom high-order stream-network snapping implementation."""

    def snap_points_to_network(
            self,
            pour_points_layer: Any,
            channel_network_layer: Any,
            output_path: str,
            order_field_name: str = 'Order',
            high_order_distance: float = 1000.0,
            max_snap_distance: float = 5000.0) -> str | None:
        """
        Snaps pour points to a channel network based on stream order and proximity.

        For each pour point, it first searches within `high_order_distance` for the
        channel segment with the highest order. If found, it snaps to the closest
        point on that segment. If multiple segments share the highest order, it
        chooses the closest one.

        If no channels are found in the first search, it expands the search to
        `max_snap_distance` and snaps to the geometrically closest channel segment,
        regardless of order.

        Args:
            pour_points_layer (QgsVectorLayer): Input layer of pour points.
            channel_network_layer (QgsVectorLayer): Input layer of the channel network.
            output_path (str): The file path for the output snapped points layer.
            order_field_name (str): The name of the stream order attribute.
            high_order_distance (float): Search radius for the high-order snap.
            max_snap_distance (float): Maximum search radius for the fallback snap.

        Returns:
            The output file path if successful, otherwise None.
        """
        self.log_message("--- Starting custom snap points procedure ---")

        # 1. Validate inputs and check for order field
        if order_field_name not in channel_network_layer.fields().names():
            self.log_message(
                f"ERROR: Order field '{order_field_name}' not found in channel network layer.")
            QMessageBox.critical(
                self.dialog, "Missing Field", f"Could not find order field '{order_field_name}' in the channel network layer.")
            return None

        # 2. Prepare the output layer
        source_fields = pour_points_layer.fields()
        output_fields = QgsFields()
        for field in source_fields:
            output_fields.append(field)
        output_fields.append(qgs_field("snap_status", "String"))
        output_fields.append(qgs_field("snap_dist", "Double"))
        output_fields.append(qgs_field("snapped_order", "Int"))

        # If the file already exists, remove it before creating a new one
        if os.path.exists(output_path):
            self._remove_vector_dataset(output_path)

        writer = create_vector_file_writer(
            output_path,
            output_fields,
            QgsWkbTypes.Point,
            self.dialog.get_crs(),
        )
        if writer.hasError():
            self.log_message(
                f"ERROR creating output file: {writer.errorMessage()}")
            return None

        # 3. Create a spatial index for fast searching
        self.log_message("Creating spatial index for channel network...")
        network_index = QgsSpatialIndex(channel_network_layer.getFeatures())

        # 4. Process each pour point
        self.log_message("Snapping pour points...")
        total_points = pour_points_layer.featureCount()
        for i, point_feat in enumerate(pour_points_layer.getFeatures()):
            if (i + 1) % 20 == 0:
                self.log_message(
                    f"  ...processing point {i + 1} of {total_points}")

            original_geom = point_feat.geometry()
            original_point = original_geom.asPoint()
            new_feat = QgsFeature(output_fields)
            new_feat.setAttributes(point_feat.attributes())
            snapped = False

            # --- Stage 1: High-order snap search ---
            search_rect = original_geom.buffer(
                high_order_distance, 5).boundingBox()
            candidate_ids = network_index.intersects(search_rect)

            best_feat_id_high_order = -1
            max_order = -1
            min_dist_at_max_order = float('inf')

            if candidate_ids:
                candidate_features = {
                    f.id(): f for f in channel_network_layer.getFeatures(candidate_ids)}
                for feat_id in candidate_ids:
                    feat = candidate_features.get(feat_id)
                    if not feat:
                        continue
                    dist = feat.geometry().distance(original_geom)
                    if dist > high_order_distance:
                        continue

                    order_val = feat.attribute(order_field_name)

                    if order_val > max_order:
                        max_order = order_val
                        min_dist_at_max_order = dist
                        best_feat_id_high_order = feat_id
                    elif order_val == max_order and dist < min_dist_at_max_order:
                        min_dist_at_max_order = dist
                        best_feat_id_high_order = feat_id

            if best_feat_id_high_order != -1:
                target_feat = candidate_features[best_feat_id_high_order]
                closest_point = target_feat.geometry().closestSegmentWithContext(original_point)
                # closestSegmentWithContext returns (sqrDist, closestPoint, afterVertex, sqrDistToSegment)
                new_geom = QgsGeometry.fromPointXY(closest_point[1])

                new_feat.setGeometry(new_geom)
                new_feat.setAttribute("snap_status", "high_order")
                new_feat.setAttribute(
                    "snap_dist", new_geom.distance(original_geom))
                new_feat.setAttribute("snapped_order", max_order)
                writer.addFeature(new_feat)
                snapped = True

            # --- Stage 2: Fallback to closest snap if high-order search failed ---
            if not snapped:
                nearest_ids = network_index.nearestNeighbor(original_point, 1)
                if nearest_ids:
                    target_feat = channel_network_layer.getFeature(
                        nearest_ids[0])
                    dist_to_geom = target_feat.geometry().distance(original_geom)

                    if dist_to_geom <= max_snap_distance:
                        closest_point = target_feat.geometry().closestSegmentWithContext(original_point)
                        new_geom = QgsGeometry.fromPointXY(closest_point[1])

                        new_feat.setGeometry(new_geom)
                        new_feat.setAttribute("snap_status", "closest")
                        new_feat.setAttribute("snap_dist", dist_to_geom)
                        order_val = target_feat.attribute(
                            order_field_name)
                        new_feat.setAttribute("snapped_order", order_val)
                        writer.addFeature(new_feat)
                        snapped = True

            # --- Stage 3: If no snap possible, write original point with 'failed' status ---
            if not snapped:
                new_feat.setGeometry(original_geom)
                new_feat.setAttribute("snap_status", "failed")
                new_feat.setAttribute("snap_dist", 0.0)
                writer.addFeature(new_feat)

        del writer  # Finalize writing to the file
        self.log_message(f"Snapping complete. Output saved to {output_path}")
        self.mark_output_prepared(
            output_path,
            name=os.path.basename(output_path),
            loaded=False
        )
        return output_path

    # --- Hydrological Processing Methods ---
