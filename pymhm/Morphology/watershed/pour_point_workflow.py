# -*- coding: utf-8 -*-
"""Pour-point snap workflow orchestration."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    QgsVectorLayer,
)
from ..core.layer_preparation import LayerPreparationMixin
from .channel_network import ChannelNetworkMixin
from .network_snapper import NetworkSnapperMixin


class PourPointWorkflowMixin(
        LayerPreparationMixin,
        ChannelNetworkMixin,
        NetworkSnapperMixin):
    """Pour-point snap workflow orchestration."""

    def snap_points(self) -> None:
        """Step 3: Snap Pour Points to the nearest high-order channel segment."""
        self.log_message(
            "\n--- Starting Geometry Step 3: Snap Pour Points ---")
        if not self.check_prerequisites(needs_pour_points=True):
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        self.snapped_points_path = os.path.join(
            geometry_folder, "2_pour_points_snapped.shp")

        # Check if snapped points already exist (from load_project_state or previous run)
        if self.snapped_points_path and os.path.exists(self.snapped_points_path):
            self.log_message(
                "Snapped points already exist. Loading existing file...")
            self.load_layer(self.snapped_points_path,
                            "2_pour_points_snapped", is_raster=False)
            return

        if not self.channel_network_vector_path or not os.path.exists(self.channel_network_vector_path):
            if not self._ensure_channel_network(
                    self.process_channel_network,
                    self.process_flow_accumulation,
                    self.fill_dem):
                return

        pour_points_layer = self.dialog.mMapLayerComboBox_pour_points.currentLayer()
        temp_reprojected_path = None

        # Check and reproject pour points layer if needed
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            self.log_message(
                "WARNING: Input CRS is not valid. Using pour points layer CRS as-is.")
        else:
            pour_points_crs = pour_points_layer.crs()
            if pour_points_crs.isValid():
                if pour_points_crs.authid() != input_crs.authid():
                    self.log_message(
                        f"Pour points CRS ({pour_points_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                    temp_reprojected_path = os.path.join(
                        geometry_folder, "2_pour_points_reprojected.shp")

                    reprojected_layer = self.reproject_vector_layer(
                        pour_points_layer, input_crs, temp_reprojected_path)
                    if reprojected_layer:
                        pour_points_layer = reprojected_layer
                        self.log_message(
                            f"Pour points reprojected successfully to {input_crs.authid()}")
                    else:
                        self.log_message(
                            "WARNING: Reprojection failed. Using original pour points layer.")
                        temp_reprojected_path = None  # Don't try to delete if reprojection failed
                else:
                    self.log_message(
                        f"Pour points CRS ({pour_points_crs.authid()}) matches input CRS. No reprojection needed.")
            else:
                self.log_message(
                    "WARNING: Pour points CRS is not valid. Using as-is.")

        channel_network_layer = QgsVectorLayer(
            self.channel_network_vector_path, "Channel Network", "ogr")
        if not channel_network_layer.isValid():
            self.log_message(
                f"ERROR: Could not load channel network layer from {self.channel_network_vector_path}")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    os.remove(temp_reprojected_path)
                except:
                    pass
            return

        # Call the new custom snapping function
        result_path = self.snap_points_to_network(
            pour_points_layer=pour_points_layer,
            channel_network_layer=channel_network_layer,
            output_path=self.snapped_points_path,
            order_field_name='Order',
            high_order_distance=1000.0,
            max_snap_distance=5000.0
        )

        # Clean up temporary reprojected file if it was created
        if temp_reprojected_path and os.path.exists(temp_reprojected_path):
            try:
                os.remove(temp_reprojected_path)
                self.log_message(
                    "Cleaned up temporary reprojected pour points file.")
            except Exception as e:
                self.log_message(
                    f"WARNING: Could not delete temporary file: {e}")

        if result_path:
            self.mark_output_prepared(
                result_path,
                name="2_pour_points_snapped",
                loaded=False
            )
            self.load_layer(result_path,
                            "2_pour_points_snapped", is_raster=False)
        else:
            self.log_message("ERROR: Custom snap points procedure failed.")
            self.snapped_points_path = None
