# -*- coding: utf-8 -*-
"""
Morphology/Geometry processing module for pymhm
Contains all geometry-related processing methods
"""
import os
from qgis.PyQt.QtWidgets import QMessageBox, QFileDialog
from qgis.core import (
    QgsVectorLayer,
    QgsRasterLayer,
    QgsApplication,
    QgsFeature,
    QgsGeometry,
    QgsSpatialIndex,
    QgsVectorFileWriter,
    QgsField,
    QgsFields,
    QgsWkbTypes,
    QgsProject,
    QgsCoordinateReferenceSystem
)
from qgis.PyQt.QtCore import QVariant, NULL
import processing
from .utils import DialogUtils


class MorphologyProcessor(DialogUtils):
    """
    Handles all morphology/geometry processing functionality.
    Expects the dialog instance to have DialogUtils attributes.
    """
    
    def __init__(self, dialog):
        """
        Initialize the morphology processor with reference to the dialog.
        
        Args:
            dialog: The pymhmDialog instance (needs DialogUtils methods)
        """
        self.dialog = dialog
        # Copy DialogUtils methods to self
        self.log_message = dialog.log_message
        self.check_prerequisites = dialog.check_prerequisites
        self.run_processing_algorithm = dialog.run_processing_algorithm
        self.load_layer = dialog.load_layer
        self.get_dem_extent_and_resolution = dialog.get_dem_extent_and_resolution
        
        # DEM layer reference (reprojected if needed)
        self.dem_layer = None  # Will store the reprojected DEM layer if CRS differs
        
        # Processed layer references
        self.land_use_layer = None  # Final land use layer with reclassified values (1, 2, or 3)
        
        # Paths for geometry processing outputs
        self.filled_dem_path = None
        self.flow_dir_path = None
        self.channel_network_vector_path = None
        self.snapped_points_path = None
        self.watershed_raster_path = None
        self.watershed_vector_path = None
        
        # Hydrological processing paths
        self.aspect_path = None
        self.slope_path = None
        self.flow_direction_path = None
        self.flow_accumulation_path = None
        self.gauge_position_path = None

    def load_project_state(self):
        """Checks for existing output files and sets instance attributes to resume work."""
        self.log_message("\n--- Checking for existing project files... ---")

        files_to_check = {
            # Original processing files
            'filled_dem_path': "1_dem_filled.sdat",
            'flow_dir_path': "2_flow_direction.sdat",
            'channel_network_vector_path': "2_channel_network.shp",
            'snapped_points_path': "2_pour_points_snapped.shp",
            'watershed_raster_path': "1_watershed_raster.sdat",
            'watershed_vector_path': "1_watershed_final.shp",

            # New hydrological processing files
            'aspect_path': "1_dem_aspect.sdat",
            'slope_path': "1_dem_slope.sdat",
            'flow_direction_path': "2_flow_direction.sdat",
            'flow_accumulation_path': "2_flow_accumulation.sdat",
            'gauge_position_path': "2_gauge_position.sdat"
        }

        found_any = False
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        for attr, filename in files_to_check.items():
            expected_path = os.path.join(geometry_folder, filename)
            if os.path.exists(expected_path):
                setattr(self, attr, expected_path)
                self.log_message(f"Found existing file: {filename}")
                found_any = True
            else:
                setattr(self, attr, None)

        # Check for final land use layer separately
        land_use_path = os.path.join(geometry_folder, "3_land_use.tif")
        if os.path.exists(land_use_path):
            self.land_use_layer = land_use_path
            self.log_message(f"Found existing file: 3_land_use.tif")
            found_any = True
        
        # Check for reprojected DEM and load it into self.dem_layer (but don't add to project)
        reprojected_dem_path = os.path.join(geometry_folder, "0_dem_reprojected.tif")
        if os.path.exists(reprojected_dem_path):
            reprojected_layer = QgsRasterLayer(reprojected_dem_path, "1_DEM_Reprojected")
            if reprojected_layer.isValid():
                self.dem_layer = reprojected_layer
                self.log_message(f"Found existing reprojected DEM: 0_dem_reprojected.tif")
                found_any = True
            else:
                self.log_message(f"WARNING: Reprojected DEM file exists but is not valid.")
        
        # Check for filled DEM and load it to QGIS interface
        if self.filled_dem_path and os.path.exists(self.filled_dem_path):
            self.log_message(f"Loading existing filled DEM to QGIS interface...")
            self.load_layer(self.filled_dem_path, "1_DEM_Filled")
        
        if found_any:
            self.log_message(
                "Project state loaded. You may skip completed steps.")
        else:
            self.log_message(
                "No existing project files found in 'Geometry' folder.")

    # --- Geometry Processing Methods ---

    def check_and_reproject_dem(self, dem_layer):
        """
        Check if DEM CRS matches input CRS, and reproject if needed.
        Stores the reprojected layer in self.dem_layer attribute.
        
        Args:
            dem_layer: The DEM raster layer
            
        Returns:
            The DEM layer (original or reprojected) to use for processing
        """
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            self.log_message("WARNING: Input CRS is not valid. Using DEM CRS as-is.")
            self.dem_layer = dem_layer
            return dem_layer
        
        dem_crs = dem_layer.crs()
        if not dem_crs.isValid():
            self.log_message("WARNING: DEM CRS is not valid. Using DEM as-is.")
            self.dem_layer = dem_layer
            return dem_layer
        
        # Check if CRS matches
        if dem_crs.authid() == input_crs.authid():
            self.log_message(f"DEM CRS ({dem_crs.authid()}) matches input CRS. No reprojection needed.")
            self.dem_layer = dem_layer
            return dem_layer
        
        # Need to reproject - check if reprojected file already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        reprojected_dem_path = os.path.join(geometry_folder, "0_dem_reprojected.tif")
        
        if os.path.exists(reprojected_dem_path):
            self.log_message(f"Found existing reprojected DEM. Loading it...")
            reprojected_layer = QgsRasterLayer(reprojected_dem_path, "1_DEM_Reprojected")
            if reprojected_layer.isValid():
                self.dem_layer = reprojected_layer
                return reprojected_layer
            else:
                self.log_message("WARNING: Existing reprojected DEM is not valid. Reprojecting again...")
        
        # Reproject DEM
        self.log_message(f"DEM CRS ({dem_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
        
        params_warp = {
            'INPUT': dem_layer.source(),
            'SOURCE_CRS': None,  # Auto-detect from input
            'TARGET_CRS': input_crs,
            'RESAMPLING': 0,  # Nearest neighbor
            'NODATA': None,
            'TARGET_RESOLUTION': None,
            'OPTIONS': None,
            'DATA_TYPE': 0,  # Use input data type
            'TARGET_EXTENT': None,
            'TARGET_EXTENT_CRS': None,
            'MULTITHREADING': False,
            'EXTRA': '',
            'OUTPUT': reprojected_dem_path
        }
        
        result = self.run_processing_algorithm("gdal:warpreproject", params_warp)
        if result and os.path.exists(reprojected_dem_path):
            reprojected_layer = QgsRasterLayer(reprojected_dem_path, "1_DEM_Reprojected")
            if reprojected_layer.isValid():
                self.log_message(f"DEM reprojected successfully to {input_crs.authid()}")
                # Do not add reprojected layer to QGIS project - it's only used for processing
                self.dem_layer = reprojected_layer
                return reprojected_layer
            else:
                self.log_message("ERROR: Reprojected DEM layer is not valid. Using original DEM.")
                self.dem_layer = dem_layer
                return dem_layer
        else:
            self.log_message("ERROR: DEM reprojection failed. Using original DEM.")
            self.dem_layer = dem_layer
            return dem_layer

    def fill_dem(self):
        """Step 1: Fill Sinks in DEM using SAGA NG's Wang & Liu algorithm."""
        self.log_message("\n--- Starting Geometry Step 1: Fill DEM ---")
        if not self.check_prerequisites():
            return

        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.filled_dem_path = os.path.join(
            geometry_folder, "1_dem_filled.sdat")
        
        # Check if filled DEM already exists
        if self.filled_dem_path and os.path.exists(self.filled_dem_path):
            self.log_message("Filled DEM already exists. Loading existing file...")
            self.load_layer(self.filled_dem_path, "1_DEM_Filled")
            return

        original_dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        self.log_message(f"Input DEM: {original_dem_layer.name()}")

        # Check and reproject DEM if needed (stores in self.dem_layer)
        dem_layer = self.check_and_reproject_dem(original_dem_layer)
        
        # Ensure self.dem_layer is set
        if not self.dem_layer:
            self.dem_layer = dem_layer

        params = {
            'ELEV': self.dem_layer,  # Use reprojected DEM layer
            'MINSLOPE': 0.01,
            'FILLED': self.filled_dem_path,
            'FDIR': 'TEMPORARY_OUTPUT',
            'WSHED': 'TEMPORARY_OUTPUT'
        }
        result = self.run_processing_algorithm(
            "sagang:fillsinkswangliu", params)
        if result:
            self.load_layer(result['FILLED'], "1_DEM_Filled")
        else:
            self.filled_dem_path = None

    def process_channel_network(self):
        """Process Flow Direction and Channel Network."""
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 1 (Fill DEM) successfully first.")
            return

        # Ensure flow accumulation is processed first (required for channel network)
        self.process_flow_accumulation()

        # Check if flow accumulation exists (required dependency)
        if not self.flow_accumulation_path or not os.path.exists(self.flow_accumulation_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Flow accumulation processing failed or file not found.")
            return

        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.channel_network_vector_path = os.path.join(
            geometry_folder, "2_channel_network.shp")

        # Check if Channel Network already exists, otherwise process it
        if self.channel_network_vector_path and os.path.exists(self.channel_network_vector_path):
            self.log_message("Channel Network already exists. Loading existing file...")
            self.load_layer(self.channel_network_vector_path,
                            "2_Channel_Network", is_raster=False)
        else:
            self.log_message("Processing Channel Network...")
            params_chnl = {
                'ELEVATION': self.filled_dem_path,
                'SINKROUTE': None,
                'CHNLNTWRK': 'TEMPORARY_OUTPUT',
                'CHNLROUTE': 'TEMPORARY_OUTPUT',
                'SHAPES': self.channel_network_vector_path,
                'INIT_GRID': self.flow_accumulation_path,
                'INIT_METHOD': 2,       # Greater than
                # Threshold for channel initiation (cell count)
                'INIT_VALUE': 1000000.0,
                'DIV_GRID': None,
                'DIV_CELLS': 5.0,
                'TRACE_WEIGHT': None,
                'MINLEN': 10.0
            }
            result_chnl = self.run_processing_algorithm(
                "sagang:channelnetwork", params_chnl)
            if result_chnl:
                self.load_layer(self.channel_network_vector_path,
                                "2_Channel_Network", is_raster=False)
                self.log_message("Channel Network processing completed successfully.")
            else:
                self.log_message("Channel Network processing failed.")
                self.channel_network_vector_path = None

    def snap_points(self):
        """Step 3: Snap Pour Points to the nearest high-order channel segment."""
        self.log_message(
            "\n--- Starting Geometry Step 3: Snap Pour Points ---")
        if not self.check_prerequisites(needs_pour_points=True):
            return
        
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.snapped_points_path = os.path.join(
            geometry_folder, "2_pour_points_snapped.shp")
        
        # Check if snapped points already exist (from load_project_state or previous run)
        if self.snapped_points_path and os.path.exists(self.snapped_points_path):
            self.log_message("Snapped points already exist. Loading existing file...")
            self.load_layer(self.snapped_points_path, "2_pour_points_snapped", is_raster=False)
            return
        
        if not self.channel_network_vector_path or not os.path.exists(self.channel_network_vector_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 2 (Create Network) successfully first.")
            return

        pour_points_layer = self.dialog.mMapLayerComboBox_pour_points.currentLayer()
        temp_reprojected_path = None
        
        # Check and reproject pour points layer if needed
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            self.log_message("WARNING: Input CRS is not valid. Using pour points layer CRS as-is.")
        else:
            pour_points_crs = pour_points_layer.crs()
            if pour_points_crs.isValid():
                if pour_points_crs.authid() != input_crs.authid():
                    self.log_message(f"Pour points CRS ({pour_points_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                    temp_reprojected_path = os.path.join(geometry_folder, "2_pour_points_reprojected.shp")
                    
                    reprojected_layer = self.reproject_vector_layer(pour_points_layer, input_crs, temp_reprojected_path)
                    if reprojected_layer:
                        pour_points_layer = reprojected_layer
                        self.log_message(f"Pour points reprojected successfully to {input_crs.authid()}")
                    else:
                        self.log_message("WARNING: Reprojection failed. Using original pour points layer.")
                        temp_reprojected_path = None  # Don't try to delete if reprojection failed
                else:
                    self.log_message(f"Pour points CRS ({pour_points_crs.authid()}) matches input CRS. No reprojection needed.")
            else:
                self.log_message("WARNING: Pour points CRS is not valid. Using as-is.")
        
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
            order_field_name='Order',  # SAGA Channel Network tool typically uses 'ORDER'
            high_order_distance=1000.0,
            max_snap_distance=5000.0
        )

        # Clean up temporary reprojected file if it was created
        if temp_reprojected_path and os.path.exists(temp_reprojected_path):
            try:
                os.remove(temp_reprojected_path)
                self.log_message("Cleaned up temporary reprojected pour points file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary file: {e}")

        if result_path:
            self.load_layer(result_path,
                            "2_pour_points_snapped", is_raster=False)
        else:
            self.log_message("ERROR: Custom snap points procedure failed.")
            self.snapped_points_path = None

    def delineate_watershed(self):
        """Step 4: Delineate watershed using the snapped points and flow direction."""
        self.log_message(
            "\n--- Starting Geometry Step 4: Delineate Watershed ---")
        if not self.snapped_points_path or not os.path.exists(self.snapped_points_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 3 (Snap Points) successfully first.")
            return
        if not self.flow_dir_path or not os.path.exists(self.flow_dir_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 2 (Create Network) successfully first.")
            return

        # Sub-step 4a: Delineate watersheds for each snapped point
        self.log_message("Delineating watersheds for each snapped point...")

        # Load the snapped points layer
        snapped_points_layer = QgsVectorLayer(
            self.snapped_points_path, "Snapped Points", "ogr")
        if not snapped_points_layer.isValid():
            QMessageBox.warning(
                self.dialog, "Error", "Could not load snapped points layer.")
            return

        # Get all features from snapped points
        features = list(snapped_points_layer.getFeatures())
        if not features:
            QMessageBox.warning(self.dialog, "Error", "No snapped points found.")
            return

        # Create a list to store all watershed outputs
        watershed_outputs = []

        # Loop over each snapped point
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        for i, feature in enumerate(features):
            # Get the geometry and attributes
            geom = feature.geometry()
            point = geom.asPoint()

            # Get the Name attribute (use index if Name is not available)
            name_attr = feature.attribute("Name")
            if not name_attr or name_attr == NULL:
                name_attr = f"Watershed_{i+1}"

            # Clean the name for filename (remove invalid characters)
            clean_name = "".join(c for c in str(
                name_attr) if c.isalnum() or c in (' ', '-', '_')).rstrip()
            clean_name = clean_name.replace(' ', '_')

            self.log_message(f"Processing watershed for point: {name_attr}")

            # Create output path for this watershed
            watershed_raster_path = os.path.join(
                geometry_folder, f"4_watershed_{clean_name}.sdat")

            # Set up parameters for this specific point
            params_ws = {
                'TARGET': None,
                'TARGET_PT_X': point.x(),
                'TARGET_PT_Y': point.y(),
                'ELEVATION': self.filled_dem_path,
                'SINKROUTE': None,
                'AREA': watershed_raster_path,
                'METHOD': 0,
                'CONVERGE': 1.1,
                'MFD_CONTOUR': False
            }

            # Run watershed delineation for this point
            result_ws = self.run_processing_algorithm(
                "sagang:upslopearea", params_ws)

            if result_ws:
                watershed_outputs.append({
                    'name': name_attr,
                    'clean_name': clean_name,
                    'raster_path': watershed_raster_path,
                    'point': point
                })
                self.load_layer(watershed_raster_path,
                                f"4_Watershed_{clean_name}")
            else:
                self.log_message(
                    f"Failed to create watershed for point: {name_attr}")

        if not watershed_outputs:
            QMessageBox.warning(
                self.dialog, "Error", "No watersheds were successfully created.")
            return

        # Sub-step 4b: Combine all watersheds into a single vector layer
        self.log_message("Combining all watersheds into final vector layer...")
        self.watershed_vector_path = os.path.join(
            geometry_folder, "1_watershed_final.shp")

        # Create a combined watershed vector layer
        self.create_combined_watershed_layer(
            watershed_outputs, self.watershed_vector_path)

        if os.path.exists(self.watershed_vector_path):
            self.load_layer(self.watershed_vector_path,
                            "1_watershed_final", is_raster=False)
        else:
            self.watershed_vector_path = None

    def snap_points_to_network(self, pour_points_layer, channel_network_layer, output_path,
                               order_field_name='Order', high_order_distance=1000.0, max_snap_distance=5000.0):
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
        output_fields.append(QgsField("snap_status", QVariant.String))
        output_fields.append(QgsField("snap_dist", QVariant.Double))
        output_fields.append(QgsField("snapped_order", QVariant.Int))

        # If the file already exists, remove it before creating a new one
        if os.path.exists(output_path):
            try:
                os.remove(output_path)
            except Exception as e:
                self.log_message(f"ERROR deleting existing output file: {e}")
                return None

        writer = QgsVectorFileWriter(output_path, "UTF-8", output_fields,
                                     QgsWkbTypes.Point, self.dialog.get_crs(), "ESRI Shapefile")
        if writer.hasError() != QgsVectorFileWriter.NoError:
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
        return output_path

    def create_combined_watershed_layer(self, watershed_outputs, output_path):
        """
        Create a combined vector layer from multiple watershed rasters.

        Args:
            watershed_outputs: List of dictionaries containing watershed information
            output_path: Path for the output vector layer
        """
        try:
            # Create a new vector layer for the combined watersheds
            fields = QgsFields()
            fields.append(QgsField("id", QVariant.Int))
            fields.append(QgsField("name", QVariant.String))
            fields.append(QgsField("clean_name", QVariant.String))
            fields.append(QgsField("area_km2", QVariant.Double))
            fields.append(QgsField("perimeter_km", QVariant.Double))

            # Create the writer
            writer = QgsVectorFileWriter(
                output_path,
                "UTF-8",
                fields,
                QgsWkbTypes.Polygon,
                QgsProject.instance().crs(),
                "GPKG"
            )

            if writer.hasError() != QgsVectorFileWriter.NoError:
                self.log_message(
                    f"Error creating vector file: {writer.errorMessage()}")
                return False

            # Process each watershed
            for i, watershed_info in enumerate(watershed_outputs):
                raster_path = watershed_info['raster_path']

                if not os.path.exists(raster_path):
                    self.log_message(f"Raster file not found: {raster_path}")
                    continue

                # Convert raster to vector for this watershed
                temp_vector_path = os.path.join(
                    os.path.dirname(output_path),
                    f"temp_watershed_{watershed_info['clean_name']}.gpkg"
                )

                params_poly = {
                    'INPUT': raster_path,
                    'BAND': 1,
                    'FIELD': 'DN',
                    'EIGHT_CONNECTEDNESS': False,
                    'OUTPUT': temp_vector_path
                }

                result_poly = self.run_processing_algorithm(
                    "gdal:polygonize", params_poly)

                if result_poly and os.path.exists(temp_vector_path):
                    # Load the temporary vector layer
                    temp_layer = QgsVectorLayer(
                        temp_vector_path, "Temp Watershed", "ogr")

                    if temp_layer.isValid():
                        # Get the largest polygon (main watershed)
                        max_area = 0
                        main_feature = None

                        for feature in temp_layer.getFeatures():
                            geom = feature.geometry()
                            area = geom.area()
                            if area > max_area:
                                max_area = area
                                main_feature = feature

                        if main_feature:
                            # Create new feature with attributes
                            new_feature = QgsFeature()
                            new_feature.setGeometry(main_feature.geometry())
                            new_feature.setAttribute("id", i + 1)
                            new_feature.setAttribute(
                                "name", watershed_info['name'])
                            new_feature.setAttribute(
                                "clean_name", watershed_info['clean_name'])

                            # Calculate area in km²
                            area_km2 = main_feature.geometry().area() / 1000000  # Convert from m² to km²
                            new_feature.setAttribute(
                                "area_km2", round(area_km2, 2))

                            # Calculate perimeter in km
                            perimeter_km = main_feature.geometry().length() / 1000  # Convert from m to km
                            new_feature.setAttribute(
                                "perimeter_km", round(perimeter_km, 2))

                            writer.addFeature(new_feature)
                            self.log_message(
                                f"Added watershed: {watershed_info['name']} (Area: {area_km2:.2f} km²)")

                    # Clean up temporary file
                    try:
                        os.remove(temp_vector_path)
                    except:
                        pass

            del writer
            self.log_message(
                f"Combined watershed layer saved to: {output_path}")
            return True

        except Exception as e:
            self.log_message(
                f"Error creating combined watershed layer: {str(e)}")
            return False

    # --- Hydrological Processing Methods ---

    def process_aspect(self):
        """Process Aspect from input DEM"""
        # Ensure DEM is checked and reprojected if needed
        if not self.dem_layer:
            input_dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
            if not input_dem_layer:
                QMessageBox.warning(self.dialog, "Input Error",
                                    "Please select an input DEM layer.")
                return
            self.check_and_reproject_dem(input_dem_layer)
        
        # Use reprojected DEM layer
        dem_layer = self.dem_layer if self.dem_layer else self.dialog.mMapLayerComboBox_dem.currentLayer()
        input_dem_path = dem_layer.source()
        
        if not os.path.exists(input_dem_path):
            QMessageBox.warning(self.dialog, "Input Error",
                                "Input DEM file not found.")
            return

        # Check if aspect already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.aspect_path = os.path.join(
            geometry_folder, "1_dem_aspect.sdat")
        if self.aspect_path and os.path.exists(self.aspect_path):
            self.log_message("Aspect already exists. Loading existing file...")
            self.load_layer(self.aspect_path, "1_DEM_Aspect")
            return

        self.log_message("Processing Aspect...")

        params_aspect = {
            'INPUT': input_dem_path,
            'BAND': 1,
            'TRIG_ANGLE': False,
            'ZERO_FLAT': False,
            'COMPUTE_EDGES': False,
            'ZEVENBERGEN': False,
            'OPTIONS': None,
            'EXTRA': '',
            'OUTPUT': self.aspect_path
        }

        result = self.run_processing_algorithm("gdal:aspect", params_aspect)
        if result:
            self.load_layer(self.aspect_path, "1_DEM_Aspect")
            self.log_message("Aspect processing completed successfully.")
        else:
            self.log_message("Aspect processing failed.")

    def process_slope(self):
        """Process Slope from input DEM"""
        # Ensure DEM is checked and reprojected if needed
        if not self.dem_layer:
            input_dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
            if not input_dem_layer:
                QMessageBox.warning(self.dialog, "Input Error",
                                    "Please select an input DEM layer.")
                return
            self.check_and_reproject_dem(input_dem_layer)
        
        # Use reprojected DEM layer
        dem_layer = self.dem_layer if self.dem_layer else self.dialog.mMapLayerComboBox_dem.currentLayer()
        input_dem_path = dem_layer.source()
        
        if not os.path.exists(input_dem_path):
            QMessageBox.warning(self.dialog, "Input Error",
                                "Input DEM file not found.")
            return

        # Check if slope already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.slope_path = os.path.join(
            geometry_folder, "1_dem_slope.sdat")
        if self.slope_path and os.path.exists(self.slope_path):
            self.log_message("Slope already exists. Loading existing file...")
            self.load_layer(self.slope_path, "1_DEM_Slope")
            return

        self.log_message("Processing Slope...")

        params_slope = {
            'INPUT': input_dem_path,
            'BAND': 1,
            'SCALE': 1,
            'AS_PERCENT': False,
            'COMPUTE_EDGES': False,
            'ZEVENBERGEN': False,
            'OPTIONS': None,
            'EXTRA': '',
            'OUTPUT': self.slope_path
        }

        result = self.run_processing_algorithm("gdal:slope", params_slope)
        if result:
            self.load_layer(self.slope_path, "1_DEM_Slope")
            self.log_message("Slope processing completed successfully.")
        else:
            self.log_message("Slope processing failed.")

    def process_flow_direction(self):
        """Process Flow Direction using Channel Network and Drainage Basins algorithm"""
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 1 (Fill DEM) successfully first.")
            return

        # Check if flow direction already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.flow_direction_path = os.path.join(
            geometry_folder, "2_flow_direction.sdat")
        if self.flow_direction_path and os.path.exists(self.flow_direction_path):
            self.log_message(
                "Flow Direction already exists. Loading existing file...")
            self.load_layer(self.flow_direction_path, "2_Flow_Direction")
            return

        self.log_message("Processing Flow Direction using Channel Network and Drainage Basins...")

        params_flow_dir = {
            'DEM': self.filled_dem_path,
            'DIRECTION': self.flow_direction_path,
            'CONNECTION': 'TEMPORARY_OUTPUT',
            'ORDER': 'TEMPORARY_OUTPUT',
            'BASIN': 'TEMPORARY_OUTPUT',
            'SEGMENTS': 'TEMPORARY_OUTPUT',
            'BASINS': 'TEMPORARY_OUTPUT',
            'NODES': 'TEMPORARY_OUTPUT',
            'THRESHOLD': 5,
            'SUBBASINS': False
        }

        result = self.run_processing_algorithm(
            "sagang:channelnetworkanddrainagebasins", params_flow_dir)
        if result:
            self.load_layer(self.flow_direction_path, "2_Flow_Direction")
            self.log_message(
                "Flow Direction processing completed successfully.")
        else:
            self.log_message("Flow Direction processing failed.")

    def process_flow_accumulation(self):
        """Process Flow Accumulation using D8 method"""
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 1 (Fill DEM) successfully first.")
            return

        # Check if flow accumulation D8 already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.flow_accumulation_path = os.path.join(
            geometry_folder, "2_flow_accumulation.sdat")
        if self.flow_accumulation_path and os.path.exists(self.flow_accumulation_path):
            self.log_message(
                "Flow Accumulation (D8) already exists. Loading existing file...")
            self.load_layer(self.flow_accumulation_path,
                            "2_Flow_Accumulation")
            return

        self.log_message("Processing Flow Accumulation (D8)...")

        params_acc = {
            'DEM': self.filled_dem_path,
            'PREPROCESSING': 1,  # Fill Sinks
            'FLOW_ROUTING': 4,  # D8
            'TCA': self.flow_accumulation_path,
            'SCA': 'TEMPORARY_OUTPUT',
            'FLOW_PATH_LENGTH': 'TEMPORARY_OUTPUT'
        }
        result = self.run_processing_algorithm(
            "sagang:flowaccumulationonestep", params_acc)
        if result:
            self.load_layer(self.flow_accumulation_path,
                            "2_Flow_Accumulation")
            self.log_message(
                "Flow Accumulation (D8) processing completed successfully.")
        else:
            self.log_message("Flow Accumulation (D8) processing failed.")

    def process_gauge_position(self):
        """Process Gauge Position from pour points"""
        # Check if gauge position already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.gauge_position_path = os.path.join(
            geometry_folder, "2_gauge_position.sdat")
        if self.gauge_position_path and os.path.exists(self.gauge_position_path):
            self.log_message(
                "Gauge Position already exists. Loading existing file...")
            self.load_layer(self.gauge_position_path, "2_Gauge_Position")
            return

        self.log_message("Processing Gauge Position...")

        params_gauge = {
            'INPUT': self.snapped_points_path,
            'FIELD': '',
            'BURN': 1,
            'USE_Z': False,
            'UNITS': 1,
            'WIDTH': 0,
            'HEIGHT': 0,
            'EXTENT': None,
            'NODATA': 0,
            'OPTIONS': None,
            'DATA_TYPE': 5,
            'INIT': None,
            'INVERT': False,
            'EXTRA': '',
            'OUTPUT': self.gauge_position_path
        }

        result = self.run_processing_algorithm("gdal:rasterize", params_gauge)
        if result:
            self.load_layer(self.gauge_position_path, "2_Gauge_Position")
            self.log_message(
                "Gauge Position processing completed successfully.")
        else:
            self.log_message("Gauge Position processing failed.")

    def reproject_vector_layer(self, vector_layer, target_crs, output_path):
        """
        Reproject a vector layer to target CRS.
        
        Args:
            vector_layer: QgsVectorLayer to reproject
            target_crs: QgsCoordinateReferenceSystem target CRS
            output_path: Path for reprojected output
            
        Returns:
            QgsVectorLayer: Reprojected layer or None if failed
        """
        params_reproject = {
            'INPUT': vector_layer,
            'SOURCE_CRS': None,  # Auto-detect
            'TARGET_CRS': target_crs,
            'OPERATION': '',  # Use default
            'TARGET_EXTENT': None,
            'TARGET_EXTENT_CRS': None,
            'CONCAT_PREFIX': '',
            'NO_DATA': None,
            'COPY_SUBDATASETS': False,
            'DATA_TYPE': 0,
            'GPKG_LAYER_NAME': '',
            'OUTPUT': output_path
        }
        
        result = self.run_processing_algorithm("native:reprojectlayer", params_reproject)
        if result and os.path.exists(output_path):
            reprojected_layer = QgsVectorLayer(output_path, f"{vector_layer.name()}_reprojected", "ogr")
            if reprojected_layer.isValid():
                self.log_message(f"Vector layer reprojected successfully to {target_crs.authid()}")
                return reprojected_layer
            else:
                self.log_message("ERROR: Reprojected vector layer is not valid.")
                return None
        else:
            self.log_message("ERROR: Vector layer reprojection failed.")
            return None

    def process_layer_with_dem_clipping(self, layer, layer_name, output_filename, display=True):
        """
        Process a layer (vector or raster) by rasterizing if needed and clipping to DEM extent
        
        Args:
            layer: QgsMapLayer (vector or raster)
            layer_name: Name for logging
            output_filename: Output filename in geometry folder
            display: bool, whether to load the processed layer to QGIS (default: True)
            
        Returns:
            str: Path to processed layer or None if failed
        """
        if not layer:
            QMessageBox.warning(self.dialog, "Input Error", f"Please select a {layer_name} layer.")
            return None
        
        if not self.check_prerequisites():
            return None
        
        # Get input CRS and DEM info
        input_crs = self.dialog.get_crs()
        
        # Use reprojected DEM layer if available, otherwise get from combo box
        if not self.dem_layer:
            dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
            # Check and reproject if needed
            self.check_and_reproject_dem(dem_layer)
        
        # Use the reprojected DEM layer stored in self.dem_layer
        dem_layer = self.dem_layer if self.dem_layer else self.dialog.mMapLayerComboBox_dem.currentLayer()
        
        # Get DEM extent and resolution from reprojected DEM
        extent_str, pixel_size_x, pixel_size_y = self.get_dem_extent_and_resolution()
        if not extent_str:
            return None
        
        # Get reprojected DEM extent with CRS for clipping
        dem_extent = self.dem_layer.extent()
        dem_crs = self.dem_layer.crs()
        # Format: xmin,xmax,ymin,ymax [EPSG:xxxxx]
        extent_projwin = f"{dem_extent.xMinimum()},{dem_extent.xMaximum()},{dem_extent.yMinimum()},{dem_extent.yMaximum()} [{dem_crs.authid()}]"
        
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        output_path = os.path.join(geometry_folder, output_filename)
        
        # Check if already processed
        if os.path.exists(output_path):
            self.log_message(f"{layer_name} already processed.")
            if display:
                self.load_layer(output_path, layer_name)
            return output_path
        
        self.log_message(f"Processing {layer_name}...")
        
        # Determine if layer is vector or raster
        if isinstance(layer, QgsVectorLayer):
            self.log_message(f"Vector layer detected. Checking CRS for {layer_name}...")
            
            # Check if vector layer CRS matches input CRS
            layer_crs = layer.crs()
            if input_crs.isValid() and layer_crs.isValid():
                if layer_crs.authid() != input_crs.authid():
                    self.log_message(f"{layer_name} CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                    temp_reprojected_path = os.path.join(geometry_folder, f"temp_{layer_name}_reprojected.gpkg")
                    reprojected_layer = self.reproject_vector_layer(layer, input_crs, temp_reprojected_path)
                    if reprojected_layer:
                        layer = reprojected_layer
                    else:
                        self.log_message(f"WARNING: Reprojection failed. Using original {layer_name} layer.")
                else:
                    self.log_message(f"{layer_name} CRS matches input CRS. No reprojection needed.")
            
            # Rasterize vector layer using DEM's georeferenced units
            self.log_message(f"Rasterizing {layer_name}...")
            
            # Get the first field for rasterization (or use a default)
            fields = layer.fields()
            field_name = fields[0].name() if fields else 'id'
            
            # Calculate width and height from reprojected DEM pixel size
            dem_extent = self.dem_layer.extent()
            width = int((dem_extent.xMaximum() - dem_extent.xMinimum()) / pixel_size_x)
            height = int((dem_extent.yMaximum() - dem_extent.yMinimum()) / pixel_size_y)
            
            params_rasterize = {
                'INPUT': layer,
                'FIELD': field_name,
                'BURN': 0,
                'USE_Z': False,
                'UNITS': 1,  # Georeferenced units
                'WIDTH': width,
                'HEIGHT': height,
                'EXTENT': extent_str,
                'NODATA': 0,
                'OPTIONS': None,
                'DATA_TYPE': 5,  # Float32
                'INIT': None,
                'INVERT': False,
                'EXTRA': '',
                'OUTPUT': output_path
            }
            
            result = self.run_processing_algorithm("gdal:rasterize", params_rasterize)
            
        elif isinstance(layer, QgsRasterLayer):
            self.log_message(f"Raster layer detected. Clipping {layer_name} to DEM extent...")
            
            # Check if raster layer CRS matches input CRS
            layer_crs = layer.crs()
            if input_crs.isValid() and layer_crs.isValid():
                if layer_crs.authid() != input_crs.authid():
                    self.log_message(f"{layer_name} CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                    temp_reprojected_path = os.path.join(geometry_folder, f"temp_{layer_name}_reprojected.tif")
                    
                    params_warp = {
                        'INPUT': layer.source(),
                        'SOURCE_CRS': None,
                        'TARGET_CRS': input_crs,
                        'RESAMPLING': 0,  # Nearest neighbor
                        'NODATA': None,
                        'TARGET_RESOLUTION': None,
                        'OPTIONS': None,
                        'DATA_TYPE': 0,
                        'TARGET_EXTENT': None,
                        'TARGET_EXTENT_CRS': None,
                        'MULTITHREADING': False,
                        'EXTRA': '',
                        'OUTPUT': temp_reprojected_path
                    }
                    
                    warp_result = self.run_processing_algorithm("gdal:warpreproject", params_warp)
                    if warp_result and os.path.exists(temp_reprojected_path):
                        layer = QgsRasterLayer(temp_reprojected_path, f"{layer_name}_reprojected")
                        self.log_message(f"{layer_name} reprojected successfully.")
                    else:
                        self.log_message(f"WARNING: Reprojection failed. Using original {layer_name} layer.")
                else:
                    self.log_message(f"{layer_name} CRS matches input CRS. No reprojection needed.")
            
            # Clip raster to DEM extent using reprojected DEM extent
            params_clip = {
                'INPUT': layer,
                'PROJWIN': extent_projwin,  # Extent from reprojected DEM with CRS
                'OVERCRS': False,
                'NODATA': None,
                'OPTIONS': None,
                'DATA_TYPE': 0,
                'EXTRA': '',
                'OUTPUT': output_path
            }
            
            result = self.run_processing_algorithm("gdal:cliprasterbyextent", params_clip)
        
        else:
            self.log_message(f"ERROR: Unsupported layer type for {layer_name}")
            return None
        
        if result:
            if display:
                self.load_layer(output_path, layer_name)
            self.log_message(f"{layer_name} processing completed successfully.")
            return output_path
        else:
            self.log_message(f"{layer_name} processing failed.")
            return None

    def load_land_cover_lookup_table(self):
        """
        Load land cover lookup table from file or use default mapping.
        
        Returns:
            dict: Mapping of grid values to type_int (1=Forest, 2=Impervious, 3=Pervious)
        """
        lookup_file = self.dialog.lineEdit_land_cover_lookup.text()
        
        # Default lookup table
        default_lookup = {
            1: 3,   # Waterbody -> Pervious
            2: 2,   # Glacier -> Impervious
            3: 3,   # Snow -> Pervious
            4: 1,   # Forest -> Forest
            5: 3,   # Riverbed -> Pervious
            6: 2,   # Built-up area -> Impervious
            7: 3,   # Cropland -> Pervious
            8: 3,   # Bare soil -> Pervious
            9: 2,   # Bare rock -> Impervious
            10: 3,  # Grassland -> Pervious
            11: 3   # Other wooded land -> Pervious
        }
        
        if not lookup_file or not os.path.exists(lookup_file):
            self.log_message("No lookup table provided or file not found. Using default mapping.")
            return default_lookup
        
        # Try to read lookup table from file (CSV format expected)
        try:
            import csv
            lookup_mapping = {}
            with open(lookup_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        grid_value = int(row.get('Grid value', 0))
                        type_str = row.get('Type', '').strip()
                        
                        # Map type string to type_int
                        if type_str == 'Forest':
                            type_int = 1
                        elif type_str == 'Impervious':
                            type_int = 2
                        elif type_str == 'Pervious':
                            type_int = 3
                        else:
                            self.log_message(f"WARNING: Unknown type '{type_str}' in lookup table. Skipping.")
                            continue
                        
                        lookup_mapping[grid_value] = type_int
                    except (ValueError, KeyError) as e:
                        self.log_message(f"WARNING: Error parsing lookup table row: {e}")
                        continue
            
            if lookup_mapping:
                self.log_message(f"Loaded {len(lookup_mapping)} entries from lookup table.")
                return lookup_mapping
            else:
                self.log_message("Lookup table is empty. Using default mapping.")
                return default_lookup
        except Exception as e:
            self.log_message(f"ERROR reading lookup table: {e}. Using default mapping.")
            return default_lookup

    def reclassify_land_use_raster(self, input_raster_path, output_path, lookup_mapping):
        """
        Reclassify land use raster values based on lookup table.
        
        Args:
            input_raster_path: Path to clipped land use raster
            output_path: Path for reclassified output
            lookup_mapping: Dictionary mapping grid values to type_int (1, 2, or 3)
            
        Returns:
            bool: True if successful, False otherwise
        """
        self.log_message("Reclassifying land use raster based on lookup table...")
        
        # Build reclassification table
        # Format: [min1, max1, value1, min2, max2, value2, ...]
        reclass_table = []
        for grid_value, type_int in sorted(lookup_mapping.items()):
            reclass_table.extend([str(grid_value-0.5), str(grid_value+0.5), str(type_int)])
        
        self.log_message(f"Reclassification table: {reclass_table}")
        
        params_reclass = {
            'INPUT_RASTER': input_raster_path,
            'RASTER_BAND': 1,
            'TABLE': reclass_table,
            'NO_DATA': -9999,
            'RANGE_BOUNDARIES': 0,  # min <= value <= max
            'NODATA_FOR_MISSING': False,
            'DATA_TYPE': 3,  # Int16 (since values are 1, 2, or 3)
            'CREATE_OPTIONS': None,
            'OUTPUT': output_path
        }
        
        result = self.run_processing_algorithm("native:reclassifybytable", params_reclass)
        
        if result and os.path.exists(output_path):
            self.log_message("Land use reclassification completed successfully.")
            return True
        else:
            self.log_message("ERROR: Land use reclassification failed.")
            return False

    def process_land_use(self):
        """Process Land Use layer with clipping and reclassification"""
        layer = self.dialog.mMapLayerComboBox_land_cover.currentLayer()
        if not layer:
            QMessageBox.warning(self.dialog, "Input Error", "Please select a Land Cover layer.")
            return
        
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        clipped_path = os.path.join(geometry_folder, "3_land_use_clipped.tif")
        final_path = os.path.join(geometry_folder, "3_land_use.tif")
        
        # Check if final land use layer already exists
        if os.path.exists(final_path):
            self.log_message("Final land use layer already exists. Loading existing file...")
            self.land_use_layer = final_path
            self.load_layer(final_path, "3_Land_Use")
            return
        
        # Step 1: Clip land use layer to DEM extent (display=False for intermediate layer)
        self.log_message("Clipping land use layer to DEM extent...")
        result = self.process_layer_with_dem_clipping(
            layer, "3_Land_Use_Clip", "3_land_use_clipped.tif", display=False)
        
        if not result:
            self.log_message("ERROR: Land use clipping failed.")
            return
        
        # Step 2: Load lookup table
        lookup_mapping = self.load_land_cover_lookup_table()
        
        # Step 3: Reclassify clipped raster
        success = self.reclassify_land_use_raster(result, final_path, lookup_mapping)
        
        if success:
            self.land_use_layer = final_path
            self.load_layer(final_path, "3_Land_Use")
            self.log_message("Land Use layer processed, clipped, and reclassified successfully.")
            
            # Clean up intermediate clipped file
            try:
                if os.path.exists(clipped_path):
                    os.remove(clipped_path)
            except:
                pass
        else:
            self.log_message("ERROR: Land use reclassification failed.")

    def process_soil(self):
        """Process Soil layer"""
        layer = self.dialog.mMapLayerComboBox_soil.currentLayer()
        result = self.process_layer_with_dem_clipping(
            layer, "3_Soil", "3_soil.tif")
        if result:
            self.log_message("Soil layer processed and clipped to DEM extent.")

    def process_geology(self):
        """Process Geology layer"""
        layer = self.dialog.mMapLayerComboBox_geology.currentLayer()
        result = self.process_layer_with_dem_clipping(
            layer, "Geology", "geology_processed.tif")
        if result:
            self.log_message("Geology layer processed and clipped to DEM extent.")

