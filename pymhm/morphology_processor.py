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
    QgsProject
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
        
        # Paths for geometry processing outputs
        self.filled_dem_path = None
        self.flow_acc_path = None
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
            'flow_acc_path': "2_flow_accumulation.sdat",
            'flow_dir_path': "2_flow_direction.sdat",
            'channel_network_vector_path': "2_channel_network.shp",
            'snapped_points_path': "3_pour_points_snapped.gpkg",
            'watershed_raster_path': "4_watershed_raster.sdat",
            'watershed_vector_path': "4_watershed_final.gpkg",

            # New hydrological processing files
            'aspect_path': "1_dem_aspect.sdat",
            'slope_path': "1_dem_slope.sdat",
            'flow_direction_path': "2_flow_direction.sdat",
            'flow_accumulation_path': "2_flow_accumulation.sdat",
            'gauge_position_path': "3_gauge_position.sdat"
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

        if found_any:
            self.log_message(
                "Project state loaded. You may skip completed steps.")
        else:
            self.log_message(
                "No existing project files found in 'Geometry' folder.")

    # --- Geometry Processing Methods ---

    def fill_dem(self):
        """Step 1: Fill Sinks in DEM using SAGA NG's Wang & Liu algorithm."""
        self.log_message("\n--- Starting Geometry Step 1: Fill DEM ---")
        if not self.check_prerequisites():
            return

        dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        self.log_message(f"Input DEM: {dem_layer.name()}")

        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.filled_dem_path = os.path.join(
            geometry_folder, "1_dem_filled.sdat")
        params = {
            'ELEV': dem_layer,
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

    def create_network(self):
        """Step 2: Create Flow Accumulation, Flow Direction, and Channel Network."""
        self.log_message("\n--- Starting Geometry Step 2: Create Network ---")
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 1 (Fill DEM) successfully first.")
            return

        # Sub-step 2a: Create Flow Accumulation and Flow Direction
        self.log_message(
            "Sub-step 2a: Creating Flow Accumulation and Flow Direction...")
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.flow_acc_path = os.path.join(
            geometry_folder, "2_flow_accumulation.sdat")
        params_acc = {
            'DEM': self.filled_dem_path,
            'PREPROCESSING': 1,  # Fill Sinks
            'FLOW_ROUTING': 4,  # D8
            'TCA': self.flow_acc_path,
            'SCA': 'TEMPORARY_OUTPUT',
            'FLOW_PATH_LENGTH': 'TEMPORARY_OUTPUT'
        }
        result_acc = self.run_processing_algorithm(
            "sagang:flowaccumulationonestep", params_acc)
        if not result_acc:
            self.flow_acc_path = None
            self.flow_dir_path = None
            return
        self.load_layer(self.flow_acc_path, "2_Flow_Accumulation")

        # Sub-step 2b: Create vector Channel Network
        self.log_message("Sub-step 2b: Creating vector Channel Network...")
        self.flow_dir_path = os.path.join(
            geometry_folder, "2_flow_direction.sdat")
        self.channel_network_vector_path = os.path.join(
            geometry_folder, "2_channel_network.shp")
        params_chnl = {
            'ELEVATION': self.filled_dem_path,
            'SINKROUTE': None,
            'CHNLNTWRK': 'TEMPORARY_OUTPUT',
            'CHNLROUTE': self.flow_dir_path,
            'SHAPES': self.channel_network_vector_path,
            'INIT_GRID': self.flow_acc_path,
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
        if not result_chnl:
            self.channel_network_vector_path = None
            return
        self.load_layer(self.channel_network_vector_path,
                        "2_Channel_Network", is_raster=False)

    def snap_points(self):
        """Step 3: Snap Pour Points to the nearest high-order channel segment."""
        self.log_message(
            "\n--- Starting Geometry Step 3: Snap Pour Points ---")
        if not self.check_prerequisites(needs_pour_points=True):
            return
        if not self.channel_network_vector_path or not os.path.exists(self.channel_network_vector_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 2 (Create Network) successfully first.")
            return

        pour_points_layer = self.dialog.mMapLayerComboBox_pour_points.currentLayer()
        channel_network_layer = QgsVectorLayer(
            self.channel_network_vector_path, "Channel Network", "ogr")
        if not channel_network_layer.isValid():
            self.log_message(
                f"ERROR: Could not load channel network layer from {self.channel_network_vector_path}")
            return

        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.snapped_points_path = os.path.join(
            geometry_folder, "3_pour_points_snapped.gpkg")

        # Call the new custom snapping function
        result_path = self.snap_points_to_network(
            pour_points_layer=pour_points_layer,
            channel_network_layer=channel_network_layer,
            output_path=self.snapped_points_path,
            order_field_name='Order',  # SAGA Channel Network tool typically uses 'ORDER'
            high_order_distance=1000.0,
            max_snap_distance=5000.0
        )

        if result_path:
            self.load_layer(result_path,
                            "3_Pour_Points_Snapped", is_raster=False)
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
            geometry_folder, "4_watershed_final.gpkg")

        # Create a combined watershed vector layer
        self.create_combined_watershed_layer(
            watershed_outputs, self.watershed_vector_path)

        if os.path.exists(self.watershed_vector_path):
            self.load_layer(self.watershed_vector_path,
                            "4_Watershed_Final", is_raster=False)
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
                                     QgsWkbTypes.Point, pour_points_layer.crs(), "GPKG")
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
        # Get the input DEM from the combo box
        input_dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        if not input_dem_layer:
            QMessageBox.warning(self.dialog, "Input Error",
                                "Please select an input DEM layer.")
            return

        input_dem_path = input_dem_layer.source()
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
        # Get the input DEM from the combo box
        input_dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        if not input_dem_layer:
            QMessageBox.warning(self.dialog, "Input Error",
                                "Please select an input DEM layer.")
            return

        input_dem_path = input_dem_layer.source()
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
        """Process Flow Direction using D8 method"""
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 1 (Fill DEM) successfully first.")
            return

        # Check if flow direction D8 already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.flow_direction_path = os.path.join(
            geometry_folder, "2_flow_direction.sdat")
        if self.flow_direction_path and os.path.exists(self.flow_direction_path):
            self.log_message(
                "Flow Direction (D8) already exists. Loading existing file...")
            self.load_layer(self.flow_direction_path, "2_Flow_Direction")
            return

        self.log_message("Processing Flow Direction (D8)...")

        params_flow_dir = {
            'ELEVATION': self.filled_dem_path,
            'FLOW': self.flow_direction_path,
            'METHOD': 0  # D8
        }

        result = self.run_processing_algorithm(
            "sagang:flowdirection", params_flow_dir)
        if result:
            self.load_layer(self.flow_direction_path, "2_Flow_Direction")
            self.log_message(
                "Flow Direction (D8) processing completed successfully.")
        else:
            self.log_message("Flow Direction (D8) processing failed.")

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

        params_flow_acc = {
            'ELEVATION': self.filled_dem_path,
            'FLOW': self.flow_accumulation_path,
            'METHOD': 0,  # D8
            'CONVERGE': 1.1
        }

        result = self.run_processing_algorithm(
            "sagang:flowaccumulation", params_flow_acc)
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
            geometry_folder, "3_gauge_position.sdat")
        if self.gauge_position_path and os.path.exists(self.gauge_position_path):
            self.log_message(
                "Gauge Position already exists. Loading existing file...")
            self.load_layer(self.gauge_position_path, "3_Gauge_Position")
            return

        # This method should be called with a pour points layer path
        # For now, we'll use a placeholder path that should be provided by the user
        pour_points_path = QFileDialog.getOpenFileName(
            self.dialog, "Select Pour Points Layer", "", "GeoPackage (*.gpkg);;Shapefile (*.shp)")[0]

        if not pour_points_path:
            return

        self.log_message("Processing Gauge Position...")

        params_gauge = {
            'INPUT': pour_points_path,
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
            self.load_layer(self.gauge_position_path, "3_Gauge_Position")
            self.log_message(
                "Gauge Position processing completed successfully.")
        else:
            self.log_message("Gauge Position processing failed.")

    def process_layer_with_dem_clipping(self, layer, layer_name, output_filename):
        """
        Process a layer (vector or raster) by rasterizing if needed and clipping to DEM extent
        
        Args:
            layer: QgsMapLayer (vector or raster)
            layer_name: Name for logging
            output_filename: Output filename in geometry folder
            
        Returns:
            str: Path to processed layer or None if failed
        """
        if not layer:
            QMessageBox.warning(self.dialog, "Input Error", f"Please select a {layer_name} layer.")
            return None
        
        if not self.check_prerequisites():
            return None
        
        # Get DEM extent and resolution
        extent_str, pixel_size_x, pixel_size_y = self.get_dem_extent_and_resolution()
        if not extent_str:
            return None
        
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        output_path = os.path.join(geometry_folder, output_filename)
        
        # Check if already processed
        if os.path.exists(output_path):
            self.log_message(f"{layer_name} already processed. Loading existing file...")
            self.load_layer(output_path, f"Processed_{layer_name}")
            return output_path
        
        self.log_message(f"Processing {layer_name}...")
        
        # Determine if layer is vector or raster
        if isinstance(layer, QgsVectorLayer):
            self.log_message(f"Vector layer detected. Rasterizing {layer_name}...")
            
            # Get the first field for rasterization (or use a default)
            fields = layer.fields()
            field_name = fields[0].name() if fields else 'id'
            
            params_rasterize = {
                'INPUT': layer,
                'FIELD': field_name,
                'BURN': 0,
                'USE_Z': False,
                'UNITS': 1,  # Georeferenced units
                'WIDTH': 0,  # Will be calculated from pixel size
                'HEIGHT': 0,  # Will be calculated from pixel size
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
            
            # Use DEM layer directly as mask
            dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
            
            params_clip = {
                'INPUT': layer,
                'MASK': dem_layer,
                'SOURCE_CRS': None,
                'TARGET_CRS': None,
                'TARGET_EXTENT': extent_str,
                'NODATA': None,
                'ALPHA_BAND': False,
                'CROP_TO_CUTLINE': True,
                'KEEP_RESOLUTION': False,
                'SET_RESOLUTION': False,
                'X_RESOLUTION': None,
                'Y_RESOLUTION': None,
                'MULTITHREADING': False,
                'OPTIONS': None,
                'DATA_TYPE': 0,
                'EXTRA': '',
                'OUTPUT': output_path
            }
            
            result = self.run_processing_algorithm("gdal:cliprasterbymasklayer", params_clip)
        
        else:
            self.log_message(f"ERROR: Unsupported layer type for {layer_name}")
            return None
        
        if result:
            self.load_layer(output_path, f"Processed_{layer_name}")
            self.log_message(f"{layer_name} processing completed successfully.")
            return output_path
        else:
            self.log_message(f"{layer_name} processing failed.")
            return None

    def process_land_use(self):
        """Process Land Use layer"""
        layer = self.dialog.mMapLayerComboBox_land_cover.currentLayer()
        result = self.process_layer_with_dem_clipping(
            layer, "Land Use", "land_use_processed.tif")
        if result:
            self.log_message("Land Use layer processed and clipped to DEM extent.")

    def process_soil(self):
        """Process Soil layer"""
        layer = self.dialog.mMapLayerComboBox_soil.currentLayer()
        result = self.process_layer_with_dem_clipping(
            layer, "Soil", "soil_processed.tif")
        if result:
            self.log_message("Soil layer processed and clipped to DEM extent.")

    def process_geology(self):
        """Process Geology layer"""
        layer = self.dialog.mMapLayerComboBox_geology.currentLayer()
        result = self.process_layer_with_dem_clipping(
            layer, "Geology", "geology_processed.tif")
        if result:
            self.log_message("Geology layer processed and clipped to DEM extent.")

