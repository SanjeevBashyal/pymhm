# -*- coding: utf-8 -*-
"""
Morphology/Geometry processing module for pymhm
Contains all geometry-related processing methods
"""
import os
import math
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
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsRectangle,
    QgsPointXY
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
        # Note: load_layer is overridden below to support skip_loading flag
        self.get_dem_extent_and_resolution = dialog.get_dem_extent_and_resolution

        # DEM layer reference (reprojected if needed)
        self.dem_layer = None  # Will store the reprojected DEM layer if CRS differs

        # Processed layer references
        # Final land use layer with reclassified values (1, 2, or 3)
        self.land_use_layer = None

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
        self.flow_accumulation_area_path = None
        self.gauge_position_path = None
        
        # Layer processing paths
        self.geology_path = None
        
        # Flag to skip loading layers (used in execute_all_processing)
        self.skip_loading = False
        
        # Lat/Lon header information (L0, L1, L11, L2)
        self.L0 = None  # Dictionary with ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value
        self.L1 = None
        self.L11 = None
        self.L2 = None

    def load_layer(self, path, name, is_raster=True):
        """
        Override load_layer to check skip_loading flag.
        If skip_loading is True, don't load the layer.
        """
        if self.skip_loading:
            self.log_message(f"Skipping loading layer: {name} (execute_all mode)")
            return
        # Call the parent load_layer method
        return self.dialog.load_layer(path, name, is_raster)

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
            'flow_direction_path': "2_flow_direction.tif",
            'flow_accumulation_path': "2_flow_accumulation.tif",
            'flow_accumulation_area_path': "2_flow_accumulation_area.sdat",
            'gauge_position_path': "2_gauge_position.tif"
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
        reprojected_dem_path = os.path.join(
            geometry_folder, "0_dem_reprojected.tif")
        if os.path.exists(reprojected_dem_path):
            reprojected_layer = QgsRasterLayer(
                reprojected_dem_path, "1_DEM_Reprojected")
            if reprojected_layer.isValid():
                self.dem_layer = reprojected_layer
                self.log_message(
                    f"Found existing reprojected DEM: 0_dem_reprojected.tif")
                found_any = True
            else:
                self.log_message(
                    f"WARNING: Reprojected DEM file exists but is not valid.")

        # Check for filled DEM and load it to QGIS interface
        if self.filled_dem_path and os.path.exists(self.filled_dem_path):
            self.log_message(
                f"Loading existing filled DEM to QGIS interface...")
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
            self.log_message(
                "WARNING: Input CRS is not valid. Using DEM CRS as-is.")
            self.dem_layer = dem_layer
            return dem_layer

        dem_crs = dem_layer.crs()
        if not dem_crs.isValid():
            self.log_message("WARNING: DEM CRS is not valid. Using DEM as-is.")
            self.dem_layer = dem_layer
            return dem_layer

        # Check if CRS matches
        if dem_crs.authid() == input_crs.authid():
            self.log_message(
                f"DEM CRS ({dem_crs.authid()}) matches input CRS. No reprojection needed.")
            self.dem_layer = dem_layer
            return dem_layer

        # Need to reproject - check if reprojected file already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        reprojected_dem_path = os.path.join(
            geometry_folder, "0_dem_reprojected.tif")

        if os.path.exists(reprojected_dem_path):
            self.log_message(f"Found existing reprojected DEM. Loading it...")
            reprojected_layer = QgsRasterLayer(
                reprojected_dem_path, "1_DEM_Reprojected")
            if reprojected_layer.isValid():
                self.dem_layer = reprojected_layer
                return reprojected_layer
            else:
                self.log_message(
                    "WARNING: Existing reprojected DEM is not valid. Reprojecting again...")

        # Reproject DEM
        self.log_message(
            f"DEM CRS ({dem_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")

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

        result = self.run_processing_algorithm(
            "gdal:warpreproject", params_warp)
        if result and os.path.exists(reprojected_dem_path):
            reprojected_layer = QgsRasterLayer(
                reprojected_dem_path, "1_DEM_Reprojected")
            if reprojected_layer.isValid():
                self.log_message(
                    f"DEM reprojected successfully to {input_crs.authid()}")
                # Do not add reprojected layer to QGIS project - it's only used for processing
                self.dem_layer = reprojected_layer
                return reprojected_layer
            else:
                self.log_message(
                    "ERROR: Reprojected DEM layer is not valid. Using original DEM.")
                self.dem_layer = dem_layer
                return dem_layer
        else:
            self.log_message(
                "ERROR: DEM reprojection failed. Using original DEM.")
            self.dem_layer = dem_layer
            return dem_layer

    def convert_dem_to_asc(self):
        """Convert DEM from TIF to ASC format. First reprojects to input CRS, then converts to ASC."""
        self.log_message("\n--- Converting DEM to ASC format ---")
        if not self.check_prerequisites():
            return

        # Output folder - save to input/morph folder
        input_folder = os.path.join(self.dialog.project_folder, "input")
        morph_folder = os.path.join(input_folder, "morph")
        os.makedirs(morph_folder, exist_ok=True)
        dem_asc_path = os.path.join(morph_folder, "dem.asc")
        
        # Geometry folder for temporary files
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")

        # Check if ASC file already exists
        if os.path.exists(dem_asc_path):
            self.log_message(
                "DEM ASC file already exists. Loading existing file...")
            self.load_layer(dem_asc_path, "DEM_ASC", is_raster=True)
            return

        # Get DEM layer from combo box
        dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        if not dem_layer:
            QMessageBox.critical(
                self.dialog, "Missing Input",
                "Please select a DEM Raster Layer.")
            return

        # Get input CRS from CRS widget
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.critical(
                self.dialog, "Missing Input",
                "Please select a valid Coordinate Reference System.")
            return

        self.log_message(f"Input DEM: {dem_layer.name()}")

        # Get DEM CRS
        dem_crs = dem_layer.crs()
        if not dem_crs.isValid():
            self.log_message(
                "WARNING: DEM CRS is not valid. Proceeding with conversion...")

        self.log_message(
            f"DEM CRS: {dem_crs.authid() if dem_crs.isValid() else 'Unknown'}")
        self.log_message(f"Target CRS: {input_crs.authid()}")

        # Get DEM source path
        dem_source = dem_layer.source()
        if not dem_source or not os.path.exists(dem_source):
            self.log_message(f"ERROR: DEM source file not found: {dem_source}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"DEM source file not found: {dem_source}")
            return

        # Step 1: Reproject DEM to input CRS (if needed)
        reprojected_dem_path = os.path.join(
            geometry_folder, "dem_reprojected.tif")
        dem_source_for_conversion = dem_source

        # Check if reprojection is needed
        needs_reprojection = False
        if dem_crs.isValid() and input_crs.isValid():
            if dem_crs.authid() != input_crs.authid():
                needs_reprojection = True

        if needs_reprojection:
            self.log_message(
                f"Reprojecting DEM from {dem_crs.authid()} to {input_crs.authid()}...")

            # Check if reprojected file already exists
            if os.path.exists(reprojected_dem_path):
                self.log_message("Found existing reprojected DEM. Using it...")
                dem_source_for_conversion = reprojected_dem_path
            else:
                # Reproject DEM
                params_warp = {
                    'INPUT': dem_source,
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

                result_warp = self.run_processing_algorithm(
                    "gdal:warpreproject", params_warp)
                if result_warp and os.path.exists(reprojected_dem_path):
                    self.log_message("DEM reprojected successfully!")
                    dem_source_for_conversion = reprojected_dem_path
                else:
                    self.log_message(
                        "ERROR: DEM reprojection failed. Using original DEM.")
                    QMessageBox.warning(
                        self.dialog, "Warning",
                        "DEM reprojection failed. Converting original DEM to ASC.")
        else:
            self.log_message(
                "DEM CRS matches target CRS. No reprojection needed.")

        # Step 2: Convert to ASC format
        self.log_message(f"Converting DEM to ASC format...")
        self.log_message(f"Output file: {dem_asc_path}")

        # Parameters for GDAL translate
        params_translate = {
            'INPUT': dem_source_for_conversion,
            'TARGET_CRS': None,  # Already in target CRS from reprojection
            'NODATA': -9999,
            'COPY_SUBDATASETS': False,
            'OPTIONS': None,
            'EXTRA': '',
            'DATA_TYPE': 0,  # Use input data type
            'OUTPUT': dem_asc_path
        }

        result = self.run_processing_algorithm(
            "gdal:translate", params_translate)
        if result and os.path.exists(dem_asc_path):
            self.load_layer(dem_asc_path, "DEM_ASC", is_raster=True)
            self.log_message("DEM converted to ASC format successfully!")
        else:
            self.log_message("ERROR: Failed to convert DEM to ASC format.")

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
            self.log_message(
                "Filled DEM already exists. Loading existing file...")
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

        # Check if flow accumulation area exists (required dependency for channel network)
        if not self.flow_accumulation_area_path or not os.path.exists(self.flow_accumulation_area_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Flow accumulation area processing failed or file not found.")
            return

        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.channel_network_vector_path = os.path.join(
            geometry_folder, "2_channel_network.shp")

        # Check if Channel Network already exists, otherwise process it
        if self.channel_network_vector_path and os.path.exists(self.channel_network_vector_path):
            self.log_message(
                "Channel Network already exists. Loading existing file...")
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
                'INIT_GRID': self.flow_accumulation_area_path,
                'INIT_METHOD': 2,       # Greater than
                # Threshold for channel initiation (area-based, e.g., m²)
                'INIT_VALUE': 10000000.0,
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
                self.log_message(
                    "Channel Network processing completed successfully.")
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
            self.log_message(
                "Snapped points already exist. Loading existing file...")
            self.load_layer(self.snapped_points_path,
                            "2_pour_points_snapped", is_raster=False)
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
            order_field_name='Order',  # SAGA Channel Network tool typically uses 'ORDER'
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
        if not self.flow_direction_path or not os.path.exists(self.flow_direction_path):
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
            QMessageBox.warning(self.dialog, "Error",
                                "No snapped points found.")
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

            if result_ws and os.path.exists(watershed_raster_path):
                # Polygonize the watershed raster and save as shapefile
                watershed_vector_path = os.path.join(
                    geometry_folder, f"4_watershed_{clean_name}.shp")

                self.log_message(
                    f"Polygonizing watershed raster for: {name_attr}")

                params_poly = {
                    'INPUT': watershed_raster_path,
                    'BAND': 1,
                    'FIELD': 'DN',
                    'EIGHT_CONNECTEDNESS': False,
                    'EXTRA': '',
                    'OUTPUT': watershed_vector_path
                }

                result_poly = self.run_processing_algorithm(
                    "gdal:polygonize", params_poly)

                if result_poly and os.path.exists(watershed_vector_path):
                    self.log_message(
                        f"Watershed polygon saved: {watershed_vector_path}")

                    watershed_outputs.append({
                        'name': name_attr,
                        'clean_name': clean_name,
                        'raster_path': watershed_raster_path,
                        'vector_path': watershed_vector_path,
                        'point': point
                    })
                else:
                    self.log_message(
                        f"Warning: Failed to polygonize watershed for: {name_attr}")
                    # Still add to outputs even if polygonization failed
                    watershed_outputs.append({
                        'name': name_attr,
                        'clean_name': clean_name,
                        'raster_path': watershed_raster_path,
                        'vector_path': None,
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

        # Merge all polygonized watershed shapefiles into one
        self.log_message("Merging all watershed vector layers...")
        self.merged_watershed_path = os.path.join(
            geometry_folder, "4_watershed_merged_vector.shp")

        # Collect all vector layer paths that exist
        vector_layer_paths = []
        for watershed_info in watershed_outputs:
            if 'vector_path' in watershed_info and watershed_info['vector_path']:
                if os.path.exists(watershed_info['vector_path']):
                    vector_layer_paths.append(watershed_info['vector_path'])

        if vector_layer_paths:
            self.log_message(
                f"Merging {len(vector_layer_paths)} watershed vector layers...")

            params_merge = {
                'LAYERS': vector_layer_paths,
                'CRS': None,
                'OUTPUT': self.merged_watershed_path
            }

            result_merge = self.run_processing_algorithm(
                "native:mergevectorlayers", params_merge)

            if result_merge and os.path.exists(self.merged_watershed_path):
                self.log_message(
                    f"Merged watershed vector saved: {self.merged_watershed_path}")
            else:
                self.log_message(
                    "Warning: Failed to merge watershed vector layers.")
        else:
            self.log_message("Warning: No valid vector layers found to merge.")

        if os.path.exists(self.merged_watershed_path):
            self.load_layer(self.merged_watershed_path,
                            "4_watershed_merged", is_raster=False)
        else:
            self.merged_watershed_path = None

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
            geometry_folder, "2_flow_direction.tif")
        if self.flow_direction_path and os.path.exists(self.flow_direction_path):
            self.log_message(
                "Flow Direction already exists. Loading existing file...")
            self.load_layer(self.flow_direction_path, "2_Flow_Direction")
            return

        self.log_message(
            "Processing Flow Direction using Channel Network and Drainage Basins...")

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
        """Process Flow Accumulation using D8 method. Converts area to pixels and saves as integer raster."""
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error", "Please run Step 1 (Fill DEM) successfully first.")
            return

        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        flow_accumulation_area_path = os.path.join(
            geometry_folder, "2_flow_accumulation_area.sdat")
        flow_accumulation_pixels_float_path = os.path.join(
            geometry_folder, "2_flow_accumulation_pixels_float.tif")
        self.flow_accumulation_path = os.path.join(
            geometry_folder, "2_flow_accumulation.tif")

        # Check if final flow accumulation already exists
        if self.flow_accumulation_path and os.path.exists(self.flow_accumulation_path):
            self.log_message(
                "Flow Accumulation (pixels, integer) already exists. Loading existing file...")
            self.load_layer(self.flow_accumulation_path,
                            "2_Flow_Accumulation")
            return

        self.log_message("Processing Flow Accumulation (Area)...")

        # Step 1: Calculate flow accumulation area
        params_acc = {
            'DEM': self.filled_dem_path,
            'PREPROCESSING': 0,  # No preprocessing (already filled)
            'FLOW_ROUTING': 0,  # D8
            'TCA': flow_accumulation_area_path,
            'SCA': 'TEMPORARY_OUTPUT',
            'FLOW_PATH_LENGTH': 'TEMPORARY_OUTPUT'
        }
        result_acc = self.run_processing_algorithm(
            "sagang:flowaccumulationonestep", params_acc)

        if not result_acc or not os.path.exists(flow_accumulation_area_path):
            self.log_message(
                "ERROR: Flow Accumulation (Area) processing failed.")
            return

        self.log_message("Flow Accumulation (Area) calculated successfully.")
        
        # Store flow accumulation area path in self attribute
        self.flow_accumulation_area_path = flow_accumulation_area_path

        # Step 2: Get cell size from filled DEM
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM to get cell size.")
            return

        # Get pixel size (cell size)
        raster_extent = filled_dem_layer.extent()
        width = filled_dem_layer.width()
        height = filled_dem_layer.height()
        cell_size_x = (raster_extent.xMaximum() -
                       raster_extent.xMinimum()) / width
        cell_size_y = (raster_extent.yMaximum() -
                       raster_extent.yMinimum()) / height

        # Use average cell size if they differ
        cell_size = (cell_size_x + cell_size_y) / 2.0
        self.log_message(f"Cell size: {cell_size:.2f} meters")

        # Step 3: Convert area to pixels using raster calculator
        self.log_message("Converting flow accumulation from area to pixels...")

        # Load the area raster temporarily to get proper layer reference
        area_layer = QgsRasterLayer(
            flow_accumulation_area_path, "2_Flow_Accumulation_Area")
        if not area_layer.isValid():
            self.log_message(
                "ERROR: Cannot load flow accumulation area layer.")
            return

        # Get base name for layer reference (QGIS uses filename without extension)
        base_name = os.path.splitext(
            os.path.basename(flow_accumulation_area_path))[0]
        # Replace spaces and special characters
        layer_name = base_name.replace(" ", "_")

        # Expression: round("2_flow_accumulation_area@1" / (cellsize * cellsize))
        expression = f'"{layer_name}@1" / ({cell_size} * {cell_size})'

        params_calc = {
            'LAYERS': [flow_accumulation_area_path],
            'EXPRESSION': expression,
            'EXTENT': None,
            'CELL_SIZE': None,
            'CRS': None,
            'OUTPUT': flow_accumulation_pixels_float_path
        }

        result_calc = self.run_processing_algorithm(
            "native:rastercalc", params_calc)

        if not result_calc or not os.path.exists(flow_accumulation_pixels_float_path):
            self.log_message("ERROR: Raster calculator processing failed.")
            return

        self.log_message("Flow accumulation converted to pixels (float).")

        # Step 4: Convert float to integer using GDAL translate
        self.log_message(
            "Converting flow accumulation from float to integer...")

        params_translate = {
            'INPUT': flow_accumulation_pixels_float_path,
            'TARGET_CRS': self.dialog.get_crs(),
            'NODATA': None,
            'COPY_SUBDATASETS': False,
            'OPTIONS': None,
            'EXTRA': '',  # Specify Int32 output type
            'DATA_TYPE': 4,  # Int32 = 4
            'OUTPUT': self.flow_accumulation_path
        }

        result_translate = self.run_processing_algorithm(
            "gdal:translate", params_translate)

        if result_translate and os.path.exists(self.flow_accumulation_path):
            self.load_layer(self.flow_accumulation_path, "2_Flow_Accumulation")
            self.log_message(
                "Flow Accumulation (pixels, integer) processing completed successfully.")

            # Clean up temporary files
            try:
                if os.path.exists(flow_accumulation_pixels_float_path):
                    os.remove(flow_accumulation_pixels_float_path)
                    self.log_message("Cleaned up temporary float raster file.")
            except Exception as e:
                self.log_message(
                    f"Warning: Could not remove temporary file: {e}")
        else:
            self.log_message(
                "ERROR: Failed to convert flow accumulation to integer.")

    def process_gauge_position(self):
        """Process Gauge Position from pour points"""
        # Check prerequisites
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error",
                "Please run Step 1 (Fill DEM) successfully first.")
            return

        if not self.snapped_points_path or not os.path.exists(self.snapped_points_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error",
                "Please run Snap Points successfully first.")
            return

        # Check if gauge position already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.gauge_position_path = os.path.join(
            geometry_folder, "2_gauge_position.tif")
        if self.gauge_position_path and os.path.exists(self.gauge_position_path):
            self.log_message(
                "Gauge Position already exists. Loading existing file...")
            self.load_layer(self.gauge_position_path, "2_Gauge_Position")
            return

        self.log_message("Processing Gauge Position...")

        # Get cell size from filled DEM
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM to get cell size.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM to get cell size.")
            return

        # Calculate cell size (width and height) and get extent
        raster_extent = filled_dem_layer.extent()
        width = filled_dem_layer.width()
        height = filled_dem_layer.height()
        cell_size_width = (raster_extent.xMaximum() -
                           raster_extent.xMinimum()) / width
        cell_size_height = (raster_extent.yMaximum() -
                            raster_extent.yMinimum()) / height

        # Format extent as string: "xmin,xmax,ymin,ymax"
        extent_str = f"{raster_extent.xMinimum()},{raster_extent.xMaximum()},{raster_extent.yMinimum()},{raster_extent.yMaximum()}"

        self.log_message(
            f"Cell size from filled DEM - Width: {cell_size_width:.2f}m, Height: {cell_size_height:.2f}m")
        self.log_message(f"Using extent from filled DEM: {extent_str}")

        params_gauge = {
            'INPUT': self.snapped_points_path,
            'FIELD': 'Name',
            'BURN': 0,
            'USE_Z': False,
            'UNITS': 1,
            'WIDTH': cell_size_width,
            'HEIGHT': cell_size_height,
            'EXTENT': extent_str,
            'NODATA': -9999,
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

        result = self.run_processing_algorithm(
            "native:reprojectlayer", params_reproject)
        if result and os.path.exists(output_path):
            reprojected_layer = QgsVectorLayer(
                output_path, f"{vector_layer.name()}_reprojected", "ogr")
            if reprojected_layer.isValid():
                self.log_message(
                    f"Vector layer reprojected successfully to {target_crs.authid()}")
                return reprojected_layer
            else:
                self.log_message(
                    "ERROR: Reprojected vector layer is not valid.")
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
            QMessageBox.warning(self.dialog, "Input Error",
                                f"Please select a {layer_name} layer.")
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
            self.log_message(
                f"Vector layer detected. Checking CRS for {layer_name}...")

            # Check if vector layer CRS matches input CRS
            layer_crs = layer.crs()
            if input_crs.isValid() and layer_crs.isValid():
                if layer_crs.authid() != input_crs.authid():
                    self.log_message(
                        f"{layer_name} CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                    temp_reprojected_path = os.path.join(
                        geometry_folder, f"temp_{layer_name}_reprojected.gpkg")
                    reprojected_layer = self.reproject_vector_layer(
                        layer, input_crs, temp_reprojected_path)
                    if reprojected_layer:
                        layer = reprojected_layer
                    else:
                        self.log_message(
                            f"WARNING: Reprojection failed. Using original {layer_name} layer.")
                else:
                    self.log_message(
                        f"{layer_name} CRS matches input CRS. No reprojection needed.")

            # Rasterize vector layer using DEM's georeferenced units
            self.log_message(f"Rasterizing {layer_name}...")

            # Get the first field for rasterization (or use a default)
            fields = layer.fields()
            field_name = fields[0].name() if fields else 'id'

            # Calculate width and height from reprojected DEM pixel size
            dem_extent = self.dem_layer.extent()
            width = int((dem_extent.xMaximum() -
                        dem_extent.xMinimum()) / pixel_size_x)
            height = int((dem_extent.yMaximum() -
                         dem_extent.yMinimum()) / pixel_size_y)

            params_rasterize = {
                'INPUT': layer,
                'FIELD': field_name,
                'BURN': 0,
                'USE_Z': False,
                'UNITS': 1,  # Georeferenced units
                'WIDTH': width,
                'HEIGHT': height,
                'EXTENT': extent_str,
                'NODATA': -9999,
                'OPTIONS': None,
                'DATA_TYPE': 5,  # Float32
                'INIT': None,
                'INVERT': False,
                'EXTRA': '',
                'OUTPUT': output_path
            }

            result = self.run_processing_algorithm(
                "gdal:rasterize", params_rasterize)

        elif isinstance(layer, QgsRasterLayer):
            self.log_message(
                f"Raster layer detected. Clipping {layer_name} to DEM extent...")

            # Check if raster layer CRS matches input CRS
            layer_crs = layer.crs()
            if input_crs.isValid() and layer_crs.isValid():
                if layer_crs.authid() != input_crs.authid():
                    self.log_message(
                        f"{layer_name} CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                    temp_reprojected_path = os.path.join(
                        geometry_folder, f"temp_{layer_name}_reprojected.tif")

                    params_warp = {
                        'INPUT': layer.source(),
                        'SOURCE_CRS': None,
                        'TARGET_CRS': input_crs,
                        'RESAMPLING': 0,  # Nearest neighbor
                        'NODATA': -9999,
                        'TARGET_RESOLUTION': None,
                        'OPTIONS': None,
                        'DATA_TYPE': 0,
                        'TARGET_EXTENT': None,
                        'TARGET_EXTENT_CRS': None,
                        'MULTITHREADING': False,
                        'EXTRA': '',
                        'OUTPUT': temp_reprojected_path
                    }

                    warp_result = self.run_processing_algorithm(
                        "gdal:warpreproject", params_warp)
                    if warp_result and os.path.exists(temp_reprojected_path):
                        layer = QgsRasterLayer(
                            temp_reprojected_path, f"{layer_name}_reprojected")
                        self.log_message(
                            f"{layer_name} reprojected successfully.")
                    else:
                        self.log_message(
                            f"WARNING: Reprojection failed. Using original {layer_name} layer.")
                else:
                    self.log_message(
                        f"{layer_name} CRS matches input CRS. No reprojection needed.")

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

            result = self.run_processing_algorithm(
                "gdal:cliprasterbyextent", params_clip)

        else:
            self.log_message(f"ERROR: Unsupported layer type for {layer_name}")
            return None

        if result:
            if display:
                self.load_layer(output_path, layer_name)
            self.log_message(
                f"{layer_name} processing completed successfully.")
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
            self.log_message(
                "No lookup table provided or file not found. Using default mapping.")
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
                            self.log_message(
                                f"WARNING: Unknown type '{type_str}' in lookup table. Skipping.")
                            continue

                        lookup_mapping[grid_value] = type_int
                    except (ValueError, KeyError) as e:
                        self.log_message(
                            f"WARNING: Error parsing lookup table row: {e}")
                        continue

            if lookup_mapping:
                self.log_message(
                    f"Loaded {len(lookup_mapping)} entries from lookup table.")
                return lookup_mapping
            else:
                self.log_message(
                    "Lookup table is empty. Using default mapping.")
                return default_lookup
        except Exception as e:
            self.log_message(
                f"ERROR reading lookup table: {e}. Using default mapping.")
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
        self.log_message(
            "Reclassifying land use raster based on lookup table...")

        # Build reclassification table
        # Format: [min1, max1, value1, min2, max2, value2, ...]
        reclass_table = []
        for grid_value, type_int in sorted(lookup_mapping.items()):
            reclass_table.extend(
                [str(grid_value-0.5), str(grid_value+0.5), str(type_int)])

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

        result = self.run_processing_algorithm(
            "native:reclassifybytable", params_reclass)

        if result and os.path.exists(output_path):
            self.log_message(
                "Land use reclassification completed successfully.")
            return True
        else:
            self.log_message("ERROR: Land use reclassification failed.")
            return False

    def process_land_use(self):
        """Process Land Use layer with clipping and reclassification"""
        layer = self.dialog.mMapLayerComboBox_land_cover.currentLayer()
        if not layer:
            QMessageBox.warning(self.dialog, "Input Error",
                                "Please select a Land Cover layer.")
            return

        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        clipped_path = os.path.join(geometry_folder, "3_land_use_clipped.tif")
        final_path = os.path.join(geometry_folder, "3_land_use.tif")

        # Check if final land use layer already exists
        if os.path.exists(final_path):
            self.log_message(
                "Final land use layer already exists. Loading existing file...")
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
        success = self.reclassify_land_use_raster(
            result, final_path, lookup_mapping)

        if success:
            self.land_use_layer = final_path
            self.load_layer(final_path, "3_Land_Use")
            self.log_message(
                "Land Use layer processed, clipped, and reclassified successfully.")

            # Clean up intermediate clipped file
            try:
                if os.path.exists(clipped_path):
                    os.remove(clipped_path)
            except:
                pass
        else:
            self.log_message("ERROR: Land use reclassification failed.")

    def load_soil_lookup_table(self):
        """
        Load soil lookup table from CSV file using pandas.
        Maps Dominant_S to CLASS (soil class code).

        Returns:
            dict: Mapping of Dominant_S to CLASS, or None if failed
        """
        lookup_file = self.dialog.lineEdit_soil_lookup.text()

        if not lookup_file or not os.path.exists(lookup_file):
            self.log_message(
                "ERROR: Soil lookup table file not found or not specified.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Soil lookup table file.")
            return None

        try:
            import pandas as pd
            
            # Try different encodings
            encodings = ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252']
            df = None
            used_encoding = None
            
            for encoding in encodings:
                try:
                    df = pd.read_csv(lookup_file, encoding=encoding)
                    used_encoding = encoding
                    self.log_message(f"Successfully read file with encoding: {encoding}")
                    break
                except (UnicodeDecodeError, pd.errors.EmptyDataError) as e:
                    continue
            
            if df is None or df.empty:
                self.log_message("ERROR: Could not read file with any encoding or file is empty.")
                QMessageBox.critical(
                    self.dialog, "Error",
                    "Could not read soil lookup table file. Encoding issue or file is empty.")
                return None
            
            # Log column names found
            self.log_message(f"CSV columns found: {', '.join(df.columns.tolist())}")
            
            # Check if required columns exist (case-insensitive)
            column_names_lower = {col.lower(): col for col in df.columns}
            dominant_s_col = column_names_lower.get('dominant_s')
            class_col = column_names_lower.get('class')
            
            if not dominant_s_col:
                self.log_message("ERROR: 'Dominant_S' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'Dominant_S' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            if not class_col:
                self.log_message("ERROR: 'CLASS' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'CLASS' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            self.log_message(f"Using column '{dominant_s_col}' for Dominant_S and '{class_col}' for CLASS")
            
            # Create lookup mapping
            lookup_mapping = {}
            
            # Remove rows with missing values in required columns
            df_clean = df[[dominant_s_col, class_col]].dropna()
            
            # Process each row
            for idx, row in df_clean.iterrows():
                try:
                    # Get Dominant_S value (strip whitespace and convert to string)
                    dominant_s = str(row[dominant_s_col]).strip() if pd.notna(row[dominant_s_col]) else ''
                    
                    # Get CLASS value
                    class_value = row[class_col]
                    
                    if not dominant_s:
                        self.log_message(f"WARNING: Row {idx + 2} (header + 1-based) has empty Dominant_S value. Skipping.")
                        continue
                    
                    if pd.isna(class_value):
                        self.log_message(f"WARNING: Row {idx + 2} has empty CLASS value. Skipping.")
                        continue
                    
                    # Convert CLASS to integer
                    try:
                        class_code = int(float(class_value))  # Use float first to handle numeric strings
                    except (ValueError, TypeError):
                        self.log_message(
                            f"WARNING: Row {idx + 2} has invalid CLASS value '{class_value}'. Skipping.")
                        continue
                    
                    if class_code > 0:
                        lookup_mapping[dominant_s] = class_code
                        if len(lookup_mapping) <= 10:  # Only log first 10 mappings
                            self.log_message(f"  Mapped: {dominant_s} -> {class_code}")
                    else:
                        self.log_message(
                            f"WARNING: Row {idx + 2} has CLASS <= 0 ({class_code}). Skipping.")
                except Exception as e:
                    self.log_message(
                        f"WARNING: Error parsing lookup table row {idx + 2}: {e}")
                    continue

            self.log_message(f"Processed {len(df)} rows from CSV file.")
            
            if lookup_mapping:
                self.log_message(
                    f"Loaded {len(lookup_mapping)} entries from soil lookup table.")
                return lookup_mapping
            else:
                self.log_message(
                    "ERROR: Soil lookup table is empty or could not be parsed.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"Soil lookup table is empty or could not be parsed.\nProcessed {len(df)} rows but found no valid mappings.")
                return None
        except ImportError:
            self.log_message("ERROR: pandas library is not installed.")
            QMessageBox.critical(
                self.dialog, "Error",
                "pandas library is required but not installed.\nPlease install it using: pip install pandas")
            return None
        except Exception as e:
            import traceback
            error_details = traceback.format_exc()
            self.log_message(
                f"ERROR reading soil lookup table: {e}\n{error_details}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error reading soil lookup table:\n{str(e)}")
            return None

    def process_soil(self):
        """Process Soil layer by rasterizing using GDAL with filled DEM extent and resolution"""
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        # Check if filled DEM exists
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error",
                "Please run Fill DEM successfully first.")
            return
        
        # Get soil layer
        layer = self.dialog.mMapLayerComboBox_soil.currentLayer()
        if not layer:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Soil layer.")
            return
        
        # Check if it's a vector layer (rasterization is for vector layers)
        if not isinstance(layer, QgsVectorLayer):
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Soil layer must be a vector layer for rasterization.")
            return
        
        # Check if output already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        output_path = os.path.join(geometry_folder, "3_soil.tif")
        
        if os.path.exists(output_path):
            self.log_message("Soil layer already processed. Loading existing file...")
            self.load_layer(output_path, "3_Soil")
            return
        
        self.log_message("Processing Soil layer...")
        
        # Get input CRS from dialog
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog, "CRS Error",
                "Please set a valid input CRS.")
            return
        
        # Check if soil layer CRS matches input CRS
        layer_crs = layer.crs()
        temp_reprojected_path = None
        layer_to_process = layer
        
        if layer_crs.isValid():
            if layer_crs.authid() != input_crs.authid():
                self.log_message(
                    f"Soil layer CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                temp_reprojected_path = os.path.join(
                    geometry_folder, "temp_soil_reprojected.shp")
                
                reprojected_layer = self.reproject_vector_layer(
                    layer, input_crs, temp_reprojected_path)
                if reprojected_layer:
                    layer_to_process = reprojected_layer
                    self.log_message(
                        f"Soil layer reprojected successfully to {input_crs.authid()}")
                else:
                    self.log_message(
                        "WARNING: Reprojection failed. Using original soil layer.")
                    temp_reprojected_path = None
            else:
                self.log_message(
                    f"Soil layer CRS ({layer_crs.authid()}) matches input CRS. No reprojection needed.")
        else:
            self.log_message(
                "WARNING: Soil layer CRS is not valid. Proceeding with processing...")
        
        # Load soil lookup table
        lookup_mapping = self.load_soil_lookup_table()
        if not lookup_mapping:
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            return
        
        # Check if Dominant_S field exists in the layer
        fields = layer_to_process.fields()
        dominant_s_field = None
        for field_name in ['Dominant_S', 'dominant_s', 'DOMINANT_S']:
            if field_name in [f.name() for f in fields]:
                dominant_s_field = field_name
                break
        
        if not dominant_s_field:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Soil layer must have a 'Dominant_S' attribute field.")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            return
        
        self.log_message(f"Using field '{dominant_s_field}' for soil class lookup")
        
        # Create a temporary layer with CLASS field
        temp_soil_with_class_path = os.path.join(
            geometry_folder, "temp_soil_with_class.shp")
        
        # Add CLASS field and populate it
        self.log_message("Adding CLASS field to soil layer...")
        
        # Create output fields (copy all existing fields + add CLASS)
        output_fields = QgsFields()
        for field in fields:
            output_fields.append(field)
        output_fields.append(QgsField("CLASS", QVariant.Int))
        
        # Remove existing file if it exists
        if os.path.exists(temp_soil_with_class_path):
            try:
                base_name = os.path.splitext(temp_soil_with_class_path)[0]
                for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                    file_to_remove = base_name + ext
                    if os.path.exists(file_to_remove):
                        os.remove(file_to_remove)
            except:
                pass
        
        # Get geometry type from input layer
        geometry_type = layer_to_process.wkbType()
        
        writer = QgsVectorFileWriter(
            temp_soil_with_class_path, "UTF-8", output_fields,
            geometry_type, input_crs, "ESRI Shapefile")
        
        if writer.hasError() != QgsVectorFileWriter.NoError:
            self.log_message(
                f"ERROR creating temporary soil layer: {writer.errorMessage()}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error creating temporary soil layer:\n{writer.errorMessage()}")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            return
        
        # Process features and add CLASS value
        features_processed = 0
        features_mapped = 0
        features_unmapped = 0
        
        for feature in layer_to_process.getFeatures():
            new_feat = QgsFeature(output_fields)
            new_feat.setGeometry(feature.geometry())
            
            # Copy all attributes
            attrs = feature.attributes()
            for i, attr in enumerate(attrs):
                if i < len(output_fields) - 1:  # Exclude the last field (CLASS)
                    new_feat.setAttribute(i, attr)
            
            # Get Dominant_S value and lookup CLASS
            dominant_s_value = feature.attribute(dominant_s_field)
            class_code = None
            
            if dominant_s_value:
                dominant_s_str = str(dominant_s_value).strip()
                class_code = lookup_mapping.get(dominant_s_str)
            
            if class_code is not None:
                new_feat.setAttribute("CLASS", class_code)
                features_mapped += 1
            else:
                # If no mapping found, use 0 or a default value
                self.log_message(
                    f"WARNING: No CLASS found for Dominant_S='{dominant_s_value}'. Using 0.")
                new_feat.setAttribute("CLASS", 0)
                features_unmapped += 1
            
            writer.addFeature(new_feat)
            features_processed += 1
        
        del writer
        
        self.log_message(
            f"Processed {features_processed} features. "
            f"Mapped: {features_mapped}, Unmapped: {features_unmapped}")
        
        if features_mapped == 0:
            QMessageBox.warning(
                self.dialog, "Mapping Error",
                "No features were successfully mapped to CLASS values. Please check your lookup table.")
            # Clean up temporary files
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            if os.path.exists(temp_soil_with_class_path):
                try:
                    base_name = os.path.splitext(temp_soil_with_class_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            return
        
        # Get filled DEM layer for extent and resolution
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM layer.")
            # Clean up temporary files
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            if os.path.exists(temp_soil_with_class_path):
                try:
                    base_name = os.path.splitext(temp_soil_with_class_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            return
        
        # Get extent and dimensions from filled DEM
        raster_extent = filled_dem_layer.extent()
        width = filled_dem_layer.width()  # Pixel width
        height = filled_dem_layer.height()  # Pixel height

        cell_size_width = (raster_extent.xMaximum() -
                           raster_extent.xMinimum()) / width
        cell_size_height = (raster_extent.yMaximum() -
                            raster_extent.yMinimum()) / height
        
        # Format extent as string: "xmin,xmax,ymin,ymax [EPSG:xxxxx]"
        extent_str = f"{raster_extent.xMinimum()},{raster_extent.xMaximum()},{raster_extent.yMinimum()},{raster_extent.yMaximum()} [{input_crs.authid()}]"
        
        self.log_message(f"Using filled DEM extent: {extent_str}")
        self.log_message(f"Using filled DEM resolution: {width} x {height} pixels")
        
        # Rasterize using GDAL with CLASS field
        params_rasterize = {
            'INPUT': temp_soil_with_class_path,
            'FIELD': 'CLASS',
            'BURN': 0,
            'USE_Z': False,
            'UNITS': 1,  # Georeferenced units
            'WIDTH': cell_size_width,
            'HEIGHT': cell_size_height,
            'EXTENT': extent_str,   
            'NODATA': -9999,
            'OPTIONS': None,
            'DATA_TYPE': 5,  # Int32
            'INIT': None,
            'INVERT': False,
            'EXTRA': '',
            'OUTPUT': output_path
        }
        
        result = self.run_processing_algorithm("gdal:rasterize", params_rasterize)
        if result:
            self.load_layer(output_path, "3_Soil")
            self.log_message("Soil layer rasterized successfully.")
        else:
            self.log_message("ERROR: Soil rasterization failed.")
        
        # Clean up temporary files
        if temp_reprojected_path and os.path.exists(temp_reprojected_path):
            try:
                base_name = os.path.splitext(temp_reprojected_path)[0]
                for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                    file_to_remove = base_name + ext
                    if os.path.exists(file_to_remove):
                        os.remove(file_to_remove)
                self.log_message("Cleaned up temporary reprojected soil file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary reprojected file: {e}")
        
        if os.path.exists(temp_soil_with_class_path):
            try:
                base_name = os.path.splitext(temp_soil_with_class_path)[0]
                for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                    file_to_remove = base_name + ext
                    if os.path.exists(file_to_remove):
                        os.remove(file_to_remove)
                self.log_message("Cleaned up temporary soil with class file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary soil with class file: {e}")

    def process_geology(self):
        """Process Geology layer by rasterizing using GDAL with filled DEM extent and resolution"""
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        # Check if filled DEM exists
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error",
                "Please run Fill DEM successfully first.")
            return
        
        # Get geology layer
        layer = self.dialog.mMapLayerComboBox_geology.currentLayer()
        if not layer:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Geology layer.")
            return
        
        # Check if it's a vector layer (rasterization is for vector layers)
        if not isinstance(layer, QgsVectorLayer):
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Geology layer must be a vector layer for rasterization.")
            return
        
        # Check if output already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.geology_path = os.path.join(geometry_folder, "3_geology_processed.tif")
        
        if os.path.exists(self.geology_path):
            self.log_message("Geology layer already processed. Loading existing file...")
            self.load_layer(self.geology_path, "Geology")
            return
        
        self.log_message("Processing Geology layer...")
        
        # Get input CRS from dialog
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog, "CRS Error",
                "Please set a valid input CRS.")
            return
        
        # Check if geology layer CRS matches input CRS
        layer_crs = layer.crs()
        temp_reprojected_path = None
        layer_to_rasterize = layer
        
        if layer_crs.isValid():
            if layer_crs.authid() != input_crs.authid():
                self.log_message(
                    f"Geology layer CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                temp_reprojected_path = os.path.join(
                    geometry_folder, "temp_geology_reprojected.shp")
                
                reprojected_layer = self.reproject_vector_layer(
                    layer, input_crs, temp_reprojected_path)
                if reprojected_layer:
                    layer_to_rasterize = reprojected_layer
                    self.log_message(
                        f"Geology layer reprojected successfully to {input_crs.authid()}")
                else:
                    self.log_message(
                        "WARNING: Reprojection failed. Using original geology layer.")
                    temp_reprojected_path = None
            else:
                self.log_message(
                    f"Geology layer CRS ({layer_crs.authid()}) matches input CRS. No reprojection needed.")
        else:
            self.log_message(
                "WARNING: Geology layer CRS is not valid. Proceeding with rasterization...")
        
        # Get filled DEM layer for extent and resolution
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM layer.")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    os.remove(temp_reprojected_path)
                except:
                    pass
            return
        
        # Get extent and dimensions from filled DEM
        raster_extent = filled_dem_layer.extent()
        width = filled_dem_layer.width()  # Pixel width
        height = filled_dem_layer.height()  # Pixel height

        cell_size_width = (raster_extent.xMaximum() -
                           raster_extent.xMinimum()) / width
        cell_size_height = (raster_extent.yMaximum() -
                            raster_extent.yMinimum()) / height
        
        # Format extent as string: "xmin,xmax,ymin,ymax [EPSG:xxxxx]"
        extent_str = f"{raster_extent.xMinimum()},{raster_extent.xMaximum()},{raster_extent.yMinimum()},{raster_extent.yMaximum()} [{input_crs.authid()}]"
        
        self.log_message(f"Using filled DEM extent: {extent_str}")
        self.log_message(f"Using filled DEM resolution: {width} x {height} pixels")
        
        # Get field name for rasterization (try common field names, or use first field)
        fields = layer_to_rasterize.fields()
        field_name = None
        # Try common geology field names
        for field_candidate in ['GEO_CLASS', 'geology', 'Geology', 'GEO', 'CLASS', 'class', 'id']:
            if field_candidate in [f.name() for f in fields]:
                field_name = field_candidate
                break
        
        # If no common field found, use first field
        if not field_name and fields:
            field_name = fields[0].name()
        
        if not field_name:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Geology layer must have at least one attribute field for rasterization.")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    os.remove(temp_reprojected_path)
                except:
                    pass
            return
        
        self.log_message(f"Using field '{field_name}' for rasterization")
        
        # Get layer source path (for vector layers)
        if temp_reprojected_path and os.path.exists(temp_reprojected_path):
            # Use the reprojected layer path
            layer_input = temp_reprojected_path
        else:
            layer_source = layer_to_rasterize.source()
            if not layer_source or not os.path.exists(layer_source):
                # If source is not a file path, use the layer directly
                layer_input = layer_to_rasterize
            else:
                layer_input = layer_source
        
        # Rasterize using GDAL
        params_rasterize = {
            'INPUT': layer_input,
            'FIELD': field_name,
            'BURN': 0,
            'USE_Z': False,
            'UNITS': 1,  # Georeferenced units
            'WIDTH': cell_size_width,
            'HEIGHT': cell_size_height,
            'EXTENT': extent_str,   
            'NODATA': -9999,
            'OPTIONS': None,
            'DATA_TYPE': 5,  # Int32
            'INIT': None,
            'INVERT': False,
            'EXTRA': '',
            'OUTPUT': self.geology_path
        }
        
        result = self.run_processing_algorithm("gdal:rasterize", params_rasterize)
        if result:
            self.load_layer(self.geology_path, "Geology")
            self.log_message("Geology layer rasterized successfully.")
        else:
            self.log_message("ERROR: Geology rasterization failed.")
            self.geology_path = None
        
        # Clean up temporary reprojected file if it was created
        if temp_reprojected_path and os.path.exists(temp_reprojected_path):
            try:
                # Remove shapefile and associated files (.shp, .shx, .dbf, .prj, etc.)
                base_name = os.path.splitext(temp_reprojected_path)[0]
                for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                    file_to_remove = base_name + ext
                    if os.path.exists(file_to_remove):
                        os.remove(file_to_remove)
                self.log_message("Cleaned up temporary reprojected geology file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary file: {e}")
    
    def mask_all_layers(self):
        """
        Mask all raster layers using the merged watershed vector.
        Reads L0, L1, L2 values from UI and calculates extent based on L2.
        """
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error",
                "Please run Fill DEM successfully first.")
            return
        
        # Check for merged watershed mask
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        merged_watershed_path = os.path.join(
            geometry_folder, "4_watershed_merged_vector.shp")
        
        if not os.path.exists(merged_watershed_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error",
                "Please run Watershed Delineation successfully first to create merged watershed layer.")
            return
        
        # Read L0, L1, L2 values and units from UI
        try:
            l0_value = float(self.dialog.lineEdit_L0.text())
            l0_unit = self.dialog.comboBox_L0.currentText()
            l1_value = float(self.dialog.lineEdit_L1.text())
            l1_unit = self.dialog.comboBox_L1.currentText()
            l2_value = float(self.dialog.lineEdit_L2.text())
            l2_unit = self.dialog.comboBox_L2.currentText()
        except (ValueError, AttributeError) as e:
            QMessageBox.warning(
                self.dialog, "Input Error",
                f"Please enter valid numeric values for L0, L1, and L2.\nError: {e}")
            return
        
        if not l0_unit or not l1_unit or not l2_unit:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select units (m or °) for L0, L1, and L2.")
            return
        
        self.log_message(f"\n--- Masking all layers ---")
        self.log_message(f"L0: {l0_value} {l0_unit}")
        self.log_message(f"L1: {l1_value} {l1_unit}")
        self.log_message(f"L2: {l2_value} {l2_unit}")
        
        # Get filled DEM layer
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM layer.")
            return
        
        dem_extent = filled_dem_layer.extent()
        dem_crs = filled_dem_layer.crs()
        input_crs = self.dialog.get_crs()
        
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog, "CRS Error",
                "Please set a valid input CRS.")
            return
        
        self.log_message(f"Filled DEM CRS: {dem_crs.authid()}")
        self.log_message(
            f"Original DEM extent: {dem_extent.xMinimum():.2f}, "
            f"{dem_extent.xMaximum():.2f}, {dem_extent.yMinimum():.2f}, "
            f"{dem_extent.yMaximum():.2f}")
        
        # Calculate extent based on L2
        # For ERA5 grid snapping: if L2 is in degrees, do all operations in WGS84 degrees
        if l2_unit == "°":
            # Transform extent to WGS84 to snap to ERA5 grid
            wgs84_crs = QgsCoordinateReferenceSystem('EPSG:4326')
            
            # Transform extent to WGS84
            transform_to_wgs84 = QgsCoordinateTransform(
                dem_crs, wgs84_crs, QgsProject.instance())
            transform_to_wgs84.setBallparkTransformsAreAppropriate(True)
            dem_extent_wgs84 = transform_to_wgs84.transform(dem_extent)
            
            self.log_message(
                f"DEM extent in WGS84: {dem_extent_wgs84.xMinimum():.6f}, "
                f"{dem_extent_wgs84.xMaximum():.6f}, {dem_extent_wgs84.yMinimum():.6f}, "
                f"{dem_extent_wgs84.yMaximum():.6f}")
            
            # Round to L2 grid in degrees (snap to ERA5 grid)
            xmin_rounded = math.floor(dem_extent_wgs84.xMinimum() / l2_value) * l2_value
            xmax_rounded = math.ceil(dem_extent_wgs84.xMaximum() / l2_value) * l2_value
            ymin_rounded = math.floor(dem_extent_wgs84.yMinimum() / l2_value) * l2_value
            ymax_rounded = math.ceil(dem_extent_wgs84.yMaximum() / l2_value) * l2_value
            
            # Shift L2/2 inside on all 4 directions (in degrees)
            half_l2 = l2_value / 2.0
            xmin_shifted = xmin_rounded + half_l2
            xmax_shifted = xmax_rounded - half_l2
            ymin_shifted = ymin_rounded + half_l2
            ymax_shifted = ymax_rounded - half_l2
            
            # Validate shifted extent
            if xmin_shifted >= xmax_shifted or ymin_shifted >= ymax_shifted:
                QMessageBox.warning(
                    self.dialog, "Extent Error",
                    f"Invalid extent after shifting. L2 ({l2_value}°) may be too large relative to DEM extent.")
                return
            
            # Create shifted extent rectangle in WGS84
            shifted_extent_wgs84 = QgsRectangle(
                xmin_shifted, ymin_shifted, xmax_shifted, ymax_shifted)
            
            self.log_message(
                f"Extent rounded to L2 grid ({l2_value}°) and shifted L2/2 in WGS84: "
                f"{xmin_shifted:.6f}, {xmax_shifted:.6f}, {ymin_shifted:.6f}, {ymax_shifted:.6f}")
            
            # Transform back to input CRS
            transform_back = QgsCoordinateTransform(
                wgs84_crs, input_crs, QgsProject.instance())
            transform_back.setBallparkTransformsAreAppropriate(True)
            
            # Transform the four corners individually to ensure valid transformation
            corners = [
                QgsPointXY(xmin_shifted, ymin_shifted),  # Bottom-left
                QgsPointXY(xmax_shifted, ymin_shifted),  # Bottom-right
                QgsPointXY(xmax_shifted, ymax_shifted),  # Top-right
                QgsPointXY(xmin_shifted, ymax_shifted)   # Top-left
            ]
            
            transformed_corners = []
            for corner in corners:
                try:
                    transformed_corner = transform_back.transform(corner)
                    transformed_corners.append(transformed_corner)
                except Exception as e:
                    self.log_message(f"WARNING: Corner transformation failed: {e}")
                    # Fallback: use geometry transform
                    geom = QgsGeometry.fromPointXY(corner)
                    geom.transform(transform_back)
                    transformed_corners.append(geom.asPoint())
            
            # Create extent from transformed corners
            x_coords = [p.x() for p in transformed_corners]
            y_coords = [p.y() for p in transformed_corners]
            xmin_final = min(x_coords)
            xmax_final = max(x_coords)
            ymin_final = min(y_coords)
            ymax_final = max(y_coords)
            
            extent_str = f"{xmin_final},{xmax_final},{ymin_final},{ymax_final} [{input_crs.authid()}]"
            
            self.log_message(
                f"Extent transformed back to input CRS ({input_crs.authid()}): {extent_str}")
        else:
            # L2 is in meters - round directly in input CRS
            # Round to L2 value
            xmin_rounded = math.floor(dem_extent.xMinimum() / l2_value) * l2_value
            xmax_rounded = math.ceil(dem_extent.xMaximum() / l2_value) * l2_value
            ymin_rounded = math.floor(dem_extent.yMinimum() / l2_value) * l2_value
            ymax_rounded = math.ceil(dem_extent.yMaximum() / l2_value) * l2_value
            
            # Shift L2/2 inside on all 4 directions
            half_l2 = l2_value / 2.0
            xmin_shifted = xmin_rounded + half_l2
            xmax_shifted = xmax_rounded - half_l2
            ymin_shifted = ymin_rounded + half_l2
            ymax_shifted = ymax_rounded - half_l2
            
            # Validate shifted extent
            if xmin_shifted >= xmax_shifted or ymin_shifted >= ymax_shifted:
                QMessageBox.warning(
                    self.dialog, "Extent Error",
                    f"Invalid extent after shifting. L2 ({l2_value} m) may be too large relative to DEM extent.")
                return
            
            extent_str = f"{xmin_shifted},{xmax_shifted},{ymin_shifted},{ymax_shifted} [{input_crs.authid()}]"
            
            self.log_message(
                f"Extent rounded to L2 ({l2_value} m) and shifted L2/2 "
                f"({half_l2:.2f} m) inside in input CRS: {extent_str}")
        
        # Convert L0 to meters if in degrees
        l0_meters = l0_value
        if l0_unit == "°":
            # Convert degrees to meters (approximate: 1 degree ≈ 111,000 meters)
            # For more accuracy, could use cos(lat) for longitude, but using average
            center_lat = (dem_extent.yMinimum() + dem_extent.yMaximum()) / 2.0
            l0_meters = l0_value * 111000
            self.log_message(
                f"L0 converted from degrees to meters: {l0_meters:.2f} m "
                f"(using center latitude: {center_lat:.4f})")
        
        x_resolution = l0_meters
        y_resolution = l0_meters
        
        self.log_message(f"X Resolution (L0): {x_resolution:.2f} m")
        self.log_message(f"Y Resolution (L0): {y_resolution:.2f} m")
        
        # Get list of raster layers to mask
        layers_to_mask = []
        
        # Check for filled DEM
        if self.filled_dem_path and os.path.exists(self.filled_dem_path):
            layers_to_mask.append({
                'input': self.filled_dem_path,
                'output': os.path.join(geometry_folder, "1_dem_filled_masked.tif"),
                'name': "1_DEM_Filled_Masked"
            })
        
        # Check for other processed layers
        layer_mapping = [
            ('aspect_path', '1_dem_aspect_masked.tif', '1_DEM_Aspect_Masked'),
            ('slope_path', '1_dem_slope_masked.tif', '1_DEM_Slope_Masked'),
            ('flow_accumulation_path', '2_flow_accumulation_masked.tif', '2_Flow_Accumulation_Masked'),
            ('flow_direction_path', '2_flow_direction_masked.tif', '2_Flow_Direction_Masked'),
            ('gauge_position_path', '2_gauge_position_masked.tif', '2_Gauge_Position_Masked'),
        ]
        
        for attr_name, output_filename, layer_name in layer_mapping:
            layer_path = getattr(self, attr_name, None)
            if layer_path and os.path.exists(layer_path):
                layers_to_mask.append({
                    'input': layer_path,
                    'output': os.path.join(geometry_folder, output_filename),
                    'name': layer_name
                })
        
        # Check for land cover layer
        if self.land_use_layer and os.path.exists(self.land_use_layer):
            layers_to_mask.append({
                'input': self.land_use_layer,
                'output': os.path.join(geometry_folder, "3_land_use_masked.tif"),
                'name': "3_Land_Use_Masked"
            })
        
        # Check for soil layer
        soil_path = os.path.join(geometry_folder, "3_soil.tif")
        if os.path.exists(soil_path):
            layers_to_mask.append({
                'input': soil_path,
                'output': os.path.join(geometry_folder, "3_soil_masked.tif"),
                'name': "3_Soil_Masked"
            })
        
        # Check for geology layer
        if self.geology_path and os.path.exists(self.geology_path):
            layers_to_mask.append({
                'input': self.geology_path,
                'output': os.path.join(geometry_folder, "3_geology_processed_masked.tif"),
                'name': "3_Geology_Masked"
            })
        
        if not layers_to_mask:
            QMessageBox.warning(
                self.dialog, "No Layers",
                "No raster layers found to mask. Please process at least the filled DEM first.")
            return
        
        self.log_message(f"Masking {len(layers_to_mask)} raster layer(s)...")
        
        # Mask each layer
        for layer_info in layers_to_mask:
            input_path = layer_info['input']
            output_path = layer_info['output']
            layer_name = layer_info['name']
            
            self.log_message(f"Masking {layer_name}...")
            
            # Build params - only include TARGET_EXTENT if we have a valid extent string
            params_mask = {
                'INPUT': input_path,
                'MASK': merged_watershed_path,
                'SOURCE_CRS': input_crs,
                'TARGET_CRS': input_crs,
                'NODATA': None,
                'ALPHA_BAND': False,
                'CROP_TO_CUTLINE': False,
                'KEEP_RESOLUTION': False,
                'SET_RESOLUTION': True,
                'X_RESOLUTION': x_resolution,
                'Y_RESOLUTION': y_resolution,
                'MULTITHREADING': False,
                'OPTIONS': None,
                'DATA_TYPE': 0,
                'EXTRA': '',
                'OUTPUT': output_path
            }
            
            # Only add TARGET_EXTENT if we have a valid extent string
            # GDAL will use TARGET_CRS for the extent automatically
            if extent_str:
                params_mask['TARGET_EXTENT'] = extent_str
            
            result = self.run_processing_algorithm(
                "gdal:cliprasterbymasklayer", params_mask)
            
            if result and os.path.exists(output_path):
                self.load_layer(output_path, layer_name)
                self.log_message(f"{layer_name} masked successfully.")
            else:
                self.log_message(f"ERROR: Failed to mask {layer_name}.")
        
        self.log_message("Masking process completed.")
    
    def write_header_file(self, L_info, output_path):
        """
        Helper function to write header.txt file for a given L level.
        
        Args:
            L_info: Dictionary with keys: ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value
            output_path: Directory path where header.txt should be created
        """
        header_file = os.path.join(output_path, "header.txt")
        
        try:
            with open(header_file, 'w', encoding='utf-8') as f:
                f.write(f"ncols\t{L_info['ncols']}\n")
                f.write(f"nrows\t{L_info['nrows']}\n")
                f.write(f"xllcorner\t{L_info['xllcorner']}\n")
                f.write(f"yllcorner\t{L_info['yllcorner']}\n")
                f.write(f"cellsize\t{L_info['cellsize']}\n")
                f.write(f"NODATA_value\t{L_info['NODATA_value']}\n")
                f.write("\n")
            
            self.log_message(f"Header file written: {header_file}")
            return True
        except Exception as e:
            self.log_message(f"ERROR: Failed to write header file {header_file}: {e}")
            return False
    
    def process_lat_lon(self):
        """
        Process lat/lon information by reading dem_masked_layer and creating header files
        for L0, L1, L11, and L2 levels. This prepares data for latlon.nc file creation.
        """
        self.log_message("\n--- Processing Lat/Lon Headers ---")
        
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        dem_masked_path = os.path.join(geometry_folder, "1_dem_filled_masked.tif")
        
        if not os.path.exists(dem_masked_path):
            QMessageBox.warning(
                self.dialog, "Dependency Error",
                "Please run 'Mask All Layers' successfully first to create masked DEM.")
            return
        
        # Read L0, L1, L2 values and units from UI
        try:
            l0_value = float(self.dialog.lineEdit_L0.text())
            l0_unit = self.dialog.comboBox_L0.currentText()
            l1_value = float(self.dialog.lineEdit_L1.text())
            l1_unit = self.dialog.comboBox_L1.currentText()
            l2_value = float(self.dialog.lineEdit_L2.text())
            l2_unit = self.dialog.comboBox_L2.currentText()
        except (ValueError, AttributeError) as e:
            QMessageBox.warning(
                self.dialog, "Input Error",
                f"Please enter valid numeric values for L0, L1, and L2.\nError: {e}")
            return
        
        if not l0_unit or not l1_unit or not l2_unit:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select units (m or °) for L0, L1, and L2.")
            return
        
        # Read the masked DEM layer
        dem_masked_layer = QgsRasterLayer(dem_masked_path, "DEM_Masked")
        if not dem_masked_layer.isValid():
            self.log_message("ERROR: Cannot read masked DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read masked DEM layer.")
            return
        
        # Get extent and dimensions from masked DEM
        raster_extent = dem_masked_layer.extent()
        width = dem_masked_layer.width()  # Number of columns
        height = dem_masked_layer.height()  # Number of rows
        
        # Calculate cell size (assuming square cells)
        cell_size_x = (raster_extent.xMaximum() - raster_extent.xMinimum()) / width
        cell_size_y = (raster_extent.yMaximum() - raster_extent.yMinimum()) / height
        cell_size = (cell_size_x + cell_size_y) / 2.0  # Average cell size
        
        # Get corner coordinates
        xllcorner = raster_extent.xMinimum()
        yllcorner = raster_extent.yMinimum()
        
        # L0: Use values from masked DEM
        self.L0 = {
            'ncols': width,
            'nrows': height,
            'xllcorner': xllcorner,
            'yllcorner': yllcorner,
            'cellsize': cell_size,
            'NODATA_value': -9999
        }
        
        self.log_message(f"L0 information extracted from masked DEM:")
        self.log_message(f"  ncols: {self.L0['ncols']}, nrows: {self.L0['nrows']}")
        self.log_message(f"  xllcorner: {self.L0['xllcorner']}, yllcorner: {self.L0['yllcorner']}")
        self.log_message(f"  cellsize: {self.L0['cellsize']:.2f}")
        
        # Convert L1 cell size to meters if needed
        l1_cellsize = l1_value
        if l1_unit == "°":
            # Convert degrees to meters (approximate: 1 degree ≈ 111,000 meters)
            center_lat = (raster_extent.yMinimum() + raster_extent.yMaximum()) / 2.0
            l1_cellsize = l1_value * 111000
            self.log_message(f"L1 cell size converted from {l1_value}° to {l1_cellsize:.2f} m")
        
        # L1: Calculate new ncols and nrows based on L1 cell size
        # Keep same xllcorner and yllcorner
        x_range = raster_extent.xMaximum() - raster_extent.xMinimum()
        y_range = raster_extent.yMaximum() - raster_extent.yMinimum()
        
        l1_ncols = int(round(x_range / l1_cellsize))
        l1_nrows = int(round(y_range / l1_cellsize))
        
        self.L1 = {
            'ncols': l1_ncols,
            'nrows': l1_nrows,
            'xllcorner': xllcorner,
            'yllcorner': yllcorner,
            'cellsize': l1_cellsize,
            'NODATA_value': -9999
        }
        
        self.log_message(f"L1 information calculated:")
        self.log_message(f"  ncols: {self.L1['ncols']}, nrows: {self.L1['nrows']}")
        self.log_message(f"  xllcorner: {self.L1['xllcorner']}, yllcorner: {self.L1['yllcorner']}")
        self.log_message(f"  cellsize: {self.L1['cellsize']:.2f}")
        
        # L11: Exactly same as L1
        self.L11 = {
            'ncols': l1_ncols,
            'nrows': l1_nrows,
            'xllcorner': xllcorner,
            'yllcorner': yllcorner,
            'cellsize': l1_cellsize,
            'NODATA_value': -9999
        }
        
        self.log_message(f"L11 information (same as L1):")
        self.log_message(f"  ncols: {self.L11['ncols']}, nrows: {self.L11['nrows']}")
        self.log_message(f"  cellsize: {self.L11['cellsize']:.2f}")
        
        # Convert L2 cell size to meters if needed
        l2_cellsize = l2_value
        if l2_unit == "°":
            # Convert degrees to meters
            l2_cellsize = l2_value * 111000
            self.log_message(f"L2 cell size converted from {l2_value}° to {l2_cellsize:.2f} m")
        
        # L2: Calculate new ncols and nrows based on L2 cell size
        l2_ncols = int(round(x_range / l2_cellsize))
        l2_nrows = int(round(y_range / l2_cellsize))
        
        self.L2 = {
            'ncols': l2_ncols,
            'nrows': l2_nrows,
            'xllcorner': xllcorner,
            'yllcorner': yllcorner,
            'cellsize': l2_cellsize,
            'NODATA_value': -9999
        }
        
        self.log_message(f"L2 information calculated:")
        self.log_message(f"  ncols: {self.L2['ncols']}, nrows: {self.L2['nrows']}")
        self.log_message(f"  xllcorner: {self.L2['xllcorner']}, yllcorner: {self.L2['yllcorner']}")
        self.log_message(f"  cellsize: {self.L2['cellsize']:.2f}")
        
        # Create header files in input/latlon folder structure
        input_folder = os.path.join(self.dialog.project_folder, "input")
        latlon_folder = os.path.join(input_folder, "latlon")
        os.makedirs(latlon_folder, exist_ok=True)
        
        # Create NetCDF file
        self.log_message("\nCreating latlon.nc file...")
        self.create_latlon_nc_file(latlon_folder)
        
        self.log_message("Lat/Lon header processing completed successfully.")
        QMessageBox.information(
            self.dialog, "Success",
            "Lat/Lon headers and NetCDF file created successfully in input/latlon folder.")
    
    def create_latlon_nc_file(self, latlon_folder):
        """
        Create latlon.nc NetCDF file with L0, L1, and L11 information.
        Includes xc/yc coordinates (1D arrays in UTM meters) and lat/lon values (2D arrays) for each level.
        """
        try:
            import numpy as np
            import xarray as xr
            from pyproj import Transformer
            import time
            
            # Get input CRS for coordinate transformation
            input_crs = self.dialog.get_crs()
            if not input_crs.isValid():
                self.log_message("ERROR: Invalid input CRS. Cannot create latlon.nc file.")
                return False
            
            # Get CRS string for pyproj
            crs_string = None
            if input_crs.postgisSrid():
                crs_string = f"EPSG:{input_crs.postgisSrid()}"
            else:
                crs_string = input_crs.authid()
            
            if not crs_string:
                self.log_message("ERROR: Could not determine CRS string. Cannot create latlon.nc file.")
                return False
            
            # Create transformer from input CRS to WGS84 (EPSG:4326)
            transformer = Transformer.from_crs(crs_string, "EPSG:4326", always_xy=True)
            
            # Helper function to calculate xc/yc arrays (1D) for a given level
            def calculate_coordinate_arrays(L_info):
                """
                Calculate xc and yc 1D arrays for grid cell centers in UTM meters.
                xc represents column centers, yc represents row centers.
                """
                ncols = L_info['ncols']
                nrows = L_info['nrows']
                xllcorner = L_info['xllcorner']
                yllcorner = L_info['yllcorner']
                cellsize = L_info['cellsize']
                
                # Calculate xc (column centers) - 1D array
                # xc = xllcorner + (column_index + 0.5) * cellsize
                xc = np.array([xllcorner + (i + 0.5) * cellsize for i in range(ncols)], dtype=np.float32)
                
                # Calculate yc (row centers) - 1D array
                # yc = yllcorner + (row_index + 0.5) * cellsize
                # Note: yllcorner is bottom-left, so we go from bottom to top
                yc = np.array([yllcorner + (i + 0.5) * cellsize for i in range(nrows)], dtype=np.float32)
                
                return xc, yc
            
            # Helper function to calculate lat/lon 2D arrays
            def calculate_latlon_arrays(L_info, xc, yc):
                """
                Calculate lat and lon 2D arrays from xc/yc coordinates.
                Creates a meshgrid and transforms UTM coordinates to lat/lon.
                """
                ncols = L_info['ncols']
                nrows = L_info['nrows']
                
                # Create meshgrid of xc and yc
                x_grid, y_grid = np.meshgrid(xc, yc)
                
                # Transform to lat/lon (transformer expects x, y and returns lon, lat)
                lon_array, lat_array = transformer.transform(x_grid.flatten(), y_grid.flatten())
                
                # Reshape to 2D (nrows, ncols)
                lon_2d = np.array(lon_array, dtype=np.float64).reshape(nrows, ncols)
                lat_2d = np.array(lat_array, dtype=np.float64).reshape(nrows, ncols)
                
                return lat_2d, lon_2d
            
            # Calculate coordinates for each level
            self.log_message("Calculating coordinate arrays for L0...")
            xc_l0, yc_l0 = calculate_coordinate_arrays(self.L0)
            lat_l0, lon_l0 = calculate_latlon_arrays(self.L0, xc_l0, yc_l0)
            
            self.log_message("Calculating coordinate arrays for L1...")
            xc_l1, yc_l1 = calculate_coordinate_arrays(self.L1)
            lat_l1, lon_l1 = calculate_latlon_arrays(self.L1, xc_l1, yc_l1)
            
            self.log_message("Calculating coordinate arrays for L11...")
            xc_l11, yc_l11 = calculate_coordinate_arrays(self.L11)
            lat_l11, lon_l11 = calculate_latlon_arrays(self.L11, xc_l11, yc_l11)
            
            # Create xarray Dataset
            self.log_message("Creating NetCDF dataset...")
            
            # Create dimensions and coordinate variables
            ds = xr.Dataset(
                {
                    # L0 variables
                    'xc_l0': (['xc_l0'], xc_l0),
                    'yc_l0': (['yc_l0'], yc_l0),
                    'lon_l0': (['yc_l0', 'xc_l0'], lon_l0),
                    'lat_l0': (['yc_l0', 'xc_l0'], lat_l0),
                    
                    # L1 variables
                    'xc': (['xc'], xc_l1),
                    'yc': (['yc'], yc_l1),
                    'lon': (['yc', 'xc'], lon_l1),
                    'lat': (['yc', 'xc'], lat_l1),
                    
                    # L11 variables (same as L1)
                    'xc_l11': (['xc_l11'], xc_l11),
                    'yc_l11': (['yc_l11'], yc_l11),
                    'lon_l11': (['yc_l11', 'xc_l11'], lon_l11),
                    'lat_l11': (['yc_l11', 'xc_l11'], lat_l11),
                },
                attrs={
                    'description': 'lat lon file',
                    'projection': crs_string.lower(),
                    'history': f'Created {time.ctime()}'
                }
            )
            
            # Set variable attributes
            # L0 attributes
            ds['xc_l0'].attrs['axis'] = 'X'
            ds['yc_l0'].attrs['axis'] = 'Y'
            ds['lon_l0'].attrs['units'] = 'degrees_east'
            ds['lon_l0'].attrs['long_name'] = 'longitude at level 0'
            ds['lat_l0'].attrs['units'] = 'degrees_north'
            ds['lat_l0'].attrs['long_name'] = 'latitude at level 0'
            
            # L1 attributes
            ds['xc'].attrs['axis'] = 'X'
            ds['yc'].attrs['axis'] = 'Y'
            ds['lon'].attrs['units'] = 'degrees_east'
            ds['lon'].attrs['long_name'] = 'longitude at level 1'
            ds['lat'].attrs['units'] = 'degrees_north'
            ds['lat'].attrs['long_name'] = 'latitude at level 1'
            
            # L11 attributes
            ds['xc_l11'].attrs['axis'] = 'X'
            ds['yc_l11'].attrs['axis'] = 'Y'
            ds['lon_l11'].attrs['units'] = 'degrees_east'
            ds['lon_l11'].attrs['long_name'] = 'longitude at level 11'
            ds['lat_l11'].attrs['units'] = 'degrees_north'
            ds['lat_l11'].attrs['long_name'] = 'latitude at level 11'
            
            # Set encoding for compression
            encoding = {
                'lon_l0': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L0['nrows'], self.L0['ncols'])
                },
                'lat_l0': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L0['nrows'], self.L0['ncols'])
                },
                'lon': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L1['nrows'], self.L1['ncols'])
                },
                'lat': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L1['nrows'], self.L1['ncols'])
                },
                'lon_l11': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L11['nrows'], self.L11['ncols'])
                },
                'lat_l11': {
                    'zlib': True,
                    'complevel': 9,
                    '_FillValue': -9999.0,
                    'chunksizes': (self.L11['nrows'], self.L11['ncols'])
                },
            }
            
            # Save to NetCDF file
            output_file = os.path.join(latlon_folder, "latlon_1.nc")
            ds.to_netcdf(output_file, encoding=encoding)
            
            self.log_message(f"NetCDF file created successfully: {output_file}")
            return True
            
        except ImportError as e:
            self.log_message(f"ERROR: Required library not available: {e}")
            QMessageBox.warning(
                self.dialog, "Missing Dependency",
                f"Required library not available: {e}\n"
                "Please install: pip install xarray netCDF4 pyproj numpy")
            return False
        except Exception as e:
            import traceback
            error_details = traceback.format_exc()
            self.log_message(f"ERROR creating latlon.nc file: {e}\n{error_details}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error creating latlon.nc file:\n{str(e)}")
            return False
    
    def write_all_layers(self):
        """
        Convert all masked layers to ASCII format using gdal:translate.
        Finds all masked .tif files in the geometry folder and converts them to .asc files
        with standardized names in the input/morph folder.
        """
        self.log_message("\n--- Converting all masked layers to ASCII format ---")
        
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        
        if not os.path.exists(geometry_folder):
            QMessageBox.warning(
                self.dialog, "Error",
                "Geometry folder does not exist. Please process layers first.")
            return
        
        # Output folder - save to input/morph folder
        input_folder = os.path.join(self.dialog.project_folder, "input")
        morph_folder = os.path.join(input_folder, "morph")
        os.makedirs(morph_folder, exist_ok=True)
        
        # Mapping from masked layer filenames to output ASCII filenames
        filename_mapping = {
            "1_dem_filled_masked.tif": "dem.asc",
            "1_dem_slope_masked.tif": "slope.asc",
            "1_dem_aspect_masked.tif": "aspect.asc",
            "2_flow_accumulation_masked.tif": "facc.asc",
            "2_flow_direction_masked.tif": "fdir.asc",
            "2_gauge_position_masked.tif": "idgauges.asc",
            "3_land_use_masked.tif": "luse.asc",
            "3_soil_masked.tif": "soil_class.asc",
            "3_geology_processed_masked.tif": "geology_class.asc"
        }
        
        # Find all masked .tif files in the geometry folder
        masked_layers = []
        for filename in os.listdir(geometry_folder):
            if filename.endswith("_masked.tif"):
                input_path = os.path.join(geometry_folder, filename)
                
                # Get output filename from mapping, or use default if not found
                output_filename = filename_mapping.get(filename)
                if not output_filename:
                    # Default: remove "_masked" and change extension to .asc
                    base_name = filename.replace("_masked.tif", "")
                    output_filename = f"{base_name}.asc"
                    self.log_message(f"WARNING: No mapping found for {filename}, using default name: {output_filename}")
                
                output_path = os.path.join(morph_folder, output_filename)
                layer_name = filename.replace("_masked.tif", "").replace("_", " ").title()
                masked_layers.append({
                    'input': input_path,
                    'output': output_path,
                    'name': layer_name,
                    'output_filename': output_filename
                })
        
        if not masked_layers:
            QMessageBox.warning(
                self.dialog, "No Layers",
                "No masked layers found. Please run 'Mask All Layers' first.")
            return
        
        self.log_message(f"Found {len(masked_layers)} masked layer(s) to convert...")
        
        # Convert each masked layer to ASCII
        success_count = 0
        failed_count = 0
        
        for layer_info in masked_layers:
            input_path = layer_info['input']
            output_path = layer_info['output']
            layer_name = layer_info['name']
            output_filename = layer_info['output_filename']
            
            # Skip if output already exists
            if os.path.exists(output_path):
                self.log_message(f"{layer_name} ASCII file ({output_filename}) already exists. Skipping...")
                success_count += 1
                continue
            
            self.log_message(f"Converting {layer_name} to {output_filename}...")
            
            # Use gdal:translate to convert to ASCII
            params_translate = {
                'INPUT': input_path,
                'TARGET_CRS': None,
                'NODATA': -9999,
                'COPY_SUBDATASETS': False,
                'OPTIONS': None,
                'EXTRA': '',
                'DATA_TYPE': 0,  # Use input data type
                'OUTPUT': output_path
            }
            
            result = self.run_processing_algorithm("gdal:translate", params_translate)
            
            if result and os.path.exists(output_path):
                self.log_message(f"{layer_name} converted to {output_filename} successfully.")
                success_count += 1
            else:
                self.log_message(f"ERROR: Failed to convert {layer_name} to {output_filename}.")
                failed_count += 1
        
        # Summary
        self.log_message(f"\n--- Conversion Summary ---")
        self.log_message(f"Successfully converted: {success_count} layer(s)")
        if failed_count > 0:
            self.log_message(f"Failed to convert: {failed_count} layer(s)")
        
        if success_count > 0:
            self.log_message(f"ASCII files saved to: {morph_folder}")
            self.log_message("ASCII conversion process completed.")
        else:
            QMessageBox.warning(
                self.dialog, "Conversion Error",
                "No layers were successfully converted to ASCII format.")
    
    def geology_classification_writer(self):
        """
        Read geological lookup file and produce geology_classdefinition.txt.
        
        Input CSV format:
        Geo-Class,Karstic,Default_Parameter_Value
        1,0,100
        2,0,100
        ...
        
        Output format:
        nGeo_Formations  <count>
        GeoParam(i)   ClassUnit     Karstic      Description
             1	          1           0      GeoUnit-1
             2     	      2           0      GeoUnit-2
        ...
        !<-END
        !***********************************
        ! NOTES
        !***********************************
        1 = Karstic
        0 = Non-karstic
        IMPORTANT ::
           Ordering has to be according to the ordering in mhm_parameter.nml
           (namelist: geoparameter)
        """
        self.log_message("\n--- Writing Geology Classification Definition File ---")
        
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        # Get geology lookup file path
        lookup_file = self.dialog.lineEdit_geology_lookup.text()
        
        if not lookup_file or not os.path.exists(lookup_file):
            self.log_message(
                "ERROR: Geology lookup table file not found or not specified.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Geology lookup table file.")
            return
        
        # Output file path - save to input/morph folder
        input_folder = os.path.join(self.dialog.project_folder, "input")
        morph_folder = os.path.join(input_folder, "morph")
        os.makedirs(morph_folder, exist_ok=True)
        output_file = os.path.join(morph_folder, "geology_classdefinition.txt")
        
        try:
            import pandas as pd
            
            # Try different encodings
            encodings = ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252']
            df = None
            used_encoding = None
            
            for encoding in encodings:
                try:
                    df = pd.read_csv(lookup_file, encoding=encoding)
                    used_encoding = encoding
                    self.log_message(f"Successfully read file with encoding: {encoding}")
                    break
                except (UnicodeDecodeError, pd.errors.EmptyDataError) as e:
                    continue
            
            if df is None or df.empty:
                self.log_message("ERROR: Could not read file with any encoding or file is empty.")
                QMessageBox.critical(
                    self.dialog, "Error",
                    "Could not read geology lookup table file. Encoding issue or file is empty.")
                return None
            
            # Log column names found
            self.log_message(f"CSV columns found: {', '.join(df.columns.tolist())}")
            
            # Check if required columns exist (case-insensitive)
            column_names_lower = {col.lower(): col for col in df.columns}
            geo_class_col = column_names_lower.get('geo-class') or column_names_lower.get('geo_class')
            karstic_col = column_names_lower.get('karstic')
            default_param_col = column_names_lower.get('default_parameter_value') or column_names_lower.get('default parameter value')
            
            if not geo_class_col:
                self.log_message("ERROR: 'Geo-Class' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'Geo-Class' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            if not karstic_col:
                self.log_message("ERROR: 'Karstic' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'Karstic' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            if not default_param_col:
                self.log_message("ERROR: 'Default_Parameter_Value' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'Default_Parameter_Value' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            self.log_message(f"Using columns: '{geo_class_col}', '{karstic_col}', '{default_param_col}'")
            
            # Clean and sort data
            df_clean = df[[geo_class_col, karstic_col, default_param_col]].dropna()
            
            # Sort by Geo-Class to ensure proper ordering
            df_clean = df_clean.sort_values(by=geo_class_col)
            
            # Get count of formations
            n_geo_formations = len(df_clean)
            
            self.log_message(f"Found {n_geo_formations} geological formations")
            
            # Write output file
            with open(output_file, 'w', encoding='utf-8') as f:
                # Write header
                f.write(f"nGeo_Formations  {n_geo_formations}\n")
                f.write("\n")
                f.write("GeoParam(i)   ClassUnit     Karstic      Description\n")
                
                # Write data rows
                for idx, row in df_clean.iterrows():
                    geo_class = int(float(row[geo_class_col])) if pd.notna(row[geo_class_col]) else 0
                    karstic = int(float(row[karstic_col])) if pd.notna(row[karstic_col]) else 0
                    default_param = row[default_param_col] if pd.notna(row[default_param_col]) else 100
                    description = f"GeoUnit-{geo_class}"
                    
                    # Format: GeoParam(i)   ClassUnit     Karstic      Description
                    # Match the exact format from the example
                    #         1	          1           0      GeoUnit-1
                    f.write(f"{geo_class:10d}\t{geo_class:10d}     {karstic:10d}      {description}\n")
                
                # Write footer
                f.write("!<-END\n")
                f.write("!***********************************\n")
                f.write("! NOTES\n")
                f.write("!***********************************\n")
                f.write("1 = Karstic\n")
                f.write("0 = Non-karstic\n")
                f.write("IMPORTANT ::\n")
                f.write("   Ordering has to be according to the ordering in mhm_parameter.nml\n")
                f.write("   (namelist: geoparameter)\n")
            
            self.log_message(f"Geology classification definition file written successfully: {output_file}")
            QMessageBox.information(
                self.dialog, "Success",
                f"Geology classification definition file created successfully.\n{output_file}")
            return output_file
            
        except ImportError:
            self.log_message("ERROR: pandas library is not installed.")
            QMessageBox.critical(
                self.dialog, "Error",
                "pandas library is required but not installed.\nPlease install it using: pip install pandas")
            return None
        except Exception as e:
            import traceback
            error_details = traceback.format_exc()
            self.log_message(
                f"ERROR writing geology classification file: {e}\n{error_details}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error writing geology classification file:\n{str(e)}")
            return None
    
    def soil_classdefinition_writer(self):
        """
        Read soil lookup file and produce soil_classdefinition.txt.
        
        Input CSV format:
        CLASS,Dominant_S,SNAM,NLAYERS,...,SOL_Z1,SOL_BD1,CLAY1,SAND1,...,SOL_Z2,SOL_BD2,CLAY2,SAND2,...
        
        Output format:
        nSoil_Types  <total_rows>
        MU_GLOBAL	HORIZON	UD[mm]	LD[mm]	CLAY[%]	SAND[%]	BD[gcm-3]
        1	1	0	300	50.0	25.0	1.75
        1	2	300	1000	50.0	25.0	1.75
        2	1	0	300	19.85	53.45	1.41
        2	2	300	1000	31.33	45.22	1.5
        ...
        """
        self.log_message("\n--- Writing Soil Classification Definition File ---")
        
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        # Get soil lookup file path
        lookup_file = self.dialog.lineEdit_soil_lookup.text()
        
        if not lookup_file or not os.path.exists(lookup_file):
            self.log_message(
                "ERROR: Soil lookup table file not found or not specified.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Soil lookup table file.")
            return
        
        # Output file path - save to input/morph folder
        input_folder = os.path.join(self.dialog.project_folder, "input")
        morph_folder = os.path.join(input_folder, "morph")
        os.makedirs(morph_folder, exist_ok=True)
        output_file = os.path.join(morph_folder, "soil_classdefinition.txt")
        
        try:
            import pandas as pd
            
            # Try different encodings
            encodings = ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252']
            df = None
            used_encoding = None
            
            for encoding in encodings:
                try:
                    df = pd.read_csv(lookup_file, encoding=encoding)
                    used_encoding = encoding
                    self.log_message(f"Successfully read file with encoding: {encoding}")
                    break
                except (UnicodeDecodeError, pd.errors.EmptyDataError) as e:
                    continue
            
            if df is None or df.empty:
                self.log_message("ERROR: Could not read file with any encoding or file is empty.")
                QMessageBox.critical(
                    self.dialog, "Error",
                    "Could not read soil lookup table file. Encoding issue or file is empty.")
                return None
            
            # Log column names found
            self.log_message(f"CSV columns found: {', '.join(df.columns.tolist())}")
            
            # Check if required columns exist (case-insensitive)
            column_names_lower = {col.lower(): col for col in df.columns}
            class_col = column_names_lower.get('class')
            nlayers_col = column_names_lower.get('nlayers')
            
            if not class_col:
                self.log_message("ERROR: 'CLASS' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'CLASS' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            if not nlayers_col:
                self.log_message("ERROR: 'NLAYERS' column not found in CSV file.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    f"CSV file must contain a 'NLAYERS' column.\nFound columns: {', '.join(df.columns.tolist())}")
                return None
            
            # Check for horizon columns (SOL_Z1, CLAY1, SAND1, SOL_BD1, etc.)
            horizon_columns = {}
            for i in range(1, 11):  # Check horizons 1-10
                sol_z_col = column_names_lower.get(f'sol_z{i}')
                clay_col = column_names_lower.get(f'clay{i}')
                sand_col = column_names_lower.get(f'sand{i}')
                sol_bd_col = column_names_lower.get(f'sol_bd{i}')
                
                if sol_z_col and clay_col and sand_col and sol_bd_col:
                    horizon_columns[i] = {
                        'sol_z': sol_z_col,
                        'clay': clay_col,
                        'sand': sand_col,
                        'sol_bd': sol_bd_col
                    }
            
            if not horizon_columns:
                self.log_message("ERROR: No horizon columns found (SOL_Z1, CLAY1, SAND1, SOL_BD1, etc.).")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    "CSV file must contain horizon columns (SOL_Z1, CLAY1, SAND1, SOL_BD1, etc.).")
                return None
            
            self.log_message(f"Found horizon columns for horizons: {sorted(horizon_columns.keys())}")
            
            # Clean data
            df_clean = df[[class_col, nlayers_col] + 
                         [col for h in horizon_columns.values() for col in h.values()]].dropna(subset=[class_col, nlayers_col])
            
            # Sort by CLASS to ensure proper ordering
            df_clean = df_clean.sort_values(by=class_col)
            
            # Build output rows
            output_rows = []
            
            for idx, row in df_clean.iterrows():
                mu_global = int(float(row[class_col])) if pd.notna(row[class_col]) else 0
                nlayers = int(float(row[nlayers_col])) if pd.notna(row[nlayers_col]) else 1
                
                # Limit nlayers to available horizons
                nlayers = min(nlayers, max(horizon_columns.keys()))
                
                prev_depth = 0  # Upper depth starts at 0 for first horizon
                
                for horizon in range(1, nlayers + 1):
                    if horizon not in horizon_columns:
                        self.log_message(f"WARNING: Horizon {horizon} data not available for CLASS {mu_global}. Skipping.")
                        continue
                    
                    cols = horizon_columns[horizon]
                    
                    # Get values
                    sol_z = float(row[cols['sol_z']]) if pd.notna(row[cols['sol_z']]) else 0
                    clay = float(row[cols['clay']]) if pd.notna(row[cols['clay']]) else 0
                    sand = float(row[cols['sand']]) if pd.notna(row[cols['sand']]) else 0
                    sol_bd = float(row[cols['sol_bd']]) if pd.notna(row[cols['sol_bd']]) else 0
                    
                    # Upper depth is previous horizon's lower depth (or 0 for first horizon)
                    ud = prev_depth
                    # Lower depth is SOL_Z for this horizon
                    ld = sol_z
                    
                    output_rows.append({
                        'mu_global': mu_global,
                        'horizon': horizon,
                        'ud': ud,
                        'ld': ld,
                        'clay': clay,
                        'sand': sand,
                        'bd': sol_bd
                    })
                    
                    # Update previous depth for next horizon
                    prev_depth = ld
            
            if not output_rows:
                self.log_message("ERROR: No valid soil horizon data found.")
                QMessageBox.warning(
                    self.dialog, "Input Error",
                    "No valid soil horizon data found in the lookup table.")
                return None
            
            n_soil_types = len(output_rows)
            self.log_message(f"Found {n_soil_types} soil horizon entries")
            
            # Write output file
            with open(output_file, 'w', encoding='utf-8') as f:
                # Write header
                f.write(f"nSoil_Types  {n_soil_types}\n")
                f.write("\n")
                f.write("MU_GLOBAL\tHORIZON\tUD[mm]\tLD[mm]\tCLAY[%]\tSAND[%]\tBD[gcm-3]\n")
                
                # Write data rows
                for row_data in output_rows:
                    f.write(f"{row_data['mu_global']}\t"
                           f"{row_data['horizon']}\t"
                           f"{row_data['ud']:.0f}\t"
                           f"{row_data['ld']:.0f}\t"
                           f"{row_data['clay']:.2f}\t"
                           f"{row_data['sand']:.2f}\t"
                           f"{row_data['bd']:.2f}\n")
            
            self.log_message(f"Soil classification definition file written successfully: {output_file}")
            QMessageBox.information(
                self.dialog, "Success",
                f"Soil classification definition file created successfully.\n{output_file}")
            return output_file
            
        except ImportError:
            self.log_message("ERROR: pandas library is not installed.")
            QMessageBox.critical(
                self.dialog, "Error",
                "pandas library is required but not installed.\nPlease install it using: pip install pandas")
            return None
        except Exception as e:
            import traceback
            error_details = traceback.format_exc()
            self.log_message(
                f"ERROR writing soil classification file: {e}\n{error_details}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error writing soil classification file:\n{str(e)}")
            return None
    
    def resetGeometry(self):
        """
        Reset geometry processing by:
        1. Removing all layers from QGIS interface that were created by the plugin
        2. Deleting all files in the geometry folder
        3. Resetting all path attributes
        """
        from qgis.PyQt.QtWidgets import QMessageBox
        
        # Confirm with user
        reply = QMessageBox.question(
            self.dialog, "Reset Geometry",
            "This will remove all geometry layers from QGIS and delete all files in the Geometry folder.\n\n"
            "Are you sure you want to continue?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No)
        
        if reply != QMessageBox.Yes:
            self.log_message("Reset cancelled by user.")
            return
        
        self.log_message("\n--- Resetting Geometry Processing ---")
        
        # Step 1: Remove layers from QGIS interface
        project = QgsProject.instance()
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        
        if not os.path.exists(geometry_folder):
            self.log_message("Geometry folder does not exist. Nothing to reset.")
            return
        
        # Get all layers and identify those from geometry folder
        layers_to_remove = []
        for layer_id, layer in project.mapLayers().items():
            source = layer.source()
            # Check if layer source is in geometry folder
            if geometry_folder in source.replace('\\', '/'):
                layers_to_remove.append(layer_id)
                self.log_message(f"Removing layer from QGIS: {layer.name()}")
        
        # Remove layers
        for layer_id in layers_to_remove:
            project.removeMapLayer(layer_id)
        
        if layers_to_remove:
            self.log_message(f"Removed {len(layers_to_remove)} layer(s) from QGIS interface.")
        else:
            self.log_message("No layers found in QGIS interface to remove.")
        
        # Step 2: Delete all files in geometry folder
        deleted_count = 0
        try:
            for root, dirs, files in os.walk(geometry_folder):
                for file in files:
                    file_path = os.path.join(root, file)
                    try:
                        os.remove(file_path)
                        deleted_count += 1
                        self.log_message(f"Deleted: {os.path.basename(file_path)}")
                    except Exception as e:
                        self.log_message(f"WARNING: Could not delete {file_path}: {e}")
            
            # Also try to remove empty subdirectories (optional, but clean)
            for root, dirs, files in os.walk(geometry_folder, topdown=False):
                for dir_name in dirs:
                    dir_path = os.path.join(root, dir_name)
                    try:
                        os.rmdir(dir_path)
                        self.log_message(f"Removed empty directory: {dir_name}")
                    except:
                        pass  # Directory not empty or error, skip
        
        except Exception as e:
            self.log_message(f"ERROR while deleting files: {e}")
            QMessageBox.warning(
                self.dialog, "Reset Error",
                f"Some files could not be deleted. Check the log for details.\nError: {e}")
        
        if deleted_count > 0:
            self.log_message(f"Deleted {deleted_count} file(s) from Geometry folder.")
        else:
            self.log_message("No files found to delete in Geometry folder.")
        
        # Step 3: Reset all path attributes
        self.filled_dem_path = None
        self.flow_dir_path = None
        self.channel_network_vector_path = None
        self.snapped_points_path = None
        self.watershed_raster_path = None
        self.watershed_vector_path = None
        self.aspect_path = None
        self.slope_path = None
        self.flow_direction_path = None
        self.flow_accumulation_path = None
        self.flow_accumulation_area_path = None
        self.gauge_position_path = None
        self.geology_path = None
        self.land_use_layer = None
        self.dem_layer = None
        
        self.log_message("All geometry processing attributes reset.")
        self.log_message("Geometry reset completed successfully.")

    def execute_all_processing(self):
        """
        Execute all processing steps sequentially in the following order:
        1. filled dem
        2. slope
        3. aspect
        4. land cover
        5. soil
        6. geology
        7. flow_accumulation
        8. flow_direction
        9. idgauges (will be processed after snap points due to dependency)
        10. channel network
        11. snap
        12. upslope area
        13. mask all layers
        14. process lat/lon headers
        15. write all layers (convert to ASCII)
        """
        self.log_message("\n=== Starting Execute All Processing ===")
        
        # Set flag to skip loading layers during batch processing
        self.skip_loading = True
        
        # Check prerequisites first
        if not self.check_prerequisites():
            self.log_message("ERROR: Prerequisites check failed. Aborting Execute All.")
            self.skip_loading = False  # Reset flag
            return
        
        try:
            # Step 1: Fill DEM
            self.log_message("\n--- Step 1/15: Fill DEM ---")
            self.fill_dem()
            if not self.filled_dem_path or not os.path.exists(self.filled_dem_path):
                self.log_message("ERROR: Fill DEM failed. Aborting Execute All.")
                return
            
            # Step 2: Slope
            self.log_message("\n--- Step 2/15: Process Slope ---")
            self.process_slope()
            
            # Step 3: Aspect
            self.log_message("\n--- Step 3/15: Process Aspect ---")
            self.process_aspect()
            
            # Step 4: Land Cover
            self.log_message("\n--- Step 4/15: Process Land Cover ---")
            self.process_land_use()
            
            # Step 5: Soil
            self.log_message("\n--- Step 5/15: Process Soil ---")
            self.process_soil()
            
            # Step 6: Geology
            self.log_message("\n--- Step 6/15: Process Geology ---")
            self.process_geology()
            
            # Step 7: Flow Accumulation
            self.log_message("\n--- Step 7/15: Process Flow Accumulation ---")
            self.process_flow_accumulation()
            if not self.flow_accumulation_path or not os.path.exists(self.flow_accumulation_path):
                self.log_message("ERROR: Flow Accumulation failed. Aborting Execute All.")
                return
            
            # Step 8: Flow Direction
            self.log_message("\n--- Step 8/15: Process Flow Direction ---")
            self.process_flow_direction()
            if not self.flow_direction_path or not os.path.exists(self.flow_direction_path):
                self.log_message("ERROR: Flow Direction failed. Aborting Execute All.")
                return
            
            # Step 9: ID Gauges (Gauge Position) - Note: requires snap points, will be processed after step 11
            self.log_message("\n--- Step 9/15: Process ID Gauges (deferred until after snap points) ---")
            # This will be processed after snap points in step 11
            
            # Step 10: Channel Network
            self.log_message("\n--- Step 10/15: Process Channel Network ---")
            self.process_channel_network()
            if not self.channel_network_vector_path or not os.path.exists(self.channel_network_vector_path):
                self.log_message("ERROR: Channel Network failed. Aborting Execute All.")
                return
            
            # Step 11: Snap Points
            self.log_message("\n--- Step 11/15: Snap Points ---")
            if not self.check_prerequisites(needs_pour_points=True):
                self.log_message("WARNING: Pour points not available. Skipping Snap Points step.")
            else:
                self.snap_points()
                # Now process ID Gauges (Step 9) after snap points are created
                if self.snapped_points_path and os.path.exists(self.snapped_points_path):
                    self.log_message("\n--- Processing Step 9/15: ID Gauges (now that snap points are available) ---")
                    self.process_gauge_position()
            
            # Step 12: Upslope Area (Delineate Watershed)
            self.log_message("\n--- Step 12/15: Delineate Watershed (Upslope Area) ---")
            if not self.snapped_points_path or not os.path.exists(self.snapped_points_path):
                self.log_message("WARNING: Snapped points not available. Skipping Upslope Area step.")
            else:
                self.delineate_watershed()
            
            # Step 13: Mask All Layers
            self.log_message("\n--- Step 13/15: Mask All Layers ---")
            if not self.merged_watershed_path or not os.path.exists(self.merged_watershed_path):
                self.log_message("WARNING: Merged watershed not available. Skipping Mask All Layers step.")
            else:
                self.mask_all_layers()
            
            # Step 14: Process Lat/Lon Headers
            self.log_message("\n--- Step 14/15: Process Lat/Lon Headers ---")
            self.process_lat_lon()
            
            # Step 15: Write All Layers (Convert to ASCII)
            self.log_message("\n--- Step 15/15: Write All Layers (Convert to ASCII) ---")
            self.write_all_layers()
            
            self.log_message("\n=== Execute All Processing Completed Successfully ===")
            
        except Exception as e:
            self.log_message(f"\nERROR: Execute All Processing failed with exception: {str(e)}")
            import traceback
            self.log_message(f"Traceback: {traceback.format_exc()}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Execute All Processing failed:\n{str(e)}")
        finally:
            # Always reset the skip_loading flag
            self.skip_loading = False
            self.log_message("Note: Raster files have been prepared but not loaded. Use individual buttons to load them.")
