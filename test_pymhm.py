# E:\0 Python\pymhm\pymhm\test\test_dialog_standalone.py

import os
import sys

# --- IMPORTANT: CONFIGURATION SECTION ---
# You MUST provide paths to some real test data here for the dialog to work.
# Replace these placeholder paths with actual paths to a DEM and a points shapefile on your system.
TEST_DATA_DIR = r"E:\8 WRE Research\5 mHM\QGIS Plugin Test\Geometry"  # <--- EDIT THIS: A folder with your test data
DEM_PATH = r"C:\Users\Ripple\Downloads\mHM workshop\clip_utm_44_project.tif"       # <--- EDIT THIS: Path to a test DEM
POUR_POINTS_PATH = r"C:\Users\Ripple\Downloads\mHM workshop\pour_points.gpkg" # <--- EDIT THIS: Path to a test pour points shapefile

# This ensures that the script can find the 'pymhm' package and its modules
# It assumes this script is in pymhm/test/ and adds the parent directory (pymhm/) to the path.
plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(plugin_dir))
# ----------------------------------------

from qgis.core import QgsApplication, QgsProject, QgsVectorLayer, QgsRasterLayer
from qgis.PyQt.QtWidgets import QApplication

# Now that the path is set up, we can import the dialog from your plugin package
from pymhm.pymhm_dialog import pymhmDialog

def main():
    """
    Main function to initialize QGIS, load test data, and run the dialog.
    """
    print("--- Initializing QGIS Application for standalone testing ---")
    
    # 1. Initialize the QGIS Application. This is crucial for any QGIS-related code.
    # We use QgsApplication.instance() which might return an existing instance if run inside QGIS,
    # or we create a new one. The 'True' indicates it's a GUI-enabled application.
    qgs_app = QgsApplication.instance()
    if qgs_app is None:
        qgs_app = QgsApplication([], True)

    # Use the QGIS_PREFIX_PATH from your launch.json to find QGIS resources
    qgis_prefix = os.environ.get('QGIS_PREFIX_PATH')
    if not qgis_prefix:
        raise RuntimeError("QGIS_PREFIX_PATH environment variable not set. Ensure you are running this with the correct debugger configuration.")
        
    QgsApplication.setPrefixPath(qgis_prefix, True)
    qgs_app.initQgis()
    print(f"QGIS Prefix Path set to: {QgsApplication.prefixPath()}")
    print("QGIS Initialized.")

    # 2. Get a reference to the project instance. This is where layers are managed.
    project = QgsProject.instance()

    # 3. Load dummy/test layers into the project.
    # This will populate the QgsMapLayerComboBox widgets in your dialog.
    print(f"Loading test DEM: {DEM_PATH}")
    dem_layer = QgsRasterLayer(DEM_PATH, "Test DEM")
    if not dem_layer.isValid():
        print(f"ERROR: Failed to load test DEM layer. Please check the path: {DEM_PATH}")
    else:
        project.addMapLayer(dem_layer)
        print("Test DEM loaded successfully.")

    print(f"Loading test Pour Points: {POUR_POINTS_PATH}")
    points_layer = QgsVectorLayer(POUR_POINTS_PATH, "Test Pour Points", "ogr")
    if not points_layer.isValid():
        print(f"ERROR: Failed to load test pour points layer. Please check the path: {POUR_POINTS_PATH}")
    else:
        project.addMapLayer(points_layer)
        print("Test Pour Points loaded successfully.")

    print("\n--- Starting the pymhmDialog ---")
    # 4. Create an instance of your dialog and show it.
    dialog = pymhmDialog()
    dialog.show()

    # 5. Run the Qt event loop. This makes the dialog interactive.
    # The script will pause here until you close the dialog.
    qgs_app.exec_()

    print("\n--- Dialog closed, cleaning up. ---")
    # 6. Clean up the QGIS application environment.
    qgs_app.exitQgis()
    print("QGIS environment closed.")


if __name__ == "__main__":
    # Ensure a QApplication exists. It's needed for any widgets.
    # This is a bit of boilerplate for running PyQt apps standalone.
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    main()