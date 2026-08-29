pyuic5 -x src/mhm_qgis/ui/mhm_qgis_dialog_base.ui -o src/mhm_qgis/pyui/ui_mhm_qgis_dialog_base.py
pyrcc5 src/mhm_qgis/ui/resources.qrc -o src/mhm_qgis/resources_rc.py

pyuic5 -x src/mhm_qgis/ui/mhm_ready_dialog.ui -o src/mhm_qgis/pyui/ui_mhm_ready_dialog.py
pyuic5 -x src/mhm_qgis/ui/single_layer_input_with_lookup_config_dialog.ui -o src/mhm_qgis/pyui/ui_single_layer_input_with_lookup_config_dialog.py
pyuic5 -x src/mhm_qgis/ui/land_use_historical_input.ui -o src/mhm_qgis/pyui/ui_land_use_historical_input.py
pyuic5 -x src/mhm_qgis/ui/soil_multi_horizon_input.ui -o src/mhm_qgis/pyui/ui_soil_multi_horizon_input.py
pyuic5 -x src/mhm_qgis/ui/lai_netcdf_input_dialog.ui -o src/mhm_qgis/pyui/ui_lai_netcdf_input_dialog.py


pyuic5 -x src/mhm_qgis/ui/domain_delineator_dialog.ui -o src/mhm_qgis/pyui/ui_domain_delineator_dialog.py
pyuic5 -x src/mhm_qgis/ui/domain_and_discharge_table_assignment_dialog.ui -o src/mhm_qgis/pyui/ui_domain_and_discharge_table_assignment_dialog.py

pyuic5 -x src/mhm_qgis/ui/discharge_table_assignment_dialog.ui -o src/mhm_qgis/pyui/ui_discharge_table_assignment_dialog.py
pyuic5 -x src/mhm_qgis/ui/elevation_band_dialog.ui -o src/mhm_qgis/pyui/ui_elevation_band_dialog.py

pyuic5 -x src/mhm_qgis/ui/project_terminal_dialog.ui -o src/mhm_qgis/pyui/ui_project_terminal_dialog.py
pyuic5 -x src/mhm_qgis/ui/thread_display_dialog.ui -o src/mhm_qgis/pyui/ui_thread_display_dialog.py
