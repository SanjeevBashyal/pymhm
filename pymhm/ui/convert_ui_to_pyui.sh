pyuic5 -x pymhm/ui/pymhm_dialog_base.ui -o pymhm/pyui/ui_pymhm_dialog_base.py
pyrcc5 pymhm/ui/resources.qrc -o pymhm/resources_rc.py

pyuic5 -x pymhm/ui/mhm_ready_dialog.ui -o pymhm/pyui/ui_mhm_ready_dialog.py
pyuic5 -x pymhm/ui/single_layer_input_with_lookup_config_dialog.ui -o pymhm/pyui/ui_single_layer_input_with_lookup_config_dialog.py
pyuic5 -x pymhm/ui/land_use_historical_input.ui -o pymhm/pyui/ui_land_use_historical_input.py
pyuic5 -x pymhm/ui/soil_multi_horizon_input.ui -o pymhm/pyui/ui_soil_multi_horizon_input.py
pyuic5 -x pymhm/ui/lai_netcdf_input_dialog.ui -o pymhm/pyui/ui_lai_netcdf_input_dialog.py


pyuic5 -x pymhm/ui/domain_delineator_dialog.ui -o pymhm/pyui/ui_domain_delineator_dialog.py
pyuic5 -x pymhm/ui/domain_and_discharge_table_assignment_dialog.ui -o pymhm/pyui/ui_domain_and_discharge_table_assignment_dialog.py

pyuic5 -x pymhm/ui/discharge_table_assignment_dialog.ui -o pymhm/pyui/ui_discharge_table_assignment_dialog.py
pyuic5 -x pymhm/ui/elevation_band_dialog.ui -o pymhm/pyui/ui_elevation_band_dialog.py

pyuic5 -x pymhm/ui/project_terminal_dialog.ui -o pymhm/pyui/ui_project_terminal_dialog.py
