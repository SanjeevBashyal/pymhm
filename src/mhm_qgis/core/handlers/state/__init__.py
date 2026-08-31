# -*- coding: utf-8 -*-
"""The project's own persisted state, one module per store.

Every JSON file the plugin keeps inside a user project is read and written from
here. Qt-free by contract, so the state can be exercised and migrated without a
dialog: the modules take a project folder, not a dialog.

- `input_state`      mhm_qgis_input_state.json
- `processing_state` mhm_qgis_processing_state.json (three writers -- see below)
- `domain_state`     mhm_qgis_domain_delineation_state.json
- `nml_settings`     nml-settings.json, the handoff to nml-tools
- `cache`            fingerprinted stage reuse, inside the processing state
- `metadata`         meteo grid and geology class metadata

`mhm_qgis_processing_state.json` has three independent writers. Each owns named
sections and must merge rather than overwrite: dumping a whole in-memory copy
erases what the others added since it was loaded, which silently disables stage
reuse. That has broken once already.
"""
