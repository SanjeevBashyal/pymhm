# -*- coding: utf-8 -*-
"""Constants for schema-driven namelist configuration."""

CONFIG_ORDER = [
    "config_project",
    "config_processes",
    "config_optimize",
    "config_time",
    "config_input",
    "config_coupling",
    "config_meteo",
    "config_mpr",
    "config_mhm",
    "config_mrm",
    "config_observations",
]

OUTPUT_ORDER = [
    "output_mhm",
    "output_mrm",
]

PROCESS_PARAMETER_MAP = {
    "interception": {1: "interception1"},
    "snow": {1: "snow1"},
    "soil_moisture": {
        1: "soilmoisture1",
        2: "soilmoisture2",
        3: "soilmoisture3",
        4: "soilmoisture4",
    },
    "direct_runoff": {1: "directrunoff1"},
    "pet": {
        -2: "petm2",
        -1: "petm1",
        1: "pet1",
        2: "pet2",
        3: "pet3",
    },
    "interflow": {1: "interflow1"},
    "percolation": {1: "percolation1"},
    "routing": {
        1: "routing1",
        2: "routing2",
        3: "routing3",
    },
    "neutrons": {
        1: "neutrons1",
        2: "neutrons2",
    },
    "temperature_routing": {1: "rivertemp1"},
}

PARAMETER_ALWAYS_BLOCKS = [
    "geoparameter",
]

OUTPUT_FILES = {
    "mhm": "mhm.nml",
    "parameters": "mhm_parameters.nml",
    "outputs": "mhm_outputs.nml",
}

OUTPUT_FILE_OVERRIDES = {
    "v5.13": {
        "parameters": "mhm_parameter.nml",
    },
}

TEMPLATE_NAMES = {
    "mhm": ["mhm-template.nml"],
    "parameters": [
        "mhm-parameters-template.nml",
        "mhm-parameter-template.nml",
    ],
    "outputs": ["mhm-outputs-template.nml"],
}

STATUS_LABELS = {
    "mhm": "label_mHMConfigurationStatus",
    "parameters": "label_parametersConfigurationStatus",
    "outputs": "label_outputConfigurationStatus",
}

STATUS_TITLES = {
    "mhm": "mHM",
    "parameters": "Parameters",
    "outputs": "Outputs",
}
