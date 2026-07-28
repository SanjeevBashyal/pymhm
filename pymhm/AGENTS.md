# AGENTS.md

Guidance for coding agents working on the `pymhm` QGIS plugin.

## Project Shape

`pymhm` is a QGIS 3 plugin that prepares mHM model inputs, namelists, and supporting morphology/meteorology data.

Key entry points:

- `__init__.py`: QGIS plugin factory.
- `pymhm.py`: plugin lifecycle, toolbar/menu action, reload-safe import handling.
- `pymhm_dialog.py`: dialog wiring, signal connections, UI state persistence.
- `ui_pymhm_dialog_base.py` and `pymhm_dialog_base.ui`: generated/base UI files. Avoid touching generated UI unless the task is explicitly UI-related.
- `Morphology/processor.py`: morphology processor aggregate.
- `Meteorology/processor.py`: meteorology processor aggregate.
- `Configuration/processor.py`: schema-driven mHM namelist UI, rendering, and mHM run command.

Important packages:

- `Morphology/`: DEM, watershed, soil, geology, LAI, crop/mask/write-all, latlon, observations.
- `Meteorology/`: ERA5-Land to mHM forcing preparation.
- `Configuration/`: namelist schemas/templates/state/rendering/version compatibility.
- `mhm_tools_adapter.py`: small QGIS-free wrappers around the externally
  installed `mhm_tools` package for categorical data and `latlon.nc`.
- `Morphology/latlon/ascii_morphology.py`: UI-free morphology ASCII alignment
  and writing.
- `Meteorology/ERA5Land/mhm/`: ERA5-Land forcing preparation.
- `nml-schemas/` and `nml-templates/`: versioned by mHM tool version, currently `v5.13` and `v6`.
- `project-template/`: versioned project folder skeleton.

## Current Conventions

- The selected mHM version comes from `comboBox_mHMversion`.
- v5.13 uses `mhm_parameter.nml`; v6 uses `mhm_parameters.nml`.
- v5.13 schemas/templates live under `nml-schemas/v5.13` and `nml-templates/v5.13`; v6 under `v6`.
- v5.13 multi-domain template values should be rendered as one-shot arrays where possible, for example `resolution_Hydrology = 12000, 24000`, not repeated `resolution_Hydrology(1)` / `resolution_Hydrology(2)` lines.
- `mhm_tools` should be called through Python functions/modules, not shelling
  out to the `mhm-tools` CLI, unless explicitly requested. The dependency floor
  is currently `mhm-tools>=0.2.2` in both `pyproject.toml` and
  `requirements.txt`.
- Import the installed `mhm_tools` package directly. This repository does not
  bundle a second `pymhm/mhm_tools` source tree.
- UI-facing logs should go through `self.log_message(...)`.
- Success popups are not wanted for routine generated classdefinition files; log silently unless an error needs user action.

## Morphology Rules

- Filled DEM is an internal prerequisite for many workflows. Requirement calls should prepare/reuse it without loading it into QGIS. The user-facing `pushButton_fillDem` action should load the filled DEM layer.
- Temporary/intermediate geometry outputs belong in `Z Temp/Geometry`.
- Final mHM static morphology ASCII/classdefinition outputs belong in `data/static/morph`.
- Crop/mask outputs in `Z Temp/Geometry` use `_crop` and `_masked` suffixes. Do not automatically add crop/mask outputs to QGIS; groupBox processing buttons are used to show them.
- `pushButton_writeAll` writes masked rasters to ASCII files using mHM-compatible headers.
- Check raster dimensions across levels: L0 must be an integer multiple of L1, and L1/L11 must be compatible with L2 as required by mHM.

## Categorical Morphology Inputs

- Plain input combo boxes are exposed through compatibility adapters under the
  former `mMapLayerComboBox_*` names. Loaded QGIS layers are always listed;
  folder searching optionally adds relative project paths.
- Land cover, soil, and geology use either `mHM_ready` or `Lookup Table` mode.
  Lookup mode stores an explicit CSV or delimited TXT lookup file, mapping
  field, and class field.
- `mhm_tools.pre.rasterize_map_data`, `format_soil_data`,
  `format_geology_data`, and `format_lc_data` perform all categorical raster
  computation. Do not restore local QGIS reclassification or guessed mappings.
- Vector lookup processing burns the selected class field. Its follow-up
  formatter must therefore map `class_field` to itself; do not map the original
  source field a second time.
- The filled DEM defines the exact target grid and CRS. Original input and DEM
  CRS arguments are missing-metadata fallbacks, and must not override
  conflicting embedded metadata.
- Lookup mode retains `3_land_use.tif`, `3_soil.tif`, and
  `3_geology_processed.tif` for the existing crop/mask/write workflow.
- Soil and geology classdefinitions are produced once by their formatter and
  moved to `data/static/morph`. Only geology `PARAMETER_VALUE` metadata remains
  plugin-specific in `Z Temp/Geometry/geology_class_metadata.json`.
- `mHM_ready` accepts local `.asc`, `.nc`, or `.tif` rasters, copies them to
  `data/static/morph`, and bypasses categorical crop/mask/write processing.
  Soil and geology require their existing classdefinition in that folder;
  geology also requires its plugin metadata in `Z Temp/Geometry`.
- Prefer direct local paths. Materialize provider sublayers and non-file-backed
  QGIS layers only at the plugin boundary.
- Do not add Shapely validity checks or `make_valid` calls before rasterizing.

## LAI

- Currently implemented LAI workflow is input type item 0: long-term mean monthly gridded data from NetCDF.
- For NetCDF LAI, disable lookup table and lookup field controls.
- Clip/resample LAI to the filled DEM L0 grid. The filled DEM may be projected while LAI may be lat/lon.
- LAI is a 3D monthly array: 12 months plus spatial dimensions.
- Crop/mask operations must apply to all months.
- LAI output NetCDF data must be double precision (`float64`).
- Write the L0 header for LAI to `data/lai/header.txt`.

## Configuration/Namelists

- Schema-driven pages are built from `Configuration/schema_loader.py`.
- Template rendering lives in `Configuration/namelist.py`.
- Version compatibility, especially v5.13 schema-to-template mapping, lives in `Configuration/version_compat.py`.
- Path defaults live in `Configuration/path_defaults.py`.
- Preserve full v5.13 parameter templates so blocks like `&PET0` remain present.
- Be careful with Fortran namelist keys containing indices or derived-type fields, for example `eval_Per%yStart`, `GeoParam(:,1)`, `GeoParam(1,:)`.
- If changing schemas, keep templates and compatibility mappings in sync.

## Development Workflow

- Use `rg` / `rg --files` for search.
- Use `apply_patch` for manual edits.
- The worktree may already contain user changes. Do not revert unrelated modifications.
- Many modules import QGIS/PyQt and cannot be fully imported in plain system Python.
- Useful low-risk checks:
  - `python -m py_compile <changed pure-python files>`
  - `python3 -m pytest -q` for the categorical adapter, input-selection,
    processing, and logging-capture tests under `tests/`.
  - Install the standalone shim before importing `MorphologyProcessor` in a
    plain-Python MRO smoke test.
  - Targeted render smoke tests for `Configuration/namelist.py` and `Configuration/version_compat.py` using a lightweight fake `pymhm` package import harness.
- The standalone shell Python may not have `PyYAML`; QGIS Python may. If schema loading fails only because `yaml` is missing in shell Python, note it rather than treating it as a plugin failure.
- When QGIS behavior matters, reason from code and keep changes scoped; do not claim full runtime verification unless tested inside QGIS.

## Import And Reload Pitfalls

- QGIS reload can import `pymhm.py` as top-level `pymhm`; keep reload-safe import handling intact.
- Avoid introducing package-relative imports into files that may be loaded top-level unless existing code already handles that context.
- Circular mixin inheritance can break plugin reload with MRO errors. Prefer linear aggregate mixins and small shared helper modules.

## Output Paths To Remember

- Temporary geometry: `Z Temp/Geometry`
- Static morphology: `data/static/morph`
- LAI: `data/lai`
- Meteorology: `data/meteo`
- Root namelists: project root (`mhm.nml`, parameter namelist, `mhm_outputs.nml`)
- mHM outputs: `output`
- Restarts: `restart`

## Style

- Keep UI code focused on UI state, validation, and logging.
- Keep the `mhm_tools_adapter.py` boundary QGIS-free; put generally reusable
  raster computation in the external `mhm-tools` project and expose it through
  a public Python function plus its CLI command where required.
- Do not copy external `mhm_tools` source back into this repository.
- Prefer clear errors over silent fallbacks when the user is expected to provide fields/layers.
- Keep changes narrow and aligned with existing mixin/module boundaries.
