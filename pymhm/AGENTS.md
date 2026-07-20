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
- `mhm_tools_to_integrate/`: path-oriented adapters and UI-free helpers. Soil
  and geology adapters under `setup_creation/` delegate computation to the
  externally installed `mhm_tools` package.
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
- Despite its legacy name,
  `mhm_tools_to_integrate._bundled.ensure_bundled_mhm_tools()` currently falls
  through to the installed package because this tree no longer contains
  `pymhm/mhm_tools`.
- UI-facing logs should go through `self.log_message(...)`.
- Success popups are not wanted for routine generated classdefinition files; log silently unless an error needs user action.

## Morphology Rules

- Filled DEM is an internal prerequisite for many workflows. Requirement calls should prepare/reuse it without loading it into QGIS. The user-facing `pushButton_fillDem` action should load the filled DEM layer.
- Temporary/intermediate geometry outputs belong in `Z Temp/Geometry`.
- Final mHM static morphology ASCII/classdefinition outputs belong in `data/static/morph`.
- Crop/mask outputs in `Z Temp/Geometry` use `_crop` and `_masked` suffixes. Do not automatically add crop/mask outputs to QGIS; groupBox processing buttons are used to show them.
- `pushButton_writeAll` writes masked rasters to ASCII files using mHM-compatible headers.
- Check raster dimensions across levels: L0 must be an integer multiple of L1, and L1/L11 must be compatible with L2 as required by mHM.

## Soil And Geology

- Soil and geology lookup tables are selected from QGIS with `QgsMapLayerComboBox`, not file-only inputs.
- The lookup field combo boxes beside each lookup table must be populated from the selected lookup layer fields.
- The selected lookup field is mandatory. Do not use fallback field names such as `Dominant_S` or `GEO_CLASS`; if the required user-selected mapping cannot be applied, raise/log a clear error.
- Soil inputs may be vector or raster. The lookup table must include `SOIL_CLASS`; output raster values come from `SOIL_CLASS`.
- Geology inputs may be vector or raster. The lookup table must include `GEOLOGY_CLASS`; output raster values come from `GEOLOGY_CLASS`.
- The filled DEM is required for both vector and raster inputs. It defines the
  exact output CRS, affine transform, extent, width, and height.
- Vector inputs are handled by `mhm_tools.pre.rasterize_map_data` through
  `rasterize_soil_map` or `rasterize_geology_map`. Pass the selected source
  field as the vector mapping field, the selected lookup field as the lookup
  mapping field, and burn `SOIL_CLASS` or `GEOLOGY_CLASS` respectively.
- The external `data-converter rasterize-map` CLI requires `-i/--input-file`,
  `-d/--dem-file`, `-o/--output-file`, and `-b/--burn-field`. Lookup mode adds
  the paired `-l/--lookup-table` and `-m/--mapping-field` options; direct mode
  burns the numeric vector field named by `-b`.
- Do not add Shapely validity checks or `make_valid` calls before rasterizing.
  The current contract passes geometries directly to Rasterio so invalid or
  empty features follow Rasterio's normal handling.
- Raster inputs are handled by `mhm_tools.pre.format_soil_data` or
  `format_geology_data`. Input and DEM files may be `.asc`, `.tif`/`.tiff`, or
  `.nc`; categorical alignment must use nearest-neighbour resampling.
- The corresponding formatter CLIs require `-i`, `-d`, `-o`, `-l`, and `-m`;
  `-t/--output-type` selects `nc` or `asc`. Across all three commands,
  `-s/--input-crs` and `-r/--dem-crs` are metadata fallbacks.
- Forward the QGIS layer CRS as `input_crs` or `dem_crs` only as a fallback for
  files without embedded metadata. A supplied CRS must not silently override a
  conflicting embedded CRS. ASCII CRS metadata uses a neighboring `.prj` file.
- Prefer a direct local provider path. Materialize provider sublayers or
  non-file-backed vector layers to temporary GeoPackage files and raster layers
  to temporary GeoTIFF files.
- `format_soil_data` and `format_geology_data` intentionally expose only `nc`
  and `asc` output types. The shared `_categorical.format_categorical_file`
  adapter uses a `TemporaryDirectory`, requests temporary NetCDF output, then
  writes the plugin's final GeoTIFF. Keep that temporary directory until both
  the raster and classdefinition have been moved successfully.
- `soil_classdefinition.txt` and `geology_classdefinition.txt` are written to `data/static/morph`.
- Classdefinition rendering belongs to `mhm_tools`; capture its log records
  through `mhm_tools_to_integrate.logging.capture_messages` and emit one
  completion message per file, never one message per definition row.
- Geology metadata for configuration/parameters is plugin-specific and is
  written to `Z Temp/Geometry/geology_class_metadata.json`. Keep
  `GEOLOGY_CLASS` and `PARAMETER_VALUE` validation in `pymhm`; standalone
  `mhm_tools` geology classdefinition output does not use `PARAMETER_VALUE`.
- A failed geology refresh must not destroy a previously valid
  classdefinition or metadata file. Restore prior contents and remove any
  partial raster output.
- Existing `3_soil.tif` and `3_geology_processed.tif` files are currently
  reused based on file existence, without input/lookup provenance checks. When
  changing cache behavior, prevent an old raster from being paired with a
  classdefinition generated from a newer lookup table.
- Prefer GeoPackage (`.gpkg`) for temporary vector outputs; shapefile paths with spaces/backslashes have caused QGIS Processing creation errors.

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
  - `python3 -m pytest -q` for `tests/test_soil_mhm_tools_adapter.py`,
    `tests/test_geology_mhm_tools_adapter.py`, and logging-capture tests.
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
- Keep QGIS adapters in `mhm_tools_to_integrate`; put generally reusable raster
  computation in the external `mhm-tools` project and expose it through a
  public Python function plus its CLI command where required.
- Do not copy external `mhm_tools` source back into this repository.
- Prefer clear errors over silent fallbacks when the user is expected to provide fields/layers.
- Keep changes narrow and aligned with existing mixin/module boundaries.
