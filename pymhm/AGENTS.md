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

## Common L0/L2 Extent

- The common extent is built once from the merged active-domain polygons in
  `grid_resolution.aligned_l0_l2_headers()`. L2 is snapped outward to the L0
  anchor with floor/ceil; L0 is then derived from L2 (`ncols*n`, `nrows*n`,
  same corners), so `n x n` L0 cells per L2 cell holds by construction.
  `validate_l0_l2_alignment()` guards every header pair.
- **Never round a grid-defining cell size.** `floor_cellsize()` / `ceil_cellsize()`
  are for display and for the resolution text saved in the UI state, nothing
  else. Two separate failures came from rounding them:
  - Flooring a repeating L0 cell size such as 1/1200 degrees to 8 places drifts
    a fraction of a cell per row; over a few thousand rows the L0 window copy
    rejects every raster as misaligned.
  - `120 x` a 1/1200-degree cell size is `0.09999999999999999`, one ulp under
    `0.1`. Flooring *that* to 8 places gives `0.09999999`, a 1e-8 error that
    fails `validate_l0_l2_alignment` and makes `possible_resolutions()` return
    nothing.
  Grid construction uses `exact_resolution` / `header["exact_cellsize"]` from
  `raster_resolution_info()` and `current_l0_resolution()`. The L2 cell size is
  carried verbatim from `aligned_l0_l2_headers()` through
  `update_l2_resolution_from_metadata()`, `grid_level_headers()`,
  `read_header_file()`, `_standardize_header()`, `_resolution_ratio()`, and the
  `header_for_*` builders. `_target_l0_header()` prefers the saved exact L0
  header; never rebuild it from a rounded resolution.
- Never shift or shrink a misaligned boundary; expand it outward only.
- Prepared L0 rasters are never resampled. `file_tasks.crop_aligned_l0_raster()`
  and `file_tasks.mask_aligned_l0_raster()` copy integer windows, pad outside
  the source with nodata, and raise on a misaligned or wrong-CRS input. Do not
  reintroduce `gdal:warpreproject` or `gdal:cliprasterbymasklayer` here.
- Advanced historical land cover and multi-horizon soil are formatted on the
  filled-DEM grid and then placed on the common L0 header by
  `Morphology/layers/advanced_l0.py`, which uses
  `ascii_morphology.pad_l0_file_to_header()` (cell lookup plus nodata padding,
  never interpolation). This is Morphology Setup step 5/5.
- **ASCII grids are padded by streaming text**, in
  `Morphology/latlon/ascii_pad.py`. Loading one to pad it costs gigabytes:
  measured at 3238 MiB for a 13201 x 6001 soil class raster, which is why soil
  failed inside QGIS once crop and mask had already claimed memory. The
  streaming version is 39 MiB and 1 s for the same file. Only the soil *class*
  raster is padded; the horizon inputs (clay, sand, silt, bulk density) are
  never touched.
- Each Morphology Setup step records its own status in `workflows` as
  `<workflow>_crop`, `_mask`, `_latlon`, `_write`, `_publish`, so a failure
  points at the step that failed.
- L2 meteorology inputs already on the target grid are sliced and padded by
  `Meteorology/l2_grid.py`. Only misaligned or reprojected inputs are resampled,
  and then with nearest-neighbour only. Every written file and header is checked
  against the saved L2 header before publication.
- ERA5-Land temperature files are read once for `tavg`, `tmin`, and `tmax`
  through `iter_daily_datasets()`, but only when the estimated hold fits
  `PYMHM_METEO_MAX_BYTES` (2 GiB default). A long record needs all three series
  in memory at once, which triples the peak, so it falls back to one variable
  at a time. Consume the iterator and release each dataset after writing it;
  `build_daily_datasets()` collects them all and is for tests and short records.
- `Meteorology/reuse.py` skips meteorology preparation when every required
  variable is recorded in `pymhm_processing_state.json`, still on disk, and its
  `header.txt` matches the current L2 header. A changed L2 grid always rebuilds.
  `pet` is only required when a PET folder is selected.
- The validated L0 header, L2 header, and multiplier are mirrored into the
  `grid` section of `pymhm_processing_state.json` for resumes.
- Model-ready shared inputs live below `data/master`: static morphology in
  `data/master/static/morph`, forcing in `data/master/meteo`, LAI in
  `data/master/lai`, lat/lon at `data/master/latlon.nc`, and observations in
  `data/master/observation/streamflow`.
- A checked domain outlet creates a watershed-masked and bounding-box-cropped
  L0 DEM at `data/<outlet-id>/dem.asc`. Domain directories are created from
  saved outlet IDs, not from literal template placeholders.
- Keep the Domain Outlet label and checkbox enabled regardless of their current
  value. An unchanged outlet location is a valid domain/gauge location; a map
  pick overrides it.

## Staging And Publication

- Execute All writes land cover, soil, and geology outputs into
  `Z Temp/Morphology/morph` (`morph_staging_folder`), never straight into
  `data/master`. Nothing is model-ready until the common L0 extent is known, and
  that only happens once meteorology has fixed the L2 grid.
- Morphology Setup step 5/5 (`align_advanced_inputs_to_l0` ->
  `advanced_l0.publish_model_inputs`) pads each staged raster onto the L0 header,
  moves the whole staged set into `data/master/static/morph` (class definitions
  and `.prj` sidecars included), and rewrites the `output_path` /
  `classdefinition_path` / `scenes[].output_path` entries in `nml-settings.json`
  to the published locations. Those values feed `vSpecific/v5_13` and
  `vSpecific/v6`, so they must name files that exist.
- Repointing also fixes paths left pointing at staging when the file is already
  in master, so a reused stage still ends up resolvable.
- The step finishes by calling `missing_model_inputs()` and fails when anything
  the handoff names is absent. That is the guarantee that a completed
  Morphology Setup leaves `nml-settings.json` usable by nml-tools.
- `mHM_ready` is the exception: it bypasses crop/mask/write, so its raster and
  class definition are already final and are written directly to `data/master`.
- Execute All verifies its categorical outputs in **staging**; only mHM-ready is
  verified in `data/master`.

## Fingerprinted Reuse

- `state_cache.py` records, in `pymhm_processing_state.json`, a fingerprint of
  each stage's inputs next to the outputs it produced. A stage reruns only when
  an input's size or modification time changed, its configuration changed, or one
  of its recorded outputs is gone. Content is never hashed, so the check stays
  cheap on multi-gigabyte inputs.
- Sections in use: `stages` (per Execute All stage) and `meteo_inspection`
  (per meteorology folder). Keep entries additive; `store_payload` merges into
  the existing state rather than replacing it.
- The categorical lookup stages and the LAI stage are gated on it. Geology used
  to rerun every time because only land cover and soil passed `reuse_existing`,
  and that check was existence-only rather than input-based. Land cover and soil
  keep `reuse_existing=True` because their advanced mode has no lookup job to
  fingerprint; removing it would make advanced runs rebuild every time.
- Outputs already on disk with no recorded fingerprint are **adopted** when they
  postdate every input (`outputs_newer_than_inputs`), so a project prepared
  before the cache existed avoids one pointless rebuild. An output older than an
  input is certainly stale and is rebuilt instead.
- `ProcessingStateMixin.save_processing_state` merges: it rewrites only
  `OWNED_SECTIONS` and preserves the rest. Dumping its in-memory copy wholesale
  erased the fingerprints another writer had added since it loaded, which
  silently disabled all reuse.
- `Meteorology/inspection_cache.py` wraps `inspect_meteo_folder`, which opens
  every NetCDF in a folder **twice**. On a 552-file ERA5-Land record that froze
  project load for ~22 s per folder; reusing an unchanged inspection brings the
  whole `load_input_state` from 22.6 s to 0.26 s. Do not call
  `inspect_meteo_folder` directly from the dialog.

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
- **Execute All has two implementations. The button uses the bridge.**
  `pushButton_executeAllMorphology` goes to
  `MorphologyTaskBridge.start_execute_all()`, a `QgsTask` stage list
  (fill, terrain, land cover, soil, geology, LAI, hydrology).
  `ExecuteAllMixin.execute_all_processing()` is the synchronous twin and is not
  reachable from the UI. A stage added only to the mixin silently never runs --
  that is exactly how the LAI stage was missed. Add stages to both.
- LAI follows the same two-stage path as every other layer:
  - **Execute All step 7/15** (`resample_lai_to_dem_grid`, bridge stage
    `start_lai`) resamples the source
    onto the filled DEM grid with **bilinear** interpolation and stages it at
    `Z Temp/Morphology/lai/lai_dem.nc`. The result has exactly the DEM's size
    and extent.
  - **Morphology Setup** then reaches the common L0 extent by integer **window
    copy** with nodata padding (`crop_lai_netcdf_to_l0`), and masks the same way
    (`mask_lai_netcdf_to_l0`). Neither stage resamples.
  Keep the resampling in Execute All: changing the L2 multiplier must not force
  a re-resample, only a re-crop.
- `Morphology/layers/lai_source.py` is the QGIS-free LAI pipeline over paths and
  primitives, so it runs inside a `QgsTask`. `lai.py` is the thin QGIS boundary:
  `lai_task_options()` snapshots the primitives on the main thread and
  `run_lai_resample()` is the worker entry point. Never call a processor method
  or touch a widget from the worker.
- `Morphology/layers/lai_l0.py` streams every stage through `stream_lai_grid()`,
  one block of target rows and one time step at a time. Never materialise the
  whole cube or the full WGS84 target mesh — 468 monthly steps on a
  13320 x 6120 grid is 284 GiB, and the two 2-D coordinate arrays alone are
  652 MiB each. Measured on that grid: ~370 MiB peak, about 70:1 compression.
- The stages differ only in their sampler: `ResampleSampler` (bilinear or
  nearest, by lon/lat lookup) and `WindowSampler` (integer window, nodata
  padding). `lai_window_offsets()` validates the cell size and alignment first
  and raises rather than silently shifting the data.
- Bilinear renormalises over the finite corners, so one missing source cell
  does not poison its neighbours.
- `_lazy_target_grid()` returns a `row_lonlat(start, stop)` generator instead of
  a mesh, with an identity fast path for a geographic project CRS.
  `separable_axes()` then collapses a geographic block to two 1-D index arrays
  so the gather uses `np.ix_`; keep both paths, projected grids need the 2-D one.
- The LAI variable is written at zlib level 1 on purpose: level 4 compresses an
  upsampled grid no better (78:1 vs 76:1) and takes nearly twice as long.
- Bilinear output barely compresses. Measured on real GIMMS LAI upsampled 100x:
  nearest-neighbour reached 70:1 on its runs of repeated values, bilinear only
  1.80:1 on the DEM grid and 1.85:1 on the L0 grid because every cell differs.
  That is why `LAI_ASSUMED_COMPRESSION` is 1.8,
  and why a long record at L0 resolution is expensive on disk: 468 monthly steps
  on a 13320 x 6120 grid is ~158 GiB per file, and the pipeline writes three
  (`lai_dem.nc`, `lai_crop.nc`, `lai.nc`). Check the numbers before assuming a
  long record is feasible; the 12-step long-term-mean option is ~4 GiB per file.
- `lai.nc` is the published masked file; `lai_masked.nc` is not written, because
  a second full-size copy costs disk for no benefit.
- `assert_lai_output_fits()` is a disk check, not a memory check: memory is
  bounded by the writer. `PYMHM_LAI_MAX_BYTES` overrides it.
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
- `pymhm_domain_delineation_state.json` remains canonical. Project domain
  directories, `dem.asc` paths, inverse domain-to-gauge membership, and gauge
  metadata into `nml-settings.json` before opening nml-tools.
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
- Shared static morphology: `data/master/static/morph`
- Domain DEM: `data/<outlet-id>/dem.asc`
- LAI: `data/master/lai`
- Meteorology: `data/master/meteo`
- Streamflow: `data/master/observation/streamflow`
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

## Background Processing

- Widgets, `QgsProject`, `QgsMapLayer`, and map-canvas operations are main-thread
  only. Snapshot paths, CRS text, numbers, and plain containers before starting
  work.
- Use `TaskCoordinator`/`QgsTask` for path-only computation and return paths or
  serializable results to a main-thread completion callback.
- **All heavy categorical formatting runs in `pymhm.native_worker`**, a child
  process, dispatched through `_run_in_worker`: advanced land cover and soil, and
  the lookup mode used by geology. The lookup path used to run inside QGIS and
  OOM-killed it -- `anon-rss:5756792kB`, i.e. 5.5 GiB, formatting a vector input
  onto a 13201 x 6001 grid. The worker sets `oom_score_adj=500` so the kernel
  takes the child first. Cancellation must terminate that child; a nonzero or
  signal exit is reported as a task failure rather than crashing QGIS.
- `categorical_lookup.py` and `geology_metadata.py` sit at the top level, not
  under `pymhm/Morphology/`, because importing anything from that package runs
  `Morphology/__init__` -> `common` -> `import processing`, which does not exist
  outside QGIS. Keep worker-side modules clear of that package.
- Keep large categorical raster algorithms bounded by raster windows and write
  one historical period before allocating the next.
