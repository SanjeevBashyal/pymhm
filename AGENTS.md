# AGENTS.md

Guidance for coding agents working at the repository root of `mhm_qgis`.

## Repository Shape

This repository publishes one Python distribution and also contains the QGIS
plugin folder:

- `src/mhm_qgis/`: QGIS plugin folder and Python package.
- `src/mhm_qgis/metadata.txt`, `src/mhm_qgis/__init__.py`, `src/mhm_qgis/mhm_qgis.py`: QGIS plugin
  entry points.
- `src/mhm_qgis/standalone/cli.py`: PyPI console entry point. The command is
  `mhm-qgis` and it opens the standalone GUI.
- `src/mhm_qgis/standalone/qgis_shim.py`: Qt-backed fallback for opening the plugin
  dialog without QGIS installed, installed via `mhm_qgis.standalone.install()`.
- `src/mhm_qgis/applications/`: one handler per external package -- everything the
  plugin asks of `mhm_tools` and `nml_tools` goes through these two modules.
- `src/mhm_qgis/qt/`: Qt code split by role -- `ui/forms`, `ui/pyui`, `bindings`,
  `controllers`. See `src/mhm_qgis/AGENTS.md` for the rule per layer.
- `src/mhm_qgis/qgis_bridge/`: drives the QGIS application (canvas, layers).
- `src/mhm_qgis/core/`: data handling that knows nothing about Qt or QGIS.
- `src/mhm_qgis/Morphology/latlon/ascii_morphology.py`: reusable morphology ASCII
  alignment and writing.
- `src/mhm_qgis/Meteorology/ERA5Land/mhm/`: local ERA5-Land forcing preparation.
- `setup.py`, `pyproject.toml`, `MANIFEST.in`: PyPI package metadata.

Read `src/mhm_qgis/AGENTS.md` before making plugin-internal changes. It contains the
current morphology, meteorology, configuration, and QGIS reload conventions.

## Dual Runtime Contract

The project supports two runtimes:

- QGIS plugin runtime: QGIS loads `src/mhm_qgis/__init__.py` and `src/mhm_qgis/mhm_qgis.py`.
  Real QGIS APIs, `QgsMapLayerComboBox`, QGIS Processing, and the QGIS project
  are available.
- Standalone PyPI runtime: the `mhm-qgis` command installs
  `mhm_qgis.standalone` before importing the dialog. Layer selectors become
  path selector widgets backed by lightweight layer objects.

Keep this rule intact:

```text
QGIS plugin code may use QGIS APIs.
Standalone startup must not require QGIS to be installed.
Reusable computation should move toward QGIS-free file/path/data APIs.
```

## Standalone GUI Notes

- `mhm_qgis.cli:main` is the console script target.
- The CLI intentionally calls `standalone_qgis.install(force=True)` so the
  command-line GUI uses file-path selectors even on machines where QGIS Python
  libraries are importable outside a real QGIS session.
- The standalone shim is for UI startup and input selection. It does not provide
  a real QGIS Processing backend. Processing actions that still call
  `processing.run(...)` need pure-Python or external-tool backends before they
  can be considered fully cluster-ready.
- Do not import `mhm_qgis.mhm_qgis_main` in standalone paths before installing the
  shim.

## mHM-Tools Integration Snapshot

- `mhm-tools>=0.3.0` is an external runtime dependency. There is no bundled
  `src/mhm_qgis/mhm_tools` source tree in the current project.
- The plugin calls exports from `mhm_tools.pre` through
  `src/mhm_qgis/applications/mhm_tools_handler.py`; it does not shell out to the
  `mhm-tools` CLI.
- `rasterize_map_data` burns vector attributes directly or maps them through a
  lookup table, always using the filled DEM's exact grid.
- `format_soil_data`, `format_geology_data`, `format_lai_data`, and
  `format_lc_data` accept
  categorical and DEM rasters in ASCII, GeoTIFF, or NetCDF form and align
  categories to the DEM with nearest-neighbour resampling.
- Shared raster file I/O and alignment live in `mhm_tools.common.file_handler`;
  `mhm_tools.common.rasterize` is intentionally limited to `rasterize_vector`.
- The corresponding external CLI commands are `data-converter rasterize-map`
  and `data-converter format-data`. Categorical inputs require lookup fields;
  gridded `--type lai` NetCDF input instead infers its source cadence and is
  warped directly onto the DEM grid.
- QGIS-specific layer selection, materialization, logging, output placement,
  and geology parameter metadata remain in `mhm_qgis`. See `src/mhm_qgis/AGENTS.md` for
  the detailed contract and current caveats.

## Domain And Worker Contracts

- Shared model inputs are published below `mhm-plugin/data/master/`. Static
  morphology, meteorology, LAI, lat/lon, and streamflow observations must not
  be written back to the former direct `data/*` folders.
- Each selected domain outlet owns exactly one physical DEM at
  `mhm-plugin/data/<outlet-id>/dem.asc`. Outlet IDs used as directory names
  must reject path separators and traversal components.
- `mhm_qgis_domain_delineation_state.json` is the canonical detailed state;
  `nml-settings.json` is the normalized handoff containing domain directories,
  DEM paths, contained gauge IDs, and shared `data/master/` paths.
- Qt and QGIS objects remain on the main thread. `QgsTask` jobs receive only
  copied primitive values and filesystem paths. Large advanced land-cover and
  soil formatting runs through `mhm_qgis.native_worker` in a child process so a
  GDAL/NumPy native failure cannot terminate QGIS.
- Historical v5.13 land-cover formatting is windowed and period-by-period in
  `mhm-tools`; never restore an eager list of full aligned period arrays.
- One common L0/L2 extent covers every active domain. Prepared L0 rasters and
  published advanced categorical inputs are placed on it by integer window copy
  and nodata padding, never resampled. See `src/mhm_qgis/AGENTS.md`.

## Packaging Rules

- Keep package discovery recursive. The project contains nested packages under
  `src/mhm_qgis/Morphology`, `src/mhm_qgis/Meteorology`, `src/mhm_qgis/qt`,
  `src/mhm_qgis/qt/ui`, `src/mhm_qgis/qt/bindings`, `src/mhm_qgis/qt/controllers`,
  `src/mhm_qgis/core`, `src/mhm_qgis/applications`, `src/mhm_qgis/qgis_bridge`,
  `src/mhm_qgis/standalone` and `src/mhm_qgis/vSpecific`. Every one needs an
  `__init__.py` or setuptools drops it from the wheel.
- Keep the `mhm-tools` dependency floor synchronized between `pyproject.toml`
  and `requirements.txt` when an adapter starts using a newer public API.
- Include plugin assets, schemas, templates, and project-template files in the
  PyPI build.
- If adding data files under `src/mhm_qgis/`, update `MANIFEST.in` and package-data
  patterns if the extension is new.
- Do not add `qgis` as a PyPI dependency. QGIS is provided by QGIS itself.

## Development Workflow

- Use `rg` / `rg --files` for search.
- Use `apply_patch` for manual edits.
- Do not touch unrelated user changes.
- Useful checks for this root-level packaging/runtime work:
  - `python3 -m py_compile src/mhm_qgis/standalone/cli.py src/mhm_qgis/standalone/qgis_shim.py`
  - `QT_QPA_PLATFORM=offscreen python3 -m mhm_qgis.cli --info`
  - A standalone dialog import/instantiation smoke test with the shim installed.
  - `python3 -m pytest -q` for the focused adapter and logging tests in `tests/`.
