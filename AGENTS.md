# AGENTS.md

Guidance for coding agents working at the repository root of `pymhm`.

## Repository Shape

This repository publishes one Python distribution and also contains the QGIS
plugin folder:

- `pymhm/`: QGIS plugin folder and Python package.
- `pymhm/metadata.txt`, `pymhm/__init__.py`, `pymhm/pymhm.py`: QGIS plugin
  entry points.
- `pymhm/cli.py`: PyPI console entry point; `pymhm` opens the standalone GUI.
- `pymhm/standalone_qgis.py`: Qt-backed fallback for opening the plugin dialog
  without QGIS installed.
- `setup.py`, `pyproject.toml`, `MANIFEST.in`: PyPI package metadata.

Read `pymhm/AGENTS.md` before making plugin-internal changes. It contains the
current morphology, meteorology, configuration, and QGIS reload conventions.

## Dual Runtime Contract

The project supports two runtimes:

- QGIS plugin runtime: QGIS loads `pymhm/__init__.py` and `pymhm/pymhm.py`.
  Real QGIS APIs, `QgsMapLayerComboBox`, QGIS Processing, and the QGIS project
  are available.
- Standalone PyPI runtime: the `pymhm` command installs
  `pymhm.standalone_qgis` before importing the dialog. Layer selectors become
  path selector widgets backed by lightweight layer objects.

Keep this rule intact:

```text
QGIS plugin code may use QGIS APIs.
Standalone startup must not require QGIS to be installed.
Reusable computation should move toward QGIS-free file/path/data APIs.
```

## Standalone GUI Notes

- `pymhm.cli:main` is the console script target.
- The CLI intentionally calls `standalone_qgis.install(force=True)` so the
  command-line GUI uses file-path selectors even on machines where QGIS Python
  libraries are importable outside a real QGIS session.
- The standalone shim is for UI startup and input selection. It does not provide
  a real QGIS Processing backend. Processing actions that still call
  `processing.run(...)` need pure-Python or external-tool backends before they
  can be considered fully cluster-ready.
- Do not import `pymhm.pymhm_dialog` in standalone paths before installing the
  shim.

## Packaging Rules

- Keep package discovery recursive. The project contains nested packages under
  `pymhm/Morphology`, `pymhm/Meteorology`, `pymhm/Configuration`, `pymhm/mhm_tools`,
  and `pymhm/mhm_tools_to_integrate`.
- Include plugin assets, schemas, templates, and project-template files in the
  PyPI build.
- If adding data files under `pymhm/`, update `MANIFEST.in` and package-data
  patterns if the extension is new.
- Do not add `qgis` as a PyPI dependency. QGIS is provided by QGIS itself.

## Development Workflow

- Use `rg` / `rg --files` for search.
- Use `apply_patch` for manual edits.
- Do not touch unrelated user changes.
- Useful checks for this root-level packaging/runtime work:
  - `python3 -m py_compile pymhm/cli.py pymhm/standalone_qgis.py`
  - `QT_QPA_PLATFORM=offscreen python3 -m pymhm.cli --info`
  - A standalone dialog import/instantiation smoke test with the shim installed.
