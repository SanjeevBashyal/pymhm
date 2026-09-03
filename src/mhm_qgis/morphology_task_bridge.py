"""Main-thread boundary for background morphology file jobs."""
from __future__ import annotations

import os
import json
import shutil
import subprocess
import sys
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory

from qgis.PyQt import QtCore
from qgis.PyQt.QtWidgets import QMessageBox
from qgis.core import QgsRasterLayer

from .core.morphology.soil_sources import local_layer_source
from .core.executions import morpho
from .core.handlers.raster.tasks import (dem_derivative_files, fill_dem_file,
                                    hydrology_files, terrain_files)
from .core.morphology.hydrology.outlets import configured_gauged_outlet_ids
from .qgis_bridge.morphology.layers.categorical import _SPECS
from .applications.mhm_tools_handler import prepare_lai_file
from .core.handlers.state.nml_settings import load_settings
from .core.handlers.state.cache import (
    cached_payload,
    fingerprint,
    outputs_newer_than_inputs,
    outputs_present,
    store_payload,
)
from .core.handlers.store.paths import geometry_folder, lai_dem_staging_path, morph_staging_folder, workspace_folder


def _saved_categorical_outputs(project_folder, kind):
    """Return complete configured outputs suitable for execute-all reuse."""
    section = "land_cover" if kind == "lc" else kind
    configured = load_settings(project_folder).get(section, {})
    if not isinstance(configured, dict):
        return ()
    values = []
    if kind == "lc":
        values.extend(
            scene.get("output_path")
            for scene in configured.get("scenes", ())
            if isinstance(scene, dict)
        )
    values.extend(
        configured.get(name) for name in ("output_path", "classdefinition_path")
    )
    root = Path(workspace_folder(project_folder))
    outputs = []
    for value in values:
        if not value:
            continue
        path = Path(value)
        path = path if path.is_absolute() else root / path
        if path not in outputs:
            outputs.append(path)
    return tuple(outputs) if outputs and all(path.is_file() for path in outputs) else ()


def _run_fill(task, options):
    return fill_dem_file(task=task, **options)


def _run_dem_derivatives(task, options):
    """Derive in the disposable worker: the pass peaks at several gigabytes."""
    def compute(prepared, staging, _log):
        result = _run_in_worker(
            task,
            {
                "kind": "dem-derivatives",
                "input_file": str(prepared),
                "output_folder": str(staging),
            },
            str(Path(staging).parent),
        )
        return {name: Path(path) for name, path in result["outputs"].items()}

    return dem_derivative_files(task=task, compute=compute, **options)


def _run_hydrology(task, options):
    return hydrology_files(task=task, **options)


def _run_terrain(task, options):
    return terrain_files(task=task, **options)


def _run_lai(task, options, messages):
    """Resample LAI onto the filled DEM grid inside the worker thread."""
    return prepare_lai_file(
        options["source_path"],
        options["filled_dem"],
        options["output_path"],
        output_temporal_resolution=options["target_timestep"],
        source_variable=options.get("source_variable"),
        dem_crs=options.get("dem_crs") or None,
        resampling=options.get("method") or "bilinear",
        task=task,
        log=messages.append,
    )


def _run_lookup(task, job):
    """Run one categorical lookup job in the disposable child process."""
    result = _run_in_worker(
        task,
        {"mode": "lookup", "kind": job["kind"], "job": job},
        job["geometry"],
    )
    return {
        "kind": result["kind"],
        "output": result["output"],
        "definition": result.get("definition", ""),
        "metadata": result.get("metadata", ""),
    }


def _run_advanced(task, job):
    payload = dict(job)
    payload["value"] = job["value"].as_dict()
    result = _run_in_worker(
        task, payload, geometry_folder(job["project_folder"]))
    return {"kind": result["kind"], "outputs": tuple(result["outputs"])}


def _run_in_worker(task, payload, working_folder):
    """Run one job in `mhm_qgis.native_worker` so a native failure cannot kill QGIS."""
    if task.isCanceled():
        raise RuntimeError("Task cancelled.")
    task.setProgress(5)
    with TemporaryDirectory(
            prefix="mhm_qgis_native_", dir=str(working_folder)) as temporary:
        temporary = Path(temporary)
        job_file = temporary / "job.json"
        result_file = temporary / "result.json"
        job_file.write_text(json.dumps(payload), encoding="utf-8")
        environment = os.environ.copy()
        package_parent = str(Path(__file__).resolve().parent.parent)
        environment["PYTHONPATH"] = os.pathsep.join(
            item
            for item in (package_parent, environment.get("PYTHONPATH", ""))
            if item
        )
        log_file = temporary / "worker.log"
        with log_file.open("w", encoding="utf-8") as log_stream:
            process = subprocess.Popen(
                [
                    sys.executable,
                    "-m",
                    "mhm_qgis.native_worker",
                    str(job_file),
                    str(result_file),
                ],
                stdout=log_stream,
                stderr=subprocess.STDOUT,
                text=True,
                env=environment,
            )
            while process.poll() is None:
                if task.isCanceled():
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                    raise RuntimeError("Task cancelled.")
                try:
                    process.wait(timeout=0.2)
                except subprocess.TimeoutExpired:
                    continue
        log = log_file.read_text(encoding="utf-8", errors="replace").strip()
        if not result_file.is_file():
            detail = log or f"exit status {process.returncode}"
            raise RuntimeError(f"Native raster worker failed: {detail}")
        payload = json.loads(result_file.read_text(encoding="utf-8"))
        if process.returncode or payload.get("error"):
            raise RuntimeError(
                payload.get("traceback") or payload.get("error") or log
            )
        result = payload["result"]
    task.setProgress(100)
    return result


class MorphologyTaskBridge(QtCore.QObject):
    """Schedule heavy work without passing Qt or QGIS objects to workers."""

    def __init__(self, dialog):
        super().__init__(dialog)
        self.dialog = dialog
        self.processor = dialog.morphology_processor
        self.coordinator = dialog.task_coordinator
        self.execute_all_active = False

    def _controls(self, controls):
        controls = tuple(controls)
        if controls:
            controls += (self.dialog.pushButton_BrowseProjectFolder,)
        return controls

    def _error(self, label, message, callback=None):
        detail = str(message).split("\n", 1)[0]
        self.dialog.log_message(f"ERROR: {label}: {detail}")
        if callback is not None:
            callback(detail)
        else:
            QMessageBox.critical(self.dialog, label, detail)

    def _dem_options(self):
        if not self.dialog.check_prerequisites():
            raise ValueError("Select a project folder and DEM first.")
        layer = self.dialog.input_combo("dem").currentLayer()
        source = str(layer.source() or "")
        if not source:
            raise ValueError("The selected DEM has no readable source.")
        source_crs = layer.crs()
        target_crs = self.dialog.get_crs()
        folder = Path(geometry_folder(self.dialog.project_folder))
        return {
            "source": source,
            "output_path": str(folder / "1_dem_filled.tif"),
            "source_crs": (
                source_crs.authid() or source_crs.toWkt()
                if source_crs is not None and source_crs.isValid()
                else ""
            ),
            "target_crs": (
                target_crs.authid() or target_crs.toWkt()
                if target_crs is not None and target_crs.isValid()
                else ""
            ),
            "reprojected_path": str(folder / "0_dem_reprojected.tif"),
        }

    def start_fill(
        self,
        *,
        key="fill-dem",
        controls=(),
        load=True,
        done=None,
        failed=None,
    ):
        try:
            options = self._dem_options()
        except Exception as error:
            self._error("Fill DEM", error, failed)
            return False
        started = self.coordinator.submit(
            key,
            "Fill DEM",
            partial(_run_fill, options=options),
            controls=self._controls(controls),
            resource="morphology-files",
            task_aware=True,
            on_success=lambda result: self._filled(result, load, done),
            on_error=lambda message: self._error("Fill DEM", message, failed),
        )
        return started

    def _filled(self, result, load, done):
        path = str(result["filled_dem"])
        self.processor.filled_dem_path = path
        self.processor.dem_layer = QgsRasterLayer(result["dem_source"], "Processing DEM")
        self.processor.mark_output_prepared(path, name="1_DEM_Filled", loaded=False)
        if load:
            self.processor.load_layer(path, "1_DEM_Filled")
        if not result.get("reused"):
            self.dialog.log_message(
                f'DEM filled successfully ({result.get("filled_cells", 0)} adjusted cells).'
            )
        self.dialog.update_l0_resolution_from_dem()
        if done is not None:
            done(result)

    def _hydrology_options(
        self,
        *,
        channel_path=None,
        threshold=None,
        direction=True,
        include_channel=True,
        write_flow=True,
        accumulation=True,
    ):
        folder = Path(geometry_folder(self.dialog.project_folder))
        return {
            "filled_dem": self.processor.filled_dem_path,
            "accumulation_path": (
                str(folder / "2_flow_accumulation.tif")
                if write_flow and accumulation
                else ""
            ),
            "area_path": (
                str(folder / "2_flow_accumulation_area.tif") if write_flow else ""
            ),
            "direction_path": str(folder / "2_flow_direction.tif") if direction else "",
            "channel_path": (
                str(channel_path or folder / "2_channel_network.shp")
                if include_channel
                else ""
            ),
            "threshold_cells": threshold,
        }

    def start_hydrology(
        self,
        *,
        key="prepare-hydrology",
        controls=(),
        load="",
        channel_path=None,
        threshold=None,
        direction=True,
        include_channel=True,
        write_flow=True,
        accumulation=True,
        done=None,
        failed=None,
    ):
        if not self.processor.filled_dem_path or not os.path.isfile(
            self.processor.filled_dem_path
        ):
            return self.start_fill(
                key=f"{key}-dem",
                controls=controls,
                load=False,
                done=lambda _result: self.start_hydrology(
                    key=key,
                    controls=controls,
                    load=load,
                    channel_path=channel_path,
                    threshold=threshold,
                    direction=direction,
                    include_channel=include_channel,
                    write_flow=write_flow,
                    accumulation=accumulation,
                    done=done,
                    failed=failed,
                ),
                failed=failed,
            )
        options = self._hydrology_options(
            channel_path=channel_path,
            threshold=threshold,
            direction=direction,
            include_channel=include_channel,
            write_flow=write_flow,
            accumulation=accumulation,
        )
        publish_channel = not channel_path
        return self.coordinator.submit(
            key,
            "Prepare flow and channel files",
            partial(_run_hydrology, options=options),
            controls=self._controls(controls),
            resource="morphology-files",
            task_aware=True,
            on_success=lambda result: self._hydrology_ready(
                result, load, done, publish_channel
            ),
            on_error=lambda message: self._error("Hydrology", message, failed),
        )

    def _hydrology_ready(self, result, load, done, publish_channel=True):
        mapping = {
            "flow_accumulation": ("flow_accumulation_path", "2_Flow_Accumulation", True),
            "flow_accumulation_area": ("flow_accumulation_area_path", "2_Flow_Accumulation_Area", True),
            "flow_direction": ("flow_direction_path", "2_Flow_Direction", True),
            "channel_network": ("channel_network_vector_path", "2_Channel_Network", False),
        }
        for key, path in result.items():
            if key == "channel_network" and not publish_channel:
                continue
            attribute, name, is_raster = mapping[key]
            setattr(self.processor, attribute, str(path))
            self.processor.mark_output_prepared(str(path), name=name, loaded=False)
            if load == key:
                self.processor.load_layer(str(path), name, is_raster=is_raster)
        if done is not None:
            done(result)

    def start_dem_derivatives(
        self,
        *,
        key="execute-all-dem-derivatives",
        controls=(),
        load=False,
        done=None,
        failed=None,
    ):
        """Fill the DEM and derive slope, aspect, facc and fdir in one pass."""
        try:
            options = self._dem_options()
        except Exception as error:
            self._error("DEM derivatives", error, failed)
            return False
        options.pop("output_path", None)
        options["output_folder"] = str(geometry_folder(self.dialog.project_folder))
        options["log"] = self.dialog.log_message
        return self.coordinator.submit(
            key,
            "Prepare DEM derivatives",
            partial(_run_dem_derivatives, options=options),
            controls=self._controls(controls),
            resource="morphology-files",
            task_aware=True,
            on_success=lambda result: self._dem_derivatives_ready(result, load, done),
            on_error=lambda message: self._error("DEM derivatives", message, failed),
        )

    def _dem_derivatives_ready(self, result, load, done):
        self._filled(result, load, None)
        self.processor.slope_path = result["slope"]
        self.processor.aspect_path = result["aspect"]
        self.processor.flow_accumulation_path = result["flow_accumulation"]
        self.processor.flow_direction_path = result["flow_direction"]
        for key in ("slope", "aspect", "flow_accumulation", "flow_direction"):
            self.processor.mark_output_prepared(
                result[key], name=Path(result[key]).name, loaded=False
            )
        if done is not None:
            done(result)

    def start_terrain(
        self,
        *,
        key="execute-all-terrain",
        controls=(),
        load="",
        done=None,
        failed=None,
    ):
        folder = Path(geometry_folder(self.dialog.project_folder))
        layer = self.dialog.input_combo("dem").currentLayer()
        scale = 111200.0 if layer.crs().isValid() and layer.crs().isGeographic() else 1.0
        options = {
            "filled_dem": self.processor.filled_dem_path,
            "slope_path": str(folder / "1_dem_slope.tif"),
            "aspect_path": str(folder / "1_dem_aspect.tif"),
            "scale": scale,
        }
        return self.coordinator.submit(
            key,
            "Prepare slope and aspect",
            partial(_run_terrain, options=options),
            controls=self._controls(controls),
            resource="morphology-files",
            task_aware=True,
            on_success=lambda result: self._terrain_ready(result, load, done),
            on_error=lambda message: self._error("Terrain", message, failed),
        )

    def _terrain_ready(self, result, load, done):
        self.processor.slope_path = result["slope"]
        self.processor.aspect_path = result["aspect"]
        for path in result.values():
            self.processor.mark_output_prepared(path, name=Path(path).name, loaded=False)
        if load in result:
            self.processor.load_layer(
                result[load], "1_DEM_Slope" if load == "slope" else "1_DEM_Aspect"
            )
        if done is not None:
            done(result)

    def start_lai(
        self,
        *,
        key="execute-all-lai",
        controls=(),
        done=None,
        failed=None,
    ):
        """Resample the selected LAI NetCDF onto the filled DEM grid."""
        processor = self.processor
        if not processor._is_lai_long_term_monthly_netcdf_selected():
            self.dialog.log_message(
                "No LAI NetCDF input is selected. Skipping LAI.")
            if done is not None:
                done(None)
            return True
        options = processor.lai_task_options(
            lai_dem_staging_path(self.dialog.project_folder)
        )
        if options is None:
            self.dialog.log_message(
                "LAI NetCDF input is not selected. Skipping LAI.")
            if done is not None:
                done(None)
            return True
        os.makedirs(os.path.dirname(options["output_path"]), exist_ok=True)
        digest = self._stage_digest(
            (options["source_path"], options["filled_dem"]),
            {
                "source_variable": options["source_variable"],
                "target_timestep": options["target_timestep"],
                "method": options.get("method") or "bilinear",
                "output_path": options["output_path"],
            },
        )
        if self._reusable_stage(
                "lai",
                digest,
                inputs=(options["source_path"], options["filled_dem"]),
                outputs=(options["output_path"],),
        ) is not None:
            self.dialog.log_message(
                "LAI inputs are unchanged. Reusing the staged DEM-grid LAI."
            )
            self.processor.mark_output_prepared(
                options["output_path"],
                name=Path(options["output_path"]).name,
                loaded=False,
            )
            if done is not None:
                QtCore.QTimer.singleShot(
                    0, lambda: done(options["output_path"]))
            return True
        self.dialog.log_message(
            "Resampling LAI to the filled DEM grid with bilinear interpolation."
        )
        messages = []
        return self.coordinator.submit(
            key,
            "Resample LAI to the filled DEM grid",
            partial(_run_lai, options=options, messages=messages),
            controls=self._controls(controls),
            resource="morphology-files",
            task_aware=True,
            on_success=lambda result: self._lai_ready(
                result, messages, done, digest=digest),
            on_error=lambda message: self._error("LAI", message, failed),
        )

    def _lai_ready(self, result, messages, done, digest=None):
        for message in messages:
            self.dialog.log_message(message)
        if result:
            self.processor.mark_output_prepared(
                result, name=Path(result).name, loaded=False)
            if digest:
                self._record_stage("lai", digest, (result,))
            self.dialog.log_message(
                f"LAI resampled to the filled DEM grid: {result}")
        if done is not None:
            done(result)

    def _stage_digest(self, inputs, config):
        """Return the fingerprint of one stage's inputs and configuration."""
        return fingerprint(inputs, config)

    def _reusable_stage(self, stage, digest, inputs=(), outputs=()):
        """Return the outputs to reuse, recording them when newly adopted.

        A recorded fingerprint that still matches is reused directly. Failing
        that, outputs already on disk are adopted when they postdate every input,
        which is how a project prepared before the cache existed avoids one
        pointless rebuild.
        """
        if not self.dialog.project_folder:
            return None
        payload = cached_payload(
            self.dialog.project_folder, "stages", stage, digest)
        if outputs_present(payload):
            return payload
        candidates = [str(path) for path in outputs if path]
        if not candidates or not all(Path(path).is_file() for path in candidates):
            return None
        if not outputs_newer_than_inputs(inputs, candidates):
            return None
        self._record_stage(stage, digest, candidates)
        self.dialog.log_message(
            f"Adopted the existing {stage} output into the project cache."
        )
        return {"outputs": candidates}

    def _record_stage(self, stage, digest, outputs):
        """Record which outputs a stage produced for the given inputs."""
        if not self.dialog.project_folder:
            return
        paths = [str(path) for path in outputs if path and Path(path).is_file()]
        if paths:
            store_payload(
                self.dialog.project_folder, "stages", stage, digest,
                {"outputs": paths},
            )

    def start_categorical(
        self,
        kind,
        *,
        key=None,
        controls=(),
        load=False,
        reuse_existing=False,
        done=None,
        failed=None,
    ):
        if not self.processor.filled_dem_path or not os.path.isfile(
            self.processor.filled_dem_path
        ):
            return self.start_fill(
                key=f"{key or kind}-dem",
                controls=controls,
                load=False,
                done=lambda _result: self.start_categorical(
                    kind,
                    key=key,
                    controls=controls,
                    load=load,
                    reuse_existing=reuse_existing,
                    done=done,
                    failed=failed,
                ),
                failed=failed,
            )
        if reuse_existing:
            outputs = _saved_categorical_outputs(self.dialog.project_folder, kind)
            if outputs:
                self.dialog.log_message(
                    f"Reusing prepared {_SPECS[kind]['label']} output."
                )
                self._categorical_complete(kind, outputs, reused=True)
                if done is not None:
                    QtCore.QTimer.singleShot(
                        0,
                        lambda: done(
                            {"kind": kind, "outputs": tuple(map(str, outputs))}
                        ),
                    )
                return True
        try:
            mode = self.processor._categorical_mode(kind)
            if mode == "mhm_ready":
                ok = self.processor._process_categorical(kind)
                if not ok:
                    raise RuntimeError(f"{_SPECS[kind]['label']} processing failed.")
                self._categorical_complete(
                    kind,
                    _saved_categorical_outputs(self.dialog.project_folder, kind),
                )
                if done is not None:
                    QtCore.QTimer.singleShot(0, lambda: done({"kind": kind}))
                return True
            if self.dialog.uses_advanced_categorical_input(kind):
                value = (
                    self.dialog._land_use_input_value(
                        self.dialog._advanced_inputs.get("land_cover")
                    )
                    if kind == "lc"
                    else self.dialog._soil_input_value(
                        self.dialog._advanced_inputs.get("soil")
                    )
                )
                job = {
                    "kind": kind,
                    "project_folder": str(self.dialog.project_folder),
                    "version": self.dialog.comboBox_mHMversion.currentText().strip(),
                    "value": value,
                    "filled_dem": self.processor.filled_dem_path,
                }
                runner = partial(_run_advanced, job=job)
                finish = lambda result: self._advanced_ready(result, done)
            else:
                if mode != "lookup_table":
                    raise ValueError(
                        f"Select an input type for {_SPECS[kind]['label']}."
                    )
                job = self._lookup_job(kind)
                digest = self._stage_digest(
                    (job["input_file"], job["filled_dem"], job["lookup_table"]),
                    {
                        "kind": kind,
                        "mapping_field": job["mapping_field"],
                        "class_field": job["class_field"],
                        "is_vector": job["is_vector"],
                        "input_crs": job["input_crs"],
                        "dem_crs": job["dem_crs"],
                    },
                )
                reusable = self._reusable_stage(
                    kind,
                    digest,
                    inputs=(
                        job["input_file"], job["filled_dem"], job["lookup_table"],
                    ),
                    outputs=(
                        job["output"], job["definition"], job["metadata"],
                    ),
                )
                if reusable is not None:
                    self.dialog.log_message(
                        f'{_SPECS[kind]["label"]} inputs are unchanged. '
                        "Reusing the prepared output."
                    )
                    self._lookup_ready(
                        {
                            "kind": kind,
                            "output": job["output"],
                            "definition": job["definition"],
                            "metadata": job["metadata"],
                        },
                        load,
                        None,
                        digest=digest,
                        reused=True,
                    )
                    if done is not None:
                        QtCore.QTimer.singleShot(0, lambda: done({"kind": kind}))
                    return True
                runner = partial(_run_lookup, job=job)
                finish = lambda result: self._lookup_ready(
                    result, load, done, digest=digest)
        except Exception as error:
            self._error(_SPECS[kind]["label"], error, failed)
            return False
        return self.coordinator.submit(
            key or f"categorical-{kind}",
            f'Prepare {_SPECS[kind]["label"]}',
            runner,
            controls=self._controls(controls),
            resource="morphology-files",
            task_aware=True,
            on_success=finish,
            on_error=lambda message: self._error(_SPECS[kind]["label"], message, failed),
        )

    def _lookup_job(self, kind):
        spec = _SPECS[kind]
        layer, source, is_vector = self.processor._categorical_input(kind, spec)
        config = self.processor._categorical_lookup(kind) or {}
        if not source or any(not config.get(key) for key in ("lookup_table", "mapping_field", "class_field")):
            raise ValueError(f"Configure the {spec['label']} input and lookup table first.")
        local = local_layer_source(layer) if layer is not None else None
        if local is None:
            candidate = Path(str(source).split("|", 1)[0])
            if candidate.is_file():
                local = str(candidate)
        if local is None:
            with TemporaryDirectory(prefix=f"mhm_qgis_{kind}_boundary_") as temporary:
                materialized = self.processor._categorical_input_file(
                    layer, source, is_vector, Path(temporary)
                )
                permanent = Path(geometry_folder(self.dialog.project_folder)) / f"_{kind}_input{Path(materialized).suffix}"
                shutil.copy2(materialized, permanent)
                local = str(permanent)
        geometry = Path(geometry_folder(self.dialog.project_folder))
        staging = Path(morph_staging_folder(self.dialog.project_folder))
        staging.mkdir(parents=True, exist_ok=True)
        output = geometry / spec["geometry"]
        definition = staging / spec["definition"] if spec.get("definition") else None
        metadata = geometry / "geology_class_metadata.json" if kind == "geology" else None
        removals = self.processor._ready_output_paths(kind)
        removals.extend(self.processor._categorical_intermediate_paths(kind)[1:])
        return {
            "kind": kind,
            "geometry": str(geometry),
            "input_file": str(local),
            "filled_dem": self.processor.filled_dem_path,
            "output": str(output),
            "definition": str(definition) if definition else "",
            "metadata": str(metadata) if metadata else "",
            "lookup_table": config["lookup_table"],
            "mapping_field": config["mapping_field"],
            "class_field": config["class_field"],
            "is_vector": bool(is_vector),
            "input_crs": self.processor._input_crs(layer),
            "dem_crs": self.processor._dem_crs(),
            "removals": tuple(map(str, removals)),
        }

    def _lookup_ready(self, result, load, done, digest=None, reused=False):
        kind = result["kind"]
        output = Path(result["output"])
        definition = Path(result["definition"]) if result["definition"] else None
        if kind == "lc":
            self.processor.land_use_layer = str(output)
        elif kind == "geology":
            self.processor.geology_path = str(output)
        self.processor.mark_output_prepared(str(output), name=output.name, loaded=False)
        if definition:
            self.processor.mark_output_prepared(str(definition), name=definition.name, loaded=False)
        if load:
            self.processor.load_layer(str(output), _SPECS[kind]["layer_name"])
        if kind == "soil":
            self.dialog.record_standard_soil_output(None, str(definition) if definition else None)
        self.processor._record_categorical_nml(kind, output, definition)
        outputs = [output]
        if definition:
            outputs.append(definition)
        if result.get("metadata"):
            outputs.append(Path(result["metadata"]))
        self._categorical_complete(kind, outputs, reused=reused)
        if digest and not reused:
            self._record_stage(kind, digest, outputs)
        if not reused:
            self.dialog.log_message(
                f'{_SPECS[kind]["label"]} data prepared successfully.')
        if done is not None:
            done(result)

    def _advanced_ready(self, result, done):
        for output in result["outputs"]:
            self.processor.mark_output_prepared(output, name=Path(output).name, loaded=False)
        self.dialog.log_message(
            f'Advanced {"land-cover" if result["kind"] == "lc" else "soil"} data prepared.'
        )
        self._categorical_complete(result["kind"], result["outputs"])
        if done is not None:
            done(result)

    def _categorical_complete(self, kind, outputs, reused=False):
        paths = [str(path) for path in outputs if Path(path).is_file()]
        for path in paths:
            self.processor.mark_output_prepared(
                path, name=Path(path).name, loaded=False
            )
        self.processor.mark_workflow_status(
            f"{'land_cover' if kind == 'lc' else kind}_processing",
            "completed",
            outputs=[self.processor.output_state_key(path) for path in paths],
            reused=bool(reused),
        )

    def start_domain_preflight(self, controls, done, failed=None):
        try:
            self._dem_options()
        except Exception as error:
            self._error("Domain Delineator", error, failed)
            return False
        folder = Path(geometry_folder(self.dialog.project_folder))
        self.processor.load_project_state()
        paths = (
            folder / "1_dem_filled.tif",
            folder / "2_flow_accumulation.tif",
            folder / "2_channel_network.shp",
        )
        if all(path.is_file() for path in paths):
            self.processor.filled_dem_path = str(paths[0])
            self.processor.flow_accumulation_path = str(paths[1])
            self.processor.channel_network_vector_path = str(paths[2])
            QtCore.QTimer.singleShot(0, lambda: done({"prepared": True}))
            return True
        return self.start_hydrology(
            key="domain-delineator-preflight",
            controls=controls,
            load="",
            direction=False,
            done=lambda result: done({"prepared": True, **result}),
            failed=failed,
        )

    def start_execute_all(self):
        if self.execute_all_active:
            return False
        if self.dialog.running_morphology_workflow_key() is not None or self.coordinator.resource_busy(
            "morphology-files"
        ):
            QMessageBox.information(
                self.dialog,
                "Execute All Processing",
                "Another preprocessing task is currently running.",
            )
            return False
        if not self.dialog.check_prerequisites():
            return False
        self.execute_all_active = True
        self.dialog.set_meteo_setup_controls_enabled(False)
        self.processor.skip_loading = True
        self.processor.load_processing_state()
        self.processor.mark_workflow_status("execute_all", "running")
        self.dialog.set_morphology_workflow_button_state("execute_all", "running")

        def start(stage, advance):
            starter = getattr(self, f"start_{stage.starter}")
            return starter(
                *stage.args,
                done=lambda _result: advance(),
                failed=self._fail_execute_all,
                **stage.options,
            )

        return morpho.run_stages(
            morpho.EXECUTE_ALL_STAGES,
            start,
            on_done=self._finish_execute_all,
            on_fail=self._fail_execute_all,
        )

    def _finish_execute_all(self):
        try:
            if self.dialog.input_combo("pour_points").currentLayer() is not None:
                self.processor.snap_points()
                if self.processor.snapped_points_path and os.path.isfile(
                    self.processor.snapped_points_path
                ):
                    configured = configured_gauged_outlet_ids(self.dialog.project_folder)
                    if configured != []:
                        self.processor.process_gauge_position()
            self.dialog.finish_morphology_workflow("execute_all", True, "")
        except Exception as error:
            self._fail_execute_all(str(error))
            return
        finally:
            self.processor.skip_loading = False
            self.dialog.set_meteo_setup_controls_enabled(True)
        self.execute_all_active = False

    def _fail_execute_all(self, message):
        self.processor.skip_loading = False
        self.dialog.set_meteo_setup_controls_enabled(True)
        if not self.execute_all_active:
            return
        self.execute_all_active = False
        self.dialog.finish_morphology_workflow("execute_all", False, str(message))

__all__ = ["MorphologyTaskBridge"]
