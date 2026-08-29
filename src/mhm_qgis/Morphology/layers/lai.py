# -*- coding: utf-8 -*-
"""LAI NetCDF preparation."""
from __future__ import annotations

import shutil
from pathlib import Path

from ..common import (
    os,
    QMessageBox,
    QgsRasterLayer,
)
from ..watershed.dem_fill import DemFillMixin
from ..core.predecessors import PredecessorMixin
from ..core.layer_preparation import LayerPreparationMixin
from ...grid_resolution import header_bounds, write_header_file
from ...project_layout import (
    geometry_folder,
    lai_dem_staging_path,
    lai_folder,
)
from ...lai_temporal import lai_time_step, prepare_lai_temporal
from .lai_source import (
    LAI_ASSUMED_COMPRESSION,
    LAI_MAX_BYTES_ENV,
    assert_lai_output_fits,
    lai_grid_byte_size,
    resample_lai_file_to_grid,
    run_lai_resample,
    target_grid_from_header,
    target_grid_from_raster,
    window_copy_lai_file,
)


class LaiProcessingMixin(
        LayerPreparationMixin,
        DemFillMixin,
        PredecessorMixin):
    """Prepare gridded LAI inputs for mHM."""

    def process_lai(self) -> bool:
        """Process the selected LAI input according to the selected input type."""
        input_type = self._selected_lai_input_type()
        if "netcdf" not in input_type:
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Categorical LAI processing will be added through mHM-tools later.")
            return False

        return self.process_lai_long_term_monthly_netcdf()

    def process_lai_long_term_monthly_netcdf(self) -> bool:
        """Resample the selected LAI NetCDF onto the filled DEM grid.

        This is the Execute All stage. The result has exactly the filled DEM's
        size and extent, like every other prepared L0 layer, so Morphology Setup
        can crop it to the common extent by window copy instead of resampling.
        """
        if not self.check_prerequisites():
            return False

        if not self._ensure_filled_dem(self.fill_dem):
            return False

        source = self._selected_lai_netcdf_source()
        if source is None:
            QMessageBox.warning(self.dialog, "Input Error", "Select a LAI NetCDF file.")
            return False
        source_path, source_variable = source
        if not source_path:
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Could not determine the selected LAI NetCDF file path.")
            return False

        if not os.path.exists(source_path):
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                f"Selected LAI NetCDF file was not found:\n{source_path}")
            return False

        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog,
                "Error",
                "Cannot read filled DEM layer.")
            return False

        output_path = lai_dem_staging_path(self.dialog.project_folder)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        self.log_message("\n--- Resampling LAI NetCDF to the filled DEM grid ---")
        self.log_message(f"Input LAI NetCDF: {source_path}")
        self.log_message(
            f"Filled DEM grid: {filled_dem_layer.width()} x "
            f"{filled_dem_layer.height()} pixels")

        try:
            self._write_resampled_lai_netcdf(
                source_path,
                source_variable,
                filled_dem_layer,
                output_path,
            )
        except ImportError as e:
            self.log_message(f"ERROR: Required library not available: {e}")
            QMessageBox.warning(
                self.dialog,
                "Missing Dependency",
                f"Required library not available: {e}\n"
                "Please install: xarray, netCDF4, pyproj, and numpy.")
            return False
        except Exception as e:
            import traceback
            self.log_message(
                f"ERROR: LAI NetCDF preparation failed: {e}\n{traceback.format_exc()}")
            QMessageBox.critical(
                self.dialog,
                "LAI Processing Error",
                f"LAI NetCDF preparation failed:\n{e}")
            return False

        self.mark_output_prepared(output_path, name="lai_dem.nc", loaded=False)
        self.log_message(f"LAI resampled to the filled DEM grid: {output_path}")
        return True

    def resample_lai_to_dem_grid(self, show_error_dialog=True) -> bool:
        """Execute All hook: stage LAI on the filled DEM grid when selected."""
        if not self._is_lai_long_term_monthly_netcdf_selected():
            self.log_message("No LAI NetCDF input is selected. Skipping LAI.")
            return True
        if self._selected_lai_netcdf_source() is None:
            self.log_message("LAI NetCDF input is not selected. Skipping LAI.")
            return True
        return self.process_lai_long_term_monthly_netcdf()

    def crop_lai_netcdf_to_l0(self, l0_header: dict, input_crs) -> bool:
        """Window-copy the staged DEM-grid LAI onto the common L0 extent."""
        if not self._is_lai_long_term_monthly_netcdf_selected():
            return False

        if self._selected_lai_netcdf_source() is None:
            self.log_message("LAI NetCDF input is not selected. Skipping LAI crop.")
            return False

        staged_path = lai_dem_staging_path(self.dialog.project_folder)
        if not os.path.exists(staged_path):
            self.log_message(
                "ERROR: LAI has not been resampled to the filled DEM grid. "
                "Run Execute All before Morphology Setup."
            )
            return False

        output_folder = lai_folder(self.dialog.project_folder)
        os.makedirs(output_folder, exist_ok=True)
        crop_path = os.path.join(output_folder, "lai_crop.nc")
        masked_path = os.path.join(output_folder, "lai_masked.nc")
        final_path = os.path.join(output_folder, "lai.nc")

        self.log_message("\n--- Cropping LAI NetCDF to L0 grid ---")
        self.log_message(f"Staged DEM-grid LAI: {staged_path}")
        self.log_message(
            f"LAI L0 grid: {int(l0_header['ncols'])} x {int(l0_header['nrows'])} cells")

        try:
            window_copy_lai_file(
                staged_path,
                crop_path,
                l0_header,
                self._qgis_crs_to_pyproj(input_crs) or "",
                "Monthly LAI cropped to the common L0 model extent",
                log=self.log_message,
            )
            for stale_path in (masked_path, final_path):
                if os.path.exists(stale_path):
                    os.remove(stale_path)
            self._write_lai_l0_header(l0_header)
        except ImportError as e:
            self.log_message(f"ERROR: Required library not available for LAI crop: {e}")
            QMessageBox.warning(
                self.dialog,
                "Missing Dependency",
                f"Required library not available for LAI processing: {e}\n"
                "Please install: xarray, netCDF4, pyproj, and numpy.")
            return False
        except Exception as e:
            import traceback
            self.log_message(f"ERROR: LAI crop failed: {e}\n{traceback.format_exc()}")
            QMessageBox.warning(
                self.dialog,
                "LAI Crop Error",
                f"LAI NetCDF crop failed:\n{e}")
            return False

        self.mark_output_prepared(crop_path, name="lai_crop.nc", loaded=False)
        self.log_message(f"LAI NetCDF cropped to L0 grid: {crop_path}")
        return True

    def mask_lai_netcdf_to_l0(
            self,
            l0_header: dict,
            input_crs,
            merged_watershed_path: str) -> bool:
        """Publish the cropped LAI cube on the L0 grid without masking it.

        LAI is not cut to the merged watershed: masking would write nodata
        into cells mHM still reads, so the whole model extent keeps a usable
        leaf area. ``merged_watershed_path`` is accepted so this stage matches
        the other mask hooks, but it is deliberately unused.
        """
        if not self._is_lai_long_term_monthly_netcdf_selected():
            return False

        output_folder = lai_folder(self.dialog.project_folder)
        os.makedirs(output_folder, exist_ok=True)
        crop_path = os.path.join(output_folder, "lai_crop.nc")
        masked_path = os.path.join(output_folder, "lai_masked.nc")
        final_path = os.path.join(output_folder, "lai.nc")

        if not os.path.exists(crop_path):
            if not self.crop_lai_netcdf_to_l0(l0_header, input_crs):
                return False

        self.log_message("\n--- Writing LAI NetCDF on the L0 grid (unmasked) ---")
        try:
            window_copy_lai_file(
                crop_path,
                final_path,
                l0_header,
                self._qgis_crs_to_pyproj(input_crs) or "",
                "Monthly LAI on the common L0 model extent",
                log=self.log_message,
            )
            # The published cube is written straight from the crop; a separate
            # intermediate copy would double the disk cost for no benefit.
            if os.path.exists(masked_path):
                os.remove(masked_path)
            self._write_lai_l0_header(l0_header)
        except ImportError as e:
            self.log_message(f"ERROR: Required library not available for LAI: {e}")
            QMessageBox.warning(
                self.dialog,
                "Missing Dependency",
                f"Required library not available for LAI processing: {e}\n"
                "Please install: xarray, netCDF4, pyproj, and numpy.")
            return False
        except Exception as e:
            import traceback
            self.log_message(f"ERROR: LAI write failed: {e}\n{traceback.format_exc()}")
            QMessageBox.warning(
                self.dialog,
                "LAI Write Error",
                f"Writing the LAI NetCDF on the L0 grid failed:\n{e}")
            return False

        self.mark_output_prepared(final_path, name="lai.nc", loaded=False)
        source = self._selected_lai_netcdf_source()
        if source is not None:
            self._record_lai_nml(final_path, source[0], source[1])
        self.log_message(f"LAI NetCDF written for all months: {final_path}")
        return True

    def _selected_lai_input_type(self) -> str:
        """Return the normalized LAI input type selection."""
        combo_box = getattr(self.dialog, "comboBox_lai_inputType", None)
        if combo_box is None:
            return ""
        return (combo_box.currentText() or "").strip().lower()

    def _is_lai_long_term_monthly_netcdf_selected(self) -> bool:
        """Return True when the active LAI workflow uses NetCDF."""
        return "netcdf" in self._selected_lai_input_type()

    def _selected_lai_layer(self):
        """Return the selected LAI layer from the LAI input combo box."""
        layer_combo = self.dialog.input_combo("lai")
        if layer_combo is None:
            return None
        return layer_combo.currentLayer()

    def _selected_lai_netcdf_source(self) -> tuple[str, str | None] | None:
        """Return the selected LAI NetCDF path and variable, if available."""
        config_getter = getattr(self.dialog, "lai_netcdf_config", None)
        config = config_getter() if callable(config_getter) else {}
        configured = str(config.get("input_path", "") or "")
        if configured:
            if not os.path.isfile(configured):
                self.log_message(f"Selected LAI NetCDF file was not found: {configured}")
                return None
            return os.path.normpath(configured), None
        lai_layer = self._selected_lai_layer()
        if lai_layer is None:
            return None

        if not isinstance(lai_layer, QgsRasterLayer):
            self.log_message("Selected LAI layer is not a NetCDF raster layer.")
            return None

        source_path, source_variable = self._parse_lai_netcdf_source(lai_layer.source())
        if not source_path:
            self.log_message("Could not determine the selected LAI NetCDF file path.")
            return None

        if not os.path.exists(source_path):
            self.log_message(f"Selected LAI NetCDF file was not found: {source_path}")
            return None

        return source_path, source_variable

    def _parse_lai_netcdf_source(self, source: str) -> tuple[str | None, str | None]:
        """Extract a NetCDF path and variable name from a QGIS/GDAL source URI."""
        if not source:
            return None, None

        base_source = source.split("|")[0].strip()
        if base_source.upper().startswith("NETCDF:"):
            netcdf_source = base_source[len("NETCDF:"):]
            if netcdf_source.startswith(("\"", "'")):
                quote = netcdf_source[0]
                end_quote = netcdf_source.find(quote, 1)
                if end_quote > 0:
                    path = netcdf_source[1:end_quote]
                    tail = netcdf_source[end_quote + 1:]
                    variable = tail[1:] if tail.startswith(":") else None
                    return os.path.normpath(path), self._clean_lai_variable(variable)

            if os.path.exists(netcdf_source):
                return os.path.normpath(netcdf_source), None

            if ":" in netcdf_source:
                path, variable = netcdf_source.rsplit(":", 1)
                return os.path.normpath(path), self._clean_lai_variable(variable)

            return os.path.normpath(netcdf_source), None

        return os.path.normpath(base_source), None

    def _clean_lai_variable(self, variable: str | None) -> str | None:
        """Normalize a variable name parsed from a NetCDF source URI."""
        if not variable:
            return None
        variable = variable.strip().strip("\"'")
        return variable or None

    def lai_task_options(self, output_path: str) -> dict | None:
        """Snapshot the primitives one LAI resample task needs.

        Called on the main thread. Everything a worker needs is a path, a string
        or a number, so the task never touches a widget or a QGIS layer.
        """
        source = self._selected_lai_netcdf_source()
        if source is None:
            return None
        source_path, source_variable = source
        getter = getattr(self.dialog, "lai_netcdf_config", None)
        config = getter() if callable(getter) else {}
        return {
            "source_path": str(source_path),
            "source_variable": source_variable,
            "output_path": str(output_path),
            "filled_dem": str(self.filled_dem_path),
            "crs_string": self._qgis_crs_to_pyproj(self.dialog.get_crs()) or "",
            "input_resolution": config.get("input_resolution") or "Monthly",
            "target_timestep": (
                config.get("target_timestep")
                or "Long Term Mean Monthly Gridded Data"
            ),
        }

    def _write_resampled_lai_netcdf(
            self,
            source_path: str,
            source_variable: str | None,
            filled_dem_layer,
            output_path: str) -> None:
        """Resample LAI onto the filled DEM grid through the QGIS-free pipeline."""
        options = self.lai_task_options(output_path)
        if options is None:
            raise ValueError("LAI NetCDF input is not selected.")
        run_lai_resample(options, log=self.log_message)

    def _find_lai_month_dimension(self, data_array, lat_dim: str, lon_dim: str) -> str | None:
        """Find the 12-value monthly dimension in the LAI variable."""
        non_spatial_dims = [
            dim for dim in data_array.dims if dim not in (lat_dim, lon_dim)
        ]
        preferred_dims = [
            dim for dim in non_spatial_dims
            if dim.lower() in ("time", "month", "months", "mon", "t")
        ]

        for dim in preferred_dims + non_spatial_dims:
            if data_array.sizes.get(dim) == 12:
                return dim
        return None

    def _qgis_crs_to_pyproj(self, crs) -> str | None:
        """Return a pyproj-readable CRS string from a QGIS CRS."""
        if crs is None or not crs.isValid():
            return None
        if crs.postgisSrid():
            return f"EPSG:{crs.postgisSrid()}"
        if crs.authid():
            return crs.authid()
        if hasattr(crs, "toWkt"):
            return crs.toWkt()
        return None

    def _lai_time_attrs(time_step: int) -> dict:
        """Return the time-axis attributes for a prepared LAI time step."""
        if time_step == 1:
            return {"long_name": "month", "units": "month"}
        return {"standard_name": "time", "axis": "T"}

    def _record_lai_nml(self, output_path, source_path, source_variable):
        from ...nml_settings import relative_workspace_path, update_section

        getter = getattr(self.dialog, "lai_netcdf_config", None)
        config = getter() if callable(getter) else {}
        target = config.get("target_timestep") or "Long Term Mean Monthly Gridded Data"
        update_section(
            self.dialog.project_folder,
            "lai",
            {
                "mode": "netcdf",
                "source_path": str(Path(source_path).resolve()),
                "source_variable": str(source_variable or ""),
                "input_resolution": config.get("input_resolution") or "Monthly",
                "target_timestep": target,
                "time_step": lai_time_step(target),
                "output_path": relative_workspace_path(
                    self.dialog.project_folder, output_path
                ),
                "variable": "lai",
            },
        )

    def _lai_watershed_mask_array(
            self,
            l0_header: dict,
            input_crs,
            merged_watershed_path: str):
        """Rasterize the merged watershed to a boolean array on the L0 grid."""
        import numpy as np
        from osgeo import gdal, osr

        ncols = int(l0_header["ncols"])
        nrows = int(l0_header["nrows"])
        cellsize = float(l0_header["cellsize"])
        xmin, xmax, ymin, ymax = header_bounds(l0_header)
        mask_path = os.path.join(
            geometry_folder(self.dialog.project_folder),
            "lai_watershed_mask.tif",
        )
        os.makedirs(os.path.dirname(mask_path), exist_ok=True)
        if os.path.exists(mask_path):
            os.remove(mask_path)

        driver = gdal.GetDriverByName("GTiff")
        dataset = driver.Create(
            mask_path,
            ncols,
            nrows,
            1,
            gdal.GDT_Byte,
            options=["COMPRESS=LZW"],
        )
        if dataset is None:
            raise RuntimeError(f"Could not create LAI watershed mask: {mask_path}")

        dataset.SetGeoTransform((xmin, cellsize, 0.0, ymax, 0.0, -cellsize))
        spatial_ref = osr.SpatialReference()
        if input_crs is not None and input_crs.isValid():
            if input_crs.postgisSrid():
                spatial_ref.ImportFromEPSG(int(input_crs.postgisSrid()))
            elif input_crs.authid().upper().startswith("EPSG:"):
                spatial_ref.ImportFromEPSG(int(input_crs.authid().split(":")[1]))
            elif hasattr(input_crs, "toWkt"):
                spatial_ref.ImportFromWkt(input_crs.toWkt())
            dataset.SetProjection(spatial_ref.ExportToWkt())

        band = dataset.GetRasterBand(1)
        band.SetNoDataValue(0)
        band.Fill(0)

        vector_dataset = gdal.OpenEx(merged_watershed_path, gdal.OF_VECTOR)
        if vector_dataset is None:
            dataset = None
            raise RuntimeError(
                f"Could not open merged watershed for LAI masking: {merged_watershed_path}")
        vector_layer = vector_dataset.GetLayer(0)
        if vector_layer is None:
            dataset = None
            vector_dataset = None
            raise RuntimeError("Merged watershed has no readable vector layer.")

        error_code = gdal.RasterizeLayer(
            dataset,
            [1],
            vector_layer,
            burn_values=[1],
        )
        if error_code != 0:
            dataset = None
            vector_dataset = None
            raise RuntimeError("GDAL failed to rasterize merged watershed for LAI masking.")

        band.FlushCache()
        mask = band.ReadAsArray()
        band = None
        dataset = None
        vector_dataset = None

        if mask is None:
            raise RuntimeError("Could not read LAI watershed mask array.")
        if mask.shape != (nrows, ncols):
            raise ValueError(
                f"LAI mask shape {mask.shape} does not match L0 shape {(nrows, ncols)}.")
        return np.asarray(mask) > 0

    def _write_lai_l0_header(self, l0_header: dict) -> None:
        """Write the LAI L0 header required by mHM."""
        output_folder = lai_folder(self.dialog.project_folder)
        os.makedirs(output_folder, exist_ok=True)
        header_path = os.path.join(output_folder, "header.txt")
        header = dict(l0_header)
        header["nodata_value"] = -9999.0
        write_header_file(header_path, header)
        self.mark_output_prepared(header_path, name="header.txt", loaded=False)

