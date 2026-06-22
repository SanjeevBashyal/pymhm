# -*- coding: utf-8 -*-
"""LAI NetCDF preparation."""
from __future__ import annotations

import shutil

from ..common import (
    os,
    QMessageBox,
    QgsRasterLayer,
)
from ..watershed.dem_fill import DemFillMixin
from ..core.predecessors import PredecessorMixin
from ..core.layer_preparation import LayerPreparationMixin
from ...grid_resolution import header_bounds, write_header_file
from ...project_layout import geometry_folder, lai_folder


class LaiProcessingMixin(
        LayerPreparationMixin,
        DemFillMixin,
        PredecessorMixin):
    """Prepare gridded LAI inputs for mHM."""

    def process_lai(self) -> bool:
        """Process the selected LAI input according to the selected input type."""
        input_type = self._selected_lai_input_type()
        if input_type != "long term mean monthly gridded data (netcdf)":
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Only Long Term Mean Monthly Gridded Data (NetCDF) is implemented for LAI at the moment.")
            return False

        return self.process_lai_long_term_monthly_netcdf()

    def process_lai_long_term_monthly_netcdf(self) -> bool:
        """Clip and resample 12-month LAI NetCDF data to the filled DEM grid."""
        if not self.check_prerequisites():
            return False

        if not self._ensure_filled_dem(self.fill_dem):
            return False

        lai_layer = self._selected_lai_layer()
        if lai_layer is None:
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Please select a LAI NetCDF layer.")
            return False

        if not isinstance(lai_layer, QgsRasterLayer):
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Long term mean monthly LAI input must be a NetCDF raster layer.")
            return False

        source_path, source_variable = self._parse_lai_netcdf_source(lai_layer.source())
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

        output_folder = lai_folder(self.dialog.project_folder)
        os.makedirs(output_folder, exist_ok=True)
        output_path = os.path.join(output_folder, "lai.nc")

        self.log_message("\n--- Processing LAI NetCDF ---")
        self.log_message(f"Input LAI NetCDF: {source_path}")
        self.log_message(
            f"Filled DEM grid: {filled_dem_layer.width()} x {filled_dem_layer.height()} pixels")

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

        self.mark_output_prepared(output_path, name="lai.nc", loaded=False)
        self.log_message(f"LAI NetCDF prepared: {output_path}")
        return True

    def crop_lai_netcdf_to_l0(self, l0_header: dict, input_crs) -> bool:
        """Crop/resample the selected monthly LAI NetCDF to the L0 model grid."""
        if not self._is_lai_long_term_monthly_netcdf_selected():
            return False

        source = self._selected_lai_netcdf_source()
        if source is None:
            self.log_message("LAI NetCDF input is not selected. Skipping LAI crop.")
            return False
        source_path, source_variable = source

        output_folder = lai_folder(self.dialog.project_folder)
        os.makedirs(output_folder, exist_ok=True)
        crop_path = os.path.join(output_folder, "lai_crop.nc")
        masked_path = os.path.join(output_folder, "lai_masked.nc")
        final_path = os.path.join(output_folder, "lai.nc")

        self.log_message("\n--- Cropping LAI NetCDF to L0 grid ---")
        self.log_message(f"Input LAI NetCDF: {source_path}")
        self.log_message(
            f"LAI L0 grid: {int(l0_header['ncols'])} x {int(l0_header['nrows'])} cells")

        try:
            self._write_resampled_lai_netcdf_to_l0_header(
                source_path,
                source_variable,
                l0_header,
                input_crs,
                crop_path,
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
        """Apply the merged watershed mask to every monthly LAI slice."""
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

        self.log_message("\n--- Masking LAI NetCDF by merged watershed ---")
        try:
            mask = self._lai_watershed_mask_array(
                l0_header,
                input_crs,
                merged_watershed_path,
            )
            self._write_masked_lai_netcdf(crop_path, masked_path, final_path, mask)
            self._write_lai_l0_header(l0_header)
        except ImportError as e:
            self.log_message(f"ERROR: Required library not available for LAI mask: {e}")
            QMessageBox.warning(
                self.dialog,
                "Missing Dependency",
                f"Required library not available for LAI processing: {e}\n"
                "Please install: xarray, netCDF4, pyproj, and numpy.")
            return False
        except Exception as e:
            import traceback
            self.log_message(f"ERROR: LAI mask failed: {e}\n{traceback.format_exc()}")
            QMessageBox.warning(
                self.dialog,
                "LAI Mask Error",
                f"LAI NetCDF mask failed:\n{e}")
            return False

        self.mark_output_prepared(masked_path, name="lai_masked.nc", loaded=False)
        self.mark_output_prepared(final_path, name="lai.nc", loaded=False)
        self.log_message(f"LAI NetCDF masked for all months: {final_path}")
        return True

    def _selected_lai_input_type(self) -> str:
        """Return the normalized LAI input type selection."""
        combo_box = getattr(self.dialog, "comboBox_laiInputType", None)
        if combo_box is None:
            return ""
        return (combo_box.currentText() or "").strip().lower()

    def _is_lai_long_term_monthly_netcdf_selected(self) -> bool:
        """Return True when the active LAI workflow is the 12-month NetCDF input."""
        return (
            self._selected_lai_input_type()
            == "long term mean monthly gridded data (netcdf)"
        )

    def _selected_lai_layer(self):
        """Return the selected LAI layer from the LAI input combo box."""
        layer_combo = getattr(self.dialog, "mMapLayerComboBox_LAI_Class", None)
        if layer_combo is None:
            return None
        return layer_combo.currentLayer()

    def _selected_lai_netcdf_source(self) -> tuple[str, str | None] | None:
        """Return the selected LAI NetCDF path and variable, if available."""
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

    def _write_resampled_lai_netcdf(
            self,
            source_path: str,
            source_variable: str | None,
            filled_dem_layer,
            output_path: str) -> None:
        """Read, clip, resample, and write LAI data as a 12-month NetCDF cube."""
        target_grid = self._filled_dem_target_grid(filled_dem_layer)
        crs = filled_dem_layer.crs()
        if not crs.isValid():
            crs = self.dialog.get_crs()
        self._write_resampled_lai_netcdf_to_grid(
            source_path,
            source_variable,
            target_grid,
            crs,
            output_path,
            "Long term mean monthly LAI resampled to the filled DEM grid",
        )

    def _write_resampled_lai_netcdf_to_l0_header(
            self,
            source_path: str,
            source_variable: str | None,
            l0_header: dict,
            input_crs,
            output_path: str) -> None:
        """Read, crop, resample, and write LAI data on the L0 model grid."""
        target_grid = self._l0_header_target_grid(l0_header, input_crs)
        self._write_resampled_lai_netcdf_to_grid(
            source_path,
            source_variable,
            target_grid,
            input_crs,
            output_path,
            "Long term mean monthly LAI cropped and resampled to the L0 model grid",
        )

    def _write_resampled_lai_netcdf_to_grid(
            self,
            source_path: str,
            source_variable: str | None,
            target_grid,
            crs,
            output_path: str,
            description: str) -> None:
        """Read, clip, resample, and write LAI data for a target grid."""
        import numpy as np
        import xarray as xr

        dataset = xr.open_dataset(source_path)
        try:
            lat_coord = self._find_lai_coordinate(dataset, "lat")
            lon_coord = self._find_lai_coordinate(dataset, "lon")
            if lat_coord is None or lon_coord is None:
                raise ValueError("LAI NetCDF must contain 1D latitude and longitude coordinates.")

            coordinate_variables = [
                name for name in (lat_coord, lon_coord)
                if name in dataset.data_vars
            ]
            if coordinate_variables:
                dataset = dataset.set_coords(coordinate_variables)

            lat_dim = dataset[lat_coord].dims[0]
            lon_dim = dataset[lon_coord].dims[0]
            lai_data = self._find_lai_variable(
                dataset,
                source_variable,
                lat_dim,
                lon_dim,
            )
            lai_data, lat_dim, lon_dim = self._use_lai_coordinate_dimensions(
                lai_data,
                lat_coord,
                lon_coord,
                lat_dim,
                lon_dim,
            )
            month_dim = self._find_lai_month_dimension(lai_data, lat_dim, lon_dim)
            if month_dim is None:
                raise ValueError("LAI variable must contain one non-spatial dimension with exactly 12 months.")

            if lai_data.sizes[month_dim] != 12:
                raise ValueError(
                    f"LAI monthly dimension '{month_dim}' has {lai_data.sizes[month_dim]} values; expected 12.")
            extra_dims = [
                dim for dim in lai_data.dims
                if dim not in (month_dim, lat_dim, lon_dim)
            ]
            if extra_dims:
                raise ValueError(
                    "LAI variable has unsupported extra dimension(s): "
                    f"{', '.join(extra_dims)}.")

            lai_data = lai_data.sortby(lat_coord).sortby(lon_coord)
            x_centers, y_centers, target_lon, target_lat = target_grid
            sample_lon = self._match_lai_longitude_range(
                target_lon,
                lai_data[lon_coord].values,
            )

            subset = self._clip_lai_to_target_bounds(
                lai_data,
                lat_coord,
                lon_coord,
                target_lat,
                sample_lon,
            )
            if subset.sizes.get(lat_dim, 0) < 2 or subset.sizes.get(lon_dim, 0) < 2:
                raise ValueError("Selected LAI NetCDF does not cover the filled DEM extent.")

            target_lat_da = xr.DataArray(target_lat, dims=("yc", "xc"))
            target_lon_da = xr.DataArray(sample_lon, dims=("yc", "xc"))
            resampled = subset.interp(
                {lat_coord: target_lat_da, lon_coord: target_lon_da},
                method="linear",
            )
            nearest = subset.interp(
                {lat_coord: target_lat_da, lon_coord: target_lon_da},
                method="nearest",
            )
            resampled = resampled.where(~np.isnan(resampled), nearest)
            resampled = resampled.transpose(month_dim, "yc", "xc")
            lai_values = resampled.values.astype(np.float32)

            output_dataset = self._build_lai_output_dataset(
                lai_values,
                x_centers,
                y_centers,
                target_lat,
                target_lon,
                crs,
                lai_data,
                description,
            )
            self._save_lai_output_dataset(output_dataset, output_path)
        finally:
            dataset.close()

    def _find_lai_coordinate(self, dataset, coordinate_type: str) -> str | None:
        """Find a 1D latitude or longitude coordinate in a NetCDF dataset."""
        if coordinate_type == "lat":
            candidates = ("lat", "latitude", "y")
            standard_name = "latitude"
            axis = "Y"
        else:
            candidates = ("lon", "longitude", "x")
            standard_name = "longitude"
            axis = "X"

        for name in candidates:
            if name in dataset.variables and dataset[name].ndim == 1:
                return name

        for name in dataset.variables:
            variable = dataset[name]
            if variable.ndim != 1:
                continue
            attrs = variable.attrs
            if str(attrs.get("standard_name", "")).lower() == standard_name:
                return name
            if str(attrs.get("axis", "")).upper() == axis:
                return name

        return None

    def _find_lai_variable(
            self,
            dataset,
            source_variable: str | None,
            lat_dim: str,
            lon_dim: str):
        """Find the LAI data variable to process."""
        source_variable_candidates = []
        if source_variable:
            source_variable_candidates.append(source_variable)
            source_variable_basename = os.path.basename(source_variable)
            if source_variable_basename not in source_variable_candidates:
                source_variable_candidates.append(source_variable_basename)

        for candidate in source_variable_candidates:
            if candidate in dataset.data_vars:
                data_array = dataset[candidate]
                if lat_dim in data_array.dims and lon_dim in data_array.dims:
                    return data_array
                raise ValueError(
                    f"Selected NetCDF variable '{candidate}' does not use the detected latitude/longitude dimensions.")

        for candidate in ("lai", "LAI", "leaf_area_index", "Leaf_Area_Index"):
            if candidate in dataset.data_vars:
                data_array = dataset[candidate]
                if lat_dim in data_array.dims and lon_dim in data_array.dims:
                    return data_array

        for name, data_array in dataset.data_vars.items():
            if lat_dim not in data_array.dims or lon_dim not in data_array.dims:
                continue
            if self._find_lai_month_dimension(data_array, lat_dim, lon_dim):
                self.log_message(f"Using LAI variable '{name}'.")
                return data_array

        raise ValueError(
            "Could not find a LAI variable with latitude, longitude, and 12 monthly values.")

    def _use_lai_coordinate_dimensions(
            self,
            data_array,
            lat_coord: str,
            lon_coord: str,
            lat_dim: str,
            lon_dim: str):
        """Use latitude and longitude coordinates as interpolation dimensions."""
        swap_dims = {}
        if lat_coord != lat_dim:
            swap_dims[lat_dim] = lat_coord
        if lon_coord != lon_dim:
            swap_dims[lon_dim] = lon_coord

        if swap_dims:
            data_array = data_array.swap_dims(swap_dims)
            lat_dim = lat_coord
            lon_dim = lon_coord

        return data_array, lat_dim, lon_dim

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

    def _filled_dem_target_grid(self, filled_dem_layer):
        """Return filled DEM cell-center coordinates and their WGS84 lon/lat values."""
        import numpy as np
        from pyproj import Transformer

        extent = filled_dem_layer.extent()
        width = filled_dem_layer.width()
        height = filled_dem_layer.height()
        if width <= 0 or height <= 0:
            raise ValueError("Filled DEM has invalid dimensions.")

        cell_width = (extent.xMaximum() - extent.xMinimum()) / width
        cell_height = (extent.yMaximum() - extent.yMinimum()) / height
        x_centers = extent.xMinimum() + (np.arange(width, dtype=np.float64) + 0.5) * cell_width
        y_centers = extent.yMaximum() - (np.arange(height, dtype=np.float64) + 0.5) * cell_height

        crs = filled_dem_layer.crs()
        if not crs.isValid():
            crs = self.dialog.get_crs()
        crs_string = self._qgis_crs_to_pyproj(crs)
        if not crs_string:
            raise ValueError("Could not determine the filled DEM CRS.")

        x_grid, y_grid = np.meshgrid(x_centers, y_centers)
        transformer = Transformer.from_crs(crs_string, "EPSG:4326", always_xy=True)
        lon_array, lat_array = transformer.transform(x_grid, y_grid)

        return (
            x_centers.astype(np.float64),
            y_centers.astype(np.float64),
            np.asarray(lon_array, dtype=np.float64),
            np.asarray(lat_array, dtype=np.float64),
        )

    def _l0_header_target_grid(self, l0_header: dict, input_crs):
        """Return L0 cell-center coordinates and their WGS84 lon/lat values."""
        import numpy as np
        from pyproj import Transformer

        ncols = int(l0_header["ncols"])
        nrows = int(l0_header["nrows"])
        cellsize = float(l0_header["cellsize"])
        if ncols <= 0 or nrows <= 0 or cellsize <= 0:
            raise ValueError("L0 header has invalid dimensions or cell size.")

        xmin, xmax, ymin, ymax = header_bounds(l0_header)
        x_centers = xmin + (np.arange(ncols, dtype=np.float64) + 0.5) * cellsize
        y_centers = ymax - (np.arange(nrows, dtype=np.float64) + 0.5) * cellsize

        crs_string = self._qgis_crs_to_pyproj(input_crs)
        if not crs_string:
            raise ValueError("Could not determine the L0 grid CRS.")

        x_grid, y_grid = np.meshgrid(x_centers, y_centers)
        transformer = Transformer.from_crs(crs_string, "EPSG:4326", always_xy=True)
        lon_array, lat_array = transformer.transform(x_grid, y_grid)

        return (
            x_centers.astype(np.float64),
            y_centers.astype(np.float64),
            np.asarray(lon_array, dtype=np.float64),
            np.asarray(lat_array, dtype=np.float64),
        )

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

    def _match_lai_longitude_range(self, target_lon, source_lon_values):
        """Convert target longitudes to the source longitude convention."""
        import numpy as np

        source_lon = np.asarray(source_lon_values, dtype=np.float64)
        if source_lon.size and np.nanmin(source_lon) >= 0 and np.nanmax(source_lon) > 180:
            return np.mod(target_lon, 360.0)
        return target_lon

    def _clip_lai_to_target_bounds(
            self,
            lai_data,
            lat_coord: str,
            lon_coord: str,
            target_lat,
            target_lon):
        """Clip LAI data to the target bounds with a small buffer for interpolation."""
        import numpy as np

        lat_values = np.asarray(lai_data[lat_coord].values, dtype=np.float64)
        lon_values = np.asarray(lai_data[lon_coord].values, dtype=np.float64)
        lat_buffer = self._coordinate_buffer(lat_values)
        lon_buffer = self._coordinate_buffer(lon_values)

        lat_min = float(np.nanmin(target_lat)) - lat_buffer
        lat_max = float(np.nanmax(target_lat)) + lat_buffer
        lon_min = float(np.nanmin(target_lon)) - lon_buffer
        lon_max = float(np.nanmax(target_lon)) + lon_buffer

        return lai_data.sel({
            lat_coord: slice(lat_min, lat_max),
            lon_coord: slice(lon_min, lon_max),
        })

    def _coordinate_buffer(self, values) -> float:
        """Return a two-cell coordinate buffer for clipping before interpolation."""
        import numpy as np

        unique_values = np.unique(np.asarray(values, dtype=np.float64))
        if unique_values.size < 2:
            return 0.0
        diffs = np.diff(unique_values)
        diffs = np.abs(diffs[np.isfinite(diffs)])
        if diffs.size == 0:
            return 0.0
        return float(np.nanmedian(diffs) * 2.0)

    def _build_lai_output_dataset(
            self,
            lai_values,
            x_centers,
            y_centers,
            target_lat,
            target_lon,
            crs,
            source_lai_data,
            description: str):
        """Create the mHM-ready LAI output dataset."""
        import numpy as np
        import xarray as xr

        if crs is None or not crs.isValid():
            crs = self.dialog.get_crs()
        crs_string = self._qgis_crs_to_pyproj(crs) or ""

        dataset = xr.Dataset(
            data_vars={
                "lai": (("time", "yc", "xc"), lai_values),
                "lat": (("yc", "xc"), target_lat.astype(np.float64)),
                "lon": (("yc", "xc"), target_lon.astype(np.float64)),
            },
            coords={
                "time": np.arange(1, 13, dtype=np.int32),
                "yc": y_centers,
                "xc": x_centers,
            },
            attrs={
                "description": description,
                "projection": crs_string.lower(),
            },
        )

        dataset["time"].attrs.update({
            "long_name": "month",
            "units": "month",
        })
        dataset["yc"].attrs.update({
            "axis": "Y",
            "units": self._projected_axis_units(crs),
        })
        dataset["xc"].attrs.update({
            "axis": "X",
            "units": self._projected_axis_units(crs),
        })
        dataset["lat"].attrs.update({
            "units": "degrees_north",
            "long_name": "latitude",
        })
        dataset["lon"].attrs.update({
            "units": "degrees_east",
            "long_name": "longitude",
        })

        source_attrs = dict(source_lai_data.attrs)
        for reserved_attr in (
                "_FillValue",
                "missing_value",
                "scale_factor",
                "add_offset",
                "coordinates",
                "grid_mapping"):
            source_attrs.pop(reserved_attr, None)
        dataset["lai"].attrs.update(source_attrs)
        dataset["lai"].attrs.setdefault("long_name", "leaf area index")
        dataset["lai"].attrs.setdefault("units", "1")

        return dataset

    def _projected_axis_units(self, crs) -> str:
        """Return a simple axis unit label for output projected coordinates."""
        if crs is not None and crs.isGeographic():
            return "degrees"
        return "m"

    def _write_masked_lai_netcdf(
            self,
            crop_path: str,
            masked_path: str,
            final_path: str,
            watershed_mask) -> None:
        """Write LAI NetCDF outputs with the watershed mask applied to every month."""
        import numpy as np
        import xarray as xr

        with xr.open_dataset(crop_path) as source_dataset:
            dataset = source_dataset.load()

        if "lai" not in dataset:
            raise ValueError("Cropped LAI NetCDF does not contain variable 'lai'.")
        if dataset["lai"].sizes.get("yc") != watershed_mask.shape[0]:
            raise ValueError("LAI rows do not match the watershed mask rows.")
        if dataset["lai"].sizes.get("xc") != watershed_mask.shape[1]:
            raise ValueError("LAI columns do not match the watershed mask columns.")

        mask_da = xr.DataArray(
            np.asarray(watershed_mask, dtype=bool),
            dims=("yc", "xc"),
            coords={
                "yc": dataset["yc"],
                "xc": dataset["xc"],
            },
        )
        dataset["lai"] = dataset["lai"].where(mask_da, -9999.0).astype("float32")
        dataset["lai"].attrs.pop("_FillValue", None)
        dataset["lai"].attrs.pop("missing_value", None)
        dataset["lai"].attrs["nodata_value"] = -9999.0
        dataset.attrs["description"] = (
            "Long term mean monthly LAI cropped to L0 and masked by merged watershed"
        )

        self._save_lai_output_dataset(dataset, masked_path)
        if os.path.exists(final_path):
            os.remove(final_path)
        shutil.copyfile(masked_path, final_path)

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

    def _save_lai_output_dataset(self, dataset, output_path: str) -> None:
        """Write the LAI dataset to disk, falling back if compression is unavailable."""
        output_folder = os.path.dirname(output_path)
        tmp_path = os.path.join(output_folder, "lai_tmp.nc")
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

        encoding = {
            "lai": {
                "zlib": True,
                "complevel": 4,
                "_FillValue": -9999.0,
                "chunksizes": (1, dataset.sizes["yc"], dataset.sizes["xc"]),
            },
            "lat": {
                "zlib": True,
                "complevel": 4,
                "_FillValue": -9999.0,
                "chunksizes": (dataset.sizes["yc"], dataset.sizes["xc"]),
            },
            "lon": {
                "zlib": True,
                "complevel": 4,
                "_FillValue": -9999.0,
                "chunksizes": (dataset.sizes["yc"], dataset.sizes["xc"]),
            },
        }

        try:
            dataset.to_netcdf(tmp_path, encoding=encoding)
        except (TypeError, ValueError):
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            dataset.to_netcdf(tmp_path)

        if os.path.exists(output_path):
            os.remove(output_path)
        os.replace(tmp_path, output_path)
