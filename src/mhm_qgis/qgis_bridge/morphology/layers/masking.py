# -*- coding: utf-8 -*-
"""Cropping and watershed masking for prepared morphology rasters."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
)
from ..watershed.watershed_delineation import WatershedDelineationMixin
from ....grid_resolution import (
    CATEGORICAL_PAD_VALUE,
    header_bounds,
    l0_header_from_l2,
    validate_l0_l2_alignment,
)
from ....core.handlers.raster.tasks import crop_aligned_l0_raster, mask_aligned_l0_raster


class MaskingMixin(WatershedDelineationMixin):
    """Cropping and watershed masking for prepared morphology rasters."""

    def crop_all_layers(self, show_error_dialog=True) -> bool:
        """Crop prepared morphology rasters to the common L0/L2 model extent."""
        if not self.check_prerequisites():
            return False

        if not self._ensure_filled_dem(self.fill_dem):
            return False

        geometry_folder = project_geometry_folder(self.session.project_folder)
        try:
            l0_header = self._target_l0_header()
        except Exception as e:
            self.log_message(f"ERROR: Cannot crop morphology rasters: {e}")
            if show_error_dialog:
                self.warn("Grid Configuration Error", str(e))
            return False

        input_crs = self.session.crs
        if not input_crs.isValid():
            if show_error_dialog:
                self.warn(
                                        "CRS Error",
                    "Please set a valid input CRS.")
            return False

        layers_to_crop = self._collect_layers_to_crop(geometry_folder)
        if not layers_to_crop:
            lai_cropped = self._crop_lai_netcdf_if_available(l0_header, input_crs)
            if lai_cropped:
                self.log_message("Crop all process completed.")
                return True
            if show_error_dialog:
                self.warn(
                                        "No Layers",
                    "No raster layers found to crop. Please process at least the filled DEM first.")
            return False

        extent_str = self._extent_from_header(l0_header, input_crs)
        resolution = float(l0_header["cellsize"])

        self.log_message("\n--- Cropping all morphology rasters ---")
        self.log_message(f"Crop extent: {extent_str}")
        self.log_message(f"Crop resolution: {resolution:.6f}")
        self.log_message(f"Cropping {len(layers_to_crop)} raster layer(s)...")

        all_success = True
        for layer_info in layers_to_crop:
            input_path = layer_info["input"]
            output_path = layer_info["crop"]
            masked_path = layer_info["masked"]
            layer_name = layer_info["name"]

            self.log_message(f"Cropping {layer_name}...")
            self._remove_existing_raster(output_path)
            self._remove_existing_raster(masked_path)
            try:
                result = crop_aligned_l0_raster(
                    input_path,
                    output_path,
                    l0_header,
                    reference_path=self.filled_dem_path,
                    pad_value=layer_info.get("pad"),
                )
            except Exception as error:
                self.log_message(f"ERROR: Failed to crop {layer_name}: {error}")
                result = None

            if result and os.path.exists(output_path):
                self.mark_output_prepared(
                    output_path,
                    name=os.path.basename(output_path),
                    loaded=False,
                    algorithm="aligned L0 window copy",
                )
                self.log_message(f"{layer_name} cropped successfully.")
            else:
                self.log_message(f"ERROR: Failed to crop {layer_name}.")
                all_success = False

        lai_selected = bool(
            getattr(self, "_is_lai_long_term_monthly_netcdf_selected", lambda: False)()
        )
        lai_cropped = self._crop_lai_netcdf_if_available(l0_header, input_crs)
        if lai_selected and not lai_cropped:
            all_success = False
        self.log_message("Crop all process completed.")
        return all_success

    def mask_all_layers(self, show_error_dialog=True) -> bool:
        """Mask the cropped DEM derivatives by merged watershed, -9999 outside.

        Land cover, soil, geology, and LAI are not masked: they are written to
        their `_masked` output as a straight copy of the crop, so they keep a
        usable value across the whole model extent.
        """
        if not self.check_prerequisites():
            return False

        geometry_folder = project_geometry_folder(self.session.project_folder)
        try:
            l0_header = self._target_l0_header()
        except Exception as e:
            self.log_message(f"ERROR: Cannot mask morphology rasters: {e}")
            if show_error_dialog:
                self.warn("Grid Configuration Error", str(e))
            return False

        input_crs = self.session.crs
        if not input_crs.isValid():
            if show_error_dialog:
                self.warn(
                                        "CRS Error",
                    "Please set a valid input CRS.")
            return False

        merged_watershed_path = self._restore_existing_path(
            "merged_watershed_path",
            os.path.join("Watersheds", "4_watershed_merged_vector.shp"),
            "4_watershed_merged_vector.shp",
        )

        if not merged_watershed_path or not os.path.exists(merged_watershed_path):
            if not self._ensure_merged_watershed(
                    self.delineate_watershed,
                    self.snap_points,
                    self.process_channel_network,
                    self.process_flow_accumulation,
                    self.fill_dem):
                return False
            merged_watershed_path = self.merged_watershed_path

        if not merged_watershed_path or not os.path.exists(merged_watershed_path):
            return False

        layers_to_mask = self._collect_layers_to_mask(geometry_folder)
        if not layers_to_mask:
            lai_masked = self._mask_lai_netcdf_if_available(
                l0_header,
                input_crs,
                merged_watershed_path,
            )
            if lai_masked:
                self.log_message("Masking process completed.")
                return True
            if show_error_dialog:
                self.warn(
                                        "No Cropped Layers",
                    "No cropped raster layers found. Please run Crop All first.")
            return False

        extent_str = self._extent_from_header(l0_header, input_crs)
        resolution = float(l0_header["cellsize"])

        self.log_message("\n--- Masking all cropped morphology rasters ---")
        self.log_message(f"Mask extent: {extent_str}")
        self.log_message(f"Mask resolution: {resolution:.6f}")
        self.log_message(f"Masking {len(layers_to_mask)} raster layer(s)...")

        all_success = True
        for layer_info in layers_to_mask:
            input_path = layer_info["input"]
            output_path = layer_info["output"]
            layer_name = layer_info["name"]
            watershed_masked = layer_info.get("watershed_masked", True)
            verb = "Masking" if watershed_masked else "Copying unmasked"

            self.log_message(f"{verb} {layer_name}...")
            self._remove_existing_raster(output_path)
            try:
                if watershed_masked:
                    result = mask_aligned_l0_raster(
                        input_path,
                        output_path,
                        l0_header,
                        merged_watershed_path,
                        reference_path=self.filled_dem_path,
                        pad_value=layer_info.get("pad"),
                    )
                else:
                    # Class layers keep the full model extent: the watershed
                    # boundary would only carve nodata into cells mHM reads.
                    result = crop_aligned_l0_raster(
                        input_path,
                        output_path,
                        l0_header,
                        reference_path=self.filled_dem_path,
                        pad_value=layer_info.get("pad"),
                    )
            except Exception as error:
                self.log_message(f"ERROR: Failed to write {layer_name}: {error}")
                result = None

            if result and os.path.exists(output_path):
                self.mark_output_prepared(
                    output_path,
                    name=os.path.basename(output_path),
                    loaded=False,
                    algorithm=(
                        "aligned L0 window mask" if watershed_masked
                        else "aligned L0 window copy"
                    ),
                )
                self.log_message(
                    f"{layer_name} "
                    f"{'masked' if watershed_masked else 'written unmasked'} "
                    "successfully."
                )
            else:
                self.log_message(f"ERROR: Failed to write {layer_name}.")
                all_success = False

        lai_selected = bool(
            getattr(self, "_is_lai_long_term_monthly_netcdf_selected", lambda: False)()
        )
        lai_masked = self._mask_lai_netcdf_if_available(
            l0_header,
            input_crs,
            merged_watershed_path,
        )
        if lai_selected and not lai_masked:
            all_success = False
        self.log_message("Masking process completed.")
        return all_success

    def align_advanced_inputs_to_l0(self, show_error_dialog=True) -> bool:
        """Publish the staged morphology inputs into data/master on the L0 grid."""
        from ....core.morphology.layers.advanced_l0 import missing_model_inputs, publish_model_inputs

        try:
            l0_header = self._target_l0_header()
            published = publish_model_inputs(
                self.session.project_folder,
                l0_header,
                log=self.log_message,
            )
        except Exception as error:
            self.log_message(f"ERROR: Cannot publish model inputs: {error}")
            if show_error_dialog:
                self.warn("Publish Error", str(error))
            return False

        for path in published:
            self.mark_output_prepared(
                str(path), name=os.path.basename(str(path)), loaded=False)
        if not published:
            self.log_message("No staged morphology inputs were waiting to publish.")

        # Everything the namelist handoff names must now exist under data/master.
        missing = missing_model_inputs(self.session.project_folder)
        if missing:
            message = (
                "These model inputs are recorded in nml-settings.json but are "
                "missing on disk: " + ", ".join(sorted(missing))
            )
            self.log_message(f"ERROR: {message}")
            if show_error_dialog:
                self.warn("Missing Model Inputs", message)
            return False
        self.log_message(
            "All model inputs recorded in nml-settings.json are present.")
        return True

    def write_domain_dems_to_l0(self, show_error_dialog=True) -> bool:
        """Write every domain DEM by masking the cropped L0 DEM."""
        from ....core.morphology.layers.domain_dem import domain_dem_plan, write_domain_dems

        try:
            l0_header = self._target_l0_header()
            cropped = self._raster_variant_path(self.filled_dem_path, "_crop")
            plan = domain_dem_plan(self.session.project_folder)
            self.save_domain_plan(plan)
            written = write_domain_dems(
                self.session.project_folder,
                cropped,
                l0_header,
                reference_path=self.filled_dem_path,
                log=self.log_message,
            )
        except Exception as error:
            self.log_message(f"ERROR: Cannot write the domain DEMs: {error}")
            if show_error_dialog:
                self.warn("Domain DEM Error", str(error))
            return False

        for entry in written:
            self.mark_output_prepared(
                entry["dem_path"], name="dem.asc", loaded=False)
            for path in entry.get("copied", ()):
                self.mark_output_prepared(
                    path, name=os.path.basename(path), loaded=False)
        if not written:
            self.log_message("No active domains needed a DEM.")
        return True

    def _target_l0_header(self):
        """Return the L0 grid header derived from the configured model extent."""
        l2_header = getattr(self.dialog, "_grid_l2_header", None)
        if not l2_header and hasattr(self.dialog, "update_l2_resolution_from_metadata"):
            self.dialog.update_l2_resolution_from_metadata()
            l2_header = getattr(self.dialog, "_grid_l2_header", None)
        if not l2_header:
            saved = self.saved_grid_contract()
            if saved:
                return dict(saved["l0_header"])
            raise ValueError("L2 grid is not available. Run and save meteorology data first.")

        l0_resolution = None
        if hasattr(self.dialog, "current_l0_resolution"):
            l0_resolution = self.session.l0_resolution
        if not l0_resolution:
            l0_info = getattr(self.dialog, "_grid_l0_info", None)
            if not l0_info and hasattr(self.dialog, "filled_dem_resolution_info"):
                l0_info = self.session.dem_resolution_info
            if l0_info:
                l0_resolution = float(l0_info["resolution"])
        if not l0_resolution:
            raise ValueError("L0 resolution is not available. Select or prepare the filled DEM first.")

        metadata = getattr(self.dialog, "_grid_l2_metadata", None) or {}
        ratio = int(
            metadata.get("l2_ratio_to_l0")
            or round(float(l2_header["cellsize"]) / float(l0_resolution))
        )
        # The saved L0 header carries the filled DEM's unrounded cell size.
        # Rebuilding it from the rounded UI resolution would drift off the
        # source grid, so prefer the saved header whenever it still validates.
        saved_l0_header = metadata.get("l0_header")
        l0_header = None
        if isinstance(saved_l0_header, dict):
            try:
                validate_l0_l2_alignment(saved_l0_header, l2_header, ratio)
                l0_header = dict(saved_l0_header)
            except (TypeError, ValueError) as error:
                self.log_message(
                    f"WARNING: Saved L0 header does not match the L2 grid: {error}"
                )
        if l0_header is None:
            l0_header = l0_header_from_l2(l2_header, float(l0_resolution), ratio)
        self.save_grid_contract(l0_header, l2_header, ratio)
        return l0_header

    def _extent_from_header(self, header, input_crs):
        """Return a QGIS Processing extent string from a grid header."""
        xmin, xmax, ymin, ymax = header_bounds(header)
        return f"{xmin},{xmax},{ymin},{ymax} [{input_crs.authid()}]"

    def _crop_lai_netcdf_if_available(self, l0_header, input_crs):
        """Crop LAI NetCDF data when the long-term monthly LAI input is selected."""
        crop_lai = getattr(self, "crop_lai_netcdf_to_l0", None)
        if crop_lai is None:
            return False
        return bool(crop_lai(l0_header, input_crs))

    def _mask_lai_netcdf_if_available(
            self,
            l0_header,
            input_crs,
            merged_watershed_path):
        """Mask LAI NetCDF data when the long-term monthly LAI input is selected."""
        mask_lai = getattr(self, "mask_lai_netcdf_to_l0", None)
        if mask_lai is None:
            return False
        return bool(mask_lai(l0_header, input_crs, merged_watershed_path))

    def _collect_layers_to_crop(self, geometry_folder):
        """Return prepared base morphology rasters that should be cropped."""
        # The last two fields are the value written where the layer does not
        # reach the model extent, and whether the watershed mask applies.
        # Only the DEM derivatives are masked. Class layers pad with a valid
        # class and keep the full extent, because mHM needs a class in every
        # cell it reads and nodata there would be a hole, not a boundary.
        layer_specs = [
            ("filled_dem_path", "1_dem_filled.tif", "1_DEM_Filled", None, True),
            ("aspect_path", "1_dem_aspect.tif", "1_DEM_Aspect", None, True),
            ("slope_path", "1_dem_slope.tif", "1_DEM_Slope", None, True),
            ("flow_accumulation_path", "2_flow_accumulation.tif", "2_Flow_Accumulation", None, True),
            ("flow_direction_path", "2_flow_direction.tif", "2_Flow_Direction", None, True),
            ("gauge_position_path", "2_gauge_position.tif", "2_Gauge_Position", None, True),
            ("land_use_layer", "3_land_use.tif", "3_Land_Use", CATEGORICAL_PAD_VALUE, False),
            (None, "3_soil.tif", "3_Soil", CATEGORICAL_PAD_VALUE, False),
            ("geology_path", "3_geology_processed.tif", "3_Geology", CATEGORICAL_PAD_VALUE, False),
        ]

        layers_to_crop = []
        seen = set()
        for attr_name, filename, layer_name, pad_value, watershed_masked in layer_specs:
            expected_path = os.path.join(geometry_folder, filename)
            layer_path = getattr(self, attr_name, None) if attr_name else None
            if not layer_path or not os.path.exists(layer_path):
                layer_path = expected_path
            if not layer_path or not os.path.exists(layer_path):
                continue

            normalized = os.path.abspath(layer_path)
            if normalized in seen:
                continue
            seen.add(normalized)

            layers_to_crop.append({
                "input": layer_path,
                "crop": self._raster_variant_path(layer_path, "_crop"),
                "masked": self._raster_variant_path(layer_path, "_masked"),
                "name": layer_name,
                "pad": pad_value,
                "watershed_masked": watershed_masked,
            })
        return layers_to_crop

    def _collect_layers_to_mask(self, geometry_folder):
        """Return cropped raster layers that should be watershed-masked."""
        layers_to_mask = []
        for layer_info in self._collect_layers_to_crop(geometry_folder):
            crop_path = layer_info["crop"]
            if not crop_path or not os.path.exists(crop_path):
                continue
            layers_to_mask.append({
                "input": crop_path,
                "output": layer_info["masked"],
                "name": f"{layer_info['name']}_Masked",
                "pad": layer_info.get("pad"),
                "watershed_masked": layer_info.get("watershed_masked", True),
            })
        return layers_to_mask

    def _raster_variant_path(self, raster_path, suffix):
        """Return a sidecar raster path with the requested suffix."""
        folder, filename = os.path.split(raster_path)
        stem, ext = os.path.splitext(filename)
        for existing_suffix in ("_crop", "_masked"):
            if stem.endswith(existing_suffix):
                stem = stem[: -len(existing_suffix)]
        return os.path.join(folder, f"{stem}{suffix}{ext or '.tif'}")

    def _remove_existing_raster(self, raster_path):
        """Remove an old raster output before overwriting it."""
        if not raster_path or not os.path.exists(raster_path):
            return
        try:
            os.remove(raster_path)
        except Exception as e:
            self.log_message(f"WARNING: Could not remove existing raster {raster_path}: {e}")

    def _preferred_display_raster_path(self, raster_path):
        """Return masked/cropped raster variant if it exists for inspection."""
        if not raster_path:
            return raster_path
        stem = os.path.splitext(os.path.basename(raster_path))[0]
        if stem.endswith("_masked") or stem.endswith("_crop"):
            return raster_path
        masked_path = self._raster_variant_path(raster_path, "_masked")
        if os.path.exists(masked_path):
            return masked_path
        crop_path = self._raster_variant_path(raster_path, "_crop")
        if os.path.exists(crop_path):
            return crop_path
        return raster_path

    def _preferred_display_layer_name(self, layer_name, raster_path):
        """Return a layer name matching the preferred display raster variant."""
        stem = os.path.splitext(os.path.basename(raster_path or ""))[0]
        if stem.endswith("_masked") and not layer_name.endswith("_Masked"):
            return f"{layer_name}_Masked"
        if stem.endswith("_crop") and not layer_name.endswith("_Crop"):
            return f"{layer_name}_Crop"
        return layer_name
