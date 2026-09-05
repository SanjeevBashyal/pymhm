"""Path-only Show/Save commands; reuse outlet masks independently of domain IDs."""
from pathlib import Path

from ...handlers.raster import tasks as raster
from ...handlers.state.domain_state import cached_delineation, delineation_fingerprint


def _outlet(project_folder, filled_dem, outlet, context, task):
    raster._cancelled(task)
    x, y = outlet["x"], outlet["y"]
    result = cached_delineation(project_folder, filled_dem, x, y, outlet.get("cached"))
    if result is None:
        if context is None:
            context = raster._flow_context(filled_dem)
        digest = delineation_fingerprint(filled_dem, x, y)
        path = Path(outlet["raster_path"])
        # A preview must not overwrite the last committed domain mask.
        path = path.with_name(f"{path.stem}_{digest[:16]}{path.suffix}")
        result = raster._delineate(context, x, y, str(path), outlet["basin_id"])
        result.pop("mask")
        result.update(fingerprint=digest, mask_value=outlet["basin_id"])
    result["picked"] = dict(outlet["picked"])
    for key in ("confidence", "snap_count", "snap_distance_m"):
        result[key] = outlet.get(key)
    return result, context


def show(task, *, project_folder, filled_dem, outlet):
    """Delineate only on explicit request, or return a valid cached mask."""
    result, _context = _outlet(project_folder, filled_dem, outlet, None, task)
    raster._progress(task, 100)
    return result


def save(task, *, project_folder, filled_dem, outlets, merged_path, dem_domain=None):
    """Prepare every listed outlet, merging only the active domain masks."""
    context = None
    results, masks = {}, []
    for index, outlet in enumerate(outlets):
        result, context = _outlet(project_folder, filled_dem, outlet, context, task)
        results[outlet["outlet_id"]] = result
        if outlet.get("domain_id"):
            masks.append((outlet["domain_id"], result["raster_path"]))
        raster._progress(task, 80 * (index + 1) / max(1, len(outlets)))
    context = None  # Release the flow grid before reading masks for the merge.
    dem_mask = None
    if dem_domain:
        raster._cancelled(task)
        domain_id, dem_mask, _dem_path = dem_domain
        # DEM extent needs only its valid-cell mask, not a flow-grid calculation.
        reference = raster._read_raster(filled_dem, as_float=True)
        _dem, invalid = raster.pyflwdir_handler.normalise_dem(reference["array"], reference["nodata"])
        raster._dem_domain({"reference": reference, "invalid": invalid}, domain_id, dem_mask)
        masks.append((domain_id, dem_mask))
        del reference, _dem, invalid
    raster._cancelled(task)
    merged = raster.merge_domain_masks(masks, merged_path)
    raster._progress(task, 100)
    return {"outlets": results, "merged_path": merged, "dem_mask_path": dem_mask,
            "dem_domain_path": dem_domain[2] if dem_domain else None}
