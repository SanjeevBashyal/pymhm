# -*- coding: utf-8 -*-
"""Cache the meteorology folder inspection between sessions.

`inspect_meteo_folder` opens every NetCDF in a folder twice: once for its
spatial signature and once to validate its time axis. On a multi-decade
ERA5-Land record that is over a thousand file opens per folder, which froze the
plugin for around half a minute on every project load. The result only depends
on the files themselves, so it is fingerprinted and reused.
"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from .forcing import (
    MeteoFolderSpec,
    SpatialMetadata,
    _files_for_spec,
    inspect_meteo_folder,
)
from ..state_cache import cached_payload, fingerprint, store_payload


SECTION = "meteo_inspection"


def inspection_fingerprint(spec: MeteoFolderSpec) -> tuple[str, list[Path]]:
    """Return the input fingerprint and the file list it was taken from.

    Listing the folder is cheap; opening the files is not, so the fingerprint is
    built from the listing plus each file's size and modification time.
    """
    files = _files_for_spec(spec)
    digest = fingerprint(
        files,
        {
            "kind": spec.kind,
            "source": spec.source,
            "crs": spec.crs,
            "folder": str(spec.folder),
        },
    )
    return digest, files


def _as_payload(metadata: SpatialMetadata) -> dict:
    payload = asdict(metadata)
    payload["files"] = [str(path) for path in metadata.files]
    payload["shape"] = list(metadata.shape)
    payload["bounds"] = list(metadata.bounds)
    return payload


def _from_payload(payload) -> SpatialMetadata | None:
    if not isinstance(payload, dict):
        return None
    try:
        values = dict(payload)
        values["files"] = tuple(Path(path) for path in values.get("files", ()))
        values["shape"] = tuple(values["shape"])
        values["bounds"] = tuple(values["bounds"])
        return SpatialMetadata(**values)
    except (KeyError, TypeError, ValueError):
        return None


def inspect_meteo_folder_cached(
        project_folder,
        folder,
        kind: str,
        source: str,
        crs_fallback: str | None = None,
        log=None) -> SpatialMetadata:
    """Inspect a meteorology folder, reusing an unchanged previous inspection."""
    spec = MeteoFolderSpec(kind, folder, source, crs_fallback).normalized()
    digest, files = inspection_fingerprint(spec)

    if project_folder:
        metadata = _from_payload(
            cached_payload(project_folder, SECTION, kind, digest))
        if metadata is not None and len(metadata.files) == len(files):
            if log:
                log(
                    f"{kind.title()} metadata reused for {len(files)} unchanged "
                    "NetCDF file(s)."
                )
            return metadata

    metadata = inspect_meteo_folder(folder, kind, source, crs_fallback)
    if project_folder:
        store_payload(
            project_folder, SECTION, kind, digest, _as_payload(metadata))
    return metadata


__all__ = ["inspect_meteo_folder_cached", "inspection_fingerprint"]
