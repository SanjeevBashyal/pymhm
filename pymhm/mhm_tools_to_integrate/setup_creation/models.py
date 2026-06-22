# -*- coding: utf-8 -*-
"""Plain data model for regions, domains, watersheds, and outlets."""
from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .latlon import create_latlon_file


Bounds = tuple[float, float, float, float]
Header = dict[str, float | int]
LogCallback = Callable[[str], None]


def _merge_bounds(left: Bounds | None, right: Bounds | None) -> Bounds | None:
    """Return a bounding box covering both inputs."""
    if left is None:
        return right
    if right is None:
        return left
    return (
        min(left[0], right[0]),
        max(left[1], right[1]),
        min(left[2], right[2]),
        max(left[3], right[3]),
    )


@dataclass
class Outlet:
    """A model outlet point independent of GIS layer/UI objects."""

    outlet_id: str
    x: float
    y: float
    crs: str | None = None
    is_domain: bool = False
    attributes: dict[str, Any] = field(default_factory=dict)
    snapped_x: float | None = None
    snapped_y: float | None = None

    @property
    def snapped(self) -> bool:
        """Return True when snap coordinates are available."""
        return self.snapped_x is not None and self.snapped_y is not None

    def coordinate(self, snapped: bool = True) -> tuple[float, float]:
        """Return snapped or original coordinate."""
        if snapped and self.snapped_x is not None and self.snapped_y is not None:
            return self.snapped_x, self.snapped_y
        return self.x, self.y


@dataclass
class Watershed:
    """Watershed data and per-watershed processing inputs."""

    watershed_id: str
    outlet: Outlet
    bounds: Bounds | None = None
    raster_path: Path | None = None
    vector_path: Path | None = None
    elevation_band_paths: list[Path] = field(default_factory=list)
    children: list["Watershed"] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def add_child(self, watershed: "Watershed") -> None:
        """Attach a nested watershed to this watershed."""
        self.children.append(watershed)
        self.bounds = _merge_bounds(self.bounds, watershed.bounds)

    def snap_outlet(
        self, snapper: Callable[[Outlet], Outlet], log: LogCallback | None = None
    ) -> Outlet:
        """Snap this watershed outlet using a caller-provided algorithm."""
        if log:
            log(f"Snapping outlet {self.outlet.outlet_id}.")
        self.outlet = snapper(self.outlet)
        return self.outlet

    def create_elevation_bands(
        self,
        algorithm: Callable[["Watershed"], Iterable[Path]],
        log: LogCallback | None = None,
    ) -> list[Path]:
        """Create elevation bands using a caller-provided algorithm."""
        if log:
            log(f"Creating elevation bands for watershed {self.watershed_id}.")
        self.elevation_band_paths = list(algorithm(self))
        return self.elevation_band_paths


@dataclass
class Domain(Watershed):
    """Domain watershed containing all nested watersheds draining inside it."""

    watersheds: list[Watershed] = field(default_factory=list)

    def add_watershed(self, watershed: Watershed) -> None:
        """Attach a watershed to this domain."""
        self.watersheds.append(watershed)
        self.bounds = _merge_bounds(self.bounds, watershed.bounds)

    @classmethod
    def from_domain_outlet(
        cls, outlet: Outlet, watershed_id: str | None = None, **kwargs: Any
    ) -> "Domain":
        """Create a domain from an outlet marked as a domain outlet."""
        outlet.is_domain = True
        return cls(
            watershed_id=watershed_id or outlet.outlet_id,
            outlet=outlet,
            **kwargs,
        )


@dataclass
class Region:
    """A region containing one or more independent domains."""

    region_id: str
    dem_path: Path | None = None
    crs: str | None = None
    bounds: Bounds | None = None
    domains: list[Domain] = field(default_factory=list)
    outputs: dict[str, Path] = field(default_factory=dict)
    grid_headers: dict[str, Header] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def add_domain(self, domain: Domain) -> None:
        """Attach a domain and expand region bounds."""
        self.domains.append(domain)
        self.bounds = _merge_bounds(self.bounds, domain.bounds)

    def add_domains_from_outlets(
        self,
        outlets: Iterable[Outlet],
        watershed_factory: Callable[[Outlet], Watershed | Domain],
    ) -> None:
        """Build domain/watershed hierarchy from outlets using plain data."""
        domain_by_id: dict[str, Domain] = {}
        watersheds: list[Watershed] = []
        for outlet in outlets:
            watershed = watershed_factory(outlet)
            if outlet.is_domain:
                domain = (
                    watershed
                    if isinstance(watershed, Domain)
                    else Domain.from_domain_outlet(outlet, watershed.watershed_id)
                )
                self.add_domain(domain)
                domain_by_id[outlet.outlet_id] = domain
            else:
                watersheds.append(watershed)

        for watershed in watersheds:
            for domain in self.domains:
                if self._watershed_belongs_to_domain(watershed, domain):
                    domain.add_watershed(watershed)
                    break

    def create_terrain_products(
        self,
        algorithm: Callable[["Region"], dict[str, Path]],
        log: LogCallback | None = None,
    ) -> dict[str, Path]:
        """Create slope, aspect, flow direction, and flow accumulation."""
        if log:
            log(f"Creating terrain products for region {self.region_id}.")
        products = algorithm(self)
        self.outputs.update(products)
        return products

    def create_latlon(
        self,
        out_file: Path,
        level0: Header | Path,
        level1: Header | float | Path,
        level11: Header | float | Path | None = None,
        level2: Header | float | Path | None = None,
        log: LogCallback | None = None,
    ) -> Path:
        """Create latlon.nc for this region through mhm_tools."""
        result = create_latlon_file(
            out_file=out_file,
            level0=level0,
            level1=level1,
            level11=level11,
            level2=level2,
            crs=self.crs,
            log=log,
        )
        self.outputs["latlon"] = result
        return result

    def _watershed_belongs_to_domain(
        self, watershed: Watershed, domain: Domain
    ) -> bool:
        """Return True when watershed outlet lies inside domain bounds."""
        if domain.bounds is None:
            return False
        x, y = watershed.outlet.coordinate()
        xmin, xmax, ymin, ymax = domain.bounds
        return xmin <= x <= xmax and ymin <= y <= ymax
