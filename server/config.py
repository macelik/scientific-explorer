"""Workspace/data-root resolution. Nothing here hardcodes historical paths."""
from __future__ import annotations

import os
from dataclasses import dataclass


HERE = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(HERE)                 # interactive-explorer/
DEFAULT_ROOT = os.path.dirname(APP_DIR)         # the planning workspace root


@dataclass
class Config:
    data_root: str
    h5ad: str
    cohort_dir: str
    package_dir: str
    production_breakpoints: str
    clusters_json: str
    cache_dir: str
    sessions_dir: str
    integration_dir: str = ""
    h5mu: str = ""
    prototype_dir: str = ""
    derived_dir: str = ""
    cluster_breakpoints: str = ""

    @staticmethod
    def from_env(data_root: str | None = None) -> "Config":
        root = os.path.abspath(data_root or os.environ.get("PYEPI_SCIENTIFIC_ROOT") or DEFAULT_ROOT)
        cohort = os.environ.get("PYEPI_SCIENTIFIC_COHORT", os.path.join(root, "default_300cells"))
        h5ad = os.environ.get("PYEPI_SCIENTIFIC_H5AD", os.path.join(root, "pyepi_results", "count_matrix.h5ad"))
        pkg = os.environ.get("PYEPI_SCIENTIFIC_PACKAGE", os.path.join(root, "pyEpiAneufinder"))
        cache = os.environ.get("PYEPI_SCIENTIFIC_CACHE", os.path.join(APP_DIR, "cache"))
        sessions = os.environ.get("PYEPI_SCIENTIFIC_SESSIONS", os.path.join(APP_DIR, "sessions"))
        integration_dir = os.environ.get("PYEPI_SCIENTIFIC_INTEGRATION", os.path.join(root, "integration-story"))
        prototype = os.path.join(integration_dir, "data", "prototype")
        reference = os.path.join(root, "scientific-explorer", "reference", "integration-prototype", "regenerated")
        if not os.path.isdir(prototype) and os.path.isfile(os.path.join(reference, "verification.json")):
            prototype = reference
        return Config(
            data_root=root,
            h5ad=h5ad,
            cohort_dir=cohort,
            package_dir=pkg,
            production_breakpoints=os.path.join(root, "pyepi_results", "breakpoints.csv"),
            clusters_json=os.path.join(root, "pyepi_results", "clusters.json"),
            cache_dir=cache,
            sessions_dir=sessions,
            integration_dir=integration_dir,
            h5mu=os.environ.get("PYEPI_SCIENTIFIC_H5MU", ""),
            prototype_dir=os.environ.get("PYEPI_SCIENTIFIC_PROTOTYPE", prototype),
            derived_dir=os.environ.get("PYEPI_SCIENTIFIC_DERIVED", os.path.join(integration_dir, "data", "derived")),
            cluster_breakpoints=os.environ.get("PYEPI_SCIENTIFIC_CLUSTER_BREAKPOINTS", os.path.join(integration_dir, "data", "raw_accepted_breakpoints.tsv")),
        )

    def validate(self) -> None:
        for label, p in [("h5ad", self.h5ad), ("cohort dir", self.cohort_dir), ("package", self.package_dir)]:
            if not os.path.exists(p):
                raise FileNotFoundError(f"{label} not found: {p}")
        os.makedirs(self.cache_dir, exist_ok=True)
        os.makedirs(self.sessions_dir, exist_ok=True)
