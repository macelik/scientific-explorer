"""Thin adapter around the supplied pyEpiAneufinder package.

The package's ``__init__`` eagerly imports plotting/postprocessing modules
(matplotlib, seaborn, ...). The app only needs the scientific functions, so
we register the package directory as a namespace-style package *without*
executing ``__init__.py`` and import the three scientific modules directly.
The functions themselves are used unmodified.
"""
from __future__ import annotations

import os
import sys
import types


def _register_package(pkg_dir: str) -> None:
    pkg_dir = os.path.abspath(pkg_dir)
    if not os.path.isfile(os.path.join(pkg_dir, "distance_statistics.py")):
        raise FileNotFoundError(f"pyEpiAneufinder package not found at {pkg_dir}")
    existing = sys.modules.get("pyEpiAneufinder")
    if existing is not None and getattr(existing, "__path__", None) == [pkg_dir]:
        return
    pkg = types.ModuleType("pyEpiAneufinder")
    pkg.__path__ = [pkg_dir]  # type: ignore[attr-defined]
    pkg.__file__ = os.path.join(pkg_dir, "__init__.py")
    sys.modules["pyEpiAneufinder"] = pkg


_loaded = False


def load(pkg_dir: str):
    """Register the package and return the scientific functions."""
    global _loaded, dist_ad, seq_dist_ad, permutation_test_ad, global_permutation_test_ad
    global trimmed_mean_iqr, assign_gainloss_new, weighted_scale_search
    _register_package(pkg_dir)
    from pyEpiAneufinder.distance_statistics import dist_ad as _dist_ad, seq_dist_ad as _seq_dist_ad
    from pyEpiAneufinder.get_breakpoints import (
        permutation_test_ad as _perm,
        global_permutation_test_ad as _gperm,
    )
    from pyEpiAneufinder.assign_somy import (
        trimmed_mean_iqr as _tmi,
        assign_gainloss_new as _agn,
        weighted_scale_search as _wss,
    )
    dist_ad = _dist_ad
    seq_dist_ad = _seq_dist_ad
    permutation_test_ad = _perm
    global_permutation_test_ad = _gperm
    trimmed_mean_iqr = _tmi
    assign_gainloss_new = _agn
    weighted_scale_search = _wss
    _loaded = True
    return types.SimpleNamespace(
        dist_ad=dist_ad,
        seq_dist_ad=seq_dist_ad,
        permutation_test_ad=permutation_test_ad,
        global_permutation_test_ad=global_permutation_test_ad,
        trimmed_mean_iqr=trimmed_mean_iqr,
        assign_gainloss_new=assign_gainloss_new,
        weighted_scale_search=weighted_scale_search,
    )
