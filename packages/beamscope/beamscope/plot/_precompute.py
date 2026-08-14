"""Shared pre-computation layer for all plot modules.

Provides a single, consistent function that maps particle coordinates
to display-friendly units (mm, mrad, MeV/c, ps, nC, %).

All plot functions should import `precompute` and `clip_percentile`
from this module rather than duplicating the logic.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from beamscope.distribution import Distribution

# ---------------------------------------------------------------------------
# Variable definitions — single source of truth for axis labels and units
# ---------------------------------------------------------------------------

_VARIABLE_DEFS: dict[str, tuple[str, str]] = {
    "x":       ("x [mm]",       "x"),
    "y":       ("y [mm]",       "y"),
    "z":       ("z [mm]",       "z"),
    "xp":      ("x' [mrad]",    "xp"),
    "yp":      ("y' [mrad]",    "yp"),
    "dp":      ("δp/p [%]",     "dp"),
    "E":       ("E [MeV]",      "E"),
    "clock":   ("clock [ps]",   "clock"),
    "r":       ("r [mm]",       "r"),
    "rp":      ("r' [mrad]",    "rp"),
    "charge":  ("Q [nC]",       "charge"),
}

# 3×2 grid layout for overview plot
_OVERVIEW_PANELS: list[tuple[int, int, str, str, str, str, str]] = [
    (0, 0, "x",  "xp",  "x–x'",    "x [mm]",    "x' [mrad]"),
    (1, 0, "y",  "yp",  "y–y'",    "y [mm]",    "y' [mrad]"),
    (2, 0, "z",  "dp",  "z–δp/p",  "z [mm]",    "δp/p [%]"),
    (0, 1, "x",  "y",   "x–y",     "x [mm]",    "y [mm]"),
    (1, 1, "z",  "x",   "z–x",     "z [mm]",    "x [mm]"),
    (2, 1, "z",  "y",   "z–y",     "z [mm]",    "y [mm]"),
]


# ---------------------------------------------------------------------------
# Core precompute function
# ---------------------------------------------------------------------------

def precompute(dist: Distribution) -> dict[str, np.ndarray]:
    """Pre-compute all display-friendly quantities for a Distribution.

    All quantities are centroid-subtracted and unit-converted:
      - positions: m → mm
      - angles: rad → mrad
      - energy: eV → MeV
      - time: s → ps
      - charge: nC (native)
      - δp/p: %

    Returns a dict mapping variable key → numpy array.
    Only active particles are included.
    """
    mask = dist.active
    x = dist.x[mask]
    y = dist.y[mask]
    z = dist.z[mask]
    pz = dist.pz[mask]
    px = dist.px[mask]
    py = dist.py[mask]
    clock = dist.clock[mask]
    charge = dist.charge[mask]

    pz_abs = np.abs(pz)

    mean_px = float(np.mean(px))
    mean_py = float(np.mean(py))
    mean_pz = float(np.mean(pz))
    mean_x = float(np.mean(x))
    mean_y = float(np.mean(y))
    mean_z = float(np.mean(z))

    return {
        "x":      (x - mean_x) * 1e3,
        "y":      (y - mean_y) * 1e3,
        "z":      (z - mean_z) * 1e3,
        "xp":     (px - mean_px) / pz_abs * 1e3,
        "yp":     (py - mean_py) / pz_abs * 1e3,
        "dp":     (pz - mean_pz) / mean_pz * 100,
        "E":      pz * 1e-6,
        "clock":  clock * 1e12,
        "r":      np.sqrt((x - mean_x)**2 + (y - mean_y)**2) * 1e3,
        "rp":     np.sqrt(
            ((px - mean_px) / pz_abs)**2 + ((py - mean_py) / pz_abs)**2
        ) * 1e3,
        "charge": charge,
    }


# ---------------------------------------------------------------------------
# Axis range helpers
# ---------------------------------------------------------------------------

def clip_percentile(
    data: np.ndarray,
    lo: float = 2.0,
    hi: float = 98.0,
    n_min: int = 50,
) -> tuple[float, float]:
    """Return (vmin, vmax) bounds at given percentiles for robust axis ranges.

    Default lo=2, hi=98 clips ~4% of particles to focus on the core
    distribution. For datasets with N < n_min, falls back to min/max.

    Args:
        data: 1D data array.
        lo: Lower percentile (default 2%).
        hi: Upper percentile (default 98%).
        n_min: Minimum data size for percentile clipping.

    Returns:
        (vmin, vmax) tuple.
    """
    if len(data) == 0:
        return 0.0, 1.0
    if len(data) < n_min:
        return float(np.min(data)), float(np.max(data))
    return float(np.percentile(data, lo)), float(np.percentile(data, hi))


def get_variable_label(key: str) -> str:
    """Get display label for a variable key."""
    info = _VARIABLE_DEFS.get(key)
    return info[0] if info else key


def get_overview_panels() -> list[tuple[int, int, str, str, str, str, str]]:
    """Get the 3×2 overview panel definitions."""
    return _OVERVIEW_PANELS


def get_variable_defs() -> dict[str, tuple[str, str]]:
    """Get the full variable definitions dict."""
    return dict(_VARIABLE_DEFS)
