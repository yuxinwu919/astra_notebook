"""2D density estimation for phase-space plots.

Unified Gaussian-kernel KDE computed on a grid: histogram + Gaussian
smoothing (mathematically the density estimate of a Gaussian kernel at
fixed bandwidth, visually identical to scipy gaussian_kde contour
plots). Benchmark: 1e6 particles ~0.1 s, 1e7 ~1 s — no subsampling
needed for typical ASTRA bunches.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def density2d(
    x: np.ndarray,
    y: np.ndarray,
    bins: int = 160,
    range_xy=None,
    weights=None,
    bandwidth_pix: float = 3.0,
):
    """2D Gaussian-kernel density estimate on a grid.

    Args:
        x, y: coordinates.
        bins: grid resolution (square).
        range_xy: ((xmin, xmax), (ymin, ymax)); default = 0.5-99.5
            percentile clip of the data (outliers excluded from the
            display range; see clip_percentile).
        weights: optional weights (macro-particle charge).
        bandwidth_pix: Gaussian kernel sigma in pixels (KDE bandwidth).

    Returns:
        (z, xe, ye): density grid and bin edges.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        raise ValueError("need at least 3 points for a density plot")
    if range_xy is None:
        range_xy = (clip_percentile(x), clip_percentile(y))

    h, xe, ye = np.histogram2d(x, y, bins=bins, range=range_xy, weights=weights)
    # Gaussian-kernel smoothing (KDE on the grid)
    if bandwidth_pix > 0:
        z = ndimage.gaussian_filter(h, sigma=bandwidth_pix)
    else:
        z = h
    # Normalize to probability density
    dx = xe[1] - xe[0]
    dy = ye[1] - ye[0]
    total = float(np.sum(z) * dx * dy)
    if total > 0:
        z = z / total
    return z, xe, ye


def clip_percentile(a, q: float = 0.5):
    """Robust display range: [q, 100-q] percentile.

    Prevents a few outliers from collapsing the bulk of the data into a
    corner (the classic 'everything clumped' plot). Returns (lo, hi);
    the fraction of points outside is available for annotation.
    """
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    if len(a) == 0:
        return 0.0, 1.0
    lo, hi = np.percentile(a, [q, 100.0 - q])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < 1e-30:
        pad = max(abs(hi), abs(lo), 1.0) * 1e-6 + 1e-30
        lo, hi = lo - pad, hi + pad
    return float(lo), float(hi)


def outside_fraction(a, lo, hi):
    """Fraction of points outside [lo, hi] (for annotation)."""
    a = np.asarray(a)
    return float(np.mean((a < lo) | (a > hi)))
