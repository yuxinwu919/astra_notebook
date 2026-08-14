"""Custom colormaps for accelerator physics visualization.

The canonical 348-color SLAC-DESY beam density colormap is provided by
the embedded ``_plotting`` core. This module re-exports it along with
beamscope-specific colormap utilities.
"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from .._plotting.cosmetics import SLAC_DESY_CMAP

BEAM_COLORMAP = SLAC_DESY_CMAP
"""Default beam density colormap (SLAC-DESY, 348 colors)."""


def get_beam_colormap() -> LinearSegmentedColormap:
    """Return the default beam density colormap."""
    return BEAM_COLORMAP


def make_transparent_cmap(
    base_cmap: LinearSegmentedColormap,
    alpha_min: float = 0.0,
    alpha_max: float = 1.0,
) -> LinearSegmentedColormap:
    """Create a version of a colormap with variable alpha.

    Low values → transparent, high values → opaque.
    Useful for overlaying density plots.

    Args:
        base_cmap: Base colormap.
        alpha_min: Alpha at the low end.
        alpha_max: Alpha at the high end.

    Returns:
        New colormap with modified alpha.
    """
    base_cmap = base_cmap.copy()
    n = base_cmap.N
    alphas = np.linspace(alpha_min, alpha_max, n)
    base_cmap._init()
    base_cmap._lut[:, -1] = alphas
    return base_cmap
