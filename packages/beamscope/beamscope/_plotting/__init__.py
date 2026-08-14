"""
Embedded ASTRA plotting core — battle-tested plotting functions for beam dynamics.

This subpackage contains the core plotting and styling modules ported from the
astra_plotter library. All functions work with raw numpy arrays.

Public API
----------
Plotter:
    DensityPlot, DensityPlot_w_Hproj, DensityPlot_w_proj,
    DensityplotwProjec2x2, PlotBunchFormFactor, PlotEigenEmits,
    PlotEmit1plt, PlotEnergy1plt, PlotSize1plt, PlotSize1pltLat,
    PlotSliceParameters, PlotTransSize1plt, PlotTransSize1pltMag,
    histogram0
Cosmetics:
    SLAC_DESY_CMAP, BeamColorMap, FormatLabelSci, PrettyPlot
"""

from .plotter import (
    DensityPlot,
    DensityPlot_w_Hproj,
    DensityPlot_w_proj,
    DensityplotwProjec2x2,
    PlotBunchFormFactor,
    PlotEigenEmits,
    PlotEmit1plt,
    PlotEnergy1plt,
    PlotSize1plt,
    PlotSize1pltLat,
    PlotSliceParameters,
    PlotTransSize1plt,
    PlotTransSize1pltMag,
)
from .cosmetics import (
    BeamColorMap,
    FormatLabelSci,
    PrettyPlot,
    SLAC_DESY_CMAP,
)

__all__ = [
    # Plotter
    'DensityPlot', 'DensityPlot_w_Hproj', 'DensityPlot_w_proj',
    'DensityplotwProjec2x2', 'PlotBunchFormFactor', 'PlotEigenEmits',
    'PlotEmit1plt', 'PlotEnergy1plt', 'PlotSize1plt', 'PlotSize1pltLat',
    'PlotSliceParameters', 'PlotTransSize1plt', 'PlotTransSize1pltMag',
    # Cosmetics
    'BeamColorMap', 'FormatLabelSci', 'PrettyPlot', 'SLAC_DESY_CMAP',
]
