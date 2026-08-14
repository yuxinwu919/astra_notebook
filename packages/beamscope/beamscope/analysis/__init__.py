"""Beam statistical analysis and phase space characterization.

All functions operate on the unified Distribution data model.
Formulae verified against:
  - ASTRA Manual V3.2, Tables 1 & 4
  - M. Reiser, "Theory and Design of Charged Particle Beams", Wiley (2008)
  - K. Floettmann, "Some basic features of the beam emittance", PRSTAB 6, 034202 (2003)

Key fixes (2026-07):
  - γ computed from momentum (not kinetic energy) per ASTRA convention
  - RMS uses sample std (ddof=1) for unbiased estimation
  - SVD-based emittance for numerical stability
  - Charge-weighted statistics support
"""

from .statistics import BeamStatistics, compute_statistics, print_statistics
from .emittance import (
    compute_geometric_emittance,
    compute_normalized_emittance,
    compute_twiss_parameters,
    compute_emittance_ellipse_params,
)
from .slices import SliceAnalysis, compute_slice_analysis
from .bff import BFFResult, compute_bff

__all__ = [
    "BeamStatistics",
    "compute_statistics",
    "print_statistics",
    "compute_geometric_emittance",
    "compute_normalized_emittance",
    "compute_twiss_parameters",
    "compute_emittance_ellipse_params",
    "SliceAnalysis",
    "compute_slice_analysis",
    "BFFResult",
    "compute_bff",
]
