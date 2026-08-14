"""Plotting functions for phase space, distributions, comparison, and dashboards.

All plot functions follow the convention:
  - Accept a Distribution (or dict of Distributions) as first argument
  - Accept an optional ``ax`` or ``axes`` parameter for subplot embedding
  - Return a matplotlib Figure
  - Never call plt.show() internally

Modules:
  - _precompute: shared pre-computation layer (unit conversions)
  - _artists: reusable matplotlib artist elements
  - overview: 6D phase space overview (3×2 grid)
  - phase_space: single-plane phase space plots
  - detail: detailed single-dimension view with marginals
  - distributions: 1D projection histograms
  - comparison: multi-distribution comparison
  - dashboard: comprehensive multi-panel dashboard
  - slice_plots: slice analysis visualization
  - bff_plots: bunch form factor visualization
  - emit_plots: ASTRA emit file evolution plots (delegates to _plotting core)
"""

from .phase_space import plot_phase_space, plot_transverse_phase_space
from .distributions import plot_distributions
from .comparison import plot_comparison
from .dashboard import plot_dashboard
from .overview import plot_overview
from .detail import plot_detail
from .slice_plots import plot_slice_dashboard, plot_current_profile, plot_energy_chirp, plot_slice_emittance
from .bff_plots import plot_bff, plot_bff_with_amplitude
from .emit_plots import (
    plot_envelope_evolution, plot_emittance_evolution, plot_energy_evolution,
    plot_eigen_emittances, plot_ref_trajectory, plot_emit_dashboard,
    plot_transverse_size, plot_size_with_lattice, plot_size_with_magnets,
)

__all__ = [
    "plot_phase_space",
    "plot_transverse_phase_space",
    "plot_distributions",
    "plot_comparison",
    "plot_dashboard",
    "plot_overview",
    "plot_detail",
    "plot_slice_emittance",
    "plot_current_profile",
    "plot_energy_chirp",
    "plot_slice_dashboard",
    "plot_bff",
    "plot_bff_with_amplitude",
    "plot_envelope_evolution",
    "plot_emittance_evolution",
    "plot_energy_evolution",
    "plot_eigen_emittances",
    "plot_ref_trajectory",
    "plot_emit_dashboard",
    "plot_transverse_size",
    "plot_size_with_lattice",
    "plot_size_with_magnets",
]
