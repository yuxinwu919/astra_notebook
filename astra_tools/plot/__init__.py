"""Plotting layer: postpro / lineplot / fieldplot replacements.

All plots: unified 2D Gaussian-kernel KDE density, 0.5-99.5 percentile
range clipping (outlier-proof), SI input / display-unit output with
complete legends. Call plot.style.set_style() once before plotting.
"""

from .style import set_style, reset_style
from .phase_space import plot_phase_space, plot_transverse_phase_space
from .distributions import plot_distributions, plot_energy_distribution
from .overview import plot_overview, plot_transverse_profile
from .emit_plots import (
    plot_envelope_evolution,
    plot_emittance_evolution,
    plot_energy_evolution,
    plot_ref_trajectory,
    plot_eigen_emittances,
    plot_emit_dashboard,
)
from .slice_plots import (
    plot_current_profile,
    plot_slice_emittance,
    plot_energy_chirp,
    plot_slice_dashboard,
)
from .bff_plots import plot_bff, plot_bff_with_amplitude
from .field_plots import plot_cavity_field, plot_solenoid_field

__all__ = [
    "set_style", "reset_style",
    "plot_phase_space", "plot_transverse_phase_space",
    "plot_distributions", "plot_energy_distribution",
    "plot_overview", "plot_transverse_profile",
    "plot_envelope_evolution", "plot_emittance_evolution",
    "plot_energy_evolution", "plot_ref_trajectory",
    "plot_eigen_emittances", "plot_emit_dashboard",
    "plot_current_profile", "plot_slice_emittance",
    "plot_energy_chirp", "plot_slice_dashboard",
    "plot_bff", "plot_bff_with_amplitude",
    "plot_cavity_field", "plot_solenoid_field",
]

from .advanced_plots import (
    plot_losses, plot_beam_loading, plot_beta_alpha, plot_phase_advance,
    plot_coherence_length, plot_phase_scan, plot_scan_fom, plot_error_hist,
    plot_reduced_emittance, plot_trace_emittance, plot_core_emittance,
    plot_larmor, plot_probe_trajectories, plot_space_charge_fields,
    plot_cathode_emission, plot_slice_mismatch, plot_3d_map_slices,
    slice_mismatch,
)

from .advanced_plots import (
    plot_pscan_dedz, plot_pscan_compression, plot_tcheck_scaling,
    plot_z_plot, plot_field_profile, plot_curved_cathode_contour,
    plot_core_brightness, plot_slice_ellipses_3d,
)
