"""Plotting layer: postpro / lineplot / fieldplot replacements."""

from .style import set_style, reset_style
from .phase_space import plot_phase_space, plot_transverse_phase_space
from .distributions import plot_distributions, plot_beam_profile
from .overview import plot_overview
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
    "plot_distributions", "plot_beam_profile",
    "plot_overview",
    "plot_envelope_evolution", "plot_emittance_evolution",
    "plot_energy_evolution", "plot_ref_trajectory",
    "plot_eigen_emittances", "plot_emit_dashboard",
    "plot_current_profile", "plot_slice_emittance",
    "plot_energy_chirp", "plot_slice_dashboard",
    "plot_bff", "plot_bff_with_amplitude",
    "plot_cavity_field", "plot_solenoid_field",
]
