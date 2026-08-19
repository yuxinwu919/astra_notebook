"""Plotting layer: postpro / lineplot / fieldplot replacements.

All plots: plain 2D scatter with deterministic subsampling, 0.5-99.5 percentile
range clipping (outlier-proof), SI input / display-unit output with
complete legends. Call plot.style.set_style() once before plotting.
"""

from .style import set_style, reset_style
from .phase_space import plot_phase_space, plot_transverse_phase_space
from .distributions import plot_distributions, plot_energy_distribution
from .overview import plot_overview, plot_transverse_profile
from .emit_plots import (
    plot_envelope_evolution,
    plot_correlated_energy_spread,
    plot_ref_momentum,
    plot_divergence_evolution,
    plot_emittance_evolution,
    plot_energy_evolution,
    plot_bunch_length_evolution,
    plot_energy_spread_evolution,
    plot_ref_trajectory,
    plot_velocity_evolution,
    plot_step_size_evolution,
    plot_eigen_emittances,
    plot_emit_dashboard,
    plot_lineplot_overview,
)
from .slice_plots import (
    plot_current_profile,
    plot_slice_emittance,
    plot_energy_chirp,
    plot_slice_dashboard,
)
from .bff_plots import plot_bff, plot_bff_with_amplitude
from .field_plots import (plot_cavity_field, plot_solenoid_field,
                          plot_solenoid_components, plot_te_field,
                          plot_field_expansion_radius,
                           plot_3d_field_map, plot_3d_field_slices,
                           resolve_3d_plane, plane_fixed_axis)

__all__ = [
    "set_style", "reset_style",
    "plot_phase_space", "plot_transverse_phase_space",
    "plot_distributions", "plot_energy_distribution",
    "plot_overview", "plot_transverse_profile",
    "plot_envelope_evolution", "plot_divergence_evolution",
    "plot_emittance_evolution", "plot_energy_evolution",
    "plot_bunch_length_evolution", "plot_energy_spread_evolution",
    "plot_ref_trajectory", "plot_velocity_evolution", "plot_step_size_evolution",
    "plot_correlated_energy_spread", "plot_ref_momentum",
    "plot_eigen_emittances",
    "plot_emit_dashboard", "plot_lineplot_overview",
    "plot_current_profile", "plot_slice_emittance",
    "plot_energy_chirp", "plot_slice_dashboard",
    "plot_bff", "plot_bff_with_amplitude",
    "plot_cavity_field", "plot_solenoid_field",
    "plot_solenoid_components", "plot_te_field", "plot_field_expansion_radius",
     "plot_3d_field_map", "plot_3d_field_slices",
    "resolve_3d_plane", "plane_fixed_axis",
    # advanced_plots (postpro 5.6 / lineplot 菜单2/3/4 / fieldplot 扩展)
    "plot_losses", "plot_beam_loading", "plot_beta_alpha",
    "plot_phase_advance", "plot_coherence_length", "plot_phase_scan",
    "plot_scan_fom", "plot_error_hist", "plot_reduced_emittance",
    "plot_trace_emittance", "plot_core_emittance", "plot_larmor",
    "plot_probe_trajectories", "plot_space_charge_fields",
    "plot_cathode_emission", "plot_slice_mismatch",
    "slice_mismatch",
    "plot_pscan_dedz", "plot_pscan_compression", "plot_pscan_compression_time",
    "plot_scan_position", "plot_tcheck_scaling", "plot_tcheck_counter",
    "plot_cr_emit",
    "plot_z_plot", "plot_field_profile", "plot_curved_cathode_contour",
    "plot_core_brightness", "plot_slice_ellipses_3d",
    # Batch C: 手册第5章末梢 (t 轴变体在 emit_plots 的 x_axis 参数)
    "plot_envelope_with_aperture", "aperture_elements",
    "plot_laser_on_axis", "plot_plasma_profile",
    "plot_core_fraction_curves",
]

from .advanced_plots import (
    plot_losses, plot_beam_loading, plot_beta_alpha, plot_phase_advance,
    plot_coherence_length, plot_phase_scan, plot_scan_fom, plot_error_hist,
    plot_reduced_emittance, plot_trace_emittance, plot_core_emittance,
    plot_larmor, plot_probe_trajectories, plot_space_charge_fields,
    plot_cathode_emission, plot_slice_mismatch,
    slice_mismatch,
)

from .advanced_plots import (
    plot_pscan_dedz, plot_pscan_compression, plot_pscan_compression_time,
    plot_scan_position, plot_tcheck_scaling, plot_tcheck_counter,
    plot_cr_emit,
    plot_z_plot, plot_field_profile, plot_curved_cathode_contour,
    plot_core_brightness, plot_slice_ellipses_3d,
    plot_envelope_with_aperture, aperture_elements,
    plot_laser_on_axis, plot_plasma_profile, plot_core_fraction_curves,
)
