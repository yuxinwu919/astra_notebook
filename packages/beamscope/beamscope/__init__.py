"""beamscope — Accelerator particle distribution visualization toolkit.

.. code-block:: python

    import beamscope
    from beamscope.io import read_distribution

    # Load a distribution
    dist = read_distribution("bunch.ini")

    # Quick stats
    stats = beamscope.compute_statistics(dist)
    beamscope.print_statistics(stats)

    # Plot phase space
    fig = beamscope.plot_phase_space(dist, plane="x")
    fig.savefig("phase_xy.png")
"""

__version__ = "0.2.0"

from .distribution import Distribution

# Lazy imports for heavy subpackages
__all__ = [
    "Distribution",
    "read_distribution",
    "compute_statistics",
    "print_statistics",
    "compute_slice_analysis",
    "compute_geometric_emittance",
    "compute_normalized_emittance",
    "compute_twiss_parameters",
    "compute_emittance_ellipse_params",
    "compute_bff",
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
    # Emit evolution plots (powered by _plotting core)
    "plot_envelope_evolution",
    "plot_emittance_evolution",
    "plot_energy_evolution",
    "plot_eigen_emittances",
    "plot_ref_trajectory",
    "plot_emit_dashboard",
    "plot_transverse_size",
    "plot_size_with_lattice",
    "plot_size_with_magnets",
    "set_publication_style",
]


def __getattr__(name: str):
    """Lazy import to avoid loading matplotlib until needed."""
    if name == "read_distribution":
        from .io import read_distribution as _fn
        return _fn
    if name == "compute_statistics":
        from .analysis.statistics import compute_statistics as _fn
        return _fn
    if name == "print_statistics":
        from .analysis.statistics import print_statistics as _fn
        return _fn
    if name == "compute_slice_analysis":
        from .analysis.slices import compute_slice_analysis as _fn
        return _fn
    if name == "compute_geometric_emittance":
        from .analysis.emittance import compute_geometric_emittance as _fn
        return _fn
    if name == "compute_normalized_emittance":
        from .analysis.emittance import compute_normalized_emittance as _fn
        return _fn
    if name == "compute_twiss_parameters":
        from .analysis.emittance import compute_twiss_parameters as _fn
        return _fn
    if name == "compute_emittance_ellipse_params":
        from .analysis.emittance import compute_emittance_ellipse_params as _fn
        return _fn
    if name == "compute_bff":
        from .analysis.bff import compute_bff as _fn
        return _fn
    if name == "plot_phase_space":
        from .plot.phase_space import plot_phase_space as _fn
        return _fn
    if name == "plot_transverse_phase_space":
        from .plot.phase_space import plot_transverse_phase_space as _fn
        return _fn
    if name == "plot_distributions":
        from .plot.distributions import plot_distributions as _fn
        return _fn
    if name == "plot_comparison":
        from .plot.comparison import plot_comparison as _fn
        return _fn
    if name == "plot_dashboard":
        from .plot.dashboard import plot_dashboard as _fn
        return _fn
    if name == "plot_overview":
        from .plot.overview import plot_overview as _fn
        return _fn
    if name == "plot_detail":
        from .plot.detail import plot_detail as _fn
        return _fn
    if name in ("plot_slice_emittance", "plot_current_profile",
                "plot_energy_chirp", "plot_slice_dashboard"):
        import importlib
        mod = importlib.import_module(".plot.slice_plots", package="beamscope")
        return getattr(mod, name)
    if name in ("plot_bff", "plot_bff_with_amplitude"):
        import importlib
        mod = importlib.import_module(".plot.bff_plots", package="beamscope")
        return getattr(mod, name)
    if name in ("plot_envelope_evolution", "plot_emittance_evolution",
                "plot_energy_evolution", "plot_eigen_emittances",
                "plot_ref_trajectory", "plot_emit_dashboard",
                "plot_transverse_size", "plot_size_with_lattice",
                "plot_size_with_magnets"):
        import importlib
        mod = importlib.import_module(".plot.emit_plots", package="beamscope")
        return getattr(mod, name)
    if name == "set_publication_style":
        from .style.rcparams import set_publication_style as _fn
        return _fn
    raise AttributeError(f"module 'beamscope' has no attribute '{name}'")
