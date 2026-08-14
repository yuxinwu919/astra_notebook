"""astra-tools - accelerator physics pre/post-processing for ASTRA & Generator.

Jupyter-first toolkit: read ASTRA distributions and evolution files,
compute audited beam statistics, and produce publication-quality plots
replacing the official postpro / lineplot / fieldplot tools.

    import astra_tools as at
    dist = at.read_distribution("bunch.ini")
    stats = at.compute_statistics(dist)
"""

__version__ = "0.1.0"

from .distribution import Distribution

__all__ = [
    "Distribution",
    "read_distribution",
    "compute_statistics",
    "print_statistics",
]


def __getattr__(name):
    if name == "read_distribution":
        from .io import read_distribution as _fn
        return _fn
    if name == "compute_statistics":
        from .analysis.statistics import compute_statistics as _fn
        return _fn
    if name == "print_statistics":
        from .analysis.statistics import print_statistics as _fn
        return _fn
    raise AttributeError("module 'astra_tools' has no attribute " + repr(name))
