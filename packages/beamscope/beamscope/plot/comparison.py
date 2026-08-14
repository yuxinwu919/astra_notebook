"""Multi-distribution comparison plots.

Compare beam properties before/after an element, or across parameter scans.
"""

from __future__ import annotations

from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.statistics import compute_statistics
from ..distribution import Distribution

ComparisonType = Literal["phase_space", "projections", "statistics"]


def plot_comparison(
    distributions: dict[str, Distribution],
    plot_type: ComparisonType = "phase_space",
    plane: str = "x",
    bins: int = 60,
    figsize: tuple[float, float] = (12, 8),
) -> plt.Figure:
    """Compare multiple distributions side by side.

    Args:
        distributions: Dict mapping labels to Distribution objects.
                       e.g. {"before": dist1, "after": dist2}
        plot_type: Type of comparison plot.
            - 'phase_space': Side-by-side phase space plots.
            - 'projections': Overlaid projection histograms.
            - 'statistics': Bar chart comparing key beam parameters.
        plane: Phase plane for 'phase_space' type ('x', 'y', 'z').
        bins: Number of bins for histograms.
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    labels = list(distributions.keys())
    n = len(labels)

    if n < 2:
        raise ValueError("Need at least 2 distributions for comparison.")

    if plot_type == "phase_space":
        from .phase_space import plot_phase_space

        fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
        if n == 1:
            axes = [axes]

        for ax, (label, dist) in zip(axes, distributions.items()):
            plot_phase_space(
                dist, plane=plane, kind="density", bins=bins,
                ax=ax, title=label, colorbar=(ax == axes[-1]),
            )

        fig.suptitle(f"Phase Space Comparison — {plane}-{plane}'", fontsize=14, fontweight="bold")
        fig.tight_layout()
        return fig

    elif plot_type == "projections":
        fig, axes = plt.subplots(1, 3, figsize=figsize)

        for ax, (dim_label, scale) in zip(
            axes,
            [("x", 1e3), ("y", 1e3), ("z", 1e3)],
        ):
            for label, dist in distributions.items():
                data = getattr(dist, dim_label)[dist.active] * scale
                ax.hist(
                    data, bins=bins, density=True,
                    alpha=0.5, label=label, edgecolor="white",
                )
            ax.set_xlabel(f"{dim_label} [mm]")
            ax.set_ylabel("Density")
            ax.set_title(f"{dim_label} Projection")
            ax.legend(fontsize=8)

        fig.suptitle("Beam Distribution Comparison", fontsize=14, fontweight="bold")
        fig.tight_layout()
        return fig

    elif plot_type == "statistics":
        stats_list = {
            label: compute_statistics(dist, label=label)
            for label, dist in distributions.items()
        }

        # Key parameters to compare
        params = [
            ("sig_x_mm", "σ_x [mm]"),
            ("sig_y_mm", "σ_y [mm]"),
            ("sig_z_mm", "σ_z [mm]"),
            ("sig_E_over_E_pct", "σ_E/E [%]"),
            ("emit_x_norm_um", "ε_nx [μm·rad]"),
            ("emit_y_norm_um", "ε_ny [μm·rad]"),
        ]

        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()

        for ax, (key, ylabel) in zip(axes, params):
            values = [s.to_dict()[key] for s in stats_list.values()]
            colors = plt.cm.Set2(np.linspace(0, 1, n))

            bars = ax.bar(range(n), values, color=colors, edgecolor="white")
            ax.set_xticks(range(n))
            ax.set_xticklabels(labels, fontsize=9)
            ax.set_ylabel(ylabel)
            ax.set_title(ylabel)

            # Add value labels on bars
            for bar, val in zip(bars, values):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height(),
                    f"{val:.3f}",
                    ha="center", va="bottom", fontsize=8,
                )

        if n > 6:
            # Remove extra axes
            for ax in axes[len(params):]:
                ax.set_visible(False)

        fig.suptitle("Beam Parameter Comparison", fontsize=14, fontweight="bold")
        fig.tight_layout()
        return fig

    else:
        raise ValueError(f"Unknown plot_type: {plot_type}")
