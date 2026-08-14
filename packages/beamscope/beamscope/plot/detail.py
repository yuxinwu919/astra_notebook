"""Detailed single-dimension phase space plot with rich information.

Features: 2x2 layout with 2D density + marginal projections + statistics panel.
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy import stats as sp_stats

from beamscope.distribution import Distribution
from beamscope.plot._precompute import precompute, clip_percentile, get_variable_label
from beamscope.plot._artists import (
    draw_emittance_ellipse, draw_reference_lines, add_colorbar,
)
from beamscope._plotting.cosmetics import SLAC_DESY_CMAP


def _get_stats_for_panel(stats, x_key: str, y_key: str) -> dict:
    """Extract relevant statistics for the current phase-space panel."""
    s = {}

    # Particle count
    n = getattr(stats, 'n_particle', None)
    s["N"] = f"{n:,}" if n else "—"

    # Per-panel parameters
    if x_key == "x" and y_key == "xp":
        s["ε_nx"] = f"{stats.emit_x_norm * 1e6:.3f} μm·rad"
        s["σ_x"] = f"{stats.sig_x * 1e3:.4f} mm"
        s["σ_x'"] = f"{stats.sig_px * 1e3:.4f} mrad"
        if hasattr(stats, 'beta_x'):
            s["β_x"] = f"{stats.beta_x:.2f} m"
            s["α_x"] = f"{stats.alpha_x:.3f}"
        s["⟨x⟩"] = f"{stats.mean_x * 1e3:.4f} mm"
    elif x_key == "y" and y_key == "yp":
        s["ε_ny"] = f"{stats.emit_y_norm * 1e6:.3f} μm·rad"
        s["σ_y"] = f"{stats.sig_y * 1e3:.4f} mm"
        s["σ_y'"] = f"{stats.sig_py * 1e3:.4f} mrad"
        if hasattr(stats, 'beta_y'):
            s["β_y"] = f"{stats.beta_y:.2f} m"
            s["α_y"] = f"{stats.alpha_y:.3f}"
        s["⟨y⟩"] = f"{stats.mean_y * 1e3:.4f} mm"
    elif x_key == "z" and y_key in ("dp", "E"):
        s["σ_z"] = f"{stats.sig_z * 1e3:.4f} mm"
        if hasattr(stats, 'sig_E_over_E'):
            s["σ_E/E"] = f"{stats.sig_E_over_E * 100:.3f} %"
        if hasattr(stats, 'ref_kinetic_energy_eV'):
            s["⟨E_kin⟩"] = f"{stats.ref_kinetic_energy_eV * 1e-6:.2f} MeV"
        s["⟨z⟩"] = f"{stats.mean_z * 1e3:.4f} mm"
    else:
        # Generic fallback
        s["σ_x"] = f"{stats.sig_x * 1e3:.4f} mm"
        s["σ_y"] = f"{stats.sig_y * 1e3:.4f} mm"
        s["σ_z"] = f"{stats.sig_z * 1e3:.4f} mm"

    if hasattr(stats, 'label') and stats.label:
        s["File"] = stats.label

    return s


def plot_detail(
    fig: plt.Figure,
    dist: Distribution,
    x_key: str = "x",
    y_key: str = "xp",
    title: str = "",
    bins: int = 60,
    cmap=None,
    show_density: bool = True,
    show_scatter: bool = False,
    show_contour: bool = False,
    show_ellipse: bool = True,
    show_marginals: bool = True,
    use_weights: bool = False,
    x_label: Optional[str] = None,
    y_label: Optional[str] = None,
) -> plt.Figure:
    """Plot a detailed single-dimension phase space view.

    2x2 layout:
        ┌───────────────────┬──────────────┐
        │                   │  Y marginal  │
        │  2D Density       │  (histogram  │
        │  (SLAC-DESY       │   + Gaussian)│
        │   + LogNorm)      ├──────────────┤
        │                   │  Statistics  │
        ├───────────────────│  panel       │
        │  X marginal       │              │
        │  (histogram       │              │
        │   + Gaussian)     │              │
        └───────────────────┴──────────────┘

    Args:
        fig: matplotlib Figure to draw on.
        dist: Particle distribution.
        x_key, y_key: Variable keys (see ``_precompute.VARIABLE_DEFS``).
        title: Plot title.
        bins: Number of bins.
        cmap: Colormap (default: SLAC-DESY beam map).
        show_density: Overlay 2D density histogram.
        show_scatter: Overlay scatter plot.
        show_contour: Overlay contour lines.
        show_ellipse: Overlay RMS emittance ellipse.
        show_marginals: Show 1D projection histograms with Gaussian fits.
        use_weights: Weight by macro-particle charge.
        x_label, y_label: Optional axis labels (auto-derived if None).

    Returns:
        The Figure (same as input).
    """
    if cmap is None:
        cmap = SLAC_DESY_CMAP

    fig.clear()

    data = precompute(dist)
    x_data = data.get(x_key, np.array([]))
    y_data = data.get(y_key, np.array([]))
    x_label = x_label or get_variable_label(x_key)
    y_label = y_label or get_variable_label(y_key)

    if use_weights:
        mask = dist.active
        weights = dist.charge[mask]
    else:
        weights = None

    if len(x_data) == 0:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No active particles", transform=ax.transAxes,
                ha="center", va="center", fontsize=12)
        return fig

    # ── 2x2 GridSpec layout ──
    # Left column: density (2 rows tall) + x-marginal
    # Right column: y-marginal + stats panel
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[3, 1.2],
        height_ratios=[3, 1],
        hspace=0.3, wspace=0.3,
    )

    ax_density = fig.add_subplot(gs[0, 0])
    ax_ymarg = fig.add_subplot(gs[0, 1])
    ax_xmarg = fig.add_subplot(gs[1, 0], sharex=ax_density)
    ax_stats = fig.add_subplot(gs[1, 1])

    # ── 2D Density (top-left) ──
    if show_density:
        vmin_x, vmax_x = clip_percentile(x_data)
        vmin_y, vmax_y = clip_percentile(y_data)
        h = ax_density.hist2d(
            x_data, y_data, bins=bins, cmap=cmap, norm=LogNorm(),
            range=[[vmin_x, vmax_x], [vmin_y, vmax_y]],
            weights=weights,
        )
        add_colorbar(fig, ax_density, h[3], fraction=0.08, pad=0.02)

    if show_scatter:
        n = min(len(x_data), 3000)
        rng = np.random.default_rng(42)
        idx = rng.choice(len(x_data), n, replace=False)
        ax_density.scatter(x_data[idx], y_data[idx], s=1, alpha=0.3,
                           color="white")

    if show_contour:
        try:
            counts, xedges, yedges = np.histogram2d(
                x_data, y_data, bins=bins, weights=weights,
            )
            xc = 0.5 * (xedges[:-1] + xedges[1:])
            yc = 0.5 * (yedges[:-1] + yedges[1:])
            X, Y = np.meshgrid(xc, yc)
            ax_density.contour(X, Y, counts.T, levels=4, colors="white",
                               linewidths=0.5, alpha=0.6)
        except Exception:
            pass

    if show_ellipse:
        x_centered = x_data - np.mean(x_data)
        y_centered = y_data - np.mean(y_data)
        draw_emittance_ellipse(
            ax_density, x_centered, y_centered,
            n_sigma=1.0, weights=weights,
        )

    draw_reference_lines(ax_density, color="white")
    ax_density.set_xlabel(x_label)
    ax_density.set_ylabel(y_label)
    if title:
        ax_density.set_title(title, fontweight="bold")
    else:
        ax_density.set_title(f"{x_key}–{y_key}", fontweight="bold")

    # ── X Marginal (bottom-left) ──
    if show_marginals:
        ax_xmarg.hist(
            x_data, bins=bins, density=True,
            color="steelblue", alpha=0.7, edgecolor="white",
        )
        mu_x, sig_x = float(np.mean(x_data)), float(np.std(x_data, ddof=1))
        xf = np.linspace(x_data.min(), x_data.max(), 200)
        ax_xmarg.plot(xf, sp_stats.norm.pdf(xf, mu_x, sig_x), "r-", lw=1.5)
        ax_xmarg.set_xlabel(x_label)
        ax_xmarg.set_ylabel("Density", fontsize=10)
        ax_xmarg.set_title(f"μ={mu_x:.3f}, σ={sig_x:.3f}", fontsize=10)
        ax_xmarg.tick_params(labelsize=9)
        ax_xmarg.grid(True, alpha=0.3)
    else:
        ax_xmarg.axis("off")

    # ── Y Marginal (top-right) ──
    if show_marginals:
        ax_ymarg.hist(
            y_data, bins=bins, density=True, orientation="horizontal",
            color="steelblue", alpha=0.7, edgecolor="white",
        )
        mu_y, sig_y = float(np.mean(y_data)), float(np.std(y_data, ddof=1))
        yf = np.linspace(y_data.min(), y_data.max(), 200)
        ax_ymarg.plot(sp_stats.norm.pdf(yf, mu_y, sig_y), yf, "r-", lw=1.5)
        ax_ymarg.set_ylabel(y_label)
        ax_ymarg.set_xlabel("Density", fontsize=10)
        ax_ymarg.set_title(f"μ={mu_y:.3f}\nσ={sig_y:.3f}",
                           fontsize=10, loc="center")
        ax_ymarg.tick_params(labelsize=9)
        ax_ymarg.grid(True, alpha=0.3)
    else:
        ax_ymarg.axis("off")

    # ── Statistics Panel (bottom-right) ──
    try:
        from beamscope.analysis.statistics import compute_statistics
        stats = compute_statistics(dist)
        stats_dict = _get_stats_for_panel(stats, x_key, y_key)

        lines = []
        for k, v in stats_dict.items():
            lines.append(f"{k}: {v}")
        textstr = "\n".join(lines)

        ax_stats.text(
            0.05, 0.95, textstr,
            transform=ax_stats.transAxes,
            fontsize=9, verticalalignment="top", family="monospace",
            bbox=dict(boxstyle="round", facecolor="lightyellow",
                      alpha=0.9, edgecolor="gray"),
        )
        ax_stats.set_title("Statistics", fontsize=11, fontweight="bold")
    except Exception:
        ax_stats.text(0.5, 0.5, "Stats unavailable", transform=ax_stats.transAxes,
                      ha="center", va="center", fontsize=10, color="gray")
    ax_stats.set_xticks([])
    ax_stats.set_yticks([])

    return fig
