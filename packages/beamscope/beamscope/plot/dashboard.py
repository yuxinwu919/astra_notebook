"""Beam statistics dashboard — comprehensive multi-panel view."""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as sp_stats

from ..analysis.statistics import compute_statistics
from ..distribution import Distribution
from ._precompute import precompute, clip_percentile
from ._artists import add_colorbar, draw_reference_lines


def plot_dashboard(
    distributions: dict[str, Distribution],
    figsize: tuple[float, float] = (16, 12),
    title: Optional[str] = None,
) -> plt.Figure:
    """Generate a comprehensive multi-panel beam dashboard.

    Includes phase space (x-x', y-y', z-dp), projections, statistics table,
    and transverse beam profile.
    """
    labels = list(distributions.keys())
    main_label = labels[0]
    main_dist = distributions[main_label]
    data = precompute(main_dist)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.35)

    # ── Row 1: Phase space triple (x-x', y-y', z-dp) ──
    for i, (x_key, y_key, xl, yl, t) in enumerate([
        ("x", "xp", "x [mm]", "x' [mrad]", "x–x'"),
        ("y", "yp", "y [mm]", "y' [mrad]", "y–y'"),
        ("z", "dp", "z [mm]", "δp/p [%]", "Longitudinal"),
    ]):
        ax = fig.add_subplot(gs[0, i])
        xd = data[x_key]; yd = data[y_key]
        vmin_x, vmax_x = clip_percentile(xd)
        vmin_y, vmax_y = clip_percentile(yd)
        h = ax.hist2d(xd, yd, bins=60, cmap="viridis",
                      range=[[vmin_x, vmax_x], [vmin_y, vmax_y]])
        add_colorbar(fig, ax, h[3])
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(t)
        draw_reference_lines(ax)

    # Transverse profile (x-y)
    ax_xy = fig.add_subplot(gs[0, 3])
    vmin_x2, vmax_x2 = clip_percentile(data["x"])
    vmin_y2, vmax_y2 = clip_percentile(data["y"])
    h = ax_xy.hist2d(data["x"], data["y"], bins=60, cmap="inferno",
                     range=[[vmin_x2, vmax_x2], [vmin_y2, vmax_y2]])
    add_colorbar(fig, ax_xy, h[3])
    ax_xy.set_xlabel("x [mm]"); ax_xy.set_ylabel("y [mm]")
    ax_xy.set_title("Transverse Profile"); ax_xy.set_aspect("equal")

    # ── Row 2: Projection histograms ──
    for i, (key, dim_label, color) in enumerate([
        ("x", "x [mm]", "steelblue"),
        ("y", "y [mm]", "seagreen"),
        ("z", "z [mm]", "darkorange"),
    ]):
        ax = fig.add_subplot(gs[1, i])
        d = data[key]
        ax.hist(d, bins=60, density=True, alpha=0.7, color=color, edgecolor="white")
        mu, sigma = float(np.mean(d)), float(np.std(d))
        x_fit = np.linspace(d.min(), d.max(), 200)
        ax.plot(x_fit, sp_stats.norm.pdf(x_fit, mu, sigma), "r-", lw=2, label=f"σ={sigma:.3f}")
        ax.set_xlabel(dim_label); ax.set_ylabel("Density")
        ax.set_title(f"{dim_label.split('[')[0].strip()} Distribution")
        ax.legend(fontsize=9)

    # Statistics table
    ax_table = fig.add_subplot(gs[1, 3])
    ax_table.axis("off")
    _draw_stats_table(ax_table, main_dist, distributions, labels)

    # ── Row 3: Current profile + Info ──
    ax_curr = fig.add_subplot(gs[2, :2])
    curr_bins, z_edges = np.histogram(data["z"], bins=100)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    ax_curr.fill_between(z_centers, curr_bins, alpha=0.6, color="steelblue")
    ax_curr.set_xlabel("z [mm]"); ax_curr.set_ylabel("Count (~Current)")
    ax_curr.set_title("Longitudinal Profile")

    ax_info = fig.add_subplot(gs[2, 2:])
    ax_info.axis("off")
    info_lines = [
        f"Source: {main_dist.source or 'N/A'}",
        f"Format: {main_dist.format or 'N/A'}",
        f"Particles: {main_dist.n_particle} total, {main_dist.n_active} active",
        f"Total charge: {main_dist.active_charge_nC:.4f} nC",
        f"Ref. momentum: {main_dist.ref_momentum_eVc * 1e-6:.2f} MeV/c",
        f"Ref. kinetic E: {main_dist.reference_kinetic_energy_eV * 1e-6:.2f} MeV",
    ]
    for i, line in enumerate(info_lines):
        ax_info.text(0.05, 0.9 - i * 0.15, line, transform=ax_info.transAxes,
                     fontsize=10, fontfamily="monospace")

    fig.suptitle(title or "Beam Dashboard", fontsize=16, fontweight="bold", y=1.01)
    try:
        fig.tight_layout()
    except Exception:
        pass
    return fig


def _draw_stats_table(ax, main_dist, distributions, labels):
    """Draw statistics table on the given axes."""
    if len(distributions) == 1:
        stats = compute_statistics(main_dist)
        rows = [
            ["Parameter", "Value"],
            ["N_active", str(stats.n_active)],
            ["sig_x [mm]", f"{stats.sig_x * 1e3:.4f}"],
            ["sig_y [mm]", f"{stats.sig_y * 1e3:.4f}"],
            ["sig_z [mm]", f"{stats.sig_z * 1e3:.4f}"],
            ["sig_E/E [%]", f"{stats.sig_E_over_E * 100:.4f}"],
            ["eps_nx [um]", f"{stats.emit_x_norm * 1e6:.2f}"],
            ["eps_ny [um]", f"{stats.emit_y_norm * 1e6:.2f}"],
            ["beta_x [m]", f"{stats.beta_x:.4f}"],
            ["alpha_x", f"{stats.alpha_x:.4f}"],
        ]
    else:
        stats_list = [compute_statistics(d, label=l) for l, d in distributions.items()]
        rows = [["Parameter"] + labels]
        for key, name in [
            ("sig_x", "sig_x [mm]"), ("sig_y", "sig_y [mm]"),
            ("sig_z", "sig_z [mm]"), ("sig_E_over_E", "sig_E/E [%]"),
            ("emit_x_norm", "eps_nx [um]"), ("emit_y_norm", "eps_ny [um]"),
            ("n_active", "N_active"),
        ]:
            row = [name]
            for s in stats_list:
                val = getattr(s, key)
                if key in ("sig_x", "sig_y", "sig_z"):
                    row.append(f"{val * 1e3:.4f}")
                elif key == "sig_E_over_E":
                    row.append(f"{val * 100:.4f}")
                elif key in ("emit_x_norm", "emit_y_norm"):
                    row.append(f"{val * 1e6:.2f}")
                else:
                    row.append(f"{val}")
            rows.append(row)

    table = ax.table(
        cellText=rows, cellLoc="center", loc="center",
        colWidths=[0.25] + [0.25] * max(1, len(distributions)),
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.4)
    for j in range(len(rows[0])):
        table[(0, j)].set_facecolor("#455A64")
        table[(0, j)].set_text_props(weight="bold", color="white")
