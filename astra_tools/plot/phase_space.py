"""Phase-space plots: x-x', y-y' (postpro style) and z-dp/p.

Physics conventions (audited against the ASTRA Manual V3.2):
  * x' = p~x / p_ref with the canonical momentum p~x = px + c Bz y / 2
    (manual 4.13.1). Pass bz_on_axis_T for bunches inside solenoid
    fields; otherwise the trace-space divergence is used.
  * dp/p = (pz - <pz>)/<pz>, pz in absolute eV/c (the reader converts
    the file's relative pz to absolute).
  * positive z = ahead of the reference particle (bunch head).
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from ..analysis.emittance import canonical_divergence, compute_emittance_ellipse_params
from ..distribution import Distribution

PLANE = ("x", "y", "z")


def _clip_percentile(a, q=0.5):
    lo, hi = np.percentile(a, [q, 100.0 - q])
    if hi - lo < 1e-30:
        pad = max(abs(hi) * 1e-6, 1e-30)
        lo, hi = lo - pad, hi + pad
    return lo, hi


def plot_phase_space(
    dist: Distribution,
    plane: str = "x",
    kind: str = "density",
    bins: int = 80,
    ax=None,
    figsize=(6, 5),
    title: Optional[str] = None,
    cmap=None,
    colorbar: bool = True,
    show_ellipse: bool = False,
    bz_on_axis_T: float = 0.0,
    use_weights: bool = False,
) -> plt.Figure:
    """Plot one 2D phase-space projection.

    Args:
        dist: Distribution.
        plane: 'x' (x-x'), 'y' (y-y') or 'z' (z-dp/p).
        kind: 'density' (2D histogram, log scale) or 'scatter'.
        bins: histogram bins.
        ax: optional existing Axes.
        show_ellipse: overlay the RMS emittance ellipse ('x'/'y' only).
        bz_on_axis_T: solenoid on-axis field at the bunch center [T]
            for the canonical divergence (manual 4.13.1).
        use_weights: weight the histogram by macro-particle charge.
    """
    if plane not in PLANE:
        raise ValueError("plane must be one of " + str(PLANE))

    mask = dist.active
    pz_abs = np.abs(dist.pz[mask])
    w = dist.charge[mask] if use_weights else None

    if plane == "x":
        ptx = canonical_divergence(dist.px[mask], dist.y[mask], bz_on_axis_T, +1.0)
        x_data = dist.x[mask] * 1e3          # mm
        y_data = (ptx - np.mean(ptx)) / dist.ref_momentum_eVc * 1e3  # mrad
        xlabel, ylabel = "x [mm]", "x' [mrad]"
        default_title = "x-x' phase space"
    elif plane == "y":
        pty = canonical_divergence(dist.py[mask], dist.x[mask], bz_on_axis_T, -1.0)
        x_data = dist.y[mask] * 1e3
        y_data = (pty - np.mean(pty)) / dist.ref_momentum_eVc * 1e3
        xlabel, ylabel = "y [mm]", "y' [mrad]"
        default_title = "y-y' phase space"
    else:
        x_data = (dist.z[mask] - np.mean(dist.z[mask])) * 1e3
        mean_pz = np.mean(dist.pz[mask])
        y_data = (dist.pz[mask] - mean_pz) / mean_pz * 100  # %
        xlabel = "z [mm]  (positive = ahead of reference)"
        ylabel = "dp/p [%]"
        default_title = "longitudinal phase space"

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure

    vx0, vx1 = _clip_percentile(x_data)
    vy0, vy1 = _clip_percentile(y_data)

    if kind == "density":
        h = ax.hist2d(
            x_data, y_data, bins=bins,
            range=[[vx0, vx1], [vy0, vy1]],
            cmap=cmap or "viridis", norm=LogNorm(),
            weights=w,
        )
        if colorbar:
            fig.colorbar(h[3], ax=ax, label="counts")
    elif kind == "scatter":
        ax.scatter(x_data, y_data, s=1, alpha=0.5)
    else:
        raise ValueError("kind must be 'density' or 'scatter'")

    ax.axhline(0, color="w", lw=0.5, ls="--")
    ax.axvline(0, color="w", lw=0.5, ls="--")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title or default_title)

    if show_ellipse and plane in ("x", "y"):
        if plane == "x":
            u = dist.x[mask] - np.mean(dist.x[mask])
            up = (ptx - np.mean(ptx)) / dist.ref_momentum_eVc
        else:
            u = dist.y[mask] - np.mean(dist.y[mask])
            up = (pty - np.mean(pty)) / dist.ref_momentum_eVc
        params = compute_emittance_ellipse_params(u, up, n_sigma=1.0, weights=w)
        theta = np.linspace(0, 2 * np.pi, 200)
        xe = params["a"] * np.cos(theta) * 1e3
        ye = params["b"] * np.sin(theta) * 1e3
        xr = xe * np.cos(params["theta"]) - ye * np.sin(params["theta"])
        yr = xe * np.sin(params["theta"]) + ye * np.cos(params["theta"])
        ax.plot(xr, yr, "r-", lw=1.5, label="1-RMS ellipse")
        ax.legend(loc="best", fontsize=8)

    fig.tight_layout()
    return fig


def plot_transverse_phase_space(
    dist: Distribution,
    figsize=(12, 5),
    title_prefix: str = "",
    **kwargs,
) -> plt.Figure:
    """Both transverse phase spaces side by side (x-x', y-y')."""
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    plot_phase_space(dist, plane="x", ax=axes[0],
                     title=(title_prefix + " x-x'").strip(), **kwargs)
    plot_phase_space(dist, plane="y", ax=axes[1],
                     title=(title_prefix + " y-y'").strip(), **kwargs)
    return fig
