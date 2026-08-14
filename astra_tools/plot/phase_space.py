"""Phase-space density plots (postpro style, modern KDE rendering).

Physics conventions (audited against ASTRA Manual V3.2 and validated
against ASTRA's own Xemit output):
  * x' = p~x / p_ref with canonical momentum p~x = px + c Bz y / 2
    (manual 4.13.1); pass bz_on_axis_T for bunches inside solenoids
  * dp/p from absolute pz (the reader converts the file's relative pz)
  * positive z = ahead of the reference particle (bunch head)

Density rendering: unified 2D Gaussian-kernel KDE (plot._density);
display range auto-clipped to the 0.5-99.5 percentile so outliers can
never collapse the bulk of the distribution.
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from ..analysis.emittance import canonical_divergence, compute_emittance_ellipse_params
from ..distribution import Distribution
from ._density import clip_percentile, density2d, outside_fraction

PLANES = ("x", "y", "z")


def plot_phase_space(
    dist: Distribution,
    plane: str = "x",
    bins: int = 160,
    ax=None,
    figsize=(6.5, 5),
    title: Optional[str] = None,
    cmap: Optional[str] = None,
    colorbar: bool = True,
    show_ellipse: bool = False,
    bz_on_axis_T: float = 0.0,
    use_weights: bool = False,
    clip_q: float = 0.5,
) -> plt.Figure:
    """KDE density plot of one 2D phase-space projection.

    Args:
        dist: Distribution.
        plane: 'x' (x-x'), 'y' (y-y') or 'z' (z-dp/p).
        bins: KDE grid resolution.
        show_ellipse: overlay the RMS emittance ellipse ('x'/'y').
        bz_on_axis_T: on-axis solenoid field at bunch center [T].
        use_weights: weight by macro-particle charge.
        clip_q: display-range percentile clip [%] on each side.
    """
    if plane not in PLANES:
        raise ValueError("plane must be one of " + str(PLANES))

    mask = dist.active
    w = dist.charge[mask] if use_weights else None

    # Reference momentum guard: fall back to mean |pz| for synthetic
    # distributions without a header reference.
    p_ref = dist.ref_momentum_eVc
    if p_ref <= 0:
        p_ref = float(np.mean(np.abs(dist.pz[mask])))
    if p_ref <= 0:
        raise ValueError("reference momentum is zero; cannot form x'")

    if plane == "x":
        ptx = canonical_divergence(dist.px[mask], dist.y[mask], bz_on_axis_T, +1.0)
        x_data = dist.x[mask] * 1e3
        y_data = (ptx - np.mean(ptx)) / p_ref * 1e3
        xlabel, ylabel = "x [mm]", "x' [mrad]"
        default_title = "x-x' phase space"
        u, up = dist.x[mask] - np.mean(dist.x[mask]), (ptx - np.mean(ptx)) / p_ref
    elif plane == "y":
        pty = canonical_divergence(dist.py[mask], dist.x[mask], bz_on_axis_T, -1.0)
        x_data = dist.y[mask] * 1e3
        y_data = (pty - np.mean(pty)) / p_ref * 1e3
        xlabel, ylabel = "y [mm]", "y' [mrad]"
        default_title = "y-y' phase space"
        u, up = dist.y[mask] - np.mean(dist.y[mask]), (pty - np.mean(pty)) / p_ref
    else:
        x_data = (dist.z[mask] - np.mean(dist.z[mask])) * 1e3
        mean_pz = np.mean(dist.pz[mask])
        y_data = (dist.pz[mask] - mean_pz) / mean_pz * 100
        xlabel = "z [mm]  (positive = ahead of reference)"
        ylabel = "dp/p [%]"
        default_title = "longitudinal phase space"
        u = up = None

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure

    xr = clip_percentile(x_data, clip_q)
    yr = clip_percentile(y_data, clip_q)
    z, xe, ye = density2d(x_data, y_data, bins=bins, range_xy=(xr, yr), weights=w)

    z_plot = np.where(z > 0, z, np.nan)
    im = ax.pcolormesh(xe, ye, z_plot.T, cmap=cmap, norm=LogNorm(),
                       shading="auto", rasterized=True)
    if colorbar:
        cb = fig.colorbar(im, ax=ax)
        cb.set_label("probability density [1/mm/mrad]" if plane in ("x", "y")
                     else "probability density")

    ax.axhline(0, color="white", lw=0.6, ls="--", alpha=0.7)
    ax.axvline(0, color="white", lw=0.6, ls="--", alpha=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title or default_title)

    if show_ellipse and plane in ("x", "y") and u is not None:
        params = compute_emittance_ellipse_params(u, up, n_sigma=1.0, weights=w)
        th = np.linspace(0, 2 * np.pi, 240)
        xe_ = params["a"] * np.cos(th) * 1e3
        ye_ = params["b"] * np.sin(th) * 1e3
        xr_ = xe_ * np.cos(params["theta"]) - ye_ * np.sin(params["theta"])
        yr_ = xe_ * np.sin(params["theta"]) + ye_ * np.cos(params["theta"])
        ax.plot(xr_, yr_, color="#CC3311", lw=1.8, ls="-", label="1-RMS ellipse")
        ax.legend(loc="best")

    # 离群点裁剪说明
    f_out = outside_fraction(x_data, *xr) + outside_fraction(y_data, *yr)
    if f_out > 0.01:
        ax.text(0.99, 0.02, "range clip: %.1f%% of points outside" % (100 * f_out),
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color="0.35")

    fig.tight_layout()
    return fig


def plot_transverse_phase_space(
    dist: Distribution,
    figsize=(13, 5),
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
