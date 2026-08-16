"""Field-map plots (fieldplot replacement).

Cavity: on-axis Ez(z) and 2D (r, z) maps of Ez, Er, Bphi reconstructed
from the axis field via the off-axis expansion (manual chapter 8).
Solenoid: on-axis Bz(z) and 2D (r, z) map of Bz with Br arrows.
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from ..io.field_map import CavityField, SolenoidField


def plot_cavity_field(
    field: CavityField,
    omega: float = 0.0,
    rmax: Optional[float] = None,
    figsize=(14, 8),
    title: Optional[str] = None,
    maxE_MVpm: Optional[float] = None,
) -> plt.Figure:
    """Cavity field: Ez(r,z), Er(r,z), Bphi(r,z) maps + on-axis Ez(z).

    Args:
        field: CavityField (on-axis table).
        omega: RF angular frequency [rad/s] for Bphi; 0 = static.
        rmax: max radius for the map [m]; default 10% of cavity length.
        maxE_MVpm: scale the field to this peak [MV/m] (like ASTRA's
            MaxE(n)); None = raw file values.
    """
    # 批 3: ez0 为原始任意单位 (手册 6.9); 给 maxE 时按峰值缩放并
    # 以 MV/m 显示, 否则按任意单位显示 (不再假装 V/m)。
    peak = float(np.max(np.abs(field.ez0)))
    if maxE_MVpm is not None:
        scale = maxE_MVpm / peak if peak > 0 else 0.0
        ez_unit, b_unit = "MV/m", "T"
    else:
        scale = 1.0
        ez_unit, b_unit = "arb. units", "arb. units"

    z0, z1 = float(field.z.min()), float(field.z.max())
    if rmax is None:
        rmax = 0.1 * (z1 - z0)
    nz, nr = 400, 60
    zz = np.linspace(z0, z1, nz)
    rr = np.linspace(0, rmax, nr)
    ZZ, RR = np.meshgrid(zz, rr)
    ez, er, bphi = field.field_at(RR, ZZ, omega)
    ez = ez * scale
    er = er * scale
    bphi = bphi * scale

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    maps = [(ez, "Ez [%s]" % ez_unit), (er, "Er [%s]" % ez_unit),
            (bphi, "Bphi [%s]" % b_unit)]
    for ax, (m, label) in zip(axes.flat[:3], maps):
        vmax = float(np.max(np.abs(m)))
        if vmax == 0:
            vmax = 1.0
        im = ax.pcolormesh(ZZ * 1e3, RR * 1e3, m, cmap="RdBu_r",
                           vmin=-vmax, vmax=vmax, shading="auto")
        fig.colorbar(im, ax=ax, label=label)
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("r [mm]")

    ax = axes[1, 1]
    ax.plot(field.z * 1e3, field.ez0 * scale, color="C0")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Ez on axis [%s]" % ez_unit)
    ax.set_title("on-axis field")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_solenoid_field(
    field: SolenoidField,
    rmax: Optional[float] = None,
    figsize=(12, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """Solenoid: Bz(r,z) map and on-axis Bz(z).

    Args:
        field: SolenoidField, already scaled to Tesla (use .scaled(maxB)).
    """
    z0, z1 = float(field.z.min()), float(field.z.max())
    if rmax is None:
        rmax = 0.05 * (z1 - z0)
    nz, nr = 400, 40
    zz = np.linspace(z0, z1, nz)
    rr = np.linspace(0, rmax, nr)
    ZZ, RR = np.meshgrid(zz, rr)
    br, bz = field.field_at(RR, ZZ)

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    ax = axes[0]
    vmax = float(np.max(np.abs(bz))) or 1.0
    im = ax.pcolormesh(ZZ * 1e3, RR * 1e3, bz, cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, shading="auto")
    fig.colorbar(im, ax=ax, label="Bz [T]")
    # Br arrows on a coarse grid
    stride = 12
    ax.quiver(ZZ[::stride, ::stride] * 1e3, RR[::stride, ::stride] * 1e3,
              np.zeros_like(ZZ[::stride, ::stride]), br[::stride, ::stride],
              color="k", scale=abs(br).max() * 2 or 1.0, width=0.002)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title("Bz(r,z)")

    ax = axes[1]
    ax.plot(field.z * 1e3, field.bz0, color="C0")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Bz on axis [T]")
    ax.set_title("on-axis field")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig
