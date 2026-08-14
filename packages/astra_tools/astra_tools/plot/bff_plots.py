"""Bunch form factor plots (CSR standard display).

BFF(k) and |F(k)| on a log-log scale; the wavelength lambda = 2 pi / k
is shown on the top axis (standard for CSR/ISR discussions).
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.bff import BFFResult


def plot_bff(
    bff: BFFResult,
    figsize=(9, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """BFF(k) with wavelength top axis."""
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.loglog(bff.k, bff.bff, lw=1.5, color="C0")
    ax.set_xlabel("k [1/m]")
    ax.set_ylabel("|F(k)|^2")
    ax.set_title(title or "bunch form factor")

    ax2 = ax.twiny()
    k_decades = np.logspace(np.log10(bff.k[0]), np.log10(bff.k[-1]), 6)
    ax2.set_xscale("log")
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(k_decades)
    ax2.set_xticklabels(["%.3g" % (2 * np.pi / k) for k in k_decades])
    ax2.set_xlabel("wavelength lambda [m]")
    fig.tight_layout()
    return fig


def plot_bff_with_amplitude(
    bff: BFFResult,
    figsize=(9, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """|F(k)| and BFF(k) together."""
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.loglog(bff.k, bff.bff, lw=1.5, color="C0", label="|F(k)|^2")
    ax.loglog(bff.k, bff.bff_amplitude, lw=1.2, color="C1", label="|F(k)|")
    if bff.csr_cutoff_k:
        ax.axvline(bff.csr_cutoff_k, color="r", ls="--", lw=1,
                   label="cutoff (BFF=0.5)")
    ax.set_xlabel("k [1/m]")
    ax.set_ylabel("form factor")
    ax.set_title(title or "bunch form factor")
    ax.legend()
    fig.tight_layout()
    return fig
