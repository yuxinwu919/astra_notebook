"""Custom axis mapping canvas: free selection of X and Y variables."""

from __future__ import annotations

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QSizePolicy

from beamscope.plot._precompute import precompute, get_variable_defs, clip_percentile
from beamscope._plotting.cosmetics import SLAC_DESY_CMAP


class CustomCanvas(FigureCanvasQTAgg):
    """Free-form axis mapping canvas.

    Uses the shared precompute() function from beamscope.plot.overview
    for all variable computations — no duplicated logic.
    """

    def __init__(self, plot_panel=None, parent=None):
        self._fig = Figure(figsize=(8, 6), dpi=120)
        super().__init__(self._fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self._plot_panel = plot_panel
        self._distributions: dict = {}
        self._data_cache: dict = {}
        self._x_var = "x"
        self._y_var = "xp"
        self._opts: dict = {}

        self._resize_timer = QTimer()
        self._resize_timer.setSingleShot(True)
        self._resize_timer.setInterval(150)
        self._resize_timer.timeout.connect(self._delayed_render)

    def render(self, distributions: dict, x_var: str = "x", y_var: str = "xp",
               opts: dict = None):
        self._distributions = distributions
        self._x_var = x_var
        self._y_var = y_var
        self._opts = opts or {}
        # Precompute data for all distributions
        self._data_cache.clear()
        for label, info in distributions.items():
            self._data_cache[label] = precompute(info["dist"])
        self._delayed_render()

    def _delayed_render(self):
        if not self._distributions:
            return
        var_defs = get_variable_defs()
        if self._x_var not in var_defs or self._y_var not in var_defs:
            return

        self._fig.clear()
        ax = self._fig.add_subplot(111)

        x_fmt, x_key = var_defs[self._x_var]
        y_fmt, y_key = var_defs[self._y_var]

        bins = self._opts.get("bins", 60)
        cmap = self._opts.get("cmap", SLAC_DESY_CMAP)
        last_mappable = None

        for label, info in self._distributions.items():
            data = self._data_cache.get(label, {})
            x_data = data.get(x_key)
            y_data = data.get(y_key)
            if x_data is None or y_data is None:
                continue

            if self._opts.get("density", True):
                vmin_x, vmax_x = clip_percentile(x_data)
                vmin_y, vmax_y = clip_percentile(y_data)
                if vmax_x <= vmin_x:
                    vmin_x, vmax_x = float(np.min(x_data)), float(np.max(x_data))
                if vmax_y <= vmin_y:
                    vmin_y, vmax_y = float(np.min(y_data)), float(np.max(y_data))
                h = ax.hist2d(
                    x_data, y_data, bins=bins, cmap=cmap, norm=LogNorm(),
                    range=[[vmin_x, vmax_x], [vmin_y, vmax_y]],
                )
                last_mappable = h[3]

            if self._opts.get("scatter", False):
                n_sample = min(len(x_data), 5000)
                idx = np.random.choice(len(x_data), n_sample, replace=False)
                color = info.get("color", None)
                ax.scatter(x_data[idx], y_data[idx], s=1, alpha=0.3, color=color)

            if self._opts.get("contour", False):
                try:
                    counts, xedges, yedges = np.histogram2d(x_data, y_data, bins=bins)
                    xc = 0.5 * (xedges[:-1] + xedges[1:])
                    yc = 0.5 * (yedges[:-1] + yedges[1:])
                    X, Y = np.meshgrid(xc, yc)
                    ax.contour(X, Y, counts.T, levels=4, colors="white",
                               linewidths=0.5, alpha=0.6)
                except Exception:
                    pass

        ax.set_xlabel(x_fmt)
        ax.set_ylabel(y_fmt)
        ax.axhline(0, color="gray", ls="--", lw=0.8)
        ax.axvline(0, color="gray", ls="--", lw=0.8)

        # Single shared colorbar
        if last_mappable is not None:
            self._fig.colorbar(last_mappable, ax=ax, label="Count")

        self._fig.tight_layout()
        self.draw_idle()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._distributions:
            self._resize_timer.start()


