"""3×2 overview canvas showing six-dimensional phase space.

Displays: x-x', y-y', z-dE (col 1), x-y, z-x, z-y (col 2).
Click on any subplot to open a detailed view.
"""

from __future__ import annotations

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from PySide6.QtCore import QTimer, Qt
from PySide6.QtWidgets import QSizePolicy

from beamscope.distribution import Distribution
from beamscope.plot.overview import find_panel
from beamscope.plot._precompute import precompute, clip_percentile, get_overview_panels
from beamscope._plotting.cosmetics import SLAC_DESY_CMAP


class OverviewCanvas(FigureCanvasQTAgg):
    """3×2 grid of phase space projections."""

    def __init__(self, plot_panel=None, parent=None):
        self._fig = Figure(figsize=(10, 8), dpi=120)
        super().__init__(self._fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setMinimumSize(400, 300)

        self._plot_panel = plot_panel
        self._distributions: dict = {}
        self._opts: dict = {}
        self._axes: list = []
        self._data_cache: dict = {}  # label -> precomputed data dict

        # Resize debounce timer
        self._resize_timer = QTimer()
        self._resize_timer.setSingleShot(True)
        self._resize_timer.setInterval(150)
        self._resize_timer.timeout.connect(self._delayed_render)

        # Click connection for drill-down
        self.mpl_connect("button_press_event", self._on_click)

        self._create_empty_overview()

    # ── Public ───────────────────────────────────────────────────

    def render(self, distributions: dict, opts: dict = None):
        """Render the overview with loaded distributions.

        Args:
            distributions: dict of {label: {'dist': Distribution, 'color': str}}
            opts: optional style dict (bins, cmap)
        """
        self._distributions = distributions
        self._opts = opts or {}
        self._data_cache.clear()
        # Precompute data for all distributions
        for label, info in distributions.items():
            self._data_cache[label] = precompute(info["dist"])
        self._delayed_render()

    def clear(self):
        self._distributions.clear()
        self._data_cache.clear()
        self._create_empty_overview()

    # ── Internals ────────────────────────────────────────────────

    def _create_empty_overview(self):
        """Create empty placeholder subplots."""
        self._fig.clear()
        self._axes = []
        for row in range(3):
            row_axes = []
            for col in range(2):
                ax = self._fig.add_subplot(3, 2, row * 2 + col + 1)
                ax.text(0.5, 0.5, "Load a file to begin\n\nFile → Open or drag .ini/.001 file",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=10, color="gray")
                ax.set_xticks([])
                ax.set_yticks([])
                row_axes.append(ax)
            self._axes.append(row_axes)
        self._fig.tight_layout(pad=2.0)
        self.draw_idle()

    def _delayed_render(self):
        """Render after resize settles."""
        if not self._distributions:
            self._create_empty_overview()
            return
        self._render_overview()

    def _active_label(self):
        """Return the label of the distribution to use as primary for overview."""
        labels = list(self._distributions.keys())
        if not labels:
            return None
        # Use the distribution whose color key has 'active' marker, or first
        for label, info in self._distributions.items():
            if info.get("active", False):
                return label
        return labels[0]

    def _render_overview(self):
        """Render the 3×2 overview Figure."""
        self._fig.clear()
        self._axes = []

        label = self._active_label()
        if label is None:
            self._create_empty_overview()
            return

        data = self._data_cache.get(label, {})
        cmap = self._opts.get("cmap", SLAC_DESY_CMAP)
        bins = self._opts.get("bins", 60)

        # Track the last mappable for a shared colorbar
        last_mappable = None

        for row, col, x_key, y_key, title, x_label, y_label in get_overview_panels():
            ax = self._fig.add_subplot(3, 2, row * 2 + col + 1)

            x_data = data.get(x_key)
            y_data = data.get(y_key)

            if x_data is not None and y_data is not None and len(x_data) > 0:
                vmin_x, vmax_x = clip_percentile(x_data)
                vmin_y, vmax_y = clip_percentile(y_data)

                # Avoid degenerate range (all values identical after clipping)
                if vmax_x <= vmin_x:
                    vmin_x, vmax_x = float(np.min(x_data)), float(np.max(x_data))
                if vmax_y <= vmin_y:
                    vmin_y, vmax_y = float(np.min(y_data)), float(np.max(y_data))
                if vmax_x <= vmin_x:
                    vmin_x, vmax_x = vmin_x - 0.1, vmax_x + 0.1
                if vmax_y <= vmin_y:
                    vmin_y, vmax_y = vmin_y - 0.1, vmax_y + 0.1

                h = ax.hist2d(
                    x_data, y_data, bins=bins, cmap=cmap, norm=LogNorm(),
                    range=[[vmin_x, vmax_x], [vmin_y, vmax_y]],
                )
                last_mappable = h[3]
                ax.axhline(0, color="white", ls="--", lw=0.8)
                ax.axvline(0, color="white", ls="--", lw=0.8)

            ax.set_xlabel(x_label, fontsize=10)
            ax.set_ylabel(y_label, fontsize=10)
            ax.set_title(title, fontsize=12, fontweight="bold")
            ax.tick_params(labelsize=9)
            # Suppress large offset text (e.g., "1e-30") on axes
            ax.ticklabel_format(useOffset=False, style='plain')

            if len(self._axes) <= row:
                self._axes.append([])
            self._axes[row].append(ax)

        # Single shared colorbar for all subplots
        if last_mappable is not None:
            self._fig.colorbar(last_mappable, ax=self._axes, fraction=0.03, pad=0.02,
                              label="Count")

        try:
            self._fig.tight_layout(pad=2.0, h_pad=1.5, w_pad=1.5)
        except Exception:
            pass
        self.draw_idle()

    def _on_click(self, event):
        """Handle click on subplot → open Detail tab."""
        if event.inaxes is None or not self._distributions:
            return

        # Find which subplot was clicked by checking axes positions
        for row in range(3):
            for col in range(2):
                if row < len(self._axes) and col < len(self._axes[row]):
                    if self._axes[row][col] == event.inaxes:
                        x_key, y_key, title, x_label, y_label = find_panel(row, col)

                        if event.button == 3:  # Right click → new window
                            self._open_detached(row, col)
                        else:
                            self._open_detail(row, col)
                        return

    def _open_detail(self, row: int, col: int):
        """Open Detail tab for the clicked subplot."""
        if self._plot_panel is None:
            return
        x_key, y_key, title, x_label, y_label = find_panel(row, col)
        self._plot_panel.open_detail(
            x_key, y_key, x_label, y_label,
            self._distributions, title=title,
        )

    def _open_detached(self, row: int, col: int):
        """Open the subplot in a detached window."""
        if self._plot_panel is None:
            return
        self._open_detail(row, col)
        self._plot_panel.detach_current_tab()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._distributions:
            self._resize_timer.start()


