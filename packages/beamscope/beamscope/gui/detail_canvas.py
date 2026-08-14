"""Detail, Slice, Emit, and Log canvas views."""

from __future__ import annotations

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QSizePolicy, QWidget

from beamscope._plotting.cosmetics import SLAC_DESY_CMAP


class DetailCanvas(FigureCanvasQTAgg):
    """Single-dimension detailed phase space with marginals and overlays."""

    def __init__(self, plot_panel=None, parent=None):
        self._fig = Figure(figsize=(8, 7), dpi=120)
        super().__init__(self._fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self._plot_panel = plot_panel
        self._distributions: dict = {}
        self._x_key = "x"
        self._y_key = "xp"
        self._x_label = "x [mm]"
        self._y_label = "x' [mrad]"
        self._title = ""
        self._opts: dict = {}

        self._resize_timer = QTimer()
        self._resize_timer.setSingleShot(True)
        self._resize_timer.setInterval(150)
        self._resize_timer.timeout.connect(self._delayed_render)

    def update_opts(self, opts: dict):
        self._opts = opts

    def render(self, distributions: dict, x_key: str, y_key: str,
               x_label: str, y_label: str, title: str = "", opts: dict = None):
        self._distributions = distributions
        self._x_key = x_key
        self._y_key = y_key
        self._x_label = x_label
        self._y_label = y_label
        self._title = title
        if opts:
            self._opts = opts
        self._delayed_render()

    def _delayed_render(self):
        if not self._distributions:
            return
        from beamscope.plot.detail import plot_detail

        opts = self._opts
        self._fig.clear()
        labels = list(self._distributions.keys())
        main_label = labels[0]
        main_dist = self._distributions[main_label]["dist"]

        plot_detail(
            self._fig, main_dist,
            x_key=self._x_key, y_key=self._y_key,
            x_label=self._x_label, y_label=self._y_label,
            title=self._title,
            bins=opts.get("bins", 60),
            cmap=opts.get("cmap", SLAC_DESY_CMAP),
            show_density=opts.get("density", True),
            show_scatter=opts.get("scatter", False),
            show_contour=opts.get("contour", False),
            show_ellipse=opts.get("ellipse", True),
            show_marginals=opts.get("marginals", True),
        )
        try:
            self._fig.tight_layout()
        except Exception:
            pass
        self.draw_idle()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._distributions:
            self._resize_timer.start()


class SliceCanvas(FigureCanvasQTAgg):
    """Slice analysis visualization canvas."""

    def __init__(self, plot_panel=None, parent=None):
        self._fig = Figure(figsize=(10, 8), dpi=120)
        super().__init__(self._fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def render(self, slice_result, label: str = ""):
        """Render slice analysis in a 3x2 dashboard."""
        import numpy as np
        sa = slice_result
        self._fig.clear()

        # Convert to sliceMatrix format for PlotSliceParameters
        n = len(sa.z_centers)
        sm = np.zeros((n, 14))
        sm[:, 0] = sa.z_centers
        sm[:, 2] = sa.n_particles
        sm[:, 3] = sa.charge
        sm[:, 4] = sa.current
        sm[:, 6] = sa.mean_x
        sm[:, 7] = sa.mean_y
        sm[:, 8] = sa.emit_x_norm
        sm[:, 9] = sa.emit_y_norm
        sm[:, 10] = sa.sig_E_over_E

        from beamscope._plotting.plotter import PlotSliceParameters
        # PlotSliceParameters creates its own figure; copy axes to ours
        tmp_fig = PlotSliceParameters(sm, figsize=(10, 8))
        for ax in tmp_fig.axes:
            ax.figure = self._fig
            ax.set_figure(self._fig)
            self._fig.axes.append(ax)
        import matplotlib.pyplot as plt
        plt.close(tmp_fig)

        self._fig.suptitle(f"Slice Analysis — {label}" if label else "Slice Analysis",
                           fontweight="bold")
        try:
            self._fig.tight_layout()
        except Exception:
            pass
        self.draw_idle()


class EmitCanvas(FigureCanvasQTAgg):
    """Canvas for emit/sigma/ref evolution plots."""

    def __init__(self, plot_panel=None, parent=None):
        self._fig = Figure(figsize=(12, 8), dpi=120)
        super().__init__(self._fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def render_emit_dashboard(self, emit_set, sigma=None):
        from beamscope.plot.emit_plots import plot_emit_dashboard
        plot_emit_dashboard(emit_set, sigma, fig=self._fig)
        try:
            self._fig.tight_layout()
        except Exception:
            pass
        self.draw_idle()

    def render_emit_plot(self, emit_set, kind: str):
        from beamscope.plot.emit_plots import (
            plot_envelope_evolution, plot_emittance_evolution, plot_energy_evolution,
            plot_transverse_size,
        )
        fn = {"envelope": plot_envelope_evolution,
              "emittance": plot_emittance_evolution,
              "energy": plot_energy_evolution,
              "transverse": plot_transverse_size}[kind]
        fn(emit_set, fig=self._fig)
        try:
            self._fig.tight_layout()
        except Exception:
            pass
        self.draw_idle()

    def render_eigen(self, sigma):
        from beamscope.plot.emit_plots import plot_eigen_emittances
        plot_eigen_emittances(sigma, fig=self._fig)
        try:
            self._fig.tight_layout()
        except Exception:
            pass
        self.draw_idle()

    def render_ref(self, ref):
        from beamscope.plot.emit_plots import plot_ref_trajectory
        plot_ref_trajectory(ref, fig=self._fig)
        try:
            self._fig.tight_layout()
        except Exception:
            pass
        self.draw_idle()


class LogCanvas(QWidget):
    """Simple text display for ASTRA log files."""

    def __init__(self, plot_panel=None, parent=None):
        super().__init__(parent)
        from PySide6.QtWidgets import QVBoxLayout, QTextEdit
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self._text = QTextEdit()
        self._text.setReadOnly(True)
        self._text.setFontFamily("Menlo, Monaco, monospace")
        self._text.setFontPointSize(10)
        layout.addWidget(self._text)

    def render(self, text: str):
        self._text.setPlainText(text)


class BFFCanvas(FigureCanvasQTAgg):
    """Bunch Form Factor visualization canvas."""

    def __init__(self, plot_panel=None, parent=None):
        self._fig = Figure(figsize=(10, 5), dpi=120)
        super().__init__(self._fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def render(self, bff_result, label: str = ""):
        from beamscope.plot.bff_plots import plot_bff
        self._fig.clear()
        plot_bff(bff_result, figsize=(10, 5))
        # Transfer axes from the new figure
        import matplotlib.pyplot as plt
        src_fig = plt.gcf()
        if src_fig is not self._fig:
            for ax in src_fig.axes:
                ax.figure = self._fig
                ax.set_figure(self._fig)
                self._fig.axes.append(ax)
            plt.close(src_fig)
        self._fig.suptitle(f"Bunch Form Factor — {label}" if label else "Bunch Form Factor",
                           fontweight="bold")
        try:
            self._fig.tight_layout()
        except Exception:
            pass
        self.draw_idle()
