"""Central plot panel: manages Overview, Detail, and Custom tabs."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QTabWidget, QToolButton,
)
from PySide6.QtCore import Qt, Signal

from beamscope.gui.overview_canvas import OverviewCanvas
from beamscope.gui.detail_canvas import DetailCanvas, EmitCanvas, LogCanvas
from beamscope.gui.custom_canvas import CustomCanvas


class PlotPanel(QWidget):
    """Tab-based plot panel managing Overview, Detail, and Custom views."""

    def __init__(self, main_window=None, parent=None):
        super().__init__(parent)
        self._main_window = main_window

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._tabs = QTabWidget()
        self._tabs.setTabsClosable(True)
        self._tabs.setMovable(True)
        self._tabs.tabCloseRequested.connect(self._on_tab_close)
        layout.addWidget(self._tabs)

        # Overview tab (always present, not closable)
        self._overview = OverviewCanvas(self)
        self._tabs.addTab(self._overview, "● Overview")
        from PySide6.QtWidgets import QTabBar
        self._tabs.tabBar().setTabButton(0, QTabBar.ButtonPosition.RightSide, None)

        # Detach button
        detach_btn = QToolButton()
        detach_btn.setText("↗")
        detach_btn.setToolTip("Detach current tab to new window (Ctrl+D)")
        detach_btn.clicked.connect(self._on_detach_current)
        self._tabs.setCornerWidget(detach_btn, Qt.TopRightCorner)

        self._detached_windows: list = []

    # ── Public API ─────────────────────────────────────────────────

    def refresh_current(self, distributions: dict, opts: dict = None):
        """Refresh whatever tab is currently active."""
        w = self._tabs.currentWidget()
        opts = opts or {}

        if isinstance(w, OverviewCanvas):
            w.render(distributions)
        elif isinstance(w, DetailCanvas):
            w.update_opts(opts)
            w.render(distributions, w._x_key, w._y_key, w._x_label, w._y_label, w._title)
        elif isinstance(w, CustomCanvas):
            w.render(distributions, opts.get("x_axis", "x"), opts.get("y_axis", "xp"), opts)
        else:
            # Also update overview in background
            self._overview.render(distributions)

    def update_overview(self, distributions: dict, opts: dict = None):
        self._overview.render(distributions, opts or {})

    def open_detail(self, x_key: str, y_key: str,
                    x_label: str, y_label: str,
                    distributions: dict, title: str = "", opts: dict = None):
        detail = DetailCanvas(self)
        detail.render(distributions, x_key, y_key, x_label, y_label, title, opts or {})
        idx = self._tabs.addTab(detail, title or f"{x_key}-{y_key}")
        self._tabs.setCurrentIndex(idx)

    def open_custom(self, distributions: dict, opts: dict = None):
        opts = opts or {}
        for i in range(self._tabs.count()):
            if isinstance(self._tabs.widget(i), CustomCanvas):
                c = self._tabs.widget(i)
                c.render(distributions, opts.get("x_axis", "x"), opts.get("y_axis", "xp"), opts)
                self._tabs.setCurrentIndex(i)
                return
        custom = CustomCanvas(self)
        custom.render(distributions, opts.get("x_axis", "x"), opts.get("y_axis", "xp"), opts)
        idx = self._tabs.addTab(custom, "Custom")
        self._tabs.setCurrentIndex(idx)

    def open_slice_plot(self, slice_result, label: str):
        """Open a slice analysis tab."""
        from beamscope.gui.detail_canvas import SliceCanvas
        canvas = SliceCanvas(self)
        canvas.render(slice_result, label)
        idx = self._tabs.addTab(canvas, f"Slices — {label}")
        self._tabs.setCurrentIndex(idx)

    def open_bff_plot(self, bff_result, label: str):
        """Open a Bunch Form Factor tab."""
        from beamscope.gui.detail_canvas import BFFCanvas
        canvas = BFFCanvas(self)
        canvas.render(bff_result, label)
        idx = self._tabs.addTab(canvas, f"BFF — {label}")
        self._tabs.setCurrentIndex(idx)

    def load_emit_set(self, sim_set):
        """Load a SimSet: create Emit/Sigma/Ref/Log tabs."""
        from beamscope.io.astra_emit import SimSet
        from beamscope.gui.detail_canvas import EmitCanvas, LogCanvas

        # Remove old emit tabs (keep Overview)
        for i in range(self._tabs.count() - 1, 0, -1):
            w = self._tabs.widget(i)
            if isinstance(w, (EmitCanvas, LogCanvas)):
                self._tabs.removeTab(i)

        if sim_set.emit is not None:
            # Emit dashboard
            canvas = EmitCanvas(self)
            canvas.render_emit_dashboard(sim_set.emit, sim_set.sigma)
            self._tabs.addTab(canvas, "Emit Dashboard")
            self._tabs.setCurrentIndex(self._tabs.count() - 1)

            # Individual emit plots
            for name, plane_key in [("Envelope", "envelope"),
                                     ("Emittance", "emittance"),
                                     ("Energy", "energy"),
                                     ("Transverse", "transverse")]:
                c = EmitCanvas(self)
                c.render_emit_plot(sim_set.emit, plane_key)
                self._tabs.addTab(c, name)

        if sim_set.sigma is not None:
            canvas = EmitCanvas(self)
            canvas.render_eigen(sim_set.sigma)
            self._tabs.addTab(canvas, "Eigen-Emittances")

        if sim_set.ref is not None:
            canvas = EmitCanvas(self)
            canvas.render_ref(sim_set.ref)
            self._tabs.addTab(canvas, "Ref Trajectory")

        if sim_set.log_text is not None:
            canvas = LogCanvas(self)
            canvas.render(sim_set.log_text)
            self._tabs.addTab(canvas, "Log")

    def current_canvas(self):
        w = self._tabs.currentWidget()
        return w if w is not None else None

    def clear(self):
        while self._tabs.count() > 1:
            self._tabs.removeTab(1)
        self._overview.clear()

    def detach_current_tab(self):
        idx = self._tabs.currentIndex()
        if idx == 0:
            return
        w = self._tabs.widget(idx)
        title = self._tabs.tabText(idx)
        self._tabs.removeTab(idx)

        from PySide6.QtWidgets import QMainWindow, QVBoxLayout, QWidget
        win = QMainWindow(self)
        win.setWindowTitle(f"{title} — beamscope")
        win.resize(800, 700)
        container = QWidget()
        lo = QVBoxLayout(container); lo.setContentsMargins(0, 0, 0, 0); lo.addWidget(w)
        win.setCentralWidget(container)
        win.setAttribute(Qt.WA_DeleteOnClose)
        win.show()
        self._detached_windows.append(win)

    def close_current_tab(self):
        """Close the currently active tab (if not the Overview tab)."""
        idx = self._tabs.currentIndex()
        if idx > 0:  # tab 0 is Overview, not closable
            self._tabs.removeTab(idx)

    # ── Slots ──────────────────────────────────────────────────────

    def _on_tab_close(self, index: int):
        if index == 0:
            return
        self._tabs.removeTab(index)

    def _on_detach_current(self):
        self.detach_current_tab()
