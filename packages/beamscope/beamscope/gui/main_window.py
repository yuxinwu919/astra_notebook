"""Main window with three-panel layout: file browser, plot area, properties."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QMainWindow, QSplitter, QDockWidget, QStatusBar,
    QWidget, QVBoxLayout, QMessageBox, QFileDialog, QTextEdit,
    QApplication,
)
from PySide6.QtCore import Qt, QSettings
from PySide6.QtGui import QAction, QKeySequence

from beamscope.gui.file_browser import FileBrowser
from beamscope.gui.plot_panel import PlotPanel
from beamscope.gui.properties_panel import PropertiesPanel


class MainWindow(QMainWindow):
    """beamscope main application window."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("beamscope — Particle Beam Visualization")
        self.setMinimumSize(1000, 600)
        self.setAcceptDrops(True)  # Enable drag-drop file loading

        self._loaded_distributions: dict[str, dict] = {}  # label -> {dist, color}
        self._active_label: str | None = None

        # ── Central: Plot Panel + Log ──
        self._plot_panel = PlotPanel(self)
        self._log_widget = self._create_log_widget()

        central_splitter = QSplitter(Qt.Vertical)
        central_splitter.addWidget(self._plot_panel)
        central_splitter.addWidget(self._log_widget)
        central_splitter.setStretchFactor(0, 4)
        central_splitter.setStretchFactor(1, 1)
        self.setCentralWidget(central_splitter)

        # ── Dock widgets ──
        self._file_browser = FileBrowser(self)
        self._properties_panel = PropertiesPanel(self)

        self._file_dock = QDockWidget("File Browser", self)
        self._file_dock.setWidget(self._file_browser)
        self._file_dock.setObjectName("FileBrowserDock")
        self.addDockWidget(Qt.LeftDockWidgetArea, self._file_dock)

        self._props_dock = QDockWidget("Properties", self)
        self._props_dock.setWidget(self._properties_panel)
        self._props_dock.setObjectName("PropertiesDock")
        self.addDockWidget(Qt.RightDockWidgetArea, self._props_dock)

        # ── Settings ── (must be before _create_menus)
        self._settings = QSettings("beamscope", "beamscope")

        # ── Menu ──
        self._create_menus()

        # ── Status bar ──
        self._status_bar = QStatusBar()
        self.setStatusBar(self._status_bar)
        self._status_bar.showMessage("Ready  |  File → Open to load a distribution")

        # ── Signals ──
        self._file_browser.file_loaded.connect(self._on_file_loaded)
        self._file_browser.sim_set_loaded.connect(self._on_sim_set_loaded)
        self._file_browser.active_file_changed.connect(self._on_active_file_changed)
        self._file_browser.files_cleared.connect(self._on_files_cleared)
        self._properties_panel.plot_requested.connect(self._on_update_plot)
        self._properties_panel.view_mode_changed.connect(self._on_view_mode_changed)
        self._properties_panel.export_requested.connect(self._on_export)

        self._file_dock.visibilityChanged.connect(self._on_dock_visibility_changed)
        self._props_dock.visibilityChanged.connect(self._on_dock_visibility_changed)

        # ── Restore saved state ──
        self._restore_settings()

        # ── Auto-load on startup ──
        from PySide6.QtCore import QTimer
        QTimer.singleShot(100, self._auto_load_startup)

    # ── Log ────────────────────────────────────────────────────────

    def _create_log_widget(self):
        log = QTextEdit()
        log.setReadOnly(True)
        log.setFontFamily("Menlo, Monaco, monospace")
        log.setFontPointSize(11)
        log.setPlaceholderText("Statistics and log output will appear here...\n"
                               "Tip: Click a file in the left panel to make it the active view.")
        return log

    def log_stats(self, text: str):
        self._log_widget.clear()
        self._log_widget.setPlainText(text)

    # ── Core logic ─────────────────────────────────────────────────

    def _get_active_dist(self):
        """Return (label, dist_info) for the currently active file, or (None, None)."""
        if self._active_label and self._active_label in self._loaded_distributions:
            return self._active_label, self._loaded_distributions[self._active_label]
        # Fallback: first loaded
        if self._loaded_distributions:
            label = list(self._loaded_distributions.keys())[0]
            return label, self._loaded_distributions[label]
        return None, None

    def _get_checked_dists(self) -> dict:
        """Return all distributions whose checkbox is checked (for comparison)."""
        checked = self._file_browser.get_selected_labels()
        if not checked and self._loaded_distributions:
            return dict(self._loaded_distributions)
        return {k: v for k, v in self._loaded_distributions.items() if k in checked}

    def _refresh_overview(self):
        """Refresh Overview with the active distribution."""
        label, info = self._get_active_dist()
        if label is None:
            return
        opts = self._properties_panel.get_opts()
        self._plot_panel.update_overview({label: info}, opts)

    def _refresh_current_tab(self):
        """Refresh the currently active tab."""
        dists = self._get_checked_dists()
        if not dists:
            return
        opts = self._properties_panel.get_opts()
        self._plot_panel.refresh_current(dists, opts)

    # ── Slots ──────────────────────────────────────────────────────

    def _on_file_loaded(self, label: str, dist, color: str):
        self._loaded_distributions[label] = {"dist": dist, "color": color}
        # active_file_changed is emitted after, which triggers _on_active_file_changed

        from beamscope.analysis.statistics import compute_statistics, print_statistics
        from io import StringIO
        import sys
        stats = compute_statistics(dist, label=label)
        old = sys.stdout; sys.stdout = buf = StringIO()
        print_statistics(stats, title=f"Beam Statistics — {label}")
        sys.stdout = old
        self.log_stats(buf.getvalue())
        self._status_bar.showMessage(
            f"Loaded: {label} | N={dist.n_particle}, active={dist.n_active} | "
            f"Click a file to switch view"
        )

    def _on_sim_set_loaded(self, base_name: str, sim_set):
        """Handle a complete simulation set loaded."""
        # Load emit/sigma/ref/Log into the plot panel
        self._plot_panel.load_emit_set(sim_set)
        self._status_bar.showMessage(
            f"Simulation set loaded: {base_name} | "
            f"Emit: {sim_set.emit is not None}, "
            f"Sigma: {sim_set.sigma is not None}, "
            f"Ref: {sim_set.ref is not None}"
        )

    def _on_active_file_changed(self, label: str):
        """Called when user clicks a file in the loaded list."""
        self._active_label = label
        self._refresh_overview()
        if label in self._loaded_distributions:
            d = self._loaded_distributions[label]["dist"]
            self._status_bar.showMessage(
                f"Active: {label} | N={d.n_particle}, active={d.n_active} | "
                f"{d.n_particle} particles"
            )

    def _on_files_cleared(self):
        self._loaded_distributions.clear()
        self._active_label = None
        self._plot_panel.clear()
        self._log_widget.clear()
        self._status_bar.showMessage("Ready  |  File → Open to load a distribution")

    def _on_update_plot(self):
        """Update Plot button clicked."""
        label, _ = self._get_active_dist()
        if label:
            self._refresh_overview()
        self._refresh_current_tab()

    def _on_view_mode_changed(self, mode: str):
        """View Mode dropdown changed in properties panel."""
        if mode == "overview":
            self._plot_panel._tabs.setCurrentIndex(0)
            self._refresh_overview()
        elif mode == "custom":
            dists = self._get_checked_dists()
            if not dists:
                return
            opts = self._properties_panel.get_opts()
            x_key = opts.get("x_axis", "x")
            y_key = opts.get("y_axis", "xp")
            self._plot_panel.open_custom(dists, opts)

    def _on_export(self, export_type: str):
        if export_type == "figure":
            self._export_figure()
        elif export_type == "csv":
            self._export_csv()
        elif export_type == "npz":
            self._export_npz()

    def _export_figure(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Figure", "", "PNG (*.png);;PDF (*.pdf);;SVG (*.svg);;EPS (*.eps)"
        )
        if not path:
            return
        canvas = self._plot_panel.current_canvas()
        if canvas is not None:
            canvas._fig.savefig(path, dpi=150, bbox_inches="tight")
        self._status_bar.showMessage(f"Figure saved: {path}")

    def _export_csv(self):
        label, info = self._get_active_dist()
        if label is None:
            QMessageBox.warning(self, "Export", "No data loaded.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Statistics CSV", "", "CSV (*.csv)")
        if not path:
            return
        from beamscope.analysis.statistics import compute_statistics
        import csv
        stats = compute_statistics(info["dist"], label=label)
        d = stats.to_dict()
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(d.keys())
            w.writerow(d.values())
        self._status_bar.showMessage(f"CSV exported: {path}")

    def _export_npz(self):
        label, info = self._get_active_dist()
        if label is None:
            QMessageBox.warning(self, "Export", "No data loaded.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Raw Data NPZ", "", "NPZ (*.npz)")
        if not path:
            return
        import numpy as np
        dist = info["dist"]; m = dist.active
        np.savez(path, x=dist.x[m], y=dist.y[m], z=dist.z[m],
                 px=dist.px[m], py=dist.py[m], pz=dist.pz[m],
                 clock=dist.clock[m], charge=dist.charge[m], status=dist.status[m])
        self._status_bar.showMessage(f"Data exported: {path}")

    # ── Menu ───────────────────────────────────────────────────────

    def _create_menus(self):
        mb = self.menuBar()

        # ── File ──
        fm = mb.addMenu("&File")
        a = QAction("&Open File...", self); a.setShortcut(QKeySequence.Open)
        a.triggered.connect(self._file_browser._on_load_clicked); fm.addAction(a)
        a = QAction("Open &Directory...", self); a.setShortcut(QKeySequence("Ctrl+Shift+O"))
        a.triggered.connect(self._file_browser._on_batch_clicked); fm.addAction(a)
        fm.addSeparator()
        self._recent_menu = fm.addMenu("&Recent Files")
        self._recent_actions: list[QAction] = []
        self._update_recent_menu()
        fm.addSeparator()
        a = QAction("&Quit", self); a.setShortcut(QKeySequence.Quit)
        a.triggered.connect(self.close); fm.addAction(a)

        # ── View ──
        vm = mb.addMenu("&View")
        self._view_file_action = QAction("File &Browser", self, checkable=True, checked=True)
        self._view_file_action.triggered.connect(lambda c: self._file_dock.setVisible(c))
        vm.addAction(self._view_file_action)
        self._view_props_action = QAction("&Properties", self, checkable=True, checked=True)
        self._view_props_action.triggered.connect(lambda c: self._props_dock.setVisible(c))
        vm.addAction(self._view_props_action)
        self._view_log_action = QAction("&Log Output", self, checkable=True, checked=True)
        self._view_log_action.triggered.connect(lambda c: self._log_widget.setVisible(c))
        vm.addAction(self._view_log_action)
        vm.addSeparator()
        a = QAction("&Dark Theme", self, checkable=True)
        a.triggered.connect(self._toggle_theme)
        vm.addAction(a)
        self._dark_theme_action = a
        a = QAction("&Full Screen", self); a.setShortcut(QKeySequence("F11"))
        a.triggered.connect(lambda: self.showFullScreen() if not self.isFullScreen() else self.showNormal())
        vm.addAction(a)

        # ── Analysis ──
        am = mb.addMenu("&Analysis")
        a = QAction("&Beam Statistics", self); a.setShortcut(QKeySequence("Ctrl+T"))
        a.triggered.connect(self._on_analyze_stats); am.addAction(a)
        a = QAction("&Slice Analysis...", self); a.triggered.connect(self._on_analyze_slices)
        am.addAction(a)
        a = QAction("&Compare Checked Files", self); a.setShortcut(QKeySequence("Ctrl+Shift+C"))
        a.triggered.connect(self._on_compare_checked); am.addAction(a)
        am.addSeparator()
        a = QAction("&Bunch Form Factor", self); a.setShortcut(QKeySequence("Ctrl+B"))
        a.triggered.connect(self._on_analyze_bff); am.addAction(a)

        # ── Export ──
        em = mb.addMenu("&Export")
        a = QAction("Save &Figure...", self); a.setShortcut(QKeySequence("Ctrl+S"))
        a.triggered.connect(lambda: self._on_export("figure")); em.addAction(a)
        a = QAction("Export Statistics &CSV...", self)
        a.triggered.connect(lambda: self._on_export("csv")); em.addAction(a)
        a = QAction("Export Raw Data &NPZ...", self)
        a.triggered.connect(lambda: self._on_export("npz")); em.addAction(a)

        # ── Window ──
        wm = mb.addMenu("&Window")
        a = QAction("Close Current &Tab", self); a.setShortcut(QKeySequence("Ctrl+W"))
        a.triggered.connect(self._close_current_tab); wm.addAction(a)
        a = QAction("&Detach Current Tab", self); a.setShortcut(QKeySequence("Ctrl+D"))
        a.triggered.connect(self._plot_panel.detach_current_tab); wm.addAction(a)

        # ── Help ──
        hm = mb.addMenu("&Help")
        a = QAction("&About beamscope", self); a.triggered.connect(self._show_about)
        hm.addAction(a)

    def _close_current_tab(self):
        """Close the currently active tab (if closable)."""
        self._plot_panel.close_current_tab()

    def _add_recent_file(self, path: str):
        """Add a file path to the recent files list (persisted in QSettings)."""
        recent = self._settings.value("recent_files", []) or []
        if path in recent:
            recent.remove(path)
        recent.insert(0, path)
        recent = recent[:10]  # keep last 10
        self._settings.setValue("recent_files", recent)
        self._update_recent_menu()

    def _update_recent_menu(self):
        """Rebuild the Recent Files submenu from QSettings."""
        self._recent_menu.clear()
        self._recent_actions.clear()
        recent = self._settings.value("recent_files", []) or []
        if not recent:
            a = QAction("(none)", self); a.setEnabled(False)
            self._recent_menu.addAction(a)
            return
        for p in recent:
            a = QAction(p, self)
            a.triggered.connect(lambda checked, path=p: self._open_recent(path))
            self._recent_menu.addAction(a)
            self._recent_actions.append(a)
        self._recent_menu.addSeparator()
        a = QAction("Clear Recent Files", self)
        a.triggered.connect(self._clear_recent_files)
        self._recent_menu.addAction(a)

    def _open_recent(self, path: str):
        """Open a file from the recent files list."""
        from pathlib import Path
        self._file_browser._load_file(Path(path))
        self._add_recent_file(path)

    def _clear_recent_files(self):
        self._settings.setValue("recent_files", [])
        self._update_recent_menu()

    def _on_dock_visibility_changed(self, visible):
        s = self.sender()
        if s == self._file_dock: self._view_file_action.setChecked(visible)
        elif s == self._props_dock: self._view_props_action.setChecked(visible)

    def _on_analyze_stats(self):
        label, info = self._get_active_dist()
        if label is None: return
        from beamscope.analysis.statistics import compute_statistics, print_statistics
        from io import StringIO; import sys
        old = sys.stdout; sys.stdout = buf = StringIO()
        print_statistics(compute_statistics(info["dist"], label=label))
        sys.stdout = old; self.log_stats(buf.getvalue())

    def _on_analyze_slices(self):
        label, info = self._get_active_dist()
        if label is None: return
        from PySide6.QtWidgets import QInputDialog
        n, ok = QInputDialog.getInt(self, "Slice Analysis", "Number of slices:", 20, 3, 200)
        if not ok: return
        from beamscope.analysis.slices import compute_slice_analysis
        slices = compute_slice_analysis(info["dist"], n_slices=n)
        self._plot_panel.open_slice_plot(slices, label)

    def _on_analyze_bff(self):
        label, info = self._get_active_dist()
        if label is None: return
        import numpy as np
        from beamscope.analysis.bff import compute_bff
        dist = info["dist"]
        active = dist.active
        z_data = dist.z[active] - np.mean(dist.z[active])
        charge_arr = dist.charge[active]
        bff_result = compute_bff(z_data, charge_arr)
        self._plot_panel.open_bff_plot(bff_result, label)

    def _on_compare_checked(self):
        dists = self._get_checked_dists()
        if len(dists) < 2:
            QMessageBox.information(self, "Compare", "Check at least 2 files for comparison.")
            return
        from beamscope.plot.comparison import plot_comparison
        import matplotlib.pyplot as plt
        fig = plot_comparison(
            {l: d["dist"] for l, d in dists.items()},
            plot_type="statistics"
        )
        # Show in a detached window
        from PySide6.QtWidgets import QMainWindow as QMW
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        win = QMW(self); win.setWindowTitle("Comparison — Statistics")
        canvas = FigureCanvasQTAgg(fig); win.setCentralWidget(canvas)
        win.resize(900, 600); win.setAttribute(Qt.WA_DeleteOnClose); win.show()

    def _toggle_theme(self, dark: bool):
        """Toggle between light and dark themes."""
        if dark:
            self.setStyleSheet("""
                QMainWindow { background-color: #2b2b2b; }
                QDockWidget { background-color: #333; color: #eee; }
                QLabel { color: #eee; }
                QGroupBox { color: #ddd; border: 1px solid #555; margin-top: 8px; padding-top: 8px; }
                QGroupBox::title { color: #ddd; }
                QPushButton { background-color: #444; color: #eee; border: 1px solid #555; padding: 4px 12px; }
                QPushButton:hover { background-color: #555; }
                QComboBox { background-color: #444; color: #eee; border: 1px solid #555; }
                QSpinBox { background-color: #444; color: #eee; border: 1px solid #555; }
                QCheckBox { color: #eee; }
                QLineEdit { background-color: #444; color: #eee; border: 1px solid #555; }
                QTabWidget::pane { background-color: #2b2b2b; }
                QTabBar::tab { background-color: #333; color: #aaa; padding: 6px 12px; }
                QTabBar::tab:selected { background-color: #444; color: #fff; }
                QTextEdit { background-color: #1e1e1e; color: #ddd; }
                QListWidget { background-color: #333; color: #eee; }
                QTreeView { background-color: #333; color: #eee; }
                QStatusBar { background-color: #222; color: #aaa; }
                QMenuBar { background-color: #333; color: #eee; }
                QMenuBar::item:selected { background-color: #555; }
                QMenu { background-color: #333; color: #eee; }
                QMenu::item:selected { background-color: #555; }
                QSplitter::handle { background-color: #555; }
            """)
        else:
            self.setStyleSheet("")

    def _show_about(self):
        from beamscope import __version__ as version
        QMessageBox.about(self, "About beamscope",
            f"<b>beamscope v{version}</b><br><br>"
            f"Accelerator particle distribution visualization toolkit.<br>"
            f"Six-dimensional phase space analysis for ASTRA simulations.<br><br>"
            f"Built with PySide6 + Matplotlib + NumPy + SciPy.<br>"
            f"SLAC-DESY beam colormap: C. Behrens (DESY) / T. Maxwell (SLAC).<br><br>"
            f"<a href='https://github.com'>GitHub</a> | MIT License")

    # ── Settings ───────────────────────────────────────────────────

    def _restore_settings(self):
        g = self._settings.value("window/geometry")
        if g: self.restoreGeometry(g)
        s = self._settings.value("window/state")
        if s: self.restoreState(s)

    def _auto_load_startup(self):
        """Auto-load simulation files on startup."""
        loaded = self._file_browser.auto_load_simulation()
        if loaded:
            self._status_bar.showMessage(f"Auto-loaded {loaded} files from simulation_files/")
        else:
            self._status_bar.showMessage(
                "Ready  |  File → Open to load a distribution  |  "
                "Tip: place .ini files in simulation_files/ for auto-load"
            )

    # ── Drag & Drop ───────────────────────────────────────────────

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            path = Path(url.toLocalFile())
            if path.is_dir():
                self._file_browser._on_batch_clicked()
            elif path.suffix.lower() in (".ini", ".dat", ".001", ".inp", ".ast"):
                self._file_browser._load_file(path)
            # Emit/sigma/ref/log files — route through _load_file which detects type
            elif len(path.suffixes) >= 2 and path.suffix.lstrip(".").isdigit():
                self._file_browser._load_file(path)
        event.acceptProposedAction()

    # ── Settings ───────────────────────────────────────────────────

    def closeEvent(self, event):
        self._settings.setValue("window/geometry", self.saveGeometry())
        self._settings.setValue("window/state", self.saveState())
        super().closeEvent(event)
