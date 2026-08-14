"""Right-side properties and controls panel."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QComboBox, QSpinBox,
    QCheckBox, QPushButton, QGroupBox, QLineEdit, QHBoxLayout, QFrame,
)
from PySide6.QtCore import Signal


class PropertiesPanel(QWidget):
    """Control panel for plot parameters, style, and export."""

    plot_requested = Signal()
    view_mode_changed = Signal(str)  # 'overview', 'detail', 'custom'
    export_requested = Signal(str)

    def __init__(self, main_window=None, parent=None):
        super().__init__(parent)
        self._main_window = main_window

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # ── View Mode ──
        view_group = QGroupBox("View Mode")
        vlayout = QVBoxLayout(view_group)
        self._view_mode = QComboBox()
        self._view_mode.addItems(["Overview", "Custom"])
        self._view_mode.currentTextChanged.connect(
            lambda t: self.view_mode_changed.emit(t.lower()))
        vlayout.addWidget(QLabel("Mode:"))
        vlayout.addWidget(self._view_mode)
        vlayout.addWidget(QLabel("(Detail: click subplot in Overview)",
                                 styleSheet="color: gray; font-size: 10px;"))
        layout.addWidget(view_group)

        # ── Axis Mapping (Custom mode) ──
        axis_group = QGroupBox("Axis Mapping (Custom mode)")
        alayout = QVBoxLayout(axis_group)

        vars_list = ["x", "y", "z", "xp", "yp", "dp", "E", "clock", "r", "rp", "charge"]
        self._x_axis = QComboBox()
        self._x_axis.addItems(vars_list)
        self._x_axis.setCurrentText("x")
        self._y_axis = QComboBox()
        self._y_axis.addItems(vars_list)
        self._y_axis.setCurrentText("xp")

        # Auto-update Custom tab when axis changes
        self._x_axis.currentTextChanged.connect(lambda: self._on_axis_changed())
        self._y_axis.currentTextChanged.connect(lambda: self._on_axis_changed())

        alayout.addWidget(QLabel("X Axis:"))
        alayout.addWidget(self._x_axis)
        alayout.addWidget(QLabel("Y Axis:"))
        alayout.addWidget(self._y_axis)
        layout.addWidget(axis_group)

        # ── Style ──
        style_group = QGroupBox("Style")
        slayout = QVBoxLayout(style_group)

        cmap_layout = QHBoxLayout()
        cmap_layout.addWidget(QLabel("Colormap:"))
        self._cmap = QComboBox()
        self._cmap.addItems([
            "viridis", "inferno", "plasma", "magma",
            "Blues", "Reds", "Greens", "hot", "jet", "turbo",
            "slac_desy_beam",
        ])
        self._cmap.setCurrentText("viridis")
        cmap_layout.addWidget(self._cmap)
        slayout.addLayout(cmap_layout)

        bins_layout = QHBoxLayout()
        bins_layout.addWidget(QLabel("Bins:"))
        self._bins = QSpinBox()
        self._bins.setRange(10, 300)
        self._bins.setValue(60)
        bins_layout.addWidget(self._bins)
        slayout.addLayout(bins_layout)

        self._density = QCheckBox("Density (hist2d)")
        self._density.setChecked(True)
        slayout.addWidget(self._density)

        self._scatter = QCheckBox("Scatter overlay")
        self._scatter.setChecked(False)
        slayout.addWidget(self._scatter)

        self._contour = QCheckBox("Contour overlay")
        self._contour.setChecked(False)
        slayout.addWidget(self._contour)

        self._ellipse = QCheckBox("RMS Emittance Ellipse")
        self._ellipse.setChecked(True)
        slayout.addWidget(self._ellipse)

        self._marginals = QCheckBox("Marginal Projections")
        self._marginals.setChecked(True)
        slayout.addWidget(self._marginals)

        layout.addWidget(style_group)

        # ── Reference ──
        ref_group = QGroupBox("Reference")
        rlayout = QVBoxLayout(ref_group)
        rlayout.addWidget(QLabel("Ref Energy [MeV]:"))
        self._ref_energy = QLineEdit()
        self._ref_energy.setPlaceholderText("auto")
        rlayout.addWidget(self._ref_energy)
        self._overlay_stats = QCheckBox("Overlay Stats Text")
        self._overlay_stats.setChecked(True)
        rlayout.addWidget(self._overlay_stats)
        layout.addWidget(ref_group)

        # ── Buttons ──
        update_btn = QPushButton("Update Plot")
        update_btn.clicked.connect(self.plot_requested.emit)
        layout.addWidget(update_btn)

        reset_btn = QPushButton("Reset Defaults")
        reset_btn.clicked.connect(self._reset_defaults)
        layout.addWidget(reset_btn)

        # ── Separator ──
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        layout.addWidget(line)

        # ── Export ──
        export_group = QGroupBox("Export")
        elayout = QVBoxLayout(export_group)
        save_fig_btn = QPushButton("Save Figure")
        save_fig_btn.clicked.connect(lambda: self.export_requested.emit("figure"))
        elayout.addWidget(save_fig_btn)
        export_csv_btn = QPushButton("Export CSV")
        export_csv_btn.clicked.connect(lambda: self.export_requested.emit("csv"))
        elayout.addWidget(export_csv_btn)
        export_npz_btn = QPushButton("Export NPZ")
        export_npz_btn.clicked.connect(lambda: self.export_requested.emit("npz"))
        elayout.addWidget(export_npz_btn)
        layout.addWidget(export_group)

        layout.addStretch()

    # ── Public ───────────────────────────────────────────────────

    def get_opts(self) -> dict:
        return {
            "bins": self._bins.value(),
            "cmap": self._cmap.currentText(),
            "density": self._density.isChecked(),
            "scatter": self._scatter.isChecked(),
            "contour": self._contour.isChecked(),
            "ellipse": self._ellipse.isChecked(),
            "marginals": self._marginals.isChecked(),
            "ref_energy": self._ref_energy.text() or None,
            "view_mode": self._view_mode.currentText().lower(),
            "x_axis": self._x_axis.currentText(),
            "y_axis": self._y_axis.currentText(),
        }

    def _on_axis_changed(self):
        """When axis dropdowns change in Custom mode, auto-update."""
        if self._view_mode.currentText() == "Custom":
            self.plot_requested.emit()

    def _reset_defaults(self):
        self._bins.setValue(60)
        self._cmap.setCurrentText("viridis")
        self._density.setChecked(True)
        self._scatter.setChecked(False)
        self._contour.setChecked(False)
        self._ellipse.setChecked(True)
        self._marginals.setChecked(True)
        self._ref_energy.clear()
        self._overlay_stats.setChecked(True)
        self.plot_requested.emit()
