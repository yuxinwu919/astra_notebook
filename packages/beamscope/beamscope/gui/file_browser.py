"""File browser panel: tree view, load buttons, and loaded files list.

Interaction model:
  - Single click → select as ACTIVE file (bold highlight), updates Overview immediately
  - Checkbox → include in multi-file comparison / batch statistics
  - First loaded file is auto-selected as active
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QTreeView, QPushButton, QListWidget,
    QListWidgetItem, QFileDialog, QLabel, QHBoxLayout, QFileSystemModel,
    QAbstractItemView, QMenu,
)
from PySide6.QtCore import QDir, Signal, Qt
from PySide6.QtGui import QAction, QFont, QBrush, QColor

from beamscope.distribution import Distribution


_COLORS = [
    "#2196F3", "#F44336", "#4CAF50", "#FF9800", "#9C27B0",
    "#00BCD4", "#E91E63", "#3F51B5", "#CDDC39", "#795548",
]


class FileBrowser(QWidget):
    """Left-side panel with file system browser and loaded files list."""

    file_loaded = Signal(str, Distribution, str)      # label, dist, color
    sim_set_loaded = Signal(str, object)               # base_name, SimSet
    active_file_changed = Signal(str)                   # label of newly active file
    files_cleared = Signal()

    def __init__(self, main_window=None, parent=None):
        super().__init__(parent)
        self._main_window = main_window
        self._color_idx = 0
        self._active_label: str | None = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        # ── File system tree ──
        layout.addWidget(QLabel("File Browser", styleSheet="font-weight: bold;"))

        self._model = QFileSystemModel()
        self._model.setNameFilters([
            "*.ini", "*.dat", "*.001", "*.ast", "*.inp",
            "*.Xemit.*", "*.Yemit.*", "*.Zemit.*", "*.Sigma.*",
        ])
        self._model.setNameFilterDisables(False)

        # Smart default directory: simulation_files/ near cwd, then home
        start_dir = self._find_simulation_dir()

        self._tree = QTreeView()
        self._tree.setModel(self._model)
        self._tree.setRootIndex(self._model.index(str(start_dir)))
        self._tree.setColumnWidth(0, 180)
        for i in range(1, 4):
            self._tree.setColumnHidden(i, True)
        self._tree.setAnimated(False)
        self._tree.setIndentation(12)
        self._tree.doubleClicked.connect(self._on_tree_double_click)
        layout.addWidget(self._tree)

        # ── Buttons ──
        btn_layout = QHBoxLayout()
        for text, slot in [("Load", self._on_load_clicked),
                           ("Batch", self._on_batch_clicked),
                           ("Sim Set", self._on_load_sim_set),
                           ("Clear", self._on_clear_clicked)]:
            btn = QPushButton(text)
            btn.clicked.connect(slot)
            btn_layout.addWidget(btn)
        layout.addLayout(btn_layout)

        # ── Loaded files list ──
        hdr = QHBoxLayout()
        hdr.addWidget(QLabel("Loaded Files", styleSheet="font-weight: bold;"))
        hdr.addStretch()
        hdr.addWidget(QLabel("☑=compare", styleSheet="color: gray; font-size: 10px;"))
        layout.addLayout(hdr)

        self._list = QListWidget()
        self._list.setSelectionMode(QAbstractItemView.SingleSelection)
        self._list.setContextMenuPolicy(Qt.CustomContextMenu)
        self._list.customContextMenuRequested.connect(self._on_context_menu)
        self._list.itemClicked.connect(self._on_item_clicked)
        layout.addWidget(self._list)

        # Bold font for active item
        self._bold_font = QFont()
        self._bold_font.setBold(True)

    # ── Public ─────────────────────────────────────────────────────

    def get_selected_labels(self) -> list[str]:
        """Return labels of checked items (for comparison)."""
        labels = []
        for i in range(self._list.count()):
            item = self._list.item(i)
            if item.checkState() == Qt.Checked:
                labels.append(item.text())
        return labels

    def get_active_label(self) -> str | None:
        return self._active_label

    def _update_active_highlight(self):
        """Bold the active item, unbold others."""
        for i in range(self._list.count()):
            item = self._list.item(i)
            is_active = item.text() == self._active_label
            f = QFont()
            f.setBold(is_active)
            item.setFont(f)
            if is_active:
                item.setForeground(QBrush(QColor("#1565C0")))
            else:
                item.setForeground(QBrush(QColor(0, 0, 0)))

    # ── Slots ──────────────────────────────────────────────────────

    def _on_item_clicked(self, item: QListWidgetItem):
        """Single click → select as active file, update Overview."""
        self._active_label = item.text()
        self._update_active_highlight()
        self.active_file_changed.emit(self._active_label)

    def _on_load_clicked(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Load Distribution File(s)",
            str(Path.cwd()),
            "ASTRA Files (*.ini *.dat *.001 *.ast *.inp *.Xemit.* *.Yemit.* *.Zemit.* *.Sigma.*);;All Files (*)"
        )
        for p in paths:
            self._load_file(Path(p))

    def _on_batch_clicked(self):
        dir_path = QFileDialog.getExistingDirectory(
            self, "Select Directory", str(Path.cwd())
        )
        if not dir_path:
            return
        d = Path(dir_path)
        for pattern in ["*.ini", "*.001", "*.dat"]:
            for f in sorted(d.glob(pattern)):
                self._load_file(f)

    def _on_load_sim_set(self):
        """Load a complete simulation set (emit/sigma/ref/log) from a directory."""
        from PySide6.QtWidgets import QMessageBox
        sim_dir = Path(self._tree.model().filePath(self._tree.rootIndex()))
        try:
            from beamscope.io.astra_emit import discover_sim_files
            sim_set = discover_sim_files(sim_dir)
            base_name = sim_set.base_name
            self.sim_set_loaded.emit(base_name, sim_set)
        except Exception as e:
            QMessageBox.warning(self, "Load Error",
                f"Could not discover simulation files in:\n{sim_dir}\n\n{e}")

    def _on_clear_clicked(self):
        self._list.clear()
        self._color_idx = 0
        self._active_label = None
        self.files_cleared.emit()

    def _on_tree_double_click(self, index):
        path = Path(self._model.filePath(index))
        if path.is_file():
            self._load_file(path)

    def _on_context_menu(self, pos):
        item = self._list.itemAt(pos)
        if item is None:
            return
        menu = QMenu(self)
        menu.addAction("Set as Active", lambda: self._on_item_clicked(item))
        menu.addAction("Rename Label...", lambda: self._rename_item(item))
        menu.addSeparator()
        menu.addAction("Remove", lambda: self._remove_item(item))
        menu.exec(self._list.mapToGlobal(pos))

    def _remove_item(self, item):
        was_active = item.text() == self._active_label
        row = self._list.row(item)
        self._list.takeItem(row)
        if was_active and self._list.count() > 0:
            new_item = self._list.item(0)
            self._on_item_clicked(new_item)
        elif self._list.count() == 0:
            self._active_label = None

    def _rename_item(self, item):
        from PySide6.QtWidgets import QInputDialog
        old_name = item.text()
        new_name, ok = QInputDialog.getText(
            self, "Rename", "New label:", text=old_name
        )
        if ok and new_name:
            item.setText(new_name)
            if self._active_label == old_name:
                self._active_label = new_name

    # ── Helpers ────────────────────────────────────────────────────

    def _find_simulation_dir(self) -> Path:
        """Find a sensible default directory for the file browser.

        Priority:
          1. ./simulation_files (relative to cwd)
          2. ../simulation_files (one level up, for running from packages/beamscope/)
          3. QSettings last directory
          4. User home
        """
        cwd = Path.cwd()
        candidates = [
            cwd / "simulation_files",
            cwd.parent / "simulation_files",
        ]
        # Also check if running from PyInstaller bundle
        import sys
        if getattr(sys, 'frozen', False):
            bundle_dir = Path(sys._MEIPASS).parent.parent  # .app/Contents/MacOS/ → .app/
            candidates.append(bundle_dir.parent / "simulation_files")

        for d in candidates:
            if d.exists() and d.is_dir():
                return d

        return Path.home()

    def auto_load_simulation(self):
        """Auto-discover and load files from the current directory."""
        root_path = Path(self._model.filePath(self._tree.rootIndex()))
        sim_dir = root_path if root_path.name == "simulation_files" else root_path / "simulation_files"
        if not sim_dir.exists():
            return

        # Load distribution files
        loaded = 0
        for pattern in ["*.ini", "*.001", "*.dat"]:
            for f in sorted(sim_dir.glob(pattern)):
                if f.stat().st_size > 100:  # skip tiny/empty files
                    self._load_file(f)
                    loaded += 1

        # Load simulation set (emit/sigma/ref) if available
        try:
            from beamscope.io.astra_emit import discover_sim_files
            sim_set = discover_sim_files(sim_dir)
            self.sim_set_loaded.emit(sim_set.base_name, sim_set)
        except Exception:
            pass

        return loaded

    # File name patterns that indicate emit/sigma/ref/log files
    _EMIT_EXTENSIONS = {'.xemit', '.yemit', '.zemit', '.cemit'}
    _SIGMA_EXTENSION = '.sigma'
    _REF_EXTENSION = '.ref'
    _LOG_EXTENSION = '.log'

    @classmethod
    def _is_emit_file(cls, path: Path) -> bool:
        """Check if path is an ASTRA emit/sigma/ref/log file (not a distribution)."""
        suffixes_lower = {s.lower() for s in path.suffixes}
        return bool(suffixes_lower & (cls._EMIT_EXTENSIONS | {cls._SIGMA_EXTENSION,
                                                                cls._REF_EXTENSION,
                                                                cls._LOG_EXTENSION}))

    def _load_file(self, path: Path):
        if not path.exists():
            return

        # Route emit/sigma/ref/log files to the simulation-set pipeline
        if self._is_emit_file(path):
            self._load_emit_file(path)
            return

        from beamscope.io import read_distribution

        try:
            dist = read_distribution(path)
        except Exception as e:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Load Error", f"Failed to load {path.name}:\n{e}")
            return

        color = _COLORS[self._color_idx % len(_COLORS)]
        self._color_idx += 1

        # Z-position from ASTRA tracking filename
        label = path.stem
        parts = path.name.split(".")
        if len(parts) == 3 and parts[1].isdigit():
            z_pos = float(parts[1]) / 100
            if z_pos < 100:
                label = f"z={z_pos:.2f}m"

        # Emit before adding to list so main_window can store dist
        self.file_loaded.emit(label, dist, color)

        item = QListWidgetItem(label)
        item.setCheckState(Qt.Checked)
        item.setData(Qt.UserRole, dist)  # store ref
        self._list.addItem(item)

        # First loaded file → auto-select as active
        if self._active_label is None:
            self._active_label = label
        self._update_active_highlight()
        self.active_file_changed.emit(self._active_label)

        # Track in recent files
        if self._main_window is not None:
            self._main_window._add_recent_file(str(path.resolve()))

    def _load_emit_file(self, path: Path):
        """Load an emit/sigma/ref/log file via the simulation-set pipeline.

        Extracts the base name and run number from the filename, then uses
        ``discover_sim_files`` to load the entire simulation set (emit, sigma,
        ref, log). Results are forwarded to the main window for display in
        Emit/Sigma/Ref/Log tabs.
        """
        from beamscope.io.astra_emit import discover_sim_files

        parent_dir = path.parent

        # Try to determine base name from the filename
        # e.g., 'Example.Xemit.001' → base='Example'
        name = path.name
        suffixes_lower = {s.lower() for s in path.suffixes}
        base_candidates = []
        for ext in self._EMIT_EXTENSIONS | {self._SIGMA_EXTENSION,
                                              self._REF_EXTENSION,
                                              self._LOG_EXTENSION}:
            lower_name = name.lower()
            idx = lower_name.find(ext)
            if idx > 0:
                base_candidates.append(name[:idx])

        if not base_candidates:
            # Fallback: use stem of the first part
            base_candidates.append(path.name.split(".")[0])

        base_name = min(base_candidates, key=len)  # shortest match is most likely

        try:
            sim_set = discover_sim_files(parent_dir, rootname=base_name)
            self.sim_set_loaded.emit(base_name, sim_set)
        except Exception as e:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(
                self, "Load Error",
                f"Failed to load simulation set from {path.name}:\n{e}"
            )
