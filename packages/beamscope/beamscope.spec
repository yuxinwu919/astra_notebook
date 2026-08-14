# -*- mode: python ; coding: utf-8 -*-
# beamscope v0.2.0 — PyInstaller spec


a = Analysis(
    ['beamscope/gui/app.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=[
        # Qt binding
        'PySide6', 'PySide6.QtWidgets', 'PySide6.QtCore', 'PySide6.QtGui',
        # Scientific
        'scipy.special', 'scipy.stats', 'scipy.linalg',
        # Matplotlib backends
        'matplotlib.backends.backend_qtagg', 'matplotlib.backends.backend_agg',
        # beamscope analysis
        'beamscope.analysis.statistics', 'beamscope.analysis.emittance',
        'beamscope.analysis.slices', 'beamscope.analysis.bff',
        # beamscope io
        'beamscope.io', 'beamscope.io.astra', 'beamscope.io.astra_emit',
        # beamscope _plotting core (embedded astra_plotter)
        'beamscope._plotting',
        'beamscope._plotting.cosmetics', 'beamscope._plotting.plotter',
        # beamscope plot
        'beamscope.plot', 'beamscope.plot._precompute', 'beamscope.plot._artists',
        'beamscope.plot.overview', 'beamscope.plot.phase_space',
        'beamscope.plot.detail', 'beamscope.plot.distributions',
        'beamscope.plot.comparison', 'beamscope.plot.dashboard',
        'beamscope.plot.emit_plots', 'beamscope.plot.slice_plots',
        'beamscope.plot.bff_plots',
        # beamscope GUI
        'beamscope.gui.main_window', 'beamscope.gui.file_browser',
        'beamscope.gui.plot_panel', 'beamscope.gui.overview_canvas',
        'beamscope.gui.detail_canvas', 'beamscope.gui.custom_canvas',
        'beamscope.gui.properties_panel',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tensorflow', 'torch', 'pandas', 'tkinter', 'jupyter', 'notebook',
        'IPython', 'sphinx', 'pytest', 'PyQt5', 'PyQt6',
        # HoloViz / web dashboarding ecosystem
        'panel', 'bokeh', 'holoviews', 'hvplot', 'plotly', 'altair',
        'vega_datasets', 'xyzservices',
        # AWS / cloud SDK
        'botocore', 'boto3', 's3fs', 'aiobotocore',
        # JIT / ML (ocelot transitive, not used by beamscope)
        'numba', 'llvmlite', 'sklearn', 'skimage',
        # Astronomy (ocelot transitive)
        'astropy', 'astropy_iers_data', 'healpy', 'pyerfa',
        # Data processing (ocelot transitive)
        'pyarrow', 'statsmodels', 'h5py',
        # Build-time tools (should NOT be in runtime bundle)
        'mypy', 'lief', 'setuptools', 'setuptools_scm', 'wheel', 'pip',
        # Network / compute infrastructure (not needed)
        'distributed', 'dask', 'mpi4py', 'sqlalchemy',
        'aiohttp', 'uvloop', 'cryptography', 'zstandard',
    ],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='beamscope',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='beamscope',
)
app = BUNDLE(
    coll,
    name='beamscope.app',
    icon=None,
    bundle_identifier='com.beamscope.app',
)
