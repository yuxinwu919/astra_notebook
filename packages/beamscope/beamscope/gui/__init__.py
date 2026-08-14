"""beamscope GUI — PySide6-based desktop application.

Usage:
    python -m beamscope.gui.app
    # or
    from beamscope.gui import main; main()
"""

__all__ = ["main"]


def __getattr__(name: str):
    """Lazy import to avoid loading PySide6 + matplotlib QtAgg until needed."""
    if name == "main":
        from beamscope.gui.app import main as _main
        return _main
    raise AttributeError(f"module 'beamscope.gui' has no attribute '{name}'")
