"""I/O layer: readers for various simulation codes."""

from pathlib import Path
from typing import Optional

from ..distribution import Distribution
from .astra import AstraReader

# Registry of format probes
_READERS: list[type] = [AstraReader]


def read_distribution(path: Path, fmt: str = "auto") -> Distribution:
    """Read a particle distribution file, auto-detecting the format.

    Args:
        path: Path to distribution file.
        fmt: Format hint. Use 'auto' to probe all registered readers.
             Use 'astra', 'elegant', 'cst', 'echo2d' to force a specific reader.

    Returns:
        Distribution instance.

    Raises:
        ValueError: If the format cannot be determined or the file cannot be read.
    """
    if fmt != "auto":
        for reader_cls in _READERS:
            if reader_cls.format_name == fmt:
                return reader_cls().read(path)
        raise ValueError(f"Unknown format '{fmt}'. Available: {[r.format_name for r in _READERS]}")

    # Auto-detect
    errors = []
    for reader_cls in _READERS:
        reader = reader_cls()
        if reader.probe(path):
            try:
                return reader.read(path)
            except Exception as e:
                errors.append(f"{reader_cls.format_name}: {e}")
                continue

    if errors:
        raise ValueError(f"Failed to read '{path}': {'; '.join(errors)}")
    raise ValueError(f"Cannot determine format of '{path}'. Tried: {[r.format_name for r in _READERS]}")


def register_reader(reader_cls: type) -> None:
    """Register a custom reader class."""
    _READERS.append(reader_cls)
