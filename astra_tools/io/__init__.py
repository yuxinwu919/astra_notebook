"""I/O layer: ASTRA distribution and evolution-file readers."""

from .astra_dist import AstraDistributionReader


def read_distribution(path, fmt: str = "auto"):
    """Read a particle distribution file, auto-detecting the format.

    Currently only the ASTRA format (binary or ASCII) is supported.

    Args:
        path: path to the distribution file.
        fmt: format hint; 'auto' (default) or 'astra'.

    Returns:
        Distribution instance.
    """
    from pathlib import Path

    path = Path(path)
    if fmt not in ("auto", "astra"):
        raise ValueError("unknown format " + repr(fmt))
    reader = AstraDistributionReader()
    if fmt == "astra" or reader.probe(path):
        return reader.read(path)
    raise ValueError("cannot determine format of " + str(path))


__all__ = ["read_distribution", "AstraDistributionReader"]
