"""CST Particle Studio export format reader (framework).

CST can export particle distributions in various ASCII formats.
This is a placeholder for future implementation.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..distribution import Distribution

logger = logging.getLogger(__name__)


class CstReader:
    """Reader for CST Particle Studio export files.

    Currently a framework — CST export format varies by version.
    """

    format_name = "cst"

    def probe(self, path: Path) -> bool:
        """Check if file appears to be CST export format."""
        try:
            with open(path) as f:
                first_line = f.readline()
            return "CST" in first_line.upper() or "PARTICLE" in first_line.upper()
        except Exception:
            return False

    def read(self, path: Path) -> Distribution:
        """Read a CST particle export file.

        Raises:
            NotImplementedError: CST parsing not yet implemented.
        """
        raise NotImplementedError(
            "CST reader is not yet implemented. "
            "Export particles as ASCII from CST and convert using a custom script."
        )
