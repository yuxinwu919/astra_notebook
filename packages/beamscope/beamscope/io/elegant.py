"""ELEGANT SDDS format reader (framework).

ELEGANT uses SDDS (Self Describing Data Sets) format for phase space output.
Full support requires the SDDS Python bindings or sdds2plain conversion.
This is a placeholder for future implementation.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..distribution import Distribution

logger = logging.getLogger(__name__)


class ElegantReader:
    """Reader for ELEGANT phase space files (SDDS format).

    Currently a framework — requires SDDS parsing infrastructure.
    """

    format_name = "elegant"

    def probe(self, path: Path) -> bool:
        """Check if file appears to be SDDS format (starts with 'SDDS')."""
        try:
            with open(path, "rb") as f:
                header = f.read(4)
            return header == b"SDDS"
        except Exception:
            return False

    def read(self, path: Path) -> Distribution:
        """Read an ELEGANT SDDS phase space file.

        Raises:
            NotImplementedError: SDDS parsing not yet implemented.
        """
        raise NotImplementedError(
            "ELEGANT SDDS reader is not yet implemented. "
            "Use sdds2plain to convert to ASCII first, then read as ASTRA format."
        )
