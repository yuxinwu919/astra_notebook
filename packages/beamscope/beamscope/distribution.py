"""Unified particle distribution data model.

All I/O readers convert their native formats into this standard representation.
Analysis and plotting functions operate on this type exclusively.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np


@dataclass
class Distribution:
    """Standardized particle bunch distribution in 6D phase space.

    All array attributes have shape (N,) for N particles.
    Coordinate conventions follow ASTRA:
      - x, y    : transverse position [m]
      - z       : longitudinal position [m], relative to reference particle
      - px, py  : transverse momentum [eV/c]
      - pz      : longitudinal momentum [eV/c]
      - clock   : arrival time [s] (ASTRA convention: stored in ns, converted to s)
      - charge  : macro-particle charge [nC]
      - status  : particle status flag (0=active, >0=lost, <0=not started)

    The ``active`` property returns a boolean mask for status == 0 particles.
    All analysis functions filter to active particles automatically.
    """

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    px: np.ndarray
    py: np.ndarray
    pz: np.ndarray
    clock: np.ndarray
    charge: np.ndarray
    status: np.ndarray

    # Metadata (ASTRA Manual V3.2 Table 1 header)
    ref_time_ns: float = 0.0       # header[0]: reference time [ns]
    ref_momentum_eVc: float = 0.0  # header[1]: reference momentum p_ref [eV/c]
    total_charge_nC: float = 0.0   # header[2]: total bunch charge [nC]
    ref_x_m: float = 0.0           # header[3]: reference x offset [m]
    ref_y_m: float = 0.0           # header[4]: reference y offset [m]
    source: str = ""
    format: str = ""  # e.g. 'astra_binary', 'astra_ascii', 'elegant_sdds'
    attrs: dict[str, Any] = field(default_factory=dict)

    # Backward-compatible alias (deprecated — use ref_momentum_eVc)
    @property
    def ref_energy_eV(self) -> float:
        """DEPRECATED: use ref_momentum_eVc.

        Historically misnamed — ASTRA stores reference *momentum* (eV/c),
        not kinetic energy, in header[1]. Kept for backward compatibility.
        """
        return self.ref_momentum_eVc

    @ref_energy_eV.setter
    def ref_energy_eV(self, value: float) -> None:
        self.ref_momentum_eVc = value

    # ── Derived properties ────────────────────────────────────────

    @property
    def active(self) -> np.ndarray:
        """Boolean mask: True for active (status == 0) particles."""
        return self.status == 0

    @property
    def not_started(self) -> np.ndarray:
        """Boolean mask: True for not-yet-started (status < 0) particles.

        ASTRA Manual V3.2: status_flag < 0 means particle not yet started
        (e.g., cathode electrons not yet emitted at this z-position).
        """
        return self.status < 0

    @property
    def lost(self) -> np.ndarray:
        """Boolean mask: True for lost (status > 0) particles.

        ASTRA Manual V3.2: status_flag > 0 encodes loss reason.
        """
        return self.status > 0

    @property
    def reference_kinetic_energy_eV(self) -> float:
        """Reference kinetic energy [eV], correctly computed from momentum.

        E_kin = sqrt(p_ref² + m_e²) - m_e

        This is the correct conversion — unlike the old code which
        treated momentum as kinetic energy directly.
        """
        if self.ref_momentum_eVc <= 0:
            return 0.0
        from .constants import kinetic_energy_from_momentum
        return kinetic_energy_from_momentum(self.ref_momentum_eVc)

    @property
    def n_particle(self) -> int:
        """Total number of particles (including lost)."""
        return len(self.x)

    @property
    def n_active(self) -> int:
        """Number of active particles."""
        return int(np.sum(self.active))

    @property
    def active_charge_nC(self) -> float:
        """Total charge of active particles [nC]."""
        return float(np.sum(self.charge[self.active]))

    @property
    def mean_energy_eV(self) -> float:
        """Mean longitudinal momentum of active particles [eV/c]."""
        return float(np.mean(self.pz[self.active]))

    @property
    def is_valid(self) -> bool:
        """Check that all arrays have consistent lengths."""
        n = len(self.x)
        return all(
            len(arr) == n
            for arr in [self.y, self.z, self.px, self.py, self.pz, self.clock, self.charge, self.status]
        )

    # ── Factory methods ────────────────────────────────────────────

    @classmethod
    def from_dict(cls, data: dict[str, Any], source: str = "", fmt: str = "") -> Distribution:
        """Create Distribution from a dictionary (compatible with utils.py output).

        Args:
            data: Dictionary with keys matching the legacy utils.read_astra_distribution format.
            source: Source file path.
            fmt: Original format identifier.

        Returns:
            Distribution instance.
        """
        header = data.get("header", np.zeros(5))
        return cls(
            x=np.asarray(data.get("x", [])),
            y=np.asarray(data.get("y", [])),
            z=np.asarray(data.get("z", [])),
            px=np.asarray(data.get("px", [])),
            py=np.asarray(data.get("py", [])),
            pz=np.asarray(data.get("pz", [])),
            clock=np.asarray(data.get("clock", [])) * 1e-9,  # ns → s
            charge=np.asarray(data.get("macro_charge", [])),
            status=np.asarray(data.get("status_flag", []), dtype=np.int32),
            ref_time_ns=float(header[0]),
            ref_momentum_eVc=float(header[1]),
            total_charge_nC=float(header[2]),
            ref_x_m=float(header[3]) if len(header) > 3 else 0.0,
            ref_y_m=float(header[4]) if len(header) > 4 else 0.0,
            source=source,
            format=fmt,
        )

    @classmethod
    def from_arrays(
        cls,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        px: np.ndarray,
        py: np.ndarray,
        pz: np.ndarray,
        clock: np.ndarray,
        charge: np.ndarray,
        status: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> Distribution:
        """Create Distribution from raw arrays.

        Units: x/y/z [m], px/py/pz [eV/c], clock [s], charge [nC].
        """
        if status is None:
            status = np.zeros(len(x), dtype=np.int32)
        return cls(
            x=np.asarray(x),
            y=np.asarray(y),
            z=np.asarray(z),
            px=np.asarray(px),
            py=np.asarray(py),
            pz=np.asarray(pz),
            clock=np.asarray(clock),
            charge=np.asarray(charge),
            status=np.asarray(status, dtype=np.int32),
            **kwargs,
        )

    # ── Slicing / filtering ────────────────────────────────────────

    def filter_active(self) -> Distribution:
        """Return a new Distribution with only active particles."""
        mask = self.active
        return Distribution(
            x=self.x[mask].copy(),
            y=self.y[mask].copy(),
            z=self.z[mask].copy(),
            px=self.px[mask].copy(),
            py=self.py[mask].copy(),
            pz=self.pz[mask].copy(),
            clock=self.clock[mask].copy(),
            charge=self.charge[mask].copy(),
            status=self.status[mask].copy(),
            ref_time_ns=self.ref_time_ns,
            ref_momentum_eVc=self.ref_momentum_eVc,
            total_charge_nC=self.active_charge_nC,
            ref_x_m=self.ref_x_m,
            ref_y_m=self.ref_y_m,
            source=self.source,
            format=self.format,
            attrs=dict(self.attrs),
        )

    def sample(self, n: int, seed: int = 42) -> Distribution:
        """Randomly sample n particles (without replacement)."""
        rng = np.random.default_rng(seed)
        active_mask = self.active
        active_indices = np.where(active_mask)[0]
        n_sample = min(n, len(active_indices))
        chosen = rng.choice(active_indices, size=n_sample, replace=False)

        return Distribution(
            x=self.x[chosen].copy(),
            y=self.y[chosen].copy(),
            z=self.z[chosen].copy(),
            px=self.px[chosen].copy(),
            py=self.py[chosen].copy(),
            pz=self.pz[chosen].copy(),
            clock=self.clock[chosen].copy(),
            charge=self.charge[chosen].copy(),
            status=self.status[chosen].copy(),
            ref_time_ns=self.ref_time_ns,
            ref_momentum_eVc=self.ref_momentum_eVc,
            total_charge_nC=float(np.sum(self.charge[chosen])),
            ref_x_m=self.ref_x_m,
            ref_y_m=self.ref_y_m,
            source=f"{self.source} (sampled n={n_sample})",
            format=self.format,
            attrs=dict(self.attrs),
        )

    def __repr__(self) -> str:
        return (
            f"Distribution(n={self.n_particle}, active={self.n_active}, "
            f"Q={self.total_charge_nC:.4g} nC, "
            f"source='{Path(self.source).name if self.source else 'N/A'}')"
        )

    def summary(self) -> str:
        """Return a multi-line summary string."""
        from .constants import EV_TO_MEV

        active_charge = self.active_charge_nC
        e_kin = self.reference_kinetic_energy_eV
        return (
            f"Distribution Summary\n"
            f"  Source:         {self.source or 'N/A'}\n"
            f"  Format:         {self.format or 'N/A'}\n"
            f"  Particles:      {self.n_particle} total, "
            f"{self.n_active} active, "
            f"{int(np.sum(self.lost))} lost, "
            f"{int(np.sum(self.not_started))} not started\n"
            f"  Total charge:   {active_charge:.4f} nC\n"
            f"  Ref. time:      {self.ref_time_ns:.4f} ns\n"
            f"  Ref. momentum:  {self.ref_momentum_eVc * EV_TO_MEV:.4f} MeV/c\n"
            f"  Ref. kinetic E: {e_kin * EV_TO_MEV:.4f} MeV\n"
            f"  Ref. position:  ({self.ref_x_m*1e3:.4f}, {self.ref_y_m*1e3:.4f}) mm\n"
        )
