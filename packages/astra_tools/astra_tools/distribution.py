"""Unified particle distribution data model.

All I/O readers convert native formats into this representation and all
analysis/plotting functions operate on it exclusively.

Coordinate and unit conventions follow the ASTRA Manual V3.2, Table 1:

    x, y, z   : position [m], z relative to the reference particle
    px, py, pz: momentum [eV/c]
    clock     : arrival time [s] (file stores ns; converted on read)
    charge    : macro-particle charge [nC]
    index     : particle species index (1=e-, 2=e+, 3=p, 4=H+); 0 if unknown
    status    : status flag, see Table 2 and the mask properties below

Header (5 values, binary files only):
    ref_time_ns, ref_momentum_eVc, total_charge_nC, ref_x_m, ref_y_m

Status semantics per ASTRA Manual V3.2 (Table 2 + section 4.13):
    -6 .. -1 : at the cathode, not yet started
    < -6     : lost (aperture, backwards, discarded, ...)
    0, 1     : passive (probe) - tracked but EXCLUDED from statistics
    > 1      : taken into account for emittance / statistics (manual 4.13)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np

from .constants import (
    kinetic_energy_from_momentum,
    gamma_from_momentum,
    beta_from_gamma,
)


@dataclass
class Distribution:
    """Standardized particle bunch distribution in 6D phase space."""

    x: np.ndarray      # [m]
    y: np.ndarray      # [m]
    z: np.ndarray      # [m], relative to reference particle
    px: np.ndarray     # [eV/c]
    py: np.ndarray     # [eV/c]
    pz: np.ndarray     # [eV/c]
    clock: np.ndarray  # [s]
    charge: np.ndarray  # [nC]
    status: np.ndarray  # int32 status flag
    index: Optional[np.ndarray] = None  # species index (10-col files)

    # Header / metadata (ASTRA Manual V3.2 Table 1)
    ref_time_ns: float = 0.0
    ref_momentum_eVc: float = 0.0
    total_charge_nC: float = 0.0
    ref_x_m: float = 0.0
    ref_y_m: float = 0.0
    ref_z_m: float = 0.0   # absolute reference z; 0 when unknown (binary)
    source: str = ""
    format: str = ""
    attrs: dict[str, Any] = field(default_factory=dict)

    # -- Status masks (ASTRA Manual V3.2, Table 2 / section 4.13) ------

    @property
    def active(self) -> np.ndarray:
        """Mask of particles used for statistics: status > 1 (manual 4.13)."""
        return self.status > 1

    @property
    def passive(self) -> np.ndarray:
        """Passive (probe) particles: status 0 or 1. Tracked but excluded."""
        return (self.status == 0) | (self.status == 1)

    @property
    def not_started(self) -> np.ndarray:
        """Particles still at the cathode: status -6 .. -1."""
        return (self.status >= -6) & (self.status <= -1)

    @property
    def lost(self) -> np.ndarray:
        """Particles lost by some mechanism: status < -6."""
        return self.status < -6

    @property
    def tracked(self) -> np.ndarray:
        """Particles currently tracked (active or passive): status >= 0."""
        return self.status >= 0

    # -- Derived scalar quantities ------------------------------------

    @property
    def n_particle(self) -> int:
        return len(self.x)

    @property
    def n_active(self) -> int:
        return int(np.sum(self.active))

    @property
    def active_charge_nC(self) -> float:
        return float(np.sum(self.charge[self.active]))

    @property
    def mean_pz_eVc(self) -> float:
        """Mean longitudinal momentum of active particles [eV/c]."""
        return float(np.mean(self.pz[self.active]))

    @property
    def reference_kinetic_energy_eV(self) -> float:
        """Reference kinetic energy [eV] from the reference momentum."""
        if self.ref_momentum_eVc <= 0:
            return 0.0
        return kinetic_energy_from_momentum(self.ref_momentum_eVc)

    @property
    def gamma_ref(self) -> float:
        if self.ref_momentum_eVc <= 0:
            return 1.0
        return gamma_from_momentum(self.ref_momentum_eVc)

    @property
    def beta_ref(self) -> float:
        return beta_from_gamma(self.gamma_ref)

    @property
    def is_valid(self) -> bool:
        n = len(self.x)
        arrays = [self.y, self.z, self.px, self.py, self.pz,
                  self.clock, self.charge, self.status]
        if self.index is not None:
            arrays.append(self.index)
        return all(len(a) == n for a in arrays)

    # -- Factory methods ----------------------------------------------

    @classmethod
    def from_arrays(
        cls,
        x, y, z, px, py, pz,
        clock, charge,
        status=None, index=None,
        **kwargs,
    ) -> "Distribution":
        """Create from raw arrays.

        Units: x/y/z [m], px/py/pz [eV/c], clock [s], charge [nC].
        """
        n = len(np.asarray(x))
        if status is None:
            status = np.full(n, 5, dtype=np.int32)  # standard particle
        return cls(
            x=np.asarray(x, dtype=float),
            y=np.asarray(y, dtype=float),
            z=np.asarray(z, dtype=float),
            px=np.asarray(px, dtype=float),
            py=np.asarray(py, dtype=float),
            pz=np.asarray(pz, dtype=float),
            clock=np.asarray(clock, dtype=float),
            charge=np.asarray(charge, dtype=float),
            status=np.asarray(status, dtype=np.int32),
            index=None if index is None else np.asarray(index, dtype=np.int32),
            **kwargs,
        )

    # -- Slicing / filtering ------------------------------------------

    def _subset(self, mask: np.ndarray, source_note: str = "") -> "Distribution":
        mask = np.asarray(mask, dtype=bool)
        return Distribution(
            x=self.x[mask].copy(), y=self.y[mask].copy(), z=self.z[mask].copy(),
            px=self.px[mask].copy(), py=self.py[mask].copy(), pz=self.pz[mask].copy(),
            clock=self.clock[mask].copy(), charge=self.charge[mask].copy(),
            status=self.status[mask].copy(),
            index=None if self.index is None else self.index[mask].copy(),
            ref_time_ns=self.ref_time_ns,
            ref_momentum_eVc=self.ref_momentum_eVc,
            total_charge_nC=float(np.sum(self.charge[mask])),
            ref_x_m=self.ref_x_m, ref_y_m=self.ref_y_m,
            ref_z_m=self.ref_z_m,
            source=self.source + source_note,
            format=self.format,
            attrs=dict(self.attrs),
        )

    def filter_active(self) -> "Distribution":
        """Return a new Distribution containing only active particles."""
        return self._subset(self.active, " (active only)")

    def sample(self, n: int, seed: int = 42) -> "Distribution":
        """Randomly sample n active particles without replacement."""
        rng = np.random.default_rng(seed)
        idx_active = np.where(self.active)[0]
        n_sample = min(n, len(idx_active))
        chosen = rng.choice(idx_active, size=n_sample, replace=False)
        keep = np.zeros(self.n_particle, dtype=bool)
        keep[chosen] = True
        return self._subset(keep, " (sampled n=" + str(n_sample) + ")")

    # -- Display ------------------------------------------------------

    def __repr__(self) -> str:
        name = Path(self.source).name if self.source else "N/A"
        return (
            "Distribution(n=" + str(self.n_particle)
            + ", active=" + str(self.n_active)
            + ", lost=" + str(int(np.sum(self.lost)))
            + ", Q=" + format(self.active_charge_nC, ".4g")
            + " nC, source='" + name + "')"
        )

    def summary(self) -> str:
        """Multi-line human-readable summary."""
        return (
            "Distribution Summary\n"
            "  Source:          " + (self.source or "N/A") + "\n"
            "  Format:          " + (self.format or "N/A") + "\n"
            "  Particles:       " + str(self.n_particle) + " total, "
            + str(self.n_active) + " active, "
            + str(int(np.sum(self.passive))) + " passive, "
            + str(int(np.sum(self.lost))) + " lost, "
            + str(int(np.sum(self.not_started))) + " not started\n"
            "  Total charge:    " + format(self.active_charge_nC, ".4f") + " nC\n"
            "  Ref. time:       " + format(self.ref_time_ns, ".4f") + " ns\n"
            "  Ref. momentum:   " + format(self.ref_momentum_eVc * 1e-6, ".4f") + " MeV/c\n"
            "  Ref. kinetic E:  " + format(self.reference_kinetic_energy_eV * 1e-6, ".4f") + " MeV\n"
            "  Ref. z:          " + format(self.ref_z_m, ".4f") + " m\n"
            "  gamma:           " + format(self.gamma_ref, ".4f") + "\n"
        )
