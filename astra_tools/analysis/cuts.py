"""Phase space cuts and rotations (postpro 5.6.4).

Apply transverse/longitudinal/energy windows or a radial aperture to a
Distribution, or rotate the x-y plane; returns a new Distribution ready
for further tracking or export.
"""

from __future__ import annotations

import numpy as np

from ..constants import kinetic_energy_from_momentum
from ..distribution import Distribution


def cut_distribution(
    dist: Distribution,
    x_range=None,
    y_range=None,
    z_range=None,
    e_range=None,
    r_aperture=None,
):
    """Cut a distribution by windows (active particles only; cut
    particles are relabelled status -31, 'discarded by user', per the
    ASTRA Manual Table 2).

    Args:
        dist: input distribution.
        x_range / y_range: (min, max) [m] transverse windows.
        z_range: (min, max) [m] longitudinal window (absolute z).
        e_range: (min, max) [eV] kinetic-energy window.
        r_aperture: radial aperture radius [m], sqrt(x^2+y^2) < r.

    Returns:
        (dist_cut, mask) - mask marks the removed particles.
    """
    cut = np.zeros(dist.n_particle, dtype=bool)
    if x_range is not None:
        cut |= (dist.x < x_range[0]) | (dist.x > x_range[1])
    if y_range is not None:
        cut |= (dist.y < y_range[0]) | (dist.y > y_range[1])
    if z_range is not None:
        cut |= (dist.z < z_range[0]) | (dist.z > z_range[1])
    if e_range is not None:
        e = kinetic_energy_from_momentum(dist.pz)
        cut |= (e < e_range[0]) | (e > e_range[1])
    if r_aperture is not None:
        r = np.sqrt(dist.x**2 + dist.y**2)
        cut |= r > r_aperture

    cut &= dist.active
    status = dist.status.copy()
    status[cut] = -31
    out = Distribution(
        x=dist.x.copy(), y=dist.y.copy(), z=dist.z.copy(),
        px=dist.px.copy(), py=dist.py.copy(), pz=dist.pz.copy(),
        clock=dist.clock.copy(), charge=dist.charge.copy(),
        status=status,
        index=None if dist.index is None else dist.index.copy(),
        ref_time_ns=dist.ref_time_ns, ref_momentum_eVc=dist.ref_momentum_eVc,
        total_charge_nC=dist.total_charge_nC,
        ref_x_m=dist.ref_x_m, ref_y_m=dist.ref_y_m, ref_z_m=dist.ref_z_m,
        source=dist.source + " (cut)", format=dist.format, attrs=dict(dist.attrs),
    )
    return out, cut


def rotate_phase_space(dist: Distribution, angle_deg: float):
    """Rotate the x-y plane (coordinates and momenta together)."""
    th = np.radians(angle_deg)
    c, s = np.cos(th), np.sin(th)
    out = Distribution(
        x=c * dist.x + s * dist.y,
        y=-s * dist.x + c * dist.y,
        z=dist.z.copy(),
        px=c * dist.px + s * dist.py,
        py=-s * dist.px + c * dist.py,
        pz=dist.pz.copy(), clock=dist.clock.copy(), charge=dist.charge.copy(),
        status=dist.status.copy(),
        index=None if dist.index is None else dist.index.copy(),
        ref_time_ns=dist.ref_time_ns, ref_momentum_eVc=dist.ref_momentum_eVc,
        total_charge_nC=dist.total_charge_nC,
        ref_x_m=dist.ref_x_m, ref_y_m=dist.ref_y_m, ref_z_m=dist.ref_z_m,
        source=dist.source + " (rotated %.1f deg)" % angle_deg,
        format=dist.format, attrs=dict(dist.attrs),
    )
    return out
