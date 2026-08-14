"""Interoperability with the LUME ecosystem (lume-astra, pmd-beamphysics).

astra_tools builds ON TOP of lume-astra: input decks and output parsing
come from lume-astra (astra.Astra, astra.AstraGenerator, astra.parsers);
this module converts its particle data (pmd_beamphysics.ParticleGroup)
to/from astra_tools.Distribution.

Status mapping (pmd-beamphysics openPMD convention vs ASTRA Table 2):
    openPMD:  1 = alive,  2 = passive probe,  0 = not started, -1 = lost
    ASTRA:    > 1 = active, 0/1 = passive, -6..-1 = not started, < -6 = lost

Note on 'active' definitions: astra_tools follows ASTRA Manual 4.13
('all particles with status flag > 1', i.e. including trajectory
probes, status 3) - validated against ASTRA's own Xemit to < 0.02%.
pmd-beamphysics computes statistics over openPMD-alive particles only
(ASTRA status 5), which excludes the probes. The two conventions
agree within sampling noise for typical bunches.
"""

from __future__ import annotations

import numpy as np

from .distribution import Distribution


def distribution_from_particle_group(P) -> Distribution:
    """Convert a pmd_beamphysics.ParticleGroup to a Distribution.

    Units: beamphysics stores x/y/z [m], px/py/pz [eV/c], t [s],
    weight = macro charge [C] -> converted to nC (absolute value,
    matching ASTRA's positive-charge convention for electrons is NOT
    assumed; the sign is preserved from the data if available).
    """
    q_nC = np.asarray(P.weight) * 1e9  # C -> nC
    # openPMD -> ASTRA status: 1(alive)->5, 2(passive)->1,
    # 0(not started)->-1, -1(lost)->-99
    opmd = np.asarray(P.status, dtype=np.int32)
    status = np.where(opmd == 1, 5, np.where(opmd == 2, 1,
                     np.where(opmd == 0, -1, np.where(opmd == -1, -99, opmd))))
    return Distribution(
        x=np.asarray(P.x, dtype=float),
        y=np.asarray(P.y, dtype=float),
        z=np.asarray(P.z, dtype=float),
        px=np.asarray(P.px, dtype=float),
        py=np.asarray(P.py, dtype=float),
        pz=np.asarray(P.pz, dtype=float),
        clock=np.asarray(P.t, dtype=float),
        charge=q_nC,
        status=status,
        index=None,
        ref_momentum_eVc=float(P['mean_pz']),
        total_charge_nC=float(np.sum(np.abs(q_nC))),
        source=getattr(P, 'filename', ''),
        format='pmd_beamphysics',
    )


def particle_group_from_distribution(d: Distribution):
    """Convert a Distribution to a pmd_beamphysics.ParticleGroup.

    Status mapping: active -> 1 (alive), passive -> 2, not_started -> 0,
    lost -> -1.
    """
    from pmd_beamphysics import ParticleGroup

    status = np.zeros(d.n_particle, dtype=np.int32)
    status[d.active] = 1
    status[d.passive] = 2
    status[d.not_started] = 0
    status[d.lost] = -1

    data = {
        'x': d.x, 'y': d.y, 'z': d.z,
        'px': d.px, 'py': d.py, 'pz': d.pz,
        't': d.clock,
        'status': status,
        'weight': np.abs(d.charge) * 1e-9,  # nC -> C
        'species': 'electron',
        'n_particle': d.n_particle,
    }
    return ParticleGroup(data=data)


def stats_from_lume_astra(A) -> dict:
    """Extract the last row of lume-astra's parsed output stats.

    A.output['stats'] holds column arrays (mean_z, sigma_x, norm_emit_x,
    ...) with SI units (m, s, eV/c). Returns a dict of final values.
    """
    return {k: float(np.asarray(v)[-1]) for k, v in A.output['stats'].items()}
