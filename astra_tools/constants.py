"""Physical constants and unit conversions for beam physics.

All internal library quantities are in SI (m, s, eV/c, nC).
Constants follow CODATA 2018; file formats follow the ASTRA Manual V3.2.

References
----------
* CODATA 2018 recommended values
* ASTRA Manual V3.2, Table 1  (particle distribution file)
* ASTRA Manual V3.2, Table 4  (Xemit/Yemit/Zemit/ref/Sigma/Cemit files)
"""

from __future__ import annotations

import numpy as np

# ============================================================
# Fundamental constants (CODATA 2018)
# ============================================================

M_E_C2_EV = 0.510998950e6   # electron rest energy [eV]
C_LIGHT = 2.99792458e8      # speed of light [m/s]
E_CHARGE = 1.602176634e-19  # elementary charge [C]
M_E_KG = 9.1093837015e-31   # electron rest mass [kg]
EPS0 = 8.8541878128e-12     # vacuum permittivity [F/m] (CODATA 2018)

# ============================================================
# Relativistic kinematics
# ============================================================
# NOTE (physics audit): ASTRA stores the reference *momentum* p_ref in
# header[1] of distribution files, in eV/c — NOT the kinetic energy.
# All relativistic factors must be computed from momentum via
# gamma_from_momentum(). Treating p_ref as kinetic energy is the classic
# error that silently zeroes the normalized emittance below p = mc^2.

def _as_scalar_or_array(x):
    """Return a plain float for scalars, keep ndarray unchanged."""
    arr = np.asarray(x)
    if arr.ndim == 0:
        return float(arr)
    return arr


def kinetic_energy_from_momentum(momentum_eVc, mass_eV=M_E_C2_EV):
    """Kinetic energy [eV] from momentum [eV/c]: E_kin = sqrt(p^2c^2 + m^2c^4) - mc^2."""
    return _as_scalar_or_array(np.sqrt(np.asarray(momentum_eVc) ** 2 + mass_eV**2) - mass_eV)


def momentum_from_kinetic_energy(kinetic_energy_eV, mass_eV=M_E_C2_EV):
    """Momentum [eV/c] from kinetic energy [eV]."""
    total_eV = np.asarray(kinetic_energy_eV) + mass_eV
    return _as_scalar_or_array(np.sqrt(np.maximum(total_eV**2 - mass_eV**2, 0.0)))


def gamma_from_momentum(momentum_eVc, mass_eV=M_E_C2_EV):
    """Relativistic gamma from momentum [eV/c]: gamma = sqrt(1 + (p/mc)^2)."""
    return _as_scalar_or_array(np.sqrt(1.0 + (np.asarray(momentum_eVc) / mass_eV) ** 2))


def gamma_from_kinetic_energy(energy_eV, mass_eV=M_E_C2_EV):
    """Relativistic gamma from kinetic energy [eV]: gamma = 1 + E_kin/mc^2."""
    return _as_scalar_or_array(1.0 + np.asarray(energy_eV) / mass_eV)


def beta_from_gamma(gamma):
    """Relativistic beta from gamma (0 if gamma <= 1)."""
    g = np.asarray(gamma, dtype=float)
    return _as_scalar_or_array(np.sqrt(np.maximum(1.0 - 1.0 / g**2, 0.0)))


def beta_gamma(gamma):
    """beta*gamma from gamma: sqrt(gamma^2 - 1)."""
    g = np.asarray(gamma, dtype=float)
    return _as_scalar_or_array(np.sqrt(np.maximum(g**2 - 1.0, 0.0)))


def gamma_from_total_energy(total_eV, mass_eV=M_E_C2_EV):
    """Relativistic gamma from total energy [eV]."""
    return _as_scalar_or_array(np.asarray(total_eV) / mass_eV)

# ============================================================
# Unit conversion factors
# ============================================================

# Length
M_TO_MM = 1e3
M_TO_UM = 1e6
MM_TO_M = 1e-3
UM_TO_M = 1e-6

# Angle
RAD_TO_MRAD = 1e3
MRAD_TO_RAD = 1e-3

# Energy / momentum
EV_TO_MEV = 1e-6
EV_TO_KEV = 1e-3
MEV_TO_EV = 1e6
KEV_TO_EV = 1e3

# Emittance units.
#
# ASTRA prints emittances in 'pi mm mrad' / 'pi keV mm'. The pi marks
# the quantity as the AREA of the RMS phase-space ellipse
# (area = pi*a*b, while the RMS statistical emittance from particle data
# is eps_rms = a*b). Numerically 'pi mm mrad' and 'mm mrad' therefore
# denote the SAME value in ASTRA files: the stored number equals
# eps_rms in mm.mrad. To convert to SI multiply by 1e-6 -> m.rad.
# Do NOT multiply by pi. This matches lume-astra's parsers
# (unit 'mm-mrad', factor 1e-6) and pmd-beamphysics
# (norm_emit = sqrt(det cov(x,px)) / mc2, in meters).
EMIT_M_TO_UM = 1e6                # m.rad -> um.rad
EMIT_M_TO_MM_MRAD = 1e6           # m.rad -> mm.mrad (numerically == pi-units)
MM_MRAD_TO_M_RAD = 1e-6           # mm.mrad -> m.rad

# Charge
C_TO_NC = 1e9
NC_TO_C = 1e-9

# Time
S_TO_NS = 1e9
NS_TO_S = 1e-9
S_TO_PS = 1e12

# ============================================================
# ASTRA file-format conventions (Manual V3.2)
# ============================================================

# Distribution-file particle record: x y z px py pz clock q index status
#   x,y,z [m]; px,py,pz [eV/c]; clock [ns]; macro charge [nC];
#   particle index (1=e-, 2=e+, 3=p, 4=H+); status flag (Table 2).
ASTRA_PARTICLE_DTYPE_9 = np.dtype([
    ("x", "f8"), ("y", "f8"), ("z", "f8"),
    ("px", "f8"), ("py", "f8"), ("pz", "f8"),
    ("clock_ns", "f8"), ("charge_nC", "f8"), ("status", "f8"),
])

ASTRA_PARTICLE_DTYPE_10 = np.dtype([
    ("x", "f8"), ("y", "f8"), ("z", "f8"),
    ("px", "f8"), ("py", "f8"), ("pz", "f8"),
    ("clock_ns", "f8"), ("charge_nC", "f8"), ("index", "f8"), ("status", "f8"),
])

# Status-flag semantics (ASTRA Manual V3.2, Table 2).
#   -6 .. -1   : at the cathode, not yet started
#   < -6       : lost (aperture, backwards, discarded, ...)
#   0, 1       : passive (probe) particles — tracked but EXCLUDED from
#                emittance / statistics and from space-charge sources
#   2, 3, 4, 5 : tracking particles (probe, crossover, standard)
#   6, 8, 9+   : secondary electrons of generation 1..n — tracked
STATUS_NOT_STARTED_RANGE = (-6, -1)   # inclusive
STATUS_PASSIVE = (0, 1)
STATUS_TRACKING = (2, 3, 4, 5)

# Species mapping for the particle-index column (Table 1)
PARTICLE_SPECIES = {
    1: "electron",
    2: "positron",
    3: "proton",
    4: "hydrogen_ion",
}
