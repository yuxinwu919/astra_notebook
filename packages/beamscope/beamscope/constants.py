"""Physical constants and unit conversions for particle beam physics.

All values in SI base units unless otherwise noted.

References
----------
* CODATA 2018 recommended values
* ASTRA Manual V3.2, Table 1 (particle distribution format)
* ASTRA Manual V3.2, Table 4 (emit file columns)
"""

import numpy as np

# ============================================================
# Fundamental physical constants (CODATA 2018)
# ============================================================

# Electron rest mass energy [eV]
M_E_C2_EV = 0.510998950e6  # m_e * c^2 / e

# Speed of light [m/s]
C_LIGHT = 2.99792458e8

# Electron charge [C]
E_CHARGE = 1.602176634e-19

# Electron rest mass [kg]
M_E_KG = 9.1093837015e-31

# ============================================================
# Relativistic kinematics
# ============================================================

def kinetic_energy_from_momentum(
    momentum_eVc: float,
    mass_eV: float = M_E_C2_EV,
) -> float:
    """Compute kinetic energy [eV] from momentum [eV/c].

    E_kin = sqrt(p²c² + m²c⁴) - mc²

    This is the CORRECT conversion for ASTRA data, where the binary header
    stores reference *momentum* p_ref [eV/c], NOT kinetic energy.

    Args:
        momentum_eVc: Particle momentum in eV/c.
        mass_eV: Rest mass energy in eV (default: electron).

    Returns:
        Kinetic energy in eV.
    """
    total_eV = np.sqrt(momentum_eVc**2 + mass_eV**2)
    return total_eV - mass_eV


def momentum_from_kinetic_energy(
    kinetic_energy_eV: float,
    mass_eV: float = M_E_C2_EV,
) -> float:
    """Compute momentum [eV/c] from kinetic energy [eV].

    p = sqrt((E_kin + mc²)² - m²c⁴) / c

    Args:
        kinetic_energy_eV: Kinetic energy in eV.
        mass_eV: Rest mass energy in eV (default: electron).

    Returns:
        Momentum in eV/c.
    """
    total_eV = kinetic_energy_eV + mass_eV
    return np.sqrt(total_eV**2 - mass_eV**2)


def gamma_from_momentum(
    momentum_eVc: float,
    mass_eV: float = M_E_C2_EV,
) -> float:
    """Compute relativistic gamma from momentum [eV/c].

    γ = sqrt(p² + m²) / m = sqrt(1 + (p/m)²)

    This is the preferred function for ASTRA data since the reference
    quantity stored in binary files is momentum, not kinetic energy.

    Args:
        momentum_eVc: Particle momentum in eV/c.
        mass_eV: Rest mass energy in eV (default: electron).

    Returns:
        Relativistic gamma factor.
    """
    return np.sqrt(1.0 + (momentum_eVc / mass_eV) ** 2)


def gamma_from_energy(energy_eV: float, mass_eV: float = M_E_C2_EV) -> float:
    """Compute relativistic gamma from kinetic energy [eV].

    γ = 1 + E_kin / (m·c²)

    NOTE: Prefer gamma_from_momentum() for ASTRA data since ASTRA stores
    reference *momentum* (eV/c), not kinetic energy.

    Args:
        energy_eV: Kinetic energy in eV.
        mass_eV: Rest mass energy in eV (default: electron).

    Returns:
        Relativistic gamma factor.
    """
    return 1.0 + energy_eV / mass_eV


def beta_from_gamma(gamma: float) -> float:
    """Compute relativistic beta from gamma."""
    if gamma <= 1.0:
        return 0.0
    return np.sqrt(1.0 - 1.0 / gamma**2)


def beta_gamma(gamma: float) -> float:
    """Compute beta*gamma product from gamma.

    βγ = sqrt(γ² - 1)
    """
    if gamma <= 1.0:
        return 0.0
    return np.sqrt(gamma**2 - 1.0)


def momentum_from_energy(energy_eV: float, mass_eV: float = M_E_C2_EV) -> float:
    """Compute momentum [eV/c] from kinetic energy [eV].

    Deprecated: use momentum_from_kinetic_energy() for clarity.
    """
    return momentum_from_kinetic_energy(energy_eV, mass_eV)


# ============================================================
# Unit conversion factors
# ============================================================

# Length
M_TO_MM = 1e3
M_TO_UM = 1e6
MM_TO_M = 1e-3

# Angle
RAD_TO_MRAD = 1e3

# Energy
EV_TO_MEV = 1e-6
EV_TO_KEV = 1e-3
MEV_TO_EV = 1e6

# Emittance: m·rad → μm·rad, π·m·rad → π·mm·mrad
EMIT_M_TO_UM = 1e6
EMIT_M_TO_MM_MRAD = 1e6

# Charge
C_TO_NC = 1e9
NC_TO_C = 1e-9

# Time
S_TO_NS = 1e9
NS_TO_S = 1e-9
S_TO_PS = 1e12
