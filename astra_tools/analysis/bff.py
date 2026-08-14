"""Bunch Form Factor (BFF) computation.

    F(k) = (1/Q_total) * sum_j q_j exp(i k z_j)
    BFF(k) = |F(k)|^2

The bunch form factor is the squared magnitude of the Fourier transform
of the longitudinal charge distribution; it describes the coherent
radiation spectrum of the bunch (CSR, ISR, FEL) and enters wakefield
convolutions (ASTRA Manual V3.2, section 6.8).

References
----------
* ASTRA Manual V3.2, section 6.8 (wakefield convolution)
* M. Dohlus et al., CSR in bunch compressors, DESY 12-012
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import C_LIGHT


@dataclass
class BFFResult:
    """Bunch form factor results."""

    k: np.ndarray            # [1/m]
    bff: np.ndarray          # |F(k)|^2
    bff_amplitude: np.ndarray  # |F(k)|
    wavelength: np.ndarray   # [m]
    frequency: np.ndarray    # [Hz]

    # CSR characteristic points (detect_features=True)
    csr_critical_k: float = 0.0
    csr_cutoff_k: float = 0.0
    peak_k: float = 0.0


def compute_bff(
    z: np.ndarray,
    charge: np.ndarray,
    kmin: float = 1.0,
    kmax: float = 1e5,
    nk: int = 200,
    log_spaced: bool = True,
    detect_features: bool = False,
) -> BFFResult:
    """Compute the bunch form factor from longitudinal particle data.

    Args:
        z: longitudinal positions [m] (active particles).
        charge: macro-particle charges [nC].
        kmin, kmax: wavenumber range [1/m].
        nk: number of k points.
        log_spaced: log-spaced k grid (recommended for CSR).
        detect_features: compute CSR characteristic k values.
    """
    z = np.asarray(z, dtype=float)
    charge = np.asarray(charge, dtype=float)

    if log_spaced:
        k = np.logspace(np.log10(kmin), np.log10(kmax), nk)
    else:
        k = np.linspace(kmin, kmax, nk)

    q_total = float(np.sum(charge))
    if q_total == 0 or len(z) == 0:
        return BFFResult(
            k=k, bff=np.zeros(nk), bff_amplitude=np.zeros(nk),
            wavelength=2.0 * np.pi / k,
            frequency=k * C_LIGHT / (2.0 * np.pi),
        )

    # Direct summation, memory-friendly for large N
    f = np.zeros(nk, dtype=complex)
    for zi, qi in zip(z, charge):
        f += qi * np.exp(1j * k * zi)
    f /= q_total

    bff = np.abs(f) ** 2
    bff_amplitude = np.abs(f)

    with np.errstate(divide="ignore"):
        wavelength = 2.0 * np.pi / k
        wavelength[k == 0] = np.inf
    frequency = k * C_LIGHT / (2.0 * np.pi)

    csr_critical_k = 0.0
    csr_cutoff_k = 0.0
    peak_k = 0.0

    if detect_features and len(k) > 2:
        peak_k = float(k[int(np.argmax(bff))])
        below = np.where(bff <= 0.9)[0]
        csr_critical_k = float(k[below[0]]) if len(below) else float(kmax)
        csr_cutoff_k = float(k[int(np.argmin(np.abs(bff - 0.5)))])

    return BFFResult(
        k=k, bff=bff, bff_amplitude=bff_amplitude,
        wavelength=wavelength, frequency=frequency,
        csr_critical_k=csr_critical_k,
        csr_cutoff_k=csr_cutoff_k,
        peak_k=peak_k,
    )
