"""Bunch Form Factor (BFF) computation.

The bunch form factor is the squared magnitude of the Fourier transform
of the longitudinal charge distribution, used for coherent radiation
calculations (CSR, ISR, FEL).

  F(k) = (1/Q_total) · Σ q_j · exp(i · k · z_j)
  BFF(k) = |F(k)|²

References:
  - ASTRA Manual V3.2, §6.8 (wakefield convolution)
  - M. Dohlus et al., "CSR in bunch compressors", DESY 12-012
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class BFFResult:
    """Bunch form factor results."""

    k: np.ndarray          # [1/m] wavenumber array
    bff: np.ndarray        # |F(k)|² — bunch form factor
    bff_amplitude: np.ndarray  # |F(k)| — amplitude

    # Derived quantities
    wavelength: np.ndarray  # [m] λ = 2π/k
    frequency: np.ndarray   # [Hz] f = k·c/(2π)

    # CSR characteristic points (computed if detect_features=True)
    csr_critical_k: float = 0.0      # k where BFF starts to drop
    csr_cutoff_k: float = 0.0        # k where BFF = 0.5
    peak_k: float = 0.0              # k of maximum BFF


def compute_bff(
    z: np.ndarray,
    charge: np.ndarray,
    kmin: float = 1.0,
    kmax: float = 1e5,
    nk: int = 200,
    log_spaced: bool = True,
    detect_features: bool = False,
    use_vectorized: bool = True,
) -> BFFResult:
    """Compute the bunch form factor from longitudinal particle data.

    F(k) = (1/Q) Σ q_j · exp(i · k · z_j)
    BFF(k) = |F(k)|²

    For large particle counts, the vectorized method uses loop-based
    computation with broadcasting which is more memory-efficient than
    a full (nk × N) complex matrix.

    Args:
        z: Longitudinal positions [m] (active particles only).
        charge: Macro-particle charges [nC].
        kmin: Minimum wavenumber [1/m].
        kmax: Maximum wavenumber [1/m].
        nk: Number of k points.
        log_spaced: If True, use log-spaced k array (recommended for CSR).
        detect_features: If True, compute CSR characteristic k values.
        use_vectorized: If True, use loop-based vectorized computation
            (memory-friendly). If False, use full broadcasting (faster but
            uses O(nk·N) memory).

    Returns:
        BFFResult with k, bff, amplitude, and derived quantities.
    """
    if log_spaced:
        k = np.logspace(np.log10(kmin), np.log10(kmax), nk)
    else:
        k = np.linspace(kmin, kmax, nk)

    q_total = np.sum(charge)
    if q_total == 0:
        return BFFResult(
            k=k, bff=np.zeros(nk), bff_amplitude=np.zeros(nk),
            wavelength=2.0 * np.pi / k,
            frequency=k * 2.99792458e8 / (2.0 * np.pi),
        )

    if use_vectorized:
        # Loop over particles — memory-friendly for large N
        F = np.zeros(nk, dtype=complex)
        for i in range(len(z)):
            F += charge[i] * np.exp(1j * k * z[i])
        F /= q_total
    else:
        # Full broadcasting — fast but O(nk·N) memory
        F = np.sum(
            charge[np.newaxis, :] * np.exp(1j * np.outer(k, z)),
            axis=1,
        ) / q_total

    bff = np.abs(F) ** 2
    bff_amplitude = np.abs(F)

    # Derived quantities
    with np.errstate(divide="ignore"):
        wavelength = 2.0 * np.pi / k
        wavelength[k == 0] = np.inf
    frequency = k * 2.99792458e8 / (2.0 * np.pi)

    # CSR feature detection
    csr_critical_k = 0.0
    csr_cutoff_k = 0.0
    peak_k = 0.0

    if detect_features and len(k) > 2:
        # Peak k: where BFF is maximum
        peak_idx = int(np.argmax(bff))
        peak_k = float(k[peak_idx])

        # Critical k: where BFF first drops below 0.9 (from low-k side)
        above_threshold = bff > 0.9
        falloff_indices = np.where(~above_threshold)[0]
        if len(falloff_indices) > 0:
            csr_critical_k = float(k[falloff_indices[0]])
        else:
            csr_critical_k = kmax

        # Cutoff k: where BFF = 0.5
        try:
            csr_cutoff_k = float(k[np.argmin(np.abs(bff - 0.5))])
        except (ValueError, IndexError):
            csr_cutoff_k = 0.0

    return BFFResult(
        k=k,
        bff=bff,
        bff_amplitude=bff_amplitude,
        wavelength=wavelength,
        frequency=frequency,
        csr_critical_k=csr_critical_k,
        csr_cutoff_k=csr_cutoff_k,
        peak_k=peak_k,
    )
