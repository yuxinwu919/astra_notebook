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
    method: str = "auto",
) -> BFFResult:
    """Compute the bunch form factor from longitudinal particle data.

    Args:
        z: longitudinal positions [m] (active particles).
        charge: macro-particle charges [nC].
        kmin, kmax: wavenumber range [1/m].
        nk: number of k points.
        log_spaced: log-spaced k grid (recommended for CSR).
        detect_features: compute CSR characteristic k values.
        method: 'direct' (exact O(N*nk) summation), 'fft' (binned FFT,
            ~0.03% accurate at kmax, O(N + n log n)) or 'auto'
            (fft when N*nk > 2e6).
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

    use_fft = method == "fft" or (method == "auto" and len(z) * nk > 2_000_000)
    if use_fft:
        # |sum_j q_j exp(i k z_j)| on the requested k grid
        f_amp = _bff_fft(z, charge, k)
        bff = (f_amp / q_total) ** 2
        bff_amplitude = f_amp / abs(q_total)
    else:
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


def _bff_fft(z: np.ndarray, charge: np.ndarray, k: np.ndarray) -> np.ndarray:
    """FFT-based |F(k)| = |sum_j q_j exp(i k z_j)| on the requested k grid.

    Only the amplitude is needed downstream, so the phase is never
    interpolated (a rotating phase between FFT samples would make a
    linear chord cut inside the unit circle and underestimate |F|).
    The smooth, non-negative quantity |F|^2 is interpolated instead.

    Accuracy design:
      * bin width dz <= pi/(256 kmax): per-particle phase error inside
        a bin <= pi/512 ~ 0.006 rad (absolute |F| floor ~1e-4)
      * zero-padding >= 16x: the FFT k-spacing dk is fine enough that
        the |F|^2 interpolation error ~ (dk sigma_z)^2/8 stays small;
        quadratic 3-point interpolation reproduces deep nulls exactly

    O(N + n_fft log n_fft) instead of O(N * nk).
    """
    zc = float(np.mean(z))
    zc_shifted = z - zc
    half = float(np.max(np.abs(zc_shifted)))
    half = half * 1.1 + 1e-6          # 10% guard band

    kmax = float(np.max(k))
    dz0 = np.pi / (256.0 * kmax)
    n_grid = max(int(np.ceil(2.0 * half / dz0)) + 1, 16)
    dz_bin = 2.0 * half / n_grid      # <= dz0
    rho, _ = np.histogram(zc_shifted, bins=n_grid,
                          range=(-half, half), weights=charge)

    n_fft = 1 << (16 * n_grid - 1).bit_length()   # >= 16x zero padding
    f = np.fft.rfft(rho, n=n_fft)
    k_fft = 2.0 * np.pi * np.arange(len(f)) / (n_fft * dz_bin)

    bff_grid = np.abs(f) ** 2
    # quadratic (3-point) interpolation of |F|^2: reproduces the sharp
    # minima (deep nulls) exactly, unlike a linear chord
    j = np.clip(np.searchsorted(k_fft, k, side="right") - 1,
                1, len(k_fft) - 2)
    k1 = k_fft[j]
    dk_grid = k_fft[j + 1] - k_fft[j - 1]
    t = (k - k1) / (dk_grid * 0.5)
    y0 = bff_grid[j - 1]
    y1 = bff_grid[j]
    y2 = bff_grid[j + 1]
    bff_interp = y1 + 0.5 * t * (y2 - y0) + 0.5 * t**2 * (y2 - 2.0 * y1 + y0)
    return np.sqrt(np.maximum(bff_interp, 0.0))
