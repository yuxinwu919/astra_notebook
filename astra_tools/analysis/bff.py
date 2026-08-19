"""Bunch Form Factor (BFF) computation.

    F(k) = (1/Q_total) * sum_j q_j exp(i k z_j)
    BFF(k) = |F(k)|^2

The bunch form factor is the squared magnitude of the Fourier transform
of the longitudinal charge distribution; it describes the coherent
radiation spectrum of the bunch (CSR, ISR, FEL) and enters wakefield
convolutions.

近中性束团 (Σq ≈ 0): 1/Q_total 归一化发散, 回退到 |q| 归一化的结构
因子 F̃(k) = Σ|q_j| e^{ikz_j} / Σ|q_j| (F̃(0)=1) 并告警 — k≠0 处的
结构因子仍有物理意义, 全零会丢失信息 (2026-08 审计 P3-1)。

References
----------
* M. Dohlus et al., CSR in bunch compressors, DESY 12-012 (批 5: 手册 6.8 实为 WAKE namelist 章, 不定义 BFF; 原引用不实已删)
"""

from __future__ import annotations

from dataclasses import dataclass
import warnings

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

    # CSR characteristic points (detect_features=True) — 启发式定义
    # (BFF=0.9/0.5 交叉点, 2026-08 审计 P3-2: 非 Dohlus 等文献的
    # 临界波数定义; 单调递减 BFF 时 peak_k 恒等于 kmin, 仅作参考)
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
        kmin: minimum wavenumber [1/m].
        kmax: maximum wavenumber [1/m].
        nk: number of k points.
        log_spaced: log-spaced k grid (recommended for CSR).
        detect_features: compute CSR characteristic k values.
        method: 'direct' (exact O(N*nk) summation), 'fft' (binned FFT,
            ~0.03% accurate at kmax, O(N + n log n)) or 'auto'
            (fft when N*nk > 2e6).
    """
    z = np.asarray(z, dtype=float)
    charge = np.asarray(charge, dtype=float)

    k: np.ndarray
    if log_spaced:
        if kmin <= 0:
            raise ValueError("kmin must be > 0 for a log-spaced k grid")
        k = np.asarray(np.logspace(np.log10(kmin), np.log10(kmax), nk))
    else:
        k = np.asarray(np.linspace(kmin, kmax, nk))

    q_total = float(np.sum(charge))
    sum_abs_q = float(np.sum(np.abs(charge)))
    # 全零电荷 (Σ|q|==0) 或无粒子: 归一化发散 (0/0), 返回全零
    if sum_abs_q == 0 or len(z) == 0:
        return BFFResult(
            k=k, bff=np.zeros(nk), bff_amplitude=np.zeros(nk),
            wavelength=2.0 * np.pi / k,
            frequency=k * C_LIGHT / (2.0 * np.pi),
        )
    # 近中性束团 (|Σq| << Σ|q|): 1/Q_total 奇异 -> |q| 归一化结构因子
    # (2026-08 审计 P3-1: 旧实现静默归零, 丢失 k≠0 结构信息)
    use_q_abs = abs(q_total) < 1e-12 * sum_abs_q
    if use_q_abs:
        warnings.warn(
            "近中性束团 (|Σq| = %.3g nC << Σ|q| = %.3g nC): BFF 改用 "
            "|q| 归一化的结构因子 F̃(k)=Σ|q_j|e^{ikz_j}/Σ|q_j| (F̃(0)=1)"
            % (abs(q_total), sum_abs_q), UserWarning, stacklevel=2)
    w = np.abs(charge) if use_q_abs else charge
    denom = sum_abs_q if use_q_abs else q_total

    use_fft = method == "fft" or (method == "auto" and len(z) * nk > 2_000_000)
    bff: np.ndarray
    bff_amplitude: np.ndarray
    if use_fft:
        # |sum_j w_j exp(i k z_j)| on the requested k grid
        f_amp: np.ndarray = _bff_fft(z, w, k)
        bff = np.asarray((f_amp / denom) ** 2)
        bff_amplitude = np.asarray(f_amp / abs(denom))
    else:
        # Direct summation, memory-friendly for large N
        f: np.ndarray = np.zeros(nk, dtype=complex)
        for zi, qi in zip(z, w):
            f += np.asarray(qi * np.exp(1j * k * zi))
        f /= denom
        bff = np.asarray(np.abs(f) ** 2)
        bff_amplitude = np.asarray(np.abs(f))

    with np.errstate(divide="ignore"):
        wavelength: np.ndarray = np.asarray(2.0 * np.pi / k)
        wavelength[k == 0] = np.inf
    frequency: np.ndarray = np.asarray(k * C_LIGHT / (2.0 * np.pi))

    csr_critical_k = 0.0
    csr_cutoff_k = 0.0
    peak_k = 0.0

    # 特征点启发式 (P3-2): peak_k = argmax (单调 BFF 时= kmin, 无信息);
    # csr_critical_k = BFF 首次 ≤ 0.9 的 k; csr_cutoff_k = BFF 最接近
    # 0.5 的 k。均为绘图辅助, 非文献定义的 CSR 临界波数。
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

    Accuracy design (2026-08 审计 P3-3 修正误差机理描述):
      * binning 是 boxcar 卷积, 频谱乘以 sinc(k·dz/2), 相对误差
        (k·dz)²/24; 取 dz <= pi/(256·kmax) 后 <= (pi/256)²/24 ≈ 6e-6。
        此前的"每粒子相位误差 pi/512 -> |F| 地板 1e-4"描述不准确:
        实测深零点处绝对误差 ~5e-5, 主要来自采样噪声与插值离散化,
        而非相位地板。
      * zero-padding >= 16x: FFT 的 k 间距 dk 足够细, |F|² 插值误差
        ~ (dk·σz)²/8; 二次三点插值保持深零点不劣化 (实测 5.3e-5),
        但不是"精确复现"。

    内存提示 (2026-08 审计 P3-5): n_fft ≈ 16·n_grid, 长束团 (half~1 m)
    且 kmax~1e5 时 n_fft 可达 ~2^28, 复数工作集数 GB 量级; 需要更大
    范围时请用 direct 法或降低 kmax。超过 2^26 时本函数发告警。

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
    if n_fft > 2**26:
        warnings.warn(
            "BFF FFT 路径 n_fft = %d (> 2^26): 复数工作集约 %.1f GB; "
            "长束团/大 kmax 建议用 direct 法或降低 kmax"
            % (n_fft, n_fft * 16 / 1e9), UserWarning, stacklevel=2)
    f = np.fft.rfft(rho, n=n_fft)
    k_fft: np.ndarray = np.asarray(
        2.0 * np.pi * np.arange(len(f)) / (n_fft * dz_bin)
    )

    bff_grid: np.ndarray = np.asarray(np.abs(f) ** 2)
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
