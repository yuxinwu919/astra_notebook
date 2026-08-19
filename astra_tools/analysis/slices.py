"""Longitudinal slice analysis.

Computes slice-by-slice beam parameters (current, charge, emittance,
energy chirp) for studying slice emittance and current profiles.

Binning strategies follow the ASTRA Manual V3.2, section 6.8:
  * 'equi_spaced': equal z-width bins
  * 'equi_charge' : equal charge per bin
  * 'equi_energy': equal charge per bin sorted by kinetic energy
    (postpro 5.6.3 items 8/9: slices w.r.t. the bunch energy)

Physics notes (audited):
  * slices are taken from active particles (status > 1) only
  * slice divergence uses the canonical momentum (manual 4.13.1)
    normalized by the global reference momentum, u' = p~u / p_ref, and
    the normalized emittance is beta*gamma(p_ref) * eps_geom - the
    same convention as the Xemit/Yemit files
  * current: I = Q/dt with dt = dz/(beta*c); the per-slice beta/gamma
    used for the current come from the slice mean momentum
  * per-slice energy: E_kin = sqrt((px^2+py^2+pz^2) + m^2c^4) - mc^2
    (FULL momentum per particle, 2026-08 audit P2-2)
  * equi_charge binning uses |q| (mixed-sign safe); slice charge stays
    signed internally (sign is conventional, |Q| at display)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from ..constants import C_LIGHT, M_E_C2_EV, NC_TO_C, \
    kinetic_energy_from_momentum_vector, gamma_from_momentum
from ..distribution import Distribution
from .emittance import compute_geometric_emittance, canonical_signs


@dataclass
class SliceAnalysis:
    """Results of longitudinal slice analysis (SI units)."""

    n_slices: int
    binning_mode: str

    z_centers: np.ndarray     # [m]; equi_energy 分箱时为动能 [eV]
    z_edges: np.ndarray       # [m] (n_slices+1); equi_energy 时为 [eV]
    current: np.ndarray       # [A]
    charge: np.ndarray        # [nC]
    n_particles: np.ndarray   # int

    mean_x: np.ndarray        # [m]
    mean_y: np.ndarray        # [m]
    mean_pz: np.ndarray       # [eV/c]
    mean_kinetic_energy_eV: np.ndarray  # [eV]

    sig_x: np.ndarray         # [m]
    sig_y: np.ndarray         # [m]
    sig_E_over_E: np.ndarray  # relative energy spread

    mean_xp: np.ndarray       # [rad] slice 内 x' 平均 (发散角, 5.6.3 项 5)
    mean_yp: np.ndarray       # [rad]
    sig_xp: np.ndarray        # [rad] slice 内 x' rms
    sig_yp: np.ndarray        # [rad]

    emit_x_norm: np.ndarray   # [m.rad]
    emit_y_norm: np.ndarray   # [m.rad]

    gamma_per_slice: np.ndarray
    beta_per_slice: np.ndarray


def compute_slice_analysis(
    dist: Distribution,
    n_slices: int = 20,
    binning: str = "equi_spaced",
    ref_momentum_eVc: Optional[float] = None,
    bz_on_axis_T: float = 0.0,
) -> SliceAnalysis:
    """Compute slice-by-slice beam parameters.

    Args:
        dist: particle distribution (active particles used).
        n_slices: number of longitudinal slices.
        binning: 'equi_spaced' (等 z 宽), 'equi_charge' (按 z 等电荷),
            'equi_energy' (按动能等电荷, 手册 5.6.3 项 8/9)。
        ref_momentum_eVc: global reference momentum [eV/c] used for the
            slice divergences and normalized emittances (Xemit convention).
        bz_on_axis_T: on-axis solenoid field at the bunch center [T] for
            the canonical momentum (manual 4.13.1).
    """
    if binning not in ("equi_spaced", "equi_charge", "equi_energy"):
        raise ValueError(
            "binning must be 'equi_spaced'/'equi_charge'/'equi_energy'")

    d = dist.filter_active()
    z = d.z
    n = d.n_active
    q = d.charge

    if n < 6:
        raise ValueError("too few active particles for slice analysis")
    if n < n_slices:
        n_slices = max(n // 3, 1)

    if ref_momentum_eVc is None:
        ref_momentum_eVc = d.ref_momentum_or_mean()
    if ref_momentum_eVc <= 0:
        raise ValueError(
            "reference momentum is zero/negative (beam at rest?); "
            "slice divergences and normalized emittances are undefined")

    # -- Binning --
    # 通用结构: bin_vals = 分箱变量, edges = 边界 (与 bin_vals 同单位)。
    # equi_spaced/equi_charge 按 z 分箱; equi_energy 按动能 E 分箱
    # (2026-08 审计 P1: 旧实现把 equi_energy 静默退化为 z 等宽)。
    if binning == "equi_charge":
        sort_idx = np.argsort(z)
        z_sorted = z[sort_idx]
        q_sorted = q[sort_idx]
        q_cumsum = np.cumsum(np.abs(q_sorted))  # |q|: sign cannot break bins
        q_total = float(q_cumsum[-1])
        if q_total == 0:
            # degenerate: fall back to equi-spaced
            edges = np.linspace(float(np.min(z)), float(np.max(z)), n_slices + 1)
        else:
            q_per_slice = q_total / n_slices
            edges = np.empty(n_slices + 1)
            edges[0] = z_sorted[0]
            for i in range(1, n_slices):
                idx = min(np.searchsorted(q_cumsum, i * q_per_slice), n - 1)
                edges[i] = float(z_sorted[idx])
            edges[-1] = z_sorted[-1]
            # 重复 z 值会使边界塌缩 (尤其最后一条边); 用"平均箱宽
            # 的 1%" 修复全部 n_slices+1 条边并保持严格递增。delta 函数
            # 式电荷的真实峰值电流为无穷大, 这里用箱宽正则化给出有界
            # 数值 (避免 1e-12 m 级别的假 dz 造出 ~1e12 A 的假电流)
            span = float(z_sorted[-1] - z_sorted[0])
            eps = max(span / (100.0 * n_slices), 1e-12)
            for i in range(1, n_slices + 1):
                if edges[i] <= edges[i - 1]:
                    edges[i] = edges[i - 1] + eps
        bin_vals = z
    elif binning == "equi_energy":
        e_all = kinetic_energy_from_momentum_vector(d.px, d.py, d.pz)
        sort_idx = np.argsort(e_all)
        e_sorted = e_all[sort_idx]
        q_sorted = q[sort_idx]
        q_cumsum = np.cumsum(np.abs(q_sorted))
        q_total = float(q_cumsum[-1])
        if q_total == 0:
            edges = np.linspace(float(e_sorted[0]), float(e_sorted[-1]),
                                n_slices + 1)
        else:
            q_per_slice = q_total / n_slices
            edges = np.empty(n_slices + 1)
            edges[0] = e_sorted[0]
            for i in range(1, n_slices):
                idx = min(np.searchsorted(q_cumsum, i * q_per_slice), n - 1)
                edges[i] = float(e_sorted[idx])
            edges[-1] = e_sorted[-1]
            # 同上: 修复塌缩边界
            span = float(e_sorted[-1] - e_sorted[0])
            eps = max(span / (100.0 * n_slices), 1e-6)
            for i in range(1, n_slices + 1):
                if edges[i] <= edges[i - 1]:
                    edges[i] = edges[i - 1] + eps
        bin_vals = e_all
    else:
        edges = np.linspace(float(np.min(z)), float(np.max(z)), n_slices + 1)
        bin_vals = z

    z_centers = 0.5 * (edges[:-1] + edges[1:])   # equi_energy 时为能量中心
    dz = np.diff(edges)

    # -- Pre-allocate --
    n_part = np.zeros(n_slices, dtype=int)
    charge_s = np.zeros(n_slices)
    mean_x = np.zeros(n_slices)
    mean_y = np.zeros(n_slices)
    sig_x = np.zeros(n_slices)
    sig_y = np.zeros(n_slices)
    mean_xp = np.zeros(n_slices)
    mean_yp = np.zeros(n_slices)
    sig_xp = np.zeros(n_slices)
    sig_yp = np.zeros(n_slices)
    mean_pz = np.zeros(n_slices)
    mean_e = np.zeros(n_slices)
    sig_ee = np.zeros(n_slices)
    emit_xn = np.zeros(n_slices)
    emit_yn = np.zeros(n_slices)
    gamma_s = np.ones(n_slices)
    beta_s = np.zeros(n_slices)
    zspan = np.zeros(n_slices)     # equi_energy 电流用 (真实 z 跨度)

    s_can_all = canonical_signs(d)
    for i in range(n_slices):
        if i == n_slices - 1:
            mask = (bin_vals >= edges[i]) & (bin_vals <= edges[i + 1])
        else:
            mask = (bin_vals >= edges[i]) & (bin_vals < edges[i + 1])
        idx = np.where(mask)[0]
        n_part[i] = len(idx)
        if n_part[i] < 3:
            continue

        xi, yi = d.x[idx], d.y[idx]
        pxi, pyi, pzi = d.px[idx], d.py[idx], d.pz[idx]
        qi = q[idx]
        s_i = s_can_all[idx]
        zspan[i] = float(np.max(d.z[idx]) - np.min(d.z[idx]))

        charge_s[i] = float(np.sum(qi))   # signed (internal); display uses |Q|
        mean_x[i] = float(np.mean(xi))
        mean_y[i] = float(np.mean(yi))
        sig_x[i] = float(np.std(xi - mean_x[i]))   # ddof=0, matches ASTRA
        sig_y[i] = float(np.std(yi - mean_y[i]))
        mean_pz[i] = float(np.mean(pzi))
        # 全动量动能 (2026-08 审计 P2-2: pz-only 对大发散束低估 ~60%)
        e_i = kinetic_energy_from_momentum_vector(pxi, pyi, pzi)
        mean_e[i] = float(np.mean(e_i))
        sig_ee[i] = float(np.std(e_i) / mean_e[i]) if mean_e[i] else 0.0

        # Canonical momentum (manual 4.13.1), centered per slice, divided
        # by the global reference momentum - the Xemit convention.
        # 种类感知符号 (2026-08 审计 F4): 正电荷种类 s=-1。
        ptx = pxi + s_i * 0.5 * C_LIGHT * bz_on_axis_T * yi
        pty = pyi - s_i * 0.5 * C_LIGHT * bz_on_axis_T * xi
        xp = (ptx - np.mean(ptx)) / ref_momentum_eVc
        yp = (pty - np.mean(pty)) / ref_momentum_eVc
        # slice 发散角统计 (手册 5.6.3 项 5: px/pz rms 与 avr vs z)
        mean_xp[i] = float(np.mean(xp))
        mean_yp[i] = float(np.mean(yp))
        sig_xp[i] = float(np.std(xp))
        sig_yp[i] = float(np.std(yp))
        # 群体矩 (无加权), 与 compute_statistics 默认及 ASTRA 一致;
        # |q| 仅用于分箱
        ex = compute_geometric_emittance(xi - mean_x[i], xp)
        ey = compute_geometric_emittance(yi - mean_y[i], yp)

        # Normalized emittance w.r.t. the reference momentum (Xemit
        # convention); the per-slice gamma below is only for the current.
        g_ref = gamma_from_momentum(ref_momentum_eVc)
        bg = float(np.sqrt(max(g_ref**2 - 1.0, 0.0)))
        emit_xn[i] = bg * ex
        emit_yn[i] = bg * ey

        # slice 纵向速度: v_z/c = pz / sqrt(|p|^2 + m^2c^4) 的平均
        # (2026-08: 大发散束下不能用 beta(gamma(mean pz)) 近似)
        p2_i = pxi**2 + pyi**2 + pzi**2
        e_tot_i = np.sqrt(p2_i + M_E_C2_EV**2)
        beta_s[i] = float(np.mean(pzi / e_tot_i))
        p_rms = float(np.sqrt(np.mean(p2_i)))
        gamma_s[i] = gamma_from_momentum(p_rms)

    # -- Current: I = Q/dt, dt = dz/(beta c) --
    # equi_energy 分箱时边界单位是 eV, 电流用箱内粒子的真实 z 跨度
    # (2026-08: dz 为能量宽度时不可直接当长度用)。
    dz_use = zspan if binning == "equi_energy" else dz
    current = np.zeros(n_slices)
    for i in range(n_slices):
        v = beta_s[i] * C_LIGHT if beta_s[i] > 0 else C_LIGHT
        dt = dz_use[i] / v if v > 0 and dz_use[i] > 0 else 1.0
        current[i] = charge_s[i] * NC_TO_C / dt

    return SliceAnalysis(
        n_slices=n_slices, binning_mode=binning,
        z_centers=z_centers, z_edges=edges,
        current=current, charge=charge_s, n_particles=n_part,
        mean_x=mean_x, mean_y=mean_y,
        mean_pz=mean_pz, mean_kinetic_energy_eV=mean_e,
        sig_x=sig_x, sig_y=sig_y, sig_E_over_E=sig_ee,
        mean_xp=mean_xp, mean_yp=mean_yp,
        sig_xp=sig_xp, sig_yp=sig_yp,
        emit_x_norm=emit_xn, emit_y_norm=emit_yn,
        gamma_per_slice=gamma_s, beta_per_slice=beta_s,
    )
