"""高级绘图 (lineplot 菜单 2/3/4 与 postpro 扩展).

对照 ASTRA 手册 5.5.2-5.5.4 / 5.6.3:
  * 损失与能量沉积 (LandF), 束载
  * beta/alpha 函数、相位推进、相干长度
  * 相位扫描 (PScan), 参数/误差扫描 (Scan/Error)
  * 缩减发射度 (Xemit2), trace-space (TRemit), 核心发射度 (Cemit),
    拉莫尔角 (Larmor), 粒子密度 (Density)
  * 探针轨迹与空间电荷场 (track), 阴极发射过程 (Cathode)
  * slice 失配参数与核心亮度 (postpro)

发射度显示统一 [pi mm mrad] (数值与 ASTRA 打印一致)。
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ..constants import C_LIGHT, M_E_C2_EV
from ..io.astra_misc import read_track_file, read_cathode_file

HBAR_EVS = 6.582119569e-16   # hbar [eV.s]
M0_C2_EV = 0.510998950e6


def _ax(ax, figsize):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    return fig, ax


def plot_losses(landf, ax=None, figsize=(8, 4), title=None):
    """粒子损失数/m 与能量沉积 (LandF 文件, lineplot 菜单 2 项 5/6)."""
    fig, ax = _ax(ax, figsize)
    ax.plot(landf["landf_z"], landf["landf_n_lost"], label="lost particles / m")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("lost particles per m")
    ax2 = ax.twinx()
    ax2.plot(landf["landf_z"], landf["landf_energy_deposited"],
             color="C1", label="deposited energy")
    ax2.set_ylabel("deposited energy [J/m]")
    ax.set_title(title or "particle loss and energy deposition")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], fontsize=9)
    fig.tight_layout()
    return fig


def plot_beam_loading(landf, ax=None, figsize=(8, 4), title=None):
    """束载: 束团与场交换的能量 (LandF 列 6, 菜单 2 项 7)."""
    fig, ax = _ax(ax, figsize)
    ax.plot(landf["landf_z"], landf["landf_energy_exchange"], label="energy exchange")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("energy exchange [J/m]")
    ax.set_title(title or "beam loading")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_beta_alpha(emit, ref=None, ax=None, figsize=(8, 5), title=None):
    """光学 beta/alpha 函数 (菜单 2 项 12/13):
    beta = sigma_u^2 / eps_geom,  alpha = -cov(u,u')/eps_geom。

    手册 5.5.2 项 12/13 用几何 RMS 发射度 (项 15 相干长度才用
    归一化 ε_n): eps_geom = eps_n / (beta*gamma), 其中 beta*gamma
    由参考粒子动量轨迹 ref.pz [eV/c] 给出 (批 2 修复: 此前误用
    ε_n 使 β/α 小 βγ≈1958 倍)。
    """
    if ref is None:
        raise ValueError(
            "plot_beta_alpha 需要参考粒子轨迹 (ref) 以计算 beta*gamma")
    fig, ax = _ax(ax, figsize)
    ax2 = ax.twinx()   # 只建一次右轴, 避免两个 twinx 重叠
    # beta*gamma = p/mc (由 gamma=sqrt(1+(p/mc)^2) 严格推出; 复核纠正)
    bg = np.interp(emit.x.z, ref.z,
                   np.maximum(ref.pz / M_E_C2_EV, 0.0))
    for e, lbl in ((emit.x, "x"), (emit.y, "y")):
        eps_geom = np.maximum(e.emit / bg, 1e-30)
        beta = e.rms**2 / eps_geom
        alpha = -e.corr / np.maximum(e.rms, 1e-30) * beta
        ax.plot(e.z, beta, label="$\\beta_%s$ [m]" % lbl)
        ax2.plot(e.z, alpha, ls="--", label="$\\alpha_%s$" % lbl)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("beta function [m]")
    ax2.set_ylabel("alpha")
    ax.set_title(title or "optical functions (from Xemit)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=9)
    fig.tight_layout()
    return fig


def plot_phase_advance(emit, ref=None, ax=None, figsize=(8, 4), title=None):
    """相位推进 theta = int(1/beta) dz (菜单 2 项 14).

    与 plot_beta_alpha 同口径: beta 由几何发射度给出 (批 2 修复)。
    """
    if ref is None:
        raise ValueError(
            "plot_phase_advance 需要参考粒子轨迹 (ref) 以计算 beta*gamma")
    fig, ax = _ax(ax, figsize)
    # beta*gamma = p/mc (由 gamma=sqrt(1+(p/mc)^2) 严格推出; 复核纠正)
    bg = np.interp(emit.x.z, ref.z,
                   np.maximum(ref.pz / M_E_C2_EV, 0.0))
    for e, lbl in ((emit.x, "x"), (emit.y, "y")):
        eps_geom = np.maximum(e.emit / bg, 1e-30)
        beta = e.rms**2 / eps_geom
        theta = np.cumsum(np.gradient(e.z) / np.maximum(beta, 1e-12))
        ax.plot(e.z, theta, label="$\\theta_%s$ [rad]" % lbl)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("phase advance [rad]")
    ax.set_title(title or "phase advance")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_coherence_length(emit, ax=None, figsize=(8, 4), title=None):
    """相干长度 Lc = hbar * sigma_x / (m0 c * eps_n) (菜单 2 项 15)."""
    fig, ax = _ax(ax, figsize)
    e = emit.x
    lc = HBAR_EVS * e.rms / (M0_C2_EV * np.maximum(e.emit, 1e-30)) * C_LIGHT
    ax.plot(e.z, lc, label="$L_c$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("coherence length [m]")
    ax.set_title(title or "coherence length (electrons)")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_phase_scan(pscan, ax=None, figsize=(8, 5), title=None):
    """能量增益 vs 腔相位 (PScan, 菜单 2 项 1)."""
    fig, ax = _ax(ax, figsize)
    ax.plot(pscan["phase_deg"], pscan["E_kin_eV"] * 1e-6, label="$E_{kin}$")
    ax.set_xlabel("RF phase [deg]")
    ax.set_ylabel("kinetic energy [MeV]")
    ax.set_title(title or "energy gain vs phase")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_scan_fom(scan, i=0, ax=None, figsize=(8, 4), title=None, lab=None):
    """参数扫描 FOM(i) (Scan 文件, 菜单 3).

    lab: read_lab_file 的输出; 有定义时用真实轴名/标题。
    """
    fig, ax = _ax(ax, figsize)
    ax.plot(scan["para"], scan["FOM"][:, i], label="FOM(%d)" % (i + 1))
    if lab is not None and i < len(lab["xlabel"]):
        xl = lab["xlabel"][i]
        yl = lab["ylabel"][i]
        tl = lab["title"][i]
        if xl and xl.lower() != "no entry":
            ax.set_xlabel(xl)
        else:
            ax.set_xlabel("scan parameter")
        if yl and yl.lower() != "no entry":
            ax.set_ylabel(yl)
        else:
            ax.set_ylabel("FOM(%d)" % (i + 1))
        ax.set_title(title or (tl if tl and tl.lower() != "no entry" else "parameter scan"))
    else:
        ax.set_xlabel("scan parameter")
        ax.set_ylabel("FOM(%d)" % (i + 1))
        ax.set_title(title or "parameter scan")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_error_hist(err, i=0, bins=20, ax=None, figsize=(8, 4), title=None,
                     lab=None):
    """误差扫描 FOM(i) 直方图 (Error 文件, 菜单 3).

    lab: read_lab_file 的输出; 有定义时用真实轴名/标题。
    """
    fig, ax = _ax(ax, figsize)
    ax.hist(err["FOM"][:, i], bins=bins, alpha=0.75, label="FOM(%d)" % (i + 1))
    if lab is not None and i < len(lab["xlabel"]):
        xl = lab["xlabel"][i]
        tl = lab["title"][i]
        ax.set_xlabel(xl if xl and xl.lower() != "no entry" else "FOM(%d)" % (i + 1))
        ax.set_title(title or (tl if tl and tl.lower() != "no entry" else "error scan histogram"))
    else:
        ax.set_xlabel("FOM(%d)" % (i + 1))
        ax.set_title(title or "error scan histogram")
    ax.set_ylabel("counts")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_reduced_emittance(x2, y2=None, ax=None, figsize=(8, 5), title=None):
    """缩减发射度 (Xemit2, 菜单 4 项 1-3): eps_red_z 与 eps_red_zE vs z."""
    fig, ax = _ax(ax, figsize)
    ax.plot(x2["z"], x2["eps_red_z"] * 1e6, label="$\\varepsilon_{red,z}$ [x]")
    ax.plot(x2["z"], x2["eps_red_zE"] * 1e6, label="$\\varepsilon_{red,zE}$ [x]")
    if y2 is not None:
        ax.plot(y2["z"], y2["eps_red_z"] * 1e6, ls="--", label="$\\varepsilon_{red,z}$ [y]")
        ax.plot(y2["z"], y2["eps_red_zE"] * 1e6, ls="--", label="$\\varepsilon_{red,zE}$ [y]")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("reduced emittance [$\\pi$ mm mrad]")
    ax.set_title(title or "reduced emittance (correlations removed)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def plot_trace_emittance(tr, ax=None, figsize=(8, 5), title=None):
    """trace-space 发射度 (TRemit, 菜单 4 项 7/8)."""
    fig, ax = _ax(ax, figsize)
    ax.plot(tr["z"], tr["eps_tr_x"] * 1e6, label="$\\varepsilon_{tr,x}$")
    ax.plot(tr["z"], tr["eps_tr_y"] * 1e6, label="$\\varepsilon_{tr,y}$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("trace-space emittance [$\\pi$ mm mrad]")
    ax.set_title(title or "trace-space emittance")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_larmor(lm, ax=None, figsize=(8, 4), title=None):
    """拉莫尔角平均与 RMS (Larmor 文件)."""
    fig, ax = _ax(ax, figsize)
    ax.plot(lm["z"], lm["avr"], label="average Larmor angle")
    ax.plot(lm["z"], lm["rms"], label="rms Larmor angle")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("angle [rad]")
    ax.set_title(title or "Larmor angle")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_probe_trajectories(path, ax=None, figsize=(8, 5), title=None):
    """探针轨迹 x(z), y(z) (track 文件, 菜单 1 项 14/15)."""
    tr = read_track_file(path) if not isinstance(path, dict) else path
    fig, ax = _ax(ax, figsize)
    for seq in np.unique(tr["seq"]):
        m = tr["seq"] == seq
        ax.plot(tr["z"][m], tr["x"][m] * 1e3, lw=0.8)
        ax.plot(tr["z"][m], tr["y"][m] * 1e3, lw=0.8, ls="--")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x / y [mm]")
    ax.set_title(title or "probe trajectories (x solid, y dashed)")
    fig.tight_layout()
    return fig


def plot_space_charge_fields(path, component="Ez", ax=None, figsize=(8, 5), title=None):
    """探针上的空间电荷场 Ez(z)/Er(z) (track 文件, fieldplot 子菜单近似)."""
    tr = read_track_file(path) if not isinstance(path, dict) else path
    fig, ax = _ax(ax, figsize)
    key = "Ez" if component == "Ez" else "Er"
    for seq in np.unique(tr["seq"]):
        m = tr["seq"] == seq
        ax.plot(tr["z"][m], tr[key][m], lw=0.8)
    ax.set_xlabel("z [m]")
    ax.set_ylabel(component + " [V/m]")
    ax.set_title(title or ("space charge field %s on probe particles" % component))
    fig.tight_layout()
    return fig


def plot_cathode_emission(path, ax=None, figsize=(8, 5), title=None,
                           include_spch: bool = True):
    """阴极发射: 加速场 + 空间电荷场(阴极表面) 与发射电荷 vs 时间
    (Cathode 文件, 菜单 4 项 15/16 与 fieldplot 阴极表面场)."""
    c = read_cathode_file(path) if not isinstance(path, dict) else path
    fig, ax = _ax(ax, figsize)
    ax.plot(c["t"] * 1e9, c["E_acc"], label="$E_{acc}$ on cathode")
    if include_spch and "E_spch" in c:
        ax.plot(c["t"] * 1e9, c["E_spch"], color="C2", ls="--",
                label="$E_{spch}$ on cathode")
    ax.set_xlabel("t [ns]")
    ax.set_ylabel("cathode field [V/m]")
    ax2 = ax.twinx()
    ax2.plot(c["t"] * 1e9, c["q"], color="C1", label="emitted charge")
    ax2.set_ylabel("charge [nC]")
    ax.set_title(title or "cathode emission")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=9)
    fig.tight_layout()
    return fig


def slice_mismatch(dist, n_slices=20, bz_on_axis_T: float = 0.0):
    """slice 失配参数 (postpro 5.6.3 项 2):

    zeta_i = 0.5 * (beta0*gamma_i - 2*alpha0*alpha_i + gamma0*beta_i) >= 1

    投影束团 (下标 0) 与各纵向 slice (下标 i) 的 Courant-Snyder 参数。
    散角用正则动量 (manual 4.13.1, bz_on_axis_T)。
    返回 (z_centers [m], zeta_x, zeta_y)。
    """
    from ..analysis.emittance import compute_twiss_parameters
    from ..analysis.slices import compute_slice_analysis

    d = dist.filter_active()
    p_ref = dist.ref_momentum_or_mean()
    sa = compute_slice_analysis(dist, n_slices=n_slices, bz_on_axis_T=bz_on_axis_T)

    # 投影 Twiss (x / y), 使用正则散角 x' = p~x/p_ref
    ptx = d.px + 0.5 * C_LIGHT * bz_on_axis_T * d.y
    pty = d.py - 0.5 * C_LIGHT * bz_on_axis_T * d.x
    xp = (ptx - np.mean(ptx)) / p_ref
    yp = (pty - np.mean(pty)) / p_ref
    b0x, a0x, g0x = compute_twiss_parameters(d.x - np.mean(d.x), xp)
    b0y, a0y, g0y = compute_twiss_parameters(d.y - np.mean(d.y), yp)

    zeta_x = np.full(sa.n_slices, np.nan)
    zeta_y = np.full(sa.n_slices, np.nan)
    z = sa.z_centers
    edges = sa.z_edges
    for i in range(sa.n_slices):
        if sa.n_particles[i] < 3:
            continue
        mask = (d.z >= edges[i]) & (d.z < edges[i + 1])
        if i == sa.n_slices - 1:
            mask = (d.z >= edges[i]) & (d.z <= edges[i + 1])
        xi = d.x[mask]
        yi = d.y[mask]
        if len(xi) < 3:
            continue
        # 与投影分支同口径: slice 散角也用正则动量 (手册 4.13.1, 批 2 修复)
        pxi = (d.px[mask] + 0.5 * C_LIGHT * bz_on_axis_T * yi)
        pyi = (d.py[mask] - 0.5 * C_LIGHT * bz_on_axis_T * xi)
        bxi, axi, gxi = compute_twiss_parameters(xi - np.mean(xi), (pxi - np.mean(pxi)) / p_ref)
        byi, ayi, gyi = compute_twiss_parameters(yi - np.mean(yi), (pyi - np.mean(pyi)) / p_ref)
        zeta_x[i] = 0.5 * (b0x * gxi - 2 * a0x * axi + g0x * bxi)
        zeta_y[i] = 0.5 * (b0y * gyi - 2 * a0y * ayi + g0y * byi)
    return z, zeta_x, zeta_y


def plot_slice_mismatch(dist, n_slices=20, ax=None, figsize=(8, 4), title=None):
    """slice 失配参数图 (postpro 5.6.3 项 2)。"""
    fig, ax = _ax(ax, figsize)
    z, zx, zy = slice_mismatch(dist, n_slices=n_slices)
    ax.plot(z * 1e3, zx, label="$\\zeta_x$")
    ax.plot(z * 1e3, zy, label="$\\zeta_y$")
    ax.axhline(1.0, color="r", ls=":", lw=1, label="matched")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("mismatch parameter $\\zeta$")
    ax.set_title(title or "slice mismatch parameter")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_3d_map_slices(path, axis="z", n_slices=3, figsize=(13, 4), title=None,
                       unit=""):
    """3D 场图截面 (fieldplot 菜单 2): 沿 axis 取 n_slices 个切面.

    Args:
        path: 3D 场图文件 (3D_test.ex / 3D_Dipole.bx ...).
        axis: 切片方向 'x'/'y'/'z' (垂直于切面的轴).
        n_slices: 切面数 (等距).
        unit: 色条单位 (如 'V/m', 'T').
    """
    from ..io.field_map import read_3d_field_map
    x, y, z, f = read_3d_field_map(path)
    shape = f.shape
    fig, axs = plt.subplots(1, n_slices, figsize=figsize, squeeze=False)
    axs = axs[0]
    vmax = float(np.max(np.abs(f)))
    if vmax == 0:
        vmax = 1.0
    for k in range(n_slices):
        i = int(round(k * (shape[{"x": 0, "y": 1, "z": 2}[axis]] - 1) / max(n_slices - 1, 1)))
        ax = axs[k]
        if axis == "z":
            im = ax.pcolormesh(x * 1e3, y * 1e3, f[:, :, i].T, cmap="RdBu_r",
                               vmin=-vmax, vmax=vmax, shading="auto")
            ax.set_xlabel("x [mm]")
            ax.set_ylabel("y [mm]")
            ax.set_title("z = %.4g m" % z[i])
        elif axis == "y":
            im = ax.pcolormesh(x * 1e3, z * 1e3, f[:, i, :].T, cmap="RdBu_r",
                               vmin=-vmax, vmax=vmax, shading="auto")
            ax.set_xlabel("x [mm]")
            ax.set_ylabel("z [mm]")
            ax.set_title("y = %.4g mm" % (y[i] * 1e3))
        else:
            im = ax.pcolormesh(y * 1e3, z * 1e3, f[i, :, :].T, cmap="RdBu_r",
                               vmin=-vmax, vmax=vmax, shading="auto")
            ax.set_xlabel("y [mm]")
            ax.set_ylabel("z [mm]")
            ax.set_title("x = %.4g mm" % (x[i] * 1e3))
    fig.colorbar(im, ax=axs, label=unit or "field")
    if title:
        fig.suptitle(title)
    fig.subplots_adjust(wspace=0.35, bottom=0.15, top=0.88)
    return fig


def plot_pscan_dedz(pscan, ax=None, figsize=(8, 4), title=None):
    """dE/dz vs 相位 (PScan 数值微分, 菜单 2 项 2; 正比于关联能散)."""
    fig, ax = _ax(ax, figsize)
    dedz = np.gradient(pscan["E_kin_eV"], np.radians(pscan["phase_deg"]))
    ax.plot(pscan["phase_deg"], dedz * 1e-6, label="dE/dphase")
    ax.set_xlabel("RF phase [deg]")
    ax.set_ylabel("dE/dphase [MeV/rad]")
    ax.set_title(title or "dE/dz vs phase (proportional to corr. energy spread)")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_pscan_compression(pscan, ax=None, figsize=(8, 4), title=None):
    """压缩因子 vs 相位 (PScan 列 3, 菜单 2 项 3)."""
    fig, ax = _ax(ax, figsize)
    ax.plot(pscan["phase_deg"], pscan["compression"], label="compression factor")
    ax.set_xlabel("RF phase [deg]")
    ax.set_ylabel("compression factor")
    ax.set_title(title or "bunch compression vs phase")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_tcheck_scaling(tc, ax=None, figsize=(8, 5), title=None):
    """空间电荷缩放因子 vs z (tcheck 文件, 菜单 2 项 9/10)."""
    fig, ax = _ax(ax, figsize)
    names = ["nr(r)", "nr(z)", "nr(gamma)", "nz(r)", "nz(gamma*z)"]
    for i, nm in enumerate(names):
        ax.plot(tc["z"], tc["scaling"][:, i], label=nm)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("space charge scaling factor")
    ax.set_title(title or "space charge scaling factors")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_z_plot(dist, ax=None, figsize=(8, 5), title=None):
    """z-plot (postpro 5.6.1 项 10): 所有粒子 (含丢失) 沿束线的位置.

    按粒子序号着色: active 蓝、lost 红、passive 灰。
    """
    fig, ax = _ax(ax, figsize)
    idx = np.arange(dist.n_particle)
    ax.scatter(idx[dist.lost], dist.z[dist.lost] * 1e3, s=2, c="#CC3311",
               label="lost (%d)" % int(np.sum(dist.lost)))
    ax.scatter(idx[dist.passive], dist.z[dist.passive] * 1e3, s=2, c="0.6",
               label="passive (%d)" % int(np.sum(dist.passive)))
    ax.scatter(idx[dist.active], dist.z[dist.active] * 1e3, s=2, c="#0077BB",
               label="active (%d)" % int(np.sum(dist.active)))
    ax.set_xlabel("particle index")
    ax.set_ylabel("z [mm]")
    ax.set_title(title or "z-plot (all particles along the beamline)")
    ax.legend(markerscale=3, fontsize=9)
    fig.tight_layout()
    return fig


def plot_field_profile(path, label="field", unit="", scale=1.0, ax=None,
                       figsize=(8, 4), title=None):
    """一维场剖面 (z, F) 表: 四极梯度/二极场等 (fieldplot 菜单 1)."""
    data = np.loadtxt(path)
    z, f = data[:, 0], data[:, 1]
    fig, ax = _ax(ax, figsize)
    ax.plot(z * 1e3, f * scale, label=label)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel(label + (" [" + unit + "]" if unit else ""))
    ax.set_title(title or (label + " profile"))
    ax.legend()
    fig.tight_layout()
    return fig


def plot_curved_cathode_contour(path, ax=None, figsize=(8, 5), title=None):
    """弯曲阴极轮廓 (Contour.dat, fieldplot 菜单 1 项 8).

    表格式 (x, y, z, R) 沿轮廓; 画 R(z) 曲线。
    """
    data = np.loadtxt(path)
    fig, ax = _ax(ax, figsize)
    ax.plot(data[:, 2] * 1e3, data[:, 3] * 1e3, label="cathode contour")
    ax.fill_between(data[:, 2] * 1e3, 0, data[:, 3] * 1e3, alpha=0.2)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_aspect("equal")
    ax.set_title(title or "curved cathode contour")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_core_brightness(ce, landf=None, ax=None, figsize=(8, 4), title=None):
    """横向核心亮度 B = Q / (eps_nx * eps_ny) (postpro 5.6.1 项 8, 近似).

    Q 取 LandF 的束团电荷 (同 z 处); 无 LandF 时用 ce 的归一化亮度。
    """
    fig, ax = _ax(ax, figsize)
    z = ce["mean_z"]
    epsn = np.maximum(ce["norm_emit_x"] * ce["norm_emit_y"], 1e-30)
    if landf is not None and "landf_total_charge" in landf:
        # landf 来自 parse_output_file, 电荷已是 SI (C)
        q = np.interp(z, landf["landf_z"], landf["landf_total_charge"])
        b = np.abs(q) / epsn
        ylabel = "brightness B = Q/(eps_nx eps_ny) [C/m^2]"
    else:
        b = 1.0 / epsn
        ylabel = "1/(eps_nx eps_ny) [1/m^2]"
    ax.plot(z, b, label="brightness")
    ax.set_xlabel("z [m]")
    ax.set_ylabel(ylabel)
    ax.set_title(title or "transverse core brightness")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_slice_ellipses_3d(dist, n_slices=10, figsize=(9, 6), title=None,
                            bz_on_axis_T: float = 0.0):
    """3D RMS slice 椭圆 (postpro 5.6.3 项 6, mplot3d 静态).

    散角用正则动量 (manual 4.13.1, bz_on_axis_T)。
    """
    from ..analysis.emittance import compute_emittance_ellipse_params
    from ..analysis.slices import compute_slice_analysis

    sa = compute_slice_analysis(dist, n_slices=n_slices, bz_on_axis_T=bz_on_axis_T)
    d = dist.filter_active()
    p_ref = dist.ref_momentum_or_mean()
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    for i in range(sa.n_slices):
        if sa.n_particles[i] < 3:
            continue
        mask = (d.z >= sa.z_edges[i]) & (d.z < sa.z_edges[i + 1])
        if i == sa.n_slices - 1:
            mask = (d.z >= sa.z_edges[i]) & (d.z <= sa.z_edges[i + 1])
        xi = d.x[mask]
        yi = d.y[mask]
        ptx = d.px[mask] + 0.5 * C_LIGHT * bz_on_axis_T * yi
        xp = (ptx - np.mean(ptx)) / p_ref
        par = compute_emittance_ellipse_params(xi - np.mean(xi), xp)
        th = np.linspace(0, 2 * np.pi, 60)
        xe = par["a"] * np.cos(th) * 1e3
        ye = par["b"] * np.sin(th) * 1e3
        xr = xe * np.cos(par["theta"]) - ye * np.sin(par["theta"])
        yr = xe * np.sin(par["theta"]) + ye * np.cos(par["theta"])
        zc = sa.z_centers[i] * 1e3
        ax.plot(np.full_like(th, zc), xr, yr, lw=1.2)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("x [mm]")
    ax.set_zlabel("x' [mrad]")
    ax.set_title(title or "3D slice emittance ellipses")
    fig.tight_layout()
    return fig
def aperture_elements(namelist: dict, base_dir=None) -> list:
    """从 &APERTURE namelist (parse_namelists 输出) 提取孔径元素列表.

    返回 [{'z1','z2','r','xoff','yoff','file'}] (单位 m)。Ap_R 单位 mm
    (手册 6.2), 负值 = 圆柱堵块 (beam stop); 文件型孔径取几何表
    (z, R[mm]) 的最大半径。
    """
    import re as _re
    from pathlib import Path as _Path
    z1 = np.atleast_1d(namelist.get("Ap_Z1", []))
    z2 = np.atleast_1d(namelist.get("Ap_Z2", []))
    r = np.atleast_1d(namelist.get("Ap_R", []))
    xoff = np.atleast_1d(namelist.get("A_xoff", np.zeros(len(z1))))
    yoff = np.atleast_1d(namelist.get("A_yoff", np.zeros(len(z1))))
    if len(xoff) < len(z1):
        xoff = np.pad(xoff, (0, len(z1) - len(xoff)), constant_values=0.0)
    if len(yoff) < len(z1):
        yoff = np.pad(yoff, (0, len(z1) - len(yoff)), constant_values=0.0)
    files = namelist.get("File_Aperture", [])
    if isinstance(files, str):
        files = [files]
    els = []
    for i in range(len(z1)):
        el = {
            "z1": float(z1[i]), "z2": float(z2[i]),
            "r": float(r[i]) * 1e-3,          # mm -> m
            "xoff": float(xoff[i]),           # m
            "yoff": float(yoff[i]),
            "file": str(files[i]) if i < len(files) else "",
        }
        if el["file"] and not _re.fullmatch(r"(?i)rad|cir|scr_x|scr_y|col_x|col_y", el["file"]):
            if base_dir is not None:
                gp = _Path(base_dir) / el["file"]
                if gp.exists():
                    g = np.loadtxt(gp)
                    el["r"] = float(np.max(np.abs(g[:, 1]))) * 1e-3  # mm -> m
        els.append(el)
    return els


def plot_envelope_with_aperture(emit, apertures, ax=None, figsize=(9, 5),
                                title=None, plane: str = "x"):
    """束包络 + 孔径几何叠加 (postpro 含孔径显示).

    Args:
        emit: EmitSet (Xemit/Yemit)。
        apertures: aperture_elements() 输出的元素列表 (SI)。
        plane: 'x' 或 'y' (y 视图下孔径带按 yoff 平移)。
    """
    from matplotlib.patches import Rectangle
    fig, ax = _ax(ax, figsize)
    e = emit.x if plane == "x" else emit.y
    ax.plot(e.z, e.rms * 1e3, label="$\\sigma_%s$" % plane)
    ax.plot(e.z, -e.rms * 1e3, alpha=0.35, ls="--", label="-$\\sigma_%s$" % plane)
    for ap in apertures:
        z1, z2, r = ap["z1"], ap["z2"], ap["r"]
        if r < 0:  # beam stop / plug
            ax.axvspan(z1, z2, color="k", alpha=0.25, label="beam stop")
            continue
        off = ap["xoff"] if plane == "x" else ap["yoff"]
        ax.add_patch(Rectangle(
            (z1, (off - r) * 1e3), z2 - z1, 2 * r * 1e3,
            fill=False, edgecolor="r", lw=1.5, label="aperture"))
    ax.set_xlabel("z [m]")
    ax.set_ylabel("beam size / aperture [mm]")
    ax.set_title(title or "beam envelope with aperture geometry")
    handles, labels = ax.get_legend_handles_labels()
    by = dict(zip(labels, handles))
    ax.legend(by.values(), by.keys(), fontsize=9)
    fig.tight_layout()
    return fig


def plot_laser_on_axis(path, unit="a.u.", figsize=(11, 4), title=None):
    """激光 3D 图 (File_A0 格式) 的轴上剖面 vs z 与 vs t (fieldplot 8 章).

    取 x=0, y=0 网格点; 右图为随光运动的观察者时间 t = (z - z0)/c。
    数值原样 (归一化到 E_a0 由 ASTRA 负责)。
    """
    from ..io.field_map import read_3d_field_map
    x, y, z, f = read_3d_field_map(path)
    ix0 = int(np.argmin(np.abs(x)))
    iy0 = int(np.argmin(np.abs(y)))
    onax = f[ix0, iy0, :]
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    axes[0].plot(z * 1e3, onax)
    axes[0].set_xlabel("z [mm]")
    axes[0].set_ylabel("on-axis value [%s]" % unit)
    axes[0].set_title("on-axis profile vs z")
    t_ps = (z - z[0]) / C_LIGHT * 1e12
    axes[1].plot(t_ps, onax)
    axes[1].set_xlabel("t = (z - z0)/c [ps]")
    axes[1].set_ylabel("on-axis value [%s]" % unit)
    axes[1].set_title("co-moving time profile")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_plasma_profile(path, peak_density_cm3=None, ax=None, figsize=(8, 5),
                        title=None):
    """等离子体密度剖面 (File_Efield='Plasma...' 两列表, 手册 6.7).

    文件: z [m], n [arb. u.]; ASTRA 归一化到峰值 P_n。若给出
    peak_density_cm3 (=P_n), 右轴显示物理密度 [cm^-3]。
    """
    fig, ax = _ax(ax, figsize)
    d = np.loadtxt(path, ndmin=2)
    z, n = d[:, 0], d[:, 1]
    ax.plot(z * 1e3, n, label="profile (normalized)")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("plasma density [arb. u.]")
    ax.set_title(title or "plasma density profile")
    if peak_density_cm3:
        ax2 = ax.twinx()
        ax2.plot(z * 1e3, n / np.max(np.abs(n)) * peak_density_cm3,
                 color="C1", ls="--", label="scaled to P_n")
        ax2.set_ylabel("density [cm$^{-3}$]", color="C1")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=9)
    else:
        ax.legend()
    fig.tight_layout()
    return fig


def plot_core_fraction_curves(dist, fractions=(0.1, 0.25, 0.5, 0.75, 0.9, 1.0),
                              figsize=(10, 4), title=None):
    """核心束长/发射度 vs 电荷分数 (postpro 核心曲线, 自算)."""
    from ..analysis.core import compute_core_fraction_curves
    c = compute_core_fraction_curves(dist, fractions=fractions)
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    axes[0].plot(c["fractions"] * 100, c["sig_z"] * 1e3, "o-", label="$\\sigma_z$")
    axes[0].plot(c["fractions"] * 100, c["sig_x"] * 1e3, "s-", label="$\\sigma_x$")
    axes[0].plot(c["fractions"] * 100, c["sig_y"] * 1e3, "^-", label="$\\sigma_y$")
    axes[0].set_xlabel("core charge fraction [%]")
    axes[0].set_ylabel("RMS size [mm]")
    axes[0].set_title("core bunch length / sizes")
    axes[0].legend(fontsize=9)
    axes[1].plot(c["fractions"] * 100, c["emit_xn"] * 1e6, "o-",
                 label="$\\varepsilon_{nx}$")
    axes[1].plot(c["fractions"] * 100, c["emit_yn"] * 1e6, "s-",
                 label="$\\varepsilon_{ny}$")
    axes[1].set_xlabel("core charge fraction [%]")
    axes[1].set_ylabel("norm. emittance [$\\pi$ mm mrad]")
    axes[1].set_title("core emittance")
    axes[1].legend(fontsize=9)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig

def plot_pscan_compression_time(pscan, ax=None, figsize=(8, 4), title=None):
    """压缩因子 (时间, lineplot 菜单 2 项 4): PScan 第 4 列 beta/beta0。"""
    fig, ax = _ax(ax, figsize)
    ax.plot(pscan["phase_deg"], pscan["beta_ratio"], label=r"$\beta/\beta_0$")
    ax.set_xlabel("RF phase [deg]")
    ax.set_ylabel(r"velocity ratio $\beta/\beta_0$")
    ax.set_title(title or "compression factor (time)")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_scan_position(scan, ax=None, figsize=(8, 4), title=None):
    """FOM 保存位置 (lineplot 菜单 3 项 11): z vs 扫描参数。"""
    fig, ax = _ax(ax, figsize)
    ax.plot(scan["para"], scan["z"], "o-", label="FOM position")
    ax.set_xlabel("scan parameter")
    ax.set_ylabel("z [m]")
    ax.set_title(title or "FOM evaluation position")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_tcheck_counter(tc, ax=None, figsize=(8, 4), title=None):
    """空间电荷缩放计数器 (lineplot 菜单 2 项 10)。"""
    fig, ax = _ax(ax, figsize)
    ax.plot(tc["z"], tc["counter"], label="scaling counter")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("counter")
    ax.set_title(title or "space charge scaling counter")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_core_emittance(ce, ax=None, figsize=(8, 5), title=None,
                        plane: str = "x"):
    """核心发射度 (Cemit, 菜单 4 项 9-11): eps_n + C95/C90/C80。

    plane: 'x' / 'y' / 'z' (z 平面纵向, 单位 keV.mm)。
    """
    fig, ax = _ax(ax, figsize)
    z = ce["mean_z"]
    # Cemit z 列内部单位 eV·m (OUTPUT_TABLES 因子 1), 1 keV·mm ≡ 1 eV·m,
    # 显示 keV·mm 的缩放系数为 1 (此前 1e-3 使值小 1000 倍 — 批 2 修复)。
    scale = 1.0 if plane == "z" else 1e6
    unit = "keV mm" if plane == "z" else r"$\pi$ mm mrad"
    ax.plot(z, ce["norm_emit_" + plane] * scale, label=r"$\varepsilon_{n%s}$" % plane)
    ax.plot(z, ce["core_emit_95percent_" + plane] * scale, ls="--",
            label="C%s_95" % plane)
    ax.plot(z, ce["core_emit_90percent_" + plane] * scale, ls=":",
            label="C%s_90" % plane)
    ax.plot(z, ce["core_emit_80percent_" + plane] * scale, ls="-.",
            label="C%s_80" % plane)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("core emittance [%s]" % unit)
    ax.set_title(title or ("core emittance (%s)" % plane))
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def plot_cr_emit(cr, ax=None, figsize=(8, 5), title=None):
    """交叉粒子 (Cr_emit, 菜单 4 项 12-14): 发射度与剩余/交叉电荷。"""
    fig, ax = _ax(ax, figsize)
    ax.plot(cr["z"], cr["eps_x"] * 1e6, label=r"$\varepsilon_x$")
    ax.plot(cr["z"], cr["eps_y"] * 1e6, label=r"$\varepsilon_y$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel(r"emittance [$\pi$ mm mrad]")
    ax.set_title(title or "cross-over particle emittance")
    ax2 = ax.twinx()
    ax2.plot(cr["z"], cr["q_rest"], color="C2", ls="--", label="rest charge")
    ax2.plot(cr["z"], cr["q_cross"], color="C3", ls=":", label="cross charge")
    ax2.set_ylabel("charge [nC]")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=9)
    fig.tight_layout()
    return fig
