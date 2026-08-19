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

import warnings

import matplotlib.pyplot as plt
import numpy as np

from ..constants import C_LIGHT, M_E_C2_EV, E_CHARGE
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
    for e, lbl in ((emit.x, "x"), (emit.y, "y")):
        # beta*gamma 逐平面按各自 z 网格插值 (X/Yemit 网格可能不一致)
        bg = np.interp(e.z, ref.z, np.maximum(ref.pz / M_E_C2_EV, 0.0))
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
    for e, lbl in ((emit.x, "x"), (emit.y, "y")):
        # beta*gamma 逐平面按各自 z 网格插值 (X/Yemit 网格可能不一致)
        bg = np.interp(e.z, ref.z, np.maximum(ref.pz / M_E_C2_EV, 0.0))
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


def plot_emittance_difference(emit, x2, y2=None, ax=None, figsize=(8, 4),
                              title=None):
    """标准与缩减横向发射度之差 (菜单 4 项 3, 手册 4.13.6).

    diff = eps_std - eps_red_z (z 相关缩减)。emit: Xemit 文件
    (norm_emit_x/y [m.rad]); x2/y2: Xemit2/Yemit2 (eps_red_z [m.rad])。
    按 z 线性插值到 emit 网格。y 平面缺 Yemit2 时用 x 近似 (虚线标注 est.)。
    """
    from numpy import interp
    # emit: EmitSet (Xemit/Yemit/Zemit) 或 dict 兼容
    if hasattr(emit, "x") and hasattr(emit.x, "emit"):
        z = np.asarray(emit.x.z)
        ex = np.asarray(emit.x.emit)
        ey = np.asarray(emit.y.emit)
    else:
        z = np.asarray(emit["z"])
        ex = np.asarray(emit["norm_emit_x"])
        ey = np.asarray(emit["norm_emit_y"])
    fig, ax = _ax(ax, figsize)
    red_x = interp(z, x2["z"], x2["eps_red_z"])
    ax.plot(z, (ex - red_x) * 1e6, label=r"$\Delta\varepsilon_x$")
    if y2 is not None:
        red_y = interp(z, y2["z"], y2["eps_red_z"])
        ax.plot(z, (ey - red_y) * 1e6, ls="--",
                label=r"$\Delta\varepsilon_y$")
    else:
        ax.plot(z, (ey - red_x) * 1e6, ls="--",
                label=r"$\Delta\varepsilon_y$ (est.)")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("emittance difference [$\\pi$ mm mrad]")
    ax.set_title(title or "standard - reduced emittance (transverse)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def plot_correlated_emittance_contributions(x2, y2=None, ax=None,
                                            figsize=(8, 4), title=None):
    """相关发射度贡献 (Eq. 4.4 K 项, 菜单 4 项 4/5).

    画 K2z/K3z/K2E/K3E 绝对值 vs z [pi mm mrad], 对应纵向位置与动能
    相关的线性/二次相关贡献 (手册 4.13.6 Eq. 4.4)。x2=Xemit2 (x),
    y2=Yemit2 (y, 虚线)。
    """
    fig, ax = _ax(ax, figsize)
    for key, lab, c in [("K2z", r"$K_{2,Z}$", "C0"),
                        ("K3z", r"$K_{3,Z}$", "C1"),
                        ("K2E", r"$K_{2,E}$", "C2"),
                        ("K3E", r"$K_{3,E}$", "C3")]:
        ax.plot(x2["z"], np.abs(x2[key]) * 1e6, label=lab + " [x]", color=c)
        if y2 is not None:
            ax.plot(y2["z"], np.abs(y2[key]) * 1e6, ls="--", color=c,
                    label=lab + " [y]")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("correlated contribution [$\\pi$ mm mrad]")
    ax.set_title(title or "correlated emittance contributions (Eq. 4.4)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def plot_reduced_longitudinal_emittance(dist, ax=None, figsize=(8, 5),
                                        title=None):
    """缩减纵向发射度 (菜单 4 项 6, 手册 4.13.6).

    纵向发射度减去 2nd/3rd 阶 z-pz 多项式相关, 连同正则相空间发射度。
    显示 z-(pz/p_ref) 相空间散点 + 3 阶相关拟合线, 框内标注两个值:
    eps_z (正则) 与 eps_z_red (去 2nd/3rd 阶相关) [mm mrad]。
    """
    from ..analysis.emittance import compute_geometric_emittance
    act = dist.filter_active()
    z = act.z
    pz = act.pz                                      # eV/c
    pref = float(np.median(np.abs(pz))) or 1.0
    u = z - np.mean(z)
    up = (pz - np.mean(pz)) / pref                   # 相对动量发散 [rad]
    w = np.abs(act.charge)
    eps_z = compute_geometric_emittance(u, up, w)
    coeffs = np.polyfit(z, up, 3)                    # 2nd/3rd 阶相关
    up_red = up - np.polyval(coeffs, z)
    eps_z_red = compute_geometric_emittance(u, up_red, w)
    fig, ax = _ax(ax, figsize)
    ax.scatter(z * 1e3, up * 1e3, s=4, label="pz/p_ref")
    zs = np.linspace(z.min(), z.max(), 100)
    ax.plot(zs * 1e3, np.polyval(coeffs, zs) * 1e3, color="C3", lw=1.2,
            label="3rd-order correlation fit")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("(pz - <pz>) / p_ref [mrad]")
    ax.set_title(title or "longitudinal emittance & reduced "
                          "(2nd/3rd-order corr. removed)")
    ax.text(0.02, 0.98,
            "$\\varepsilon_z$ = %.3f mm.mrad\n"
            "$\\varepsilon_{z,red}$ = %.3f mm.mrad" % (eps_z * 1e6,
                                                       eps_z_red * 1e6),
            transform=ax.transAxes, va="top", fontsize=10,
            bbox=dict(boxstyle="round", fc="white", alpha=0.85))
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_trace_emittance(tr, ax=None, figsize=(8, 5), title=None):
    """trace-space 发射度 (TRemit, 菜单 4 项 7/8).

    横向 eps_tr_x/y [pi mm mrad] (左轴) + 纵向 eps_tr_z [pi um] (右轴,
    TRemit 纵向列, 覆盖项 8)。
    """
    fig, ax = _ax(ax, figsize)
    ax.plot(tr["z"], tr["eps_tr_x"] * 1e6, label="$\\varepsilon_{tr,x}$")
    ax.plot(tr["z"], tr["eps_tr_y"] * 1e6, label="$\\varepsilon_{tr,y}$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("trace-space emittance [$\\pi$ mm mrad]")
    ax.set_title(title or "trace-space emittance")
    if "eps_tr_z" in tr:
        ax2 = ax.twinx()
        ax2.plot(tr["z"], tr["eps_tr_z"] * 1e6, ls="--", color="C3",
                 label="$\\varepsilon_{tr,z}$")
        ax2.set_ylabel("long. trace-space emit [$\\pi$ um]", color="0.3")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=9)
    else:
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


def plot_probe_trajectories(path, mode="cartesian", ax=None, figsize=None,
                            title=None):
    """探针轨迹 (track 文件, 菜单 1 项 14/15).

    mode='cartesian': 单面板 x(z)/y(z) (实/虚线, 项 14)。
    mode='cylindrical': 双面板 r(z)=sqrt(x²+y²) 与 x/y 投影 (项 15)。
    """
    tr = read_track_file(path) if not isinstance(path, dict) else path
    if mode == "cylindrical":
        fig, axes = plt.subplots(1, 2, figsize=figsize or (11, 4))
        for seq in np.unique(tr["seq"]):
            m = tr["seq"] == seq
            r = np.hypot(tr["x"][m], tr["y"][m])
            axes[0].plot(tr["z"][m], r * 1e3, lw=0.8)
            axes[1].plot(tr["z"][m], tr["x"][m] * 1e3, lw=0.8)
            axes[1].plot(tr["z"][m], tr["y"][m] * 1e3, lw=0.8, ls="--")
        axes[0].set_xlabel("z [m]")
        axes[0].set_ylabel("r [mm]")
        axes[0].set_title("cylindrical radius r(z)")
        axes[1].set_xlabel("z [m]")
        axes[1].set_ylabel("x / y [mm]")
        axes[1].set_title("x / y projections")
        if title:
            fig.suptitle(title)
        fig.tight_layout()
        return fig
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


def plot_quadrupole_field(path, figsize=(12, 8), title=None):
    """四极场 (fieldplot 菜单 5 + next page).

    主图: 水平/垂直梯度 Gx(z), Gy(z) (三列文件 z, Gx, Gy; 两列文件
    z, G 用理想四极 Gy=-Gx)。
    next page: 纵向磁场 Bz(z) + 综合图。ASTRA 四极通常用解析聚焦
    (无场表), 纵向边缘场未建模: 文件提供第 4 列才画 Bz, 否则标注
    理想四极 Bz=0。2x2: [Gx/Gy 主图, Bz(z) / 理想标注],
                        [综合图 (各分量叠加), 说明文本]。
    """
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] < 2:
        raise ValueError("quadrupole field file needs >= 2 columns")
    z = data[:, 0]
    if data.shape[1] >= 3:
        gx, gy = data[:, 1], data[:, 2]
    else:
        g = data[:, 1]
        gx, gy = g, -g
    bz = data[:, 3] if data.shape[1] >= 4 else None

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    # 主图: 水平/垂直梯度
    axes[0, 0].plot(z * 1e3, gx, label="Gx (hor.)")
    axes[0, 0].plot(z * 1e3, gy, label="Gy (ver.)")
    axes[0, 0].set_xlabel("z [mm]")
    axes[0, 0].set_ylabel("gradient [T/m]")
    axes[0, 0].set_title("quadrupole gradients")
    axes[0, 0].legend(fontsize=9)
    # next page: Bz
    if bz is not None:
        axes[0, 1].plot(z * 1e3, bz, color="C3")
        axes[0, 1].set_ylabel("Bz [T]")
        axes[0, 1].set_title("longitudinal field Bz")
    else:
        axes[0, 1].text(0.5, 0.5, "ideal quadrupole: Bz = 0\n"
                        "(longitudinal edge field not modeled)",
                        ha="center", va="center", transform=axes[0, 1].transAxes)
        axes[0, 1].set_title("longitudinal field Bz (ideal)")
        axes[0, 1].set_xticks([])
        axes[0, 1].set_yticks([])
    axes[0, 1].set_xlabel("z [mm]")
    # next page: 综合图
    axes[1, 0].plot(z * 1e3, gx, label="Gx")
    axes[1, 0].plot(z * 1e3, gy, label="Gy")
    if bz is not None:
        axes[1, 0].plot(z * 1e3, bz, color="C3", label="Bz")
    axes[1, 0].set_xlabel("z [mm]")
    axes[1, 0].set_ylabel("field [T/m or T]")
    axes[1, 0].set_title("combined plot of all components")
    axes[1, 0].legend(fontsize=9)
    axes[1, 1].axis("off")
    axes[1, 1].text(0.02, 0.98,
                    "main view: horizontal & vertical gradient\n"
                    "next page: longitudinal field Bz + combined plot\n"
                    "quadrupole edge (fringe) fields are not stored in\n"
                    "the deck; provide a 4-column field table to show Bz.",
                    va="top", fontsize=9, transform=axes[1, 1].transAxes)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_curved_cathode_contour(path, ax=None, figsize=(8, 5), title=None,
                                show_rings=False, ring_offset=None):
    """弯曲阴极轮廓 (Contour.dat, fieldplot 菜单 1 项 8).

    手册 4.4.5: 前两列为纵向 z 与径向 R 坐标, 第三/四列为阴极在该点
    的切向单位矢量分量 (t_z, t_R)。这里画 z(R) 截面轮廓。

    show_rings: 同时画出电荷环位置。每个表点后方一点放置一个电荷环
    (用于修正非平面阴极的镜象电荷场, 手册 4.4.5)。环位于表面法向
    (指向阴极背面) 偏移 ring_offset [m] 处 (默认 10% 纵向跨度);
    图为示意图, 偏移量可调。
    """
    data = np.loadtxt(path)
    z, r = data[:, 0], data[:, 1]
    tz, tr = data[:, 2], data[:, 3]
    fig, ax = _ax(ax, figsize)
    ax.plot(z * 1e3, r * 1e3, label="cathode contour")
    ax.fill_between(z * 1e3, 0, r * 1e3, alpha=0.2)
    if show_rings:
        if ring_offset is None:
            ring_offset = 0.1 * float(np.ptp(z))
        # 发射侧法向 = (tr, -tz); 阴极背面 = (-tr, tz) (环所在)
        nz, nr = -tr, tz
        zr = z + ring_offset * nz
        rr = r + ring_offset * nr
        ax.plot(zr * 1e3, rr * 1e3, "o", ms=3, color="C3",
                label="charge rings")
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
                            bz_on_axis_T: float = 0.0, plane="xxp",
                            subtract_corr=False):
    """3D RMS slice 椭圆 (postpro 5.6.3 项 6, mplot3d 静态).

    plane: 'xxp' (x-x', 默认) / 'yyp' (y-y') / 'xyp' (x-y') /
        'yxp' (y-x') (项 7 投影切换)。
    subtract_corr: 从散角减去对位置的线性相关 (项 11)。
    散角用正则动量 (manual 4.13.1, bz_on_axis_T)。
    """
    from ..analysis.emittance import compute_emittance_ellipse_params
    from ..analysis.slices import compute_slice_analysis
    from .arbitrary_phase_space import subtract_linear_corr

    if plane not in ("xxp", "yyp", "xyp", "yxp"):
        raise ValueError("plane 必须为 'xxp'/'yyp'/'xyp'/'yxp': %r" % (plane,))
    sa = compute_slice_analysis(dist, n_slices=n_slices, bz_on_axis_T=bz_on_axis_T)
    d = dist.filter_active()
    p_ref = dist.ref_momentum_or_mean()
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    lab = {"xxp": ("x [mm]", "x' [mrad]"),
           "yyp": ("y [mm]", "y' [mrad]"),
           "xyp": ("x [mm]", "y' [mrad]"),
           "yxp": ("y [mm]", "x' [mrad]")}[plane]
    for i in range(sa.n_slices):
        if sa.n_particles[i] < 3:
            continue
        mask = (d.z >= sa.z_edges[i]) & (d.z < sa.z_edges[i + 1])
        if i == sa.n_slices - 1:
            mask = (d.z >= sa.z_edges[i]) & (d.z <= sa.z_edges[i + 1])
        xi = d.x[mask]
        yi = d.y[mask]
        ptx = d.px[mask] + 0.5 * C_LIGHT * bz_on_axis_T * yi
        pty = d.py[mask] - 0.5 * C_LIGHT * bz_on_axis_T * xi
        xp = (ptx - np.mean(ptx)) / p_ref
        yp = (pty - np.mean(pty)) / p_ref
        if plane == "xxp":
            u, up = xi, xp
        elif plane == "yyp":
            u, up = yi, yp
        elif plane == "xyp":
            u, up = xi, yp
        else:
            u, up = yi, xp
        if subtract_corr:
            up = subtract_linear_corr(u - np.mean(u), up)
        par = compute_emittance_ellipse_params(u - np.mean(u), up)
        th = np.linspace(0, 2 * np.pi, 60)
        xe = par["a"] * np.cos(th) * 1e3
        ye = par["b"] * np.sin(th) * 1e3
        xr = xe * np.cos(par["theta"]) - ye * np.sin(par["theta"])
        yr = xe * np.sin(par["theta"]) + ye * np.cos(par["theta"])
        zc = sa.z_centers[i] * 1e3
        ax.plot(np.full_like(th, zc), xr, yr, lw=1.2)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel(lab[0])
    ax.set_zlabel(lab[1])
    ax.set_title(title or "3D slice emittance ellipses (plane=%s%s)"
                 % (plane, ", corr removed" if subtract_corr else ""))
    fig.tight_layout()
    return fig


def plot_slice_ellipses_2d(dist, n_slices=10, figsize=(8, 6), title=None,
                           bz_on_axis_T: float = 0.0, plane="xxp",
                           subtract_corr=False):
    """投影 rms slice 发射度椭圆 (postpro 5.6.3 项 4, 2D).

    plane: 'xxp' (x-x', 默认) / 'yyp' (y-y') / 'xyp' (x-y') /
        'yxp' (y-x') (项 7 投影切换)。
    subtract_corr: 从散角减去对位置的线性相关 (项 11)。
    散角用正则动量 (manual 4.13.1, bz_on_axis_T)。粒子散点 (灰) +
    各 slice 的 RMS 椭圆 (按 slice 索引着色)。
    """
    from ..analysis.emittance import compute_emittance_ellipse_params
    from ..analysis.slices import compute_slice_analysis
    from .arbitrary_phase_space import subtract_linear_corr

    if plane not in ("xxp", "yyp", "xyp", "yxp"):
        raise ValueError("plane 必须为 'xxp'/'yyp'/'xyp'/'yxp': %r" % (plane,))
    sa = compute_slice_analysis(dist, n_slices=n_slices, bz_on_axis_T=bz_on_axis_T)
    d = dist.filter_active()
    p_ref = dist.ref_momentum_or_mean()
    lab = {"xxp": ("x [mm]", "x' [mrad]"),
           "yyp": ("y [mm]", "y' [mrad]"),
           "xyp": ("x [mm]", "y' [mrad]"),
           "yxp": ("y [mm]", "x' [mrad]")}[plane]
    xi, yi = d.x, d.y
    ptx = d.px + 0.5 * C_LIGHT * bz_on_axis_T * yi
    pty = d.py - 0.5 * C_LIGHT * bz_on_axis_T * xi
    xp = (ptx - np.mean(ptx)) / p_ref
    yp = (pty - np.mean(pty)) / p_ref
    if plane == "xxp":
        u, up = xi, xp
    elif plane == "yyp":
        u, up = yi, yp
    elif plane == "xyp":
        u, up = xi, yp
    else:
        u, up = yi, xp
    if subtract_corr:
        up = subtract_linear_corr(u - np.mean(u), up)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(u * 1e3, up * 1e3, s=3, alpha=0.25, color="0.6",
               label="particles")
    cmap = plt.get_cmap("viridis")
    for i in range(sa.n_slices):
        if sa.n_particles[i] < 3:
            continue
        mask = (d.z >= sa.z_edges[i]) & (d.z < sa.z_edges[i + 1])
        if i == sa.n_slices - 1:
            mask = (d.z >= sa.z_edges[i]) & (d.z <= sa.z_edges[i + 1])
        xi2 = d.x[mask]
        yi2 = d.y[mask]
        ptx2 = d.px[mask] + 0.5 * C_LIGHT * bz_on_axis_T * yi2
        pty2 = d.py[mask] - 0.5 * C_LIGHT * bz_on_axis_T * xi2
        xp2 = (ptx2 - np.mean(ptx2)) / p_ref
        yp2 = (pty2 - np.mean(pty2)) / p_ref
        if plane == "xxp":
            ui, upi = xi2, xp2
        elif plane == "yyp":
            ui, upi = yi2, yp2
        elif plane == "xyp":
            ui, upi = xi2, yp2
        else:
            ui, upi = yi2, xp2
        if subtract_corr:
            upi = subtract_linear_corr(ui - np.mean(ui), upi)
        par = compute_emittance_ellipse_params(ui - np.mean(ui), upi)
        th = np.linspace(0, 2 * np.pi, 60)
        xe = par["a"] * np.cos(th) * 1e3
        ye = par["b"] * np.sin(th) * 1e3
        xr = xe * np.cos(par["theta"]) - ye * np.sin(par["theta"])
        yr = xe * np.sin(par["theta"]) + ye * np.cos(par["theta"])
        ax.plot(xr, yr, color=cmap(i / max(sa.n_slices - 1, 1)), lw=1.3)
    ax.set_xlabel(lab[0])
    ax.set_ylabel(lab[1])
    ax.set_title(title or "projected slice emittance ellipses (plane=%s%s)"
                 % (plane, ", corr removed" if subtract_corr else ""))
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(0, max(sa.n_slices - 1, 1)))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="slice index")
    ax.legend(fontsize=9)
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


def plot_laser_envelope(path, figsize=(8, 5), title=None):
    """激光 3D 图 (File_A0 格式) 的 rms 横向束包络 + 焦点位置 (5.7.3).

    每个 z 切片以幅值 a₀ = √f 为权重计算 rms 半径 σx(z), σy(z) ——
    laser.dat 存的是归一化矢势平方 a₀⊥² (手册 File_A0 语义,
    create_lasermap.m: asq = a0² w0²/w² exp(-2r²/w²)), 物理包络权重是
    a₀ 而非 f² (2026-08 审计 P1: f² 权重使 rms 系统性偏小 1/√2 ≈ 29%)。
    高斯束下每轴 rms = w/√2 (w = 1/e² 幅值半径), 焦点处
    √(σx²+σy²) = w0; 图中竖线标注焦点。
    """
    from ..io.field_map import read_3d_field_map
    x, y, z, f = read_3d_field_map(path)
    sx = np.empty(len(z))
    sy = np.empty(len(z))
    X, Y = np.meshgrid(x, y, indexing="ij")
    for k in range(len(z)):
        sl = f[:, :, k]
        w = np.sqrt(np.abs(sl))     # a₀ 幅值权重 (File_A0 存 a₀²)
        wsum = float(w.sum())
        if wsum <= 0:
            sx[k] = sy[k] = np.nan
            continue
        xm = float((w * X).sum() / wsum)
        ym = float((w * Y).sum() / wsum)
        sx[k] = float(np.sqrt((w * (X - xm) ** 2).sum() / wsum))
        sy[k] = float(np.sqrt((w * (Y - ym) ** 2).sum() / wsum))
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(z * 1e3, sx * 1e3, label="rms x envelope")
    ax.plot(z * 1e3, sy * 1e3, label="rms y envelope")
    sig = np.sqrt(sx ** 2 + sy ** 2)
    ok = np.isfinite(sig)
    if ok.any():
        kf = int(np.nanargmin(sig))
        ax.axvline(z[kf] * 1e3, color="C2", ls="--", lw=1)
        ax.annotate("focus: z=%.2f mm, sigma=%.3f mm"
                    % (z[kf] * 1e3, sig[kf] * 1e3),
                    xy=(z[kf] * 1e3, sig[kf] * 1e3),
                    xytext=(6, 6), textcoords="offset points", fontsize=9)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("rms spot size [mm]")
    ax.set_title(title or "laser rms beam envelope & focus position")
    ax.legend(fontsize=9)
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
        nmax = np.max(np.abs(n))
        if nmax == 0:
            nmax = 1.0
        ax2.plot(z * 1e3, n / nmax * peak_density_cm3,
                 color="C1", ls="--", label="scaled to P_n")
        ax2.set_ylabel("density [cm$^{-3}$]", color="C1")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=9)
    else:
        ax.legend()
    fig.tight_layout()
    return fig


def plot_plasma_fields(path, peak_density_cm3=None, a0=1.0, sigma_z_m=30e-6,
                       vs="z", figsize=(8, 5), title=None):
    """等离子体场 vs z / vs zeta (fieldplot 5.7.3 项 8/9, 手册 6.7 + 附录).

    读两列表 (z [m], n 归一化密度)。用线性等离子体尾场解析公式
    (手册附录, 轴上 r=0, 简化幅度) 重建纵向场:

      Ez(ζ) = kp² m_e c² a0² / (2e) · sqrt(π/2) · σz
              · exp(-kp² σz² / 2) · cos(kp ζ)

    kp = sqrt(n_peak e²/(ε0 m_e c²)) 由峰值密度 (peak_density_cm3)
    确定; ζ = z - c t 为共动参数 (粒子随光速运动看到)。

    vs='z': 横轴 z [mm], 左轴 Ez [a.u.], 右轴归一化密度 (项 8)。
    vs='zeta': 横轴 ζ = z - <z> (共动, 项 9)。

    注: 波形为解析模型演示 (非 ASTRA 内部数值), 幅度归一化到峰值。
    """
    from ..constants import EPS0, M_E_KG
    d = np.loadtxt(path, ndmin=2)
    z, n = d[:, 0], d[:, 1]
    if peak_density_cm3:
        n_peak = peak_density_cm3 * 1e6           # cm^-3 -> m^-3
        kp = float(np.sqrt(n_peak * E_CHARGE ** 2
                           / (EPS0 * M_E_KG * C_LIGHT ** 2)))
    else:
        kp = 1.0                                  # 任意单位演示
    zeta = z - float(np.mean(z))
    # 轴上纵向尾场 (手册附录线性等离子体尾场, 幅度峰值归一)
    ez = (kp ** 2 * M_E_KG * C_LIGHT ** 2 * a0 ** 2 / (2 * E_CHARGE)
          * np.sqrt(np.pi / 2) * sigma_z_m
          * np.exp(-kp ** 2 * sigma_z_m ** 2 / 2) * np.cos(kp * zeta))
    amp = float(np.max(np.abs(ez)))
    if amp > 0:
        ez = ez / amp

    fig, ax = plt.subplots(figsize=figsize)
    if vs == "zeta":
        x, xl = zeta * 1e6, "zeta = z - <z> [um]"
    else:
        x, xl = z * 1e3, "z [mm]"
    ax.plot(x, ez, color="C0", label="Ez (on axis, linear wake model)")
    ax.set_xlabel(xl)
    ax.set_ylabel("Ez [a.u.]")
    ax2 = ax.twinx()
    nmax = float(np.max(np.abs(n)))
    if nmax == 0:
        nmax = 1.0
    ax2.plot(x, n / nmax, color="C1", ls="--", label="density (normalized)")
    ax2.set_ylabel("plasma density [arb. u.]", color="C1")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=9, loc="best")
    ax.set_title(title or ("plasma fields vs " + vs))
    fig.tight_layout()
    return fig


def plot_core_fraction_curves(*args, **kwargs):
    """已弃用别名 (批 3): 请用 plot_central_charge_fraction_curves。"""
    import warnings
    warnings.warn("plot_core_fraction_curves 已更名为 "
                  "plot_central_charge_fraction_curves", DeprecationWarning, stacklevel=2)
    return plot_central_charge_fraction_curves(*args, **kwargs)


def plot_central_charge_fraction_curves(dist, fractions=(0.1, 0.25, 0.5, 0.75, 0.9, 1.0),
                              figsize=(10, 4), title=None):
    """束长/发射度 vs 纵向中心电荷分数 (批 3 更名; 非 ASTRA Cemit 口径)."""
    from ..analysis.core import compute_central_charge_fraction_curves
    c = compute_central_charge_fraction_curves(dist, fractions=fractions)
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


def plot_core_emittance_curve(dist, figsize=(9, 4), title=None,
                              bz_on_axis_T: float = 0.0):
    """核心发射度 vs 粒子百分比 (手册 4.13.5, postpro 5.6.1 项 6).

    按单粒子发射度不变量排序取核心, 画 x/y/z 三平面归一化发射度随
    粒子百分比的变化; 100% 处 = 标准 rms 发射度。
    显示单位: 横向 [pi mm mrad] (数值 x1e6), 纵向 [keV mm] (数值 x1)。
    """
    from ..analysis.core_emit import compute_core_emittance_curves
    curves = compute_core_emittance_curves(dist, bz_on_axis_T=bz_on_axis_T)
    fig, ax = _ax(None, figsize)
    specs = (("x", r"$\varepsilon_{nx}$", r"$\pi$ mm mrad"),
             ("y", r"$\varepsilon_{ny}$", r"$\pi$ mm mrad"),
             ("z", r"$\varepsilon_z$", "keV mm"))
    for plane, lab, unit in specs:
        fracs = sorted(curves[plane])
        vals = np.array([curves[plane][f] for f in fracs])
        scale = 1.0 if plane == "z" else 1e6
        ax.plot(np.asarray(fracs) * 100, vals * scale, "o-", label=lab)
    ax.set_xlabel("fraction of particles [%]")
    ax.set_ylabel("core emittance [%s]" % unit)
    ax.set_title(title or "core emittance vs particle fraction (4.13.5)")
    ax.legend(fontsize=9)
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


def plot_cr_emit(cr, ax=None, figsize=(12, 5), title=None):
    """交叉粒子 (Cr_emit, 菜单 4 项 12-14).

    左: 发射度 eps_x/y [pi mm mrad] + 剩余/交叉电荷 [nC];
    右: 排除 cross-over 的束斑 x/y rms [mm] (项 13, 需 x_rms/y_rms 字段)。
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize, squeeze=False)
    ax = axes[0, 0]
    ax.plot(cr["z"], cr["eps_x"] * 1e6, label=r"$\varepsilon_x$")
    ax.plot(cr["z"], cr["eps_y"] * 1e6, label=r"$\varepsilon_y$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel(r"emittance [$\pi$ mm mrad]")
    ax.set_title("emittance (w/o cross-over particles)")
    ax2 = ax.twinx()
    ax2.plot(cr["z"], cr["q_rest"], color="C2", ls="--", label="rest charge")
    ax2.plot(cr["z"], cr["q_cross"], color="C3", ls=":", label="cross charge")
    ax2.set_ylabel("charge [nC]")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=9)

    axb = axes[0, 1]
    if "x_rms" in cr and "y_rms" in cr:
        axb.plot(cr["z"], cr["x_rms"] * 1e3, label=r"$\sigma_x$")
        axb.plot(cr["z"], cr["y_rms"] * 1e3, label=r"$\sigma_y$")
        axb.set_xlabel("z [m]")
        axb.set_ylabel("beam size [mm]")
        axb.set_title("beam size (w/o cross-over particles)")
        axb.legend(fontsize=9)
    else:
        axb.text(0.5, 0.5, "no x_rms/y_rms data", ha="center",
                 va="center", transform=axb.transAxes, color="0.4")
        axb.set_xlabel("z [m]")
        axb.set_title("beam size (w/o cross-over particles)")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig
