"""其余 ASTRA 输出文件读取器 (手册 Table 4 全量覆盖).

单位换算与 Xemit 同约定: 手册 'pi' 单位数值上即 SI 前缀单位
(pi mrad mm == mm.mrad, pi um == um), 见 docs/physics_notes/01。

文件类型 (手册 Table 3/4):
  track   探针轨迹:  seq stat z[m] x[mm] y[mm] Ez[V/m] Er[V/m]
  Cathode 阴极发射:  z t[ns] sp-ch场 acc场 q[nC] 网格边界 发射标志
  Xemit2  缩减发射度: z K2z K3z eps_red_z K2E K3E eps_red_zE
  TRemit  trace发射度: z t eps_tr_x eps_tr_y eps_tr_z
  Cr_emit 交叉粒子:   z t x_rms y_rms eps_x eps_y 剩余电荷 交叉电荷
  Larmor  拉莫尔角:   z t avr rms
  PScan   相位扫描:   phase[deg] E_kin[MeV] 压缩因子 beta/beta0
  Scan    参数扫描:   para z FOM(1..10)
  Error   误差扫描:   run# z FOM(1..10)
  Density 粒子密度:   z t N(i) dens(i) i=1..5
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

MM_TO_M = 1e-3
NS_TO_S = 1e-9
MEV_TO_EV = 1e6


def _load(path, ncols):
    d = np.loadtxt(path, ndmin=2)
    if d.shape[1] < ncols:
        raise ValueError("%s: 需要 >= %d 列, 实际 %d" % (path, ncols, d.shape[1]))
    return d


def read_track_file(path) -> dict:
    """探针轨迹 (track 文件): 每行 = 一个探针在一个积分步的状态."""
    d = _load(path, 7)
    return {
        "seq": d[:, 0].astype(int), "status": d[:, 1].astype(int),
        "z": d[:, 2], "x": d[:, 3] * MM_TO_M, "y": d[:, 4] * MM_TO_M,
        "Ez": d[:, 5], "Er": d[:, 6],
    }


def read_cathode_file(path) -> dict:
    """阴极发射过程 (Cathode 文件)."""
    d = _load(path, 8)
    return {
        "z": d[:, 0], "t": d[:, 1] * NS_TO_S,
        "E_spch": d[:, 2], "E_acc": d[:, 3], "q": d[:, 4],
        "grid_min": d[:, 5], "grid_max": d[:, 6], "flag": d[:, 7],
    }


def read_xemit2(path) -> dict:
    """缩减发射度 (Xemit2/Yemit2, Lsub_cor=T 时生成)."""
    d = _load(path, 7)
    return {
        "z": d[:, 0],
        "K2z": d[:, 1], "K3z": d[:, 2],          # pi rad m (数值 rad.m)
        "eps_red_z": d[:, 3] * 1e-6,               # pi mrad mm -> m.rad
        "K2E": d[:, 4], "K3E": d[:, 5],
        "eps_red_zE": d[:, 6] * 1e-6,
    }


def read_tremit(path) -> dict:
    """trace-space 发射度 (TRemit, TR_EmitS=T 时生成)."""
    d = _load(path, 5)
    return {
        "z": d[:, 0], "t": d[:, 1] * NS_TO_S,
        "eps_tr_x": d[:, 2] * 1e-6,                # pi mrad mm -> m.rad
        "eps_tr_y": d[:, 3] * 1e-6,
        "eps_tr_z": d[:, 4] * 1e-6,                # pi um -> m
    }


def read_cr_emit(path) -> dict:
    """交叉粒子发射度 (Cr_emit, Cross_start != Cross_end 时生成)."""
    d = _load(path, 8)
    return {
        "z": d[:, 0], "t": d[:, 1] * NS_TO_S,
        "x_rms": d[:, 2] * MM_TO_M, "y_rms": d[:, 3] * MM_TO_M,
        "eps_x": d[:, 4] * 1e-6, "eps_y": d[:, 5] * 1e-6,
        "q_rest": d[:, 6], "q_cross": d[:, 7],
    }


def read_larmor(path) -> dict:
    """拉莫尔角 (Larmor 文件)."""
    d = _load(path, 4)
    return {"z": d[:, 0], "t": d[:, 1] * NS_TO_S,
            "avr": d[:, 2], "rms": d[:, 3]}


def read_pscan(path) -> dict:
    """相位扫描 (PScan 文件, Phase_Scan=T 时生成)."""
    d = _load(path, 4)
    return {"phase_deg": d[:, 0], "E_kin_eV": d[:, 1] * MEV_TO_EV,
            "compression": d[:, 2], "beta_ratio": d[:, 3]}


def read_scan(path) -> dict:
    """参数扫描 (Scan 文件, LScan=T 时生成): para, z, FOM(1..10)."""
    d = _load(path, 12)
    return {"para": d[:, 0], "z": d[:, 1], "FOM": d[:, 2:]}


def read_error(path) -> dict:
    """误差扫描 (Error 文件): run#, z, FOM(1..10)."""
    d = _load(path, 12)
    return {"run": d[:, 0].astype(int), "z": d[:, 1], "FOM": d[:, 2:]}


def read_density(path) -> dict:
    """粒子密度 (Density 文件, DensityS=T 时生成)."""
    d = _load(path, 12)
    return {"z": d[:, 0], "t": d[:, 1] * NS_TO_S,
            "N": d[:, 2:7], "dens": d[:, 7:12]}
