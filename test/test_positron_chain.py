"""第三阶段 Task 3: 正电子全链口径验证 (纯内部验证, 不跑 ASTRA).

五条消费者路径口径一致 (2026-08 F4 种类感知正则动量符号 canonical_signs):
  1. statistics.compute_statistics
  2. slices.compute_slice_analysis
  3. core_emit._plane_coords / single_particle_amplitudes / 核心发射度
  4. plot.phase_space.plot_phase_space (x-x', y-y' 绘图数据)
  5. plot.arbitrary_phase_space.param_columns (xp/yp 列)

物理口径 (emittance.py 4.13.1):
    p~x = px + s*(c/2)*Bz*y,  p~y = py - s*(c/2)*Bz*x
    电子 (species=1) s=+1; 正电荷种类 (2=正电子) s=-1。
  正则动量 P = p + qA (q 带符号), A_x = -Bz*y/2, A_y = +Bz*x/2:
    - 电子零正则动量束:  px = -(c/2)Bz*y, py = +(c/2)Bz*x
    - 正电子零正则动量束: px = +(c/2)Bz*y, py = -(c/2)Bz*x
  任务简报字面给出的 "px = -c Bz y/2, py = +c Bz x/2" 是**电子**的零
  正则动量构造 (brief 括注判据 p~x = px - c Bz y/2 = 0 推出正电子应为
  px = +c Bz y/2 — 以物理判据为准, 两种构造均覆盖):
    - 正电子 + 零正则构造 (px=+cBz y/2): 五条路径 sig_xp ≈ 0 仅当符号正确
      (若某路径仍硬编码电子符号 s=+1, 会得到 ±c*Bz*y 的非零散角);
    - 正电子 + 电子式构造 (px=-cBz y/2): 正确实现给出 sig_xp ≈ c*Bz*σy/p_ref
      (符号误用的反向探测器, 与 test_slice_species_eigen_fixes 一致);
    - 电子 + 零正则构造: 全部路径行为不回归。

关键数值 (Bz=0.5 T, σx=σy=1e-3 m, p_ref=1e9 eV/c):
    符号误用判别尺度 c*Bz*σy/p_ref ≈ 1.49e-4 rad (0.149 mrad);
    断言容差 1e-8 rad, 低于判别尺度 4 个数量级, 高于浮点噪声。
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.constants import C_LIGHT
from astra_tools.distribution import Distribution
from astra_tools.analysis.statistics import compute_statistics
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.analysis.core_emit import (
    _plane_coords,
    single_particle_amplitudes,
    compute_core_emittance_by_fraction,
)
from astra_tools.plot.phase_space import plot_phase_space
from astra_tools.plot.arbitrary_phase_space import param_columns, plot_arbitrary
import matplotlib.pyplot as plt

BZ = 0.5          # 螺线管轴上场 [T]
P_REF = 1.0e9     # 参考动量 [eV/c] = 1 GeV/c
SIG_XY = 1.0e-3   # 横向均方根尺寸 [m]
SIG_UPC = 1.0e5   # 非零正则散角动量的均方根 [eV/c]
# 符号误用判别尺度: c*Bz*σy/p_ref ≈ 1.5e-4 rad
# 零正则断言容差: 低于判别尺度 4 个数量级, 高于浮点噪声
TOL_RAD = 1e-8
TOL_MRAD = 1e-6


def _canonical_species_sign(index: np.ndarray) -> np.ndarray:
    """每粒子正则动量符号 (与 canonical_signs 口径一致): 正电荷种类 -1."""
    return np.where(index == 2, -1.0, 1.0)


def _beam(species, n=10000, seed=9, px_can=0.0, py_can=0.0,
          index=None, z_uniform=True, charge_sign=None):
    """解析测试束构造 (正/电子零正则动量或指定正则相空间).

    px_can/py_can: 每粒子的正则动量散角部分 [eV/c] (标量或数组)。
    物理动量: px = px_can - s*(c/2)*Bz*y, py = py_can + s*(c/2)*Bz*x,
    使代码口径 p~ = px + s*(c/2)*Bz*y 精确回到 px_can。
    """
    rng = np.random.default_rng(seed)
    n = int(n)
    x = rng.normal(0.0, SIG_XY, n)
    y = rng.normal(0.0, SIG_XY, n)
    if z_uniform:
        z = rng.uniform(-3e-3, 3e-3, n) + rng.normal(0.0, 1e-6, n)
    else:
        z = rng.normal(0.0, 1e-3, n)
    if index is None:
        index = np.full(n, species, dtype=np.int32)
    index = np.asarray(index, dtype=np.int32)
    s = _canonical_species_sign(index)
    if np.isscalar(px_can):
        px_can = np.full(n, float(px_can))
    if np.isscalar(py_can):
        py_can = np.full(n, float(py_can))
    t = 0.5 * C_LIGHT * BZ * y
    u = 0.5 * C_LIGHT * BZ * x
    px = np.asarray(px_can, dtype=float) - s * t
    py = np.asarray(py_can, dtype=float) + s * u
    if charge_sign is None:
        # 电子负宏电荷, 正电子正宏电荷 (ASTRA Table 1 惯例)
        charge = np.where(index == 2, 2e-3, -2e-3)
    else:
        charge = np.full(n, float(charge_sign))
    return Distribution(
        x=x, y=y, z=z, px=px, py=py, pz=np.full(n, P_REF),
        clock=np.zeros(n), charge=charge,
        status=np.full(n, 5), index=index,
        ref_momentum_eVc=P_REF,
    )


def _zero_can_beam(species, **kw):
    """零正则动量束: 正电子 px=+cBz y/2, 电子 px=-cBz y/2 (见模块 docstring)."""
    return _beam(species, **kw)


# ---------------------------------------------------------------- statistics
def test_statistics_positron_zero_canonical():
    """正电子零正则动量束: sig_xp/sig_yp ≈ 0 (符号正确才可能)."""
    dist = _zero_can_beam(2, n=10000, seed=11)
    st = compute_statistics(dist, bz_on_axis_T=BZ)
    assert st.sig_xp < TOL_RAD, "正电子 sig_xp 应为 0, 得到 %g" % st.sig_xp
    assert st.sig_yp < TOL_RAD, "正电子 sig_yp 应为 0, 得到 %g" % st.sig_yp
    # 零发散束几何发射度恒为零
    assert st.emit_x_geom == pytest.approx(0.0, abs=1e-20)
    assert st.emit_y_geom == pytest.approx(0.0, abs=1e-20)
    # 电荷符号为正 (正电子)
    assert st.total_charge_nC > 0


def test_statistics_electron_control_no_regression():
    """电子对照束 (同构造, 零正则动量): 行为不回归."""
    dist = _zero_can_beam(1, n=10000, seed=11)
    st = compute_statistics(dist, bz_on_axis_T=BZ)
    assert st.sig_xp < TOL_RAD
    assert st.sig_yp < TOL_RAD
    assert st.total_charge_nC < 0


def test_statistics_positron_with_electron_style_momenta_stays_nonzero():
    """正电子用电子式动量 (brief 字面 px=-cBz y/2): 散角必须是非零的
    c*Bz*σy/p_ref — 若实现把正电子当电子处理 (硬编码 s=+1), 此值会错误
    归零. 反向探测器, 与 test_slice_species_eigen_fixes.py 一致."""
    # 显式构造: px = -c Bz y/2 (电子式动量), 正电子代码口径下 ptx = -c Bz y
    rng = np.random.default_rng(13)
    n = 20000
    y = rng.normal(0.0, SIG_XY, n)
    x = rng.normal(0.0, SIG_XY, n)
    z = rng.uniform(-3e-3, 3e-3, n) + rng.normal(0.0, 1e-6, n)
    t = 0.5 * C_LIGHT * BZ * y
    u = 0.5 * C_LIGHT * BZ * x
    dist = Distribution(
        x=x, y=y, z=z, px=-t, py=+u, pz=np.full(n, P_REF),
        clock=np.zeros(n), charge=np.full(n, 2e-3),
        status=np.full(n, 5), index=np.full(n, 2, dtype=np.int32),
        ref_momentum_eVc=P_REF)
    st = compute_statistics(dist, bz_on_axis_T=BZ)
    expected = C_LIGHT * BZ * float(np.std(y)) / P_REF
    assert st.sig_xp == pytest.approx(expected, rel=1e-6)
    assert st.sig_xp > 1e-5   # 远离零 (判别尺度 ~1.5e-4 rad)


def test_statistics_nonzero_canonical_matches_electron():
    """非零正则相空间 (px_can~N(0, 1e5)): 正电子与电子对照束逐项一致,
    且 sig_xp = σ(px_can)/p_ref (p_ref 口径一致性)."""
    rng = np.random.default_rng(17)
    px_can = rng.normal(0.0, SIG_UPC, 10000)
    py_can = rng.normal(0.0, SIG_UPC, 10000)
    # 注意: 正则散角随机流与 _beam 内部几何随机流必须不同种子,
    # 否则 x 与 px_can 完全相关 (同一 RNG 流的同比例缩放) -> 束矩阵退化
    pe = _beam(1, n=10000, seed=17001, px_can=px_can, py_can=py_can)
    pp = _beam(2, n=10000, seed=17001, px_can=px_can, py_can=py_can)
    se = compute_statistics(pe, bz_on_axis_T=BZ)
    sp = compute_statistics(pp, bz_on_axis_T=BZ)
    for key in ("sig_xp", "sig_yp", "emit_x_geom", "emit_y_geom",
                "emit_x_norm", "emit_y_norm", "beta_x", "alpha_x",
                "beta_y", "alpha_y"):
        assert getattr(sp, key) == pytest.approx(getattr(se, key), rel=1e-9), key
    assert sp.sig_xp == pytest.approx(float(np.std(px_can)) / P_REF, rel=1e-6)
    assert sp.sig_yp == pytest.approx(float(np.std(py_can)) / P_REF, rel=1e-6)
    assert sp.sig_xp > 1e-5   # 非退化: 判别尺度之上
    # 非退化束发射度: eps_geom ~ sigma_x*sigma_x' ~ 1e-7, 归一化 ~ 2e-4 m.rad
    assert sp.emit_x_norm > 1e-6 and sp.emit_y_norm > 1e-6


def test_statistics_mixed_beam_per_particle_signs():
    """混合束 (半电子半正电子, 逐粒子零正则动量): 逐粒子符号必须正确,
    全局 sig_xp ≈ 0 — 若用单一全局符号, 一半粒子散角 ±c*Bz*y, sig_xp
    ≈ c*Bz*σy/(√2 p_ref) ≈ 1e-4 rad 远大于容差."""
    n = 10000
    rng = np.random.default_rng(23)
    index = np.concatenate([np.full(n // 2, 1, dtype=np.int32),
                            np.full(n - n // 2, 2, dtype=np.int32)])
    rng.shuffle(index)
    # 逐粒子零正则构造 (两种物种各用各的符号)
    x = rng.normal(0.0, SIG_XY, n)
    y = rng.normal(0.0, SIG_XY, n)
    z = rng.uniform(-3e-3, 3e-3, n) + rng.normal(0.0, 1e-6, n)
    s = _canonical_species_sign(index)
    t = 0.5 * C_LIGHT * BZ * y
    u = 0.5 * C_LIGHT * BZ * x
    dist = Distribution(
        x=x, y=y, z=z, px=-s * t, py=+s * u, pz=np.full(n, P_REF),
        clock=np.zeros(n), charge=np.where(index == 2, 2e-3, -2e-3),
        status=np.full(n, 5), index=index, ref_momentum_eVc=P_REF)
    st = compute_statistics(dist, bz_on_axis_T=BZ)
    assert st.sig_xp < TOL_RAD, "混合束 sig_xp 应为 0, 得到 %g" % st.sig_xp
    assert st.sig_yp < TOL_RAD
    # 加权模式 (|q| 均匀 -> 与无加权一致)
    sw = compute_statistics(dist, bz_on_axis_T=BZ, use_weights=True)
    assert sw.sig_xp < TOL_RAD
    assert sw.sig_xp == pytest.approx(st.sig_xp, abs=1e-20)


# -------------------------------------------------------------------- slices
def _slice_ok(sa, tol=TOL_RAD):
    """全部有粒子 slice 的 sig_xp/sig_yp 均小于容差."""
    nz = sa.n_particles > 0
    assert np.all(nz), "存在空 slice (测试束构造问题)"
    return (np.max(np.abs(sa.sig_xp[nz])) < tol
            and np.max(np.abs(sa.sig_yp[nz])) < tol)


def test_slices_positron_zero_canonical():
    """正电子: 各 slice 正则散角 ≈ 0."""
    dist = _zero_can_beam(2, n=6000, seed=29)
    sa = compute_slice_analysis(dist, n_slices=8, bz_on_axis_T=BZ)
    assert _slice_ok(sa), "slice sig_xp/sig_yp 应为 0"
    assert np.all(sa.sig_xp < TOL_RAD)
    assert np.all(sa.sig_yp < TOL_RAD)


def test_slices_electron_control_no_regression():
    """电子对照束: 各 slice 行为不回归."""
    dist = _zero_can_beam(1, n=6000, seed=29)
    sa = compute_slice_analysis(dist, n_slices=8, bz_on_axis_T=BZ)
    assert _slice_ok(sa)


def test_slices_mixed_beam_per_slice_signs():
    """混合束: 逐 slice 逐粒子符号正确 -> 各 slice sig_xp ≈ 0."""
    n = 6000
    rng = np.random.default_rng(31)
    index = np.concatenate([np.full(n // 2, 1, dtype=np.int32),
                            np.full(n - n // 2, 2, dtype=np.int32)])
    rng.shuffle(index)
    x = rng.normal(0.0, SIG_XY, n)
    y = rng.normal(0.0, SIG_XY, n)
    z = rng.uniform(-3e-3, 3e-3, n) + rng.normal(0.0, 1e-6, n)
    s = _canonical_species_sign(index)
    t = 0.5 * C_LIGHT * BZ * y
    u = 0.5 * C_LIGHT * BZ * x
    dist = Distribution(
        x=x, y=y, z=z, px=-s * t, py=+s * u, pz=np.full(n, P_REF),
        clock=np.zeros(n), charge=np.where(index == 2, 2e-3, -2e-3),
        status=np.full(n, 5), index=index, ref_momentum_eVc=P_REF)
    sa = compute_slice_analysis(dist, n_slices=8, bz_on_axis_T=BZ)
    assert _slice_ok(sa), "混合束 slice sig_xp/sig_yp 应为 0"


# ---------------------------------------------------------------- core_emit
def test_core_emit_zero_canonical_amplitudes_vanish():
    """零正则动量束: 单粒子振幅恒为零 (零发射度), 正/电子一致."""
    pe = _zero_can_beam(1, n=4000, seed=37)
    pp = _zero_can_beam(2, n=4000, seed=37)
    for plane in ("x", "y"):
        je = single_particle_amplitudes(pe, plane=plane, bz_on_axis_T=BZ)
        jp = single_particle_amplitudes(pp, plane=plane, bz_on_axis_T=BZ)
        assert np.all(je == 0.0)
        assert np.all(jp == 0.0)
        ce = compute_core_emittance_by_fraction(
            pe, plane=plane, fractions=(0.8, 0.95, 1.0), bz_on_axis_T=BZ)
        cp = compute_core_emittance_by_fraction(
            pp, plane=plane, fractions=(0.8, 0.95, 1.0), bz_on_axis_T=BZ)
        for f in ce:
            assert ce[f] == 0.0 and cp[f] == 0.0


def test_core_emit_amplitudes_match_electron():
    """非零正则相空间: 单粒子振幅逐粒子与电子对照束一致 (横向)."""
    rng = np.random.default_rng(41)
    px_can = rng.normal(0.0, SIG_UPC, 4000)
    py_can = rng.normal(0.0, SIG_UPC, 4000)
    pe = _beam(1, n=4000, seed=41001, px_can=px_can, py_can=py_can)
    pp = _beam(2, n=4000, seed=41001, px_can=px_can, py_can=py_can)
    for plane in ("x", "y"):
        je = single_particle_amplitudes(pe, plane=plane, bz_on_axis_T=BZ)
        jp = single_particle_amplitudes(pp, plane=plane, bz_on_axis_T=BZ)
        assert np.all(je > 0), "非退化束振幅应 > 0"
        np.testing.assert_allclose(jp, je, rtol=1e-9, atol=0.0)


def test_core_emit_fractions_match_electron_and_statistics():
    """核心发射度曲线与电子对照束一致; f=1.0 与统计模块归一化发射度一致."""
    rng = np.random.default_rng(43)
    px_can = rng.normal(0.0, SIG_UPC, 4000)
    py_can = rng.normal(0.0, SIG_UPC, 4000)
    pe = _beam(1, n=4000, seed=43001, px_can=px_can, py_can=py_can)
    pp = _beam(2, n=4000, seed=43001, px_can=px_can, py_can=py_can)
    se = compute_statistics(pe, bz_on_axis_T=BZ)
    sp = compute_statistics(pp, bz_on_axis_T=BZ)
    for plane, key in (("x", "emit_x_norm"), ("y", "emit_y_norm")):
        ce = compute_core_emittance_by_fraction(
            pe, plane=plane, fractions=(0.8, 0.9, 0.95, 1.0), bz_on_axis_T=BZ)
        cp = compute_core_emittance_by_fraction(
            pp, plane=plane, fractions=(0.8, 0.9, 0.95, 1.0), bz_on_axis_T=BZ)
        for f in ce:
            assert cp[f] == pytest.approx(ce[f], rel=1e-9), "f=%s" % f
        assert cp[1.0] == pytest.approx(getattr(sp, key), rel=1e-9)
        assert ce[1.0] == pytest.approx(getattr(se, key), rel=1e-9)
        # 核心发射度随分数递减 (非退化束)
        vals = [ce[f] for f in (0.8, 0.9, 0.95, 1.0)]
        assert vals[0] < vals[1] < vals[2] < vals[3]


def test_plane_coords_direct():
    """_plane_coords: 零正则 -> up 恒零; 非零正则 -> 与电子逐元素一致."""
    pe0 = _zero_can_beam(1, n=2000, seed=47)
    pp0 = _zero_can_beam(2, n=2000, seed=47)
    for plane in ("x", "y"):
        u, up, w = _plane_coords(pp0, plane, BZ, P_REF, pp0.charge)
        assert np.all(up == 0.0), "零正则束 up 应恒为零"
        ue, upe, we = _plane_coords(pe0, plane, BZ, P_REF, pe0.charge)
        np.testing.assert_allclose(u, ue, rtol=0.0, atol=0.0)
    rng = np.random.default_rng(47)
    px_can = rng.normal(0.0, SIG_UPC, 2000)
    py_can = rng.normal(0.0, SIG_UPC, 2000)
    pe1 = _beam(1, n=2000, seed=47001, px_can=px_can, py_can=py_can)
    pp1 = _beam(2, n=2000, seed=47001, px_can=px_can, py_can=py_can)
    for plane in ("x", "y"):
        up_p, up_e = (_plane_coords(d, plane, BZ, P_REF, d.charge)[1]
                      for d in (pp1, pe1))
        np.testing.assert_allclose(up_p, up_e, rtol=1e-9, atol=0.0)


# ------------------------------------------------------------ plot paths
def _plotted_divergence(fig):
    """散点集合的 y 坐标 (显示单位 mrad)."""
    return np.asarray(fig.axes[0].collections[0].get_offsets())[:, 1]


def test_phase_space_plot_positron_divergence_zero():
    """plot.phase_space 数据路径: 正电子零正则束 x'/y' 散角 ≈ 0."""
    dist = _zero_can_beam(2, n=10000, seed=53)
    for plane in ("x", "y"):
        fig = plot_phase_space(dist, plane=plane, bz_on_axis_T=BZ)
        try:
            yp = _plotted_divergence(fig)
            assert np.max(np.abs(yp)) < TOL_MRAD, \
                "plane %s 绘图散角应 ≈ 0, 得到 max|y'|=%g mrad" % (plane, np.max(np.abs(yp)))
        finally:
            plt.close(fig)


def test_phase_space_plot_electron_control_no_regression():
    """电子对照束绘图路径不回归."""
    dist = _zero_can_beam(1, n=10000, seed=53)
    for plane in ("x", "y"):
        fig = plot_phase_space(dist, plane=plane, bz_on_axis_T=BZ)
        try:
            yp = _plotted_divergence(fig)
            assert np.max(np.abs(yp)) < TOL_MRAD
        finally:
            plt.close(fig)


def test_phase_space_plot_positron_electron_style_nonzero():
    """绘图路径反向探测器: 正电子用电子式动量 -> 绘图散角 ≠ 0
    (~c*Bz*σy/p_ref ≈ 0.15 mrad), 若硬编码电子符号则错误归零."""
    n = 10000
    rng = np.random.default_rng(59)
    y = rng.normal(0.0, SIG_XY, n)
    x = rng.normal(0.0, SIG_XY, n)
    z = rng.uniform(-3e-3, 3e-3, n) + rng.normal(0.0, 1e-6, n)
    t = 0.5 * C_LIGHT * BZ * y
    u = 0.5 * C_LIGHT * BZ * x
    dist = Distribution(
        x=x, y=y, z=z, px=-t, py=+u, pz=np.full(n, P_REF),
        clock=np.zeros(n), charge=np.full(n, 2e-3),
        status=np.full(n, 5), index=np.full(n, 2, dtype=np.int32),
        ref_momentum_eVc=P_REF)
    fig = plot_phase_space(dist, plane="x", bz_on_axis_T=BZ)
    try:
        yp = _plotted_divergence(fig)
        assert float(np.std(yp)) > 0.1, "判别尺度 0.15 mrad 之上"
    finally:
        plt.close(fig)


def test_param_columns_positron_zero():
    """arbitrary_phase_space.param_columns: 正电子 xp/yp ≈ 0,
    且与 statistics.sig_xp 跨路径一致 (p_ref 口径相同)."""
    dist = _zero_can_beam(2, n=10000, seed=61)
    cols = param_columns(dist, bz_on_axis_T=BZ)
    for name in ("xp", "yp"):
        vals = np.asarray(cols[name][0])
        assert np.max(np.abs(vals)) < TOL_MRAD, "%s 列应 ≈ 0" % name
    st = compute_statistics(dist, bz_on_axis_T=BZ)
    # mrad -> rad: 与统计模块 sig_xp 一致 (同一 p_ref 口径)
    assert float(np.std(cols["xp"][0])) / 1e3 == pytest.approx(st.sig_xp, abs=1e-12)
    assert float(np.std(cols["yp"][0])) / 1e3 == pytest.approx(st.sig_yp, abs=1e-12)
    # x/y 位置列不受影响
    assert np.all(cols["x"][0] == dist.x * 1e3)
    assert np.all(cols["y"][0] == dist.y * 1e3)


def test_param_columns_electron_control_no_regression():
    """电子对照束 param_columns 不回归."""
    dist = _zero_can_beam(1, n=10000, seed=61)
    cols = param_columns(dist, bz_on_axis_T=BZ)
    for name in ("xp", "yp"):
        assert np.max(np.abs(cols[name][0])) < TOL_MRAD


def test_param_columns_nonzero_canonical_matches_electron():
    """非零正则: xp/yp 列与电子对照束逐元素一致."""
    rng = np.random.default_rng(67)
    px_can = rng.normal(0.0, SIG_UPC, 10000)
    py_can = rng.normal(0.0, SIG_UPC, 10000)
    pe = _beam(1, n=10000, seed=67001, px_can=px_can, py_can=py_can)
    pp = _beam(2, n=10000, seed=67001, px_can=px_can, py_can=py_can)
    cpe = param_columns(pe, bz_on_axis_T=BZ)
    cpp = param_columns(pp, bz_on_axis_T=BZ)
    for name in ("xp", "yp"):
        np.testing.assert_allclose(cpp[name][0], cpe[name][0], rtol=1e-9, atol=0.0)
    # 与统计模块 sig_xp 一致
    sp = compute_statistics(pp, bz_on_axis_T=BZ)
    assert float(np.std(cpp["xp"][0])) / 1e3 == pytest.approx(sp.sig_xp, rel=1e-9)
    assert float(np.std(cpp["yp"][0])) / 1e3 == pytest.approx(sp.sig_yp, rel=1e-9)


def test_arbitrary_plot_positron_uses_zero_divergence():
    """plot_arbitrary (x-xp) 消费 param_columns 的数据路径: 散角 ≈ 0."""
    dist = _zero_can_beam(2, n=10000, seed=71)
    fig = plot_arbitrary(dist, "x", "xp", bz_on_axis_T=BZ)
    try:
        yp = _plotted_divergence(fig)
        assert np.max(np.abs(yp)) < TOL_MRAD
    finally:
        plt.close(fig)
