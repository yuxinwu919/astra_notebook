"""代码审查修复的回归测试 (2026-08 审查轮)."""
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.distribution import Distribution
from astra_tools.analysis.statistics import compute_statistics
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.analysis.bff import compute_bff
from astra_tools.analysis.emittance import compute_emittance_ellipse_params


def test_equi_charge_duplicate_z_no_fake_current():
    """重复 z 值时 equi_charge 分箱不得产生 ~1e12 A 假电流。"""
    n = 60
    z = np.concatenate([np.full(20, 1e-3), np.linspace(0, 0.9e-3, 40)])
    d = Distribution.from_arrays(
        x=np.zeros(n), y=np.zeros(n), z=z,
        px=np.zeros(n), py=np.zeros(n), pz=np.full(n, 5e6),
        clock=np.zeros(n), charge=np.full(n, 1e-3))   # 0.1 pC/粒子
    sa = compute_slice_analysis(d, n_slices=10, binning="equi_charge")
    assert np.all(np.diff(sa.z_edges) > 0)          # 严格递增
    assert np.all(np.isfinite(sa.current))
    # delta 电荷用箱宽正则化后电流有界 (~1e3 A 量级), 不再是 1e12 A
    assert float(np.max(np.abs(sa.current))) < 1e5


def test_zero_momentum_raises():
    """零动量束团: 发散角/发射度无定义, 应明确报错而非 NaN。"""
    n = 50
    d = Distribution.from_arrays(
        x=np.zeros(n), y=np.zeros(n), z=np.zeros(n),
        px=np.zeros(n), py=np.zeros(n), pz=np.zeros(n),
        clock=np.zeros(n), charge=np.ones(n))
    with pytest.raises(ValueError, match="zero/negative"):
        compute_statistics(d)
    with pytest.raises(ValueError, match="zero/negative"):
        compute_slice_analysis(d, n_slices=4)


def test_ellipse_theta_is_major_axis():
    """1-RMS 椭圆主轴角: 正相关数据 -> +45 度 (不是 -45 短轴)。"""
    rng = np.random.default_rng(0)
    base = rng.normal(size=(2, 50000))
    rho = 0.6
    L = np.array([[1.0, rho], [rho, 1.0]])
    u, up = L @ base
    par = compute_emittance_ellipse_params(u, up)
    assert np.degrees(par["theta"]) == pytest.approx(45.0, abs=2.0)
    # 与协方差矩阵主特征向量一致
    S = np.cov(u, up)
    w, v = np.linalg.eigh(S)
    main = np.arctan2(v[1, np.argmax(w)], v[0, np.argmax(w)])
    assert par["theta"] == pytest.approx(main, abs=0.02) or \
        abs(abs(par["theta"] - main) - np.pi) < 0.02


def test_bff_near_neutral_bunch_returns_zero():
    """近中性束团 (|Σq| << Σ|q|) 的 BFF 归一化发散, 返回零而非爆炸。"""
    rng = np.random.default_rng(5)
    z = rng.normal(0, 1e-3, 3000)
    q = np.where(np.arange(3000) % 2 == 0, 1.0, -1.0)
    b = compute_bff(z, q, kmin=1, kmax=1e4, nk=50)
    assert np.all(b.bff == 0.0)


def test_slice_emittance_unweighted_matches_statistics_convention():
    """slice 矩为群体矩: 与 compute_statistics 约定一致 (均匀电荷)。"""
    rng = np.random.default_rng(42)
    n = 4000
    d = Distribution.from_arrays(
        x=rng.normal(0, 2e-4, n), y=rng.normal(0, 2e-4, n),
        z=rng.uniform(-0.5e-3, 0.5e-3, n),
        px=rng.normal(0, 100.0, n), py=rng.normal(0, 100.0, n),
        pz=np.full(n, 5e6), clock=np.zeros(n), charge=np.ones(n))
    sa = compute_slice_analysis(d, n_slices=20, ref_momentum_eVc=5e6)
    st = compute_statistics(d)
    # 均匀束团: 各 slice 发射度均值应接近全束团值
    mask = sa.n_particles >= 3
    assert np.median(sa.emit_x_norm[mask]) == pytest.approx(
        st.emit_x_norm, rel=0.10)
# ---------------------------------------------------------------
# 审查轮 2: io / namelist / forms / exec 修复
# ---------------------------------------------------------------
def test_binary_9_10_col_ambiguity_uses_index_column(tmp_path):
    """10 列二进制文件在 N%9==0 时必须按 index 列 (第 8 列) 判 10 列."""
    from astra_tools.io.astra_dist import AstraDistributionReader
    rng = np.random.default_rng(0)
    n = 9  # 9 % 9 == 0, 触发歧义分支
    rows = np.column_stack([
        rng.normal(0, 1e-4, n), rng.normal(0, 1e-4, n), np.zeros(n),
        rng.normal(0, 50.0, n), rng.normal(0, 50.0, n),
        np.full(n, 1e6), np.zeros(n), np.full(n, 1e-3),
        np.arange(1, n + 1, dtype=float),   # index 列 (第 8 列, 顺序)
        np.full(n, 5.0),                    # status 列 (恒值)
    ])
    header = np.array([0.0, 1e6, n * 1e-3, 0.0, 0.0])
    p = tmp_path / "test.ini"
    np.concatenate([header, rows.flatten()]).astype(np.float64).tofile(p)
    d = AstraDistributionReader().read(p)
    assert d.n_particle == n
    assert d.index is not None
    assert list(d.index) == list(range(1, n + 1))


def test_track_file_8_columns_keeps_ey(tmp_path):
    """track 文件 8 列: Ez, (Er 或 Ex), (0 或 Ey) 全部保留."""
    from astra_tools.io.astra_misc import read_track_file
    rows = np.column_stack([
        np.arange(1, 11), np.full(10, 5), np.linspace(0, 1.5, 10),
        np.zeros(10), np.zeros(10), np.linspace(0, 1e6, 10),
        np.full(10, 2e5), np.linspace(3e5, 4e5, 10),
    ])
    p = tmp_path / "test.track.001"
    np.savetxt(p, rows)
    d = read_track_file(p)
    assert "Ey" in d and np.allclose(d["Ey"], rows[:, 7])
    assert np.allclose(d["Er"], rows[:, 6])


def test_namelist_quotes_round_trip(tmp_path):
    """引号内 ! / 逗号 / 撇号转义 的写-读往返对称."""
    from astra_tools.namelist.write import write_input_deck
    from astra_tools.namelist.parse import parse_namelists
    deck = {
        "NEWRUN": {
            "Head": "test ! bang , comma ' quote",
            "RUN": 1,
            "Distribution": "bunch.ini",
        },
        "CAVITY": {"MaxE": [-40.0, -20.0], "Nue": [2.857, 1.3]},
    }
    p = tmp_path / "rt.in"
    write_input_deck(deck, p)
    d = parse_namelists(p.read_text())
    assert d["NEWRUN"]["Head"] == "test ! bang , comma ' quote"
    assert d["NEWRUN"]["RUN"] == 1
    assert d["NEWRUN"]["Distribution"] == "bunch.ini"
    assert d["CAVITY"]["MaxE"] == [-40.0, -20.0]


def test_namelist_uppercase_true_false():
    from astra_tools.namelist.parse import parse_namelists
    d = parse_namelists("&NEWRUN\nA=TRUE, B=FALSE, C=.TRUE.\n/\n")["NEWRUN"]
    assert d["A"] is True and d["B"] is False and d["C"] is True


def test_forms_array_params_are_numbers():
    """表单数组参数 (如 Nue/MaxE) 输出数值列表, 不是带引号字符串."""
    from astra_tools.widgets import forms as F
    wmap, getter2 = F.namelist_form("CAVITY", only=["Nue"], show=False)
    wmap["Nue"].value = "2.857, 1.3"
    out2 = getter2()
    assert out2["Nue"] == [2.857, 1.3]


def test_discover_outputs_yemit_zemit_landf(tmp_path):
    from astra_tools.run.exec import discover_outputs
    for name in ("t.Xemit.001", "t.Yemit.001", "t.Zemit.001",
                 "t.Sigma.001", "t.ref.001", "t.Log.001", "t.LandF.001",
                 "t.Cemit.001", "t.0100.001"):
        (tmp_path / name).write_text("0 0 0 0 0 0 0\n")
    out = discover_outputs(tmp_path, "t", run="001")
    assert out["yemit"].name == "t.Yemit.001"
    assert out["zemit"].name == "t.Zemit.001"
    assert out["landf"].name == "t.LandF.001"
    assert [f.name for f in out["phase"]] == ["t.0100.001"]

