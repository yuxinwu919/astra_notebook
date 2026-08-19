"""R5 真实 ErrorS 输出 golden 验证 (Task 1 第三阶段).

金样: examples/Manual_Example/golden/Example.Error.001 — 本地 ASTRA
V4.0 真跑产生 (deck = Manual_Example + &ERROR Loop=T LError=T
ErrorS=T Err_MaxB(1)=0.02T, FOM(1..3)='mean beam energy' /
'horizontal rms emittance' / 'rms bunch length', NEWRUN NLoop=5)。

文件列语义 (手册 Table 4): run#, z, FOM(1..10), 1P,12E17.8。
FOM 数值单位与对应输出文件一致 (MeV / mm mrad / mm)。
read_error (astra_tools/io/astra_misc.py) 此前仅由合成数据测试,
这里用真实 ASTRA 输出验证列数/列语义, 并抽查与名义值运行
(Xemit/Zemit 末行) 的一致性: 误差只加在螺线管 MaxB 上,
能量/束长不受影响, 发射度随 run 小量散布并包围名义值。
"""

from pathlib import Path

import numpy as np
import pytest

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io.astra_misc import read_error

DATA = PROJECT_ROOT / "examples" / "Manual_Example"
ERROR_GOLDEN = DATA / "golden" / "Example.Error.001"


def test_read_error_shape_and_run_column():
    """列数/行数与手册 Table 4 一致: run#, z, FOM(1..10)."""
    err = read_error(ERROR_GOLDEN)
    assert err["run"].shape == (5,)                    # deck: NLoop=5
    assert err["z"].shape == (5,)
    assert err["FOM"].shape == (5, 10)                 # FOM(1..10)
    assert np.array_equal(err["run"], [1, 2, 3, 4, 5])


def test_z_column_is_run_end_position():
    """z 列 = run 末端束团位置 (ZSTOP=1.5 m)."""
    err = read_error(ERROR_GOLDEN)
    assert np.allclose(err["z"], 1.5, atol=1e-6)


def test_fom_columns_semantics_match_manual_table4():
    """FOM 语义与手册一致: 能量[MeV], 横向 rms 发射度[mm mrad],
    束长[mm]; 未用 FOM(4..10) 为 0."""
    err = read_error(ERROR_GOLDEN)
    # FOM(1) = mean beam energy ≈ 999.99 MeV (Zemit 末行 <E_kin>)
    assert err["FOM"][:, 0] == pytest.approx(9.999891e2, rel=1e-6)
    # FOM(2) = horizontal rms emittance: 5 个误差 run 在名义值附近散布
    assert err["FOM"][:, 1] == pytest.approx(1.0003, rel=5e-3)
    # FOM(3) = rms bunch length (螺线管误差不影响, 逐 run 恒定)
    assert err["FOM"][:, 2] == pytest.approx(0.6595227, rel=1e-5)
    # 未指定的 FOM(4..10) 写 0
    assert np.all(err["FOM"][:, 3:] == 0.0)


def test_nominal_consistency_with_xemit_zemit():
    """与名义值运行 (Xemit/Zemit 末行) 的一致性抽查:
    误差 run 的 FOM 均值应接近名义值, 且名义值在散布范围内."""
    xemit = np.loadtxt(DATA / "Example.Xemit.001")
    zemit = np.loadtxt(DATA / "Example.Zemit.001")
    err = read_error(ERROR_GOLDEN)

    e_nominal = zemit[-1, 2]           # <E_kin> = 999.99 MeV
    eps_nominal = xemit[-1, 5]         # eps_n = 1.0003 mm mrad
    zrms_nominal = zemit[-1, 3]        # sig_z = 0.6595 mm

    assert np.mean(err["FOM"][:, 0]) == pytest.approx(e_nominal, rel=1e-4)
    assert np.mean(err["FOM"][:, 1]) == pytest.approx(eps_nominal, rel=5e-3)
    assert np.mean(err["FOM"][:, 2]) == pytest.approx(zrms_nominal, rel=1e-3)
    # 散布包含名义值 (MaxB 误差 ±0.02T, 发射度响应 ±0.03%)
    assert np.min(err["FOM"][:, 1]) <= eps_nominal <= np.max(err["FOM"][:, 1])


def test_frozen_golden_values():
    """归档金样数值冻结 (回归钉住, 任何解析改动立即失败)."""
    err = read_error(ERROR_GOLDEN)
    assert err["FOM"][:, 0] == pytest.approx(
        [9.999891e2] * 5, rel=1e-9)
    assert err["FOM"][:, 1] == pytest.approx(
        [1.000259e0, 1.000645e0, 1.000134e0, 1.000318e0, 1.000665e0],
        rel=1e-9)
    assert err["FOM"][:, 2] == pytest.approx([6.595227e-1] * 5, rel=1e-9)
