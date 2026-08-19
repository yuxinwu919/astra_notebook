"""TWS 周期展开 expand_tws_field_map 回归修复 (task 4)。

背景: 2026-08 第一轮审计把 expand_tws_field_map 当死代码删除 (当时
无调用方, 且 n/m 语义有误), 但 examples/Cavity_Example.ipynb cell 6
一直用它做 TWS_Sband.dat 的 9 单元展开演示绘图 (调用:
expand_tws_field_map(data[:,0], data[:,1], attrs["z1"], attrs["z2"],
m_cells_in_body=attrs["m"], n_cell=9)), e2e 19 本 notebook 中
Cavity_Example 因此 FAIL (ImportError)。

本文件按手册 6.9 语义重新锁定行为 (与 2026-08 审计结论一致: 仅轴上
场演示, 一阶横向展开不在本函数范围, C_numb·n/m 需为整数):
  * 布局 |Entrance| Cells×(n_cell/m) | Exit|; 入口/出口段原表采样
    原样保留 (不动), 周期段 (z1..z2, 含 m 个单元) 按单元平移拼接,
    场值用 np.interp 从原表取 (表格必须覆盖一个完整周期段);
  * n_cell % m != 0 -> ValueError (不再静默截断);
  * 周期段 f(z+λ)=f(z) 逐点误差 0.0, λ = z2-z1 = 104.926 mm
    (与 2026-08 审计验证一致)。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import astra_tools.io as io
from astra_tools.io import field_map as fm
from astra_tools.io.field_map import parse_field_map_file

CAV = PROJECT_ROOT / "examples" / "Cavity_Example"
TWS_FILE = CAV / "TWS_Sband.dat"

# 与 2026-08 审计验证一致: TWS_Sband.dat 头 (z1 z2 n m) = 0.052464
# 0.15739 1.0 3.0, 周期段 1λ = 104.926 mm, C_numb = 9。


def _load_tws():
    """按 notebook 方式读取真实 TWS_Sband.dat (头 4 值 + 401 行表)。"""
    attrs, data = parse_field_map_file(TWS_FILE)
    assert data.shape[0] == 401 and data.shape[1] == 2
    return attrs, data


# ---------------- 导入面回归 (e2e ImportError 的直接原因) ----------------

def test_export_surface_restored():
    """notebook cell 6 的导入路径与 io.__all__ 导出均恢复。"""
    # 与 Cavity_Example.ipynb cell 6 完全相同的导入方式 (导入即验证)
    from astra_tools.io.field_map import read_cavity_field, parse_field_map_file, expand_tws_field_map
    assert callable(read_cavity_field)
    assert callable(parse_field_map_file)
    assert callable(expand_tws_field_map)
    assert "expand_tws_field_map" in io.__all__
    from astra_tools.io import expand_tws_field_map as f1
    assert f1 is expand_tws_field_map


# ---------------- 真实 TWS_Sband.dat (notebook 调用口径) ----------------

def test_real_tws_sband_notebook_expansion():
    """按 notebook 调用展开: 长度增大 / z 单调 / 周期段逐点 0.0 误差。"""
    attrs, data = _load_tws()
    z1, z2, m = attrs["z1"], attrs["z2"], attrs["m"]
    assert (z1, z2, attrs["n"], m) == (0.052464, 0.15739, 1, 3)
    lam = z2 - z1
    assert lam == pytest.approx(104.926e-3)   # 与 2026-08 审计验证一致

    zf, ff = fm.expand_tws_field_map(data[:, 0], data[:, 1], z1, z2,
                                     m_cells_in_body=m, n_cell=9)
    assert len(zf) == len(ff) and len(zf) > len(data)
    assert np.all(np.diff(zf) > 0)            # z 严格单调

    # 入口/出口段与原表一致 (原样保留, 仅周期段重复)
    n_ent = int(np.count_nonzero(data[:, 0] < z1))
    n_exit = int(np.count_nonzero(data[:, 0] >= z2))
    assert n_ent == 100 and n_exit == 101     # 表结构: 入口100 + 周期201 + 出口100
    assert np.array_equal(zf[:n_ent], data[:n_ent, 0])
    assert np.array_equal(ff[:n_ent], data[:n_ent, 1])
    assert np.array_equal(zf[-n_exit:], data[-n_exit:, 0] + (9 // m - 1) * lam)
    assert np.array_equal(ff[-n_exit:], data[-n_exit:, 1])

    # 周期段 f(z+λ) = f(z) 逐点误差 0.0 (n_cell=9, m=3 -> 3 次重复)
    mask0 = (zf >= z1) & (zf < z2)
    zq = zf[mask0] + lam
    idx = np.searchsorted(zf, zq)
    assert np.all(zf[idx] == zq)              # λ 平移后的点精确存在于表内
    assert np.all(ff[idx] == ff[mask0])       # 场值逐点相等
    assert np.max(np.abs(ff[idx] - ff[mask0])) == 0.0


def test_real_tws_repeat_count_layout():
    """n_cell=9, m=3 -> 周期段重复 3 次; 出口平移 (n_repeat-1)λ。"""
    attrs, data = _load_tws()
    z1, z2, m = attrs["z1"], attrs["z2"], attrs["m"]
    lam = z2 - z1
    n_repeat = 9 // m                          # = 3
    n_ent = int(np.count_nonzero(data[:, 0] < z1))
    n_exit = int(np.count_nonzero(data[:, 0] >= z2))
    n_body = len(data) - n_ent - n_exit

    zf, _ = fm.expand_tws_field_map(data[:, 0], data[:, 1], z1, z2,
                                    m_cells_in_body=m, n_cell=9)
    assert len(zf) == n_ent + n_repeat * n_body + n_exit
    assert zf[n_ent] == z1                     # 第一段周期起点 = z1
    assert zf[-n_exit] == z2 + (n_repeat - 1) * lam
    assert zf[-1] == data[-1, 0] + (n_repeat - 1) * lam


@pytest.mark.parametrize("n_cell,m", [(10, 3), (9, 5), (7, 3), (1, 3)])
def test_not_divisible_raises_value_error(n_cell, m):
    """n_cell % m != 0 -> ValueError, 不再静默截断。"""
    attrs, data = _load_tws()
    with pytest.raises(ValueError, match="n_cell="):
        fm.expand_tws_field_map(data[:, 0], data[:, 1], attrs["z1"],
                                attrs["z2"], m_cells_in_body=m, n_cell=n_cell)


# ---------------- 合成表 (解析式场 sin) 拼接正确性 ----------------

def test_synthetic_sin_stitching():
    """均匀合成表: 入口/出口原样, 周期段平移拼接, 场值符合解析式。"""
    z0 = np.linspace(0.0, 1.0, 101)
    z1, z2 = z0[40], z0[60]                    # 取表内精确采样点
    lam = z2 - z1
    f0 = np.sin(2 * np.pi * z0 / lam)          # 周期 λ 的解析场
    n_ent, n_exit, n_body = 40, 41, 20
    n_repeat = 2

    zf, ff = fm.expand_tws_field_map(z0, f0, z1, z2,
                                     m_cells_in_body=2, n_cell=4)
    assert len(zf) == n_ent + n_repeat * n_body + n_exit
    assert np.all(np.diff(zf) > 0)

    # 入口/出口段与原表位对位一致 (出口整体平移 (n_repeat-1)λ)
    assert np.array_equal(zf[:n_ent], z0[:n_ent])
    assert np.array_equal(ff[:n_ent], f0[:n_ent])
    assert np.array_equal(zf[-n_exit:], z0[-n_exit:] + lam)
    assert np.array_equal(ff[-n_exit:], f0[-n_exit:])

    # 周期段: n_repeat 次重复逐位一致, 行间 z 差恰为 λ
    seg = zf[n_ent:n_ent + n_repeat * n_body].reshape(n_repeat, n_body)
    fsg = ff[n_ent:n_ent + n_repeat * n_body].reshape(n_repeat, n_body)
    assert np.array_equal(fsg[1], fsg[0])
    assert np.array_equal(seg[1], seg[0] + lam)

    # 解析一致性: 周期段场值 ≈ sin(2π z/λ)
    zper = zf[n_ent:n_ent + n_repeat * n_body]
    fper = ff[n_ent:n_ent + n_repeat * n_body]
    assert np.allclose(fper, np.sin(2 * np.pi * zper / lam), atol=1e-10)


# ---------------- 输入护栏 ----------------

def test_invalid_inputs_raise():
    z0 = np.linspace(0.0, 1.0, 11)
    f0 = np.sin(z0)
    with pytest.raises(ValueError):
        fm.expand_tws_field_map(z0, f0, 0.6, 0.4, 1, 2)        # z2 <= z1
    with pytest.raises(ValueError):
        fm.expand_tws_field_map(z0[::-1], f0, 0.4, 0.6, 1, 2)  # 非递增 z
    with pytest.raises(ValueError):
        fm.expand_tws_field_map(z0, f0, 0.4, 1.2, 1, 2)        # 周期段越界
    with pytest.raises(ValueError):
        fm.expand_tws_field_map(z0, f0, 0.4, 0.6, 0, 2)        # m < 1
