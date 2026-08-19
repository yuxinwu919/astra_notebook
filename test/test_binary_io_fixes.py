"""2026-08 对抗性审计修复:二进制/ASCII 分布 IO (astra_dist.py).

每个测试的输入文件都用 numpy 独立构造 (不经本项目 writer),
避免循环验证。对应审计发现:
  P1-1 二进制 9/10 列消歧前提错误 (第 9 列 = 粒子种类, 恒为 1)
  P1-2 写入头带符号总电荷 -> 应为 |Q|
  P2-1 ASCII 写入格式应为参考粒子首行 (10 列绝对坐标)
  P2-3 index=None 时种类列应写 1 (电子) 而非 1..N
  P3   字节序探测 / 小文件二进制误判为 ASCII
"""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution, write_distribution
from astra_tools.distribution import Distribution


def _make_10col_binary(path, header, n=90, species=1, status=5, dtype="<f8"):
    """独立构造 ASTRA 10 列二进制分布文件 (不经项目 writer)."""
    rng = np.random.default_rng(7)
    x = rng.normal(0, 1e-4, n)
    y = rng.normal(0, 1e-4, n)
    z = rng.normal(0, 1e-4, n)          # relative z
    px = rng.normal(0, 1e3, n)
    py = rng.normal(0, 1e3, n)
    pz = rng.normal(0, 1e5, n)          # relative pz [eV/c]
    clock = rng.normal(0, 1e-4, n)      # relative clock [ns]
    charge = np.full(n, -2e-3)
    spec = np.full(n, species)
    stat = np.full(n, status)
    body = np.column_stack([x, y, z, px, py, pz, clock, charge, spec, stat])
    hdr = np.asarray(header, dtype=np.float64)
    with open(path, "wb") as fh:
        hdr.astype(dtype).tofile(fh)
        body.astype(dtype).tofile(fh)
    return x, y, z, px, py, pz, clock


def test_binary_10col_species_all_ones_reads_10_cols(tmp_path):
    """N=90 (9/10 列歧义区): 第 9 列恒 1 (电子种类) 的真实 10 列文件
    必须按 10 列读取, status 列不被种类列顶替 (P1-1)."""
    p = tmp_path / "b10.001"
    header = [5.0035, 1.0005e9, 1.0, -6.8e-13, -3.9e-13]  # Q=+|Q|
    _make_10col_binary(p, header, n=90, species=1, status=5)
    dist = read_distribution(p)
    assert dist.n_particle == 90
    assert dist.n_active == 90                      # 全部 status=5
    assert np.all(dist.status == 5)                 # status 列正确
    assert np.all(dist.index == 1)                  # 种类列正确
    assert dist.ref_momentum_eVc == pytest.approx(1.0005e9)
    assert dist.total_charge_nC == pytest.approx(1.0)
    # z 保持相对 (binary 头无 z), pz/clock 相对 -> 绝对
    assert dist.ref_z_m == 0.0


def test_binary_9col_legacy_still_readable(tmp_path):
    """9 列遗留格式 (非 ASTRA 标准) 回退读取 + 告警 (回归保护)."""
    p = tmp_path / "b9.001"
    rng = np.random.default_rng(11)
    n = 10
    x = rng.normal(0, 1e-4, n)
    y = rng.normal(0, 1e-4, n)
    z = rng.normal(0, 1e-4, n)
    px = rng.normal(0, 1e3, n)
    py = rng.normal(0, 1e3, n)
    pz = rng.normal(0, 1e5, n)
    clock = rng.normal(0, 1e-4, n)
    charge = np.full(n, -2e-3)
    stat = np.full(n, 5)
    body = np.column_stack([x, y, z, px, py, pz, clock, charge, stat])
    hdr = np.array([5.0, 1.0e9, 0.02, 1e-3, -2e-3])
    with open(p, "wb") as fh:
        hdr.tofile(fh)
        body.astype(np.float64).tofile(fh)
    dist = read_distribution(p)
    assert dist.n_particle == n
    assert np.all(dist.status == 5)
    assert dist.n_active == n


def test_write_binary_header_charge_positive(tmp_path):
    """写入器文件头 Q 必须为 |Q| (P1-2): 电子束 (负宏电荷) 头值 > 0."""
    n = 20
    dist = Distribution(
        x=np.zeros(n), y=np.zeros(n), z=np.linspace(-1e-3, 1e-3, n),
        px=np.zeros(n), py=np.zeros(n), pz=np.full(n, 1.0e9),
        clock=np.zeros(n), charge=np.full(n, -0.05),   # 总和 -1.0 nC
        status=np.full(n, 5),
        ref_time_ns=5.0, ref_momentum_eVc=1.0e9,
        ref_x_m=0.0, ref_y_m=0.0, ref_z_m=1.5,
    )
    p = tmp_path / "neg.001"
    write_distribution(dist, p, format="binary", include_index=True)
    raw = np.fromfile(p, dtype=np.float64)
    header_q = raw[2]
    assert header_q > 0
    assert header_q == pytest.approx(1.0)             # sum(|q|) = 1.0 nC


def test_write_binary_species_default_ones(tmp_path):
    """index=None 时种类列写 1 (电子), 非 1..N 序列 (P2-3)."""
    n = 12
    dist = Distribution(
        x=np.zeros(n), y=np.zeros(n), z=np.linspace(-1e-3, 1e-3, n),
        px=np.zeros(n), py=np.zeros(n), pz=np.full(n, 1.0e9),
        clock=np.zeros(n), charge=np.full(n, -0.05),
        status=np.full(n, 5),
        ref_time_ns=5.0, ref_momentum_eVc=1.0e9,
    )
    p = tmp_path / "spec.001"
    write_distribution(dist, p, format="binary", include_index=True)
    d2 = read_distribution(p)
    assert d2.index is not None
    assert np.all(d2.index == 1)                      # 全部电子
    assert d2.n_active == n


def test_write_ascii_first_row_is_reference_particle(tmp_path):
    """ASCII 写入首行必须是 10 列参考粒子行 (绝对坐标), 非 5 值头 (P2-1).
    用真实 ASTRA 文件 (Example.0150.001) 做输入, 首行应与粒子 0 一致."""
    dist = read_distribution(
        PROJECT_ROOT / "examples" / "Manual_Example" / "Example.0150.001")
    p = tmp_path / "out.txt"
    write_distribution(dist, p, format="ascii")
    with open(p) as fh:
        first = fh.readline().split()
        second = fh.readline().split()
    assert len(first) == 10                            # 参考粒子行
    assert len(second) == 10
    i0 = 0  # 粒子 0 = 参考粒子 (绝对坐标)
    ref = [dist.x[i0], dist.y[i0], dist.z[i0],
           dist.px[i0], dist.py[i0], dist.pz[i0],
           dist.clock[i0] * 1e9, dist.charge[i0],
           dist.index[i0], dist.status[i0]]
    got = [float(v) for v in first]
    assert got == pytest.approx(ref, rel=1e-9)
    # 第二行 z/pz/clock 相对参考粒子
    assert float(second[2]) == pytest.approx(dist.z[1] - dist.z[0], abs=1e-9)
    assert float(second[5]) == pytest.approx(dist.pz[1] - dist.pz[0], abs=1e-3)
    # 回读一致
    d2 = read_distribution(p)
    assert d2.n_particle == dist.n_particle
    assert np.allclose(d2.z, dist.z, rtol=1e-9, atol=1e-12)


def test_big_endian_binary_file_reads(tmp_path):
    """大端字节序二进制文件应被自动探测并正确读取 (P3 字节序)."""
    p = tmp_path / "be.001"
    header = [5.0035, 1.0005e9, 1.0, -6.8e-13, -3.9e-13]
    x, y, z, px, py, pz, clock = _make_10col_binary(
        p, header, n=30, species=1, status=5, dtype=">f8")
    dist = read_distribution(p)
    assert dist.n_particle == 30
    assert dist.n_active == 30
    assert dist.ref_momentum_eVc == pytest.approx(1.0005e9)
    assert np.allclose(dist.x, x, rtol=1e-12)
    assert np.allclose(dist.pz, pz + 1.0005e9, rtol=1e-9)


def _nozero_variants(n):
    """构造无 0x00 字节的二进制 9 列文件体.

    1.2345678901234567 = 3F F3 C0 CA 42 8C 59 FB (无 0x00 字节);
    乘 2**k 仅改指数、尾数不变, 仍无 0x00。任何 int32 范围内的整数
    float64 都含 0x00 字节, 故 status 列用非整数无 0x00 值 (本测试
    只考察格式检测的字节级解析, 不考察 status 语义)。
    """
    v = 1.2345678901234567
    body = np.empty((n, 9))
    for k in range(9):
        body[:, k] = v * (2.0 ** (k + 1)) * (1.0 + 1e-12 * np.arange(n))
    return body, v


def test_tiny_binary_file_without_zero_bytes_detected(tmp_path):
    """≤10 粒子且全文件不含 0x00 字节的二进制文件不得误判为 ASCII
    (P3 小文件误判; 值精心选择避免 0x00 字节)."""
    p = tmp_path / "tiny.001"
    n = 10
    body, v = _nozero_variants(n)
    # 头 5 值同样取无 0x00 字节的值: [t_ns, p_ref, Q, x_ref, y_ref]
    hdr = np.array([v, v * 2.0 ** 40, v * 2.0 ** 10, -v * 2.0 ** -20,
                    v * 2.0 ** -20])
    with open(p, "wb") as fh:
        hdr.astype(np.float64).tofile(fh)
        body.astype(np.float64).tofile(fh)
    raw = open(p, "rb").read()
    assert b"\x00" not in raw                     # 前提: 全文件无 0x00 字节
    dist = read_distribution(p)
    assert dist.n_particle == n
    assert dist.ref_momentum_eVc == pytest.approx(v * 2.0 ** 40)
    assert np.allclose(dist.x, body[:, 0], rtol=1e-12)


def test_binary_dat_suffix_probe_accepts(tmp_path):
    """probe 对无 0x00 字节的二进制 .dat 文件应识别为分布文件."""
    p = tmp_path / "tiny.dat"
    n = 10
    body, v = _nozero_variants(n)
    hdr = np.array([v, v * 2.0 ** 40, v * 2.0 ** 10, -v * 2.0 ** -20,
                    v * 2.0 ** -20])
    with open(p, "wb") as fh:
        hdr.astype(np.float64).tofile(fh)
        body.astype(np.float64).tofile(fh)
    assert b"\x00" not in open(p, "rb").read()
    dist = read_distribution(p)                    # probe 必须接受
    assert dist.n_particle == n
