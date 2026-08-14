#!/usr/bin/env python3
"""beamscope 快速测试脚本 — 载入仿真数据并生成验证图表。

Usage:
    cd /Users/yuxinwu/astra_notebook
    python packages/beamscope/tests/quick_test.py

若无命令行参数，自动搜索以下仿真文件：
  1. test/diwakecyl_astra_test/astra.0210.001 (二进制输出, ~2000 particles)
  2. examples/Manual_Example/Example.ini (ASCII 输入, ~200 particles)
"""

import sys
from pathlib import Path
import warnings

# Allow running from any directory
_PROJ_ROOT = Path(__file__).resolve().parent.parent  # packages/beamscope
_NOTEBOOK_ROOT = _PROJ_ROOT.parent.parent  # astra_notebook root
sys.path.insert(0, str(_PROJ_ROOT))

import matplotlib
matplotlib.use("Agg")  # 无 GUI 后端
import matplotlib.pyplot as plt
import numpy as np

from beamscope.io import read_distribution
from beamscope.distribution import Distribution
from beamscope.analysis.statistics import compute_statistics, print_statistics
from beamscope.plot.overview import plot_overview
from beamscope.plot._precompute import precompute, clip_percentile, get_overview_panels
from beamscope.plot.detail import plot_detail
from beamscope.plot.dashboard import plot_dashboard
from beamscope.plot.phase_space import plot_phase_space, plot_transverse_phase_space
from beamscope.style.rcparams import set_publication_style

# 输出目录
_OUTDIR = _PROJ_ROOT / "tests" / "output"
_OUTDIR.mkdir(parents=True, exist_ok=True)

# 查找仿真文件的候选路径
_FILE_CANDIDATES = [
    # 二进制 ASTRA 输出（较大，~2000-10000 particles）
    _NOTEBOOK_ROOT / "test" / "diwakecyl_astra_test" / "astra.0210.001",
    # ASCII .ini 文件（较小，~200 particles）
    _NOTEBOOK_ROOT / "examples" / "Manual_Example" / "Example.ini",
    # 另一个二进制文件（Wake test）
    _NOTEBOOK_ROOT / "test" / "diwakecyl_astra_test" / "astra_tesla.0210.002",
    # chirper test
    _NOTEBOOK_ROOT / "test" / "dechirper_test" / "astra.0210.001",
    # Manual Example 二进制输出
    _NOTEBOOK_ROOT / "examples" / "Manual_Example" / "Example.0150.001",
]


def find_file():
    """Search for a loadable simulation file."""
    for p in _FILE_CANDIDATES:
        if p.exists():
            print(f"✓ 找到文件: {p}")
            return p
    return None


def run_all_tests(dist: Distribution, label: str):
    """Run every plot function and save output."""
    set_publication_style(use_tex=False)
    print(f"\n{'='*60}")
    print(f"  测试 Distribution: {label}")
    print(f"  Particles: {dist.n_particle} total, {dist.n_active} active")
    print(f"  Total charge: {dist.active_charge_nC:.4f} nC")
    print(f"{'='*60}\n")

    tests_passed = 0
    tests_failed = []

    # ── Test 1: clip_percentile ──
    print("Test 1: clip_percentile ...", end=" ")
    try:
        data = precompute(dist)
        for _, _, x_key, y_key, *_ in get_overview_panels():
            x_data = data.get(x_key, np.array([]))
            y_data = data.get(y_key, np.array([]))
            if len(x_data) > 0:
                vx = clip_percentile(x_data)
                vy = clip_percentile(y_data)
                assert vx[0] < vx[1], f"Invalid x range: {vx}"
                assert vy[0] < vy[1], f"Invalid y range: {vy}"
        print("PASS")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: {e}")
        tests_failed.append(("clip_percentile", str(e)))

    # ── Test 2: plot_overview ──
    print("Test 2: plot_overview (3×2 grid) ...", end=" ")
    try:
        fig, axes = plot_overview(dist, title=f"Overview — {label}")
        fig.savefig(_OUTDIR / f"overview_{label}.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"PASS → overview_{label}.png")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: {e}")
        tests_failed.append(("plot_overview", str(e)))

    # ── Test 3: plot_detail (x-x') ──
    print("Test 3: plot_detail (x-x' with ellipse) ...", end=" ")
    try:
        fig = plt.figure(figsize=(8, 7))
        plot_detail(fig, dist, x_key="x", y_key="xp",
                     title=f"Detail x-x' — {label}",
                     show_ellipse=True, show_marginals=True)
        fig.savefig(_OUTDIR / f"detail_x_xp_{label}.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"PASS → detail_x_xp_{label}.png")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: {e}")
        tests_failed.append(("plot_detail_x_xp", str(e)))

    # ── Test 4: plot_detail (z-dp) ──
    print("Test 4: plot_detail (z-dp, longitudinal) ...", end=" ")
    try:
        fig = plt.figure(figsize=(8, 7))
        plot_detail(fig, dist, x_key="z", y_key="dp",
                     title=f"Detail z-dp — {label}",
                     show_ellipse=True, show_marginals=True)
        fig.savefig(_OUTDIR / f"detail_z_dp_{label}.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"PASS → detail_z_dp_{label}.png")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: {e}")
        tests_failed.append(("plot_detail_z_dp", str(e)))

    # ── Test 5: plot_dashboard ──
    print("Test 5: plot_dashboard (multi-panel) ...", end=" ")
    try:
        fig = plot_dashboard({label: dist}, title=f"Dashboard — {label}")
        fig.savefig(_OUTDIR / f"dashboard_{label}.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"PASS → dashboard_{label}.png")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: {e}")
        tests_failed.append(("plot_dashboard", str(e)))

    # ── Test 6: plot_phase_space (x) ──
    print("Test 6: plot_phase_space (x, density) ...", end=" ")
    try:
        fig = plot_phase_space(dist, plane="x", kind="density")
        fig.savefig(_OUTDIR / f"phase_x_{label}.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"PASS → phase_x_{label}.png")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: {e}")
        tests_failed.append(("plot_phase_space_x", str(e)))

    # ── Test 7: plot_phase_space (z) ──
    print("Test 7: plot_phase_space (z, density) ...", end=" ")
    try:
        fig = plot_phase_space(dist, plane="z", kind="density")
        fig.savefig(_OUTDIR / f"phase_z_{label}.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"PASS → phase_z_{label}.png")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: {e}")
        tests_failed.append(("plot_phase_space_z", str(e)))

    # ── Test 8: plot_transverse_phase_space ──
    print("Test 8: plot_transverse_phase_space ...", end=" ")
    try:
        fig = plot_transverse_phase_space(dist)
        fig.savefig(_OUTDIR / f"transverse_{label}.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"PASS → transverse_{label}.png")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: {e}")
        tests_failed.append(("plot_transverse", str(e)))

    # ── Test 9: Statistics ──
    print("Test 9: compute_statistics ...", end=" ")
    try:
        stats = compute_statistics(dist)
        assert stats.n_active > 0, "No active particles"
        assert stats.sig_x >= 0, "Negative sigma_x"
        assert stats.sig_y >= 0, "Negative sigma_y"
        d = stats.to_dict()
        assert "sig_x_mm" in d, "Missing sig_x_mm"
        print_statistics(stats, title=f"Stats — {label}")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: {e}")
        tests_failed.append(("statistics", str(e)))

    # ── Summary ──
    print(f"\n{'─'*40}")
    print(f"  Results: {tests_passed}/9 passed")
    if tests_failed:
        print(f"  Failed: {[t[0] for t in tests_failed]}")
        for name, err in tests_failed:
            print(f"    - {name}: {err[:80]}")
    else:
        print("  ✅ All tests passed!")
    print(f"  Output: {_OUTDIR.resolve()}")
    print(f"{'─'*40}\n")
    return tests_passed, tests_failed


def main():
    warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    filepath = None
    if len(sys.argv) > 1:
        filepath = Path(sys.argv[1])
        if not filepath.exists():
            print(f"✗ 文件不存在: {filepath}")
            sys.exit(1)
    else:
        filepath = find_file()

    if filepath is None:
        print("✗ 未找到仿真文件。请将 .ini 或 .001 文件作为参数传入。")
        print("  Usage: python quick_test.py <path_to_file>")
        print(f"  Searched: {[_FILE_CANDIDATES]}")
        sys.exit(1)

    # Load
    print(f"\n正在载入: {filepath}")
    try:
        dist = read_distribution(filepath)
    except Exception as e:
        print(f"✗ 载入失败: {e}")
        import traceback; traceback.print_exc()
        sys.exit(1)

    label = filepath.stem
    n_pass, n_fail = run_all_tests(dist, label)

    # If more files are available, also load the secondary file
    for extra in _FILE_CANDIDATES:
        if extra == filepath or not extra.exists():
            continue
        label2 = extra.stem
        if label2 == label:
            continue
        try:
            dist2 = read_distribution(extra)
            run_all_tests(dist2, label2)
        except Exception:
            pass
        break  # Only test one extra at most

    sys.exit(0 if not n_fail else 1)


if __name__ == "__main__":
    main()