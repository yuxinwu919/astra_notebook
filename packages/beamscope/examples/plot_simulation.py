#!/usr/bin/env python3
"""Plot ASTRA simulation results from simulation_files/.

Usage:
    cd /path/to/astra_notebook
    python packages/beamscope/examples/plot_simulation.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

from beamscope.io import read_distribution
from beamscope.analysis.statistics import compute_statistics, print_statistics
from beamscope.plot.comparison import plot_comparison
from beamscope.plot.dashboard import plot_dashboard
from beamscope.plot.phase_space import plot_transverse_phase_space, plot_phase_space
from beamscope.style.rcparams import set_publication_style

# ── Paths ──
SIM_DIR = Path(__file__).resolve().parent.parent.parent.parent / "simulation_files"
OUT_DIR = SIM_DIR / "plots"
OUT_DIR.mkdir(exist_ok=True)

set_publication_style(use_tex=False)

# ── 1. 加载所有相空间文件 ──
print("=" * 60)
print("  加载 ASTRA 模拟结果")
print("=" * 60)

# 初始分布
dist_ini = read_distribution(SIM_DIR / "Example.ini")
print(f"\n📥 初始分布: {dist_ini}")

# 追踪步骤输出 (Example.0050.001, Example.0100.001, Example.0150.001)
tracking_files = sorted(SIM_DIR.glob("Example.????.001"))
dists = {"initial": dist_ini}

for f in tracking_files:
    # 从文件名提取 z 位置：Example.0050.001 → z=0.50m
    z_pos = float(f.name.split(".")[1]) / 100
    label = f"z={z_pos:.2f}m"
    dist = read_distribution(f)
    dists[label] = dist
    print(f"📥 {label}: {dist}")

# ── 2. 统计量对比 ──
print("\n" + "=" * 60)
print("  束团统计量对比")
print("=" * 60)

for label, dist in dists.items():
    stats = compute_statistics(dist, label=label)
    print_statistics(stats, title=f"Beam Statistics — {label}")

# ── 3. 初始 vs 最终对比 ──
print("生成对比图...")

# 相空间对比
final_label = list(dists.keys())[-1]
fig = plot_comparison(
    {"initial": dists["initial"], final_label: dists[final_label]},
    plot_type="phase_space",
    plane="x",
)
fig.savefig(OUT_DIR / "comparison_phase_x.png", dpi=150)
print(f"  ✓ {OUT_DIR / 'comparison_phase_x.png'}")

fig = plot_comparison(
    {"initial": dists["initial"], final_label: dists[final_label]},
    plot_type="phase_space",
    plane="z",
)
fig.savefig(OUT_DIR / "comparison_phase_z.png", dpi=150)
print(f"  ✓ {OUT_DIR / 'comparison_phase_z.png'}")

# 统计量柱状图对比
fig = plot_comparison(dists, plot_type="statistics")
fig.savefig(OUT_DIR / "comparison_statistics.png", dpi=150)
print(f"  ✓ {OUT_DIR / 'comparison_statistics.png'}")

# ── 4. 初始分布综合仪表盘 ──
fig = plot_dashboard(dists)
fig.savefig(OUT_DIR / "dashboard.png", dpi=150)
print(f"  ✓ {OUT_DIR / 'dashboard.png'}")

# ── 5. 每个追踪步骤的横向相空间 ──
for label, dist in dists.items():
    safe_name = label.replace("=", "_").replace(".", "p")
    fig = plot_transverse_phase_space(dist, title_prefix=f"{label} — ")
    fig.savefig(OUT_DIR / f"phase_space_{safe_name}.png", dpi=150)
    print(f"  ✓ {OUT_DIR / f'phase_space_{safe_name}.png'}")

# ── 6. 纵向相空间演化 ──
for label, dist in dists.items():
    safe_name = label.replace("=", "_").replace(".", "p")
    fig = plot_phase_space(dist, plane="z", title=f"Longitudinal Phase Space — {label}")
    fig.savefig(OUT_DIR / f"longitudinal_{safe_name}.png", dpi=150)
    print(f"  ✓ {OUT_DIR / f'longitudinal_{safe_name}.png'}")

print(f"\n✅ 全部完成！图片保存在: {OUT_DIR.resolve()}")
