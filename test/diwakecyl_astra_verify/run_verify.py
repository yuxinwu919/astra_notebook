#!/usr/bin/env python3
"""
DiWakeCyl × ASTRA 尾场验证测试 (v2)
=====================================
验证 DiWakeCyl 生成的介质波导尾场文件能被 ASTRA 正确读取。

- 使用 Taylor_Method_F 格式（ASTRA 已验证兼容）
- 测试 m=0 (纵向) 和 m=1 (横向) 尾场效应
- 关注：纵向相空间、能散、发射度变化

用法:
  python run_verify.py
"""

import sys, os, time, subprocess, shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- 路径设置 ----
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
TEST_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "Wakefield" / "DiWakeCyl"))
sys.path.insert(0, str(PROJECT_ROOT))
import utils
from DWFA_cyl_func_Ng import *

ASTRA_EXE = str((PROJECT_ROOT / "astra" / "astra").resolve())
GENERATOR_EXE = str((PROJECT_ROOT / "astra" / "generator").resolve())

# ============================================================
# 参数定义
# ============================================================
BEAM_E      = 200.0       # MeV
BEAM_Q      = 0.05        # nC (50 pC)
BEAM_SIGZ   = 0.1         # mm RMS
BEAM_SIGE   = 10.0        # keV RMS
BEAM_EMIT_N = 1.0         # π·mrad·mm
BEAM_NPART  = 500

DLW_B       = 0.5e-3      # 内径 [m]
DLW_A       = 1.0e-3      # 外径 [m]
DLW_EPS     = 4.41        # 介电常数
DLW_MU      = 1.0         # 磁导率
DLW_NMODE   = 2           # 模式数（用主导模式）

ASTRA_ZSTART = 0.0
ASTRA_ZSTOP  = 0.1
ASTRA_WK_Z   = 0.05
ASTRA_N_BIN  = 80

# DiWakeCyl 网格
SIGMAZ_M   = BEAM_SIGZ * 1e-3
NZ         = 500         # 减少网格点数
ZMIN       = 0
ZMAX       = 16 * SIGMAZ_M

# 测试算例（仅 Taylor 格式）
CASES = [
    {
        "name": "case_longitudinal",
        "dir": TEST_DIR / "case_longitudinal",
        "m": 0,
        "desc": "纵向尾场 m=0 (Taylor)",
    },
    {
        "name": "case_transverse",
        "dir": TEST_DIR / "case_transverse",
        "m": 1,
        "desc": "横向尾场 m=1 (Taylor)",
    },
]

print("=" * 70)
print("  DiWakeCyl × ASTRA 尾场验证测试 v2")
print("=" * 70)
print(f"  束团: E={BEAM_E} MeV, Q={BEAM_Q} nC, σz={BEAM_SIGZ} mm")
print(f"  波导: b={DLW_B*1e3:.1f} mm, a={DLW_A*1e3:.1f} mm, ε={DLW_EPS}")
print(f"  ASTRA: {ASTRA_ZSTART}→{ASTRA_ZSTOP} m, Wk@z={ASTRA_WK_Z} m")
print(f"  DiWakeCyl: Nmode={DLW_NMODE}, Nz={NZ}")

# ============================================================
# Step 1: 生成初始束团
# ============================================================
print("\n" + "=" * 50)
print("📐 Step 1: 生成初始束团分布")
print("-" * 50)

gen_params = {
    "FNAME": "'test.ini'",
    "IPart": BEAM_NPART,
    "Species": "'electrons'",
    "Q_total": BEAM_Q,
    "Probe": True,
    "Noise_reduc": True,
    "Cathode": False,
    "Add": False, "N_add": 0,
    "Ref_Ekin": BEAM_E,
    "Ref_zpos": 0.0,
    "Dist_z": "'gauss'",
    "sig_z": BEAM_SIGZ,
    "C_sig_z": 3.0,
    "Dist_pz": "'g'",
    "sig_Ekin": BEAM_SIGE,
    "cor_Ekin": 0.0,
    "Dist_x": "'gauss'", "sig_x": 0.05,
    "Dist_y": "'gauss'", "sig_y": 0.05,
    "Dist_px": "'g'", "Nemit_x": BEAM_EMIT_N, "cor_px": 0.0,
    "Dist_py": "'g'", "Nemit_y": BEAM_EMIT_N, "cor_py": 0.0,
}

gen_in = TEST_DIR / "generator.in"
utils.write_namelist("INPUT", gen_params, gen_in)

result = subprocess.run(
    [GENERATOR_EXE, "generator.in"],
    cwd=str(TEST_DIR), capture_output=True, text=True, timeout=60
)
if result.returncode != 0:
    print(f"Generator stderr:\n{result.stderr}")
    raise RuntimeError("Generator 运行失败")

test_ini = TEST_DIR / "test.ini"
if not test_ini.exists():
    raise FileNotFoundError(f"Generator 输出未找到: {test_ini}")

dist_initial = utils.read_astra_distribution(test_ini)
print(f"  ✅ 初始束团: {dist_initial['n_particle']} particles")
print(f"     Q={dist_initial['total_charge']:.4f} nC")

# ============================================================
# Step 2: DiWakeCyl 计算尾场 + 无尾场对照组
# ============================================================
print("\n" + "=" * 50)
print("📐 Step 2: DiWakeCyl 计算格林函数")
print("-" * 50)

# 对照组：无尾场（纯漂移）
print("  设置对照组 (无尾场)...")
ctrl_dir = TEST_DIR / "case_control"
ctrl_dir.mkdir(parents=True, exist_ok=True)
shutil.copy2(test_ini, ctrl_dir / "test.ini")
(ctrl_dir / "astra.in").write_text(f"""&NEWRUN
  Head='Control (no wakefield)',
  Distribution='test.ini',
  RUN=1,
  H_max=0.0001,
  Track_All=T,
  Auto_Phase=F,
 /
&OUTPUT
  ZSTART={ASTRA_ZSTART},
  ZSTOP={ASTRA_ZSTOP},
  High_res=T,
  EmitS=T,
  Zemit=10,
  PhaseS=T,
  Zphase=1,
  RefS=T,
 /
""")

# m=0 纵场
print(f"  计算 m=0 (单极模式)...")
RootAmplit0, RootWavVec0 = FindMode(DLW_A, DLW_B, 0, DLW_MU, DLW_EPS, DLW_NMODE, 0.3)
r0_val = DLW_B
zGreen, WlGreen = Long_GreenFunction(
    RootAmplit0, RootWavVec0, r0_val, r0_val, DLW_B, DLW_A, 0,
    ZMIN, ZMAX, NZ, DLW_MU, DLW_EPS
)
print(f"    {len(RootAmplit0)} modes, Wl(0)={WlGreen[0]:.3e} V/m/C")

# m=1 横场
print(f"  计算 m=1 (偶极模式)...")
RootAmplit1, RootWavVec1 = FindMode(DLW_A, DLW_B, 1, DLW_MU, DLW_EPS, DLW_NMODE, 0.3)
zGreen_t, WtGreen = Trans_GreenFunction(
    RootAmplit1, RootWavVec1, r0_val, r0_val, DLW_B, DLW_A, 1,
    ZMIN, ZMAX, NZ, DLW_MU, DLW_EPS
)
print(f"    {len(RootAmplit1)} modes, Wt(max)={np.max(np.abs(WtGreen)):.3e} V/m/C")

# ============================================================
# Step 3: 绘制尾场图
# ============================================================
print("\n" + "=" * 50)
print("📐 Step 3: 绘制格林函数")
print("-" * 50)

fig_wake, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
ax1.plot(zGreen * 1e3, WlGreen * 1e-16, 'b-', linewidth=1.5)
ax1.set_xlabel('ζ [mm]')
ax1.set_ylabel('W_long [×10^16 V/m/C]')
ax1.set_title(f'Longitudinal Wake (m=0), Nmode={DLW_NMODE}')
ax1.grid(True, alpha=0.3)

ax2.plot(zGreen_t * 1e3, WtGreen * 1e-18, 'r-', linewidth=1.5)
ax2.set_xlabel('ζ [mm]')
ax2.set_ylabel('W_trans [×10^18 V/m/C]')
ax2.set_title(f'Transverse Wake (m=1), Nmode={DLW_NMODE}')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
fig_wake.savefig(TEST_DIR / "wakefield_green.png", dpi=150, bbox_inches='tight')
plt.close(fig_wake)
print(f"  ✅ 尾场图: {TEST_DIR / 'wakefield_green.png'}")

# ============================================================
# Step 4: 生成尾场文件和 astra.in (Taylor 格式)
# ============================================================
print("\n" + "=" * 50)
print("📐 Step 4: 生成尾场文件 & ASTRA 输入")
print("-" * 50)

ds = zGreen[1] - zGreen[0]

wake_data = {
    0: {"dir": TEST_DIR / "case_longitudinal", "data": (zGreen, WlGreen)},
    1: {"dir": TEST_DIR / "case_transverse", "data": (zGreen_t, WtGreen)},
}

for case in CASES:
    case_dir = case["dir"]
    case_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(test_ini, case_dir / "test.ini")

    m_val = case["m"]
    zg, wf = wake_data[m_val]["data"]

    wake_fname = f"wake_taylor.dat"
    wake_path = case_dir / wake_fname

    # Taylor 格式: 3 项展开
    ds_local = zg[1] - zg[0]
    Nterms = 3
    dw = np.gradient(wf, ds_local)
    d2w = np.gradient(dw, ds_local)

    with open(wake_path, 'w') as f:
        f.write(f'                    {Nterms}                   0\n')
        f.write(f'                 {NZ}                   0\n')
        f.write('                    0                   0\n')
        f.write('                    0                   0\n')
        for i in range(NZ):
            f.write(f'{zg[i]:25.15E}{wf[i]:25.15E}\n')
        # 一阶导数段
        f.write(f'                    0                 {NZ}\n')
        f.write('                    0                    0\n')
        f.write('                    0                   13\n')
        f.write('                    0                    0\n')
        for i in range(1, NZ):
            f.write(f'{zg[i]:25.15E}{dw[i]:25.15E}\n')
        # 二阶导数段
        f.write(f'                    0                 {NZ}\n')
        f.write('                    0                    0\n')
        f.write('                    0                   24\n')
        f.write('                    0                    0\n')
        for i in range(1, NZ):
            f.write(f'{zg[i]:25.15E}{d2w[i]:25.15E}\n')

    # 生成 astra.in
    if m_val == 0:
        wk_type = "'Taylor_Method_F'"
    else:
        wk_type = "'Taylor_Method_F'"

    astra_in = case_dir / "astra.in"
    astra_in.write_text(f"""&NEWRUN
  Head='{case["desc"]}',
  Distribution='test.ini',
  RUN=1,
  H_max=0.0001,
  Track_All=T,
  Auto_Phase=F,
 /
&OUTPUT
  ZSTART={ASTRA_ZSTART},
  ZSTOP={ASTRA_ZSTOP},
  High_res=T,
  EmitS=T,
  Zemit=10,
  PhaseS=T,
  Zphase=1,
  RefS=T,
 /
&WAKE
  LWake=T,
  WK_Z(1)={ASTRA_WK_Z},
  WK_EQUI_GRID(1)=1.0,
  WK_N_BIN(1)={ASTRA_N_BIN},
  WK_TYPE(1)={wk_type},
  WK_FILENAME(1)='{wake_fname}',
 /
""")

    print(f"  ✅ {case['name']}: {wake_fname} ({NZ} points, {Nterms} terms)")

# ============================================================
# Step 5: 运行 ASTRA
# ============================================================
print("\n" + "=" * 50)
print("📐 Step 5: 运行 ASTRA")
print("-" * 50)

all_runs = [
    ("control", ctrl_dir),
] + [(c["name"], c["dir"]) for c in CASES]

results = {}

for name, run_dir in all_runs:
    print(f"  {name}...", end=" ", flush=True)
    t0 = time.perf_counter()

    result = subprocess.run(
        [ASTRA_EXE, "astra.in"],
        cwd=str(run_dir), capture_output=True, text=True, timeout=120
    )
    elapsed = time.perf_counter() - t0

    if result.returncode != 0:
        print(f"❌ (返回码={result.returncode})")
        if result.stderr:
            # 只显示关键错误
            err_lines = result.stderr.strip().split('\n')
            for line in err_lines[:3]:
                print(f"     {line}")
        results[name] = {"status": "failed"}
        continue

    # 查找输出
    phase_files = sorted(run_dir.glob("astra.*.001"))
    phase_files = [f for f in phase_files if f.name.split('.')[1].isdigit()]
    if not phase_files:
        print("❌ 无输出")
        results[name] = {"status": "no_output"}
        continue

    astra_out = max(phase_files, key=lambda p: p.stat().st_size)
    try:
        dist_after = utils.read_astra_distribution(astra_out)
    except Exception as e:
        print(f"❌ 读取失败: {e}")
        results[name] = {"status": "read_error"}
        continue

    results[name] = {
        "status": "ok",
        "dist": dist_after,
        "elapsed": elapsed,
        "output": astra_out,
    }
    print(f"✅ ({elapsed:.2f}s, {dist_after['n_active']} active)")

# ============================================================
# Step 6: 分析对比
# ============================================================
print("\n" + "=" * 50)
print("📐 Step 6: 分析与绘图")
print("-" * 50)

ref_eV = BEAM_E * 1e6
stats_initial = utils.compute_beam_statistics(dist_initial, ref_energy_eV=ref_eV)

# 对照组统计
if "control" in results and results["control"]["status"] == "ok":
    stats_control = utils.compute_beam_statistics(results["control"]["dist"], ref_energy_eV=ref_eV)
else:
    stats_control = None

print("\n" + "=" * 100)
print("  Beam Statistics Comparison")
print("=" * 100)
hdr = f"  {'Case':<28s} {'σ(δp/p) [%]':>14s} {'ε_nx [μm]':>12s} {'ε_ny [μm]':>12s} {'σ_z [mm]':>10s} {'Mean pz [MeV]':>14s}"
print(hdr)
print("  " + "-" * 96)
print(f"  {'Initial':<28s} {stats_initial['sig_E_over_E']*100:>14.4f} {stats_initial['emit_x_norm']*1e6:>12.4f} {stats_initial['emit_y_norm']*1e6:>12.4f} {stats_initial['sig_z']*1e3:>10.4f} {stats_initial['mean_pz']/1e6:>14.4f}")

for name in ["control"] + [c["name"] for c in CASES]:
    if name in results and results[name]["status"] == "ok":
        s = utils.compute_beam_statistics(results[name]["dist"], ref_energy_eV=ref_eV)
        print(f"  {name:<28s} {s['sig_E_over_E']*100:>14.4f} {s['emit_x_norm']*1e6:>12.4f} {s['emit_y_norm']*1e6:>12.4f} {s['sig_z']*1e3:>10.4f} {s['mean_pz']/1e6:>14.4f}")
print("=" * 100 + "\n")

# ---- 逐个算例绘图 ----
for name in ["control"] + [c["name"] for c in CASES]:
    if name not in results or results[name]["status"] != "ok":
        print(f"  ⊘ {name}: 跳过")
        continue

    case_dir = ctrl_dir if name == "control" else {c["name"]: c["dir"] for c in CASES}[name]
    r = results[name]
    dist_a = r["dist"]
    s = utils.compute_beam_statistics(dist_a, ref_energy_eV=ref_eV)

    mask_i = dist_initial['active_mask']
    mask_a = dist_a['active_mask']

    pz_i = dist_initial['pz'][mask_i]
    pz_a = dist_a['pz'][mask_a]
    mean_pz = np.mean(pz_i)
    dp_i = (pz_i - mean_pz) / mean_pz * 100
    dp_a = (pz_a - mean_pz) / mean_pz * 100
    z_i = dist_initial['z'][mask_i] * 1e3
    z_a = dist_a['z'][mask_a] * 1e3
    x_i = dist_initial['x'][mask_i] * 1e3
    y_i = dist_initial['y'][mask_i] * 1e3
    x_a = dist_a['x'][mask_a] * 1e3
    y_a = dist_a['y'][mask_a] * 1e3

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    z_lim = max(np.percentile(np.abs(np.concatenate([z_i, z_a])), 99), 0.02) * 1.5
    dp_lim = max(np.percentile(np.abs(np.concatenate([dp_i, dp_a])), 99), 0.002) * 1.2

    # Before
    ax = axes[0, 0]
    h = ax.hist2d(z_i, dp_i, bins=(80, 80), cmap='Reds',
                  range=[[-z_lim, z_lim], [-dp_lim, dp_lim]])
    plt.colorbar(h[3], ax=ax, label='Count')
    ax.set_xlabel('z [mm]'); ax.set_ylabel('δp/p [%]')
    ax.set_title('Before: z vs δp/p', fontweight='bold')
    ax.axhline(0, color='gray', ls='--', lw=0.5)
    ax.axvline(0, color='gray', ls='--', lw=0.5)
    ax.text(0.03, 0.97, f'σ={np.std(dp_i):.4f}%', transform=ax.transAxes,
            fontsize=8, va='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # After
    ax = axes[0, 1]
    h = ax.hist2d(z_a, dp_a, bins=(80, 80), cmap='Blues',
                  range=[[-z_lim, z_lim], [-dp_lim, dp_lim]])
    plt.colorbar(h[3], ax=ax, label='Count')
    ax.set_xlabel('z [mm]'); ax.set_ylabel('δp/p [%]')
    ax.set_title(f'After ({name}): z vs δp/p', fontweight='bold')
    ax.axhline(0, color='gray', ls='--', lw=0.5)
    ax.axvline(0, color='gray', ls='--', lw=0.5)
    ax.text(0.03, 0.97, f'σ={np.std(dp_a):.4f}%', transform=ax.transAxes,
            fontsize=8, va='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # δp/p 分布对比
    ax = axes[0, 2]
    ax.hist(dp_i, bins=50, density=True, alpha=0.6, color='#F44336',
            edgecolor='white', label=f'Before (σ={np.std(dp_i):.4f}%)', histtype='stepfilled')
    ax.hist(dp_a, bins=50, density=True, alpha=0.5, color='#2196F3',
            edgecolor='white', label=f'After  (σ={np.std(dp_a):.4f}%)', histtype='stepfilled')
    ax.set_xlabel('δp/p [%]'); ax.set_ylabel('Density')
    ax.set_title('Energy Spread'); ax.legend(fontsize=9)

    # x 分布
    ax = axes[1, 0]
    ax.hist(x_i, bins=60, density=True, alpha=0.6, color='#F44336',
            edgecolor='white', label=f'Before (σ={np.std(x_i):.4f}mm)', histtype='stepfilled')
    ax.hist(x_a, bins=60, density=True, alpha=0.5, color='#2196F3',
            edgecolor='white', label=f'After  (σ={np.std(x_a):.4f}mm)', histtype='stepfilled')
    ax.set_xlabel('x [mm]'); ax.set_ylabel('Density')
    ax.set_title('x Distribution'); ax.legend(fontsize=9)

    # y 分布
    ax = axes[1, 1]
    ax.hist(y_i, bins=60, density=True, alpha=0.6, color='#F44336',
            edgecolor='white', label=f'Before (σ={np.std(y_i):.4f}mm)', histtype='stepfilled')
    ax.hist(y_a, bins=60, density=True, alpha=0.5, color='#2196F3',
            edgecolor='white', label=f'After  (σ={np.std(y_a):.4f}mm)', histtype='stepfilled')
    ax.set_xlabel('y [mm]'); ax.set_ylabel('Density')
    ax.set_title('y Distribution'); ax.legend(fontsize=9)

    # 统计
    ax = axes[1, 2]
    ax.axis('off')
    lines = [
        f"Case: {name}",
        f"Particles: {dist_a['n_active']}/{dist_a['n_particle']}",
        f"Elapsed: {r['elapsed']:.2f}s",
        "",
        "Parameter      Before → After",
        f"σ(δp/p) [%]   {np.std(dp_i):.4f} → {np.std(dp_a):.4f}",
        f"ε_nx [μm]      {stats_initial['emit_x_norm']*1e6:.2f} → {s['emit_x_norm']*1e6:.2f}",
        f"ε_ny [μm]      {stats_initial['emit_y_norm']*1e6:.2f} → {s['emit_y_norm']*1e6:.2f}",
        f"σ_z [mm]       {stats_initial['sig_z']*1e3:.4f} → {s['sig_z']*1e3:.4f}",
        f"Mean pz [MeV]  {stats_initial['mean_pz']/1e6:.2f} → {s['mean_pz']/1e6:.2f}",
    ]
    for i, txt in enumerate(lines):
        ax.text(0.05, 0.95 - i * 0.06, txt, transform=ax.transAxes, fontsize=8, family='monospace')

    plt.suptitle(f"{name} — Before vs After", fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig_path = case_dir / "result_before_after.png"
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✅ {name}: {fig_path}")

# ---- 汇总图 ----
ok_cases = [(n, results[n]) for n in ["control"] + [c["name"] for c in CASES]
            if n in results and results[n]["status"] == "ok"]
if len(ok_cases) >= 2:
    ncols = len(ok_cases)
    fig_sum, axes_sum = plt.subplots(1, ncols, figsize=(5 * ncols, 5))
    if ncols == 1:
        axes_sum = [axes_sum]

    for ax, (name, r) in zip(axes_sum, ok_cases):
        dist_d = r["dist"]
        mask = dist_d['active_mask']
        dp_d = (dist_d['pz'][mask] - mean_pz) / mean_pz * 100
        z_d = dist_d['z'][mask] * 1e3
        h = ax.hist2d(z_d, dp_d, bins=(80, 80), cmap='viridis',
                      range=[[-z_lim, z_lim], [-dp_lim, dp_lim]])
        plt.colorbar(h[3], ax=ax, label='Count')
        ax.set_xlabel('z [mm]'); ax.set_ylabel('δp/p [%]')
        ax.set_title(f'{name}', fontweight='bold')
        ax.axhline(0, color='white', ls='--', lw=0.5)
        ax.axvline(0, color='white', ls='--', lw=0.5)

    plt.suptitle('All Cases — Longitudinal Phase Space', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig_sum.savefig(TEST_DIR / "summary_all_cases.png", dpi=150, bbox_inches='tight')
    plt.close(fig_sum)
    print(f"  ✅ 汇总图: {TEST_DIR / 'summary_all_cases.png'}")

# ============================================================
# Step 7: 生成 README
# ============================================================
print("\n" + "=" * 50)
print("📐 Step 7: 生成 README.md")
print("-" * 50)

ok_names = [n for n in ["control"] + [c["name"] for c in CASES]
            if n in results and results[n]["status"] == "ok"]

readme = [
    "# DiWakeCyl × ASTRA 尾场验证测试",
    "",
    f"**测试时间**: {time.strftime('%Y-%m-%d %H:%M:%S')}",
    "",
    "## 参数",
    "",
    "| 参数 | 值 |",
    "|------|------|",
    f"| 束团能量 | {BEAM_E} MeV |",
    f"| 束团电荷 | {BEAM_Q} nC (50 pC) |",
    f"| RMS 长度 σz | {BEAM_SIGZ} mm |",
    f"| 能散 σE | {BEAM_SIGE} keV RMS |",
    f"| 归一化发射度 | {BEAM_EMIT_N} π·mrad·mm |",
    f"| 波导内径 b | {DLW_B*1e3:.1f} mm |",
    f"| 波导外径 a | {DLW_A*1e3:.1f} mm |",
    f"| 介电常数 ε | {DLW_EPS} |",
    f"| 模式数 | {DLW_NMODE} |",
    f"| ASTRA 追踪 | {ASTRA_ZSTART} → {ASTRA_ZSTOP} m |",
    f"| 尾场位置 Wk_z | {ASTRA_WK_Z} m |",
    f"| Bin 数 | {ASTRA_N_BIN} |",
    f"| 尾场格式 | Taylor_Method_F, 3 terms |",
    "",
    "## 结果",
    "",
    f"成功: {len(ok_names)} 个",
    "",
    "### 统计对比",
    "",
    "| Case | σ(δp/p) [%] | ε_nx [μm] | ε_ny [μm] | σ_z [mm] | Mean pz [MeV] |",
    "|------|-------------|-----------|-----------|---------|---------------|",
]

for name in ["control"] + [c["name"] for c in CASES]:
    if name in results and results[name]["status"] == "ok":
        s = utils.compute_beam_statistics(results[name]["dist"], ref_energy_eV=ref_eV)
        readme.append(
            f"| {name} | {s['sig_E_over_E']*100:.4f} | {s['emit_x_norm']*1e6:.2f} | {s['emit_y_norm']*1e6:.2f} | {s['sig_z']*1e3:.4f} | {s['mean_pz']/1e6:.2f} |"
        )
    else:
        readme.append(f"| {name} | N/A | N/A | N/A | N/A | N/A |")

readme += [
    "",
    "### 注意事项",
    "",
    "- DiWakeCyl 计算的是介质的格林函数（点电荷 δ-响应），而非束团总尾场势",
    "- ASTRA 内部会自行完成卷积（binning → 离散卷积 → 粒子 kick）",
    "- 介质波导的尾场格林函数值量级约 10^16 V/m/C，这是由物理过程决定的",
    "- Taylor_Method_F 格式已在实际测试中验证兼容",
    "",
    "## 文件结构",
    "",
    "```",
    "test/diwakecyl_astra_verify/",
    "├── README.md",
    "├── run_verify.py",
    "├── wakefield_green.png",
    "├── summary_all_cases.png",
    "├── case_control/          # 对照组 (无尾场)",
    "├── case_longitudinal/     # m=0 纵向尾场",
    "└── case_transverse/       # m=1 横向尾场",
    "```",
]

(TEST_DIR / "README.md").write_text("\n".join(readme) + "\n")
print(f"  ✅ README: {TEST_DIR / 'README.md'}")

# ============================================================
# 完成
# ============================================================
print("\n" + "=" * 70)
print(f"  ✅ 测试完成")
print(f"  成功: {len(ok_names)}/{len(all_runs)}")
for name in all_runs:
    s = "✅" if name[0] in ok_names else "❌"
    print(f"    {s} {name[0]}")
print(f"  工作目录: {TEST_DIR}")
print("=" * 70)