#!/usr/bin/env python3
"""测试 examples/ 下三个算例，分别用模式 A（导入 .in）和模式 B（字典生成）运行 ASTRA。"""

import shutil, subprocess, sys, time
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
SIM_DIR = PROJECT_DIR / "simulation_files"
ASTRA_EXE = str(PROJECT_DIR / "ASTRA" / "astra")
GEN_EXE = str(PROJECT_DIR / "ASTRA" / "generator")

def run_cmd(exe, work_dir, input_file=None):
    cmd = [exe, input_file] if input_file else [exe]
    t0 = time.perf_counter()
    r = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True, timeout=120)
    dt = time.perf_counter() - t0
    if r.returncode != 0:
        print(f"  ✗ 失败 (rc={r.returncode})")
        if r.stderr: print(f"  stderr: {r.stderr[:500]}")
        return None
    print(f"  ✓ 成功 ({dt:.1f}s)")
    return r.stdout

def save_output(src_dir, dst_dir):
    dst_dir.mkdir(parents=True, exist_ok=True)
    for f in src_dir.iterdir():
        if f.is_file():
            shutil.copy2(f, dst_dir / f.name)
    print(f"  → 输出保存到 {dst_dir}")

# ============================================================
# 1. Manual_Example
# ============================================================
print("="*60)
print("1. Manual_Example")
print("="*60)

ex_dir = PROJECT_DIR / "examples/Manual_Example"
ex_in = ex_dir / "Example.in"
gen_in = ex_dir / "generator.in"

# Step: 生成 Example.ini
SIM_DIR.mkdir(parents=True, exist_ok=True)
for f in SIM_DIR.iterdir(): f.unlink()
shutil.copy2(gen_in, SIM_DIR / "generator.in")
print("  运行 Generator...")
run_cmd(GEN_EXE, SIM_DIR, input_file=str(SIM_DIR / "generator.in"))

# -- 模式 A：直接导入 Example.in --
print("\n  --- 模式 A：导入 Example.in ---")
for f in SIM_DIR.iterdir(): 
    if f.name != "Example.ini": f.unlink()
shutil.copy2(ex_in, SIM_DIR / "astra.in")
run_cmd(ASTRA_EXE, SIM_DIR, input_file="astra.in")
save_output(SIM_DIR, ex_dir / "output_A")

# -- 模式 B：字典生成等效参数 --
print("\n  --- 模式 B：字典生成 ---")
for f in SIM_DIR.iterdir(): 
    if f.name != "Example.ini": f.unlink()

# 用 Python 生成完全等效的 astra.in
import sys; sys.path.insert(0, str(PROJECT_DIR))
import utils

params = {
    "NEWRUN": {
        "Head": "' Example of ASTRA users manual'",
        "RUN": 1,
        "Distribution": "'Example.ini'",
        "Xoff": 0.0, "Yoff": 0.0,
        "Track_All": True, "Auto_Phase": True,
        "H_max": 0.001, "H_min": 0.0,
        "Loop": False, "Nloop": 1, "Phase_Scan": False, "check_ref_part": False,
    },
    "OUTPUT": {
        "ZSTART": 0.0, "ZSTOP": 1.5,
        "Zemit": 500, "Zphase": 1,
        "RefS": True, "EmitS": True, "PhaseS": True,
    },
    "CHARGE": {
        "LSPCH": False, "Nrad": 10, "Cell_var": 2.0, "Nlong_in": 10,
        "min_grid": 0.0, "Max_Scale": 0.05,
    },
    "CAVITY": {
        "LEfield": True,
        "File_Efield(1)": "'3_cell_L-Band.dat'", "C_pos(1)": 0.3,
        "Nue(1)": 1.3, "MaxE(1)": 40.0, "Phi(1)": 0.0,
    },
    "SOLENOID": {
        "LBField": True,
        "File_Bfield(1)": "'Solenoid.dat'", "S_pos(1)": 1.2,
        "MaxB(1)": 0.35, "S_smooth(1)": 10,
    },
}

lines = []
for name, p in params.items():
    lines.append(utils.write_namelist(name, p, None))
with open(SIM_DIR / "astra.in", "w") as f:
    f.write("\n".join(lines))
print("  astra.in 已生成")
run_cmd(ASTRA_EXE, SIM_DIR, input_file="astra.in")
save_output(SIM_DIR, ex_dir / "output_B")

# ============================================================
# 2. Wake
# ============================================================
print("\n" + "="*60)
print("2. Wake")
print("="*60)

wake_dir = PROJECT_DIR / "examples/Wake/Wake_Files"
wake_in = wake_dir / "Wake.in"
ini_file = wake_dir / "test.ini"

SIM_DIR.mkdir(parents=True, exist_ok=True)
for f in SIM_DIR.iterdir(): f.unlink()
shutil.copy2(ini_file, SIM_DIR / "test.ini")

# -- 模式 A --
print("\n  --- 模式 A：导入 Wake.in ---")
shutil.copy2(wake_in, SIM_DIR / "astra.in")
run_cmd(ASTRA_EXE, SIM_DIR, input_file="astra.in")
save_output(SIM_DIR, PROJECT_DIR / "examples/Wake/output_A")

# -- 模式 B --
print("\n  --- 模式 B：字典生成 ---")
for f in SIM_DIR.iterdir(): 
    if f.name != "test.ini": f.unlink()

params = {
    "NEWRUN": {
        "Head": "'DRIFT to z=2.1 with WAKE at 2.0m (from z=1.9)'",
        "Distribution": "'test.ini'", "RUN": 1,
        "H_max": 0.00333,
        "Loop": False, "Nloop": 1, "Xoff": 0.0, "Yoff": 0.0,
        "Track_All": True, "Phase_Scan": False, "Auto_Phase": True,
        "check_ref_part": False, "H_min": 0.0,
    },
    "OUTPUT": {
        "ZSTART": 1.9, "ZSTOP": 2.1,
        "High_res": True, "EmitS": True, "Zemit": 10,
        "PhaseS": True, "Zphase": 1, "RefS": True,
    },
    "WAKE": {
        "LWake": True,
        "Wk_z(1)": 2.0, "Wk_equi_grid(1)": 1.0,
        "Wk_N_bin(1)": 50, "Wk_Type(1)": "'taylor_method_F'",
        "Wk_filename(1)": "'TESLA_MODULE_WAKE_TAYLOR.dat'",
        "Wk_testfile(1)": "'test.dat'",
    },
}

lines = []
for name, p in params.items():
    lines.append(utils.write_namelist(name, p, None))
with open(SIM_DIR / "astra.in", "w") as f:
    f.write("\n".join(lines))
print("  astra.in 已生成")
run_cmd(ASTRA_EXE, SIM_DIR, input_file="astra.in")
save_output(SIM_DIR, PROJECT_DIR / "examples/Wake/output_B")

# ============================================================
# 3. Aperture
# ============================================================
print("\n" + "="*60)
print("3. Aperture")
print("="*60)
print("  ⚠ Aperture 算例仅有 &APERTURE 块，无 NEWRUN/OUTPUT，无法独立运行 ASTRA。")
print("  跳过。")

print("\n" + "="*60)
print("全部测试完成！")
print("="*60)
