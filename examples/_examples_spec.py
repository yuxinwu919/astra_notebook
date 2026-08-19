"""官方算例共享规格: 文件清单/运行步骤/黄金比对目标。

供 examples/<name>/<name>.ipynb 单算例示例 notebook 共用 - 单一
数据源, 避免 9 份拷贝漂移 (批 5: 06_examples 汇总本已删除)。
"""

import json
import shutil
from pathlib import Path

import numpy as np

from astra_tools.run import run_program
from astra_tools.io.astra_emit import parse_output_file
from astra_tools.io.field_map import fix_laser_map_header


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


PROJECT_ROOT = project_root()
EXAMPLES_DIR = PROJECT_ROOT / "examples"
SIM_DIR = PROJECT_ROOT / "data" / "workspace"

GOLDEN_EXPECTED = json.loads((EXAMPLES_DIR / "golden_expected.json").read_text())

# 每个算例: 输入文件清单 / 运行步骤 / 黄金比对目标 / 一句话简介
EXAMPLES = {
    "Manual_Example": dict(
        title="手册主算例: 光阴极枪 + 螺线管, 发射度补偿",
        stems=["Example"],
        copy=["generator.in", "Example.in", "3_cell_L-Band.dat", "Solenoid.dat"],
        patch={"Example.in": (
            "/Users/yuxinwu/astra_notebook/simulation_files/Example.ini",
            "Example.ini")},
        steps=[("generator", "generator.in"), ("astra", "Example.in")],
        golden_xemit=EXAMPLES_DIR / "Manual_Example/Example.Xemit.001"),
    "Aperture": dict(
        title="孔径: 圆孔 + 圆柱堵块",
        stems=["astra"],
        copy=["astra.in", "aperture.in", "Geometry.dat", "test.ini"],
        steps=[("astra", "astra.in")],
        golden_xemit=EXAMPLES_DIR / "Aperture/golden/astra.Xemit.001"),
    "Wake": dict(
        title="尾场: TESLA 模块尾场表",
        stems=["Wake"],
        copy=["Wake.in", "test.ini", "TESLA_MODULE_WAKE_TAYLOR.dat", "test.dat"],
        src_dir="Wake/Wake_Files",
        steps=[("astra", "Wake.in")],
        golden_xemit=EXAMPLES_DIR / "Wake/golden/Wake.Xemit.001"),
    "Cavity_Example": dict(
        title="腔与 TWS: 1D 腔表/TWS/3D 场图/TDS",
        stems=["astra"],
        copy=["generator.in", "astra.in", "TWS_Sband.dat", "3_cell_L-Band.dat",
              "dcfield.dat", "3D_test.bx", "3D_test.by", "3D_test.bz",
              "3D_test.ex", "3D_test.ey", "3D_test.ez"],
        steps=[("generator", "generator.in"), ("astra", "astra.in")],
        golden_xemit=EXAMPLES_DIR / "Cavity_Example/golden/astra.Xemit.001"),
    "Curved_Cathode_Example": dict(
        title="弯曲阴极",
        stems=["astra"],
        copy=["generator.in", "astra.in", "Contour.dat", "efld.dat"],
        steps=[("generator", "generator.in"), ("astra", "astra.in")],
        golden_xemit=EXAMPLES_DIR / "Curved_Cathode_Example/golden/astra.Xemit.001"),
    "90deg_bend_Example": dict(
        title="90° 弯转: 3D 二极场, 分段追踪",
        stems=["Section1"],
        copy=["Section1.in", "Section2.in", "test.ini",
              "3D_Dipole.bx", "3D_Dipole.by", "3D_Dipole.bz"],
        patch={"Section2.in": ("Section1_n.0100.001", "Section1.0100.001")},
        steps=[("astra", "Section1.in"), ("astra", "Section2.in")],
        golden_xemit=EXAMPLES_DIR / "90deg_bend_Example/golden/Section2.Log.001",
        compare_mode="log"),
    "Plasma_Example_1": dict(
        title="束驱动等离子体尾场",
        stems=["plasma"],
        copy=["phsp.in", "plasma.in", "PLASMA_flattop.txt"],
        steps=[("generator", "phsp.in"), ("astra", "plasma.in")],
        golden_xemit=EXAMPLES_DIR / "Plasma_Example_1/golden/plasma.Xemit.001"),
    "Plasma_Example_2": dict(
        title="激光驱动等离子体尾场 (需本地 laser.dat)",
        stems=["plasma"],
        copy=["phsp.in", "plasma.in", "PLASMA_flattop.txt", "laser.dat"],
        laser_fix=True,
        steps=[("generator", "phsp.in"), ("astra", "plasma.in")],
        golden_xemit=EXAMPLES_DIR / "Plasma_Example_2/golden/plasma.Xemit.001"),
}


def stage_files(name):
    """复制官方算例输入到工作目录 (不运行), 返回 (work, spec)。

    先清空工作目录: 上次运行的陈旧输出 (如 ZSTOP 改版前遗留的
    相空间 dump) 会污染 phase_files/统计/比对 (二轮审计 R2-3-1
    实测: 残留 astra.0190.001 使 ph[-1] 取到旧 dump, lost 误判为 0)。
    跨步骤算例 (90deg_bend: Section1 -> Section2) 在同一目录内连续
    运行, 清空只发生在 stage 时, 不影响步骤间产物。
    """
    spec = EXAMPLES[name]
    src = EXAMPLES_DIR / spec.get("src_dir", name)
    work = SIM_DIR / name
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    for f in spec["copy"]:
        s = src / f
        if not s.exists():
            raise FileNotFoundError("%s 缺失: %s" % (f, s))
        shutil.copy2(s, work / f)
    for deck, (old, new) in spec.get("patch", {}).items():
        p = work / deck
        p.write_text(p.read_text().replace(old, new))
    if spec.get("laser_fix"):
        # R2-3-2: 仓库 laser.dat = DESY 原版 (MD5 68d016175859b20c1e2ccee5057c2d46,
        # 浮点计数头是原版特征); macOS Apple Silicon 构建直接运行会报
        # "Error while reading file: laser.dat"。此处仅把头 3 行计数取整
        # (数据体不动), 整数头重跑可字节复现 golden。
        fix_laser_map_header(work / "laser.dat")
        print("  laser.dat 图头计数已转为整数形式")
    return work, spec


def run_example(name):
    """备文件 + 跑完 spec 里全部步骤 (generator/astra)。"""
    work, spec = stage_files(name)
    for kind, deck in spec["steps"]:
        exe = check(kind)
        run_program(exe, work, input_file=deck)
    return work


def check(kind):
    """定位可执行文件 (PATH -> 项目 ASTRA/ 目录)。"""
    from astra_tools.run import check_executable
    exe = check_executable(kind, project_dir=PROJECT_ROOT)
    return exe


def phase_files(work, stem):
    """工作目录里按 z 排序的相空间文件列表。"""
    return sorted(p for p in work.glob(stem + ".*.001")
                  if p.name.split(".")[1].lstrip("-").isdigit())


def _deck_zstop(work, deck):
    """从输入卡解析 OUTPUT/ZSTOP [m]; 解析不到返回 None。

    compare_xemit 护栏用: 断言新运行末行 z 到达 ZSTOP, 防止
    空跑算例 (参考粒子落后于 ZSTOP -> 0 迭代) 生成空洞 golden。
    """
    from astra_tools.namelist.parse import parse_namelists
    try:
        z = parse_namelists(str(work / deck))["OUTPUT"]["ZSTOP"]
    except Exception:
        return None
    arr = np.asarray(z, dtype=float).ravel()
    return None if arr.size == 0 else float(arr[0])


def compare_xemit(name, work):
    """新运行 vs 归档 golden 的末行比对 (rel < 0.5% 判 OK)。

    护栏 (二轮审计 R2-3-1): 空洞比对防护。Xemit 必须多行
    (> 1 行, 新运行与归档 golden 都检查) 且新运行末行 z 到达
    输入卡 ZSTOP (相对容差 0.5%), 否则判不通过 — 旧算例 0 迭代
    时 Xemit 仅 1 行且数值 = 输入统计原样, rel=0.0000% 恒真。

    Returns:
        True 全部通过; False 任一护栏或末行比对不通过。
    """
    spec = EXAMPLES[name]
    golden = spec["golden_xemit"]
    new_file = work / golden.name
    if not new_file.exists():
        print("  (无 %s 输出, 跳过比对)" % golden.name)
        return False
    if spec.get("compare_mode") == "log":
        ok_new = "finished" in new_file.read_text()
        ok_ref = "finished" in golden.read_text()
        print("  末行比对 (new vs golden, log):")
        print("    finished        %-10s %-10s %s"
              % (ok_new, ok_ref, "OK" if ok_new and ok_ref else "MISMATCH"))
        return ok_new and ok_ref
    new = parse_output_file(new_file)
    ref = parse_output_file(golden)
    ok = True

    # ---- 护栏 1: Xemit 行数 > 1 (新与 golden 都检查) ----
    n_new = len(np.asarray(new["mean_z"]))
    n_ref = len(np.asarray(ref["mean_z"]))
    print("  护栏: Xemit 行数 new=%d golden=%d" % (n_new, n_ref))
    if n_new <= 1 or n_ref <= 1:
        print("  护栏失败: Xemit 行数必须 > 1 (空跑算例仅 1 行,"
              " 数值=输入统计原样, 比对无意义)")
        ok = False

    # ---- 护栏 2: 新运行末行 z 到达输入卡 ZSTOP ----
    zstop = _deck_zstop(work, spec["steps"][-1][1])
    if zstop is not None:
        z_last = float(np.asarray(new["mean_z"])[-1])
        rel_z = abs(z_last - zstop) / zstop * 100
        print("  护栏: 末行 z=%.6g vs ZSTOP=%.6g (rel=%.4f%%)"
              % (z_last, zstop, rel_z))
        if rel_z > 0.5:
            print("  护栏失败: 末行 z 未到达 ZSTOP"
                  " (束团可能在孔径/堵块处提前全部损失)")
            ok = False

    print("  末行比对 (new vs golden):")
    for key in ("norm_emit_x", "sigma_x", "mean_z"):
        a = float(np.asarray(new[key])[-1])
        b = float(np.asarray(ref[key])[-1])
        rel = abs(a - b) / abs(b) * 100 if b else float("inf")
        passed = rel < 0.5
        ok = ok and passed
        print("    %-14s %-10.6g %-10.6g rel=%.4f%% %s"
              % (key, a, b, rel, "OK" if passed else "MISMATCH"))
    print("  总体: %s" % ("OK" if ok else "不通过 (见上方护栏/比对标记)"))
    return ok
