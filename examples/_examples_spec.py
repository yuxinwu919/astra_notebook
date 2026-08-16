"""官方算例共享规格: 文件清单/运行步骤/黄金比对目标。

供 examples/<name>/<name>.ipynb 单算例示例 notebook 共用 - 单一
数据源, 避免 9 份拷贝漂移 (批 5: 06_examples 汇总本已删除)。
"""

import json
import shutil
from pathlib import Path

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
    """复制官方算例输入到工作目录 (不运行), 返回 (work, spec)。"""
    spec = EXAMPLES[name]
    src = EXAMPLES_DIR / spec.get("src_dir", name)
    work = SIM_DIR / name
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


def compare_xemit(name, work):
    """新运行 vs 归档 golden 的末行比对 (rel < 0.5% 判 OK)。"""
    spec = EXAMPLES[name]
    golden = spec["golden_xemit"]
    new_file = work / golden.name
    if not new_file.exists():
        print("  (无 %s 输出, 跳过比对)" % golden.name)
        return
    if spec.get("compare_mode") == "log":
        ok_new = "finished" in new_file.read_text()
        ok_ref = "finished" in golden.read_text()
        print("  末行比对 (new vs golden, log):")
        print("    finished        %-10s %-10s %s"
              % (ok_new, ok_ref, "OK" if ok_new and ok_ref else "MISMATCH"))
        return
    new = parse_output_file(new_file)
    ref = parse_output_file(golden)
    print("  末行比对 (new vs golden):")
    for key in ("norm_emit_x", "sigma_x", "mean_z"):
        a = float(__import__("numpy").asarray(new[key])[-1])
        b = float(__import__("numpy").asarray(ref[key])[-1])
        rel = abs(a - b) / abs(b) * 100
        print("    %-14s %-10.6g %-10.6g rel=%.4f%% %s"
              % (key, a, b, rel, "OK" if rel < 0.5 else "MISMATCH"))
