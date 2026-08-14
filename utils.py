# ============================================================
# utils.py — ASTRA / Generator 共享工具模块
# ============================================================
# 本模块提供 Namelist 生成、二进制文件解析、外部程序调用、
# 束团统计计算和相空间绘图等通用功能。
# 供 generator_interface.ipynb 和 astra_interface.ipynb 导入使用。
# ============================================================

import os
import sys
import subprocess
import shutil
import logging
from pathlib import Path
from typing import Any, Callable, Optional, Tuple, Union
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, constants

# ============================================================
# 全局配置
# ============================================================

# Matplotlib 中文字体设置
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'Heiti SC', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 150

# Logging 配置
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# ============================================================
# 解析器检测
# ============================================================

HAS_LUME_ASTRA = False
HAS_OCELOT = False

try:
    import ocelot
    HAS_OCELOT = True
    logger.info("✓ ocelot 已加载（首选解析器）")
except ImportError:
    logger.warning("✗ ocelot 未安装。将尝试 lume-astra。")
    logger.warning("  安装方法: pip install ocelot")

try:
    import lume_astra
    HAS_LUME_ASTRA = True
    logger.info("✓ lume-astra 已加载（备选解析器）")
except ImportError:
    logger.warning("✗ lume-astra 未安装。将使用自定义二进制解析器。")
    logger.warning("  安装方法: pip install lume-astra")

logger.info("所有基础库导入完成。")


# ============================================================
# 可执行文件检测
# ============================================================

def check_executable(name: str, fallback_dir: Optional[Path] = None) -> str:
    """检查可执行文件是否可用，返回完整路径或抛出异常。

    查找顺序：1) 系统 PATH → 2) fallback_dir/astra/ 子目录

    Args:
        name: 可执行文件名称（如 'astra', 'generator'）
        fallback_dir: 项目根目录（通常为 PROJECT_DIR）

    Returns:
        可执行文件的完整路径

    Raises:
        FileNotFoundError: 若在所有位置都未找到
    """
    # 1. 查系统 PATH
    exe_path = shutil.which(name)
    if exe_path is not None:
        logger.info(f"✓ 在 PATH 中找到 {name}: {exe_path}")
        return exe_path

    # 2. 查项目 astra/ 子目录
    if fallback_dir is not None:
        local_path = fallback_dir / "astra" / name
        if local_path.is_file() and os.access(local_path, os.X_OK):
            logger.info(f"✓ 在项目目录中找到 {name}: {local_path}")
            return str(local_path)

    raise FileNotFoundError(
        f"未找到可执行文件 '{name}'。请通过以下方式之一获取：\n"
        f"  1) 将 {name} 安装到系统 PATH，或\n"
        f"  2) 将 {name} 可执行文件放入项目 ({fallback_dir / 'astra'}) 目录，或\n"
        f"  3) 从 DESY 下载 ASTRA 程序包: https://www.desy.de/~mpyflo/"
    )


def get_version(exe_path: str, work_dir: Path, timeout: int = 5) -> str:
    """运行可执行文件获取版本信息（从前几行输出中提取）。

    Args:
        exe_path: 可执行文件路径
        work_dir: 工作目录
        timeout: 超时时间（秒）

    Returns:
        版本字符串，或 "未知版本" / "无法检测"
    """
    try:
        result = subprocess.run(
            [exe_path],
            capture_output=True, text=True, timeout=timeout,
            cwd=work_dir,
            input="",
        )
        for line in result.stdout.split("\n"):
            if "version" in line.lower():
                return line.strip()
        return "未知版本"
    except (subprocess.TimeoutExpired, OSError, ValueError) as e:
        logger.warning(f"版本检测失败: {e}")
        return "无法检测"


# ============================================================
# Namelist 文件生成
# ============================================================

def write_namelist(
    namelist_name: str,
    params: dict[str, Any],
    filepath: Optional[Path] = None,
) -> Optional[str]:
    """将 Python 字典转换为 ASTRA/Generator 风格的 namelist 格式。

    处理规则（适配 ASTRA 自定义解析器）：
    - str 值：直接输出（用户需自行添加引号，如 "'bunch.ini'"）
    - bool 值：转换为 T / F（ASTRA 风格）
    - int / float：直接输出
    - list / tuple：转换为逗号分隔的值列表

    Args:
        namelist_name: Namelist 名称（如 'NEWRUN', 'INPUT', 'OUTPUT'）
        params: 参数字典
        filepath: 输出文件路径。若为 None，则返回格式化字符串而不写文件。

    Returns:
        若 filepath 为 None，返回格式化后的字符串；否则返回 None。
    """
    lines = [f"&{namelist_name}"]

    for key, value in params.items():
        # Skip empty lists, empty strings, and None values (optional parameters)
        if value is None:
            continue
        if isinstance(value, str) and value == "":
            continue
        if isinstance(value, (list, tuple)) and len(value) == 0:
            continue

        if isinstance(value, bool):
            formatted = "T" if value else "F"
        elif isinstance(value, str):
            formatted = value
        elif isinstance(value, (list, tuple, np.ndarray)):
            arr = np.asarray(value).flatten()
            str_values = []
            for v in arr:
                if isinstance(v, (bool, np.bool_)):
                    str_values.append("T" if v else "F")
                elif isinstance(v, str):
                    str_values.append(v)
                elif isinstance(v, int):
                    str_values.append(str(v))
                else:
                    str_values.append(f"{v:.12g}")
            formatted = ", ".join(str_values)
        elif isinstance(value, int):
            formatted = str(value)
        elif isinstance(value, float):
            formatted = f"{value:.12g}"
        else:
            formatted = str(value)

        lines.append(f"  {key}={formatted},")

    lines.append(" /")

    result = "\n".join(lines) + "\n"

    if filepath is None:
        return result

    with open(filepath, "w") as f:
        f.write(result)

    logger.info(f"Namelist 文件已写入: {filepath}")
    return None


# ============================================================
# 二进制粒子分布文件读取
# ============================================================

def read_astra_binary(filepath: Path) -> Tuple[np.ndarray, np.ndarray]:
    """读取 ASTRA / Generator 格式的二进制粒子分布文件。

    文件格式（标准 ASTRA 分布格式，参考手册 Table 1）：
        - 前 5 个 float64：文件头信息
          [0]: 参考时间 (ns)
          [1]: 参考动量/能量 (eV/c)
          [2]: 束团总电荷 (nC)
          [3]: 保留位
          [4]: 保留位
        - 后续每 9 个 float64 描述一个粒子（二进制格式）：
          [0] x  [m]       横向位置
          [1] y  [m]       垂直位置
          [2] z  [m]       纵向位置（相对于参考粒子）
          [3] px [eV/c]    横向动量
          [4] py [eV/c]    垂直动量
          [5] pz [eV/c]    纵向动量
          [6] clock [ns]   到达时间（注意：单位是 ns！）
          [7] macro_charge [nC]  宏粒子电荷
          [8] status_flag        状态标志（0=活跃, >0=丢失, <0=未开始/丢失）

    Args:
        filepath: 二进制文件路径

    Returns:
        header: np.ndarray, shape (5,), float64
        particles: np.ndarray, shape (N, 9), float64
    """
    data = np.fromfile(filepath, dtype=np.float64)

    if len(data) < 5:
        raise ValueError(f"文件 {filepath} 数据不足：仅读取到 {len(data)} 个值")

    header = data[:5].copy()

    # 尝试 9 列或 10 列格式
    n_particles_9 = (len(data) - 5) // 9
    n_particles_10 = (len(data) - 5) // 10

    remainder_9 = (len(data) - 5) % 9
    remainder_10 = (len(data) - 5) % 10

    if remainder_9 == 0 and remainder_10 == 0:
        logger.info("检测到可能包含 particle_index 的 10 列格式，尝试按 9 列读取...")
        n_particles = n_particles_9
    elif remainder_9 == 0:
        n_particles = n_particles_9
    elif remainder_10 == 0:
        n_particles = n_particles_10
        logger.warning(
            f"检测到 10 列格式（含 particle_index），"
            f"将忽略第 10 列。N={n_particles}"
        )
    else:
        logger.warning(
            f"文件大小与预期不符：剩余 {len(data) - 5} 个值，"
            f"9列余数={remainder_9}, 10列余数={remainder_10}。"
            f"尝试按 {n_particles_9} 个粒子（9列）读取。"
        )
        n_particles = n_particles_9

    particles = data[5:5 + n_particles * 9].reshape(n_particles, 9)

    logger.info(
        f"从 {filepath.name} 读取："
        f"ref_time={header[0]:.4e} ns, ref_energy={header[1]:.4e} eV, "
        f"Q_total={header[2]:.4f} nC, N_particles={n_particles}"
    )

    return header, particles


def _array_to_dict(header: np.ndarray, particles: np.ndarray, filepath: Path) -> dict[str, Any]:
    """将 numpy 数组转换为标准化的字典格式。"""
    n_particle = len(particles)
    active = particles[:, 8] == 0

    return {
        'filepath': filepath,
        'header': header,
        'x': particles[:, 0],
        'y': particles[:, 1],
        'z': particles[:, 2],
        'px': particles[:, 3],
        'py': particles[:, 4],
        'pz': particles[:, 5],
        'clock': particles[:, 6],
        'macro_charge': particles[:, 7],
        'status_flag': particles[:, 8],
        'n_particle': n_particle,
        'n_active': int(np.sum(active)),
        'total_charge': float(np.sum(particles[active, 7])),
        'active_mask': active,
    }


def _normalize_distribution(data: Any, filepath: Path) -> dict[str, Any]:
    """将 lume-astra 的输出标准化为统一字典格式。"""
    if isinstance(data, dict):
        result = {'filepath': filepath}
        key_map = {
            'x': 'x', 'y': 'y', 'z': 'z',
            'px': 'px', 'py': 'py', 'pz': 'pz',
            't': 'clock', 'clock': 'clock',
            'macro_charge': 'macro_charge', 'charge': 'macro_charge',
            'status': 'status_flag', 'flag': 'status_flag',
        }
        for src, dst in key_map.items():
            if src in data:
                result[dst] = np.asarray(data[src])
        result['n_particle'] = len(result.get('x', []))
        result['total_charge'] = float(np.sum(result.get('macro_charge', [0])))
        result['active_mask'] = np.ones(result['n_particle'], dtype=bool)
        result['header'] = np.zeros(5)
        return result
    else:
        return _array_to_dict(np.zeros(5), np.asarray(data), filepath)


def _ocelot_to_dict(p_array: Any, filepath: Path) -> dict[str, Any]:
    """将 ocelot ParticleArray 转换为标准化字典格式。

    ocelot ParticleArray 的坐标约定：
    - .x(), .y(): 横向位置 [m]
    - .tau(): 纵向位置 τ = z - ct [m]
    - .px(), .py(): 横向散角 x'=px/p0, y'=py/p0 [rad]（无量纲）
    - .p(): 相对能量偏离 δ = (p-p0)/p0 ≈ (E-E0)/E0
    - .E: 参考能量 [GeV]
    - .q_array: 宏粒子电荷 [C]
    """
    n = p_array.n
    p_ref = p_array.E * 1e9  # GeV → eV/c 近似

    xp = p_array.px()
    yp = p_array.py()
    delta = p_array.p()
    pz = p_ref * (1.0 + delta)
    px = xp * pz
    py = yp * pz

    q_nC = p_array.q_array * 1e9

    return {
        'filepath': filepath,
        'header': np.array([0.0, p_ref, q_nC.sum(), 0.0, 0.0]),
        'x': p_array.x(),
        'y': p_array.y(),
        'z': p_array.tau(),
        'px': px,
        'py': py,
        'pz': pz,
        'clock': np.zeros(n),
        'macro_charge': q_nC,
        'status_flag': np.zeros(n),
        'n_particle': n,
        'n_active': n,
        'total_charge': float(q_nC.sum()),
        'active_mask': np.ones(n, dtype=bool),
    }


def read_astra_distribution(filepath: Path) -> dict[str, Any]:
    """读取粒子分布文件（自动选择最优解析器）。

    优先使用 ocelot，其次 lume-astra，最后使用自定义二进制解析。

    Args:
        filepath: 分布文件路径

    Returns:
        dict，包含以下键：
        - 'header': np.ndarray (5,)
        - 'x', 'y', 'z': 位置数组 [m]
        - 'px', 'py', 'pz': 动量数组 [eV/c]
        - 'clock': 时间数组 [s]
        - 'macro_charge': 宏粒子电荷 [nC]
        - 'status_flag': 状态标志
        - 'n_particle': 粒子总数
        - 'total_charge': 总电荷 [nC]
    """
    if HAS_OCELOT:
        try:
            logger.info("使用 ocelot 读取分布文件...")
            from ocelot.adaptors.astra2ocelot import astraBeam2particleArray
            p_array = astraBeam2particleArray(str(filepath))
            logger.info("ocelot 读取成功。")
            return _ocelot_to_dict(p_array, filepath)
        except Exception as e:
            logger.warning(f"ocelot 读取失败: {e}，回退到 lume-astra。")

    if HAS_LUME_ASTRA:
        try:
            logger.info("使用 lume-astra 读取分布文件...")
            data_dict = lume_astra.read_astra(str(filepath))
            logger.info("lume-astra 读取成功。")
            return _normalize_distribution(data_dict, filepath)
        except Exception as e:
            logger.warning(f"lume-astra 读取失败: {e}，回退到自定义解析器。")

    # Fallback: 自定义二进制解析
    logger.info("使用自定义二进制解析器读取分布文件...")
    header, particles = read_astra_binary(filepath)
    return _array_to_dict(header, particles, filepath)


# ============================================================
# 外部程序调用
# ============================================================

def run_exe(
    exe_path: str,
    work_dir: Path,
    input_file: Optional[str] = None,
    timeout: int = 300,
) -> subprocess.CompletedProcess:
    """在指定工作目录中运行可执行文件，捕获并打印输出。

    Args:
        exe_path: 可执行文件路径
        work_dir: 工作目录
        input_file: 可选，输入文件路径（作为命令行参数传递）
        timeout: 超时时间（秒）

    Returns:
        subprocess.CompletedProcess 对象

    Raises:
        RuntimeError: 若程序返回码非零
    """
    exe_name = Path(exe_path).name
    logger.info(f"运行 {exe_name} (工作目录: {work_dir})...")

    try:
        if input_file is not None:
            result = subprocess.run(
                [exe_path, input_file],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        else:
            result = subprocess.run(
                [exe_path],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=timeout,
            )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"{exe_name} 执行超时（>{timeout}s）")

    if result.stdout:
        print(result.stdout)

    if result.returncode != 0:
        error_msg = f"{exe_name} 返回非零退出码 {result.returncode}"
        if result.stderr:
            error_msg += f"\nstderr:\n{result.stderr}"
        raise RuntimeError(error_msg)

    logger.info(f"{exe_name} 执行成功。（返回码: {result.returncode}）")
    return result


# ============================================================
# 束团统计计算
# ============================================================

def compute_beam_statistics(dist: dict[str, Any], ref_energy_eV: Optional[float] = None) -> dict[str, Any]:
    """计算束团的关键统计量。

    Args:
        dist: read_astra_distribution 返回的字典
        ref_energy_eV: 参考能量 [eV]，用于计算归一化发射度。
                       若为 None，则从 header[1] 或 pz 均值估计。

    Returns:
        dict，包含束团统计量
    """
    mask = dist['active_mask']
    x, y, z = dist['x'][mask], dist['y'][mask], dist['z'][mask]
    px, py, pz = dist['px'][mask], dist['py'][mask], dist['pz'][mask]

    if ref_energy_eV is None:
        ref_energy_eV = dist['header'][1] if dist['header'][1] != 0 else np.mean(pz)

    # ---- 位置统计 ----
    mean_x, mean_y, mean_z = np.mean(x), np.mean(y), np.mean(z)
    sig_x = np.std(x - mean_x)
    sig_y = np.std(y - mean_y)
    sig_z = np.std(z - mean_z)

    # ---- 动量统计 ----
    mean_px, mean_py, mean_pz = np.mean(px), np.mean(py), np.mean(pz)
    sig_px = np.std(px - mean_px)
    sig_py = np.std(py - mean_py)
    sig_pz = np.std(pz - mean_pz)

    sig_E_over_E = sig_pz / mean_pz if mean_pz != 0 else 0.0

    # ---- 横向散角 ----
    xp = (px - mean_px) / np.abs(pz)
    yp = (py - mean_py) / np.abs(pz)
    x_centered = x - mean_x
    y_centered = y - mean_y

    # ---- 几何发射度 (RMS) ----
    x2 = np.mean(x_centered ** 2)
    xp2 = np.mean(xp ** 2)
    xxp = np.mean(x_centered * xp)
    emit_x_geom = np.sqrt(max(x2 * xp2 - xxp ** 2, 0))

    y2 = np.mean(y_centered ** 2)
    yp2 = np.mean(yp ** 2)
    yyp = np.mean(y_centered * yp)
    emit_y_geom = np.sqrt(max(y2 * yp2 - yyp ** 2, 0))

    # ---- 归一化发射度 ----
    m_e_c2_eV = constants.m_e * constants.c**2 / constants.e
    gamma = ref_energy_eV / m_e_c2_eV
    beta_rel = np.sqrt(1 - 1 / gamma**2) if gamma > 1 else 0.0
    bg = beta_rel * gamma

    emit_x_norm = bg * emit_x_geom if bg > 0 else 0.0
    emit_y_norm = bg * emit_y_geom if bg > 0 else 0.0

    # ---- Twiss 参数 ----
    beta_x = x2 / emit_x_geom if emit_x_geom > 0 else 0.0
    alpha_x = -xxp / emit_x_geom if emit_x_geom > 0 else 0.0
    beta_y = y2 / emit_y_geom if emit_y_geom > 0 else 0.0
    alpha_y = -yyp / emit_y_geom if emit_y_geom > 0 else 0.0

    return {
        'n_particle': dist['n_particle'],
        'n_active': int(dist['n_active']),
        'total_charge_nC': dist['total_charge'],
        'mean_x': mean_x, 'mean_y': mean_y, 'mean_z': mean_z,
        'sig_x': sig_x, 'sig_y': sig_y, 'sig_z': sig_z,
        'mean_px': mean_px, 'mean_py': mean_py, 'mean_pz': mean_pz,
        'sig_px': sig_px, 'sig_py': sig_py, 'sig_pz': sig_pz,
        'sig_E_over_E': sig_E_over_E,
        'emit_x_geom': emit_x_geom, 'emit_y_geom': emit_y_geom,
        'emit_x_norm': emit_x_norm, 'emit_y_norm': emit_y_norm,
        'beta_x': beta_x, 'alpha_x': alpha_x,
        'beta_y': beta_y, 'alpha_y': alpha_y,
        'ref_energy_eV': ref_energy_eV,
        'gamma': gamma,
        'beta_rel': beta_rel,
    }


def print_statistics(stats: dict[str, Any], title: str = "束团统计") -> None:
    """格式化打印统计结果。"""
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")
    print(f"  粒子总数:        {stats['n_particle']:>10d}")
    print(f"  活跃粒子数:      {stats['n_active']:>10d}")
    print(f"  总电荷:          {stats['total_charge_nC']:>10.4f} nC")
    print(f"{'-'*40}")
    print(f"  参考能量:        {stats['ref_energy_eV']/1e6:>10.4f} MeV")
    print(f"  相对论 γ:        {stats['gamma']:>10.4f}")
    print(f"  相对论 β:        {stats['beta_rel']:>10.6f}")
    print(f"{'-'*40}")
    print(f"  质心 x:          {stats['mean_x']*1e3:>10.4f} mm")
    print(f"  质心 y:          {stats['mean_y']*1e3:>10.4f} mm")
    print(f"  质心 z:          {stats['mean_z']*1e3:>10.4f} mm")
    print(f"  σ_x:             {stats['sig_x']*1e3:>10.4f} mm")
    print(f"  σ_y:             {stats['sig_y']*1e3:>10.4f} mm")
    print(f"  σ_z:             {stats['sig_z']*1e3:>10.4f} mm")
    print(f"{'-'*40}")
    print(f"  平均 pz:         {stats['mean_pz']/1e6:>10.4f} MeV/c")
    print(f"  σ_pz:            {stats['sig_pz']/1e6:>10.4f} MeV/c")
    print(f"  σ_E / E:         {stats['sig_E_over_E']:>10.6f}")
    print(f"{'-'*40}")
    print(f"  几何发射度 ε_x:  {stats['emit_x_geom']*1e6:>10.4f} μm·rad")
    print(f"  几何发射度 ε_y:  {stats['emit_y_geom']*1e6:>10.4f} μm·rad")
    print(f"  归一化发射度 ε_nx:{stats['emit_x_norm']*1e6:>10.4f} μm·rad")
    print(f"  归一化发射度 ε_ny:{stats['emit_y_norm']*1e6:>10.4f} μm·rad")
    print(f"{'='*60}\n")


# ============================================================
# 绘图函数
# ============================================================

def plot_transverse_phase_space(
    dist: dict[str, Any],
    title_prefix: str = "",
    figsize: Tuple[int, int] = (12, 5),
) -> plt.Figure:
    """绘制横向相空间散点图 (x-x', y-y')。"""
    mask = dist['active_mask']
    x = dist['x'][mask] * 1e3
    y = dist['y'][mask] * 1e3
    pz = np.abs(dist['pz'][mask])
    px = dist['px'][mask]
    py = dist['py'][mask]

    xp = px / pz * 1e3
    yp = py / pz * 1e3

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    ax = axes[0]
    counts, xedges, yedges, im = ax.hist2d(
        x, xp, bins=80, cmap='viridis',
        range=[[x.min(), x.max()], [xp.min(), xp.max()]]
    )
    plt.colorbar(im, ax=ax, label='Counts')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel("x' [mrad]")
    ax.set_title(f'{title_prefix} x–x\' Phase Space')
    ax.axhline(0, color='white', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='white', linewidth=0.5, linestyle='--')

    ax = axes[1]
    counts, xedges, yedges, im = ax.hist2d(
        y, yp, bins=80, cmap='viridis',
        range=[[y.min(), y.max()], [yp.min(), yp.max()]]
    )
    plt.colorbar(im, ax=ax, label='Counts')
    ax.set_xlabel('y [mm]')
    ax.set_ylabel("y' [mrad]")
    ax.set_title(f'{title_prefix} y–y\' Phase Space')
    ax.axhline(0, color='white', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='white', linewidth=0.5, linestyle='--')

    plt.tight_layout()
    return fig


def plot_longitudinal_phase_space(
    dist: dict[str, Any],
    title_prefix: str = "",
    figsize: Tuple[int, int] = (10, 5),
) -> plt.Figure:
    """绘制纵向相空间散点图 (z vs δp/p)。"""
    mask = dist['active_mask']
    z = dist['z'][mask] * 1e3
    pz = dist['pz'][mask]
    mean_pz = np.mean(pz)
    delta_p = (pz - mean_pz) / mean_pz * 100

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    counts, xedges, yedges, im = ax.hist2d(z, delta_p, bins=100, cmap='viridis')
    plt.colorbar(im, ax=ax, label='Counts')
    ax.set_xlabel('z [mm]（粒子前进方向 +z）')
    ax.set_ylabel('δp/p [%]')
    ax.set_title(f'{title_prefix} Longitudinal Phase Space (z vs δp/p)')
    ax.axhline(0, color='white', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='white', linewidth=0.5, linestyle='--')

    plt.tight_layout()
    return fig


def plot_spatial_distributions(
    dist: dict[str, Any],
    title_prefix: str = "",
    figsize: Tuple[int, int] = (14, 4),
) -> plt.Figure:
    """绘制 x, y, z 的投影分布直方图。"""
    mask = dist['active_mask']
    x = dist['x'][mask] * 1e3
    y = dist['y'][mask] * 1e3
    z = dist['z'][mask] * 1e3

    fig, axes = plt.subplots(1, 3, figsize=figsize)

    for ax, data, label in zip(axes, [x, y, z], ['x', 'y', 'z']):
        ax.hist(data, bins=80, density=True, alpha=0.7,
                color='steelblue', edgecolor='white')
        mu, sigma = np.mean(data), np.std(data)
        x_fit = np.linspace(data.min(), data.max(), 200)
        ax.plot(x_fit, stats.norm.pdf(x_fit, mu, sigma),
                'r-', linewidth=2, label=f'σ={sigma:.3f} mm')
        ax.set_xlabel(f'{label} [mm]')
        ax.set_ylabel('Probability Density')
        ax.set_title(f'{title_prefix} {label} Distribution')
        ax.legend()

    plt.tight_layout()
    return fig


def plot_transverse_scatter(
    dist: dict[str, Any],
    title_prefix: str = "",
    figsize: Tuple[int, int] = (7, 6),
) -> plt.Figure:
    """绘制 x-y 横向散点图。"""
    mask = dist['active_mask']
    x = dist['x'][mask] * 1e3
    y = dist['y'][mask] * 1e3

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    counts, xedges, yedges, im = ax.hist2d(x, y, bins=80, cmap='inferno')
    plt.colorbar(im, ax=ax, label='Counts')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_title(f'{title_prefix} Transverse Beam Profile (x–y)')
    ax.set_aspect('equal')
    ax.axhline(0, color='white', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='white', linewidth=0.5, linestyle='--')

    plt.tight_layout()
    return fig


# ============================================================
# 环境信息展示
# ============================================================

def _display_width(s: str) -> int:
    """计算字符串的显示宽度（中文字符计为 2 宽度）。"""
    w = 0
    for ch in s:
        w += 2 if '\u4e00' <= ch <= '\u9fff' or '\u3000' <= ch <= '\u303f' or '\uff00' <= ch <= '\uffef' else 1
    return w


def print_environment_header(
    tool_name: str,
    tool_version: str,
    sim_dir: Path,
    os_version: str,
) -> None:
    """打印统一的环境初始化横幅。

    Args:
        tool_name: 工具名称（如 'Generator', 'ASTRA'）
        tool_version: 版本字符串
        sim_dir: 工作目录
        os_version: 操作系统版本字符串
    """
    from datetime import datetime

    def _p(label, value, target_width=20):
        pad = target_width - _display_width(label)
        print(f"  {label}{' ' * max(pad, 1)} {value}")

    preferred = 'ocelot' if HAS_OCELOT else 'lume-astra' if HAS_LUME_ASTRA else '自定义 np.fromfile'

    print("\n" + "=" * 60)
    print(f"  {tool_name} 初始化完成 — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    _p("操作系统",     os_version)
    _p(tool_name,   tool_version)
    _p("工作目录",     sim_dir)
    _p("首选解析器",   preferred)
    print("=" * 60 + "\n")
