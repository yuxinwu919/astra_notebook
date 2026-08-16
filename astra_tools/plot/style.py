"""Matplotlib style presets for modern, publication-quality figures.

No import-time side effects; call set_style() once per notebook.

字体策略 (2026-08): 正文统一 Times New Roman (用户要求), 按
font.family 列表逐字形回退: 缺字形时依次尝试 STIXGeneral (数学
符号) -> DejaVu Serif -> 各中文宋/黑体, 保证 π/′/中文等符号在任何
平台都不出现方框 (matplotlib >= 3.6 的字体回退机制)。数学模式
(mathtext) 用 'stix' 字库, 与 Times 风格一致。
"""

from __future__ import annotations

from matplotlib import font_manager, rcParams

# 正文字体回退链 (候选池): 拉丁 -> 数学符号 -> 通用衬线 -> 中文。
# 运行时只保留本机真实存在的字体: 若把不存在的家族留在
# font.family 里, matplotlib 对每个文字对象都会为每个缺失家族
# 发一条 findfont 警告 —— 一张图数千条警告, Jupyter 前端渲染
# 这些消息流会把 notebook 拖到"几分钟跑不完"。
FONT_FAMILY = [
    "Times New Roman",
    "STIXGeneral",
    "DejaVu Serif",
    "Liberation Serif",
    "Songti SC",      # macOS 宋体
    "STSong",
    "SimSun",         # Windows 宋体
    "Noto Serif CJK SC",
    "WenQuanYi Zen Hei",
    "Microsoft YaHei",
    "PingFang SC",
    "Heiti SC",
    "Arial Unicode MS",
    "DejaVu Sans",    # 兜底: 覆盖面最广
]

_AVAILABLE_FONTS: list | None = None


def available_font_family() -> list:
    """FONT_FAMILY 中本机真实存在的字体 (按原顺序), 末尾兜底 DejaVu Sans.

    matplotlib 自带的 DejaVu Sans 必然存在; 扫描结果缓存一次。
    """
    global _AVAILABLE_FONTS
    if _AVAILABLE_FONTS is None:
        installed = {f.name for f in font_manager.fontManager.ttflist}
        keep = [f for f in FONT_FAMILY if f in installed]
        if "DejaVu Sans" not in keep:
            keep.append("DejaVu Sans")
        _AVAILABLE_FONTS = keep
    return _AVAILABLE_FONTS


# Colorblind-friendly qualitative palette (Paul Tol)
COLORS = ["#0077BB", "#EE7733", "#009988", "#CC3311", "#33BBEE", "#EE3377"]

# Density colormap (default: viridis; alternative SLAC-DESY beam map)
DEFAULT_CMAP = "viridis"


def set_style(
    font_size: int = 12,
    dpi: int = 120,
    fig_width_inches: float = 6.0,
    fig_height_inches: float = 4.5,
) -> None:
    """Apply the astra-notebook figure theme."""
    rcParams.update({
        "figure.dpi": dpi,
        "savefig.dpi": 150,
        "figure.figsize": (fig_width_inches, fig_height_inches),
        "savefig.bbox": "tight",
        "font.size": font_size,
        # 统一 Times New Roman; 列表即 matplotlib 的字形回退链,
        # 缺字形 (π/′/中文等) 自动落到后续字体, 不再出现方框。
        # 只含本机存在的字体, 避免缺失家族引发 findfont 警告洪流
        # (数千条 IOPub 消息会把 Jupyter 单元拖到几分钟不结束)。
        "font.family": available_font_family(),
        "font.serif": available_font_family(),
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        # mathtext 用 Times 风格的 STIX 字库 ($\beta$ 等数学符号)
        "mathtext.fontset": "stix",
        "axes.labelsize": font_size + 1,
        "axes.titlesize": font_size + 2,
        "axes.linewidth": 1.0,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
        "axes.unicode_minus": False,
        "axes.prop_cycle": __import__("matplotlib").rcsetup.cycler("color", COLORS),
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.fontsize": font_size - 1,
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.edgecolor": "0.8",
        "lines.linewidth": 1.6,
        "lines.markersize": 5,
        "image.cmap": DEFAULT_CMAP,
        "text.usetex": False,
    })


def reset_style() -> None:
    """Reset rcParams to matplotlib defaults."""
    import matplotlib.pyplot as plt

    rcParams.update(plt.rcParamsDefault)
