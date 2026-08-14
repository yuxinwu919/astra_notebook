"""Jupyter 前端组件: 选择器 / 参数表单 / 统计面板 / 运行面板.

需要 ipywidgets; 无 Jupyter 环境下可导入但组件不显示。
"""

from .selectors import discover_sim_runs, run_selector, phase_selector
from .forms import namelist_form, form_values
from .panels import stats_table_html, distribution_summary_html

__all__ = [
    "discover_sim_runs", "run_selector", "phase_selector",
    "namelist_form", "form_values",
    "stats_table_html", "distribution_summary_html",
]
