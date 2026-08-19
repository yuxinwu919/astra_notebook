"""HTML 统计面板 (Jupyter 富呈现).
"""

from __future__ import annotations

from IPython.display import HTML


def stats_table_html(stats) -> HTML:
    """束团统计两列表格 (量/值/单位)."""
    d = stats.to_dict()
    rows = []
    for k, v in d.items():
        unit = _UNITS.get(k, "")
        if k == "total_charge_nC":
            v = "%.6g (|Q|)" % abs(v)
        elif isinstance(v, float):
            v = "%.6g" % v
        rows.append("<tr><td>%s</td><td align='right'>%s</td><td>%s</td></tr>"
                    % (k, v, unit))
    body = "\n".join(rows)
    return HTML(
        "<table style='border-collapse:collapse;font-family:monospace'>"
        "<tr><th>quantity</th><th>value</th><th>unit</th></tr>" + body + "</table>"
    )


def distribution_summary_html(dist) -> HTML:
    """分布概要面板 (批 5: HTML 转义 source 路径)。"""
    import html
    return HTML("<pre>" + html.escape(dist.summary()) + "</pre>")


def display_bz_warning(sim_dir) -> bool:
    """deck 含螺线管且统计未传 bz 时的黄色告警 (批 2, P0-4 兜底).

    从 astra.in 解析 &SOLENOID: LBField=T 即认为束线含螺线管场;
    此时默认 compute_statistics(bz=0) 的结果与 ASTRA Xemit 口径
    不一致 (正则动量, 手册 4.13.1)。不自动改数值, 只提醒。
    返回是否有告警; 找不到 deck 时静默返回 False。
    """
    from pathlib import Path
    from IPython.display import display
    from ..namelist.parse import get_ci, iter_namelist_blocks, parse_namelists
    deck = Path(sim_dir) / "astra.in"
    if not deck.exists():
        return False
    blocks = parse_namelists(deck)
    sol_blocks = list(iter_namelist_blocks(blocks, "SOLENOID"))
    if not sol_blocks:
        return False
    sol = sol_blocks[-1]  # 重复块取最后一块 (Fortran 读取语义)
    if not get_ci(sol, "LBfield", False):
        return False
    display(HTML(
        "<div style='background:#fff8e1;border:1px solid #e0c36a;"
        "padding:6px 10px;margin:4px 0'>注意: 该 deck 含螺线管场, "
        "以下统计未应用正则动量 (bz=0), 发射度/散角与 ASTRA Xemit "
        "口径可能不一致 (手册 4.13.1); 参见 stats_validation_demo "
        "或显式传入 bz_on_axis_T。</div>"))
    return True


_UNITS = {
    "total_charge_nC": "nC",
    "mean_x_mm": "mm", "mean_y_mm": "mm", "mean_z_mm": "mm",
    "sig_x_mm": "mm", "sig_y_mm": "mm", "sig_z_mm": "mm",
    "mean_pz_MeVc": "MeV/c", "sig_p_over_p_pct": "%",
    "mean_E_kin_MeV": "MeV", "sig_E_keV": "keV", "sig_E_over_E_pct": "%",
    "sig_xp_mrad": "mrad", "sig_yp_mrad": "mrad",
    "emit_x_geom_um": "um.rad", "emit_y_geom_um": "um.rad",
    "emit_x_norm_um": "um.rad", "emit_y_norm_um": "um.rad",
    "emit_x_norm_mm_mrad": "mm.mrad", "emit_y_norm_mm_mrad": "mm.mrad",
    "beta_x_m": "m", "alpha_x": "1", "gamma_t_x": "1/m",
    "beta_y_m": "m", "alpha_y": "1", "gamma_t_y": "1/m",
    "ref_momentum_MeVc": "MeV/c", "ref_kinetic_energy_MeV": "MeV",
    "gamma": "1", "beta_rel": "1",
}
