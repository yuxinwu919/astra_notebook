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
        if isinstance(v, float):
            v = "%.6g" % v
        rows.append("<tr><td>%s</td><td align='right'>%s</td><td>%s</td></tr>"
                    % (k, v, unit))
    body = "\n".join(rows)
    return HTML(
        "<table style='border-collapse:collapse;font-family:monospace'>"
        "<tr><th>quantity</th><th>value</th><th>unit</th></tr>" + body + "</table>"
    )


def distribution_summary_html(dist) -> HTML:
    """分布概要面板."""
    return HTML("<pre>" + dist.summary() + "</pre>")


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
