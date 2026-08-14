"""人性化数据导出: CSV (带单位表头注释) 与 .npz (原始数组).

目的: notebook 内置图不满意时, 用户可用导出的原始数据自行绘图。
CSV 文件头含单位与物理量注释, 开箱即懂。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def _write_csv_units(arrays: dict, path: Path, units: dict, note: str = ""):
    """CSV 导出: 每列一个物理量, 表头注释含单位与说明."""
    import pandas as pd

    df = pd.DataFrame(arrays)
    path = Path(path)
    with open(path, "w") as f:
        if note:
            for line in note.splitlines():
                f.write("# " + line + "\n")
        f.write("# columns (unit): ")
        f.write("; ".join("%s [%s]" % (k, units.get(k, "-")) for k in arrays))
        f.write("\n")
        df.to_csv(f, index=False)
    return path


def export_distribution(dist, out_dir, stem: str = "distribution"):
    """导出相空间分布: <stem>.npz (原始数组) + <stem>.csv (表格).

    CSV 列: x[m] y[m] z[m] px[eV/c] py[eV/c] pz[eV/c] t[s]
            q[nC] status index(可选)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    arrays = {
        "x": dist.x, "y": dist.y, "z": dist.z,
        "px": dist.px, "py": dist.py, "pz": dist.pz,
        "t": dist.clock, "q": dist.charge, "status": dist.status,
    }
    if dist.index is not None:
        arrays["index"] = dist.index
    units = {"x": "m", "y": "m", "z": "m", "px": "eV/c", "py": "eV/c",
             "pz": "eV/c", "t": "s", "q": "nC", "status": "1", "index": "1"}

    npz_path = out_dir / (stem + ".npz")
    np.savez(npz_path, **arrays,
             ref_time_ns=dist.ref_time_ns, ref_momentum_eVc=dist.ref_momentum_eVc,
             ref_z_m=dist.ref_z_m, total_charge_nC=dist.total_charge_nC,
             source=dist.source, format=dist.format)

    csv_path = _write_csv_units(
        arrays, out_dir / (stem + ".csv"), units,
        note=("ASTRA particle distribution exported by astra-notebook.\n"
              "reference: time=%.6g ns, p=%.6g eV/c, z=%.6g m, Q=%.6g nC"
              % (dist.ref_time_ns, dist.ref_momentum_eVc, dist.ref_z_m,
                 dist.total_charge_nC)))
    return {"npz": npz_path, "csv": csv_path}


def export_statistics(stats, out_dir, stem: str = "statistics"):
    """导出束团统计表: <stem>.csv (显示单位)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    d = stats.to_dict()
    df_rows = [{"quantity": k, "value": v, "unit": _UNIT_HINTS.get(k, "")}
               for k, v in d.items()]
    import pandas as pd
    path = out_dir / (stem + ".csv")
    pd.DataFrame(df_rows).to_csv(path, index=False)
    return path


_UNIT_HINTS = {
    "n_particle": "1", "n_active": "1", "total_charge_nC": "nC",
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


def export_emit(emit, out_dir, stem: str = "emit"):
    """导出 Xemit/Yemit/Zemit 演化数据: <stem>_{x,y,z}.csv (SI 列)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    units = {"z": "m", "t": "s", "avg": "m/eV", "rms": "m",
             "rmsprime": "rad/eV", "emit": "m.rad / eV.m", "corr": "-"}
    written = {}
    for key in ("x", "y", "z"):
        e = getattr(emit, key)
        arrays = {"z": e.z, "t": e.t, "avg": e.avg, "rms": e.rms,
                  "rmsprime": e.rmsprime, "emit": e.emit, "corr": e.corr}
        written[key] = _write_csv_units(
            arrays, out_dir / ("%s_%s.csv" % (stem, key)), units,
            note="ASTRA %semit data exported by astra-notebook (SI units)." % key.upper())
    return written
