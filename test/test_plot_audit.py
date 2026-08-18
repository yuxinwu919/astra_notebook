"""第 4 层扩展: 全部绘图函数逐一渲染审计 (无异常/标签/有限性).

渲染全部公开绘图函数并检查:
  * 不抛异常
  * 每个含曲线的轴都有 x/y 标签 (twinx 共享 x 轴除外)
  * 所有线/集合/轴范围数据有限
  * 字体: 任何 findfont / 缺字形警告一律视为错误 (字体回退链必须
    覆盖全部使用到的字形; 历史事故: 缺失字体家族产生上万条警告
    把 notebook 拖到几分钟不结束)
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 字体熔断: 本模块渲染全部绘图函数, findfont 或缺字形警告一旦出现
# 立即失败 (而不是让 Jupyter 前端去消化几千条警告消息)。
warnings.filterwarnings(
    "error", message=r".*(findfont|missing from font).*",
    category=UserWarning)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.plot.style import set_style
set_style()

from astra_tools.io import read_distribution, cp_index_colors
from astra_tools.io.astra_emit import (read_emit_files, read_ref_file,
                                        read_sigma_file, parse_output_file)
from astra_tools.io.field_map import (read_cavity_field,
                                      read_solenoid_field, TEField)
from astra_tools.io.astra_misc import read_pscan, read_scan
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.analysis.bff import compute_bff
from astra_tools.namelist.parse import parse_namelists

DATA = PROJECT_ROOT / "examples/Manual_Example"
CAV = PROJECT_ROOT / "examples/Cavity_Example"
DIPOLE = PROJECT_ROOT / "examples/90deg_bend_Example/3D_Dipole"


def _is_twinx(ax):
    """轴与同图其他轴共享 x (twinx) -> 不要求其 xlabel。"""
    try:
        return len(ax.get_shared_x_axes().get_siblings(ax)) > 1
    except Exception:
        return False


@pytest.fixture(scope="module")
def F():
    d = read_distribution(DATA / "Example.0150.001")
    emit = read_emit_files(str(DATA / "Example"))
    ref = read_ref_file(str(DATA / "Example"))
    sigma = read_sigma_file(str(DATA / "Example"))
    sa = compute_slice_analysis(d, n_slices=20)
    bff = compute_bff(d.filter_active().z, d.filter_active().charge,
                      kmin=10, kmax=1e5, nk=150, detect_features=True)
    cav = read_cavity_field(DATA / "3_cell_L-Band.dat")
    sol = read_solenoid_field(DATA / "Solenoid.dat").scaled(0.35)
    _tez = np.linspace(0, 0.3, 201)
    te = TEField(z=_tez, bz0=np.exp(-((_tez - 0.15) / 0.05) ** 2))
    pscan = read_pscan(CAV / "golden/astra.PScan.001")
    scan = read_scan(DATA / "Example.Scan.001")
    z = np.linspace(0, 1.5, 30)
    landf = dict(landf_z=z, landf_n_particles=np.full(30, 500.0),
                 landf_total_charge=np.full(30, 1e-9),
                 landf_n_lost=np.zeros(30),
                 landf_energy_deposited=np.linspace(0, 1e-4, 30),
                 landf_energy_exchange=np.linspace(0, 2e-5, 30))
    track = dict(seq=np.repeat([1, 2, 3], 20), status=np.ones(60, int),
                 z=np.tile(np.linspace(0, 1.5, 20), 3),
                 x=np.tile(np.linspace(0, 1e-3, 20), 3),
                 y=np.tile(np.linspace(0, -0.5e-3, 20), 3),
                 Ez=np.linspace(0, 1e6, 60), Er=np.linspace(1e5, 0, 60))
    cathode = dict(t=np.linspace(0, 1e-9, 20),
                   E_acc=np.linspace(5e6, 3e6, 20),
                   E_spch=np.linspace(-2e6, 0, 20),
                   q=np.linspace(0, 1e-3, 20))
    x2 = dict(z=np.linspace(0, 1.5, 20), K2z=np.zeros(20), K3z=np.zeros(20),
              eps_red_z=np.linspace(0.9e-6, 1e-6, 20),
              K2E=np.zeros(20), K3E=np.zeros(20),
              eps_red_zE=np.linspace(1.2e-6, 1.1e-6, 20))
    tr = dict(z=np.linspace(0, 1.5, 20), t=np.linspace(0, 5e-9, 20),
              eps_tr_x=np.linspace(1e-6, 1.1e-6, 20),
              eps_tr_y=np.linspace(1e-6, 0.9e-6, 20),
              eps_tr_z=np.linspace(1e-6, 1.2e-6, 20))
    lm = dict(z=np.linspace(0, 1.5, 20), avr=np.linspace(0, 0.5, 20),
              rms=np.linspace(0.1, 0.3, 20))
    tc = dict(z=np.linspace(0, 1.5, 20),
              scaling=np.random.default_rng(0).uniform(0.5, 1.5, (20, 5)),
              counter=np.arange(20))
    err = dict(run=np.arange(30), z=np.full(30, 1.5),
               FOM=np.random.default_rng(1).normal(1e-6, 1e-8, (30, 10)))
    cz = np.linspace(0, 1.5, 10)
    crows = np.column_stack([cz] + [
        np.full(10, 1.0), np.full(10, 0.95), np.full(10, 0.90),
        np.full(10, 0.85), np.full(10, 1.0), np.full(10, 0.95),
        np.full(10, 0.90), np.full(10, 0.85), np.full(10, 5.0),
        np.full(10, 4.8), np.full(10, 4.5), np.full(10, 4.2)])
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".Cemit.001", delete=False) as tf:
        np.savetxt(tf.name, crows)
        cpath = tf.name
    ce = parse_output_file(cpath)
    ap = parse_namelists(PROJECT_ROOT / "examples/Aperture/astra.in")["APERTURE"]
    with tempfile.NamedTemporaryFile(suffix=".quad.dat", delete=False) as tf:
        _qz = np.linspace(0, 1.0, 50)
        np.savetxt(tf.name, np.column_stack([_qz, 5 * np.ones(50),
                                             -5 * np.ones(50)]))
        qpath = tf.name
    return dict(dist=d, emit=emit, ref=ref, sigma=sigma, sa=sa, bff=bff,
                cav=cav, sol=sol, pscan=pscan, scan=scan, landf=landf,
                track=track, cathode=cathode, x2=x2, tr=tr, lm=lm, tc=tc,
                err=err, ce=ce, ap=ap, te=te, qpath=qpath)


from astra_tools.plot import (phase_space as _ps, overview as _ov,
    distributions as _ds, emit_plots as _ep, slice_plots as _sp,
    bff_plots as _bp, field_plots as _fp, advanced_plots as _ap,
    arbitrary_phase_space as _arb)


def _overlay_case(d):
    om = _arb.OverlayManager()
    om.add(d, "x", "xp")
    om.add(d, "y", "yp")
    return om.plot()


def _cases(F):
    d, emit, ref, sigma = F["dist"], F["emit"], F["ref"], F["sigma"]
    return [
        ("phase_space_x", lambda: _ps.plot_phase_space(d, plane="x")),
        ("phase_space_x_norm", lambda: _ps.plot_phase_space(d, plane="x", normalize=True)),
        ("phase_space_x_weighted", lambda: _ps.plot_phase_space(d, plane="x", use_weights=True)),
        ("phase_space_y", lambda: _ps.plot_phase_space(d, plane="y")),
        ("phase_space_z", lambda: _ps.plot_phase_space(d, plane="z")),
        ("transverse", lambda: _ps.plot_transverse_phase_space(d)),
        ("overview", lambda: _ov.plot_overview(d)),
        ("transverse_profile", lambda: _ov.plot_transverse_profile(d)),
        ("distributions", lambda: _ds.plot_distributions(d)),
        ("energy_dist", lambda: _ds.plot_energy_distribution(d)),
        ("envelope", lambda: _ep.plot_envelope_evolution(emit)),
        ("envelope_t", lambda: _ep.plot_envelope_evolution(emit, x_axis="t")),
        ("divergence", lambda: _ep.plot_divergence_evolution(emit)),
        ("emittance", lambda: _ep.plot_emittance_evolution(emit)),
        ("energy", lambda: _ep.plot_energy_evolution(emit)),
        ("bunch_length", lambda: _ep.plot_bunch_length_evolution(emit)),
        ("energy_spread", lambda: _ep.plot_energy_spread_evolution(emit)),
        ("ref_traj", lambda: _ep.plot_ref_trajectory(ref)),
        ("velocity", lambda: _ep.plot_velocity_evolution(ref)),
        ("step_size", lambda: _ep.plot_step_size_evolution(ref)),
        ("eigen", lambda: _ep.plot_eigen_emittances(sigma)),
        ("emit_dashboard", lambda: _ep.plot_emit_dashboard(emit, sigma)),
        ("lineplot_overview", lambda: _ep.plot_lineplot_overview(emit)),
        ("current", lambda: _sp.plot_current_profile(F["sa"])),
        ("slice_emit", lambda: _sp.plot_slice_emittance(F["sa"])),
        ("slice_sizes", lambda: _sp.plot_slice_sizes(F["sa"])),
        ("chirp", lambda: _sp.plot_energy_chirp(F["sa"])),
        ("slice_dashboard", lambda: _sp.plot_slice_dashboard(F["sa"])),
        ("bff", lambda: _bp.plot_bff(F["bff"])),
        ("bff_amp", lambda: _bp.plot_bff_with_amplitude(F["bff"])),
        ("cavity", lambda: _fp.plot_cavity_field(F["cav"], omega=2 * np.pi * 1.3e9)),
        ("solenoid", lambda: _fp.plot_solenoid_field(F["sol"])),
        ("te_field", lambda: _fp.plot_te_field(F["te"], omega=2 * np.pi * 1.3e9)),
        ("r3rd_tm", lambda: _fp.plot_field_expansion_radius(F["cav"], omega=2 * np.pi * 1.3e9)),
        ("r3rd_te", lambda: _fp.plot_field_expansion_radius(F["te"], omega=2 * np.pi * 1.3e9)),
        ("r3rd_solenoid", lambda: _fp.plot_field_expansion_radius(F["sol"])),
        ("solenoid_components", lambda: _fp.plot_solenoid_components(F["sol"])),
        ("quadrupole", lambda: _ap.plot_quadrupole_field(F["qpath"])),
        ("laser_envelope", lambda: _ap.plot_laser_envelope(str(CAV / "3D_test.ex"))),
        ("plasma_z", lambda: _ap.plot_plasma_fields(str(PROJECT_ROOT / "examples/Plasma_Example_1/PLASMA_flattop.txt"), peak_density_cm3=1e17, vs="z")),
        ("plasma_zeta", lambda: _ap.plot_plasma_fields(str(PROJECT_ROOT / "examples/Plasma_Example_1/PLASMA_flattop.txt"), peak_density_cm3=1e17, vs="zeta")),
        ("cathode_rings", lambda: _ap.plot_curved_cathode_contour(str(PROJECT_ROOT / "examples/Curved_Cathode_Example/Contour.dat"), show_rings=True)),
        ("losses", lambda: _ap.plot_losses(F["landf"])),
        ("beam_loading", lambda: _ap.plot_beam_loading(F["landf"])),
        ("beta_alpha", lambda: _ap.plot_beta_alpha(emit, ref=ref)),
        ("phase_advance", lambda: _ap.plot_phase_advance(emit, ref=ref)),
        ("coherence", lambda: _ap.plot_coherence_length(emit)),
        ("phase_scan", lambda: _ap.plot_phase_scan(F["pscan"])),
        ("pscan_dedz", lambda: _ap.plot_pscan_dedz(F["pscan"])),
        ("pscan_comp", lambda: _ap.plot_pscan_compression(F["pscan"])),
        ("scan_fom", lambda: _ap.plot_scan_fom(F["scan"], i=0)),
        ("corr_energy_spread", lambda: _ep.plot_correlated_energy_spread(emit)),
        ("ref_momentum", lambda: _ep.plot_ref_momentum(ref)),
        ("pscan_comp_time", lambda: _ap.plot_pscan_compression_time(F["pscan"])),
        ("scan_position", lambda: _ap.plot_scan_position(F["scan"])),
        ("tcheck_counter", lambda: _ap.plot_tcheck_counter(F["tc"])),
        ("core_emit_z", lambda: _ap.plot_core_emittance(F["ce"], plane="z")),
        ("cr_emit", lambda: _ap.plot_cr_emit(dict(z=np.linspace(0, 1.5, 20),
            eps_x=np.full(20, 1e-6), eps_y=np.full(20, 1e-6),
            q_rest=np.linspace(1, 0.8, 20), q_cross=np.linspace(0, 0.2, 20),
            x_rms=np.linspace(1e-3, 2e-3, 20),
            y_rms=np.linspace(1e-3, 1.5e-3, 20)))),
        ("error_hist", lambda: _ap.plot_error_hist(F["err"], i=0)),
        ("reduced", lambda: _ap.plot_reduced_emittance(F["x2"], F["x2"])),
        ("emit_diff", lambda: _ap.plot_emittance_difference(F["emit"], F["x2"])),
        ("corr_contrib", lambda: _ap.plot_correlated_emittance_contributions(F["x2"])),
        ("red_long", lambda: _ap.plot_reduced_longitudinal_emittance(d)),
        ("trace", lambda: _ap.plot_trace_emittance(F["tr"])),
        ("core_emit", lambda: _ap.plot_core_emittance(F["ce"])),
        ("core_brightness", lambda: _ap.plot_core_brightness(F["ce"], F["landf"])),
        ("larmor", lambda: _ap.plot_larmor(F["lm"])),
        ("tcheck", lambda: _ap.plot_tcheck_scaling(F["tc"])),
        ("z_plot", lambda: _ap.plot_z_plot(d)),
        ("probe_traj", lambda: _ap.plot_probe_trajectories(F["track"])),
        ("probe_traj_cyl", lambda: _ap.plot_probe_trajectories(F["track"], mode="cylindrical")),
        ("sc_fields", lambda: _ap.plot_space_charge_fields(F["track"])),
        ("cathode", lambda: _ap.plot_cathode_emission(F["cathode"])),
        ("slice_mismatch", lambda: _ap.plot_slice_mismatch(d, n_slices=10)),
        ("3d_map_slices", lambda: _fp.plot_3d_field_map(str(CAV / "3D_test.ex"), view="slices", component="z", n_slices=2)),
        ("3d_vector_slices", lambda: _fp.plot_3d_field_map(str(DIPOLE), view="slices", n_slices=2)),
        ("3d_quiver", lambda: _fp.plot_3d_field_quiver(str(DIPOLE), plane="xy", index=22)),
        ("3d_quiver_mag", lambda: _fp.plot_3d_field_quiver(str(DIPOLE), plane="xz", position=0.0, color_by="magnitude")),
        ("3d_contour", lambda: _fp.plot_3d_field_map(str(DIPOLE), view="slices", kind="contour", n_slices=2)),
        ("3d_stack3d", lambda: _fp.plot_3d_field_map(str(DIPOLE), view="stack3d", n_slices=3)),
        ("3d_scalar_slices", lambda: _fp.plot_3d_field_map(str(DIPOLE), view="slices", component="y", n_slices=2)),
        ("field_profile", lambda: _ap.plot_field_profile(str(DATA / "3_cell_L-Band.dat"), label="Ez", unit="MV/m")),
        ("curved_cathode", lambda: _ap.plot_curved_cathode_contour(str(PROJECT_ROOT / "examples/Curved_Cathode_Example/Contour.dat"))),
        ("laser_on_axis", lambda: _ap.plot_laser_on_axis(str(CAV / "3D_test.ex"), unit="V/m")),
        ("plasma_profile", lambda: _ap.plot_plasma_profile(str(PROJECT_ROOT / "examples/Plasma_Example_1/PLASMA_flattop.txt"), peak_density_cm3=1e17)),
        ("envelope_aperture", lambda: _ap.plot_envelope_with_aperture(emit, _ap.aperture_elements(F["ap"]))),
        ("core_fraction", lambda: _ap.plot_central_charge_fraction_curves(d)),
        ("slice_ellipses_3d", lambda: _ap.plot_slice_ellipses_3d(d, n_slices=6)),
        ("slice_ellipses_2d", lambda: _ap.plot_slice_ellipses_2d(d, n_slices=6)),
        ("slice_ellipses_2d_corr", lambda: _ap.plot_slice_ellipses_2d(d, n_slices=6, subtract_corr=True)),
        ("phase_space_cpcolors", lambda: _ps.plot_phase_space(d, plane="x", colors=cp_index_colors(d.index, {1: (1, 0, 0), 2: (0, 1, 0)}))),
        # 2026-08 全覆盖新增
        ("phase_space_t", lambda: _ps.plot_phase_space(d, plane="t")),
        ("phase_space_status", lambda: _ps.plot_phase_space(d, plane="x", color_by_status=True)),
        ("overview_time", lambda: _ov.plot_overview(d, time=True)),
        ("core_emit_curve", lambda: _ap.plot_core_emittance_curve(d)),
        ("arbitrary", lambda: _arb.plot_arbitrary(d, "x", "xp")),
        ("arbitrary_corr", lambda: _arb.plot_arbitrary(d, "x", "xp", subtract_corr=True)),
        ("arbitrary_proj", lambda: _arb.plot_arbitrary(d, "x", "xp", add_proj=True)),
        ("arbitrary_status", lambda: _arb.plot_arbitrary(d, "x", "xp", color_by_status=True)),
        ("slice_sizes_div", lambda: _sp.plot_slice_sizes(F["sa"], divergences=True)),
        ("slice_ellipses_yyp", lambda: _ap.plot_slice_ellipses_3d(d, plane="yyp")),
        ("slice_ellipses_corr", lambda: _ap.plot_slice_ellipses_3d(d, subtract_corr=True)),
        ("overlay", lambda: _overlay_case(d)),
    ]



# 72 个绘图用例逐一参数化: 任一失败可精确定位到具体图 (批 1c)。
AUDIT_CASE_IDS = ["phase_space_x","phase_space_x_norm","phase_space_x_weighted","phase_space_y","phase_space_z","transverse","overview","transverse_profile","distributions","energy_dist","envelope","envelope_t","divergence","emittance","energy","bunch_length","energy_spread","ref_traj","velocity","step_size","eigen","emit_dashboard","lineplot_overview","current","slice_emit","slice_sizes","chirp","slice_dashboard","bff","bff_amp","cavity","solenoid","te_field","r3rd_tm","r3rd_te","r3rd_solenoid","solenoid_components","quadrupole","laser_envelope","plasma_z","plasma_zeta","cathode_rings","losses","beam_loading","beta_alpha","phase_advance","coherence","phase_scan","pscan_dedz","pscan_comp","scan_fom","corr_energy_spread","ref_momentum","pscan_comp_time","scan_position","tcheck_counter","core_emit_z","cr_emit","error_hist","reduced","emit_diff","corr_contrib","red_long","trace","core_emit","core_brightness","larmor","tcheck","z_plot","probe_traj","probe_traj_cyl","sc_fields","cathode","slice_mismatch","3d_map_slices","3d_vector_slices","3d_quiver","3d_quiver_mag","3d_contour","3d_stack3d","3d_scalar_slices","field_profile","curved_cathode","laser_on_axis","plasma_profile","envelope_aperture","core_fraction","slice_ellipses_3d","slice_ellipses_2d","slice_ellipses_2d_corr","phase_space_cpcolors","phase_space_t","phase_space_status","overview_time","core_emit_curve","arbitrary","arbitrary_corr","arbitrary_proj","arbitrary_status","slice_sizes_div","slice_ellipses_yyp","slice_ellipses_corr","overlay"]


@pytest.mark.parametrize("case_name", AUDIT_CASE_IDS)
def test_all_plots_render_cleanly(F, case_name):
    """每个绘图函数: 不抛异常、标签齐全、数据有限。"""
    fn = dict(_cases(F))[case_name]
    fig = fn()
    try:
        if not isinstance(fig, plt.Figure):
            fig = fig[0]
        for ax in fig.axes:
            if ax.get_label() == "<colorbar>":
                continue
            has_lines = len(ax.lines) > 0 or len(ax.collections) > 0
            if not has_lines:
                continue
            if not _is_twinx(ax):
                assert ax.get_xlabel().strip(), case_name + ": missing xlabel"
            assert ax.get_ylabel().strip(), case_name + ": missing ylabel"
            for ln in ax.lines:
                assert np.all(np.isfinite(ln.get_xdata())), case_name + ": non-finite x"
                assert np.all(np.isfinite(ln.get_ydata())), case_name + ": non-finite y"
            for coll in ax.collections:
                arr = coll.get_array()
                if arr is not None and len(arr):
                    assert np.all(np.isfinite(arr)), case_name + ": non-finite collection"
            assert np.isfinite(ax.get_xlim()).all() and np.isfinite(ax.get_ylim()).all(), case_name + ": bad limits"
    finally:
        plt.close("all")


def test_audit_case_ids_complete(F):
    """ID 清单与 _cases 同步防漂移。"""
    assert [n for n, _ in _cases(F)] == AUDIT_CASE_IDS
