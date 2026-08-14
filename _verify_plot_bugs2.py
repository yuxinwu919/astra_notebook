import sys
sys.path.insert(0, "/Users/yuxinwu/my_projects/astra_notebook")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from astra_tools.io import read_distribution
from astra_tools.io.astra_emit import read_emit_files
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.plot.overview import plot_overview
from astra_tools.plot.slice_plots import plot_current_profile
from astra_tools.plot.advanced_plots import plot_core_brightness

DATA = Path("/Users/yuxinwu/my_projects/astra_notebook/examples/Manual_Example")
dist = read_distribution(DATA / "Example.0150.001")
emit = read_emit_files(str(DATA / "Example"))
sa = compute_slice_analysis(dist, n_slices=20)

# --- Bug A: overview use_weights with negative charge
print("=== overview use_weights=True (electron bunch, negative charge) ===")
try:
    fig, axes = plot_overview(dist, use_weights=True)
    print("NO CRASH")
    plt.close(fig)
except Exception as e:
    print("CRASH: %s: %s" % (type(e).__name__, e))

# --- Bug C: current profile Q sign
fig = plot_current_profile(sa)
txt = [t.get_text() for t in fig.axes[0].texts]
print("\n=== current profile text ===", txt)
print("sa.charge sum (signed, nC):", np.sum(sa.charge))
print("|sum| =", abs(np.sum(sa.charge)))
plt.close(fig)

# --- Bug D: core brightness 1e-9 factor
z = np.linspace(0, 1.5, 10)
landf = dict(landf_z=z, landf_total_charge=np.full(10, 1e-9))  # 1 nC in C
ce = dict(mean_z=z, norm_emit_x=np.full(10, 1e-6), norm_emit_y=np.full(10, 1e-6))
fig = plot_core_brightness(ce, landf)
line = fig.axes[0].lines[0]
print("\n=== core brightness ===")
print("computed b[0] =", line.get_ydata()[0])
print("expected Q/eps^2 = 1e-9 / 1e-12 =", 1e-9/1e-12)
print("ylabel:", fig.axes[0].get_ylabel())
plt.close(fig)
