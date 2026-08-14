import sys
sys.path.insert(0, "/Users/yuxinwu/my_projects/astra_notebook")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from astra_tools.io import read_distribution
from astra_tools.plot.phase_space import plot_phase_space

DATA = Path("/Users/yuxinwu/my_projects/astra_notebook/examples/Manual_Example")
dist = read_distribution(DATA / "Example.0150.001")
m = dist.active
print("mean x [mm]:", np.mean(dist.x[m])*1e3)
print("mean x' [mrad]:", np.mean(dist.px[m])/dist.ref_momentum_eVc*1e3)

# Build a clearly off-center beam to expose ellipse misplacement
from astra_tools.distribution import Distribution
n = 2000
rng = np.random.default_rng(0)
x = rng.normal(3e-3, 1e-4, n)   # centroid at 3 mm
px = rng.normal(0, 1e-3, n)
d = Distribution.from_arrays(
    x, rng.normal(0, 1e-4, n), rng.normal(0, 1e-4, n),
    px, rng.normal(0, 1e-3, n), np.full(n, 1e6),
    np.zeros(n), np.full(n, 1e-3),
    status=np.full(n, 5), ref_momentum_eVc=1e6)

fig = plot_phase_space(d, plane="x", show_ellipse=True)
ax = fig.axes[0]
# find the ellipse line (the one with label '1-RMS ellipse')
for ln in ax.lines:
    if ln.get_label() == "1-RMS ellipse":
        xs = ln.get_xdata()
        print("ellipse x-centroid [mm]:", (xs.min()+xs.max())/2)
        print("ellipse x-range:", xs.min(), xs.max())
        print("beam density centroid should be ~3 mm; ellipse drawn at 0 -> misplaced")
plt.close(fig)
