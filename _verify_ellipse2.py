import sys
sys.path.insert(0, "/Users/yuxinwu/my_projects/astra_notebook")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astra_tools.distribution import Distribution
from astra_tools.plot.phase_space import plot_phase_space
from astra_tools.analysis.emittance import compute_emittance_ellipse_params

n = 2000
rng = np.random.default_rng(0)
x = rng.normal(3e-3, 1e-4, n)
px = rng.normal(0, 1e-3, n)
p_ref = 1e6
d = Distribution.from_arrays(
    x, rng.normal(0, 1e-4, n), rng.normal(0, 1e-4, n),
    px, rng.normal(0, 1e-3, n), np.full(n, 1e6),
    np.zeros(n), np.full(n, 1e-3),
    status=np.full(n, 5), ref_momentum_eVc=1e6)

u = x - np.mean(x)
up = (px - np.mean(px))/p_ref
par = compute_emittance_ellipse_params(u, up)
print("direct: a*1e3 =", par["a"]*1e3, " b*1e3 =", par["b"]*1e3)

fig = plot_phase_space(d, plane="x", show_ellipse=True)
ax = fig.axes[0]
print("n lines:", len(ax.lines))
for ln in ax.lines:
    xs = np.asarray(ln.get_xdata())
    print("  label=%r  xmin=%g xmax=%g  ymin=%g ymax=%g" % (ln.get_label(), xs.min(), xs.max(), np.asarray(ln.get_ydata()).min(), np.asarray(ln.get_ydata()).max()))
plt.close(fig)
