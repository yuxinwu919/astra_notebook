import sys
sys.path.insert(0, "/Users/yuxinwu/my_projects/astra_notebook")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astra_tools.plot.advanced_plots import plot_envelope_with_aperture, aperture_elements
from astra_tools.distribution import Distribution

# --- Bug E: aperture offset swap
aps = [{"z1": 0.0, "z2": 1.0, "r": 0.01, "xoff": 0.005, "yoff": -0.003, "file": ""}]
# fake emit
class E:
    z = np.linspace(0, 1.5, 10)
    rms = np.full(10, 0.001)
class ES:
    x = E(); y = E()
emit = ES()

fig = plot_envelope_with_aperture(emit, aps, plane="x")
ax = fig.axes[0]
# find the aperture Rectangle patch
from matplotlib.patches import Rectangle
for p in ax.patches:
    if isinstance(p, Rectangle):
        print("plane=x aperture rectangle y-extent:", p.get_y(), "..", p.get_y()+p.get_height())
        print("  (expect centered at xoff=0.005 -> [0.005-0.01, 0.005+0.01]*1e3 =", (-0.005*1e3, 0.015*1e3), ")")
plt.close(fig)

fig = plot_envelope_with_aperture(emit, aps, plane="y")
ax = fig.axes[0]
for p in ax.patches:
    if isinstance(p, Rectangle):
        print("plane=y aperture rectangle y-extent:", p.get_y(), "..", p.get_y()+p.get_height())
        print("  (expect centered at yoff=-0.003 -> [(-0.003-0.01), (-0.003+0.01)]*1e3 =", (-0.013*1e3, 0.007*1e3), ")")
plt.close(fig)

# --- zero-momentum z phase space
from astra_tools.plot.phase_space import plot_phase_space
n = 100
d0 = Distribution.from_arrays(
    np.random.default_rng(0).normal(0, 1e-3, n),
    np.random.default_rng(1).normal(0, 1e-3, n),
    np.random.default_rng(2).normal(0, 1e-3, n),
    np.zeros(n), np.zeros(n), np.zeros(n),  # pz = 0
    np.zeros(n), np.full(n, 1e-3),
    status=np.full(n, 5), ref_momentum_eVc=0.0)
print("\n=== plot_phase_space z with pz=0 ===")
try:
    fig = plot_phase_space(d0, plane="z")
    print("NO CRASH")
    plt.close(fig)
except Exception as e:
    print("CRASH/ERR: %s: %s" % (type(e).__name__, e))
