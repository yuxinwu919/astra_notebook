import sys
sys.path.insert(0, "/Users/yuxinwu/my_projects/astra_notebook")
import numpy as np
from astra_tools.distribution import Distribution
from astra_tools.analysis.emittance import compute_emittance_ellipse_params

n = 2000
rng = np.random.default_rng(0)
x = rng.normal(3e-3, 1e-4, n)
px = rng.normal(0, 1e-3, n)
p_ref = 1e6
u = x - np.mean(x)
up = (px - np.mean(px)) / p_ref
print("sigma_u [m]:", np.std(u))
print("sigma_up [rad]:", np.std(up))
par = compute_emittance_ellipse_params(u, up)
print("a [m]:", par["a"], " b:", par["b"], " theta:", par["theta"], " eps:", par["eps"])
print("a*1e3 [mm]:", par["a"]*1e3)
print("b*1e3 [mrad]:", par["b"]*1e3)
