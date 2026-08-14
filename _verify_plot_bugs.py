import sys
sys.path.insert(0, "/Users/yuxinwu/my_projects/astra_notebook")
import numpy as np
from pathlib import Path
from astra_tools.io import read_distribution
from astra_tools.io.astra_emit import read_emit_files

DATA = Path("/Users/yuxinwu/my_projects/astra_notebook/examples/Manual_Example")

zemit = np.loadtxt(DATA / "Example.Zemit.001")
last = zemit[-1]
print("Zemit last row:", last)
print("file col5 (eps_zn, keV.mm) =", last[5])
print("sig_E[keV]*sig_z[mm] =", (last[4]) * (last[3]))

emit = read_emit_files(str(DATA / "Example"))
print("emit.z.emit[-1] (eV.m internal) =", emit.z.emit[-1])
print("emit.z.emit[-1] * 1e-3 (as plotted) =", emit.z.emit[-1]*1e-3)
print("expected display 'keV mm' =", last[5])
print("ratio plotted/expected =", (emit.z.emit[-1]*1e-3)/last[5])

dist = read_distribution(DATA / "Example.0150.001")
m = dist.active
print("charge min/max:", dist.charge[m].min(), dist.charge[m].max())
print("all charge negative?", np.all(dist.charge[m] < 0))
print("ref_momentum_eVc =", dist.ref_momentum_eVc)
print("mean pz =", np.mean(dist.pz[m]))
