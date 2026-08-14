# beamscope — Accelerator Particle Distribution Visualization

`beamscope` is a Python toolkit for loading, analyzing, and visualizing
particle beam distributions from accelerator simulation codes.

## Features

- **Multi-format I/O** — Read ASTRA binary/ASCII distributions (ELEGANT, CST, ECHO2D planned)
- **Beam statistics** — RMS sizes, geometric/normalized emittance, Twiss parameters, energy spread
- **Phase space plots** — Transverse (x-x', y-y') and longitudinal (z-δp/p) 2D histograms
- **Projection distributions** — x, y, z histograms with Gaussian fits
- **Multi-file comparison** — Side-by-side phase space, overlaid projections, parameter bar charts
- **Dashboard** — Comprehensive multi-panel beam overview
- **Slice analysis** — Longitudinal slice-by-slice emittance and current profile
- **Publication-quality styling** — SLAC-DESY beam colormap, LaTeX optional

## Installation

```bash
pip install -e packages/beamscope
```

For GUI support:
```bash
pip install -e "packages/beamscope[gui]"
```

## Quick Start

```python
from beamscope.io import read_distribution
from beamscope.analysis.statistics import compute_statistics, print_statistics
from beamscope.plot.dashboard import plot_dashboard

# Load
dist = read_distribution("bunch.ini")
print(dist.summary())

# Statistics
stats = compute_statistics(dist)
print_statistics(stats)

# Dashboard
fig = plot_dashboard({"beam": dist})
fig.savefig("dashboard.png")
```

## CLI Usage

```bash
beamscope stats bunch.ini
beamscope plot bunch.ini -o output.png
beamscope compare before.ini after.ini --type phase_space
```

## License

MIT
