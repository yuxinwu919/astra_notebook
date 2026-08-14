#!/usr/bin/env python3
"""beamscope quickstart example.

Demonstrates the core workflow:
  1. Load a distribution file
  2. Compute and print statistics
  3. Generate plots (phase space, projections, dashboard)

Usage:
    python quickstart.py [path_to_ini_file]
"""

import sys
from pathlib import Path

# If beamscope is not installed, add the parent directory to the path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend for script runs

from beamscope.io import read_distribution
from beamscope.analysis.statistics import compute_statistics, print_statistics
from beamscope.plot.phase_space import plot_transverse_phase_space, plot_phase_space
from beamscope.plot.distributions import plot_distributions
from beamscope.plot.dashboard import plot_dashboard
from beamscope.style.rcparams import set_publication_style
from beamscope.style.colormaps import BEAM_COLORMAP


def main():
    # Use command-line argument or default path
    if len(sys.argv) > 1:
        filepath = Path(sys.argv[1])
    else:
        # Try common locations
        candidates = [
            Path("../../examples/Manual_Example/Example.ini"),
            Path("../../simulation_files/bunch.ini"),
        ]
        filepath = None
        for c in candidates:
            if c.exists():
                filepath = c.resolve()
                break
        if filepath is None:
            print("No distribution file found. Usage: python quickstart.py <path_to_ini>")
            print("Tried:", [str(c) for c in candidates])
            return

    # ── Setup publication style ──
    set_publication_style(use_tex=False)

    # ── Load distribution ──
    print(f"Loading: {filepath}")
    dist = read_distribution(filepath)
    print(dist.summary())

    # ── Statistics ──
    stats = compute_statistics(dist)
    print_statistics(stats)

    # ── Phase space plots ──
    print("Generating phase space plots...")
    fig = plot_transverse_phase_space(dist, cmap=BEAM_COLORMAP.name or "viridis")
    outdir = Path("output")
    outdir.mkdir(exist_ok=True)
    fig.savefig(outdir / "phase_space_transverse.png", dpi=150)
    print(f"  → {outdir / 'phase_space_transverse.png'}")

    fig = plot_phase_space(dist, plane="z", cmap=BEAM_COLORMAP.name or "viridis")
    fig.savefig(outdir / "phase_space_longitudinal.png", dpi=150)
    print(f"  → {outdir / 'phase_space_longitudinal.png'}")

    # ── Projection distributions ──
    fig = plot_distributions(dist)
    fig.savefig(outdir / "projections.png", dpi=150)
    print(f"  → {outdir / 'projections.png'}")

    # ── Dashboard ──
    print("Generating dashboard...")
    fig = plot_dashboard({"beam": dist})
    fig.savefig(outdir / "dashboard.png", dpi=150)
    print(f"  → {outdir / 'dashboard.png'}")

    print("\nDone! All plots saved to output/ directory.")


if __name__ == "__main__":
    main()
