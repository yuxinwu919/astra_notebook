#!/usr/bin/env python3
"""beamscope v0.2.0 — Accelerator Particle Distribution Visualization & Analysis

Double-click to launch GUI, or run with subcommands from terminal:

    python beamscope_app.py                        # Launch GUI (default)
    python beamscope_app.py stats bunch.ini        # Print beam statistics
    python beamscope_app.py plot bunch.ini          # Quick 6D overview plot
    python beamscope_app.py slice bunch.ini         # Slice analysis dashboard
    python beamscope_app.py bff bunch.ini           # Bunch form factor plot
    python beamscope_app.py compare A.ini B.ini     # Side-by-side comparison
    python beamscope_app.py emit Example --run 001  # Emittance evolution

Requires: beamscope package (pip install -e packages/beamscope)
"""

import sys
from pathlib import Path

# Ensure beamscope is importable
_project_root = Path(__file__).resolve().parent
_beamscope_path = _project_root / "packages" / "beamscope"
if str(_beamscope_path) not in sys.path:
    sys.path.insert(0, str(_beamscope_path))


# ---------------------------------------------------------------------------
# CLI mode
# ---------------------------------------------------------------------------
def _cli_stats(dist_path: str) -> None:
    """Print comprehensive beam statistics."""
    from beamscope.io import read_distribution
    from beamscope.analysis.statistics import compute_statistics, print_statistics

    dist = read_distribution(Path(dist_path))
    stats = compute_statistics(dist)
    print_statistics(stats)


def _cli_plot(dist_path: str, output: str | None = None, kind: str = "overview") -> None:
    """Generate a quick plot of a distribution."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from beamscope.io import read_distribution
    from beamscope.style.rcparams import set_publication_style

    set_publication_style(use_tex=False, dpi=150)
    dist = read_distribution(Path(dist_path))
    label = Path(dist_path).stem

    if kind == "overview":
        from beamscope.plot.overview import plot_overview
        fig, _ = plot_overview(dist, title=f"6D Overview — {label}")
    elif kind == "phase_x":
        from beamscope.plot.phase_space import plot_phase_space
        fig = plot_phase_space(dist, plane="x", show_ellipse=True)
    elif kind == "phase_y":
        from beamscope.plot.phase_space import plot_phase_space
        fig = plot_phase_space(dist, plane="y", show_ellipse=True)
    elif kind == "phase_z":
        from beamscope.plot.phase_space import plot_phase_space
        fig = plot_phase_space(dist, plane="z")
    elif kind == "detail_x":
        from beamscope.plot.detail import plot_detail
        fig = plt.figure(figsize=(8, 7))
        plot_detail(fig, dist, x_key="x", y_key="xp", show_ellipse=True)
    elif kind == "dashboard":
        from beamscope.plot.dashboard import plot_dashboard
        fig = plot_dashboard({label: dist}, title=f"Dashboard — {label}")
    elif kind == "distributions":
        from beamscope.plot.distributions import plot_distributions
        fig = plot_distributions(dist)
    else:
        raise ValueError(f"Unknown plot kind: {kind}")

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved: {output}")
    else:
        out_path = _project_root / f"beamscope_{kind}_{label}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")
    plt.close(fig)


def _cli_slice(dist_path: str, n_slices: int = 20, output: str | None = None) -> None:
    """Generate slice analysis dashboard."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from beamscope.io import read_distribution
    from beamscope.analysis.slices import compute_slice_analysis
    from beamscope.plot.slice_plots import plot_slice_dashboard
    from beamscope.style.rcparams import set_publication_style

    set_publication_style(use_tex=False, dpi=150)
    dist = read_distribution(Path(dist_path))
    label = Path(dist_path).stem

    sa = compute_slice_analysis(dist, n_slices=n_slices)
    fig = plot_slice_dashboard(sa, title=f"Slice Analysis — {label} ({n_slices} slices)")

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved: {output}")
    else:
        out_path = _project_root / f"beamscope_slice_{label}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")
    plt.close(fig)


def _cli_bff(dist_path: str, output: str | None = None) -> None:
    """Generate bunch form factor plot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from beamscope.io import read_distribution
    from beamscope.analysis.bff import compute_bff
    from beamscope.plot.bff_plots import plot_bff_with_amplitude
    from beamscope.style.rcparams import set_publication_style

    set_publication_style(use_tex=False, dpi=150)
    dist = read_distribution(Path(dist_path))
    label = Path(dist_path).stem

    d = dist.filter_active()
    bff_result = compute_bff(d.z, d.charge, detect_features=True)
    fig = plot_bff_with_amplitude(bff_result, title=f"BFF — {label}")

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved: {output}")
    else:
        out_path = _project_root / f"beamscope_bff_{label}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")
    plt.close(fig)


def _cli_compare(file1: str, file2: str, output: str | None = None) -> None:
    """Generate side-by-side comparison of two distributions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from beamscope.io import read_distribution
    from beamscope.plot.comparison import plot_comparison
    from beamscope.style.rcparams import set_publication_style

    set_publication_style(use_tex=False, dpi=150)

    dist1 = read_distribution(Path(file1))
    dist2 = read_distribution(Path(file2))
    label1, label2 = Path(file1).stem, Path(file2).stem

    fig = plot_comparison(
        {label1: dist1, label2: dist2},
        plot_type="phase_space",
    )

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved: {output}")
    else:
        out_path = _project_root / f"beamscope_compare_{label1}_vs_{label2}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Usage helper
# ---------------------------------------------------------------------------
_USAGE = """
beamscope v0.2.0 — Particle Beam Visualization & Analysis
==========================================================

Usage:
    python beamscope_app.py                              Launch GUI
    python beamscope_app.py stats <file>                 Beam statistics
    python beamscope_app.py plot <file> [--kind ...]     Quick plot
    python beamscope_app.py slice <file> [--n 20]        Slice analysis
    python beamscope_app.py bff <file>                   Bunch form factor
    python beamscope_app.py compare <f1> <f2>            Side-by-side compare

Plot kinds (for 'plot' subcommand):
    overview (default), phase_x, phase_y, phase_z, detail_x, dashboard, distributions

Examples:
    python beamscope_app.py stats simulation_files/Example.ini
    python beamscope_app.py plot simulation_files/Example.ini --kind dashboard
    python beamscope_app.py slice simulation_files/Example.ini --n 30
    python beamscope_app.py bff simulation_files/Example.ini
    python beamscope_app.py compare before.ini after.ini
"""


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 1:
        # No arguments → launch GUI
        from beamscope.gui.app import main
        main()
    elif sys.argv[1] in ("-h", "--help", "help"):
        print(_USAGE)
    elif sys.argv[1] == "stats" and len(sys.argv) >= 3:
        _cli_stats(sys.argv[2])
    elif sys.argv[1] == "plot" and len(sys.argv) >= 3:
        kind = "overview"
        output = None
        remaining = sys.argv[2:]
        i = 0
        while i < len(remaining):
            if remaining[i] == "--kind" and i + 1 < len(remaining):
                kind = remaining[i + 1]
                i += 2
            elif remaining[i] == "-o" and i + 1 < len(remaining):
                output = remaining[i + 1]
                i += 2
            else:
                filepath = remaining[i]
                i += 1
        _cli_plot(filepath, output=output, kind=kind)
    elif sys.argv[1] == "slice" and len(sys.argv) >= 3:
        n_slices = 20
        output = None
        remaining = sys.argv[2:]
        i = 0
        while i < len(remaining):
            if remaining[i] == "--n" and i + 1 < len(remaining):
                n_slices = int(remaining[i + 1])
                i += 2
            elif remaining[i] == "-o" and i + 1 < len(remaining):
                output = remaining[i + 1]
                i += 2
            else:
                filepath = remaining[i]
                i += 1
        _cli_slice(filepath, n_slices=n_slices, output=output)
    elif sys.argv[1] == "bff" and len(sys.argv) >= 3:
        output = None
        remaining = sys.argv[2:]
        i = 0
        while i < len(remaining):
            if remaining[i] == "-o" and i + 1 < len(remaining):
                output = remaining[i + 1]
                i += 2
            else:
                filepath = remaining[i]
                i += 1
        _cli_bff(filepath, output=output)
    elif sys.argv[1] == "compare" and len(sys.argv) >= 4:
        output = None
        remaining = sys.argv[2:]
        i = 0
        files = []
        while i < len(remaining):
            if remaining[i] == "-o" and i + 1 < len(remaining):
                output = remaining[i + 1]
                i += 2
            else:
                files.append(remaining[i])
                i += 1
        if len(files) >= 2:
            _cli_compare(files[0], files[1], output=output)
        else:
            print("Need two files to compare.")
    else:
        print(f"Unknown command: {' '.join(sys.argv[1:])}")
        print(_USAGE)
        sys.exit(1)
