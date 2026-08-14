"""beamscope CLI v0.2.0 — quick plotting and analysis from the command line.

Usage:
    beamscope plot bunch.ini                  # Quick overview dashboard
    beamscope compare before.ini after.ini    # Side-by-side comparison
    beamscope stats bunch.ini                 # Print beam statistics
    beamscope slice bunch.ini --n 20          # Slice analysis
    beamscope bff bunch.ini                   # Bunch form factor
    beamscope emit Example --type dashboard   # Emit/sigma evolution plots
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="beamscope v0.2.0 — Accelerator particle distribution visualization",
        prog="beamscope",
    )
    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    # plot
    plot_parser = subparsers.add_parser("plot", help="Plot overview dashboard")
    plot_parser.add_argument("file", type=str, help="Distribution file to plot")
    plot_parser.add_argument("-o", "--output", type=str, help="Output image file")
    plot_parser.add_argument("--kind", type=str, default="overview",
                             choices=["overview", "phase_x", "phase_y", "phase_z",
                                      "detail_x", "dashboard", "distributions"],
                             help="Plot type (default: overview)")
    plot_parser.add_argument("--no-tex", action="store_true")

    # compare
    cmp_parser = subparsers.add_parser("compare", help="Compare two distributions")
    cmp_parser.add_argument("file1", type=str)
    cmp_parser.add_argument("file2", type=str)
    cmp_parser.add_argument("-o", "--output", type=str)
    cmp_parser.add_argument("--type", choices=["phase_space", "projections", "statistics"],
                            default="phase_space")

    # stats
    stats_parser = subparsers.add_parser("stats", help="Print beam statistics")
    stats_parser.add_argument("file", type=str)
    stats_parser.add_argument("--weighted", action="store_true")

    # slice (NEW)
    slice_parser = subparsers.add_parser("slice", help="Slice analysis dashboard")
    slice_parser.add_argument("file", type=str)
    slice_parser.add_argument("--n", type=int, default=20, dest="n_slices")
    slice_parser.add_argument("--binning", choices=["equi_spaced", "equi_charge"],
                              default="equi_spaced")
    slice_parser.add_argument("-o", "--output", type=str)

    # bff (NEW)
    bff_parser = subparsers.add_parser("bff", help="Bunch form factor plot")
    bff_parser.add_argument("file", type=str)
    bff_parser.add_argument("-o", "--output", type=str)
    bff_parser.add_argument("--kmin", type=float, default=1.0)
    bff_parser.add_argument("--kmax", type=float, default=1e5)

    # emit
    emit_parser = subparsers.add_parser("emit", help="Plot emittance/sigma evolution")
    emit_parser.add_argument("rootname", type=str)
    emit_parser.add_argument("--run", type=str, default="001")
    emit_parser.add_argument("-o", "--output", type=str)
    emit_parser.add_argument("--type", choices=["envelope", "emittance", "energy",
                             "eigen", "ref", "dashboard", "transverse"],
                             default="dashboard")
    emit_parser.add_argument("--no-tex", action="store_true")

    return parser.parse_args()


def _save_or_show(fig, args) -> None:
    output = getattr(args, 'output', None)
    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved: {output}")
    else:
        import matplotlib.pyplot as plt
        plt.show()


def _run_emit_command(args) -> None:
    from beamscope.io.astra_emit import read_emit_files, read_sigma_file, read_ref_file
    from beamscope.style.rcparams import set_publication_style
    set_publication_style(use_tex=not getattr(args, 'no_tex', True))
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rootname = args.rootname
    emit = read_emit_files(rootname, args.run)
    plot_type = args.type

    if plot_type == "dashboard":
        sigma = None
        try: sigma = read_sigma_file(rootname, args.run)
        except Exception: pass
        from beamscope.plot.emit_plots import plot_emit_dashboard
        fig = plot_emit_dashboard(emit, sigma)
    elif plot_type == "envelope":
        from beamscope.plot.emit_plots import plot_envelope_evolution
        fig = plot_envelope_evolution(emit)
    elif plot_type == "emittance":
        from beamscope.plot.emit_plots import plot_emittance_evolution
        fig = plot_emittance_evolution(emit)
    elif plot_type == "energy":
        from beamscope.plot.emit_plots import plot_energy_evolution
        fig = plot_energy_evolution(emit)
    elif plot_type == "eigen":
        sigma = read_sigma_file(rootname, args.run)
        from beamscope.plot.emit_plots import plot_eigen_emittances
        fig = plot_eigen_emittances(sigma)
    elif plot_type == "ref":
        ref = read_ref_file(rootname, args.run)
        from beamscope.plot.emit_plots import plot_ref_trajectory
        fig = plot_ref_trajectory(ref)
    elif plot_type == "transverse":
        from beamscope.plot.emit_plots import plot_transverse_size
        fig = plot_transverse_size(emit)
    else:
        raise ValueError(f"Unknown emit plot type: {plot_type}")
    _save_or_show(fig, args)


def main() -> None:
    args = _parse_args()
    if args.command is None:
        print("beamscope v0.2.0 — subcommands: plot, compare, stats, slice, bff, emit")
        sys.exit(0)

    try:
        from beamscope.io import read_distribution
    except ImportError as e:
        print(f"Error: Cannot import beamscope. ({e})")
        sys.exit(1)

    # stats
    if args.command == "stats":
        dist = read_distribution(Path(args.file))
        from beamscope.analysis.statistics import compute_statistics, print_statistics
        s = compute_statistics(dist, use_weights=getattr(args, 'weighted', False))
        print_statistics(s)
        return

    # emit
    if args.command == "emit":
        _run_emit_command(args)
        return

    # All below need matplotlib
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from beamscope.style.rcparams import set_publication_style
    set_publication_style(use_tex=not getattr(args, 'no_tex', True))

    # slice
    if args.command == "slice":
        dist = read_distribution(Path(args.file))
        from beamscope.analysis.slices import compute_slice_analysis
        from beamscope.plot.slice_plots import plot_slice_dashboard
        sa = compute_slice_analysis(dist, n_slices=args.n_slices, binning=args.binning)
        fig = plot_slice_dashboard(sa, title=f"Slice — {Path(args.file).stem}")
        _save_or_show(fig, args)
        return

    # bff
    if args.command == "bff":
        dist = read_distribution(Path(args.file))
        d = dist.filter_active()
        from beamscope.analysis.bff import compute_bff
        from beamscope.plot.bff_plots import plot_bff_with_amplitude
        bff_r = compute_bff(d.z, d.charge, kmin=args.kmin, kmax=args.kmax, detect_features=True)
        fig = plot_bff_with_amplitude(bff_r, title=f"BFF — {Path(args.file).stem}")
        _save_or_show(fig, args)
        return

    # plot
    if args.command == "plot":
        dist = read_distribution(Path(args.file))
        kind = getattr(args, 'kind', 'overview')
        label = Path(args.file).stem
        if kind == "overview":
            from beamscope.plot.overview import plot_overview
            fig, _ = plot_overview(dist, title=f"Overview — {label}")
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
            fig = plot_dashboard({label: dist})
        elif kind == "distributions":
            from beamscope.plot.distributions import plot_distributions
            fig = plot_distributions(dist)
        else:
            raise ValueError(f"Unknown plot kind: {kind}")
        _save_or_show(fig, args)
        return

    # compare
    if args.command == "compare":
        d1 = read_distribution(Path(args.file1))
        d2 = read_distribution(Path(args.file2))
        from beamscope.plot.comparison import plot_comparison
        fig = plot_comparison(
            {Path(args.file1).stem: d1, Path(args.file2).stem: d2},
            plot_type=args.type)
        _save_or_show(fig, args)
        return


if __name__ == "__main__":
    main()
