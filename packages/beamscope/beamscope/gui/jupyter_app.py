"""Interactive Jupyter notebook GUI using ipywidgets.

Usage in a Jupyter notebook:
    from beamscope.gui.jupyter_app import BeamExplorer
    explorer = BeamExplorer()
    explorer.display()
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from ..distribution import Distribution


class BeamExplorer:
    """Interactive beam data explorer for Jupyter notebooks.

    Requires ipywidgets to be installed.
    """

    def __init__(self):
        self.distribution: Optional[Distribution] = None

    def load(self, path: str | Path) -> Distribution:
        """Load a distribution file and store it."""
        from ..io import read_distribution

        self.distribution = read_distribution(Path(path))
        print(self.distribution.summary())
        return self.distribution

    def stats(self) -> None:
        """Print beam statistics."""
        if self.distribution is None:
            print("No distribution loaded. Call .load() first.")
            return
        from ..analysis.statistics import compute_statistics, print_statistics

        stats = compute_statistics(self.distribution)
        print_statistics(stats)

    def plot_phase_space(self, plane: str = "x") -> None:
        """Plot phase space."""
        if self.distribution is None:
            print("No distribution loaded. Call .load() first.")
            return
        from ..plot.phase_space import plot_phase_space
        import matplotlib.pyplot as plt

        fig = plot_phase_space(self.distribution, plane=plane)
        plt.show()

    def plot_dashboard(self) -> None:
        """Display the full dashboard."""
        if self.distribution is None:
            print("No distribution loaded. Call .load() first.")
            return
        from ..plot.dashboard import plot_dashboard
        import matplotlib.pyplot as plt

        fig = plot_dashboard({"beam": self.distribution})
        plt.show()

    def display(self) -> None:
        """Display the interactive widget interface (requires ipywidgets)."""
        try:
            import ipywidgets as widgets
            from IPython.display import display
        except ImportError:
            print("ipywidgets is required for interactive GUI. Install with: pip install ipywidgets")
            print("Use the method-based API instead: .load(), .stats(), .plot_dashboard()")
            return

        file_input = widgets.Text(
            placeholder="Path to .ini file",
            description="File:",
            layout=widgets.Layout(width="400px"),
        )
        load_btn = widgets.Button(description="Load")
        output = widgets.Output()

        plot_type = widgets.Dropdown(
            options=["dashboard", "phase_x", "phase_y", "phase_z", "projections"],
            description="Plot:",
        )
        plot_btn = widgets.Button(description="Plot")
        stats_btn = widgets.Button(description="Stats")

        def on_load(_):
            with output:
                output.clear_output()
                try:
                    self.load(file_input.value)
                    print("✓ Loaded successfully.")
                except Exception as e:
                    print(f"✗ Error: {e}")

        def on_plot(_):
            with output:
                output.clear_output()
                if self.distribution is None:
                    print("Load a file first.")
                    return
                import matplotlib.pyplot as plt
                if plot_type.value == "dashboard":
                    self.plot_dashboard()
                elif plot_type.value == "projections":
                    from ..plot.distributions import plot_distributions
                    plot_distributions(self.distribution)
                    plt.show()
                else:
                    plane = plot_type.value.split("_")[1]
                    self.plot_phase_space(plane)

        def on_stats(_):
            with output:
                output.clear_output()
                self.stats()

        load_btn.on_click(on_load)
        plot_btn.on_click(on_plot)
        stats_btn.on_click(on_stats)

        display(
            widgets.HBox([file_input, load_btn]),
            widgets.HBox([plot_type, plot_btn, stats_btn]),
            output,
        )
