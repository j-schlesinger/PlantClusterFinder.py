"""Plugin manager for PlantClusterFinder.

This module is responsible for:
- Initializing the Pluggy plugin manager
- Registering built-in and external plugins
- Executing plugin hooks in a controlled manner
"""

from pathlib import Path

import polars as pl

from plantclusterfinder_py.plugin_arch.context import RunContext
from plantclusterfinder_py.plugin_arch.enzyme_stats import EnzymeStatsPlugin
from plantclusterfinder_py.plugin_arch.example import ExamplePlugin
from plantclusterfinder_py.plugin_arch.manager import PluginManager


def get_plugin_manager():
    """Get the plugin manager for PlantClusterFinder."""
    pm = PluginManager()
    # register plugins
    pm.register(ExamplePlugin())
    pm.register(EnzymeStatsPlugin())

    return pm


def run_plugins(df: pl.DataFrame):
    """Run the plugins defined by the input DataFrame."""
    pm = get_plugin_manager()

    return pm.run(df, RunContext(run_name="TestPlugins", output_dir=Path(), config={}))


def main():
    """Main function execution: test execution of plugins."""
    df = pl.DataFrame(
        {
            "Gene Name": ["g1", "g2", "g3", "g4"],
            "Protein Name": ["PAL", "C4H", "4CL", "G4PTN"],
            "EC-NUMBER": [["1.1.1.2"], None, ["3.1.2"], ["3.1.2", "3.1.2.2"]],
        }
    )

    results = run_plugins(df)

    print(results)


if __name__ == "__main__":
    main()
