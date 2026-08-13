"""Plugin manager for PlantClusterFinder.

This module is responsible for:
- Initializing the Pluggy plugin manager
- Registering built-in and external plugins
- Executing plugin hooks in a controlled manner

This script is explicitly for testing the plugins
"""

import sys
from pathlib import Path

import polars as pl

from plantclusterfinder_py.plugin_arch.context import RunContext
from plantclusterfinder_py.plugin_arch.enzyme_stats import EnzymeStatsPlugin
from plantclusterfinder_py.plugin_arch.example import ExamplePlugin
from plantclusterfinder_py.plugin_arch.manager import PluginManager
from plantclusterfinder_py.plugin_arch.overlap_clusters import OverlapClustersPlugin
from plantclusterfinder_py.plugin_arch.visualize import VisualizationPlugin


def get_plugin_manager():
    """Get the plugin manager for PlantClusterFinder."""
    pm = PluginManager()
    # register plugins
    pm.register(ExamplePlugin())
    pm.register(EnzymeStatsPlugin())
    # pm.register(OverlapClustersPlugin())
    # pm.register(VisualizationPlugin())

    return pm


def run_plugins(df: pl.DataFrame, output_dir: Path, genome_dir: Path):
    """Run the plugins defined by the input DataFrame."""
    pm = get_plugin_manager()
    run_context = RunContext(
        run_name="TestPlugins", output_dir=output_dir, genome_file=genome_dir, config={}
    )
    return pm.run(df, run_context)


def main(argv: list[str] | None = None):
    """Main function execution: test execution of plugins."""
    if argv is None:
        argv = []
    gene_count = 4
    df = pl.DataFrame(
        {
            "True Chromosome": ["scaffold_1"] * gene_count,
            "Chromosome": ["scaffold_1"] * gene_count,
            "Gene Name": ["g1", "g2", "g3", "g4"],
            "Protein Name": ["PAL", "C4H", "4CL", "G4PTN"],
            "Start bp": [10000 * x + 50000 for x in range(gene_count)],
            "End bp": [10000 * x + 58000 for x in range(gene_count)],
            "EC-NUMBER": [["1.1.1.2"], None, ["3.1.2"], ["3.1.2", "3.1.2.2"]],
        }
    )
    if len(argv) < 2:
        print("""ERROR: wrong number of arguments;
    plugin_main usage:
        plugin_main genome_dir output_dir
        """)
        raise Exception("Cannot run plugin main due to input errors")
    output_dir = Path(argv[-1])
    genome_dir = Path(argv[-2])

    results = run_plugins(df, output_dir, genome_dir)

    print(results)


if __name__ == "__main__":
    main(sys.argv)
