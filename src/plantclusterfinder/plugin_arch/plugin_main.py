"""Plugin manager for PlantClusterFinder.

This module is responsible for:
- Initializing the Pluggy plugin manager
- Registering built-in and external plugins
- Executing plugin hooks in a controlled manner

This script is explicitly for testing the plugins
"""

import logging
from pathlib import Path

import polars as pl

from plantclusterfinder.plugin_arch.context import RunContext
from plantclusterfinder.plugin_arch.enzyme_stats import EnzymeStatsPlugin
from plantclusterfinder.plugin_arch.example import ExamplePlugin
from plantclusterfinder.plugin_arch.manager import PluginManager
from plantclusterfinder.plugin_arch.overlap_clusters import OverlapClustersPlugin
from plantclusterfinder.plugin_arch.visualize import VisualizationPlugin

_log = logging.getLogger()


def get_plugin_manager():
    """Get the plugin manager for PlantClusterFinder."""
    pm = PluginManager()
    # register plugins
    pm.register(ExamplePlugin())
    pm.register(EnzymeStatsPlugin())
    pm.register(OverlapClustersPlugin())
    pm.register(VisualizationPlugin())

    return pm


def run_plugins(df: pl.DataFrame, output_dir: Path, genome_file: Path, cluster_filename: Path):
    """Run the plugins defined by the input DataFrame."""
    pm = get_plugin_manager()
    output_dict = {"cluster_filename": cluster_filename}
    run_context = RunContext(
        run_name="TestPlugins",
        output_dir=output_dir,
        genome_file=genome_file,
        config={},
        output_dict=output_dict,
    )
    return pm.run(df, run_context)


def test_plugins(output_dir: Path, genome_file: Path, cluster_filename: Path) -> dict[str, Path]:
    """Test plugins on a sample dataset."""
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
    _log.debug("Testing plugins with the GeneOutput dataframe: %s", df)
    return run_plugins(df, output_dir, genome_file, cluster_filename)


def main(argv: list[str] | None = None):
    """Main function execution: test execution of plugins."""
    if argv is None:
        argv = []
    if len(argv) < 3:
        print("""ERROR: wrong number of arguments;
    plugin_main usage:
        plugin_main cluster_filename genome_file output_dir
        """)
        raise Exception("Cannot run plugin main due to input errors")
    output_dir = Path(argv[-1])
    genome_file = Path(argv[-2])
    cluster_filename = Path(argv[-3])

    results = ", ".join(
        [
            f"{name} at {filename}"
            for name, filename in test_plugins(output_dir, genome_file, cluster_filename).items()
        ]
    )

    print(results)


if __name__ == "__main__":
    import sys

    main(sys.argv)
