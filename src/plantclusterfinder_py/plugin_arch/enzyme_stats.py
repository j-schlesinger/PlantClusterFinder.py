"""Built-in plugin for computing enzyme statistics.

This plugin demonstrates how to consume PlantClusterFinder results
and produce a derived analysis without modifying core pipeline behavior.
"""

import logging

import polars as pl

from plantclusterfinder_py.plugin_arch.context import ReturnContext, RunContext
from plantclusterfinder_py.plugin_arch.hookspecs import hookimpl

_log = logging.getLogger(__name__)


class EnzymeStatsPlugin:
    """Enzyme Statistics: calculates enzyme stats from results."""

    @hookimpl
    def process_results(self, results: pl.DataFrame, context: RunContext) -> ReturnContext:
        """Find the number of ___ sharing an enzyme."""
        stats = (
            results.explode("EC-NUMBER")
            .group_by("EC-NUMBER")
            .agg(count=pl.len())
            # .select(pl.col("EC-NUMBER"), pl.col("len").alias("count"))
            .sort("count", descending=True)
        )
        out = context.output_dir / f"enzyme_stats.{context.run_name}.csv"
        try:
            stats.write_csv(out)
        except pl.exceptions.ComputeError:
            _log.error("Cannot write stats DF to CSV: %s", stats)
        return ReturnContext(plugin_name="enzyme_stats", output_filename=out)
