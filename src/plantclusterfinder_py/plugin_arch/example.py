"""Built-in plugin for computing enzyme statistics.

This plugin demonstrates how to consume PlantClusterFinder results
and produce a derived analysis without modifying core pipeline behavior.
"""

import logging

import polars as pl

from plantclusterfinder_py.plugin_arch.context import ReturnContext, RunContext
from plantclusterfinder_py.plugin_arch.hookspecs import hookimpl

_log = logging.getLogger(__name__)


class ExamplePlugin:
    """A hook specification namespace."""

    @hookimpl
    def process_results(self, results: pl.DataFrame, context: RunContext) -> ReturnContext:
        """Example plugin logic: process DataFrame to get non-null genes."""
        processed_df = results.filter(pl.col("Gene Name").is_not_null())
        output_filename = context.output_dir / f"example_results.{context.run_name}.csv"
        _log.info("Writing output to %s", output_filename)
        processed_df.write_csv(output_filename)
        return ReturnContext(plugin_name="example_plugin", output_filename=output_filename)
