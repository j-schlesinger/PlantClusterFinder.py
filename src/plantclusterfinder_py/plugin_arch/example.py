"""Built-in plugin for computing enzyme statistics.

This plugin demonstrates how to consume PlantClusterFinder results
and produce a derived analysis without modifying core pipeline behavior.
"""

import polars as pl

from plantclusterfinder_py.plugin_arch.hookspecs import hookimpl


class ExamplePlugin:
    """A hook specification namespace."""

    @hookimpl
    def process_dataframe(self, df: pl.DataFrame):
        """Example plugin logic: process dataframe to get non-null genes."""
        processed_df = df.filter(pl.col("gene").is_not_null())

        print("Processing Dataframe")

        return processed_df
