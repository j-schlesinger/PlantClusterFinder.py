"""
Built-in plugin for computing enzyme statistics.

This plugin demonstrates how to consume PlantClusterFinder results
and produce a derived analysis without modifying core pipeline behavior.
"""

import polars as pl
from hookspecs import hookspec
import pluggy

hookimpl = pluggy.HookimplMarker("plantclusterfinder")


class EnzymeStatsPlugin:
    
    
    @hookimpl
    def process_results(self, results: pl.DataFrame, context):
        stats = (
            results
            .group_by("enzyme")
            .count()
            .sort("count", descending=True)
        )
        out = context.output_dir / "enzyme_stats.csv"
        stats.write_csv(out)
        return {"enzyme_stats": out}
