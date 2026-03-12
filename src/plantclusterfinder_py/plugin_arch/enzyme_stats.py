"""Built-in plugin for computing enzyme statistics.

This plugin demonstrates how to consume PlantClusterFinder results
and produce a derived analysis without modifying core pipeline behavior.
"""

from pathlib import Path

import pluggy
import polars as pl

hookimpl = pluggy.HookimplMarker("plantclusterfinder")


class EnzymeStatsPlugin:
    """Enzyme Statistics: calculates enzyme stats from results."""

    @hookimpl
    def process_results(self, results: pl.DataFrame, context) -> dict[str, Path]:
        """Find the number of ___ sharing an enzyme."""
        stats = results.group_by("enzyme").count().sort("count", descending=True)
        out = context.output_dir / "enzyme_stats.csv"
        stats.write_csv(out)
        return {"enzyme_stats": out}
