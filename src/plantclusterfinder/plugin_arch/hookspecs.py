"""Hook specifications for the PlantClusterFinder plugin system.

This module defines the public extension points exposed by the framework.
"""

import pluggy
import polars as pl

from plantclusterfinder.plugin_arch.context import ReturnContext, RunContext

hookspec = pluggy.HookspecMarker("plantclusterfinder")
hookimpl = pluggy.HookimplMarker("plantclusterfinder")


class PlantClusterSpecs:
    """Hook specifications for PlantClusterFinder plugins."""

    expected_outputs: None | list[str]

    @hookspec
    def required_outputs(self) -> list[str]:
        """Return a list with necessary `output_dict` keys populated."""
        return []

    @hookspec
    def process_results(self, results: pl.DataFrame, context: RunContext) -> ReturnContext:
        """Receive a Polars DataFrame from main.py and return any result."""
        raise NotImplementedError("This is an abstract class for plugins")
