"""
Hook specifications for the PlantClusterFinder plugin system.

This module defines the public extension points exposed by the framework.
"""

import pluggy
import polars as pl
from context import RunContext

hookspec = pluggy.HookspecMarker("plantclusterfinder")


class PluginSpecs:
     """
    Hook specifications describing how plugins can interact with
    PlantClusterFinder results.
    """
     
    @hookspec
    def process_results(
        self,
        results: pl.DataFrame,
        context: RunContext,
    ):
        """
        Process final clustering results.

        Should return:
        - None (side effects only)
        - OR a named result object
        """
