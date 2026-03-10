"""
Hook specifications for the PlantClusterFinder plugin system.

This module defines the public extension points exposed by the framework.
"""


import pluggy
import polars as pl


hookspec = pluggy.HookspecMarker("plantclusterfinder")
hookimpl = pluggy.HookimplMarker("plantclusterfinder")


class PLantCLusterSpecs:
    """
    Hook specifications for PlantClusterFinder plugins.
    """

    @hookspec
    def process_dataframe(self, df: pl.DataFrame):
        """
        Receive a Polars DataFrame from main.py
        and return any result.
        """

