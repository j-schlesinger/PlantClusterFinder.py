"""Plugin manager for PlantClusterFinder.

This module is responsible for:
- Initializing the Pluggy plugin manager
- Registering built-in and external plugins
- Executing plugin hooks in a controlled manner
"""

import pluggy
import polars as pl

from plantclusterfinder_py.plugin_arch.example import ExamplePlugin
from plantclusterfinder_py.plugin_arch.hookspecs import PlantClusterSpecs


def get_plugin_manager():
    """Get the plugin manager for PlantClusterFinder."""
    pm = pluggy.PluginManager("plantclusterfinder")

    # register hook specifications
    pm.add_hookspecs(PlantClusterSpecs)

    # register plugins
    pm.register(ExamplePlugin())

    return pm


def run_plugins(df: pl.DataFrame):
    """Run the plugins defined by the input DataFrame."""
    pm = get_plugin_manager()

    return pm.hook.process_dataframe(df=df)


if __name__ == "__main__":
    df = pl.DataFrame({"gene": ["g1", None, "g3"], "enzyme": ["PAL", "C4H", "4CL"]})

    results = run_plugins(df)

    print(results)
