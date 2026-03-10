"""
Plugin manager for PlantClusterFinder.

This module is responsible for:
- Initializing the Pluggy plugin manager
- Registering built-in and external plugins
- Executing plugin hooks in a controlled manner
"""

import pluggy
import polars as pl
from plugin_arch.hookspecs import PLantCLusterSpecs
from plugin_arch.example import ExamplePlugin

def get_plugin_manager():
    pm = pluggy.PluginManager("plantclusterfinder")
    
    # register hook specifications
    pm.add_hookspecs(PLantCLusterSpecs)

    # register plugins
    pm.register(ExamplePlugin())
    #pm.register(Plugin_2())
    
    return pm


def run_plugins(df: pl.DataFrame):

    pm = get_plugin_manager()

    results = pm.hook.process_dataframe(df=df)

    return results


if __name__ == "__main__":
    df = pl.DataFrame(
        {
            "gene": ["g1", None, "g3"],
            "enzyme": ["PAL", "C4H", "4CL"]
        }
    )

    results = run_plugins(df)

    print(results)