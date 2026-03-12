"""Plugin manager for PlantClusterFinder.

This module is responsible for:
- Initializing the Pluggy plugin manager
- Registering built-in and external plugins
- Executing plugin hooks in a controlled manner
"""

import pluggy

from plantclusterfinder_py.plugin_arch.context import RunContext
from plantclusterfinder_py.plugin_arch.hookspecs import PlantClusterSpecs


class PluginManager:
    """Central manager for loading and executing PlantClusterFinder plugins.

    This class abstracts Pluggy internals and provides a stable interface
    for the core pipeline to interact with plugins.
    """

    def __init__(self):
        """Create new PluginManager instance for PlantClusterFinder."""
        self.pm = pluggy.PluginManager("plantclusterfinder")
        self.pm.add_hookspecs(PlantClusterSpecs)

    def register(self, plugin):
        """Register a plugin implementation.

        Args:
            plugin: An object or class decorated with @hookimpl
        """
        self.pm.register(plugin)

    def load_entrypoints(self):
        """Load and register plugins exposed via setuptools entry points."""
        self.pm.load_setuptools_entrypoints("plantclusterfinder.plugins")

    def run(self, results, context: RunContext):
        """Run all registered plugins (process_results).

        Args:
            results: Results DataFrame as input for each plugin
            context: Context for run of all plugins
        """
        outputs = []
        for out in self.pm.hook.process_results(
            results=results,
            context=context,
        ):
            if out is not None:
                outputs.append(out)
        return outputs
