"""Plugin manager for PlantClusterFinder.

This module is responsible for:
- Initializing the Pluggy plugin manager
- Registering built-in and external plugins
- Executing plugin hooks in a controlled manner
"""

import logging
from pathlib import Path

import pluggy

from plantclusterfinder_py.plugin_arch.context import RunContext
from plantclusterfinder_py.plugin_arch.hookspecs import PlantClusterSpecs


class PluginManager(pluggy.PluginManager):
    """Central manager for loading and executing PlantClusterFinder plugins.

    This class abstracts Pluggy internals and provides a stable interface
    for the core pipeline to interact with plugins.
    """

    def __init__(self):
        """Create new PluginManager instance for PlantClusterFinder."""
        super().__init__("plantclusterfinder")
        self.add_hookspecs(PlantClusterSpecs)
        self._log = logging.getLogger()

    def load_entrypoints(self):
        """Load and register plugins exposed via setuptools entry points."""
        super().load_setuptools_entrypoints("plantclusterfinder.plugins")

    def run(self, results, context: RunContext) -> dict[str, Path]:
        """Run all registered plugins (process_results).

        Args:
            results: Results DataFrame as input for each plugin
            context: Context for run of all plugins
        """
        for out in self.hook.process_results(
            results=results,
            context=context,
        ):
            if out is None:
                self._log.warning("Output not found for a plugin")
                continue
            k = out.plugin_name
            v = out.output_filename
            if k in context.output_dict:
                self._log.warning(
                    "Overwriting plugin output value for %s, was %s", k, context.output_dict[k]
                )
                context.output_dict[k + "_old"] = context.output_dict[k]
            context.output_dict[k] = v
        return context.output_dict
