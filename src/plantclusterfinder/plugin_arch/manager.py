"""Plugin manager for PlantClusterFinder.

This module is responsible for:
- Initializing the Pluggy plugin manager
- Registering built-in and external plugins
- Executing plugin hooks in a controlled manner
"""

import logging
from pathlib import Path

import pluggy
import polars as pl

from plantclusterfinder.plugin_arch.context import RunContext
from plantclusterfinder.plugin_arch.hookspecs import PlantClusterSpecs


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

    def run(self, results: pl.DataFrame, context: RunContext) -> dict[str, Path]:
        """Run all registered plugins (process_results).

        Args:
            results: Results DataFrame as input for each plugin
            context: Context for run of all plugins
        """
        name_plugins = self.list_name_plugin()
        names = [n[0] for n in name_plugins]
        for required_outputs in self.hook.required_outputs():
            print(required_outputs)
            if any(req in names for req in required_outputs):
                raise NotImplementedError("Currently plugin heirarchy is not supported.")
            if any(req not in context.output_dict for req in required_outputs):
                raise pluggy.PluginValidationError(
                    self,
                    "Plugin requires output not available in the output dictionary: "
                    + str(required_outputs),
                )
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
