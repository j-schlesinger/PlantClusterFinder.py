"""Run Context object for Plugins."""

from pathlib import Path

from pydantic import BaseModel


class RunContext(BaseModel):
    """Run Context dataclass for plugin execution.

    Args:
        run_name: Unique name of PlantClusterFinder execution
        output_dir: Output directory for all plugins
        config: configuration dictionary for execution
        genome_file: FASTA Genome file (optional)
        annotation_file: GFF3 gene annotation file (optional)
        output_dict: dictionary of output files
    """

    run_name: str
    output_dir: Path
    config: dict
    genome_file: Path | None = None
    annotation_file: Path | None = None
    output_dict: dict = {}


class ReturnContext(BaseModel):
    """Return Context object for plugin execution.

    Args:
        plugin_name: Name of the plugin output
        output_filename: location of filename of output
    """

    plugin_name: str
    output_filename: Path
