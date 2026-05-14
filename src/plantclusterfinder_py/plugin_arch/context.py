"""Run Context object for Plugins."""

from pathlib import Path

from pydantic import BaseModel


class RunContext(BaseModel):
    """Run Context dataclass for plugin execution."""

    """Unique name of PlantClusterFinder execution"""
    run_name: str
    """Output directory for all plugins"""
    output_dir: Path
    """configuration dictionary for execution"""
    config: dict
    """FASTA Genome file (optional)"""
    genome_file: Path | None = None
    """GFF3 gene annotation file (optional)"""
    annotation_file: Path | None = None
    """dictionary of output files"""
    output_dict: dict = {}


class ReturnContext(BaseModel):
    """Return Context object for plugin execution."""

    """Name of the plugin output."""
    plugin_name: str
    """Location of the output file."""
    output_filename: Path
