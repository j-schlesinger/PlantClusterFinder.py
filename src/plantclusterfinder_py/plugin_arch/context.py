"""Run Context object for Plugins."""

from pathlib import Path

from pydantic import BaseModel


class RunContext(BaseModel):
    """Run Context dataclass for plugin execution.

    Args:
        output_dir: Output directory for all plugins
        config: configuration dictionary for execution
        genome_file: FASTA Genome file
        annotation_file: GFF3 gene annotation file
    """

    output_dir: Path
    config: dict
    genome_file: Path | None = None
    annotation_file: Path | None = None
