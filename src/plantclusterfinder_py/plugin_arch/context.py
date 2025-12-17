from dataclasses import dataclass
from pathlib import Path

@dataclass
class RunContext:
    output_dir: Path
    config: dict
    genome_file: Path | None = None
    annotation_file: Path | None = None
