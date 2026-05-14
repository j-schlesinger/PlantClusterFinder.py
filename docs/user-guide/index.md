# Using PlantClusterFinder.py

## Quick Start
If you don't have MCL and Diamond already installed, the fastest way to get started is using Conda:
```bash
# mamba or pixi are also ok
conda env create bgc bioconda::plantclusterfinder_py
```

## Data Prep
PlantClusterFinder.py expects three sources:
1. The genome data as a FASTA
2. the gene annotation as GFF3
3. the enzyme annotation in PGDB

The genome data can be obtained from [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) using [datasets](https://www.ncbi.nlm.nih.gov/datasets/), for example for Arabidopsis:
```bash
datasets download genome accession GCF_000001735 --include genome,gff3
```

The PGDB data can be obtained from [Plant Metabolic Network (PMN)](https://pmn.plantcyc.org/downloads), or created using Pathway-Tools Pathologic.

## Execution
PlantClusterFinder.py creates a `plantclusterfinder` executable in your path, which is enabled when you activate the bgc conda environment.
