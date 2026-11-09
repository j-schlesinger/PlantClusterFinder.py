# Using PlantClusterFinder.next

## Quick Start

If you don't have MCL and Diamond already installed, the fastest way to get started is using Conda:

```bash
# mamba or pixi are also ok
conda env create bgc bioconda::plantclusterfinder_py
conda activate bgc
plantclusterfinder --test-install
```

## Data Prep

PlantClusterFinder.next expects three sources:

1. The genome data as a FASTA
2. the gene annotation as GFF3
3. the enzyme annotation in PGDB

The genome data can be obtained from [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) using [datasets](https://www.ncbi.nlm.nih.gov/datasets/), for example for Arabidopsis:

```bash
# Download the genome FASTA and annotation GFF3 for Arabidopsis thaliana
datasets download genome accession GCF_000001735 --include genome,gff3
```

The PGDB data can be obtained from [Plant Metabolic Network (PMN)](https://pmn.plantcyc.org/downloads), or created using Pathway-Tools Pathologic.

## Execution

PlantClusterFinder.next creates a `plantclusterfinder` executable in your path, which is enabled when you activate the bgc conda environment.

For details on the arguments that you need to predict plant BGCs, see the [Running PlantClusterFinder.next](arguments.md) section.
