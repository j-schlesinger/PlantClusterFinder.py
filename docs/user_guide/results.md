# Examining PlantClusterFinder.py Results

There are two main files that PlantClusterFinder.py produces, the gene list, and the cluster list.
These are tab separated value files which have the data about the organism genes, and predicted BGCs respectively.

Additionally, these may be in one of three different "flavors" old, simple, and verbose.
By default, PlantClusterFinder.py produces "simple" output.

!!! danger
    The "old" flavor of output is for feature parity with PlantClusterFinder MATLAB **only**, and should be considered deprecated.

## Gene Output
The gene output file has a row for all annotated genes provided to PlantClusterFinder.py.
All "flavors" have the following output:

| Column | Type | Data |
|----|----|----|
| Gene Name | String | Name of annotated gene |
| Gene Start (bp) | Int | Start position of gene |
| Gene End (bp) | Int | End position of gene |
| Chromosome Name | String | Name of Chromosome or Scaffold |
| Enzyme | Boolean | Whether or not gene is an enzyme |
| EC-NUMBER | String | Comma separated list of EC numbers (empty if not enzyme) |
| RXN | String | Reaction ID from PGDB/PMN |
| Locally Duplicated | Boolean | True if gene is considered locally duplicated |
| LD_ClusterIds | String | Comma separated list of IDs used by MCL to determine local duplication |
| ... | ... | Section of flavor-dependent columns |
| ... | Boolean | Column for every selected metabolic domain, true if enzyme is in metabolic domain |
| Signature Tailoring Enzyme Class | String | Comma separated list of if the enzyme is signature, tailoring, or both |
| Signature or Tailoring | String | Comma Sperated list of if the enzyme is signature or tailoring |

### Flavor-Dependent Columns
The "simple" flavor adds the Boolean columns "Stepsize (finalStepSize)", and "Top (top percent clusters)%" which are true if the gene is picked up in the final stepsize or in the top percent clusters.

The "verbose" flavor additionally adds boolean columns for all step sizes.

The "old" flavor has a column for every combination of top percent clusters and step sizes.

## Cluster Output

| Column | Type | Data |
|----|----|----|
| Stepsize | Int | Step size used when finding the cluster |
| Cluster ID | Int | Cluster identifier from this run -- not stable across runs |
| Top n% cluster | Boolean | |
| Chromosome | String | Chromosome or scaffold name |
| Start bp | Int | Start of first gene on cluster |
| End bp | Int | End of last gene on cluster |
| Physical size bp | Int | Distance between the start and end genes |
| Number of Genes | Int | Number of genes in cluster |
| Number of Hypothetical Genes | Int | Number of hypothetical genes that make up cluster |
| Number of Enzymes | Int | Number of enzymes in cluster |
| Start Gene | String | First gene in cluster |
| End Gene | String | Last gene in cluster |
| Gene IDs | String | Comma separated list of all genes in cluster |

### Flavor-Dependent Rows
The "simple" flavor returns clusters where the Stepsize is the "FinalStepSize" and "Top n% cluster" is true.

The "verbose" flavor and "old" flavor return clusters where Stepsize is less than or equal to the "FinalStepSize"
