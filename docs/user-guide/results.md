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

## GeneLocation Parquet Output
