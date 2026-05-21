# PlantClusterFinder
---
PlantClusterFinder (PCF) detects metabolic gene clusters in a sequenced genome.
It uses a gene location file provided by the user (see below) and a PGDB created with Pathway Tools as well as further information (see below) to identify enzyme-coding genes (metabolic genes) located together on a chromosome.
Initially only continuous stretches of metabolic genes lying directly next to each other are allowed. This condition is relaxed by iteratively increasing the intervening (non-metabolic) gene size by one. Several criteria to select for clusters are provided. In addition to this, clusters can be prevented from forming by a section of criteria.
Details of PCF (version 1.0) can be found in PMID: 28228535.

To get started, please visit the [user-guide](user-guide/install.md).
